/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* HDR composition of registered linear images taken with different exposures.
 *
 * Each image is mapped into the units of the longest exposure by a linear
 * relation measured on the pixels well exposed in both members of each
 * consecutive pair, then the images are averaged with noise-optimal weights
 * derived from a Poisson-Gaussian noise model estimated on each image,
 * excluding saturated and non-linear pixels through feathered masks.
 *
 * References:
 * - M. Granados, B. Ajdin, M. Wand, C. Theobalt, H.-P. Seidel, H. Lensch,
 *   "Optimal HDR reconstruction with linear digital cameras", CVPR 2010.
 * - C. Aguerrebere, J. Delon, Y. Gousseau, P. Musé, "Best algorithms for HDR
 *   image generation. A study of performance bounds", SIAM J. Imaging Sci. 2014.
 * - P. Hanji, F. Zhong, R. K. Mantiuk, "Noise-aware merging of high dynamic
 *   range image stacks without camera calibration", ECCV Workshops 2020.
 * - A. Foi, M. Trimeche, V. Katkovnik, K. Egiazarian, "Practical
 *   Poissonian-Gaussian noise modeling and fitting for single-image raw-data",
 *   IEEE Trans. Image Process. 2008.
 * - M. Granados, K. I. Kim, J. Tompkin, C. Theobalt, "Automatic noise
 *   modeling for ghost-free HDR reconstruction", ACM Trans. Graph. 2013.
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "algos/PSF.h"
#include "algos/sorting.h"
#include "compositing/hdr.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"

#define NOISE_BLOCK 16
#define FIT_MAX_SAMPLES 2000000
#define FIT_MIN_SAMPLES 500
#define FIT_LOW_SIGMA 5.0	// lower bound of the fit window, above the background
#define SAT_LEVEL 0.95		// fraction of each channel maximum considered saturated
#define FIT_HIGH 0.85		// upper bound of the fit window, fraction of the saturation level
#define TUKEY_C 4.685
#define DEFICIT_SIGMA 5.0	// significance of a deficit against the shorter exposure
#define DEFICIT_LEVEL 0.5	// deficits are only looked for above this fraction of the maximum
#define REJECT_SIGMA 5.0
#define RATIO_WARN 0.1

typedef struct {
	fits *fit;
	const char *name;
	double exposure;	// seconds, 0 if unknown
	double eff;		// exposure times ISO when available, used for ordering and ratios
	double signal;		// sum of mean - background over channels, ordering fallback
	double bg[3], max[3], sat[3];
	double s2[3], k[3];	// variance = s2 + k * max(x - bg, 0)
	double a, b[3];		// y = a * x + b, in units of the longest exposure
	double c, d[3];		// x = c * x_prev + d, against the previous (longer) image
	uint8_t *mask;
	float *h;		// validity weight in [0, 1], NULL means always valid
} hdr_img;

void hdr_params_init(struct hdr_params *params) {
	params->feather = HDR_DEFAULT_FEATHER;
}

void free_hdr_args(struct hdr_args *args) {
	if (!args)
		return;
	g_slist_free_full(args->files, g_free);
	free(args);
}

static double sampled_median(const float *x, size_t n) {
	size_t stride = n > 1000000 ? n / 1000000 : 1;
	size_t m = n / stride;
	float *tmp = malloc(m * sizeof(float));
	if (!tmp)
		return 0.0;
	for (size_t i = 0; i < m; i++)
		tmp[i] = x[i * stride];
	double med = quickmedian_float(tmp, m);
	free(tmp);
	return med;
}

static void channel_stats(hdr_img *im, int nch, size_t npix) {
	im->signal = 0.0;
	for (int ch = 0; ch < nch; ch++) {
		const float *x = im->fit->fpdata[ch];
		double sum = 0.0;
		float mx = -INFINITY;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(+:sum) reduction(max:mx)
#endif
		for (size_t p = 0; p < npix; p++) {
			sum += x[p];
			if (x[p] > mx)
				mx = x[p];
		}
		im->bg[ch] = sampled_median(x, npix);
		im->max[ch] = mx;
		im->signal += sum / npix - im->bg[ch];
	}
}

/* Poisson-Gaussian noise model of a channel, in the spirit of Foi et al.:
 * local variances of a high-pass residual on small blocks, grouped by local
 * mean, then a robust fit of variance = s2 + k * (mean - bg). */
struct nblock {
	float mu, var;
};

static int cmp_nblock(const void *a, const void *b) {
	float x = ((const struct nblock *)a)->mu, y = ((const struct nblock *)b)->mu;
	return (x > y) - (x < y);
}

static int cmp_float(const void *a, const void *b) {
	float x = *(const float *)a, y = *(const float *)b;
	return (x > y) - (x < y);
}

static void fit_noise_model(const double *u, const double *v, int ng, double *s2_out, double *k_out) {
	int nlow = max(1, ng / 4);
	double *tmp = malloc(nlow * sizeof(double));
	memcpy(tmp, v, nlow * sizeof(double));
	double s2 = quickmedian_double(tmp, nlow), k = 0.0;
	free(tmp);
	for (int it = 0; it < 30; it++) {
		double sw = 0.0, swu = 0.0, swv = 0.0, swuu = 0.0, swuv = 0.0;
		for (int g = 0; g < ng; g++) {
			double pred = s2 + k * u[g];
			if (pred <= 0.0)
				continue;
			if (it > 0 && fabs(v[g] / pred - 1.0) > 1.0)
				continue;
			double w = 1.0 / (pred * pred);
			sw += w;
			swu += w * u[g];
			swv += w * v[g];
			swuu += w * u[g] * u[g];
			swuv += w * u[g] * v[g];
		}
		double det = sw * swuu - swu * swu;
		if (sw <= 0.0)
			break;
		double nk = det > 0.0 ? (sw * swuv - swu * swv) / det : 0.0;
		double ns2;
		if (nk < 0.0) {
			nk = 0.0;
			ns2 = swv / sw;
		} else {
			ns2 = (swv - nk * swu) / sw;
		}
		gboolean done = fabs(ns2 - s2) <= 1e-6 * fabs(s2) && fabs(nk - k) <= 1e-6 * fabs(k) + 1e-15;
		s2 = ns2;
		k = nk;
		if (done)
			break;
	}
	if (s2 <= 0.0) {
		s2 = INFINITY;
		for (int g = 0; g < ng; g++)
			if (v[g] > 0.0 && v[g] < s2)
				s2 = v[g];
		if (!isfinite(s2))
			s2 = 1e-12;
	}
	*s2_out = s2;
	*k_out = k;
}

static void estimate_noise(hdr_img *im, int ch, int rx, int ry) {
	const float *x = im->fit->fpdata[ch];
	const float sat = (float) im->sat[ch];
	int bw = (rx - 2) / NOISE_BLOCK, bh = (ry - 2) / NOISE_BLOCK;
	size_t nb = (size_t) bw * bh;
	im->s2[ch] = 1e-12;
	im->k[ch] = 0.0;
	if (nb < 16)
		return;
	struct nblock *blocks = malloc(nb * sizeof(struct nblock));
	if (!blocks)
		return;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (size_t b = 0; b < nb; b++) {
		float d[NOISE_BLOCK * NOISE_BLOCK];
		int x0 = 1 + (int) (b % bw) * NOISE_BLOCK, y0 = 1 + (int) (b / bw) * NOISE_BLOCK;
		double sum = 0.0;
		int n = 0;
		gboolean valid = TRUE;
		for (int yy = y0; yy < y0 + NOISE_BLOCK && valid; yy++) {
			for (int xx = x0; xx < x0 + NOISE_BLOCK; xx++) {
				const float *c = x + (size_t) yy * rx + xx;
				if (*c >= sat) {
					valid = FALSE;
					break;
				}
				float s9 = c[-rx - 1] + c[-rx] + c[-rx + 1] + c[-1] + c[0] + c[1] + c[rx - 1] + c[rx] + c[rx + 1];
				d[n++] = *c - s9 / 9.f;
				sum += *c;
			}
		}
		if (!valid) {
			blocks[b].var = -1.f;
			continue;
		}
		float med = quickmedian_float(d, n);
		for (int i = 0; i < n; i++)
			d[i] = fabsf(d[i] - med);
		/* the 3x3 high-pass keeps 8/9 of the variance of white noise */
		double mad = 1.4826 * quickmedian_float(d, n);
		blocks[b].var = (float) (mad * mad * 9.0 / 8.0);
		blocks[b].mu = (float) (sum / n);
	}

	size_t nvalid = 0;
	for (size_t b = 0; b < nb; b++)
		if (blocks[b].var >= 0.f)
			blocks[nvalid++] = blocks[b];
	if (nvalid < 16) {
		free(blocks);
		return;
	}
	qsort(blocks, nvalid, sizeof(struct nblock), cmp_nblock);

	size_t gsize = max(32, nvalid / 256);
	int ng = (int) (nvalid / gsize);
	if (ng < 1) ng = 1;
	double *u = malloc(ng * sizeof(double));
	double *v = malloc(ng * sizeof(double));
	float *tmp = malloc(gsize * sizeof(float));
	for (int g = 0; g < ng; g++) {
		struct nblock *grp = blocks + g * gsize;
		for (size_t i = 0; i < gsize; i++)
			tmp[i] = grp[i].var;
		qsort(tmp, gsize, sizeof(float), cmp_float);
		/* lower quartile: blocks with real structure only inflate the variance */
		v[g] = tmp[gsize / 4];
		u[g] = max(0.0, grp[gsize / 2].mu - im->bg[ch]);
	}
	fit_noise_model(u, v, ng, &im->s2[ch], &im->k[ch]);
	free(tmp);
	free(u);
	free(v);
	free(blocks);
}

static void dilate_line(const uint8_t *in, uint8_t *out, int n, size_t stride, int r) {
	int count = 0;
	for (int i = 0; i < min(r, n); i++)
		count += in[i * stride] != 0;
	for (int i = 0; i < n; i++) {
		if (i + r < n)
			count += in[(i + r) * stride] != 0;
		if (i - r - 1 >= 0)
			count -= in[(i - r - 1) * stride] != 0;
		out[i * stride] = count > 0;
	}
}

static int dilate_mask(uint8_t *mask, int rx, int ry, int r) {
	uint8_t *tmp = malloc((size_t) rx * ry);
	if (!tmp)
		return 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int y = 0; y < ry; y++)
		dilate_line(mask + (size_t) y * rx, tmp + (size_t) y * rx, rx, 1, r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int x = 0; x < rx; x++)
		dilate_line(tmp + x, mask + x, ry, rx, r);
	free(tmp);
	return 0;
}

static void blur_line(const float *in, float *out, int n, size_t stride, int r) {
	double sum = 0.0;
	int count = 0;
	for (int i = 0; i < min(r, n); i++) {
		sum += in[i * stride];
		count++;
	}
	for (int i = 0; i < n; i++) {
		if (i + r < n) {
			sum += in[(i + r) * stride];
			count++;
		}
		if (i - r - 1 >= 0) {
			sum -= in[(i - r - 1) * stride];
			count--;
		}
		out[i * stride] = (float) (sum / count);
	}
}

static int box_blur(float *f, int rx, int ry, int r) {
	float *tmp = malloc((size_t) rx * ry * sizeof(float));
	if (!tmp)
		return 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int y = 0; y < ry; y++)
		blur_line(f + (size_t) y * rx, tmp + (size_t) y * rx, rx, 1, r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int x = 0; x < rx; x++)
		blur_line(tmp + x, f + x, ry, rx, r);
	free(tmp);
	return 0;
}

/* The mask is grown by r pixels, then smoothed over r pixels so that the
 * originally masked pixels keep a null weight. */
static float *feather_mask(uint8_t *mask, int rx, int ry, int r) {
	size_t npix = (size_t) rx * ry;
	float *h = malloc(npix * sizeof(float));
	if (!h)
		return NULL;
	if (r > 0 && dilate_mask(mask, rx, ry, r)) {
		free(h);
		return NULL;
	}
	for (size_t p = 0; p < npix; p++)
		h[p] = mask[p] ? 1.f : 0.f;
	int r1 = (r + 1) / 2, r2 = r / 2;
	if ((r1 > 0 && box_blur(h, rx, ry, r1)) || (r2 > 0 && box_blur(h, rx, ry, r2))) {
		free(h);
		return NULL;
	}
	for (size_t p = 0; p < npix; p++)
		h[p] = 1.f - h[p];
	return h;
}

static float channel_median(float *v, const uint8_t *chan, size_t n, int ch, float *tmp) {
	size_t m = 0;
	for (size_t i = 0; i < n; i++)
		if (chan[i] == ch)
			tmp[m++] = v[i];
	return m ? (float) quickmedian_float(tmp, m) : 0.f;
}

/* Robust fit of S = c * L + d_ch on the pixels where L is well exposed.
 * Regressing the short exposure on the long one keeps the fit unbiased:
 * the selection is made on the regressor, which is also the less noisy. */
static int fit_pair(const hdr_img *L, const hdr_img *S, const uint8_t *excl, int nch, size_t npix,
		double c0, double *c_out, double *d_out, size_t *nused) {
	size_t stride = npix * nch / FIT_MAX_SAMPLES + 1;
	size_t cap = (npix / stride + 1) * nch;
	float *xs = malloc(cap * sizeof(float));
	float *ys = malloc(cap * sizeof(float));
	float *tmp = malloc(cap * sizeof(float));
	float *w = malloc(cap * sizeof(float));
	uint8_t *chan = malloc(cap);
	double lo[3], hi[3];
	size_t n = 0;
	int retval = 1;
	*nused = 0;
	if (!xs || !ys || !tmp || !w || !chan)
		goto end;

	for (int ch = 0; ch < nch; ch++) {
		lo[ch] = L->bg[ch] + FIT_LOW_SIGMA * sqrt(L->s2[ch]);
		hi[ch] = FIT_HIGH * L->sat[ch];
	}
	for (size_t p = 0; p < npix; p += stride) {
		if (excl[p])
			continue;
		for (int ch = 0; ch < nch; ch++) {
			float l = L->fit->fpdata[ch][p];
			if (l > lo[ch] && l < hi[ch]) {
				xs[n] = l;
				ys[n] = S->fit->fpdata[ch][p];
				chan[n++] = ch;
			}
		}
	}
	if (n < FIT_MIN_SAMPLES)
		goto end;

	double c = c0, d[3] = { 0.0 };
	if (c <= 0.0) {
		size_t m = 0;
		for (size_t i = 0; i < n; i++)
			tmp[m++] = (ys[i] - S->bg[chan[i]]) / (xs[i] - L->bg[chan[i]]);
		c = quickmedian_float(tmp, m);
	}
	for (int ch = 0; ch < nch; ch++) {
		for (size_t i = 0; i < n; i++)
			w[i] = ys[i] - c * xs[i];
		d[ch] = channel_median(w, chan, n, ch, tmp);
	}

	for (int it = 0; it < 50; it++) {
		double scale[3];
		for (size_t i = 0; i < n; i++)
			w[i] = fabsf(ys[i] - (float) (c * xs[i] + d[chan[i]]));
		for (int ch = 0; ch < nch; ch++) {
			scale[ch] = 1.4826 * channel_median(w, chan, n, ch, tmp);
			if (scale[ch] <= 0.0)
				scale[ch] = 1e-9;
		}
		double sw[3] = { 0.0 }, sx[3] = { 0.0 }, sy[3] = { 0.0 };
		for (size_t i = 0; i < n; i++) {
			double r = w[i] / (TUKEY_C * scale[chan[i]]);
			w[i] = r < 1.0 ? (float) ((1.0 - r * r) * (1.0 - r * r)) : 0.f;
			sw[chan[i]] += w[i];
			sx[chan[i]] += w[i] * xs[i];
			sy[chan[i]] += w[i] * ys[i];
		}
		double xm[3], ym[3];
		for (int ch = 0; ch < nch; ch++) {
			if (sw[ch] <= 0.0)
				goto end;
			xm[ch] = sx[ch] / sw[ch];
			ym[ch] = sy[ch] / sw[ch];
		}
		double num = 0.0, den = 0.0;
		for (size_t i = 0; i < n; i++) {
			double dx = xs[i] - xm[chan[i]];
			num += w[i] * dx * (ys[i] - ym[chan[i]]);
			den += w[i] * dx * dx;
		}
		if (den <= 0.0)
			goto end;
		double nc = num / den;
		for (int ch = 0; ch < nch; ch++)
			d[ch] = ym[ch] - nc * xm[ch];
		gboolean done = fabs(nc - c) <= 1e-7 * fabs(c);
		c = nc;
		if (done && it > 0)
			break;
	}
	if (!(c > 0.0) || !isfinite(c))
		goto end;
	for (size_t i = 0; i < n; i++)
		if (w[i] > 0.f)
			(*nused)++;
	if (*nused < FIT_MIN_SAMPLES)
		goto end;
	*c_out = c;
	for (int ch = 0; ch < nch; ch++)
		d_out[ch] = d[ch];
	retval = 0;

end:
	free(xs);
	free(ys);
	free(tmp);
	free(w);
	free(chan);
	return retval;
}

static int cmp_hdr_img(const void *a, const void *b) {
	const hdr_img *x = a, *y = b;
	if (x->eff != y->eff)
		return x->eff < y->eff ? 1 : -1;
	return (x->signal < y->signal) - (x->signal > y->signal);
}

static inline double var_model(const hdr_img *im, int ch, double x) {
	return im->s2[ch] + im->k[ch] * max(x - im->bg[ch], 0.0);
}

/* Pixels of a longer exposure fall short of what the next shorter exposure
 * predicts when they are saturated, non-linear or bloomed. */
static void flag_deficits(hdr_img *L, const hdr_img *S, int nch, size_t npix) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (size_t p = 0; p < npix; p++) {
		if (L->mask[p])
			continue;
		for (int ch = 0; ch < nch; ch++) {
			double xs = S->fit->fpdata[ch][p];
			double pred = (xs - S->d[ch]) / S->c;
			if (pred < DEFICIT_LEVEL * L->max[ch])
				continue;
			double v = var_model(L, ch, pred) + var_model(S, ch, xs) / (S->c * S->c);
			if (L->fit->fpdata[ch][p] - pred < -DEFICIT_SIGMA * sqrt(v)) {
				L->mask[p] = 1;
				break;
			}
		}
	}
}

/* Relative sensitivity of each image, in ADU per unit of exposure time. The
 * GAIN keyword has vendor-specific units and is never converted: only the
 * e-/ADU value can account for different gains. */
static gboolean sensitivity_factors(fits **images, int nb, double *factor) {
	gboolean iso = TRUE, cvf = TRUE, gains_differ = FALSE, same_cvf = TRUE;
	guint gain = DEFAULT_UINT_VALUE;
	for (int i = 0; i < nb; i++) {
		const fits *f = images[i];
		if (f->keywords.iso_speed <= 0.0)
			iso = FALSE;
		if (f->keywords.cvf <= 0.0)
			cvf = FALSE;
		else if (f->keywords.cvf != images[0]->keywords.cvf)
			same_cvf = FALSE;
		if (f->keywords.key_gain != DEFAULT_UINT_VALUE) {
			if (gain == DEFAULT_UINT_VALUE)
				gain = f->keywords.key_gain;
			else if (f->keywords.key_gain != gain)
				gains_differ = TRUE;
		}
	}
	/* some drivers write the same e-/ADU value whatever the gain */
	if (cvf && !(gains_differ && same_cvf)) {
		for (int i = 0; i < nb; i++)
			factor[i] = 1.0 / images[i]->keywords.cvf;
		if (gains_differ)
			siril_log_message(_("Images taken with different gains: using their e-/ADU values\n"));
		return TRUE;
	}
	if (!gains_differ) {
		for (int i = 0; i < nb; i++)
			factor[i] = iso ? images[i]->keywords.iso_speed : 1.0;
		return TRUE;
	}
	siril_log_warning(_("Images taken with different gains and without usable e-/ADU values: images are ordered by "
				"signal level and the exposure ratios rely on the measurements only\n"));
	return FALSE;
}

static void free_hdr_imgs(hdr_img *im, int nb) {
	for (int i = 0; i < nb; i++) {
		free(im[i].mask);
		free(im[i].h);
	}
	free(im);
}

fits *hdr_compose(fits **images, const char **names, int nb, const struct hdr_params *params) {
	if (nb < 2) {
		siril_log_error(_("HDR composition needs at least two images\n"));
		return NULL;
	}
	int rx = images[0]->rx, ry = images[0]->ry, nch = (int) images[0]->naxes[2];
	for (int i = 0; i < nb; i++) {
		if (images[i]->type != DATA_FLOAT) {
			siril_log_error(_("HDR composition requires 32-bit images\n"));
			return NULL;
		}
		if ((int) images[i]->rx != rx || (int) images[i]->ry != ry || images[i]->naxes[2] != nch) {
			siril_log_error(_("HDR composition: all images must have the same size and number of channels\n"));
			return NULL;
		}
	}
	size_t npix = (size_t) rx * ry;
	fits *out = NULL;
	hdr_img *im = calloc(nb, sizeof(hdr_img));
	if (!im) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	gboolean known = TRUE;
	for (int i = 0; i < nb; i++) {
		im[i].fit = images[i];
		im[i].name = names ? names[i] : "";
		im[i].exposure = images[i]->keywords.exposure > 0.0 ? images[i]->keywords.exposure : 0.0;
		if (im[i].exposure <= 0.0)
			known = FALSE;
		channel_stats(&im[i], nch, npix);
	}
	double *factor = malloc(nb * sizeof(double));
	if (!factor) {
		PRINT_ALLOC_ERR;
		free(im);
		return NULL;
	}
	if (!known)
		siril_log_warning(_("Exposure time missing from at least one header: images are ordered by signal level\n"));
	else known = sensitivity_factors(images, nb, factor);
	for (int i = 0; i < nb; i++)
		im[i].eff = known ? im[i].exposure * factor[i] : 0.0;
	free(factor);
	qsort(im, nb, sizeof(hdr_img), cmp_hdr_img);
	if (images[0]->keywords.bayer_pattern[0] != '\0' && nch == 1)
		siril_log_warning(_("HDR composition of undebayered CFA images: masks and noise estimates mix the colour channels\n"));

	gui_iface.set_progress(PROGRESS_RESET, _("HDR composition: noise estimation"));
	for (int i = 0; i < nb; i++) {
		double noise = 0.0, k = 0.0;
		for (int ch = 0; ch < nch; ch++) {
			im[i].sat[ch] = SAT_LEVEL * im[i].max[ch];
			estimate_noise(&im[i], ch, rx, ry);
			noise += sqrt(im[i].s2[ch]) / nch;
			k += im[i].k[ch] / nch;
		}
		if (im[i].exposure > 0.0)
			siril_log_message(_("HDR input %d: %s, exposure %g s, noise %.3e, shot-noise factor %.3e\n"),
					i + 1, im[i].name, im[i].exposure, noise, k);
		else siril_log_message(_("HDR input %d: %s, unknown exposure, noise %.3e, shot-noise factor %.3e\n"),
					i + 1, im[i].name, noise, k);
		gui_iface.set_progress(0.3 * (i + 1) / nb, NULL);
	}

	/* the shortest exposure is the last resort and is never masked */
	for (int i = 0; i < nb - 1; i++) {
		im[i].mask = calloc(npix, 1);
		if (!im[i].mask) {
			PRINT_ALLOC_ERR;
			goto end;
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t p = 0; p < npix; p++)
			for (int ch = 0; ch < nch; ch++)
				if (im[i].fit->fpdata[ch][p] >= im[i].sat[ch]) {
					im[i].mask[p] = 1;
					break;
				}
	}

	gui_iface.set_progress(0.3, _("HDR composition: exposure ratios"));
	im[0].a = 1.0;
	for (int i = 1; i < nb; i++) {
		hdr_img *L = &im[i - 1], *S = &im[i];
		double c0 = known ? S->eff / L->eff : 0.0;
		size_t nused = 0;
		uint8_t *excl = malloc(npix);
		if (!excl) {
			PRINT_ALLOC_ERR;
			goto end;
		}
		memcpy(excl, L->mask, npix);
		int err = dilate_mask(excl, rx, ry, max(params->feather, 1));
		if (!err)
			err = fit_pair(L, S, excl, nch, npix, c0, &S->c, S->d, &nused);
		free(excl);
		if (err) {
			if (!known) {
				siril_log_error(_("HDR composition: not enough well exposed pixels in common between %s and %s "
							"to measure their exposure ratio\n"), L->name, S->name);
				goto end;
			}
			siril_log_warning(_("Not enough well exposed pixels in common between %s and %s, "
						"using the header exposure ratio and the background levels\n"), L->name, S->name);
			S->c = c0;
			for (int ch = 0; ch < nch; ch++)
				S->d[ch] = S->bg[ch] - c0 * L->bg[ch];
		} else if (known) {
			siril_log_message(_("Exposure ratio %s / %s: %.4f measured on %zu pixels, %.4f from headers\n"),
					L->name, S->name, 1.0 / S->c, nused, 1.0 / c0);
			if (fabs(S->c / c0 - 1.0) > RATIO_WARN)
				siril_log_warning(_("The measured ratio differs from the header ratio by more than %d%%: "
							"check the exposure times, gains and sky transparency\n"), (int) (RATIO_WARN * 100.0));
		} else {
			siril_log_message(_("Exposure ratio %s / %s: %.4f measured on %zu pixels\n"),
					L->name, S->name, 1.0 / S->c, nused);
		}
		S->a = L->a / S->c;
		for (int ch = 0; ch < nch; ch++)
			S->b[ch] = L->b[ch] - L->a * S->d[ch] / S->c;
	}

	gui_iface.set_progress(0.5, _("HDR composition: saturation masks"));
	for (int i = 0; i < nb - 1; i++) {
		flag_deficits(&im[i], &im[i + 1], nch, npix);
		size_t nmasked = 0;
		for (size_t p = 0; p < npix; p++)
			nmasked += im[i].mask[p];
		siril_log_message(ngettext("%s: %zu pixel (%.3g%%) is saturated or non-linear\n",
					"%s: %zu pixels (%.3g%%) are saturated or non-linear\n", nmasked),
				im[i].name, nmasked, 100.0 * nmasked / npix);
		im[i].h = feather_mask(im[i].mask, rx, ry, params->feather);
		if (!im[i].h) {
			PRINT_ALLOC_ERR;
			goto end;
		}
	}

	gui_iface.set_progress(0.7, _("HDR composition: merging"));
	if (new_fit_image(&out, rx, ry, nch, DATA_FLOAT)) {
		out = NULL;
		goto end;
	}
	float maxv = 0.f;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:maxv)
#endif
	for (size_t p = 0; p < npix; p++) {
		for (int ch = 0; ch < nch; ch++) {
			double num = 0.0, den = 0.0;
			for (int i = 0; i < nb; i++) {
				double hw = im[i].h ? im[i].h[p] : 1.0;
				if (hw <= 0.0)
					continue;
				double y = im[i].a * im[i].fit->fpdata[ch][p] + im[i].b[ch];
				double w = hw / (im[i].a * im[i].a * max(im[i].s2[ch], 1e-12));
				num += w * y;
				den += w;
			}
			double y0 = num / den;
			num = den = 0.0;
			/* weights from the estimated radiance rather than from each sample */
			for (int i = 0; i < nb; i++) {
				double hw = im[i].h ? im[i].h[p] : 1.0;
				if (hw <= 0.0)
					continue;
				double a = im[i].a;
				double y = a * im[i].fit->fpdata[ch][p] + im[i].b[ch];
				double v = a * a * max(var_model(&im[i], ch, (y0 - im[i].b[ch]) / a), 1e-12);
				if (fabs(y - y0) > REJECT_SIGMA * sqrt(v))
					continue;
				num += hw / v * y;
				den += hw / v;
			}
			float yf = (float) (den > 0.0 ? num / den : y0);
			out->fpdata[ch][p] = yf;
			if (yf > maxv)
				maxv = yf;
		}
	}
	if (!(maxv > 0.f)) {
		siril_log_error(_("HDR composition: empty result\n"));
		clearfits(out);
		free(out);
		out = NULL;
		goto end;
	}
	float norm = 1.f / maxv;
	for (size_t p = 0; p < npix * nch; p++)
		out->fdata[p] = max(out->fdata[p] * norm, 0.f);

	fits **sorted = malloc((nb + 1) * sizeof(fits *));
	for (int i = 0; i < nb; i++)
		sorted[i] = im[i].fit;
	sorted[nb] = NULL;
	merge_fits_headers_to_result2(out, sorted, TRUE);
	free(sorted);
	/* pixel values now correspond to this exposure of the longest image */
	out->keywords.exposure = im[0].exposure > 0.0 ? im[0].exposure / maxv : 0.0;
	GString *hist = g_string_new("HDR composition of");
	for (int i = 0; i < nb; i++) {
		if (im[i].exposure > 0.0)
			g_string_append_printf(hist, " %gs", im[i].exposure);
		else g_string_append(hist, " ?s");
	}
	out->history = g_slist_append(out->history, g_string_free(hist, FALSE));
	if (out->keywords.exposure > 0.0)
		siril_log_message(_("HDR composition done: dynamic range extended by a factor %.1f, "
					"equivalent exposure %g s\n"), im[nb - 1].a, out->keywords.exposure);
	else siril_log_message(_("HDR composition done: dynamic range extended by a factor %.1f\n"), im[nb - 1].a);
	gui_iface.set_progress(PROGRESS_DONE, _("HDR composition done"));

end:
	free_hdr_imgs(im, nb);
	if (!out)
		gui_iface.set_progress(PROGRESS_DONE, _("HDR composition failed"));
	return out;
}

/* gfit must only be replaced on the main thread, see install_new_single_image() */
static gboolean install_hdr_result(gpointer p) {
	clear_stars_list(TRUE);
	if (!install_new_single_image(p, strdup(_("Unsaved HDR composition"))) && com.uniq)
		com.uniq->comment = strdup(_("HDR composition"));
	return FALSE;
}

gpointer hdr_worker(gpointer p) {
	struct hdr_args *args = (struct hdr_args *) p;
	int nb = g_slist_length(args->files);
	fits **images = calloc(nb, sizeof(fits *));
	gchar **names = calloc(nb + 1, sizeof(gchar *));
	fits *out = NULL;
	int retval = 1, i = 0;

	for (GSList *l = args->files; l; l = l->next, i++) {
		images[i] = calloc(1, sizeof(fits));
		if (!images[i] || readfits(l->data, images[i], NULL, TRUE)) {
			siril_log_error(_("Could not open %s\n"), (char *) l->data);
			goto end;
		}
		names[i] = g_path_get_basename(l->data);
	}
	out = hdr_compose(images, (const char **) names, nb, &args->params);
	if (!out)
		goto end;

	update_fits_header(out);
	gui_iface.execute_idle_sync(install_hdr_result, out);
	out = NULL;
	retval = 0;

end:
	for (i = 0; i < nb; i++) {
		if (images[i]) {
			clearfits(images[i]);
			free(images[i]);
		}
	}
	free(images);
	g_strfreev(names);
	if (out) {
		clearfits(out);
		free(out);
	}
	free_hdr_args(args);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}
