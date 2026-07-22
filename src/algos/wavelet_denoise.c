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

/* Per-scale wavelet denoising — noise foundation (see wavelet_denoise.h). */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <glib.h>

#include "core/siril.h"
#include "algos/Def_Wavelet.h"
#include "algos/sorting.h"
#include "algos/wavelet_denoise.h"

/* 1 / 0.6744897501960817 — scales the MAD to a Gaussian standard deviation. */
#define MAD_TO_SIGMA 1.482602218505602

/* Cache of per-scale propagation factors. The factors depend only on the
 * transform type and number of planes, so they are computed once and reused.
 * Indexed by transform type (TO_PAVE_* are small positive ints). */
typedef struct {
	gboolean valid;
	int nbr_plan;
	double e[WD_MAX_PLAN];
} wd_factor_cache;

#define WD_TYPE_SLOTS 9 /* covers TO_PAVE_LINEAR(1)..TO_MALLAT_BARLAUD(8) */
static wd_factor_cache wd_cache[WD_TYPE_SLOTS];
static GMutex wd_cache_mutex;

int wavelet_noise_factors(int type, int nbr_plan, double *e_out) {
	if (!e_out)
		return 1;
	if (type != TO_PAVE_LINEAR && type != TO_PAVE_BSPLINE)
		return 1;
	if (nbr_plan < 2 || nbr_plan > WD_MAX_PLAN)
		return 1;

	const int ndetail = nbr_plan - 1; /* detail scales 0..nbr_plan-2 */

	g_mutex_lock(&wd_cache_mutex);
	wd_factor_cache *c = &wd_cache[type];
	if (c->valid && c->nbr_plan == nbr_plan) {
		memcpy(e_out, c->e, ndetail * sizeof(double));
		g_mutex_unlock(&wd_cache_mutex);
		return 0;
	}
	g_mutex_unlock(&wd_cache_mutex);

	/* Transform a unit impulse through the actual kernels. The detail plane j
	 * then holds the equivalent detail filter's impulse response h_j, whose L2
	 * norm is exactly the std of that scale for unit-variance white noise
	 * (Var = sum of squared filter taps). The image is sized to satisfy
	 * wavelet_transform_data's minimum (Min >= 2^(nbr_plan+2)), which also
	 * leaves several times the filter support as margin so the compact-support
	 * response is captured exactly and border handling never intrudes. */
	const int N = 1 << (nbr_plan + 2);
	const size_t npix = (size_t) N * (size_t) N;
	float *impulse = calloc(npix, sizeof(float));
	if (!impulse)
		return 1;
	impulse[(size_t) (N / 2) * N + (N / 2)] = 1.0f;

	wave_transf_des wave = { 0 };
	if (wavelet_transform_data(impulse, N, N, &wave, type, nbr_plan)) {
		free(impulse);
		return 1;
	}

	for (int j = 0; j < ndetail; j++) {
		const float *plane = wave.Pave.Data + npix * (size_t) j;
		double ss = 0.0;
		for (size_t i = 0; i < npix; i++)
			ss += (double) plane[i] * (double) plane[i];
		e_out[j] = sqrt(ss);
	}

	wave_io_free(&wave);
	free(impulse);

	g_mutex_lock(&wd_cache_mutex);
	c->valid = TRUE;
	c->nbr_plan = nbr_plan;
	memcpy(c->e, e_out, ndetail * sizeof(double));
	g_mutex_unlock(&wd_cache_mutex);
	return 0;
}

double wavelet_mad_sigma_float(const float *band, size_t n) {
	if (!band || n == 0)
		return 0.0;
	float *buf = malloc(n * sizeof(float));
	if (!buf)
		return -1.0;

	memcpy(buf, band, n * sizeof(float));
	const float med = (float) quickmedian_float(buf, n);

	for (size_t i = 0; i < n; i++)
		buf[i] = fabsf(band[i] - med);
	const double mad = quickmedian_float(buf, n);

	free(buf);
	return mad * MAD_TO_SIGMA;
}

double wavelet_estimate_noise_float(const float *band0, size_t n, double e1) {
	if (e1 <= 0.0)
		return -1.0;
	const double s = wavelet_mad_sigma_float(band0, n);
	if (s < 0.0)
		return -1.0;
	return s / e1;
}

int wavelet_sigma_from_file(const char *filename, double *sigma_out) {
	if (!filename || !sigma_out)
		return 1;
	wave_transf_des wave = { 0 };
	if (wave_io_read((char *) filename, &wave))
		return 1;
	int ret = 1;
	if (wave.Nbr_Plan >= 2) {
		double e[WD_MAX_PLAN];
		if (!wavelet_noise_factors(wave.Type_Wave_Transform, wave.Nbr_Plan, e)) {
			const size_t npix = (size_t) wave.Nbr_Ligne * (size_t) wave.Nbr_Col;
			const double s = wavelet_estimate_noise_float(wave.Pave.Data, npix, e[0]);
			if (s >= 0.0) {
				*sigma_out = s;
				ret = 0;
			}
		}
	}
	wave_io_free(&wave);
	return ret;
}

void denoise_params_init(struct denoise_params *dp) {
	if (!dp)
		return;
	dp->enabled = FALSE;
	dp->method = WD_BISHRINK;
	dp->k = 3.0f;
	dp->sigma_source = WD_SIGMA_PROPAGATED;
	dp->soft = TRUE;
	dp->anscombe = FALSE;
	for (int i = 0; i < WD_MAX_PLAN; i++)
		dp->f[i] = 1.0f;
}

void anscombe_forward(float *data, size_t n, double scale) {
	if (!data || scale <= 0.0)
		return;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		double x = (double) data[i] * scale;
		if (x < 0.0)
			x = 0.0;
		data[i] = (float) (2.0 * sqrt(x + 0.375));
	}
}

void anscombe_inverse(float *data, size_t n, double scale) {
	if (!data || scale <= 0.0)
		return;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		const double y = (double) data[i];
		data[i] = (float) ((y * y * 0.25 - 0.375) / scale);
	}
}

/* Half-width of the local-variance window used by the bivariate shrinkage. */
#define WD_BISHRINK_RADIUS 3

static inline int clampi(int v, int lo, int hi) {
	return v < lo ? lo : (v > hi ? hi : v);
}

/* In-place soft/hard thresholding of one detail plane at threshold t. */
static void threshold_plane(float *plane, size_t n, float t, gboolean soft) {
	if (t <= 0.f)
		return; /* factor 0 (or no noise) -> leave this scale untouched */
	if (soft) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			const float w = plane[i];
			const float a = fabsf(w) - t;
			plane[i] = (a > 0.f) ? copysignf(a, w) : 0.f;
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			if (fabsf(plane[i]) <= t)
				plane[i] = 0.f;
		}
	}
}

/* Local mean of src^2 over a (2r+1)x(2r+1) window (separable, replicated
 * borders). dst receives the result; tmp is scratch of npix floats. */
static void box_mean_sq(const float *src, float *dst, float *tmp, int Nl,
		int Nc, int r) {
	const double norm = 1.0 / ((double) (2 * r + 1) * (2 * r + 1));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int y = 0; y < Nl; y++) {
		for (int x = 0; x < Nc; x++) {
			double s = 0.0;
			for (int dx = -r; dx <= r; dx++) {
				const float v = src[y * Nc + clampi(x + dx, 0, Nc - 1)];
				s += (double) v * v;
			}
			tmp[y * Nc + x] = (float) s;
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int y = 0; y < Nl; y++) {
		for (int x = 0; x < Nc; x++) {
			double s = 0.0;
			for (int dy = -r; dy <= r; dy++)
				s += tmp[clampi(y + dy, 0, Nl - 1) * Nc + x];
			dst[y * Nc + x] = (float) (s * norm);
		}
	}
}

/* Sendur-Selesnick bivariate (parent + local-variance) shrinkage of one detail
 * plane, in place. parent may be NULL (coarsest detail scale) -> reduces to a
 * locally-adaptive soft threshold. locms is the local mean of child^2; sigma_n
 * is the (effective) noise std of this scale. */
static void bishrink_plane(float *child, const float *parent,
		const float *locms, size_t n, double sigma_n) {
	const double sn2 = sigma_n * sigma_n;
	const double sqrt3 = 1.7320508075688772;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		const double c = child[i];
		const double p = parent ? (double) parent[i] : 0.0;
		const double R = sqrt(c * c + p * p);
		if (R <= 0.0) {
			child[i] = 0.f;
			continue;
		}
		/* marginal signal std from the local variance of the noisy band */
		const double sig = sqrt(fmax(0.0, (double) locms[i] - sn2));
		if (sig < 1e-12) {
			child[i] = 0.f; /* locally indistinguishable from noise */
			continue;
		}
		const double T = sqrt3 * sn2 / sig;
		double fac = (R - T) / R;
		if (fac < 0.0)
			fac = 0.0;
		child[i] = (float) (c * fac);
	}
}

int wavelet_denoise_planes(float *pave_data, int type, int nbr_plan, int Nl,
		int Nc, const struct denoise_params *dp) {
	if (!pave_data || !dp || !dp->enabled)
		return 0;
	if (nbr_plan < 2)
		return 0;

	const size_t npix = (size_t) Nl * (size_t) Nc;
	const int ndetail = nbr_plan - 1; /* planes 0..nbr_plan-2; last is residual */

	double e[WD_MAX_PLAN];
	if (wavelet_noise_factors(type, nbr_plan, e))
		return 1;

	/* Global noise from the finest detail plane (plane 0), unless each band is
	 * measured independently. */
	double sigma_g = 0.0;
	if (dp->sigma_source != WD_SIGMA_PER_BAND) {
		sigma_g = wavelet_estimate_noise_float(pave_data, npix, e[0]);
		if (sigma_g < 0.0)
			return 1;
	}

	const float k = (dp->k > 0.f) ? dp->k : 3.0f;

	/* Bivariate shrinkage needs scratch buffers and the parent (next coarser)
	 * band. Process finest -> coarsest so each plane's parent is still the
	 * original (noisy) band when used. */
	float *locms = NULL, *tmp = NULL;
	if (dp->method == WD_BISHRINK) {
		locms = malloc(npix * sizeof(float));
		tmp = malloc(npix * sizeof(float));
		if (!locms || !tmp) {
			free(locms);
			free(tmp);
			return 1;
		}
	}

	for (int j = 0; j < ndetail; j++) {
		float *plane = pave_data + npix * (size_t) j;
		const double sigma_j = (dp->sigma_source == WD_SIGMA_PER_BAND)
				? wavelet_mad_sigma_float(plane, npix)
				: sigma_g * e[j];

		if (dp->method == WD_BISHRINK) {
			/* k=3 is the neutral point: sigma_n = (k/3)*f_j*sigma_j, so the
			 * default uses the physical per-scale noise and k scales
			 * aggressiveness consistently with threshold mode. */
			const double sigma_n = (k / 3.0) * dp->f[j] * sigma_j;
			if (sigma_n <= 0.0)
				continue; /* factor 0 -> leave this scale untouched */
			const float *parent = (j + 1 < ndetail)
					? pave_data + npix * (size_t) (j + 1) : NULL;
			box_mean_sq(plane, locms, tmp, Nl, Nc, WD_BISHRINK_RADIUS);
			bishrink_plane(plane, parent, locms, npix, sigma_n);
		} else {
			const float t = (float) (k * dp->f[j] * sigma_j);
			threshold_plane(plane, npix, t, dp->soft);
		}
	}

	free(locms);
	free(tmp);
	return 0;
}
