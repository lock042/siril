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

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_version.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/siril_random.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "filters/mtf.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "core/gui_iface.h"
#include "opencv/opencv.h"
#include "background_extraction.h"

#define NPARAM_POLY4 15		// Number of parameters used with 4rd order
#define NPARAM_POLY3 10		// Number of parameters used with 3rd order
#define NPARAM_POLY2 6		// Number of parameters used with 2nd order
#define NPARAM_POLY1 3		// Number of parameters used with 1nd order

//C contains background function
#define C(i) (gsl_vector_get(c,(i)))

static double poly_4(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + (x * x) * C(3) + (x * y) * C(4)
			+ (y * y) * C(5) + (x * x) * x * C(6) + (x * x) * y * C(7)
			+ x * (y * y) * C(8)+ (y * y) * y * C(9) + (x * x) * (x * x) * C(10)
			+ (x * x) * (x * y) * C(11) + (x * x) * (y * y) * C(12)
			+ (x * y) * (y * y) * C(13) + (y * y) * (y * y) * C(14);

	return (value);
}

static double poly_3(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + (x * x) * C(3) + (y * x) * C(4)
			+ (y * y) * C(5) + (x * x) * x * C(6) + (x * x) * y * C(7) + x * (y * y) * C(8)
			+ (y * y) * y * C(9);

	return (value);
}

static double poly_2(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + C(3) * (x * x) + C(4) * (y * x)
			+ C(5) * (y * y);

	return (value);
}

static double poly_1(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y;

	return (value);
}

/* this function come from GSL 2.7 */
static double siril_gsl_vector_sum(const gsl_vector *v) {
	size_t i;
	double sum = 0;

	for (i = 0; i < v->size; ++i) {
		double vi = gsl_vector_get(v, i);
		sum += vi;
	}

	return sum;
}

static void cfachans_cleanup(fits **cfachans) {
	for (int i = 0 ; i < 4 ; i++) {
		clearfits(cfachans[i]);
		free(cfachans[i]);
	}
	free(cfachans);
}

static gboolean computeBackground_RBF(GSList *list, double *background, int channel, unsigned int width, unsigned int height, double smoothing, gchar **err, int threads) {
	/* Implementation of RBF interpolation with a thin-plate Kernel k(r) = r^2 * log(r)

	References:
	----------
	[1] G. B. Wright. Radial Basis Function Interpolation: Numerical and
		Analytical Developments. PhD thesis, University of Colorado, 2003.

	[2] SciPy version 1.1.0
		https://github.com/scipy/scipy/blob/v1.1.0/scipy/interpolate/rbf.py

	[3] S. Rippa. An algorithm for selecting a good value for the parameter
		c in radial basis function interpolation. Advances in Computational
		Mathematics, 11(2-3):193–210, 1999.

	[4] J. D. Martin and T. W. Simpson. Use of Kriging Models to Approximate
		Deterministic Computer Models. AIAA journal, 43(4):853–863, 2005.
	*/

	gsl_matrix *K;
	gsl_vector *f, *coef;
	GSList *l_i;
	guint n = g_slist_length(list);
	double mean = 0.0;
	double (*list_array)[3];
	siril_log_debug("RBF: %d samples, channel %d, w %u, h %u, smooth %f, %d threads\n", n, channel, width, height, smoothing, threads);

	K = gsl_matrix_calloc(n + 1, n + 1);
	f = gsl_vector_calloc(n + 1);
	coef = gsl_vector_calloc(n + 1);

	/* Scaling */
	int scaling_factor = 4;
	int width_scaled = round_to_int(width / scaling_factor);
	int height_scaled = round_to_int(height / scaling_factor);
	if (width_scaled <= 0 || height_scaled <= 0)
		return FALSE;
	double *background_scaled = calloc(width_scaled * height_scaled, sizeof(double));
	double *kernel_scaled = calloc(width_scaled * height_scaled, sizeof(double));
	double x_scaling = (double)height_scaled / (double)height;
	double y_scaling = (double)width_scaled / (double)width;

	/* Copy linked list into array */
	list_array = malloc(n * sizeof(double[3]));
	l_i = list;
	for (int i = 0; i < n; i++) {
		background_sample *sample_i = (background_sample*) l_i->data;
		list_array[i][0] = round(sample_i->position.x * x_scaling);
		list_array[i][1] = round(sample_i->position.y * y_scaling);
		list_array[i][2] = sample_i->median[channel];
		l_i = l_i->next;
	}

	/* Setup Kernel matrix K with K_ij = k(r_ij) and vector f with f_i = median(sample_i) */

	for (int i = 0; i < n; i++) {
		double x_i = list_array[i][0];
		double y_i = list_array[i][1];

		gsl_matrix_set(K, i, n, 1.0);
		gsl_matrix_set(K, n, i, 1.0);
		gsl_vector_set(f, i, list_array[i][2]);

		for (int j = 0; j < n; j++) {
			double x_j = list_array[j][0];
			double y_j = list_array[j][1];
			double distance = pow(x_i - x_j, 2) + pow(y_i - y_j, 2); // we can use r^2 directly in kernel calc
			double kernel = distance * log(distance) * 0.5;
			if (distance <= 1e-5) {
				kernel = 0.0;
			}

			gsl_matrix_set(K, i, j, kernel);
			mean += kernel / n;
		}
	}

	/* Smoothing */
	smoothing = 1e-4 * pow(10.0, (smoothing-0.5) * 3);

	for (int i = 0; i < n; i++) {
		gsl_matrix_set(K, i, i, smoothing * mean);
	}

	gsl_matrix_set(K, n, n, 0.0);
	gsl_vector_set(f, n, 0.0);

	/* Solve K*coef = f for coef */

	int s, status;
	gsl_permutation *p = gsl_permutation_alloc(n + 1);
	status = gsl_linalg_LU_decomp(K, p, &s);
	if (status) {
		siril_log_error("Error in RBF algorithm: %s\n", gsl_strerror(status));
		gsl_permutation_free(p);
		gsl_matrix_free(K);
		gsl_vector_free(f);
		gsl_vector_free(coef);
		free(background_scaled);
		free(kernel_scaled);
		free(list_array);
		return FALSE;
	}
	status = gsl_linalg_LU_solve(K, p, f, coef);
	if (status) {
		siril_log_error("Error in RBF algorithm: %s\n", gsl_strerror(status));
		gsl_permutation_free(p);
		gsl_matrix_free(K);
		gsl_vector_free(f);
		gsl_vector_free(coef);
		free(background_scaled);
		free(kernel_scaled);
		free(list_array);
		return FALSE;
	}

	/* precompute the kernel in the x>0, y>0 plane */
	for (int i = 0; i < height_scaled; i++) {
		for (int j = 0; j < width_scaled; j++) {
			if ((i == 0) && (j == 0)) {
				kernel_scaled[j + i * width_scaled] = 0.;
			} else {
				double distance = pow((double)j, 2.) + pow((double)i, 2.); // we can use r^2 directly in kernel calc
				kernel_scaled[j + i * width_scaled] = distance * log(distance) * 0.5;
			}
		}
	}

	/* Calculate background from coefficients coef */
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
	{
		gsl_vector *A = gsl_vector_calloc(n + 1);
#ifdef _OPENMP
#pragma omp for
#endif
		for (int i = 0; i < height_scaled; i++) {
			for (int j = 0; j < width_scaled; j++) {
				for (int k = 0; k < n; k++) {
					int deltax = abs(j - (int)list_array[k][0]);
					int deltay = abs(i - (int)list_array[k][1]);
					gsl_vector_set(A, k, kernel_scaled[deltax + deltay * width_scaled]);
				}

				gsl_vector_set(A, n, 1.0);
				gsl_vector_mul(A, coef);

				double pixel = siril_gsl_vector_sum(A);
				background_scaled[j + i * width_scaled] = pixel;
			}
		}

		gsl_vector_free(A);
	}


	cvResizeArray(background_scaled, background, height_scaled, width_scaled, height, width);

	gsl_permutation_free(p);
	gsl_matrix_free(K);
	gsl_vector_free(f);
	gsl_vector_free(coef);
	free(background_scaled);
	free(kernel_scaled);
	free(list_array);
	return TRUE;
}

static gboolean computeBackground_Polynom(GSList *list, double *background, int channel, unsigned int width, unsigned int height, poly_order order, gchar **err) {
	size_t k = 0;
	double chisq, pixel;
	gsl_matrix *J, *cov;
	gsl_vector *y, *w, *c;
	GSList *l;

	guint n = g_slist_length(list);

	int nbParam;
	switch (order) {
	case BACKGROUND_POLY_1:
		nbParam = NPARAM_POLY1;
		break;
	case BACKGROUND_POLY_2:
		nbParam = NPARAM_POLY2;
		break;
	case BACKGROUND_POLY_3:
		nbParam = NPARAM_POLY3;
		break;
	case BACKGROUND_POLY_4:
	default:
		nbParam = NPARAM_POLY4;
	}

	if (n < nbParam) {
		*err = siril_log_message(_("There are not enough background samples. "
				"The background to be extracted cannot be computed.\n"));
		return FALSE;
	}

	// J is the Jacobian
	// y contains data (pixel intensity)
	J = gsl_matrix_calloc(n, nbParam);
	y = gsl_vector_calloc(n);
	w = gsl_vector_calloc(n);
	c = gsl_vector_calloc(nbParam);
	cov = gsl_matrix_calloc(nbParam, nbParam);

	for (l = list; l; l = l->next) {
		background_sample *sample = (background_sample *) l->data;

		double col = sample->position.x;
		double row = sample->position.y;
		pixel = sample->median[channel];
		// here, it is a bit sketchy in the sense that if there is no
		// value to report in a box (because the threshold is too low
		// for example), then I just skip the initialization of J and y.
		// gsl automatically discards the unassigned values during the
		// minimization. I tested it with Matlab and it works fine. The
		// results agree.
		if (pixel < 0)
			continue;

		gsl_matrix_set(J, k, 0, 1.0);
		gsl_matrix_set(J, k, 1, col);
		gsl_matrix_set(J, k, 2, row);

		if (order != BACKGROUND_POLY_1) {
			gsl_matrix_set(J, k, 3, col * col);
			gsl_matrix_set(J, k, 4, col * row);
			gsl_matrix_set(J, k, 5, row * row);
		}

		if (order == BACKGROUND_POLY_3 || order == BACKGROUND_POLY_4) {
			gsl_matrix_set(J, k, 6, col * col * col);
			gsl_matrix_set(J, k, 7, col * col * row);
			gsl_matrix_set(J, k, 8, col * row * row);
			gsl_matrix_set(J, k, 9, row * row * row);
		}

		if (order == BACKGROUND_POLY_4) {
			gsl_matrix_set(J, k, 10, col * col * col * col);
			gsl_matrix_set(J, k, 11, col * col * col * row);
			gsl_matrix_set(J, k, 12, col * col * row * row);
			gsl_matrix_set(J, k, 13, col * row * row * row);
			gsl_matrix_set(J, k, 14, row * row * row * row);
		}

		gsl_vector_set(y, k, pixel);
		gsl_vector_set(w, k, 1.0);

		k++;
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nbParam);
	int status = gsl_multifit_wlinear(J, w, y, c, cov, &chisq, work);
	if (status != GSL_SUCCESS) {
		*err = siril_log_error("GSL multifit error: %s\n", gsl_strerror(status));
		gsl_matrix_free(J);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);
		return FALSE;
	}

	// Calculation of the background with the same dimension that the input matrix.
	for (unsigned int i = 0; i < height; i++) {
		for (unsigned int j = 0; j < width; j++) {
			switch (order) {
			case BACKGROUND_POLY_1:
				pixel = poly_1(c, (double) j, (double) i);
				break;
			case BACKGROUND_POLY_2:
				pixel = poly_2(c, (double) j, (double) i);
				break;
			case BACKGROUND_POLY_3:
				pixel = poly_3(c, (double) j, (double) i);
				break;
			default:
			case BACKGROUND_POLY_4:
				pixel = poly_4(c, (double) j, (double) i);
			}
			background[j + i * width] = pixel;
		}
	}

	/* free memory */
	gsl_multifit_linear_free(work);
	gsl_matrix_free(J);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

	return TRUE;
}

static background_sample *get_sample(float *buf, const int xx,
		const int yy, const int w, const int h) {
	size_t size = SAMPLE_SIZE * SAMPLE_SIZE;
	int radius = SAMPLE_SIZE / 2;
	background_sample *sample = malloc(sizeof(background_sample));
	if (!sample) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	double *data = malloc(size * sizeof(double));
	if (!data) {
		free(sample);
		PRINT_ALLOC_ERR;
		return NULL;
	}
	int n = 0;
	for (int y = yy - radius; y <= yy + radius; y ++) {
		for (int x = xx - radius; x <= xx + radius; x ++) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					data[n++] = (double)buf[x + y * w];
				}
			}
		}
	}
	if (n != size) {
		siril_log_debug("sample did not have the expected size (on border?)\n");
		free(sample);
		free(data);
		return NULL;
	}
	gsl_stats_minmax(&sample->min, &sample->max, data, 1, size);
	sample->mean = gsl_stats_mean(data, 1, size);
	sample->median[RLAYER] = quickmedian_double(data, size);
	sample->position.x = xx;
	sample->position.y = yy;
	sample->size = SAMPLE_SIZE;
	sample->valid = TRUE;

	free(data);
	return sample;
}

static double get_sample_median(const double *buf, const int xx,
		const int yy, const int w, const int h) {
	size_t size = SAMPLE_SIZE * SAMPLE_SIZE;
	int radius = SAMPLE_SIZE / 2;

	int n = 0;
	double *data = calloc(size, sizeof(double));
	if (!data) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}
	for (int y = yy - radius; y <= yy + radius; y ++) {
		for (int x = xx - radius; x <= xx + radius; x ++) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					data[n++] = buf[x + y * w];
				}
			}
		}
	}
	double median = quickmedian_double(data, size);

	free(data);
	return median;
}

static double get_background_mean(GSList *list, int num_channels) {
	/* Estimate background mean for all channels from samples */
	GSList *l;
	guint n = g_slist_length(list);
	double mean = 0.0;

	for (int channel = 0; channel < num_channels; channel++) {
		for (l = list; l; l = l->next) {
			background_sample *sample = (background_sample *) l->data;
			mean += sample->median[channel];
		}
	}

	return mean / n / num_channels;
}

static gboolean convert_fits_to_img(fits *fit, double *image, int channel, gboolean add_dither) {

	double invnorm = 1.0 / USHRT_MAX;
	const int height = fit->ry;
	const int width = fit->rx;
	if (fit->type == DATA_USHORT) {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				image[y * width + x] = fit->pdata[channel][(height - y - 1) * width + x] * invnorm;
				if (add_dither) {
					/* add dithering in order to avoid colour banding */
					image[y * width + x] += (siril_random_uint() & 0xFFFFF) * (1.0f / (1ULL << 36));
				}
			}
		}
	} else {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				image[y * width + x] = fit->fpdata[channel][(height - y - 1) * width + x];
				if (add_dither) {
					/* add dithering in order to avoid colour banding */
					image[y * width + x] += (siril_random_uint() & 0xFFFFF) * (1.0f / (1ULL << 36));
				}
			}
		}
	}
	return TRUE;
}

static float* convert_fits_to_luminance(fits *fit, threading_type threads) {
	g_assert(fit->type == DATA_USHORT || fit->type == DATA_FLOAT);
	const size_t n = fit->naxes[0] * fit->naxes[1];
	float invnorm = (float)(1.0 / USHRT_MAX);
	/* allocating memory to image */
	float *image = malloc(n * sizeof(float));
	if (!image) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	const int height = fit->ry;
	const int width = fit->rx;

#ifdef _OPENMP
	limit_threading(&threads, 200000/width, height);
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for (int y = 0; y < height; ++y) {
		size_t in_idx = (height - y - 1) * width;
		size_t out_idx = y * width;
		for (int x = 0; x < width; ++x) {
			if (fit->naxes[2] > 1) {
				float r, g, b;
				if (fit->type == DATA_USHORT) {
					r = fit->pdata[RLAYER][in_idx] * invnorm;
					g = fit->pdata[GLAYER][in_idx] * invnorm;
					b = fit->pdata[BLAYER][in_idx] * invnorm;
				} else {
					r = fit->fpdata[RLAYER][in_idx];
					g = fit->fpdata[GLAYER][in_idx];
					b = fit->fpdata[BLAYER][in_idx];
				}
				image[out_idx] = 0.3333f * r + 0.3333f * g + 0.3333f * b;
			} else {
				if (fit->type == DATA_USHORT) {
					image[out_idx] = fit->pdata[RLAYER][in_idx] * invnorm;
				} else if (fit->type == DATA_FLOAT) {
					image[out_idx] = fit->fpdata[RLAYER][in_idx];
				}
			}
			in_idx++;
			out_idx++;
		}
	}

	return image;
}

static void convert_img_to_fits(double *image, fits *fit, int channel) {
	const int height = fit->ry;
	const int width = fit->rx;
	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[channel];
		size_t out_idx = 0;
		for (int y = 0; y < height; ++y) {
			size_t in_idx = (height - y - 1) * width;
			for (int x = 0; x < width; ++x) {
				buf[out_idx] = round_to_WORD(image[in_idx] * USHRT_MAX);
				in_idx++;
				out_idx++;
			}
		}
	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[channel];
		size_t out_idx = 0;
		for (int y = 0; y < height; ++y) {
			size_t in_idx = (height - y - 1) * width;
			for (int x = 0; x < width; ++x) {
				buf[out_idx] = (float) image[in_idx];
				in_idx++;
				out_idx++;
			}
		}
	}
}

/* Remove samples whose centres are within SAMPLE_SIZE pixels of one already
 * accepted, keeping the first occurrence.  The caller must NOT use the
 * original list after this call — its nodes are consumed or freed. */
static GSList *deduplicate_background_samples(GSList *list) {
	GSList *out = NULL;
	for (GSList *l = list; l; l = l->next) {
		background_sample *s = (background_sample *)l->data;
		gboolean too_close = FALSE;
		for (GSList *o = out; o && !too_close; o = o->next) {
			background_sample *q = (background_sample *)o->data;
			int dx = (int)s->position.x - (int)q->position.x;
			int dy = (int)s->position.y - (int)q->position.y;
			if (dx * dx + dy * dy < SAMPLE_SIZE * SAMPLE_SIZE)
				too_close = TRUE;
		}
		if (too_close)
			free(s);
		else
			out = g_slist_prepend(out, s);
	}
	g_slist_free(list);
	return g_slist_reverse(out);
}

/* ---- Gradient-descent helpers ---- */

/* Compute median of the SAMPLE_SIZE x SAMPLE_SIZE patch centred on (x,y).
 * Uses a stack buffer since patch area is compile-time constant. */
static double gd_patch_median(const float *image, int x, int y, int w, int h) {
	double buf[SAMPLE_SIZE * SAMPLE_SIZE];
	int half = SAMPLE_SIZE / 2;
	int n = 0;
	for (int py = y - half; py <= y + half; py++) {
		for (int px = x - half; px <= x + half; px++) {
			if (px >= 0 && px < w && py >= 0 && py < h)
				buf[n++] = (double)image[py * w + px];
		}
	}
	return n > 0 ? quickmedian_double(buf, n) : 0.0;
}

/* Move (x,y) one pixel at a time toward the darkest 8-connected neighbour
 * (as judged by patch median) until a local minimum is found or the iteration
 * limit is reached.  The result is clamped so that get_sample() will succeed.
 * This optimization is based on the approach used by SetiAstro's AutoDBE, which
 * is (c) Franklin Marek and licensed as GPL v3.0-or-later.
 */
static void gradient_descent_to_dim_spot(const float *image, int *px, int *py, int w, int h) {
	int half = SAMPLE_SIZE / 2;
	int x = CLAMP(*px, half, w - half - 1);
	int y = CLAMP(*py, half, h - half - 1);

	for (int iter = 0; iter < 100; iter++) {
		double cur = gd_patch_median(image, x, y, w, h);
		int bx = x, by = y;
		double best = cur;

		for (int dy = -1; dy <= 1; dy++) {
			for (int dx = -1; dx <= 1; dx++) {
				if (dx == 0 && dy == 0) continue;
				int nx = x + dx, ny = y + dy;
				if (nx < half || nx >= w - half || ny < half || ny >= h - half)
					continue;
				double m = gd_patch_median(image, nx, ny, w, h);
				if (m < best) {
					best = m;
					bx = nx;
					by = ny;
				}
			}
		}
		if (bx == x && by == y)
			break;
		x = bx;
		y = by;
	}
	*px = x;
	*py = y;
}

/* ---- Random dark-area sample generation ---- */

/* Generate nb_samples interior samples from darker areas of each image
 * quadrant, plus fixed border/corner points.  Optionally applies
 * gradient descent to move each candidate to a local dark minimum.
 * This sample generation function is inspired by SetiAstro's AutoDBE
 * script which is (c) Franklin Marek and licensed as GPL-v3.0-or-later.
 */
static GSList *generate_samples_random(fits *fit, int nb_samples, int size,
		gboolean grad_descent, const char **error, threading_type threads,
		const rectangle *bbox) {
	int nx = fit->rx;
	int ny = fit->ry;
	GSList *list = NULL;
	int radius = size / 2;

	float *image = convert_fits_to_luminance(fit, threads);
	if (!image) {
		if (error) *error = "out of memory";
		return NULL;
	}

	/* Autostretch buffer */
	fits *tmp = NULL;
	new_fit_image_with_data(&tmp, fit->rx, fit->ry, 1, DATA_FLOAT, (void*) image);
	struct mtf_params params;
	find_unlinked_midtones_balance_default(tmp, &params);
	apply_unlinked_mtf_to_fits(tmp, tmp, &params);
	tmp->fdata = NULL;
	clearfits(tmp);
	free(tmp);

	/* Determine bounding box for sample placement */
	gboolean use_bbox = (bbox && bbox->w > 0 && bbox->h > 0);
	int sel_x0 = use_bbox ? bbox->x : 0;
	int sel_y0 = use_bbox ? bbox->y : 0;
	int sel_w  = use_bbox ? bbox->w : nx;
	int sel_h  = use_bbox ? bbox->h : ny;

	/* Clamp the selection to the image so all index arithmetic below stays
	 * within the nx*ny luminance buffer: a bbox that overruns the image
	 * (negative origin, or origin+extent past nx/ny) would otherwise drive
	 * the threshold-sampling read image[ry * nx + rx] out of bounds. */
	if (sel_x0 < 0) { sel_w += sel_x0; sel_x0 = 0; }
	if (sel_y0 < 0) { sel_h += sel_y0; sel_y0 = 0; }
	if (sel_x0 + sel_w > nx) sel_w = nx - sel_x0;
	if (sel_y0 + sel_h > ny) sel_h = ny - sel_y0;
	if (sel_w <= 0 || sel_h <= 0) {
		free(image);
		if (error) *error = "selection does not overlap the image";
		return NULL;
	}

	/* minimum coordinate so get_sample() succeeds */
	int margin = radius + 1;

/* Helper: accept a sample only if it contains no zero pixels */
#define ACCEPT_SAMPLE(s, list) \
	do { \
		if ((s) && (s)->min > 0.0) \
			(list) = g_slist_prepend((list), (s)); \
		else if (s) \
			free(s); \
	} while (0)

	/* ---- Border points ---- */
	/* 4 corners */
	int corners[4][2] = {
		{sel_x0 + margin,             sel_y0 + margin},
		{sel_x0 + sel_w - margin - 1, sel_y0 + margin},
		{sel_x0 + margin,             sel_y0 + sel_h - margin - 1},
		{sel_x0 + sel_w - margin - 1, sel_y0 + sel_h - margin - 1},
	};
	for (int i = 0; i < 4; i++) {
		int x = corners[i][0], y = corners[i][1];
		if (grad_descent)
			gradient_descent_to_dim_spot(image, &x, &y, nx, ny);
		if (!use_bbox || (x >= sel_x0 && x < sel_x0 + sel_w &&
				y >= sel_y0 && y < sel_y0 + sel_h)) {
			background_sample *s = get_sample(image, x, y, nx, ny);
			ACCEPT_SAMPLE(s, list);
		}
	}

	/* 5 evenly-spaced points along top and bottom edges */
	for (int k = 0; k < 5; k++) {
		int x = (sel_x0 + margin) + k * (sel_w - 2 * margin) / 4;
		/* top */
		{
			int xt = x, yt = sel_y0 + margin;
			if (grad_descent)
				gradient_descent_to_dim_spot(image, &xt, &yt, nx, ny);
			if (!use_bbox || (xt >= sel_x0 && xt < sel_x0 + sel_w &&
					yt >= sel_y0 && yt < sel_y0 + sel_h)) {
				background_sample *st = get_sample(image, xt, yt, nx, ny);
				ACCEPT_SAMPLE(st, list);
			}
		}
		/* bottom */
		{
			int xb = x, yb = sel_y0 + sel_h - margin - 1;
			if (grad_descent)
				gradient_descent_to_dim_spot(image, &xb, &yb, nx, ny);
			if (!use_bbox || (xb >= sel_x0 && xb < sel_x0 + sel_w &&
					yb >= sel_y0 && yb < sel_y0 + sel_h)) {
				background_sample *sb = get_sample(image, xb, yb, nx, ny);
				ACCEPT_SAMPLE(sb, list);
			}
		}
	}

	/* 5 evenly-spaced points along left and right edges */
	for (int k = 0; k < 5; k++) {
		int y = (sel_y0 + margin) + k * (sel_h - 2 * margin) / 4;
		/* left */
		{
			int xl = sel_x0 + margin, yl = y;
			if (grad_descent)
				gradient_descent_to_dim_spot(image, &xl, &yl, nx, ny);
			if (!use_bbox || (xl >= sel_x0 && xl < sel_x0 + sel_w &&
					yl >= sel_y0 && yl < sel_y0 + sel_h)) {
				background_sample *sl = get_sample(image, xl, yl, nx, ny);
				ACCEPT_SAMPLE(sl, list);
			}
		}
		/* right */
		{
			int xr = sel_x0 + sel_w - margin - 1, yr = y;
			if (grad_descent)
				gradient_descent_to_dim_spot(image, &xr, &yr, nx, ny);
			if (!use_bbox || (xr >= sel_x0 && xr < sel_x0 + sel_w &&
					yr >= sel_y0 && yr < sel_y0 + sel_h)) {
				background_sample *sr = get_sample(image, xr, yr, nx, ny);
				ACCEPT_SAMPLE(sr, list);
			}
		}
	}

	/* ---- Random interior points from 4 quadrants ---- */
	int pts_per_quad = MAX(1, nb_samples / 4);
	int half_nx = sel_x0 + sel_w / 2, half_ny = sel_y0 + sel_h / 2;

	int qxs[2] = {sel_x0, half_nx}, qxe[2] = {half_nx, sel_x0 + sel_w};
	int qys[2] = {sel_y0, half_ny}, qye[2] = {half_ny, sel_y0 + sel_h};

	for (int qi = 0; qi < 4; qi++) {
		int xs = qxs[qi % 2], xe = qxe[qi % 2];
		int ys = qys[qi / 2], ye = qye[qi / 2];
		int qw = xe - xs, qh = ye - ys;

		/* Estimate the 50th-percentile brightness by sampling ~10 000 pixels */
		int n_thresh = MIN(10000, qw * qh);
		float *tbuf = malloc(n_thresh * sizeof(float));
		if (!tbuf) continue;
		for (int i = 0; i < n_thresh; i++) {
			int rx = xs + (int)(siril_random_double() * qw);
			int ry = ys + (int)(siril_random_double() * qh);
			rx = CLAMP(rx, xs, xe - 1);
			ry = CLAMP(ry, ys, ye - 1);
			tbuf[i] = image[ry * nx + rx];
		}
		double threshold = histogram_median_float(tbuf, n_thresh, FALSE);
		free(tbuf);

		int min_x = MAX(xs, radius) + 1;
		int max_x = MIN(xe, nx - radius) - 1;
		int min_y = MAX(ys, radius) + 1;
		int max_y = MIN(ye, ny - radius) - 1;
		if (min_x >= max_x || min_y >= max_y) continue;

		/* Rejection sampling: pick random points below the threshold */
		int count = 0;
		int max_attempts = pts_per_quad * 200;
		for (int attempt = 0; attempt < max_attempts && count < pts_per_quad; attempt++) {
			int x = min_x + (int)(siril_random_double() * (max_x - min_x));
			int y = min_y + (int)(siril_random_double() * (max_y - min_y));
			x = CLAMP(x, min_x, max_x);
			y = CLAMP(y, min_y, max_y);
			if ((double)image[y * nx + x] < threshold) {
				if (grad_descent) {
					gradient_descent_to_dim_spot(image, &x, &y, nx, ny);
					if (use_bbox && (x < sel_x0 || x >= sel_x0 + sel_w ||
							y < sel_y0 || y >= sel_y0 + sel_h))
						continue;
				}
				background_sample *s = get_sample(image, x, y, nx, ny);
				if (s && s->min > 0.0) {
					list = g_slist_prepend(list, s);
					count++;
				} else if (s) {
					free(s);
				}
			}
		}
	}
#undef ACCEPT_SAMPLE

	list = deduplicate_background_samples(list);
	list = g_slist_reverse(list);
	free(image);
	if (!list && error)
		*error = "none of the samples matched the thresholds";
	return list;
}

GSList *generate_samples(fits *fit, int nb_per_line, double tolerance, int size, gboolean grad_descent, const char **error, threading_type threads, const rectangle *bbox) {
	int nx = fit->rx;
	int ny = fit->ry;
	size_t n = fit->naxes[0] * fit->naxes[1];
	GSList *list = NULL;
	float *image = convert_fits_to_luminance(fit, threads);	// upside down
	if (!image) {
		if (error)
			*error = "out of memory";
		return NULL;
	}
	float median = (float)histogram_median_float(image, n, threads);
	if (median <= 0.0f) {
		if (error)
			*error = "removing the gradient on negative images is not supported";
		free(image);
		return NULL;
	}

	gboolean use_bbox = (bbox && bbox->w > 0 && bbox->h > 0);
	int sel_x0 = use_bbox ? bbox->x : 0;
	int sel_y0 = use_bbox ? bbox->y : 0;
	int sel_w  = use_bbox ? bbox->w : nx;
	int sel_h  = use_bbox ? bbox->h : ny;

	int boxes_width = nb_per_line * size + 2;	// leave a margin of 1 px on the sides
	float spacing = (sel_w - boxes_width) / (float)(nb_per_line-1);
	int radius = size / 2;

	// Calculate nb_per_column using the same spacing as x-axis
	int nb_per_column = 1;
	while (nb_per_column * size + round_to_int((nb_per_column - 1) * spacing) < (sel_h - 2))
		nb_per_column++;
	nb_per_column--;
	if (nb_per_column == 0) {
		if (error)
			*error = "image is smaller than the sample size of the background extraction";
		free(image);
		return NULL;
	}

	// Calculate symmetric placement for y-axis within the region
	int total_grid_height = nb_per_column * size + (nb_per_column - 1) * round_to_int(spacing);
	int available_height = sel_h - 2;  // subtract margins
	int y_offset = (available_height - total_grid_height) / 2 + 1;  // +1 for the margin

	guint nb = nb_per_line * nb_per_column;
	float *mad = malloc(nb * sizeof(float));
	int k = 0;
	for (int i = 0; i < nb_per_line; i++) {
		for (int j = 0; j < nb_per_column; j++) {
			int x = sel_x0 + round_to_int(i * (spacing + size)) + radius + 1;
			int y = sel_y0 + y_offset + round_to_int(j * (spacing + size)) + radius;
			if (grad_descent) {
				gradient_descent_to_dim_spot(image, &x, &y, nx, ny);
				if (use_bbox && (x < sel_x0 || x >= sel_x0 + sel_w ||
						y < sel_y0 || y >= sel_y0 + sel_h))
					continue;
			}
			background_sample *sample = get_sample(image, x, y, nx, ny);
			if (sample) {
				mad[k++] = fabs(sample->median[RLAYER] - median);
				list = g_slist_prepend(list, sample);
			}
		}
	}
	/* compute mad */
	double mad0 = histogram_median_float(mad, k, TRUE);
	free(mad);
	double threshold = median + mad0 * tolerance;
	if (tolerance < 0.0)
		siril_log_debug("Background gradient: %d samples per line, no threshold\n", nb_per_line);
	else siril_log_debug("Background gradient: %d samples per line, threshold %f\n", nb_per_line, threshold);

	/* remove bad samples */
	GSList *l = list;
	while (l != NULL) {
		background_sample *sample = (background_sample*) l->data;
		/* Store next element's pointer before removing it */
		GSList *next = g_slist_next(l);
		if (sample->median[RLAYER] <= 0.0 || (tolerance > 0.0 && sample->median[RLAYER] >= threshold)) {
			free(sample);
			list = g_slist_delete_link(list, l);
		}
		l = next;
	}
	list = deduplicate_background_samples(list);
	list = g_slist_reverse(list);
	free(image);
	if (!list && error)
		*error = "none of the samples matched the thresholds";
	return list;
}

static GSList *update_median_samples(GSList *orig, fits *fit) {
	const int nx = fit->rx;
	const int ny = fit->ry;

	const size_t n = fit->naxes[0] * fit->naxes[1];
	double *channelData = malloc(n * sizeof(double));
	if (!channelData) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		convert_fits_to_img(fit, channelData, channel, FALSE);

		for (GSList *list = orig; list; list = list->next) {
			background_sample *sample = (background_sample*) list->data;
			sample->median[channel] = get_sample_median(channelData,
					sample->position.x, sample->position.y, nx, ny);
		}
	}

	free(channelData);
	return orig;
}

static void remove_gradient(double *img, const double *background, double background_mean, size_t ndata, background_correction type, threading_type threads) {
	switch (type) {
		default:
		case BACKGROUND_CORRECTION_SUBTRACT:
#ifdef _OPENMP
			limit_threading(&threads, 300000, ndata);
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
			for (size_t i = 0; i < ndata; i++) {
				img[i] -= background[i];
				img[i] += background_mean;
			}
			break;
		case BACKGROUND_CORRECTION_DIVIDE:
			{
				double mean = gsl_stats_mean(img, 1, ndata);
#ifdef _OPENMP
				limit_threading(&threads, 300000, ndata);
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
				for (size_t i = 0; i < ndata; i++) {
					img[i] /= background[i];
					img[i] *= mean;
				}
			}
	}
}

/************* PUBLIC FUNCTIONS *************/

static GMutex bgsamples_mutex = { 0 };

void sample_mutex_lock() {
	g_mutex_lock(&bgsamples_mutex);
}

void sample_mutex_unlock() {
	g_mutex_unlock(&bgsamples_mutex);
}

int get_background_sample_radius() {
       return SAMPLE_SIZE / 2;
}

void free_background_sample_list(GSList *list) {
	if (list == NULL) return;
	g_slist_free_full(list, free);
}

// This behaves as per add_background_sample but add a whole list of sample points
// at once, for use with the python interface to add points defined in python
// without incurring the cost of convert_fits_to_luminance for each one

GSList *add_background_samples(GSList *orig, fits *fit, GSList *pts) {
	GSList *list;
	int nx = fit->rx;
	int ny = fit->ry;
	float *image;
	image = convert_fits_to_luminance(fit, MULTI_THREADED);
	list = orig;
	for (GSList *iter = pts ; iter ; iter = iter->next) {
		point *pt = (point*) iter->data;
		background_sample *sample = get_sample(image, pt->x, pt->y, nx, ny);
		if (sample)
			list = g_slist_append(list, sample);
	}
	if (fit->naxes[2] > 1) {
		list = update_median_samples(list, fit);
	}
	free(image);

	return list;
}

GSList *add_background_sample(GSList *orig, fits *fit, point pt, gboolean grad_descent) {
	GSList *list;
	int nx = fit->rx;
	int ny = fit->ry;
	float *image;

	image = convert_fits_to_luminance(fit, MULTI_THREADED);

	list = orig;

	if (grad_descent) {
		int x = (int)pt.x, y = (int)pt.y;
		gradient_descent_to_dim_spot(image, &x, &y, nx, ny);
		pt.x = x;
		pt.y = y;
	}

	background_sample *sample = get_sample(image, pt.x, pt.y, nx, ny);
	list = g_slist_append(list, sample);

	if (fit->naxes[2] > 1) {
		list = update_median_samples(list, fit);
	}

	free(image);

	return list;
}

GSList *remove_background_sample(GSList *orig, fits *fit, point pt) {
	GSList *list;
	float *image;
	double min_radius = DBL_MAX;

	image = convert_fits_to_luminance(fit, MULTI_THREADED);

	/* search for the min radius vale */
	for (list = orig; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		double dx = pt.x - sample->position.x;
		double dy = pt.y - sample->position.y;
		double radius = sqrt(dx * dx + dy * dy);

		min_radius = min(min_radius, radius);
	}
	/* remove this value */
	for (list = orig; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		double dx = pt.x - sample->position.x;
		double dy = pt.y - sample->position.y;
		double radius = sqrt(dx * dx + dy * dy);

		if (radius == min_radius) {
			orig = g_slist_remove(orig, sample);
			g_free((background_sample *) sample);
			break;
		}
	}
	free(image);

	return orig;
}

/* generates samples and stores them in com.grad_samples */
int generate_background_samples(int nb_of_samples, double tolerance, gboolean randomize, gboolean grad_descent, const rectangle *override_bbox) {
	g_mutex_lock(&bgsamples_mutex);
	free_background_sample_list(com.grad_samples);
	const char *err = NULL;
	g_rw_lock_reader_lock(&gfit->rwlock);
	const rectangle *sel = (override_bbox && override_bbox->w > 0 && override_bbox->h > 0) ? override_bbox : NULL;
	if (sel)
		siril_log_debug("BGE: constraining sample placement to region (%d,%d %dx%d)\n",
			sel->x, sel->y, sel->w, sel->h);
	if (randomize) {
		com.grad_samples = generate_samples_random(gfit, nb_of_samples, SAMPLE_SIZE, grad_descent, &err, MULTI_THREADED, sel);
	} else {
		com.grad_samples = generate_samples(gfit, nb_of_samples, tolerance, SAMPLE_SIZE, grad_descent, &err, MULTI_THREADED, sel);
	}
	if (!com.grad_samples) {
		siril_log_error(_("Failed to generate background samples for image: %s\n"), err ? _(err) : _("unknown error"));
		g_rw_lock_reader_unlock(&gfit->rwlock);
		g_mutex_unlock(&bgsamples_mutex);
		return 1;
	}

	if (com.grad_samples && gfit->naxes[2] > 1) {
		/* If RGB we need to update all local median, not only the first one */
		com.grad_samples = update_median_samples(com.grad_samples, gfit);
	}
	g_rw_lock_reader_unlock(&gfit->rwlock);
	g_mutex_unlock(&bgsamples_mutex);
	return 0;
}

gboolean end_background(gpointer p);	// in gui/background_extraction.c

/* uses samples from com.grad_samples */
#if 0 /* dead code — no call sites; superseded by remove_gradient_image_hook */
gpointer remove_gradient_from_image(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	gchar *error = NULL;
	double *background = malloc(gfit->ry * gfit->rx * sizeof(double));

	if (!background) {
		PRINT_ALLOC_ERR;
		if (!com.script) {
			gui_iface.set_busy(FALSE);
		}
		return GINT_TO_POINTER(1);
	}

	const size_t n = gfit->naxes[0] * gfit->naxes[1];
	double *image = malloc(n * sizeof(double));
	if (!image) {
		free(background);
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(1);
	}

	/* Make sure to update local median. Useful if undo is pressed */
	g_mutex_lock(&bgsamples_mutex);
	update_median_samples(com.grad_samples, gfit);
	g_mutex_unlock(&bgsamples_mutex);

	double background_mean = get_background_mean(com.grad_samples, gfit->naxes[2]);
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	for (int channel = 0; channel < gfit->naxes[2]; channel++) {
		/* compute background */
		gboolean interpolation_worked = TRUE;
		if (args->interpolation_method == BACKGROUND_INTER_POLY) {
			interpolation_worked = computeBackground_Polynom(com.grad_samples, background, channel,
					gfit->rx, gfit->ry, args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(com.grad_samples, background, channel,
					gfit->rx, gfit->ry, args->smoothing, &error, args->threads);
		}

		if (!interpolation_worked) {
			free(image);
			free(background);
			gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Not enough samples."), error ? error : _("Insufficient samples"));
			if (!args->from_ui) {
				g_mutex_lock(&bgsamples_mutex);
				free_background_sample_list(com.grad_samples);
				com.grad_samples = NULL;
				g_mutex_unlock(&bgsamples_mutex);
			}
			free(args);
			notify_gfit_data_modified();
			siril_add_idle(end_background, NULL);
			return GINT_TO_POINTER(1);
		}
		/* remove background */
		const char *c_name;
		if (gfit->naxes[2] > 1)
			c_name = channel_number_to_name(channel);
		else
			c_name = _("monochrome");
		siril_log_message(_("Background extraction from %s channel.\n"), c_name);
		convert_fits_to_img(gfit, image, channel, args->dither);
		remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
		convert_img_to_fits(image, gfit, channel);

	}
	siril_log_message(_("Background with %s interpolation computed.\n"),
			(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	/* free memory */
	free(image);
	free(background);
	invalidate_stats_from_fit(gfit);
	if (!args->from_ui) {
		g_mutex_lock(&bgsamples_mutex);
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
		g_mutex_unlock(&bgsamples_mutex);
	}
	notify_gfit_data_modified();
	siril_add_idle(end_background, args);
	return GINT_TO_POINTER(0);
}
#endif /* dead code */

static GSList* rescale_sample_list_for_cfa(GSList *original_list, fits *fit) {
	GSList *new_list = NULL;
	GSList *current = original_list;
	int nx = fit->rx;
	int ny = fit->ry;
	int radius = SAMPLE_SIZE / 2;

	float *image = convert_fits_to_luminance(fit, MULTI_THREADED);	// upside down
	if (!image) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	// Traverse the original list in reverse order because g_slist_prepend is more efficient
	while (current != NULL) {
		background_sample *original_sample = (background_sample *)current->data;
		int x = min(nx - 1 - radius, max(radius, round_to_int((float) original_sample->position.x / 2.f)));
		int y = min(ny - 1 - radius, max(radius, round_to_int((float) original_sample->position.y / 2.f)));
		// Generate a new sample at half the original coordinates
		background_sample *new_sample = get_sample(image, x, y, nx, ny);

		// Add the new sample to the new list
		new_list = g_slist_prepend(new_list, new_sample);

		current = current->next;
	}

	free(image);

	// Reverse the list to maintain the original order
	new_list = g_slist_reverse(new_list);

	return new_list;
}

/* uses samples from com.grad_samples */
#if 0 /* dead code — no call sites; superseded by remove_gradient_image_hook */
gpointer remove_gradient_from_cfa_image(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	sensor_pattern pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
		siril_log_error(_("Error: unsupported CFA pattern for this operation.\n"));
		return GINT_TO_POINTER(1);
	}
	gchar *error = NULL;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);

	// Allocate an array of fits* and the fits objects themselves. These will hold the subchannels
	fits** cfachans = calloc(4, sizeof(fits*));
	for (int i = 0 ; i < 4 ; i++) {
		cfachans[i] = calloc(1, sizeof(fits));
	}

	// Split the FITS into the 4 subchannels
	int ret = 1;
	if (gfit->type == DATA_USHORT) {
		ret = split_cfa_ushort(gfit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	else if (gfit->type == DATA_FLOAT) {
		ret = split_cfa_float(gfit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	if (ret) {
		siril_log_error(_("Error splitting into CFA subchannels, aborting...\n"));
		cfachans_cleanup(cfachans);
		return GINT_TO_POINTER(1);
	}

	// Check subchannel medians are OK
	for (int i = 0 ; i < 4 ; i++) {
		imstats* stat = statistics(NULL, -1, cfachans[i], 0, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_error(_("Error: statistics computation failed.\n"));
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}
		float median = (float) stat->median;
		free_stats(stat);

		if (median <= 0.0f) {
			siril_log_error(_("Subchannel with negative median detected: removing the gradient on negative images is not supported\n"));
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}
	}

	for (int i = 0; i < 4; i++) {
		fits *subchannel = cfachans[i];

		GSList *samples = rescale_sample_list_for_cfa(com.grad_samples, subchannel);

		if (!samples) {
			siril_log_error(_("Failed to adapt background samples for CFA image\n"));
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}

		double *background = (double*)malloc(subchannel->naxes[0] * subchannel->naxes[1] * sizeof(double));
		if (!background) {
			PRINT_ALLOC_ERR;
			siril_log_error(_("Out of memory - aborting"));
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}

		const size_t n = subchannel->naxes[0] * subchannel->naxes[1];
		double *image = malloc(n * sizeof(double));
		if (!image) {
			free(background);
			free_background_sample_list(samples);
			PRINT_ALLOC_ERR;
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}

		double background_mean = get_background_mean(samples, 1);
		/* compute background */
		gboolean interpolation_worked = TRUE;
		if (args->interpolation_method == BACKGROUND_INTER_POLY) {
			interpolation_worked = computeBackground_Polynom(samples, background, 0,
					subchannel->rx, subchannel->ry, args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(samples, background, 0,
					subchannel->rx, subchannel->ry, args->smoothing, &error, args->threads);
		}

		if (!interpolation_worked) {
			free(image);
			free(background);
			gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Not enough samples."), error);
			if (!args->from_ui) {
				g_mutex_lock(&bgsamples_mutex);
				free_background_sample_list(com.grad_samples);
				com.grad_samples = NULL;
				g_mutex_unlock(&bgsamples_mutex);
			}
			cfachans_cleanup(cfachans);
			free(args);
			notify_gfit_data_modified();
			siril_add_idle(end_background, NULL);
			return GINT_TO_POINTER(1);
		}
		/* remove background */
		convert_fits_to_img(subchannel, image, 0, args->dither);
		remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
		convert_img_to_fits(image, subchannel, 0);
		free(image);
		free(background);
		free_background_sample_list(samples);

	}
	fits *out = merge_cfa(cfachans[0], cfachans[1], cfachans[2], cfachans[3], pattern);
	fits_swap_image_data(out, gfit); // Efficiently move the merged pixeldata from out to gfit
	clearfits(out);
	free(out);
	siril_log_message(_("Background with %s interpolation computed for CFA image.\n"),
			(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	/* free memory */
	cfachans_cleanup(cfachans);
	invalidate_stats_from_fit(gfit);
	if (!args->from_ui) {
		g_mutex_lock(&bgsamples_mutex);
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
		g_mutex_unlock(&bgsamples_mutex);
	}
	notify_gfit_data_modified();
	siril_add_idle(end_background, args);
	return GINT_TO_POINTER(0);
}
#endif /* dead code */

void free_background_data(void *p) {
	struct background_data *args = (struct background_data *)p;
	if (!args) return;
	free(args->seqEntry);
	free(args);
}

/* Single-image hook for generic_image_worker.
 * Handles both CFA and non-CFA paths based on args->is_cfa.
 * Uses fit parameter instead of gfit; no notify/idle calls. */
gchar *remove_gradient_log_hook(gpointer p, log_hook_detail detail) {
	struct background_data *args = (struct background_data *)p;
	if (args->method == BACKGROUND_METHOD_AUTO) {
		const gchar *model = args->autograd.simplified ? _("simplified") : _("multiscale");
		return g_strdup_printf(_("Automatic gradient removal (%s)"), model);
	}
	const gchar *interp = (args->interpolation_method == BACKGROUND_INTER_POLY) ? _("polynomial") : _("RBF");
	if (args->is_cfa)
		return g_strdup_printf(_("Background extraction (%s, CFA)"), interp);
	return g_strdup_printf(_("Background extraction (%s)"), interp);
}

/* Automatic (sample-free) background model
 *
 * The background is fitted on every pixel that survives an iterative robust rejection of
 * structures (stars, nebulae). Two models: a multiscale smooth surface with
 * structure protection (default) or a stiff low-degree polynomial (simplified).
 * The estimation path works on planar float buffers of size width*height and
 * uses RawTherapee's fast separable Gaussian (rt/gauss.cc) for every low-pass. */

#define AG_HIGH_K 2.0f
#define AG_LOW_K 4.0f
#define AG_N_ITER 20
#define AG_PASSES 3

/* Gaussian sigma that variance-matches `passes` box blurs of radius r: this is
 * how the old repeated-box low-pass is mapped onto the single RT Gaussian, so
 * the effective smoothing scale (and hence the model behaviour) is preserved. */
static double ag_sigma(int r, int passes) {
	if (r < 1) return 0.0;
	return sqrt((double)passes * r * (r + 1) / 3.0);
}

/* One separable box blur of radius r along a single line of `n` float samples
 * with the given stride, via a running sum. Edge pixels are replicated (numpy
 * "edge" padding), so the divisor is a constant 2r+1. */
static void ag_box1d_line(const float *in, float *out, int n, int stride, int r) {
	const int w = 2 * r + 1;
	const float inv = 1.0f / w;
	float s = 0.0f;
	for (int k = -r; k <= r; k++) {
		int idx = k < 0 ? 0 : (k >= n ? n - 1 : k);
		s += in[(size_t)idx * stride];
	}
	out[0] = s * inv;
	for (int x = 1; x < n; x++) {
		int leave = x - 1 - r; leave = leave < 0 ? 0 : (leave >= n ? n - 1 : leave);
		int enter = x + r;     enter = enter < 0 ? 0 : (enter >= n ? n - 1 : enter);
		s += in[(size_t)enter * stride] - in[(size_t)leave * stride];
		out[(size_t)x * stride] = s * inv;
	}
}

/* Separable low-pass of a planar float buffer approximating a Gaussian of the
 * given sigma with AG_PASSES running-sum box blurs (variance-matched radius).
 * O(n) per pass and independent of radius; `in` is preserved unless out == in. */
static void ag_gauss(const float *in, float *out, int w, int h, double sigma, int threads) {
	const size_t n = (size_t)w * h;
	/* invert ag_sigma(): sigma^2 = passes * r(r+1)/3  ->  r */
	int r = (int)lround((-1.0 + sqrt(1.0 + 12.0 * sigma * sigma / AG_PASSES)) / 2.0);
	if (sigma < 0.25 || r < 1) {
		if (out != in) memcpy(out, in, n * sizeof(float));
		return;
	}
	float *tmp = malloc(n * sizeof(float));
	if (!tmp) {
		if (out != in) memcpy(out, in, n * sizeof(float));
		return;
	}
	memcpy(out, in, n * sizeof(float));
	for (int p = 0; p < AG_PASSES; p++) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(h > 64)
#endif
		for (int y = 0; y < h; y++)                       /* horizontal: out -> tmp */
			ag_box1d_line(out + (size_t)y * w, tmp + (size_t)y * w, w, 1, r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(w > 64)
#endif
		for (int x = 0; x < w; x++)                       /* vertical: tmp -> out */
			ag_box1d_line(tmp + x, out + x, h, w, r);
	}
	free(tmp);
}

/* Robust (median, sigma) of a float buffer via Siril's quickmedian_float
 * (quickselect) applied twice: once for the median, once for the median of the
 * absolute deviations (MAD); sigma = 1.4826*MAD. `x` is preserved, `buf`
 * (n floats) is scratch. */
static void ag_mad_sigma(const float *x, size_t n, float *buf, float *med_out,
		float *sigma_out) {
	if (n == 0) { *med_out = 0.0f; *sigma_out = 1e-12f; return; }
	memcpy(buf, x, n * sizeof(float));
	float med = (float)quickmedian_float(buf, n);
	double mad = siril_stats_float_mad(x, n, med, MULTI_THREADED, buf);
	*med_out = med;
	*sigma_out = (float)(1.4826 * mad + 1e-12);
}

/* Robust smooth fill that bridges arbitrarily large rejected holes.
 * Harmonic-style inpainting: repeatedly low-pass then restore the kept pixels.
 * `out` receives the filled+smoothed model; `f0` and `f1` are scratch buffers
 * of size w*h. */
static void ag_inpaint_lowpass(const float *img, const gboolean *mask, int w,
		int h, int radius, int n_fill, int passes, float *out, float *f0,
		float *f1, int threads) {
	const size_t n = (size_t)w * h;
	const double sigma = ag_sigma(radius, passes);
	size_t nkept = 0;
	double sum = 0.0;
	for (size_t i = 0; i < n; i++)
		if (mask[i]) { nkept++; sum += img[i]; }
	if (nkept == 0) {
		memset(out, 0, n * sizeof(float));
		return;
	}
	if (nkept == n) {                        /* no holes: plain low-pass */
		ag_gauss(img, out, w, h, sigma, threads);
		return;
	}
	const float known_mean = (float)(sum / nkept);
	float *filled = f0;
	for (size_t i = 0; i < n; i++)
		filled[i] = mask[i] ? img[i] : known_mean;
	for (int it = 0; it < n_fill; it++) {
		ag_gauss(filled, f1, w, h, sigma, threads);   /* sm = lowpass(filled) */
		for (size_t i = 0; i < n; i++)
			filled[i] = mask[i] ? img[i] : f1[i];
	}
	ag_gauss(filled, out, w, h, sigma, threads);
}

/* Spatially-coherent mask of extended bright structures to protect.
 * `struct_mask` receives TRUE where a structure is to be excluded from the
 * fit. `det` and `grown` are scratch buffers of size w*h. */
static void ag_structure_mask(const float *residual, int w, int h,
		int model_radius, float protect_threshold, float protect_amount,
		gboolean *struct_mask, float *det, float *grown, int threads) {
	const size_t n = (size_t)w * h;
	gboolean any = FALSE;
	for (size_t i = 0; i < n; i++) {
		det[i] = (residual[i] > protect_threshold) ? 1.0f : 0.0f;
		if (det[i] > 0.0f) any = TRUE;
	}
	if (!any) {
		memset(struct_mask, 0, n * sizeof(gboolean));
		return;
	}
	int grow_r = (int)lround(model_radius * (0.5 + protect_amount));
	if (grow_r < 1) grow_r = 1;
	ag_gauss(det, grown, w, h, ag_sigma(grow_r, 2), threads);
	const float cutoff = (1.0f - protect_amount) * 0.5f + 1e-3f;
	for (size_t i = 0; i < n; i++)
		struct_mask[i] = grown[i] > cutoff;
}

/* Number of terms x^i*y^j with i+j <= degree. */
static int ag_poly_nterms(int degree) {
	return (degree + 1) * (degree + 2) / 2;
}

/* Least-squares fit of the polynomial basis over masked pixels via the normal
 * equations (small T*T system), then evaluate the model over all pixels.
 * `xn`/`yn` are the normalized [-1,1] coordinate axes. Returns FALSE on a
 * singular system. */
static gboolean ag_poly_fit(const float *ch, const gboolean *mask, int w, int h,
		int degree, const double *xn, const double *yn, float *model) {
	const int T = ag_poly_nterms(degree);
	int *ei = malloc(T * sizeof(int)), *ej = malloc(T * sizeof(int));
	if (!ei || !ej) { free(ei); free(ej); return FALSE; }
	int t = 0;
	for (int i = 0; i <= degree; i++)
		for (int j = 0; j <= degree - i; j++) { ei[t] = i; ej[t] = j; t++; }

	gsl_matrix *AtA = gsl_matrix_calloc(T, T);
	gsl_vector *Atb = gsl_vector_calloc(T);
	gsl_vector *coef = gsl_vector_calloc(T);
	double *tv = malloc(T * sizeof(double));
	if (!AtA || !Atb || !coef || !tv) {
		free(ei); free(ej); free(tv);
		if (AtA) gsl_matrix_free(AtA);
		if (Atb) gsl_vector_free(Atb);
		if (coef) gsl_vector_free(coef);
		return FALSE;
	}

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			size_t idx = (size_t)y * w + x;
			if (!mask[idx]) continue;
			for (int k = 0; k < T; k++)
				tv[k] = pow(xn[x], ei[k]) * pow(yn[y], ej[k]);
			const double b = ch[idx];
			for (int a = 0; a < T; a++) {
				double *row = gsl_matrix_ptr(AtA, a, 0);
				const double ta = tv[a];
				for (int c = 0; c < T; c++)
					row[c] += ta * tv[c];
				*gsl_vector_ptr(Atb, a) += ta * b;
			}
		}
	}

	/* Tiny ridge term keeps the normal equations well-conditioned. */
	for (int a = 0; a < T; a++)
		*gsl_matrix_ptr(AtA, a, a) += 1e-9;

	gsl_error_handler_t *old = gsl_set_error_handler_off();
	int status = gsl_linalg_cholesky_decomp1(AtA);
	if (status == GSL_SUCCESS)
		status = gsl_linalg_cholesky_solve(AtA, Atb, coef);
	gsl_set_error_handler(old);
	if (status != GSL_SUCCESS) {
		gsl_matrix_free(AtA); gsl_vector_free(Atb); gsl_vector_free(coef);
		free(ei); free(ej); free(tv);
		return FALSE;
	}

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			double v = 0.0;
			for (int k = 0; k < T; k++)
				v += gsl_vector_get(coef, k) * pow(xn[x], ei[k]) * pow(yn[y], ej[k]);
			model[(size_t)y * w + x] = (float)v;
		}
	}

	gsl_matrix_free(AtA); gsl_vector_free(Atb); gsl_vector_free(coef);
	free(ei); free(ej); free(tv);
	return TRUE;
}

/* Iterative robust background model for one channel of size w*h. `model`
 * receives the fitted background. Returns FALSE on allocation/fit failure. */
static gboolean ag_estimate_background(const float *ch, int w, int h, int radius,
		const struct autograd_data *p, float *model, int threads,
		double prog_base, double prog_span, const char *prog_msg) {
	const size_t n = (size_t)w * h;
	if (radius < 1) radius = 1;

	gboolean ok = TRUE;
	gboolean *keep = malloc(n * sizeof(gboolean));
	gboolean *new_keep = malloc(n * sizeof(gboolean));
	gboolean *smask = malloc(n * sizeof(gboolean));
	float *residual = malloc(n * sizeof(float));
	float *ref = malloc(n * sizeof(float));
	float *sortbuf = malloc(n * sizeof(float));
	/* scratch for inpaint (f0,f1) / structure mask (det,grown) / residual-med */
	float *s0 = malloc(n * sizeof(float));
	float *s1 = malloc(n * sizeof(float));
	float *s2 = malloc(n * sizeof(float));
	double *xn = NULL, *yn = NULL;
	if (!keep || !new_keep || !smask || !residual || !ref || !sortbuf ||
			!s0 || !s1 || !s2) { ok = FALSE; goto cleanup; }

	if (p->simplified) {
		xn = malloc(w * sizeof(double));
		yn = malloc(h * sizeof(double));
		if (!xn || !yn) { ok = FALSE; goto cleanup; }
		for (int x = 0; x < w; x++) xn[x] = (double)x / (w > 1 ? w - 1 : 1) * 2.0 - 1.0;
		for (int y = 0; y < h; y++) yn[y] = (double)y / (h > 1 ? h - 1 : 1) * 2.0 - 1.0;
	}

#define AG_FIT(mask)                                                        \
	do {                                                                    \
		if (p->simplified) {                                                \
			if (!ag_poly_fit(ch, (mask), w, h, p->degree, xn, yn, model)) { \
				ok = FALSE; goto cleanup;                                   \
			}                                                               \
		} else {                                                            \
			ag_inpaint_lowpass(ch, (mask), w, h, radius, 10, 2, model,      \
					s0, s1, threads);                                       \
		}                                                                   \
	} while (0)

	for (size_t i = 0; i < n; i++) keep[i] = TRUE;
	AG_FIT(keep);
	size_t prev = n;
	size_t min_keep = (size_t)(0.02 * n);
	if (min_keep < 16) min_keep = 16;

	for (int it = 0; it < AG_N_ITER; it++) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(n > 100000)
#endif
		for (size_t i = 0; i < n; i++)
			residual[i] = ch[i] - model[i];
		size_t nref = 0;
		for (size_t i = 0; i < n; i++)
			if (keep[i]) ref[nref++] = residual[i];
		if (nref == 0)
			for (size_t i = 0; i < n; i++) ref[i] = residual[i], nref = n;
		float med, sigma;
		ag_mad_sigma(ref, nref, sortbuf, &med, &sigma);
		const float hi = med + AG_HIGH_K * sigma;
		const float lo = med - AG_LOW_K * sigma;
		for (size_t i = 0; i < n; i++)
			new_keep[i] = (residual[i] <= hi) && (residual[i] >= lo);
		if (p->protect) {
			/* structure_mask is fed residual - med */
			for (size_t i = 0; i < n; i++) s2[i] = residual[i] - med;
			ag_structure_mask(s2, w, h, radius, (float)p->protect_threshold,
					(float)p->protect_amount, smask, s0, s1, threads);
			for (size_t i = 0; i < n; i++)
				if (smask[i]) new_keep[i] = FALSE;
		}
		size_t kept = 0;
		for (size_t i = 0; i < n; i++) if (new_keep[i]) kept++;
		if (kept < min_keep) {                 /* never empty the fit set (rare) */
			memcpy(sortbuf, residual, n * sizeof(float));
			quicksort_f(sortbuf, n);
			size_t rank = min_keep >= n ? n - 1 : min_keep;
			float thr = sortbuf[rank];
			kept = 0;
			for (size_t i = 0; i < n; i++) {
				new_keep[i] = residual[i] <= thr;
				if (new_keep[i]) kept++;
			}
		}
		AG_FIT(new_keep);
		double change = (double)(kept > prev ? kept - prev : prev - kept) / n;
		memcpy(keep, new_keep, n * sizeof(gboolean));
		prev = kept;
		/* report progress within this channel's slice; converging early just
		 * makes the bar reach the slice end sooner (skipped on sequences) */
		if (prog_span > 0.0)
			gui_iface.set_progress(prog_base + prog_span * (it + 1) / (double)AG_N_ITER, prog_msg);
		if (it > 0 && change < 1e-4)
			break;
	}
#undef AG_FIT

	if (p->smoothness > 0.0) {
		int sr = (int)lround(radius * p->smoothness);
		if (sr < 1) sr = 1;
		ag_gauss(model, s0, w, h, ag_sigma(sr, AG_PASSES), threads);
		memcpy(model, s0, n * sizeof(float));
	}

cleanup:
	free(keep); free(new_keep); free(smask); free(residual); free(ref);
	free(sortbuf); free(s0); free(s1); free(s2);
	free(xn); free(yn);
	return ok;
}

/* Area-average downsample by integer factor f (block mean) into `out` of size
 * (w/f)*(h/f). Returns the small dimensions via *sw,*sh. */
static void ag_downsample(const double *img, int w, int h, int f, float *out,
		int *sw, int *sh) {
	int ow = w / f, oh = h / f;
	*sw = ow; *sh = oh;
	for (int y = 0; y < oh; y++) {
		for (int x = 0; x < ow; x++) {
			double s = 0.0;
			for (int by = 0; by < f; by++)
				for (int bx = 0; bx < f; bx++)
					s += img[(size_t)(y * f + by) * w + (x * f + bx)];
			out[(size_t)y * ow + x] = (float)(s / (f * f));
		}
	}
}

/* Bilinear resize of `img` (sw*sh) up to (ow*oh) into `out`. */
static void ag_resize_bilinear(const float *img, int sw, int sh, int ow, int oh,
		float *out, int threads) {
	if (sw == ow && sh == oh) {
		memcpy(out, img, (size_t)ow * oh * sizeof(float));
		return;
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(oh > 64)
#endif
	for (int y = 0; y < oh; y++) {
		double fy = (oh > 1) ? (double)y * (sh - 1) / (oh - 1) : 0.0;
		int y0 = (int)fy; int y1 = y0 + 1 < sh ? y0 + 1 : sh - 1;
		double wy = fy - y0;
		for (int x = 0; x < ow; x++) {
			double fx = (ow > 1) ? (double)x * (sw - 1) / (ow - 1) : 0.0;
			int x0 = (int)fx; int x1 = x0 + 1 < sw ? x0 + 1 : sw - 1;
			double wx = fx - x0;
			double Ia = img[(size_t)y0 * sw + x0], Ib = img[(size_t)y0 * sw + x1];
			double Ic = img[(size_t)y1 * sw + x0], Id = img[(size_t)y1 * sw + x1];
			double top = Ia * (1 - wx) + Ib * wx;
			double bot = Ic * (1 - wx) + Id * wx;
			out[(size_t)y * ow + x] = (float)(top * (1 - wy) + bot * wy);
		}
	}
}

/* Correct one channel in place. `image` (w*h, double for FITS I/O precision) is
 * both input and output; the heavy estimation runs on float buffers. Returns
 * FALSE on failure. */
static gboolean auto_gradient_channel(double *image, int w, int h,
		const struct autograd_data *p, background_correction mode, int threads,
		double prog_base, double prog_span, const char *prog_msg) {
	int f = p->downsample;
	if (f < 1 || f > w || f > h) f = 1;
	int sw, sh;
	float *small = malloc((size_t)(w / f + 1) * (h / f + 1) * sizeof(float));
	if (!small) { PRINT_ALLOC_ERR; return FALSE; }
	ag_downsample(image, w, h, f, small, &sw, &sh);
	if (sw < 2 || sh < 2) { free(small); return FALSE; }

	int mind = sw < sh ? sw : sh;
	int radius = (int)lround(p->scale / 100.0 * mind);
	if (radius < 1) radius = 1;

	float *model_small = malloc((size_t)sw * sh * sizeof(float));
	float *bg = malloc((size_t)w * h * sizeof(float));
	float *sortbuf = malloc((size_t)w * h * sizeof(float));
	if (!model_small || !bg || !sortbuf) {
		free(small); free(model_small); free(bg); free(sortbuf);
		PRINT_ALLOC_ERR; return FALSE;
	}

	/* keep a little of the slice for the final resize/correction step */
	gboolean ok = ag_estimate_background(small, sw, sh, radius, p, model_small,
			threads, prog_base, prog_span * 0.9, prog_msg);
	if (ok) {
		ag_resize_bilinear(model_small, sw, sh, w, h, bg, threads);
		const size_t n = (size_t)w * h;
		memcpy(sortbuf, bg, n * sizeof(float));         /* preserve bg for the correction */
		double level = quickmedian_float(sortbuf, n);
		if (mode == BACKGROUND_CORRECTION_DIVIDE) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(n > 100000)
#endif
			for (size_t i = 0; i < n; i++)
				image[i] = image[i] / (bg[i] > 1e-6f ? bg[i] : 1e-6f) * level;
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(n > 100000)
#endif
			for (size_t i = 0; i < n; i++)
				image[i] = image[i] - bg[i] + level;
		}
	}
	free(small); free(model_small); free(bg); free(sortbuf);
	return ok;
}

/* Automatic (sample-free) path of the background extraction image hook.
 * Processes each channel of `fit` independently; CFA images are treated as a
 * plain mono frame (no per-subchannel split). */
static int auto_gradient_image_hook(struct background_data *args, fits *fit, int threads) {
	const size_t n = fit->naxes[0] * fit->naxes[1];
	const int nchan = fit->naxes[2];
	/* On a sequence the per-frame progress bar is driven by the sequence worker,
	 * so only report fine-grained progress when processing a single image. */
	const gboolean report = (args->seq == NULL);
	double *image = malloc(n * sizeof(double));
	if (!image) { PRINT_ALLOC_ERR; return 1; }
	for (int channel = 0; channel < nchan; channel++) {
		const char *c_name = nchan > 1 ? channel_number_to_name(channel) : _("monochrome");
		gchar *prog_msg = g_strdup_printf(_("Automatic gradient removal (%s)"), c_name);
		/* each channel owns an equal slice of the [0,1] progress bar */
		double base = (double)channel / nchan;
		double span = 1.0 / nchan;
		if (report)
			gui_iface.set_progress(base, prog_msg);
		convert_fits_to_img(fit, image, channel, args->dither);
		if (!auto_gradient_channel(image, fit->rx, fit->ry, &args->autograd,
				args->correction, threads,
				report ? base : PROGRESS_NONE, report ? span : 0.0, prog_msg)) {
			siril_log_error(_("Automatic gradient removal failed.\n"));
			g_free(prog_msg);
			free(image);
			return 1;
		}
		siril_log_message(_("Automatic gradient removed from %s channel.\n"), c_name);
		convert_img_to_fits(image, fit, channel);
		g_free(prog_msg);
	}
	free(image);
	invalidate_stats_from_fit(fit);
	if (report)
		gui_iface.set_progress(PROGRESS_DONE, _("Automatic gradient removal done"));
	return 0;
}

int remove_gradient_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct background_data *args = (struct background_data *)gargs->user;

	if (args->method == BACKGROUND_METHOD_AUTO)
		return auto_gradient_image_hook(args, fit, threads);

	gchar *error = NULL;

	if (args->is_cfa) {
		sensor_pattern pattern = get_validated_cfa_pattern(fit, FALSE, FALSE);
		if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
			siril_log_error(_("Error: unsupported CFA pattern for this operation.\n"));
			return 1;
		}

		fits **cfachans = calloc(4, sizeof(fits *));
		for (int i = 0 ; i < 4 ; i++)
			cfachans[i] = calloc(1, sizeof(fits));

		int ret = 1;
		if (fit->type == DATA_USHORT)
			ret = split_cfa_ushort(fit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
		else if (fit->type == DATA_FLOAT)
			ret = split_cfa_float(fit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
		if (ret) {
			siril_log_error(_("Error splitting into CFA subchannels, aborting...\n"));
			cfachans_cleanup(cfachans);
			return 1;
		}

		for (int i = 0 ; i < 4 ; i++) {
			imstats *stat = statistics(NULL, -1, cfachans[i], 0, NULL, STATS_BASIC, MULTI_THREADED);
			if (!stat) {
				siril_log_error(_("Error: statistics computation failed.\n"));
				cfachans_cleanup(cfachans);
				return 1;
			}
			float median = (float)stat->median;
			free_stats(stat);
			if (median <= 0.0f) {
				siril_log_error(_("Subchannel with negative median detected: removing the gradient on negative images is not supported\n"));
				cfachans_cleanup(cfachans);
				return 1;
			}
		}

		for (int i = 0 ; i < 4 ; i++) {
			fits *subchannel = cfachans[i];
			GSList *samples = rescale_sample_list_for_cfa(com.grad_samples, subchannel);
			if (!samples) {
				siril_log_error(_("Failed to adapt background samples for CFA image\n"));
				cfachans_cleanup(cfachans);
				return 1;
			}
			double *background = malloc(subchannel->naxes[0] * subchannel->naxes[1] * sizeof(double));
			if (!background) {
				PRINT_ALLOC_ERR;
				free_background_sample_list(samples);
				cfachans_cleanup(cfachans);
				return 1;
			}
			const size_t n = subchannel->naxes[0] * subchannel->naxes[1];
			double *image = malloc(n * sizeof(double));
			if (!image) {
				free(background);
				free_background_sample_list(samples);
				cfachans_cleanup(cfachans);
				PRINT_ALLOC_ERR;
				return 1;
			}
			double background_mean = get_background_mean(samples, 1);
			gboolean interpolation_worked = TRUE;
			if (args->interpolation_method == BACKGROUND_INTER_POLY) {
				interpolation_worked = computeBackground_Polynom(samples, background, 0,
						subchannel->rx, subchannel->ry, args->degree, &error);
			} else {
				interpolation_worked = computeBackground_RBF(samples, background, 0,
						subchannel->rx, subchannel->ry, args->smoothing, &error, args->threads);
			}
			if (!interpolation_worked) {
				free(image);
				free(background);
				free_background_sample_list(samples);
				if (!args->from_ui) {
					g_mutex_lock(&bgsamples_mutex);
					free_background_sample_list(com.grad_samples);
					com.grad_samples = NULL;
					g_mutex_unlock(&bgsamples_mutex);
				}
				cfachans_cleanup(cfachans);
				gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Not enough samples."), error);
				return 1;
			}
			convert_fits_to_img(subchannel, image, 0, args->dither);
			remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
			convert_img_to_fits(image, subchannel, 0);
			free(image);
			free(background);
			free_background_sample_list(samples);
		}
		fits *out = merge_cfa(cfachans[0], cfachans[1], cfachans[2], cfachans[3], pattern);
		fits_swap_image_data(out, fit);
		clearfits(out);
		free(out);
		siril_log_message(_("Background with %s interpolation computed for CFA image.\n"),
				(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
		cfachans_cleanup(cfachans);
		invalidate_stats_from_fit(fit);
		if (!args->from_ui) {
			g_mutex_lock(&bgsamples_mutex);
			free_background_sample_list(com.grad_samples);
			com.grad_samples = NULL;
			g_mutex_unlock(&bgsamples_mutex);
		}
	} else {
		double *background = malloc(fit->ry * fit->rx * sizeof(double));
		if (!background) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		const size_t n = fit->naxes[0] * fit->naxes[1];
		double *image = malloc(n * sizeof(double));
		if (!image) {
			free(background);
			PRINT_ALLOC_ERR;
			return 1;
		}
		g_mutex_lock(&bgsamples_mutex);
		update_median_samples(com.grad_samples, fit);
		g_mutex_unlock(&bgsamples_mutex);
		double background_mean = get_background_mean(com.grad_samples, fit->naxes[2]);
		for (int channel = 0; channel < fit->naxes[2]; channel++) {
			gboolean interpolation_worked = TRUE;
			if (args->interpolation_method == BACKGROUND_INTER_POLY) {
				interpolation_worked = computeBackground_Polynom(com.grad_samples, background, channel,
						fit->rx, fit->ry, args->degree, &error);
			} else {
				interpolation_worked = computeBackground_RBF(com.grad_samples, background, channel,
						fit->rx, fit->ry, args->smoothing, &error, args->threads);
			}
			if (!interpolation_worked) {
				free(image);
				free(background);
				gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Not enough samples."), error ? error : _("Insufficient samples"));
				if (!args->from_ui) {
					g_mutex_lock(&bgsamples_mutex);
					free_background_sample_list(com.grad_samples);
					com.grad_samples = NULL;
					g_mutex_unlock(&bgsamples_mutex);
				}
				return 1;
			}
			const char *c_name = fit->naxes[2] > 1 ? channel_number_to_name(channel) : _("monochrome");
			siril_log_message(_("Background extraction from %s channel.\n"), c_name);
			convert_fits_to_img(fit, image, channel, args->dither);
			remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
			convert_img_to_fits(image, fit, channel);
		}
		siril_log_message(_("Background with %s interpolation computed.\n"),
				(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
		free(image);
		free(background);
		invalidate_stats_from_fit(fit);
		if (!args->from_ui) {
			g_mutex_lock(&bgsamples_mutex);
			free_background_sample_list(com.grad_samples);
			com.grad_samples = NULL;
			g_mutex_unlock(&bgsamples_mutex);
		}
	}
	return 0;
}

/** Apply for sequence **/

/* Compute the inset bbox obtained by excluding a border strip from all sides.
 * pixel_scale allows halving the pixel border for CFA subchannels (pass 0.5). */
static rectangle compute_border_bbox(int rx, int ry, double border_value, gboolean border_is_percent, double pixel_scale) {
	int bx, by;
	if (border_is_percent) {
		bx = (int)(rx * border_value / 100.0 + 0.5);
		by = (int)(ry * border_value / 100.0 + 0.5);
	} else {
		bx = (int)(border_value * pixel_scale + 0.5);
		by = (int)(border_value * pixel_scale + 0.5);
	}
	if (bx >= rx / 2) bx = rx / 2 - 1;
	if (by >= ry / 2) by = ry / 2 - 1;
	rectangle bbox = { bx, by, rx - 2 * bx, ry - 2 * by };
	return bbox;
}

static int background_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct background_data *b_args = (struct background_data*) args->user;

	if (b_args->method == BACKGROUND_METHOD_AUTO)
		return auto_gradient_image_hook(b_args, fit, threads);

	double *background = (double*)malloc(fit->naxes[0] * fit->naxes[1] * sizeof(double));
	if (!background) {
		PRINT_ALLOC_ERR;
		siril_log_error(_("Out of memory - aborting"));
		return 1;
	}

	rectangle border_bbox;
	const rectangle *bbox = NULL;
	if (b_args->border_value > 0.0) {
		border_bbox = compute_border_bbox(fit->rx, fit->ry, b_args->border_value, b_args->border_is_percent, 1.0);
		bbox = &border_bbox;
	}

	const char *err;
	GSList *samples;
	if (b_args->randomize) {
		samples = generate_samples_random(fit, b_args->nb_of_samples, SAMPLE_SIZE, b_args->grad_descent, &err, (threading_type)threads, bbox);
	} else {
		samples = generate_samples(fit, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE, b_args->grad_descent, &err, (threading_type)threads, bbox);
	}
	if (!samples) {
		siril_log_error(_("Failed to generate background samples for image %d: %s\n"), i, _(err));
		free(background);
		return 1;
	}

	/* If RGB we need to update all local median, not only the first one */
	if (fit->naxes[2] > 1) {
		samples = update_median_samples(samples, fit);
	}

	const size_t n = fit->naxes[0] * fit->naxes[1];
	double *image = malloc(n * sizeof(double));
	if (!image) {
		free(background);
		free_background_sample_list(samples);
		PRINT_ALLOC_ERR;
		return 1;
	}

	double background_mean = get_background_mean(samples, fit->naxes[2]);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		/* compute background */
		gboolean interpolation_worked = TRUE;
		gchar *error = NULL;
		if (b_args->interpolation_method == BACKGROUND_INTER_POLY){
			interpolation_worked = computeBackground_Polynom(samples, background, channel, fit->rx, fit->ry, b_args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(samples, background, channel, fit->rx, fit->ry, b_args->smoothing, &error, threads);
		}

		if (!interpolation_worked) {
			if (error) {
				siril_log_message(error);
			}
			free(image);
			free(background);
			free_background_sample_list(samples);
			return 1;
		}
		/* remove background */
		convert_fits_to_img(fit, image, channel, b_args->dither);
		remove_gradient(image, background, background_mean, fit->naxes[0] * fit->naxes[1], b_args->correction, (threading_type)threads);
		convert_img_to_fits(image, fit, channel);
	}
	/* free memory */
	free(image);
	free(background);
	free_background_sample_list(samples);

	return 0;
}

/*******************************************************************************************
 * In the case of Bayer CFA images, we can consider the image as 4 interleaved sub-images: *
 * we split the image into these 4 sub-images and apply background removal to each         *
 * independently. We allow for different gradients in each sub-image as this may           *
 * realistically model a frequency dependence of the observed gradient, for example        *
 * moonlight.                                                                              *
 * On completion of gradient removal from each sub-image the sub-images are reassembled    *
 * ready for drizzle.                                                                      *
 ******************************************************************************************/

static int bgcfa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct background_data *b_args = (struct background_data*) args->user;
	if (b_args->interpolation_method == BACKGROUND_INTER_RBF) {
		siril_log_warning(_("Warning: RBF background removal is not recommended for CFA images. Only linear background removal is recommended.\n"));
	} else if (b_args->degree > 1) {
		siril_log_warning(_("Warning: polynomial background removal order > 1 is not recommended for CFA images. Only linear background removal is recommended.\n"));
	}
	sensor_pattern pattern = get_validated_cfa_pattern(fit, FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
		siril_log_error(_("Error: unsupported CFA pattern for this operation.\n"));
		return 1;
	}

	// Allocate an array of fits* and the fits objects themselves. These will hold the subchannels
	fits** cfachans = calloc(4, sizeof(fits*));
	for (int i = 0 ; i < 4 ; i++) {
		cfachans[i] = calloc(1, sizeof(fits));
	}

	// Split the FITS into the 4 subchannels
	int ret = 1;
	if (fit->type == DATA_USHORT) {
		ret = split_cfa_ushort(fit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	else if (fit->type == DATA_FLOAT) {
		ret = split_cfa_float(fit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	if (ret) {
		siril_log_error(_("Error splitting into CFA subchannels, aborting...\n"));
		cfachans_cleanup(cfachans);
		return 1;
	}

	// Clear the original fit image data, we will put the result back into it later
	free(fit->data);
	free(fit->fdata);
	fit->data = fit->pdata[0] = fit->pdata[1] = fit->pdata[2] = NULL;
	fit->fdata = fit->fpdata[0] = fit->fpdata[1] = fit->fpdata[2] = NULL;


	// Carry out background removal on each subchannel in turn
	for (int i = 0 ; i < 4 ; i++) {
		fits *subchannel = cfachans[i];
		double *background = (double*)malloc(subchannel->naxes[0] * subchannel->naxes[1] * sizeof(double));
		if (!background) {
			PRINT_ALLOC_ERR;
			siril_log_error(_("Out of memory - aborting"));
			cfachans_cleanup(cfachans);
			return 1;
		}

		rectangle sub_border_bbox;
		const rectangle *sub_bbox = NULL;
		if (b_args->border_value > 0.0) {
			/* CFA subchannels are half the original image size in each dimension;
			 * pixel borders must be scaled accordingly */
			sub_border_bbox = compute_border_bbox(subchannel->rx, subchannel->ry,
				b_args->border_value, b_args->border_is_percent,
				b_args->border_is_percent ? 1.0 : 0.5);
			sub_bbox = &sub_border_bbox;
		}

		const char *err;
		GSList *samples;
		if (b_args->randomize) {
			samples = generate_samples_random(subchannel, b_args->nb_of_samples, SAMPLE_SIZE, b_args->grad_descent, &err, (threading_type)threads, sub_bbox);
		} else {
			samples = generate_samples(subchannel, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE, b_args->grad_descent, &err, (threading_type)threads, sub_bbox);
		}
		if (!samples) {
			siril_log_error(_("Failed to generate background samples for image %d: %s\n"), i, _(err));
			free(background);
			cfachans_cleanup(cfachans);
			return 1;
		}

		const size_t n = subchannel->naxes[0] * subchannel->naxes[1];
		double *image = malloc(n * sizeof(double));
		if (!image) {
			free(background);
			free_background_sample_list(samples);
			PRINT_ALLOC_ERR;
			cfachans_cleanup(cfachans);
			return 1;
		}

		double background_mean = get_background_mean(samples, subchannel->naxes[2]);
		/* compute background */
		gboolean interpolation_worked = TRUE;
		gchar *error = NULL;
		if (b_args->interpolation_method == BACKGROUND_INTER_POLY){
			interpolation_worked = computeBackground_Polynom(samples, background, 0, subchannel->rx, subchannel->ry, b_args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(samples, background, 0, subchannel->rx, subchannel->ry, b_args->smoothing, &error, threads);
		}

		if (!interpolation_worked) {
			if (error) {
				siril_log_message(error);
			}
			free(image);
			free(background);
			free_background_sample_list(samples);
			cfachans_cleanup(cfachans);
			return 1;
		}
		/* remove background */
		convert_fits_to_img(subchannel, image, 0, b_args->dither);
		remove_gradient(image, background, background_mean, subchannel->naxes[0] * subchannel->naxes[1], b_args->correction, (threading_type)threads);
		convert_img_to_fits(image, subchannel, 0);
		/* free memory */
		free(image);
		free(background);
		free_background_sample_list(samples);
	}
	fits *out = merge_cfa(cfachans[0], cfachans[1], cfachans[2], cfachans[3], pattern);
	fits_swap_image_data(out, fit); // Efficiently move the merged pixeldata from out to fit
	clearfits(out);
	free(out);
	for (int i = 0 ; i < 4 ; i++) {
		free(cfachans[i]); // no need to use cfachans_cleanup here as the FITS have already
				// been cleared in merge_cfa()
	}
	free(cfachans);
	return 0;
}

static int background_mem_limits_hook(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail;

	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	unsigned int required = MB_per_image;
	if (limit > 0) {
		/* allocations:
		 * generate_samples convert_fits_to_luminance allocates         rx * ry * sizeof(float)
		 * generate_samples allocates a buffer for MAD computation      rx * ry * sizeof(double)
		 * both are freed at generate_samples exit
		 * for color images:
		 *	update_median_for_rgb_samples allocates for median      rx * ry * sizeof(double)
		 * freed at update_median_for_rgb_samples exit
		 * remove_gradient_from_image allocates the background image to rx * ry * sizeof(double)
		 * remove_gradient_from_image allocates the image            to rx * ry * sizeof(double)
		 *
		 * so at maximum, ignoring the samples, we need 2 times the double channel size.
		 *
		 * Despite creating the CFA subsequences, peak memory use does not increase when handling
		 * a CFA sequence as we free the pixel data in fit before starting to process each of the
		 * subchannels (essentially we are doing the same operation on 4 FITS each 1/4 the size)
		 *
		 */
		uint64_t double_channel_size = args->seq->rx * args->seq->ry * sizeof(double);
		unsigned int double_channel_size_MB = double_channel_size / BYTES_IN_A_MB;
		if (double_channel_size_MB == 0)
			double_channel_size_MB = 1;
		required = MB_per_image + double_channel_size_MB * 2;
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_image;
		} else limit = thread_limit;

	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_error(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"), args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_log_debug("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int bg_extract_finalize_hook(struct generic_seq_args *args) {
	struct background_data *data = (struct background_data *) args->user;
	int retval = seq_finalize_hook(args);
	free(data->seqEntry);
	free(data);
	return retval;
}

void apply_background_extraction_to_sequence(struct background_data *background_args) {
	fits metadata = { 0 };
	struct generic_seq_args *args = create_default_seqargs(background_args->seq);
	if (seq_read_frame_metadata(background_args->seq, sequence_find_refimage(background_args->seq), &metadata)) {
		siril_log_error(_("Error reading reference metadata.\n"));
		free(background_args->seqEntry);
		free(background_args);
		free_generic_seq_args(args, TRUE);
		return;
	}
	sensor_pattern pattern = get_cfa_pattern_index_from_string(metadata.keywords.bayer_pattern);
	background_args->is_cfa = background_args->seq->nb_layers == 1 && pattern >= BAYER_FILTER_MIN && pattern <= BAYER_FILTER_MAX;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = background_args->seq->selnum;
	args->compute_mem_limits_hook = background_mem_limits_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = bg_extract_finalize_hook;
	/* the automatic model treats CFA frames as a plain mono image (no split) */
	args->image_hook = (background_args->is_cfa && background_args->method != BACKGROUND_METHOD_AUTO)
			? bgcfa_image_hook : background_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Background Extraction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(background_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = background_args;

	background_args->fit = NULL;	// not used here

	if(!start_in_new_thread(generic_sequence_worker, args)) {
		free(background_args->seqEntry);
		free(background_args);
		free_generic_seq_args(args, TRUE);
	}
}

/**** getters ***/

gboolean background_sample_is_valid(background_sample *sample) {
	return sample->valid;
}

gdouble background_sample_get_size(background_sample *sample) {
	return sample->size;
}

point background_sample_get_position(background_sample *sample) {
	return sample->position;
}
