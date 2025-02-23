/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "gui/message_dialog.h"
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
	siril_debug_print("RBF: %d samples, channel %d, w %u, h %u, smooth %f, %d threads\n", n, channel, width, height, smoothing, threads);

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
		siril_log_color_message("Error in RBF algorithm: %s\n", "red", gsl_strerror(status));
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
		siril_log_color_message("Error in RBF algorithm: %s\n", "red", gsl_strerror(status));
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
		*err = siril_log_message("GSL multifit error: %s\n", gsl_strerror(status));
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
		siril_debug_print("sample did not have the expected size (on border?)\n");
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

GSList *generate_samples(fits *fit, int nb_per_line, double tolerance, int size, const char **error, threading_type threads) {
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

	int boxes_width = nb_per_line * size + 2;	// leave a margin of 1 px on the sides
	float spacing = (nx - boxes_width) / (float)(nb_per_line-1);
	int radius = size / 2;
	// solving iteratively n for: ny - 2 = n * size + (n-1) * spacing;
	int nb_per_column = 1;
	while (nb_per_column * size + round_to_int((nb_per_column - 1) * spacing) <= (ny - 2))
		nb_per_column++;
	nb_per_column--;
	if (nb_per_column == 0) {
		if (error)
			*error = "image is smaller than the sample size of the background extraction";
		free(image);
		return NULL;
	}

	guint nb = nb_per_line * nb_per_column;
	float *mad = malloc(nb * sizeof(float));
	int k = 0;

	for (int i = 0; i < nb_per_line; i++) {
		for (int j = 0; j < nb_per_column; j++) {
			int x = round_to_int(i * (spacing + size)) + radius + 1;
			int y = round_to_int(j * (spacing + size)) + radius + 1;
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
	siril_debug_print("Background gradient: %d samples per line, threshold %f\n", nb_per_line, threshold);

	/* remove bad samples */
	GSList *l = list;
	while (l != NULL) {
		background_sample *sample = (background_sample*) l->data;
		/* Store next element's pointer before removing it */
		GSList *next = g_slist_next(l);
		if (sample->median[RLAYER] <= 0.0 || sample->median[RLAYER] >= threshold) {
			free(sample);
			list = g_slist_delete_link(list, l);
		}
		l = next;
	}

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

GSList *add_background_sample(GSList *orig, fits *fit, point pt) {
	GSList *list;
	int nx = fit->rx;
	int ny = fit->ry;
	float *image;

	image = convert_fits_to_luminance(fit, MULTI_THREADED);

	list = orig;

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
int generate_background_samples(int nb_of_samples, double tolerance) {
	free_background_sample_list(com.grad_samples);
	const char *err;
	com.grad_samples = generate_samples(&gfit, nb_of_samples, tolerance, SAMPLE_SIZE, &err, MULTI_THREADED);
	if (!com.grad_samples) {
		siril_log_color_message(_("Failed to generate background samples for image: %s\n"), "red", _(err));
		return 1;
	}

	if (com.grad_samples && gfit.naxes[2] > 1) {
		/* If RGB we need to update all local median, not only the first one */
		com.grad_samples = update_median_samples(com.grad_samples, &gfit);
	}
	return 0;
}

gboolean end_background(gpointer p);	// in gui/background_extraction.c

/* uses samples from com.grad_samples */
gpointer remove_gradient_from_image(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	gchar *error = NULL;
	double *background = malloc(gfit.ry * gfit.rx * sizeof(double));

	if (!background) {
		PRINT_ALLOC_ERR;
		if (!com.script) {
			set_cursor_waiting(FALSE);
		}
		return GINT_TO_POINTER(1);
	}

	const size_t n = gfit.naxes[0] * gfit.naxes[1];
	double *image = malloc(n * sizeof(double));
	if (!image) {
		free(background);
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(1);
	}

	/* Make sure to update local median. Useful if undo is pressed */
	update_median_samples(com.grad_samples, &gfit);

	double background_mean = get_background_mean(com.grad_samples, gfit.naxes[2]);
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	for (int channel = 0; channel < gfit.naxes[2]; channel++) {
		/* compute background */
		gboolean interpolation_worked = TRUE;
		if (args->interpolation_method == BACKGROUND_INTER_POLY) {
			interpolation_worked = computeBackground_Polynom(com.grad_samples, background, channel,
					gfit.rx, gfit.ry, args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(com.grad_samples, background, channel,
					gfit.rx, gfit.ry, args->smoothing, &error, args->threads);
		}

		if (!interpolation_worked) {
			free(image);
			free(background);
			queue_error_message_dialog(_("Not enough samples."), error ? error : _("Insufficient samples"));
			if (!args->from_ui) {
				free_background_sample_list(com.grad_samples);
				com.grad_samples = NULL;
			}
			free(args);
			siril_add_idle(end_background, NULL);
			return GINT_TO_POINTER(1);
		}
		/* remove background */
		const char *c_name;
		if (gfit.naxes[2] > 1)
			c_name = channel_number_to_name(channel);
		else
			c_name = _("monochrome");
		siril_log_message(_("Background extraction from %s channel.\n"), c_name);
		convert_fits_to_img(&gfit, image, channel, args->dither);
		remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
		convert_img_to_fits(image, &gfit, channel);

	}
	siril_log_message(_("Background with %s interpolation computed.\n"),
			(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	/* free memory */
	free(image);
	free(background);
	invalidate_stats_from_fit(&gfit);
	if (!args->from_ui) {
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
	}
	siril_add_idle(end_background, args);
	return GINT_TO_POINTER(0);
}

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
gpointer remove_gradient_from_cfa_image(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	sensor_pattern pattern;
	if (!strncmp(gfit.keywords.bayer_pattern, "RGGB", 4)) {
		pattern = BAYER_FILTER_RGGB;
	} else if (!strncmp(gfit.keywords.bayer_pattern, "BGGR", 4)) {
		pattern = BAYER_FILTER_BGGR;
	} else if (!strncmp(gfit.keywords.bayer_pattern, "GBRG", 4)) {
		pattern = BAYER_FILTER_GBRG;
	} else if (!strncmp(gfit.keywords.bayer_pattern, "GRBG", 4)) {
		pattern = BAYER_FILTER_GBRG;
	} else {
		siril_log_color_message(_("Error: unsupported CFA pattern for this operation.\n"), "red");
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
	if (gfit.type == DATA_USHORT) {
		ret = split_cfa_ushort(&gfit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	else if (gfit.type == DATA_FLOAT) {
		ret = split_cfa_float(&gfit, cfachans[0], cfachans[1], cfachans[2], cfachans[3]);
	}
	if (ret) {
		siril_log_color_message(_("Error splitting into CFA subcannels, aborting...\n"), "red");
		cfachans_cleanup(cfachans);
		return GINT_TO_POINTER(1);
	}

	// Check subchannel medians are OK
	for (int i = 0 ; i < 4 ; i++) {
		size_t n = cfachans[i]->naxes[0] * cfachans[i]->naxes[1];
		float median = (float)histogram_median_float(cfachans[i], n, MULTI_THREADED);
		if (median <= 0.0f) {
			siril_log_color_message(_("Subchannel with negative median detected: removing the gradient on negative images is not supported\n"), "red");
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}
	}

	for (int i = 0; i < 4; i++) {
		fits *subchannel = cfachans[i];

		GSList *samples = rescale_sample_list_for_cfa(com.grad_samples, subchannel);

		if (!samples) {
			siril_log_color_message(_("Failed to adapt background samples for CFA image\n"), "red");
			cfachans_cleanup(cfachans);
			return GINT_TO_POINTER(1);
		}

		double *background = (double*)malloc(subchannel->naxes[0] * subchannel->naxes[1] * sizeof(double));
		if (!background) {
			PRINT_ALLOC_ERR;
			siril_log_message(_("Out of memory - aborting"));
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
			queue_error_message_dialog(_("Not enough samples."), error);
			if (!args->from_ui) {
				free_background_sample_list(com.grad_samples);
				com.grad_samples = NULL;
			}
			cfachans_cleanup(cfachans);
			free(args);
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
	fits_swap_image_data(out, &gfit); // Efficiently move the merged pixeldata from out to gfit
	clearfits(out);
	free(out);
	siril_log_message(_("Background with %s interpolation computed for CFA image.\n"),
			(args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	/* free memory */
	cfachans_cleanup(cfachans);
	invalidate_stats_from_fit(&gfit);
	if (!args->from_ui) {
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
	}
	siril_add_idle(end_background, args);
	return GINT_TO_POINTER(0);
}

/** Apply for sequence **/

static int background_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct background_data *b_args = (struct background_data*) args->user;

	double *background = (double*)malloc(fit->naxes[0] * fit->naxes[1] * sizeof(double));
	if (!background) {
		PRINT_ALLOC_ERR;
		siril_log_message(_("Out of memory - aborting"));
		return 1;
	}

	const char *err;
	GSList *samples = generate_samples(fit, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE, &err, (threading_type)threads);
	if (!samples) {
		siril_log_color_message(_("Failed to generate background samples for image %d: %s\n"), "red", i, _(err));
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
		siril_log_color_message(_("Warning: RBF background removal is not recommended for CFA images. Only linear background removal is recommended.\n"), "salmon");
	} else if (b_args->degree > 1) {
		siril_log_color_message(_("Warning: polynomial background removal order > 1 is not recommended for CFA images. Only linear background removal is recommended.\n"), "salmon");
	}
	// Obtain CFA pattern as a sensor_pattern
	fits metadata = { 0 };
	if (seq_read_frame_metadata(args->seq, i, &metadata)) {
		siril_log_color_message(_("Error reading metadata.\n"), "red");
		return 1;
	}
	sensor_pattern pattern;
	if (!strncmp(metadata.keywords.bayer_pattern, "RGGB", 4)) {
		pattern = BAYER_FILTER_RGGB;
	} else if (!strncmp(metadata.keywords.bayer_pattern, "BGGR", 4)) {
		pattern = BAYER_FILTER_BGGR;
	} else if (!strncmp(metadata.keywords.bayer_pattern, "GBRG", 4)) {
		pattern = BAYER_FILTER_GBRG;
	} else if (!strncmp(metadata.keywords.bayer_pattern, "GRBG", 4)) {
		pattern = BAYER_FILTER_GBRG;
	} else {
		siril_log_color_message(_("Error: unsupported CFA pattern for this operation.\n"), "red");
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
		siril_log_color_message(_("Error splitting into CFA subcannels, aborting...\n"), "red");
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
			siril_log_message(_("Out of memory - aborting"));
			cfachans_cleanup(cfachans);
			return 1;
		}

		const char *err;
		GSList *samples = generate_samples(subchannel, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE, &err, (threading_type)threads);
		if (!samples) {
			siril_log_color_message(_("Failed to generate background samples for image %d: %s\n"), "red", i, _(err));
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

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
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
		siril_log_color_message(_("Error reading reference metadata.\n"), "red");
		free(background_args->seqEntry);
		free(background_args);
		free_generic_seq_args(args, TRUE);
		return;
	}
	background_args->is_cfa = background_args->seq->nb_layers == 1 && (!strncmp(metadata.keywords.bayer_pattern, "RGGB", 4) ||
							!strncmp(metadata.keywords.bayer_pattern, "BGGR", 4) ||
							!strncmp(metadata.keywords.bayer_pattern, "GBRG", 4) ||
							!strncmp(metadata.keywords.bayer_pattern, "GRBG", 4));
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = background_args->seq->selnum;
	args->compute_mem_limits_hook = background_mem_limits_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = bg_extract_finalize_hook;
	args->image_hook = background_args->is_cfa ? bgcfa_image_hook : background_image_hook;
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
