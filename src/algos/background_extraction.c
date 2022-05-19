/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "algos/statistics.h"
#include "algos/geometry.h"
#include "algos/sorting.h"
#include "opencv/opencv.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "registration/registration.h"	// for mouse_status
#include "background_extraction.h"

static gboolean background_computed = FALSE;

static void background_startup() {
	copy_gfit_to_backup();
}

#define NPARAM_POLY4 15		// Number of parameters used with 4rd order
#define NPARAM_POLY3 10		// Number of parameters used with 3rd order
#define NPARAM_POLY2 6		// Number of parameters used with 2nd order
#define NPARAM_POLY1 3		// Number of parameters used with 1nd order

#define SAMPLE_SIZE 25		// must be odd to compute a radius

struct sample {
	double median[3]; // median of each channel of the sample (if color)
	double mean; // mean of the 3 channel of the sample (if color)
	double min, max;
	size_t size;
	point position;
	gboolean valid;
};

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

	char *msg = siril_log_color_message(_("RBF Extraction: processing channel %d...\n"), "green", channel);
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);

	double pixel;
	gsl_matrix *K;
	gsl_vector *f, *coef, *A;
	GSList *l_i;
	guint n = g_slist_length(list);
	double mean = 0.0;
	double (*list_array)[3];

	K = gsl_matrix_calloc(n + 1, n + 1);
	f = gsl_vector_calloc(n + 1);
	coef = gsl_vector_calloc(n + 1);

	/* Scaling */
	int scaling_factor = 4;
	int width_scaled = round_to_int(width / scaling_factor);
	int height_scaled = round_to_int(height / scaling_factor);
	
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

	int s;
	gsl_permutation *p = gsl_permutation_alloc(n + 1);
	gsl_linalg_LU_decomp(K, p, &s);
	gsl_linalg_LU_solve(K, p, f, coef);

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
	double total = height_scaled * width_scaled;
	int progress = 0;

#pragma omp parallel shared(background_scaled) private(A) num_threads(threads)
	{
		A = gsl_vector_calloc(n + 1);
#pragma omp for
		for (int i = 0; i < height_scaled; i++) {
			for (int j = 0; j < width_scaled; j++) {
				for (int k = 0; k < n; k++) {
					int deltax = abs(j - (int)list_array[k][0]);
					int deltay = abs(i - (int)list_array[k][1]);
					gsl_vector_set(A, k, kernel_scaled[deltax + deltay * width_scaled]);
				}

				gsl_vector_set(A, n, 1.0);
				gsl_vector_mul(A, coef);

				pixel = siril_gsl_vector_sum(A);
				background_scaled[j + i * width_scaled] = pixel;

#ifdef _OPENMP
#pragma omp critical
#endif
				{
					++progress;
					if (!(progress % 32)) {
						set_progress_bar_data(NULL, (double) progress / total);
					}
				}

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

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

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

	// Must turn off error handler or it aborts on error
	gsl_set_error_handler_off();

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
	background_sample *sample = (background_sample *) g_malloc(sizeof(background_sample));
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

static unsigned int _rand(guint64 *const p_rng) {
	*p_rng = *p_rng * 1103515245 + 12345U;
	return (unsigned int) *p_rng;
}

static gboolean convert_fits_to_img(fits *fit, double *image, int channel, gboolean add_dither) {
	guint64 seed = time(NULL);

	double invnorm = 1.0 / USHRT_MAX;
	const int height = fit->ry;
	const int width = fit->rx;
	if (fit->type == DATA_USHORT) {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				image[y * width + x] = fit->pdata[channel][(height - y - 1) * width + x] * invnorm;
				if (add_dither) {
					/* add dithering in order to avoid colour banding */
					image[y * width + x] += (_rand(&seed) % 1048576) * 0.000000000095367431640625f;
				}
			}
		}
	} else {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				image[y * width + x] = fit->fpdata[channel][(height - y - 1) * width + x];
				if (add_dither) {
					/* add dithering in order to avoid colour banding */
					image[y * width + x] += (_rand(&seed) % 1048576) * 0.000000000095367431640625f;
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

static GSList *generate_samples(fits *fit, int nb_per_line, double tolerance, int size, threading_type threads) {
	int nx = fit->rx;
	int ny = fit->ry;
	size_t n = fit->naxes[0] * fit->naxes[1];
	GSList *list = NULL;

	float *image = convert_fits_to_luminance(fit, threads);	// upside down
	if (!image) return NULL;
	float median = (float)histogram_median_float(image, n, threads);

	int boxes_width = nb_per_line * size + 2;	// leave a margin of 1 px on the sides
	float spacing = (nx - boxes_width) / (float)(nb_per_line-1);
	int radius = size / 2;
	// solving iteratively n for: ny - 2 = n * size + (n-1) * spacing;
	int nb_per_column = 1;
	while (nb_per_column * size + round_to_int((nb_per_column - 1) * spacing) <= (ny - 2))
		nb_per_column++;
	nb_per_column--;
	if (nb_per_column == 0) {
		siril_debug_print("image is smaller than the sample size");
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
			g_free(sample);
			list = g_slist_delete_link(list, l);
		}
		l = next;
	}

	list = g_slist_reverse(list);
	free(image);
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

static poly_order get_poly_order() {
	GtkComboBox *combo_box_poly_order = GTK_COMBO_BOX(lookup_widget("box_background_order"));
	return gtk_combo_box_get_active(combo_box_poly_order);
}

static background_correction get_correction_type() {
	GtkComboBox *combo_box_correction = GTK_COMBO_BOX(lookup_widget("box_background_correction"));
	return gtk_combo_box_get_active(combo_box_correction);
}

static int get_nb_samples_per_line() {
	GtkSpinButton *nb_samples = GTK_SPIN_BUTTON(lookup_widget("spin_background_nb_samples"));
	return gtk_spin_button_get_value_as_int(nb_samples);
}

static double get_tolerance_value() {
	GtkRange *tol = GTK_RANGE(lookup_widget("scale_background_tolerance"));
	return gtk_range_get_value(tol);
}

static background_interpolation get_interpolation_method() {
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("background_extraction_combo"));
	return gtk_combo_box_get_active(combo);
}

static double get_smoothing_parameter() {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_background_smoothing"));
	return gtk_spin_button_get_value(spin);
}

static void remove_gradient(double *img, const double *background, double background_mean, size_t ndata, background_correction type, threading_type threads) {
	size_t i;
	double mean;

	switch (type) {
		default:
		case BACKGROUND_CORRECTION_SUBTRACT:
#ifdef _OPENMP
			limit_threading(&threads, 300000, ndata);
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				img[i] -= background[i];
				img[i] += background_mean;
			}
			break;
		case BACKGROUND_CORRECTION_DIVIDE:
			mean = gsl_stats_mean(img, 1, ndata);
#ifdef _OPENMP
			limit_threading(&threads, 300000, ndata);
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				img[i] /= background[i];
				img[i] *= mean;
			}
	}
}

static gboolean is_dither_checked() {
	return (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bkg_dither_button"))));
}

/************* PUBLIC FUNCTIONS *************/

int get_background_sample_radius() {
       return SAMPLE_SIZE / 2;
}

void free_background_sample_list(GSList *list) {
	if (list == NULL) return;
	g_slist_free_full(list, g_free);
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
		}
	}
	free(image);

	return orig;
}

/* generates samples and stores them in com.grad_samples */
void generate_background_samples(int nb_of_samples, double tolerance) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = generate_samples(&gfit, nb_of_samples, tolerance, SAMPLE_SIZE, MULTI_THREADED);

	if (com.grad_samples && gfit.naxes[2] > 1) {
		/* If RGB we need to update all local median, not only the first one */
		com.grad_samples = update_median_samples(com.grad_samples, &gfit);
	}
	redraw(REDRAW_OVERLAY);
}

static gboolean end_background(gpointer p) {
	stop_processing_thread();
	struct background_data *args = (struct background_data *)p;
	invalidate_stats_from_fit(args->fit);
	background_computed = TRUE;
	if (!args->from_ui) {
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
	}

	notify_gfit_modified();
	free(args);
	return FALSE;
}


/* uses samples from com.grad_samples */
gpointer remove_gradient_from_image(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	gchar *error;
	double *background = malloc(gfit.ry * gfit.rx * sizeof(double));
	
	if (!background && !com.script) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return NULL;
	}

	const size_t n = gfit.naxes[0] * gfit.naxes[1];
	double *image = malloc(n * sizeof(double));
	if (!image) {
		free(background);
		PRINT_ALLOC_ERR;
		return NULL;
	}
	
	/* Make sure to update local median. Useful if undo is pressed */
	update_median_samples(com.grad_samples, &gfit);

	double background_mean = get_background_mean(com.grad_samples, gfit.naxes[2]);
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	for (int channel = 0; channel < gfit.naxes[2]; channel++) {
		/* compute background */
		gboolean interpolation_worked = TRUE;
		if (args->interpolation_method == BACKGROUND_INTER_POLY){
			interpolation_worked = computeBackground_Polynom(com.grad_samples, background, channel, 
									gfit.rx, gfit.ry, args->degree, &error);
		} else {
			interpolation_worked = computeBackground_RBF(com.grad_samples, background, channel, 
									gfit.rx, gfit.ry, args->smoothing, &error, args->threads);
		}
		
		if (!interpolation_worked) {
			free(image);
			free(background);
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Not enough samples."),
					error);
			set_cursor_waiting(FALSE);
			return NULL;
		}
		/* remove background */
		const char *c_name = vport_number_to_name(channel);
		siril_log_message(_("Background extraction from %s channel.\n"), c_name);
		convert_fits_to_img(&gfit, image, channel, args->dither);
		remove_gradient(image, background, background_mean, n, args->correction, MULTI_THREADED);
		convert_img_to_fits(image, &gfit, channel);
	}
	siril_log_message(_("Background with %s interpolation computed.\n"), (args->interpolation_method == BACKGROUND_INTER_POLY) ? "polynomial" : "RBF");
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	/* free memory */
	free(image);
	free(background);
	siril_add_idle(end_background, args);
	return args;
}

/** Apply for sequence **/

static int background_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct background_data *b_args = (struct background_data*) args->user;

	gchar *error;
	double *background = (double*)malloc(fit->naxes[0] * fit->naxes[1] * sizeof(double));
	if (!background) {
		PRINT_ALLOC_ERR;
		error = _("Out of memory - aborting");
		siril_log_message(error);
		return 1;
	}

	GSList *samples = generate_samples(fit, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE, (threading_type)threads);
	if (!samples) {
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
		 */
		uint64_t double_channel_size = args->seq->rx * args->seq->ry * sizeof(double);
		unsigned int double_channel_size_MB = double_channel_size / BYTES_IN_A_MB;
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

void apply_background_extraction_to_sequence(struct background_data *background_args) {
	struct generic_seq_args *args = create_default_seqargs(background_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = background_args->seq->selnum;
	args->compute_mem_limits_hook = background_mem_limits_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = background_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Background Extraction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = background_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->user = background_args;

	background_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

/**** getter and setter ***/

gboolean background_sample_is_valid(background_sample *sample) {
	return sample->valid;
}

gdouble background_sample_get_size(background_sample *sample) {
	return sample->size;
}

point background_sample_get_position(background_sample *sample) {
	return sample->position;
}

/************* CALLBACKS *************/

void on_background_generate_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	int nb_of_samples;
	double tolerance;

	nb_of_samples = get_nb_samples_per_line();
	tolerance = get_tolerance_value();

	generate_background_samples(nb_of_samples, tolerance);
	set_cursor_waiting(FALSE);
}

void on_background_clear_all_clicked(GtkButton *button, gpointer user_data) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;

	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
}

void on_bkg_compute_bkg_clicked(GtkButton *button, gpointer user_data) {
	if (com.grad_samples == NULL) {
		return;
	}
	set_cursor_waiting(TRUE);
	copy_backup_to_gfit();

	background_correction correction = get_correction_type();
	poly_order degree = get_poly_order();
	gboolean use_dither = is_dither_checked();
	double smoothing = get_smoothing_parameter();
	background_interpolation interpolation_method = get_interpolation_method();

	struct background_data *args = malloc(sizeof(struct background_data));
	args->threads = com.max_thread;
	args->from_ui = TRUE;
	args->correction = correction;
	args->interpolation_method = interpolation_method;
	args->degree = (poly_order) (degree - 1);
	args->smoothing = smoothing;
	args->dither = use_dither;
	args->fit = &gfit;

	start_in_new_thread(remove_gradient_from_image, args);
}

void on_background_ok_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkBkgSeq"));
	if (gtk_toggle_button_get_active(seq_button) && sequence_is_loaded()) {
		struct background_data *args = malloc(sizeof(struct background_data));
		args->nb_of_samples = get_nb_samples_per_line();
		args->tolerance = get_tolerance_value();
		args->correction = get_correction_type();
		args->degree = get_poly_order();
		args->smoothing = get_smoothing_parameter();
		args->dither = is_dither_checked();
		args->interpolation_method = get_interpolation_method();
		
		if (args->interpolation_method == BACKGROUND_INTER_POLY && args->degree > BACKGROUND_POLY_1) {
			int confirm = siril_confirm_dialog(_("Polynomial order seems too high."),
					_("You are about to process a sequence of preprocessed files with "
							"a polynomial degree greater than 1. This is unlikely because such "
							"gradients are often linear and a correction with a polynomial "
							"function of degree 1 is probably enough."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		} else if (args->interpolation_method != BACKGROUND_INTER_POLY) {
			int confirm = siril_confirm_dialog(_("Using wrong interpolation method"),
					_("You are about to process a sequence of preprocessed files with an RBF algorithm. "
							"This algorithm may not be very well suited for automated processing "
							"and we advise you to use the polynomial algorithm with a "
							"degree order of 1."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		}

		set_cursor_waiting(TRUE);

		args->seqEntry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryBkgSeq")));
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "bkg_";
		args->seq = &com.seq;
		/* now we uncheck the button */
		gtk_toggle_button_set_active(seq_button, FALSE);
		apply_background_extraction_to_sequence(args);
	} else {
		if (background_computed) {
			background_correction correction = get_correction_type();
			undo_save_state(get_preview_gfit_backup(), _("Background extraction (Correction: %s)"),
					correction == BACKGROUND_CORRECTION_DIVIDE ? "Division" : "Subtraction");
			background_computed = FALSE;
			clear_backup();
			siril_close_dialog("background_extraction_dialog");
		} else {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("No Background model computed"),
					_("You must first compute the background model."));
		}
	}
}

void apply_background_cancel() {
	siril_close_dialog("background_extraction_dialog");
}

void on_background_close_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("background_extraction_dialog");
}

void on_background_extraction_dialog_hide(GtkWidget *widget, gpointer user_data) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(REDRAW_OVERLAY);

	if (background_computed) {
		siril_preview_hide();
	} else {
		clear_backup();
	}
}

void on_background_extraction_dialog_show(GtkWidget *widget, gpointer user_data) {
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	background_startup();
}

void on_background_extraction_combo_changed(GtkComboBox *combo, gpointer user_data) {
	GtkNotebook *notebook = GTK_NOTEBOOK(lookup_widget("bkg_notebook_inter"));
	gtk_notebook_set_current_page(notebook, gtk_combo_box_get_active(combo));
}

