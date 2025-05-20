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

/* HOW STATISTICS WORK
 * Stats for an image are computed on request and are carried in an imstats
 * structure. Depending on the current mode, this structure is stored for
 * caching in the fits->stats or in the sequence->stats.
 * If it is stored in the fits, when it is disposed, it is copied in the
 * sequence if it belongs to one.
 * If it is stored in the sequence, when the sequence is disposed, it is saved
 * in the seqfile, in `M' fields.
 * When a sequence is loaded, previously computed stats are recovered that way.
 * All operations that need to access stats should do it with the statistics()
 * function at the bottom of this file which provides this abstraction.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "sorting.h"
#include "statistics.h"
#include "demosaicing.h"
#include "gui/progress_and_log.h"

// uncomment to debug statistics
#undef siril_debug_print
#define siril_debug_print(fmt, ...) { }

// copies the area of an image into the memory buffer data
static void select_area_float(fits *fit, float *data, int layer, const rectangle *bounds) {
	int i, j, k = 0;

	const float *from = fit->fpdata[layer] +
		(fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
	int stridefrom = fit->rx - bounds->w;

	for (i = 0; i < bounds->h; ++i) {
		for (j = 0; j < bounds->w; ++j) {
			data[k++] = *from++;
		}
		from += stridefrom;
	}
}

// in this case, N is the number of frames, so int is fine
float siril_stats_float_sd(const float data[], const int N, float *m) {
	double accumulator = 0.0; // accumulating in double precision is important for accuracy
	for (int i = 0; i < N; ++i) {
		accumulator += data[i];
	}
	float mean = (float)(accumulator / N);
	accumulator = 0.0;
	for (int i = 0; i < N; ++i)
		accumulator += (data[i] - mean) * (data[i] - mean);

	if (m) *m = mean;

	return sqrtf((float)(accumulator / (N - 1)));
}

/* For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median
 * of the absolute deviations from the data's median:
 *  MAD = median (| Xi − median(X) |)
 */
double siril_stats_float_mad(const float *data, const size_t n, const double m, threading_type threads, float *buffer) {
	double mad;
	const float median = (float)m;
	float *tmp = buffer ? buffer : malloc(n * sizeof(float));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		return 0.0;
	}

#ifdef _OPENMP
	threads = limit_threading(&threads, 400000, n);
#pragma omp parallel for num_threads(threads) if(threads>1) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		tmp[i] = fabsf(data[i] - median);
	}

	mad = histogram_median_float(tmp, n, threads);
	if (!buffer) {
	    free(tmp);
	}
	return mad;
}

static double siril_stats_float_bwmv(const float* data, const size_t n,
		const float mad, const float median, threading_type threads) {
	double bwmv = 0.0;
	double up = 0.0, down = 0.0;

	if (mad > 0.f) {
        const float factor = 1.f / (9.f * mad);
#ifdef _OPENMP
	threads = limit_threading(&threads, 150000, n);
#pragma omp parallel for num_threads(threads) if(threads>1) schedule(static) reduction(+:up,down)
#endif
		for (size_t i = 0; i < n; i++) {
			const float i_med = data[i] - median;

			const float yi = i_med * factor;
			const float yi2 = fabsf(yi) < 1.f ? yi * yi : 1.f;

			up += SQR(i_med * SQR (1 - yi2));
			down += (1 - yi2) * (1 - 5 * yi2);
		}
		bwmv = down ? n * (up / (down * down)) : 0.0;
	}

	return bwmv;
}

float siril_stats_trmean_from_sorted_data(const float trim,
		const float sorted_data[], const size_t stride, const size_t size) {
	if (trim >= 0.5f) {
		return (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	} else {
		size_t ilow = (size_t) floorf(trim * size);
		size_t ihigh = size - ilow - 1;
		float mean = 0.f;
		float k = 0.f;

		/* compute mean of middle samples in [ilow,ihigh] */
		for (size_t i = ilow; i <= ihigh; ++i) {
			float delta = sorted_data[i * stride] - mean;
			k += 1.f;
			mean += delta / k;
		}
		return mean;
	}
}

#if 0
int IKSS(float *data, size_t n, double *location, double *scale, gboolean multithread) {
	size_t i, j;
	double mad, s, s0, m;
	float xlow, xhigh;

	quicksort_f(data, n);	// this sort is mandatory
	i = 0;
	j = n;
	s0 = 1;
	float *buffer = malloc(n * sizeof(float));
	if (!buffer) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	for (;;) {
		if (j - i < 1) {
			*location = *scale = 0;
			break;
		}
		m = gsl_stats_float_median_from_sorted_data(data + i, 1, j - i);
		mad = siril_stats_float_mad(data + i, j - i, m, multithread, buffer);
		if (mad == 0.0f) {
			free(buffer);
			return 1;
		}
		s = sqrt(siril_stats_float_bwmv(data + i, j - i, mad, m, multithread));
		if (s < 2E-23) {
			*location = m;
			*scale = 0;
			break;
		}
		if (((s0 - s) / s) < 10E-6) {
			*location = m;
			*scale = 0.991 * s;
			break;
		}
		s0 = s;
		xlow = m - 4 * s;
		xhigh = m + 4 * s;
		while (data[i] < xlow)
			i++;
		while (data[j - 1] > xhigh)
			j--;
	}
	free(buffer);
	return 0;
}
#endif

int IKSSlite(float *data, size_t n, const float median, float mad, double *location, double *scale, threading_type threads) {
	size_t i, kept;
	float xlow, xhigh;

	xlow = median - 6.0 * mad;
	xhigh = median + 6.0 * mad;

	kept = 0;
	// removing pixels outside of +/- 6mad from median
	for (i = 0; i < n; i++) {
		if ((data[i] >= xlow) && (data[i] <= xhigh)) {
			if (i != kept) {
				data[kept] = data[i];
			}
			kept++;
		}
	}
	if (kept == 0)
		return 1;

	*location = histogram_median_float(data, kept, threads);
	mad = siril_stats_float_mad(data, kept, *location, threads, NULL);
	if (mad == 0.0f) {
		siril_log_color_message(_("MAD is null. Statistics cannot be computed.\n"), "red");
		return 1;
	}

	*scale = sqrt(siril_stats_float_bwmv(data, kept, mad, *location, threads)) *.991;
	/* 0.991 factor is to keep consistency with IKSS scale */
	return 0;
}

static float* reassign_to_non_null_data_float(float *data, size_t inputlen, size_t outputlen, int free_input) {
	size_t i, j = 0;
	float *ndata = malloc(outputlen * sizeof(float));
	if (!ndata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (i = 0; i < inputlen; i++) {
		if (data[i] != 0.f && !isnan(data[i])) {
			if (j >= outputlen) {
				fprintf(stderr, "\n- stats MISMATCH in sizes (in: %zu, out: %zu), THIS IS A BUG: seqfile is wrong *********\n\n", inputlen, outputlen);
				break;
			}
			ndata[j] = data[i];
			j++;
		}
	}
	if (free_input)
		free(data);
	return ndata;
}

static void siril_stats_float_minmax(float *min_out, float *max_out,
		const float data[], const size_t n, threading_type threads) {
	/* finds the smallest and largest members of a dataset */

	if (n > 0 && data) {
		float min = data[0];
		float max = data[0];
#ifdef _OPENMP
	threads = limit_threading(&threads, 400000, n);
#pragma omp parallel for num_threads(threads) schedule(static) if(threads>1) reduction(max:max) reduction(min:min)
#endif
		for (size_t i = 0; i < n; i++) {
			const float xi = data[i];
			if (!isnan(xi)) {
				min = min(min, xi);
				max = max(max, xi);
			}
		}

		*min_out = min;
		*max_out = max;
	}
}


/* this function tries to get the requested stats from the passed stats,
 * computes them and stores them in it if they have not already been */
imstats* statistics_internal_float(fits *fit, int layer, rectangle *selection, int option, imstats *stats, int bitpix, threading_type threads) {
	int nx = 0, ny = 0;
	float *data = NULL;
	int stat_is_local = 0, free_data = 0;
	imstats* stat = stats;
	// median is included in STATS_BASIC but required to compute other data
	int compute_median = (option & STATS_BASIC) || (option & STATS_AVGDEV) ||
		(option & STATS_MAD) || (option & STATS_BWMV) || (option & STATS_IKSS);

	gboolean valid_selection = selection && selection->h > 0 && selection->w > 0;
	if (!fit && (layer < 0 || valid_selection))
		return NULL;	// not in cache, don't compute

	if (!stat) {
		allocate_stats(&stat);
		if (!stat) return NULL;
		stat_is_local = 1;
	}

	if (fit) {
		if (valid_selection) {
			nx = selection->w;
			ny = selection->h;
			if (layer < 0) {
				size_t newsz;
				data = extract_CFA_buffer_area_float(fit, -layer - 1, selection, &newsz);
				if (!data || newsz == 0) {
					siril_log_message(_("Failed to compute CFA statistics for channel %d\n"), -layer-1);
					return NULL;
				}
				nx = newsz;
				ny = 1;
			} else {
				data = malloc(nx * ny * sizeof(float));
				if (!data) {
					PRINT_ALLOC_ERR;
					if (stat_is_local) free(stat);
					return NULL;
				}
				g_assert(layer < fit->naxes[2]);
				select_area_float(fit, data, layer, selection);
			}
			free_data = 1;
		} else {
			if (layer >= 0) {
				g_assert(layer < fit->naxes[2]);
				nx = fit->rx;
				ny = fit->ry;
				data = fit->fpdata[layer];
			} else {
				/* we just create a buffer containing all pixels that have the
				 * filter number -layer, it's not a real image but ok for stats
				 */
				size_t newsz;
				data = extract_CFA_buffer_float(fit, -layer - 1, &newsz);
				if (!data) {
					siril_log_color_message(_("Failed to compute CFA statistics\n"), "red");
					return NULL;
				}
				nx = newsz;
				ny = 1;
				free_data = 1;
			}
		}
		stat->total = nx * ny;
		if (stat->total == 0L) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			return NULL;
		}
	}

	/* We need to convert stats to the bitdepth of the orginal images
	 * to make sure they are logged in a consistent manner and no mixing 16b/32b can occur in the seq file
	 * - if from fit, the original fit bitpix is passed as argument
	 * - if from cache, so from seq file, we take it from seq->bitpix and pass it as argument
	 * This is really important because we can compute float statistics from ushort images.
	 * Fixes https://gitlab.com/free-astro/siril/-/issues/700
	 */

	if (stat->normValue == NULL_STATS) {
		if (bitpix == FLOAT_IMG) stat->normValue = 1.f;
		else stat->normValue = (bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	}

	/* Calculation of min and max */
	if ((option & (STATS_MINMAX | STATS_BASIC)) && (stat->min == NULL_STATS || stat->max == NULL_STATS)) {
		float min = 0, max = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing minmax\n", stat, fit, layer);
		siril_stats_float_minmax(&min, &max, data, stat->total, threads);
		stat->min = (double)min * stat->normValue;
		stat->max = (double)max * stat->normValue;
	}

	/* Calculation of ngoodpix, mean, sigma and background noise */
	if ((option & (STATS_SIGMEAN | STATS_BASIC)) && (stat->ngoodpix <= 0L || stat->mean == NULL_STATS ||
			stat->sigma == NULL_STATS || stat->bgnoise == NULL_STATS)) {
		int status = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing basic\n", stat, fit, layer);
		siril_fits_img_stats_float(data, nx, ny, &stat->ngoodpix,
				NULL, NULL, &stat->mean, &stat->sigma, &stat->bgnoise,
				NULL, NULL, NULL, threads, &status);
		if (status) {
			if (free_data) free(data);
			if (stat_is_local) free(stat);
			siril_log_message("fits_img_stats_float failed\n");
			return NULL;
		}
		stat->mean *= stat->normValue;
		stat->sigma *= stat->normValue;
		stat->bgnoise *= stat->normValue;
	}
	if (stat->ngoodpix == 0L) {
		if (free_data) free(data);
		if (stat_is_local) free(stat);
		return NULL;
	}


	/* we exclude 0 and nans if some computations remain to be done or copy data if
	 * median has to be computed */
	if (fit && compute_median && stat->total != stat->ngoodpix) {
		data = reassign_to_non_null_data_float(data, stat->total, stat->ngoodpix, free_data);
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;
		}
		free_data = 1;
	}


	/* Calculation of median */
	if (compute_median && stat->median == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing median\n", stat, fit, layer);
		stat->median = histogram_median_float(data, stat->ngoodpix, threads) * stat->normValue;
	}

	/* Calculation of average absolute deviation from the median */
	if ((option & STATS_AVGDEV) && stat->avgDev == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing absdev\n", stat, fit, layer);
		stat->avgDev = gsl_stats_float_absdev_m(data, 1, stat->ngoodpix, stat->median / stat->normValue) * stat->normValue;
	}

	/* Calculation of median absolute deviation */
	if (((option & STATS_MAD) || (option & STATS_BWMV) || (option & STATS_IKSS)) && stat->mad == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing mad\n", stat, fit, layer);
		stat->mad = siril_stats_float_mad(data, stat->ngoodpix, stat->median / stat->normValue, threads, NULL) * stat->normValue;
	}

	/* Calculation of Bidweight Midvariance */
	if ((option & STATS_BWMV) && stat->sqrtbwmv == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing bimid\n", stat, fit, layer);
		double bwmv = siril_stats_float_bwmv(data, stat->ngoodpix, stat->mad / stat->normValue, stat->median / stat->normValue, threads);
		stat->sqrtbwmv = sqrt(bwmv) * stat->normValue;
	}

	/* Calculation of IKSS. Only used for stacking normalization */
	if ((option & STATS_IKSS) && (stat->location == NULL_STATS || stat->scale == NULL_STATS)) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing ikss\n", stat, fit, layer);
		if (IKSSlite(data, stat->ngoodpix, stat->median / stat->normValue, stat->mad / stat->normValue, &stat->location, &stat->scale, threads)) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			return
			 NULL;
		}
		stat->location *= stat->normValue;
		stat->scale *= stat->normValue;
	}

	if (free_data) free(data);
	return stat;
}

int compute_means_from_flat_cfa_float(const fits *fit, double mean[36]) {
	int row, col, c, i[36] = {0};
	float *data;
	unsigned int width, height;
	unsigned int startx, starty;

	data = fit->fdata;
	width = fit->rx;
	height = fit->ry;

	/* due to vignetting it is better to take an area in the
	 * center of the flat image
	 */
	startx = width / 3;
	starty = height / 3;

	siril_debug_print("Computing stat in (%d, %d, %d, %d)\n", startx, starty,
			width - 1 - startx, height - 1 - starty);

	/* compute mean of each element in 6x6 blocks */
	for (row = starty; row < height - starty; row++) {
		for (col = startx; col < width - startx; col++) {
			mean[(col % 6) + (row % 6) * 6] += data[col + row * width];
			i[(col % 6) + (row % 6) * 6]++;
		}
	}

	for (c = 0; c < 36; c++) {
		mean[c] /= (double) i[c];
	}
	return 0;
}

/****************** ROBUST MEAN FOR PHOTOMETRY *******************/

#define hampel_a   1.7
#define hampel_b   3.4
#define hampel_c   8.5
#define sign(x,y)  ((y)>=0?fabs(x):-fabs(x))
#define epsilon(x) 0.00000001
#define maxit      50

static double hampel(double x) {
	if (x >= 0) {
		if (x < hampel_a)
			return x;
		if (x < hampel_b)
			return hampel_a;
		if (x < hampel_c)
			return hampel_a * (x - hampel_c) / (hampel_b - hampel_c);
	} else {
		if (x > -hampel_a)
			return x;
		if (x > -hampel_b)
			return -hampel_a;
		if (x > -hampel_c)
			return hampel_a * (x + hampel_c) / (hampel_b - hampel_c);
	}
	return 0.0;
}

static double dhampel(double x) {
	if (x >= 0) {
		if (x < hampel_a)
			return 1;
		if (x < hampel_b)
			return 0;
		if (x < hampel_c)
			return hampel_a / (hampel_b - hampel_c);
	} else {
		if (x > -hampel_a)
			return 1;
		if (x > -hampel_b)
			return 0;
		if (x > -hampel_c)
			return -hampel_a / (hampel_b - hampel_c);
	}
	return 0.0;
}

static double qmedD(int n, double *a)
	/* Vypocet medianu algoritmem Quick Median (Wirth) */
{
	double w;
	int k = ((n & 1) ? (n / 2) : ((n / 2) - 1));
	int l = 0;
	int r = n - 1;

	while (l < r) {
		double x = a[k];
		int i = l;
		int j = r;
		do {
			while (a[i] < x)
				i++;
			while (x < a[j])
				j--;
			if (i <= j) {
				w = a[i];
				a[i] = a[j];
				a[j] = w;
				i++;
				j--;
			}
		} while (i <= j);
		if (j < k)
			l = i;
		if (k < i)
			r = j;
	}
	return a[k];
}

static float Qn0(const float sorted_data[], const size_t stride, const size_t n) {
	const size_t wsize = n * (n - 1) / 2;
	const size_t n_2 = n / 2;
	const size_t k = ((n_2 + 1) * n_2) / 2;
	size_t idx = 0;

	if (n < 2)
		return (0.0);

	float *work = malloc(wsize * sizeof(float));
	if (!work) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j)
			work[idx++] = fabsf(sorted_data[i] - sorted_data[j]);
	}

	quicksort_f(work, idx);
	float Qn = work[k - 1];

	free(work);
	return Qn;
}

/* a robust mean from sorted data, see also robust_mean in sorting.c */
float siril_stats_robust_mean(const float sorted_data[],
		const size_t stride, const size_t size, double *deviation) {
	float mx = (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	float qn0 = Qn0(sorted_data, 1, size);
	if (qn0 < 0)
		return -1.0f;
	float sx = 2.2219f * qn0;
	float *x, mean;
	int i, j;

	x = malloc(size * sizeof(float));
	if (!x) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (i = 0, j = 0; i < size; ++i) {
		if (fabsf(sorted_data[i] - mx) <= 3.f * sx) {
			x[j++] = sorted_data[i];
		}
	}
	siril_debug_print("keeping %d samples on %zu for the robust mean (mx: %f, sx: %f)\n", j, size, mx, sx);
	/* not enough stars, try something anyway */
	if (j < 5) {
		mean = siril_stats_trmean_from_sorted_data(0.3f, sorted_data, stride, size);
	} else {
		mean = (float) gsl_stats_float_mean(x, stride, j);
	}
	free(x);

	/* compute the deviation of the mean against the values */
	if (deviation) {
		double dev = 0.0;
		int inliers = 0;
		for (i = 0, j = 0; i < size; ++i) {
			if (fabsf(sorted_data[i] - mx) <= 3.f * sx) {
				dev += fabsf(sorted_data[i] - mean);
				inliers++;
			}
		}
		if (inliers)
			*deviation = dev / inliers;
		else *deviation = 1.0;
	}

	return mean;
}

/************* another robust mean ***************/

int robustmean(int n, const double *x, double *mean, double *stdev)
	/* Newton's iterations */
{
	int i, it;
	double a, c, dt, r, s, psir;
	double *xx;

	if (n < 1) {
		if (mean)
			*mean = 0.0; /* a few data */
		if (stdev)
			*stdev = -1.0;
		return 1;
	}
	if (n == 1) { /* only one point, but correct case */
		if (mean)
			*mean = x[0];
		if (stdev)
			*stdev = 0.0;
		return 0;
	}

	/* initial values:
	   - median is the first approximation of location
	   - MAD/0.6745 is the first approximation of scale */
	xx = malloc(n * sizeof(double));
	if (!xx) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	memcpy(xx, x, n * sizeof(double));
	a = qmedD(n, xx);
	for (i = 0; i < n; i++)
		xx[i] = fabs(x[i] - a);
	s = qmedD(n, xx) / 0.6745;
	free(xx);

	/* almost identical points on input */
	if (fabs(s) < epsilon(s)) {
		if (mean)
			*mean = a;
		if (stdev) {
			double sum = 0.0;
			for (i = 0; i < n; i++)
				sum += (x[i] - a) * (x[i] - a);
			*stdev = sqrt(sum / n);
		}
		return 0;
	}

	/* corrector's estimation */
	dt = 0;
	c = s * s * n * n / (n - 1);
	for (it = 1; it <= maxit; it++) {
		double sum1, sum2, sum3;
		sum1 = sum2 = sum3 = 0.0;
		for (i = 0; i < n; i++) {
			r = (x[i] - a) / s;
			psir = hampel(r);
			sum1 += psir;
			sum2 += dhampel(r);
			sum3 += psir * psir;
		}
		if (fabs(sum2) < epsilon(sum2))
			break;
		double d = s * sum1 / sum2;
		a = a + d;
		dt = c * sum3 / (sum2 * sum2);
		if ((it > 2) && ((d * d < 1e-4 * dt) || (fabs(d) < 10.0 * epsilon(d))))
			break;
	}
	if (mean)
		*mean = a;
	if (stdev)
		*stdev = (dt > 0 ? sqrt(dt) : 0);
	return 0;
}

double robust_median_f(fits *fit, rectangle *area, int chan, float lower, float upper) {
	uint32_t x0, y0, x1,y1;
	if (area) {
		x0 = area->x;
		y0 = area->y;
		x1 = area->x + area->w;
		y1 = area->y + area->h;
	} else {
		x0 = y0 = 0;
		x1 = gfit.rx;
		y1 = gfit.ry;
	}
	size_t npixels = (x1 - x0) * (y1 - y0);
	float *data = fit->fpdata[chan];
	float *filtered_data = malloc(npixels * sizeof(float));
	size_t count = 0;
	for (uint32_t y = y0 ; y < y1 ; y++) {
		size_t j = y * fit->rx;
		for (uint32_t x = x0 ; x < x1 ; x++) {
			size_t i = x + j;
			if (data[i] >= lower && data[i] <= upper) {
				filtered_data[count++] = data[i];
			}
		}
	}
	// Check if there are any elements in the specified range
	if (count == 0) {
		free(filtered_data);
		return 0.0; // No elements in the range, return 0 as median
	}
	// Sort the filtered data
	// use histogram_median_float here instead of quickmedian for speed (see #1458)
	double retval = histogram_median_float(filtered_data, count, MULTI_THREADED);

	// Free the allocated memory for filtered_data
	free(filtered_data);

	return retval;
}
