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
#include <gsl/gsl_statistics.h>
#include "algos/sorting.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "statistics.h"
#include "statistics_float.h"
#include "demosaicing.h"
#include "core/OS_utils.h"

// comment to debug statistics
#undef siril_debug_print
#define siril_debug_print(fmt, ...) { }

static void stats_set_default_values(imstats *stat);

// copies the area of an image into the memory buffer data
static void select_area_ushort(fits *fit, WORD *data, int layer, rectangle *bounds) {
	int i, j, k = 0;

	WORD *from = fit->pdata[layer] +
		(fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
	int stridefrom = fit->rx - bounds->w;

	for (i = 0; i < bounds->h; ++i) {
		for (j = 0; j < bounds->w; ++j) {
			data[k++] = *from++;
		}
		from += stridefrom;
	}
}

// sd computation for stacking.
// In this case, N is the number of frames, so int is fine
float siril_stats_ushort_sd_64(const WORD data[], const int N) {
	guint64 intsum = 0;
	for (int i = 0; i < N; ++i) {
		intsum += data[i];
	}
	float mean = (float)(((double)intsum) / ((double)N));
	double accumulator = 0.0;
	for (int i = 0; i < N; ++i) {
		float pixel = (float)data[i];
		accumulator += (pixel - mean) * (pixel - mean);
	}
	return sqrtf((float) (accumulator / (N - 1)));
}

// 32 bits sum version, enough if there are less than 64k images
float siril_stats_ushort_sd_32(const WORD data[], const int N) {
	guint32 intsum = 0;
	for (int i = 0; i < N; ++i) {
		intsum += data[i];
	}
	float mean = (float)(((double)intsum) / ((double)N));
	double accumulator = 0.0;
	for (int i = 0; i < N; ++i) {
		float pixel = (float)data[i];
		accumulator += (pixel - mean) * (pixel - mean);
	}
	return sqrtf((float) (accumulator / (N - 1)));
}

/* For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median
 * of the absolute deviations from the data's median:
 *  MAD = median (| Xi − median(X) |)
 */
float siril_stats_ushort_mad(const WORD* data, const size_t n, const double m,
		threading_type threads) {
	float mad;
	int median = round_to_int(m);	// we use it on integer data anyway
	WORD *tmp = malloc(n * sizeof(WORD));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		return 0.0f;
	}

#ifdef _OPENMP
	threads = limit_threading(&threads, 400000, n);
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		tmp[i] = (WORD)abs(data[i] - median);
	}

	mad = (float) histogram_median(tmp, n, threads);
	free(tmp);
	return mad;
}

static double siril_stats_ushort_bwmv(const WORD* data, const size_t n,
		const double mad, const double median, threading_type threads) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;

	if (mad > 0.0) {
#ifdef _OPENMP
	threads = limit_threading(&threads, 150000, n);
#pragma omp parallel for num_threads(threads) schedule(static) reduction(+:up,down) if(threads>1)
#endif
		for (size_t i = 0; i < n; i++) {
			double yi, ai, yi2;

			yi = ((double) data[i] - median) / (9 * mad);
			yi2 = yi * yi;
			ai = (fabs(yi) < 1.0) ? 1.0 : 0.0;

			up += ai * SQR((double ) data[i] - median) * SQR(SQR (1 - yi2));
			down += (ai * (1 - yi2) * (1 - 5 * yi2));
		}
		bwmv = down ? n * (up / (down * down)) : 0.0;
	}

	return bwmv;
}

static WORD* reassign_to_non_null_data_ushort(WORD *data, size_t inputlen, size_t outputlen, int free_input) {
	size_t i, j = 0;
	WORD *ndata = malloc(outputlen * sizeof(WORD));
	if (!ndata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (i = 0; i < inputlen; i++) {
		if (data[i] > 0) {
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

static void siril_stats_ushort_minmax(WORD *min_out, WORD *max_out,
		const WORD data[], const size_t n, threading_type threads) {
	/* finds the smallest and largest members of a dataset */

	if (n > 0 && data) {
		WORD min = data[0];
		WORD max = data[0];

#ifdef _OPENMP
	threads = limit_threading(&threads, 400000, n);
#pragma omp parallel for num_threads(threads) schedule(static) if(threads>1) reduction(max:max) reduction(min:min)
#endif
		for (size_t i = 0; i < n; i++) {
			WORD xi = data[i];
			if (xi < min)
				min = xi;
			if (xi > max)
				max = xi;
		}

		*min_out = min;
		*max_out = max;
	}
}

/* this function tries to get the requested stats from the passed stats,
 * computes them and stores them in it if they have not already been */
static imstats* statistics_internal_ushort(fits *fit, int layer, rectangle *selection,
		int option, imstats *stats, int bitpix, threading_type threads) {
	int nx = 0, ny = 0;
	WORD *data = NULL;
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
				data = extract_CFA_buffer_area_ushort(fit, -layer - 1, selection, &newsz);
				if (!data || newsz == 0) {
					siril_log_message(_("Failed to compute CFA statistics for channel %d\n"), -layer-1);
					return NULL;
				}
				nx = newsz;
				ny = 1;
			} else {
				data = malloc(nx * ny * sizeof(WORD));
				if (!data) {
					PRINT_ALLOC_ERR;
					if (stat_is_local) free(stat);
					return NULL;
				}
				g_assert(layer < fit->naxes[2]);
				select_area_ushort(fit, data, layer, selection);
			}
			free_data = 1;
		} else {
			if (layer >= 0) {
				g_assert(layer < fit->naxes[2]);
				nx = fit->rx;
				ny = fit->ry;
				data = fit->pdata[layer];
			} else {
				/* we just create a buffer containing all pixels that have the
				 * filter number -layer, it's not a real image but ok for stats
				 */
				size_t newsz;
				data = extract_CFA_buffer_ushort(fit, -layer - 1, &newsz);
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

	if (stat->normValue == NULL_STATS) {
		stat->normValue = (bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	}

	/* Calculation of min and max */
	if ((option & (STATS_MINMAX | STATS_BASIC)) && (stat->min == NULL_STATS || stat->max == NULL_STATS)) {
		WORD min = 0, max = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing minmax\n", stat, fit, layer);
		siril_stats_ushort_minmax(&min, &max, data, stat->total, threads);
		stat->min = (double) min;
		stat->max = (double) max;
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
		siril_fits_img_stats_ushort(data, nx, ny, &stat->ngoodpix,
				NULL, NULL, &stat->mean, &stat->sigma, &stat->bgnoise,
				NULL, NULL, NULL, threads, &status);
		if (status) {
			if (free_data) free(data);
			if (stat_is_local) free(stat);
			return NULL;
		}
	}

	if (stat->ngoodpix == 0L) {
		if (free_data) free(data);
		if (stat_is_local) free(stat);
		return NULL;
	}

	/* we exclude 0 if some computations remain to be done or copy data if
	 * median has to be computed (this is deactivated in the ngoodpix computation) */
	if (fit && compute_median && stat->total != stat->ngoodpix) {
		data = reassign_to_non_null_data_ushort(data, stat->total, stat->ngoodpix, free_data);
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
		stat->median = histogram_median(data, stat->ngoodpix, threads);
	}

	/* Calculation of average absolute deviation from the median */
	if ((option & STATS_AVGDEV) && stat->avgDev == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing absdev\n", stat, fit, layer);
		stat->avgDev = gsl_stats_ushort_absdev_m(data, 1, stat->ngoodpix, stat->median);
	}

	/* Calculation of median absolute deviation */
	if (((option & STATS_MAD) || (option & STATS_BWMV) || (option & STATS_IKSS)) && stat->mad == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing mad\n", stat, fit, layer);
		stat->mad = siril_stats_ushort_mad(data, stat->ngoodpix, stat->median, threads);
	}

	/* Calculation of Bidweight Midvariance */
	if ((option & STATS_BWMV) && stat->sqrtbwmv == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing bimid\n", stat, fit, layer);
		double bwmv = siril_stats_ushort_bwmv(data, stat->ngoodpix, stat->mad, stat->median, threads);
		stat->sqrtbwmv = sqrt(bwmv);
	}


	/* Calculation of IKSS. Only used for stacking normalization */
	if ((option & STATS_IKSS) && (stat->location == NULL_STATS || stat->scale == NULL_STATS)) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing ikss\n", stat, fit, layer);
		float *newdata = malloc(stat->ngoodpix * sizeof(float));
		if (!newdata) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			PRINT_ALLOC_ERR;
			return NULL;
		}
		double normValue = (fit->bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
		/* we convert in the [0, 1] range */
		float invertNormValue = (float)(1.0 / normValue);
#ifdef _OPENMP
		int loopthreads = limit_threading(&threads, 400000, stat->ngoodpix);
#pragma omp parallel for num_threads(loopthreads) if (loopthreads>1) schedule(static)
#endif
		for (size_t i = 0; i < stat->ngoodpix; i++) {
			newdata[i] = (float) data[i] * invertNormValue;
		}

		float med = (float)(stat->median) * invertNormValue;
		float mad = (float)(stat->mad) * invertNormValue;
		if (IKSSlite(newdata, stat->ngoodpix, med, mad, &stat->location, &stat->scale, threads)) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			free(newdata);
			return NULL;
		}
		/* go back to the original range */
		stat->location *= normValue;
		stat->scale *= normValue;
		free(newdata);
	}

	if (free_data) free(data);
	return stat;
}

static imstats* statistics_internal(fits *fit, int layer, rectangle *selection, int option, imstats *stats, int bitpix, threading_type threads) {
	if (fit) {
		if (fit->type == DATA_USHORT)
			return statistics_internal_ushort(fit, layer, selection, option, stats, bitpix, threads);
		if (fit->type == DATA_FLOAT)
			return statistics_internal_float(fit, layer, selection, option, stats, bitpix, threads);
		return NULL;
	}
	if (bitpix == FLOAT_IMG)
		return statistics_internal_float(fit, layer, selection, option, stats, bitpix, threads);
	return statistics_internal_ushort(fit, layer, selection, option, stats, bitpix, threads);
}

/* Computes statistics on the given layer of the given opened image.
 * Mean, sigma and noise are computed with a cfitsio function rewritten here.
 * Min and max value, average deviation, MAD, Bidweight Midvariance and IKSS
 * are computed with gsl stats.
 *
 * layer can be negative, this is a special trick to get statistics on a CFA
 * channel, -1 being red, -2 green and -3 blue. These computations are not
 * stored.
 * If the selection is not null or empty, computed data is not stored and seq
 * is not used.
 * If seq is null (single image processing), image_index is ignored, data is
 * stored in the fit, which cannot be NULL.
 * If seq is non-null, fit can be null to check for cached data.
 * The return value, if non-null, may be freed only using the free_stats()
 * function because of the special rule of this object that has a reference
 * counter because it can be referenced in 3 different places.
 */
imstats* statistics(sequence *seq, int image_index, fits *fit, int super_layer, rectangle *selection, int option, threading_type threads) {
	imstats *oldstat = NULL, *stat;
	check_threading(&threads);
	int layer = abs(super_layer);
	if (selection && selection->h > 0 && selection->w > 0) {
		// we have a selection, don't store anything
		if (!fit) return NULL;
		return statistics_internal(fit, super_layer, selection, option, NULL, fit->bitpix, threads);
	} else if (super_layer < 0) {
		// we are computing stats per filter on a CFA image, don't store anything
		if (!fit) return NULL;
		return statistics_internal(fit, super_layer, NULL, option, NULL, fit->bitpix, threads);
	} else if (!seq || image_index < 0) {
		// we have a single image, store in the fits
		g_assert(layer < fit->naxes[2]);
		if (fit->stats && fit->stats[layer]) {
			oldstat = fit->stats[layer];
			g_atomic_int_inc(&oldstat->_nb_refs);
		}
		stat = statistics_internal(fit, layer, NULL, option, oldstat, fit->bitpix, threads);
		if (!stat) {
			fprintf(stderr, "- stats failed for fit %p (%d)\n", fit, layer);
			if (oldstat) {
				stats_set_default_values(oldstat);
				g_atomic_int_dec_and_test(&oldstat->_nb_refs);
			}
			return NULL;
		}
		if (!oldstat)
			add_stats_to_fit(fit, layer, stat);
		return stat;
	} else {
		// we have sequence data, store in the sequence
		if (seq->stats && seq->stats[layer]) {
			oldstat = seq->stats[layer][image_index];
			if (oldstat)	// can be NULL here
				g_atomic_int_inc(&oldstat->_nb_refs);
		}
		stat = statistics_internal(fit, layer, NULL, option, oldstat, seq->bitpix, threads);
		if (!stat) {
			if (fit)
				fprintf(stderr, "- stats failed for %d in seq (%d)\n",
						image_index, layer);
			if (oldstat) {
				stats_set_default_values(oldstat);
				g_atomic_int_dec_and_test(&oldstat->_nb_refs);
			}
			return NULL;
		}
		if (!oldstat)
			add_stats_to_seq(seq, image_index, layer, stat);
		if (fit)
			add_stats_to_fit(fit, layer, stat);	// can be useful too
		return stat;
	}
}

int compute_means_from_flat_cfa_ushort(fits *fit, double mean[36]) {
	int row, col, c, i[36] = {0};
	WORD *data;
	unsigned int width, height;
	unsigned int startx, starty;

	data = fit->data;
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
			mean[(col % 6) + (row % 6) * 6] += (double) data[col + row * width];
			i[(col % 6) + (row % 6) * 6]++;
		}
	}

	for (c = 0; c < 36; c++) {
		mean[c] /= (double) i[c];
	}
	return 0;
}

/*  the mean is computed by adding data over and over in a long loop;
    but if the mean variable is a single-precision float, this results in
    the accumulation of rounding errors. */
int compute_means_from_flat_cfa(fits *fit, double mean[36]) {
	if (fit->type == DATA_USHORT)
		return compute_means_from_flat_cfa_ushort(fit, mean);
	if (fit->type == DATA_FLOAT)
		return compute_means_from_flat_cfa_float(fit, mean);
	return -1;
}

/****************** statistics caching and data management *****************/

/* reference an imstats struct to a fits, creates the stats array if needed */
void add_stats_to_fit(fits *fit, int layer, imstats *stat) {
	if (!fit->stats) {
		fit->stats = calloc(fit->naxes[2], sizeof(imstats *));
		if (!fit->stats) {
			PRINT_ALLOC_ERR;
			return;
		}
	}
	if (fit->stats[layer]) {
		if (fit->stats[layer] != stat) {
			siril_debug_print("- stats %p in fit %p (%d) is being replaced\n", fit->stats[layer], fit, layer);
			free_stats(fit->stats[layer]);
		} else return;
	}
	fit->stats[layer] = stat;
	g_atomic_int_inc(&stat->_nb_refs);
	siril_debug_print("- stats %p saved to fit %p (%d)\n", stat, fit, layer);
}

static void add_stats_to_stats(sequence *seq, int nb_layers, imstats ****stats, int image_index, int layer, imstats *stat) {
	if (!*stats) {
		*stats = calloc(nb_layers, sizeof(imstats **));
		if (!*stats) {
			PRINT_ALLOC_ERR;
			return;
		}
	}
	if (!(*stats)[layer]) {
		(*stats)[layer] = calloc(seq->number, sizeof(imstats *));
		if (!(*stats)[layer]) {
			PRINT_ALLOC_ERR;
			return;
		}
	}

	if ((*stats)[layer][image_index]) {
		if ((*stats)[layer][image_index] != stat) {
			siril_debug_print("- stats %p, %d in seq (%d) is being replaced\n", (*stats)[layer][image_index], image_index, layer);
			free_stats((*stats)[layer][image_index]);
		} else return;
	}
	siril_debug_print("- stats %p, %d in seq (%d): saving data\n", stat, image_index, layer);
	(*stats)[layer][image_index] = stat;
	seq->needs_saving = TRUE;
	g_atomic_int_inc(&stat->_nb_refs);
}

void add_stats_to_seq(sequence *seq, int image_index, int layer, imstats *stat) {
	add_stats_to_stats(seq, seq->nb_layers, &seq->stats, image_index, layer, stat);
}

void add_stats_to_seq_backup(sequence *seq, int image_index, int layer, imstats *stat) {
	add_stats_to_stats(seq, 3, &seq->stats_bkp, image_index, layer, stat);
}

/* saves cached stats from the fits to its sequence, and clears the cache of the fits */
void save_stats_from_fit(fits *fit, sequence *seq, int index) {
	int layer;
	if (!fit || !fit->stats || !seq || index < 0) return;
	for (layer = 0; layer < fit->naxes[2]; layer++) {
		if (fit->stats[layer])
			add_stats_to_seq(seq, index, layer, fit->stats[layer]);
		free_stats(fit->stats[layer]);
		fit->stats[layer] = NULL;
	}
}

/* fit must be already read from disk or have naxes set at least */
void copy_seq_stats_to_fit(sequence *seq, int index, fits *fit) {
	if (seq->stats) {
		int layer;
		for (layer = 0; layer < fit->naxes[2]; layer++) {
			if (seq->stats[layer] && seq->stats[layer][index]) {
				add_stats_to_fit(fit, layer, seq->stats[layer][index]);
				siril_debug_print("- stats %p, copied from seq %d (%d)\n", fit->stats[layer], index, layer);
			}
		}
	}
}

/* if image data has changed, use this to force recomputation of the stats */
void invalidate_stats_from_fit(fits *fit) {
	if (fit->stats) {
		int layer;
		for (layer = 0; layer < fit->naxes[2]; layer++) {
			siril_debug_print("- stats %p cleared from fit (%d)\n", fit->stats[layer], layer);
			free_stats(fit->stats[layer]);
			fit->stats[layer] = NULL;
		}
	}
	fit->maxi = -1;
}

/* if image data and image structure has changed, invalidate the complete stats data structure */
void full_stats_invalidation_from_fit(fits *fit) {
	if (fit->stats) {
		invalidate_stats_from_fit(fit);
		free(fit->stats);
		fit->stats = NULL;
	}
}

static void stats_set_default_values(imstats *stat) {
	stat->total = -1L;
	stat->ngoodpix = -1L;
	stat->mean = stat->avgDev = stat->median = stat->sigma = stat->bgnoise = stat->min = stat->max = stat->normValue = stat->mad = stat->sqrtbwmv = stat->location = stat->scale = NULL_STATS;
}

/* allocates an imstat structure and initializes it with default values that
 * are used by the statistics() function.
 * Only use free_stats() to free the return value.
 * Increment the _nb_refs if a new reference to the struct's address is made. */
void allocate_stats(imstats **stat) {
	if (stat) {
		if (!*stat)
			*stat = malloc(sizeof(imstats));
		if (!*stat) { PRINT_ALLOC_ERR; return; } // OOM
		stats_set_default_values(*stat);
		(*stat)->_nb_refs = 1;
		siril_debug_print("- stats %p allocated\n", *stat);
	}
}

/* frees an imstats struct if there are no more references to it.
 * returns NULL if it was freed, the argument otherwise. */
imstats* free_stats(imstats *stat) {
	if (stat && (g_atomic_int_dec_and_test(&stat->_nb_refs)) == TRUE) {
		siril_debug_print("- stats %p has no more refs, freed\n", stat);
		free(stat);
		return NULL;
	}
	if (stat)
		siril_debug_print("- stats %p has refs (%d)\n", stat, n);
	return stat;
}

/* calls free_stats on all stats of a sequence */
void clear_stats(sequence *seq, int layer) {
	if (seq->stats && seq->stats[layer]) {
		int i;
		for (i = 0; i < seq->number; i++) {
			if (seq->stats[layer][i]) {
				siril_debug_print("- stats %p freed from seq %d (%d)\n", seq->stats[layer][i], i, layer);
				free_stats(seq->stats[layer][i]);
				seq->stats[layer][i] = NULL;
			}
		}
	}
}

/* calls free_stats on all stats of a sequence */
void clear_stats_bkp(sequence *seq, int layer) {
	if (seq->stats_bkp && seq->stats_bkp[layer]) {
		int i;
		for (i = 0; i < seq->number; i++) {
			if (seq->stats_bkp[layer][i]) {
				siril_debug_print("- stats %p freed from seq %d (%d)\n", seq->stats_bkp[layer][i], i, layer);
				free_stats(seq->stats_bkp[layer][i]);
				seq->stats_bkp[layer][i] = NULL;
			}
		}
	}
}

/** generic function for sequences */

static void free_stat_list(gchar **list, int nb) {
	for (int i = 0; i < nb; i++) {
		g_free(list[i]);
	}
	free(list);
}

static int stat_prepare_hook(struct generic_seq_args *args) {
	struct stat_data *s_args = (struct stat_data*) args->user;
	if (s_args->option != STATS_BASIC && s_args->option != STATS_MAIN &&
			s_args->option != (STATS_NORM | STATS_MAIN)) {
		siril_log_color_message(_("Bad argument to stats option\n"), "red");
		return 1;
	}
	int nb_layers = s_args->cfa ? 3 : s_args->seq->nb_layers;
	// cfa may be set to TRUE for a non CFA sequence, but we don't know yet, so we still alloc for 3
	s_args->list = calloc(args->nb_filtered_images * nb_layers, sizeof(char*));
	return 0;
}

static int stat_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct stat_data *s_args = (struct stat_data*) args->user;

	gboolean is_cfa = fit->keywords.bayer_pattern[0] != '\0' && s_args->cfa;
	int nb_image_layers = (int)fit->naxes[2];
	if (is_cfa)
		nb_image_layers = 3;
	int nb_data_layers = s_args->cfa ? 3 : s_args->seq->nb_layers;

	for (int layer = 0; layer < nb_image_layers; layer++) {
		/* we first check for data in cache */
		int super_layer = is_cfa ? -layer - 1 : layer;
		imstats* stat = statistics(args->seq, i, NULL, super_layer, &s_args->selection, s_args->option, SINGLE_THREADED);
		if (!stat) {
			/* if no cache */
			stat = statistics(args->seq, i, fit, super_layer, &s_args->selection, s_args->option, SINGLE_THREADED);
			if (!stat) {
				siril_log_message(_("Error: statistics computation failed.\n"));
				return 1;
			}
		}

		int new_index = o * nb_data_layers;
		if (s_args->option == STATS_BASIC) {
			s_args->list[new_index + layer] = g_strdup_printf("%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n",
					i + 1,
					layer,
					stat->mean,
					stat->median,
					stat->sigma,
					stat->min,
					stat->max,
					stat->bgnoise
			);
		} else if (s_args->option == STATS_MAIN){
			s_args->list[new_index + layer] = g_strdup_printf("%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
					i + 1,
					layer,
					stat->mean,
					stat->median,
					stat->sigma,
					stat->min,
					stat->max,
					stat->bgnoise,
					stat->avgDev,
					stat->mad,
					stat->sqrtbwmv
			);
		} else {
			s_args->list[new_index + layer] = g_strdup_printf("%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
					i + 1,
					layer,
					stat->mean,
					stat->median,
					stat->sigma,
					stat->min,
					stat->max,
					stat->bgnoise,
					stat->avgDev,
					stat->mad,
					stat->sqrtbwmv,
					stat->location,
					stat->scale
			);
		}
		if (free_stats(stat))
			siril_debug_print("Error freeing stats in seqstat...\n");
	}
	return 0;
}

static int stat_finalize_hook(struct generic_seq_args *args) {
	GError *error = NULL;
	struct stat_data *s_args = (struct stat_data*) args->user;
	if (!s_args->list) {
		free(s_args);
		return 1;
	}

	int nb_data_layers = s_args->cfa ? 3 : s_args->seq->nb_layers;
	int size = nb_data_layers * args->nb_filtered_images;
	GFile *file = g_file_new_for_path(s_args->csv_name);
	GOutputStream* output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
	g_free(s_args->csv_name);
	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "Cannot save histo\n");
		}
		g_object_unref(file);
		free_stat_list(s_args->list, size);
		free(s_args);
		return 1;
	}
	const gchar *header;
	if (s_args->option == STATS_BASIC) {
		header = "image\tchan\tmean\tmedian\tsigma\tmin\tmax\tnoise\n";
	} else if (s_args->option == STATS_MAIN){
		header = "image\tchan\tmean\tmedian\tsigma\tmin\tmax\tnoise\tavgDev\tmad\tsqrtbwmv\n";
	} else {
		header = "image\tchan\tmean\tmedian\tsigma\tmin\tmax\tnoise\tavgDev\tmad\tsqrtbwmv\tlocation\tscale\n";
	}
	if (!g_output_stream_write_all(output_stream, header, strlen(header), NULL, NULL, &error)) {
		g_warning("%s\n", error->message);
		g_clear_error(&error);
		g_object_unref(output_stream);
		g_object_unref(file);
		free_stat_list(s_args->list, size);
		free(s_args);
		return 1;
	}

	for (int i = 0; i < args->nb_filtered_images * nb_data_layers; i++) {
		if (!s_args->list[i]) continue; //stats can fail
		if (!g_output_stream_write_all(output_stream, s_args->list[i], strlen(s_args->list[i]), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			free_stat_list(s_args->list, size);
			free(s_args);
			return 1;
		}
	}

	siril_log_message(_("Statistic file %s was successfully created.\n"), g_file_peek_path(file));
	writeseqfile(args->seq);
	g_object_unref(output_stream);
	g_object_unref(file);
	free_stat_list(s_args->list, size);
	free(s_args);

	return 0;
}

static int stat_compute_mem_limit(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 0.0, FALSE, &MB_per_image, NULL, &MB_avail); // no output - output_scale is set to 0.

	int is_color = args->seq->nb_layers == 3;
	required = is_color ? MB_per_image * 4 / 3 : MB_per_image * 2;

	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
                        thread_limit = com.max_thread;
		limit = thread_limit;

		/* should not happen */
		if (for_writer)
			return 1;
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
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
#endif
	}
	return limit;
}


void apply_stats_to_sequence(struct stat_data *stat_args) {
	struct generic_seq_args *args = create_default_seqargs(stat_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = stat_args->seq->selnum;
	args->compute_mem_limits_hook = stat_compute_mem_limit;
	args->prepare_hook = stat_prepare_hook;
	args->finalize_hook = stat_finalize_hook;
	args->image_hook = stat_image_hook;
	args->description = _("Statistics");
	args->has_output = FALSE;
	args->new_seq_prefix = NULL;
	args->user = stat_args;

	stat_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		int nb_data_layers = stat_args->cfa ? 3 : stat_args->seq->nb_layers;
		int size = nb_data_layers * args->nb_filtered_images;
		free_stat_list(stat_args->list, size);
		free (stat_args);
		free_generic_seq_args(args, TRUE);
	}
}

/**** callbacks ****/

void on_menu_gray_stat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	set_cursor_waiting(FALSE);
}

/**** threading helpers ****/

/* compute statistics for all channels of an image from a sequence, only on full image, make sure the result (stats) is allocated */
int compute_all_channels_statistics_seqimage(sequence *seq, int image_index, fits *fit, int option,
		threading_type threading, int image_thread_id, imstats **stats) {
	int required_computations = 0;
	int stats_to_compute[3];
	int nb_layers = seq->nb_layers;
	int retval = 0;
	fits local_fit = { 0 };

	// try from the cache first to see if image loading is required
	for (int layer = 0; layer < nb_layers; ++layer) {
		// try with no fit passed: fails if data is needed because data is not cached
		stats[layer] = statistics(seq, image_index, NULL, layer, NULL, option, SINGLE_THREADED);
		if (!stats[layer])
			stats_to_compute[required_computations++] = layer;
	}
	siril_debug_print("%d stats to recompute for image %d\n", required_computations, image_index);

	if (required_computations) {
		if (!fit) {
			if (seq_read_frame(seq, image_index, &local_fit, TRUE, image_thread_id)) {
				return -1;
			}

			if (nb_layers != (int)local_fit.naxes[2]) {
				siril_log_color_message(_("It looks like your sequence contains a mix of monochrome and RGB images.\n"), "red");
				clearfits(&local_fit);
				return -1;
			}
			fit = &local_fit;
		}

#ifdef _OPENMP
		/* we have 'threading' threads and 'required_computations' channels to
		 * compute stats on. The most efficient approach is to process the
		 * channels statistics in parallel, because the statistics of each channel
		 * can be a small operation that would suffer from the threading overhead.
		 * Consequently, if we have enough threads available, we parallelize first
		 * the channels and if there are more than enough, we parallelize the
		 * stats computation. The stats computation then decides if it will use
		 * all the threads it was given or not, depending on data size.
		 */
		int threads = check_threading(&threading);
		int channels_per_thread = required_computations;
		int *threads_per_thread = NULL;
		if (threads > 1) {
			if (threads >= required_computations) {
				channels_per_thread = 1;
				threads_per_thread = compute_thread_distribution(required_computations, threading);
				omp_set_max_active_levels(INT_MAX);	// to be done at each level
				if ((omp_get_max_active_levels() < 2) && threads_per_thread[0] > 1)
					siril_log_message(_("Threading statistics computation per channel, but enabling threading in channels too is not supported.\n"));
				else siril_debug_print("threading statistics computation per channel (at most %d threads each)\n", threads_per_thread[0]);
			}
			else {
				threads_per_thread = compute_thread_distribution(1, threading);
				siril_debug_print("threading statistics computation of channels (%d threads)\n", threads_per_thread[0]);
			}
		}

#pragma omp parallel for num_threads(required_computations) schedule(static) if(required_computations > 1 && channels_per_thread == 1)
#endif
		for (int layer_index = 0; layer_index < required_computations; ++layer_index) {
			threading_type subthreads = SINGLE_THREADED;
			int layer = stats_to_compute[layer_index];
#ifdef _OPENMP
			//siril_debug_print("actual number of threads in channel thread: %d (level %d)\n", omp_get_num_threads(), omp_get_level());
			if (threads_per_thread) {
				int thread_id = omp_get_thread_num();
				subthreads = threads_per_thread[thread_id];
			}
			siril_debug_print("requesting stats for normalization of channel %d with %d threads\n", layer, subthreads);
#endif
			stats[layer] = statistics(seq, image_index, fit, layer, NULL, option, subthreads);
			if (!stats[layer]) {
				retval = -1;
				continue;
			}
		}

#ifdef _OPENMP
		free(threads_per_thread);
#endif

		if (seq->type != SEQ_INTERNAL && fit == &local_fit)
			clearfits(&local_fit);
	}
	return retval;
}

/* compute statistics for all channels of a single image, only on full image, make sure the result (stats) is allocated */
int compute_all_channels_statistics_single_image(fits *fit, int option,
		threading_type threading, imstats **stats) {
	g_assert(fit);
	gboolean cfa = (option & STATS_FOR_CFA) && fit->keywords.bayer_pattern[0] != '\0';
	int required_computations = cfa ? 3 : (int)fit->naxes[2];
	if (required_computations == 0) {
		stats[0] = NULL;
		return -1;
	}
	int retval = 0;

#ifdef _OPENMP
	/* we have 'threading' threads and 'required_computations' channels to
	 * compute stats on. The most efficient approach is to process the
	 * channels statistics in parallel, because the statistics of each channel
	 * can be a small operation that would suffer from the threading overhead.
	 * Consequently, if we have enough threads available, we parallelize first
	 * the channels and if there are more than enough, we parallelize the
	 * stats computation. The stats computation then decides if it will use
	 * all the threads it was given or not, depending on data size.
	 */
	int threads = check_threading(&threading);
	int channels_per_thread = required_computations;
	int *threads_per_thread = NULL;
	if (threads > 1) {
		if (threads >= required_computations) {
			channels_per_thread = 1;
			threads_per_thread = compute_thread_distribution(required_computations, threading);
			omp_set_max_active_levels(INT_MAX);	// to be done at each level
			if ((omp_get_max_active_levels() < 2) && threads_per_thread[0] > 1)
				siril_log_message(_("Threading statistics computation per channel, but enabling threading in channels too is not supported.\n"));
			else siril_debug_print("threading statistics computation per channel (at most %d threads each)\n", threads_per_thread[0]);
		}
		else {
			threads_per_thread = compute_thread_distribution(1, threading);
			siril_debug_print("threading statistics computation of channels (%d threads)\n", threads_per_thread[0]);
		}
	}

#pragma omp parallel for num_threads(required_computations) schedule(static) if(required_computations > 1 && channels_per_thread == 1)
#endif
	for (int layer = 0; layer < required_computations; ++layer) {
		threading_type subthreads = SINGLE_THREADED;
#ifdef _OPENMP
		//siril_debug_print("actual number of threads in channel thread: %d (level %d)\n", omp_get_num_threads(), omp_get_level());
		if (threads_per_thread) {
			int thread_id = omp_get_thread_num();
			subthreads = threads_per_thread[thread_id];
		}
		siril_debug_print("requesting stats for normalization of channel %d with %d threads\n", layer, subthreads);
#endif
		int layer_arg = cfa ? -layer - 1 : layer;
		stats[layer] = statistics(NULL, -1, fit, layer_arg, NULL, option, subthreads);
		if (!stats[layer])
			retval = -1;
	}

#ifdef _OPENMP
	free(threads_per_thread);
#endif

	return retval;
}

/* get a reference to stats of an image, only if already computed
 * make sure channels is allocated to the number of channels
 * use free_stats() to release the returned stats
 */
int copy_cached_stats_for_image(sequence *seq, int image, imstats **channels) {
	if (!seq->stats)
		return 1;
	int all_copied = 1;
	for (int i = 0; i < seq->nb_layers; i++) {
		if (all_copied && seq->stats[i] && seq->stats[i][image]) {
			imstats *stats = seq->stats[i][image];
			g_atomic_int_inc(&stats->_nb_refs);
			channels[i] = stats;
		} else {
			all_copied = 0;
		}
	}
	if (!all_copied) {
		for (int i = 0; i < seq->nb_layers; i++) {
			if (channels[i])
				free_stats(channels[i]);
			channels[i] = NULL;
		}
	}
	return !all_copied;
}

int sos_update_noise_float(float *array, long nx, long ny, long nchans, double *noise) {
	int status, retval = 0;
	float *colarray[3];
	double fSigma = 0.0;
	if (nchans == 1) {
		retval = siril_fits_img_stats_float(array, nx, ny, NULL, NULL, NULL,
		NULL, NULL, noise, NULL, NULL, NULL, MULTI_THREADED, &status);
		return retval;
	} else {
		colarray[0] = array;
		colarray[1] = array + (nx * ny);
		colarray[2] = array + 2 * (nx * ny);
		for (unsigned i = 0 ; i < nchans ; i++) {
			retval += siril_fits_img_stats_float(colarray[i], nx, ny, NULL, NULL, NULL,
						NULL, NULL, noise, NULL, NULL, NULL, MULTI_THREADED, &status);
			fSigma += *noise;
		}
		*noise = fSigma / nchans;
	}
	return retval;
}

double robust_median_w(fits *fit, rectangle *area, int chan, float lower, float upper) {
	uint32_t x0, y0, x1,y1;
	if (area) {
		x0 = area->x;
		y0 = area->y;
		x1 = area->x + area->w;
		y1 = area->y + area->h;
	} else {
		x0 = y0 = 0;
		x1 = fit->rx;
		y1 = fit->ry;
	}
	size_t npixels = (x1 - x0) * (y1 - y0);
	WORD *data = fit->pdata[chan];
	WORD *filtered_data = malloc(npixels * sizeof(WORD));
	size_t count = 0;
	WORD lowerw = roundf_to_WORD(lower);
	WORD upperw = roundf_to_WORD(upper);
	for (uint32_t y = y0 ; y < y1 ; y++) {
		size_t j = y * fit->rx;
		for (uint32_t x = x0 ; x < x1 ; x++) {
			size_t i = x + j;
			if (data[i] >= lowerw && data[i] <= upperw) {
				filtered_data[count++] = data[i];
			}
		}
	}

	// Check if there are any elements in the specified range
	if (count == 0) {
		free(filtered_data);
		return 0.0; // No elements in the range, return 0 as median
	}

	// use histogram_median here instead of quickmedian for speed (see #1458)
	double retval = histogram_median(filtered_data, count, MULTI_THREADED);

	// Free the allocated memory for filtered_data
	free(filtered_data);

	return retval;
}

// Function to quickly compute min and max values
int quick_minmax(fits *fit, double *minval, double *maxval) {
    imstats *stats[3] = { NULL };
    int retval = compute_all_channels_statistics_single_image(fit, STATS_MINMAX, MULTI_THREADED, stats);

    if (retval) {
        siril_log_color_message(_("Error: statistics computation failed. Unable to check for out-of-range values.\n"), "red");
    } else {
        if (fit->naxes[2] == 1) {
            *maxval = stats[0]->max;
            *minval = stats[0]->min;
        } else {
            *maxval = max(max(stats[RLAYER]->max, stats[GLAYER]->max), stats[BLAYER]->max);
            *minval = min(min(stats[RLAYER]->min, stats[GLAYER]->min), stats[BLAYER]->min);
        }
        for (int i = 0; i < fit->naxes[2]; i++) {
            free_stats(stats[i]);
        }
    }
    return retval;
}
