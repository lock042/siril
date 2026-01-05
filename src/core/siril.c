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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/arithm.h"
#include "core/siril_log.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "algos/demosaicing.h"
#include "opencv/opencv.h"
#include "rt/gauss.h"

int threshlo(fits *fit, WORD level) {
	size_t i, n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->data;
		for (i = 0; i < n; ++i) {
			buf[i] = max(level, buf[i]);
		}
	} else if (fit->type == DATA_FLOAT) {
		float l = (float) level / USHRT_MAX_SINGLE;
		float *buf = fit->fdata;
		for (i = 0; i < n; ++i) {
			buf[i] = max(l, buf[i]);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int threshhi(fits *fit, WORD level) {
	size_t i, n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->data;
		for (i = 0; i < n; ++i) {
			buf[i] = min(level, buf[i]);
		}
	} else if (fit->type == DATA_FLOAT) {
		float l = (float) level / USHRT_MAX_SINGLE;
		float *buf = fit->fdata;
		for (i = 0; i < n; ++i) {
			buf[i] = min(l, buf[i]);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

// level is for ushort data, adapted automatically in case of float data
int nozero(fits *fit, WORD level) {
	size_t i, n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->data;
		for (i = 0; i < n; ++i) {
			if (buf[i] == 0)
				buf[i] = level;
		}
	} else if (fit->type == DATA_FLOAT) {
		float l = (float) level / USHRT_MAX_SINGLE;
		float *buf = fit->fdata;
		for (i = 0; i < n; ++i) {
			if (buf[i] <= 0.0)
				buf[i] = l;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

// in-place Gaussian blur with RawTherapee's implementation (SSE/vectorization and OpenMP)
// only implemented for float and monochrome images, as support for other one-channel algorithms
int gaussian_blur_RT(fits *fit, double sigma, int threads) {
	g_assert(fit->naxes[2] == 1);
	if (fit->type == DATA_FLOAT) {
		siril_debug_print("Using RawTherapee in-place Gaussian blur with sigma=%f and %d threads\n", sigma, threads);
		// RawTherapee gaussianBlur (mono only)
		int rx = (int)fit->naxes[0];
		int ry = (int)fit->naxes[1];
		float **src = malloc(ry * sizeof(float *));
		if (!src) { PRINT_ALLOC_ERR; return 1; }
		for (int k = 0; k < ry; k++) {
			src[k] = fit->fdata + k * rx;
		}

		gaussianBlurC(src, src, rx, ry, sigma, threads);
		free(src);
		return 0;
	}
	else {
		// OpenCV GaussianBlur
		return cvUnsharpFilter(fit, sigma, 1.0);
	}
}

// using a temporary buffer
int gaussian_blur_RT2(fits *fit, double sigma, int threads) {
	if (fit->type == DATA_FLOAT) {
		siril_debug_print("Using RawTherapee out-of-place Gaussian blur with sigma=%f and %d threads\n", sigma, threads);
		// RawTherapee gaussianBlur (mono only)
		size_t n = fit->naxes[0] * fit->naxes[1];
		int rx = (int)fit->naxes[0];
		int ry = (int)fit->naxes[1];
		float *result = malloc(n * sizeof(float));
		if (!result) { PRINT_ALLOC_ERR; return 1; }
		float **src = malloc(ry * sizeof(float *));
		if (!src) { PRINT_ALLOC_ERR; free(result); return 1; }
		float **dst = malloc(ry * sizeof(float *));
		if (!dst) { PRINT_ALLOC_ERR; free(src); free(result); return 1; }
		for (int k = 0; k < ry; k++) {
			src[k] = fit->fdata + k * rx;
			dst[k] = result + k * rx;
		}

		gaussianBlurC(src, dst, rx, ry, sigma, threads);
		free(src);
		free(dst);
		float *olddata = gfit->fdata;
		gfit->fdata = result;
		gfit->fpdata[RLAYER] = gfit->fdata;
		if (gfit->naxis == 3) {
			gfit->fpdata[GLAYER] = gfit->fdata + n;
			gfit->fpdata[BLAYER] = gfit->fdata + n * 2;
		} else {
			gfit->fpdata[GLAYER] = gfit->fdata;
			gfit->fpdata[BLAYER] = gfit->fdata;
		}
		free(olddata);
		return 0;
	}
	else {
		// OpenCV GaussianBlur
		return cvUnsharpFilter(fit, sigma, 1.0);
	}
}

int unsharp(fits *fit, double sigma, double amount, gboolean verbose) {
	struct timeval t_start, t_end;

	if (sigma <= 0.0)
		return 1;
	if (verbose) {
		siril_log_color_message(_("Unsharp: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	cvUnsharpFilter(fit, sigma, amount);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	return 0;
}

/* This entropy function computes the entropy for the image in gfit for its
 * layer 'layer', in the area designated by area which can be NULL.
 * An optional imstats parameter can be used to provide the background and
 * sigma value, and when it is given, the entropy will only be computed for
 * pixels with values above background + 1 * sigma. It must be NULL otherwise.
 */
float entropy(fits *fit, int layer, rectangle *area, const imstats *opt_stats) {
	float e = 0.f;
	gsl_histogram *histo;

	if (area == NULL)
		histo = computeHisto(fit, layer);
	else
		histo = computeHisto_Selection(fit, layer, area);

	size_t n = fit->naxes[0] * fit->naxes[1];
	g_assert(n > 0);
	size_t size = gsl_histogram_bins(histo);

	for (size_t i = 0; i < size; i++) {
		double p = gsl_histogram_get(histo, i) / n;
		if (p > 0)
			e -= p * log(p);
	}

	gsl_histogram_free(histo);

	return e;
}

static int loglut_ushort(fits *fit) {
	// This function maps fit with a log LUT
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double norm = USHRT_MAX_DOUBLE / log(USHRT_MAX_DOUBLE);

	for (int layer = 0; layer < fit->naxes[2]; ++layer) {
		imstats *stat = statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX, MULTI_THREADED);
		double min = stat->min;
		double wd = stat->max - stat->min;
		size_t i, n = fit->naxes[0] * fit->naxes[1];
		for (i = 0; i < n; i++) {
			float px = (float)buf[layer][i];
			buf[layer][i] = round_to_WORD(log1pf((px - min) / wd) * norm);
		}
		free_stats(stat);
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int loglut_float(fits *fit) {
	// This function maps fit with a log LUT
	float *buf[3] = { fit->fpdata[RLAYER],
			fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	for (int layer = 0; layer < fit->naxes[2]; ++layer) {
		imstats *stat = statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX, MULTI_THREADED);
		double min = stat->min;
		double wd = stat->max - stat->min;
		size_t i, n = fit->naxes[0] * fit->naxes[1];
		for (i = 0; i < n; i++) {
			float px = buf[layer][i];
			buf[layer][i] = log1pf((px - min) / wd);
		}
		free_stats(stat);
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int loglut(fits *fit) {
	if (fit->type == DATA_USHORT) {
		return loglut_ushort(fit);
	} else if (fit->type == DATA_FLOAT) {
		return loglut_float(fit);
	}
	return -1;
}

int ddp(fits *a, float level, float coeff, float sigma) {
	fits fit = { 0 };
	if (a->orig_bitpix == BYTE_IMG) {
		siril_log_color_message(_("This process cannot be applied to 8b images\n"), "red");
		return 1;
	}
	if (level < 0.f || level > USHRT_MAX_SINGLE) {
		siril_log_color_message(_("ddp level argument must be [0, 65535]\n"), "green");
		return 1;
	}
	if (level < 1.f && a->type == DATA_FLOAT)
		level *= USHRT_MAX_SINGLE;
	float l = ushort_to_float_range((WORD) level);

	int ret = copyfits(a, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (!ret) ret = unsharp(&fit, sigma, 0, FALSE);
	if (!ret) ret = soper(&fit, l, OPER_ADD, TRUE);
	if (!ret) ret = nozero(&fit, 1);
	if (!ret) ret = siril_fdiv(a, &fit, l, TRUE);
	if (!ret) ret = soper(a, coeff, OPER_MUL, TRUE);
	clearfits(&fit);
	invalidate_stats_from_fit(a);
	if (!ret) {
		char log[90];
		sprintf(log, "DDP stretch, threshold: %.2f, multiplier: %.2f, sigma: %.1f", level, coeff, sigma);
		a->history = g_slist_append(a->history, strdup(log));
	}
	return ret;
}

int visu(fits *fit, int low, int high) {
	if (low < 0 || low > USHRT_MAX || high < 1 || high > USHRT_MAX)
		return 1;
	if (!single_image_is_loaded() && !sequence_is_loaded())
		return 1;
	gui.lo = low;
	gui.hi = high;
	set_cutoff_sliders_values();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	return 0;
}

/* fill an image or selection with the value 'level' */
int fill(fits *fit, int level, const rectangle *arearg) {
	rectangle area;

	if (arearg) {
		memcpy(&area, arearg, sizeof(rectangle));
	} else {
		if (com.selection.h && com.selection.w) {
			memcpy(&area, &com.selection, sizeof(rectangle));
		} else {
			area.w = fit->rx;
			area.h = fit->ry;
			area.x = 0;
			area.y = 0;
		}
	}
	for (int layer = 0; layer < fit->naxes[2]; ++layer) {
		if (fit->type == DATA_USHORT) {
			int maxlevel = (fit->orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
			if ((level > maxlevel) || (level < 0)) {
				siril_log_message(_("Fill value must be in the range [0,%d]\n"), maxlevel);
				return 1;
			}
			WORD *buf = fit->pdata[layer]
					+ (fit->ry - area.y - area.h) * fit->rx + area.x;
			int stridebuf = fit->rx - area.w;
			for (int i = 0; i < area.h; ++i) {
				for (int j = 0; j < area.w; ++j) {
					*buf++ = level;
				}
				buf += stridebuf;
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer]
					+ (fit->ry - area.y - area.h) * fit->rx + area.x;
			int stridebuf = fit->rx - area.w;
			float flevel = level * INV_USHRT_MAX_SINGLE;
			for (int i = 0; i < area.h; ++i) {
				for (int j = 0; j < area.w; ++j) {
					*buf++ = flevel;
				}
				buf += stridebuf;
			}
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int off_uchar(fits *fit, float level) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER],
			fit->pdata[BLAYER] };
	g_assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -UCHAR_MAX_SINGLE)
		level = -UCHAR_MAX_SINGLE;
	else if (level > UCHAR_MAX_SINGLE)
		level = UCHAR_MAX_SINGLE;
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	for (i = 0; i < n; ++i) {
		for (int layer = 0; layer < fit->naxes[2]; ++layer) {
			float val = (float)buf[layer][i];
			buf[layer][i] = (WORD)roundf_to_BYTE(val + level);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int off_ushort(fits *fit, float level) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER],
			fit->pdata[BLAYER] };
	g_assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -USHRT_MAX_SINGLE)
		level = -USHRT_MAX_SINGLE;
	else if (level > USHRT_MAX_SINGLE)
		level = USHRT_MAX_SINGLE;
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	for (i = 0; i < n; ++i) {
		for (int layer = 0; layer < fit->naxes[2]; ++layer) {
			float val = (float)buf[layer][i];
			buf[layer][i] = roundf_to_WORD(val + level);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int off_float(fits *fit, float level) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER],
			fit->fpdata[BLAYER] };
	g_assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -1.f)
		level = -1.f;
	else if (level > 1.f)
		level = 1.f;
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	for (i = 0; i < n; ++i) {
		for (int layer = 0; layer < fit->naxes[2]; ++layer) {
			float val = buf[layer][i];
			buf[layer][i] = set_float_in_interval(val + level, 0.f, 1.f);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int off(fits *fit, float level) {
	if (fit->type == DATA_USHORT) {
		if (fit->orig_bitpix == BYTE_IMG) return off_uchar(fit, level);
		return off_ushort(fit, level);
	} else if (fit->type == DATA_FLOAT) {
		level /= USHRT_MAX_SINGLE;
		return off_float(fit, level);
	}
	return -1;
}

/* computes the background value using the histogram and/or median value.
 * The argument layer can be -1 for automatic setting (= green for RGB) */
double background(fits* fit, int reqlayer, rectangle *selection, threading_type threading) {
	int layer = RLAYER;

	if (reqlayer >= 0)
		layer = reqlayer;
	else if (isrgb(gfit))
		layer = GLAYER;		//GLAYER is better to evaluate background

	imstats* stat = statistics(NULL, -1, fit, layer, selection, STATS_BASIC, threading);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		return -1.0;
	}
	double bg = stat->median;
	free_stats(stat);
	return bg;
}

void compute_grey_flat(fits *fit) {
	float coeff1, coeff2;
	double mean[36] = {0};
	double green_var[num_filter_patterns];
	double ch_mean[num_filter_patterns][3];
	double delta;
	unsigned int cell, pat, pat_width, pat_cell, ch_n[3];
	unsigned int guessed_pat = 0;

	/* compute mean of each element in 6x6 blocks */
	compute_means_from_flat_cfa(fit, mean);

	/* compute coefficients */

	for (pat = 0; pat < num_filter_patterns; pat++) {
		green_var[pat] = 0;
		ch_mean[pat][0] = ch_mean[pat][1] = ch_mean[pat][2] = 0;
		ch_n[0] = ch_n[1] = ch_n[2] = 0;
		/* compute width of the (square) CFA pattern */
		/* added 0.1 in case the result of sqrt is something like 5.999999 and it gets casted to int as 5 */
		pat_width = (unsigned int) (sqrt(strlen(filter_pattern[pat]))+0.1);

		for (cell = 0; cell < 36; cell++) {
			/* convert 6 x 6 block coordinates to CFA pattern local coordinates */
			pat_cell = ( cell / 6 % pat_width) * pat_width + cell % pat_width;

			switch(filter_pattern[pat][pat_cell]) {
				case 'G':
					delta = mean[cell] - ch_mean[pat][1];
					ch_n[1]++;
					ch_mean[pat][1] += delta/ch_n[1];
					/* we use welford's algortihm to compute the variance */
					green_var[pat] += delta*(mean[cell]-ch_mean[pat][1]);
					break;
				case 'R':
					delta = mean[cell] - ch_mean[pat][0];
					ch_n[0]++;
					ch_mean[pat][0] += delta/ch_n[0];
					break;
				case 'B':
					delta = mean[cell] - ch_mean[pat][2];
					ch_n[2]++;
					ch_mean[pat][2] += delta/ch_n[2];
					break;
			}
		}
		/* we use bessel's correction, so n-1 */
		green_var[pat] /= (double) (ch_n[1]-1);
		if (green_var[pat] < green_var[guessed_pat]) {
			guessed_pat = pat;
		}
	}

	siril_debug_print("Guessed pattern: #%d (%s)\n", guessed_pat, filter_pattern[guessed_pat]);

	coeff1 = ch_mean[guessed_pat][0]/ch_mean[guessed_pat][1];
	coeff2 = ch_mean[guessed_pat][2]/ch_mean[guessed_pat][1];
	siril_debug_print("coeff1: %.5f, coeff2: %.5f\n", coeff1, coeff2);

	/* applies coefficients to cfa image */
	equalize_cfa_fit_with_coeffs(fit, coeff1, coeff2, filter_pattern[guessed_pat]);
}
