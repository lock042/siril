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

#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "core/arithm.h"
#include "algos/astrometry_solver.h"
#include "algos/demosaicing.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "core/processing.h"
#include "opencv/opencv.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"

#include "geometry.h"

/* this method rotates the image 180 degrees, useful after german mount flip.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
static void fits_rotate_pi_ushort(fits *fit) {
	int i, line, axis;
	WORD *line1 = NULL, *line2 = NULL, *src, *dst, swap;

	size_t line_size = fit->rx * sizeof(WORD);
	line1 = malloc(line_size);
	line2 = malloc(line_size);
	if (!line1 || !line2) {
		PRINT_ALLOC_ERR;
		if (line1) free(line1);
		if (line2) free(line2);
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(line1, src, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line1[i];
				line1[i] = line1[fit->rx - i - 1];
				line1[fit->rx - i - 1] = swap;
			}
			memcpy(line2, dst, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line2[i];
				line2[i] = line2[fit->rx - i - 1];
				line2[fit->rx - i - 1] = swap;
			}
			memcpy(src, line2, line_size);
			memcpy(dst, line1, line_size);
		}
		if (fit->ry & 1) {
			/* swap the middle line */
			src = fit->pdata[axis] + line * fit->rx;
			for (i = 0; i < fit->rx / 2; i++) {
				swap = src[i];
				src[i] = src[fit->rx - i - 1];
				src[fit->rx - i - 1] = swap;
			}
		}
	}
	free(line1);
	free(line2);
}

static void fits_rotate_pi_float(fits *fit) {
	int i, line, axis;
	float *line1 = NULL, *line2 = NULL, *src, *dst, swap;

	size_t line_size = fit->rx * sizeof(float);
	line1 = malloc(line_size);
	line2 = malloc(line_size);
	if (!line1 || !line2) {
		PRINT_ALLOC_ERR;
		if (line1) free(line1);
		if (line2) free(line2);
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->fpdata[axis] + line * fit->rx;
			dst = fit->fpdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(line1, src, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line1[i];
				line1[i] = line1[fit->rx - i - 1];
				line1[fit->rx - i - 1] = swap;
			}
			memcpy(line2, dst, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line2[i];
				line2[i] = line2[fit->rx - i - 1];
				line2[fit->rx - i - 1] = swap;
			}
			memcpy(src, line2, line_size);
			memcpy(dst, line1, line_size);
		}
		if (fit->ry & 1) {
			/* swap the middle line */
			src = fit->fpdata[axis] + line * fit->rx;
			for (i = 0; i < fit->rx / 2; i++) {
				swap = src[i];
				src[i] = src[fit->rx - i - 1];
				src[fit->rx - i - 1] = swap;
			}
		}
	}
	free(line1);
	free(line2);
}

static void fits_rotate_pi(fits *fit) {
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	if (fit->type == DATA_USHORT) {
		fits_rotate_pi_ushort(fit);
	} else if (fit->type == DATA_FLOAT) {
		fits_rotate_pi_float(fit);
	}
}

static void fit_update_buffer(fits *fit, void *newbuf, int width, int height, int bin_factor) {
	size_t nbdata = width * height;

	full_stats_invalidation_from_fit(fit);

	if (fit->type == DATA_USHORT) {
		if (fit->data)
			free(fit->data);
		fit->data = (WORD *)newbuf;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->naxes[2] == 3 ? fit->data + nbdata : fit->data;
		fit->pdata[BLAYER] = fit->naxes[2] == 3 ? fit->data + nbdata * 2 : fit->data;
	}
	else if (fit->type == DATA_FLOAT) {
		if (fit->fdata)
			free(fit->fdata);
		fit->fdata = (float *)newbuf;
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->naxes[2] == 3 ? fit->fdata + nbdata : fit->fdata;
		fit->fpdata[BLAYER] = fit->naxes[2] == 3 ? fit->fdata + nbdata * 2 : fit->fdata;
	}
	/* update size */
	fit->naxes[0] = width;
	fit->naxes[1] = height;
	fit->rx = width;
	fit->ry = height;

	if (fit->keywords.binning_x == 0 || fit->keywords.binning_x == 1) {
		fit->keywords.binning_x = bin_factor;
		fit->keywords.binning_y = bin_factor;
	} else {
		fit->keywords.binning_x *= bin_factor;
		fit->keywords.binning_y *= bin_factor;
	}
	if (!com.pref.binning_update) {
		fit->keywords.pixel_size_x *= bin_factor;
		fit->keywords.pixel_size_y *= bin_factor;
	}
}

static void fits_binning_float(fits *fit, int bin_factor, gboolean mean) {
	if (bin_factor == 0) // Bin 0 would be nonsensical.
		return;
	int width = fit->rx;
	int height = fit->ry;
	int new_width = width / bin_factor;
	int new_height = height / bin_factor;

	size_t npixels = new_width * new_height;

	float *newbuf = malloc(npixels * fit->naxes[2] * sizeof(float));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		const float *buf = fit->fdata + (width * height) * channel;

		long k = 0 + channel * npixels;
		for (int row = 0, nrow = 0; row < height - bin_factor + 1; row += bin_factor, nrow++) {
			for (int col = 0, ncol = 0; col < width - bin_factor + 1; col += bin_factor, ncol++) {
				int c = 0;
				newbuf[k] = 0;
				for (int i = 0; i < bin_factor; i++) {
					for (int j = 0; j < bin_factor; j++) {
						newbuf[k] += buf[i + col + (j + row) * width];
						c++;
					}
				}
				if (mean) newbuf[k] /= c;
				k++;
			}
		}
	}
	fit_update_buffer(fit, newbuf, new_width, new_height, bin_factor);
	return;
}

static void fits_binning_ushort(fits *fit, int bin_factor, gboolean mean) {
	if (bin_factor == 0) // Bin 0 would be nonsensical.
		return;
	int width = fit->rx;
	int height = fit->ry;
	int new_width = width / bin_factor;
	int new_height = height / bin_factor;

	size_t npixels = new_width * new_height;

	WORD *newbuf = malloc(npixels * fit->naxes[2] * sizeof(WORD));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		const WORD *buf = fit->data + (width * height) * channel;

		long k = 0 + channel * npixels;
		for (int row = 0, nrow = 0; row < height - bin_factor + 1; row += bin_factor, nrow++) {
			for (int col = 0, ncol = 0; col < width - bin_factor + 1; col += bin_factor, ncol++) {
				int c = 0;
				int tmp = 0;
				for (int i = 0; i < bin_factor; i++) {
					for (int j = 0; j < bin_factor; j++) {
						tmp += (buf[i + col + (j + row) * width]);
						c++;
					}
				}
				if (mean) tmp /= c;
				newbuf[k] = truncate_to_WORD(tmp);
				k++;
			}
		}
	}
	fit_update_buffer(fit, newbuf, new_width, new_height, bin_factor);
}

int fits_binning(fits *fit, int factor, gboolean mean) {
	struct timeval t_start, t_end;

	siril_log_color_message(_("Binning x%d: processing...\n"), "green", factor);
	gettimeofday(&t_start, NULL);
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	if (fit->type == DATA_USHORT) {
		fits_binning_ushort(fit, factor, mean);
	} else if (fit->type == DATA_FLOAT) {
		fits_binning_float(fit, factor, mean);
	}

	free_wcs(fit);
	reset_wcsdata(fit);
	refresh_annotations(TRUE);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	char log[90];
	sprintf(log, "Binned x%d (%s)", factor, mean ? "mean" : "sum");
	fit->history = g_slist_append(fit->history, strdup(log));

	siril_log_message(_("New image size: %dx%d pixels.\n"), fit->rx, fit->ry);

	return 0;
}

const char *interp_to_str(opencv_interpolation interpolation) {
	switch (interpolation) {
		case OPENCV_NEAREST:
			return _("Nearest-Neighbor");
		default:
		case OPENCV_LINEAR:
			return _("Bilinear");
		case OPENCV_AREA:
			return _("Pixel Area Relation");
		case OPENCV_CUBIC:
			return _("Bicubic");
		case OPENCV_LANCZOS4:
			return _("Lanczos4");
	}
}

/* These functions do not more than resize_gaussian and rotate_image
 * except for console outputs.
 * Indeed, siril_log_message seems not working in a cpp file */
int verbose_resize_gaussian(fits *image, int toX, int toY, opencv_interpolation interpolation, gboolean clamp) {
	int retvalue;
	struct timeval t_start, t_end;
	float factor_X = (float)image->rx / (float)toX;
	float factor_Y = (float)image->ry / (float)toY;

	siril_log_color_message(_("Resample (%s interpolation): processing...\n"),
			"green", interp_to_str(interpolation));

	gettimeofday(&t_start, NULL);
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	retvalue = cvResizeGaussian(image, toX, toY, interpolation, clamp);
	if (image->keywords.pixel_size_x > 0) image->keywords.pixel_size_x *= factor_X;
	if (image->keywords.pixel_size_y > 0) image->keywords.pixel_size_y *= factor_Y;
	free_wcs(image);
	reset_wcsdata(image);
	refresh_annotations(TRUE);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return retvalue;
}

// computes H matrix for rotation and crop
static void GetMatrixReframe(fits *image, rectangle area, double angle, int cropped, int *target_rx, int *target_ry, Homography *H) {
	*target_rx = area.w;
	*target_ry = area.h;
	double orig_x = (double)area.x;
	double orig_y = (double)area.y;
	if (!cropped) {
		point center = {orig_x + (double)*target_rx * 0.5, orig_y + (double)*target_ry * 0.5};
		cvGetBoundingRectSize(image, center, angle, target_rx, target_ry);
		orig_x = (double)((int)image->rx - *target_rx) * 0.5;
		orig_y = (double)((int)image->ry - *target_ry) * 0.5;
	}
	cvGetMatrixReframe(orig_x, orig_y, *target_rx, *target_ry, angle, H);
}

// wraps cvRotateImage to update WCS data as well
int verbose_rotate_fast(fits *image, int angle) {
	if (angle % 90 != 0) return 1;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	siril_log_color_message(
			_("Rotation (%s interpolation, angle=%g): processing...\n"), "green",
			_("No"), (double)angle);

	int orig_ry = image->ry; // required to compute flips afterwards
	int target_rx, target_ry;
	Homography H = { 0 };
	rectangle area = {0, 0, image->rx, image->ry};
	GetMatrixReframe(image, area, (double)angle, 0, &target_rx, &target_ry, &H);

	if (cvRotateImage(image, angle)) return 1;

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	if (has_wcs(image)) {
		cvApplyFlips(&H, orig_ry, target_ry);
		reframe_astrometry_data(image, &H);
		update_wcsdata_from_wcs(image);
		update_fits_header(image);
		refresh_annotations(FALSE);
	}
	return 0;
}

int verbose_rotate_image(fits *image, rectangle area, double angle, int interpolation,
		int cropped, gboolean clamp) {
	const char *str_inter;
	struct timeval t_start, t_end;

	switch (interpolation) {
		case -1:
			str_inter = _("No");
			break;
		case OPENCV_NEAREST:
			str_inter = _("Nearest-Neighbor");
			break;
		default:
		case OPENCV_LINEAR:
			str_inter = _("Bilinear");
			break;
		case OPENCV_AREA:
			str_inter = _("Pixel Area Relation");
			break;
		case OPENCV_CUBIC:
			str_inter = _("Bicubic");
			break;
		case OPENCV_LANCZOS4:
			str_inter = _("Lanczos4");
			break;
	}

	siril_log_color_message(
			_("Rotation (%s interpolation, angle=%g): processing...\n"), "green",
			str_inter, angle);
	gettimeofday(&t_start, NULL);
	on_clear_roi(); // ROI is cleared on geometry-altering operations

	int orig_ry = image->ry; // required to compute flips afterwards
	int target_rx, target_ry;
	Homography H = { 0 }, Hocv = { 0 };
	GetMatrixReframe(image, area, angle, cropped, &target_rx, &target_ry, &H);
	// The matrix has been written in display convention (i.e, siril flipped)
	// We need to convert to opencv convention (pixel-based) before applying to the image
	// The original matrix will still be used to reframe astrometry data
	Hocv = H;
	cvdisplay2ocv(&Hocv);
	if (cvTransformImage(image, target_rx, target_ry, Hocv, 1.f, interpolation, clamp, NULL)) return 1;

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	if (has_wcs(image)) {
		cvApplyFlips(&H, orig_ry, target_ry);
		reframe_astrometry_data(image, &H);
		update_wcsdata_from_wcs(image);
		update_fits_header(image);
		refresh_annotations(FALSE);
	}
	return 0;
}

static void mirrorx_ushort(fits *fit, gboolean verbose) {
	int line, axis;
	WORD *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	size_t line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);
	if (!swapline) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

static void mirrorx_float(fits *fit, gboolean verbose) {
	int line, axis;
	float *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	size_t line_size = fit->rx * sizeof(float);
	swapline = malloc(line_size);
	if (!swapline) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->fpdata[axis] + line * fit->rx;
			dst = fit->fpdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

void mirrorx(fits *fit, gboolean verbose) {
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	if (fit->type == DATA_USHORT) {
		mirrorx_ushort(fit, verbose);
	} else if (fit->type == DATA_FLOAT) {
		mirrorx_float(fit, verbose);
	}
	if (!strcmp(fit->keywords.row_order, "BOTTOM-UP"))
		sprintf(fit->keywords.row_order, "TOP-DOWN");
	else { //if (!strcmp(fit->row_order, "TOP-DOWN"))
		// let's create the keyword in all cases
		sprintf(fit->keywords.row_order, "BOTTOM-UP");
	}
	fit->history = g_slist_append(fit->history, strdup("TOP-DOWN mirror"));
	if (has_wcs(fit)) {
		Homography H = { 0 };
		cvGetEye(&H);
		H.h11 = -1.;
		H.h12 = (double)fit->ry;
		reframe_astrometry_data(fit, &H);
		update_wcsdata_from_wcs(fit);
		update_fits_header(fit);
		refresh_annotations(FALSE);
	}
}

void mirrory(fits *fit, gboolean verbose) {
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Vertical mirror: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	fits_flip_top_to_bottom(fit);
	fits_rotate_pi(fit);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	fit->history = g_slist_append(fit->history, strdup("Left-right mirror"));
	if (has_wcs(fit)) {
		Homography H = { 0 };
		cvGetEye(&H);
		H.h00 = -1.;
		H.h02 = (double)fit->rx;
		reframe_astrometry_data(fit, &H);
		update_wcsdata_from_wcs(fit);
		update_fits_header(fit);
		refresh_annotations(FALSE);
	}
}

/* inplace cropping of the image in fit
 * fit->data is not realloc, only fit->pdata points to a different area and
 * data is correctly written to this new area, which makes this function
 * quite dangerous to use when fit is used for something else afterwards.
 */
static int crop_ushort(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	rectangle bounds_cpy = { 0 };

	bounds_cpy.x = bounds->x;
	bounds_cpy.y = bounds->y;
	bounds_cpy.w = (bounds->x + bounds->w > fit->rx) ? fit->rx - bounds->x : bounds->w;
	bounds_cpy.h = (bounds->y + bounds->h > fit->ry) ? fit->ry - bounds->y : bounds->h;

	newnbdata = bounds_cpy.w * bounds_cpy.h;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *from = fit->pdata[layer]
			+ (fit->ry - bounds_cpy.y - bounds_cpy.h) * fit->rx + bounds_cpy.x;
		fit->pdata[layer] = fit->data + layer * newnbdata;
		WORD *to = fit->pdata[layer];
		int stridefrom = fit->rx - bounds_cpy.w;

		for (i = 0; i < bounds_cpy.h; ++i) {
			for (j = 0; j < bounds_cpy.w; ++j) {
				*to++ = *from++;
			}
			from += stridefrom;
		}
	}
	fit->rx = fit->naxes[0] = bounds_cpy.w;
	fit->ry = fit->naxes[1] = bounds_cpy.h;
	return 0;
}

static int crop_float(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	rectangle bounds_cpy = { 0 };

	bounds_cpy.x = bounds->x;
	bounds_cpy.y = bounds->y;
	bounds_cpy.w = (bounds->x + bounds->w > fit->rx) ? fit->rx - bounds->x : bounds->w;
	bounds_cpy.h = (bounds->y + bounds->h > fit->ry) ? fit->ry - bounds->y : bounds->h;

	newnbdata = bounds_cpy.w * bounds_cpy.h;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		float *from = fit->fpdata[layer]
			+ (fit->ry - bounds_cpy.y - bounds_cpy.h) * fit->rx + bounds_cpy.x;
		fit->fpdata[layer] = fit->fdata + layer * newnbdata;
		float *to = fit->fpdata[layer];
		int stridefrom = fit->rx - bounds_cpy.w;

		for (i = 0; i < bounds_cpy.h; ++i) {
			for (j = 0; j < bounds_cpy.w; ++j) {
				*to++ = *from++;
			}
			from += stridefrom;
		}
	}
	fit->rx = fit->naxes[0] = bounds_cpy.w;
	fit->ry = fit->naxes[1] = bounds_cpy.h;
	return 0;
}

int crop(fits *fit, rectangle *bounds) {
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	if (bounds->w <= 0 || bounds->h <= 0 || bounds->x < 0 || bounds->y < 0) return -1;
	if (bounds->x + bounds->w > fit->rx) return -1;
	if (bounds->y + bounds->h > fit->ry) return -1;
	int cfa = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern); // we don't need the validated value here because we just want to know if it's CFA, XTRANS or NONE
	switch (cfa) {
		case BAYER_FILTER_NONE:
			break;
		case BAYER_FILTER_RGGB: // Fallthrough intentional
		case BAYER_FILTER_BGGR:
		case BAYER_FILTER_GBRG:
		case BAYER_FILTER_GRBG:
			bounds->x &= ~1;
			bounds->y &= ~1;
			bounds->w &= ~1;
			bounds->h &= ~1;
			if(bounds->w < 2)
				bounds->w = 2;
			if (bounds->h < 2)
				bounds->h = 2;
			siril_log_message(_("Rounding selection to match CFA pattern\n"));
			break;
		default: // X-trans
			bounds->x = round_down_to_multiple(bounds->x, 6);
			bounds->y = round_down_to_multiple(bounds->y, 6);
			bounds->w = round_down_to_multiple(bounds->w, 6);
			bounds->h = round_down_to_multiple(bounds->h, 6);
			if(bounds->w < 2)
				bounds->w = 2;
			if (bounds->h < 2)
				bounds->h = 2;
			siril_log_message(_("Rounding selection to match CFA pattern\n"));
	}
	int orig_ry = fit->ry; // required to compute flips afterwards
	int target_rx, target_ry;
	Homography H = { 0 };
	gboolean wcs = has_wcs(fit);
	if (wcs)
		GetMatrixReframe(fit, *bounds, 0., 1, &target_rx, &target_ry, &H);

	if (fit->type == DATA_USHORT) {
		crop_ushort(fit, bounds);
	} else if (fit->type == DATA_FLOAT) {
		crop_float(fit, bounds);
	} else {
		return -1;
	}

	invalidate_stats_from_fit(fit);
	if (wcs) {
		cvApplyFlips(&H, orig_ry, target_ry);
		reframe_astrometry_data(fit, &H);
		update_wcsdata_from_wcs(fit);
		update_fits_header(fit);
		refresh_annotations(FALSE);
	}
	return 0;
}

gint64 crop_compute_size_hook(struct generic_seq_args *args, int nb_frames) {
	struct crop_sequence_data *c_args = (struct crop_sequence_data*) args->user;
	double ratio = (double)(c_args->area.h * c_args->area.w) / (double)(args->seq->rx * args->seq->ry);
	double fullseqsize = seq_compute_size(args->seq, nb_frames, args->output_type);
	return (gint64)(fullseqsize * ratio);
}

int crop_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct crop_sequence_data *c_args = (struct crop_sequence_data*) args->user;

	int ret = crop(fit, &(c_args->area));

	if (!ret) {
		char log[90];
		sprintf(log, _("Crop (x=%d, y=%d, w=%d, h=%d)"),
				c_args->area.x, c_args->area.y, c_args->area.w, c_args->area.h);
		fit->history = g_slist_append(fit->history, strdup(log));
	}
	return ret;
}

int crop_finalize_hook(struct generic_seq_args *args) {
	struct crop_sequence_data *data = (struct crop_sequence_data *) args->user;
	int retval = seq_finalize_hook(args);
	free(data);
	return retval;
}

int scale_compute_mem_limits_hook(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;

	/* The transformation memory consumption is:
		* the original image
		* the transformed image, including upscale if required (4x)
		*/
	unsigned int required = MB_per_orig_image + MB_per_scaled_image;
	// If interpolation clamping is set, 2x additional Mats of the same format
	// as the original image are required
	struct scale_sequence_data *s_args = args->user;
	if (s_args->clamp && (s_args->interpolation == OPENCV_CUBIC ||
			s_args->interpolation == OPENCV_LANCZOS4)) {
		float factor = (is_float) ? 0.25 : 0.5;
		required += (1 + factor) * MB_per_scaled_image;
	}

	if (limit > 0) {

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
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
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

gint64 scale_compute_size_hook(struct generic_seq_args *args, int nb_frames) {
	struct scale_sequence_data *s_args = (struct scale_sequence_data*) args->user;
	double ratio = s_args->scale * s_args->scale;
	double fullseqsize = seq_compute_size(args->seq, nb_frames, args->output_type);
	return (gint64)(fullseqsize * ratio);
}

int scale_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct scale_sequence_data *s_args = (struct scale_sequence_data*) args->user;
	int toX = fit->rx * s_args->scale;
	int toY = fit->ry * s_args->scale;
	s_args->retvalue = verbose_resize_gaussian(fit, toX, toY, s_args->interpolation, s_args->clamp);
	return s_args->retvalue;
}

int scale_finalize_hook(struct generic_seq_args *args) {
	struct scale_sequence_data *data = (struct scale_sequence_data *) args->user;
	int retval = seq_finalize_hook(args);
	free(data);
	return retval;
}

/* TODO: should we use the partial image? */
gpointer crop_sequence(struct crop_sequence_data *crop_sequence_data) {
	struct generic_seq_args *args = create_default_seqargs(crop_sequence_data->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = crop_sequence_data->seq->selnum;
	args->compute_size_hook = crop_compute_size_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = crop_finalize_hook;
	args->image_hook = crop_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Crop Sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(crop_sequence_data->prefix);
	args->load_new_sequence = TRUE;
	args->user = crop_sequence_data;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(crop_sequence_data->prefix);
		free(crop_sequence_data);
		free_generic_seq_args(args, TRUE);
		return GINT_TO_POINTER(1);
	}

	return GINT_TO_POINTER(0);
}

gpointer scale_sequence(struct scale_sequence_data *scale_sequence_data) {
	struct generic_seq_args *args = create_default_seqargs(scale_sequence_data->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = scale_sequence_data->seq->selnum;
	args->compute_mem_limits_hook = scale_compute_mem_limits_hook;
	args->compute_size_hook = scale_compute_size_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = scale_finalize_hook;
	args->image_hook = scale_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Scale Sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(scale_sequence_data->prefix);
	args->load_new_sequence = TRUE;
	args->user = scale_sequence_data;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(scale_sequence_data->prefix);
		free(scale_sequence_data);
		free_generic_seq_args(args, TRUE);
		return GINT_TO_POINTER(1);
	}

	return GINT_TO_POINTER(0);
}
