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

#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "core/arithm.h"
#include "core/masks.h"
#include "algos/astrometry_solver.h"
#include "algos/demosaicing.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "core/processing.h"
#include "opencv/opencv.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"

#include "geometry.h"

/* Mask helper functions */
// Helper function to resize mask
static int resize_mask(mask_t *mask, int old_rx, int old_ry, int new_rx, int new_ry, opencv_interpolation interpolation) {
	if (!mask || !mask->data) {
		return 0; // No mask to resize, not an error
	}

	// Determine element size based on bitpix
	size_t elem_size;
	switch (mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			return -1; // Invalid bitpix value
	}

	int newnbdata = new_rx * new_ry;
	void *newdata = malloc((size_t)newnbdata * elem_size);
	if (!newdata) {
		return -1;
	}

	// Create temporary fits structures for mask data
	fits temp_in = { 0 };
	temp_in.rx = old_rx;
	temp_in.ry = old_ry;
	temp_in.naxes[0] = old_rx;
	temp_in.naxes[1] = old_ry;
	temp_in.naxes[2] = 1;

	fits temp_out = { 0 };
	temp_out.rx = new_rx;
	temp_out.ry = new_ry;
	temp_out.naxes[0] = new_rx;
	temp_out.naxes[1] = new_ry;
	temp_out.naxes[2] = 1;

	if (mask->bitpix == 32) {
		temp_in.type = DATA_FLOAT;
		temp_in.fdata = (float *)mask->data;
		temp_in.fpdata[0] = temp_in.fdata;
		temp_in.fpdata[1] = temp_in.fdata;
		temp_in.fpdata[2] = temp_in.fdata;

		temp_out.type = DATA_FLOAT;
		temp_out.fdata = (float *)newdata;
		temp_out.fpdata[0] = temp_out.fdata;
		temp_out.fpdata[1] = temp_out.fdata;
		temp_out.fpdata[2] = temp_out.fdata;
	} else {
		temp_in.type = DATA_USHORT;
		temp_in.data = (WORD *)mask->data;
		temp_in.pdata[0] = temp_in.data;
		temp_in.pdata[1] = temp_in.data;
		temp_in.pdata[2] = temp_in.data;

		temp_out.type = DATA_USHORT;
		temp_out.data = (WORD *)newdata;
		temp_out.pdata[0] = temp_out.data;
		temp_out.pdata[1] = temp_out.data;
		temp_out.pdata[2] = temp_out.data;
	}

	if (cvResizeGaussian(&temp_in, new_rx, new_ry, interpolation, FALSE)) {
		free(newdata);
		return -1;
	}

	// Replace old mask data
	void *tmp = mask->data;
	mask->data = newdata;
	free(tmp);

	return 0;
}

// Helper function to bin mask
static int bin_mask(mask_t *mask, int old_rx, int old_ry, int bin_factor, gboolean mean) {
	if (!mask || !mask->data || bin_factor <= 0) {
		return 0; // No mask to bin, not an error
	}

	int new_rx = old_rx / bin_factor;
	int new_ry = old_ry / bin_factor;

	// Determine element size based on bitpix
	size_t elem_size;
	switch (mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			return -1;
	}

	int newnbdata = new_rx * new_ry;
	void *newdata = malloc((size_t)newnbdata * elem_size);
	if (!newdata) {
		return -1;
	}

	if (mask->bitpix == 32) {
		float *buf = (float *)mask->data;
		float *new_buf = (float *)newdata;

		long k = 0;
		for (int row = 0; row < old_ry - bin_factor + 1; row += bin_factor) {
			for (int col = 0; col < old_rx - bin_factor + 1; col += bin_factor) {
				int c = 0;
				new_buf[k] = 0;
				for (int i = 0; i < bin_factor; i++) {
					for (int j = 0; j < bin_factor; j++) {
						new_buf[k] += buf[i + col + (j + row) * old_rx];
						c++;
					}
				}
				if (mean) new_buf[k] /= c;
				k++;
			}
		}
	} else {
		WORD *buf = (WORD *)mask->data;
		WORD *new_buf = (WORD *)newdata;

		long k = 0;
		for (int row = 0; row < old_ry - bin_factor + 1; row += bin_factor) {
			for (int col = 0; col < old_rx - bin_factor + 1; col += bin_factor) {
				int c = 0;
				int tmp = 0;
				for (int i = 0; i < bin_factor; i++) {
					for (int j = 0; j < bin_factor; j++) {
						tmp += buf[i + col + (j + row) * old_rx];
						c++;
					}
				}
				if (mean) tmp /= c;
				new_buf[k] = (mask->bitpix == 8) ? (uint8_t)tmp : truncate_to_WORD(tmp);
				k++;
			}
		}
	}

	void *tmp = mask->data;
	mask->data = newdata;
	free(tmp);

	return 0;
}

// Helper function to rotate mask 180 degrees
static int rotate_mask_pi(mask_t *mask, int rx, int ry) {
	if (!mask || !mask->data) {
		return 0;
	}

	size_t elem_size;
	switch (mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			return -1;
	}

	size_t line_size = rx * elem_size;
	void *line1 = malloc(line_size);
	void *line2 = malloc(line_size);
	if (!line1 || !line2) {
		if (line1) free(line1);
		if (line2) free(line2);
		return -1;
	}

	uint8_t *data = (uint8_t *)mask->data;

	for (int line = 0; line < ry / 2; line++) {
		uint8_t *src = data + line * rx * elem_size;
		uint8_t *dst = data + (ry - line - 1) * rx * elem_size;

		memcpy(line1, src, line_size);
		// Reverse line1
		for (int i = 0; i < rx / 2; i++) {
			memcpy(line2, line1 + i * elem_size, elem_size);
			memcpy(line1 + i * elem_size, line1 + (rx - i - 1) * elem_size, elem_size);
			memcpy(line1 + (rx - i - 1) * elem_size, line2, elem_size);
		}

		memcpy(line2, dst, line_size);
		// Reverse line2
		for (int i = 0; i < rx / 2; i++) {
			uint8_t temp[32]; // Large enough for any elem_size we support
			memcpy(temp, line2 + i * elem_size, elem_size);
			memcpy(line2 + i * elem_size, line2 + (rx - i - 1) * elem_size, elem_size);
			memcpy(line2 + (rx - i - 1) * elem_size, temp, elem_size);
		}

		memcpy(src, line2, line_size);
		memcpy(dst, line1, line_size);
	}

	if (ry & 1) {
		// Reverse middle line
		uint8_t *src = data + (ry / 2) * rx * elem_size;
		for (int i = 0; i < rx / 2; i++) {
			memcpy(line1, src + i * elem_size, elem_size);
			memcpy(src + i * elem_size, src + (rx - i - 1) * elem_size, elem_size);
			memcpy(src + (rx - i - 1) * elem_size, line1, elem_size);
		}
	}

	free(line1);
	free(line2);
	return 0;
}

// Helper function to mirror mask horizontally
static int mirrorx_mask(mask_t *mask, int rx, int ry) {
	if (!mask || !mask->data) {
		return 0;
	}

	size_t elem_size;
	switch (mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			return -1;
	}

	size_t line_size = rx * elem_size;
	void *swapline = malloc(line_size);
	if (!swapline) {
		return -1;
	}

	uint8_t *data = (uint8_t *)mask->data;

	for (int line = 0; line < ry / 2; line++) {
		uint8_t *src = data + line * rx * elem_size;
		uint8_t *dst = data + (ry - line - 1) * rx * elem_size;

		memcpy(swapline, src, line_size);
		memcpy(src, dst, line_size);
		memcpy(dst, swapline, line_size);
	}

	free(swapline);
	return 0;
}

// Helper function to transform mask with homography
static int transform_mask(mask_t *mask, int old_rx, int old_ry, int new_rx, int new_ry,
                         Homography H, opencv_interpolation interpolation) {
	if (!mask || !mask->data) {
		return 0;
	}

	int old_nbdata = old_rx * old_ry;
	int new_nbdata = new_rx * new_ry;

	// Convert mask data to USHORT for processing if it's 8-bit
	WORD *input_ushort = NULL;
	float *input_float = NULL;
	gboolean needs_conversion = FALSE;

	if (mask->bitpix == 32) {
		// Float mask - copy to float array
		input_float = malloc(old_nbdata * sizeof(float));
		if (!input_float) return -1;
		memcpy(input_float, mask->data, old_nbdata * sizeof(float));
	} else {
		// 8-bit or 16-bit mask - convert to USHORT
		input_ushort = malloc(old_nbdata * sizeof(WORD));
		if (!input_ushort) return -1;

		if (mask->bitpix == 8) {
			needs_conversion = TRUE;
			uint8_t *src = (uint8_t *)mask->data;
			for (int i = 0; i < old_nbdata; i++) {
				input_ushort[i] = (WORD)src[i];
			}
		} else {
			memcpy(input_ushort, mask->data, old_nbdata * sizeof(WORD));
		}
	}

	// Create temporary fits structure
	fits temp_in = { 0 };
	temp_in.rx = old_rx;
	temp_in.ry = old_ry;
	temp_in.naxes[0] = old_rx;
	temp_in.naxes[1] = old_ry;
	temp_in.naxes[2] = 1;

	if (mask->bitpix == 32) {
		temp_in.type = DATA_FLOAT;
		temp_in.fdata = input_float;
		temp_in.fpdata[0] = temp_in.fdata;
		temp_in.fpdata[1] = temp_in.fdata;
		temp_in.fpdata[2] = temp_in.fdata;
	} else {
		temp_in.type = DATA_USHORT;
		temp_in.data = input_ushort;
		temp_in.pdata[0] = temp_in.data;
		temp_in.pdata[1] = temp_in.data;
		temp_in.pdata[2] = temp_in.data;
	}

	// Convert H to OpenCV convention
	Homography Hocv = H;
	cvdisplay2ocv(&Hocv);

	// cvTransformImage will free the old data and allocate new data
	if (cvTransformImage(&temp_in, new_rx, new_ry, Hocv, 1.f, interpolation, FALSE, NULL)) {
		if (input_ushort) free(input_ushort);
		if (input_float) free(input_float);
		return -1;
	}

	// cvTransformImage has replaced the data pointer with transformed data
	// Now convert back to original format if needed
	free(mask->data);

	if (mask->bitpix == 32) {
		mask->data = temp_in.fdata;
	} else if (mask->bitpix == 8 && needs_conversion) {
		// Convert USHORT back to uint8_t
		uint8_t *new_data = malloc(new_nbdata * sizeof(uint8_t));
		if (!new_data) {
			free(temp_in.data);
			return -1;
		}
		for (int i = 0; i < new_nbdata; i++) {
			new_data[i] = (uint8_t)(temp_in.data[i] > 255 ? 255 : temp_in.data[i]);
		}
		free(temp_in.data);
		mask->data = new_data;
	} else {
		// bitpix == 16
		mask->data = temp_in.data;
	}

	return 0;
}


static int crop_mask(mask_t *mask, rectangle *bounds, int rx, int ry) {
	if (!mask || !mask->data || !bounds) {
		return 0; // No mask to crop, not an error
	}

	// Determine element size based on bitpix
	size_t elem_size;
	switch (mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			return -1; // Invalid bitpix value
	}

	// Clamp bounds to image dimensions
	rectangle bounds_cpy = { 0 };
	bounds_cpy.x = bounds->x;
	bounds_cpy.y = bounds->y;
	bounds_cpy.w = (bounds->x + bounds->w > rx) ? rx - bounds->x : bounds->w;
	bounds_cpy.h = (bounds->y + bounds->h > ry) ? ry - bounds->y : bounds->h;

	// Validate resulting crop region
	if (bounds_cpy.w <= 0 || bounds_cpy.h <= 0 ||
	    bounds_cpy.x < 0 || bounds_cpy.y < 0 ||
	    bounds_cpy.x >= rx || bounds_cpy.y >= ry) {
		return -1;
	}

	int newnbdata = bounds_cpy.w * bounds_cpy.h;

	// Allocate new buffer for cropped mask
	void *newdata = malloc((size_t)newnbdata * elem_size);
	if (!newdata) {
		return -1;
	}

	// Source pointer: bottom-left origin (FITS coordinate system)
	uint8_t *from = (uint8_t *)mask->data +
	                (ry - bounds_cpy.y - bounds_cpy.h) * rx * elem_size +
	                bounds_cpy.x * elem_size;

	// Destination pointer
	uint8_t *to = (uint8_t *)newdata;

	size_t row_bytes = bounds_cpy.w * elem_size;
	size_t stride_bytes = rx * elem_size;

	// Copy row by row
	for (int i = 0; i < bounds_cpy.h; ++i) {
		memcpy(to, from, row_bytes);
		to += row_bytes;
		from += stride_bytes;
	}

	// Replace old mask data
	void *tmp = mask->data;
	mask->data = newdata;
	free(tmp);

	return 0;
}

/* Image geometry functions */

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
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	gboolean tmp_mask_active = fit->mask_active;
	set_mask_active(fit, FALSE);

	int old_rx = fit->rx;
	int old_ry = fit->ry;

	if (fit->type == DATA_USHORT) {
		fits_binning_ushort(fit, factor, mean);
	} else if (fit->type == DATA_FLOAT) {
		fits_binning_float(fit, factor, mean);
	}

	if (fit->mask) {
		if (bin_mask(fit->mask, old_rx, old_ry, factor, mean)) {
			siril_log_color_message(_("Error binning mask\n"), "red");
			free_mask(fit->mask);
			fit->mask = NULL;
			show_or_hide_mask_tab();
			return -1;
		} else {
			set_mask_active(fit, tmp_mask_active);
		}
	}

	free_wcs(fit);
	reset_wcsdata(fit);
	refresh_annotations(TRUE);

	return 0;
}

const char* interp_to_str(int interpolation) {
	const char *str_inter;
	switch (interpolation) {
		case -1:
		case OPENCV_NONE:
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
	return str_inter;
}

/* These functions do not more than resize_gaussian and rotate_image
 * except for console outputs.
 * Indeed, siril_log_message seems not working in a cpp file */
int verbose_resize_gaussian(fits *image, int toX, int toY, opencv_interpolation interpolation, gboolean clamp) {
	int retvalue;
	float factor_X = (float)image->rx / (float)toX;
	float factor_Y = (float)image->ry / (float)toY;
	gboolean tmp_mask_active = FALSE;
	if (image->mask) {
		tmp_mask_active = image->mask_active;
		set_mask_active(image, FALSE);
	}

	int old_rx = image->rx;
	int old_ry = image->ry;

	on_clear_roi(); // ROI is cleared on geometry-altering operations
	retvalue = cvResizeGaussian(image, toX, toY, interpolation, clamp);

	if (retvalue == 0 && image->mask) {
		if (resize_mask(image->mask, old_rx, old_ry, toX, toY, interpolation)) {
			siril_log_color_message(_("Error resizing mask\n"), "red");
			free_mask(image->mask);
			image->mask = NULL;
			show_or_hide_mask_tab();
		} else {
			set_mask_active(image, tmp_mask_active);
		}
	}

	if (image->keywords.pixel_size_x > 0) image->keywords.pixel_size_x *= factor_X;
	if (image->keywords.pixel_size_y > 0) image->keywords.pixel_size_y *= factor_Y;
	free_wcs(image);
	reset_wcsdata(image);
	refresh_annotations(TRUE);

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
	if (angle % 90 != 0) return 1; // only for multiples of 90 \deg
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	gboolean tmp_mask_active = FALSE;
	if (image->mask) {
		tmp_mask_active = image->mask_active;
		set_mask_active(image, FALSE);
	}

	int orig_ry = image->ry;
	int orig_rx = image->rx;
	int target_rx, target_ry;
	Homography H = { 0 };
	rectangle area = {0, 0, image->rx, image->ry};
	GetMatrixReframe(image, area, (double)angle, 0, &target_rx, &target_ry, &H);

	if (cvRotateImage(image, angle)) return 1;

	if (image->mask) {
		// OPENCV_NEAREST is fine because we are rotating by a multiple of 90 \deg
		if (transform_mask(image->mask, orig_rx, orig_ry, target_rx, target_ry, H, OPENCV_NEAREST)) {
			siril_log_color_message(_("Error rotating mask\n"), "red");
			free_mask(image->mask);
			image->mask = NULL;
			set_mask_active(image, FALSE);
			show_or_hide_mask_tab();
		} else {
			set_mask_active(image, tmp_mask_active);
		}
	}

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
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	gboolean tmp_mask_active = FALSE;
	if (image->mask) {
		tmp_mask_active = image->mask_active;
		set_mask_active(image, FALSE);
	}

	int orig_ry = image->ry;
	int orig_rx = image->rx;
	int target_rx, target_ry;
	Homography H = { 0 }, Hocv = { 0 };
	GetMatrixReframe(image, area, angle, cropped, &target_rx, &target_ry, &H);

	Hocv = H;
	cvdisplay2ocv(&Hocv);
	if (cvTransformImage(image, target_rx, target_ry, Hocv, 1.f, interpolation, clamp, NULL)) return 1;

	if (image->mask) {
		if (transform_mask(image->mask, orig_rx, orig_ry, target_rx, target_ry, H, OPENCV_CUBIC)) {
			siril_log_color_message(_("Error rotating mask\n"), "red");
			free_mask(image->mask);
			image->mask = NULL;
			show_or_hide_mask_tab();
			set_mask_active(image, FALSE);
		} else {
			set_mask_active(image, tmp_mask_active);
		}
	}

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
	gboolean tmp_mask_active = fit->mask_active;
	set_mask_active(fit, FALSE);

	if (fit->type == DATA_USHORT) {
		mirrorx_ushort(fit, verbose);
	} else if (fit->type == DATA_FLOAT) {
		mirrorx_float(fit, verbose);
	}

	if (fit->mask) {
		if (mirrorx_mask(fit->mask, fit->rx, fit->ry)) {
			siril_log_color_message(_("Error mirroring mask\n"), "red");
			free_mask(fit->mask);
			fit->mask = NULL;
			show_or_hide_mask_tab();
		} else {
			set_mask_active(fit, tmp_mask_active);
		}
	}

	if (!strcmp(fit->keywords.row_order, "BOTTOM-UP"))
		sprintf(fit->keywords.row_order, "TOP-DOWN");
	else {
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
	gboolean tmp_mask_active = fit->mask_active;
	set_mask_active(fit, FALSE);

	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Vertical mirror: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	fits_flip_top_to_bottom(fit);
	fits_rotate_pi(fit);

	if (fit->mask) {
		// For vertical mirror: flip top-to-bottom then rotate 180
		if (mirrorx_mask(fit->mask, fit->rx, fit->ry) ||
		    rotate_mask_pi(fit->mask, fit->rx, fit->ry)) {
			siril_log_color_message(_("Error mirroring mask\n"), "red");
			free_mask(fit->mask);
			fit->mask = NULL;
			show_or_hide_mask_tab();
		} else {
			set_mask_active(fit, tmp_mask_active);
		}
	}

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

static int crop_ushort(fits *fit, rectangle *bounds) {
	if (!fit || !bounds) {
		return -1;
	}

	// Validate and clamp bounds to image dimensions
	rectangle bounds_cpy = { 0 };
	bounds_cpy.x = bounds->x;
	bounds_cpy.y = bounds->y;
	bounds_cpy.w = (bounds->x + bounds->w > fit->rx) ? fit->rx - bounds->x : bounds->w;
	bounds_cpy.h = (bounds->y + bounds->h > fit->ry) ? fit->ry - bounds->y : bounds->h;

	// Validate the resulting crop region
	if (bounds_cpy.w <= 0 || bounds_cpy.h <= 0 ||
	    bounds_cpy.x < 0 || bounds_cpy.y < 0 ||
	    bounds_cpy.x >= fit->rx || bounds_cpy.y >= fit->ry) {
		return -1;
	}

	int newnbdata = bounds_cpy.w * bounds_cpy.h;
	int nlayers = fit->naxes[2];

	// Allocate new buffer for cropped data
	WORD *newdata = malloc(newnbdata * nlayers * sizeof(WORD));
	if (!newdata) {
		return -1;
	}

	// Copy cropped data from each layer
	for (int layer = 0; layer < nlayers; ++layer) {
		// Source pointer: start at bottom-left of crop region (FITS coordinate system)
		WORD *from = fit->pdata[layer] +
		             (fit->ry - bounds_cpy.y - bounds_cpy.h) * fit->rx + bounds_cpy.x;

		// Destination pointer in new buffer
		WORD *to = newdata + layer * newnbdata;

		// Copy row by row
		for (int i = 0; i < bounds_cpy.h; ++i) {
			memcpy(to, from, bounds_cpy.w * sizeof(WORD));
			to += bounds_cpy.w;
			from += fit->rx;
		}
	}

	// Free old data and update pointers
	free(fit->data);
	fit->data = newdata;

	// Update layer pointers
	for (int layer = 0; layer < nlayers; ++layer) {
		fit->pdata[layer] = fit->data + layer * newnbdata;
	}

	// For single-layer images, set all pointers to the same data
	if (nlayers == 1) {
		fit->pdata[1] = fit->data;
		fit->pdata[2] = fit->data;
	}

	// Update dimensions
	fit->rx = fit->naxes[0] = bounds_cpy.w;
	fit->ry = fit->naxes[1] = bounds_cpy.h;

	return 0;
}

static int crop_float(fits *fit, rectangle *bounds) {
	if (!fit || !bounds) {
		return -1;
	}

	// Validate and clamp bounds to image dimensions
	rectangle bounds_cpy = { 0 };
	bounds_cpy.x = bounds->x;
	bounds_cpy.y = bounds->y;
	bounds_cpy.w = (bounds->x + bounds->w > fit->rx) ? fit->rx - bounds->x : bounds->w;
	bounds_cpy.h = (bounds->y + bounds->h > fit->ry) ? fit->ry - bounds->y : bounds->h;

	// Validate the resulting crop region
	if (bounds_cpy.w <= 0 || bounds_cpy.h <= 0 ||
	    bounds_cpy.x < 0 || bounds_cpy.y < 0 ||
	    bounds_cpy.x >= fit->rx || bounds_cpy.y >= fit->ry) {
		return -1;
	}

	int newnbdata = bounds_cpy.w * bounds_cpy.h;
	int nlayers = fit->naxes[2];

	// Allocate new buffer for cropped data
	float *newfdata = malloc(newnbdata * nlayers * sizeof(float));
	if (!newfdata) {
		return -1;
	}

	// Copy cropped data from each layer
	for (int layer = 0; layer < nlayers; ++layer) {
		// Source pointer: start at bottom-left of crop region (FITS coordinate system)
		float *from = fit->fpdata[layer] +
		              (fit->ry - bounds_cpy.y - bounds_cpy.h) * fit->rx + bounds_cpy.x;

		// Destination pointer in new buffer
		float *to = newfdata + layer * newnbdata;

		// Copy row by row
		for (int i = 0; i < bounds_cpy.h; ++i) {
			memcpy(to, from, bounds_cpy.w * sizeof(float));
			to += bounds_cpy.w;
			from += fit->rx;
		}
	}

	// Free old data and update pointers
	free(fit->fdata);
	fit->fdata = newfdata;

	// Update layer pointers
	for (int layer = 0; layer < nlayers; ++layer) {
		fit->fpdata[layer] = fit->fdata + layer * newnbdata;
	}

	// For single-layer images, set all pointers to the same data
	if (nlayers == 1) {
		fit->fpdata[1] = fit->fdata;
		fit->fpdata[2] = fit->fdata;
	}

	// Update dimensions
	fit->rx = fit->naxes[0] = bounds_cpy.w;
	fit->ry = fit->naxes[1] = bounds_cpy.h;

	return 0;
}

int crop(fits *fit, rectangle *bounds) {
	on_clear_roi(); // ROI is cleared on geometry-altering operations
	gboolean tmp_mask_active = fit->mask_active;
	set_mask_active(fit, FALSE);
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
	int orig_rx = fit->rx; // required to compute flips afterwards
	int orig_ry = fit->ry; // required to compute flips afterwards
	int target_rx, target_ry;
	Homography H = { 0 };
	gboolean wcs = has_wcs(fit);
	if (wcs)
		GetMatrixReframe(fit, *bounds, 0., 1, &target_rx, &target_ry, &H);

	if (fit->type == DATA_USHORT) {
		if (crop_ushort(fit, bounds)) return -1;
	} else if (fit->type == DATA_FLOAT) {
		if (crop_float(fit, bounds)) return -1;
	} else {
		return -1;
	}

	if (fit->mask) {
		if (crop_mask(fit->mask, bounds, orig_rx, orig_ry)) {
			siril_log_color_message(_("Error cropping mask\n"), "red");
			free_mask(fit->mask);
			fit->mask = NULL;
			show_or_hide_mask_tab();
			return -1;
		} else {
			set_mask_active(fit, tmp_mask_active);
		}
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
	if (args->seq->type == SEQ_INTERNAL) {
		// For SEQ_INTERNAL we update the sequence in place, we must update the metadata too
		args->seq->imgparam[o].rx = c_args->area.w;
		args->seq->imgparam[o].ry = c_args->area.h;
	}
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
	args->already_in_a_thread = args->seq->type == SEQ_INTERNAL;
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

/*****************************************************************************
 *      G E O M E T R Y   A R G S   A L L O C A T O R S   A N D   D E S T R U C T O R S
 ****************************************************************************/

struct binning_args *new_binning_args() {
	struct binning_args *args = calloc(1, sizeof(struct binning_args));
	if (args) {
		args->destroy_fn = free_binning_args;
	}
	return args;
}

void free_binning_args(void *ptr) {
	if (!ptr) return;
	free(ptr);
}

struct crop_args *new_crop_args() {
	struct crop_args *args = calloc(1, sizeof(struct crop_args));
	if (args) {
		args->destroy_fn = free_crop_args;
	}
	return args;
}

void free_crop_args(void *ptr) {
	if (!ptr) return;
	free(ptr);
}

struct mirror_args *new_mirror_args() {
	struct mirror_args *args = calloc(1, sizeof(struct mirror_args));
	if (args) {
		args->destroy_fn = free_mirror_args;
	}
	return args;
}

void free_mirror_args(void *ptr) {
	if (!ptr) return;
	free(ptr);
}

struct resample_args *new_resample_args() {
	struct resample_args *args = calloc(1, sizeof(struct resample_args));
	if (args) {
		args->destroy_fn = free_resample_args;
	}
	return args;
}

void free_resample_args(void *ptr) {
	if (!ptr) return;
	free(ptr);
}

struct rotation_args *new_rotation_args() {
	struct rotation_args *args = calloc(1, sizeof(struct rotation_args));
	if (args) {
		args->destroy_fn = free_rotation_args;
	}
	return args;
}

void free_rotation_args(void *ptr) {
	if (!ptr) return;
	free(ptr);
}

/**********************************************************************
 *      G E O M E T R Y   L O G   H O O K S   F O R   W O R K E R
 **********************************************************************/

gchar *crop_log_hook(gpointer p, log_hook_detail detail) {
	struct crop_args *params = (struct crop_args*) p;
	gchar *msg = g_strdup_printf(_("Crop (x=%d, y=%d, w=%d, h=%d)"),
			params->area.x, params->area.y, params->area.w,
			params->area.h);
	return msg;
}

gchar *resample_log_hook(gpointer p, log_hook_detail detail) {
	struct resample_args *params = (struct resample_args *)p;

	gchar *msg = g_strdup_printf(_("Resampling to %d x %d pixels with %s interpolation"),
			params->toX, params->toY, interp_to_str(params->interpolation));
	return msg;
}

gchar *binning_log_hook(gpointer p, log_hook_detail detail) {
	struct binning_args *params = (struct binning_args *) p;
	gchar *msg = g_strdup_printf(_("Binning x%d (%s)"), params->factor, params->mean ? _("average") : _("sum"));
	return msg;
}

gchar *rotation_log_hook(gpointer p, log_hook_detail detail) {
	struct rotation_args *params = (struct rotation_args *) p;
	gchar *msg = NULL;
	if (detail == SUMMARY) {
		msg = g_strdup_printf(_("Rotation (%.1fdeg)"), params->angle);
	} else {
		msg = g_strdup_printf(_("Rotation (%.1fdeg, cropped=%s, clamped=%s, %s interpolation)"), params->angle,
				params->cropped ? _("TRUE") : _("FALSE"), params->clamp ? _("TRUE") : _("FALSE"),
				interp_to_str(params->interpolation));
	}
	return msg;
}

/*****************************************************************************
 *      G E O M E T R Y   I M A G E   H O O K S   F O R   W O R K E R
 ****************************************************************************/

int binning_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct binning_args *params = (struct binning_args *)args->user;
	if (!params)
		return 1;

	return fits_binning(fit, params->factor, params->mean);
}

static gboolean crop_gui_updates(gpointer user) {
	clear_stars_list(TRUE);
	delete_selected_area();
	reset_display_offset();
	update_zoom_label();
	return FALSE;
}

int crop_image_hook_single(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct crop_args *params = (struct crop_args *)args->user;
	if (!params)
		return 1;
	int retval = crop(fit, &params->area);
	gui_function(crop_gui_updates, NULL);
	return retval;
}

int mirrorx_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	mirrorx(fit, FALSE);
	return 0;
}

int mirrory_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	mirrory(fit, FALSE);
	return 0;
}

int resample_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct resample_args *params = (struct resample_args *)args->user;
	if (!params)
		return 1;

	int retval = verbose_resize_gaussian(fit, params->toX, params->toY,
	                                params->interpolation, params->clamp);
	gui_function(update_MenuItem, NULL);
	return retval;
}

int rotation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct rotation_args *params = (struct rotation_args *)args->user;
	if (!params)
		return 1;

	// Check for fast rotation (90 degree increments)
	int angle_int = (int)params->angle;
	int retval;
	if (params->angle == angle_int && angle_int % 90 == 0 &&
	    params->area.x == 0 && params->area.y == 0 &&
	    params->area.w == fit->rx && params->area.h == fit->ry) {
		retval = verbose_rotate_fast(fit, angle_int);
	} else {
		retval = verbose_rotate_image(fit, params->area, params->angle,
	                             params->interpolation, params->cropped,
	                             params->clamp);
	}

	// If a selection is set, we set it to the entire image
	if (com.selection.w > 0 && com.selection.h > 0) {
		com.selection = (rectangle){ 0, 0, gfit->rx, gfit->ry };
		gui_function(new_selection_zone, NULL);
	}
	update_zoom_label();
	return retval;
}
