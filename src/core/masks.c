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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "core/siril_log.h"
#include "filters/synthstar.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/histogram.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/user_polygons.h"
#include "io/image_format_fits.h"
#include "masks.h"
#include "opencv/opencv.h" // for mask functions that use OpenCV
// (feathering, Gaussian blur and set_poly_in_mask())

void free_mask(mask_t* mask) {
	free(mask->data);
	free(mask);
}

gpointer clear_mask_worker(gpointer p) {
	if (gfit && gfit->mask) {
		set_mask_active(gfit, FALSE);
		free_mask(gfit->mask);
		gfit->mask = NULL;
		show_or_hide_mask_tab();
		siril_log_message(_("Mask cleared\n"));
	}
	siril_add_idle(end_generic, NULL);
	return FALSE;
}

gpointer autostretch_mask_worker(gpointer p) {
	if (gfit && gfit->mask) {
		mask_autostretch(gfit);
		queue_redraw_mask();
	}
	siril_add_idle(end_generic, NULL);
	return FALSE;
}

gpointer binarize_mask_from_gui_worker(gpointer p) {
	if (gfit && gfit->mask) {
		gchar *msg = g_strdup_printf(_("This will binarize the mask based on the range slider values. Values between %f and %f will be masked TRUE, other values will be masked FALSE"), (float) gui.lo, (float) gui.hi);
		if (siril_confirm_dialog(_("Binarize mask"), msg, _("Proceed"))) {
			mask_binarize(gfit, (float) gui.lo, (float) gui.hi);
			queue_redraw_mask();
			siril_log_message(_("Mask binarized (TRUE between %f and %f)\n"), (float) gui.lo, (float) gui.hi);
		}
	}
	siril_add_idle(end_generic, NULL);
	return FALSE;
}

gpointer invert_mask_worker(gpointer p) {
	if (gfit && gfit->mask) {
		mask_invert(gfit);
		queue_redraw_mask();
		siril_log_message(_("Mask inverted\n"));
	}
	siril_add_idle(end_generic, NULL);
	return FALSE;
}

gboolean set_mask_active_idle(gpointer p) {
	gboolean state = GPOINTER_TO_INT(p);
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("mask_active_check"));
	g_signal_handlers_block_by_func(GTK_TOGGLE_BUTTON(button), on_mask_active_toggled, NULL);
	gtk_toggle_button_set_active(button, state);
	g_signal_handlers_unblock_by_func(GTK_TOGGLE_BUTTON(button), on_mask_active_toggled, NULL);
	return FALSE;
}

void set_mask_active(fits *fit, gboolean state) {
	fit->mask_active = state;
	if (fit == gfit && !com.headless) {
		gui_function(set_mask_active_idle, GINT_TO_POINTER(state));
	}
}

// Create a test mask : the left half of the image has mask value 255,
// the right half has mask value 0. Any existing mask is
// freed and remade. Return 0 on success.

FAST_MATH_PUSH
int mask_create_test(fits *fit, uint8_t bitpix) {
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32))
		return 1;

	if (fit->mask)
		free_mask(fit->mask);

	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fit->mask->bitpix = bitpix;
	fit->mask->data = malloc(fit->rx * fit->ry * bitpix);
	if (!fit->mask->data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	size_t rx = fit->rx, ry = fit->ry, hrx = rx / 2;

	switch (bitpix) {
		case 8: {
			uint8_t* restrict m = (uint8_t*) fit->mask->data;
			for (size_t y = 0 ; y < ry ; y++) {
				for (size_t x = 0 ; x < rx ; x++) {
					m[y * rx + x] = x < hrx ? 255 : 0;
				}
			}
			break;
		}
		case 16: {
			uint16_t* restrict m = (uint16_t*) fit->mask->data;
			for (size_t y = 0 ; y < ry ; y++) {
				for (size_t x = 0 ; x < rx ; x++) {
					m[y * rx + x] = x < hrx ? 65535 : 0;
				}
			}
			break;
		}
		case 32: {
			float* restrict m = (float*) fit->mask->data;
			for (size_t y = 0 ; y < ry ; y++) {
				for (size_t x = 0 ; x < rx ; x++) {
					m[y * rx + x] = x < hrx ? 1.f : 0;
				}
			}
			break;
		}
	}
	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP

// Create a mask with ones-like property (i.e. full applicability: we
// use 8-bit masks so all values are set to 255). Any existing mask is
// freed and remade. Return 0 on success.

FAST_MATH_PUSH
int mask_create_ones_like(fits *fit, uint8_t bitpix) {
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32))
		return 1;

	if (fit->mask)
		free_mask(fit->mask);

	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fit->mask->bitpix = bitpix;
	fit->mask->data = malloc(fit->rx * fit->ry * bitpix);
	if (!fit->mask->data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	size_t rx = fit->rx, ry = fit->ry;

	switch (bitpix) {
		case 8: {
			uint8_t* m = (uint8_t*) fit->mask->data;
			memset(m, 255, rx * ry);
			break;
		}
		case 16: {
			uint16_t* m = (uint16_t*) fit->mask->data;
			memset(m, 0xFF, rx * ry * sizeof(uint16_t));
			break;
		}
		case 32: {
			float* restrict m = (float*) fit->mask->data;
			size_t n = rx * ry;
			for (size_t i = 0; i < n; ++i)
				m[i] = 1.0f;
			break;
		}
		default:
			siril_debug_print("Error! Unhandled bitpix in mask_create_ones_like\n");
			free_mask(fit->mask);
			fit->mask = NULL;
			return 1;
	}
	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP

// Create a mask with zeroes-like property. Any existing mask is freed
// and remade. Returns 0 on success.

int mask_create_zeroes_like(fits *fit, uint8_t bitpix) {
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32))
		return 1;

	if (fit->mask)
		free_mask(fit->mask);

	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fit->mask->bitpix = bitpix;
	fit->mask->data = calloc(fit->rx * fit->ry, bitpix);
	if (!fit->mask->data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}

// TODO: should the following functions take account of orig_bitpix or something to identify 16-bit FITS
// where the data is actually 8-bit (ie has  a max value of 255)?

FAST_MATH_PUSH
int mask_create_from_channel(fits *fit, fits *source, int chan, uint8_t bitpix) {
	if (!fit) return 1; // no FITS struct
	if (!source) return 1; // no source FITS struct
	if (chan >= source->naxes[2]) return 1; // channel out of range

	set_mask_active(fit, FALSE);
	if (fit->mask)
		free_mask(fit->mask);

	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fit->mask->bitpix = bitpix;
	fit->mask->data = malloc(fit->rx * fit->ry * bitpix);
	if (!fit->mask->data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	size_t rx = fit->rx, ry = fit->ry, npixels = rx * ry;

	switch (bitpix) {
		case 8: {
			uint8_t* m = (uint8_t*) fit->mask->data;
			if (source->type == DATA_USHORT) {
				WORD* c = source->pdata[chan];
				if (source->orig_bitpix == BYTE_IMG) {
					for (size_t i = 0 ; i < npixels ;i++) {
						m[i] = (uint8_t) c[i];
					}
				} else {
					for (size_t i = 0 ; i < npixels ;i++) {
						m[i] = (uint8_t) ((c[i] * 255 + 32895) >> 16);
					}
				}
			} else {
				float *c = source->fpdata[chan];
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = roundf_to_BYTE(c[i] * UCHAR_MAX_SINGLE);
				}
			}
			break;
		}
		case 16: {
			uint16_t* m = (uint16_t*) fit->mask->data;
			if (source->type == DATA_USHORT) {
				uint16_t* c =source->pdata[chan];
				if (source->orig_bitpix == BYTE_IMG) {
					for (size_t i = 0 ; i < npixels ;i++) {
						m[i] = (c[i] << 8) | c[i];
					}
				} else {
					memcpy(fit->mask->data, source->pdata[chan], npixels * sizeof(WORD));
				}
			} else {
				float* c = source->fpdata[chan];
				for (size_t i = 0 ; i < npixels ; i++) {
					m[i] = roundf_to_WORD(c[i] * USHRT_MAX_SINGLE);
				}
			}
			break;
		}
		case 32: {
			if (source->type == DATA_FLOAT) {
				memcpy(fit->mask->data, source->fpdata[chan], npixels * sizeof(float));
			} else {
				float* m = (float*) fit->mask->data;
				WORD* c = source->pdata[chan];
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = c[i] * INV_USHRT_MAX_SINGLE;
				}
			}
			break;
		}
	}
	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP

FAST_MATH_PUSH
int mask_create_from_luminance(fits *fit, fits *source, float rw, float gw, float bw, uint8_t bitpix) {
	if (!fit || !source) return 1;

	if (source->naxes[2] == 1) {
		siril_debug_print("mask_create_from_luminance called on mono image, using mono channel as luminance\n");
		return mask_create_from_channel(fit, source, 0, bitpix);
	}

	set_mask_active(fit, FALSE);
	if (fit->mask) free_mask(fit->mask);

	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) return 1;

	fit->mask->bitpix = bitpix;
	size_t npixels = (size_t)fit->rx * fit->ry;

	// Correct byte allocation: bitpix / 8 gives bytes per pixel
	size_t bytes_per_pixel = bitpix / 8;
	fit->mask->data = malloc(npixels * bytes_per_pixel);

	if (!fit->mask->data) {
		free(fit->mask);
		fit->mask = NULL;
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Pre-calculate scaling factors
	// If weights rw+gw+bw = 1.0, no extra division is needed for float.
	// For USHORT (0-65535) to Float (0-1.0), factor is 1/65535.
	float ushort_to_float = 1.0f / 65535.0f;

	switch (bitpix) {
		case 8: {
			uint8_t* m = (uint8_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float lum;
				if (source->type == DATA_USHORT) {
					// Convert USHORT sum to 0.0-255.0 range
					lum = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw) / 257.0f;
				} else {
					// Float (0.0-1.0) to 0.0-255.0 range
					lum = (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw) * 255.0f;
				}
				m[i] = roundf_to_BYTE(lum);
			}
			break;
		}
		case 16: {
			uint16_t* m = (uint16_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float lum;
				if (source->type == DATA_USHORT) {
					// Source is already 0-65535, just apply weights
					lum = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw);
				} else {
					// Float 0.0-1.0 to 0-65535
					lum = (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw) * 65535.0f;
				}
				m[i] = roundf_to_WORD(lum);
			}
			break;
		}
		case 32: {
			float* m = (float*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				if (source->type == DATA_USHORT) {
					// USHORT to 0.0-1.0
					m[i] = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw) * ushort_to_float;
				} else {
					// Float to Float
					m[i] = (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw);
				}
			}
			break;
		}
	}

	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP

int mask_create_from_luminance_even(fits *fit, fits *source, uint8_t bitpix) {
	float third = 1.f / 3.f;
	return mask_create_from_luminance(fit, source, third, third, third, bitpix);
}

int mask_create_from_luminance_human(fits *fit, fits *source, uint8_t bitpix) {
	return mask_create_from_luminance(fit, source, 0.2126, 0.7152, 0.0722, bitpix);
}


/**
 * Creates a mask based on chromaticity (pure color) and luminance (brightness) ranges
 * This separates color selection from brightness, allowing you to select colors
 * across all brightness levels independently.
 *
 * Chromaticity represents the pure color ratios (r/(r+g+b), g/(r+g+b), b/(r+g+b))
 * Luminance is calculated using standard photometric weights
 *
 * @param fit Output mask structure
 * @param source Source image to analyze
 * @param chrom_center_r Red component of target chromaticity (0.0 - 1.0)
 * @param chrom_center_g Green component of target chromaticity (0.0 - 1.0)
 * @param chrom_center_b Blue component of target chromaticity (0.0 - 1.0)
 *                       Note: these should sum to 1.0 for a valid chromaticity
 * @param chrom_tolerance Maximum distance in chromaticity space to accept (0.0 - 1.0)
 *                        Smaller = more selective, larger = broader color range
 * @param lum_min Minimum luminance (0.0 - 1.0)
 * @param lum_max Maximum luminance (0.0 - 1.0)
 * @param feather_radius Feathering radius in pixels (0 = no feathering)
 * @param invert If TRUE, inverts the mask
 * @param bitpix Output bit depth (8, 16, or 32)
 * @return 0 on success, 1 on failure
 *
 * Example use cases:
 * - Select all red stars regardless of magnitude:
 *   chrom_center=(0.6, 0.2, 0.2), tolerance=0.15, lum_min=0.01, lum_max=1.0
 * - Select bright blue nebula regions:
 *   chrom_center=(0.2, 0.3, 0.5), tolerance=0.2, lum_min=0.3, lum_max=1.0
 * - Select faint H-alpha emission:
 *   chrom_center=(0.7, 0.2, 0.1), tolerance=0.1, lum_min=0.0, lum_max=0.2
 */
FAST_MATH_PUSH
int mask_create_from_chromaticity_luminance(fits *fit, fits *source,
                                            float chrom_center_r, float chrom_center_g, float chrom_center_b,
                                            float chrom_tolerance,
                                            float lum_min, float lum_max,
                                            int feather_radius, gboolean invert,
                                            uint8_t bitpix) {
	if (!fit || !source) return 1;
	if (source->naxes[2] == 1) {
		siril_log_color_message(_("Color mask requires RGB image\n"), "red");
		return 1;
	}

	// Normalize chromaticity center (in case user didn't)
	float chrom_sum = chrom_center_r + chrom_center_g + chrom_center_b;
	if (chrom_sum < 0.001f) {
		siril_log_color_message(_("Invalid chromaticity center (sum near zero)\n"), "red");
		return 1;
	}
	chrom_center_r /= chrom_sum;
	chrom_center_g /= chrom_sum;
	chrom_center_b /= chrom_sum;

	set_mask_active(fit, FALSE);
	if (fit->mask) free_mask(fit->mask);
	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) return 1;
	fit->mask->bitpix = bitpix;

	size_t npixels = (size_t)fit->rx * fit->ry;
	size_t bytes_per_pixel = bitpix / 8;
	fit->mask->data = malloc(npixels * bytes_per_pixel);
	if (!fit->mask->data) {
		free(fit->mask);
		fit->mask = NULL;
		PRINT_ALLOC_ERR;
		return 1;
	}

	float *temp_mask = malloc(npixels * sizeof(float));
	if (!temp_mask) {
		free(fit->mask->data);
		free(fit->mask);
		fit->mask = NULL;
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Standard photometric luminance weights (ITU-R BT.709)
	const float lum_r = 0.2126f;
	const float lum_g = 0.7152f;
	const float lum_b = 0.0722f;

	// Create initial mask based on chromaticity and luminance ranges
	for (size_t i = 0; i < npixels; i++) {
		float r, g, b;

		// Normalize RGB to 0.0-1.0
		if (source->type == DATA_USHORT) {
			r = source->pdata[RLAYER][i] / 65535.0f;
			g = source->pdata[GLAYER][i] / 65535.0f;
			b = source->pdata[BLAYER][i] / 65535.0f;
		} else {
			r = source->fpdata[RLAYER][i];
			g = source->fpdata[GLAYER][i];
			b = source->fpdata[BLAYER][i];
		}

		// Calculate luminance (brightness independent of color)
		float luminance = lum_r * r + lum_g * g + lum_b * b;

		// Calculate chromaticity (color independent of brightness)
		float rgb_sum = r + g + b;
		float chrom_r, chrom_g, chrom_b;

		if (rgb_sum > 0.0001f) {
			// Normal case: pixel has some brightness
			chrom_r = r / rgb_sum;
			chrom_g = g / rgb_sum;
			chrom_b = b / rgb_sum;
		} else {
			// Very dark pixel: chromaticity is undefined
			// Set to neutral gray to avoid matching
			chrom_r = 1.0f / 3.0f;
			chrom_g = 1.0f / 3.0f;
			chrom_b = 1.0f / 3.0f;
		}

		// Calculate chromaticity distance (Euclidean distance in RGB chromaticity space)
		float dr = chrom_r - chrom_center_r;
		float dg = chrom_g - chrom_center_g;
		float db = chrom_b - chrom_center_b;
		float chrom_distance = sqrtf(dr*dr + dg*dg + db*db);

		// Check if pixel is within chromaticity tolerance and luminance range
		gboolean in_chrom_range = (chrom_distance <= chrom_tolerance);
		gboolean in_lum_range = (luminance >= lum_min && luminance <= lum_max);

		// For very dark pixels, we can optionally be more lenient
		// since chromaticity is unreliable near black
		if (rgb_sum < 0.001f && in_lum_range) {
			// Pixel is essentially black - exclude it regardless of "chromaticity"
			temp_mask[i] = 0.0f;
		} else {
			temp_mask[i] = (in_chrom_range && in_lum_range) ? 1.0f : 0.0f;
		}
	}

	// Convert to output format with optional inversion
	switch (bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float val = invert ? (1.0f - temp_mask[i]) : temp_mask[i];
				m[i] = roundf_to_BYTE(val * 255.0f);
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float val = invert ? (1.0f - temp_mask[i]) : temp_mask[i];
				m[i] = roundf_to_WORD(val * 65535.0f);
			}
			break;
		}
		case 32: {
			float *m = (float*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = invert ? (1.0f - temp_mask[i]) : temp_mask[i];
			}
			break;
		}
	}

	// Apply feathering if requested
	if (feather_radius > 0) {
			mask_feather(fit, feather_radius, FEATHER_OUTER);
	}

	free(temp_mask);
	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP

/**
* mask_create_from_image:
* @fit: The target fits structure where the mask will be stored
* @filename: Path to the source image file
* @chan: Channel to use (-1 for luminance, 0 for R, 1 for G, 2 for B)
* @bitpix: Output bit depth (8, 16, or 32)
* @weight_r: Red channel weight for luminance calculation (ignored if chan != -1)
* @weight_g: Green channel weight for luminance calculation (ignored if chan != -1)
* @weight_b: Blue channel weight for luminance calculation (ignored if chan != -1)
*
* Creates a mask from an external image file. If chan is -1, uses luminance
* with the specified weights. Otherwise, uses the specified channel.
*
* Returns: 0 on success, 1 on error
*/
int mask_create_from_image(fits *fit, gchar *filename, int chan, uint8_t bitpix,
						double weight_r, double weight_g, double weight_b) {
	if (!fit || !filename) {
		siril_debug_print("mask_create_from_image: invalid parameters\n");
		return 1;
	}
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32)) {
		siril_debug_print("mask_create_from_image: bitpix must be 8, 16, or 32\n");
		return 1;
	}
	if (chan < -1 || chan > 2) {
		siril_debug_print("mask_create_from_image: chan must be -1 (luminance), 0 (R), 1 (G), or 2 (B)\n");
		return 1;
	}

	// Validate luminance weights if using luminance mode
	if (chan == -1) {
		double weight_sum = weight_r + weight_g + weight_b;
		if (weight_sum < 0.99 || weight_sum > 1.01) {
			siril_log_color_message(_("Warning: luminance weights sum to %.3f (should be 1.0), normalizing\n"),
				"salmon", weight_sum);
			// Normalize weights
			if (weight_sum > 0.0) {
				weight_r /= weight_sum;
				weight_g /= weight_sum;
				weight_b /= weight_sum;
			} else {
				// Use even weights as fallback
				weight_r = weight_g = weight_b = 0.333;
			}
		}
	}

	// Load the source image
	fits *source = calloc(1, sizeof(fits));
	if (!source) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	int retval = readfits(filename, source, FALSE, FALSE);
	if (retval) {
		siril_log_color_message(_("Failed to load mask image: %s\n"), "red", filename);
		free(source);
		return 1;
	}

	// Check dimensions match
	if (fit->rx != source->rx || fit->ry != source->ry) {
		siril_log_color_message(_("Mask image dimensions (%ux%u) do not match target image (%ux%u)\n"),
			"red", source->rx, source->ry, fit->rx, fit->ry);
		clearfits(source);
		free(source);
		return 1;
	}

	// Check channel validity for the source image
	if (chan >= 0 && chan >= source->naxes[2]) {
		siril_log_color_message(_("Channel %d not available in source image (has %u channels)\n"),
			"red", chan, source->naxes[2]);
		clearfits(source);
		free(source);
		return 1;
	}

	if (chan == -1) {
		// Use luminance with specified weights
		siril_log_message(_("Creating mask from luminance of %s (weights: R=%.3f, G=%.3f, B=%.3f)\n"),
			filename, weight_r, weight_g, weight_b);
		retval = mask_create_from_luminance(fit, source, weight_r, weight_g, weight_b, bitpix);
	} else {
		// Use specific channel
		const char *chan_names[] = {"red", "green", "blue"};
		const char *chan_name = (source->naxes[2] == 1) ? "mono" : chan_names[chan];
		siril_log_message(_("Creating mask from %s channel of %s\n"), chan_name, filename);
		retval = mask_create_from_channel(fit, source, chan, bitpix);
	}

	// Clean up
	clearfits(source);
	free(source);

	if (retval == 0) {
		siril_log_message(_("Mask created successfully from image\n"));
	} else {
		siril_log_color_message(_("Failed to create mask from image\n"), "red");
	}

	return retval;
}

#undef max
#define max(a,b) \
({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; })

#undef min
#define min(a,b) \
({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a < _b ? _a : _b; })

FAST_MATH_PUSH
int mask_create_from_stars(fits *fit, float n_fwhm, uint8_t bitpix) {
	if (!fit) return 1;
	if (n_fwhm <= 0.f) {
		siril_debug_print("mask_create_from_stars: n_fwhm must be positive\n");
		return 1;
	}
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32)) {
		siril_debug_print("mask_create_from_stars: bitpix must be 8, 16, or 32\n");
		return 1;
	}

	psf_star **stars = NULL;
	int nb_stars = 0;
	gboolean stars_needs_freeing = FALSE;

	// Check if we already have stars in com.stars
	if (starcount(com.stars) < 1) {
		// Need to detect stars
		struct starfinder_data *sf_data = calloc(1, sizeof(struct starfinder_data));
		if (!sf_data) {
			siril_log_color_message(_("Memory allocation failed\n"), "red");
			return 1;
		}

		sf_data->im.fit = fit;
		sf_data->im.from_seq = NULL;
		sf_data->im.index_in_seq = -1;
		sf_data->layer = (fit->naxes[2] == 1) ? 0 : 1;
		sf_data->max_stars_fitted = MAX_STARS;
		sf_data->selection = (rectangle){0, 0, 0, 0};
		sf_data->save_eqcoords = FALSE;
		sf_data->ref_wcs = NULL;
		sf_data->stars = &stars;
		sf_data->nb_stars = &nb_stars;
		sf_data->threading = MULTI_THREADED;
		sf_data->update_GUI = FALSE;
		sf_data->process_all_images = FALSE;
		sf_data->already_in_thread = FALSE;
		sf_data->keep_stars = FALSE;

		int retval = GPOINTER_TO_INT(findstar_worker(sf_data));
		free(sf_data);

		if (retval != 0 || !stars) {
			siril_log_color_message(_("Star detection failed\n"), "red");
			if (stars)
				free_fitted_stars(stars);
			return 1;
		}
		stars_needs_freeing = TRUE;
	} else {
		stars = com.stars;
		nb_stars = starcount(com.stars);
	}

	if (nb_stars < 1 || !stars) {
		siril_log_color_message(_("No stars detected in the image.\n"), "red");
		if (stars_needs_freeing)
			free_fitted_stars(stars);
		return 1;
	}

	siril_log_message(_("Creating mask from %d stars (n_fwhm = %.2f)...\n"), nb_stars, n_fwhm);

	// Free existing mask if present
	if (fit->mask) {
		set_mask_active(fit, FALSE);
		free_mask(fit->mask);
	}

	// Create new mask
	fit->mask = calloc(1, sizeof(mask_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		if (stars_needs_freeing)
			free_fitted_stars(stars);
		return 1;
	}

	fit->mask->bitpix = bitpix;
	size_t npixels = fit->rx * fit->ry;
	fit->mask->data = calloc(npixels, bitpix);
	if (!fit->mask->data) {
		PRINT_ALLOC_ERR;
		free(fit->mask);
		fit->mask = NULL;
		if (stars_needs_freeing)
			free_fitted_stars(stars);
		return 1;
	}

	int dimx = fit->rx;
	int dimy = fit->ry;

	// For each star, create a disc-shaped mask
	for (int n = 0; n < nb_stars; n++) {
		// Use the average of fwhmx and fwhmy
		float fwhm = (stars[n]->fwhmx + stars[n]->fwhmy) / 2.0f;
		float radius = n_fwhm * fwhm / 2.0f;
		int x = (int)roundf(stars[n]->xpos);
		int y = (int)roundf(stars[n]->ypos);

		// Calculate bounding box
		int xmin = max(0, (int)floorf(x - radius));
		int xmax = min(dimx - 1, (int)ceilf(x + radius));
		int ymin = max(0, (int)floorf(y - radius));
		int ymax = min(dimy - 1, (int)ceilf(y + radius));

		float radius_sq = radius * radius;

		// Fill the disc
		for (int yy = ymin; yy <= ymax; yy++) {
			for (int xx = xmin; xx <= xmax; xx++) {
				float dx = xx - stars[n]->xpos;
				float dy = yy - stars[n]->ypos;
				float dist_sq = dx * dx + dy * dy;

				if (dist_sq <= radius_sq) {
					// Use inverted y-coordinate like synthstar functions do
					size_t idx = xx + (dimy - yy) * dimx;
					switch (bitpix) {
						case 8: {
							uint8_t *m = (uint8_t*)fit->mask->data;
							m[idx] = 255;
							break;
						}
						case 16: {
							uint16_t *m = (uint16_t*)fit->mask->data;
							m[idx] = 65535;
							break;
						}
						case 32: {
							float *m = (float*)fit->mask->data;
							m[idx] = 1.0f;
							break;
						}
					}
				}
			}
		}
	}

	if (stars_needs_freeing)
		free_fitted_stars(stars);

	set_mask_active(fit, TRUE);
	show_or_hide_mask_tab();
	siril_log_message(_("Star mask created successfully.\n"));
	return 0;
}
FAST_MATH_POP

int mask_autostretch(fits *fit) {
	struct mtf_data *data = create_mtf_data();
	if (!data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fits *mfit = mask_to_fits(fit);

	data->fit = mfit;
	data->auto_display_compensation = FALSE;
	data->is_preview = FALSE;
	data->linked = TRUE;

	// Create generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		destroy_mtf_data(data);
		return 1;
	}
	// Compute the autostretch parameters
	find_linked_midtones_balance(mfit, AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, &data->params);
	data->params.do_red = data->params.do_green = data->params.do_blue = TRUE;

	args->fit = mfit;
	args->mem_ratio = 1.0f;
	args->image_hook = mtf_single_image_hook;
	args->description = _("Autostretch mask");
	args->updates_mask = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	// Run worker synchronously - cleanup happens via destructor
	gpointer result = generic_image_worker(args);
	if (result)
		siril_log_color_message(_("Warning: mask autostretch encountered a problem.\n"), "salmon");
	return 0;
}

int get_default_mask_bitpix() {
	return com.pref.default_mask_bitpix > 0 ? com.pref.default_mask_bitpix : gfit->type == DATA_FLOAT ? 32 : gfit->orig_bitpix == BYTE_IMG ? 8 : 16;
}

// Creates a fits struct representing the mask in fit. This allows for masks to be
// processed with the generic_image_worker (using a special idle to re-update the
// mask with the result)
// 8-bit masks are scaled to 16-bit range for better interoperability with processing
// functions, 16-bit and 32-bit masks are converted directly into data or fdata
// These fits structs are always mono (because masks are always mono)

fits *mask_to_fits(fits *fit) {
	if (!fit || !fit->mask || !fit->mask->data)
		return NULL;

	size_t npixels = fit->rx * fit->ry;
	fits *mfit = calloc(1, sizeof(fits));
	copyfits(fit, mfit, CP_FORMAT, -1);
	mfit->naxes[2] = 1;
	mfit->data = NULL;
	mfit->fdata = NULL;
	switch (fit->mask->bitpix) {
		case 8: {
			mfit->data = malloc(fit->rx * fit->ry * sizeof(WORD));
			uint8_t *m = (uint8_t*) fit->mask->data;
			for (size_t i = 0 ; i < npixels ; i++) {
				mfit->data[i] = ((uint16_t) m[i] << 8) | m[i]; // exact scaling to WORD range
			}
			mfit->orig_bitpix = BYTE_IMG; // We set this so we can revert it when converting back
			mfit->bitpix = USHORT_IMG;
			mfit->type = DATA_USHORT;
			break;
		}
		case 16: {
			mfit->data = malloc(fit->rx * fit->ry * sizeof(WORD));
			memcpy(mfit->data, fit->mask->data, npixels * sizeof(WORD));
			mfit->orig_bitpix = USHORT_IMG; // So we know this was always 16-bit when converting back
			mfit->bitpix = USHORT_IMG;
			mfit->type = DATA_USHORT;
			break;
		}
		case 32: {
			mfit->fdata = malloc(fit->rx * fit->ry * sizeof(float));
			memcpy(mfit->fdata, fit->mask->data, npixels * sizeof(float));
			mfit->orig_bitpix = FLOAT_IMG;
			mfit->bitpix = FLOAT_IMG;
			mfit->type = DATA_FLOAT;
			break;
		}
		default: {
			free(mfit);
			return NULL;
		}
	}
	mfit->pdata[0] = mfit->pdata[1] = mfit->pdata[2] = mfit->data;
	mfit->fpdata[0] = mfit->fpdata[1] = mfit->fpdata[2] = mfit->fdata;
	invalidate_stats_from_fit(mfit);
	return mfit;
}

mask_t *fits_to_mask(fits *mfit) {
	if (!mfit || (!mfit->data && !mfit->fdata))
		return NULL;

	size_t npixels = mfit->rx * mfit->ry;

	mask_t *mask = calloc(1, sizeof(mask_t));
	if (!mask)
		return NULL;

	switch (mfit->orig_bitpix) {

		/* Was originally 8-bit, expanded to 16-bit */
		case BYTE_IMG: {
			mask->bitpix = 8;
			mask->data = malloc(npixels * sizeof(uint8_t));
			if (!mask->data) {
				free(mask);
				return NULL;
			}

			uint8_t  *dst = (uint8_t *)mask->data;
			uint16_t *src = (uint16_t *)mfit->data;

			for (size_t i = 0; i < npixels; i++) {
				dst[i] = (uint8_t)((uint32_t) src[i] * 255 + 32895) >> 16;
			}
			break;
		}

		/* Was originally 16-bit */
		case USHORT_IMG: {
			mask->bitpix = 16;
			mask->data = malloc(npixels * sizeof(uint16_t));
			if (!mask->data) {
				free(mask);
				return NULL;
			}

			memcpy(mask->data, mfit->data, npixels * sizeof(uint16_t));
			break;
		}

		/* Was originally 32-bit float */
		case FLOAT_IMG: {
			mask->bitpix = 32;
			mask->data = malloc(npixels * sizeof(float));
			if (!mask->data) {
				free(mask);
				return NULL;
			}

			memcpy(mask->data, mfit->fdata, npixels * sizeof(float));
			break;
		}

		default:
			free(mask);
			return NULL;
	}

	return mask;
}

// Binarize the mask based on min-max range
int mask_binarize(fits *fit, float min_val, float max_val) {
	if (!fit || !fit->mask || !fit->mask->data) {
		siril_debug_print("mask_binarize: invalid mask\n");
		return 1;
	}

	if (min_val > max_val) {
		siril_debug_print("mask_binarize: min_val must be <= max_val\n");
		return 1;
	}

	size_t npixels = fit->rx * fit->ry;
	float actual_min = min_val;
	float actual_max = max_val;

	// Scale range values if they are normalized (< 1) but mask is 8 or 16 bit
	if (min_val < 1.f && max_val < 1.f) {
		if (fit->mask->bitpix == 8) {
			actual_min = min_val * UCHAR_MAX;
			actual_max = max_val * UCHAR_MAX;
		} else if (fit->mask->bitpix == 16) {
			actual_min = min_val * USHRT_MAX;
			actual_max = max_val * USHRT_MAX;
		}
	}

	switch (fit->mask->bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = (m[i] >= actual_min && m[i] <= actual_max) ? UCHAR_MAX : 0;
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = (m[i] >= actual_min && m[i] <= actual_max) ? USHRT_MAX : 0;
			}
			break;
		}
		case 32: {
			float *m = (float*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = (m[i] >= actual_min && m[i] <= actual_max) ? 1.f : 0.f;
			}
			break;
		}
		default:
			siril_debug_print("mask_binarize: unsupported bitpix %d\n", fit->mask->bitpix);
			return 1;
	}

	return 0;
}

// Invert the mask. Returns 0 on success.
int mask_invert(fits *fit) {
	if (!fit->mask || !fit->mask->data)
		return 1;
	uint8_t bitpix = fit->mask->bitpix;
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32))
		return 1;
	size_t npixels = fit->rx * fit->ry;
	switch(bitpix) {
		case 8: {
			uint8_t* restrict m = (uint8_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = ~m[i];
			}
			break;
		}
		case 16: {
			uint16_t* restrict m = (uint16_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				m[i] = ~m[i];
			}
			break;
		}
		case 32: {
			float* restrict m = (float*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float v = 1.f - m[i];
				if (v < 0.f) v = 0.f;
				if (v > 1.f) v = 1.f;
				m[i] = v;
			}
			break;
		}
		default:
			siril_debug_print("Error! Unhandled bitpix in mask_create_ones_like\n");
			return 1;
	}
	return 0;
}

// Scale the mask by factor f. Returns 0 on success.
int mask_scale(fits *fit, float f) {
	if (!fit->mask || !fit->mask->data)
		return 1;
	uint8_t bitpix = fit->mask->bitpix;
	if (!(bitpix == 8 || bitpix == 16 || bitpix == 32))
		return 1;
	size_t npixels = fit->rx * fit->ry;
	switch(bitpix) {
		case 8: {
			uint8_t* restrict m = (uint8_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float v = m[i] * f;
				m[i] = roundf_to_BYTE(v);
			}
			break;
		}
		case 16: {
			uint16_t* restrict m = (uint16_t*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float v = m[i] * f;
				m[i] = roundf_to_WORD(v);
			}
			break;
		}
		case 32: {
			float* restrict m = (float*) fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float v = m[i] * f;
				if (v < 0.f) v = 0.f;
				if (v > 1.f) v = 1.f;
				m[i] = v;
			}
			break;
		}
		default:
			siril_debug_print("Error! Unhandled bitpix in mask_scale\n");
			return 1;
	}
	return 0;
}

int mask_change_bitpix(fits* fit, uint8_t new_bitpix) {
	if (!fit->mask || !fit->mask->data) return 1;

	uint8_t old_bitpix = fit->mask->bitpix;
	size_t npixels = fit->rx * fit->ry;
	if (new_bitpix == old_bitpix) return 0; // nothing to do, but not an error

	switch (old_bitpix) {
		case 8: {
			uint8_t* old_data = (uint8_t*) fit->mask->data;
			switch (new_bitpix) {
				case 16: {
					uint16_t *new_data = malloc(npixels * sizeof(uint16_t));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (uint16_t) (old_data[i] << 8) | old_data[i];
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 16;
					break;
				}
				case 32: {
					float *new_data = malloc(npixels * sizeof(float));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (float) old_data[i] * INV_UCHAR_MAX_SINGLE;
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 32;
					break;
				}
				default: {
					return 1;
				}
			}
			break;
		}
		case 16: {
			uint16_t* old_data = (uint16_t*) fit->mask->data;
			switch (new_bitpix) {
				case 8: {
					uint8_t *new_data = malloc(npixels * sizeof(uint8_t));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (uint8_t)((uint32_t) old_data[i] * 255 + 32895) >> 16;
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 8;
					break;
				}
				case 32: {
					float *new_data = malloc(npixels * sizeof(float));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (float) old_data[i] * INV_USHRT_MAX_SINGLE;
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 32;
					break;
				}
				default: {
					return 1;
				}
			}
			break;
		}
		case 32: {
			float* old_data = (float*) fit->mask->data;
			switch (new_bitpix) {
				case 8: {
					uint8_t *new_data = malloc(npixels * sizeof(uint16_t));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (uint8_t) roundf_to_BYTE(old_data[i] * UCHAR_MAX_SINGLE);
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 8;
					break;
				}
				case 16: {
					float *new_data = malloc(npixels * sizeof(float));
					if (!new_data) {
						PRINT_ALLOC_ERR;
						return 1;
					}
					for (size_t i = 0 ; i < npixels ; i++) {
						new_data[i] = (uint16_t) roundf_to_WORD(old_data[i] * USHRT_MAX_SINGLE);
					}
					free(old_data);
					fit->mask->data = (void*) new_data;
					fit->mask->bitpix = 16;
					break;
				}
				default: {
					return 1;
				}
			}
			break;
		}
		default: {
			return 1;
		}
	}
	return 0;
}
