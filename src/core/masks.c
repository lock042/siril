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
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "core/siril_log.h"
#include "filters/synthstar.h"
#include "gui/callbacks.h"
#include "io/image_format_fits.h"
#include "masks.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include <opencv2/imgproc/imgproc_c.h>
#pragma GCC diagnostic pop

void free_mask(mask_t* mask) {
	free(mask->data);
	free(mask);
	show_or_hide_mask_tab();
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
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = (uint8_t) ((c[i] * 255 + 32895) >> 16);
				}
			} else {
				float *c = source->fpdata[chan];
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = roundf_to_BYTE(c[i]);
				}
			}
			break;
		}
		case 16: {
			if (source->type == DATA_USHORT) {
				memcpy(fit->mask->data, source->pdata[chan], npixels * sizeof(WORD));
			} else {
				uint16_t* m = (uint16_t*) fit->mask->data;
				WORD* c = source->pdata[chan];
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = roundf_to_WORD(c[i]);
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
	show_or_hide_mask_tab();
	return 0;
}
FAST_MATH_POP
int mask_create_from_luminance(fits *fit, fits *source, float rw, float gw, float bw, uint8_t bitpix) {
	if (!fit) return 1; // no FITS struct
	if (!source) return 1; // no source FITS struct

	// Handle mono
	if (source->naxes[2] == 1) {
		siril_debug_print("mask_create_from_luminance called on mono image, using mono channel as luminance\n");
		return mask_create_from_channel(fit, source, 0, bitpix);
	}

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
				float factor = 1.f / (source->naxes[2] * USHRT_MAX_SINGLE);
				for (size_t i = 0 ; i < npixels ;i++) {
					float lum = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw) / factor;
					m[i] = roundf_to_BYTE(lum);
				}
			} else {
				float third = 1.f / 3.f;
				for (size_t i = 0 ; i < npixels ;i++) {
					float lum = (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw) * third;
					m[i] = roundf_to_BYTE(lum);
				}
			}
			break;
		}
		case 16: {
			uint16_t* m = (uint16_t*) fit->mask->data;
			if (fit->type == DATA_USHORT) {
				float factor = 1.f / (source->naxes[2] * USHRT_MAX_SINGLE);
				for (size_t i = 0 ; i < npixels ;i++) {
					float lum = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw) * factor;
					m[i] = roundf_to_WORD(lum);
				}
			} else {
				float third = 1.f / 3.f;
				for (size_t i = 0 ; i < npixels ;i++) {
					float lum = (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw) * third;
					m[i] = roundf_to_WORD(lum);
				}
			}
			break;
		}
		case 32: {
			float* m = (float*) fit->mask->data;
			if (fit->type == DATA_USHORT) {
				float factor = 1.f / (source->naxes[2] * USHRT_MAX_SINGLE);
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i] = (source->pdata[RLAYER][i] * rw + source->pdata[GLAYER][i] * gw + source->pdata[BLAYER][i] * bw) * factor;
				}
			} else {
				float third = 1.f / 3.f;
				for (size_t i = 0 ; i < npixels ;i++) {
					m[i]= (source->fpdata[RLAYER][i] * rw + source->fpdata[GLAYER][i] * gw + source->fpdata[BLAYER][i] * bw) * third;
				}
			}
			break;
		}
	}
	show_or_hide_mask_tab();
	return 0;
}

int mask_create_from_luminance_even(fits *fit, fits *source, uint8_t bitpix) {
	float third = 1.f / 3.f;
	return mask_create_from_luminance(fit, source, third, third, third, bitpix);
}

int mask_create_from_luminance_human(fits *fit, fits *source, uint8_t bitpix) {
	return mask_create_from_luminance(fit, source, 0.2126, 0.7152, 0.0722, bitpix);
}

// Apply Gaussian blur to the mask
int mask_apply_gaussian_blur(fits *fit, float radius) {
	if (!fit || !fit->mask || !fit->mask->data) {
		siril_debug_print("mask_apply_gaussian_blur: invalid mask\n");
		return 1;
	}

	if (radius <= 0.f) {
		siril_debug_print("mask_apply_gaussian_blur: radius must be positive\n");
		return 1;
	}

	size_t rx = fit->rx, ry = fit->ry;
	int cv_type;

	// Determine OpenCV type based on bitpix
	switch (fit->mask->bitpix) {
		case 8:
			cv_type = CV_8UC1;
			break;
		case 16:
			cv_type = CV_16UC1;
			break;
		case 32:
			cv_type = CV_32FC1;
			break;
		default:
			siril_debug_print("mask_apply_gaussian_blur: unsupported bitpix %d\n", fit->mask->bitpix);
			return 1;
	}

	// Create OpenCV Mat headers for source and destination
	CvMat *src = cvCreateMatHeader(ry, rx, cv_type);
	CvMat *dst = cvCreateMatHeader(ry, rx, cv_type);

	if (!src || !dst) {
		if (src) cvReleaseMat(&src);
		if (dst) cvReleaseMat(&dst);
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Allocate destination data
	void *blur_data = malloc(rx * ry * fit->mask->bitpix);
	if (!blur_data) {
		cvReleaseMat(&src);
		cvReleaseMat(&dst);
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Set data pointers
	cvSetData(src, fit->mask->data, rx * (fit->mask->bitpix / 8));
	cvSetData(dst, blur_data, rx * (fit->mask->bitpix / 8));

	// Calculate kernel size from radius (should be odd)
	// Using approximation: kernel_size = 2 * ceil(3 * sigma) + 1, where sigma = radius
	int kernel_size = 2 * (int)ceilf(3.f * radius) + 1;

	// Apply Gaussian blur
	cvSmooth(src, dst, CV_GAUSSIAN, kernel_size, kernel_size, radius, radius);

	// Replace original data with blurred data
	free(fit->mask->data);
	fit->mask->data = blur_data;

	cvReleaseMat(&src);
	cvReleaseMat(&dst);

	return 0;
}

// Create a mask from an external image file
// filename: path to the image file to load
// chan: channel to use (0, 1, 2 for R, G, B or mono; -1 for luminance)
// bitpix: 8, 16, or 32 for mask bit depth
int mask_create_from_image(fits *fit, gchar *filename, int chan, uint8_t bitpix) {
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

	// Load the source image
	fits *source = NULL;
	int retval = readfits(filename, source, FALSE, FALSE);
	if (retval) {
		siril_log_color_message(_("Failed to load mask image: %s\n"), "red", filename);
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
		// Use luminance
		siril_log_message(_("Creating mask from luminance of %s\n"), filename);
		retval = mask_create_from_luminance_human(fit, source, bitpix);
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
	if (fit->mask)
		free_mask(fit->mask);

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

	show_or_hide_mask_tab();
	siril_log_message(_("Star mask created successfully.\n"));
	return 0;
}
FAST_MATH_POP

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

// Feather (soften edges) of a mask using distance transform
// feather_dist: distance in pixels over which to feather from the edge
// mode: FEATHER_INNER, FEATHER_OUTER, or FEATHER_EDGE
// Returns 0 on success
int mask_feather(fits *fit, float feather_dist, feather_mode mode) {
	if (!fit || !fit->mask || !fit->mask->data) {
		siril_debug_print("mask_feather: invalid mask\n");
		return 1;
	}

	if (feather_dist <= 0.f) {
		siril_debug_print("mask_feather: feather_dist must be positive\n");
		return 1;
	}

	size_t rx = fit->rx, ry = fit->ry;
	size_t npixels = rx * ry;
	uint8_t bitpix = fit->mask->bitpix;

	// Convert mask to binary uint8 for OpenCV processing
	uint8_t *binary = malloc(npixels * sizeof(uint8_t));
	if (!binary) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Convert mask to binary uint8
	switch (bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary[i] = (m[i] > 127) ? 255 : 0;
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary[i] = (m[i] > 32767) ? 255 : 0;
			}
			break;
		}
		case 32: {
			float *m = (float*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary[i] = (m[i] > 0.5f) ? 255 : 0;
			}
			break;
		}
	}

	// Determine which distance transforms we need
	gboolean need_inside = (mode == FEATHER_INNER || mode == FEATHER_EDGE);
	gboolean need_outside = (mode == FEATHER_OUTER || mode == FEATHER_EDGE);

	// Adjust feather distance for EDGE mode
	float effective_dist = (mode == FEATHER_EDGE) ? feather_dist / 2.0f : feather_dist;

	// Allocate distance maps
	float *dist_inside = NULL;
	float *dist_outside = NULL;

	if (need_inside) {
		dist_inside = malloc(npixels * sizeof(float));
		if (!dist_inside) {
			free(binary);
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	if (need_outside) {
		dist_outside = malloc(npixels * sizeof(float));
		if (!dist_outside) {
			free(binary);
			if (dist_inside) free(dist_inside);
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	// Create OpenCV matrices
	CvMat *src_inside = NULL, *dst_inside = NULL;
	CvMat *src_outside = NULL, *dst_outside = NULL;

	// Compute inside distance transform if needed
	if (need_inside) {
		src_inside = cvCreateMatHeader(ry, rx, CV_8UC1);
		dst_inside = cvCreateMatHeader(ry, rx, CV_32FC1);

		if (!src_inside || !dst_inside) {
			if (src_inside) cvReleaseMat(&src_inside);
			if (dst_inside) cvReleaseMat(&dst_inside);
			free(binary);
			if (dist_inside) free(dist_inside);
			if (dist_outside) free(dist_outside);
			PRINT_ALLOC_ERR;
			return 1;
		}

		cvSetData(src_inside, binary, rx);
		cvSetData(dst_inside, dist_inside, rx * sizeof(float));

		// Distance transform on the mask (inside pixels)
		// cvDistTransform signature: (src, dst, distance_type, mask_size, mask, labels, labelType)
		cvDistTransform(src_inside, dst_inside, CV_DIST_L2, CV_DIST_MASK_PRECISE, NULL, NULL, CV_DIST_LABEL_CCOMP);

		cvReleaseMat(&src_inside);
		cvReleaseMat(&dst_inside);
	}

	// Compute outside distance transform if needed
	if (need_outside) {
		// Invert binary mask for outside distance
		uint8_t *binary_inv = malloc(npixels * sizeof(uint8_t));
		if (!binary_inv) {
			free(binary);
			if (dist_inside) free(dist_inside);
			if (dist_outside) free(dist_outside);
			PRINT_ALLOC_ERR;
			return 1;
		}

		for (size_t i = 0; i < npixels; i++) {
			binary_inv[i] = (binary[i] == 0) ? 255 : 0;
		}

		src_outside = cvCreateMatHeader(ry, rx, CV_8UC1);
		dst_outside = cvCreateMatHeader(ry, rx, CV_32FC1);

		if (!src_outside || !dst_outside) {
			if (src_outside) cvReleaseMat(&src_outside);
			if (dst_outside) cvReleaseMat(&dst_outside);
			free(binary);
			free(binary_inv);
			if (dist_inside) free(dist_inside);
			if (dist_outside) free(dist_outside);
			PRINT_ALLOC_ERR;
			return 1;
		}

		cvSetData(src_outside, binary_inv, rx);
		cvSetData(dst_outside, dist_outside, rx * sizeof(float));

		// Distance transform on inverted mask (outside pixels)
		// cvDistTransform signature: (src, dst, distance_type, mask_size, mask, labels, labelType)
		cvDistTransform(src_outside, dst_outside, CV_DIST_L2, CV_DIST_MASK_PRECISE, NULL, NULL, CV_DIST_LABEL_CCOMP);

		cvReleaseMat(&src_outside);
		cvReleaseMat(&dst_outside);
		free(binary_inv);
	}

	// Apply feathering based on mode
	switch (bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside[i] >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside[i] / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary[i] > 127) {
						value = 1.0f;
					} else if (dist_outside[i] >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside[i] / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary[i] > 127) {
						// Inside: feather inward
						if (dist_inside[i] >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside[i] / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside[i] >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside[i] / feather_dist);
						}
					}
				}
				m[i] = (uint8_t)roundf(value * 255.0f);
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside[i] >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside[i] / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary[i] > 127) {
						value = 1.0f;
					} else if (dist_outside[i] >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside[i] / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary[i] > 127) {
						// Inside: feather inward
						if (dist_inside[i] >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside[i] / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside[i] >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside[i] / feather_dist);
						}
					}
				}
				m[i] = (uint16_t)roundf(value * 65535.0f);
			}
			break;
		}
		case 32: {
			float *m = (float*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside[i] >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside[i] / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary[i] > 127) {
						value = 1.0f;
					} else if (dist_outside[i] >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside[i] / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary[i] > 127) {
						// Inside: feather inward
						if (dist_inside[i] >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside[i] / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside[i] >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside[i] / feather_dist);
						}
					}
				}
				m[i] = value;
			}
			break;
		}
	}

	// Cleanup
	free(binary);
	if (dist_inside) free(dist_inside);
	if (dist_outside) free(dist_outside);

	const char *mode_str = (mode == FEATHER_INNER) ? "inward" :
	                       (mode == FEATHER_OUTER) ? "outward" : "on edge";
	siril_log_message(_("Mask feathered %s with distance %.1f pixels\n"), mode_str, feather_dist);
	return 0;
}
