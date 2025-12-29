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
#include "siril.h"


void free_mask(mask_t* mask) {
	free(mask->data);
	free(mask);
}

// Create a test mask : the left half of the image has mask value 255,
// the right half has mask value 0. Any existing mask is
// freed and remade. Return 0 on success.

int mask_create_test(fits *fit, uint8_t bitpix) {
WITH_FAST_MATH
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
	return 0;
}

// Create a mask with ones-like property (i.e. full applicability: we
// use 8-bit masks so all values are set to 255). Any existing mask is
// freed and remade. Return 0 on success.

int mask_create_ones_like(fits *fit, uint8_t bitpix) {
WITH_FAST_MATH
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
	return 0;
}

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
