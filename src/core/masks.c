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

// Create a mask with ones-like property (i.e. full applicability: we
// use 8-bit masks so all values are set to 255). Any existing mask is
// freed and remade. Return 0 on success.

int mask_create_ones_like(fits *fit) {
	if (fit->mask)
		free(fit->mask);

	size_t npixels = fit->rx * fit->ry;
	fit->mask = malloc(fit->rx * fit->ry * sizeof(uint8_t));
	uint8_t *m = fit->mask;
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
	for (size_t i = 0 ; i < npixels ; i++)) {
		m[i] = 255;
	}
	return m;
}

// Create a mask with zeroes-like property. Any existing mask is freed
// and remade. Returns 0 on success.

int mask_create_zeroes_like(fits *fit) {
	if (fit->mask)
		free(fit->mask);
	fit->mask = calloc(fit->rx * fit->ry, sizeof(uint8_t));
	if (!fit->mask) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}

// Invert the mask. Returns 0 on success.

int mask_invert(fits *fit) {
	if (!fit->mask)
		return 1;
	size_t npixels = fit->rx * fit->ry;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
	for (size_t i = 0; i < npixels; i++) {
		fit->mask[i] = ~fit->mask[i];
	}
	return 0;
}
