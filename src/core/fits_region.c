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
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"

#include "fits_region.h"

/* Bytes per mask element for a given mask bitpix, 0 for an unsupported one. */
static size_t mask_elem_size(uint8_t bitpix) {
	switch (bitpix) {
	case 8:  return sizeof(uint8_t);
	case 16: return sizeof(uint16_t);
	case 32: return sizeof(float);
	default: return 0;
	}
}

/* TRUE when area is non-empty and wholly inside a rx x ry image. */
static gboolean area_within(const rectangle *area, uint32_t rx, uint32_t ry) {
	if (!area || area->w <= 0 || area->h <= 0)
		return FALSE;
	if (area->x < 0 || area->y < 0)
		return FALSE;
	return (uint32_t)(area->x + area->w) <= rx
	    && (uint32_t)(area->y + area->h) <= ry;
}

/* Row y of the region (0 = its top row, display convention) maps to this row
 * of a bottom-up fits of height ry.  The single place the flip is expressed;
 * both directions call it, which is what keeps them exact mirrors. */
static inline size_t src_row_of(uint32_t ry, int area_y, uint32_t y) {
	return (size_t)ry - y - (size_t)(area_y + 1);
}

int crop_fits_region(fits *src, const rectangle *area, fits *dst) {
	if (!src || !dst || !area)
		return 1;
	if (!area_within(area, src->rx, src->ry)) {
		siril_log_debug("crop_fits_region: area outside the image\n");
		return 1;
	}
	if (!src->data && !src->fdata)
		return 1;

	const size_t nchans = src->naxes[2];
	if (nchans != 1 && nchans != 3)
		return 1;
	const gboolean rgb = (nchans == 3);
	const size_t npix_dst = (size_t)area->w * area->h;
	const size_t npix_src = (size_t)src->rx * src->ry;

	clearfits(dst);
	copyfits(src, dst, CP_FORMAT, -1);

	dst->rx = dst->naxes[0] = area->w;
	dst->ry = dst->naxes[1] = area->h;
	dst->naxes[2] = nchans;
	dst->naxis = (nchans == 1 ? 2 : 3);

	if (dst->type == DATA_FLOAT) {
		dst->fdata = malloc(npix_dst * nchans * sizeof(float));
		if (!dst->fdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		dst->fpdata[0] = dst->fdata;
		dst->fpdata[1] = rgb ? dst->fdata + npix_dst : dst->fdata;
		dst->fpdata[2] = rgb ? dst->fdata + 2 * npix_dst : dst->fdata;

		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const float *from = src->fdata + (npix_src * c)
					+ (src_row_of(src->ry, area->y, y) * src->rx) + area->x;
				float *to = dst->fdata + (npix_dst * c) + ((size_t)dst->rx * y);
				memcpy(to, from, (size_t)area->w * sizeof(float));
			}
		}
	} else {
		dst->data = malloc(npix_dst * nchans * sizeof(WORD));
		if (!dst->data) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		dst->pdata[0] = dst->data;
		dst->pdata[1] = rgb ? dst->data + npix_dst : dst->data;
		dst->pdata[2] = rgb ? dst->data + 2 * npix_dst : dst->data;

		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const WORD *from = src->data + (npix_src * c)
					+ (src_row_of(src->ry, area->y, y) * src->rx) + area->x;
				WORD *to = dst->data + (npix_dst * c) + ((size_t)dst->rx * y);
				memcpy(to, from, (size_t)area->w * sizeof(WORD));
			}
		}
	}

	/* Mask.  Absent or of an unsupported depth leaves dst unmasked rather
	 * than failing the crop: the pixels are still valid, and that is what the
	 * ROI path did before this was factored out. */
	dst->mask_active = FALSE;
	dst->mask = NULL;

	if (src->mask && src->mask->data) {
		const size_t elem = mask_elem_size(src->mask->bitpix);
		if (!elem)
			return 0;

		dst->mask = malloc(sizeof(mask_t));
		if (!dst->mask) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		dst->mask->bitpix = src->mask->bitpix;
		dst->mask->data = malloc(npix_dst * elem);
		if (!dst->mask->data) {
			PRINT_ALLOC_ERR;
			free(dst->mask);
			dst->mask = NULL;
			return 1;
		}

		for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
			const uint8_t *from = (const uint8_t *)src->mask->data
				+ (src_row_of(src->ry, area->y, y) * src->rx + area->x) * elem;
			uint8_t *to = (uint8_t *)dst->mask->data + ((size_t)area->w * y) * elem;
			memcpy(to, from, (size_t)area->w * elem);
		}
		dst->mask_active = src->mask_active;
	}

	return 0;
}

int copy_fits_region(const fits *src, fits *dst, const rectangle *area) {
	if (!src || !dst || !area)
		return 1;
	if (src->rx != dst->rx || src->ry != dst->ry)
		return 1;
	if (src->type != dst->type || src->naxes[2] != dst->naxes[2])
		return 1;
	if (!area_within(area, dst->rx, dst->ry)) {
		siril_log_debug("copy_fits_region: area outside the image\n");
		return 1;
	}

	const size_t nchans = src->naxes[2];
	const size_t npix = (size_t)src->rx * src->ry;

	if (src->type == DATA_FLOAT) {
		if (!src->fdata || !dst->fdata)
			return 1;
		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const size_t off = (npix * c)
					+ (src_row_of(src->ry, area->y, y) * src->rx) + area->x;
				memcpy(dst->fdata + off, src->fdata + off,
				       (size_t)area->w * sizeof(float));
			}
		}
	} else {
		if (!src->data || !dst->data)
			return 1;
		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const size_t off = (npix * c)
					+ (src_row_of(src->ry, area->y, y) * src->rx) + area->x;
				memcpy(dst->data + off, src->data + off,
				       (size_t)area->w * sizeof(WORD));
			}
		}
	}
	return 0;
}

int paste_fits_region(const fits *src, fits *dst, const rectangle *area) {
	if (!src || !dst || !area)
		return 1;
	if (!area_within(area, dst->rx, dst->ry)) {
		siril_log_debug("paste_fits_region: area outside the image\n");
		return 1;
	}
	if (src->rx != (uint32_t)area->w || src->ry != (uint32_t)area->h) {
		siril_log_debug("paste_fits_region: source does not match the area\n");
		return 1;
	}
	if (src->type != dst->type)
		return 1;
	if (src->naxes[2] != dst->naxes[2])
		return 1;

	const size_t nchans = src->naxes[2];
	const size_t npix_src = (size_t)src->rx * src->ry;
	const size_t npix_dst = (size_t)dst->rx * dst->ry;

	if (src->type == DATA_FLOAT) {
		if (!src->fdata || !dst->fdata)
			return 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const float *from = src->fdata + ((size_t)src->rx * y) + (npix_src * c);
				float *to = dst->fdata + (npix_dst * c)
					+ (src_row_of(dst->ry, area->y, y) * dst->rx) + area->x;
				memcpy(to, from, (size_t)area->w * sizeof(float));
			}
		}
	} else {
		if (!src->data || !dst->data)
			return 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (size_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < (uint32_t)area->h; y++) {
				const WORD *from = src->data + ((size_t)src->rx * y) + (npix_src * c);
				WORD *to = dst->data + (npix_dst * c)
					+ (src_row_of(dst->ry, area->y, y) * dst->rx) + area->x;
				memcpy(to, from, (size_t)area->w * sizeof(WORD));
			}
		}
	}

	return 0;
}
