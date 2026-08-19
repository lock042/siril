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
#ifndef _FITS_REGION_H_
#define _FITS_REGION_H_

/**
 * \file fits_region.h
 * \brief Out-of-place rectangular region copy between two fits.
 *
 * The pair is the primitive under the region-preview path in
 * generic_image_worker, and is intended to serve region-scoped processing more
 * generally — see roi-nde-plan.md.
 *
 * COORDINATE CONVENTION.  \c area is in DISPLAY convention: x/y measured from
 * the TOP-left, which is what gui.roi.selection holds once
 * flis_display_to_active_layer_rect has translated it to layer-local
 * coordinates.  FITS pixel storage is bottom-up, so both functions apply the
 * same vertical flip internally.  They are exact mirrors: for any area within
 * bounds,
 *
 *     crop_fits_region(img, area, tmp);
 *     paste_fits_region(tmp, img, area);
 *
 * leaves img bit-identical.  Get the flip wrong in one of them only and the
 * pair still "works" for a vertically centred region, which is why the mirror
 * property has a unit test rather than a comment.
 *
 * Not to be confused with crop() in algos/geometry.c, which is in-place,
 * GUI-coupled (it clears the ROI and logs), and rounds the rectangle to the
 * CFA pattern.
 */

#include "core/siril.h"   /* fits, rectangle */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * crop_fits_region:
 * @src: source image; not modified.
 * @area: region to extract, display convention, must lie within @src.
 * @dst: destination, cleared and (re)allocated to @area's size.  May be a
 *       freshly zeroed fits or one holding a previous result.
 *
 * Copies scalar metadata (CP_FORMAT), the pixel data of every channel, and the
 * mask if @src carries one (with @src's mask_active).  @dst takes @src's data
 * type; no conversion is performed.
 *
 * Returns: 0 on success, non-zero on bad arguments or allocation failure.
 */
int crop_fits_region(fits *src, const rectangle *area, fits *dst);

/**
 * copy_fits_region:
 * @src: source image; not modified.
 * @dst: destination image, written in place at @area.  Must have the same
 *       dimensions, channel count and pixel type as @src.
 * @area: region to copy, display convention, must lie within both.
 *
 * Copies just @area between two images of the SAME geometry — the same rows
 * crop_fits_region would take, written to the same place they came from.  This
 * is the whole-to-whole case the crop/paste pair cannot express without a
 * region-sized buffer in between; it is what lets a preview restore touch only
 * the rectangle a region preview modified.  Pixels only: no mask, no metadata.
 *
 * Returns: 0 on success, non-zero on bad or mismatched arguments.
 */
int copy_fits_region(const fits *src, fits *dst, const rectangle *area);

/**
 * paste_fits_region:
 * @src: region-sized image, whose dimensions must equal @area's.
 * @dst: destination image, written in place at @area.
 * @area: destination region, display convention, must lie within @dst.
 *
 * Pixel data only — the mask is NOT written back, matching the ROI path this
 * was factored out of: image operations do not edit the mask, so pasting one
 * back would at best be a no-op and at worst overwrite a mask the user changed
 * while a preview was in flight.
 *
 * @src and @dst must share a data type; the caller converts first if not (the
 * ROI path does this itself, since it must convert the preview backup in step).
 * Takes no locks: the caller owns @dst's rwlock.
 *
 * Returns: 0 on success, non-zero on bad arguments.
 */
int paste_fits_region(const fits *src, fits *dst, const rectangle *area);

#ifdef __cplusplus
}
#endif

#endif /* _FITS_REGION_H_ */
