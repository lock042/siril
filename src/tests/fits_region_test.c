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

/* Mirror property of the region copy pair (core/fits_region.c).
 *
 * The point of these tests is the display-to-FITS vertical flip.  Both
 * directions apply it, so a bug that flips one of them the wrong way still
 * round-trips for a region whose rows happen to be symmetric about the image
 * centre — hence the deliberately off-centre rectangles below, and the
 * separate "crop reads the rows we expect" test which pins the direction
 * rather than just the symmetry. */

#include <criterion/criterion.h>
#include <stdlib.h>

#include "core/siril.h"
#include "core/fits_region.h"
#include "io/image_format_fits.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

/* Fill every channel with a value that identifies its (channel, row, column)
 * uniquely, so a misplaced row or plane cannot pass unnoticed. */
static fits *make_image(uint32_t rx, uint32_t ry, int nchans, data_type type) {
	fits *f = calloc(1, sizeof(fits));
	cr_assert_not_null(f);
	if (new_fit_image(&f, rx, ry, nchans, type))
		cr_assert(0, "new_fit_image failed");

	const size_t npix = (size_t)rx * ry;
	for (int c = 0; c < nchans; c++) {
		for (uint32_t y = 0; y < ry; y++) {
			for (uint32_t x = 0; x < rx; x++) {
				const size_t i = (size_t)c * npix + (size_t)y * rx + x;
				const unsigned v = (unsigned)(c * 100000u + y * 311u + x * 7u + 1u);
				if (type == DATA_FLOAT)
					f->fdata[i] = (float)(v % 65535u) / 65535.0f;
				else
					f->data[i] = (WORD)(v % 65535u);
			}
		}
	}
	return f;
}

static void attach_mask(fits *f, uint8_t bitpix) {
	const size_t npix = (size_t)f->rx * f->ry;
	size_t elem = (bitpix == 8) ? 1 : (bitpix == 16) ? 2 : 4;
	f->mask = malloc(sizeof(mask_t));
	cr_assert_not_null(f->mask);
	f->mask->bitpix = bitpix;
	f->mask->data = malloc(npix * elem);
	cr_assert_not_null(f->mask->data);
	for (size_t i = 0; i < npix; i++) {
		if (bitpix == 8)
			((uint8_t *)f->mask->data)[i] = (uint8_t)(i % 251);
		else if (bitpix == 16)
			((uint16_t *)f->mask->data)[i] = (uint16_t)(i % 65521);
		else
			((float *)f->mask->data)[i] = (float)(i % 1009) / 1009.0f;
	}
	f->mask_active = TRUE;
}

static void free_image(fits *f) {
	clearfits(f);
	free(f);
}

/* Byte-compare the pixel buffers of two same-shaped images. */
static void assert_pixels_equal(const fits *a, const fits *b, const char *what) {
	cr_assert_eq(a->rx, b->rx, "%s: rx", what);
	cr_assert_eq(a->ry, b->ry, "%s: ry", what);
	cr_assert_eq(a->naxes[2], b->naxes[2], "%s: channels", what);
	const size_t n = (size_t)a->rx * a->ry * a->naxes[2];
	if (a->type == DATA_FLOAT) {
		for (size_t i = 0; i < n; i++)
			cr_assert_float_eq(a->fdata[i], b->fdata[i], 0.0f,
					   "%s: float pixel %zu", what, i);
	} else {
		for (size_t i = 0; i < n; i++)
			cr_assert_eq(a->data[i], b->data[i],
				     "%s: ushort pixel %zu", what, i);
	}
}

/* crop then paste back into the same place must leave the image untouched. */
static void check_roundtrip(uint32_t rx, uint32_t ry, int nchans,
			    data_type type, rectangle area, gboolean with_mask) {
	fits *img = make_image(rx, ry, nchans, type);
	fits *ref = make_image(rx, ry, nchans, type);
	if (with_mask)
		attach_mask(img, 16);

	fits region = { 0 };
	cr_assert_eq(crop_fits_region(img, &area, &region), 0, "crop failed");
	cr_assert_eq(region.rx, (uint32_t)area.w, "region width");
	cr_assert_eq(region.ry, (uint32_t)area.h, "region height");
	cr_assert_eq(region.naxes[2], (size_t)nchans, "region channels");
	if (with_mask) {
		cr_assert_not_null(region.mask, "region mask");
		cr_assert_eq(region.mask->bitpix, 16, "region mask depth");
		cr_assert(region.mask_active, "region mask_active");
	}

	cr_assert_eq(paste_fits_region(&region, img, &area), 0, "paste failed");
	assert_pixels_equal(img, ref, "roundtrip");

	clearfits(&region);
	free_image(img);
	free_image(ref);
}

Test(fits_region, roundtrip_is_the_identity_ushort_mono) {
	check_roundtrip(64, 48, 1, DATA_USHORT, (rectangle){ 5, 3, 20, 11 }, FALSE);
}

Test(fits_region, roundtrip_is_the_identity_float_rgb) {
	check_roundtrip(64, 48, 3, DATA_FLOAT, (rectangle){ 5, 3, 20, 11 }, FALSE);
}

Test(fits_region, roundtrip_is_the_identity_with_a_mask) {
	check_roundtrip(64, 48, 3, DATA_USHORT, (rectangle){ 5, 3, 20, 11 }, TRUE);
}

/* Deliberately hard against every edge: a region touching the top-left of the
 * DISPLAY is the bottom-left in storage, which is where an inverted flip shows
 * up first. */
Test(fits_region, roundtrip_is_the_identity_at_the_display_corners) {
	check_roundtrip(37, 29, 1, DATA_USHORT, (rectangle){ 0, 0, 6, 5 }, FALSE);
	check_roundtrip(37, 29, 1, DATA_USHORT, (rectangle){ 31, 24, 6, 5 }, FALSE);
	check_roundtrip(37, 29, 3, DATA_FLOAT, (rectangle){ 0, 24, 6, 5 }, FALSE);
	check_roundtrip(37, 29, 3, DATA_FLOAT, (rectangle){ 31, 0, 6, 5 }, FALSE);
}

Test(fits_region, roundtrip_is_the_identity_for_the_whole_image) {
	check_roundtrip(16, 16, 1, DATA_USHORT, (rectangle){ 0, 0, 16, 16 }, FALSE);
}

/* Pins the flip DIRECTION, which the round-trip alone cannot: the top display
 * row of the region must be the row at (ry - 1 - area.y) in storage. */
Test(fits_region, crop_reads_display_rows_top_down) {
	const uint32_t rx = 24, ry = 20;
	fits *img = make_image(rx, ry, 1, DATA_USHORT);
	rectangle area = { 3, 2, 8, 4 };

	fits region = { 0 };
	cr_assert_eq(crop_fits_region(img, &area, &region), 0);

	for (uint32_t y = 0; y < (uint32_t)area.h; y++) {
		const size_t storage_row = (size_t)ry - y - (size_t)(area.y + 1);
		for (uint32_t x = 0; x < (uint32_t)area.w; x++) {
			cr_assert_eq(region.data[(size_t)y * area.w + x],
				     img->data[storage_row * rx + area.x + x],
				     "row %u col %u", y, x);
		}
	}

	clearfits(&region);
	free_image(img);
}

/* A partially-overlapping or oversized rectangle is a caller error, not
 * something to clamp silently: the ROI path clips through
 * flis_display_to_active_layer_rect before it gets here. */
Test(fits_region, out_of_bounds_areas_are_refused) {
	fits *img = make_image(16, 16, 1, DATA_USHORT);
	fits region = { 0 };

	cr_assert_neq(crop_fits_region(img, &(rectangle){ 10, 0, 10, 4 }, &region), 0,
		      "overhanging right edge accepted");
	cr_assert_neq(crop_fits_region(img, &(rectangle){ 0, 10, 4, 10 }, &region), 0,
		      "overhanging bottom edge accepted");
	cr_assert_neq(crop_fits_region(img, &(rectangle){ -1, 0, 4, 4 }, &region), 0,
		      "negative origin accepted");
	cr_assert_neq(crop_fits_region(img, &(rectangle){ 0, 0, 0, 4 }, &region), 0,
		      "empty area accepted");

	free_image(img);
}

/* The paste side must reject a source that is not the size of the area, rather
 * than reading past the end of it. */
Test(fits_region, paste_refuses_a_mismatched_source) {
	fits *img = make_image(16, 16, 1, DATA_USHORT);
	rectangle area = { 2, 2, 8, 8 };

	fits region = { 0 };
	cr_assert_eq(crop_fits_region(img, &area, &region), 0);

	rectangle bigger = { 2, 2, 9, 8 };
	cr_assert_neq(paste_fits_region(&region, img, &bigger), 0,
		      "size mismatch accepted");

	clearfits(&region);
	free_image(img);
}

/* copy_fits_region is the whole-to-whole case: it must touch exactly the
 * rectangle crop_fits_region would take and leave every other pixel alone.
 *
 * Checking "the region now matches the source" alone would pass for a routine
 * that copied the whole image, and checking "the rest is untouched" alone would
 * pass for one that copied nothing — so both halves are asserted, against a
 * source deliberately made to differ everywhere. */
Test(fits_region, copy_writes_the_region_and_only_the_region) {
	fits *dst = make_image(23, 17, 3, DATA_USHORT);
	fits *src = make_image(23, 17, 3, DATA_USHORT);
	fits *before = make_image(23, 17, 3, DATA_USHORT);
	const size_t n = (size_t)23 * 17 * 3;

	/* make src differ from dst at every pixel, and keep dst's original */
	for (size_t i = 0; i < n; i++)
		src->data[i] = (WORD)((dst->data[i] + 12345u) % 65535u);
	memcpy(before->data, dst->data, n * sizeof(WORD));

	const rectangle area = { 5, 3, 9, 6 };
	cr_assert_eq(copy_fits_region(src, dst, &area), 0);

	/* the rectangle took src's pixels: compare against an independent crop,
	 * so an error in the flip cannot cancel itself out */
	fits want = { 0 }, got = { 0 };
	cr_assert_eq(crop_fits_region(src, &area, &want), 0);
	cr_assert_eq(crop_fits_region(dst, &area, &got), 0);
	cr_assert_eq(memcmp(want.data, got.data, (size_t)9 * 6 * 3 * sizeof(WORD)), 0,
		     "region does not match the source");

	/* everything outside it is untouched */
	size_t changed_outside = 0;
	for (size_t c = 0; c < 3; c++) {
		for (uint32_t fy = 0; fy < 17; fy++) {
			for (uint32_t fx = 0; fx < 23; fx++) {
				/* display row of this FITS row */
				const uint32_t dy = 17 - 1 - fy;
				const gboolean inside = (fx >= 5 && fx < 5 + 9
							 && dy >= 3 && dy < 3 + 6);
				if (inside)
					continue;
				const size_t i = c * (size_t)23 * 17 + (size_t)fy * 23 + fx;
				if (dst->data[i] != before->data[i])
					changed_outside++;
			}
		}
	}
	cr_assert_eq(changed_outside, 0, "%zu pixels outside the area were modified",
		     changed_outside);

	clearfits(&want);
	clearfits(&got);
	free_image(before);
	free_image(src);
	free_image(dst);
}

/* Mismatched geometry is a caller error: this variant has no way to resize. */
Test(fits_region, copy_refuses_mismatched_images) {
	fits *a = make_image(16, 16, 1, DATA_USHORT);
	fits *b = make_image(16, 15, 1, DATA_USHORT);
	fits *c = make_image(16, 16, 1, DATA_FLOAT);
	const rectangle area = { 1, 1, 4, 4 };

	cr_assert_neq(copy_fits_region(a, b, &area), 0, "size mismatch accepted");
	cr_assert_neq(copy_fits_region(a, c, &area), 0, "type mismatch accepted");
	cr_assert_neq(copy_fits_region(a, a, &(rectangle){ 14, 0, 4, 4 }), 0,
		      "out-of-bounds area accepted");

	free_image(a);
	free_image(b);
	free_image(c);
}
