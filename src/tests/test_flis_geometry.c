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
 */

/*
 * test_flis_geometry — exercises the FLIS geometry-helper functions
 * that adjust layer offsets after geometry operations on gfit:
 *   - flis_update_layer_offset_after_crop
 *   - flis_update_layer_offset_after_resize
 *   - flis_update_layer_offset_after_rotate
 *   - flis_update_all_layer_offsets_after_rotate
 *
 * These functions don't touch pixels — they update each non-active
 * layer's position_x/y so that after a geometry op on the active layer
 * (or canvas), the remaining layers stay in the right place.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_geometry, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Crop of the BASE layer: every non-base layer's position should be
 * reduced by the crop selection origin (so layers that were at canvas
 * (10, 10) become (10-sel_x, 10-sel_y) on the new smaller canvas). */
Test(flis_geometry, crop_base_layer_shifts_non_base_offsets) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top1 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch1");
	flis_layer_t *top2 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 1.0f), "patch2");
	top1->position_x = 10; top1->position_y = 15;
	top2->position_x = 50; top2->position_y = 60;

	/* simulate cropping the base by selection (5, 8, …) — only the offset
	 * matters for this helper.  The active layer here is the most recently
	 * added (top2), so make sure we switch to the base before calling. */
	uniq_set_active_layer(com.uniq, 0);
	flis_update_layer_offset_after_crop(5, 8);

	cr_assert_eq(top1->position_x, 10 - 5);
	cr_assert_eq(top1->position_y, 15 - 8);
	cr_assert_eq(top2->position_x, 50 - 5);
	cr_assert_eq(top2->position_y, 60 - 8);
}

/* Crop of a NON-BASE layer: only that layer's position changes (to the
 * absolute selection origin in canvas space). */
Test(flis_geometry, crop_non_base_sets_only_that_layers_offset) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top1 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch1");
	flis_layer_t *top2 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 1.0f), "patch2");
	top1->position_x = 10; top1->position_y = 15;
	top2->position_x = 50; top2->position_y = 60;

	/* top2 is the active layer (added last) — crop it */
	uniq_set_active_layer(com.uniq, 2);
	flis_update_layer_offset_after_crop(3, 7);

	cr_assert_eq(top2->position_x, 3);
	cr_assert_eq(top2->position_y, 7);
	/* top1 untouched */
	cr_assert_eq(top1->position_x, 10);
	cr_assert_eq(top1->position_y, 15);
}

/* Resize of the BASE layer: non-base layer centres scale proportionally.
 * A layer whose centre was at (cx, cy) on an old_rx × old_ry canvas
 * becomes centred at (cx * new_rx / old_rx, cy * new_ry / old_ry). */
Test(flis_geometry, resize_base_scales_non_base_centres) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;  /* centre at (50, 40) */
	uniq_set_active_layer(com.uniq, 0);
	flis_update_layer_offset_after_resize(100, 100, 200, 200);

	/* New centre should be at (100, 80); layer is still 20x20 so the
	 * new position is centre - (10, 10) = (90, 70). */
	cr_assert_eq(top->position_x, 90);
	cr_assert_eq(top->position_y, 70);
}

/* Rotate of the BASE layer by 180°: every non-base centre is flipped
 * around the canvas centre.  A layer whose centre was at (25, 25) on a
 * 100x100 canvas becomes centred at (75, 75) on the new 100x100 canvas. */
Test(flis_geometry, rotate_base_180_flips_offsets_through_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 20;  /* centre at (25, 25) */
	uniq_set_active_layer(com.uniq, 0);
	flis_update_layer_offset_after_rotate(100, 100, 100, 100, 180.0);

	/* (25, 25) rotated 180° around (50, 50) → (75, 75); offset is (70, 70). */
	cr_assert_eq(top->position_x, 70);
	cr_assert_eq(top->position_y, 70);
}

/* Group-level rotation: every non-base layer is rotated, regardless of
 * which one is currently active.  Verifies the _all_ variant. */
Test(flis_geometry, rotate_all_layers_90deg) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top1 = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.5f), "p1");
	flis_layer_t *top2 = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.7f), "p2");
	top1->position_x = 20; top1->position_y = 20;
	top2->position_x = 70; top2->position_y = 20;

	/* The implementation applies a standard rotation matrix
	 *   (new_dx, new_dy) = (dx cos − dy sin, dx sin + dy cos)
	 * to each layer's centre offset from the canvas centre.  At +90°
	 * (cos=0, sin=1) the matrix is (−dy, dx), so:
	 *   top1 centre (25,25): dx=−25, dy=−25 → new_d=(25, −25) → new centre (75, 25), pos (70, 20)
	 *   top2 centre (75,25): dx=+25, dy=−25 → new_d=(25, +25) → new centre (75, 75), pos (70, 70)
	 * (Positive angle is math-CCW in coordinate-system convention, which
	 * corresponds to clockwise in image space where y grows downward.) */
	flis_update_all_layer_offsets_after_rotate(100, 100, 100, 100, 90.0);

	cr_assert_eq(top1->position_x, 70);
	cr_assert_eq(top1->position_y, 20);
	cr_assert_eq(top2->position_x, 70);
	cr_assert_eq(top2->position_y, 70);
}

/* Sanity guard: helpers are no-ops when no FLIS is loaded. */
Test(flis_geometry, helpers_are_noops_without_flis) {
	/* com.uniq->layers is NULL — should not crash */
	flis_update_layer_offset_after_crop(10, 10);
	flis_update_layer_offset_after_resize(100, 100, 200, 200);
	flis_update_layer_offset_after_rotate(100, 100, 100, 100, 45.0);
	flis_update_all_layer_offsets_after_rotate(100, 100, 100, 100, 45.0);
}
