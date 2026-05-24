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
 * test_flis_geometry — exercises the FLIS per-layer geometry helpers
 * that adjust the active layer's canvas offset after a geometry
 * operation has been applied in place to its fits:
 *   - flis_update_layer_offset_after_crop
 *   - flis_update_layer_offset_after_resize
 *   - flis_update_layer_offset_after_rotate
 *   - flis_update_layer_offset_after_mirrorx / _mirrory
 *
 * §7 canvas decoupling: every layer (including the base / bottom-of-
 * stack) behaves identically.  Only the active layer's position is
 * updated.  Canvas-scoped equivalents — which DO touch all layers —
 * live in flis_canvas_* and are covered by their own test file.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_geometry, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Crop on the bottom (formerly "base") layer only repositions that layer;
 * other layers are untouched.  Pre-§7 this case shifted every layer; the
 * new contract is "every layer behaves like any other layer", and canvas-
 * scoped shifts are a separate operation. */
Test(flis_geometry, crop_bottom_layer_only_moves_itself) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;

	uniq_set_active_layer(com.uniq, 0);
	flis_update_layer_offset_after_crop(5, 8);

	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 5,  "cropped layer should land at the selection top-left");
	cr_assert_eq(base->position_y, 8);
	cr_assert_eq(top->position_x, 40, "non-active layer must not move");
	cr_assert_eq(top->position_y, 30);
}

/* Crop on a top layer: only that layer's position changes (to the
 * selection top-left in canvas coords). */
Test(flis_geometry, crop_top_layer_only_moves_itself) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top1 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch1");
	flis_layer_t *top2 = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 1.0f), "patch2");
	top1->position_x = 10; top1->position_y = 15;
	top2->position_x = 50; top2->position_y = 60;

	uniq_set_active_layer(com.uniq, 2);
	flis_update_layer_offset_after_crop(3, 7);

	cr_assert_eq(top2->position_x, 3);
	cr_assert_eq(top2->position_y, 7);
	cr_assert_eq(top1->position_x, 10);
	cr_assert_eq(top1->position_y, 15);
}

/* Resize preserves the active layer's centre on the canvas. */
Test(flis_geometry, resize_preserves_active_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;  /* centre at (50, 40) */
	uniq_set_active_layer(com.uniq, 1);

	/* Layer's fit goes from 20x20 to 40x40 — half-delta is (20-40)/2 = -10
	 * on each axis, so position becomes (40-10, 30-10) = (30, 20).  Centre
	 * stays at (30+20, 20+20) = (50, 40). */
	flis_update_layer_offset_after_resize(20, 20, 40, 40);

	cr_assert_eq(top->position_x, 30);
	cr_assert_eq(top->position_y, 20);
}

/* Rotation preserves the active layer's centre on the canvas. */
Test(flis_geometry, rotate_preserves_active_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 20;  /* centre (25, 25) */
	uniq_set_active_layer(com.uniq, 1);

	/* Layer rotated 180° (no dim change); centre is preserved → position
	 * stays at (20, 20).  The pixel-data rotation has already happened
	 * in the in-place op; only the canvas anchor matters here. */
	flis_update_layer_offset_after_rotate(10, 10, 10, 10, 180.0);

	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 20);
}

/* Group-rotation variant still rotates every non-base layer's centre
 * around the canvas centre.  This is the path the layer-group-rotate
 * UI will use; it's separate from the per-layer helper. */
Test(flis_geometry, rotate_all_layers_90deg) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top1 = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.5f), "p1");
	flis_layer_t *top2 = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.7f), "p2");
	top1->position_x = 20; top1->position_y = 20;
	top2->position_x = 70; top2->position_y = 20;

	flis_update_all_layer_offsets_after_rotate(100, 100, 100, 100, 90.0);

	cr_assert_eq(top1->position_x, 70);
	cr_assert_eq(top1->position_y, 20);
	cr_assert_eq(top2->position_x, 70);
	cr_assert_eq(top2->position_y, 70);
}

/* Mirroring a layer about its own centre is position-invariant. */
Test(flis_geometry, mirror_does_not_move_active) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 1);

	flis_update_layer_offset_after_mirrorx();
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 30);

	flis_update_layer_offset_after_mirrory();
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 30);
}

/* Sanity guard: helpers are no-ops when no FLIS is loaded. */
Test(flis_geometry, helpers_are_noops_without_flis) {
	flis_update_layer_offset_after_crop(10, 10);
	flis_update_layer_offset_after_resize(100, 100, 200, 200);
	flis_update_layer_offset_after_rotate(100, 100, 100, 100, 45.0);
	flis_update_all_layer_offsets_after_rotate(100, 100, 100, 100, 45.0);
	flis_update_layer_offset_after_mirrorx();
	flis_update_layer_offset_after_mirrory();
}
