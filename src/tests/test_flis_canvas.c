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
 * test_flis_canvas — exercises the §7 canvas-scoped helpers
 * (flis_canvas_resize / _fit_to_layers / _rotate / _mirrorx / _mirrory)
 * and the canvas state fields on com.uniq populated by
 * flis_promote_from_gfit + load_flis (round-trip + missing-keys
 * leniency).
 *
 * Helpers act on com.uniq->canvas_w/h AND on every layer's position;
 * pixel data is never modified.  The save/load round-trip is covered
 * separately by test_flis_roundtrip (which already loads/saves files);
 * here we focus on the canvas-specific bookkeeping.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_canvas, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Helper: build a 2-layer fixture with known canvas + positions. */
static void make_two_layer_fixture(guint canvas_w, guint canvas_h) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;
	com.uniq->canvas_w = canvas_w;
	com.uniq->canvas_h = canvas_h;
}

Test(flis_canvas, resize_changes_dims_and_shifts_layers) {
	make_two_layer_fixture(100, 100);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *top  = (flis_layer_t *)com.uniq->layers->next->data;

	cr_assert_eq(flis_canvas_resize(200, 150, 10, 20), 0);
	cr_assert_eq(com.uniq->canvas_w, 200u);
	cr_assert_eq(com.uniq->canvas_h, 150u);
	cr_assert_eq(base->position_x, 10);   /* was 0 + dx */
	cr_assert_eq(base->position_y, 20);
	cr_assert_eq(top->position_x, 50);    /* was 40 + dx */
	cr_assert_eq(top->position_y, 50);    /* was 30 + dy */
}

Test(flis_canvas, resize_rejects_zero_dims) {
	make_two_layer_fixture(100, 100);
	cr_assert_neq(flis_canvas_resize(0, 100, 0, 0), 0);
	cr_assert_neq(flis_canvas_resize(100, 0, 0, 0), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 100u);
}

Test(flis_canvas, fit_to_layers_shrinks_to_bbox) {
	make_two_layer_fixture(500, 500);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *top  = (flis_layer_t *)com.uniq->layers->next->data;
	/* base 100x100 at (0,0); top 20x20 at (40,30).  bbox = [0,60)x[0,100). */

	cr_assert_eq(flis_canvas_fit_to_layers(FALSE), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 100u);
	cr_assert_eq(base->position_x, 0);
	cr_assert_eq(base->position_y, 0);
	cr_assert_eq(top->position_x, 40);
	cr_assert_eq(top->position_y, 30);
}

Test(flis_canvas, fit_to_layers_handles_layers_extending_past_canvas) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(40, 40, 0.5f), "patch");
	top->position_x = -10; top->position_y = -20;  /* extends past 0,0 */
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	cr_assert_eq(flis_canvas_fit_to_layers(FALSE), 0);
	/* bbox: x [-10, 100), y [-20, 100) → 110 x 120; layers shift by (+10,+20). */
	cr_assert_eq(com.uniq->canvas_w, 110u);
	cr_assert_eq(com.uniq->canvas_h, 120u);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 10);
	cr_assert_eq(base->position_y, 20);
	cr_assert_eq(top->position_x, 0);
	cr_assert_eq(top->position_y, 0);
}

Test(flis_canvas, fit_to_layers_skips_invisible_by_default) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *huge = flis_test_add_layer(
	    flis_test_make_mono_fits(500, 500, 0.5f), "huge_invisible");
	huge->position_x = 0; huge->position_y = 0;
	huge->visible = FALSE;
	(void)base;
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	cr_assert_eq(flis_canvas_fit_to_layers(FALSE), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);  /* huge layer ignored */
	cr_assert_eq(com.uniq->canvas_h, 100u);

	cr_assert_eq(flis_canvas_fit_to_layers(TRUE), 0);
	cr_assert_eq(com.uniq->canvas_w, 500u);
	cr_assert_eq(com.uniq->canvas_h, 500u);
}

Test(flis_canvas, rotate_90_swaps_dims_and_rotates_layers) {
	flis_test_add_layer(flis_test_make_mono_fits(40, 40, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "p");
	top->position_x = 60; top->position_y = 20;  /* centre (65, 25) */
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	cr_assert_eq(flis_canvas_rotate(90.0), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);  /* square stays square */
	cr_assert_eq(com.uniq->canvas_h, 100u);
	/* (65,25) → centre (50,50) + R90·(15,-25) = (50,50) + (25,15) = (75,65);
	 * pos = (75-5, 65-5) = (70, 60). */
	cr_assert_eq(top->position_x, 70);
	cr_assert_eq(top->position_y, 60);
}

Test(flis_canvas, rotate_180_keeps_dims) {
	flis_test_add_layer(flis_test_make_mono_fits(40, 40, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "p");
	top->position_x = 20; top->position_y = 20;  /* centre (25, 25) */
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 80;

	cr_assert_eq(flis_canvas_rotate(180.0), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 80u);
	/* (25,25) flipped through (50, 40) → (75, 55); pos = (70, 50). */
	cr_assert_eq(top->position_x, 70);
	cr_assert_eq(top->position_y, 50);
}

Test(flis_canvas, mirrorx_flips_layer_y) {
	flis_test_add_layer(flis_test_make_mono_fits(40, 40, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "p");
	top->position_x = 20; top->position_y = 30;
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	cr_assert_eq(flis_canvas_mirrorx(), 0);
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 100 - 30 - 10);  /* 60 */
}

Test(flis_canvas, mirrory_flips_layer_x) {
	flis_test_add_layer(flis_test_make_mono_fits(40, 40, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "p");
	top->position_x = 20; top->position_y = 30;
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	cr_assert_eq(flis_canvas_mirrory(), 0);
	cr_assert_eq(top->position_x, 100 - 20 - 10);  /* 70 */
	cr_assert_eq(top->position_y, 30);
}

/* No-FLIS sanity: every helper must refuse cleanly. */
Test(flis_canvas, helpers_refuse_when_no_flis) {
	cr_assert_neq(flis_canvas_resize(100, 100, 0, 0), 0);
	cr_assert_neq(flis_canvas_fit_to_layers(FALSE), 0);
	cr_assert_neq(flis_canvas_rotate(45.0), 0);
	cr_assert_neq(flis_canvas_mirrorx(), 0);
	cr_assert_neq(flis_canvas_mirrory(), 0);
}

/* Promote initialises canvas fields from the promoted fits. */
Test(flis_canvas, promote_sets_canvas_from_fits) {
	fits *f = flis_test_make_mono_fits(123, 456, 0.5f);
	gfit = f;
	com.uniq->fit = f;
	com.uniq->chans = 1;

	cr_assert_eq(flis_promote_from_gfit("base"), 0);
	cr_assert_eq(com.uniq->canvas_w, 123u);
	cr_assert_eq(com.uniq->canvas_h, 456u);
	cr_assert_eq(com.uniq->canvas_bg_r, 0.0);
	cr_assert_eq(com.uniq->canvas_bg_g, 0.0);
	cr_assert_eq(com.uniq->canvas_bg_b, 0.0);
	/* flis_canvas_rx/ry now read from com.uniq. */
	cr_assert_eq(flis_canvas_rx(), 123u);
	cr_assert_eq(flis_canvas_ry(), 456u);
}
