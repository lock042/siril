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
 * test_flis_sparse — §6.1 sparse-layer correctness audit.
 *
 * Layers may extend past the canvas in any direction (FLIS spec
 * §6.2 / §12.2 + §7 canvas decoupling).  These tests exercise the
 * key code paths (compositor, coordinate transform, canvas helpers)
 * with deliberately out-of-canvas positions to confirm everything
 * clips rather than crashes / produces garbage.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"
#include "io/flis_compose.h"

cominfo com;
fits *gfit;

TestSuite(flis_sparse, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Compositing a layer that pokes off the left edge clips correctly:
 * the in-canvas portion appears at the expected place; out-of-canvas
 * pixels never trigger memory access. */
Test(flis_sparse, compose_layer_off_left_edge) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 1.0f, 1.0f, 1.0f), "patch");
	/* Top layer at x=-2: rightmost 2 columns of the layer overlap the
	 * canvas's first 2 columns. */
	top->position_x = -2;
	top->position_y = 0;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out, "compose returned NULL on sparse layer");
	cr_assert_eq(out->rx, 8u);
	cr_assert_eq(out->ry, 8u);
	/* FITS storage is bottom-up: row index N corresponds to display row
	 * (canvas_h - 1 - N).  Top layer at position_y=0 covers display rows
	 * [0, 4) = FITS rows 4..7.  Columns: position_x=-2 + lW=4 → covers
	 * canvas columns 0, 1 only.  Sample the top-left of the canvas:
	 * display (0,0) = FITS (7,0) = index 56. */
	cr_assert_float_eq(out->fpdata[0][7 * 8 + 0], 1.0f, 1e-5);
	cr_assert_float_eq(out->fpdata[0][7 * 8 + 1], 1.0f, 1e-5);
	cr_assert_float_eq(out->fpdata[0][7 * 8 + 2], 0.0f, 1e-5);
	cr_assert_float_eq(out->fpdata[0][7 * 8 + 7], 0.0f, 1e-5);
	/* Bottom row of canvas (FITS row 0) is base only — top layer doesn't
	 * reach down there. */
	cr_assert_float_eq(out->fpdata[0][0 * 8 + 0], 0.0f, 1e-5);
	free(out->fdata); free(out);
}

/* Layer extending past the right + bottom: the portion outside the
 * canvas is discarded; the in-canvas portion blends normally. */
Test(flis_sparse, compose_layer_off_right_and_bottom) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(6, 6, 1.0f, 1.0f, 1.0f), "patch");
	top->position_x = 5;  /* extends to x=11 — clipped at x=8 */
	top->position_y = 5;  /* extends to y=11 — clipped at y=8 */
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);
	/* Top layer covers display rows [5, 8) (= FITS rows 0..2) and
	 * columns [5, 8).  Sample bottom-right of canvas (display (7,7) =
	 * FITS (0,7)): covered by top → white. */
	cr_assert_float_eq(out->fpdata[0][0 * 8 + 7], 1.0f, 1e-5);
	/* Top-left of canvas (display (0,0) = FITS (7,0)): base only. */
	cr_assert_float_eq(out->fpdata[0][7 * 8 + 0], 0.0f, 1e-5);
	/* Just outside the in-canvas portion of top: display (4,4) = FITS
	 * (3, 4); top doesn't reach here. */
	cr_assert_float_eq(out->fpdata[0][3 * 8 + 4], 0.0f, 1e-5);
	free(out->fdata); free(out);
}

/* A layer entirely outside the canvas leaves the composite untouched. */
Test(flis_sparse, compose_layer_entirely_off_canvas) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.3f, 0.3f, 0.3f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 1.0f, 1.0f, 1.0f), "way_outside");
	top->position_x = -100;
	top->position_y = -100;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);
	/* Every canvas pixel is the base colour (0.3). */
	for (size_t i = 0; i < 64; i++) {
		cr_assert_float_eq(out->fpdata[0][i], 0.3f, 1e-5);
	}
	free(out->fdata); free(out);
}

/* canvas_to_pixel_index returns FALSE for clicks outside the active
 * layer (including outside-canvas clicks when the active layer is
 * sparse). */
Test(flis_sparse, pixel_index_returns_false_outside_active_layer) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.5f), "patch");
	top->position_x = 0;
	top->position_y = 0;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();

	size_t idx;
	/* Click inside top: returns TRUE */
	cr_assert(flis_canvas_to_pixel_index(2, 2, 8, &idx));
	/* Click outside top but inside canvas: returns FALSE (top is sparse) */
	cr_assert_not(flis_canvas_to_pixel_index(6, 6, 8, &idx));
}

/* Canvas resize accepting layers that extend past current canvas dims
 * — flis_canvas_fit_to_layers grows the canvas to fit and shifts every
 * layer accordingly. */
Test(flis_sparse, fit_to_layers_grows_canvas_to_include_sparse) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.5f), "patch");
	top->position_x = -3;     /* extends past left edge */
	top->position_y = 9;      /* extends past bottom edge */
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	cr_assert_eq(flis_canvas_fit_to_layers(FALSE), 0);
	/* bbox: x [-3, 8) → width 11.  y [0, 13) → height 13.
	 * After shift (+3, 0): base (3, 0), top (0, 9). */
	cr_assert_eq(com.uniq->canvas_w, 11u);
	cr_assert_eq(com.uniq->canvas_h, 13u);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 3);
	cr_assert_eq(base->position_y, 0);
	cr_assert_eq(top->position_x, 0);
	cr_assert_eq(top->position_y, 9);
}

/* Canvas mirror with a sparse layer: position flips correctly even
 * for negative starting positions. */
Test(flis_sparse, mirrorx_flips_sparse_layer_position) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.5f), "patch");
	top->position_x = -2;
	top->position_y = -1;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	cr_assert_eq(flis_canvas_mirrorx(), 0);
	/* mirrorx flips y: new_y = ch - y - layer_h = 8 - (-1) - 4 = 5. */
	cr_assert_eq(top->position_y, 5);
	cr_assert_eq(top->position_x, -2);  /* X unchanged */
}

/* Canvas rotation handles sparse layers — bbox of layer centre is
 * rotated cleanly regardless of starting position. */
Test(flis_sparse, rotate_handles_sparse_positions) {
	flis_test_add_layer(flis_test_make_mono_fits(40, 40, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "p");
	top->position_x = -20;   /* off-canvas top */
	top->position_y = -10;
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;

	/* 180° rotation just flips through canvas centre. */
	cr_assert_eq(flis_canvas_rotate(180.0), 0);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 100u);
	/* old centre (-15, -5) → new centre (115, 105); pos (110, 100). */
	cr_assert_eq(top->position_x, 110);
	cr_assert_eq(top->position_y, 100);
}
