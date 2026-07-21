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

/* ---- USHORT sparse-layer compose (C5 regression) -------------------- */

/* Build a constant-value USHORT fits, mono or RGB. */
static fits *make_ushort_fits(int rx, int ry, int chans, WORD value) {
	fits *f = NULL;
	if (new_fit_image(&f, rx, ry, chans, DATA_USHORT)) return NULL;
	size_t n = (size_t)rx * ry * chans;
	for (size_t i = 0; i < n; i++) f->data[i] = value;
	return f;
}

/* A USHORT layer larger than the canvas: the 16-bit→float conversion
 * buffer must be sized and strided by the LAYER pixel count, not the
 * canvas's.  Before the fix the blend loop read past the conversion
 * buffer for layer-local indices beyond canvas_w*canvas_h, compositing
 * heap garbage. */
Test(flis_sparse, ushort_layer_larger_than_canvas_composes_exactly) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    make_ushort_fits(12, 12, 1, USHRT_MAX), "big16");
	top->position_x = -2;   /* layer covers the whole canvas */
	top->position_y = -2;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);
	/* Every canvas pixel is covered by the (white) top layer; layer-local
	 * indices reach 9*12+9 = 117, past the canvas pixel count of 64. */
	for (guint y = 0; y < 8; y++)
		for (guint x = 0; x < 8; x++)
			cr_assert_float_eq(out->fpdata[0][y * 8 + x], 1.0f, 1e-4,
			                   "pixel (%u,%u) not composed from the layer buffer", x, y);
	free(out->fdata); free(out);
}

/* RGB variant: the per-plane stride of the conversion buffer must be the
 * layer pixel count.  Before the fix pre_g/pre_b pointed at canvas-count
 * offsets, so high layer-local indices read the wrong plane. */
Test(flis_sparse, ushort_rgb_layer_larger_than_canvas_composes_exactly) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f), "base");
	fits *tf = NULL;
	cr_assert_eq(new_fit_image(&tf, 12, 12, 3, DATA_USHORT), 0);
	size_t n = (size_t)12 * 12;
	for (size_t i = 0; i < n; i++) {
		tf->pdata[0][i] = USHRT_MAX;        /* R = 1.0 */
		tf->pdata[1][i] = USHRT_MAX / 2;    /* G ≈ 0.5 */
		tf->pdata[2][i] = 0;                /* B = 0.0 */
	}
	flis_layer_t *top = flis_test_add_layer(tf, "rgb16");
	top->position_x = -2;
	top->position_y = -2;
	com.uniq->canvas_w = 8;
	com.uniq->canvas_h = 8;

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);
	for (guint y = 0; y < 8; y++)
		for (guint x = 0; x < 8; x++) {
			size_t i = y * 8 + x;
			cr_assert_float_eq(out->fpdata[0][i], 1.0f, 1e-4);
			cr_assert_float_eq(out->fpdata[1][i], 0.5f, 1e-3,
			                   "G plane misread at (%u,%u) — conversion stride wrong", x, y);
			cr_assert_float_eq(out->fpdata[2][i], 0.0f, 1e-4);
		}
	free(out->fdata); free(out);
}

/* ---- displayed-mask dimension gate (C6 regression) ------------------ */

/* flis_get_displayed_mask hands its buffer to callers that index it with
 * gfit-linear indices and mask_t carries no dimensions: it must refuse to
 * return a layer mask whose dims differ from the displayed image (e.g.
 * a sparse layer's lmask while the canvas-sized composite is swapped
 * into gfit). */
Test(flis_sparse, displayed_mask_refuses_dim_mismatch) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 1.0f), "small");
	layermask_t *lm = flis_test_make_const_lmask(4, 4, 8, 1.0);
	cr_assert_eq(flis_layer_set_lmask(top, lm), 0);
	top->lmask_active = TRUE;
	com.uniq->flis_mask_view = 1;   /* LAYER mask view */
	uniq_set_active_layer(com.uniq, 1);   /* top active; gfit = 4x4 layer fit */

	mask_t out = { 0 };
	cr_assert(flis_get_displayed_mask(&out),
	          "matching dims must return the layer mask");
	cr_assert_eq(out.data, top->lmask->data);

	/* Simulate display time: the canvas-sized composite is swapped into
	 * gfit while the active layer (and its mask) stays 4x4. */
	fits *canvas = flis_test_make_rgb_fits(8, 8, 0.0f, 0.0f, 0.0f);
	fits *saved = gfit;
	gfit = canvas;
	mask_t out2 = { 0 };
	gboolean got = flis_get_displayed_mask(&out2);
	cr_assert(!got || out2.data != top->lmask->data,
	          "layer-sized mask must not be returned against a canvas-sized gfit");
	gfit = saved;
	clearfits(canvas);
	free(canvas);
}

/* ---- coordinate-model helpers (display ↔ active-layer-local) -------- */

Test(flis_sparse, display_to_layer_pt_handles_offset_and_bounds) {
	flis_test_add_layer(flis_test_make_rgb_fits(16, 16, 0.f, 0.f, 0.f), "base");
	flis_layer_t *sp = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 1.f), "sparse");
	sp->position_x = 6;
	sp->position_y = 8;
	com.uniq->canvas_w = 16;
	com.uniq->canvas_h = 16;
	uniq_set_active_layer(com.uniq, 1);   /* sparse active, gfit = 4x4 */

	pointi out;
	/* On-layer point translates to layer-local. */
	cr_assert(flis_display_to_active_layer_pt((pointi){ 7, 9 }, &out));
	cr_assert_eq(out.x, 1);
	cr_assert_eq(out.y, 1);
	/* Canvas point outside the layer: FALSE, but translation still
	 * reported. */
	cr_assert_not(flis_display_to_active_layer_pt((pointi){ 0, 0 }, &out));
	cr_assert_eq(out.x, -6);
	cr_assert_eq(out.y, -8);
	/* Just past the far corner. */
	cr_assert_not(flis_display_to_active_layer_pt((pointi){ 10, 12 }, NULL));
	cr_assert(flis_display_to_active_layer_pt((pointi){ 9, 11 }, NULL));

	/* Base (non-offset) layer behaves as identity. */
	uniq_set_active_layer(com.uniq, 0);
	cr_assert(flis_display_to_active_layer_pt((pointi){ 15, 15 }, &out));
	cr_assert_eq(out.x, 15);
	cr_assert_eq(out.y, 15);
	cr_assert_not(flis_display_to_active_layer_pt((pointi){ 16, 0 }, NULL));
}

Test(flis_sparse, display_to_layer_rect_partial_clips_full_requires_containment) {
	flis_test_add_layer(flis_test_make_rgb_fits(16, 16, 0.f, 0.f, 0.f), "base");
	flis_layer_t *sp = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 1.f), "sparse");
	sp->position_x = 6;
	sp->position_y = 8;
	com.uniq->canvas_w = 16;
	com.uniq->canvas_h = 16;
	uniq_set_active_layer(com.uniq, 1);

	rectangle out;
	/* Selection partially overlapping the sparse layer: partial mode
	 * clips to the intersection (the user-visible ROI-on-sparse-layer
	 * behaviour); full mode refuses. */
	rectangle sel = { 4, 6, 4, 4 };   /* overlaps layer rows/cols 0..1 */
	cr_assert(flis_display_to_active_layer_rect(&sel, &out, TRUE));
	cr_assert_eq(out.x, 0);
	cr_assert_eq(out.y, 0);
	cr_assert_eq(out.w, 2);
	cr_assert_eq(out.h, 2);
	cr_assert_not(flis_display_to_active_layer_rect(&sel, &out, FALSE));

	/* Fully inside. */
	rectangle sel2 = { 7, 9, 2, 2 };
	cr_assert(flis_display_to_active_layer_rect(&sel2, &out, FALSE));
	cr_assert_eq(out.x, 1);
	cr_assert_eq(out.y, 1);
	cr_assert_eq(out.w, 2);
	cr_assert_eq(out.h, 2);

	/* No overlap at all. */
	rectangle sel3 = { 0, 0, 4, 4 };
	cr_assert_not(flis_display_to_active_layer_rect(&sel3, &out, TRUE));
}
