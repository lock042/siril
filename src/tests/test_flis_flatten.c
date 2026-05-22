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
 * test_flis_flatten — black-box tests of the flatten and merge-down
 * pipelines.  These exercise the same compositing arithmetic as
 * test_flis_compose, but through the user-visible flis_flatten_all
 * and flis_merge_down_layer entry points — so they catch bugs in the
 * flatten/merge plumbing (layer deletion, base preservation, mask
 * consumption, undo purge) on top of the kernel.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_flatten, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* W3C reference formulas — duplicated from test_flis_compose so the
 * flatten tests are self-contained. */
static float ref_normal(float s, float d)    { (void)d; return s; }
static float ref_multiply(float s, float d)  { return s * d; }
static float ref_screen(float s, float d)    { return s + d - s * d; }
static float ref_darken(float s, float d)    { return fminf(s, d); }
static float ref_lighten(float s, float d)   { return fmaxf(s, d); }
static float ref_difference(float s, float d){ return fabsf(d - s); }
static float ref_exclusion(float s, float d) { return d + s - 2.f * d * s; }

/* Build a 2-layer FLIS with constant base + top, run flis_flatten_all,
 * and return the (now-1-layer) com.uniq's base layer's fits.  Caller
 * does NOT free — flis_test_cleanup_com handles it. */
static fits *do_flatten(float br, float bg, float bb,
                        float tr, float tg, float tb,
                        flis_blend_mode_t mode, float opacity) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, br, bg, bb), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, tr, tg, tb), "top");
	flis_layer_set_blend_mode(top, mode);
	flis_layer_set_opacity(top, opacity);
	int rv = flis_flatten_all();
	cr_assert_eq(rv, 0, "flis_flatten_all returned %d", rv);
	cr_assert_eq(flis_layer_count(), 1, "flatten should leave exactly 1 layer");
	return ((flis_layer_t *)com.uniq->layers->data)->fit;
}

/* Per-blend-mode flatten test: pixels match the reference formula. */
#define FLATTEN_BLEND(name, mode_enum, ref_fn)                                 \
Test(flis_flatten, mode_##name) {                                              \
	float br = 0.2f, bg = 0.4f, bb = 0.6f;                                     \
	float tr = 0.7f, tg = 0.3f, tb = 0.5f;                                     \
	fits *fit = do_flatten(br, bg, bb, tr, tg, tb, mode_enum, 1.0f);           \
	cr_assert_float_eq(fit->fpdata[0][0], ref_fn(tr, br), 1e-4);               \
	cr_assert_float_eq(fit->fpdata[1][0], ref_fn(tg, bg), 1e-4);               \
	cr_assert_float_eq(fit->fpdata[2][0], ref_fn(tb, bb), 1e-4);               \
}

FLATTEN_BLEND(normal,     FLIS_BLEND_NORMAL,     ref_normal)
FLATTEN_BLEND(multiply,   FLIS_BLEND_MULTIPLY,   ref_multiply)
FLATTEN_BLEND(screen,     FLIS_BLEND_SCREEN,     ref_screen)
FLATTEN_BLEND(darken,     FLIS_BLEND_DARKEN,     ref_darken)
FLATTEN_BLEND(lighten,    FLIS_BLEND_LIGHTEN,    ref_lighten)
FLATTEN_BLEND(difference, FLIS_BLEND_DIFFERENCE, ref_difference)
FLATTEN_BLEND(exclusion,  FLIS_BLEND_EXCLUSION,  ref_exclusion)

/* Opacity-scan flatten: result is a Porter-Duff interpolation between
 * base and full-blend at intermediate opacities. */
Test(flis_flatten, multiply_opacity_scan) {
	for (float opa = 0.f; opa <= 1.0001f; opa += 0.25f) {
		fits *fit = do_flatten(0.3f, 0.3f, 0.3f, 0.9f, 0.9f, 0.9f,
		                       FLIS_BLEND_MULTIPLY, opa);
		float blended = ref_multiply(0.9f, 0.3f);  /* = 0.27 */
		float expected = blended * opa + 0.3f * (1.f - opa);
		cr_assert_float_eq(fit->fpdata[0][0], expected, 1e-4,
		                   "opacity %g: got %g expected %g",
		                   opa, fit->fpdata[0][0], expected);
		flis_free_layers(com.uniq);
		com.uniq->next_item_id = 1;
	}
}

/* Invisible top is skipped — base preserved unchanged. */
Test(flis_flatten, invisible_layer_skipped) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.4f, 0.5f, 0.6f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 1.f, 1.f, 1.f), "top");
	flis_layer_set_visible(top, FALSE);
	cr_assert_eq(flis_flatten_all(), 0);
	cr_assert_eq(flis_layer_count(), 1);
	fits *fit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	cr_assert_float_eq(fit->fpdata[0][0], 0.4f, 1e-5);
	cr_assert_float_eq(fit->fpdata[1][0], 0.5f, 1e-5);
	cr_assert_float_eq(fit->fpdata[2][0], 0.6f, 1e-5);
}

/* Tinted mono top in Screen mode over RGB base — exercises the
 * mono-broadcast-with-tint path through the flatten pipeline. */
Test(flis_flatten, tinted_mono_screen) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.1f, 0.1f, 0.1f), "base");
	flis_layer_t *ha = flis_test_add_layer(flis_test_make_mono_fits(2, 2, 0.5f), "Ha");
	flis_layer_set_blend_mode(ha, FLIS_BLEND_SCREEN);
	flis_layer_set_tint(ha, 1.0, 0.2, 0.1);

	cr_assert_eq(flis_flatten_all(), 0);
	fits *fit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	/* effective source: (0.5*1.0, 0.5*0.2, 0.5*0.1) = (0.5, 0.1, 0.05) */
	cr_assert_float_eq(fit->fpdata[0][0], ref_screen(0.5f,  0.1f), 1e-4);
	cr_assert_float_eq(fit->fpdata[1][0], ref_screen(0.1f,  0.1f), 1e-4);
	cr_assert_float_eq(fit->fpdata[2][0], ref_screen(0.05f, 0.1f), 1e-4);
}

/* Post-flatten state: exactly one layer remains, base has NORMAL blend /
 * full opacity / visible, masks cleared.  This is the "did flatten clean
 * up correctly" sanity check. */
Test(flis_flatten, post_state_single_layer_reset) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 1.0f, 0.0f, 0.0f), "red");
	flis_layer_set_blend_mode(top, FLIS_BLEND_SCREEN);
	flis_layer_set_opacity(top, 0.5f);

	cr_assert_eq(flis_flatten_all(), 0);
	cr_assert_eq(flis_layer_count(), 1);

	flis_layer_t *result = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(result->blend_mode, FLIS_BLEND_NORMAL,
	             "post-flatten base blend mode should be NORMAL");
	cr_assert_float_eq(result->opacity, 1.0f, 1e-6,
	                   "post-flatten base opacity should be 1.0");
	cr_assert(result->visible, "post-flatten base must be visible");
	cr_assert_null(result->lmask, "post-flatten base lmask cleared");
	cr_assert(!result->fit->mask || !result->fit->mask->data,
	          "post-flatten base pmask cleared");
}

/* Merge-down between two layers: behaves identically to a 2-layer
 * flatten for the colour math, and consumes the top layer. */
Test(flis_flatten, merge_down_per_blend_mode_normal) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.2f, 0.4f, 0.6f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.7f, 0.3f, 0.5f), "top");
	flis_layer_set_blend_mode(top, FLIS_BLEND_NORMAL);

	cr_assert_eq(flis_merge_down_layer(top), 0);
	cr_assert_eq(flis_layer_count(), 1);
	fits *fit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	cr_assert_float_eq(fit->fpdata[0][0], 0.7f, 1e-4);
	cr_assert_float_eq(fit->fpdata[1][0], 0.3f, 1e-4);
	cr_assert_float_eq(fit->fpdata[2][0], 0.5f, 1e-4);
}

Test(flis_flatten, merge_down_per_blend_mode_multiply) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.5f, 0.5f, 0.5f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.4f, 0.6f, 0.8f), "top");
	flis_layer_set_blend_mode(top, FLIS_BLEND_MULTIPLY);

	cr_assert_eq(flis_merge_down_layer(top), 0);
	cr_assert_eq(flis_layer_count(), 1);
	fits *fit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	cr_assert_float_eq(fit->fpdata[0][0], 0.5f * 0.4f, 1e-4);
	cr_assert_float_eq(fit->fpdata[1][0], 0.5f * 0.6f, 1e-4);
	cr_assert_float_eq(fit->fpdata[2][0], 0.5f * 0.8f, 1e-4);
}

/* Sparse layer flatten: top covers only a subrect of the canvas.  Pixels
 * inside the rect blend per mode; pixels outside are base unchanged.
 *
 * Note on coordinates: position_y is in FITS convention (origin
 * bottom-left), so the kernel maps a patch with position_y=2 + lH=2 on
 * an H=8 canvas into top-down canvas rows oy=H−py−lH=4, covering rows
 * [4, 6) — not rows [2, 4) as a naive reading would suggest. */
Test(flis_flatten, sparse_layer_outside_extent_preserved) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.1f, 0.1f, 0.1f), "canvas");
	flis_layer_t *patch = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 1.0f, 1.0f, 1.0f), "patch");
	patch->position_x = 2;
	patch->position_y = 2;   /* maps to top-down rows [4, 6) */
	flis_layer_set_blend_mode(patch, FLIS_BLEND_NORMAL);

	cr_assert_eq(flis_flatten_all(), 0);
	fits *fit = ((flis_layer_t *)com.uniq->layers->data)->fit;

	/* Row 0 col 0: outside patch — base */
	cr_assert_float_eq(fit->fpdata[0][0 * 8 + 0], 0.1f, 1e-5);
	/* Row 4 col 2: inside patch — fully replaced with top */
	cr_assert_float_eq(fit->fpdata[0][4 * 8 + 2], 1.0f, 1e-5);
	cr_assert_float_eq(fit->fpdata[0][5 * 8 + 3], 1.0f, 1e-5);
	/* Row 7 col 7: outside — base */
	cr_assert_float_eq(fit->fpdata[0][7 * 8 + 7], 0.1f, 1e-5);
	/* Row 3 col 3: just outside patch in y — base */
	cr_assert_float_eq(fit->fpdata[0][3 * 8 + 3], 0.1f, 1e-5);
}
