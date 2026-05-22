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
 * test_flis_compose — white-box tests of the compositing kernel
 * (flis_render_layers).  For each spec §5 blend mode, feed two
 * constant-colour layers and compare the output to a W3C reference
 * formula evaluated at the test colours.  This is the canonical
 * correctness oracle for the FLIS compositing arithmetic; the §1.2
 * flatten / merge-down tests exercise the same arithmetic through
 * the user-visible pipelines, and the §3.5 display test exercises
 * it through the GPU snapshot path — all three must agree.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_compose, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* W3C compositing reference formulas (Compositing and Blending Level 1,
 * <https://www.w3.org/TR/compositing-1/>).  Inputs are linear-light [0,1]. */
static float ref_normal(float s, float d)     { (void)d; return s; }
static float ref_multiply(float s, float d)   { return s * d; }
static float ref_screen(float s, float d)     { return s + d - s * d; }
static float ref_overlay(float s, float d)    {
	return (d <= 0.5f) ? 2.f * s * d : 1.f - 2.f * (1.f - s) * (1.f - d);
}
static float ref_darken(float s, float d)     { return fminf(s, d); }
static float ref_lighten(float s, float d)    { return fmaxf(s, d); }
static float ref_color_dodge(float s, float d) {
	if (d <= 0.f) return 0.f;
	if (s >= 1.f) return 1.f;
	return fminf(d / (1.f - s), 1.f);
}
static float ref_color_burn(float s, float d) {
	if (d >= 1.f) return 1.f;
	if (s <= 0.f) return 0.f;
	return 1.f - fminf((1.f - d) / s, 1.f);
}
static float ref_hard_light(float s, float d) {
	return (s <= 0.5f) ? 2.f * s * d : 1.f - 2.f * (1.f - s) * (1.f - d);
}
static float ref_soft_light(float s, float d) {
	float D = (d <= 0.25f) ? ((16.f * d - 12.f) * d + 4.f) * d : sqrtf(d);
	return (s <= 0.5f) ? d - (1.f - 2.f * s) * d * (1.f - d)
	                   : d + (2.f * s - 1.f) * (D - d);
}
static float ref_difference(float s, float d) { return fabsf(d - s); }
static float ref_exclusion(float s, float d)  { return d + s - 2.f * d * s; }

/* Helper: build a constant-colour FLIS with two layers and a chosen blend
 * mode + opacity on top, run the kernel, return the float RGB output.
 * Caller frees the returned fits.  Base and top are RGB (so per-channel
 * blend modes apply uniformly). */
static fits *run_compose_two(float br, float bg, float bb,
                             float tr, float tg, float tb,
                             flis_blend_mode_t mode, float opacity) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, br, bg, bb), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, tr, tg, tb), "top");
	flis_layer_set_blend_mode(top, mode);
	flis_layer_set_opacity(top, opacity);
	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *result = flis_render_layers(sorted);
	g_slist_free(sorted);
	return result;
}

/* Per-blend-mode parameterised test: one Test() per mode, asserting all
 * three channels of the (2x2) composite match the reference formula. */
#define BLEND_TEST(name, mode_enum, ref_fn)                                      \
Test(flis_compose, blend_##name) {                                               \
	float br = 0.2f, bg = 0.4f, bb = 0.6f;                                       \
	float tr = 0.7f, tg = 0.3f, tb = 0.5f;                                       \
	fits *out = run_compose_two(br, bg, bb, tr, tg, tb, mode_enum, 1.0f);        \
	cr_assert_not_null(out, "compose returned NULL");                            \
	cr_assert_eq(out->naxes[2], 3);                                              \
	float exp_r = ref_fn(tr, br);                                                \
	float exp_g = ref_fn(tg, bg);                                                \
	float exp_b = ref_fn(tb, bb);                                                \
	cr_assert_float_eq(out->fpdata[0][0], exp_r, 1e-4,                           \
	                   "R: got %g, expected %g", out->fpdata[0][0], exp_r);      \
	cr_assert_float_eq(out->fpdata[1][0], exp_g, 1e-4,                           \
	                   "G: got %g, expected %g", out->fpdata[1][0], exp_g);     \
	cr_assert_float_eq(out->fpdata[2][0], exp_b, 1e-4,                           \
	                   "B: got %g, expected %g", out->fpdata[2][0], exp_b);     \
	clearfits(out); free(out);                                                   \
}

BLEND_TEST(normal,      FLIS_BLEND_NORMAL,      ref_normal)
BLEND_TEST(multiply,    FLIS_BLEND_MULTIPLY,    ref_multiply)
BLEND_TEST(screen,      FLIS_BLEND_SCREEN,      ref_screen)
BLEND_TEST(overlay,     FLIS_BLEND_OVERLAY,     ref_overlay)
BLEND_TEST(darken,      FLIS_BLEND_DARKEN,      ref_darken)
BLEND_TEST(lighten,     FLIS_BLEND_LIGHTEN,     ref_lighten)
BLEND_TEST(color_dodge, FLIS_BLEND_COLOR_DODGE, ref_color_dodge)
BLEND_TEST(color_burn,  FLIS_BLEND_COLOR_BURN,  ref_color_burn)
BLEND_TEST(hard_light,  FLIS_BLEND_HARD_LIGHT,  ref_hard_light)
BLEND_TEST(soft_light,  FLIS_BLEND_SOFT_LIGHT,  ref_soft_light)
BLEND_TEST(difference,  FLIS_BLEND_DIFFERENCE,  ref_difference)
BLEND_TEST(exclusion,   FLIS_BLEND_EXCLUSION,   ref_exclusion)

/* Opacity test: at opacity=0 the composite equals the base; at opacity=1
 * it equals the per-mode blend formula; intermediate values are a linear
 * Porter-Duff interpolation between the two. */
Test(flis_compose, opacity_interpolates_blend_with_base) {
	float br = 0.3f, bg = 0.3f, bb = 0.3f;
	float tr = 0.9f, tg = 0.9f, tb = 0.9f;

	for (float opa = 0.f; opa <= 1.0001f; opa += 0.25f) {
		fits *out = run_compose_two(br, bg, bb, tr, tg, tb, FLIS_BLEND_MULTIPLY, opa);
		cr_assert_not_null(out);
		float blended = ref_multiply(tr, br);  /* same per channel: 0.27 */
		float expected = blended * opa + br * (1.f - opa);
		cr_assert_float_eq(out->fpdata[0][0], expected, 1e-4,
		                   "opacity=%g: got %g, expected %g",
		                   opa, out->fpdata[0][0], expected);
		clearfits(out); free(out);
		flis_free_layers(com.uniq);
	}
}

/* Invisible top layer should be skipped entirely, regardless of blend mode. */
Test(flis_compose, invisible_top_layer_skipped) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.4f, 0.4f, 0.4f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 1.f, 0.f, 0.f), "top");
	flis_layer_set_visible(top, FALSE);

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);
	cr_assert_float_eq(out->fpdata[0][0], 0.4f, 1e-5);
	cr_assert_float_eq(out->fpdata[1][0], 0.4f, 1e-5);
	cr_assert_float_eq(out->fpdata[2][0], 0.4f, 1e-5);
	clearfits(out); free(out);
}

/* Mono layer with no tint: pixel value broadcasts uniformly to all three
 * RGB channels (spec §6.5 default behaviour). */
Test(flis_compose, mono_top_no_tint_broadcasts_to_rgb) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.1f, 0.2f, 0.3f), "base");
	flis_layer_t *mono = flis_test_add_layer(flis_test_make_mono_fits(2, 2, 0.5f), "mono");
	flis_layer_set_blend_mode(mono, FLIS_BLEND_SCREEN);

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);

	/* screen(0.5, base[c]) per channel — same source 0.5 broadcast everywhere. */
	cr_assert_float_eq(out->fpdata[0][0], ref_screen(0.5f, 0.1f), 1e-4);
	cr_assert_float_eq(out->fpdata[1][0], ref_screen(0.5f, 0.2f), 1e-4);
	cr_assert_float_eq(out->fpdata[2][0], ref_screen(0.5f, 0.3f), 1e-4);
	clearfits(out); free(out);
}

/* Mono with LAYER_COLOR tint: pixel value is multiplied by the tint vector
 * before blending (spec §6.5).  Ha red tint (1.0, 0.2, 0.1) on a mono
 * value of 0.5 produces RGB (0.5, 0.1, 0.05) as the effective source. */
Test(flis_compose, mono_top_with_layer_color_tint) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *ha = flis_test_add_layer(flis_test_make_mono_fits(2, 2, 0.5f), "Ha");
	flis_layer_set_blend_mode(ha, FLIS_BLEND_NORMAL);
	flis_layer_set_tint(ha, 1.0, 0.2, 0.1);

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);

	cr_assert_float_eq(out->fpdata[0][0], 0.5f, 1e-5);
	cr_assert_float_eq(out->fpdata[1][0], 0.1f, 1e-5);
	cr_assert_float_eq(out->fpdata[2][0], 0.05f, 1e-5);
	clearfits(out); free(out);
}

/* Layer mask test: an all-white lmask should produce the same result as
 * no mask; an all-mid-grey lmask is equivalent to opacity 0.5. */
Test(flis_compose, lmask_constant_grey_equivalent_to_opacity) {
	flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 0.0f, 0.0f, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_rgb_fits(2, 2, 1.0f, 1.0f, 1.0f), "top");
	flis_layer_set_blend_mode(top, FLIS_BLEND_NORMAL);
	flis_layer_set_lmask(top, flis_test_make_const_lmask(2, 2, 8, 0.5));

	GSList *sorted = g_slist_copy(com.uniq->layers);
	fits *out = flis_render_layers(sorted);
	g_slist_free(sorted);
	cr_assert_not_null(out);

	/* effective alpha = opacity (1.0) * mask (~0.5) → result ≈ 0.5
	 * The 0.5 mask comes from value*255 rounded; 0.5*255=127.5→128, mask
	 * value back to float = 128/255 ≈ 0.50196.  Tolerance accounts for it. */
	cr_assert_float_eq(out->fpdata[0][0], 0.5f, 5e-3);
	clearfits(out); free(out);
}
