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
 * test_flis_register_layers — surface tests for the §5.6
 * `flis_register_layers` command and the underlying primitive in
 * src/registration/flis_register.c.
 *
 * Pixel-level equivalence (offset detection on real star fixtures) is
 * deferred to fixture-driven integration tests — the DFT cross-
 * correlation primitive itself is already covered by the existing
 * registration test suite.  This file focuses on:
 *
 *   • argument parsing (-ref=, -interp=, -noclamp);
 *   • refusal paths (no FLIS, < 2 layers, unknown options);
 *   • the primitive accepts a heterogeneous layer list and runs to
 *     completion against trivial pixel content.
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "registration/registration.h"
#include "registration/flis_register.h"

cominfo com;
fits *gfit;
char *word[MAX_COMMAND_WORDS];

extern int process_flis_register_layers(int nb);

TestSuite(flis_register_layers, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

Test(flis_register_layers, refuses_when_no_flis) {
	word[0] = "flis_register_layers";
	word[1] = NULL;
	cr_assert_eq(process_flis_register_layers(1), CMD_GENERIC_ERROR);
}

Test(flis_register_layers, refuses_single_layer) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "only");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_register_layers";
	word[1] = NULL;
	cr_assert_eq(process_flis_register_layers(1), CMD_GENERIC_ERROR);
}

Test(flis_register_layers, rejects_unknown_ref_layer) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "top");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_register_layers";
	word[1] = "-ref=NoSuchLayer";
	word[2] = NULL;
	cr_assert_eq(process_flis_register_layers(2), CMD_ARG_ERROR);
}

Test(flis_register_layers, rejects_unknown_interp) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "top");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_register_layers";
	word[1] = "-interp=fantasymode";
	word[2] = NULL;
	cr_assert_eq(process_flis_register_layers(2), CMD_ARG_ERROR);
}

/* "accepts_known_interp_names" was intentionally omitted: a successful
 * -interp= parse leads into register_shift_dft / register_apply_reg which
 * are not designed to run on the constant 16x16 mono fixtures available
 * to a unit test (no DFT peak; the pipeline may dereference NULL stats
 * or hit an internal assert).  Coverage for the parsing surface comes
 * from the rejection tests above — a successful match path is exercised
 * by integration tests against real star fixtures. */

Test(flis_register_layers, rejects_unknown_option) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "top");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_register_layers";
	word[1] = "-bogus=1";
	word[2] = NULL;
	cr_assert_eq(process_flis_register_layers(2), CMD_ARG_ERROR);
}

/* =====================================================================
 * End-to-end alignment regression test.
 *
 * Guards against the flis-gtk4 double-apply regression: the GLOBAL
 * method used to resolve to single-pass register_star_alignment, which
 * resamples the layer pixels inline; the register_apply_reg() pass in
 * flis_register_layers then warped them a second time and every layer
 * came out misplaced (in both axes once rotation was involved).
 *
 * Uses the DFT shift method on a synthetic gaussian star field: fully
 * deterministic, no star matching.  The assertion is content-anchored
 * (composite peak positions), so it holds regardless of which position
 * convention flis_register_layers stores internally.
 * ===================================================================== */

/* Deterministic star field: gaussian spots from an LCG, content offset
 * by (sx, sy) in FITS buffer coordinates. */
static fits *make_star_field(int w, int h, double sx, double sy) {
	fits *f = flis_test_make_mono_fits(w, h, 0.02f);
	if (!f) return NULL;
	guint32 lcg = 42;
	for (int s = 0; s < 40; s++) {
		lcg = lcg * 1103515245u + 12345u;
		double cx = 30.0 + (lcg >> 8) % (w - 60) + sx;
		lcg = lcg * 1103515245u + 12345u;
		double cy = 30.0 + (lcg >> 8) % (h - 60) + sy;
		lcg = lcg * 1103515245u + 12345u;
		double amp = 0.3 + 0.6 * ((lcg >> 8) % 1000) / 1000.0;
		for (int y = (int)cy - 6; y <= (int)cy + 6; y++) {
			if (y < 0 || y >= h) continue;
			for (int x = (int)cx - 6; x <= (int)cx + 6; x++) {
				if (x < 0 || x >= w) continue;
				double d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
				float v = f->fdata[(size_t)y * w + x]
				          + (float)(amp * exp(-d2 / (2.0 * 2.0 * 2.0)));
				f->fdata[(size_t)y * w + x] = v > 1.f ? 1.f : v;
			}
		}
	}
	return f;
}

/* Composite with only @keep visible; return the buffer index of the
 * brightest pixel of channel 0, or -1. */
static gssize solo_composite_peak(flis_layer_t *keep, guint *out_w) {
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (lay) lay->visible = (lay == keep);
	}
	fits *comp = flis_render_layers(com.uniq->layers);
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (lay) lay->visible = TRUE;
	}
	if (!comp) return -1;
	size_t n = (size_t)comp->rx * comp->ry;
	size_t best = 0;
	for (size_t i = 1; i < n; i++)
		if (comp->fpdata[0][i] > comp->fpdata[0][best]) best = i;
	if (out_w) *out_w = comp->rx;
	gssize rv = (gssize)best;
	clearfits(comp);
	free(comp);
	return rv;
}

Test(flis_register_layers, dft_aligns_content_shifted_layer) {
	const int W = 400, H = 300;
	flis_layer_t *base = flis_test_add_layer(make_star_field(W, H, 0, 0), "base");
	flis_layer_t *top  = flis_test_add_layer(make_star_field(W, H, 15, 9), "shifted");
	cr_assert_not_null(base);
	cr_assert_not_null(top);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	/* The registration pipeline consults the memory limits. */
	com.pref.mem_mode = RATIO;
	com.pref.memory_ratio = 0.9;

	/* DFT phase correlation on a square selection. */
	com.selection = (rectangle){ 100, 50, 128, 128 };
	selection_type sel;
	transformation_type tx;
	registration_function method =
	    flis_register_resolve_method(FLIS_REG_DFT, &sel, &tx);
	cr_assert_not_null(method);

	int rv = flis_register_layers(NULL, NULL, method, sel, tx,
	                              OPENCV_LANCZOS4, TRUE);
	cr_assert_eq(rv, 0, "flis_register_layers failed (%d)", rv);

	/* Content-anchored check: the brightest star must land on the same
	 * canvas pixel whichever layer is composited. */
	guint cw1 = 0, cw2 = 0;
	gssize p_base = solo_composite_peak(base, &cw1);
	gssize p_top  = solo_composite_peak(top, &cw2);
	cr_assert(p_base >= 0 && p_top >= 0, "composite failed");
	cr_assert_eq(cw1, cw2);
	gint bx = (gint)(p_base % cw1), by = (gint)(p_base / cw1);
	gint tx_ = (gint)(p_top % cw2),  ty = (gint)(p_top / cw2);
	cr_assert(ABS(bx - tx_) <= 1 && ABS(by - ty) <= 1,
	          "layers misaligned: base peak (%d,%d) vs shifted peak (%d,%d)",
	          bx, by, tx_, ty);

	/* Documented position convention for a pure (+15,+9) buffer-space
	 * content shift; if these change, the compose convention changed
	 * with them and the content check above is what matters. */
	cr_assert_eq(top->position_x, -15, "position_x = %d", top->position_x);
	cr_assert_eq(top->position_y, 9, "position_y = %d", top->position_y);
}
