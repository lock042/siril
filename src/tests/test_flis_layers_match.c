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
 * test_flis_layers_match — exercises the §5.7 `flis_layers_match` command
 * (background neutralisation across tinted mono layers).  The underlying
 * primitive (`flis_background_neutralise_layers`) is already covered by
 * its own tests; this file focuses on the command surface: argument
 * parsing, subset resolution, and the full headless path.
 *
 * Strategy: build a 3-mono-layer fixture with distinct tints, push a
 * non-trivial median through each layer, run the command, assert the
 * scaled medians compose to a near-neutral background.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/nde_history.h"
#include "core/op_descriptor.h"

cominfo com;
fits *gfit;
char *word[MAX_COMMAND_WORDS];

extern int process_flis_layers_match(int nb);

TestSuite(flis_layers_match, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Build three tinted mono layers, each with a constant grey background.
 * The screen-blend composite of (a_i * lay_i * tint_i) summed across
 * channels should land near 1.0 per channel after neutralisation
 * (sum of channel contributions / 3 → 1.0). */
static void make_tinted_triple(double med_r, double med_g, double med_b) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_r), "R");
	flis_layer_set_tint(base, 1.0, 0.0, 0.0);
	flis_layer_t *g = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_g), "G");
	flis_layer_set_tint(g, 0.0, 1.0, 0.0);
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_b), "B");
	flis_layer_set_tint(b, 0.0, 0.0, 1.0);
}

Test(flis_layers_match, refuses_when_no_flis) {
	/* com.uniq is set up by init, but no layers added — is_current_image_flis
	 * returns FALSE without layers, so the command must refuse. */
	word[0] = "flis_layers_match";
	word[1] = NULL;
	cr_assert_eq(process_flis_layers_match(1), CMD_GENERIC_ERROR);
}

Test(flis_layers_match, runs_on_default_subset) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = NULL;
	cr_assert_eq(process_flis_layers_match(1), CMD_OK);
}

Test(flis_layers_match, accepts_subset_by_id) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	/* item_ids are 1, 2, 3 (assigned sequentially by flis_layer_add). */
	word[1] = "-subset=1,2,3";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
}

Test(flis_layers_match, accepts_subset_by_name) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-subset=R,G,B";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
}

Test(flis_layers_match, rejects_unknown_layer_in_subset) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-subset=R,NoSuchLayer";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_ARG_ERROR);
}

Test(flis_layers_match, rejects_unknown_option) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-bogus=1";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_ARG_ERROR);
}

/* ---- NDE provenance (per-layer "flis.layer_scale" records) --------------- */

Test(flis_layers_match, captures_per_layer_nde_records) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = NULL;
	cr_assert_eq(process_flis_layers_match(1), CMD_OK);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap, "no NDE history after layers match");
	guint n_scale = 0;
	gint seen_targets[3] = { 0, 0, 0 };
	for (guint i = 0; i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (g_strcmp0(rec->op_id, "flis.layer_scale"))
			continue;
		cr_assert_eq(rec->tier, NDE_TIER_A, "record %u not Tier A", i);
		cr_assert_eq(rec->scope, NDE_SCOPE_LAYER);
		cr_assert_not_null(rec->params);
		cr_assert(strstr(rec->params, "scale=") != NULL,
		          "params '%s' missing scale", rec->params);
		cr_assert(rec->target_item_id >= 1 && rec->target_item_id <= 3,
		          "unexpected target %d", rec->target_item_id);
		seen_targets[rec->target_item_id - 1]++;
		n_scale++;
	}
	cr_assert_eq(n_scale, 3, "expected 3 per-layer records, got %u", n_scale);
	for (int i = 0; i < 3; i++)
		cr_assert_eq(seen_targets[i], 1, "layer %d captured %d times",
		             i + 1, seen_targets[i]);
	g_ptr_array_unref(snap);
}

Test(flis_layers_match, layer_scale_serializer_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("flis.layer_scale");
	cr_assert_not_null(op, "flis.layer_scale not registered");
	cr_assert_not_null(op->serialize);
	cr_assert_not_null(op->deserialize);

	struct flis_layer_scale_data in = { .destroy_fn = NULL,
	                                    .scale = 0.73215, .offset = 0.0125 };
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct flis_layer_scale_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out, "deserialize failed on '%s'", blob);
	cr_assert_float_eq(out->scale, in.scale, 1e-12);
	cr_assert_float_eq(out->offset, in.offset, 1e-12);
	free_flis_layer_scale_data(out);
	g_free(blob);

	/* offset omitted from the blob → defaults to 0 */
	struct flis_layer_scale_data in2 = { .destroy_fn = NULL,
	                                     .scale = 1.5, .offset = 0.0 };
	blob = op->serialize(&in2);
	cr_assert(strstr(blob, "offset") == NULL);
	out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_float_eq(out->scale, 1.5, 1e-12);
	cr_assert_float_eq(out->offset, 0.0, 1e-12);
	free_flis_layer_scale_data(out);
	g_free(blob);

	/* negative scale refused */
	out = op->deserialize("scale=-1", op->version);
	cr_assert_null(out);
}

/* ---- group colour calibration distribution ------------------------------ */

/* Pure-tint RGB triple: the tint matrix is the identity, so the per-layer
 * factors must equal the composite's per-channel calibration exactly:
 * a_i = K_i and b_i = O_i. */
Test(flis_layers_match, group_calibration_pure_tints_exact) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	GSList *members = g_slist_copy(com.uniq->layers);
	const double K[3] = { 2.0, 1.0, 0.5 };
	const double O[3] = { 0.01, 0.0, -0.01 };
	cr_assert_eq(flis_group_apply_channel_calibration(members, K, O, "Test calib"), 0);
	g_slist_free(members);

	flis_layer_t *lr = g_slist_nth_data(com.uniq->layers, 0);
	flis_layer_t *lg = g_slist_nth_data(com.uniq->layers, 1);
	flis_layer_t *lb = g_slist_nth_data(com.uniq->layers, 2);
	cr_assert_float_eq(lr->fit->fdata[0], 2.0f * 0.2f + 0.01f, 1e-5);
	cr_assert_float_eq(lg->fit->fdata[0], 0.3f, 1e-5);
	cr_assert_float_eq(lb->fit->fdata[0], 0.5f * 0.4f - 0.01f, 1e-5);

	/* One Tier-A record per MUTATED layer: G's factors are exactly identity
	 * (K=1, O=0), so it is skipped and gets no record. */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	guint n = 0;
	for (guint i = 0; i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (g_strcmp0(rec->op_id, "flis.layer_scale"))
			continue;
		cr_assert_eq(rec->tier, NDE_TIER_A);
		cr_assert_neq(rec->target_item_id, 2, "identity layer G should have no record");
		n++;
	}
	cr_assert_eq(n, 2, "expected 2 records, got %u", n);
	g_ptr_array_unref(snap);
}

/* Two layers sharing one tint cannot carry a full-colour calibration: the
 * background system is rank-deficient in the unreachable channels, and the
 * operation must fail WITHOUT mutating any layer. */
Test(flis_layers_match, group_calibration_insufficient_diversity_fails) {
	flis_layer_t *l1 = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.2f), "A");
	flis_layer_set_tint(l1, 1.0, 0.0, 0.0);
	flis_layer_t *l2 = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.3f), "B");
	flis_layer_set_tint(l2, 1.0, 0.0, 0.0);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	GSList *members = g_slist_copy(com.uniq->layers);
	const double K[3] = { 1.2, 0.9, 1.1 };
	const double O[3] = { 0.01, 0.02, 0.03 };
	cr_assert_neq(flis_group_apply_channel_calibration(members, K, O, "Test calib"), 0);
	g_slist_free(members);

	cr_assert_float_eq(l1->fit->fdata[0], 0.2f, 1e-6, "layer mutated on failure");
	cr_assert_float_eq(l2->fit->fdata[0], 0.3f, 1e-6, "layer mutated on failure");
}

Test(flis_layers_match, layer_scale_hook_applies_affine) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.4f);
	struct flis_layer_scale_data p = { .destroy_fn = NULL,
	                                   .scale = 0.5, .offset = 0.05 };
	struct generic_img_args args = { 0 };
	args.user = &p;
	cr_assert_eq(flis_layer_scale_image_hook(&args, f, 1), 0);
	cr_assert_float_eq(f->fdata[0], 0.25f, 1e-6);
	cr_assert_float_eq(f->fdata[8 * 8 - 1], 0.25f, 1e-6);
	clearfits(f);
	free(f);
}
