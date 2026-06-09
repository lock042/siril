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
 * test_flis_stress — §6.4 stress tests.
 *
 * Unit-test-scale stress: many small layers + every blend mode + every
 * group configuration.  We don't measure wall-time here (the stress
 * tests at production resolutions are run manually); the goal is to
 * confirm correctness at scale, not perf.  VRAM-pressure tests against
 * the GPU compose path are deferred — they need a GPU-capable harness.
 */

#include <criterion/criterion.h>
#include <glib/gstdio.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

static char *tmpdir = NULL;
static char *tmppath = NULL;

static void setup(void) {
	flis_test_init_com();
	tmpdir = g_dir_make_tmp("flis-stress-XXXXXX", NULL);
	tmppath = g_build_filename(tmpdir, "stress.flis", NULL);
}

static void teardown(void) {
	if (tmppath) { g_unlink(tmppath); g_free(tmppath); tmppath = NULL; }
	if (tmpdir)  { g_rmdir(tmpdir);   g_free(tmpdir);   tmpdir  = NULL; }
	flis_test_cleanup_com();
}

TestSuite(flis_stress, .init = setup, .fini = teardown);

/* 16-layer FLIS round-trip — verifies the meta table handles many
 * layers and the layer-order resort logic is stable at scale.
 * Resolution kept tiny so the test runs fast; the structural
 * complexity is the load-bearing part. */
Test(flis_stress, sixteen_layers_round_trip) {
	const int n = 16;
	for (int i = 0; i < n; i++) {
		float v = 0.05f * (float)i;
		gchar *name = g_strdup_printf("L%02d", i);
		flis_layer_t *lay = flis_test_add_layer(
		    flis_test_make_mono_fits(4, 4, v), name);
		g_free(name);
		cr_assert_not_null(lay);
		lay->opacity = 1.0f - 0.02f * (float)i;
		lay->visible = (i % 3 != 0);
	}
	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_eq(flis_layer_count(), n);

	int i = 0;
	for (GSList *l = com.uniq->layers; l; l = l->next, i++) {
		flis_layer_t *r = (flis_layer_t *)l->data;
		gchar *expected = g_strdup_printf("L%02d", i);
		cr_assert_str_eq(r->layer_name, expected);
		cr_assert_float_eq(r->fit->fdata[0], 0.05f * (float)i, 1e-5);
		cr_assert_float_eq(r->opacity, 1.0f - 0.02f * (float)i, 1e-5);
		cr_assert_eq(r->visible, (i % 3 != 0));
		g_free(expected);
	}
}

/* Round-trip every blend mode — confirms the BLEND_TABLE mapping
 * round-trips losslessly via the 16A string column.  Existing per-mode
 * kernel tests verify the math; this verifies the metadata. */
Test(flis_stress, every_blend_mode_round_trips) {
	const flis_blend_mode_t modes[] = {
		FLIS_BLEND_NORMAL, FLIS_BLEND_MULTIPLY, FLIS_BLEND_SCREEN,
		FLIS_BLEND_OVERLAY, FLIS_BLEND_SOFT_LIGHT, FLIS_BLEND_HARD_LIGHT,
		FLIS_BLEND_COLOR_DODGE, FLIS_BLEND_COLOR_BURN, FLIS_BLEND_DARKEN,
		FLIS_BLEND_LIGHTEN, FLIS_BLEND_DIFFERENCE, FLIS_BLEND_EXCLUSION,
		FLIS_BLEND_HUE, FLIS_BLEND_SATURATION, FLIS_BLEND_COLOR,
		FLIS_BLEND_LUMINOSITY, FLIS_BLEND_CHROMA,
	};
	const int n_modes = G_N_ELEMENTS(modes);

	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	for (int i = 0; i < n_modes; i++) {
		gchar *name = g_strdup_printf("m%d", i);
		flis_layer_t *lay = flis_test_add_layer(
		    flis_test_make_mono_fits(4, 4, 0.5f), name);
		g_free(name);
		lay->blend_mode = modes[i];
	}

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	GSList *l = com.uniq->layers->next;  /* skip base */
	for (int i = 0; i < n_modes && l; i++, l = l->next) {
		flis_layer_t *r = (flis_layer_t *)l->data;
		cr_assert_eq(r->blend_mode, modes[i],
		    "blend mode %d round-trip mismatch: wrote %d, read %d",
		    i, (int)modes[i], (int)r->blend_mode);
	}
}

/* Group config round-trip: PASS_THROUGH and NORMAL group blend modes,
 * with members on either side of the membership.  Verifies that group
 * blend_mode, opacity, visibility, and member assignment all survive
 * a save/load cycle. */
Test(flis_stress, group_configurations_round_trip) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	flis_layer_t *l1 = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.25f), "in_pt");
	flis_layer_t *l2 = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.50f), "in_normal");
	flis_layer_t *l3 = flis_test_add_layer(
	    flis_test_make_mono_fits(4, 4, 0.75f), "ungrouped");

	flis_group_t *g_pt = flis_group_add("g_pt");
	flis_group_t *g_n  = flis_group_add("g_n");
	g_n->blend_mode = FLIS_BLEND_NORMAL;
	g_n->opacity = 0.6f;
	g_n->visible = FALSE;

	l1->group_id = g_pt->item_id;
	l2->group_id = g_n->item_id;
	(void)l3;

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	flis_free_groups(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert_eq(g_slist_length(com.uniq->groups), 2u);
	cr_assert_eq(flis_layer_count(), 4);

	flis_group_t *r_pt = (flis_group_t *)com.uniq->groups->data;
	flis_group_t *r_n  = (flis_group_t *)com.uniq->groups->next->data;
	cr_assert_eq(r_pt->blend_mode, FLIS_BLEND_PASS_THROUGH);
	cr_assert_eq(r_n->blend_mode, FLIS_BLEND_NORMAL);
	cr_assert_float_eq(r_n->opacity, 0.6f, 1e-5);
	cr_assert_not(r_n->visible);

	/* Member assignment preserved */
	flis_layer_t *r_l1 = flis_layer_get_by_name("in_pt");
	flis_layer_t *r_l2 = flis_layer_get_by_name("in_normal");
	flis_layer_t *r_l3 = flis_layer_get_by_name("ungrouped");
	cr_assert_eq(r_l1->group_id, r_pt->item_id);
	cr_assert_eq(r_l2->group_id, r_n->item_id);
	cr_assert_eq(r_l3->group_id, 0);
}

/* Stack-mutation stress: add and remove 32 layers in alternation,
 * confirming the item_id allocator never reuses an id and that the
 * layer-order steps remain stable. */
Test(flis_stress, alternating_add_remove_keeps_ids_unique) {
	GHashTable *seen_ids = g_hash_table_new(g_direct_hash, g_direct_equal);
	for (int i = 0; i < 32; i++) {
		gchar *name = g_strdup_printf("L%d", i);
		flis_layer_t *lay = flis_test_add_layer(
		    flis_test_make_mono_fits(2, 2, 0.1f), name);
		g_free(name);
		cr_assert_not_null(lay);
		cr_assert(!g_hash_table_contains(seen_ids,
		    GINT_TO_POINTER(lay->item_id)),
		    "item_id %d reused — must be globally unique", lay->item_id);
		g_hash_table_insert(seen_ids, GINT_TO_POINTER(lay->item_id), lay);
		/* Every 4th iteration, remove the most-recently-added layer. */
		if (i > 0 && i % 4 == 0)
			flis_layer_remove(lay);
	}
	g_hash_table_destroy(seen_ids);
}
