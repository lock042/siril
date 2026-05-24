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
 * test_flis_roundtrip — full save/load round-trip of a synthetic FLIS,
 * verifying that metadata, pixel data, layer ordering, and FLIS-specific
 * keys (LAYER_COLOR tint, blend modes, opacities, group membership,
 * unknown METADATA entries) are preserved byte-for-byte through one
 * save_flis → load_flis cycle.
 */

#include <criterion/criterion.h>
#include <unistd.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

static char *tmpdir = NULL;
static char *tmppath = NULL;

static void setup(void) {
	flis_test_init_com();
	tmpdir = g_dir_make_tmp("flis-rt-XXXXXX", NULL);
	tmppath = g_build_filename(tmpdir, "test.flis", NULL);
}

static void teardown(void) {
	if (tmppath) { g_unlink(tmppath); g_free(tmppath); tmppath = NULL; }
	if (tmpdir)  { g_rmdir(tmpdir);   g_free(tmpdir);   tmpdir  = NULL; }
	flis_test_cleanup_com();
}

TestSuite(flis_roundtrip, .init = setup, .fini = teardown);

/* Smallest possible round-trip: 1 layer with default everything.  Validates
 * the core save → load path independently of any FLIS-specific feature. */
Test(flis_roundtrip, single_layer_minimal) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	cr_assert_not_null(l);
	cr_assert_eq(save_flis(tmppath), 0, "save_flis failed");

	/* Tear down the in-memory state and reload from disk */
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0, "load_flis failed");
	cr_assert_eq(flis_layer_count(), 1);

	flis_layer_t *reloaded = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_str_eq(reloaded->layer_name, "base");
	cr_assert_eq(reloaded->fit->rx, 8);
	cr_assert_eq(reloaded->fit->ry, 8);
	cr_assert_eq(reloaded->fit->naxes[2], 1);
	cr_assert_float_eq(reloaded->fit->fdata[0], 0.25f, 1e-5);
}

/* Two-layer round-trip with non-default blend mode and opacity on top.
 * Verifies that blend_mode + opacity + visibility serialise correctly. */
Test(flis_roundtrip, two_layers_with_blend_and_opacity) {
	flis_layer_t *base = flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.1f, 0.2f, 0.3f), "OSC");
	flis_layer_t *top  = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.7f), "Ha");
	cr_assert_eq(flis_layer_set_blend_mode(top, FLIS_BLEND_SCREEN), 0);
	cr_assert_eq(flis_layer_set_opacity(top, 0.8f), 0);

	/* Snapshot IDs before save/free — `base` and `top` become dangling
	 * pointers after flis_free_layers and must not be dereferenced. */
	gint base_id = base->item_id;
	gint top_id  = top->item_id;

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_eq(flis_layer_count(), 2);

	flis_layer_t *r_base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *r_top  = (flis_layer_t *)com.uniq->layers->next->data;

	cr_assert_str_eq(r_base->layer_name, "OSC");
	cr_assert_eq(r_base->fit->naxes[2], 3);
	cr_assert_float_eq(r_base->fit->fpdata[0][0], 0.1f, 1e-5);
	cr_assert_float_eq(r_base->fit->fpdata[1][0], 0.2f, 1e-5);
	cr_assert_float_eq(r_base->fit->fpdata[2][0], 0.3f, 1e-5);

	cr_assert_str_eq(r_top->layer_name, "Ha");
	cr_assert_eq(r_top->fit->naxes[2], 1);
	cr_assert_float_eq(r_top->fit->fdata[0], 0.7f, 1e-5);
	cr_assert_eq(r_top->blend_mode, FLIS_BLEND_SCREEN);
	cr_assert_float_eq(r_top->opacity, 0.8f, 1e-5);

	/* base layer item_id should be preserved */
	cr_assert_eq(r_base->item_id, base_id);
	cr_assert_eq(r_top->item_id, top_id);
}

Test(flis_roundtrip, mono_layer_tint_round_trips) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	flis_layer_t *ha = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "Ha");
	cr_assert_eq(flis_layer_set_tint(ha, 1.0, 0.2, 0.1), 0);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	flis_layer_t *r_ha = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert(r_ha->has_tint, "tint should survive round-trip");
	cr_assert_float_eq(r_ha->layer_tint.r, 1.0, 1e-5);
	cr_assert_float_eq(r_ha->layer_tint.g, 0.2, 1e-5);
	cr_assert_float_eq(r_ha->layer_tint.b, 0.1, 1e-5);
}

Test(flis_roundtrip, locked_flag_survives) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	cr_assert_eq(flis_layer_set_locked(l, TRUE), 0);
	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);
	flis_layer_t *r = (flis_layer_t *)com.uniq->layers->data;
	cr_assert(r->locked);
}

Test(flis_roundtrip, group_membership_survives) {
	flis_layer_t *base = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	flis_layer_t *top  = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "narrow");
	flis_group_t *g = flis_group_add("Narrowband");
	cr_assert_not_null(g);
	cr_assert_eq(flis_layer_set_group(top, g->item_id), 0);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	flis_free_groups(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert_eq(flis_group_count(), 1);
	flis_layer_t *r_base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *r_top  = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(r_base->group_id, 0);
	cr_assert_neq(r_top->group_id, 0);

	flis_group_t *r_g = (flis_group_t *)com.uniq->groups->data;
	cr_assert_str_eq(r_g->name, "Narrowband");
	cr_assert_eq(r_top->group_id, r_g->item_id);
	(void)base;  /* unused after free */
}

Test(flis_roundtrip, layer_order_preserved) {
	/* Add three layers, then save/load: the order they come back in must
	 * match (the lowest layer_order is at the head of com.uniq->layers). */
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.1f), "first");
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.2f), "second");
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.3f), "third");

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	GSList *l = com.uniq->layers;
	cr_assert_str_eq(((flis_layer_t *)l->data)->layer_name, "first");
	l = l->next;
	cr_assert_str_eq(((flis_layer_t *)l->data)->layer_name, "second");
	l = l->next;
	cr_assert_str_eq(((flis_layer_t *)l->data)->layer_name, "third");
}

Test(flis_roundtrip, position_offset_survives) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "canvas");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "patch");
	top->position_x = 4;
	top->position_y = 3;

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	flis_layer_t *r_top = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(r_top->position_x, 4);
	cr_assert_eq(r_top->position_y, 3);
}

/* §6.3 forward-compat: any METADATA key the loader does not recognise
 * (e.g. attributes added by a future Siril release) is stashed on the
 * layer / group and re-emitted on save so the round-trip is lossless. */
Test(flis_roundtrip, unknown_layer_metadata_round_trips) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "base");
	/* Simulate a future-version key that this build doesn't understand. */
	base->unknown_metadata = g_strdup("FUTURE_KEY=hello world;ANOTHER=42");

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	flis_layer_t *r = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_not_null(r->unknown_metadata,
	    "unknown metadata must be preserved across save/load");
	/* The exact ordering matches the build/parse contract: known keys
	 * first, then unknown.  Pair order within unknown is preserved. */
	cr_assert(strstr(r->unknown_metadata, "FUTURE_KEY=hello world") != NULL,
	    "FUTURE_KEY value preserved (got: %s)", r->unknown_metadata);
	cr_assert(strstr(r->unknown_metadata, "ANOTHER=42") != NULL,
	    "ANOTHER value preserved (got: %s)", r->unknown_metadata);
}

Test(flis_roundtrip, unknown_group_metadata_round_trips) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "base");
	flis_group_t *grp = flis_group_add("g1");
	cr_assert_not_null(grp);
	grp->unknown_metadata = g_strdup("FUTURE_GROUP_KEY=alpha;BLEND_PROFILE=lab");

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	flis_free_groups(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert_not_null(com.uniq->groups);
	flis_group_t *r_grp = (flis_group_t *)com.uniq->groups->data;
	cr_assert_not_null(r_grp->unknown_metadata);
	cr_assert(strstr(r_grp->unknown_metadata, "FUTURE_GROUP_KEY=alpha") != NULL);
	cr_assert(strstr(r_grp->unknown_metadata, "BLEND_PROFILE=lab") != NULL);
}
