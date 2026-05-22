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
 * test_flis_undo — exercises the FLIS undo machinery:
 *   - the headless gate (com.script): undo_save_flis_* must return 0
 *     immediately without pushing anything to the stack
 *   - the GUI-mode push path: undo_save_flis_layer_props /
 *     undo_save_flis_lmask / undo_save_flis_layer_reorder /
 *     undo_save_flis_multi_layer_props all add to com.undo_stack and
 *     clear com.redo_stack
 *   - flis_undo_purge_layer removes entries belonging to the deleted
 *     layer from both stacks
 *   - compound (multi-layer) entries shrink in place when one of their
 *     sub-entries matches the deleted layer
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/undo.h"

cominfo com;
fits *gfit;

/* Override the default `com.script = TRUE` from flis_test_init_com — the
 * undo save paths use that as their headless gate.  GUI-mode tests need
 * com.script == FALSE so the save paths actually push entries. */
static void setup_gui_mode(void) {
	flis_test_init_com();
	com.script = FALSE;
}
static void setup_headless_mode(void) {
	flis_test_init_com();
	com.script = TRUE;
}

TestSuite(flis_undo_headless, .init = setup_headless_mode, .fini = flis_test_cleanup_com);
TestSuite(flis_undo_gui,      .init = setup_gui_mode,      .fini = flis_test_cleanup_com);

/* ----- headless: every save function returns 0 and pushes nothing ----- */

Test(flis_undo_headless, props_save_is_noop_when_script) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	cr_assert_eq(undo_save_flis_layer_props(l, "test"), 0);
	cr_assert_null(com.undo_stack, "headless save must not push to undo_stack");
}

Test(flis_undo_headless, lmask_save_is_noop_when_script) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	cr_assert_eq(undo_save_flis_lmask(l, "test"), 0);
	cr_assert_null(com.undo_stack);
}

Test(flis_undo_headless, reorder_save_is_noop_when_script) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	cr_assert_eq(undo_save_flis_layer_reorder(a, b, "test"), 0);
	cr_assert_null(com.undo_stack);
}

Test(flis_undo_headless, multi_layer_save_is_noop_when_script) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	GSList *layers = g_slist_append(NULL, a);
	layers = g_slist_append(layers, b);
	cr_assert_eq(undo_save_flis_multi_layer_props(layers, "test"), 0);
	cr_assert_null(com.undo_stack);
	g_slist_free(layers);
}

/* ----- GUI-mode: each save function pushes one entry and clears redo --- */

Test(flis_undo_gui, props_save_pushes_one_entry) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	cr_assert_eq(undo_save_flis_layer_props(l, "change opacity"), 0);
	cr_assert_eq(g_list_length(com.undo_stack), 1);
	historic *h = (historic *)com.undo_stack->data;
	cr_assert_not_null(h->layer_props);
	cr_assert_eq(h->flis_layer_id, l->item_id);
	cr_assert_str_eq(h->history, "change opacity");
}

Test(flis_undo_gui, save_clears_redo_stack) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	/* Pre-populate redo with a fake entry (fd/mask_fd = -1 to look like a
	 * filename-mode entry with no actual file; otherwise undo_free_item
	 * would try to g_close(0) and emit a GLib warning). */
	historic *fake = g_new0(historic, 1);
	fake->fd = -1;
	fake->mask_fd = -1;
	com.redo_stack = g_list_prepend(com.redo_stack, fake);
	cr_assert_eq(g_list_length(com.redo_stack), 1);

	cr_assert_eq(undo_save_flis_layer_props(l, "x"), 0);
	cr_assert_null(com.redo_stack, "new save must invalidate redo branch");
	cr_assert_eq(g_list_length(com.undo_stack), 1);
}

Test(flis_undo_gui, reorder_save_records_both_layers) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	gint a_order = a->layer_order;
	gint b_order = b->layer_order;
	cr_assert_eq(undo_save_flis_layer_reorder(a, b, "swap"), 0);
	historic *h = (historic *)com.undo_stack->data;
	cr_assert_eq(h->reorder_layer_a_id, a->item_id);
	cr_assert_eq(h->reorder_layer_a_order, a_order);
	cr_assert_eq(h->reorder_layer_b_id, b->item_id);
	cr_assert_eq(h->reorder_layer_b_order, b_order);
}

Test(flis_undo_gui, multi_layer_props_save_holds_one_entry_per_layer) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	flis_layer_t *c = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "c");
	GSList *layers = g_slist_append(NULL, a);
	layers = g_slist_append(layers, b);
	layers = g_slist_append(layers, c);

	cr_assert_eq(undo_save_flis_multi_layer_props(layers, "group drag"), 0);
	historic *h = (historic *)com.undo_stack->data;
	cr_assert_eq(h->n_multi_entries, 3);
	cr_assert_eq(h->multi_entries[0].flis_layer_id, a->item_id);
	cr_assert_eq(h->multi_entries[1].flis_layer_id, b->item_id);
	cr_assert_eq(h->multi_entries[2].flis_layer_id, c->item_id);
	for (guint i = 0; i < 3; i++) {
		cr_assert(h->multi_entries[i].props_only,
		          "every entry should be props-only");
		cr_assert_not_null(h->multi_entries[i].layer_props);
	}
	g_slist_free(layers);
}

/* ----- purge_layer removes entries pointing at the deleted layer ------- */

Test(flis_undo_gui, purge_removes_matching_single_entries) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");

	/* Two undo entries: one for a, one for b */
	undo_save_flis_layer_props(a, "change a");
	undo_save_flis_layer_props(b, "change b");
	cr_assert_eq(g_list_length(com.undo_stack), 2);

	flis_undo_purge_layer(b->item_id);
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "purge must remove the entry pointing at b");
	historic *survivor = (historic *)com.undo_stack->data;
	cr_assert_eq(survivor->flis_layer_id, a->item_id,
	             "surviving entry must belong to a");
}

Test(flis_undo_gui, purge_compacts_compound_entries) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	flis_layer_t *c = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "c");
	GSList *layers = g_slist_append(NULL, a);
	layers = g_slist_append(layers, b);
	layers = g_slist_append(layers, c);
	undo_save_flis_multi_layer_props(layers, "group drag");
	g_slist_free(layers);

	historic *h = (historic *)com.undo_stack->data;
	cr_assert_eq(h->n_multi_entries, 3);

	flis_undo_purge_layer(b->item_id);
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "compound entry should survive (still has a and c)");
	cr_assert_eq(h->n_multi_entries, 2,
	             "compound entry should shrink from 3 to 2");
	/* a and c remain; b's sub-entry is gone */
	cr_assert_eq(h->multi_entries[0].flis_layer_id, a->item_id);
	cr_assert_eq(h->multi_entries[1].flis_layer_id, c->item_id);
}

Test(flis_undo_gui, purge_drops_compound_when_all_match) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	GSList *layers = g_slist_append(NULL, a);
	layers = g_slist_append(layers, b);
	undo_save_flis_multi_layer_props(layers, "group drag");
	g_slist_free(layers);

	/* Purge both layers in turn */
	flis_undo_purge_layer(a->item_id);
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "still 1 entry — compound has b left");
	flis_undo_purge_layer(b->item_id);
	cr_assert_null(com.undo_stack,
	               "compound entry should be removed once all sub-entries match");
}

Test(flis_undo_gui, purge_walks_both_stacks) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	undo_save_flis_layer_props(a, "step1");
	undo_save_flis_layer_props(a, "step2");

	/* Manually move one entry from undo to redo */
	historic *moved = (historic *)com.undo_stack->data;
	com.undo_stack = g_list_remove(com.undo_stack, moved);
	com.redo_stack = g_list_prepend(com.redo_stack, moved);
	cr_assert_eq(g_list_length(com.undo_stack), 1);
	cr_assert_eq(g_list_length(com.redo_stack), 1);

	flis_undo_purge_layer(a->item_id);
	cr_assert_null(com.undo_stack, "undo entries for a should be gone");
	cr_assert_null(com.redo_stack, "redo entries for a should be gone");
}

Test(flis_undo_gui, purge_with_layer_none_is_noop) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	undo_save_flis_layer_props(a, "step1");
	cr_assert_eq(g_list_length(com.undo_stack), 1);

	flis_undo_purge_layer(FLIS_UNDO_LAYER_NONE);
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "purge with sentinel must touch nothing");
}
