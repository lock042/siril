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
#include "core/nde_history.h"
#include "core/nde_snapstore.h"

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
	fake->mask_fd = -1;   /* snap NULL via g_new0 — no pixels */
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

/* ----- restore path: undo/redo round-trips (C3/C4 regressions) -------- */

/* Redo of a props-only entry must reapply the property change.  Before
 * the counterpart-capture fix, undo_display_data pushed a flavour-blind
 * pixel snapshot to the redo stack, so redo restored identical pixels
 * and the property change silently vanished. */
Test(flis_undo_gui, props_undo_redo_round_trip) {
	com.headless = TRUE;
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	uniq_set_active_layer(com.uniq, 0);

	l->opacity = 0.3f;
	cr_assert_eq(undo_save_flis_layer_props(l, "opacity"), 0);
	l->opacity = 0.9f;

	cr_assert_eq(undo_display_data(UNDO), 0);
	cr_assert_float_eq(l->opacity, 0.3f, 1e-6, "undo must restore the property");

	cr_assert_eq(g_list_length(com.redo_stack), 1);
	historic *h = (historic *)com.redo_stack->data;
	cr_assert_not_null(h->layer_props,
	                   "redo counterpart must be props-flavoured, not a pixel snapshot");
	cr_assert_null(h->snap, "props counterpart must not carry a pixel snapshot");

	cr_assert_eq(undo_display_data(REDO), 0);
	cr_assert_float_eq(l->opacity, 0.9f, 1e-6, "redo must reapply the property change");
	/* And the undo stack got a fresh props counterpart for another undo. */
	cr_assert_eq(g_list_length(com.undo_stack), 1);
	cr_assert_eq(undo_display_data(UNDO), 0);
	cr_assert_float_eq(l->opacity, 0.3f, 1e-6, "second undo must work after redo");
}

/* A plain pixel entry saved while layer A was active must restore into
 * layer A even if the active layer has changed since (C4: entries now
 * record the layer identity; before, undo wrote A's pixels into whatever
 * layer was active, resizing it to A's dims). */
Test(flis_undo_gui, plain_undo_restores_into_saved_layer) {
	com.headless = TRUE;
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.25f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(6, 6, 0.75f), "b");
	uniq_set_active_layer(com.uniq, 0);   /* a active, gfit = a->fit */

	cr_assert_eq(undo_save_state(gfit, "stretch"), 0);
	cr_assert_eq(g_list_length(com.undo_stack), 1);
	historic *h = (historic *)com.undo_stack->data;
	cr_assert_eq(h->flis_layer_id, a->item_id,
	             "plain entry must record the layer it was saved from");

	for (size_t i = 0; i < 16; i++)
		a->fit->fdata[i] = 0.6f;

	uniq_set_active_layer(com.uniq, 1);   /* b now active, gfit = b->fit */
	cr_assert_eq(undo_display_data(UNDO), 0);

	cr_assert_float_eq(a->fit->fdata[0], 0.25f, 1e-6,
	                   "undo must restore the layer the state was saved from");
	cr_assert_eq(a->fit->rx, 4u);
	cr_assert_float_eq(b->fit->fdata[0], 0.75f, 1e-6,
	                   "undo must not clobber the now-active layer");
	cr_assert_eq(b->fit->rx, 6u, "active layer must not be resized by the restore");

	/* Redo counterpart must also target layer a. */
	cr_assert_eq(g_list_length(com.redo_stack), 1);
	historic *r = (historic *)com.redo_stack->data;
	cr_assert_eq(r->flis_layer_id, a->item_id);
	cr_assert_eq(undo_display_data(REDO), 0);
	cr_assert_float_eq(a->fit->fdata[0], 0.6f, 1e-6, "redo must reapply into layer a");
	cr_assert_float_eq(b->fit->fdata[0], 0.75f, 1e-6);
}

/* ----- NDE provenance coupling (nde sketch §13.3) ----- */

Test(flis_undo_gui, nde_tag_top_entry) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	undo_tag_top_nde_record(7);           /* empty stack: no-op, no crash */
	cr_assert_eq(undo_save_flis_layer_props(l, "op"), 0);
	historic *top = (historic *)com.undo_stack->data;
	cr_assert_eq(top->nde_record_id, 0, "fresh entries start uncoupled");
	undo_tag_top_nde_record(0);           /* id 0: no-op */
	cr_assert_eq(top->nde_record_id, 0);
	undo_tag_top_nde_record(42);
	cr_assert_eq(top->nde_record_id, 42);
	nde_history_attach(NULL);
}

Test(flis_undo_gui, nde_undo_redo_moves_live_count) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	nde_record *r1 = nde_record_new();
	r1->op_id = g_strdup("a.a");
	nde_history_append(r1);
	nde_record *r2 = nde_record_new();
	r2->op_id = g_strdup("b.b");
	gint64 id2 = nde_history_append(r2);
	cr_assert_eq(nde_history_live_count(), 2);

	cr_assert_eq(undo_save_flis_layer_props(l, "op"), 0);
	undo_tag_top_nde_record(id2);

	cr_assert_eq(undo_display_data(UNDO), 0);
	cr_assert_eq(nde_history_live_count(), 1, "undo must retire the coupled record");
	/* the counterpart pushed to the redo stack inherited the coupling */
	cr_assert_eq(((historic *)com.redo_stack->data)->nde_record_id, id2);

	cr_assert_eq(undo_display_data(REDO), 0);
	cr_assert_eq(nde_history_live_count(), 2, "redo must revive the coupled record");
	nde_history_attach(NULL);
}

/* ----- C2: pixel snapshots live in the shared store (convergence) ----- */

Test(flis_undo_gui, snapstore_pre_post_tags_follow_undo_redo) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.25f), "a");
	gfit = a->fit;
	uniq_set_active_layer(com.uniq, 0);

	/* a pixel op: record + coupled pixel undo entry (the worker's order:
	 * capture, save undo, tag) */
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("geometry.mirrorx");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("x_axis=1");
	gint64 id = nde_history_append(rec);
	cr_assert(id > 0);
	cr_assert_eq(undo_save_state(gfit, "op"), 0);
	undo_tag_top_nde_record(id);

	gint item = a->item_id;
	cr_assert(nde_snapstore_has(item, id, FALSE),
	          "a tagged undo entry must register its snapshot as PRE(record)");
	cr_assert(!nde_snapstore_has(item, id, TRUE));

	/* mutate then undo: the popped entry dies (PRE deregisters) and the
	 * redo counterpart registers the current state as POST(record) */
	gfit->fdata[0] = 0.9f;
	cr_assert_eq(undo_display_data(UNDO), 0);
	cr_assert(!nde_snapstore_has(item, id, FALSE),
	          "the freed undo entry's PRE tag must vanish with it");
	cr_assert(nde_snapstore_has(item, id, TRUE),
	          "the redo counterpart must register POST(record)");

	/* the POST snapshot holds the mutated (post-op) pixels */
	fits *post = nde_snapstore_lookup(item, id, TRUE);
	cr_assert_not_null(post);
	cr_assert_float_eq(post->fdata[0], 0.9f, 1e-6,
	                   "POST must hold the state after the record");
	clearfits(post); free(post);

	/* redo: POST dies with the popped redo entry; the fresh undo
	 * counterpart re-registers PRE */
	cr_assert_eq(undo_display_data(REDO), 0);
	cr_assert(!nde_snapstore_has(item, id, TRUE));
	cr_assert(nde_snapstore_has(item, id, FALSE),
	          "redo's undo-stack counterpart must register PRE(record)");
	fits *pre = nde_snapstore_lookup(item, id, FALSE);
	cr_assert_not_null(pre);
	cr_assert_float_eq(pre->fdata[0], 0.25f, 1e-6,
	                   "PRE must hold the state before the record");
	clearfits(pre); free(pre);

	/* flush releases everything — the no-meta-undo invalidation path */
	undo_flush();
	cr_assert(!nde_snapstore_has(item, id, FALSE));
	cr_assert(!nde_snapstore_has(item, id, TRUE));

	nde_history_attach(NULL);
	gfit = NULL;
}
