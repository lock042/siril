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
 * test_flis_mask_route — exercises the §5.2 FLIS layermask routing in
 * generic_mask_worker.  When args->target_layer_id is set and the
 * mask_hook produces a mask on args->fit->mask, the worker moves that
 * mask onto the named layer's lmask slot and clears args->fit->mask.
 *
 * The mask creation path itself (mask_from_channel_hook etc.) is
 * tested elsewhere; this file uses a trivial test hook that fills a
 * known mask so the assertions can focus on routing semantics.
 *
 * The worker is called directly on the test thread (com.script +
 * args->command run the synchronous-completion / no-undo path).
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/masks.h"
#include "core/undo.h"
#include "core/op_descriptor.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_replay.h"
#include "core/nde/nde_script_scope.h"
#include "core/processing_thread.h"

cominfo com;
fits *gfit;

TestSuite(flis_mask_route, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Test hook: creates a constant 8-bit mask on the fit. */
static int fill_mask_hook(struct generic_mask_args *args) {
	fits *f = args->fit;
	size_t n = (size_t)f->rx * f->ry;
	f->mask = calloc(1, sizeof(mask_t));
	f->mask->bitpix = 8;
	f->mask->data = malloc(n);
	memset(f->mask->data, 0xAA, n);
	return 0;
}

/* Hook that always fails — used to verify a failed op records nothing. */
static int failing_hook(struct generic_mask_args *args) {
	(void)args;
	return 1;
}

/* Hook that does nothing — used to verify "no routing without -layermask=". */
static int noop_hook(struct generic_mask_args *args) {
	(void)args;
	return 0;
}

/* Routing happens: mask hook produces fit->mask, worker moves it to the
 * target layer's lmask and clears the processing mask. */
Test(flis_mask_route, target_id_set_routes_to_layer_lmask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.5f), "patch");
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();
	cr_assert_not_null(gfit);

	cr_assert_null(top->lmask, "precondition: top has no lmask");

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test routing to top";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = top->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	/* Mask now lives on top->lmask; gfit->mask cleared. */
	cr_assert_null(gfit->mask, "processing mask should be cleared");
	cr_assert_not_null(top->lmask, "target layer lmask should be populated");
	cr_assert_eq(top->lmask->w, 16u);
	cr_assert_eq(top->lmask->h, 16u);
	cr_assert_eq(top->lmask->bitpix, 8);
	const uint8_t *p = (const uint8_t *)top->lmask->data;
	cr_assert_eq(p[0], 0xAA);
	cr_assert_eq(p[16 * 16 - 1], 0xAA);
}

/* Routing target a different (base) layer than the active layer:
 * routing key is item_id, not "active layer". */
Test(flis_mask_route, target_id_can_be_non_active_layer) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.5f), "patch");
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test routing to base";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = base->item_id;  /* not the active one */
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask);
	cr_assert_not_null(base->lmask);
	cr_assert_null(top->lmask,  "non-target layer must remain untouched");
}

/* When target_layer_id == 0 (the default), the worker leaves the mask on
 * the processing slot and activates it (mask_creation path). */
Test(flis_mask_route, target_id_zero_keeps_processing_mask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test default routing";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = 0;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_not_null(gfit->mask, "processing mask should remain");
	cr_assert_eq(gfit->mask->bitpix, 8);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_null(base->lmask, "no lmask routing should have happened");
}

/* Dimension mismatch: the target layer is smaller than args->fit.  The
 * worker's pre-check refuses the whole operation up front — no mask is
 * generated at all.  (The old behaviour generated the mask and left it
 * on the processing slot, which surprised users with an error message
 * AND a stray mask they never asked for.) */
Test(flis_mask_route, dim_mismatch_refuses_without_creating_mask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *small = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "small");
	uniq_set_active_layer(com.uniq, 0);  /* base active (16x16) */
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "mismatch test";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = small->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask,
	               "refused op must not leave a stray processing mask");
	cr_assert_null(small->lmask, "mismatch layer must not receive a mask");
}

/* No mask produced by the hook (hook returned success but didn't set
 * fit->mask): routing block is skipped, no crash. */
Test(flis_mask_route, no_mask_produced_no_op) {
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "x");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = noop_hook;
	args->description      = "no mask test";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = top->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask);
	cr_assert_null(top->lmask);
}

/* Routed lmask creation must be undoable: the worker saves an
 * lmask-flavoured entry (not a flavour-blind pixel snapshot, which left
 * the freshly routed mask surviving undo), and redo reinstates it. */
Test(flis_mask_route, routed_lmask_undo_removes_mask) {
	com.script = FALSE;    /* GUI-mode: worker-side undo saves are live */
	com.headless = TRUE;
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "routed undo test";
	args->command          = FALSE;   /* GUI path: worker saves undo */
	args->mask_creation    = TRUE;
	args->target_layer_id  = base->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_not_null(base->lmask, "mask must be routed to the lmask");
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "routed op must push exactly one (lmask) undo entry");

	cr_assert_eq(undo_display_data(UNDO), 0);
	cr_assert_null(base->lmask, "undo must remove the routed lmask");

	cr_assert_eq(undo_display_data(REDO), 0);
	cr_assert_not_null(base->lmask, "redo must reinstate the lmask");
}

/* A refused routed op (dimension mismatch) must leave no undo entry —
 * the pre-check fires before any undo state is saved. */
Test(flis_mask_route, refused_routing_leaves_no_undo) {
	com.script = FALSE;
	com.headless = TRUE;
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *small = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "small");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "refused undo test";
	args->command          = FALSE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = small->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask, "refused op must create nothing");
	cr_assert_null(small->lmask);
	cr_assert_null(com.undo_stack,
	               "refused op must leave no undo entry");
}

/* ---------------- NDE provenance for mask operations ---------------- */

/* A descriptor so the worker takes the capture path; no serializer, which is
 * the state of every real mask op today (Tier B). */
static const op_descriptor op_desc_test_mask = {
	.id = "mask.test_fill", .version = 1, .mask_hook = fill_mask_hook,
	.description = "Test mask fill", .mem_ratio = 0.0f, .flags = 0,
};

static const nde_record *only_record(void) {
	static GPtrArray *snap = NULL;
	if (snap) g_ptr_array_unref(snap);
	snap = nde_history_snapshot(NULL);
	return (snap && snap->len) ? g_ptr_array_index(snap, snap->len - 1) : NULL;
}

static void run_mask_op(fits *fit, gint target_layer_id) {
	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit             = fit;
	args->op              = &op_desc_test_mask;
	args->command         = TRUE;
	args->mask_creation   = TRUE;
	args->target_layer_id = target_layer_id;
	args->max_threads     = 1;
	generic_mask_worker(args);
}

/* A mask op records against the MASK, not the layer whose pixels it will
 * later modulate — that separation is what step 4's input pins need. */
Test(flis_mask_route, mask_op_records_against_the_processing_mask_item) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	run_mask_op(gfit, 0);
	cr_assert_not_null(gfit->mask, "precondition: the op left a processing mask");

	gint pmask = flis_layer_pmask_id(base);
	cr_assert_neq(pmask, 0);
	cr_assert_neq(pmask, base->item_id, "the mask is its own item");
	const nde_record *rec = only_record();
	cr_assert_not_null(rec, "a mask op must leave provenance");
	cr_assert_str_eq(rec->op_id, "mask.test_fill");
	cr_assert_eq(rec->target_item_id, pmask);
	cr_assert_eq(rec->scope, NDE_SCOPE_LAYER);
	cr_assert_eq(rec->tier, NDE_TIER_B, "no mask op has a serializer yet");
	cr_assert(!rec->mask_active, "a mask op has no mask input of its own");
}

Test(flis_mask_route, routed_mask_op_records_against_the_layer_mask_item) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.5f), "patch");
	uniq_set_active_layer(com.uniq, 1);
	gfit = flis_active_layer_fit();

	run_mask_op(gfit, top->item_id);
	cr_assert_not_null(top->lmask, "precondition: the mask was routed");

	gint lmask = flis_layer_lmask_id(top);
	const nde_record *rec = only_record();
	cr_assert_not_null(rec);
	cr_assert_eq(rec->target_item_id, lmask,
	             "a routed op belongs to the layer mask it produced, "
	             "not the processing slot it passed through");
}

/* The mask's lineage is separate from the layer's: a mask op must not join
 * the layer's replay chain, or every masked workflow would drag mask edits
 * into the pixel history. */
Test(flis_mask_route, mask_records_stay_out_of_the_layer_chain) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	run_mask_op(gfit, 0);

	nde_chain *layer_chain = nde_chain_build(base->item_id);
	cr_assert_eq(layer_chain->records->len, 0,
	             "the mask op is not part of the layer's pixel lineage");
	nde_chain_free(layer_chain);

	nde_chain *mask_chain = nde_chain_build(flis_layer_pmask_id(base));
	cr_assert_eq(mask_chain->records->len, 1,
	             "but it IS the mask's own chain");
	nde_chain_free(mask_chain);
}

/* On a plain image the processing mask has no layer to hang off, so it uses
 * the NDE_ITEM_PLAIN_MASK sentinel — the mask counterpart of the -1 that
 * names the image itself. */
Test(flis_mask_route, plain_image_mask_op_uses_the_sentinel_item) {
	fits *f = flis_test_make_mono_fits(16, 16, 0.25f);
	gfit = f;
	/* no com.uniq->layers: this is a plain single image */
	cr_assert(!is_current_image_flis());

	run_mask_op(gfit, 0);
	const nde_record *rec = only_record();
	cr_assert_not_null(rec);
	cr_assert_eq(rec->target_item_id, NDE_ITEM_PLAIN_MASK);
	cr_assert_neq(rec->target_item_id, NDE_ITEM_IMAGE,
	              "the mask is not the image");
	gfit = NULL;
	clearfits(f);
	free(f);
}

/* A failed hook must not record: the log describes changes that happened. */
Test(flis_mask_route, a_failed_mask_op_records_nothing) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	static const op_descriptor failing = {
		.id = "mask.test_fail", .version = 1, .mask_hook = failing_hook,
		.description = "Failing mask op", .mem_ratio = 0.0f, .flags = 0,
	};
	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit           = gfit;
	args->op            = &failing;
	args->command       = TRUE;
	args->mask_creation = TRUE;
	args->max_threads   = 1;
	generic_mask_worker(args);

	cr_assert_eq(nde_history_live_count(), 0);
}

/* A Tier-C script re-run's mask commands travel through this worker under
 * SLOT_REPLAY; they reproduce an existing record and must not append a
 * duplicate on every replay (the guard the pixel-op capture always had). */
Test(flis_mask_route, replay_rerun_mask_op_captures_no_record) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	GPtrArray *before = nde_history_snapshot(NULL);
	guint n_before = before ? before->len : 0;
	if (before) g_ptr_array_unref(before);

	cr_assert(replay_reserve_slot(), "precondition: replay slot free");
	run_mask_op(gfit, 0);
	replay_release_slot();

	cr_assert_not_null(gfit->mask, "the op itself must still run under replay");
	GPtrArray *after = nde_history_snapshot(NULL);
	guint n_after = after ? after->len : 0;
	if (after) g_ptr_array_unref(after);
	cr_assert_eq(n_after, n_before,
	             "a replayed mask op must not append a duplicate record");
}

/* A resident script's open provenance scope must swallow only the script's
 * own ops.  This op runs from the test thread — a stand-in for a user-run
 * op, which is not a python-comm job — so it must capture per-op provenance
 * even while a scope is open. */
Test(flis_mask_route, resident_scope_does_not_swallow_user_mask_op) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	nde_script_scope_begin(NULL);
	run_mask_op(gfit, 0);
	nde_script_scope_end();

	const nde_record *rec = only_record();
	cr_assert_not_null(rec, "user op must leave provenance despite the scope");
	cr_assert_str_eq(rec->op_id, "mask.test_fill");
	cr_assert_eq(rec->target_item_id, flis_layer_pmask_id(base));
}

/* A command-context mask op saves no undo entry (the save is gated on
 * !args->command), so it must not re-tag whatever unrelated entry sits on
 * top of the undo stack — that would re-couple the previous operation's
 * snapshot to the mask record and desynchronise undo from the history. */
Test(flis_mask_route, command_mask_op_does_not_retag_foreign_undo_entry) {
	com.script = FALSE;    /* GUI-mode: the undo stack is live */
	com.headless = TRUE;
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	/* A previous operation's undo entry, coupled to its own record. */
	undo_save_state(gfit, "previous op");
	cr_assert_eq(g_list_length(com.undo_stack), 1);
	historic *top = (historic *)com.undo_stack->data;
	top->nde_record_id = 4242;

	run_mask_op(gfit, 0);   /* command context: saves no undo entry */

	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "a command mask op must not push an undo entry");
	cr_assert_eq(top->nde_record_id, 4242,
	             "…and must not re-tag the previous op's entry");
}

/* Guard-order regression: a Tier-C replay re-run while an unrelated resident
 * script's scope is open must not mark that scope dirty — or the innocent
 * script's exit would emit a spurious opaque barrier record for pixels the
 * replay legitimately reproduced.  The replay guard runs BEFORE the scope
 * check at every capture site. */
Test(flis_mask_route, replay_rerun_does_not_mark_stranger_scope) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	GPtrArray *before = nde_history_snapshot(NULL);
	guint n_before = before ? before->len : 0;
	if (before) g_ptr_array_unref(before);

	nde_script_scope_begin(NULL);           /* the resident stranger */
	cr_assert(replay_reserve_slot());
	run_mask_op(gfit, 0);                   /* the replay's mask command */
	replay_release_slot();
	nde_script_scope_end();                 /* stranger exits, wrote nothing */

	GPtrArray *after = nde_history_snapshot(NULL);
	guint n_after = after ? after->len : 0;
	if (after) g_ptr_array_unref(after);
	cr_assert_eq(n_after, n_before,
	             "neither the replayed op nor the innocent scope's exit "
	             "may append a record");
}
