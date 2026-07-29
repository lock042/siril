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
 * test_nde_composite — merge-down as a replayable composite node
 * (design note §3 and §9.0).
 *
 * The bug that motivated this: duplicate a layer, edit each copy, merge down,
 * and the pre-merge steps are visible but inert.  They were never lost — the
 * merged-away layer keeps its baseline and its records — what was missing was
 * the EDGE from those records to the merged result, and a merge that could
 * re-run.  These tests pin down both, and the honest refusal that must remain
 * for records written before either existed.
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/op_descriptors.h"
#include "core/op_descriptor.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_composite.h"
#include "core/nde_replay.h"
#include "core/nde_graph.h"
#include "core/masks.h"
#include "filters/asinh.h"

cominfo com;
fits *gfit;

TestSuite(nde_composite, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* ---- fixtures ---------------------------------------------------------- */

static asinh_params *asinh_beta(float beta) {
	asinh_params *p = calloc(1, sizeof(*p));
	p->beta = beta;
	p->clip_mode = RESCALE;
	return p;
}

/* Apply an op to the active layer through the real capture path, so the
 * record, its baseline and its checkpoint all come from production code. */
static int run_op_on_active(const op_descriptor *op, gpointer user) {
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->op = op;
	args->user = user;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;
	args->max_threads = 1;
	args->mem_ratio = -1.0f;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_image_worker(args));
	com.headless = prev;
	return rc;
}

static void select_layer(flis_layer_t *lay) {
	gint idx = flis_layer_get_index(lay);
	cr_assert(idx >= 0);
	uniq_set_active_layer(com.uniq, idx);
	gfit = flis_active_layer_fit();
}

/* Two RGB layers on a canvas, each carrying one asinh record of its own, then
 * merged.  Returns the surviving layer; @out_bottom_item / @out_top_item get
 * the two item ids (the top's outlives its layer). */
static flis_layer_t *two_edited_layers_merged(gint *out_bottom_item,
                                              gint *out_top_item,
                                              flis_blend_mode_t top_blend,
                                              float top_opacity) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "top");
	cr_assert_eq(flis_layer_set_blend_mode(top, top_blend), 0);
	cr_assert_eq(flis_layer_set_opacity(top, top_opacity), 0);

	select_layer(bottom);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(5.f)), 0);
	select_layer(top);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(30.f)), 0);

	if (out_bottom_item) *out_bottom_item = bottom->item_id;
	if (out_top_item)    *out_top_item    = top->item_id;

	cr_assert_eq(flis_merge_down_layer(top), 0);
	cr_assert_eq(flis_layer_count(), 1);
	flis_layer_t *merged = (flis_layer_t *)com.uniq->layers->data;
	gfit = merged->fit;
	return merged;
}

static void done(void) {
	nde_history_attach(NULL);
	gfit = NULL;
}

static const nde_record *find_composite_record(void) {
	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *found = NULL;
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (nde_composite_is_op(r->op_id))
			found = r;
	}
	/* The snapshot owns the records; keep it alive for the caller by leaking
	 * the array is not acceptable, so copy what the tests need instead. */
	static nde_record *copy = NULL;
	if (copy) { nde_record_free(copy); copy = NULL; }
	if (found) copy = nde_record_copy(found);
	if (snap) g_ptr_array_unref(snap);
	return copy;
}

/* ---- what the merge records ------------------------------------------- */

/* The merge is the moment at which both inputs still exist, so it is the only
 * moment at which their state can be written down.  Everything the replay
 * needs comes from here. */
Test(nde_composite, the_merge_records_both_inputs) {
	gint bottom_item = 0, top_item = 0;
	two_edited_layers_merged(&bottom_item, &top_item, FLIS_BLEND_NORMAL, 0.5f);

	const nde_record *merge = find_composite_record();
	cr_assert_not_null(merge, "the merge must be recorded");

	const nde_input_pin *base = nde_record_input(merge, NDE_COMPOSITE_ROLE_BASE);
	const nde_input_pin *over = nde_record_input(merge, NDE_COMPOSITE_ROLE_OVERLAY);
	cr_assert_not_null(base, "the surviving layer is an input, not just a target");
	cr_assert_not_null(over, "the merged-away layer is an input");
	cr_assert_eq(base->src_item_id, bottom_item);
	cr_assert_eq(over->src_item_id, top_item);
	cr_assert_neq(base->src_record_id, 0,
	              "each pin names the input's state AT THE MERGE, not its baseline");
	cr_assert_neq(over->src_record_id, 0);

	nde_composite_state *st = nde_composite_state_parse(merge->params);
	cr_assert_not_null(st, "the per-input compositing state must round-trip");
	cr_assert(st->raw_first, "merge-down paints its bottom input raw");
	cr_assert_eq(st->inputs->len, 2);
	const nde_composite_input *bs = &g_array_index(st->inputs, nde_composite_input, 0);
	const nde_composite_input *ov = &g_array_index(st->inputs, nde_composite_input, 1);
	cr_assert_eq(ov->item_id, top_item);
	cr_assert_str_eq(ov->name, "top");
	cr_assert_eq(ov->blend_mode, FLIS_BLEND_NORMAL);
	cr_assert_float_eq(ov->opacity, 0.5, 1e-6,
	                   "opacity dies with the layer, so it must be recorded");
	cr_assert(ov->visible);
	cr_assert_eq(bs->item_id, bottom_item);
	nde_composite_state_free(st);
	done();
}

/* A layer merged without ever being edited has no baseline yet, and after the
 * merge its pre-merge pixels exist nowhere else.  The capture site is the last
 * place that can take one. */
Test(nde_composite, merging_an_unedited_layer_still_leaves_it_rebuildable) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "top");
	gint bottom_item = bottom->item_id, top_item = top->item_id;
	cr_assert(!nde_checkpoint_baseline_exists(top_item),
	          "precondition: nothing has been recorded against either layer");

	cr_assert_eq(flis_merge_down_layer(top), 0);
	cr_assert(nde_checkpoint_baseline_exists(top_item),
	          "the merged-away layer's pixels must survive the merge");
	cr_assert(nde_checkpoint_baseline_exists(bottom_item),
	          "and so must the surviving layer's PRE-merge pixels");
	done();
}

/* ---- the chain ---------------------------------------------------------- */

/* The freeze this closes: before step 7 the merge was a wall, and every step
 * the surviving layer had taken before it was locked behind that wall. */
Test(nde_composite, the_merge_no_longer_freezes_the_layer_underneath) {
	gint bottom_item = 0;
	two_edited_layers_merged(&bottom_item, NULL, FLIS_BLEND_NORMAL, 1.0f);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable, "the merge must be replayable, not a blocker: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(chain->has_composite, "and the chain must know it consumes another item");
	cr_assert_eq(chain->records->len, 2,
	             "the layer's own step plus the merge");
	cr_assert_eq(chain->tail_start, 0,
	             "nothing is frozen: the whole chain is editable");
	nde_chain_free(chain);
	done();
}

/* Replaying the surviving layer's chain reproduces the merge exactly: the
 * composite is re-run from the two inputs rather than read back from a stored
 * copy of its own result. */
Test(nde_composite, replaying_through_the_merge_reproduces_it) {
	gint bottom_item = 0;
	flis_layer_t *merged = two_edited_layers_merged(&bottom_item, NULL,
	                                                FLIS_BLEND_MULTIPLY, 1.0f);
	fits expected = { 0 };
	copyfits(merged->fit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_eq(result->rx, expected.rx);
	cr_assert_eq(result->ry, expected.ry);
	cr_assert_eq(result->naxes[2], expected.naxes[2]);
	size_t n = (size_t)expected.rx * expected.ry * expected.naxes[2];
	float maxdev = 0.f;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(result->fdata[i] - expected.fdata[i]);
		if (d > maxdev) maxdev = d;
	}
	cr_assert(maxdev <= 1e-6f, "replayed composite deviates by %.3g", (double)maxdev);

	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	done();
}

/* ---- amending across the merge ----------------------------------------- */

/* Find the asinh record targeting @item. */
static gint64 record_for_item(gint item) {
	GPtrArray *snap = nde_history_snapshot(NULL);
	gint64 id = 0;
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (r->target_item_id == item && !g_strcmp0(r->op_id, "stretch.asinh"))
			id = r->record_id;
	}
	if (snap) g_ptr_array_unref(snap);
	return id;
}

static float first_pixel(void) {
	flis_layer_t *merged = (flis_layer_t *)com.uniq->layers->data;
	return merged->fit->fdata[0];
}

/* The surviving layer's own pre-merge step: it used to be locked behind the
 * merge, and now it is not. */
Test(nde_composite, amending_the_surviving_layers_pre_merge_step_takes_effect) {
	gint bottom_item = 0;
	two_edited_layers_merged(&bottom_item, NULL, FLIS_BLEND_NORMAL, 0.5f);
	gint64 rid = record_for_item(bottom_item);
	cr_assert_neq(rid, 0);
	float before = first_pixel();

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(rid, "beta=200.000000;offset=0.000000;human=0;clip_mode=0", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);

	cr_assert_neq(first_pixel(), before,
	              "the merged result must reflect the amended step");
	done();
}

/* The report that started this (del.png): the merged-away layer's steps were
 * listed, looked editable, and did nothing.  Amending one now re-runs the
 * merge, because the merge knows that layer is an input. */
Test(nde_composite, amending_the_merged_away_layers_step_takes_effect) {
	gint top_item = 0;
	two_edited_layers_merged(NULL, &top_item, FLIS_BLEND_NORMAL, 0.5f);
	gint64 rid = record_for_item(top_item);
	cr_assert_neq(rid, 0, "the merged-away layer's records are still in the log");
	cr_assert_null(flis_layer_get_by_id(top_item),
	               "and its layer really is gone");
	float before = first_pixel();

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(rid, "beta=1.000000;offset=0.000000;human=0;clip_mode=0", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);

	cr_assert_neq(first_pixel(), before,
	              "amending a step upstream of the merge must reach the merged image");
	done();
}

/* A pin names a record in the item's HISTORY, which is a wider set than its
 * pixel chain: setting the opacity just before merging makes the last record
 * against that layer a compositing-state one, which the chain does not contain.
 * Reported as "record N is no longer part of item M's history" — the merge
 * still worked, but every edit upstream of it announced a failure. */
Test(nde_composite, a_pin_may_name_a_record_the_pixel_chain_does_not_contain) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "top");
	select_layer(bottom);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(5.f)), 0);
	select_layer(top);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(30.f)), 0);
	gint top_item = top->item_id;
	/* the last word on this layer is not a pixel step (captured as the GUI
	 * and the flis_layer_opacity command do, outside the setter) */
	cr_assert_eq(flis_layer_set_opacity(top, 0.5f), 0);
	{
		GString *kv = nde_kv_start();
		nde_kv_add_float(kv, "opacity", 0.5f);
		cr_assert_neq(nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
		                                     top_item, nde_kv_end(kv), "Set opacity"), 0);
	}
	cr_assert_eq(flis_merge_down_layer(top), 0);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;

	const nde_record *merge = find_composite_record();
	const nde_input_pin *over = nde_record_input(merge, NDE_COMPOSITE_ROLE_OVERLAY);
	cr_assert_not_null(over);
	cr_assert_neq(over->src_record_id, record_for_item(top_item),
	              "precondition: the pin names the opacity record, not the asinh");

	float before = first_pixel();
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(record_for_item(top_item),
	                            "beta=1.000000;offset=0.000000;human=0;clip_mode=0", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);
	cr_assert_neq(first_pixel(), before);
	done();
}

/* ---- flatten: the same node, N inputs ----------------------------------- */

/* Three RGB layers, one step each, then flattened.  The middle one blends and
 * the top one is half-transparent, so the flattened pixels depend on state that
 * nothing but the record can carry once the layers are gone. */
static void three_layers_flattened(gint items[3]) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *l[3] = {
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "middle"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.1f, 0.1f, 0.1f), "top"),
	};
	com.uniq->canvas_w = 4;   /* the document paths set this; flis_layer_add does not */
	com.uniq->canvas_h = 4;
	cr_assert_eq(flis_layer_set_blend_mode(l[1], FLIS_BLEND_MULTIPLY), 0);
	cr_assert_eq(flis_layer_set_opacity(l[2], 0.5f), 0);
	const float betas[3] = { 5.f, 30.f, 12.f };
	for (int i = 0; i < 3; i++) {
		items[i] = l[i]->item_id;
		select_layer(l[i]);
		cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(betas[i])), 0);
	}
	cr_assert_eq(flis_flatten_all(), 0);
	cr_assert_eq(flis_layer_count(), 1);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;
}

/* Flatten consumes the whole stack, so the whole stack is recorded — and unlike
 * merge-down it blends its bottom input too, over the canvas background. */
Test(nde_composite, the_flatten_records_every_layer) {
	gint items[3] = { 0 };
	three_layers_flattened(items);

	const nde_record *flat = find_composite_record();
	cr_assert_not_null(flat);
	cr_assert_str_eq(flat->op_id, "document.flatten");

	nde_composite_state *st = nde_composite_state_parse(flat->params);
	cr_assert_not_null(st, "the state of every input must round-trip");
	cr_assert(!st->raw_first, "flatten blends its bottom input, it does not paint it raw");
	cr_assert_eq(st->inputs->len, 3);
	cr_assert_eq(st->canvas_w, 4, "the canvas it was composited against is recorded");
	cr_assert_eq(st->canvas_h, 4);
	for (guint i = 0; i < 3; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		cr_assert_eq(in->item_id, items[i]);
		const nde_input_pin *pin = nde_record_input_by_item(flat, items[i]);
		cr_assert_not_null(pin, "every input is pinned, not just the consumed ones");
		cr_assert_neq(pin->src_record_id, 0);
	}
	const nde_composite_input *mid = &g_array_index(st->inputs, nde_composite_input, 1);
	cr_assert_eq(mid->blend_mode, FLIS_BLEND_MULTIPLY);
	const nde_composite_input *top = &g_array_index(st->inputs, nde_composite_input, 2);
	cr_assert_float_eq(top->opacity, 0.5, 1e-6);
	nde_composite_state_free(st);
	done();
}

Test(nde_composite, the_flatten_no_longer_freezes_the_layer_it_survives_as) {
	gint items[3] = { 0 };
	three_layers_flattened(items);

	nde_chain *chain = nde_chain_build(items[0]);
	cr_assert(chain->replayable, "the flatten must be replayable, not a blocker: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(chain->has_composite);
	cr_assert_eq(chain->records->len, 2, "the layer's own step plus the flatten");
	cr_assert_eq(chain->tail_start, 0);
	nde_chain_free(chain);
	done();
}

Test(nde_composite, replaying_through_the_flatten_reproduces_it) {
	gint items[3] = { 0 };
	three_layers_flattened(items);
	flis_layer_t *flat = (flis_layer_t *)com.uniq->layers->data;
	fits expected = { 0 };
	copyfits(flat->fit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);

	nde_chain *chain = nde_chain_build(items[0]);
	cr_assert(chain->replayable);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_eq(result->rx, expected.rx);
	cr_assert_eq(result->naxes[2], expected.naxes[2]);
	size_t n = (size_t)expected.rx * expected.ry * expected.naxes[2];
	float maxdev = 0.f;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(result->fdata[i] - expected.fdata[i]);
		if (d > maxdev) maxdev = d;
	}
	cr_assert(maxdev <= 1e-6f, "replayed flatten deviates by %.3g", (double)maxdev);

	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	done();
}

/* The del.png bug, for flatten: a layer the flatten consumed is gone from the
 * document but its steps are live, and amending one re-runs the flatten. */
Test(nde_composite, amending_a_flattened_away_layers_step_takes_effect) {
	gint items[3] = { 0 };
	three_layers_flattened(items);
	cr_assert_null(flis_layer_get_by_id(items[1]), "the middle layer really is gone");
	gint64 rid = record_for_item(items[1]);
	cr_assert_neq(rid, 0);
	float before = first_pixel();

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(rid, "beta=1.000000;offset=0.000000;human=0;clip_mode=0", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);

	cr_assert_neq(first_pixel(), before,
	              "amending a step upstream of the flatten must reach the flattened image");
	done();
}

/* The reason the compositor needed an explicit context.  A group with a real
 * blend mode is pre-composited before it reaches the main composite, and after
 * the flatten the group has no members left — so a replay that consulted the
 * live document would quietly render something else. */
Test(nde_composite, a_flattened_group_replays_from_the_recorded_group_state) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *l[3] = {
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "g1"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.6f, 0.6f, 0.6f), "g2"),
	};
	flis_group_t *grp = flis_group_add("stack");
	cr_assert_not_null(grp);
	grp->blend_mode = FLIS_BLEND_MULTIPLY;   /* not PASS_THROUGH: pre-composited */
	grp->opacity    = 0.6f;
	grp->visible    = TRUE;
	cr_assert_eq(flis_layer_set_group(l[1], grp->item_id), 0);
	cr_assert_eq(flis_layer_set_group(l[2], grp->item_id), 0);
	gint base_item = l[0]->item_id;
	select_layer(l[0]);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(5.f)), 0);

	cr_assert_eq(flis_flatten_all(), 0);
	fits expected = { 0 };
	flis_layer_t *flat = (flis_layer_t *)com.uniq->layers->data;
	copyfits(flat->fit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);
	gfit = flat->fit;

	const nde_record *rec = find_composite_record();
	nde_composite_state *st = nde_composite_state_parse(rec->params);
	cr_assert_not_null(st);
	cr_assert_eq(st->groups->len, 1, "the group is recorded, not looked up later");
	cr_assert_eq(g_array_index(st->groups, nde_composite_group, 0).blend_mode,
	             FLIS_BLEND_MULTIPLY);
	nde_composite_state_free(st);

	nde_chain *chain = nde_chain_build(base_item);
	cr_assert(chain->replayable, "%s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	g_free(err);
	size_t n = (size_t)expected.rx * expected.ry * expected.naxes[2];
	float maxdev = 0.f;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(result->fdata[i] - expected.fdata[i]);
		if (d > maxdev) maxdev = d;
	}
	cr_assert(maxdev <= 1e-6f, "replayed group flatten deviates by %.3g", (double)maxdev);

	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	done();
}

/* Merging a group is a run of merge_downs, so the composites nest: the second
 * merge's overlay is an item whose own chain contains the first.  Resolving it
 * recurses one level and terminates — composites form a DAG — and an amend at
 * the top of the stack still reaches the bottom of it. */
Test(nde_composite, merges_of_merges_resolve_through_each_other) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *l[3] = {
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "middle"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.1f, 0.1f, 0.1f), "top"),
	};
	const float betas[3] = { 5.f, 30.f, 12.f };
	gint items[3];
	for (int i = 0; i < 3; i++) {
		items[i] = l[i]->item_id;
		select_layer(l[i]);
		cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(betas[i])), 0);
	}
	cr_assert_eq(flis_merge_down_layer(l[2]), 0);   /* top   -> middle */
	cr_assert_eq(flis_merge_down_layer(l[1]), 0);   /* middle -> bottom */
	cr_assert_eq(flis_layer_count(), 1);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;

	nde_chain *chain = nde_chain_build(items[0]);
	cr_assert(chain->replayable, "%s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);

	float before = first_pixel();
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(record_for_item(items[2]),
	                            "beta=1.000000;offset=0.000000;human=0;clip_mode=0", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);
	cr_assert_neq(first_pixel(), before,
	              "an amend two merges upstream must still reach the result");
	done();
}

/* ---- amending the composite itself -------------------------------------- */

/* @params with @key's value replaced.  The blob is one flat kv list, so this
 * is a substring swap rather than a re-encode. */
static gchar *replace_kv(const char *params, const char *key, const char *value) {
	gchar *needle = g_strdup_printf("%s=", key);
	const char *at = strstr(params, needle);
	cr_assert_not_null(at, "no key '%s' in: %s", key, params);
	const char *end = strchr(at, ';');
	cr_assert_not_null(end);
	gchar *head = g_strndup(params, (gsize)(at - params));
	gchar *out = g_strdup_printf("%s%s%s%s", head, needle, value, end);
	g_free(head);
	g_free(needle);
	return out;
}

/* The compositing state of a merged-away layer is recorded, which makes it the
 * one thing about that layer a user can still reconsider afterwards.  The
 * result must equal the merge that would have happened at that opacity —
 * "it changed" would not prove it changed to the right thing. */
Test(nde_composite, the_opacity_a_layer_was_merged_at_can_be_amended) {
	two_edited_layers_merged(NULL, NULL, FLIS_BLEND_NORMAL, 1.0f);
	fits native = { 0 };
	copyfits(((flis_layer_t *)com.uniq->layers->data)->fit, &native,
	         CP_DEEPCOPY | CP_ALLOC, -1);
	done();
	flis_free_layers(com.uniq);

	two_edited_layers_merged(NULL, NULL, FLIS_BLEND_NORMAL, 0.5f);
	const nde_record *merge = find_composite_record();
	gint64 rid = merge->record_id;
	gchar *amended = replace_kv(merge->params, "i1_opacity", "1.000000");
	cr_assert_neq(first_pixel(), native.fdata[0], "fixture: 0.5 and 1.0 must differ");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(rid, amended, &err), "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);
	g_free(amended);

	fits *now = ((flis_layer_t *)com.uniq->layers->data)->fit;
	size_t n = (size_t)native.rx * native.ry * native.naxes[2];
	float maxdev = 0.f;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(now->fdata[i] - native.fdata[i]);
		if (d > maxdev) maxdev = d;
	}
	cr_assert(maxdev <= 1e-6f,
	          "re-merging at opacity 1 must match having merged at 1 (deviation %.3g)",
	          (double)maxdev);
	clearfits(&native);
	done();
}

/* A group's properties modify every member at once, and after the flatten the
 * group is gone from the Layers panel along with its members — so the record is
 * the only place left that can reach them.  Its name is recorded for the same
 * reason each input's is: nothing else names it once it is gone. */
Test(nde_composite, the_opacity_a_group_was_flattened_at_can_be_amended) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *l[3] = {
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "g1"),
		flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.6f, 0.6f, 0.6f), "g2"),
	};
	flis_group_t *grp = flis_group_add("stack");
	cr_assert_not_null(grp);
	grp->blend_mode = FLIS_BLEND_MULTIPLY;   /* not PASS_THROUGH: pre-composited */
	grp->opacity    = 0.5f;
	grp->visible    = TRUE;
	cr_assert_eq(flis_layer_set_group(l[1], grp->item_id), 0);
	cr_assert_eq(flis_layer_set_group(l[2], grp->item_id), 0);
	select_layer(l[0]);

	cr_assert_eq(flis_flatten_all(), 0);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;

	const nde_record *rec = find_composite_record();
	nde_composite_state *st = nde_composite_state_parse(rec->params);
	cr_assert_not_null(st);
	cr_assert_eq(st->groups->len, 1);
	cr_assert_str_eq(g_array_index(st->groups, nde_composite_group, 0).name, "stack");
	nde_composite_state_free(st);

	float before = first_pixel();
	gint64 rid = rec->record_id;
	gchar *amended = replace_kv(rec->params, "g0_opacity", "1.000000");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(rid, amended, &err), "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);
	g_free(amended);

	cr_assert(fabsf(first_pixel() - before) > 1e-5f,
	          "the group's opacity must reach the flattened pixels (%.6f unchanged)",
	          (double)before);
	done();
}

/* An edit changes how a step composited, never what it composited: the graph
 * is not rewired by editing parameters (design note §5.4). */
Test(nde_composite, a_composite_cannot_be_rewired_by_an_amend) {
	gint top_item = 0;
	two_edited_layers_merged(NULL, &top_item, FLIS_BLEND_NORMAL, 0.5f);
	const nde_record *merge = find_composite_record();
	gint64 rid = merge->record_id;
	gchar *rewired = replace_kv(merge->params, "i1_item", "99");
	gchar *silly   = replace_kv(merge->params, "i1_opacity", "3.000000");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(!nde_amend_execute(rid, rewired, &err),
	          "pointing an input at another item must be refused");
	cr_assert_not_null(err);
	g_free(err);
	err = NULL;
	cr_assert(!nde_amend_execute(rid, silly, &err),
	          "an opacity outside 0..1 must be refused");
	cr_assert_not_null(err);
	unreserve_thread();
	g_free(err);
	g_free(rewired);
	g_free(silly);

	/* and the record is untouched, so the image still replays */
	nde_chain *chain = nde_chain_build(((flis_layer_t *)com.uniq->layers->data)->item_id);
	cr_assert(chain->replayable);
	nde_chain_free(chain);
	done();
}

/* ---- honest refusal ----------------------------------------------------- */

/* A merge recorded before step 7 carries neither pins nor per-input state.
 * There is nothing to re-run it from, and guessing at an opacity nobody wrote
 * down would be worse than saying so. */
Test(nde_composite, a_merge_without_recorded_inputs_still_blocks) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom");
	gint bottom_item = bottom->item_id;
	select_layer(bottom);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(5.f)), 0);

	/* Exactly what the pre-step-7 capture site wrote: the op id, the top
	 * layer's id and name, and nothing else. */
	nde_record *old = nde_record_new();
	old->op_id          = g_strdup("layer.merge_down");
	old->op_version     = 1;
	old->scope          = NDE_SCOPE_DOCUMENT;
	old->target_item_id = bottom_item;
	old->tier           = NDE_TIER_A;
	old->params         = g_strdup("top_item=2;top_name=top");
	old->summary        = g_strdup("Merge down");
	cr_assert_neq(nde_history_append(old), 0);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(!chain->replayable,
	          "a merge with no recorded inputs cannot be re-run");
	cr_assert(chain->reasons->len > 0, "and must say why");
	nde_chain_free(chain);
	done();
}

/* ---- masked inputs ------------------------------------------------------ */

/* Two layers, the top one half-masked, flattened.  Returns the surviving
 * layer's item id; @out_mask_item gets the LMASK item the mask was stored
 * under. */
static gint masked_flatten(gint *out_mask_item) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.5f, 0.5f, 0.5f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(4, 4, 0.25f, 0.25f, 0.25f), "top");
	cr_assert_eq(flis_layer_set_lmask(top, flis_test_make_const_lmask(4, 4, 8, 0.5)), 0);
	top->lmask_active = TRUE;
	gint bottom_item = bottom->item_id;
	select_layer(bottom);
	cr_assert_eq(run_op_on_active(&op_desc_asinh, asinh_beta(5.f)), 0);

	cr_assert_eq(flis_flatten_all(), 0);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	if (out_mask_item) {
		const nde_record *rec = find_composite_record();
		nde_composite_state *st = nde_composite_state_parse(rec->params);
		cr_assert_not_null(st);
		*out_mask_item = g_array_index(st->inputs, nde_composite_input, 1).mask_item_id;
		cr_assert_neq(*out_mask_item, 0, "the masked input records its LMASK item");
		nde_composite_state_free(st);
	}
	return bottom_item;
}

/* A layer mask can be painted or loaded from a file, so it may have no chain to
 * replay at all.  It is kept as a stored copy instead — and that is enough to
 * make a composite that applied one reproducible. */
Test(nde_composite, a_flatten_over_a_masked_layer_replays) {
	gint mask_item = 0;
	gint bottom_item = masked_flatten(&mask_item);
	fits expected = { 0 };
	copyfits(((flis_layer_t *)com.uniq->layers->data)->fit, &expected,
	         CP_DEEPCOPY | CP_ALLOC, -1);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable, "%s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	g_free(err);
	size_t n = (size_t)expected.rx * expected.ry * expected.naxes[2];
	float maxdev = 0.f;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(result->fdata[i] - expected.fdata[i]);
		if (d > maxdev) maxdev = d;
	}
	/* If the mask were dropped the top layer would land at full strength, so
	 * this comparison is what proves it was applied. */
	cr_assert(maxdev <= 1e-6f, "replayed masked flatten deviates by %.3g", (double)maxdev);

	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	done();
}

/* The stored copy lives in a cache under a budget.  If it goes, the composite
 * says so rather than compositing without the mask and calling the result a
 * reproduction. */
Test(nde_composite, a_masked_composite_blocks_once_its_mask_is_gone) {
	gint mask_item = 0;
	gint bottom_item = masked_flatten(&mask_item);
	nde_chain *ok = nde_chain_build(bottom_item);
	cr_assert(ok->replayable, "precondition: it replays while the mask is stored");
	nde_chain_free(ok);

	nde_checkpoint_drop(mask_item);
	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(!chain->replayable, "with the mask gone it cannot be reproduced");
	cr_assert(chain->reasons->len > 0, "and must say why");
	nde_chain_free(chain);
	done();
}

/* The first merge-down format, before flatten joined this node and the params
 * became indexed.  Documents saved with it must keep replaying rather than
 * silently reverting to blockers, so the parser still accepts it. */
Test(nde_composite, the_first_merge_format_still_decodes) {
	nde_composite_state *st = nde_composite_state_parse(
	    "base_item=1;base_x=0;base_y=0;base_tinted=0;base_tint_r=1.000000;"
	    "base_tint_g=1.000000;base_tint_b=1.000000;top_item=2;top_name=top;"
	    "top_blend=3;top_opacity=0.500000;top_x=1;top_y=2;top_visible=1;"
	    "top_tinted=0;top_tint_r=1.000000;top_tint_g=1.000000;"
	    "top_tint_b=1.000000;top_masked=0");
	cr_assert_not_null(st, "a merge recorded by the shipped format must still parse");
	cr_assert(st->raw_first, "it is a merge-down, so its base is painted raw");
	cr_assert_eq(st->inputs->len, 2);
	const nde_composite_input *ov = &g_array_index(st->inputs, nde_composite_input, 1);
	cr_assert_eq(ov->item_id, 2);
	cr_assert_str_eq(ov->name, "top");
	cr_assert_eq(ov->blend_mode, 3);
	cr_assert_float_eq(ov->opacity, 0.5, 1e-6);
	cr_assert_eq(ov->position_x, 1);
	cr_assert_eq(ov->position_y, 2);
	cr_assert(!ov->was_masked);
	nde_composite_state_free(st);
}

/* Reported (weird.png): the mask item a flatten consumed appeared as "Layer 3,
 * no longer in the document" — a layer the user never had.  Item ids are one
 * sequence across layers, masks and groups, so the fallback label has to say
 * which KIND it lost. */
Test(nde_composite, a_consumed_layer_mask_is_labelled_as_one) {
	gint mask_item = 0;
	masked_flatten(&mask_item);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *n = nde_graph_node_for(g, mask_item);
	cr_assert_not_null(n, "the mask a composite kept is still a node");
	cr_assert_null(strstr(n->label, "Layer 3"),
	               "it is not a layer: got '%s'", n->label);
	cr_assert_not_null(strstr(n->label, "mask"),
	                   "and it must say so: got '%s'", n->label);
	cr_assert_not_null(strstr(n->label, "top"),
	                   "named after the layer it masked: got '%s'", n->label);
	nde_graph_free(g);
	done();
}

/* ---- persistence -------------------------------------------------------- */

/* A merge survives a save/load: the pins and the per-input state ride in the
 * history table, and the merged-away layer's baseline has to be written even
 * though the document no longer lists the layer. */
Test(nde_composite, a_merge_is_still_replayable_after_a_round_trip) {
	gint bottom_item = 0, top_item = 0;
	two_edited_layers_merged(&bottom_item, &top_item, FLIS_BLEND_NORMAL, 0.5f);

	gchar *dir = g_dir_make_tmp("nde-composite-XXXXXX", NULL);
	cr_assert_not_null(dir);
	gchar *path = g_build_filename(dir, "merged.flis", NULL);
	cr_assert_eq(save_flis(path), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(path), 0);

	cr_assert(nde_checkpoint_baseline_exists(top_item),
	          "the merged-away layer's baseline must be written and read back");
	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable, "the reloaded merge must still replay: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(chain->has_composite);
	nde_chain_free(chain);

	g_unlink(path);
	g_free(path);
	g_rmdir(dir);
	g_free(dir);
	done();
}

/* Route a mask op onto @target's LMASK slot, as the -layermask= command path
 * does, so the mask gets a chain of its own. */
static int run_mask_op_on(flis_layer_t *target, const op_descriptor *op,
                          gpointer user_heap) {
	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit             = target->fit;
	args->op              = op;
	args->user            = user_heap;
	args->command         = TRUE;
	args->mask_creation   = (op->flags & OP_MASK_FROM_IMAGE) != 0;
	args->target_layer_id = target->item_id;
	args->max_threads     = 1;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_mask_worker(args));
	com.headless = prev;
	return rc;
}

static mask_from_channel_data *from_channel_8bit(int chan) {
	mask_from_channel_data *d = calloc(1, sizeof(*d));
	d->channel = chan;
	d->bitpix = 8;
	return d;
}

/* A stored copy is not a dead end when the mask DOES have a chain: the mask
 * cascade already refreshes stored copies, and it now matches them by source
 * rather than by role, so a composite's per-input mask pin refreshes with the
 * rest.  Amending the step that built the mask must therefore reach the
 * flattened image — through a layer that no longer exists. */
Test(nde_composite, amending_a_mask_a_composite_used_reaches_the_result) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *bottom = flis_test_add_layer(
	    flis_test_make_rgb_fits(8, 8, 0.5f, 0.5f, 0.5f), "bottom");
	/* channels deliberately unequal: the mask must differ per channel, or
	 * amending which one it came from would prove nothing. */
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(8, 8, 0.2f, 0.5f, 0.9f), "top");
	gint bottom_item = bottom->item_id, top_item = top->item_id;
	select_layer(top);
	cr_assert_eq(run_mask_op_on(top, &op_desc_mask_from_channel, from_channel_8bit(0)), 0);
	cr_assert_not_null(top->lmask, "the op must have routed to the layer mask");
	top->lmask_active = TRUE;

	gint64 mask_rec = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (!g_strcmp0(r->op_id, "mask.from_channel"))
			mask_rec = r->record_id;
	}
	g_ptr_array_unref(snap);
	cr_assert_neq(mask_rec, 0);

	cr_assert_eq(flis_flatten_all(), 0);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	float before = first_pixel();
	cr_assert_null(flis_layer_get_by_id(top_item));

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(mask_rec, "channel=2;autostretch=0;invert=0;bitpix=8", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable, "the stored mask must have been refreshed, not dropped: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);
	cr_assert_neq(first_pixel(), before,
	              "a mask taken from a different channel must change the flattened image");
	done();
}

/* Reported (weird.png): create a layer mask, flatten, then try to delete the
 * mask — "failed to load the baseline checkpoint", about an item the user had
 * never seen.  Deleting the only step of a MASK item took the image path,
 * because a trial chain with that step excluded is empty and the mask test
 * reads the first member's op.  Deleting it is also the only way to replay a
 * composite WITHOUT a mask that applied at the time, so it has to work rather
 * than merely fail better. */
Test(nde_composite, deleting_the_step_that_made_a_mask_unmasks_the_composite) {
	com.pref.nde_cache_mb = 256;
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.5f, 0.5f, 0.5f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_rgb_fits(8, 8, 0.2f, 0.5f, 0.9f), "top");
	select_layer(top);
	cr_assert_eq(run_mask_op_on(top, &op_desc_mask_from_channel, from_channel_8bit(0)), 0);
	cr_assert_not_null(top->lmask);
	top->lmask_active = TRUE;
	gint64 mask_rec = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (!g_strcmp0(r->op_id, "mask.from_channel"))
			mask_rec = r->record_id;
	}
	g_ptr_array_unref(snap);
	cr_assert_neq(mask_rec, 0);
	cr_assert_eq(flis_flatten_all(), 0);
	gfit = ((flis_layer_t *)com.uniq->layers->data)->fit;
	gint mask_item = 0;
	{
		const nde_record *rec = find_composite_record();
		nde_composite_state *st = nde_composite_state_parse(rec->params);
		cr_assert_not_null(st);
		mask_item = g_array_index(st->inputs, nde_composite_input, 1).mask_item_id;
		cr_assert_neq(mask_item, 0);
		nde_composite_state_free(st);
	}

	float masked_result = first_pixel();
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_delete_execute(mask_rec, &err),
	          "deleting the step that built the mask must work: %s", err ? err : "?");
	unreserve_thread();
	g_free(err);

	/* The composite kept a copy of the mask, and the point of deleting the
	 * step is to composite without it — so the flattened image must change,
	 * and the record must say it no longer applies one. */
	cr_assert_neq(first_pixel(), masked_result,
	              "the flatten must be re-run without the mask");
	const nde_record *flat = find_composite_record();
	nde_composite_state *st = nde_composite_state_parse(flat->params);
	cr_assert_not_null(st);
	const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, 1);
	cr_assert(!in->was_masked, "and must record that no mask applies");
	cr_assert_eq(in->mask_item_id, 0);
	nde_composite_state_free(st);
	cr_assert_null(nde_record_input_by_item(flat, mask_item),
	               "the pin to the now-historyless mask must go too");

	nde_chain *chain = nde_chain_build(
	    ((flis_layer_t *)com.uniq->layers->data)->item_id);
	cr_assert(chain->replayable, "and the flatten must still replay: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);
	done();
}

/* The stored mask is not the layer's any more — the layer is gone — so it has
 * to be written under its own item id and read back with it. */
Test(nde_composite, a_masked_composite_is_still_replayable_after_a_round_trip) {
	gint mask_item = 0;
	gint bottom_item = masked_flatten(&mask_item);

	gchar *dir = g_dir_make_tmp("nde-composite-XXXXXX", NULL);
	cr_assert_not_null(dir);
	gchar *path = g_build_filename(dir, "masked.flis", NULL);
	cr_assert_eq(save_flis(path), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	nde_checkpoint_purge();
	cr_assert_eq(load_flis(path), 0);

	nde_chain *chain = nde_chain_build(bottom_item);
	cr_assert(chain->replayable, "the reloaded masked flatten must still replay: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);

	g_unlink(path);
	g_free(path);
	g_rmdir(dir);
	g_free(dir);
	done();
}

/* ---- the graph ---------------------------------------------------------- */

/* The pins are what turn two disconnected lists into one graph: the
 * merged-away layer is still a node, and now there is an edge from it to the
 * layer that consumed it. */
Test(nde_composite, the_merged_away_layer_is_joined_to_its_consumer) {
	gint bottom_item = 0, top_item = 0;
	two_edited_layers_merged(&bottom_item, &top_item, FLIS_BLEND_NORMAL, 1.0f);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *top_node = nde_graph_node_for(g, top_item);
	cr_assert_not_null(top_node, "the merged-away layer is still a node");
	cr_assert(top_node->orphan, "the Layers panel no longer lists it");
	cr_assert_not_null(strstr(top_node->label, "top"),
	                   "and it is named, not numbered: got '%s'", top_node->label);

	GPtrArray *out = nde_graph_consumers_of(g, top_item);
	cr_assert_eq(out->len, 1, "exactly one node consumes it");
	const nde_graph_edge *e = g_ptr_array_index(out, 0);
	cr_assert_eq(e->dst_item_id, bottom_item);
	cr_assert_str_eq(e->role, NDE_COMPOSITE_ROLE_OVERLAY);
	g_ptr_array_unref(out);
	nde_graph_free(g);
	done();
}
