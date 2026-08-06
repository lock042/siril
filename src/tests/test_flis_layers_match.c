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
#include "core/processing_thread.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "core/nde_graph.h"
#include "core/nde_joint.h"
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

/* ---- NDE provenance (ONE joint "flis.layers_match" record) --------------- */

/* The joint record in the live log, or 0. */
static gint64 find_joint_record(void) {
	gint64 id = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "flis.layers_match"))
			id = rec->record_id;
	}
	if (snap)
		g_ptr_array_unref(snap);
	return id;
}

/* Live tint + its compositing-state record, the way the GUI records one. */
static void set_and_record_tint(flis_layer_t *lay, double r, double g, double b) {
	flis_layer_set_tint(lay, r, g, b);
	gchar *params = g_strdup_printf("tinted=1;r=%.17g;g=%.17g;b=%.17g", r, g, b);
	nde_capture_structural("layer.set_tint", NDE_SCOPE_LAYER, lay->item_id,
	                       params, "Layer tint");
}

/* Three tinted mono layers whose tints are ALSO in the log, so the joint
 * replay's positional fold sees them.  Factors are all non-identity for
 * medians 0.2 / 0.25 / 0.4: old_bg = 0.28333, s = (1.41667, 1.13333,
 * 0.70833), every post-match median 0.28333. */
static void make_recorded_triple(void) {
	flis_layer_t *r = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.2f), "R");
	set_and_record_tint(r, 1.0, 0.0, 0.0);
	flis_layer_t *g = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.25f), "G");
	set_and_record_tint(g, 0.0, 1.0, 0.0);
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.4f), "B");
	set_and_record_tint(b, 0.0, 0.0, 1.0);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
}

static int run_layers_match(void) {
	word[0] = "flis_layers_match";
	word[1] = NULL;
	return process_flis_layers_match(1);
}

static float layer_value(guint nth) {
	flis_layer_t *lay = g_slist_nth_data(com.uniq->layers, nth);
	return lay->fit->fdata[0];
}

/* The whole operation is ONE record, not one per layer: its meaning is the
 * intent (neutralise the composite background given these layers and their
 * tints), and only a single record can recompute that at replay.  It targets
 * the first participant as its anchor and pins every participant, in order,
 * under the roles the graph reads. */
Test(flis_layers_match, captures_one_joint_record) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap, "no NDE history after layers match");
	guint n_joint = 0, n_scale = 0;
	for (guint i = 0; i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "flis.layer_scale"))
			n_scale++;
		if (g_strcmp0(rec->op_id, "flis.layers_match"))
			continue;
		n_joint++;
		cr_assert_eq(rec->tier, NDE_TIER_A);
		cr_assert_eq(rec->scope, NDE_SCOPE_LAYER);
		cr_assert_eq(rec->target_item_id, 1,
		             "anchor must be the first participant");
		cr_assert_not_null(rec->params);
		cr_assert(strstr(rec->params, "n=3") != NULL,
		          "params '%s' missing participant count", rec->params);

		cr_assert_not_null(rec->inputs);
		cr_assert_eq(rec->inputs->len, 3, "every participant must be pinned");
		for (guint k = 0; k < rec->inputs->len; k++) {
			const nde_input_pin *pin = g_ptr_array_index(rec->inputs, k);
			char role[16];
			g_snprintf(role, sizeof(role), "in%u", k);
			cr_assert_str_eq(pin->role, role);
			cr_assert_eq(pin->src_item_id, (gint)k + 1);
		}
	}
	cr_assert_eq(n_joint, 1, "expected ONE joint record, got %u", n_joint);
	cr_assert_eq(n_scale, 0,
	             "layers match must not capture per-layer records any more");
	g_ptr_array_unref(snap);

	/* And it did the arithmetic: every median lands on the old background. */
	cr_assert_float_eq(layer_value(0), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.28333f, 1e-4);
}

/* One record at one log position, shared by N chains: the chain-membership
 * extension is what makes an undo across it one step and what lets any
 * participant's chain replay through it. */
Test(flis_layers_match, joint_record_is_member_of_every_participant_chain) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);

	for (gint item = 1; item <= 3; item++) {
		nde_chain *chain = nde_chain_build(item);
		gboolean member = FALSE;
		for (guint i = 0; i < chain->records->len; i++) {
			const nde_record *rec = g_ptr_array_index(chain->records, i);
			member = member || rec->record_id == jid;
		}
		cr_assert(member, "joint record missing from item %d's chain", item);
		cr_assert(chain->replayable, "item %d's chain must replay", item);
		nde_chain_free(chain);
	}
}

/* Deleting the one record undoes the whole operation: each participant
 * replays its own chain without it and returns to its original pixels. */
Test(flis_layers_match, delete_reverts_every_participant) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);

	gchar *err = NULL;
	cr_assert(nde_delete_execute(jid, &err), "delete failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_float_eq(layer_value(0), 0.2f, 1e-5);
	cr_assert_float_eq(layer_value(1), 0.25f, 1e-5);
	cr_assert_float_eq(layer_value(2), 0.4f, 1e-5);
	cr_assert_eq(find_joint_record(), 0);
}

/* The point of recording the intent rather than the factors: re-assign one
 * layer's tint UPSTREAM of the match and the analysis re-runs against the
 * new tint matrix, so all three participants take fresh factors.  B's tint
 * halved to (0,0,0.5) drops the composite background to 0.21667 and leaves
 * B itself at twice that. */
Test(flis_layers_match, tint_amend_recomputes_all_participants) {
	make_recorded_triple();

	gint64 b_tint_rec = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "layer.set_tint") && rec->target_item_id == 3)
			b_tint_rec = rec->record_id;
	}
	if (snap)
		g_ptr_array_unref(snap);
	cr_assert(b_tint_rec > 0);

	cr_assert_eq(run_layers_match(), CMD_OK);
	cr_assert_float_eq(layer_value(2), 0.28333f, 1e-4);

	gchar *err = NULL;
	gboolean ok = nde_amend_execute(b_tint_rec, "tinted=1;r=0;g=0;b=0.5", &err);
	cr_assert(ok, "tint amend failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_float_eq(layer_value(0), 0.21667f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.21667f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.43333f, 1e-4);
}

/* The mirror image: a tint recorded AFTER the match is not upstream of it,
 * so the positional fold must not see it and nothing rescales.  (Amending
 * the joint record itself with its own params re-runs the analysis, which
 * is what makes this a real test rather than a no-op.) */
Test(flis_layers_match, tint_appended_after_match_does_not_rescale) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);

	flis_layer_t *b = g_slist_nth_data(com.uniq->layers, 2);
	set_and_record_tint(b, 0.0, 0.0, 0.5);

	gchar *params = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid)
			params = g_strdup(rec->params);
	}
	if (snap)
		g_ptr_array_unref(snap);
	cr_assert_not_null(params);

	gchar *err = NULL;
	gboolean ok = nde_amend_execute(jid, params, &err);
	cr_assert(ok, "joint self-amend failed: %s", err ? err : "?");
	g_free(err);
	g_free(params);

	cr_assert_float_eq(layer_value(0), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.28333f, 1e-4);
}

/* The factors are cached per (record, history generation), so one edit costs
 * ONE analysis however many participants replay through the record. */
Test(flis_layers_match, one_analysis_per_edit) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);

	gint64 b_tint_rec = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "layer.set_tint") && rec->target_item_id == 3)
			b_tint_rec = rec->record_id;
	}
	if (snap)
		g_ptr_array_unref(snap);

	guint before = nde_joint_analysis_runs();
	gchar *err = NULL;
	cr_assert(nde_amend_execute(b_tint_rec, "tinted=1;r=0;g=0;b=0.5", &err),
	          "%s", err ? err : "?");
	g_free(err);
	guint delta = nde_joint_analysis_runs() - before;
	cr_assert_eq(delta, 1, "3 participants replayed but the analysis ran %u times",
	             delta);
}

/* Editing the joint record recomputes every participant, so the freeze rule
 * has to hold on each of them: an opaque step recorded on a SIBLING's chain
 * locks the record just as one on the anchor's would. */
Test(flis_layers_match, opaque_step_on_sibling_locks_joint_record) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();

	flis_layer_t *g = g_slist_nth_data(com.uniq->layers, 1);
	cr_assert(nde_capture_opaque("opaque.test", NDE_SCOPE_LAYER, 2,
	                             "Opaque step", g->fit) > 0);

	gchar *params = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid)
			params = g_strdup(rec->params);
	}
	if (snap)
		g_ptr_array_unref(snap);

	gchar *err = NULL;
	cr_assert(!nde_amend_execute(jid, params, &err),
	          "amend must refuse across the sibling's opaque step");
	cr_assert_not_null(err);
	cr_assert(strstr(err, "locked by a later opaque step") != NULL,
	          "unexpected error: %s", err);
	g_free(err);
	g_free(params);
}

/* A change strictly upstream of the match on ONE participant re-derives the
 * factors for ALL of them.  R starts at 0.4 and is halved to 0.2 before the
 * match (old_bg 0.28333); amending that step to scale=1 puts R back at 0.4,
 * and the re-run analysis lands every layer on 0.35. */
Test(flis_layers_match, upstream_amend_recomputes_siblings) {
	flis_layer_t *r = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.4f), "R");
	set_and_record_tint(r, 1.0, 0.0, 0.0);
	flis_layer_t *g = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.25f), "G");
	set_and_record_tint(g, 0.0, 1.0, 0.0);
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.4f), "B");
	set_and_record_tint(b, 0.0, 0.0, 1.0);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct flis_layer_scale_data *p = new_flis_layer_scale_data();
	nde_checkpoint_baseline_ensure(r->fit, r->item_id);
	flis_affine_layer_pixels(r->fit, 0.5, 0.0);
	p->scale = 0.5;
	gint64 scale_rec = nde_capture_from_descriptor_for_item(
	    &op_desc_flis_layer_scale, p, "Halve R", r->fit, r->item_id);
	free_flis_layer_scale_data(p);
	cr_assert(scale_rec > 0);

	cr_assert_eq(run_layers_match(), CMD_OK);
	cr_assert_float_eq(layer_value(0), 0.28333f, 1e-4);

	gchar *err = NULL;
	cr_assert(nde_amend_execute(scale_rec, "scale=1", &err), "%s", err ? err : "?");
	g_free(err);

	cr_assert_float_eq(layer_value(0), 0.35f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.35f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.35f, 1e-4);
}

/* In the graph the joint record is not a step on any one column: it is a
 * BAND drawn alone across its participants, with no edges of its own, and
 * every participant's later steps continue in a segment below it. */
Test(flis_layers_match, graph_shows_one_spanning_node) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	gint pid = nde_graph_joint_pseudo_item(jid);

	/* A step recorded AFTER the match on one participant. */
	flis_layer_t *g_lay = g_slist_nth_data(com.uniq->layers, 1);
	flis_affine_layer_pixels(g_lay->fit, 0.5, 0.0);
	struct flis_layer_scale_data *p = new_flis_layer_scale_data();
	p->scale = 0.5;
	gint64 after = nde_capture_from_descriptor_for_item(
	    &op_desc_flis_layer_scale, p, "Halve G", g_lay->fit, g_lay->item_id);
	free_flis_layer_scale_data(p);
	cr_assert(after > 0);

	nde_graph *graph = nde_graph_build();
	const nde_graph_node *jn = nde_graph_node_for(graph, pid);
	cr_assert_not_null(jn, "no dedicated node for the joint record");
	cr_assert_eq(jn->kind, NDE_NODE_JOINT);
	cr_assert(jn->spanning, "the joint record must be a spanning band");
	cr_assert_eq(jn->record_ids->len, 1);

	for (guint i = 0; i < graph->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(graph->edges, i);
		cr_assert(e->src_item_id != pid && e->dst_item_id != pid,
		          "a joint band must draw no edges");
	}

	/* The participants sit above it and no longer list its record, nor the
	 * step that came after it. */
	for (gint item = 1; item <= 3; item++) {
		const nde_graph_node *n = nde_graph_node_for(graph, item);
		cr_assert_not_null(n);
		cr_assert(n->level < jn->level, "participants must sit above the band");
		for (guint i = 0; i < n->record_ids->len; i++) {
			gint64 rid = *(gint64 *)g_ptr_array_index(n->record_ids, i);
			cr_assert_neq(rid, jid, "joint record inside item %d's node", item);
			cr_assert_neq(rid, after,
			              "a post-match step must not draw above the band");
		}
	}

	/* It draws in a continuation segment below the band instead. */
	gint seg = nde_graph_route_for_record(graph, after);
	cr_assert(seg != 0 && seg != g_lay->item_id,
	          "the post-match step must route to a continuation segment");
	const nde_graph_node *sn = nde_graph_node_for(graph, seg);
	cr_assert_not_null(sn);
	cr_assert_eq(sn->real_item, g_lay->item_id);
	cr_assert(sn->level > jn->level, "the segment must sit below the band");
	cr_assert_str_eq(sn->label, "G", "the segment keeps its layer's name");
	cr_assert_eq(nde_graph_route_for_record(graph, jid), pid);
	nde_graph_free(graph);
}

/* A subset amend rewrites the participant list, so the list has to be
 * validated against the document before anything is replayed: a layer the
 * document does not have must refuse, not fail deep inside the replay. */
Test(flis_layers_match, subset_amend_refuses_unknown_layer) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();

	gchar *params = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid)
			params = g_strdup(rec->params);
	}
	if (snap)
		g_ptr_array_unref(snap);
	cr_assert_not_null(params);

	GString *edited = g_string_new(params);
	g_string_replace(edited, "i2_item=3", "i2_item=7", 0);

	gchar *err = NULL;
	cr_assert(!nde_amend_execute(jid, edited->str, &err),
	          "a participant the document lacks must refuse");
	cr_assert_not_null(err);
	cr_assert(strstr(err, "not in the document") != NULL,
	          "unexpected error: %s", err);
	g_free(err);
	g_string_free(edited, TRUE);
	g_free(params);
}

/* Insertion BEFORE the joint record: the edit-at session must show the
 * pre-match state (the insertion is upstream of the match, so the match's
 * effect must be peeled back for it), and finishing re-runs the analysis
 * against the inserted step — halving R before the match lands all three
 * layers on 0.25 instead of 0.28333. */
Test(flis_layers_match, edit_at_insert_before_joint_rescales_all) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert_float_eq(layer_value(1), 0.28333f, 1e-4);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	gboolean opened = nde_edit_at_begin_execute(jid, &err);
	unreserve_thread();
	cr_assert(opened, "insertion before the joint record must open: %s",
	          err ? err : "?");
	g_free(err);
	cr_assert(nde_edit_at_active());

	err = NULL;
	flis_layer_t *r = g_slist_nth_data(com.uniq->layers, 0);
	cr_assert_float_eq(r->fit->fdata[0], 0.2f, 1e-5,
	                   "the insertion must show the pre-match state");

	flis_affine_layer_pixels(r->fit, 0.5, 0.0);
	struct flis_layer_scale_data *p = new_flis_layer_scale_data();
	p->scale = 0.5;
	gint64 ins = nde_capture_from_descriptor_for_item(
	    &op_desc_flis_layer_scale, p, "Halve R", r->fit, r->item_id);
	free_flis_layer_scale_data(p);
	cr_assert(ins > 0);

	cr_assert(reserve_thread());
	gboolean done = nde_edit_at_end_execute(TRUE, &err);
	unreserve_thread();
	cr_assert(done, "finish failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_float_eq(layer_value(0), 0.25f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.25f, 1e-4,
	                   "siblings must rescale when the insertion changed the solve");
	cr_assert_float_eq(layer_value(2), 0.25f, 1e-4);
}

/* The joint record's params, deserialized — the subset-amend tests edit the
 * participant list through the op's own structs rather than by hand. */
static struct nde_joint_layers_match_data *joint_record_data(gint64 jid) {
	struct nde_joint_layers_match_data *out = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid)
			out = op_desc_flis_layers_match.deserialize(rec->params,
			                                            op_desc_flis_layers_match.version);
	}
	if (snap)
		g_ptr_array_unref(snap);
	return out;
}

static void copy_participant(nde_joint_participant *dst,
                             const nde_joint_participant *src) {
	*dst = *src;
	dst->name = g_strdup(src->name);
}

/* How many input pins the joint record carries — the participant list and
 * the pins have to stay in step across a subset amend. */
static guint joint_record_pin_count(gint64 jid) {
	guint n = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid && rec->inputs)
			n = rec->inputs->len;
	}
	if (snap)
		g_ptr_array_unref(snap);
	return n;
}

/* Dropping a participant re-solves for the ones that remain (R and G alone
 * neutralise at 0.15) AND reverts the removed layer to its own pixels — it
 * is no longer a member of the record's chain, so the record's effect on it
 * has to be undone, not merely stopped. */
Test(flis_layers_match, subset_amend_removes_participant) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();

	struct nde_joint_layers_match_data *old = joint_record_data(jid);
	cr_assert_not_null(old);
	struct nde_joint_layers_match_data *edited = nde_joint_layers_match_data_new(2);
	for (guint k = 0; k < 2; k++)
		copy_participant(&edited->parts[k], &old->parts[k]);
	gchar *blob = op_desc_flis_layers_match.serialize(edited);
	nde_joint_layers_match_data_free(edited);
	nde_joint_layers_match_data_free(old);

	gchar *err = NULL;
	cr_assert(nde_amend_execute(jid, blob, &err), "%s", err ? err : "?");
	g_free(err);
	g_free(blob);

	cr_assert_float_eq(layer_value(0), 0.15f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.15f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.4f, 1e-4,
	                   "the removed layer must revert to its original pixels");

	cr_assert_eq(joint_record_pin_count(jid), 2, "pins must follow the new list");
}

/* And the other direction: a layer that was never a participant joins, takes
 * a baseline from its current pixels, and is scaled by the new solve — all
 * three land on 0.28333 as if the match had been run on them from the
 * start. */
Test(flis_layers_match, subset_amend_adds_participant) {
	make_recorded_triple();
	word[0] = "flis_layers_match";
	word[1] = "-subset=1,2";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);
	cr_assert_float_eq(layer_value(0), 0.15f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.4f, 1e-5, "B must start untouched");
	cr_assert_eq(joint_record_pin_count(jid), 2);

	struct nde_joint_layers_match_data *old = joint_record_data(jid);
	cr_assert_not_null(old);
	cr_assert_eq(old->n, 2);
	struct nde_joint_layers_match_data *edited = nde_joint_layers_match_data_new(3);
	for (guint k = 0; k < 2; k++)
		copy_participant(&edited->parts[k], &old->parts[k]);
	edited->parts[2].item_id = 3;
	edited->parts[2].name = g_strdup("B");
	edited->parts[2].tinted = TRUE;
	edited->parts[2].tint[0] = 0.0;
	edited->parts[2].tint[1] = 0.0;
	edited->parts[2].tint[2] = 1.0;
	edited->parts[2].diag_scale = 1.0;
	edited->parts[2].diag_offset = 0.0;
	gchar *blob = op_desc_flis_layers_match.serialize(edited);
	nde_joint_layers_match_data_free(edited);
	nde_joint_layers_match_data_free(old);

	gchar *err = NULL;
	cr_assert(nde_amend_execute(jid, blob, &err), "%s", err ? err : "?");
	g_free(err);
	g_free(blob);

	cr_assert_float_eq(layer_value(0), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(1), 0.28333f, 1e-4);
	cr_assert_float_eq(layer_value(2), 0.28333f, 1e-4,
	                   "the added layer must be scaled by the new solve");
	cr_assert_eq(joint_record_pin_count(jid), 3, "pins must follow the new list");
}

/* The anchor is the layer whose chain hosts the record, so it cannot be
 * dropped from the list: the record would have no chain to live on.  The
 * honest answer is to refuse and say what to do instead. */
Test(flis_layers_match, subset_amend_refuses_dropping_anchor) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();

	struct nde_joint_layers_match_data *old = joint_record_data(jid);
	cr_assert_not_null(old);
	struct nde_joint_layers_match_data *edited = nde_joint_layers_match_data_new(2);
	for (guint k = 0; k < 2; k++)
		copy_participant(&edited->parts[k], &old->parts[k + 1]);
	gchar *blob = op_desc_flis_layers_match.serialize(edited);
	nde_joint_layers_match_data_free(edited);
	nde_joint_layers_match_data_free(old);

	gchar *err = NULL;
	cr_assert(!nde_amend_execute(jid, blob, &err));
	cr_assert_not_null(err);
	cr_assert(strstr(err, "must stay a participant") != NULL,
	          "unexpected error: %s", err);
	g_free(err);
	g_free(blob);
}

/* The badgraph.png refinements: the band spans only its participants (a
 * luminance layer outside the group stays outside it), a composite/mask
 * captured AFTER the joint pins each participant at the JOINT record (the
 * replay-correctness half — without it the pin resolved the pre-calibration
 * state), and the graph binds such a pin to the band aligned on the
 * participant's column while a non-participant pin keeps its own node. */
Test(flis_layers_match, graph_binds_joint_consumers_to_the_band) {
	make_recorded_triple();
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.3f), "Lum");
	word[0] = "flis_layers_match";
	word[1] = "-subset=1,2,3";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);

	/* Pin fix: the joint record IS every participant's latest state. */
	for (gint item = 1; item <= 3; item++)
		cr_assert_eq(nde_history_last_record_for_item(item), jid,
		             "participant %d must pin at the joint record", item);
	cr_assert_neq(nde_history_last_record_for_item(4), jid,
	              "a non-participant must not pin at the joint record");

	/* A consumer capturing pins the way a composite does. */
	nde_pin_spec pins[2] = {
		{ "inG", 2, nde_history_last_record_for_item(2) },
		{ "inL", 4, nde_history_last_record_for_item(4) },
	};
	gint64 cons = nde_capture_structural_pinned("test.consumer", NDE_SCOPE_LAYER,
			1, NULL, "Consumer", pins, 2);
	cr_assert(cons > 0);

	nde_graph *graph = nde_graph_build();
	gint pid = nde_graph_joint_pseudo_item(jid);
	const nde_graph_node *jn = nde_graph_node_for(graph, pid);
	cr_assert_not_null(jn);
	cr_assert_eq(jn->n_span_items, 3, "the band spans its participants only");
	for (guint k = 0; k < 3; k++)
		cr_assert_eq(jn->span_items[k], (gint)k + 1);

	gboolean band_edge = FALSE, lum_edge = FALSE;
	for (guint i = 0; i < graph->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(graph->edges, i);
		if (e->dst_record_id != cons)
			continue;
		if (e->src_item_id == pid) {
			cr_assert_eq(e->src_align_item, 2,
			             "the band edge must align on its channel's column");
			band_edge = TRUE;
		} else if (e->src_item_id == 4) {
			cr_assert_eq(e->src_align_item, 0);
			lum_edge = TRUE;
		}
	}
	cr_assert(band_edge, "the participant pin must bind to the band");
	cr_assert(lum_edge, "the non-participant pin must bind to its own node");
	nde_graph_free(graph);
}

/* Reordering is defined on ONE chain — move a step past its neighbour — and
 * a record that is a member of several chains has no such neighbour.  It
 * must refuse rather than pick one chain's ordering and corrupt the rest. */
Test(flis_layers_match, joint_record_refuses_reorder) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 first = find_joint_record();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 second = find_joint_record();
	cr_assert(second > first);

	gchar *err = NULL;
	cr_assert(!nde_reorder_execute(second, first, FALSE, &err));
	cr_assert_not_null(err);
	cr_assert(strstr(err, "spans several layers") != NULL,
	          "unexpected error: %s", err);
	g_free(err);
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

	/* ONE joint "flis.group_calibration" record for the whole distribution,
	 * carrying the composite affine it was asked for — not the per-layer
	 * factors it happened to derive, which an upstream tint edit changes. */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	guint n = 0;
	for (guint i = 0; i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		cr_assert(g_strcmp0(rec->op_id, "flis.layer_scale"),
		          "group calibration must not capture per-layer records any more");
		if (g_strcmp0(rec->op_id, "flis.group_calibration"))
			continue;
		n++;
		cr_assert_eq(rec->tier, NDE_TIER_A);
		cr_assert_eq(rec->scope, NDE_SCOPE_LAYER);
		cr_assert_not_null(rec->params);
		cr_assert(strstr(rec->params, "n=3") != NULL);
		cr_assert(strstr(rec->params, "kind=3") != NULL,
		          "direct-affine capture must record kind DIRECT: %s", rec->params);
		cr_assert(strstr(rec->params, "dK0=2") != NULL,
		          "composite affine missing from params: %s", rec->params);
		cr_assert_not_null(rec->inputs);
		cr_assert_eq(rec->inputs->len, 3);
	}
	cr_assert_eq(n, 1, "expected ONE joint record, got %u", n);
	g_ptr_array_unref(snap);
}

/* A direct-affine group calibration keeps its K/O as the parameter, so
 * amending a tint upstream of it re-derives the DISTRIBUTION: R's tint
 * halved needs twice the offset to reach the same calibrated background,
 * and the layers that kept their tints are unmoved. */
Test(flis_layers_match, group_direct_tint_amend_redistributes) {
	flis_layer_t *r = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.2f), "R");
	set_and_record_tint(r, 1.0, 0.0, 0.0);
	flis_layer_t *g = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.3f), "G");
	set_and_record_tint(g, 0.0, 1.0, 0.0);
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.4f), "B");
	set_and_record_tint(b, 0.0, 0.0, 1.0);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	gint64 r_tint_rec = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "layer.set_tint") && rec->target_item_id == 1)
			r_tint_rec = rec->record_id;
	}
	if (snap)
		g_ptr_array_unref(snap);
	cr_assert(r_tint_rec > 0);

	GSList *members = g_slist_copy(com.uniq->layers);
	const double K[3] = { 2.0, 1.0, 0.5 };
	const double O[3] = { 0.01, 0.0, -0.01 };
	cr_assert_eq(flis_group_apply_channel_calibration(members, K, O, "Test calib"), 0);
	g_slist_free(members);
	cr_assert_float_eq(r->fit->fdata[0], 0.41f, 1e-5);

	gchar *err = NULL;
	gboolean ok = nde_amend_execute(r_tint_rec, "tinted=1;r=0.5;g=0;b=0", &err);
	cr_assert(ok, "tint amend failed: %s", err ? err : "?");
	g_free(err);

	cr_assert_float_eq(r->fit->fdata[0], 0.42f, 1e-4);
	cr_assert_float_eq(g->fit->fdata[0], 0.3f, 1e-4);
	cr_assert_float_eq(b->fit->fdata[0], 0.19f, 1e-4);
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

/* spcc-edit.png: after a FLATTEN, clicking "edit parameters" on the group
 * calibration step opened the raw key/value grid instead of the real SPCC
 * dialog.  The record is fine — the panel never asks for its editor:
 * on_hist_edit_clicked gates the native editor on
 * !nde_item_is_retained_input(target), and a flatten turns every participant
 * into a retained input, the anchor included.  That gate exists for the
 * LIVE-PREVIEW editors, which have nothing to preview against once the layer
 * is gone; the joint and photometric editors are apply-on-OK and do not care.
 * This pins the precondition down so the gate's fix cannot silently regress. */
Test(flis_layers_match, a_flatten_makes_the_joint_anchor_a_retained_input) {
	make_recorded_triple();
	cr_assert_eq(run_layers_match(), CMD_OK);
	gint64 jid = find_joint_record();
	cr_assert(jid > 0);

	gint anchor = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == jid)
			anchor = rec->target_item_id;
	}
	if (snap) g_ptr_array_unref(snap);
	cr_assert_eq(anchor, 1, "the anchor is the first participant");
	cr_assert(!nde_item_is_retained_input(anchor),
	          "before the flatten the anchor is an ordinary layer");

	cr_assert_eq(flis_flatten_all(), 0);

	cr_assert(nde_item_is_retained_input(anchor),
	          "after the flatten the anchor IS a retained input — this is the "
	          "condition that sent the SPCC step to the kv grid");
}
