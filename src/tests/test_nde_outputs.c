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
 * test_nde_outputs — which items a record WRITES (nde-simplify-plan step 5):
 *   - the invariant every older reader depends on: output 0 is the target,
 *     however the producer ordered its outputs
 *   - a record that never said keeps costing nothing and writes its target
 *   - a joint record's participants become its outputs when it is appended,
 *     and when a file written before OUTPUTS existed is attached
 *   - the list survives a FLIS save and load
 *   - a secondary output gets its own graph column and one edge into it
 *   - and, until the replay can carry N results, a multi-output record is
 *     refused from a chain and from the edit path rather than half-done
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_joint.h"
#include "core/nde/nde_graph.h"
#include "core/nde/nde_replay.h"
#include "core/nde/nde_checkpoint.h"
#include "core/op_descriptor.h"
#include "core/processing.h"

cominfo com;
fits *gfit;

static char *tmpdir = NULL;
static char *tmppath = NULL;

static void setup(void) {
	flis_test_init_com();
	tmpdir = g_dir_make_tmp("nde-outputs-XXXXXX", NULL);
	tmppath = g_build_filename(tmpdir, "test.flis", NULL);
}

static void teardown(void) {
	nde_history_attach(NULL);
	if (tmppath) { g_unlink(tmppath); g_free(tmppath); tmppath = NULL; }
	if (tmpdir)  { g_rmdir(tmpdir);   g_free(tmpdir);   tmpdir  = NULL; }
	flis_test_cleanup_com();
}

TestSuite(nde_outputs, .init = setup, .fini = teardown);

/* ---- the model, on a bare record --------------------------------------- */

/* Output 0 is target_item_id even when the producer names a secondary first:
 * the TARGET column, chain membership and the graph's node placement all read
 * it, and none of them should depend on the order a capture site happened to
 * list its outputs in. */
Test(nde_outputs, an_output_list_always_begins_with_the_target) {
	nde_record *rec = nde_record_new();
	rec->target_item_id = 7;
	nde_record_add_output(rec, "starmask", 9);
	nde_record_add_output(rec, "starless", 7);

	cr_assert_eq(nde_record_output_count(rec), 2);
	cr_assert_eq(nde_record_output(rec, 0)->item_id, 7);
	/* ...and naming it afterwards relabelled the seeded primary rather than
	 * listing item 7 twice. */
	cr_assert_str_eq(nde_record_output(rec, 0)->role, "starless");
	cr_assert_eq(nde_record_output(rec, 1)->item_id, 9);
	nde_record_free(rec);
}

/* The 95% case: no list, no memory, and it still answers the only question
 * anybody asks. */
Test(nde_outputs, a_record_writes_its_target_without_saying_so) {
	nde_record *rec = nde_record_new();
	rec->target_item_id = 3;
	cr_assert_null(rec->outputs);
	cr_assert_eq(nde_record_output_count(rec), 1);
	cr_assert_null(nde_record_output(rec, 0),
	               "there is no pin to describe an output nobody wrote down");
	cr_assert(nde_record_writes_item(rec, 3));
	cr_assert_not(nde_record_writes_item(rec, 4));
	nde_record_free(rec);
}

Test(nde_outputs, writes_item_covers_every_output) {
	nde_record *rec = nde_record_new();
	rec->target_item_id = 3;
	nde_record_add_output(rec, "out", 3);
	nde_record_add_output(rec, "extra", 8);
	cr_assert(nde_record_writes_item(rec, 3));
	cr_assert(nde_record_writes_item(rec, 8));
	cr_assert_not(nde_record_writes_item(rec, 4));
	nde_record_free(rec);
}

Test(nde_outputs, a_copied_record_carries_its_outputs) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("test.two_outputs");
	rec->target_item_id = 3;
	nde_record_add_output(rec, "out", 3);
	nde_record_add_output(rec, "extra", 8);
	nde_record *copy = nde_record_copy(rec);
	cr_assert_eq(nde_record_output_count(copy), 2);
	cr_assert(nde_record_writes_item(copy, 8));
	cr_assert_neq(copy->outputs, rec->outputs, "a copy, not a shared list");
	nde_record_free(copy);
	nde_record_free(rec);
}

/* ---- joint records ------------------------------------------------------ */

/* Params naming @n participants, in the shape nde_joint_params_participants
 * reads.  Only the participant list matters here. */
static gchar *joint_params(const gint *items, guint n) {
	GString *s = g_string_new(NULL);
	g_string_append_printf(s, "n=%u", n);
	for (guint k = 0; k < n; k++)
		g_string_append_printf(s, ";i%u_item=%d;i%u_s=1", k, items[k], k);
	return g_string_free(s, FALSE);
}

static gint64 append_joint(const gint *items, guint n, gint target) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("flis.layers_match");
	rec->op_version = 1;
	rec->params = joint_params(items, n);
	rec->summary = g_strdup("Layers match");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = target;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	return nde_history_append(rec);
}

/* A joint record writes every layer it read.  The capture sites say so once,
 * in the params they were already writing; the edges are derived on append, so
 * nothing downstream parses that blob to find out. */
Test(nde_outputs, a_joint_record_writes_all_its_participants) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "a");
	gint items[3] = { 1, 2, 3 };
	gint64 id = append_joint(items, 3, 2);

	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert_eq(rec->record_id, id);
	cr_assert_eq(nde_record_output_count(rec), 3);
	/* Target first, as ever — even though it is participant 1 in the params. */
	cr_assert_eq(nde_record_output(rec, 0)->item_id, 2);
	for (int i = 1; i <= 3; i++)
		cr_assert(nde_record_writes_item(rec, i),
		          "participant %d is written by the joint record", i);
	cr_assert_not(nde_record_writes_item(rec, 4));
	g_ptr_array_unref(snap);
}

/* A file written before the OUTPUTS column says who it wrote in its params and
 * nowhere else.  Attaching derives the edges once, so the runtime never meets a
 * joint record without them and the legacy case is not a second code path. */
Test(nde_outputs, a_legacy_joint_record_gains_outputs_when_attached) {
	/* Built the way the FLIS loader builds one, so that the record reaches
	 * attach exactly as it would off disk: no outputs, only params. */
	nde_history *h = g_new0(nde_history, 1);
	h->records = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
	h->next_record_id = 2;
	nde_record *rec = nde_record_new();
	rec->record_id = 1;
	rec->op_id = g_strdup("flis.layers_match");
	gint items[2] = { 5, 6 };
	rec->params = joint_params(items, 2);
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = 5;
	g_ptr_array_add(h->records, rec);
	h->live_count = 1;
	cr_assert_null(rec->outputs, "as it came off disk");

	nde_history_attach(h);
	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *loaded = g_ptr_array_index(snap, 0);
	cr_assert_eq(nde_record_output_count(loaded), 2);
	cr_assert(nde_record_writes_item(loaded, 6));
	g_ptr_array_unref(snap);
}

/* ---- persistence -------------------------------------------------------- */

Test(nde_outputs, outputs_survive_a_save_and_a_load) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;   /* @l does not survive flis_free_layers below */
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("test.two_outputs");
	rec->op_version = 1;
	rec->summary = g_strdup("Two images");
	rec->tier = NDE_TIER_B;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = lid;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	nde_record_add_output(rec, "starless", lid);
	nde_record_add_output(rec, "starmask", lid + 100);
	nde_history_append(rec);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	const nde_record *back = g_ptr_array_index(snap, 0);
	cr_assert_eq(nde_record_output_count(back), 2);
	cr_assert_eq(nde_record_output(back, 0)->item_id, lid,
	             "output 0 is still the target after a round trip");
	cr_assert_str_eq(nde_record_output(back, 0)->role, "starless");
	cr_assert_eq(nde_record_output(back, 1)->item_id, lid + 100);
	cr_assert_str_eq(nde_record_output(back, 1)->role, "starmask");
	g_ptr_array_unref(snap);
}

/* An ordinary record writes an empty cell, and reads back as what it is: a
 * single-output record with no list at all.  This is also what a v2 file looks
 * like to a v3 reader — the column is simply absent — so the degradation path
 * is the one every save already exercises. */
Test(nde_outputs, a_single_output_record_stores_no_list) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;   /* @l does not survive flis_free_layers below */
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("stretch.mtf");
	rec->op_version = 1;
	rec->params = g_strdup("linked=1");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = lid;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	nde_history_append(rec);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *back = g_ptr_array_index(snap, 0);
	cr_assert_null(back->outputs);
	cr_assert_eq(nde_record_output_count(back), 1);
	cr_assert(nde_record_writes_item(back, lid));
	g_ptr_array_unref(snap);
}

/* ---- the graph ---------------------------------------------------------- */

/* A secondary output is an ordinary outgoing edge into a column of its own —
 * not a spanning band.  A band is one op across columns that already existed
 * and continue afterwards; a multi-output op CREATES the column. */
Test(nde_outputs, the_graph_draws_one_edge_per_secondary_output) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint produced = l->item_id + 100;
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("test.two_outputs");
	rec->tier = NDE_TIER_B;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = l->item_id;
	rec->timestamp = nde_iso8601_now();
	nde_record_add_output(rec, "starless", l->item_id);
	nde_record_add_output(rec, "starmask", produced);
	gint64 id = nde_history_append(rec);

	nde_graph *g = nde_graph_build();
	guint found = 0;
	for (guint i = 0; i < g->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(g->edges, i);
		if (e->dst_item_id != produced)
			continue;
		found++;
		cr_assert_eq(e->src_item_id, l->item_id);
		cr_assert_eq(e->src_record_id, id);
		cr_assert_eq(e->dst_record_id, id);
		cr_assert_str_eq(e->role, "starmask");
	}
	cr_assert_eq(found, 1, "one edge per secondary output, no more");
	cr_assert_not_null(nde_graph_node_for(g, produced),
	                   "the produced item gets a column of its own");
	nde_graph_free(g);
}

/* ---- the refusals, until the replay can carry N results ----------------- */

/* An ordinary op's hook returns ONE image, so replaying a record that wrote
 * two would reproduce one of them and claim to have reproduced the step.  The
 * chain says so rather than finding out halfway through. */
Test(nde_outputs, a_multi_output_record_is_not_replayable) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("stretch.mtf");   /* a real, replayable descriptor */
	rec->op_version = 1;
	rec->params = g_strdup("linked=1;lo=0;mid=0.5;hi=1");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = l->item_id;
	rec->timestamp = nde_iso8601_now();
	nde_record_add_output(rec, "out", l->item_id);
	nde_record_add_output(rec, "extra", l->item_id + 100);
	nde_history_append(rec);

	nde_chain *c = nde_chain_build(l->item_id);
	cr_assert_not(c->replayable);
	cr_assert_gt(c->reasons->len, 0);
	cr_assert(strstr(g_ptr_array_index(c->reasons, 0), "more than one image"),
	          "the reason names the actual problem: %s",
	          (const char *)g_ptr_array_index(c->reasons, 0));
	nde_chain_free(c);
}

/* Amending one would have to replay and commit N chains atomically.  Refuse it
 * in words the user can act on rather than editing one output of two. */
Test(nde_outputs, editing_a_multi_output_record_is_refused) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("stretch.mtf");
	rec->op_version = 1;
	rec->params = g_strdup("linked=1;lo=0;mid=0.5;hi=1");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = l->item_id;
	rec->timestamp = nde_iso8601_now();
	nde_record_add_output(rec, "out", l->item_id);
	nde_record_add_output(rec, "extra", l->item_id + 100);
	gint64 id = nde_history_append(rec);

	gchar *err = NULL;
	cr_assert_not(nde_amend_execute(id, "linked=1;lo=0;mid=0.6;hi=1", &err));
	cr_assert_not_null(err);
	cr_assert(strstr(err, "more than one image"), "got: %s", err);
	g_free(err);
}

/* The guard above must not catch the family that already works: a joint record
 * writes all its participants and keeps its own multi-participant edit path. */
Test(nde_outputs, a_joint_record_is_not_refused_for_having_many_outputs) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.50f), "b");
	gint items[2] = { a->item_id, b->item_id };
	gint64 id = append_joint(items, 2, a->item_id);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_gt(nde_record_output_count(g_ptr_array_index(snap, 0)), 1);
	g_ptr_array_unref(snap);

	gchar *err = NULL;
	gchar *params = joint_params(items, 2);
	nde_amend_execute(id, params, &err);
	/* It may well fail for want of baselines in this fixture; what it must not
	 * do is refuse on the grounds of having several outputs. */
	if (err)
		cr_assert_not(strstr(err, "more than one image"),
		              "a joint record was refused as multi-output: %s", err);
	g_free(params);
	g_free(err);
}
