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
 * test_flis_register_nde — the JOINT "flis.register" NDE record.
 *
 * Multi-layer registration is a multi-layer OPERATION: every participant's
 * transform is solved against the same reference, so recording it as N
 * independent steps would both misdescribe it and make undo N steps.  It
 * captures ONE joint record spanning every participant, exactly as layers
 * match and group calibration do (test_flis_layers_match.c is the model for
 * this harness).
 *
 * The load-bearing test here is `replay_reproduces_live_pixels`: it proves
 * the transform and the output size the record stores really are the ones
 * register_apply_reg applied.  Nothing else in the file means anything if
 * that one is weakened.
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
#include "core/op_descriptors.h"
#include "algos/geometry.h"
#include "registration/registration.h"
#include "registration/flis_register.h"

cominfo com;
fits *gfit;
char *word[MAX_COMMAND_WORDS];

TestSuite(flis_register_nde, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* ---- fixture ------------------------------------------------------------ */

/* Deterministic star field: gaussian spots from an LCG, content offset by
 * (sx, sy) in FITS buffer coordinates.  Same generator as
 * test_flis_register_layers.c — the DFT shift method needs real structure to
 * correlate on, and constant fixtures give it none. */
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

#define FIX_W 400
#define FIX_H 300

/* Two star fields, the second content-shifted, plus the memory-mode and
 * selection state the DFT path consults. */
static void make_shifted_pair(void) {
	flis_test_add_layer(make_star_field(FIX_W, FIX_H, 0, 0), "base");
	flis_test_add_layer(make_star_field(FIX_W, FIX_H, 15, 9), "shifted");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	com.uniq->canvas_w = FIX_W;
	com.uniq->canvas_h = FIX_H;
	com.pref.mem_mode = RATIO;
	com.pref.memory_ratio = 0.9;
	com.selection = (rectangle){ 100, 50, 128, 128 };
}

static void make_shifted_triple(void) {
	flis_test_add_layer(make_star_field(FIX_W, FIX_H, 0, 0), "base");
	flis_test_add_layer(make_star_field(FIX_W, FIX_H, 15, 9), "shifted");
	flis_test_add_layer(make_star_field(FIX_W, FIX_H, -7, 11), "other");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	com.uniq->canvas_w = FIX_W;
	com.uniq->canvas_h = FIX_H;
	com.pref.mem_mode = RATIO;
	com.pref.memory_ratio = 0.9;
	com.selection = (rectangle){ 100, 50, 128, 128 };
}

static int run_register(void) {
	selection_type sel;
	transformation_type tx;
	registration_function method =
	    flis_register_resolve_method(FLIS_REG_DFT, &sel, &tx);
	cr_assert_not_null(method);
	return flis_register_layers(NULL, NULL, method, sel, tx,
	                            OPENCV_LANCZOS4, TRUE, FLIS_REG_DFT);
}

/* The joint record in the live log, or 0. */
static gint64 find_register_record(void) {
	gint64 id = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (!g_strcmp0(rec->op_id, "flis.register"))
			id = rec->record_id;
	}
	if (snap)
		g_ptr_array_unref(snap);
	return id;
}

static gchar *record_params(gint64 rid) {
	gchar *out = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == rid)
			out = g_strdup(rec->params);
	}
	if (snap)
		g_ptr_array_unref(snap);
	return out;
}

static flis_layer_t *layer_nth(guint n) {
	return g_slist_nth_data(com.uniq->layers, n);
}

/* ---- 1. capture --------------------------------------------------------- */

/* ONE record for the whole registration: Tier A, scope LAYER, anchored on the
 * reference layer, one input pin per participant under the roles the graph
 * reads.  And no per-layer geometry records alongside it — that would make
 * undo N steps and describe an operation nobody performed. */
Test(flis_register_nde, captures_one_joint_record) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap, "no NDE history after registration");
	guint n_joint = 0;
	for (guint i = 0; i < snap->len; i++) {
		nde_record *rec = g_ptr_array_index(snap, i);
		if (g_strcmp0(rec->op_id, "flis.register"))
			continue;
		n_joint++;
		cr_assert_eq(rec->tier, NDE_TIER_A);
		cr_assert_eq(rec->scope, NDE_SCOPE_LAYER);
		cr_assert_eq(rec->target_item_id, 1,
		             "the anchor must be the reference layer");
		cr_assert_not_null(rec->params);
		cr_assert(strstr(rec->params, "n=2") != NULL,
		          "params '%s' missing participant count", rec->params);
		cr_assert(strstr(rec->params, "ref=1") != NULL,
		          "params '%s' missing the reference", rec->params);

		cr_assert_not_null(rec->inputs);
		cr_assert_eq(rec->inputs->len, 2, "every participant must be pinned");
		for (guint k = 0; k < rec->inputs->len; k++) {
			const nde_input_pin *pin = g_ptr_array_index(rec->inputs, k);
			char role[16];
			g_snprintf(role, sizeof(role), "in%u", k);
			cr_assert_str_eq(pin->role, role);
			cr_assert_eq(pin->src_item_id, (gint)k + 1);
		}
	}
	cr_assert_eq(n_joint, 1, "expected ONE joint record, got %u", n_joint);
	g_ptr_array_unref(snap);
}

/* The record's op must be registered in the descriptor table and round-trip
 * its own params — a joint record that cannot re-read its blob replays
 * nothing at all. */
Test(flis_register_nde, serializer_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	cr_assert_not_null(op, "flis.register not registered");
	cr_assert_not_null(op->serialize);
	cr_assert_not_null(op->deserialize);
	cr_assert(op->flags & OP_GEOMETRY_CHANGING,
	          "registration must count as geometry-changing");

	struct nde_joint_register_data *in = nde_joint_register_data_new(2);
	in->method = 2;
	in->tx_type = 1;
	in->interpolation = 4;
	in->clamp = TRUE;
	in->ref_item = 1;
	in->selection = (rectangle){ 3, 4, 5, 6 };
	in->canvas_w = 400;
	in->canvas_h = 300;
	for (guint k = 0; k < 2; k++) {
		in->parts[k].item_id = (gint)k + 1;
		in->parts[k].name = g_strdup(k ? "shifted" : "base");
		for (int c = 0; c < 9; c++)
			in->parts[k].H[c] = 0.125 * (c + 1) + k;
		in->parts[k].pos_x = -15 * (gint)k;
		in->parts[k].pos_y = 9 * (gint)k;
		in->parts[k].out_rx = 411;
		in->parts[k].out_ry = 312;
		g_free(in->parts[k].geom_sig);
		in->parts[k].geom_sig = g_strdup(k ? "cafe" : "");
	}
	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct nde_joint_register_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out, "deserialize failed on '%s'", blob);
	cr_assert_eq(out->n, 2u);
	cr_assert_eq(out->method, 2);
	cr_assert_eq(out->interpolation, 4);
	cr_assert_eq(out->clamp, TRUE);
	cr_assert_eq(out->ref_item, 1);
	cr_assert_eq(out->canvas_w, 400);
	cr_assert_eq(out->selection.w, 5);
	for (guint k = 0; k < 2; k++) {
		cr_assert_eq(out->parts[k].item_id, in->parts[k].item_id);
		cr_assert_str_eq(out->parts[k].name, in->parts[k].name);
		for (int c = 0; c < 9; c++)
			cr_assert_float_eq(out->parts[k].H[c], in->parts[k].H[c], 1e-12);
		cr_assert_eq(out->parts[k].out_rx, 411);
		cr_assert_eq(out->parts[k].pos_x, in->parts[k].pos_x);
		cr_assert_str_eq(out->parts[k].geom_sig, in->parts[k].geom_sig);
	}
	nde_joint_register_data_free(out);
	nde_joint_register_data_free(in);
	g_free(blob);

	/* Defensive: a blob missing a participant's output size is not a
	 * registration record, and must be rejected whole rather than replayed
	 * with a zero-sized warp. */
	cr_assert_null(op->deserialize("n=2;i0_item=1;i1_item=2", op->version));
	/* A layer listed twice would be warped twice. */
	struct nde_joint_register_data *dup = nde_joint_register_data_new(2);
	dup->ref_item = 1;
	for (guint k = 0; k < 2; k++) {
		dup->parts[k].item_id = 1;
		dup->parts[k].name = g_strdup("x");
		dup->parts[k].out_rx = dup->parts[k].out_ry = 16;
	}
	gchar *dupblob = op->serialize(dup);
	cr_assert_null(op->deserialize(dupblob, op->version));
	g_free(dupblob);
	nde_joint_register_data_free(dup);
}

/* ---- 2. chain membership ------------------------------------------------ */

/* One record at one log position, shared by N chains: what makes an undo
 * across it one step and lets any participant's chain replay through it. */
Test(flis_register_nde, record_is_member_of_every_participant_chain) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);

	for (gint item = 1; item <= 2; item++) {
		nde_chain *chain = nde_chain_build(item);
		gboolean member = FALSE;
		for (guint i = 0; i < chain->records->len; i++) {
			const nde_record *rec = g_ptr_array_index(chain->records, i);
			member = member || rec->record_id == rid;
		}
		cr_assert(member, "the record is missing from item %d's chain", item);
		cr_assert(chain->replayable, "item %d's chain must replay: %s", item,
		          chain->reasons->len ?
		              (const char *)g_ptr_array_index(chain->reasons, 0) : "?");
		cr_assert(chain->has_geometry,
		          "registration moves the layer — item %d's chain must carry "
		          "geometry so the replay threads its position", item);
		nde_chain_free(chain);
	}
}

/* ---- 3. THE L1 ORACLE (the one that matters) ---------------------------- */

/* Replaying a participant's chain must reproduce the live pixels.  This is
 * the whole justification for storing a homography rather than re-solving:
 * if the stored H is not the post-framing transform register_apply_reg
 * actually applied, or out_rx/out_ry are not the box it warped into, the
 * replayed layer comes out shifted, cropped, or both — and every later edit
 * silently inherits that.
 *
 * Nothing upstream has changed here, so this exercises L1 specifically. */
Test(flis_register_nde, replay_reproduces_live_pixels) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);
	cr_assert(find_register_record() > 0);

	for (guint n = 0; n < 2; n++) {
		flis_layer_t *lay = layer_nth(n);
		cr_assert_not_null(lay);
		cr_assert_not_null(lay->fit);

		nde_chain *chain = nde_chain_build(lay->item_id);
		cr_assert(chain->replayable, "item %d not replayable", lay->item_id);
		gchar *err = NULL;
		fits *replayed = nde_chain_replay(chain, &err);
		nde_chain_free(chain);
		cr_assert_not_null(replayed, "replay of layer %u failed: %s", n,
		                   err ? err : "?");
		g_free(err);

		cr_assert_eq(replayed->rx, lay->fit->rx,
		             "layer %u replayed %ux%u, live is %ux%u — the stored "
		             "output size is not the framed box", n,
		             replayed->rx, replayed->ry, lay->fit->rx, lay->fit->ry);
		cr_assert_eq(replayed->ry, lay->fit->ry);
		cr_assert_eq(replayed->naxes[2], lay->fit->naxes[2]);
		cr_assert_eq(replayed->type, DATA_FLOAT);

		size_t npix = (size_t)replayed->rx * replayed->ry * replayed->naxes[2];
		double worst = 0.0;
		size_t worst_i = 0;
		for (size_t i = 0; i < npix; i++) {
			double d = fabs((double)replayed->fdata[i] - (double)lay->fit->fdata[i]);
			if (d > worst) { worst = d; worst_i = i; }
		}
		/* "bit-for-bit-ish": the same H through the same interpolator on the
		 * same pixels is deterministic, so anything above float noise means
		 * a DIFFERENT transform was applied, not a rounding difference. */
		cr_assert(worst <= 1e-6,
		          "layer %u: replayed pixels differ from live by %g at index "
		          "%zu — the stored transform is not the one that was applied",
		          n, worst, worst_i);

		clearfits(replayed);
		free(replayed);
	}
}

/* The other half of a layer's value: the record has to put the layer back
 * where the live registration put it, not merely reproduce its pixels. */
Test(flis_register_nde, replay_restores_positions) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);

	gint live_x[2], live_y[2];
	for (guint n = 0; n < 2; n++) {
		live_x[n] = layer_nth(n)->position_x;
		live_y[n] = layer_nth(n)->position_y;
	}
	/* Scramble the offsets, then make the engine replay: a self-amend
	 * re-runs the record on every participant's chain. */
	for (guint n = 0; n < 2; n++) {
		layer_nth(n)->position_x = 12345;
		layer_nth(n)->position_y = -999;
	}

	gchar *params = record_params(rid);
	cr_assert_not_null(params);
	gchar *err = NULL;
	cr_assert(nde_amend_execute(rid, params, &err), "%s", err ? err : "?");
	g_free(err);
	g_free(params);

	for (guint n = 0; n < 2; n++) {
		cr_assert_eq(layer_nth(n)->position_x, live_x[n],
		             "layer %u position_x not restored by the replay", n);
		cr_assert_eq(layer_nth(n)->position_y, live_y[n],
		             "layer %u position_y not restored by the replay", n);
	}
}

/* ---- 4. delete ---------------------------------------------------------- */

/* Deleting the one record undoes the whole registration: each participant
 * replays its chain without it and returns to its original pixels, size AND
 * canvas offset. */
Test(flis_register_nde, delete_reverts_every_participant) {
	make_shifted_pair();
	fits *orig0 = NULL, *orig1 = NULL;
	orig0 = calloc(1, sizeof(fits));
	orig1 = calloc(1, sizeof(fits));
	cr_assert_eq(copyfits(layer_nth(0)->fit, orig0, CP_ALLOC | CP_COPYA | CP_FORMAT, -1), 0);
	cr_assert_eq(copyfits(layer_nth(1)->fit, orig1, CP_ALLOC | CP_COPYA | CP_FORMAT, -1), 0);

	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);
	cr_assert_neq(layer_nth(1)->position_x, 0,
	              "the fixture must actually move a layer for this to test "
	              "anything");

	gchar *err = NULL;
	cr_assert(nde_delete_execute(rid, &err), "delete failed: %s", err ? err : "?");
	g_free(err);
	cr_assert_eq(find_register_record(), 0);

	fits *orig[2] = { orig0, orig1 };
	for (guint n = 0; n < 2; n++) {
		flis_layer_t *lay = layer_nth(n);
		cr_assert_eq(lay->fit->rx, orig[n]->rx, "layer %u kept its framed width", n);
		cr_assert_eq(lay->fit->ry, orig[n]->ry);
		size_t npix = (size_t)orig[n]->rx * orig[n]->ry;
		double worst = 0.0;
		for (size_t i = 0; i < npix; i++)
			worst = MAX(worst, fabs((double)lay->fit->fdata[i] -
			                        (double)orig[n]->fdata[i]));
		cr_assert(worst <= 1e-6, "layer %u pixels not reverted (worst %g)", n, worst);
		cr_assert_eq(lay->position_x, 0, "layer %u offset not reverted", n);
		cr_assert_eq(lay->position_y, 0);
		clearfits(orig[n]);
		free(orig[n]);
	}
}

/* ---- 5. the graph ------------------------------------------------------- */

/* In the graph the record is not a step on any one column: it is a BAND drawn
 * across its participants, with no edges of its own. */
Test(flis_register_nde, graph_shows_one_spanning_node) {
	make_shifted_triple();
	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);
	gint pid = nde_graph_joint_pseudo_item(rid);

	nde_graph *graph = nde_graph_build();
	const nde_graph_node *jn = nde_graph_node_for(graph, pid);
	cr_assert_not_null(jn, "no dedicated node for the register record");
	cr_assert_eq(jn->kind, NDE_NODE_JOINT);
	cr_assert(jn->spanning, "the register record must be a spanning band");
	cr_assert_eq(jn->record_ids->len, 1);
	cr_assert_eq(jn->n_span_items, 3, "the band spans its participants");
	for (guint k = 0; k < 3; k++)
		cr_assert_eq(jn->span_items[k], (gint)k + 1);

	for (guint i = 0; i < graph->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(graph->edges, i);
		cr_assert(e->src_item_id != pid && e->dst_item_id != pid,
		          "a joint band must draw no edges");
	}
	nde_graph_free(graph);

	/* And the record IS every participant's latest state, so a later
	 * consumer pins there rather than at the pre-registration pixels. */
	for (gint item = 1; item <= 3; item++)
		cr_assert_eq(nde_history_last_record_for_item(item), rid,
		             "participant %d must pin at the register record", item);
}

/* Reordering is defined on ONE chain, and a record that is a member of
 * several has no single neighbour to move past. */
Test(flis_register_nde, record_refuses_reorder) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);
	gint64 first = find_register_record();
	cr_assert(first > 0);
	cr_assert_eq(run_register(), 0);
	gint64 second = find_register_record();
	cr_assert(second > first);

	gchar *err = NULL;
	cr_assert(!nde_reorder_execute(second, first, FALSE, &err));
	cr_assert_not_null(err);
	cr_assert(strstr(err, "spans several layers") != NULL,
	          "unexpected error: %s", err);
	g_free(err);
}

/* ---- 6. the geometry signature (the L1 / L2 discriminator) --------------- */

/* Nothing upstream: every participant's signature is empty, so nothing can
 * ever look like a change. */
Test(flis_register_nde, signature_empty_without_geometry) {
	make_shifted_pair();
	gchar *s1 = nde_joint_geometry_signature(1, 0);
	gchar *s2 = nde_joint_geometry_signature(2, 0);
	cr_assert_not_null(s1);
	cr_assert_str_eq(s1, "");
	cr_assert_str_eq(s2, "");
	g_free(s1);
	g_free(s2);
}

/* A geometry op upstream changes the signature; a NON-geometry op does not.
 * That asymmetry is the entire L1/L2 decision, so it is worth pinning
 * directly rather than only through a replay. */
Test(flis_register_nde, signature_follows_geometry_only) {
	make_shifted_pair();
	flis_layer_t *lay = layer_nth(1);

	struct flis_layer_scale_data *p = new_flis_layer_scale_data();
	p->scale = 0.5;
	flis_affine_layer_pixels(lay->fit, 0.5, 0.0);
	nde_checkpoint_baseline_ensure(lay->fit, lay->item_id);
	gint64 stretch = nde_capture_from_descriptor_for_item(
	    &op_desc_flis_layer_scale, p, "Halve", lay->fit, lay->item_id);
	free_flis_layer_scale_data(p);
	cr_assert(stretch > 0);

	gchar *after_stretch = nde_joint_geometry_signature(lay->item_id, 0);
	cr_assert_str_eq(after_stretch, "",
	                 "a pixel-only op must not enter the geometry signature");
	g_free(after_stretch);

	/* A geometry record on the same layer.  Captured directly (the record is
	 * what the signature reads; running the op itself is beside the point). */
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("geometry.rotation");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("angle=90;cropped=0;interpolation=4;clamp=1");
	rec->scope = NDE_SCOPE_CANVAS;
	rec->target_item_id = lay->item_id;
	rec->summary = g_strdup("Rotation");
	gint64 rot = nde_history_append(rec);
	cr_assert(rot > 0);

	gchar *after_rot = nde_joint_geometry_signature(lay->item_id, 0);
	cr_assert_str_neq(after_rot, "",
	                  "a geometry op must enter the geometry signature");
	/* And it is POSITIONAL: asking for the signature as of the rotation's own
	 * position must not include the rotation. */
	gchar *before_rot = nde_joint_geometry_signature(lay->item_id, rot);
	cr_assert_str_eq(before_rot, "");
	/* The other layer is unaffected — signatures are per-layer. */
	gchar *sibling = nde_joint_geometry_signature(1, 0);
	cr_assert_str_eq(sibling, "");
	g_free(after_rot);
	g_free(before_rot);
	g_free(sibling);
}

/* The record captures a signature per participant, so the replay has
 * something to compare against at all. */
Test(flis_register_nde, record_carries_participant_signatures) {
	make_shifted_pair();
	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	gchar *params = record_params(rid);
	cr_assert_not_null(params);
	cr_assert(strstr(params, "i0_geom") != NULL,
	          "no geometry signature recorded: %s", params);
	cr_assert(strstr(params, "i1_geom") != NULL);
	g_free(params);
}

/* The intent half (spec item 6): an upstream geometry AMEND flips a
 * participant's signature, so the record takes L2 and RE-SOLVES the
 * registration instead of replaying a transform that no longer describes the
 * pixels it is applied to.  An upstream NON-geometry amend leaves the
 * signature alone and stays on L1.
 *
 * Counted through nde_joint_analysis_runs(), which the register re-solve
 * increments exactly as the scalar joint analyses do. */
Test(flis_register_nde, geometry_amend_takes_l2_and_stretch_stays_l1) {
	make_shifted_pair();
	flis_layer_t *lay = layer_nth(1);
	nde_checkpoint_baseline_ensure(lay->fit, lay->item_id);
	nde_checkpoint_baseline_set_offset(lay->item_id, lay->position_x,
	                                   lay->position_y);

	/* A pixel-only step and a (currently no-op, full-frame) geometry step,
	 * both UPSTREAM of the registration. */
	flis_affine_layer_pixels(lay->fit, 0.9, 0.0);
	struct flis_layer_scale_data *p = new_flis_layer_scale_data();
	p->scale = 0.9;
	gint64 stretch = nde_capture_from_descriptor_for_item(
	    &op_desc_flis_layer_scale, p, "Dim", lay->fit, lay->item_id);
	free_flis_layer_scale_data(p);
	cr_assert(stretch > 0);

	struct crop_args *ca = new_crop_args();
	ca->area = (rectangle){ 0, 0, FIX_W, FIX_H };
	gint64 crop = nde_capture_from_descriptor_for_item(
	    &op_desc_crop, ca, "Crop", lay->fit, lay->item_id);
	if (ca->destroy_fn) ca->destroy_fn(ca);
	cr_assert(crop > 0);

	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);

	/* L1: the stretch is not geometry, so amending it must not re-solve —
	 * the transform still describes the pixels it is applied to. */
	guint before = nde_joint_analysis_runs();
	gchar *err = NULL;
	cr_assert(nde_amend_execute(stretch, "scale=0.8", &err), "%s", err ? err : "?");
	g_free(err);
	cr_assert_eq(nde_joint_analysis_runs(), before,
	             "a non-geometry amend must not re-run the registration");

	/* L2: shrinking the upstream crop changes the geometry the layer reaches
	 * the registration with, so the stored transform is stale and the solve
	 * must run again. */
	before = nde_joint_analysis_runs();
	err = NULL;
	gchar *newcrop = g_strdup_printf("x=0;y=0;w=%d;h=%d", FIX_W - 40, FIX_H - 30);
	gboolean amended = nde_amend_execute(crop, newcrop, &err);
	g_free(newcrop);
	cr_assert(amended, "upstream crop amend failed: %s", err ? err : "?");
	g_free(err);
	cr_assert(nde_joint_analysis_runs() > before,
	          "a changed geometry signature must re-run the registration");
}

/* One solve per edit however many participants replay: the register solution
 * is generation-cached for exactly the reason the scalar factors are. */
Test(flis_register_nde, one_solve_per_edit) {
	make_shifted_triple();
	flis_layer_t *lay = layer_nth(2);
	nde_checkpoint_baseline_ensure(lay->fit, lay->item_id);
	nde_checkpoint_baseline_set_offset(lay->item_id, lay->position_x,
	                                   lay->position_y);

	cr_assert_eq(run_register(), 0);
	gint64 rid = find_register_record();
	cr_assert(rid > 0);

	/* Force L2 by planting a geometry record upstream of the registration on
	 * one participant — its signature no longer matches the stored one. */
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("geometry.crop");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("x=0;y=0;w=400;h=300");
	rec->scope = NDE_SCOPE_CANVAS;
	rec->target_item_id = lay->item_id;
	rec->summary = g_strdup("Crop");
	cr_assert(nde_history_append(rec) > 0);

	guint before = nde_joint_analysis_runs();
	gchar *params = record_params(rid);
	gchar *err = NULL;
	gboolean amended = nde_amend_execute(rid, params, &err);
	g_free(params);
	g_free(err);
	guint delta = nde_joint_analysis_runs() - before;
	if (amended)
		cr_assert(delta <= 1,
		          "3 participants replayed but the registration solved %u times",
		          delta);
}

/* ---- 6. the amend dialog's params rewrite ------------------------------- */

/* A two-participant record with known settings and non-trivial stored state,
 * round-tripped through the descriptor so the tests below exercise the same
 * blob the editor sees.  @sig is what both participants' geom_sig start as. */
static struct nde_joint_register_data *make_settings_record(const char *sig) {
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	struct nde_joint_register_data *in = nde_joint_register_data_new(2);
	in->method = FLIS_REG_GLOBAL;
	in->tx_type = HOMOGRAPHY_TRANSFORMATION;
	in->interpolation = OPENCV_LANCZOS4;
	in->clamp = TRUE;
	in->ref_item = 1;
	in->selection = (rectangle){ 7, 8, 9, 10 };
	in->canvas_w = 512;
	in->canvas_h = 384;
	for (guint k = 0; k < 2; k++) {
		in->parts[k].item_id = (gint)k + 1;
		in->parts[k].name = g_strdup(k ? "green" : "red");
		for (int c = 0; c < 9; c++)
			in->parts[k].H[c] = 0.25 * (c + 1) - k;
		in->parts[k].pos_x = 11 * (gint)k;
		in->parts[k].pos_y = -7 * (gint)k;
		in->parts[k].out_rx = 520;
		in->parts[k].out_ry = 390;
		g_free(in->parts[k].geom_sig);
		in->parts[k].geom_sig = g_strdup(sig);
	}
	gchar *blob = op->serialize(in);
	nde_joint_register_data_free(in);
	struct nde_joint_register_data *out = op->deserialize(blob, op->version);
	g_free(blob);
	cr_assert_not_null(out);
	return out;
}

/* Changing a SETTING makes the stored transforms the answer to a question no
 * longer being asked, so the amend must stop replay reusing them.  It says so
 * by poisoning every participant's geometry signature: replay's L1 test then
 * cannot match and the alignment is solved again with the new settings.
 * Everything else in the blob must survive verbatim — the transforms stay as
 * the L3 fallback, and the framing, positions, selection and canvas size are
 * not the dialog's to touch. */
Test(flis_register_nde, amend_settings_change_invalidates_stored_solve) {
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	struct nde_joint_register_data *before = make_settings_record("deadbeef");
	struct nde_joint_register_data *p = make_settings_record("deadbeef");

	cr_assert(nde_joint_register_apply_settings(p, FLIS_REG_2PASS,
	                                            SHIFT_TRANSFORMATION,
	                                            OPENCV_LINEAR, FALSE, 2, NULL),
	          "a changed setting must report the solve as invalidated");

	/* Round-trip: the editor amends with serialize(), so the invalidation has
	 * to survive the blob, not merely the in-memory struct. */
	gchar *blob = op->serialize(p);
	struct nde_joint_register_data *out = op->deserialize(blob, op->version);
	g_free(blob);
	cr_assert_not_null(out, "the rewritten blob must still deserialize");

	cr_assert_eq(out->method, FLIS_REG_2PASS);
	cr_assert_eq(out->tx_type, SHIFT_TRANSFORMATION);
	cr_assert_eq(out->interpolation, OPENCV_LINEAR);
	cr_assert_eq(out->clamp, FALSE);
	cr_assert_eq(out->ref_item, 2);

	cr_assert_eq(out->n, before->n);
	cr_assert_eq(out->selection.x, before->selection.x);
	cr_assert_eq(out->selection.w, before->selection.w);
	cr_assert_eq(out->canvas_w, before->canvas_w);
	cr_assert_eq(out->canvas_h, before->canvas_h);
	for (guint k = 0; k < out->n; k++) {
		cr_assert_str_neq(out->parts[k].geom_sig, before->parts[k].geom_sig,
		                  "participant %u kept a signature that would send "
		                  "replay down L1 with a stale transform", k);
		/* The machine-derived state is preserved verbatim. */
		cr_assert_eq(out->parts[k].item_id, before->parts[k].item_id);
		cr_assert_str_eq(out->parts[k].name, before->parts[k].name);
		for (int c = 0; c < 9; c++)
			cr_assert_float_eq(out->parts[k].H[c], before->parts[k].H[c], 1e-12);
		cr_assert_eq(out->parts[k].pos_x, before->parts[k].pos_x);
		cr_assert_eq(out->parts[k].pos_y, before->parts[k].pos_y);
		cr_assert_eq(out->parts[k].out_rx, before->parts[k].out_rx);
		cr_assert_eq(out->parts[k].out_ry, before->parts[k].out_ry);
	}
	nde_joint_register_data_free(out);
	nde_joint_register_data_free(p);
	nde_joint_register_data_free(before);
}

/* The signature a layer with NO geometry upstream of the record computes is
 * the EMPTY STRING, and that is the ordinary case for a first registration.
 * So "clear the signature" cannot mean "set it to empty": an empty stored
 * signature compares EQUAL to a freshly computed empty one, replay takes L1,
 * and the stale transform is applied — the precise failure the invalidation
 * exists to prevent.  The marker must be something no signature can be. */
Test(flis_register_nde, invalidation_marker_cannot_be_a_real_signature) {
	struct nde_joint_register_data *p = make_settings_record("");
	/* KOMBAT rather than DFT: the fixture's selection is 9x10 and DFT wants a
	 * square one, so DFT would now be refused for a reason this test is not
	 * about. */
	cr_assert(nde_joint_register_apply_settings(p, FLIS_REG_KOMBAT,
	                                            SHIFT_TRANSFORMATION,
	                                            OPENCV_NEAREST, FALSE, 1, NULL));
	for (guint k = 0; k < p->n; k++) {
		cr_assert_str_neq(p->parts[k].geom_sig, "",
		                  "an empty marker IS a valid signature and would "
		                  "compare equal for an ungeometried layer");
		/* A real signature is empty or SHA256 hex; the marker must be neither. */
		gboolean hexish = strlen(p->parts[k].geom_sig) == 64;
		cr_assert(!hexish, "the marker must not look like a real digest");
	}
	nde_joint_register_data_free(p);

	/* And it really does drive replay to L2: a layer with no upstream geometry
	 * signs as "", which the poisoned value must not match. */
	gchar *live = nde_joint_geometry_signature(1, 0);
	cr_assert_str_eq(live, "", "no geometry upstream should sign as empty");
	cr_assert_str_neq(live, NDE_JOINT_GEOM_SIG_STALE);
	g_free(live);
}

/* The other half of the contract: an accidental OK must be a no-op.  When
 * every setting comes back unchanged the struct is left alone, so the
 * re-serialized blob is byte-identical and the amend is the plain "re-run"
 * verb rather than a silent forced re-solve. */
Test(flis_register_nde, amend_without_changes_preserves_the_blob) {
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	gchar *original = op->serialize(p);

	gchar *why = NULL;
	cr_assert_not(nde_joint_register_apply_settings(p, FLIS_REG_GLOBAL,
	                                                HOMOGRAPHY_TRANSFORMATION,
	                                                OPENCV_LANCZOS4, TRUE, 1,
	                                                &why),
	              "re-applying the recorded settings is not a change");
	cr_assert_null(why, "an unchanged amend is not an error: %s", why ? why : "");

	gchar *after = op->serialize(p);
	cr_assert_str_eq(after, original,
	                 "an unchanged amend must re-serialize byte-identically");
	for (guint k = 0; k < p->n; k++)
		cr_assert_str_eq(p->parts[k].geom_sig, "deadbeef",
		                 "an unchanged amend must keep the signatures so "
		                 "replay can still reuse the stored solve");
	g_free(after);
	g_free(original);
	nde_joint_register_data_free(p);
}

/* The reference must name one of the participants — the re-solve builds its
 * sequence around it, and the deserializer rejects a blob whose reference is
 * outside the list.  A refused rewrite must leave the struct untouched rather
 * than half-applied. */
Test(flis_register_nde, amend_rejects_a_reference_outside_the_participants) {
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	gchar *why = NULL;
	cr_assert_not(nde_joint_register_apply_settings(p, FLIS_REG_2PASS,
	                                                SHIFT_TRANSFORMATION,
	                                                OPENCV_LINEAR, FALSE, 99,
	                                                &why));
	cr_assert_not_null(why, "a refusal must say why, so the dialog can too");
	g_free(why);
	cr_assert_eq(p->method, FLIS_REG_GLOBAL, "a refused amend changed method");
	cr_assert_eq(p->ref_item, 1, "a refused amend changed the reference");
	cr_assert_eq(p->clamp, TRUE, "a refused amend changed clamp");
	for (guint k = 0; k < p->n; k++)
		cr_assert_str_eq(p->parts[k].geom_sig, "deadbeef",
		                 "a refused amend poisoned the signatures");
	nde_joint_register_data_free(p);
}

/* ===================================================================== */
/* The selection must travel with the method                             */
/*                                                                       */
/* A record written by a global star alignment stores NO selection — that */
/* method needs none.  Change its method to KOMBAT in the amend dialog    */
/* and the re-solve inherits {0,0,0,0}, which reaches cv::matchTemplate   */
/* as a 0x0 template.  That is an ASSERT, not an error return: OpenCV     */
/* throws, the exception crosses C frames, and the process aborts.        */
/* Nothing downstream can catch it, so it must never be reachable.        */
/* ===================================================================== */

/* Nothing to align to, no way to say so politely: refuse the amend. */
Test(flis_register_nde, amend_to_kombat_without_a_selection_is_refused) {
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	p->selection = (rectangle){ 0, 0, 0, 0 };   /* as a global-stars record has */
	memset(&com.selection, 0, sizeof com.selection);

	gchar *why = NULL;
	cr_assert_not(nde_joint_register_apply_settings(p, FLIS_REG_KOMBAT,
	                                               SHIFT_TRANSFORMATION,
	                                               OPENCV_LINEAR, FALSE, 1,
	                                               &why),
	              "a selectionless KOMBAT amend must be refused, not stored");
	cr_assert_not_null(why, "the refusal must be explained to the user");
	cr_assert_eq(p->method, FLIS_REG_GLOBAL,
	             "a refused amend must leave the method alone");
	g_free(why);
	nde_joint_register_data_free(p);
}

/* The selection the user has drawn is what they expect to be used, and it has
 * to end up IN THE RECORD: replay can run later, headless, with com.selection
 * empty or somewhere else entirely. */
Test(flis_register_nde, amend_to_kombat_captures_the_live_selection) {
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	p->selection = (rectangle){ 0, 0, 0, 0 };
	com.selection = (rectangle){ 40, 50, 64, 48 };

	gchar *why = NULL;
	cr_assert(nde_joint_register_apply_settings(p, FLIS_REG_KOMBAT,
	                                            SHIFT_TRANSFORMATION,
	                                            OPENCV_LINEAR, FALSE, 1, &why),
	          "a drawn selection makes the method change legal: %s",
	          why ? why : "?");

	gchar *blob = op->serialize(p);
	struct nde_joint_register_data *out = op->deserialize(blob, op->version);
	g_free(blob);
	cr_assert_not_null(out);
	cr_assert_eq(out->selection.x, 40, "the live selection must be recorded");
	cr_assert_eq(out->selection.y, 50);
	cr_assert_eq(out->selection.w, 64);
	cr_assert_eq(out->selection.h, 48);
	nde_joint_register_data_free(out);
	nde_joint_register_data_free(p);
	memset(&com.selection, 0, sizeof com.selection);
}

/* ...but re-opening a KOMBAT record to change something else must not drag its
 * selection off to wherever the cursor last left one. */
Test(flis_register_nde, an_unrelated_amend_keeps_the_recorded_selection) {
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	p->method = FLIS_REG_KOMBAT;
	p->selection = (rectangle){ 7, 8, 9, 10 };
	memset(&com.selection, 0, sizeof com.selection);

	cr_assert(nde_joint_register_apply_settings(p, FLIS_REG_KOMBAT,
	                                            SHIFT_TRANSFORMATION,
	                                            OPENCV_NEAREST, FALSE, 1, NULL),
	          "the interpolation moved, so this is a change");
	cr_assert_eq(p->selection.x, 7, "an unrelated amend moved the selection");
	cr_assert_eq(p->selection.w, 9);
	cr_assert_eq(p->selection.h, 10);
	nde_joint_register_data_free(p);
}

/* DFT is the stricter case: it needs a SQUARE selection, and a rectangular one
 * is refused rather than quietly squared off. */
Test(flis_register_nde, amend_to_dft_refuses_a_non_square_selection) {
	struct nde_joint_register_data *p = make_settings_record("deadbeef");
	com.selection = (rectangle){ 10, 10, 64, 48 };

	gchar *why = NULL;
	cr_assert_not(nde_joint_register_apply_settings(p, FLIS_REG_DFT,
	                                               SHIFT_TRANSFORMATION,
	                                               OPENCV_LINEAR, FALSE, 1,
	                                               &why));
	cr_assert_not_null(why);
	g_free(why);

	com.selection = (rectangle){ 10, 10, 64, 64 };   /* square: accepted */
	cr_assert(nde_joint_register_apply_settings(p, FLIS_REG_DFT,
	                                            SHIFT_TRANSFORMATION,
	                                            OPENCV_LINEAR, FALSE, 1, &why),
	          "%s", why ? why : "?");
	cr_assert_eq(p->selection.w, 64);
	cr_assert_eq(p->selection.h, 64);
	nde_joint_register_data_free(p);
	memset(&com.selection, 0, sizeof com.selection);
}

/* The validator itself, including the bound the crash needed: a template
 * larger than the image it is matched against fails the same assertion an
 * empty one does. */
Test(flis_register_nde, selection_validator_covers_every_fatal_shape) {
	rectangle empty = { 0, 0, 0, 0 };
	rectangle ok    = { 10, 10, 32, 32 };
	rectangle oblong= { 10, 10, 32, 16 };
	rectangle over  = { 10, 10, 200, 32 };
	rectangle neg   = { -5, 10, 32, 32 };
	gchar *err = NULL;

	cr_assert(flis_register_selection_ok(REQUIRES_NO_SELECTION, &empty, 0, 0, &err),
	          "a method needing no selection accepts an empty one");
	cr_assert_null(err);

	cr_assert_not(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &empty, 0, 0, &err));
	cr_assert_not_null(err, "an empty selection is the crash case; say so");
	g_free(err); err = NULL;

	cr_assert_not(flis_register_selection_ok(REQUIRES_ANY_SELECTION, NULL, 0, 0, &err),
	              "a NULL selection is treated as empty, not dereferenced");
	g_free(err); err = NULL;

	cr_assert(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &oblong, 0, 0, NULL),
	          "KOMBAT takes any rectangle");
	cr_assert_not(flis_register_selection_ok(REQUIRES_SQUARED_SELECTION, &oblong, 0, 0, NULL),
	              "DFT does not");
	cr_assert(flis_register_selection_ok(REQUIRES_SQUARED_SELECTION, &ok, 0, 0, NULL));

	/* Bounds: only checked when the image size is supplied. */
	cr_assert(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &over, 0, 0, NULL));
	cr_assert_not(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &over, 128, 128, &err),
	              "a template wider than the image trips the same assertion "
	              "an empty one does");
	g_free(err); err = NULL;
	cr_assert_not(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &neg, 128, 128, NULL),
	              "a selection starting off the left edge is out of bounds too");
	cr_assert(flis_register_selection_ok(REQUIRES_ANY_SELECTION, &ok, 128, 128, NULL),
	          "a selection inside the image is fine");
}

/* The load-bearing one: drive the solve exactly as the replay did when it
 * aborted — KOMBAT with an empty selection — and require an error return.
 * Before the guard this call did not fail, it killed the process, so a green
 * result here IS the regression check. */
Test(flis_register_nde, an_empty_selection_cannot_reach_matchtemplate) {
	fits *a = make_star_field(128, 128, 0.0, 0.0);
	fits *b = make_star_field(128, 128, 3.0, -2.0);
	cr_assert_not_null(a);
	cr_assert_not_null(b);
	fits *in[2] = { a, b };
	flis_reg_solution sol[2] = { 0 };

	selection_type sel = REQUIRES_NO_SELECTION;
	transformation_type tx = HOMOGRAPHY_TRANSFORMATION;
	registration_function kombat =
			flis_register_resolve_method(FLIS_REG_KOMBAT, &sel, &tx);
	cr_assert_not_null(kombat);
	cr_assert_eq(sel, REQUIRES_ANY_SELECTION);

	rectangle none = { 0, 0, 0, 0 };
	cr_assert_neq(flis_register_solve(in, 2, 0, kombat, sel, tx, none,
	                                  -1, 0, 0, sol), 0,
	              "KOMBAT with no selection must return an error — reaching "
	              "cv::matchTemplate with a 0x0 template aborts the process");

	/* And one that is present but too big for the image is just as fatal. */
	rectangle huge = { 0, 0, 512, 512 };
	cr_assert_neq(flis_register_solve(in, 2, 0, kombat, sel, tx, huge,
	                                  -1, 0, 0, sol), 0,
	              "a template larger than the image must be refused too");

	clearfits(a); free(a);
	clearfits(b); free(b);
}
