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
 * test_nde_graph — the history read as the DAG it is (design note §9.1):
 *   - a LINEAR document is one node, which is why it draws as today's list
 *   - a mask is its own node, joined to its consumer by a real edge
 *   - records stay in LOG ORDER inside a node (order is meaning there)
 *   - levels put a source before what consumes it, over forward edges only,
 *     because at ITEM granularity the graph really does contain cycles
 *   - items the Layers panel does not list are marked orphan, not "deleted"
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/op_descriptors.h"
#include "core/op_descriptor.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_graph.h"
#include "core/masks.h"
#include "filters/asinh.h"

cominfo com;
fits *gfit;

TestSuite(nde_graph, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* ---- fixtures (mirroring test_nde_replay's real-capture helpers) -------- */

static void fill_gradient(fits *f) {
	size_t n = (size_t)f->rx * f->ry;
	for (size_t i = 0; i < n; i++)
		f->fdata[i] = (float)i / (float)n;
}

static asinh_params *asinh_beta(float beta) {
	asinh_params *p = calloc(1, sizeof(*p));
	p->beta = beta;
	p->clip_mode = RESCALE;
	return p;
}

static mask_from_channel_data *from_channel(int chan) {
	mask_from_channel_data *d = calloc(1, sizeof(*d));
	d->channel = chan;
	d->bitpix = 8;
	return d;
}

static mask_blur_data *blur_mask(float radius) {
	mask_blur_data *d = calloc(1, sizeof(*d));
	d->radius = radius;
	return d;
}

static int run_op(const op_descriptor *op, gpointer user, gboolean masked) {
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->op = op;
	args->user = user;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;
	args->mask_aware = masked;
	args->max_threads = com.max_thread;
	args->mem_ratio = -1.0f;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_image_worker(args));
	com.headless = prev;
	return rc;
}

static int run_mask_op(const op_descriptor *op, gpointer user) {
	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->op = op;
	args->user = user;
	args->max_threads = com.max_thread;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_mask_worker(args));
	com.headless = prev;
	return rc;
}

static fits *fresh_image(void) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_gradient(f);
	gfit = f;
	return f;
}

static void done(fits *f) {
	nde_history_attach(NULL);
	clearfits(f); free(f);
	gfit = NULL;
}

/* ---- the linear case --------------------------------------------------- */

/* The property the whole design rests on: with nothing but a chain of image
 * operations the graph is a SINGLE node holding a single list — so a linear
 * document renders exactly as today's History window, and nobody pays for
 * structure they do not have. */
Test(nde_graph, a_linear_document_is_one_node) {
	fits *f = fresh_image();
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(10.f), FALSE), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), FALSE), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(30.f), FALSE), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 1, "one item, one node");
	cr_assert_eq(g->edges->len, 0, "and nothing to draw between");
	cr_assert_eq(g->depth, 1, "a single band");

	const nde_graph_node *n = g_ptr_array_index(g->nodes, 0);
	cr_assert_eq(n->item_id, NDE_ITEM_IMAGE);
	cr_assert_eq(n->kind, NDE_NODE_IMAGE);
	cr_assert(!n->orphan, "the image is what the Layers panel shows");
	cr_assert_eq(n->record_ids->len, 3, "all three steps live in it");

	nde_graph_free(g);
	done(f);
}

/* The user-visible half of the origin record: an image that has had nothing
 * done to it still has a node, so adding a second layer to a single-layer
 * document shows BOTH — the panel used to show the new layer and not the one
 * that was already there. */
Test(nde_graph, an_image_has_a_node_as_soon_as_it_says_where_it_came_from) {
	cr_assert_neq(nde_capture_image_origin("file", "x.fit", "Opened 'x.fit'"), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 1, "the image is in the graph before it is edited");
	const nde_graph_node *n = g_ptr_array_index(g->nodes, 0);
	cr_assert_eq(n->item_id, NDE_ITEM_IMAGE);
	cr_assert_eq(n->record_ids->len, 1, "and its origin is the step it shows");
	nde_graph_free(g);
}

/* And it travels: promoting a plain image to a layered one rebinds the whole
 * plain-image history onto the base layer, origin record included, so the
 * layer that was already there is a node from the moment the document has
 * two.  This is the case the report named. */
Test(nde_graph, a_promoted_image_takes_its_origin_with_it) {
	fits *f = fresh_image();
	com.uniq->fit = f;
	com.uniq->chans = 1;
	cr_assert_neq(nde_capture_image_origin("file", "x.fit", "Opened 'x.fit'"), 0);
	cr_assert_eq(flis_promote_from_gfit("Background"), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 1);
	const nde_graph_node *n = g_ptr_array_index(g->nodes, 0);
	cr_assert_neq(n->item_id, NDE_ITEM_IMAGE, "rebound onto the base layer");
	cr_assert_eq(n->kind, NDE_NODE_LAYER);
	cr_assert_str_eq(n->label, "Background");
	cr_assert_eq(n->record_ids->len, 1);
	nde_graph_free(g);
	/* the base layer owns f now; the suite teardown frees it */
	nde_history_attach(NULL);
}

Test(nde_graph, an_empty_history_is_an_empty_graph) {
	fits *f = fresh_image();
	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 0);
	cr_assert_eq(g->edges->len, 0);
	cr_assert_eq(g->depth, 1);
	nde_graph_free(g);
	done(f);
}

/* ---- the case the list could not show ---------------------------------- */

/* A mask is a second node joined to its consumer by a real edge — the thing
 * the flat list had to smuggle into a target column. */
Test(nde_graph, a_mask_is_its_own_node_with_an_edge_to_its_consumer) {
	fits *f = fresh_image();
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), TRUE), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 2, "image and mask");

	const nde_graph_node *mask = nde_graph_node_for(g, NDE_ITEM_PLAIN_MASK);
	cr_assert_not_null(mask);
	cr_assert_eq(mask->kind, NDE_NODE_MASK);
	cr_assert(mask->orphan, "no Layers row lists a mask");

	/* BOTH directions are real, and together they form a cycle at item
	 * granularity: the mask was built FROM the image, and the image then
	 * consumed the mask.  Neither edge is spurious — the cycle is what
	 * collapsing records into items costs, and is why levels are computed
	 * over forward edges only. */
	cr_assert_eq(g->edges->len, 2);
	gboolean saw_mask_edge = FALSE, saw_image_edge = FALSE;
	for (guint i = 0; i < g->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(g->edges, i);
		if (e->src_item_id == NDE_ITEM_PLAIN_MASK && e->dst_item_id == NDE_ITEM_IMAGE) {
			saw_mask_edge = TRUE;
			cr_assert_str_eq(e->role, "mask");
		}
		if (e->src_item_id == NDE_ITEM_IMAGE && e->dst_item_id == NDE_ITEM_PLAIN_MASK) {
			saw_image_edge = TRUE;
			cr_assert_str_eq(e->role, "image");
			cr_assert_eq(e->src_record_id, 0, "built from the image's baseline");
		}
	}
	cr_assert(saw_mask_edge, "the mask supplies the masked step");
	cr_assert(saw_image_edge, "and was itself built from the image");

	/* The mask's chain starts first here, so it ranks first and the image
	 * sits one band below it — a top-to-bottom layout reads that off. */
	cr_assert_lt(mask->level, nde_graph_node_for(g, NDE_ITEM_IMAGE)->level);
	cr_assert_eq(g->depth, 2);

	nde_graph_free(g);
	done(f);
}

/* Order is meaning inside a node: it is what the user drags to reorder, so
 * the graph must not sort, dedupe or regroup a node's records. */
/* A node with no records of its own — an image consumed only as a mask's
 * source, a layer known only from a composite's pins — has no log position to
 * rank by, and discovery order appends it AFTER the record that named it.
 * Ranked that way, its derivation edge read as feedback and the mask sat
 * BESIDE its source with the arrow drawn backwards into empty space.  A
 * record-less node is a pure source (every edge destination is some record's
 * target), so it ranks before the whole log: what derives from it sits below
 * it, and the edge is forward. */
Test(nde_graph, a_recordless_source_ranks_before_what_derives_from_it) {
	fits *f = fresh_image();
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_eq(g->nodes->len, 2, "image and mask");
	const nde_graph_node *img  = nde_graph_node_for(g, NDE_ITEM_IMAGE);
	const nde_graph_node *mask = nde_graph_node_for(g, NDE_ITEM_PLAIN_MASK);
	cr_assert_not_null(img);
	cr_assert_not_null(mask);
	cr_assert_eq(img->record_ids->len, 0u,
	             "precondition: the image has no records of its own");
	cr_assert_lt(img->level, mask->level,
	             "the source sits above what was derived from it");
	cr_assert_eq(g->depth, 2);

	cr_assert_eq(g->edges->len, 1);
	const nde_graph_edge *e = g_ptr_array_index(g->edges, 0);
	cr_assert_eq(e->src_item_id, NDE_ITEM_IMAGE);
	cr_assert_eq(e->dst_item_id, NDE_ITEM_PLAIN_MASK);
	cr_assert(!nde_graph_edge_is_feedback(g, e),
	          "derivation is a forward edge, not feedback");
	nde_graph_free(g);
	done(f);
}

/* graph2.png: a processing mask created on a layer, used by two steps, then
 * CLEARED — its item resolves to nothing, and the fallback called it
 * "Layer 5, no longer in the document": a layer the user never had.  The
 * graph knows better: the node's records are mask ops and a record consumed
 * it through a "mask" pin, so it is named as the mask of what it masked. */
Test(nde_graph, a_cleared_mask_item_is_labelled_as_a_mask) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 12, 0.f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	fill_gradient(gfit);

	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	gint mask_item = flis_layer_pmask_id(base);
	cr_assert_neq(mask_item, 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), TRUE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_clear, NULL), 0);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *mn = nde_graph_node_for(g, mask_item);
	cr_assert_not_null(mn);
	cr_assert_eq(mn->kind, NDE_NODE_MASK,
	             "a cleared mask item is still a mask, not an unknown layer");
	cr_assert_not_null(strstr(mn->label, "Mask of base"),
	                   "…and says whose it was, got: %s", mn->label);
	nde_graph_free(g);
	nde_history_attach(NULL);
	gfit = NULL;
}

/* procmasksnag: a mask cleared without EVER having been used leaves nothing
 * behind — no step consumed it, so its records describe an item that
 * influenced no pixels, and their only future was an orphan node.  (A USED
 * mask keeps everything and goes dormant instead — the test above.) */
Test(nde_graph, clearing_a_never_used_mask_forgets_it_entirely) {
	fits *f = fresh_image();
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(10.f), FALSE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_clear, NULL), 0);

	nde_graph *g = nde_graph_build();
	cr_assert_null(nde_graph_node_for(g, NDE_ITEM_PLAIN_MASK),
	               "nothing consumed the mask, so its episode is forgotten");
	const nde_graph_node *img = nde_graph_node_for(g, NDE_ITEM_IMAGE);
	cr_assert_not_null(img);
	cr_assert_eq(img->record_ids->len, 1, "the image's own history is untouched");
	nde_graph_free(g);
	done(f);
}

Test(nde_graph, records_keep_log_order_inside_their_node) {
	fits *f = fresh_image();
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(10.f), FALSE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), TRUE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_blur, blur_mask(3.f)), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(30.f), TRUE), 0);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *img  = nde_graph_node_for(g, NDE_ITEM_IMAGE);
	const nde_graph_node *mask = nde_graph_node_for(g, NDE_ITEM_PLAIN_MASK);
	cr_assert_eq(img->record_ids->len, 3, "three image steps");
	cr_assert_eq(mask->record_ids->len, 2, "two mask steps");

	for (guint i = 1; i < img->record_ids->len; i++)
		cr_assert_lt(*(gint64 *)g_ptr_array_index(img->record_ids, i - 1),
		             *(gint64 *)g_ptr_array_index(img->record_ids, i),
		             "ascending, as logged");
	for (guint i = 1; i < mask->record_ids->len; i++)
		cr_assert_lt(*(gint64 *)g_ptr_array_index(mask->record_ids, i - 1),
		             *(gint64 *)g_ptr_array_index(mask->record_ids, i));

	/* The two chains INTERLEAVE in the log, which is exactly what a single
	 * flat list showed and a per-node view deliberately separates. */
	cr_assert_lt(*(gint64 *)g_ptr_array_index(mask->record_ids, 0),
	             *(gint64 *)g_ptr_array_index(img->record_ids, 1));

	nde_graph_free(g);
	done(f);
}

/* ---- the queries the GUI needs ----------------------------------------- */

/* Consumers-of is reverse invalidation made visible: these are the nodes an
 * amend would re-run, and the graph can say so before the user commits. */
Test(nde_graph, consumers_and_inputs_are_both_reachable) {
	fits *f = fresh_image();
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), TRUE), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(30.f), TRUE), 0);

	nde_graph *g = nde_graph_build();

	GPtrArray *consumers = nde_graph_consumers_of(g, NDE_ITEM_PLAIN_MASK);
	cr_assert_eq(consumers->len, 2, "both masked steps consume the mask");
	g_ptr_array_unref(consumers);

	GPtrArray *inputs = nde_graph_inputs_of(g, NDE_ITEM_IMAGE);
	cr_assert_eq(inputs->len, 2);
	for (guint i = 0; i < inputs->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(inputs, i);
		cr_assert_eq(e->src_item_id, NDE_ITEM_PLAIN_MASK);
	}
	g_ptr_array_unref(inputs);

	/* The image has a consumer too — the mask was derived from it.  That is
	 * the edge which, together with the one above, makes the item-level graph
	 * cyclic. */
	GPtrArray *from_image = nde_graph_consumers_of(g, NDE_ITEM_IMAGE);
	cr_assert_eq(from_image->len, 1, "the mask derives from the image");
	cr_assert_str_eq(((nde_graph_edge *)g_ptr_array_index(from_image, 0))->role, "image");
	g_ptr_array_unref(from_image);

	nde_graph_free(g);
	done(f);
}

/* Node order must not reshuffle as the history grows, or the layout would
 * jump under the user between one operation and the next. */
Test(nde_graph, node_order_is_stable_as_records_are_added) {
	fits *f = fresh_image();
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(10.f), FALSE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);

	nde_graph *g1 = nde_graph_build();
	gint first_a = ((nde_graph_node *)g_ptr_array_index(g1->nodes, 0))->item_id;
	gint second_a = ((nde_graph_node *)g_ptr_array_index(g1->nodes, 1))->item_id;
	nde_graph_free(g1);

	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(20.f), TRUE), 0);
	cr_assert_eq(run_mask_op(&op_desc_mask_blur, blur_mask(3.f)), 0);

	nde_graph *g2 = nde_graph_build();
	cr_assert_eq(g2->nodes->len, 2);
	cr_assert_eq(((nde_graph_node *)g_ptr_array_index(g2->nodes, 0))->item_id, first_a);
	cr_assert_eq(((nde_graph_node *)g_ptr_array_index(g2->nodes, 1))->item_id, second_a);
	nde_graph_free(g2);

	done(f);
}

/* ---- layout ------------------------------------------------------------
 * Pure: no history, no widgets.  The view measures its children and hands
 * the sizes here, which is the whole reason this is a function and not a
 * size_allocate() implementation.
 */

static GArray *boxes_new(void) {
	return g_array_new(FALSE, TRUE, sizeof(nde_graph_box));
}

static void add_box(GArray *a, gint item, gint level, gint w, gint h) {
	nde_graph_box b = { .item_id = item, .level = level, .w = w, .h = h };
	g_array_append_val(a, b);
}

static void add_sat(GArray *a, gint item, gint host, gint w, gint h) {
	nde_graph_box b = { .item_id = item, .host_item = host, .w = w, .h = h };
	g_array_append_val(a, b);
}

static void add_span(GArray *a, gint item, gint level, gint w, gint h,
                     const gint *items, guint n_items) {
	nde_graph_box b = { .item_id = item, .level = level, .w = w, .h = h,
	                    .spanning = TRUE, .span_items = items,
	                    .n_span_items = n_items };
	g_array_append_val(a, b);
}

static const nde_graph_place *place_for(GArray *places, gint item) {
	for (guint i = 0; i < places->len; i++) {
		const nde_graph_place *p = &g_array_index(places, nde_graph_place, i);
		if (p->item_id == item)
			return p;
	}
	return NULL;
}

/* Derivation runs downward, so nodes that derive from nothing share the top
 * band and sit side by side. */
Test(nde_graph, one_level_lays_out_as_a_single_band) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_box(b, 2, 0, 240, 60);
	gint w = 0, h = 0;
	GArray *p = nde_graph_layout(b, 20, 6, &w, &h);

	cr_assert_eq(p->len, 2);
	cr_assert_eq(place_for(p, 1)->y, 0);
	cr_assert_eq(place_for(p, 2)->y, 0, "one level is one band");
	cr_assert_eq(place_for(p, 1)->x, 0);
	cr_assert_eq(place_for(p, 2)->x, 220, "side by side with the column gap between");
	cr_assert_eq(h, 100, "the band is as tall as its tallest node");
	cr_assert_eq(w, 460);
	g_array_unref(p);
	g_array_unref(b);
}

/* A band is levelled to its tallest member so that every edge leaving it
 * leaves at one y — the reason places carry a height at all. */
Test(nde_graph, a_band_is_as_tall_as_its_tallest_node) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 50, 100);
	add_box(b, 2, 0, 50, 180);
	add_box(b, 3, 1, 50, 90);
	gint h = 0;
	GArray *p = nde_graph_layout(b, 20, 6, NULL, &h);

	cr_assert_eq(place_for(p, 1)->h, 180);
	cr_assert_eq(place_for(p, 2)->h, 180);
	cr_assert_eq(place_for(p, 3)->y, 186, "180 tall plus the 6 gap");
	cr_assert_eq(h, 276);
	g_array_unref(p);
	g_array_unref(b);
}

/* Nodes pack tightly WITHIN a band, but a column belongs to the item that
 * claimed it, so a band's nodes start beyond any column another item already
 * holds above them.  (This test used to assert the opposite — bands filled
 * independently and reused each other's columns — which is what put a
 * newly-added layer underneath Background's history.) */
Test(nde_graph, a_band_packs_tightly_but_clears_columns_already_held) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 300, 100);
	add_box(b, 2, 1, 40, 100);
	add_box(b, 3, 1, 40, 100);
	gint w = 0;
	GArray *p = nde_graph_layout(b, 5, 10, &w, NULL);

	cr_assert_eq(place_for(p, 2)->x, 305,
	             "past item 1's column, which it does not own");
	cr_assert_eq(place_for(p, 3)->x, 350, "then tight against its band-mate");
	cr_assert_eq(w, 390, "the widest band sets the width");
	g_array_unref(p);
	g_array_unref(b);
}

/* A SPANNING band (a joint record, nde_joint.h) stretches over exactly its
 * participants' columns — a layer the operation never read stays outside it
 * — and its vertical whitespace is the COLUMN gap, not the band gap. */
Test(nde_graph, a_spanning_band_covers_its_participants_columns) {
	GArray *b = boxes_new();
	/* Band 0: luminance (1) then R/G/B (2,3,4), each 100 wide, gap 20. */
	add_box(b, 1, 0, 100, 50);   /* x 0   */
	add_box(b, 2, 0, 100, 50);   /* x 120 */
	add_box(b, 3, 0, 100, 50);   /* x 240 */
	add_box(b, 4, 0, 100, 50);   /* x 360 */
	const gint parts[3] = { 2, 3, 4 };
	add_span(b, -1001, 1, 80, 30, parts, 3);
	add_box(b, 5, 2, 100, 40);   /* something after the band */
	gint w = 0, h = 0;
	GArray *p = nde_graph_layout(b, 20, 36, &w, &h);

	const nde_graph_place *band = place_for(p, -1001);
	cr_assert_eq(band->x, 120, "the band starts at its first participant");
	cr_assert_eq(band->w, 340, "and ends at its last participant's right edge");
	/* Whitespace: 50-tall band 0, then the COLUMN gap (20), not the row
	 * gap (36) — on both sides of the spanning band. */
	cr_assert_eq(band->y, 70);
	cr_assert_eq(place_for(p, 5)->y, 70 + 30 + 20);
	cr_assert_eq(h, 70 + 30 + 20 + 40);
	g_array_unref(p);
	g_array_unref(b);
}

/* With no participant columns placed (all merged away, or no list at all)
 * the band degrades to the full layout width. */
Test(nde_graph, a_spanning_band_without_participants_spans_everything) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 100, 50);
	add_box(b, 2, 0, 100, 50);
	const gint gone[1] = { 99 };
	add_span(b, -1001, 1, 80, 30, gone, 1);
	gint w = 0;
	GArray *p = nde_graph_layout(b, 20, 36, &w, NULL);

	cr_assert_eq(place_for(p, -1001)->x, 0);
	cr_assert_eq(place_for(p, -1001)->w, w);
	g_array_unref(p);
	g_array_unref(b);
}

/* Node order in equals node order out, so the layout is a function of the
 * graph's stable order and nothing else. */
Test(nde_graph, layout_preserves_the_order_it_was_given) {
	GArray *b = boxes_new();
	add_box(b, 7, 1, 10, 10);
	add_box(b, 3, 0, 10, 10);
	add_box(b, 5, 1, 10, 10);
	GArray *p = nde_graph_layout(b, 4, 4, NULL, NULL);

	cr_assert_eq(g_array_index(p, nde_graph_place, 0).item_id, 7);
	cr_assert_eq(g_array_index(p, nde_graph_place, 1).item_id, 3);
	cr_assert_eq(g_array_index(p, nde_graph_place, 2).item_id, 5);
	/* Item 3 is alone in band 0 and takes the first column; the band-1 pair
	 * clear it and then pack against each other, in list order. */
	cr_assert_eq(g_array_index(p, nde_graph_place, 1).x, 0);
	cr_assert_eq(g_array_index(p, nde_graph_place, 0).x, 14);
	cr_assert_eq(g_array_index(p, nde_graph_place, 2).x, 28,
	             "second in its band, whatever its position in the list");
	g_array_unref(p);
	g_array_unref(b);
}

/* The defect wrongorder.png showed: a mask laid out after every band node had
 * its connectors crossing the layers in between to get back to its host.  A
 * satellite goes beside its host, and the band's next node beyond it. */
Test(nde_graph, a_satellite_sits_between_its_host_and_the_next_node) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 300);      /* the host layer */
	add_box(b, 2, 0, 200, 100);      /* the next layer along */
	add_sat(b, 9, 1, 150, 60);       /* a mask of layer 1 */
	gint w = 0;
	GArray *p = nde_graph_layout(b, 20, 6, &w, NULL);

	cr_assert_eq(place_for(p, 1)->x, 0);
	cr_assert_eq(place_for(p, 9)->x, 220, "beside its host, not after the band");
	cr_assert_eq(place_for(p, 2)->x, 390, "the next layer starts beyond the mask");
	cr_assert_eq(place_for(p, 9)->y, place_for(p, 1)->y, "in its host's band");
	cr_assert_eq(w, 590);
	g_array_unref(p);
	g_array_unref(b);
}

/* Several masks of one layer are a column, and the band has to be tall enough
 * to hold it — otherwise the last of them overhangs the band below. */
Test(nde_graph, satellites_stack_and_the_band_grows_to_hold_them) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 50);
	add_sat(b, 9, 1, 100, 60);
	add_sat(b, 8, 1, 140, 70);
	add_box(b, 2, 1, 200, 40);
	gint h = 0;
	GArray *p = nde_graph_layout(b, 20, 6, NULL, &h);

	cr_assert_eq(place_for(p, 9)->y, 0);
	cr_assert_eq(place_for(p, 8)->y, 80, "60 tall plus the 20 column gap");
	cr_assert_eq(place_for(p, 1)->h, 150, "the band covers the whole column");
	cr_assert_eq(place_for(p, 2)->y, 156, "the next band clears it");
	cr_assert_eq(h, 196);
	g_array_unref(p);
	g_array_unref(b);
}

/* A satellite may be slid down to meet the host row it masks, but only as far
 * as leaves room for the satellites still below it in the column. */
Test(nde_graph, a_satellites_slack_stops_short_of_the_ones_below_it) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 400);
	add_sat(b, 9, 1, 100, 60);
	add_sat(b, 8, 1, 100, 70);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, 9)->y_min, 0);
	cr_assert_eq(place_for(p, 9)->y_max, 400 - 60 - 90,
	             "leaves the 70-tall mask below it, and the 20 gap, in the band");
	cr_assert_eq(place_for(p, 8)->y_max, 400 - 70);
	cr_assert_eq(place_for(p, 1)->y_min, place_for(p, 1)->y_max,
	             "a band node has no slack at all");
	g_array_unref(p);
	g_array_unref(b);
}

/* Levelling a band gives every node in it the same bottom edge, which is where
 * the old code started the outgoing line — far below a short node's last step,
 * with a gap of nothing in between.  The node's own height is kept for that. */
Test(nde_graph, a_place_remembers_the_height_the_node_measured) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 50, 100);
	add_box(b, 2, 0, 50, 180);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, 1)->h, 180, "levelled to the band");
	cr_assert_eq(place_for(p, 1)->content_h, 100, "but its own history ends here");
	cr_assert_eq(place_for(p, 2)->content_h, 180);
	g_array_unref(p);
	g_array_unref(b);
}

/* A satellite whose host is not in the layout would otherwise go unplaced.
 * Wrong-looking beats invisible. */
Test(nde_graph, a_satellite_with_no_host_lays_out_as_an_ordinary_node) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_sat(b, 9, 77, 150, 60);
	gint w = 0;
	GArray *p = nde_graph_layout(b, 20, 6, &w, NULL);

	cr_assert_not_null(place_for(p, 9));
	cr_assert_eq(place_for(p, 9)->x, 220, "second in the band, on its own account");
	cr_assert_eq(w, 370);
	g_array_unref(p);
	g_array_unref(b);
}

Test(nde_graph, an_empty_graph_lays_out_to_nothing) {
	GArray *b = boxes_new();
	gint w = -1, h = -1;
	GArray *p = nde_graph_layout(b, 8, 8, &w, &h);
	cr_assert_eq(p->len, 0);
	cr_assert_eq(w, 0);
	cr_assert_eq(h, 0);
	g_array_unref(p);
	g_array_unref(b);
}

/* The edge that makes the item-level graph cyclic must be identifiable, so
 * that it can be drawn differently and kept out of the layout. */
Test(nde_graph, the_mask_feedback_edge_is_identified_as_one) {
	fits *f = fresh_image();
	cr_assert_eq(run_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(run_op(&op_desc_asinh, asinh_beta(10.f), TRUE), 0);

	nde_graph *g = nde_graph_build();
	guint forward = 0, feedback = 0;
	for (guint i = 0; i < g->edges->len; i++) {
		const nde_graph_edge *e = g_ptr_array_index(g->edges, i);
		if (nde_graph_edge_is_feedback(g, e))
			feedback++;
		else
			forward++;
	}
	cr_assert_gt(forward, 0, "the image feeds the mask");
	cr_assert_gt(feedback, 0, "and the mask feeds the image back");
	nde_graph_free(g);
	done(f);
}

/* ===================================================================== */
/* A segment stays in its item's column                                  */
/*                                                                       */
/* badorder.png: an SPCC over layers 2 and 3, then a median filter on    */
/* layer 2.  The median is a SEGMENT — layer 2's history continued below */
/* the joint band — and it was drawn under layer 1, because bands are    */
/* left-aligned in isolation and it was the only node in its own band.   */
/* A step of layer 2 must sit in layer 2's column, or the history reads  */
/* as belonging to whichever layer happens to be leftmost.               */
/* ===================================================================== */

static void add_segment(GArray *a, gint item, gint level, gint w, gint h,
                        gint align_item) {
	nde_graph_box b = { .item_id = item, .level = level, .w = w, .h = h,
	                    .align_item = align_item };
	g_array_append_val(a, b);
}

Test(nde_graph, a_segment_sits_in_its_own_items_column) {
	const gint participants[] = { 2, 3 };
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);                       /* Background */
	add_box(b, 2, 0, 200, 100);                       /* r_layers_00002 */
	add_box(b, 3, 0, 200, 100);                       /* r_layers_00003 */
	add_span(b, -1, 1, 200, 30, participants, 2);     /* the SPCC band */
	add_segment(b, -2, 2, 200, 40, 2);                /* #5 Median on layer 2 */
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, 2)->x, 220, "layer 2 is the middle column");
	cr_assert_eq(place_for(p, -2)->x, place_for(p, 2)->x,
	             "the median filter belongs to layer 2, so it must sit in "
	             "layer 2's column — not under layer 1 at x=0");
	cr_assert_gt(place_for(p, -2)->y, place_for(p, -1)->y,
	             "and still below the joint band that split it off");
	g_array_unref(p);
	g_array_unref(b);
}

/* Two participants both continuing below the band each keep their own
 * column, rather than being packed left to right in list order. */
Test(nde_graph, sibling_segments_keep_their_separate_columns) {
	const gint participants[] = { 2, 3 };
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_box(b, 2, 0, 200, 100);
	add_box(b, 3, 0, 200, 100);
	add_span(b, -1, 1, 200, 30, participants, 2);
	add_segment(b, -2, 2, 200, 40, 2);
	add_segment(b, -3, 2, 200, 40, 3);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, -2)->x, place_for(p, 2)->x);
	cr_assert_eq(place_for(p, -3)->x, place_for(p, 3)->x);
	cr_assert_neq(place_for(p, -2)->x, place_for(p, -3)->x,
	              "two segments must not be stacked into one column");
	g_array_unref(p);
	g_array_unref(b);
}

/* An item whose history STARTS after the joint band has no column to
 * inherit, so it flows as an ordinary node — but it must not land on top of
 * a segment that has already claimed a column in that band. */
Test(nde_graph, a_new_item_below_a_band_does_not_collide_with_a_segment) {
	const gint participants[] = { 1, 2 };
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_box(b, 2, 0, 200, 100);
	add_span(b, -1, 1, 200, 30, participants, 2);
	add_segment(b, -2, 2, 200, 40, 1);   /* claims column 0 */
	add_box(b, 4, 2, 200, 40);           /* a layer added after the SPCC */
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, -2)->x, 0, "the segment keeps layer 1's column");
	cr_assert_geq(place_for(p, 4)->x, 220,
	              "a new node must start beyond the segment, not under it");
	g_array_unref(p);
	g_array_unref(b);
}

/* A segment whose anchor is missing (its item's first node filtered out of
 * the view) must still be placed: visible in the wrong column beats
 * unplaced. */
Test(nde_graph, a_segment_without_its_anchor_falls_back_to_the_flow) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_segment(b, -2, 1, 200, 40, 77);   /* no box for item 77 */
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_not_null(place_for(p, -2), "an unresolved segment is still placed");
	cr_assert_eq(place_for(p, -2)->x, 220,
	             "with no column to inherit it flows as an ordinary node, "
	             "which means clearing item 1's column rather than landing "
	             "on top of it");
	g_array_unref(p);
	g_array_unref(b);
}

/* ===================================================================== */
/* A column belongs to the item that claimed it                          */
/*                                                                       */
/* Bands used to fill independently, so an item whose history STARTS      */
/* below another item's band was left-aligned into that item's column.    */
/* A layer added after a joint band landed under Background and read as   */
/* part of its history.  A column is the width of one item's story down   */
/* the whole page, so once claimed it is not available to anyone else.    */
/* ===================================================================== */

static void add_retired(GArray *a, gint item, gint level, gint w, gint h) {
	nde_graph_box b = { .item_id = item, .level = level, .w = w, .h = h,
	                    .retired = TRUE };
	g_array_append_val(a, b);
}

Test(nde_graph, a_new_item_lower_down_does_not_reuse_an_occupied_column) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);   /* Background */
	add_box(b, 2, 0, 200, 100);   /* L2 */
	add_box(b, 3, 2, 200, 100);   /* L3, added after a joint band */
	gint w = 0;
	GArray *p = nde_graph_layout(b, 20, 6, &w, NULL);

	cr_assert_eq(place_for(p, 1)->x, 0);
	cr_assert_eq(place_for(p, 2)->x, 220);
	cr_assert_eq(place_for(p, 3)->x, 440,
	             "L3 starts its own column rather than sitting under "
	             "Background, whose column is occupied further up");
	cr_assert_eq(w, 640);
	g_array_unref(p);
	g_array_unref(b);
}

/* The claim is per ITEM, not per band, so a node may reuse the column of the
 * item it continues — that is what a segment does — while still being kept
 * out of anyone else's. */
Test(nde_graph, an_item_may_reoccupy_its_own_column_lower_down) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_box(b, 2, 0, 200, 100);
	add_segment(b, -2, 1, 200, 40, 2);   /* layer 2, continued */
	add_box(b, 3, 1, 200, 100);          /* a genuinely new item */
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, -2)->x, 220, "its own column, reoccupied");
	cr_assert_eq(place_for(p, 3)->x, 440, "the new item takes a fresh one");
	g_array_unref(p);
	g_array_unref(b);
}

/* A satellite column is claimed too: it is part of its host's column group,
 * and a later node walking left to right must clear the whole group. */
Test(nde_graph, a_later_node_clears_a_hosts_satellite_column_too) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_sat(b, 9, 1, 150, 60);
	add_box(b, 2, 1, 200, 100);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, 9)->x, 220, "the mask sits beside its host");
	cr_assert_eq(place_for(p, 2)->x, 390,
	             "the next item starts beyond the mask's column, not under it");
	g_array_unref(p);
	g_array_unref(b);
}

/* ===================================================================== */
/* Removed layers get out of the way                                     */
/*                                                                       */
/* A layer the document no longer holds keeps its node — the records are  */
/* real and still listed — but it is no longer part of the live story, so */
/* it goes to the far right where it does not push the live columns       */
/* around or sit between them.                                           */
/* ===================================================================== */

Test(nde_graph, a_removed_layer_goes_to_the_right_of_every_live_one) {
	GArray *b = boxes_new();
	add_retired(b, 9, 0, 200, 100);   /* first in the order, but gone */
	add_box(b, 1, 0, 200, 100);
	add_box(b, 2, 1, 200, 100);
	gint w = 0;
	GArray *p = nde_graph_layout(b, 20, 6, &w, NULL);

	cr_assert_eq(place_for(p, 1)->x, 0, "the live layer keeps the left edge");
	cr_assert_eq(place_for(p, 2)->x, 220,
	             "and the second live one takes the next column");
	cr_assert_gt(place_for(p, 9)->x, place_for(p, 2)->x,
	             "the removed layer is out of the way on the right, however "
	             "early it appears in the order");
	cr_assert_eq(w, 640);
	g_array_unref(p);
	g_array_unref(b);
}

/* Removed layers do not pile into one column either — they are out of the
 * way, not merged. */
Test(nde_graph, removed_layers_keep_separate_columns_from_each_other) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_retired(b, 8, 0, 200, 100);
	add_retired(b, 9, 1, 200, 100);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, 1)->x, 0);
	cr_assert_gt(place_for(p, 8)->x, 0);
	cr_assert_neq(place_for(p, 9)->x, place_for(p, 8)->x);
	g_array_unref(p);
	g_array_unref(b);
}

/* A segment of a removed layer follows its anchor to the right, rather than
 * being placed among the live columns on its own account. */
Test(nde_graph, a_removed_layers_segment_follows_it_right) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 200, 100);
	add_retired(b, 8, 0, 200, 100);
	add_segment(b, -8, 1, 200, 40, 8);
	GArray *p = nde_graph_layout(b, 20, 6, NULL, NULL);

	cr_assert_eq(place_for(p, -8)->x, place_for(p, 8)->x);
	cr_assert_gt(place_for(p, -8)->x, place_for(p, 1)->x);
	g_array_unref(p);
	g_array_unref(b);
}

/* ===================================================================== */
/* Consumed is not deleted                                               */
/*                                                                       */
/* A layer flattened away and a layer the user deleted both leave the     */
/* document, so describe_item() calls both NDE_NODE_UNKNOWN.  They are    */
/* not the same thing to read: the consumed one is where the surviving    */
/* image CAME FROM and belongs among the live columns with the composite  */
/* below it, while the deleted one is a dead end kept only so it can be   */
/* undone.  Telling them apart by "no longer in the document" put the     */
/* flatten's inputs out on the right and left its result stranded.        */
/* ===================================================================== */

Test(nde_graph, a_flattened_away_layer_is_not_a_deleted_one) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.3f), "A");
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.4f), "B");
	cr_assert_eq(flis_flatten_all(), 0);

	nde_graph *g = nde_graph_build();
	guint consumed = 0;
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->kind != NDE_NODE_UNKNOWN)
			continue;
		consumed++;
		cr_assert_not(n->deleted,
		              "'%s' was consumed by the flatten, not deleted — it "
		              "must stay among the live columns", n->label);
	}
	cr_assert_gt(consumed, 0, "the flatten should have consumed its inputs");
	nde_graph_free(g);
	nde_history_attach(NULL);
}

/* The flatten's result reads as the place its inputs arrive at, so it sits
 * under them instead of claiming a column of its own to the right. */
Test(nde_graph, a_flatten_result_sits_under_the_layers_it_consumed) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.3f), "A");
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.4f), "B");
	cr_assert_eq(flis_flatten_all(), 0);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *result = NULL;
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->kind == NDE_NODE_LAYER && n->level > 0)
			result = n;
	}
	cr_assert_not_null(result, "the flatten should leave one live layer below");
	cr_assert_neq(result->column_item, 0,
	              "a derived node must name the column it belongs in, or the "
	              "layout will start it a column to the right of its inputs");

	/* And that column is one it actually consumed, at a lower level. */
	const nde_graph_node *anchor = nde_graph_node_for(g, result->column_item);
	cr_assert_not_null(anchor);
	cr_assert_lt(anchor->level, result->level);
	gboolean is_input = FALSE;
	for (guint e = 0; e < g->edges->len && !is_input; e++) {
		const nde_graph_edge *ed = g_ptr_array_index(g->edges, e);
		is_input = ed->dst_item_id == result->item_id &&
		           ed->src_item_id == anchor->item_id;
	}
	cr_assert(is_input, "the column it takes must be one of its own inputs");
	nde_graph_free(g);
	nde_history_attach(NULL);
}

/* A deleted layer, by contrast, IS flagged — it keeps its node so the
 * deletion can be reversed, but it is out of the live story. */
Test(nde_graph, a_removed_layer_is_flagged_as_deleted) {
	flis_layer_t *a = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.3f), "A");
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.4f), "B");
	cr_assert_not_null(a);
	const gint gone = b->item_id;

	GString *kv = nde_kv_start();
	nde_kv_add_str(kv, "name", "B");
	cr_assert_neq(nde_capture_structural("layer.remove", NDE_SCOPE_DOCUMENT,
	                                     gone, nde_kv_end(kv), "Remove layer"), 0);
	cr_assert_eq(flis_layer_remove(b), 0);

	nde_graph *g = nde_graph_build();
	const nde_graph_node *n = nde_graph_node_for(g, gone);
	cr_assert_not_null(n, "a deleted layer keeps its node so the deletion "
	                      "can be reversed");
	cr_assert(n->deleted, "…and is marked as deleted, not merely gone");
	nde_graph_free(g);
	nde_history_attach(NULL);
}
