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

/* Each band fills independently: a wide node in one band must not push the
 * nodes of another band along to meet it. */
Test(nde_graph, bands_fill_independently) {
	GArray *b = boxes_new();
	add_box(b, 1, 0, 300, 100);
	add_box(b, 2, 1, 40, 100);
	add_box(b, 3, 1, 40, 100);
	gint w = 0;
	GArray *p = nde_graph_layout(b, 5, 10, &w, NULL);

	cr_assert_eq(place_for(p, 2)->x, 0);
	cr_assert_eq(place_for(p, 3)->x, 45);
	cr_assert_eq(w, 300, "the widest band sets the width");
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
	cr_assert_eq(g_array_index(p, nde_graph_place, 0).x, 0);
	cr_assert_eq(g_array_index(p, nde_graph_place, 2).x, 14,
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
