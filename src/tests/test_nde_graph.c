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
	cr_assert_eq(g->depth, 1, "a single column");

	const nde_graph_node *n = g_ptr_array_index(g->nodes, 0);
	cr_assert_eq(n->item_id, NDE_ITEM_IMAGE);
	cr_assert_eq(n->kind, NDE_NODE_IMAGE);
	cr_assert(!n->orphan, "the image is what the Layers panel shows");
	cr_assert_eq(n->record_ids->len, 3, "all three steps live in it");

	nde_graph_free(g);
	done(f);
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
	 * sits one level past it — a left-to-right layout reads that off. */
	cr_assert_lt(mask->level, nde_graph_node_for(g, NDE_ITEM_IMAGE)->level);
	cr_assert_eq(g->depth, 2);

	nde_graph_free(g);
	done(f);
}

/* Order is meaning inside a node: it is what the user drags to reorder, so
 * the graph must not sort, dedupe or regroup a node's records. */
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
