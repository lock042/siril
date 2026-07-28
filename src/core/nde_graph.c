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
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* See nde_graph.h for what this derives and why it is separate from the GUI. */

#include "core/siril.h"
#include "core/nde_history.h"
#include "core/nde_graph.h"
#include "io/image_format_flis.h"

static void node_free(gpointer p) {
	nde_graph_node *n = p;
	if (!n)
		return;
	g_free(n->label);
	if (n->record_ids)
		g_ptr_array_unref(n->record_ids);
	g_free(n);
}

static void edge_free(gpointer p) {
	nde_graph_edge *e = p;
	if (!e)
		return;
	g_free(e->role);
	g_free(e);
}

void nde_graph_free(nde_graph *g) {
	if (!g)
		return;
	g_ptr_array_unref(g->nodes);
	g_ptr_array_unref(g->edges);
	g_free(g);
}

/* Name and classify an item.  The two sentinels come first because a plain
 * image has no FLIS document to ask. */
static void describe_item(gint item_id, nde_node_kind *kind, gchar **label,
                          gboolean *orphan) {
	if (item_id == NDE_ITEM_IMAGE) {
		*kind = NDE_NODE_IMAGE;
		*label = g_strdup(_("Image"));
		*orphan = FALSE;
		return;
	}
	if (item_id == NDE_ITEM_PLAIN_MASK) {
		*kind = NDE_NODE_MASK;
		*label = g_strdup(_("Mask"));
		*orphan = TRUE;
		return;
	}

	flis_layer_t *owner = NULL;
	switch (flis_item_lookup(item_id, &owner)) {
	case FLIS_ITEM_LAYER:
		*kind = NDE_NODE_LAYER;
		*label = g_strdup(owner && owner->layer_name ? owner->layer_name : _("Layer"));
		*orphan = FALSE;
		return;
	case FLIS_ITEM_MASK:
		*kind = NDE_NODE_MASK;
		/* Named after what it masks: "Mask" alone is ambiguous once a
		 * document has several. */
		*label = owner && owner->layer_name
		         ? g_strdup_printf(_("Mask of %s"), owner->layer_name)
		         : g_strdup(_("Mask"));
		*orphan = TRUE;
		return;
	case FLIS_ITEM_LMASK:
		*kind = NDE_NODE_LAYERMASK;
		*label = owner && owner->layer_name
		         ? g_strdup_printf(_("Layer mask of %s"), owner->layer_name)
		         : g_strdup(_("Layer mask"));
		*orphan = TRUE;
		return;
	default:
		break;
	}

	/* An id the document no longer claims — most often a layer consumed by a
	 * merge or flatten.  Its records are real and still listed, so the node
	 * is real; what is gone is the layer, and the name went with it.  Say
	 * exactly that: "deleted layer" was a guess dressed as a fact, and
	 * "Item 2" tells the user nothing they can act on. */
	*kind = NDE_NODE_UNKNOWN;
	*label = g_strdup_printf(_("Layer %d, no longer in the document"), item_id);
	*orphan = TRUE;
}

static nde_graph_node *node_ensure(nde_graph *g, GHashTable *by_id, gint item_id) {
	nde_graph_node *n = g_hash_table_lookup(by_id, GINT_TO_POINTER(item_id));
	if (n)
		return n;
	n = g_new0(nde_graph_node, 1);
	n->item_id = item_id;
	n->record_ids = g_ptr_array_new_with_free_func(g_free);
	describe_item(item_id, &n->kind, &n->label, &n->orphan);
	/* Appended in first-record order, which is what keeps node order stable
	 * as the history grows. */
	g_ptr_array_add(g->nodes, n);
	g_hash_table_insert(by_id, GINT_TO_POINTER(item_id), n);
	return n;
}

/* AT ITEM GRANULARITY THE GRAPH HAS CYCLES, and a layout has to survive that.
 * A mask derives from the image and the image then consumes the mask, so the
 * two nodes point at each other.  Nothing is wrong: the cycle is an artefact
 * of collapsing records into items, and it vanishes at record granularity
 * (the mask reads the image's BASELINE and feeds a LATER record of it).  A
 * longest-path level would simply never terminate on that.
 *
 * The log breaks the tie, which is the design's own thesis — the flat history
 * is already a valid topological order.  Rank each node by where its FIRST
 * record sits, then compute levels over the forward edges only: those whose
 * source ranks earlier than their destination.  Ranks strictly increase along
 * such edges, so that subset is acyclic by construction.  The remaining edges
 * are still real and still returned; they just do not get a vote on layout.
 */
static void assign_levels(nde_graph *g, GHashTable *by_id) {
	GHashTable *rank = g_hash_table_new(g_direct_hash, g_direct_equal);
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		/* Nodes were appended in first-record order, so the index IS the rank. */
		g_hash_table_insert(rank, GINT_TO_POINTER(n->item_id), GINT_TO_POINTER(i));
	}

	gboolean changed = TRUE;
	guint guard = 0;
	while (changed && guard++ <= g->nodes->len) {
		changed = FALSE;
		for (guint i = 0; i < g->edges->len; i++) {
			const nde_graph_edge *e = g_ptr_array_index(g->edges, i);
			nde_graph_node *src = g_hash_table_lookup(by_id, GINT_TO_POINTER(e->src_item_id));
			nde_graph_node *dst = g_hash_table_lookup(by_id, GINT_TO_POINTER(e->dst_item_id));
			if (!src || !dst || src == dst)
				continue;
			gint rs = GPOINTER_TO_INT(g_hash_table_lookup(rank, GINT_TO_POINTER(src->item_id)));
			gint rd = GPOINTER_TO_INT(g_hash_table_lookup(rank, GINT_TO_POINTER(dst->item_id)));
			if (rs >= rd)
				continue;   /* a feedback edge — drawn, but not laid out on */
			if (dst->level < src->level + 1) {
				dst->level = src->level + 1;
				changed = TRUE;
			}
		}
	}
	g_hash_table_destroy(rank);
	g->depth = 1;
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->level + 1 > g->depth)
			g->depth = n->level + 1;
	}
}

nde_graph *nde_graph_build(void) {
	nde_graph *g = g_new0(nde_graph, 1);
	g->nodes = g_ptr_array_new_with_free_func(node_free);
	g->edges = g_ptr_array_new_with_free_func(edge_free);
	g->depth = 1;

	GPtrArray *snap = nde_history_snapshot(NULL);
	if (!snap)
		return g;

	GHashTable *by_id = g_hash_table_new(g_direct_hash, g_direct_equal);

	for (guint i = 0; i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		nde_graph_node *n = node_ensure(g, by_id, rec->target_item_id);
		gint64 *rid = g_new(gint64, 1);
		*rid = rec->record_id;
		g_ptr_array_add(n->record_ids, rid);

		if (!rec->inputs)
			continue;
		for (guint p = 0; p < rec->inputs->len; p++) {
			const nde_input_pin *pin = g_ptr_array_index(rec->inputs, p);
			/* A pin to the node's own chain is the implicit image edge
			 * restated, not a second node: drawing it would put a
			 * self-loop on every mask-creation step. */
			if (pin->src_item_id == rec->target_item_id)
				continue;
			node_ensure(g, by_id, pin->src_item_id);
			nde_graph_edge *e = g_new0(nde_graph_edge, 1);
			e->src_item_id   = pin->src_item_id;
			e->src_record_id = pin->src_record_id;
			e->dst_item_id   = rec->target_item_id;
			e->dst_record_id = rec->record_id;
			e->role          = g_strdup(pin->role ? pin->role : "");
			g_ptr_array_add(g->edges, e);
		}
	}

	/* An item the document no longer claims was usually consumed by a merge
	 * or flatten — and that operation recorded its NAME on the way past.  Use
	 * it: the layer is gone, but what it was called is not, and "Layer 2" is
	 * of no use to anyone.  (Once step 7 retains merge inputs as real items
	 * this fallback stops firing for merges, because flis_item_lookup will
	 * find them again.) */
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->kind != NDE_NODE_UNKNOWN)
			continue;
		for (guint r = 0; r < snap->len; r++) {
			const nde_record *rec = g_ptr_array_index(snap, r);
			if (!rec->params || g_strcmp0(rec->op_id, "layer.merge_down"))
				continue;
			GHashTable *kv = nde_kv_parse(rec->params);
			gint64 top = 0;
			const char *name = NULL;
			if (nde_kv_get_int(kv, "top_item", &top) && top == n->item_id)
				name = nde_kv_get_str(kv, "top_name");
			if (name && *name) {
				g_free(n->label);
				n->label = g_strdup_printf(_("%s, merged away"), name);
			}
			g_hash_table_unref(kv);
			if (name && *name)
				break;
		}
	}

	assign_levels(g, by_id);
	g_hash_table_destroy(by_id);
	g_ptr_array_unref(snap);
	return g;
}

const nde_graph_node *nde_graph_node_for(const nde_graph *g, gint item_id) {
	if (!g)
		return NULL;
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->item_id == item_id)
			return n;
	}
	return NULL;
}

static GPtrArray *edges_matching(const nde_graph *g, gint item_id, gboolean as_source) {
	GPtrArray *out = g_ptr_array_new();
	if (!g)
		return out;
	for (guint i = 0; i < g->edges->len; i++) {
		nde_graph_edge *e = g_ptr_array_index(g->edges, i);
		gint side = as_source ? e->src_item_id : e->dst_item_id;
		if (side == item_id)
			g_ptr_array_add(out, e);
	}
	return out;
}

GPtrArray *nde_graph_inputs_of(const nde_graph *g, gint item_id) {
	return edges_matching(g, item_id, FALSE);
}

GPtrArray *nde_graph_consumers_of(const nde_graph *g, gint item_id) {
	return edges_matching(g, item_id, TRUE);
}
