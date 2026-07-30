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
#include "core/nde_composite.h"
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
                          gboolean *orphan, gint *owner_item) {
	*owner_item = 0;
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
		*owner_item = NDE_ITEM_IMAGE;
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
		*owner_item = owner ? owner->item_id : 0;
		return;
	case FLIS_ITEM_LMASK:
		*kind = NDE_NODE_LAYERMASK;
		*label = owner && owner->layer_name
		         ? g_strdup_printf(_("Layer mask of %s"), owner->layer_name)
		         : g_strdup(_("Layer mask"));
		*orphan = TRUE;
		*owner_item = owner ? owner->item_id : 0;
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
	describe_item(item_id, &n->kind, &n->label, &n->orphan, &n->owner_item);
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
		/* Nodes were appended in first-record order, so the index IS the rank
		 * — except for a node with NO records, which has no log position to
		 * rank by.  Discovery order appended it AFTER the record that named
		 * it, which misread its derivation edge as feedback: a mask's source
		 * image sat BESIDE the mask, and the arrow between them was drawn
		 * backwards into empty space.  Such a node is a pure source (every
		 * edge destination is some record's target, so nothing can point into
		 * it), so ranking it before the whole log is acyclic by construction
		 * and puts what derives from it below it, where derivation reads. */
		gint r = n->record_ids->len ? (gint)i : -1;
		g_hash_table_insert(rank, GINT_TO_POINTER(n->item_id), GINT_TO_POINTER(r));
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
		gboolean named = FALSE;
		for (guint r = 0; r < snap->len && !named; r++) {
			const nde_record *rec = g_ptr_array_index(snap, r);
			if (!nde_composite_is_op(rec->op_id))
				continue;
			nde_composite_state *st = nde_composite_state_parse(rec->params);
			if (!st)
				continue;
			for (guint k = 0; k < st->inputs->len && !named; k++) {
				const nde_composite_input *in =
						&g_array_index(st->inputs, nde_composite_input, k);
				if (!in->name || !*in->name)
					continue;
				/* A composite records its inputs' LAYER MASK items too, and
				 * one of those going unnamed is worse than an unnamed layer:
				 * "Layer 3, no longer in the document" describes a layer the
				 * user never had. */
				if (in->mask_item_id == n->item_id) {
					g_free(n->label);
					n->label = g_strdup_printf(_("Layer mask of %s, no longer in the document"),
					                           in->name);
					n->kind = NDE_NODE_LAYERMASK;
					/* The layer is a retained input and still a node, so the
					 * mask can still be drawn beside the layer it masked. */
					n->owner_item = in->item_id;
					named = TRUE;
					break;
				}
				if (in->item_id != n->item_id)
					continue;
				g_free(n->label);
				n->label = g_strdup_printf(st->raw_first ? _("%s, merged away")
				                                         : _("%s, flattened"),
				                           in->name);
				named = TRUE;
			}
			nde_composite_state_free(st);
		}
	}

	/* A node still unknown after the composite pass, but that some record
	 * consumed through a "mask" pin — or whose own records are all mask ops —
	 * was a processing or layer mask whose slot has since been cleared.
	 * "Layer N, no longer in the document" describes a layer the user never
	 * had; say what it was, and whose.  Not "no longer in the document",
	 * either: a used mask goes DORMANT when cleared — its records and the
	 * pinned states its consumers read all persist, and editing them still
	 * works (nde_replay.c) — so the label says only that it was cleared. */
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->kind != NDE_NODE_UNKNOWN)
			continue;
		const nde_graph_node *consumer = NULL;
		gboolean was_mask = FALSE;
		for (guint e = 0; e < g->edges->len && !consumer; e++) {
			const nde_graph_edge *ed = g_ptr_array_index(g->edges, e);
			if (ed->src_item_id == n->item_id && !g_strcmp0(ed->role, "mask")) {
				was_mask = TRUE;
				consumer = g_hash_table_lookup(by_id,
				                               GINT_TO_POINTER(ed->dst_item_id));
			}
		}
		if (!was_mask && n->record_ids->len && snap) {
			was_mask = TRUE;
			for (guint r = 0; r < snap->len && was_mask; r++) {
				const nde_record *rec = g_ptr_array_index(snap, r);
				if (rec->target_item_id == n->item_id &&
				    !g_str_has_prefix(rec->op_id, "mask."))
					was_mask = FALSE;
			}
		}
		if (!was_mask)
			continue;
		g_free(n->label);
		n->label = (consumer && consumer->label)
				? g_strdup_printf(_("Mask of %s (cleared)"), consumer->label)
				: g_strdup(_("Mask (cleared)"));
		n->kind = NDE_NODE_MASK;
		n->owner_item = consumer ? consumer->item_id : 0;
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

gboolean nde_graph_edge_is_feedback(const nde_graph *g, const nde_graph_edge *e) {
	if (!g || !e)
		return FALSE;
	const nde_graph_node *src = nde_graph_node_for(g, e->src_item_id);
	const nde_graph_node *dst = nde_graph_node_for(g, e->dst_item_id);
	if (!src || !dst)
		return FALSE;
	/* Levels were assigned over the forward edges only, so an edge that does
	 * not increase the level is exactly one that was left out of them. */
	return dst->level <= src->level;
}

/* ---- layout ------------------------------------------------------------- */

GArray *nde_graph_layout(const GArray *boxes, gint col_gap, gint row_gap,
                         gint *total_w, gint *total_h) {
	GArray *out = g_array_new(FALSE, TRUE, sizeof(nde_graph_place));
	if (total_w) *total_w = 0;
	if (total_h) *total_h = 0;
	if (!boxes || !boxes->len)
		return out;

	const guint n = boxes->len;
	#define BOX(i) (&g_array_index(boxes, nde_graph_box, (i)))

	/* Resolve each satellite to its host's INDEX.  A host that is not here, or
	 * that is a satellite itself, leaves the box hostless: it then lays out as
	 * an ordinary node at its own level, which is wrong-looking but visible —
	 * and visible beats unplaced. */
	gint *host_of = g_new0(gint, n);
	for (guint i = 0; i < n; i++) {
		host_of[i] = -1;
		if (!BOX(i)->host_item)
			continue;
		for (guint j = 0; j < n; j++) {
			if (j != i && !BOX(j)->host_item &&
			    BOX(j)->item_id == BOX(i)->host_item) {
				host_of[i] = (gint)j;
				break;
			}
		}
	}
	/* A satellite belongs to its host's band, whatever level it was given. */
	gint *level_of = g_new0(gint, n);
	gint n_bands = 0;
	for (guint i = 0; i < n; i++) {
		level_of[i] = host_of[i] >= 0 ? BOX(host_of[i])->level : BOX(i)->level;
		if (level_of[i] + 1 > n_bands)
			n_bands = level_of[i] + 1;
	}

	/* Each host's satellite column, indexed by the HOST's box index: how wide
	 * it is, how tall the stack comes to, and — per satellite — how much of
	 * that stack still follows it, which is what bounds how far down the
	 * satellite may be nudged to meet its anchor row. */
	gint *col_w = g_new0(gint, n);
	gint *col_h = g_new0(gint, n);
	gint *sat_y = g_new0(gint, n);   /* offset into its host's column */
	gint *sat_tail = g_new0(gint, n);
	for (guint i = 0; i < n; i++) {
		gint hi = host_of[i];
		if (hi < 0)
			continue;
		if (BOX(i)->w > col_w[hi])
			col_w[hi] = BOX(i)->w;
		/* Satellites are a column, so what separates them vertically is the
		 * column gap, not the gap between bands. */
		sat_y[i] = col_h[hi] ? col_h[hi] + col_gap : 0;
		col_h[hi] = sat_y[i] + BOX(i)->h;
	}
	for (guint i = 0; i < n; i++) {
		if (host_of[i] >= 0)
			sat_tail[i] = col_h[host_of[i]] - (sat_y[i] + BOX(i)->h);
	}

	/* A level with no node in it still costs a band: the gap says the
	 * derivation skipped a rank rather than that two nodes are siblings. */
	gint *band_h = g_new0(gint, (gsize)n_bands);
	for (guint i = 0; i < n; i++) {
		if (host_of[i] >= 0)
			continue;   /* counted through its host's column below */
		const gint own = MAX(BOX(i)->h, col_h[i]);
		if (own > band_h[level_of[i]])
			band_h[level_of[i]] = own;
	}
	gint *band_y = g_new0(gint, (gsize)n_bands);
	gint y = 0;
	for (gint r = 0; r < n_bands; r++) {
		band_y[r] = y;
		y += band_h[r] + row_gap;
	}

	/* Place hosts left to right, each followed by its own column, so that the
	 * band's next node starts beyond the satellites rather than before them. */
	gint *band_x = g_new0(gint, (gsize)n_bands);
	gint *px = g_new0(gint, n);
	for (guint i = 0; i < n; i++) {
		if (host_of[i] >= 0)
			continue;
		const gint lvl = level_of[i];
		px[i] = band_x[lvl];
		band_x[lvl] += BOX(i)->w + col_gap;
		if (col_w[i]) {
			for (guint j = 0; j < n; j++)
				if (host_of[j] == (gint)i)
					px[j] = band_x[lvl];
			band_x[lvl] += col_w[i] + col_gap;
		}
	}

	gint used_w = 0;
	for (gint r = 0; r < n_bands; r++)
		if (band_x[r] - col_gap > used_w)
			used_w = band_x[r] - col_gap;

	for (guint i = 0; i < n; i++) {
		const nde_graph_box *b = BOX(i);
		const gint top = band_y[level_of[i]];
		nde_graph_place p = {
			.item_id   = b->item_id,
			.x         = px[i],
			.w         = b->w,
			.content_h = b->h,
		};
		if (host_of[i] >= 0) {
			p.y     = top + sat_y[i];
			p.h     = b->h;   /* nothing levels a satellite */
			p.y_min = p.y;
			p.y_max = MAX(p.y, top + band_h[level_of[i]] - b->h - sat_tail[i]);
		} else {
			p.y     = top;
			p.h     = band_h[level_of[i]];
			p.y_min = p.y_max = p.y;
		}
		g_array_append_val(out, p);
	}

	if (total_w) *total_w = used_w;
	if (total_h) *total_h = y > 0 ? y - row_gap : 0;
	g_free(host_of);
	g_free(level_of);
	g_free(col_w);
	g_free(col_h);
	g_free(sat_y);
	g_free(sat_tail);
	g_free(band_h);
	g_free(band_x);
	g_free(band_y);
	g_free(px);
	#undef BOX
	return out;
}
