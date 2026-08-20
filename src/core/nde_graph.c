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
#include "core/nde_op_class.h"
#include "core/nde_graph.h"
#include "core/nde_composite.h"
#include "core/nde_joint.h"
#include "io/image_format_flis.h"

static void node_free(gpointer p) {
	nde_graph_node *n = p;
	if (!n)
		return;
	g_free(n->label);
	g_free(n->span_items);
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

/* A node for @item_id's steps in chronological @stage.  The item's FIRST
 * node keeps the real item id (edges, satellites and by-id lookups find it);
 * a later SEGMENT — its steps below a joint band — takes a fresh pseudo id
 * and inherits the first node's identity for labelling. */
static nde_graph_node *node_new_for_item(nde_graph *g, GHashTable *by_id,
                                         GHashTable *first_of, gint item_id,
                                         gint stage, gint *next_seg) {
	nde_graph_node *first = g_hash_table_lookup(first_of, GINT_TO_POINTER(item_id));
	nde_graph_node *n = g_new0(nde_graph_node, 1);
	n->real_item = item_id;
	n->stage = stage;
	n->record_ids = g_ptr_array_new_with_free_func(g_free);
	if (!first) {
		n->item_id = item_id;
		describe_item(item_id, &n->kind, &n->label, &n->orphan, &n->owner_item);
		g_hash_table_insert(first_of, GINT_TO_POINTER(item_id), n);
	} else {
		n->item_id = (*next_seg)--;
		n->kind = first->kind;
		n->label = g_strdup(first->label);
		n->orphan = first->orphan;
		n->owner_item = first->owner_item;
	}
	/* Appended in first-record order, which is what keeps node order stable
	 * as the history grows. */
	g_ptr_array_add(g->nodes, n);
	g_hash_table_insert(by_id, GINT_TO_POINTER(n->item_id), n);
	return n;
}

/* An edge as the scan found it, before its endpoints are bound to nodes:
 * with segments an item may have several, and the edge belongs to the one
 * holding the referenced record. */
typedef struct {
	gint    src_item;
	gint64  src_rec;
	gint    dst_node;    /* already a node id — the consuming node existed */
	gint64  dst_rec;
	gchar  *role;
} raw_edge;

static nde_graph_node *node_holding_record(nde_graph *g, gint item_id,
                                           gint64 record_id) {
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->real_item != item_id || n->spanning)
			continue;
		for (guint r = 0; r < n->record_ids->len; r++)
			if (*(gint64 *)g_ptr_array_index(n->record_ids, r) == record_id)
				return n;
	}
	return NULL;
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

	/* Joint records are chronological boundaries (nde_graph.h): the final
	 * level orders by STAGE first, edge level second, densified so empty
	 * ranks cost nothing.  With no joint records every node has stage 0 and
	 * the levels are the pure edge levels, unchanged. */
	GArray *pairs = g_array_new(FALSE, FALSE, sizeof(gint64));
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		gint64 key = ((gint64)n->stage << 32) | (guint32)n->level;
		gboolean seen = FALSE;
		for (guint j = 0; j < pairs->len && !seen; j++)
			seen = g_array_index(pairs, gint64, j) == key;
		if (!seen)
			g_array_append_val(pairs, key);
	}
	for (guint i = 1; i < pairs->len; i++) {   /* insertion sort: pairs are few */
		gint64 k = g_array_index(pairs, gint64, i);
		gint j = (gint)i - 1;
		while (j >= 0 && g_array_index(pairs, gint64, j) > k) {
			g_array_index(pairs, gint64, j + 1) = g_array_index(pairs, gint64, j);
			j--;
		}
		g_array_index(pairs, gint64, j + 1) = k;
	}
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		gint64 key = ((gint64)n->stage << 32) | (guint32)n->level;
		for (guint j = 0; j < pairs->len; j++)
			if (g_array_index(pairs, gint64, j) == key) {
				n->level = (gint)j;
				break;
			}
	}
	g_array_unref(pairs);

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
	GHashTable *first_of = g_hash_table_new(g_direct_hash, g_direct_equal);
	/* item -> its node in the CURRENT chronological stage; cleared at every
	 * joint boundary so the next record of the item opens a segment below. */
	GHashTable *cur_of = g_hash_table_new(g_direct_hash, g_direct_equal);
	GArray *raw = g_array_new(FALSE, TRUE, sizeof(raw_edge));
	gint stage = 0;
	gint next_seg = NDE_GRAPH_SEG_ITEM_BASE;

	for (guint i = 0; i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		/* A JOINT record (nde_joint.h) is one operation across several
		 * layers, and it is drawn like one: a single band SPANNING the
		 * graph, no edges — the layers above it are what it read, the
		 * segments below are what follows it.  It is a chronological
		 * boundary: every record after it, on any item, draws below. */
		if (nde_joint_is_op(rec->op_id)) {
			nde_graph_node *jn = g_new0(nde_graph_node, 1);
			jn->item_id = nde_graph_joint_pseudo_item(rec->record_id);
			jn->kind = NDE_NODE_JOINT;
			jn->spanning = TRUE;
			/* The band spans its PARTICIPANTS' columns, not the whole
			 * graph — a layer the operation never read stays visually
			 * outside it. */
			jn->span_items = nde_joint_record_participants(rec, &jn->n_span_items);
			jn->real_item = rec->target_item_id;
			jn->stage = ++stage;
			jn->label = g_strdup(rec->summary && *rec->summary ?
			                     rec->summary : _("Joint calibration"));
			jn->record_ids = g_ptr_array_new_with_free_func(g_free);
			jn->orphan = TRUE;   /* the Layers panel lists no such item */
			gint64 *jrid = g_new(gint64, 1);
			*jrid = rec->record_id;
			g_ptr_array_add(jn->record_ids, jrid);
			g_ptr_array_add(g->nodes, jn);
			g_hash_table_insert(by_id, GINT_TO_POINTER(jn->item_id), jn);
			stage++;
			g_hash_table_remove_all(cur_of);
			continue;
		}
		nde_graph_node *n = g_hash_table_lookup(cur_of,
				GINT_TO_POINTER(rec->target_item_id));
		if (!n) {
			n = node_new_for_item(g, by_id, first_of, rec->target_item_id,
			                      stage, &next_seg);
			g_hash_table_insert(cur_of, GINT_TO_POINTER(rec->target_item_id), n);
		}
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
			raw_edge re = {
				.src_item = pin->src_item_id,
				.src_rec  = pin->src_record_id,
				.dst_node = n->item_id,
				.dst_rec  = rec->record_id,
				.role     = g_strdup(pin->role ? pin->role : ""),
			};
			g_array_append_val(raw, re);
		}
	}

	/* Bind the collected edges to nodes: the source is the segment holding
	 * the pinned record (a pin into an earlier segment must not point at the
	 * continuation), falling back to the item's first node — created here
	 * when the item has no records at all (an input-only source; rank -1 in
	 * assign_levels puts it above everything that derives from it).
	 *
	 * A pin naming a JOINT record (a composite consuming a participant's
	 * calibrated state — nde_history_last_record_for_item points there now)
	 * binds to the BAND, aligned on the participant's column: the edge
	 * leaves the band's bottom under that channel. */
	for (guint i = 0; i < raw->len; i++) {
		raw_edge *re = &g_array_index(raw, raw_edge, i);
		gint align = 0;
		nde_graph_node *src = NULL;
		if (re->src_rec) {
			for (guint j = 0; j < g->nodes->len && !src; j++) {
				nde_graph_node *n = g_ptr_array_index(g->nodes, j);
				if (!n->spanning)
					continue;
				for (guint r = 0; r < n->record_ids->len && !src; r++) {
					if (*(gint64 *)g_ptr_array_index(n->record_ids, r)
					    == re->src_rec) {
						src = n;
						align = re->src_item;
					}
				}
			}
			if (!src)
				src = node_holding_record(g, re->src_item, re->src_rec);
		}
		if (!src)
			src = g_hash_table_lookup(first_of, GINT_TO_POINTER(re->src_item));
		if (!src)
			src = node_new_for_item(g, by_id, first_of, re->src_item, 0,
			                        &next_seg);
		nde_graph_edge *e = g_new0(nde_graph_edge, 1);
		e->src_item_id    = src->item_id;
		e->src_record_id  = re->src_rec;
		e->dst_item_id    = re->dst_node;
		e->dst_record_id  = re->dst_rec;
		e->src_align_item = align;
		e->role           = re->role;   /* ownership moves */
		g_ptr_array_add(g->edges, e);
	}
	g_array_unref(raw);
	g_hash_table_destroy(cur_of);
	g_hash_table_destroy(first_of);

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
				if (in->mask_item_id == n->real_item) {
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
				if (in->item_id != n->real_item)
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
				if (rec->target_item_id == n->real_item &&
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

	/* DELETED, as opposed to merely gone.  Both a deleted layer and one
	 * consumed by a flatten leave the document, so both describe_item() as
	 * NDE_NODE_UNKNOWN — but they are not the same thing to a reader.  A
	 * consumed layer is where the surviving image CAME FROM and belongs among
	 * the live columns with the composite below it; a deleted one is a dead
	 * end, kept only so the deletion can be undone.  The record says which:
	 * NDE_OPT_DELETES_ITEM is written for a deletion and for nothing else. */
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		for (guint r = 0; r < snap->len && !n->deleted; r++) {
			const nde_record *rec = g_ptr_array_index(snap, r);
			n->deleted = (nde_op_class_for(rec->op_id)->traits & NDE_OPT_DELETES_ITEM) &&
			             rec->target_item_id == n->real_item;
		}
	}

	/* Which item's COLUMN each node belongs in.  Two cases, one field:
	 *
	 *  - a SEGMENT continues its own item's history below a joint band, and
	 *    its pseudo id hides that;
	 *  - a node DERIVED from others — a flatten or merge result — reads as
	 *    the place its inputs arrive at, so it belongs under them rather than
	 *    off in a column of its own.  Its inputs' columns are ordered as the
	 *    nodes are, so the first input in node order is the leftmost.
	 *
	 * Both are "this box does not own a fresh column", which is exactly what
	 * the layout needs to know; deciding it here keeps the policy next to the
	 * edges it is derived from, and testable without a widget. */
	for (guint i = 0; i < g->nodes->len; i++) {
		nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		if (n->spanning)
			continue;
		if (n->real_item != n->item_id) {
			n->column_item = n->real_item;
			continue;
		}
		for (guint j = 0; j < g->nodes->len && !n->column_item; j++) {
			const nde_graph_node *cand = g_ptr_array_index(g->nodes, j);
			if (cand == n || cand->spanning || cand->level >= n->level)
				continue;
			for (guint e = 0; e < g->edges->len; e++) {
				const nde_graph_edge *ed = g_ptr_array_index(g->edges, e);
				if (ed->dst_item_id != n->item_id)
					continue;
				/* An input reached THROUGH a joint band still arrives from a
				 * column, and the edge already says which: src_align_item is
				 * the participant it leaves under.  Without this the base
				 * layer of a flatten is invisible as an input — its edge
				 * comes from the band, not from itself — and the result
				 * anchors to whichever input happens to have a direct edge,
				 * which is the far right of the group rather than its base. */
				const nde_graph_node *src =
						nde_graph_node_for(g, ed->src_item_id);
				const gint from = (src && src->spanning && ed->src_align_item)
						? ed->src_align_item : ed->src_item_id;
				if (from == cand->item_id) {
					n->column_item = cand->item_id;
					break;
				}
			}
		}
	}

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

gint nde_graph_route_for_record(const nde_graph *g, gint64 record_id) {
	if (!g)
		return 0;
	for (guint i = 0; i < g->nodes->len; i++) {
		const nde_graph_node *n = g_ptr_array_index(g->nodes, i);
		for (guint r = 0; r < n->record_ids->len; r++)
			if (*(gint64 *)g_ptr_array_index(n->record_ids, r) == record_id)
				return n->item_id;
	}
	return 0;
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

/* One item's hold on a horizontal range, for the whole height of the layout.
 * @owner is the ITEM, not the box: a segment continues its anchor's column and
 * so shares the claim, and a host's satellite column is part of the same one. */
typedef struct { gint x0, x1, owner; } column_claim;

static void claim_add(GArray *cl, gint x0, gint x1, gint owner) {
	column_claim k = { x0, x1, owner };
	g_array_append_val(cl, k);
}

/* The first x at or after @from where a box of width @w belonging to @owner
 * fits without landing on somebody else's column.  Claims number one or two
 * per node and this runs once per box, so the repeated scan is not worth a
 * cleverer structure; it re-scans after each move because stepping past one
 * claim can land on another. */
static gint claim_free_x(const GArray *cl, gint from, gint w, gint owner,
                         gint gap) {
	gint x = from;
	gboolean moved = TRUE;
	while (moved) {
		moved = FALSE;
		for (guint c = 0; c < cl->len; c++) {
			const column_claim *k = &g_array_index(cl, column_claim, c);
			if (k->owner == owner || x >= k->x1 || x + w <= k->x0)
				continue;
			x = k->x1 + gap;
			moved = TRUE;
		}
	}
	return x;
}

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
	gboolean *band_span = g_new0(gboolean, (gsize)n_bands);
	for (guint i = 0; i < n; i++) {
		if (host_of[i] >= 0)
			continue;   /* counted through its host's column below */
		const gint own = MAX(BOX(i)->h, col_h[i]);
		if (own > band_h[level_of[i]])
			band_h[level_of[i]] = own;
		if (BOX(i)->spanning)
			band_span[level_of[i]] = TRUE;
	}
	/* A spanning band keeps the COLUMN gap as its vertical whitespace —
	 * enough to read as its own stage without opening a chasm.  Ordinary
	 * band boundaries keep the row gap. */
	gint *band_y = g_new0(gint, (gsize)n_bands);
	gint y = 0;
	gint last_gap = 0;
	for (gint r = 0; r < n_bands; r++) {
		band_y[r] = y;
		const gboolean at_span = band_span[r] ||
				(r + 1 < n_bands && band_span[r + 1]);
		last_gap = at_span ? col_gap : row_gap;
		y += band_h[r] + last_gap;
	}

	/* Resolve each segment to the box whose column it continues.  The anchor
	 * has to be an ordinary box in an EARLIER band — a segment is a later
	 * chronological stage than the node it continues, so that always holds for
	 * a real one, and demanding it keeps a malformed input from making the
	 * placement order circular.  Unresolved boxes flow normally. */
	gint *anchor_of = g_new0(gint, n);
	for (guint i = 0; i < n; i++) {
		anchor_of[i] = -1;
		if (!BOX(i)->align_item || host_of[i] >= 0 || BOX(i)->spanning)
			continue;
		for (guint j = 0; j < n; j++) {
			if (j == i || host_of[j] >= 0 || BOX(j)->spanning ||
			    BOX(j)->align_item)
				continue;
			if (BOX(j)->item_id == BOX(i)->align_item &&
			    level_of[j] < level_of[i]) {
				anchor_of[i] = (gint)j;
				break;
			}
		}
	}

	/* A COLUMN BELONGS TO THE ITEM THAT CLAIMED IT, for the whole height of
	 * the layout.  Bands used to fill independently, which meant an item whose
	 * history starts lower down was left-aligned into a column another item
	 * already occupied above — a layer added after a joint band was drawn
	 * under Background and read as part of its history.  A column is one
	 * item's story down the page, so a claim is a claim.
	 *
	 * Claims are keyed by the item, not the box: a segment continues its
	 * anchor's column and so shares the claim, while everyone else is kept
	 * out.  A host's satellite column is part of the same claim — a later node
	 * walking left to right has to clear the whole group. */
	GArray *claims = g_array_new(FALSE, FALSE, sizeof(column_claim));

	/* Place hosts left to right, each followed by its own column, so that the
	 * band's next node starts beyond the satellites rather than before them.
	 * A spanning box takes no x slot: its width is the whole layout's, set
	 * once that width is known below.
	 *
	 * Band by band, because a segment takes its x from a box in an earlier
	 * band and so cannot be placed until that one has been.  Within a band the
	 * segments go FIRST: they have no choice of column, so the ordinary nodes
	 * — which do — are the ones that must give way.
	 *
	 * Then the whole thing twice: live boxes, then RETIRED ones.  A removed
	 * layer's node is real and still listed, but it is no longer part of the
	 * live story, so it belongs to the right of every live column instead of
	 * between them.  That cannot be decided band by band — a retired node in
	 * the top band has to clear a live one three bands down — so it takes a
	 * second pass, once the live extent is known. */
	/* A segment shares its anchor's column, so it shares its anchor's PASS
	 * too: a removed layer's continuation goes right with it rather than
	 * being placed among the live columns while the anchor is still
	 * unplaced. */
	gboolean *retired_of = g_new0(gboolean, n);
	for (guint i = 0; i < n; i++)
		retired_of[i] = anchor_of[i] >= 0 ? BOX(anchor_of[i])->retired
		                                  : BOX(i)->retired;

	gint *band_x = g_new0(gint, (gsize)n_bands);
	gint *px = g_new0(gint, n);
	gint live_right = 0;
	for (int pass = 0; pass < 2; pass++) {
		const gboolean want_retired = (pass == 1);
		if (want_retired)
			for (guint c = 0; c < claims->len; c++)
				live_right = MAX(live_right,
				                 g_array_index(claims, column_claim, c).x1 + col_gap);
		for (gint r = 0; r < n_bands; r++) {
			for (guint i = 0; i < n; i++) {
				if (level_of[i] != r || anchor_of[i] < 0 ||
				    (!retired_of[i]) == want_retired)
					continue;
				/* No choice of column: it IS the anchor's. */
				px[i] = px[anchor_of[i]];
				const gint owner = BOX(anchor_of[i])->item_id;
				gint end = px[i] + BOX(i)->w;
				claim_add(claims, px[i], end, owner);
				end += col_gap;
				if (col_w[i]) {
					for (guint j = 0; j < n; j++)
						if (host_of[j] == (gint)i)
							px[j] = end;
					claim_add(claims, end, end + col_w[i], owner);
					end += col_w[i] + col_gap;
				}
				band_x[r] = MAX(band_x[r], end);
			}
			for (guint i = 0; i < n; i++) {
				if (level_of[i] != r || host_of[i] >= 0 || BOX(i)->spanning ||
				    anchor_of[i] >= 0 || (!retired_of[i]) == want_retired)
					continue;
				const gint owner = BOX(i)->item_id;
				const gint floor_x = MAX(band_x[r], want_retired ? live_right : 0);
				const gint total = col_w[i] ? BOX(i)->w + col_gap + col_w[i]
				                            : BOX(i)->w;
				px[i] = claim_free_x(claims, floor_x, total, owner, col_gap);
				claim_add(claims, px[i], px[i] + BOX(i)->w, owner);
				band_x[r] = px[i] + BOX(i)->w + col_gap;
				if (col_w[i]) {
					for (guint j = 0; j < n; j++)
						if (host_of[j] == (gint)i)
							px[j] = band_x[r];
					claim_add(claims, band_x[r], band_x[r] + col_w[i], owner);
					band_x[r] += col_w[i] + col_gap;
				}
			}
		}
	}
	g_array_unref(claims);

	gint used_w = 0;
	for (gint r = 0; r < n_bands; r++)
		if (band_x[r] - col_gap > used_w)
			used_w = band_x[r] - col_gap;

	for (guint i = 0; i < n; i++) {
		const nde_graph_box *b = BOX(i);
		const gint top = band_y[level_of[i]];
		gint bx = px[i], bw = b->w;
		if (b->spanning) {
			/* Stretch over the participants' columns — the layers the
			 * operation actually read.  When none of them is present
			 * (or no list was given), fall back to the full width. */
			gint sx0 = G_MAXINT, sx1 = G_MININT;
			for (guint j = 0; j < n; j++) {
				if (j == i || host_of[j] >= 0 || BOX(j)->spanning)
					continue;
				gboolean member = FALSE;
				for (guint k = 0; k < b->n_span_items && !member; k++)
					member = b->span_items[k] == BOX(j)->item_id;
				if (!member)
					continue;
				if (px[j] < sx0)
					sx0 = px[j];
				if (px[j] + BOX(j)->w > sx1)
					sx1 = px[j] + BOX(j)->w;
			}
			if (sx1 > sx0) {
				bx = sx0;
				bw = sx1 - sx0;
			} else {
				bx = 0;
				bw = MAX(used_w, b->w);
			}
		}
		nde_graph_place p = {
			.item_id   = b->item_id,
			.x         = bx,
			.w         = bw,
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
	if (total_h) *total_h = y > 0 ? y - last_gap : 0;
	g_free(host_of);
	g_free(anchor_of);
	g_free(retired_of);
	g_free(level_of);
	g_free(col_w);
	g_free(col_h);
	g_free(sat_y);
	g_free(sat_tail);
	g_free(band_h);
	g_free(band_span);
	g_free(band_x);
	g_free(band_y);
	g_free(px);
	#undef BOX
	return out;
}
