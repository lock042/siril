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
#ifndef _NDE_GRAPH_H_
#define _NDE_GRAPH_H_

/**
 * \file nde_graph.h
 * \brief The history read as the DAG it has been since input pins landed.
 *
 * The log is flat and totally ordered; that ordering is a valid topological
 * order of a graph whose NODES are items (the image, its masks, each layer)
 * and whose EDGES are the input pins on records.  This module derives that
 * view.  It computes nothing new — every node and every edge is already in
 * the history — it just stops projecting the graph onto a single list.
 *
 * Why it exists as a separate, GUI-free module: the graph view is the core
 * History UI (design note §9.1), and the part worth testing is the shape,
 * not the widgets.  A node's records stay in log order, because ORDER IS
 * MEANING inside a node: position, not id, is what the user reorders.
 *
 * Threading: call from the main thread.  Builds from a history snapshot and
 * owns everything it returns.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** What an item is, for labelling and for deciding what may be done to it. */
typedef enum {
	NDE_NODE_IMAGE,       /* the plain single image, or a FLIS layer */
	NDE_NODE_LAYER,
	NDE_NODE_MASK,        /* an operation mask (fit->mask) */
	NDE_NODE_LAYERMASK,   /* a layer mask (layer->lmask) */
	NDE_NODE_UNKNOWN,     /* an item id nothing in the document claims */
	NDE_NODE_JOINT,       /* a JOINT multi-layer record (nde_joint.h): one
	                       * record shown as a SPANNING band across the whole
	                       * graph — a normal operation that happens to affect
	                       * every layer above it, so it draws no edges */
} nde_node_kind;

/* Joint records are not items, so their nodes carry PSEUDO item ids derived
 * from the record id, well below the real sentinels (-1 image, -2 its mask).
 * SEGMENT nodes (an item's steps AFTER a joint record — the joint is a
 * chronological boundary, so what follows it is drawn below it) draw their
 * pseudo ids sequentially from a second, deeper range. */
#define NDE_GRAPH_JOINT_ITEM_BASE (-1000)
#define NDE_GRAPH_SEG_ITEM_BASE   (-1000000)

static inline gint nde_graph_joint_pseudo_item(gint64 record_id) {
	return NDE_GRAPH_JOINT_ITEM_BASE - (gint)record_id;
}
static inline gboolean nde_graph_item_is_joint_node(gint item_id) {
	return item_id <= NDE_GRAPH_JOINT_ITEM_BASE && item_id > NDE_GRAPH_SEG_ITEM_BASE;
}
static inline gint64 nde_graph_joint_node_record(gint item_id) {
	return (gint64)(NDE_GRAPH_JOINT_ITEM_BASE - item_id);
}

typedef struct {
	gint            item_id;
	nde_node_kind   kind;
	gchar          *label;        /* owned; user-facing, already translated */
	GPtrArray      *record_ids;   /* gint64* (owned); THIS ITEM'S records, in log order */
	/**
	 * Layout depth, 0 = derives from nothing earlier.  Computed over the
	 * FORWARD edges only (source's first record precedes the destination's),
	 * because at item granularity the graph genuinely contains cycles — a
	 * mask derives from the image and the image then consumes the mask.  See
	 * assign_levels() in nde_graph.c.
	 */
	gint            level;
	/**
	 * The item this one BELONGS to, or 0.  A mask is owned by the layer whose
	 * pixels it modulates — the document says so directly, whatever the edges
	 * happen to be — and that is what lets a layer mask be drawn beside its
	 * layer even though the only edge it has runs to a composite further down.
	 */
	gint            owner_item;
	/**
	 * TRUE when the Layers panel does not list this item: masks always, and
	 * (from step 7) inputs retained behind a composite.  These are exactly
	 * the items the flat list had nowhere to put — it rendered them
	 * "deleted layer", which for a live, editable item is simply false.
	 */
	gboolean        orphan;
	/**
	 * The DOCUMENT item this node's steps belong to.  == item_id for an
	 * ordinary node; a SEGMENT (a continuation of an item's history below a
	 * joint band) carries a pseudo item_id and names its item here, and a
	 * joint node names its anchor.
	 */
	gint            real_item;
	/**
	 * Chronological stage: joint records are band boundaries — everything
	 * recorded after one draws below it, whatever the derivation edges say
	 * — and this is which inter-joint interval the node's steps fall in.
	 * The final @level is (stage, edge level) densified; with no joint
	 * records every node has stage 0 and levels are pure edge levels.
	 */
	gint            stage;
	/** TRUE for a joint node: laid out alone in its band, spanning its
	 *  participants' columns, with its own whitespace above and below. */
	gboolean        spanning;
	/** The participants a joint node spans (owned; NULL otherwise).  The
	 *  layout stretches the band over exactly these items' columns. */
	gint           *span_items;
	guint           n_span_items;
} nde_graph_node;

typedef struct {
	gint    src_item_id;
	gint64  src_record_id;   /* 0 = the source item's baseline */
	gint    dst_item_id;
	gint64  dst_record_id;   /* the consuming record */
	gchar  *role;            /* owned; "mask", later "base", "overlay"… */
	/**
	 * When the source is a JOINT band, the participant whose state this
	 * edge carries (0 otherwise).  The view draws the edge leaving the
	 * band's bottom at that participant's column, so a flatten consuming
	 * three calibrated channels shows three departures, each under its
	 * channel, instead of one line from the band's middle.
	 */
	gint    src_align_item;
} nde_graph_edge;

typedef struct {
	GPtrArray *nodes;   /* nde_graph_node* */
	GPtrArray *edges;   /* nde_graph_edge* */
	gint       depth;   /* 1 + the highest level present; 1 for a linear document */
} nde_graph;

/**
 * Build the graph for the current history.  Never NULL: a document with no
 * records yields an empty graph, and one with a single chain yields exactly
 * one node — which is why a linear document draws as today's list.
 *
 * Node order is stable: by the position of each item's FIRST record, so the
 * graph does not reshuffle as records are added.
 */
nde_graph *nde_graph_build(void);

void nde_graph_free(nde_graph *g);

/** The node for @item_id, or NULL. */
const nde_graph_node *nde_graph_node_for(const nde_graph *g, gint item_id);

/** The node whose steps include LIVE record @record_id — how the History
 *  panel routes a row when segments split an item across bands.  0 when no
 *  node lists it (a dead-tail record: route to its target item). */
gint nde_graph_route_for_record(const nde_graph *g, gint64 record_id);

/** Edges INTO @item_id (what it consumes).  Caller g_ptr_array_unref()s the
 *  array; the edges themselves stay owned by @g. */
GPtrArray *nde_graph_inputs_of(const nde_graph *g, gint item_id);

/** Edges OUT of @item_id (what consumes it) — the reverse-invalidation set
 *  made visible: these are the nodes an amend here would re-run. */
GPtrArray *nde_graph_consumers_of(const nde_graph *g, gint item_id);

/**
 * TRUE for an edge that runs against the levels — the destination sits at or
 * before the source.  These are the cycle at item granularity (§9.1.1): real,
 * drawable, and deliberately given no vote on where anything sits.
 */
gboolean nde_graph_edge_is_feedback(const nde_graph *g, const nde_graph_edge *e);

/* ---- layout -------------------------------------------------------------
 * Kept here, and kept a pure function of measured sizes, because the part
 * worth testing is the shape and not the widgets: the view measures its
 * children, calls this, and allocates what it gets back.
 */

/** One node's identity and measured size, as the caller found it. */
typedef struct {
	gint item_id;
	gint level;
	/**
	 * 0 for an ordinary band node.  Otherwise the item this box is a
	 * SATELLITE of: it is laid out immediately beside that box, in its band,
	 * and the band's next node starts to the right of the pair.  @level is
	 * then ignored — a satellite sits in its host's band by definition.
	 * A satellite whose host is absent, or is itself a satellite, degrades to
	 * an ordinary node at its own level rather than going unplaced.
	 */
	gint host_item;
	gint w, h;
	/**
	 * TRUE for a joint node's box: laid out alone in its band, stretched
	 * over its participants' columns (@span_items; the whole layout when
	 * none of them are present), with the COLUMN gap as its vertical
	 * whitespace — a normal operation that spans the channels it read.
	 */
	gboolean spanning;
	const gint *span_items;   /* borrowed for the duration of the call */
	guint       n_span_items;
} nde_graph_box;

/** Where that node goes.  Origin is the top-left of the whole layout. */
typedef struct {
	gint item_id;
	gint x, y, w;
	/**
	 * @h is the BAND's height for an ordinary node — every node in a band is
	 * levelled to it, so the band has one bottom edge and satellites have one
	 * range to be anchored in.  @content_h is what the node actually measured,
	 * and is where its outgoing edges leave: a short node in a tall band drew
	 * its line starting at the band's bottom, which read as a line that had
	 * come adrift from the end of its history.  For a satellite the two are
	 * equal — it is not levelled to anything.
	 */
	gint h, content_h;
	gint y_min, y_max;   /* the range a satellite's host allows it; else = y */
} nde_graph_place;

/**
 * Place @boxes in horizontal bands by level, left-aligned, in the order given
 * — which is nde_graph_build()'s stable node order, so the layout does not
 * reshuffle as records are added.  A band is as tall as its tallest node, so
 * every edge leaves a band at the same y.
 *
 * Derivation runs DOWNWARD and siblings sit side by side: three layers across
 * the top, the merge that consumed them below.  Steps already run top to
 * bottom inside a node, and running derivation left to right as well made the
 * two axes fight — a layer's history read down the page while the history of
 * the document read across it.
 *
 * A box with a host_item is a SATELLITE: it goes in a column immediately to
 * the right of its host, inside the host's band, and the band's next node
 * starts beyond that column.  Masks are laid out this way, so a mask sits
 * beside the layer it applies to instead of after every layer — where its
 * connecting lines had to cross the layers in between to get back.  Several
 * satellites of one host stack down the column; the band grows to hold them.
 *
 * Returns a GArray of nde_graph_place, one per box and in the same order; free
 * with g_array_unref().  @total_w / @total_h may be NULL.  No boxes yields an
 * empty array and a zero size, and one level yields one band — a single row of
 * nodes, which for a linear document is the one node it has.
 */
GArray *nde_graph_layout(const GArray *boxes, gint col_gap, gint row_gap,
                         gint *total_w, gint *total_h);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_GRAPH_H_ */
