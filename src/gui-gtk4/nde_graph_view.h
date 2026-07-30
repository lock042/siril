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
#ifndef _NDE_GRAPH_VIEW_H_
#define _NDE_GRAPH_VIEW_H_

/**
 * \file nde_graph_view.h
 * \brief The History window's node container: places its children by graph
 *        level and paints the edges between them.
 *
 * It draws the EDGES and nothing else.  The nodes are ordinary widgets — the
 * frame, expander and GtkListView the History window already built — so the
 * row factory, drag-to-reorder, the popover, the focus chain and the whole
 * accessible tree come along unchanged.  A node graph that drew its own nodes
 * would have had to reimplement every one of those.
 *
 * Layout itself is not here: nde_graph_layout() is a pure function of the
 * measured sizes, tested without GTK.  This widget measures, calls it, and
 * allocates what comes back.
 */

#include <gtk/gtk.h>

G_BEGIN_DECLS

#define NDE_TYPE_GRAPH_VIEW (nde_graph_view_get_type())
G_DECLARE_FINAL_TYPE(NdeGraphView, nde_graph_view, NDE, GRAPH_VIEW, GtkWidget)

GtkWidget *nde_graph_view_new(void);

/** Drop every node and edge, unparenting the node widgets.  Popovers parented
 *  here for positioning are left alone: they are not nodes. */
void nde_graph_view_reset(NdeGraphView *self);

/** Add @child as the node for @item_id at @level.  Takes the parent ref. */
void nde_graph_view_add_node(NdeGraphView *self, gint item_id, gint level,
                             GtkWidget *child);

/** Draw an edge from @src_item to @dst_item.  Duplicates are dropped: several
 *  records may consume the same item, and one line says it once.  @feedback
 *  edges are the ones that run against the levels (nde_graph.h). */
void nde_graph_view_add_edge(NdeGraphView *self, gint src_item, gint dst_item,
                             gboolean feedback);

/**
 * Add @child as a SATELLITE node: a mask belonging to @host_item.  It is laid
 * out in a column IMMEDIATELY BESIDE the host, inside the host's band, with
 * the band's next node starting beyond that column — a mask is an aside to one
 * item's history, not a stage of the document's, and putting it after every
 * band node left its connectors crossing the layers in between.
 *
 * When @first_use_record is non-zero the mask masked a run of the host's own
 * steps, and the relationship is drawn as a bracket around that run rather
 * than as node-to-node edges: a solid connector leaving the host just above
 * the row of @first_use_record (the step its pixels first modulated) and a
 * dashed one returning just below the row of @last_use_record.  Rows are
 * found by the "nde-record-id" data the History panel sets on each row
 * widget; when a row is missing (collapsed node) the host's frame edges
 * stand in.
 *
 * Pass 0 for both records when there is no such run — a LAYER mask applies
 * where its layer is composited, not to any stretch of the layer's own steps.
 * A single connector then says whose mask it is, and the mask's ordinary edge
 * to the composite says where it is used.  Takes the parent ref, like add_node.
 */
void nde_graph_view_add_satellite(NdeGraphView *self, gint item_id, gint host_item,
                                  gint64 first_use_record, gint64 last_use_record,
                                  GtkWidget *child);

/**
 * Add a further bracket to the satellite @item_id (added above).  A mask can
 * be DISABLED for a step and re-enabled after it (missedasinh.png): one
 * bracket around the whole range would claim the mask covered the gap, so
 * each contiguous run of masked steps gets its own.  The satellite itself is
 * placed against the FIRST span; later spans only add connectors.
 */
void nde_graph_view_satellite_add_span(NdeGraphView *self, gint item_id,
                                       gint64 first_use_record,
                                       gint64 last_use_record);

G_END_DECLS

#endif /* _NDE_GRAPH_VIEW_H_ */
