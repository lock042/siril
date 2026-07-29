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

#include <math.h>

#include "core/nde_graph.h"
#include "gui-gtk4/nde_graph_view.h"

/* Gaps between columns and between the nodes stacked within one. */
#define COL_GAP 28
#define ROW_GAP 8

/* A step summary is an ellipsized label, so its NATURAL width is the whole
 * untruncated text — under the box that used to hold these, the popover
 * clamped it, and a free layout would not.  Columns are held between these
 * instead, which also makes them uniform enough to read. */
#define NODE_MIN_W 220
#define NODE_MAX_W 420

/* How far a feedback edge bows out of the corridor its forward twin uses.
 * The image feeds the mask and the mask feeds the image back, so the two run
 * between the same pair of column edges; without the bow they would be one
 * line drawn twice. */
#define FEEDBACK_BOW 22.0

typedef struct {
	gint       item_id;
	gint       level;
	GtkWidget *child;    /* borrowed; we hold the parent ref */
} view_node;

typedef struct {
	gint     src;
	gint     dst;
	gboolean feedback;
} view_edge;

struct _NdeGraphView {
	GtkWidget  parent_instance;
	GArray    *nodes;    /* view_node */
	GArray    *edges;    /* view_edge */
	GArray    *places;   /* nde_graph_place, from the last allocation */
};

G_DEFINE_FINAL_TYPE(NdeGraphView, nde_graph_view, GTK_TYPE_WIDGET)

/* ---- geometry ----------------------------------------------------------- */

/* Measure every node and hand the sizes to the layout.  Children are measured
 * at their natural size: a node is as wide as its widest step summary and as
 * tall as its step count, and nothing here wants to negotiate that. */
static GArray *layout_now(NdeGraphView *self, gint *w, gint *h) {
	GArray *boxes = g_array_new(FALSE, TRUE, sizeof(nde_graph_box));
	for (guint i = 0; i < self->nodes->len; i++) {
		const view_node *n = &g_array_index(self->nodes, view_node, i);
		if (!gtk_widget_should_layout(n->child))
			continue;
		int min_w, nat_w, min_h, nat_h, ignored;
		gtk_widget_measure(n->child, GTK_ORIENTATION_HORIZONTAL, -1,
		                   &min_w, &nat_w, &ignored, &ignored);
		const int w = MAX(CLAMP(nat_w, NODE_MIN_W, NODE_MAX_W), min_w);
		gtk_widget_measure(n->child, GTK_ORIENTATION_VERTICAL, w,
		                   &min_h, &nat_h, &ignored, &ignored);
		nde_graph_box b = {
			.item_id = n->item_id,
			.level   = n->level,
			.w       = w,
			.h       = MAX(nat_h, min_h),
		};
		g_array_append_val(boxes, b);
	}
	GArray *places = nde_graph_layout(boxes, COL_GAP, ROW_GAP, w, h);
	g_array_unref(boxes);
	return places;
}

static const nde_graph_place *place_of(const GArray *places, gint item_id) {
	if (!places)
		return NULL;
	for (guint i = 0; i < places->len; i++) {
		const nde_graph_place *p = &g_array_index(places, nde_graph_place, i);
		if (p->item_id == item_id)
			return p;
	}
	return NULL;
}

static void nde_graph_view_measure(GtkWidget *widget, GtkOrientation orientation,
                                   int for_size, int *minimum, int *natural,
                                   int *minimum_baseline, int *natural_baseline) {
	(void)for_size;
	NdeGraphView *self = NDE_GRAPH_VIEW(widget);
	gint w = 0, h = 0;
	GArray *places = layout_now(self, &w, &h);
	g_array_unref(places);

	*minimum = *natural = (orientation == GTK_ORIENTATION_HORIZONTAL) ? w : h;
	*minimum_baseline = *natural_baseline = -1;
}

static void nde_graph_view_size_allocate(GtkWidget *widget, int width, int height,
                                         int baseline) {
	(void)width; (void)height; (void)baseline;
	NdeGraphView *self = NDE_GRAPH_VIEW(widget);

	g_clear_pointer(&self->places, g_array_unref);
	self->places = layout_now(self, NULL, NULL);

	for (guint i = 0; i < self->nodes->len; i++) {
		const view_node *n = &g_array_index(self->nodes, view_node, i);
		if (!gtk_widget_should_layout(n->child))
			continue;
		const nde_graph_place *p = place_of(self->places, n->item_id);
		if (!p)
			continue;
		GskTransform *t = gsk_transform_translate(
				NULL, &GRAPHENE_POINT_INIT((float)p->x, (float)p->y));
		gtk_widget_allocate(n->child, p->w, p->h, -1, t);
	}

	/* A popover parented here for positioning is not a node and is not laid
	 * out, but GTK still requires it be presented from size_allocate. */
	for (GtkWidget *c = gtk_widget_get_first_child(widget); c;
	     c = gtk_widget_get_next_sibling(c)) {
		if (GTK_IS_POPOVER(c))
			gtk_popover_present(GTK_POPOVER(c));
	}
}

/* ---- edges -------------------------------------------------------------- */

static void arrowhead(cairo_t *cr, double x, double y, double from_x, double from_y) {
	const double dx = x - from_x, dy = y - from_y;
	const double len = hypot(dx, dy);
	if (len < 1e-3)
		return;
	const double ux = dx / len, uy = dy / len;
	const double size = 7.0, half = 3.5;
	cairo_move_to(cr, x, y);
	cairo_line_to(cr, x - size * ux + half * uy, y - size * uy - half * ux);
	cairo_line_to(cr, x - size * ux - half * uy, y - size * uy + half * ux);
	cairo_close_path(cr);
	cairo_fill(cr);
}

/* Forward edges leave the source's bottom edge and enter the destination's
 * top: the bands run down the page, so that is the direction derivation runs.
 * Feedback edges take the same corridor backwards, bowed clear of it. */
static void draw_edge(cairo_t *cr, const nde_graph_place *s,
                      const nde_graph_place *d, gboolean feedback) {
	double x0, y0, x1, y1, c0x, c0y, c1x, c1y;

	if (!feedback) {
		x0 = s->x + s->w / 2.0;      y0 = s->y + s->h;
		x1 = d->x + d->w / 2.0;      y1 = d->y;
		const double reach = MAX(fabs(y1 - y0) * 0.5, 16.0);
		c0x = x0;                    c0y = y0 + reach;
		c1x = x1;                    c1y = y1 - reach;
	} else {
		x0 = s->x + s->w / 2.0 + 10.0;   y0 = s->y;
		x1 = d->x + d->w / 2.0 + 10.0;   y1 = d->y + d->h;
		const double reach = MAX(fabs(y1 - y0) * 0.5, 16.0);
		c0x = x0 + FEEDBACK_BOW;     c0y = y0 - reach;
		c1x = x1 + FEEDBACK_BOW;     c1y = y1 + reach;
	}

	cairo_move_to(cr, x0, y0);
	cairo_curve_to(cr, c0x, c0y, c1x, c1y, x1, y1);
	cairo_stroke(cr);
	cairo_set_dash(cr, NULL, 0, 0);   /* the head is solid even when dashed */
	arrowhead(cr, x1, y1, c1x, c1y);
}

/* Several inputs converging on one node — a merge or a flatten — draw as a
 * BUS: a drop from each input to a shared bar, and one drop from the bar into
 * the result.  Three curves crossing each other say the same thing and say it
 * worse; the bar is also how the operation reads on paper. */
static void draw_join(cairo_t *cr, const nde_graph_place **srcs, guint n,
                      const nde_graph_place *d) {
	double bar_y = d->y - ROW_GAP / 2.0;
	double min_x = G_MAXDOUBLE, max_x = -G_MAXDOUBLE;
	for (guint i = 0; i < n; i++) {
		const double cx = srcs[i]->x + srcs[i]->w / 2.0;
		const double bottom = srcs[i]->y + srcs[i]->h;
		/* An input from further up the page than the band immediately above
		 * would otherwise drop THROUGH the nodes between; keep the bar below
		 * every source it serves. */
		if (bottom + 4.0 > bar_y)
			bar_y = bottom + 4.0;
		if (cx < min_x) min_x = cx;
		if (cx > max_x) max_x = cx;
	}
	const double drop_x = d->x + d->w / 2.0;
	if (drop_x < min_x) min_x = drop_x;
	if (drop_x > max_x) max_x = drop_x;

	for (guint i = 0; i < n; i++) {
		const double cx = srcs[i]->x + srcs[i]->w / 2.0;
		cairo_move_to(cr, cx, srcs[i]->y + srcs[i]->h);
		cairo_line_to(cr, cx, bar_y);
	}
	cairo_move_to(cr, min_x, bar_y);
	cairo_line_to(cr, max_x, bar_y);
	cairo_move_to(cr, drop_x, bar_y);
	cairo_line_to(cr, drop_x, d->y);
	cairo_stroke(cr);
	arrowhead(cr, drop_x, d->y, drop_x, bar_y);
}

static void nde_graph_view_snapshot(GtkWidget *widget, GtkSnapshot *snapshot) {
	NdeGraphView *self = NDE_GRAPH_VIEW(widget);
	const int w = gtk_widget_get_width(widget);
	const int h = gtk_widget_get_height(widget);

	if (self->edges->len && self->places && w > 0 && h > 0) {
		GdkRGBA fg;
		gtk_widget_get_color(widget, &fg);
		cairo_t *cr = gtk_snapshot_append_cairo(
				snapshot, &GRAPHENE_RECT_INIT(0, 0, (float)w, (float)h));
		cairo_set_line_width(cr, 1.5);
		cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
		const double dashes[] = { 4.0, 3.0 };
		/* Forward edges are grouped by destination first, so that the several
		 * inputs of one merge become a single join rather than a sheaf. */
		GHashTable *drawn = g_hash_table_new(g_direct_hash, g_direct_equal);
		for (guint i = 0; i < self->edges->len; i++) {
			const view_edge *e = &g_array_index(self->edges, view_edge, i);
			if (e->feedback ||
			    g_hash_table_contains(drawn, GINT_TO_POINTER(e->dst)))
				continue;
			g_hash_table_add(drawn, GINT_TO_POINTER(e->dst));
			const nde_graph_place *d = place_of(self->places, e->dst);
			if (!d)
				continue;
			const nde_graph_place *srcs[64];
			guint n = 0;
			for (guint j = 0; j < self->edges->len && n < G_N_ELEMENTS(srcs); j++) {
				const view_edge *f = &g_array_index(self->edges, view_edge, j);
				if (f->feedback || f->dst != e->dst)
					continue;
				const nde_graph_place *s = place_of(self->places, f->src);
				if (s)
					srcs[n++] = s;
			}
			if (!n)
				continue;
			cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.55);
			cairo_set_dash(cr, NULL, 0, 0);
			if (n == 1)
				draw_edge(cr, srcs[0], d, FALSE);
			else
				draw_join(cr, srcs, n, d);
		}
		g_hash_table_destroy(drawn);
		for (guint i = 0; i < self->edges->len; i++) {
			const view_edge *e = &g_array_index(self->edges, view_edge, i);
			if (!e->feedback)
				continue;
			const nde_graph_place *s = place_of(self->places, e->src);
			const nde_graph_place *d = place_of(self->places, e->dst);
			if (!s || !d)
				continue;
			cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.30);
			cairo_set_dash(cr, dashes, 2, 0);
			draw_edge(cr, s, d, TRUE);
		}
		cairo_destroy(cr);
	}

	/* Chain up last: GTK draws the children, so the nodes cover the edges
	 * rather than the other way round. */
	GTK_WIDGET_CLASS(nde_graph_view_parent_class)->snapshot(widget, snapshot);
}

/* ---- object ------------------------------------------------------------- */

static void nde_graph_view_dispose(GObject *object) {
	NdeGraphView *self = NDE_GRAPH_VIEW(object);
	GtkWidget *child;
	while ((child = gtk_widget_get_first_child(GTK_WIDGET(self))))
		gtk_widget_unparent(child);
	g_clear_pointer(&self->nodes, g_array_unref);
	g_clear_pointer(&self->edges, g_array_unref);
	g_clear_pointer(&self->places, g_array_unref);
	G_OBJECT_CLASS(nde_graph_view_parent_class)->dispose(object);
}

static void nde_graph_view_class_init(NdeGraphViewClass *klass) {
	GObjectClass   *object_class = G_OBJECT_CLASS(klass);
	GtkWidgetClass *widget_class = GTK_WIDGET_CLASS(klass);

	object_class->dispose        = nde_graph_view_dispose;
	widget_class->measure        = nde_graph_view_measure;
	widget_class->size_allocate  = nde_graph_view_size_allocate;
	widget_class->snapshot       = nde_graph_view_snapshot;
	gtk_widget_class_set_css_name(widget_class, "ndegraphview");
	gtk_widget_class_set_accessible_role(widget_class, GTK_ACCESSIBLE_ROLE_GROUP);
}

static void nde_graph_view_init(NdeGraphView *self) {
	self->nodes = g_array_new(FALSE, TRUE, sizeof(view_node));
	self->edges = g_array_new(FALSE, TRUE, sizeof(view_edge));
}

/* ---- API ---------------------------------------------------------------- */

GtkWidget *nde_graph_view_new(void) {
	return g_object_new(NDE_TYPE_GRAPH_VIEW, NULL);
}

void nde_graph_view_reset(NdeGraphView *self) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	for (guint i = 0; i < self->nodes->len; i++) {
		view_node *n = &g_array_index(self->nodes, view_node, i);
		gtk_widget_unparent(n->child);
	}
	g_array_set_size(self->nodes, 0);
	g_array_set_size(self->edges, 0);
	g_clear_pointer(&self->places, g_array_unref);
	gtk_widget_queue_resize(GTK_WIDGET(self));
}

void nde_graph_view_add_node(NdeGraphView *self, gint item_id, gint level,
                             GtkWidget *child) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	g_return_if_fail(GTK_IS_WIDGET(child));
	view_node n = { .item_id = item_id, .level = level, .child = child };
	g_array_append_val(self->nodes, n);
	gtk_widget_set_parent(child, GTK_WIDGET(self));
	gtk_widget_queue_resize(GTK_WIDGET(self));
}

void nde_graph_view_add_edge(NdeGraphView *self, gint src_item, gint dst_item,
                             gboolean feedback) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	if (src_item == dst_item)
		return;
	for (guint i = 0; i < self->edges->len; i++) {
		const view_edge *e = &g_array_index(self->edges, view_edge, i);
		if (e->src == src_item && e->dst == dst_item)
			return;
	}
	view_edge e = { .src = src_item, .dst = dst_item, .feedback = feedback };
	g_array_append_val(self->edges, e);
	gtk_widget_queue_draw(GTK_WIDGET(self));
}
