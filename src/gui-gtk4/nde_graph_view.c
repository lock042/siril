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

/* COL_GAP separates siblings across a band; ROW_GAP separates the bands, and
 * has to be deep enough for a join to be drawn in — the bar, the drop below
 * it and the arrowhead all live in that space. */
#define COL_GAP 20
#define ROW_GAP 36

/* How far above the destination the join bar sits.  It must exceed the
 * arrowhead's own length by a clear margin: at ROW_GAP/2 the drop was 4px, the
 * head was 7, and the head came out sitting across the bar itself. */
#define JOIN_DROP 18.0

/* A step summary is an ellipsized label, so its NATURAL width is the whole
 * untruncated text — under the box that used to hold these, the popover
 * clamped it, and a free layout would not.  Columns are held between these
 * instead, which also makes them uniform enough to read. */
#define NODE_MIN_W 220
#define NODE_MAX_W 420

/* How long a straight run an edge makes out of its source and into its
 * arrowhead.  A cubic's tangent at its endpoint is exact, but exactness is not
 * legibility: the edge from the first layer to a mask beside the last one is
 * 590px wide and 26px tall, so the vertical arrival the head is drawn for
 * occupies about a pixel of the sweep.  The line then reads as horizontal and
 * the head reads as a wedge stuck on sideways.  A leg gives the head something
 * to stand on, and says what the head says for long enough to be seen. */
#define ARROW_LEG 14.0

/* How far a feedback edge bows out of the corridor its forward twin uses.
 * The image feeds the mask and the mask feeds the image back, so the two run
 * between the same pair of column edges; without the bow they would be one
 * line drawn twice. */
#define FEEDBACK_BOW 22.0

typedef struct {
	gint       item_id;
	gint       level;
	gboolean   spanning;   /* a joint band: spans its participants' columns */
	gint      *span_items; /* owned; the participants a spanning band covers */
	guint      n_span_items;
	GtkWidget *child;    /* borrowed; we hold the parent ref */
} view_node;

typedef struct {
	gint     src;
	gint     dst;
	gboolean feedback;
	/* When the source is a joint band: the participant column the edge
	 * leaves under (0 = plain source).  Part of the dedup identity — one
	 * band feeds a merge once PER CHANNEL. */
	gint     src_align;
} view_edge;

/* One bracket: the host rows it opens above and closes below. */
typedef struct {
	gint64 first_rec;
	gint64 last_rec;
} view_span;

/* A mask beside the item it belongs to (see the header).  spans[0] anchors the
 * satellite's placement; a mask disabled and re-enabled mid-history carries one
 * span per contiguous masked run.  A first_rec of 0 means there is no run to
 * bracket at all — a layer mask acts on the composite, not on any stretch of
 * its layer's own steps — and then adjacency alone says whose it is. */
typedef struct {
	gint       item_id;
	gint       host;
	GArray    *spans;    /* view_span, owned */
	GtkWidget *child;    /* borrowed; we hold the parent ref */
} view_sat;

static void view_sat_clear(gpointer p) {
	view_sat *s = p;
	g_array_unref(s->spans);
}

struct _NdeGraphView {
	GtkWidget  parent_instance;
	GArray    *nodes;    /* view_node */
	GArray    *edges;    /* view_edge */
	GArray    *sats;     /* view_sat */
	GArray    *places;   /* nde_graph_place: bands, then satellites */
};

G_DEFINE_FINAL_TYPE(NdeGraphView, nde_graph_view, GTK_TYPE_WIDGET)

/* ---- geometry ----------------------------------------------------------- */

/* Measure one child the way this view wants it: as wide as its widest step
 * summary, held between the column bounds, and as tall as that width makes it. */
static void measure_child(GtkWidget *child, gint *w, gint *h) {
	int min_w, nat_w, min_h, nat_h, ignored;
	gtk_widget_measure(child, GTK_ORIENTATION_HORIZONTAL, -1,
	                   &min_w, &nat_w, &ignored, &ignored);
	*w = MAX(CLAMP(nat_w, NODE_MIN_W, NODE_MAX_W), min_w);
	gtk_widget_measure(child, GTK_ORIENTATION_VERTICAL, *w,
	                   &min_h, &nat_h, &ignored, &ignored);
	*h = MAX(nat_h, min_h);
}

/* Measure every node AND satellite and hand the sizes to the layout.  Both go
 * in together because a satellite's column is what the band's next node has to
 * start beyond — reserving that space is the layout's job, not a second pass
 * bolted on afterwards. */
static GArray *layout_now(NdeGraphView *self, gint *w, gint *h) {
	GArray *boxes = g_array_new(FALSE, TRUE, sizeof(nde_graph_box));
	for (guint i = 0; i < self->nodes->len; i++) {
		const view_node *n = &g_array_index(self->nodes, view_node, i);
		if (!gtk_widget_should_layout(n->child))
			continue;
		nde_graph_box b = { .item_id = n->item_id, .level = n->level,
		                    .spanning = n->spanning,
		                    .span_items = n->span_items,
		                    .n_span_items = n->n_span_items };
		measure_child(n->child, &b.w, &b.h);
		g_array_append_val(boxes, b);
	}
	for (guint i = 0; i < self->sats->len; i++) {
		const view_sat *s = &g_array_index(self->sats, view_sat, i);
		if (!gtk_widget_should_layout(s->child))
			continue;
		nde_graph_box b = { .item_id = s->item_id, .host_item = s->host };
		measure_child(s->child, &b.w, &b.h);
		g_array_append_val(boxes, b);
	}
	GArray *places = nde_graph_layout(boxes, COL_GAP, ROW_GAP, w, h);
	g_array_unref(boxes);
	return places;
}

static nde_graph_place *place_of_mut(GArray *places, gint item_id) {
	if (!places)
		return NULL;
	for (guint i = 0; i < places->len; i++) {
		nde_graph_place *p = &g_array_index(places, nde_graph_place, i);
		if (p->item_id == item_id)
			return p;
	}
	return NULL;
}

static const nde_graph_place *place_of(const GArray *places, gint item_id) {
	return place_of_mut((GArray *)places, item_id);
}

/* The row widget carrying "nde-record-id" == @rec, anywhere under @root.  The
 * History panel stamps its rows with it (on_hist_row_bind); the internal
 * wrappers a GtkListView puts between us and them are why this walks. */
static GtkWidget *row_for_record(GtkWidget *root, gint64 rec) {
	for (GtkWidget *c = gtk_widget_get_first_child(root); c;
	     c = gtk_widget_get_next_sibling(c)) {
		const gint64 *id = g_object_get_data(G_OBJECT(c), "nde-record-id");
		if (id && *id == rec)
			return c;
		GtkWidget *r = row_for_record(c, rec);
		if (r)
			return r;
	}
	return NULL;
}

/* Where a host row's top and bottom sit in @relative_to's coordinates.  FALSE
 * when the row does not exist (node collapsed, record gone). */
static gboolean row_bounds(GtkWidget *host_child, gint64 rec, GtkWidget *relative_to,
                           double *top, double *bottom) {
	GtkWidget *row = row_for_record(host_child, rec);
	graphene_rect_t b;
	if (!row || !gtk_widget_compute_bounds(row, relative_to, &b))
		return FALSE;
	*top = b.origin.y;
	*bottom = b.origin.y + b.size.height;
	return TRUE;
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

/* The layout stacked each host's satellites from the top of its band; slide
 * them down to meet the first host row they masked, within the range the
 * layout reserved for them.  Nudging only: the column's width and the space
 * after it are already fixed, so this cannot disturb anything to the right.
 * A satellite with no span (a layer mask, whose effect is on the composite and
 * not on any run of its layer's steps) stays where it was put. */
static void sat_anchor(NdeGraphView *self) {
	GHashTable *bottom = g_hash_table_new(g_direct_hash, g_direct_equal);
	for (guint i = 0; i < self->sats->len; i++) {
		const view_sat *s = &g_array_index(self->sats, view_sat, i);
		if (!gtk_widget_should_layout(s->child))
			continue;
		nde_graph_place *p = place_of_mut(self->places, s->item_id);
		const nde_graph_place *host = place_of(self->places, s->host);
		if (!p || !host)
			continue;
		gint y = p->y;
		const view_span *sp0 = s->spans->len
				? &g_array_index(s->spans, view_span, 0) : NULL;
		if (sp0 && sp0->first_rec) {
			/* Align the satellite's title with the first row it masked.  The
			 * bounds are relative to the HOST CHILD, whose place gives the
			 * offset into the view — the view's own subtree is what was just
			 * allocated. */
			const view_node *hn = NULL;
			for (guint k = 0; k < self->nodes->len; k++) {
				const view_node *n = &g_array_index(self->nodes, view_node, k);
				if (n->item_id == s->host) {
					hn = n;
					break;
				}
			}
			double rt, rb;
			if (hn && row_bounds(hn->child, sp0->first_rec, hn->child, &rt, &rb))
				y = host->y + (gint)rt - 4;
		}
		/* Not above its own slot in the column, not past what the layout left
		 * for the satellites still below it, and never overlapping the one
		 * before it. */
		gint floor_y = p->y_min;
		gpointer prev = g_hash_table_lookup(bottom, GINT_TO_POINTER(s->host));
		if (prev && GPOINTER_TO_INT(prev) + COL_GAP > floor_y)
			floor_y = GPOINTER_TO_INT(prev) + COL_GAP;
		p->y = CLAMP(y, floor_y, MAX(floor_y, p->y_max));
		g_hash_table_insert(bottom, GINT_TO_POINTER(s->host),
		                    GINT_TO_POINTER(p->y + p->h));
	}
	g_hash_table_destroy(bottom);
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

	/* Satellites AFTER the band nodes: their anchors are rows inside a host,
	 * and a row has a position only once its host's subtree is allocated. */
	sat_anchor(self);
	for (guint i = 0; i < self->sats->len; i++) {
		const view_sat *s = &g_array_index(self->sats, view_sat, i);
		if (!gtk_widget_should_layout(s->child))
			continue;
		const nde_graph_place *p = place_of(self->places, s->item_id);
		if (!p)
			continue;
		GskTransform *t = gsk_transform_translate(
				NULL, &GRAPHENE_POINT_INIT((float)p->x, (float)p->y));
		gtk_widget_allocate(s->child, p->w, p->h, -1, t);
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
	/* The direction the edge leaves by and the direction it arrives by: the
	 * legs run along these, and the head points along the second.  They are
	 * stated rather than read back off the control points because a control
	 * point tells you the tangent, and on a long shallow edge the tangent is
	 * not what the shape looks like. */
	double ux0, uy0, ux1, uy1;

	if (!feedback) {
		x0 = s->x + s->w / 2.0;      y0 = s->y + s->content_h;
		x1 = d->x + d->w / 2.0;      y1 = d->y;
		const double reach = MAX(fabs(y1 - y0) * 0.5, 16.0);
		c0x = x0;                    c0y = y0 + reach;
		c1x = x1;                    c1y = y1 - reach;
		ux0 = ux1 = 0.0;             uy0 = uy1 = 1.0;
	} else if (s->y == d->y) {
		/* Same band (band nodes are levelled to it, so equal y IS equal
		 * level): there is no "up the page" for the corridor geometry to run
		 * along — it used to leave the source's top and dive to the
		 * destination's band bottom, an arrow into empty space.  Run between
		 * the facing side edges at title height instead, with a slight sag so
		 * the line reads as a connector rather than a node border. */
		const gboolean rightward = d->x > s->x;
		x0 = rightward ? s->x + s->w : s->x;
		x1 = rightward ? d->x : d->x + d->w;
		y0 = y1 = s->y + 16.0;
		const double sag = 8.0;
		c0x = x0 + (x1 - x0) / 3.0;        c0y = y0 + sag;
		c1x = x0 + 2.0 * (x1 - x0) / 3.0;  c1y = y1 + sag;
		/* Both ends face along the run here, which is the one arrangement a leg
		 * cannot rescue — a horizontal stub on a horizontal line is the line.
		 * The legs are pitched at 45 instead, dipping out of the source and
		 * climbing into the destination, which is the direction the sag was
		 * already taking the middle. */
		const double diag = 0.70710678;
		ux0 = rightward ? diag : -diag;    uy0 =  diag;
		ux1 = ux0;                         uy1 = -diag;
	} else {
		x0 = s->x + s->w / 2.0 + 10.0;   y0 = s->y;
		x1 = d->x + d->w / 2.0 + 10.0;   y1 = d->y + d->content_h;
		const double reach = MAX(fabs(y1 - y0) * 0.5, 16.0);
		c0x = x0 + FEEDBACK_BOW;     c0y = y0 - reach;
		c1x = x1 + FEEDBACK_BOW;     c1y = y1 + reach;
		ux0 = ux1 = 0.0;             uy0 = uy1 = -1.0;
	}

	/* A leg never eats more than half of what the edge spans along its own
	 * axis, so the two can meet but never cross and double back. */
	const double dx = x1 - x0, dy = y1 - y0;
	const double leg0 = MIN(ARROW_LEG, fabs(dx * ux0 + dy * uy0) / 2.0);
	const double leg1 = MIN(ARROW_LEG, fabs(dx * ux1 + dy * uy1) / 2.0);

	cairo_move_to(cr, x0, y0);
	cairo_line_to(cr, x0 + leg0 * ux0, y0 + leg0 * uy0);
	cairo_curve_to(cr, c0x, c0y, c1x, c1y,
	               x1 - leg1 * ux1, y1 - leg1 * uy1);
	cairo_line_to(cr, x1, y1);
	cairo_stroke(cr);
	cairo_set_dash(cr, NULL, 0, 0);   /* the head is solid even when dashed */
	arrowhead(cr, x1, y1, x1 - ux1, y1 - uy1);
}

/* Several inputs converging on one node — a merge or a flatten — draw as a
 * BUS: a drop from each input to a shared bar, and one drop from the bar into
 * the result.  Three curves crossing each other say the same thing and say it
 * worse; the bar is also how the operation reads on paper. */
static void draw_join(cairo_t *cr, const nde_graph_place **srcs, guint n,
                      const nde_graph_place *d) {
	double bar_y = d->y - JOIN_DROP;
	double min_x = G_MAXDOUBLE, max_x = -G_MAXDOUBLE;
	for (guint i = 0; i < n; i++) {
		const double cx = srcs[i]->x + srcs[i]->w / 2.0;
		/* The END OF ITS HISTORY, not the bottom of its band: a short node in
		 * a band levelled to a taller one used to start its drop far below its
		 * own last step, with a gap of nothing in between. */
		const double bottom = srcs[i]->y + srcs[i]->content_h;
		/* An input from further up the page than the band immediately above
		 * would otherwise drop THROUGH the nodes between; keep the bar below
		 * every source it serves.  Never so close to the destination that the
		 * arrowhead has no drop to sit on, though — that is the defect this
		 * clamp exists for. */
		if (bottom + 4.0 > bar_y)
			bar_y = MIN(bottom + 4.0, d->y - 10.0);
		if (cx < min_x) min_x = cx;
		if (cx > max_x) max_x = cx;
	}
	const double drop_x = d->x + d->w / 2.0;
	if (drop_x < min_x) min_x = drop_x;
	if (drop_x > max_x) max_x = drop_x;

	for (guint i = 0; i < n; i++) {
		const double cx = srcs[i]->x + srcs[i]->w / 2.0;
		cairo_move_to(cr, cx, srcs[i]->y + srcs[i]->content_h);
		cairo_line_to(cr, cx, bar_y);
	}
	cairo_move_to(cr, min_x, bar_y);
	cairo_line_to(cr, max_x, bar_y);
	cairo_move_to(cr, drop_x, bar_y);
	cairo_line_to(cr, drop_x, d->y);
	cairo_stroke(cr);
	arrowhead(cr, drop_x, d->y, drop_x, bar_y);
}

/* The place an edge LEAVES from.  Ordinarily the source node's own place; an
 * edge out of a joint band aligned on a participant (view_edge.src_align)
 * synthesizes one from the band's y/height and the participant column's
 * x/width — clamped inside the band — so the departure sits under the
 * channel whose state the edge carries.  @storage backs the synthetic case
 * and must outlive the returned pointer's use. */
static const nde_graph_place *edge_src_place(const NdeGraphView *self,
                                             const view_edge *e,
                                             nde_graph_place *storage) {
	const nde_graph_place *s = place_of(self->places, e->src);
	if (!s || !e->src_align)
		return s;
	const nde_graph_place *a = place_of(self->places, e->src_align);
	if (!a)
		return s;
	*storage = *s;
	storage->w = MIN(a->w, s->w);
	storage->x = CLAMP(a->x, s->x, s->x + s->w - storage->w);
	return storage;
}

static void nde_graph_view_snapshot(GtkWidget *widget, GtkSnapshot *snapshot) {
	NdeGraphView *self = NDE_GRAPH_VIEW(widget);
	const int w = gtk_widget_get_width(widget);
	const int h = gtk_widget_get_height(widget);

	if ((self->edges->len || self->sats->len) && self->places && w > 0 && h > 0) {
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
			nde_graph_place srcs_store[64];
			guint n = 0;
			for (guint j = 0; j < self->edges->len && n < G_N_ELEMENTS(srcs); j++) {
				const view_edge *f = &g_array_index(self->edges, view_edge, j);
				if (f->feedback || f->dst != e->dst)
					continue;
				const nde_graph_place *s = edge_src_place(self, f,
						&srcs_store[n]);
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
			nde_graph_place fb_store;
			const nde_graph_place *s = edge_src_place(self, e, &fb_store);
			const nde_graph_place *d = place_of(self->places, e->dst);
			if (!s || !d)
				continue;
			cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.30);
			cairo_set_dash(cr, dashes, 2, 0);
			draw_edge(cr, s, d, TRUE);
		}

		/* Satellite brackets: a solid connector opening just above the first
		 * row the mask modulated, a dashed one closing just below the last.
		 * The span they enclose IS the mask's period of effect, which is what
		 * a per-row glyph used to approximate. */
		for (guint i = 0; i < self->sats->len; i++) {
			const view_sat *sat = &g_array_index(self->sats, view_sat, i);
			const nde_graph_place *sp = place_of(self->places, sat->item_id);
			const nde_graph_place *hp = place_of(self->places, sat->host);
			if (!sp || !hp)
				continue;
			const view_node *hn = NULL;
			for (guint k = 0; k < self->nodes->len; k++) {
				const view_node *n = &g_array_index(self->nodes, view_node, k);
				if (n->item_id == sat->host) {
					hn = n;
					break;
				}
			}
			const double hx = hp->x + hp->w;
			const double mid = hx + (sp->x - hx) / 2.0;

			/* No run of host steps to bracket — a layer mask, which acts where
			 * the layer is composited rather than on any of its steps.  One
			 * connector at title height says whose mask it is, and its real
			 * edge to the composite says where it is used. */
			if (!sat->spans->len ||
			    !g_array_index(sat->spans, view_span, 0).first_rec) {
				const double ty = sp->y + 16.0;
				/* Leaving the host's side, not its empty band space below. */
				const double hy = MIN(ty, hp->y + hp->content_h - 4.0);
				cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.55);
				cairo_set_dash(cr, NULL, 0, 0);
				cairo_move_to(cr, hx, hy);
				cairo_curve_to(cr, mid, hy, mid, ty, sp->x, ty);
				cairo_stroke(cr);
				arrowhead(cr, sp->x, ty, mid, ty);
				continue;
			}

			for (guint k = 0; k < sat->spans->len; k++) {
				const view_span *span = &g_array_index(sat->spans, view_span, k);
				double top = hp->y + 24.0;
				double bottom = hp->y + hp->content_h - 6.0;
				if (hn) {
					double rt, rb;
					if (row_bounds(hn->child, span->first_rec, widget, &rt, &rb))
						top = rt - 1.0;
					if (row_bounds(hn->child, span->last_rec, widget, &rt, &rb))
						bottom = rb + 1.0;
				}

				/* Opening line: host → mask, at the boundary above the first
				 * masked step, arriving at the satellite's title (first span)
				 * or at the nearest point on its edge (a later span, after
				 * the mask was re-enabled). */
				const double sy = k == 0 ? sp->y + 16.0 :
						CLAMP(top, sp->y + 16.0, sp->y + sp->h - 12.0);
				cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.55);
				cairo_set_dash(cr, NULL, 0, 0);
				cairo_move_to(cr, hx, top);
				cairo_curve_to(cr, mid, top, mid, sy, sp->x, sy);
				cairo_stroke(cr);
				arrowhead(cr, sp->x, sy, mid, sy);

				/* Closing line: mask → host, at the boundary below the last
				 * masked step.  Dashed, like every edge that runs back up the
				 * derivation. */
				const double ey = CLAMP(bottom, sy + 8.0, sp->y + sp->h - 4.0);
				cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha * 0.30);
				cairo_set_dash(cr, dashes, 2, 0);
				cairo_move_to(cr, sp->x, ey);
				cairo_curve_to(cr, mid, ey, mid, bottom, hx, bottom);
				cairo_stroke(cr);
				cairo_set_dash(cr, NULL, 0, 0);
				arrowhead(cr, hx, bottom, mid, bottom);
			}
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
	g_clear_pointer(&self->sats, g_array_unref);
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
	self->sats  = g_array_new(FALSE, TRUE, sizeof(view_sat));
	g_array_set_clear_func(self->sats, view_sat_clear);
}

/* ---- API ---------------------------------------------------------------- */

GtkWidget *nde_graph_view_new(void) {
	return g_object_new(NDE_TYPE_GRAPH_VIEW, NULL);
}

void nde_graph_view_reset(NdeGraphView *self) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	for (guint i = 0; i < self->nodes->len; i++) {
		view_node *n = &g_array_index(self->nodes, view_node, i);
		g_clear_pointer(&n->span_items, g_free);
		gtk_widget_unparent(n->child);
	}
	for (guint i = 0; i < self->sats->len; i++) {
		view_sat *s = &g_array_index(self->sats, view_sat, i);
		gtk_widget_unparent(s->child);
	}
	g_array_set_size(self->nodes, 0);
	g_array_set_size(self->edges, 0);
	g_array_set_size(self->sats, 0);
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

void nde_graph_view_add_span_node(NdeGraphView *self, gint item_id, gint level,
                                  const gint *items, guint n_items,
                                  GtkWidget *child) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	g_return_if_fail(GTK_IS_WIDGET(child));
	view_node n = { .item_id = item_id, .level = level, .spanning = TRUE,
	                .child = child };
	if (items && n_items) {
		n.span_items = g_new(gint, n_items);
		memcpy(n.span_items, items, n_items * sizeof(gint));
		n.n_span_items = n_items;
	}
	g_array_append_val(self->nodes, n);
	gtk_widget_set_parent(child, GTK_WIDGET(self));
	gtk_widget_queue_resize(GTK_WIDGET(self));
}

void nde_graph_view_add_edge(NdeGraphView *self, gint src_item, gint dst_item,
                             gboolean feedback) {
	nde_graph_view_add_edge_aligned(self, src_item, dst_item, feedback, 0);
}

void nde_graph_view_add_edge_aligned(NdeGraphView *self, gint src_item,
                                     gint dst_item, gboolean feedback,
                                     gint align_item) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	if (src_item == dst_item)
		return;
	for (guint i = 0; i < self->edges->len; i++) {
		const view_edge *e = &g_array_index(self->edges, view_edge, i);
		if (e->src == src_item && e->dst == dst_item &&
		    e->src_align == align_item)
			return;
	}
	view_edge e = { .src = src_item, .dst = dst_item, .feedback = feedback,
	                .src_align = align_item };
	g_array_append_val(self->edges, e);
	gtk_widget_queue_draw(GTK_WIDGET(self));
}

void nde_graph_view_add_satellite(NdeGraphView *self, gint item_id, gint host_item,
                                  gint64 first_use_record, gint64 last_use_record,
                                  GtkWidget *child) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	g_return_if_fail(GTK_IS_WIDGET(child));
	view_sat s = { .item_id = item_id, .host = host_item,
	               .spans = g_array_new(FALSE, TRUE, sizeof(view_span)),
	               .child = child };
	/* 0 = no run of host steps to bracket; adjacency alone says whose it is. */
	if (first_use_record) {
		view_span sp = { .first_rec = first_use_record,
		                 .last_rec = last_use_record };
		g_array_append_val(s.spans, sp);
	}
	g_array_append_val(self->sats, s);
	gtk_widget_set_parent(child, GTK_WIDGET(self));
	gtk_widget_queue_resize(GTK_WIDGET(self));
}

void nde_graph_view_satellite_add_span(NdeGraphView *self, gint item_id,
                                       gint64 first_use_record,
                                       gint64 last_use_record) {
	g_return_if_fail(NDE_IS_GRAPH_VIEW(self));
	for (guint i = 0; i < self->sats->len; i++) {
		view_sat *s = &g_array_index(self->sats, view_sat, i);
		if (s->item_id != item_id)
			continue;
		view_span sp = { .first_rec = first_use_record,
		                 .last_rec = last_use_record };
		g_array_append_val(s->spans, sp);
		gtk_widget_queue_draw(GTK_WIDGET(self));
		return;
	}
}
