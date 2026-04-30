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

/* GTK draw callbacks and aberration-inspector UI for the CCD inspector.
 * Processing logic lives in algos/ccd-inspector.c. */

#include <gtk/gtk.h>
#include <cairo/cairo.h>

#include "core/siril.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/dialogs.h"
#include "gui/utils.h"

static cairo_surface_t *edge_surface = NULL;
static int image_width = -1;
static int image_height = -1;

static char *edge_w[] = {
		"left_top",
		"center_top",
		"right_top",
		"left_center",
		"center_center",
		"right_center",
		"left_bottom",
		"center_bottom",
		"right_bottom"
};

static GtkWidget *edge_panels[9] = { NULL };
static GtkWidget *edge_dialog_widget = NULL;

static void ccd_inspector_init_statics(void) {
	if (edge_panels[0]) return;
	for (int i = 0; i < 9; i++)
		edge_panels[i] = GTK_WIDGET(gtk_builder_get_object(gui.builder, edge_w[i]));
	edge_dialog_widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "edge_dialog"));
}

static void set_edge_square(gchar **panel) {
	int cvport = gfit->naxes[2] > 1 ? RGB_VPORT : RED_VPORT;

	struct image_view *view = &gui.view[cvport];

	siril_open_dialog("edge_dialog");
	if (edge_surface)
		cairo_surface_destroy(edge_surface);

	image_width = gfit->rx;
	image_height = gfit->ry;
	/* New surface as we modify it */
	edge_surface = cairo_image_surface_create_for_data(view->buf, CAIRO_FORMAT_RGB24, image_width, image_height, view->full_surface_stride);

	if (cairo_surface_status(edge_surface) != CAIRO_STATUS_SUCCESS) {
		cairo_surface_destroy(edge_surface);
		edge_surface = NULL;
		return;
	}
	int widget_size = com.pref.analysis.mosaic_window / 3;
	double scale = (double) com.pref.analysis.mosaic_panel / widget_size;
	if (scale < 1.0) scale = 1.0;
	cairo_surface_set_device_scale(edge_surface, scale, scale);
	image_width = (int) ((double)image_width / scale);
	image_height = (int) ((double) image_height / scale);

	ccd_inspector_init_statics();
	for (int i = 0; i < G_N_ELEMENTS(edge_w); i++)
		gtk_widget_queue_draw(edge_panels[i]);
}

gboolean on_left_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_left_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_left_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width(widget);
	int area_height = gtk_widget_get_allocated_height(widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

void compute_aberration_inspector(void) {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		ccd_inspector_init_statics();
		int widget_size = com.pref.analysis.mosaic_window / 3;
		for (int i = 0; i < 9; i++)
			gtk_widget_set_size_request(edge_panels[i], widget_size, widget_size);
		set_edge_square(edge_w);
	}
}

void redraw_aberration_inspector(void) {
	ccd_inspector_init_statics();
	if (!gtk_widget_is_visible(edge_dialog_widget)) return;
	compute_aberration_inspector();
}
