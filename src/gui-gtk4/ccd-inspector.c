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
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/utils.h"

/* One small surface per panel, each holding only the corner/edge/centre
 * region that panel displays — never a full-image surface.  Built in
 * set_edge_square() via siril_image_view_copy_region(), so arbitrarily
 * large images (and lazy-mode images with no contiguous view->buf) work. */
static cairo_surface_t *edge_surfaces[9] = { NULL };

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

/* GTK4 draw_func for each of the 9 inspector panels.  Each panel owns a
 * surface holding exactly the region it shows (set in set_edge_square), so
 * it just paints that surface at the origin.  The panel index 0..8 is
 * passed via user_data.  Replaces the nine GTK3 on_*_draw handlers —
 * Phase 18 removed the .ui "draw" signals but never wired up replacements,
 * leaving the inspector blank. */
static void inspector_panel_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                                    int width, int height, gpointer user_data) {
	int idx = GPOINTER_TO_INT(user_data);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);
	if (idx >= 0 && idx < 9 && edge_surfaces[idx]) {
		cairo_set_source_surface(cr, edge_surfaces[idx], 0, 0);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
		cairo_paint(cr);
	}
}

static void ccd_inspector_init_statics(void) {
	if (edge_panels[0]) return;
	for (int i = 0; i < 9; i++) {
		edge_panels[i] = GTK_WIDGET(gtk_builder_get_object(gui.builder, edge_w[i]));
		if (edge_panels[i]) {
			gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(edge_panels[i]),
			                               inspector_panel_draw_cb,
			                               GINT_TO_POINTER(i), NULL);
		}
	}
	edge_dialog_widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "edge_dialog"));
}

static void set_edge_square(gchar **panel) {
	int cvport = gfit->naxes[2] > 1 ? RGB_VPORT : RED_VPORT;

	siril_open_dialog("edge_dialog");

	const int rx = (int) gfit->rx;
	const int ry = (int) gfit->ry;
	const int widget_size = com.pref.analysis.mosaic_window / 3;
	double scale = (double) com.pref.analysis.mosaic_panel / widget_size;
	if (scale < 1.0) scale = 1.0;
	/* Number of native image pixels each panel shows; displayed shrunk to
	 * widget_size via the surface device scale. */
	const int region = (int) (widget_size * scale);
	if (region <= 0)
		return;

	ccd_inspector_init_statics();

	for (int i = 0; i < 9; i++) {
		const int col = i % 3;   /* 0=left, 1=centre, 2=right  */
		const int row = i / 3;   /* 0=top,  1=centre, 2=bottom */
		int src_x = (col == 0) ? 0 : (col == 1) ? (rx - region) / 2 : rx - region;
		int src_y = (row == 0) ? 0 : (row == 1) ? (ry - region) / 2 : ry - region;

		if (edge_surfaces[i])
			cairo_surface_destroy(edge_surfaces[i]);
		cairo_surface_t *surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24,
		                                                   region, region);
		if (cairo_surface_status(surf) != CAIRO_STATUS_SUCCESS) {
			cairo_surface_destroy(surf);
			edge_surfaces[i] = NULL;
			continue;
		}
		cairo_surface_flush(surf);
		guchar *data = cairo_image_surface_get_data(surf);
		const int stride = cairo_image_surface_get_stride(surf);
		/* Pre-clear: regions of the panel past the image edge stay black,
		 * and small images are centred by copy_region's negative-origin
		 * handling. */
		memset(data, 0, (size_t) stride * region);
		siril_image_view_copy_region(cvport, src_x, src_y, region, region,
		                             data, stride);
		cairo_surface_mark_dirty(surf);
		cairo_surface_set_device_scale(surf, scale, scale);
		edge_surfaces[i] = surf;
	}

	for (int i = 0; i < G_N_ELEMENTS(edge_w); i++)
		gtk_widget_queue_draw(edge_panels[i]);
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
