/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include "io/siril_plot.h"

#include <cairo.h>
#include "core/siril_log.h"
#include "gui/progress_and_log.h"

// utilities
static gboolean spl_data_has_any_plot(siril_plot_data *spl_data) {
	if (!spl_data)
		return FALSE;
	// siril_debug_print("Plot: %d\n", g_list_length(spl_data->plot));
	// siril_debug_print("Plots: %d\n", g_list_length(spl_data->plots));
	return (g_list_length(spl_data->plot) + g_list_length(spl_data->plots) > 0);
}

static gboolean is_inside_grid(double x, double y, plot_draw_data_t pdd) {
	if (x <= pdd.offset.x + pdd.range.x &&
		x >= pdd.offset.x &&
		y <= pdd.offset.y + pdd.range.y &&
		y >= pdd.offset.y)
			return TRUE;
	return FALSE;
}

// callbacks
static gboolean on_siril_plot_window_closed(GtkWidget *widget, GdkEvent *event, gpointer data) {
	siril_debug_print("Freeing siril_plot data and closing\n");
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(widget), "spl_data");
	free_siril_plot_data(spl_data);
	gtk_widget_destroy(widget);
	return TRUE;
}

static gboolean on_siril_plot_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
	// retrieve the parent window and its attached spl_data
	GtkWidget *window = gtk_widget_get_toplevel(widget);
	if (!window)
		return TRUE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	if (!spl_data)
		return TRUE;

	double width =  gtk_widget_get_allocated_width(widget);
	double height = gtk_widget_get_allocated_height(widget);
	if (!siril_plot_draw(cr, spl_data, width, height))
		siril_debug_print("Problem while creating siril_plot\n");
	return TRUE;
}

static gboolean on_siril_plot_enter_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor("tcross");
	return TRUE;
}

static gboolean on_siril_plot_leave_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor_waiting(FALSE);
	return TRUE;
}

static gboolean on_siril_plot_motion_notify_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data) {
	GtkWidget *window = gtk_widget_get_toplevel(widget);
	if (!window)
		return FALSE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *label = (GtkWidget *)g_object_get_data(G_OBJECT(window), "display_label_handle");
	if (!spl_data || !label)
		return TRUE;
	double x = (double)event->x;
	double y = (double)event->y;
	if (is_inside_grid(x, y, spl_data->pdd)) {
		if (!spl_data->revertX)
			x = spl_data->pdd.datamin.x + (x - spl_data->pdd.offset.x) / spl_data->pdd.range.x * (spl_data->pdd.datamax.x - spl_data->pdd.datamin.x);
		else
			x = spl_data->pdd.datamax.x - (x - spl_data->pdd.offset.x) / spl_data->pdd.range.x * (spl_data->pdd.datamax.x - spl_data->pdd.datamin.x);
		if (!spl_data->revertY)
			y = spl_data->pdd.datamax.y - (y - spl_data->pdd.offset.y) / spl_data->pdd.range.y * (spl_data->pdd.datamax.y - spl_data->pdd.datamin.y);
		else
			y = spl_data->pdd.datamin.y + (y - spl_data->pdd.offset.y) / spl_data->pdd.range.y * (spl_data->pdd.datamax.y - spl_data->pdd.datamin.y);
		gchar *labeltext = g_strdup_printf("%10g ; %10g", x, y);
		if (labeltext) {
			gtk_label_set_text(GTK_LABEL(label), labeltext);
		}
		g_free(labeltext);
	} else {
		gtk_label_set_text(GTK_LABEL(label), "-");
	}
	return TRUE;
}

gboolean create_new_siril_plot_window(gpointer p) {
	GtkWidget *window;
	GtkWidget *vbox;
	GtkWidget *da;
	GtkWidget *label;
	siril_plot_data *spl_data = (siril_plot_data *)p;

	if (!spl_data_has_any_plot(spl_data)) {
		siril_debug_print("Trying to display plot that contains no data, freeing and aborting\n");
		free_siril_plot_data(spl_data);
		spl_data = NULL;
		return FALSE;
	}

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "Siril plot");
	// attaching the spl_data to the window widget
	g_object_set_data(G_OBJECT(window), "spl_data", spl_data);

	// add css data to uniformize all backgrounds
	GtkCssProvider *cssProvider = gtk_css_provider_new();
	gchar *data = "window {color: grey; background: white;}\0";
	gtk_css_provider_load_from_data(cssProvider, data, -1, NULL);
	GtkStyleContext *styleContext = gtk_widget_get_style_context(window);
	gtk_style_context_add_provider(styleContext, GTK_STYLE_PROVIDER(cssProvider), GTK_STYLE_PROVIDER_PRIORITY_USER);

	// connect the delete-event signal, triggered when the window is closed
	// the callback frees the attached spl_data
	g_signal_connect(G_OBJECT(window), "delete-event", G_CALLBACK(on_siril_plot_window_closed), NULL);
	gtk_container_set_border_width(GTK_CONTAINER(window), 5);

	// add a vertical box
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_container_add(GTK_CONTAINER(window), vbox);

	// add the label
	label = gtk_label_new("-");
	gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
	gtk_widget_set_halign(label, GTK_ALIGN_END);
	// and cache its handle
	g_object_set_data(G_OBJECT(window), "display_label_handle", label);

	// and finally add the drawing area...
	da = gtk_drawing_area_new();
	gtk_widget_set_size_request(da, SIRIL_PLOT_DISPLAY_WIDTH, SIRIL_PLOT_DISPLAY_HEIGHT);
	gtk_box_pack_end(GTK_BOX(vbox), da, TRUE, TRUE, 0);
	gtk_widget_add_events(da, GDK_POINTER_MOTION_MASK | GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK);
	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), NULL);
	g_signal_connect(G_OBJECT(da), "enter-notify-event", G_CALLBACK(on_siril_plot_enter_notify_event), NULL);
	g_signal_connect(G_OBJECT(da), "leave-notify-event", G_CALLBACK(on_siril_plot_leave_notify_event), NULL);
	g_signal_connect(G_OBJECT(da), "motion-notify-event", G_CALLBACK(on_siril_plot_motion_notify_event), NULL);
	gtk_widget_show_all(window);
	return FALSE;
}










