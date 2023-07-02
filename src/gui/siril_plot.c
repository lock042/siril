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

// utilities
static gboolean spl_data_has_any_plot(siril_plot_data *spl_data) {
	if (!spl_data)
		return FALSE;
	return (g_list_length(spl_data->plot) + g_list_length(spl_data->plots) > 0);
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
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	if (!spl_data || !widget)
		return FALSE;

	double width =  gtk_widget_get_allocated_width(widget);
	double height = gtk_widget_get_allocated_height(widget);
	if (!siril_plot_draw(cr, spl_data, width, height))
		siril_debug_print("Problem while creating siril_plot\n");
	return FALSE;
}

gboolean create_new_siril_plot_window(gpointer p) {
	GtkWidget *window;
	GtkWidget *da;
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

	// connecting the delete-event signal, triggered when the window is closed
	// the callback frees the attached spl_data
	g_signal_connect(G_OBJECT(window), "delete-event", G_CALLBACK(on_siril_plot_window_closed), NULL);

	da = gtk_drawing_area_new();
	gtk_widget_set_size_request(da, SIRIL_PLOT_DISPLAY_WIDTH, SIRIL_PLOT_DISPLAY_HEIGHT);
	gtk_container_add(GTK_CONTAINER(window), da);

	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), NULL);
	gtk_widget_show_all(window);
	return FALSE;
}










