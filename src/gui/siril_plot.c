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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gui/utils.h"
#include "core/siril_log.h"

static gboolean on_siril_plot_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {

	siril_plot_data *spl_data = (siril_plot_data *)user_data;
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

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "Siril plot");
	// // g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK(clear_siril_plot_data), spl_data);

	da = gtk_drawing_area_new();
	gtk_widget_set_size_request(da, SIRIL_PLOT_DISPLAY_WIDTH, SIRIL_PLOT_DISPLAY_HEIGHT);
	gtk_container_add(GTK_CONTAINER(window), da);

	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), spl_data);
	gtk_widget_show_all(window);
	return FALSE;
}










