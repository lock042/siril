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
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/single_image.h"

// utilities
static gboolean spl_data_has_any_plot(siril_plot_data *spl_data) {
	if (!spl_data)
		return FALSE;
	// siril_debug_print("Plot: %d\n", g_list_length(spl_data->plot));
	// siril_debug_print("Plots: %d\n", g_list_length(spl_data->plots));
	return (g_list_length(spl_data->plot) + g_list_length(spl_data->plots) > 0);
}

static gboolean is_inside_grid(double x, double y, plot_draw_data_t *pdd) {
	if (x <= pdd->offset.x + pdd->range.x &&
		x >= pdd->offset.x &&
		y <= pdd->offset.y + pdd->range.y &&
		y >= pdd->offset.y)
			return TRUE;
	return FALSE;
}

static gchar* build_save_filename(gchar *prepend, gchar *ext, gboolean forsequence, gboolean add_time_stamp){
	gchar *temp = NULL, *timestamp = NULL;
	GString *filename = NULL;
	
	filename = g_string_new(prepend);

	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		temp = g_path_get_basename(com.uniq->filename);
	} else if (sequence_is_loaded() && !forsequence) {
		char* seq_image_canonical_name = calloc(256, 1);
		seq_get_image_filename(&com.seq, com.seq.current, seq_image_canonical_name);
		temp = g_strdup(seq_image_canonical_name);
		free(seq_image_canonical_name);
	}
	if (temp) {
		gchar *tmp = remove_ext_from_filename(temp);
		g_string_append_printf(filename, "_%s", tmp);
		g_free(temp);
		g_free(tmp);
	}

		timestamp = build_timestamp_filename();
	if (add_time_stamp) {
		g_string_append_printf(filename, "_%s", timestamp);
		g_free(timestamp);
	}
	g_string_append_printf(filename, "%s", ext);
	return g_string_free(filename, FALSE);
}

// callbacks
static gboolean on_siril_plot_window_closed(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	siril_debug_print("Freeing siril_plot data and closing\n");
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(widget), "spl_data");
	free_siril_plot_data(spl_data);
	gtk_widget_destroy(widget);
	return TRUE;
}

static gboolean on_siril_plot_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
	// retrieve the parent window and its attached spl_data
	GtkWidget *window = (GtkWidget *)(user_data);
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
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return FALSE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *label = (GtkWidget *)g_object_get_data(G_OBJECT(window), "display_label_handle");
	if (!spl_data || !label)
		return TRUE;
	double x = (double)event->x;
	double y = (double)event->y;
	if (is_inside_grid(x, y, &spl_data->pdd)) {
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
	}
	return TRUE;
}

static gboolean on_siril_plot_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return FALSE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *menu = (GtkWidget *)g_object_get_data(G_OBJECT(window), "menu_handle");
	if (!spl_data || !menu)
		return TRUE;
	if (event->button == GDK_BUTTON_SECONDARY) {
		gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
	}
	return TRUE;
}

static void on_siril_plot_grid_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !da)
		return;
	spl_data->cfgplot.grid = (gtk_check_menu_item_get_active(checkmenuitem)) ? GRID_ALL : 0;
	gtk_widget_queue_draw(da);
}

static void on_siril_plot_legend_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !da)
		return;
	spl_data->show_legend = gtk_check_menu_item_get_active(checkmenuitem);
	gtk_widget_queue_draw(da);
	return;
}

static void on_siril_plot_save_png_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	if (!spl_data)
		return;
	gchar *filename = build_save_filename(spl_data->savename, ".png", spl_data->forsequence, TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_plot_save_png(spl_data, filename);
	g_free(filename);
}

static void on_siril_plot_save_dat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	if (!spl_data)
		return;
	gchar *filename = build_save_filename(spl_data->savename, ".dat", spl_data->forsequence, TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_plot_save_dat(spl_data, filename);
	g_free(filename);
}

gboolean create_new_siril_plot_window(gpointer p) {
	GtkWidget *window, *vbox, *da, *label, *menu;
	GtkWidget *spl_menu_grid, *spl_menu_legend;
	GtkWidget *spl_menu_save_png, *spl_menu_save_dat, *spl_menu_sep;
	siril_plot_data *spl_data = (siril_plot_data *)p;

	// sanity checks
	if (!spl_data) {
		siril_debug_print("Passed an empty spl_data structure\n");
		return FALSE;
	}
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
	gchar *data = "window {color: grey; background: white; font-size: 12px}\0";
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

	// add the drawing area
	da = gtk_drawing_area_new();
	gtk_widget_set_size_request(da, SIRIL_PLOT_DISPLAY_WIDTH, SIRIL_PLOT_DISPLAY_HEIGHT);
	gtk_box_pack_start(GTK_BOX(vbox), da, TRUE, TRUE, 0);
	gtk_widget_add_events(da, GDK_POINTER_MOTION_MASK | GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK | GDK_BUTTON_MOTION_MASK | GDK_BUTTON1_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK);
	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), window);
	g_signal_connect(G_OBJECT(da), "enter-notify-event", G_CALLBACK(on_siril_plot_enter_notify_event), window);
	g_signal_connect(G_OBJECT(da), "leave-notify-event", G_CALLBACK(on_siril_plot_leave_notify_event), window);
	g_signal_connect(G_OBJECT(da), "motion-notify-event", G_CALLBACK(on_siril_plot_motion_notify_event), window);
	g_signal_connect(G_OBJECT(da), "button-press-event", G_CALLBACK(on_siril_plot_button_press_event), window);
	// and cache its handle
	g_object_set_data(G_OBJECT(window), "drawing_area_handle", da);

	// add the label
	label = gtk_label_new("0;0");
	gtk_box_pack_end(GTK_BOX(vbox), label, FALSE, FALSE, 0);
	gtk_widget_set_halign(label, GTK_ALIGN_START);
	// and cache its handle
	g_object_set_data(G_OBJECT(window), "display_label_handle", label);

	// add the pop-up menu
	menu = gtk_menu_new();

	// grid
	gtk_menu_attach_to_widget(GTK_MENU(menu), window, NULL);
	spl_menu_grid = gtk_check_menu_item_new_with_label("Grid");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(spl_menu_grid), TRUE);
	g_signal_connect(G_OBJECT(spl_menu_grid), "toggled", G_CALLBACK(on_siril_plot_grid_toggled), window);

	// legend
	spl_menu_legend = gtk_check_menu_item_new_with_label("Legend");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(spl_menu_legend), TRUE);
	g_signal_connect(G_OBJECT(spl_menu_legend), "toggled", G_CALLBACK(on_siril_plot_legend_toggled), window);

	// sep
	spl_menu_sep = gtk_separator_menu_item_new();

	// save as png
	spl_menu_save_png = gtk_menu_item_new_with_label("Save view as PNG");
	g_signal_connect(G_OBJECT(spl_menu_save_png), "activate", G_CALLBACK(on_siril_plot_save_png_activate), window);

	// save as dat
	spl_menu_save_dat = gtk_menu_item_new_with_label("Export dat file");
	g_signal_connect(G_OBJECT(spl_menu_save_dat), "activate", G_CALLBACK(on_siril_plot_save_dat_activate), window);

	//fill the menu
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_grid);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_legend);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_sep);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_png);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_dat);

	gtk_widget_show_all(menu);
	// and cache its handle
	g_object_set_data(G_OBJECT(window), "menu_handle", menu);

	gtk_widget_show_all(window);
	return FALSE;
}










