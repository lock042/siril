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

#include "io/siril_plot.h"

#include <cairo.h>
#include <math.h>
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/single_image.h"

#define SIRIL_PLOT_ZOOM_OUT 1.5
#define SIRIL_PLOT_ZOOM_IN 1. / SIRIL_PLOT_ZOOM_OUT
#define SIRIL_PLOT_MAX_ZOOM 1000.
#define SIRIL_PLOT_MIN_ZOOM 0.5

// utilities
static gboolean spl_data_has_any_plot(siril_plot_data *spl_data) {
	if (!spl_data)
		return FALSE;
	// siril_debug_print("Plot: %d\n", g_list_length(spl_data->plot));
	// siril_debug_print("Plots: %d\n", g_list_length(spl_data->plots));
	return (g_list_length(spl_data->plot) + g_list_length(spl_data->plots) > 0);
}

gchar* build_save_filename(gchar *prepend, gchar *ext, gboolean forsequence, gboolean add_time_stamp){
	gchar *temp = NULL, *timestamp = NULL;
	GString *filename = NULL;
	
	if (!prepend)
		return NULL;
	filename = g_string_new(prepend);

	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		temp = g_path_get_basename(com.uniq->filename);
	} else if (sequence_is_loaded() && !forsequence) {
		char seq_image_canonical_name[256] = "";
		seq_get_image_filename(&com.seq, com.seq.current, seq_image_canonical_name);
		temp = g_strdup(seq_image_canonical_name);
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

static void convert_surface_to_plot(siril_plot_data *spl_data, double x, double y, double *xpl, double *ypl) {
	if (!spl_data->revertX)
		*xpl = spl_data->pdd.pdatamin.x + (x - spl_data->pdd.offset.x) / spl_data->pdd.range.x * (spl_data->pdd.pdatamax.x - spl_data->pdd.pdatamin.x);
	else
		*xpl = spl_data->pdd.pdatamax.x - (x - spl_data->pdd.offset.x) / spl_data->pdd.range.x * (spl_data->pdd.pdatamax.x - spl_data->pdd.pdatamin.x);
	if (!spl_data->revertY)
		*ypl = spl_data->pdd.pdatamax.y - (y - spl_data->pdd.offset.y) / spl_data->pdd.range.y * (spl_data->pdd.pdatamax.y - spl_data->pdd.pdatamin.y);
	else
		*ypl = spl_data->pdd.pdatamin.y + (y - spl_data->pdd.offset.y) / spl_data->pdd.range.y * (spl_data->pdd.pdatamax.y - spl_data->pdd.pdatamin.y);
}

static void reset_selection(plot_draw_data_t *pdd) {
	pdd->action = SELACTION_NONE;
	pdd->selection = (rectangled){0., 0., 0., 0.};
}

static void reset_zoom(siril_plot_data *spl_data) {
	spl_data->pdd.datamin = spl_data->datamin;
	spl_data->pdd.datamax = spl_data->datamax;
}

static gboolean update_zoom(siril_plot_data *spl_data, double x, double y, double scale) {
	if (!spl_data)
		return FALSE;
	double x1, y1;
	convert_surface_to_plot(spl_data, x, y, &x1, &y1);
	double xrangeorig = spl_data->datamax.x - spl_data->datamin.x;
	double xrangep = (spl_data->pdd.pdatamax.x - x1) * scale;
	double xrangem = (x1 - spl_data->pdd.pdatamin.x) * scale;
	double yrangep = (spl_data->pdd.pdatamax.y - y1) * scale;
	double yrangem = (y1 - spl_data->pdd.pdatamin.y) * scale;
	double xrangezoomratio = xrangeorig/ (xrangep + xrangem);
	if (xrangezoomratio > SIRIL_PLOT_MAX_ZOOM || xrangezoomratio < SIRIL_PLOT_MIN_ZOOM)
		return FALSE;

	spl_data->pdd.datamin.x = x1 - xrangem;
	spl_data->pdd.datamax.x = x1 + xrangep;
	spl_data->pdd.datamin.y = y1 - yrangem;
	spl_data->pdd.datamax.y = y1 + yrangep;
	return TRUE;
}

static void set_filter(GtkFileChooser *dialog, const gchar *name, const gchar *pattern) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, name);
	gtk_file_filter_add_pattern(f, pattern);
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

static gchar* save_siril_plot_dialog(GtkWindow *parent, const gchar *defaultfilename, const gchar *filter_name, const gchar *filter_pattern) {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	gint res;
	gchar *savefilename = NULL;

	widgetdialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, defaultfilename);
	gtk_file_chooser_set_local_only(dialog, FALSE);
	set_filter(dialog, filter_name, filter_pattern);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		savefilename = siril_file_chooser_get_filename(dialog);
	}
	siril_widget_destroy(widgetdialog);
	return savefilename;
}

gboolean save_siril_plot_to_clipboard(siril_plot_data *spl_data, int width, int height) {
	if (!spl_data)
		return TRUE;

	cairo_surface_t *surface = siril_plot_draw_to_image_surface(spl_data, width, height);
	if (!surface)
		return TRUE;

	GdkPixbuf *pixbuf = gdk_pixbuf_get_from_surface(surface, 0, 0, width, height);
	cairo_surface_destroy(surface);

	if (pixbuf) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
		gtk_clipboard_set_image(cb, pixbuf);
#if !defined _WIN32
		gtk_clipboard_store(cb);
#endif
		siril_log_message(_("Snapshot was saved into the clipboard.\n"));
		g_object_unref(pixbuf);
	} else {
		siril_log_message(_("Could not copy the snapshot into the clipboard.\n"));
	}
	return TRUE;

}

static GdkModifierType get_primary() {
	return gdk_keymap_get_modifier_mask(
			gdk_keymap_get_for_display(gdk_display_get_default()),
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);
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
	if (!siril_plot_draw(cr, spl_data, width, height, FALSE))
		siril_debug_print("Problem while creating siril_plot\n");
	if (spl_data->pdd.action == SELACTION_SELECTING) {
		// drawing the selection box
		if (fabs(spl_data->pdd.selection.w) > 1. || fabs(spl_data->pdd.selection.h) > 1.) {
			double xc, yc, w, h;
			if (spl_data->pdd.selection.w > 0) {
				xc = spl_data->pdd.selection.x;
				w = spl_data->pdd.selection.w;
			} else {
				w = -spl_data->pdd.selection.w;
				xc = spl_data->pdd.selection.x - w;
			}
			if (spl_data->pdd.selection.h > 0) {
				yc = spl_data->pdd.selection.y;
				h = spl_data->pdd.selection.h;
			} else {
				h = -spl_data->pdd.selection.h;
				yc = spl_data->pdd.selection.y - h;
			}
			cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
			cairo_set_line_width(cr, 1.);
			cairo_rectangle(cr, xc, yc, w, h);
			cairo_stroke(cr);
		}
	}
	return TRUE;
}

static gboolean on_siril_plot_enter_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor("crosshair");
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
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !label || !da)
		return TRUE;
	double x = (double)event->x;
	double y = (double)event->y;
	double xpos, ypos;
	
	if (is_inside_grid(x, y, &spl_data->pdd)) {
		// display cursor position in bottom left label
		convert_surface_to_plot(spl_data, x, y, &xpos, &ypos);
		gchar *labeltext = g_strdup_printf("%10g ; %10g", xpos, ypos);
		if (labeltext) {
			gtk_label_set_text(GTK_LABEL(label), labeltext);
		}
		g_free(labeltext);
		// display selection box
		if (spl_data->pdd.action == SELACTION_SELECTING && spl_data->zoomable) {
			spl_data->pdd.selection.w = x - spl_data->pdd.start.x;
			spl_data->pdd.selection.h = y - spl_data->pdd.start.y;
			gtk_widget_queue_draw(da);
			return TRUE;
		}
		if (spl_data->pdd.action == SELACTION_MOVING && spl_data->zoomable) {
			double x1, x2, y1, y2;
			convert_surface_to_plot(spl_data, spl_data->pdd.start.x, spl_data->pdd.start.y, &x1, &y1);
			convert_surface_to_plot(spl_data, x, y, &x2, &y2);
			spl_data->pdd.datamin.x += x1 - x2;
			spl_data->pdd.datamax.x += x1 - x2;
			spl_data->pdd.datamin.y += y1 - y2;
			spl_data->pdd.datamax.y += y1 - y2;
			spl_data->pdd.start = (point){x, y};
			gtk_widget_queue_draw(da);
			return TRUE;
		}
	}
	return TRUE;
}

static gboolean on_siril_plot_scroll_event(GtkWidget *widget, GdkEventScroll *event, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return TRUE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !spl_data->zoomable || !da)
		return TRUE;
	gboolean handled = FALSE;
	if (event->state & get_primary()) {
		point delta;
		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
			if (delta.y < 0) {
				handled = update_zoom(spl_data, event->x, event->y, SIRIL_PLOT_ZOOM_IN);
			}
			if (delta.y > 0) {
				handled = update_zoom(spl_data, event->x, event->y, SIRIL_PLOT_ZOOM_OUT);
			}
			break;
		case GDK_SCROLL_DOWN:
			handled = update_zoom(spl_data, event->x, event->y, SIRIL_PLOT_ZOOM_OUT);
			break;
		case GDK_SCROLL_UP:
			handled = update_zoom(spl_data, event->x, event->y, SIRIL_PLOT_ZOOM_IN);
			break;
		default:
			handled = FALSE;
		}
	}
	if (handled) {
		spl_data->autotic = TRUE;
		gtk_widget_queue_draw(da);
	}
	return TRUE;
}

static gboolean on_siril_plot_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return TRUE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *menu = (GtkWidget *)g_object_get_data(G_OBJECT(window), "menu_handle");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !menu || !da)
		return TRUE;

	// right-click pops-up menu
	if (event->button == GDK_BUTTON_SECONDARY) {
		gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
		return TRUE;
	}
	if(!spl_data->zoomable)
		return TRUE;

	// single left-click draws selection (terminates when released)
	// double left-click resets zoom
	double x = (double)event->x;
	double y = (double)event->y;
	if (event->button == GDK_BUTTON_PRIMARY) {
		if (event->type == GDK_DOUBLE_BUTTON_PRESS) {  // double-click resets zoom
			reset_zoom(spl_data);
			reset_selection(&spl_data->pdd);
			spl_data->autotic = TRUE;
			gtk_widget_queue_draw(da);
			return TRUE;
		} else if (spl_data->pdd.action == SELACTION_NONE && is_inside_grid(x, y, &spl_data->pdd)) { // start drawing selection
			if (event->state & get_primary()) {
				spl_data->pdd.action = SELACTION_MOVING; // pan start
				spl_data->autotic = FALSE;
				set_cursor("all-scroll");
			} else {
				spl_data->pdd.action = SELACTION_SELECTING; // selection
				spl_data->pdd.selection = (rectangled){x, y, 0., 0.};
			}
			spl_data->pdd.start = (point){x, y};
			return TRUE;
		}
	}
	return TRUE;
}

static gboolean on_siril_plot_button_release_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return FALSE;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !spl_data->zoomable || !da)
		return TRUE;
	if (event->button == GDK_BUTTON_PRIMARY && spl_data->pdd.action == SELACTION_SELECTING) {
		if (fabs(spl_data->pdd.selection.w) > 1. && fabs(spl_data->pdd.selection.h) > 1.) {
			double x1, x2, y1, y2;
			convert_surface_to_plot(spl_data, spl_data->pdd.selection.x, spl_data->pdd.selection.y, &x1, &y1);
			convert_surface_to_plot(spl_data, spl_data->pdd.selection.x + spl_data->pdd.selection.w, spl_data->pdd.selection.y + spl_data->pdd.selection.h, &x2, &y2);
			spl_data->pdd.datamin.x = max(min(x1, x2), spl_data->datamin.x);
			spl_data->pdd.datamax.x = min(max(x1, x2), spl_data->datamax.x);
			spl_data->pdd.datamin.y = max(min(y1, y2), spl_data->datamin.y);
			spl_data->pdd.datamax.y = min(max(y1, y2), spl_data->datamax.y);
		}
		reset_selection(&spl_data->pdd);
		gtk_widget_queue_draw(da);
	}
	if (event->button == GDK_BUTTON_PRIMARY && spl_data->pdd.action == SELACTION_MOVING) {
		reset_selection(&spl_data->pdd);
		spl_data->autotic = TRUE;
		set_cursor("crosshair");
		gtk_widget_queue_draw(da);
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
}

// generic save callback
// the name of the caller is retrieved with gtk_widget_get_name
// the appropriate save procedure is then called accordingly
static void on_siril_plot_save_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gchar *filename = NULL, *outname = NULL;
	GtkWidget *window = (GtkWidget *)(user_data);
	if (!window)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(window), "spl_data");
	GtkWidget *da = (GtkWidget *)g_object_get_data(G_OBJECT(window), "drawing_area_handle");
	if (!spl_data || !da)
		return;
	int width =  gtk_widget_get_allocated_width(da);
	int height = gtk_widget_get_allocated_height(da);
	
	const gchar *widgetname = gtk_widget_get_name(GTK_WIDGET(menuitem));
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!g_strcmp0(widgetname, "png")) {
		filename = build_save_filename(spl_data->savename, ".png", spl_data->forsequence, TRUE);
		if ((outname = save_siril_plot_dialog(GTK_WINDOW(window), filename, _("PNG files (*.png)"), "*.png")))
			siril_plot_save_png(spl_data, outname, width, height);
	} else if (!g_strcmp0(widgetname, "dat")) {
		filename = build_save_filename(spl_data->savename, ".dat", spl_data->forsequence, FALSE);
		if ((outname = save_siril_plot_dialog(GTK_WINDOW(window), filename, _("DAT files (*.dat)"), "*.dat")))
			siril_plot_save_dat(spl_data, outname, FALSE);
	} else if (!g_strcmp0(widgetname, "cb")) {
		save_siril_plot_to_clipboard(spl_data, width, height);
	}
#ifdef CAIRO_HAS_SVG_SURFACE
	else if (!g_strcmp0(widgetname, "svg")) {
		filename = build_save_filename(spl_data->savename, ".svg", spl_data->forsequence, TRUE);
		if ((outname = save_siril_plot_dialog(GTK_WINDOW(window), filename, _("SVG files (*.svg)"), "*.svg")))
			siril_plot_save_svg(spl_data, outname, width, height);
	}
#endif
	g_free(filename);
	g_free(outname);
}

gboolean create_new_siril_plot_window(gpointer p) {
	GtkWidget *window, *vbox, *da, *label, *menu;
	GtkWidget *spl_menu_grid, *spl_menu_legend;
	GtkWidget *spl_menu_save_cb, *spl_menu_save_png, *spl_menu_save_dat, *spl_menu_sep;
#ifdef CAIRO_HAS_SVG_SURFACE
	GtkWidget *spl_menu_save_svg;
#endif
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

	//prepare interactivity for spl_data
	reset_zoom(spl_data);
	reset_selection(&spl_data->pdd);
	spl_data->zoomable = spl_data->bkg == NULL;

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
	int width = (!spl_data->width) ? SIRIL_PLOT_DISPLAY_WIDTH : spl_data->width;
	int height = (!spl_data->height) ? SIRIL_PLOT_DISPLAY_WIDTH : spl_data->height;
	gtk_widget_set_size_request(da, width, height);
	gtk_box_pack_start(GTK_BOX(vbox), da, TRUE, TRUE, 0);
	gtk_widget_add_events(da, GDK_POINTER_MOTION_MASK | GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK |
	GDK_BUTTON_MOTION_MASK | GDK_BUTTON1_MOTION_MASK | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | 
	GDK_SCROLL_MASK | GDK_SMOOTH_SCROLL_MASK | GDK_TOUCH_MASK );
	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), window);
	g_signal_connect(G_OBJECT(da), "enter-notify-event", G_CALLBACK(on_siril_plot_enter_notify_event), window);
	g_signal_connect(G_OBJECT(da), "leave-notify-event", G_CALLBACK(on_siril_plot_leave_notify_event), window);
	g_signal_connect(G_OBJECT(da), "motion-notify-event", G_CALLBACK(on_siril_plot_motion_notify_event), window);
	g_signal_connect(G_OBJECT(da), "button-press-event", G_CALLBACK(on_siril_plot_button_press_event), window);
	g_signal_connect(G_OBJECT(da), "button-release-event", G_CALLBACK(on_siril_plot_button_release_event), window);
	g_signal_connect(G_OBJECT(da), "scroll-event", G_CALLBACK(on_siril_plot_scroll_event), window);	
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

	// save to clipboard
	spl_menu_save_cb = gtk_menu_item_new_with_label("Save view to Clipboard");
	gtk_widget_set_name(spl_menu_save_cb, "cb");
	g_signal_connect(G_OBJECT(spl_menu_save_cb), "activate", G_CALLBACK(on_siril_plot_save_activate), window);

	// save as png
	spl_menu_save_png = gtk_menu_item_new_with_label("Save view as PNG");
	gtk_widget_set_name(spl_menu_save_png, "png");
	g_signal_connect(G_OBJECT(spl_menu_save_png), "activate", G_CALLBACK(on_siril_plot_save_activate), window);

#ifdef CAIRO_HAS_SVG_SURFACE
	// save as svg
	spl_menu_save_svg = gtk_menu_item_new_with_label("Save view as SVG");
	gtk_widget_set_name(spl_menu_save_svg, "svg");
	g_signal_connect(G_OBJECT(spl_menu_save_svg), "activate", G_CALLBACK(on_siril_plot_save_activate), window);
#endif

	// save as dat
	spl_menu_save_dat = gtk_menu_item_new_with_label("Export dat file");
	gtk_widget_set_name(spl_menu_save_dat, "dat");
	g_signal_connect(G_OBJECT(spl_menu_save_dat), "activate", G_CALLBACK(on_siril_plot_save_activate), window);

	//fill the menu
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_grid);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_legend);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_sep);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_cb);
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_png);
#ifdef CAIRO_HAS_SVG_SURFACE
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_svg);
#endif
	gtk_container_add(GTK_CONTAINER(menu), spl_menu_save_dat);
	gtk_widget_show_all(menu);
	// and cache its handle
	g_object_set_data(G_OBJECT(window), "menu_handle", menu);

	gtk_window_present(GTK_WINDOW(window));

	gtk_widget_show_all(window);
	return FALSE;
}










