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
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
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
	// siril_log_debug("Plot: %d\n", g_list_length(spl_data->plot));
	// siril_log_debug("Plots: %d\n", g_list_length(spl_data->plots));
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

static void set_filter(SirilFileChooser *fc, const gchar *name, const gchar *pattern) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, name);
	gtk_file_filter_add_pattern(f, pattern);
	siril_fc_add_filter(fc, f, TRUE);
	g_object_unref(f);
}

static gchar* save_siril_plot_dialog(GtkWindow *parent, const gchar *defaultfilename, const gchar *filter_name, const gchar *filter_pattern) {
	SirilFileChooser *fc = siril_fc_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_set_current_name(fc, defaultfilename);
	set_filter(fc, filter_name, filter_pattern);

	gchar *savefilename = NULL;
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT)
		savefilename = siril_fc_get_filename(fc);
	siril_fc_destroy(fc);
	return savefilename;
}

gboolean save_siril_plot_to_clipboard(siril_plot_data *spl_data, int width, int height) {
	if (!spl_data)
		return TRUE;

	cairo_surface_t *surface = siril_plot_draw_to_image_surface(spl_data, width, height);
	if (!surface)
		return TRUE;

	GdkTexture *tex = siril_texture_from_cairo_surface(surface);
	cairo_surface_destroy(surface);
	if (tex) {
		GdkClipboard *cb = gdk_display_get_clipboard(gdk_display_get_default());
		gdk_clipboard_set_texture(cb, tex);
		g_object_unref(tex);
		siril_log_message(_("Snapshot was saved into the clipboard.\n"));
	} else {
		siril_log_message(_("Could not copy the snapshot into the clipboard.\n"));
	}
	return TRUE;

}

static GdkModifierType get_primary() {
	/* GTK4: GdkKeymap is gone.  The "primary" accelerator modifier is
	 * Cmd (META) on macOS and Ctrl elsewhere. */
#ifdef OS_OSX
	return GDK_META_MASK;
#else
	return GDK_CONTROL_MASK;
#endif
}

// callbacks

// closes a siril-plot window: frees every spl_data displayed in it (one for a
// single-plot window, several for a grouped one) and unparents every
// right-click popover (they are attached via gtk_widget_set_parent rather
// than through a container; GTK4 will warn ("still has children left") if
// still parented when the window is finalized).
static gboolean on_siril_plot_window_closed(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	siril_log_debug("Freeing siril_plot data and closing\n");
	GList *spl_data_list = (GList *)g_object_get_data(G_OBJECT(widget), "spl_data_list");
	g_list_free_full(spl_data_list, (GDestroyNotify)free_siril_plot_data);
	GList *menu_list = (GList *)g_object_get_data(G_OBJECT(widget), "menu_list");
	for (GList *l = menu_list; l; l = l->next)
		gtk_widget_unparent(GTK_WIDGET(l->data));
	g_list_free(menu_list);
	g_list_free((GList *)g_object_get_data(G_OBJECT(widget), "da_list")); // widgets go with the window
	gtk_window_destroy(GTK_WINDOW(widget));
	return TRUE;
}

static void on_siril_plot_draw(GtkDrawingArea *area, cairo_t *cr,
                                int width_i, int height_i, gpointer user_data) {
	// each pane's drawing area carries its own spl_data
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(area), "spl_data");
	if (!spl_data)
		return;

	double width = (double)width_i;
	double height = (double)height_i;
	// the header (caption, title, metrics) is drawn as widgets around the
	// canvas; only the exports draw it into the image itself
	if (!siril_plot_draw_with(cr, spl_data, width, height, FALSE, SPL_DRAW_CANVAS))
		siril_log_debug("Problem while creating siril_plot\n");
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
}

/* GTK4 enter / leave: cursor management.  Replaces "enter-notify-event"
 * and "leave-notify-event"; both signals route through GtkEventControllerMotion. */
static void on_siril_plot_enter(GtkEventControllerMotion *controller,
		double x, double y, gpointer user_data) {
	(void)controller; (void)x; (void)y; (void)user_data;
	set_cursor("crosshair");
}

static void on_siril_plot_leave(GtkEventControllerMotion *controller,
		gpointer user_data) {
	(void)controller; (void)user_data;
	set_cursor_waiting(FALSE);
}

/* GTK4 motion: track cursor + drive selection/pan drag.  Replaces
 * "motion-notify-event"; (x, y) come straight from the controller. */
static void siril_plot_handle_pointer(GtkWidget *da, double x, double y) {
	if (!da) return;
	siril_plot_data *spl_data = (siril_plot_data *) g_object_get_data(G_OBJECT(da), "spl_data");
	GtkWidget *label = (GtkWidget *) g_object_get_data(G_OBJECT(da), "display_label_handle");
	if (!spl_data || !label) return;
	if (!is_inside_grid(x, y, &spl_data->pdd)) return;

	/* Show cursor coordinates in the bottom-left label. */
	double xpos, ypos;
	convert_surface_to_plot(spl_data, x, y, &xpos, &ypos);
	if (spl_data->logy) // the plot space is log10, the readout is not
		ypos = pow(10., ypos);
	gchar *labeltext = g_strdup_printf("%10g ; %10g", xpos, ypos);
	if (labeltext) gtk_label_set_text(GTK_LABEL(label), labeltext);
	g_free(labeltext);

	/* Update selection rectangle / pan offset while dragging. */
	if (spl_data->pdd.action == SELACTION_SELECTING && spl_data->zoomable) {
		spl_data->pdd.selection.w = x - spl_data->pdd.start.x;
		spl_data->pdd.selection.h = y - spl_data->pdd.start.y;
		gtk_widget_queue_draw(da);
	} else if (spl_data->pdd.action == SELACTION_MOVING && spl_data->zoomable) {
		double x1, x2, y1, y2;
		convert_surface_to_plot(spl_data, spl_data->pdd.start.x, spl_data->pdd.start.y, &x1, &y1);
		convert_surface_to_plot(spl_data, x, y, &x2, &y2);
		spl_data->pdd.datamin.x += x1 - x2;
		spl_data->pdd.datamax.x += x1 - x2;
		spl_data->pdd.datamin.y += y1 - y2;
		spl_data->pdd.datamax.y += y1 - y2;
		spl_data->pdd.start = (point){x, y};
		gtk_widget_queue_draw(da);
	}
}

/* Hover motion: bail out while any button is held so the drag gesture
 * is the single source of truth for motion-while-pressed (drag-to-select
 * and pan rely on motion-while-pressed, which macOS routes exclusively
 * to GtkGestureDrag — the motion controller never fires there during a
 * press). */
static void on_siril_plot_motion(GtkEventControllerMotion *controller,
		double x, double y, gpointer user_data) {
	GdkModifierType state = gtk_event_controller_get_current_event_state(
	    GTK_EVENT_CONTROLLER(controller));
	if (state & (GDK_BUTTON1_MASK | GDK_BUTTON2_MASK | GDK_BUTTON3_MASK))
		return;
	siril_plot_handle_pointer((GtkWidget *)user_data, x, y);
}

/* Drag-update: only path that fires on macOS while the button is held. */
static void on_siril_plot_drag_update(GtkGestureDrag *g,
		double offset_x, double offset_y, gpointer user_data) {
	double sx, sy;
	gtk_gesture_drag_get_start_point(g, &sx, &sy);
	siril_plot_handle_pointer((GtkWidget *)user_data, sx + offset_x, sy + offset_y);
}

/* GTK4 scroll: Cmd/Ctrl + scroll zooms in/out around the cursor.
 * Replaces "scroll-event"; (dx, dy) come straight from the controller. */
static gboolean on_siril_plot_scroll(GtkEventControllerScroll *controller,
		double dx, double dy, gpointer user_data) {
	(void)dx;
	GtkWidget *da = (GtkWidget *) user_data;
	if (!da) return FALSE;
	siril_plot_data *spl_data = (siril_plot_data *) g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data || !spl_data->zoomable) return FALSE;

	GdkModifierType mods = gtk_event_controller_get_current_event_state(
	    GTK_EVENT_CONTROLLER(controller));
	if (!(mods & get_primary())) return FALSE;

	/* GtkEventControllerScroll only knows about deltas; for the
	 * cursor-anchored zoom we need the cursor (x, y) as well — query
	 * it from the underlying GdkEvent. */
	GdkEvent *ev = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(controller));
	double sx = 0, sy = 0;
	if (ev) gdk_event_get_position(ev, &sx, &sy);

	gboolean handled = FALSE;
	if (dy > 0) handled = update_zoom(spl_data, sx, sy, SIRIL_PLOT_ZOOM_OUT);
	else if (dy < 0) handled = update_zoom(spl_data, sx, sy, SIRIL_PLOT_ZOOM_IN);
	if (handled) {
		spl_data->autotic = TRUE;
		gtk_widget_queue_draw(da);
	}
	return handled;
}

/* GTK4 press: secondary opens the popover menu; primary starts a
 * selection rectangle (or, with the primary modifier held, starts a
 * pan); double-click resets the zoom. */
static void on_siril_plot_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	GtkWidget *da = (GtkWidget *) user_data;
	if (!da) return;
	siril_plot_data *spl_data = (siril_plot_data *) g_object_get_data(G_OBJECT(da), "spl_data");
	GtkWidget *menu = (GtkWidget *) g_object_get_data(G_OBJECT(da), "menu_handle");
	if (!spl_data || !menu) return;

	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	if (button == GDK_BUTTON_SECONDARY) {
		GdkRectangle rect = { (int)x, (int)y, 0, 0 };
		gtk_popover_set_pointing_to(GTK_POPOVER(menu), &rect);
		gtk_popover_popup(GTK_POPOVER(menu));
		gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
		return;
	}
	if (!spl_data->zoomable) return;

	if (button == GDK_BUTTON_PRIMARY) {
		if (n_press == 2) {
			reset_zoom(spl_data);
			reset_selection(&spl_data->pdd);
			spl_data->autotic = TRUE;
			gtk_widget_queue_draw(da);
			gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
		} else if (spl_data->pdd.action == SELACTION_NONE
		           && is_inside_grid(x, y, &spl_data->pdd)) {
			GdkModifierType mods = gtk_event_controller_get_current_event_state(
			    GTK_EVENT_CONTROLLER(gesture));
			if (mods & get_primary()) {
				spl_data->pdd.action = SELACTION_MOVING;
				spl_data->autotic = FALSE;
				set_cursor("all-scroll");
			} else {
				spl_data->pdd.action = SELACTION_SELECTING;
				spl_data->pdd.selection = (rectangled){x, y, 0., 0.};
			}
			spl_data->pdd.start = (point){x, y};
		}
	}
}

/* GTK4 release: commit the in-progress selection/pan drag. */
static void on_siril_plot_released(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)n_press; (void)x; (void)y;
	GtkWidget *da = (GtkWidget *) user_data;
	if (!da) return;
	siril_plot_data *spl_data = (siril_plot_data *) g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data || !spl_data->zoomable) return;
	if (spl_data->pdd.action == SELACTION_SELECTING) {
		if (fabs(spl_data->pdd.selection.w) > 1. && fabs(spl_data->pdd.selection.h) > 1.) {
			double x1, x2, y1, y2;
			convert_surface_to_plot(spl_data, spl_data->pdd.selection.x,
			                        spl_data->pdd.selection.y, &x1, &y1);
			convert_surface_to_plot(spl_data,
			                        spl_data->pdd.selection.x + spl_data->pdd.selection.w,
			                        spl_data->pdd.selection.y + spl_data->pdd.selection.h, &x2, &y2);
			spl_data->pdd.datamin.x = max(min(x1, x2), spl_data->datamin.x);
			spl_data->pdd.datamax.x = min(max(x1, x2), spl_data->datamax.x);
			spl_data->pdd.datamin.y = max(min(y1, y2), spl_data->datamin.y);
			spl_data->pdd.datamax.y = min(max(y1, y2), spl_data->datamax.y);
		}
		reset_selection(&spl_data->pdd);
		gtk_widget_queue_draw(da);
	} else if (spl_data->pdd.action == SELACTION_MOVING) {
		reset_selection(&spl_data->pdd);
		spl_data->autotic = TRUE;
		set_cursor("crosshair");
		gtk_widget_queue_draw(da);
	}
}

static void on_siril_plot_grid_toggled(GtkCheckButton *checkmenuitem, gpointer user_data) {
	GtkWidget *da = (GtkWidget *)(user_data);
	if (!da)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data)
		return;
	spl_data->cfgplot.grid = (siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(checkmenuitem)))) ? GRID_ALL : 0;
	gtk_widget_queue_draw(da);
}

static void on_siril_plot_legend_toggled(GtkCheckButton *checkmenuitem, gpointer user_data) {
	GtkWidget *da = (GtkWidget *)(user_data);
	if (!da)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data)
		return;
	spl_data->show_legend = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(checkmenuitem)));
	gtk_widget_queue_draw(da);
}

// generic save callback
// the name of the caller is retrieved with gtk_widget_get_name
// the appropriate save procedure is then called accordingly
static void on_siril_plot_save_activate(GtkWidget *menuitem, gpointer user_data) {
	gchar *filename = NULL, *outname = NULL;
	GtkWidget *da = (GtkWidget *)(user_data);
	if (!da)
		return;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data)
		return;
	// the save dialog needs the top-level window, which may host several panes
	GtkWindow *window = GTK_WINDOW(gtk_widget_get_root(da));
	int width =  gtk_widget_get_width(da);
	int height = gtk_widget_get_height(da);

	const gchar *widgetname = gtk_widget_get_name(GTK_WIDGET(menuitem));
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!g_strcmp0(widgetname, "png")) {
		filename = build_save_filename(spl_data->savename, ".png", spl_data->forsequence, TRUE);
		if ((outname = save_siril_plot_dialog(window, filename, _("PNG files (*.png)"), "*.png")))
			siril_plot_save_png(spl_data, outname, width, height);
	} else if (!g_strcmp0(widgetname, "dat")) {
		filename = build_save_filename(spl_data->savename, ".dat", spl_data->forsequence, FALSE);
		if ((outname = save_siril_plot_dialog(window, filename, _("DAT files (*.dat)"), "*.dat")))
			siril_plot_save_dat(spl_data, outname, FALSE);
	} else if (!g_strcmp0(widgetname, "cb")) {
		save_siril_plot_to_clipboard(spl_data, width, height);
	}
#ifdef CAIRO_HAS_SVG_SURFACE
	else if (!g_strcmp0(widgetname, "svg")) {
		filename = build_save_filename(spl_data->savename, ".svg", spl_data->forsequence, TRUE);
		if ((outname = save_siril_plot_dialog(window, filename, _("SVG files (*.svg)"), "*.svg")))
			siril_plot_save_svg(spl_data, outname, width, height);
	}
#endif
	g_free(filename);
	g_free(outname);
}

// A GtkButton centers its label by default; for a menu-style entry inside the
// popover we want the text flush left, so left-align the internal label.
static GtkWidget *spl_menu_button_new(const char *text) {
	GtkWidget *button = gtk_button_new_with_label(text);
	GtkWidget *child = gtk_button_get_child(GTK_BUTTON(button));
	if (GTK_IS_LABEL(child)) {
		gtk_label_set_xalign(GTK_LABEL(child), 0.0);
		gtk_widget_set_halign(child, GTK_ALIGN_FILL);
	}
	return button;
}

// ---- summary widgets: tiles, metrics, pane header ------------------------

// a plain-text label, left-aligned, carrying a CSS class
static GtkWidget *spl_text_label(const gchar *text, const char *css_class) {
	GtkWidget *label = gtk_label_new((text) ? text : "");
	gtk_label_set_xalign(GTK_LABEL(label), 0.);
	gtk_widget_set_halign(label, GTK_ALIGN_START);
	if (css_class)
		gtk_widget_add_css_class(label, css_class);
	return label;
}

// a label carrying pango markup, or NULL when there is nothing to show
static GtkWidget *spl_markup_label(const gchar *markup, const char *css_class, gboolean centered) {
	if (!markup || !*markup)
		return NULL;
	GtkWidget *label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(label), markup);
	gtk_label_set_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_justify(GTK_LABEL(label), (centered) ? GTK_JUSTIFY_CENTER : GTK_JUSTIFY_LEFT);
	gtk_label_set_xalign(GTK_LABEL(label), (centered) ? 0.5 : 0.);
	gtk_widget_set_halign(label, (centered) ? GTK_ALIGN_CENTER : GTK_ALIGN_START);
	if (css_class)
		gtk_widget_add_css_class(label, css_class);
	return label;
}

// the strip of tiles summarizing a whole group, above its plots. A flow box so
// that the tiles wrap instead of forcing the window wider than the screen
static GtkWidget *build_tiles_strip(GList *tiles) {
	if (!tiles)
		return NULL;
	GtkIconTheme *theme = gtk_icon_theme_get_for_display(gdk_display_get_default());
	GtkWidget *flow = gtk_flow_box_new();
	gtk_flow_box_set_selection_mode(GTK_FLOW_BOX(flow), GTK_SELECTION_NONE);
	gtk_flow_box_set_homogeneous(GTK_FLOW_BOX(flow), TRUE);
	gtk_flow_box_set_max_children_per_line(GTK_FLOW_BOX(flow), g_list_length(tiles));
	gtk_flow_box_set_column_spacing(GTK_FLOW_BOX(flow), 8);
	gtk_flow_box_set_row_spacing(GTK_FLOW_BOX(flow), 8);
	for (GList *l = tiles; l; l = l->next) {
		spltile *tile = (spltile *)l->data;
		GtkWidget *card = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
		gtk_widget_add_css_class(card, "siril-plot-tile");
		// an unknown icon name would show up as a broken-image placeholder
		if (tile->icon && gtk_icon_theme_has_icon(theme, tile->icon)) {
			GtkWidget *icon = gtk_image_new_from_icon_name(tile->icon);
			gtk_widget_set_valign(icon, GTK_ALIGN_START);
			gtk_box_append(GTK_BOX(card), icon);
		}
		GtkWidget *col = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
		gtk_widget_set_hexpand(col, TRUE);
		gtk_box_append(GTK_BOX(col), spl_text_label(tile->label, "siril-plot-tile-label"));
		gtk_box_append(GTK_BOX(col), spl_text_label(tile->value, "siril-plot-tile-value"));
		if (tile->sublabel) {
			GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
			gtk_widget_set_margin_top(sep, 4);
			gtk_widget_set_margin_bottom(sep, 4);
			gtk_box_append(GTK_BOX(col), sep);
			gtk_box_append(GTK_BOX(col), spl_text_label(tile->sublabel, "siril-plot-tile-label"));
			gtk_box_append(GTK_BOX(col), spl_text_label(tile->subvalue, "siril-plot-tile-subvalue"));
		}
		gtk_box_append(GTK_BOX(card), col);
		gtk_flow_box_append(GTK_FLOW_BOX(flow), card);
	}
	return flow;
}

static void on_metrics_copy_clicked(GtkButton *button, gpointer user_data) {
	siril_plot_data *spl_data = (siril_plot_data *)user_data;
	gchar *text = siril_plot_metrics_to_string(spl_data, "\t");
	if (!text)
		return;
	gdk_clipboard_set_text(gtk_widget_get_clipboard(GTK_WIDGET(button)), text);
	siril_log_message(_("Plot figures were copied to the clipboard.\n"));
	g_free(text);
}

// the strip of key figures below a plot, with a button to copy them all
static GtkWidget *build_metrics_strip(siril_plot_data *spl_data) {
	if (!spl_data->metrics)
		return NULL;
	GtkWidget *strip = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
	gtk_widget_add_css_class(strip, "siril-plot-metrics");
	for (GList *l = spl_data->metrics; l; l = l->next) {
		splmetric *metric = (splmetric *)l->data;
		GtkWidget *item = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
		gtk_box_append(GTK_BOX(item), spl_text_label(metric->label, "siril-plot-metric-label"));
		gtk_box_append(GTK_BOX(item), spl_text_label(metric->value, "siril-plot-metric-value"));
		gtk_box_append(GTK_BOX(strip), item);
	}
	GtkWidget *spacer = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_set_hexpand(spacer, TRUE);
	gtk_box_append(GTK_BOX(strip), spacer);
	GtkWidget *copy = gtk_button_new_from_icon_name("edit-copy-symbolic");
	gtk_widget_set_tooltip_text(copy, _("Copy these figures to the clipboard"));
	gtk_widget_add_css_class(copy, "flat");
	g_signal_connect(copy, "clicked", G_CALLBACK(on_metrics_copy_clicked), spl_data);
	gtk_box_append(GTK_BOX(strip), copy);
	return strip;
}

static void on_siril_plot_log_toggled(GObject *sw, GParamSpec *pspec, gpointer user_data) {
	(void)pspec;
	GtkWidget *da = (GtkWidget *)user_data;
	siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(da), "spl_data");
	if (!spl_data)
		return;
	siril_plot_set_logscale(spl_data, gtk_switch_get_active(GTK_SWITCH(sw)));
	gtk_widget_queue_draw(da);
}

// solo mode: hide every sibling pane (and the separators between them) so this
// one gets the whole window
static void on_siril_plot_expand_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *card = (GtkWidget *)user_data;
	GtkWidget *row = gtk_widget_get_parent(card);
	if (!row)
		return;
	gboolean solo = gtk_toggle_button_get_active(button);
	for (GtkWidget *child = gtk_widget_get_first_child(row); child; child = gtk_widget_get_next_sibling(child))
		gtk_widget_set_visible(child, !solo || child == card);
	gtk_button_set_icon_name(GTK_BUTTON(button), (solo) ? "view-restore-symbolic" : "view-fullscreen-symbolic");
}

// the pane toolbar button opens the same popover as a right click does
static void on_siril_plot_menu_clicked(GtkButton *button, gpointer user_data) {
	(void)button;
	GtkWidget *da = (GtkWidget *)user_data;
	GtkWidget *menu = (GtkWidget *)g_object_get_data(G_OBJECT(da), "menu_handle");
	if (!menu)
		return;
	// the popover is parented to the canvas, so point it at its top right
	// corner, just below the button that opened it
	GdkRectangle rect = { gtk_widget_get_width(da) - 1, 0, 1, 1 };
	gtk_popover_set_pointing_to(GTK_POPOVER(menu), &rect);
	gtk_popover_popup(GTK_POPOVER(menu));
}

// what the plot is (caption, title, subtitle) on the left, its tools on the right
static GtkWidget *build_pane_header(siril_plot_data *spl_data, GtkWidget *da, GtkWidget *card,
		gboolean with_caption, gboolean with_expand) {
	GtkWidget *header = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *titles = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_widget_set_hexpand(titles, TRUE);

	/* Producers that predate the subtitle pack every line into the title:
	 * show its first line as the name and the rest as details. Skipped when
	 * the title holds markup, as a tag could span the line break. */
	const gchar *title = spl_data->title, *subtitle = spl_data->subtitle;
	gchar *first_line = NULL, *other_lines = NULL;
	if (title && !subtitle && !g_strstr_len(title, -1, "<")) {
		const gchar *newline = g_strstr_len(title, -1, "\n");
		if (newline) {
			title = first_line = g_strndup(title, newline - title);
			subtitle = other_lines = g_strdup(newline + 1);
		}
	}

	GtkWidget *label;
	if (with_caption && (label = spl_markup_label(spl_data->caption, "siril-plot-caption", FALSE)))
		gtk_box_append(GTK_BOX(titles), label);
	if ((label = spl_markup_label(title, "siril-plot-pane-title", FALSE)))
		gtk_box_append(GTK_BOX(titles), label);
	if ((label = spl_markup_label(subtitle, "siril-plot-pane-subtitle", FALSE)))
		gtk_box_append(GTK_BOX(titles), label);
	g_free(first_line);
	g_free(other_lines);
	gtk_box_append(GTK_BOX(header), titles);

	GtkWidget *tools = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_valign(tools, GTK_ALIGN_START);
	if (siril_plot_can_logscale(spl_data)) {
		GtkWidget *sw = gtk_switch_new();
		gtk_switch_set_active(GTK_SWITCH(sw), spl_data->logy);
		gtk_widget_set_valign(sw, GTK_ALIGN_CENTER);
		gtk_widget_set_tooltip_text(sw, _("Display the y axis in log scale"));
		g_signal_connect(sw, "notify::active", G_CALLBACK(on_siril_plot_log_toggled), da);
		gtk_box_append(GTK_BOX(tools), sw);
		gtk_box_append(GTK_BOX(tools), spl_text_label(_("Log scale"), NULL));
	}
	if (with_expand) {
		GtkWidget *expand = gtk_toggle_button_new();
		gtk_button_set_icon_name(GTK_BUTTON(expand), "view-fullscreen-symbolic");
		gtk_widget_set_tooltip_text(expand, _("Show this plot alone"));
		gtk_widget_add_css_class(expand, "flat");
		g_signal_connect(expand, "toggled", G_CALLBACK(on_siril_plot_expand_toggled), card);
		gtk_box_append(GTK_BOX(tools), expand);
	}
	GtkWidget *menu_button = gtk_button_new_from_icon_name("view-more-symbolic");
	gtk_widget_set_tooltip_text(menu_button, _("Plot options and exports"));
	gtk_widget_add_css_class(menu_button, "flat");
	g_signal_connect(menu_button, "clicked", G_CALLBACK(on_siril_plot_menu_clicked), da);
	gtk_box_append(GTK_BOX(tools), menu_button);
	gtk_box_append(GTK_BOX(header), tools);
	return header;
}

// ---- window-level exports ------------------------------------------------

// renders a widget exactly as displayed into a png file
static gboolean snapshot_widget_to_png(GtkWidget *widget, const char *filename) {
	GtkWidget *parent = gtk_widget_get_parent(widget);
	int w = gtk_widget_get_width(widget), h = gtk_widget_get_height(widget);
	graphene_rect_t bounds;
	if (!parent || w <= 0 || h <= 0 || !gtk_widget_compute_bounds(widget, parent, &bounds))
		return FALSE;
	/* A GtkWidgetPaintable of a widget deeper than the window's own child
	 * snapshots empty, so the parent is asked to render it instead. The node
	 * it hands back sits in the parent's coordinates. */
	GtkSnapshot *snapshot = gtk_snapshot_new();
	gtk_widget_snapshot_child(parent, widget, snapshot);
	GskRenderNode *node = gtk_snapshot_free_to_node(snapshot);
	if (!node)
		return FALSE;

	gboolean success = FALSE;
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, w, h);
	if (!cairo_surface_status(surface)) {
		cairo_t *cr = cairo_create(surface);
		// the window background belongs to the window node, not to this child
		cairo_set_source_rgb(cr, 1., 1., 1.);
		cairo_paint(cr);
		cairo_translate(cr, -bounds.origin.x, -bounds.origin.y);
		gsk_render_node_draw(node, cr);
		cairo_destroy(cr);
		success = (cairo_surface_write_to_png(surface, filename) == CAIRO_STATUS_SUCCESS);
	}
	cairo_surface_destroy(surface);
	gsk_render_node_unref(node);
	return success;
}

static void on_export_snapshot_clicked(GtkWidget *button, gpointer user_data) {
	(void)button;
	GtkWindow *window = GTK_WINDOW(user_data);
	// the action bar itself has no place in a report
	GtkWidget *report = (GtkWidget *)g_object_get_data(G_OBJECT(window), "report_widget");
	if (!report)
		return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	gchar *filename = build_save_filename("plots", ".png", FALSE, TRUE);
	gchar *outname = save_siril_plot_dialog(window, filename, _("PNG files (*.png)"), "*.png");
	if (outname) {
		if (snapshot_widget_to_png(report, outname))
			siril_log_message(_("%s has been saved.\n"), outname);
		else
			siril_log_message(_("Could not save the report snapshot.\n"));
	}
	g_free(filename);
	g_free(outname);
}

// two panes of a group may share a savename, and the timestamp only goes down
// to the second: make sure a plot never overwrites the one saved just before
static gchar *unique_plot_path(const gchar *folder, const gchar *name) {
	gchar *path = g_build_filename(folder, name, NULL);
	if (!g_file_test(path, G_FILE_TEST_EXISTS))
		return path;
	gchar *base = g_strdup(name);
	gchar *dot = g_strrstr(base, ".");
	gchar *ext = g_strdup((dot) ? dot : "");
	if (dot)
		*dot = '\0';
	for (int i = 2; i < 100; i++) {
		g_free(path);
		gchar *candidate = g_strdup_printf("%s_%d%s", base, i, ext);
		path = g_build_filename(folder, candidate, NULL);
		g_free(candidate);
		if (!g_file_test(path, G_FILE_TEST_EXISTS))
			break;
	}
	g_free(base);
	g_free(ext);
	return path;
}

// saves every pane of the window as a png in a folder chosen by the user
static void on_save_plots_clicked(GtkButton *button, gpointer user_data) {
	(void)button;
	GtkWindow *window = GTK_WINDOW(user_data);
	GList *das = (GList *)g_object_get_data(G_OBJECT(window), "da_list");
	if (!das)
		return;

	SirilFileChooser *fc = siril_fc_open(window, GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
	siril_fc_set_current_folder_path(fc, com.wd);
	gchar *folder = NULL;
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT)
		folder = siril_fc_get_filename(fc);
	siril_fc_destroy(fc);
	if (!folder)
		return;

	control_window_switch_to_tab(OUTPUT_LOGS);
	int saved = 0;
	for (GList *l = das; l; l = l->next) {
		GtkWidget *da = (GtkWidget *)l->data;
		siril_plot_data *spl_data = (siril_plot_data *)g_object_get_data(G_OBJECT(da), "spl_data");
		if (!spl_data)
			continue;
		gchar *name = build_save_filename((spl_data->savename) ? spl_data->savename : "plot",
				".png", spl_data->forsequence, TRUE);
		gchar *path = unique_plot_path(folder, name);
		if (siril_plot_save_png(spl_data, path, gtk_widget_get_width(da), gtk_widget_get_height(da)))
			saved++;
		g_free(name);
		g_free(path);
	}
	siril_log_message(_("%d plot(s) saved in %s\n"), saved, folder);
	g_free(folder);
}

// the bar at the bottom of every siril-plot window
static GtkWidget *build_action_bar(GtkWidget *window) {
	GtkWidget *bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_widget_add_css_class(bar, "siril-plot-actionbar");

	// exports the summary and every plot as one image, as displayed
	GtkWidget *export = gtk_button_new_with_label(_("Export PNG"));
	gtk_widget_set_tooltip_text(export, _("Save the whole window as a PNG image"));
	g_signal_connect(export, "clicked", G_CALLBACK(on_export_snapshot_clicked), window);
	gtk_box_append(GTK_BOX(bar), export);

	GtkWidget *spacer = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_set_hexpand(spacer, TRUE);
	gtk_box_append(GTK_BOX(bar), spacer);

	GtkWidget *save = gtk_button_new_with_label(_("Save plots"));
	g_signal_connect(save, "clicked", G_CALLBACK(on_save_plots_clicked), window);
	gtk_box_append(GTK_BOX(bar), save);

	GtkWidget *close = gtk_button_new_with_label(_("Close"));
	g_signal_connect_swapped(close, "clicked", G_CALLBACK(gtk_window_close), window);
	gtk_box_append(GTK_BOX(bar), close);
	return bar;
}

// creates the bare siril-plot window (chrome only, no plot pane yet): title,
// CSS class, close-request handler and margins. Shared by the single-plot
// and grouped-plot window builders below.
static GtkWidget *create_siril_plot_window_shell(const gchar *title) {
	GtkWidget *window = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(window), (title) ? title : "Siril plot");

	/* Uniform plot-window background.  Class-scoped to avoid styling
	 * other windows; provider registered display-wide once. */
	/* The window is a light grey board on which the plots sit as white
	 * cards; sizes are in pt to match the text drawn inside the plots,
	 * which pango renders in points too. */
	siril_register_css_for_display("siril-plot-window",
		"window.siril-plot { color: #303030; background: #f2f2f4; font-size: 12px; }"
		"window.siril-plot label.siril-plot-caption { font-size: 13pt; font-weight: bold; color: #303030; }"
		"window.siril-plot .siril-plot-card,"
		"window.siril-plot .siril-plot-tile { background: white; border: 1px solid #dcdcdc; border-radius: 8px; padding: 8px; }"
		"window.siril-plot .siril-plot-pane-title { font-size: 12pt; font-weight: bold; color: #303030; }"
		"window.siril-plot .siril-plot-pane-subtitle { font-size: 10pt; color: #707070; }"
		"window.siril-plot .siril-plot-tile-label { font-size: 9pt; color: #808080; }"
		"window.siril-plot .siril-plot-tile-value { font-size: 13pt; font-weight: bold; color: #303030; }"
		"window.siril-plot .siril-plot-tile-subvalue { color: #505050; }"
		"window.siril-plot .siril-plot-metrics { border-top: 1px solid #ececec; padding-top: 4px; }"
		"window.siril-plot .siril-plot-metric-label { font-size: 9pt; color: #808080; }"
		"window.siril-plot .siril-plot-metric-value { font-weight: bold; color: #303030; }"
		"window.siril-plot .siril-plot-actionbar { padding-top: 4px; }"
		/* The window forces a light background, so the controls have to be
		 * dressed for it too: left to the theme they would come out dark. */
		"window.siril-plot button { background-image: none; background-color: #fbfbfb; color: #303030; border: 1px solid #d0d0d0; }"
		"window.siril-plot button:hover { background-color: #f0f0f0; }"
		"window.siril-plot button:active, window.siril-plot button:checked { background-color: #e2e2e6; }"
		"window.siril-plot button.flat { background-color: transparent; border-color: transparent; }"
		"window.siril-plot button.flat:hover { background-color: #ececee; }"
		"window.siril-plot switch { background-color: #d0d0d6; border: 1px solid #c0c0c6; }"
		"window.siril-plot switch:checked { background-color: #3584e4; border-color: #2b6fc4; }"
		"window.siril-plot switch > slider { background-color: white; border: 1px solid #c8c8ce; }"
		"window.siril-plot popover > contents { background-color: white; color: #303030; }"
		"window.siril-plot popover > arrow { background-color: white; border: 1px solid #dcdcdc; }"
		"window.siril-plot check { background-color: white; border: 1px solid #b4b4ba; color: #303030; }"
		"window.siril-plot check:checked { background-color: #3584e4; border-color: #2b6fc4; color: white; }");
	gtk_widget_add_css_class(GTK_WIDGET(window), "siril-plot");
	// connect the delete-event signal, triggered when the window is closed
	// the callback frees every spl_data displayed in the window
	g_signal_connect(G_OBJECT(window), "close-request", G_CALLBACK(on_siril_plot_window_closed), NULL);
	gtk_widget_set_margin_start(GTK_WIDGET(window), 5); gtk_widget_set_margin_end(GTK_WIDGET(window), 5); gtk_widget_set_margin_top(GTK_WIDGET(window), 5); gtk_widget_set_margin_bottom(GTK_WIDGET(window), 5);
	return window;
}

// builds one interactive plot pane (drawing area + coordinate label + its
// right-click save/grid/legend popover) bound to spl_data, and registers
// spl_data and the popover on `window` so on_siril_plot_window_closed can
// free/unparent them regardless of how many panes the window ends up
// holding. Returns the pane's root widget, to be packed by the caller.
// `with_caption` asks the pane to show spl_data->caption itself, which a group
// window does not want as it displays it once for all its panes; `with_expand`
// adds the button that gives this pane the whole window.
static GtkWidget *build_siril_plot_pane(GtkWidget *window, siril_plot_data *spl_data,
		gboolean with_caption, gboolean with_expand) {
	GtkWidget *vbox, *da, *label, *menu;
	GtkWidget *spl_menu_grid, *spl_menu_legend;
	GtkWidget *spl_menu_save_cb, *spl_menu_save_png, *spl_menu_save_dat, *spl_menu_sep;
#ifdef CAIRO_HAS_SVG_SURFACE
	GtkWidget *spl_menu_save_svg;
#endif

	//prepare interactivity for spl_data
	reset_zoom(spl_data);
	reset_selection(&spl_data->pdd);
	spl_data->zoomable = spl_data->bkg == NULL;

	// register spl_data on the window for freeing when it is closed
	GList *spl_data_list = (GList *)g_object_get_data(G_OBJECT(window), "spl_data_list");
	spl_data_list = g_list_append(spl_data_list, spl_data);
	g_object_set_data(G_OBJECT(window), "spl_data_list", spl_data_list);

	// the pane is a card: header, canvas, cursor readout, key figures
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_widget_add_css_class(vbox, "siril-plot-card");

	// add the drawing area
	da = gtk_drawing_area_new();
	int width = (!spl_data->width) ? SIRIL_PLOT_DISPLAY_WIDTH : spl_data->width;
	int height = (!spl_data->height) ? SIRIL_PLOT_DISPLAY_WIDTH : spl_data->height;
	gtk_widget_set_size_request(da, width, height);
	gtk_widget_set_hexpand(da, TRUE);
	gtk_widget_set_vexpand(da, TRUE);
	gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(da), on_siril_plot_draw, da, NULL);

	// attach the spl_data to its own drawing area, so several panes can coexist in one window
	g_object_set_data(G_OBJECT(da), "spl_data", spl_data);
	// the window keeps its panes in order, for "save plots" and the reports
	GList *da_list = (GList *)g_object_get_data(G_OBJECT(window), "da_list");
	da_list = g_list_append(da_list, da);
	g_object_set_data(G_OBJECT(window), "da_list", da_list);

	gtk_box_append(GTK_BOX(vbox), build_pane_header(spl_data, da, vbox, with_caption, with_expand));
	gtk_box_append(GTK_BOX(vbox), da);

	/* GTK4 event-controllers on the drawing area: replace the GTK3
	 * enter-notify, leave-notify, motion-notify, button-press,
	 * button-release, and scroll-event signals. */
	GtkEventController *motion_ctl = gtk_event_controller_motion_new();
	g_signal_connect(motion_ctl, "enter",  G_CALLBACK(on_siril_plot_enter),  da);
	g_signal_connect(motion_ctl, "leave",  G_CALLBACK(on_siril_plot_leave),  da);
	g_signal_connect(motion_ctl, "motion", G_CALLBACK(on_siril_plot_motion), da);
	gtk_widget_add_controller(da, motion_ctl);

	GtkGesture *click_gst = gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_gst), 0); /* any button */
	g_signal_connect(click_gst, "pressed",  G_CALLBACK(on_siril_plot_pressed),  da);
	g_signal_connect(click_gst, "released", G_CALLBACK(on_siril_plot_released), da);
	gtk_widget_add_controller(da, GTK_EVENT_CONTROLLER(click_gst));

	/* Drag: primary-only (selection rectangle, pan); see also the
	 * motion controller above for hover-only label updates. */
	GtkGesture *drag_gst = gtk_gesture_drag_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_gst),
	                              GDK_BUTTON_PRIMARY);
	g_signal_connect(drag_gst, "drag-update",
	                 G_CALLBACK(on_siril_plot_drag_update), da);
	gtk_widget_add_controller(da, GTK_EVENT_CONTROLLER(drag_gst));

	/* Group click + drag so neither claim denies the other. */
	gtk_gesture_group(click_gst, drag_gst);

	GtkEventController *scroll_ctl = gtk_event_controller_scroll_new(
	    GTK_EVENT_CONTROLLER_SCROLL_VERTICAL);
	g_signal_connect(scroll_ctl, "scroll", G_CALLBACK(on_siril_plot_scroll), da);
	gtk_widget_add_controller(da, scroll_ctl);

	// add the label
	label = gtk_label_new("0;0");
	gtk_widget_set_halign(label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(label, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(vbox), label);
	gtk_widget_set_halign(label, GTK_ALIGN_START);
	// and cache its handle
	g_object_set_data(G_OBJECT(da), "display_label_handle", label);

	// the key figures of this plot, when its producer supplied any
	GtkWidget *metrics = build_metrics_strip(spl_data);
	if (metrics)
		gtk_box_append(GTK_BOX(vbox), metrics);

	/* Phase 15: GtkMenu / GtkMenuItem / gtk_separator_menu_item_new are
	 * gone in GTK4.  Build the right-click context menu as a GtkPopover
	 * with a vertical GtkBox of GtkCheckButtons (grid/legend) and
	 * GtkButtons (save actions); the existing toggled / clicked handlers
	 * keep working because they already accept GtkCheckButton* /
	 * GtkWidget* respectively. Parented to `da` (not the window) so its
	 * pointing-to rectangle, expressed in da's own coordinates, lines up
	 * even when the window hosts several panes. */
	menu = gtk_popover_new();
	gtk_widget_set_parent(menu, da);
	GtkWidget *menu_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_widget_set_margin_start(menu_box, 6);
	gtk_widget_set_margin_end(menu_box, 6);
	gtk_widget_set_margin_top(menu_box, 6);
	gtk_widget_set_margin_bottom(menu_box, 6);
	gtk_popover_set_child(GTK_POPOVER(menu), menu_box);

	// grid
	spl_menu_grid = gtk_check_button_new_with_label("Grid");
	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(spl_menu_grid)), TRUE);
	g_signal_connect(G_OBJECT(spl_menu_grid), "toggled", G_CALLBACK(on_siril_plot_grid_toggled), da);

	// legend
	spl_menu_legend = gtk_check_button_new_with_label("Legend");
	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(spl_menu_legend)), TRUE);
	g_signal_connect(G_OBJECT(spl_menu_legend), "toggled", G_CALLBACK(on_siril_plot_legend_toggled), da);

	// sep
	spl_menu_sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);

	// save to clipboard
	spl_menu_save_cb = spl_menu_button_new("Save view to Clipboard");
	gtk_widget_set_name(spl_menu_save_cb, "cb");
	g_signal_connect(G_OBJECT(spl_menu_save_cb), "clicked", G_CALLBACK(on_siril_plot_save_activate), da);

	// save as png
	spl_menu_save_png = spl_menu_button_new("Save view as PNG");
	gtk_widget_set_name(spl_menu_save_png, "png");
	g_signal_connect(G_OBJECT(spl_menu_save_png), "clicked", G_CALLBACK(on_siril_plot_save_activate), da);

#ifdef CAIRO_HAS_SVG_SURFACE
	// save as svg
	spl_menu_save_svg = spl_menu_button_new("Save view as SVG");
	gtk_widget_set_name(spl_menu_save_svg, "svg");
	g_signal_connect(G_OBJECT(spl_menu_save_svg), "clicked", G_CALLBACK(on_siril_plot_save_activate), da);
#endif

	// save as dat
	spl_menu_save_dat = spl_menu_button_new("Export dat file");
	gtk_widget_set_name(spl_menu_save_dat, "dat");
	g_signal_connect(G_OBJECT(spl_menu_save_dat), "clicked", G_CALLBACK(on_siril_plot_save_activate), da);

	gtk_box_append(GTK_BOX(menu_box), spl_menu_grid);
	gtk_box_append(GTK_BOX(menu_box), spl_menu_legend);
	gtk_box_append(GTK_BOX(menu_box), spl_menu_sep);
	gtk_box_append(GTK_BOX(menu_box), spl_menu_save_cb);
	gtk_box_append(GTK_BOX(menu_box), spl_menu_save_png);
#ifdef CAIRO_HAS_SVG_SURFACE
	gtk_box_append(GTK_BOX(menu_box), spl_menu_save_svg);
#endif
	gtk_box_append(GTK_BOX(menu_box), spl_menu_save_dat);
	// and cache its handle
	g_object_set_data(G_OBJECT(da), "menu_handle", menu);

	// register the popover on the window for unparenting when it is closed
	GList *menu_list = (GList *)g_object_get_data(G_OBJECT(window), "menu_list");
	menu_list = g_list_append(menu_list, menu);
	g_object_set_data(G_OBJECT(window), "menu_list", menu_list);

	return vbox;
}

gboolean create_new_siril_plot_window(gpointer p) {
	siril_plot_data *spl_data = (siril_plot_data *)p;

	// sanity checks
	if (!spl_data) {
		siril_log_debug("Passed an empty spl_data structure\n");
		return FALSE;
	}
	if (!spl_data_has_any_plot(spl_data)) {
		siril_log_debug("Trying to display plot that contains no data, freeing and aborting\n");
		free_siril_plot_data(spl_data);
		return FALSE;
	}

	GtkWidget *window = create_siril_plot_window_shell(NULL);
	GtkWidget *content = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	GtkWidget *pane = build_siril_plot_pane(window, spl_data, TRUE, FALSE);
	gtk_box_append(GTK_BOX(content), pane);
	// what the snapshot export captures: everything but the action bar
	g_object_set_data(G_OBJECT(window), "report_widget", pane);
	gtk_box_append(GTK_BOX(content), build_action_bar(window));
	gtk_window_set_child(GTK_WINDOW(window), content);

	gtk_window_present(GTK_WINDOW(window));
	gtk_widget_set_visible(window, TRUE);
	return FALSE;
}

// displays every siril_plot_data of the group side by side in a single
// window, each pane keeping its own independent zoom/pan/save interactivity.
// Takes ownership of `p` (a siril_plot_group*): the group container is
// consumed here, its items are handed off to the window for lifetime
// management (freed by on_siril_plot_window_closed when the window closes).
gboolean create_new_siril_plot_group_window(gpointer p) {
	siril_plot_group *grp = (siril_plot_group *)p;
	if (!grp) {
		siril_log_debug("Passed an empty siril_plot_group structure\n");
		return FALSE;
	}

	// keep only the plots that actually hold data, freeing the empty ones
	GList *valid = NULL;
	for (GList *l = grp->items; l; l = l->next) {
		siril_plot_data *spl_data = (siril_plot_data *)l->data;
		if (spl_data_has_any_plot(spl_data))
			valid = g_list_append(valid, spl_data);
		else {
			siril_log_debug("Trying to display plot that contains no data, freeing and skipping\n");
			free_siril_plot_data(spl_data);
		}
	}
	g_list_free(grp->items);
	grp->items = NULL;

	if (!valid) {
		siril_log_debug("No valid plots in siril_plot_group, aborting\n");
		siril_plot_free_tiles(grp->tiles);
		g_free(grp->title);
		free(grp);
		return FALSE;
	}

	GtkWidget *window = create_siril_plot_window_shell(grp->title);

	/* When all the plots carry the same caption, it is displayed once above
	 * them instead of being repeated in each pane. The panes then draw only
	 * their own title, but their png/svg exports still include the caption,
	 * as a saved plot has to stand on its own. */
	const gchar *shared_caption = ((siril_plot_data *)valid->data)->caption;
	for (GList *l = valid->next; shared_caption && l; l = l->next) {
		const gchar *caption = ((siril_plot_data *)l->data)->caption;
		if (g_strcmp0(caption, shared_caption))
			shared_caption = NULL;
	}

	GtkWidget *scroller = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller), GTK_POLICY_AUTOMATIC, GTK_POLICY_NEVER);
	gtk_widget_set_hexpand(scroller, TRUE);
	gtk_widget_set_vexpand(scroller, TRUE);
	/* A GtkScrolledWindow normally reports a small minimum/natural size of
	 * its own regardless of its child (that's how it lets arbitrarily large
	 * content scroll); left alone, the window opens too small to show every
	 * pane and the user has to resize it manually. Propagating the child's
	 * natural size makes the window open wide enough for all panes (up to
	 * what the screen allows), while still falling back to scrolling if the
	 * group holds more plots than can fit. */
	gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(scroller), TRUE);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(scroller), TRUE);

	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gboolean multi = valid->next != NULL;
	gboolean first = TRUE;
	for (GList *l = valid; l; l = l->next) {
		if (!first)
			gtk_box_append(GTK_BOX(hbox), gtk_separator_new(GTK_ORIENTATION_VERTICAL));
		GtkWidget *pane = build_siril_plot_pane(window, (siril_plot_data *)l->data, FALSE, multi);
		gtk_box_append(GTK_BOX(hbox), pane);
		first = FALSE;
	}
	g_list_free(valid);

	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroller), hbox);

	/* The report holds everything a saved report should show: the shared
	 * summary, the tiles and the plots, but not the action bar. */
	GtkWidget *content = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	GtkWidget *report = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	GtkWidget *caption = spl_markup_label(shared_caption, "siril-plot-caption", TRUE);
	if (caption) {
		gtk_widget_set_margin_top(caption, SIRIL_PLOT_MARGIN);
		gtk_box_append(GTK_BOX(report), caption);
	}
	GtkWidget *tiles = build_tiles_strip(grp->tiles);
	if (tiles)
		gtk_box_append(GTK_BOX(report), tiles);
	siril_plot_free_tiles(grp->tiles); // the widgets hold their own copies
	grp->tiles = NULL;
	gtk_box_append(GTK_BOX(report), scroller);
	gtk_box_append(GTK_BOX(content), report);
	gtk_box_append(GTK_BOX(content), build_action_bar(window));
	gtk_window_set_child(GTK_WINDOW(window), content);

	// what the snapshot export captures
	g_object_set_data(G_OBJECT(window), "report_widget", report);

	g_free(grp->title);
	free(grp); // ownership of its items was transferred to the window

	gtk_window_present(GTK_WINDOW(window));
	gtk_widget_set_visible(window, TRUE);
	return FALSE;
}
