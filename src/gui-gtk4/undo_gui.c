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

/*
 * GTK4 equivalent of src/gui/undo_gui.c.  Secondary-click (right-click) on
 * the header undo / redo button pops up a list of entries from
 * com.undo_stack / com.redo_stack; activating an entry steps undo/redo
 * that many times in one go.
 *
 * Differences from the GTK3 version:
 *   - "button-press-event" / GdkEventButton no longer exist; instead a
 *     GtkGestureClick with the secondary button bound is attached to each
 *     header button.
 *   - GtkContainer is gone: boxes use gtk_box_append, GtkListBoxRow uses
 *     gtk_list_box_row_set_child, popovers are parented to their anchor
 *     via gtk_popover_set_parent.
 *   - gtk_widget_show_all / gtk_widget_destroy are gone: widgets are
 *     visible by default, and popovers must be explicitly unparented
 *     (deferred via an idle so the "closed" emission can finish before
 *     the popover is removed from its parent).
 */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_language.h"
#include "core/undo.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/mpp_ap_editor.h"
#include "undo_gui.h"

static gboolean unparent_popover_idle(gpointer data) {
	GtkWidget *popover = GTK_WIDGET(data);
	if (gtk_widget_get_parent(popover))
		gtk_widget_unparent(popover);
	g_object_unref(popover);
	return G_SOURCE_REMOVE;
}

static void on_undo_popover_closed(GtkPopover *popover, gpointer user_data) {
	/* Defer unparent: gtk_popover_popdown emits "closed" while GTK still
	 * holds internal state on the widget, and unparenting here would
	 * re-enter GTK during dispose. */
	g_object_ref(popover);
	g_idle_add(unparent_popover_idle, popover);
}

static void on_undo_popover_row_activated(GtkListBox *box, GtkListBoxRow *row,
                                          gpointer user_data) {
	GtkPopover *popover = GTK_POPOVER(g_object_get_data(G_OBJECT(box), "popover"));
	int dir   = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "dir"));
	int level = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "level"));
	gtk_popover_popdown(popover);
	set_cursor_waiting(TRUE);
	for (int i = 0; i < level; i++) {
		if (dir == UNDO && !is_undo_available()) break;
		if (dir == REDO && !is_redo_available()) break;
		undo_display_data(dir);
	}
	set_cursor_waiting(FALSE);
}

static void show_undo_history_popover(GtkWidget *button, int dir) {
	GList *stack = (dir == UNDO) ? com.undo_stack : com.redo_stack;
	if (!stack) return;

	GtkWidget *popover = gtk_popover_new();
	/* GtkPopover in GTK4 is parented like any other widget; there's no
	 * dedicated set_parent() on GtkPopover, only gtk_widget_set_parent. */
	gtk_widget_set_parent(popover, button);
	gtk_popover_set_position(GTK_POPOVER(popover), GTK_POS_BOTTOM);
	g_signal_connect(popover, "closed", G_CALLBACK(on_undo_popover_closed), NULL);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_margin_start(vbox, 4);
	gtk_widget_set_margin_end(vbox, 4);
	gtk_widget_set_margin_top(vbox, 6);
	gtk_widget_set_margin_bottom(vbox, 4);

	const gchar *title = (dir == UNDO) ? _("Undo history") : _("Redo history");
	GtkWidget *heading = gtk_label_new(title);
	PangoAttrList *attrs = pango_attr_list_new();
	pango_attr_list_insert(attrs, pango_attr_weight_new(PANGO_WEIGHT_BOLD));
	gtk_label_set_attributes(GTK_LABEL(heading), attrs);
	pango_attr_list_unref(attrs);
	gtk_widget_set_margin_bottom(heading, 4);
	gtk_box_append(GTK_BOX(vbox), heading);

	GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_margin_bottom(sep, 2);
	gtk_box_append(GTK_BOX(vbox), sep);

	GtkWidget *scroll = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_max_content_height(GTK_SCROLLED_WINDOW(scroll), 300);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(scroll), TRUE);

	GtkWidget *listbox = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(listbox), GTK_SELECTION_NONE);
	g_object_set_data(G_OBJECT(listbox), "dir",     GINT_TO_POINTER(dir));
	g_object_set_data(G_OBJECT(listbox), "popover", popover);
	g_signal_connect(listbox, "row-activated",
	                 G_CALLBACK(on_undo_popover_row_activated), NULL);

	int n = 0;
	for (GList *l = stack; l; l = l->next, n++) {
		historic *h = (historic *)l->data;
		const gchar *label = (h->history[0] != '\0') ? h->history : _("(unnamed)");

		GtkWidget *row_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		gtk_widget_set_margin_start(row_box, 6);
		gtk_widget_set_margin_end(row_box, 6);
		gtk_widget_set_margin_top(row_box, 3);
		gtk_widget_set_margin_bottom(row_box, 3);

		GtkWidget *lbl = gtk_label_new(label);
		gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
		gtk_widget_set_hexpand(lbl, TRUE);
		gtk_box_append(GTK_BOX(row_box), lbl);

		GtkWidget *list_row = gtk_list_box_row_new();
		gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(list_row), row_box);
		g_object_set_data(G_OBJECT(list_row), "dir",   GINT_TO_POINTER(dir));
		g_object_set_data(G_OBJECT(list_row), "level", GINT_TO_POINTER(n + 1));
		gtk_list_box_append(GTK_LIST_BOX(listbox), list_row);
	}

	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroll), listbox);
	gtk_widget_set_vexpand(scroll, TRUE);
	gtk_box_append(GTK_BOX(vbox), scroll);
	gtk_popover_set_child(GTK_POPOVER(popover), vbox);
	gtk_popover_popup(GTK_POPOVER(popover));
}

static void on_undo_secondary_pressed(GtkGestureClick *gesture, int n_press,
                                      double x, double y, gpointer user_data) {
	(void) n_press; (void) x; (void) y;
	GtkWidget *button = GTK_WIDGET(user_data);
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
	/* The history popover lists the image-processing history; it isn't
	 * meaningful while Undo/Redo are diverted to the AP editor. */
	if (mpp_ap_editor_is_open())
		return;
	show_undo_history_popover(button, UNDO);
}

static void on_redo_secondary_pressed(GtkGestureClick *gesture, int n_press,
                                      double x, double y, gpointer user_data) {
	(void) n_press; (void) x; (void) y;
	GtkWidget *button = GTK_WIDGET(user_data);
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
	if (mpp_ap_editor_is_open())
		return;
	show_undo_history_popover(button, REDO);
}

static void attach_secondary_click(GtkWidget *button, GCallback handler) {
	GtkGesture *click = gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), GDK_BUTTON_SECONDARY);
	/* CAPTURE so we run before the button's own primary-button gesture
	 * gets a chance to consume the event. */
	gtk_event_controller_set_propagation_phase(GTK_EVENT_CONTROLLER(click),
	                                           GTK_PHASE_CAPTURE);
	g_signal_connect(click, "pressed", handler, button);
	gtk_widget_add_controller(button, GTK_EVENT_CONTROLLER(click));
}

void setup_undo_redo_long_press(void) {
	GtkWidget *undo_btn = lookup_widget("header_undo_button");
	GtkWidget *redo_btn = lookup_widget("header_redo_button");

	attach_secondary_click(undo_btn, G_CALLBACK(on_undo_secondary_pressed));
	attach_secondary_click(redo_btn, G_CALLBACK(on_redo_secondary_pressed));
}
