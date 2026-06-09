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

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_language.h"
#include "core/undo.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "undo_gui.h"

static gboolean destroy_widget_idle(gpointer data) {
	gtk_widget_destroy(GTK_WIDGET(data));
	g_object_unref(G_OBJECT(data));
	return G_SOURCE_REMOVE;
}

static void on_undo_popover_closed(GtkPopover *popover, gpointer user_data) {
	/* Schedule destruction via idle to avoid re-entrancy: gtk_popover_popdown
	 * emits "closed" while still holding internal GTK state on the widget. */
	g_object_ref(popover);
	g_idle_add(destroy_widget_idle, popover);
}

static void on_undo_popover_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer user_data) {
	GtkWidget *popover = GTK_WIDGET(g_object_get_data(G_OBJECT(box), "popover"));
	int dir   = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "dir"));
	int level = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "level"));
	gtk_widget_hide(GTK_WIDGET(popover));
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

	GtkWidget *popover = gtk_popover_new(button);
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
	gtk_box_pack_start(GTK_BOX(vbox), heading, FALSE, FALSE, 0);

	GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_margin_bottom(sep, 2);
	gtk_box_pack_start(GTK_BOX(vbox), sep, FALSE, FALSE, 0);

	GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_max_content_height(GTK_SCROLLED_WINDOW(scroll), 300);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(scroll), TRUE);

	GtkWidget *listbox = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(listbox), GTK_SELECTION_NONE);
	g_object_set_data(G_OBJECT(listbox), "dir",     GINT_TO_POINTER(dir));
	g_object_set_data(G_OBJECT(listbox), "popover", popover);
	g_signal_connect(listbox, "row-activated", G_CALLBACK(on_undo_popover_row_activated), NULL);

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
		gtk_box_pack_start(GTK_BOX(row_box), lbl, TRUE, TRUE, 0);

		GtkWidget *list_row = gtk_list_box_row_new();
		gtk_container_add(GTK_CONTAINER(list_row), row_box);
		g_object_set_data(G_OBJECT(list_row), "dir",   GINT_TO_POINTER(dir));
		g_object_set_data(G_OBJECT(list_row), "level", GINT_TO_POINTER(n + 1));
		gtk_list_box_insert(GTK_LIST_BOX(listbox), list_row, -1);
	}

	gtk_container_add(GTK_CONTAINER(scroll), listbox);
	gtk_box_pack_start(GTK_BOX(vbox), scroll, TRUE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(popover), vbox);
	gtk_widget_show_all(popover);
	gtk_popover_popup(GTK_POPOVER(popover));
}

static gboolean on_undo_button_press(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	if (event->button == GDK_BUTTON_SECONDARY) {
		show_undo_history_popover(widget, UNDO);
		return TRUE;
	}
	return FALSE;
}

static gboolean on_redo_button_press(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	if (event->button == GDK_BUTTON_SECONDARY) {
		show_undo_history_popover(widget, REDO);
		return TRUE;
	}
	return FALSE;
}

void setup_undo_redo_long_press(void) {
	GtkWidget *undo_btn = lookup_widget("header_undo_button");
	GtkWidget *redo_btn = lookup_widget("header_redo_button");

	g_signal_connect(undo_btn, "button-press-event", G_CALLBACK(on_undo_button_press), NULL);
	g_signal_connect(redo_btn, "button-press-event", G_CALLBACK(on_redo_button_press), NULL);
}
