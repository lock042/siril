/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"

#include "core/siril_cmd_help.h"

static GtkWidget *shortcuts_window = NULL;

static void
on_shortcuts_window_destroyed(gpointer data, GObject *where_the_object_was)
{
    GtkWidget **pwindow = data;
    *pwindow = NULL;
}

void siril_cmd_help_keyboard_shortcuts(GtkWindow *window) {
    if (shortcuts_window == NULL) {
        GtkBuilder *s_builder;

        s_builder = gtk_builder_new_from_resource("/org/siril/ui/siril-shortcuts.ui");
        shortcuts_window = GTK_WIDGET(gtk_builder_get_object(s_builder, "shortcuts-siril"));

        g_object_weak_ref(G_OBJECT(shortcuts_window), on_shortcuts_window_destroyed, &shortcuts_window);

        g_object_unref(s_builder);
    }

    if (GTK_WINDOW(window) != gtk_window_get_transient_for(GTK_WINDOW(shortcuts_window))) {
        gtk_window_set_transient_for(GTK_WINDOW(shortcuts_window), GTK_WINDOW(window));
    }

    gtk_widget_show_all(shortcuts_window);
    gtk_window_present(GTK_WINDOW(shortcuts_window));
}

