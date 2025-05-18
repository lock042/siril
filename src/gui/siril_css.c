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

#include "siril_css.h"

#define CSS_FILENAME "siril.css"

/**
 * Loads the css sheet
 */
void load_css_style_sheet () {
	GError *error;
	GBytes *css_buffer = g_resource_lookup_data(com.resource, "/org/siril/ui/siril.css", G_RESOURCE_LOOKUP_FLAGS_NONE , &error);
	if (css_buffer) {
		/* make sure that scale is ok */
		if (com.pref.gui.font_scale < 70.0) com.pref.gui.font_scale = 100;

		GString *string = g_string_new_len(NULL, 0);
		string = g_string_append_len(string, g_bytes_get_data(css_buffer, NULL), g_bytes_get_size(css_buffer));

		gchar *first_line = g_strdup_printf("* { font-size: %lfem; -gtk-icon-style: %s; }",
				1.0 + (com.pref.gui.font_scale - 100.0) / 1000.0, com.pref.gui.icon_symbolic ? "symbolic" : "regular");

		g_string_replace(string, "* { font-size: 1.0em; -gtk-icon-style: regular; }", first_line, 1);
		gchar *updated_css = g_string_free(string, FALSE);

		GtkCssProvider *css_provider = gtk_css_provider_new();
		GdkDisplay *display = gdk_display_get_default();
		GdkScreen *screen = gdk_display_get_default_screen(display);
		gtk_style_context_add_provider_for_screen(screen,
				GTK_STYLE_PROVIDER(css_provider),
				GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

		gtk_css_provider_load_from_data(css_provider, updated_css, -1, NULL);

		siril_debug_print("Successfully loaded /org/siril/ui/siril.css\n");
		g_free(first_line);
		g_bytes_unref(css_buffer);
		g_free(updated_css);
		g_object_unref(css_provider);
	} else {
		g_print(_("Failed to load CSS sheet: %s.\n"), error->message);
	}
}

