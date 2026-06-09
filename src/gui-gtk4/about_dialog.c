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

#include "core/siril.h"
#include "core/proto.h"
#include "gui-gtk4/utils.h"
#include "git-version.h"

#include "about_dialog.h"

static gchar **authors = (gchar *[] ) {
				"Vincent Hourdin <vh@free-astro.vinvin.tf>",
				"Cyril Richard <cyril@free-astro.org>",
				"Cécile Melis <cissou8@gmail.com>",
				"Adrian Knagg-Baugh <aje.baugh@gmail.com>",
				"François Meyer",
				NULL
};

static gchar **documenters = (gchar *[] ) {
				"Laurent Rogé <l.roge@siril.org> (2016-2022)",
				"Team free-astro (2022-"SIRIL_GIT_LAST_COMMIT_YEAR")",
				NULL
};

static gchar **artists = (gchar *[] ) {
				"Maxime Oudoux <max.oudoux@gmail.com>",
				"Tobias Bernard <tbernard@gnome.org>",
				"Cyril Richard <cyril@free-astro.org>",
				NULL
};

void siril_show_about_dialog() {
	GtkWindow *parent;
	gchar *copyright;
	gchar *version;
#ifdef SIRIL_UNSTABLE
	version = g_strdup_printf(_("%s\nThis is an unstable development release\n"
					"commit %s\nGTK %d.%d.%d\n"), VERSION, SIRIL_GIT_VERSION_ABBREV,
			gtk_get_major_version(), gtk_get_minor_version(), gtk_get_micro_version());
#else
	version = g_strdup_printf("%s\ncommit %s\nGTK %d.%d.%d", VERSION, SIRIL_GIT_VERSION_ABBREV,
			gtk_get_major_version(), gtk_get_minor_version(), gtk_get_micro_version());
#endif
	copyright = g_strdup_printf("Copyright © 2004-2011 François Meyer\n"
			"Copyright © 2012-%s Team free-astro", SIRIL_GIT_LAST_COMMIT_YEAR);

	parent = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	/* "logo-icon-name" resolves through the icon theme: main.c registers
	 * /org/siril/ui/pixmaps as a resource path so "siril" picks up
	 * siril.svg via GTK4's native paintable infrastructure, avoiding the
	 * deprecated "logo" pixbuf path. */

	GtkAboutDialog *dialog = GTK_ABOUT_DIALOG(gtk_about_dialog_new());
	g_object_set(G_OBJECT(dialog),
			"program-name",      PACKAGE,
			"title",             _("About Siril"),
			"logo-icon-name",    "siril",
			"version",           version,
			"copyright",         copyright,
			"authors",           authors,
			"documenters",       documenters,
			"artists",           artists,
			"comments",          _("Astronomical image (pre-)processing program"),
			"translator-credits", _("translator-credits"),
			"website",           PACKAGE_URL,
			"website-label",     _("Visit the Siril website"),
			"license-type",      GTK_LICENSE_GPL_3_0,
			NULL);
	gtk_window_set_transient_for(GTK_WINDOW(dialog), parent);

#if defined(OS_OSX) && GTK_CHECK_VERSION(4, 18, 0)
	/* Enable macOS traffic-light window controls on the About dialog header bar. */
	GtkWidget *titlebar = gtk_window_get_titlebar(GTK_WINDOW(dialog));
	if (GTK_IS_HEADER_BAR(titlebar))
		gtk_header_bar_set_use_native_controls(GTK_HEADER_BAR(titlebar), TRUE);
#endif

	gtk_window_set_hide_on_close(GTK_WINDOW(dialog), FALSE);
	gtk_window_present(GTK_WINDOW(dialog));

	g_free(copyright);
	g_free(version);
}
