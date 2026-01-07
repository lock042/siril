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
#include "splashscreen.h"
#include "core/siril.h"

static GtkWidget *splash_window = NULL;
static GtkWidget *splash_progress = NULL;
static GtkWidget *splash_label = NULL;
static gboolean splash_is_active = FALSE;

/* Create and show the splash screen */
void show_splash_screen() {
	GtkWidget *overlay;
	GtkWidget *vbox;
	GtkWidget *image;
	GdkPixbuf *pixbuf;
	GError *error = NULL;

	splash_is_active = TRUE;

	/* Create a window without decorations */
	splash_window = gtk_window_new(GTK_WINDOW_POPUP);
	gtk_window_set_position(GTK_WINDOW(splash_window), GTK_WIN_POS_CENTER);
	gtk_window_set_resizable(GTK_WINDOW(splash_window), FALSE);
	gtk_window_set_keep_above(GTK_WINDOW(splash_window), TRUE);
	gtk_widget_set_app_paintable(splash_window, TRUE);

	/* Try to load the splash image from resources */
	pixbuf = gdk_pixbuf_new_from_resource("/org/siril/ui/pixmaps/splash.png", &error);
	if (!pixbuf) {
		/* If no splash image available, create a placeholder */
		pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, 600, 400);
		gdk_pixbuf_fill(pixbuf, 0x1a1a1aff); /* Dark gray background */

		if (error) {
			g_warning("Failed to load splash image: %s", error->message);
			g_clear_error(&error);
		}
	}

	/* Create the base image */
	image = gtk_image_new_from_pixbuf(pixbuf);

	/* Get image dimensions to size the window correctly */
	int img_width = gdk_pixbuf_get_width(pixbuf);
	g_object_unref(pixbuf);

	/* Create an overlay to put text and progress bar over the image */
	overlay = gtk_overlay_new();
	gtk_container_add(GTK_CONTAINER(splash_window), overlay);
	gtk_container_add(GTK_CONTAINER(overlay), image);

	/* Create a vbox for all overlay elements */
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_halign(vbox, GTK_ALIGN_FILL);
	gtk_widget_set_valign(vbox, GTK_ALIGN_FILL);
	gtk_overlay_add_overlay(GTK_OVERLAY(overlay), vbox);

	/* Add spacer to push title towards center/top */
	GtkWidget *top_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_vexpand(top_spacer, TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), top_spacer, TRUE, TRUE, 0);

	/* Add title and subtitle labels */
	GtkWidget *title_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 20);
	gtk_widget_set_halign(title_vbox, GTK_ALIGN_START);
	gtk_widget_set_margin_top(title_vbox, 20);
	gtk_widget_set_margin_start(title_vbox, 20);
	gtk_box_pack_start(GTK_BOX(vbox), title_vbox, FALSE, FALSE, 0);

	/* Title label with version */
	GtkWidget *title_label = gtk_label_new(NULL);
	gchar *package_str = g_strdup(PACKAGE_STRING);
	if (package_str && package_str[0]) {
		package_str[0] = g_ascii_toupper(package_str[0]);
	}
	gchar *title_markup = g_strdup_printf("<span size='30000' weight='bold' foreground='white'>%s</span>", package_str);
	gtk_label_set_markup(GTK_LABEL(title_label), title_markup);
	g_free(title_markup);
	g_free(package_str);
	gtk_label_set_xalign(GTK_LABEL(title_label), 0.0);
	gtk_box_pack_start(GTK_BOX(title_vbox), title_label, FALSE, FALSE, 0);

	/* Subtitle label */
	GtkWidget *subtitle_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(subtitle_label),
		"<span style='italic' size='large' foreground='white'>Astronomical Image Processing</span>");
	gtk_label_set_xalign(GTK_LABEL(subtitle_label), 0.0);
	gtk_box_pack_start(GTK_BOX(title_vbox), subtitle_label, FALSE, FALSE, 0);

	/* Add spacer to push progress bar to bottom */
	GtkWidget *middle_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_vexpand(middle_spacer, TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), middle_spacer, TRUE, TRUE, 0);

	/* Create a vertical box for progress elements at the bottom */
	GtkWidget *progress_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
	gtk_widget_set_margin_start(progress_vbox, 10);
	gtk_widget_set_margin_end(progress_vbox, 10);
	gtk_widget_set_margin_bottom(progress_vbox, 10);
	gtk_box_pack_start(GTK_BOX(vbox), progress_vbox, FALSE, FALSE, 0);

	/* Create the progress bar */
	splash_progress = gtk_progress_bar_new();
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(splash_progress), FALSE);
	gtk_widget_set_size_request(splash_progress, img_width - 20, 10);
	gtk_box_pack_start(GTK_BOX(progress_vbox), splash_progress, FALSE, FALSE, 0);

	/* Create the label for status messages */
	splash_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(splash_label), "<span foreground='white'>Starting Siril...</span>");
	gtk_label_set_xalign(GTK_LABEL(splash_label), 0.0);
	gtk_widget_set_margin_start(splash_label, 5);
	gtk_widget_set_margin_end(splash_label, 5);
	gtk_box_pack_start(GTK_BOX(progress_vbox), splash_label, FALSE, FALSE, 0);

	/* Apply CSS for styling - only to this window */
	GtkCssProvider *css_provider = gtk_css_provider_new();
	const gchar *css_data =
		"window.splash-screen { background-color: #1a1a1a; border: 1px solid #333333; }"
		"window.splash-screen progressbar { min-height: 10px; }"
		"window.splash-screen progressbar trough { background-color: rgba(42, 42, 42, 0.8); border-radius: 5px; }"
		"window.splash-screen progressbar progress { background-color: #4a9eff; border-radius: 5px; }"
		"window.splash-screen label { color: #ffffff; font-size: 11px; }";

	gtk_css_provider_load_from_data(css_provider, css_data, -1, NULL);

	/* Apply CSS only to this window's context */
	GtkStyleContext *context = gtk_widget_get_style_context(splash_window);
	gtk_style_context_add_provider(context,
		GTK_STYLE_PROVIDER(css_provider),
		GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
	gtk_style_context_add_class(context, "splash-screen");

	g_object_unref(css_provider);

	/* Show all widgets */
	gtk_widget_show_all(splash_window);

	/* Process events to make the window appear */
	while (gtk_events_pending())
		gtk_main_iteration();
}

/* Update the progress bar and message */
void update_splash_progress(const gchar *message, gdouble fraction) {
	if (!splash_window || !splash_progress || !splash_label)
		return;

	/* Clamp fraction between 0.0 and 1.0 */
	if (fraction < 0.0)
		fraction = 0.0;
	if (fraction > 1.0)
		fraction = 1.0;

	/* Update the progress bar */
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(splash_progress), fraction);

	/* Update the message if provided */
	if (message && *message) {
		gchar *markup = g_strdup_printf("<span foreground='white'>%s</span>", message);
		gtk_label_set_markup(GTK_LABEL(splash_label), markup);
		g_free(markup);
	}

	/* Force the splash screen to stay on top */
//	gtk_window_present(GTK_WINDOW(splash_window));

	/* Process events to update the display */
	while (gtk_events_pending())
		gtk_main_iteration();
}

/* Close and destroy the splash screen */
void close_splash_screen() {
	if (splash_window) {
		gtk_widget_destroy(splash_window);
		splash_window = NULL;
		splash_progress = NULL;
		splash_label = NULL;
		splash_is_active = FALSE;

		/* Process remaining events */
		while (gtk_events_pending())
			gtk_main_iteration();
	}
}

/* Check if splash screen is currently active */
gboolean is_splash_screen_active() {
	return splash_is_active;
}