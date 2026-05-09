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
#include "gui-gtk4/utils.h"

static GtkWidget *splash_window = NULL;
static GtkWidget *splash_progress = NULL;
static GtkWidget *splash_label = NULL;
static gboolean splash_is_active = FALSE;

/* Create and show the splash screen */
void show_splash_screen() {
	GtkWidget *overlay;
	GtkWidget *vbox;
	GtkWidget *image;
	GError *error = NULL;

	splash_is_active = TRUE;

	/* Phase 14: GTK4 dropped gtk_window_set_type_hint, gtk_window_set_position,
	 * gtk_widget_set_app_paintable, and the GTK_WINDOW_TOPLEVEL constructor
	 * argument.  The window manager places splash-style windows correctly
	 * when they are decorated=FALSE + transient_for unset. */
	splash_window = gtk_window_new();
	gtk_window_set_resizable(GTK_WINDOW(splash_window), FALSE);
	gtk_window_set_decorated(GTK_WINDOW(splash_window), FALSE);
	gtk_window_set_title(GTK_WINDOW(splash_window), "Siril Startup");

	/* Use GtkPicture rather than GtkImage for the splash background.
	 * GtkImage is sized by icon-size (small) and won't render the
	 * full-resolution PNG; GtkPicture takes the paintable's natural
	 * size when can-shrink is FALSE. */
	GdkTexture *splash_tex = gdk_texture_new_from_resource("/org/siril/ui/pixmaps/splash.png");
	int img_width = 600;
	int img_height = 400;
	if (splash_tex) {
		image = gtk_picture_new_for_paintable(GDK_PAINTABLE(splash_tex));
		gtk_picture_set_can_shrink(GTK_PICTURE(image), FALSE);
		gtk_picture_set_content_fit(GTK_PICTURE(image), GTK_CONTENT_FIT_FILL);
		img_width = gdk_texture_get_width(splash_tex);
		img_height = gdk_texture_get_height(splash_tex);
		g_object_unref(splash_tex);
	} else {
		/* Placeholder dark window if the splash resource is missing. */
		if (error) {
			g_warning("Failed to load splash image: %s", error->message);
			g_clear_error(&error);
		}
		image = gtk_drawing_area_new();
		gtk_widget_set_size_request(image, img_width, img_height);
	}
	gtk_window_set_default_size(GTK_WINDOW(splash_window), img_width, img_height);

	/* Create an overlay to put text and progress bar over the image */
	overlay = gtk_overlay_new();
	gtk_window_set_child(GTK_WINDOW(splash_window), overlay);
	gtk_overlay_set_child(GTK_OVERLAY(overlay), image);

	/* Create a vbox for all overlay elements */
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_halign(vbox, GTK_ALIGN_FILL);
	gtk_widget_set_valign(vbox, GTK_ALIGN_FILL);
	gtk_overlay_add_overlay(GTK_OVERLAY(overlay), vbox);

	/* Add spacer to push title towards center/top */
	GtkWidget *top_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_vexpand(top_spacer, TRUE);
	gtk_widget_set_hexpand(top_spacer, TRUE);
	gtk_widget_set_vexpand(top_spacer, TRUE);
	gtk_box_append(GTK_BOX(vbox), top_spacer);

	/* Add title and subtitle labels */
	GtkWidget *title_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 20);
	gtk_widget_set_halign(title_vbox, GTK_ALIGN_START);
	gtk_widget_set_margin_top(title_vbox, 20);
	gtk_widget_set_margin_start(title_vbox, 20);
	gtk_widget_set_halign(title_vbox, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(title_vbox, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(vbox), title_vbox);

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
	gtk_widget_set_halign(title_label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(title_label, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(title_vbox), title_label);

	/* Subtitle label */
	GtkWidget *subtitle_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(subtitle_label),
		"<span style='italic' size='large' foreground='white'>Astronomical Image Processing</span>");
	gtk_label_set_xalign(GTK_LABEL(subtitle_label), 0.0);
	gtk_widget_set_halign(subtitle_label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(subtitle_label, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(title_vbox), subtitle_label);

	/* Add spacer to push progress bar to bottom */
	GtkWidget *middle_spacer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_vexpand(middle_spacer, TRUE);
	gtk_widget_set_hexpand(middle_spacer, TRUE);
	gtk_widget_set_vexpand(middle_spacer, TRUE);
	gtk_box_append(GTK_BOX(vbox), middle_spacer);

	/* Create a vertical box for progress elements at the bottom */
	GtkWidget *progress_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
	gtk_widget_set_margin_start(progress_vbox, 10);
	gtk_widget_set_margin_end(progress_vbox, 10);
	gtk_widget_set_margin_bottom(progress_vbox, 10);
	gtk_widget_set_halign(progress_vbox, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(progress_vbox, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(vbox), progress_vbox);

	/* Create the progress bar */
	splash_progress = gtk_progress_bar_new();
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(splash_progress), FALSE);
	gtk_widget_set_size_request(splash_progress, img_width - 20, 10);
	gtk_widget_set_halign(splash_progress, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(splash_progress, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(progress_vbox), splash_progress);

	/* Create the label for status messages */
	splash_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(splash_label), "<span foreground='white'>Starting Siril...</span>");
	gtk_label_set_xalign(GTK_LABEL(splash_label), 0.0);
	gtk_widget_set_margin_start(splash_label, 5);
	gtk_widget_set_margin_end(splash_label, 5);
	gtk_widget_set_halign(splash_label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(splash_label, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(progress_vbox), splash_label);

	/* Apply CSS for styling — class selectors keep the rules scoped to
	 * .splash-screen even though the provider is registered display-wide
	 * (gtk_style_context_add_provider on a GtkStyleContext is deprecated
	 * in GTK 4.10+). */
	siril_register_css_for_display("siril-splash-screen",
		"window.splash-screen { background-color: #1a1a1a; border: 1px solid #333333; }"
		"window.splash-screen progressbar { min-height: 10px; }"
		"window.splash-screen progressbar trough { background-color: rgba(42, 42, 42, 0.8); border-radius: 5px; }"
		"window.splash-screen progressbar progress { background-color: #4a9eff; border-radius: 5px; }"
		"window.splash-screen label { color: #ffffff; font-size: 11px; }");
	gtk_widget_add_css_class(splash_window, "splash-screen");

	/* Show all widgets */
	gtk_widget_set_visible(splash_window, TRUE);

	/* Process events to make the window appear */
	while (g_main_context_pending(NULL))
		g_main_context_iteration(NULL, FALSE);
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

	/* Process events to update the display */
	while (g_main_context_pending(NULL))
		g_main_context_iteration(NULL, FALSE);
}

/* Close and destroy the splash screen */
void close_splash_screen() {
	if (splash_window) {
		gtk_window_destroy(GTK_WINDOW(splash_window));
		splash_window = NULL;
		splash_progress = NULL;
		splash_label = NULL;
		splash_is_active = FALSE;

		/* Process remaining events */
		while (g_main_context_pending(NULL))
			g_main_context_iteration(NULL, FALSE);
	}
}

/* Check if splash screen is currently active */
gboolean is_splash_screen_active() {
	return splash_is_active;
}
