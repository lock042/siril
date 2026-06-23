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

void show_splash_screen() {
	GtkWidget *overlay;
	GtkWidget *image;

	splash_is_active = TRUE;

	/* GTK4 dropped gtk_window_set_type_hint, gtk_window_set_position,
	 * gtk_widget_set_app_paintable, and the GTK_WINDOW_TOPLEVEL constructor
	 * argument.  The window manager places splash-style windows correctly
	 * when they are decorated=FALSE + transient_for unset. */
	splash_window = gtk_window_new();
	gtk_window_set_resizable(GTK_WINDOW(splash_window), FALSE);
	gtk_window_set_decorated(GTK_WINDOW(splash_window), FALSE);
	gtk_window_set_title(GTK_WINDOW(splash_window), "Siril Startup");

	/* Use GtkPicture rather than GtkImage: GtkImage is sized by icon-size
	 * and won't render the full-resolution PNG; GtkPicture takes the
	 * paintable's natural size when can-shrink is FALSE. */
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
		image = gtk_drawing_area_new();
		gtk_widget_set_size_request(image, img_width, img_height);
	}
	gtk_window_set_default_size(GTK_WINDOW(splash_window), img_width, img_height);

	overlay = gtk_overlay_new();
	gtk_window_set_child(GTK_WINDOW(splash_window), overlay);
	gtk_overlay_set_child(GTK_OVERLAY(overlay), image);

	/* Outer column anchored to the bottom: branding on top, progress strip below.
	 * Only the progress strip gets a gradient; the branding floats over the image
	 * with no backdrop (the image bottom is dark enough). */
	GtkWidget *bottom_col = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_halign(bottom_col, GTK_ALIGN_FILL);
	gtk_widget_set_valign(bottom_col, GTK_ALIGN_END);
	gtk_overlay_add_overlay(GTK_OVERLAY(overlay), bottom_col);

	/* Branding — no background, relies on the naturally dark lower portion of
	 * the image for contrast. */
	GtkWidget *branding_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_widget_set_margin_start(branding_box, 20);
	gtk_widget_set_margin_end(branding_box, 20);
	gtk_widget_set_margin_bottom(branding_box, 10);
	gtk_box_append(GTK_BOX(bottom_col), branding_box);

	GtkWidget *title_label = gtk_label_new(NULL);
	gchar *package_str = g_strdup(PACKAGE_STRING);
	if (package_str && package_str[0])
		package_str[0] = g_ascii_toupper(package_str[0]);
	gchar *title_markup = g_strdup_printf(
		"<span size='28pt' weight='bold' foreground='white'>%s</span>", package_str);
	gtk_label_set_markup(GTK_LABEL(title_label), title_markup);
	g_free(title_markup);
	g_free(package_str);
	gtk_label_set_xalign(GTK_LABEL(title_label), 0.0);
	gtk_box_append(GTK_BOX(branding_box), title_label);

	GtkWidget *subtitle_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(subtitle_label),
		"<span style='italic' size='12pt' foreground='#cccccc'>Astronomical Image Processing</span>");
	gtk_label_set_xalign(GTK_LABEL(subtitle_label), 0.0);
	gtk_box_append(GTK_BOX(branding_box), subtitle_label);

	/* Progress strip — gradient scoped only to this small band */
	GtkWidget *progress_strip = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_widget_set_margin_start(progress_strip, 20);
	gtk_widget_set_margin_end(progress_strip, 20);
	gtk_widget_set_margin_bottom(progress_strip, 14);
	gtk_widget_add_css_class(progress_strip, "splash-progress");
	gtk_box_append(GTK_BOX(bottom_col), progress_strip);

	splash_progress = gtk_progress_bar_new();
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(splash_progress), FALSE);
	gtk_widget_set_hexpand(splash_progress, TRUE);
	gtk_box_append(GTK_BOX(progress_strip), splash_progress);

	splash_label = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(splash_label),
		"<span foreground='#aaaaaa'>Starting Siril…</span>");
	gtk_label_set_xalign(GTK_LABEL(splash_label), 0.0);
	gtk_box_append(GTK_BOX(progress_strip), splash_label);

	/* Class selectors keep the rules scoped to .splash-screen even though
	 * the provider is registered display-wide (gtk_style_context_add_provider
	 * on a GtkStyleContext is deprecated in GTK 4.10+). */
	siril_register_css_for_display("siril-splash-screen",
		"window.splash-screen { background-color: #0d0d0d; }"
		".splash-progress {"
		"  padding: 6px 0 4px 0;"
		"  background: linear-gradient(to bottom,"
		"    rgba(0,0,0,0) 0%,"
		"    rgba(0,0,0,0.50) 50%,"
		"    rgba(0,0,0,0.65) 100%);"
		"}"
		"window.splash-screen progressbar { min-height: 5px; }"
		"window.splash-screen progressbar trough {"
		"  background-color: rgba(255,255,255,0.18);"
		"  border-radius: 3px;"
		"}"
		"window.splash-screen progressbar progress {"
		"  background-color: #4a9eff;"
		"  border-radius: 3px;"
		"}"
		"window.splash-screen label { font-size: 11pt; }");
	gtk_widget_add_css_class(splash_window, "splash-screen");

	gtk_widget_set_visible(splash_window, TRUE);

	while (g_main_context_pending(NULL))
		g_main_context_iteration(NULL, FALSE);
}

void update_splash_progress(const gchar *message, gdouble fraction) {
	if (!splash_window || !splash_progress || !splash_label)
		return;

	fraction = CLAMP(fraction, 0.0, 1.0);
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(splash_progress), fraction);

	if (message && *message) {
		gchar *markup = g_strdup_printf("<span foreground='#aaaaaa'>%s</span>", message);
		gtk_label_set_markup(GTK_LABEL(splash_label), markup);
		g_free(markup);
	}

	while (g_main_context_pending(NULL))
		g_main_context_iteration(NULL, FALSE);
}

void close_splash_screen() {
	if (splash_window) {
		gtk_window_destroy(GTK_WINDOW(splash_window));
		splash_window = NULL;
		splash_progress = NULL;
		splash_label = NULL;
		splash_is_active = FALSE;

		while (g_main_context_pending(NULL))
			g_main_context_iteration(NULL, FALSE);
	}
}

gboolean is_splash_screen_active() {
	return splash_is_active;
}
