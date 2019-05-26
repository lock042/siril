/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
 *
 *
 * This file contains functions related GUI initialization and mode change.
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include "core/siril.h"
#include "gui/callbacks.h"
#include "gui/gui.h"
#include "core/proto.h"
#include "gui/script_menu.h"
#include "gui/image_display.h"
#include "core/command.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "io/sequence.h"
#include "algos/star_finder.h"

#define GLADE_FILE "siril3.glade"
#define PLANETARY_GLADE_FILE "siril_planetary3.glade"

char *siril_sources[] = {
#ifdef _WIN32
	"../share/siril",
#elif (defined(__APPLE__) && defined(__MACH__))
	"/tmp/siril/Contents/Resources/share/siril/",
#endif
	PACKAGE_DATA_DIR"/",
	"/usr/share/siril/",
	"/usr/local/share/siril/",
	""
};

GtkWidget* lookup_widget(const gchar *widget_name) {
	if (!builder)
		fprintf(stderr, "builder is null, GUI not loaded\n");
	GObject *object = gtk_builder_get_object(builder, widget_name);
	if (!object)
		fprintf(stderr, "Could not find widget %s in builder\n", widget_name);
	return GTK_WIDGET(object);
}

static const gchar *checking_css_filename() {
	printf(_("Checking GTK version ... GTK-%d.%d\n"), GTK_MAJOR_VERSION, GTK_MINOR_VERSION);
	if ((GTK_MAJOR_VERSION >= 3) && (GTK_MINOR_VERSION >= 20))
		return "gtk.css";
	else if ((GTK_MAJOR_VERSION >= 3) && (GTK_MINOR_VERSION < 20))
		return "gtk_old.css";
	else return NULL;
}

/**
 * Loads the css sheet
 * @param path path of the file being loaded
 */
static void load_css_style_sheet (char *path) {
	GtkCssProvider *css_provider;
	GdkDisplay *display;
	GdkScreen *screen;
	gchar *CSSFile;
	const gchar *css_filename;

	css_filename = checking_css_filename();
	if (css_filename == NULL) {
		printf(_("The version of GTK does not match requirements: (GTK-%d.%d)\n"), GTK_MAJOR_VERSION, GTK_MINOR_VERSION);
		exit(1);
	}

	CSSFile = g_build_filename (path, css_filename, NULL);
	if (!g_file_test (CSSFile, G_FILE_TEST_EXISTS)) {
		g_error (_("Unable to load CSS style sheet file: %s. Please reinstall %s\n"), CSSFile, PACKAGE);
	}
	else {
		css_provider = gtk_css_provider_new();
		display = gdk_display_get_default();
		screen = gdk_display_get_default_screen(display);
		gtk_style_context_add_provider_for_screen(screen,
				GTK_STYLE_PROVIDER(css_provider),
				GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
		gtk_css_provider_load_from_path(css_provider, CSSFile, NULL);
		fprintf(stdout, _("Successfully loaded '%s'\n"), CSSFile);
		g_object_unref (css_provider);
	}
	g_free(CSSFile);
}


static void uninit_gui() {
	if (!builder) return;
	GSList *list = gtk_builder_get_objects(builder);
	while (list) {
		gtk_widget_destroy(list->data);
		list = g_slist_next(list);
	}
	g_slist_free(list);
	g_object_unref(builder);
	builder = NULL;
}

char *load_glade_file(char *start_cwd, enum _siril_mode mode) {
	int i = 0;
	const char *glade_file_name = mode == MODE_PLANETARY ? PLANETARY_GLADE_FILE : GLADE_FILE;

	/* try to load the glade file, from the sources defined above */
	builder = gtk_builder_new();

	do {
		GError *err = NULL;
		gchar *gladefile;

		gladefile = g_build_filename (siril_sources[i], glade_file_name, NULL);
		if (gtk_builder_add_from_file (builder, gladefile, &err)) {
			fprintf(stdout, _("Successfully loaded '%s'\n"), gladefile);
			g_free(gladefile);
			break;
		}
		fprintf (stderr, _("%s. Looking into another directory...\n"), err->message);
		g_error_free(err);
		g_free(gladefile);
		i++;
	} while (i < G_N_ELEMENTS(siril_sources));
	if (i == G_N_ELEMENTS(siril_sources)) {
		fprintf(stderr, _("%s was not found or contains errors, cannot render GUI. Exiting.\n"), glade_file_name);
		exit(EXIT_FAILURE);
	}
	/* get back to the saved working directory */
	return siril_sources[i];
}

void show_supported_files(gchar *supported_files) {
	GtkLabel *label_supported = GTK_LABEL(lookup_widget("label_supported_types"));
	gtk_label_set_text(label_supported, supported_files);
}

/* load new GUI mode: planetary or deep-sky.
 * releases the old GUI before loading the new mode.
 */
void init_gui(enum _siril_mode new_mode, char *start_cwd) {
	gchar *siril_path = NULL;
	int i = 0;

	if (new_mode == com.siril_mode)
		return;

	uninit_gui();	// release all widgets and the builder
	com.siril_mode = new_mode;

	if (new_mode == MODE_NO_GUI)
		return;

	if (new_mode == MODE_PLANETARY)
		fprintf(stdout, _("Initializing Siril Planetary interface\n"));
	else fprintf(stdout, _("Initializing Siril Deep-Sky interface\n"));

	load_prefered_theme(com.combo_theme);
	siril_path = load_glade_file(start_cwd, new_mode);
	load_css_style_sheet(siril_path);

	gtk_builder_connect_signals (builder, NULL);
	initialize_all_GUI();

	/* handling OS-X integration */
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	set_osx_integration(osx_app, siril_path);
	g_object_unref(osx_app);
#endif

	update_used_memory();
}

/* the trigger for mode change */
void on_menuitemplanetary_toggled(GtkCheckMenuItem *checkmenuitem,
		gpointer user_data) {
	/* this is not at all that simple.
	 * a master thread will have to manage the two sets of widgets, because
	 * with this technique, a gtk_main being disposed is calling his own
	 * destroying and re-creating a new builder. gtk_main has to stop
	 * naturally, and a new builder then can be loaded. */
	/*if (gtk_check_menu_item_get_active(checkmenuitem))
		init_gui(MODE_PLANETARY);
	else init_gui(MODE_DEEP_SKY);*/
	com.requested_mode =
		gtk_check_menu_item_get_active(checkmenuitem) ? MODE_PLANETARY : MODE_DEEP_SKY;
	gtk_main_quit();	// the main will handle the change, maybe
}

void load_prefered_theme(gint theme) {
	GtkSettings *settings;

	settings = gtk_settings_get_default();
	g_object_get(settings, "gtk-application-prefer-dark-theme", &com.have_dark_theme,
				NULL);

	if ((theme == 1) || (com.have_dark_theme && theme == 0)) {
		com.want_dark = TRUE;
	} else {
		com.want_dark = FALSE;
	}

	g_object_set(settings, "gtk-application-prefer-dark-theme", com.want_dark, NULL);
}

#ifdef MAC_INTEGRATION

static gboolean osx_open_file(GtkosxApplication *osx_app, gchar *path, gpointer data){
	if (path != NULL) {
		open_single_image(path);
		return FALSE;
	}
	return TRUE;
}

static void gui_add_osx_to_app_menu(GtkosxApplication *osx_app, const gchar *item_name, gint index) {
	GtkWidget *item;

	item = lookup_widget(item_name);
	if (GTK_IS_MENU_ITEM(item))
		gtkosx_application_insert_app_menu_item(osx_app, GTK_WIDGET(item), index);
}

static void set_osx_integration(GtkosxApplication *osx_app, gchar *siril_path) {
	GtkWidget *menubar = lookup_widget("menubar1");
	GtkWidget *file_quit_menu_item = lookup_widget("exit");
	GtkWidget *help_menu = lookup_widget("help1");
	GtkWidget *window_menu = lookup_widget("menuitemWindows");
	GtkWidget *sep;
	GdkPixbuf *icon;
	gchar *icon_path;
	
	g_signal_connect(osx_app, "NSApplicationOpenFile", G_CALLBACK(osx_open_file), NULL);

	gtk_widget_hide(menubar);

	gtkosx_application_set_menu_bar(osx_app, GTK_MENU_SHELL(menubar));
	gtkosx_application_set_window_menu(osx_app, GTK_MENU_ITEM(window_menu));

	gui_add_osx_to_app_menu(osx_app, "help_item1", 0);
	gui_add_osx_to_app_menu(osx_app, "help_get_scripts", 1);
	gui_add_osx_to_app_menu(osx_app, "help_update", 2);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 2);
	gui_add_osx_to_app_menu(osx_app, "settings", 3);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 4);

	gtk_widget_hide(file_quit_menu_item);
	gtk_widget_hide(help_menu);
	
	icon_path = g_build_filename(siril_path, "pixmaps/siril.svg", NULL);
	icon = gdk_pixbuf_new_from_file(icon_path, NULL);
	gtkosx_application_set_dock_icon_pixbuf(osx_app, icon);
		
	gtkosx_application_ready(osx_app);
	g_free(icon_path);
}

#endif

