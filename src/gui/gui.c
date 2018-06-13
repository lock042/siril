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

/* load new GUI mode: TRUE for planetary, FALSE for the deep-sky.
 * frees everything of the previous mode stored in com.planetary before loading
 * the new mode.
 */
void init_gui(gchar *startup_dir, enum _siril_mode new_mode) {
	gchar *siril_path = NULL;
	int i = 0;
	const char *glade_file_name = new_mode == PLANETARY ? PLANETARY_GLADE_FILE : GLADE_FILE;

	if (new_mode == com.siril_mode)
		return;

	uninit_gui();	// release all widgets and the builder
	com.siril_mode = new_mode;

	if (new_mode == NO_GUI)
		return;

	/* try to load the glade file, from the sources defined above */
	builder = gtk_builder_new();
	do {
		GError *err = NULL;
		gchar *gladefile;

		if (siril_sources[i][0] != '\0')
			gladefile = g_build_filename (siril_sources[i], glade_file_name, NULL);
		else gladefile = g_build_filename (startup_dir, glade_file_name, NULL);
		if (gtk_builder_add_from_file (builder, gladefile, &err)) {
			fprintf(stdout, _("Successfully loaded '%s'\n"), gladefile);
			g_free(gladefile);
			break;
		}
		fprintf(stderr, _("%s. Looking into another directory...\n"), err->message);
		g_error_free(err);
		g_free(gladefile);
		i++;
	} while (i < sizeof(siril_sources)/sizeof(char *));

	if (i == sizeof(siril_sources) / sizeof(char *)) {
		fprintf(stderr, _("%s was not found or contains errors, cannot render GUI. Exiting.\n"),
				glade_file_name);
		exit(EXIT_FAILURE);
	}
	siril_path = siril_sources[i];
	if (siril_path[0] == '\0')
		siril_path = startup_dir;

	gtk_builder_connect_signals(builder, NULL);

	com.vport[RED_VPORT] = lookup_widget("drawingarear");
	com.vport[GREEN_VPORT] = lookup_widget("drawingareag");
	com.vport[BLUE_VPORT] = lookup_widget("drawingareab");
	com.vport[RGB_VPORT] = lookup_widget("drawingareargb");
	if (com.siril_mode == DEEP_SKY) {
		com.preview_area[0] = lookup_widget("drawingarea_preview1");
		com.preview_area[1] = lookup_widget("drawingarea_preview2");
	}

	initialize_remap();
	initialize_scrollbars();
	init_mouse();

	/* Keybord Shortcuts */
	initialize_shortcuts();

	if (com.siril_mode == DEEP_SKY) {
		/* Select combo boxes that trigger some text display or other things */
		gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods")), 0);
		gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstacksel")), 0);

		/* initialize comboboxs of extraction background */
		GtkComboBox *order = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_polyorder"));
		gtk_combo_box_set_active(order, POLY_4);
		GtkComboBox *grad_inter = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combobox_gradient_inter"));
		gtk_combo_box_set_active(grad_inter, 0);

		/* Create tags associated with the buffer for the output text. */
		GtkTextBuffer *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(gtk_builder_get_object(builder, "output")));
		/* Tag with weight bold and tag name "bold" . */
		gtk_text_buffer_create_tag (tbuf, "bold", "weight", PANGO_WEIGHT_BOLD, NULL);
		/* Tag with style normal */
		gtk_text_buffer_create_tag (tbuf, "normal", "weight", PANGO_WEIGHT_NORMAL, NULL);
		/* Couleur Tags */
		gtk_text_buffer_create_tag (tbuf, "red", "foreground", "#e72828", NULL);
		gtk_text_buffer_create_tag (tbuf, "salmon", "foreground", "#ff9898", NULL);
		gtk_text_buffer_create_tag (tbuf, "green", "foreground", "#01b301", NULL);
		gtk_text_buffer_create_tag (tbuf, "blue", "foreground", "#7a7af8", NULL);
		gtk_text_buffer_create_tag (tbuf, "plum", "foreground", "#8e4585", NULL);
	}

	/* needs com.zoom_value */
	zoomcombo_update_display_for_zoom();

	/* needs com.seq init */
	adjust_sellabel();

	/* load the css sheet for general style */
	load_css_style_sheet (siril_path);

	/* initialize menu gui */
	update_MenuItem();
	if (com.siril_mode == DEEP_SKY) {
		initialize_script_menu();

		/* initialize command completion */
		init_completion_command();

		/* initialize registration methods */
		initialize_registration_methods();

		/* initialize stacking methods */
		initialize_stacking_methods();

		// This crashes for an unknown reason
		//init_peaker_GUI();
	}

	/* initialize preprocessing */
	initialize_preprocessing();

	/* register some callbacks */
	register_selection_update_callback(update_export_crop_label);

	/* initialization of the binning parameters */
	GtkComboBox *binning = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combobinning"));
	gtk_combo_box_set_active(binning, 0);

	/* initialization of some paths */
	initialize_path_directory();

	/* initialization of default FITS extension */
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("combobox_ext"));
	gtk_combo_box_set_active_id(box, com.ext);
	initialize_FITS_name_entries();

	set_GUI_CWD();	// to call after readinitfile
	set_GUI_misc();	// to call after readinitfile

#ifdef HAVE_LIBRAW
	set_GUI_LIBRAW();
#endif
	set_GUI_photometry();

	update_spinCPU(com.max_thread);

	if (com.have_dark_theme) {
		if (com.siril_mode == DEEP_SKY) {
			/* Put dark icons */
			printf("Loading dark theme...\n");
			gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("rotate90_anticlock_button")), lookup_widget("rotate90-acw_dark"));
			gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("rotate90_clock_button")), lookup_widget("rotate90-cw_dark"));
			gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("mirrorx_button")), lookup_widget("image_mirrorx_dark"));
			gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("mirrory_button")), lookup_widget("image_mirrory_dark"));
			gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("seqlist_button")), lookup_widget("image_seqlist_dark"));
		}
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget("histogram_button")), lookup_widget("image_histogram_dark"));
	}

	update_used_memory();
}
