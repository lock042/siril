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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "algos/comparison_stars.h"
#include "nina_light_curve.h"
#include "io/sequence.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/plot.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *file_chooser = NULL;
static GtkWidget *use_comp1 = NULL;
static GtkWidget *use_comp2 = NULL;
static GtkWidget *display_curve = NULL;
static GtkWidget *sep = NULL;
static GtkWidget *apert = NULL;
static GtkWidget *apert_value = NULL;
static GtkWidget *inner_value = NULL;
static GtkWidget *outer_value = NULL;
static double radius_value = 0.;
const char *radius_label = NULL;

static void on_nina_lc_response(GtkDialog* self, gint response_id, gpointer user_data);

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Automated light curve"), NULL,
			0, _("_Close"), GTK_RESPONSE_REJECT, _("_OK"), GTK_RESPONSE_ACCEPT, NULL);
	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(siril_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_nina_lc_response), NULL);

	file_chooser = gtk_file_chooser_button_new (_("Select the comparison star list file"),
			GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(file_chooser), com.wd);
	g_object_set(G_OBJECT(file_chooser), "margin", 15, NULL);
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("CSV file (*.csv)"));
	gtk_file_filter_add_pattern(f, "*.csv");
	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(file_chooser), f);

	// Sets the "OK" button as default one
	GtkWidget *OK_button = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
	gtk_widget_set_can_default(OK_button, TRUE);
	gtk_widget_grab_default(OK_button);
	gtk_widget_grab_focus(OK_button);
	gtk_style_context_add_class(gtk_widget_get_style_context(OK_button), "suggested-action");

	// Definition of all the graphical items
	// Top label
	GtkWidget *win_label = gtk_label_new(_("Process a sequence to get a light curve on a star using the list\nof reference stars created by Siril or the NINA exoplanet plugin"));
	gtk_label_set_line_wrap(GTK_LABEL(win_label), TRUE);
	g_object_set(G_OBJECT(win_label), "margin", 15, NULL);

	/* Defines the options box */
	GtkWidget *options_box;
	options_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(options_box), TRUE);
	// Comp1 radio
	use_comp1 = gtk_check_button_new_with_label(_("Use comparative stars selected for their color"));
	gtk_widget_set_tooltip_text(use_comp1, _("Color similar to the target mean they will get extincted the same way by the changing atmosphere"));
	g_object_set(G_OBJECT(use_comp1), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(use_comp1), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(use_comp1), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(use_comp1), "margin-bottom", 0, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(use_comp1), TRUE);
	gtk_container_add(GTK_CONTAINER(options_box), use_comp1);
	// Comp2 radio
	use_comp2 = gtk_check_button_new_with_label(_("Use comparative stars selected by the AAVSO"));
	gtk_widget_set_tooltip_text(use_comp2, _("The AAVSO gives stars that are known to not be variable"));
	g_object_set(G_OBJECT(use_comp2), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(use_comp2), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(use_comp2), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(use_comp2), "margin-bottom", 0, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(use_comp2), TRUE);
	gtk_container_add(GTK_CONTAINER(options_box), use_comp2);
	// Display radio
	display_curve = gtk_check_button_new_with_label(_("Display the light curve"));
	gtk_widget_set_tooltip_text(display_curve, _("if not checked, a PNG image of the graph will be generated instead"));
	g_object_set(G_OBJECT(display_curve), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(display_curve), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(display_curve), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(display_curve), "margin-bottom", 0, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(display_curve), TRUE);
	gtk_container_add(GTK_CONTAINER(options_box), display_curve);

	// Horizontal Separator
	sep = gtk_separator_new (GTK_ORIENTATION_HORIZONTAL);
	g_object_set(G_OBJECT(sep), "margin-top", 5, NULL);
	g_object_set(G_OBJECT(sep), "margin-bottom", 5, NULL);

	// Aperture Label
	GtkWidget *aperture_label = gtk_label_new(_("The current aperture parameters are as follow:\n(they can be managed in the Preferences/Photometry tab)"));
	gtk_label_set_line_wrap(GTK_LABEL(aperture_label), TRUE);
	g_object_set(G_OBJECT(aperture_label), "margin", 15, NULL);
	g_object_set(G_OBJECT(aperture_label), "margin-top", 0, NULL);
	gtk_widget_set_halign(aperture_label, GTK_ALIGN_START);

	// Aperture parameters
	/* Defines the 3 boxes */
	GtkWidget *aperture_box, *inner_box, *outer_box;
	aperture_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(aperture_box), TRUE);
	inner_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(inner_box), TRUE);
	outer_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(outer_box), TRUE);

	/* Fill the 1st box */
	/* The label */
	radius_label = !com.pref.phot_set.force_radius ? "Radius/half-FWHM ratio:" : "Aperture radius (px):";
	apert = gtk_label_new(radius_label);	// The label depends on the "force_radius" value
	gtk_label_set_line_wrap(GTK_LABEL(apert), TRUE);
	g_object_set(G_OBJECT(apert), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(apert), "margin-top", 0, NULL);	
	gtk_widget_set_halign(apert, GTK_ALIGN_START);
	gtk_container_add(GTK_CONTAINER(aperture_box), apert);
	/* The entry */
	apert_value = gtk_entry_new();
	gtk_widget_set_sensitive(apert_value, FALSE);
	g_object_set(G_OBJECT(apert_value), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(apert_value), "margin-top", 0, NULL);
	gtk_widget_set_halign(apert_value, GTK_ALIGN_START);
	radius_value = !com.pref.phot_set.force_radius ? com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
	gtk_entry_set_text(GTK_ENTRY(apert_value), g_strdup_printf("%1.2lf", radius_value));
	gtk_container_add(GTK_CONTAINER(aperture_box), apert_value);

	/* Fill the 2nd box */
	/* The label */
	GtkWidget *inn = gtk_label_new(_("Inner radius of the annulus (px):"));
	gtk_label_set_line_wrap(GTK_LABEL(inn), TRUE);
	g_object_set(G_OBJECT(inn), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(inn), "margin-top", 0, NULL);
	gtk_widget_set_halign(inn, GTK_ALIGN_START);
	gtk_container_add(GTK_CONTAINER(inner_box), inn);
	/* The entry */
	inner_value = gtk_entry_new();
	gtk_widget_set_sensitive(inner_value, FALSE);
	g_object_set(G_OBJECT(inner_value), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(inner_value), "margin-top", 0, NULL);
	gtk_widget_set_halign(inner_value, GTK_ALIGN_START);
	gtk_entry_set_text(GTK_ENTRY(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
	gtk_container_add(GTK_CONTAINER(inner_box), inner_value);

	/* Fill the 3rd box */
	/* The label */
	GtkWidget *oute = gtk_label_new(_("Outer radius of the annulus (px):"));
	gtk_label_set_line_wrap(GTK_LABEL(oute), TRUE);
	g_object_set(G_OBJECT(oute), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(oute), "margin-top", 0, NULL);
	gtk_widget_set_halign(oute, GTK_ALIGN_START);
	gtk_container_add(GTK_CONTAINER(outer_box), oute);
	/* The entry */
	outer_value = gtk_entry_new();
	gtk_widget_set_sensitive(outer_value, FALSE);
	g_object_set(G_OBJECT(outer_value), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(outer_value), "margin-top", 0, NULL);
	gtk_widget_set_halign(outer_value, GTK_ALIGN_START);
	gtk_entry_set_text(GTK_ENTRY(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));
	gtk_container_add(GTK_CONTAINER(outer_box), outer_value);

	// Gather the graphical elements //
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 20);
	gtk_container_add(GTK_CONTAINER(content_area), win_label);
	gtk_container_add(GTK_CONTAINER(content_area), file_chooser);
	gtk_container_add(GTK_CONTAINER(content_area), options_box);
	gtk_container_add(GTK_CONTAINER(content_area), sep);
	gtk_container_add(GTK_CONTAINER(content_area), aperture_label);
	gtk_container_add(GTK_CONTAINER(content_area), aperture_box);
	gtk_container_add(GTK_CONTAINER(content_area), inner_box);
	gtk_container_add(GTK_CONTAINER(content_area), outer_box);
	gtk_widget_show_all(GTK_WIDGET(content_area));
}

// the public getter
GtkWidget *get_nina_lc_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_nina_lc_response(GtkDialog* self, gint response_id, gpointer user_data) {
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		radius_label = !com.pref.phot_set.force_radius ? "Radius/half-FWHM ratio:" : "Aperture radius (px):";
		gtk_label_set_text(GTK_LABEL(apert), radius_label);
		radius_value = !com.pref.phot_set.force_radius ? com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
		gtk_entry_set_text(GTK_ENTRY(apert_value), g_strdup_printf("%1.2lf", radius_value));
		gtk_entry_set_text(GTK_ENTRY(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
		gtk_entry_set_text(GTK_ENTRY(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));
		gtk_widget_hide(dialog);
		return;
	}
	if (!sequence_is_loaded()) {	// Tests if a valid sequence is loaded
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No sequence loaded"));
		return;
	}
	gchar *nina_file = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
	if (!nina_file){	// Any file for processing?
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Choose a comparison stars file"));
		g_free(nina_file);
		return;
	}
	gchar *dirname = g_path_get_dirname(nina_file);
	if (g_strcmp0(dirname, com.wd) != 0) {	// Tests if the file is in the CWD
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The current comparison stars file is not located in the CWD"));
		gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(file_chooser), com.wd);
		gtk_file_chooser_unselect_all (GTK_FILE_CHOOSER(file_chooser));
		g_free(dirname);
		g_free(nina_file);
		return;
	}
	g_free(dirname);
	if (!has_wcs(&gfit)) {	// Is the current image properly plate solved
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The currently loaded image must be plate solved"));
		return;
	}
	purge_user_catalogue(CAT_AN_USER_TEMP);
	int layer = -1;
	if (com.seq.regparam) {
		for (int i = 0; i < com.seq.nb_layers; i++) {
			if (com.seq.regparam[i]) {
				layer = i;
				break;
			}
		}
	}
	if (layer == -1) {
		siril_debug_print("unregistered sequence\n");
		//siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The sequence must be registered"));
		if (com.seq.nb_layers == 3)
			layer = 1;
		else layer = 0;
	}

	if (sequence_drifts(&com.seq, layer, com.seq.rx / 4)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), _("The sequence seems to have a heavy drift, the computation of a light curve may not be accurate or possible"));
		// TODO: if clicked on close, do not continue
	}

	siril_log_message(_("Using preconfigured inner and outer photometry ring radii of %.1f and %.1f\n"),
			com.pref.phot_set.inner, com.pref.phot_set.outer);

	clear_all_photometry_and_plot();
	init_plot_colors();

	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	gboolean use_c1 = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(use_comp1));
	gboolean use_c2 = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(use_comp2));

	if (parse_nina_stars_file_using_WCS(args, nina_file, use_c1, use_c2, &gfit)) {
		// fail
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Something went wrong while saving plot"));
		free(args);
		g_free(nina_file);
		return;
	}
	g_free(nina_file);
	args->seq = &com.seq;
	args->layer = layer;
	args->display_graph = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(display_curve));
	siril_debug_print("starting PSF analysis of %d stars\n", args->nb);

	if (!start_in_new_thread(light_curve_worker, args)) {
		free(args);
		g_free(nina_file);
		return;
	}
	// Update of the UI
	radius_label = !com.pref.phot_set.force_radius ? "Radius/half-FWHM ratio:" : "Aperture radius (px):";
	gtk_label_set_text(GTK_LABEL(apert), radius_label);
	radius_value = !com.pref.phot_set.force_radius ? com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
	gtk_entry_set_text(GTK_ENTRY(apert_value), g_strdup_printf("%1.2lf", radius_value));
	gtk_entry_set_text(GTK_ENTRY(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
	gtk_entry_set_text(GTK_ENTRY(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));

}
