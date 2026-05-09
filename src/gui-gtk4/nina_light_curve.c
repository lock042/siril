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
#include "core/siril_log.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "algos/comparison_stars.h"
#include "nina_light_curve.h"
#include "io/sequence.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/plot.h"
#include "gui-gtk4/utils.h"

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

static void on_nina_lc_response(GtkWindow *self, gint response_id, gpointer user_data);

static GtkWidget *nina_ok_button = NULL;

/* Phase 14G.4: GtkDialog → GtkWindow.  Per-button "clicked" → response shape. */
static void nina_btn_ok_clicked(GtkButton *btn, gpointer ud)    { (void)btn; on_nina_lc_response(GTK_WINDOW(dialog), GTK_RESPONSE_ACCEPT, ud); }
static void nina_btn_close_clicked(GtkButton *btn, gpointer ud) { (void)btn; on_nina_lc_response(GTK_WINDOW(dialog), GTK_RESPONSE_REJECT, ud); }

static void build_the_dialog() {
	/* Phase 14G.4: GtkDialog → GtkWindow. */
	dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Automated light curve"));
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_hide_on_close(GTK_WINDOW(dialog), TRUE);
	g_signal_connect(G_OBJECT(dialog), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	/* Phase 14: GtkFileChooserButton removed in GTK4.  Use a placeholder
	 * GtkButton until Phase 18 routes this through GtkFileDialog. */
	file_chooser = gtk_button_new_with_label(_("Select the comparison star list file"));
	gtk_widget_set_margin_start(GTK_WIDGET(file_chooser), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(file_chooser), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(file_chooser), 15);
	gtk_widget_set_margin_bottom(GTK_WIDGET(file_chooser), 15);
	// Definition of all the graphical items
	// Top label
	GtkWidget *win_label = gtk_label_new(_("Process a sequence to get a light curve on a star using the list\nof reference stars created by Siril or the NINA exoplanet plugin"));
	gtk_label_set_wrap(GTK_LABEL(win_label), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(win_label), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(win_label), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(win_label), 15);
	gtk_widget_set_margin_bottom(GTK_WIDGET(win_label), 15);
	/* Defines the options box */
	GtkWidget *options_box;
	options_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(options_box), TRUE);
	// Comp1 radio
	use_comp1 = gtk_check_button_new_with_label(_("Use comparative stars selected for their color"));
	gtk_widget_set_tooltip_text(use_comp1, _("Color similar to the target mean they will get extincted the same way by the changing atmosphere"));
	gtk_widget_set_margin_start(GTK_WIDGET(use_comp1), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(use_comp1), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(use_comp1), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(use_comp1), 0);
	siril_toggle_set_active(GTK_WIDGET(use_comp1), TRUE);
	gtk_box_append(GTK_BOX(options_box), use_comp1);
	// Comp2 radio
	use_comp2 = gtk_check_button_new_with_label(_("Use comparative stars selected by the AAVSO"));
	gtk_widget_set_tooltip_text(use_comp2, _("The AAVSO gives stars that are known to not be variable"));
	gtk_widget_set_margin_start(GTK_WIDGET(use_comp2), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(use_comp2), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(use_comp2), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(use_comp2), 0);
	siril_toggle_set_active(GTK_WIDGET(use_comp2), TRUE);
	gtk_box_append(GTK_BOX(options_box), use_comp2);
	// Display radio
	display_curve = gtk_check_button_new_with_label(_("Display the light curve"));
	gtk_widget_set_tooltip_text(display_curve, _("if not checked, a PNG image of the graph will be generated instead"));
	gtk_widget_set_margin_start(GTK_WIDGET(display_curve), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(display_curve), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(display_curve), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(display_curve), 0);
	siril_toggle_set_active(GTK_WIDGET(display_curve), TRUE);
	gtk_box_append(GTK_BOX(options_box), display_curve);

	// Horizontal Separator
	sep = gtk_separator_new (GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_margin_top(GTK_WIDGET(sep), 5);
	gtk_widget_set_margin_bottom(GTK_WIDGET(sep), 5);
	// Aperture Label
	GtkWidget *aperture_label = gtk_label_new(_("The current aperture parameters are as follow:\n(they can be managed in the Preferences/Photometry tab)"));
	gtk_label_set_wrap(GTK_LABEL(aperture_label), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(aperture_label), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(aperture_label), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(aperture_label), 15);
	gtk_widget_set_margin_bottom(GTK_WIDGET(aperture_label), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(aperture_label), 0);
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
	gtk_label_set_wrap(GTK_LABEL(apert), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(apert), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(apert), 0);
	gtk_widget_set_halign(apert, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(aperture_box), apert);
	/* The entry */
	apert_value = gtk_entry_new();
	gtk_widget_set_sensitive(apert_value, FALSE);
	gtk_widget_set_margin_start(GTK_WIDGET(apert_value), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(apert_value), 0);
	gtk_widget_set_halign(apert_value, GTK_ALIGN_START);
	radius_value = !com.pref.phot_set.force_radius ? com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
	gtk_editable_set_text(GTK_EDITABLE(apert_value), g_strdup_printf("%1.2lf", radius_value));
	gtk_box_append(GTK_BOX(aperture_box), apert_value);

	/* Fill the 2nd box */
	/* The label */
	GtkWidget *inn = gtk_label_new(_("Inner radius of the annulus (px):"));
	gtk_label_set_wrap(GTK_LABEL(inn), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(inn), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(inn), 0);
	gtk_widget_set_halign(inn, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(inner_box), inn);
	/* The entry */
	inner_value = gtk_entry_new();
	gtk_widget_set_sensitive(inner_value, FALSE);
	gtk_widget_set_margin_start(GTK_WIDGET(inner_value), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(inner_value), 0);
	gtk_widget_set_halign(inner_value, GTK_ALIGN_START);
	gtk_editable_set_text(GTK_EDITABLE(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
	gtk_box_append(GTK_BOX(inner_box), inner_value);

	/* Fill the 3rd box */
	/* The label */
	GtkWidget *oute = gtk_label_new(_("Outer radius of the annulus (px):"));
	gtk_label_set_wrap(GTK_LABEL(oute), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(oute), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(oute), 0);
	gtk_widget_set_halign(oute, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(outer_box), oute);
	/* The entry */
	outer_value = gtk_entry_new();
	gtk_widget_set_sensitive(outer_value, FALSE);
	gtk_widget_set_margin_start(GTK_WIDGET(outer_value), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(outer_value), 0);
	gtk_widget_set_halign(outer_value, GTK_ALIGN_START);
	gtk_editable_set_text(GTK_EDITABLE(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));
	gtk_box_append(GTK_BOX(outer_box), outer_value);

	// Gather the graphical elements //
	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 20);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_box_append(GTK_BOX(content_area), win_label);
	gtk_box_append(GTK_BOX(content_area), file_chooser);
	gtk_box_append(GTK_BOX(content_area), options_box);
	gtk_box_append(GTK_BOX(content_area), sep);
	gtk_box_append(GTK_BOX(content_area), aperture_label);
	gtk_box_append(GTK_BOX(content_area), aperture_box);
	gtk_box_append(GTK_BOX(content_area), inner_box);
	gtk_box_append(GTK_BOX(content_area), outer_box);

	/* Action area */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *btn_close = gtk_button_new_with_mnemonic(_("_Close"));
	g_signal_connect(btn_close, "clicked", G_CALLBACK(nina_btn_close_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), btn_close);
	nina_ok_button = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(nina_ok_button, "suggested-action");
	g_signal_connect(nina_ok_button, "clicked", G_CALLBACK(nina_btn_ok_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), nina_ok_button);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), nina_ok_button);

	gtk_window_set_child(GTK_WINDOW(dialog), content_area);
}

// the public getter
GtkWidget *get_nina_lc_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_nina_lc_response(GtkWindow *self, gint response_id, gpointer user_data) {
	(void)self; (void)user_data;
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		radius_label = !com.pref.phot_set.force_radius ? "Radius/half-FWHM ratio:" : "Aperture radius (px):";
		gtk_label_set_text(GTK_LABEL(apert), radius_label);
		radius_value = !com.pref.phot_set.force_radius ? com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
		gtk_editable_set_text(GTK_EDITABLE(apert_value), g_strdup_printf("%1.2lf", radius_value));
		gtk_editable_set_text(GTK_EDITABLE(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
		gtk_editable_set_text(GTK_EDITABLE(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));
		gtk_widget_set_visible(dialog, FALSE);
		return;
	}
	if (!sequence_is_loaded()) {	// Tests if a valid sequence is loaded
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No sequence loaded"));
		return;
	}
	gchar *nina_file = siril_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
	if (!nina_file){	// Any file for processing?
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Choose a comparison stars file"));
		g_free(nina_file);
		return;
	}
	gchar *dirname = g_path_get_dirname(nina_file);
	if (g_strcmp0(dirname, com.wd) != 0) {	// Tests if the file is in the CWD
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The current comparison stars file is not located in the CWD"));
		siril_file_chooser_set_current_folder_path(GTK_FILE_CHOOSER(file_chooser), com.wd);
		/* GTK4: gtk_file_chooser_unselect_all removed */;
		g_free(dirname);
		g_free(nina_file);
		return;
	}
	g_free(dirname);
	if (!has_wcs(gfit)) {	// Is the current image properly plate solved
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
	gboolean use_c1 = siril_toggle_get_active(GTK_WIDGET(use_comp1));
	gboolean use_c2 = siril_toggle_get_active(GTK_WIDGET(use_comp2));

	if (parse_nina_stars_file_using_WCS(args, nina_file, use_c1, use_c2, gfit)) {
		// fail
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Something went wrong while saving plot"));
		free(args);
		g_free(nina_file);
		return;
	}
	g_free(nina_file);
	args->seq = &com.seq;
	args->layer = layer;
	args->display_graph = siril_toggle_get_active(GTK_WIDGET(display_curve));
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
	gtk_editable_set_text(GTK_EDITABLE(apert_value), g_strdup_printf("%1.2lf", radius_value));
	gtk_editable_set_text(GTK_EDITABLE(inner_value), g_strdup_printf("%1.2lf", com.pref.phot_set.inner));
	gtk_editable_set_text(GTK_EDITABLE(outer_value), g_strdup_printf("%1.2lf", com.pref.phot_set.outer));

}
