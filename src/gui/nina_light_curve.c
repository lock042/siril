/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include "gui/plot.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *file_chooser = NULL;
static GtkWidget *use_comp1 = NULL;
static GtkWidget *use_comp2 = NULL;
static GtkWidget *display_curve = NULL;

static void on_nina_lc_response(GtkDialog* self, gint response_id, gpointer user_data);

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Light curve with NINA star list"), NULL,
			0, _("_OK"), GTK_RESPONSE_ACCEPT, "_Cancel", GTK_RESPONSE_REJECT, NULL);
	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_nina_lc_response), NULL);

	file_chooser = gtk_file_chooser_button_new (_("Select the NINA star list file"),
			GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(file_chooser), com.wd);
	g_object_set(G_OBJECT(file_chooser), "margin", 15, NULL);
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("CSV file (*.csv)"));
	gtk_file_filter_add_pattern(f, "*.csv");
	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(file_chooser), f);

	GtkWidget *label = gtk_label_new(_("Process a sequence to get a light curve on a star using the list of reference stars created by the NINA exoplanet plugin"));
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	g_object_set(G_OBJECT(label), "margin", 15, NULL);

	use_comp1 = gtk_check_button_new_with_label(_("Use comparative stars selected for their color"));
	use_comp2 = gtk_check_button_new_with_label(_("Use comparative stars selected by the AAVSO"));
	display_curve = gtk_check_button_new_with_label(_("Display the light curve"));
	g_object_set(G_OBJECT(use_comp1), "margin-left", 6, NULL);
	g_object_set(G_OBJECT(use_comp1), "margin-right", 6, NULL);
	g_object_set(G_OBJECT(use_comp2), "margin-left", 6, NULL);
	g_object_set(G_OBJECT(use_comp2), "margin-right", 6, NULL);
	g_object_set(G_OBJECT(display_curve), "margin", 6, NULL);
	gtk_widget_set_tooltip_text(use_comp1, _("Color similar to the target mean they will get extincted the same way by the changing atmosphere"));
	gtk_widget_set_tooltip_text(use_comp2, _("The AAVSO gives stars that are know to not be variable"));
	gtk_widget_set_tooltip_text(display_curve, _("if not checked, a PNG image of the graph will be generated instead"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(use_comp1), TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(use_comp2), TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(display_curve), TRUE);

	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 20);
	gtk_container_add(GTK_CONTAINER(content_area), label);
	gtk_container_add(GTK_CONTAINER(content_area), file_chooser);
	gtk_container_add(GTK_CONTAINER(content_area), use_comp1);
	gtk_container_add(GTK_CONTAINER(content_area), use_comp2);
	gtk_container_add(GTK_CONTAINER(content_area), display_curve);
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
	gtk_widget_hide(dialog);
	if (response_id != GTK_RESPONSE_ACCEPT) {
		return;
	}
	gchar *nina_file = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
	if (!nina_file)
		return;

	if (!has_wcs(&gfit)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The currently loaded image must be plate solved"));
		return;
	}
	purge_temp_user_catalogue();
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
		// TODO: if clicked on cancel, do not continue
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
		return;
	}

	args->seq = &com.seq;
	args->layer = layer;
	args->display_graph = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(display_curve));
	siril_debug_print("starting PSF analysis of %d stars\n", args->nb);

	start_in_new_thread(light_curve_worker, args);
}
