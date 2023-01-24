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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/photometric_cc.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/registration_preview.h"
#include "io/catalogues.h"
#include "photometric_cc.h"

static rectangle get_bkg_selection();
void on_combophoto_catalog_changed(GtkComboBox *combo, gpointer user_data);

static void start_photometric_cc() {
	GtkToggleButton *auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));
	GtkToggleButton *force_platesolve_button = GTK_TOGGLE_BUTTON(lookup_widget("force_astrometry_button"));
	gboolean plate_solve;

	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("There is no valid WCS information in the header. Let's make a plate solving.\n"), "salmon");
		plate_solve = TRUE;
	} else {
		plate_solve = gtk_toggle_button_get_active(force_platesolve_button);
		if (plate_solve)
			siril_log_message(_("Plate solving will be recomputed for image\n"));
		else siril_log_message(_("Existing plate solve (WCS information) will be resused for image\n"));
	}

	struct astrometry_data *args = NULL;
	struct photometric_cc_data *pcc_args = calloc(1, sizeof(struct photometric_cc_data));
	pcc_args->catalog = get_photometry_catalog_from_GUI();
	pcc_args->use_local_cat = FALSE;
	if (local_catalogues_available()) {
		if (pcc_args->catalog == CAT_NOMAD) {
			siril_debug_print("using local star catalogues\n");
			pcc_args->use_local_cat = TRUE;
		}
		else siril_log_message(_("Using remote APASS instead of local NOMAD catalogue\n"));
	}
	if (plate_solve) {
		args = calloc(1, sizeof(struct astrometry_data));
		args->fit = &gfit;

		args->for_photometry_cc = TRUE;
		args->use_local_cat = pcc_args->use_local_cat;

		args->pcc = pcc_args;
		args->pcc->fit = &gfit;
		args->pcc->bg_auto = gtk_toggle_button_get_active(auto_bkg);
		args->pcc->bg_area = get_bkg_selection();
	}

	pcc_args->fit = &gfit;
	pcc_args->bg_auto = gtk_toggle_button_get_active(auto_bkg);
	pcc_args->mag_mode = LIMIT_MAG_AUTO;

	set_cursor_waiting(TRUE);

	if (plate_solve) {
		if (!fill_plate_solver_structure_from_GUI(args)) {
			pcc_args->mag_mode = args->mag_mode;
			pcc_args->magnitude_arg = args->magnitude_arg;
			start_in_new_thread(plate_solver, args);
		}
	} else {
		get_limit_mag_from_GUI(&pcc_args->mag_mode, &pcc_args->magnitude_arg);
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_in_new_thread(photometric_cc_standalone, pcc_args);
	}
}

static rectangle get_bkg_selection() {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	rectangle black_selection;
	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);
	return black_selection;
}

static gboolean is_selection_ok() {
	rectangle selection = get_bkg_selection();
	return selection.w != 0 && selection.h != 0;
}

/****
 * PUBLIC FUNCTIONS
 */

void initialize_photometric_cc_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg,
			*catalog_label_pcc, *force_platesolve, *lasnet;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	force_platesolve = lookup_widget("force_astrometry_button");
	lasnet = lookup_widget("localasnet_check_button");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_ips, FALSE);
	gtk_widget_set_visible(catalog_box_pcc, TRUE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(force_platesolve, TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lasnet), FALSE);
	gtk_widget_set_visible(lasnet, FALSE);
	gtk_widget_grab_focus(button_cc_ok);

	gtk_window_set_title(parent, _("Photometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

	on_combophoto_catalog_changed(GTK_COMBO_BOX(catalog_box_pcc), NULL);

	// not sure about this one. Fails for a lot of images
	//gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("downsample_ips_button")), gfit.rx > 6000);
}

int get_photometry_catalog_from_GUI() {
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxPCCCatalog"));
	if (gtk_combo_box_get_active(box) == 1)
		return CAT_APASS;
	return CAT_NOMAD;
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_button_cc_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_bkg;

	auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	if ((!gtk_toggle_button_get_active(auto_bkg)) && (!is_selection_ok())) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
	}
	else start_photometric_cc();
}

void on_button_cc_bkg_auto_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *box_cc_manual_bkg;

	box_cc_manual_bkg = lookup_widget("box_cc_manual_bkg");
	gtk_widget_set_sensitive(box_cc_manual_bkg, !gtk_toggle_button_get_active(button));
}

void on_button_cc_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_cc_bkg_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_cc_bkg_value[0]) {
		selection_cc_bkg_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_cc_bkg_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_cc_bkg_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_cc_bkg_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	gtk_spin_button_set_value(selection_cc_bkg_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_cc_bkg_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_cc_bkg_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_cc_bkg_value[3], com.selection.h);

	delete_selected_area(); // needed because we don't want the selection being used for astrometry
}

void on_combophoto_catalog_changed(GtkComboBox *combo, gpointer user_data) {
	static gboolean have_local_cat = FALSE;
	static GtkLabel *photocat_label = NULL;
	if (!photocat_label) {
		photocat_label = GTK_LABEL(lookup_widget("photometric_catalog_label"));
		have_local_cat = local_catalogues_available();
	}
	if (gtk_combo_box_get_active(combo) == 1 || !have_local_cat)
		gtk_label_set_text(photocat_label, _("(online catalogue)"));
	else gtk_label_set_text(photocat_label, _("(local catalogue)"));
}

