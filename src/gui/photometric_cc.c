/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/siril_plot.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/registration_preview.h"
#include "io/remote_catalogues.h"
#include "io/local_catalogues.h"
#include "photometric_cc.h"
#include "io/image_format_fits.h"

static const cmsCIEXYZ D65 = {0.95045471, 1.0, 1.08905029};
static const cmsCIEXYZ D50 = {0.964199999, 1.000000000, 0.824899998};
static const cmsCIExyY D58 = {0.344994428, 0.35152261, 1.0}; // Sun as white point, modelled as a Black Body
// (note using the Black Body locus gives a slightly different result than using the Daylight locus but
// is *probably* what is wanted here.
static gboolean spcc_filters_initialized = FALSE;
static rectangle get_bkg_selection();
void on_combophoto_catalog_changed(GtkComboBox *combo, gpointer user_data);
void set_spcc_args(struct photometric_cc_data *args);
void get_whitepoint_from_ui(struct photometric_cc_data *args);
void populate_spcc_combos();
void on_spcc_toggle_sensor_type_toggled(GtkToggleButton *button, gpointer user_data);

static void start_photometric_cc(gboolean spcc) {
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
	if (spcc) {
		pcc_args->catalog = CAT_GAIADR3_DIRECT;
		siril_log_message(_("Using Gaia DR3 for SPCC\n"));
		pcc_args->spcc = TRUE;
		set_spcc_args(pcc_args);
		get_whitepoint_from_ui(pcc_args);
	} else {
		pcc_args->catalog = get_photometry_catalog_from_GUI();
		pcc_args->spcc = FALSE;
	}
	if (local_catalogues_available()) {
		if (pcc_args->catalog == CAT_NOMAD) {
			pcc_args->catalog = CAT_LOCAL;
			siril_debug_print("using local star catalogues\n");
		}
		else siril_log_message(_("Using remote APASS instead of local NOMAD catalogue\n"));
	}
	if (plate_solve) {
		args = calloc(1, sizeof(struct astrometry_data));
		args->fit = &gfit;
		args->for_photometry_cc = !spcc;
		args->for_photometry_spcc = spcc;


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
		get_mag_settings_from_GUI(&pcc_args->mag_mode, &pcc_args->magnitude_arg);
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
	GtkWidget *button_ips_ok, *button_cc_ok, *button_spcc_ok, *catalog_label,
			*astrometry_catalog_label, *pcc_catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *stardet,
			*catalog_label_pcc, *force_platesolve, *lasnet, *spcc_options,
			*labelIPScatparams;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	button_spcc_ok = lookup_widget("button_spcc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	astrometry_catalog_label = lookup_widget("astrometry_catalog_label");
	pcc_catalog_label = lookup_widget("photometric_catalog_label");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	force_platesolve = lookup_widget("force_astrometry_button");
	lasnet = lookup_widget("localasnet_check_button");
	spcc_options = lookup_widget("spcc_options");
	stardet = lookup_widget("Frame_IPS_star_detection");
	labelIPScatparams = lookup_widget("labelIPSCatalogParameters");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(button_spcc_ok, FALSE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(astrometry_catalog_label, TRUE);
	gtk_widget_set_visible(pcc_catalog_label, TRUE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(stardet, TRUE);
	gtk_widget_set_visible(catalog_box_ips, FALSE);
	gtk_widget_set_visible(catalog_box_pcc, TRUE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(force_platesolve, TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lasnet), FALSE);
	gtk_widget_set_visible(lasnet, FALSE);
	gtk_widget_set_visible(spcc_options, FALSE);
	gtk_widget_grab_focus(button_cc_ok);
	gtk_expander_set_expanded(GTK_EXPANDER(labelIPScatparams), TRUE);

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
	gtk_label_set_text(GTK_LABEL(lookup_widget("astrometry_catalog_label")), "");

	// not sure about this one. Fails for a lot of images
	//gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("downsample_ips_button")), gfit.rx > 6000);
}

void initialize_spectrophotometric_cc_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *button_spcc_ok, *catalog_label,
			*astrometry_catalog_label, *pcc_catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *stardet,
			*catalog_label_pcc, *force_platesolve, *lasnet, *spcc_options,
			*labelIPScatparams;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	button_spcc_ok = lookup_widget("button_spcc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	astrometry_catalog_label = lookup_widget("astrometry_catalog_label");
	pcc_catalog_label = lookup_widget("photometric_catalog_label");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	force_platesolve = lookup_widget("force_astrometry_button");
	lasnet = lookup_widget("localasnet_check_button");
	spcc_options = lookup_widget("spcc_options");
	stardet = lookup_widget("Frame_IPS_star_detection");
	labelIPScatparams = lookup_widget("labelIPSCatalogParameters");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, FALSE);
	gtk_widget_set_visible(button_spcc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(astrometry_catalog_label, FALSE);
	gtk_widget_set_visible(pcc_catalog_label, FALSE);
	gtk_widget_set_visible(catalog_label_pcc, FALSE);
	gtk_widget_set_visible(stardet, FALSE);
	gtk_widget_set_visible(catalog_box_ips, FALSE);
	gtk_widget_set_visible(catalog_box_pcc, FALSE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(force_platesolve, TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lasnet), FALSE);
	gtk_widget_set_visible(lasnet, FALSE);
	gtk_widget_set_visible(spcc_options, TRUE);
	gtk_widget_grab_focus(button_cc_ok);
	gtk_expander_set_expanded(GTK_EXPANDER(labelIPScatparams), FALSE);
	gtk_expander_set_expanded(GTK_EXPANDER(spcc_options), TRUE);

	gtk_window_set_title(parent, _("Spectrophotometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

	on_combophoto_catalog_changed(GTK_COMBO_BOX(catalog_box_pcc), NULL);
	gtk_label_set_text(GTK_LABEL(lookup_widget("astrometry_catalog_label")), "");
	on_spcc_toggle_sensor_type_toggled(GTK_TOGGLE_BUTTON(lookup_widget("spcc_toggle_sensor_type")), NULL);
	populate_spcc_combos();
}

int get_photometry_catalog_from_GUI() {
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxPCCCatalog"));
	if (gtk_combo_box_get_active(box) == 2)
		return CAT_GAIADR3;
	else if (gtk_combo_box_get_active(box) == 1)
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
	else start_photometric_cc(FALSE);
}

void on_button_spcc_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_bkg;

	auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	if ((!gtk_toggle_button_get_active(auto_bkg)) && (!is_selection_ok())) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
	}
	else start_photometric_cc(TRUE);
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
	if (gtk_combo_box_get_active(combo) > 0 || !have_local_cat) // 1 = APASS
		gtk_label_set_text(photocat_label, _("(online catalogue)"));
	else gtk_label_set_text(photocat_label, _("(local catalogue)"));
}

void set_spcc_args(struct photometric_cc_data *args) {
	GtkWidget *mono_sensor_check = lookup_widget("spcc_toggle_sensor_type");
	args->spcc_mono_sensor = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mono_sensor_check));
	GtkWidget *monosensor = lookup_widget("combo_spcc_sensors_mono");
	GtkWidget *oscsensor = lookup_widget("combo_spcc_sensors_osc");
	GtkWidget *filters_r = lookup_widget("combo_spcc_filters_r");
	GtkWidget *filters_g = lookup_widget("combo_spcc_filters_g");
	GtkWidget *filters_b = lookup_widget("combo_spcc_filters_b");
	GtkWidget *filters_osc = lookup_widget("combo_spcc_filters_osc");
	GtkWidget *osc_filters_enable = lookup_widget("osc_filters_enable");
	GtkWidget *max_stars_spin = lookup_widget("SPCC_max_stars");

	args->selected_sensor_m = gtk_combo_box_get_active(GTK_COMBO_BOX(monosensor));
	args->selected_sensor_osc = gtk_combo_box_get_active(GTK_COMBO_BOX(oscsensor));
	args->selected_filter_r = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_r));
	args->selected_filter_g = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_g));
	args->selected_filter_b = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_b));
	args->selected_filter_osc = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_osc));
	args->use_osc_filter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(osc_filters_enable));
	args->max_spcc_stars = gtk_spin_button_get_value(GTK_SPIN_BUTTON(max_stars_spin));

}

void get_whitepoint_from_ui(struct photometric_cc_data *args) {
	GtkWidget *wp = lookup_widget("combo_spcc_whitepoint");
	wp_t wp_index = gtk_combo_box_get_active(GTK_COMBO_BOX(wp));
	switch (wp_index) {
		case WP_D50:
			memcpy(&args->whitepoint, &D50, sizeof(cmsCIExyY));
			break;
		case WP_D65:
			memcpy(&args->whitepoint, &D65, sizeof(cmsCIExyY));
			break;
		case WP_SOL:
			memcpy(&args->whitepoint, &D58, sizeof(cmsCIExyY));
			break;
		case WP_GAL_AVGSPIRAL:
		case WP_GAL_AVGELLIPTICAL:
			// TODO: This currently just returns D50, update this based on data
			memcpy(&args->whitepoint, &D50, sizeof(cmsCIExyY));
			break;
	}
}

void fill_combo_from_glist(gchar *comboname, GList *list, int channel) {
	GtkComboBox *combo;
	GtkListStore *store;
	GtkTreeIter iter;

	combo=GTK_COMBO_BOX(lookup_widget(comboname));
	// Clear the model
	gtk_combo_box_set_model(combo, NULL);

	GtkTreeModel *model;

	GtkCellRenderer *renderer=gtk_cell_renderer_text_new();

	model=GTK_TREE_MODEL((store=gtk_list_store_new(1,G_TYPE_STRING)));
	gtk_combo_box_set_model(combo, model);
	gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo),renderer,TRUE);
	gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo),renderer,"text",0,NULL);

	GList *iterator = list;

	if (list == com.spcc_data.osc_sensors) {
		while (iterator) {
			osc_sensor *object = (osc_sensor*) iterator->data;
			gtk_list_store_append(store, &iter);
			gtk_list_store_set(store,&iter,0,object->channel[0].model,-1);
			iterator = iterator->next;
		}

} else {
		while (iterator) {
			// Easier, just add all objects by name
			spcc_object *object = (spcc_object*) iterator->data;
			gtk_list_store_append(store,&iter);
			gtk_list_store_set(store,&iter,0,object->name,-1);
			iterator = iterator->next;
		}
	}
	gtk_combo_box_set_active(combo, 0);
}

void populate_spcc_combos() {
	// Initialize filters if required
	if (!spcc_filters_initialized) {
		load_all_spcc_metadata();
		fill_combo_from_glist("combo_spcc_filters_r", com.spcc_data.mono_filters, RED);
		fill_combo_from_glist("combo_spcc_filters_g", com.spcc_data.mono_filters, GREEN);
		fill_combo_from_glist("combo_spcc_filters_b", com.spcc_data.mono_filters, BLUE);
		fill_combo_from_glist("combo_spcc_filters_osc", com.spcc_data.osc_filters, -1);
		fill_combo_from_glist("combo_spcc_sensors_mono", com.spcc_data.mono_sensors, -1);
		fill_combo_from_glist("combo_spcc_sensors_osc", com.spcc_data.osc_sensors, -1);
		spcc_filters_initialized = TRUE;
	}
}

void on_spcc_toggle_sensor_type_toggled(GtkToggleButton *button, gpointer user_data) {
	int state = gtk_toggle_button_get_active(button);
	GtkWidget *widget;
	gtk_button_set_label(GTK_BUTTON(button), state ? _("Mono Sensor") : _("OSC Sensor"));
	widget = lookup_widget("label_spcc_sensors_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("combo_spcc_sensors_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("details_spcc_sensors_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("label_spcc_filters_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("combo_spcc_filters_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("details_spcc_filters_osc");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("osc_filters_enable");
	gtk_widget_set_visible(widget, !state);
	widget = lookup_widget("combo_spcc_sensors_mono");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("label_spcc_sensors_mono");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("details_spcc_sensors_mono");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("combo_spcc_filters_r");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("label_spcc_filters_r");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("details_spcc_filters_r");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("combo_spcc_filters_g");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("label_spcc_filters_g");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("details_spcc_filters_g");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("combo_spcc_filters_b");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("label_spcc_filters_b");
	gtk_widget_set_visible(widget, state);
	widget = lookup_widget("details_spcc_filters_b");
	gtk_widget_set_visible(widget, state);
}

void on_osc_filter_enable_toggled(GtkToggleButton *button, gpointer user_data) {
	int state = gtk_toggle_button_get_active(button);
	GtkWidget *widget = lookup_widget("combo_spcc_filters_osc");
	gtk_widget_set_sensitive(widget, state);
}

static spcc_object *cbdata = NULL;

void on_spcc_details_plot_clicked(GtkButton *button, gpointer user_data);

void on_spcc_details_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("spcc_details");
	gtk_widget_set_visible(win, FALSE);
	GtkWidget *widget = GTK_WIDGET(button);
	GtkComboBox *combo = NULL;
	int n;
	GList *list = NULL;
	spcc_object *object = NULL;
	if (widget == lookup_widget("details_spcc_sensors_osc")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_osc"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.osc_sensors, n);
	} else if (widget == lookup_widget("details_spcc_sensors_mono")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_mono"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.mono_sensors, n);
	} else if (widget == lookup_widget("details_spcc_filters_osc")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_osc"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.osc_filters, n);
	} else if (widget == lookup_widget("details_spcc_filters_r")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_r"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.mono_filters, n);
	} else if (widget == lookup_widget("details_spcc_filters_g")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_g"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.mono_filters, n);
	} else if (widget == lookup_widget("details_spcc_filters_b")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_b"));
		n = gtk_combo_box_get_active(combo);
		list = g_list_nth(com.spcc_data.mono_filters, n);
	}
	// For OSC sensors which use the osc_sensor data structure this is a bit cheeky
	// but it works because the first element of the struct is a spcc_object, and
	// saves handling them differently.
	object = (spcc_object*) list->data;

	GtkLabel *label = GTK_LABEL(lookup_widget("spcc_details_name"));
	gtk_label_set_text(label, object->type == 2 ? object->model : object->name);
	label = GTK_LABEL(lookup_widget("spcc_details_mfr"));
	gtk_label_set_text(label, object->manufacturer);
	label = GTK_LABEL(lookup_widget("spcc_details_version"));
	gchar *version_text = g_strdup_printf("%d", object->version);
	gtk_label_set_text(label, version_text);
	g_free(version_text);
	label = GTK_LABEL(lookup_widget("spcc_details_source"));
	gtk_label_set_text(label, object->source);
	label = GTK_LABEL(lookup_widget("spcc_details_nsamples"));
	gchar *nsamples_text = g_strdup_printf("%d", object->n);
	gtk_label_set_text(label, nsamples_text);
	g_free(nsamples_text);
	cbdata = object;

	gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(lookup_widget("IPS_dialog")));
	/* Here this is wanted that we do not use siril_open_dialog */
	gtk_widget_show(win);
}

void on_spcc_details_plot_clicked(GtkButton *button, gpointer user_data) {
	siril_plot_data *spl_data = NULL;
	gboolean is_osc_sensor = (cbdata->type == 2);
	spl_data = malloc(sizeof(siril_plot_data));
	init_siril_plot_data(spl_data);
	siril_plot_set_xlabel(spl_data, _("Wavelength / nm"));
	siril_plot_set_savename(spl_data, "SPCC_data");
	if (is_osc_sensor) {
		osc_sensor *osc = (osc_sensor*) cbdata;
		gchar *title = g_strdup_printf(_("SPCC Data: %s"), osc->channel[0].model);
		siril_plot_set_title(spl_data, title);
		siril_plot_set_ylabel(spl_data, _("Quantum Efficiency"));
		for (int i = 0 ; i <3 ; i++) {
			load_spcc_object_arrays(&osc->channel[i]);
			gchar *spl_legend = g_strdup(osc->channel[i].name);
			siril_plot_add_xydata(spl_data, spl_legend, osc->channel[i].n, osc->channel[i].x, osc->channel[i].y, NULL, NULL);
			siril_plot_set_nth_color(spl_data, i+1, (double[3]){(double) i == 0, (double) i == 1, (double) i == 2});
			g_free(spl_legend);
			spcc_object_free_arrays(&osc->channel[i]);
		}
		g_free(title);
	} else {
		spcc_object *object = (spcc_object*) cbdata;
		load_spcc_object_arrays(object);
		gchar *title = g_strdup_printf(_("SPCC Data: %s"), object->name);
		siril_plot_set_title(spl_data, title);
		g_free(title);
		if (object->type == 1 || object->type == 2)
			siril_plot_set_ylabel(spl_data, _("Quantum Efficiency"));
		else
			siril_plot_set_ylabel(spl_data, _("Transmittance"));

		gchar *spl_legend = g_strdup(object->name);

		siril_plot_add_xydata(spl_data, spl_legend, object->n, object->x, object->y, NULL, NULL);
		spcc_object_free_arrays(object);
	}
	cbdata = NULL;
	siril_add_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}

void on_spcc_details_close_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("spcc_details");
	gtk_widget_hide(win);
}
