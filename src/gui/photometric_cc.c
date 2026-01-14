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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "core/siril_log.h"
#include "core/processing.h"
#include "core/siril_networking.h"
#include "algos/photometric_cc.h"
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/siril_plot.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "io/local_catalogues.h"
#include "io/healpix/healpix_cat.h"
#include "photometric_cc.h"

#define MIN_PLOT 336.0
#define MAX_PLOT 1020.0

// Uncomment the next line to enable more verbose debug reporting
//#define SPCC_DEBUG_TEST

#ifdef SPCC_DEBUG_TEST
#define DEBUG_SPCC 1
#else
#define DEBUG_SPCC 0
#endif
#define spcc_debug_print(fmt, ...) \
	do { if (DEBUG_TEST && DEBUG_SPCC) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

// Array of primary + mirrors for the remote catalogue
gchar **spcc_mirrors = NULL;

static gboolean spcc_filters_initialized = FALSE;
static int get_spcc_catalog_from_GUI();
static rectangle get_bkg_selection();
void on_combophoto_catalog_changed(GtkComboBox *combo, gpointer user_data);
void on_ComboBoxSPCCCatalog_changed(GtkComboBox *combo, gpointer user_data);
static void set_bg_sigma(struct photometric_cc_data *args);
static int set_spcc_args(struct photometric_cc_data *args);
gboolean populate_spcc_combos(gpointer user_data);
void on_spcc_toggle_nb_toggled(GtkToggleButton *button, gpointer user_data);
void on_spcc_sensor_switch_state_set(GtkSwitch *widget, gboolean state, gpointer user_data);

void initialize_spcc_mirrors() {
	if (spcc_mirrors)
		g_strfreev(spcc_mirrors);

	spcc_mirrors = g_new(gchar*, 2);
	spcc_mirrors[0] = g_strdup("https://zenodo.org/records/17988559/files");
	spcc_mirrors[1] = NULL;
}

void reset_spcc_filters() {
	spcc_filters_initialized = FALSE;
}

void on_buttonPCC_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("s_pcc_dialog");
}

gboolean s_pcc_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("s_pcc_dialog");
	return TRUE;
}

static void get_mag_settings_from_GUI(limit_mag_mode *mag_mode, double *magnitude_arg) {
	GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("S_PCC_Mag_Limit"));
	gboolean autob = gtk_toggle_button_get_active(autobutton);
	if (autob)
		*mag_mode = LIMIT_MAG_AUTO;
	else {
		GtkSpinButton *magButton = GTK_SPIN_BUTTON(lookup_widget("GtkSpinPCC_Mag_Limit"));
		*magnitude_arg = gtk_spin_button_get_value(magButton);
		*mag_mode = LIMIT_MAG_ABSOLUTE;
	}
}

void on_S_PCC_Mag_Limit_toggled(GtkToggleButton *button, gpointer user) {
	GtkWidget *spinmag;

	spinmag = lookup_widget("GtkSpinPCC_Mag_Limit");
	gtk_widget_set_sensitive(spinmag, !gtk_toggle_button_get_active(button));
}

static gboolean end_gaiacheck_idle(gpointer p) {
	GtkWidget *image = lookup_widget("gaia_status_widget");
	gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_grey.svg");
	gtk_widget_set_tooltip_text(image, _("Checking SPCC remote catalogue status..."));
	gtk_widget_show(image);
	int resptime = GPOINTER_TO_INT(p);;
	gchar *text = NULL, *colortext = NULL;

	if (resptime == -1) {
		// Failed to fetch status
		text = g_strdup(_("The Gaia remote SPCC catalogue is not responding."));
		colortext = "red";
		gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_red.svg");
	} else {
		// Check if response contains "true" or "false"
		if (resptime < 500) {
			text = g_strdup_printf(_("Gaia remote SPCC catalogue available, response %d ms"), resptime);
			colortext = "green";
			gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_green.svg");
		} else if (resptime < 1000) {
			text = g_strdup_printf(_("Gaia remote SPCC catalogue slow, response %d ms"), resptime);
			colortext = "salmon";
			gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_yellow.svg");
		} else {
			text = g_strdup_printf(_("Gaia remote SPCC catalogue very slow, response %d ms"), resptime);
			colortext = "red";
			gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_red.svg");
		}
	}
	gtk_widget_show(image);
	gtk_widget_set_tooltip_text(image, text);
	siril_log_color_message("%s\n", colortext, text);
	g_free(text);
	return FALSE;
}

gpointer gaia_check(gpointer user_data) {
    int best_mirror_index = -1;
    int best_responsetime = INT32_MAX;
    int num_mirrors = 0;

    // Count mirrors
    while (spcc_mirrors[num_mirrors]) num_mirrors++;

    // Arrays to store results from parallel checks
    int *response_times = g_new(int, num_mirrors);

    // Parallel loop to check all mirrors
    #pragma omp parallel for
    for (int i = 0; i < num_mirrors; i++) {
        gchar *mirror_url = g_strdup_printf("%s/siril_cat1_healpix8_xpsamp_0.dat",
                                           spcc_mirrors[i]);
        response_times[i] = http_check(mirror_url);
        g_free(mirror_url);
    }
    #pragma omp barrier

    // Log results for all mirrors for transparency
    int working_mirrors = 0;
    for (int i = 0; i < num_mirrors; i++) {
        if (response_times[i] > -1) {
            working_mirrors++;
            siril_debug_print("Mirror %s: %d ms\n", spcc_mirrors[i], response_times[i]);
            if (response_times[i] < best_responsetime) {
                best_mirror_index = i;
                best_responsetime = response_times[i];
            }
        } else {
            siril_log_color_message(_("Mirror %s: Not responding\n"), "salmon", spcc_mirrors[i]);
        }
    }

    g_free(response_times);

    if (best_mirror_index != -1) {
        g_free(com.spcc_remote_catalogue);
        com.spcc_remote_catalogue = g_strdup(spcc_mirrors[best_mirror_index]);
        siril_log_message(_("Primary SPCC catalogue set to: %s (%d working mirrors available)\n"),
                               com.spcc_remote_catalogue, working_mirrors);
    }

    execute_idle_and_wait_for_it(end_gaiacheck_idle, GINT_TO_POINTER(best_responsetime));
    return NULL;
}

void check_gaia_archive_status() {
	if (!is_online()) {
		GtkWidget *image = lookup_widget("gaia_status_widget");
		gtk_image_set_from_resource(GTK_IMAGE(image), "/org/siril/ui/pixmaps/status_red.svg");
		const gchar *text = N_("Siril is offline or built without networking. The remote Gaia catalogue is unavailable.\n");
		gtk_widget_set_tooltip_text(image, text);
		siril_log_color_message("%s", "red", text);
		return;
	}
	g_thread_unref(g_thread_new("gaia-status-check", gaia_check, NULL));
}

static void start_photometric_cc(gboolean spcc) {
	GtkToggleButton *auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	if (!has_wcs(gfit)) {
		siril_log_color_message(_("There is no valid WCS information in the header. Please platesolve first.\n"), "red");
		return;
	}

	struct photometric_cc_data *pcc_args = calloc(1, sizeof(struct photometric_cc_data));
	set_bg_sigma(pcc_args);
	if (spcc) {
		pcc_args->catalog = get_spcc_catalog_from_GUI();
		siril_debug_print(_("Using Gaia DR3 for SPCC\n"));
		pcc_args->spcc = TRUE;
		if (set_spcc_args(pcc_args)) {
			free(pcc_args);
			return;
		}
	} else {
		pcc_args->catalog = get_photometry_catalog_from_GUI();
		pcc_args->spcc = FALSE;
	}

	pcc_args->fit = gfit;
	pcc_args->bg_auto = gtk_toggle_button_get_active(auto_bkg);
	pcc_args->bg_area = get_bkg_selection();
	pcc_args->mag_mode = LIMIT_MAG_AUTO;

	set_cursor_waiting(TRUE);

	get_mag_settings_from_GUI(&pcc_args->mag_mode, &pcc_args->magnitude_arg);
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!start_in_new_thread(photometric_cc_standalone, pcc_args)) {
		g_free(pcc_args->datalink_path);
		free(pcc_args);
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
	GtkWidget *button_cc_ok, *button_spcc_ok, *catalog_label,
			*pcc_catalog_label, *catalog_box_pcc, *catalog_box_spcc, *frame_cc_bkg,
			*catalog_label_pcc, *spcc_options, *spcc_do_plot, *spcc_nb_controls,
			*spcc_toggle_nb, *gaia_status_check;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_cc_ok = lookup_widget("button_cc_ok");
	button_spcc_ok = lookup_widget("button_spcc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	pcc_catalog_label = lookup_widget("photometric_catalog_label");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_box_spcc = lookup_widget("ComboBoxSPCCCatalog");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	spcc_options = lookup_widget("spcc_options");
	spcc_do_plot = lookup_widget("spcc_plot_fits");
	spcc_nb_controls = lookup_widget("spcc_nb_controls");
	spcc_toggle_nb = lookup_widget("spcc_toggle_nb");
	gaia_status_check = lookup_widget("button_gaia_status_check");

	parent = GTK_WINDOW(lookup_widget("s_pcc_dialog"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(button_spcc_ok, FALSE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(pcc_catalog_label, TRUE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_spcc, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(spcc_options, FALSE);
	gtk_widget_set_visible(spcc_do_plot, FALSE);
	gtk_widget_set_visible(spcc_nb_controls, FALSE);
	gtk_widget_set_visible(spcc_toggle_nb, FALSE);
	gtk_widget_set_visible(gaia_status_check, FALSE);
	gtk_widget_grab_focus(button_cc_ok);

	gtk_window_set_title(parent, _("Photometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

	gtk_combo_box_set_active(GTK_COMBO_BOX(catalog_box_pcc), local_gaia_available() ? 2 : (local_kstars_available() ? 0 : 2));
	gtk_label_set_text(GTK_LABEL(lookup_widget("astrometry_catalog_label")), "");
}

void populate_nb_spinbuttons() {
	GtkSpinButton *button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_r_wl"));
	gtk_spin_button_set_value(button, com.pref.spcc.red_wl);
	button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_g_wl"));
	gtk_spin_button_set_value(button, com.pref.spcc.green_wl);
	button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_b_wl"));
	gtk_spin_button_set_value(button, com.pref.spcc.blue_wl);
	button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_r_bw"));
	gtk_spin_button_set_value(button, com.pref.spcc.red_bw);
	button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_g_bw"));
	gtk_spin_button_set_value(button, com.pref.spcc.green_bw);
	button = GTK_SPIN_BUTTON(lookup_widget("spcc_nb_b_bw"));
	gtk_spin_button_set_value(button, com.pref.spcc.blue_bw);
}

void initialize_spectrophotometric_cc_dialog() {
	check_gaia_archive_status();
	GtkWidget *button_cc_ok, *button_spcc_ok, *catalog_label,
			*pcc_catalog_label, *catalog_box_pcc, *catalog_box_spcc, *frame_cc_bkg,
			*catalog_label_pcc, *spcc_options, *spcc_do_plot, *spcc_nb_controls,
			*spcc_toggle_nb, *gaia_status_check;
	GtkSwitch *monoselector;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_cc_ok = lookup_widget("button_cc_ok");
	button_spcc_ok = lookup_widget("button_spcc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	pcc_catalog_label = lookup_widget("photometric_catalog_label");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_box_spcc = lookup_widget("ComboBoxSPCCCatalog");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	spcc_options = lookup_widget("spcc_options");
	spcc_do_plot = lookup_widget("spcc_plot_fits");
	spcc_nb_controls = lookup_widget("spcc_nb_controls");
	spcc_toggle_nb = lookup_widget("spcc_toggle_nb");
	monoselector = GTK_SWITCH(lookup_widget("spcc_sensor_switch"));
	gaia_status_check = lookup_widget("button_gaia_status_check");
	GtkWidget *spcc_airmass_entry = lookup_widget("spcc_airmass_entry");
	GtkSpinButton *spcc_height = GTK_SPIN_BUTTON(lookup_widget("spcc_height"));

	parent = GTK_WINDOW(lookup_widget("s_pcc_dialog"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_cc_ok, FALSE);
	gtk_widget_set_visible(button_spcc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(pcc_catalog_label, TRUE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_pcc, FALSE);
	gtk_widget_set_visible(catalog_box_spcc, TRUE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(spcc_options, TRUE);
	gtk_widget_set_visible(spcc_do_plot, TRUE);
	gtk_widget_grab_focus(button_cc_ok);
	gtk_widget_set_visible(gaia_status_check, TRUE);
	gtk_widget_set_visible(spcc_nb_controls, com.pref.spcc.nb_mode);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(spcc_toggle_nb), com.pref.spcc.nb_mode);
	gtk_widget_set_visible(spcc_toggle_nb, TRUE);
	gtk_expander_set_expanded(GTK_EXPANDER(spcc_options), TRUE);

	gtk_window_set_title(parent, _("Spectrophotometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

	gtk_combo_box_set_active(GTK_COMBO_BOX(catalog_box_spcc), local_gaia_xpsamp_available() ? 1 : 0);

	gchar *tooltip = NULL;
	double airmass, centalt = gfit->keywords.centalt, height = gfit->keywords.siteelev;
	if (gfit->keywords.airmass > 0.0) {
		airmass = gfit->keywords.airmass;
		tooltip = g_strdup(N_("Airmass read from FITS header AIRMASS card"));
	} else if (centalt > 0.0 && centalt < 90.0) {
		airmass = compute_airmass(90.0 - centalt);
		tooltip = g_strdup(N_("Airmass computed from FITS header CENTALT card"));
	} else {
		airmass = compute_airmass(41.9);
		tooltip = g_strdup(N_("No airmass data available or computable. Estimating airmass based on average observation height of 48.1°"));
	}
	if (height != DEFAULT_DOUBLE_VALUE) {
		gtk_spin_button_set_value(spcc_height, height);
		gtk_widget_set_tooltip_text(GTK_WIDGET(spcc_height), N_("Observer height read from FITS header SITEELEV card"));
		gtk_widget_set_sensitive(GTK_WIDGET(spcc_height), FALSE);
	} else {
		gtk_spin_button_set_value(spcc_height, 10.0);
		gtk_widget_set_tooltip_text(GTK_WIDGET(spcc_height), "");
		gtk_widget_set_sensitive(GTK_WIDGET(spcc_height), TRUE);
	}
	gchar *txt = g_strdup_printf("%.3f", airmass);
	gtk_entry_set_text(GTK_ENTRY(spcc_airmass_entry), txt);
	gtk_widget_set_tooltip_text(spcc_airmass_entry, tooltip);
	g_free(txt);
	g_free(tooltip);

	on_ComboBoxSPCCCatalog_changed(GTK_COMBO_BOX(catalog_box_spcc), NULL);
	gtk_label_set_text(GTK_LABEL(lookup_widget("astrometry_catalog_label")), "");
	g_signal_handlers_block_by_func(G_OBJECT(monoselector), on_spcc_sensor_switch_state_set, NULL);
	gtk_switch_set_active(monoselector, com.pref.spcc.is_mono);
	g_signal_handlers_unblock_by_func(G_OBJECT(monoselector), on_spcc_sensor_switch_state_set, NULL);
	on_spcc_sensor_switch_state_set(monoselector, com.pref.spcc.is_mono, lookup_widget("spcc_label_switch"));
	g_signal_handlers_block_by_func(G_OBJECT(spcc_toggle_nb), on_spcc_toggle_nb_toggled, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(spcc_toggle_nb), com.pref.spcc.nb_mode);
	g_signal_handlers_unblock_by_func(G_OBJECT(spcc_toggle_nb), on_spcc_toggle_nb_toggled, NULL);
	populate_nb_spinbuttons();
}

int get_photometry_catalog_from_GUI() {
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxPCCCatalog"));
	if (gtk_combo_box_get_active(box) == 2) {
		if (local_gaia_available())
			return CAT_LOCAL_GAIA_ASTRO;
		else
			return CAT_GAIADR3;
	}
	else if (gtk_combo_box_get_active(box) == 1)
		return CAT_APASS;
	if (local_kstars_available())
		return CAT_LOCAL_KSTARS;
	return CAT_NOMAD;
	
}

int get_spcc_catalog_from_GUI() {
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxSPCCCatalog"));
	if (gtk_combo_box_get_active(box) == 1) {
		return CAT_LOCAL_GAIA_XPSAMP;
	} else {
		return CAT_REMOTE_GAIA_XPSAMP;
	}
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
	static GtkLabel *photocat_label = NULL;
	if (!photocat_label) {
		photocat_label = GTK_LABEL(lookup_widget("photometric_catalog_label"));
	}
	switch (gtk_combo_box_get_active(combo)) {
		case 0: // NOMAD
			gtk_label_set_text(photocat_label, local_kstars_available() ? _("(local catalogue)") : _("(online catalogue)"));
			break;
		case 1: // APASS
			gtk_label_set_text(photocat_label, _("(online catalogue)"));
			break;
		case 2: // GAIA
			gtk_label_set_text(photocat_label, local_gaia_available() ? ("(local catalogue)"): _("(online catalogue)"));
			break;
		default:
			break;
	}
}

void on_ComboBoxSPCCCatalog_changed(GtkComboBox *combo, gpointer user_data) {
	static GtkLabel *photocat_label = NULL;
	if (!photocat_label) {
		photocat_label = GTK_LABEL(lookup_widget("photometric_catalog_label"));
	}
	switch (gtk_combo_box_get_active(combo)) {
		case 0: // Gaia archive
			gtk_label_set_text(photocat_label, _("(online catalogue)"));
			break;
		case 1: // Local Gaia xp_sampled
			gtk_label_set_text(photocat_label, _("(local catalogue)"));
			break;
		default:
			break;
	}
}

static void set_bg_sigma(struct photometric_cc_data *args) {
	GtkWidget *lower_widget = lookup_widget("bg_tol_lower");
	GtkWidget *upper_widget = lookup_widget("bg_tol_upper");
	args->t0 = 0.f - gtk_spin_button_get_value(GTK_SPIN_BUTTON(lower_widget));
	args->t1 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(upper_widget));
}

static int set_spcc_args(struct photometric_cc_data *args) {
	GtkWidget *mono_sensor_check = lookup_widget("spcc_sensor_switch");
	GtkWidget *whiteref = lookup_widget("combo_spcc_whitepoint");
	GtkWidget *monosensor = lookup_widget("combo_spcc_sensors_mono");
	GtkWidget *oscsensor = lookup_widget("combo_spcc_sensors_osc");
	GtkWidget *filters_r = lookup_widget("combo_spcc_filters_r");
	GtkWidget *filters_g = lookup_widget("combo_spcc_filters_g");
	GtkWidget *filters_b = lookup_widget("combo_spcc_filters_b");
	GtkWidget *filters_osc = lookup_widget("combo_spcc_filters_osc");
	GtkWidget *filters_lpf = lookup_widget("combo_spcc_filters_lpf");
	GtkWidget *spcc_plot = lookup_widget("spcc_plot_fits");
	GtkWidget *nb_mode = lookup_widget("spcc_toggle_nb");
	GtkWidget *nb_r_wl = lookup_widget("spcc_nb_r_wl");
	GtkWidget *nb_g_wl = lookup_widget("spcc_nb_g_wl");
	GtkWidget *nb_b_wl = lookup_widget("spcc_nb_b_wl");
	GtkWidget *nb_r_bw = lookup_widget("spcc_nb_r_bw");
	GtkWidget *nb_g_bw = lookup_widget("spcc_nb_g_bw");
	GtkWidget *nb_b_bw = lookup_widget("spcc_nb_b_bw");
	GtkWidget *spcc_atmos_corr = lookup_widget("spcc_atmos_corr");
	GtkWidget *spcc_height = lookup_widget("spcc_height");
	GtkWidget *spcc_pressure = lookup_widget("spcc_pressure");
	GtkWidget *spcc_pressure_is_slp = lookup_widget("spcc_pressure_is_slp");

	args->atmos_corr = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(spcc_atmos_corr));
	args->atmos_obs_height = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spcc_height));
	args->atmos_pressure = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spcc_pressure));
	args->atmos_pressure_is_slp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(spcc_pressure_is_slp));
	args->spcc_mono_sensor = gtk_switch_get_active(GTK_SWITCH(mono_sensor_check));
	args->selected_white_ref = gtk_combo_box_get_active(GTK_COMBO_BOX(whiteref));
	args->selected_sensor_m = gtk_combo_box_get_active(GTK_COMBO_BOX(monosensor));
	args->selected_sensor_osc = gtk_combo_box_get_active(GTK_COMBO_BOX(oscsensor));
	args->selected_filter_r = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_r));
	args->selected_filter_g = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_g));
	args->selected_filter_b = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_b));
	args->selected_filter_osc = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_osc));
	GList *osc = g_list_nth(com.spcc_data.osc_sensors, args->selected_sensor_osc);
	if (osc) {
		osc_sensor *oscsen = (osc_sensor*) osc->data;
		args->is_dslr = oscsen->channel[0].is_dslr;
	} else {
		args->is_dslr = com.pref.spcc.is_dslr;
	}

	args->selected_filter_lpf = gtk_combo_box_get_active(GTK_COMBO_BOX(filters_lpf));
	args->do_plot = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(spcc_plot));
	args->nb_mode = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(nb_mode));
	args->nb_center[RLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_r_wl));
	args->nb_center[GLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_g_wl));
	args->nb_center[BLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_b_wl));
	args->nb_bandwidth[RLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_r_bw));
	args->nb_bandwidth[GLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_g_bw));
	args->nb_bandwidth[BLAYER] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(nb_b_bw));
	// Check for problems
	if (args->spcc_mono_sensor) {
		if (args->selected_sensor_m == -1 || args->selected_filter_r == -1 ||
						args->selected_filter_g == -1 || args->selected_filter_b == -1 ||
						args->selected_white_ref == -1) {
			siril_log_color_message(_("Error: ensure mono sensor, R, G and B filters and white reference are selected.\n"), "red");
			return 1;
		}
	} else {
		if (args->selected_sensor_osc == -1 || (args->selected_filter_osc == -1 && !args->nb_mode) ||
					(args->selected_filter_lpf == -1 && args->is_dslr == TRUE) ||
					args->selected_white_ref == -1) {
			siril_log_color_message(_("Error: ensure OSC sensor, OSC filter or NB filter values, and white reference are selected.\n"), "red");
			return 1;
		}
	}
	return 0;
}

int get_favourite_spccobject(GList *list, const gchar *favourite) {
	if (!list)
		return 0;

	GList *current = list;
	while (current != NULL) {
		spcc_object *haystack = current->data;
		if (haystack && g_strcmp0(haystack->name, favourite) == 0) {
			return g_list_position(list, current);  // Found a match, return the GList node
		}
		current = current->next;
	}
	return -1;  // No match found
}

int get_favourite_oscsensor(GList *list, const gchar *favourite) {
	if (!list)
		return 0;

	GList *current = list;
	while (current != NULL) {
		osc_sensor *haystack = current->data;
		if (g_strcmp0(haystack->channel[0].model, favourite) == 0) {
			return g_list_position(list, current);  // Found a match, return the GList node
		}
		current = current->next;
	}
	return -1;  // No match found
}

void on_spcc_combo_changed(GtkComboBox *combo, gpointer user_data);

/**
 * Validates a string for safe use in GTK widgets
 * @param str The string to validate
 * @return TRUE if string is valid, FALSE otherwise
 */
static gboolean is_valid_string(const char *str) {
	// Check for NULL pointer
	if (!str) {
		return FALSE;
	}

	// Check for empty string
	if (strlen(str) == 0) {
		return FALSE;
	}

	// Check for reasonable length (prevent buffer overflows/memory issues)
	if (strlen(str) > 1000) {
		return FALSE;
	}

	return TRUE;
}

void fill_combo_from_glist(GtkWidget *widget, GList *list, int channel, const gchar *favourite) {
	// Input validation
	if (!widget) {
		siril_debug_print("fill_combo_from_glist: widget is NULL\n");
		return;
	}

	if (!GTK_IS_COMBO_BOX_TEXT(widget)) {
		siril_debug_print("fill_combo_from_glist: widget is not a combo box text widget\n");
		return;
	}

	// Critical: Check if list is NULL
	if (!list) {
		siril_debug_print("fill_combo_from_glist: list is NULL, cannot populate combo\n");
		return;
	}

	GtkComboBoxText *combo = GTK_COMBO_BOX_TEXT(widget);

	// Validate that the combo box is still valid
	if (!GTK_IS_COMBO_BOX_TEXT(combo)) {
		siril_debug_print("fill_combo_from_glist: combo cast failed\n");
		return;
	}

	// Block signal handler so we don't call the callback while messing about with the combobox
	// Find the specific handler ID instead of using function-based blocking
	guint signal_id = g_signal_lookup("changed", GTK_TYPE_COMBO_BOX);
	if (signal_id == 0) {
		siril_debug_print("fill_combo_from_glist: failed to lookup 'changed' signal\n");
		return;
	}

	guint n_blocked = g_signal_handlers_block_matched(G_OBJECT(combo),
		G_SIGNAL_MATCH_ID,
		signal_id,
		0, NULL, NULL, NULL);

	// Clear the model
	gtk_combo_box_text_remove_all(combo);

	GList *iterator = list;
	int active_index = -1;
	int item_count = 0;

	// Validate list length to prevent infinite loops
	int max_items = 10000; // reasonable upper limit
	GList *temp_iter = list;
	int list_length = 0;
	while (temp_iter && list_length < max_items) {
		list_length++;
		temp_iter = temp_iter->next;
	}

	if (list_length >= max_items) {
		siril_debug_print("fill_combo_from_glist: list appears to be corrupted (too long or circular)\n");
		goto cleanup;
	}

	if (list == com.spcc_data.osc_sensors) {
		spcc_debug_print("fill_combo_from_glist: populating OSC sensors (%d items)\n", list_length);

		while (iterator && item_count < max_items) {
			// Validate iterator and its data
			if (!iterator->data) {
				siril_debug_print("fill_combo_from_glist: NULL data in osc_sensors list at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			osc_sensor *object = (osc_sensor *)iterator->data;

			// Comprehensive validation of osc_sensor object
			if (!object) {
				siril_debug_print("fill_combo_from_glist: NULL osc_sensor object at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			if (!object->channel[0].model) {
				siril_debug_print("fill_combo_from_glist: NULL model string in osc_sensor at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			// Validate string
			if (!is_valid_string(object->channel[0].model)) {
				siril_debug_print("fill_combo_from_glist: invalid string in osc_sensor at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			// Attempt to add to combo box
			gtk_combo_box_text_append_text(combo, object->channel[0].model);
			iterator = iterator->next;
			item_count++;
		}

		// Get favourite index with validation
		if (favourite) {
			active_index = get_favourite_oscsensor(list, favourite);
		}

	} else {
		spcc_debug_print("fill_combo_from_glist: populating SPCC objects (%d items)\n", list_length);

		while (iterator && item_count < max_items) {
			// Validate iterator and its data
			if (!iterator->data) {
				gchar *listname = NULL;
				if (list == com.spcc_data.osc_sensors)
					listname = g_strdup("ocs_sensors");
				else if (list == com.spcc_data.mono_sensors)
					listname = g_strdup("mono_sensors");
				else if (list == com.spcc_data.osc_filters)
					listname = g_strdup("osc_filters");
				else if (list == com.spcc_data.wb_ref)
					listname = g_strdup("wb_ref");
				else if (list == com.spcc_data.osc_lpf)
					listname = g_strdup("osc_lpf");
				else if (list == com.spcc_data.mono_filters[RLAYER])
					listname = g_strdup("mono filters (red)");
				else if (list == com.spcc_data.mono_filters[GLAYER])
					listname = g_strdup("mono filters (green)");
				else if (list == com.spcc_data.mono_filters[BLAYER])
					listname = g_strdup("mono filters (blue)");
				siril_debug_print("fill_combo_from_glist: NULL data in spcc_object list %s at item %d\n", listname, item_count);
				g_free(listname);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			spcc_object *object = (spcc_object *)iterator->data;

			// Comprehensive validation of spcc_object
			if (!object) {
				siril_debug_print("fill_combo_from_glist: NULL spcc_object at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			if (!object->name) {
				siril_debug_print("fill_combo_from_glist: NULL name in spcc_object at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			// Validate string
			if (!is_valid_string(object->name)) {
				siril_debug_print("fill_combo_from_glist: invalid string in spcc_object at item %d\n", item_count);
				iterator = iterator->next;
				item_count++;
				continue;
			}

			// Attempt to add to combo box
			gtk_combo_box_text_append_text(combo, object->name);
			iterator = iterator->next;
			item_count++;
		}

		// Get favourite index with validation
		if (favourite) {
			active_index = get_favourite_spccobject(list, favourite);
		}
	}

	spcc_debug_print("fill_combo_from_glist: successfully added %d items to combo\n", item_count);

cleanup:
	// Unblock handlers we blocked
	if (n_blocked > 0) {
		g_signal_handlers_unblock_matched(G_OBJECT(combo),
			G_SIGNAL_MATCH_ID,
			signal_id,
			0, NULL, NULL, NULL);
	}

	// Set active index with validation
	if (active_index >= 0 && active_index < item_count) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(combo), active_index);
		spcc_debug_print("fill_combo_from_glist: set active index to %d\n", active_index);
	} else if (active_index >= 0) {
		siril_debug_print("fill_combo_from_glist: favourite index %d is out of range (max: %d)\n", active_index, item_count - 1);
	}
}

gboolean populate_spcc_combos(gpointer user_data) {
	siril_debug_print("populate_spcc_combos: starting\n");

	// Validate we're on the main thread
	if (!g_main_context_is_owner(g_main_context_default())) {
		siril_debug_print("populate_spcc_combos: ERROR - not called from main thread!\n");
		goto end;
	}

	// Comprehensive validation of critical data structures before use
	if (!com.spcc_data.osc_sensors) {
		siril_debug_print("populate_spcc_combos: ERROR - osc_sensors is NULL\n");
		goto end;
	}
	if (!com.spcc_data.mono_sensors) {
		siril_debug_print("populate_spcc_combos: ERROR - mono_sensors is NULL\n");
		goto end;
	}
	if (!com.spcc_data.osc_filters) {
		siril_debug_print("populate_spcc_combos: ERROR - osc_filters is NULL\n");
		goto end;
	}
	if (!com.spcc_data.wb_ref) {
		siril_debug_print("populate_spcc_combos: ERROR - wb_ref is NULL\n");
		goto end;
	}
	if (!com.spcc_data.osc_lpf) {
		siril_debug_print("populate_spcc_combos: ERROR - osc_lpf is NULL\n");
		goto end;
	}
	// Check individual filter arrays
	if (!com.spcc_data.mono_filters[RLAYER]) {
		siril_debug_print("populate_spcc_combos: ERROR - red filters (RLAYER) is NULL\n");
		goto end;
	}
	if (!com.spcc_data.mono_filters[GLAYER]) {
		siril_debug_print("populate_spcc_combos: ERROR - green filters (GLAYER) is NULL\n");
		goto end;
	}
	if (!com.spcc_data.mono_filters[BLAYER]) {
		siril_debug_print("populate_spcc_combos: ERROR - blue filters (BLAYER) is NULL\n");
		goto end;
	}
	// Look up all the widgets here with validation
	GtkWidget *oscfilters = lookup_widget("combo_spcc_filters_osc");
	if (!oscfilters) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_filters_osc widget\n");
		goto end;
	}
	GtkWidget *rfilters = lookup_widget("combo_spcc_filters_r");
	if (!rfilters) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_filters_r widget\n");
		goto end;
	}
	GtkWidget *gfilters = lookup_widget("combo_spcc_filters_g");
	if (!gfilters) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_filters_g widget\n");
		goto end;
	}
	GtkWidget *bfilters = lookup_widget("combo_spcc_filters_b");
	if (!bfilters) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_filters_b widget\n");
		goto end;
	}
	GtkWidget *lpfilters = lookup_widget("combo_spcc_filters_lpf");
	if (!lpfilters) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_filters_lpf widget\n");
		goto end;
	}
	GtkWidget *monosensors = lookup_widget("combo_spcc_sensors_mono");
	if (!monosensors) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_sensors_mono widget\n");
		goto end;
	}
	GtkWidget *oscsensors = lookup_widget("combo_spcc_sensors_osc");
	if (!oscsensors) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_sensors_osc widget\n");
		goto end;
	}
	GtkWidget *whitepoint = lookup_widget("combo_spcc_whitepoint");
	if (!whitepoint) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find combo_spcc_whitepoint widget\n");
		goto end;
	}
	GtkWidget *switch_widget_lookup = lookup_widget("spcc_sensor_switch");
	if (!switch_widget_lookup) {
		siril_debug_print("populate_spcc_combos: ERROR - could not find spcc_sensor_switch widget\n");
		goto end;
	}
	GtkSwitch *switch_widget = GTK_SWITCH(switch_widget_lookup);
	if (!GTK_IS_SWITCH(switch_widget)) {
		siril_debug_print("populate_spcc_combos: ERROR - spcc_sensor_switch is not a GtkSwitch\n");
		goto end;
	}
	// Validate preferences exist (basic check)
	if (!com.pref.spcc.oscfilterpref) {
		siril_debug_print("populate_spcc_combos: WARNING - oscfilterpref is NULL\n");
	}

	// Initialize filters with error handling
	spcc_debug_print("populate_spcc_combos: populating OSC filters\n");
	fill_combo_from_glist(oscfilters, com.spcc_data.osc_filters, -1, com.pref.spcc.oscfilterpref);

	spcc_debug_print("populate_spcc_combos: populating red filters\n");
	fill_combo_from_glist(rfilters, com.spcc_data.mono_filters[RLAYER], RLAYER, com.pref.spcc.redpref);

	spcc_debug_print("populate_spcc_combos: populating green filters\n");
	fill_combo_from_glist(gfilters, com.spcc_data.mono_filters[GLAYER], GLAYER, com.pref.spcc.greenpref);

	spcc_debug_print("populate_spcc_combos: populating blue filters\n");
	fill_combo_from_glist(bfilters, com.spcc_data.mono_filters[BLAYER], BLAYER, com.pref.spcc.bluepref);

	spcc_debug_print("populate_spcc_combos: populating LP filters\n");
	fill_combo_from_glist(lpfilters, com.spcc_data.osc_lpf, -1, com.pref.spcc.lpfpref);

	spcc_debug_print("populate_spcc_combos: populating mono sensors\n");
	fill_combo_from_glist(monosensors, com.spcc_data.mono_sensors, -1, com.pref.spcc.monosensorpref);

	spcc_debug_print("populate_spcc_combos: populating OSC sensors\n");
	fill_combo_from_glist(oscsensors, com.spcc_data.osc_sensors, -1, com.pref.spcc.oscsensorpref);

	spcc_debug_print("populate_spcc_combos: populating white point\n");
	fill_combo_from_glist(whitepoint, com.spcc_data.wb_ref, -1, "Average Spiral Galaxy");

	spcc_debug_print("populate_spcc_combos: setting switch state\n");
	gtk_switch_set_active(switch_widget, com.pref.spcc.is_mono);

	siril_debug_print("populate_spcc_combos: completed successfully\n");
end:
	return FALSE;
}

gpointer populate_spcc_combos_async(gpointer user_data) {
	load_spcc_metadata_if_needed();
	// update combos back in the GTK thread
	gui_mutex_lock(); // force deconfliction with other GUI updates (specifically the script menu)
	execute_idle_and_wait_for_it(populate_spcc_combos, NULL);
	gui_mutex_unlock();
	return GINT_TO_POINTER(0);
}

static void set_osc_lpf_visibility() {
	gboolean sensor_is_osc = !gtk_switch_get_active(GTK_SWITCH(lookup_widget("spcc_sensor_switch")));
	gboolean is_dslr;
	GList *osc = g_list_nth(com.spcc_data.osc_sensors, gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_osc"))));
	if (osc) {
		osc_sensor *oscsensor = (osc_sensor*) osc->data;
		is_dslr = oscsensor->channel[0].is_dslr;
	} else {
		is_dslr = com.pref.spcc.is_dslr;
	}
	GtkWidget *widget = lookup_widget("label_spcc_filters_lpf");
	gtk_widget_set_visible(widget, (is_dslr && sensor_is_osc));
	widget = lookup_widget("combo_spcc_filters_lpf");
	gtk_widget_set_visible(widget, (is_dslr && sensor_is_osc));
	widget = lookup_widget("details_spcc_filters_lpf");
	gtk_widget_set_visible(widget, (is_dslr && sensor_is_osc));
}

void on_spcc_sensor_switch_state_set(GtkSwitch *widget, gboolean state, gpointer user_data) {
	com.pref.spcc.is_mono = state;
	GtkLabel *label = GTK_LABEL(user_data);
    GtkWidget *w;

    state ? gtk_label_set_text(label, _("mono")) : gtk_label_set_text(label, _("one-shot color (OSC)"));

	w = lookup_widget("label_spcc_sensors_osc");
	gtk_widget_set_visible(w, !state);
	w = lookup_widget("combo_spcc_sensors_osc");
	gtk_widget_set_visible(w, !state);
	w = lookup_widget("details_spcc_sensors_osc");
	gtk_widget_set_visible(w, !state);
	w = lookup_widget("label_spcc_filters_osc");
	gtk_widget_set_visible(w, !state && !com.pref.spcc.nb_mode);
	w = lookup_widget("combo_spcc_filters_osc");
	gtk_widget_set_visible(w, !state && !com.pref.spcc.nb_mode);
	w = lookup_widget("details_spcc_filters_osc");
	gtk_widget_set_visible(w, !state && !com.pref.spcc.nb_mode);
	w = lookup_widget("combo_spcc_sensors_mono");
	gtk_widget_set_visible(w, state);
	w = lookup_widget("label_spcc_sensors_mono");
	gtk_widget_set_visible(w, state);
	w = lookup_widget("details_spcc_sensors_mono");
	gtk_widget_set_visible(w, state);
	w = lookup_widget("combo_spcc_filters_r");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("label_spcc_filters_r");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("details_spcc_filters_r");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("combo_spcc_filters_g");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("label_spcc_filters_g");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("details_spcc_filters_g");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("combo_spcc_filters_b");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("label_spcc_filters_b");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	w = lookup_widget("details_spcc_filters_b");
	gtk_widget_set_visible(w, state && !com.pref.spcc.nb_mode);
	set_osc_lpf_visibility();
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
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.osc_sensors, n);
	} else if (widget == lookup_widget("details_spcc_sensors_mono")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_mono"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.mono_sensors, n);
	} else if (widget == lookup_widget("details_spcc_filters_osc")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_osc"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.osc_filters, n);
	} else if (widget == lookup_widget("details_spcc_filters_r")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_r"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.mono_filters[RLAYER], n);
	} else if (widget == lookup_widget("details_spcc_filters_g")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_g"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.mono_filters[GLAYER], n);
	} else if (widget == lookup_widget("details_spcc_filters_b")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_b"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.mono_filters[BLAYER], n);
	} else if (widget == lookup_widget("details_spcc_filters_lpf")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_lpf"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.osc_lpf, n);
	} else if (widget == lookup_widget("details_spcc_wb")) {
		combo = GTK_COMBO_BOX(lookup_widget("combo_spcc_whitepoint"));
		n = gtk_combo_box_get_active(combo);
		if (n == -1)
			goto no_selection;
		list = g_list_nth(com.spcc_data.wb_ref, n);
	}
	if (!list)
		goto no_selection;
	// For OSC sensors which use the osc_sensor data structure this is a bit cheeky
	// but it works because the first element of the struct is a spcc_object, and
	// saves handling them differently.
	object = (spcc_object*) list->data;

	GtkLabel *label = GTK_LABEL(lookup_widget("spcc_details_name"));
	gtk_label_set_text(label, object->type == OSC_SENSORS ? object->model : object->name);
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
	label = GTK_LABEL(lookup_widget("spcc_details_comment"));
	gboolean comment_visible;
	if (object->comment != NULL) {
		gtk_label_set_text(label, object->comment);
		comment_visible = TRUE;
	} else {
		comment_visible = FALSE;
	}
	gtk_widget_set_visible(GTK_WIDGET(label), comment_visible);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spcc_details_comment_label")), comment_visible);
	cbdata = object;

	gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(lookup_widget("s_pcc_dialog")));
	/* Here this is wanted that we do not use siril_open_dialog */
	gtk_widget_show(win);
	return;
	no_selection:
	siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), _("No selection made: no details to show."));
	return;
}

static void normalize_y(spcc_object *object, double min_wl, double max_wl) {
	int i = 0, j = 0;
	double maximum = -DBL_MAX;
	if (min_wl < object->x[0] || max_wl > object->x[object->n-1])
		return;

	while (object->x[i] < min_wl)
		i++;
	while (object->x[j] < max_wl)
		j++;
	for (int k = i ; k < j ; k++) {
		if (object->y[k] > maximum)
			maximum = object->y[k];
	}
	for (int k = 0 ; k < object->n ; k++)
		object->y[k] /= maximum;
}

void on_spcc_plot_all_clicked(GtkButton *button, gpointer user_data) {
	struct photometric_cc_data args = { 0 };
	if (set_spcc_args(&args))
		return;
	GList *sensor_list = args.spcc_mono_sensor ? com.spcc_data.mono_sensors : com.spcc_data.osc_sensors;
	GList *filter_list_r = com.spcc_data.mono_filters[RLAYER];
	GList *filter_list_g = com.spcc_data.mono_filters[GLAYER];
	GList *filter_list_b = com.spcc_data.mono_filters[BLAYER];
	GList *filter_list_lpf = com.spcc_data.osc_lpf;
	GList *filter_list_osc = com.spcc_data.osc_filters;
	GList *whiteref_list = com.spcc_data.wb_ref;
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data)
		return;
	siril_plot_set_xlabel(spl_data, _("Wavelength / nm"));
	siril_plot_set_savename(spl_data, "SPCC_data");
	siril_plot_set_title(spl_data, _("SPCC Data"));
	siril_plot_set_ylabel(spl_data, _("Quantum Efficiency / Transmittance / Rel. Photon Count"));

	int i = 0;
	if (args.spcc_mono_sensor) {
		spcc_object *sensor = NULL, *filter_r = NULL, *filter_g = NULL, *filter_b = NULL, *whiteref = NULL;
		if (args.selected_sensor_m >= 0 && args.selected_sensor_m < g_list_length (sensor_list))
			sensor = (spcc_object*) g_list_nth(sensor_list, args.selected_sensor_m)->data;
		if (!args.nb_mode) {
			if (args.selected_filter_r >= 0 && args.selected_filter_r < g_list_length (filter_list_r))
				filter_r = (spcc_object*) g_list_nth(filter_list_r, args.selected_filter_r)->data;
			if (args.selected_filter_g >= 0 && args.selected_filter_g < g_list_length (filter_list_g))
				filter_g = (spcc_object*) g_list_nth(filter_list_g, args.selected_filter_g)->data;
			if (args.selected_filter_b >= 0 && args.selected_filter_b < g_list_length (filter_list_b))
				filter_b = (spcc_object*) g_list_nth(filter_list_b, args.selected_filter_b)->data;
		}
		if (args.selected_white_ref >= 0 && args.selected_white_ref < g_list_length (whiteref_list))
			whiteref = (spcc_object*) g_list_nth(whiteref_list, args.selected_white_ref)->data;
		// there must not be any NULL spcc_object*s with non-NULL ones after them
		// (these should all be populated anyway)
		spcc_object* structs[6] = {  whiteref, sensor, filter_r, filter_g, filter_b, NULL };
		while (structs[i]) {
			load_spcc_object_arrays(structs[i]);
			if (structs[i] == whiteref)
				normalize_y(structs[i], 400, 700);
			gchar *spl_legend = g_strdup(structs[i]->name);
			siril_plot_add_xydata(spl_data, spl_legend, structs[i]->n, structs[i]->x, structs[i]->y, NULL, NULL);
			siril_plot_set_nth_color(spl_data, i+1, (double[3]){(double) (i == 2), (double) ((i == 3) + ((i == 1) * 0.5)), (double) ((i == 4) + ((i == 1) * 0.5)) });
			g_free(spl_legend);
			spcc_object_free_arrays(structs[i]);
			i++;
		}
		if (args.nb_mode) {
			double nbx[14] = {400.0,
				args.nb_center[0] - 0.5*args.nb_bandwidth[0], args.nb_center[0] - 0.5*args.nb_bandwidth[0],
				args.nb_center[0] + 0.5*args.nb_bandwidth[0], args.nb_center[0] + 0.5*args.nb_bandwidth[0],
				args.nb_center[1] - 0.5*args.nb_bandwidth[1], args.nb_center[1] - 0.5*args.nb_bandwidth[1],
				args.nb_center[1] + 0.5*args.nb_bandwidth[1], args.nb_center[1] + 0.5*args.nb_bandwidth[1],
				args.nb_center[2] - 0.5*args.nb_bandwidth[2], args.nb_center[2] - 0.5*args.nb_bandwidth[2],
				args.nb_center[2] + 0.5*args.nb_bandwidth[2], args.nb_center[2] + 0.5*args.nb_bandwidth[2], 1000.0};
			double nby[14] = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
			siril_plot_add_xydata(spl_data, _("Narrowband filters"), 14, nbx, nby, NULL, NULL);
			siril_plot_set_nth_color(spl_data, i+1, (double[3]){0.8, 0.8, 0.0});
		}
	} else {
		spcc_object *sensor_r = NULL, *sensor_g = NULL, *sensor_b = NULL, *filter_osc = NULL, *filter_lpf = NULL, *whiteref = NULL;
		if (args.selected_sensor_osc >= 0 && args.selected_sensor_osc < g_list_length (sensor_list)) {
			osc_sensor *osc = (osc_sensor*) g_list_nth(sensor_list, args.selected_sensor_osc)->data;
			sensor_r = &osc->channel[RLAYER];
			sensor_g = &osc->channel[GLAYER];
			sensor_b = &osc->channel[BLAYER];
		}
		if (!args.nb_mode) {
			if (args.selected_filter_osc >= 0 && args.selected_filter_osc < g_list_length (filter_list_osc))
				filter_osc = (spcc_object*) g_list_nth(filter_list_osc, args.selected_filter_osc)->data;
		}
		if (args.selected_filter_lpf >= 0 && args.selected_filter_lpf < g_list_length (filter_list_lpf))
			filter_lpf = (spcc_object*) g_list_nth(filter_list_lpf, args.selected_filter_lpf)->data;
		if (args.selected_white_ref >= 0 && args.selected_white_ref < g_list_length (whiteref_list))
			whiteref = (spcc_object*) g_list_nth(whiteref_list, args.selected_white_ref)->data;
		// there must not be any NULL spcc_object*s with non-NULL ones after them
		spcc_object* structs[7] = { whiteref, sensor_r, sensor_g, sensor_b, filter_osc, filter_lpf, NULL };
		while (structs[i]) {
			load_spcc_object_arrays(structs[i]);
			if (structs[i] == whiteref)
				normalize_y(structs[i], 400, 700);
			gchar *spl_legend = g_strdup(structs[i]->name);
			siril_plot_add_xydata(spl_data, spl_legend, structs[i]->n, structs[i]->x, structs[i]->y, NULL, NULL);
			siril_plot_set_nth_color(spl_data, i+1, (double[3]){(double) (i == 1 || i == 4), (double) ((i == 2) + ((i == 5) * 0.5)), (double) ((i == 3 || i == 4) + ((i == 5) * 0.5)) });
			g_free(spl_legend);
			spcc_object_free_arrays(structs[i]);
			i++;
		}
		if (args.nb_mode) {
			double nbx[14] = {400.0,
				args.nb_center[0] - 0.5*args.nb_bandwidth[0], args.nb_center[0] - 0.5*args.nb_bandwidth[0],
				args.nb_center[0] + 0.5*args.nb_bandwidth[0], args.nb_center[0] + 0.5*args.nb_bandwidth[0],
				args.nb_center[1] - 0.5*args.nb_bandwidth[1], args.nb_center[1] - 0.5*args.nb_bandwidth[1],
				args.nb_center[1] + 0.5*args.nb_bandwidth[1], args.nb_center[1] + 0.5*args.nb_bandwidth[1],
				args.nb_center[2] - 0.5*args.nb_bandwidth[2], args.nb_center[2] - 0.5*args.nb_bandwidth[2],
				args.nb_center[2] + 0.5*args.nb_bandwidth[2], args.nb_center[2] + 0.5*args.nb_bandwidth[2], 1000.0};
			double nby[14] = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
			siril_plot_add_xydata(spl_data, _("Narrowband filters"), 14, nbx, nby, NULL, NULL);
			siril_plot_set_nth_color(spl_data, i+1, (double[3]){0.8, 0.8, 0.0});
		}
	}
	if (args.atmos_corr) {
		if (!args.fit) args.fit = gfit;
		xpsampled atmos = init_xpsampled();
		fill_xpsampled_from_atmos_model(&atmos, &args);
		gchar *spl_legend = g_strdup(_("Atmosphere model"));
		siril_plot_add_xydata(spl_data, spl_legend, XPSAMPLED_LEN, atmos.x, atmos.y, NULL, NULL);
		siril_plot_set_nth_color(spl_data, i, (double[3]){ 0.67, 0.67, 1.0 } );
		g_free(spl_legend);
	}
	spl_data->datamin.x = MIN_PLOT;
	spl_data->datamax.x = MAX_PLOT;
	spl_data->datamax.y = min(spl_data->datamax.y, 1.2);
	spl_data->cfgdata.line.sz = 2;
	siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}

void on_spcc_details_plot_clicked(GtkButton *button, gpointer user_data) {
	gboolean is_osc_sensor = (cbdata->type == OSC_SENSORS);
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data)
		return;
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
		siril_plot_set_ylabel(spl_data, (object->type == MONO_SENSORS || object->type == OSC_SENSORS) ? _("Quantum Efficiency") : object->type == WB_REFS ? _("Relative photon count") : _("Transmittance"));

		siril_plot_add_xydata(spl_data, object->name, object->n, object->x, object->y, NULL, NULL);
		spcc_object_free_arrays(object);
	}
	spl_data->datamin.x = MIN_PLOT;
	if (cbdata->type != WB_REFS)
		spl_data->datamax.x = MAX_PLOT;
	spl_data->cfgdata.line.sz = 2;
	siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}

void on_spcc_details_close_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("spcc_details");
	gtk_widget_hide(win);
}

void on_spcc_combo_changed(GtkComboBox *combo, gpointer user_data) {
	if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_osc"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.oscsensorpref) {
			g_free(com.pref.spcc.oscsensorpref);
			com.pref.spcc.oscsensorpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.osc_sensors, index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.oscsensorpref = g_strdup(object->model);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_sensors_mono"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.monosensorpref) {
			g_free(com.pref.spcc.monosensorpref);
			com.pref.spcc.monosensorpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.mono_sensors, index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.monosensorpref = g_strdup(object->name);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_osc"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.oscfilterpref) {
			g_free(com.pref.spcc.oscfilterpref);
			com.pref.spcc.oscfilterpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.osc_filters, index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.oscfilterpref = g_strdup(object->name);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_r"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.redpref) {
			g_free(com.pref.spcc.redpref);
			com.pref.spcc.redpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.mono_filters[RLAYER], index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.redpref = g_strdup(object->name);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_g"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.greenpref) {
			g_free(com.pref.spcc.greenpref);
			com.pref.spcc.greenpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.mono_filters[GLAYER], index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.greenpref = g_strdup(object->name);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_b"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.bluepref) {
			g_free(com.pref.spcc.bluepref);
			com.pref.spcc.bluepref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.mono_filters[BLAYER], index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.bluepref = g_strdup(object->name);
		}
	} else if (combo == GTK_COMBO_BOX(lookup_widget("combo_spcc_filters_lpf"))) {
		int index = gtk_combo_box_get_active(combo);
		if (com.pref.spcc.lpfpref) {
			g_free(com.pref.spcc.lpfpref);
			com.pref.spcc.lpfpref = NULL;
		}
		GList *list = g_list_nth(com.spcc_data.osc_lpf, index);
		if (list) {
			spcc_object* object = (spcc_object*) list->data;
			com.pref.spcc.lpfpref = g_strdup(object->name);
		}
	}
	set_osc_lpf_visibility();
}

void on_spcc_toggle_nb_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean state = gtk_toggle_button_get_active(button);
	com.pref.spcc.nb_mode = state;
	gtk_widget_set_visible(lookup_widget("spcc_nb_controls"), state);
	gtk_widget_set_visible(lookup_widget("label_spcc_filters_r"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("combo_spcc_filters_r"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("details_spcc_filters_r"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("label_spcc_filters_g"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("combo_spcc_filters_g"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("details_spcc_filters_g"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("label_spcc_filters_b"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("combo_spcc_filters_b"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("details_spcc_filters_b"), !state && com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("label_spcc_filters_osc"), !state && !com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("combo_spcc_filters_osc"), !state && !com.pref.spcc.is_mono);
	gtk_widget_set_visible(lookup_widget("details_spcc_filters_osc"), !state && !com.pref.spcc.is_mono);
}

void on_nb_spin_changed(GtkSpinButton *button, gpointer user_data) {
	if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_r_wl"))) {
		com.pref.spcc.red_wl = gtk_spin_button_get_value(button);
	} else if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_g_wl"))) {
		com.pref.spcc.green_wl = gtk_spin_button_get_value(button);
	} else if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_b_wl"))) {
		com.pref.spcc.blue_wl = gtk_spin_button_get_value(button);
	} else if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_r_bw"))) {
		com.pref.spcc.red_bw = gtk_spin_button_get_value(button);
	} else if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_g_bw"))) {
		com.pref.spcc.green_bw = gtk_spin_button_get_value(button);
	} else if (button == GTK_SPIN_BUTTON(lookup_widget("spcc_nb_b_bw"))) {
		com.pref.spcc.blue_bw = gtk_spin_button_get_value(button);
	}
}

void on_button_gaia_status_check_clicked(GtkButton *button, gpointer user_data) {
	check_gaia_archive_status();
}

void on_spcc_atmos_corr_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_expander_set_expanded(GTK_EXPANDER(lookup_widget("spcc_atmos_widgets")), gtk_toggle_button_get_active(button));
}

void on_spcc_plot_atmos_clicked(GtkButton* button, gpointer user_data) {
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data)
		return;

	siril_plot_set_xlabel(spl_data, _("Wavelength / nm"));
	siril_plot_set_savename(spl_data, "SPCC_data");
	gchar *title = g_strdup_printf(_("SPCC Data: Atmospheric model"));
	siril_plot_set_title(spl_data, title);
	g_free(title);
	siril_plot_set_ylabel(spl_data, _("Transmittance"));

	// Generate data
	xpsampled data = init_xpsampled();
	struct photometric_cc_data args = { 0 };
	set_spcc_args(&args);
	args.fit = gfit;
	fill_xpsampled_from_atmos_model(&data, &args);
	gchar *spl_legend = g_strdup_printf(_("Atmospheric model\nHeight: %d m\nPressure: %.2f hPa (%s)"), (int) args.atmos_obs_height, args.atmos_pressure, args.atmos_pressure_is_slp ? N_("sea level") : N_("local"));
	siril_plot_add_xydata(spl_data, spl_legend, XPSAMPLED_LEN, data.x, data.y, NULL, NULL);
	g_free(spl_legend);
	spl_data->datamin.x = MIN_PLOT;
	spl_data->datamax.x = MAX_PLOT;
	spl_data->cfgdata.line.sz = 2;
	siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}
