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

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/siril_app_dirs.h"
#include "core/siril_language.h"
#include "core/settings.h"
#include "core/siril_log.h"
#include "algos/astrometry_solver.h"
#include "io/annotation_catalogues.h"
#include "io/siril_pythonmodule.h"
#include "gui/annotations_pref.h"
#include "gui/callbacks.h"
#include "gui/icc_profile.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/script_menu.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"
#include "gui/photometric_cc.h"
#include "gui/python_gui.h"
#include "gui/registration.h"
#include "gui/siril_intro.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "stacking/stacking.h"
#include "io/siril_git.h"

#include "preferences.h"
#include "filters/starnet.h"

#ifndef W_OK
#define W_OK 2
#endif

static gchar *sw_dir = NULL;
static gchar *st_weights = NULL;
static starnet_version st_version = NIL;
static gboolean update_custom_gamut = FALSE;
void on_working_gamut_changed(GtkComboBox *combo, gpointer user_data);

static gboolean scripts_updated = FALSE;

void notify_script_update() {
	scripts_updated = TRUE;
}

static void reset_swapdir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	const gchar *dir;

	dir = g_get_tmp_dir();

	if (g_strcmp0(dir, com.pref.swap_dir)) {
		g_free(com.pref.swap_dir);
		com.pref.swap_dir = g_strdup(dir);
		gtk_file_chooser_set_filename(swap_dir, dir);
	}
}

static void update_debayer_preferences() {
	com.pref.debayer.use_bayer_header = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_use_header")));
	com.pref.debayer.orientation = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_roworder")));
	com.pref.debayer.xbayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")));
	com.pref.debayer.ybayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")));
	com.pref.debayer.bayer_pattern = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")));
	com.pref.debayer.bayer_inter = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")));
	com.pref.debayer.xtrans_passes = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xtranspass_spin")));
}

static void update_astrometry_preferences() {
	get_astrometry_catalogue_values();

	com.pref.gui.position_compass = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("compass_combobox")));
	com.pref.wcs_formalism = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("wcs_formalism_combobox")));
	gchar *newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path1")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[0]);
		com.pref.catalogue_paths[0] = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path2")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[1]);
		com.pref.catalogue_paths[1] = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path3")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[2]);
		com.pref.catalogue_paths[2] = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path4")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[3]);
		com.pref.catalogue_paths[3] = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path5")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[4]);
		com.pref.catalogue_paths[4] = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path6")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.catalogue_paths[5]);
		com.pref.catalogue_paths[5] = newpath;
	}

	com.pref.astrometry.sip_correction_order = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_sip_order")));
	com.pref.astrometry.percent_scale_range = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("spin_sampling_tolerance")));
	com.pref.astrometry.radius_degrees = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_radius")));
	com.pref.astrometry.keep_xyls_files = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_xyls")));
	com.pref.astrometry.keep_wcs_files = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_wcs")));
	com.pref.astrometry.max_seconds_run = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_max_sec")));
	com.pref.astrometry.update_default_scale = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("astrometry_update_fields")));
	com.pref.astrometry.show_asnet_output = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_show_output")));
	if (com.pref.astrometry.default_obscode)
		g_free(com.pref.astrometry.default_obscode);
	com.pref.astrometry.default_obscode = g_strdup(gtk_entry_get_text(GTK_ENTRY(lookup_widget("obscode_entry"))));
	if (strlen(com.pref.astrometry.default_obscode) != 3 && strlen(com.pref.astrometry.default_obscode) != 0) {
		g_free(com.pref.astrometry.default_obscode);
		com.pref.astrometry.default_obscode = NULL;
		siril_log_color_message(_("Error: invalid IAU observatory code read from preferences file. Code must be a 3-character code.\n"), "red");
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("obscode_entry")), "");
	}
	if (com.pref.astrometry.default_obscode && strlen(com.pref.astrometry.default_obscode) == 0) {
		g_free(com.pref.astrometry.default_obscode);
		com.pref.astrometry.default_obscode = NULL;
	}
	// In the prefs structure, the dir is stored alongside starnet, not in astrometry
	com.pref.asnet_dir = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_asnet")));
	reset_astrometry_checks();
}

static void update_prepro_preferences() {

	if (com.pref.prepro.bias_lib) {
		g_free(com.pref.prepro.bias_lib);
		com.pref.prepro.bias_lib = NULL;
	}
	const gchar *biasentry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("biaslib_entry")));
	if (biasentry) {
		com.pref.prepro.bias_lib = g_strdup(biasentry);
		com.pref.prepro.use_bias_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias")));
	}

	if (com.pref.prepro.dark_lib) {
		g_free(com.pref.prepro.dark_lib);
		com.pref.prepro.dark_lib = NULL;
	}
	const gchar *darkentry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("darklib_entry")));
	if (darkentry) {
		com.pref.prepro.dark_lib = g_strdup(darkentry);
		com.pref.prepro.use_dark_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark")));
	}

	if (com.pref.prepro.flat_lib) {
		g_free(com.pref.prepro.flat_lib);
		com.pref.prepro.flat_lib = NULL;
	}
	const gchar *flatentry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("flatlib_entry")));
	if (flatentry) {
		com.pref.prepro.flat_lib = g_strdup(flatentry);
		com.pref.prepro.use_flat_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat")));
	}

	if (com.pref.prepro.disto_lib) {
		g_free(com.pref.prepro.disto_lib);
		com.pref.prepro.disto_lib = NULL;
	}
	const gchar *distoentry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("distolib_entry")));
	if (distoentry) {
		com.pref.prepro.disto_lib = g_strdup(distoentry);
		com.pref.prepro.use_disto_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_disto")));
	}

	if (com.pref.prepro.stack_default) {
		g_free(com.pref.prepro.stack_default);
		com.pref.prepro.stack_default = NULL;
	}
	const gchar *stackentry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("stack_default_entry")));
	if (stackentry && stackentry[0] != '\0') {
		com.pref.prepro.stack_default = g_strdup(stackentry);
		com.pref.prepro.use_stack_default = (com.pref.prepro.stack_default) && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_stack")));
	}

	com.pref.prepro.xtrans_af.x = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_x"))), NULL, 10);
	com.pref.prepro.xtrans_af.y = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_y"))), NULL, 10);
	com.pref.prepro.xtrans_af.w = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_w"))), NULL, 10);
	com.pref.prepro.xtrans_af.h = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_h"))), NULL, 10);

	com.pref.prepro.xtrans_sample.x = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_x"))), NULL, 10);
	com.pref.prepro.xtrans_sample.y = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_y"))), NULL, 10);
	com.pref.prepro.xtrans_sample.w = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_w"))), NULL, 10);
	com.pref.prepro.xtrans_sample.h = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_h"))), NULL, 10);
}

static void update_photometry_preferences() {
	com.pref.phot_set.gain = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")));
	com.pref.phot_set.inner = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")));
	com.pref.phot_set.outer = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")));
	com.pref.phot_set.aperture = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinAperture")));
	com.pref.phot_set.force_radius = !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("photometry_force_radius_button")));
	com.pref.phot_set.auto_aperture_factor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinRadRatio")));
	com.pref.phot_set.minval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")));
	com.pref.phot_set.maxval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")));
}

static void update_analysis_preferences() {
	com.pref.analysis.mosaic_panel = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinpanel")));
	com.pref.analysis.mosaic_window = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinwindow")));
}

static void update_scripts_preferences() {
	if (com.pref.gui.script_path)
		g_slist_free_full(g_steal_pointer(&com.pref.gui.script_path), g_free);
	com.pref.gui.script_path = get_list_from_preferences_dialog();
	com.pref.gui.warn_scripts_run = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")));
	com.pref.script_check_requires = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")));
#ifdef HAVE_LIBGIT2
	com.pref.use_scripts_repository = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_use_gitscripts")));
	com.pref.auto_script_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_script_automatic_updates")));
	com.pref.spcc.use_spcc_repository = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("spcc_repo_enable")));
	com.pref.spcc.auto_spcc_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("spcc_repo_sync_at_startup")));
#else
	com.pref.use_scripts_repository = FALSE;
	com.pref.spcc.use_spcc_repository = FALSE;
#endif
}

static void update_user_interface_preferences() {
	com.pref.lang = get_interface_language();
	int theme = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_theme")));
	if (theme != com.pref.gui.combo_theme) {
		com.pref.gui.combo_theme = theme;
		siril_set_theme(theme);
	}
	com.pref.gui.font_scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("pref_fontsize")));
	com.pref.gui.icon_symbolic = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_iconstyle")));
	com.pref.gui.remember_windows = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck")));
	com.pref.gui.show_thumbnails = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button")));
	com.pref.gui.thumbnail_size = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("thumbnails_box_size"))) == 1 ? 256 : 128;
	com.pref.gui.default_rendering_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_default_stf")));
	com.pref.gui.display_histogram_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_default_histo_mode")));
	com.pref.gui.roi_mode = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_ui_roimode")));
	com.pref.gui.mmb_action = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_ui_mmb_zoom")));
	com.pref.gui.mouse_speed_limit = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_mouse_speed_limit")));
	update_roi_config();
	/* Configure colors */
	GdkRGBA color;

	gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_bkg")), &color);
	g_free(com.pref.gui.config_colors.color_bkg_samples);
	com.pref.gui.config_colors.color_bkg_samples = gdk_rgba_to_string(&color);

	gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_std_annot")), &color);
	g_free(com.pref.gui.config_colors.color_std_annotations);
	com.pref.gui.config_colors.color_std_annotations = gdk_rgba_to_string(&color);

	gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_dso_annot")), &color);
	g_free(com.pref.gui.config_colors.color_dso_annotations);
	com.pref.gui.config_colors.color_dso_annotations = gdk_rgba_to_string(&color);

	gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_sso_annot")), &color);
	g_free(com.pref.gui.config_colors.color_sso_annotations);
	com.pref.gui.config_colors.color_sso_annotations = gdk_rgba_to_string(&color);

	gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_tmp_annot")), &color);
	g_free(com.pref.gui.config_colors.color_tmp_annotations);
	com.pref.gui.config_colors.color_tmp_annotations = gdk_rgba_to_string(&color);
}

static void update_color_management_preferences() {
	gchar *newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("pref_custom_monitor_profile")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.icc.icc_path_monitor);
		com.pref.icc.icc_path_monitor = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("pref_soft_proofing_profile")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.icc.icc_path_soft_proof);
		com.pref.icc.icc_path_soft_proof = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("custom_gray_icc_matching_trc")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.icc.custom_icc_gray);
		com.pref.icc.custom_icc_gray = newpath;
	}
	newpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("custom_icc_standard_trc")));
	if (newpath && newpath[0] != '\0') {
		g_free(com.pref.icc.custom_icc_trc);
		com.pref.icc.custom_icc_trc = newpath;
	}
	com.pref.icc.rendering_intent = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_rendering_intent")));
	com.pref.icc.export_intent = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_export_intent")));
	com.pref.icc.default_to_srgb = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("default_to_srgb")));
	com.pref.icc.working_gamut = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("working_gamut")));
	com.pref.icc.export_8bit_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("export_profile_8bit")));
	com.pref.icc.export_16bit_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("export_profile_16bit")));
	com.pref.icc.cmf = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_cmf_pref")));
	com.pref.icc.custom_monitor_profile_active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("custom_monitor_profile_active")));
	com.pref.icc.soft_proofing_profile_active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("custom_proofing_profile_active")));
	com.pref.icc.rendering_bpc = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_rendering_bpc")));
	com.pref.icc.pedantic_linear = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_pedantic_linear")));
	com.icc.rendering_flags = ((com.pref.icc.rendering_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.rendering_intent == INTENT_ABSOLUTE_COLORIMETRIC));
	gui.icc.proofing_flags = ((com.pref.icc.rendering_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.rendering_intent == INTENT_ABSOLUTE_COLORIMETRIC)) | cmsFLAGS_SOFTPROOFING;
	com.pref.icc.autoassignment = 	((gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_load"))) * ICC_ASSIGN_ON_LOAD) +
									(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_stack"))) * ICC_ASSIGN_ON_STACK) +
									(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_stretch"))) * ICC_ASSIGN_ON_STRETCH) +
									(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_composition"))) * ICC_ASSIGN_ON_COMPOSITION));
	com.pref.icc.autoconversion = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_icc_autoconversion")));
}

static void update_FITS_options_preferences() {
	com.pref.comp.fits_enabled = !(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio"))));
	com.pref.comp.fits_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")));
	com.pref.comp.fits_quantization = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")));
	com.pref.comp.fits_hcompress_scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")));

	com.pref.rgb_aladin = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_aladin")));
	com.pref.use_checksum = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_fits_chksum")));
	com.pref.binning_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_binned_update")));

	const gchar *ext = gtk_combo_box_get_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")));
	com.pref.ext = g_strdup(ext);

	com.pref.force_16bit = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_type"))) == 0;
	com.pref.fits_save_icc = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_fits_save_icc")));

	com.pref.allow_heterogeneous_fitseq = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbox_heterogeneous_fitseq")));
}

static void update_performances_preferences() {
	GSList *tmp_list = gtk_radio_button_get_group (GTK_RADIO_BUTTON(lookup_widget("memfreeratio_radio")));
	GtkWidget *amount;
	GtkToggleButton *tmp_button = NULL;//Create a temp toggle button.

	amount = lookup_widget("memfixed_radio");
	while (tmp_list) {
		tmp_button = tmp_list->data;
		tmp_list = tmp_list->next;

		if (gtk_toggle_button_get_active(tmp_button))
			break;

		tmp_button = NULL; //We've enumerated all of them, and none of them is active.
	}
	if (tmp_button) {
		if (GTK_WIDGET(tmp_button) == amount) {
			com.pref.mem_mode = AMOUNT;
		} else {
			com.pref.mem_mode = RATIO;
		}
	}

	com.pref.memory_ratio = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio")));
	com.pref.memory_amount = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount")));
	com.pref.fftw_conf.timelimit = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("pref_fftw_plan_timelimit")));
	com.pref.fftw_conf.strategy = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_fftw_plan_strategy")));
	com.pref.fftw_conf.multithreaded = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_fftw_multithreaded")));
	com.pref.fftw_conf.fft_cutoff = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("pref_conv_min_fft")));
	int max_slice_size = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("pref_max_slice_size")));
	com.pref.max_slice_size = max_slice_size == 0 ? 32769 : 1 << (max_slice_size + 8);
	int bitdepth = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_hd_bitdepth")));
	com.pref.hd_bitdepth = bitdepth;
}

static void update_misc_preferences() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	GtkFileChooser *starnet_exe = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet"));
	GtkFileChooser *starnet_weights = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet_weights"));
	GtkFileChooser *graxpert_exe = GTK_FILE_CHOOSER(lookup_widget("filechooser_graxpert"));

	com.pref.swap_dir = gtk_file_chooser_get_filename(swap_dir);

	com.pref.starnet_exe = gtk_file_chooser_get_filename(starnet_exe);
	com.pref.starnet_weights = gtk_file_chooser_get_filename(starnet_weights);
	com.pref.graxpert_path = gtk_file_chooser_get_filename(graxpert_exe);
#ifdef OS_OSX
	if (g_str_has_suffix(com.pref.graxpert_path, ".app"))
		str_append(&com.pref.graxpert_path, "/Contents/MacOS/GraXpert");
#endif

	com.pref.gui.silent_quit = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")));
	com.pref.gui.silent_linear = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskSave")));

	const gchar *copy = gtk_entry_get_text(GTK_ENTRY(lookup_widget("miscCopyright")));
	com.pref.copyright = g_strdup(copy);

	com.pref.check_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")));
	com.pref.gui.enable_roi_warning = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscHideInfoROI"))) ? FALSE : TRUE;
}

void on_checkbutton_use_header_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = !gtk_toggle_button_get_active(button);
	GtkWidget *combo = lookup_widget("comboBayer_pattern");
	GtkWidget *spin1 = lookup_widget("xbayeroff_spin");
	GtkWidget *spin2 = lookup_widget("ybayeroff_spin");

	gtk_widget_set_sensitive(combo, active);
	gtk_widget_set_sensitive(spin1, active);
	gtk_widget_set_sensitive(spin2, active);
}

void on_photometry_force_radius_button_toggled(GtkToggleButton *button, gpointer user_data) {
//	GtkWidget *spin = (GtkWidget *)user_data;
	GtkWidget *spin1 = lookup_widget("spinAperture");
	GtkWidget *spin2 = lookup_widget("spinRadRatio");
	gtk_widget_set_sensitive(spin1, !gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(spin2, gtk_toggle_button_get_active(button));
}

void initialize_path_directory(const gchar *path) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (swap_dir, path);
	} else {
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
	}
}

void initialize_graxpert_executable(gchar *path) {
	GtkFileChooser *graxpert_exe = GTK_FILE_CHOOSER(lookup_widget("filechooser_graxpert"));
#ifdef OS_OSX
	const gchar *suffix = "/Contents/MacOS/GraXpert";
	if (path && g_str_has_suffix(path, suffix)) {
		gsize len = strlen(path) - strlen(suffix);
		gchar *new_path = g_strndup(path, len);
		gtk_file_chooser_set_filename(graxpert_exe, new_path);
		g_free(new_path);
		return;
	}
#endif

	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (graxpert_exe, path);
	}
}

void initialize_starnet_executable(gchar *path) {
#ifdef HAVE_LIBTIFF
	GtkFileChooser *starnet_exe = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet"));
	GtkWidget *starnet_weights_reset = GTK_WIDGET(lookup_widget("starnet_weights_clear"));
	GtkWidget *starnet_weights = GTK_WIDGET(lookup_widget("filechooser_starnet_weights"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (starnet_exe, path);
		if (starnet_executablecheck(path) & TORCH) {
			gtk_widget_set_sensitive(starnet_weights, TRUE);
			gtk_widget_set_sensitive(starnet_weights_reset, TRUE);
		} else {
			gtk_widget_set_sensitive(starnet_weights, FALSE);
			gtk_widget_set_sensitive(starnet_weights_reset, FALSE);
		}
	} else {
		gtk_widget_set_sensitive(starnet_weights, FALSE);
		gtk_widget_set_sensitive(starnet_weights_reset, FALSE);
	}
#endif
}

void initialize_starnet_weights(gchar *path) {
	GtkFileChooser *starnet_weights = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet_weights"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (starnet_weights, path);
	}
}

void initialize_asnet_directory(const gchar *path) {
	GtkFileChooser *asnet_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_asnet"));
	reset_asnet_version();
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (asnet_dir, path);
	}
}

void on_asnet_clear_clicked(GtkButton *button, gpointer user_data) {
	GtkFileChooser *asnet_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_asnet"));
	gtk_file_chooser_set_filename(asnet_dir, "\0");
	reset_asnet_version();
}

void on_button_reset_swap_clicked(GtkButton *button, gpointer user_data) {
	reset_swapdir();
}

void on_combobox_comp_fits_method_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *hcompress_scale_spin = lookup_widget("spinbutton_comp_fits_hcompress_scale");
	GtkSpinButton *button = (GtkSpinButton *)user_data;
	gint method = gtk_combo_box_get_active(GTK_COMBO_BOX(box));
	if (gtk_spin_button_get_value(button) == 0.0) {
		gtk_spin_button_set_value(button, 16.0);
	}
	gtk_widget_set_sensitive(hcompress_scale_spin, (method == HCOMPRESS_COMP) ? TRUE : FALSE);
}

void on_mem_radio_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")),
			*amount = GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio"));
	GtkWidget *ratio_spin = lookup_widget("spinbutton_mem_ratio"),
		  *amount_spin = lookup_widget("spinbutton_mem_amount");
	if (!gtk_toggle_button_get_active(togglebutton)) return;

	gtk_widget_set_sensitive(ratio_spin, togglebutton == ratio);
	gtk_widget_set_sensitive(amount_spin, togglebutton == amount);
}

void on_comp_fits_radio_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *disabled = GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio"));
	GtkWidget *method_box = lookup_widget("combobox_comp_fits_method"),
		*quantization_spin = lookup_widget("spinbutton_comp_fits_quantization"),
		*tilex_spin = lookup_widget("spinbutton_comp_fits_tileX"),
		*tiley_spin = lookup_widget("spinbutton_comp_fits_tileY"),
		*hcompress_scale_spin = lookup_widget("spinbutton_comp_fits_hcompress_scale");
	if (!gtk_toggle_button_get_active(togglebutton)) return;

	gint method = gtk_combo_box_get_active(GTK_COMBO_BOX(method_box));
	gtk_widget_set_sensitive(method_box, togglebutton != disabled);
	gtk_widget_set_sensitive(quantization_spin, togglebutton != disabled);
	gtk_widget_set_sensitive(tilex_spin, FALSE);
	gtk_widget_set_sensitive(tiley_spin, FALSE);
	gtk_widget_set_sensitive(hcompress_scale_spin, (method == HCOMPRESS_COMP) ? togglebutton != disabled : FALSE);
}

void on_filechooser_swap_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(fileChooser);
	gchar *dir;

	dir = gtk_file_chooser_get_filename (swap_dir);

	if (g_access(dir, W_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to write in this directory: %s\n"), "red", dir);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(swap_dir, com.pref.swap_dir);
		return;
	}
	if (sw_dir) {
		g_free(sw_dir);
		sw_dir = gtk_file_chooser_get_filename(swap_dir);
	}
}

void on_filechooser_starnet_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
#ifdef HAVE_LIBTIFF
	GtkFileChooser *starnet_exe = GTK_FILE_CHOOSER(fileChooser);
	gchar *path;

	path = gtk_file_chooser_get_filename (starnet_exe);
	printf("%s\n", path);
	if (g_access(path, X_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to execute this file: %s\n"), "red", path);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(starnet_exe, com.pref.starnet_exe);
		return;
	}
	st_version = starnet_executablecheck(path);
	GtkWidget *starnet_weights_reset = GTK_WIDGET(lookup_widget("starnet_weights_clear"));
	GtkWidget *starnet_weights = GTK_WIDGET(lookup_widget("filechooser_starnet_weights"));
	if (st_version & TORCH) {
		gtk_widget_set_sensitive(starnet_weights, TRUE);
		gtk_widget_set_sensitive(starnet_weights_reset, TRUE);
	} else {
		gtk_widget_set_sensitive(starnet_weights, FALSE);
		gtk_widget_set_sensitive(starnet_weights_reset, FALSE);
	}
#endif
}

void on_starnet_weights_clear_clicked(GtkButton *button, gpointer user_data) {
	GtkFileChooser *starnet_weights = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet_weights"));
	if (st_weights) {
		g_free(st_weights);
	}
	st_weights = g_strdup("\0");
	gtk_file_chooser_set_filename(starnet_weights, st_weights);
}

void on_filechooser_starnet_weights_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
	GtkFileChooser *starnet_weights = GTK_FILE_CHOOSER(fileChooser);
	gchar *path;

	path = gtk_file_chooser_get_filename (starnet_weights);

	if (g_access(path, R_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to read this file: %s\n"), "red", path);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(starnet_weights, com.pref.starnet_weights);
		return;
	}
	if (st_weights) {
		g_free(st_weights);
	}
	st_weights = gtk_file_chooser_get_filename(starnet_weights);
}

void on_show_preview_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *label = lookup_widget("thumbnails_label_size");
	GtkWidget *box = lookup_widget("thumbnails_box_size");

	gtk_widget_set_sensitive(label, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(box, gtk_toggle_button_get_active(togglebutton));
}

void on_play_introduction_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("settings_window");
	start_intro_script();
}

void on_reload_script_button_clicked(GtkButton *button, gpointer user_data) {
	g_thread_unref(g_thread_new("refresh_scripts", refresh_scripts_in_thread, NULL));
}

void on_check_button_pref_bias_bis_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *bias_button = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias"));

	if (gtk_toggle_button_get_active(bias_button)) {
		g_signal_handlers_block_by_func(bias_button, on_check_button_pref_bias_toggled, NULL);
		gtk_toggle_button_set_active(bias_button, FALSE);
		g_signal_handlers_unblock_by_func(bias_button, on_check_button_pref_bias_toggled, NULL);
		gtk_toggle_button_set_active(togglebutton, TRUE);
	}
}

void on_check_button_pref_bias_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *bias_button_bis = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias_bis"));

	if (gtk_toggle_button_get_active(bias_button_bis)) {
		g_signal_handlers_block_by_func(bias_button_bis, on_check_button_pref_bias_bis_toggled, NULL);
		gtk_toggle_button_set_active(bias_button_bis, FALSE);
		g_signal_handlers_unblock_by_func(bias_button_bis, on_check_button_pref_bias_bis_toggled, NULL);
		gtk_toggle_button_set_active(togglebutton, TRUE);
	}
}

//static gboolean from_prefs_init = FALSE;	// NOT USED, SHOULD BE DELETED

void update_preferences_from_model() {
	siril_debug_print("updating preferences GUI from settings data\n");
	preferences *pref = &com.pref;
	/* tab FITS/SER Debayer */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_use_header")), pref->debayer.use_bayer_header);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")), pref->debayer.bayer_pattern);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")), pref->debayer.bayer_inter);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")), pref->debayer.xbayeroff);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")), pref->debayer.ybayeroff);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_roworder")),  pref->debayer.orientation);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("xtranspass_spin")), pref->debayer.xtrans_passes);

	/* tab FITS Options */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_aladin")), pref->rgb_aladin);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_fits_chksum")), pref->use_checksum);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_binned_update")), pref->binning_update);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio")), !pref->comp.fits_enabled);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_enabled_radio")), pref->comp.fits_enabled);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")), pref->comp.fits_method);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")), pref->comp.fits_quantization);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")), pref->comp.fits_hcompress_scale);
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")), pref->ext == NULL ? ".fit" : pref->ext);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_type")), pref->force_16bit ? 0 : 1);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_fits_save_icc")), pref->fits_save_icc);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbox_heterogeneous_fitseq")), pref->allow_heterogeneous_fitseq);

	/* tab Astrometry */
	fill_astrometry_catalogue(pref->gui.catalog);

	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("compass_combobox")), pref->gui.position_compass);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("wcs_formalism_combobox")), pref->wcs_formalism);
	if (pref->catalogue_paths[0] && (g_file_test(pref->catalogue_paths[0], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path1"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[0]);
	}
	if (pref->catalogue_paths[1] && (g_file_test(pref->catalogue_paths[1], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path2"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[1]);
	}
	if (pref->catalogue_paths[2] && (g_file_test(pref->catalogue_paths[2], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path3"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[2]);
	}
	if (pref->catalogue_paths[3] && (g_file_test(pref->catalogue_paths[3], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path4"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[3]);
	}
	if (pref->catalogue_paths[4] && (g_file_test(pref->catalogue_paths[4], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path5"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[4]);
	}
	if (pref->catalogue_paths[5] && (g_file_test(pref->catalogue_paths[5], G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path6"));
		gtk_file_chooser_set_filename(button, pref->catalogue_paths[5]);
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_sip_order")), pref->astrometry.sip_correction_order);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_sampling_tolerance")), pref->astrometry.percent_scale_range);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_radius")), pref->astrometry.radius_degrees);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_xyls")), pref->astrometry.keep_xyls_files);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_wcs")), pref->astrometry.keep_wcs_files);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_asnet_max_sec")), pref->astrometry.max_seconds_run);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("astrometry_update_fields")), pref->astrometry.update_default_scale);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_asnet_show_output")), pref->astrometry.show_asnet_output);
	if (pref->astrometry.default_obscode)
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("obscode_entry")), pref->astrometry.default_obscode);

	/* tab Pre-processing */
	if (pref->prepro.bias_lib) {
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("biaslib_entry")), pref->prepro.bias_lib);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias")), pref->prepro.use_bias_lib);
	}

	if (pref->prepro.dark_lib) {
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("darklib_entry")), pref->prepro.dark_lib);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark")), pref->prepro.use_dark_lib);
	}

	if (pref->prepro.flat_lib) {
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("flatlib_entry")), pref->prepro.flat_lib);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat")), pref->prepro.use_flat_lib);
	}

	if (pref->prepro.disto_lib) {
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("distolib_entry")), pref->prepro.disto_lib);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_disto")), pref->prepro.use_disto_lib);
	}

	if (pref->prepro.stack_default) {
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("stack_default_entry")), pref->prepro.stack_default);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_stack")), pref->prepro.use_stack_default);
	} else {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_stack")), FALSE);
	}

	gchar tmp[256];
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_af.x);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_x")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_af.y);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_y")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_af.w);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_w")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_af.h);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_h")), tmp);

	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_sample.x);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_x")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_sample.y);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_y")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_sample.w);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_w")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->prepro.xtrans_sample.h);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_h")), tmp);

	/* tab Photometry */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")), pref->phot_set.outer);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")), pref->phot_set.inner);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("photometry_force_radius_button")), !pref->phot_set.force_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinRadRatio")), pref->phot_set.auto_aperture_factor);
	gtk_widget_set_sensitive(lookup_widget("spinRadRatio"), !pref->phot_set.force_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinAperture")), pref->phot_set.aperture);
	gtk_widget_set_sensitive(lookup_widget("spinAperture"), pref->phot_set.force_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")), pref->phot_set.gain);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")), pref->phot_set.minval);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")), pref->phot_set.maxval);

	/* tab Analysis Tools */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinpanel")), pref->analysis.mosaic_panel);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinwindow")), pref->analysis.mosaic_window);

	/* tab User Interface */
	siril_language_fill_combo(pref->lang);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_theme")), pref->gui.combo_theme);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("pref_fontsize")), pref->gui.font_scale);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_iconstyle")), pref->gui.icon_symbolic);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck")), pref->gui.remember_windows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button")), pref->gui.show_thumbnails);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("thumbnails_box_size")), pref->gui.thumbnail_size == 256 ? 1 : 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_default_stf")), pref->gui.default_rendering_mode);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_default_histo_mode")), pref->gui.display_histogram_mode);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_ui_roimode")), pref->gui.roi_mode);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_ui_mmb_zoom")), pref->gui.mmb_action);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_mouse_speed_limit")), pref->gui.mouse_speed_limit);

	GdkRGBA color;
	gdk_rgba_parse(&color, pref->gui.config_colors.color_bkg_samples);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_bkg")), &color);

	gdk_rgba_parse(&color, pref->gui.config_colors.color_std_annotations);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_std_annot")), &color);

	gdk_rgba_parse(&color, pref->gui.config_colors.color_dso_annotations);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_dso_annot")), &color);

	gdk_rgba_parse(&color, pref->gui.config_colors.color_sso_annotations);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_sso_annot")), &color);

	gdk_rgba_parse(&color, pref->gui.config_colors.color_tmp_annotations);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(lookup_widget("color_button_tmp_annot")), &color);

	/* tab Color Management */
	GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("pref_custom_monitor_profile"));
	if (pref->icc.icc_path_monitor && (g_file_test(pref->icc.icc_path_monitor, G_FILE_TEST_EXISTS))) {
		gtk_file_chooser_set_filename(button, pref->icc.icc_path_monitor);
	} else {
		gtk_file_chooser_unselect_all(button);
	}
	button = GTK_FILE_CHOOSER(lookup_widget("pref_soft_proofing_profile"));
	if (pref->icc.icc_path_soft_proof && (g_file_test(pref->icc.icc_path_soft_proof, G_FILE_TEST_EXISTS))) {
		gtk_file_chooser_set_filename(button, pref->icc.icc_path_soft_proof);
	} else {
		gtk_file_chooser_unselect_all(button);
	}
	button = GTK_FILE_CHOOSER(lookup_widget("custom_icc_standard_trc"));
	if (pref->icc.custom_icc_trc && (g_file_test(pref->icc.custom_icc_trc, G_FILE_TEST_EXISTS))) {
		gtk_file_chooser_set_filename(button, pref->icc.custom_icc_trc);
	} else {
		gtk_file_chooser_unselect_all(button);
	}
	button = GTK_FILE_CHOOSER(lookup_widget("custom_gray_icc_matching_trc"));
	if (pref->icc.custom_icc_gray && (g_file_test(pref->icc.custom_icc_gray, G_FILE_TEST_EXISTS))) {
		gtk_file_chooser_set_filename(button, pref->icc.custom_icc_gray);
	} else {
		gtk_file_chooser_unselect_all(button);
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_rendering_intent")), pref->icc.rendering_intent);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_export_intent")), pref->icc.export_intent);
	initialize_icc_preferences_widgets();
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("working_gamut")), pref->icc.working_gamut);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("export_profile_8bit")), pref->icc.export_8bit_method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("export_profile_16bit")), pref->icc.export_16bit_method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_cmf_pref")), pref->icc.cmf);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("custom_monitor_profile_active")), pref->icc.custom_monitor_profile_active);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("custom_proofing_profile_active")), pref->icc.soft_proofing_profile_active);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("default_to_srgb")), pref->icc.default_to_srgb);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_rendering_bpc")), pref->icc.rendering_bpc);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_pedantic_linear")), pref->icc.pedantic_linear);
	icc_assign_type autoassign = pref->icc.autoassignment;
	gboolean val = autoassign & ICC_ASSIGN_ON_LOAD;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_load")), val);
	val = autoassign & ICC_ASSIGN_ON_STACK;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_stack")), val);
	val = autoassign & ICC_ASSIGN_ON_STRETCH;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_stretch")), val);
	val = autoassign & ICC_ASSIGN_ON_COMPOSITION;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_on_composition")), val);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_icc_assign_never")), !pref->icc.autoassignment);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_icc_autoconversion")), pref->icc.autoconversion);

	/* tab Scripts */
	pref->gui.script_path = set_list_to_preferences_dialog(pref->gui.script_path);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")), pref->gui.warn_scripts_run);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")), pref->script_check_requires);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_use_gitscripts")), pref->use_scripts_repository);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_script_automatic_updates")), pref->auto_script_update);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("spcc_repo_enable")), pref->spcc.use_spcc_repository);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("spcc_repo_sync_at_startup")), pref->spcc.auto_spcc_update);

	/* tab Performances */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")), pref->mem_mode == RATIO);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")), pref->mem_mode == AMOUNT);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio")), pref->memory_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount")), pref->memory_amount);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_hd_bitdepth")), pref->hd_bitdepth);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("pref_fftw_plan_timelimit")), pref->fftw_conf.timelimit);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_fftw_plan_strategy")), pref->fftw_conf.strategy);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_fftw_multithreaded")), pref->fftw_conf.multithreaded);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("pref_conv_min_fft")), pref->fftw_conf.fft_cutoff);
	int n = pref->max_slice_size, max_slice_size = 0;
	// Shift n to the right until it becomes 1
	if (n > 0) {
		while (n > 1) {
			n >>= 1;
			max_slice_size++;
		}
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("pref_max_slice_size")), pref->max_slice_size > 32768 || pref->max_slice_size <= 0 ? 0 : max_slice_size - 8);

	/* tab Miscellaneous */
	initialize_path_directory(pref->swap_dir);
	initialize_graxpert_executable(pref->graxpert_path);
	initialize_starnet_executable(pref->starnet_exe);
	initialize_starnet_weights(pref->starnet_weights);
	initialize_asnet_directory(pref->asnet_dir);

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscHideInfoROI")), pref->gui.enable_roi_warning ? FALSE : TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")), pref->gui.silent_quit);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskSave")), pref->gui.silent_linear);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("miscCopyright")), pref->copyright == NULL ? "" : pref->copyright);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")), pref->check_update);
}

static void set_icc_filechooser_directories() {
	GtkFileChooser* fc2 = (GtkFileChooser*) lookup_widget("custom_icc_standard_trc");
	GtkFileChooser* fc3 = (GtkFileChooser*) lookup_widget("custom_gray_icc_matching_trc");
	GtkFileChooser* fc4 = (GtkFileChooser*) lookup_widget("pref_custom_monitor_profile");
	GtkFileChooser* fc5 = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");

	gtk_file_chooser_set_current_folder(fc2, default_system_icc_path());
	gtk_file_chooser_set_current_folder(fc3, default_system_icc_path());
	gtk_file_chooser_set_current_folder(fc4, default_system_icc_path());
	gtk_file_chooser_set_current_folder(fc5, default_system_icc_path());
	GtkComboBox* combo = (GtkComboBox*) lookup_widget("working_gamut");
	on_working_gamut_changed(combo, NULL);
}

static void dump_ui_to_global_var() {
	siril_debug_print("updating settings from preferences GUI\n");
	/* tab FITS/SER Debayer */
	update_debayer_preferences();
	/* tab FITS Options */
	update_FITS_options_preferences();
	/* tab Astrometry */
	update_astrometry_preferences();
	/* tab Pre-processing */
	update_prepro_preferences();
	/* tab Analysis Tools */
	update_analysis_preferences();
	/* tab Photometry */
	update_photometry_preferences();
	/* tab Scripts */
	update_scripts_preferences();
	/* tab User Interface */
	update_user_interface_preferences();
	/* tab Color Management */
	update_color_management_preferences();
	/* tab Performances */
	update_performances_preferences();
	/* tab Miscellaneous */
	update_misc_preferences();
}

void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data);
void on_spcc_repo_enable_toggled(GtkToggleButton *button, gpointer user_data);

void on_settings_window_show(GtkWidget *widget, gpointer user_data) {
	siril_debug_print("show preferences window: updating it\n");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path1")), "filter_namedstars", "Catalogue");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path2")), "filter_unnamedstars", "Catalogue");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path3")), "filter_deepstars", "Catalogue");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("localcatalogue_path4")), "filter_USNO-NOMAD-1e8", "Catalogue");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("custom_icc_standard_trc")), "icc_filter", "ICC files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("custom_gray_icc_matching_trc")), "icc_filter", "ICC files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("pref_custom_monitor_profile")), "icc_filter", "ICC files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("pref_soft_proofing_profile")), "icc_filter", "ICC files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("filechooser_graxpert")), "all_files", "All supported files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet")), "all_files", "All supported files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet_weights")), "all_files", "All supported files");
	GtkLabel* spcc_path_label = GTK_LABEL(lookup_widget("label_spcc_repo_path"));
	gtk_label_set_text(spcc_path_label, siril_get_spcc_repo_path());

	set_icc_filechooser_directories();
	g_signal_handlers_block_by_func(G_OBJECT(lookup_widget("spcc_repo_enable")), on_spcc_repo_enable_toggled, NULL);
	g_signal_handlers_block_by_func(G_OBJECT(lookup_widget("pref_use_gitscripts")), on_pref_use_gitscripts_toggled, NULL);
	update_preferences_from_model();
	g_signal_handlers_unblock_by_func(G_OBJECT(lookup_widget("spcc_repo_enable")), on_spcc_repo_enable_toggled, NULL);
	g_signal_handlers_unblock_by_func(G_OBJECT(lookup_widget("pref_use_gitscripts")), on_pref_use_gitscripts_toggled, NULL);
	scripts_updated = FALSE;
#ifndef HAVE_LIBGIT2
	hide_git_widgets();
#else
	fill_script_repo_tree(FALSE);
#endif
}

gboolean check_pref_sanity(gchar **error) {
	gchar *err = NULL;
	/* check for photometry pref */
	double in = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")));
	double out = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")));

	if (in >= out) {
		err = _("Inner radius value must be less than outer. Please change the value in Photometry tab.");
		*error = err;
		return FALSE;
	}
	/* check for FITS pref */
	gdouble quantization = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")));
	if (quantization == 0.0) {
		GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method"));
		if (gtk_combo_box_get_active(combo) != GZIP1_COMP && gtk_combo_box_get_active(combo) != GZIP2_COMP) {
			err = _("Setting quantization to 0 has only a sense with a GZIP compression "
									"and GZIP 2 often produces better compression of floating­point images.");
			*error = err;
			return FALSE;
		}
	}
	return TRUE;
}

void on_apply_settings_button_clicked(GtkButton *button, gpointer user_data) {
	gchar *err;
	if (!check_pref_sanity(&err)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong value"), err);
	} else {
		free_preferences(&com.pref);
		dump_ui_to_global_var();
#ifdef HAVE_LIBGIT2
		if (!com.pref.use_scripts_repository)
			on_disable_gitscripts();
		else
#endif
			g_thread_unref(g_thread_new("refresh_script_menu", refresh_script_menu_in_thread, GINT_TO_POINTER((int) scripts_updated)));
			// To update the UI with scripts from the repo
			// Note this line is part of the if/else with the #ifdef and always runs
			// otherwise. This is intentional.
		scripts_updated = FALSE;
		initialize_FITS_name_entries();	// To update UI with new preferences
		refresh_star_list();		// To update star list with new preferences
		if (com.found_object)
			refresh_found_objects();
		save_main_window_state(NULL);
		writeinitfile();
		validate_custom_profiles(); // Validate and load custom ICC profiles
		if (update_custom_gamut) {
			update_profiles_after_gamut_change();
			update_custom_gamut = FALSE;
		}
		refresh_icc_transforms();
		if (single_image_is_loaded() || sequence_is_loaded()) {
			color_manage(&gfit, gfit.color_managed);
			notify_gfit_modified();
		}
		siril_close_dialog("settings_window");
		update_reg_interface(TRUE); // To update UI with new preferences
	}
	// Update the GtkSourceView theme, if the widget has been created, in case
	// the Siril theme has been changed
	if (code_view_exists())
		set_code_view_theme();
}

void on_cancel_settings_button_clicked(GtkButton *button, gpointer user_data) {
	update_custom_gamut = FALSE;
	siril_close_dialog("settings_window");
}

void on_reset_settings_button_clicked(GtkButton *button, gpointer user_data) {
	int confirm = siril_confirm_dialog(_("Reset all preferences"),
			_("Do you really want to reset all preferences to default value?"),
			_("Reset Preferences"));
	if (confirm) {
		initialize_default_settings();
		update_preferences_from_model();
		update_custom_gamut = FALSE;
	}
}

void on_settings_window_hide(GtkWidget *widget, gpointer user_data) {
	update_custom_gamut = FALSE;
}

void on_external_preferred_profile_set(GtkFileChooser *chooser, gpointer user_data) {
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("working_gamut"));
	if (gtk_combo_box_get_active(combo) != 2)
		siril_message_dialog(GTK_MESSAGE_INFO, _("External ICC Profiles"), _("In order to "
		"use ICC profiles loaded from files as your preferred color space you must set "
		"the Preferred Color Space drop-down to \"From files\"."));
}

void on_working_gamut_changed(GtkComboBox *combo, gpointer user_data) {
	update_custom_gamut = TRUE;
}

void on_combo_rendering_intent_changed(GtkComboBox *combo, gpointer user_data) {
	if (gtk_combo_box_get_active(combo) == INTENT_ABSOLUTE_COLORIMETRIC) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_rendering_bpc")), FALSE);
	}
}

gchar *get_swap_dir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));

	if (sw_dir == NULL) {
		sw_dir = gtk_file_chooser_get_filename(swap_dir);
	}
	return sw_dir;
}

void on_button_python_reset_venv_clicked(gpointer user_data) {
	if (siril_confirm_dialog(_("WARNING!"), _("This will kill all running python processes "
			"and delete and reinstall the python venv directory. This is not normally "
			"necessary for upgrades or similar as Siril will try to manage the venv "
			"automatically. However the option is provided as a last resort bug mitigation "
			"to reset the venv to a known good state. All modules installed by scripts "
			"will be deleted and will require reinstallation."), _("Proceed"))) {
		rebuild_venv();
	}
}

/* these one are not on the preference dialog */

void on_cosmCFACheck_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.cfa = gtk_toggle_button_get_active(button);
}

void on_checkbutton_equalize_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.equalize_cfa = gtk_toggle_button_get_active(button);
}

void on_fix_xtrans_af_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.fix_xtrans = gtk_toggle_button_get_active(button);
}

/* ********************************** */
