/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/siril_language.h"
#include "core/settings.h"
#include "algos/photometry.h"
#include "algos/astrometry_solver.h"
#include "algos/annotate.h"
#include "gui/annotations_pref.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/script_menu.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"
#include "gui/siril_intro.h"
#include "gui/fix_xtrans_af.h"
#include "stacking/stacking.h"

#include "preferences.h"

#ifndef W_OK
#define W_OK 2
#endif

static gchar *sw_dir = NULL;
static gchar *st_dir = NULL;

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
	com.pref.debayer.top_down = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")));
	com.pref.debayer.xbayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")));
	com.pref.debayer.ybayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")));
	com.pref.debayer.bayer_pattern = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")));
	com.pref.debayer.bayer_inter = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")));
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
}

static void update_prepro_preferences() {
	if (com.pref.prepro.bias_lib)
		g_free(com.pref.prepro.bias_lib);

	com.pref.prepro.bias_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_bias_lib")));
	com.pref.prepro.use_bias_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias")));

	if (com.pref.prepro.bias_synth)
		g_free(com.pref.prepro.bias_synth);
	const gchar *entry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("bias_synth_entry")));

	if (entry) {
		com.pref.prepro.bias_synth = g_strdup(entry);
		com.pref.prepro.use_bias_synth = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias_bis")));
	}

	if (com.pref.prepro.dark_lib)
		g_free(com.pref.prepro.dark_lib);

	com.pref.prepro.dark_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_dark_lib")));
	com.pref.prepro.use_dark_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark")));

	if (com.pref.prepro.flat_lib)
		g_free(com.pref.prepro.flat_lib);

	com.pref.prepro.flat_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_flat_lib")));
	com.pref.prepro.use_flat_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat")));

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
	com.pref.phot_set.minval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")));
	com.pref.phot_set.maxval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")));
}

static void update_analysis_preferences() {
	com.pref.analysis.mosaic_panel = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinpanel")));
}

static void update_scripts_preferences() {
	com.pref.gui.script_path = get_list_from_preferences_dialog();
	com.pref.gui.warn_script_run = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")));
	com.pref.script_check_requires = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")));
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
}

static void update_FITS_options_preferences() {
	com.pref.comp.fits_enabled = !(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio"))));
	com.pref.comp.fits_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")));
	com.pref.comp.fits_quantization = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")));
	com.pref.comp.fits_hcompress_scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")));

	com.pref.rgb_aladin = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_aladin")));

	const gchar *ext = gtk_combo_box_get_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")));
	com.pref.ext = g_strdup(ext);

	com.pref.force_16bit = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_type"))) == 0;

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
}

static void update_misc_preferences() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	GtkFileChooser *starnet_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet"));

	com.pref.swap_dir = gtk_file_chooser_get_filename(swap_dir);

	com.pref.starnet_dir = gtk_file_chooser_get_filename(starnet_dir);

	com.pref.gui.silent_quit = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")));
	com.pref.gui.silent_linear = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskSave")));

	const gchar *copy = gtk_entry_get_text(GTK_ENTRY(lookup_widget("miscCopyright")));
	com.pref.copyright = g_strdup(copy);

	com.pref.check_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")));
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
	GtkWidget *spin = (GtkWidget *)user_data;
	gtk_widget_set_sensitive(spin, !gtk_toggle_button_get_active(button));
}

void initialize_path_directory(const gchar *path) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (swap_dir, path);
	} else {
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
	}
}

void initialize_starnet_directory(const gchar *path) {
	GtkFileChooser *starnet_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (starnet_dir, path);
	} else {
		gtk_file_chooser_set_filename (starnet_dir, g_get_tmp_dir());
	}
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
	GtkFileChooser *starnet_dir = GTK_FILE_CHOOSER(fileChooser);
	gchar *dir;

	dir = gtk_file_chooser_get_filename (starnet_dir);

	if (g_access(dir, W_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to write in this directory: %s\n"), "red", dir);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(starnet_dir, com.pref.starnet_dir);
		return;
	}
	if (st_dir) {
		g_free(st_dir);
		st_dir = gtk_file_chooser_get_filename(starnet_dir);
	}
}

void on_show_preview_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *label = lookup_widget("thumbnails_label_size");
	GtkWidget *box = lookup_widget("thumbnails_box_size");

	gtk_widget_set_sensitive(label, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(box, gtk_toggle_button_get_active(togglebutton));
}

void bias_text_insert(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);

	if (*position == 0 && text[0] != '=') {
		g_signal_handlers_block_by_func(G_OBJECT (editable), G_CALLBACK (bias_text_insert), data);
		gchar *new_str = g_strdup_printf("=%s", text);
		gtk_editable_insert_text(editable, new_str, length + 1, position);
		g_free(new_str);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable), G_CALLBACK (bias_text_insert), data);
		g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");
	}

}

void on_clear_bias_entry_clicked(GtkButton *button, gpointer user_data) {
	gtk_file_chooser_unselect_all((GtkFileChooser *)user_data);
}

void on_clear_dark_entry_clicked(GtkButton *button, gpointer user_data) {
	gtk_file_chooser_unselect_all((GtkFileChooser *)user_data);
}

void on_clear_flat_entry_clicked(GtkButton *button, gpointer user_data) {
	gtk_file_chooser_unselect_all((GtkFileChooser *)user_data);
}

void on_play_introduction_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("settings_window");
	start_intro_script();
}

void on_reload_script_button_clicked(GtkButton *button, gpointer user_data) {
	gchar *error;
	int retval = refresh_scripts(FALSE, &error);

	if (retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Cannot refresh script list"), error);
	}
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

static gboolean from_prefs_init = FALSE;

void update_preferences_from_model() {
	siril_debug_print("updating preferences GUI from settings data\n");
	preferences *pref = &com.pref;
	/* tab 1 */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_use_header")), pref->debayer.use_bayer_header);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")), pref->debayer.bayer_pattern);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")), pref->debayer.bayer_inter);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")), pref->debayer.xbayeroff);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")), pref->debayer.ybayeroff);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")), pref->debayer.top_down);

	/* tab 2*/
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_aladin")), pref->rgb_aladin);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio")), !pref->comp.fits_enabled);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_enabled_radio")), pref->comp.fits_enabled);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")), pref->comp.fits_method);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")), pref->comp.fits_quantization);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")), pref->comp.fits_hcompress_scale);
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")), pref->ext == NULL ? ".fit" : pref->ext);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_type")), pref->force_16bit ? 0 : 1);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbox_heterogeneous_fitseq")), pref->allow_heterogeneous_fitseq);

	/* tab 3 */
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

	/* tab 4 */
	if (pref->prepro.bias_lib && (g_file_test(pref->prepro.bias_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_bias_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias"));

		gtk_file_chooser_set_filename(button, pref->prepro.bias_lib);
		gtk_toggle_button_set_active(toggle, pref->prepro.use_bias_lib);
	}

	if (pref->prepro.bias_synth) {
		GtkEntry *entry = GTK_ENTRY(lookup_widget("bias_synth_entry"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias_bis"));

		gtk_entry_set_text(entry, pref->prepro.bias_synth);
		gtk_toggle_button_set_active(toggle, pref->prepro.use_bias_synth);
	}

	if (pref->prepro.dark_lib && (g_file_test(pref->prepro.dark_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_dark_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark"));

		gtk_file_chooser_set_filename(button, pref->prepro.dark_lib);
		gtk_toggle_button_set_active(toggle, pref->prepro.use_dark_lib);
	}

	if (pref->prepro.flat_lib && (g_file_test(pref->prepro.flat_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_flat_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat"));

		gtk_file_chooser_set_filename(button, pref->prepro.flat_lib);
		gtk_toggle_button_set_active(toggle, pref->prepro.use_flat_lib);
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

	/* tab 5 */
	from_prefs_init = TRUE;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")), pref->phot_set.outer);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")), pref->phot_set.inner);
	from_prefs_init = FALSE;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinAperture")), pref->phot_set.aperture);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("photometry_force_radius_button")), !pref->phot_set.force_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")), pref->phot_set.gain);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")), pref->phot_set.minval);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")), pref->phot_set.maxval);

	/* tab 6 */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinpanel")), pref->analysis.mosaic_panel);

	/* tab 7 */
	pref->gui.script_path = set_list_to_preferences_dialog(pref->gui.script_path);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")), pref->gui.warn_script_run);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")), pref->script_check_requires);

	/* tab 8 */
	siril_language_fill_combo(pref->lang);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_theme")), pref->gui.combo_theme);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("pref_fontsize")), pref->gui.font_scale);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("pref_iconstyle")), pref->gui.icon_symbolic);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck")), pref->gui.remember_windows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button")), pref->gui.show_thumbnails);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("thumbnails_box_size")), pref->gui.thumbnail_size == 256 ? 1 : 0);

	/* tab 9 */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")), pref->mem_mode == RATIO);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")), pref->mem_mode == AMOUNT);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio")), pref->memory_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount")), pref->memory_amount);

	/* tab 10 */
	initialize_path_directory(pref->swap_dir);
	initialize_starnet_directory(pref->starnet_dir);

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")), pref->gui.silent_quit);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskSave")), pref->gui.silent_linear);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("miscCopyright")), pref->copyright == NULL ? "" : pref->copyright);
#ifdef HAVE_JSON_GLIB
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")), pref->check_update);
#else
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")), FALSE);
	gtk_widget_set_sensitive(lookup_widget("miscAskUpdateStartup"), FALSE);
#endif
}

static void dump_ui_to_global_var() {
	siril_debug_print("updating settings from preferences GUI\n");
	/* tab 1 */
	update_debayer_preferences();
	/* tab 2 */
	update_FITS_options_preferences();
	/* tab 3 */
	update_astrometry_preferences();
	/* tab 4 */
	update_prepro_preferences();
	/* tab 5 */
	update_analysis_preferences();
	/* tab 6 */
	update_photometry_preferences();
	/* tab 7 */
	update_scripts_preferences();
	/* tab 8 */
	update_user_interface_preferences();
	/* tab 9 */
	update_performances_preferences();
	/* tab 10 */
	update_misc_preferences();
}

void on_settings_window_show(GtkWidget *widget, gpointer user_data) {
	siril_debug_print("show preferences window: updating it\n");
	update_preferences_from_model();
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

		initialize_FITS_name_entries(); // To update UI with new preferences
		refresh_star_list(com.stars); // To update star list with new preferences
		if (com.found_object)
			force_to_refresh_catalogue_list();
		save_main_window_state();
		writeinitfile();

		siril_close_dialog("settings_window");
	}
}

void on_cancel_settings_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("settings_window");
}

void on_reset_settings_button_clicked(GtkButton *button, gpointer user_data) {
	int confirm = siril_confirm_dialog(_("Reset all preferences"),
			_("Do you really want to reset all preferences to default value?"),
			_("Reset Preferences"));
	if (confirm) {
		initialize_default_settings();
		update_preferences_from_model();
	}
}

gchar *get_swap_dir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));

	if (sw_dir == NULL) {
		sw_dir = gtk_file_chooser_get_filename(swap_dir);
	}
	return sw_dir;
}

gchar *get_starnet_dir() {
	GtkFileChooser *starnet_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_starnet"));
	if (st_dir)
		st_dir = gtk_file_chooser_get_filename(starnet_dir);

	return sw_dir;
}

/* these one are not on the preference dialog */

void on_cosmCFACheck_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_checkbutton_equalize_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.equalize_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_fix_xtrans_af_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro.fix_xtrans = gtk_toggle_button_get_active(button);
	writeinitfile();
}

/* ********************************** */
