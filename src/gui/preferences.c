/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "algos/photometry.h"
#include "gui/callbacks.h"
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

static preferences pref_init = {
		TRUE, // first_start
		TRUE, // remember_windows
		{ 0, 0, 0, 0 }, // main_w_pos
		FALSE, // is_maximized
		FALSE, // prepro_cfa
		TRUE, // prepro_equalize_cfa
		FALSE, // fix_xtrans
		{ 0, 0, 0, 0 }, // xtrans_af
		{ 0, 0, 0, 0 }, // xtrans_sample
		NULL, // prepro_bias_lib
		FALSE, // use_bias_lib
		NULL, // prepro_dark_lib
		FALSE, // use_dark_lib
		NULL, // prepro_flat_lib
		FALSE, // use_flat_lib
		{ // save structure
				TRUE, // show quit
				TRUE, // show scripts message
		},
		TRUE, // show_thumbnails
		256, //thumbnail_size
		TRUE, // check_update
		TRUE, // check_script_version
		0, // combo_theme
		100, // font_scale
		0, // combo_lang
		NULL, // ext
		NULL, // swap_dir
		NULL, // script_path
		{
				{1.0, 1.0, 1.0}, // mul[3]
				1.0, // bright
				1, //auto_mul
				0, // use_camera_wb
				0, // use_auto_wb
				2, // user_qual
				0, // user_black
				{ 1.0, 1.0 }, //gamm[2]
		},
		{
				FALSE, // open_debayer
				TRUE, // use_bayer_header
				BAYER_FILTER_RGGB, // bayer_pattern
				BAYER_RCD, // bayer_inter
				TRUE, // top_down
				0, // xbayeroff
				0, // ybayeroff
		},
		{
				2.3, // gain
				20.0, // inner
				30.0, // outer
				0, // minval
				60000, // maxval
		},
		{
				0, // stack method
				ADDITIVE_SCALING, // normalisation
				WINSORIZED, // type of rejection
				3.0, 3.0, // sigma low sigma high
				5.0, 5.0, // linear low linear high
				3.0, 3.0, // percentile low percentile high
				RATIO, // memory management
				0.9, // memory ratio
				10, // memory amount
		},
		{
				FALSE, // fits_enabled
				0, // fits_method
				16.0, // fits_quantization
				4.0, // fits_hcompress_scale
		},
		FALSE, // force_to_16bit
		0, // selection_guides
		NULL, // copyright
};

void on_cosmCFACheck_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_checkbutton_equalize_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.prepro_equalize_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_fix_xtrans_af_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.fix_xtrans = gtk_toggle_button_get_active(button);
	writeinitfile();
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

static void update_libraw_preferences() {
	/**********COLOR ADJUSTEMENT**************/
	com.pref.raw_set.bright = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")));
	com.pref.raw_set.mul[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")));
	com.pref.raw_set.mul[2] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")));

	com.pref.raw_set.auto_mul = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")));
	com.pref.raw_set.user_black = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")));

	/**************WHITE BALANCE**************/
	com.pref.raw_set.use_camera_wb = (int) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")));
	com.pref.raw_set.use_auto_wb = (int) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")));

	/********MATRIX INTERPOLATION**************/
	com.pref.raw_set.user_qual = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")));

	/********GAMMA CORRECTION**************/
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")))) {
		/* Linear Gamma Curve */
		com.pref.raw_set.gamm[0] = 1.0;
		com.pref.raw_set.gamm[1] = 1.0;
	} else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")))) {
		/* BT.709 Gamma curve */
		com.pref.raw_set.gamm[0] = 2.222;
		com.pref.raw_set.gamm[1] = 4.5;
	} else {
		/* sRGB Gamma curve */
		com.pref.raw_set.gamm[0] = 2.40;
		com.pref.raw_set.gamm[1] = 12.92;
	}
}

void update_debayer_preferences() {
	com.pref.debayer.use_bayer_header = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header")));
	com.pref.debayer.top_down = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")));
	com.pref.debayer.xbayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")));
	com.pref.debayer.ybayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")));
	com.pref.debayer.bayer_pattern = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")));
	com.pref.debayer.bayer_inter = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")));
}

void update_prepro_preferences() {
	if (com.pref.prepro_bias_lib) {
		g_free(com.pref.prepro_bias_lib);
		com.pref.prepro_bias_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_bias_lib")));
	}
	com.pref.use_bias_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias")));
	if (com.pref.prepro_dark_lib) {
		g_free(com.pref.prepro_dark_lib);
		com.pref.prepro_dark_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_dark_lib")));
	}
	com.pref.use_dark_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark")));
	if (com.pref.prepro_flat_lib) {
		g_free(com.pref.prepro_flat_lib);
		com.pref.prepro_flat_lib = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(lookup_widget("filechooser_flat_lib")));
	}
	com.pref.use_flat_lib = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat")));

	com.pref.xtrans_af.x = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_x"))), NULL, 10);
	com.pref.xtrans_af.y = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_y"))), NULL, 10);
	com.pref.xtrans_af.w = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_w"))), NULL, 10);
	com.pref.xtrans_af.h = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_af_h"))), NULL, 10);

	com.pref.xtrans_sample.x = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_x"))), NULL, 10);
	com.pref.xtrans_sample.y = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_y"))), NULL, 10);
	com.pref.xtrans_sample.w = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_w"))), NULL, 10);
	com.pref.xtrans_sample.h = g_ascii_strtoull(gtk_entry_get_text(GTK_ENTRY(lookup_widget("xtrans_sample_h"))), NULL, 10);
}

static void update_photometry_preferences() {
	com.pref.phot_set.gain = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")));
	com.pref.phot_set.inner = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")));
	com.pref.phot_set.outer = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")));
	com.pref.phot_set.minval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")));
	com.pref.phot_set.maxval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")));
}

static void update_scripts_preferences() {
	com.pref.script_path = get_list_from_preferences_dialog();
	com.pref.save.quit = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")));
	com.pref.save.warn_script = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")));
}

static void update_user_interface_preferences() {
	com.pref.combo_lang = get_interface_language();
	com.pref.combo_theme = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_theme")));
	com.pref.font_scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("pref_fontsize")));
	com.pref.remember_windows = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck")));
	com.pref.show_thumbnails = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button")));
	com.pref.thumbnail_size = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("thumbnails_box_size"))) == 1 ? 256 : 128;
}

static void update_performances_preferences() {
	GSList *tmp_list = gtk_radio_button_get_group (GTK_RADIO_BUTTON(lookup_widget("memfreeratio_radio")));
	GtkWidget *ratio, *amount;
	GtkToggleButton *tmp_button = NULL;//Create a temp toggle button.

	ratio = lookup_widget("memfreeratio_radio");
	amount = lookup_widget("memfixed_radio");
	while (tmp_list) {
		tmp_button = tmp_list->data;
		tmp_list = tmp_list->next;

		if (gtk_toggle_button_get_active(tmp_button))
			break;

		tmp_button = NULL; //We've enumerated all of them, and none of them is active.
	}
	if (tmp_button) {
		if (GTK_WIDGET(tmp_button) == ratio) {
			com.pref.stack.mem_mode = RATIO;
		} else if (GTK_WIDGET(tmp_button) == amount) {
			com.pref.stack.mem_mode = AMOUNT;
		} else {
			com.pref.stack.mem_mode = UNLIMITED;
		}
	}

	com.pref.stack.memory_ratio = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio")));
	com.pref.stack.memory_amount = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount")));

	com.pref.comp.fits_enabled = !(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio"))));
	com.pref.comp.fits_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")));
	com.pref.comp.fits_quantization = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")));
	com.pref.comp.fits_hcompress_scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")));
}

static void update_misc_preferences() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));

	com.pref.swap_dir = gtk_file_chooser_get_filename(swap_dir);

	const gchar *ext = gtk_combo_box_get_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")));
	com.pref.ext = g_strdup(ext);

	com.pref.force_to_16bit = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_type"))) == 0;

	com.pref.save.quit = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")));

	const gchar *copy = gtk_entry_get_text(GTK_ENTRY(lookup_widget("miscCopyright")));
	com.pref.copyright = g_strdup(copy);

	com.pref.check_update = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")));
}

void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(auto_button)) {
		g_signal_handlers_block_by_func(auto_button, on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, FALSE);
		g_signal_handlers_unblock_by_func(auto_button, on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, TRUE);
	}
}

void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(cam_button)) {
		g_signal_handlers_block_by_func(cam_button, on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, FALSE);
		g_signal_handlers_unblock_by_func(cam_button, on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, TRUE);
	}
}

void on_checkbutton_auto_evaluate_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *entry = (GtkWidget *)user_data;
	gtk_widget_set_sensitive(entry, !gtk_toggle_button_get_active(button));
}

void on_checkbutton_multipliers_toggled(GtkToggleButton *button,
		gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(lookup_widget("hbox8"), !active);
	gtk_widget_set_sensitive(lookup_widget("hbox11"), !active);
	if (active) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), 1.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), 1.0);
	}
}

void set_GUI_photometry() {
	if (gfit.cvf > 0.0) {
		com.pref.phot_set.gain = gfit.cvf;
	}
	if (com.pref.phot_set.gain > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")), com.pref.phot_set.gain);
	}
	if (com.pref.phot_set.inner > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")), com.pref.phot_set.inner);
	}
	if (com.pref.phot_set.outer > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")), com.pref.phot_set.outer);
	}
	if (com.pref.phot_set.minval >= 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")), (gdouble) com.pref.phot_set.minval);
	}
	if (com.pref.phot_set.maxval >= 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")), (gdouble) com.pref.phot_set.maxval);
	}
}

void initialize_path_directory(const gchar *path) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	if (path && path[0] != '\0') {
		gtk_file_chooser_set_filename (swap_dir, path);
	} else {
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
	}
}

void set_libraw_settings_menu_available(gboolean activate) {
	if (!com.script) {
		gtk_widget_set_visible(lookup_widget("box_stack_p1"), activate);
	}
}

void on_button_reset_swap_clicked(GtkButton *button, gpointer user_data) {
	reset_swapdir();
}


void on_spinbutton_comp_fits_quantization_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble quantization = gtk_spin_button_get_value(button);
	if (quantization == 0.0) {
		GtkComboBox *combo = (GtkComboBox *)user_data;
		if (gtk_combo_box_get_active(combo) != GZIP1_COMP && gtk_combo_box_get_active(combo) != GZIP2_COMP) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Incorrect parameters detected"),
								"Setting quantization to 0 has only a sense with a GZIP compression "
								"and GZIP 2 often produces better compression of floating­point images.");
			gtk_spin_button_set_value(button, com.pref.comp.fits_quantization != 0.0 ? com.pref.comp.fits_quantization : 16.0);
		}
	}
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

void set_GUI_compression() {
	if (com.pref.comp.fits_enabled) {
		GtkToggleButton *enabled = GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_enabled_radio"));
		gtk_toggle_button_set_active(enabled, com.pref.comp.fits_enabled);
		GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method"));
		GtkSpinButton *quantiz = GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization"));
		GtkSpinButton *hscale = GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale"));

		gtk_combo_box_set_active(box, com.pref.comp.fits_method);
		gtk_spin_button_set_value(quantiz, com.pref.comp.fits_quantization);
		if (com.pref.comp.fits_method == HCOMPRESS_COMP) {
			gtk_spin_button_set_value(hscale, com.pref.comp.fits_hcompress_scale);
		}
		siril_log_message(_("FITS compression enabled\n"), com.pref.ext);
	} else {
		GtkToggleButton *disabled = GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio"));
		gtk_toggle_button_set_active(disabled, !com.pref.comp.fits_enabled);
		siril_log_message(_("FITS compression disabled\n"));
	}
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
}


void on_show_preview_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *label = lookup_widget("thumbnails_label_size");
	GtkWidget *box = lookup_widget("thumbnails_box_size");

	gtk_widget_set_sensitive(label, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(box, gtk_toggle_button_get_active(togglebutton));
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

void on_combo_theme_changed(GtkComboBox *box, gpointer user_data) {
	GtkSettings *settings;

	int active = gtk_combo_box_get_active(box);

	settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", active == 0, NULL);
	update_icons_to_theme(active == 0);
}

void on_reload_script_button_clicked(GtkButton *button, gpointer user_data) {
	gchar *error;
	int retval = refresh_scripts(FALSE, &error);

	if (retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Cannot refresh script list"), error);
	}
}

void on_spinInner_value_changed(GtkSpinButton *inner, gpointer user_data) {
	GtkSpinButton *outer;
	double in, out;

	outer = GTK_SPIN_BUTTON(lookup_widget("spinOuter"));
	in = gtk_spin_button_get_value(inner);
	out = gtk_spin_button_get_value(outer);

	if (in >= out) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong value"),
				_("Inner radius value must be less than outer. Please change the value."));
	}
}

void on_spinOuter_value_changed(GtkSpinButton *outer, gpointer user_data) {
	GtkSpinButton *inner;
	double in, out;

	inner = GTK_SPIN_BUTTON(lookup_widget("spinInner"));
	in = gtk_spin_button_get_value(inner);
	out = gtk_spin_button_get_value(outer);

	if (in >= out) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong value"),
				_("Inner radius value must be less than outer. Please change the value."));
	}
}

static void set_preferences_ui(preferences *pref) {
	/* tab 1 */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")), pref->raw_set.auto_mul);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), pref->raw_set.mul[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), pref->raw_set.mul[2]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")), pref->raw_set.bright);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")), pref->raw_set.user_black);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")), pref->raw_set.gamm[0]);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")), pref->raw_set.user_qual);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")), pref->raw_set.use_camera_wb);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")), pref->raw_set.use_auto_wb);

	/* tab 2 */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header")), pref->debayer.use_bayer_header);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_pattern")), pref->debayer.bayer_pattern);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboBayer_inter")), pref->debayer.bayer_inter);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")), pref->debayer.xbayeroff);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")), pref->debayer.ybayeroff);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")), pref->debayer.top_down);

	/* tab 3 */
	if (pref->prepro_bias_lib && (g_file_test(pref->prepro_bias_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_bias_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_bias"));

		gtk_file_chooser_set_filename(button, pref->prepro_bias_lib);
		gtk_toggle_button_set_active(toggle, pref->use_bias_lib);
	}

	if (pref->prepro_dark_lib && (g_file_test(pref->prepro_dark_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_dark_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_dark"));

		gtk_file_chooser_set_filename(button, pref->prepro_dark_lib);
		gtk_toggle_button_set_active(toggle, pref->use_dark_lib);
	}

	if (pref->prepro_flat_lib && (g_file_test(pref->prepro_flat_lib, G_FILE_TEST_EXISTS))) {
		GtkFileChooser *button = GTK_FILE_CHOOSER(lookup_widget("filechooser_flat_lib"));
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("check_button_pref_flat"));

		gtk_file_chooser_set_filename(button, pref->prepro_flat_lib);
		gtk_toggle_button_set_active(toggle, pref->use_flat_lib);
	}

	gchar tmp[256];
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_af.x);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_x")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_af.y);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_y")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_af.w);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_w")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_af.h);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_af_h")), tmp);

	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_sample.x);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_x")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_sample.y);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_y")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_sample.w);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_w")), tmp);
	g_snprintf(tmp, sizeof(tmp), "%d", pref->xtrans_sample.h);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("xtrans_sample_h")), tmp);

	/* tab 4 */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")), pref->phot_set.inner);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")), pref->phot_set.outer);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")), pref->phot_set.gain);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")), pref->phot_set.minval);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")), pref->phot_set.maxval);

	/* tab 5 */
	pref->script_path = set_list_to_preferences_dialog(pref->script_path);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")), pref->save.warn_script);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("script_check_version")), pref->check_update);

	/* tab 6 */
	siril_language_fill_combo(pref->combo_lang);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_theme")), pref->combo_theme);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("pref_fontsize")), pref->font_scale);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck")), pref->remember_windows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button")), pref->show_thumbnails);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("thumbnails_box_size")), pref->thumbnail_size == 256 ? 1 : 0);

	/* tab 7 */
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")), pref->stack.mem_mode == RATIO);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")), pref->stack.mem_mode == AMOUNT);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("memunlimited_radio")), pref->stack.mem_mode == UNLIMITED);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio")), pref->stack.memory_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount")), pref->stack.memory_amount);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("comp_fits_disabled_radio")), !pref->comp.fits_enabled);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_comp_fits_method")), pref->comp.fits_method);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_quantization")), pref->comp.fits_quantization);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_comp_fits_hcompress_scale")), pref->comp.fits_hcompress_scale);

	/* tab 8 */
	initialize_path_directory(pref->swap_dir);
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(lookup_widget("combobox_ext")), pref->ext == NULL ? ".fit" : pref->ext);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_type")), pref->force_to_16bit ? 0 : 1);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit")), pref->save.quit);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("miscCopyright")), pref->copyright == NULL ? "" : pref->copyright);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup")), pref->check_update);
}

static void dump_ui_to_global_var() {
	/* tab 1 */
	update_libraw_preferences();
	/* tab 2 */
	update_debayer_preferences();
	/* tab 3 */
	update_prepro_preferences();
	/* tab 4 */
	update_photometry_preferences();
	/* tab 5 */
	update_scripts_preferences();
	/* tab 6 */
	update_user_interface_preferences();
	/* tab 7 */
	update_performances_preferences();
	/* tab 8 */
	update_misc_preferences();
}

static void reset_preferences() {
	memcpy(&com.pref, &pref_init, sizeof(preferences));
}

static void free_preferences(preferences *pref) {
	g_free(pref->ext);
	pref->ext = NULL;
	g_free(pref->swap_dir);
	pref->swap_dir = NULL;
	g_free(pref->copyright);
	pref->copyright = NULL;
	g_free(pref->combo_lang);
	pref->combo_lang = NULL;
	g_slist_free_full(pref->script_path, g_free);
	pref->script_path = NULL;
}

void initialize_default_preferences() {
	reset_preferences();
}

void on_apply_settings_button_clicked(GtkButton *button, gpointer user_data) {
	free_preferences(&com.pref);
	dump_ui_to_global_var();

	initialize_FITS_name_entries(); // To update UI with new preferences
	refresh_star_list(com.stars); // To update star list with new preferences
	save_main_window_state();
	siril_close_dialog("settings_window");
	writeinitfile();
}

void on_cancel_settings_button_clicked(GtkButton *button, gpointer user_data) {
	set_preferences_ui(&com.pref);
	siril_close_dialog("settings_window");
}

void on_reset_settings_button_clicked(GtkButton *button, gpointer user_data) {
	int confirm = siril_confirm_dialog(_("Reset all preferences"), _("Do you really want to reset all preferences to default value?"));
	if (confirm) {
		set_preferences_ui(&pref_init);
	}
}

gchar *get_swap_dir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));

	return gtk_file_chooser_get_filename(swap_dir);
}

void set_preferences_dialog_from_global() {
	set_preferences_ui(&com.pref);
}
