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

#include <glib.h>
#include "core/siril.h"
#include "algos/colors.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/undo.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/icc_profile.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/image_interactions.h"
#include "gui/siril-window.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/siril_plot.h"
#include "gui/siril_preview.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/siril_plot.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "core/proto.h"

static cmsHPROFILE target = NULL; // Target profile for the GUI tool

static void export_elle_stone_profiles() {
	cmsHPROFILE profile;
	control_window_switch_to_tab(OUTPUT_LOGS);

	profile = srgb_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = srgb_trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = srgb_trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_srgbtrc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_srgbtrcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_rec709trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_rec709trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	// Also we use this function to export the ICC sRGB perceptual display profile,
	// but use a confirmation dialog to confirm the ICC's terms of use.
	if (siril_confirm_dialog(_("Terms of Use"), _("To export a copy of the sRGB "
			"monitor profile with perceptual intent tables, please accept the ICC's terms "
			"of use:\n\nTo anyone who acknowledges that the file \"sRGB_v4_ICC_preference.icc\" "
			"is provided \"AS IS\" WITH NO EXPRESS OR IMPLIED WARRANTY, permission to use, "
			"copy and distribute this file for any purpose is hereby granted without fee, "
			"provided that the file is not changed including the ICC copyright notice tag, "
			"and that the name of ICC shall not be used in advertising or publicity "
			"pertaining to distribution of the software without specific, written prior "
			"permission. ICC makes no representations about the suitability of this software "
			"for any purpose."), _("Accept"))) {
		profile = srgb_monitor_perceptual();
		export_profile(profile, "sRGB_v4_ICC_preference.icc");
		cmsCloseProfile(profile);
	}
}

void set_source_information() {
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_copyright_label");
	if (!gfit.color_managed) {
		siril_debug_print("Target image is not color managed\n");
		gtk_label_set_text(label, _("No ICC profile"));
		gtk_label_set_text(mfr_label, "");
		gtk_label_set_text(copyright_label, "");
		return;
	}
	if (!gfit.icc_profile) {
		siril_debug_print("Target profile is NULL\n");
		gtk_label_set_text(label, _("No ICC profile"));
		gtk_label_set_text(mfr_label, "");
		gtk_label_set_text(copyright_label, "");
		return;
	}
	// Set description
	gchar *buffer = siril_color_profile_get_description(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(label, buffer);
	free(buffer);

	// Set manufacturer
	buffer = siril_color_profile_get_manufacturer(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(mfr_label, buffer);
	free(buffer);

	// Set copyright
	buffer = siril_color_profile_get_copyright(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(copyright_label, buffer);
	free(buffer);
}

void set_icc_description_in_TIFF() {
	// Set description
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_save_label");
	gchar *buffer = NULL;
	if (gfit.icc_profile) {
		gtk_widget_set_tooltip_text((GtkWidget*) label, "");
		int length = cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", NULL, 0);
		if (length) {
			buffer = (char*) g_malloc(length * sizeof(char));
			cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", buffer, length);
		}
	} else {
			gtk_widget_set_tooltip_text((GtkWidget*) label, _("To write an ICC profile, assign a profile to the image using the Color Management dialog."));
			buffer = g_strdup(_("Image is not color managed. Will not write an ICC profile."));
	}
	gtk_label_set_text(label, buffer);
	g_free(buffer);
}

static void set_target_information() {
	if (!target) {
		return;
	}
	// Set description
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_target_mfr_label");
	GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_target_copyright_label");
	int length = cmsGetProfileInfoASCII(target, cmsInfoDescription, "en", "US", NULL, 0);
	char *buffer = NULL;
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoDescription, "en", "US", buffer, length);
		gtk_label_set_text(label, buffer);
		free(buffer);
	}

	// Set manufacturer
	length = cmsGetProfileInfoASCII(target, cmsInfoManufacturer, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoManufacturer, "en", "US", buffer, length);
		gtk_label_set_text(mfr_label, buffer);
		free(buffer);
	}

	// Set copyright
	length = cmsGetProfileInfoASCII(target, cmsInfoCopyright, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoCopyright, "en", "US", buffer, length);
		gtk_label_set_text(copyright_label, buffer);
		free(buffer);
	}
}

void initialize_icc_preferences_widgets() {
	GtkToggleButton *monitortogglebutton = (GtkToggleButton*) lookup_widget("custom_monitor_profile_active");
	GtkFileChooser *monitorfilechooser = (GtkFileChooser*) lookup_widget("pref_custom_monitor_profile");

	GtkToggleButton *proofingtogglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	GtkFileChooser *proofingfilechooser = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");

	if (!gtk_file_chooser_get_filename(monitorfilechooser)) {
		gtk_toggle_button_set_active(monitortogglebutton, FALSE);
		gtk_widget_set_sensitive((GtkWidget*) monitortogglebutton, FALSE);
	} else {
		gtk_toggle_button_set_active(monitortogglebutton, TRUE);
		gtk_widget_set_sensitive((GtkWidget*) monitortogglebutton, TRUE);
	}

	gtk_toggle_button_set_active(proofingtogglebutton, (gui.icc.soft_proof != NULL));
	if (!gtk_file_chooser_get_filename(proofingfilechooser)) {
		gtk_widget_set_sensitive((GtkWidget*) proofingtogglebutton, FALSE);
	} else {
		gtk_widget_set_sensitive((GtkWidget*) proofingtogglebutton, TRUE);
	}
}

void on_pref_custom_monitor_profile_file_set(GtkFileChooser* filechooser, gpointer user_data) {
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_monitor_profile_active");
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		gtk_widget_set_sensitive((GtkWidget*) togglebutton, TRUE);
	}
}

void on_pref_soft_proofing_profile_file_set(GtkFileChooser* filechooser, gpointer user_data) {
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		gtk_widget_set_sensitive((GtkWidget*) togglebutton, TRUE);
	}
}

void on_pref_icc_assign_never_toggled(GtkToggleButton *button, gpointer user_data);

void on_pref_icc_assign_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *never = (GtkToggleButton*) lookup_widget("pref_icc_assign_never");
	if (gtk_toggle_button_get_active(button)) {
		g_signal_handlers_block_by_func(never, on_pref_icc_assign_never_toggled, NULL);
		gtk_toggle_button_set_active(never, FALSE);
		g_signal_handlers_unblock_by_func(never, on_pref_icc_assign_never_toggled, NULL);
	}
}

void on_pref_icc_assign_never_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *load = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_load");
	GtkToggleButton *stack = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_stack");
	GtkToggleButton *stretch = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_stretch");
	GtkToggleButton *composition = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_composition");
	if (gtk_toggle_button_get_active(button)) {
		g_signal_handlers_block_by_func(load, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(stack, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(stretch, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(composition, on_pref_icc_assign_toggled, NULL);
		gtk_toggle_button_set_active(load, FALSE);
		gtk_toggle_button_set_active(stack, FALSE);
		gtk_toggle_button_set_active(stretch, FALSE);
		gtk_toggle_button_set_active(composition, FALSE);
		g_signal_handlers_unblock_by_func(load, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(stack, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(stretch, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(composition, on_pref_icc_assign_toggled, NULL);
	}
}
//////// GUI callbacks for the color management dialog

void on_icc_cancel_clicked(GtkButton* button, gpointer* user_data) {
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	GtkLabel* label2 = (GtkLabel*) lookup_widget("icc_target_mfr_label");
	GtkLabel* label3 = (GtkLabel*) lookup_widget("icc_target_copyright_label");
	GtkLabel* label4 = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* label5 = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* label6 = (GtkLabel*) lookup_widget("icc_copyright_label");
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	gtk_combo_box_set_active(target_combo, 0);
	gtk_file_chooser_unselect_all(filechooser);
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	gtk_label_set_text(label, "");
	gtk_label_set_text(label2, "");
	gtk_label_set_text(label3, "");
	gtk_label_set_text(label4, "");
	gtk_label_set_text(label5, "");
	gtk_label_set_text(label6, "");
siril_close_dialog("icc_dialog");
}

void on_icc_assign_clicked(GtkButton* button, gpointer* user_data) {
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No color profile chosen, nothing to assign."));
		return;
	}
	on_clear_roi();
	// We save the undo state as dealing with gfit
	undo_save_state(&gfit, _("Color profile assignment"));

	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);

	// Handle initial assignment of an ICC profile
	if (!gfit.color_managed || !gfit.icc_profile) {
		if (gfit.icc_profile) {
			cmsCloseProfile(gfit.icc_profile);
			gfit.icc_profile = NULL;
		}
		if (target_colorspace_channels > gfit.naxes[2]) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space has incorrect channels"), _("Mismatch in number of channels between the current image and the ICC profile. You cannot assign a RGB ICC profile to a mono image."));
			return;
		}
		goto FINISH;
	}

	cmsUInt32Number gfit_colorspace = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number gfit_colorspace_channels = cmsChannelsOf(gfit_colorspace);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}
	if (gfit_colorspace_channels != target_colorspace_channels) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Transform not supported"), _("Image cannot be assigned a color profile with a different number of channels to its current color profile"));
		return;
	}
	if (gfit.icc_profile) {
		cmsCloseProfile(gfit.icc_profile);
		gfit.icc_profile = NULL;
	}
FINISH:
	gfit.icc_profile = copyICCProfile(target);
	if (gfit.icc_profile)
		color_manage(&gfit, TRUE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	notify_gfit_modified();
}

void on_icc_remove_clicked(GtkButton* button, gpointer* user_data) {
	on_clear_roi();
	// We save the undo state as dealing with gfit
	undo_save_state(&gfit, _("Color profile removal"));
	if (gfit.icc_profile) {
		cmsCloseProfile(gfit.icc_profile);
		gfit.icc_profile = NULL;
	}
	color_manage(&gfit, FALSE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	notify_gfit_modified();

}

void on_icc_convertto_clicked(GtkButton* button, gpointer* user_data) {
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No color profile chosen, nothing to assign."));
		return;
	}
	on_clear_roi();
	if (!gfit.color_managed || !gfit.icc_profile) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No color profile set"), _("The current image has no color profile. You need to assign one first."));
		return;
	}

	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}

	// Do the transform
	// We save the undo state if dealing with gfit
	undo_save_state(&gfit, _("Color profile conversion"));
	cmsUInt32Number temp_intent = com.pref.icc.processing_intent;
	com.pref.icc.processing_intent = com.pref.icc.export_intent;
	siril_colorspace_transform(&gfit, target);
	com.pref.icc.processing_intent = temp_intent;

	// Assign the new color space to gfit
	if (gfit.icc_profile)
		cmsCloseProfile(gfit.icc_profile);
	gfit.icc_profile = copyICCProfile(target);
	if (gfit.icc_profile)
		color_manage(&gfit, TRUE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	close_tab();
	init_right_tab();
	notify_gfit_modified();
}

void on_icc_target_combo_changed(GtkComboBox* combo, gpointer* user_data) {
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	internal_icc target_index = gtk_combo_box_get_active(combo);
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	gtk_file_chooser_unselect_all(filechooser);
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	GtkLabel *mfr_label = NULL, *copyright_label = NULL;
	switch (target_index) {
		case NONE:
			mfr_label = (GtkLabel*) lookup_widget("icc_target_mfr_label");
			copyright_label = (GtkLabel*) lookup_widget("icc_target_copyright_label");
			gtk_label_set_text(label, "");
			gtk_label_set_text(copyright_label, "");
			gtk_label_set_text(mfr_label, "");
			return;
		case SRGB_LINEAR:
			target = srgb_linear();
			break;
		case SRGB_TRC:
			target = srgb_trc();
			break;
		case REC2020_LINEAR:
			target = rec2020_linear();
			break;
		case REC2020_TRC:
			target = rec2020_trc();
			break;
		case GRAY_LINEAR:
			target = gray_linear();
			break;
		case GRAY_SRGBTRC:
			target = gray_srgbtrc();
			break;
		case GRAY_REC709TRC:
			target = gray_rec709trc();
			break;
	}
	set_target_information();
}

void on_icc_target_filechooser_file_set(GtkFileChooser* filechooser, gpointer* user_data) {
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	g_signal_handlers_block_by_func(target_combo, on_icc_target_combo_changed, NULL);
	gtk_combo_box_set_active(target_combo, 0);
	g_signal_handlers_unblock_by_func(target_combo, on_icc_target_combo_changed, NULL);
	gchar *filename = siril_file_chooser_get_filename(filechooser);
	if (filename) {
		target = cmsOpenProfileFromFile(filename, "r");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error getting filename from widget."));
		gtk_file_chooser_unselect_all(filechooser);
		return;
	}
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error loading selected file or it does not contain a valid ICC profile."));
		gtk_file_chooser_unselect_all(filechooser);
	} else {
		cmsColorSpaceSignature target_signature = cmsGetColorSpace(target);
		if (target_signature == cmsSigRgbData) {
			set_target_information();
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Profile error"), _("Error: profile does not describe a RGB color space. Using non-RGB color spaces (e.g. CIE La*b*, XYZ, HSL) as image working spaces is not supported, though these color spaces may be used internally for some operations."));
			cmsCloseProfile(target);
			target = NULL;
			gtk_file_chooser_unselect_all(filechooser);
		}
	}
	g_free(filename);
}

gboolean on_icc_main_window_button_clicked(GtkWidget *btn, GdkEventButton *event, gpointer userdata) {
	if (!(single_image_is_loaded() || sequence_is_loaded()))
		return FALSE;
	if (event->type == GDK_BUTTON_PRESS  &&  event->button == 3) {
		// Right mouse button press
        if (gui.icc.iso12646) {
			siril_debug_print("Disabling approximate ISO12646 viewing conditions\n");
			disable_iso12646_conditions(TRUE, TRUE, TRUE);
		} else {
			siril_debug_print("Enabling approximate ISO12646 viewing conditions\n");
			enable_iso12646_conditions();
		}
        return TRUE;
	}
	if (event->type == GDK_BUTTON_PRESS  &&  event->button == 1) {
		// Left mouse button press
		siril_open_dialog("icc_dialog");
		return TRUE;
	}
	return FALSE;
}

void on_icc_dialog_show(GtkWidget *dialog, gpointer user_data) {
	set_source_information();
	set_target_information();
	GtkFileChooser* fc = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	gtk_file_chooser_set_current_folder(fc, default_system_icc_path());
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
}

void on_icc_export_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!gfit.icc_profile) {
		siril_log_color_message(_("Error: no ICC profile associated with the current image. Cannot export.\n"), "red");
		return;
	}
	char *filename = NULL;
	if (com.uniq->filename && com.uniq->filename[0] != '\0')
		filename = remove_ext_from_filename(com.uniq->filename);
	else
		filename = strdup("image");
	gchar* temp = g_strdup_printf("%s.icc", filename);
	free(filename);
	filename = strdup(temp);
	g_free(temp);
	if (cmsSaveProfileToFile(gfit.icc_profile, filename))
		siril_log_color_message(_("Exported ICC profile to %s\n"), "green", filename);
	else
		siril_log_color_message(_("Failed to export ICC profile to %s\n"), "red", filename);
	free(filename);
}

void on_icc_export_builtin_clicked(GtkButton *button, gpointer user_data) {
	export_elle_stone_profiles();
}

void on_icc_plot_clicked(GtkButton *button, gpointer user_data) {
	if (gfit.icc_profile && siril_color_profile_is_rgb (gfit.icc_profile)) {
		siril_plot_colorspace(gfit.icc_profile, TRUE);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Chromaticity plot only works with RGB color profiles"));
	}
}

static gboolean colorspace_comparison_image_set = FALSE;

void on_icc_gamut_visualisation_clicked() {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(lookup_widget("settings_window")));
	if (!colorspace_comparison_image_set) {
		GError *error = NULL;
		GdkPixbuf *pixbuf = gdk_pixbuf_new_from_resource("/org/siril/ui/pixmaps/CIE1931.svg", &error);
		if (error) {
			siril_debug_print("Error: %s\n", error->message);
			g_error_free(error);
		}
		gtk_image_set_from_pixbuf(GTK_IMAGE(lookup_widget("colorspace_comparison")), pixbuf);
		colorspace_comparison_image_set = TRUE;
	}
	/* Here this is wanted that we do not use siril_open_dialog */
	gtk_widget_show(win);
}

void on_icc_gamut_close_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_widget_hide(win);
}

gboolean iso_12646_draw_event(GtkWidget *widget, cairo_t *cr, gpointer data) {
	if (!(single_image_is_loaded() || (sequence_is_loaded())))
		return FALSE;
	GtkWidget *parent_widget = GTK_WIDGET(data);
	GtkWidget *red = lookup_widget("vbox_r");
	GtkWidget *green = lookup_widget("vbox_g");
	GtkWidget *blue = lookup_widget("vbox_b");
	GtkWidget *rgb = lookup_widget("vbox_rgb");
	GtkWidget *child_widget;
	if (parent_widget == red) {
		child_widget = lookup_widget("drawingarear");
	} else if (parent_widget == green) {
		child_widget = lookup_widget("drawingareag");
	} else if (parent_widget == blue) {
		child_widget = lookup_widget("drawingareab");
	} else if (parent_widget == rgb) {
		child_widget = lookup_widget("drawingareargb");
	} else {
		return FALSE;
	}

	// Get child widget's allocation
	GtkAllocation child_allocation;
	gtk_widget_get_allocation(child_widget, &child_allocation);

	// Translate child widget's coordinates to parent widget's drawing context
	gint white_border_width = 36;
	gdouble z = get_zoom_val();
	gdouble x0 = 0.0, y0 = 0.0;
	gint child_x0, child_y0, child_x1, child_y1;
	cairo_matrix_transform_point(&gui.display_matrix, &x0, &y0);
	gtk_widget_translate_coordinates(child_widget, parent_widget, (gint) x0, (gint) y0, &child_x0, &child_y0);
	gtk_widget_translate_coordinates(child_widget, parent_widget, (gint) ((z * gfit.rx) + white_border_width), (gint) ((z * gfit.ry) + white_border_width), &child_x1, &child_y1);

	// Get parent widget's allocation
	GtkAllocation parent_allocation;
	gtk_widget_get_allocation(parent_widget, &parent_allocation);

	// Fill parent widget's background with mid grey
	cairo_set_source_rgb(cr, 0.7647, 0.7647, 0.7647); // D50 Gray as recommended by ISO 12646
	cairo_rectangle(cr, 0, 0, parent_allocation.width, parent_allocation.height);
	cairo_fill(cr);

	// Draw white rectangle around the child widget with a margin of {white_border_width} pixels
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); // White color
	cairo_rectangle(cr, child_x0 - white_border_width, child_y0 - white_border_width,
					child_x1 + white_border_width, child_y1 + white_border_width - 16);
	cairo_fill(cr);
	// Return FALSE to propagate the draw event further
	return FALSE;
}

gboolean on_iso12646_panel_hide_completed(GtkWidget *widget, gpointer data) {
	gboolean *remap = (gboolean *) data;
	// Set the zoom value to 80% of the zoom-to-fit value
	update_zoom_fit_button();
	gui.zoom_value = -1.0;
	double z = 0.8 * get_zoom_val();
	gui.zoom_value = z;
	// Set the display offset. This has to be calculated from the width of the control
	// window less that of the pane button because the width of the drawarea won't
	// have updated in time if we have to hide the side pane.
	int window_width = gtk_widget_get_allocated_width(gui.view[RED_VPORT].drawarea);
	int window_height = gtk_widget_get_allocated_height(gui.view[RED_VPORT].drawarea);
	gui.display_offset.x = (window_width / 2 - (gfit.rx / 2) * z);
	gui.display_offset.y = (window_height / 2) - (gfit.ry / 2 * z);
	adjust_vport_size_to_image();
	queue_redraw(remap ? REMAP_ALL : REDRAW_IMAGE); // Has to do the remap in case the sliders have been changed
	gtk_widget_queue_draw(lookup_widget("vbox_r"));
	gtk_widget_queue_draw(lookup_widget("vbox_g"));
	gtk_widget_queue_draw(lookup_widget("vbox_b"));
	gtk_widget_queue_draw(lookup_widget("vbox_rgb"));
	return FALSE;
}

static gboolean panel_state = TRUE;
static double prior_zoom = -1;
static point prior_offset = { 0.0 , 0.0 };
static sliders_mode prior_sliders = MINMAX;
static WORD prior_lo = 0, prior_hi = 65535;
static display_mode prior_rendering_mode = LINEAR_DISPLAY;
static gboolean mode_changed = FALSE;

// This function overrides any GTK theme to assign a white border and D50 Gray
// background to the 4 vports to approximate ISO 12646 viewing conditions.
// It is recommended in conjunction with the excellent Equilux GTK theme.
void enable_iso12646_conditions() {
	// Cache prior state
	prior_zoom = get_zoom_val();
	memcpy(&prior_offset, &gui.display_offset, sizeof(point));
	prior_sliders = gui.sliders;
	prior_lo = gui.lo;
	prior_hi = gui.hi;
	prior_rendering_mode = gui.rendering_mode;
	if (prior_rendering_mode != LINEAR_DISPLAY)
		mode_changed = TRUE;
	// Add draw callbacks
	GtkWidget *parent_widget = lookup_widget("vbox_rgb");
	gui.icc.sh_rgb = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_b");
	gui.icc.sh_b = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_g");
	gui.icc.sh_g = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_r");
	gui.icc.sh_r = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	// Set the text color of the labels within the vbox
	PangoAttrList *attrs = pango_attr_list_new();
	PangoAttribute *attr = pango_attr_foreground_new(15 * PANGO_SCALE, 15 * PANGO_SCALE, 15 * PANGO_SCALE);
	pango_attr_list_insert(attrs, attr);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_red")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_green")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_blue")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_rgb")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_red")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_green")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_blue")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_rgb")), attrs);
	// Get the panel out of the way
	GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
	gtk_image_set_from_icon_name(image, "pan-start-symbolic", GTK_ICON_SIZE_BUTTON);
	GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
	GtkWidget *widget = gtk_paned_get_child2(paned);
	panel_state = gtk_widget_get_visible(widget);
	if (panel_state)
		gtk_widget_set_visible(widget, FALSE);
	// Set the sliders to min/max
	gboolean is_8bit = gfit.orig_bitpix == BYTE_IMG;
	gboolean remap = ((gui.lo == 0 && gui.hi == 65535) || (is_8bit && (gui.lo == 0 && gui.hi == 255))) || mode_changed;
	gui.sliders = USER;
	gui.lo = 0;
	gui.hi = is_8bit ? 255 : 65535;
	gui.rendering_mode = LINEAR_DISPLAY;
	set_display_mode();
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_user")), TRUE);
	set_cutoff_sliders_values(); // The redraw will happen in the idle
	gui.icc.iso12646 = TRUE;
	g_idle_add((GSourceFunc)on_iso12646_panel_hide_completed, &remap);
}

void disable_iso12646_conditions(gboolean revert_zoom, gboolean revert_panel, gboolean revert_rendering_mode) {
	GtkWidget *parent_widget = lookup_widget("vbox_rgb");
	if (gui.icc.sh_rgb)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_rgb);
	gui.icc.sh_rgb = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_r");
	if (gui.icc.sh_r)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_r);
	gui.icc.sh_r = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_g");
	if (gui.icc.sh_g)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_g);
	gui.icc.sh_g = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_b");
	if (gui.icc.sh_b)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_b);
	gui.icc.sh_b = 0;
	// Revert the text color of the labels within the vbox
    PangoAttrList *attrs = NULL;
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_red")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_green")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_blue")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_rgb")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_red")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_green")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_blue")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_rgb")), attrs);
	gui.icc.iso12646 = FALSE;
	if (revert_panel) {
		// Return the panel to its previous state
		GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
		if (!panel_state)
			gtk_image_set_from_icon_name(image, "pan-end-symbolic", GTK_ICON_SIZE_BUTTON);
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		GtkWidget *widget = gtk_paned_get_child2(paned);
		if (panel_state)
			gtk_widget_set_visible(widget, TRUE);
	}
	if (revert_zoom)
		gui.zoom_value = prior_zoom;
	memcpy(&gui.display_offset, &prior_offset, sizeof(point));
	gui.sliders = prior_sliders;
	gui.lo = prior_lo;
	gui.hi = prior_hi;
	if (gui.sliders == MINMAX)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_minmax")), TRUE);
	else if (gui.sliders == MIPSLOHI)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_hilo")), TRUE);
	set_cutoff_sliders_values();
	if (revert_rendering_mode) {
		gui.rendering_mode = prior_rendering_mode;
		set_display_mode();
	}
	if (mode_changed)
		redraw(REMAP_ALL);
	gtk_widget_queue_draw(lookup_widget("control_window"));
}
