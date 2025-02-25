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
#include "algos/geometry.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "algos/astrometry_solver.h"

#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/registration_preview.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"
#include "menu_gray_geometry.h"

/**
 *  ROTATION
 */
static void rotate_gui(fits *fit) {
	if (!check_ok_if_cfa())
		return;
	if (com.selection.w == 0 || com.selection.h == 0) return;
	static GtkToggleButton *crop_rotation = NULL;
	double angle = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")));
	int interpolation = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_interpolation_rotation")));
	int cropped;

	if (crop_rotation == NULL) {
		crop_rotation = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_rotation_crop"));
	}
	cropped = gtk_toggle_button_get_active(crop_rotation);
	if ((!cropped) & (com.selection.w < gfit.rx || com.selection.h < gfit.ry)) {
		cropped = siril_confirm_dialog(_("Crop confirmation"), ("A selection is active and its size is smaller than the original image. Do you want to crop to current selection?"), _("Crop"));
		if (cropped)
			gtk_toggle_button_set_active(crop_rotation, TRUE);
	}
	gboolean clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_rot_clamp")));
	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Rotation (%.1lfdeg, cropped=%s, clamped=%s)"), angle,
			cropped ? "TRUE" : "FALSE", clamp ? "TRUE" : "FALSE");
	verbose_rotate_image(fit, com.selection, angle, interpolation, cropped, clamp);

	// the UI is still opened, need to reset selection
	// to current image size and reset rotation
	rectangle area = {0, 0, fit->rx, fit->ry};
	memcpy(&com.selection, &area, sizeof(rectangle));
	gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")), 0.);
	gui_function(new_selection_zone, NULL);

	update_zoom_label();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}
void siril_rotate90() {
	set_cursor_waiting(TRUE);
	undo_save_state(&gfit, _("Rotation (90 deg)"));
	verbose_rotate_fast(&gfit, 90); // fast rotation, no interpolation, no crop
	update_zoom_label();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}

void siril_rotate270() {
	set_cursor_waiting(TRUE);
	undo_save_state(&gfit, _("Rotation (-90 deg)"));
	verbose_rotate_fast(&gfit, -90); // fast rotation, no interpolation, no crop
	update_zoom_label();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}

void on_button_rotation_close_clicked(GtkButton *button, gpointer user_data) {
	delete_selected_area();
	siril_close_dialog("rotation_dialog");
}

void on_button_rotation_ok_clicked(GtkButton *button, gpointer user_data) {
	rotate_gui(&gfit);
}

void on_spin_rotation_value_changed(GtkSpinButton *button, gpointer user_data) {
	if (com.selection.w != 0 && com.selection.h != 0) {
		gui.rotation = gtk_spin_button_get_value(button);
		redraw(REDRAW_OVERLAY);
	}
}

void on_checkbutton_rotation_crop_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(button)) {
		rectangle area = {0, 0, gfit.rx, gfit.ry};
		memcpy(&com.selection, &area, sizeof(rectangle));
		gui_function(new_selection_zone, NULL);
	}
}

void on_combo_interpolation_rotation_changed(GtkComboBox *combo_box, gpointer user_data) {
	gint idx = gtk_combo_box_get_active(combo_box);

	gtk_widget_set_sensitive(lookup_widget("toggle_rot_clamp"), idx == OPENCV_CUBIC || idx == OPENCV_LANCZOS4);
}

/******
 * MIRROR
 */

void on_menuitem_mirrorx_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void on_mirrorx_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void on_menuitem_mirrory_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrory_gui(&gfit);
}

void on_mirrory_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrory_gui(&gfit);
}

void mirrorx_gui(fits *fit) {
	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Mirror X"));
	mirrorx(fit, TRUE);
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}

void mirrory_gui(fits *fit) {
	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Mirror Y"));
	mirrory(fit, TRUE);
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}

/*************
 * BINNING
 */

void on_button_binning_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (confirm_delete_wcs_keywords(&gfit)) {
		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

		gboolean mean = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_binning"))) == 0;
		int factor = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("spinbutton_binning")));

		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, _("Binning x%d (%s)"), factor, mean ? _("average") : _("sum"));
		fits_binning(&gfit, factor, mean);

		gui_function(update_MenuItem, NULL); // WCS not available anymore
		notify_gfit_modified();
	}
}

void on_button_binning_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("binxy_dialog");
}

/*************
 * RESAMPLE
 */
void on_button_resample_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (confirm_delete_wcs_keywords(&gfit)) {
		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

		double sample[2];
		sample[0] = gtk_spin_button_get_value( GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));
		sample[1] = gtk_spin_button_get_value( GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));
		int interpolation = gtk_combo_box_get_active( GTK_COMBO_BOX(lookup_widget("combo_interpolation")));
		gboolean clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_rot_clamp")));

		set_cursor_waiting(TRUE);
		int toX = round_to_int((sample[0] / 100.0) * gfit.rx);
		int toY = round_to_int((sample[1] / 100.0) * gfit.ry);
		undo_save_state(&gfit, _("Resample (%g - %g)"), sample[0] / 100.0, sample[1] / 100.0);
		verbose_resize_gaussian(&gfit, toX, toY, interpolation, clamp);

		gui_function(update_MenuItem, NULL); // WCS not available anymore
		notify_gfit_modified();
	}
}

void on_button_resample_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("resample_dialog");
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton, gpointer user_data);
void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton, gpointer user_data);
void on_spinbutton_resample_X_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data);
void on_spinbutton_resample_Y_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data);

static GtkWidget *spinbutton_resample_X = NULL, *spinbutton_resample_Y = NULL,
	*spinbutton_resample_X_px = NULL, *spinbutton_resample_Y_px = NULL,
	*button_sample_ratio = NULL, *combo_interpolation = NULL, *toggle_scale_clamp = NULL;

static void initialize_resample_widgets_if_needed() {
	if (!spinbutton_resample_X) {
		spinbutton_resample_X = lookup_widget("spinbutton_resample_X");
		spinbutton_resample_Y = lookup_widget("spinbutton_resample_Y");
		spinbutton_resample_X_px = lookup_widget("spinbutton_resample_X_px");
		spinbutton_resample_Y_px = lookup_widget("spinbutton_resample_Y_px");
		button_sample_ratio = lookup_widget("button_sample_ratio");
		combo_interpolation = lookup_widget("combo_interpolation");
		toggle_scale_clamp = lookup_widget("toggle_scale_clamp");
	}
}

static void pause_resample_signal_handlers() {
	g_signal_handlers_block_by_func(spinbutton_resample_X, on_spinbutton_resample_X_value_changed, NULL);
	g_signal_handlers_block_by_func(spinbutton_resample_Y, on_spinbutton_resample_Y_value_changed, NULL);
	g_signal_handlers_block_by_func(spinbutton_resample_X_px, on_spinbutton_resample_X_px_value_changed, NULL);
	g_signal_handlers_block_by_func(spinbutton_resample_Y_px, on_spinbutton_resample_Y_px_value_changed, NULL);
}

static void restart_resample_signal_handlers() {
	g_signal_handlers_unblock_by_func(spinbutton_resample_X, on_spinbutton_resample_X_value_changed, NULL);
	g_signal_handlers_unblock_by_func(spinbutton_resample_Y, on_spinbutton_resample_Y_value_changed, NULL);
	g_signal_handlers_unblock_by_func(spinbutton_resample_X_px, on_spinbutton_resample_X_px_value_changed, NULL);
	g_signal_handlers_unblock_by_func(spinbutton_resample_Y_px, on_spinbutton_resample_Y_px_value_changed, NULL);
}

void on_resample_dialog_show(GtkWidget *dialog, gpointer user_data) {
	pause_resample_signal_handlers();
	initialize_resample_widgets_if_needed();
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), gfit.rx);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), gfit.ry);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), 100.0);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), 100.0);
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double xvalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_X));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), round_to_int(gfit.rx * xvalue / 100.0));

	if (gtk_toggle_button_get_active(ratio)) {
		double yvalue = xvalue;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), yvalue);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), round_to_int(gfit.ry * yvalue / 100.0));
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double yvalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_Y));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), round_to_int(gfit.ry * yvalue / 100.0));

	if (gtk_toggle_button_get_active(ratio)) {
		double xvalue = yvalue;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), xvalue);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), round_to_int(gfit.rx * xvalue / 100.0));
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_X_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio_button = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double xpix = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px));
	double ratio = xpix / gfit.rx;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), ratio * 100.0);

	if (gtk_toggle_button_get_active(ratio_button)) {
		double ypix = round_to_int(gfit.ry * ratio);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), ratio * 100.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), ypix);
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_Y_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio_button = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double ypix = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px));
	double ratio = ypix / gfit.ry;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), ratio * 100.0);

	if (gtk_toggle_button_get_active(ratio_button)) {
		double xpix = round_to_int(gfit.rx * ratio);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), ratio * 100.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), xpix);
	}
	restart_resample_signal_handlers();
}

void on_button_sample_ratio_toggled(GtkToggleButton *button, gpointer user_data) {
	double xvalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_X));

	if (gtk_toggle_button_get_active(button))
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), xvalue);
}

void on_combo_interpolation_changed(GtkComboBox *combo_box, gpointer user_data) {
	gint idx = gtk_combo_box_get_active(combo_box);

	gtk_widget_set_sensitive(toggle_scale_clamp, idx == OPENCV_CUBIC || idx == OPENCV_LANCZOS4);
}

/**************
 * CROP
 */
void siril_crop() {
	undo_save_state(&gfit, _("Crop (x=%d, y=%d, w=%d, h=%d)"),
			com.selection.x, com.selection.y, com.selection.w,
			com.selection.h);
	if (is_preview_active()) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Image backup is active"),
				_("It is impossible to crop the image when the image backup is active. "
				"This occurs when a filter has a live preview active or when a filter with "
				"on-demand ROI preview is open. Please close the filter dialog before "
				"cropping."));
		return;
	}
	clear_stars_list(TRUE);
	crop(&gfit, &com.selection);

	char log[90];
	sprintf(log, _("Crop (x=%d, y=%d, w=%d, h=%d)"), com.selection.x,
			com.selection.y, com.selection.w, com.selection.h);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	delete_selected_area();
	reset_display_offset();
	update_zoom_label();
	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}

void on_crop_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

#ifdef HAVE_FFMS2
	if (com.seq.type == SEQ_AVI) {
		siril_log_message(_("Crop does not work with "
				"avi film. Please, convert your file to SER first.\n"));
		return;
	}
#endif
	if (com.seq.type == SEQ_INTERNAL) {
		siril_log_message(_("Not a valid sequence for cropping.\n"));
	}

	struct crop_sequence_data *args = calloc(1, sizeof(struct crop_sequence_data));

	GtkEntry *cropped_entry = GTK_ENTRY(lookup_widget("cropped_entry"));

	args->seq = &com.seq;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	args->prefix = strdup(gtk_entry_get_text(cropped_entry));

	set_cursor_waiting(TRUE);
	crop_sequence(args);
}

void on_crop_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("crop_dialog");
}

