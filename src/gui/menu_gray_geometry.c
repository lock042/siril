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
	double clamping_factor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_rot_clamp")));

	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Rotation (%.1lfdeg, cropped=%s, clamped=%.1lf)"), angle,
			cropped ? "TRUE" : "FALSE", clamp ? clamping_factor : 0.0);
	verbose_rotate_image(fit, com.selection, angle, interpolation, cropped, clamp, clamping_factor);

	// the UI is still opened, need to reset selection
	// to current image size and reset rotation
	rectangle area = {0, 0, fit->rx, fit->ry};
	memcpy(&com.selection, &area, sizeof(rectangle));
	gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")), 0.);
	new_selection_zone();

	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}
void siril_rotate90() {
	set_cursor_waiting(TRUE);
	undo_save_state(&gfit, _("Rotation (90 deg)"));
	verbose_rotate_fast(&gfit, 90); // fast rotation, no interpolation, no crop
	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void siril_rotate270() {
	set_cursor_waiting(TRUE);
	undo_save_state(&gfit, _("Rotation (-90 deg)"));
	verbose_rotate_fast(&gfit, -90); // fast rotation, no interpolation, no crop
	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
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
		new_selection_zone();
	}
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
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void mirrory_gui(fits *fit) {
	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Mirror Y"));
	mirrory(fit, TRUE);
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/*************
 * RESAMPLE
 */
void on_button_resample_ok_clicked(GtkButton *button, gpointer user_data) {
	if (confirm_delete_wcs_keywords(&gfit)) {
		double sample[2];
		sample[0] = gtk_spin_button_get_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));
		sample[1] = gtk_spin_button_get_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));
		int interpolation = gtk_combo_box_get_active(
				GTK_COMBO_BOX(lookup_widget("combo_interpolation")));
		gboolean clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_rot_clamp")));
		double clamping_factor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_rot_clamp")));

		set_cursor_waiting(TRUE);
		int toX = round_to_int((sample[0] / 100.0) * gfit.rx);
		int toY = round_to_int((sample[1] / 100.0) * gfit.ry);
		undo_save_state(&gfit, _("Resample (%g - %g)"), sample[0] / 100.0,
				sample[1] / 100.0);
		verbose_resize_gaussian(&gfit, toX, toY, interpolation, clamp, clamping_factor);

		redraw(REMAP_ALL);
		redraw_previews();
		update_MenuItem();
		set_cursor_waiting(FALSE);
	}
}

void on_button_resample_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("resample_dialog");
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}

void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double yvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")),
				yvalue);
}

void on_button_sample_ratio_toggled(GtkToggleButton *button, gpointer user_data) {
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(button))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}
/**************
 * CROP
 */
void siril_crop() {
	undo_save_state(&gfit, _("Crop (x=%d, y=%d, w=%d, h=%d)"),
			com.selection.x, com.selection.y, com.selection.w,
			com.selection.h);
	if (is_preview_active()) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("A live preview session is active"),
				_("It is impossible to crop the image when a filter with preview session is active. "
						"Please consider to close the filter dialog first."));
		return;
	}
	crop(&gfit, &com.selection);
	delete_selected_area();
	reset_display_offset();
	update_zoom_label();
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
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

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	GtkEntry *cropped_entry = GTK_ENTRY(lookup_widget("cropped_entry"));

	args->seq = &com.seq;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	args->prefix = gtk_entry_get_text(cropped_entry);

	set_cursor_waiting(TRUE);
	crop_sequence(args);
}

void on_crop_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("crop_dialog");
}

