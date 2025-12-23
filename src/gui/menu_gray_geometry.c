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

/* Idle function for rotation completion */
static gboolean rotation_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		// Reset selection to current image size and reset rotation
		rectangle area = {0, 0, gfit->rx, gfit->ry};
		memcpy(&com.selection, &area, sizeof(rectangle));
		gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")), 0.);
		gui_function(new_selection_zone, NULL);

		update_zoom_label();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

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
	if ((!cropped) & (com.selection.w < gfit->rx || com.selection.h < gfit->ry)) {
		cropped = siril_confirm_dialog(_("Crop confirmation"),
			("A selection is active and its size is smaller than the original image. Do you want to crop to current selection?"),
			_("Crop"));
		if (cropped)
			gtk_toggle_button_set_active(crop_rotation, TRUE);
	}
	gboolean clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_rot_clamp")));

	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct rotation_args *params = new_rotation_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->area = com.selection;
	params->angle = angle;
	params->interpolation = interpolation;
	params->cropped = cropped;
	params->clamp = clamp;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rotation_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = fit;
	args->mem_ratio = 2.0f;  // Rotation needs space for transformation
	args->image_hook = rotation_image_hook;
	args->mask_hook = rotation_mask_hook;
	args->log_hook = rotation_log_hook;
	args->idle_function = rotation_idle;
	args->description = _("Rotation");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
		return;
	}
}

/* Idle function for fast rotation */
static gboolean fast_rotation_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		update_zoom_label();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

void siril_rotate90() {
	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct rotation_args *params = new_rotation_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->area = (rectangle){0, 0, gfit->rx, gfit->ry};
	params->angle = 90.0;
	params->interpolation = -1;  // Fast rotation
	params->cropped = 0;
	params->clamp = FALSE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rotation_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = rotation_image_hook;
	args->mask_hook = rotation_mask_hook;
	args->log_hook = rotation_log_hook;
	args->idle_function = fast_rotation_idle;
	args->description = _("Rotation 90°");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = 1;  // Fast rotation doesn't benefit from threading

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void siril_rotate270() {
	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct rotation_args *params = new_rotation_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->area = (rectangle){0, 0, gfit->rx, gfit->ry};
	params->angle = -90.0;
	params->interpolation = -1;  // Fast rotation
	params->cropped = 0;
	params->clamp = FALSE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rotation_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = rotation_image_hook;
	args->mask_hook = rotation_mask_hook;
	args->log_hook = rotation_log_hook;
	args->idle_function = fast_rotation_idle;
	args->description = _("Rotation -90°");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = 1;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void on_button_rotation_close_clicked(GtkButton *button, gpointer user_data) {
	delete_selected_area();
	siril_close_dialog("rotation_dialog");
}

gboolean rotation_hide_on_delete(GtkWidget *widget) {
	delete_selected_area();
	siril_close_dialog("rotation_dialog");
	return TRUE;
}

void on_button_rotation_ok_clicked(GtkButton *button, gpointer user_data) {
	rotate_gui(gfit);
}

void on_spin_rotation_value_changed(GtkSpinButton *button, gpointer user_data) {
	if (com.selection.w != 0 && com.selection.h != 0) {
		gui.rotation = gtk_spin_button_get_value(button);
		redraw(REDRAW_OVERLAY);
	}
}

void on_checkbutton_rotation_crop_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(button)) {
		rectangle area = {0, 0, gfit->rx, gfit->ry};
		memcpy(&com.selection, &area, sizeof(rectangle));
		gui_function(new_selection_zone, NULL);
	}
}

void on_combo_interpolation_rotation_changed(GtkComboBox *combo_box, gpointer user_data) {
	gint idx = gtk_combo_box_get_active(combo_box);
	gtk_widget_set_sensitive(lookup_widget("toggle_rot_clamp"),
	                         idx == OPENCV_CUBIC || idx == OPENCV_LANCZOS4);
}

/******
 * MIRROR
 */

/* Idle function for mirror operations */
static gboolean mirror_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

void on_menuitem_mirrorx_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrorx_gui(gfit);
}

void on_mirrorx_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrorx_gui(gfit);
}

void on_menuitem_mirrory_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrory_gui(gfit);
}

void on_mirrory_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrory_gui(gfit);
}

void mirrorx_gui(fits *fit) {
	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct mirror_args *params = new_mirror_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->x_axis = TRUE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_mirror_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = fit;
	args->mem_ratio = 1.0f;  // Mirror needs minimal extra memory
	args->image_hook = mirrorx_image_hook;
	args->mask_hook = mirrorx_mask_hook;
	args->idle_function = mirror_idle;
	args->description = _("Mirror X");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void mirrory_gui(fits *fit) {
	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct mirror_args *params = new_mirror_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->x_axis = FALSE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_mirror_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = fit;
	args->mem_ratio = 1.0f;
	args->image_hook = mirrory_image_hook;
	args->mask_hook = mirrory_mask_hook;
	args->idle_function = mirror_idle;
	args->description = _("Mirror Y");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

/*************
 * BINNING
 */

/* Idle function for binning */
static gboolean binning_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		gui_function(update_MenuItem, NULL); // WCS not available anymore
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

void on_button_binning_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!confirm_delete_wcs_keywords(gfit))
		return;

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	gboolean mean = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combobox_binning"))) == 0;
	int factor = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("spinbutton_binning")));

	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct binning_args *params = new_binning_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->factor = factor;
	params->mean = mean;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_binning_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;  // Binning needs space for the new buffer
	args->image_hook = binning_image_hook;
	args->mask_hook = binning_mask_hook;
	args->log_hook = binning_log_hook;
	args->idle_function = binning_idle;
	args->description = _("Binning");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void on_button_binning_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("binxy_dialog");
}

gboolean binxy_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("binxy_dialog");
	return TRUE;
}

/*************
 * RESAMPLE
 */

/* Idle function for resample */
static gboolean resample_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		gui_function(update_MenuItem, NULL); // WCS not available anymore
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

void on_button_resample_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!confirm_delete_wcs_keywords(gfit))
		return;

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	double sample[2];
	sample[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));
	sample[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));
	int interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_interpolation")));
	gboolean clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_scale_clamp")));

	set_cursor_waiting(TRUE);
	int toX = round_to_int((sample[0] / 100.0) * gfit->rx);
	int toY = round_to_int((sample[1] / 100.0) * gfit->ry);

	// Allocate parameters
	struct resample_args *params = new_resample_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->toX = toX;
	params->toY = toY;
	params->interpolation = interpolation;
	params->clamp = clamp;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_resample_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;  // Resample needs space for transformation
	args->image_hook = resample_image_hook;
	args->mask_hook = resample_mask_hook;
	args->idle_function = resample_idle;
	args->description = _("Resample");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void on_button_resample_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("resample_dialog");
}

gboolean resample_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("resample_dialog");
	return TRUE;
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
	initialize_resample_widgets_if_needed();
	pause_resample_signal_handlers();
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), gfit->rx);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), gfit->ry);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), 100.0);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), 100.0);
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double xvalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_X));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), round_to_int(gfit->rx * xvalue / 100.0));

	if (gtk_toggle_button_get_active(ratio)) {
		double yvalue = xvalue;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), yvalue);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), round_to_int(gfit->ry * yvalue / 100.0));
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double yvalue = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_Y));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), round_to_int(gfit->ry * yvalue / 100.0));

	if (gtk_toggle_button_get_active(ratio)) {
		double xvalue = yvalue;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), xvalue);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px), round_to_int(gfit->rx * xvalue / 100.0));
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_X_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio_button = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double xpix = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_X_px));
	double ratio = xpix / gfit->rx;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_X), ratio * 100.0);

	if (gtk_toggle_button_get_active(ratio_button)) {
		double ypix = round_to_int(gfit->ry * ratio);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), ratio * 100.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px), ypix);
	}
	restart_resample_signal_handlers();
}

void on_spinbutton_resample_Y_px_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	pause_resample_signal_handlers();
	GtkToggleButton *ratio_button = GTK_TOGGLE_BUTTON(button_sample_ratio);
	double ypix = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spinbutton_resample_Y_px));
	double ratio = ypix / gfit->ry;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spinbutton_resample_Y), ratio * 100.0);

	if (gtk_toggle_button_get_active(ratio_button)) {
		double xpix = round_to_int(gfit->rx * ratio);
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

/* Idle function for crop */
static gboolean crop_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	stop_processing_thread();

	if (args->retval == 0) {
		delete_selected_area();
		reset_display_offset();
		update_zoom_label();
		notify_gfit_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}

	free_generic_img_args(args);
	return FALSE;
}

void siril_crop() {
	if (is_preview_active()) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Image backup is active"),
				_("It is impossible to crop the image when the image backup is active. "
				"This occurs when a filter has a live preview active or when a filter with "
				"on-demand ROI preview is open. Please close the filter dialog before "
				"cropping."));
		return;
	}
	clear_stars_list(TRUE);

	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct crop_args *params = new_crop_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	params->area = com.selection;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_crop_args(params);
		set_cursor_waiting(FALSE);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;  // Crop is in-place
	args->image_hook = crop_image_hook_single;
	args->mask_hook = crop_mask_hook;
	args->idle_function = crop_idle;
	args->description = _("Crop");
	args->log_hook = crop_log_hook;
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = 1;  // Crop doesn't benefit from threading

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}
