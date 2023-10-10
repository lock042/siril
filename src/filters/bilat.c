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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "bilat.h"

static float bilat_d_value = 0.0f, bilat_sigma_col_value = 0.0f, bilat_sigma_spatial_value = 0.0f;

static int bilat_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))))
		copy_backup_to_gfit();
	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;
	bilateral(fit, bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value, FALSE);
	notify_gfit_modified();
	return 0;
}

void bilat_change_between_roi_and_image() {
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

static void bilat_startup() {
	add_roi_callback(bilat_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
}

static void bilat_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Bilateral filter: (d=%2.2f, sigma_col=%2.2f, sigma_spatial=%2.2f)"),
				bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value);
	}
	backup_roi();
	roi_supported(FALSE);
	remove_roi_callback(bilat_change_between_roi_and_image);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int bilat_process_all() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))))
		copy_backup_to_gfit();
	bilateral(&gfit, bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value, FALSE);
	populate_roi();
	notify_gfit_modified();
	return 0;
}

static void apply_bilat_changes() {
	gboolean status = (bilat_sigma_col_value != 0.0f) || (bilat_sigma_spatial_value != 0.0f);
	bilat_close(!status);
}

void apply_bilat_cancel() {
	bilat_close(TRUE);
	siril_close_dialog("bilat_dialog");
}

/*** callbacks **/

void on_bilat_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_d"));
	GtkSpinButton *spin_bilat_sigma_col = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_col"));
	GtkSpinButton *spin_bilat_sigma_spatial = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_spatial"));

	if (gui.rendering_mode == LINEAR_DISPLAY)
		setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535

	bilat_startup();
	bilat_d_value = 0.0f;
	bilat_sigma_col_value = 0.0f;
	bilat_sigma_spatial_value = 0.0f;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_d, bilat_d_value);
	gtk_spin_button_set_value(spin_bilat_sigma_spatial, bilat_sigma_spatial_value);
	gtk_spin_button_set_value(spin_bilat_sigma_col, bilat_sigma_col_value);
	set_notify_block(FALSE);

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_bilat_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_bilat_cancel();
}

void on_bilat_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))) || gui.roi.active) {
		bilat_process_all();
	}

	apply_bilat_changes();
	siril_close_dialog("bilat_dialog");
}

void on_bilat_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_bilat_changes();
}

void on_bilat_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_d"));
	GtkSpinButton *spin_bilat_sigma_col = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_col"));
	GtkSpinButton *spin_bilat_sigma_spatial = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_spatial"));
	bilat_d_value = 0.0;
	bilat_sigma_col_value = 0.0;
	bilat_sigma_spatial_value = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_d, bilat_d_value);
	gtk_spin_button_set_value(spin_bilat_sigma_col, bilat_sigma_col_value);
	gtk_spin_button_set_value(spin_bilat_sigma_spatial, bilat_sigma_spatial_value);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_spin_bilat_d_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_d_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_spin_bilat_sigma_col_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_sigma_col_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_spin_bilat_sigma_spatial_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_sigma_spatial_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_bilat_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")))) {
		copy_backup_to_gfit();
		redraw(REMAP_ALL);
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = bilat_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
