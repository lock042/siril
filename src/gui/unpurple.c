/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"
#include "opencv/opencv.h"
#include "filters/synthstar.h"
#include "filters/unpurple.h"

static double mod_b = 1.0, thresh = 0.0, old_thresh = 0.0;
static fits starmask = {0};
static gboolean is_roi = FALSE;

/* Helper function to get current widget values */
static void get_unpurple_values(double *mod_b_out, double *thresh_out, gboolean *withstarmask_out) {
	if (mod_b_out)
		*mod_b_out = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_mod_b")));
	if (thresh_out)
		*thresh_out = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_thresh")));
	if (withstarmask_out)
		*withstarmask_out = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_stars")));
}

/* Create and launch unpurple processing */
static int unpurple_process_with_worker(gboolean for_preview, gboolean for_roi) {
	// Allocate parameters
	struct unpurpleargs *params = new_unpurple_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Get current values from widgets
	gboolean withstarmask;
	get_unpurple_values(&params->mod_b, &params->thresh, &withstarmask);
	params->withstarmask = withstarmask;
	params->fit = for_roi ? &gui.roi.fit : gfit;

	// Set up starmask if needed
	params->starmask_needs_freeing = FALSE;
	if (withstarmask) {
		// For preview, reuse existing starmask if valid
		if (for_preview && starmask.naxis != 0 && gui.roi.active == is_roi && old_thresh == params->thresh) {
			params->starmask = &starmask;
		} else {
			// Need to create new starmask
			params->starmask = calloc(1, sizeof(fits));
			if (!params->starmask) {
				PRINT_ALLOC_ERR;
				free_unpurple_args(params);
				return 1;
			}

			if (generate_binary_starmask(params->fit, &params->starmask, params->thresh)) {
				free_unpurple_args(params);
				return 1;
			}
			params->starmask_needs_freeing = TRUE;

			// Update static starmask for preview reuse
			if (for_preview) {
				clearfits(&starmask);
				copyfits(params->starmask, &starmask, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
				is_roi = gui.roi.active;
				old_thresh = params->thresh;
			}
		}
	} else {
		params->starmask = NULL;
	}

	params->verbose = !for_preview;
	params->applying = !for_preview;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_unpurple_args(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = for_roi ? &gui.roi.fit : gfit;
	args->mem_ratio = 2.0f; // unpurple needs some extra memory
	args->image_hook = unpurple_image_hook;
	args->idle_function = for_preview ? unpurple_preview_idle : unpurple_apply_idle;
	args->description = _("Unpurple Filter");
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = for_roi;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}
	return 0;
}

static int unpurple_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")))) {
		copy_backup_to_gfit();
		return unpurple_process_with_worker(TRUE, gui.roi.active);
	}
	return 0;
}

void unpurple_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

static void unpurple_startup() {
	copy_gfit_to_backup();
	add_roi_callback(unpurple_change_between_roi_and_image);
	roi_supported(TRUE);
}

static void unpurple_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copy_backup_to_gfit();
		notify_gfit_modified();
	} else {
		invalidate_stats_from_fit(gfit);

		double mod_b_val, thresh_val;
		get_unpurple_values(&mod_b_val, &thresh_val, NULL);

		undo_save_state(get_preview_gfit_backup(),
				_("Unpurple filter: (thresh=%2.2lf, mod_b=%2.2lf)"),
				thresh_val, mod_b_val);
	}
	roi_supported(FALSE);
	remove_roi_callback(unpurple_change_between_roi_and_image);
	clearfits(&starmask);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static void apply_unpurple_changes() {
	double mod_b_val, thresh_val;
	get_unpurple_values(&mod_b_val, &thresh_val, NULL);
	gboolean status = (mod_b_val != 1.0) || (thresh_val != 0.0);
	unpurple_close(!status);
}

void apply_unpurple_cancel() {
	unpurple_close(TRUE);
	siril_close_dialog("unpurple_dialog");
}

/*** callbacks **/

void on_unpurple_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_unpurple_mod_b = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_mod_b"));
	GtkSpinButton *spin_unpurple_thresh = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_thresh"));

	unpurple_startup();
	clearfits(&starmask);

	mod_b = 1.0;
	thresh = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_unpurple_mod_b, mod_b);
	gtk_spin_button_set_value(spin_unpurple_thresh, thresh);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")), FALSE);
	set_notify_block(FALSE);
}

void on_unpurple_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_unpurple_cancel();
}

void on_unpurple_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	// If preview is on, need to copy backup to gfit first
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")))) {
		copy_backup_to_gfit();
	}

	// Always process full image when Apply is clicked
	unpurple_process_with_worker(FALSE, FALSE);

	apply_unpurple_changes();
	siril_close_dialog("unpurple_dialog");
}

void on_unpurple_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_unpurple_changes();
}

void on_unpurple_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_unpurple_mod_b = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_mod_b"));
	GtkSpinButton *spin_unpurple_thresh = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_thresh"));

	mod_b = 1.0;
	thresh = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_unpurple_mod_b, mod_b);
	gtk_spin_button_set_value(spin_unpurple_thresh, thresh);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

/*** adjusters **/
void on_spin_unpurple_mod_b_value_changed(GtkSpinButton *button, gpointer user_data) {
	mod_b = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

void on_spin_unpurple_thresh_value_changed(GtkSpinButton *button, gpointer user_data) {
	thresh = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

void on_unpurple_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")))) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = unpurple_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer)param);
	}
}

void on_unpurple_stars_toggled(GtkToggleButton *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}
