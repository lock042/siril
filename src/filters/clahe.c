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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/siril_log.h"
#include "opencv/opencv.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "gui/siril_preview.h"

#include "clahe.h"

/* The actual CLAHE processing hook */
int clahe_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	clahe_params *params = (clahe_params *)args->user;
	if (!params)
		return 1;

	return cvClahe(fit, params->clip, params->tileSize);
}

gchar *clahe_log_hook(gpointer p, log_hook_detail detail) {
	clahe_params *params = (clahe_params *) p;
	if (!params) return NULL;
	gchar *message = g_strdup_printf(_("CLAHE (size=%d, clip=%.2f)"), params->tileSize, params->clip);
	return message;
}

/* Idle function for preview / final application updates */
static gboolean clahe_worker_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();

	if (args->retval == 0) {
		notify_gfit_modified();
	}

	free_generic_img_args(args);
	return FALSE;
}

/* Helper function to get current widget values */
static void get_clahe_values(double *clip, int *tileSize) {
	if (clip)
		*clip = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")));
	if (tileSize)
		*tileSize = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")));
}

/* Create and launch CLAHE processing */
static int clahe_process_with_worker(gboolean for_preview) {
	// Allocate parameters
	clahe_params *params = calloc(1, sizeof(clahe_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Get current values from widgets
	get_clahe_values(&params->clip, &params->tileSize);

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f; // CLAHE needs additional memory for tile processing
	args->image_hook = clahe_image_hook;
	args->idle_function = clahe_worker_idle;
	args->description = _("CLAHE");
	args->verbose = !for_preview; // Only verbose for final application
	args->user = params;
	args->log_hook = clahe_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = FALSE; // CLAHE always operates on full image (no ROI support)

	start_in_new_thread(generic_image_worker, args);
	return 0;
}

/* Update preview using the worker */
static int clahe_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")))) {
		copy_backup_to_gfit();
		return clahe_process_with_worker(TRUE);
	}
	return 0;
}

static void clahe_startup() {
	copy_gfit_to_backup();
}

static void clahe_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(gfit);

		double clip;
		int tileSize;
		get_clahe_values(&clip, &tileSize);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

void apply_clahe_cancel() {
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
}

/** callbacks **/

gboolean on_clahe_cancel_clicked(GtkMenuItem *menuitem, gpointer user_data) {
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
	return FALSE;
}

void on_CLAHE_dialog_close(GtkDialog *dialog, gpointer user_data) {
	clahe_close(TRUE);
}

void on_clahe_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")), 2);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")), 8);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	// If preview is on, need to copy backup to gfit first
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")))) {
		copy_backup_to_gfit();
	}

	// Allocate parameters for full image processing
	clahe_params *params = calloc(1, sizeof(clahe_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return;
	}

	// Get current values from widgets
	get_clahe_values(&params->clip, &params->tileSize);

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return;
	}

	// Always process the full gfit
	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = clahe_image_hook;
	args->idle_function = clahe_worker_idle;
	args->description = _("CLAHE");
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = clahe_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	start_in_new_thread(generic_image_worker, args);

	clahe_close(FALSE);
	siril_close_dialog("CLAHE_dialog");
}

void on_CLAHE_dialog_show(GtkWidget *widget, gpointer user_data) {
	clahe_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")), 8);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")), 2.0);
	set_notify_block(FALSE);

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

/** adjusters **/
void on_spin_clahe_value_changed(GtkSpinButton *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

void on_clahe_tiles_size_spin_value_changed(GtkSpinButton *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

void on_clahe_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = clahe_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
