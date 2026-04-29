/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

/* GTK callbacks and preview logic for the CLAHE dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "filters/clahe.h"
#include "gui/clahe.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"

static void get_clahe_values(double *clip, int *tileSize) {
	if (clip)
		*clip = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")));
	if (tileSize)
		*tileSize = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")));
}

static int clahe_process_with_worker(gboolean for_preview) {
	clahe_params *params = calloc(1, sizeof(clahe_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	get_clahe_values(&params->clip, &params->tileSize);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = clahe_image_hook;
	args->idle_function = NULL;
	args->description = _("CLAHE");
	args->verbose = !for_preview;
	args->user = params;
	args->log_hook = clahe_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = FALSE;

	start_in_new_thread(generic_image_worker, args);
	return 0;
}

static int clahe_update_preview(void) {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")))) {
		copy_backup_to_gfit();
		return clahe_process_with_worker(TRUE);
	}
	return 0;
}

static void clahe_startup(void) {
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

void apply_clahe_cancel(void) {
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
}

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

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")))) {
		copy_backup_to_gfit();
	}

	clahe_params *params = calloc(1, sizeof(clahe_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return;
	}

	get_clahe_values(&params->clip, &params->tileSize);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = clahe_image_hook;
	args->idle_function = NULL;
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

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));
	notify_update((gpointer) param);
}

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
