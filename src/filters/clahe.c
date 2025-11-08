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
#include "opencv/opencv.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "gui/siril_preview.h"

#include "clahe.h"

static double clahe_limit_value;
static int clahe_tile_size;
static gboolean clahe_show_preview;

static void clahe_startup() {
	copy_gfit_to_backup();
	clahe_limit_value = 2.0;
	clahe_tile_size = 8;
}

static void clahe_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("CLAHE (size=%d, clip=%.2f)"), clahe_tile_size, clahe_limit_value);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int clahe_update_preview() {
	copy_backup_to_gfit();

	struct CLAHE_data *args = calloc(1, sizeof(struct CLAHE_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->clip = clahe_limit_value;
	args->tileSize = clahe_tile_size;

	if (!start_in_new_thread(clahe, args))
		free(args);
	return 0;
}

static gboolean end_clahe(gpointer p) {
	struct CLAHE_data *args = (struct CLAHE_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);

	free(args);

	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	return FALSE;
}

gpointer clahe(gpointer p) {
	struct CLAHE_data *args = (struct CLAHE_data*) p;

	cvClahe(args->fit, args->clip, args->tileSize);

	siril_add_idle(end_clahe, args);
	return GINT_TO_POINTER(0);
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
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (clahe_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = clahe_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	clahe_close(FALSE);
	siril_close_dialog("CLAHE_dialog");
}

void on_CLAHE_dialog_show(GtkWidget *widget, gpointer user_data) {
	clahe_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")), clahe_tile_size);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")),	clahe_limit_value);
	set_notify_block(FALSE);

	clahe_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

/** adjusters **/
void on_spin_clahe_value_changed(GtkSpinButton *button, gpointer user_data) {
	clahe_limit_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_tiles_size_spin_value_changed(GtkSpinButton *button, gpointer user_data) {
	clahe_tile_size = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	clahe_show_preview = gtk_toggle_button_get_active(button);
	if (!clahe_show_preview) {
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

