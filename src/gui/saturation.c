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

/* GTK callbacks and preview logic for the saturation dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "filters/saturation.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/saturation.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"

static GtkToggleButton *satu_preview_btn = NULL;
static GtkComboBox *satu_combo = NULL;
static GtkRange *satu_scale = NULL, *satu_scale_bkg = NULL;

static void satu_dialog_init_statics(void) {
	if (satu_preview_btn) return;
	satu_preview_btn = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "satu_preview"));
	satu_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_saturation"));
	satu_scale = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_satu"));
	satu_scale_bkg = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_satu_bkg"));
}

static double satu_amount, background_factor;
static int satu_hue_type;
static gboolean satu_show_preview;

static int satu_update_preview(void);

static gboolean satu_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		gfit_modified_update_gui();
	}
	free_generic_img_args(args);
	clear_backup();
	return FALSE;
}

static int satu_process_with_worker(gboolean for_preview) {
	saturation_params *params = calloc(1, sizeof(saturation_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	params->coeff = satu_amount;
	params->background_factor = background_factor;
	satu_set_hues_from_types(params, satu_hue_type);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	args->fit = gui.roi.active ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = saturation_image_hook;
	args->idle_function = for_preview ? NULL : satu_apply_idle;
	args->description = _("Saturation");
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;
	args->mask_aware = TRUE;
	if (!for_preview) {
		args->log_hook = satu_log_hook;
	}

	if (for_preview)
		generic_image_worker(args);
	else
		start_in_new_thread(generic_image_worker, args);
	return 0;
}

static int satu_update_preview(void) {
	if (gtk_toggle_button_get_active(satu_preview_btn)) {
		copy_backup_to_gfit();
		return satu_process_with_worker(TRUE);
	}
	return 0;
}

void satu_change_between_roi_and_image(void) {
	gui.roi.operation_supports_roi = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

static void satu_startup(void) {
	roi_supported(TRUE);
	add_roi_callback(satu_change_between_roi_and_image);
	copy_gfit_to_backup();
	satu_amount = 0.0;
	satu_hue_type = 6;
}

static void satu_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		if (satu_amount != 0.0) {
			copy_backup_to_gfit();
			notify_gfit_data_modified();
			gfit_modified_update_gui();
		}
	}
	roi_supported(FALSE);
	remove_roi_callback(satu_change_between_roi_and_image);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static void apply_satu_changes(void) {
	gboolean status = satu_amount != 0.0;
	satu_close(!status);
}

gboolean on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
	return FALSE;
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(satu_preview_btn))
		copy_backup_to_gfit();

	satu_process_with_worker(FALSE);
	siril_close_dialog("satu_dialog");
}

void on_satu_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_satu_changes();
}

void on_satu_dialog_show(GtkWidget *widget, gpointer user_data) {
	satu_dialog_init_statics();
	satu_startup();
	satu_amount = 0.0;
	satu_hue_type = 6;
	background_factor = 1.0;

	set_notify_block(TRUE);
	gtk_combo_box_set_active(satu_combo, satu_hue_type);
	gtk_range_set_value(satu_scale, satu_amount);
	gtk_range_set_value(satu_scale_bkg, background_factor);
	set_notify_block(FALSE);

	satu_show_preview = gtk_toggle_button_get_active(satu_preview_btn);
}

void on_combo_saturation_changed(GtkComboBox* box, gpointer user_data) {
	satu_hue_type = gtk_combo_box_get_active(box);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_undo_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	double prev_satu = satu_amount;

	set_notify_block(TRUE);
	gtk_range_set_value(satu_scale, 0);
	gtk_range_set_value(satu_scale_bkg, 1);
	gtk_toggle_button_set_active(satu_preview_btn, TRUE);
	gtk_combo_box_set_active(satu_combo, 6);
	set_notify_block(FALSE);

	if (prev_satu != 0.0) {
		copy_backup_to_gfit();
		gfit_modified_update_gui();
		set_cursor_waiting(FALSE);
	}
}

void apply_satu_cancel(void) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
}

void on_spin_satu_value_changed(GtkSpinButton *button, gpointer user_data) {
	satu_amount = gtk_spin_button_get_value(button);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_spin_satu_bkg_value_changed(GtkSpinButton *button, gpointer user_data) {
	background_factor = gtk_spin_button_get_value(button);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	satu_show_preview = gtk_toggle_button_get_active(button);
	if (!satu_show_preview) {
		cancel_and_wait_for_preview();
		copy_backup_to_gfit();
		notify_gfit_data_modified();
		redraw(REDRAW_ALL);
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = satu_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
