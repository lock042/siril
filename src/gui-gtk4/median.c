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

/* GTK callbacks, ROI callback and idle function for the Median Filter dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "filters/median.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/median.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"

static GtkWidget *median_roi_preview_btn = NULL;

static void median_dialog_init_statics(void) {
	if (median_roi_preview_btn) return;
	median_roi_preview_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "Median_roi_preview"));
}

static void fill_median_params_from_gui(struct median_filter_data *params, gboolean for_preview) {
	if (!params)
		return;

	int combo_size = gtk_drop_down_get_selected(
			GTK_DROP_DOWN(
				gtk_builder_get_object(gui.builder, "combo_ksize_median")));
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_median")));
	int iterations = round_to_int(gtk_spin_button_get_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "median_button_iterations"))));

	switch (combo_size) {
		default:
		case 0:
			params->ksize = 3;
			break;
		case 1:
			params->ksize = 5;
			break;
		case 2:
			params->ksize = 7;
			break;
		case 3:
			params->ksize = 9;
			break;
		case 4:
			params->ksize = 11;
			break;
		case 5:
			params->ksize = 13;
			break;
		case 6:
			params->ksize = 15;
			break;
	}
	params->fit = for_preview && gui.roi.active && !com.headless ? &gui.roi.fit : gfit;
	params->amount = amount;
	params->iterations = iterations;
}

void median_roi_callback(void) {
	gui.roi.operation_supports_roi = TRUE;
	gtk_widget_set_visible(median_roi_preview_btn, gui.roi.active);
	copy_backup_to_gfit();
	gfit_modified_update_gui();
}

static gboolean median_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		copy_gfit_to_backup();
		gfit_modified_update_gui();
	}
	free_generic_img_args(args);
	median_close();
	return FALSE;
}

void median_close(void) {
	siril_preview_hide();
	roi_supported(FALSE);
	remove_roi_callback(median_roi_callback);
	siril_close_dialog("Median_dialog");
}

void on_Median_dialog_show(GtkWidget *widget, gpointer user_data) {
	median_dialog_init_statics();
	roi_supported(TRUE);
	gtk_widget_set_visible(median_roi_preview_btn, gui.roi.active);
	copy_gfit_to_backup();
	add_roi_callback(median_roi_callback);
}

void on_Median_cancel_clicked(GtkButton *button, gpointer user_data) {
	median_close();
}

void on_Median_Apply_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!check_ok_if_cfa())
		return;
	set_cursor_waiting(TRUE);

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	copy_backup_to_gfit();

	struct median_filter_data *params = calloc(1, sizeof(struct median_filter_data));
	gboolean for_preview = ((GtkWidget*) button == median_roi_preview_btn);
	fill_median_params_from_gui(params, for_preview);

	if (!for_preview && !com.script) {
		undo_save_state(gfit, _("Median Filter (filter=%dx%d px, iters=%d), mod=%.3lf"),
			params->ksize, params->ksize, params->iterations, params->amount);
		siril_log_info(_("Median Filter (filter=%dx%d px, iterations=%d, modulation=%.3lf)\n"),
			params->ksize, params->ksize, params->iterations, params->amount);
	}

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return;
	}

	args->fit = gui.roi.active ? &gui.roi.fit : gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = median_image_hook;
	args->idle_function = for_preview ? NULL : median_apply_idle;
	args->description = _("Median filter");
	args->command_updates_gfit = TRUE;
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;

	generic_image_worker(args);
}
