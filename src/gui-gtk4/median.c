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
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "filters/median.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/median.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"

static GtkWidget *median_roi_preview_btn = NULL;
static GtkWidget *median_amend_note = NULL;

/* Amend mode (convergence C5): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start; Apply/
 * Cancel route through nde_amend_preview_end instead of a worker run + undo
 * save. */
static gboolean median_amend_mode = FALSE;

static void median_dialog_init_statics(void) {
	if (median_roi_preview_btn) return;
	median_roi_preview_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "Median_roi_preview"));
	median_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "median_amend_note"));
}

/* Reverse of the combo-index → kernel-size map in fill_median_params_from_gui. */
static int median_ksize_to_combo(int ksize) {
	switch (ksize) {
		case 3:  return 0;
		case 5:  return 1;
		case 7:  return 2;
		case 9:  return 3;
		case 11: return 4;
		case 13: return 5;
		case 15: return 6;
		default: return 0;
	}
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply builds. */
static void median_prefill_from_amend(void) {
	struct median_filter_data *p = op_desc_median.deserialize(
			nde_amend_preview_params(), nde_amend_preview_op_version());
	if (!p)
		return;
	gtk_drop_down_set_selected(
			GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_ksize_median")),
			median_ksize_to_combo(p->ksize));
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_median")), p->amount);
	gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "median_button_iterations")),
			p->iterations);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
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
	params->fit = gfit;
	params->amount = amount;
	params->iterations = iterations;
}

void median_roi_callback(void) {
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
	if (median_amend_mode) {
		/* Leave amend mode: drop the backup and restore the true pixels
		 * (nothing committed).  Reached via the dialog registry's cancel
		 * hook (Esc / window close). */
		clear_backup();
		nde_amend_preview_end(FALSE, NULL);
		median_amend_mode = FALSE;
		siril_close_dialog("Median_dialog");
		return;
	}
	siril_preview_hide();
	roi_declare_op(NULL);
	remove_roi_callback(median_roi_callback);
	siril_close_dialog("Median_dialog");
}

void on_Median_dialog_show(GtkWidget *widget, gpointer user_data) {
	median_dialog_init_statics();
	gtk_widget_set_visible(median_amend_note, median_amend_mode);

	if (median_amend_mode) {
		/* gfit already shows the pre-record state.  No ROI (the amend
		 * display is the whole pre-record image); the backup is armed with
		 * that state so a ROI preview would still work if drawn. */
		gtk_widget_set_visible(median_roi_preview_btn, FALSE);
		if (gui.roi.active)
			on_clear_roi();
		copy_gfit_to_backup();   /* backup := pre-record state */
		set_notify_block(TRUE);
		median_prefill_from_amend();
		set_notify_block(FALSE);
		return;
	}

	roi_declare_op(&op_desc_median);
	gtk_widget_set_visible(median_roi_preview_btn, gui.roi.active);
	copy_gfit_to_backup();
	add_roi_callback(median_roi_callback);
}

void on_Median_cancel_clicked(GtkButton *button, gpointer user_data) {
	median_close();   /* handles amend mode internally */
}

void on_Median_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (median_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path. */
		struct median_filter_data applied = { 0 };
		fill_median_params_from_gui(&applied, FALSE);
		gchar *blob = op_desc_median.serialize(&applied);
		clear_backup();
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		median_amend_mode = FALSE;
		siril_close_dialog("Median_dialog");
		return;
	}

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

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return;
	}

	/* Always gfit.  A ROI is a preview scope, never an edit scope: the
	 * worker decides on its own whether this run is region-scoped, from
	 * for_preview (set by the dedicated preview button) and the descriptor.
	 *
	 * Apply therefore always reaches the worker with fit == gfit (use_swap),
	 * where the worker itself saves undo, captures the NDE record and logs
	 * through median_log_hook.  The hand-rolled block that used to stand here
	 * existed only for the ROI-scoped apply; on the full-image path it was
	 * duplicating the worker's undo and log. */
	args->fit = params->fit;
	args->op = &op_desc_median;
	args->idle_function = for_preview ? NULL : median_apply_idle;
	args->command_updates_gfit = TRUE;
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;

	generic_image_worker(args);
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void median_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	median_amend_mode = TRUE;
	siril_open_dialog("Median_dialog");
}

gboolean median_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, median_amend_ready, NULL);
	return TRUE;
}
