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
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "core/nde/nde_replay.h"
#include "filters/clahe.h"
#include "gui-gtk4/clahe.h"
#include "gui-gtk4/nde_editors.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"

static GtkSpinButton *clahe_spin = NULL, *clahe_tile_spin = NULL;
static GtkCheckButton *clahe_preview_btn = NULL;
static GtkWidget *clahe_amend_note = NULL;

/* Amend mode (convergence C5b): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start;
 * previews run against it through the ordinary backup machinery, and Apply/
 * Cancel route through nde_amend_preview_end instead of a worker run + undo
 * save. */
static gboolean clahe_amend_mode = FALSE;

static void clahe_dialog_init_statics(void) {
	if (clahe_spin) return;
	clahe_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_clahe"));
	clahe_tile_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "clahe_tiles_size_spin"));
	clahe_preview_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "clahe_preview"));
	clahe_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "clahe_amend_note"));
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply builds. */
static void clahe_prefill_from_amend(void) {
	clahe_params *p = op_desc_clahe.deserialize(nde_amend_preview_params(),
	                                            nde_amend_preview_op_version());
	if (!p)
		return;
	gtk_spin_button_set_value(clahe_spin, p->clip);
	gtk_spin_button_set_value(clahe_tile_spin, p->tileSize);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

static void get_clahe_values(double *clip, int *tileSize) {
	if (clip)
		*clip = gtk_spin_button_get_value(clahe_spin);
	if (tileSize)
		*tileSize = gtk_spin_button_get_value(clahe_tile_spin);
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
	args->op = &op_desc_clahe;
	args->idle_function = NULL;
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;

	start_in_new_thread(generic_image_worker, args);
	return 0;
}

static int clahe_update_preview(void) {
	if (siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn))) {
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
	if (clahe_amend_mode) {
		/* Leave amend mode: drop the pre-record backup and restore the true
		 * pixels (nothing was committed). */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		clear_backup();
		nde_amend_preview_end(FALSE, NULL);
		clahe_amend_mode = FALSE;
		siril_close_dialog("CLAHE_dialog");
		return;
	}
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
}

gboolean on_clahe_cancel_clicked(GtkWidget *menuitem, gpointer user_data) {
	apply_clahe_cancel();   /* handles amend mode internally */
	return FALSE;
}

/* close-request handler: return FALSE so the window proceeds to close
 * after we've committed pending changes. */
gboolean on_CLAHE_dialog_close(GtkWindow *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	apply_clahe_cancel();   /* handles amend mode internally */
	return FALSE;
}

void on_clahe_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	if (clahe_amend_mode) {
		/* Reset means "back to the record's recorded parameters" here. */
		clahe_prefill_from_amend();
	} else {
		gtk_spin_button_set_value(clahe_spin, 2);
		gtk_spin_button_set_value(clahe_tile_spin, 8);
	}
	siril_toggle_set_active(GTK_WIDGET(clahe_preview_btn), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn));
	notify_update((gpointer) param);
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (clahe_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path. */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		clahe_params applied = { 0 };
		get_clahe_values(&applied.clip, &applied.tileSize);
		gchar *blob = op_desc_clahe.serialize(&applied);
		clear_backup();
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		clahe_amend_mode = FALSE;
		siril_close_dialog("CLAHE_dialog");
		return;
	}

	if (!check_ok_if_cfa())
		return;

	if (siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn))) {
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
	args->op = &op_desc_clahe;
	args->idle_function = NULL;
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;

	start_in_new_thread(generic_image_worker, args);

	clahe_close(FALSE);
	siril_close_dialog("CLAHE_dialog");
}

void on_CLAHE_dialog_show(GtkWidget *widget, gpointer user_data) {
	clahe_dialog_init_statics();
	nde_amend_note_update(clahe_amend_note, clahe_amend_mode,
	                      NULL);

	if (clahe_amend_mode) {
		/* gfit already shows the pre-record state.  No ICC/ROI juggling. */
		if (gui.roi.active)
			on_clear_roi();
		copy_gfit_to_backup();   /* backup := pre-record state */

		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(clahe_preview_btn), TRUE);
		clahe_prefill_from_amend();
		set_notify_block(FALSE);

		/* One preview tick so the dialog opens showing the record applied
		 * with its current parameters. */
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = clahe_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
		return;
	}

	clahe_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(clahe_tile_spin, 8);
	gtk_spin_button_set_value(clahe_spin, 2.0);
	set_notify_block(FALSE);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn));
	notify_update((gpointer) param);
}

void on_spin_clahe_value_changed(GtkSpinButton *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn));
	notify_update((gpointer) param);
}

void on_clahe_tiles_size_spin_value_changed(GtkSpinButton *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(clahe_preview_btn));
	notify_update((gpointer) param);
}

void on_clahe_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(button))) {
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

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void clahe_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	clahe_amend_mode = TRUE;
	siril_open_dialog("CLAHE_dialog");
}

gboolean clahe_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, clahe_amend_ready, NULL);
	return TRUE;
}
