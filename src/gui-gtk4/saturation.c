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
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "core/nde_replay.h"
#include "filters/saturation.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/saturation.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"

static GtkCheckButton *satu_preview_btn = NULL;
static GtkDropDown *satu_combo = NULL;
static GtkRange *satu_scale = NULL, *satu_scale_bkg = NULL;
static GtkWidget *satu_amend_note = NULL;

/* Amend mode (convergence C5): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start;
 * previews run against it through the ordinary backup machinery, and Apply/
 * Cancel route through nde_amend_preview_end instead of a worker run + undo
 * save. */
static gboolean satu_amend_mode = FALSE;

static void satu_dialog_init_statics(void) {
	if (satu_preview_btn) return;
	satu_preview_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "satu_preview"));
	satu_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_saturation"));
	satu_scale = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_satu"));
	satu_scale_bkg = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_satu_bkg"));
	satu_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "satu_amend_note"));
}

static double satu_amount, background_factor;
static int satu_hue_type;
static gboolean satu_show_preview;

/* Reverse of satu_set_hues_from_types(): map a recorded (h_min,h_max) degree
 * range back to the hue-type combo index.  Falls back to Global (6) for any
 * range that does not match a preset (e.g. a script-set custom range). */
static int satu_hue_type_from_hues(double h_min, double h_max) {
	for (int t = 0; t <= 6; t++) {
		saturation_params probe = { 0 };
		satu_set_hues_from_types(&probe, t);
		if (probe.h_min == h_min && probe.h_max == h_max)
			return t;
	}
	return 6;
}

/* Fill the widgets/statics from the amended record's current parameters
 * through the op's own deserializer — the same fields the normal apply reads. */
static void satu_prefill_from_amend(void) {
	saturation_params *p = op_desc_saturation.deserialize(
			nde_amend_preview_params(), nde_amend_preview_op_version());
	if (!p)
		return;
	satu_amount = p->coeff;
	background_factor = p->background_factor;
	satu_hue_type = satu_hue_type_from_hues(p->h_min, p->h_max);
	gtk_drop_down_set_selected(satu_combo, satu_hue_type);
	gtk_range_set_value(satu_scale, satu_amount);
	gtk_range_set_value(satu_scale_bkg, background_factor);
	free(p);
}

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
	args->op = &op_desc_saturation;
	args->idle_function = for_preview ? NULL : satu_apply_idle;
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;
	args->mask_aware = TRUE;

	if (for_preview)
		generic_image_worker(args);
	else
		start_in_new_thread(generic_image_worker, args);
	return 0;
}

static int satu_update_preview(void) {
	if (siril_toggle_get_active(GTK_WIDGET(satu_preview_btn))) {
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

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void satu_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	satu_amend_mode = TRUE;
	siril_open_dialog("satu_dialog");
}

void satu_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, satu_amend_ready, NULL);
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

/* Leave amend mode: drop the backup, restore the true pixels (nothing
 * committed) and close.  Returns TRUE when it handled the amend case. */
static gboolean satu_amend_cancel(void) {
	if (!satu_amend_mode)
		return FALSE;
	cancel_pending_update();
	cancel_and_wait_for_preview();
	clear_backup();
	nde_amend_preview_end(FALSE, NULL);
	satu_amend_mode = FALSE;
	siril_close_dialog("satu_dialog");
	return TRUE;
}

gboolean on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	if (satu_amend_cancel())
		return FALSE;
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
	return FALSE;
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	if (satu_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses (satu_set_hues_from_types maps the hue-type
		 * combo to h_min/h_max exactly as satu_process_with_worker does),
		 * then route it through the amend path. */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		saturation_params applied = { 0 };
		applied.coeff = satu_amount;
		applied.background_factor = background_factor;
		satu_set_hues_from_types(&applied, satu_hue_type);
		gchar *blob = op_desc_saturation.serialize(&applied);
		clear_backup();
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		satu_amend_mode = FALSE;
		siril_close_dialog("satu_dialog");
		return;
	}

	if (siril_toggle_get_active(GTK_WIDGET(satu_preview_btn)))
		copy_backup_to_gfit();

	satu_process_with_worker(FALSE);
	siril_close_dialog("satu_dialog");
}

/* close-request handler: return FALSE so the window proceeds to close
 * after we've committed any pending changes. */
gboolean on_satu_dialog_close(GtkWindow *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	if (satu_amend_cancel())
		return FALSE;
	apply_satu_changes();
	return FALSE;
}

void on_satu_dialog_show(GtkWidget *widget, gpointer user_data) {
	satu_dialog_init_statics();
	gtk_widget_set_visible(satu_amend_note, satu_amend_mode);

	if (satu_amend_mode) {
		/* gfit already shows the pre-record state.  No ROI (previews are
		 * full-image against the pre-record backup). */
		if (gui.roi.active)
			on_clear_roi();
		copy_gfit_to_backup();   /* backup := pre-record state */

		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(satu_preview_btn), TRUE);
		satu_prefill_from_amend();
		set_notify_block(FALSE);
		satu_show_preview = TRUE;

		/* One preview tick so the dialog opens showing the record applied
		 * with its current parameters. */
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = satu_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
		return;
	}

	satu_startup();
	satu_amount = 0.0;
	satu_hue_type = 6;
	background_factor = 1.0;

	set_notify_block(TRUE);
	gtk_drop_down_set_selected(satu_combo, satu_hue_type);
	gtk_range_set_value(satu_scale, satu_amount);
	gtk_range_set_value(satu_scale_bkg, background_factor);
	set_notify_block(FALSE);

	satu_show_preview = siril_toggle_get_active(GTK_WIDGET(satu_preview_btn));
}

void on_combo_saturation_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	satu_hue_type = gtk_drop_down_get_selected(box);

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
	siril_toggle_set_active(GTK_WIDGET(satu_preview_btn), TRUE);
	gtk_drop_down_set_selected(satu_combo, 6);
	set_notify_block(FALSE);

	if (prev_satu != 0.0) {
		copy_backup_to_gfit();
		gfit_modified_update_gui();
		set_cursor_waiting(FALSE);
	}
}

void apply_satu_cancel(void) {
	if (satu_amend_cancel())
		return;
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

void on_satu_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	satu_show_preview = siril_toggle_get_active(GTK_WIDGET(button));
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
