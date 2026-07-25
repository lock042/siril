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

/* GTK callbacks and preview logic for the SCNR dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/nde_replay.h"
#include "filters/scnr.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/scnr.h"
#include "io/single_image.h"

static GtkToggleButton *scnr_roi_preview = NULL;
static GtkDropDown *scnr_combo = NULL;
static GtkScale *scnr_scale = NULL;
static GtkLabel *scnr_label56 = NULL;
static GtkSpinButton *scnr_spin = NULL;
static GtkCheckButton *scnr_preserve_light = NULL;
static GtkWidget *scnr_amend_note = NULL;

/* Amend mode (convergence C5): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start.  SCNR
 * has no preview machinery (dialogs.c has_preview=FALSE), so there is no
 * backup to arm or clear: the display simply shows the installed pre-record
 * state while the form is open, and Apply/Cancel route through
 * nde_amend_preview_end. */
static gboolean scnr_amend_mode = FALSE;

static void scnr_dialog_init_statics(void) {
	if (scnr_roi_preview) return;
	scnr_roi_preview = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "SCNR_roi_preview"));
	scnr_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_scnr"));
	scnr_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "scale_scnr"));
	scnr_label56 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label56"));
	scnr_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_scnr"));
	scnr_preserve_light = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "preserve_light"));
	scnr_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "scnr_amend_note"));
}

void on_combo_scnr_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply builds. */
static void scnr_prefill_from_amend(void) {
	struct scnr_data *p = op_desc_scnr.deserialize(nde_amend_preview_params(),
	                                               nde_amend_preview_op_version());
	if (!p)
		return;
	gtk_drop_down_set_selected(scnr_combo, p->type);
	gtk_range_set_value(GTK_RANGE(scnr_scale), p->amount);
	siril_toggle_set_active(GTK_WIDGET(scnr_preserve_light), p->preserve);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

static int scnr_process_with_worker(scnr_type type, double amount, gboolean preserve,
                                    gboolean for_preview, gboolean for_roi) {
	struct scnr_data *params = new_scnr_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	params->type = type;
	params->amount = amount;
	params->preserve = preserve;
	params->verbose = !for_preview;
	params->applying = !for_preview;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_scnr_data(params);
		return 1;
	}

	args->fit = for_roi ? &gui.roi.fit : gfit;
	args->op = &op_desc_scnr;
	args->idle_function = NULL;
	args->description = _("Subtractive Chromatic Noise Reduction");
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

static int scnr_update_preview(void) {
	if (siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview))) {
		int type = gtk_drop_down_get_selected(scnr_combo);
		gboolean preserve = siril_toggle_get_active(GTK_WIDGET(scnr_preserve_light));
		double amount = gtk_range_get_value(GTK_RANGE(scnr_scale));
		copy_backup_to_gfit();
		return scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	}
	return 0;
}

void scnr_change_between_roi_and_image(void) {
	gui.roi.operation_supports_roi = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview));
	notify_update((gpointer) param);
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	scnr_dialog_init_statics();
	gtk_widget_set_visible(scnr_amend_note, scnr_amend_mode);

	if (scnr_amend_mode) {
		/* gfit already shows the pre-record state.  SCNR has no preview,
		 * so there is no backup and no ROI preview in amend mode: hide the
		 * ROI preview button and clear any ROI (the amend display is the
		 * whole pre-record image). */
		gtk_widget_set_visible(GTK_WIDGET(scnr_roi_preview), FALSE);
		if (gui.roi.active)
			on_clear_roi();
		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(scnr_roi_preview), FALSE);
		scnr_prefill_from_amend();
		set_notify_block(FALSE);
		/* Enable the amount widgets according to the prefilled type. */
		on_combo_scnr_changed(G_OBJECT(scnr_combo), NULL, NULL);
		return;
	}

	roi_supported(TRUE);
	gtk_widget_set_visible(GTK_WIDGET(scnr_roi_preview), gui.roi.active);

	copy_gfit_to_backup();
	add_roi_callback(scnr_change_between_roi_and_image);

	if (gui.roi.active) {
		scnr_change_between_roi_and_image();
	}

	int type = gtk_drop_down_get_selected(scnr_combo);
	if (type == -1)
		gtk_drop_down_set_selected(scnr_combo, 0);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	if (scnr_amend_mode) {
		/* Leave amend mode: restore the true pixels (nothing committed).
		 * No backup/preview machinery to tear down. */
		nde_amend_preview_end(FALSE, NULL);
		scnr_amend_mode = FALSE;
		siril_close_dialog("SCNR_dialog");
		return;
	}
	roi_supported(FALSE);

	if (siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview))) {
		copy_backup_to_gfit();
		gfit_modified_update_gui();
	}

	clear_backup();
	remove_roi_callback(scnr_change_between_roi_and_image);
	siril_close_dialog("SCNR_dialog");
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (scnr_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path.
		 * preview_end restores the true pixels first, then replays the
		 * edited history from the pre-record state this session deposited. */
		struct scnr_data applied = { 0 };
		applied.type = gtk_drop_down_get_selected(scnr_combo);
		applied.amount = gtk_range_get_value(GTK_RANGE(scnr_scale));
		applied.preserve = siril_toggle_get_active(GTK_WIDGET(scnr_preserve_light));
		gchar *blob = op_desc_scnr.serialize(&applied);
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		scnr_amend_mode = FALSE;
		siril_close_dialog("SCNR_dialog");
		return;
	}

	int type = gtk_drop_down_get_selected(scnr_combo);
	gboolean preserve = siril_toggle_get_active(GTK_WIDGET(scnr_preserve_light));
	double amount = gtk_range_get_value(GTK_RANGE(scnr_scale));

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	gboolean is_preview = ((GtkWidget*) button == GTK_WIDGET(scnr_roi_preview));

	if (is_preview) {
		scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	} else {
		if (siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview))) {
			copy_backup_to_gfit();
		}

		set_cursor_waiting(TRUE);
		scnr_process_with_worker(type, amount, preserve, FALSE, FALSE);

		clear_backup();
		remove_roi_callback(scnr_change_between_roi_and_image);
		roi_supported(FALSE);
		siril_close_dialog("SCNR_dialog");
	}
}

void on_combo_scnr_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	int type = gtk_drop_down_get_selected(scnr_combo);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_label56), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_spin), type > 1);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview));
	notify_update((gpointer) param);
}

void on_SCNR_parameter_changed(GtkWidget *widget, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(scnr_roi_preview));
	notify_update((gpointer) param);
}

void on_SCNR_roi_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(button))) {
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = scnr_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void scnr_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	scnr_amend_mode = TRUE;
	siril_open_dialog("SCNR_dialog");
}

void scnr_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, scnr_amend_ready, NULL);
}
