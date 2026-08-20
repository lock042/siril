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

/* GTK callbacks and idle function for the Canon Banding Reduction dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/nde/nde_replay.h"
#include "filters/banding.h"
#include "gui-gtk4/banding.h"
#include "gui-gtk4/nde_editors.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

static GtkRange *banding_scale_amount = NULL, *banding_scale_invsigma = NULL;
static GtkCheckButton *banding_protect_highlights = NULL, *banding_vertical = NULL, *banding_seq = NULL;
static GtkEntry *banding_seq_entry = NULL;
static GtkWidget *banding_spin_invsigma = NULL;
static GtkWidget *banding_amend_note = NULL;

/* Amend mode (convergence C5b): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start.
 * Banding has no preview machinery (dialogs.c has_preview=FALSE), so there is
 * no backup to arm or clear: the display shows the installed pre-record state
 * while the form is open, and Apply/Cancel route through
 * nde_amend_preview_end. */
static gboolean banding_amend_mode = FALSE;

/* Exposed so the activate handler can cache the widget pointers before
 * the dialog is shown — previously init only ran from Apply/processing
 * handlers. */
void banding_dialog_init_statics(void);
void banding_dialog_init_statics(void) {
	if (banding_scale_amount) return;
	banding_scale_amount = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_fixbanding_amount"));
	banding_scale_invsigma = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_fixbanding_invsigma"));
	banding_protect_highlights = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_fixbanding"));
	banding_vertical = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkBandingVertical"));
	banding_seq = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkBandingSeq"));
	banding_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryBandingSeq"));
	banding_spin_invsigma = GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_fixbanding_invsigma"));
	banding_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "banding_amend_note"));
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply builds.  The
 * scale value IS sigma (apply reads sigma = invsigma scale). */
static void banding_prefill_from_amend(void) {
	struct banding_data *p = op_desc_banding.deserialize(nde_amend_preview_params(),
	                                                    nde_amend_preview_op_version());
	if (!p)
		return;
	gtk_range_set_value(banding_scale_amount, p->amount);
	gtk_range_set_value(banding_scale_invsigma, p->sigma);
	siril_toggle_set_active(GTK_WIDGET(banding_protect_highlights), p->protect_highlights);
	siril_toggle_set_active(GTK_WIDGET(banding_vertical), p->applyRotation);
	/* mirror the highlight-protect sensitivity wiring */
	gtk_widget_set_sensitive(GTK_WIDGET(banding_scale_invsigma), p->protect_highlights);
	gtk_widget_set_sensitive(banding_spin_invsigma, p->protect_highlights);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

/* Common amend-mode exit: leave amend mode (preview_end restores the true
 * pixels wholesale), then close.  No backup/preview to tear down. */
static void banding_amend_exit(gboolean apply, gchar *blob) {
	nde_amend_preview_end(apply, blob);
	banding_amend_mode = FALSE;
	siril_close_dialog("canon_fixbanding_dialog");
}

void on_canon_fixbanding_dialog_show(GtkWidget *widget, gpointer user_data) {
	banding_dialog_init_statics();
	nde_amend_note_update(banding_amend_note, banding_amend_mode,
	                      NULL);
	if (banding_amend_mode) {
		/* gfit already shows the pre-record state.  No ROI (banding is
		 * whole-image); the sequence controls are meaningless in amend mode. */
		if (gui.roi.active)
			on_clear_roi();
		gtk_widget_set_visible(GTK_WIDGET(banding_seq), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(banding_seq_entry), FALSE);
		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(banding_seq), FALSE);
		banding_prefill_from_amend();
		set_notify_block(FALSE);
	}
}

static gboolean banding_single_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();

	if (args->retval == 0) {
		gfit_modified_update_gui();
	}

	free_generic_img_args(args);
	return FALSE;
}

void on_button_ok_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	if (banding_amend_mode) {
		/* The Close button does not apply — in amend mode that means leave
		 * without committing (restore the true pixels). */
		banding_amend_exit(FALSE, NULL);
		return;
	}
	siril_close_dialog("canon_fixbanding_dialog");
}

gboolean banding_hide_on_delete(GtkWidget *widget) {
	if (banding_amend_mode) {
		nde_amend_preview_end(FALSE, NULL);
		banding_amend_mode = FALSE;
	}
	siril_close_dialog("canon_fixbanding_dialog");
	return TRUE;
}

void on_button_apply_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	if (banding_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path.  The
		 * sequence path is disabled in amend mode. */
		banding_dialog_init_statics();
		struct banding_data applied = { 0 };
		applied.amount = gtk_range_get_value(banding_scale_amount);
		applied.sigma = gtk_range_get_value(banding_scale_invsigma);
		applied.protect_highlights = siril_toggle_get_active(GTK_WIDGET(banding_protect_highlights));
		applied.applyRotation = siril_toggle_get_active(GTK_WIDGET(banding_vertical));
		gchar *blob = op_desc_banding.serialize(&applied);
		banding_amend_exit(TRUE, blob);
		g_free(blob);
		return;
	}

	if (!check_ok_if_cfa())
		return;
	double amount, invsigma;
	gboolean protect_highlights;

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	banding_dialog_init_statics();
	amount = gtk_range_get_value(banding_scale_amount);
	invsigma = gtk_range_get_value(banding_scale_invsigma);
	protect_highlights = siril_toggle_get_active(GTK_WIDGET(banding_protect_highlights));
	gboolean applyRotation = siril_toggle_get_active(GTK_WIDGET(banding_vertical));

	set_cursor_waiting(TRUE);

	if (siril_toggle_get_active(GTK_WIDGET(banding_seq)) && sequence_is_loaded()) {
		struct banding_data *seq_args = new_banding_data();
		if (!seq_args) {
			PRINT_ALLOC_ERR;
			set_cursor_waiting(FALSE);
			return;
		}

		const char *entry_text = gtk_editable_get_text(GTK_EDITABLE(banding_seq_entry));
		seq_args->seqEntry = strdup((entry_text && entry_text[0] != '\0') ? entry_text : "unband_");
		seq_args->protect_highlights = protect_highlights;
		seq_args->amount = amount;
		seq_args->sigma = invsigma;
		seq_args->applyRotation = applyRotation;
		seq_args->seq = &com.seq;
		seq_args->fit = NULL;

		siril_toggle_set_active(GTK_WIDGET(banding_seq), FALSE);
		apply_banding_to_sequence(seq_args);
	} else {
		struct banding_data *params = new_banding_data();
		if (!params) {
			PRINT_ALLOC_ERR;
			set_cursor_waiting(FALSE);
			return;
		}

		params->protect_highlights = protect_highlights;
		params->amount = amount;
		params->sigma = invsigma;
		params->applyRotation = applyRotation;
		params->seqEntry = NULL;
		params->seq = NULL;
		params->fit = NULL;

		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			PRINT_ALLOC_ERR;
			free_banding_data(params);
			free(params);
			set_cursor_waiting(FALSE);
			return;
		}

		args->fit = gfit;
		args->op = &op_desc_banding;
		args->idle_function = banding_single_idle;
		args->verbose = TRUE;
		args->user = params;
		args->max_threads = com.max_thread;
		args->for_preview = FALSE;

		if (!start_in_new_thread(generic_image_worker, args)) {
			free_banding_data(params);
			free(params);
			free(args);
			set_cursor_waiting(FALSE);
		}
	}
}

void on_checkbutton_fixbanding_toggled(GtkCheckButton *togglebutton,
		gpointer user_data) {
	banding_dialog_init_statics();
	gboolean is_active = siril_toggle_get_active(GTK_WIDGET(togglebutton));
	gtk_widget_set_sensitive(GTK_WIDGET(banding_scale_invsigma), is_active);
	gtk_widget_set_sensitive(banding_spin_invsigma, is_active);
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void banding_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	banding_amend_mode = TRUE;
	siril_open_dialog("canon_fixbanding_dialog");
}

gboolean banding_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, banding_amend_ready, NULL);
	return TRUE;
}
