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
#include "core/processing.h"
#include "core/siril_log.h"
#include "filters/banding.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

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
	siril_close_dialog("canon_fixbanding_dialog");
}

gboolean banding_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("canon_fixbanding_dialog");
	return TRUE;
}

void on_button_apply_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	static GtkRange *range_amount = NULL;
	static GtkRange *range_invsigma = NULL;
	static GtkToggleButton *toggle_protect_highlights_banding = NULL,
		*vertical = NULL, *seq = NULL;
	static GtkEntry *bandingSeqEntry = NULL;
	double amount, invsigma;
	gboolean protect_highlights;

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (range_amount == NULL) {
		range_amount = GTK_RANGE(lookup_widget("scale_fixbanding_amount"));
		range_invsigma = GTK_RANGE(lookup_widget("scale_fixbanding_invsigma"));
		toggle_protect_highlights_banding = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_fixbanding"));
		vertical = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingVertical"));
		seq = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq"));
		bandingSeqEntry = GTK_ENTRY(lookup_widget("entryBandingSeq"));
	}
	amount = gtk_range_get_value(range_amount);
	invsigma = gtk_range_get_value(range_invsigma);
	protect_highlights = gtk_toggle_button_get_active(
			toggle_protect_highlights_banding);
	gboolean applyRotation = gtk_toggle_button_get_active(vertical);

	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		struct banding_data *seq_args = new_banding_data();
		if (!seq_args) {
			PRINT_ALLOC_ERR;
			set_cursor_waiting(FALSE);
			return;
		}

		const char *entry_text = gtk_entry_get_text(bandingSeqEntry);
		seq_args->seqEntry = strdup((entry_text && entry_text[0] != '\0') ? entry_text : "unband_");
		seq_args->protect_highlights = protect_highlights;
		seq_args->amount = amount;
		seq_args->sigma = invsigma;
		seq_args->applyRotation = applyRotation;
		seq_args->seq = &com.seq;
		seq_args->fit = NULL;

		gtk_toggle_button_set_active(seq, FALSE);
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
		args->mem_ratio = 2.0f;
		args->image_hook = banding_single_image_hook;
		args->idle_function = banding_single_idle;
		args->description = _("Canon Banding Reduction");
		args->verbose = TRUE;
		args->user = params;
		args->log_hook = banding_log_hook;
		args->max_threads = com.max_thread;
		args->for_preview = FALSE;
		args->for_roi = FALSE;

		if (!start_in_new_thread(generic_image_worker, args)) {
			free_banding_data(params);
			free(params);
			free(args);
			set_cursor_waiting(FALSE);
		}
	}
}

void on_checkbutton_fixbanding_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	static GtkWidget *scalebandingHighlightBox = NULL;
	static GtkWidget *spinbandingHighlightBox = NULL;
	gboolean is_active;

	if (scalebandingHighlightBox == NULL) {
		scalebandingHighlightBox = lookup_widget("scale_fixbanding_invsigma");
		spinbandingHighlightBox = lookup_widget("spin_fixbanding_invsigma");
	}

	is_active = gtk_toggle_button_get_active(togglebutton);
	gtk_widget_set_sensitive(scalebandingHighlightBox, is_active);
	gtk_widget_set_sensitive(spinbandingHighlightBox, is_active);
}
