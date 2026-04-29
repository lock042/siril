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

/* GTK callbacks and idle functions for the Cosmetic Correction dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "filters/cosmetic_correction.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

/* Idle function for the old autoDetectThreaded path */
gboolean end_autoDetect(gpointer p) {
	stop_processing_thread();
	gfit_modified_update_gui();
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
	return FALSE;
}

/* Idle function for the generic_image_worker cosmetic-correction path */
static gboolean cosme_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();

	if (args->retval == 0) {
		gfit_modified_update_gui();
	}

	free_generic_img_args(args);
	return FALSE;
}

void on_button_cosmetic_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cosmetic_dialog");
}

gboolean cosmetic_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("cosmetic_dialog");
	return TRUE;
}

void on_checkSigCosme_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkWidget *cosmeticApply = NULL;
	static GtkToggleButton *checkCosmeSigCold = NULL;
	static GtkToggleButton *checkCosmeSigHot = NULL;
	gboolean checkCold, checkHot;

	if (cosmeticApply == NULL) {
		cosmeticApply = lookup_widget("button_cosmetic_ok");
		checkCosmeSigCold = GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"));
		checkCosmeSigHot = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"));
	}
	checkCold = gtk_toggle_button_get_active(checkCosmeSigCold);
	checkHot = gtk_toggle_button_get_active(checkCosmeSigHot);
	gtk_widget_set_sensitive(cosmeticApply, checkCold || checkHot);
}

void on_button_cosmetic_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *cosmeticSeqEntry;
	GtkToggleButton *CFA, *seq;
	GtkSpinButton *sigma[2];
	GtkAdjustment *adjCosmeAmount;

	CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheckBox"));
	sigma[0] = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox"));
	sigma[1] = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox"));
	seq = GTK_TOGGLE_BUTTON(lookup_widget("checkCosmeticSeq"));
	cosmeticSeqEntry = GTK_ENTRY(lookup_widget("entryCosmeticSeq"));
	adjCosmeAmount = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjCosmeAmount"));

	struct cosmetic_data *args = calloc(1, sizeof(struct cosmetic_data));

	if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"))))
		args->sigma[0] = gtk_spin_button_get_value(sigma[0]);
	else
		args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"))))
		args->sigma[1] = gtk_spin_button_get_value(sigma[1]);
	else
		args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->amount = gtk_adjustment_get_value(adjCosmeAmount);

	args->fit = gfit;
	args->seqEntry = strdup(gtk_entry_get_text(cosmeticSeqEntry));
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("cc_");
		}
		args->seq = &com.seq;
		args->threading = SINGLE_THREADED;
		gtk_toggle_button_set_active(seq, FALSE);
		apply_cosmetic_to_sequence(args);
	} else {
		args->threading = MULTI_THREADED;

		struct generic_img_args *ga = calloc(1, sizeof(struct generic_img_args));
		if (!ga) {
			PRINT_ALLOC_ERR;
			free_cosmetic_data(args);
			return;
		}

		ga->fit = gfit;
		ga->mem_ratio = 2.f;
		ga->image_hook = cosmetic_image_hook_generic;
		ga->idle_function = cosme_idle;
		ga->description = _("Cosmetic Correction");
		ga->verbose = TRUE;
		ga->user = args;
		args->fit = gfit;
		ga->log_hook = cosmetic_log_hook;
		ga->max_threads = com.max_thread;

		if (!start_in_new_thread(generic_image_worker, ga)) {
			free_generic_img_args(ga);
			return;
		}
	}
}
