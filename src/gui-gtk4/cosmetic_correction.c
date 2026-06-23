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
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

static GtkWidget *cosme_apply_btn = NULL;
static GtkCheckButton *cosme_check_cold = NULL, *cosme_check_hot = NULL;
static GtkCheckButton *cosme_cfa = NULL, *cosme_seq = NULL;
static GtkSpinButton *cosme_sigma_cold = NULL, *cosme_sigma_hot = NULL;
static GtkEntry *cosme_seq_entry = NULL;

/* Exposed so the activate handler can cache the widget pointers before
 * the dialog is shown — previously init only ran from the toggle and
 * Apply handlers, so widget pointers were NULL until the user touched
 * the dialog. */
void cosmetic_dialog_init_statics(void);
void cosmetic_dialog_init_statics(void) {
	if (cosme_apply_btn) return;
	cosme_apply_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_cosmetic_ok"));
	cosme_check_cold = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkSigColdBox"));
	cosme_check_hot = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkSigHotBox"));
	cosme_cfa = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "cosmCFACheckBox"));
	cosme_sigma_cold = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinSigCosmeColdBox"));
	cosme_sigma_hot = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinSigCosmeHotBox"));
	cosme_seq = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkCosmeticSeq"));
	cosme_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryCosmeticSeq"));
}

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

void on_checkSigCosme_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	cosmetic_dialog_init_statics();
	gboolean checkCold = siril_toggle_get_active(GTK_WIDGET(cosme_check_cold));
	gboolean checkHot = siril_toggle_get_active(GTK_WIDGET(cosme_check_hot));
	gtk_widget_set_sensitive(cosme_apply_btn, checkCold || checkHot);
}

void on_button_cosmetic_ok_clicked(GtkButton *button, gpointer user_data) {
	cosmetic_dialog_init_statics();
	GtkAdjustment *adjCosmeAmount = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjCosmeAmount"));

	struct cosmetic_data *args = calloc(1, sizeof(struct cosmetic_data));

	if (siril_toggle_get_active(GTK_WIDGET(cosme_check_cold)))
		args->sigma[0] = gtk_spin_button_get_value(cosme_sigma_cold);
	else
		args->sigma[0] = -1.0;

	if (siril_toggle_get_active(GTK_WIDGET(cosme_check_hot)))
		args->sigma[1] = gtk_spin_button_get_value(cosme_sigma_hot);
	else
		args->sigma[1] = -1.0;

	args->is_cfa = siril_toggle_get_active(GTK_WIDGET(cosme_cfa));
	args->amount = gtk_adjustment_get_value(adjCosmeAmount);

	args->fit = gfit;
	args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(cosme_seq_entry)));
	set_cursor_waiting(TRUE);

	if (siril_toggle_get_active(GTK_WIDGET(cosme_seq)) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("cc_");
		}
		args->seq = &com.seq;
		args->threading = SINGLE_THREADED;
		siril_toggle_set_active(GTK_WIDGET(cosme_seq), FALSE);
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
