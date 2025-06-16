/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/undo.h"
#include "algos/background_extraction.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"

static poly_order get_poly_order() {
	GtkComboBox *combo_box_poly_order = GTK_COMBO_BOX(lookup_widget("box_background_order"));
	return gtk_combo_box_get_active(combo_box_poly_order);
}

static background_correction get_correction_type() {
	GtkComboBox *combo_box_correction = GTK_COMBO_BOX(lookup_widget("box_background_correction"));
	return gtk_combo_box_get_active(combo_box_correction);
}

static int get_nb_samples_per_line() {
	GtkSpinButton *nb_samples = GTK_SPIN_BUTTON(lookup_widget("spin_background_nb_samples"));
	return gtk_spin_button_get_value_as_int(nb_samples);
}

static double get_tolerance_value() {
	GtkRange *tol = GTK_RANGE(lookup_widget("scale_background_tolerance"));
	return gtk_range_get_value(tol);
}

static background_interpolation get_interpolation_method() {
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("background_extraction_combo"));
	return gtk_combo_box_get_active(combo);
}

static double get_smoothing_parameter() {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_background_smoothing"));
	return gtk_spin_button_get_value(spin);
}

static gboolean is_dither_checked() {
	return (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bkg_dither_button"))));
}

/* management of the graphical state and backup image */
static fits background_backup;
static gboolean background_computed = FALSE;

static void background_startup() {
	copy_gfit_to_backup();
}

static void copy_gfit_to_bkg_backup() {
	if (!background_computed) return;
	if (copyfits(&gfit, &background_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_debug_print("Image copy error in previews\n");
		return;
	}
}

static int copy_bkg_backup_to_gfit() {
	if (!background_computed) return 0;
	int retval = 0;
	if (!gfit.data && !gfit.fdata)
		retval = 1;
	else if (copyfits(&background_backup, &gfit, CP_COPYA, -1)) {
		siril_debug_print("Image copy error in previews\n");
		retval = 1;
	}
	return retval;
}

gboolean end_background(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	stop_processing_thread();
	if (args) {
		background_computed = TRUE;
		notify_gfit_modified();
		gtk_widget_set_sensitive(lookup_widget("background_ok_button"), TRUE);
		gtk_widget_set_sensitive(lookup_widget("bkg_show_original"), TRUE);
		free(args);
	}
	set_cursor_waiting(FALSE);
	return FALSE;
}

/************* CALLBACKS *************/

void on_background_generate_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	int nb_of_samples;
	double tolerance;
	GtkToggleButton* keep_all_button = (GtkToggleButton*) lookup_widget("subsky_keep_samples");
	gboolean keep_all = gtk_toggle_button_get_active(keep_all_button);
	nb_of_samples = get_nb_samples_per_line();
	tolerance = keep_all ? -1. : get_tolerance_value();

	if (generate_background_samples(nb_of_samples, tolerance))
		control_window_switch_to_tab(OUTPUT_LOGS);
	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
}

void on_background_clear_all_clicked(GtkButton *button, gpointer user_data) {
	sample_mutex_lock();
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	sample_mutex_unlock();

	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
}

void on_bkg_compute_bkg_clicked(GtkButton *button, gpointer user_data) {
	if (com.grad_samples == NULL) {
		return;
	}
	set_cursor_waiting(TRUE);
	copy_backup_to_gfit();

	background_correction correction = get_correction_type();
	poly_order degree = get_poly_order();
	gboolean use_dither = is_dither_checked();
	double smoothing = get_smoothing_parameter();
	background_interpolation interpolation_method = get_interpolation_method();

	struct background_data *args = calloc(1, sizeof(struct background_data));
	args->threads = com.max_thread;
	args->from_ui = TRUE;
	args->correction = correction;
	args->interpolation_method = interpolation_method;
	args->degree = (poly_order) degree;
	args->smoothing = smoothing;
	args->dither = use_dither;
	args->fit = &gfit;

	// Check if the image has a Bayer CFA pattern
	gboolean is_cfa = gfit.naxes[2] == 1 && (!strncmp(gfit.keywords.bayer_pattern, "RGGB", 4) ||
					  !strncmp(gfit.keywords.bayer_pattern, "BGGR", 4) ||
					  !strncmp(gfit.keywords.bayer_pattern, "GBRG", 4) ||
					  !strncmp(gfit.keywords.bayer_pattern, "GRBG", 4));
	if (!start_in_new_thread(is_cfa ? remove_gradient_from_cfa_image :
						remove_gradient_from_image, args)) {
		free(args->seqEntry);
		free(args);
	}
}

void on_background_ok_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkBkgSeq"));
	if (gtk_toggle_button_get_active(seq_button) && sequence_is_loaded()) {
		struct background_data *args = calloc(1, sizeof(struct background_data));
		args->nb_of_samples = get_nb_samples_per_line();
		args->tolerance = get_tolerance_value();
		args->correction = get_correction_type();
		args->degree = get_poly_order();
		args->smoothing = get_smoothing_parameter();
		args->dither = is_dither_checked();
		args->interpolation_method = get_interpolation_method();

		if (args->interpolation_method == BACKGROUND_INTER_POLY && args->degree > BACKGROUND_POLY_1) {
			int confirm = siril_confirm_dialog(_("Polynomial order seems too high."),
					_("You are about to process a sequence of preprocessed files with "
						"a polynomial degree greater than 1. This is unlikely because such "
						"gradients are often linear and a correction with a polynomial "
						"function of degree 1 is probably enough."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		} else if (args->interpolation_method != BACKGROUND_INTER_POLY) {
			int confirm = siril_confirm_dialog(_("Using wrong interpolation method"),
					_("You are about to process a sequence of preprocessed files with an RBF algorithm. "
						"This algorithm may not be very well suited for automated processing "
						"and we advise you to use the polynomial algorithm with a "
						"degree order of 1."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		}

		set_cursor_waiting(TRUE);

		args->seqEntry = strdup( gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryBkgSeq"))));
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("bkg_");
		}
		args->seq = &com.seq;
		/* now we uncheck the button */
		gtk_toggle_button_set_active(seq_button, FALSE);
		apply_background_extraction_to_sequence(args);
	} else {
		if (background_computed) {
			background_correction correction = get_correction_type();
			undo_save_state(get_preview_gfit_backup(), _("Background extraction (Correction: %s)"),
					correction == BACKGROUND_CORRECTION_DIVIDE ? "Division" : "Subtraction");
			background_computed = FALSE;
			clear_backup();
			siril_close_dialog("background_extraction_dialog");
		} else {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("No Background model computed"),
					_("You must first compute the background model."));
		}
	}
}

void apply_background_cancel() {
	siril_close_dialog("background_extraction_dialog");
}

void on_background_close_button_clicked(GtkButton *button, gpointer user_data) {
	apply_background_cancel();
}

gboolean bge_hide_on_delete(GtkWidget *widget) {
	apply_background_cancel();
	return TRUE;
}

void on_background_extraction_dialog_hide(GtkWidget *widget, gpointer user_data) {
	sample_mutex_lock();
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	sample_mutex_unlock();
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(REDRAW_OVERLAY);

	if (background_computed) {
		siril_preview_hide();
		background_computed = FALSE;
	} else {
		clear_backup();
	}
	gtk_widget_set_sensitive(lookup_widget("background_ok_button"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("bkg_show_original"), FALSE);
}

void on_background_extraction_dialog_show(GtkWidget *widget, gpointer user_data) {
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	background_startup();
}

void on_background_extraction_combo_changed(GtkComboBox *combo, gpointer user_data) {
	GtkNotebook *notebook = GTK_NOTEBOOK(lookup_widget("bkg_notebook_inter"));
	gtk_notebook_set_current_page(notebook, gtk_combo_box_get_active(combo));
}

static gboolean pressed = FALSE;

void on_bkg_show_original_button_press_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	GtkStateFlags new_state;
	new_state = gtk_widget_get_state_flags(widget)
		& ~(GTK_STATE_FLAG_PRELIGHT | GTK_STATE_FLAG_ACTIVE);

	new_state |= GTK_STATE_FLAG_PRELIGHT;
	new_state |= GTK_STATE_FLAG_ACTIVE;

	gtk_widget_set_state_flags(widget, new_state, TRUE);
	pressed = TRUE;

	copy_gfit_to_bkg_backup();
	copy_backup_to_gfit();
	notify_gfit_modified();
}

void on_bkg_show_original_button_release_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	GtkStateFlags new_state;
	new_state = gtk_widget_get_state_flags(widget)
		& ~(GTK_STATE_FLAG_PRELIGHT | GTK_STATE_FLAG_ACTIVE);

	new_state |= GTK_STATE_FLAG_PRELIGHT;

	gtk_widget_set_state_flags(widget, new_state, TRUE);
	pressed = FALSE;

	copy_bkg_backup_to_gfit();
	clearfits(&background_backup);
	notify_gfit_modified();
}

gboolean on_bkg_show_original_enter_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	GtkStateFlags new_state;
	new_state = gtk_widget_get_state_flags(widget)
		& ~(GTK_STATE_FLAG_PRELIGHT | GTK_STATE_FLAG_ACTIVE);

	new_state |= GTK_STATE_FLAG_PRELIGHT;
	if (pressed) {
		new_state |= GTK_STATE_FLAG_ACTIVE;

	}
	gtk_widget_set_state_flags(widget, new_state, TRUE);
	return TRUE;
}

void on_checkBkgSeq_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *ok = lookup_widget("background_ok_button");
	if (gtk_toggle_button_get_active(button)) {
		gtk_widget_set_sensitive(ok, TRUE);
	} else {
		gtk_widget_set_sensitive(ok, (com.grad_samples != NULL) && background_computed);
	}
}
