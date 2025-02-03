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
#include "core/proto.h"
#include "core/undo.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/siril_preview.h" // For copy_gfit_to_backup() etc
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"

#include "scnr.h"

const char *scnr_type_to_string(scnr_type t) {
	switch (t) {
		default:
		case SCNR_AVERAGE_NEUTRAL:
			return _("average neutral");
		case SCNR_MAXIMUM_NEUTRAL:
			return _("maximum neutral");
		case SCNR_MAXIMUM_MASK:
			return _("maximum mask");
		case SCNR_ADDITIVE_MASK:
			return _("additive mask");
	}
}

/* Subtractive Chromatic Noise Reduction */
gpointer scnr(gpointer p) {
	lock_roi_mutex(); // this prevents changes to the ROI occurring while the thread
	// is running as this could rip out data structures from under us

	// Ensure we are starting with a fresh copy of gfit with no previous
	// ROI preview data
	copy_backup_to_gfit();
	struct scnr_data *args = (struct scnr_data *) p;
	g_assert(args->fit->type == DATA_USHORT || args->fit->type == DATA_FLOAT);
	size_t i, nbdata = args->fit->naxes[0] * args->fit->naxes[1];
	gint nb_above_1 = 0;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	double norm = get_normalized_value(args->fit);
	double invnorm = 1.0 / norm;
	// If we aren't being run from a command or ROI previewing, set undo state
	if (!com.script && !args->previewing)
		undo_save_state(&gfit, _("SCNR (type=%s, amount=%0.2lf, preserve=%s)"),
				scnr_type_to_string(args->type), args->amount, args->preserve ? "true" : "false");

	siril_log_color_message(_("SCNR: processing with %s algorithm%s...\n"),
			"green", scnr_type_to_string(args->type),
			args->preserve ? _(", preserving lightness") : "");

	int error = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; i++) {
		double red, green, blue;
		switch (args->fit->type) {
			case DATA_USHORT:
				red = args->fit->pdata[RLAYER][i] * invnorm;
				green = args->fit->pdata[GLAYER][i] * invnorm;
				blue = args->fit->pdata[BLAYER][i] * invnorm;
				break;
			case DATA_FLOAT:
				red = (double)args->fit->fpdata[RLAYER][i];
				green = (double)args->fit->fpdata[GLAYER][i];
				blue = (double)args->fit->fpdata[BLAYER][i];
				break;
			default: // Default needs to be included in the switch to avoid warning about omitting DATA_UNSUPPORTED
				break;
		}

		double x, y, z, L, a, b, m;
		if (args->preserve) {
			linrgb_to_xyz(red, green, blue, &x, &y, &z, TRUE);
			xyz_to_LAB(x, y, z, &L, &a, &b);
		}

		switch (args->type) {
			case SCNR_AVERAGE_NEUTRAL:
				m = 0.5 * (red + blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_NEUTRAL:
				m = max(red, blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_MASK:
				m = max(red, blue);
				green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
				break;
			case SCNR_ADDITIVE_MASK:
				m = min(1.0, red + blue);
				green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
		}

		if (args->preserve) {
			double tmp;
			linrgb_to_xyz(red, green, blue, &x, &y, &z, TRUE);
			xyz_to_LAB(x, y, z, &tmp, &a, &b);
			LAB_to_xyz(L, a, b, &x, &y, &z);
			xyz_to_linrgb(x, y, z, &red, &green, &blue, TRUE);
			if (red > 1.000001 || green > 1.000001 || blue > 1.000001)
				g_atomic_int_inc(&nb_above_1);
		}

		if (args->fit->type == DATA_USHORT) {
			if (args->fit->orig_bitpix == BYTE_IMG) {
				args->fit->pdata[RLAYER][i] = round_to_BYTE(red * norm);
				args->fit->pdata[GLAYER][i] = round_to_BYTE(green * norm);
				args->fit->pdata[BLAYER][i] = round_to_BYTE(blue * norm);
			} else {
				args->fit->pdata[RLAYER][i] = round_to_WORD(red * norm);
				args->fit->pdata[GLAYER][i] = round_to_WORD(green * norm);
				args->fit->pdata[BLAYER][i] = round_to_WORD(blue * norm);
			}
		}
		else if (args->fit->type == DATA_FLOAT) {
			args->fit->fpdata[RLAYER][i] = set_float_in_interval(red, 0.0f, 1.0f);
			args->fit->fpdata[GLAYER][i] = set_float_in_interval(green, 0.0f, 1.0f);
			args->fit->fpdata[BLAYER][i] = set_float_in_interval(blue, 0.0f, 1.0f);
		}
	}

	/* normalize in case of preserve, it can under/overshoot */
	if (args->preserve && nb_above_1)
		siril_log_message("%d pixels were truncated to a maximum value of 1\n", nb_above_1);

	// If we are not previewing, we now update the backup *and the ROI*
	// This must be done before notify_gfit_modified() is called or the ROI will
	// be drawn with the old data
	if (!args->previewing) {
		copy_gfit_to_backup();
		populate_roi();
	}

	free(args);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	unlock_roi_mutex(); // We're done, so unlock the ROI mutex to allow ROI changes to happen again
	notify_gfit_modified();
	return GINT_TO_POINTER(error);
}

void scnr_roi_callback() {
	gui.roi.operation_supports_roi = TRUE;
	// Sets the visbility of the ROI preview button
	gtk_widget_set_visible(lookup_widget("SCNR_roi_preview"), gui.roi.active);
	// The next 2 lines ensure gfit looks fresh again and there isn't the
	// remains of an old preview in the wrong place
	copy_backup_to_gfit();
	notify_gfit_modified();
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	// Notify the overlay that this dialog supports ROI processing
	roi_supported(TRUE);
	gtk_widget_set_visible(lookup_widget("SCNR_roi_preview"), gui.roi.active);
	if (gui.roi.active) {
		copy_gfit_to_backup();
		// Call this directly on startup to set the ROI preview
		scnr_roi_callback();
	}
	add_roi_callback(scnr_roi_callback);
	// Set up the backup (this is used to go back to a fresh copy of the
	// image after abandoning a ROI preview)
	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(gui.builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	// Notify the overlay that we are leaving a dialog that supports ROI
	roi_supported(FALSE);
	// Get rid of the gfit backup
	if (is_preview_active())
		siril_preview_hide();
	// Remove the callback
	remove_roi_callback(scnr_roi_callback);
	siril_close_dialog("SCNR_dialog");
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_scnr")));
	GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(gui.builder, "preserve_light"));
	gboolean preserve = gtk_toggle_button_get_active(light_button);
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_scnr")));

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	// Tell the threaded function if this is a preview or for real
	args->previewing = ((GtkWidget*) button == lookup_widget("SCNR_roi_preview"));
	// undo_save_state(...) // We don't do this here, it has to be done
	// in the thread after copy_gfit_to_backup

	// The fit to be operated on depends on if we are previewing with a ROI
	// active or carrying out the final operation on gfit
	args->fit = (args->previewing && gui.roi.active) ? &gui.roi.fit : &gfit;
	args->type = type;
	args->amount = amount;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(scnr, args))
		free(args);
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));
	GtkSpinButton *spinButton = GTK_SPIN_BUTTON(lookup_widget("spin_scnr"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(spinButton), type > 1);
}
