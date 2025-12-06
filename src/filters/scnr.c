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
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
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

/*****************************************************************************
 *      S C N R      A L L O C A T O R   A N D   D E S T R U C T O R        *
 ****************************************************************************/

/* Allocator for scnr_data */
struct scnr_data *new_scnr_data() {
	struct scnr_data *args = calloc(1, sizeof(struct scnr_data));
	if (args) {
		args->destroy_fn = free_scnr_data;
	}
	return args;
}

/* Destructor for scnr_data */
void free_scnr_data(void *ptr) {
	struct scnr_data *args = (struct scnr_data *)ptr;
	if (!args)
		return;
	free(ptr);
}

gchar *scnr_log_hook(gpointer p, log_hook_detail detail) {
	struct scnr_data* args = (struct scnr_data*) p;
	return g_strdup_printf(_("SCNR: %s algorithm%s..."),
				scnr_type_to_string(args->type),
				args->preserve ? _(", preserving lightness") : "");
}

/* Subtractive Chromatic Noise Reduction - core processing function */
static int scnr_process(struct scnr_data *args, fits *fit) {
	g_assert(fit->type == DATA_USHORT || fit->type == DATA_FLOAT);
	size_t i, nbdata = fit->naxes[0] * fit->naxes[1];
	gint nb_above_1 = 0;

	gchar *msg = scnr_log_hook(args, SUMMARY);
	set_progress_bar_data(msg, PROGRESS_PULSATE);
	g_free(msg);

	double norm = get_normalized_value(fit);
	double invnorm = 1.0 / norm;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; i++) {
		double red, green, blue;
		switch (fit->type) {
			case DATA_USHORT:
				red = fit->pdata[RLAYER][i] * invnorm;
				green = fit->pdata[GLAYER][i] * invnorm;
				blue = fit->pdata[BLAYER][i] * invnorm;
				break;
			case DATA_FLOAT:
				red = (double)fit->fpdata[RLAYER][i];
				green = (double)fit->fpdata[GLAYER][i];
				blue = (double)fit->fpdata[BLAYER][i];
				break;
			default:
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

		if (fit->type == DATA_USHORT) {
			if (fit->orig_bitpix == BYTE_IMG) {
				fit->pdata[RLAYER][i] = round_to_BYTE(red * norm);
				fit->pdata[GLAYER][i] = round_to_BYTE(green * norm);
				fit->pdata[BLAYER][i] = round_to_BYTE(blue * norm);
			} else {
				fit->pdata[RLAYER][i] = round_to_WORD(red * norm);
				fit->pdata[GLAYER][i] = round_to_WORD(green * norm);
				fit->pdata[BLAYER][i] = round_to_WORD(blue * norm);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			fit->fpdata[RLAYER][i] = set_float_in_interval(red, 0.0f, 1.0f);
			fit->fpdata[GLAYER][i] = set_float_in_interval(green, 0.0f, 1.0f);
			fit->fpdata[BLAYER][i] = set_float_in_interval(blue, 0.0f, 1.0f);
		}
	}

	/* normalize in case of preserve, it can under/overshoot */
	if (args->preserve && nb_above_1)
		siril_log_message("%d pixels were truncated to a maximum value of 1\n", nb_above_1);

	if (fit == gfit && args->applying) {
		populate_roi();
	}

	return 0;
}

/* The actual SCNR processing hook for generic_image_worker */
int scnr_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct scnr_data *params = (struct scnr_data *)args->user;
	if (!params)
		return 1;
	return scnr_process(params, fit);
}

/* Idle function for preview updates */
gboolean scnr_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Idle function for final application */
gboolean scnr_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Create and launch SCNR processing */
static int scnr_process_with_worker(scnr_type type, double amount, gboolean preserve,
                                     gboolean for_preview, gboolean for_roi) {
	// Allocate parameters
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

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_scnr_data(params);
		free(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = for_roi ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.5f; // SCNR needs minimal extra memory
	args->image_hook = scnr_image_hook;
	args->idle_function = for_preview ? scnr_preview_idle : scnr_apply_idle;
	args->description = _("Subtractive Chromatic Noise Reduction");
	args->verbose = !for_preview;
	args->user = params;
	args->log_hook = scnr_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = for_roi;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}
	return 0;
}

/* Update preview */
static int scnr_update_preview() {
	GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
	if (gtk_toggle_button_get_active(preview_button)) {
		int type = gtk_combo_box_get_active(
				GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_scnr")));
		GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(gui.builder, "preserve_light"));
		gboolean preserve = gtk_toggle_button_get_active(light_button);
		double amount = gtk_range_get_value(
				GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_scnr")));

		copy_backup_to_gfit();
		return scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	}
	return 0;
}

void scnr_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
	param->show_preview = gtk_toggle_button_get_active(preview_button);
	notify_update((gpointer) param);
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	// Notify the overlay that this dialog supports ROI processing
	roi_supported(TRUE);
	gtk_widget_set_visible(lookup_widget("SCNR_roi_preview"), gui.roi.active);

	copy_gfit_to_backup();
	add_roi_callback(scnr_change_between_roi_and_image);

	if (gui.roi.active) {
		// Call this directly on startup to set the ROI preview
		scnr_change_between_roi_and_image();
	}

	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(gui.builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	// Notify the overlay that we are leaving a dialog that supports ROI
	roi_supported(FALSE);

	// Revert to backup if preview was active
	GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
	if (gtk_toggle_button_get_active(preview_button)) {
		copy_backup_to_gfit();
		notify_gfit_modified();
	}

	clear_backup();
	remove_roi_callback(scnr_change_between_roi_and_image);
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

	// Check if this is a preview click or apply click
	gboolean is_preview = ((GtkWidget*) button == lookup_widget("SCNR_roi_preview"));

	if (is_preview) {
		// For preview with ROI
		scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	} else {
		// For final apply - always process full image
		// If preview was on, restore backup first
		GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
		if (gtk_toggle_button_get_active(preview_button)) {
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

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));
	GtkSpinButton *spinButton = GTK_SPIN_BUTTON(lookup_widget("spin_scnr"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(spinButton), type > 1);

	// Update preview if active
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
	param->show_preview = gtk_toggle_button_get_active(preview_button);
	notify_update((gpointer) param);
}

void on_SCNR_parameter_changed(GtkWidget *widget, gpointer user_data) {
	// Update preview if active
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	GtkToggleButton *preview_button = GTK_TOGGLE_BUTTON(lookup_widget("SCNR_roi_preview"));
	param->show_preview = gtk_toggle_button_get_active(preview_button);
	notify_update((gpointer) param);
}

void on_SCNR_roi_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		waiting_for_thread();
		siril_preview_hide();
		copy_backup_to_gfit();
		notify_gfit_modified();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = scnr_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
