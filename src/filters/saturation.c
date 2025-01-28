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

#include <stdlib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"

#include "saturation.h"

static double satu_amount, background_factor;
static int satu_hue_type;
static gboolean satu_show_preview;
static int satu_update_preview();

void satu_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

static void satu_startup() {
	roi_supported(TRUE);
	add_roi_callback(satu_change_between_roi_and_image);
	copy_gfit_to_backup();
	satu_amount = 0.0;
	satu_hue_type = 6;
}

static void satu_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		if (satu_amount != 0.0) {
			siril_preview_hide();
		} else {
			clear_backup();
		}
	} else {
		undo_save_state(get_preview_gfit_backup(),
				_("Saturation enhancement (amount=%4.2lf)"), satu_amount);
	}
	roi_supported(FALSE);
	remove_roi_callback(satu_change_between_roi_and_image);
	clear_backup();
}

static void apply_satu_changes() {
	gboolean status = satu_amount != 0.0;
	satu_close(!status);
}

void satu_set_hues_from_types(struct enhance_saturation_data *args, int type) {
	switch (type) {
		case 0:		// Pink-Red to Red-Orange
			args->h_min = 346.0;
			args->h_max = 20.0;
			break;
		case 1:		// Orange-Brown to Yellow
			args->h_min = 21.0;
			args->h_max = 60.0;
			break;
		case 2:		// Yellow-Green to Green-Cyan
			args->h_min = 61.0;
			args->h_max = 200.0;
			break;
		case 3:		// Cyan
			args->h_min = 170.0;
			args->h_max = 200.0;
			break;
		case 4:		// Cyan-Blue to Blue-Magenta
			args->h_min = 201.0;
			args->h_max = 280.0;
			break;
		case 5:		// Magenta to Pink
			args->h_min = 281.0;
			args->h_max = 345.0;
			break;
		default:
		case 6:		// Global
			args->h_min = 0.0;
			args->h_max = 360.0;
	}
}

static int satu_process_all() {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	set_cursor_waiting(TRUE);
	if (satu_show_preview)
		copy_backup_to_gfit();
	else if (gui.roi.active)
		restore_roi();

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	satu_set_hues_from_types(args, satu_hue_type);

	args->input = &gfit;
	args->output = &gfit;
	args->coeff = satu_amount;
	args->background_factor = background_factor;
	args->for_preview = TRUE;
	args->for_final = TRUE;

	if (!start_in_new_thread(enhance_saturation, args))
		free(args);

	return 0;
}

static int satu_update_preview() {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	set_cursor_waiting(TRUE);
	if (satu_show_preview)
		copy_backup_to_gfit();
	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	satu_set_hues_from_types(args, satu_hue_type);

	args->input = fit;
	args->output = fit;
	args->coeff = satu_amount;
	args->background_factor = background_factor;
	args->for_preview = TRUE;
	args->for_final = FALSE;

	if (!start_in_new_thread(enhance_saturation, args))
		free(args);

	return 0;
}

void on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	if (satu_show_preview == FALSE || gui.roi.active) {
		satu_process_all();
	}

	apply_satu_changes();
	siril_close_dialog("satu_dialog");
}

void on_satu_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_satu_changes();
}

static int enhance_saturation_ushort(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	double bg = 0.0;

	WORD *in[3] = { args->input->pdata[RLAYER], args->input->pdata[GLAYER],
		args->input->pdata[BLAYER] };
	WORD *out[3] = { args->output->pdata[RLAYER], args->output->pdata[GLAYER],
		args->output->pdata[BLAYER] };

	args->h_min /= 360.0;
	args->h_max /= 360.0;
	if (args->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, args->input, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (stat->median + stat->sigma) * args->background_factor;
		bg /= stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = args->h_min > args->h_max;
	double s_mult = 1.0 + args->coeff;
	float h_min = args->h_min;
	float h_max = args->h_max;
	double norm = args->input->bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	double invnorm = 1.0 / norm;
	size_t i, n = args->input->naxes[0] * args->input->naxes[1];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < n; i++) {
		double h, s, l;
		double r = in[RLAYER][i] * invnorm;
		double g = in[GLAYER][i] * invnorm;
		double b = in[BLAYER][i] * invnorm;
		rgb_to_hsl(r, g, b, &h, &s, &l);
		if (l > bg) {
			if (loop_range) {
				if (h >= h_min || h <= h_max)
					s *= s_mult;
			} else {
				if (h >= h_min && h <= h_max)
					s *= s_mult;
			}
			if (s < 0.0) s = 0.0;
			else if (s > 1.0) s = 1.0;

			hsl_to_rgb(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = round_to_WORD(r * norm);
		out[GLAYER][i] = round_to_WORD(g * norm);
		out[BLAYER][i] = round_to_WORD(b * norm);
	}
	return 0;
}

static int enhance_saturation_float(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	float bg = 0.0f;

	float *in[3] = { args->input->fpdata[RLAYER], args->input->fpdata[GLAYER],
		args->input->fpdata[BLAYER] };
	float *out[3] = { args->output->fpdata[RLAYER], args->output->fpdata[GLAYER],
		args->output->fpdata[BLAYER] };

	args->h_min /= 60.0;
	args->h_max /= 60.0;
	if (args->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, args->input, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (stat->median + stat->sigma) * args->background_factor;
		bg /= stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = args->h_min > args->h_max;
	float s_mult = 1.f + args->coeff;
	float h_min = args->h_min;
	float h_max = args->h_max;

	size_t i, n = args->input->naxes[0] * args->input->naxes[1];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(dynamic, args->input->rx * 16)
#endif
	for (i = 0; i < n; i++) {
		float h, s, l;
		float r = in[RLAYER][i];
		float g = in[GLAYER][i];
		float b = in[BLAYER][i];
		rgb_to_hsl_float_sat(r, g, b, bg, &h, &s, &l);
		if (l > bg) {
			if (loop_range) {
				if (h >= h_min || h <= h_max)
					s *= s_mult;
			} else {
				if (h >= h_min && h <= h_max)
					s *= s_mult;
			}
			if (s < 0.f) s = 0.f;
			else if (s > 1.f) s = 1.f;

			hsl_to_rgb_float_sat(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = r;
		out[GLAYER][i] = g;
		out[BLAYER][i] = b;
	}
	return 0;
}

gpointer enhance_saturation(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;

	int retval = -1;
	if (args->input->type == DATA_USHORT) {
		retval = enhance_saturation_ushort(args);
	} else if (args->input->type == DATA_FLOAT) {
		retval = enhance_saturation_float(args);
	}

	if (!args->for_preview) {
		char log[90];
		sprintf(log, "Color saturation %d%%, threshold %.2f",
			round_to_int(args->coeff * 100.0), args->background_factor);
		args->output->history = g_slist_append(args->output->history, strdup(log));

	}
	if (args->for_final)
		populate_roi();
	notify_gfit_modified();

	free(args);
	return GINT_TO_POINTER(retval);
}

/** callbacks **/

void on_satu_dialog_show(GtkWidget *widget, gpointer user_data) {
	satu_startup();
	satu_amount = 0.0;
	satu_hue_type = 6;
	background_factor = 1.0;

	set_notify_block(TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_saturation")), satu_hue_type);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), satu_amount);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu_bkg")), background_factor);
	set_notify_block(FALSE);

	satu_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("satu_preview")));
}

void on_combo_saturation_changed(GtkComboBox* box, gpointer user_data) {
	satu_hue_type = gtk_combo_box_get_active(box);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_undo_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	double prev_satu = satu_amount;

	set_notify_block(TRUE);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), 0);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu_bkg")), 1);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("satu_preview")), TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_saturation")), 6);
	set_notify_block(FALSE);

	// Update preview only if required
	if (prev_satu != 0.0) {
		copy_backup_to_gfit();
		notify_gfit_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
		set_cursor_waiting(FALSE);
	}
}

void apply_satu_cancel() {
	satu_close(TRUE);
}

/*** adjusters **/
void on_spin_satu_value_changed(GtkSpinButton *button, gpointer user_data) {
	satu_amount = gtk_spin_button_get_value(button);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = 	satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_spin_satu_bkg_value_changed(GtkSpinButton *button, gpointer user_data) {
	background_factor = gtk_spin_button_get_value(button);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = 	satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	satu_show_preview = gtk_toggle_button_get_active(button);
	if (!satu_show_preview) {
		/* if user click very fast */
		waiting_for_thread();
		copy_backup_to_gfit();
		redraw(REMAP_ALL);
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = satu_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
