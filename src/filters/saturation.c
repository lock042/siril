/*
 * Refactored saturation using generic_image_worker
 */

#include <stdlib.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"

#include "saturation.h"

static double satu_amount, background_factor;
static int satu_hue_type;
static gboolean satu_show_preview;
static int satu_update_preview();

/* Helper to map hue types to degree ranges */
void satu_set_hues_from_types(saturation_params *args, int type) {
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

/* Core Algorithm: USHORT */
static int enhance_saturation_ushort(fits *fit, saturation_params *params) {
	double bg = 0.0;
	WORD *in[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	// Operate in-place
	WORD *out[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double h_min = params->h_min / 360.0;
	double h_max = params->h_max / 360.0;

	if (params->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, fit, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (stat->median + stat->sigma) * params->background_factor;
		bg /= stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = h_min > h_max;
	double s_mult = 1.0 + params->coeff;
	double norm = fit->bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	double invnorm = 1.0 / norm;
	size_t i, n = fit->naxes[0] * fit->naxes[1];

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
			s = (s < 0.0) ? 0.0 : (s > 1.0) ? 1.0 : s;
			hsl_to_rgb(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = round_to_WORD(r * norm);
		out[GLAYER][i] = round_to_WORD(g * norm);
		out[BLAYER][i] = round_to_WORD(b * norm);
	}
	return 0;
}

/* Core Algorithm: FLOAT */
static int enhance_saturation_float(fits *fit, saturation_params *params) {
	float bg = 0.0f;
	float *in[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	// Operate in-place
	float *out[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	float h_min = (float)params->h_min / 60.0f;
	float h_max = (float)params->h_max / 60.0f;

	if (params->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, fit, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (float)((stat->median + stat->sigma) * params->background_factor);
		bg /= (float)stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = h_min > h_max;
	float s_mult = 1.f + (float)params->coeff;

	size_t i, n = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(dynamic, fit->rx * 16)
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
			s = (s < 0.f) ? 0.f : (s > 1.f) ? 1.f : s;
			hsl_to_rgb_float_sat(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = r;
		out[GLAYER][i] = g;
		out[BLAYER][i] = b;
	}
	return 0;
}

/* The Generic Processing Hook */
int saturation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	saturation_params *params = (saturation_params *)args->user;
	if (!params)
		return 1;

	if (fit->type == DATA_USHORT) {
		return enhance_saturation_ushort(fit, params);
	} else if (fit->type == DATA_FLOAT) {
		return enhance_saturation_float(fit, params);
	}
	return 1;
}

/* Idle function: PREVIEW */
static gboolean satu_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();

	if (args->retval == 0) {
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Idle function: APPLY/FINAL */
static gboolean satu_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();

	if (args->retval == 0) {
		saturation_params *params = (saturation_params *)args->user;
		// Save Undo
		undo_save_state(get_preview_gfit_backup(),
				_("Saturation enhancement (amount=%4.2lf)"), params->coeff);

		// Add to history
		char log[90];
		sprintf(log, "Color saturation %d%%, threshold %.2f",
			round_to_int(params->coeff * 100.0), params->background_factor);
		gfit->history = g_slist_append(gfit->history, strdup(log));

		populate_roi();
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Helper to launch the worker */
static int satu_process_with_worker(gboolean for_preview) {
	saturation_params *params = calloc(1, sizeof(saturation_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Destructor required by generic_img_args
	params->free = free;
	params->coeff = satu_amount;
	params->background_factor = background_factor;
	satu_set_hues_from_types(params, satu_hue_type);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	// Logic from source 1: processing.h, source 26: processing.c
	args->fit = gui.roi.active ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = saturation_image_hook;
	args->idle_function = for_preview ? satu_preview_idle : satu_apply_idle;
	args->description = _("Saturation");
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;

	if (for_preview)
		generic_image_worker(args);
	else
		start_in_new_thread(generic_image_worker, args);
	return 0;
}

static int satu_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("satu_preview")))) {
		copy_backup_to_gfit();
		return satu_process_with_worker(TRUE);
	}
	return 0;
}

void satu_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
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
			copy_backup_to_gfit();
			notify_gfit_modified();
		} else {
			clear_backup();
		}
	} else {
		// Undo is now saved in the worker idle function (satu_apply_idle)
	}
	roi_supported(FALSE);
	remove_roi_callback(satu_change_between_roi_and_image);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static void apply_satu_changes() {
	gboolean status = satu_amount != 0.0;
	satu_close(!status);
}

gboolean on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
	return FALSE;
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("satu_preview")))) {
		copy_backup_to_gfit();
	}

	// Launch worker for Final Application (FALSE = not preview)
	satu_process_with_worker(FALSE);

	// Cleanup happens in idle function, close dialog now
	siril_close_dialog("satu_dialog");
}

void on_satu_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_satu_changes();
}

/** callbacks - kept mostly the same, just update preview logic **/

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
		set_cursor_waiting(FALSE);
	}
}

void apply_satu_cancel() {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
}

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
