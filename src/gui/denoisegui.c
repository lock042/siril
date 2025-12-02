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

#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/dialogs.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"

// Widgets
GtkAdjustment *denoise_modulation_adjustment = NULL, *denoise_rho_adjustment = NULL, *denoise_sos_adjustment = NULL;
GtkButton *denoise_cancel = NULL, *denoise_roi_preview = NULL, *denoise_apply = NULL;
GtkDialog *denoise_dialog = NULL;
GtkFrame *denoise_artefact_control = NULL, *sos_advanced_options = NULL;
GtkRadioButton *radio_denoise_nosecondary = NULL, *radio_denoise_vst = NULL, *radio_denoise_da3d = NULL, *radio_denoise_sos = NULL;
GtkScale *slide_denoise_modulation = NULL;
GtkSpinButton *spin_sos_iters = NULL, *spin_rho = NULL, *spin_denoise_modulation = NULL;
GtkToggleButton *check_denoise_cosmetic = NULL, *check_denoise_suppress_artefacts = NULL, *denoise_preview_toggle = NULL;

// Statics init
static void denoise_dialog_init_statics() {
	if (denoise_modulation_adjustment == NULL) {
		// GtkAdjustment
		denoise_modulation_adjustment = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "denoise-modulation-adjustment"));
		denoise_rho_adjustment = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "denoise-rho-adjustment"));
		denoise_sos_adjustment = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "denoise_sos_adjustment"));
		// GtkButton
		denoise_cancel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "denoise_cancel"));
		denoise_roi_preview = GTK_BUTTON(gtk_builder_get_object(gui.builder, "denoise_roi_preview"));
		denoise_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "denoise_apply"));
		// GtkDialog
		denoise_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "denoise_dialog"));
		// GtkFrame
		denoise_artefact_control = GTK_FRAME(gtk_builder_get_object(gui.builder, "denoise_artefact_control"));
		sos_advanced_options = GTK_FRAME(gtk_builder_get_object(gui.builder, "sos_advanced_options"));
		// GtkRadioButton
		radio_denoise_nosecondary = GTK_RADIO_BUTTON(gtk_builder_get_object(gui.builder, "radio_denoise_nosecondary"));
		radio_denoise_vst = GTK_RADIO_BUTTON(gtk_builder_get_object(gui.builder, "radio_denoise_vst"));
		radio_denoise_da3d = GTK_RADIO_BUTTON(gtk_builder_get_object(gui.builder, "radio_denoise_da3d"));
		radio_denoise_sos = GTK_RADIO_BUTTON(gtk_builder_get_object(gui.builder, "radio_denoise_sos"));
		// GtkScale
		slide_denoise_modulation = GTK_SCALE(gtk_builder_get_object(gui.builder, "slide_denoise_modulation"));
		// GtkSpinButton
		spin_sos_iters = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_sos_iters"));
		spin_rho = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_rho"));
		spin_denoise_modulation = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_denoise_modulation"));
		// GtkToggleButton
		check_denoise_cosmetic = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_denoise_cosmetic"));
		check_denoise_suppress_artefacts = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_denoise_suppress_artefacts"));
		denoise_preview_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "denoise_preview"));
	}
}

/* Helper function to get current widget values */
static void get_denoise_values(struct denoise_args *params) {
	if (!params)
		return;

	params->modulation = (float)gtk_spin_button_get_value(spin_denoise_modulation);
	params->da3d = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_da3d));
	params->sos = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_sos)) ?
	              gtk_spin_button_get_value(spin_sos_iters) : 1;
	params->rho = gtk_spin_button_get_value(spin_rho);
	params->do_anscombe = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_vst));
	params->do_cosme = gtk_toggle_button_get_active(check_denoise_cosmetic);
	params->suppress_artefacts = gtk_toggle_button_get_active(check_denoise_suppress_artefacts);
}

/* Validate parameter combinations */
static gboolean validate_denoise_params(struct denoise_args *params) {
	if (params->modulation <= 0.f || params->modulation > 1.f) {
		siril_log_message(_("Error: modulation must be > 0.0 and <= 1.0.\n"));
		return FALSE;
	}

	if (params->rho == 1.f) {
		siril_log_message(_("Warning: rho = 1 will cause SOS to do nothing. Adjusting to 0.99...\n"));
		params->rho = 0.99f;
		gtk_spin_button_set_value(spin_rho, 0.99);
	}

	if (params->rho == 0.f) {
		siril_log_message(_("Warning: rho = 0 means SOS will never converge. Adjusting to 0.01...\n"));
		params->rho = 0.01f;
		gtk_spin_button_set_value(spin_rho, 0.01);
	}

	if (params->do_anscombe && (params->sos != 1 || params->da3d)) {
		siril_log_color_message(_("Error: will not carry out DA3D or SOS iterations with Anscombe transform VST selected.\n"), "red");
		return FALSE;
	}

	if (params->da3d && params->sos != 1) {
		siril_log_message(_("Will not carry out both DA3D and SOS. SOS iterations set to 1.\n"));
		params->sos = 1;
	}

	return TRUE;
}

/* Create and launch denoising processing */
static int denoise_process_with_worker(gboolean for_preview, gboolean for_roi) {
	// Allocate parameters
	struct denoise_args *params = new_denoise_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Get current values from widgets
	get_denoise_values(params);

	// Validate parameters
	if (!validate_denoise_params(params)) {
		free_denoise_args(params);
		free(params);
		return 1;
	}

	params->previewing = for_preview;

	// Log what we're doing
	if (!for_preview) {
		if (params->do_anscombe)
			siril_log_message(_("Will apply generalised Anscombe variance stabilising transform.\n"));
		if (params->da3d)
			siril_log_message(_("Will carry out final stage DA3D denoising.\n"));
	}

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_denoise_args(params);
		free(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = for_roi ? &gui.roi.fit : gfit;
	args->mem_ratio = 3.0f; // Denoising needs extra memory
	args->image_hook = denoise_image_hook;
	args->idle_function = for_preview ? denoise_preview_idle : denoise_apply_idle;
	args->description = _("NL-Bayes Denoising");
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

/* Update preview using the worker */
static int denoise_update_preview() {
	if (denoise_preview_toggle && gtk_toggle_button_get_active(denoise_preview_toggle)) {
		copy_backup_to_gfit();
		return denoise_process_with_worker(TRUE, gui.roi.active);
	}
	return 0;
}

void denoise_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// Restore original image first
	copy_backup_to_gfit();
	// If we are showing the preview, update it after the ROI change.
	restore_roi();
	notify_gfit_modified();
}

static void denoise_startup() {
	copy_gfit_to_backup();
	add_roi_callback(denoise_change_between_roi_and_image);
	roi_supported(TRUE);
}

static void denoise_close_internal(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copy_backup_to_gfit();
		notify_gfit_modified();
	} else {
		invalidate_stats_from_fit(gfit);
	}
	roi_supported(FALSE);
	remove_roi_callback(denoise_change_between_roi_and_image);
	clear_backup();
	set_cursor_waiting(FALSE);
}

void denoise_roi_callback() {
	gui.roi.operation_supports_roi = TRUE;
	if (denoise_roi_preview)
		gtk_widget_set_visible(GTK_WIDGET(denoise_roi_preview), gui.roi.active);
}

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	denoise_dialog_init_statics();

	denoise_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_denoise_modulation, 1.0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(radio_denoise_nosecondary), TRUE);
	gtk_toggle_button_set_active(check_denoise_cosmetic, TRUE);
	gtk_toggle_button_set_active(check_denoise_suppress_artefacts, FALSE);
	gtk_spin_button_set_value(spin_rho, 0.2);
	gtk_spin_button_set_value(spin_sos_iters, 1);
	if (denoise_preview_toggle)
		gtk_toggle_button_set_active(denoise_preview_toggle, FALSE);
	set_notify_block(FALSE);

	gtk_widget_set_visible(GTK_WIDGET(denoise_artefact_control), (gfit->naxes[2] == 3));
	gtk_widget_set_visible(GTK_WIDGET(sos_advanced_options), FALSE);
}

void close_denoise() {
	denoise_close_internal(FALSE);
	siril_close_dialog("denoise_dialog");
}

void on_denoise_cancel_clicked(GtkButton *button, gpointer user_data) {
	denoise_close_internal(TRUE);
	siril_close_dialog("denoise_dialog");
}

void update_sos(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_sos)))
		gtk_widget_set_visible(GTK_WIDGET(sos_advanced_options), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(sos_advanced_options), FALSE);

	// Update preview if active
	if (denoise_preview_toggle) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = denoise_update_preview;
		param->show_preview = gtk_toggle_button_get_active(denoise_preview_toggle);
		notify_update((gpointer) param);
	}
}

/*****************************************************************************
 *      I D L E   F U N C T I O N S                                         *
 ****************************************************************************/

/* Idle function for preview updates */
gboolean denoise_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Idle function for final application */
gboolean denoise_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	free_generic_img_args(args);
	// Close dialog after successful apply
	denoise_close_internal(FALSE);
	siril_close_dialog("denoise_dialog");
	return FALSE;
}

void on_denoise_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	// Check if modulation is zero
	float mod = (float)gtk_spin_button_get_value(spin_denoise_modulation);
	if (mod == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		return;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);

	// Determine if this is a preview operation
	gboolean is_preview = (button == denoise_roi_preview);

	if (is_preview) {
		// For ROI preview, just process the ROI
		if (denoise_process_with_worker(TRUE, gui.roi.active)) {
			return; // Error occurred
		}
	} else {
		// For Apply, weset the undo save state
		float mod = (float)gtk_spin_button_get_value(spin_denoise_modulation);
		undo_save_state(get_preview_gfit_backup(),
				_("NL-Bayes denoising: (modulation=%.2f)"), mod);

		// If ROI is active, restore to full image first
		if (gui.roi.active)
			restore_roi();

		// Process the full image (not ROI, not preview)
		if (denoise_process_with_worker(FALSE, FALSE)) {
			return; // Error occurred
		}
	}
}

void on_denoise_parameter_changed(GtkWidget *widget, gpointer user_data) {
	if (denoise_preview_toggle) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = denoise_update_preview;
		param->show_preview = gtk_toggle_button_get_active(denoise_preview_toggle);
		notify_update((gpointer) param);
	}
}

void on_denoise_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
		copy_backup_to_gfit();
		notify_gfit_modified();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = denoise_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
