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
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/dialogs.h"
#include "gui/callbacks.h"
#include "core/processing.h"
#include "core/command.h"
#include "io/single_image.h"

// Widgets
GtkAdjustment *denoise_modulation_adjustment = NULL, *denoise_rho_adjustment = NULL, *denoise_sos_adjustment = NULL;
GtkButton *denoise_cancel = NULL, *denoise_roi_preview = NULL, *denoise_apply = NULL;
GtkDialog *denoise_dialog = NULL;
GtkFrame *denoise_artefact_control = NULL, *sos_advanced_options = NULL;
GtkRadioButton *radio_denoise_nosecondary = NULL, *radio_denoise_vst = NULL, *radio_denoise_da3d = NULL, *radio_denoise_sos = NULL;
GtkScale *slide_denoise_modulation = NULL;
GtkSpinButton *spin_sos_iters = NULL, *spin_rho = NULL, *spin_denoise_modulation = NULL;
GtkToggleButton *check_denoise_cosmetic = NULL, *check_denoise_suppress_artefacts = NULL;

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
	}
}

void denoise_roi_callback() {
	gui.roi.operation_supports_roi = TRUE;
	gtk_widget_set_visible(lookup_widget("denoise_roi_preview"), gui.roi.active);
	copy_backup_to_gfit();
	notify_gfit_modified();
}

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	denoise_dialog_init_statics();
	copy_gfit_to_backup();
	roi_supported(TRUE);
	denoise_roi_callback();
	add_roi_callback(denoise_roi_callback);
	gtk_spin_button_set_value(spin_denoise_modulation, 1.0);
	gtk_widget_set_visible(GTK_WIDGET(denoise_artefact_control), (gfit->naxes[2] == 3));
}

void close_denoise() {
	roi_supported(FALSE);
	siril_preview_hide();
	remove_roi_callback(denoise_roi_callback);
	siril_close_dialog("denoise_dialog");
}

void on_denoise_cancel_clicked(GtkButton *button, gpointer user_data) {
	close_denoise();
}

void update_sos(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_sos)))
		gtk_widget_set_visible(GTK_WIDGET(sos_advanced_options), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(sos_advanced_options), FALSE);
}

void on_denoise_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (gui.roi.active)
		restore_roi();
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->modulation = (float)gtk_spin_button_get_value(spin_denoise_modulation);
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		free(args);
		return;
	}
	control_window_switch_to_tab(OUTPUT_LOGS);
	args->previewing = (button == denoise_roi_preview);
	siril_debug_print("Previewing: %d\n", args->previewing);
	args->da3d = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_da3d));
//	args->sos = 1;
	args->sos = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_sos));
	args->rho = gtk_spin_button_get_value(spin_rho);
	if (args->rho == 1.f) {
		siril_log_message("Warning: rho = 1 will cause SOS to do nothing. Adjusting to 0.99...\n");
		args->rho = 0.99f;
		gtk_spin_button_set_value(spin_rho, 0.99);
	}
	if (args->rho == 0.f) {
		siril_log_message("Warning: rho = 0 means SOS will never converge. Adjusting to 0.01...\n");
		args->rho = 0.01f;
		gtk_spin_button_set_value(spin_rho, 0.01);
	}
	args->do_anscombe = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio_denoise_vst));
	args->do_cosme = gtk_toggle_button_get_active(check_denoise_cosmetic);
	args->suppress_artefacts = gtk_toggle_button_get_active(check_denoise_suppress_artefacts);
	if (args->sos)
		args->sos = gtk_spin_button_get_value(spin_sos_iters);
	args->fit = (gui.roi.active && args->previewing) ? &gui.roi.fit : gfit;
	if (!start_in_new_thread(run_nlbayes_on_fit, args)) {
		free(args);
	}
}
