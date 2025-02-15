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
#include "algos/statistics.h"
#include "filters/median.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/dialogs.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/command.h"
#include "io/single_image.h"

// Callbacks
float denoise_modulation;
int da3d = 0;
int sos = 0;
int sos_iters = 3;
float sos_rho = 0.2f;
gboolean do_anscombe = FALSE;
gboolean do_cosme = TRUE;
gboolean suppress_artefacts = FALSE;

void denoise_roi_callback() {
	gui.roi.operation_supports_roi = TRUE;
	gtk_widget_set_visible(lookup_widget("denoise_roi_preview"), gui.roi.active);
	copy_backup_to_gfit();
	notify_gfit_modified();
}

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	copy_gfit_to_backup();
	roi_supported(TRUE);
	denoise_roi_callback();
	add_roi_callback(denoise_roi_callback);
	da3d = 0;
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = 1.f;
	gtk_spin_button_set_value(spin_denoise_modulation, denoise_modulation);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("denoise_artefact_control")), (gfit.naxes[2] == 3));
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

void on_spin_sos_iters_value_changed(GtkSpinButton *button, gpointer user_data) {
	sos_iters = (float)gtk_spin_button_get_value(button);
}
void on_spin_rho_value_changed(GtkSpinButton *button, gpointer user_data) {
	sos_rho = (float)gtk_spin_button_get_value(button);
	if (sos_rho == 1.f) {
		siril_log_message("Warning: rho = 1 will cause SOS to do nothing. Adjusting to 0.99...\n");
		sos_rho = 0.99f;
		gtk_spin_button_set_value(button, 0.99);
	}
	if (sos_rho == 0.f) {
		siril_log_message("Warning: rho = 0 means SOS will never converge. Adjusting to 0.01...\n");
		sos_rho = 0.01f;
		gtk_spin_button_set_value(button, 0.01);
	}
}

void on_spin_denoise_modulation_value_changed(GtkSpinButton *button, gpointer user_data) {
	denoise_modulation = (float)gtk_spin_button_get_value(button);
}

void on_radio_denoise_nosecondary_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_vst = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_vst"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	do_anscombe = (gtk_toggle_button_get_active(toggle_vst) ? TRUE : FALSE);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
	if (sos == 1)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), FALSE);
}

void on_radio_denoise_da3d_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_vst = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_vst"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(button) ? 1 : 0);
	do_anscombe = (gtk_toggle_button_get_active(toggle_vst) ? TRUE : FALSE);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
	if (sos == 1)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), FALSE);
}

void on_radio_denoise_sos_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_vst = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_vst"));
	sos = (gtk_toggle_button_get_active(button) ? 1 : 0);
	do_anscombe = (gtk_toggle_button_get_active(toggle_vst) ? TRUE : FALSE);
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	if (sos == 1)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), FALSE);
}

void on_check_denoise_cosmetic_toggled(GtkToggleButton *button, gpointer user_data) {
	do_cosme = (gtk_toggle_button_get_active(button) ? TRUE : FALSE);
}

void on_radio_denoise_vst_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	do_anscombe = (gtk_toggle_button_get_active(button) ? TRUE : FALSE);
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
	if (sos == 1)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), FALSE);
}

void on_denoise_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = (float)gtk_spin_button_get_value(spin_denoise_modulation);
	gboolean suppress_artefacts = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_denoise_suppress_artefacts")));
	if (gui.roi.active)
		restore_roi();
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->modulation = denoise_modulation;
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		free(args);
		return;
	}
	control_window_switch_to_tab(OUTPUT_LOGS);
	args->previewing = ((GtkWidget*) button == lookup_widget("denoise_roi_preview"));
	siril_debug_print("Previewing: %d\n", args->previewing);
	args->da3d = da3d;
	args->sos = 1;
	args->rho = sos_rho;
	args->do_anscombe = do_anscombe;
	args->do_cosme = do_cosme;
	args->suppress_artefacts = suppress_artefacts;
	if (sos == 1)
		args->sos = sos_iters;
	args->fit = (gui.roi.active && args->previewing) ? &gui.roi.fit : &gfit;
	if (!start_in_new_thread(run_nlbayes_on_fit, args)) {
		free(args);
	}
}
