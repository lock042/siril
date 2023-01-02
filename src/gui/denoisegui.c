/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	da3d = 0;
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = 1.f;
	gtk_spin_button_set_value(spin_denoise_modulation, denoise_modulation);
	if (gfit.naxes[2] == 3)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("denoise_artefact_control")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("denoise_artefact_control")), FALSE);
}

void on_denoise_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("denoise_dialog");
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
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = (float)gtk_spin_button_get_value(spin_denoise_modulation);
	//	copy_gfit_to_backup();
	gboolean suppress_artefacts = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("check_denoise_suppress_artefacts")));
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->fit = &gfit;
	args->da3d = da3d;
	args->sos = 1;
	args->rho = sos_rho;
	args->do_anscombe = do_anscombe;
	args->do_cosme = do_cosme;
	args->suppress_artefacts = suppress_artefacts;
	if (sos == 1)
		args->sos = sos_iters;
	args->modulation = denoise_modulation;
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		free(args);
		return;
	}
	char *msg1 = NULL, *msg2 = NULL, *msg3 = NULL, *log_msg = NULL;
	int n = 0, m = 0, q = 0;
	n = snprintf(NULL, 0, _("NL-Bayes denoise (mod=%.3f"), args->modulation);
	msg1 = malloc(n + 1);
	n = snprintf(msg1, n + 1, _("NL-Bayes denoise (mod=%.3f"), args->modulation);
	if (args->da3d) {
		m = snprintf(NULL, 0, _(", DA3D enabled"));
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", DA3D enabled"));
	} else if (args->sos > 1) {
		m = snprintf(NULL, 0, _(", SOS enabled (iters=%d, rho=%.3f)"), args->sos, args->rho);
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", SOS enabled (iters=%d, rho=%.3f)"), args->sos, args->rho);
	} else if (args->do_anscombe) {
		m = snprintf(NULL, 0, _(", VST enabled"));
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", VST enabled"));
	}
	if (args->do_cosme) {
		q = snprintf(NULL, 0, _(", CC enabled)"));
		msg3 = malloc(q + 1);
		q = snprintf(msg3, q + 1, _(", CC enabled)"));
	} else {
		q = 1;
		msg3 = malloc(q + 1);
		q = snprintf(msg3, q + 1, _(")"));
	}
	log_msg = malloc(n + m + q + 1);
	if (m == 0 && q == 0)
		snprintf(log_msg, n + 1, "%s", msg1);
	else if (m > 0 && q == 0)
		snprintf(log_msg, n + m + 1, "%s%s", msg1, msg2);
	else if (m == 0 && q > 0)
		snprintf(log_msg, n + q + 1, "%s%s", msg1, msg3);
	else if (m > 0 && q > 0)
		snprintf(log_msg, n + m + q + 1, "%s%s%s", msg1, msg2, msg3);
	else
		snprintf(log_msg, 26, "Error, this can't happen!");

	if (msg1) free(msg1);
	if (msg2) free(msg2);
	if (msg3) free(msg3);
	undo_save_state(&gfit, "%s", log_msg);
	free(log_msg);

	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(run_nlbayes_on_fit, args);
}
