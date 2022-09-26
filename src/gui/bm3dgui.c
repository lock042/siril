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

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	da3d = 0;
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = 1.f;
	gtk_spin_button_set_value(spin_denoise_modulation, denoise_modulation);
}

void on_denoise_cancel_clicked(GtkButton *button, gpointer user_data) {
  siril_close_dialog("denoise_dialog");
}

void on_spin_sos_iters_value_changed(GtkSpinButton *button, gpointer user_data) {
  sos_iters = (float) gtk_spin_button_get_value(button);
}
void on_spin_rho_value_changed(GtkSpinButton *button, gpointer user_data) {
  sos_rho = (float) gtk_spin_button_get_value(button);
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
  denoise_modulation = (float) gtk_spin_button_get_value(button);
}
/*
void on_radio_denoise_nosecondary_group_changed(GtkWidget *widget, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
	if (sos ==1)
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("sos_advanced_options")), FALSE);
}
*/
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
	denoise_modulation = (float) gtk_spin_button_get_value(spin_denoise_modulation);
//	copy_gfit_to_backup();
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->fit = &gfit;
	args->da3d = da3d;
	args->sos = 1;
	args->rho = sos_rho;
	args->do_anscombe = do_anscombe;
	args->do_cosme = do_cosme;
	if (sos == 1)
		args->sos = sos_iters;
	args->modulation = denoise_modulation;
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		free(args);
		return;
	}
	siril_log_message(_("Modulation: %f\n"),args->modulation);
	if (args->do_cosme)
		siril_log_message(_("Will apply salt and pepper noise removal.\n"));
	else
		siril_log_message(_("Salt and pepper noise removal disabled.\n"));
	if (args->do_anscombe)
		siril_log_message(_("Will apply Anscombe variance stabilising transform.\n"));
	else
		siril_log_message(_("Anscombe variance stabilising transform disabled.\n"));
	if (args->da3d)
		siril_log_message(_("Will apply final stage DA3D denoising.\n"));
	else
		siril_log_message(_("Final stage DA3D denoising disabled.\n"));
	if (sos)
		siril_log_message(_("Will apply SOS iterative denoise booster.\n"));
	else
		siril_log_message(_("SOS denoise booster disabled.\n"));
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(run_nlbayes_on_fit, args);
	siril_close_dialog("denoise_dialog");
}
