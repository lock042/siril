#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
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
int da3d;
int sos = 1;
int sos_iters;
float sos_rho;

void on_denoise_dialog_show(GtkWidget *widget, gpointer user_data) {
	da3d = 0;
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = 1.f;
	GtkToggleButton *toggle_bm3d_da3d = GTK_TOGGLE_BUTTON(lookup_widget("toggle_denoise_da3d"));
	gtk_spin_button_set_value(spin_denoise_modulation, denoise_modulation);
}

void on_bm3d_cancel_clicked(GtkButton *button, gpointer user_data) {
  siril_close_dialog("bm3d_dialog");
}

void on_spin_sos_iters_value_changed(GtkSpinButton *button, gpointer user_data) {
  sos_iters = (float) gtk_spin_button_get_value(button);
}
void on_spin_sos_1mrho_value_changed(GtkSpinButton *button, gpointer user_data) {
  sos_rho = 1.f - (float) gtk_spin_button_get_value(button);
}

void on_spin_denoise_modulation_value_changed(GtkSpinButton *button, gpointer user_data) {
  denoise_modulation = (float) gtk_spin_button_get_value(button);
}

void on_radio_denoise_nosecondary_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
}

void on_radio_denoise_da3d_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
}
void on_radio_denoise_sos_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *toggle_da3d = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_da3d"));
	GtkToggleButton *toggle_sos = GTK_TOGGLE_BUTTON(lookup_widget("radio_denoise_sos"));
	da3d = (gtk_toggle_button_get_active(toggle_da3d) ? 1 : 0);
	sos = (gtk_toggle_button_get_active(toggle_sos) ? 1 : 0);
}
void on_denoise_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_denoise_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_denoise_modulation"));
	denoise_modulation = (float) gtk_spin_button_get_value(spin_denoise_modulation);
//	copy_gfit_to_backup();
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->fit = &gfit;
	args->da3d = da3d;
	args->sos = 1;
	if (sos == 1)
		args->sos = sos_iters;
	args->modulation = denoise_modulation;
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		return;
	}
	siril_log_message(_("Modulation: %f\n"),args->modulation);
	if (args->da3d)
		siril_log_message(_("Will carry out final stage DA3D denoising.\n"));
	else
		siril_log_message(_("Final stage DA3D denoising disabled.\n"));

	start_in_new_thread(run_nlbayes_on_fit, args);
	siril_close_dialog("denoise_dialog");
}
