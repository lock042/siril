#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
#include "filters/median.h"
#include "filters/bm3d/call_bm3d.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/dialogs.h"
#include "core/processing.h"
#include "core/command.h"
#include "io/single_image.h"


// Callbacks
float bm3d_modulation;
int da3d;

void on_bm3d_dialog_show(GtkWidget *widget, gpointer user_data) {
	da3d = 0;
	GtkSpinButton *spin_bm3d_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_bm3d_modulation"));
	bm3d_modulation = 1.f;
	GtkToggleButton *toggle_bm3d_da3d = GTK_TOGGLE_BUTTON(lookup_widget("toggle_bm3d_da3d"));
	gtk_spin_button_set_value(spin_bm3d_modulation, bm3d_modulation);
	gtk_toggle_button_set_active(toggle_bm3d_da3d, da3d);
}

void on_bm3d_cancel_clicked(GtkButton *button, gpointer user_data) {
  siril_close_dialog("bm3d_dialog");
}

void on_spin_bm3d_modulation_value_changed(GtkSpinButton *button, gpointer user_data) {
  bm3d_modulation = (float) gtk_spin_button_get_value(button);
}

void on_toggle_bm3d_da3d_toggled(GtkToggleButton *button, gpointer user_data) {
	da3d = gtk_toggle_button_get_active(button);
}

void on_bm3d_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_bm3d_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_bm3d_modulation"));
	bm3d_modulation = (float) gtk_spin_button_get_value(spin_bm3d_modulation);
//	copy_gfit_to_backup();
	bm3d_args *args = calloc(1, sizeof(bm3d_args));
	args->fit = &gfit;
	args->da3d = da3d;
	args->modulation = bm3d_modulation;
	if (args->modulation == 0.f) {
		siril_log_message(_("Modulation is zero: doing nothing.\n"));
		return;
	}
	unsigned npixels = (unsigned) gfit.naxes[0] * gfit.naxes[1];
	float memGB = (float) (get_available_memory() / 1000000000);
	float imgmemMpix = (float) npixels / 1000000.f;
	unsigned numchunks = (unsigned) (imgmemMpix / (memGB / 5));
	if (numchunks < 1)
		numchunks = 1;
	siril_log_message(_("Available memory: %f GB, processing in %u chunks.\n"), memGB, numchunks);
	siril_log_message(_("Modulation: %f\n"),args->modulation);
	if (args->da3d)
		siril_log_message(_("Will carry out final stage DA3D denoising.\n"));
	else
		siril_log_message(_("Final stage DA3D denoising disabled.\n"));

	start_in_new_thread(run_bm3d_on_fit, args);
	siril_close_dialog("bm3d_dialog");
}
