#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
#include "filters/median.h"
#include "filters/bm3d/call_bm3d.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "core/processing.h"
#include "core/command.h"
#include "io/single_image.h"


// Callbacks
float bm3d_modulation;

void on_bm3d_dialog_show(GtkWidget *widget, gpointer user_data) {
  GtkSpinButton *spin_bm3d_modulation = GTK_SPIN_BUTTON(lookup_widget("spin_bm3d_modulation"));
  bm3d_modulation = 1.f;
  gtk_spin_button_set_value(spin_bm3d_modulation, bm3d_modulation);
}

void on_bm3d_cancel_clicked(GtkButton *button, gpointer user_data) {
  siril_close_dialog("bm3d_dialog");
}

void on_spin_bm3d_modulation_value_changed(GtkSpinButton *button, gpointer user_data) {
  bm3d_modulation = (float) gtk_spin_button_get_value(button);
}

void on_bm3d_apply_clicked(GtkButton *button, gpointer user_data) {
  bm3d_args *args;
  args = calloc(1, sizeof(bm3d_args));
  args->fit = &gfit;
  args->modulation = bm3d_modulation;
  if (get_thread_run()) {
    PRINT_ANOTHER_THREAD_RUNNING;
    return;
  }
  start_in_new_thread(run_bm3d_on_fit, args);
  siril_add_idle(end_generic, NULL);
}
