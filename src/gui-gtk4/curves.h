#ifndef SIRIL_CURVES_H
#include <gtk/gtk.h>
#define SIRIL_CURVES_H

#include "filters/curve_transform.h"

struct curve_data {
	fits *fit;
	sequence *seq;
	struct curve_params params;
	char *seq_entry;
};

// Snapshot structure for Undo Stack
typedef struct {
	curve_channel_config channels[CHAN_COUNT];
	enum curve_algorithm algorithm;
} curve_state_snapshot;

void curves_histogram_change_between_roi_and_image();
void update_gfit_curves_histogram_if_needed();
void curves_reset_after_undo();
void apply_curves_cancel();
void toggle_curves_window_visibility();
void on_curves_close_button_clicked(GtkButton *button, gpointer user_data);
void on_curves_display_toggle(GtkToggleButton *togglebutton, gpointer user_data);
void apply_curve_to_sequence(struct curve_data *curve_args);

void on_curve_check_range_button_toggled(GtkToggleButton *button, gpointer user_data);
void on_curves_range_value_changed(GtkRange *range, gpointer user_data);
void on_curves_feather_value_changed(GtkRange *range, gpointer user_data);

void curves_handle_pipette_click(int x, int y);

#endif
