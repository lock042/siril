#ifndef SIRIL_CURVES_H
#define SIRIL_CURVES_H

#include "filters/curve_transform.h"

struct curve_data {
	fits *fit;
	sequence *seq;
	struct curve_params params;
	char *seq_entry;
};

void curves_histogram_change_between_roi_and_image();

void update_gfit_curves_histogram_if_needed();

void apply_curves_cancel();

void toggle_curves_window_visibility();

void on_curves_close_button_clicked(GtkButton *button, gpointer user_data); // callback needed
void on_curves_display_toggle(GtkToggleButton *togglebutton, gpointer user_data);

void apply_curve_to_sequence(struct curve_data *curve_args);

void erase_curves_histogram_display(cairo_t *cr, int width, int height);

#endif