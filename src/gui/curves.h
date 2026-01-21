/*
* This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifndef SIRIL_CURVES_H
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

void on_curves_channel_combo_changed(GtkComboBox *widget, gpointer user_data);
void on_curve_check_range_button_toggled(GtkToggleButton *button, gpointer user_data);
void on_curves_range_value_changed(GtkRange *range, gpointer user_data);
void on_curves_feather_value_changed(GtkRange *range, gpointer user_data);

#endif