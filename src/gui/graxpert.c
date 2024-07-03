/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/icc_profile.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/background_extraction.h"
#include "filters/graxpert.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"

// Statics declarations

GtkButton *button_graxpert_cancel = NULL, *button_graxpert_apply = NULL, *graxpert_clear_samples = NULL, *graxpert_generate_samples = NULL;
GtkComboBox *combo_graxpert_operation = NULL, *combo_graxpert_algorithm = NULL, *combo_graxpert_correction = NULL, *combo_graxpert_kernel = NULL;
GtkDialog *graxpert_dialog = NULL;
GtkExpander *graxpert_bg_settings = NULL, *graxpert_denoise_settings = NULL;
GtkSpinButton *spin_graxpert_smoothing = NULL, *graxpert_spin_bgtol = NULL, *graxpert_spin_nb_samples = NULL, *spin_graxpert_strength = NULL, *graxpert_spin_sample_size = NULL, *graxpert_spin_spline_order = NULL;
GtkToggleButton *toggle_graxpert_gpu = NULL, *graxpert_toggle_keep_background = NULL, *graxpert_toggle_apply_to_sequence = NULL;

void initialize_graxpert_widgets_if_needed() {
	if (button_graxpert_cancel == NULL) {
		// GtkButton
		button_graxpert_cancel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_cancel"));
		button_graxpert_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_apply"));
		graxpert_clear_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_clear_samples"));
		graxpert_generate_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_generate_samples"));
		// GtkComboBoxText
		combo_graxpert_operation = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_operation"));
		combo_graxpert_algorithm = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_algorithm"));
		combo_graxpert_correction = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_correction"));
		combo_graxpert_kernel = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_kernel"));
		// GtkDialog
		graxpert_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "graxpert_dialog"));
		// GtkExpander
		graxpert_bg_settings = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "graxpert_bg_settings"));
		graxpert_denoise_settings = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "graxpert_denoise_settings"));
		// GtkSpinButton
		spin_graxpert_smoothing = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_graxpert_smoothing"));
		graxpert_spin_bgtol = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_bgtol"));
		graxpert_spin_spline_order = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_spline_order"));
		graxpert_spin_nb_samples = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_nb_samples"));
		graxpert_spin_sample_size = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_sample_size"));
		spin_graxpert_strength = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_graxpert_strength"));
		// GtkToggleButton
		toggle_graxpert_gpu = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_graxpert_gpu"));
		graxpert_toggle_keep_background = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_toggle_keep_background"));
		graxpert_toggle_apply_to_sequence = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_toggle_apply_to_sequence"));
	}
}

graxpert_data* fill_graxpert_data_from_gui() {
	graxpert_data *p = calloc(1, sizeof(graxpert_data));
	p->fit = &gfit;
	p->operation = (graxpert_operation) gtk_combo_box_get_active(combo_graxpert_operation);
	p->bg_smoothing = gtk_spin_button_get_value(spin_graxpert_smoothing);
	p->bg_algo = (graxpert_bg_algo) gtk_combo_box_get_active(combo_graxpert_algorithm);
	p->bg_mode = (graxpert_bg_mode) gtk_combo_box_get_active(combo_graxpert_correction);
	p->kernel = (graxpert_rbf_kernel) gtk_combo_box_get_active(combo_graxpert_kernel);
	p->stretch_option = STRETCH_OPTION_NONE; // Doesn't really matter, this is only for the GUI
	p->sample_size = gtk_spin_button_get_value(graxpert_spin_sample_size);
	p->spline_order = gtk_spin_button_get_value(graxpert_spin_spline_order);
	p->bg_tol_option = gtk_spin_button_get_value(graxpert_spin_bgtol);
	p->keep_bg = gtk_toggle_button_get_active(graxpert_toggle_keep_background);
	p->denoise_strength = gtk_spin_button_get_value(spin_graxpert_strength);
	p->use_gpu = gtk_toggle_button_get_active(toggle_graxpert_gpu);
	p->ai_batch_size = 4;
	p->bg_pts_option = gtk_spin_button_get_value(graxpert_spin_nb_samples);
	return p;
}

void on_graxpert_generate_samples_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	int nb_of_samples = (int) gtk_spin_button_get_value(graxpert_spin_nb_samples);
	int sample_size = (int) gtk_spin_button_get_value(graxpert_spin_sample_size);
	double tolerance = gtk_spin_button_get_value(graxpert_spin_bgtol);

	// Use the sample size here. We don't care about the statistics (so we don't
	// call update_median_samples), but it gives a visual indication of the size
	// of the samples that GraXpert will use.
	free_background_sample_list(com.grad_samples);
	const char *err;
	com.grad_samples = generate_samples(&gfit, nb_of_samples, tolerance, sample_size, &err, TRUE);
	if (!com.grad_samples) {
		siril_log_color_message(_("Failed to generate background samples for image: %s\n"), "red", _(err));
		return;
	}
	control_window_switch_to_tab(OUTPUT_LOGS);
	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
}

void on_graxpert_dialog_show(GtkWidget *widget, gpointer user_data) {
	initialize_graxpert_widgets_if_needed();
}

void on_graxpert_dialog_hide(GtkWidget *widget, gpointer user_data) {
	siril_close_dialog("graxpert_dialog");
};

void on_combo_graxpert_operation_changed(GtkComboBox *combo, gpointer user_data) {
	graxpert_operation operation = (graxpert_operation) gtk_combo_box_get_active(combo);
	gtk_expander_set_expanded(graxpert_bg_settings, (operation == GRAXPERT_BG));
	gtk_expander_set_expanded(graxpert_denoise_settings, (operation == GRAXPERT_DENOISE));
}

void on_button_graxpert_apply_clicked(GtkWidget *widget, gpointer user_data) {
	graxpert_data *data = fill_graxpert_data_from_gui();
	// Undo is not possible, as the result is read as a new image
	if (gtk_toggle_button_get_active(graxpert_toggle_apply_to_sequence)) {
		data->seq = &com.seq;
		apply_graxpert_to_sequence(data);
	} else {
		start_in_new_thread(do_graxpert, data);
	}
}

void on_button_graxpert_cancel_clicked(GtkWidget *widget, gpointer user_data) {
	siril_close_dialog("graxpert_dialog");
}

void on_graxpert_clear_samples_clicked(GtkWidget *widget, gpointer user_data) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	redraw(REDRAW_OVERLAY);
}

void on_graxpert_spin_sample_size_value_changed(GtkSpinButton *button, gpointer user_data) {
	int sample_size = (int) gtk_spin_button_get_value(button);
	if (!(sample_size % 2)) // Must be odd
		sample_size++;
	g_signal_handlers_block_by_func(button, on_graxpert_spin_sample_size_value_changed, NULL);
	gtk_spin_button_set_value(button, sample_size);
	g_signal_handlers_unblock_by_func(button, on_graxpert_spin_sample_size_value_changed, NULL);
}
