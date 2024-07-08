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
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"

// Statics declarations

static GtkButton *button_graxpert_cancel = NULL, *button_graxpert_apply = NULL, *graxpert_clear_samples = NULL, *graxpert_generate_samples = NULL, *button_graxpert_roipreview = NULL;
static GtkComboBox *combo_graxpert_operation = NULL, *combo_graxpert_algorithm = NULL, *combo_graxpert_correction = NULL, *combo_graxpert_kernel = NULL;
static GtkDialog *graxpert_dialog = NULL;
static GtkExpander *graxpert_bg_settings = NULL, *graxpert_denoise_settings = NULL;
static GtkSpinButton *spin_graxpert_smoothing = NULL, *graxpert_spin_bgtol = NULL, *graxpert_spin_nb_samples = NULL, *spin_graxpert_strength = NULL, *graxpert_spin_sample_size = NULL, *graxpert_spin_spline_order = NULL;
static GtkToggleButton *toggle_graxpert_gpu = NULL, *graxpert_toggle_keep_background = NULL, *graxpert_toggle_apply_to_sequence = NULL;

static gboolean is_bg = TRUE;
static graxpert_operation previous_operation = GRAXPERT_BG;

void initialize_graxpert_widgets_if_needed() {
	if (button_graxpert_cancel == NULL) {
		// GtkButton
		button_graxpert_cancel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_cancel"));
		button_graxpert_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_apply"));
		graxpert_clear_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_clear_samples"));
		graxpert_generate_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_generate_samples"));
		button_graxpert_roipreview = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_roipreview"));
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

graxpert_data* fill_graxpert_data_from_gui(gboolean previewing) {
	graxpert_data *p = calloc(1, sizeof(graxpert_data));
	p->previewing = previewing;
	p->fit = previewing ? &gui.roi.fit : &gfit;
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

void graxpert_roi_callback() {
	// ROI not supported for GraXpert background removal
	gui.roi.operation_supports_roi = !is_bg;
	gtk_widget_set_visible(GTK_WIDGET(button_graxpert_roipreview), (!is_bg && gui.roi.active));
	copy_backup_to_gfit();
	notify_gfit_modified();
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

void configure_graxpert_dialog_for_roi() {
	if (!is_bg) {
		roi_supported(TRUE);
		graxpert_roi_callback();
		add_roi_callback(graxpert_roi_callback);
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
		copy_gfit_to_backup();
	} else {
		roi_supported(FALSE);
		siril_preview_hide();
		remove_roi_callback(graxpert_roi_callback);
		mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	}
}

void on_graxpert_dialog_show(GtkWidget *widget, gpointer user_data) {
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	initialize_graxpert_widgets_if_needed();
	configure_graxpert_dialog_for_roi();
}

void on_graxpert_dialog_hide(GtkWidget *widget, gpointer user_data) {
	roi_supported(FALSE);
	siril_preview_hide();
	remove_roi_callback(graxpert_roi_callback);
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	siril_close_dialog("graxpert_dialog");
};

void on_combo_graxpert_operation_changed(GtkComboBox *combo, gpointer user_data) {
	gboolean user_cancelled = FALSE;
	if ((com.child_is_running == EXT_NONE) || (com.child_is_running == EXT_GRAXPERT &&
		(user_cancelled = siril_confirm_dialog(_("Warning!"), _("GraXpert is running. Changing the GraXpert operation will cancel the current GraXpert process. Proceed?"), _("Yes"))))) {
		kill_child_process(FALSE);
		if (user_cancelled) {
			siril_log_color_message(_("GraXpert operation cancelled by user\n"), "red");
		}
		graxpert_operation operation = (graxpert_operation) gtk_combo_box_get_active(combo);
		previous_operation = operation;
		gtk_expander_set_expanded(graxpert_bg_settings, (operation == GRAXPERT_BG));
		gtk_widget_set_sensitive(GTK_WIDGET(graxpert_bg_settings), (operation == GRAXPERT_BG));
		gtk_expander_set_expanded(graxpert_denoise_settings, (operation == GRAXPERT_DENOISE));
		gtk_widget_set_sensitive(GTK_WIDGET(graxpert_denoise_settings), (operation == GRAXPERT_DENOISE));
		is_bg = (operation == GRAXPERT_BG);
		gtk_widget_set_visible(GTK_WIDGET(button_graxpert_roipreview), (!is_bg && gui.roi.active));
		configure_graxpert_dialog_for_roi();
		redraw(REDRAW_OVERLAY);
	} else {
		g_signal_handlers_block_by_func(combo, on_combo_graxpert_operation_changed, NULL);
		gtk_combo_box_set_active(combo, (gint) previous_operation);
		g_signal_handlers_unblock_by_func(combo, on_combo_graxpert_operation_changed, NULL);
	}
}

void on_button_graxpert_apply_clicked(GtkButton *button, gpointer user_data) {
	gboolean previewing = (button == button_graxpert_roipreview);
	graxpert_data *data = fill_graxpert_data_from_gui(previewing);
	if (com.grad_samples) {
		data->bg_samples = com.grad_samples;
		com.grad_samples = NULL;
	}
	if (gtk_toggle_button_get_active(graxpert_toggle_apply_to_sequence)) {
		data->seq = &com.seq;
		apply_graxpert_to_sequence(data);
	} else {
		start_in_new_thread(do_graxpert, data);
	}
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
