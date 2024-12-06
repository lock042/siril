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
#include "gui/graxpert.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

// Statics declarations

static GtkButton *button_graxpert_cancel = NULL, *button_graxpert_apply = NULL, *graxpert_clear_samples = NULL, *graxpert_generate_samples = NULL, *button_graxpert_roipreview = NULL;
static GtkNotebook *notebook_graxpert_operation = NULL;
static GtkComboBox *combo_graxpert_algorithm = NULL, *combo_graxpert_correction = NULL, *combo_graxpert_kernel = NULL, *combo_graxpert_ai_models_bg = NULL, *combo_graxpert_ai_models_denoise = NULL, *combo_graxpert_ai_models_deconv = NULL;
static GtkDialog *graxpert_dialog = NULL;
static GtkSpinButton *spin_graxpert_smoothing = NULL, *graxpert_spin_bgtol = NULL, *graxpert_spin_nb_samples = NULL, *spin_graxpert_strength = NULL, *spin_graxpert_deconv_strength = NULL, *spin_graxpert_deconv_blur_psf_size = NULL, *graxpert_spin_sample_size = NULL, *graxpert_spin_spline_order = NULL;
static GtkToggleButton *toggle_graxpert_gpu = NULL, *graxpert_toggle_keep_background = NULL, *graxpert_toggle_apply_to_sequence = NULL;
static GtkLabel *graxpert_available = NULL;
static GtkWidget *graxpert_ai_settings = NULL, *graxpert_classical_settings = NULL, *graxpert_samples_controls = NULL, *graxpert_rbf_settings = NULL, *graxpert_spline_settings = NULL, *ai_model_settings_bg = NULL, *ai_model_settings_denoise = NULL;

static gboolean is_bg = TRUE;

gboolean initialize_graxpert_widgets_if_needed(gpointer user_data) {
	int populate_ai_combos = GPOINTER_TO_INT(user_data);
	if (button_graxpert_cancel == NULL) {
		notebook_graxpert_operation = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_graxpert_operation"));
		// GtkButton
		button_graxpert_cancel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_cancel"));
		button_graxpert_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_apply"));
		graxpert_clear_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_clear_samples"));
		graxpert_generate_samples = GTK_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_generate_samples"));
		button_graxpert_roipreview = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_graxpert_roipreview"));
		// GtkComboBoxText
		combo_graxpert_algorithm = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_algorithm"));
		combo_graxpert_correction = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_correction"));
		combo_graxpert_kernel = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_kernel"));
		combo_graxpert_ai_models_bg = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_ai_models_bg"));
		combo_graxpert_ai_models_denoise = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_ai_models_denoise"));
		combo_graxpert_ai_models_deconv = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_graxpert_ai_models_deconv"));
		// GtkDialog
		graxpert_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "graxpert_dialog"));
		// GtkSpinButton
		spin_graxpert_smoothing = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_graxpert_smoothing"));
		graxpert_spin_bgtol = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_bgtol"));
		graxpert_spin_spline_order = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_spline_order"));
		graxpert_spin_nb_samples = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_nb_samples"));
		graxpert_spin_sample_size = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_spin_sample_size"));
		spin_graxpert_strength = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_graxpert_strength"));
		spin_graxpert_deconv_blur_psf_size = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_deconv_graxpert_psf"));
		spin_graxpert_deconv_strength = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_deconv_graxpert_strength"));
		// GtkToggleButton
		toggle_graxpert_gpu = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_graxpert_gpu"));
		graxpert_toggle_keep_background = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_toggle_keep_background"));
		graxpert_toggle_apply_to_sequence = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "graxpert_toggle_apply_to_sequence"));
		// GtkLabel
		graxpert_available = GTK_LABEL(gtk_builder_get_object(gui.builder, "graxpert_available"));
		// GtkWidget
		graxpert_samples_controls = GTK_WIDGET(gtk_builder_get_object(gui.builder, "graxpert_samples_controls"));
		graxpert_ai_settings = GTK_WIDGET(gtk_builder_get_object(gui.builder, "graxpert_ai_settings"));
		graxpert_classical_settings = GTK_WIDGET(gtk_builder_get_object(gui.builder, "graxpert_classical_settings"));
		graxpert_rbf_settings = GTK_WIDGET(gtk_builder_get_object(gui.builder, "graxpert_rbf_settings"));
		graxpert_spline_settings = GTK_WIDGET(gtk_builder_get_object(gui.builder, "graxpert_spline_settings"));
		ai_model_settings_bg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "ai_model_settings_bg"));
		ai_model_settings_denoise = GTK_WIDGET(gtk_builder_get_object(gui.builder, "ai_model_settings_denoise"));
	}
	if (populate_ai_combos) {
		populate_graxpert_ai_combos(NULL);
	}
	return FALSE;
}

graxpert_data* fill_graxpert_data_from_gui(gboolean previewing) {
	graxpert_data *p = calloc(1, sizeof(graxpert_data));
	p->previewing = previewing;
	p->fit = previewing ? &gui.roi.fit : &gfit;
	p->operation = (graxpert_operation) gtk_notebook_get_current_page(notebook_graxpert_operation);
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
	p->deconv_strength = gtk_spin_button_get_value(spin_graxpert_deconv_strength);
	p->deconv_blur_psf_size = gtk_spin_button_get_value(spin_graxpert_deconv_blur_psf_size);
	p->use_gpu = gtk_toggle_button_get_active(toggle_graxpert_gpu);
	p->ai_batch_size = 4;
	p->bg_pts_option = gtk_spin_button_get_value(graxpert_spin_nb_samples);
	if (p->operation == GRAXPERT_DENOISE) {
		int n = gtk_combo_box_get_active(combo_graxpert_ai_models_denoise);
		const gchar **ai_models = get_ai_models(GRAXPERT_DENOISE);
		if (ai_models) {
			int num_models = g_strv_length((gchar**) ai_models);
			if (n >= 0 && n < num_models) {
				p->ai_version = g_strdup(ai_models[n]);
			}
		}
	} else if (p->operation == GRAXPERT_BG && p->bg_algo == GRAXPERT_BG_AI) {
		int n = gtk_combo_box_get_active(combo_graxpert_ai_models_bg);
		const gchar **ai_models = get_ai_models(GRAXPERT_BG);
		if (ai_models) {
			int num_models = g_strv_length((gchar**) ai_models);
			if (n >= 0 && n < num_models) {
				p->ai_version = g_strdup(ai_models[n]);
			}
		}
	} else if (p->operation == GRAXPERT_DECONV) {
		int n = gtk_combo_box_get_active(combo_graxpert_ai_models_deconv);
		const gchar **ai_models = get_ai_models(GRAXPERT_DECONV);
		if (ai_models) {
			int num_models = g_strv_length((gchar**) ai_models);
			if (n >= 0 && n < num_models) {
				p->ai_version = g_strdup(ai_models[n]);
			}
		}
	}
	return p;
}

void graxpert_roi_callback() {
	// ROI not supported for GraXpert background removal
	gui.roi.operation_supports_roi = !is_bg;
	gtk_widget_set_visible(GTK_WIDGET(button_graxpert_roipreview), (!is_bg && gui.roi.active));
	if (is_preview_active())
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
		if (is_preview_active())
			siril_preview_hide();
		remove_roi_callback(graxpert_roi_callback);
		mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	}
}

static void confirm_availability(guint page_num) {
	graxpert_operation operation = (graxpert_operation) page_num;

	gboolean available = graxpert_executablecheck(com.pref.graxpert_path, operation);
	gtk_widget_set_sensitive(GTK_WIDGET(button_graxpert_apply), available);
	if (available) {
		gchar *txt = g_strdup_printf(_("GraXpert available and supports %s."),
				operation == GRAXPERT_BG ? _("background extraction") :
				operation == GRAXPERT_DENOISE ? _("denoising") :
				operation == GRAXPERT_DECONV ? _("deconvolution") :
		    _("GUI"));
		gtk_label_set_text(graxpert_available, txt);
		g_free(txt);
	} else {
		gtk_label_set_markup(graxpert_available, _("<span foreground=\"red\">No suitable version of GraXpert is available.\n"
				"Configure the executable in Preferences -> Miscellaneous.</span>"));
	}
}

static void populate_combo_box(GtkComboBoxText *combo, const gchar **models) {
	int i;

	// Clear existing entries
	gtk_combo_box_text_remove_all(combo);

	// Add entries from the models array
	for (i = 0; models[i] != NULL; ++i) {
		gtk_combo_box_text_append_text(combo, models[i]);
	}

	// Add "latest" entry
	gtk_combo_box_text_append_text(combo, "latest");

	// Set "latest" as the active (default) item
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), i);  // i is now the index of "latest"
}

gboolean populate_graxpert_ai_combos(gpointer user_data) {
	const gchar** ai_models_bg = get_ai_models(GRAXPERT_BG);
	if (combo_graxpert_ai_models_bg && ai_models_bg)
		populate_combo_box(GTK_COMBO_BOX_TEXT(combo_graxpert_ai_models_bg), ai_models_bg);
	const gchar** ai_models_denoise = get_ai_models(GRAXPERT_DENOISE);
	if (combo_graxpert_ai_models_denoise && ai_models_denoise)
		populate_combo_box(GTK_COMBO_BOX_TEXT(combo_graxpert_ai_models_denoise), get_ai_models(GRAXPERT_DENOISE));
	const gchar** ai_models_deconv = get_ai_models(GRAXPERT_DECONV);
	if (combo_graxpert_ai_models_deconv && ai_models_deconv)
		populate_combo_box(GTK_COMBO_BOX_TEXT(combo_graxpert_ai_models_deconv), get_ai_models(GRAXPERT_DECONV));
	return FALSE;
}

static void set_widgets() {
	graxpert_operation operation = (graxpert_operation) gtk_notebook_get_current_page(notebook_graxpert_operation);
	graxpert_bg_algo algorithm = gtk_combo_box_get_active(combo_graxpert_algorithm);
	gtk_widget_set_visible(graxpert_ai_settings, algorithm == GRAXPERT_BG_AI || operation == GRAXPERT_DENOISE || operation == GRAXPERT_DECONV);
	gtk_widget_set_visible(graxpert_classical_settings, algorithm != GRAXPERT_BG_AI);
	gtk_widget_set_visible(graxpert_samples_controls, algorithm != GRAXPERT_BG_AI);
	gtk_widget_set_visible(graxpert_rbf_settings, algorithm == GRAXPERT_BG_RBF ||  algorithm == GRAXPERT_BG_KRIGING);
	gtk_widget_set_visible(graxpert_spline_settings, algorithm == GRAXPERT_BG_SPLINE);
	is_bg = (operation == GRAXPERT_BG);
	gtk_widget_set_visible(GTK_WIDGET(button_graxpert_roipreview), (!is_bg && gui.roi.active));
	configure_graxpert_dialog_for_roi();
	redraw(REDRAW_OVERLAY);
}

void on_graxpert_dialog_show(GtkWidget *widget, gpointer user_data) {
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	initialize_graxpert_widgets_if_needed(GINT_TO_POINTER(0));
	set_widgets();
	confirm_availability((guint) GRAXPERT_BG);
	clear_backup();
}

void on_graxpert_dialog_hide(GtkWidget *widget, gpointer user_data) {
	roi_supported(FALSE);
	if (is_preview_active())
		siril_preview_hide();
	remove_roi_callback(graxpert_roi_callback);
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	siril_close_dialog("graxpert_dialog");
};

void on_notebook_graxpert_operation_switch_page(GtkNotebook *notebook, GtkWidget *page, guint page_num, gpointer user_data) {
	gboolean user_cancelled = FALSE;
	GPid pid = get_running_graxpert_pid();
	if (pid != (GPid) -1 && (user_cancelled = siril_confirm_dialog(_("Warning!"), _("GraXpert is running. Changing the GraXpert "
						"operation will cancel the current GraXpert process. Proceed?"), _("Yes")))) {
		kill_child_process(pid, FALSE);
		if (user_cancelled) {
			siril_log_color_message(_("GraXpert operation cancelled by user\n"), "red");
		}
	}
	set_widgets();
	confirm_availability(page_num);
}

void on_combo_graxpert_algorithm_changed(GtkComboBox *combo, gpointer user_data) {
	set_widgets();
}

void on_button_graxpert_apply_clicked(GtkButton *button, gpointer user_data) {
	gboolean previewing = (button == button_graxpert_roipreview);
	graxpert_data *data = fill_graxpert_data_from_gui(previewing);
	if (com.grad_samples) {
		data->bg_samples = com.grad_samples;
		com.grad_samples = NULL;
	}
	if (gtk_toggle_button_get_active(graxpert_toggle_apply_to_sequence)) {
		if (sequence_is_loaded()) {
			data->seq = &com.seq;
			apply_graxpert_to_sequence(data);
		} else {
			siril_log_color_message(_("Error: no sequence loaded.\n"), "red");
		}
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

void on_graxpert_deconv_switch_state_set(GtkSwitch *widget, gboolean state, gpointer user_data) {
	gtk_label_set_text((GtkLabel *) user_data, state ? _("Stellar") : _("Objects"));
}
