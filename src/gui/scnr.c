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

/* GTK callbacks and preview logic for the SCNR dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "filters/scnr.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "gui/utils.h"
#include "io/single_image.h"

static GtkToggleButton *scnr_roi_preview = NULL;
static GtkComboBox *scnr_combo = NULL;
static GtkScale *scnr_scale = NULL;
static GtkLabel *scnr_label56 = NULL;
static GtkSpinButton *scnr_spin = NULL;
static GtkToggleButton *scnr_preserve_light = NULL;

static void scnr_dialog_init_statics(void) {
	if (scnr_roi_preview) return;
	scnr_roi_preview = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "SCNR_roi_preview"));
	scnr_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_scnr"));
	scnr_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "scale_scnr"));
	scnr_label56 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label56"));
	scnr_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_scnr"));
	scnr_preserve_light = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "preserve_light"));
}

static int scnr_process_with_worker(scnr_type type, double amount, gboolean preserve,
                                    gboolean for_preview, gboolean for_roi) {
	struct scnr_data *params = new_scnr_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	params->type = type;
	params->amount = amount;
	params->preserve = preserve;
	params->verbose = !for_preview;
	params->applying = !for_preview;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_scnr_data(params);
		return 1;
	}

	args->fit = for_roi ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = scnr_image_hook;
	args->idle_function = NULL;
	args->description = _("Subtractive Chromatic Noise Reduction");
	args->verbose = !for_preview;
	args->user = params;
	args->log_hook = scnr_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = for_roi;
	args->populate_roi_on_complete = !for_preview;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}
	return 0;
}

static int scnr_update_preview(void) {
	if (gtk_toggle_button_get_active(scnr_roi_preview)) {
		int type = gtk_combo_box_get_active(scnr_combo);
		gboolean preserve = gtk_toggle_button_get_active(scnr_preserve_light);
		double amount = gtk_range_get_value(GTK_RANGE(scnr_scale));
		copy_backup_to_gfit();
		return scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	}
	return 0;
}

void scnr_change_between_roi_and_image(void) {
	gui.roi.operation_supports_roi = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = gtk_toggle_button_get_active(scnr_roi_preview);
	notify_update((gpointer) param);
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	scnr_dialog_init_statics();
	roi_supported(TRUE);
	gtk_widget_set_visible(GTK_WIDGET(scnr_roi_preview), gui.roi.active);

	copy_gfit_to_backup();
	add_roi_callback(scnr_change_between_roi_and_image);

	if (gui.roi.active) {
		scnr_change_between_roi_and_image();
	}

	int type = gtk_combo_box_get_active(scnr_combo);
	if (type == -1)
		gtk_combo_box_set_active(scnr_combo, 0);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	roi_supported(FALSE);

	if (gtk_toggle_button_get_active(scnr_roi_preview)) {
		copy_backup_to_gfit();
		gfit_modified_update_gui();
	}

	clear_backup();
	remove_roi_callback(scnr_change_between_roi_and_image);
	siril_close_dialog("SCNR_dialog");
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	int type = gtk_combo_box_get_active(scnr_combo);
	gboolean preserve = gtk_toggle_button_get_active(scnr_preserve_light);
	double amount = gtk_range_get_value(GTK_RANGE(scnr_scale));

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	gboolean is_preview = ((GtkWidget*) button == GTK_WIDGET(scnr_roi_preview));

	if (is_preview) {
		scnr_process_with_worker(type, amount, preserve, TRUE, gui.roi.active);
	} else {
		if (gtk_toggle_button_get_active(scnr_roi_preview)) {
			copy_backup_to_gfit();
		}

		set_cursor_waiting(TRUE);
		scnr_process_with_worker(type, amount, preserve, FALSE, FALSE);

		clear_backup();
		remove_roi_callback(scnr_change_between_roi_and_image);
		roi_supported(FALSE);
		siril_close_dialog("SCNR_dialog");
	}
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(scnr_combo);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_label56), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(scnr_spin), type > 1);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = gtk_toggle_button_get_active(scnr_roi_preview);
	notify_update((gpointer) param);
}

void on_SCNR_parameter_changed(GtkWidget *widget, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = scnr_update_preview;
	param->show_preview = gtk_toggle_button_get_active(scnr_roi_preview);
	notify_update((gpointer) param);
}

void on_SCNR_roi_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		cancel_and_wait_for_preview();
		siril_preview_hide();
		copy_backup_to_gfit();
		gfit_modified_update_gui();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = scnr_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
