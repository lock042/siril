/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/epf.h"
#include "gui/epf.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"
#include "opencv/opencv.h"

static float epf_d_value = 0.0f, epf_sigma_col_value = 11.0f, epf_sigma_spatial_value = 11.0f, mod = 1.f;
static ep_filter_t filter_type = EP_BILATERAL;
static fits *guide = NULL, loaded_fit = { 0 };

static int epf_update_preview() {
	gboolean guide_needs_freeing = FALSE;
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview"))))
		copy_backup_to_gfit();
	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;
	if (filter_type == EP_GUIDED) {
		if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")))) {
			guide = fit;
		} else {
			if (loaded_fit.rx != 0) {
				guide = &loaded_fit;
			} else {
				return 1;
			}
		}
	} else {
		guide = NULL;
	}
	struct epfargs *args = calloc(1, sizeof(struct epfargs));
	*args = (struct epfargs) {.fit = fit, .guidefit = guide, .d = epf_d_value, .sigma_col = epf_sigma_col_value,
								.sigma_space = epf_sigma_spatial_value, .mod = mod, .filter = filter_type,
								.guide_needs_freeing = guide_needs_freeing, .verbose = FALSE };
	set_cursor_waiting(TRUE);
	// We call epf_filter here as update_preview already handles the ROI mutex lock
	start_in_new_thread(epf_filter, args);
	return 0;
}

void epf_change_between_roi_and_image() {
	// If we are showing the preview, update it after the ROI change.
	roi_supported(TRUE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

static void epf_startup() {
	copy_gfit_to_backup();
	add_roi_callback(epf_change_between_roi_and_image);
	roi_supported(TRUE);
}

static void epf_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Bilateral filter: (d=%2.2f, sigma_col=%2.2f, sigma_spatial=%2.2f)"),
				epf_d_value, epf_sigma_col_value, epf_sigma_spatial_value);
	}
	backup_roi();
	roi_supported(FALSE);
	remove_roi_callback(epf_change_between_roi_and_image);
	clearfits(&loaded_fit);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int epf_process_all() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview"))))
		copy_backup_to_gfit();
	fits *fit = &gfit;
	gboolean guide_needs_freeing = FALSE;
	struct epfargs *args = calloc(1, sizeof(struct epfargs));
	if (filter_type != EP_BILATERAL) {
		if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")))) {
			guide = fit;
		} else {
			if (loaded_fit.rx != 0) {
				guide = &loaded_fit;
				guide_needs_freeing = TRUE;
			} else {
				siril_log_color_message(_("Error, no guide image loaded."), "red");
				free(args);
				return 1;
			}
		}
	} else {
		guide = NULL;
	}
	*args = (struct epfargs) {	.fit = fit, .guidefit = guide, .d = epf_d_value, .sigma_col = epf_sigma_col_value,
								.sigma_space = epf_sigma_spatial_value, .mod = mod, .filter = filter_type,
								.guide_needs_freeing = guide_needs_freeing, .verbose = FALSE };
	// We call epfhandler here as we need to take care of the ROI mutex lock
	start_in_new_thread(epfhandler, args);
	populate_roi();
	notify_gfit_modified();
	return 0;
}

static void apply_epf_changes() {
	gboolean status = (epf_sigma_col_value != 0.0f) || (epf_sigma_spatial_value != 0.0f);
	epf_close(!status);
}

void apply_epf_cancel() {
	epf_close(TRUE);
	siril_close_dialog("epf_dialog");
}

/*** callbacks **/

void on_epf_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_epf_d"));
	GtkSpinButton *spin_epf_sigma_col = GTK_SPIN_BUTTON(lookup_widget("spin_epf_sigma_col"));
	GtkSpinButton *spin_epf_sigma_spatial = GTK_SPIN_BUTTON(lookup_widget("spin_epf_sigma_spatial"));

//	if (gui.rendering_mode == LINEAR_DISPLAY)
//		setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535

	epf_startup();
	epf_d_value = 0.0f;
	epf_sigma_col_value = 11.0f;
	epf_sigma_spatial_value = 11.0f;
	filter_type = EP_BILATERAL;
	clearfits(&loaded_fit);

	set_notify_block(TRUE);
	gtk_widget_set_visible(lookup_widget("guide_image_widgets"), FALSE);
	gtk_widget_set_visible(lookup_widget("epf_sigma_spatial_settings"), TRUE);
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("guided_filter_guideimage")));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")), TRUE);
	gtk_widget_set_sensitive(lookup_widget("guided_filter_selfguide"), FALSE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("ep_filter_type")), EP_BILATERAL);
	gtk_spin_button_set_value(spin_d, epf_d_value);
	gtk_spin_button_set_value(spin_epf_sigma_spatial, epf_sigma_spatial_value);
	gtk_spin_button_set_value(spin_epf_sigma_col, epf_sigma_col_value);
	set_notify_block(FALSE);

	// Default parameters transform the image, so update the preview if toggle is active
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);

}

void on_epf_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_epf_cancel();
}

void on_epf_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview"))) || gui.roi.active) {
		epf_process_all();
	}

	apply_epf_changes();

	siril_close_dialog("epf_dialog");
}

void on_epf_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_epf_changes();
}

void on_epf_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_epf_d")), 0);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_epf_sigma_col")), 11);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_epf_sigma_spatial")), 11);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_epf_mod")), 1);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("ep_filter_type")), EP_BILATERAL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_ep_filter_type_changed(GtkComboBox *combo, gpointer user_data) {
	filter_type = gtk_combo_box_get_active(combo);
	gtk_widget_set_visible(lookup_widget("guide_image_widgets"), filter_type != EP_BILATERAL);
	gtk_widget_set_visible(lookup_widget("epf_sigma_spatial_settings"), filter_type == EP_BILATERAL);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_guided_filter_selfguide_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);
	if (active) {
		clearfits(&loaded_fit);
		gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("guided_filter_guideimage")));
		gtk_widget_set_sensitive(lookup_widget("guided_filter_selfguide"), FALSE);

	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_guided_filter_guideimage_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	gchar *filename = siril_file_chooser_get_filename(filechooser);
	clearfits(&loaded_fit);
	if (readfits(filename, &loaded_fit, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		clearfits(&loaded_fit);
		return;
	}
	if (loaded_fit.rx != gfit.rx || loaded_fit.ry != gfit.ry) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image dimensions do not match"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		clearfits(&loaded_fit);
		return;
	}
	gtk_widget_set_sensitive(lookup_widget("guided_filter_selfguide"), TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")), FALSE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_spin_epf_d_value_changed(GtkSpinButton *button, gpointer user_data) {
	epf_d_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_spin_epf_sigma_col_value_changed(GtkSpinButton *button, gpointer user_data) {
	epf_sigma_col_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_spin_epf_mod_value_changed(GtkSpinButton *button, gpointer user_data) {
	mod = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_spin_epf_sigma_spatial_value_changed(GtkSpinButton *button, gpointer user_data) {
	epf_sigma_spatial_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")));
	notify_update((gpointer) param);
}

void on_epf_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("epf_preview")))) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = epf_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
