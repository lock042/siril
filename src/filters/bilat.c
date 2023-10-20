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
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"
#include "opencv/opencv.h"

#include "bilat.h"

static float bilat_d_value = 0.0f, bilat_sigma_col_value = 0.0f, bilat_sigma_spatial_value = 0.0f;
static ep_filter_t filter_type = EP_BILATERAL;
static fits *guide = NULL, *loaded_fit = NULL;

int match_guide_to_roi(fits *guide, fits *guide_roi) {
	int retval = 0;
	if (!gui.roi.active)
		return 1;
	uint32_t nchans = guide->naxes[2];
	size_t npixels = guide->rx * guide->ry;
	gboolean rgb = (nchans == 3);
	size_t npixels_roi = gui.roi.fit.rx * gui.roi.fit.ry;
	copyfits(guide, guide_roi, CP_FORMAT, -1);
	guide_roi->rx = guide_roi->naxes[0] = gui.roi.selection.w;
	guide_roi->ry = guide_roi->naxes[1] = gui.roi.selection.h;
	guide_roi->naxes[2] = nchans;
	guide_roi->naxis = nchans == 1 ? 2 : 3;
	if (guide->type == DATA_FLOAT) {
		guide_roi->fdata = malloc(npixels_roi * nchans * sizeof(float));
		if (!guide_roi->fdata)
			retval = 1;
		guide_roi->fpdata[0] = gui.roi.fit.fdata;
		guide_roi->fpdata[1] = rgb? guide_roi->fdata + npixels_roi : guide_roi->fdata;
		guide_roi->fpdata[2] = rgb? guide_roi->fdata + 2 * npixels_roi : guide_roi->fdata;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				float *srcindex = guide->fdata + (npixels * c) + ((guide->ry - y - gui.roi.selection.y) * guide->rx) + gui.roi.selection.x;
				float *destindex = guide_roi->fdata + (npixels_roi * c) + (guide_roi->rx * y);
				memcpy(destindex, srcindex, (gui.roi.selection.w) * sizeof(float));
			}
		}
	} else {
		guide_roi->data = malloc(npixels_roi * nchans * sizeof(WORD));
		if (!guide_roi->data)
			retval = 1;
		guide_roi->pdata[0] = gui.roi.fit.data;
		guide_roi->pdata[1] = rgb? guide_roi->data + npixels_roi : guide_roi->data;
		guide_roi->pdata[2] = rgb? guide_roi->data + 2 * npixels_roi : guide_roi->data;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				WORD *srcindex = guide->data + (npixels * c) + ((guide->ry - y - gui.roi.selection.y) * guide->rx) + gui.roi.selection.x;
				WORD *destindex = guide_roi->data + (npixels_roi * c) + (guide_roi->rx * y);
				memcpy(destindex, srcindex, (gui.roi.selection.w) * sizeof(WORD));
			}
		}
	}
	return retval;
}

int edge_preserving_filter(fits *fit, fits *guide, double d, double sigma_col, double sigma_space, ep_filter_t filter_type, gboolean verbose) {
	struct timeval t_start, t_end;
	if (sigma_col <= 0.0 || (sigma_space <= 0.0 && filter_type == EP_BILATERAL))
		return 1;
	if (verbose) {
		siril_log_color_message(_("Bilateral filter: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}
	sigma_col = sigma_col / 100.0;

	// cv::BilateralFilter() only works on 8u and 32f images, so we convert 16-bit to 32-bit
	size_t ndata;
	data_type orig_type = fit->type;
	if (orig_type == DATA_USHORT) {
		ndata = fit->rx * fit->ry * fit->naxes[2];
		fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
	}
	double eps = sigma_col * sigma_col;
	fits *guide_roi = NULL, *guidance = NULL;
	gboolean roi_fitting_method;
	switch (filter_type) {
		case EP_BILATERAL:
			cvBilateralFilter(fit, d, eps, sigma_space);
			break;
		case EP_GUIDED:
			guide_roi = malloc(sizeof(fits));
			roi_fitting_needed = (fit == &gui.roi.fit && guide != &gui.roi.fit && gui.roi.active);
			if (roi_fitting_needed)
				match_guide_to_roi(guide, guide_roi);
			guidance = roi_fitting_needed ? guide_roi : guide;
			cvGuidedFilter(fit, guidance, d, eps);
			break;
	}

	if (orig_type == DATA_USHORT) {
		fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
	}

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	return 0;
}

static int bilat_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))))
		copy_backup_to_gfit();
	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;
	if (filter_type == EP_GUIDED) {
		if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")))) {
			guide = fit;
		} else {
			if (loaded_fit) {
				guide = loaded_fit;
			} else {
				return 1;
			}
		}
	} else {
		guide = NULL;
	}
	edge_preserving_filter(fit, guide, bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value, filter_type, FALSE);
	notify_gfit_modified();
	return 0;
}

void bilat_change_between_roi_and_image() {
	// If we are showing the preview, update it after the ROI change.
	roi_supported(TRUE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

static void bilat_startup() {
	copy_gfit_to_backup();
	add_roi_callback(bilat_change_between_roi_and_image);
	roi_supported(TRUE);
}

static void bilat_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Bilateral filter: (d=%2.2f, sigma_col=%2.2f, sigma_spatial=%2.2f)"),
				bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value);
	}
	backup_roi();
	roi_supported(FALSE);
	remove_roi_callback(bilat_change_between_roi_and_image);
	if (loaded_fit) {
		clearfits(loaded_fit);
		free(loaded_fit);
		loaded_fit = NULL;
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int bilat_process_all() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))))
		copy_backup_to_gfit();
	fits *fit = &gfit;
	if (filter_type != EP_BILATERAL) {
		if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")))) {
			guide = fit;
		} else {
			if (loaded_fit) {
				guide = loaded_fit;
			} else {
				siril_log_color_message(_("Error, no guide image loaded."), "red");
				return 1;
			}
		}
	} else {
		guide = NULL;
	}
	edge_preserving_filter(fit, guide, bilat_d_value, bilat_sigma_col_value, bilat_sigma_spatial_value, filter_type, FALSE);
	populate_roi();
	notify_gfit_modified();
	return 0;
}

static void apply_bilat_changes() {
	gboolean status = (bilat_sigma_col_value != 0.0f) || (bilat_sigma_spatial_value != 0.0f);
	bilat_close(!status);
}

void apply_bilat_cancel() {
	bilat_close(TRUE);
	siril_close_dialog("bilat_dialog");
}

/*** callbacks **/

void on_bilat_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_d"));
	GtkSpinButton *spin_bilat_sigma_col = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_col"));
	GtkSpinButton *spin_bilat_sigma_spatial = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_spatial"));

	if (gui.rendering_mode == LINEAR_DISPLAY)
		setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535

	bilat_startup();
	bilat_d_value = 0.0f;
	bilat_sigma_col_value = 0.0f;
	bilat_sigma_spatial_value = 0.0f;
	filter_type = EP_BILATERAL;
	if (loaded_fit) {
		clearfits(loaded_fit);
		free(loaded_fit);
		loaded_fit = NULL;
	}

	set_notify_block(TRUE);
	gtk_widget_set_visible(lookup_widget("guide_image_widgets"), FALSE);
	gtk_widget_set_visible(lookup_widget("bilat_sigma_spatial_settings"), TRUE);
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("guided_filter_guideimage")));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")), TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("ep_filter_type")), EP_BILATERAL);
	gtk_spin_button_set_value(spin_d, bilat_d_value);
	gtk_spin_button_set_value(spin_bilat_sigma_spatial, bilat_sigma_spatial_value);
	gtk_spin_button_set_value(spin_bilat_sigma_col, bilat_sigma_col_value);
	set_notify_block(FALSE);

	/* default parameters don't transform image, no need to update preview */
}

void on_bilat_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_bilat_cancel();
}

void on_bilat_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview"))) || gui.roi.active) {
		bilat_process_all();
	}

	apply_bilat_changes();

	siril_close_dialog("bilat_dialog");
}

void on_bilat_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_bilat_changes();
}

void on_bilat_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_d"));
	GtkSpinButton *spin_bilat_sigma_col = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_col"));
	GtkSpinButton *spin_bilat_sigma_spatial = GTK_SPIN_BUTTON(lookup_widget("spin_bilat_sigma_spatial"));
	bilat_d_value = 0.0;
	bilat_sigma_col_value = 0.0;
	bilat_sigma_spatial_value = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_d, bilat_d_value);
	gtk_spin_button_set_value(spin_bilat_sigma_col, bilat_sigma_col_value);
	gtk_spin_button_set_value(spin_bilat_sigma_spatial, bilat_sigma_spatial_value);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_ep_filter_type_changed(GtkComboBox *combo, gpointer user_data) {
	filter_type = gtk_combo_box_get_active(combo);
	gtk_widget_set_visible(lookup_widget("guide_image_widgets"), filter_type != EP_BILATERAL);
	gtk_widget_set_visible(lookup_widget("bilat_sigma_spatial_settings"), filter_type == EP_BILATERAL);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_guided_filter_selfguide_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);
	if (active) {
		if (loaded_fit) {
			clearfits(loaded_fit);
			free(loaded_fit);
			loaded_fit = NULL;
		}
		gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("guided_filter_guideimage")));
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);

}

void on_guided_filter_guideimage_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	gchar *filename = siril_file_chooser_get_filename(filechooser);
	if (!loaded_fit) {
		loaded_fit = calloc(1, sizeof(fits));
	} else {
		clearfits(loaded_fit);
		free(loaded_fit);
		loaded_fit = NULL;
	}
	if (readfits(filename, loaded_fit, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		return;
	}
	if (loaded_fit->rx != gfit.rx || loaded_fit->ry != gfit.ry) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image dimensions do not match"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		clearfits(loaded_fit);
		free(loaded_fit);
		loaded_fit = NULL;
		return;
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("guided_filter_selfguide")), FALSE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);

}

void on_spin_bilat_d_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_d_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_spin_bilat_sigma_col_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_sigma_col_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_spin_bilat_sigma_spatial_value_changed(GtkSpinButton *button, gpointer user_data) {
	bilat_sigma_spatial_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = bilat_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")));
	notify_update((gpointer) param);
}

void on_bilat_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bilat_preview")))) {
		copy_backup_to_gfit();
		redraw(REMAP_ALL);
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = bilat_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
