/*
* This file is part of Siril, an astronomy image processor.
* Copyright (C) 2005-2023 Free Software Foundation, Inc.
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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/masks.h"
#include "core/processing.h"
#include "core/undo.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/histogram.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "io/image_format_fits.h"
#include "masks_gui.h"
#include "opencv/opencv.h"

/* Static widget pointers */
static GtkWidget *combo_mask_from_image_type = NULL;
static GtkWidget *combo_mask_luminance_type = NULL;
static GtkWidget *combo_mask_from_image_bitdepth = NULL;
static GtkWidget *mask_luminance_grid = NULL;
static GtkWidget *mask_channel_grid = NULL;
static GtkWidget *mask_file_chooser_box = NULL;
static GtkWidget *mask_method_box = NULL;
static GtkWidget *entry_mask_filename = NULL;
static GtkSpinButton *spin_mask_lr = NULL;
static GtkSpinButton *spin_mask_lg = NULL;
static GtkSpinButton *spin_mask_lb = NULL;
static GtkSpinButton *spin_mask_channel = NULL;
static GtkWidget *toggle_mask_autostretch = NULL;
static GtkWidget *toggle_mask_invert = NULL;

/* Static widget pointers for mask from stars dialog */
static GtkWidget *combo_mask_from_stars_bitdepth = NULL;
static GtkSpinButton *spin_mask_star_radius = NULL;
static GtkSpinButton *spin_mask_feather = NULL;
static GtkWidget *toggle_mask_stars_invert = NULL;

/* State for mask from image dialog */
static gboolean mask_from_file_mode = FALSE;
static gchar *selected_mask_filename = NULL;

/* Static widget pointers for color mask dialog */
static GtkWidget *combo_color_mask_bitdepth = NULL;
static GtkWidget *drawing_area_color_display = NULL;
static GtkScale *scale_color_tolerance = NULL;
static GtkScale *scale_lum_min = NULL;
static GtkScale *scale_lum_max = NULL;
static GtkSpinButton *spin_color_feather = NULL;
static GtkWidget *toggle_color_mask_invert = NULL;
static GtkWidget *toggle_color_mask_cleanup = NULL;

/* Store the selected color RGB values */
static float selected_color_r = 1.0f;
static float selected_color_g = 0.39f;
static float selected_color_b = 0.39f;

/* Static widget pointers for mask thresholding */
static GtkWidget *threshold_mask_close = NULL;
static GtkWidget *threshold_mask_apply = NULL;
static GtkSpinButton *spin_mask_threshold_range = NULL;
static GtkSpinButton *spin_mask_threshold_min = NULL;
static GtkSpinButton *spin_mask_threshold_max = NULL;
static GtkAdjustment *adj_mask_thresh_min_float = NULL;
static GtkAdjustment *adj_mask_thresh_max_float = NULL;
static GtkAdjustment *adj_mask_thresh_min_8bit = NULL;
static GtkAdjustment *adj_mask_thresh_max_8bit = NULL;
static GtkAdjustment *adj_mask_thresh_min_16bit = NULL;
static GtkAdjustment *adj_mask_thresh_max_16bit = NULL;
static GtkAdjustment *adj_thresh_feather_float = NULL;
static GtkAdjustment *adj_thresh_feather_8bit = NULL;
static GtkAdjustment *adj_thresh_feather_16bit = NULL;

/**
* on_mask_from_image_dialog_show:
* @widget: The dialog widget
* @user_data: User data passed to the callback
*
* Handler for the dialog show event.
* Initializes static widget pointers on first show.
*/
void on_mask_from_image_dialog_show(GtkWidget *widget, gpointer user_data) {
	static gboolean widgets_initialized = FALSE;

	if (!widgets_initialized) {
		combo_mask_from_image_type = lookup_widget("combo_mask_from_image_type");
		combo_mask_luminance_type = lookup_widget("combo_mask_luminance_type");
		combo_mask_from_image_bitdepth = lookup_widget("combo_mask_from_image_bitdepth");
		mask_luminance_grid = lookup_widget("mask_luminance_grid");
		mask_channel_grid = lookup_widget("mask_channel_grid");
		mask_file_chooser_box = lookup_widget("mask_file_chooser_box");
		mask_method_box = lookup_widget("mask_method_box");
		entry_mask_filename = lookup_widget("entry_mask_filename");
		spin_mask_lr = GTK_SPIN_BUTTON(lookup_widget("spin_mask_lr"));
		spin_mask_lg = GTK_SPIN_BUTTON(lookup_widget("spin_mask_lg"));
		spin_mask_lb = GTK_SPIN_BUTTON(lookup_widget("spin_mask_lb"));
		spin_mask_channel = GTK_SPIN_BUTTON(lookup_widget("spin_mask_channel"));
		toggle_mask_autostretch = lookup_widget("toggle_mask_autostretch");
		toggle_mask_invert = lookup_widget("toggle_mask_invert");

		widgets_initialized = TRUE;
	}

	/* Configure dialog based on mode */
	if (mask_from_file_mode) {
		/* Show file chooser and keep method box visible */
		if (mask_file_chooser_box) gtk_widget_show(mask_file_chooser_box);
		if (mask_method_box) gtk_widget_show(mask_method_box);

		/* Clear previous filename */
		if (selected_mask_filename) {
			g_free(selected_mask_filename);
			selected_mask_filename = NULL;
		}
		if (entry_mask_filename) {
			gtk_entry_set_text(GTK_ENTRY(entry_mask_filename), "");
		}

		/* Set initial visibility based on type combo */
		gint active = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_type));
		if (active == 0) {
			if (mask_luminance_grid) gtk_widget_show(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_hide(mask_channel_grid);
		} else {
			if (mask_luminance_grid) gtk_widget_hide(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_show(mask_channel_grid);
		}
	} else {
		/* Hide file chooser, show method box */
		if (mask_file_chooser_box) gtk_widget_hide(mask_file_chooser_box);
		if (mask_method_box) gtk_widget_show(mask_method_box);

		/* Set initial visibility based on type combo */
		gint active = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_type));
		if (active == 0) {
			if (mask_luminance_grid) gtk_widget_show(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_hide(mask_channel_grid);
		} else {
			if (mask_luminance_grid) gtk_widget_hide(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_show(mask_channel_grid);
		}
	}

	// TODO: initialize the bitdepth combo box based on the preference, when added
}

/**
* mask_from_image_dialog_set_file_mode:
* @file_mode: TRUE to show file chooser, FALSE for normal mode
*
* Sets the dialog mode before showing it.
*/
void mask_from_image_dialog_set_file_mode(gboolean file_mode) {
	mask_from_file_mode = file_mode;
}

/**
* on_mask_from_image_close_clicked:
* @button: The button that was clicked
* @user_data: User data passed to the callback
*
* Handler for the Close button click event.
* Closes the mask from image dialog.
*/
void on_mask_from_image_close_clicked(GtkButton *button, gpointer user_data) {
	/* Clean up filename if in file mode */
	if (mask_from_file_mode && selected_mask_filename) {
		g_free(selected_mask_filename);
		selected_mask_filename = NULL;
	}
	mask_from_file_mode = FALSE;
	siril_close_dialog("mask_from_image_dialog");
}

/**
* on_mask_file_chooser_clicked:
* @button: The button that was clicked
* @user_data: User data passed to the callback
*
* Handler for the file chooser button click event.
* Opens a file chooser dialog to select an image file.
*/
void on_mask_file_chooser_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *dialog;
	GtkFileFilter *filter;
	gint response;

	dialog = gtk_file_chooser_dialog_new("Select Mask Image File",
	                                      GTK_WINDOW(lookup_widget("mask_from_image_dialog")),
	                                      GTK_FILE_CHOOSER_ACTION_OPEN,
	                                      "_Cancel", GTK_RESPONSE_CANCEL,
	                                      "_Open", GTK_RESPONSE_ACCEPT,
	                                      NULL);

	/* Create file filter for supported formats */
	filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, "Image files (*.fits, *.fit, *.fts, *.tif)");
	gtk_file_filter_add_pattern(filter, "*.fits");
	gtk_file_filter_add_pattern(filter, "*.fit");
	gtk_file_filter_add_pattern(filter, "*.fts");
	gtk_file_filter_add_pattern(filter, "*.tif");
	gtk_file_filter_add_pattern(filter, "*.FITS");
	gtk_file_filter_add_pattern(filter, "*.FIT");
	gtk_file_filter_add_pattern(filter, "*.FTS");
	gtk_file_filter_add_pattern(filter, "*.TIF");
	gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

	/* Add all files filter */
	filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, "All files");
	gtk_file_filter_add_pattern(filter, "*");
	gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

	response = gtk_dialog_run(GTK_DIALOG(dialog));

	if (response == GTK_RESPONSE_ACCEPT) {
		gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

		if (filename) {
			/* Store the selected filename */
			if (selected_mask_filename) {
				g_free(selected_mask_filename);
			}
			selected_mask_filename = filename;

			/* Update the entry widget to show the filename */
			if (entry_mask_filename) {
				gtk_entry_set_text(GTK_ENTRY(entry_mask_filename), filename);
			}
		}
	}

	gtk_widget_destroy(dialog);
}

/**
* on_mask_from_image_apply_clicked:
* @button: The button that was clicked
* @user_data: User data passed to the callback
*
* Handler for the Apply button click event.
* Applies the mask creation operation based on the current settings.
*/
void on_mask_from_image_apply_clicked(GtkButton *button, gpointer user_data) {
	/* Get option flags */
	gboolean autostretch = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_mask_autostretch));
	gboolean invert = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_mask_invert));

	/* Get the bitdepth */
	gint bitdepth_index = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_bitdepth));
	uint8_t bitdepth = (uint8_t) bitdepth_index == 0 ? 8 : bitdepth_index == 1 ? 16 : 32;

	/* Variables for mask parameters */
	gdouble weight_r = 0.0, weight_g = 0.0, weight_b = 0.0;
	gint channel = 0;

	/* Get mask type to determine channel or luminance mode */
	gint mask_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_type));

	if (mask_from_file_mode) {
		/* File mode - create mask from external file */
		if (!selected_mask_filename || strlen(selected_mask_filename) == 0) {
			queue_error_message_dialog(_("No file selected"),
			                     _("Please select an image file to create the mask from."));
			return;
		}

		if (mask_type == 0) {
			/* Luminance mode */
			weight_r = gtk_spin_button_get_value(spin_mask_lr);
			weight_g = gtk_spin_button_get_value(spin_mask_lg);
			weight_b = gtk_spin_button_get_value(spin_mask_lb);

			mask_from_lum_data *data = calloc(1, sizeof(mask_from_lum_data));
			data->rw = weight_r;
			data->gw = weight_g;
			data->bw = weight_b;
			data->autostretch = autostretch;
			data->invert = invert;
			data->use_human = FALSE;
			data->use_even = FALSE;
			data->bitpix = bitdepth;
			data->filename = g_strdup(selected_mask_filename);

			struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
			args->fit =  gfit;
			args->mem_ratio = 1.0f;
			args->mask_hook = mask_from_lum_hook;
			args->log_hook = mask_from_lum_log;
			args->description = _("Mask from luminance");
			args->verbose = TRUE;
			args->user = data;
			args->mask_creation = TRUE;
			args->max_threads = com.max_thread;

			start_in_new_thread(generic_mask_worker, args);
		} else {
			/* Channel mode */
			channel = gtk_spin_button_get_value_as_int(spin_mask_channel);

			mask_from_channel_data *data = calloc(1, sizeof(mask_from_channel_data));
			data->channel = channel;
			data->autostretch = autostretch;
			data->invert = invert;
			data->bitpix = bitdepth;
			data->filename = g_strdup(selected_mask_filename);

			struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
			args->fit =  gfit;
			args->mem_ratio = 1.0f;
			args->mask_hook = mask_from_channel_hook;
			args->log_hook = mask_from_channel_log;
			args->description = _("Mask from channel");
			args->verbose = TRUE;
			args->user = data;
			args->mask_creation = TRUE;
			args->max_threads = com.max_thread;

			start_in_new_thread(generic_mask_worker, args);
		}
	} else {
		/* Normal mode - create mask from current image */

		if (mask_type == 0) {
			/* Luminance method */
			weight_r = gtk_spin_button_get_value(spin_mask_lr);
			weight_g = gtk_spin_button_get_value(spin_mask_lg);
			weight_b = gtk_spin_button_get_value(spin_mask_lb);

			mask_from_lum_data *data = calloc(1, sizeof(mask_from_lum_data));
			data->rw = weight_r;
			data->gw = weight_g;
			data->bw = weight_b;
			data->autostretch = autostretch;
			data->invert = invert;
			data->use_human = FALSE;
			data->use_even = FALSE;
			data->bitpix = bitdepth;
			data->filename = NULL;

			struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
			args->fit =  gfit;
			args->mem_ratio = 1.0f;
			args->mask_hook = mask_from_lum_hook;
			args->log_hook = mask_from_lum_log;
			args->description = _("Mask from luminance");
			args->verbose = TRUE;
			args->user = data;
			args->mask_creation = TRUE;
			args->max_threads = com.max_thread;

			start_in_new_thread(generic_mask_worker, args);
		} else if (mask_type == 1) {
			/* Channel method */
			channel = gtk_spin_button_get_value_as_int(spin_mask_channel);

			mask_from_channel_data *data = calloc(1, sizeof(mask_from_channel_data));
			data->channel = channel;
			data->autostretch = autostretch;
			data->invert = invert;
			data->bitpix = bitdepth;
			data->filename = NULL;

			struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
			args->fit =  gfit;
			args->mem_ratio = 1.0f;
			args->mask_hook = mask_from_channel_hook;
			args->log_hook = mask_from_channel_log;
			args->description = _("Mask from channel");
			args->verbose = TRUE;
			args->user = data;
			args->mask_creation = TRUE;
			args->max_threads = com.max_thread;

			start_in_new_thread(generic_mask_worker, args);
		}
	}

	siril_close_dialog("mask_from_image_dialog");
}

/**
* on_combo_mask_from_image_type_changed:
* @combo: The combo box that changed
* @user_data: User data passed to the callback
*
* Handler for mask type combo box change event.
* Shows/hides the appropriate grid based on whether Luminance or Channel is selected.
*/
void on_combo_mask_from_image_type_changed(GtkComboBox *combo, gpointer user_data) {
	gint active = gtk_combo_box_get_active(combo);

	switch (active) {
		case 0: /* Luminance */
			if (mask_luminance_grid) gtk_widget_show(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_hide(mask_channel_grid);
			break;
		case 1: /* Channel */
			if (mask_luminance_grid) gtk_widget_hide(mask_luminance_grid);
			if (mask_channel_grid) gtk_widget_show(mask_channel_grid);
			break;
		default:
			break;
	}
}

/**
* on_combo_mask_luminance_type_changed:
* @combo: The combo box that changed
* @user_data: User data passed to the callback
*
* Handler for luminance type combo box change event.
* Enables/disables the custom luminance weight spin buttons based on selection
* and sets preset values for Even-weighted and Human-weighted options.
*/
void on_combo_mask_luminance_type_changed(GtkComboBox *combo, gpointer user_data) {
	gint active = gtk_combo_box_get_active(combo);

	gboolean sensitive = (active == 2); /* Custom option */

	if (spin_mask_lr) gtk_widget_set_sensitive(GTK_WIDGET(spin_mask_lr), sensitive);
	if (spin_mask_lg) gtk_widget_set_sensitive(GTK_WIDGET(spin_mask_lg), sensitive);
	if (spin_mask_lb) gtk_widget_set_sensitive(GTK_WIDGET(spin_mask_lb), sensitive);

	/* Set preset values based on selection */
	switch (active) {
		case 0: /* Even-weighted */
			if (spin_mask_lr) gtk_spin_button_set_value(spin_mask_lr, 0.333);
			if (spin_mask_lg) gtk_spin_button_set_value(spin_mask_lg, 0.333);
			if (spin_mask_lb) gtk_spin_button_set_value(spin_mask_lb, 0.333);
			break;
		case 1: /* Human-weighted */
			if (spin_mask_lr) gtk_spin_button_set_value(spin_mask_lr, 0.2126);
			if (spin_mask_lg) gtk_spin_button_set_value(spin_mask_lg, 0.7152);
			if (spin_mask_lb) gtk_spin_button_set_value(spin_mask_lb, 0.0722);
			break;
		case 2: /* Custom */
			/* Leave values as-is, user will set them */
			break;
		default:
			break;
	}
}

/**
* on_mask_from_stars_dialog_show:
* @widget: The dialog widget
* @user_data: User data passed to the callback
*
* Handler for the dialog show event.
* Initializes static widget pointers on first show.
*/
void on_mask_from_stars_dialog_show(GtkWidget *widget, gpointer user_data) {
	static gboolean widgets_initialized = FALSE;

	if (!widgets_initialized) {
		combo_mask_from_stars_bitdepth = lookup_widget("combo_mask_from_stars_bitdepth");
		spin_mask_star_radius = GTK_SPIN_BUTTON(lookup_widget("spin_mask_star_radius"));
		spin_mask_feather = GTK_SPIN_BUTTON(lookup_widget("spin_mask_feather"));
		toggle_mask_stars_invert = lookup_widget("toggle_mask_stars_invert");

		widgets_initialized = TRUE;
	}
	// TODO: initialize the bitdepth combo box based on the preference, when added
}

/**
* on_mask_from_stars_close_clicked:
* @button: The button that was clicked
* @user_data: User data passed to the callback
*
* Handler for the Close button click event.
* Closes the mask from stars dialog.
*/
void on_mask_from_stars_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_from_stars_dialog");
}

/**
* on_mask_from_stars_apply_clicked:
* @button: The button that was clicked
* @user_data: User data passed to the callback
*
* Handler for the Apply button click event.
* Applies the mask creation operation from detected stars.
*/
void on_mask_from_stars_apply_clicked(GtkButton *button, gpointer user_data) {
	if (gfit->mask) {
		set_mask_active( gfit, FALSE);
		free_mask(gfit->mask);
		gfit->mask = NULL;
	}

	/* Get star radius (FWHM multiplier) */
	gfloat star_radius = (gfloat) gtk_spin_button_get_value(spin_mask_star_radius);

	/* Get the bitdepth */
	gint bitdepth_index = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_stars_bitdepth));
	uint8_t bitdepth = (uint8_t) bitdepth_index == 0 ? 8 : bitdepth_index == 1 ? 16 : 32;

	/* Get feather radius */
	gdouble feather_radius = gtk_spin_button_get_value(spin_mask_feather);

	/* Get invert option */
	gboolean invert = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_mask_stars_invert));

	undo_save_state( gfit, _("Create mask from stars"));

	mask_from_stars_data *data = calloc(1, sizeof(mask_from_stars_data));
	data->r = star_radius;
	data->feather = feather_radius;
	data->invert = invert;
	data->bitdepth = bitdepth;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_from_stars_hook;
	args->log_hook = mask_from_stars_log;
	args->description = _("Mask from stars");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);

	siril_close_dialog("mask_from_stars_dialog");
}

void on_blur_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_blur_radius"));
	float radius = gtk_spin_button_get_value(spin);

	mask_blur_data *data = calloc(1, sizeof(mask_blur_data));
	data->radius = radius;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 2.0f;
	args->mask_hook = mask_blur_hook;
	args->log_hook = mask_blur_log;
	args->description = _("Blur mask");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void on_feather_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_feather_distance"));
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("combo_mask_feather_type"));
	float distance = gtk_spin_button_get_value(spin);
	feather_mode mode = (feather_mode) gtk_combo_box_get_active(combo);

	mask_feather_data *data = calloc(1, sizeof(mask_feather_data));
	data->distance = distance;
	data->mode = mode;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 2.0f;
	args->mask_hook = mask_feather_hook;
	args->log_hook = mask_feather_log;
	args->description = _("Feather mask");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void on_multiply_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_multiply_factor"));
	float factor = gtk_spin_button_get_value(spin);

	mask_fmul_data *data = calloc(1, sizeof(mask_fmul_data));
	data->factor = factor;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 0.0f;
	args->mask_hook = mask_fmul_hook;
	args->log_hook = mask_fmul_log;
	args->description = _("Multiply mask");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void on_blur_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_blur_dialog");
}

void on_feather_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_feather_dialog");
}

void on_multiply_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_scale_dialog");
}

void on_button_pick_from_image_clicked(GtkButton *button, gpointer user_data) {
	mouse_status = MOUSE_ACTION_SAMPLE_MASK_COLOR;
//	siril_message_dialog(GTK_MESSAGE_INFO, _("Pick a color"),
//	                     _("Click on the image to select a color. The chromaticity will be extracted from the pixel you click."));
}

/**
 * on_color_display_draw:
 * @widget: The drawing area widget
 * @cr: The cairo context
 * @user_data: User data
 *
 * Draw callback for the color display area.
 */
gboolean on_color_display_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
	GtkAllocation allocation;
	gtk_widget_get_allocation(widget, &allocation);

	/* Set the color */
	cairo_set_source_rgb(cr, selected_color_r, selected_color_g, selected_color_b);

	/* Fill the entire area */
	cairo_rectangle(cr, 0, 0, allocation.width, allocation.height);
	cairo_fill(cr);

	/* Draw a border */
	cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	cairo_set_line_width(cr, 1.0);
	cairo_rectangle(cr, 0.5, 0.5, allocation.width - 1, allocation.height - 1);
	cairo_stroke(cr);

	return FALSE;
}

void on_mask_from_color_dialog_show(GtkWidget *widget, gpointer user_data) {
	static gboolean widgets_initialized = FALSE;

	if (!widgets_initialized) {
		combo_color_mask_bitdepth = lookup_widget("combo_color_mask_bitdepth");
		drawing_area_color_display = lookup_widget("drawing_area_color_display");
		scale_color_tolerance = GTK_SCALE(lookup_widget("scale_color_tolerance"));
		scale_lum_min = GTK_SCALE(lookup_widget("scale_lum_min"));
		scale_lum_max = GTK_SCALE(lookup_widget("scale_lum_max"));
		spin_color_feather = GTK_SPIN_BUTTON(lookup_widget("spin_color_feather"));
		toggle_color_mask_invert = lookup_widget("toggle_color_mask_invert");
		toggle_color_mask_cleanup = lookup_widget("toggle_color_mask_cleanup");

		widgets_initialized = TRUE;
	}

	/* Trigger initial draw */
	if (drawing_area_color_display) {
		gtk_widget_queue_draw(drawing_area_color_display);
	}
}

void mask_color_handle_image_click(int x, int y) {
	if (! gfit) return;

	if (x < 0 || x >= gfit->rx || y < 0 || y >= gfit->ry) return;

	size_t pixel_index = (gfit->ry - y) * gfit->rx + x;
	float r, g, b;

	if (gfit->type == DATA_USHORT) {
		r = gfit->pdata[RLAYER][pixel_index] / 65535.0f;
		g = gfit->pdata[GLAYER][pixel_index] / 65535.0f;
		b = gfit->pdata[BLAYER][pixel_index] / 65535.0f;
	} else {
		r = gfit->fpdata[RLAYER][pixel_index];
		g = gfit->fpdata[GLAYER][pixel_index];
		b = gfit->fpdata[BLAYER][pixel_index];
	}

	/* Store the color values */
	selected_color_r = CLAMP(r, 0.0f, 1.0f);
	selected_color_g = CLAMP(g, 0.0f, 1.0f);
	selected_color_b = CLAMP(b, 0.0f, 1.0f);
	siril_debug_print("Selected color R: %f, G: %f, B: %f\n", selected_color_r, selected_color_g, selected_color_b);

	/* Redraw the color display */
	if (drawing_area_color_display) {
		gtk_widget_queue_draw(drawing_area_color_display);
	}
}

void on_mask_color_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_from_color_dialog");
}

void on_mask_color_apply_clicked(GtkButton *button, gpointer user_data) {
	if (! gfit) return;

	/* Check that widgets are initialized */
	if (!combo_color_mask_bitdepth || !drawing_area_color_display ||
	    !scale_color_tolerance || !scale_lum_min || !scale_lum_max ||
	    !spin_color_feather || !toggle_color_mask_invert) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Initialization error"),
		                     _("Color mask dialog widgets not properly initialized."));
		return;
	}

	gint bitdepth_index = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_color_mask_bitdepth));
	uint8_t bitdepth = (uint8_t)(bitdepth_index == 0 ? 8 : bitdepth_index == 1 ? 16 : 32);
	gint feather_radius = gtk_spin_button_get_value_as_int(spin_color_feather);
	gboolean invert = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_color_mask_invert));
	gboolean cleanup = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_color_mask_cleanup));

	float r = selected_color_r;
	float g = selected_color_g;
	float b = selected_color_b;

	/* Ensure we have non-zero values */
	if (r < 0.0001f) r = 0.0001f;
	if (g < 0.0001f) g = 0.0001f;
	if (b < 0.0001f) b = 0.0001f;

	float sum = r + g + b;

	/* This check should now be much more lenient */
	if (sum < 0.0003f) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid color"),
		                     _("Selected color is too close to black. Please choose a color with some brightness."));
		return;
	}

	float chrom_r = r / sum;
	float chrom_g = g / sum;
	float chrom_b = b / sum;

	float tolerance = gtk_range_get_value(GTK_RANGE(scale_color_tolerance));
	float lum_min = gtk_range_get_value(GTK_RANGE(scale_lum_min));
	float lum_max = gtk_range_get_value(GTK_RANGE(scale_lum_max));

	mask_from_color_data *data = calloc(1, sizeof(mask_from_color_data));
	data->chrom_center_r = chrom_r;
	data->chrom_center_g = chrom_g;
	data->chrom_center_b = chrom_b;
	data->chrom_tolerance = tolerance;
	data->lum_min = lum_min;
	data->lum_max = lum_max;
	data->feather_radius = feather_radius;
	data->invert = invert;
	data->bitpix = bitdepth;
	data->cleanup = cleanup;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.5f;
	args->mask_hook = mask_from_color_hook;
	args->log_hook = mask_from_color_log;
	args->description = _("Mask from color");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);

	// For now, don't close the dialog as this mask typically takes a bit of refinement.
//	siril_close_dialog("mask_from_color_dialog");
}

void on_threshold_mask_show(GtkWidget *widget, gpointer user_data) {
	static gboolean widgets_initialized = FALSE;

	if (!widgets_initialized) {
		threshold_mask_close = lookup_widget("threshold_mask_close");
		threshold_mask_apply = lookup_widget("threshold_mask_apply");
		spin_mask_threshold_range = GTK_SPIN_BUTTON(lookup_widget("spin_mask_threshold_range"));
		spin_mask_threshold_min = GTK_SPIN_BUTTON(lookup_widget("spin_mask_threshold_min"));
		spin_mask_threshold_max = GTK_SPIN_BUTTON(lookup_widget("spin_mask_threshold_max"));
		adj_mask_thresh_min_float = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_min_float"));
		adj_mask_thresh_max_float = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_max_float"));
		adj_mask_thresh_min_8bit = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_min_8bit"));
		adj_mask_thresh_max_8bit = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_max_8bit"));
		adj_mask_thresh_min_16bit = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_min_16bit"));
		adj_mask_thresh_max_16bit = GTK_ADJUSTMENT(lookup_gobject("adj_mask_thresh_max_16bit"));
		adj_thresh_feather_float = GTK_ADJUSTMENT(lookup_gobject("adj_thresh_feather_float"));;
		adj_thresh_feather_8bit = GTK_ADJUSTMENT(lookup_gobject("adj_thresh_feather_8bit"));;
		adj_thresh_feather_16bit = GTK_ADJUSTMENT(lookup_gobject("adj_thresh_feather_16bit"));;
		widgets_initialized = TRUE;
	}
	if (gfit->mask) {
		switch(gfit->mask->bitpix) {
			case 8: {
				gtk_spin_button_set_adjustment(spin_mask_threshold_min, adj_mask_thresh_min_8bit);
				gtk_spin_button_set_adjustment(spin_mask_threshold_max, adj_mask_thresh_max_8bit);
				gtk_spin_button_set_adjustment(spin_mask_threshold_range, adj_thresh_feather_8bit);
				gtk_spin_button_set_digits(spin_mask_threshold_min, 0);
				gtk_spin_button_set_digits(spin_mask_threshold_max, 0);
				gtk_spin_button_set_digits(spin_mask_threshold_range, 0);
				break;
			}
			case 16: {
				gtk_spin_button_set_adjustment(spin_mask_threshold_min, adj_mask_thresh_min_16bit);
				gtk_spin_button_set_adjustment(spin_mask_threshold_max, adj_mask_thresh_max_16bit);
				gtk_spin_button_set_adjustment(spin_mask_threshold_range, adj_thresh_feather_16bit);
				gtk_spin_button_set_digits(spin_mask_threshold_min, 0);
				gtk_spin_button_set_digits(spin_mask_threshold_max, 0);
				gtk_spin_button_set_digits(spin_mask_threshold_range, 0);
				break;
			}
			default: {
				gtk_spin_button_set_adjustment(spin_mask_threshold_min, adj_mask_thresh_min_float);
				gtk_spin_button_set_adjustment(spin_mask_threshold_max, adj_mask_thresh_max_float);
				gtk_spin_button_set_adjustment(spin_mask_threshold_range, adj_thresh_feather_float);
				gtk_spin_button_set_digits(spin_mask_threshold_min, 3);
				gtk_spin_button_set_digits(spin_mask_threshold_max, 3);
				gtk_spin_button_set_digits(spin_mask_threshold_range, 3);
				break;
			}
		}
	}
}

void on_threshold_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_thresholds_dialog");
}

void on_threshold_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	float min_val = gtk_spin_button_get_value(spin_mask_threshold_min);
	float max_val = gtk_spin_button_get_value(spin_mask_threshold_max);
	float range = gtk_spin_button_get_value(spin_mask_threshold_range);

	mask_thresh_data *data = calloc(1, sizeof(mask_thresh_data));
	data->min_val = min_val;
	data->max_val = max_val;
	data->range = range;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_thresh_hook;
	args->log_hook = mask_thresh_log;
	args->description = _("Mask intensity thresholding");
	args->verbose = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}
