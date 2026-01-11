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
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/histogram.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "io/image_format_fits.h"
#include "masks_gui.h"

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
	int result = 0;

	if (mask_from_file_mode) {
		/* File mode - create mask from external file */
		if (!selected_mask_filename || strlen(selected_mask_filename) == 0) {
			queue_error_message_dialog(_("No file selected"),
			                     _("Please select an image file to create the mask from."));
			return;
		}

		/* Get mask type to determine channel or luminance mode */
		gint mask_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_type));

		if (mask_type == 0) {
			/* Luminance mode: channel = -1, pass weights */
			weight_r = gtk_spin_button_get_value(spin_mask_lr);
			weight_g = gtk_spin_button_get_value(spin_mask_lg);
			weight_b = gtk_spin_button_get_value(spin_mask_lb);
			result = mask_create_from_image(gfit, selected_mask_filename, -1, bitdepth, weight_r, weight_g, weight_b);
		} else {
			/* Channel mode */
			channel = gtk_spin_button_get_value_as_int(spin_mask_channel);
			result = mask_create_from_image(gfit, selected_mask_filename, channel, bitdepth, 0.0, 0.0, 0.0);
		}

		if (result != 0) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Mask creation failed"),
			                     _("Failed to create mask from the selected file."));
			return;
		}
	} else {
		/* Normal mode - create mask from current image */
		gint mask_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_mask_from_image_type));

		// We have to autostretch first, otherwise we get signal crushing
		fits *fit = gfit;
		if (autostretch) {
			fit = calloc(1, sizeof(fits));
			copyfits(gfit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
			struct mtf_data *data = create_mtf_data();
			if (!data) {
				PRINT_ALLOC_ERR;
				clearfits(fit);
				free(fit);
				return;
			}
			data->fit = fit;
			data->auto_display_compensation = FALSE;
			data->is_preview = FALSE;
			data->linked = TRUE;

			// Create generic_img_args
			struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				PRINT_ALLOC_ERR;
				destroy_mtf_data(data);
				clearfits(fit);
				free(fit);
				return;
			}
			// Compute the autostretch parameters
			find_linked_midtones_balance(fit, AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, &data->params);
			data->params.do_red = data->params.do_green = data->params.do_blue = TRUE;

			args->fit = fit;
			args->mem_ratio = 1.0f;
			args->image_hook = mtf_single_image_hook;
			args->description = _("Autostretch mask");
			args->idle_function = end_generic_image;
			args->user = data;
			args->max_threads = com.max_thread;

			// Run worker synchronously - cleanup happens via destructor
			gpointer result = generic_image_worker(args);

			if (result) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Mask creation failed"),
									_("Failed to create mask."));
				destroy_mtf_data(data);
				clearfits(fit);
				free(fit);
				return;
			}
		}
		if (mask_type == 0) {
			/* Luminance method */
			weight_r = gtk_spin_button_get_value(spin_mask_lr);
			weight_g = gtk_spin_button_get_value(spin_mask_lg);
			weight_b = gtk_spin_button_get_value(spin_mask_lb);

			mask_create_from_luminance(gfit, fit, weight_r, weight_g, weight_b, bitdepth);
		} else if (mask_type == 1) {
			/* Channel method */
			channel = gtk_spin_button_get_value_as_int(spin_mask_channel);
			mask_create_from_channel(gfit, fit, channel, bitdepth);
		}
		if (fit != gfit) { // autostretch - we have to clear up our temporary fits!
			clearfits(fit);
			free(fit);
		}
	}

	if (invert) {
		mask_invert(gfit);
	}
	queue_redraw_mask();
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
	if (gfit && gfit->mask) {
		set_mask_active(gfit, FALSE);
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

	/* Create mask from stars */
	int result = mask_create_from_stars(gfit, star_radius, bitdepth);

	if (result != 0) {
		/* Error occurred during mask creation */
		return;
	}

	/* Apply feathering if radius > 0 */
	if (feather_radius > 0.0) {
		mask_feather(gfit, feather_radius, FEATHER_OUTER);
	}

	/* Apply inversion if requested */
	if (invert) {
		mask_invert(gfit);
	}

	queue_redraw_mask();
	siril_close_dialog("mask_from_image_dialog");
}

void on_blur_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_blur_radius"));
	float radius = gtk_spin_button_get_value(spin);
	mask_apply_gaussian_blur(gfit, radius);
	redraw_mask_idle(NULL);
}

void on_feather_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_feather_distance"));
	GtkComboBox *combo = GTK_COMBO_BOX(lookup_widget("combo_mask_feather_type"));
	float distance = gtk_spin_button_get_value(spin);
	feather_mode mode = (feather_mode) gtk_combo_box_get_active(combo);
	mask_feather(gfit, distance, mode);
	redraw_mask_idle(NULL);
}

void on_multiply_mask_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin = GTK_SPIN_BUTTON(lookup_widget("spin_mask_multiply_factor"));
	float factor = gtk_spin_button_get_value(spin);
	mask_scale(gfit, factor);
	redraw_mask_idle(NULL);
}

void on_blur_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_blur_dialog");
}

void on_feather_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_feather_dialog");
}

void on_multiply_mask_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mask_fmul_dialog");
}
