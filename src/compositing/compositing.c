/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/command.h" // process_close
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/astrometry_solver.h"
#include "algos/colors.h"
#include "algos/fitting.h"
#include "algos/siril_wcs.h"
#include "filters/linear_match.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/PSF_list.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/photometric_cc.h"
#include "gui/progress_and_log.h"
#include "gui/sequence_list.h"
#include "gui/colors.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "opencv/opencv.h"

#include "compositing.h"
#include "filters.h"

#undef DEBUG

static int compositing_loaded = 0;

typedef enum {
	HSL,
	HSV,
	CIELAB
} coloring_type_enum;

static coloring_type_enum coloring_type = HSL;

/* The result is stored in gfit.
 * gfit.rx and gfit.ry are the reference 1x1 binning output image size. */

// Absolute maximum possible layers
#define MAX_LAYERS 8
// Working limit (dependent on memory check)
static int maximum_layers = MAX_LAYERS;

/* the list of layers. It is dynamic in content but fixed in size.
 * The first element is reserved for luminance and cannot be removed. */
static layer *layers[MAX_LAYERS+1];	// NULL terminated
static int layers_count = 1;		// the glade has only luminance
static unsigned int orig_rx[MAX_LAYERS] = { 0 }; // these are used to hold the original layer image
static unsigned int orig_ry[MAX_LAYERS] = { 0 }; // dimensions before any upscaling

static int luminance_mode = 0;		// 0 if luminance is not used

static gboolean timespan_warning_given = FALSE;

static struct registration_method *reg_methods[5];

static sequence *seq = NULL;		// the sequence of layers, for alignments and normalization
static norm_coeff *coeff = NULL;	// the normalization coefficients
static transformation_type the_type = HOMOGRAPHY_TRANSFORMATION; // always HOMOGRAPHY_TRANSFORMATION unless using the shift-only method

static GdkRGBA list_of_12_palette_colors[12];
static const char *list_of_12_color_names[12] = {
	"#ff0000", "#7f0000",
	"#00ff00", "#007f00",
	"#0000ff", "#00007f",
	"#ffff00", "#7f7f00",
	"#ff00ff", "#7f007f",
	"#00ffff", "#007f7f"
};

static GtkColorChooserDialog *color_dialog = NULL;
static int current_layer_color_choosing = 0;
static int color_quick_edit = 0;
static GdkRGBA qe_ref_color;
static GtkEntry *wl_entry = NULL;
static GtkComboBoxText *box = NULL;

static GtkGrid *grid_layers = NULL;

static GtkButton *add_button = NULL;
static cmsHPROFILE reference = NULL;
/******* internal functions *******/
static void remove_layer(int layer);
static void add_the_layer_add_button();
static void grid_add_row(int layer, int index, int first_time);
static void grid_remove_row(int layer, int free_the_row);
static int has_fit(int layer);
static int number_of_images_loaded();
static void update_compositing_registration_interface();
static float get_composition_pixel_value(int fits_index, int reg_layer, int x, int y);
static void increment_pixel_components_from_layer_value(int fits_index, GdkRGBA *rgbpixel, float vlayer_pixel_value);
static void increment_pixel_components_from_layer_saturated_value(int fits_index, GdkRGBA *rgbpixel, float layer_pixel_value);
static void colors_align_and_compose();		// the rgb procedure
static void luminance_and_colors_align_and_compose();	// the lrgb procedure
static void color_has_been_updated(int layer);
static void update_color_from_saturation(int layer, double newl);
static void rgb_pixel_limiter(GdkRGBA *pixel);
static void clear_pixel(GdkRGBA *pixel);
static void update_result(int assume_gfit_is_ours);
static void populate_filter_lists();
static void coeff_clear();
int register_manual(struct registration_args *regargs);

/* callbacks for programatic GTK */
void on_layer_remove(const GtkButton *button, gpointer user_data);
gboolean on_color_button_press_event(const GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data);
gboolean on_color_button_release_event(const GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data);
gboolean on_color_button_motion_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data);
gboolean draw_layer_color(GtkDrawingArea *widget, cairo_t *cr, gpointer data);
void on_filechooser_file_set(GtkFileChooserButton *widget, gpointer user_data);
void on_centerbutton_toggled(GtkToggleButton *button, gpointer user_data);


/********************************************************/

/* mem_limits function */
/* A bit like the sequence compute_mem_limits functions except this one is used
 * to set a maximum on the number of layers that can be added, based on the size
 * of the first one loaded. */
static void compute_compositor_mem_limits(fits* fit) {
	unsigned int MB_per_image, MB_per_scaled_image, MB_avail, required;
	float overlap_allowance = 1.5f;
	int limit = compute_nb_images_fit_memory_from_fit(fit, 1.0, FALSE, &MB_per_image, &MB_per_scaled_image, &MB_avail);
	if (limit > 0) {
		uint64_t float_channel_size = fit->rx * fit->ry * sizeof(float) / BYTES_IN_A_MB;
		if (float_channel_size == 0)
			float_channel_size = 1;
		required = float_channel_size * 2;
		limit = (MB_avail - (3.f * float_channel_size) * overlap_allowance) / required;
		// UI limitations make it practically difficult to deal with > 10 images
		// regardless of memory constraints, so we retain 10 as an upper bound.
		if (limit > 8) {
			limit = 8;
		}
		if (limit < 3) {
			siril_log_color_message(_("Warning: memory limit check determines that fewer than 3 images of this size can be composited. RGB composition is not possible without freeing up more memory.\n"), "salmon");
		} else if (limit < 4) {
			siril_log_color_message(_("Warning: memory limit check determines that fewer than 4 images of this size can be composited. LRGB composition is not possible without freeing up more memory.\n"), "salmon");
		} else {
			siril_log_message(_("Based on the available memory and initial image dimensions, up to %d images may be composited (the maximum can never exceed 8).\n"), limit);
		}
	}
	maximum_layers = limit;
	// Ensure layers exceeding the limit are removed from the GUI
	for (int layer = MAX_LAYERS ; layer > limit ; layer--) {
		remove_layer(layer);
	}
	int count = 0;
	while (layers[count])
		count++;
	if (count == maximum_layers)
		gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));

}

/* the compositing menu callback */

/* creates a new row with all widgets and bindings at the row index in the
 * layers grid. Indices start at 0, but row 0 holds only one label, and row 1 is
 * reserved to the luminance layer. */
layer *create_layer(int index) {
	layer *ret = malloc(sizeof(layer));
	/* create the widgets and set properties and signals */
	ret->remove_button = GTK_BUTTON(gtk_button_new());
	gtk_button_set_image(ret->remove_button,
			gtk_image_new_from_icon_name("list-remove", GTK_ICON_SIZE_BUTTON));
	g_object_ref(G_OBJECT(ret->remove_button));	// don't destroy it on removal from grid
	g_signal_connect(ret->remove_button, "clicked", G_CALLBACK(on_layer_remove), NULL);

	ret->color_w = GTK_DRAWING_AREA(gtk_drawing_area_new());
	gtk_widget_set_events(GTK_WIDGET(ret->color_w),
			GDK_BUTTON_MOTION_MASK | GDK_BUTTON_PRESS_MASK |
			GDK_BUTTON_RELEASE_MASK);
	g_signal_connect(GTK_WIDGET(ret->color_w), "button-release-event", G_CALLBACK(on_color_button_release_event), NULL);
	g_signal_connect(GTK_WIDGET(ret->color_w), "button-press-event", G_CALLBACK(on_color_button_press_event), NULL);
	g_signal_connect(GTK_WIDGET(ret->color_w), "motion-notify-event", G_CALLBACK(on_color_button_motion_event), NULL);
	g_signal_connect(ret->color_w, "draw", G_CALLBACK(draw_layer_color), NULL);
	g_object_ref(G_OBJECT(ret->color_w));	// don't destroy it on removal from grid

	ret->chooser = GTK_FILE_CHOOSER_BUTTON(gtk_file_chooser_button_new("Select source image", GTK_FILE_CHOOSER_ACTION_OPEN));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(ret->chooser), com.wd);
	siril_set_file_filter(GTK_FILE_CHOOSER(ret->chooser), "filefilter_comp", "FITS and TIFF files");

	gtk_file_chooser_button_set_width_chars(ret->chooser, 16);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(ret->chooser), FALSE);
	g_signal_connect(ret->chooser, "file-set", G_CALLBACK(on_filechooser_file_set), NULL);
	g_object_ref(G_OBJECT(ret->chooser));	// don't destroy it on removal from grid

	ret->label = GTK_LABEL(gtk_label_new(_("not loaded")));
	gtk_widget_set_tooltip_text(GTK_WIDGET(ret->label), _("not loaded"));
	g_object_ref(G_OBJECT(ret->label));	// don't destroy it on removal from grid

	ret->spinbutton_x = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(-10000.0, 10000.0, 0.1));
	gtk_spin_button_set_value(ret->spinbutton_x, 0.0);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->spinbutton_x), FALSE);
	g_object_ref(G_OBJECT(ret->spinbutton_x));	// don't destroy it on removal from grid

	ret->spinbutton_y = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(-10000.0, 10000.0, 0.1));
	gtk_spin_button_set_value(ret->spinbutton_y, 0.0);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->spinbutton_y), FALSE);
	g_object_ref(G_OBJECT(ret->spinbutton_y));	// don't destroy it on removal from grid

	ret->spinbutton_r = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(-360.0, 360.0, 0.001));
	gtk_spin_button_set_value(ret->spinbutton_r, 0.0);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->spinbutton_r), FALSE);
	g_object_ref(G_OBJECT(ret->spinbutton_r));	// don't destroy it on removal from grid

	ret->centerbutton = GTK_TOGGLE_BUTTON(gtk_toggle_button_new());
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ret->centerbutton), FALSE);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->centerbutton), FALSE);
	gtk_widget_set_tooltip_text(GTK_WIDGET(ret->centerbutton), _("Set rotation center for this layer"));
	gtk_button_set_image(GTK_BUTTON(ret->centerbutton), gtk_image_new_from_icon_name("gtk-cdrom", GTK_ICON_SIZE_BUTTON));
	g_signal_connect(ret->centerbutton, "toggled", G_CALLBACK(on_centerbutton_toggled), NULL);
	g_object_ref(G_OBJECT(ret->centerbutton));	// don't destroy it on removal from grid

	/* set other layer data */
	if (number_of_images_loaded() > 0) {
		// If we have images loaded, set the layer center to the center of the image
		for (int layer = 0 ; layer < maximum_layers ; layer++) {
			if (layers[layer]->the_fit.rx != 0) {
				ret->center.x = layers[layer]->the_fit.rx / 2.0;
				ret->center.y = layers[layer]->the_fit.ry / 2.0;
				break;
			}
		}
	} else {
		// Otherwise initialize it to zero
		memset(&ret->center, 0, sizeof(point));
	}
	memset(&ret->the_fit, 0, sizeof(fits));
	assert(index >= 2);	// 1 is luminance
	if (index <= 7) {		// copy default RGB colours
		memcpy(&ret->color,
				&list_of_12_palette_colors[(index-2)*2],
				sizeof(GdkRGBA));
		memcpy(&ret->display_color, &ret->color, sizeof(GdkRGBA));
	} else {
		clear_pixel(&ret->color);
	}
	clear_pixel(&ret->saturated_color);
	return ret;
}

/* callback of the '+' button that is clicked to add a layer in the list */
void on_layer_add(GtkButton *button, gpointer user_data) {
	++layers_count;

	/* move down the plus button */
	gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));
	if (layers_count < maximum_layers)
		add_the_layer_add_button();

	/* add the new layer */
	layers[layers_count-1] = create_layer(layers_count);
	layers[layers_count] = NULL;

	grid_add_row(layers_count-1, layers_count, 1);

	color_has_been_updated(layers_count-1);

	coeff_clear();
}

/* adds the '+' button at the bottom of the list. Creates the trailing grid row too */
static void add_the_layer_add_button() {
	int first_time = 0;
	if (!add_button) {
		first_time = 1;
		add_button = GTK_BUTTON(gtk_button_new());
		gtk_button_set_image(add_button,
				gtk_image_new_from_icon_name("list-add", GTK_ICON_SIZE_BUTTON));
		g_object_ref(G_OBJECT(add_button));	// don't destroy it on removal from grid
		g_signal_connect(add_button, "clicked", G_CALLBACK(on_layer_add), NULL);
	}

	gtk_grid_attach(grid_layers, GTK_WIDGET(add_button), 0, layers_count+1, 1, 1);
	if (first_time)
		gtk_widget_show(GTK_WIDGET(add_button));
}

static void remove_layer(int layer) {
	// Abort if the layer we are asked to remove does not exist
	if (!layers[layer])
		return;
	int refresh = 0;
	if (layer != maximum_layers-1) {
		// the add button is not present if we're at the maximum number of layers
		gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));
	}
	if (has_fit(layer)) {
		clearfits(&layers[layer]->the_fit);
		refresh = 1;
	}
	grid_remove_row(layer, 1);	// remove these widgets for good
	free(layers[layer]);

	do {
		layers[layer] = layers[layer+1];
		grid_remove_row(layer, 0);		// switch rows
		grid_add_row(layer, layer+1, 0);
		layer++;
	} while (layers[layer]) ;

	--layers_count;

	add_the_layer_add_button();

	coeff_clear();

	if (refresh)
		update_result(0);
}

/* callback of the '-' button that is clicked to remove a layer in the list */
void on_layer_remove(const GtkButton *button, gpointer user_data) {
	int layer;
	for (layer = 1; layers[layer]; layer++)
		if (layers[layer]->remove_button == button)
			break;
	if (!layers[layer]) return;

	remove_layer(layer);
}

static void grid_remove_row(int layer, int free_the_row) {
	GtkContainer *cont = GTK_CONTAINER(grid_layers);
	if (!layers[layer]) return;
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->remove_button));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->color_w));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->chooser));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->label));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->spinbutton_x));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->spinbutton_y));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->spinbutton_r));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->centerbutton));
	if (free_the_row) {
		g_object_unref(G_OBJECT(layers[layer]->remove_button));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->color_w));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->chooser));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->label));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->spinbutton_x));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->spinbutton_y));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->spinbutton_r));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->centerbutton));	// don't destroy it on removal from grid
	}
}

static void grid_add_row(int layer, int index, int first_time) {
	if (!layers[layer]) return;
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->remove_button),	0, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->color_w),	1, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->chooser),	2, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->label),		3, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->spinbutton_x),	4, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->spinbutton_y),	5, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->spinbutton_r),	6, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->centerbutton),	7, index, 1, 1);

	if (first_time) {
		gtk_widget_show(GTK_WIDGET(layers[layer]->remove_button));
		gtk_widget_show(GTK_WIDGET(layers[layer]->color_w));
		gtk_widget_show(GTK_WIDGET(layers[layer]->chooser));
		gtk_widget_show(GTK_WIDGET(layers[layer]->label));
		gtk_widget_show(GTK_WIDGET(layers[layer]->spinbutton_x));
		gtk_widget_show(GTK_WIDGET(layers[layer]->spinbutton_y));
		gtk_widget_show(GTK_WIDGET(layers[layer]->spinbutton_r));
		gtk_widget_show(GTK_WIDGET(layers[layer]->centerbutton));
	}
}

/* load all glade data, connect signals, configure the dynamic objects of the
 * composition window and make it visible */
void open_compositing_window() {
	int i;
	GtkWidget *button = lookup_widget("demosaicingButton");
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
		siril_log_color_message(_("Disabling debayer-on-open setting: this must be unset in order to open monochrome images in the compositing tool.\n"), "salmon");
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), FALSE);
	}

	if (!compositing_loaded) {
		register_selection_update_callback(update_compositing_registration_interface);

		gtk_builder_connect_signals(gui.builder, NULL);
		siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("filechooser_lum")), "filefilter_comp", "FITS and TIFF files");

		/* parse the default palette */
		for (i=0; i<sizeof(list_of_12_color_names)/sizeof(const char*); i++)
			gdk_rgba_parse(&list_of_12_palette_colors[i], list_of_12_color_names[i]);
		color_dialog = GTK_COLOR_CHOOSER_DIALOG(gtk_builder_get_object(gui.builder, "colorchooserdialog"));
		gtk_color_chooser_add_palette(GTK_COLOR_CHOOSER(color_dialog),
				GTK_ORIENTATION_VERTICAL, 2, 12, list_of_12_palette_colors);
		populate_filter_lists();
		wl_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_wavelength"));
		box = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboboxtext_filters"));


		/* allocate default layers and populate widget data */
		grid_layers = GTK_GRID(gtk_builder_get_object(gui.builder, "grid_layers"));
		add_the_layer_add_button();

		layers[0] = calloc(1, sizeof(layer));
		layers[0]->chooser = GTK_FILE_CHOOSER_BUTTON(gtk_builder_get_object(gui.builder, "filechooser_lum"));
		gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(layers[0]->chooser), com.wd);
		layers[0]->label = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_lum"));
		layers[0]->spinbutton_x = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbutton_lum_x"));
		layers[0]->spinbutton_y = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbutton_lum_y"));
		layers[0]->spinbutton_r = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbutton_lum_r"));
		layers[0]->centerbutton = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "centerbutton_lum"));

		for (i=1; i<4; i++) {
			/* Create the three default layers */
			on_layer_add(NULL, NULL);
		}
		layers[i] = NULL;

		/* the list below depends on the content of the glade file. It
		 * should be done in the same way as in registration.c, but it
		 * woud be easier if the two glades are merged. */
		reg_methods[0] = new_reg_method(_("Deep Sky (two-step global star registration)"), &register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
		reg_methods[1] = new_reg_method(_("Planetary (DFT image pattern alignment)"), &register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
		reg_methods[2] = new_reg_method(_("Planetary (KOMBAT image pattern alignment)"), &register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_PLANETARY);
		reg_methods[3] = new_reg_method(_("Manual alignment"), &register_manual, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);

		reg_methods[4] = NULL;
		update_compositing_registration_interface();
		/* fill compositing_align_method_combo */
		GtkComboBoxText *aligncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "compositing_align_method_combo"));
		gtk_combo_box_text_remove_all(aligncombo);
		i = 0;
		while (reg_methods[i] != NULL) {
			gtk_combo_box_text_append_text(aligncombo, reg_methods[i]->name);
			i++;
		}
		if (i > 0) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(aligncombo), 0);
		}
		update_compositing_registration_interface();
		compositing_loaded = 1;
	} else {
		/* not the first load, update the CWD just in case it changed in the meantime */
		i = 0;

		do {
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(layers[i]->chooser), com.wd);
			if (!reference)
				reference = copyICCProfile(layers[i]->the_fit.icc_profile);
			i++;
		} while (layers[i]);
		gui_function(update_MenuItem, NULL);
	}
	if (compositing_loaded == 1)
		siril_open_dialog("composition_dialog");
}

/* returns true if the layer number layer has a loaded FITS image */
static int has_fit(int layer) {
	return (layers[layer] && layers[layer]->the_fit.rx != 0);
}

/* returns the number of images loaded */
static int number_of_images_loaded() {
	int i, count = 0;
	for (i=0; layers[i]; i++)
		if (has_fit(i))
			count++;
	return count;
}

/* returns true if none of the colour layers have an image loaded */
static int no_color_available() {	// don't test luminance
	int i;
	for (i=1; layers[i]; i++) {
		if (has_fit(i))
			return 0;
	}
	return 1;
}

/* the 'enable luminance' checkbox callback */
void on_composition_use_lum_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	luminance_mode = gtk_toggle_button_get_active(togglebutton);

	// Reset the luminance file if it is toggled off when the maximum number
	// of files are loaded.
	if (!luminance_mode && layers[0]->the_fit.rx != 0 && number_of_images_loaded() == maximum_layers) {
		clearfits(&layers[0]->the_fit);
		gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(layers[0]->chooser));
		gtk_label_set_text((GtkLabel*)layers[0]->label, _("not loaded"));
	}

	// Update the result if necessary
	if (has_fit(0) && number_of_images_loaded() >= 1)
		update_result(0);
}

static void check_gfit_is_ours() {
	if (number_of_images_loaded() < 1)
		return;
	gboolean update_needed = FALSE;
	int update_from_layer = -1;
	if (luminance_mode) {
		if (has_fit(0) &&
				(layers[0]->the_fit.rx != gfit.rx ||
				 layers[0]->the_fit.ry != gfit.ry ||
				 !gfit.fdata)) {
			update_needed = TRUE;
			update_from_layer = 0;
		}
	} else {
		for (int i = 1; layers[i]; i++)
			if (has_fit(i) &&
					(layers[i]->the_fit.rx != gfit.rx ||
					 layers[i]->the_fit.ry != gfit.ry ||
					 !gfit.fdata)) {
				update_needed = TRUE;
				update_from_layer = i;
				break;
			}
	}

	if (!update_needed)
		return;
	/* create the new result image if it's the first opened image */
	close_single_image();
	if (copyfits(&layers[update_from_layer]->the_fit, &gfit, CP_ALLOC | CP_FORMAT | CP_INIT | CP_EXPAND, -1)) {
		clearfits(&layers[update_from_layer]->the_fit);
		siril_log_color_message(_("Could not display image, unloading it\n"), "red");
		return;
	}
	icc_auto_assign(&gfit, ICC_ASSIGN_ON_COMPOSITION);
	/* open the single image.
	 * code taken from stacking.c:start_stacking() and read_single_image() */
	clear_stars_list(TRUE);
	com.seq.current = UNRELATED_IMAGE;
	char *temp = strdup(_("Unsaved compositing result"));
	if (!create_uniq_from_gfit(temp, FALSE))
		com.uniq->comment = strdup(_("Compositing result image"));
	else
		free(temp);
	initialize_display_mode();
	update_zoom_label();
	display_filename();
	gui_function(set_precision_switch, NULL);
	sliders_mode_set_state(gui.sliders);

	init_layers_hi_and_lo_values(MIPSLOHI);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();
	set_display_mode();
	gui_function(update_MenuItem, NULL);
	redraw(REMAP_ALL);

	sequence_list_change_current();
}

// Called from the filechooser
static void update_metadata(gboolean do_sum) {
	int nb = number_of_images_loaded();
	fits **f = malloc((nb + 1) * sizeof(fits *));
	int j = 0;
	int firstlayer = -1;
	for (int i = 0; layers[i] ; i++)
		if (has_fit(i)) {
			if (firstlayer == -1)
				firstlayer = i;
			f[j++] = &layers[i]->the_fit;
		}
	f[j] = NULL;

	merge_fits_headers_to_result2(&gfit, f, do_sum);
	update_fits_header(&gfit);
	free(f);
}

// Called after alignment
static void update_comp_metadata(fits *fit, gboolean do_sum) {
	int nb = number_of_images_loaded();
	fits **f = malloc((nb + 1) * sizeof(fits *));
	int j = 0;
	for (int i = 0; i < nb ; i++)
		if (seq->internal_fits[i])
			f[j++] = seq->internal_fits[i];
	f[j] = NULL;

	merge_fits_headers_to_result2(&gfit, f, do_sum);
	update_fits_header(&gfit);
	free(f);
}

/* callback for the file chooser's file selection: try to load the pointed file, allocate the
 * destination image if this is the first image, and update the result. */
void on_filechooser_file_set(GtkFileChooserButton *chooser, gpointer user_data) {
	int layer, retval;
	char buf[48], *filename;

	if (!number_of_images_loaded()
			&& (single_image_is_loaded() || (sequence_is_loaded()))) {
		process_close(0);
	}
	check_gfit_is_ours();

	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->chooser == chooser)
			break;
	if (!layers[layer]) return;	// not found
	/* layer is the number of the opened layer */

	// If this layer doesn't already have a fit loaded, and
	// we already have the maximum number of images loaded,
	// error message and return. This may happen because the
	// memory check function can't be too strict about removing
	// layers, it has to allow for either luminance or non-
	// luminance compositions
	if (layers[layer]->the_fit.rx == 0 && number_of_images_loaded() == maximum_layers) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error: image could not be loaded"), _("The maximum number of images of this size has been reached based on available memory limits."));
		return;
	}

	filename = siril_file_chooser_get_filename(GTK_FILE_CHOOSER(chooser));
	if (!filename) return;
	if (layers[layer]->the_fit.rx != 0) {	// already loaded image
		clearfits(&layers[layer]->the_fit);
	}
	if ((retval = read_single_image(filename, &layers[layer]->the_fit,
					NULL, FALSE, NULL, FALSE, TRUE))) {
		gtk_label_set_markup(layers[layer]->label, _("<span foreground=\"red\">ERROR</span>"));
		gtk_widget_set_tooltip_text(GTK_WIDGET(layers[layer]->label), _("Cannot load the file, See the log for more information."));
	} else {
		/* first we want test that we load a single-channel image */
		if (layers[layer]->the_fit.naxes[2] > 1) {
			gtk_label_set_markup(layers[layer]->label, _("<span foreground=\"red\">ERROR</span>"));
			gtk_widget_set_tooltip_text(GTK_WIDGET(layers[layer]->label), _("Only single channel images can be loaded"));
			retval = 1;
		} else {
			/* Force first tab to be Red and not B&W if an image was already loaded */
			GtkNotebook* Color_Layers = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook1"));
			GtkWidget *page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
			gtk_notebook_set_tab_label_text(GTK_NOTEBOOK(Color_Layers), page, _("Red"));
			gui_function(close_tab, NULL);

			if (number_of_images_loaded() == 1) { // set orig_rx and orig_ry, and check memory limits
				orig_rx[layer] = layers[layer]->the_fit.rx;
				orig_ry[layer] = layers[layer]->the_fit.ry;
				compute_compositor_mem_limits(&layers[layer]->the_fit);
				if (reference) {
					cmsCloseProfile(reference);
					reference = NULL;
				}
				if (layers[layer]->the_fit.icc_profile)
					reference = copyICCProfile(layers[layer]->the_fit.icc_profile);
			}
			if (number_of_images_loaded() > 1 && !profiles_identical(reference,
							layers[layer]->the_fit.icc_profile)) {
				if (reference) {
					siril_log_color_message(_("ICC profile differs to that of the first image loaded. "
								"Converting this image to match the first one loaded.\n"), "salmon");
					siril_colorspace_transform(&layers[layer]->the_fit, reference);
				} else {
					siril_log_color_message(_("Input images have inconsistent ICC profiles. First image "
								"had no ICC profile. All input layers will be treated as raw data.\n"), "salmon");
					cmsCloseProfile(layers[layer]->the_fit.icc_profile);
					layers[layer]->the_fit.icc_profile = NULL;
					color_manage(&layers[layer]->the_fit, FALSE);
				}
			}
			if (number_of_images_loaded() > 1 &&
					(gfit.rx != layers[layer]->the_fit.rx ||
							gfit.ry != layers[layer]->the_fit.ry)) {
				if (gfit.rx < layers[layer]->the_fit.rx ||
						gfit.ry < layers[layer]->the_fit.ry) {
					siril_log_message(_("The first loaded image should have the greatest sizes for now\n"));
					sprintf(buf, _("NOT OK %ux%u"), layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
					gtk_label_set_text(layers[layer]->label, buf);
					gtk_widget_set_tooltip_text(GTK_WIDGET(layers[layer]->label), _("The first loaded image should have the greatest sizes for now"));
					retval = 1;
				} else {
					siril_log_message(_("Resizing the loaded image from %dx%d to %dx%d\n"),
							layers[layer]->the_fit.rx,
							layers[layer]->the_fit.ry, gfit.rx, gfit.ry);
						sprintf(buf, _("OK upscaled from %ux%u"),
								layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
						cvResizeGaussian(&layers[layer]->the_fit, gfit.rx, gfit.ry, OPENCV_LANCZOS4, TRUE);
						gtk_label_set_text(layers[layer]->label, buf);
						layers[layer]->center.x = layers[layer]->the_fit.rx / 2.0;
						layers[layer]->center.y = layers[layer]->the_fit.ry / 2.0;
						gtk_widget_set_tooltip_text(GTK_WIDGET(layers[layer]->label), _("Image loaded, and upscaled"));
				}
			}
			else if (!retval) {
				layers[layer]->center.x = layers[layer]->the_fit.rx / 2.0;
				layers[layer]->center.y = layers[layer]->the_fit.ry / 2.0;
				sprintf(buf, _("OK %ux%u"), layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
				gtk_label_set_text(layers[layer]->label, buf);
				gtk_widget_set_tooltip_text(GTK_WIDGET(layers[layer]->label), _("Image loaded"));
			}
		}
	}
	g_free(filename);

	/* special case of luminance selected */
	if (layer == 0) {
		GtkToggleButton *lum_button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "composition_use_lum"));
		g_signal_handlers_block_by_func(lum_button, on_composition_use_lum_toggled, NULL);
		gtk_toggle_button_set_active(lum_button, !retval);
		g_signal_handlers_unblock_by_func(lum_button, on_composition_use_lum_toggled, NULL);
		luminance_mode = !retval;
	}
	if (retval) {
		clearfits(&layers[layer]->the_fit);
		orig_rx[layer] = 0;
		orig_ry[layer] = 0;
		return;
	}
	// Check if we have a color image. If so, show the RGB viewport. But set a variable so if the
	// user switches back to a mono viewport and then loads another layer it won't change again -
	// at that point it's up to them.
	if (number_of_images_loaded() > 1) {
		GtkNotebook* notebook = (GtkNotebook*) lookup_widget("notebook1");
		gtk_notebook_set_current_page(notebook, 3);
		gui.cvport = 3;
		redraw(REMAP_ALL);
		update_display_selection();	// update the dimensions of the selection when switching page
		update_display_fwhm();
	}

	update_compositing_registration_interface();

	update_result(0);
	update_metadata(TRUE);
	gui_function(update_MenuItem, NULL);
}

void create_the_internal_sequence() {
	int i, j, nb_layers;
	if (seq) free_sequence(seq, TRUE);

	nb_layers = number_of_images_loaded();
	if (nb_layers == 0 || nb_layers == 1) {
		char *msg = siril_log_message(_("You must at least load two layers before!\n"));
		siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
		seq = NULL;
		return;
	}
	if (luminance_mode) {
		seq = create_internal_sequence(nb_layers);
		for (j = 0, i=0; i<layers_count; i++) {
			if (has_fit(i)) {
				internal_sequence_set(seq, j, &layers[i]->the_fit);
				j++;
			}
		}
	} else {
		if (has_fit(0))
			nb_layers--;
		seq = create_internal_sequence(nb_layers);
		for (j=0, i=1; i<layers_count; i++) {
			if (has_fit(i)) {
				internal_sequence_set(seq, j, &layers[i]->the_fit);
				j++;
			}
		}
	}
	seq->bitpix = gfit.bitpix;
	seq->rx = gfit.rx;
	seq->ry = gfit.ry;
}

void on_centerbutton_toggled(GtkToggleButton *button, gpointer user_data) {
	// First, get the state of this centerbutton
	gboolean state = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
	// If state == active, ensure all the other centerbuttons are inactivated
	if (state) {
		for (int layer = 0 ; layer < maximum_layers ; layer++) {
			if (layers[layer]) {
				if (layers[layer]->centerbutton != button) {
					// Suppress the callback while deactivating each of the other buttons
					g_signal_handlers_block_by_func(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), on_centerbutton_toggled, NULL);
					gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), FALSE);
					g_signal_handlers_unblock_by_func(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), on_centerbutton_toggled, NULL);
				} else {
					// Set gui.comp_layer_centering to the current layer
					gui.comp_layer_centering = layers[layer];
					redraw(REDRAW_OVERLAY);
				}
			}
		}
		// Set the mouse_status to allow the image_interactions callbacks to set the center point
		mouse_status = MOUSE_ACTION_GET_COMP_CENTER_COORDINATE;
	} else {
		// Reset the mouse status
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	}
}


/* start aligning the layers: create an 'internal' sequence and run the selected method on it */
void on_button_align_clicked(GtkButton *button, gpointer user_data) {
	// Set gui.comp_layer_centering to NULL so that center points cease to be drawn
	gui.comp_layer_centering = NULL;
	// Ensure the transformation type is set correctly
	struct registration_method *method;
	framing_type framing;
	GtkComboBox *regcombo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "compositing_align_method_combo"));
	method = reg_methods[gtk_combo_box_get_active(regcombo)];
	GtkComboBox *framingcombo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "compositing_align_framing_combo"));
	int ft = gtk_combo_box_get_active(framingcombo);
	framing = ft == 0 ? FRAMING_CURRENT : ft == 1 ? FRAMING_MIN : FRAMING_COG;
	if (method->method_ptr == register_shift_fwhm || method->method_ptr == register_shift_dft)
		the_type = SHIFT_TRANSFORMATION;
	else
		the_type = HOMOGRAPHY_TRANSFORMATION;

	gboolean two_step = (method->method_ptr == register_multi_step_global ||
		method->method_ptr == register_kombat || method->method_ptr == register_manual) ? TRUE : FALSE;

	gboolean do_sum = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("cumulate_rgb_button")));

	// Avoid crash if gfit has been closed since populating the layers
	if (!gfit.data && !gfit.fdata) {
		update_result(0);
	} else {
		// Avoid crash if the size of gfit is larger than the input images
		// This can happen if compositing is run twice in a row with a
		// framing type that increases its dimensions
		check_gfit_is_ours();
	}

	// Check the input images are still all the same size. Size changes
	// can happen if a prior composition failed to align some images.
	gboolean variable = FALSE;
	int start = luminance_mode ? 0 : 1;
	int rx = layers[start]->the_fit.rx;
	int ry = layers[start]->the_fit.ry;
	for (int i = start ; i < layers_count ; i++) {
		if (layers[i]->the_fit.rx != 0 && (rx != layers[i]->the_fit.rx || ry != layers[i]->the_fit.ry)) {
			char buf[48];
			sprintf(buf, _("NOT OK: %ux%u"), layers[i]->the_fit.rx, layers[i]->the_fit.ry);
			gtk_label_set_text(layers[i]->label, buf);
			variable = TRUE;
		}
	}
	if (variable) {
		siril_log_color_message(_("The image sizes are not consistent. This may happen as the result of previous alignment operations. Please reload the input images.\n"), "red");
		return;
	}

	int i = 0;
	struct registration_args regargs = { 0 };
	char *msg;

	create_the_internal_sequence();
	/* align it */

	regargs.seq = seq;
	regargs.no_output = FALSE;
	get_the_registration_area(&regargs, method);
	regargs.layer = 0;
	seq->reference_image = 0;
	regargs.max_stars_candidates = MAX_STARS_FITTED;
	regargs.run_in_thread = FALSE;
	regargs.interpolation = OPENCV_LANCZOS4;
	regargs.output_scale = 1.f;
	regargs.clamp = TRUE;
	regargs.framing = framing;
	regargs.percent_moved = 0.50f; // Only needed for KOMBAT
	regargs.two_pass = (method->method_ptr == register_multi_step_global &&
						framing != FRAMING_CURRENT) ? TRUE : FALSE;
	regargs.type = HOMOGRAPHY_TRANSFORMATION;
	com.run_thread = TRUE;	// fix for the cancelling check in processing

	// Update the spinbutton values if we are doing manual reg
	// This avoids GTK calls from threads
	if (method->method_ptr == register_manual) {
		for (int layer = 0 ; layer < maximum_layers ; layer++) {
			if (layers[layer]) {
				layers[layer]->spinbutton_x_value = gtk_spin_button_get_value(layers[layer]->spinbutton_x);
				layers[layer]->spinbutton_y_value = gtk_spin_button_get_value(layers[layer]->spinbutton_y);
				layers[layer]->spinbutton_r_value = gtk_spin_button_get_value(layers[layer]->spinbutton_r);
			}
		}
	}

	msg = siril_log_message(_("Starting registration using method: %s\n"), method->name);
	msg[strlen(msg)-1] = '\0';
	set_cursor_waiting(TRUE);
	set_progress_bar_data(msg, PROGRESS_RESET);
	int ret1 = method->method_ptr(&regargs);
	free(regargs.imgparam);
	regargs.imgparam = NULL;
	free(regargs.regparam);
	regargs.regparam = NULL;
	if (ret1) {
		set_progress_bar_data(_("Error in layers alignment."), PROGRESS_DONE);
		set_cursor_waiting(FALSE);
		com.run_thread = FALSE;	// fix for the cancelling check in processing
		return;
	}
	if (two_step) {
		int count_A = 0;
		int count_B = 0;
		for (int index = 0 ; index < layers_count - start ; index++) {
			int layer = index + start;
			if (layers[layer]->the_fit.rx != 0) {
				count_B++;
				if (seq->imgparam[index].incl) count_A++;
			}
		}
		if (count_A != count_B) {
			if (!siril_confirm_dialog(_("Incomplete alignment"),
					_("Some images did not align correctly. Proceed to see the "
					"partially aligned result? (This may alter image dimensions "
					"in which case the images must be re-loaded to retry "
					"alignment.)"), _("Proceed"))) {
				set_cursor_waiting(FALSE);
				com.run_thread = FALSE;	// fix for the cancelling check in processing
				return;
			}
		}
		if (luminance_mode) {
			// Set the reference image to 0 (i.e. Luminance)
			regargs.seq->reference_image = 0;
			siril_log_message(_("Using Luminance channel as reference.\n"));
		}
		int ret2 = register_apply_reg(&regargs);
		free(regargs.imgparam);
		regargs.imgparam = NULL;
		free(regargs.regparam);
		regargs.regparam = NULL;
		if (ret2) {
			set_progress_bar_data(_("Error in layers alignment."), PROGRESS_DONE);
			set_cursor_waiting(FALSE);
			com.run_thread = FALSE;	// fix for the cancelling check in processing
			return;
		}
	}
	// update WCS etc.
	update_comp_metadata(seq->internal_fits[seq->reference_image], do_sum);
	set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	com.run_thread = FALSE;	// fix for the cancelling check in processing

	/* display the values */
	if (method->method_ptr != register_manual) {
		if (!luminance_mode)
			i = 1;

		for (int j=0; i<layers_count; i++) {
			if (has_fit(i)) {
				double dx, dy, rotation;
				translation_from_H(seq->regparam[0][j].H, &dx, &dy);
				rotation = atan2(seq->regparam[0][j].H.h01, seq->regparam[0][j].H.h00) * 180 / M_PI;
				gtk_spin_button_set_value(layers[i]->spinbutton_x, dx);
				gtk_spin_button_set_value(layers[i]->spinbutton_y, dy);
				gtk_spin_button_set_value(layers[i]->spinbutton_r, rotation);
				j++;
			}
		}
	}
	/* align the image and display it.
	 * Layers are aligned against the reference layer, with zeros where there is not data */
	update_result(1);
	// reset the transformation type so that it is always in this state by default
	the_type = HOMOGRAPHY_TRANSFORMATION;
	// Reset rotation centers: owing to the change of framing the previous rotation centers
	// cannot be relied on. Note this means only one rotation should be carried out at a time.
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layers[layer]) {
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), FALSE);
			layers[layer]->center.x = gfit.rx / 2.0;
			layers[layer]->center.y = gfit.ry / 2.0;
		}
	}
}

float get_normalized_pixel_value(int fits_index, float layer_pixel_value) {
	double tmp = (double)layer_pixel_value;
	if (!has_fit(0))
		fits_index--;

	tmp *= coeff->scale[fits_index];
	tmp -= coeff->offset[fits_index];
	return tmp;
}

/* get the pixel value at coordinates x,y for image in layers[fits_index]->the_fit.
 * x and y are given in buffer coordinates, not image coordinates.
 * Handles (not yet - binning and) registration offset */
static float get_composition_pixel_value(int fits_index, int reg_layer, int x, int y) {
	int realX = x, realY = y;
	if (seq && seq->regparam && reg_layer < seq->number && reg_layer >= 0) {
		double dx = 0.0, dy = 0.0;
		// Not needed except for shift transformation
		if (the_type == SHIFT_TRANSFORMATION)
			translation_from_H(seq->regparam[0][reg_layer].H, &dx, &dy);
		// all images have one layer, hence the 0 below
		realX = x - round_to_int(dx);
		if (realX < 0 || realX >= layers[fits_index]->the_fit.rx) return 0.0f;
		realY = y - round_to_int(dy);
		if (realY < 0 || realY >= layers[fits_index]->the_fit.ry) return 0.0f;
	}
	float pixel_value;
	if (layers[fits_index]->the_fit.type == DATA_FLOAT)
		pixel_value = layers[fits_index]->the_fit.fpdata[0][realX + realY * layers[fits_index]->the_fit.rx];
	else if (layers[fits_index]->the_fit.type == DATA_USHORT)
		pixel_value = (float) layers[fits_index]->the_fit.pdata[0][realX + realY * layers[fits_index]->the_fit.rx] / USHRT_MAX_SINGLE;
	else
		pixel_value = 0.f;
	if (coeff) {
		// normalization
		pixel_value = get_normalized_pixel_value(fits_index, pixel_value);
	}
	return pixel_value;
}

/* increments the color values in rgbpixel from the pixel value for a particular
 * layer. GdkRGBA values are stored in the [0, 1] interval. */
static void increment_pixel_components_from_layer_value(int fits_index, GdkRGBA *rgbpixel, float layer_pixel_value) {
	GdkRGBA *layer_color = &layers[fits_index]->color;
	rgbpixel->red += layer_color->red * layer_pixel_value;
	rgbpixel->green += layer_color->green * layer_pixel_value;
	rgbpixel->blue += layer_color->blue * layer_pixel_value;
}

/* increments the color values in rgbpixel from the saturated pixel value for a
 * particular layer. GdkRGBA values are stored in the [0, 1] interval. */
static void increment_pixel_components_from_layer_saturated_value(int fits_index, GdkRGBA *rgbpixel, float layer_pixel_value) {
	GdkRGBA *layer_color = &layers[fits_index]->saturated_color;
	if (layer_pixel_value > 1.0f) {
		/* images could have pixel values above 1, especially when
		 * demosaicing is used, we shouldn't count them as overflow
		 * here */
		layer_pixel_value = 1.0f;
	}
	rgbpixel->red += layer_color->red * layer_pixel_value;
	rgbpixel->green += layer_color->green * layer_pixel_value;
	rgbpixel->blue += layer_color->blue * layer_pixel_value;
}

/* called when selection changed */
static void update_compositing_registration_interface() {
	if (!gui.builder) return;
	GtkLabel *label = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_msg"));
	GtkComboBox *combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "compositing_align_method_combo"));
	int sel_method = gtk_combo_box_get_active(combo);

	if (com.selection.w <= 0 && com.selection.h <= 0 && (sel_method == 1 || sel_method == 2)) {
		// DFT shift ad KOMBAT require a selection to be made
		gtk_label_set_text(label, _("An image area must be selected for align"));
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);
	/*} else if (ref_layer == -1 || (!luminance_mode && ref_layer == 0)) {
		gtk_label_set_text(label, "A reference layer must be selected for align");
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);*/
	} else if (number_of_images_loaded() < 2) {
		gtk_label_set_text(label, _("At least 2 channels must be loaded for align"));
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);
	} else {
		gtk_label_set_text(label, "");
		gtk_widget_set_sensitive(lookup_widget("button_align"), TRUE);
	}
	gui_function(update_MenuItem, NULL);
}

/* callback for changes of the selected reference layer */

void on_compositing_align_method_combo_changed(GtkComboBox *widget, gpointer user_data) {
	update_compositing_registration_interface();
	int sel_method = gtk_combo_box_get_active(widget);
	if (sel_method == 3) { // Prepare for manual alignment
		for (int layer = 0 ; layer < maximum_layers ; layer++) {
			if (layers[layer]) {
				// Make the spinbuttons and centerbutton sensitive to input
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_x), TRUE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_y), TRUE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_r), TRUE);
				// When changing method we always want to set all the center buttons to be inactive initially
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), FALSE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->centerbutton), TRUE);
				// Zero the spinbuttons
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(layers[layer]->spinbutton_x), 0.0);
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(layers[layer]->spinbutton_y), 0.0);
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(layers[layer]->spinbutton_r), 0.0);
				layers[layer]->center.x = layers[layer]->the_fit.rx / 2.0;
				layers[layer]->center.y = layers[layer]->the_fit.ry  /2.0;
			}
		}
	} else {
		gui.comp_layer_centering = NULL;
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		for (int layer = 0 ; layer < maximum_layers ; layer++) {
			if (layers[layer]) {
				// Make the spinbuttons and centerbutton insensitive to input
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_x), FALSE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_y), FALSE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->spinbutton_r), FALSE);
				gtk_widget_set_sensitive(GTK_WIDGET(layers[layer]->centerbutton), FALSE);
				// When changing method we always want to set all the center buttons to be inactive initially
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(layers[layer]->centerbutton), FALSE);
			}
		}
	}
}

void on_composition_combo_coloringtype_changed(GtkComboBox *widget, gpointer user_data) {
	coloring_type = gtk_combo_box_get_active(widget);
	update_result(0);
}

/* Image composition without luminance. Used for RGB composition for example.
 * Result is in gfit. */
static void colors_align_and_compose() {
	int x, y;
	if (no_color_available()) return;
	// Sort the date_obs and pic the earliest one
	GList *date_obs_list = NULL;
	for (int layer = 1 ; layers[layer]; layer++) {
		if (has_fit(layer) && layers[layer]->the_fit.keywords.date_obs) {
			date_obs_list = g_list_append(date_obs_list, layers[layer]->the_fit.keywords.date_obs);
		}
	}
	int len = g_list_length(date_obs_list);
	if (len > 0) {
		date_obs_list = g_list_sort(date_obs_list, g_date_time_compare);
		GDateTime *earliest = (GDateTime*) g_list_nth(date_obs_list, 0)->data;
		GDateTime *latest = (GDateTime*) g_list_nth(date_obs_list, len-1)->data;
		GTimeSpan timespan = g_date_time_difference(latest, earliest);
		timespan /= 3600000000;
		if (timespan > 24 && !timespan_warning_given) {
			siril_log_message(_("Attention: channels are dated more than 24 hours apart. DATE_OBS "
								"set to the earliest observation start date but for "
								"some purposes this field may be of limited use (e.g. solar "
								"system objects may have moved significantly between channels).\n"));
			timespan_warning_given = TRUE;
		} else {
			timespan_warning_given = FALSE;
		}
		if (gfit.keywords.date_obs) {
			g_date_time_unref(gfit.keywords.date_obs);
		}
		gfit.keywords.date_obs = g_date_time_ref(earliest);
	}
	g_list_free(date_obs_list);
	fprintf(stdout, "colour layers only composition\n");
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y,x) schedule(static)
#endif
	for (y = 0; y < gfit.ry; ++y) {
		for (x = 0; x < gfit.rx; ++x) {
			int layer;
			GdkRGBA pixel;
			clear_pixel(&pixel);
			for (layer = 1; layers[layer]; layer++) {
				if (has_fit(layer)) {
					int reg_layer = seq ? internal_sequence_find_index(seq, &layers[layer]->the_fit) : -1;
					float layer_value = get_composition_pixel_value(layer, reg_layer, x, y);
					if (layer_value > 0.0f)
						increment_pixel_components_from_layer_value(layer, &pixel, layer_value);
				}
			}

			rgb_pixel_limiter(&pixel);
			size_t dst_index = y * gfit.rx + x;
			gfit.fpdata[RLAYER][dst_index] = pixel.red;
			gfit.fpdata[GLAYER][dst_index] = pixel.green;
			gfit.fpdata[BLAYER][dst_index] = pixel.blue;
		}
	}
}

/* This function fills the data in the gfit image with LRGB information from
 * layers[*]->the_fit images. Layers are aligned with registration data, no
 * binning yet. */
static void luminance_and_colors_align_and_compose() {
	/* Each pixel is transformed from RGB to HSI, I is replaced by the
	 * luminance layer's value and transformed back to RGB. */
	guint x, y;
	assert(has_fit(0));
	// Copy the date_obs field from the luminance layer
	if (layers[0]->the_fit.keywords.date_obs) {
		if (gfit.keywords.date_obs) {
			g_date_time_unref(gfit.keywords.date_obs);
		}
		g_date_time_ref(layers[0]->the_fit.keywords.date_obs);
		gfit.keywords.date_obs = layers[0]->the_fit.keywords.date_obs;
	}
	if (no_color_available()) {
		/* luminance only: we copy its data to all result layers */
		int i;
		size_t nbdata = gfit.rx * gfit.ry;
		fprintf(stdout, "luminance-only, no composition\n");
		for (i=0; i<3; i++)
			memcpy(gfit.fpdata[i], layers[0]->the_fit.fdata, nbdata*sizeof(float));
		return;
	}
	fprintf(stdout, "luminance-enabled composition\n");

	image_find_minmax(&layers[0]->the_fit);
	double norm = (double)(layers[0]->the_fit.maxi);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y,x) schedule(static)
#endif
	for (y = 0; y < gfit.ry; y++) {
		for (x = 0; x < gfit.rx; x++) {
			int layer;
			gdouble h, s, i;
			gdouble X, Y, Z;
			gdouble a, b;
			/* get color information */
			GdkRGBA pixel;
			clear_pixel(&pixel);
			for (layer = 1; layers[layer]; layer++) {
				if (has_fit(layer)) {
					float layer_value = get_composition_pixel_value(layer, layer, x, y);
					if (layer_value > 0.0f)
						increment_pixel_components_from_layer_value(layer, &pixel, layer_value);
				}
			}
			rgb_pixel_limiter(&pixel);

			switch (coloring_type) {
			case HSL:
				rgb_to_hsl(pixel.red, pixel.green, pixel.blue, &h, &s, &i);
				/* add luminance by replacing it in the HSI */
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				/* converting back to RGB */
				hsl_to_rgb(h, s, i, &pixel.red, &pixel.green, &pixel.blue);
				break;
			case HSV:
				rgb_to_hsv(pixel.red, pixel.green, pixel.blue, &h, &s, &i);
				/* add luminance by replacing it in the HSI */
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				/* converting back to RGB */
				hsv_to_rgb(h, s, i, &pixel.red, &pixel.green, &pixel.blue);
				break;
			case CIELAB:
				rgb_to_xyz(pixel.red, pixel.green, pixel.blue, &X, &Y, &Z);
				xyz_to_LAB(X, Y, Z, &i, &a, &b);
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				i *= 100.0;		// 0 < L < 100
				LAB_to_xyz(i, a, b, &X, &Y, &Z);
				xyz_to_rgb(X, Y, Z, &pixel.red, &pixel.green, &pixel.blue);
				break;
			}

			rgb_pixel_limiter(&pixel);

			/* and store in gfit */
			size_t dst_index = y * gfit.rx + x;
			gfit.fpdata[RLAYER][dst_index] = pixel.red;
			gfit.fpdata[GLAYER][dst_index] = pixel.green;
			gfit.fpdata[BLAYER][dst_index] = pixel.blue;
		}
	}
}

void on_compositing_cancel_clicked(GtkButton *button, gpointer user_data){
	gui.comp_layer_centering = NULL;
	reset_compositing_module();
	siril_close_dialog("composition_dialog");
}

/* When summing all layers to get the RGB values for one pixel, it may overflow.
 * This procedure defines what happens in that case. */
static void rgb_pixel_limiter(GdkRGBA *pixel) {
#ifdef DEBUG
	if (pixel->red > 1.2f || pixel->green > 1.2f || pixel->blue > 1.2f)
		fprintf(stdout, "large overflow %g,%g,%g\n", pixel->red,
				pixel->green, pixel->blue);
#endif
	if (pixel->red >= 1.0f)
		pixel->red = 1.0f;
	if (pixel->green >= 1.0f)
		pixel->green = 1.0f;
	if (pixel->blue >= 1.0f)
		pixel->blue = 1.0f;
}

/* initializes a GdkRGBA to black */
static void clear_pixel(GdkRGBA *pixel) {
	pixel->red = 0.0f;
	pixel->green = 0.0f;
	pixel->blue = 0.0f;
	pixel->alpha = 1.0f;
}

/* recompute the layer composition and optionnally refresh the displayed result image */
static void update_result(int assume_gfit_is_ours) {
	icc_auto_assign(&gfit, ICC_ASSIGN_ON_COMPOSITION);

	if (!assume_gfit_is_ours)
		check_gfit_is_ours();
	if (luminance_mode && has_fit(0)) {
		luminance_and_colors_align_and_compose();
	} else {
		colors_align_and_compose();
	}
	if (number_of_images_loaded() > 0) {
		notify_gfit_modified();
		redraw(REMAP_ALL);
	}
}

/****************** colour management ******************/

// update the saturated color from the new real colour
static void color_has_been_updated(int layer) {
	double h, s, v;
	GdkRGBA *real = &layers[layer]->color;
	GdkRGBA *satu = &layers[layer]->saturated_color;
	rgb_to_hsv(real->red, real->green, real->blue, &h,&s,&v);
	printf("%d: saturation: %g, light: %g\n", layer, s, v);
	/* in HSL, the actual saturated pure colour happens at l=0.5 and s=1 */
	s = 1.0; v = 1.0;
	hsv_to_rgb(h,s,v, &satu->red, &satu->green, &satu->blue);
	printf("%d: r: %g, g: %g, b: %g\n", layer, satu->red, satu->green, satu->blue);
}

// update a real colour from the saturated colour with a new lightness
static void update_color_from_saturation(int layer, double newl) {
	double h, s, v;
	GdkRGBA *real = &layers[layer]->color;
	GdkRGBA *satu = &layers[layer]->saturated_color;
	rgb_to_hsv(satu->red, satu->green, satu->blue, &h,&s,&v);
	hsv_to_rgb(h,s,newl, &real->red, &real->green, &real->blue);
}

void on_colordialog_response(GtkColorChooserDialog *chooser, gint response_id, gpointer user_data) {
	/* this callback is called on any action of the dialog, and must be
	 * filtered according to the response_id. List is here:
	 * https://developer.gnome.org/gtk3/stable/GtkDialog.html#GTK-RESPONSE-NONE:CAPS
	 */
	if (response_id == GTK_RESPONSE_DELETE_EVENT ||
			response_id == GTK_RESPONSE_CANCEL ||
			response_id == GTK_RESPONSE_CLOSE) {
		current_layer_color_choosing = 0;
		gtk_widget_hide(GTK_WIDGET(chooser));
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		return;
	}

	if (current_layer_color_choosing > 0 && layers[current_layer_color_choosing]) {
		gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(chooser), &layers[current_layer_color_choosing]->display_color);

		// The color chooser returns a RGBA in the monitor color profile. We need to change this to
		// the image color profile to set it accurately for the composition
		float disp[3] = { (float) layers[current_layer_color_choosing]->display_color.red,
						  (float) layers[current_layer_color_choosing]->display_color.green,
						  (float) layers[current_layer_color_choosing]->display_color.blue };
		float img[3];
		cmsHPROFILE image_profile = NULL;
		if (gfit.icc_profile)
			image_profile = copyICCProfile(gfit.icc_profile);
		else if ((com.pref.icc.autoassignment & ICC_ASSIGN_ON_COMPOSITION) && com.icc.working_standard)
			image_profile = copyICCProfile(com.icc.working_standard);
		if (image_profile) {
			// We use the NONEGATIVES flag to bound the transform, otherwise the
			// negative components play havoc with the compositing.
			cmsHTRANSFORM transform = cmsCreateTransformTHR(com.icc.context_single,
															gui.icc.monitor,
															TYPE_RGB_FLT_PLANAR,
															image_profile,
															TYPE_RGB_FLT_PLANAR,
															INTENT_RELATIVE_COLORIMETRIC,
															cmsFLAGS_NONEGATIVES);
			cmsCloseProfile(image_profile);
			if (transform) {
				cmsDoTransform(transform, disp, img, 1);
				cmsDeleteTransform(transform);
				siril_debug_print("Color picker transform (display to image): R: %f -> %f, G: %f -> %f, B: %f -> %f\n", disp[0], img[0], disp[1], img[1], disp[2], img[2]);
			} else {
				memcpy(&img, &disp, 3 * sizeof(float));
				siril_debug_print("Unable to complete color picker transform\n");
			}
		} else {
			memcpy(&img, &disp, 3 * sizeof(float));
		}
		layers[current_layer_color_choosing]->color.red = img[0];
		layers[current_layer_color_choosing]->color.green = img[1];
		layers[current_layer_color_choosing]->color.blue = img[2];

		color_has_been_updated(current_layer_color_choosing);
		gtk_widget_queue_draw(GTK_WIDGET(layers[current_layer_color_choosing]->color_w));
		gtk_widget_hide(GTK_WIDGET(chooser));
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		if (has_fit(current_layer_color_choosing))
			update_result(0);
	}
}

/* managing the drawing area that displays the color and opens the color dialog when clicked */
gboolean draw_layer_color(GtkDrawingArea *widget, cairo_t *cr, gpointer data) {
	int layer, w, h;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;

	w = gtk_widget_get_allocated_width(GTK_WIDGET(widget));
	h = gtk_widget_get_allocated_height(GTK_WIDGET(widget));
	cairo_set_source_rgb(cr, layers[layer]->display_color.red,
			layers[layer]->display_color.green, layers[layer]->display_color.blue);
	//cairo_rectangle(cr, (double)w*0.33, 1, (double)w*0.33, (double)h-2.0);
	cairo_rectangle(cr, 1, 1, w-2.0, h-2.0);
	cairo_fill(cr);
	return FALSE;
}

/* click on the colored area: button press, only configure for quick color edit with
 * right button */
gboolean on_color_button_press_event(const GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data) {
	int layer;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;
	if (event->button == GDK_BUTTON_SECONDARY) {	// right click
		current_layer_color_choosing = layer;
		color_quick_edit = 1;
		memcpy(&qe_ref_color, &layers[layer]->color, sizeof(GdkRGBA));
		return TRUE;
	}
	return FALSE;
}

/* click on the colored area: open the color chooser dialog */
gboolean on_color_button_release_event(const GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data) {
	int layer;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;
	if (event->button == GDK_BUTTON_PRIMARY) {	// left click
		current_layer_color_choosing = layer;
		gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(color_dialog), &layers[layer]->color);
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		gtk_widget_show(GTK_WIDGET(color_dialog));
	} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
		if (has_fit(current_layer_color_choosing))
			update_result(0);
		current_layer_color_choosing = 0;
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
	}
	color_quick_edit = 0;
	return TRUE;
}

gboolean on_color_button_motion_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data) {
	if (color_quick_edit) {	// right click
		double h, s, v;
		//fprintf(stdout, "%g\n", event->x);
		rgb_to_hsv(qe_ref_color.red, qe_ref_color.green,
				qe_ref_color.blue, &h,&s,&v);
		h += event->x / 600.0;
		v -= event->y / 600.0;
		while (h < 0.0) h += 1.0;
		while (h > 1.0) h -= 1.0;
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		hsv_to_rgb(h,s,v, &layers[current_layer_color_choosing]->color.red,
				&layers[current_layer_color_choosing]->color.green,
				&layers[current_layer_color_choosing]->color.blue);
		color_has_been_updated(current_layer_color_choosing);
		gtk_widget_queue_draw(GTK_WIDGET(layers[current_layer_color_choosing]->color_w));
	}
	return FALSE;
}
/*******************************************************/

/* fill the combo box containing filter names */
static void populate_filter_lists() {
	int i, nb_filters;
	GtkComboBoxText *cbox = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboboxtext_filters"));
	nb_filters = get_nb_narrow_filters();
	gtk_combo_box_text_remove_all(cbox);
	for (i=0; i<nb_filters; i++)
		gtk_combo_box_text_append_text(cbox, narrow_band_filters[i].name);
	//gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
}

/* the combo box containing filter names has one item selected: update colours */
void on_filter_changed(GtkComboBox *widget, gpointer user_data) {
	gint active = gtk_combo_box_get_active(widget);
	char wl_text[20];
	if (active == -1) return;
	sprintf(wl_text, "%g", narrow_band_filters[active].wavelength);
	gtk_entry_set_text(wl_entry, wl_text);
}

void on_wavelength_changed(GtkEditable *editable, gpointer user_data){
	GdkRGBA color;
	double wavelength = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(editable)), NULL);
	if (wavelength < 380.0 || wavelength > 780.0) return;
	wavelength_to_display_RGB(wavelength, &color);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(color_dialog), &color);
}

void reset_compositing_module() {
	timespan_warning_given = FALSE;
	gui.comp_layer_centering = NULL;
	if (!compositing_loaded)
		return;

	if (has_fit(0))
		clearfits(&layers[0]->the_fit);
	for (int i = 1; layers[i]; i++) {
		if (has_fit(i))
			clearfits(&layers[i]->the_fit);
		grid_remove_row(i, 1);
		free(layers[i]);
		layers[i] = NULL;
		layers_count--;
	}

	if (reference)
		cmsCloseProfile(reference);
	reference = NULL;

	/* Reset GtkFileChooserButton Luminance */
	GtkFileChooserButton *lum = GTK_FILE_CHOOSER_BUTTON(gtk_builder_get_object(gui.builder, "filechooser_lum"));
	gtk_file_chooser_unselect_all (GTK_FILE_CHOOSER (lum));

	if (seq) {
		free_sequence(seq, TRUE);
		seq = NULL;
	}

	int i;
	for (i = 1; i < 4; i++) {
		/* Create the three default layers */
		on_layer_add(NULL, NULL);
	}
	layers[i] = NULL;

	gtk_widget_hide(GTK_WIDGET(color_dialog));
	current_layer_color_choosing = 0;

	luminance_mode = 0;
	GtkToggleButton *lum_button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "composition_use_lum"));
	gtk_toggle_button_set_active(lum_button, 0);
	gtk_label_set_text(layers[0]->label, _("not loaded"));
	gtk_widget_set_tooltip_text(GTK_WIDGET(layers[0]->label), _("not loaded"));

	update_compositing_registration_interface();
	on_compositing_align_method_combo_changed((GtkComboBox*)lookup_widget("compositing_align_method_combo"), NULL);
}

void on_compositing_reset_clicked(GtkButton *button, gpointer user_data){
	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;
	process_close(0);
	reset_compositing_module();
	open_compositing_window();	// update the CWD just in case
}

/* Reduce brightness of colours associated to layers so that they never overflow on composition.
 * Algorithm: take the max of the composition and normalize brightness with
 * this max. It has to be done three times because it's the same layers that
 * can act on the resulting RGB channels.
 * This algorithm doesn't give the optimal answer, which could be found
 * iteratively, but it should never give an overflow.
 */
void autoadjust(int force_redraw) {
	int layer, nb_images_red = 0, nb_images_green = 0, nb_images_blue = 0;
	GdkRGBA max_pixel;

	set_cursor_waiting(TRUE);
	clear_pixel(&max_pixel);
	/* sum the max per channel */
	for (layer = 1; layers[layer]; layer++) {
		if (has_fit(layer)) {
			image_find_minmax(&layers[layer]->the_fit);
			double max_value = layers[layer]->the_fit.maxi;
			if (coeff)
				max_value = get_normalized_pixel_value(layer, max_value);
			increment_pixel_components_from_layer_saturated_value(
					layer, &max_pixel, max_value);

			if (layers[layer]->color.red > 0.0f) nb_images_red++;
			if (layers[layer]->color.green > 0.0f) nb_images_green++;
			if (layers[layer]->color.blue > 0.0f) nb_images_blue++;
		}
	}

	if (max_pixel.red <= 1.0f && max_pixel.green <= 1.0f && max_pixel.blue <= 1.0f) {
		if (force_redraw) {
			siril_log_message(_("No overflow with the current colours, redrawing only\n"));
			update_result(0);
		} else {
			siril_log_message(_("Nothing to adjust, no overflow\n"));
			set_cursor_waiting(FALSE);
		}
		return;
	}

	/* update the real colours of layers from their saturated colour, based
	 * on how much each colour of the composition overflows */
	// amounts of normalization to be done on each layer's image for each channel
	double to_redistribute_red = nb_images_red == 0 ? 0.0 : (max_pixel.red - 1.0) / (double)nb_images_red;
	double to_redistribute_green = nb_images_green == 0 ? 0.0 : (max_pixel.green - 1.0) / (double)nb_images_green;
	double to_redistribute_blue = nb_images_blue == 0 ? 0.0 : (max_pixel.blue - 1.0) / (double)nb_images_blue;
	for (layer = 1; layers[layer]; layer++) {
		if (has_fit(layer)) {
			double to_redistribute = 0.0;	// for this layer

			if (layers[layer]->color.red > 0.0f && to_redistribute_red > 0.0f) {
				to_redistribute = to_redistribute_red;
			}
			/* for each layer, we check if a channel requires the
			 * current layer to be readjusted and we take the more
			 * severe value of all requirements */

			if (layers[layer]->color.green > 0.0f && to_redistribute_green > 0.0) {
				if (to_redistribute_green > to_redistribute)
					to_redistribute = to_redistribute_green;
			}

			if (layers[layer]->color.blue > 0.0f && to_redistribute_blue > 0.0) {
				if (to_redistribute_blue > to_redistribute)
					to_redistribute = to_redistribute_blue;
			}

			siril_log_message(_("Readjusting layer %d to %g times bright\n"),
					layer, 1.0-to_redistribute);
			/* to_redistribute here is the maximum reduction we
			 * need to give to the layer */
			update_color_from_saturation(layer, 1.0 - to_redistribute);
		}
	}

	/* redraw colours and composition */
	for (layer = 1; layers[layer]; layer++) {
		gtk_widget_queue_draw(GTK_WIDGET(layers[layer]->color_w));
	}
	update_result(0);
	set_cursor_waiting(FALSE);
}

void on_compositing_autoadjust_clicked(GtkButton *button, gpointer user_data){
	autoadjust(0);
}

/* Normalization functions, not used anymore */

static void coeff_clear() {
	if (coeff) {
		free(coeff->offset);
		free(coeff->scale);
		free(coeff->mul);
		free(coeff);
		coeff = NULL;
	}
}

void on_compositing_reload_all_clicked(GtkButton *button, gpointer user_data) {
	if (number_of_images_loaded() < 1)
		return;

	// Clear the image data and the sequence, if populated
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layers[layer] && layers[layer]->the_fit.rx != 0)
			clearfits(&layers[layer]->the_fit);
	}
	if (seq) {
		free_sequence(seq, TRUE);
		seq = NULL;
	}
	close_single_image();

	// Ensure the first image that was loaded before is reloaded first, as the first image to be
	// loaded must be the largest dimensions
	int biggest_layer;
	for (biggest_layer = 0 ; biggest_layer < maximum_layers ; biggest_layer++) {
		if (layers[biggest_layer])
			break;
	}
	for (int layer = 1 ; layer < maximum_layers ; layer++) {
		if (orig_rx[layer] > layers[biggest_layer]->the_fit.rx && orig_ry[layer] > layers[biggest_layer]->the_fit.ry)
			biggest_layer = layer;
	}
	on_filechooser_file_set(layers[biggest_layer]->chooser, NULL);
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layers[layer] && layer != biggest_layer) {
			on_filechooser_file_set(layers[layer]->chooser, NULL);
		}
	}
	siril_log_message(_("All images reloaded successfully.\n"));
	update_result(0);
}

void on_compositing_save_all_clicked(GtkButton *button, gpointer user_data) {
	if (number_of_images_loaded() == 0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No layers loaded: no files to save"));
		return;
	}
	int retval = 0;
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layers[layer] && layers[layer]->the_fit.rx != 0) {
			gchar *filename = siril_file_chooser_get_filename(GTK_FILE_CHOOSER(layers[layer]->chooser));
			gchar *basename = g_path_get_basename(filename);
			gchar *prepended_filename = g_strdup_printf("comp_%s", basename);
			retval += savefits(prepended_filename, &layers[layer]->the_fit);
			g_free(prepended_filename);
			g_free(basename);
			g_free(filename);
		}
	}
	if (!retval) {
		siril_log_color_message(_("All layer images saved correctly.\n"), "green");
	} else {
		siril_log_color_message(_("Error encountered in saving one or more files. Check previous log messages for details.\n"), "red");
	}
}

void on_compositing_linear_match_clicked(GtkButton *button, gpointer *user_data) {
	int error = 0;
	if (number_of_images_loaded() < 2) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("At least 2 images must be loaded to perform linear matching."));
		return;
	}
	int ref_layer = -1;
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layers[layer] && layers[layer]->the_fit.rx != 0) {
			ref_layer = layer;
			break;
		}
	}
	if (ref_layer < 0) return;
	fits *ref = &layers[ref_layer]->the_fit;
	for (int layer = 0 ; layer < maximum_layers ; layer++) {
		if (layer == ref_layer)
			continue;
		if (layers[layer] && layers[layer]->the_fit.rx != 0) {
			double a[3] = { 0.0 }, b[3] = { 0.0 };
			double low = 1.e-7;
			double high = 0.875;
			if (!find_linear_coeff(&layers[layer]->the_fit, ref, low, high, a, b, NULL)) {
				apply_linear_to_fits(&layers[layer]->the_fit, a, b);
			} else {
				siril_log_message(_("Unable to find linear match coefficients for layer %d. Skipping...\n"), layer);
				error++;
			}
		}
		set_progress_bar_data(_("Linear match layers"), (double) layer / (double) maximum_layers);
	}
	set_progress_bar_data(_("Linear match layers complete"), 1.0);
	if (error)
		siril_log_color_message(_("Error: failed to match %d layers with the reference layer.\n"), "red", error);
	update_result(0);
}

/*********************************
 *                               *
 * Manual registration functions *
 *                               *
 *********************************
 * This registration function takes the shift and rotation values entered in the
 * spinbutton boxes, generates a homography matrix from them and applies it to the
 * set of images.
 */

int manual_align_prepare_results(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	if (!regargs->no_output) {
		// allocate destination sequence data
		regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
		regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
		if (!regargs->imgparam  || !regargs->regparam) {
			PRINT_ALLOC_ERR;
			return 1;
		}

		if (seq_prepare_hook(args))
			return 1;
	}

	sadata->success = calloc(args->nb_filtered_images, sizeof(BYTE));
	if (!sadata->success) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}

int manual_align_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	fits fit = { 0 };
	/* read the reference frame to get sadata->ref.x and ref.y, and check
	   it isn't a CFA sequence */
	if (seq_read_frame(args->seq, regargs->reference_image, &fit, FALSE, -1)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		return 1;
	}
	if (fit.naxes[2] == 1 && fit.keywords.bayer_pattern[0] != '\0')
		siril_log_color_message(_("Registering a sequence opened as CFA is a bad idea.\n"), "red");

	siril_log_color_message(_("Reference Image:\n"), "green");
	sadata->ref.x = fit.rx;
	sadata->ref.y = fit.ry;
	// For internal sequences the data / fdata pointer still
	// points to the original memory in seq->internal_fits.
	// It must not be freed by clearfits here so we set the
	// pointers in fit to NULL
	if (args->seq->type == SEQ_INTERNAL) {
		fit.data = NULL;
		fit.fdata = NULL;
	}
	clearfits(&fit);

	int retval = manual_align_prepare_results(args);
	if (!retval)
		sadata->current_regdata = regargs->regparam;
	return retval;
}

Homography H_from_translation_and_rotation(double dx, double dy, double dr, point center) {
	double alpha = cos(dr);
	double beta = sin(dr);
	Homography H = { alpha,	beta,	dx + (1 - alpha) * center.x - beta * center.y,
					 -beta,	alpha,	dy + beta * center.x + (1 - alpha) * center.y,
					 0,		0,		1 };
	return H;
}

int manual_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int offset = luminance_mode ? 0 : 1;
	if (regargs->no_output) {
		/* if "save transformation only", we choose to initialize all frames
		 * to exclude status. If registration is ok, the status is
		 * set to include */
		args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
	}
	if (layers[in_index + offset] && layers[in_index + offset]->the_fit.rx != 0) {
		double dx = layers[in_index + offset]->spinbutton_x_value;
		double dy = layers[in_index + offset]->spinbutton_y_value;
		double dr = layers[in_index + offset]->spinbutton_r_value;
		// Convert dr to radians
		dr *= M_PI;
		dr /= 180.0;
		point* center = &layers[in_index + offset]->center;
		if (dx == 0.0 && dy == 0.0 && dr == 0.0) {
			regargs->reference_image = in_index;
		}
		regargs->seq->regparam[regargs->layer][in_index].H = H_from_translation_and_rotation(dx, dy, dr, *center);
	} else {
		args->seq->imgparam[in_index].incl = FALSE;
		return 1;
	}
	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	regargs->imgparam[out_index].rx = sadata->ref.x;
	regargs->imgparam[out_index].ry = sadata->ref.y;
	cvGetEye(&regargs->regparam[out_index].H);

	sadata->success[out_index] = 1;
	return 0;
}

int manual_align_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int failed = 0;
	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);
	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;

		if (!regargs->no_output) {
			if (failed) {
				// regargs->imgparam and regargs->regparam may have holes caused by images
				// that failed to be registered - compact them
				for (int i = 0, j = 0; i < regargs->new_total; i++, j++) {
					while (!sadata->success[j] && j < args->nb_filtered_images) j++;
					g_assert(sadata->success[j]);
					if (i != j) {
						regargs->imgparam[i] = regargs->imgparam[j];
						regargs->regparam[i] = regargs->regparam[j];
					}
				}
			}

			seq_finalize_hook(args);
		}
	} else {
		regargs->new_total = 0;
		free(args->seq->regparam[regargs->layer]);
		args->seq->regparam[regargs->layer] = NULL;

		if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		}
	}

	if (sadata->success) free(sadata->success);
	free(sadata);
	args->user = NULL;
	if (!args->retval) {
		siril_log_message(_("Registration finished.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->new_total);

		g_free(str);
	}
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	return regargs->new_total == 0;
}

int manual_align_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	if (limit > com.max_thread)
		limit = com.max_thread;
	return limit;
}

int register_manual(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = manual_align_compute_mem_limits;
	args->prepare_hook = manual_align_prepare_hook;
	args->image_hook = manual_align_image_hook;
	args->finalize_hook = manual_align_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Manual registration");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = 1.0;
	args->new_seq_prefix = strdupnullok(regargs->prefix);
	args->load_new_sequence = !regargs->no_output;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free_generic_seq_args(args, FALSE);
		return -1;
	}
	sadata->regargs = regargs;
	args->user = sadata;

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free_generic_seq_args(args, FALSE);
	return regargs->retval;
}

