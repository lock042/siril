/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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

#include <gtk/gtk.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "stacking/stacking.h"
#include "algos/planetary.h"
#include "plot.h"
#include "gui.h"
#include "callbacks.h"
#include "image_display.h"
#include "planetary_callbacks.h"

static gboolean seqimage_range_mouse_pressed = FALSE;
static gboolean add_zones_mode = FALSE;
static gboolean remove_zones_mode = FALSE;
static double lowest_accepted_quality = 0.0;

/* the planetary mode uses a slider to change displayed image whereas the
 * deep-sky mode uses a spin button (on_imagenumberspin_output) */
void on_seqimage_changed(GtkRange *range, gpointer user_data) {
	if (!sequence_is_loaded() || seqimage_range_mouse_pressed) return;
	int index = (int)gtk_range_get_value(range);
	seq_load_image(&com.seq, index, TRUE);
}

gboolean on_seqimage_button_press(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	seqimage_range_mouse_pressed = TRUE;
	return FALSE;
}

gboolean on_seqimage_button_release(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	seqimage_range_mouse_pressed = FALSE;
	on_seqimage_changed(GTK_RANGE(widget), NULL);
	return FALSE;
}

void on_bestimage_changed(GtkRange *range, gpointer user_data) {
	if (!sequence_is_loaded()) return;
	lowest_accepted_quality = compute_lowest_accepted_quality(gtk_range_get_value(range));
	plot_set_filtering_threshold(lowest_accepted_quality);
}


void on_add_zones_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	add_zones_mode = gtk_toggle_button_get_active(togglebutton);
	if (add_zones_mode)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("remove_zones")), FALSE);
}

void on_remove_zones_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	remove_zones_mode = gtk_toggle_button_get_active(togglebutton);
	if (remove_zones_mode)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("add_zones")), FALSE);
}

static void activate_mpp_processing_button() {
	static GtkWidget *stack_button = NULL;
	int status = sequence_is_loaded() && refimage_is_set() &&
		com.stacking_zones && com.stacking_zones[0].centre.x >= 0.0;
	// we should also check that there are at least two images to stack
	if (!stack_button)
	       stack_button = lookup_widget("gostack_button");
	gtk_widget_set_sensitive(stack_button, status);
}

void add_stacking_zone(double x, double y, double half_side) {
	if (!com.stacking_zones_size || !com.stacking_zones) {
		com.stacking_zones_size = 40;
		com.stacking_zones = malloc(com.stacking_zones_size * sizeof(stacking_zone));
		com.stacking_zones[0].centre.x = -1.0; 
	}

	int i = 0;
	while (com.stacking_zones[i].centre.x >= 0.0 && i < com.stacking_zones_size - 1)
		i++;
	if (i == com.stacking_zones_size - 1) {
		com.stacking_zones_size *= 2;
		com.stacking_zones = realloc(com.stacking_zones,
				com.stacking_zones_size * sizeof(stacking_zone));
	}

	stacking_zone *zone = &com.stacking_zones[i];
	zone->centre.x = x;
	zone->centre.y = y;
	zone->half_side = half_side;

	com.stacking_zones[i+1].centre.x = -1.0; 

	activate_mpp_processing_button();
}

gboolean on_remove_all_zones_clicked(GtkButton *button, gpointer user_data) {
	if (com.stacking_zones) {
		int i = 0;
		while (com.stacking_zones[i].centre.x >= 0.0) {
			com.stacking_zones[i].centre.x = -1.0;
			i++;
		}
	}
	redraw(com.cvport, REMAP_NONE);
	return FALSE;
}

void planetary_click_in_image(double x, double y) {
	int i = 0;
	if (add_zones_mode) {
		GtkAdjustment *sizeadj = GTK_ADJUSTMENT(lookup_widget("adjustment_zonesize"));
		double size = gtk_adjustment_get_value(sizeadj);
		if (size > 0.0)
			add_stacking_zone(x, y, size);
	}
	else if (remove_zones_mode) {
		double closest_distance = DBL_MAX;
		int closest_zone = -1;
		if (com.stacking_zones) {
			while (com.stacking_zones[i].centre.x >= 0.0) {
				stacking_zone *zone = &com.stacking_zones[i];
				double xdist = zone->centre.x - x;
				double ydist = zone->centre.y - y;
				double distance = sqrt(xdist * xdist + ydist * ydist);
				// distance is good for circular zones, these are square
				// but since a point can be in several zones, distance is
				// a good first indicator
				if (distance > zone->half_side * 1.414) {
					i++;
					continue;
				}
				if (distance < closest_distance &&
						point_is_inside_zone((int)x, (int)y, zone)) {
					closest_distance = distance;
					closest_zone = i;
				}
				i++;
			}

			if (closest_zone != -1) {
				// replace by last zone
				memcpy(&com.stacking_zones[closest_zone],
						&com.stacking_zones[i-1],
						sizeof(stacking_zone));
				com.stacking_zones[i-1].centre.x = -1.0;

				activate_mpp_processing_button();
			}
		}
	}
}

gboolean on_planetary_processing_button_clicked(GtkButton *button, gpointer user_data) {
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	GtkEntry *output_file = GTK_ENTRY(gtk_builder_get_object(builder, "entryresultfile"));
	GtkToggleButton *overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutoverwrite"));

	struct mpr_args *args = malloc(sizeof(struct mpr_args));
	args->seq = &com.seq;
	args->layer = gtk_combo_box_get_active(cbbt_layers);
	args->filtering_criterion = stack_filter_quality;
	args->filtering_parameter = lowest_accepted_quality;
	args->nb_closest_AP = 5;
	args->max_distance = 350.0;
	args->own_distance_f = 0.5;
	args->output_filename = strdup(gtk_entry_get_text(output_file));
	args->output_overwrite = gtk_toggle_button_get_active(overwrite);
	start_in_new_thread(the_multipoint_processing, args);
	return FALSE;
}

/* when the reglayer is changed we update the status of the reference image */
void on_comboboxreglayer_planetary_changed(GtkComboBox *widget, gpointer user_data) {
	update_refimage_on_layer_change(&com.seq, gtk_combo_box_get_active(widget));
	activate_mpp_processing_button();
	display_refimage_if_needed();
}

void on_showref_check_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	int condition = refimage_is_set() && gtk_toggle_button_get_active(togglebutton);
	if (condition) {
		if (single_image_is_loaded()) return;
		clearfits(&gfit);
		const fits *refimage = get_refimage();
		copyfits(refimage, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		// code from end_stacking
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup(_("Stacking result image"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		com.uniq->filename = strdup(get_refimage_filename());

		com.seq.current = RESULT_IMAGE;
		adjust_cutoff_from_updated_gfit();
		set_sliders_value_to_gfit();
		initialize_display_mode();

		sliders_mode_set_state(com.sliders);
		set_cutoff_sliders_max_values();

		set_display_mode();
		display_filename();
	} else {
		if (!single_image_is_loaded() || !sequence_is_loaded()) return;
		if (com.seq.current != RESULT_IMAGE) return;
		// close_single_image does not display a sequence image, it had
		// to be modified as below
		if (com.uniq->filename)
			free(com.uniq->filename);
		if (com.uniq->comment)
			free(com.uniq->comment);
		free(com.uniq);
		com.uniq = NULL;

		int image_to_load = sequence_find_refimage(&com.seq);
		seq_load_image(&com.seq, image_to_load, TRUE);
	}
	redraw(com.cvport, REMAP_ALL);
}

void display_refimage_if_needed() {
	on_showref_check_toggled(GTK_TOGGLE_BUTTON(lookup_widget("showref_check")), NULL);
}
