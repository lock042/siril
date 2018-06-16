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
#include "io/sequence.h"
#include "stacking/stacking.h"
#include "plot.h"
#include "gui.h"

static gboolean seqimage_range_mouse_pressed = FALSE;
static gboolean add_zones_mode = FALSE;
static gboolean remove_zones_mode = FALSE;

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
	double qual = compute_lowest_accepted_quality(gtk_range_get_value(range));
	plot_set_filtering_threshold(qual);
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
				if (distance > zone->half_side) {
					i++;
					continue;
				}
				if (distance < closest_distance) {
					closest_distance = distance;
					closest_zone = i;
				}
				i++;
			}

			if (closest_zone != -1) {
				// move last here
				memcpy(&com.stacking_zones[closest_zone],
						&com.stacking_zones[i-1],
						sizeof(stacking_zone));
				com.stacking_zones[i-1].centre.x = -1.0;
			}
		}
	}
}
