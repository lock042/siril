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

/* This file manages small selection areas called zones, used in multi-point
 * planetary processing and in background extraction.
 * Zones are stored in com.stacking_zones. */

#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "gui.h"
#include "zones.h"
#include "image_display.h"
#include "planetary_callbacks.h"

static gboolean add_zones_mode = FALSE;
static gboolean remove_zones_mode = FALSE;

int get_number_of_zones() {
	int i = 0;
	if (!com.stacking_zones)
		return 0;
	while (com.stacking_zones[i].centre.x >= 0.0) i++;
	return i;
}

int get_side(const stacking_zone *zone) {
	return round_to_int(zone->half_side * 2.0);
}

int point_is_inside_zone(int px, int py, const stacking_zone *zone) {
	int side = round_to_int(zone->half_side * 2.0);
	int startx = round_to_int(zone->centre.x - zone->half_side);
	int starty = round_to_int(zone->centre.y - zone->half_side);
	return px > startx && px < startx + side && py > starty && py < starty + side;
}

/***************** GUI stuff ********************/

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

void remove_stacking_zones() {
	if (com.stacking_zones) {
		int i = 0;
		while (com.stacking_zones[i].centre.x >= 0.0) {
			com.stacking_zones[i].centre.x = -1.0;
			i++;
		}
	}
}

/* add a stacking zone centred around x,y with the given half side */
void add_stacking_zone(double x, double y, double half_side) {
	if (!com.stacking_zones_size || !com.stacking_zones) {
		com.stacking_zones_size = 40;
		com.stacking_zones = malloc(com.stacking_zones_size * sizeof(stacking_zone));
		com.stacking_zones[0].centre.x = -1.0;
		com.stacking_zones[0].mpregparam = NULL;
	}

	int i = 0;
	while (com.stacking_zones && com.stacking_zones[i].centre.x >= 0.0 && i < com.stacking_zones_size - 1)
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
	zone->mpregparam = NULL;
	fprintf(stdout, "Added stacking zone %d at %4g,%4g. Half side: %g\n", i, x, y, half_side);

	com.stacking_zones[i+1].centre.x = -1.0;
	com.stacking_zones[i+1].mpregparam = NULL;

	activate_mpp_processing_button();
}

gboolean on_remove_all_zones_clicked(GtkButton *button, gpointer user_data) {
	remove_stacking_zones();
	update_zones_list();
	redraw(com.cvport, REMAP_NONE);
	return FALSE;
}

void planetary_click_in_image(double x, double y) {
	int i = 0;
	if (add_zones_mode) {
		GtkAdjustment *sizeadj = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_zonesize"));
		double size = gtk_adjustment_get_value(sizeadj);
		if (size > 0.0)
			add_stacking_zone(x, y, size*0.5);
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
				if (com.stacking_zones[closest_zone].mpregparam) {
					free(com.stacking_zones[closest_zone].mpregparam);
					com.stacking_zones[closest_zone].mpregparam = NULL;
				}
				// replace by last zone
				memcpy(&com.stacking_zones[closest_zone],
						&com.stacking_zones[i-1],
						sizeof(stacking_zone));
				com.stacking_zones[i-1].centre.x = -1.0;

				activate_mpp_processing_button();
				update_zones_list();
			}
		}
	}
}

/* more than 50% overlapping */
gboolean zone_is_too_close(int x, int y, int size) {
	if (com.stacking_zones) {
		int i = 0;
		while (com.stacking_zones[i].centre.x >= 0.0) {
			if ((fabs(com.stacking_zones[i].centre.x - x) < size)
					&& (fabs(com.stacking_zones[i].centre.y - y) < size))
				return TRUE;
			i++;
		}
	}
	return FALSE;
}

