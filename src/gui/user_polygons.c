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

#include "core/siril.h"
#include "gui/user_polygons.h"
#include "gui/image_display.h"

#define FROM_BE64_INTO(dest, val, type) \
do { \
	union { type v; uint64_t i; } conv; \
	memcpy(&conv.i, &val, sizeof(type)); \
	conv.i = GUINT64_FROM_BE(conv.i); \
	(dest) = conv.v; \
} while(0)

#define MAX_LEGEND_LENGTH 4095  // Allows for +1 null terminator without exceeding 4096

static gint unused_polygon_id = 0;

UserPolygon *find_polygon_by_id(int id) {
	GSList *node;

	// Special case: if id is -1, find most recently created polygon
	if (id == -1) {
		int highest_id = g_atomic_int_get(&unused_polygon_id) - 1;
		if (highest_id < 0) return NULL; // No polygons exist

		for (node = gui.user_polygons; node != NULL; node = node->next) {
			UserPolygon *polygon = (UserPolygon *)node->data;
			if (polygon->id == highest_id) {
				return polygon;
			}
		}
		return NULL; // Polygon with highest ID was deleted
	}

	// Normal case: find polygon with matching id
	for (node = gui.user_polygons; node != NULL; node = node->next) {
		UserPolygon *polygon = (UserPolygon *)node->data;
		if (polygon->id == id) {
			return polygon;
		}
	}
	return NULL;
}

int get_unused_polygon_id(void) {
	return g_atomic_int_add(&unused_polygon_id, 1);
}

UserPolygon* create_user_polygon_from_points(GSList *point_list) {
	if (!point_list) {
		return NULL;  // Handle empty list
	}

	// Allocate memory for the UserPolygon
	UserPolygon *polygon = g_malloc(sizeof(UserPolygon));
	if (!polygon) {
		return NULL;  // Memory allocation failed
	}

	// Count the number of points in the list
	int n_points = g_slist_length(point_list);

	// Allocate memory for the points array
	polygon->points = g_malloc(n_points * sizeof(point));
	if (!polygon->points) {
		g_free(polygon);
		return NULL;  // Memory allocation failed
	}

	// Copy points from GSList to the array
	GSList *current = point_list;
	for (int i = 0; i < n_points; i++) {
		polygon->points[i] = *(point*)current->data;
		current = current->next;
	}

	// Initialize the struct fields
	polygon->id = -1;
	polygon->n_points = n_points;

	// Set color to 0xFFFFFF40 (green with alpha)
	polygon->color.red = 0.0;
	polygon->color.green = 1.0;
	polygon->color.blue = 0.0;
	polygon->color.alpha = 0.25;

	polygon->fill = TRUE;
	polygon->legend = NULL;

	return polygon;
}

int add_user_polygon(point *points, int n_points, const GdkRGBA *color, gboolean fill) {
	int id = get_unused_polygon_id();

	UserPolygon *polygon = g_new(UserPolygon, 1);
	polygon->id = id;
	polygon->n_points = n_points;
	polygon->points = g_new(point, n_points);
	polygon->fill = fill;
	polygon->color = *color;

	for (int i = 0; i < n_points; i++) {
		polygon->points[i] = points[i];
	}

	gui.user_polygons = g_slist_prepend(gui.user_polygons, polygon);

	return id; // Return the generated ID to the caller
}

int add_existing_polygon(UserPolygon* polygon, const GdkRGBA *color, gboolean fill) {
	int id = get_unused_polygon_id();
	polygon->id = id;
	polygon->fill = fill;
	if (color) {
		polygon->color = *color;
	}
	gui.user_polygons = g_slist_prepend(gui.user_polygons, polygon);
	return id; // Return the generated ID to the caller
}


void free_user_polygon(gpointer data) {
	UserPolygon *polygon = (UserPolygon *)data;
	if (polygon) {
		g_free(polygon->points);
		g_free(polygon->legend);
		g_free(polygon);
	}
}

gboolean delete_user_polygon(int id) {
	GSList *l;
	for (l = gui.user_polygons; l != NULL; l = l->next) {
		UserPolygon *polygon = (UserPolygon *)l->data;
		if (polygon->id == id) {
			gui.user_polygons = g_slist_delete_link(gui.user_polygons, l);
			free_user_polygon(polygon);
			return TRUE;
		}
	}
	return FALSE;
}

void clear_user_polygons(void) {
	g_slist_free_full(gui.user_polygons, free_user_polygon);
	gui.user_polygons = NULL;
	g_atomic_int_set(&unused_polygon_id, 0);
	queue_redraw(REDRAW_OVERLAY);
}

UserPolygon* deserialize_polygon(const uint8_t *data, size_t size) {
	if (size < 17) {  // 4 + 4 + 4 + 1 + 4 bytes (minimum with string length)
		fprintf(stderr, "Invalid data size\n");
		return NULL;
	}
	const uint8_t *ptr = data;
	UserPolygon *polygon = g_malloc0(sizeof(UserPolygon));
	if (!polygon) {
		fprintf(stderr, "Memory allocation failed for polygon\n");
		return NULL;
	}

	// Read ID
	polygon->id = (int32_t)g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	// Read number of points
	polygon->n_points = (int32_t)g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	// Read RGBA (packed into uint32_t)
	uint32_t packed_color = g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);
	polygon->color.red   = ((packed_color >> 24) & 0xFF) / 255.0;
	polygon->color.green = ((packed_color >> 16) & 0xFF) / 255.0;
	polygon->color.blue  = ((packed_color >> 8)  & 0xFF) / 255.0;
	polygon->color.alpha = (packed_color & 0xFF) / 255.0;

	// Read fill flag (1 byte boolean)
	polygon->fill = (gboolean)(*ptr);
	ptr += 1;  // Move pointer past boolean

	if (polygon->n_points <= 0 || polygon->n_points > MAX_POLYGON_POINTS) {
		fprintf(stderr, "Invalid number of points: %d\n", polygon->n_points);
		g_free(polygon);
		return NULL;
	}

	// Allocate memory for points
	polygon->points = g_malloc0(polygon->n_points * sizeof(point));
	if (!polygon->points) {
		fprintf(stderr, "Memory allocation for points failed\n");
		g_free(polygon);
		return NULL;
	}

	// Ensure enough data for points
	size_t points_size = polygon->n_points * (2 * sizeof(double));
	size_t required_size = 13 + points_size + 4;  // Basic data + points + legend length field
	if (size < required_size) {
		fprintf(stderr, "Not enough data for points and legend length\n");
		g_free(polygon->points);
		g_free(polygon);
		return NULL;
	}

	// Read points
	for (int i = 0; i < polygon->n_points; i++) {
		double x_BE, y_BE;
		memcpy(&x_BE, ptr, sizeof(double));
		ptr += sizeof(double);
		FROM_BE64_INTO(polygon->points[i].x, x_BE, double);
		memcpy(&y_BE, ptr, sizeof(double));
		ptr += sizeof(double);
		FROM_BE64_INTO(polygon->points[i].y, y_BE, double);
	}

	// Read legend string length
	int32_t legend_length = (int32_t)g_ntohl(*(uint32_t *)ptr);
	ptr += sizeof(uint32_t);

	// Check if we have enough data for the legend string
	if (size < required_size + legend_length) {
		fprintf(stderr, "Not enough data for legend string\n");
		g_free(polygon->points);
		g_free(polygon);
		return NULL;
	}

	// Handle the legend string
	if (legend_length > 0) {
		if (legend_length > MAX_LEGEND_LENGTH) {
			fprintf(stderr, "Legend length exceeds maximum allowed\n");
			g_free(polygon->points);
			g_free(polygon);
			return NULL;
		}
		polygon->legend = g_malloc0(legend_length + 1);  // +1 for null terminator
		if (!polygon->legend) {
			fprintf(stderr, "Memory allocation for legend failed\n");
			g_free(polygon->points);
			g_free(polygon);
			return NULL;
		}
		memcpy(polygon->legend, ptr, legend_length);
		polygon->legend[legend_length] = '\0';  // Ensure null termination
	} else {
		polygon->legend = NULL;
	}

	return polygon;
}

uint8_t* serialize_polygon(UserPolygon *polygon, size_t *size) {
	if (!polygon || !size) {
		fprintf(stderr, "Invalid arguments\n");
		return NULL;
	}

	// Compute the required buffer size
	size_t points_size = polygon->n_points * (2 * sizeof(double));
	size_t legend_size = polygon->legend ? strlen(polygon->legend) : 0;
	*size = 4 + 4 + 4 + 1 + points_size + 4 + legend_size; // id + n_points + color + fill + points + legend_len + legend

	uint8_t *buffer = g_malloc0(*size);
	if (!buffer) {
		fprintf(stderr, "Memory allocation for serialization failed\n");
		return NULL;
	}

	uint8_t *ptr = buffer;

	// Write ID
	*(uint32_t *)ptr = g_htonl((uint32_t)polygon->id);
	ptr += sizeof(uint32_t);

	// Write number of points
	*(uint32_t *)ptr = g_htonl((uint32_t)polygon->n_points);
	ptr += sizeof(uint32_t);

	// Write color as packed RGBA
	uint32_t packed_color =
	((uint32_t)(polygon->color.red * 255) << 24) |
	((uint32_t)(polygon->color.green * 255) << 16) |
	((uint32_t)(polygon->color.blue * 255) << 8) |
	((uint32_t)(polygon->color.alpha * 255));
	*(uint32_t *)ptr = g_htonl(packed_color);
	ptr += sizeof(uint32_t);

	// Write fill flag
	*ptr = (uint8_t)polygon->fill;
	ptr += 1;

	// Write points
	for (int i = 0; i < polygon->n_points; i++) {
		double x_BE, y_BE;
		FROM_BE64_INTO(x_BE, polygon->points[i].x, double);
		FROM_BE64_INTO(y_BE, polygon->points[i].y, double);
		memcpy(ptr, &x_BE, sizeof(double));
		ptr += sizeof(double);
		memcpy(ptr, &y_BE, sizeof(double));
		ptr += sizeof(double);
	}

	// Write legend length and string
	*(uint32_t *)ptr = g_htonl((uint32_t)legend_size);
	ptr += sizeof(uint32_t);

	if (legend_size > 0) {
		memcpy(ptr, polygon->legend, legend_size);
	}

	return buffer;
}

uint8_t* serialize_polygon_list(GSList *polygons, size_t *out_size) {
	if (!polygons) {
		*out_size = 0;
		return NULL;
	}

	size_t total_size = sizeof(uint32_t); // Space for number of polygons
	GSList *node;
	for (node = polygons; node != NULL; node = node->next) {
		UserPolygon *polygon = (UserPolygon *)node->data;
		size_t polygon_size;
		uint8_t *serialized = serialize_polygon(polygon, &polygon_size);
		if (!serialized) {
			continue;
		}
		total_size += polygon_size;
		g_free(serialized);
	}

	uint8_t *buffer = g_malloc0(total_size);
	if (!buffer) {
		fprintf(stderr, "Failed to allocate memory for serialized polygon list\n");
		*out_size = 0;
		return NULL;
	}

	uint8_t *ptr = buffer;

	// Write number of polygons
	*(uint32_t *)ptr = g_htonl(g_slist_length(polygons));
	ptr += sizeof(uint32_t);

	for (node = polygons; node != NULL; node = node->next) {
		UserPolygon *polygon = (UserPolygon *)node->data;
		size_t polygon_size;
		uint8_t *serialized = serialize_polygon(polygon, &polygon_size);
		if (!serialized) {
			g_free(buffer);
			return NULL;
		}
		memcpy(ptr, serialized, polygon_size);
		ptr += polygon_size;
		g_free(serialized);
	}

	*out_size = total_size;
	return buffer;
}
