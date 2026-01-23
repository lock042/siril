/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2025 Free Software Foundation, Inc.
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

#ifndef _HISTO_DISPLAY_H
#define _HISTO_DISPLAY_H

#include <gtk/gtk.h>

/* Resize edge flags */
typedef enum {
	RESIZE_NONE = 0,
	RESIZE_LEFT = 1 << 0,
	RESIZE_RIGHT = 1 << 1,
	RESIZE_TOP = 1 << 2,
	RESIZE_BOTTOM = 1 << 3
} ResizeEdge;

/* Structure to store histogram overlay display state */
typedef struct {
	gboolean show_histo;      /* Toggle state */
	gint x;                   /* Position X */
	gint y;                   /* Position Y */
	gint width;               /* Histogram width in pixels */
	gint height;              /* Histogram height in pixels */
	gint min_width;           /* Minimum width */
	gint min_height;          /* Minimum height */
	gint max_width;           /* Maximum width */
	gint max_height;          /* Maximum height */
	gdouble opacity;          /* Overlay opacity (0.0-1.0) */
	GtkWidget *rgb_area;      /* Reference to RGB drawing area */

	/* Display mode */
	gboolean logarithmic;     /* TRUE = log scale, FALSE = linear */

	/* Channel visibility */
	gboolean show_red;        /* Show red channel */
	gboolean show_green;      /* Show green channel */
	gboolean show_blue;       /* Show blue channel */

	/* Drag state */
	gboolean is_dragging;
	gint drag_start_x;
	gint drag_start_y;
	gint drag_offset_x;
	gint drag_offset_y;

	/* Resize state */
	gboolean is_resizing;
	ResizeEdge resize_edge;
	gint resize_start_x;
	gint resize_start_y;
	gint resize_start_width;
	gint resize_start_height;
	gint resize_start_pos_x;
	gint resize_start_pos_y;

	/* UI elements */
	gint header_height;       /* Height of drag header */
	gint resize_border;       /* Width of resize border area */
} histo_overlay_state;

/* Global overlay state */
extern histo_overlay_state histo_state;

/* Initialize histogram overlay */
void init_histogram_overlay(void);

/* Show or hide the histogram overlay */
void set_histogram_overlay_visible(gboolean visible);

/* Invalidate histogram cache */
void invalidate_histogram_cache(void);

/* Update histogram display */
void update_histogram_display(void);

#endif /* _HISTO_DISPLAY_H */