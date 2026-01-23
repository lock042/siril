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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "histo_display.h"
#include <math.h>

/* Global histogram overlay state */
histo_overlay_state histo_state = {
	.show_histo = FALSE,
	.x = 20,
	.y = 0,  /* Will be calculated in init to be bottom-left */
	.width = 256,
	.height = 195,  /* Increased for stats panel */
	.min_width = 150,
	.min_height = 140,  /* Minimum to show stats properly */
	.max_width = 600,
	.max_height = 400,
	.opacity = 0.85,
	.rgb_area = NULL,
	.logarithmic = FALSE,  /* Start in linear mode */
	.is_dragging = FALSE,
	.is_resizing = FALSE,
	.resize_edge = RESIZE_NONE,
	.header_height = 20,
	.resize_border = 8
};

/* Cursor tracking for interactive histogram */
static struct {
	gboolean active;
	gint mouse_x;
	gint mouse_y;
} cursor_state = { .active = FALSE, .mouse_x = 0, .mouse_y = 0 };

/* Structure for histogram cache */
typedef struct {
	guint32 histo_r[UCHAR_MAX + 1];
	guint32 histo_g[UCHAR_MAX + 1];
	guint32 histo_b[UCHAR_MAX + 1];
	guint32 max_r, max_g, max_b;
	gboolean is_valid;
	guint64 image_checksum;

	/* Statistics */
	double mean_r, mean_g, mean_b;
	double median_r, median_g, median_b;
	double stddev_r, stddev_g, stddev_b;
} histo_cache;

static histo_cache cache = { .is_valid = FALSE, .image_checksum = 0 };

/* Compute histogram for a single channel */
static void compute_histogram_for_channel(int channel, guint32 *histo, guint32 *max_val) {
	int i;
	guint32 max = 0;

	for (i = 0; i < UCHAR_MAX + 1; i++) {
		histo[i] = 0;
	}

	if (!single_image_is_loaded()) {
		*max_val = 0;
		return;
	}

	if (gfit->type == DATA_USHORT) {
		WORD *buf = gfit->pdata[channel];
		guint32 nb_pixels = gfit->rx * gfit->ry;

		for (i = 0; i < nb_pixels; i++) {
			int bin = buf[i] >> 8;
			histo[bin]++;
			if (histo[bin] > max)
				max = histo[bin];
		}
	} else if (gfit->type == DATA_FLOAT) {
		float *buf = gfit->fpdata[channel];
		guint32 nb_pixels = gfit->rx * gfit->ry;

		for (i = 0; i < nb_pixels; i++) {
			int bin = (int)(buf[i] * UCHAR_MAX);
			if (bin < 0) bin = 0;
			if (bin > UCHAR_MAX) bin = UCHAR_MAX;
			histo[bin]++;
			if (histo[bin] > max)
				max = histo[bin];
		}
	}

	*max_val = max;
}

/* Compute image checksum */
static guint64 compute_image_checksum(void) {
	if (!single_image_is_loaded())
		return 0;

	guint64 checksum = 0;
	int step = (gfit->rx * gfit->ry) / 1000;
	if (step < 1) step = 1;

	if (gfit->type == DATA_USHORT) {
		for (int c = 0; c < gfit->naxes[2]; c++) {
			WORD *buf = gfit->pdata[c];
			for (int i = 0; i < gfit->rx * gfit->ry; i += step) {
				checksum = checksum * 31 + buf[i];
			}
		}
	} else if (gfit->type == DATA_FLOAT) {
		for (int c = 0; c < gfit->naxes[2]; c++) {
			float *buf = gfit->fpdata[c];
			for (int i = 0; i < gfit->rx * gfit->ry; i += step) {
				checksum = checksum * 31 + (guint64)(buf[i] * 65535.0f);
			}
		}
	}

	return checksum;
}

/* Update histogram cache if necessary */
static void update_histogram_cache(void) {
	if (!single_image_is_loaded()) {
		cache.is_valid = FALSE;
		return;
	}

	guint64 current_checksum = compute_image_checksum();

	if (cache.is_valid && cache.image_checksum == current_checksum)
		return;

	int nchannels = gfit->naxes[2];

	compute_histogram_for_channel(RLAYER, cache.histo_r, &cache.max_r);

	if (nchannels == 3) {
		compute_histogram_for_channel(GLAYER, cache.histo_g, &cache.max_g);
		compute_histogram_for_channel(BLAYER, cache.histo_b, &cache.max_b);
	} else {
		memcpy(cache.histo_g, cache.histo_r, sizeof(cache.histo_r));
		memcpy(cache.histo_b, cache.histo_r, sizeof(cache.histo_r));
		cache.max_g = cache.max_r;
		cache.max_b = cache.max_r;
	}

	cache.image_checksum = current_checksum;
	cache.is_valid = TRUE;

	/* Compute statistics */
	guint64 total_pixels = (guint64)gfit->rx * gfit->ry;

	/* Mean, median, stddev for each channel */
	for (int ch = 0; ch < 3; ch++) {
		guint32 *histo = (ch == 0) ? cache.histo_r : (ch == 1) ? cache.histo_g : cache.histo_b;
		double mean = 0.0;
		double variance = 0.0;
		double median = 0.0;

		/* Calculate mean */
		for (int i = 0; i <= UCHAR_MAX; i++) {
			mean += (i / 255.0) * histo[i];
		}
		mean /= total_pixels;

		/* Calculate variance */
		for (int i = 0; i <= UCHAR_MAX; i++) {
			double diff = (i / 255.0) - mean;
			variance += diff * diff * histo[i];
		}
		variance /= total_pixels;

		/* Calculate median (50th percentile) */
		guint64 cumulative = 0;
		guint64 half = total_pixels / 2;
		for (int i = 0; i <= UCHAR_MAX; i++) {
			cumulative += histo[i];
			if (cumulative >= half) {
				median = i / 255.0;
				break;
			}
		}

		/* Store results */
		if (ch == 0) {
			cache.mean_r = mean;
			cache.median_r = median;
			cache.stddev_r = sqrt(variance);
		} else if (ch == 1) {
			cache.mean_g = mean;
			cache.median_g = median;
			cache.stddev_g = sqrt(variance);
		} else {
			cache.mean_b = mean;
			cache.median_b = median;
			cache.stddev_b = sqrt(variance);
		}
	}
}

void invalidate_histogram_cache(void) {
	cache.is_valid = FALSE;
	cache.image_checksum = 0;
}

/* Draw tonal zones background (shadows, midtones, highlights) */
static void draw_tonal_zones(cairo_t *cr, double x, double y, double width, double height) {
	/* Shadows zone (0-85 / 0-33%) - dark background */
	cairo_rectangle(cr, x, y, width * 0.33, height);
	cairo_set_source_rgba(cr, 0.08, 0.08, 0.12, 0.3);
	cairo_fill(cr);

	/* Midtones zone (85-170 / 33%-66%) - medium background */
	cairo_rectangle(cr, x + width * 0.33, y, width * 0.33, height);
	cairo_set_source_rgba(cr, 0.12, 0.12, 0.15, 0.3);
	cairo_fill(cr);

	/* Highlights zone (170-255 / 66%-100%) - lighter background */
	cairo_rectangle(cr, x + width * 0.66, y, width * 0.34, height);
	cairo_set_source_rgba(cr, 0.15, 0.15, 0.18, 0.3);
	cairo_fill(cr);
}

/* Draw vertical grid lines */
static void draw_grid(cairo_t *cr, double x, double y, double width, double height) {
	cairo_set_source_rgba(cr, 0.3, 0.3, 0.3, 0.4);
	cairo_set_line_width(cr, 1.0);

	/* Vertical lines at 25%, 50%, 75% */
	for (int i = 1; i < 4; i++) {
		double grid_x = x + (width * i / 4.0);
		cairo_move_to(cr, grid_x, y);
		cairo_line_to(cr, grid_x, y + height);
	}
	cairo_stroke(cr);
}

/* Draw filled histogram with gradient for a channel */
static void draw_filled_channel(cairo_t *cr, double x, double y, double width, double height,
                                 guint32 *histo, guint32 max_val,
                                 double r, double g, double b) {
	if (max_val == 0)
		return;

	double display_values[UCHAR_MAX + 1];
	double display_max = 0.0;

	for (int i = 0; i <= UCHAR_MAX; i++) {
		double bin_val = (double)histo[i];

		/* Apply log to bin value if enabled (like Siril does) */
		if (histo_state.logarithmic && bin_val != 0.0) {
			bin_val = log(bin_val);  /* Natural log */
		}

		display_values[i] = bin_val;
		if (bin_val > display_max)
			display_max = bin_val;
	}

	if (display_max == 0.0)
		return;

	/* Create path for filled area */
	cairo_move_to(cr, x, y + height);

	for (int i = 0; i <= UCHAR_MAX; i++) {
		double bin_x = x + (i * width / (double)UCHAR_MAX);
		double normalized_height = display_values[i] / display_max;
		double bin_height = normalized_height * height;
		double bin_y = y + height - bin_height;
		cairo_line_to(cr, bin_x, bin_y);
	}

	cairo_line_to(cr, x + width, y + height);
	cairo_close_path(cr);

	/* Create gradient pattern - vertical gradient from transparent to color */
	cairo_pattern_t *gradient = cairo_pattern_create_linear(0, y, 0, y + height);
	cairo_pattern_add_color_stop_rgba(gradient, 0.0, r, g, b, 0.5);  /* Top: semi-transparent */
	cairo_pattern_add_color_stop_rgba(gradient, 1.0, r, g, b, 0.1);  /* Bottom: very transparent */

	cairo_set_source(cr, gradient);
	cairo_fill_preserve(cr);
	cairo_pattern_destroy(gradient);

	/* Draw outline */
	cairo_set_source_rgba(cr, r, g, b, 0.9);
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}

/* Draw statistics panel */
static void draw_statistics(cairo_t *cr, double x, double y, double width) {
	if (!cache.is_valid)
		return;

	cairo_set_source_rgba(cr, 0.15, 0.15, 0.15, 0.9);
	cairo_rectangle(cr, x, y, width, 25);
	cairo_fill(cr);

	/* Separator line */
	cairo_set_source_rgba(cr, 0.4, 0.4, 0.4, 0.8);
	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, x, y);
	cairo_line_to(cr, x + width, y);
	cairo_stroke(cr);

	cairo_set_source_rgba(cr, 0.8, 0.8, 0.8, 0.9);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr, 9.0);

	char stats_text[256];

	if (gfit->naxes[2] == 3) {
		/* RGB stats - show average across channels */
		double mean_avg = (cache.mean_r + cache.mean_g + cache.mean_b) / 3.0;
		double median_avg = (cache.median_r + cache.median_g + cache.median_b) / 3.0;
		double stddev_avg = (cache.stddev_r + cache.stddev_g + cache.stddev_b) / 3.0;

		snprintf(stats_text, sizeof(stats_text),
		         "μ:%.3f  Med:%.3f  σ:%.3f",
		         mean_avg, median_avg, stddev_avg);
	} else {
		/* Mono stats */
		snprintf(stats_text, sizeof(stats_text),
		         "μ:%.3f  Med:%.3f  σ:%.3f",
		         cache.mean_r, cache.median_r, cache.stddev_r);
	}

	cairo_move_to(cr, x + 5, y + 16);
	cairo_show_text(cr, stats_text);
}

/* Draw interactive cursor */
static void draw_interactive_cursor(cairo_t *cr, double x, double y, double width, double height) {
	if (!cursor_state.active || !cache.is_valid)
		return;

	/* Check if mouse is within histogram area */
	if (cursor_state.mouse_x < x || cursor_state.mouse_x > x + width ||
	    cursor_state.mouse_y < y || cursor_state.mouse_y > y + height)
		return;

	/* Calculate bin index from mouse position */
	double relative_x = cursor_state.mouse_x - x;
	int bin = (int)((relative_x / width) * UCHAR_MAX);
	if (bin < 0) bin = 0;
	if (bin > UCHAR_MAX) bin = UCHAR_MAX;

	/* Draw vertical cursor line */
	cairo_set_source_rgba(cr, 0.9, 0.9, 0.9, 0.7);
	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, cursor_state.mouse_x, y);
	cairo_line_to(cr, cursor_state.mouse_x, y + height);
	cairo_stroke(cr);

	/* Draw tooltip with value and pixel count */
	char tooltip[128];
	double value = bin / 255.0;

	if (gfit->naxes[2] == 3) {
		snprintf(tooltip, sizeof(tooltip), "%.2f: R:%u G:%u B:%u",
		         value, cache.histo_r[bin], cache.histo_g[bin], cache.histo_b[bin]);
	} else {
		snprintf(tooltip, sizeof(tooltip), "%.2f: %u px",
		         value, cache.histo_r[bin]);
	}

	/* Draw tooltip background */
	cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 0.9);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr, 9.0);

	cairo_text_extents_t extents;
	cairo_text_extents(cr, tooltip, &extents);

	/* Position tooltip above cursor, or below if near top */
	double tooltip_x = cursor_state.mouse_x - extents.width / 2.0;
	double tooltip_y = (cursor_state.mouse_y < y + 30) ? cursor_state.mouse_y + 20 : cursor_state.mouse_y - 10;

	/* Keep tooltip within bounds */
	if (tooltip_x < x + 2) tooltip_x = x + 2;
	if (tooltip_x + extents.width > x + width - 2) tooltip_x = x + width - extents.width - 2;

	cairo_rectangle(cr, tooltip_x - 3, tooltip_y - extents.height - 3,
	                extents.width + 6, extents.height + 6);
	cairo_fill(cr);

	/* Draw tooltip text */
	cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 0.95);
	cairo_move_to(cr, tooltip_x, tooltip_y);
	cairo_show_text(cr, tooltip);
}

/* Draw clipping indicators */
static void draw_clipping_indicators(cairo_t *cr, double x, double y, double width, double height) {
	if (!cache.is_valid)
		return;

	/* Count clipped pixels (first and last bins) */
	guint32 black_clipped = cache.histo_r[0] + cache.histo_g[0] + cache.histo_b[0];
	guint32 white_clipped = cache.histo_r[255] + cache.histo_g[255] + cache.histo_b[255];

	/* Threshold: show warning if more than 0.1% of pixels are clipped */
	guint32 total_pixels = gfit->rx * gfit->ry * gfit->naxes[2];
	guint32 clip_threshold = total_pixels / 1000;

	if (black_clipped > clip_threshold) {
		/* Left clipping indicator - red triangle */
		cairo_set_source_rgba(cr, 1.0, 0.2, 0.2, 0.8);
		cairo_move_to(cr, x, y);
		cairo_line_to(cr, x + 8, y);
		cairo_line_to(cr, x + 4, y + 6);
		cairo_close_path(cr);
		cairo_fill(cr);
	}

	if (white_clipped > clip_threshold) {
		/* Right clipping indicator - blue triangle */
		cairo_set_source_rgba(cr, 0.3, 0.5, 1.0, 0.8);
		cairo_move_to(cr, x + width, y);
		cairo_line_to(cr, x + width - 8, y);
		cairo_line_to(cr, x + width - 4, y + 6);
		cairo_close_path(cr);
		cairo_fill(cr);
	}
}

/* Determine which edge/area is under cursor */
static ResizeEdge get_resize_edge(gint mouse_x, gint mouse_y, gboolean *in_header) {
	gint x = histo_state.x;
	gint y = histo_state.y;
	gint w = histo_state.width;
	gint h = histo_state.height;
	gint border = histo_state.resize_border;
	ResizeEdge edge = RESIZE_NONE;

	/* Check if in header (for dragging) */
	if (mouse_x >= x && mouse_x <= x + w &&
	    mouse_y >= y && mouse_y <= y + histo_state.header_height) {
		*in_header = TRUE;
		return RESIZE_NONE;
	}
	*in_header = FALSE;

	/* Check resize borders */
	if (mouse_x >= x - border && mouse_x <= x + border)
		edge |= RESIZE_LEFT;
	else if (mouse_x >= x + w - border && mouse_x <= x + w + border)
		edge |= RESIZE_RIGHT;

	if (mouse_y >= y - border && mouse_y <= y + border)
		edge |= RESIZE_TOP;
	else if (mouse_y >= y + h - border && mouse_y <= y + h + border)
		edge |= RESIZE_BOTTOM;

	return edge;
}

/* Check if mouse is over Lin/Log toggle button */
static gboolean is_over_log_button(gint mouse_x, gint mouse_y) {
	gint button_width = 35;
	gint button_height = 16;
	gint button_x = histo_state.x + histo_state.width - button_width - 5;
	gint button_y = histo_state.y + 2;

	return (mouse_x >= button_x && mouse_x <= button_x + button_width &&
	        mouse_y >= button_y && mouse_y <= button_y + button_height);
}

/* Draw Lin/Log toggle button in header */
static void draw_log_toggle_button(cairo_t *cr, double x, double y, double width) {
	gint button_width = 35;
	gint button_height = 16;
	gint button_x = x + width - button_width - 5;
	gint button_y = y + 2;

	/* Draw button background */
	cairo_rectangle(cr, button_x, button_y, button_width, button_height);
	if (histo_state.logarithmic) {
		cairo_set_source_rgba(cr, 0.3, 0.5, 0.7, 0.9);  /* Highlighted when log */
	} else {
		cairo_set_source_rgba(cr, 0.25, 0.25, 0.25, 0.9);  /* Darker when linear */
	}
	cairo_fill_preserve(cr);

	/* Draw button border */
	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.9);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);

	/* Draw button text */
	cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 0.95);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_set_font_size(cr, 10.0);

	const char *text = histo_state.logarithmic ? "Log" : "Lin";
	cairo_text_extents_t extents;
	cairo_text_extents(cr, text, &extents);

	double text_x = button_x + (button_width - extents.width) / 2.0;
	double text_y = button_y + (button_height + extents.height) / 2.0;

	cairo_move_to(cr, text_x, text_y);
	cairo_show_text(cr, text);
}

/* Update cursor based on mouse position */
static void update_cursor(GtkWidget *widget, gint mouse_x, gint mouse_y) {
	GdkWindow *window = gtk_widget_get_window(widget);
	if (!window)
		return;

	GdkCursor *cursor = NULL;
	GdkDisplay *display = gdk_display_get_default();

	/* Check if over Lin/Log button first */
	if (is_over_log_button(mouse_x, mouse_y)) {
		cursor = gdk_cursor_new_for_display(display, GDK_HAND2);
		gdk_window_set_cursor(window, cursor);
		if (cursor)
			g_object_unref(cursor);
		return;
	}

	gboolean in_header = FALSE;
	ResizeEdge edge = get_resize_edge(mouse_x, mouse_y, &in_header);

	if (in_header) {
		cursor = gdk_cursor_new_for_display(display, GDK_FLEUR);
	} else if (edge == (RESIZE_LEFT | RESIZE_TOP) || edge == (RESIZE_RIGHT | RESIZE_BOTTOM)) {
		cursor = gdk_cursor_new_for_display(display, GDK_TOP_LEFT_CORNER);
	} else if (edge == (RESIZE_RIGHT | RESIZE_TOP) || edge == (RESIZE_LEFT | RESIZE_BOTTOM)) {
		cursor = gdk_cursor_new_for_display(display, GDK_TOP_RIGHT_CORNER);
	} else if (edge & RESIZE_LEFT || edge & RESIZE_RIGHT) {
		cursor = gdk_cursor_new_for_display(display, GDK_SB_H_DOUBLE_ARROW);
	} else if (edge & RESIZE_TOP || edge & RESIZE_BOTTOM) {
		cursor = gdk_cursor_new_for_display(display, GDK_SB_V_DOUBLE_ARROW);
	}

	gdk_window_set_cursor(window, cursor);
	if (cursor)
		g_object_unref(cursor);
}

/* Mouse button press event */
static gboolean on_histogram_button_press(GtkWidget *widget, GdkEventButton *event, gpointer data) {
	if (!histo_state.show_histo)
		return FALSE;

	if (event->button != 1) /* Left button only */
		return FALSE;

	/* Check if clicking on Lin/Log toggle button */
	if (is_over_log_button(event->x, event->y)) {
		histo_state.logarithmic = !histo_state.logarithmic;
		gtk_widget_queue_draw(widget);
		return TRUE;
	}

	gboolean in_header = FALSE;
	ResizeEdge edge = get_resize_edge(event->x, event->y, &in_header);

	if (in_header) {
		/* Start dragging */
		histo_state.is_dragging = TRUE;
		histo_state.drag_start_x = event->x;
		histo_state.drag_start_y = event->y;
		histo_state.drag_offset_x = event->x - histo_state.x;
		histo_state.drag_offset_y = event->y - histo_state.y;
		return TRUE;
	} else if (edge != RESIZE_NONE) {
		/* Start resizing */
		histo_state.is_resizing = TRUE;
		histo_state.resize_edge = edge;
		histo_state.resize_start_x = event->x;
		histo_state.resize_start_y = event->y;
		histo_state.resize_start_width = histo_state.width;
		histo_state.resize_start_height = histo_state.height;
		histo_state.resize_start_pos_x = histo_state.x;
		histo_state.resize_start_pos_y = histo_state.y;
		return TRUE;
	}

	return FALSE;
}

/* Mouse button release event */
static gboolean on_histogram_button_release(GtkWidget *widget, GdkEventButton *event, gpointer data) {
	if (event->button != 1)
		return FALSE;

	histo_state.is_dragging = FALSE;
	histo_state.is_resizing = FALSE;
	histo_state.resize_edge = RESIZE_NONE;

	return FALSE;
}

/* Mouse motion event */
static gboolean on_histogram_motion_notify(GtkWidget *widget, GdkEventMotion *event, gpointer data) {
	if (!histo_state.show_histo) {
		cursor_state.active = FALSE;
		return FALSE;
	}

	/* Update cursor position for interactive cursor */
	cursor_state.mouse_x = event->x;
	cursor_state.mouse_y = event->y;

	if (histo_state.is_dragging) {
		cursor_state.active = FALSE;  /* Disable cursor while dragging */

		/* Update position */
		histo_state.x = event->x - histo_state.drag_offset_x;
		histo_state.y = event->y - histo_state.drag_offset_y;

		/* Clamp to widget bounds */
		gint widget_width = gtk_widget_get_allocated_width(widget);
		gint widget_height = gtk_widget_get_allocated_height(widget);

		if (histo_state.x < 0)
			histo_state.x = 0;
		if (histo_state.y < 0)
			histo_state.y = 0;
		if (histo_state.x + histo_state.width > widget_width)
			histo_state.x = widget_width - histo_state.width;
		if (histo_state.y + histo_state.height > widget_height)
			histo_state.y = widget_height - histo_state.height;

		gtk_widget_queue_draw(widget);
		return TRUE;
	} else if (histo_state.is_resizing) {
		cursor_state.active = FALSE;  /* Disable cursor while resizing */

		/* Calculate size changes */
		gint dx = event->x - histo_state.resize_start_x;
		gint dy = event->y - histo_state.resize_start_y;

		gint new_width = histo_state.resize_start_width;
		gint new_height = histo_state.resize_start_height;
		gint new_x = histo_state.resize_start_pos_x;
		gint new_y = histo_state.resize_start_pos_y;

		if (histo_state.resize_edge & RESIZE_LEFT) {
			new_width -= dx;
			new_x += dx;
		} else if (histo_state.resize_edge & RESIZE_RIGHT) {
			new_width += dx;
		}

		if (histo_state.resize_edge & RESIZE_TOP) {
			new_height -= dy;
			new_y += dy;
		} else if (histo_state.resize_edge & RESIZE_BOTTOM) {
			new_height += dy;
		}

		/* Enforce size limits */
		if (new_width < histo_state.min_width) {
			if (histo_state.resize_edge & RESIZE_LEFT)
				new_x = histo_state.resize_start_pos_x + histo_state.resize_start_width - histo_state.min_width;
			new_width = histo_state.min_width;
		}
		if (new_width > histo_state.max_width) {
			if (histo_state.resize_edge & RESIZE_LEFT)
				new_x = histo_state.resize_start_pos_x + histo_state.resize_start_width - histo_state.max_width;
			new_width = histo_state.max_width;
		}

		if (new_height < histo_state.min_height) {
			if (histo_state.resize_edge & RESIZE_TOP)
				new_y = histo_state.resize_start_pos_y + histo_state.resize_start_height - histo_state.min_height;
			new_height = histo_state.min_height;
		}
		if (new_height > histo_state.max_height) {
			if (histo_state.resize_edge & RESIZE_TOP)
				new_y = histo_state.resize_start_pos_y + histo_state.resize_start_height - histo_state.max_height;
			new_height = histo_state.max_height;
		}

		histo_state.width = new_width;
		histo_state.height = new_height;
		histo_state.x = new_x;
		histo_state.y = new_y;

		gtk_widget_queue_draw(widget);
		return TRUE;
	} else {
		/* Not dragging or resizing - check if over histogram */
		double x = histo_state.x;
		double y = histo_state.y + histo_state.header_height;  /* Below header */
		double width = histo_state.width;
		double height = histo_state.height - histo_state.header_height - 25;  /* Above stats */

		if (event->x >= x && event->x <= x + width &&
		    event->y >= y && event->y <= y + height) {
			/* Mouse is over histogram area - activate cursor */
			cursor_state.active = TRUE;
			gtk_widget_queue_draw(widget);
		} else {
			/* Mouse outside histogram - check for cursor change */
			if (cursor_state.active) {
				cursor_state.active = FALSE;
				gtk_widget_queue_draw(widget);
			}

			/* Update cursor based on position */
			update_cursor(widget, event->x, event->y);
		}
	}

	return FALSE;
}

/* Draw histogram overlay */
static gboolean on_histogram_overlay_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	if (!histo_state.show_histo || !single_image_is_loaded())
		return FALSE;

	update_histogram_cache();

	if (!cache.is_valid)
		return FALSE;

	/* Auto-position to bottom-left on first draw if widget wasn't allocated during init */
	static gboolean positioned = FALSE;
	if (!positioned) {
		gint widget_height = gtk_widget_get_allocated_height(widget);
		if (widget_height > histo_state.height + 40) {
			histo_state.x = 20;
			histo_state.y = widget_height - histo_state.height - 20;
			positioned = TRUE;
		}
	}

	double x = histo_state.x;
	double y = histo_state.y;
	double width = histo_state.width;
	double height = histo_state.height;
	double header_height = histo_state.header_height;

	cairo_save(cr);

	/* Draw outer frame with rounded corners */
	double corner_radius = 8.0;
	cairo_new_sub_path(cr);
	cairo_arc(cr, x + corner_radius, y + corner_radius, corner_radius, M_PI, 3 * M_PI / 2);
	cairo_arc(cr, x + width - corner_radius, y + corner_radius, corner_radius, 3 * M_PI / 2, 2 * M_PI);
	cairo_arc(cr, x + width - corner_radius, y + height - corner_radius, corner_radius, 0, M_PI / 2);
	cairo_arc(cr, x + corner_radius, y + height - corner_radius, corner_radius, M_PI / 2, M_PI);
	cairo_close_path(cr);

	cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, histo_state.opacity);
	cairo_fill_preserve(cr);

	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, histo_state.opacity);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);

	/* Draw header (drag area) */
	cairo_rectangle(cr, x, y, width, header_height);
	cairo_set_source_rgba(cr, 0.2, 0.2, 0.2, histo_state.opacity);
	cairo_fill(cr);

	/* Draw header separator line */
	cairo_move_to(cr, x, y + header_height);
	cairo_line_to(cr, x + width, y + header_height);
	cairo_set_source_rgba(cr, 0.4, 0.4, 0.4, histo_state.opacity);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);

	/* Draw "grip" dots in header */
	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, histo_state.opacity);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			double dot_x = x + width/2 - 6 + i * 6;
			double dot_y = y + header_height/2 - 2 + j * 4;
			cairo_arc(cr, dot_x, dot_y, 1.5, 0, 2 * M_PI);
			cairo_fill(cr);
		}
	}

	/* Draw Lin/Log toggle button */
	draw_log_toggle_button(cr, x, y, width);

	/* Define histogram content area (below header, with space for stats panel) */
	double histo_y = y + header_height;
	double stats_height = 25;  /* Height of statistics panel */
	double histo_height = height - header_height - stats_height;

	/* Clip to histogram content area for cleaner drawing */
	cairo_rectangle(cr, x, histo_y, width, histo_height);
	cairo_clip(cr);

	/* Draw tonal zones background */
	draw_tonal_zones(cr, x, histo_y, width, histo_height);

	/* Draw vertical grid */
	draw_grid(cr, x, histo_y, width, histo_height);

	/* Draw filled histogram curves */
	if (gfit->naxes[2] == 3) {
		/* Draw filled areas for each channel */
		draw_filled_channel(cr, x, histo_y, width, histo_height,
		                    cache.histo_b, cache.max_b, 0.3, 0.5, 1.0);  /* Blue */
		draw_filled_channel(cr, x, histo_y, width, histo_height,
		                    cache.histo_g, cache.max_g, 0.0, 1.0, 0.0);  /* Green */
		draw_filled_channel(cr, x, histo_y, width, histo_height,
		                    cache.histo_r, cache.max_r, 1.0, 0.0, 0.0);  /* Red */
	} else {
		/* Mono - white/gray */
		draw_filled_channel(cr, x, histo_y, width, histo_height,
		                    cache.histo_r, cache.max_r, 0.8, 0.8, 0.8);
	}

	/* Reset clip */
	cairo_reset_clip(cr);

	/* Draw statistics panel at bottom */
	draw_statistics(cr, x, y + height - stats_height, width);

	/* Draw clipping indicators */
	draw_clipping_indicators(cr, x, histo_y, width, histo_height);

	/* Draw interactive cursor if active */
	draw_interactive_cursor(cr, x, histo_y, width, histo_height);

	cairo_restore(cr);

	return FALSE;
}

/* Initialize histogram overlay */
void init_histogram_overlay(void) {
	siril_log_color_message("Initializing histogram overlay...\n", "green");

	histo_state.rgb_area = lookup_widget("drawingareargb");
	if (!histo_state.rgb_area) {
		siril_log_color_message("ERROR: Cannot find drawingareargb widget!\n", "red");
		return;
	}

	/* Calculate initial position at bottom-left */
	gint widget_height = gtk_widget_get_allocated_height(histo_state.rgb_area);
	if (widget_height > 0) {
		/* Position at bottom-left with margin */
		histo_state.x = 20;
		histo_state.y = widget_height - histo_state.height - 20;
		siril_log_message("Initial histogram position: %dx%d at (%d,%d)\n",
		                  histo_state.width, histo_state.height,
		                  histo_state.x, histo_state.y);
	} else {
		/* Widget not yet allocated, use default that will be adjusted on first draw */
		histo_state.x = 20;
		histo_state.y = 300;
	}

	/* Connect draw signal */
	g_signal_connect_after(histo_state.rgb_area, "draw",
	                        G_CALLBACK(on_histogram_overlay_draw), NULL);

	/* Connect mouse events for drag and resize */
	gtk_widget_add_events(histo_state.rgb_area,
	                      GDK_BUTTON_PRESS_MASK |
	                      GDK_BUTTON_RELEASE_MASK |
	                      GDK_POINTER_MOTION_MASK);

	g_signal_connect(histo_state.rgb_area, "button-press-event",
	                 G_CALLBACK(on_histogram_button_press), NULL);
	g_signal_connect(histo_state.rgb_area, "button-release-event",
	                 G_CALLBACK(on_histogram_button_release), NULL);
	g_signal_connect(histo_state.rgb_area, "motion-notify-event",
	                 G_CALLBACK(on_histogram_motion_notify), NULL);

	siril_log_color_message("Histogram overlay initialized successfully!\n", "green");
}

void set_histogram_overlay_visible(gboolean visible) {
	histo_state.show_histo = visible;

	if (!histo_state.rgb_area) {
		siril_log_color_message("ERROR: rgb_area is NULL!\n", "red");
		return;
	}

	gtk_widget_queue_draw(histo_state.rgb_area);
}

void update_histogram_display(void) {
	if (!histo_state.rgb_area || !histo_state.show_histo)
		return;

	invalidate_histogram_cache();
	update_histogram_cache();
	gtk_widget_queue_draw(histo_state.rgb_area);
}