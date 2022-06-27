/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/histogram.h"

#include "global_histogram.h"

#include <math.h>

#define PU(V, params) (MIN((V), (params->bins_count - 1)))

static uint32_t *hist = NULL;

typedef enum global_histogram_scale_t {
	GLOBAL_HISTOGRAM_SCALE_LOGARITHMIC = 0,
	GLOBAL_HISTOGRAM_SCALE_LINEAR,
	GLOBAL_HISTOGRAM_SCALE_N
} global_histogram_scale_t;

typedef struct global_histogram_t {
	// histogram for display
	uint32_t *histogram;
	uint32_t histogram_max;

	global_histogram_scale_t histogram_scale;

} global_histogram_t;

static float histogram_max = 0;

GdkRGBA histo_color[] = {
		{ 1.0, 0.0, 0.0, 1.0 },
		{ 0.0, 1.0, 0.0, 1.0 },
		{ 0.0, 0.0, 1.0, 1.0 }
};

static void global_draw_line(cairo_t *cr, float left, float top, float right, float bottom) {
	cairo_move_to(cr, left, top);
	cairo_line_to(cr, right, bottom);
}

static void global_draw_grid(cairo_t *cr, const int num, const int left, const int top, const int right, const int bottom) {
	float width = right - left;
	float height = bottom - top;

	for (int k = 1; k < num; k++) {
		global_draw_line(cr, left + k / (float) num * width, top, left + k / (float) num * width, bottom);
		cairo_stroke(cr);
		global_draw_line(cr, left, top + k / (float) num * height, right, top + k / (float) num * height);
		cairo_stroke(cr);
	}
}

static void global_draw_histogram_16_linxliny(cairo_t *cr, const uint32_t *hist, int32_t channels, int32_t channel) {
	cairo_move_to(cr, 0, 0);
	for (int k = 0; k < HISTOGRAM_BINS; k++)
		cairo_line_to(cr, k, hist[channels * k + channel]);
	cairo_line_to(cr, 65535, 0);
	cairo_close_path(cr);
	cairo_fill(cr);
}

// linear x log y
static void global_draw_histogram_16_linxlogy(cairo_t *cr, const uint32_t *hist, int32_t channels, int32_t channel) {
	cairo_move_to(cr, 0, 0);
	for (int k = 0; k < HISTOGRAM_BINS; k++)
		cairo_line_to(cr, k, logf(1.0 + hist[channels * k + channel]));
	cairo_line_to(cr, 65535, 0);
	cairo_close_path(cr);
	cairo_fill(cr);
}

// linear x
static void global_draw_histogram_16(cairo_t *cr, const uint32_t *hist, int32_t channels, int32_t channel, const gboolean linear) {
	if (linear) // linear y
		global_draw_histogram_16_linxliny(cr, hist, channels, channel);
	else
		// log y
		global_draw_histogram_16_linxlogy(cr, hist, channels, channel);
}

static uint32_t *compute_histo() {
	uint32_t *hist = calloc(3 * HISTOGRAM_BINS, sizeof(uint32_t));

	for (int c = 0; c < gfit.naxes[2]; c++) {
		gsl_histogram *histo = computeHisto(&gfit, c);

		double tmp = gsl_histogram_max_val(histo);
		histogram_max = MAX(histogram_max, tmp);

		for (size_t i = 0; i < USHRT_MAX + 1; i++) {
			uint32_t bin = (int) gsl_histogram_get (histo, i);
			hist[i * 3 + c] = bin;
		}
	}

	return hist;
}

static void global_histogram_draw_histogram(/*dt_lib_histogram_t *d, */cairo_t *cr,
		int width, int height) {
//	if (!d->histogram_max)
//		return;

	const float hist_max = histogram_max;
		//	d->histogram_scale == DT_LIB_HISTOGRAM_SCALE_LINEAR ?
			//		d->histogram_max : logf(1.0 + d->histogram_max);
	cairo_save(cr);
	cairo_push_group_with_content(cr, CAIRO_CONTENT_COLOR);
	cairo_translate(cr, 0, height);
	cairo_scale(cr, width / 65535.0, -(height - 10) / hist_max);
	cairo_set_operator(cr, CAIRO_OPERATOR_ADD);
	cairo_set_line_width(cr, 1.);
	for (int k = 0; k < 3; k++) {
		cairo_set_source_rgba(cr, histo_color[k].red, histo_color[k].green, histo_color[k].blue, histo_color[k].alpha);
		global_draw_histogram_16(cr, hist, 3, k, TRUE);
	}
	cairo_pop_group_to_source(cr);
	cairo_set_operator(cr, CAIRO_OPERATOR_ADD);
	cairo_paint_with_alpha(cr, 0.5);
	cairo_restore(cr);
}

gboolean on_drawingarea_global_histograms_draw(GtkWidget *widget, cairo_t *crf, gpointer user_data) {
	global_histogram_t *d = (global_histogram_t*) user_data;

	GtkAllocation allocation;
	gtk_widget_get_allocation(widget, &allocation);
	const int width = allocation.width, height = allocation.height;

	cairo_surface_t *cst = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
			width, height);
	cairo_t *cr = cairo_create(cst);

	gtk_render_background(gtk_widget_get_style_context(widget), cr, 0, 0, width,
			height);
	cairo_set_line_width(cr, 0.5); // borders width

	// Draw frame and background
	cairo_save(cr);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 1.0);
	cairo_fill(cr);
	cairo_restore(cr);

	// draw grid
	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 1.0);
	global_draw_grid(cr, 4, 0, 0, width, height);

	g_mutex_lock(&com.mutex);

	global_histogram_draw_histogram(cr, width, height);

	g_mutex_unlock(&com.mutex);

	// finally a thin border
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_set_source_rgba(cr, 0.9, 0.9, 0.9, 1.0);
	cairo_stroke(cr);

	cairo_destroy(cr);
	cairo_set_source_surface(crf, cst, 0, 0);
	cairo_paint(crf);
	cairo_surface_destroy(cst);

	return TRUE;
}

void compute_global_histo() {
	if (hist == NULL) {
		hist = compute_histo(&gfit);
	}
}

void free_global_histo() {
	free(hist);
	hist = NULL;
}

void gui_global_histo_init() {
	/* initialize ui widgets */
	global_histogram_t *d = malloc(sizeof(global_histogram_t));

}
