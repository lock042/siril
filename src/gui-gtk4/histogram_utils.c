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

#include <gsl/gsl_histogram.h>
#include <string.h>
#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "histogram_utils.h"

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// Colors of layer histograms: R, G, B, RGB (last entry unused for display_histo)
static const double histo_color_r[] = { 1.0, 0.0, 0.0, 0.0 };
static const double histo_color_g[] = { 0.0, 1.0, 0.0, 0.0 };
static const double histo_color_b[] = { 0.0, 0.0, 1.0, 0.0 };

size_t get_histo_size(fits *fit) {
	if (fit->type == DATA_USHORT) {
		if (fit->orig_bitpix == BYTE_IMG)
			return UCHAR_MAX;
	}
	return (size_t)USHRT_MAX;
}

// Create a new histogram object for the passed fit and layer
gsl_histogram *computeHisto(fits *fit, int layer) {
	g_assert(layer < 3);
	size_t i, ndata, size;

	size = get_histo_size(fit);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);
	ndata = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		gsl_histogram *histo_thr = gsl_histogram_alloc(size + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);

		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				if (buf[i] == 0)
					continue;
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				if (buf[i] == 0.f)
					continue;
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			gsl_histogram_add(histo, histo_thr);
		}
		gsl_histogram_free(histo_thr);
	}

	return histo;
}

gsl_histogram *computeHisto_Selection(fits *fit, int layer, rectangle *selection) {
	g_assert(layer < 3);

	size_t size = get_histo_size(fit);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 : size);
	size_t stridefrom = fit->rx - selection->w;

	if (fit->type == DATA_USHORT) {
		WORD *from = fit->pdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	} else if (fit->type == DATA_FLOAT) {
		float *from = fit->fpdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	return histo;
}

void set_histogram(gsl_histogram *histo, int layer) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

void set_sat_histogram(gsl_histogram *histo) {
	if (com.sat_hist)
		gsl_histogram_free(com.sat_hist);
	com.sat_hist = histo;
}

/* call from main thread */
void compute_histo_for_fit(fits *thefit) {
	int nb_layers = 3;
	if (thefit->naxis == 2)
		nb_layers = 1;
	for (int i = 0; i < nb_layers; i++) {
		if (!com.layers_hist[i])
			set_histogram(computeHisto(thefit, i), i);
	}
}

/* call from any thread */
void invalidate_gfit_histogram() {
	for (int layer = 0; layer < MAXVPORT; layer++) {
		set_histogram(NULL, layer);
	}
	set_sat_histogram(NULL);
}

void fill_histo_background(cairo_t *cr, int width, int height) {
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);
}

void draw_grid(cairo_t *cr, int width, int height) {
	double dash_format[] = { 1.0, 1.0 };

	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
	// quarters in solid, eighths in dashed line
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_move_to(cr, width * 0.25, 0);
	cairo_line_to(cr, width * 0.25, height);
	cairo_move_to(cr, width * 0.5, 0);
	cairo_line_to(cr, width * 0.5, height);
	cairo_move_to(cr, width * 0.75, 0);
	cairo_line_to(cr, width * 0.75, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.25);
	cairo_line_to(cr, width, height * 0.25);
	cairo_move_to(cr, 0, height * 0.5);
	cairo_line_to(cr, width, height * 0.5);
	cairo_move_to(cr, 0, height * 0.75);
	cairo_line_to(cr, width, height * 0.75);

	cairo_stroke(cr);

	cairo_set_line_width(cr, 1.0);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_move_to(cr, width * 0.125, 0);
	cairo_line_to(cr, width * 0.125, height);
	cairo_move_to(cr, width * 0.375, 0);
	cairo_line_to(cr, width * 0.375, height);
	cairo_move_to(cr, width * 0.625, 0);
	cairo_line_to(cr, width * 0.625, height);
	cairo_move_to(cr, width * 0.875, 0);
	cairo_line_to(cr, width * 0.875, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.125);
	cairo_line_to(cr, width, height * 0.125);
	cairo_move_to(cr, 0, height * 0.375);
	cairo_line_to(cr, width, height * 0.375);
	cairo_move_to(cr, 0, height * 0.625);
	cairo_line_to(cr, width, height * 0.625);
	cairo_move_to(cr, 0, height * 0.875);
	cairo_line_to(cr, width, height * 0.875);
	cairo_stroke(cr);
}

void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV, gboolean isOrig, gboolean is_log,
		gboolean is_mono) {
	if (!histo) return;
	if (width <= 0) return;
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

	float vals_per_px = (float)norm / (float)width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static gfloat *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		gfloat *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(gfloat));
		if (!tmp) {
			if (displayed_values != NULL) {
				g_free(displayed_values);
				displayed_values = NULL;
			}
			PRINT_ALLOC_ERR;
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}
	if (is_mono)
		cairo_set_source_rgb(cr, 255.0, 255.0, 255.0);
	else if (layer < 0)
		cairo_set_source_rgb(cr, 255.0, 255.0, 0.0);
	else
		cairo_set_source_rgb(cr, histo_color_r[layer], histo_color_g[layer],
				histo_color_b[layer]);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);
	if (isOrig || layer == -2)
		cairo_set_line_width(cr, 0.5);

	// first loop builds the bins and finds the maximum
	i = 0;
	double graph_height = 0.f;
	current_bin = 0;
	do {
		double bin_val = 0.f;
		while (i < nb_orig_bins
				&& (double)i / vals_per_px <= current_bin + 0.5f) {
			bin_val += gsl_histogram_get(histo, i);
			i++;
		}
		if (is_log && bin_val != 0.f) {
			bin_val = logf(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);
	if (!graph_height)
		return;
	for (i = 0; i < nb_bins_allocated; i++) {
		double bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}
