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

#ifndef _HISTOGRAM_UTILS_H_
#define _HISTOGRAM_UTILS_H_

#include <gsl/gsl_histogram.h>
#include "core/siril.h"

/* Histogram size */
size_t get_histo_size(fits *fit);

/* Histogram computation */
gsl_histogram *computeHisto(fits *fit, int layer);
gsl_histogram *computeHisto_Selection(fits *fit, int layer, rectangle *selection);

/* Global histogram cache management (com.layers_hist[]) */
void set_histogram(gsl_histogram *histo, int layer);
void set_sat_histogram(gsl_histogram *histo);
void compute_histo_for_fit(fits *thefit);
void invalidate_gfit_histogram(void);

/* Cairo rendering primitives shared between histogram dialog, curves dialog, and remixer */
void fill_histo_background(cairo_t *cr, int width, int height);
void draw_grid(cairo_t *cr, int width, int height);
void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width, int height,
                   double zoomH, double zoomV, gboolean isOrig, gboolean is_log, gboolean is_mono);

#endif /* _HISTOGRAM_UTILS_H_ */
