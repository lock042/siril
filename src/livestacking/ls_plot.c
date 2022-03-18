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

#include <math.h>

#include "core/siril.h"
#include "kplot.h"

static void format(double v, char *buf, size_t bufsz) {
	snprintf(buf, bufsz, "%0.1lf", v);
}

static void format_void(double v, char *buf, size_t bufsz) {
	snprintf(buf, bufsz, "%s", "");
}

gboolean on_ls_plot_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	struct kplotcfg cfgplot;
	struct kpair points1[50], points2[50];
	struct kdata *d1;
	struct kplot *p;
	size_t i;

	kplotcfg_defaults(&cfgplot);

	cfgplot.yticlabelfmt = format;
	cfgplot.xticlabelfmt = format_void;

	for (i = 0; i < 50; i++) {
		points1[i].x = points2[i].x = i;
		points1[i].y = log((i + 1) / 50.0);
		points2[i].y = -log((i + 1) / 50.0) + points1[0].y;
	}
	d1 = kdata_array_alloc(points1, 50);
	p = kplot_alloc(&cfgplot);
	kplot_attach_data(p, d1, KPLOT_LINES, NULL);

	gint width = gtk_widget_get_allocated_width(widget);
	gint height = gtk_widget_get_allocated_height(widget);

	kplot_draw(p, width, height, cr);

	kdata_destroy(d1);
	kplot_free(p);


	return FALSE;
}
