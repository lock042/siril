/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include "io/siril_plot.h"

#include <cairo.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gui/utils.h"
#include "core/siril_log.h"

// TODO: these variables should not be static
static char *fmtx = NULL, *fmty = NULL;

static void spl_formatX(double v, char *buf, size_t bufsz) {
	snprintf(buf, 128, fmtx, v);
}

static void spl_formatY(double v, char *buf, size_t bufsz) {
	snprintf(buf, 128, fmty, v);
}

// static void (*spl_format(char *fmt, ...))(double, char *, size_t){
// 	// format(fmt);
// 	va_list parametersInfos;
// 	/* Initialize the va_list structure */
// 	va_start(parametersInfos, fmt);
// 	double v = (double) va_arg(parametersInfos, double);
// 	char *buf = (char *) va_arg(parametersInfos, char *);
// 	// size_t bufsz = (size_t) va_arg(parametersInfos, size_t);
// 	va_end(parametersInfos); 
// 	snprintf(buf, 128, fmt, v);
// }

static gboolean siril_plot_draw(cairo_t *cr, siril_plot_data *spl_data, double width, double height) {
	struct kdata *d1 = NULL, *d2 = NULL;
	splxydata *plot = spl_data->plot;
	splxydata *plots = spl_data->plots;
	double color = 1.0;
	if (spl_data->xlabel)
		spl_data->kplotcfg.xaxislabel = spl_data->xlabel;
	if (spl_data->ylabel)
		spl_data->kplotcfg.yaxislabel = spl_data->ylabel;
	if (spl_data->xfmt) {
		if (fmtx)
			g_free(fmtx);
		fmtx = g_strdup(spl_data->xfmt);
		spl_data->kplotcfg.xticlabelfmt = spl_formatX;
	}
	if (spl_data->yfmt) {
		if (fmty)
			g_free(fmty);
		fmty = g_strdup(spl_data->yfmt);
		spl_data->kplotcfg.yticlabelfmt = spl_formatY;
	}

	// computing the tics spacing and bounds
	double xmin, xmax, ymin, ymax;
	int nbticX,nbticY;
	if (siril_plot_autotic(spl_data->datamin.x, spl_data->datamax.x, &nbticX, &xmin, &xmax) &&
		siril_plot_autotic(spl_data->datamin.y, spl_data->datamax.y, &nbticY, &ymin, &ymax)) {
		spl_data->kplotcfg.extrema = 0x0F;
		spl_data->kplotcfg.extrema_xmin = xmin;
		spl_data->kplotcfg.extrema_xmax = xmax;
		spl_data->kplotcfg.extrema_ymin = ymin;
		spl_data->kplotcfg.extrema_ymax = ymax;
		spl_data->kplotcfg.xtics = nbticX;
		spl_data->kplotcfg.ytics = nbticY;
	}

	struct kplot *p = kplot_alloc(&spl_data->kplotcfg);
	struct kdatacfg cfgdata;
	kdatacfg_defaults(&cfgdata);
	cfgdata.line.sz = 0.5;

	// data plots
	int nb_graphs = 0;

	// xylines
	while (plot) {
		d1 = kdata_array_alloc(plot->data, plot->nb);
		kplot_attach_data(p, d1, spl_data->plottype, &cfgdata);
		plot = plot->next;
		kdata_destroy(d1);
		d1 = NULL;
		nb_graphs++;
	}
	// xy points with y error bars
	while (plots) {
		d2 = kdata_array_alloc(plot->data, plot->nb);
		kplot_attach_data(p, d1, spl_data->plottype, NULL);
		plots = plots->nextplots;
		kdata_destroy(d2);
		d2 = NULL;
		nb_graphs++;
	}

	// preparing the surfaces
	double drawwidth = width;
	double drawheight = height;
	double top = 0.;
	
	// booking space for title
	if (spl_data->title) {
		top = SPL_TITLE_RATIO * height;
		drawheight = height - top;
	}

	// painting the whole surface white
	cairo_set_source_rgb(cr, color, color, color);
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_fill(cr);

	// creating a surface to draw the plot (accounting for title reserved space if required)
	cairo_surface_t *draw_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, (int)drawwidth, (int)drawheight);
	cairo_t *draw_cr = cairo_create(draw_surface);
	kplot_draw(p, drawwidth, drawheight, draw_cr);
	cairo_set_source_surface(cr, draw_surface, 0., top);
	cairo_paint(cr);
	cairo_surface_destroy(draw_surface);
	kplot_free(p);

	// writing the title if any
	if (spl_data->title) {
		cairo_save(cr); // save the orginal context
		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		cairo_select_font_face(cr, spl_data->kplotcfg.axislabelfont.family, 
		spl_data->kplotcfg.axislabelfont.slant, spl_data->kplotcfg.axislabelfont.weight);
		cairo_set_font_size(cr, spl_data->kplotcfg.axislabelfont.sz * 1.5);
		cairo_text_extents_t te1;
		cairo_text_extents(cr, spl_data->title, &te1); // getting the dimensions of the textbox
		cairo_translate(cr, (width - te1.width) * 0.5, top - te1.height * 0.5);
		cairo_show_text(cr, spl_data->title);
		cairo_stroke(cr);
		cairo_restore(cr); // restore the orginal context
	}

	if (spl_data->pngfilename) {
		cairo_surface_t *png_surface = cairo_get_target(cr);
		cairo_surface_write_to_png(png_surface, spl_data->pngfilename);
		siril_log_message(_("%s has been saved.\n"), spl_data->pngfilename);
		g_free(spl_data->pngfilename);
		spl_data->pngfilename = NULL;
	}
	return TRUE;
}

static gboolean on_siril_plot_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data) {

	siril_plot_data *spl_data = (siril_plot_data *)user_data;
	if (!spl_data || !widget)
		return FALSE;

	double width =  gtk_widget_get_allocated_width(widget);
	double height = gtk_widget_get_allocated_height(widget);
	if (!siril_plot_draw(cr, spl_data, width, height))
		siril_debug_print("Problem while creating siril_plot\n");
	return FALSE;
}

gboolean create_new_siril_plot_window(gpointer p) {
	GtkWidget *window;
	GtkWidget *da;
	siril_plot_data *spl_data = (siril_plot_data *)p;

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "Siril plot");
	// // g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK(clear_siril_plot_data), spl_data);

	da = gtk_drawing_area_new();
	gtk_widget_set_size_request(da, SIRIL_PLOT_WIDTH, SIRIL_PLOT_HEIGHT);
	gtk_container_add(GTK_CONTAINER(window), da);

	g_signal_connect(G_OBJECT(da), "draw", G_CALLBACK(on_siril_plot_draw), spl_data);
	gtk_widget_show_all(window);
	return FALSE;
}










