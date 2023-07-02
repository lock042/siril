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

#include "siril_plot.h"

#include <cairo.h>
#include <math.h>
#include "core/siril_log.h"

// static variables

#define GUIDE 20 // used to determine number of tics and spacing

// TODO: these variables should not be static
static char *fmtx = NULL, *fmty = NULL;

static void spl_formatX(double v, char *buf, size_t bufsz) {
	snprintf(buf, 128, fmtx, v);
}

static void spl_formatY(double v, char *buf, size_t bufsz) {
	snprintf(buf, 128, fmty, v);
}

// static functions

static void free_xyplot_data(splxydata *plot) {
	free(plot->data);
	g_free(plot->label);
	plot = NULL;
}

static void free_xyerrplot_data(splxyerrdata *plots) {
	g_free(plots->label);
	for (int i = 0; i < 3; i++)
		free_xyplot_data(plots->plots[i]);
	plots = NULL;
}

static void free_list_plot(gpointer data) {
	splxydata *item = (splxydata *)data;
	free_xyplot_data(data);
	g_slice_free(splxydata, item);
}

static void free_list_plots(gpointer data) {
	splxyerrdata *item = (splxyerrdata *)data;
	free_xyerrplot_data(data);
	g_slice_free(splxyerrdata, item);
}

// allocate a simple xy data structure
static splxydata *alloc_xyplot_data(int nb) {
	splxydata *plot = malloc(sizeof(splxydata));
	plot->data = calloc(nb, sizeof(struct kpair));
	if (!plot->data) {
		PRINT_ALLOC_ERR;
		free(plot);
		return NULL;
	}
	plot->nb = nb;
	plot->label = NULL;
	return plot;
}

// allocate a xyerr data structure to hold data, errp and errm
static splxyerrdata *alloc_xyerrplot_data(int nb) {
	splxyerrdata *plots = malloc(sizeof(splxyerrdata));
	plots->label = NULL;
	plots->nb = nb;
	gboolean ok = TRUE;
	for (int i = 0; i < 3; i++) {
		plots->plots[i] = alloc_xyplot_data(nb);
		if (!plots->plots[i])
			ok = FALSE;
	}
	if (!ok) // one of the 3 allocations failed
		free_xyerrplot_data(plots);
	return plots;
}

// init/free spl_data

void init_siril_plot_data(siril_plot_data *spl_data) {
	spl_data->plot = NULL;
	spl_data->plots = NULL;
	spl_data->title = NULL;
	spl_data->xlabel = NULL;
	spl_data->ylabel = NULL;
	spl_data->xfmt = NULL;
	spl_data->yfmt = NULL;
	spl_data->plottype = KPLOT_LINES;
	spl_data->plotstype = KPLOTS_YERRORBAR;
	spl_data->datamin = (point){ DBL_MAX, DBL_MAX};
	spl_data->datamax = (point){ -DBL_MAX, -DBL_MAX};

	// initializing kplot cfg structs
	kplotcfg_defaults(&spl_data->cfgplot);
	kdatacfg_defaults(&spl_data->cfgdata);
	spl_data->cfgplot.ticlabel = TICLABEL_LEFT | TICLABEL_BOTTOM;
	spl_data->cfgplot.border = BORDER_ALL;
	spl_data->cfgplot.borderline.clr.type = KPLOTCTYPE_RGBA;
	spl_data->cfgplot.borderline.clr.rgba[0] = 0.5;
	spl_data->cfgplot.borderline.clr.rgba[1] = 0.5;
	spl_data->cfgplot.borderline.clr.rgba[2] = 0.5;
	spl_data->cfgplot.borderline.clr.rgba[3] = 1.0;
	spl_data->cfgplot.yaxislabelrot = M_PI_2 * 3.0;
	spl_data->cfgdata.line.sz = 0.5;
}

void free_siril_plot_data(siril_plot_data *spl_data) {
	// freeing gchars
	g_free(spl_data->title);
	g_free(spl_data->xlabel);
	g_free(spl_data->ylabel);
	g_free(spl_data->xfmt);
	g_free(spl_data->yfmt);
	// freeing the xy and xyerr plots
	g_list_free_full(spl_data->plot, (GDestroyNotify)free_list_plot);
	g_list_free_full(spl_data->plots, (GDestroyNotify)free_list_plots);
	//freeing kplot cfg structures
	// TODO: check if smthg is required

}

// setters
void siril_plot_set_title(siril_plot_data *spl_data, const gchar *title) {
	if (spl_data->title)
		g_free(spl_data->title);
	spl_data->title = g_strdup(title);
}

void siril_plot_set_xlabel(siril_plot_data *spl_data, const gchar *xlabel) {
	if (spl_data->xlabel)
		g_free(spl_data->xlabel);
	spl_data->xlabel = g_strdup(xlabel);
}

void siril_plot_set_ylabel(siril_plot_data *spl_data, const gchar *ylabel) {
	if (spl_data->ylabel)
		g_free(spl_data->ylabel);
	spl_data->ylabel = g_strdup(ylabel);
}

void siril_plot_set_xfmt(siril_plot_data *spl_data, const gchar *xfmt) {
	if (spl_data->xfmt)
		g_free(spl_data->xfmt);
	spl_data->xfmt = g_strdup(xfmt);
}

void siril_plot_set_yfmt(siril_plot_data *spl_data, const gchar *yfmt) {
	if (spl_data->yfmt)
		g_free(spl_data->yfmt);
	spl_data->yfmt = g_strdup(yfmt);
}

// utilities
static gboolean siril_plot_autotic(double vmin, double vmax, int *nbtics, double *tmin, double *tmax) {
	double extent = vmax - vmin;
	if (extent <= 0.)
		return FALSE;
	double power = pow(10, floor(log10(extent)));
	double xnorm = extent / power;
	double posns = (double)GUIDE / xnorm;
	double tics;
	if (posns > 40)
		tics = 0.05;
	else if (posns > 20)
		tics = 0.1;
	else if (posns > 10)
		tics = 0.2;
	else if (posns > 4)
		tics = 0.5;
	else if (posns > 2)
		tics = 1;
	else if (posns > 0.5)
		tics = 2;
	else
		tics = ceil(xnorm);
	tics *= power;
	*tmin = floor(vmin / tics) * tics;
	*tmax = ceil(vmax / tics) * tics;
	*nbtics = (int)((*tmax - *tmin) / tics) + 1;
	// siril_debug_print("autotic:\t%g\t%g=>%d\t%g\t%g\n", vmin, vmax, *nbtics, *tmin, *tmax);
	return TRUE;
}

// data assignment
gboolean siril_plot_add_xydata(siril_plot_data *spl_data, gchar *label, size_t nb, double *x, double *y, double *errp, double *errm) {
	// single plot case
	if (!errp) {
		// allocate data
		splxydata *plot = alloc_xyplot_data(nb);
		if (!plot) {
			siril_debug_print("Could not allocate plot data\n");
			return FALSE;
		}
		// fill and update spl_data bounds
		for (int i = 0; i < nb; i++) {
			plot->data[i].x = x[i];
			plot->data[i].y = y[i];
			if (x[i] < spl_data->datamin.x) spl_data->datamin.x = x[i];
			if (y[i] < spl_data->datamin.y) spl_data->datamin.y = y[i];
			if (x[i] > spl_data->datamax.x) spl_data->datamax.x = x[i];
			if (y[i] > spl_data->datamax.y) spl_data->datamax.y = y[i];
		}
		if (label)
			plot->label = g_strdup(label);
		// and append to plot GList
		spl_data->plot = g_list_append(spl_data->plot, plot);
		return TRUE;
	}
	// xyerror plot case
	splxyerrdata *plots = alloc_xyerrplot_data(nb);
	if (!plots) {
		siril_debug_print("Could not allocate plots data\n");
		return FALSE;
	}
	// if no errm is passed, we assume it is the same as errp
	if (!errm)
		errm = errp;
	// fill and update spl_data bounds
	for (int i = 0; i < nb; i++) {
		plots->plots[0]->data[i].x = x[i];
		plots->plots[0]->data[i].y = y[i];
		plots->plots[1]->data[i].x = x[i];
		plots->plots[1]->data[i].y = errp[i];
		plots->plots[2]->data[i].x = x[i];
		plots->plots[2]->data[i].y = errm[i];
		if (x[i] < spl_data->datamin.x) spl_data->datamin.x = x[i];
		if (y[i] - errm[i] < spl_data->datamin.y) spl_data->datamin.y = y[i] - errm[i];
		if (x[i] > spl_data->datamax.x) spl_data->datamax.x = x[i];
		if (y[i] + errp[i] > spl_data->datamax.y) spl_data->datamax.y = y[i] + errp[i];
	}
	if (label)
		plots->label = g_strdup(label);
	// and append to plots GList
	spl_data->plots = g_list_append(spl_data->plots, plots);
	return TRUE;
}

// draw the data contained in spl_data to the cairo context cr
gboolean siril_plot_draw(cairo_t *cr, siril_plot_data *spl_data, double width, double height) {
	struct kdata *d1 = NULL, *d2 = NULL;
	double color = 1.0;
	if (spl_data->xlabel)
		spl_data->cfgplot.xaxislabel = spl_data->xlabel;
	if (spl_data->ylabel)
		spl_data->cfgplot.yaxislabel = spl_data->ylabel;
	if (spl_data->xfmt) {
		if (fmtx)
			g_free(fmtx);
		fmtx = g_strdup(spl_data->xfmt);
		spl_data->cfgplot.xticlabelfmt = spl_formatX;
	}
	if (spl_data->yfmt) {
		if (fmty)
			g_free(fmty);
		fmty = g_strdup(spl_data->yfmt);
		spl_data->cfgplot.yticlabelfmt = spl_formatY;
	}

	// computing the tics spacing and bounds
	double xmin, xmax, ymin, ymax;
	int nbticX,nbticY;
	if (siril_plot_autotic(spl_data->datamin.x, spl_data->datamax.x, &nbticX, &xmin, &xmax) &&
		siril_plot_autotic(spl_data->datamin.y, spl_data->datamax.y, &nbticY, &ymin, &ymax)) {
		spl_data->cfgplot.extrema = 0x0F;
		spl_data->cfgplot.extrema_xmin = xmin;
		spl_data->cfgplot.extrema_xmax = xmax;
		spl_data->cfgplot.extrema_ymin = ymin;
		spl_data->cfgplot.extrema_ymax = ymax;
		spl_data->cfgplot.xtics = nbticX;
		spl_data->cfgplot.ytics = nbticY;
	}

	struct kplot *p = kplot_alloc(&spl_data->cfgplot);

	// data plots
	int nb_graphs = 0;

	// xylines
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
		d1 = kdata_array_alloc(plot->data, plot->nb);
		kplot_attach_data(p, d1, spl_data->plottype, &spl_data->cfgdata);
		kdata_destroy(d1);
		d1 = NULL;
		nb_graphs++;
	}
	// // xy points with y error bars
	// while (plots) {
	// 	d2 = kdata_array_alloc(plot->data, plot->nb);
	// 	kplot_attach_data(p, d1, spl_data->plotstype, &spl_data->cfgdata);
	// 	plots = plots->nextplots;
	// 	kdata_destroy(d2);
	// 	d2 = NULL;
	// 	nb_graphs++;
	// }

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

	// creating a surface to draw the plot (accounting for title-reserved space if required)
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
		cairo_select_font_face(cr, spl_data->cfgplot.axislabelfont.family, 
		spl_data->cfgplot.axislabelfont.slant, spl_data->cfgplot.axislabelfont.weight);
		cairo_set_font_size(cr, spl_data->cfgplot.axislabelfont.sz * 1.2);
		cairo_text_extents_t te1;
		cairo_text_extents(cr, spl_data->title, &te1); // getting the dimensions of the textbox
		cairo_translate(cr, (int)((width - te1.width) * 0.5), (int)(top - te1.height * 0.5));
		cairo_show_text(cr, spl_data->title);
		cairo_stroke(cr);
		cairo_restore(cr); // restore the orginal context
	}
	return TRUE;
}

// draw the data contained in spl_data and saves as png file
gboolean siril_plot_save_png(siril_plot_data *spl_data, char *pngfilename) {
	gboolean success = TRUE;
	cairo_t *png_cr = NULL;
	//create the surface
	cairo_surface_t *png_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, SIRIL_PLOT_PNG_WIDTH, SIRIL_PLOT_PNG_HEIGHT);
	if (cairo_surface_status(png_surface)) {
		siril_debug_print("Could not create png surface\n");
		success = FALSE;
	}
	//create the context
	if (success) {
		png_cr = cairo_create(png_surface);
		if (cairo_status(png_cr)) {
			siril_debug_print("Could not create png context\n");
			success = FALSE;
		}
	}
	// draw the plot and save the surface to png
	if (success && siril_plot_draw(png_cr, spl_data, (double)SIRIL_PLOT_PNG_WIDTH, (double)SIRIL_PLOT_PNG_HEIGHT)) {
		siril_debug_print("Successfully created png plot\n");
		if (!cairo_surface_write_to_png(png_surface, pngfilename))
			siril_log_message(_("%s has been saved.\n"), pngfilename);
		else
			success = FALSE;
	}
	if (png_cr)
		cairo_destroy(png_cr);
	cairo_surface_destroy(png_surface);
	return success;
}









