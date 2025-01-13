/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include "siril_plot.h"

#include <cairo.h>
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo/cairo-svg.h>
#endif
#include <pango/pangocairo.h>
#include <math.h>
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/plot.h"


// static variables

#define GUIDE 20 // used to determine number of tics and spacing


// static functions

static void free_xyplot_data(splxydata *plot) {
	free(plot->data);
	g_free(plot->label);
}

static void free_xyerrplot_data(splxyerrdata *plots) {
	g_free(plots->label);
	for (int i = 0; i < 3; i++)
		free_xyplot_data(plots->plots[i]);
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

static void free_bkg(splbkg *bkg) {
	if (!bkg)
		return;
	g_free(bkg->bkgfilepath);
	if (bkg->img)
		g_object_unref(bkg->img);
	free(bkg);
}

// allocate a simple xy data structure
static splxydata *alloc_xyplot_data(int nb) {
	splxydata *plot = g_slice_new(splxydata);
	plot->data = calloc(nb, sizeof(struct kpair));
	if (!plot->data) {
		PRINT_ALLOC_ERR;
		g_slice_free(splxydata, plot);
		return NULL;
	}
	plot->nb = nb;
	plot->label = NULL;
	plot->pl_type = KPLOT_UNDEFINED;
	return plot;
}

// allocate a xyerr data structure to hold data, errp and errm
static splxyerrdata *alloc_xyerrplot_data(int nb) {
	splxyerrdata *plots = g_slice_new(splxyerrdata);
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

static void free_legend(void *data) {
	spllegend *legend = (spllegend*) data;
	g_slice_free(spllegend, legend);
}

// allocate a legend entry
static spllegend *new_legend_entry(spl_type type, const double color[3]) {
	spllegend *legend = g_slice_new(spllegend);
	legend->type = type;
	memcpy(legend->color, color, 3 * sizeof(double));
	return legend;
}

// sort kpairs by ascending x
static int comparex(const void *a, const void *b) {
	struct kpair datax_a = *((struct kpair *) a);
	struct kpair datax_b = *((struct kpair *) b);
	if (datax_a.x > datax_b.x) return 1;
	if (datax_a.x < datax_b.x) return -1;
	return 0;
}

// TODO later (zoomable background)
// get subpixmap to display background
// static GdkPixbuf *extract_sub_bkg(siril_plot_data *spl_data, double xmin, double xmax, double ymin, double ymax, int *offsx, int *offsy) {
// 	*offsx = 0;
// 	*offsy = 0;
// 	double deltax = spl_data->datamax.x - spl_data->datamin.x;
// 	double deltay = spl_data->datamax.y - spl_data->datamin.y;
// 	double xminp = xmin * (double)spl_data->bkg->width / deltax;
// 	double xmaxp = xmax * (double)spl_data->bkg->width / deltax;
// 	double yminp = (spl_data->datamax.y - ymin) * (double)spl_data->bkg->width / deltax;// the display has y down
// 	double ymaxp = (spl_data->datamax.y - ymax) * (double)spl_data->bkg->width / deltax;
// 	double rangex = (xmax - xmin) * (double)spl_data->bkg->width / deltax;
// 	double rangey = (ymax - ymin) * (double)spl_data->bkg->height / deltay;
// 	if (xminp < 0) {
// 		*offsx = (int)xminp;
// 		xminp = 0.;
// 	}
// 	if (ymaxp < 0) {
// 		*offsy = (int)ymaxp;
// 		ymaxp = 0.;
// 	}
// 	if (rangex > (double)spl_data->bkg->width) {
// 		rangex = (double)spl_data->bkg->width;
// 	}
// 	if (rangey > (double)spl_data->bkg->height) {
// 		rangey = (double)spl_data->bkg->height;
// 	}
// 	GdkPixbuf *subbkg = gdk_pixbuf_new_subpixbuf(spl_data->bkg->img, (int)xminp, (int)ymaxp, (int)(rangex), (int)(rangey));
// 	return subbkg;
// }

// init/free spl_data

void init_siril_plot_data(siril_plot_data *spl_data) {
	spl_data->plot = NULL;
	spl_data->plots = NULL;
	spl_data->title = NULL;
	spl_data->xlabel = NULL;
	spl_data->ylabel = NULL;
	spl_data->xfmt = NULL;
	spl_data->yfmt = NULL;
	spl_data->savename = NULL;
	spl_data->forsequence = FALSE;
	spl_data->plottype = KPLOT_LINES;
	spl_data->plotstype = KPLOTS_YERRORBAR;
	spl_data->plotstypes[0] = KPLOT_POINTS;
	spl_data->plotstypes[1] = KPLOT_HYPHENS;
	spl_data->plotstypes[2] = KPLOT_HYPHENS;
	spl_data->datamin = (point){ DBL_MAX, DBL_MAX};
	spl_data->datamax = (point){ -DBL_MAX, -DBL_MAX};
	spl_data->show_legend = TRUE;
	spl_data->autotic = TRUE;
	spl_data->revertX = FALSE;
	spl_data->revertY = FALSE;
	spl_data->zoomable = FALSE;
	spl_data->bkg = NULL;
	spl_data->width = 0;
	spl_data->height = 0;

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
	spl_data->cfgplot.ticlabelfont.family = SIRIL_PLOT_FONT_FAMILY;
	spl_data->cfgplot.axislabelfont.family = SIRIL_PLOT_FONT_FAMILY;
	spl_data->cfgdata.line.sz = 0.5;
	spl_data->cfgdata.point.sz = 0.5;
	size_t clrsz;

	// we init the colors here to ease creating the legend entries
	struct kplotccfg *cfgcolor;
	kplotcfg_default_palette(&cfgcolor, &clrsz);
	spl_data->cfgplot.clrsz = clrsz;
	spl_data->cfgplot.clrs = cfgcolor;

	// initializing the plot_draw_data
	memset(&spl_data->pdd, 0, sizeof(plot_draw_data_t));
}

void free_siril_plot_data(siril_plot_data *spl_data) {
	if (!spl_data)
		return;
	// freeing gchars
	g_free(spl_data->title);
	g_free(spl_data->xlabel);
	g_free(spl_data->ylabel);
	g_free(spl_data->xfmt);
	g_free(spl_data->yfmt);
	g_free(spl_data->savename);
	// freeing the xy and xyerr plots
	g_list_free_full(spl_data->plot, (GDestroyNotify)free_list_plot);
	g_list_free_full(spl_data->plots, (GDestroyNotify)free_list_plots);
	//freeing kplot cfg structures
	free(spl_data->cfgplot.clrs);
	free_bkg(spl_data->bkg);
	free(spl_data);

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

void siril_plot_set_savename(siril_plot_data *spl_data, const gchar *savename) {
	if (spl_data->savename)
		g_free(spl_data->savename);
	spl_data->savename = g_strdup(savename);
}

// set the color of the nth plot (n is one-based)
// e.g.: siril_plot_set_nth_color(spl_data, 1, (double[3]){1., 0., 0.});
// sets the first series color to red
void siril_plot_set_nth_color(siril_plot_data *spl_data, int n, double color[3]) {
	if (n > spl_data->cfgplot.clrsz) {
		siril_debug_print("can't add color out of palette size (%lu)\n", spl_data->cfgplot.clrsz);
		return;
	}
	memcpy(spl_data->cfgplot.clrs[n - 1].rgba, color, 3 * sizeof(double));
}

// set the type of the nth plot (n is one-based)
// Note: spl_data->plot member is the list of xydata (w/o error bars)
// e.g.: siril_plot_set_nth_plot_type(spl_data, 1, KPLOT_LINESMARKS);
// sets the first series type to line with cross markers
// This overiddes the setting defined in spl_data->plottype which applies
// if no specific type is passed for a series
void siril_plot_set_nth_plot_type(siril_plot_data *spl_data, int n, enum kplottype pl_type) {
	GList *current_entry = g_list_nth(spl_data->plot, n - 1);
	if (!current_entry) {
		siril_debug_print("can't add plot type out of plot list size\n");
		return;
	}
	splxydata *plot = (splxydata *)current_entry->data;
	plot->pl_type = pl_type;
}

// set an image to be used as background
// loads the image as a GdkPixBuf and gets its dimensions
// `bkgfilename` is the name of the bkg file which should be added in
// `pixmaps/plot_background` folder and added to siril_resource.xml
gboolean siril_plot_set_background(siril_plot_data *spl_data, const gchar *bkgfilename) {
	if (spl_data->bkg) {
		free_bkg(spl_data->bkg);
	}
	GError *error = NULL;
	spl_data->bkg = calloc(1, sizeof(splbkg));
	spl_data->bkg->bkgfilepath = g_build_filename("/org/siril/ui/pixmaps/plot_background", bkgfilename, NULL);
	spl_data->bkg->img = gdk_pixbuf_new_from_resource(spl_data->bkg->bkgfilepath, &error);
	if (error) {
		free_bkg(spl_data->bkg);
		siril_debug_print("can't load background image %s (Error: %s)", spl_data->bkg->bkgfilepath, error->message);
		g_error_free(error);
		return FALSE;
	}
	spl_data->bkg->height = gdk_pixbuf_get_height(spl_data->bkg->img);
	spl_data->bkg->width = gdk_pixbuf_get_width(spl_data->bkg->img);
	return TRUE;
}

// set the types of the nth plots (n is one-based)
// Note: spl_data->plots member is the list of xydata + error bars
// e.g.: siril_plot_set_nth_plot_types(spl_data, 1, (enum kplottype[3]){KPLOT_MARKS, KPLOT_POINTS, KPLOT_POINTS});
// sets the first series with error bars to crosses with circle error bar markers
// This overiddes the setting defined in spl_data->plotstypes which applies
// if no specific type is passed for a series
void siril_plot_set_nth_plots_types(siril_plot_data *spl_data, int n, enum kplottype pl_type[3]) {
	GList *current_entry = g_list_nth(spl_data->plots, n - 1);
	if (!current_entry) {
		siril_debug_print("can't add plots types out of plots list size\n");
		return;
	}
	splxyerrdata *plots = (splxyerrdata *)current_entry->data;
	for (int i = 0; i < 3; i++) {
		plots->plots[i]->pl_type = pl_type[i];
	}
}

// utilities
static gboolean siril_plot_autotic(double vmin, double vmax, int *nbtics, double *tmin, double *tmax, int *sig) {
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
	*nbtics = (int)round(((*tmax - *tmin) / tics)) + 1;
	//computing number of decimals
	double logtics = log10(tics);
	*sig = abs((int)floor(min(0., logtics)));
	// siril_debug_print("autotic:\t%g\t%g=>%d\t%g\t%g\n", vmin, vmax, *nbtics, *tmin, *tmax);
	return TRUE;
}

// data assignment
gboolean siril_plot_add_xydata(siril_plot_data *spl_data, const gchar *label, size_t nb, const double *x, const double *y, const double *errp, const double *errm) {
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

// sort all plots by ascending x
void siril_plot_sort_x(siril_plot_data *spl_data) {
	if (!spl_data)
		return;
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
		qsort(plot->data, plot->nb, sizeof(struct kpair), comparex);
	}
	for (GList *list = spl_data->plots; list; list = list->next) {
		splxyerrdata *plots = (splxyerrdata *)list->data;
		for (int i = 0; i < 3; i++) {
			splxydata *plot = plots->plots[i];
			qsort(plot->data, plot->nb, sizeof(struct kpair), comparex);
		}
	}
}

// draw the data contained in spl_data to the cairo context cr
gboolean siril_plot_draw(cairo_t *cr, siril_plot_data *spl_data, double width, double height, gboolean for_svg) {
	struct kdata *d1 = NULL, *d2[3];
	double color = 1.0;
	if (spl_data->xlabel)
		spl_data->cfgplot.xaxislabel = spl_data->xlabel;
	if (spl_data->ylabel)
		spl_data->cfgplot.yaxislabel = spl_data->ylabel;
	spl_data->cfgplot.xaxisrevert = (spl_data->revertX) ? 1 : 0;
	spl_data->cfgplot.yaxisrevert = (spl_data->revertY) ? 1 : 0;

	// computing the tics spacing and bounds
	double x1 = (!spl_data->zoomable) ? spl_data->datamin.x : spl_data->pdd.datamin.x;
	double x2 = (!spl_data->zoomable) ? spl_data->datamax.x : spl_data->pdd.datamax.x;
	double y1 = (!spl_data->zoomable) ? spl_data->datamin.y : spl_data->pdd.datamin.y;
	double y2 = (!spl_data->zoomable) ? spl_data->datamax.y : spl_data->pdd.datamax.y;
	double xmin, xmax, ymin, ymax;
	int nbticX, nbticY, sigX, sigY;
	if (spl_data->autotic &&
		siril_plot_autotic(x1, x2, &nbticX, &xmin, &xmax, &sigX) &&
		siril_plot_autotic(y1, y2, &nbticY, &ymin, &ymax, &sigY)) {
		spl_data->cfgplot.extrema = 0x0F;
		spl_data->cfgplot.extrema_xmin = xmin;
		spl_data->cfgplot.extrema_xmax = xmax;
		spl_data->cfgplot.extrema_ymin = ymin;
		spl_data->cfgplot.extrema_ymax = ymax;
		spl_data->cfgplot.xtics = nbticX;
		spl_data->cfgplot.ytics = nbticY;
		spl_data->pdd.pdatamin = (point){xmin, ymin};
		spl_data->pdd.pdatamax = (point){xmax, ymax};
		// if the formats are not forced by caller, they are adjusted
		if (!spl_data->xfmt) {
			g_free(spl_data->cfgplot.xticlabelfmtstr);
			spl_data->cfgplot.xticlabelfmtstr = g_strdup_printf("%%.%df", sigX);
		}
		if (!spl_data->yfmt) {
			g_free(spl_data->cfgplot.yticlabelfmtstr);
			spl_data->cfgplot.yticlabelfmtstr = g_strdup_printf("%%.%df", sigY);
		}
	} else {  // fallback
		spl_data->cfgplot.extrema = 0x0F;
		spl_data->cfgplot.extrema_xmin = x1;
		spl_data->cfgplot.extrema_xmax = x2;
		spl_data->cfgplot.extrema_ymin = y1;
		spl_data->cfgplot.extrema_ymax = y2;
		spl_data->pdd.pdatamin = (point){x1, y1};
		spl_data->pdd.pdatamax = (point){x2, y2};
		if (spl_data->autotic) {
			spl_data->cfgplot.xtics = 5;
			spl_data->cfgplot.ytics = 5;
			if (!spl_data->xfmt) {
				g_free(spl_data->cfgplot.xticlabelfmtstr);
				spl_data->cfgplot.xticlabelfmtstr = g_strdup("%g");
			}
			if (!spl_data->yfmt) {
				g_free(spl_data->cfgplot.yticlabelfmtstr);
				spl_data->cfgplot.yticlabelfmtstr = g_strdup("%g");
			}
		}
	}
	// if the formats are forced by caller, they are passed
	if (spl_data->xfmt) {
		g_free(spl_data->cfgplot.xticlabelfmtstr);
		spl_data->cfgplot.xticlabelfmtstr = g_strdup(spl_data->xfmt);
	}
	if (spl_data->yfmt) {
		g_free(spl_data->cfgplot.yticlabelfmtstr);
		spl_data->cfgplot.yticlabelfmtstr = g_strdup(spl_data->yfmt);
	}

	struct kplot *p = kplot_alloc(&spl_data->cfgplot);

	// data plots
	int nb_graphs = 0, nb_xygraphs = 0;
	GList *legend = NULL;
	GString *legend_text = NULL;

	// xylines
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
		d1 = kdata_array_alloc(plot->data, plot->nb);
		enum kplottype plottype = (plot->pl_type == KPLOT_UNDEFINED) ? spl_data->plottype : plot->pl_type;
		kplot_attach_data(p, d1, plottype, &spl_data->cfgdata);
		if (spl_data->show_legend) {
			int index = nb_graphs % spl_data->cfgplot.clrsz;
			legend = g_list_append(legend, new_legend_entry(SIRIL_PLOT_XY, spl_data->cfgplot.clrs[index].rgba));
			if (!nb_graphs)
				legend_text = g_string_new((!plot->label) ? "\n" : plot->label); // in case the first label is empty
			else
				g_string_append_printf(legend_text, "\n%s", plot->label);
		}
		kdata_destroy(d1);
		nb_graphs++;
		nb_xygraphs++;
	}
	d1 = NULL;
	// xy points with y error bars
	for (GList *list = spl_data->plots; list; list = list->next) {
		splxyerrdata *plots = (splxyerrdata *)list->data;
		const struct kdatacfg *cfgs[3];
		enum kplottype plotstypes[3];
		for (int i = 0; i < 3; i++) {
			d2[i] = kdata_array_alloc(plots->plots[i]->data, plots->nb);
			cfgs[i] = &spl_data->cfgdata;
			plotstypes[i] = (plots->plots[i]->pl_type == KPLOT_UNDEFINED) ? spl_data->plotstypes[i] : plots->plots[i]->pl_type;
		}
		kplot_attach_datas(p, 3, d2, plotstypes, cfgs, spl_data->plotstype);
		if (spl_data->show_legend) {
			int index = nb_graphs % spl_data->cfgplot.clrsz;
			legend = g_list_append(legend, new_legend_entry(SIRIL_PLOT_XYERR, spl_data->cfgplot.clrs[index].rgba));
			if (!nb_graphs)
				legend_text = g_string_new((!plots->label) ? "\n" : plots->label); // in case the first label is empty
			else
				g_string_append_printf(legend_text, "\n%s", plots->label);
		}
		for (int i = 0; i < 3; i++) {
			kdata_destroy(d2[i]);
			d2[i] = NULL;
		}
		nb_graphs++;
	}

	// preparing the surfaces
	double drawwidth = width;
	double drawheight = height;
	double top = 0.;

	// painting the whole surface white
	cairo_set_source_rgb(cr, color, color, color);
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_fill(cr);

	// writing the title if any and booking space
	if (spl_data->title) {
		cairo_save(cr);
		PangoLayout *layout;
		PangoFontDescription *desc;
		int pw, ph;
		layout = pango_cairo_create_layout(cr);
		pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);
		pango_layout_set_markup(layout, spl_data->title, -1);
		// set max width to wrap title if required
		pango_layout_set_width(layout, (width - 2 * SIRIL_PLOT_MARGIN) * PANGO_SCALE);
		pango_layout_set_wrap(layout,PANGO_WRAP_WORD);
		desc = pango_font_description_new();
		pango_font_description_set_size(desc, SIRIL_PLOT_TITLE_SIZE * PANGO_SCALE);
		pango_font_description_set_family(desc, SIRIL_PLOT_FONT_FAMILY);
		pango_layout_set_font_description(layout, desc);
		pango_font_description_free(desc);
		pango_layout_get_size(layout, &pw, &ph);
		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		cairo_move_to(cr, (double)SIRIL_PLOT_MARGIN, (double)SIRIL_PLOT_MARGIN);
		pango_cairo_show_layout(cr, layout);
		g_object_unref(layout);
		cairo_restore(cr); // restore the orginal context
		top = (double)SIRIL_PLOT_MARGIN + (double)ph / PANGO_SCALE;
		drawheight = height - top;
	}

	// creating a surface to draw the plot (accounting for title-reserved space if required)
	cairo_surface_t *draw_surface = NULL;
#ifdef CAIRO_HAS_SVG_SURFACE
	if (!for_svg)
		draw_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, (int)drawwidth, (int)drawheight);
	else
		draw_surface = cairo_svg_surface_create_for_stream(NULL, NULL, drawwidth, drawheight);
#else
	draw_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, (int)drawwidth, (int)drawheight);
#endif
	struct kplotctx ctx = { 0 };
	cairo_t *draw_cr = cairo_create(draw_surface);
	kplot_draw(p, drawwidth, drawheight, draw_cr, &ctx);
	spl_data->pdd.range = (point){ ctx.dims.x,  ctx.dims.y};
	spl_data->pdd.offset = (point){ ctx.offs.x,  ctx.offs.y + top};
	if (spl_data->bkg) {
		// TODO later zoomable bkg
		// int offsx, offsy;
		// GdkPixbuf *subbkg = extract_sub_bkg(spl_data, xmin, xmax, ymin, ymax, &offsx, &offsy);
		GdkPixbuf *bkg = gdk_pixbuf_scale_simple(spl_data->bkg->img, (int)ctx.dims.x, (int)ctx.dims.y, GDK_INTERP_BILINEAR);
		gdk_cairo_set_source_pixbuf(cr, bkg, ctx.offs.x, ctx.offs.y + top);
		cairo_paint(cr);
		cairo_fill(cr);
		g_object_unref(bkg);
		// g_object_unref(subbkg);
	}
	cairo_set_source_surface(cr, draw_surface, 0., (int)top);
	cairo_paint(cr);
	cairo_destroy(draw_cr);
	cairo_surface_destroy(draw_surface);
	kplot_free(p);

	if (spl_data->show_legend) {
		// creating the legend
		PangoLayout *layout;
		PangoFontDescription *desc;
		int pw, ph;
		layout = pango_cairo_create_layout(cr);
		desc = pango_font_description_new();
		pango_font_description_set_size(desc, SIRIL_PLOT_LEGEND_SIZE * PANGO_SCALE);
		pango_font_description_set_family(desc, SIRIL_PLOT_FONT_FAMILY);
		pango_layout_set_font_description(layout, desc);
		pango_font_description_free(desc);

		pango_layout_set_markup(layout, legend_text->str, -1);
		pango_layout_get_size(layout, &pw, &ph);
		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		double px0 = spl_data->pdd.offset.x - (double)SIRIL_PLOT_MARGIN + spl_data->pdd.range.x - (double)pw / PANGO_SCALE;
		double py0 = spl_data->pdd.offset.y + (double)SIRIL_PLOT_MARGIN;
		cairo_move_to(cr, px0, py0);
		pango_cairo_show_layout(cr, layout);
		cairo_stroke(cr);

		px0 -= 6. * SIRIL_PLOT_MARGIN;
		PangoLayoutIter *iter = pango_layout_get_iter(layout);
		int y0;
		guint index = 0;
		do {
			GList *current_entry = g_list_nth(legend, index);
			double color3[3];
			memcpy(color3, ((spllegend *)current_entry->data)->color, 3 * sizeof(double));
			if (iter == NULL)
				break;
			y0 = pango_layout_iter_get_baseline(iter);
			double dy = (double)y0 / PANGO_SCALE - 0.5 * SIRIL_PLOT_LEGEND_SIZE;
			cairo_set_source_rgb(cr, color3[0], color3[1], color3[2]);
			cairo_set_line_width(cr, 1.);
			if (index < nb_xygraphs) {
				cairo_move_to(cr, px0, py0 + dy);
				cairo_rel_line_to(cr, 4. * SIRIL_PLOT_MARGIN, 0.);
			} else {
				cairo_arc(cr, px0 + 3. * SIRIL_PLOT_MARGIN, py0 + dy, 0.5 * SIRIL_PLOT_LEGEND_SIZE, 0., 2. * M_PI);
			}
			cairo_stroke(cr);
			index++;
		} while (pango_layout_iter_next_line(iter) && index < nb_graphs);
		g_object_unref(layout);

		// freeing
		g_string_free(legend_text, TRUE);
		g_list_free_full(legend, (GDestroyNotify) free_legend);
	}

	return TRUE;
}

cairo_surface_t *siril_plot_draw_to_image_surface(siril_plot_data *spl_data, int width, int height) {
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
	if (cairo_surface_status(surface)) {
		siril_debug_print("Could not create cairo surface\n");
		return NULL;
	}
	cairo_t *cr = cairo_create(surface);
	if (cairo_status(cr)) {
		cairo_surface_destroy(surface);
		siril_debug_print("Could not create cairo context\n");
		return NULL;
	}
	if (!siril_plot_draw(cr, spl_data, (double)width, (double)height, FALSE)) {
		siril_debug_print("Could not draw to cairo context\n");
		cairo_surface_destroy(surface);
		surface = NULL;
	}
	cairo_destroy(cr);
	return surface;
}
// draw the data contained in spl_data and saves as png file
gboolean siril_plot_save_png(siril_plot_data *spl_data, char *pngfilename, int width, int height) {
	gboolean success = TRUE;
	cairo_surface_t *png_surface = siril_plot_draw_to_image_surface(spl_data, (width) ? width : SIRIL_PLOT_PNG_WIDTH, (height) ? height : SIRIL_PLOT_PNG_HEIGHT);
	if (!png_surface)
		return FALSE;

	siril_debug_print("Successfully created png plot\n");
	if (!cairo_surface_write_to_png(png_surface, pngfilename))
		siril_log_message(_("%s has been saved.\n"), pngfilename);
	else
		success = FALSE;
	cairo_surface_destroy(png_surface);
	return success;
}

#ifdef CAIRO_HAS_SVG_SURFACE
// draw the data contained in spl_data and saves as svg file
gboolean siril_plot_save_svg(siril_plot_data *spl_data, char *svgfilename, int width, int height) {
	gboolean success = TRUE;
	cairo_t *svg_cr = NULL;
	//create the surface
	cairo_surface_t *svg_surface = cairo_svg_surface_create(svgfilename, (width) ? width : SIRIL_PLOT_PNG_WIDTH, (height) ? height : SIRIL_PLOT_PNG_HEIGHT);
	if (cairo_surface_status(svg_surface)) {
		siril_debug_print("Could not create svg surface\n");
		success = FALSE;
	}
	//create the context
	if (success) {
		cairo_svg_surface_set_document_unit(svg_surface, CAIRO_SVG_UNIT_PX);
		svg_cr = cairo_create(svg_surface);
		if (cairo_status(svg_cr)) {
			siril_debug_print("Could not create svg context\n");
			success = FALSE;
		}
	}
	// draw the plot and save the surface to svg
	if (success && siril_plot_draw(svg_cr, spl_data, (width) ? width : SIRIL_PLOT_PNG_WIDTH, (height) ? height : SIRIL_PLOT_PNG_HEIGHT, TRUE))
		siril_log_message(_("%s has been saved.\n"), svgfilename);
	else {
		success = FALSE;
		siril_debug_print("Could not draw to svg context\n");
	}

	if (svg_cr)
		cairo_destroy(svg_cr);
	cairo_surface_destroy(svg_surface);
	return success;
}
#endif

// save the data to a dat file
gboolean siril_plot_save_dat(siril_plot_data *spl_data, const char *datfilename, gboolean add_title) {
	GString *header = NULL;
	FILE* fileout = NULL;
	gboolean retval = TRUE;
	double *data = NULL;
	char* newfilename = NULL;

	// TODO: for the time-being, we will assume that x of the first series
	// is valid for all other series. May need to complexify this in the future
	int nbpoints = 0, nbcols = 1, nbgraphs = 0, j = 0;
	if (add_title && spl_data->title) {
		// spl_data->title is assumed to have the # signs at each line start as necessary
		// and to finish by a \n character
		header = g_string_new(spl_data->title);
		g_string_append_printf(header, "#x");
	} else
		header = g_string_new("#x");
	// xylines
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
		if (nbpoints == 0)
			nbpoints = plot->nb;
		else if (plot->nb != nbpoints) {
			siril_debug_print("Cannot export to *.dat series of different length, skipping\n");
			continue;
		}
		gchar *label = (plot->label) ? g_strdup(plot->label) : g_strdup_printf("Series_%02d", nbgraphs + 1);
		replace_spaces_from_str(label, '_');
		g_string_append_printf(header, " %s", label);
		g_free(label);
		nbgraphs++;
		nbcols++;
	}
	// xy points with y error bars
	for (GList *list = spl_data->plots; list; list = list->next) {
		splxyerrdata *plots = (splxyerrdata *)list->data;
		if (nbpoints == 0)
			nbpoints = plots->nb;
		else if (plots->nb != nbpoints) {
			siril_debug_print("Cannot export to *.dat series of different length, skipping\n");
			continue;
		}
		gchar *label = (plots->label) ? g_strdup(plots->label) : g_strdup_printf("Series_%02d", nbgraphs + 1);
		replace_spaces_from_str(label, '_');
		g_string_append_printf(header, " %s %s_err+ %s_err-", label, label, label);
		g_free(label);
		nbgraphs++;
		nbcols += 3;
	}

	// gathering all the data
	data = calloc(nbpoints, nbcols * sizeof(double)); // Initialize data to avoid possible use of uninitalized var in the loop that writes to the file
	if (!data) {
		PRINT_ALLOC_ERR;
		retval = FALSE;
		goto clean_and_exit;
	}
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
		if (plot->nb != nbpoints)
			continue;
		// adding x if none is present
		if (j == 0) {
			int index = j;
			for (int i = 0; i < nbpoints; i++) {
				data[index] = plot->data[i].x;
				index += nbcols;
			}
			j++;
		}
		// adding y
		int index = j;
		for (int i = 0; i < nbpoints; i++) {
			data[index] = plot->data[i].y;
			index += nbcols;
		}
		j++;
	}
	for (GList *list = spl_data->plots; list; list = list->next) {
		splxyerrdata *plots = (splxyerrdata *)list->data;
		if (plots->nb != nbpoints)
			continue;
		// adding x if none is present
		if (j == 0) {
			int index = j;
			for (int i = 0; i < nbpoints; i++) {
				data[index] = plots->plots[0]->data[i].x;
				index += nbcols;
			}
			j++;
		}
		// adding ys
		int index = j;
		for (int i = 0; i < nbpoints; i++) {
			for (int k = 0; k < 3; k++)
				data[index + k] = plots->plots[k]->data[i].y;
			index += nbcols;
		}
		j += 3;
	}
	newfilename = strdup(datfilename);
	if (!g_str_has_suffix(newfilename, ".dat")) {
		str_append(&newfilename, ".dat");
	}
	fileout = g_fopen(newfilename, "w");
	if (fileout == NULL) {
		siril_log_message(_("Could not create %s, aborting\n"));
		retval = FALSE;
		goto clean_and_exit;
	}
	fprintf(fileout, "%s", header->str);
	int index = 0;
	for (int r = 0 ; r < nbpoints ; r++) {
		fprintf(fileout, "\n%g", data[index++]); // print newline and x
		for (int c = 1 ; c < nbcols ; c++)
			fprintf(fileout, " %g", data[index++]);
	}
	fclose(fileout);
	siril_log_message(_("%s has been saved.\n"), newfilename);

clean_and_exit:
	g_string_free(header, TRUE);
	free(data);
	free(newfilename);
	return retval;
}
