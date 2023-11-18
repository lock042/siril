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
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo/cairo-svg.h>
#endif
#include <pango/pangocairo.h>
#include <math.h>
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/plot.h"


// static variables

#define GUIDE 20 // used to determine number of tics and spacing


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
	splxydata *plot = g_slice_new(splxydata);
	plot->data = calloc(nb, sizeof(struct kpair));
	if (!plot->data) {
		PRINT_ALLOC_ERR;
		g_slice_free(splxydata, plot);
		return NULL;
	}
	plot->nb = nb;
	plot->label = NULL;
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

// allocate a legend entry
static spllegend *new_legend_entry(spl_type type, double color[3]) {
	spllegend *legend = g_slice_new(spllegend);
	legend->type = type;
	memcpy(legend->color, color, 3 * sizeof(double));
	return legend;
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
	spl_data->interactive = FALSE;

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
	free(spl_data);
	spl_data = NULL;
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
	double x1 = (!spl_data->interactive) ? spl_data->datamin.x : spl_data->pdd.datamin.x;
	double x2 = (!spl_data->interactive) ? spl_data->datamax.x : spl_data->pdd.datamax.x;
	double y1 = (!spl_data->interactive) ? spl_data->datamin.y : spl_data->pdd.datamin.y;
	double y2 = (!spl_data->interactive) ? spl_data->datamax.y : spl_data->pdd.datamax.y;
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
		kplot_attach_data(p, d1, spl_data->plottype, &spl_data->cfgdata);
		if (spl_data->show_legend) {
			int index = nb_graphs % spl_data->cfgplot.clrsz;
			legend = g_list_append(legend, new_legend_entry(SIRIL_PLOT_XY, spl_data->cfgplot.clrs[index].rgba));
			if (!nb_graphs)
				legend_text = g_string_new((!plot->label) ? "\n" : plot->label); // in case the first label is empty
			else
				g_string_append_printf(legend_text, "\n%s", plot->label);
		}
		kdata_destroy(d1);
		d1 = NULL;
		nb_graphs++;
		nb_xygraphs++;
	}
	// xy points with y error bars
	for (GList *list = spl_data->plots; list; list = list->next) {
		splxyerrdata *plots = (splxyerrdata *)list->data;
		const struct kdatacfg *cfgs[3];
		for (int i = 0; i < 3; i++) {
			d2[i] = kdata_array_alloc(plots->plots[i]->data, plots->nb);
			cfgs[i] = &spl_data->cfgdata;
		}
		kplot_attach_datas(p, 3, d2, spl_data->plotstypes, cfgs, spl_data->plotstype);
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
	cairo_set_source_surface(cr, draw_surface, 0., (int)top);
	cairo_paint(cr);
	cairo_destroy(draw_cr);
	cairo_surface_destroy(draw_surface);
	kplot_free(p);
	spl_data->pdd.range = (point){ ctx.dims.x,  ctx.dims.y};
	spl_data->pdd.offset = (point){ ctx.offs.x,  ctx.offs.y + top};

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
			double color[3];
			memcpy(color, ((spllegend *)current_entry->data)->color, 3 * sizeof(double));
			if (iter == NULL)
				break;
			y0 = pango_layout_iter_get_baseline(iter);
			double dy = (double)y0 / PANGO_SCALE - 0.5 * SIRIL_PLOT_LEGEND_SIZE;
			cairo_set_source_rgb(cr, color[0], color[1], color[2]);
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
		g_list_free(legend);
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
	if (success && siril_plot_draw(svg_cr, spl_data, (double)SIRIL_PLOT_PNG_WIDTH, (double)SIRIL_PLOT_PNG_HEIGHT, TRUE))
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
	data = malloc(nbpoints * nbcols * sizeof(double));
	if (!data) {
		PRINT_ALLOC_ERR;
		retval = FALSE;
		goto clean_and_exit;
	}
	for (GList *list = spl_data->plot; list; list = list->next) {
		splxydata *plot = (splxydata *)list->data;
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
	fileout = g_fopen(datfilename, "w");
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
	siril_log_message(_("%s has been saved.\n"), datfilename);

clean_and_exit:
	g_string_free(header, TRUE);
	free(data);
	return retval;
}








