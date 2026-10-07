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

#ifndef SRC_IO_SIRIL_PLOT_H_
#define SRC_IO_SIRIL_PLOT_H_

#include "core/siril.h"
#include "core/gtk_forward_decls.h"
#include "gui-gtk4/plot.h"
#include "kplot.h"
#include "extern.h"

#define SIRIL_PLOT_DISPLAY_WIDTH 600
#define SIRIL_PLOT_DISPLAY_HEIGHT 400
#define SIRIL_PLOT_PNG_WIDTH 800
#define SIRIL_PLOT_PNG_HEIGHT 600
#define SIRIL_PLOT_TITLE_SIZE 12
#define SIRIL_PLOT_MARGIN 5
#define SIRIL_PLOT_LEGEND_SIZE 10.
#define SIRIL_PLOT_FONT_FAMILY "Sans"

typedef enum {
	SIRIL_PLOT_XY,
	SIRIL_PLOT_XYERR
} spl_type;
typedef struct siril_plot_xydata_struct {
	struct kpair *data; // data
	int nb;		// number of points in the plot
	gchar *label; // the name of the series
	double x_offset; // the offset to be added in x
	enum kplottype pl_type; // the plot type (-1 if undefined)
} splxydata;

typedef struct siril_plot_xyerrdata_struct {
	splxydata *plots[3]; // data
	int nb;		// number of points in the plot
	gchar *label; // the name of the series
} splxyerrdata;

typedef struct siril_plot_legend_struct {
	spl_type type;
	double color[3];
} spllegend;

// a key figure displayed in the strip below a plot (label "R²", value "0.76")
typedef struct siril_plot_metric_struct {
	gchar *label;
	gchar *value;
} splmetric;

// a tile of the summary strip displayed above a group of plots
typedef struct siril_plot_tile_struct {
	gchar *icon;     // themed icon name, may be NULL
	gchar *label;
	gchar *value;
	gchar *sublabel; // secondary line, may be NULL
	gchar *subvalue;
} spltile;

// what siril_plot_draw_with() renders on top of the plot itself. The GUI draws
// the header as widgets around the canvas and asks for the canvas alone;
// exports ask for the header too, so that a saved plot names what it shows
typedef enum {
	SPL_DRAW_CANVAS = 0,
	SPL_DRAW_HEADER = 1 << 0 // caption, title and subtitle
} spl_draw_flags;

typedef struct siril_plot_bkg_struct {
	gchar *bkgfilepath;
	cairo_surface_t *img;
	double width, height;
} splbkg;

typedef struct siril_plot_data_struct {
	GList *plot; // a list of splxydata structures to hold simple data plot (only data)
	GList *plots; // a list of splxydata structures to hold data plots (data and errors)
	gchar *title; // title
	gchar *caption; // optional header (markup) shared by several plots, drawn above the title
	gchar *subtitle; // optional detail lines (markup) below the title
	GList *metrics; // a list of splmetric structures, the key figures of the plot
	gboolean logy; // y axis displayed in log10 scale
	gchar *xlabel; //xlabel
	gchar *ylabel; //ylabel
	gchar *xfmt; // x axis number format
	gchar *yfmt; // y axis number format
	gchar *savename; // chain to be prepended when interactively saving (png, dat)
	gboolean forsequence; // using for saving
	enum kplottype plottype;
	enum kplotstype plotstype;
	enum kplottype plotstypes[3];
	struct kplotcfg cfgplot;
	struct kdatacfg cfgdata;
	point datamin; // min x/y of data
	point datamax; // max x/y of data
	gboolean autotic;
	gboolean show_legend;
	gboolean revertX;
	gboolean revertY;
	plot_draw_data_t pdd; // data for display interaction
	gboolean zoomable; // true if GUI display (unless there is a backgroung set)
	splbkg *bkg; // background image structure
	int width; // width of the display area (if 0, uses SIRIL_PLOT_DISPLAY_WIDTH)
	int height; // height of the display area (if 0, uses SIRIL_PLOT_DISPLAY_HEIGHT)
} siril_plot_data;

siril_plot_data* init_siril_plot_data();
void free_siril_plot_data(siril_plot_data *spl_data);
void siril_plot_sort_x(siril_plot_data *spl_data);

// a group of siril_plot_data to be displayed together, side by side, in a single window
typedef struct siril_plot_group_struct {
	GList *items; // list of siril_plot_data* to be displayed together
	gchar *title; // optional title for the window holding the group
	GList *tiles; // a list of spltile structures, summarizing the whole group
} siril_plot_group;

siril_plot_group *siril_plot_group_new();
void siril_plot_group_add(siril_plot_group *grp, siril_plot_data *spl_data);
void siril_plot_group_set_title(siril_plot_group *grp, const gchar *title);
void siril_plot_group_add_tile(siril_plot_group *grp, const gchar *icon, const gchar *label,
		const gchar *value, const gchar *sublabel, const gchar *subvalue);
void siril_plot_free_tiles(GList *tiles);
void free_siril_plot_group(siril_plot_group *grp);

void siril_plot_set_title(siril_plot_data *spl_data, const gchar *title);
void siril_plot_set_caption(siril_plot_data *spl_data, const gchar *caption);
void siril_plot_set_subtitle(siril_plot_data *spl_data, const gchar *subtitle);
void siril_plot_add_metric(siril_plot_data *spl_data, const gchar *label, const gchar *value);
// joins the metrics into one line, e.g. for the clipboard or a plot header
gchar *siril_plot_metrics_to_string(siril_plot_data *spl_data, const gchar *separator);
// log scale only makes sense for strictly positive data, and is not supported
// for series carrying error bars
gboolean siril_plot_can_logscale(siril_plot_data *spl_data);
// switches the y axis to/from log10, converting the stored bounds so that
// zoom and pan keep working in the displayed space
void siril_plot_set_logscale(siril_plot_data *spl_data, gboolean logy);
void siril_plot_set_xlabel(siril_plot_data *spl_data, const gchar *xlabel);
void siril_plot_set_ylabel(siril_plot_data *spl_data, const gchar *ylabel);
void siril_plot_set_xfmt(siril_plot_data *spl_data, const gchar *xfmt);
void siril_plot_set_yfmt(siril_plot_data *spl_data, const gchar *yfmt);
void siril_plot_set_savename(siril_plot_data *spl_data, const gchar *savename);
void siril_plot_set_nth_color(siril_plot_data *spl_data, int n, double color[3]);
void siril_plot_set_nth_plot_type(siril_plot_data *spl_data, int n, enum kplottype pl_type);
void siril_plot_set_nth_plots_types(siril_plot_data *spl_data, int n, enum kplottype pl_type[3]);
gboolean siril_plot_set_background(siril_plot_data *spl_data, const gchar *bkgfilepath);

gboolean siril_plot_add_xydata(siril_plot_data *spl_data, const gchar *label, size_t nb, const double *x, const double *y, const double *errp, const double *errm);
gboolean siril_plot_draw(cairo_t *cr, siril_plot_data *spl_data, double width, double height, gboolean for_svg);
gboolean siril_plot_draw_with(cairo_t *cr, siril_plot_data *spl_data, double width, double height, gboolean for_svg, spl_draw_flags flags);
cairo_surface_t *siril_plot_draw_to_image_surface(siril_plot_data *spl_data, int width, int height);
gboolean siril_plot_save_png(siril_plot_data *spl_data, char *pngfilename, int width, int height);
gboolean siril_plot_save_svg(siril_plot_data *spl_data, char *svgfilename, int width, int height);
gboolean siril_plot_save_dat(siril_plot_data *spl_data, const char *datfilename, gboolean add_title);

#endif /* SRC_IO_PLOT_H_ */
