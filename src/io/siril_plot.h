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

#ifndef SRC_IO_SIRIL_PLOT_H_
#define SRC_IO_SIRIL_PLOT_H_

#include "core/siril.h"
#include "gui/plot.h"
#include "kplot.h"

#define SIRIL_PLOT_WIDTH 800
#define SIRIL_PLOT_HEIGHT 600
#define SPL_TITLE_RATIO 0.1 // veritical proportion of the drawing area booked for title (if any)

typedef struct siril_plot_xydata_struct {
	struct kpair *data; // data
	int nb;		// number of points in the plot
	struct siril_plot_xydata_struct *next; // pointer to the next plot
	struct siril_plot_xydata_struct *nextplots; // pointer to the next plots (used for xyerr plots only)
} splxydata;

typedef struct siril_plot_data_struct {
	splxydata *plot; // a splxydata structure to hold simple data plot (only data)
	splxydata *plots; // a splxydata structure to hold data plots (data and errors)
	gchar *title; // title
	gchar *xlabel; //xlabel
	gchar *ylabel; //ylabel
	gchar *xfmt; // x axis number formatting
	gchar *yfmt; // y axis number formatting
	gchar *datfilename; // name of the datfile to be saved
	gchar *pngfilename; // name of the pngfile to be saved
	enum kplottype plottype;
	enum kplotstype plotstype;
	struct kplotcfg kplotcfg;
	point datamin; // min x/y of data
	point datamax; // max x/y of data
} siril_plot_data;

void init_siril_plot_data(siril_plot_data *spl_data);
void clear_siril_plot_data(GtkWidget *widget, gpointer user_data);

void siril_plot_set_title(siril_plot_data *spl_data, const gchar *title);
void siril_plot_set_xlabel(siril_plot_data *spl_data, const gchar *xlabel);
void siril_plot_set_ylabel(siril_plot_data *spl_data, const gchar *ylabel);
void siril_plot_set_xfmt(siril_plot_data *spl_data, const gchar *xfmt);
void siril_plot_set_yfmt(siril_plot_data *spl_data, const gchar *yfmt);
void siril_plot_set_datfilename(siril_plot_data *spl_data, const gchar *datfilename);
void siril_plot_set_pngfilename(siril_plot_data *spl_data, const gchar *pngfilename);

gboolean siril_plot_autotic(double vmin, double vmax, int *nbtics, double *tmin, double *tmax);
gboolean siril_plot_add_xydata(siril_plot_data *spl_data, size_t nb, double *x, double *y, double *errp, double *errm);

#endif /* SRC_IO_PLOT_H_ */