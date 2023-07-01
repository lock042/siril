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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gui/utils.h"

// static variables

#define GUIDE 20

// static functions

static void freexyplot(splxydata *plot) {
	while (plot) {
		splxydata *next = plot->next;
		free(plot->data);
		free(plot);
		plot = next;
	}
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
	plot->next = NULL;
	plot->nextplots = NULL;
	return plot;
}

// allocate 3 xy data structures to hold data, errp and errm
static splxydata *alloc_xyplots_data(int nb) {
	splxydata *plots = alloc_xyplot_data(nb);
	if (!plots)
		return NULL;
	plots->next = alloc_xyplot_data(nb);
	if (!plots->next) {
		freexyplot(plots);
		return NULL;
	}
	plots->next->next = alloc_xyplot_data(nb);
	if (!plots->next->next) {
		freexyplot(plots);
		return NULL;
	}
	return plots;
}

// init/destroy

void init_siril_plot_data(siril_plot_data *spl_data) {
	spl_data->plot = NULL;
	spl_data->plots = NULL;
	spl_data->title = NULL;
	spl_data->xlabel = NULL;
	spl_data->ylabel = NULL;
	spl_data->xfmt = NULL;
	spl_data->yfmt = NULL;
	spl_data->datfilename = NULL;
	spl_data->pngfilename = NULL;
	spl_data->plottype = KPLOT_LINES;
	spl_data->plotstype = KPLOTS_YERRORBAR;
	spl_data->datamin = (point){ DBL_MAX, DBL_MAX};
	spl_data->datamax = (point){ -DBL_MAX, -DBL_MAX};

	// initializing kplotcfg
	kplotcfg_defaults(&spl_data->kplotcfg);
	// set_colors(&spl_data->kplotcfg);
	spl_data->kplotcfg.ticlabel = TICLABEL_LEFT | TICLABEL_BOTTOM;
	spl_data->kplotcfg.border = BORDER_ALL;
	spl_data->kplotcfg.borderline.clr.type = KPLOTCTYPE_RGBA;
	spl_data->kplotcfg.borderline.clr.rgba[0] = 0.5;
	spl_data->kplotcfg.borderline.clr.rgba[1] = 0.5;
	spl_data->kplotcfg.borderline.clr.rgba[2] = 0.5;
	spl_data->kplotcfg.borderline.clr.rgba[3] = 1.0;
	spl_data->kplotcfg.yaxislabelrot = M_PI_2 * 3.0;
}

void clear_siril_plot_data(GtkWidget *widget, gpointer user_data) {
	siril_plot_data *spl_data = (siril_plot_data *)user_data;

	g_free(spl_data->title);
	g_free(spl_data->xlabel);
	g_free(spl_data->ylabel);
	g_free(spl_data->xfmt);
	g_free(spl_data->yfmt);
	g_free(spl_data->datfilename);
	g_free(spl_data->pngfilename);

	// freeing the xy plots
	freexyplot(spl_data->plot);

	// freeing the xyerr plots
	splxydata *plots = spl_data->plots;
	while (plots) {
		splxydata *nexts = plots->nextplots;
		freexyplot(plots);
		plots = nexts;
	}
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

void siril_plot_set_datfilename(siril_plot_data *spl_data, const gchar *datfilename) {
	if (spl_data->datfilename)
		g_free(spl_data->datfilename);
	spl_data->datfilename = g_strdup(datfilename);
}

void siril_plot_set_pngfilename(siril_plot_data *spl_data, const gchar *pngfilename) {
	if (spl_data->pngfilename)
		g_free(spl_data->pngfilename);
	spl_data->pngfilename = g_strdup(pngfilename);
}

// utilities
gboolean siril_plot_autotic(double vmin, double vmax, int *nbtics, double *tmin, double *tmax) {
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
gboolean siril_plot_add_xydata(siril_plot_data *spl_data, size_t nb, double *x, double *y, double *errp, double *errm) {
	// single plot case
	if (!errp) {
		splxydata *plot = spl_data->plot;
		if (!plot) {
			spl_data->plot = alloc_xyplot_data(nb);
			plot = spl_data->plot;
		} else {
			while (plot->next) {
				plot = plot->next;
			}
			plot->next = alloc_xyplot_data(nb);
			plot = plot->next;
		}
		if (!plot) {
			siril_debug_print("Could not allocate plot data\n");
			return FALSE;
		}
		for (int i = 0; i < nb; i++) {
			plot->data[i].x = x[i];
			plot->data[i].y = y[i];
			if (x[i] < spl_data->datamin.x) spl_data->datamin.x = x[i];
			if (y[i] < spl_data->datamin.y) spl_data->datamin.y = y[i];
			if (x[i] > spl_data->datamax.x) spl_data->datamax.x = x[i];
			if (y[i] > spl_data->datamax.y) spl_data->datamax.y = y[i];
		}
		return TRUE;
	}
	// error plot case
	splxydata *plots = spl_data->plots;
	if (!plots) {
		spl_data->plots = alloc_xyplots_data(nb);
		plots = spl_data->plots;
	} else {
		while (plots->nextplots) {
			plots = plots->nextplots;
		}
		plots->nextplots = alloc_xyplots_data(nb);
		plots = plots->nextplots;
	}
	if (!plots) {
		siril_debug_print("Could not allocate plots data\n");
		return FALSE;
	}
	for (int i = 0; i < nb; i++) {
		plots->data[i].x = x[i];
		plots->data[i].y = y[i];
		plots->next->data[i].x = x[i];
		plots->next->data[i].y = errp[i];
		plots->next->next->data[i].x = x[i];
		plots->next->next->data[i].y = errm[i];
		if (x[i] < spl_data->datamin.x) spl_data->datamin.x = x[i];
		if (y[i] < spl_data->datamin.y) spl_data->datamin.y = y[i];
		if (x[i] > spl_data->datamax.x) spl_data->datamax.x = x[i];
		if (y[i] > spl_data->datamax.y) spl_data->datamax.y = y[i];
		if (y[i] < spl_data->datamin.y) spl_data->datamin.y = y[i];
		if (y[i] > spl_data->datamax.y) spl_data->datamax.y = y[i];
		if (y[i] < spl_data->datamin.y) spl_data->datamin.y = y[i];
		if (y[i] > spl_data->datamax.y) spl_data->datamax.y = y[i];
	}
	return TRUE;

}










