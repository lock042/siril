/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#ifndef SRC_GUI_PLOT_H_
#define SRC_GUI_PLOT_H_

#include "core/siril.h"
#include "kplot.h"

void clear_all_photometry_and_plot();
void reset_plot();
void drawPlot();
void notify_new_photometry();
void init_plot_colors();

/* for one curve */
typedef struct plot_data_struct {
	double *frame;
	struct kpair *data, *err;
	int nb;		// number of points in the plot (= number of images)
	struct plot_data_struct *next;
} pldata;

typedef struct phot_data_struct {
	int frame;
	struct kpair data, err;
} photdata_t;

/* has to be the same as in the glade file */
enum photometry_source {
	FWHM,
	ROUNDNESS,
	AMPLITUDE,
	MAGNITUDE,
	BACKGROUND,
	X_POSITION,
	Y_POSITION,
	SNR
};

/* has to be the same as in the glade file */
enum registration_source {
	r_FWHM,
	r_ROUNDNESS,
	r_WFWHM,
	r_BACKGROUND,
	r_NBSTARS,
	r_X_POSITION,
	r_Y_POSITION,
	r_QUALITY,
	r_FRAME
};

enum marker_type {
	MARKER_NONE = -1,
	MARKER_X_MIN = 0,
	MARKER_X_MAX = 1,
	MARKER_Y_MIN = 2,
	MARKER_Y_MAX = 3
};

enum slider_type {
	SLIDER_NONE = -1,
	SLIDER_X = 0,
	SLIDER_Y = 1
};

enum border_type {
	SELBORDER_NONE = -1,
	SELBORDER_TOP = 0,
	SELBORDER_BOTTOM = 1,
	SELBORDER_LEFT = 2,
	SELBORDER_RIGHT = 3
};

enum selaction_type {
	SELACTION_NONE = -1,
	SELACTION_SELECTING = 0,
	SELACTION_RESIZING = 1,
	SELACTION_MOVING = 2
};

// same as rectangle but double to avoid recasting
typedef struct {
	double x, y, w, h;
} rectangled;

typedef struct plot_draw_data {
	point datamin; // coordinates of the min (x,y) data values in data units
	point datamax; // coordinates of the max (x,y) data values in data units
	point pdatamin; // coordinates of the plotted min (x,y) data values in data units (accounting for the sliders)
	point pdatamax; // coordinates of the plotted max (x,y) data values in data units (accounting for the sliders)
	point range; // coordinates of the extent of (x,y) axes in pixel units
	point scale; // scales on x and y in data unit/pixel
	point offset; // coordinates of the topleft corner (x,y) axes in pixel units
	double surf_w; // x size of the cairosurface in pixel
	double surf_h; // y size of the cairosurface in pixel
	double xrange[2]; // pair between 0 and 1 giving the extent of plotted x values in the datamin,datamax range
	double yrange[2]; // pair between 0 and 1 giving the extent of plotted y values in the datamin,datamax range
	enum marker_type marker_grabbed; // flag telling which slider marker is being grabbed (0 to 3 from X_MIN to Y_MAX)
	enum slider_type slider_grabbed; // flag telling which slider is being grabbed (0 to 3 from SLIDER_X to SLIDER_Y)
	rectangled selection; // area selected in pixel units (x, y, w, h)
	point start; // coordinates x,y of the initial selection point in pixels
	enum selaction_type action; // selection action being performed
	enum border_type border_grabbed; // the border being dragged
	gboolean *selected; // array same length as data array stating if point is within selection
	int nbselected;
} plot_draw_data_t;

gboolean is_inside_grid(double x, double y, plot_draw_data_t *pddstruct);

#endif /* SRC_GUI_PLOT_H_ */
