/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#ifndef SRC_GUI_PLOT_H_
#define SRC_GUI_PLOT_H_

#include "core/siril.h"

void reset_plot();
void drawPlot();
void notify_new_photometry();
void free_photometry_set(sequence *seq, int set);

/* for one curve */
typedef struct plot_data_struct {
	double *julian, *frame;
	struct kpair *data, *err;
	int nb;
	struct plot_data_struct *next;
} pldata;


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

/* for headless light curve generation */
void generate_magnitude_data(sequence *seq, int dataset, int ref_image, pldata *plot);
int light_curve(pldata *plot, sequence *seq, gchar *filename);

pldata *alloc_plot_data(int size);
void free_plot_data(pldata *plot);

#endif /* SRC_GUI_PLOT_H_ */
