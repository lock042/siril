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
#ifndef SRC_ALGOS_CCD_INSPECTOR_H_
#define SRC_ALGOS_CCD_INSPECTOR_H_

struct tilt_struct {
	point pt[4];
	double fwhm[4];
	double fwhm_centre;
};

struct tilt_data {
	fits *fit;
	sequence *seq;
	gboolean draw_polygon;
	int nbstars;
	float m, m1, m2, m3, m4, mr1, mr2;
};

void clear_sensor_tilt();
int draw_sensor_tilt(fits *fit);

void apply_tilt_to_sequence(struct tilt_data *tilt_args);

void compute_aberration_inspector();
void redraw_aberration_inspector();

#endif /* SRC_ALGOS_CCD_INSPECTOR_H_ */
