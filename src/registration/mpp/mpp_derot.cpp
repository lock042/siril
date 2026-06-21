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

#include <cmath>

#include "mpp_derot.h"
#include "mpp_derot_geom.h"

#define DEG (M_PI / 180.0)

void mpp_derot_frame_map(const mpp_derot_t *d, int frame,
                         int out_w, int out_h,
                         double org_x, double org_y, double scale,
                         cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu) {
	derot_diskfit_t disk;
	disk.cx = (d->cx - org_x) * scale;
	disk.cy = (d->cy - org_y) * scale;
	disk.r_eq = d->r_eq * scale;
	disk.flattening = d->flattening;
	disk.parity = d->parity;

	derot_geom_t epoch;
	epoch.sub_obs_lat = d->epoch_sub_obs_lat * DEG;
	epoch.cm = d->epoch_cm * DEG;
	epoch.pole_angle = d->pole_angle_epoch;

	derot_geom_t fr;
	fr.sub_obs_lat = d->sub_obs_lat[frame] * DEG;
	fr.cm = d->cm[frame] * DEG;
	/* the fitted in-image pole angle plus the slow ephemeris pole-PA drift */
	fr.pole_angle = d->pole_angle_epoch
	              + (d->pole_pa[frame] - d->epoch_pole_pa) * DEG;

	cv::Mat valid;   /* mu already encodes validity (0 off-disk) */
	derot_build_map(out_w, out_h, disk, epoch, fr, mapx, mapy, valid, mu);
}
