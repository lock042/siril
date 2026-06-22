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

/* disk fit of `d` placed into a canvas whose pixel (0,0) is original-frame
 * coordinate (org_x, org_y), scaled by `scale`. */
static derot_diskfit_t place_disk(const mpp_derot_t *d,
                                  double org_x, double org_y, double scale) {
	derot_diskfit_t disk;
	disk.cx = (d->cx - org_x) * scale;
	disk.cy = (d->cy - org_y) * scale;
	disk.r_eq = d->r_eq * scale;
	disk.flattening = d->flattening;
	disk.parity = d->parity;
	return disk;
}

void mpp_derot_frame_map_ms(const mpp_derot_t *out_d, const mpp_derot_t *src_d,
                            int frame, int out_w, int out_h,
                            double out_org_x, double out_org_y, double out_scale,
                            double src_org_x, double src_org_y, double src_scale,
                            cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu) {
	/* output globe = reference sequence's epoch disk in the common canvas */
	derot_diskfit_t out_disk = place_disk(out_d, out_org_x, out_org_y, out_scale);
	/* source globe = this sequence's own fitted disk in its source frame */
	derot_diskfit_t src_disk = place_disk(src_d, src_org_x, src_org_y, src_scale);

	/* epoch (output) orientation from the reference sequence */
	derot_geom_t epoch;
	epoch.sub_obs_lat = out_d->epoch_sub_obs_lat * DEG;
	epoch.cm = out_d->epoch_cm * DEG;
	epoch.pole_angle = out_d->pole_angle_epoch;

	/* frame (source) orientation from this sequence: its fitted in-image pole
	 * angle plus its own slow ephemeris pole-PA drift to the frame time */
	derot_geom_t fr;
	fr.sub_obs_lat = src_d->sub_obs_lat[frame] * DEG;
	fr.cm = src_d->cm[frame] * DEG;
	fr.pole_angle = src_d->pole_angle_epoch
	              + (src_d->pole_pa[frame] - src_d->epoch_pole_pa) * DEG;

	cv::Mat valid;   /* mu already encodes validity (0 off-disk) */
	derot_build_map_ms(out_w, out_h, out_disk, src_disk, epoch, fr,
	                   mapx, mapy, valid, mu);
}

void mpp_derot_frame_map(const mpp_derot_t *d, int frame,
                         int out_w, int out_h,
                         double org_x, double org_y, double scale,
                         cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu) {
	/* single-source: output and source globe coincide (same disk, same canvas) */
	mpp_derot_frame_map_ms(d, d, frame, out_w, out_h,
	                       org_x, org_y, scale, org_x, org_y, scale,
	                       mapx, mapy, mu);
}
