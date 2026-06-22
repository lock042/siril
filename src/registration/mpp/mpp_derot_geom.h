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

#ifndef SRC_REGISTRATION_MPP_MPP_DEROT_GEOM_H_
#define SRC_REGISTRATION_MPP_MPP_DEROT_GEOM_H_

#include <opencv2/core.hpp>

/* Planetary derotation geometry.
 *
 * Models the planet as an oblate spheroid viewed orthographically (the planet
 * is effectively at infinity at these angular scales). All projection is done
 * in a normalised "north-up" frame where the equatorial radius is 1; the pixel
 * mapping (centre, scale, pole position angle, mirror parity) is a separable
 * affine applied at the boundary.
 *
 * Derotation maps a pixel of the reference-epoch canvas to the pixel of a
 * frame's (MPP-aligned) image showing the same surface point: unproject with
 * the epoch geometry -> rotate the body in longitude (carried implicitly by the
 * differing central-meridian longitudes) -> reproject with the frame geometry.
 * Since MPP global+AP alignment already co-registers each frame's disk to the
 * reference (same centre/scale), only the disk *orientation* (B, CM, pole PA)
 * differs between epoch and frame here; translation/scale ride along in the
 * warp engine's other map terms.
 */

/* Disk fit in canvas pixels (same for epoch input and frame output, because the
 * frames are MPP-aligned to the reference). */
typedef struct {
	double cx, cy;       /* planet centre, canvas pixels */
	double r_eq;         /* apparent equatorial radius, pixels */
	double flattening;   /* geometric flattening f = (Req-Rpol)/Req */
	double parity;       /* +1 normal, -1 mirrored (star-diagonal etc.) */
} derot_diskfit_t;

/* Per-instant disk orientation (radians). */
typedef struct {
	double sub_obs_lat;  /* B: planetocentric sub-observer latitude */
	double cm;           /* central-meridian longitude */
	double pole_angle;   /* theta: position angle of projected north in image */
} derot_geom_t;

/* Forward projection of a body surface point (planetographic latitude `lat_g`,
 * longitude `lon`, radians) to normalised north-up coordinates (u right,
 * v up; equatorial radii). Returns true and fills *u,*v,*mu (cosine of the
 * emission angle, in (0,1]) when the point is on the visible hemisphere; false
 * otherwise. `f` is the flattening. */
bool derot_project(double lat_g, double lon, double B, double CM, double f,
                   double *u, double *v, double *mu);

/* Inverse projection: normalised north-up (u,v) -> visible-side surface point
 * (*lat_g,*lon, radians). Returns false when (u,v) falls outside the limb. */
bool derot_unproject(double u, double v, double B, double CM, double f,
                     double *lat_g, double *lon);

/* Build the dense derotation maps for one frame, sized out_h x out_w (CV_32F
 * for mapx/mapy/mu, CV_8U for valid). Only the globe is derotated; pixels
 * outside the fitted globe (rings/sky) pass through as identity. For each
 * output (epoch-canvas) pixel:
 *   - mapx/mapy: source pixel to sample — the derotated globe point, the
 *                identity outside the globe, or -1 where masked.
 *   - valid:     1 where the map is usable; 0 only where a globe point rotated
 *                behind the limb (no source data).
 *   - mu:        1 where usable, 0 where masked — a hard visible/hidden mask
 *                the warp engine folds into its weight (no brightness taper).
 */
void derot_build_map(int out_w, int out_h,
                     const derot_diskfit_t &disk,
                     const derot_geom_t &epoch,
                     const derot_geom_t &frame,
                     cv::Mat &mapx, cv::Mat &mapy,
                     cv::Mat &valid, cv::Mat &mu);

/* Multi-source generalisation: the output (epoch) globe and the source (frame)
 * globe may sit at different places. `out_disk` (the reference sequence's globe
 * in the output canvas) drives the epoch unprojection; `src_disk` (this
 * sequence's own fitted globe, in its source frame) drives the reprojection.
 * The single resample therefore both derotates the globe to the epoch AND
 * relocates/rescales/re-orients it from this sequence's frame into the common
 * reference canvas. With out_disk == src_disk this is identical to the
 * single-disk derot_build_map above. The flattening (same body) is taken from
 * out_disk. Outside the output globe, pixels pass through as identity (rings/
 * sky of the reference canvas stack normally). */
void derot_build_map_ms(int out_w, int out_h,
                        const derot_diskfit_t &out_disk,
                        const derot_diskfit_t &src_disk,
                        const derot_geom_t &epoch,
                        const derot_geom_t &frame,
                        cv::Mat &mapx, cv::Mat &mapy,
                        cv::Mat &valid, cv::Mat &mu);

#endif /* SRC_REGISTRATION_MPP_MPP_DEROT_GEOM_H_ */
