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

/* Oblate-spheroid orthographic projection for planetary derotation — see
 * mpp_derot_geom.h. Closed-form both ways:
 *
 *   Screen basis from the observer direction O = (cosB cosCM, cosB sinCM, sinB):
 *     r = (-sinCM,  cosCM, 0)               (image +u, "east-ish")
 *     u = (-sinB cosCM, -sinB sinCM, cosB)  (projected north)
 *     O                                     (toward observer, +depth)
 *
 *   Forward: surface point S(lat_g,lon) on the ellipsoid (a=1, c=1-f),
 *     normal n(lat_g,lon); project u=S.r, v=S.u, emission cos mu=n.O.
 *   Inverse: ray (u r + v u + d O) intersected with the ellipsoid, near root.
 */

#include <cmath>

#include "mpp_derot_geom.h"

bool derot_project(double lat_g, double lon, double B, double CM, double f,
                   double *u, double *v, double *mu) {
	const double cf = 1.0 - f;
	const double c2 = cf * cf;           /* = 1 - e^2 */
	const double e2 = 1.0 - c2;
	const double sphi = sin(lat_g), cphi = cos(lat_g);
	const double dlam = lon - CM;
	const double cB = cos(B), sB = sin(B);

	const double m = cphi * cB * cos(dlam) + sphi * sB;  /* n . O */
	if (m <= 0.0)
		return false;                    /* far side / on the limb */

	const double N = 1.0 / sqrt(1.0 - e2 * sphi * sphi);
	*u  = N * cphi * sin(dlam);
	*v  = -N * cphi * sB * cos(dlam) + N * c2 * sphi * cB;
	*mu = m;
	return true;
}

bool derot_unproject(double u, double v, double B, double CM, double f,
                     double *lat_g, double *lon) {
	const double cf = 1.0 - f;
	const double c2 = cf * cf;
	const double cB = cos(B), sB = sin(B), cC = cos(CM), sC = sin(CM);

	/* screen basis */
	const double rx = -sC,        ry = cC,         rz = 0.0;
	const double ux = -sB * cC,   uy = -sB * sC,   uz = cB;
	const double ox = cB * cC,    oy = cB * sC,    oz = sB;

	/* R = u*r + v*u (the d-independent part of the line of sight) */
	const double Rx = u * rx + v * ux;
	const double Ry = u * ry + v * uy;
	const double Rz = u * rz + v * uz;

	/* ellipsoid X^2 + Y^2 + Z^2/c2 = 1, point = R + d*O, solve for d */
	const double A = ox * ox + oy * oy + oz * oz / c2;
	const double Bq = Rx * ox + Ry * oy + Rz * oz / c2;
	const double Cq = Rx * Rx + Ry * Ry + Rz * Rz / c2 - 1.0;
	const double disc = Bq * Bq - A * Cq;
	if (disc < 0.0)
		return false;                    /* outside the limb */
	const double d = (-Bq + sqrt(disc)) / A;   /* near side: larger depth */

	const double Px = Rx + d * ox;
	const double Py = Ry + d * oy;
	const double Pz = Rz + d * oz;

	*lon = atan2(Py, Px);
	const double nz = Pz / c2;            /* surface normal (un-normalised) */
	*lat_g = atan2(nz, sqrt(Px * Px + Py * Py));
	return true;
}

/* normalised north-up (u,v) -> canvas pixel, applying mirror, pole rotation,
 * scale and centre (image rows increase downward, so v subtracts). */
static inline void norm_to_pixel(double u, double v, const derot_diskfit_t &d,
                                 double theta, double *ix, double *iy) {
	const double up = d.parity * u;
	const double ct = cos(theta), st = sin(theta);
	const double u2 = up * ct - v * st;
	const double v2 = up * st + v * ct;
	*ix = d.cx + d.r_eq * u2;
	*iy = d.cy - d.r_eq * v2;
}

static inline void pixel_to_norm(double ix, double iy, const derot_diskfit_t &d,
                                 double theta, double *u, double *v) {
	const double u2 = (ix - d.cx) / d.r_eq;
	const double v2 = -(iy - d.cy) / d.r_eq;
	const double ct = cos(theta), st = sin(theta);
	const double up =  u2 * ct + v2 * st;
	const double vv = -u2 * st + v2 * ct;
	*u = up * d.parity;                  /* parity is +/-1 -> own inverse */
	*v = vv;
}

void derot_build_map_ms(int out_w, int out_h,
                        const derot_diskfit_t &out_disk,
                        const derot_diskfit_t &src_disk,
                        const derot_geom_t &epoch,
                        const derot_geom_t &frame,
                        cv::Mat &mapx, cv::Mat &mapy,
                        cv::Mat &valid, cv::Mat &mu) {
	mapx.create(out_h, out_w, CV_32F);
	mapy.create(out_h, out_w, CV_32F);
	valid.create(out_h, out_w, CV_8U);
	mu.create(out_h, out_w, CV_32F);

	const double f = out_disk.flattening;   /* same body on both sides */

	for (int y = 0; y < out_h; ++y) {
		float *mx = mapx.ptr<float>(y);
		float *my = mapy.ptr<float>(y);
		uint8_t *vp = valid.ptr<uint8_t>(y);
		float *mp = mu.ptr<float>(y);
		for (int x = 0; x < out_w; ++x) {
			double u0, v0, lat, lon, uf, vf, m;
			pixel_to_norm(x, y, out_disk, epoch.pole_angle, &u0, &v0);
			if (!derot_unproject(u0, v0, epoch.sub_obs_lat, epoch.cm,
			                     f, &lat, &lon)) {
				/* Outside the fitted globe (rings, moons, sky): pass the pixel
				 * through unchanged so it is aligned/stacked normally — only the
				 * globe is derotated. */
				mx[x] = (float) x;
				my[x] = (float) y;
				vp[x] = 1;
				mp[x] = 1.0f;
				continue;
			}
			if (!derot_project(lat, lon, frame.sub_obs_lat, frame.cm,
			                   f, &uf, &vf, &m)) {
				/* On the globe but rotated behind the limb at this frame's time:
				 * no source data, so mask it out. */
				mx[x] = -1.0f;
				my[x] = -1.0f;
				vp[x] = 0;
				mp[x] = 0.0f;
				continue;
			}
			/* Visible globe point. Full weight: the emission-angle (mu = m)
			 * down-weighting darkened the limb via the background blend, so the
			 * map carries only a hard visible/hidden mask, not a brightness
			 * taper. */
			double sx, sy;
			norm_to_pixel(uf, vf, src_disk, frame.pole_angle, &sx, &sy);
			mx[x] = (float) sx;
			my[x] = (float) sy;
			vp[x] = 1;
			mp[x] = 1.0f;
			(void) m;
		}
	}
}

void derot_build_map(int out_w, int out_h,
                     const derot_diskfit_t &disk,
                     const derot_geom_t &epoch,
                     const derot_geom_t &frame,
                     cv::Mat &mapx, cv::Mat &mapy,
                     cv::Mat &valid, cv::Mat &mu) {
	/* single-source: the output and source globe coincide */
	derot_build_map_ms(out_w, out_h, disk, disk, epoch, frame,
	                   mapx, mapy, valid, mu);
}
