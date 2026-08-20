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

#ifndef SRC_ALGOS_ZPN_FIT_H_
#define SRC_ALGOS_ZPN_FIT_H_

#include "core/siril.h"

#define ZPN_MAX_PV 7 // highest PV term index supported (odd radial terms)

typedef struct {
	double scale;             // central plate scale, deg/px (goes into CDELT/CD)
	double pv[ZPN_MAX_PV + 1];// radial polynomial: pv[1]=1, pv[3], pv[5], pv[7]; rest 0
	int max_order;            // highest odd order actually fit
	double rms_px;            // radial residual of the fit, in pixels
	double max_px;            // worst radial residual, in pixels
	gboolean monotonic;       // R(theta) strictly increasing over the field
} zpn_fit_result;

/* Fit the ZPN radial law R(theta_zd) = sum_m pv[m] * theta_zd^m (degrees, odd m,
 * pv[1] fixed to 1 with the plate scale carried separately) from matched
 * pixel<->sky correspondences about a projection center (ra0,dec0) at pixel
 * (crpix0,crpix1). Rotation/parity are irrelevant here (pixel radius and zenith
 * distance are both rotation-invariant). max_order is the highest odd term to
 * fit (3, 5 or 7). Returns 0 on success. */
int fit_zpn_radial(const double *px, const double *py,
		const double *ra, const double *dec, int n,
		double ra0, double dec0, double crpix0, double crpix1,
		int max_order, zpn_fit_result *out);

/* Fill an (uninitialised) wcsprm with a ZPN solution: CTYPE RA---ZPN/DEC--ZPN,
 * the given center (crval, deg) and reference pixel (crpix, FITS 1-based), the
 * linear core cd[2][2] (deg/px, decomposed into CDELT+PC), and the radial PV
 * terms from the fit. The cd scale must match fit->scale (PV_1 = 1 convention).
 * Returns 0 on success (wcsset succeeded). */
int build_zpn_wcs(struct wcsprm *prm, double crval_ra, double crval_dec,
		double crpix0, double crpix1, double cd[2][2],
		const zpn_fit_result *fit);

#endif /* SRC_ALGOS_ZPN_FIT_H_ */
