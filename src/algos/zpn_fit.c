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

/*
 * ZPN (zenithal polynomial) radial-distortion fit for ultra-wide / fisheye
 * plate solving. ZPN models the radial mapping R(theta_zd) = sum_m PV_m *
 * theta_zd^m (R and the zenith distance theta_zd both in degrees), which stays
 * well-behaved across very wide fields where TAN (R = tan theta_zd) diverges.
 *
 * Given matched image<->catalogue stars and a projection center, the zenith
 * distance comes from the sky positions and the projection-plane radius is the
 * plate scale times the pixel radius from the center:  s * rho = R(theta_zd).
 * Fixing PV_1 = 1 (the plate scale s carries the linear term, as in the FITS
 * convention) leaves a linear least-squares problem in (s, PV_3, PV_5, ...).
 * Rotation and parity are irrelevant: pixel radius and zenith distance are both
 * rotation-invariant, so the orientation/CD is determined elsewhere.
 */

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "algos/siril_wcs.h"
#include "zpn_fit.h"

static double gc_dist_deg(double ra1, double dec1, double ra2, double dec2) {
	const double d2r = G_PI / 180.0;
	double dra = (ra2 - ra1) * d2r, ddec = (dec2 - dec1) * d2r;
	double a = sin(ddec * 0.5) * sin(ddec * 0.5)
		+ cos(dec1 * d2r) * cos(dec2 * d2r) * sin(dra * 0.5) * sin(dra * 0.5);
	if (a < 0.0) a = 0.0;
	if (a > 1.0) a = 1.0;
	return 2.0 * atan2(sqrt(a), sqrt(1.0 - a)) / d2r;
}

/* Solve the m x m linear system A x = b in place by Gauss elimination with
 * partial pivoting. Returns 0 on success, 1 if singular. m <= 8. */
static int solve_linear(double A[8][8], double b[8], double x[8], int m) {
	for (int col = 0; col < m; col++) {
		int piv = col;
		double best = fabs(A[col][col]);
		for (int r = col + 1; r < m; r++) {
			if (fabs(A[r][col]) > best) { best = fabs(A[r][col]); piv = r; }
		}
		if (best < 1e-300)
			return 1;
		if (piv != col) {
			for (int c = 0; c < m; c++) { double t = A[col][c]; A[col][c] = A[piv][c]; A[piv][c] = t; }
			double t = b[col]; b[col] = b[piv]; b[piv] = t;
		}
		for (int r = col + 1; r < m; r++) {
			double f = A[r][col] / A[col][col];
			for (int c = col; c < m; c++)
				A[r][c] -= f * A[col][c];
			b[r] -= f * b[col];
		}
	}
	for (int r = m - 1; r >= 0; r--) {
		double s = b[r];
		for (int c = r + 1; c < m; c++)
			s -= A[r][c] * x[c];
		x[r] = s / A[r][r];
	}
	return 0;
}

int fit_zpn_radial(const double *px, const double *py,
		const double *ra, const double *dec, int n,
		double ra0, double dec0, double crpix0, double crpix1,
		int max_order, zpn_fit_result *out) {
	if (!px || !py || !ra || !dec || !out || n < 3)
		return 1;
	if (max_order < 3) max_order = 3;
	if (max_order > ZPN_MAX_PV) max_order = ZPN_MAX_PV;
	if ((max_order & 1) == 0) max_order--; // odd terms only

	/* Choose which odd radial terms to fit from the field extent: a small theta
	 * span cannot constrain high powers (the theta^5 column is ~1e6x the rho
	 * column), so only add PV_5/PV_7 when the field is wide enough. */
	double theta_max = 0.0;
	for (int i = 0; i < n; i++) {
		double th = gc_dist_deg(ra0, dec0, ra[i], dec[i]);
		if (th > theta_max) theta_max = th;
	}
	int orders[4]; // odd radial terms (beyond PV_1) actually fitted
	int n_extra = 0;
	orders[n_extra++] = 3;
	if (theta_max > 25.0 && max_order >= 5) orders[n_extra++] = 5;
	if (theta_max > 50.0 && max_order >= 7) orders[n_extra++] = 7;
	int m = 1 + n_extra; // unknowns: u[0] = s (scale), then the PV terms

	/* model:  s * rho_i  =  theta_i + sum_k PV_(orders[k]) * theta_i^orders[k]
	 * design row d_i = [ rho_i, -th^o1, -th^o2, ... ], target t_i = theta_i */
	double AtA[8][8] = {{0}}, Atb[8] = {0};
	int used = 0;
	for (int i = 0; i < n; i++) {
		double th = gc_dist_deg(ra0, dec0, ra[i], dec[i]); // zenith distance, deg
		double dx = px[i] - crpix0, dy = py[i] - crpix1;
		double rho = sqrt(dx * dx + dy * dy);              // pixel radius
		if (rho < 1e-6)
			continue;                                      // skip the center point
		double d[8];
		d[0] = rho;
		for (int k = 0; k < n_extra; k++)
			d[1 + k] = -pow(th, orders[k]);
		for (int a = 0; a < m; a++) {
			for (int b = 0; b < m; b++)
				AtA[a][b] += d[a] * d[b];
			Atb[a] += d[a] * th;
		}
		used++;
	}
	if (used < m)
		return 1;

	/* Jacobi (column) equilibration of the normal equations to tame the huge
	 * dynamic range between the rho and theta^m columns, then solve and undo. */
	double cs[8];
	for (int a = 0; a < m; a++) {
		cs[a] = sqrt(AtA[a][a]);
		if (!(cs[a] > 0.0))
			return 1;
	}
	double AtAn[8][8] = {{0}}, Atbn[8] = {0};
	for (int a = 0; a < m; a++) {
		for (int b = 0; b < m; b++)
			AtAn[a][b] = AtA[a][b] / (cs[a] * cs[b]);
		Atbn[a] = Atb[a] / cs[a];
	}
	double un[8] = {0}, u[8] = {0};
	if (solve_linear(AtAn, Atbn, un, m))
		return 1;
	for (int a = 0; a < m; a++)
		u[a] = un[a] / cs[a];

	double s = u[0];
	if (!(s > 0.0) || !isfinite(s))
		return 1;

	/* q[] holds the radial coefficients in the natural "degrees" convention
	 * R(deg) = theta_deg + q_3 theta_deg^3 + ... ; q_1 = 1. The internal fit,
	 * residual and monotonicity all use q. */
	double q[ZPN_MAX_PV + 1] = {0};
	q[1] = 1.0;
	for (int k = 0; k < n_extra; k++)
		q[orders[k]] = u[1 + k];

	/* residuals in pixels: rho_model(theta) = R(theta)/s */
	double sse = 0.0, worst = 0.0;
	int cnt = 0;
	for (int i = 0; i < n; i++) {
		double th = gc_dist_deg(ra0, dec0, ra[i], dec[i]);
		double dx = px[i] - crpix0, dy = py[i] - crpix1;
		double rho = sqrt(dx * dx + dy * dy);
		if (rho < 1e-6)
			continue;
		double R = th; // q_1 = 1
		double thp = th * th * th;
		for (int k = 0; k < n_extra; k++) { R += q[orders[k]] * thp; thp *= th * th; }
		double e = fabs(R / s - rho);
		sse += e * e;
		if (e > worst) worst = e;
		cnt++;
	}

	memset(out, 0, sizeof(*out));
	out->scale = s;
	out->max_order = max_order;
	out->rms_px = cnt ? sqrt(sse / cnt) : 0.0;
	out->max_px = worst;

	/* monotonicity over the field: R'(theta) = 1 + 3 q_3 th^2 + 5 q_5 th^4 + ...
	 * must stay positive so the wcslib inverse converges. */
	out->monotonic = TRUE;
	int steps = 64;
	for (int j = 1; j <= steps; j++) {
		double th = theta_max * j / steps;
		double deriv = 1.0; // d/dtheta of q_1*theta
		for (int k = 0; k < n_extra; k++) {
			int ord = orders[k];
			deriv += ord * q[ord] * pow(th, ord - 1);
		}
		if (deriv <= 0.0) { out->monotonic = FALSE; break; }
	}

	/* Convert to FITS/wcslib ZPN convention: wcslib evaluates
	 * R(deg) = r0 * sum PV_m * theta_rad^m with r0 = 180/pi, theta in radians.
	 * Matching the degrees form gives PV_m = q_m * (180/pi)^(m-1); PV_1 = 1. */
	const double rad2deg = 180.0 / G_PI;
	out->pv[1] = 1.0;
	for (int k = 0; k < n_extra; k++) {
		int ord = orders[k];
		out->pv[ord] = q[ord] * pow(rad2deg, ord - 1);
	}
	return 0;
}

int build_zpn_wcs(struct wcsprm *prm, double crval_ra, double crval_dec,
		double crpix0, double crpix1, double cd[2][2],
		const zpn_fit_result *fit) {
	if (!prm || !fit)
		return 1;
	prm->flag = -1;
	if (wcsinit(1, 2, prm, ZPN_MAX_PV + 2, 0, 0))
		return 1;
	prm->equinox = 2000.0;
	prm->crpix[0] = crpix0;
	prm->crpix[1] = crpix1;
	prm->crval[0] = crval_ra;
	prm->crval[1] = crval_dec;
	prm->lonpole = 180.0;
	const char ctype[2][9] = { "RA---ZPN", "DEC--ZPN" };
	for (int i = 0; i < 2; i++) {
		strncpy(prm->cunit[i], "deg", 71);
		strncpy(prm->ctype[i], ctype[i], 71);
	}
	wcs_decompose_cd(prm, cd); // fills cdelt + pc

	/* ZPN projection parameters live on the latitude axis (axis 2) */
	int k = 0;
	prm->pv[k].i = 2; prm->pv[k].m = 0; prm->pv[k].value = 0.0; k++;
	prm->pv[k].i = 2; prm->pv[k].m = 1; prm->pv[k].value = 1.0; k++;
	for (int m = 3; m <= fit->max_order; m += 2) {
		prm->pv[k].i = 2; prm->pv[k].m = m; prm->pv[k].value = fit->pv[m]; k++;
	}
	prm->npv = k;
	return wcsset(prm);
}
