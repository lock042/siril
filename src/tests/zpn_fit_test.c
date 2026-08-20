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
 * along with Siril. If not, see <http://www.gnu.org/licenses/ >.
 */

/*
 * Unit test for the ZPN radial fitter (zpn_fit.c). Ground truth is wcslib
 * itself: build a ZPN wcsprm with known PV_3/PV_5, generate synthetic
 * pixel<->sky correspondences across a wide field via wcsp2s, then assert
 * fit_zpn_radial recovers the plate scale and PV terms.
 */

#include <criterion/criterion.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <wcslib.h>

#include "core/siril.h"
#include "algos/zpn_fit.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image (now a pointer)

#define RX 4080
#define RY 3072

static void build_zpn(struct wcsprm *w, double s, double pv3, double pv5) {
	w->flag = -1;
	wcsinit(1, 2, w, 16, 0, 0);
	w->crpix[0] = RX * 0.5 + 0.5;
	w->crpix[1] = RY * 0.5 + 0.5;
	w->crval[0] = 180.0;
	w->crval[1] = 60.0;
	w->cdelt[0] = -s;
	w->cdelt[1] = s;
	w->pc[0] = 1.0; w->pc[1] = 0.0; w->pc[2] = 0.0; w->pc[3] = 1.0;
	strcpy(w->cunit[0], "deg");
	strcpy(w->cunit[1], "deg");
	strcpy(w->ctype[0], "RA---ZPN");
	strcpy(w->ctype[1], "DEC--ZPN");
	w->npv = 4;
	w->pv[0].i = 2; w->pv[0].m = 0; w->pv[0].value = 0.0;
	w->pv[1].i = 2; w->pv[1].m = 1; w->pv[1].value = 1.0;
	w->pv[2].i = 2; w->pv[2].m = 3; w->pv[2].value = pv3;
	w->pv[3].i = 2; w->pv[3].m = 5; w->pv[3].value = pv5;
	int st = wcsset(w);
	cr_assert(st == 0, "wcsset failed: %d", st);
}

Test(zpn_fit, recovers_scale_and_pv) {
	const double s = 0.02, pv3 = 3.0e-5, pv5 = 1.0e-8; // ~73 arcsec/px, ~40 deg field
	struct wcsprm w;
	memset(&w, 0, sizeof w);
	build_zpn(&w, s, pv3, pv5);

	int GX = 50, GY = 38, n = 0;
	double *px = malloc(GX * GY * sizeof(double));
	double *py = malloc(GX * GY * sizeof(double));
	double *ra = malloc(GX * GY * sizeof(double));
	double *dec = malloc(GX * GY * sizeof(double));
	for (int gy = 0; gy < GY; gy++)
		for (int gx = 0; gx < GX; gx++) {
			double pix[2] = { 1.0 + (RX - 1.0) * gx / (GX - 1),
					  1.0 + (RY - 1.0) * gy / (GY - 1) };
			double img[2], world[2], phi, theta; int st;
			if (wcsp2s(&w, 1, 2, pix, img, &phi, &theta, world, &st))
				continue;
			px[n] = pix[0]; py[n] = pix[1];
			ra[n] = world[0]; dec[n] = world[1];
			n++;
		}
	cr_assert(n > 100, "too few synthetic points: %d", n);

	zpn_fit_result r;
	int ret = fit_zpn_radial(px, py, ra, dec, n, w.crval[0], w.crval[1],
			w.crpix[0], w.crpix[1], 5, &r);
	cr_assert(ret == 0, "fit_zpn_radial failed");

	cr_assert(fabs(r.scale - s) / s < 1e-3, "scale off: %.6f vs %.6f", r.scale, s);
	cr_assert(fabs(r.pv[3] - pv3) / pv3 < 0.05, "PV_3 off: %.4e vs %.4e", r.pv[3], pv3);
	cr_assert(fabs(r.pv[5] - pv5) < 5.0e-9, "PV_5 off: %.4e vs %.4e", r.pv[5], pv5);
	cr_assert(r.monotonic, "fit should be monotonic");
	cr_assert(r.rms_px < 0.05, "radial residual too large: %.4f px", r.rms_px);

	free(px); free(py); free(ra); free(dec);
	wcsfree(&w);
}

Test(zpn_fit, build_wcs_reproduces_correspondences) {
	const double s = 0.02, pv3 = 3.0e-5, pv5 = 1.0e-8;
	struct wcsprm w;
	memset(&w, 0, sizeof w);
	build_zpn(&w, s, pv3, pv5);

	int GX = 50, GY = 38, n = 0;
	double *px = malloc(GX * GY * sizeof(double));
	double *py = malloc(GX * GY * sizeof(double));
	double *ra = malloc(GX * GY * sizeof(double));
	double *dec = malloc(GX * GY * sizeof(double));
	for (int gy = 0; gy < GY; gy++)
		for (int gx = 0; gx < GX; gx++) {
			double pix[2] = { 1.0 + (RX - 1.0) * gx / (GX - 1),
					  1.0 + (RY - 1.0) * gy / (GY - 1) };
			double img[2], world[2], phi, theta; int st;
			if (wcsp2s(&w, 1, 2, pix, img, &phi, &theta, world, &st))
				continue;
			px[n] = pix[0]; py[n] = pix[1];
			ra[n] = world[0]; dec[n] = world[1];
			n++;
		}

	zpn_fit_result r;
	int ret = fit_zpn_radial(px, py, ra, dec, n, w.crval[0], w.crval[1],
			w.crpix[0], w.crpix[1], 5, &r);
	cr_assert(ret == 0);

	/* CD = scale * (PC of the truth); for this fixture PC is identity with the
	 * RA axis negated. Use the fitted scale so PV_1 = 1 is consistent. */
	double cd[2][2] = {{ -r.scale, 0.0 }, { 0.0, r.scale }};
	struct wcsprm built;
	memset(&built, 0, sizeof built);
	int bret = build_zpn_wcs(&built, w.crval[0], w.crval[1], w.crpix[0], w.crpix[1], cd, &r);
	cr_assert(bret == 0, "build_zpn_wcs/wcsset failed: %d", bret);

	double maxe = 0.0;
	for (int i = 0; i < n; i++) {
		double world[2] = { ra[i], dec[i] }, img[2], pix2[2], phi, theta; int st;
		if (wcss2p(&built, 1, 2, world, &phi, &theta, img, pix2, &st))
			continue;
		double e = hypot(pix2[0] - px[i], pix2[1] - py[i]);
		if (e > maxe) maxe = e;
	}
	cr_assert(maxe < 0.1, "built ZPN WCS does not reproduce correspondences: %.4f px", maxe);

	free(px); free(py); free(ra); free(dec);
	wcsfree(&w); wcsfree(&built);
}

Test(zpn_fit, tan_like_field_is_near_linear) {
	// A near-TAN field (small cubic) should fit with PV_3 ~ tan expansion and
	// tiny residual; mostly a sanity check that the fitter is stable when the
	// distortion is small.
	const double s = 0.0005, pv3 = 1.0e-6, pv5 = 0.0; // ~1.8 arcsec/px, ~1 deg
	struct wcsprm w;
	memset(&w, 0, sizeof w);
	build_zpn(&w, s, pv3, pv5);

	int GX = 40, GY = 30, n = 0;
	double *px = malloc(GX * GY * sizeof(double));
	double *py = malloc(GX * GY * sizeof(double));
	double *ra = malloc(GX * GY * sizeof(double));
	double *dec = malloc(GX * GY * sizeof(double));
	for (int gy = 0; gy < GY; gy++)
		for (int gx = 0; gx < GX; gx++) {
			double pix[2] = { 1.0 + (RX - 1.0) * gx / (GX - 1),
					  1.0 + (RY - 1.0) * gy / (GY - 1) };
			double img[2], world[2], phi, theta; int st;
			if (wcsp2s(&w, 1, 2, pix, img, &phi, &theta, world, &st))
				continue;
			px[n] = pix[0]; py[n] = pix[1];
			ra[n] = world[0]; dec[n] = world[1];
			n++;
		}
	zpn_fit_result r;
	int ret = fit_zpn_radial(px, py, ra, dec, n, w.crval[0], w.crval[1],
			w.crpix[0], w.crpix[1], 5, &r);
	cr_assert(ret == 0);
	cr_assert(fabs(r.scale - s) / s < 1e-3, "scale off: %.6f vs %.6f", r.scale, s);
	cr_assert(r.rms_px < 0.05, "residual too large: %.4f px", r.rms_px);

	free(px); free(py); free(ra); free(dec);
	wcsfree(&w);
}
