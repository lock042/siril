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

/* Built-in analytic planetary ephemeris — see planet_ephem.h.
 *
 * Pipeline (all validated against IMCCE Miriade, geocentre, IAU 2015 model):
 *   1. VSOP87A heliocentric rectangular X,Y,Z (ecliptic/equinox J2000) for the
 *      Earth and the target body.
 *   2. Geocentric vector with light-time iteration (planet retarded, Earth at t)
 *      — this yields the astrometric J2000 direction Miriade reports.
 *   3. Rotate ecliptic->equatorial by the J2000 mean obliquity for RA/Dec and
 *      the IAU pole frame.
 *   4. IAU/WGCCRE pole (alpha0, delta0) and prime-meridian W(t) -> sub-observer
 *      latitude B, central-meridian longitude CM (per rotation system) and
 *      pole position angle P.
 */

#include <math.h>
#include <string.h>

#include "planet_ephem.h"
#include "planet_ephem_vsop_data.h"

#define DEG       (M_PI / 180.0)
#define AU_KM     1.495978707e8
#define C_AUD     173.1446326847   /* speed of light, AU per day */
#define EPS0_DEG  23.4392911       /* J2000 mean obliquity of the ecliptic */

static double wrap360(double x) {
	x = fmod(x, 360.0);
	if (x < 0.0) x += 360.0;
	return x;
}

/* ---- VSOP87A series evaluation ----------------------------------------- */

/* One rectangular coordinate (AU) of a body at millennia T from J2000. */
static double vsop_coord_eval(const vsop_coord *c, double T) {
	double total = 0.0, Tp = 1.0;
	for (int p = 0; p < 6; ++p) {
		double s = 0.0;
		for (int i = c->off[p]; i < c->off[p + 1]; ++i)
			s += c->t[i].a * cos(c->t[i].b + c->t[i].c * T);
		total += s * Tp;
		Tp *= T;
	}
	return total;
}

static void vsop_xyz(const vsop_body *b, double jd_tdb, double xyz[3]) {
	const double T = (jd_tdb - 2451545.0) / 365250.0;
	for (int k = 0; k < 3; ++k)
		xyz[k] = vsop_coord_eval(&b->c[k], T);
}

/* ---- delta-T (TT - UTC), seconds --------------------------------------- */

/* Contemporary yearly values (IERS observed / IERS Bulletin-A predictions),
 * adequate to ~1 s. Derotation depends only on CM differences, where any
 * delta-T error cancels; absolute CM shifts by 0.01 deg per second of error. */
static const double DELTAT_Y0 = 2000.0;
static const double DELTAT[] = {
	63.83, 64.09, 64.30, 64.47, 64.57, 64.69, 64.85, 65.15, 65.46, 65.78, /*2009*/
	66.07, 66.32, 66.60, 66.91, 67.28, 67.64, 68.10, 68.59, 68.97, 69.22, /*2019*/
	69.36, 69.36, 69.29, 69.20, 69.10, 69.00, 68.90, 68.80, 68.70, 68.60, /*2029*/
	68.50, 68.40                                                          /*2031*/
};

double planet_ephem_delta_t(double jd_utc) {
	const int n = (int) (sizeof(DELTAT) / sizeof(DELTAT[0]));
	const double year = 2000.0 + (jd_utc - 2451545.0) / 365.25;
	double f = year - DELTAT_Y0;
	if (f <= 0.0) return DELTAT[0];
	if (f >= n - 1) {
		/* Espenak & Meeus long-range fit beyond the table (years > ~2031). */
		double t = year - 2000.0;
		return 62.92 + 0.32217 * t + 0.005589 * t * t;
	}
	int i = (int) f;
	double frac = f - i;
	return DELTAT[i] + (DELTAT[i + 1] - DELTAT[i]) * frac;
}

/* ---- body physical data ------------------------------------------------ */

typedef struct {
	const vsop_body *vsop;
	double req_km, flattening;
	double a0, a0_dot;   /* pole RA: deg and deg/century */
	double d0, d0_dot;   /* pole Dec: deg and deg/century */
	int jupiter_periodic;
	double w0[3], wdot[3];  /* prime meridian per system: deg, deg/day */
} body_phys_t;

static const body_phys_t BODIES[PLANET_BODY_COUNT] = {
	[PLANET_JUPITER] = {
		&VSOP_JUPITER, 71492.0, 0.06487439,
		268.056595, -0.006499, 64.495303, 0.002413, 1,
		{ 67.1, 43.3, 284.95 }, { 877.900, 870.270, 870.5360000 }
	},
	[PLANET_SATURN] = {
		&VSOP_SATURN, 60268.0, 0.09796243,
		40.589, -0.036, 83.537, -0.004, 0,
		{ 38.90, 38.90, 38.90 }, { 810.7939024, 810.7939024, 810.7939024 }
	},
	[PLANET_MARS] = {
		&VSOP_MARS, 3396.19, 0.005886007,
		317.68143, -0.1061, 52.88650, -0.0609, 0,
		{ 176.630, 176.630, 176.630 },
		{ 350.89198226, 350.89198226, 350.89198226 }
	},
};

/* Jupiter periodic pole terms (IAU 2015): refine alpha0/delta0 by ~0.002 deg.
 * T in centuries from J2000. */
static void jupiter_pole_periodic(double T, double *da0, double *dd0) {
	const double Ja = (99.360714  + 4850.4046 * T) * DEG;
	const double Jb = (175.895369 + 1191.9605 * T) * DEG;
	const double Jc = (300.323162 +  262.5475 * T) * DEG;
	const double Jd = (114.012305 + 6070.2476 * T) * DEG;
	const double Je = (49.511251  +   64.3000 * T) * DEG;
	*da0 = 0.000117 * sin(Ja) + 0.000938 * sin(Jb) + 0.001432 * sin(Jc)
	     + 0.000030 * sin(Jd) + 0.002150 * sin(Je);
	*dd0 = 0.000050 * cos(Ja) + 0.000404 * cos(Jb) + 0.000617 * cos(Jc)
	     - 0.000013 * cos(Jd) + 0.000926 * cos(Je);
}

/* ---- main ephemeris ---------------------------------------------------- */

int planet_ephemeris(planet_body_t body, double jd_utc,
                     double obs_lat, double obs_lon, double obs_elev,
                     planet_geom_t *out) {
	(void) obs_lat; (void) obs_lon; (void) obs_elev; /* topocentric: TODO refinement */
	if (!out || body < 0 || body >= PLANET_BODY_COUNT)
		return 1;
	const body_phys_t *bp = &BODIES[body];

	const double jd_tdb = jd_utc + planet_ephem_delta_t(jd_utc) / 86400.0;

	/* Geocentric astrometric vector (ecliptic J2000), light-time corrected. */
	double earth[3], plan[3], g[3];
	vsop_xyz(&VSOP_EARTH, jd_tdb, earth);
	vsop_xyz(bp->vsop, jd_tdb, plan);
	for (int k = 0; k < 3; ++k) g[k] = plan[k] - earth[k];
	double dist = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
	double jd_emit = jd_tdb;
	for (int it = 0; it < 3; ++it) {
		jd_emit = jd_tdb - dist / C_AUD;
		vsop_xyz(bp->vsop, jd_emit, plan);
		for (int k = 0; k < 3; ++k) g[k] = plan[k] - earth[k];
		dist = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
	}

	/* Ecliptic -> equatorial (J2000). */
	const double eps = EPS0_DEG * DEG, ce = cos(eps), se = sin(eps);
	const double xe = g[0];
	const double ye = g[1] * ce - g[2] * se;
	const double ze = g[1] * se + g[2] * ce;
	const double ra = atan2(ye, xe);
	const double dec = asin(ze / dist);

	/* IAU pole (equatorial J2000). */
	const double Tcen = (jd_emit - 2451545.0) / 36525.0;
	double a0 = bp->a0 + bp->a0_dot * Tcen;
	double d0 = bp->d0 + bp->d0_dot * Tcen;
	if (bp->jupiter_periodic) {
		double da0, dd0;
		jupiter_pole_periodic(Tcen, &da0, &dd0);
		a0 += da0; d0 += dd0;
	}
	a0 *= DEG; d0 *= DEG;

	const double npole[3] = { cos(d0) * cos(a0), cos(d0) * sin(a0), sin(d0) };
	/* Unit vector planet -> Earth, equatorial J2000. */
	const double ux = -xe / dist, uy = -ye / dist, uz = -ze / dist;

	/* Sub-observer (planetocentric) latitude B. */
	const double sinB = npole[0]*ux + npole[1]*uy + npole[2]*uz;
	out->sub_obs_lat = asin(sinB) / DEG;

	/* Pole position angle P (north pole, from celestial north, eastward),
	 * [0,360). Referred to the J2000 equator; the IMCCE physical-ephemeris
	 * convention differs by the precession angle (~0.1 deg at present), which is
	 * constant over a derotation span and so cancels in CM/P differences. */
	out->pole_pa = wrap360(atan2(cos(d0) * sin(a0 - ra),
	                             sin(d0) * cos(dec) - cos(d0) * sin(dec) * cos(a0 - ra))
	                       / DEG);

	/* Body-fixed frame: z = pole, x = ascending node of the equator on the
	 * J2000 equator (direction alpha0 + 90 deg), y = z x x. Inertial longitude
	 * of the sub-observer point, measured prograde from the node. */
	const double nx = -sin(a0), ny = cos(a0), nz = 0.0;
	const double ey0 = npole[1]*nz - npole[2]*ny;
	const double ey1 = npole[2]*nx - npole[0]*nz;
	const double ey2 = npole[0]*ny - npole[1]*nx;
	const double lon = atan2(ux*ey0 + uy*ey1 + uz*ey2, ux*nx + uy*ny + uz*nz) / DEG;

	const double d_days = jd_emit - 2451545.0;
	for (int s = 0; s < 3; ++s) {
		const double W = bp->w0[s] + bp->wdot[s] * d_days;
		out->cm[s] = wrap360(W - lon);
	}

	out->ra = wrap360(ra / DEG);
	out->dec = dec / DEG;
	out->dist_au = dist;
	out->flattening = bp->flattening;
	out->ang_diam_eq = 2.0 * asin((bp->req_km / AU_KM) / dist) / DEG * 3600.0;
	return 0;
}
