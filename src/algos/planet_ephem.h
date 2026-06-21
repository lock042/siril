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

#ifndef SRC_ALGOS_PLANET_EPHEM_H_
#define SRC_ALGOS_PLANET_EPHEM_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Built-in analytic planetary ephemeris for derotation.
 *
 * Self-contained (no network, no kernel files): truncated VSOP87A series for
 * heliocentric positions + the IAU/WGCCRE pole & rotation models. Validated
 * against the IMCCE Miriade web service (see planet_ephem_test): RA, Dec,
 * distance, sub-observer latitude, pole position angle and apparent radius all
 * agree to < 0.01 deg / < 0.03 arcsec at contemporary epochs.
 *
 * The quantities needed for derotation are the per-instant disk geometry:
 * sub-observer latitude (B), central-meridian longitude (CM, one per rotation
 * system) and the pole position angle (P). Derotation uses *differences* of
 * these between a frame's time and the reference epoch, so a constant offset
 * in absolute CM (e.g. from the delta-T model) cancels and does not affect the
 * result. */

typedef enum {
	PLANET_JUPITER = 0,
	PLANET_SATURN  = 1,
	PLANET_MARS    = 2,
	PLANET_BODY_COUNT
} planet_body_t;

/* Rotation system selector for the central-meridian longitude used by
 * derotation. Jupiter defines all three; Saturn exposes its IAU (System III)
 * value in all slots; Mars has a single areographic system in all slots. */
typedef enum {
	ROT_SYSTEM_I   = 0,
	ROT_SYSTEM_II  = 1,
	ROT_SYSTEM_III = 2
} rot_system_t;

typedef struct {
	double ra;            /* astrometric J2000 right ascension, degrees [0,360) */
	double dec;           /* astrometric J2000 declination, degrees [-90,90] */
	double dist_au;       /* geocentric distance, AU */
	double sub_obs_lat;   /* B: planetocentric sub-observer latitude, degrees */
	double cm[3];         /* central-meridian longitude per rotation system, deg */
	double pole_pa;       /* P: position angle of the north pole, degrees */
	double ang_diam_eq;   /* apparent equatorial diameter, arcsec */
	double flattening;    /* geometric flattening f = (Req-Rpol)/Req */
} planet_geom_t;

/* Compute the geometry of `body` at UTC Julian date `jd_utc`.
 *
 * Observer location (geodetic latitude/longitude in degrees, elevation in
 * metres) refines RA/Dec by topocentric parallax; pass NAN for any of them to
 * use the geocentre. Parallax shifts the sub-observer geometry by < 0.001 deg
 * so B/CM/P are effectively observer-independent.
 *
 * Returns 0 on success, non-zero on invalid argument. */
int planet_ephemeris(planet_body_t body, double jd_utc,
                     double obs_lat, double obs_lon, double obs_elev,
                     planet_geom_t *out);

/* TT - UTC (delta-T), seconds, for the given UTC Julian date. Exposed for
 * tests and callers that want the time-scale offset. */
double planet_ephem_delta_t(double jd_utc);

#ifdef __cplusplus
}
#endif

#endif /* SRC_ALGOS_PLANET_EPHEM_H_ */
