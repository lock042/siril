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

#include <math.h>

#include <glib.h>

#include "algos/barycentric.h"

#define DEG2RAD (G_PI / 180.0)
#define AU_LIGHT_TIME 499.004783836	// seconds
#define J2000 2451545.0
#define OBLIQUITY_J2000 (23.43928 * DEG2RAD)

/* TAI-UTC steps (IERS Bulletin C), by MJD of the change. To be extended when
 * a new leap second is announced. */
static const struct { int mjd; int tai_utc; } leap_seconds[] = {
	{ 41317, 10 }, { 41499, 11 }, { 41683, 12 }, { 42048, 13 }, { 42413, 14 },
	{ 42778, 15 }, { 43144, 16 }, { 43509, 17 }, { 43874, 18 }, { 44239, 19 },
	{ 44786, 20 }, { 45151, 21 }, { 45516, 22 }, { 46247, 23 }, { 47161, 24 },
	{ 47892, 25 }, { 48257, 26 }, { 48804, 27 }, { 49169, 28 }, { 49534, 29 },
	{ 50083, 30 }, { 50630, 31 }, { 51179, 32 }, { 53736, 33 }, { 54832, 34 },
	{ 56109, 35 }, { 57204, 36 }, { 57754, 37 }
};

/* Keplerian elements and rates per century, J2000 ecliptic, valid 1800-2050
 * (Standish, JPL "Approximate Positions of the Planets", table 1):
 * a (au), e, I, L, long. perihelion, long. node (degrees),
 * and the planet's mass as GM_sun / GM_planet */
static const struct {
	double el[6], rate[6];
	double mass_ratio;
} planets[] = {
	{ { 0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593 },
	  { 0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081 }, 6023600.0 },
	{ { 0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255 },
	  { 0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418 }, 408523.71 },
	{ { 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0 },
	  { 0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0 }, 328900.56 },
	{ { 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891 },
	  { 0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343 }, 3098708.0 },
	{ { 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909 },
	  { -0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106 }, 1047.3486 },
	{ { 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448 },
	  { -0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794 }, 3497.898 },
	{ { 19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503 },
	  { -0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589 }, 22902.98 },
	{ { 30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574 },
	  { 0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.01262724 }, 19412.24 }
};
#define EARTH_MOON_INDEX 2

static double tai_minus_utc(double jd_utc) {
	double mjd = jd_utc - 2400000.5;
	int dat = 10;
	for (size_t i = 0; i < G_N_ELEMENTS(leap_seconds) && mjd >= leap_seconds[i].mjd; i++)
		dat = leap_seconds[i].tai_utc;
	return (double)dat;
}

double julian_utc_to_tdb(double jd_utc) {
	double jd_tt = jd_utc + (tai_minus_utc(jd_utc) + 32.184) / 86400.0;
	// periodic TDB-TT term, below 2 ms
	double g = (357.53 + 0.98560028 * (jd_tt - J2000)) * DEG2RAD;
	return jd_tt + (0.001657 * sin(g) + 0.000014 * sin(2.0 * g)) / 86400.0;
}

/* heliocentric position of a planet, equatorial J2000, au */
static void planet_position(int p, double T, double pos[3]) {
	double el[6];
	for (int k = 0; k < 6; k++)
		el[k] = planets[p].el[k] + planets[p].rate[k] * T;
	double a = el[0], e = el[1], I = el[2] * DEG2RAD;
	double node = el[5] * DEG2RAD;
	double omega = (el[4] - el[5]) * DEG2RAD;
	double M = fmod((el[3] - el[4]) * DEG2RAD, 2.0 * G_PI);

	double E = M + e * sin(M);
	for (int i = 0; i < 10; i++) {
		double dE = (E - e * sin(E) - M) / (1.0 - e * cos(E));
		E -= dE;
		if (fabs(dE) < 1e-12) break;
	}
	double xp = a * (cos(E) - e);
	double yp = a * sqrt(1.0 - e * e) * sin(E);

	double co = cos(omega), so = sin(omega), cn = cos(node), sn = sin(node);
	double ci = cos(I), si = sin(I);
	double x = (co * cn - so * sn * ci) * xp + (-so * cn - co * sn * ci) * yp;
	double y = (co * sn + so * cn * ci) * xp + (-so * sn + co * cn * ci) * yp;
	double z = (so * si) * xp + (co * si) * yp;

	double ce = cos(OBLIQUITY_J2000), se = sin(OBLIQUITY_J2000);
	pos[0] = x;
	pos[1] = y * ce - z * se;
	pos[2] = y * se + z * ce;
}

double julian_utc_to_bjd_tdb(double jd_utc, double ra, double dec) {
	double jd_tdb = julian_utc_to_tdb(jd_utc);
	double T = (jd_tdb - J2000) / 36525.0;

	// barycentre relative to the Sun, the Earth approximated by the Earth-Moon barycentre
	double bary[3] = { 0.0 }, earth[3] = { 0.0 }, total_mass = 1.0;
	for (int p = 0; p < (int)G_N_ELEMENTS(planets); p++) {
		double pos[3], m = 1.0 / planets[p].mass_ratio;
		planet_position(p, T, pos);
		for (int k = 0; k < 3; k++)
			bary[k] += m * pos[k];
		total_mass += m;
		if (p == EARTH_MOON_INDEX)
			for (int k = 0; k < 3; k++)
				earth[k] = pos[k];
	}

	double cd = cos(dec * DEG2RAD);
	double n[3] = { cd * cos(ra * DEG2RAD), cd * sin(ra * DEG2RAD), sin(dec * DEG2RAD) };
	double proj = 0.0;
	for (int k = 0; k < 3; k++)
		proj += (earth[k] - bary[k] / total_mass) * n[k];
	return jd_tdb + proj * AU_LIGHT_TIME / 86400.0;
}
