/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include <glib.h>
#include <math.h>
#include <stdio.h>

#include "siril_world_cs.h"

struct _SirilWorldCS {
	gdouble alpha;

	gdouble delta;

	gint ref_count; /* (atomic) */
};

static SirilWorldCS* siril_world_cs_alloc() {
	SirilWorldCS *world_cs;

	world_cs = g_slice_new0(SirilWorldCS);
	world_cs->ref_count = 1;

	return world_cs;
}


SirilWorldCS *siril_world_cs_ref(SirilWorldCS *world_cs) {
	g_return_val_if_fail(world_cs != NULL, NULL);
	g_return_val_if_fail(world_cs->ref_count > 0, NULL);

	g_atomic_int_inc(&world_cs->ref_count);

	return world_cs;
}

void siril_world_cs_unref(SirilWorldCS *world_cs) {
	g_return_if_fail(world_cs != NULL);
	g_return_if_fail(world_cs->ref_count > 0);

	if (g_atomic_int_dec_and_test(&world_cs->ref_count)) {
		g_slice_free(SirilWorldCS, world_cs);
	}
}

SirilWorldCS* siril_world_cs_new_from_a_d(gdouble alpha, gdouble delta) {
	SirilWorldCS *world_cs;
	/*
	 * make sure new RA lies in range  0 < RA < 360
	 */
	if (alpha < 0) {
		alpha += 360.0;
	}
	if (alpha >= 360.0) {
		alpha -= 360.0;
	}

	/*
	 * make sure Dec lies in range  -90 < Dec < +90
	 */
	if (delta < -90) {
		delta += 180;
	}
	if (delta > 90) {
		delta -= 180;
	}

	world_cs = siril_world_cs_alloc();
	world_cs->alpha = alpha;
	world_cs->delta = delta;

	return world_cs;
}

SirilWorldCS* siril_world_cs_new_from_ra_dec(gdouble ra_h, gdouble ra_m, gdouble ra_s, gdouble dec_deg, gdouble dec_m, gdouble dec_s) {
	SirilWorldCS *world_cs;

	world_cs = siril_world_cs_alloc();

	world_cs->alpha = ra_h * 15.0 + ra_m * 15.0 / 60.0 + ra_s * 15.0 / 3600.0;
	if (dec_deg > 0) {
		world_cs->delta = ((dec_s / 3600.0) + (dec_m / 60.0) + dec_deg);
	} else {
		world_cs->delta = (-(dec_s / 3600.0) - (dec_m / 60.0) + dec_deg);
	}
	return world_cs;
}

/* parses RA and DEC as strings to create a new world_cs object
 * Format of the strings can be decimal values in degrees or
 * 'hours minutes seconds' for RA and 'degrees minutes seconds' for DEC or 
 * 'hours:minutes:seconds' for RA and 'degrees:minutes:seconds' for DEC
 */
SirilWorldCS* siril_world_cs_new_from_objct_ra_dec(gchar *objctra, gchar *objctdec) {
	int ra_h, ra_m, dec_deg, dec_m;
	gdouble ra_s, dec_s;
	gboolean south;

	if (!objctra || objctra[0] == '\0' || !objctdec || objctdec[0] == '\0')
		return NULL;
	gchar *end;
	double ra = g_ascii_strtod(objctra, &end);
	if (end - objctra != strlen(objctra)) {
		ra = NAN;
		if (sscanf(objctra, "%d %d %lf", &ra_h, &ra_m, &ra_s) == 3 ||
				sscanf(objctra, "%d:%d:%lf", &ra_h, &ra_m, &ra_s) == 3)
			ra = ra_h * 15.0 + ra_m * 15.0 / 60.0 + ra_s * 15.0 / 3600.0;
	}

	south = objctdec[0] == '-';
	double dec = g_ascii_strtod(objctdec, &end);
	if (end - objctdec != strlen(objctdec)) {
		dec = NAN;
		if (sscanf(objctdec, "%d %d %lf", &dec_deg, &dec_m, &dec_s) == 3 ||
				sscanf(objctdec, "%d:%d:%lf", &dec_deg, &dec_m, &dec_s) == 3) {
			if ((dec_deg == 0 && !south) || dec_deg > 0)
				dec = ((dec_s / 3600.0) + (dec_m / 60.0) + dec_deg);
			else dec = (-(dec_s / 3600.0) - (dec_m / 60.0) + dec_deg);
		}
	}

	if (isnan(ra) || isnan(dec))
		return NULL;
	return siril_world_cs_new_from_a_d(ra, dec);
}

gdouble siril_world_cs_get_alpha(SirilWorldCS *world_cs) {
	return world_cs->alpha;
}

gdouble siril_world_cs_get_delta(SirilWorldCS *world_cs) {
	return world_cs->delta;
}

gchar* siril_world_cs_delta_format(SirilWorldCS *world_cs, const gchar *format) {
	g_return_val_if_fail(world_cs != NULL, NULL);
	g_return_val_if_fail(format != NULL, NULL);
	g_return_val_if_fail(g_utf8_validate (format, -1, NULL), NULL);

	gchar sig = '+';
	gdouble dec = world_cs->delta;

	if (dec < 0) sig = '-';

	dec = fabs(dec);

	int degree = (int) dec;
	int min = abs((int) ((dec - degree) * 60.0));
	double sec = (fabs((dec - degree) * 60.0) - min) * 60.0;
	if (sec >= 60.) {
		sec -= 60.;
		min += 1;
	}
	if (min >= 60) {
		min -= 60;
		degree += 1;
	}
	if (degree >= 360) degree = 0;

	gchar *ptr = g_strrstr(format, "lf");
	if (ptr) { // floating point for second
		return g_strdup_printf(format, sig, degree, min, sec);
	}
	return g_strdup_printf(format, sig, degree, min, (int) round(sec));
}

gchar* siril_world_cs_alpha_format(SirilWorldCS *world_cs, const gchar *format) {
	g_return_val_if_fail(world_cs != NULL, NULL);
	g_return_val_if_fail(format != NULL, NULL);
	g_return_val_if_fail(g_utf8_validate (format, -1, NULL), NULL);

	gdouble ra = world_cs->alpha;

	ra = fabs(ra);

	int hour = (int)(ra / 15.0);
	int min = (int)(((ra / 15.0) - hour) * 60.0);
	double sec = ((((ra / 15.0) - hour) * 60.0) - min) * 60.0;
	if (sec >= 60.) {
		sec -= 60.;
		min += 1;
	}
	if (min >= 60) {
		min -= 60;
		hour += 1;
	}
	if (hour >= 24) hour = 0;

	gchar *ptr = g_strrstr(format, "lf");
	if (ptr) { // floating point for second
		return g_strdup_printf(format, hour, min, sec);
	}
	return g_strdup_printf(format, hour, min, (int) round(sec));
}

void siril_world_cs_get_ra_hour_min_sec(SirilWorldCS *world_cs, int *hour, int *min, double *sec) {
	int h, m;
	gdouble s;

	h = (int)(world_cs->alpha / 15.0);
	m = (int)(((world_cs->alpha / 15.0) - h) * 60.0);
	s = ((((world_cs->alpha / 15.0) - h) * 60.0) - m) * 60.0;

	if (hour)
		*hour = h;
	if (min)
		*min = m;
	if (sec)
		*sec = s;
}

void siril_world_cs_get_dec_deg_min_sec(SirilWorldCS *world_cs, int *deg, int *min, double *sec) {
	int d, m;
	gdouble s;

	d = (int) world_cs->delta;
	m = abs((int) ((world_cs->delta - d) * 60.0));
	s = (fabs((world_cs->delta - d) * 60.0) - m) * 60.0;

	if (deg)
		*deg = d;
	if (min)
		*min = m;
	if (sec)
		*sec = s;
}
