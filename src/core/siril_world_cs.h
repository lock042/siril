/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_CORE_SIRIL_WORLD_CS_H_
#define SRC_CORE_SIRIL_WORLD_CS_H_

#include <glib.h>

/* This is an implementation of the equatorial coordinate system, consisting of
 * a couple or RA and DEC coordinates, with their parsing and transformation
 * tools.
 * Internal data is stored in degrees.
 */

typedef struct _SirilWorldCS SirilWorldCS;

SirilWorldCS *siril_world_cs_ref(SirilWorldCS *world_cs);
SirilWorldCS *siril_world_cs_copy(SirilWorldCS *world_cs);
void siril_world_cs_unref(SirilWorldCS *world_cs);

SirilWorldCS* siril_world_cs_new_from_a_d(gdouble alpha, gdouble delta);
SirilWorldCS* siril_world_cs_new_from_ra_dec(gdouble ra_h, gdouble ra_m, gdouble ra_s, gdouble dec_deg, gdouble dec_m, gdouble dec_s);
SirilWorldCS* siril_world_cs_new_from_objct_ra_dec(gchar *objctra, gchar *objctdec);


gdouble siril_world_cs_get_alpha(SirilWorldCS *world_cs);
gdouble siril_world_cs_get_delta(SirilWorldCS *world_cs);
gchar *siril_world_cs_delta_format_from_double(gdouble dec, const gchar *format);
gchar* siril_world_cs_delta_format(SirilWorldCS *world_cs, const gchar *format);
gchar *siril_world_cs_alpha_format_from_double(gdouble ra, const gchar *format);
gchar* siril_world_cs_alpha_format(SirilWorldCS *world_cs, const gchar *format);
void siril_world_cs_get_ra_hour_min_sec(SirilWorldCS *world_cs, int *hour, int *min, double *sec);
void siril_world_cs_get_dec_deg_min_sec(SirilWorldCS *world_cs, int *deg, int *min, double *sec);

double parse_hms(const char *objctra);
double parse_dms(const char *objctdec);

#endif /* SRC_CORE_SIRIL_WORLD_CS_H_ */
