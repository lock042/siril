/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#ifndef SRC_CORE_SIRIL_DATE_H_
#define SRC_CORE_SIRIL_DATE_H_

#include <glib.h>
#include <stdint.h>

gchar *build_timestamp_filename();
GDateTime *ser_timestamp_to_date_time(guint64 timestamp);
guint64 date_time_to_ser_timestamp(GDateTime *dt);
double date_time_to_Julian(GDateTime *dt);
GDateTime* Julian_to_date_time(gdouble jd);
GDateTime *FITS_date_to_date_time(gchar *date);
gchar *date_time_to_FITS_date(GDateTime *date);
gchar *date_time_to_date(GDateTime *datetime);
gchar *date_time_to_date_time(GDateTime *datetime);
GDateTime *julian_sec_to_date(uint32_t jsecs, uint32_t ms);
double timediff_in_s(GDateTime *dt1, GDateTime *dt2);

#endif /* SRC_CORE_SIRIL_DATE_H_ */
