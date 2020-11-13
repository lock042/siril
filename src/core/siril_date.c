/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include <stdio.h>

#include "siril_date.h"

/**
 * From a datetime it computes the Julian date needed in photometry
 * (code borrowed from muniwin)
 * @param dt timestamp in datetime format
 * @return the Julian date
 */
double encode_to_Julian_date(GDateTime *dt) {
	double jd1;
	int before, d1, d2;
	int year, month, day;
	int hour, min, sec, ms;

	if (!dt) return 0;

	g_date_time_get_ymd(dt, &year, &month, &day);

	/* Check date and time */
	if (day <= 0 || year <= 0 || month <= 0)
		return 0;

	/* Compute Julian date from input citizen year, month and day. */
	/* Tested for YEAR>0 except 1582-10-07/15 */
	if (year > 1582) {
		before = 0;
	} else if (year < 1582) {
		before = 1;
	} else if (month > 10) {
		before = 0;
	} else if (month < 10) {
		before = 1;
	} else if (day >= 15) {
		before = 0;
	} else {
		before = 1;
	}
	if (month <= 2) {
		d1 = (int) (365.25 * (year - 1));
		d2 = (int) (30.6001 * (month + 13));
	} else {
		d1 = (int) (365.25 * (year));
		d2 = (int) (30.6001 * (month + 1));
	}

	hour = g_date_time_get_hour(dt);
	min = g_date_time_get_minute(dt);
	sec = g_date_time_get_second(dt);
	ms = g_date_time_get_microsecond(dt);

	jd1 = 1720994.5 + d1 + d2 + day;
	jd1 += 1.0 * hour / 24;
	jd1 += 1.0 * min / 1440.0;
	jd1 += 1.0 * sec / 86400.0;
	jd1 += 1.0 * ms / 86400000.0;

	if (before) {
		return jd1;
	} else {
		return jd1 + 2 - (year / 100) + (year / 400);
	}
}

/**
 * From a char * in FITS format to a GDateTime
 * @param date
 * @return a GDateTime
 */
GDateTime *siril_FITS_to_date_time(char *date) {
	int year = 0, month = 0, day = 0, hour = 0, min = 0;
	float sec = 0.f;

	if (date[0] == '\0')
		return 0;

	if (sscanf(date, "%04d-%02d-%02dT%02d:%02d:%f", &year, &month, &day, &hour,
			&min, &sec) != 6) {
		return 0;
	}
	GTimeZone *tz = g_time_zone_new_utc();
	GDateTime *new_date = g_date_time_new(tz, year, month, day, hour, min, sec);
	g_time_zone_unref(tz);

	return new_date;
}

/**
 *
 * @param
 * @return
 */
gchar *siril_format_date_time(GDateTime *date) {
	return g_date_time_format(date, "%Y-%m-%dT%H:%M:%S.%f");
}

/**
 *
 * @param from
 * @return
 */
GDateTime *siril_copy_date_time(GDateTime *from) {
	GDateTime *to = NULL;
	if (from) {
		GTimeZone *tz = g_date_time_get_timezone(from);

		to = g_date_time_new(tz, g_date_time_get_year(from),
				g_date_time_get_month(from),
				g_date_time_get_day_of_month(from),
				g_date_time_get_hour(from),
				g_date_time_get_minute(from),
				g_date_time_get_seconds(from));

		g_time_zone_unref(tz);
	}
	return to;
}
