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

#include <stdio.h>
#include <criterion/criterion.h>
#include "core/siril_date.h"

#define UNDER_US      G_GUINT64_CONSTANT(7)
#define INPUT_TIME    G_GUINT64_CONSTANT(637232717926133380) + UNDER_US
#define SER_TIME_1970 G_GUINT64_CONSTANT(621355968000000000) // 621.355.968.000.000.000 ticks between 1st Jan 0001 and 1st Jan 1970.

#ifdef WITH_MAIN
#include <stdio.h>
#define CHECK(cond, ...) \
	if (!(cond)) { \
		fprintf(stderr, __VA_ARGS__); \
		return 1; \
	}
#else
#define CHECK cr_expect
#endif


/**
 *  Test consistency of siril_date_time functions
 *   */
int test_ser_date() {
	GDateTime *dt1, *dt2, *dt3, *dt4;
	guint64 diff, ts;
	gchar *date_str;

	dt1 = ser_timestamp_to_date_time(INPUT_TIME);
	ts = date_time_to_ser_timestamp(dt1);

	/**
	 *  ser timestamp precision is down 0.1 microsecond while our
	 *  structure is accurate down to 1 microsecond
	 */
	diff = INPUT_TIME - ts;
	CHECK(diff == UNDER_US, "Failed with retval=%"G_GUINT64_FORMAT, diff);

	dt2 = g_date_time_new_from_iso8601 ("2016-11-30T22:10:42Z", NULL);
	ts = date_time_to_ser_timestamp(dt2);
	dt3 = ser_timestamp_to_date_time(ts);
	CHECK(g_date_time_equal(dt2, dt3), "date_time from ser are not equal");

	/**
	 *  Test FITS date time consistency
	 */
	date_str = date_time_to_FITS_date(dt2);
	dt4 = FITS_date_to_date_time(date_str);
	CHECK(g_date_time_equal(dt2, dt4), "date_time from FITS are not equal");

	g_date_time_unref(dt1);
	g_date_time_unref(dt2);
	g_date_time_unref(dt3);
	g_date_time_unref(dt4);
	return 0;
}

int test_julian_date() {
	GDateTime *dt = g_date_time_new_from_iso8601("2023-05-20T21:15:55Z", NULL);
	//printf("expected: 2460085.386053, obtained: %f\n", date_time_to_Julian(dt));
	CHECK(abs(date_time_to_Julian(dt) - 2460085.386053) < 1e-7, "Julian date conversion failed\n");
	g_date_time_unref(dt);

	// same, with centiseconds (at the limit of double precision)
	dt = g_date_time_new_from_iso8601("2023-05-20T21:15:55.2Z", NULL);
	printf("expected: 2460085.386056, obtained: %.7f\n", date_time_to_Julian(dt));
	CHECK(abs(date_time_to_Julian(dt) - 2460085.386056) < 1e-7, "Julian date conversion failed\n");
	g_date_time_unref(dt);
	return 0;
}

int test_julian_seconds() {
	GDateTime *ref = g_date_time_new_utc(1995, 10, 10, 0, 0, 0.0);
	GDateTime *test = julian_sec_to_date(0, 0);
	CHECK(g_date_time_equal(ref, test), "Julian seconds base is wrong");
	g_date_time_unref(test);
	g_date_time_unref(ref);

	// 2023-05-31T17:01:56
	ref = g_date_time_new_utc(2023, 5, 31, 17, 1, 56.42);
	test = julian_sec_to_date(872269316, 420000);
	CHECK(g_date_time_equal(ref, test), "Julian seconds conversion failed");
	printf("converted truncated Julian seconds: %s\n", g_date_time_format_iso8601(ref));
	g_date_time_unref(test);
	g_date_time_unref(ref);
	return 0;
}

#ifdef WITH_MAIN
int main() {
	int retval = test_ser_date();
	retval += test_julian_date();
	retval += test_julian_seconds();
	if (retval)
		fprintf(stderr, "TESTS FAILED\n");
	else fprintf(stderr, "ALL TESTS PASSED\n");
	return retval;
}
#else // with criterion
Test(check_date, test_ser) { cr_assert(!test_ser_date()); }
Test(check_date, test_julian_date) { cr_assert(!test_julian_date()); }
Test(check_date, test_julian_secs) { cr_assert(!test_julian_seconds()); }
#endif
