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

#include <criterion/criterion.h>
#include <math.h>
#include "core/siril_world_cs.h"


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

int test_parsing() {
	SirilWorldCS *ecs = siril_world_cs_new_from_objct_ra_dec("02:13:42.6", "-4:0:36");
	CHECK(ecs, "did not parse with colon separator\n");
	CHECK(fabs(siril_world_cs_get_alpha(ecs) - 33.4275) < 0.00001, "bad alpha parse\n");
	CHECK(fabs(siril_world_cs_get_delta(ecs) - -4.01) < 0.00001, "bad delta parse\n");
	siril_world_cs_unref(ecs);

	ecs = siril_world_cs_new_from_objct_ra_dec("02 13 42.6", "-4 0 36");
	CHECK(ecs, "did not parse with space separator\n");
	CHECK(fabs(siril_world_cs_get_alpha(ecs) - 33.4275) < 0.00001, "bad alpha parse\n");
	CHECK(fabs(siril_world_cs_get_delta(ecs) - -4.01) < 0.00001, "bad delta parse\n");
	siril_world_cs_unref(ecs);

	ecs = siril_world_cs_new_from_objct_ra_dec("33.4275", "-4.01");
	CHECK(ecs, "did not parse decimal\n");
	CHECK(fabs(siril_world_cs_get_alpha(ecs) - 33.4275) < 0.00001, "bad alpha parse\n");
	CHECK(fabs(siril_world_cs_get_delta(ecs) - -4.01) < 0.00001, "bad delta parse\n");
	siril_world_cs_unref(ecs);

	ecs = siril_world_cs_new_from_objct_ra_dec("0", "4:00:36");
	CHECK(ecs, "did not parse decimal\n");
	CHECK(fabs(siril_world_cs_get_alpha(ecs)) < 0.00000001, "bad alpha parse\n");
	CHECK(fabs(siril_world_cs_get_delta(ecs) - 4.01) < 0.00001, "bad delta parse\n");
	siril_world_cs_unref(ecs);


	/** test for RA/DEC conversion */
	ecs = siril_world_cs_new_from_a_d(274.2499, 42.9601);

	gchar *ra = siril_world_cs_alpha_format(ecs, "%02d %02d %.3lf");
	CHECK(!g_strcmp0(ra, "18 16 59.976"), "RA conversion is not good: %s", ra);
	g_free(ra);

	ra = siril_world_cs_alpha_format(ecs, "%02d %02d %02d");
	CHECK(!g_strcmp0(ra, "18 17 00"), "RA conversion is not good: %s", ra);
	g_free(ra);

	gchar *dec = siril_world_cs_delta_format(ecs, "%c%02d %02d %.3lf");
	CHECK(!g_strcmp0(dec, "+42 57 36.360"), "DEC conversion is not good: %s", dec);
	g_free(dec);

	dec = siril_world_cs_delta_format(ecs, "%c%02d %02d %02d");
	CHECK(!g_strcmp0(dec, "+42 57 36"), "DEC conversion is not good: %s", dec);
	g_free(dec);

	siril_world_cs_unref(ecs);


	return 0;
}

#ifdef WITH_MAIN
int main() {
	int retval = test_parsing();
	if (retval)
		fprintf(stderr, "TESTS FAILED\n");
	else fprintf(stderr, "ALL TESTS PASSED\n");
	return retval;
}
#else // with criterion
Test(check_date, test1) { cr_assert(!test_parsing()); }
#endif
