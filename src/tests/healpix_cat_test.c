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

/* Checks that a file which is not a Siril Healpix catalogue is rejected
 * instead of being handed to the healpix library: an out of range indexing
 * level made it throw, and the exception, crossing the extern "C" boundary,
 * terminated Siril while plate solving. */

#include <criterion/criterion.h>
#include <stdio.h>
#include <string.h>

#include "core/siril.h"
#include "io/siril_catalogues.h"
#include "io/local_catalogues.h"
#include "io/healpix/healpix_cat.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

/* writes a file of the given size, filled with a byte pattern that is not a
 * valid catalogue header, and points the preferences at it */
static gchar *make_bogus_catalogue(size_t size) {
	gchar *path = g_build_filename(g_get_tmp_dir(), "siril_bogus_catalogue.dat", NULL);
	FILE *f = g_fopen(path, "wb");
	cr_assert_not_null(f, "could not create the test file");
	for (size_t i = 0; i < size; i++)
		fputc(0xA5, f);	// gives an indexing level of 165, way out of range
	fclose(f);
	g_free(com.pref.catalogue_paths[4]);
	com.pref.catalogue_paths[4] = path;
	return path;
}

/* writes a file holding a valid, if empty, monolithic level 8 catalogue header
 * and points the preferences at it */
static void make_valid_catalogue_header() {
	guint8 header[128] = { 0 };
	strcpy((char *)header, "Siril Gaia DR3 astrometric extract");
	header[48] = 3;	// gaia_version: DR3
	header[49] = 8;	// healpix_level
	header[50] = 1;	// cat_type: astrometric
	header[51] = 0;	// not chunked
	gchar *path = g_build_filename(g_get_tmp_dir(), "siril_valid_catalogue.dat", NULL);
	FILE *f = g_fopen(path, "wb");
	cr_assert_not_null(f, "could not create the test file");
	cr_assert_eq(fwrite(header, 1, sizeof(header), f), sizeof(header));
	fclose(f);
	g_free(com.pref.catalogue_paths[4]);
	com.pref.catalogue_paths[4] = path;
}

Test(healpix_cat, valid_catalogue_is_available) {
	make_valid_catalogue_header();
	cr_expect_neq(local_gaia_astro_available(), 0,
			"a valid catalogue header must be accepted");
	g_unlink(com.pref.catalogue_paths[4]);
}

Test(healpix_cat, invalid_catalogue_is_not_available) {
	make_bogus_catalogue(4096);
	cr_expect_eq(local_gaia_astro_available(), 0,
			"a file that is not a catalogue must not be reported as available");
	cr_expect_not(local_gaia_available(),
			"local Gaia solving must not be selected for an invalid catalogue");
	g_unlink(com.pref.catalogue_paths[4]);
}

Test(healpix_cat, invalid_catalogue_query_fails_cleanly) {
	make_bogus_catalogue(4096);
	deepStarData *stars = (deepStarData*) 0x1;	// must be reset by the callee
	uint32_t nb_stars = 42;
	int retval = get_raw_stars_from_local_gaia_astro_catalogue(302.0, 60.0, 30.0,
			17.0, FALSE, &stars, &nb_stars);
	cr_expect_neq(retval, 0, "querying an invalid catalogue must report an error");
	cr_expect_null(stars);
	cr_expect_eq(nb_stars, 0);
	g_unlink(com.pref.catalogue_paths[4]);
}

Test(healpix_cat, truncated_catalogue_query_fails_cleanly) {
	make_bogus_catalogue(12);	// smaller than a header
	deepStarData *stars = (deepStarData*) 0x1;
	uint32_t nb_stars = 42;
	int retval = get_raw_stars_from_local_gaia_astro_catalogue(302.0, 60.0, 30.0,
			17.0, FALSE, &stars, &nb_stars);
	cr_expect_neq(retval, 0, "querying a truncated catalogue must report an error");
	cr_expect_null(stars);
	cr_expect_eq(nb_stars, 0);
	g_unlink(com.pref.catalogue_paths[4]);
}
