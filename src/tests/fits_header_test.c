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

/* Reading of the FITS header keywords into fit->keywords */

#include <string.h>

#include <criterion/criterion.h>
#include <fitsio.h>

#include "core/siril.h"
#include "io/image_format_fits.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

#define W 4
#define H 4

/* writes a small 16-bit image with the given raw header cards and reads it back */
static void read_header(const char *name, const char **cards, int nb_cards, fits *fit) {
	fitsfile *fptr = NULL;
	int status = 0;
	long naxes[2] = { W, H };
	WORD pixels[W * H] = { 0 };
	gchar *path = g_build_filename(g_get_tmp_dir(), name, NULL);

	g_unlink(path);
	fits_create_diskfile(&fptr, path, &status);
	fits_create_img(fptr, USHORT_IMG, 2, naxes, &status);
	fits_write_img(fptr, TUSHORT, 1, W * H, pixels, &status);
	for (int i = 0; i < nb_cards; i++)
		fits_write_record(fptr, cards[i], &status);
	fits_close_file(fptr, &status);
	cr_assert(status == 0, "could not write %s (cfitsio status %d)", path, status);

	cr_assert(readfits(path, fit, NULL, FALSE) == 0);
	g_unlink(path);
	g_free(path);
}

static gboolean unknown_key_present(const fits *fit, const char *key) {
	gchar **cards = g_strsplit(fit->unknown_keys ? fit->unknown_keys : "", "\n", -1);
	gboolean found = FALSE;
	for (int i = 0; cards[i] && !found; i++)
		found = g_str_has_prefix(cards[i], key) && cards[i][strlen(key)] == ' ';
	g_strfreev(cards);
	return found;
}

Test(fits_header, primary_keywords) {
	const char *cards[] = {
		"EXPTIME =                 30.5 / [s] Exposure time",
		"GAIN    =                  120",
		"CCD-TEMP=                -10.2",
		"OBJECT  = '  M 42   '",
		"BAYERPAT= 'RGGB    '",
		"FOCALLEN=                  530",
		"XPIXSZ  =                 3.76",
		"XBINNING=                    2"
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_primary.fit", cards, G_N_ELEMENTS(cards), &fit);

	cr_expect_float_eq(fit.keywords.exposure, 30.5, 1e-9);
	cr_expect(fit.keywords.key_gain == 120);
	cr_expect_float_eq(fit.keywords.ccd_temp, -10.2, 1e-9);
	cr_expect_str_eq(fit.keywords.object, "M 42");	// unquoted and stripped
	cr_expect_str_eq(fit.keywords.bayer_pattern, "RGGB");
	cr_expect_float_eq(fit.keywords.focal_length, 530.0, 1e-9);
	cr_expect(fit.focalkey);
	cr_expect_float_eq(fit.keywords.pixel_size_x, 3.76, 1e-9);
	cr_expect(fit.pixelkey);
	cr_expect(fit.keywords.binning_x == 2);
	clearfits(&fit);
}

/* the alternative names other software use for the same data */
Test(fits_header, secondary_keywords) {
	const char *cards[] = {
		"EXPOSURE=                  120",
		"CCDTEMP =                 -5.0",
		"NCOMBINE=                   42",
		"BLKLEVEL=                   30",
		"FRAMETYP= 'Light   '",
		"FLENGTH =                 0.53 / [m]"
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_secondary.fit", cards, G_N_ELEMENTS(cards), &fit);

	cr_expect_float_eq(fit.keywords.exposure, 120.0, 1e-9);
	cr_expect_float_eq(fit.keywords.ccd_temp, -5.0, 1e-9);
	cr_expect(fit.keywords.stackcnt == 42);
	cr_expect(fit.keywords.key_offset == 30);
	cr_expect_str_eq(fit.keywords.image_type, "Light");
	cr_expect_float_eq(fit.keywords.focal_length, 530.0, 1e-9);	// m to mm
	clearfits(&fit);
}

/* integers in scientific notation, and values the handlers fix or reject */
Test(fits_header, value_sanitizing) {
	const char *cards[] = {
		"GAIN    =            5.600E+01",
		"XBINNING=                    0",
		"BAYERPAT= 'NONE    '",
		"STACKCNT=                   -3"	// out of range, falls back to the default
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_sanitize.fit", cards, G_N_ELEMENTS(cards), &fit);

	cr_expect(fit.keywords.key_gain == 56);
	cr_expect(fit.keywords.binning_x == 1);
	cr_expect_str_eq(fit.keywords.bayer_pattern, "");
	cr_expect(fit.keywords.stackcnt == 1);
	clearfits(&fit);
}

/* DATE-OBS without a time is completed by TIME-OBS */
Test(fits_header, date_obs_with_time_obs) {
	const char *cards[] = {
		"DATE-OBS= '2024-03-15'",
		"TIME-OBS= '21:30:15.250'"
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_date.fit", cards, G_N_ELEMENTS(cards), &fit);

	GDateTime *d = fit.keywords.date_obs;
	cr_assert(d, "DATE-OBS was not read");
	cr_expect(g_date_time_get_year(d) == 2024);
	cr_expect(g_date_time_get_month(d) == 3);
	cr_expect(g_date_time_get_day_of_month(d) == 15);
	cr_expect(g_date_time_get_hour(d) == 21);
	cr_expect(g_date_time_get_minute(d) == 30);
	cr_expect_float_eq(g_date_time_get_seconds(d), 15.25, 1e-6);
	clearfits(&fit);
}

/* SITELAT, SITELONG and RA come as numbers or as strings depending on the software */
Test(fits_header, coordinates_as_string_or_number) {
	const char *cards[] = {
		"SITELAT = '+48:51:24'",
		"SITELONG=               2.3522",
		"RA      = '83.6331 '"
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_coords.fit", cards, G_N_ELEMENTS(cards), &fit);

	cr_expect_float_eq(fit.keywords.sitelat, 48.0 + 51.0 / 60.0 + 24.0 / 3600.0, 1e-6);
	cr_expect_float_eq(fit.keywords.sitelong, 2.3522, 1e-9);
	cr_expect_float_eq(fit.keywords.wcsdata.ra, 83.6331, 1e-6);
	/* and the string form is rewritten as HMS */
	cr_expect(g_str_has_prefix(fit.keywords.wcsdata.objctra, "05 34 3"),
			"OBJCTRA is '%s'", fit.keywords.wcsdata.objctra);
	clearfits(&fit);
}

/* keywords Siril doesn't know are kept verbatim to be saved again, but not
 * HISTORY or checksums, which are handled or recomputed elsewhere */
Test(fits_header, unknown_keywords_kept) {
	const char *cards[] = {
		"MYKEY   = 'hello   '",
		"SWOWNER = 'someone '",
		"HISTORY processed somewhere",
		"CHECKSUM= 'abcdefghijklmnop'",
		"DATASUM = '0       '"
	};
	fits fit = { 0 };
	read_header("siril_test_hdr_unknown.fit", cards, G_N_ELEMENTS(cards), &fit);

	cr_expect(unknown_key_present(&fit, "MYKEY"));
	cr_expect(unknown_key_present(&fit, "SWOWNER"));
	cr_expect(!unknown_key_present(&fit, "HISTORY"));
	cr_expect(!unknown_key_present(&fit, "CHECKSUM"));
	cr_expect(!unknown_key_present(&fit, "DATASUM"));
	/* and the structural keywords never end up there */
	cr_expect(!unknown_key_present(&fit, "NAXIS1"));
	cr_expect(!unknown_key_present(&fit, "BITPIX"));
	clearfits(&fit);
}
