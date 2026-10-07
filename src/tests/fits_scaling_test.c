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

/* Reading of FITS images whose integer data is a scaled representation of
 * floating point physical values (BSCALE/BZERO), as opposed to the usual
 * unsigned integer offset trick. */

#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include "core/siril.h"
#include "io/image_format_fits.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

#define W 16
#define H 8
#define NBPIX (W * H)

/* physical = BZERO + BSCALE * stored, giving values in [0.7, 1.3] */
#define BSCALE 1.5e-5
#define BZERO 1.0

static void fill_raw(short *raw) {
	for (int i = 0; i < NBPIX; i++)
		raw[i] = (short)(-20000 + i * (40000 / (NBPIX - 1)));
}

/* extra keywords added to the file to check which ones the rescaling drops.
 * str is written instead of value when it is not NULL */
struct extra_key { const char *key; double value; const char *str; };

/* writes the raw shorts, then adds the scaling keywords so that cfitsio does
 * not apply the inverse transform to the data we give it */
static gchar *write_scaled_fits_full(const char *name, const short *raw,
		gboolean with_scaling, double bscale, double bzero,
		const struct extra_key *extra, int nb_extra) {
	fitsfile *fptr = NULL;
	int status = 0;
	long naxes[2] = { W, H };
	gchar *path = g_build_filename(g_get_tmp_dir(), name, NULL);

	g_unlink(path);
	fits_create_diskfile(&fptr, path, &status);
	fits_create_img(fptr, SHORT_IMG, 2, naxes, &status);
	fits_write_img(fptr, TSHORT, 1, NBPIX, (void *)raw, &status);
	if (with_scaling) {
		fits_update_key(fptr, TDOUBLE, "BSCALE", &bscale, NULL, &status);
		fits_update_key(fptr, TDOUBLE, "BZERO", &bzero, NULL, &status);
	}
	for (int i = 0; i < nb_extra; i++) {
		if (extra[i].str) {
			fits_update_key(fptr, TSTRING, extra[i].key, (void *)extra[i].str, NULL, &status);
		} else {
			double v = extra[i].value;
			fits_update_key(fptr, TDOUBLE, extra[i].key, &v, NULL, &status);
		}
	}
	fits_close_file(fptr, &status);
	cr_assert(status == 0, "could not write %s (cfitsio status %d)", path, status);
	return path;
}

static gchar *write_scaled_fits(const char *name, const short *raw,
		gboolean with_scaling, double bscale, double bzero) {
	return write_scaled_fits_full(name, raw, with_scaling, bscale, bzero, NULL, 0);
}

/* TRUE if the keyword is still among the cards Siril kept, whatever its type */
static gboolean unknown_key_present(const fits *fit, const char *key) {
	gchar **cards = g_strsplit(fit->unknown_keys ? fit->unknown_keys : "", "\n", -1);
	gboolean found = FALSE;

	for (int i = 0; cards[i] && !found; i++) {
		char keyname[FLEN_KEYWORD];
		int length = 0, status = 0;
		if (!cards[i][0])
			continue;
		fits_get_keyname(cards[i], keyname, &length, &status);
		found = !status && !g_strcmp0(keyname, key);
	}
	g_strfreev(cards);
	return found;
}

/* looks a keyword up in the raw cards Siril kept for keywords it doesn't know */
static gboolean unknown_key_value(const fits *fit, const char *key, double *value) {
	gchar **cards = g_strsplit(fit->unknown_keys ? fit->unknown_keys : "", "\n", -1);
	gboolean found = FALSE;

	for (int i = 0; cards[i] && !found; i++) {
		char keyname[FLEN_KEYWORD], val[FLEN_VALUE], comment[FLEN_COMMENT];
		int length = 0, status = 0;
		if (!cards[i][0])
			continue;
		fits_get_keyname(cards[i], keyname, &length, &status);
		fits_parse_value(cards[i], val, comment, &status);
		if (!status && !g_strcmp0(keyname, key)) {
			*value = g_ascii_strtod(val, NULL);
			found = TRUE;
		}
	}
	g_strfreev(cards);
	return found;
}

/* the FITS convention has the first row at the bottom, Siril stores it the
 * same way for its float data, so pixel i of the file is pixel i of fdata */
Test(fits_scaling, physical_values_read_as_float) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	gchar *path = write_scaled_fits("siril_test_phys.fit", raw, TRUE, BSCALE, BZERO);

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	/* BITPIX is 16 in the file but the data is floating point */
	cr_expect(fit.type == DATA_FLOAT);
	cr_expect(fit.bitpix == FLOAT_IMG);
	cr_expect(fit.orig_bitpix == FLOAT_IMG);
	cr_assert(fit.fdata);

	/* physical values are in [0.7, 1.3]: they are all positive, so they are
	 * only divided by the maximum to fit our [0, 1] range */
	double maxi = BZERO + BSCALE * raw[NBPIX - 1];
	for (int i = 0; i < NBPIX; i++) {
		double expected = (BZERO + BSCALE * raw[i]) / maxi;
		cr_expect_float_eq(fit.fdata[i], (float)expected, 1e-5,
				"pixel %d: got %f, expected %f", i, fit.fdata[i], expected);
	}

	/* the transform has to be recorded, it is the only trace left of the
	 * original physical range once the image is saved */
	gboolean recorded = FALSE;
	for (GSList *l = fit.history; l; l = l->next)
		if (strstr((char *)l->data, "physical = pixel /"))
			recorded = TRUE;
	cr_expect(recorded, "the rescaling was not recorded in the history");
	clearfits(&fit);
}

/* a scaling that keeps the data inside [0, 1] must be applied but not rescaled */
Test(fits_scaling, physical_values_already_in_range) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	/* physical = 0.5 + 1e-5 * stored, in [0.3, 0.7] */
	gchar *path = write_scaled_fits("siril_test_phys_inrange.fit", raw, TRUE, 1e-5, 0.5);

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	cr_assert(fit.type == DATA_FLOAT);
	for (int i = 0; i < NBPIX; i++) {
		double expected = 0.5 + 1e-5 * raw[i];
		cr_expect_float_eq(fit.fdata[i], (float)expected, 1e-6,
				"pixel %d: got %f, expected %f", i, fit.fdata[i], expected);
	}

	/* nothing was rescaled here, so there is nothing to record */
	for (GSList *l = fit.history; l; l = l->next)
		cr_expect(!strstr((char *)l->data, "physical = pixel /"),
				"unexpected rescaling record: %s", (char *)l->data);
	clearfits(&fit);
}

/* the unsigned 16-bit trick (BSCALE 1, BZERO 2^15) must keep working */
Test(fits_scaling, unsigned_offset_still_read_as_ushort) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	gchar *path = write_scaled_fits("siril_test_ushort.fit", raw, TRUE, 1.0, 32768.0);

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	cr_expect(fit.type == DATA_USHORT);
	cr_assert(fit.data);
	for (int i = 0; i < NBPIX; i++)
		cr_expect(fit.data[i] == (WORD)(raw[i] + 32768),
				"pixel %d: got %u, expected %d", i, fit.data[i], raw[i] + 32768);
	clearfits(&fit);
}

/* the statistics the instrument put in the file describe the original physical
 * data: once the pixels have been rescaled they no longer apply, and Siril
 * cannot reinterpret them, so they are dropped rather than left to be believed */
Test(fits_scaling, data_keywords_stripped_after_rescaling) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	const struct extra_key extra[] = {
		{ "DATAMIN",  0.7, NULL },
		{ "DATAMED",  1.0, NULL },
		{ "DATAMEAN", 0.9, NULL },
		{ "DATARMS",  0.2, NULL },
		{ "DATASKEW", -1.5, NULL },
		{ "DATALVL",  4095, NULL },	// integer, still a pixel value
		{ "DATATYPE", 0, "science" },	// string, says nothing about the values
		{ "DATANPIX", 1234, NULL },	// SOLARNET: a count of pixels
		{ "DATAVALS", 5678, NULL },	// SDO/JSOC: a count of pixels
		{ "DATASIGN", -1, NULL },	// SDO/JSOC: a sign, not a value
		{ "MFGSMEAN", 0.8, NULL }	// not a statistic of the pixel values
	};
	gchar *path = write_scaled_fits_full("siril_test_stats.fit", raw, TRUE,
			BSCALE, BZERO, extra, G_N_ELEMENTS(extra));

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	double value;
	cr_expect(!unknown_key_value(&fit, "DATAMIN", &value), "DATAMIN was kept");
	cr_expect(!unknown_key_value(&fit, "DATAMED", &value), "DATAMED was kept");
	cr_expect(!unknown_key_value(&fit, "DATAMEAN", &value), "DATAMEAN was kept");
	cr_expect(!unknown_key_value(&fit, "DATARMS", &value), "DATARMS was kept");
	cr_expect(!unknown_key_value(&fit, "DATASKEW", &value), "DATASKEW was kept");
	/* integers are pixel values too */
	cr_expect(!unknown_key_present(&fit, "DATALVL"), "integer DATALVL was kept");

	/* a DATA* keyword that is not a number carries no pixel value, so it stays */
	cr_expect(unknown_key_present(&fit, "DATATYPE"), "string DATATYPE was removed");

	/* neither do counts and flags, which the rescaling leaves valid */
	cr_expect(unknown_key_value(&fit, "DATANPIX", &value), "DATANPIX was removed");
	cr_expect_float_eq(value, 1234.0, 1e-9);
	cr_expect(unknown_key_value(&fit, "DATAVALS", &value), "DATAVALS was removed");
	cr_expect_float_eq(value, 5678.0, 1e-9);
	cr_expect(unknown_key_value(&fit, "DATASIGN", &value), "DATASIGN was removed");
	cr_expect_float_eq(value, -1.0, 1e-9);

	/* and nothing outside the DATA* family is touched */
	cr_assert(unknown_key_value(&fit, "MFGSMEAN", &value), "MFGSMEAN was removed");
	cr_expect_float_eq(value, 0.8, 1e-9);

	/* and the removal is recorded, like the rescaling itself */
	gboolean recorded = FALSE;
	for (GSList *l = fit.history; l; l = l->next)
		if (strstr((char *)l->data, "Removed DATA*"))
			recorded = TRUE;
	cr_expect(recorded, "the removal was not recorded in the history");
	clearfits(&fit);
}

/* without a rescaling the pixels are the physical values, so the statistics
 * still describe them and must be left alone */
Test(fits_scaling, data_keywords_kept_when_not_rescaled) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	const struct extra_key extra[] = { { "DATAMED", 0.5, NULL } };
	/* physical = 0.5 + 1e-5 * stored, in [0.3, 0.7], already in our range */
	gchar *path = write_scaled_fits_full("siril_test_stats_inrange.fit", raw, TRUE,
			1e-5, 0.5, extra, G_N_ELEMENTS(extra));

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	double value;
	cr_assert(unknown_key_value(&fit, "DATAMED", &value), "DATAMED was removed");
	cr_expect_float_eq(value, 0.5, 1e-9);
	clearfits(&fit);
}

/* DATAMAX is deliberately never written: read_fits_with_convert() would divide
 * float data above 10 by USHRT_MAX on the next load */
Test(fits_scaling, datamax_not_written_on_save) {
	short raw[NBPIX];
	fits fit = { 0 };
	fitsfile *fptr = NULL;
	double datamax = 0.0;
	int status = 0;
	fill_raw(raw);
	gchar *path = write_scaled_fits("siril_test_save_src.fit", raw, FALSE, 1.0, 0.0);

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	com.pref.ext = g_strdup(".fit");	// used by savefits() to name the file
	gchar *out = g_build_filename(g_get_tmp_dir(), "siril_test_save_out", NULL);
	cr_assert(savefits(out, &fit) == 0);
	clearfits(&fit);

	gchar *outfit = g_strdup_printf("%s.fit", out);
	g_free(out);
	fits_open_diskfile(&fptr, outfit, READONLY, &status);
	cr_assert(status == 0, "could not reopen %s", outfit);
	fits_read_key(fptr, TDOUBLE, "DATAMAX", &datamax, NULL, &status);
	cr_expect(status == KEY_NO_EXIST, "DATAMAX must not be written (got %f)", datamax);
	status = 0;
	fits_close_file(fptr, &status);
	g_free(outfit);
}

/* and so must plain signed 16-bit data without any scaling keyword */
Test(fits_scaling, plain_short_still_read_as_ushort) {
	short raw[NBPIX];
	fits fit = { 0 };
	fill_raw(raw);
	gchar *path = write_scaled_fits("siril_test_short.fit", raw, FALSE, 1.0, 0.0);

	cr_assert(readfits(path, &fit, NULL, FALSE) == 0);
	g_free(path);

	cr_expect(fit.type == DATA_USHORT);
	cr_assert(fit.data);
	for (int i = 0; i < NBPIX; i++)
		cr_expect(fit.data[i] == (WORD)(raw[i] + 32768),
				"pixel %d: got %u, expected %d", i, fit.data[i], raw[i] + 32768);
	clearfits(&fit);
}
