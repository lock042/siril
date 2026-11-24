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

#ifndef WITH_MAIN
#include <criterion/criterion.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "algos/demosaicing.h"
#include "core/siril.h"
#include "core/siril_date.h"
#include "io/image_format_fits.h"
#include "io/ser.h"

#ifdef WITH_MAIN
#define CHECK(cond, ...) \
	if (!(cond)) { \
		fprintf(stderr, __VA_ARGS__); \
		return 1; \
	}
#else
#define CHECK cr_expect

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image
#endif

#ifndef _WIN32
#define TMP_FILE1 "/tmp/test_tmp1.ser"
#define TMP_FILE2 "/tmp/test_tmp2.ser"
#define TMP_FILE3 "/tmp/test_tmp3.ser"
#define TMP_FILE4 "/tmp/test_tmp4.ser"
#define TMP_FILE5 "/tmp/test_tmp5.ser"
#define TMP_FILE6 "/tmp/test_tmp6.ser"
#define TMP_FILE7 "/tmp/test_tmp7.ser"
#define TMP_FILE8 "/tmp/test_tmp8.ser"
#define TMP_FILE9 "/tmp/test_tmp9.ser"
#define TMP_FILE10 "/tmp/test_tmp10.ser"
#define TMP_FILE11 "/tmp/test_tmp11.ser"
#else
#define TMP_FILE1 ".\\test_tmp1.ser"
#define TMP_FILE2 ".\\test_tmp2.ser"
#define TMP_FILE3 ".\\test_tmp3.ser"
#define TMP_FILE4 ".\\test_tmp4.ser"
#define TMP_FILE5 ".\\test_tmp5.ser"
#define TMP_FILE6 ".\\test_tmp6.ser"
#define TMP_FILE7 ".\\test_tmp7.ser"
#define TMP_FILE8 ".\\test_tmp8.ser"
#define TMP_FILE9 ".\\test_tmp9.ser"
#define TMP_FILE10 ".\\test_tmp10.ser"
#define TMP_FILE11 ".\\test_tmp11.ser"
#endif

static fits *create_image(int w, int h, int layers) {
	fits *fit = NULL;
	if (new_fit_image(&fit, w, h, layers, DATA_USHORT))
		return NULL;
	return fit;
}

int test_ser_image_number() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE1, ser, TRUE, NULL), "create file\n");
	fits *fit = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 1), "writing image\n");
	fits *fit3 = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit3, 2), "writing image\n");
	CHECK(!ser_write_and_close(ser), "close file\n");
	CHECK(!ser_open_file(TMP_FILE1, ser), "reopen\n");
	CHECK(ser->frame_count == 3, "wrong number of frames\n");
	free(ser);
	CHECK(!unlink(TMP_FILE1), "error unlinking file " TMP_FILE1);
	fprintf(stdout, "* test 1 passed *\n\n");
	return 0;
}

int test_ser_image_overwrite() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE2, ser, TRUE, NULL), "create file\n");
	fits *fit = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 0), "writing image 2\n");
	CHECK(ser_write_and_close(ser), "close file succeeded with lost image\n");

	CHECK(ser_open_file(TMP_FILE2, ser), "reopen found a file\n");

	free(ser);
	unlink(TMP_FILE2);
	fprintf(stdout, "* test 2 passed *\n\n");
	return 0;
}

int test_ser_image_sizes() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE3, ser, TRUE, NULL), "create file\n");
	fits *fit = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(21, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 1), "writing image\n");
	CHECK(ser_write_and_close(ser), "close file succeeded for error write\n");

	CHECK(ser_open_file(TMP_FILE3, ser), "reopen found a file\n");

	free(ser);
	unlink(TMP_FILE3);
	fprintf(stdout, "* test 3 passed *\n\n");
	return 0;
}

static void set_fits_date(fits *fit, gint64 time) {
	fit->keywords.date_obs = g_date_time_new_from_unix_utc(time);
}

int test_ser_dates() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE4, ser, TRUE, NULL), "create file\n");
	fits *fit = create_image(20, 10, 1);
	set_fits_date(fit, 100);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(20, 10, 1);
	set_fits_date(fit2, 200);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 1), "writing image\n");
	fits *fit3 = create_image(20, 10, 1);
	set_fits_date(fit3, 300);
	CHECK(!ser_write_frame_from_fit(ser, fit3, 2), "writing image\n");
	CHECK(!ser_write_and_close(ser), "close file\n");

	CHECK(!ser_open_file(TMP_FILE4, ser), "reopen\n");
	CHECK(ser->frame_count == 3, "wrong number of frames\n");
	CHECK(ser->ts, "no date information in SER\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[0])) == 100,
			"first image date is wrong\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[1])) == 200,
			"second image date is wrong\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[2])) == 300,
			"third image date is wrong\n");

	free(ser);
	CHECK(!unlink(TMP_FILE4), "error unlinking file " TMP_FILE4);
	fprintf(stdout, "* test 4 passed *\n\n");
	return 0;
}

int test_ser_with_holes() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE5, ser, TRUE, NULL), "create file\n");
	fits *fit = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 1), "writing image 2\n");
	fits *fit3 = create_image(20, 10, 1);
	CHECK(!ser_write_frame_from_fit(ser, fit3, 2), "writing image 3\n");
	CHECK(!ser_write_and_close(ser), "close file\n");

	CHECK(!ser_open_file(TMP_FILE5, ser), "reopen\n");
	CHECK(ser->frame_count == 3, "wrong number of frames\n");

	free(ser);
	CHECK(!unlink(TMP_FILE5), "error unlinking file " TMP_FILE5);
	fprintf(stdout, "* test 5 passed *\n\n");
	return 0;
}

int test_ser_ooo_write() {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE6, ser, TRUE, NULL), "create file\n");
	CHECK(!ser_write_frame_from_fit(ser, NULL, 0), "writing image\n");
	fits *fit = create_image(20, 10, 1);
	set_fits_date(fit, 100);
	CHECK(!ser_write_frame_from_fit(ser, fit, 1), "writing image\n");
	fits *fit2 = create_image(20, 10, 1);
	set_fits_date(fit2, 200);
	fits *fit3 = create_image(20, 10, 1);
	set_fits_date(fit3, 300);
	CHECK(!ser_write_frame_from_fit(ser, fit3, 4), "writing image\n");
	CHECK(!ser_write_frame_from_fit(ser, fit2, 2), "writing image\n");
	CHECK(!ser_write_frame_from_fit(ser, NULL, 3), "writing image\n");
	CHECK(!ser_write_frame_from_fit(ser, NULL, 6), "writing image\n");
	CHECK(!ser_write_frame_from_fit(ser, NULL, 5), "writing image\n");
	CHECK(!ser_write_and_close(ser), "close file\n");

	CHECK(!ser_open_file(TMP_FILE6, ser), "reopen\n");
	CHECK(ser->frame_count == 3, "wrong number of frames\n");
	CHECK(ser->ts, "no date information in SER\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[0])) == 100,
			"first image date is wrong\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[1])) == 200,
			"second image date is wrong\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[2])) == 300,
			"third image date is wrong\n");

	free(ser);
	CHECK(!unlink(TMP_FILE6), "error unlinking file " TMP_FILE6);
	fprintf(stdout, "* test 6 passed *\n\n");
	return 0;
}

static struct ser_struct *create_fake_ser(int w, int h, int with_nb_frames, char *obs, guint64 utc) {
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	ser->file_id = strdup("LUCAM-RECORDER");
	ser->lu_id = 0;
	ser->color_id = SER_MONO;
	ser->little_endian = SER_LITTLE_ENDIAN;
	ser->image_width = w;
	ser->image_height = h;
	ser->bit_pixel_depth = 16;
	ser->frame_count = with_nb_frames;
	strcpy(ser->observer, obs);
	memset(ser->instrument, 0, 40);
	memset(ser->telescope, 0, 40);
	memset(&ser->date, 0, 8);
	ser->date_utc = utc;
	ser->byte_pixel_depth = SER_PIXEL_DEPTH_16;
	ser->number_of_planes = 1;
	ser->ts = malloc(with_nb_frames * sizeof(guint64));
	ser->ts_alloc = with_nb_frames;
	return ser;
}

int test_ser_create_from_copy() {
	char *observer_str = "super observer";
	guint64 utc_time = (guint64)100;
	struct ser_struct *orig = create_fake_ser(20, 10, 3, observer_str, utc_time);
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);

	CHECK(!ser_create_file(TMP_FILE7, ser, TRUE, orig), "create file\n");
	fits *fit = create_image(40, 20, 3);
	set_fits_date(fit, 100);
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	fits *fit2 = create_image(40, 20, 3);
	set_fits_date(fit2, 200);
	CHECK(!ser_write_frame_from_fit(ser, fit2, 1), "writing image\n");
	fits *fit3 = create_image(40, 20, 3);
	set_fits_date(fit3, 300);
	CHECK(!ser_write_frame_from_fit(ser, fit3, 2), "writing image\n");
	CHECK(!ser_write_and_close(ser), "close file\n");

	CHECK(!ser_open_file(TMP_FILE7, ser), "reopen\n");
	CHECK(ser->color_id == SER_RGB, "wrong image color id\n");
	CHECK(ser->image_width == 40, "wrong image width\n");
	CHECK(ser->image_height == 20, "wrong image width\n");
	CHECK(ser->frame_count == 3, "wrong number of frames\n");
	CHECK(ser->ts, "no date information in SER\n");
	CHECK(!strcmp(ser->observer, observer_str), "observer was not copied\n");
	CHECK(ser->date_utc == utc_time, "UTC time was not copied\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[0])) == 100,
			"first image date is wrong\n");
	CHECK(g_date_time_to_unix(ser_timestamp_to_date_time(ser->ts[1])) == 200,
			"second image date is wrong\n");
	free(ser);

	CHECK(!unlink(TMP_FILE7), "error unlinking file " TMP_FILE7);
	fprintf(stdout, "* test 7 passed *\n\n");
	return 0;
}

/* A CFA frame written to SER and read back must keep its data and be
 * debayered with the same pattern as the original FITS */
static int check_cfa_roundtrip(row_order_t pref, const char *row_order, int h, sensor_pattern pattern) {
	const int w = 20;
	size_t n = (size_t) w * h;
	com.pref.debayer.use_bayer_header = TRUE;
	com.pref.debayer.orientation = pref;

	fits *fit = create_image(w, h, 1);
	CHECK(fit, "image creation\n");
	strcpy(fit->keywords.bayer_pattern, filter_pattern[pattern]);
	strcpy(fit->keywords.row_order, row_order);
	for (size_t i = 0; i < n; i++)
		fit->data[i] = (WORD) i;
	sensor_pattern expected = get_validated_cfa_pattern(fit, FALSE, FALSE);
	WORD *orig = g_memdup2(fit->data, n * sizeof(WORD));

	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_create_file(TMP_FILE8, ser, TRUE, NULL), "create file\n");
	CHECK(!ser_write_frame_from_fit(ser, fit, 0), "writing image\n");
	CHECK(!ser_write_and_close(ser), "close file\n");
	CHECK(!ser_open_file(TMP_FILE8, ser), "reopen\n");
	fits *out = calloc(1, sizeof(fits));
	CHECK(!ser_read_frame(ser, 0, out, FALSE, FALSE), "reading image\n");
	sensor_pattern got = get_validated_cfa_pattern(out, FALSE, FALSE);
	int same_data = !memcmp(orig, out->data, n * sizeof(WORD));
	ser_close_file(ser);
	free(ser);
	clearfits(out);
	free(out);
	g_free(orig);
	unlink(TMP_FILE8);

	CHECK(same_data, "pref %d, ROWORDER '%s', h %d, %s: data changed\n",
			pref, row_order, h, filter_pattern[pattern]);
	CHECK(got == expected, "pref %d, ROWORDER '%s', h %d, %s: pattern %s instead of %s\n",
			pref, row_order, h, filter_pattern[pattern], filter_pattern[got], filter_pattern[expected]);
	return !(same_data && got == expected);
}

int test_ser_cfa_roundtrip() {
	const row_order_t prefs[] = { ROW_ORDER_HEADER_TOPDOWN, ROW_ORDER_HEADER_BOTTOMUP,
		ROW_ORDER_FORCE_TOPDOWN, ROW_ORDER_FORCE_BOTTOMUP };
	const char *row_orders[] = { "TOP-DOWN", "BOTTOM-UP", "" };
	const int heights[] = { 10, 11 };
	int failed = 0;
	for (int p = 0; p < G_N_ELEMENTS(prefs); p++)
		for (int r = 0; r < G_N_ELEMENTS(row_orders); r++)
			for (int h = 0; h < G_N_ELEMENTS(heights); h++)
				for (sensor_pattern pat = BAYER_FILTER_RGGB; pat <= BAYER_FILTER_GRBG; pat++)
					failed |= check_cfa_roundtrip(prefs[p], row_orders[r], heights[h], pat);
	if (!failed)
		fprintf(stdout, "* test 8 passed *\n\n");
	return failed;
}

/* rewriting a CFA SER frame to SER must keep the data and the Bayer pattern,
 * whatever the orientation preference */
static int check_cfa_ser_passthrough(row_order_t pref, int h, sensor_pattern pattern) {
	const int w = 20;
	size_t n = (size_t) w * h;
	com.pref.debayer.use_bayer_header = TRUE;
	com.pref.debayer.orientation = ROW_ORDER_HEADER_TOPDOWN;

	fits *fit = create_image(w, h, 1);
	CHECK(fit, "image creation\n");
	strcpy(fit->keywords.bayer_pattern, filter_pattern[pattern]);
	strcpy(fit->keywords.row_order, "TOP-DOWN");
	for (size_t i = 0; i < n; i++)
		fit->data[i] = (WORD) i;
	struct ser_struct *src = malloc(sizeof(struct ser_struct));
	ser_init_struct(src);
	CHECK(!ser_create_file(TMP_FILE9, src, TRUE, NULL), "create file\n");
	CHECK(!ser_write_frame_from_fit(src, fit, 0), "writing image\n");
	CHECK(!ser_write_and_close(src), "close file\n");

	com.pref.debayer.orientation = pref;
	CHECK(!ser_open_file(TMP_FILE9, src), "reopen\n");
	fits *frame = calloc(1, sizeof(fits));
	CHECK(!ser_read_frame(src, 0, frame, FALSE, FALSE), "reading image\n");
	WORD *orig = g_memdup2(frame->data, n * sizeof(WORD));
	struct ser_struct *dst = malloc(sizeof(struct ser_struct));
	ser_init_struct(dst);
	CHECK(!ser_create_file(TMP_FILE10, dst, TRUE, NULL), "create file\n");
	CHECK(!ser_write_frame_from_fit(dst, frame, 0), "writing image\n");
	CHECK(!ser_write_and_close(dst), "close file\n");
	CHECK(!ser_open_file(TMP_FILE10, dst), "reopen\n");
	fits *out = calloc(1, sizeof(fits));
	CHECK(!ser_read_frame(dst, 0, out, FALSE, FALSE), "reading image\n");
	int same_data = !memcmp(orig, out->data, n * sizeof(WORD));
	int same_color = src->color_id == dst->color_id && src->color_id != SER_MONO;
	ser_close_file(src);
	ser_close_file(dst);
	free(src);
	free(dst);
	clearfits(out);
	free(out);
	g_free(orig);
	unlink(TMP_FILE9);
	unlink(TMP_FILE10);

	CHECK(same_data, "pref %d, h %d, %s: data changed\n", pref, h, filter_pattern[pattern]);
	CHECK(same_color, "pref %d, h %d, %s: SER color id changed\n", pref, h, filter_pattern[pattern]);
	return !(same_data && same_color);
}

int test_ser_cfa_passthrough() {
	const row_order_t prefs[] = { ROW_ORDER_HEADER_TOPDOWN, ROW_ORDER_HEADER_BOTTOMUP,
		ROW_ORDER_FORCE_TOPDOWN, ROW_ORDER_FORCE_BOTTOMUP };
	const int heights[] = { 10, 11 };
	int failed = 0;
	for (int p = 0; p < G_N_ELEMENTS(prefs); p++)
		for (int h = 0; h < G_N_ELEMENTS(heights); h++)
			for (sensor_pattern pat = BAYER_FILTER_RGGB; pat <= BAYER_FILTER_GRBG; pat++)
				failed |= check_cfa_ser_passthrough(prefs[p], heights[h], pat);
	if (!failed)
		fprintf(stdout, "* test 9 passed *\n\n");
	return failed;
}

/* writes an 8-bit single frame SER whose pixel values are the colour (1 + RGB
 * layer) of the Bayer filter they sit under, rows stored top-down */
static int write_native_cfa_ser(const char *filename, int w, int h, sensor_pattern pattern) {
	const ser_color ids[] = { [BAYER_FILTER_RGGB] = SER_BAYER_RGGB, [BAYER_FILTER_BGGR] = SER_BAYER_BGGR,
		[BAYER_FILTER_GBRG] = SER_BAYER_GBRG, [BAYER_FILTER_GRBG] = SER_BAYER_GRBG };
	char header[SER_HEADER_LEN] = { 0 };
	gint32 fields[7] = { 0, ids[pattern], SER_LITTLE_ENDIAN, w, h, 8, 1 };
	memcpy(header, "LUCAM-RECORDER", 14);
	memcpy(header + 14, fields, sizeof(fields));
	BYTE *data = malloc((size_t) w * h);
	const char *cfa = filter_pattern[pattern];
	for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++) {
			char c = cfa[((y & 1) << 1) | (x & 1)];
			data[y * w + x] = 1 + (c == 'R' ? RLAYER : c == 'G' ? GLAYER : BLAYER);
		}
	FILE *f = fopen(filename, "wb");
	if (!f) {
		free(data);
		return 1;
	}
	int ret = fwrite(header, SER_HEADER_LEN, 1, f) != 1 || fwrite(data, (size_t) w * h, 1, f) != 1;
	fclose(f);
	free(data);
	return ret;
}

/* debayering a frame of a camera SER must find the right colour under each
 * pixel, whatever the orientation preference */
static int check_native_ser_colours(row_order_t pref, int h, sensor_pattern pattern) {
	const int w = 20;
	com.pref.debayer.use_bayer_header = TRUE;
	com.pref.debayer.orientation = pref;
	CHECK(!write_native_cfa_ser(TMP_FILE11, w, h, pattern), "writing SER\n");
	struct ser_struct *ser = malloc(sizeof(struct ser_struct));
	ser_init_struct(ser);
	CHECK(!ser_open_file(TMP_FILE11, ser), "open\n");
	fits *fit = calloc(1, sizeof(fits));
	CHECK(!ser_read_frame(ser, 0, fit, FALSE, FALSE), "reading image\n");
	ser_close_file(ser);
	free(ser);
	unlink(TMP_FILE11);

	BYTE cfa[36];
	int size = 0;
	int ok = !get_compiled_pattern(fit, cfa, &size, FALSE) && size == 2;
	for (int y = 0; ok && y < h; y++)
		for (int x = 0; ok && x < w; x++)
			ok = FC_array(y, x, cfa, size) + 1 == fit->data[y * w + x];
	clearfits(fit);
	free(fit);

	CHECK(ok, "pref %d, h %d, %s: wrong colours\n", pref, h, filter_pattern[pattern]);
	return !ok;
}

int test_ser_native_cfa_colours() {
	const row_order_t prefs[] = { ROW_ORDER_HEADER_TOPDOWN, ROW_ORDER_HEADER_BOTTOMUP,
		ROW_ORDER_FORCE_TOPDOWN, ROW_ORDER_FORCE_BOTTOMUP };
	const int heights[] = { 10, 11 };
	int failed = 0;
	for (int p = 0; p < G_N_ELEMENTS(prefs); p++)
		for (int h = 0; h < G_N_ELEMENTS(heights); h++)
			for (sensor_pattern pat = BAYER_FILTER_RGGB; pat <= BAYER_FILTER_GRBG; pat++)
				failed |= check_native_ser_colours(prefs[p], heights[h], pat);
	if (!failed)
		fprintf(stdout, "* test 10 passed *\n\n");
	return failed;
}

#ifdef WITH_MAIN
int main() {
	int retval = 0;
	retval |= test_ser_image_number();
	retval |= test_ser_image_overwrite();
	retval |= test_ser_image_sizes();
	retval |= test_ser_dates();
	retval |= test_ser_with_holes();
	retval |= test_ser_ooo_write();
	retval |= test_ser_create_from_copy();
	retval |= test_ser_cfa_roundtrip();
	retval |= test_ser_cfa_passthrough();
	retval |= test_ser_native_cfa_colours();
	if (retval)
		fprintf(stderr, "TESTS FAILED\n");
	else fprintf(stderr, "ALL TESTS PASSED\n");
	return retval;
}
#else // with criterion
Test(ser, image_number) { cr_assert(!test_ser_image_number()); }
Test(ser, image_overwrite) { cr_assert(!test_ser_image_overwrite()); }
Test(ser, image_sizes) { cr_assert(!test_ser_image_sizes()); }
Test(ser, dates) { cr_assert(!test_ser_dates()); }
Test(ser, holes) { cr_assert(!test_ser_with_holes()); }
Test(ser, out_of_order_write) { cr_assert(!test_ser_ooo_write()); }
Test(ser, create_from_copy) { cr_assert(!test_ser_create_from_copy()); }
Test(ser, cfa_roundtrip) { cr_assert(!test_ser_cfa_roundtrip()); }
Test(ser, cfa_passthrough) { cr_assert(!test_ser_cfa_passthrough()); }
Test(ser, native_cfa_colours) { cr_assert(!test_ser_native_cfa_colours()); }
#endif
