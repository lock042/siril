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
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "core/siril.h"
#include "core/siril_date.h"
#include "io/ser.h"
#include "io/image_format_fits.h"

#ifdef WITH_MAIN
#define CHECK(cond, ...) \
	if (!(cond)) { \
		fprintf(stderr, __VA_ARGS__); \
		return 1; \
	}
#else
#define CHECK cr_expect

cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image
#endif

#ifndef _WIN32
#define TMP_FILE1 "/tmp/test_tmp1.ser"
#define TMP_FILE2 "/tmp/test_tmp2.ser"
#define TMP_FILE3 "/tmp/test_tmp3.ser"
#define TMP_FILE4 "/tmp/test_tmp4.ser"
#define TMP_FILE5 "/tmp/test_tmp5.ser"
#define TMP_FILE6 "/tmp/test_tmp6.ser"
#define TMP_FILE7 "/tmp/test_tmp7.ser"
#else
#define TMP_FILE1 ".\\test_tmp1.ser"
#define TMP_FILE2 ".\\test_tmp2.ser"
#define TMP_FILE3 ".\\test_tmp3.ser"
#define TMP_FILE4 ".\\test_tmp4.ser"
#define TMP_FILE5 ".\\test_tmp5.ser"
#define TMP_FILE6 ".\\test_tmp6.ser"
#define TMP_FILE7 ".\\test_tmp7.ser"
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
#endif
