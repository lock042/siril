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

#include <criterion/criterion.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include "core/siril.h"
#include "core/proto.h"
#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#endif
#include "core/siril_date.h"
#include "io/gps_parser.h"
#include "io/image_format_fits.h"
#include "algos/geometry.h"

/************ download test files **************/

#define TEST_FILES_BUCKET_URL "https://siril-share-public.s3.rbx.io.cloud.ovh.net"
#define TEST_FILES_DOWNLOAD_DIR "siril/tests" // in user cache dir
#define LOCK_FILE_FOR_DOWNLOAD "/tmp/siril-test-download.lock"

#define LOCK(retval) \
	int fd_lock = open(LOCK_FILE_FOR_DOWNLOAD, O_RDWR | O_CREAT, 0640); \
	if (fd_lock < 0) { \
		siril_log_debug("could not create lock file\n"); \
		return retval; \
	} \
	flock(fd_lock, LOCK_EX)

#define UNLOCK \
	flock(fd_lock, LOCK_UN); \
	close(fd_lock);

static size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
	size_t written = fwrite(ptr, size, nmemb, stream);
	return written;
}

static gboolean directory_created = FALSE;

/* each test runs in a separate process, we cannot use a mutex to protect check_or_download_test_file(),
 * we use a file lock instead, which may not be portable. The other solution would be to make a function
 * with all the test function and call this one instead in the single Test() */
static void init_download() {
	LOCK();
	const gchar *cache_dir_path = g_get_user_cache_dir();
	gchar *test_files_dir_path = g_strdup_printf("%s/%s", cache_dir_path, TEST_FILES_DOWNLOAD_DIR);
	if (g_mkdir_with_parents(test_files_dir_path, 0700)) {
		siril_log_debug("failed to create dir %s\n", test_files_dir_path);
	} else {
		directory_created = TRUE;
	}
	g_free(test_files_dir_path);
	UNLOCK;
}

static gchar *check_or_download_test_file(const char *filename) {
	LOCK(NULL);
	const gchar *cache_dir_path = g_get_user_cache_dir();
	gchar *file_path = g_strdup_printf("%s/%s/%s", cache_dir_path, TEST_FILES_DOWNLOAD_DIR, filename);
	if (is_readable_file(file_path)) {
		UNLOCK;
		return file_path;
	}
	if (!directory_created) {
		siril_log_debug("cannot download to expected directory which creation failed\n");
		g_free(file_path);
		UNLOCK;
		return NULL;
	}
#ifdef HAVE_LIBCURL
	siril_log_debug("test file %s is missing, downloading...\n", file_path);
	FILE *dest_file = g_fopen(file_path, "wb");
	if (!dest_file) {
		perror("fopen");
		siril_log_debug("failed to create file %s for downloads\n", file_path);
		UNLOCK;
		return NULL;
	}
	gchar *file_url = g_strdup_printf("%s/%s", TEST_FILES_BUCKET_URL, filename);

	CURL *curl;
	curl = curl_easy_init();
	if (curl) {
		curl_easy_setopt(curl, CURLOPT_URL, file_url);
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, dest_file);
		curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L);
		CURLcode res = curl_easy_perform(curl);
		curl_easy_cleanup(curl);
		if (res != CURLE_OK) {
			siril_log_debug("failed to download from %s\n", file_url);
			g_unlink(file_path);
			g_free(file_path);
			file_path = NULL;
		}
	}
	fclose(dest_file);
	g_free(file_url);
	//siril_log_debug("test file %s was downloaded\n", file_path);
#else
	siril_log_debug("test file %s is missing and cannot be doawnloaded (libcurl missing)\n", file_path);
#endif
	UNLOCK;
	return file_path;
}

/***********************************************/

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

void test_rsgps_data_loading() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);
	cr_assert(fit.keywords.gps_data);
	cr_assert(!fit.keywords.date_and_exp_from_gps);
	clearfits(&fit);

	file_path = check_or_download_test_file("gps_global_shutter_not_locked.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);
	cr_assert(!fit.keywords.gps_data);
	cr_assert(!fit.keywords.date_and_exp_from_gps);
	clearfits(&fit);
}

void test_metadata_reading() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_global_shutter_not_locked.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);
	cr_assert(!fit.keywords.gps_data);
	cr_assert(!fit.keywords.date_and_exp_from_gps);
	GDateTime *original_dateobs = g_date_time_ref(fit.keywords.date_obs);
	double original_exposure = fit.keywords.exposure;
	long original_naxes1 = fit.naxes[1];
	int original_ry = fit.ry;

	int crop_rows = 1;
	struct generic_seq_args args = { 0 };
	args.user = GINT_TO_POINTER(crop_rows);
	int retval = gps_extract_image_hook(&args, 0, 0, &fit, NULL, MULTI_THREADED);
	cr_assert(!retval);

	double new_exposure = fit.keywords.exposure;
	GDateTime *new_dateobs = g_date_time_ref(fit.keywords.date_obs);
	gchar *new_dateobs_str = date_time_to_FITS_date(new_dateobs);
	siril_log_debug("exp %f -> %f\n", original_exposure, new_exposure);
	siril_log_debug("date-obs %s -> %s\n", date_time_to_FITS_date(original_dateobs), new_dateobs_str);
	cr_assert(fabs(original_exposure - new_exposure) < original_exposure * 0.001);
	cr_assert(timediff_in_s(original_dateobs, new_dateobs) < 1.0);
	cr_assert(!strcmp("2023-08-31T20:20:12.854125", new_dateobs_str));
	// image metadata has been modified using the GPS metadata
	cr_assert(fit.ry == original_ry - 1);
	cr_assert(fit.naxes[1] == original_naxes1 - 1);
	cr_assert(!fit.keywords.gps_data);
	// gps_extract_image_hook replaces data even if it's not locked
	cr_assert(fit.keywords.date_and_exp_from_gps);
	clearfits(&fit);
	// TODO: some headers problems remain after using the gps command:
	// - IMAGETYP and TELESCOP get some extra spaces
	// - GPS_EUTC and GPS_EFLG are removed because managed by the rolling shutter function
}

void test_rsgps_timestamp_computation() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	int row = 0;
	enum timestamp_type moment = EXP_START;
	GDateTime *date_start0 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	//date_str = date_time_to_FITS_date(date_start);
	//cr_assert(!strcmp("", date_str));

	moment = EXP_MIDDLE;
	GDateTime *date_middle0 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	cr_assert(fabs(timediff_in_s(date_start0, date_middle0) - fit.keywords.exposure * 0.5) < 1e-3);

	moment = EXP_END;
	GDateTime *date_end0 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	cr_assert(fabs(timediff_in_s(date_start0, date_end0) - fit.keywords.exposure) < 1e-3);

	row = 3193;	// this is on-screen row, so effectively the top of the image
	moment = EXP_START;
	GDateTime *date_start1 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	//date_str = date_time_to_FITS_date(date_start);
	//cr_assert(!strcmp("", date_str));

	moment = EXP_MIDDLE;
	GDateTime *date_middle1 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	cr_assert(fabs(timediff_in_s(date_start1, date_middle1) - fit.keywords.exposure * 0.5) < 1e-3);

	moment = EXP_END;
	GDateTime *date_end1 = get_timestamp_for_pixel(fit.keywords.gps_data, moment, 0, row);
	cr_assert(fabs(timediff_in_s(date_start1, date_end1) - fit.keywords.exposure) < 1e-3);

	double rolling_shutter_duration = timediff_in_s(date_start1, date_start0);
	siril_log_debug("rolling shutter duration: %f\n", rolling_shutter_duration);
	cr_assert(rolling_shutter_duration > 0.2);
	cr_assert(rolling_shutter_duration < 0.3);

	clearfits(&fit);
}

void test_rsgps_binning() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	int row = 500;
	GDateTime *date_before = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, row);
	cr_assert(fit.keywords.gps_data->binning == 2);

	fits_binning(&fit, 2, TRUE);

	cr_assert(fit.keywords.gps_data->binning == 4);
	GDateTime *date_after = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, row / 2);
	siril_log_debug("timestamp change after binning: %s -> %s, diff = %f\n",
			date_time_to_FITS_date(date_before), date_time_to_FITS_date(date_after),
			timediff_in_s(date_before, date_after));
	cr_assert(fabs(timediff_in_s(date_before, date_after)) < 1e-4);

	clearfits(&fit);
}

void test_rsgps_flip() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	GDateTime *date_before0 = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 0);
	GDateTime *date_before1 = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 3193);
	cr_assert(fit.keywords.gps_data->top_down);
	mirrorx(&fit, FALSE);
	cr_assert(!fit.keywords.gps_data->top_down);
	GDateTime *date_after0 = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 0);
	GDateTime *date_after1 = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 3193);
	cr_assert(fabs(timediff_in_s(date_before0, date_after1)) < 1e-6);
	cr_assert(fabs(timediff_in_s(date_before1, date_after0)) < 1e-6);

	clearfits(&fit);
}

void test_rsgps_crop() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	GDateTime *date_before = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 500);
	cr_assert(fit.keywords.gps_data->crop_offset_y == 0);
	rectangle area = { 300, 500, 200, 200 };
	crop(&fit, &area);
	siril_log_debug("new offset: %d\n", fit.keywords.gps_data->crop_offset_y);
	cr_assert(fit.keywords.gps_data->crop_offset_y == 3194 - 500 - 200);
	GDateTime *date_after = get_timestamp_for_pixel(fit.keywords.gps_data, EXP_MIDDLE, 0, 0);
	cr_assert(fabs(timediff_in_s(date_before, date_after)) < 1e-6);

	clearfits(&fit);
}

void test_non_gps_images() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("B_stacked.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	cr_assert(!fit.keywords.gps_data);
	cr_assert(!fit.keywords.date_and_exp_from_gps);

	struct _qhy_struct qhy_header = { 0 };
	int retval = parse_gps_image(&fit, &qhy_header);
	cr_assert(retval);

	retval = parse_gps_from_header(&fit, file_path, &qhy_header);
	cr_assert(retval);

	g_free(file_path);
	clearfits(&fit);
}


TestSuite(qhy_gps, .init = init_download);
Test(qhy_gps, load) { test_rsgps_data_loading(); }
Test(qhy_gps, metadata) { test_metadata_reading(); }
Test(qhy_gps, rs_timestamp) { test_rsgps_timestamp_computation(); }
Test(qhy_gps, binning) { test_rsgps_binning(); }
Test(qhy_gps, flip) { test_rsgps_flip(); }
Test(qhy_gps, crop) { test_rsgps_crop(); }
Test(qhy_gps, normal_image) { test_non_gps_images(); }

