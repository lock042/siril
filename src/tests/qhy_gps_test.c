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
#include "core/siril.h"
#include "core/proto.h"
#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#endif
//#include "core/siril_date.h"
#include "io/gps_parser.h"
#include "io/image_format_fits.h"

/************ download test files **************/

#define TEST_FILES_BUCKET_URL "https://siril-share-public.s3.rbx.io.cloud.ovh.net"
#define TEST_FILES_DOWNLOAD_DIR "siril/tests" // in user cache dir

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

gchar *check_or_download_test_file(const char *filename) {
	const gchar *cache_dir_path = g_get_user_cache_dir();
	gchar *test_files_dir_path = g_strdup_printf("%s/%s", cache_dir_path, TEST_FILES_DOWNLOAD_DIR);
	gchar *file_path = g_strdup_printf("%s/%s", test_files_dir_path, filename);
	if (is_readable_file(file_path))
		return file_path;
	if (g_mkdir_with_parents(test_files_dir_path, 0700)) {
		siril_debug_print("failed to create dir %s\n", test_files_dir_path);
		g_free(test_files_dir_path);
		g_free(file_path);
		return NULL;
	}
#ifdef HAVE_LIBCURL
	siril_debug_print("test file %s is missing, downloading...\n", file_path);
	FILE *dest_file = g_fopen(file_path, "wb");
	if (!dest_file) {
		perror("fopen");
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
		/* always cleanup */
		curl_easy_cleanup(curl);
		if (res != CURLE_OK) {
			siril_debug_print("failed to download from %s\n", file_url);
			g_unlink(file_path);
			g_free(file_path);
			file_path = NULL;
		}
	}
	fclose(dest_file);
	g_free(file_url);
#else
	siril_debug_print("test file %s is missing and cannot be doawnloaded (libcurl missing)\n", file_path);
#endif
	return file_path;
}

/***********************************************/

cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits *gfit;	// currently loaded image

void test_rsgps_data_loading() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("gps_target_rate_tracking_rolling_shutter_1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);
	cr_assert(fit.keywords.gps_data);
	clearfits(&fit);

	file_path = check_or_download_test_file("gps_global_shutter_not_locked.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);
	cr_assert(!fit.keywords.gps_data);
}

Test(qhy_gps, load) { test_rsgps_data_loading(); }

