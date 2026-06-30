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
#include <stdio.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "io/image_format_fits.h"
#include "core/command.h"
#include "algos/streaks.h"
#include "core/processing.h"
#include "download_files.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

void test_streak_detection_gps() {
	initialize_default_settings();
	processing_system_init();
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("streaks/streak1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	gchar *result_basename = g_strdup("/tmp/result");
	gchar *result_file = g_strdup_printf("%s.streaks", result_basename);
	g_unlink(result_file);

	detect_streaks_async(&fit, 200, FALSE, 0, result_basename);

	waiting_for_thread();
	siril_log_debug("wait over\n");
	cr_assert(g_file_test(result_file, G_FILE_TEST_EXISTS));

	// TODO read content or return the date struct from the processing function
}

TestSuite(streaks, .init = init_download);
Test(streaks, detect) { test_streak_detection_gps(); }
