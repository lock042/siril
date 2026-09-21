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
#include "core/siril.h"
#include "core/proto.h"
#include "io/image_format_fits.h"
#include "core/processing_thread.h"
#include "core/processing.h"
#include "filters/banding.h"
#include "filters/cosmetic_correction.h"
#include "core/op_descriptors.h"
#include "algos/statistics.h"
#include "download_files.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

static int apply_banding(fits *fit, double sigma, bool vertical) {
	struct banding_data *params = new_banding_data();
	params->protect_highlights = sigma > 0;
	params->amount = 1.0;
	params->sigma = sigma;
	params->vertical = vertical;
	params->seqEntry = NULL;
	params->seq = NULL;
	params->fit = NULL;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = fit;
	args->op = &op_desc_banding;
	args->idle_function = NULL; // Use default idle function for command-line
	args->command_updates_gfit = FALSE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = FALSE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	int retval = banding_single_image_hook(args, fit, com.max_thread);

	free_banding_data(params);
	free(args);
	return retval;
}

void test_banding_ushort() {
	initialize_default_settings();
	processing_system_init();
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("banding_example.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	imstats *stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double original_sigma = stats->sigma;
	free_stats(stats);

	cr_assert(!apply_banding(&fit, 0.0, FALSE));

	stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double after_horizontal = stats->sigma;
	free_stats(stats);
	cr_assert(after_horizontal < original_sigma);

	cr_assert(!apply_banding(&fit, 0.0, TRUE));
	stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double after_vertical = stats->sigma;
	free_stats(stats);
	cr_assert(after_vertical < after_horizontal);

	clearfits(&fit);
}

void test_banding_float() {
	initialize_default_settings();
	processing_system_init();
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("banding_example.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, TRUE) == 0);
	g_free(file_path);

	imstats *stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double original_sigma = stats->sigma;
	free_stats(stats);

	cr_assert(!apply_banding(&fit, 0.0, FALSE));

	stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double after_horizontal = stats->sigma;
	free_stats(stats);
	cr_assert(after_horizontal < original_sigma);

	cr_assert(!apply_banding(&fit, 0.0, TRUE));
	stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double after_vertical = stats->sigma;
	free_stats(stats);
	cr_assert(after_vertical < after_horizontal);

	clearfits(&fit);
}

void test_dead_column_cosmetic_correction() {
	initialize_default_settings();
	processing_system_init();
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("banding_example.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, TRUE) == 0);
	g_free(file_path);

	GFileIOStream *stream;
	GFile *file = g_file_new_tmp("cosme_XXXXXX.lst", &stream, NULL);
	cr_assert(file);
	const gchar *lst = "C 750 0\n";
	cr_assert(g_output_stream_write(g_io_stream_get_output_stream(G_IO_STREAM(stream)), lst, strlen(lst), NULL, NULL) > 0);
	cr_assert(g_io_stream_close(G_IO_STREAM(stream), NULL, NULL));

	imstats *stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	double original_sigma = stats->sigma;
	double original_mean = stats->mean;
	free_stats(stats);

	cr_assert(!apply_cosme_to_image(&fit, file, FALSE));

	stats = statistics(NULL, -1, &fit, 0, NULL, STATS_BASIC, MULTI_THREADED);
	cr_assert(stats->sigma < original_sigma);
	cr_assert(stats->mean > original_mean);
	free_stats(stats);
	g_object_unref(file);
	clearfits(&fit);
}

TestSuite(banding, .init = init_download);
Test(banding, remove_ushort) { test_banding_ushort(); }
Test(banding, remove_float) { test_banding_float(); }
TestSuite(cosme, .init = init_download);
Test(cosme, remove_column) { test_dead_column_cosmetic_correction(); }
