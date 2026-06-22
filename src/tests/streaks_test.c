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
#include "download_files.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

void test_streak_detection_gps() {
	fits fit = {0};
	gchar *file_path = check_or_download_test_file("streaks/streak1.fit");
	cr_assert(file_path);
	cr_assert(readfits(file_path, &fit, NULL, FALSE) == 0);
	g_free(file_path);

	float fwhm = 3.0f;
	simple_star_removal(&fit, 0, -0.1, &fwhm, &com.pref.starfinder_conf);

	// TODO: refactor the main function
	/*int expected_length = 250;
	struct streak_detection_conf *arg = calloc(1, sizeof(struct streak_detection_conf));
	arg->fit = &fit;
	arg->free_fit = FALSE;
	arg->im_idx = com.uniq ? 0 : com.seq.current;
	arg->filename = g_strdup("streak1");
	arg->layer = 0;
	arg->initial_segment_length = expected_length;
	arg->minimum_segment_length = round_to_int(expected_length * 0.7);
	arg->bright_target = FALSE;
	arg->display_streaks = FALSE;
	arg->use_idle = TRUE;
	arg->nb_threads = com.max_thread;
	arg->results = alloc_results(1);
	arg->fwhm = fwhm;
	int retval = GPOINTER_TO_INT(streak_detection_worker(arg));
	cr_assert(!retval);
	cr_assert(is_readable_file("streak1.streaks"));
	*/
}
