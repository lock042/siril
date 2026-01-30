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

#include <math.h>
#include <string.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

#include "noise.h"

static GThread *thread;

static gboolean end_noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);

	if (args->display_start_end) {
		struct timeval t_end;
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	free(args);
	return FALSE;
}

// returns the argument if not freed
gpointer noise_worker(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	args->mean_noise = 0.0;

	if (args->display_start_end) {
		siril_log_color_message(_("Noise standard deviation: calculating...\n"),
				"green");
		gettimeofday(&args->t_start, NULL);
	}

	imstats *stats[3];
	int retval = compute_all_channels_statistics_single_image(args->fit, STATS_SIGMEAN, MULTI_THREADED, stats);
	for (int chan = 0; chan < args->fit->naxes[2]; chan++) {
		if (!retval) {
			args->bgnoise[chan] = stats[chan]->bgnoise;
			args->mean_noise += stats[chan]->bgnoise;
			double norm = stats[chan]->normValue;

			if (args->display_results) {
				if (args->fit->type == DATA_USHORT)
					siril_log_message(
							_("Background noise value (channel: #%d): %0.3f (%.3e)\n"),
							chan, args->bgnoise[chan], args->bgnoise[chan] / norm);
				else
					siril_log_message(
							_("Background noise value (channel: #%d): %0.3f (%.3e)\n"),
							chan, args->bgnoise[chan] * USHRT_MAX_DOUBLE, args->bgnoise[chan]);
			}
		}
		if (stats[chan])
			free_stats(stats[chan]);
	}

	if (retval) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		args->mean_noise = -1.0;
	} else {
		args->mean_noise = args->mean_noise / (double)args->fit->naxes[2];
		if (args->fit->type == DATA_FLOAT)
			args->mean_noise *= USHRT_MAX_DOUBLE;
	}

	if (args->use_idle) {
		siril_add_idle(end_noise, args);
		args = NULL;
	} else {
		args->retval = retval;
	}

	return args;
}

// called by the GUI or the command, uses the processing thread
void evaluate_noise_in_image() {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	set_cursor_waiting(TRUE);

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	struct noise_data *args = calloc(1, sizeof(struct noise_data));
	args->fit = &gfit;
	args->use_idle = TRUE;
	args->display_results = TRUE;
	args->display_start_end = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	if (!start_in_new_thread(noise_worker, args)) {
		free(args);
	}
}

// called in general from another function like stacking,
// bgnoise_await() has to be called to free resources
void bgnoise_async(fits *fit, gboolean display_values) {
	if (thread) {
		siril_debug_print("bgnoise request ignored, still running\n");
		return;
	}
	struct noise_data *args = calloc(1, sizeof(struct noise_data));
	args->fit = fit;
	args->use_idle = FALSE;
	args->display_start_end = FALSE;
	args->display_results = display_values;
	memset(args->bgnoise, 0.0, sizeof(double[3]));

	thread = g_thread_new("bgnoise", noise_worker, args);
}

double bgnoise_await() {
	if (!thread)
		return -1.0;
	struct noise_data *args = g_thread_join(thread);
	thread = NULL;
	if (!args)
		return -1.0;
	double value = args->mean_noise;
	free(args);
	return value;
}

