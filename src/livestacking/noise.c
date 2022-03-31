/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include <glib.h>
#include <stdlib.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "algos/noise.h"
#include "algos/statistics.h"

static GThread *thread;

struct _noise_data {
	double noise;
};

/* starting noise in a new thread */
gpointer noise_worker(gpointer arg) {
	double noise[3];
	int nb_channels = (int)gfit.naxes[2];
	for (int chan = 0; chan < nb_channels; chan++) {
		imstats *stat = statistics(NULL, -1, &gfit, chan, NULL, STATS_SIGMEAN, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return NULL;
		}
		noise[chan] = stat->bgnoise;
		free_stats(stat);
	}

	struct _noise_data *retval = calloc(1, sizeof(struct _noise_data));
	retval->noise = noise[0];
	return retval;
}

/* starts computing bgnoise on gfit, make sure it is not modified */
void bgnoise_async() {
	if (thread) {
		siril_debug_print("bgnoise request ignored, still running\n");
		return;
	}
	thread = g_thread_new("bgnoise", noise_worker, NULL);
}

double bgnoise_await() {
	if (!thread)
		return -1.0;
	gpointer retval = g_thread_join(thread);
	thread = NULL;
	if (!retval)
		return -1.0;
	double value = ((struct _noise_data *)retval)->noise;
	free(retval);
	if (gfit.type == DATA_FLOAT)
		value *= USHRT_MAX_DOUBLE;
	return value;
}

