/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril_free-astro.org/index.php/Siril
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "gui/callbacks.h"
#include "algos/colors.h"
#include "filters/unpurple.h"

static gboolean end_unpurple(gpointer p) {
	stop_processing_thread();
	notify_gfit_modified();
	return FALSE;
}

gpointer unpurple_handler(gpointer args) {
	lock_roi_mutex();
	struct unpurpleargs *p = (struct unpurpleargs *)args;
	gpointer retval = unpurple(p);
	unlock_roi_mutex();
	if (!com.script)
		siril_add_idle(end_unpurple, NULL);
	return retval;
}

gpointer unpurple_filter(gpointer args) {
	struct unpurpleargs *p = (struct unpurpleargs *)args;
	gpointer retval = unpurple(p);
	if (!com.script)
		siril_add_idle(end_unpurple, NULL);
	return retval;
}

// TODO: Perhaps we still need a better purple detector?
static gboolean is_purple(float red, float green, float blue) {
	float h, s, v;
	rgb_to_hsvf(red, green, blue, &h, &s, &v);

	return (h >= 0.40f && h <= 0.99f) && (s >= 0.0f) && (v >= 0.0f);
}

//TODO improve this filter!
gpointer unpurple(gpointer p) {
	struct unpurpleargs *args = (struct unpurpleargs*) p;

	fits *fit = args->fit;
	fits *starmask = args->starmask;
	double mod_b = args->mod_b;
	double thresh = args->thresh;
	gboolean withstarmask = args->withstarmask;
	gboolean verbose = args->verbose;

	struct timeval t_start = { 0 }, t_end = { 0 };

	if (mod_b < 1) {
		for (size_t j = 0; j < fit->ry; j++) {
			size_t offset = j * fit->rx;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) if (com.max_thread > 1)
#endif
			for (size_t i = 0; i < fit->rx; i++) {
				float red, green, blue;
				if (fit->type == DATA_USHORT) {
					red = fit->pdata[RLAYER][i + offset] * INV_USHRT_MAX_SINGLE;
					green = fit->pdata[GLAYER][i + offset] * INV_USHRT_MAX_SINGLE;
					blue = fit->pdata[BLAYER][i + offset] * INV_USHRT_MAX_SINGLE;
				} else {
					red = fit->fpdata[RLAYER][i + offset];
					green = fit->fpdata[GLAYER][i + offset];
					blue = fit->fpdata[BLAYER][i + offset];
				}
				float luminance = 0.299f * red + 0.587f * green + 0.114f * blue;

				// Is this purple?
				if (is_purple(red, green, blue)) {

					// Only affect pixels that fall in the starmask or are greater than our background luminance threshold
					if ((withstarmask && starmask->pdata[RLAYER][i + offset] > 0) || (!withstarmask && luminance > thresh))  {

						// If we're modifying blue
						if (mod_b < 1.0) {
							// Scale colour values
							float target_blue = green * (0.587f / 0.114f);
							blue = blue * (float) mod_b + target_blue * (1.0f - (float) mod_b);

							// Scale luminance for blue channel. This is bloated from the chromatic aberration
							float luminance_factor = green;
							blue = blue * ((float) mod_b * luminance_factor + (float) mod_b);

							if (fit->type == DATA_USHORT) {
								fit->pdata[BLAYER][i + offset] = blue * USHRT_MAX_SINGLE;
							} else {
								fit->fpdata[BLAYER][i + offset] = blue;
							}
						}
					}
				}
			}
		}
	}

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	char log[90];
	sprintf(log, "Unpurple mod: %.2f, threshold: %.2f, withstarmask: %d", mod_b, thresh, withstarmask);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	if (args->for_final)
		populate_roi();
	siril_free(args);
	return GINT_TO_POINTER(0);
}
