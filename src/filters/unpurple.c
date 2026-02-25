/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "gui/callbacks.h"
#include "algos/colors.h"
#include "filters/synthstar.h"
#include "filters/unpurple.h"

/*****************************************************************************
 *      U N P U R P L E   A L L O C A T O R   A N D   D E S T R U C T O R    *
 ****************************************************************************/

/* Allocator for unpurpleargs */
struct unpurpleargs *new_unpurple_args() {
	struct unpurpleargs *args = calloc(1, sizeof(struct unpurpleargs));
	if (args) {
		args->destroy_fn = free_unpurple_args;
	}
	return args;
}

/* Destructor for unpurpleargs */
void free_unpurple_args(void *ptr) {
	struct unpurpleargs *args = (struct unpurpleargs *)ptr;
	if (!args)
		return;

	if (args->starmask_needs_freeing && args->starmask) {
		clearfits(args->starmask);
		free(args->starmask);
		args->starmask = NULL;
	}
	free(ptr);
}

// TODO: Perhaps we still need a better purple detector?
static gboolean is_purple(float red, float green, float blue) {
	float h, s, v;
	rgb_to_hsvf(red, green, blue, &h, &s, &v);

	return (h >= 0.40f && h <= 0.99f) && (s >= 0.0f) && (v >= 0.0f);
}

int generate_binary_starmask(fits *fit, fits **star_mask, double threshold) {
	gboolean stars_needs_freeing = FALSE;
	psf_star **stars = NULL;
	int channel = 1;

	int nb_stars = starcount(com.stars);
	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;

	// Do we have stars from Dynamic PSF or not?
	if (nb_stars < 1) {
		image *input_image = NULL;
		input_image = calloc(1, sizeof(image));
		input_image->fit = fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;
		stars = peaker(input_image, channel, &com.pref.starfinder_conf, &nb_stars,
						NULL, FALSE, FALSE, 0, com.pref.starfinder_conf.profile, com.max_thread);
		free(input_image);
		stars_needs_freeing = TRUE;
	} else {
		stars = com.stars;
		stars_needs_freeing = FALSE;
	}

	if (starcount(stars) < 1) {
		siril_log_color_message(_("No stars detected in the image.\n"), "red");
		return -1;
	}

	siril_log_message(_("Creating binary star mask for %d stars...\n"), nb_stars);
	if (new_fit_image(star_mask, dimx, dimy, 1, DATA_USHORT)) {
		return -1;
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) if (com.max_thread > 1)
#endif
	for (size_t i = 0; i < dimx * dimy; i++) {
		(*star_mask)->pdata[0][i] = 0;
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread)
#endif
	for (int n = 0; n < nb_stars; n++) {
		int size = (int) 2 * max(stars[n]->fwhmx, stars[n]->fwhmy);

		// The fringe factor tries to adjust the star size to allow for the purple fringe
		// The threshold slider decides the scaling factor
		double base_fringe = 10.0;
		double fringe = base_fringe + threshold * pow(size, 1.5);
		size = sqrt(pow(size, 2) + pow(fringe, 2));

		if (size % 2 == 0)
			size++;

		int x0 = (int)(stars[n]->xpos - size / 2);
		int y0 = (int)((dimy - stars[n]->ypos) - size / 2);
		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				int px = x0 + x;
				int py = y0 + y;
				if (px >= 0 && px < dimx && py >= 0 && py < dimy) {
					int idx = py * dimx + px;
					if (idx >= 0 && idx < count) {
						double dx = x - (size / 2.0);
						double dy = y - (size / 2.0);
						double distance = sqrt(dx * dx + dy * dy);
						int is_star = (distance <= (size / 2.0)) ? 1 : 0;
						(*star_mask)->pdata[0][idx] = is_star ? USHRT_MAX : (*star_mask)->pdata[0][idx];
					}
				}
			}
		}
	}

	if (stars_needs_freeing)
		free_fitted_stars(stars);

	return 0;
}

gchar *unpurple_log_hook(gpointer p, log_hook_detail detail) {
	struct unpurpleargs *args = (struct unpurpleargs*) p;
	gchar *message = g_strdup_printf(_("Unpurple mod: %.2f, threshold: %.2f, starmask: %s"),
			args->mod_b, args->thresh, args->withstarmask ? _("true"): _("false"));
	return message;
}

//TODO improve this filter!
static int unpurple_filter(struct unpurpleargs *args) {
	fits *fit = args->fit;
	fits *starmask = args->starmask;
	double mod_b = args->mod_b;
	double thresh = args->thresh;
	gboolean withstarmask = args->withstarmask;

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

	if (fit == gfit && args->applying && !com.script) {
		populate_roi();
	}

	if (fit == gfit && args->applying) {
		siril_log_color_message(_("Unpurple filter applied: mod_b=%.3f, threshold=%.3f, withstarmask=%d\n"), "green",
								args->mod_b, args->thresh, args->withstarmask);
	}

	return 0;
}

/* The actual unpurple processing hook */
int unpurple_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct unpurpleargs *params = (struct unpurpleargs *)args->user;
	if (!params)
		return 1;
	params->fit = fit;
	return unpurple_filter(params);
}

/* Idle function for preview updates */
gboolean unpurple_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	// Free using the generic cleanup which will call the destructor
	free_generic_img_args(args);
	return FALSE;
}

/* Idle function for final application */
gboolean unpurple_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		notify_gfit_modified();
	}
	// Free using the generic cleanup which will call the destructor
	free_generic_img_args(args);
	return FALSE;
}
