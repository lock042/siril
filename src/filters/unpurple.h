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

#ifndef SRC_FILTERS_CA_H_
#define SRC_FILTERS_CA_H_

#include "core/siril.h"

struct unpurpleargs {
	destructor destroy_fn;  // Must be first member
	fits *fit;  // just a reference, not freed
	fits *starmask;
	gboolean starmask_needs_freeing;
	double mod_b;
	double thresh;
	gboolean withstarmask;
	gboolean verbose;
	gboolean applying;
	/* NDE Convention 2 provenance of the mask's star set (see synthstar.h).
	 * EXPLICIT (star_auto FALSE): stars_blob pins the consumed com.stars.
	 * DELEGATED (star_auto TRUE): star_conf records the auto-detection params,
	 * re-run at replay; stars_blob unused. */
	gchar *stars_blob;
	gboolean star_auto;
	star_finder_params star_conf;
};

/* Allocator and destructor functions */
struct unpurpleargs *new_unpurple_args();
void free_unpurple_args(void *args);

/* Image processing hook */
int unpurple_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
gchar *unpurple_log_hook(gpointer p, log_hook_detail detail);

void apply_unpurple_cancel();
/* Build a binary star mask from the star set.  Provenance outputs (any may be
 * NULL) record NDE replay info: on return *auto_out tells whether the set was
 * auto-detected (DELEGATED) rather than taken from com.stars (EXPLICIT);
 * *stars_blob_out receives the effective list (EXPLICIT only, caller g_free);
 * *conf_out receives the detection params used (DELEGATED only).
 * @conf_override (non-NULL only during DELEGATED replay) forces auto-detection
 * with those params installed, ignoring com.stars. */
int generate_binary_starmask(fits *fit, fits **star_mask, double threshold,
                             gchar **stars_blob_out, gboolean *auto_out,
                             star_finder_params *conf_out,
                             const star_finder_params *conf_override);

#endif /* SRC_FILTERS_CA_H_ */
