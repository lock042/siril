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
};

/* Allocator and destructor functions */
struct unpurpleargs *new_unpurple_args();
void free_unpurple_args(void *args);

/* Image processing hook */
int unpurple_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
gchar *unpurple_log_hook(gpointer p, log_hook_detail detail);
/* Idle functions */
gboolean unpurple_preview_idle(gpointer p);
gboolean unpurple_apply_idle(gpointer p);

void apply_unpurple_cancel();
int generate_binary_starmask(fits *fit, fits **star_mask, double threshold);

#endif /* SRC_FILTERS_CA_H_ */
