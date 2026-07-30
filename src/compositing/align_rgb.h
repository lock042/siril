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

#ifndef SRC_COMPOSITING_ALIGN_RGB_H_
#define SRC_COMPOSITING_ALIGN_RGB_H_

#include "core/siril.h"          /* fits, rectangle, destructor */
#include "core/processing.h"     /* struct generic_img_args */

/* What a channel alignment needs to be repeated.  The method is the user's
 * choice; the area is not — get_the_registration_area() derives it from
 * com.selection, which is a live GUI state that no longer exists by the time
 * the history is replayed.  So the effective area is captured, in the same
 * spirit as synthstar's star list: the detection is delegated and re-runs
 * against whatever image the replay presents, but the frame it runs in is
 * pinned. */
struct rgb_align_data {
	destructor destroy_fn;
	int        method;      /* index into the internal registration method table */
	rectangle  area;        /* effective registration area */
	gboolean   have_area;
};

struct rgb_align_data *new_rgb_align_data(void);
void free_rgb_align_data(void *p);

/* Align gfit's three channels.  @area_used, when non-NULL, receives the
 * registration area actually used, for the caller to record.
 * Returns 0 on success; non-zero on failure (pixels left unchanged). */
int rgb_align(int m, rectangle *area_used);

int rgb_align_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

#endif /* SRC_COMPOSITING_ALIGN_RGB_H_ */
