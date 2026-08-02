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

/* Values are indices into the reg_methods[] table in align_rgb.c */
typedef enum {
	RGBALIGN_PSF = 0,
	RGBALIGN_DFT = 1,
	RGBALIGN_GLOBAL = 2,
	RGBALIGN_KOMBAT = 3,
} rgb_align_method;

/* What a channel alignment needs to be repeated.  The method is the user's
 * choice; the area is not — get_the_registration_area() derives it from
 * com.selection, which is a live GUI state that no longer exists by the time
 * the history is replayed.  So the effective area is captured, in the same
 * spirit as synthstar's star list: the detection is delegated and re-runs
 * against whatever image the replay presents, but the frame it runs in is
 * pinned. */
struct rgb_align_data {
	destructor destroy_fn;
	int        method;      /* an rgb_align_method value */
	rectangle  area;        /* effective registration area */
	gboolean   have_area;
};

struct rgb_align_data *new_rgb_align_data(void);
void free_rgb_align_data(void *p);

const char *rgb_align_method_name(rgb_align_method m);
gboolean rgb_align_prerequisites_met(rgb_align_method m);

/* Align gfit's three channels.  @area_used, when non-NULL, receives the
 * registration area actually used, for the caller to record.
 * Returns 0 on success; non-zero on failure (pixels left unchanged). */
int rgb_align(rgb_align_method m, rectangle *area_used);

/* rgb_align + NDE provenance capture, shared by the RGB align menu entries
 * and the rgbalign command.  @undo_err is the return of the caller's
 * undo_save_state for this operation (the capture tags that undo entry on
 * success); @summary is the same label used for the undo entry.
 * Returns rgb_align's return value. */
int rgb_align_and_capture(rgb_align_method m, int undo_err, const char *summary);

int rgb_align_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

#endif /* SRC_COMPOSITING_ALIGN_RGB_H_ */
