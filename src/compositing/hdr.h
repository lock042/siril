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

#ifndef SRC_COMPOSITING_HDR_H_
#define SRC_COMPOSITING_HDR_H_

#include <glib.h>
#include "core/siril.h"

#define HDR_DEFAULT_FEATHER 4

struct hdr_params {
	int feather;		// growth and transition width of the saturation masks, in pixels
};

struct hdr_args {
	GSList *files;		// gchar *, at least 2
	struct hdr_params params;
};

void hdr_params_init(struct hdr_params *params);
void free_hdr_args(struct hdr_args *args);

fits *hdr_compose(fits **images, const char **names, int nb, const struct hdr_params *params);
gpointer hdr_worker(gpointer p);

#endif /* SRC_COMPOSITING_HDR_H_ */
