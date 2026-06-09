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

#include "core/siril.h"
#include "opencv/opencv.h"
#include "clahe.h"

/* The actual CLAHE processing hook */
int clahe_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	clahe_params *params = (clahe_params *)args->user;
	if (!params)
		return 1;

	return cvClahe(fit, params->clip, params->tileSize);
}

gchar *clahe_log_hook(gpointer p, log_hook_detail detail) {
	clahe_params *params = (clahe_params *) p;
	if (!params) return NULL;
	gchar *message = g_strdup_printf(_("CLAHE (size=%d, clip=%.2f)"), params->tileSize, params->clip);
	return message;
}

/* GUI callbacks moved to src/gui/clahe.c */
