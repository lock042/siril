/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "core/command.h"
#include "drizzle/driz_portability.h"
#include "drizzle/cdrizzleutil.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzlebox.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

void on_apply_drizzle_clicked(GtkButton *button, gpointer user_data) {
	struct driz_args_t *driz = calloc(1, sizeof(struct driz_args_t));
	driz->seq = &com.seq;
	driz->reference_image = sequence_find_refimage(&com.seq);
	driz->is_bayer = FALSE;
	driz->use_wcs = FALSE;
	driz->scale = 1.f;
	apply_drizzle(driz);
}


