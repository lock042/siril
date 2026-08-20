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

/* The NDE unit of currency — see nde_state.h for what it is and why. */

#include "core/siril.h"
#include "core/nde_state.h"
#include "io/image_format_fits.h"   /* clearfits */

nde_state *nde_state_new(fits *pix_taken) {
	if (!pix_taken)
		return NULL;
	nde_state *s = g_new0(nde_state, 1);
	s->pix = pix_taken;
	return s;
}

nde_state *nde_state_new_at(fits *pix_taken, gint pos_x, gint pos_y) {
	nde_state *s = nde_state_new(pix_taken);
	if (s) {
		s->pos_x   = pos_x;
		s->pos_y   = pos_y;
		s->has_pos = TRUE;
	}
	return s;
}

void nde_state_free(nde_state *s) {
	if (!s)
		return;
	if (s->pix) {
		clearfits(s->pix);
		free(s->pix);
	}
	g_free(s);
}

fits *nde_state_release(nde_state *s) {
	if (!s)
		return NULL;
	fits *pix = s->pix;
	g_free(s);
	return pix;
}
