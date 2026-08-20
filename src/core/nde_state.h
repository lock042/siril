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
#ifndef _NDE_STATE_H_
#define _NDE_STATE_H_

/**
 * \file nde_state.h
 * \brief What an NDE item's value IS: one type, carried everywhere.
 *
 * A layer's value is not its pixels alone.  Geometry operations — crop,
 * resample, rotate, mirror, binning, registration — change the pixels AND
 * move the layer on the canvas, and each move is relative to the previous
 * position.  So a replay that carried only pixels could rebuild them and
 * still have nowhere to put them: it would need to know where the layer was
 * when the chain started.
 *
 * That "where" used to travel beside the pixels, in a second hash table keyed
 * the same way as the checkpoint tables, so that call sites could keep their
 * signatures.  Five places had to keep the two tables in step (store, drop,
 * purge, rebind, and the loader), the LRU restart cache could not hold a
 * position at all — which cost every layer that had ever been cropped its
 * cached restart points — and "one record, N written values" had nowhere to
 * go.  One type ends all three.
 *
 * OWNERSHIP.  A state from nde_state_new() or from a store read OWNS its
 * pixels and is released with nde_state_free().  A state declared on the
 * stack and filled in by hand BORROWS them and is simply left to go out of
 * scope — that is how a caller hands existing pixels to something that only
 * reads them.  The two are told apart by where the state came from, so do
 * not nde_state_free() one you filled in yourself.
 *
 * POSITION.  @has_pos is FALSE for anything that has no place on a canvas: a
 * plain single image, a mask, a state whose chain never moved its layer.  It
 * is not "unknown" — it means the value genuinely has no position, and a
 * consumer must leave whatever position the layer already has alone.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;
typedef struct ffit fits;

typedef struct {
	fits     *pix;            /* the pixels; owned unless the state is a borrow */
	gint      pos_x, pos_y;   /* where they sit on the canvas */
	gboolean  has_pos;        /* FALSE = this value has no position at all */
} nde_state;

/** Wrap @pix_taken (taken over) as a positionless state.  NULL in, NULL out. */
nde_state *nde_state_new(fits *pix_taken);

/** Wrap @pix_taken (taken over) as a state at (@pos_x, @pos_y). */
nde_state *nde_state_new_at(fits *pix_taken, gint pos_x, gint pos_y);

/** Free @s and the pixels it owns.  NULL-safe. */
void nde_state_free(nde_state *s);

/** Hand the pixels out and free the wrapper — for a caller that wants only
 *  the image and has no use for the position.  NULL-safe (returns NULL). */
fits *nde_state_release(nde_state *s);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_STATE_H_ */
