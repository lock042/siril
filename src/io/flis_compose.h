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

#ifndef FLIS_COMPOSE_H
#define FLIS_COMPOSE_H

#include <glib.h>
#include "core/siril.h"   /* for fits */

/**
 * flis_render_layers:
 * @layers: GSList of flis_layer_t*, sorted ascending by layer_order.
 *          The first element defines the canvas dimensions; all visible
 *          layers in the list are composited into a freshly allocated
 *          fits* in float RGB.
 *
 * Returns: a heap-allocated 3-channel float fits* containing the
 *          composited image, or NULL on error.  Caller takes ownership
 *          and must clearfits()+free() the result.
 *
 * This is the kernel called by the merge-down, flatten, thumbnail-
 * generation, and (stage 3 fallback) display paths.  Stage 3's GPU
 * compose path uses it as the correctness oracle in pixel-equivalence
 * tests.  Stage 3 itself will also add region-based and BGRA8 entry
 * points for the tiled display pipeline; those land in this same file
 * when they are needed.
 */
fits *flis_render_layers(GSList *layers);

#endif /* FLIS_COMPOSE_H */
