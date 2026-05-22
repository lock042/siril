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
 */

/*
 * flis_gpu_compose.h — per-layer GdkTexture cache + GSK push_blend
 * composition path (stage 3.2).
 *
 * For FLIS files whose layer stack consists entirely of "simple" cases
 * (no sparse positioning, no groups, no FLIS_BLEND_CHROMA), each visible
 * layer is uploaded once as a canvas-sized BGRA8 texture with stretch
 * and mono tint baked in.  The display snapshot then composites them
 * via nested gtk_snapshot_push_blend / push_opacity / push_mask nodes
 * — the GPU does the per-pixel blend arithmetic rather than the CPU.
 *
 * Property-only changes (opacity, visibility, blend mode swap on an
 * existing layer) require no texture rebuild — they just change the
 * snapshot tree on the next frame.  Pixel changes invalidate one
 * layer's texture; stretch changes invalidate all (the stretch is
 * baked in).
 *
 * Layers with active lmask are supported via a second per-layer
 * texture used as a luminance mask.
 *
 * Stage 3.3 will replace the canvas-sized texture with per-tile
 * textures that materialise lazily, mirroring the existing
 * plain-FITS tile pipeline.  The public API here stays stable across
 * that refactor.
 */

#ifndef FLIS_GPU_COMPOSE_H
#define FLIS_GPU_COMPOSE_H

#include <gtk/gtk.h>
#include <glib.h>

/* Returns TRUE if every layer in @layers can be rendered via the
 * GPU compose path.  When this returns FALSE the caller must fall
 * back to the §3.1 CPU composite.  Reasons to fall back:
 *   • Any layer uses FLIS_BLEND_CHROMA (LRGB — no GSK equivalent)
 *   • Any layer is sparse (position_x/y != 0 OR layer extent != canvas)
 *   • Any layer belongs to a group (groups are §3.3 territory)
 *   • Empty layer list */
gboolean flis_gpu_compose_compatible(GSList *layers);

/* Render @layers into @snapshot inside the rectangle @dst_rect (in
 * widget-space pixels).  The canvas-to-widget transform must already
 * be on @snapshot's stack (caller does gtk_snapshot_save /
 * push_transform around this call).
 *
 * Requires flis_gpu_compose_compatible(layers) == TRUE; callers
 * must check before invoking.  Builds and caches textures lazily.
 * Stretch lo/hi are read from gui.lo / gui.hi at call time and
 * trigger a cache rebuild when they change. */
void flis_gpu_compose_render(GtkSnapshot *snapshot,
                              GSList *layers,
                              guint canvas_w, guint canvas_h,
                              const graphene_rect_t *dst_rect,
                              GskScalingFilter filter);

/* Invalidate the cached textures belonging to one layer.  Call after
 * a layer's pixel data, tint, or lmask data changes.  No-op if the
 * layer was never cached. */
void flis_gpu_compose_invalidate_layer(gint item_id);

/* Invalidate every cached texture.  Call when the stretch state
 * (gui.lo / gui.hi) changes, the layer stack order changes, or the
 * display mode (HISTEQ / STF / LINEAR) changes.  Cheaper than
 * flis_gpu_compose_free_all because it keeps the cache slot structs
 * alive; only the underlying GdkTextures are released. */
void flis_gpu_compose_invalidate_all(void);

/* Release every GPU resource owned by the cache.  Called on image
 * close.  Idempotent. */
void flis_gpu_compose_free_all(void);

#endif /* FLIS_GPU_COMPOSE_H */
