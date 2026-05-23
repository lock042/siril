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
 * flis_gpu_compose.h — per-layer per-tile GdkTexture cache + GSK
 * push_blend composition path (stages 3.2 + 3.3).
 *
 * For FLIS files whose layer stack can be expressed in GSK
 * primitives (translatable blend modes, base-layer covers canvas),
 * each visible layer is uploaded as a grid of per-tile BGRA8
 * textures with stretch and mono tint baked in.  Tiles materialise
 * lazily based on the visible-canvas rect passed to the render call,
 * so off-screen tiles on multi-100-Mpix mosaics never hit GPU memory.
 *
 * Property-only changes (opacity, visibility, blend mode swap on an
 * existing layer) require no texture rebuild — they just change the
 * snapshot tree on the next frame.  Pixel/tint/extent changes
 * invalidate one layer's tile grid; stretch changes invalidate all
 * (the stretch is baked into each tile).
 *
 * Layers with active lmask are supported via a single per-layer
 * lmask texture used as a luminance mask source.  Sparse layers are
 * placed at their FITS position_x/y via per-tile snapshot rects;
 * groups composite into sub-trees wrapped in
 * push_opacity(group_opacity) + push_blend(group_blend_mode).
 */

#ifndef FLIS_GPU_COMPOSE_H
#define FLIS_GPU_COMPOSE_H

#include <gtk/gtk.h>
#include <glib.h>

/* Returns TRUE if every layer in @layers can be rendered via the
 * GPU compose path.  When this returns FALSE the caller must fall
 * back to the §3.1 CPU composite.  Reasons to fall back:
 *   • Any layer uses FLIS_BLEND_CHROMA (LRGB — no GSK equivalent)
 *   • The base layer is sparse, or in a non-PASS_THROUGH group
 *   • A group's blend mode has no GSK equivalent
 *   • Empty layer list */
gboolean flis_gpu_compose_compatible(GSList *layers);

/* Render @layers into @snapshot inside the rectangle @dst_rect (in
 * widget-space pixels).  The canvas-to-widget transform must already
 * be on @snapshot's stack (caller does gtk_snapshot_save /
 * push_transform around this call).
 *
 * @visible_canvas is an optional culling hint in canvas image-space
 * (display top-down): tiles whose canvas rect does not intersect
 * this rect are skipped (no texture materialised, no GdkTexture
 * appended).  Pass NULL to draw all tiles unconditionally.
 *
 * Requires flis_gpu_compose_compatible(layers) == TRUE; callers
 * must check before invoking.  Builds and caches textures lazily.
 * Stretch lo/hi are read from gui.lo / gui.hi at call time and
 * trigger a cache rebuild when they change. */
void flis_gpu_compose_render(GtkSnapshot *snapshot,
                              GSList *layers,
                              guint canvas_w, guint canvas_h,
                              const graphene_rect_t *dst_rect,
                              const graphene_rect_t *visible_canvas,
                              GskScalingFilter filter);

/* Invalidate the cached textures belonging to one layer.  Call after
 * a layer's pixel data, tint, or lmask data changes.  No-op if the
 * layer was never cached. */
void flis_gpu_compose_invalidate_layer(gint item_id);

/* Drop just the lmask texture for one layer.  Cheaper than
 * invalidate_layer when only the layer mask data changed.  No-op
 * if the layer was never cached or had no lmask. */
void flis_gpu_compose_invalidate_lmask(gint item_id);

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
