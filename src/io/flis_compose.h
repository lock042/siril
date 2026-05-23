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

/**
 * flis_compose_bake_tile_bgra8:
 * @lay:      layer to bake from (uses lay->fit, lay->has_tint, lay->layer_tint)
 * @tile_dim: tile side length; tile (tx, ty) covers layer-local display rows
 *            [ty*tile_dim, ty*tile_dim+tile_h) and cols
 *            [tx*tile_dim, tx*tile_dim+tile_w) — edge tiles may be smaller.
 * @tx:       tile column (0-based)
 * @ty:       tile row (0-based)
 * @tile_w:   actual tile width in pixels (≤ tile_dim)
 * @tile_h:   actual tile height in pixels (≤ tile_dim)
 * @mip:      output downsample factor; 1 = full resolution, 2 = half,
 *            4 = quarter.  Output texture dimensions become
 *            (tile_w / mip) × (tile_h / mip).  Stride-sampled (every
 *            mip-th source pixel) — fast but aliasing.  Must be a
 *            power of two and divide both tile_w and tile_h.
 * @guard:    output-pixel border included on every edge to feed GSK's
 *            sampler with neighbour data at tile boundaries (1 px is
 *            enough for bilinear; 2 px for trilinear seam-free mips).
 *            0 disables.  Source pixels @ guard*mip away from the
 *            canonical tile are sampled, edge-clamped at layer
 *            boundaries.  Output dimensions become
 *            (tile_w/mip + 2*guard) × (tile_h/mip + 2*guard).
 * @lut:      256K-entry uint8 lookup table mapping 16-bit pixel values
 *            to display bytes (the linked stretch LUT)
 * @out_bgra: caller-provided buffer of size
 *            (tile_w/mip + 2*guard) * (tile_h/mip + 2*guard) * 4 bytes
 *
 * Bakes one tile of @lay into BGRA8 with the supplied stretch LUT applied
 * and mono tint (if lay->has_tint) folded in.  Pixels are laid out
 * display-orientation (row 0 = layer's display top, i.e. FITS row
 * layer_h-1).  Returns TRUE on success, FALSE on bad inputs.
 *
 * Pure CPU; no GTK / GdkTexture dependency — this is the byte-baking
 * core that flis_gpu_compose.c wraps in gdk_memory_texture_new().  The
 * test harness in src/tests/test_flis_gpu_compose_bake.c exercises this
 * directly with synthetic fixtures.
 */
gboolean flis_compose_bake_tile_bgra8(const flis_layer_t *lay,
                                       int tile_dim,
                                       int tx, int ty,
                                       int tile_w, int tile_h,
                                       int mip,
                                       int guard,
                                       const BYTE *lut,
                                       uint8_t *out_bgra);

/**
 * flis_compose_bake_lmask_bgra8:
 * Bakes the layer's lmask data into a layer-sized BGRA8 buffer (grey
 * encoded as R=G=B=value, alpha=255).  @out_bgra must be sized
 * lay->lmask->w * lay->lmask->h * 4 bytes.  Returns FALSE if the layer
 * has no lmask.
 */
gboolean flis_compose_bake_lmask_bgra8(const flis_layer_t *lay,
                                        uint8_t *out_bgra);

#endif /* FLIS_COMPOSE_H */
