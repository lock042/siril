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
 * flis_gpu_compose.c — per-layer per-tile texture cache + GSK push_blend
 * compositing path (stages 3.2 + 3.3).
 *
 * Cache shape (post-§3.3 slice 3): a small fixed-size array of slots,
 * keyed by layer item_id.  Each slot holds:
 *   - tiles[] : a tile_cols × tile_rows array of GdkTexture *.  Each
 *               tile is up to FLIS_TILE_DIM × FLIS_TILE_DIM pixels (edge
 *               tiles may be smaller).  Pixels are layer-local in display
 *               orientation (row 0 = layer's display top) with stretch
 *               and (mono+tint) baked in.  Materialised lazily by the
 *               render path based on the visible-image rect — tiles
 *               outside the viewport are never built.
 *   - lmask_tex : single layer-sized BGRA8 GdkTexture, the layer mask
 *               uploaded as a luminance-mask source for
 *               GSK_MASK_MODE_LUMINANCE.  NULL if no lmask.  (Lmasks
 *               are not currently tiled — they're typically smaller
 *               than layer pixels and the mask node spans the layer.)
 *   - stretch_lo / stretch_hi / has_tint / tint_*: the LUT + bake
 *               parameters; mismatch on any → drop the entire tile
 *               grid for that layer (cheap re-build on demand).
 *
 * Snapshot composition: nested gtk_snapshot_push_blend nodes form a
 * left-leaning chain — the deepest "bottom" is the base layer, each
 * successive layer is wrapped by one push_blend whose top subtree
 * holds (and optionally masks / opacifies) the layer.  Within each
 * layer's "top subtree" we append all visible tiles for that layer.
 * Groups composite as a sub-tree of their member layers wrapped in
 * push_opacity(group_opacity) + push_blend(group_blend_mode).
 *
 * Limits (caller must check via flis_gpu_compose_compatible):
 *   - FLIS_BLEND_CHROMA has no GSK equivalent — fall back to CPU
 *   - Base layer must be canvas-sized at position (0,0) — non-sparse
 *     (FLIS spec invariant); upper layers may be sparse.
 *   - Base layer in a non-PASS_THROUGH group — that group would need
 *     to blend against empty content for its bottom, not yet expressed.
 */

#include "flis_gpu_compose.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"           /* roundf_to_WORD */
#include "core/siril_log.h"
#include "io/image_format_flis.h"
#include "gui-gtk4/gui_state.h"   /* gui.lo / gui.hi / gui.remap_index */

/* Tile side length, layer-local.  Sized to fit comfortably inside every
 * GPU's max-texture-size limit and to amortise the per-tile bookkeeping
 * cost — most small/medium FLIS layers fit in 1-4 tiles.  Matches
 * SIRIL_TILE_DIM in image_display.c for consistency. */
#define FLIS_TILE_DIM 2048

/* Capacity sized comfortably above typical FLIS layer counts (spec
 * recommendation: 2–8 layers).  Slots are O(N) scanned. */
#define FLIS_CACHE_MAX 64

struct cache_slot {
	gint        item_id;     /* 0 = empty slot */

	/* Layer pixel cache, tiled.  tiles is a tile_cols * tile_rows array
	 * of GdkTexture *; entries are NULL until materialised by the
	 * snapshot path (lazy).  All non-NULL entries are baked at the
	 * stretch + tint state recorded below and at the layer extent
	 * (layer_w, layer_h).  Any of those changing → drop the grid via
	 * slot_drop_tiles(). */
	int         tile_cols;
	int         tile_rows;
	GdkTexture **tiles;      /* tile_cols * tile_rows; NULL slots OK */
	guint64    *last_used;   /* parallel array; tile last access epoch */
	guint       layer_w;     /* layer extent the grid was sized for */
	guint       layer_h;

	/* Layer mask — single texture, layer-sized (not tiled). */
	GdkTexture *lmask_tex;
	guint       lmask_w;
	guint       lmask_h;

	/* Bake parameters for the tile contents.  Mismatch → invalidate
	 * whole grid. */
	WORD        stretch_lo;
	WORD        stretch_hi;
	gboolean    has_tint;
	double      tint_r, tint_g, tint_b;
};

static struct cache_slot g_cache[FLIS_CACHE_MAX];
static guint64           g_epoch;  /* monotonic counter for tile last_used */

static struct cache_slot *find_slot(gint item_id) {
	if (item_id <= 0) return NULL;
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		if (g_cache[i].item_id == item_id) return &g_cache[i];
	return NULL;
}

static struct cache_slot *find_or_alloc_slot(gint item_id) {
	struct cache_slot *s = find_slot(item_id);
	if (s) return s;
	for (int i = 0; i < FLIS_CACHE_MAX; i++) {
		if (g_cache[i].item_id == 0) {
			g_cache[i].item_id = item_id;
			return &g_cache[i];
		}
	}
	siril_log_warning(_("FLIS GPU cache full (>%d layers); falling back\n"),
	                  FLIS_CACHE_MAX);
	return NULL;
}

static void slot_drop_tiles(struct cache_slot *s) {
	if (s->tiles) {
		const int n = s->tile_cols * s->tile_rows;
		for (int i = 0; i < n; i++) {
			if (s->tiles[i]) g_object_unref(s->tiles[i]);
		}
		g_free(s->tiles);
		s->tiles = NULL;
	}
	if (s->last_used) {
		g_free(s->last_used);
		s->last_used = NULL;
	}
	s->tile_cols = s->tile_rows = 0;
	s->layer_w = s->layer_h = 0;
}

static void slot_drop_textures(struct cache_slot *s) {
	slot_drop_tiles(s);
	if (s->lmask_tex) { g_object_unref(s->lmask_tex); s->lmask_tex = NULL; }
	s->lmask_w = s->lmask_h = 0;
}

static void slot_clear(struct cache_slot *s) {
	slot_drop_textures(s);
	s->item_id = 0;
	s->stretch_lo = 0;
	s->stretch_hi = 0;
	s->has_tint = FALSE;
	s->tint_r = s->tint_g = s->tint_b = 0;
}

/* =====================================================================
 * Texture builders
 * ===================================================================== */

/* Materialise one tile of a layer at (tx, ty) into a GdkTexture.  The
 * tile occupies layer-local DISPLAY rows [ty*FLIS_TILE_DIM, +tile_h)
 * and cols [tx*FLIS_TILE_DIM, +tile_w) (edge tiles may be smaller).
 * Reads the layer's source fits (FITS bottom-up) for that region, applies
 * the linked stretch via gui.remap_index[0], bakes in tint for mono
 * layers with has_tint set, encodes BGRA8.  Returns NULL on failure. */
static GdkTexture *materialise_layer_tile(flis_layer_t *lay,
                                            int tx, int ty,
                                            int tile_w, int tile_h) {
	fits *fit = lay->fit;
	if (!fit) return NULL;

	const guint layer_w = (guint)fit->rx;
	const guint layer_h = (guint)fit->ry;
	const int   tile_x0 = tx * FLIS_TILE_DIM;
	const int   tile_y0 = ty * FLIS_TILE_DIM;  /* in layer-local display coords */
	const BYTE *lut = gui.remap_index[0];   /* linked stretch */
	const gboolean is_rgb = (fit->naxes[2] >= 3);

	const size_t n = (size_t)tile_w * tile_h;
	uint8_t *bgra = malloc(n * 4);
	if (!bgra) return NULL;

	const float tr = lay->has_tint ? (float)lay->layer_tint.r : 1.f;
	const float tg = lay->has_tint ? (float)lay->layer_tint.g : 1.f;
	const float tb = lay->has_tint ? (float)lay->layer_tint.b : 1.f;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int local_y = 0; local_y < tile_h; local_y++) {
		/* Display y within the layer: tile_y0 .. tile_y0+tile_h-1.
		 * FITS y within the layer (bottom-up): layer_h - 1 - disp_y. */
		const int layer_disp_y = tile_y0 + local_y;
		const int layer_fits_y = (int)layer_h - 1 - layer_disp_y;
		const size_t src_row = (size_t)layer_fits_y * layer_w + (size_t)tile_x0;
		const size_t dst_row = (size_t)local_y * tile_w * 4;
		for (int local_x = 0; local_x < tile_w; local_x++) {
			const size_t src_i = src_row + (size_t)local_x;
			const size_t dst_i = dst_row + (size_t)local_x * 4;
			BYTE r_b, g_b, b_b;
			if (is_rgb) {
				if (fit->type == DATA_USHORT) {
					r_b = lut[fit->pdata[0][src_i]];
					g_b = lut[fit->pdata[1][src_i]];
					b_b = lut[fit->pdata[2][src_i]];
				} else {
					WORD r = roundf_to_WORD(fit->fpdata[0][src_i] * USHRT_MAX_SINGLE);
					WORD g = roundf_to_WORD(fit->fpdata[1][src_i] * USHRT_MAX_SINGLE);
					WORD b = roundf_to_WORD(fit->fpdata[2][src_i] * USHRT_MAX_SINGLE);
					r_b = lut[r];
					g_b = lut[g];
					b_b = lut[b];
				}
			} else {
				BYTE v_b;
				if (fit->type == DATA_USHORT) {
					v_b = lut[fit->pdata[0][src_i]];
				} else {
					WORD v = roundf_to_WORD(fit->fpdata[0][src_i] * USHRT_MAX_SINGLE);
					v_b = lut[v];
				}
				if (lay->has_tint) {
					float vf = (float)v_b;
					r_b = (BYTE)fminf(255.f, vf * tr);
					g_b = (BYTE)fminf(255.f, vf * tg);
					b_b = (BYTE)fminf(255.f, vf * tb);
				} else {
					r_b = g_b = b_b = v_b;
				}
			}
			bgra[dst_i + 0] = b_b;
			bgra[dst_i + 1] = g_b;
			bgra[dst_i + 2] = r_b;
			bgra[dst_i + 3] = 255;
		}
	}

	GBytes *bytes = g_bytes_new_take(bgra, n * 4);
	GdkTexture *tex = gdk_memory_texture_new(tile_w, tile_h,
	                                          GDK_MEMORY_B8G8R8A8,
	                                          bytes, (gsize)tile_w * 4);
	g_bytes_unref(bytes);
	return tex;
}

/* Build a BGRA8 texture from a layer mask at the lmask's own dimensions,
 * encoded as grayscale (R=G=B=value, alpha=255).  Used as the mask
 * source for GSK_MASK_MODE_LUMINANCE.  Returns NULL if the layer has no
 * mask data or on failure.  The lmask matches the layer's extent (not
 * the canvas) for both full-canvas and sparse layers. */
static GdkTexture *build_lmask_texture(flis_layer_t *lay) {
	if (!lay->lmask || !lay->lmask->data) return NULL;

	const guint w = lay->lmask->w;
	const guint h = lay->lmask->h;
	const size_t n = (size_t)w * h;
	uint8_t *bgra = malloc(n * 4);
	if (!bgra) return NULL;
	const guint8 bp = lay->lmask->bitpix;
	const void *src = lay->lmask->data;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (guint y = 0; y < h; y++) {
		const size_t src_row = (size_t)(h - 1 - y) * w;
		const size_t dst_row = (size_t)y * w * 4;
		for (guint x = 0; x < w; x++) {
			const size_t src_i = src_row + x;
			const size_t dst_i = dst_row + (size_t)x * 4;
			BYTE v;
			switch (bp) {
				case 8:  v = ((const uint8_t *)src)[src_i]; break;
				case 16: v = (BYTE)(((const uint16_t *)src)[src_i] >> 8); break;
				case 32: {
					float f = ((const float *)src)[src_i];
					v = (BYTE)fminf(255.f, fmaxf(0.f, f * 255.f));
					break;
				}
				default: v = 0; break;
			}
			bgra[dst_i + 0] = v;
			bgra[dst_i + 1] = v;
			bgra[dst_i + 2] = v;
			bgra[dst_i + 3] = 255;
		}
	}

	GBytes *bytes = g_bytes_new_take(bgra, n * 4);
	GdkTexture *tex = gdk_memory_texture_new(w, h,
	                                          GDK_MEMORY_B8G8R8A8,
	                                          bytes, w * 4);
	g_bytes_unref(bytes);
	return tex;
}

/* Compute tile dimensions for tile (tx, ty) within a layer of (lw, lh).
 * Edge tiles on the right/bottom may be smaller than FLIS_TILE_DIM. */
static inline void tile_dim_for(guint lw, guint lh, int tx, int ty,
                                 int *out_w, int *out_h) {
	const int x0 = tx * FLIS_TILE_DIM;
	const int y0 = ty * FLIS_TILE_DIM;
	*out_w = MIN((int)lw - x0, FLIS_TILE_DIM);
	*out_h = MIN((int)lh - y0, FLIS_TILE_DIM);
}

/* Look up the slot for @lay and (re)allocate its tile grid + bake
 * parameters if any of stretch/tint/extent has changed.  Returns the
 * slot on success or NULL on alloc failure.  Does NOT materialise any
 * tile pixels — that's deferred to ensure_tile() below, called by the
 * render path for each visible tile. */
static struct cache_slot *ensure_layer_cache_ready(flis_layer_t *lay,
                                                    WORD lo, WORD hi) {
	struct cache_slot *s = find_or_alloc_slot(lay->item_id);
	if (!s) return NULL;
	if (!lay->fit) return NULL;

	const guint cur_w = (guint)lay->fit->rx;
	const guint cur_h = (guint)lay->fit->ry;

	gboolean stretch_changed = (s->stretch_lo != lo || s->stretch_hi != hi);
	gboolean tint_changed = (s->has_tint != lay->has_tint) ||
	                        (lay->has_tint && (
	                            s->tint_r != lay->layer_tint.r ||
	                            s->tint_g != lay->layer_tint.g ||
	                            s->tint_b != lay->layer_tint.b));
	gboolean extent_changed = (s->layer_w != cur_w || s->layer_h != cur_h);

	if (!s->tiles || stretch_changed || tint_changed || extent_changed) {
		slot_drop_tiles(s);
		s->layer_w = cur_w;
		s->layer_h = cur_h;
		s->tile_cols = (int)((cur_w + FLIS_TILE_DIM - 1) / FLIS_TILE_DIM);
		s->tile_rows = (int)((cur_h + FLIS_TILE_DIM - 1) / FLIS_TILE_DIM);
		if (s->tile_cols <= 0 || s->tile_rows <= 0) {
			slot_clear(s);
			return NULL;
		}
		const int n = s->tile_cols * s->tile_rows;
		s->tiles     = g_new0(GdkTexture *, n);
		s->last_used = g_new0(guint64, n);
		s->stretch_lo = lo;
		s->stretch_hi = hi;
		s->has_tint = lay->has_tint;
		s->tint_r = lay->layer_tint.r;
		s->tint_g = lay->layer_tint.g;
		s->tint_b = lay->layer_tint.b;
	}

	/* Lmask texture — single layer-sized BGRA, rebuilt only when the
	 * mask data or extent changes (caller invalidates externally). */
	gboolean want_lmask = (lay->lmask && lay->lmask_active);
	if (want_lmask) {
		const guint mw = lay->lmask->w;
		const guint mh = lay->lmask->h;
		gboolean lm_extent_changed = (s->lmask_w != mw || s->lmask_h != mh);
		if (!s->lmask_tex || lm_extent_changed) {
			if (s->lmask_tex) { g_object_unref(s->lmask_tex); s->lmask_tex = NULL; }
			s->lmask_tex = build_lmask_texture(lay);
			if (s->lmask_tex) {
				s->lmask_w = mw;
				s->lmask_h = mh;
			}
		}
	} else if (s->lmask_tex) {
		g_object_unref(s->lmask_tex);
		s->lmask_tex = NULL;
		s->lmask_w = s->lmask_h = 0;
	}

	return s;
}

/* Ensure tile (tx, ty) of a layer is materialised in @s.  Bumps the
 * tile's last_used epoch.  Returns the texture (borrowed; do not unref)
 * or NULL on failure. */
static GdkTexture *ensure_tile(struct cache_slot *s, flis_layer_t *lay,
                                int tx, int ty) {
	if (tx < 0 || ty < 0 || tx >= s->tile_cols || ty >= s->tile_rows)
		return NULL;
	const int idx = ty * s->tile_cols + tx;
	if (!s->tiles[idx]) {
		int tw, th;
		tile_dim_for(s->layer_w, s->layer_h, tx, ty, &tw, &th);
		if (tw <= 0 || th <= 0) return NULL;
		s->tiles[idx] = materialise_layer_tile(lay, tx, ty, tw, th);
	}
	if (s->tiles[idx])
		s->last_used[idx] = ++g_epoch;
	return s->tiles[idx];
}

/* =====================================================================
 * GSK blend-mode translation
 * ===================================================================== */

static gboolean translate_blend_mode(flis_blend_mode_t mode, GskBlendMode *out) {
	switch (mode) {
		case FLIS_BLEND_NORMAL:       *out = GSK_BLEND_MODE_DEFAULT;     return TRUE;
		case FLIS_BLEND_MULTIPLY:     *out = GSK_BLEND_MODE_MULTIPLY;    return TRUE;
		case FLIS_BLEND_SCREEN:       *out = GSK_BLEND_MODE_SCREEN;      return TRUE;
		case FLIS_BLEND_OVERLAY:      *out = GSK_BLEND_MODE_OVERLAY;     return TRUE;
		case FLIS_BLEND_DARKEN:       *out = GSK_BLEND_MODE_DARKEN;      return TRUE;
		case FLIS_BLEND_LIGHTEN:      *out = GSK_BLEND_MODE_LIGHTEN;     return TRUE;
		case FLIS_BLEND_COLOR_DODGE:  *out = GSK_BLEND_MODE_COLOR_DODGE; return TRUE;
		case FLIS_BLEND_COLOR_BURN:   *out = GSK_BLEND_MODE_COLOR_BURN;  return TRUE;
		case FLIS_BLEND_HARD_LIGHT:   *out = GSK_BLEND_MODE_HARD_LIGHT;  return TRUE;
		case FLIS_BLEND_SOFT_LIGHT:   *out = GSK_BLEND_MODE_SOFT_LIGHT;  return TRUE;
		case FLIS_BLEND_DIFFERENCE:   *out = GSK_BLEND_MODE_DIFFERENCE;  return TRUE;
		case FLIS_BLEND_EXCLUSION:    *out = GSK_BLEND_MODE_EXCLUSION;   return TRUE;
		case FLIS_BLEND_HUE:          *out = GSK_BLEND_MODE_HUE;         return TRUE;
		case FLIS_BLEND_SATURATION:   *out = GSK_BLEND_MODE_SATURATION;  return TRUE;
		case FLIS_BLEND_COLOR:        *out = GSK_BLEND_MODE_COLOR;       return TRUE;
		case FLIS_BLEND_LUMINOSITY:   *out = GSK_BLEND_MODE_LUMINOSITY;  return TRUE;
		/* FLIS_BLEND_CHROMA (LRGB) and FLIS_BLEND_PASS_THROUGH (groups only)
		 * have no GSK equivalent — caller falls back. */
		default: return FALSE;
	}
}

/* =====================================================================
 * Compatibility check
 * ===================================================================== */

gboolean flis_gpu_compose_compatible(GSList *layers) {
	if (!layers) return FALSE;
	guint canvas_w = flis_canvas_rx();
	guint canvas_h = flis_canvas_ry();
	gboolean is_base = TRUE;
	for (GSList *l = layers; l; l = l->next, is_base = FALSE) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay || !lay->fit) return FALSE;
		if (is_base) {
			if (lay->position_x != 0 || lay->position_y != 0) return FALSE;
			if (lay->fit->rx != canvas_w || lay->fit->ry != canvas_h) return FALSE;
			if (lay->group_id != 0) {
				flis_group_t *grp = flis_group_get_by_id(lay->group_id);
				if (grp && grp->blend_mode != FLIS_BLEND_PASS_THROUGH)
					return FALSE;
			}
		}
		GskBlendMode dummy;
		if (!translate_blend_mode(lay->blend_mode, &dummy)) return FALSE;
		if (lay->group_id != 0) {
			flis_group_t *grp = flis_group_get_by_id(lay->group_id);
			if (grp && grp->blend_mode != FLIS_BLEND_PASS_THROUGH
			    && !translate_blend_mode(grp->blend_mode, &dummy))
				return FALSE;
		}
	}
	return TRUE;
}

/* =====================================================================
 * Render-item model
 *
 * A "render item" is a single composite atom that produces one source
 * subtree in the snapshot tree.  See the comment block above
 * flis_gpu_compose_render for the SINGLE / GROUP discriminator and
 * z-ordering rules.
 * ===================================================================== */

enum render_item_kind { RENDER_ITEM_SINGLE, RENDER_ITEM_GROUP };

struct render_item {
	enum render_item_kind kind;
	flis_blend_mode_t     blend_mode;
	gfloat                opacity;
	/* SINGLE */
	flis_layer_t         *layer;
	struct cache_slot    *slot;
	/* GROUP: flat array of member layers + their slots.  Group opacity
	 * lives in `opacity` above; each member's own opacity is applied
	 * by emit_group_member. */
	int                   member_count;
	flis_layer_t        **member_layers;
	struct cache_slot   **member_slots;
};

static void render_item_free(struct render_item *it) {
	g_free(it->member_layers);
	g_free(it->member_slots);
}

/* Compute the widget-space rect a layer's whole texture would occupy
 * (origin = layer's display position; size = layer extent × xx/yy).
 * See companion comment in slice 1 commit for the position_y
 * convention (top-down display offset, matches CPU compose oracle). */
static void layer_dst_rect(const flis_layer_t *lay,
                            guint canvas_w, guint canvas_h,
                            const graphene_rect_t *canvas_dst,
                            graphene_rect_t *out) {
	const float xx = canvas_dst->size.width  / (float)canvas_w;
	const float yy = canvas_dst->size.height / (float)canvas_h;
	const guint lw = (guint)lay->fit->rx;
	const guint lh = (guint)lay->fit->ry;
	const float dx = (float)lay->position_x * xx;
	const float dy = (float)lay->position_y * yy;
	*out = (graphene_rect_t)GRAPHENE_RECT_INIT(
		canvas_dst->origin.x + dx,
		canvas_dst->origin.y + dy,
		(float)lw * xx,
		(float)lh * yy);
}

/* Walk @lay's visible tiles (intersecting visible_canvas, in canvas
 * image-space) and append each as a scaled texture inside whatever
 * snapshot frame is currently open.  Materialises tiles on demand. */
static void emit_layer_tiles(GtkSnapshot *snap,
                              flis_layer_t *lay, struct cache_slot *slot,
                              guint canvas_w, guint canvas_h,
                              const graphene_rect_t *canvas_dst,
                              const graphene_rect_t *visible_canvas,
                              GskScalingFilter filter) {
	const float xx = canvas_dst->size.width  / (float)canvas_w;
	const float yy = canvas_dst->size.height / (float)canvas_h;
	const int   layer_disp_x0 = lay->position_x;
	const int   layer_disp_y0 = lay->position_y;

	for (int ty = 0; ty < slot->tile_rows; ty++) {
		for (int tx = 0; tx < slot->tile_cols; tx++) {
			int tw, th;
			tile_dim_for(slot->layer_w, slot->layer_h, tx, ty, &tw, &th);
			if (tw <= 0 || th <= 0) continue;

			/* Tile in canvas image-space (display top-down). */
			const float cx = (float)(layer_disp_x0 + tx * FLIS_TILE_DIM);
			const float cy = (float)(layer_disp_y0 + ty * FLIS_TILE_DIM);
			const float cw = (float)tw;
			const float ch = (float)th;

			/* Visibility cull: skip tiles entirely outside the
			 * caller's visible rect (when one was supplied — NULL
			 * means "draw everything"). */
			if (visible_canvas) {
				if (cx + cw < visible_canvas->origin.x) continue;
				if (cy + ch < visible_canvas->origin.y) continue;
				if (cx > visible_canvas->origin.x + visible_canvas->size.width)
					continue;
				if (cy > visible_canvas->origin.y + visible_canvas->size.height)
					continue;
			}

			GdkTexture *tex = ensure_tile(slot, lay, tx, ty);
			if (!tex) continue;

			const graphene_rect_t r = GRAPHENE_RECT_INIT(
				canvas_dst->origin.x + cx * xx,
				canvas_dst->origin.y + cy * yy,
				cw * xx,
				ch * yy);
			gtk_snapshot_append_scaled_texture(snap, tex, filter, &r);
		}
	}
}

/* Emit a SINGLE item: wrap any visible tiles of the layer in
 * push_opacity / push_mask as needed. */
static void emit_single(GtkSnapshot *snap, const struct render_item *it,
                         guint canvas_w, guint canvas_h,
                         const graphene_rect_t *canvas_dst,
                         const graphene_rect_t *visible_canvas,
                         GskScalingFilter filter) {
	gboolean wrap_opacity = (it->opacity < 0.999f);
	if (wrap_opacity)
		gtk_snapshot_push_opacity(snap, (double)it->opacity);
	gboolean wrap_mask = (it->slot->lmask_tex != NULL);
	if (wrap_mask)
		gtk_snapshot_push_mask(snap, GSK_MASK_MODE_LUMINANCE);

	emit_layer_tiles(snap, it->layer, it->slot,
	                 canvas_w, canvas_h, canvas_dst, visible_canvas, filter);

	if (wrap_mask) {
		gtk_snapshot_pop(snap);  /* close source subtree */
		/* Lmask spans the layer's full rect. */
		graphene_rect_t lr;
		layer_dst_rect(it->layer, canvas_w, canvas_h, canvas_dst, &lr);
		gtk_snapshot_append_scaled_texture(snap, it->slot->lmask_tex, filter, &lr);
		gtk_snapshot_pop(snap);  /* close mask node */
	}
	if (wrap_opacity)
		gtk_snapshot_pop(snap);
}

/* Append a member layer inside a GROUP sub-tree (its own tiles +
 * opacity/mask).  No group-opacity multiplication here — the group
 * opacity wraps the entire sub-tree at the caller. */
static void emit_group_member(GtkSnapshot *snap,
                               flis_layer_t *lay, struct cache_slot *slot,
                               guint canvas_w, guint canvas_h,
                               const graphene_rect_t *canvas_dst,
                               const graphene_rect_t *visible_canvas,
                               GskScalingFilter filter) {
	gboolean wrap_opacity = (lay->opacity < 0.999f);
	if (wrap_opacity)
		gtk_snapshot_push_opacity(snap, (double)lay->opacity);
	gboolean wrap_mask = (slot->lmask_tex != NULL);
	if (wrap_mask)
		gtk_snapshot_push_mask(snap, GSK_MASK_MODE_LUMINANCE);

	emit_layer_tiles(snap, lay, slot,
	                 canvas_w, canvas_h, canvas_dst, visible_canvas, filter);

	if (wrap_mask) {
		gtk_snapshot_pop(snap);
		graphene_rect_t lr;
		layer_dst_rect(lay, canvas_w, canvas_h, canvas_dst, &lr);
		gtk_snapshot_append_scaled_texture(snap, slot->lmask_tex, filter, &lr);
		gtk_snapshot_pop(snap);
	}
	if (wrap_opacity)
		gtk_snapshot_pop(snap);
}

/* Emit a GROUP item's sub-tree: push_opacity(group), then a nested
 * push_blend chain of the group's members (bottom-most member is the
 * group's base, no blend wrap for it). */
static void emit_group_subtree(GtkSnapshot *snap, const struct render_item *it,
                                guint canvas_w, guint canvas_h,
                                const graphene_rect_t *canvas_dst,
                                const graphene_rect_t *visible_canvas,
                                GskScalingFilter filter) {
	gboolean wrap_opacity = (it->opacity < 0.999f);
	if (wrap_opacity)
		gtk_snapshot_push_opacity(snap, (double)it->opacity);

	const int m = it->member_count;
	for (int i = 0; i < m - 1; i++) {
		flis_layer_t *top = it->member_layers[m - 1 - i];
		GskBlendMode mode;
		translate_blend_mode(top->blend_mode, &mode);
		gtk_snapshot_push_blend(snap, mode);
	}
	emit_group_member(snap, it->member_layers[0], it->member_slots[0],
	                  canvas_w, canvas_h, canvas_dst, visible_canvas, filter);
	for (int i = 1; i < m; i++) {
		gtk_snapshot_pop(snap);
		emit_group_member(snap, it->member_layers[i], it->member_slots[i],
		                  canvas_w, canvas_h, canvas_dst, visible_canvas, filter);
		gtk_snapshot_pop(snap);
	}

	if (wrap_opacity)
		gtk_snapshot_pop(snap);
}

/* =====================================================================
 * Main render entry
 * ===================================================================== */

void flis_gpu_compose_render(GtkSnapshot *snapshot,
                              GSList *layers,
                              guint canvas_w, guint canvas_h,
                              const graphene_rect_t *dst_rect,
                              const graphene_rect_t *visible_canvas,
                              GskScalingFilter filter) {
	if (!layers) return;

	const WORD lo = gui.lo;
	const WORD hi = gui.hi;

	/* ---- Phase 1: assemble render items in z-order. ---- */
	GArray *items = g_array_new(FALSE, TRUE, sizeof(struct render_item));
	GHashTable *consumed_groups = g_hash_table_new(g_direct_hash, g_direct_equal);

	for (GSList *l = layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay || !lay->fit || !lay->visible) continue;

		flis_group_t *grp = (lay->group_id != 0)
			? flis_group_get_by_id(lay->group_id) : NULL;
		if (grp && !grp->visible) continue;

		gboolean grouped_real = grp && grp->blend_mode != FLIS_BLEND_PASS_THROUGH;
		if (grouped_real) {
			if (g_hash_table_contains(consumed_groups, GINT_TO_POINTER(grp->item_id)))
				continue;
			GSList *members = flis_group_get_layers(grp);
			GArray *mem_layers = g_array_new(FALSE, FALSE, sizeof(flis_layer_t *));
			GArray *mem_slots  = g_array_new(FALSE, FALSE, sizeof(struct cache_slot *));
			gboolean ok = TRUE;
			for (GSList *m = members; m; m = m->next) {
				flis_layer_t *mlay = (flis_layer_t *)m->data;
				if (!mlay || !mlay->fit || !mlay->visible) continue;
				struct cache_slot *mslot = ensure_layer_cache_ready(mlay, lo, hi);
				if (!mslot) { ok = FALSE; break; }
				g_array_append_val(mem_layers, mlay);
				g_array_append_val(mem_slots,  mslot);
			}
			g_slist_free(members);
			if (!ok || mem_layers->len == 0) {
				g_array_free(mem_layers, TRUE);
				g_array_free(mem_slots, TRUE);
				if (!ok) goto bail;
				continue;
			}
			struct render_item it = { 0 };
			it.kind        = RENDER_ITEM_GROUP;
			it.blend_mode  = grp->blend_mode;
			it.opacity     = grp->opacity;
			it.member_count   = mem_layers->len;
			it.member_layers  = (flis_layer_t **)     g_array_free(mem_layers, FALSE);
			it.member_slots   = (struct cache_slot **)g_array_free(mem_slots,  FALSE);
			g_array_append_val(items, it);
			g_hash_table_insert(consumed_groups,
				GINT_TO_POINTER(grp->item_id), GINT_TO_POINTER(1));
		} else {
			struct cache_slot *slot = ensure_layer_cache_ready(lay, lo, hi);
			if (!slot) goto bail;
			struct render_item it = { 0 };
			it.kind        = RENDER_ITEM_SINGLE;
			it.blend_mode  = lay->blend_mode;
			it.opacity     = lay->opacity * (grp ? grp->opacity : 1.0f);
			it.layer       = lay;
			it.slot        = slot;
			g_array_append_val(items, it);
		}
	}

	if (items->len == 0) {
		g_array_free(items, TRUE);
		g_hash_table_destroy(consumed_groups);
		return;
	}

	/* ---- Phase 2: emit the nested push_blend chain over items. ---- */
	const int n = items->len;
	for (int i = 0; i < n - 1; i++) {
		struct render_item *top = &g_array_index(items, struct render_item, n - 1 - i);
		GskBlendMode mode;
		translate_blend_mode(top->blend_mode, &mode);
		gtk_snapshot_push_blend(snapshot, mode);
	}

	{
		struct render_item *base = &g_array_index(items, struct render_item, 0);
		if (base->kind == RENDER_ITEM_SINGLE)
			emit_single(snapshot, base, canvas_w, canvas_h,
			            dst_rect, visible_canvas, filter);
		else
			emit_group_subtree(snapshot, base, canvas_w, canvas_h,
			                   dst_rect, visible_canvas, filter);
	}
	for (int i = 1; i < n; i++) {
		gtk_snapshot_pop(snapshot);
		struct render_item *cur = &g_array_index(items, struct render_item, i);
		if (cur->kind == RENDER_ITEM_SINGLE)
			emit_single(snapshot, cur, canvas_w, canvas_h,
			            dst_rect, visible_canvas, filter);
		else
			emit_group_subtree(snapshot, cur, canvas_w, canvas_h,
			                   dst_rect, visible_canvas, filter);
		gtk_snapshot_pop(snapshot);
	}

	for (guint i = 0; i < items->len; i++)
		render_item_free(&g_array_index(items, struct render_item, i));
	g_array_free(items, TRUE);
	g_hash_table_destroy(consumed_groups);
	return;

bail:
	for (guint i = 0; i < items->len; i++)
		render_item_free(&g_array_index(items, struct render_item, i));
	g_array_free(items, TRUE);
	g_hash_table_destroy(consumed_groups);
}

/* =====================================================================
 * Invalidation API
 * ===================================================================== */

void flis_gpu_compose_invalidate_layer(gint item_id) {
	struct cache_slot *s = find_slot(item_id);
	if (s) slot_drop_textures(s);
}

void flis_gpu_compose_invalidate_all(void) {
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		slot_drop_textures(&g_cache[i]);
}

void flis_gpu_compose_free_all(void) {
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		slot_clear(&g_cache[i]);
}
