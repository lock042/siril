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
#include "io/flis_compose.h"       /* flis_compose_bake_tile_bgra8 / _lmask_bgra8 */
#include "gui-gtk4/gui_state.h"   /* gui.lo / gui.hi / gui.remap_index */

/* Tile side length, layer-local.  Sized to fit comfortably inside every
 * GPU's max-texture-size limit and to amortise the per-tile bookkeeping
 * cost — most small/medium FLIS layers fit in 1-4 tiles.  Matches
 * SIRIL_TILE_DIM in image_display.c for consistency. */
#define FLIS_TILE_DIM 2048

/* Capacity sized comfortably above typical FLIS layer counts (spec
 * recommendation: 2–8 layers).  Slots are O(N) scanned. */
#define FLIS_CACHE_MAX 64

/* Guard band, in SOURCE pixels (matches SIRIL_TILE_GUARD in
 * image_display.c).  GSK's TRILINEAR filter builds a mipmap pyramid
 * per texture; at mipmap level k each texel is averaged from 2^k
 * source pixels.  For adjacent tiles' boundary mip-texels to be
 * averaged from the same source-pixel range — i.e. for the boundary
 * to look seamless at zoom 1/2^k — each tile needs at least 2^k
 * source pixels of the next tile's data on each open edge.  64
 * source pixels covers zoom levels down to 1/64 of the bake mip,
 * which past mip=4 means past 1/256 effective image-space zoom —
 * well past anything users routinely use.  Without enough guard,
 * tile boundaries appear as visible seams (dark lines) at low zoom.
 *
 * The output-pixel guard scales inversely with the bake mip stride:
 * at mip=1 we need 64 output pixels; at mip=4 only 16 output pixels
 * give the same 64-source-pixel guard.  flis_tile_guard_out() does
 * that math. */
#define FLIS_TILE_GUARD 64

static inline int flis_tile_guard_out(int mip) {
	int g = FLIS_TILE_GUARD / mip;
	return g < 1 ? 1 : g;
}

/* Tile-bytes budget across the whole cache.  When materialising a new
 * tile would push g_cache_bytes above this, we evict the oldest
 * non-NULL tile (any layer, any position) until we're back under.
 * Sourced from com.pref.gui.flis_tile_budget_mb (default 256 MiB,
 * clamped 128–4096 MiB by settings.c).  The fallback constants here
 * cover the case where preferences haven't been initialised yet
 * (early startup / headless test binaries that don't load settings). */
#define FLIS_CACHE_BUDGET_DEFAULT_MB  256
#define FLIS_CACHE_BUDGET_MIN_MB      128
#define FLIS_CACHE_BUDGET_MAX_MB     4096

static inline gsize flis_cache_budget_bytes(void) {
	int mb = com.pref.gui.flis_tile_budget_mb;
	if (mb < FLIS_CACHE_BUDGET_MIN_MB || mb > FLIS_CACHE_BUDGET_MAX_MB)
		mb = FLIS_CACHE_BUDGET_DEFAULT_MB;
	return (gsize)mb * 1024 * 1024;
}

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
	gsize      *tile_bytes;  /* parallel array; per-tile byte size */
	guint       layer_w;     /* layer extent the grid was sized for */
	guint       layer_h;

	/* Layer mask — single texture, layer-sized (not tiled). */
	GdkTexture *lmask_tex;
	gsize       lmask_bytes;
	guint       lmask_w;
	guint       lmask_h;

	/* Bake parameters for the tile contents.  Mismatch → invalidate
	 * whole grid. */
	WORD        stretch_lo;
	WORD        stretch_hi;
	gboolean    has_tint;
	double      tint_r, tint_g, tint_b;
	/* Mip stride the tile grid was baked at (1, 2, 4 …).  Higher mip
	 * = smaller, downsampled textures (used at low zoom).  When the
	 * render path's required mip is finer (smaller integer) than the
	 * cached mip, the grid is dropped and rebuilt — coarser is OK to
	 * reuse because GSK downsamples on the GPU anyway. */
	int         mip_stride;
};

static struct cache_slot g_cache[FLIS_CACHE_MAX];
static guint64           g_epoch;       /* monotonic counter for tile last_used */
static gsize             g_cache_bytes; /* total tile + lmask bytes resident */

/* Background-prefetch state.  After each render we record what we just
 * drew; an idle handler then materialises tiles in a one-ring around
 * the visible rect so smooth panning doesn't bake on the critical
 * path.  Single source ID = at most one idle queued at a time. */
static guint  g_prefetch_idle_id;
static struct {
	gboolean valid;
	guint    canvas_w, canvas_h;
	graphene_rect_t visible;
	int      desired_mip;
} g_prefetch;

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

/* Drop a single tile in-place: unref the texture, deduct its bytes from
 * the global counter, and zero the slot entry.  Safe on NULL tiles. */
static void tile_release(struct cache_slot *s, int idx) {
	if (!s->tiles || !s->tiles[idx]) return;
	g_object_unref(s->tiles[idx]);
	s->tiles[idx] = NULL;
	if (s->tile_bytes) {
		g_cache_bytes -= s->tile_bytes[idx];
		s->tile_bytes[idx] = 0;
	}
	if (s->last_used) s->last_used[idx] = 0;
}

static void slot_drop_tiles(struct cache_slot *s) {
	if (s->tiles) {
		const int n = s->tile_cols * s->tile_rows;
		for (int i = 0; i < n; i++) tile_release(s, i);
		g_free(s->tiles);
		s->tiles = NULL;
	}
	g_free(s->last_used);  s->last_used  = NULL;
	g_free(s->tile_bytes); s->tile_bytes = NULL;
	s->tile_cols = s->tile_rows = 0;
	s->layer_w = s->layer_h = 0;
}

static void slot_drop_textures(struct cache_slot *s) {
	slot_drop_tiles(s);
	if (s->lmask_tex) {
		g_object_unref(s->lmask_tex);
		s->lmask_tex = NULL;
		g_cache_bytes -= s->lmask_bytes;
		s->lmask_bytes = 0;
	}
	s->lmask_w = s->lmask_h = 0;
}

/* LRU eviction: find the oldest non-NULL tile across all slots and
 * release it.  Returns the bytes freed, or 0 if nothing was evictable
 * (which happens when the cache is empty — meaning we're already
 * under budget and the caller's about-to-bake tile just can't fit
 * alone, an unrecoverable condition the caller treats as a soft fail). */
static gsize cache_evict_one(void) {
	struct cache_slot *victim_slot = NULL;
	int                victim_idx  = -1;
	guint64            oldest      = G_MAXUINT64;
	for (int s = 0; s < FLIS_CACHE_MAX; s++) {
		struct cache_slot *cs = &g_cache[s];
		if (!cs->tiles) continue;
		const int n = cs->tile_cols * cs->tile_rows;
		for (int i = 0; i < n; i++) {
			if (!cs->tiles[i]) continue;
			if (cs->last_used[i] < oldest) {
				oldest = cs->last_used[i];
				victim_slot = cs;
				victim_idx  = i;
			}
		}
	}
	if (!victim_slot || victim_idx < 0) return 0;
	const gsize bytes = victim_slot->tile_bytes[victim_idx];
	tile_release(victim_slot, victim_idx);
	return bytes;
}

/* Evict the oldest tiles until the global byte counter would
 * accommodate @incoming_bytes within the configured cache budget. */
static void cache_make_room_for(gsize incoming_bytes) {
	const gsize budget = flis_cache_budget_bytes();
	while (g_cache_bytes + incoming_bytes > budget) {
		if (cache_evict_one() == 0) break;  /* nothing left to evict */
	}
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

/* Materialise one tile of a layer at (tx, ty) into a GdkTexture, baked
 * at the supplied @mip stride (1 = full res, 2 = half-res, …).  The
 * source tile region is layer-local DISPLAY rows
 * [ty*FLIS_TILE_DIM, +tile_h) × cols [tx*FLIS_TILE_DIM, +tile_w);
 * output texture is (tile_w / mip) × (tile_h / mip) BGRA8 pixels.
 * Delegates byte bake to flis_compose_bake_tile_bgra8 (testable from
 * the headless harness). */
static GdkTexture *materialise_layer_tile(flis_layer_t *lay,
                                            int tx, int ty,
                                            int tile_w, int tile_h,
                                            int mip) {
	if (!lay || !lay->fit) return NULL;
	const int inner_w = tile_w / mip;
	const int inner_h = tile_h / mip;
	if (inner_w <= 0 || inner_h <= 0) return NULL;
	const int g_out = flis_tile_guard_out(mip);
	const int out_w = inner_w + 2 * g_out;
	const int out_h = inner_h + 2 * g_out;
	const size_t n = (size_t)out_w * out_h;
	uint8_t *bgra = malloc(n * 4);
	if (!bgra) return NULL;
	if (!flis_compose_bake_tile_bgra8(lay, FLIS_TILE_DIM, tx, ty,
	                                   tile_w, tile_h, mip, g_out,
	                                   gui.remap_index[0], bgra)) {
		free(bgra);
		return NULL;
	}
	GBytes *bytes = g_bytes_new_take(bgra, n * 4);
	GdkTexture *tex = gdk_memory_texture_new(out_w, out_h,
	                                          GDK_MEMORY_B8G8R8A8,
	                                          bytes, (gsize)out_w * 4);
	g_bytes_unref(bytes);
	return tex;
}

/* Build a BGRA8 texture from a layer's lmask at the lmask's own
 * dimensions.  Delegates the byte bake to flis_compose_bake_lmask_bgra8
 * in flis_compose.c and wraps in a GdkTexture. */
static GdkTexture *build_lmask_texture(flis_layer_t *lay) {
	if (!lay || !lay->lmask || !lay->lmask->data) return NULL;
	const guint w = lay->lmask->w;
	const guint h = lay->lmask->h;
	const size_t n = (size_t)w * h;
	uint8_t *bgra = malloc(n * 4);
	if (!bgra) return NULL;
	if (!flis_compose_bake_lmask_bgra8(lay, bgra)) {
		free(bgra);
		return NULL;
	}
	GBytes *bytes = g_bytes_new_take(bgra, n * 4);
	GdkTexture *tex = gdk_memory_texture_new(w, h,
	                                          GDK_MEMORY_B8G8R8A8,
	                                          bytes, w * 4);
	g_bytes_unref(bytes);
	return tex;
}

/* Compute tile dimensions for tile (tx, ty) within a layer of (lw, lh),
 * rounded down to a multiple of @mip so the bake function's
 * "tile_w % mip == 0" precondition holds.  Edge tiles on the right/
 * bottom may be smaller than FLIS_TILE_DIM and may also drop a few
 * pixels at the edge to satisfy mip alignment — those edge pixels
 * are not displayed when mipped, an acceptable trade for low-zoom
 * previews. */
static inline void tile_dim_for(guint lw, guint lh, int tx, int ty,
                                 int mip, int *out_w, int *out_h) {
	const int x0 = tx * FLIS_TILE_DIM;
	const int y0 = ty * FLIS_TILE_DIM;
	int tw = MIN((int)lw - x0, FLIS_TILE_DIM);
	int th = MIN((int)lh - y0, FLIS_TILE_DIM);
	if (mip > 1) {
		tw -= tw % mip;
		th -= th % mip;
	}
	*out_w = tw;
	*out_h = th;
}

/* Look up the slot for @lay and (re)allocate its tile grid + bake
 * parameters if any of stretch/tint/extent/mip has changed.  Returns the
 * slot on success or NULL on alloc failure.  Does NOT materialise any
 * tile pixels — that's deferred to ensure_tile() below, called by the
 * render path for each visible tile.
 *
 * @desired_mip is the smallest (finest) mip the render path currently
 * needs.  If the cached grid is at a coarser mip (mip_stride > desired)
 * we drop and rebuild so subsequent ensure_tile gets pixel-perfect
 * textures.  If cached is finer (mip_stride < desired) we keep it —
 * GSK downsamples on the GPU without extra CPU cost. */
static struct cache_slot *ensure_layer_cache_ready(flis_layer_t *lay,
                                                    WORD lo, WORD hi,
                                                    int desired_mip) {
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
	gboolean mip_too_coarse = (s->tiles && s->mip_stride > desired_mip);

	if (!s->tiles || stretch_changed || tint_changed || extent_changed
	    || mip_too_coarse) {
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
		s->tiles      = g_new0(GdkTexture *, n);
		s->last_used  = g_new0(guint64, n);
		s->tile_bytes = g_new0(gsize, n);
		s->stretch_lo = lo;
		s->stretch_hi = hi;
		s->has_tint = lay->has_tint;
		s->tint_r = lay->layer_tint.r;
		s->tint_g = lay->layer_tint.g;
		s->tint_b = lay->layer_tint.b;
		s->mip_stride = desired_mip;
	}

	/* Lmask texture — single layer-sized BGRA, rebuilt only when the
	 * mask data or extent changes (caller invalidates externally). */
	gboolean want_lmask = (lay->lmask && lay->lmask_active);
	if (want_lmask) {
		const guint mw = lay->lmask->w;
		const guint mh = lay->lmask->h;
		gboolean lm_extent_changed = (s->lmask_w != mw || s->lmask_h != mh);
		if (!s->lmask_tex || lm_extent_changed) {
			if (s->lmask_tex) {
				g_object_unref(s->lmask_tex);
				g_cache_bytes -= s->lmask_bytes;
				s->lmask_tex = NULL;
				s->lmask_bytes = 0;
			}
			const gsize lm_bytes = (gsize)mw * mh * 4;
			cache_make_room_for(lm_bytes);
			s->lmask_tex = build_lmask_texture(lay);
			if (s->lmask_tex) {
				s->lmask_w = mw;
				s->lmask_h = mh;
				s->lmask_bytes = lm_bytes;
				g_cache_bytes += lm_bytes;
			}
		}
	} else if (s->lmask_tex) {
		g_object_unref(s->lmask_tex);
		g_cache_bytes -= s->lmask_bytes;
		s->lmask_tex = NULL;
		s->lmask_bytes = 0;
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
		tile_dim_for(s->layer_w, s->layer_h, tx, ty, s->mip_stride, &tw, &th);
		if (tw <= 0 || th <= 0) return NULL;
		const int inner_w = tw / s->mip_stride;
		const int inner_h = th / s->mip_stride;
		const int g_out   = flis_tile_guard_out(s->mip_stride);
		const int out_w   = inner_w + 2 * g_out;
		const int out_h   = inner_h + 2 * g_out;
		const gsize bytes = (gsize)out_w * out_h * 4;
		cache_make_room_for(bytes);
		s->tiles[idx] = materialise_layer_tile(lay, tx, ty, tw, th, s->mip_stride);
		if (s->tiles[idx]) {
			s->tile_bytes[idx] = bytes;
			g_cache_bytes += bytes;
		}
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

	/* Push a single clip at the LAYER's extent so each tile can draw
	 * its full guard-padded texture without the per-tile clip — that
	 * approach left sub-pixel gaps at the canonical boundaries where
	 * neighbouring clips meet (visible as 1-px dark lines at all
	 * zooms).  Instead, adjacent tiles' guard regions overlap and
	 * paint identical pixels (the bake samples neighbouring tile
	 * source data into each tile's guard), so the union of all tiles
	 * inside the layer clip is seamless.  The layer-extent clip then
	 * trims the edge-clamped overflow at the layer's outer boundary,
	 * which would otherwise show as duplicated edge pixels beyond
	 * the layer's canonical area. */
	graphene_rect_t layer_clip;
	layer_dst_rect(lay, canvas_w, canvas_h, canvas_dst, &layer_clip);
	gtk_snapshot_push_clip(snap, &layer_clip);

	for (int ty = 0; ty < slot->tile_rows; ty++) {
		for (int tx = 0; tx < slot->tile_cols; tx++) {
			int tw, th;
			tile_dim_for(slot->layer_w, slot->layer_h, tx, ty,
			             slot->mip_stride, &tw, &th);
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

			/* The texture extends by FLIS_TILE_GUARD source pixels on
			 * each side beyond the canonical tile rect.  Draw at the
			 * expanded extent; the outer layer clip handles edge
			 * overflow, and adjacent tiles overlap in their guard
			 * regions with identical pixel data. */
			const float gpx = (float)FLIS_TILE_GUARD;
			const graphene_rect_t expanded = GRAPHENE_RECT_INIT(
				canvas_dst->origin.x + (cx - gpx) * xx,
				canvas_dst->origin.y + (cy - gpx) * yy,
				(cw + 2.f * gpx) * xx,
				(ch + 2.f * gpx) * yy);
			gtk_snapshot_append_scaled_texture(snap, tex, filter, &expanded);
		}
	}

	gtk_snapshot_pop(snap);   /* layer_clip */
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
 * Background prefetch
 *
 * After each render we know what was visible.  The prefetch idle then
 * walks the same layer stack and materialises tiles in a one-tile-wide
 * ring around the visible rect, so panning into adjacent area doesn't
 * stall the render thread on a cold bake.  Runs at G_PRIORITY_LOW so
 * it never preempts user input or animation, and bails immediately if
 * the cache is at budget (eviction churn would defeat the purpose).
 * ===================================================================== */

static gboolean prefetch_idle_cb(gpointer user_data);

static void schedule_prefetch_idle(guint canvas_w, guint canvas_h,
                                    const graphene_rect_t *visible,
                                    int desired_mip) {
	if (!visible) return;
	g_prefetch.valid       = TRUE;
	g_prefetch.canvas_w    = canvas_w;
	g_prefetch.canvas_h    = canvas_h;
	g_prefetch.visible     = *visible;
	g_prefetch.desired_mip = desired_mip;
	if (g_prefetch_idle_id == 0) {
		g_prefetch_idle_id = g_idle_add_full(G_PRIORITY_LOW,
		                                      prefetch_idle_cb, NULL, NULL);
	}
}

static gboolean prefetch_idle_cb(gpointer user_data) {
	(void)user_data;
	g_prefetch_idle_id = 0;
	if (!g_prefetch.valid) return G_SOURCE_REMOVE;
	g_prefetch.valid = FALSE;

	/* Don't prefetch when we're already near the byte budget — every
	 * new tile would evict a currently-visible one. */
	if (g_cache_bytes >= flis_cache_budget_bytes() * 3 / 4)
		return G_SOURCE_REMOVE;

	/* Read live layer list — the snapshot taken at render time isn't
	 * pointer-safe across the idle gap (a layer may have been freed). */
	if (!com.uniq || !com.uniq->layers) return G_SOURCE_REMOVE;

	/* Extend the visible rect by one tile-side worth of canvas-pixels
	 * (FLIS_TILE_DIM at mip=1).  This produces a one-tile ring around
	 * what's currently on screen. */
	const float ring = (float)FLIS_TILE_DIM;
	const graphene_rect_t extended = GRAPHENE_RECT_INIT(
		g_prefetch.visible.origin.x - ring,
		g_prefetch.visible.origin.y - ring,
		g_prefetch.visible.size.width  + 2.f * ring,
		g_prefetch.visible.size.height + 2.f * ring);

	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay || !lay->fit || !lay->visible) continue;

		struct cache_slot *s = find_slot(lay->item_id);
		if (!s || !s->tiles) continue;       /* not yet rendered */

		const int layer_disp_x0 = lay->position_x;
		const int layer_disp_y0 = lay->position_y;
		for (int ty = 0; ty < s->tile_rows; ty++) {
			for (int tx = 0; tx < s->tile_cols; tx++) {
				const int idx = ty * s->tile_cols + tx;
				if (s->tiles[idx]) continue;  /* already cached */

				int tw, th;
				tile_dim_for(s->layer_w, s->layer_h, tx, ty,
				             s->mip_stride, &tw, &th);
				if (tw <= 0 || th <= 0) continue;

				/* Tile rect in canvas image-space. */
				const float cx = (float)(layer_disp_x0 + tx * FLIS_TILE_DIM);
				const float cy = (float)(layer_disp_y0 + ty * FLIS_TILE_DIM);
				const float cw = (float)tw;
				const float ch = (float)th;

				/* Skip if outside extended rect (we only want the ring). */
				if (cx + cw < extended.origin.x) continue;
				if (cy + ch < extended.origin.y) continue;
				if (cx > extended.origin.x + extended.size.width) continue;
				if (cy > extended.origin.y + extended.size.height) continue;
				/* Skip if inside the strictly-visible rect — already
				 * baked by the render call. */
				if (cx >= g_prefetch.visible.origin.x &&
				    cy >= g_prefetch.visible.origin.y &&
				    cx + cw <= g_prefetch.visible.origin.x + g_prefetch.visible.size.width &&
				    cy + ch <= g_prefetch.visible.origin.y + g_prefetch.visible.size.height)
					continue;

				ensure_tile(s, lay, tx, ty);

				/* Bail if a single ring tile pushed us up to budget;
				 * the next render will reschedule us. */
				if (g_cache_bytes >= flis_cache_budget_bytes() * 3 / 4)
					return G_SOURCE_REMOVE;
			}
		}
	}
	return G_SOURCE_REMOVE;
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

	/* Pick the smallest mip the GPU sampler still benefits from at the
	 * current zoom.  zoom = widget pixels per canvas pixel; below 0.5
	 * we're sampling fewer than half the source pixels per output, so
	 * baking at half-res is lossless from the screen's POV.  Capped at
	 * 8× so we don't degrade quality beyond mip3 even at extreme
	 * fit-to-window.  Powers-of-two only — the bake function asserts. */
	int desired_mip = 1;
	if (dst_rect && canvas_w > 0 && canvas_h > 0) {
		const float zoom_x = dst_rect->size.width  / (float)canvas_w;
		const float zoom_y = dst_rect->size.height / (float)canvas_h;
		const float zoom = MIN(zoom_x, zoom_y);
		if      (zoom < 0.125f) desired_mip = 8;
		else if (zoom < 0.25f)  desired_mip = 4;
		else if (zoom < 0.5f)   desired_mip = 2;
	}

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
				struct cache_slot *mslot = ensure_layer_cache_ready(mlay, lo, hi, desired_mip);
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
			struct cache_slot *slot = ensure_layer_cache_ready(lay, lo, hi, desired_mip);
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

	/* Queue a background prefetch of the surrounding tile ring so the
	 * next pan doesn't stall on a cold bake.  Skipped if the caller
	 * passed no visibility hint (visible_canvas == NULL means "draw
	 * everything", so there's nothing to prefetch). */
	if (visible_canvas)
		schedule_prefetch_idle(canvas_w, canvas_h, visible_canvas, desired_mip);
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

void flis_gpu_compose_invalidate_lmask(gint item_id) {
	struct cache_slot *s = find_slot(item_id);
	if (s && s->lmask_tex) {
		g_object_unref(s->lmask_tex);
		g_cache_bytes -= s->lmask_bytes;
		s->lmask_tex = NULL;
		s->lmask_bytes = 0;
		s->lmask_w = s->lmask_h = 0;
	}
}

void flis_gpu_compose_invalidate_all(void) {
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		slot_drop_textures(&g_cache[i]);
}

void flis_gpu_compose_free_all(void) {
	if (g_prefetch_idle_id != 0) {
		g_source_remove(g_prefetch_idle_id);
		g_prefetch_idle_id = 0;
	}
	g_prefetch.valid = FALSE;
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		slot_clear(&g_cache[i]);
}
