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
 * flis_gpu_compose.c — per-layer texture cache + GSK push_blend
 * compositing path (stage 3.2).
 *
 * Cache shape: a small fixed-size array of slots, keyed by layer
 * item_id.  Each slot holds:
 *   - layer_tex  : canvas-sized BGRA8 GdkTexture, the layer's pixels
 *                  with stretch and (mono+tint) baked in
 *   - lmask_tex  : canvas-sized BGRA8 GdkTexture, the layer mask
 *                  uploaded as a luminance-mask source for
 *                  GSK_MASK_MODE_LUMINANCE.  NULL if no lmask.
 *   - stretch_lo / stretch_hi : the gui.lo / gui.hi the layer
 *                  texture was baked at; mismatch → rebuild
 *
 * Snapshot composition: nested gtk_snapshot_push_blend nodes form a
 * left-leaning chain — the deepest "bottom" is the base layer, each
 * successive layer is wrapped by one push_blend whose top subtree
 * holds (and optionally masks / opacifies) the layer.  The order is
 * preserved by the GSL list being already sorted ascending by
 * layer_order (load_flis / flis_layer_add maintain that invariant).
 *
 * Limits (caller must check via flis_gpu_compose_compatible):
 *   - FLIS_BLEND_CHROMA has no GSK equivalent — fall back to CPU
 *   - Sparse layers (position != 0 OR extent != canvas) — fall back;
 *     §3.3 will replace canvas-sized textures with per-tile ones and
 *     handle sparse correctly via per-tile placement.
 *   - Group membership — also §3.3 territory.
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

/* Capacity sized comfortably above typical FLIS layer counts (spec
 * recommendation: 2–8 layers).  Slots are O(N) scanned. */
#define FLIS_CACHE_MAX 64

struct cache_slot {
	gint        item_id;     /* 0 = empty slot */
	GdkTexture *layer_tex;
	GdkTexture *lmask_tex;
	WORD        stretch_lo;
	WORD        stretch_hi;
	gboolean    has_tint;    /* baked-in tint state — invalidate if toggled */
	double      tint_r, tint_g, tint_b;
	guint       layer_w;     /* extent the texture was built for; rebuild on mismatch */
	guint       layer_h;
	guint       lmask_w;     /* same for the mask texture */
	guint       lmask_h;
};

static struct cache_slot g_cache[FLIS_CACHE_MAX];

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

static void slot_drop_textures(struct cache_slot *s) {
	if (s->layer_tex) { g_object_unref(s->layer_tex); s->layer_tex = NULL; }
	if (s->lmask_tex) { g_object_unref(s->lmask_tex); s->lmask_tex = NULL; }
}

static void slot_clear(struct cache_slot *s) {
	slot_drop_textures(s);
	s->item_id = 0;
	s->stretch_lo = 0;
	s->stretch_hi = 0;
	s->has_tint = FALSE;
	s->tint_r = s->tint_g = s->tint_b = 0;
	s->layer_w = s->layer_h = 0;
	s->lmask_w = s->lmask_h = 0;
}

/* =====================================================================
 * Texture builders
 * ===================================================================== */

/* Build a layer-sized BGRA8 GdkTexture from a layer's pixel data.
 * Applies the linked stretch via gui.remap_index[0]; bakes in tint
 * for mono layers with has_tint set.  Texture dimensions equal the
 * layer's fit extent (which may be smaller than the canvas for sparse
 * layers — those land in the snapshot tree at position_x/y).  Returns
 * NULL on failure. */
static GdkTexture *build_layer_texture(flis_layer_t *lay) {
	fits *fit = lay->fit;
	if (!fit) return NULL;

	const guint w = (guint)fit->rx;
	const guint h = (guint)fit->ry;
	const BYTE *lut = gui.remap_index[0];   /* linked stretch */
	const gboolean is_rgb = (fit->naxes[2] >= 3);
	const size_t n = (size_t)w * h;
	uint8_t *bgra = malloc(n * 4);
	if (!bgra) return NULL;

	const float tr = lay->has_tint ? (float)lay->layer_tint.r : 1.f;
	const float tg = lay->has_tint ? (float)lay->layer_tint.g : 1.f;
	const float tb = lay->has_tint ? (float)lay->layer_tint.b : 1.f;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (guint y = 0; y < h; y++) {
		/* FITS convention: row 0 = bottom; display y=0 = top.
		 * Source row = h - 1 - y. */
		const size_t src_row = (size_t)(h - 1 - y) * w;
		const size_t dst_row = (size_t)y * w * 4;
		for (guint x = 0; x < w; x++) {
			const size_t src_i = src_row + x;
			const size_t dst_i = dst_row + (size_t)x * 4;
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
				/* Mono: broadcast (with tint if applicable).
				 * Tint is applied to the stretched byte to preserve
				 * the colour ratio at all display intensities. */
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
	GdkTexture *tex = gdk_memory_texture_new(w, h,
	                                          GDK_MEMORY_B8G8R8A8,
	                                          bytes, w * 4);
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

/* Look up or rebuild a layer's cache slot for the given stretch.
 * Returns the slot on success (with layer_tex non-NULL) or NULL on
 * any failure (caller should fall back). */
static struct cache_slot *ensure_layer_built(flis_layer_t *lay,
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

	if (!s->layer_tex || stretch_changed || tint_changed || extent_changed) {
		if (s->layer_tex) { g_object_unref(s->layer_tex); s->layer_tex = NULL; }
		s->layer_tex = build_layer_texture(lay);
		if (!s->layer_tex) {
			slot_clear(s);
			return NULL;
		}
		s->stretch_lo = lo;
		s->stretch_hi = hi;
		s->has_tint = lay->has_tint;
		s->tint_r = lay->layer_tint.r;
		s->tint_g = lay->layer_tint.g;
		s->tint_b = lay->layer_tint.b;
		s->layer_w = cur_w;
		s->layer_h = cur_h;
	}

	/* Lmask texture is stretch-independent and only rebuilt when the
	 * mask data itself changes (signalled by external invalidation) or
	 * the layer/lmask extent changes. */
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
			/* lmask build failure is non-fatal — just skip masking */
		}
	} else if (s->lmask_tex) {
		g_object_unref(s->lmask_tex);
		s->lmask_tex = NULL;
		s->lmask_w = s->lmask_h = 0;
	}

	return s;
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
		/* Group membership not yet supported via GPU path. */
		if (lay->group_id != 0) return FALSE;
		/* Base layer must fill the canvas — the FLIS spec invariant.
		 * Upper layers may be sparse (positioned anywhere with their
		 * own extent); the snapshot path places them via per-layer
		 * dst_rects.  Sparse-layer support added in §3.3 slice 1. */
		if (is_base) {
			if (lay->position_x != 0 || lay->position_y != 0) return FALSE;
			if (lay->fit->rx != canvas_w || lay->fit->ry != canvas_h) return FALSE;
		}
		/* Blend mode must translate to a GSK mode. */
		GskBlendMode dummy;
		if (!translate_blend_mode(lay->blend_mode, &dummy)) return FALSE;
	}
	return TRUE;
}

/* =====================================================================
 * Per-layer append (texture + optional mask + optional opacity)
 * ===================================================================== */

static void append_layer_full(GtkSnapshot *snap, struct cache_slot *s,
                               flis_layer_t *lay,
                               const graphene_rect_t *dst,
                               GskScalingFilter filter) {
	gboolean wrap_opacity = (lay->opacity < 0.999f);
	gboolean wrap_mask    = (s->lmask_tex != NULL);

	if (wrap_opacity)
		gtk_snapshot_push_opacity(snap, (double)lay->opacity);
	if (wrap_mask)
		gtk_snapshot_push_mask(snap, GSK_MASK_MODE_LUMINANCE);

	/* Source subtree (will be masked if mask was pushed) */
	gtk_snapshot_append_scaled_texture(snap, s->layer_tex, filter, dst);

	if (wrap_mask) {
		gtk_snapshot_pop(snap);  /* close source */
		gtk_snapshot_append_scaled_texture(snap, s->lmask_tex, filter, dst);
		gtk_snapshot_pop(snap);  /* close mask node */
	}
	if (wrap_opacity)
		gtk_snapshot_pop(snap);  /* close opacity */
}

/* =====================================================================
 * Main render entry
 * ===================================================================== */

/* Compute the widget-space rect a layer's texture should be painted
 * into, given the canvas dst_rect.  For full-canvas layers this is
 * just dst_rect; for sparse layers it's translated to the layer's
 * position.
 *
 * Coordinate-system note: the FLIS spec describes position_y as FITS
 * (bottom-up) but the CPU compose kernel (flis_compose.c) writes the
 * layer's pixels at FITS canvas rows [H-position_y-lH, H-position_y),
 * which through the standard fpdata→display y-flip ends up displaying
 * the layer at top-down rows [position_y, position_y+lH).  In other
 * words, position_y behaves as a top-down display offset of the layer's
 * top edge from the canvas top.  Match that convention so the GPU
 * path puts the layer in the same place the CPU oracle does. */
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

void flis_gpu_compose_render(GtkSnapshot *snapshot,
                              GSList *layers,
                              guint canvas_w, guint canvas_h,
                              const graphene_rect_t *dst_rect,
                              GskScalingFilter filter) {
	if (!layers) return;

	/* Filter to visible layers only — invisible ones produce no
	 * contribution but would still push a needless blend frame. */
	GSList *visible = NULL;
	for (GSList *l = layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (lay && lay->visible) visible = g_slist_append(visible, lay);
	}
	if (!visible) return;

	/* Build all layer textures up front so we can short-circuit on
	 * any failure before pushing any snapshot nodes. */
	const WORD lo = gui.lo;
	const WORD hi = gui.hi;
	const int n = g_slist_length(visible);
	struct cache_slot **slots = g_new0(struct cache_slot *, n);
	graphene_rect_t *rects = g_new0(graphene_rect_t, n);
	GSList *it = visible;
	for (int i = 0; i < n; i++, it = it->next) {
		flis_layer_t *lay = (flis_layer_t *)it->data;
		slots[i] = ensure_layer_built(lay, lo, hi);
		if (!slots[i] || !slots[i]->layer_tex) {
			/* Build failure: bail and let caller fall back to CPU. */
			g_free(slots);
			g_free(rects);
			g_slist_free(visible);
			return;
		}
		layer_dst_rect(lay, canvas_w, canvas_h, dst_rect, &rects[i]);
	}

	/* Nested push_blend chain (see header comment).  For N visible
	 * layers we push N-1 blends, each pop emits one layer in turn. */
	GSList *iter_layers = visible;
	for (int i = 0; i < n - 1; i++) {
		/* Push blends in outer-to-inner order — outermost = top layer */
		flis_layer_t *top = (flis_layer_t *)g_slist_nth_data(visible, n - 1 - i);
		GskBlendMode mode;
		translate_blend_mode(top->blend_mode, &mode);
		gtk_snapshot_push_blend(snapshot, mode);
	}

	/* Append the base layer (deepest "bottom"). */
	flis_layer_t *base = (flis_layer_t *)iter_layers->data;
	append_layer_full(snapshot, slots[0], base, &rects[0], filter);

	/* For each subsequent layer, pop into its top subtree, emit it,
	 * pop again to finalise the blend node. */
	for (int i = 1; i < n; i++) {
		gtk_snapshot_pop(snapshot);  /* switch to top subtree of this blend */
		flis_layer_t *lay = (flis_layer_t *)g_slist_nth_data(visible, i);
		append_layer_full(snapshot, slots[i], lay, &rects[i], filter);
		gtk_snapshot_pop(snapshot);  /* finalise this blend node */
	}

	g_free(slots);
	g_free(rects);
	g_slist_free(visible);
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
	/* Preserve item_id and tint state so the next render can detect
	 * whether anything actually changed. */
}

void flis_gpu_compose_free_all(void) {
	for (int i = 0; i < FLIS_CACHE_MAX; i++)
		slot_clear(&g_cache[i]);
}
