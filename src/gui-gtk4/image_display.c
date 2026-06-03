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

#include <math.h>
#include <ctype.h>

#include "git-version.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "algos/astrometry_solver.h"
#include "algos/colors.h"
#include "algos/ccd-inspector.h"
#include "gui-gtk4/ccd-inspector.h"
#include "algos/background_extraction.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/sorting.h"
#include "io/annotation_catalogues.h"
#include "filters/mtf.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/user_polygons.h"
#include "gui-gtk4/histo_display.h"
#include "livestacking/livestacking.h"
#include "histogram.h"
#include "registration/matching/degtorad.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

#include <wcslib.h>
#include <wcsfix.h>

#include "image_display.h"

/* is gfit->icc_profile identical to the monitor profile, if so we can avoid the
 * transform */
static cmsBool identical = FALSE;

#define ANGLE_TOP 315. * DEGTORAD
#define ANGLE_BOT 45. * DEGTORAD

/* remap index data, an index for each layer */
static float last_pente;
static display_mode last_mode;

/* STF (auto-stretch) data */
static gboolean stf_computed = FALSE; // Flag to know if STF parameters are available
static struct mtf_params stf[3];

/* widgets for draw_reg_data*/
GtkDropDown *seqcombo;
GtkToggleButton *drawframe;
static GtkWidget *rotation_dlg = NULL;

/* widgets for cut tool*/
static GtkWidget *cut_dialog = NULL, *cut_cdialog = NULL, *cut_sdialog = NULL;
static GtkCheckButton *tri_cut_toggle = NULL;
static GtkSpinButton *tri_cut_spin_step = NULL;

static GtkApplicationWindow *imgdisp_app_win = NULL;
static GtkCheckButton *imgdisp_autohd_item = NULL;
static GtkWidget *imgdisp_drawing_rgb = NULL;
static GtkWidget *imgdisp_drawing_r = NULL;

static void image_display_init_statics(void) {
	if (imgdisp_app_win) return;
	imgdisp_app_win = GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	imgdisp_autohd_item = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "autohd_item"));
	imgdisp_drawing_rgb = GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawingareargb"));
	imgdisp_drawing_r = GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawingarear"));
	seqcombo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "seqlist_dialog_combo"));
	drawframe = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "drawframe_check"));
	rotation_dlg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rotation_dialog"));
	cut_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_dialog"));
	cut_cdialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_coords_dialog"));
	cut_sdialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_spectroscopy_dialog"));
	tri_cut_toggle = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "cut_tri_cut"));
	tri_cut_spin_step = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "cut_tricut_step"));
}

static void invalidate_image_render_cache(int vport);

/* Tile side length.  Each viewport's buf is exposed to GSK as a grid of
 * GdkTextures of up to SIRIL_TILE_DIM × SIRIL_TILE_DIM pixels; right- and
 * bottom-edge tiles may be smaller.  The chosen size sits comfortably
 * inside every GPU's max-texture-size limit (typically ≥ 4096 on any
 * GL/Vulkan path) and divides easily into the typical megapixel-scale
 * images Siril works with. */
#define SIRIL_TILE_DIM 2048

/* Per-tile guard band (extra pixels along the right and bottom edges).
 * GSK's TRILINEAR filter builds an independent mipmap pyramid per
 * texture; at mipmap level k each texel is averaged from 2^k source
 * pixels.  For adjacent tiles' boundary mip-texels to be averaged from
 * the *same* source-pixel range (i.e. for the boundary to look seamless
 * at zoom 1/2^k), each tile needs a guard of at least 2^k pixels of the
 * next tile's data on each open edge.  64 covers zoom levels down to
 * 1/64, which is well past what users normally zoom out to.  Costs ~6 %
 * extra texture memory for interior tiles (2112² vs 2048²). */
#define SIRIL_TILE_GUARD 64

/* Per-viewport RAM budget — serves two roles:
 *  1. Eager-vs-lazy threshold in allocate_full_surface.  If the whole
 *     image fits inside this budget at native resolution we go eager
 *     (one contiguous buf, tiles slice into it, no LRU, no mip
 *     downsampling, no glitching on zoom).  This covers almost every
 *     realistic image — a 24 Mpix RGB frame is ~96 MB.
 *  2. Soft cap on lazy-mode resident tile bytes for the pathological
 *     outliers (huge mosaics) that don't fit eager.  materialise_tile
 *     evicts the LRU tile(s) until back under budget.
 * User-configurable via Preferences → Performance (com.pref.lazy_tile_cache_mb,
 * default 128 MB, up to 4 GB): a larger value loads more images eagerly (no
 * glitching) and lets the lazy cache hold the whole visible set of a big mosaic
 * so the LRU never has to blank a visible tile.  Read live so a settings change
 * takes effect on the next image (re)allocation.  Worst case across the four
 * vports (R, G, B, RGB) is 4× this; mono images use one view. */
#define SIRIL_LAZY_BUDGET    ((gsize)com.pref.lazy_tile_cache_mb * 1024 * 1024)

/* Cap on the downsample level a tile may be materialised at.  At mip M the
 * texture is 1/(2^M)² of the source area; at mip 6 a 2048-pixel tile yields
 * a 32-pixel texture (4 KB), and there's no useful saving past that point.
 * Also keeps tex_w >> mip strictly positive on every realistic tile. */
#define SIRIL_MAX_MIP        6

/* Long-edge size (pixels) of the lazy-mode full-image proxy.  One texture of
 * at most SIRIL_PROXY_MAX on its longer axis is kept resident per lazy
 * viewport and used as the fallback backdrop while real tiles are still being
 * materialised after a fast zoom-out.  2048 keeps it inside every GPU's
 * texture-size limit and makes it pixel-adequate at zoom_fit on any normal
 * monitor; it is point-sampled (not box-averaged) so its build cost is
 * O(proxy area) — a few megapixels — regardless of how large the source is. */
#define SIRIL_PROXY_MAX      2048

/* Choose the mip level we want to materialise tiles at, given the current
 * display zoom.  At mip M each texture pixel covers 2^M × 2^M source pixels;
 * the coarsest mip that still keeps the texture at least as dense as the
 * display is M_necessary = floor(-log2(zoom)).  We deliberately materialise
 * at one level finer (M_necessary − 1) so that zooming in by one octave
 * doesn't force a re-materialise of every visible tile (which would stutter
 * on every wheel notch).  Cost: 4× the memory of M_necessary, still 4^M_target
 * less than full-res.  At zoom ≥ 1 we stay at mip 0. */
static int compute_target_mip(double zoom) {
	if (zoom >= 1.0 || !(zoom > 0.0)) return 0;
	const int necessary = (int)floor(-log2(zoom));
	int target = necessary - 1;
	if (target < 0) target = 0;
	if (target > SIRIL_MAX_MIP) target = SIRIL_MAX_MIP;
	return target;
}

/* Return the dimensions of tile (tx, ty) in image-space.  Edge tiles may
 * be smaller than tile_dim. */
static inline void tile_dims(const struct image_view *view, int tx, int ty,
                              int *out_x0, int *out_y0, int *out_w, int *out_h) {
	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	const int x0 = tx * view->tile_dim;
	const int y0 = ty * view->tile_dim;
	*out_x0 = x0;
	*out_y0 = y0;
	*out_w = MIN(view->tile_dim, img_w - x0);
	*out_h = MIN(view->tile_dim, img_h - y0);
}

/* As tile_dims, but also returns the texture's actual width and height
 * including a SIRIL_TILE_GUARD-pixel "guard band" along the right and
 * bottom edges (when there's an adjacent tile to take the guard data
 * from).  The guard duplicates the first SIRIL_TILE_GUARD rows/columns
 * of the neighbouring tile so the scaling filter — including TRILINEAR
 * mipmap pre-filtering — samples consistent data across the boundary
 * instead of seeing per-tile pyramid discontinuities. */
static inline void tile_dims_padded(const struct image_view *view,
                                     int tx, int ty,
                                     int *out_x0, int *out_y0,
                                     int *out_w, int *out_h,
                                     int *out_tex_w, int *out_tex_h) {
	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	tile_dims(view, tx, ty, out_x0, out_y0, out_w, out_h);
	/* The guard is the smaller of SIRIL_TILE_GUARD and however many real
	 * image pixels actually exist past this tile's edge — never invent
	 * data, and never read past the image. */
	const int guard_x = (tx < view->tile_cols - 1)
		? MIN(SIRIL_TILE_GUARD, img_w - (*out_x0 + *out_w))
		: 0;
	const int guard_y = (ty < view->tile_rows - 1)
		? MIN(SIRIL_TILE_GUARD, img_h - (*out_y0 + *out_h))
		: 0;
	*out_tex_w = *out_w + guard_x;
	*out_tex_h = *out_h + guard_y;
}

/* Release a lazy-mode tile.  We unref the texture but do NOT free
 * tile->data directly: the texture's GBytes was constructed with
 * g_bytes_new_with_free_func(free), so the bytes are released
 * automatically when the GBytes' last reference drops — which is
 * deferred until GSK finishes any in-flight render that referenced
 * this texture, avoiding the use-after-free that would otherwise hit
 * gdk_memory_convert.
 *
 * lazy_bytes_used tracks the "logical" budget — we decrement it here
 * even though the physical bytes may live a little longer.  Worst case
 * the budget overshoots briefly during a busy redraw; for the
 * 512 MB-default target that's a sub-tile margin.  Must be called with
 * gui.cairo_mutex held. */
static void tile_release(struct image_view *view, struct image_tile *tile) {
	if (!view->lazy) return;  /* eager-mode tiles share view->buf */
	const gsize bytes = tile->data ? tile->bytes : 0;
	if (tile->texture) {
		g_object_unref(tile->texture);
		tile->texture = NULL;
	}
	tile->data = NULL;
	tile->bytes = 0;
	if (view->lazy_bytes_used >= bytes)
		view->lazy_bytes_used -= bytes;
	else
		view->lazy_bytes_used = 0;
	tile->dirty = TRUE;
}

/* Release every tile texture (and, in lazy mode, the bytes too).  Must be
 * called with gui.cairo_mutex held. */
static void view_drop_all_tile_textures(struct image_view *view) {
	if (!view->tiles) return;
	const int n = view->tile_cols * view->tile_rows;
	for (int i = 0; i < n; i++) {
		if (view->lazy) {
			tile_release(view, &view->tiles[i]);
		} else {
			if (view->tiles[i].texture) {
				g_object_unref(view->tiles[i].texture);
				view->tiles[i].texture = NULL;
			}
			/* eager tile->data points into view->buf; not owned. */
			view->tiles[i].data = NULL;
		}
	}
}

/* Eager mode: (re)build all tile textures, each wrapping a sub-region of
 * view->buf via the buf row stride.  Each tile's GBytes is a sub-bytes
 * of the parent view->buf_gbytes, so it retains a refcount on the parent
 * — the parent's free function (which frees view->buf) won't fire until
 * every tile texture that GSK still holds is also released, preventing
 * the use-after-free that happens if realloc/free runs while the
 * renderer is mid-frame.
 *
 * Non-rightmost / non-bottommost tiles include a SIRIL_TILE_GUARD-pixel
 * guard band beyond their nominal extent — the underlying buf already
 * contains the neighbouring tile's first pixels there, so the texture
 * simply extends into them via the buf row stride.  See
 * tile_dims_padded for the rationale (mipmap-consistent boundaries).
 *
 * Must be called with gui.cairo_mutex held. */
static void view_refresh_tile_textures_eager(struct image_view *view) {
	view_drop_all_tile_textures(view);
	if (!view->buf_gbytes || !view->tiles || view->buf_stride <= 0
	    || view->buf_height <= 0)
		return;
	for (int ty = 0; ty < view->tile_rows; ty++) {
		for (int tx = 0; tx < view->tile_cols; tx++) {
			int x0, y0, tw, th, tex_w, tex_h;
			tile_dims_padded(view, tx, ty, &x0, &y0, &tw, &th,
				&tex_w, &tex_h);
			if (tex_w <= 0 || tex_h <= 0) continue;
			const gsize offset = (gsize)y0 * view->buf_stride
			                   + (gsize)x0 * 4;
			const gsize len = (gsize)(tex_h - 1) * view->buf_stride
			                + (gsize)tex_w * 4;
			GBytes *bytes = g_bytes_new_from_bytes(view->buf_gbytes,
				offset, len);
			struct image_tile *t = &view->tiles[ty * view->tile_cols + tx];
			t->data = view->buf + offset;
			t->texture = gdk_memory_texture_new(tex_w, tex_h,
				GDK_MEMORY_B8G8R8X8, bytes, view->buf_stride);
			t->mip = 0;
			t->bytes = 0;   /* eager tiles don't count against lazy_bytes_used */
			t->dirty = FALSE;
			g_bytes_unref(bytes);
		}
	}
}

/* Lazy mode: mark every tile as needing re-materialisation on next snapshot.
 *
 * We deliberately KEEP each tile's existing texture (just flag it dirty)
 * instead of dropping it.  The stale texture — built from the previous LUT —
 * is a perfectly good thing to display for the one or few frames it takes the
 * background worker to rebuild it, and showing it avoids both a blank/proxy
 * "blur flash" and, more importantly, any synchronous main-thread re-fill.
 * That matters on big lazy images where a single full-resolution tile re-fill
 * (per-channel LUT, ×3 for RGB, plus ICC) is tens to hundreds of milliseconds:
 * a remap triggered by e.g. clearing a selection used to freeze the UI while
 * the visible tiles were rebuilt inline.  Now the snapshot shows the stale
 * textures immediately and the worker refreshes them off-thread.
 *
 * Must be called with gui.cairo_mutex held. */
static void view_invalidate_tiles_lazy(struct image_view *view) {
	if (!view->tiles) return;
	const int n = view->tile_cols * view->tile_rows;
	for (int i = 0; i < n; i++)
		view->tiles[i].dirty = TRUE;
	/* The LUT changed, so the resident proxy is stale too — rebuild it on the
	 * next snapshot that actually needs it. */
	view->proxy_dirty = TRUE;
	/* Advance the invalidate sequence so a worker job that captured the old LUT
	 * and is mid-fill discards its (now stale) texture instead of clearing the
	 * dirty flag we just set.  Without this, a remap landing during a tile fill
	 * leaves that tile clean-but-stale, and zoom-out never refreshes it because
	 * its mip stays "acceptable" — only a zoom-in forces the rebuild. */
	view->invalidate_seq++;
}

/* Compatibility wrapper called from every remap exit point.  In eager
 * mode this rebuilds all tile textures from the freshly-filled view->buf;
 * in lazy mode it just marks tiles dirty so the next snapshot
 * materialises the visible ones using the up-to-date LUT. */
static void view_refresh_tile_textures(struct image_view *view) {
	if (view->lazy)
		view_invalidate_tiles_lazy(view);
	else
		view_refresh_tile_textures_eager(view);
}


/* Allocate (or reallocate) the per-viewport tile grid for the current
 * gfit dimensions.  Images that fit in a single tile (≤ SIRIL_TILE_DIM
 * on both axes) get eager mode — one contiguous buf, no LRU; anything
 * larger uses lazy mode — tile bytes materialised on demand and bounded
 * by SIRIL_LAZY_BUDGET.  Returns 0 on success. */
static int allocate_full_surface(struct image_view *view) {
	g_mutex_lock(&gui.cairo_mutex);

	const int img_w = (int)gfit->rx;
	const int img_h = (int)gfit->ry;
	if (img_w <= 0 || img_h <= 0) {
		g_mutex_unlock(&gui.cairo_mutex);
		return 1;
	}

	const int stride   = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, img_w);
	const int tile_dim = SIRIL_TILE_DIM;
	const int tile_cols = (img_w + tile_dim - 1) / tile_dim;
	const int tile_rows = (img_h + tile_dim - 1) / tile_dim;
	const gsize total_bytes = (gsize)stride * img_h;
	/* Eager when the whole image fits in one viewport's budget at native
	 * resolution: a single contiguous buf, tile textures slice into it,
	 * no LRU and no mip-driven re-materialisation — so no glitches when
	 * zoom crosses a mip boundary.  Lazy is the fallback for pathological
	 * outliers (huge mosaics) where holding the whole image is infeasible. */
	const gboolean want_lazy = (total_bytes > SIRIL_LAZY_BUDGET);

	if (stride != view->buf_stride
				|| img_h != view->buf_height
				|| tile_cols != view->tile_cols
				|| tile_rows != view->tile_rows
				|| want_lazy != view->lazy
				|| !view->tiles
				|| (!want_lazy && !view->buf_gbytes)) {
		siril_log_debug("display buffer (re-)allocation %p — %d×%d, %d×%d tiles, %s mode\n",
				view, img_w, img_h, tile_cols, tile_rows,
				want_lazy ? "lazy" : "eager");

		/* Bump the generation so any materialise-worker job that's mid-flight
		 * for the old grid discards its result instead of writing into the
		 * tile array we're about to free/replace. */
		view->generation++;

		/* Drop existing textures.  In eager mode the tile textures hold
		 * sub-refs on view->buf_gbytes; once we unref our copy of
		 * buf_gbytes below, the bytes live until every still-referenced
		 * texture is released by GSK. */
		view_drop_all_tile_textures(view);
		g_free(view->tiles);
		view->tiles = NULL;

		if (view->buf_gbytes) {
			g_bytes_unref(view->buf_gbytes);
			view->buf_gbytes = NULL;
		}
		view->buf = NULL;

		/* The image dimensions / mode are changing — drop the proxy so it's
		 * rebuilt at the new size (and never lingers into eager mode). */
		if (view->proxy) {
			g_object_unref(view->proxy);
			view->proxy = NULL;
		}
		view->proxy_dirty = TRUE;
		view->vis_valid = FALSE;

		if (!want_lazy) {
			guchar *fresh = malloc(total_bytes * sizeof(guchar));
			if (!fresh) {
				PRINT_ALLOC_ERR;
				view->buf_stride = 0;
				view->buf_height = 0;
				view->tile_dim = view->tile_cols = view->tile_rows = 0;
				view->lazy = FALSE;
				view->lazy_bytes_used = 0;
				g_mutex_unlock(&gui.cairo_mutex);
				return 1;
			}
			view->buf_gbytes = g_bytes_new_with_free_func(
				fresh, total_bytes, free, fresh);
			view->buf = fresh;
		}

		view->buf_stride = stride;
		view->buf_height = img_h;
		view->tile_dim   = tile_dim;
		view->tile_cols  = tile_cols;
		view->tile_rows  = tile_rows;
		view->lazy       = want_lazy;
		view->lazy_bytes_used = 0;
		view->tiles = g_new0(struct image_tile, tile_cols * tile_rows);
		if (view->lazy)
			for (int i = 0; i < tile_cols * tile_rows; i++)
				view->tiles[i].dirty = TRUE;
	}
	g_mutex_unlock(&gui.cairo_mutex);
	return 0;
}

/* ── Lazy-mode materialisation ────────────────────────────────────────────
 *
 * In lazy mode the remap pass only recomputes the LUTs (HISTEQ histogram,
 * STF midtones, …) and marks every tile dirty.  The per-pixel work is
 * deferred to materialise_tile, which runs at snapshot time for each tile
 * the renderer is about to draw.  Per-tile bytes are LRU-evicted once the
 * viewport's resident set exceeds SIRIL_LAZY_BUDGET.
 *
 * All functions in this block must be called with gui.cairo_mutex held. */

static void materialise_tile_evict_until(struct image_view *view, gsize headroom);

/* Apply the LUT + flags to a tile rect of a gray (R/G/B) viewport and
 * build a fresh texture from a freshly-malloc'd byte buffer.  The
 * texture's GBytes owns the bytes via a free function, so when the
 * texture is unrefed (now, or later when GSK finishes rendering) the
 * bytes are reclaimed automatically.  Crucially we don't re-use the
 * tile's previous data buffer — GSK may still be reading it for an
 * in-flight frame, and writing into it in-place would corrupt that
 * upload (the assertion that fires inside gdk_memory_convert). */
static gboolean materialise_tile_gray(struct image_view *view, int vport,
                                       int tx, int ty, int mip,
                                       const BYTE *lut, gboolean hd_mode,
                                       gboolean inverted, gboolean cut_over,
                                       WORD remap_lo, WORD remap_hi) {
	(void) remap_hi; (void) cut_over; (void) remap_lo;
	int x0, y0, tw, th, tex_w, tex_h;
	tile_dims_padded(view, tx, ty, &x0, &y0, &tw, &th, &tex_w, &tex_h);
	if (tex_w <= 0 || tex_h <= 0) return TRUE;  /* empty edge tile */

	struct image_tile *tile = &view->tiles[ty * view->tile_cols + tx];
	/* At mip M each output pixel averages a step×step block of source pixels
	 * (with the LUT applied per source pixel — averaging in display-byte
	 * space, not in raw-pixel space, so the stretch curve is preserved).  The
	 * texture is handed to GSK at this reduced size; the snapshot rect stays
	 * in image-space coords so the existing display_matrix transform works
	 * unchanged and GSK upscales the small texture into the image-space rect. */
	const int step = 1 << mip;
	const int out_w = MAX(1, tex_w >> mip);
	const int out_h = MAX(1, tex_h >> mip);
	const gsize tile_bytes = (gsize)out_w * out_h * 4;

	/* Reclaim budget for the new buffer.  If the previous buffer for this
	 * tile is still alive in a pending render, it'll be freed when GSK is
	 * done; in the meantime we'll briefly overshoot. */
	const gsize prev_bytes = (tile->data ? tile->bytes : 0);
	if (prev_bytes && view->lazy_bytes_used >= prev_bytes)
		view->lazy_bytes_used -= prev_bytes;
	materialise_tile_evict_until(view, tile_bytes);

	guchar *data = malloc(tile_bytes);
	if (!data) return FALSE;

	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	const WORD *src   = gfit->pdata[vport];
	const float *fsrc = gfit->fpdata[vport];
	const guint hd_max = gui.hd_remap_max;

	if (mip == 0) {
		for (int dy = 0; dy < tex_h; dy++) {
			/* display row (y0+dy) → image row (img_h-1-(y0+dy)).  This
			 * is the y-flip the legacy remap performs at fill time.
			 * Guard rows (dy in [th, tex_h)) sample from the next
			 * tile's image rows — same source data the neighbour would
			 * render, so the per-tile mipmap pyramids agree at the
			 * boundary. */
			const int iy = img_h - 1 - (y0 + dy);
			guint32 *row_out = (guint32 *)(data + (size_t)dy * out_w * 4);
			for (int dx = 0; dx < tex_w; dx++) {
				const int ix = x0 + dx;
				const size_t si = (size_t)iy * img_w + ix;
				BYTE val;
				if (gfit->type == DATA_USHORT) {
					const WORD sv = src[si];
					val = hd_mode ? lut[(guint)sv * hd_max / USHRT_MAX] : lut[sv];
				} else {
					val = hd_mode
						? lut[float_to_max_range(fsrc[si], hd_max)]
						: lut[roundf_to_WORD(fsrc[si] * USHRT_MAX_SINGLE)];
				}
				if (inverted) val = UCHAR_MAX - val;
				row_out[dx] = ((guint32)val << 16) | ((guint32)val << 8) | val;
			}
		}
	} else {
		/* Box-average path.  Each output pixel sums up to step*step
		 * LUT-applied source values, then divides by the actual count
		 * (smaller at the right/bottom edge when tex_w/h isn't a
		 * multiple of step). */
		for (int oy = 0; oy < out_h; oy++) {
			guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);
			for (int ox = 0; ox < out_w; ox++) {
				unsigned sum = 0;
				int count = 0;
				for (int sy = 0; sy < step; sy++) {
					const int dy = oy * step + sy;
					if (dy >= tex_h) break;
					const int iy = img_h - 1 - (y0 + dy);
					for (int sx = 0; sx < step; sx++) {
						const int dx = ox * step + sx;
						if (dx >= tex_w) break;
						const int ix = x0 + dx;
						const size_t si = (size_t)iy * img_w + ix;
						BYTE val;
						if (gfit->type == DATA_USHORT) {
							const WORD sv = src[si];
							val = hd_mode ? lut[(guint)sv * hd_max / USHRT_MAX] : lut[sv];
						} else {
							val = hd_mode
								? lut[float_to_max_range(fsrc[si], hd_max)]
								: lut[roundf_to_WORD(fsrc[si] * USHRT_MAX_SINGLE)];
						}
						if (inverted) val = UCHAR_MAX - val;
						sum += val;
						count++;
					}
				}
				const BYTE avg = count > 0 ? (BYTE)(sum / count) : 0;
				row_out[ox] = ((guint32)avg << 16) | ((guint32)avg << 8) | avg;
			}
		}
	}

	/* The new GBytes owns `data` — when the texture's last reference
	 * drops (anywhere from this scope to many frames later, once GSK is
	 * done with it), free(data) will be called automatically. */
	GBytes *bytes = g_bytes_new_with_free_func(data, tile_bytes, free, data);
	GdkTexture *new_tex = gdk_memory_texture_new(out_w, out_h, GDK_MEMORY_B8G8R8X8,
	                                              bytes, out_w * 4);
	g_bytes_unref(bytes);

	if (tile->texture)
		g_object_unref(tile->texture);  /* old bytes survive until GSK is done */
	tile->texture = new_tex;
	tile->data = data;  /* non-owning view; valid for as long as tile->texture lives */
	tile->mip = mip;
	tile->bytes = tile_bytes;
	tile->dirty = FALSE;
	tile->last_used = ++view->lazy_epoch;
	view->lazy_bytes_used += tile_bytes;
	return TRUE;
}

/* Apply the per-channel LUTs and optional ICC transform to a tile rect of
 * the RGB composite viewport.  Same fresh-buffer-each-call pattern as
 * materialise_tile_gray — see that function for the rationale.  idx[c]
 * is the LUT slot to use for channel c (the dispatcher mirrors remap()'s
 * target_index logic so HISTEQ / linked-STF / linear modes correctly
 * read slot 0 instead of stale per-channel slots). */
static gboolean materialise_tile_rgb(struct image_view *view,
                                      int tx, int ty, int mip,
                                      const BYTE * const idx[3],
                                      gboolean hd_mode,
                                      gboolean inverted,
                                      gboolean do_icc) {
	int x0, y0, tw, th, tex_w, tex_h;
	tile_dims_padded(view, tx, ty, &x0, &y0, &tw, &th, &tex_w, &tex_h);
	if (tex_w <= 0 || tex_h <= 0) return TRUE;

	struct image_tile *tile = &view->tiles[ty * view->tile_cols + tx];
	const int step = 1 << mip;
	const int out_w = MAX(1, tex_w >> mip);
	const int out_h = MAX(1, tex_h >> mip);
	const gsize tile_bytes = (gsize)out_w * out_h * 4;

	const gsize prev_bytes = (tile->data ? tile->bytes : 0);
	if (prev_bytes && view->lazy_bytes_used >= prev_bytes)
		view->lazy_bytes_used -= prev_bytes;
	materialise_tile_evict_until(view, tile_bytes);

	guchar *data = malloc(tile_bytes);
	if (!data) return FALSE;

	/* Per-row scratch for the three byte channels (optionally ICC-
	 * transformed in-place).  Sized to the OUTPUT width so the ICC
	 * transform at mip>0 runs once per averaged pixel rather than once per
	 * source pixel — averaging-then-ICC is a small approximation vs ICC-
	 * then-averaging, but the perceptual error is far smaller than the
	 * downsampling itself and saves 4^mip× ICC calls.  At mip=0 it
	 * matches the legacy semantics exactly. */
	BYTE *lb = malloc((size_t)out_w * 3);
	if (!lb) {
		free(data);
		return FALSE;
	}
	BYTE *lb_r = lb;
	BYTE *lb_g = lb + out_w;
	BYTE *lb_b = lb + 2 * out_w;

	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	const guint hd_max = gui.hd_remap_max;

	for (int oy = 0; oy < out_h; oy++) {
		guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);

		for (int ox = 0; ox < out_w; ox++) {
			unsigned sum[3] = { 0, 0, 0 };
			int count = 0;
			for (int sy = 0; sy < step; sy++) {
				const int dy = oy * step + sy;
				if (dy >= tex_h) break;
				const int iy = img_h - 1 - (y0 + dy);
				for (int sx = 0; sx < step; sx++) {
					const int dx = ox * step + sx;
					if (dx >= tex_w) break;
					const int ix = x0 + dx;
					const size_t si = (size_t)iy * img_w + ix;
					for (int c = 0; c < 3; c++) {
						/* Mirror remap()'s indexing exactly — the LUT
						 * already encodes lo/hi, cut_over, and the
						 * display curve, so just look the source pixel
						 * up directly.  No remap_lo subtraction; no
						 * cut_over branch. */
						BYTE val;
						if (gfit->type == DATA_USHORT) {
							const WORD sv = gfit->pdata[c][si];
							val = hd_mode
								? idx[c][(guint)sv * hd_max / USHRT_MAX]
								: idx[c][sv];
						} else {
							const float fv = gfit->fpdata[c][si];
							val = hd_mode
								? idx[c][float_to_max_range(fv, hd_max)]
								: idx[c][roundf_to_WORD(fv * USHRT_MAX_SINGLE)];
						}
						if (inverted) val = UCHAR_MAX - val;
						sum[c] += val;
					}
					count++;
				}
			}
			if (count > 0) {
				lb_r[ox] = (BYTE)(sum[0] / count);
				lb_g[ox] = (BYTE)(sum[1] / count);
				lb_b[ox] = (BYTE)(sum[2] / count);
			} else {
				lb_r[ox] = lb_g[ox] = lb_b[ox] = 0;
			}
		}

		if (do_icc)
			cmsDoTransformLineStride(com.gui_icc.proofing_transform,
				lb, lb, out_w, 1, out_w * 3, out_w * 3, out_w, out_w);

		for (int ox = 0; ox < out_w; ox++) {
			row_out[ox] = ((guint32)lb_r[ox] << 16)
				| ((guint32)lb_g[ox] << 8)
				|  (guint32)lb_b[ox];
		}
	}

	free(lb);

	GBytes *bytes = g_bytes_new_with_free_func(data, tile_bytes, free, data);
	GdkTexture *new_tex = gdk_memory_texture_new(out_w, out_h, GDK_MEMORY_B8G8R8X8,
	                                              bytes, out_w * 4);
	g_bytes_unref(bytes);

	if (tile->texture)
		g_object_unref(tile->texture);
	tile->texture = new_tex;
	tile->data = data;
	tile->mip = mip;
	tile->bytes = tile_bytes;
	tile->dirty = FALSE;
	tile->last_used = ++view->lazy_epoch;
	view->lazy_bytes_used += tile_bytes;
	return TRUE;
}

/* Dispatch to the right per-vport implementation.  Materialises the tile
 * if it's dirty, has no data, or its current mip differs from the one the
 * caller is asking for; no-op otherwise.  `mip` is the downsample level
 * the caller wants — the snapshot path picks it from the current zoom
 * (compute_target_mip); paths that consume tile->data directly at native
 * stride (the Cairo fallback) must pass 0.  Returns FALSE on allocation
 * failure. */
static gboolean materialise_tile(struct image_view *view, int vport,
                                  int tx, int ty, int mip) {
	struct image_tile *tile = &view->tiles[ty * view->tile_cols + tx];
	if (!tile->dirty && tile->texture && tile->mip == mip) {
		tile->last_used = ++view->lazy_epoch;  /* refresh LRU position */
		return TRUE;
	}

	const int target_index = (gui.rendering_mode == STF_DISPLAY && gui.unlink_channels)
		? vport : 0;
	const gboolean hd_mode = (gui.rendering_mode == STF_DISPLAY
		&& gui.use_hd_remap && gfit->type == DATA_FLOAT);

	WORD remap_lo, remap_hi;
	g_mutex_lock(&com.mutex);
	remap_lo = gui.lo;
	remap_hi = gui.hi;
	g_mutex_unlock(&com.mutex);

	/* Negative-view flag — sample the action state once per tile. */
	gboolean neg = FALSE;
	GAction *action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win),
	                                                  "negative-view");
	if (action_neg) {
		GVariant *st = g_action_get_state(action_neg);
		neg = g_variant_get_boolean(st);
		g_variant_unref(st);
	}

	if (vport == RGB_VPORT) {
		/* Pick the LUT slot for each channel the same way remap() picks
		 * its target_index: STF unlinked → per-channel slot; everything
		 * else (HISTEQ, linked STF, LINEAR/LOG/SQRT/SQUARED/ASINH) →
		 * slot 0.  Reading gui.remap_index[c] unconditionally was the
		 * autostretch / histogram bug: slots 1 and 2 contain stale
		 * data for non-unlinked modes. */
		const BYTE *idx[3];
		for (int c = 0; c < 3; c++) {
			const int slot = (gui.rendering_mode == STF_DISPLAY
			                  && gui.unlink_channels) ? c : 0;
			idx[c] = hd_mode ? gui.hd_remap_index[slot]
			                 : gui.remap_index[slot];
		}
		const gboolean do_icc = (gfit->color_managed
			&& com.gui_icc.proofing_transform && !identical
			&& !com.gui_icc.same_primaries);
		return materialise_tile_rgb(view, tx, ty, mip, idx, hd_mode, neg, do_icc);
	}

	if (vport >= 0 && vport <= 2) {
		const BYTE *lut = hd_mode
			? gui.hd_remap_index[target_index]
			: gui.remap_index[target_index];
		return materialise_tile_gray(view, vport, tx, ty, mip, lut, hd_mode,
			neg, gui.cut_over, remap_lo, remap_hi);
	}

	/* Mask vport doesn't yet have a lazy materialiser; fall back to a
	 * cleared tile so the renderer just shows transparent in lazy mode.
	 * In practice the mask viewport is small and rarely triggers lazy. */
	return TRUE;
}

/* ── Full-image proxy (lazy mode) ─────────────────────────────────────────
 *
 * fill_proxy_{gray,rgb} write a point-sampled downscale of the *whole* image
 * into a caller-provided BGRX buffer (out_w × out_h), taking one source pixel
 * per output pixel at stride `step`.  This mirrors the LUT / HD / inverted /
 * ICC handling of materialise_tile_{gray,rgb} (kept in sync with those
 * functions) but skips box-averaging: the proxy is a transient stand-in shown
 * only until real mipmapped tiles arrive, so decimation is visually adequate
 * and keeps the build cost proportional to the proxy size, not the image size.
 * The same bottom-up y-flip as the tile materialisers is applied so the proxy
 * lines up with the tiles and the display_matrix. */
static void fill_proxy_gray(struct image_view *view, int vport,
                            const BYTE *lut, gboolean hd_mode, gboolean inverted,
                            guchar *data, int out_w, int out_h, int step) {
	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	const WORD *src   = gfit->pdata[vport];
	const float *fsrc = gfit->fpdata[vport];
	const guint hd_max = gui.hd_remap_max;
	for (int oy = 0; oy < out_h; oy++) {
		const int iy = img_h - 1 - MIN(oy * step, img_h - 1);  /* y-flip */
		guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);
		for (int ox = 0; ox < out_w; ox++) {
			const int ix = MIN(ox * step, img_w - 1);
			const size_t si = (size_t)iy * img_w + ix;
			BYTE val;
			if (gfit->type == DATA_USHORT) {
				const WORD sv = src[si];
				val = hd_mode ? lut[(guint)sv * hd_max / USHRT_MAX] : lut[sv];
			} else {
				val = hd_mode
					? lut[float_to_max_range(fsrc[si], hd_max)]
					: lut[roundf_to_WORD(fsrc[si] * USHRT_MAX_SINGLE)];
			}
			if (inverted) val = UCHAR_MAX - val;
			row_out[ox] = ((guint32)val << 16) | ((guint32)val << 8) | val;
		}
	}
}

static void fill_proxy_rgb(struct image_view *view, const BYTE * const idx[3],
                           gboolean hd_mode, gboolean inverted, gboolean do_icc,
                           guchar *data, int out_w, int out_h, int step) {
	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	const guint hd_max = gui.hd_remap_max;
	BYTE *lb = malloc((size_t)out_w * 3);
	if (!lb) {                          /* extremely unlikely; show black proxy */
		memset(data, 0, (size_t)out_w * out_h * 4);
		return;
	}
	BYTE *lb_r = lb, *lb_g = lb + out_w, *lb_b = lb + 2 * out_w;
	for (int oy = 0; oy < out_h; oy++) {
		const int iy = img_h - 1 - MIN(oy * step, img_h - 1);  /* y-flip */
		guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);
		for (int ox = 0; ox < out_w; ox++) {
			const int ix = MIN(ox * step, img_w - 1);
			const size_t si = (size_t)iy * img_w + ix;
			BYTE v[3];
			for (int c = 0; c < 3; c++) {
				BYTE val;
				if (gfit->type == DATA_USHORT) {
					const WORD sv = gfit->pdata[c][si];
					val = hd_mode ? idx[c][(guint)sv * hd_max / USHRT_MAX]
					              : idx[c][sv];
				} else {
					const float fv = gfit->fpdata[c][si];
					val = hd_mode ? idx[c][float_to_max_range(fv, hd_max)]
					              : idx[c][roundf_to_WORD(fv * USHRT_MAX_SINGLE)];
				}
				if (inverted) val = UCHAR_MAX - val;
				v[c] = val;
			}
			lb_r[ox] = v[0]; lb_g[ox] = v[1]; lb_b[ox] = v[2];
		}
		if (do_icc)
			cmsDoTransformLineStride(com.gui_icc.proofing_transform,
				lb, lb, out_w, 1, out_w * 3, out_w * 3, out_w, out_w);
		for (int ox = 0; ox < out_w; ox++)
			row_out[ox] = ((guint32)lb_r[ox] << 16)
				| ((guint32)lb_g[ox] << 8) | (guint32)lb_b[ox];
	}
	free(lb);
}

/* (Re)build the lazy-mode proxy for `view` from the current LUTs if it's dirty
 * or missing.  No-op (and frees any stale proxy) in eager mode.  The LUT-slot
 * selection mirrors the materialise_tile dispatcher exactly so the proxy and
 * the real tiles agree on which remap_index slot each channel reads.  Must be
 * called with gui.cairo_mutex held. */
static void ensure_proxy(struct image_view *view, int vport) {
	if (!view->lazy) {
		if (view->proxy) { g_object_unref(view->proxy); view->proxy = NULL; }
		return;
	}
	if (view->proxy && !view->proxy_dirty) return;

	const int img_w = view->buf_stride / 4;
	const int img_h = view->buf_height;
	if (img_w <= 0 || img_h <= 0) return;

	int m = 0;
	while ((img_w >> m) > SIRIL_PROXY_MAX || (img_h >> m) > SIRIL_PROXY_MAX)
		m++;
	const int pw = MAX(1, img_w >> m);
	const int ph = MAX(1, img_h >> m);
	const int step = 1 << m;
	const gsize bytes_len = (gsize)pw * ph * 4;
	guchar *data = malloc(bytes_len);
	if (!data) return;

	const gboolean hd_mode = (gui.rendering_mode == STF_DISPLAY
		&& gui.use_hd_remap && gfit->type == DATA_FLOAT);
	gboolean neg = FALSE;
	GAction *action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win),
	                                                  "negative-view");
	if (action_neg) {
		GVariant *st = g_action_get_state(action_neg);
		neg = g_variant_get_boolean(st);
		g_variant_unref(st);
	}

	if (vport == RGB_VPORT) {
		const BYTE *idx[3];
		for (int c = 0; c < 3; c++) {
			const int slot = (gui.rendering_mode == STF_DISPLAY
			                  && gui.unlink_channels) ? c : 0;
			idx[c] = hd_mode ? gui.hd_remap_index[slot] : gui.remap_index[slot];
		}
		const gboolean do_icc = (gfit->color_managed
			&& com.gui_icc.proofing_transform && !identical
			&& !com.gui_icc.same_primaries);
		fill_proxy_rgb(view, idx, hd_mode, neg, do_icc, data, pw, ph, step);
	} else if (vport >= 0 && vport <= 2) {
		const int target_index = (gui.rendering_mode == STF_DISPLAY
			&& gui.unlink_channels) ? vport : 0;
		const BYTE *lut = hd_mode ? gui.hd_remap_index[target_index]
		                          : gui.remap_index[target_index];
		fill_proxy_gray(view, vport, lut, hd_mode, neg, data, pw, ph, step);
	} else {
		free(data);
		return;
	}

	GBytes *bytes = g_bytes_new_with_free_func(data, bytes_len, free, data);
	GdkTexture *tex = gdk_memory_texture_new(pw, ph, GDK_MEMORY_B8G8R8X8,
	                                          bytes, pw * 4);
	g_bytes_unref(bytes);

	if (view->proxy) g_object_unref(view->proxy);
	view->proxy = tex;
	view->proxy_mip = m;
	view->proxy_w = pw;
	view->proxy_h = ph;
	view->proxy_dirty = FALSE;
}

/* True when a tile materialised at `mip` is good enough to *display* at the
 * current target without re-materialising.  Any tile at least as dense as the
 * target (mip <= target_mip) is fine: an equal-or-finer texture just gets
 * downscaled by GSK and looks correct.  We deliberately do NOT reject
 * finer-than-needed tiles here — re-materialising a visible tile purely to
 * coarsen it (e.g. after a multi-octave zoom-out) is pure latency on the
 * snapshot/render thread for zero visual gain, and was a real source of zoom
 * jerk.  Over-fine tiles are reclaimed lazily for memory via LRU eviction in
 * materialise_tile, not synchronously on the draw path.  Only mip > target_mip
 * (too coarse / aliased) is unacceptable for display. */
static inline gboolean tile_mip_acceptable(int mip, int target_mip) {
	return mip <= target_mip;
}

/* How many levels finer than the target a tile may stay before the *background
 * worker* reclaims it by rebuilding at the target mip.  The display path keeps
 * accepting any-finer tile (tile_mip_acceptable above) — this only governs
 * memory.  A tile at mip M holds 4^(target-M) times the target's bytes, so a
 * full-resolution tile (M=0) left visible after a multi-octave zoom-out is ~17
 * MB; a handful of them blow SIRIL_LAZY_BUDGET, and the LRU then perpetually
 * evicts the tiles it can't fit (a static black block on the last-built corner,
 * un-rebuilt because their fine version won't fit and their coarse version is
 * never made).  Allowing ONE finer level keeps single-octave zoom-in prefetch
 * cheap (the reason compute_target_mip already targets one level fine) while
 * capping the resident working set near the target resolution. */
#define SIRIL_MIP_COARSEN_SLACK 1

/* True when a visible tile is so much finer than the target that the worker
 * should rebuild it coarser to reclaim memory.  Multi-octave jump only: a
 * normal single-step zoom never trips this, so it adds no churn to ordinary
 * panning/zooming. */
static inline gboolean tile_mip_too_fine(int mip, int target_mip) {
	return mip < target_mip - SIRIL_MIP_COARSEN_SLACK;
}

/* LRU eviction: free tile bytes (and their texture) starting from the
 * least-recently-used tile until the viewport has at least `headroom`
 * bytes of budget remaining.  The currently-being-materialised tile is
 * NOT in the LRU yet (caller hasn't bumped its last_used), so it can't
 * accidentally evict itself.  Called with gui.cairo_mutex held. */
static void materialise_tile_evict_until(struct image_view *view, gsize headroom) {
	if (!view->lazy) return;
	while (view->lazy_bytes_used + headroom > SIRIL_LAZY_BUDGET) {
		/* Linear scan for the oldest materialised tile.  N is small
		 * (hundreds at most for typical huge images), and eviction is
		 * a rare event compared to per-pixel work. */
		struct image_tile *victim = NULL;
		const int n = view->tile_cols * view->tile_rows;
		for (int i = 0; i < n; i++) {
			struct image_tile *t = &view->tiles[i];
			if (!t->data) continue;
			if (!victim || t->last_used < victim->last_used)
				victim = t;
		}
		if (!victim) break;  /* nothing to evict */
		tile_release(view, victim);
	}
}

/* ── Background materialise worker ────────────────────────────────────────
 *
 * The heavy per-pixel tile fill runs on a shared single-thread GThreadPool, so
 * neither the fill itself nor cairo_mutex contention can stall the GTK main
 * thread during a zoom.  A worker job for a viewport loops, materialising one
 * still-unready visible tile up to target_mip per iteration and queueing a
 * redraw, so the proxy backdrop is progressively replaced by sharp tiles — this
 * is what makes a fast zoom *resolve*.  It exits as soon as every visible tile
 * is ready, so a settled view (e.g. while only an overlay like the selection
 * rubber-band is being redrawn) schedules no worker at all.
 *
 * There is deliberately NO speculative prefetch of finer mips: stepping visible
 * tiles one level finer costs 4× the memory, which on a large zoomed-out view
 * blows SIRIL_LAZY_BUDGET and sends the LRU into evict/re-materialise thrash
 * that never settles while redraws keep coming (overlay drags).  On-demand fill
 * plus the resident proxy already make zoom smooth without it.
 *
 * Locking (critical):
 *   - Lock order is gfit->rwlock reader THEN gui.cairo_mutex, matching
 *     redraw_mask_idle; no path takes them the other way round.
 *   - The reader lock is held for the WHOLE job and keeps gfit's pixel buffers
 *     alive across the unlocked fill phase (a processing thread takes the
 *     writer lock to mutate gfit, so it blocks until we're done).
 *   - cairo_mutex is held only briefly to pick a job + snapshot its inputs
 *     (phase 1) and to assign the finished texture (phase 3).  The expensive
 *     fill (phase 2) holds NEITHER cairo_mutex — so the snapshot never waits on
 *     us — relying on the reader lock for pixels and on the generation stamp to
 *     detect a grid realloc before assigning.
 *   - LUT contents (gui.remap_index) can be rewritten by a concurrent remap;
 *     that's the same benign race the inline path already tolerates (remap
 *     builds the LUT before taking cairo_mutex), self-healed by remap's own
 *     dirty+redraw. */

static GThreadPool *materialise_pool = NULL;

/* Select the LUT pointer(s) for a viewport, mirroring the materialise_tile
 * dispatcher's slot logic exactly.  idx[0..2] for the RGB composite; idx[0]
 * only for a gray channel.  Returns hd_mode and sets *out_do_icc.  Reads gui.*
 * / gfit->* and so must be called with gui.cairo_mutex held. */
static gboolean pick_render_luts(int vport, const BYTE *idx[3],
                                 gboolean *out_do_icc) {
	const gboolean hd_mode = (gui.rendering_mode == STF_DISPLAY
		&& gui.use_hd_remap && gfit->type == DATA_FLOAT);
	idx[0] = idx[1] = idx[2] = NULL;
	if (vport == RGB_VPORT) {
		for (int c = 0; c < 3; c++) {
			const int slot = (gui.rendering_mode == STF_DISPLAY
			                  && gui.unlink_channels) ? c : 0;
			idx[c] = hd_mode ? gui.hd_remap_index[slot] : gui.remap_index[slot];
		}
		*out_do_icc = (gfit->color_managed && com.gui_icc.proofing_transform
			&& !identical && !com.gui_icc.same_primaries);
	} else {
		const int ti = (gui.rendering_mode == STF_DISPLAY
		                && gui.unlink_channels) ? vport : 0;
		idx[0] = hd_mode ? gui.hd_remap_index[ti] : gui.remap_index[ti];
		*out_do_icc = FALSE;
	}
	return hd_mode;
}

/* Pure pixel producers used by the worker's phase 2.  These hold no lock and
 * touch no view/tile/gui state — everything comes in as parameters.  They MUST
 * stay byte-for-byte equivalent to materialise_tile_gray / materialise_tile_rgb
 * (the inline path) so a tile looks identical whether the main thread or the
 * worker built it. */
static void wk_fill_gray(int gfit_type, const WORD *src, const float *fsrc,
                         int img_w, int img_h, int x0, int y0,
                         int tex_w, int tex_h, int mip, const BYTE *lut,
                         guint hd_max, gboolean hd_mode, gboolean inverted,
                         guchar *data, int out_w, int out_h) {
	const int step = 1 << mip;
	if (mip == 0) {
		for (int dy = 0; dy < tex_h; dy++) {
			const int iy = img_h - 1 - (y0 + dy);
			guint32 *row_out = (guint32 *)(data + (size_t)dy * out_w * 4);
			for (int dx = 0; dx < tex_w; dx++) {
				const int ix = x0 + dx;
				const size_t si = (size_t)iy * img_w + ix;
				BYTE val;
				if (gfit_type == DATA_USHORT) {
					const WORD sv = src[si];
					val = hd_mode ? lut[(guint)sv * hd_max / USHRT_MAX] : lut[sv];
				} else {
					val = hd_mode ? lut[float_to_max_range(fsrc[si], hd_max)]
					              : lut[roundf_to_WORD(fsrc[si] * USHRT_MAX_SINGLE)];
				}
				if (inverted) val = UCHAR_MAX - val;
				row_out[dx] = ((guint32)val << 16) | ((guint32)val << 8) | val;
			}
		}
	} else {
		for (int oy = 0; oy < out_h; oy++) {
			guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);
			for (int ox = 0; ox < out_w; ox++) {
				unsigned sum = 0;
				int count = 0;
				for (int sy = 0; sy < step; sy++) {
					const int dy = oy * step + sy;
					if (dy >= tex_h) break;
					const int iy = img_h - 1 - (y0 + dy);
					for (int sx = 0; sx < step; sx++) {
						const int dx = ox * step + sx;
						if (dx >= tex_w) break;
						const int ix = x0 + dx;
						const size_t si = (size_t)iy * img_w + ix;
						BYTE val;
						if (gfit_type == DATA_USHORT) {
							const WORD sv = src[si];
							val = hd_mode ? lut[(guint)sv * hd_max / USHRT_MAX] : lut[sv];
						} else {
							val = hd_mode ? lut[float_to_max_range(fsrc[si], hd_max)]
							              : lut[roundf_to_WORD(fsrc[si] * USHRT_MAX_SINGLE)];
						}
						if (inverted) val = UCHAR_MAX - val;
						sum += val;
						count++;
					}
				}
				const BYTE avg = count > 0 ? (BYTE)(sum / count) : 0;
				row_out[ox] = ((guint32)avg << 16) | ((guint32)avg << 8) | avg;
			}
		}
	}
}

static void wk_fill_rgb(int gfit_type, const WORD *const psrc[3],
                        const float *const fpsrc[3], int img_w, int img_h,
                        int x0, int y0, int tex_w, int tex_h, int mip,
                        const BYTE *const idx[3], guint hd_max, gboolean hd_mode,
                        gboolean inverted, gboolean do_icc,
                        guchar *data, int out_w, int out_h) {
	const int step = 1 << mip;
	BYTE *lb = malloc((size_t)out_w * 3);
	if (!lb) {
		memset(data, 0, (size_t)out_w * out_h * 4);
		return;
	}
	BYTE *lb_r = lb, *lb_g = lb + out_w, *lb_b = lb + 2 * out_w;
	for (int oy = 0; oy < out_h; oy++) {
		guint32 *row_out = (guint32 *)(data + (size_t)oy * out_w * 4);
		for (int ox = 0; ox < out_w; ox++) {
			unsigned sum[3] = { 0, 0, 0 };
			int count = 0;
			for (int sy = 0; sy < step; sy++) {
				const int dy = oy * step + sy;
				if (dy >= tex_h) break;
				const int iy = img_h - 1 - (y0 + dy);
				for (int sx = 0; sx < step; sx++) {
					const int dx = ox * step + sx;
					if (dx >= tex_w) break;
					const int ix = x0 + dx;
					const size_t si = (size_t)iy * img_w + ix;
					for (int c = 0; c < 3; c++) {
						BYTE val;
						if (gfit_type == DATA_USHORT) {
							const WORD sv = psrc[c][si];
							val = hd_mode ? idx[c][(guint)sv * hd_max / USHRT_MAX]
							              : idx[c][sv];
						} else {
							const float fv = fpsrc[c][si];
							val = hd_mode ? idx[c][float_to_max_range(fv, hd_max)]
							              : idx[c][roundf_to_WORD(fv * USHRT_MAX_SINGLE)];
						}
						if (inverted) val = UCHAR_MAX - val;
						sum[c] += val;
					}
					count++;
				}
			}
			if (count > 0) {
				lb_r[ox] = (BYTE)(sum[0] / count);
				lb_g[ox] = (BYTE)(sum[1] / count);
				lb_b[ox] = (BYTE)(sum[2] / count);
			} else {
				lb_r[ox] = lb_g[ox] = lb_b[ox] = 0;
			}
		}
		if (do_icc)
			cmsDoTransformLineStride(com.gui_icc.proofing_transform,
				lb, lb, out_w, 1, out_w * 3, out_w * 3, out_w, out_w);
		for (int ox = 0; ox < out_w; ox++)
			row_out[ox] = ((guint32)lb_r[ox] << 16)
				| ((guint32)lb_g[ox] << 8) | (guint32)lb_b[ox];
	}
	free(lb);
}

/* Assign a worker-built texture to a tile and update the LRU budget — the
 * mutating tail of materialise_tile_gray/rgb, factored out for the worker's
 * phase 3.  Must be called with gui.cairo_mutex held. */
static void wk_assign_tile(struct image_view *view, int tx, int ty, int mip,
                           GdkTexture *tex, guchar *data, gsize tile_bytes) {
	struct image_tile *t = &view->tiles[ty * view->tile_cols + tx];
	const gsize prev = (t->data ? t->bytes : 0);
	if (prev && view->lazy_bytes_used >= prev)
		view->lazy_bytes_used -= prev;
	materialise_tile_evict_until(view, tile_bytes);
	if (t->texture)
		g_object_unref(t->texture);
	t->texture = tex;
	t->data = data;
	t->mip = mip;
	t->bytes = tile_bytes;
	t->dirty = FALSE;
	t->last_used = ++view->lazy_epoch;
	view->lazy_bytes_used += tile_bytes;
}

/* Runs on the main thread (scheduled from the worker) to request a redraw of
 * the viewport whose deferred tile just landed. */
static gboolean wk_queue_draw_idle(gpointer p) {
	const int vport = GPOINTER_TO_INT(p);
	if (vport >= 0 && vport < MAXVPORT && gui.view[vport].drawarea)
		gtk_widget_queue_draw(gui.view[vport].drawarea);
	return G_SOURCE_REMOVE;
}

static void materialise_worker(gpointer data, gpointer user) {
	(void) user;
	const int vport = GPOINTER_TO_INT(data);
	/* Only the gray channels (0..2) and the RGB composite have lazy
	 * materialisers; the mask vport has none (see materialise_tile). */
	if (vport < 0 || vport > RGB_VPORT)
		return;

	for (;;) {
		if (!gfit) {
			g_mutex_lock(&gui.cairo_mutex);
			gui.view[vport].worker_active = FALSE;
			g_mutex_unlock(&gui.cairo_mutex);
			return;
		}

		g_rw_lock_reader_lock(&gfit->rwlock);
		g_mutex_lock(&gui.cairo_mutex);
		struct image_view *view = &gui.view[vport];

		if (!view->lazy || !view->tiles
		    || g_atomic_int_get(&gui.suppress_drawarea_redraw)) {
			view->worker_active = FALSE;
			g_mutex_unlock(&gui.cairo_mutex);
			g_rw_lock_reader_unlock(&gfit->rwlock);
			return;
		}

		const int cols = view->tile_cols;
		const int rows = view->tile_rows;
		const int target = view->wk_target_mip;

		/* Phase 1a — pick the next still-unready visible tile (raster order
		 * within the range the last snapshot recorded). */
		int jtx = -1, jty = -1;
		const int jmip = target;
		gboolean found = FALSE;
		if (view->vis_valid) {
			const int tx0 = CLAMP(view->vis_tx0, 0, cols - 1);
			const int tx1 = CLAMP(view->vis_tx1, 0, cols - 1);
			const int ty0 = CLAMP(view->vis_ty0, 0, rows - 1);
			const int ty1 = CLAMP(view->vis_ty1, 0, rows - 1);
			/* Pass 1 — reclaim before consume.  Coarsen an over-fine resident
			 * tile first: rebuilding a tile that is 4^n-too-large at the target
			 * mip FREES its excess bytes with no eviction (wk_assign_tile credits
			 * the old size before reserving the new).  This is what stops a
			 * budget-full view from livelocking: building a missing tile would
			 * otherwise have to evict, but the giant over-fine tiles hogging the
			 * budget are immune to LRU (visible → last_used bumped every frame),
			 * so the worker would ping-pong the few small tiles forever.  Draining
			 * the over-fine tiles first opens enough budget for every missing tile
			 * to build cleanly.
			 *
			 * The threshold adapts to memory pressure: with headroom we only
			 * reclaim tiles finer than the prefetch level (target − SLACK), so a
			 * single-octave zoom-in still finds tiles ready; when the budget is
			 * nearly full we reclaim every tile finer than the target itself, so
			 * the resident set converges to the target-resolution working set —
			 * which always fits as long as that set alone fits the budget. */
			const gboolean pressure =
				view->lazy_bytes_used > (SIRIL_LAZY_BUDGET * 3) / 4;
			const int reclaim_below =
				pressure ? target : target - SIRIL_MIP_COARSEN_SLACK;
			for (int ty = ty0; !found && ty <= ty1; ty++) {
				for (int tx = tx0; tx <= tx1; tx++) {
					const struct image_tile *t = &view->tiles[ty * cols + tx];
					if (t->texture && t->mip < reclaim_below) {
						jtx = tx; jty = ty;
						found = TRUE;
						break;
					}
				}
			}
			/* Pass 2 — otherwise the first not-ready tile (missing / dirty /
			 * too coarse), raster order. */
			for (int ty = ty0; !found && ty <= ty1; ty++) {
				for (int tx = tx0; tx <= tx1; tx++) {
					const struct image_tile *t = &view->tiles[ty * cols + tx];
					const gboolean ready = t->texture && !t->dirty
						&& tile_mip_acceptable(t->mip, target)
						&& !tile_mip_too_fine(t->mip, target);
					if (!ready) {
						jtx = tx; jty = ty;
						found = TRUE;
						break;
					}
				}
			}
		}
		if (!found) {
			view->worker_active = FALSE;
			g_mutex_unlock(&gui.cairo_mutex);
			g_rw_lock_reader_unlock(&gfit->rwlock);
			return;
		}

		/* Phase 1b — snapshot every input the fill needs, so phase 2 touches no
		 * shared mutable state.  LUT *pointers* are captured here (contents may
		 * race a remap, see the section note); gfit pixel pointers stay valid
		 * under the reader lock. */
		const guint gen = view->generation;
		const guint seq = view->invalidate_seq;
		int gx0, gy0, gtw, gth, gtex_w, gtex_h;
		tile_dims_padded(view, jtx, jty, &gx0, &gy0, &gtw, &gth, &gtex_w, &gtex_h);
		const int img_w = view->buf_stride / 4;
		const int img_h = view->buf_height;
		const int out_w = MAX(1, gtex_w >> jmip);
		const int out_h = MAX(1, gtex_h >> jmip);
		const gsize tile_bytes = (gsize)out_w * out_h * 4;
		const gboolean neg = view->render_neg;
		const int gtype = gfit->type;
		const guint hd_max = gui.hd_remap_max;
		const BYTE *idx[3];
		gboolean do_icc = FALSE;
		const gboolean hd_mode = pick_render_luts(vport, idx, &do_icc);
		const WORD *psrc[3] = { gfit->pdata[0], gfit->pdata[1], gfit->pdata[2] };
		const float *fpsrc[3] = { gfit->fpdata[0], gfit->fpdata[1], gfit->fpdata[2] };

		g_mutex_unlock(&gui.cairo_mutex);   /* keep the reader lock for the fill */

		/* Phase 2 — heavy fill, NO cairo_mutex held. */
		guchar *buf = malloc(tile_bytes);
		if (!buf) {
			g_rw_lock_reader_unlock(&gfit->rwlock);
			g_mutex_lock(&gui.cairo_mutex);
			gui.view[vport].worker_active = FALSE;
			g_mutex_unlock(&gui.cairo_mutex);
			return;
		}
		if (vport == RGB_VPORT)
			wk_fill_rgb(gtype, psrc, fpsrc, img_w, img_h, gx0, gy0,
				gtex_w, gtex_h, jmip, idx, hd_max, hd_mode, neg, do_icc,
				buf, out_w, out_h);
		else
			wk_fill_gray(gtype, psrc[vport], fpsrc[vport], img_w, img_h,
				gx0, gy0, gtex_w, gtex_h, jmip, idx[0], hd_max, hd_mode, neg,
				buf, out_w, out_h);

		GBytes *bytes = g_bytes_new_with_free_func(buf, tile_bytes, free, buf);
		GdkTexture *tex = gdk_memory_texture_new(out_w, out_h,
			GDK_MEMORY_B8G8R8X8, bytes, out_w * 4);
		g_bytes_unref(bytes);

		/* Phase 3 — assign, but only if neither the grid moved (generation) nor
		 * the LUT changed (invalidate_seq) under us.  A seq mismatch means a
		 * remap re-dirtied this tile mid-fill: drop our stale texture and leave
		 * the dirty flag set so the next loop iteration rebuilds it with the
		 * current LUT. */
		g_mutex_lock(&gui.cairo_mutex);
		gboolean assigned = FALSE;
		if (view->generation == gen && view->invalidate_seq == seq
		    && view->lazy && view->tiles
		    && view->tile_cols == cols && view->tile_rows == rows) {
			wk_assign_tile(view, jtx, jty, jmip, tex, buf, tile_bytes);
			assigned = TRUE;
		}
		g_mutex_unlock(&gui.cairo_mutex);
		g_rw_lock_reader_unlock(&gfit->rwlock);

		if (!assigned)
			g_object_unref(tex);   /* discard; frees buf via the GBytes free func */
		else
			g_idle_add(wk_queue_draw_idle, GINT_TO_POINTER(vport));

		/* Loop for the next job. */
	}
}

/* Kick the materialise worker for a viewport if one isn't already churning it
 * and there's deferred work to do (a visible tile still showing the proxy).
 * Called with gui.cairo_mutex held (from the snapshot); the worker re-checks
 * under the lock before doing anything. */
static void schedule_view_worker(struct image_view *view, int vport,
                                 gboolean have_deferred) {
	if (!view->lazy) return;
	if (vport < 0 || vport > RGB_VPORT) return;  /* mask vport has no materialiser */
	if (!have_deferred) return;
	if (view->worker_active) return;
	if (!materialise_pool) {
		/* One dedicated worker thread: the job function loops until the view is
		 * drained, and all jobs contend on cairo_mutex anyway, so extra threads
		 * would only add contention. */
		materialise_pool = g_thread_pool_new(materialise_worker, NULL,
		                                     1, TRUE, NULL);
		if (!materialise_pool) return;
	}
	view->worker_active = TRUE;
	g_thread_pool_push(materialise_pool, GINT_TO_POINTER(vport), NULL);
}

void check_gfit_profile_identical_to_monitor() {
	if (!com.headless && gfit->icc_profile && gfit->color_managed)
		identical = profiles_identical(gfit->icc_profile, com.gui_icc.monitor);
	siril_log_debug("gfit profile identical to monitor profile: %d\n", identical);
}

static void remaprgb(void) {
	guint32 *dst;
	const guint32 *bufr, *bufg, *bufb;
	gint i;
	int nbdata;

	siril_log_debug("remaprgb\n");
	if (!isrgb(gfit))
		return;

	struct image_view *rgbview = &gui.view[RGB_VPORT];
	if (allocate_full_surface(rgbview))
		return;

	g_mutex_lock(&gui.cairo_mutex);

	/* Lazy mode: no gray bufs to composite from.  materialise_tile_rgb
	 * reads gfit + the per-channel LUTs directly when it fills each
	 * RGB tile, so we just invalidate here. */
	if (rgbview->lazy) {
		view_refresh_tile_textures(rgbview);
		g_mutex_unlock(&gui.cairo_mutex);
		invalidate_image_render_cache(RGB_VPORT);
		return;
	}

	bufr = (const guint32*) gui.view[RED_VPORT].buf;
	bufg = (const guint32*) gui.view[GREEN_VPORT].buf;
	bufb = (const guint32*) gui.view[BLUE_VPORT].buf;
	if (bufr == NULL || bufg == NULL || bufb == NULL) {
		siril_log_debug("remaprgb: gray buffers not allocated for display\n");
		g_mutex_unlock(&gui.cairo_mutex);
		return;
	}
	dst = (guint32*) rgbview->buf;	// index is j
	nbdata = (rgbview->buf_stride / 4) * rgbview->buf_height;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; ++i) {
		dst[i] = (bufr[i] & 0xFF0000) | (bufg[i] & 0xFF00) | (bufb[i] & 0xFF);
	}

	view_refresh_tile_textures(rgbview);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(RGB_VPORT);
}

void allocate_hd_remap_indices() {
	gui.hd_remap_max = 1 << (guint) com.pref.hd_bitdepth;
	for (unsigned i=0; i < 3; i++) {
		if (gui.hd_remap_index[i] != NULL)
			free(gui.hd_remap_index[i]);
		gui.hd_remap_index[i] = (BYTE*) calloc(gui.hd_remap_max + 1, sizeof(BYTE));
		if (gui.hd_remap_index[i] == NULL) {
			siril_log_error(_("Error: memory allocaton failure when instantiating HD LUTs. Reverting to standard 16 bit LUTs.\n"));
			gui.use_hd_remap = FALSE;
			image_display_init_statics();
			siril_toggle_set_active(GTK_WIDGET(imgdisp_autohd_item), FALSE);
			hd_remap_indices_cleanup();
			return;
		}
	}
}

void hd_remap_indices_cleanup() {
	for (unsigned i=0 ; i < 3; i++) {
		if (gui.hd_remap_index[i] != NULL) {
			free(gui.hd_remap_index[i]);
			gui.hd_remap_index[i] = NULL;
		}
	}
}

static int make_index_for_current_display(int vport);

static int make_hd_index_for_current_display(int vport);

static int make_index_for_rainbow(BYTE index[][3]);

static void remap_mask(mask_t *mask) {
	siril_log_debug("mask remap\n");

	int vport = MASK_VPORT;
	struct image_view *view = &gui.view[vport];
	if (allocate_full_surface(view))
		return;

	g_mutex_lock(&gui.cairo_mutex);

	/* Lazy mode: no view->buf to write into.  Mask materialisation is
	 * stubbed out for now (see materialise_tile); just refresh textures
	 * (which marks tiles dirty) and bail. */
	if (view->lazy) {
		view_refresh_tile_textures(view);
		g_mutex_unlock(&gui.cairo_mutex);
		invalidate_image_render_cache(vport);
		test_and_allocate_reference_image(vport);
		return;
	}

	unsigned char *dst_base = view->buf;
	int dst_stride = view->buf_stride;

	// Mask setup
	void *src_data = mask->data;
	guint rx = gfit->rx;
	guint ry = gfit->ry;
	int bitpix = mask->bitpix;

	// We switch ONCE. This requires duplicating the loop code,
	// but ensures the tightest possible inner loops for the compiler to vectorize.
	switch (bitpix) {
		case 8: {
			uint8_t *s8 = (uint8_t *)src_data;
			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				// Pre-calculate pointers for this row
				uint8_t *src_row = s8 + (y * rx);
				// Vertical flip logic
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					uint8_t val = src_row[x];
					// Alpha is unused/ignored in Cairo RGB24 (upper 8 bits)
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		case 16: {
			uint16_t *s16 = (uint16_t *)src_data;
			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				uint16_t *src_row = s16 + (y * rx);
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					uint32_t val = ((uint32_t)src_row[x] * 255 + 32895) >> 16;
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		case 32: {
			float *sf = (float *)src_data;

			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				float *src_row = sf + (y * rx);

				// Destination pointer calculation (Vertical Flip)
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					// Scale
					float v = src_row[x] * UCHAR_MAX_SINGLE;

					// Round & Clamp
					uint8_t val = roundf_to_BYTE(v);

					// Pack to ARGB (Alpha unused)
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		default: {
			siril_log_debug("Error: invalid mask bitpix\n");
			break;
		}
	}

	view_refresh_tile_textures(view);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(vport);
	test_and_allocate_reference_image(vport);
}

/*
 * Idle wrapper for set_viewer_mode_widgets_sensitive().
 *
 * remap() and remap_all_vports() may be called from a non-GTK thread
 * (via notify_gfit_data_modified()).  Widget-sensitivity changes must only
 * happen on the GTK main thread, so we dispatch them as idle callbacks.
 * The gint value (0 or 1) is passed via GINT_TO_POINTER / GPOINTER_TO_INT.
 */
static gboolean viewer_mode_sensitive_idle(gpointer data) {
	set_viewer_mode_widgets_sensitive(GPOINTER_TO_INT(data));
	return G_SOURCE_REMOVE;
}

// remapping one vport at a time is used for DISPLAY_STF and DISPLAY_HISTEQ
static void remap(int vport) {
	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	BYTE *dst, *index, rainbow_index[UCHAR_MAX + 1][3];
	WORD *src;
	float *fsrc;
	gboolean inverted;
	siril_log_debug("HISTEQ / STF remap %d\n", vport);
	if (vport == RGB_VPORT) {
		remaprgb();
		return;
	}
	if (gfit->type == DATA_UNSUPPORTED) {
		siril_log_debug("data is not loaded yet\n");
		return;
	}

	struct image_view *view = &gui.view[vport];
	if (allocate_full_surface(view))
		return;

	/* Cache the two GAction pointers on first call.
	 * g_action_map_lookup_action() is not thread-safe; caching here (initialised
	 * from the GTK main thread on first image display) avoids calling it from
	 * worker threads.  g_action_get_state() on a GSimpleAction is thread-safe
	 * (atomic pointer read) and is therefore safe to call from any thread. */
	static GAction *action_neg_cached = NULL;
	static GAction *action_color_cached = NULL;
	if (action_neg_cached == NULL) {
		image_display_init_statics();
		action_neg_cached   = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
		action_color_cached = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "color-map");
	}
	GVariant *neg_state = g_action_get_state(action_neg_cached);
	inverted = g_variant_get_boolean(neg_state);
	g_variant_unref(neg_state);
	neg_state = NULL;

	if (gui.rendering_mode == HISTEQ_DISPLAY) {
		double hist_sum, nb_pixels;
		size_t i, hist_nb_bins;
		gsl_histogram *histo;
		compute_histo_for_fit(gfit);
		histo = com.layers_hist[vport];
		hist_nb_bins = gsl_histogram_bins(histo);
		nb_pixels = (double)(gfit->rx * gfit->ry);
		// build the remap_index
		index = gui.remap_index[0];
		index[0] = 0;
		hist_sum = gsl_histogram_get(histo, 0);
		for (i = 1; i < hist_nb_bins; i++) {
			hist_sum += gsl_histogram_get(histo, i);
			index[i] = round_to_BYTE((hist_sum / nb_pixels) * UCHAR_MAX_DOUBLE);
		}

		last_mode = gui.rendering_mode;
		histo = com.layers_hist[vport];
		/* Widget-sensitivity changes must happen on the GTK main thread. */
		siril_add_idle(viewer_mode_sensitive_idle, GINT_TO_POINTER(FALSE));
	} else {
		if (gui.rendering_mode == STF_DISPLAY && !stf_computed) {
			if (gui.unlink_channels)
				find_unlinked_midtones_balance_default(gfit, stf);
			else find_linked_midtones_balance_default(gfit, stf);
			stf_computed = TRUE;
		}
		if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type == DATA_FLOAT) {
			make_hd_index_for_current_display(vport);
		}
		else
			make_index_for_current_display(vport);
		siril_add_idle(viewer_mode_sensitive_idle,
		               GINT_TO_POINTER(gui.rendering_mode != STF_DISPLAY));
	}

	src = gfit->pdata[vport];
	fsrc = gfit->fpdata[vport];

	g_mutex_lock(&gui.cairo_mutex);

	/* Lazy mode: gui.remap_index[vport] is now populated for the new
	 * display state; mark all tiles dirty and let snapshot materialise
	 * the visible ones on demand. */
	if (view->lazy) {
		view_refresh_tile_textures(view);
		g_mutex_unlock(&gui.cairo_mutex);
		invalidate_image_render_cache(vport);
		test_and_allocate_reference_image(vport);
		return;
	}

	dst = view->buf;

	GVariant *rainbow_state = g_action_get_state(action_color_cached);
	color_map color = g_variant_get_boolean(rainbow_state);
	g_variant_unref(rainbow_state);
	rainbow_state = NULL;

	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.unlink_channels ? vport : 0;

	gboolean hd_mode = (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type == DATA_FLOAT);
	if (hd_mode) {
		index = gui.hd_remap_index[target_index];
	}
	else
		index = gui.remap_index[target_index];

	// Check if mask overlay is active
	gboolean apply_mask = (gfit->mask != NULL && gfit->mask_active && com.pref.gui.mask_tints_vports);
	uint8_t *mask_u8 = NULL;
	uint16_t *mask_u16 = NULL;
	float *mask_f32 = NULL;
	int mask_bitpix = 0;

	if (apply_mask) {
		mask_bitpix = gfit->mask->bitpix;
		switch (mask_bitpix) {
			case 8:
				mask_u8 = (uint8_t *)gfit->mask->data;
				break;
			case 16:
				mask_u16 = (uint16_t *)gfit->mask->data;
				break;
			case 32:
				mask_f32 = (float *)gfit->mask->data;
				break;
			default:
				apply_mask = FALSE; // Invalid bitpix, disable mask
				break;
		}
	}

	const guint width = gfit->rx;
	const guint height = gfit->ry;

/* Macro to compute a pixel from source index src_i and write it to dst+dst_i */
#define REMAP_WRITE_PIXEL(src_i, dst_i) \
	do { \
		BYTE dpv; \
		if (gfit->type == DATA_USHORT) { \
			const WORD sv = src[(src_i)]; \
			dpv = hd_mode ? index[sv * gui.hd_remap_max / USHRT_MAX] : index[sv]; \
		} else { \
			dpv = hd_mode ? index[float_to_max_range(fsrc[(src_i)], gui.hd_remap_max)] \
			              : index[roundf_to_WORD(fsrc[(src_i)] * USHRT_MAX_SINGLE)]; \
		} \
		if (inverted) dpv = UCHAR_MAX - dpv; \
		if (apply_mask) { \
			BYTE mv; \
			if (mask_bitpix == 8) mv = mask_u8[(src_i)]; \
			else if (mask_bitpix == 16) mv = mask_u16[(src_i)] >> 8; \
			else mv = (BYTE)(mask_f32[(src_i)] * UCHAR_MAX); \
			const uint16_t imv = UCHAR_MAX - mv; \
			if (color == RAINBOW_COLOR) { \
				const uint16_t rs = imv + rainbow_index[dpv][0]; \
				*(guint32*)(dst + (dst_i)) = \
					((rs > UCHAR_MAX ? UCHAR_MAX : (BYTE)rs) << 16) | \
					((rainbow_index[dpv][1] * mv) >> 8) << 8 | \
					((rainbow_index[dpv][2] * mv) >> 8); \
			} else { \
				const uint16_t rs = imv + dpv; \
				*(guint32*)(dst + (dst_i)) = \
					((rs > UCHAR_MAX ? UCHAR_MAX : (BYTE)rs) << 16) | \
					((dpv * mv) >> 8) << 8 | ((dpv * mv) >> 8); \
			} \
		} else { \
			if (color == RAINBOW_COLOR) \
				*(guint32*)(dst + (dst_i)) = rainbow_index[dpv][0] << 16 | \
				                             rainbow_index[dpv][1] << 8 | \
				                             rainbow_index[dpv][2]; \
			else \
				*(guint32*)(dst + (dst_i)) = dpv << 16 | dpv << 8 | dpv; \
		} \
	} while (0)

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (guint y = 0; y < height; y++) {
		const guint src_row_start = y * width;
		const guint dst_row_start = (height - 1 - y) * width * 4;
		for (guint x = 0; x < width; ++x) {
			REMAP_WRITE_PIXEL(src_row_start + x, dst_row_start + x * 4);
		}
	}
#undef REMAP_WRITE_PIXEL

	view_refresh_tile_textures(view);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(vport);
	test_and_allocate_reference_image(vport);
}

static void remap_all_vports() {
	gboolean inverted;
	/* Snapshot gui.hi/gui.lo: init_layers_hi_and_lo_values() may write them
	 * from the worker thread concurrently while this path runs from
	 * remap_all() called directly on the worker (via notify_gfit_data_modified). */
	g_mutex_lock(&com.mutex);
	WORD remap_hi = gui.hi;
	WORD remap_lo = gui.lo;
	g_mutex_unlock(&com.mutex);

	/* Cache GAction pointers on first call — same reasoning as in remap() above. */
	static GAction *action_neg_cached = NULL;
	static GAction *action_color_cached = NULL;
	if (action_neg_cached == NULL) {
		image_display_init_statics();
		action_neg_cached   = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
		action_color_cached = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "color-map");
	}
	struct image_view *view[3] = { &gui.view[0], &gui.view[1], &gui.view[2] };
	GVariant *state_neg = g_action_get_state(action_neg_cached);
	inverted = g_variant_get_boolean(state_neg);
	g_variant_unref(state_neg);
	state_neg = NULL;
	// We are now dealing with a 3-channel image

	// Check if we need a rainbow color map
	BYTE rainbow_index[UCHAR_MAX + 1][3];
	GVariant* rainbow_state = g_action_get_state(action_color_cached);
	color_map color = g_variant_get_boolean(rainbow_state);
	g_variant_unref(rainbow_state);
	rainbow_state = NULL;
	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);

	// Check if mask overlay is active
	gboolean apply_mask = (gfit->mask != NULL && gfit->mask_active && com.pref.gui.mask_tints_vports);
	uint8_t *mask_u8 = NULL;
	uint16_t *mask_u16 = NULL;
	float *mask_f32 = NULL;
	int mask_bitpix = 0;

	if (apply_mask) {
		mask_bitpix = gfit->mask->bitpix;
		switch (mask_bitpix) {
			case 8:
				mask_u8 = (uint8_t *)gfit->mask->data;
				break;
			case 16:
				mask_u16 = (uint16_t *)gfit->mask->data;
				break;
			case 32:
				mask_f32 = (float *)gfit->mask->data;
				break;
			default:
				apply_mask = FALSE; // Invalid bitpix, disable mask
				break;
		}
	}

	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	BYTE *dst[3], *index[3];
	WORD *src[3];
	float *fsrc[3];

	if (gfit->type == DATA_UNSUPPORTED) {
		siril_log_debug("data is not loaded yet\n");
		return;
	}

	lock_display_transform();
	if (gfit->color_managed) {
		// Set the transform in case it is missing
		if (!com.gui_icc.proofing_transform) {
			com.gui_icc.proofing_transform = initialize_proofing_transform();
			com.gui_icc.profile_changed = TRUE;
		}
		if (com.gui_icc.profile_changed) {
			com.gui_icc.same_primaries = same_primaries(gfit->icc_profile, com.gui_icc.monitor, (com.gui_icc.soft_proof && com.pref.icc.soft_proofing_profile_active) ? com.gui_icc.soft_proof : NULL);
//			com.gui_icc.same_primaries = FALSE;
			check_gfit_profile_identical_to_monitor();
			// Calling color_manage() like this updates the color management button tooltip
			color_manage(gfit, gfit->color_managed);
			if (is_preview_active())
				copy_gfit_icc_to_backup();
		}
	}

	make_index_for_current_display(0);
	index[0] = gui.remap_index[0];
	if (gfit->color_managed) {
		for (int i = 1 ; i < 3 ; i++) {
			make_index_for_current_display(i);
			index[i] = gui.remap_index[i];
		}
	}
	unlock_display_transform();

	com.gui_icc.profile_changed = FALSE;
	/* Widget-sensitivity changes must happen on the GTK main thread. */
	siril_add_idle(viewer_mode_sensitive_idle, GINT_TO_POINTER(TRUE));

	last_mode = gui.rendering_mode;

	// Allocate surfaces and capture src/fsrc outside the lock.
	// dst[] pointers are captured below, inside the lock, to prevent a stale-
	// pointer race: a concurrent allocate_full_surface() could realloc view->buf
	// between its internal unlock and the cairo_mutex acquisition below.
	for (int i = 0 ; i < 3 ; i++) {
		src[i] = gfit->pdata[i];
		fsrc[i] = gfit->fpdata[i];
		if (allocate_full_surface(view[i]))
			return;
	}
	g_mutex_lock(&gui.cairo_mutex);

	/* Lazy mode: per-channel LUTs (gui.remap_index[0..2]) are now valid.
	 * Mark all gray tiles dirty so the snapshot path can materialise
	 * the visible ones; skip the eager pixel pass entirely. */
	if (view[0]->lazy) {
		for (int i = 0; i < 3; i++)
			view_refresh_tile_textures(view[i]);
		g_mutex_unlock(&gui.cairo_mutex);
		for (int vport = 0; vport < 3; vport++) {
			invalidate_image_render_cache(vport);
			test_and_allocate_reference_image(vport);
		}
		return;
	}

	for (int i = 0 ; i < 3 ; i++)
		dst[i] = view[i]->buf;

	int norm = (int) get_normalized_value(gfit);
	const guint width = gfit->rx;
	const guint height = gfit->ry;

	/* Render is always at the image's natural resolution now — the
	 * downscaled path that worked around the Cairo 32767-px limit has
	 * been retired by the tile-grid texture wrapping. */
	const int loop_max = (int)height;
	const int sw = (int)width;

	// Shared flag for per-thread allocation failures; checked after the parallel block.
	gboolean alloc_error = FALSE;

	{
		siril_log_debug((com.gui_icc.proofing_transform && !identical && (!com.gui_icc.same_primaries || com.gui_icc.profile_changed)) ? "Non-identical primaries: doing expensive color transform\n" : "");
		const gboolean do_transform = (com.gui_icc.proofing_transform && !identical && (!com.gui_icc.same_primaries || com.gui_icc.profile_changed));

		if (do_transform)
			lock_display_transform();

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) shared(alloc_error)
{
		// Thread-local buffers - allocated once per thread
		WORD *pixelbuf = malloc(width * 3 * sizeof(WORD));
		WORD *linebuf[3] = { pixelbuf, (pixelbuf + width), (pixelbuf + 2 * width) };
		BYTE *pixelbuf_byte = malloc(width * 3);
		BYTE *linebuf_byte[3] = { pixelbuf_byte, (pixelbuf_byte + width), (pixelbuf_byte + 2 * width) };
		BYTE *mask_row = apply_mask ? malloc(width * sizeof(BYTE)) : NULL;

		if (!pixelbuf || !pixelbuf_byte || (apply_mask && !mask_row)) {
			PRINT_ALLOC_ERR;
#pragma omp atomic write
			alloc_error = TRUE;
		} else {
#pragma omp for schedule(static)
#else
		// Single-threaded: allocate once
		WORD *pixelbuf = malloc(width * 3 * sizeof(WORD));
		WORD *linebuf[3] = { pixelbuf, (pixelbuf + width), (pixelbuf + 2 * width) };
		BYTE *pixelbuf_byte = malloc(width * 3);
		BYTE *linebuf_byte[3] = { pixelbuf_byte, (pixelbuf_byte + width), (pixelbuf_byte + 2 * width) };
		BYTE *mask_row = apply_mask ? malloc(width * sizeof(BYTE)) : NULL;

		if (!pixelbuf || !pixelbuf_byte || (apply_mask && !mask_row)) {
			PRINT_ALLOC_ERR;
			free(pixelbuf);
			free(pixelbuf_byte);
			free(mask_row);
			alloc_error = TRUE;
		} else {
#endif
		for (int row = 0; row < loop_max; row++) {
			guint x;
			/* row == destination row; sy == source row (1:1 mapping). */
			const guint sy = (guint)row;
			const guint src_i = sy * width;

			// Precompute mask values for entire source row if mask is active
			if (apply_mask) {
				if (mask_bitpix == 8) {
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = mask_u8[src_i + x];
					}
				} else if (mask_bitpix == 16) {
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = mask_u16[src_i + x] >> 8;
					}
				} else { // 32
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = (BYTE)(mask_f32[src_i + x] * UCHAR_MAX);
					}
				}
			}

			if (gfit->type == DATA_FLOAT) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
					float *source = fsrc[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = roundf_to_WORD(source[src_i + x] * USHRT_MAX_SINGLE);
				}
			} else if (norm == UCHAR_MAX) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
					const WORD *source = src[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = source[src_i + x] << 8;
				}
			} else {
				for (int c = 0 ; c < 3 ; c++)
// No omp simd here as memcpy should already be highly optimized
					memcpy(linebuf[c], src[c] + src_i, width * sizeof(WORD));
			}
			if (gfit->type == DATA_USHORT && norm == UCHAR_MAX) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = line[x] >> 8;
				}
			}
			for (int c = 0 ; c < 3 ; c++) {
				const int cc = gfit->color_managed ? c : 0;
				for (x = 0; x < width; ++x) {
					WORD val = linebuf[c][x];
					if (gui.cut_over && val > remap_hi) {	// cut
						linebuf_byte[c][x] = 0;
					} else {
						linebuf_byte[c][x] = index[cc][val - remap_lo < 0 ? 0 : val - remap_lo];
					}
					if (inverted)
						linebuf_byte[c][x] = UCHAR_MAX - linebuf_byte[c][x];
				}
			}
			if (do_transform) {
				cmsDoTransformLineStride(com.gui_icc.proofing_transform, pixelbuf_byte, pixelbuf_byte, width, 1, width * 3, width * 3, width, width);
			}

			const guint dst_row_start = (height - 1 - sy) * (guint)width * 4;
			const int out_w = sw;  /* number of destination pixels to write */

			if (color == RAINBOW_COLOR) {
				if (apply_mask) {
					// Rainbow with mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[dx];
							const BYTE mask_val = mask_row[dx];

							const BYTE rainbow_r = rainbow_index[pixel_val][0];
							const BYTE rainbow_g = rainbow_index[pixel_val][1];
							const BYTE rainbow_b = rainbow_index[pixel_val][2];

							const uint16_t inv_mask = UCHAR_MAX - mask_val;
							const uint16_t r_sum = inv_mask + rainbow_r;
							const BYTE r_val = (r_sum > UCHAR_MAX) ? UCHAR_MAX : (BYTE)r_sum;
							const BYTE g_val = (rainbow_g * mask_val) >> 8;
							const BYTE b_val = (rainbow_b * mask_val) >> 8;

							*(guint32*)(dst[c] + dst_index) = r_val << 16 | g_val << 8 | b_val;
						}
					}
				} else {
					// Rainbow without mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[dx];
							*(guint32*)(dst[c] + dst_index) = rainbow_index[pixel_val][0] << 16 |
							                                   rainbow_index[pixel_val][1] << 8 |
							                                   rainbow_index[pixel_val][2];
						}
					}
				}
			} else {
				// NORMAL_COLOR
				if (apply_mask) {
					// Normal with mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[dx];
							const BYTE mask_val = mask_row[dx];

							const uint16_t inv_mask = UCHAR_MAX - mask_val;
							const uint16_t r_sum = inv_mask + pixel_val;
							const BYTE r_val = (r_sum > UCHAR_MAX) ? UCHAR_MAX : (BYTE)r_sum;
							const BYTE g_val = (pixel_val * mask_val) >> 8;
							const BYTE b_val = (pixel_val * mask_val) >> 8;

							*(guint32*)(dst[c] + dst_index) = r_val << 16 | g_val << 8 | b_val;
						}
					}
				} else {
					// Normal without mask — SIMD-vectorizable unit-stride
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
#pragma omp simd
						for (x = 0 ; x < (guint)out_w ; x++) {
							const guint dst_index = dst_row_start + x * 4;
							const BYTE pixel_val = line_byte[x];
							*(guint32*)(dst[c] + dst_index) = pixel_val << 16 | pixel_val << 8 | pixel_val;
						}
					}
				}
			}
		}
		// Close the else-branch of the allocation check, then free buffers.
		// In the OpenMP path free() is called by every thread for its own
		// allocations; free(NULL) is safe for the error path where some
		// allocations may have failed.
		}
#ifdef _OPENMP
		free(pixelbuf);
		free(pixelbuf_byte);
		free(mask_row);
}
#else
		free(pixelbuf);
		free(pixelbuf_byte);
		free(mask_row);
#endif
		if (do_transform)
			unlock_display_transform();
	}

	// Bail out cleanly if any thread's allocation failed
	if (alloc_error) {
		g_mutex_unlock(&gui.cairo_mutex);
		return;
	}

	for (int vport = 0 ; vport < 3 ; vport++) {
		view_refresh_tile_textures(view[vport]);
	}
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; both it and
	// test_and_allocate_reference_image must be called after our unlock
	for (int vport = 0 ; vport < 3 ; vport++) {
		invalidate_image_render_cache(vport);
		test_and_allocate_reference_image(vport);
	}
}

static int make_hd_index_for_current_display(int vport) {
	float slope;
	int i;
	BYTE *index;
	float pxl;
	// Check if the bit depth matches the LUT size, if not we need to realloc
	if (gui.hd_remap_max != 1 << com.pref.hd_bitdepth) {
		gui.hd_remap_max = 1 << com.pref.hd_bitdepth;
		allocate_hd_remap_indices();
	}

	/* initialization of data required to build the remap_index
	 *
	 * The HD remap curve is only used with STF rendering mode
	 * as it is the only mode that results in a slope steep enough
	 * to cause noticeable quantization of levels */
	slope = UCHAR_MAX_SINGLE;
	/************* Building the HD remap_index **************/
	siril_log_debug("Rebuilding HD remap_index\n");
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gui.unlink_channels ? vport : 0;
	index = gui.hd_remap_index[target_index];

	for (i = 0; i <= gui.hd_remap_max; i++) {
		pxl = (float) i / gui.hd_remap_max;
		index[i] = roundf_to_BYTE((MTFp(pxl, stf[target_index])) * slope);
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != gui.hd_remap_max + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= gui.hd_remap_max; i++)
			index[i] = UCHAR_MAX;
	}
	return 0;
}

static int make_index_for_current_display(int vport) {
	g_mutex_lock(&com.mutex);
	WORD lo = gui.lo;
	WORD hi = gui.hi;
	g_mutex_unlock(&com.mutex);
	float slope, delta = hi - lo;
	int i;
	BYTE *index;
	float pxl;
	/* initialization of data required to build the remap_index */
	switch (gui.rendering_mode) {
		case LINEAR_DISPLAY:
			slope = UCHAR_MAX_SINGLE / delta;
			break;
		case LOG_DISPLAY:
			slope = fabsf(UCHAR_MAX_SINGLE / logf(delta * 0.1f));
			break;
		case SQRT_DISPLAY:
			slope = UCHAR_MAX_SINGLE / sqrtf(delta);
			break;
		case SQUARED_DISPLAY:
			slope = UCHAR_MAX_SINGLE / SQR(delta);
			break;
		case ASINH_DISPLAY:
			slope = UCHAR_MAX_SINGLE / asinhf(delta * 0.001f);
			break;
		case STF_DISPLAY:
			slope = UCHAR_MAX_SINGLE;
			break;
		default:
			return 1;
	}
	if(!(slope == last_pente && gui.rendering_mode == last_mode))
		com.gui_icc.profile_changed = TRUE;

	if ((gui.rendering_mode != HISTEQ_DISPLAY && gui.rendering_mode != STF_DISPLAY) &&
			slope == last_pente && gui.rendering_mode == last_mode && !com.gui_icc.profile_changed) {
		siril_log_debug("Re-using previous gui.remap_index\n");
		return 0;
	}

	/************* Building the remap_index **************/
	siril_log_debug("Rebuilding gui.remap_index %d\n", vport);
	// target_index only used for STF mode
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.unlink_channels ? vport : 0;
	index = gui.remap_index[vport];
	for (i = 0; i <= USHRT_MAX; i++) {
		switch (gui.rendering_mode) {
			case LOG_DISPLAY:
				// ln(5.56*10^110) = 255
				if (i < 10)
					index[i] = 0; /* avoid null and negative values */
				else
					index[i] = roundf_to_BYTE(logf((float) i / 10.f) * slope); //10.f is arbitrary: good matching with ds9
				break;
			case SQRT_DISPLAY:
				// sqrt(2^16) = 2^8
				index[i] = roundf_to_BYTE(sqrtf((float) i) * slope);
				break;
			case SQUARED_DISPLAY:
				// pow(2^4,2) = 2^8
				index[i] = roundf_to_BYTE(SQR((float)i) * slope);
				break;
			case ASINH_DISPLAY:
				// asinh(2.78*10^110) = 255
				index[i] = roundf_to_BYTE(asinhf((float) i / 1000.f) * slope); //1000.f is arbitrary: good matching with ds9, could be asinhf(a*Q*i)/Q
				break;
			case LINEAR_DISPLAY:
				index[i] = roundf_to_BYTE((float) i * slope);
				break;
			case STF_DISPLAY:
				pxl = (gfit->orig_bitpix == BYTE_IMG ?
						(float) i / UCHAR_MAX_SINGLE :
						(float) i / USHRT_MAX_SINGLE);
				index[i] = roundf_to_BYTE((MTFp(pxl, stf[target_index])) * slope);
				break;
			default:
				return 1;
		}
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++) {
			index[i] = UCHAR_MAX;
		}
	}
	if (gfit->color_managed && com.gui_icc.same_primaries && com.gui_icc.proofing_transform && gui.rendering_mode != STF_DISPLAY)
		display_index_transform(index, vport);

	last_pente = slope;
	last_mode = gui.rendering_mode;
	return 0;
}

static int make_index_for_rainbow(BYTE index[][3]) {
	int i;
	double h, s, v, r, g, b;

	for (i = 0; i < UCHAR_MAX + 1; i++) {
		r = g = b = (double) i / UCHAR_MAX_DOUBLE;
		rgb_to_hsv(r, g, b, &h, &s, &v);
		double off = 300.0 / 360.0; /* Arbitrary: we want h from 300 to 0 deg */
		h = (off - (double) i * (off / UCHAR_MAX_DOUBLE));
		s = 1.;
		v = 1.; /* Saturation and Value are set to 100% */
		hsv_to_rgb(h, s, v, &r, &g, &b);
		index[i][0] = round_to_BYTE(r * UCHAR_MAX_DOUBLE);
		index[i][1] = round_to_BYTE(g * UCHAR_MAX_DOUBLE);
		index[i][2] = round_to_BYTE(b * UCHAR_MAX_DOUBLE);
	}
	return 0;
}

/*****************************************************************************
 * ^ ^ ^ above:     R E M A P P I N G     I M A G E     D A T A        ^ ^ ^ *
 *                                                                           *
 * v v v below:          R E D R A W I N G     W I D G E T S           v v v *
 *****************************************************************************/

typedef struct label_point_struct {
	double x, y, ra, dec, angle;
	gboolean isRA;
	int border;
} label_point;

static void request_gtk_redraw_of_cvport() {
	//siril_log_debug("image redraw requested (vport %d)\n", gui.cvport);
	GtkWidget *widget = gui.view[gui.cvport].drawarea;
	gtk_widget_queue_draw(widget);
}

/* GTK4 has no signal-block facility on the draw_func / snapshot vfunc.
 * Use an atomic flag the renderers consult to short-circuit drawing while
 * a worker repaints the buffers. */
static gint drawarea_handlers_blocked;

/* Forward declarations for overlay-drawing helpers used by both the
 * legacy redraw_drawingarea path and the new SirilImageView::snapshot path
 * below.  All are defined later in this file. */
static void draw_empty_image(const draw_data_t* dd);
static void draw_selection(const draw_data_t* dd);
static void draw_roi(const draw_data_t *dd);
static void draw_cut_line(const draw_data_t* dd);
static void draw_measurement_line(const draw_data_t* dd);
static void draw_in_progress_poly(const draw_data_t* dd);
static void draw_user_polygons(const draw_data_t* dd);
static void draw_stars(const draw_data_t* dd);
static void draw_wcs_grid(const draw_data_t* dd);
static void draw_wcs_disto(const draw_data_t *dd);
static void draw_annotates(const draw_data_t* dd);
static void draw_analysis(const draw_data_t *dd);
static void draw_brg_boxes(const draw_data_t* dd);
static void draw_regframe(const draw_data_t* dd);
static void draw_rgb_centers(const draw_data_t* dd);

/* ── SirilImageView: GTK4 snapshot-based image viewport ───────────────────
 *
 * This widget replaces GtkDrawingArea for the five image viewports.  Its
 * snapshot() vfunc renders the image via gtk_snapshot_append_scaled_texture
 * (from a GdkMemoryTexture that wraps view->buf) and overlays via
 * gtk_snapshot_append_cairo, so the image flows through GSK to the GPU
 * compositor while the existing Cairo overlay drawing code remains
 * unchanged.
 *
 * The disp_surface cache and the Cairo full_surface → disp_surface compose
 * pass that the legacy redraw_drawingarea path uses are bypassed here —
 * GSK already caches its render-node tree between frames, so panning and
 * zooming don't re-walk pixel data.  The legacy Cairo path is retained for
 * add_image_and_label_to_cairo (used by the snapshot-save feature). */

struct _SirilImageView {
	GtkWidget parent_instance;
};

G_DEFINE_FINAL_TYPE(SirilImageView, siril_image_view, GTK_TYPE_WIDGET)

static void siril_image_view_snapshot(GtkWidget *widget, GtkSnapshot *snapshot);

static void siril_image_view_class_init(SirilImageViewClass *klass) {
	GtkWidgetClass *widget_class = GTK_WIDGET_CLASS(klass);
	widget_class->snapshot = siril_image_view_snapshot;
}

static void siril_image_view_init(SirilImageView *self) {
	gtk_widget_set_focusable(GTK_WIDGET(self), TRUE);
}

void siril_image_view_register(void) {
	g_type_ensure(SIRIL_TYPE_IMAGE_VIEW);
}

/* Snapshot vfunc — see redraw_drawingarea for the legacy Cairo equivalent.
 * The two share the overlay-drawing call sequence; only the image draw
 * differs (snapshot append_texture vs Cairo set_source_surface + paint). */
static void siril_image_view_snapshot(GtkWidget *widget, GtkSnapshot *snapshot) {
	if (g_atomic_int_get(&drawarea_handlers_blocked))
		return;

	int vport = match_drawing_area_widget(widget, TRUE);
	if (vport < 0)
		return;

	const int width  = gtk_widget_get_width(widget);
	const int height = gtk_widget_get_height(widget);
	if (width <= 0 || height <= 0)
		return;

	image_display_init_statics();

	g_mutex_lock(&gui.cairo_mutex);
	const gboolean has_image = (gui.view[vport].tiles != NULL);
	g_mutex_unlock(&gui.cairo_mutex);

	/* GTK4 doesn't auto-clip a widget's snapshot tree to its allocation —
	 * append_texture with a transform that extends past the widget bounds
	 * will paint over neighbouring widgets.  Push a clip rectangle that
	 * matches our allocation, covering both the texture pass and the
	 * Cairo overlay below. */
	gtk_snapshot_push_clip(snapshot,
		&GRAPHENE_RECT_INIT(0, 0, (float)width, (float)height));

	if (!has_image) {
		/* No image loaded: route to the existing Cairo "splash" code via
		 * a snapshot Cairo subnode. */
		cairo_t *cr = gtk_snapshot_append_cairo(snapshot,
			&GRAPHENE_RECT_INIT(0, 0, width, height));
		draw_data_t dd = { 0 };
		dd.cr = cr;
		dd.vport = vport;
		dd.window_width = width;
		dd.window_height = height;
		draw_empty_image(&dd);
		cairo_destroy(cr);
		gtk_snapshot_pop(snapshot);
		return;
	}

	/* Recompute display_matrix from current zoom/offset (matches the
	 * legacy path before draw_main_image runs). */
	adjust_vport_size_to_image();

	/* Compute the image-space rectangle currently visible in the widget
	 * by inverting the display matrix on the widget corners.  Tiles
	 * outside this rect are skipped — both to avoid wasted GPU upload
	 * (eager mode) and to avoid eager materialisation (lazy mode). */
	cairo_matrix_t inv = gui.display_matrix;
	cairo_matrix_invert(&inv);
	double cx[4] = { 0, (double)width, (double)width, 0 };
	double cy[4] = { 0, 0, (double)height, (double)height };
	double vis_xmin = +1e30, vis_ymin = +1e30, vis_xmax = -1e30, vis_ymax = -1e30;
	for (int i = 0; i < 4; i++) {
		double x = cx[i], y = cy[i];
		cairo_matrix_transform_point(&inv, &x, &y);
		if (x < vis_xmin) vis_xmin = x;
		if (y < vis_ymin) vis_ymin = y;
		if (x > vis_xmax) vis_xmax = x;
		if (y > vis_ymax) vis_ymax = y;
	}

	/* Walk the tile grid under the mutex.  Materialise dirty/missing
	 * tiles in lazy mode (skipped automatically in eager mode where they
	 * were filled by view_refresh_tile_textures_eager at remap time).
	 * Take refs on every texture we intend to draw so a concurrent remap
	 * can't free them out from under the snapshot. */
	g_mutex_lock(&gui.cairo_mutex);
	struct image_view *view = &gui.view[vport];
	const int tile_dim  = view->tile_dim;
	const int tile_cols = view->tile_cols;
	const int tile_rows = view->tile_rows;
	const int img_w     = view->buf_stride / 4;
	const int img_h     = view->buf_height;
	const gboolean suppressed = g_atomic_int_get(&gui.suppress_drawarea_redraw);

	/* Zoom-aware mip selection: at strong zoom-outs we materialise tiles
	 * at a downsampled resolution so resident bytes track what's actually
	 * visible (see compute_target_mip).  Eager mode (single-tile images)
	 * stays at mip 0 — there's no memory pressure to relieve.
	 *
	 * Once the entire image fits inside the window (zoom ≤ zoom_fit) no
	 * finer mip will reveal more pixels, so cap the effective zoom at
	 * zoom_fit.  This prevents the visible "jerk" of a re-materialise on
	 * every octave of further zoom-out past the fit boundary — at that
	 * point GSK's TRILINEAR pyramid handles the remaining downscale
	 * cheaply on the GPU. */
	const double zoom_now = get_zoom_val();
	double zoom_effective = zoom_now;
	if (img_w > 0 && img_h > 0) {
		const double zoom_fit = MIN((double)width  / img_w,
		                            (double)height / img_h);
		if (zoom_effective < zoom_fit) zoom_effective = zoom_fit;
	}
	const int target_mip = view->lazy ? compute_target_mip(zoom_effective) : 0;

	/* Record the render state the materialise worker can't query itself off
	 * the main thread: the zoom-derived target mip, and the negative-view flag
	 * (GAction lookup is main-thread-only).  Captured under cairo_mutex so the
	 * worker reads a consistent pair. */
	if (view->lazy) {
		view->wk_target_mip = target_mip;
		gboolean neg_now = FALSE;
		GAction *an = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win),
		                                         "negative-view");
		if (an) {
			GVariant *st = g_action_get_state(an);
			neg_now = g_variant_get_boolean(st);
			g_variant_unref(st);
		}
		view->render_neg = neg_now;
	}

	GdkTexture **tiles_snapshot = NULL;
	GdkTexture *proxy_ref = NULL;
	gboolean draw_proxy = FALSE;
	if (view->tiles && tile_cols > 0 && tile_rows > 0) {
		tiles_snapshot = g_new0(GdkTexture *, tile_cols * tile_rows);

		/* Single pass over the tile grid.  We NEVER materialise on the render
		 * thread — a full-resolution tile re-fill is far too slow to do during
		 * a frame on a big image.  Instead we ref and display whatever texture
		 * each visible tile already has (even a stale/dirty one from the
		 * previous LUT, or a too-coarse one from a prior zoom — both beat a
		 * blank or a blurry proxy), record whether anything needs the worker,
		 * and whether any tile is entirely missing (the only case needing the
		 * proxy backdrop). */
		gboolean have_vis = FALSE;
		int vtx0 = 0, vty0 = 0, vtx1 = 0, vty1 = 0;
		gboolean any_work = FALSE;     /* a visible tile is dirty / missing / too coarse */
		gboolean any_missing = FALSE;  /* a visible tile has no texture at all */
		for (int ty = 0; ty < tile_rows; ty++) {
			for (int tx = 0; tx < tile_cols; tx++) {
				const int x0 = tx * tile_dim;
				const int y0 = ty * tile_dim;
				/* Display-space y → image-space y is flipped, but the
				 * visible rect was computed in image-space, so the tile's
				 * stored display position (x0, y0) and visible bounds
				 * compare directly. */
				if (x0 + tile_dim < vis_xmin) continue;
				if (y0 + tile_dim < vis_ymin) continue;
				if (x0 > vis_xmax) continue;
				if (y0 > vis_ymax) continue;
				if (!have_vis) {
					vtx0 = vtx1 = tx; vty0 = vty1 = ty; have_vis = TRUE;
				} else {
					if (tx < vtx0) vtx0 = tx;
					if (tx > vtx1) vtx1 = tx;
					if (ty < vty0) vty0 = ty;
					if (ty > vty1) vty1 = ty;
				}
				struct image_tile *t = &view->tiles[ty * tile_cols + tx];
				if (view->lazy && !suppressed) {
					if (t->dirty || !t->texture
					    || !tile_mip_acceptable(t->mip, target_mip)
					    || tile_mip_too_fine(t->mip, target_mip))
						any_work = TRUE;
					if (!t->texture)
						any_missing = TRUE;
				}
				if (t->texture) {
					tiles_snapshot[ty * tile_cols + tx] = t->texture;
					g_object_ref(t->texture);
					t->last_used = ++view->lazy_epoch;
				}
			}
		}

		/* Hand the visible range to the worker so it only works on-screen. */
		view->vis_tx0 = vtx0; view->vis_ty0 = vty0;
		view->vis_tx1 = vtx1; view->vis_ty1 = vty1;
		view->vis_valid = have_vis;

		/* Proxy backdrop only when a visible tile has NO texture at all (a
		 * never-materialised region just revealed by a zoom/pan).  A stale or
		 * coarse existing texture is displayed directly instead. */
		if (view->lazy && !suppressed && any_missing) {
			ensure_proxy(view, vport);
			if (view->proxy) {
				proxy_ref = g_object_ref(view->proxy);
				draw_proxy = TRUE;
			}
		}

		/* Kick the worker to (re)materialise any visible tile that's dirty,
		 * missing, or too coarse.  No-op if already running, not lazy, or
		 * nothing needs work — so a settled view (e.g. an overlay-only redraw
		 * such as moving the selection) schedules no background work. */
		if (!suppressed)
			schedule_view_worker(view, vport, any_work);
	}
	g_mutex_unlock(&gui.cairo_mutex);

	if (tiles_snapshot && img_w > 0 && img_h > 0) {
		gtk_snapshot_save(snapshot);

		const double xx = gui.display_matrix.xx;
		const double yy = gui.display_matrix.yy;

		/* We deliberately do NOT apply gtk_snapshot_scale(xx, yy) here and
		 * pass image-space rects.  Doing so makes the TextureScaleNode's
		 * rect equal to the mip-0 texture size (a no-op scale), so its
		 * GskScalingFilter is ignored — GSK then upscales to widget-space
		 * via the ambient transform with the renderer's default (linear)
		 * filter, blurring pixel edges at zoom > 1.  Instead we push the
		 * zoom into the per-tile rect dimensions, so the TextureScaleNode
		 * itself does the rasterise-and-scale and NEAREST is honoured. */
		gtk_snapshot_translate(snapshot,
			&GRAPHENE_POINT_INIT((float)gui.display_matrix.x0,
			                      (float)gui.display_matrix.y0));

		/* The remap writes rows bottom-up so that the image renders
		 * top-down on screen — except in the livestacking TOP-DOWN case,
		 * where the buffer is stored top-down and needs a y-flip.  The
		 * flip is now done in widget-space (after the per-tile rect has
		 * had `yy` folded in), so the translate distance is yy*img_h. */
		if (livestacking_is_started()
		    && !g_strcmp0(gfit->keywords.row_order, "TOP-DOWN")) {
			gtk_snapshot_translate(snapshot,
				&GRAPHENE_POINT_INIT(0.0f, (float)(yy * img_h)));
			gtk_snapshot_scale(snapshot, 1.0f, -1.0f);
		}

		const double zoom = get_zoom_val();
		/* TRILINEAR for downscale: mipmap pre-filtering gives correct
		 * antialiasing.  The wider SIRIL_TILE_GUARD masks the per-tile
		 * pyramid discontinuity that LINEAR was a band-aid for —
		 * adjacent tiles' boundary mip-texels are now built from the
		 * same source range at every mipmap level up to log2(GUARD).
		 * NEAREST at zoom ≥ 1 preserves individual-pixel display when
		 * zoomed in (matches the legacy behaviour). */
		const GskScalingFilter filter = (zoom < 1.0)
			? GSK_SCALING_FILTER_TRILINEAR
			: GSK_SCALING_FILTER_NEAREST;

		/* Proxy backdrop: a single scaled-texture node covering the whole
		 * image, drawn *under* the real tiles.  Any visible tile that hasn't
		 * been materialised yet (fast zoom-out) shows the low-res proxy here
		 * instead of a hole; the refine idle then paints sharp tiles over it
		 * frame by frame.  Same widget-space rect convention as the tiles. */
		if (draw_proxy && proxy_ref) {
			gtk_snapshot_append_scaled_texture(snapshot, proxy_ref, filter,
				&GRAPHENE_RECT_INIT(0.0f, 0.0f,
					(float)(img_w * xx), (float)(img_h * yy)));
		}

		for (int ty = 0; ty < tile_rows; ty++) {
			for (int tx = 0; tx < tile_cols; tx++) {
				GdkTexture *tile = tiles_snapshot[ty * tile_cols + tx];
				if (!tile) continue;
				/* Rect is in WIDGET-SPACE (post-translate): image-space
				 * coords multiplied by the per-axis zoom.  This makes
				 * rect_size != texture_size at any zoom ≠ 1, so GSK's
				 * TextureScaleNode rasteriser is what scales the texture
				 * — and `filter` actually applies. */
				int x0, y0, tw_unused, th_unused, tex_w_img, tex_h_img;
				tile_dims_padded(view, tx, ty, &x0, &y0, &tw_unused,
					&th_unused, &tex_w_img, &tex_h_img);
				gtk_snapshot_append_scaled_texture(snapshot, tile, filter,
					&GRAPHENE_RECT_INIT(
						(float)(x0 * xx),
						(float)(y0 * yy),
						(float)(tex_w_img * xx),
						(float)(tex_h_img * yy)));
			}
		}

		gtk_snapshot_restore(snapshot);
	}
	if (tiles_snapshot) {
		for (int i = 0; i < tile_cols * tile_rows; i++) {
			if (tiles_snapshot[i]) g_object_unref(tiles_snapshot[i]);
		}
		g_free(tiles_snapshot);
	}
	if (proxy_ref)
		g_object_unref(proxy_ref);

	/* Overlays: open a Cairo subnode covering the full widget, apply
	 * display_matrix so overlay drawing functions operate in image-space
	 * (matches what redraw_drawingarea sets up). */
	cairo_t *cr = gtk_snapshot_append_cairo(snapshot,
		&GRAPHENE_RECT_INIT(0, 0, width, height));
	cairo_transform(cr, &gui.display_matrix);

	draw_data_t dd = { 0 };
	dd.cr = cr;
	dd.vport = vport;
	dd.window_width = width;
	dd.window_height = height;
	dd.zoom = get_zoom_val();
	dd.image_width = gfit->rx;
	dd.image_height = gfit->ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_NEAREST;

	static GAction *action_neg = NULL;
	if (action_neg == NULL)
		action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
	if (action_neg) {
		GVariant *state = g_action_get_state(action_neg);
		dd.neg_view = g_variant_get_boolean(state);
		g_variant_unref(state);
	}

	/* Same overlay sequence redraw_drawingarea uses, minus draw_main_image
	 * (which we've already rendered via the texture path above). */
	draw_selection(&dd);
	draw_roi(&dd);
	draw_cut_line(&dd);
	draw_measurement_line(&dd);
	draw_in_progress_poly(&dd);
	draw_user_polygons(&dd);
	g_mutex_lock(&com.mutex);
	draw_stars(&dd);
	g_mutex_unlock(&com.mutex);
	draw_wcs_grid(&dd);
	draw_wcs_disto(&dd);
	draw_annotates(&dd);
	draw_analysis(&dd);
	draw_brg_boxes(&dd);
	draw_regframe(&dd);
	draw_rgb_centers(&dd);
	if (gui.draw_extra)
		gui.draw_extra(&dd);

	/* GTK3-equivalent: histogram_overlay_paint inherits the
	 * display_matrix transform that's been applied to cr above — this
	 * matches the legacy redraw_drawingarea behaviour and is what the
	 * widget-coords-vs-image-coords positioning expects. */
	histogram_overlay_paint(widget, cr);

	cairo_destroy(cr);

	gtk_snapshot_pop(snapshot);  /* matches push_clip at the top */
}


void block_drawarea_handlers(void) {
	g_atomic_int_set(&drawarea_handlers_blocked, 1);
}

void unblock_drawarea_handlers(void) {
	g_atomic_int_set(&drawarea_handlers_blocked, 0);
}

void install_drawarea_draw_funcs(void) {
	/* No-op in Phase 2: SirilImageView's snapshot vfunc is set at class
	 * init time, not per-instance.  Kept as an empty entry point so
	 * callbacks.c::initialize_all_GUI can keep calling it without a
	 * compile-time conditional. */
}

static void draw_empty_image(const draw_data_t* dd) {
	cairo_t *cr = dd->cr;
	guint width = dd->window_width;
	guint height = dd->window_height;
	guint pix_size = height / 3;
	guint offset;
#ifdef SIRIL_UNSTABLE
	offset = 32;
#else
	offset = 2;
#endif

	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	/* Render siril.svg via GTK4's paintable pipeline (works without the
	 * librsvg GdkPixbuf loader, which isn't installed on every system). */
	siril_cairo_paint_resource(cr, "/org/siril/ui/pixmaps/siril.svg",
	                           (width - pix_size) / 2,
	                           (height - pix_size) / offset,
	                           pix_size, pix_size);


	image_display_init_statics();
	/* Use the drawing area we're actually painting on (via dd's
	 * window_width/window_height) — the previous imgdisp_drawing_rgb
	 * lookup returned the RGB-tab drawingarea, which has no allocation
	 * when the user is on a different viewport tab, leading to
	 * scale=0 and invisible text. */
	GtkWidget *widget = imgdisp_drawing_rgb;
	(void) gtk_widget_get_state_flags(widget);
	PangoLayout *layout;
	gchar *msg;
	gdouble scale;
	GdkRGBA color;
	gint w, h;
	const gdouble alloc_w = (gdouble) width;
	const gdouble alloc_h = (gdouble) height;

	layout = gtk_widget_create_pango_layout(widget, NULL);

#ifdef SIRIL_UNSTABLE

	msg = g_strdup_printf(_("<big>Unstable Development Version</big>\n\n"
			    "<small>%c%s</small>\n"
				"<small>commit <tt>%s</tt></small>\n"
				"<small>Please test bugs against "
				"latest git master branch\n"
				"before reporting them.</small>"),
			toupper(PACKAGE_STRING[0]),
			(char *)PACKAGE_STRING + 1,
			SIRIL_GIT_VERSION_ABBREV);
	pango_layout_set_markup(layout, msg, -1);
	g_free(msg);
	pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

	pango_layout_get_pixel_size(layout, &w, &h);

	scale = MIN((alloc_w / 2.0) / (gdouble) w,
	            (alloc_h / 2.0) / (gdouble) h / 2.0);
	if (scale <= 0.0) scale = 1.0;

	gtk_widget_get_color(widget, &color);
	gdk_cairo_set_source_rgba(cr, &color);

	cairo_move_to(cr, (alloc_w - (w * scale)) / 2,
	                  (alloc_h - (h * scale)) / 2);

#else
	msg = g_strdup_printf("%c%s", toupper(PACKAGE_STRING[0]), (char *)PACKAGE_STRING + 1);

	pango_layout_set_markup(layout, msg, -1);
	g_free(msg);
	pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

	pango_layout_get_pixel_size(layout, &w, &h);

	scale = MIN((alloc_w / 4.0) / (gdouble) w,
	            (alloc_h / 4.0) / (gdouble) h / 4.0);
	if (scale <= 0.0) scale = 1.0;

	gtk_widget_get_color(widget, &color);
	gdk_cairo_set_source_rgba(cr, &color);

	cairo_move_to(cr, (alloc_w - (w * scale)) / 2,
	                  3 * (alloc_h - (h * scale)) / 4);

#endif /* SIRIL_UNSTABLE */
	cairo_scale(cr, scale, scale);

	pango_cairo_show_layout(cr, layout);

	g_object_unref(layout);
}

/* Compose the viewport image into dd->cr by painting each tile of
 * view->buf as a per-tile Cairo image surface (wrapping a sub-region of
 * the buffer via the buf row stride).  Used by add_image_and_label_to_cairo
 * for the snapshot-save / copy-to-clipboard path.  On-screen drawing
 * goes through SirilImageView::snapshot instead.
 *
 * After this returns the cr's transform has display_matrix applied so
 * subsequent overlay drawing operates in image-space — matching the
 * legacy convention. */
static void draw_vport(const draw_data_t* dd) {
	g_mutex_lock(&gui.cairo_mutex);
	struct image_view *view = &gui.view[dd->vport];
	const int img_w = view->buf_stride > 0 ? view->buf_stride / 4 : 0;
	const int img_h = view->buf_height;
	const int tile_cols = view->tile_cols;
	const int tile_rows = view->tile_rows;
	const int buf_stride = view->buf_stride;

	if (img_w <= 0 || img_h <= 0 || !view->tiles) {
		g_mutex_unlock(&gui.cairo_mutex);
		return;
	}

	cairo_save(dd->cr);
	cairo_transform(dd->cr, &gui.display_matrix);

	if (livestacking_is_started()
	    && !g_strcmp0(gfit->keywords.row_order, "TOP-DOWN")) {
		cairo_translate(dd->cr, 0.0, (double)img_h);
		cairo_scale(dd->cr, 1.0, -1.0);
	}

	/* Walk the tile grid, painting each as a separate Cairo image
	 * surface.  In eager mode the surface wraps a sub-region of
	 * view->buf via the full row stride.  In lazy mode the tile is
	 * materialised on demand, the surface wraps its compact data block
	 * with a per-tile stride.  Painting + destroying the surface
	 * synchronously means cairo never holds on to view bytes past the
	 * fill, so the LRU eviction policy is free to reclaim them. */
	for (int ty = 0; ty < tile_rows; ty++) {
		for (int tx = 0; tx < tile_cols; tx++) {
			/* Use the PADDED dims (tex_w × tex_h, which include up to
			 * SIRIL_TILE_GUARD pixels of overlap with the right /
			 * bottom neighbour) as both the source-surface size and
			 * the destination rectangle.  Adjacent tiles then overlap
			 * by `guard` pixels, and the filter at each tile's edge
			 * has real data to sample instead of transparent black
			 * past the surface bounds — the on-screen path uses the
			 * same trick at snapshot_image_view():2129–2137.  The
			 * overlap pixels are duplicates of the next tile's data so
			 * later-drawn tiles overpaint earlier ones with identical
			 * values; visually a no-op, but the filter sees continuous
			 * source data and no grid lines appear at boundaries. */
			int x0, y0, tw_unused, th_unused, tex_w, tex_h;
			tile_dims_padded(view, tx, ty, &x0, &y0,
			                 &tw_unused, &th_unused, &tex_w, &tex_h);
			(void) tw_unused; (void) th_unused;
			if (tex_w <= 0 || tex_h <= 0) continue;

			guchar *p = NULL;
			int stride;
			if (view->lazy) {
				/* The Cairo fallback consumes tile->data as a native-stride
				 * image surface; it can't accept a mip-downsampled buffer.
				 * Force mip 0 here.  This path is only used for non-snapshot
				 * destinations (export, save-as) so the occasional re-
				 * materialise vs the snapshot mip is not on the hot path. */
				if (!materialise_tile(view, dd->vport, tx, ty, 0)) continue;
				p = view->tiles[ty * tile_cols + tx].data;
				/* Lazy tiles are laid out at the padded width. */
				stride = tex_w * 4;
			} else {
				/* Eager mode shares a single buffer for the whole
				 * image; the guard region is naturally available
				 * because it's just the next tile's leftmost pixels
				 * in the same buffer. */
				p = view->buf + (size_t)y0 * buf_stride + (size_t)x0 * 4;
				stride = buf_stride;
			}
			if (!p) continue;

			cairo_surface_t *tile_surf = cairo_image_surface_create_for_data(
				p, CAIRO_FORMAT_RGB24, tex_w, tex_h, stride);
			if (cairo_surface_status(tile_surf) != CAIRO_STATUS_SUCCESS) {
				cairo_surface_destroy(tile_surf);
				continue;
			}
			cairo_set_source_surface(dd->cr, tile_surf, x0, y0);
			cairo_pattern_set_filter(cairo_get_source(dd->cr), dd->filter);
			cairo_rectangle(dd->cr, x0, y0, tex_w, tex_h);
			cairo_fill(dd->cr);
			cairo_surface_destroy(tile_surf);
		}
	}

	cairo_restore(dd->cr);
	/* Re-apply display_matrix to dd->cr (cairo_restore reset it) so
	 * subsequent overlay drawing operates in image-space, matching the
	 * legacy contract. */
	cairo_transform(dd->cr, &gui.display_matrix);
	g_mutex_unlock(&gui.cairo_mutex);
}

static void draw_main_image(const draw_data_t* dd) {
	g_mutex_lock(&gui.cairo_mutex);
	gboolean has_image = (gui.view[dd->vport].tiles != NULL);
	g_mutex_unlock(&gui.cairo_mutex);

	if (has_image) {
		draw_vport(dd);
	} else {
		draw_empty_image(dd);
	}
}

gboolean get_context_rotation_matrix(double rotation, cairo_matrix_t *transform, gboolean invert) {
	if (rotation == 0.) return FALSE;
	double dx = (double)com.selection.x + (double)com.selection.w * 0.5;
	double dy = (double)com.selection.y + (double)com.selection.h * 0.5;
	cairo_matrix_init_translate(transform, dx, dy);
	cairo_matrix_rotate(transform, rotation * DEGTORAD);
	cairo_matrix_translate(transform, -dx, -dy);
	if (invert) return (cairo_matrix_invert(transform) == CAIRO_STATUS_SUCCESS);
	return TRUE;
}

static void rotate_context(cairo_t *cr, double rotation) {
	cairo_matrix_t transform;
	if (!get_context_rotation_matrix(rotation, &transform, FALSE)) return;
	cairo_transform(cr, &transform);
}

static void draw_roi(const draw_data_t *dd) {
	double r, g, b;
	if (gui.roi.operation_supports_roi) {
		r = 0.3; g = 1.0; b = 0.3;
	} else {
		r = 1.0; g = 0.0; b = 0.0;
	}
	if (gui.roi.selection.w > 0 && gui.roi.selection.h > 0 && gui.roi.active) {
		cairo_t *cr = dd->cr;
		static double dash_format[] = { 4.0, 2.0 };
		cairo_set_line_width(cr, 1.5 / dd->zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, r, g, b);
		cairo_save(cr); // save the original transform
		cairo_rectangle(cr, (double) gui.roi.selection.x, (double) gui.roi.selection.y,
						(double) gui.roi.selection.w, (double) gui.roi.selection.h);
		cairo_stroke(cr);
		cairo_restore(cr);
	}
}

static void draw_selection(const draw_data_t* dd) {
	if (com.selection.w > 0 && com.selection.h > 0) {
		if ((com.selection.x + com.selection.w > gfit->rx) ||
		(com.selection.y + com.selection.h > gfit->ry)) {
			rectangle area = {0, 0, gfit->rx, gfit->ry};
			memcpy(&com.selection, &area, sizeof(rectangle));
		}
		if (!rotation_dlg) image_display_init_statics();
		cairo_t *cr = dd->cr;
		static double dash_format[] = { 4.0, 2.0 };
		cairo_set_line_width(cr, 1.5 / dd->zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
		cairo_save(cr); // save the original transform
		if (gtk_widget_is_visible(rotation_dlg)) {
			double dashes2[]={5.0, 5.0};
			cairo_set_dash(cr, dashes2, 2, 0);
			cairo_set_line_width(cr, 0.5 / dd->zoom);
			cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
						(double) com.selection.w, (double) com.selection.h);
			cairo_stroke(cr);
			cairo_set_line_width(cr, 3. / dd->zoom);
			cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);
			rotate_context(cr, -gui.rotation); // cairo is positive CW while opencv is positive CCW

			// draw a circle at top left corner to visualize rots larger than 90
			double size = 10. / dd->zoom;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_arc(cr, com.selection.x, com.selection.y, size * 0.5, 0., 2. * M_PI);
			cairo_stroke_preserve(cr);
			cairo_fill(cr);
			cairo_set_dash(cr, dash_format, 2, 0);
		}
		cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
						(double) com.selection.w, (double) com.selection.h);
		cairo_stroke(cr);

		// display a grid when the selection is being made / modified, when it is big enough
		if (com.pref.gui.selection_guides > 1 && gui.drawing && com.selection.w > 40 / dd->zoom && com.selection.h > 40 / dd->zoom) {
			cairo_set_line_width(cr, 0.4 / dd->zoom);
			cairo_set_dash(cr, NULL, 0, 0);
			for (int i = 1; i < com.pref.gui.selection_guides; i++) {
				int x = com.selection.x + com.selection.w * i / com.pref.gui.selection_guides;
				int y = com.selection.y + com.selection.h * i / com.pref.gui.selection_guides;
				cairo_move_to(cr, x, com.selection.y);
				cairo_line_to(cr, x, com.selection.y + com.selection.h);
				cairo_move_to(cr, com.selection.x, y);
				cairo_line_to(cr, com.selection.x + com.selection.w, y);
			}
			cairo_stroke(cr);
		}

		// display a mini cross when the selection is being dragged
		if ((gui.freezeX && gui.freezeY) || gtk_widget_is_visible(rotation_dlg)) {
			cairo_set_line_width(cr, 1.0 / dd->zoom);
			point selection_center = { com.selection.x + (double)com.selection.w / 2.,
				com.selection.y + (double)com.selection.h / 2. };
			cairo_move_to(cr, selection_center.x, selection_center.y - 5 / dd->zoom);
			cairo_line_to(cr, selection_center.x, selection_center.y + 5 / dd->zoom);
			cairo_move_to(cr, selection_center.x - 5 / dd->zoom, selection_center.y);
			cairo_line_to(cr, selection_center.x + 5 / dd->zoom, selection_center.y);
			cairo_stroke(cr);
		}
		cairo_restore(cr); // restore the original transform
	}
}

static void draw_cut_line(const draw_data_t* dd) {
	if (!cut_dialog)
		image_display_init_statics();
	if (!(gtk_widget_get_visible(cut_dialog) || gtk_widget_get_visible(cut_cdialog) || gtk_widget_get_visible(cut_sdialog)))
		return;
	if (gui.cut.cut_end.x == -1 || gui.cut.cut_end.y == -1 || gui.cut.seq)
		return;
	gboolean tri = siril_toggle_get_active(GTK_WIDGET(tri_cut_toggle));
	double offstartx, offstarty, offendx, offendy, step;

	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	static double solid_format[] = { 1.0, 0.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);

	if (tri) {
		point delta;
		delta.x = gui.cut.cut_end.x - gui.cut.cut_start.x;
		delta.y = gui.cut.cut_end.y - gui.cut.cut_start.y;
		double length = sqrt(delta.x * delta.x + delta.y * delta.y);
		if (length < 1.) return;
		int nbr_points = (int) length;
		double point_spacing_x = delta.x / nbr_points;
		double point_spacing_y = delta.y / nbr_points;
		step = gtk_spin_button_get_value(tri_cut_spin_step);
		double line_r[3] = { 0.58, 0.0, 0.34 }; // These colours match the 3 lines plotted by siril plot
		double line_g[3] = { 0.0, 0.62, 0.70 };
		double line_b[3] = { 0.83, 0.45, 0.91 };
		double arrow_length = 10 / dd->zoom;
		double arrow_angle = 0.5;
		double angle = atan2(gui.cut.cut_end.y - gui.cut.cut_start.y, gui.cut.cut_end.x - gui.cut.cut_start.x);
		for (int offset = -1 ; offset < 2 ; offset++) {
			cairo_set_dash(cr, dash_format, 2, 0);
			offstartx = gui.cut.cut_start.x + (offset * point_spacing_y * step);
			offstarty = gui.cut.cut_start.y - (offset * point_spacing_x * step);
			offendx = gui.cut.cut_end.x + (offset * point_spacing_y * step);
			offendy = gui.cut.cut_end.y - (offset * point_spacing_x * step);
			cairo_set_source_rgb(cr, line_r[offset+1], line_g[offset+1], line_b[offset+1]);
			cairo_save(cr);
			cairo_move_to(cr, offstartx + 0.5, offstarty + 0.5);
			cairo_line_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_stroke(cr);
			// Draw arrowheads at the end
			cairo_set_dash(cr, solid_format, 0, 0); // Draw the arrow heads solid
			point pt1 = { offendx + 0.5 - arrow_length * cos(angle - arrow_angle), offendy + 0.5 - arrow_length * sin(angle - arrow_angle) };
			point pt2 = { offendx + 0.5 - arrow_length * cos(angle + arrow_angle), offendy + 0.5 - arrow_length * sin(angle + arrow_angle) };
			cairo_line_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_line_to(cr, pt1.x, pt1.y);
			cairo_move_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_line_to(cr, pt2.x, pt2.y);
			cairo_stroke(cr);
			cairo_restore(cr);
		}
	} else {
		cairo_set_source_rgb(cr, 0.0, 0.62, 0.70); // This matches the single line plotted by siril plot
		cairo_save(cr);
		cairo_move_to(cr, gui.cut.cut_start.x + 0.5, gui.cut.cut_start.y + 0.5);
		cairo_line_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		// Draw an arrowhead at the end
		double arrow_length = 10 / dd->zoom;
		double arrow_angle = 0.5;
		double angle = atan2(gui.cut.cut_end.y - gui.cut.cut_start.y, gui.cut.cut_end.x - gui.cut.cut_start.x);
		cairo_stroke(cr);
		cairo_set_dash(cr, solid_format, 0, 0); // Draw the arrow heads solid
		point pt1 = { gui.cut.cut_end.x + 0.5 - arrow_length * cos(angle - arrow_angle), gui.cut.cut_end.y + 0.5 - arrow_length * sin(angle - arrow_angle) };
		point pt2 = { gui.cut.cut_end.x + 0.5 - arrow_length * cos(angle + arrow_angle), gui.cut.cut_end.y + 0.5 - arrow_length * sin(angle + arrow_angle) };
		cairo_line_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		cairo_line_to(cr, pt1.x, pt1.y);
		cairo_move_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		cairo_line_to(cr, pt2.x, pt2.y);
		cairo_stroke(cr);
		cairo_restore(cr);
	}
}

static void draw_measurement_line(const draw_data_t* dd) {
	if (gui.measure_start.x == -1)
		return;
	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
	cairo_save(cr);
	cairo_move_to(cr, gui.measure_start.x + 0.5, gui.measure_start.y + 0.5);
	cairo_line_to(cr, gui.measure_end.x + 0.5, gui.measure_end.y + 0.5);
	cairo_stroke(cr);
	cairo_restore(cr);
}

static void draw_user_polygons(const draw_data_t *dd) {
	if (gui.user_polygons == NULL)
		return;
	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);
	GSList *l;
	for (l = gui.user_polygons; l != NULL; l = l->next) {
		UserPolygon *polygon = (UserPolygon *)l->data;
		if (polygon->n_points < 2)
			continue;

		// Calculate center of polygon for legend placement
		double center_x = 0, center_y = 0;
		for (int i = 0; i < polygon->n_points; i++) {
			center_x += polygon->points[i].x;
			center_y += polygon->points[i].y;
		}
		center_x /= polygon->n_points;
		center_y /= polygon->n_points;

		if (polygon->fill) {
			cairo_save(cr);
			// Set color for filling
			cairo_set_source_rgba(cr,
								  polygon->color[0],
						 polygon->color[1],
						 polygon->color[2],
						 polygon->color[3]);
			cairo_move_to(cr, polygon->points[0].x + 0.5, polygon->points[0].y + 0.5);
			for (int i = 1; i < polygon->n_points; i++) {
				cairo_line_to(cr, polygon->points[i].x + 0.5, polygon->points[i].y + 0.5);
			}
			cairo_close_path(cr);
			cairo_fill_preserve(cr);
			/* slightly darker outline when filled */
			cairo_set_source_rgba(cr,
								  polygon->color[0] * 0.8,
						 polygon->color[1] * 0.8,
						 polygon->color[2] * 0.8,
						 polygon->color[3]);
			cairo_stroke(cr);
			cairo_restore(cr);
		} else {
			cairo_save(cr);
			cairo_set_source_rgba(cr,
								  polygon->color[0],
						 polygon->color[1],
						 polygon->color[2],
						 polygon->color[3]);
			cairo_move_to(cr, polygon->points[0].x + 0.5, polygon->points[0].y + 0.5);
			for (int i = 1; i < polygon->n_points; i++) {
				cairo_line_to(cr, polygon->points[i].x + 0.5, polygon->points[i].y + 0.5);
			}
			cairo_close_path(cr);
			cairo_stroke(cr);
			cairo_restore(cr);
		}

		// Draw legend text if it exists
		if (polygon->legend != NULL && polygon->legend[0] != '\0') {
			cairo_save(cr);

			// Set up text properties
			cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, 12.0 / dd->zoom);

			// Get text dimensions to center properly
			cairo_text_extents_t text_extents;
			cairo_text_extents(cr, polygon->legend, &text_extents);

			// Create a background for the text for better visibility
			double padding = 4.0 / dd->zoom;
			cairo_set_source_rgba(cr, 0, 0, 0, 0.5);  // Darken the background
			cairo_rectangle(cr,
							center_x - text_extents.width/2 - padding,
				   center_y - text_extents.height/2 - padding,
				   text_extents.width + 2*padding,
				   text_extents.height + 2*padding);
			cairo_fill(cr);

			// Draw the text
			cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);  // White text
			cairo_move_to(cr,
						  center_x - text_extents.width/2,
				 center_y + text_extents.height/2);
			cairo_show_text(cr, polygon->legend);
			cairo_restore(cr);
		}
	}
}

static void draw_stars(const draw_data_t* dd) {
	cairo_t *cr = dd->cr;
	int i = 0;

	if (!com.script && (single_image_is_loaded() || sequence_is_loaded()) &&
			g_rw_lock_reader_trylock(&com.stars_lock)) {
		/* com.stars is a NULL-terminated array */
		if (com.stars) {
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / dd->zoom);

		while (com.stars[i]) {
			double size = com.stars[i]->fwhmx * 2.0;
			if (size <= 0.0) size = com.pref.phot_set.aperture;
			if (i == gui.selected_star) {
				// We draw horizontal and vertical lines to show the star
				cairo_set_line_width(cr, 2.0 / dd->zoom);
				cairo_set_source_rgba(cr, 0.0, 0.4, 1.0, 0.6);

				cairo_move_to(cr, com.stars[i]->xpos, 0);
				cairo_line_to(cr, com.stars[i]->xpos, dd->image_height);
				cairo_stroke(cr);
				cairo_move_to(cr, 0, com.stars[i]->ypos);
				cairo_line_to(cr, dd->image_width, com.stars[i]->ypos);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / dd->zoom);
			}
			if (com.stars[i]->has_saturated) {
				cairo_set_source_rgba(cr, 0.75, 0.22, 1.0, 0.9);
				cairo_set_line_width(cr, 3.0 / dd->zoom);
			}
			cairo_save(cr); // save the original transform
			cairo_translate(cr, com.stars[i]->xpos, com.stars[i]->ypos);
			cairo_rotate(cr, M_PI * 0.5 + com.stars[i]->angle * M_PI / 180.);
			double r = com.stars[i]->fwhmx > 0.0 ? com.stars[i]->fwhmy / com.stars[i]->fwhmx : 1.0;
			cairo_scale(cr, r, 1);
			cairo_arc(cr, 0., 0., size, 0., 2 * M_PI);
			cairo_restore(cr); // restore the original transform
			cairo_stroke(cr);
			/* to keep  for debugging boxes adjustements */
			// if (com.stars[i]->R > 0)
			// 	cairo_rectangle(cr, com.stars[i]->xpos - (double)com.stars[i]->R, com.stars[i]->ypos - (double)com.stars[i]->R, (double)com.stars[i]->R * 2 + 1, (double)com.stars[i]->R * 2 + 1);
			// cairo_stroke(cr);
			if (com.stars[i]->has_saturated) {
				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / dd->zoom);
			}

			i++;
		}
		} // if (com.stars)
		g_rw_lock_reader_unlock(&com.stars_lock);
	}

	/* quick photometry */
	if (!com.script && gui.qphot && mouse_status == MOUSE_ACTION_PHOTOMETRY) {
		double size = (!com.pref.phot_set.force_radius) ? 0.5 * gui.qphot->fwhmx * com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
		if (size <= 0.0) size = com.pref.phot_set.aperture;

		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / dd->zoom);

		/* fwhm * 2: first circle */
		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, size, 0., 2. * M_PI);
		cairo_stroke(cr);

		/* sky annulus */
		if (dd->neg_view) {
			cairo_set_source_rgba(cr, 0.5, 0.0, 0.7, 0.9);
		} else {
			cairo_set_source_rgba(cr, 0.5, 1.0, 0.3, 0.9);
		}

		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, com.pref.phot_set.inner, 0., 2. * M_PI);
		cairo_stroke(cr);
		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, com.pref.phot_set.outer, 0., 2. * M_PI);
		cairo_stroke(cr);
		cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		cairo_set_font_size(cr, 40.0 / dd->zoom);
		cairo_move_to(cr, gui.qphot->xpos + com.pref.phot_set.outer + 5, gui.qphot->ypos);
		cairo_show_text(cr, "V");  // was missing — stroke on empty path was a no-op
		cairo_stroke(cr);
	}

	/* draw seqpsf stars */
	if (sequence_is_loaded() && com.seq.current >= 0) {
		for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
			psf_star *the_psf = com.seq.photometry[i][com.seq.current];
			if (the_psf) {
				double size = (!com.pref.phot_set.force_radius && the_psf->fwhmx > 0.0) ?
					0.5 * the_psf->fwhmx * com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
				cairo_set_dash(cr, NULL, 0, 0);
				// make the aperture slightly brighter
				cairo_set_source_rgba(cr, min(com.seq.photometry_colors[i][0] + 0.2, 1.0),
						min(com.seq.photometry_colors[i][1] + 0.2, 1.0),
						min(com.seq.photometry_colors[i][2] + 0.2, 1.0), 1.0);
				cairo_set_line_width(cr, 2.0 / dd->zoom);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, size, 0., 2. * M_PI);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, com.seq.photometry_colors[i][0],
						com.seq.photometry_colors[i][1],
						com.seq.photometry_colors[i][2], 1.0);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.inner, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.outer, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
				cairo_set_font_size(cr, 40.0 / dd->zoom);
				cairo_move_to(cr, the_psf->xpos + com.pref.phot_set.outer + 5, the_psf->ypos);
				if (i == 0) {
					cairo_show_text(cr, "V");
				} else {
					char tmp[16];
					sprintf(tmp, "%d", i);
					cairo_show_text(cr, tmp);
				}
				cairo_stroke(cr);
			}
		}

		/* draw a cross on excluded images */
		if (com.seq.imgparam && com.seq.current >= 0 &&
				!com.seq.imgparam[com.seq.current].incl) {
			int w = dd->image_width > gfit->rx ? gfit->rx : dd->image_width;
			int h = dd->image_height > gfit->ry ? gfit->ry : dd->image_height;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
			cairo_set_line_width(cr, 2.0 / dd->zoom);
			cairo_move_to(cr, 0, 0);
			cairo_line_to(cr, w, h);
			cairo_move_to(cr, 0, h);
			cairo_line_to(cr, w, 0.0);
			cairo_stroke(cr);
		}

		/* draw preview rectangles for the manual registration */
		for (i = 0; i < PREVIEW_NB; i++) {
			if (com.seq.previewX[i] >= 0) {
				int textX, textY;
				gchar *text;
				cairo_set_line_width(cr, 1.0 / dd->zoom);
				cairo_set_source_rgb(cr, 0.1, 0.6, 0.0);
				cairo_rectangle(cr,
						com.seq.previewX[i] - com.seq.previewW[i] / 2,
						com.seq.previewY[i] - com.seq.previewH[i] / 2,
						com.seq.previewW[i], com.seq.previewH[i]);
				cairo_stroke(cr);

				textX = com.seq.previewX[i] - com.seq.previewW[i] / 2;
				textX += 0.1 * com.seq.previewW[i];

				textY = com.seq.previewY[i] - com.seq.previewH[i] / 2;
				textY += 0.1 * com.seq.previewH[i];

				text = g_strdup_printf("%d", i + 1);

				cairo_set_font_size(cr, 12.0 / dd->zoom);
				cairo_move_to(cr, textX, textY);
				cairo_show_text(cr, text);
				cairo_stroke(cr);
				g_free(text);
			}
		}
	}
}

static void draw_brg_boxes(const draw_data_t* dd) {
	GSList *list;
	GdkRGBA gdk_color;
	for (list = com.grad_samples; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		if (sample && background_sample_is_valid(sample)) {
			int radius = (int) (background_sample_get_size(sample) / 2);
			point position = background_sample_get_position(sample);
			cairo_set_line_width(dd->cr, 1.5 / dd->zoom);
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_bkg_samples);
			cairo_set_source_rgba(dd->cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			cairo_rectangle(dd->cr, position.x - radius - 1, position.y - radius,
					radius * 2, radius * 2);
			cairo_stroke(dd->cr);
		}
	}
}

static void draw_in_progress_poly(const draw_data_t* dd) {
	if (!gui.drawing_polypoints) return;

	cairo_t *cr = dd->cr;
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	gdk_cairo_set_source_rgba(cr, &gui.poly_ink);

	GSList *list = gui.drawing_polypoints;
	const point* start = (point*) list->data;
	cairo_move_to(cr, start->x + 0.5, start->y + 0.5);

	// Build the complete path first
	for (list = list->next; list; list = list->next) {
		const point* position = (point*) list->data;
		cairo_line_to(cr, position->x + 0.5, position->y + 0.5);
	}

	// Now stroke the entire path at once
	cairo_stroke(cr);
}

static void draw_compass(const draw_data_t* dd) {
	int pos = com.pref.gui.position_compass;
	if (!pos) return; // User chose None
	fits *fit = gfit;
	cairo_t *cr = dd->cr;
	cairo_set_line_width(cr, 3.0 / dd->zoom);
	double ra0, dec0;
	double xN, yN, xE, yE;

	double xpos = fit->rx * 0.5;
	double ypos = fit->ry * 0.5;
	pix2wcs(fit, xpos, ypos, &ra0, &dec0);
	if (ra0 == -1) return; // checks implicitly that wcslib member exists
	if (90. - fabs(dec0) < 2.78e-3) return;// center is less than 10"off from a pole
	double len = (double)fit->ry / 20.;
	wcs2pix(fit, ra0, dec0 + 0.1, &xN, &yN);
	wcs2pix(fit, ra0 + 0.1, dec0, &xE, &yE);
	double angleN = -atan2(yN - ypos, xN - xpos);
	double angleE = -atan2(yE - ypos, xE - xpos);

	cairo_set_font_size(cr, len / 3);

	double xdraw, ydraw;
	double pos_values[5][2] = { { 0.5, 0.5 }, { 0.1, 0.1 }, { 0.9, 0.1 }, { 0.1, 0.9 }, { 0.9, 0.9 } };
	xdraw = pos_values[pos - 1][0] * fit->rx;
	ydraw = pos_values[pos - 1][1] * fit->ry;

	/* draw north line and filled-arrow*/
	cairo_set_source_rgba(cr, 1., 0., 0., 1.0);
	cairo_save(cr); // save the original transform
	cairo_translate(cr, xdraw, ydraw);
	cairo_rotate(cr, angleN);
	cairo_move_to(cr, 0., 0.);
	cairo_line_to(cr, len, 0.);
	cairo_stroke(cr);
	cairo_line_to(cr, 0.75 * len, -0.15 * len);
	cairo_line_to(cr, 0.75 * len, +0.15 * len);
	cairo_line_to(cr, len, 0.);
	cairo_fill(cr);
	cairo_move_to(cr, len * 1.3, 0.1 * len);
	cairo_rotate(cr, -angleN);
	cairo_show_text(cr, "N");
	cairo_restore(cr); // restore the original transform

	/* draw east line */
	if (dd->neg_view)
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
	else
		cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);

	cairo_save(cr); // save the original transform
	cairo_translate(cr, xdraw, ydraw);
	cairo_rotate(cr, angleE);
	cairo_move_to(cr, 0., 0.);
	cairo_line_to(cr, len / 2.0, 0.);
	cairo_stroke(cr);
	cairo_move_to(cr, (len / 2) * 2.0, -0.1 * len);
	cairo_rotate(cr, -angleE);
	cairo_show_text(cr, "E");
	cairo_stroke(cr);
	cairo_restore(cr); // restore the original transform
}

static label_point *new_label_point(double height, const double *pix1, const double *pix2, const double *world, gboolean isRA, int border) {
	label_point *pt = g_new(label_point, 1);

	pt->x = pix1[0];
	pt->y = height - pix1[1];
	pt->ra = world[0];
	pt->dec = world[1];
	pt->angle = -atan2(pix2[1] - pix1[1], pix2[0] - pix1[0]);
	pt->isRA = isRA;
	pt->border = border;

	return pt;
}

static int has_pole(fits *fit) {
	if (!wcs2pix(fit, 0., 90., NULL, NULL))
		return 1;
	if (!wcs2pix(fit, 0., -90., NULL, NULL))
		return -1;
	return 0;
}

static gboolean get_line_intersection(double p0_x, double p0_y, double p1_x,
		double p1_y, double p2_x, double p2_y, double p3_x, double p3_y,
		double *i_x, double *i_y) {
	double s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;	s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;	s2_y = p3_y - p2_y;
	double det = -s2_x * s1_y + s1_x * s2_y;
	if (fabs(det) < DBL_EPSILON)
		return FALSE;
	double s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / det;
	t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / det;
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		// Collision detected
		if (i_x != NULL)
			*i_x = p0_x + (t * s1_x);
		if (i_y != NULL)
			*i_y = p0_y + (t * s1_y);
		return TRUE;
	}
	return FALSE; // No collision
}

static gint border_compare(const label_point *a, const label_point *b) {
	if (a->border > b->border) return 1;
	if (a->border < b->border) return -1;
	return 0;
}

static double ra_values[] = { 45, 30, 15, 10, 7.5, 5, 3.75, 2.5, 1.5, 1.25, 1, 3. / 4., 1.
		/ 2., 1. / 4., 1. / 6., 1. / 8., 1. / 12., 1. / 16., 1. / 24., 1. / 40., 1. / 48. };


static void draw_wcs_grid(const draw_data_t* dd) {
	if (!gui.show_wcs_grid) return;
	fits *fit = gfit;
	if (!has_wcs(fit)) return;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1. / dd->zoom);
	cairo_set_font_size(cr, 12.0 / dd->zoom);
	double ra0, dec0;
	GList *ptlist = NULL;
	double world[2], pix[2], pix2[2], img[2];
	double phi, theta;
	int status;

	double width = (double) fit->rx;
	double height = (double) fit->ry;
	cairo_rectangle(cr, 0., 0., width, height); // to clip the grid
	cairo_clip(cr);
	/* get ra and dec of center of the image */
	center2wcs(fit, &ra0, &dec0);
	if (ra0 == -1.) return;
	dec0 *= (M_PI / 180.0);
	ra0  *= (M_PI / 180.0);
	double range = get_wcs_image_resolution(fit) * sqrt(pow((width / 2.0), 2) + pow((height / 2.0), 2)); // range in degrees, FROM CENTER
	double step;

	/* Compute borders in pixel for tags*/
	const double pixbox[5][2] = { { 0., 0. }, { width, 0. }, { width, height }, { 0., height }, { 0., 0. } };
	const double pixval[4] = { 0., width, height, 0. }; // bottom, right, top, left with ref bottom left
	int pixtype[4] = { 1, 0, 1, 0 }; // y, x, y, x
	int polesign = has_pole(fit);

	/* calculate DEC step size */
	if (range > 16.0) {
		step = 8.; //step DEC 08:00
	} else if (range > 8.0) {
		step = 4.; // step DEC 04:00
	} else if (range > 4.0) { // image FOV about >2*4/sqrt(2) so >5 degrees
		step = 2.; // step DEC 02:00
	} else if (range > 2.0) {
		step = 1.; // step DEC 01:00
	} else if (range > 1.0) {
		step = 0.5; // step DEC 00:30
	} else if (range > 0.5) {
		step = 0.25; // step DEC 00:15
	} else if (range > 0.3) {
		step = 1. / 6.; // 0.166666, step DEC 00:10
	} else {
		step = 1. / 12.; // step DEC 00:05
	}

	// calculate RA step size
	double step2 = min(45, step / (cos(dec0) + 0.000001)); // exact value for stepRA, but not well rounded
	int iter = 0;
	double stepRA;
	do { // select nice rounded values for ra_step
		stepRA = ra_values[iter];
		iter++;
	} while ((stepRA >= step2) && (iter < G_N_ELEMENTS(ra_values))); // repeat until compatible value is found in ra_values
	if (polesign) stepRA = 45.;

	// round image centers
	double centra = stepRA * round(ra0 * 180 / (M_PI * stepRA)); // rounded image centers
	double centdec = step * round(dec0 * 180 / (M_PI * step));

	// plot DEC grid
	cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);
	double di = (polesign) ? 0. : centra - 6 * stepRA;
	do { // dec lines
		double dj = max(centdec - 6 * step, -90);
		do {
			double xa, ya, xb, yb, x1, x2, y1, y2;

			wcs2pix(fit, di, dj, &xa, &ya);
			x1 = round(xa - 1);
			y1 = round(height - ya);

			wcs2pix(fit, di, (dj + step), &xb, &yb);
			x2 = round(xb - 1);
			y2 = round(height - yb);

			if (((x1 >= 0) && (y1 >= 0) && (x1 < width) && (y1 < height))
					|| ((x2 >= 0) && (y2 >= 0) && (x2 < width) && (y2 < height))) {
				cairo_move_to(cr, x1, y1);
				cairo_line_to(cr, x2, y2);
				cairo_stroke(cr);
			}
			// check crossing
			if (!(((xa >= 0) && (ya >= 0) && (xa < width) && (ya < height))
						&& ((xb >= 0) && (yb >= 0) && (xb < width) && (yb < height)))) {
				for (int k = 0; k < 4; k ++) {
					if (get_line_intersection(xa, ya, xb, yb, pixbox[k][0], pixbox[k][1], pixbox[k+1][0], pixbox[k+1][1], NULL, NULL)) {
						world[0] = di;
						pix[pixtype[k]] = pixval[k];
						double latspan[2] = {dj, dj+step};
						status = wcsmix(fit->keywords.wcslib, pixtype[k], 1, latspan, 1.0, 0, world, &phi, &theta, img, pix);
						if(!status) {
							wcs2pix(fit, world[0], world[1] + 0.1, &pix2[0], &pix2[1]);
							ptlist = g_list_append(ptlist, new_label_point(height, pix, pix2, world, TRUE, k));
						}
						break;
					}
				}
			}
			dj = dj + step;
		} while (dj <= min(centdec + 6 * step, 90.));
		di = di + stepRA;
	} while (di <= ((polesign) ? 360. : centra + 6 * stepRA));

	// plot RA grid
	cairo_set_source_rgb(cr, 0.0, 0.5, 1.0);
	double dj = max(centdec - step * 6, -90);
	do { // ra lines
		di = (polesign) ? 0. : centra - 6 * stepRA;
		do {
			double xa, ya, xb, yb, x1, x2, y1, y2;

			wcs2pix(fit, di, dj, &xa, &ya);
			x1 = round(xa - 1);
			y1 = round(height - ya);

			wcs2pix(fit, (di + step), dj, &xb, &yb);
			x2 = round(xb - 1);
			y2 = round(height - yb);

			if (((x1 >= 0) && (y1 >= 0) && (x1 < width) && (y1 < height))
					|| ((x2 >= 0) && (y2 >= 0) && (x2 < width) && (y2 < height))) {
				cairo_move_to(cr, x1, y1);
				cairo_line_to(cr, x2, y2);
				cairo_stroke(cr);
			}
				// check crossing
			if (!(((xa >= 0) && (ya >= 0) && (xa < width) && (ya < height))
				&& ((xb >= 0) && (yb >= 0) && (xb < width) && (yb < height)))) {
				for (int k = 0; k < 4; k ++) {
					if (get_line_intersection(xa, ya, xb, yb, pixbox[k][0], pixbox[k][1], pixbox[k+1][0], pixbox[k+1][1], NULL, NULL)) {
						world[1] = dj;
						pix[pixtype[k]] = pixval[k];
						double lngspan[2] = {di, di+step};
						status = wcsmix(fit->keywords.wcslib, pixtype[k], 2, lngspan, 1.0, 0, world, &phi, &theta, img, pix);
						if(!status) {
							wcs2pix(fit, world[0] + 0.1, world[1], &pix2[0], &pix2[1]);
							ptlist = g_list_append(ptlist, new_label_point(height, pix, pix2, world, FALSE, k));
						}
						break;
					}
				}
			}
			di = di + step;
		} while (di <= ((polesign) ? 360. : centra + 6 * stepRA));
		dj = dj + step;
	} while (dj <= min(centdec + step * 6, 90));

	// Add crossings labels
	ptlist = g_list_sort(ptlist, (GCompareFunc) border_compare); // sort potential tags by increasing border number
	if (dd->neg_view) {
		cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
	} else {
		cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
	}
	GSList *existingtags = NULL;
	gchar *RAfmt = (stepRA < 1./4.) ? "%02dh%02dm%02ds" : "%02dh%02dm";
	for (GList *l = ptlist; l != NULL; l = l->next) {
		// getting the label
		SirilWorldCS *world_cs;
		label_point *pt = (label_point*) l->data;
		world_cs = siril_world_cs_new_from_a_d(pt->ra, pt->dec);
		if (world_cs) {
			gchar *tag = (pt->isRA) ? siril_world_cs_alpha_format(world_cs, RAfmt) : siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'");
			siril_world_cs_unref(world_cs);
			if (!g_slist_find_custom(existingtags, tag, (GCompareFunc) strcompare)) { // this tag has already been used - skipping
				existingtags = g_slist_append(existingtags, (gpointer) tag);
				cairo_text_extents_t te1, te2;
				cairo_text_extents(cr, tag, &te1); // getting the dimensions of the textbox
				cairo_save(cr); // save the orginal transform
				cairo_translate(cr, pt->x, pt->y);
				// add pi for angles larger than +/- pi/2
				if (pt->angle > M_PI_2)
					pt->angle -= M_PI;
				if (pt->angle < -M_PI_2)
					pt->angle += M_PI;
				double dx = 0., dy = 0.;
				switch (pt->border) { // shift to get back in the image
				case 0: // bottom
					if (pt->angle > 0.)
						dx -= te1.x_advance;
					break;
				case 1: // right
					dx -= te1.x_advance;
					if (pt->angle > 0.)
						dy += te1.height;
					break;
				case 2: // top
					dy += te1.height;
					if (pt->angle < 0.)
						dx -= te1.x_advance;
					break;
				case 3: // left
					if (pt->angle < 0.)
						dy += te1.height;
					break;
				default:
					break;
				}
				cairo_rotate(cr, pt->angle);
				cairo_move_to(cr, dx, dy);
				cairo_text_extents(cr, tag, &te2);
				cairo_show_text(cr, tag);
				cairo_stroke(cr);
				cairo_restore(cr); // restore the orginal transform
			} else {
				g_free(tag);
			}
		}
	}
	g_list_free_full(ptlist, (GDestroyNotify) g_free);
	g_slist_free_full(existingtags, (GDestroyNotify) g_free);

	draw_compass(dd);
}

static void draw_wcs_disto(const draw_data_t* dd) {
	if (!gui.show_wcs_disto) return;
	fits *fit = gfit;
	if (!has_wcs(fit) || !fit->keywords.wcslib->lin.dispre) return; // no platesolve or no distortions
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 3. / dd->zoom);
	cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);

	int nbpoints = 20;
	double radius = (dd->zoom < 1. )?  3. : 3. / dd->zoom;

	double paceX = (double)fit->rx / (double)(nbpoints - 1);
	double paceY = (double)fit->ry / (double)(nbpoints - 1);
	double currX = 0.5; // coords in fits/wcs conventions start at (0.5, 0.5) for the bottom left corner
	double currY = 0.5;
	for (int i = 0; i < nbpoints; i++) {
		currY = 0.5;
		for (int j = 0; j < nbpoints; j++) {
			double rawcrd[2] = { currX, currY };
			double discrd[2] = {0., 0.};
			int status = disp2x(fit->keywords.wcslib->lin.dispre, rawcrd, discrd);
			if (!status) {
				double disX = (discrd[0] - rawcrd[0]) * 5. +  rawcrd[0];
				double disY = (discrd[1] - rawcrd[1]) * 5. +  rawcrd[1];
				double startX, startY, endX, endY;
				fits_to_display(currX, currY, &startX, &startY, fit->ry);
				fits_to_display(disX, disY, &endX, &endY, fit->ry);
				cairo_arc(cr, startX, startY, radius,  0., 2. * M_PI);
				cairo_fill(cr);
				cairo_move_to(cr, startX, startY);
				cairo_line_to(cr, endX, endY);
				cairo_stroke(cr);
			}
			currY += paceY;
		}
		currX += paceX;
	}
}

static gdouble x_circle(gdouble x, gdouble radius, gdouble angle) {
	return x + radius * cos(angle);
}

static gdouble y_circle(gdouble y, gdouble radius, gdouble angle) {
	return y + radius * sin(angle);
}

static void draw_annotates(const draw_data_t* dd) {
	if (!com.found_object) return;
	gdouble resolution = get_wcs_image_resolution(gfit);
	if (resolution <= 0) return;
	double width = (double) gfit->rx;
	double height = (double) gfit->ry;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);

	cairo_set_line_width(cr, 1.0 / dd->zoom);
	cairo_rectangle(cr, 0., 0., width, height); // to clip the grid
	cairo_clip(cr);

	for (GSList *list = com.found_object; list; list = list->next) {
		CatalogObjects *object = (CatalogObjects *)list->data;
		gdouble radius = get_catalogue_object_radius(object);
		gdouble x = get_catalogue_object_x(object);
		gdouble y = get_catalogue_object_y(object);
		gdouble x1 = get_catalogue_object_x1(object);
		gdouble y1 = get_catalogue_object_y1(object);
		gchar *code = get_catalogue_object_code_pretty(object);
		guint catalog = get_catalogue_object_cat(object);
		gboolean revert = FALSE;
		double angle = ANGLE_TOP;
		double addoffset = 0.;
		GdkRGBA gdk_color;

		switch (catalog) {
		case CAT_AN_USER_DSO:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_dso_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			break;
		case CAT_AN_USER_SSO:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_sso_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			break;
		case CAT_AN_USER_TEMP:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_tmp_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			revert = TRUE;
			angle = ANGLE_BOT;
			break;
		default:
		case 0:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_std_annotations);
			if (dd->neg_view) {
				cairo_set_source_rgba(cr, 1.0 - gdk_color.red, 1.0 - gdk_color.green, 1.0 - gdk_color.blue, gdk_color.alpha);
			} else {
				cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			}
			break;
		}

		radius = radius / resolution / 60.0;
		// radius now in pixels

		point offset = {5., revert ? 5. : -5.};
		if (catalog == CAT_AN_CONST) { // constellation line
			cairo_move_to(cr, x, y);
			cairo_line_to(cr, x1, y1);
			cairo_stroke(cr);
		} else if (radius < 0 || catalog == CAT_AN_CONST_NAME) {
			// objects we don't have an accurate location (LdN, Sh2)
		} else if (radius > 5) {
			cairo_arc(cr, x, y, radius, 0., 2. * M_PI);
			cairo_stroke(cr);
			if (code) {
				cairo_move_to(cr, x_circle(x, radius, angle), y_circle(y, radius, angle));
				offset.x = x_circle(x, radius * 1.3, angle) - x;
				offset.y = y_circle(y, radius * 1.3, angle) - y;
				cairo_line_to(cr, offset.x + x, offset.y + y);
				cairo_stroke(cr);
			}
		} else {
			/* it is punctual */
			cairo_move_to(cr, x, y - 15);
			cairo_line_to(cr, x, y - 5);
			cairo_stroke(cr);
			cairo_move_to(cr, x, y + 15);
			cairo_line_to(cr, x, y + 5);
			cairo_stroke(cr);
			cairo_move_to(cr, x - 15, y);
			cairo_line_to(cr, x - 5, y);
			cairo_stroke(cr);
			cairo_move_to(cr, x + 15, y);
			cairo_line_to(cr, x + 5, y);
			cairo_stroke(cr);
		}
		if (code) {
			gchar *name = code;
			gdouble size = 18 * (com.pref.gui.font_scale / 100.0);
			cairo_select_font_face(cr, "Liberation Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, size / dd->zoom);
			if (revert) {
				cairo_text_extents_t te;
				cairo_text_extents(cr, name, &te); // getting the dimensions of the textbox
				addoffset = te.height;
			}
			cairo_move_to(cr, x + offset.x, y + offset.y + addoffset);
			cairo_show_text(cr, name);
			cairo_stroke(cr);
		}

	}
}

static void draw_rgb_centers(const draw_data_t* dd) {
	if (!gui.comp_layer_centering) return;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);

	double red = gui.comp_layer_centering->saturated_color.red;
	double green = gui.comp_layer_centering->saturated_color.green;
	double blue = gui.comp_layer_centering->saturated_color.blue;
	cairo_set_source_rgb(cr, red, green, blue);
	cairo_set_line_width(cr, 2.0 / dd->zoom);

	double size = 10. / dd->zoom;
	cairo_arc(cr, gui.comp_layer_centering->center.x, gui.comp_layer_centering->center.y, size * 0.5, 0., 2. * M_PI);
	cairo_stroke(cr);
}

static void draw_analysis(const draw_data_t* dd) {
	if (com.tilt) {
		cairo_t *cr = dd->cr;
		cairo_set_dash(cr, NULL, 0, 0);

		cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
		cairo_set_line_width(cr, 2.0 / dd->zoom);
		cairo_move_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y);
		cairo_line_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y);
		cairo_line_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y);
		cairo_line_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y);
		cairo_line_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y);
		cairo_move_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y);
		cairo_line_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y);
		cairo_line_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y);
		cairo_stroke(cr);

		/* draw text */
		cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		int size = 20.0 / dd->zoom;
		cairo_set_font_size(cr, size);

		/* fwhm 1 */
		gchar *str = g_strdup_printf("%.2f", com.tilt->fwhm[0]);
		cairo_move_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y - size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 2 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[1]);
		cairo_move_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y - size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 3 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[2]);
		cairo_move_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y + size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 4 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[3]);
		cairo_move_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y + size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm center */
		str = g_strdup_printf("%.2f", com.tilt->fwhm_centre);
		cairo_move_to(cr, gfit->rx / 2.0, (gfit->ry / 2.0) + size);
		cairo_show_text(cr, str);
		cairo_stroke(cr);
		g_free(str);
	}
}

static void draw_regframe(const draw_data_t* dd) {
	if (com.script || com.headless) return;
	if (!sequence_is_loaded()) return;
	if (com.seq.current == RESULT_IMAGE) return;
	if (!drawframe)
		image_display_init_statics();
	if (!siril_toggle_get_active(GTK_WIDGET(drawframe))) return;
	int activelayer = gtk_drop_down_get_selected(seqcombo);
	if (!layer_has_registration(&com.seq, activelayer)) {
		activelayer = seq_has_any_regdata(&com.seq);
		if (activelayer < 0) return;
	}
	if (com.seq.reg_invalidated) return;
	transformation_type min, max;
	guess_transform_from_seq(&com.seq, activelayer, &min, &max, FALSE);
	if (max <= IDENTITY_TRANSFORMATION) return;

	if (guess_transform_from_H(com.seq.regparam[activelayer][com.seq.current].H) == NULL_TRANSFORMATION)
		return; // current image H matrix is null matrix

	// For ext_ref sequences, H matrices are absolute (relative to external ref).
	// Use identity as Href so the frame shows the external reference boundary.
	// For normal sequences, use H[reference_image] so the frame shows the sequence reference.
	Homography Href = { 0 };
	int ref_rx, ref_ry;
	if (com.seq.ext_ref) {
		cvGetEye(&Href);
		ref_rx = com.seq.ext_ref_rx;
		ref_ry = com.seq.ext_ref_ry;
	} else {
		if (guess_transform_from_H(com.seq.regparam[activelayer][com.seq.reference_image].H) == NULL_TRANSFORMATION)
			return;
		Href = com.seq.regparam[activelayer][com.seq.reference_image].H;
		ref_rx = com.seq.is_variable ? com.seq.imgparam[com.seq.reference_image].rx : com.seq.rx;
		ref_ry = com.seq.is_variable ? com.seq.imgparam[com.seq.reference_image].ry : com.seq.ry;
	}

	regframe framing = { 0 };
	framing.pt[0].x = 0.;
	framing.pt[0].y = 0.;
	framing.pt[1].x = (double)ref_rx;
	framing.pt[1].y = 0.;
	framing.pt[2].x = (double)ref_rx;
	framing.pt[2].y = (double)ref_ry;
	framing.pt[3].x = 0.;
	framing.pt[3].y = (double)ref_ry;
	double cogx = 0., cogy = 0., cx, cy;
	for (int i = 0; i < 4; i++) {
		cvTransfPoint(&framing.pt[i].x, &framing.pt[i].y, Href, com.seq.regparam[activelayer][com.seq.current].H, 1.);
		cogx += framing.pt[i].x;
		cogy += framing.pt[i].y;
	}
	cogx *= 0.25;
	cogy *= 0.25;
	cx = (com.seq.is_variable) ? (double)com.seq.imgparam[com.seq.current].rx * 0.5 : (double)com.seq.rx * 0.5;
	cy = (com.seq.is_variable) ? (double)com.seq.imgparam[com.seq.current].ry * 0.5 : (double)com.seq.ry * 0.5;

	cairo_t *cr = dd->cr;
	double size = 10. / dd->zoom;
	cairo_set_dash(cr, NULL, 0, 0);
	gboolean has_disto = seq_has_any_distortion(&com.seq);
	if (max <= SHIFT_TRANSFORMATION)
		cairo_set_source_rgb(cr, 0.0, 0.5, 1.0);
	else {
		if (has_disto)
			cairo_set_source_rgb(cr, 0., 1., 0.5);
		else
			cairo_set_source_rgb(cr, 1., 0., 0.);
	}

	cairo_set_line_width(cr, 2.0 / dd->zoom);
	// reference origin
	cairo_arc(cr, framing.pt[0].x, framing.pt[0].y, size * 0.5, 0., 2. * M_PI);
	cairo_stroke_preserve(cr);
	cairo_fill(cr);
	// reference frame
	cairo_move_to(cr, framing.pt[0].x, framing.pt[0].y);
	cairo_line_to(cr, framing.pt[1].x, framing.pt[1].y);
	cairo_line_to(cr, framing.pt[2].x, framing.pt[2].y);
	cairo_line_to(cr, framing.pt[3].x, framing.pt[3].y);
	cairo_line_to(cr, framing.pt[0].x, framing.pt[0].y);
	cairo_stroke(cr);

	// reference center
	cairo_arc(cr, cogx, cogy, size * 0.5, 0., 2. * M_PI);
	cairo_stroke(cr);
	// current center
	cairo_set_source_rgb(cr, 0., 1., 0.);
	cairo_move_to(cr, cx - size * 0.5, cy);
	cairo_rel_line_to(cr, size, 0.);
	cairo_stroke(cr);
	cairo_move_to(cr, cx, cy - size * 0.5);
	cairo_rel_line_to(cr, 0., size);
	cairo_stroke(cr);
}

void initialize_image_display() {
	int i;
	siril_log_debug("HD AutoStretch bitdepth: %d\n", com.pref.hd_bitdepth);
	gui.hd_remap_max = 1 << (guint) com.pref.hd_bitdepth;
	for (i = 0; i < MAXGRAYVPORT; i++) {
		memset(gui.remap_index[i], 0, sizeof(gui.remap_index[i]));
		last_pente = 0.f;
		last_mode = HISTEQ_DISPLAY;
		// only HISTEQ mode always computes the index, it's a good initializer here
	}
	cairo_matrix_init_identity(&gui.display_matrix);
}

/* this function calculates the "fit to window" zoom values, given the window
 * size in argument and the image size in gfit->
 * Should not be called before displaying the main gray window when using zoom to fit */
double get_zoom_val() {
	int window_width = 0, window_height = 0;
	if (gui.zoom_value > 0.)
		return gui.zoom_value;
	/* else if zoom is < 0, it means fit to window */
	/* Measure the viewport that is actually on screen.  In GTK4 a GtkNotebook
	 * only allocates its visible page, so the off-screen vports' drawareas
	 * report a zero size.  Hardcoding RED_VPORT therefore pinned colour images
	 * (displayed in RGB_VPORT) to the 1.0 fallback, breaking fit-to-window. */
	int vport = gui.cvport;
	if (vport >= 0 && vport < MAXVPORT && gui.view[vport].drawarea) {
		window_width = gtk_widget_get_width(gui.view[vport].drawarea);
		window_height = gtk_widget_get_height(gui.view[vport].drawarea);
	}
	if (window_width <= 1 || window_height <= 1) {
		/* cvport not allocated yet or stale: use whichever vport currently
		 * has a real allocation (only the visible notebook page does). */
		for (int i = 0; i < MAXVPORT; i++) {
			if (!gui.view[i].drawarea)
				continue;
			int w = gtk_widget_get_width(gui.view[i].drawarea);
			int h = gtk_widget_get_height(gui.view[i].drawarea);
			if (w > 1 && h > 1) {
				window_width = w;
				window_height = h;
				break;
			}
		}
	}
	if (gfit->rx == 0 || gfit->ry == 0 || window_height <= 1 || window_width <= 1)
		return 1.0;
	double wtmp = (double) window_width / (double) gfit->rx;
	double htmp = (double) window_height / (double) gfit->ry;
	return min(wtmp, htmp);
}

static void invalidate_image_render_cache(int vport) {
	/* Phase 2 retired the disp_surface render cache; GSK now caches the
	 * snapshot render-node tree per frame.  The cache-invalidation call
	 * sites are kept (they're harmless no-ops here) so that future
	 * additions to the render cache can hook back in without auditing
	 * every remap exit point. */
	(void) vport;
}

void adjust_vport_size_to_image() {
	if (com.script) return;
	double zoom = get_zoom_val();
	if (zoom <= 0.0) return;
	/* Init display matrix from current display state */
	cairo_matrix_t new_matrix;
	/*siril_log_debug("computing matrix for zoom %g and offset [%g, %g]\n",
			zoom, gui.display_offset.x, gui.display_offset.y);*/
	cairo_matrix_init(&new_matrix,
			zoom, 0, 0, zoom,
			gui.display_offset.x,
			gui.display_offset.y);
	if (memcmp(&new_matrix, &gui.display_matrix, sizeof(new_matrix))) {
		invalidate_image_render_cache(-1);
		gui.display_matrix = new_matrix;

		/* Compute the inverse display matrix used for coordinate transformation */
		gui.image_matrix = gui.display_matrix;
		cairo_matrix_invert(&gui.image_matrix);
		//siril_log_debug("  matrix changed\n");
	}
}

void copy_roi_into_gfit() {
	size_t npixels_roi = gui.roi.selection.w * gui.roi.selection.h;
	if (npixels_roi == 0 || com.script || com.python_command)
		return;
	g_rw_lock_writer_lock(&gfit->rwlock);
	size_t npixels_gfit = gfit->rx * gfit->ry;
	if (gui.roi.fit.type != gfit->type) {
		size_t roi_ndata = gui.roi.fit.rx * gui.roi.fit.ry * gui.roi.fit.naxes[2];
		if (gfit->type == DATA_FLOAT) {
			fit_replace_buffer(&gui.roi.fit, gui.roi.fit.bitpix == BYTE_IMG ? ushort8_buffer_to_float(gui.roi.fit.data, roi_ndata): ushort_buffer_to_float(gui.roi.fit.data, roi_ndata), DATA_FLOAT);
			if (is_preview_active()) {
				fits *roi_backup = get_roi_backup();
				fit_replace_buffer(roi_backup, roi_backup->bitpix == BYTE_IMG ? ushort8_buffer_to_float(roi_backup->data, roi_ndata) : ushort_buffer_to_float(roi_backup->data, roi_ndata), DATA_FLOAT);
			}
		} else {
			fit_replace_buffer(&gui.roi.fit, float_buffer_to_ushort(gui.roi.fit.fdata, roi_ndata), DATA_USHORT);
			if (gfit->bitpix == BYTE_IMG) {
				for (size_t i = 0 ; i < roi_ndata ; i++) {
					gui.roi.fit.data[i] >>= 8;
				}
				gui.roi.fit.bitpix = BYTE_IMG;
			}
			if (is_preview_active()) {
				fits *roi_backup = get_roi_backup();
				fit_replace_buffer(roi_backup, float_buffer_to_ushort(roi_backup->fdata, roi_ndata), DATA_USHORT);
				if (gfit->bitpix == BYTE_IMG) {
					for (size_t i = 0 ; i < roi_ndata ; i++) {
						roi_backup->data[i] >>= 8;
					}
					roi_backup->bitpix = BYTE_IMG;
				}
			}
		}
	}

	if (gui.roi.fit.type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (uint32_t c = 0 ; c < gui.roi.fit.naxes[2] ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				const float *rowindex = gui.roi.fit.fdata + (y * gui.roi.fit.rx) + (c * npixels_roi);
				float *destindex = gfit->fdata + (c * npixels_gfit) + ((gfit->ry - gui.roi.selection.y - y - 1) * gfit->rx) + gui.roi.selection.x;
				memcpy(destindex, rowindex, gui.roi.selection.w * sizeof(float));
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (uint32_t c = 0 ; c < gui.roi.fit.naxes[2] ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				const WORD *rowindex = gui.roi.fit.data + (y * gui.roi.fit.rx) + (c * npixels_roi);
				WORD *destindex = gfit->data + (npixels_gfit * c) + ((gfit->ry - gui.roi.selection.y - y - 1) * gfit->rx) + gui.roi.selection.x;
				memcpy(destindex, rowindex, gui.roi.selection.w * sizeof(WORD));
			}
		}
	}
	g_rw_lock_writer_unlock(&gfit->rwlock);
}

void remap_all() {
	stf_computed = FALSE;
	if (gui.rendering_mode == HISTEQ_DISPLAY || gui.rendering_mode == STF_DISPLAY) {
		for (int i = 0; i < gfit->naxes[2]; i++) {
			remap(i);
		}
	} else {
		remap_all_vports();
	}
	if (gfit->naxis == 3)
		remaprgb();
}

void redraw(remap_type doremap) {
	if (com.script && !com.python_script) return;
	switch (doremap) {
		case REDRAW_OVERLAY:
			break;
		case REDRAW_IMAGE:
			invalidate_image_render_cache(-1);
			break;
		case REDRAW_ALL:
			/* redraw the 9-panel mosaic dialog if needed */
			redraw_aberration_inspector();
			break;
		default:
			siril_log_debug("UNKNOWN REMAP\n\n");
	}
	request_gtk_redraw_of_cvport();
}

static gboolean redraw_idle(gpointer p) {
	redraw((remap_type)GPOINTER_TO_INT(p)); // draw stars
	return FALSE;
}

static gpointer redraw_idle_thread_func(gpointer data) {
	sample_mutex_lock();
	execute_idle_and_wait_for_it(redraw_idle, data);
	sample_mutex_unlock();
	return FALSE;
}

void queue_redraw_and_wait_for_it(remap_type doremap) {
	if (!com.script && !com.python_command && !com.headless) {
		GThread *thread = g_thread_new("redraw", redraw_idle_thread_func, GINT_TO_POINTER((int)doremap));
		g_thread_join(thread);
	}
}

gboolean redraw_mask_idle(gpointer p) {
	g_rw_lock_reader_lock(&gfit->rwlock);
	if (gfit->mask && gfit->mask->data)
		remap_mask(gfit->mask);
	g_rw_lock_reader_unlock(&gfit->rwlock);
	if (com.pref.gui.mask_tints_vports) {
		/* Reader lock released before this call: notify_gfit_data_modified()
		 * may call copy_roi_into_gfit() which acquires the writer lock. */
		notify_gfit_data_modified();
		redraw(REDRAW_ALL); // need to remap all to tint the image vports correctly
	}
	return FALSE;
}

void queue_redraw_mask() {
	siril_add_idle(redraw_mask_idle, NULL);
}

void queue_redraw(remap_type doremap) {
	// request a redraw from another thread
	siril_add_idle(redraw_idle, GINT_TO_POINTER((int)doremap));
}

/* The legacy GtkDrawingArea-driven Cairo redraw_drawingarea entry point
 * was removed in Phase 2.  All on-screen image rendering now flows
 * through SirilImageView::snapshot, and the Cairo save / clipboard
 * pathway calls draw_main_image directly via
 * add_image_and_label_to_cairo. */

point get_center_of_vport() {
	image_display_init_statics();
	GtkWidget *widget = imgdisp_drawing_r;

	guint window_width = gtk_widget_get_width(widget);
	guint window_height = gtk_widget_get_height(widget);

	point center = { window_width / 2., window_height / 2. };

	return center;
}

void add_image_and_label_to_cairo(cairo_t *cr, int vport) {
	draw_data_t dd;
	image_display_init_statics();
	GtkWidget *widget = imgdisp_drawing_r;
	GAction *action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
	GVariant *state = g_action_get_state(action_neg);

	dd.vport = vport;
	dd.cr = cr;
	dd.window_width = gtk_widget_get_width(widget);
	dd.window_height = gtk_widget_get_height(widget);
	dd.zoom = get_zoom_val();
	dd.image_width = gfit->rx;
	dd.image_height = gfit->ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_NEAREST;
	dd.neg_view = g_variant_get_boolean(state);
	g_variant_unref(state);

	/* RGB or gray images */
	draw_main_image(&dd);
	/* detected stars and highlight the selected star */
	g_mutex_lock(&com.mutex);
	draw_stars(&dd);
	g_mutex_unlock(&com.mutex);
	/* wcs_grid */
	draw_wcs_grid(&dd);
	/* detected objects */
	draw_annotates(&dd);
	/* analysis */
	draw_analysis(&dd);
	/* distortions */
	draw_wcs_disto(&dd);
}
