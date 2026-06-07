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

/*
 * GUI-only data structures and global state.
 *
 * This header defines types that carry GTK/GDK/Cairo fields and therefore
 * must never be included from non-GUI compilation units.  All processing,
 * I/O, registration, and stacking code must remain free of this header.
 *
 * Moved here from core/siril.h as part of plan step 1.2.
 */

#ifndef GUI_STATE_H
#define GUI_STATE_H

#include <gtk/gtk.h>    /* GtkWidget, GtkBuilder, GdkRGBA, cairo_*, … */
#include "core/siril.h" /* guiinfo typedef, MAXVPORT, PREVIEW_NB, and all non-GTK types */

/* Convenience alias kept in GUI translation units only. */
typedef GtkWidget SirilWidget;

/* ── compositing layer ────────────────────────────────────────────────────── */

/* One source image and its associated colour in the RGB compositor.
 * Defined here rather than in siril.h because every field that is not a
 * plain numeric value is a GTK widget pointer or GdkRGBA. */
typedef struct {
	/* widgets */
	GtkButton       *remove_button;
	GtkDrawingArea  *color_w;           /* simulated colour chooser */
	GtkButton       *chooser_button;    /* file chooser */
	gchar           *selected_filename;
	GtkLabel        *chooser_button_label;
	GtkLabel        *label;
	GtkSpinButton   *spinbutton_x;
	GtkSpinButton   *spinbutton_y;
	GtkSpinButton   *spinbutton_r;
	GtkToggleButton *centerbutton;
	double           spinbutton_x_value;
	double           spinbutton_y_value;
	double           spinbutton_r_value;
	/* data */
	GdkRGBA          color;             /* layer colour in image colorspace */
	GdkRGBA          saturated_color;   /* saturated layer colour */
	GdkRGBA          display_color;     /* layer colour in display colorspace */
	fits             the_fit;           /* the FITS for this layer */
	point            center;
} layer;

/* ── viewport rendering cache ─────────────────────────────────────────────── */

/* Rendering cache for one display viewport.  The raw pixel buffer (buf) is
 * written by worker threads; all other fields are GTK/Cairo and live only
 * on the main thread. */
/* Per-tile state.  In eager mode tile->data is a pointer into the shared
 * view->buf (no separate allocation); in lazy mode each tile owns its own
 * malloc'd byte buffer that can be freed under LRU pressure.  In both
 * modes tile->texture wraps tile->data (or its slice of buf) via a static
 * GBytes; both pointers are NULL when the tile is evicted / not yet
 * materialised. */
struct image_tile {
	guchar     *data;
	GdkTexture *texture;
	guint64     last_used;   /* lazy-mode LRU epoch */
	gboolean    dirty;       /* lazy-mode flag: LUT/data changed since last fill */
	int         mip;         /* lazy-mode downsample level: texture dim = src_dim >> mip;
	                          * 0 = full-resolution (eager mode is always 0).
	                          * Set by materialise_tile_{gray,rgb}; consumed by tile_release
	                          * for accurate budget bookkeeping and by the snapshot path to
	                          * detect "current mip differs from target mip" → re-materialise. */
	gsize       bytes;       /* size in bytes of tile->data when materialised; tracks the
	                          * actual allocation so tile_release returns the right amount
	                          * to lazy_bytes_used regardless of mip/edge size. */
};

struct image_view {
	GtkWidget        *drawarea;

	/* Image dimensions the tile grid represents. */
	int               buf_stride;        /* bytes per image row (= img_w * 4 rounded) */
	int               buf_height;        /* img_h */
	int               view_width;        /* drawing area dimensions */
	int               view_height;

	/* Tile grid.  Each tile is at most SIRIL_TILE_DIM × SIRIL_TILE_DIM pixels;
	 * right- and bottom-edge tiles may be smaller. */
	int               tile_dim;
	int               tile_cols, tile_rows;
	struct image_tile *tiles;            /* tile_cols * tile_rows entries */

	/* Eager-mode contiguous buffer: a single allocation of the whole image,
	 * with each tile->data pointing into it (each texture wraps a sub-region
	 * via buf_stride).  NULL in lazy mode.  buf is borrowed from buf_gbytes
	 * for convenient indexing — the refcounted GBytes is what actually owns
	 * the allocation, so the bytes outlive any tile texture that GSK is
	 * still rendering. */
	guchar           *buf;
	GBytes           *buf_gbytes;

	/* Background materialise worker (lazy mode).  The heavy per-pixel tile
	 * fill runs on a shared single-thread GThreadPool rather than the GTK
	 * main thread, so neither the fill nor cairo_mutex contention can stall
	 * rendering.  worker_active (set under cairo_mutex) is a dedupe flag: the
	 * snapshot pushes a job only when no worker is already churning this view.
	 * generation is bumped whenever the tile grid is reallocated; the worker
	 * stamps each job with it and discards a finished texture whose generation
	 * no longer matches (the grid moved under it).  wk_target_mip and
	 * render_neg are the bits of render state the worker can't compute itself
	 * off the main thread (zoom/widget size, GAction state); the snapshot
	 * records them under cairo_mutex each frame for the worker to read. */
	gboolean          worker_active;
	guint             generation;
	/* Bumped by every lazy invalidate (LUT change).  A worker job captures it
	 * alongside its LUT pointers and discards its finished texture if it changed
	 * during the unlocked fill — otherwise the worker would clear the tile's
	 * dirty flag that the mid-flight invalidate just set, leaving the tile
	 * clean-but-stale until a later zoom-in changes its target mip. */
	guint             invalidate_seq;
	int               wk_target_mip;
	gboolean          render_neg;

	/* Lazy mode: tile bytes are allocated on demand by materialise_tile and
	 * may be freed by LRU eviction once lazy_bytes_used exceeds the budget. */
	gboolean          lazy;
	gsize             lazy_bytes_used;
	guint64           lazy_epoch;

	/* Lazy-mode low-resolution proxy: a single point-sampled downscale of the
	 * whole image (long edge <= SIRIL_PROXY_MAX), built from the current LUTs
	 * and kept resident as an always-available fallback.  When a fast zoom-out
	 * makes many never-materialised tiles visible at once, the snapshot draws
	 * this proxy as a full-image backdrop instead of blocking to materialise
	 * each tile — the real, sharp tiles then fill in over the next few frames
	 * via the refine idle.  proxy_mip / proxy_w / proxy_h describe the current
	 * texture; proxy_dirty is set on every lazy remap and grid realloc so the
	 * next snapshot rebuilds it.  Not counted against lazy_bytes_used (small,
	 * always resident).  NULL / unused in eager mode. */
	GdkTexture       *proxy;
	int               proxy_mip;
	int               proxy_w, proxy_h;
	gboolean          proxy_dirty;

	/* Inclusive visible tile-index range recorded by the most recent snapshot,
	 * consumed by the materialise worker so it only does work for on-screen
	 * tiles.  vis_valid is FALSE until the first snapshot populates it. */
	int               vis_tx0, vis_ty0, vis_tx1, vis_ty1;
	gboolean          vis_valid;
};

/* ── draw callback context ────────────────────────────────────────────────── */

/* Bundle passed to draw_extra callbacks and overlay painters. */
typedef struct draw_data {
	cairo_t         *cr;            /* destination context */
	int              vport;         /* viewport index */
	double           zoom;          /* current zoom value */
	gboolean         neg_view;      /* negative (inverted) display */
	cairo_filter_t   filter;        /* image scaling filter */
	guint            image_width, image_height;   /* image dimensions */
	guint            window_width, window_height; /* drawing-area dimensions */
} draw_data_t;

/* ── master GUI state ─────────────────────────────────────────────────────── */

struct guiinf {
	GtkBuilder       *builder;          /* GtkBuilder for the whole interface */

	/*** rendering of the currently loaded image ***/
	struct image_view view[MAXVPORT];
	int               cvport;           /* active viewport index */

	cairo_matrix_t    display_matrix;   /* image → display coordinate transform */
	cairo_matrix_t    image_matrix;     /* display → image coordinate transform */
	double            zoom_value;       /* 1.0 = 100 %; use get_zoom_val() */
	point             display_offset;   /* image pan offset */

	gboolean          translating;      /* panning in progress */

	/*** FLIS canvas drag-layer state (§4.3 — toolbar drag toggle).
	 * Cursor start position is in UNCLAMPED image coords (point,
	 * doubles) so dragging past the canvas edge keeps producing
	 * correct deltas — FLIS supports negative position_x/y and
	 * positions beyond the canvas extent (sparse layers). ***/
	gboolean          flis_layer_dragging;
	gint              flis_drag_layer_id;
	point             flis_drag_start_image;
	pointi            flis_drag_start_layer;

	gboolean          show_excluded;    /* show excluded images in sequences */

	int               selected_star;    /* selected star in the star list */

	gboolean          show_wcs_grid;
	gboolean          show_wcs_disto;

	psf_star         *qphot;            /* quick-photometry highlight */

	point             measure_start;    /* alt-drag measurement endpoints */
	point             measure_end;

	GSList           *user_polygons;    /* user-drawn polygon overlays */

	void            (*draw_extra)(draw_data_t *dd); /* extra overlay painter */

	/*** colour mapping ***/
	WORD              lo, hi;           /* cutoff slider values */
	gboolean          cut_over;         /* show values > hi as inverted */
	sliders_mode      sliders;
	display_mode      rendering_mode;
	gboolean          unlink_channels;  /* for autostretch */
	BYTE              remap_index[3][USHRT_MAX + 1]; /* 8-bit LUT per channel */
	BYTE             *hd_remap_index[3];             /* high-precision LUT */
	guint             hd_remap_max;
	gboolean          use_hd_remap;

	/*** selection / drawing state ***/
	gboolean          drawing;
	cut_struct        cut;
	pointi            start;
	pointi            origin;
	gboolean          freezeX, freezeY;
	double            ratio;
	double            rotation;

	/*** alignment preview ***/
	cairo_surface_t  *preview_surface[PREVIEW_NB];
	GtkWidget        *preview_area[PREVIEW_NB];
	guchar           *refimage_regbuffer;
	cairo_surface_t  *refimage_surface;

	int               file_ext_filter;

	/*** command history ***/
	char            **cmd_history;
	int               cmd_hist_size;
	int               cmd_hist_current;
	int               cmd_hist_display;

	layer            *comp_layer_centering; /* RGB compositor centering target */

	roi_t             roi;
	GSList           *mouse_actions;
	GSList           *scroll_actions;
	gboolean          drawing_polygon;
	GSList           *drawing_polypoints;
	GdkRGBA           poly_ink;         /* polygon drawing colour */
	gboolean          poly_fill;

	/* cairo_mutex guards view[].buf and all Cairo display surfaces.
	 * Held by display-update code on the main thread and any worker that
	 * composites into the display buffers. */
	GMutex            cairo_mutex;

	/* Set atomically by generic_image_worker before UI updates; cleared
	 * just before the completion idle is posted.  When set,
	 * redraw_drawingarea repaints from the cached surface to avoid a
	 * grey flash while processing. */
	gint              suppress_drawarea_redraw;
};

/* The single global GUI state instance.  Defined in main.c / main-cli.c.
 * Only GUI translation units should reference this directly; processing code
 * must use gui_iface callbacks instead. */
#ifndef MAIN
extern guiinfo gui;
#endif

#endif /* GUI_STATE_H */
