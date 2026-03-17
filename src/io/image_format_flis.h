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

#ifndef IMAGE_FORMAT_FLIS_H
#define IMAGE_FORMAT_FLIS_H

#include <stdint.h>
#include <glib.h>
#include "core/siril.h"   /* fits, single, gfit, com, WORD, DATA_FLOAT etc. */

/* -----------------------------------------------------------------------
 * FLIS specification version written/expected by this implementation.
 * ----------------------------------------------------------------------- */
#define FLIS_VERSION_STRING  "1.0"
#define FLIS_THUMBNAIL_MAX   512   /* max pixels on longest side */

/* -----------------------------------------------------------------------
 * Blend modes - values are FLISBLND bitmask positions (Section 5).
 * ----------------------------------------------------------------------- */
typedef enum {
    FLIS_BLEND_NORMAL      = 1,
    FLIS_BLEND_MULTIPLY    = 2,
    FLIS_BLEND_SCREEN      = 4,
    FLIS_BLEND_OVERLAY     = 8,
    FLIS_BLEND_SOFT_LIGHT  = 16,
    FLIS_BLEND_HARD_LIGHT  = 32,
    FLIS_BLEND_COLOR_DODGE = 64,
    FLIS_BLEND_COLOR_BURN  = 128,
    FLIS_BLEND_DARKEN      = 256,
    FLIS_BLEND_LIGHTEN     = 512,
    FLIS_BLEND_DIFFERENCE  = 1024,
    FLIS_BLEND_EXCLUSION   = 2048,
    FLIS_BLEND_HUE         = 4096,
    FLIS_BLEND_SATURATION  = 8192,
    FLIS_BLEND_COLOR       = 16384,
    FLIS_BLEND_LUMINOSITY  = 32768,
} flis_blend_mode_t;

#define FLIS_BLND_ALL 65535  /* bitmask declaring all 16 modes supported */

/* -----------------------------------------------------------------------
 * Tint colour for MONO layers composited into an RGB stack.
 * Components are normalised linear-light values [0.0, 1.0].
 * Default (no tint) is {1.0, 1.0, 1.0} — neutral broadcast.
 * Example: Hα red:       {1.0, 0.2, 0.1}
 *          Oiii blue-green: {0.2, 0.8, 1.0}
 * See FLIS spec Section 6.5 and LAYER_COLOR metadata key.
 * ----------------------------------------------------------------------- */
typedef struct {
    double r, g, b;
} flis_tint_t;

/* -----------------------------------------------------------------------
 * Layer mask — lightweight structure for a greyscale opacity or processing
 * mask.  Does not carry FITS headers or science metadata; just the pixel
 * array and its geometry.
 *
 * bitpix: element size in bits (always 8 for LMASK; 32 for MASK because
 *         Siril's existing processing mask pipeline uses 32-bit float).
 * data:   row-major, bottom-to-top (FITS convention) pixel array.
 *         For bitpix=8  the element type is uint8_t.
 *         For bitpix=32 the element type is float.
 * ----------------------------------------------------------------------- */
typedef struct _layermask_t {
    size_t   w;
    size_t   h;
    uint8_t  bitpix;   /* 8 (layer mask) or 32 (processing mask) */
    void    *data;
} layermask_t;

/* -----------------------------------------------------------------------
 * Per-layer data.  The GSList in single carries flis_layer_t* entries
 * sorted ascending by layer_order.
 *
 * Ownership:
 *   fit        — owned; free with clearfits() then g_free() / free()
 *   lmask      — owned; NULL if no layer mask; free with layermask_free()
 *   Processing mask lives inside fit->mask (existing mechanism).
 *   layer_name — owned; free with g_free()
 *   created / modified — owned; NULL if unset; free with g_free()
 * ----------------------------------------------------------------------- */
typedef struct _flis_layer_t {
    /* Identity and ordering */
    gint              item_id;       /* Stable ITEM_ID from FLIS_META;
                                        never reused within a file.      */
    gint              layer_order;   /* Stacking z-order (ascending;
                                        higher = rendered on top).        */
    gchar            *layer_name;    /* User-visible name (≤32 chars for
                                        on-disk storage).                 */

    /* Image and mask data */
    fits             *fit;           /* Layer pixel data.  gfit == this->fit
                                        for the active layer.             */
    layermask_t      *lmask;         /* Layer mask (FLIS LMASK type):
                                        greyscale, same spatial size as
                                        fit.  NULL if not present.        */
    /* Processing mask (FLIS MASK type) lives in fit->mask — unchanged.  */

    /* Compositing parameters */
    flis_blend_mode_t blend_mode;    /* How this layer blends into the
                                        composite below it.               */
    gfloat            opacity;       /* 0.0 (transparent) – 1.0 (opaque)*/
    gboolean          visible;       /* FALSE: layer skipped in composite*/

    /* Mono-to-RGB tinting (FLIS Section 6.5 / LAYER_COLOR metadata key) */
    gboolean          has_tint;      /* TRUE: layer_tint overrides neutral
                                        broadcast.  Only meaningful when
                                        fit->naxes[2] == 1.              */
    flis_tint_t       layer_tint;    /* Tint applied to mono pixel value
                                        before RGB compositing.          */

    /* Sparse layer support — reserved, must be 0 until implemented */
    gint              position_x;    /* X offset relative to canvas.     */
    gint              position_y;    /* Y offset relative to canvas.     */

    /* Provenance and lock state */
    gchar            *created;       /* ISO 8601 creation timestamp or NULL */
    gchar            *modified;      /* ISO 8601 modification timestamp or NULL */
    gboolean          locked;        /* TRUE: layer locked against edits  */
} flis_layer_t;

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

/**
 * load_flis:
 * @filename: path to the .flis file to open.
 *
 * Reads a FLIS file and populates com.uniq with the layer stack.
 * On success, com.uniq->layers is populated, com.uniq->active_layer is set
 * to the index of the base layer (lowest layer_order), com.uniq->fit and
 * gfit point to the base layer's fits*, and com.uniq->chans reflects its
 * channel count.
 *
 * Returns: 0 on success, non-zero on error.
 */
int load_flis(const gchar *filename);

/**
 * save_flis:
 * @filename: output path (a .flis extension will be appended if absent).
 *
 * Writes all layers in com.uniq->layers to a new FLIS file.
 * The composite thumbnail in HDU 0 is generated from the base layer.
 * Compression follows the same CFITSIO settings as standard FITS saves.
 *
 * Returns: 0 on success, non-zero on error.
 */
int save_flis(const gchar *filename);

/**
 * flis_layer_new:
 * @fit:  Pixel data.  Ownership is transferred to the new layer.
 * @name: Human-readable layer name; copied internally.
 *
 * Allocates a flis_layer_t with sensible defaults (Normal blend, fully
 * opaque, visible, no masks, no tint, position 0,0).
 * The caller must assign item_id and layer_order explicitly.
 *
 * Returns: heap-allocated flis_layer_t*, or NULL on allocation failure.
 */
flis_layer_t *flis_layer_new(fits *fit, const gchar *name);

/**
 * flis_layer_free:
 * @layer: layer to free.  May be NULL (no-op).
 *
 * Frees the flis_layer_t and all owned resources:
 *   fit (via clearfits + free), lmask (via layermask_free),
 *   layer_name, created, modified.
 * Does NOT close the FITS file descriptor — that must be done before
 * handing the fit to this layer, or by the caller beforehand.
 */
void flis_layer_free(flis_layer_t *layer);

/**
 * layermask_free:
 * @mask: mask to free.  May be NULL (no-op).
 */
void layermask_free(layermask_t *mask);

/**
 * flis_free_layers:
 * @uniq: the single struct whose layers list is to be freed.
 *
 * Frees every flis_layer_t in uniq->layers and sets the list to NULL.
 * Also clears uniq->fit and uniq->chans.  Does not free uniq itself.
 */
void flis_free_layers(single *uniq);

/**
 * uniq_set_active_layer:
 * @uniq:  the single struct to update.
 * @index: 0-based index into uniq->layers.
 *
 * Updates uniq->active_layer, uniq->fit, uniq->chans, and gfit to
 * reflect the newly active layer.  Safe to call whenever the active
 * layer changes (e.g., user switches layer in the UI).
 */
void uniq_set_active_layer(single *uniq, gint index);

#endif /* IMAGE_FORMAT_FLIS_H */
