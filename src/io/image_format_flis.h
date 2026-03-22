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

/* flis_layer_props_t — lightweight snapshot of all non-pixel layer properties.
 * Used by the undo mechanism for property-only undo states that require no
 * swap file.  See undo_save_flis_layer_props() in undo.c. */
typedef struct flis_layer_props_t {
    flis_blend_mode_t blend_mode;
    gfloat            opacity;
    gboolean          visible;
    gboolean          locked;
    gboolean          has_tint;
    flis_tint_t       tint;
    gchar             name[33];   /* matches max-length in the UI */
} flis_layer_props_t;

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

/* -----------------------------------------------------------------------
 * Layer lookup helpers
 * ----------------------------------------------------------------------- */

void flis_sort_layer_stack(void);

/**
 * flis_active_layer:
 *
 * Returns the flis_layer_t* for the currently active layer in com.uniq,
 * or NULL if no layers are loaded.
 */
flis_layer_t *flis_active_layer(void);

/**
 * flis_active_layer_fit:
 *
 * Returns the fits* for the currently active FLIS layer, or NULL.
 * Use this instead of gfit when you need a stable pointer to the active
 * layer's pixel data that is unaffected by the gfit-to-composite swap
 * that redraw() performs during FLIS rendering.
 */
fits *flis_active_layer_fit(void);

/**
 * flis_update_layer_offset_after_crop:
 * @sel_x: x coordinate of the crop selection top-left (canvas pixels)
 * @sel_y: y coordinate of the crop selection top-left (canvas pixels)
 *
 * Updates FLIS layer offsets after a crop has been applied to gfit.
 * Call this immediately after the crop operation, passing the selection
 * origin that was in effect before the crop (i.e. com.selection.x/y).
 *
 * - Base layer cropped: all other layers' position_x/y are reduced by
 *   (sel_x, sel_y) so they remain correctly placed on the new canvas.
 * - Non-base layer cropped: that layer's position_x/y is set to (sel_x, sel_y).
 *
 * Invalidates the composite. No-op if no FLIS is loaded.
 */
void flis_update_layer_offset_after_crop(gint sel_x, gint sel_y);

/**
 * flis_update_layer_offset_after_resize:
 * @old_rx: canvas width before resize
 * @old_ry: canvas height before resize
 * @new_rx: canvas width after resize
 * @new_ry: canvas height after resize
 *
 * Updates FLIS layer offsets after a resize has been applied to gfit.
 *
 * - Base layer resized: all non-base layers' centres are scaled proportionally.
 * - Non-base layer resized: that layer's centre on the canvas is preserved.
 *
 * Invalidates the composite. No-op if no FLIS is loaded.
 */
void flis_update_layer_offset_after_resize(gint old_rx, gint old_ry,
                                           gint new_rx, gint new_ry);

/**
 * flis_update_layer_offset_after_rotate:
 * @old_rx: canvas width before rotation
 * @old_ry: canvas height before rotation
 * @new_rx: canvas width after rotation
 * @new_ry: canvas height after rotation
 * @angle:  rotation angle in degrees (positive = CCW in Siril convention)
 *
 * Updates FLIS layer offsets after a rotation has been applied to gfit.
 *
 * - Base layer rotated: all non-base layers' centres are rotated by the same
 *   angle around the canvas centre and mapped into the new canvas.
 * - Non-base layer rotated: that layer's centre on the canvas is preserved.
 *
 * Invalidates the composite. No-op if no FLIS is loaded.
 */
void flis_update_layer_offset_after_rotate(gint old_rx, gint old_ry,
                                           gint new_rx, gint new_ry,
                                           double angle);

/**
 * flis_promote_from_gfit:
 * @name: name for the base layer, or NULL to use "Background".
 *
 * Converts the currently loaded plain FITS image (com.uniq->layers == NULL)
 * into a single-layer FLIS in memory by wrapping gfit as the base layer.
 * After this call is_current_image_flis() returns TRUE and additional
 * layers can be added with flis_layer_add().
 *
 * gfit must be heap-allocated (i.e. obtained via malloc/calloc, not a
 * reference to a static struct) because flis_layer_free() will free it
 * when the layer is destroyed.
 *
 * Returns: 0 on success, non-zero on error.
 */
int flis_promote_from_gfit(const gchar *name);

/**
 * flis_layer_get_by_id:
 * @item_id: the stable ITEM_ID to look up.
 *
 * Returns the flis_layer_t* with the given item_id, or NULL if not found.
 */
flis_layer_t *flis_layer_get_by_id(gint item_id);

/**
 * flis_layer_get_index:
 * @layer: layer pointer to locate.
 *
 * Returns the 0-based index of @layer in com.uniq->layers, or -1 if not
 * found.
 */
gint flis_layer_get_index(const flis_layer_t *layer);

/**
 * flis_layer_count:
 *
 * Returns the number of layers in com.uniq->layers, or 0 if com.uniq is
 * NULL.
 */
gint flis_layer_count(void);

/**
 * is_current_image_flis:
 *
 * Returns TRUE if a FLIS file is currently loaded (com.uniq->layers is
 * non-NULL).  Returns FALSE for a regular FITS image or when no image is
 * loaded.
 *
 * This is the canonical test for "are we in FLIS mode?".  Use it
 * everywhere rather than inspecting com.uniq directly, so the definition
 * stays in one place.  In particular:
 *   - image_display.c uses it to decide whether to build a composite
 *   - UI code uses it to enable/disable layer-related widgets and menus
 *   - Tool and processing code uses it to refuse operations that are not
 *     meaningful on a layered image (e.g. direct gfit pixel modification)
 */
gboolean is_current_image_flis(void);

/* -----------------------------------------------------------------------
 * Stack management
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_add:
 * @fit:  Pixel data for the new layer.  Ownership is transferred.
 * @name: Layer name; copied internally.  May be NULL (defaults to "Layer").
 *
 * Creates a new flis_layer_t, assigns the next available item_id, places
 * it above all existing layers (highest layer_order + FLIS_ORDER_STEP),
 * appends it to com.uniq->layers (which remains sorted), and sets it as
 * the active layer.
 *
 * Returns: pointer to the newly created layer, or NULL on failure.
 */
flis_layer_t *flis_layer_add(fits *fit, const gchar *name);

/**
 * flis_layer_remove:
 * @layer: layer to remove and free.
 *
 * Removes @layer from com.uniq->layers and frees it.  If the removed layer
 * was the active layer, the active layer is moved to the layer below it
 * (or above if it was the base layer).  Fails if only one layer remains or
 * if @layer is locked.
 *
 * Returns: 0 on success, non-zero on error.
 */
int flis_layer_remove(flis_layer_t *layer);

/**
 * flis_layer_duplicate:
 * @layer: layer to duplicate.
 *
 * Creates a deep copy of @layer (pixel data, layer mask, processing mask,
 * all metadata).  The duplicate receives a new item_id and a layer_order
 * one step above @layer.  The duplicate is inserted into com.uniq->layers
 * and set as the active layer.
 *
 * Returns: pointer to the duplicate layer, or NULL on failure.
 */
flis_layer_t *flis_layer_duplicate(const flis_layer_t *layer);

/**
 * flis_layer_move_up:
 * @layer: layer to move.
 *
 * Increases this layer's stacking position by one step (swaps layer_order
 * with the next-higher neighbour).  "Up" means rendered on top of more
 * layers.  No-op if @layer is already at the top of the stack.  Fails if
 * @layer is locked.
 *
 * Returns: 0 on success, 1 if already at top, negative on error.
 */
int flis_layer_move_up(flis_layer_t *layer);

/**
 * flis_layer_move_down:
 * @layer: layer to move.
 *
 * Decreases this layer's stacking position by one step (swaps layer_order
 * with the next-lower neighbour).  No-op if @layer is the base layer.
 * Fails if @layer is locked.
 *
 * Returns: 0 on success, 1 if already at bottom, negative on error.
 */
int flis_layer_move_down(flis_layer_t *layer);

/* -----------------------------------------------------------------------
 * Compositing properties
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_set_blend_mode:
 * @layer: target layer.
 * @mode:  new blend mode.
 *
 * Returns: 0 on success, non-zero if @layer is locked.
 */
int flis_layer_set_blend_mode(flis_layer_t *layer, flis_blend_mode_t mode);

/**
 * flis_layer_set_opacity:
 * @layer:   target layer.
 * @opacity: new opacity, clamped to [0.0, 1.0].
 *
 * Returns: 0 on success, non-zero if @layer is locked.
 */
int flis_layer_set_opacity(flis_layer_t *layer, gfloat opacity);

/**
 * flis_layer_set_visible:
 * @layer:   target layer.
 * @visible: new visibility state.
 *
 * Visibility may be toggled even on a locked layer (it is a display
 * property, not a data property).
 *
 * Returns: 0 on success.
 */
int flis_layer_set_visible(flis_layer_t *layer, gboolean visible);

/* -----------------------------------------------------------------------
 * Layer name and lock
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_set_name:
 * @layer: target layer.
 * @name:  new name; copied internally.
 *
 * Returns: 0 on success, non-zero if @layer is locked.
 */
int flis_layer_set_name(flis_layer_t *layer, const gchar *name);

/**
 * flis_layer_set_locked:
 * @layer:  target layer.
 * @locked: TRUE to lock, FALSE to unlock.
 *
 * Lock state may always be changed regardless of current lock state.
 *
 * Returns: 0 on success.
 */
int flis_layer_set_locked(flis_layer_t *layer, gboolean locked);

/* -----------------------------------------------------------------------
 * Layer mask (LMASK) operations
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_set_lmask:
 * @layer: target layer.
 * @lmask: new layer mask.  Ownership is transferred.  Must have the same
 *         spatial dimensions as @layer->fit.  If NULL, any existing layer
 *         mask is removed (equivalent to flis_layer_remove_lmask()).
 *
 * Any previously attached layer mask is freed before the new one is set.
 *
 * Returns: 0 on success, non-zero if @layer is locked or dimensions
 *          do not match.
 */
int flis_layer_set_lmask(flis_layer_t *layer, layermask_t *lmask);

/**
 * flis_layer_remove_lmask:
 * @layer: target layer.
 *
 * Removes and frees the layer mask attached to @layer, if any.
 *
 * Returns: 0 on success, non-zero if @layer is locked.
 */
int flis_layer_remove_lmask(flis_layer_t *layer);

/**
 * flis_layer_move_lmask:
 * @from: source layer whose lmask will be detached.
 * @to:   destination layer that will receive the lmask.
 *
 * Transfers the layer mask from @from to @to.  The transfer is refused if:
 *   - @from has no layer mask
 *   - @from or @to is locked
 *   - the mask dimensions do not match @to->fit's spatial dimensions
 * Any existing lmask on @to is freed before the transfer.
 *
 * Returns: 0 on success, non-zero on any of the above errors.
 */
int flis_layer_move_lmask(flis_layer_t *from, flis_layer_t *to);

/* -----------------------------------------------------------------------
 * Mono layer tinting
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_set_tint:
 * @layer: target layer.  Must be a MONO layer (fit->naxes[2] == 1).
 * @r, @g, @b: tint colour components, normalised to [0.0, 1.0].
 *
 * Sets the LAYER_COLOR tint applied to this mono layer when compositing
 * into an RGB stack (FLIS spec Section 6.5).
 *
 * Returns: 0 on success, non-zero if @layer is locked or is not MONO.
 */
int flis_layer_set_tint(flis_layer_t *layer, double r, double g, double b);

/**
 * flis_layer_clear_tint:
 * @layer: target layer.
 *
 * Removes the tint from a mono layer, reverting to neutral broadcast
 * (all three RGB channels receive equal weight).
 *
 * Returns: 0 on success, non-zero if @layer is locked.
 */
int flis_layer_clear_tint(flis_layer_t *layer);

/* -----------------------------------------------------------------------
 * Timestamp utility (also useful for callers updating layer data
 * programmatically)
 * ----------------------------------------------------------------------- */

/**
 * flis_layer_touch_modified:
 * @layer: layer to timestamp.
 *
 * Sets layer->modified to the current UTC time in ISO 8601 format.
 * Called internally by all mutating layer functions.
 */
void flis_layer_touch_modified(flis_layer_t *layer);

/* -----------------------------------------------------------------------
 * Convenience: layer_order gap between adjacent layers on creation.
 * Gaps allow insertion without renumbering existing layers.
 * ----------------------------------------------------------------------- */
#define FLIS_ORDER_STEP 10

#endif /* IMAGE_FORMAT_FLIS_H */
