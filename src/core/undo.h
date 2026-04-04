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
#ifndef UNDO_H_
#define UNDO_H_

#define UNDO		-1
#define REDO		 1
#define HISTORY_SIZE 30		// size of the operations history

/* Sentinel value for historic.flis_layer_id indicating that the state was
 * saved from a plain FITS image (not a FLIS layer). */
#define FLIS_UNDO_LAYER_NONE  (-1)

/* Forward declarations for FLIS types used in function signatures below.
 * Full definitions are in io/image_format_flis.h.
 * The tag on flis_layer_props_t must match the tagged struct definition there.
 * siril.h also forward-declares flis_layer_props_t; both declarations are
 * compatible because they name the same struct tag. */
typedef struct _flis_layer_t      flis_layer_t;
typedef struct flis_layer_props_t flis_layer_props_t;

/*
 * IMPORTANT — siril.h change required
 * =====================================
 * The `historic` struct in siril.h must have the following fields added:
 *
 *   gint               flis_layer_id;   (already added)
 *   flis_layer_props_t *layer_props;    (already added)
 *   gchar             *lmask_filename;  (already added)
 *   gint               lmask_layer_id;  (already added)
 *   size_t             lmask_w;         (already added)
 *   size_t             lmask_h;         (already added)
 *   guint8             lmask_bitpix;    (already added)
 *   gint               lmask_dest_layer_id; (already added)
 *   gint               reorder_layer_a_id;  (already added)
 *   gint               reorder_layer_a_order; (already added)
 *   gint               reorder_layer_b_id;  (already added)
 *   gint               reorder_layer_b_order; (already added)
 *   gint               flis_position_x; ← ADD THESE
 *   gint               flis_position_y;
 *
 * `layer_props` is non-NULL only for property-only undo states (no swap
 * file).  When non-NULL and `filename` is NULL, undo_get_data() restores
 * only the layer properties rather than pixel data.
 *
 * Suggested placement immediately after flis_layer_id:
 *
 *   gint               flis_layer_id;
 *   flis_layer_props_t *layer_props;    ← ADD THIS
 */

gboolean is_undo_available(void);
gboolean is_redo_available(void);
int      undo_display_data(int dir);
int      undo_save_state(fits *fit, const char *message, ...);
int      undo_flush(void);
void     set_undo_redo_tooltip(void);
gboolean undo_in_thread(gpointer user_data);
gboolean redo_in_thread(gpointer user_data);

/**
 * flis_undo_purge_layer:
 * @item_id: the stable item_id of the layer being removed.
 *
 * Purges all undo/redo states belonging to the specified layer from the
 * history ring buffer and frees their swap files.  All states for other
 * layers are preserved.
 *
 * Call this from flis_gui.c (on_flis_remove_layer_clicked) immediately
 * BEFORE the layer is destroyed, so the item_id is still valid.
 *
 * Add, duplicate, and reorder operations do NOT need to call this:
 *   - Add / duplicate: the new layer has no prior states.
 *   - Reorder: layer lookup uses item_id not position, so reordering is
 *     transparent to the undo mechanism.
 */
void flis_undo_purge_layer(gint item_id);

/**
 * undo_save_flis_layer_props:
 * @layer:   the FLIS layer whose properties are to be snapshotted.
 * @message: undo label shown in the tooltip (printf-style format string).
 *
 * Saves a lightweight property-only undo state for the given layer.
 * No pixel swap file is written; only the layer's blend mode, opacity,
 * visibility, lock, tint, and name are stored.
 *
 * Call this from property-change handlers in flis_gui.c BEFORE applying
 * the change, so the saved state represents the pre-change values.
 *
 * Returns 0 on success, non-zero on failure.
 */
int undo_save_flis_layer_props(flis_layer_t *layer, const char *message, ...);

/**
 * undo_save_flis_layer_props_snapshot:
 * @item_id:  the layer's stable item_id.
 * @props:    pre-built snapshot of the layer properties at the desired point.
 * @message:  undo label.
 *
 * Variant of undo_save_flis_layer_props() that takes an already-built
 * flis_layer_props_t snapshot rather than reading from the live layer.
 * Used by the opacity drag-end path to record the pre-drag state after
 * the drag is complete, without polluting history with every tick value.
 */
int undo_save_flis_layer_props_snapshot(gint item_id,
                                        const flis_layer_props_t *props,
                                        const char *message);

/**
 * undo_save_flis_lmask:
 * @layer:   the FLIS layer whose lmask state is to be snapshotted.
 * @message: undo label.
 *
 * Saves the current layer-mask state of @layer as an undo entry.
 * Call this BEFORE any change to layer->lmask:
 *   - Before adding a mask:   saves "no mask" marker; undo removes the mask.
 *   - Before removing a mask: saves mask pixels to swap; undo restores them.
 *   - Before moving a mask:   call on both source and destination layers.
 *
 * Returns 0 on success, non-zero on failure.
 */
int undo_save_flis_lmask(flis_layer_t *layer, const char *message);

/**
 * undo_save_flis_lmask_move:
 * @source: the layer currently holding the mask (before the move).
 * @dest:   the layer that will receive the mask.
 * @message: undo label.
 *
 * Saves an atomic undo state for a mask move.  A single ring-buffer entry
 * records both layer IDs; undo reverses the move in one step with no broken
 * intermediate state.  Call BEFORE flis_layer_move_lmask(source, dest).
 */
int undo_save_flis_lmask_move(flis_layer_t *source, flis_layer_t *dest,
                               const char *message);

/**
 * undo_save_flis_layer_reorder:
 * @layer_a: the layer being moved (flis_selected).
 * @layer_b: the layer it will swap order with.
 * @message: undo label.
 *
 * Saves an atomic undo state for a layer reorder (move-up / move-down).
 * Records the current layer_order of both layers before the swap.
 * Call BEFORE flis_layer_move_up() or flis_layer_move_down().
 */
int undo_save_flis_layer_reorder(flis_layer_t *layer_a, flis_layer_t *layer_b,
                                  const char *message);

/**
 * undo_save_flis_multi_layer:
 * @layers:  GSList of flis_layer_t* to snapshot (all layers affected by the
 *           operation).  The list is not modified or freed.
 * @message: undo label (printf-style format string).
 *
 * Saves an atomic compound undo state covering every layer in @layers.
 * For each layer a pixel swap file is written (and a processing-mask swap
 * file if the layer carries a processing mask).  A full flis_layer_props_t
 * snapshot (blend mode, opacity, visibility, lock, tint, name, position) is
 * also stored per layer.
 *
 * Use this for multi-layer operations such as registration or layer matching
 * where pixel data, WCS, positions, and/or properties of several layers
 * change atomically.  A single undo step reverts all layers at once.
 *
 * Call BEFORE the operation runs, so the saved state holds the pre-operation
 * values.
 *
 * Returns 0 on success, non-zero if any swap file could not be written
 * (in which case no undo state is pushed and any already-written swap files
 * are cleaned up).
 */
int undo_save_flis_multi_layer(GSList *layers, const char *message, ...);

#endif /* UNDO_H_ */
