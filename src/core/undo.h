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

/* Sentinel value for historic_struct fields that hold a FLIS layer item_id,
 * indicating that the state was saved from a plain FITS image (not a FLIS
 * layer).  Distinguishes "no layer" from item_id 0, which would be invalid
 * but might be confused with "unset" if 0 were the sentinel. */
#define FLIS_UNDO_LAYER_NONE (-1)

gboolean is_undo_available();
gboolean is_redo_available();
int undo_display_data(int dir);
int undo_save_state(fits *fit, const char *message, ...);
int	undo_flush();
gboolean undo_in_thread(gpointer user_data);
gboolean redo_in_thread(gpointer user_data);

/* ----------------------------------------------------------------------- */
/* FLIS undo API — added at stage 1.4.  All functions below save undo      */
/* state for FLIS layer mutations.  They are no-ops when running headless  */
/* (com.script == TRUE); the gate lives in the worker layer in stage 1.5   */
/* but each save function also short-circuits if no single image is loaded.*/
/* All functions return 0 on success, non-zero on failure.                  */
/* ----------------------------------------------------------------------- */

/**
 * flis_undo_purge_layer:
 * @item_id: stable item_id of the layer being removed.
 *
 * Purges all undo/redo entries belonging to the specified layer from both
 * stacks.  For compound (multi-layer) entries, only the matching layer's
 * entry is removed; if that leaves the compound entry empty it is removed
 * altogether.  Call BEFORE the layer is destroyed so item_id is still valid.
 */
void flis_undo_purge_layer(gint item_id);

/**
 * undo_save_flis_layer_props:
 * Save a property-only undo entry for @layer.  No swap file is written;
 * only blend mode, opacity, visibility, lock, tint, position, and name are
 * stored.  Call BEFORE the mutation so the saved state is the pre-change.
 */
int undo_save_flis_layer_props(flis_layer_t *layer, const char *message, ...);

/**
 * undo_save_flis_layer_props_snapshot:
 * Variant of undo_save_flis_layer_props that takes a pre-built snapshot
 * rather than reading the live layer.  Used by the opacity-drag-end path
 * to record the pre-drag state without per-tick undo entries.
 */
int undo_save_flis_layer_props_snapshot(gint item_id,
                                        const flis_layer_props_t *props,
                                        const char *message);

/**
 * undo_save_flis_lmask:
 * Save the current layer-mask state of @layer.  Used for add/remove/move
 * operations.  Empty mask → "no mask" marker; existing mask → swap file.
 * Call BEFORE the mutation.
 */
int undo_save_flis_lmask(flis_layer_t *layer, const char *message);

/**
 * undo_save_flis_lmask_move:
 * Atomic undo entry for a mask move from @source to @dest.  A single
 * entry captures both layer IDs; undo reverses the move in one step.
 * Call BEFORE the move.
 */
int undo_save_flis_lmask_move(flis_layer_t *source, flis_layer_t *dest,
                              const char *message);

/**
 * undo_save_flis_layer_reorder:
 * Atomic undo entry for a layer-order swap of @layer_a and @layer_b.
 * Records both layers' layer_order so undo restores them simultaneously.
 * Call BEFORE the swap.
 */
int undo_save_flis_layer_reorder(flis_layer_t *layer_a, flis_layer_t *layer_b,
                                 const char *message);

/**
 * undo_save_flis_multi_layer:
 * Atomic compound undo entry covering every layer in @layers.  Each layer
 * gets a pixel + pmask + props snapshot stored in one historic entry; a
 * single undo step reverts all layers.  Used for registration and other
 * multi-layer operations.  Returns non-zero on any swap-file failure (in
 * which case nothing is pushed and partial files are cleaned up).
 */
int undo_save_flis_multi_layer(GSList *layers, const char *message, ...);

/**
 * undo_save_flis_multi_layer_props:
 * Compound props-only variant — no pixel swap files.  Used for group
 * drag-moves that update position_x/y on every group member at once.
 */
int undo_save_flis_multi_layer_props(GSList *layers, const char *message, ...);

/**
 * undo_save_processing_mask:
 * Save the processing mask state of @fit only (no pixel data).  Analogous
 * to undo_save_flis_lmask but for fit->mask.  Used by operations that
 * create / modify / delete a processing mask without touching pixels.
 */
int undo_save_processing_mask(fits *fit, const char *message, ...);

/**
 * undo_save_flis_layer_full:
 * Save a complete single-layer state: pixels (from @fit_snapshot), processing
 * mask, layer mask (@lmask_snapshot, or NULL if @lay had none), and the
 * @props_snapshot.  Used for geometry-changing operations (rotate, mirror,
 * crop, binning, resample) so all aspects of the layer are reverted together.
 * @lmask_snapshot ownership is NOT transferred — caller frees after this call.
 */
int undo_save_flis_layer_full(fits *fit_snapshot,
                              flis_layer_t *lay,
                              layermask_t *lmask_snapshot,
                              const flis_layer_props_t *props_snapshot,
                              const char *message, ...);

#endif /* UNDO_H_ */
