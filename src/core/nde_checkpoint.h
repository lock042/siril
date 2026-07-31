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
#ifndef _NDE_CHECKPOINT_H_
#define _NDE_CHECKPOINT_H_

/**
 * \file nde_checkpoint.h
 * \brief In-session baseline pixel-snapshot store for nondestructive editing.
 *
 * Phase 2 (replay) rebuilds a document's pixels from a per-item BASELINE
 * checkpoint through the recorded op chain.  This module keeps that baseline:
 * the pre-op pixels of the FIRST operation captured against a layer item_id
 * (−1 = plain single image), snapshotted to a delete-on-close swap file so
 * the working set stays out of RAM, exactly like the undo swap cache.
 *
 * Keyed by layer item_id.  A baseline is taken ONCE, at the first capture
 * touching that item — later ensure() calls are cheap no-ops, so the stored
 * pixels always reflect the pre-FIRST-op state (never a later intermediate).
 *
 * Threading: guarded by an internal GMutex that is a strict LEAF in the lock
 * order — it acquires no other lock and never blocks while held.  Fits copies
 * and swap-file I/O are performed OUTSIDE the lock (the registry only holds
 * swap-file handles, never pixel work).  It must NEVER be held together with
 * the nde history mutex.  All public functions are callable from any thread.
 *
 * Lifetime coupling (flis-nde-sketch.md §7 decision 13): the swap-file
 * storage is deliberately FACTORED (a clean write/read pair) so it can later
 * become the shared refcounted snapshot store that undo and NDE both draw on.
 * Nothing here may assume checkpoints and undo entries have independent
 * lifetimes.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The image type is `fits` (typedef of `struct ffit` in siril.h).  Forward-
 * declare the tag so this header stays light; the includers all pull in
 * siril.h for the full definition. */
struct ffit;
typedef struct ffit fits;

/**
 * Snapshot @pre (PRE-op pixels) as the baseline for @item_id if none is
 * recorded yet; a cheap no-op when one already exists.  The fits is
 * deep-copied to a delete-on-close swap file OUTSIDE the registry lock, so
 * this may be called from a worker thread while other locks are held.
 * @pre must be a fully-populated single image (DATA_USHORT or DATA_FLOAT).
 * A NULL @pre or a fits with no pixel buffer is ignored.
 */
void nde_checkpoint_baseline_ensure(const fits *pre, gint item_id);

/**
 * Load the baseline for @item_id into a freshly allocated, fully-owned fits
 * (caller clearfits()+frees it), or NULL when none is recorded.  Pixel-exact
 * with the buffer passed to ensure()/adopt() — baselines round-trip losslessly.
 */
fits *nde_checkpoint_baseline_get(gint item_id);

/** TRUE when a baseline is recorded for @item_id. */
gboolean nde_checkpoint_baseline_exists(gint item_id);

/**
 * Adopt @src as the baseline for @item_id, replacing any existing one.  Used
 * by the FLIS loader to install persisted NDE_BASE HDUs.  Deep-copies @src
 * to a swap file; @src is NOT taken over (caller keeps ownership).  Unlike
 * ensure(), adopt overwrites — the on-disk baseline is authoritative on load.
 */
void nde_checkpoint_baseline_adopt(const fits *src, gint item_id);

/** Drop the baseline AND output checkpoints for @item_id (layer removed). */
void nde_checkpoint_drop(gint item_id);

/** Drop every checkpoint (document closed / new document loaded). */
void nde_checkpoint_purge(void);

/* ---- output checkpoints (phase-4 barriers, nde-phase4-plan.md) ----------
 * A barrier record's POST-op pixels, keyed by record_id: the restart point
 * that lets everything after the barrier stay editable.  Same lossless /
 * leaf-lock discipline as baselines; a checkpoint dies with its record
 * (delete, dead-tail truncation), its layer (drop) or the document
 * (purge).                                                                */

/** Store @post as the output checkpoint of @record_id (overwrites). */
void nde_checkpoint_output_store(const fits *post, gint64 record_id, gint item_id);

/**
 * Adopt @src as the output checkpoint of @record_id (owning layer @item_id),
 * replacing any existing one.  Used by the FLIS loader to install persisted
 * NDE_CKPT HDUs.  Deep-copies @src to a swap file; @src is NOT taken over
 * (caller keeps ownership).  Mirrors nde_checkpoint_baseline_adopt.
 */
void nde_checkpoint_output_adopt(const fits *src, gint64 record_id, gint item_id);

/** Load a copy of @record_id's output checkpoint (caller clears+frees). */
fits *nde_checkpoint_output_get(gint64 record_id);

gboolean nde_checkpoint_output_exists(gint64 record_id);

/** Drop one output checkpoint (its record was deleted or truncated). */
void nde_checkpoint_output_drop(gint64 record_id);

/* ---- pin coordinates ----------------------------------------------------
 * An input pin names {item, record}, where record 0 means "that item's
 * baseline" and anything else means that record's output (nde_history.h).
 * These three wrap the choice so the convention lives in one place rather
 * than at every pin-resolving site.                                        */

void     nde_checkpoint_store_at(const fits *f, gint item_id, gint64 record_id);
fits    *nde_checkpoint_get_at(gint item_id, gint64 record_id);
gboolean nde_checkpoint_exists_at(gint item_id, gint64 record_id);

/* ---- layer offsets (graph step 5) --------------------------------------
 * A checkpoint stores pixels, but the replay unit for a FLIS layer is not
 * pixels alone: geometry operations (crop, resample, rotate, mirror, binning)
 * move the layer on the canvas as well as changing its pixels, and the
 * offset helpers evolve position_x/position_y from their previous value.  So
 * a chain that contains one cannot be replayed from pixels alone — the
 * restart point has to say where the layer WAS.
 *
 * The offset is attached beside the snapshot rather than inside it, so every
 * existing call site keeps its signature and a checkpoint with no offset
 * recorded (a plain image, or a file written before this) simply reports
 * none.  It persists on the NDE_BASE / NDE_CKPT HDU as FLIS_POSX/FLIS_POSY.
 */

/** Record the layer offset that goes with @item_id's baseline.  First value
 *  wins, exactly as the baseline itself does — the baseline is the pre-FIRST-op
 *  state, so its position is the one the layer had then.  Ignored when no
 *  baseline exists for @item_id. */
void nde_checkpoint_baseline_set_offset(gint item_id, gint pos_x, gint pos_y);

/** FALSE when no offset is recorded (plain image, or an older file). */
gboolean nde_checkpoint_baseline_get_offset(gint item_id, gint *pos_x, gint *pos_y);

/** Same, for a record's output checkpoint. */
void     nde_checkpoint_output_set_offset(gint64 record_id, gint pos_x, gint pos_y);
gboolean nde_checkpoint_output_get_offset(gint64 record_id, gint *pos_x, gint *pos_y);

/* Plain → FLIS promote: move the baseline and every output checkpoint
 * tagged for @from_item onto @to_item (see flis_promote_from_gfit). */
void nde_checkpoint_rebind_item(gint from_item, gint to_item);

/* ---- retention (design note §7; policy lives in nde_retention.c) --------
 * Checkpoints are EDITABILITY PINS, not cache: dropping one does not cost
 * time, it removes something the user could do.  So they are evicted only
 * after the cache pool is exhausted, and never silently.
 */

/** Payload bytes held by the baseline and output-checkpoint tables. */
gint64 nde_checkpoint_bytes(void);



/**
 * Active FLIS layer item_id for baseline targeting, or −1 for a plain image.
 * The single idiom shared by every capture site (is_current_image_flis() ?
 * active layer item_id : −1).
 */
gint nde_checkpoint_active_item_id(void);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_CHECKPOINT_H_ */
