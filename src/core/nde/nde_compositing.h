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
#ifndef _NDE_COMPOSITING_H_
#define _NDE_COMPOSITING_H_

/**
 * \file nde_compositing.h
 * \brief Live compositing-state history records (nde-graph-design-note.md §6.1).
 *
 * `layer.set_opacity`, `layer.set_blend` and `layer.set_visible` are NOT pixel
 * operations: they are inputs to the compositor, which reads them from live
 * FLIS layer state at render time.  They are therefore never members of a
 * replay chain (nde_chain_build skips them) — but they ARE editable, and until
 * now editing them did nothing at all.
 *
 * The model here is a FOLD, not a replay: a layer's compositing state is the
 * last value recorded for each property, and amending / deleting / reordering
 * such a record simply re-folds the log onto the layer and re-renders.  No
 * pixels move, no chain is built, nothing is replayed.
 *
 * THE LOG IS AUTHORITATIVE.  The fold starts from the state the layer was in
 * before the log's first word about it and applies every record, so a layer's
 * compositing state is exactly what its history says it is — including when
 * the last record for a property is DELETED, which must return the property
 * to where it started rather than strand it at the deleted record's value.
 * The alternative (leaving properties the log does not mention at their live
 * value) was tried and rejected: it lets the layer sit in a state no part of
 * the log describes, so a later save writes FLIS_META that disagrees with
 * FLIS_HIST.
 *
 * THE STARTING STATE is the new-layer defaults for a layer the log saw
 * created.  A layer that arrived with its own state (a file opened without a
 * history, a duplicate) starts from that state instead: see the block above
 * nde_compositing_params_with_base() below.  Only a layer with neither a
 * recorded base nor a mark — one whose records predate the base_* keys and
 * which was not loaded or duplicated by this build — still folds from the
 * defaults.
 *
 * RESET POINTS.  Two structural records restore an item to the compositing
 * defaults and so make every property "known" from that point on:
 * `document.flatten` (which explicitly resets blend/opacity/visible on the
 * surviving base) and `layer.add` (a fresh layer starts at the defaults).
 * `layer.merge_down` is deliberately NOT a reset — it preserves the surviving
 * bottom layer's blend and opacity — and `layer.duplicate` copies the source
 * layer's properties, so neither makes the defaults knowable.
 *
 * Threading: reads the history through nde_history_snapshot() and mutates
 * live FLIS layer state, so it must run on a thread that owns the processing
 * slot (the edit job) or the main thread.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** TRUE for the three compositing-state op ids. */
gboolean nde_compositing_is_op(const char *op_id);

/**
 * Validate a candidate params blob for @op_id (range- and enum-checked), the
 * compositing counterpart of an op descriptor's deserialize round-trip.
 * FALSE + heap message in @err when the blob lacks the key or is out of range.
 */
gboolean nde_compositing_validate(const char *op_id, const char *params, gchar **err);

/** TRUE for a blend mode a LAYER may carry.  PASS_THROUGH is excluded on
 *  purpose: it is a group mode.  Shared with the composite node, which range-
 *  checks the same value once per recorded input (nde_composite.h). */
gboolean nde_compositing_blend_valid(gint64 v);

/**
 * Re-fold the live history onto layer @item_id's compositing state and
 * invalidate the composite.  A no-op returning TRUE for @item_id < 0 (a plain
 * image has no compositing state).  FALSE + heap @err when the layer is gone.
 */
gboolean nde_compositing_recompute(gint item_id, gchar **err);

/* The bare fold: derive the compositing state the live history describes for
 * @item_id (creation defaults + every compositing record).  Works for items
 * with no live layer (retained composite inputs). */
void nde_compositing_fold(gint item_id, gfloat *out_opacity,
                          gint *out_blend /* flis_blend_mode_t */,
                          gboolean *out_visible,
                          gboolean *out_tinted, double *out_tint /* [3] */);

/* The POSITIONAL fold: the state as of the live log just BEFORE record
 * @upto_record_id (exclusive; 0 = the whole log, i.e. nde_compositing_fold).
 * The joint-record replay (nde_joint.c) reads participants' tints through
 * this so a tint recorded AFTER the joint record does not change what the
 * record recomputes — only amending one recorded before it does. */
void nde_compositing_fold_upto(gint item_id, gint64 upto_record_id,
                               gfloat *out_opacity,
                               gint *out_blend /* flis_blend_mode_t */,
                               gboolean *out_visible,
                               gboolean *out_tinted, double *out_tint /* [3] */);

/* TRUE when @item_id's live history contains a tint record (gates tint
 * folding for documents recorded before tint capture existed). */
gboolean nde_compositing_has_tint_record(gint item_id);

/**
 * The fold starts from the state a layer was in before the log's first word
 * about it.  For a layer the log saw created that is the new-layer defaults;
 * for one that arrived with its own blend mode, opacity or tint (a file
 * opened without a history, a duplicate) it is that state, which the layer
 * is marked with (flis_layer_mark_compositing_base) and which every
 * compositing record of it carries as base_* keys, so that it survives a
 * save and reload.
 *
 * nde_compositing_params_with_base() appends those keys to a new record's
 * @params (taken; the result is returned and owned by the caller) when the
 * starting state of @item_id is known.  nde_compositing_carry_base() returns
 * a copy of @new_params with @old_params' base keys added — an amend replaces
 * a record's params wholesale and must not drop them.
 *
 * nde_compositing_adopt_recorded_bases() is for the loader, once the file's
 * history is attached: a loaded layer is marked with the state it has NOW,
 * which is where it started only if the history never touched it.  Where a
 * record says otherwise, the record is right, and the layer's mark is
 * corrected to match — so the answer outlives the deletion of that record.
 */
void nde_compositing_adopt_recorded_bases(void);
gchar *nde_compositing_params_with_base(gchar *params, gint item_id);
gchar *nde_compositing_carry_base(const char *old_params, const char *new_params);

/* As above, bounded to records strictly before @upto_record_id (0 = whole
 * log).  The joint replay's "does the log say anything about this layer's
 * tint before me?" question — FALSE means fall back to the tint recorded in
 * the joint record's own params. */
gboolean nde_compositing_has_tint_record_upto(gint item_id, gint64 upto_record_id);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_COMPOSITING_H_ */
