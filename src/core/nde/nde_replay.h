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
#ifndef _NDE_REPLAY_H_
#define _NDE_REPLAY_H_

/**
 * \file nde_replay.h
 * \brief NDE replay engine (phase 2, nde-phase2-3-plan.md P2.D).
 *
 * Rebuilds an item's pixels headlessly: baseline checkpoint → each Tier-A
 * record deserialized and applied through generic_image_worker with the
 * nde_replay flag, always into a PRIVATE fits (never gfit).  Phase 2 uses
 * it for verification (flis_replay_check); phase 3's amend-and-replay
 * feeds an edited chain through the same driver.
 *
 * Threading: nde_chain_replay() MUST run inside a processing job (the
 * caller owns the slot — one continuous occupancy for the whole chain, so
 * python cannot interleave).  nde_chain_build() only snapshots and
 * validates; it is safe anywhere.
 *
 * REGIONS.  A replay that REBUILDS an item is never region-scoped: it applies
 * whole records to a private full-size fits, and generic_image_worker's region
 * path is gated on for_preview, which such a replay never sets.  What IS
 * region-scoped is the amend PREVIEW's tail — see the nde_region_tail block
 * further down, and the ROI model above roi_t (core/siril.h) for the scope
 * rules.  The classification rule both ends share:
 *
 *   an op is region-replayable iff its record is PARAMETER-COMPLETE — it pins
 *   the derived quantity, not the recipe for deriving it — and its spatial
 *   support is bounded and declared (OP_ROI_CAPABLE + roi_halo).
 *
 * Many "global" ops are MORE region-safe on replay than interactively, because
 * the record carries the computed result: an MTF record carries its transfer
 * function, so at replay time it is pixel-local.  Background extraction is the
 * counterexample — its replay_pre reinstalls recorded sample POSITIONS and
 * refits from the image, so on a cropped image it would fit the wrong pixels.
 * It is a barrier until its record pins the model.
 */

#include <glib.h>

#include "core/settings.h"   /* rectangle */

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;
struct op_descriptor;

/* Per-member flag bits in nde_chain.member_flags. */
#define NDE_CHAIN_MEMBER_BARRIER 0x1  /* non-replayable member (Tier B / mask) */

typedef struct {
	gint       item_id;    /* layer item_id; -1 = plain single image */
	GPtrArray *records;    /* chain members: nde_record* views into ->snapshot */
	GPtrArray *snapshot;   /* owns the records (live-prefix deep copy) */
	gboolean   replayable; /* the WHOLE chain replays from the baseline */
	GPtrArray *reasons;    /* gchar*: hard blockers (checkpoint-less barriers,
	                        * document-wide records, missing baseline...) */

	/* Phase-4 barrier model (nde-phase4-plan.md).  A barrier is a member
	 * that cannot be replayed; with an output checkpoint it is a restart
	 * point instead of a hard blocker.  records[tail_start..] is the
	 * editable tail (everything strictly after the last freeze cause). */
	gboolean   is_mask;         /* the item is a MASK: its members are mask ops
	                             * and its value is a mask, so it replays
	                             * through the mask hooks and starts from the
	                             * image its first member read, not from a
	                             * baseline of its own */
	gboolean   has_geometry;    /* a member moves the layer on the canvas, so
	                             * the replay must carry the layer's offset and
	                             * the restart point must supply one */
	gboolean   has_composite;   /* a member is a composite node (graph step 7):
	                             * the chain consumes ANOTHER item's pixels, so
	                             * replaying it replays that item's chain too */
	gboolean   from_composite;  /* the FIRST member is that composite: the item
	                             * was BORN of a merge or flatten and has no
	                             * baseline of its own, the way a mask has none.
	                             * Its origin is the composite's inputs, and
	                             * replay starts by rendering them */
	GArray    *member_flags;    /* guint8 per member (NDE_CHAIN_MEMBER_*) */
	guint      tail_start;      /* first tail member index (== len ⇒ empty tail) */
	gboolean   tail_replayable; /* the tail applies cleanly from its restart point */
	gint64     restart_ckpt_id; /* record id of the restart checkpoint; 0 = baseline */
	guint      first_block;     /* member index of the first HARD blocker (the
	                             * position the first ->reasons entry describes):
	                             * members [0..first_block) replay cleanly from
	                             * the baseline or from barrier checkpoints among
	                             * them.  G_MAXUINT = no hard blocker.  Lets a
	                             * PREFIX consumer (an input pin) succeed when
	                             * only the far end of the chain is broken. */
} nde_chain;

/**
 * Snapshot the live history and extract + validate the replay chain for
 * @item_id.  Membership: LAYER-scope records targeting the item; CANVAS-
 * scope records on plain images (the image IS the layer).  Blockers:
 * Tier-B records, records whose mask input no longer resolves, unknown ops /
 * missing deserializers / newer-version params / unparsable params,
 * canvas-scope records targeting a FLIS layer whose baseline recorded no
 * starting position (there is nothing to anchor the moves to), a flatten or an
 * un-re-runnable merge-down that replaced the item's pixels, and a missing
 * baseline checkpoint.  DOCUMENT-scope property/structure records are ignored
 * (they do not touch the item's pixels; layer.add is embodied in the
 * baseline) — EXCEPT a composite node, which targets the item it produced and
 * is an ordinary member when its inputs were recorded (nde_composite.h).
 * Never returns NULL.
 */
nde_chain *nde_chain_build(gint item_id);

void nde_chain_free(nde_chain *chain);

/* ---- policy predicates (shared by the History panel and edit_execute) ----
 * POLICY only: these answer "could this record ever be edited / deleted?"
 * from the record's own kind, so the UI can grey out buttons without a
 * replay.  They do NOT check liveness or whether the surviving trial chain
 * is actually replayable — the execute path (nde_amend_execute /
 * nde_delete_execute) is the authority on both.                             */
struct nde_record;

/** TRUE for a Tier-A record whose params can be read back, edited and
 *  round-tripped.  For an ordinary op that means a registered deserializer;
 *  the descriptor-less families check their own blobs, which the registry
 *  knows how to ask (nde_op_class_params_valid). */
gboolean nde_record_amendable(const struct nde_record *rec);

/** TRUE unless the record is structural (DOCUMENT scope) — those cannot
 *  resurrect what they destroyed and so cannot be deleted from the log. */
gboolean nde_record_deletable(const struct nde_record *rec);

/**
 * TRUE when @item_id names nothing in the document any more but is still an
 * input to a composite — a RETAINED INPUT (nde_composite.h).  Its steps stay
 * editable, and an edit reaches the image that consumed it.
 *
 * The UI needs this to know what it CANNOT do: there is no such layer to make
 * active and nothing of it on screen, so an editor that previews against the
 * displayed image has nothing to preview.
 */
gboolean nde_item_is_retained_input(gint item_id);

/**
 * Replay the WHOLE chain from the item's baseline into a freshly allocated
 * fits (caller clearfits()+frees()s).  Requires chain->replayable.  NULL on
 * failure with a heap message in @err (caller g_free()s).  Honours
 * cancellation between records.  Job-slot contract in the file header.
 */
struct ffit *nde_chain_replay(const nde_chain *chain, gchar **err);

/**
 * The state of @item_id just BEFORE record @before_record_id's position in
 * the live log (a positional, EXCLUSIVE prefix of the item's chain).  How a
 * joint record (nde_joint.h) resolves its siblings: exclusive because the
 * joint record is itself a member of the sibling's chain — resolving "up to
 * and including" would replay it and recurse — and positional so steps
 * inserted between a sibling's last record and the joint record are
 * honoured.  Caller owns the result.  Job-slot contract as nde_chain_replay.
 */
struct ffit *nde_replay_resolve_before(gint item_id, gint64 before_record_id,
                                       gchar **err);

/**
 * The id of the record the replay driver is currently applying, or 0 outside
 * one.  For hooks whose recompute reads record-scoped side data: the
 * photometric pipeline fetches its embedded star catalogue (nde_cat.h) under
 * this id instead of re-running the network conesearch.
 */
gint64 nde_replay_current_record_id(void);

/** Restore/override the current-record id — for the JOINT recompute
 *  (nde_joint.c), whose sibling resolutions recurse through the replay
 *  driver and clear the id before the joint record's own analysis runs.
 *  Not for general use. */
void nde_replay_set_current_record(gint64 record_id);

/**
 * Replay only the editable tail, starting from the chain's restart point
 * (the last barrier's output checkpoint, or the baseline when there is no
 * barrier).  Requires chain->tail_replayable.  An empty tail returns the
 * restart pixels themselves.  Same ownership/err/job contract as above.
 */
struct ffit *nde_chain_replay_tail(const nde_chain *chain, gchar **err);

/* ---- amend-and-replay commit machinery (phase 3, P3.B) ------------------
 * An amend/delete recomputes the target's pixels from the baseline through
 * the modified chain, then commits atomically: scratch replay first, a
 * single writer-locked swap into the target fits, then the log-side commit
 * (nde_history_amend/delete) and an undo_flush() — there is no meta-undo
 * (sketch §7 decision, UI confirms before calling).  Failure at any
 * pre-commit stage leaves pixels and log untouched.                        */

/**
 * Job-context implementations (the caller owns the processing slot; the
 * _start wrappers and tests call these).  @new_params is the full new kv
 * blob for the record.  FALSE + heap @err message on any failure.
 */
gboolean nde_amend_execute(gint64 record_id, const gchar *new_params, gchar **err);
gboolean nde_delete_execute(gint64 record_id, gchar **err);

/**
 * Spawn the edit as a processing job (refuses while the thread is reserved
 * for a Python script; errors are logged).  For the command layer and the
 * History panel.  @new_params NULL is not a delete — use the dedicated
 * wrapper.  Returns FALSE when the job could not start.
 */
gboolean nde_amend_start(gint64 record_id, const gchar *new_params);
gboolean nde_delete_start(gint64 record_id);

/**
 * Reorder (C3.5): move live Tier-A chain member @record_id so it sits
 * before (@after == FALSE) or after (@after == TRUE) chain member
 * @anchor_id of the SAME item.  Both the old and the new position must
 * lie in the editable tail (crossing an opaque barrier refuses with the
 * locked reason).  Replays from the earliest affected position's cached
 * restart and commits atomically like amend/delete.  A no-op move
 * returns TRUE without replaying.
 */
gboolean nde_reorder_execute(gint64 record_id, gint64 anchor_id, gboolean after, gchar **err);
gboolean nde_reorder_start(gint64 record_id, gint64 anchor_id, gboolean after);

/**
 * flis_replay_check: replay @item_id's whole chain (or its editable tail) on
 * the conductor and log the deviation against the current pixels.  Runs under
 * the SLOT_REPLAY reservation like every other replay — mandatory for chains
 * containing Tier-C members, whose script re-runs need the reservation's
 * capture/undo suppression and command admission.  Returns FALSE when the
 * slot could not be reserved.
 */
gboolean nde_replay_check_start(gint item_id);

/* ---- amend preview (convergence C4) -------------------------------------
 * While record K is being edited in its native dialog, the target image
 * temporarily shows the chain up to K−1 (downstream steps hidden); the
 * dialog previews against that state through the ordinary preview-backup
 * machinery, and Apply routes the widget state back through the normal
 * amend path.  Single instance: one step is editable at a time, and every
 * history edit (amend/delete/reorder) is refused while a preview is
 * installed.  The true pixels are swapped out wholesale (metadata
 * included) and restored bit-exactly on end — Apply restores them FIRST,
 * so a failed amend changes nothing.                                      */

/** Ready callback for nde_amend_preview_start(): fires on the main thread
 *  once the pre-K state is installed (@ok TRUE — open the dialog now) or
 *  synthesis failed (@ok FALSE — state already cleaned up). */
typedef void (*nde_amend_preview_ready_fn)(gboolean ok, gpointer user);

/**
 * Spawn the pre-K synthesis as a processing job and install the result
 * into the record's target image.  Only the record's target may be the
 * ACTIVE image (gfit): the dialog's preview pipeline works on gfit.
 * Returns FALSE when the job could not start (a preview is already
 * active, python owns the thread, ...) — @on_ready is not called then.
 */
gboolean nde_amend_preview_start(gint64 record_id,
                                 nde_amend_preview_ready_fn on_ready,
                                 gpointer user);

/**
 * Leave amend-preview mode.  Restores the true pixels first (bit-exact
 * swap-back); with @apply TRUE then runs the normal amend of the edited
 * record with @new_params (the full new kv blob) as a processing job.
 * Cancelling when no preview is active is a tolerated no-op.
 */
gboolean nde_amend_preview_end(gboolean apply, const gchar *new_params);

/** TRUE from a successful start until end (or a failed synthesis). */
gboolean nde_amend_preview_active(void);

/* Record identity/params for the dialog's pre-fill.  Valid only while
 * active, between the ready callback and end. */
gint64       nde_amend_preview_record_id(void);
const gchar *nde_amend_preview_op_id(void);
gint         nde_amend_preview_op_version(void);
const gchar *nde_amend_preview_params(void);

/**
 * Job-context implementations (the caller owns the processing slot; the
 * wrappers above and the tests call these).  FALSE + heap @err on failure.
 */
gboolean nde_amend_preview_begin_execute(gint64 record_id, gchar **err);
gboolean nde_amend_preview_end_execute(gboolean apply, const gchar *new_params, gchar **err);

/* ---- region-scoped tail replay (roi-nde-plan.md phase 9) ----------------
 * The amend preview shows the chain up to K−1 with the steps AFTER K hidden,
 * so the user tunes a parameter against a state the image will never look
 * like.  Recomputing K+1…head every slider tick is what makes it honest, and
 * it is affordable exactly when it is confined to the ROI rectangle.
 *
 * So when a region preview runs while an amend preview is installed, the
 * worker grows its crop by K's halo PLUS the tail's, runs K's hook on the
 * crop, then replays K+1…head over it before pasting the requested rectangle
 * back.  Inside the rectangle the user sees the whole chain; outside it the
 * pre-K state — which is the contract every region preview has always had
 * ("outside the rectangle, the pre-operation state"), not a new one.
 *
 * A tail is region-replayable iff every record in it is Tier-A with a
 * registered descriptor, declares OP_ROI_CAPABLE, does not change geometry,
 * and is parameter-complete.  In code the last of those reads "has no
 * replay_pre": that hook exists precisely for records that re-derive
 * something from the image they are about to run on, which is the one thing a
 * crop cannot reproduce.  Halos ADD down the tail.
 *
 * Nothing here is recorded, deposited or committed.  A region replay produces
 * preview pixels for one rectangle, and its intermediate states are worth
 * nothing to anyone — a snapstore deposit of a region-sized state would be
 * actively harmful, since a later full replay could restart from it.
 */

/** Opaque plan: a snapshot of the tail, taken once per preview run so the
 *  halo the crop is grown by and the records replayed into it are the same
 *  list.  Free with nde_region_tail_free(). */
typedef struct nde_region_tail nde_region_tail;

/**
 * Plan the tail replay for a preview about to run @editing.  NULL when there
 * is nothing to do — no amend preview installed, a different op than the one
 * being amended, or a tail that cannot be region-replayed.  On success *halo
 * (optional) receives the pixels of extra context the tail needs, which the
 * caller must add to its own op's halo before cropping.
 *
 * An EMPTY tail (the amended record is the last one) plans successfully with
 * a zero halo: there is nothing to replay and nothing hidden, which is the
 * truthful answer rather than a refusal.
 */
nde_region_tail *nde_region_tail_begin(const struct op_descriptor *editing,
                                       int *halo);

/**
 * Replay the planned tail over @region in place.  @rect is the crop @region
 * was taken from, in the amended item's coordinates — needed to crop the
 * records' pinned masks to match.  @max_threads is passed through.
 *
 * FALSE when the tail did not complete: @region then holds a partial result,
 * and the caller should paste it anyway.  That degrades the rectangle to what
 * the REST of the image is already showing (pre-K plus the step being
 * edited), which beats freezing the preview.
 *
 * Two ways not to complete, and telling them apart is load-bearing.  A record
 * that FAILS is a fault: logged once, and the plan refuses to build again for
 * the rest of this amend preview, so a broken tail costs one message rather
 * than one per slider tick.  Being CANCELLED is not: notify_update cancels
 * the running preview job every time a newer tick arrives, so it happens
 * constantly during an ordinary drag — it is silent and leaves the feature
 * armed.  Treating the two alike disarmed the tail for good the first time
 * the user moved the rectangle.
 */
gboolean nde_region_tail_apply(nde_region_tail *plan, struct ffit *region,
                               const rectangle *rect, int max_threads);

void nde_region_tail_free(nde_region_tail *plan);

/**
 * The regime, for the amend banner to state once at open: TRUE when the
 * installed amend preview's tail would be recomputed inside a region preview.
 * @halo (optional) receives the tail's halo, @why (optional) a translated
 * one-line reason when the answer is FALSE *and there is one to give* — it
 * stays NULL when no amend preview is installed at all, which is not a
 * refusal but a different situation entirely.  Caller g_free()s it.
 */
gboolean nde_region_tail_available(int *halo, gchar **why);

/* ---- edit at / insert before K (graph step 2) ---------------------------
 * The same pre-K install as the amend preview, with a different exit verb.
 * Instead of re-opening the anchor's own dialog, the user runs ORDINARY
 * operations against the pre-anchor state; nde_history_append() inserts them
 * before the anchor (nde_history.h), and finishing replays the anchor and
 * everything after it forward over the result.
 *
 * The inserted steps are never themselves replayed — their output is the
 * pixels on screen — so an opaque insert is legal; it simply becomes a
 * barrier for later edits.  What must replay is [anchor..end], which is
 * validated when the point opens, before any pixels change.  Shares the
 * single-instance lock with the amend preview: nde_amend_preview_active() is
 * TRUE for both, and every history edit is refused while either is open.  */

/** Open an insertion point before live record @record_id (which must be a
 *  chain member in the editable tail of its item).  @on_ready fires on the
 *  main thread as for nde_amend_preview_start().  Truncates the dead tail
 *  and flushes undo. */
gboolean nde_edit_at_start(gint64 record_id,
                           nde_amend_preview_ready_fn on_ready, gpointer user);

/**
 * Close the insertion point.  @apply TRUE commits: restore the true pixels,
 * replay the anchor and its successors over the inserted work, swap in.
 * @apply FALSE (or a commit that cannot proceed) discards the inserted
 * records and leaves the image exactly as it was.  Cancelling with no point
 * open is a tolerated no-op.
 */
gboolean nde_edit_at_end(gboolean apply);

/** TRUE while an insertion point is open (a strict subset of
 *  nde_amend_preview_active()); the anchor's record id, or 0. */
gboolean nde_edit_at_active(void);
gint64   nde_edit_at_record_id(void);

/* ---- pinned mask inputs (graph step 4) ----------------------------------
 * A record that ran under a mask carries a "mask" input pin naming the mask
 * item and the record after which its state is wanted (nde_history.h).  For
 * the replay to reproduce the operation rather than refuse it, that state has
 * to be retrievable — so the consumer stores it at capture, under the pin's
 * own coordinate, using the ordinary checkpoint store: a mask is a mono image
 * and needs no store of its own.  Storing at the CONSUMER rather than at every
 * mask op means one state per masked operation instead of one per mask edit,
 * and two operations sharing a mask state share the stored copy.
 *
 * This is a pin, not a live edge: the stored state is what the mask WAS.
 * Re-deriving it by replaying the mask's own chain is a later step, and it
 * changes nothing here — the pin's coordinate is what both resolve.          */

/**
 * Store the mask currently on @fit as the state named by @rec's "mask" pin.
 * No-op when @rec has no mask pin or @fit carries no mask.  Call at capture,
 * after the pin is attached and while the mask is still in place.
 */
void nde_mask_pin_store(const struct nde_record *rec, const struct ffit *fit);

/** TRUE when @rec's mask pin resolves to a stored state, so a replay of the
 *  record can reproduce its mask.  Registry lookup only — no I/O. */
gboolean nde_mask_pin_resolvable(const struct nde_record *rec);

/**
 * Gate for the operations an open insertion point cannot survive: canvas
 * geometry on a LAYERED document, flatten and merge-down.  Their pixel
 * effect is not part of the inserted lineage, and closing the point restores
 * the true image over it — so they must not run while a point is open.
 * Returns TRUE (having logged an error naming @what, e.g. "Flattening") when
 * the caller must refuse; FALSE when there is nothing to refuse.
 */
gboolean nde_edit_at_refuses_op(const char *what);

/** Job-context implementations (the caller owns the processing slot). */
gboolean nde_edit_at_begin_execute(gint64 record_id, gchar **err);

/**
 * The item a capture should be stamped with while an insertion is open on an
 * item that has no layer of its own (one a merge or flatten consumed).  0 when
 * no such insertion is open, in which case the capture uses the active layer
 * as usual.  Without this the new record carries the ACTIVE layer's id, does
 * not qualify for the armed insertion, and is appended to the end of the
 * merged result instead of landing where it was aimed.
 */
gint nde_edit_at_borrowed_item(void);

/* ---- undoing a merge or a flatten ---------------------------------------
 * Deleting a composite step is not a chain edit: the step did not change an
 * image, it CONSUMED several and produced one.  Undoing it therefore restores
 * the document — the consumed layers come back with the identities, names and
 * compositing properties the record kept — and takes everything built on the
 * result with it, because those steps describe an image that stops existing.
 * It is a bulk undo to just before the merge, and the caller must say so
 * before doing it.
 */

/** A human-readable list of the steps that would be discarded, one per line,
 *  for the confirmation the caller must show.  @n_layers receives how many
 *  layers would come back and @n_discarded how many steps would go.  NULL +
 *  @err when the undo is not possible, with the reason spelled out. */
gchar *nde_composite_undo_describe(gint64 record_id, guint *n_layers,
                                   guint *n_discarded, gchar **err);

/** Perform it.  Job-slot contract as nde_chain_replay (the caller owns the
 *  slot).  FALSE + @err leaves the document untouched. */
gboolean nde_composite_undo_execute(gint64 record_id, gchar **err);

/** Async wrapper: runs the above on the replay conductor. */
gboolean nde_composite_undo_start(gint64 record_id);
gboolean nde_edit_at_end_execute(gboolean apply, gchar **err);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_REPLAY_H_ */
