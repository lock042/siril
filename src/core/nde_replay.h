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
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;

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
	GArray    *member_flags;    /* guint8 per member (NDE_CHAIN_MEMBER_*) */
	guint      tail_start;      /* first tail member index (== len ⇒ empty tail) */
	gboolean   tail_replayable; /* the tail applies cleanly from its restart point */
	gint64     restart_ckpt_id; /* record id of the restart checkpoint; 0 = baseline */
} nde_chain;

/**
 * Snapshot the live history and extract + validate the replay chain for
 * @item_id.  Membership: LAYER-scope records targeting the item; CANVAS-
 * scope records on plain images (the image IS the layer).  Blockers:
 * Tier-B records, active-mask records, unknown ops / missing
 * deserializers / newer-version params / unparsable params, canvas-scope
 * records targeting the item on a FLIS document (their layer-offset side
 * effects cannot be verified on a scratch fits yet), merge-down/flatten
 * that replaced the item's pixels, and a missing baseline checkpoint.
 * DOCUMENT-scope property/structure records are ignored (they do not
 * touch the item's pixels; layer.add is embodied in the baseline).
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

/** TRUE for a Tier-A record whose op has a registered deserializer (so its
 *  params can be parsed, edited and round-tripped). */
gboolean nde_record_amendable(const struct nde_record *rec);

/** TRUE unless the record is structural (DOCUMENT scope) or records
 *  compositing state (layer.set_opacity/blend/visible) — those cannot be
 *  meaningfully deleted from the log. */
gboolean nde_record_deletable(const struct nde_record *rec);

/**
 * Replay the WHOLE chain from the item's baseline into a freshly allocated
 * fits (caller clearfits()+frees()s).  Requires chain->replayable.  NULL on
 * failure with a heap message in @err (caller g_free()s).  Honours
 * cancellation between records.  Job-slot contract in the file header.
 */
struct ffit *nde_chain_replay(const nde_chain *chain, gchar **err);

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

#ifdef __cplusplus
}
#endif

#endif /* _NDE_REPLAY_H_ */
