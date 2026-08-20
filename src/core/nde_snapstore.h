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
#ifndef _NDE_SNAPSTORE_H_
#define _NDE_SNAPSTORE_H_

/**
 * \file nde_snapstore.h
 * \brief Shared refcounted pixel-snapshot store (convergence phase C1,
 * nde-convergence-plan.md; sketch §7 decision 13).
 *
 * One snapshot = one delete-on-close swap file holding a fits' pixels plus
 * the geometry/WCS needed to reconstruct it.  Snapshots are refcounted;
 * owners today are the NDE checkpoint tables (baseline pinned per item,
 * barrier outputs per record), the LRU cache pool below, and — from C2 —
 * undo entries' main pixel snapshots.
 *
 * TAG REGISTRY: a snapshot may carry one position tag
 * {item_id, record_id, post}: post=FALSE means "state just BEFORE
 * record_id" (undo entries), post=TRUE "state just AFTER" (barrier
 * checkpoints, redo counterparts, replay deposits; the baseline is
 * POST(0)).  The registry holds NO references — it indexes live snapshots
 * and an entry vanishes when the last owner unrefs.  Lookups return fresh
 * fits copies, never handles.
 *
 * LRU POOL: store-owned deposits made by the replay driver so successive
 * history edits restart near the edit point instead of the baseline.
 * Bounded by com.pref.nde_cache_mb of DISK (swap files, not RAM); evicted
 * by LRU, by the invalidation rules (any edit at position M drops
 * POST(K>=M) and PRE(K>M) pool entries for that item), and wholesale on
 * document change.  The pool is pure cache: correctness never depends on
 * it.
 *
 * Threading: an internal strict-LEAF mutex guards the registry/pool
 * bookkeeping; refcounts are atomic; all swap I/O and fits copies run
 * OUTSIDE the lock.  nde_snap_unref() may take the leaf mutex (final-unref
 * deregistration), so never call it while holding another NDE leaf lock —
 * use the steal-then-unref pattern (see nde_checkpoint.c).
 */

#include <glib.h>
#include "core/nde_state.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct nde_snap nde_snap;

/* ---- snapshot lifecycle ------------------------------------------------ */

/** Deep-copy @st's pixels (+ geometry/WCS) and its canvas position into a
 *  fresh snapshot, ref 1.  NULL on error.  Swap I/O only — safe while other
 *  locks are held.  @st is only read; it is not taken over. */
nde_snap *nde_snap_create(const nde_state *st);

nde_snap *nde_snap_ref(nde_snap *s);

/** Drop one reference; at zero the snapshot deregisters and its swap file
 *  vanishes.  May take the store's leaf mutex — do not call while holding
 *  another NDE leaf lock. */
void nde_snap_unref(nde_snap *s);

/** Reconstruct the stored value (caller nde_state_free()s).  The caller must
 *  hold a reference across the call. */
nde_state *nde_snap_read(const nde_snap *s);

/**
 * Read the snapshot's pixels INTO @dst, undo-restore style: buffers are
 * realloc'd, dims/type and pdata/fpdata set, stats invalidated.  Header
 * metadata (WCS, ICC, keywords) is deliberately NOT touched — undo's
 * historic_struct keeps carrying and restoring those itself.  0 on
 * success.
 */
int nde_snap_read_into(const nde_snap *s, fits *dst);

/** On-disk payload size in bytes (for budgets/diagnostics). */
gint64 nde_snap_size(const nde_snap *s);

/* ---- tag registry ------------------------------------------------------ */

/**
 * Tag @s as the state at a timeline position and index it.  A later tag
 * with the same key replaces the index entry (latest wins); the index
 * entry disappears when @s's last reference drops.
 */
void nde_snap_set_tag(nde_snap *s, gint item_id, gint64 record_id, gboolean post);

/** TRUE when a live snapshot is indexed at this position. */
gboolean nde_snapstore_has(gint item_id, gint64 record_id, gboolean post);

/** Where the stored value sits on the canvas, without reading its pixels back.
 *  FALSE when it has none.  @pos_x / @pos_y may be NULL. */
gboolean nde_snap_position(const nde_snap *s, gint *pos_x, gint *pos_y);

/** Read @s's tag (immutable once set).  FALSE when untagged. */
gboolean nde_snap_tag_get(const nde_snap *s, gint *item_id, gint64 *record_id, gboolean *post);

/** Fresh copy of the indexed state (caller nde_state_free()s), or NULL.
 *  Counts as a lookup (and a hit) in the stats; touches pool LRU. */
nde_state *nde_snapstore_lookup(gint item_id, gint64 record_id, gboolean post);

/* ---- LRU cache pool ---------------------------------------------------- */

/** Deposit @state as POST(@record_id) for @item_id into the pool (replaces
 *  any same-tag pool entry; evicts LRU over budget; no-op when the budget
 *  pref is 0).  Failures are silent — the pool is cache, not correctness. */
void nde_snapstore_deposit(const nde_state *state, gint item_id, gint64 record_id);

/** Invalidation rule for an edit at position @record_id of @item_id:
 *  evict pool entries POST(K >= record_id) and PRE(K > record_id). */
void nde_snapstore_invalidate_from(gint item_id, gint64 record_id);

/** Evict pool entries tagged with @record_id (the record was deleted or
 *  truncated — any item, PRE and POST alike). */
void nde_snapstore_evict_record(gint64 record_id);

/** Evict the whole pool (document change). */
void nde_snapstore_pool_purge(void);

/* ---- budget accounting (nde_retention.c drives these) ------------------- */

/**
 * Payload bytes of EVERY live snapshot, whoever owns it — pool, checkpoint
 * tables, undo.  This, not the pool alone, is what the single retention
 * budget is measured against (design note §7).
 */
gint64 nde_snapstore_total_bytes(void);

/** Payload bytes held by the LRU pool alone (the freely-evictable part). */
gint64 nde_snapstore_pool_bytes(void);

/**
 * Evict pool entries LRU-first until @want bytes have been freed or the pool
 * is empty; returns the bytes actually freed.  This is the SILENT half of
 * retention — losing a pool entry costs time, never a capability, so it
 * happens before anything the user must be told about.
 */
gint64 nde_snapstore_reclaim_pool(gint64 want);

/* ---- stats (tests/diagnostics) ----------------------------------------- */

typedef struct {
	guint lookups;
	guint hits;
	guint deposits;
	guint evictions;
} nde_snapstore_stats_t;

void nde_snapstore_stats(nde_snapstore_stats_t *out);
void nde_snapstore_stats_reset(void);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_SNAPSTORE_H_ */
