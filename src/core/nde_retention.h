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
#ifndef _NDE_RETENTION_H_
#define _NDE_RETENTION_H_

/**
 * \file nde_retention.h
 * \brief One storage budget for everything nondestructive editing keeps.
 *
 * ONE BUDGET, TWO PRIORITIES (design note §7).  `com.pref.nde_cache_mb` is
 * not the cache's budget — it is the budget for ALL edit-history storage,
 * measured against nde_snapstore_total_bytes().  What differs between the
 * two kinds of stored data is not their allowance but their eviction order:
 *
 *   1. SPEED CACHE (the snapstore LRU pool) goes first, silently.  Losing a
 *      pool entry costs time: a later replay restarts further back.
 *   2. EDITABILITY PINS (baselines, output checkpoints, pinned operation
 *      masks) go only when the cache is exhausted, and are ANNOUNCED.
 *      Losing one costs a capability: steps that depended on it become
 *      fixed and can no longer be amended.
 *
 * The honesty requirement is the part not to compromise on: a silent
 * capability downgrade is worse than a refusal.  The user is told at the
 * moment it happens, not when a later amend fails.
 *
 * There is no "off".  A single always-on budget is one code path; a machine
 * short of disk sets the number small and gets the old behaviour without a
 * separate mode to test.
 *
 * Threading: callable from any thread.  Takes the snapstore and checkpoint
 * leaf locks one at a time, never together and never across I/O, so it
 * observes the same lock order as everything else in the NDE core.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Bring total NDE storage back within the budget, cache first and pins only
 * if that is not enough.  A cheap no-op when already inside it, so capture
 * sites can call it unconditionally after storing.
 *
 * @keep_record_id is the checkpoint just stored (0 for none), which is never
 * evicted: a budget too small to hold even one checkpoint must degrade into
 * "keeps the newest", not into storing and immediately dropping it.
 */
void nde_retention_enforce(gint64 keep_record_id);

/**
 * Forget that the limit has been announced (document closed or replaced), so
 * the next document explains itself once rather than staying silent because a
 * previous one already said it.
 */
void nde_retention_notice_reset(void);

/** The configured budget in bytes (0 disables all deposits, as it does now). */
gint64 nde_retention_budget_bytes(void);

/* ---- diagnostics / tests ----------------------------------------------- */

typedef struct {
	guint  cache_evictions;    /* silent pool reclaims that freed something */
	guint  pins_dropped;       /* announced editability losses */
	guint  baselines_dropped;  /* of those, whole-chain losses */
	gint64 bytes_reclaimed;
} nde_retention_stats_t;

void nde_retention_stats(nde_retention_stats_t *out);
void nde_retention_stats_reset(void);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_RETENTION_H_ */
