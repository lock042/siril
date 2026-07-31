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
 * ONE BUDGET, TWO KINDS OF DATA (design note §7).  `com.pref.nde_cache_mb` is
 * not the cache's budget — it is the budget for ALL edit-history storage,
 * measured against nde_snapstore_total_bytes().  What differs between the two
 * kinds of stored data is whether the budget may take them:
 *
 *   1. SPEED CACHE (the snapstore LRU pool) is evictable, oldest first, and
 *      goes silently.  Losing a pool entry costs only time: a later replay
 *      restarts further back and recomputes more.
 *   2. RECONSTRUCTION DATA (baselines, barrier output checkpoints, pinned
 *      operation masks) is NEVER evicted, even when that means exceeding the
 *      budget.  Losing one of these does not cost time, it costs the ability
 *      to rebuild the image at all: a baseline is where a chain starts, and a
 *      barrier's checkpoint is the only thing the steps after an
 *      unreplayable step can be rebuilt from.
 *
 * The rule that decides this: an edit history that cannot reconstruct its
 * image is not an edit history.  A budget is a request about disk; being
 * unable to re-run the user's own processing is a broken promise, and the
 * second outranks the first.  So the budget bends.  It is announced once so
 * the overshoot is not a surprise, and raising the setting buys back cache —
 * it does not change what is kept.
 *
 * There is no "off".  A single always-on budget is one code path; a machine
 * short of disk sets the number small, keeps only what reconstruction needs,
 * and recomputes the rest.
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
 * Reclaim from the speed cache until total NDE storage is inside the budget,
 * or until the cache is empty and only reconstruction data remains — in which
 * case the budget is knowingly exceeded and said so once.  A cheap no-op when
 * already inside it, so capture sites can call it unconditionally.
 */
void nde_retention_enforce(void);

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
	gint64 bytes_reclaimed;
	/* What the budget could not take because dropping it would cost the
	 * ability to rebuild an image: 0 while inside the budget, otherwise the
	 * most recent overshoot.  There is no pins_dropped counter any more —
	 * the answer is structurally zero. */
	gint64 bytes_over_budget;
} nde_retention_stats_t;

void nde_retention_stats(nde_retention_stats_t *out);
void nde_retention_stats_reset(void);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_RETENTION_H_ */
