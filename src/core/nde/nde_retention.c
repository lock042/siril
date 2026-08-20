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

/* See nde_retention.h for the policy this implements. */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/nde/nde_snapstore.h"
#include "core/nde/nde_checkpoint.h"
#include "core/nde/nde_retention.h"

static GMutex ret_mutex;                  /* serialises enforcement only */
static nde_retention_stats_t ret_stats;

gint64 nde_retention_budget_bytes(void) {
	return (gint64)com.pref.nde_cache_mb * 1024 * 1024;
}

/* Announced ONCE per document, not once per pass.  Once the budget is full,
 * every subsequent operation evicts something, so a per-pass message repeats
 * on every command and reads as a fault rather than as the policy working.
 * The ongoing signal is the History marking those steps fixed; this line
 * exists to explain the first time it happens.  Reset by
 * nde_retention_notice_reset() when the document changes. */
static gboolean warned_limit;

static void announce(gint64 over) {
	if (over <= 0 || warned_limit)
		return;
	warned_limit = TRUE;
	siril_log_message(_("Edit history storage is over the %d MB budget by %.0f MB. "
	                    "The states needed to rebuild your images are kept "
	                    "regardless — only cached intermediate results are "
	                    "discarded, so nothing becomes un-editable.\n"),
	                  com.pref.nde_cache_mb, (double)over / (1024.0 * 1024.0));
	/* Name the setting rather than a Preferences page: this one has no GUI
	 * control, so "see Preferences" would send the user looking for something
	 * that is not there.  And be accurate about what raising it buys — not
	 * editability, which is never traded away, but the cache that makes an
	 * edit fast instead of slow. */
	siril_log_message(_("Raise core.nde_cache_mb to keep more intermediate "
	                    "results and make history edits faster.\n"));
}

void nde_retention_notice_reset(void) {
	g_mutex_lock(&ret_mutex);
	warned_limit = FALSE;
	g_mutex_unlock(&ret_mutex);
}

void nde_retention_enforce(void) {
	gint64 budget = nde_retention_budget_bytes();
	if (budget <= 0)
		return;   /* deposits are already suppressed; reconstruction data is all
		           * that is left, and that is kept whatever the budget says */

	/* Serialised so two workers finishing together cannot both reclaim for the
	 * same overshoot and announce it twice. */
	g_mutex_lock(&ret_mutex);

	gint64 total = nde_snapstore_total_bytes();
	if (total <= budget) {
		ret_stats.bytes_over_budget = 0;
		g_mutex_unlock(&ret_mutex);
		return;
	}

	/* The cache is the only thing the budget may take.  Baselines, barrier
	 * checkpoints and pinned masks are what the image is rebuilt FROM; dropping
	 * one to save disk trades a capability for a number, which is the wrong way
	 * round (see the header). */
	gint64 freed = nde_snapstore_reclaim_pool(total - budget);
	if (freed) {
		ret_stats.cache_evictions++;
		ret_stats.bytes_reclaimed += freed;
	}

	/* Re-read rather than subtract: another thread may have freed something
	 * underneath us, and the overshoot we report should be the real one. */
	gint64 over = nde_snapstore_total_bytes() - budget;
	ret_stats.bytes_over_budget = over > 0 ? over : 0;
	g_mutex_unlock(&ret_mutex);

	announce(over);
}

void nde_retention_stats(nde_retention_stats_t *out) {
	g_mutex_lock(&ret_mutex);
	*out = ret_stats;
	g_mutex_unlock(&ret_mutex);
}

void nde_retention_stats_reset(void) {
	g_mutex_lock(&ret_mutex);
	memset(&ret_stats, 0, sizeof(ret_stats));
	g_mutex_unlock(&ret_mutex);
}
