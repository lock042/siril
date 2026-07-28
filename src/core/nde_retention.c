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
#include "core/nde_snapstore.h"
#include "core/nde_checkpoint.h"
#include "core/nde_retention.h"

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

static void announce(guint pins, guint baselines) {
	if (!pins)
		return;
	if (!warned_limit) {
		warned_limit = TRUE;
		siril_log_message(_("Edit history storage limit reached (%d MB). The oldest "
		                    "steps can no longer be re-run, and are marked as fixed "
		                    "in the history.\n"), com.pref.nde_cache_mb);
		/* Name the setting rather than a Preferences page: this one has no GUI
		 * control, so "see Preferences" would send the user looking for
		 * something that is not there. */
		siril_log_message(_("Raise core.nde_cache_mb to keep more of the history "
		                    "editable.\n"));
	}
	/* A lost baseline is a strictly worse event — a whole image's history stops
	 * being amendable, not just its oldest steps — and is rare enough that
	 * saying so every time cannot become noise. */
	if (baselines)
		siril_log_message(_("%u image(s) lost their original state, so none of their "
		                    "steps can be amended.\n"), baselines);
}

void nde_retention_notice_reset(void) {
	g_mutex_lock(&ret_mutex);
	warned_limit = FALSE;
	g_mutex_unlock(&ret_mutex);
}

void nde_retention_enforce(gint64 keep_record_id) {
	gint64 budget = nde_retention_budget_bytes();
	if (budget <= 0)
		return;   /* deposits are already suppressed; pins are all that is left */

	/* Serialised so two workers finishing together cannot both decide the
	 * same pins must go and announce it twice. */
	g_mutex_lock(&ret_mutex);

	gint64 total = nde_snapstore_total_bytes();
	if (total <= budget) {
		g_mutex_unlock(&ret_mutex);
		return;
	}

	/* 1. Cache first, silently. */
	gint64 freed = nde_snapstore_reclaim_pool(total - budget);
	if (freed) {
		ret_stats.cache_evictions++;
		ret_stats.bytes_reclaimed += freed;
	}

	/* 2. Only then pins, and only as many as the overshoot actually needs.
	 * Re-reading the total each round rather than subtracting keeps this
	 * honest when another thread frees something underneath us. */
	guint pins = 0, baselines = 0;
	while (nde_snapstore_total_bytes() > budget) {
		nde_checkpoint_eviction ev = { 0 };
		if (!nde_checkpoint_evict_oldest(&ev, keep_record_id))
			break;   /* nothing left that may go — see the header's floor */
		pins++;
		if (ev.is_baseline)
			baselines++;
		ret_stats.pins_dropped++;
		ret_stats.baselines_dropped += ev.is_baseline ? 1 : 0;
		ret_stats.bytes_reclaimed += ev.bytes;
	}
	g_mutex_unlock(&ret_mutex);

	announce(pins, baselines);
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
