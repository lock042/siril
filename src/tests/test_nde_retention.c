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
 */

/*
 * test_nde_retention — one storage budget over ALL nondestructive-editing
 * storage, with two eviction priorities (design note §7):
 *   - the budget is measured against every live snapshot, not the pool alone
 *   - the speed cache yields to editability pins, never the other way round
 *   - pins go oldest-first, output checkpoints before baselines
 *   - the active item's baseline is never evicted
 *   - what was lost is counted, so the announcement reports fact
 */

#include <criterion/criterion.h>
#include <glib.h>

#include "core/siril.h"
#include "io/image_format_fits.h"
#include "core/nde_history.h"
#include "core/nde_snapstore.h"
#include "core/nde_checkpoint.h"
#include "core/nde_retention.h"

cominfo com;
fits *gfit;

static void setup(void) {
	memset(&com, 0, sizeof(com));
	com.max_thread = 1;
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	com.pref.nde_cache_mb = 1;
	nde_retention_stats_reset();
	nde_snapstore_stats_reset();
}

static void teardown(void) {
	nde_checkpoint_purge();
	nde_snapstore_pool_purge();
	g_free(com.pref.swap_dir);
	com.pref.swap_dir = NULL;
}

TestSuite(nde_retention, .init = setup, .fini = teardown);

/* A float mono image of a known payload size: rx*ry*4 bytes. */
static fits *make_image(int rx, int ry, float base) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, rx, ry, 1, DATA_FLOAT), 0);
	size_t n = (size_t)rx * ry;
	for (size_t i = 0; i < n; i++)
		f->fdata[i] = base + (float)i * 1e-4f;
	return f;
}

/* 256 KB each: four of them exceed a 1 MB budget. */
#define QUARTER_MB_SIDE 256
static fits *quarter_mb(float base) { return make_image(QUARTER_MB_SIDE, QUARTER_MB_SIDE, base); }

/* ---- the budget covers everything, not just the pool -------------------- */

Test(nde_retention, the_budget_counts_pins_as_well_as_cache) {
	fits *pin = quarter_mb(0.1f);
	nde_checkpoint_output_store(pin, 1, 7);

	cr_assert_geq(nde_snapstore_total_bytes(), 256 * 256 * 4,
	              "a checkpoint counts towards the total");
	cr_assert_eq(nde_snapstore_pool_bytes(), 0, "and it is not a pool entry");
	cr_assert_eq(nde_checkpoint_bytes(), nde_snapstore_total_bytes(),
	             "with nothing else stored the two agree");
	clearfits(pin); free(pin);
}

/* The cache must shrink to make room for pins, which is the whole point of a
 * single budget: pool_bytes alone would have let the two grow independently. */
Test(nde_retention, the_cache_yields_to_pins) {
	fits *a = quarter_mb(0.1f), *b = quarter_mb(0.2f);
	nde_snapstore_deposit(a, 1, 101);
	nde_snapstore_deposit(b, 1, 102);
	cr_assert_eq(nde_snapstore_pool_bytes(), 2 * 256 * 256 * 4);

	/* Two pins push the total to 1 MB; the next one must cost the cache. */
	fits *p1 = quarter_mb(0.3f), *p2 = quarter_mb(0.4f), *p3 = quarter_mb(0.5f);
	nde_checkpoint_output_store(p1, 1, 7);
	nde_checkpoint_output_store(p2, 2, 7);
	nde_checkpoint_output_store(p3, 3, 7);

	cr_assert_leq(nde_snapstore_total_bytes(), nde_retention_budget_bytes(),
	              "enforcement must bring the total inside the budget");
	cr_assert_lt(nde_snapstore_pool_bytes(), 2 * 256 * 256 * 4,
	             "and the cache is what gave way");

	nde_retention_stats_t st;
	nde_retention_stats(&st);
	cr_assert_eq(st.pins_dropped, 0, "no pin should have been touched");

	clearfits(a); free(a); clearfits(b); free(b);
	clearfits(p1); free(p1); clearfits(p2); free(p2); clearfits(p3); free(p3);
}

/* ---- pins, oldest first, and only once the cache is gone --------------- */

Test(nde_retention, pins_go_oldest_first_once_the_cache_is_exhausted) {
	/* Five quarter-MB pins against a 1 MB budget: the oldest must go. */
	fits *p[5];
	for (int i = 0; i < 5; i++) {
		p[i] = quarter_mb(0.1f * (float)(i + 1));
		nde_checkpoint_output_store(p[i], 10 + i, 7);
	}

	cr_assert(!nde_checkpoint_output_exists(10), "record 10 was the oldest");
	cr_assert(nde_checkpoint_output_exists(14), "the newest must survive");
	cr_assert_leq(nde_snapstore_total_bytes(), nde_retention_budget_bytes());

	nde_retention_stats_t st;
	nde_retention_stats(&st);
	cr_assert_gt(st.pins_dropped, 0, "the loss must be counted, not silent");
	cr_assert_eq(st.baselines_dropped, 0, "no baseline was involved");

	for (int i = 0; i < 5; i++) { clearfits(p[i]); free(p[i]); }
}

/* procmasksnag: a used mask persists ONLY as its pinned states — clearing it
 * releases the live slot, so once the stored copy is gone every step that ran
 * under it is a permanent barrier.  And its record id is the OLDEST around (a
 * mask precedes everything it masked), which is exactly what the age rule
 * used to pick first. */
Test(nde_retention, a_pinned_mask_state_is_never_evicted) {
	com.uniq = g_new0(single, 1);
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("image.asinh");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = 7;
	rec->summary = g_strdup("masked stretch");
	rec->mask_active = TRUE;
	nde_record_add_input(rec, "mask", 9, 10);
	cr_assert(nde_history_append(rec) > 0);

	/* Five quarter-MB states against a 1 MB budget: one must go, and the
	 * oldest (record 10) is the pinned mask state. */
	fits *p[5];
	for (int i = 0; i < 5; i++) {
		p[i] = quarter_mb(0.1f * (float)(i + 1));
		nde_checkpoint_output_store(p[i], 10 + i, 9);
	}
	cr_assert(nde_checkpoint_output_exists(10),
	          "the pinned mask state must be passed over");
	cr_assert(!nde_checkpoint_output_exists(11),
	          "the oldest UNPINNED state goes instead");
	cr_assert_leq(nde_snapstore_total_bytes(), nde_retention_budget_bytes());

	for (int i = 0; i < 5; i++) { clearfits(p[i]); free(p[i]); }
	nde_history_attach(NULL);
	g_free(com.uniq);
	com.uniq = NULL;
}

/* A baseline costs an item its WHOLE chain, so it goes only after every
 * output checkpoint has already gone. */
Test(nde_retention, baselines_are_evicted_after_output_checkpoints) {
	fits *base = quarter_mb(0.1f);
	nde_checkpoint_baseline_ensure(base, 4);       /* item 4, not the active one */
	cr_assert(nde_checkpoint_baseline_exists(4));

	fits *o1 = quarter_mb(0.2f), *o2 = quarter_mb(0.3f), *o3 = quarter_mb(0.4f);
	nde_checkpoint_output_store(o1, 20, 4);
	nde_checkpoint_output_store(o2, 21, 4);
	nde_checkpoint_output_store(o3, 22, 4);   /* baseline + 3 = exactly 1 MB */
	cr_assert(nde_checkpoint_output_exists(20), "at the budget nothing goes yet");

	/* One past the budget, and it must be the oldest OUTPUT checkpoint that
	 * goes, not the baseline. */
	fits *o4 = quarter_mb(0.5f);
	nde_checkpoint_output_store(o4, 23, 4);

	cr_assert(nde_checkpoint_baseline_exists(4), "the baseline must outlast them");
	cr_assert(!nde_checkpoint_output_exists(20), "the oldest checkpoint went first");
	cr_assert(nde_checkpoint_output_exists(23), "the newest survives");

	clearfits(base); free(base);
	clearfits(o1); free(o1); clearfits(o2); free(o2);
	clearfits(o3); free(o3); clearfits(o4); free(o4);
}

/* Evicting the baseline of the item being worked on would make the user's
 * very next operation unreplayable, before they could act on the warning. */
Test(nde_retention, the_active_items_baseline_is_never_evicted) {
	/* nde_checkpoint_active_item_id() reports -1 for a plain image, which is
	 * the item every capture in a non-FLIS session targets. */
	cr_assert_eq(nde_checkpoint_active_item_id(), -1);

	fits *active = quarter_mb(0.1f);
	nde_checkpoint_baseline_ensure(active, -1);

	/* Baselines ALONE past the budget, so there is nothing cheaper to drop:
	 * five quarter-MB baselines against 1 MB. */
	fits *other[4];
	for (int i = 0; i < 4; i++) {
		other[i] = quarter_mb(0.2f + 0.1f * (float)i);
		nde_checkpoint_baseline_ensure(other[i], 9 + i);
	}

	cr_assert(nde_checkpoint_baseline_exists(-1),
	          "the active item's baseline must survive any pressure");
	cr_assert(!nde_checkpoint_baseline_exists(9),
	          "the oldest inactive baseline is what gives way");
	cr_assert(nde_checkpoint_baseline_exists(12), "the newest survives");

	nde_retention_stats_t st;
	nde_retention_stats(&st);
	cr_assert_eq(st.baselines_dropped, 1, "and the whole-chain loss is reported");

	clearfits(active); free(active);
	for (int i = 0; i < 4; i++) { clearfits(other[i]); free(other[i]); }
}

/* A budget too small to hold even one checkpoint must degrade into "keeps the
 * newest", not into storing each one and immediately dropping it — which is
 * what an unprotected oldest-first rule does when the newest IS also the
 * oldest.  Found by running it, not by reasoning about it. */
Test(nde_retention, the_checkpoint_just_stored_is_never_the_one_evicted) {
	fits *big = make_image(1024, 512, 0.1f);   /* 2 MB, twice the budget */

	nde_checkpoint_output_store(big, 80, 7);
	cr_assert(nde_checkpoint_output_exists(80),
	          "storing it and dropping it again would be pure waste");

	/* The next one displaces it, rather than being discarded on arrival. */
	nde_checkpoint_output_store(big, 81, 7);
	cr_assert(!nde_checkpoint_output_exists(80), "the older one gives way");
	cr_assert(nde_checkpoint_output_exists(81), "the newest is what is kept");

	clearfits(big); free(big);
}

/* ---- the escape hatch -------------------------------------------------- */

/* A zero budget is the "constrained machine" setting, not an off switch for
 * enforcement: deposits are already suppressed, and nothing may be silently
 * dropped on a pref value that means "store as little as possible". */
Test(nde_retention, a_zero_budget_suppresses_the_cache_and_keeps_pins) {
	com.pref.nde_cache_mb = 0;

	fits *a = quarter_mb(0.1f);
	nde_snapstore_deposit(a, 1, 101);
	cr_assert_eq(nde_snapstore_pool_bytes(), 0, "no cache at a zero budget");

	fits *pin = quarter_mb(0.2f);
	nde_checkpoint_output_store(pin, 50, 7);
	cr_assert(nde_checkpoint_output_exists(50),
	          "a pin is correctness, not cache: it is stored regardless");

	nde_retention_stats_t st;
	nde_retention_stats(&st);
	cr_assert_eq(st.pins_dropped, 0);

	clearfits(a); free(a);
	clearfits(pin); free(pin);
}

/* Enforcement must be idempotent: calling it when already inside the budget
 * costs nothing and drops nothing. */
Test(nde_retention, enforcing_inside_the_budget_is_a_no_op) {
	fits *pin = quarter_mb(0.1f);
	nde_checkpoint_output_store(pin, 60, 7);
	nde_retention_stats_reset();

	for (int i = 0; i < 5; i++)
		nde_retention_enforce(0);

	cr_assert(nde_checkpoint_output_exists(60));
	nde_retention_stats_t st;
	nde_retention_stats(&st);
	cr_assert_eq(st.pins_dropped, 0);
	cr_assert_eq(st.bytes_reclaimed, 0);

	clearfits(pin); free(pin);
}

/* The counter must fall back to zero when everything is released, or a long
 * session would drift over budget on phantom bytes. */
Test(nde_retention, released_storage_stops_counting) {
	fits *pin = quarter_mb(0.1f);
	nde_checkpoint_output_store(pin, 70, 7);
	cr_assert_gt(nde_snapstore_total_bytes(), 0);

	nde_checkpoint_output_drop(70);
	cr_assert_eq(nde_snapstore_total_bytes(), 0,
	             "dropping the last pin must return the total to zero");

	clearfits(pin); free(pin);
}
