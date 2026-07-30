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
 * test_nde_snapstore — the shared refcounted snapshot store (convergence
 * phase C1, nde-convergence-plan.md): lifecycle, tag registry ownership
 * semantics, read/read_into equivalence, LRU pool budget + eviction +
 * invalidation, stats.
 */

#include <criterion/criterion.h>
#include <string.h>
#include "flis_test_helpers.h"
#include "core/nde_snapstore.h"

cominfo com;
fits *gfit;

static void setup(void) {
	flis_test_init_com();
	com.pref.nde_cache_mb = 2048;
	nde_snapstore_stats_reset();
	nde_snapstore_pool_purge();
}

static void teardown(void) {
	nde_snapstore_pool_purge();
	flis_test_cleanup_com();
}

TestSuite(nde_snapstore, .init = setup, .fini = teardown);

static fits *make_fixture(guint w, guint h, float seed) {
	fits *f = flis_test_make_mono_fits(w, h, 0.f);
	for (guint i = 0; i < w * h; i++)
		f->fdata[i] = seed + i / (float)(w * h);
	return f;
}

Test(nde_snapstore, lifecycle_and_read_roundtrip) {
	fits *f = make_fixture(8, 6, 0.1f);
	nde_snap *s = nde_snap_create(f);
	cr_assert_not_null(s);
	cr_assert_eq(nde_snap_size(s), 8 * 6 * (gint64)sizeof(float));

	fits *back = nde_snap_read(s);
	cr_assert_not_null(back);
	cr_assert_eq(back->rx, 8);
	cr_assert_eq(back->ry, 6);
	cr_assert(memcmp(back->fdata, f->fdata, 48 * sizeof(float)) == 0,
	          "read must be bit-exact");
	clearfits(back); free(back);

	/* extra ref keeps it readable after one unref */
	nde_snap_ref(s);
	nde_snap_unref(s);
	back = nde_snap_read(s);
	cr_assert_not_null(back);
	clearfits(back); free(back);
	nde_snap_unref(s);

	clearfits(f); free(f);
}

Test(nde_snapstore, read_into_matches_read) {
	fits *f = make_fixture(10, 4, 0.3f);
	nde_snap *s = nde_snap_create(f);
	cr_assert_not_null(s);

	fits *via_read = nde_snap_read(s);
	fits *via_into = flis_test_make_mono_fits(2, 2, 0.f);  /* wrong dims: must realloc */
	cr_assert_eq(nde_snap_read_into(s, via_into), 0);

	cr_assert_eq(via_into->rx, 10);
	cr_assert_eq(via_into->ry, 4);
	cr_assert_eq(via_into->type, DATA_FLOAT);
	cr_assert(memcmp(via_into->fdata, via_read->fdata, 40 * sizeof(float)) == 0);
	cr_assert_eq(via_into->fpdata[0], via_into->fdata, "fpdata must be rewired");

	clearfits(via_read); free(via_read);
	clearfits(via_into); free(via_into);
	nde_snap_unref(s);
	clearfits(f); free(f);
}

/* A snapshot's BITPIX describes the buffer it holds, not the depth its content
 * really has — an 8-bit mask travels widened into WORDs, and only orig_bitpix
 * says so.  Dropping it silently promoted restored 8-bit masks to 16-bit. */
Test(nde_snapstore, orig_bitpix_survives_the_round_trip) {
	fits *f = flis_test_make_mono_fits(8, 6, 0.f);
	free(f->fdata); f->fdata = NULL;
	f->type = DATA_USHORT;
	f->bitpix = USHORT_IMG;
	f->orig_bitpix = BYTE_IMG;          /* 8-bit content in a 16-bit buffer */
	f->data = calloc(48, sizeof(WORD));
	f->pdata[0] = f->pdata[1] = f->pdata[2] = f->data;

	nde_snap *s = nde_snap_create(f);
	cr_assert_not_null(s);

	fits *back = nde_snap_read(s);
	cr_assert_not_null(back);
	cr_assert_eq(back->orig_bitpix, BYTE_IMG,
	             "read must not promote 8-bit content to 16-bit");
	clearfits(back); free(back);

	fits *into = flis_test_make_mono_fits(2, 2, 0.f);
	cr_assert_eq(nde_snap_read_into(s, into), 0);
	cr_assert_eq(into->orig_bitpix, BYTE_IMG, "read_into must carry it too");
	clearfits(into); free(into);

	nde_snap_unref(s);
	clearfits(f); free(f);
}

Test(nde_snapstore, registry_ownership_semantics) {
	fits *f = make_fixture(6, 6, 0.2f);
	nde_snap *s = nde_snap_create(f);
	nde_snap_set_tag(s, -1, 7, TRUE);
	cr_assert(nde_snapstore_has(-1, 7, TRUE));
	cr_assert(!nde_snapstore_has(-1, 7, FALSE));
	cr_assert(!nde_snapstore_has(0, 7, TRUE));

	gint item; gint64 rec; gboolean post;
	cr_assert(nde_snap_tag_get(s, &item, &rec, &post));
	cr_assert_eq(item, -1);
	cr_assert_eq(rec, 7);
	cr_assert(post);

	fits *got = nde_snapstore_lookup(-1, 7, TRUE);
	cr_assert_not_null(got);
	cr_assert(memcmp(got->fdata, f->fdata, 36 * sizeof(float)) == 0);
	clearfits(got); free(got);

	/* latest same-tag snapshot wins; the older one's death must not evict
	 * the newer entry */
	fits *f2 = make_fixture(6, 6, 0.9f);
	nde_snap *s2 = nde_snap_create(f2);
	nde_snap_set_tag(s2, -1, 7, TRUE);
	nde_snap_unref(s);                         /* old snap dies */
	cr_assert(nde_snapstore_has(-1, 7, TRUE), "newer same-tag entry must survive");
	got = nde_snapstore_lookup(-1, 7, TRUE);
	cr_assert(memcmp(got->fdata, f2->fdata, 36 * sizeof(float)) == 0,
	          "lookup must return the newer snapshot");
	clearfits(got); free(got);

	/* final owner death removes the index entry */
	nde_snap_unref(s2);
	cr_assert(!nde_snapstore_has(-1, 7, TRUE));
	cr_assert_null(nde_snapstore_lookup(-1, 7, TRUE));

	clearfits(f); free(f);
	clearfits(f2); free(f2);
}

Test(nde_snapstore, pool_budget_lru_and_invalidation) {
	/* Budget of 1 MB; each 128x128 float deposit is 64 KiB. */
	com.pref.nde_cache_mb = 1;
	const guint w = 128, h = 128;

	fits *f = make_fixture(w, h, 0.f);
	/* 20 deposits x 64KiB = 1.25 MB > 1 MB → evictions of the oldest */
	for (gint64 k = 1; k <= 20; k++) {
		f->fdata[0] = (float)k;
		nde_snapstore_deposit(f, -1, k);
	}
	nde_snapstore_stats_t st;
	nde_snapstore_stats(&st);
	cr_assert_eq(st.deposits, 20);
	cr_assert(st.evictions > 0, "budget must force evictions");
	cr_assert(!nde_snapstore_has(-1, 1, TRUE), "oldest deposit must be evicted");
	cr_assert(nde_snapstore_has(-1, 20, TRUE), "newest deposit must survive");

	/* LRU touch: look up an old-ish survivor, then deposit more — the
	 * touched one must outlive its untouched neighbours */
	gint64 survivor = 0;
	for (gint64 k = 20; k >= 1; k--)
		if (nde_snapstore_has(-1, k, TRUE)) survivor = k;   /* oldest live */
	cr_assert(survivor > 0);
	fits *touched = nde_snapstore_lookup(-1, survivor, TRUE);
	clearfits(touched); free(touched);
	for (gint64 k = 21; k <= 24; k++) {
		f->fdata[0] = (float)k;
		nde_snapstore_deposit(f, -1, k);
	}
	cr_assert(nde_snapstore_has(-1, survivor, TRUE),
	          "LRU-touched entry must survive newer evictions");

	/* invalidation: an edit at record 22 drops POST(K>=22), keeps older */
	nde_snapstore_invalidate_from(-1, 22);
	cr_assert(!nde_snapstore_has(-1, 22, TRUE));
	cr_assert(!nde_snapstore_has(-1, 23, TRUE));
	cr_assert(!nde_snapstore_has(-1, 24, TRUE));
	cr_assert(nde_snapstore_has(-1, survivor, TRUE),
	          "entries before the edit stay valid");
	/* other items untouched */
	nde_snapstore_deposit(f, 5, 30);
	nde_snapstore_invalidate_from(-1, 1);
	cr_assert(nde_snapstore_has(5, 30, TRUE));

	/* purge drops everything */
	nde_snapstore_pool_purge();
	cr_assert(!nde_snapstore_has(5, 30, TRUE));
	cr_assert(!nde_snapstore_has(-1, survivor, TRUE));

	/* budget 0 disables deposits entirely */
	com.pref.nde_cache_mb = 0;
	nde_snapstore_deposit(f, -1, 40);
	cr_assert(!nde_snapstore_has(-1, 40, TRUE));

	clearfits(f); free(f);
}

/* Regression (af5ceedd8): the deposit-time over-budget loop tests
 * live_bytes, but live_bytes only falls in nde_snap_unref, which runs after
 * the mutex is dropped — so while the loop runs live_bytes never changes.
 * With the naked `live_bytes > budget` condition that is loop-invariant, so
 * one byte of overshoot drains the WHOLE pool to a single entry instead of
 * trimming to fit.  Deposit exactly one snap's worth over budget and assert
 * we shed only that one, not everything. */
Test(nde_snapstore, deposit_trims_to_fit_not_to_one) {
	/* Budget 1 MiB; each 128x128 float snap is 64 KiB → budget holds 16.
	 * Depositing 17 leaves us 64 KiB (one snap) over budget. */
	com.pref.nde_cache_mb = 1;
	const guint w = 128, h = 128;
	const gint64 snap_bytes = (gint64)w * h * sizeof(float);
	const gint64 budget = 1 * 1024 * 1024;
	const gint64 fit = budget / snap_bytes;   /* 16 entries fit */
	const gint64 n = fit + 1;                 /* 17: one snap over budget */

	fits *f = make_fixture(w, h, 0.f);
	for (gint64 k = 1; k <= n; k++) {
		f->fdata[0] = (float)k;
		nde_snapstore_deposit(f, -1, k);
	}

	nde_snapstore_stats_t st;
	nde_snapstore_stats(&st);
	cr_assert_eq(st.deposits, n);
	/* The fix trims exactly one entry; the bug drained the pool to one,
	 * evicting fit (=16) entries. */
	cr_assert_eq(st.evictions, 1,
	             "over-budget deposit must trim to fit, not drain the pool");

	/* Count survivors directly: the oldest (record 1) is the sole victim,
	 * the remaining `fit` entries must all still be present. */
	gint64 survivors = 0;
	for (gint64 k = 1; k <= n; k++)
		if (nde_snapstore_has(-1, k, TRUE)) survivors++;
	cr_assert_eq(survivors, fit,
	             "only the LRU entry should be evicted, %ld should remain", (long)fit);
	cr_assert(!nde_snapstore_has(-1, 1, TRUE), "the LRU entry is the victim");
	cr_assert(nde_snapstore_has(-1, n, TRUE), "the freshest deposit survives");
	cr_assert(nde_snapstore_has(-1, 2, TRUE),
	          "the second-oldest must NOT be drained (the bug would evict it)");

	clearfits(f); free(f);
}

Test(nde_snapstore, same_tag_deposit_replaces) {
	com.pref.nde_cache_mb = 64;
	fits *f = make_fixture(16, 16, 0.f);
	f->fdata[0] = 1.f;
	nde_snapstore_deposit(f, -1, 3);
	f->fdata[0] = 2.f;
	nde_snapstore_deposit(f, -1, 3);
	fits *got = nde_snapstore_lookup(-1, 3, TRUE);
	cr_assert_not_null(got);
	cr_assert_eq(got->fdata[0], 2.f, "same-tag deposit must replace");
	clearfits(got); free(got);
	nde_snapstore_stats_t st;
	nde_snapstore_stats(&st);
	cr_assert_eq(st.evictions, 1, "the replaced entry counts as one eviction");
	clearfits(f); free(f);
}
