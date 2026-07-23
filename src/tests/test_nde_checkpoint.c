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
 * test_nde_checkpoint — the in-session baseline checkpoint store (P2.B):
 *   - ensure/get round-trip is pixel-EXACT for ushort and float
 *   - ensure is idempotent (first pre-op pixels win, later ones ignored)
 *   - exists / drop / purge behave
 *   - ensure-after-first-record ordering: the baseline reflects the pixels
 *     handed to the FIRST ensure, not a later state
 */

#include <criterion/criterion.h>
#include <glib.h>
#include <glib/gstdio.h>

#include "core/siril.h"
#include "io/image_format_fits.h"
#include "core/nde_checkpoint.h"

cominfo com;
fits *gfit;

static void setup(void) {
	memset(&com, 0, sizeof(com));
	com.max_thread = 1;
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
}

static void teardown(void) {
	nde_checkpoint_purge();
	g_free(com.pref.swap_dir);
	com.pref.swap_dir = NULL;
}

TestSuite(nde_checkpoint, .init = setup, .fini = teardown);

/* ---- fixtures ---------------------------------------------------------- */

static fits *make_float(int rx, int ry, int nch, float base) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, rx, ry, nch, DATA_FLOAT), 0);
	size_t n = (size_t)rx * ry * nch;
	for (size_t i = 0; i < n; i++)
		f->fdata[i] = base + (float)i * 1e-4f;
	return f;
}

static fits *make_ushort(int rx, int ry, int nch, WORD base) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, rx, ry, nch, DATA_USHORT), 0);
	size_t n = (size_t)rx * ry * nch;
	for (size_t i = 0; i < n; i++)
		f->data[i] = (WORD)(base + (WORD)(i * 7));
	return f;
}

static void free_fit(fits *f) {
	if (f) { clearfits(f); free(f); }
}

/* ---- tests ------------------------------------------------------------- */

Test(nde_checkpoint, roundtrip_float_exact) {
	fits *src = make_float(17, 13, 3, 0.25f);
	nde_checkpoint_baseline_ensure(src, 4);
	cr_assert(nde_checkpoint_baseline_exists(4));

	fits *got = nde_checkpoint_baseline_get(4);
	cr_assert_not_null(got);
	cr_assert_eq(got->rx, src->rx);
	cr_assert_eq(got->ry, src->ry);
	cr_assert_eq(got->naxes[2], src->naxes[2]);
	cr_assert_eq(got->type, DATA_FLOAT);
	size_t n = (size_t)src->rx * src->ry * src->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(got->fdata[i], src->fdata[i],
		             "float pixel %zu mismatch: %g != %g", i, got->fdata[i], src->fdata[i]);
	free_fit(got);
	free_fit(src);
}

Test(nde_checkpoint, roundtrip_ushort_exact) {
	fits *src = make_ushort(9, 11, 1, 100);
	nde_checkpoint_baseline_ensure(src, -1);   /* plain image key */
	cr_assert(nde_checkpoint_baseline_exists(-1));

	fits *got = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(got);
	cr_assert_eq(got->type, DATA_USHORT);
	size_t n = (size_t)src->rx * src->ry * src->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(got->data[i], src->data[i],
		             "ushort pixel %zu mismatch: %u != %u", i, got->data[i], src->data[i]);
	free_fit(got);
	free_fit(src);
}

Test(nde_checkpoint, ensure_is_idempotent) {
	fits *first = make_float(8, 8, 1, 0.1f);
	fits *later = make_float(8, 8, 1, 0.9f);
	nde_checkpoint_baseline_ensure(first, 7);   /* wins */
	nde_checkpoint_baseline_ensure(later, 7);   /* ignored */

	fits *got = nde_checkpoint_baseline_get(7);
	cr_assert_not_null(got);
	/* Baseline must reflect the FIRST ensure's pixels. */
	cr_assert_eq(got->fdata[0], first->fdata[0]);
	cr_assert_neq(got->fdata[0], later->fdata[0]);
	free_fit(got);
	free_fit(first);
	free_fit(later);
}

Test(nde_checkpoint, adopt_overwrites) {
	fits *a = make_float(8, 8, 1, 0.1f);
	fits *b = make_float(8, 8, 1, 0.9f);
	nde_checkpoint_baseline_ensure(a, 3);
	nde_checkpoint_baseline_adopt(b, 3);   /* load-path: overwrite */

	fits *got = nde_checkpoint_baseline_get(3);
	cr_assert_not_null(got);
	cr_assert_eq(got->fdata[0], b->fdata[0]);
	free_fit(got);
	free_fit(a);
	free_fit(b);
}

Test(nde_checkpoint, drop_and_purge) {
	fits *f1 = make_float(4, 4, 1, 0.2f);
	fits *f2 = make_float(4, 4, 1, 0.3f);
	nde_checkpoint_baseline_ensure(f1, 1);
	nde_checkpoint_baseline_ensure(f2, 2);
	cr_assert(nde_checkpoint_baseline_exists(1));
	cr_assert(nde_checkpoint_baseline_exists(2));

	nde_checkpoint_drop(1);
	cr_assert_not(nde_checkpoint_baseline_exists(1));
	cr_assert(nde_checkpoint_baseline_exists(2));

	nde_checkpoint_purge();
	cr_assert_not(nde_checkpoint_baseline_exists(2));
	cr_assert_null(nde_checkpoint_baseline_get(2));
	free_fit(f1);
	free_fit(f2);
}

Test(nde_checkpoint, get_absent_is_null) {
	cr_assert_not(nde_checkpoint_baseline_exists(42));
	cr_assert_null(nde_checkpoint_baseline_get(42));
}

/* Ordering: ensure captures the pre-FIRST-op pixels even though the source
 * buffer is subsequently mutated in place (mirrors the worker `orig` being
 * reused/cleared after append). */
Test(nde_checkpoint, captures_pre_first_op_state) {
	fits *pre = make_float(6, 6, 1, 0.4f);
	float saved0 = pre->fdata[0];
	nde_checkpoint_baseline_ensure(pre, 5);
	/* Mutate the source AFTER ensure — the swap copy must be independent. */
	size_t n = (size_t)pre->rx * pre->ry;
	for (size_t i = 0; i < n; i++)
		pre->fdata[i] = 0.99f;

	fits *got = nde_checkpoint_baseline_get(5);
	cr_assert_not_null(got);
	cr_assert_eq(got->fdata[0], saved0);
	cr_assert_neq(got->fdata[0], pre->fdata[0]);
	free_fit(got);
	free_fit(pre);
}
