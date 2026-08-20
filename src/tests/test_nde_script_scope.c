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
 * test_nde_script_scope — the phase-5 script provenance scope decision logic
 * (nde_script_scope.c).  Drives begin / mark / declare / end directly (no
 * Python needed) and checks the at-most-one-record net-effect decision:
 * read-only → nothing; net pixel change → Tier-B; previewed-then-reverted →
 * nothing; declared → Tier-C carrying script/sha256/args.
 */

#include <criterion/criterion.h>
#include <glib.h>
#include "flis_test_helpers.h"
#include "core/siril.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_script_scope.h"

cominfo com;
fits *gfit;

static void scope_test_init(void) {
	flis_test_init_com();
	com.pref.nde_cache_mb = 2048;   /* else the checkpoint pool is silently off */
	gfit = flis_test_make_mono_fits(8, 8, 0.25f);
	cr_assert_not_null(gfit);
}

static void scope_test_fini(void) {
	if (gfit) {
		clearfits(gfit);
		free(gfit);
		gfit = NULL;
	}
	flis_test_cleanup_com();
}

TestSuite(nde_script_scope, .init = scope_test_init, .fini = scope_test_fini);

/* Number of live records currently in the history. */
static guint record_count(void) {
	GPtrArray *snap = nde_history_snapshot(NULL);
	guint n = snap ? snap->len : 0;
	if (snap)
		g_ptr_array_unref(snap);
	return n;
}

/* The single live record (fails if there is not exactly one).  The snapshot is
 * a deep copy whose records are freed on unref, so we intentionally leak the
 * (short-lived, test-only) array to keep the returned record valid. */
static const nde_record *only_record(void) {
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	cr_assert_eq(snap->len, 1, "expected exactly one record, got %u", snap->len);
	return g_ptr_array_index(snap, 0);
}

/* ---- read-only: no record --------------------------------------------- */

Test(nde_script_scope, readonly_leaves_nothing) {
	nde_script_scope_begin(NULL);
	cr_assert(nde_script_scope_active());
	/* touch nothing */
	nde_script_scope_end();
	cr_assert(!nde_script_scope_active());
	cr_assert_eq(record_count(), 0, "a read-only script must leave no record");
}

/* ---- net pixel change: one Tier-B barrier ------------------------------ */

Test(nde_script_scope, pixel_change_makes_one_tier_b) {
	nde_script_scope_begin(NULL);
	gfit->fdata[0] = 0.9f;                 /* a real mutation */
	nde_script_scope_mark_pixels_dirty();
	nde_script_scope_end();

	cr_assert_eq(record_count(), 1);
	const nde_record *rec = only_record();
	cr_assert_str_eq(rec->op_id, "python.script");
	cr_assert_eq(rec->tier, NDE_TIER_B);
	cr_assert_null(rec->params, "a Tier-B barrier carries no params");
}

/* ---- previewed then reverted: nothing (net effect zero) ---------------- */

Test(nde_script_scope, preview_then_revert_leaves_nothing) {
	float original = gfit->fdata[0];
	nde_script_scope_begin(NULL);
	gfit->fdata[0] = 0.9f;                 /* "preview" write */
	nde_script_scope_mark_pixels_dirty();
	gfit->fdata[0] = original;             /* script reverts on cancel */
	nde_script_scope_end();
	cr_assert_eq(record_count(), 0,
	             "a previewed-then-reverted script must leave no record");
}

/* ---- non-pixel mutation (mask): forces a record ------------------------ */

Test(nde_script_scope, nonpixel_dirty_forces_record) {
	nde_script_scope_begin(NULL);
	/* no pixel change, but a mask write happened */
	nde_script_scope_mark_nonpixel_dirty();
	nde_script_scope_end();
	cr_assert_eq(record_count(), 1,
	             "a mask-only change is invisible to the pixel hash but must record");
	cr_assert_eq(only_record()->tier, NDE_TIER_B);
}

/* ---- declared + changed: one Tier-C carrying script/sha256/args -------- */

Test(nde_script_scope, declared_change_makes_tier_c) {
	/* A real file to hash so the declaration is honoured (Tier-C needs a
	 * script sha256). */
	gchar *dir = g_dir_make_tmp("nde_scope_XXXXXX", NULL);
	cr_assert_not_null(dir);
	gchar *path = g_build_filename(dir, "myscript.py", NULL);
	cr_assert(g_file_set_contents(path, "print('hi')\n", -1, NULL));

	nde_script_scope_begin(path);
	gfit->fdata[0] = 0.9f;
	nde_script_scope_mark_pixels_dirty();
	nde_script_scope_record_args("npoints=25;degree=4");
	nde_script_scope_end();

	cr_assert_eq(record_count(), 1);
	const nde_record *rec = only_record();
	cr_assert_str_eq(rec->op_id, "python.script");
	cr_assert_eq(rec->tier, NDE_TIER_C);
	cr_assert_not_null(rec->params);
	/* params carry the replay recipe */
	cr_assert(g_strstr_len(rec->params, -1, "script=") != NULL);
	cr_assert(g_strstr_len(rec->params, -1, "sha256=") != NULL);
	cr_assert(g_strstr_len(rec->params, -1, "args=") != NULL);
	/* summary is the script basename */
	cr_assert_str_eq(rec->summary, "myscript.py");

	g_free(path);
	g_free(dir);
}

/* ---- declared but NO net change: nothing (previewed, never committed) --- */

Test(nde_script_scope, declared_without_change_leaves_nothing) {
	gchar *dir = g_dir_make_tmp("nde_scope_XXXXXX", NULL);
	gchar *path = g_build_filename(dir, "s.py", NULL);
	g_file_set_contents(path, "x\n", -1, NULL);

	float original = gfit->fdata[0];
	nde_script_scope_begin(path);
	gfit->fdata[0] = 0.9f;
	nde_script_scope_mark_pixels_dirty();
	nde_script_scope_record_args("a=1");
	gfit->fdata[0] = original;             /* reverted */
	nde_script_scope_end();
	cr_assert_eq(record_count(), 0,
	             "no net change → no record even when a declaration was made");

	g_free(path);
	g_free(dir);
}

/* ---- nesting flattens to one scope ------------------------------------- */

Test(nde_script_scope, nested_begins_flatten_to_one_record) {
	nde_script_scope_begin(NULL);
	nde_script_scope_begin(NULL);          /* flattened onto the first */
	cr_assert(nde_script_scope_active());
	gfit->fdata[0] = 0.9f;
	nde_script_scope_mark_pixels_dirty();
	nde_script_scope_end();                /* inner: depth 2 → 1, still open */
	cr_assert(nde_script_scope_active());
	cr_assert_eq(record_count(), 0, "no record until the outermost scope ends");
	nde_script_scope_end();                /* outer: depth 1 → 0, decides */
	cr_assert(!nde_script_scope_active());
	cr_assert_eq(record_count(), 1);
}
