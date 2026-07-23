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
 * test_nde_replay — NDE replay engine (phase 2, nde-phase2-3-plan.md).
 * P2.A: the nde_replay flag lets the replay driver call
 * generic_image_worker() directly on a PRIVATE fits with zero side
 * effects: no undo entry, no provenance record, no HISTORY card, no
 * idles/completion machinery (the caller owns args).
 * P2.E extends this file with per-op golden round-trips.
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "algos/geometry.h"

cominfo com;
fits *gfit;

TestSuite(nde_replay, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

Test(nde_replay, worker_flag_runs_headless_with_no_side_effects) {
	fits *f = flis_test_make_mono_fits(4, 4, 0.0f);
	cr_assert_not_null(f);
	f->fdata[0] = 1.0f;

	struct mirror_args ma = { 0 };
	ma.x_axis = TRUE;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = f;                 /* private fits — never gfit in replay */
	args->op = &op_desc_mirrorx;
	args->user = &ma;
	args->nde_replay = TRUE;
	args->mem_ratio = -1.0f;       /* no memory check in the test harness */
	args->max_threads = 1;

	gpointer ret = generic_image_worker(args);
	cr_assert_eq(GPOINTER_TO_INT(ret), 0, "replay invocation failed");

	/* the op really ran on the private fits */
	cr_assert_eq(f->fdata[0], 0.0f, "marker pixel should have moved");
	float sum = 0.f;
	for (int i = 0; i < 16; i++)
		sum += f->fdata[i];
	cr_assert_float_eq(sum, 1.0f, 1e-6, "mirror must be a pure permutation");

	/* zero WORKER side effects.  (fit->history is not asserted: several
	 * hooks — the geometry ops included — append their own FITS HISTORY
	 * entry inside the algo, which is op behaviour, not worker machinery,
	 * and is correct for replay: a replayed chain rebuilds HISTORY.) */
	cr_assert_eq(nde_history_live_count(), 0, "replay must not append records");
	cr_assert_null(com.undo_stack, "replay must not save undo");

	/* the worker freed nothing on the replay path — args is still ours */
	cr_assert_eq(args->retval, 0);
	free(args);
	clearfits(f);
	free(f);
}

Test(nde_replay, worker_flag_reports_hook_failure) {
	/* crop with an empty area fails inside the hook; the flag must hand
	 * the failure back without any completion machinery. */
	fits *f = flis_test_make_mono_fits(4, 4, 0.5f);
	struct crop_args ca = { 0 };   /* zero-size area -> hook error */

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = f;
	args->op = &op_desc_crop;
	args->user = &ca;
	args->nde_replay = TRUE;
	args->mem_ratio = -1.0f;
	args->max_threads = 1;

	gpointer ret = generic_image_worker(args);
	cr_assert_neq(GPOINTER_TO_INT(ret), 0, "zero-size crop should fail");
	cr_assert_eq(nde_history_live_count(), 0);
	cr_assert_null(com.undo_stack);

	free(args);
	clearfits(f);
	free(f);
}

/* ---------------- P2.D: chain build / validate / replay ---------------- */

/* Build a Tier-A chain artificially: baseline + captured records, exactly
 * what the capture sites produce, without driving the full gfit swap path. */

Test(nde_replay, chain_replays_bit_exact) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.0f);
	for (int i = 0; i < 64; i++)
		f->fdata[i] = i / 64.0f;    /* asymmetric so mirrors matter */

	nde_checkpoint_baseline_ensure(f, -1);
	struct mirror_args ma = { 0 };
	ma.x_axis = TRUE;
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m1") > 0);
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m2") > 0);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable, "all-Tier-A chain with baseline must be replayable (first reason: %s)",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert_eq(chain->records->len, 2);

	gchar *errmsg = NULL;
	fits *result = nde_chain_replay(chain, &errmsg);
	cr_assert_not_null(result, "replay failed: %s", errmsg ? errmsg : "?");
	/* mirrorx twice is the identity: the result must equal the baseline
	 * (== f, untouched) bit-exactly */
	cr_assert_eq(result->rx, 8);
	cr_assert(memcmp(result->fdata, f->fdata, 64 * sizeof(float)) == 0,
	          "mirror+mirror must reproduce the baseline bit-exactly");
	clearfits(result); free(result);
	nde_chain_free(chain);

	/* single-record chain equals a direct application of the op */
	nde_history_on_undo(2);   /* retire m2; live chain = [m1] */
	chain = nde_chain_build(-1);
	cr_assert(chain->replayable);
	cr_assert_eq(chain->records->len, 1);
	result = nde_chain_replay(chain, &errmsg);
	cr_assert_not_null(result, "replay failed: %s", errmsg ? errmsg : "?");
	fits *expected = calloc(1, sizeof(fits));
	copyfits(f, expected, CP_DEEPCOPY | CP_ALLOC, -1);
	mirrorx(expected, FALSE);
	cr_assert(memcmp(result->fdata, expected->fdata, 64 * sizeof(float)) == 0,
	          "single-op replay must equal direct application");
	clearfits(expected); free(expected);
	clearfits(result); free(result);
	nde_chain_free(chain);

	nde_history_attach(NULL);
	clearfits(f); free(f);
}

Test(nde_replay, chain_blockers_are_reported) {
	fits *f = flis_test_make_mono_fits(4, 4, 0.5f);

	/* opaque record → not replayable */
	nde_checkpoint_baseline_ensure(f, -1);
	nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "opaque");
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert_eq(chain->reasons->len, 1);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "opaque") != NULL);
	nde_chain_free(chain);
	nde_history_attach(NULL);

	/* mask-active record → not replayable */
	nde_checkpoint_baseline_ensure(f, -1);
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("geometry.mirrorx");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("x_axis=1");
	rec->mask_active = TRUE;
	nde_history_append(rec);
	chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "mask") != NULL);
	nde_chain_free(chain);
	nde_history_attach(NULL);

	/* missing baseline → not replayable */
	struct mirror_args ma = { 0 };
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m");
	nde_checkpoint_purge();
	chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "baseline") != NULL);
	nde_chain_free(chain);
	nde_history_attach(NULL);

	/* unknown op id → not replayable */
	nde_checkpoint_baseline_ensure(f, -1);
	rec = nde_record_new();
	rec->op_id = g_strdup("future.unknown_op");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	nde_history_append(rec);
	chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "unknown") != NULL);
	nde_chain_free(chain);
	nde_history_attach(NULL);

	/* empty history: trivially "replayable", nothing to do */
	chain = nde_chain_build(-1);
	cr_assert(chain->replayable);
	cr_assert_eq(chain->records->len, 0);
	nde_chain_free(chain);

	clearfits(f); free(f);
}
