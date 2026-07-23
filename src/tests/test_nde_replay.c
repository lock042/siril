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
