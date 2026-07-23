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
#include "core/processing_thread.h"
#include "core/op_descriptors.h"
#include "core/op_descriptor.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "algos/geometry.h"
#include "filters/asinh.h"
#include "filters/mtf.h"
#include "filters/ght.h"
#include "filters/scnr.h"
#include "filters/median.h"
#include "filters/curve_transform.h"
#include "algos/colors.h"
#include "algos/background_extraction.h"

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

/* ===================================================================== *
 *  P2.E — per-op golden round-trips through the REAL capture path.       *
 *                                                                        *
 *  Each test: build a fixture with structured content → apply the op via *
 *  generic_image_worker() on gfit exactly as a headless command does     *
 *  (swap path fires the NDE capture + baseline checkpoint) → replay the  *
 *  chain from the baseline into a scratch fits → compare the replay      *
 *  result against the post-op pixels.  These run the same code on the    *
 *  same pixels in-process, so bit-exactness is expected for every        *
 *  deterministic hook; only BGE (below) is tolerance-only, and only      *
 *  because its GSL polynomial fit is not bit-reproducible.               *
 * ===================================================================== */

/* Fill a mono fits with a smooth, asymmetric gradient + a diagonal ripple
 * so that mirrors/rotations/stretches all have observable, order-sensitive
 * effects (a constant field would make several ops the identity). */
static void fill_mono_gradient(fits *f) {
	for (unsigned y = 0; y < f->ry; y++)
		for (unsigned x = 0; x < f->rx; x++) {
			float v = 0.05f + 0.6f * ((float)x / f->rx)
			               + 0.3f * ((float)y / f->ry)
			               + 0.03f * sinf(0.5f * (x + y));
			if (v < 0.f) v = 0.f;
			if (v > 1.f) v = 1.f;
			f->fdata[x + (size_t)y * f->rx] = v;
		}
}

/* Fill an RGB fits with three distinct per-channel gradients (needed by the
 * ops that mix channels — scnr, ccm). */
static void fill_rgb_gradient(fits *f) {
	size_t n = (size_t)f->rx * f->ry;
	for (unsigned y = 0; y < f->ry; y++)
		for (unsigned x = 0; x < f->rx; x++) {
			size_t i = x + (size_t)y * f->rx;
			f->fpdata[RLAYER][i] = 0.10f + 0.50f * ((float)x / f->rx);
			f->fpdata[GLAYER][i] = 0.20f + 0.40f * ((float)y / f->ry);
			f->fpdata[BLAYER][i] = 0.15f + 0.45f * ((float)(x + y) / (f->rx + f->ry));
		}
	(void)n;
}

/* Bit-exact whole-buffer compare of two fits with identical geometry. */
static void assert_pixels_bit_exact(const fits *a, const fits *b, const char *what) {
	cr_assert_eq(a->rx, b->rx, "%s: rx mismatch (%u vs %u)", what, a->rx, b->rx);
	cr_assert_eq(a->ry, b->ry, "%s: ry mismatch (%u vs %u)", what, a->ry, b->ry);
	cr_assert_eq(a->naxes[2], b->naxes[2], "%s: channel count mismatch", what);
	size_t n = (size_t)a->rx * a->ry * a->naxes[2];
	cr_assert(memcmp(a->fdata, b->fdata, n * sizeof(float)) == 0,
	          "%s: replay must reproduce the post-op pixels bit-exactly", what);
}

/* Max-abs / mean-abs deviation between two fits with identical geometry
 * (used for the tolerance-only BGE case). */
static void assert_pixels_close(const fits *a, const fits *b, float tol, const char *what) {
	cr_assert_eq(a->rx, b->rx, "%s: rx mismatch", what);
	cr_assert_eq(a->ry, b->ry, "%s: ry mismatch", what);
	cr_assert_eq(a->naxes[2], b->naxes[2], "%s: channel count mismatch", what);
	size_t n = (size_t)a->rx * a->ry * a->naxes[2];
	float maxdev = 0.f;
	double sumdev = 0.0;
	for (size_t i = 0; i < n; i++) {
		float d = fabsf(a->fdata[i] - b->fdata[i]);
		if (d > maxdev) maxdev = d;
		sumdev += d;
	}
	cr_assert(maxdev <= tol,
	          "%s: replay deviation %.3g exceeds tolerance %.3g (mean %.3g)",
	          what, (double)maxdev, (double)tol, sumdev / n);
}

/*
 * Apply @op to @gfit through the real headless-command capture path.
 *
 * @user_heap must be a heap-allocated params struct whose FIRST member is a
 * destructor pointer (destroy_fn), exactly like every command construction
 * site: the worker owns it and frees it (destructor-first) along with args.
 * We set com.headless = TRUE so the worker's completion tail runs
 * stop_processing_thread() + free_generic_img_args(args) itself, never
 * touching a GUI idle.  (stop_processing_thread() sets the sticky cancel
 * flag; callers clear it via reserve_thread() before replaying.)
 *
 * gfit is the same fits struct after the swap, now holding the post-op
 * pixels — the caller compares the replay result against it.
 */
static int apply_op_real(const op_descriptor *op, gpointer user_heap) {
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	cr_assert_not_null(args);
	args->fit = gfit;                 /* swap path → capture + baseline fire */
	args->op = op;
	args->user = user_heap;
	args->command = TRUE;             /* headless-command completion tail */
	args->command_updates_gfit = TRUE;
	args->max_threads = 1;            /* deterministic, single-threaded */
	args->mem_ratio = -1.0f;          /* skip the memory check in the harness */
	args->verbose = FALSE;
	gboolean prev_headless = com.headless;
	com.headless = TRUE;              /* worker owns + frees args */
	int rc = GPOINTER_TO_INT(generic_image_worker(args));
	com.headless = prev_headless;
	return rc;
}

/* Build the item −1 chain, replay it, and return the scratch result.  Clears
 * the sticky cancel flag left by the capture path's stop_processing_thread().
 * Asserts a replayable single-record chain unless @expect_records says
 * otherwise. */
static fits *replay_current_chain(guint expect_records) {
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable, "chain not replayable: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert_eq(chain->records->len, expect_records,
	             "expected %u records, got %u", expect_records, chain->records->len);
	/* clear the cancel_flag the capture path set (stop_processing_thread);
	 * reserve_thread() is the documented reset and also claims the slot the
	 * replay driver expects to own. */
	cr_assert(reserve_thread(), "could not reserve the processing slot");
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	g_free(err);
	nde_chain_free(chain);
	return result;
}

/* Teardown shared by every per-op golden test: drop the scratch, the chain
 * records + baselines (attach(NULL) purges checkpoints too), and gfit. */
static void golden_teardown(fits *result, fits *f) {
	if (result) { clearfits(result); free(result); }
	nde_history_attach(NULL);      /* also purges the checkpoint store */
	clearfits(f); free(f);
	gfit = NULL;
}

/* ---- geometry: the six ops ---- */

Test(nde_replay, golden_mirrorx) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct mirror_args *u = calloc(1, sizeof(*u));   /* POD → destroy_fn NULL */
	u->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "mirrorx");
	golden_teardown(result, f);
}

Test(nde_replay, golden_mirrory) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct mirror_args *u = calloc(1, sizeof(*u));
	u->x_axis = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_mirrory, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "mirrory");
	golden_teardown(result, f);
}

Test(nde_replay, golden_crop) {
	fits *f = flis_test_make_mono_fits(20, 16, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct crop_args *u = calloc(1, sizeof(*u));
	u->area = (rectangle){ .x = 3, .y = 2, .w = 12, .h = 9 };  /* non-degenerate */
	cr_assert_eq(apply_op_real(&op_desc_crop, u), 0);
	cr_assert_eq(gfit->rx, 12); cr_assert_eq(gfit->ry, 9);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "crop");
	golden_teardown(result, f);
}

Test(nde_replay, golden_binning) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct binning_args *u = calloc(1, sizeof(*u));
	u->factor = 2;
	u->mean = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_binning, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "binning");
	golden_teardown(result, f);
}

Test(nde_replay, golden_resample) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct resample_args *u = calloc(1, sizeof(*u));
	u->toX = 24; u->toY = 18;
	u->interpolation = OPENCV_CUBIC;
	u->clamp = TRUE;
	u->update_wcs = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_resample, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "resample");
	golden_teardown(result, f);
}

/* Rotation: pick 90° with area == whole image → verbose_rotate_fast(), a pure
 * lossless 90° transpose (no interpolation), so replay is BIT-EXACT.  A
 * non-axis-aligned angle would go through OpenCV interpolation which is
 * deterministic here too, but 90° removes any doubt about interpolation
 * reproducibility across the two in-process runs. */
Test(nde_replay, golden_rotation_90) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct rotation_args *u = calloc(1, sizeof(*u));
	u->area = (rectangle){ .x = 0, .y = 0, .w = 16, .h = 12 };
	u->angle = 90.0;
	u->interpolation = OPENCV_LANCZOS4;   /* ignored by the fast path */
	u->cropped = 1;
	u->clamp = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_rotation, u), 0);
	cr_assert_eq(gfit->rx, 12); cr_assert_eq(gfit->ry, 16);  /* transposed */
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "rotation 90");
	golden_teardown(result, f);
}

/* ---- stretches / filters ---- */

Test(nde_replay, golden_asinh) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	asinh_params *u = calloc(1, sizeof(*u));    /* destroy_fn (first) NULL: POD */
	u->beta = 15.0f;
	u->offset = 0.02f;
	u->human_luminance = FALSE;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "asinh");
	golden_teardown(result, f);
}

Test(nde_replay, golden_mtf) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	/* mtf owns a create/destroy pair; mirror the command site (data->fit and
	 * auto_display_compensation are runtime — the latter forced FALSE so the
	 * hook never reaches the ICC-compensation gui path). */
	struct mtf_data *u = create_mtf_data();
	cr_assert_not_null(u);
	u->linked = TRUE;
	u->params.midtones = 0.35f;
	u->params.shadows = 0.02f;
	u->params.highlights = 0.98f;
	u->params.do_red = u->params.do_green = u->params.do_blue = TRUE;
	u->auto_display_compensation = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_mtf, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "mtf");
	golden_teardown(result, f);
}

Test(nde_replay, golden_ghs) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	/* ghs: params_ght is a pointer owned by ght_data — build exactly like the
	 * command site (create_ght_data() + heap ght_params) and free via
	 * destroy_ght_data (the struct's destroy_fn). */
	ght_params *p = calloc(1, sizeof(*p));
	cr_assert_not_null(p);
	p->B = 1.0f; p->D = 1.2f; p->LP = 0.0f; p->SP = 0.4f; p->HP = 1.0f; p->BP = 0.0f;
	p->stretchtype = STRETCH_PAYNE_NORMAL;
	p->payne_colourstretchmodel = COL_INDEP;
	p->do_red = p->do_green = p->do_blue = TRUE;
	p->clip_mode = CLIP;
	struct ght_data *u = create_ght_data();
	cr_assert_not_null(u);
	u->params_ght = p;
	u->auto_display_compensation = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_ghs, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "ghs");
	golden_teardown(result, f);
}

Test(nde_replay, golden_median) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	struct median_filter_data *u = calloc(1, sizeof(*u));    /* POD */
	u->ksize = 3;
	u->amount = 0.8;
	u->iterations = 1;
	cr_assert_eq(apply_op_real(&op_desc_median, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "median");
	golden_teardown(result, f);
}

/* ---- ops that require 3 channels ---- */

Test(nde_replay, golden_scnr) {
	fits *f = flis_test_make_rgb_fits(16, 12, 0.f, 0.f, 0.f);
	fill_rgb_gradient(f);
	gfit = f;
	struct scnr_data *u = calloc(1, sizeof(*u));    /* POD */
	u->type = SCNR_AVERAGE_NEUTRAL;
	u->amount = 0.5;
	u->preserve = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_scnr, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "scnr");
	golden_teardown(result, f);
}

Test(nde_replay, golden_ccm) {
	fits *f = flis_test_make_rgb_fits(16, 12, 0.f, 0.f, 0.f);
	fill_rgb_gradient(f);
	gfit = f;
	struct ccm_data *u = calloc(1, sizeof(*u));     /* POD */
	/* a mild, invertible-ish colour matrix + gamma */
	u->matrix[0][0] = 1.05f; u->matrix[0][1] = -0.03f; u->matrix[0][2] = 0.02f;
	u->matrix[1][0] = 0.01f; u->matrix[1][1] = 0.98f;  u->matrix[1][2] = 0.04f;
	u->matrix[2][0] = -0.02f; u->matrix[2][1] = 0.03f; u->matrix[2][2] = 1.02f;
	u->power = 1.1f;
	cr_assert_eq(apply_op_real(&op_desc_ccm, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "ccm");
	golden_teardown(result, f);
}

/* ---- curves: GUI-only construction, core hook + serializer ---- */

Test(nde_replay, golden_curves) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	/* Build curve_params exactly as the GUI does: new_curve_params() sets
	 * destroy_fn = free_curve_params (which does NOT free the points list — the
	 * caller owns it), and we attach a points GList we free ourselves after the
	 * worker.  Mirrors src/gui-gtk4/curves.c curves_process_with_worker(). */
	struct curve_params *u = new_curve_params();
	cr_assert_not_null(u);
	GList *pts = NULL;
	point *p0 = g_new(point, 1); p0->x = 0.0; p0->y = 0.0; pts = g_list_append(pts, p0);
	point *p1 = g_new(point, 1); p1->x = 0.4; p1->y = 0.6; pts = g_list_append(pts, p1);
	point *p2 = g_new(point, 1); p2->x = 1.0; p2->y = 1.0; pts = g_list_append(pts, p2);
	u->points = pts;
	u->algorithm = LINEAR;
	u->do_channel[0] = u->do_channel[1] = u->do_channel[2] = TRUE;
	u->fit = gfit;
	u->verbose = FALSE;
	u->for_preview = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_curves, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "curves");
	golden_teardown(result, f);
	g_list_free_full(pts, g_free);   /* the caller owns the points list */
}

/* ---- multi-op chain: three different ops stacked ---- */

Test(nde_replay, golden_chain_multi_op) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);

	asinh_params *a = calloc(1, sizeof(*a));
	a->beta = 8.0f; a->offset = 0.01f; a->human_luminance = FALSE; a->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, a), 0);

	struct median_filter_data *md = calloc(1, sizeof(*md));
	md->ksize = 3; md->amount = 0.5; md->iterations = 1;
	cr_assert_eq(apply_op_real(&op_desc_median, md), 0);

	fits *result = replay_current_chain(3);
	assert_pixels_bit_exact(result, gfit, "chain(mirrorx+asinh+median)");
	golden_teardown(result, f);
}

/* ---- BGE: tolerance-only + replay_pre samples reinstall ----
 *
 * Background extraction is not bit-reproducible: its GSL least-squares
 * polynomial fit is only reproducible to floating-point rounding, so this is
 * the one op the plan flags as tolerance-only (<=1e-5) even with dither=0.
 * Additionally, its effective sample positions live in com.grad_samples, not
 * in the params struct — replay_pre (op_desc_remove_gradient) reinstalls the
 * recorded samples before the hook runs; we assert com.grad_samples is
 * non-NULL after replay to prove that pathway executed.
 *
 * We set from_ui = TRUE (the GUI apply path).  This matters: the sample-based
 * hook FREES com.grad_samples after running when from_ui == FALSE, which
 * happens BEFORE the worker's capture block serializes the record — so on the
 * command path (subsky, from_ui == FALSE) the samples key is never captured.
 * See the report: that command-path gap is a genuine pre-existing limitation
 * flagged for the maintainer, not fixed here. */
Test(nde_replay, golden_bge_tolerance_and_samples_reinstall) {
	/* 96x96 leaves room for full 25x25 sample windows away from the border. */
	fits *f = flis_test_make_mono_fits(96, 96, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* Place a grid of valid samples on the pre-op pixels (>=12px from every
	 * edge so get_sample() accepts the full window). */
	int coords[3] = { 24, 48, 72 };
	for (int yi = 0; yi < 3; yi++)
		for (int xi = 0; xi < 3; xi++) {
			point pt = { coords[xi], coords[yi] };
			com.grad_samples = add_background_sample(com.grad_samples, gfit, pt, FALSE);
		}
	cr_assert_not_null(com.grad_samples);

	struct background_data *u = calloc(1, sizeof(*u));
	u->destroy_fn = free_background_data;
	u->method = BACKGROUND_METHOD_SAMPLES;
	u->nb_of_samples = 9;
	u->tolerance = 2.0;
	u->correction = BACKGROUND_CORRECTION_SUBTRACT;
	u->interpolation_method = BACKGROUND_INTER_POLY;
	u->degree = BACKGROUND_POLY_1;        /* needs >=3 samples; we have 9 */
	u->smoothing = 0.5;
	u->threads = 1;
	u->dither = FALSE;                    /* deterministic input conversion */
	u->from_ui = TRUE;                    /* keep com.grad_samples for capture */
	u->randomize = FALSE;
	u->grad_descent = FALSE;

	cr_assert_eq(apply_op_real(&op_desc_remove_gradient, u), 0);

	/* The capture (from_ui == TRUE kept com.grad_samples alive) must have
	 * serialized the effective sample positions into the record — this is the
	 * external input replay_pre reinstalls.  Verify it is present. */
	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert(strstr(rec->params, "samples=") != NULL,
	          "the captured BGE record must carry the effective sample positions");
	g_ptr_array_unref(snap);

	/* Drop the live sample list so the ONLY way the replayed hook can build a
	 * non-empty model is via replay_pre reinstalling the captured positions. */
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;

	fits *result = replay_current_chain(1);
	/* Replay SUCCEEDING is itself proof the replay_pre pathway ran: with no
	 * samples in com.grad_samples the sample-based hook's polynomial fit would
	 * fail (n < nbParam) and the replay would return NULL.  The replayed
	 * record deserializes with from_ui == FALSE, so the hook frees
	 * com.grad_samples again after applying — hence we do NOT assert it is
	 * non-NULL post-replay (it is intentionally cleared by the hook); the
	 * successful bounded-deviation reconstruction below is the observable. */
	assert_pixels_close(result, gfit, 1e-5f, "bge");
	/* Tolerance 1e-5 (not bit-exact): BGE's GSL least-squares polynomial fit
	 * is reproducible only to floating-point rounding across the two runs. */

	if (com.grad_samples) {
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
	}
	golden_teardown(result, f);
}

Test(nde_replay, golden_bge_command_path_captures_samples) {
	/* Regression for the P2.E finding: with from_ui == FALSE (the command
	 * path) the hook frees com.grad_samples on success BEFORE the worker's
	 * capture serializes the record.  The hook now stashes the effective
	 * positions into args->effective_samples first, so the record still
	 * carries them and replay reproduces the pixels. */
	fits *f = flis_test_make_mono_fits(96, 96, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	int coords[3] = { 24, 48, 72 };
	for (int yi = 0; yi < 3; yi++)
		for (int xi = 0; xi < 3; xi++) {
			point pt = { coords[xi], coords[yi] };
			com.grad_samples = add_background_sample(com.grad_samples, gfit, pt, FALSE);
		}
	cr_assert_not_null(com.grad_samples);

	struct background_data *u = calloc(1, sizeof(*u));
	u->destroy_fn = free_background_data;
	u->method = BACKGROUND_METHOD_SAMPLES;
	u->nb_of_samples = 9;
	u->tolerance = 2.0;
	u->correction = BACKGROUND_CORRECTION_SUBTRACT;
	u->interpolation_method = BACKGROUND_INTER_POLY;
	u->degree = BACKGROUND_POLY_1;
	u->smoothing = 0.5;
	u->threads = 1;
	u->dither = FALSE;
	u->from_ui = FALSE;                   /* command path: hook clears samples */
	u->randomize = FALSE;
	u->grad_descent = FALSE;

	cr_assert_eq(apply_op_real(&op_desc_remove_gradient, u), 0);
	cr_assert_null(com.grad_samples, "command path must still clear the live list");

	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert(rec->params && strstr(rec->params, "samples=") != NULL,
	          "command-path BGE record must carry the effective sample positions");
	g_ptr_array_unref(snap);

	fits *result = replay_current_chain(1);
	assert_pixels_close(result, gfit, 1e-5f, "bge-cmd");

	golden_teardown(result, f);
}
