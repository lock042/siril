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
#include "filters/epf.h"
#include "filters/curve_transform.h"
#include "algos/colors.h"
#include "algos/background_extraction.h"
#include "filters/synthstar.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "io/image_format_fits.h"

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
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m1", NULL) > 0);
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m2", NULL) > 0);

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
	nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "opaque", NULL);
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
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", NULL);
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

/* ---------------- P3.B: amend / delete commit machinery ---------------- */

Test(nde_replay, amend_execute_recomputes_pixels_and_log) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u = calloc(1, sizeof(*u));
	u->beta = 15.0f;
	u->offset = 0.02f;
	u->human_luminance = FALSE;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);

	gint64 rec_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		rec_id = ((nde_record *)g_ptr_array_index(snap, 0))->record_id;
		g_ptr_array_unref(snap);
	}
	fits before = { 0 };
	copyfits(gfit, &before, CP_DEEPCOPY | CP_ALLOC, -1);

	/* a fake undo entry that the amend must flush (no meta-undo) */
	historic *fake = g_new0(historic, 1);
	fake->fd = -1;
	fake->mask_fd = -1;
	com.undo_stack = g_list_prepend(com.undo_stack, fake);

	/* amend beta 15 → 40 */
	cr_assert(reserve_thread());
	gchar *err = NULL;
	gboolean ok = nde_amend_execute(rec_id,
			"beta=40;offset=0.0199999996;human=0;clip_mode=1", &err);
	unreserve_thread();
	cr_assert(ok, "amend failed: %s", err ? err : "?");
	cr_assert_null(com.undo_stack, "amend must flush the undo stack");

	/* pixels actually changed, and the log carries the new params */
	size_t n = (size_t)gfit->rx * gfit->ry;
	gboolean changed = FALSE;
	for (size_t i = 0; i < n && !changed; i++)
		changed = (gfit->fdata[i] != before.fdata[i]);
	cr_assert(changed, "amended params must change the pixels");
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		nde_record *rec = g_ptr_array_index(snap, 0);
		cr_assert_eq(rec->record_id, rec_id);
		cr_assert(strstr(rec->params, "beta=40") != NULL, "params: %s", rec->params);
		g_ptr_array_unref(snap);
	}

	/* self-consistency: replaying the amended chain reproduces the new
	 * pixels bit-exactly */
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "amended-chain");

	clearfits(&before);
	golden_teardown(result, f);
}

Test(nde_replay, delete_execute_removes_step_from_pixels_and_log) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u = calloc(1, sizeof(*u));
	u->beta = 15.0f;
	u->offset = 0.02f;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);
	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);

	gint64 asinh_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		asinh_id = ((nde_record *)g_ptr_array_index(snap, 0))->record_id;
		g_ptr_array_unref(snap);
	}

	/* expected pixels after deleting the asinh step: mirrorx(baseline) */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	mirrorx(expected, FALSE);

	cr_assert(reserve_thread());
	gchar *err = NULL;
	gboolean ok = nde_delete_execute(asinh_id, &err);
	unreserve_thread();
	cr_assert(ok, "delete failed: %s", err ? err : "?");

	cr_assert_eq(nde_history_live_count(), 1, "one record must remain");
	assert_pixels_bit_exact(gfit, expected, "post-delete");

	clearfits(expected);
	free(expected);
	golden_teardown(NULL, f);
}

Test(nde_replay, amend_failures_leave_everything_untouched) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u = calloc(1, sizeof(*u));
	u->beta = 15.0f;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);
	gint64 rec_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		rec_id = ((nde_record *)g_ptr_array_index(snap, 0))->record_id;
		g_ptr_array_unref(snap);
	}
	fits before = { 0 };
	copyfits(gfit, &before, CP_DEEPCOPY | CP_ALLOC, -1);

	cr_assert(reserve_thread());
	gchar *err = NULL;
	/* unknown record */
	cr_assert(!nde_amend_execute(999, "beta=40", &err));
	g_clear_pointer(&err, g_free);
	/* unparsable params */
	cr_assert(!nde_amend_execute(rec_id, "nonsense=only", &err));
	cr_assert_not_null(err);
	g_clear_pointer(&err, g_free);
	/* delete of unknown record */
	cr_assert(!nde_delete_execute(999, &err));
	g_clear_pointer(&err, g_free);
	unreserve_thread();

	/* nothing changed */
	assert_pixels_bit_exact(gfit, &before, "untouched-pixels");
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 1);
	cr_assert(strstr(((nde_record *)g_ptr_array_index(snap, 0))->params, "beta=15") != NULL);
	g_ptr_array_unref(snap);

	clearfits(&before);
	golden_teardown(NULL, f);
}

Test(nde_replay, deleting_opaque_record_regains_editability) {
	/* [asinh, opaque freehand, mirrorx]: blocked; delete the opaque step
	 * and the chain becomes fully editable again, with the freehand
	 * contribution removed by construction. */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u = calloc(1, sizeof(*u));
	u->beta = 15.0f;
	u->offset = 0.02f;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);

	/* simulate a python freehand edit: poke pixels + opaque record, the
	 * same shape the set_pixeldata handler produces */
	gfit->fdata[5] = 0.999f;
	gfit->fdata[6] = 0.001f;
	gint64 opaque_id = nde_capture_opaque("python.set_pixeldata",
	                                      NDE_SCOPE_LAYER, -1, "freehand", NULL);
	cr_assert(opaque_id > 0);

	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);

	/* blocked while the opaque record lives */
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	nde_chain_free(chain);

	/* expected post-delete pixels: mirrorx(asinh(baseline)) */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	{
		asinh_params *eu = calloc(1, sizeof(*eu));
		eu->beta = 15.0f; eu->offset = 0.02f; eu->clip_mode = RESCALE;
		struct generic_img_args *eargs = calloc(1, sizeof(*eargs));
		eargs->fit = expected;
		eargs->op = &op_desc_asinh;
		eargs->user = eu;
		eargs->nde_replay = TRUE;
		eargs->mem_ratio = -1.0f;
		eargs->max_threads = 1;
		cr_assert_eq(GPOINTER_TO_INT(generic_image_worker(eargs)), 0);
		free_generic_img_args(eargs);
		mirrorx(expected, FALSE);
	}

	cr_assert(reserve_thread());
	gchar *err = NULL;
	gboolean ok = nde_delete_execute(opaque_id, &err);
	unreserve_thread();
	cr_assert(ok, "deleting the opaque record failed: %s", err ? err : "?");

	assert_pixels_bit_exact(gfit, expected, "post-opaque-delete");

	/* the chain is replayable again, and upstream records are editable */
	chain = nde_chain_build(-1);
	cr_assert(chain->replayable, "chain must be replayable after the opaque delete");
	cr_assert_eq(chain->records->len, 2);
	nde_chain_free(chain);

	gint64 asinh_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		asinh_id = ((nde_record *)g_ptr_array_index(snap, 0))->record_id;
		g_ptr_array_unref(snap);
	}
	cr_assert(reserve_thread());
	ok = nde_amend_execute(asinh_id, "beta=40;offset=0.0199999996;human=0;clip_mode=1", &err);
	unreserve_thread();
	cr_assert(ok, "amending upstream of the deleted opaque failed: %s", err ? err : "?");

	clearfits(expected);
	free(expected);
	golden_teardown(NULL, f);
}

Test(nde_replay, compositing_state_records_do_not_block_chains) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u = calloc(1, sizeof(*u));
	u->beta = 15.0f;
	u->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u), 0);

	/* a property commit as Part A captures it: Tier A, LAYER scope, no
	 * descriptor — it must be neither a member nor a blocker */
	nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER, -1,
	                       g_strdup("opacity=0.5"), "opacity");

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable,
	          "property records must not block the chain (first reason: %s)",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert_eq(chain->records->len, 1, "property records are not chain members");
	nde_chain_free(chain);

	/* and deleting a property record is refused with a clear message */
	gint64 prop_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		prop_id = ((nde_record *)g_ptr_array_index(snap, 1))->record_id;
		g_ptr_array_unref(snap);
	}
	cr_assert(reserve_thread());
	gchar *err = NULL;
	cr_assert(!nde_delete_execute(prop_id, &err));
	unreserve_thread();
	cr_assert(strstr(err, "property") != NULL, "got: %s", err);
	g_free(err);

	golden_teardown(NULL, f);
}

Test(nde_replay, document_scope_pixel_ops_block_chains) {
	/* A DOCUMENT-scope record that is not a known structural op mutated
	 * pixels document-wide (FLIS icc.convert) — it must fail closed as a
	 * chain blocker, or an amend spanning it would replay without the
	 * transform and commit wrong pixels. */
	fits *f = flis_test_make_mono_fits(8, 8, 0.5f);
	nde_checkpoint_baseline_ensure(f, -1);
	struct mirror_args ma = { 0 };
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", NULL);
	nde_capture_opaque("icc.convert", NDE_SCOPE_DOCUMENT, -1, "Converted profile", NULL);
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m2", NULL);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable,
	          "a document-wide pixel record must block the chain");
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "whole document") != NULL);
	nde_chain_free(chain);

	/* and it is NOT deletable: removing a document-wide pixel record
	 * would require recomputing EVERY layer's chain, but the delete
	 * machinery recomputes only the target item — cross-layer
	 * consistency demands the refusal (multi-layer recompute is
	 * phase-4+ territory). */
	gint64 doc_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		doc_id = ((nde_record *)g_ptr_array_index(snap, 1))->record_id;
		g_ptr_array_unref(snap);
	}
	gfit = f;
	cr_assert(reserve_thread());
	gchar *err = NULL;
	gboolean ok = nde_delete_execute(doc_id, &err);
	unreserve_thread();
	cr_assert(!ok, "document-wide records must refuse deletion");
	cr_assert_not_null(err);
	g_free(err);

	golden_teardown(NULL, f);
}

/* ---------------- P4.1: barrier model ---------------- */

/* Apply one op headlessly to @f via the replay-flag worker (no capture). */
static void apply_direct(const op_descriptor *op, gpointer user, fits *f) {
	struct generic_img_args *args = calloc(1, sizeof(*args));
	args->fit = f;
	args->op = op;
	args->user = user;
	args->nde_replay = TRUE;
	args->mem_ratio = -1.0f;
	args->max_threads = 1;
	cr_assert_eq(GPOINTER_TO_INT(generic_image_worker(args)), 0);
	free_generic_img_args(args);
}

Test(nde_replay, barrier_checkpoint_enables_tail_editing) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* A1: asinh(15) */
	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 15.0f; u1->offset = 0.02f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);

	/* B: opaque freehand + output checkpoint.  P4.3's capture wiring stores
	 * the checkpoint automatically from the POST-op fits at a barrier
	 * capture — exercise that path (pass gfit) rather than a manual store. */
	gfit->fdata[3] = 0.987f;
	gint64 b_id = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "freehand", gfit);
	cr_assert(nde_checkpoint_output_exists(b_id),
	          "barrier capture with a post fits must store an output checkpoint");

	/* A2: asinh(20) */
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->offset = 0.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);
	gint64 a1_id = b_id - 1, a2_id = b_id + 1;

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable, "the full chain crosses an opaque record");
	cr_assert(chain->tail_replayable, "the tail restarts from the barrier checkpoint");
	cr_assert_eq(chain->records->len, 3);
	cr_assert_eq(chain->tail_start, 2);
	cr_assert_eq(chain->restart_ckpt_id, b_id);
	cr_assert(g_array_index(chain->member_flags, guint8, 1) & NDE_CHAIN_MEMBER_BARRIER);
	cr_assert_eq(chain->reasons->len, 0, "a checkpointed barrier is not a blocker");

	/* tail replay reproduces the current pixels bit-exactly */
	gchar *err = NULL;
	cr_assert(reserve_thread());
	fits *result = nde_chain_replay_tail(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "tail replay failed: %s", err ? err : "?");
	assert_pixels_bit_exact(result, gfit, "tail");
	clearfits(result); free(result);
	nde_chain_free(chain);

	/* frozen prefix refuses both amend and delete */
	cr_assert(reserve_thread());
	cr_assert(!nde_amend_execute(a1_id, "beta=99;offset=0;human=0;clip_mode=1", &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	cr_assert(!nde_delete_execute(a1_id, &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);

	/* tail amend works and replays from the checkpoint, not the baseline */
	gboolean ok = nde_amend_execute(a2_id, "beta=50;offset=0;human=0;clip_mode=1", &err);
	unreserve_thread();
	cr_assert(ok, "tail amend failed: %s", err ? err : "?");
	fits *expected = nde_checkpoint_output_get(b_id);
	cr_assert_not_null(expected);
	{
		asinh_params *eu = calloc(1, sizeof(*eu));
		eu->beta = 50.0f; eu->offset = 0.0f; eu->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, eu, expected);
	}
	assert_pixels_bit_exact(gfit, expected, "tail-amend");
	clearfits(expected); free(expected);

	/* delete of the (last) barrier: replay falls back to the baseline and
	 * the freehand poke disappears */
	cr_assert(reserve_thread());
	ok = nde_delete_execute(b_id, &err);
	unreserve_thread();
	cr_assert(ok, "barrier delete failed: %s", err ? err : "?");
	expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 15.0f; e1->offset = 0.02f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, expected);
		asinh_params *e2 = calloc(1, sizeof(*e2));
		e2->beta = 50.0f; e2->offset = 0.0f; e2->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e2, expected);
	}
	assert_pixels_bit_exact(gfit, expected, "post-barrier-delete");
	cr_assert(!nde_checkpoint_output_exists(b_id),
	          "the deleted barrier's checkpoint must be dropped");
	clearfits(expected); free(expected);

	golden_teardown(NULL, f);
}

Test(nde_replay, non_last_barrier_delete_refused_and_ckpt_less_barrier_blocks) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	struct mirror_args *m1 = calloc(1, sizeof(*m1));
	m1->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m1), 0);   /* A1 */

	gfit->fdata[5] = 0.911f;
	/* b1: barrier capture wires its own output checkpoint from gfit (P4.3). */
	gint64 b1 = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "b1", gfit);
	cr_assert(nde_checkpoint_output_exists(b1));
	fits b1_pixels = { 0 };
	copyfits(gfit, &b1_pixels, CP_DEEPCOPY | CP_ALLOC, -1);

	struct mirror_args *m2 = calloc(1, sizeof(*m2));
	m2->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m2), 0);   /* A2 */

	gfit->fdata[6] = 0.077f;
	/* b2: pass NULL — simulate a pre-phase-4 capture with no checkpoint. */
	gint64 b2 = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "b2", NULL);

	struct mirror_args *m3 = calloc(1, sizeof(*m3));
	m3->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m3), 0);   /* A3 */

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(!chain->tail_replayable, "a checkpoint-less last barrier blocks the tail");
	cr_assert_eq(chain->tail_start, 4);
	cr_assert(chain->reasons->len >= 1);
	nde_chain_free(chain);

	/* give b2 its checkpoint: the tail unfreezes, but b1 stays locked */
	nde_checkpoint_output_store(gfit /* wrong pixels for realism but fine for the model */, b2, -1);
	/* recompute a CORRECT b2 checkpoint: state right after b2 = mirrorx
	 * applied since... easier: the model checks only existence; pixel
	 * correctness for tail replay uses A3-on-b2ckpt below, so store the
	 * true post-b2 state: current is A3(b2); b2 = mirrorx(current) since
	 * A3 is an involution. */
	{
		fits post_b2 = { 0 };
		copyfits(gfit, &post_b2, CP_DEEPCOPY | CP_ALLOC, -1);
		mirrorx(&post_b2, FALSE);
		nde_checkpoint_output_store(&post_b2, b2, -1);
		clearfits(&post_b2);
	}
	chain = nde_chain_build(-1);
	cr_assert(chain->tail_replayable);
	cr_assert_eq(chain->restart_ckpt_id, b2);
	nde_chain_free(chain);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	/* non-last barrier: delete refused */
	cr_assert(!nde_delete_execute(b1, &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* last barrier: delete works; pixels = A3(A2(b1_ckpt)) = b1 pixels
	 * (mirror twice), the b2 poke gone */
	gboolean ok = nde_delete_execute(b2, &err);
	unreserve_thread();
	cr_assert(ok, "last-barrier delete failed: %s", err ? err : "?");
	assert_pixels_bit_exact(gfit, &b1_pixels, "post-b2-delete");

	clearfits(&b1_pixels);
	golden_teardown(NULL, f);
}

/* ---------------- P4.3: capture-wiring stores output checkpoints -------- */

/* A barrier capture (Tier B) via nde_capture_opaque with a POST fits stores
 * an output checkpoint keyed by the returned record id; without a fits it
 * does not (the pre-phase-4 no-restart-point behaviour). */
Test(nde_replay, capture_opaque_wires_output_checkpoint) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.3f);
	gfit = f;
	nde_checkpoint_baseline_ensure(f, -1);

	/* with a post fits → checkpoint stored */
	gint64 b = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1,
	                              "freehand", f);
	cr_assert(b > 0);
	cr_assert(nde_checkpoint_output_exists(b),
	          "a Tier-B capture with a post fits must store an output checkpoint");
	/* the stored pixels are the ones we passed */
	fits *ck = nde_checkpoint_output_get(b);
	cr_assert_not_null(ck);
	assert_pixels_bit_exact(ck, f, "opaque-checkpoint");
	clearfits(ck); free(ck);

	/* without a post fits → no checkpoint (pre-phase-4 barrier) */
	gint64 b2 = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1,
	                               "freehand-nockpt", NULL);
	cr_assert(b2 > 0);
	cr_assert(!nde_checkpoint_output_exists(b2),
	          "a Tier-B capture with NULL post must NOT store a checkpoint");

	golden_teardown(NULL, f);
}

/* A Tier-A capture is not a barrier: it stores no output checkpoint even when
 * a post fits is supplied.  A Tier-A capture with a mask active IS a barrier
 * and does store one. */
Test(nde_replay, capture_tier_a_checkpoint_only_when_mask_active) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.4f);
	gfit = f;
	nde_checkpoint_baseline_ensure(f, -1);

	/* Tier A, no mask → not a barrier → no checkpoint even with a post fits */
	struct mirror_args ma = { 0 };
	ma.x_axis = TRUE;
	gint64 a = nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", f);
	cr_assert(a > 0);
	cr_assert(!nde_checkpoint_output_exists(a),
	          "a Tier-A capture without a mask is not a barrier — no checkpoint");

	/* Tier A WITH a mask active → barrier → checkpoint stored.  The
	 * descriptor path reads gfit->mask/mask_active, so install a minimal
	 * mask (a non-NULL mask_t with mask_active is all the check needs). */
	cr_assert_null(gfit->mask);
	gfit->mask = calloc(1, sizeof(mask_t));
	gfit->mask->bitpix = FLOAT_IMG;
	gfit->mask->data = calloc((size_t)gfit->rx * gfit->ry, sizeof(float));
	gfit->mask_active = TRUE;

	gint64 am = nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m-masked", f);
	cr_assert(am > 0);
	cr_assert(nde_checkpoint_output_exists(am),
	          "a Tier-A capture with a mask active is a barrier — checkpoint stored");

	gfit->mask_active = FALSE;
	free(gfit->mask->data);
	free(gfit->mask);
	gfit->mask = NULL;

	golden_teardown(NULL, f);
}

/* The generic_image_worker capture block stores an output checkpoint for a
 * serializer-less (Tier B) op applied via the real path.  (This coverage used
 * op_desc_epf before phase 4.5 upgraded EPF to Tier A — Convention 1 file
 * operands — so it now drives op_desc_bg (stats.bg): a genuinely
 * serializer-less op with no external inputs that runs headlessly.  It does not
 * mutate pixels, so the stored checkpoint equals gfit after the op, which is
 * exactly the barrier-restart invariant under test.)  bg needs a heap user with
 * a destructor; we mirror its command.c-local struct layout (destructor-first,
 * plain free destructor). */
struct t_bg_data {
	void (*destructor)(void *);
	rectangle selection;
	gboolean has_selection;
	double bg_values[3];
	guint16 us_bg_values[3];
	int nchannels;
};
Test(nde_replay, worker_block_stores_checkpoint_for_serializerless_op) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	const op_descriptor *bg = op_descriptor_by_id("stats.bg");
	cr_assert_not_null(bg);
	cr_assert_null(bg->serialize, "stats.bg must stay serializer-less for this test");
	struct t_bg_data *u = calloc(1, sizeof(*u));
	cr_assert_not_null(u);
	u->destructor = free;
	cr_assert_eq(apply_op_real(bg, u), 0);

	/* one live record, Tier B (no serializer) */
	cr_assert_eq(nde_history_live_count(), 1);
	gint64 rid;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		nde_record *rec = g_ptr_array_index(snap, 0);
		cr_assert_eq(rec->tier, NDE_TIER_B, "stats.bg is serializer-less → Tier B");
		rid = rec->record_id;
		g_ptr_array_unref(snap);
	}
	/* the worker capture block stored the post-op pixels (== gfit now) */
	cr_assert(nde_checkpoint_output_exists(rid),
	          "the worker block must checkpoint a serializer-less barrier op");
	fits *ck = nde_checkpoint_output_get(rid);
	cr_assert_not_null(ck);
	assert_pixels_bit_exact(ck, gfit, "worker-bg-checkpoint");
	clearfits(ck); free(ck);

	golden_teardown(NULL, f);
}

/* ---------------- phase 4.5: file-operand replay_pre (Convention 1) --------
 *
 * imoper (arith.imoper) pins its operand file: the hook re-reads it from the
 * recorded path and replay_pre verifies existence/size/sha before the apply.
 * These tests write a real operand FITS to a temp file, drive the op through
 * the real capture path (apply_op_real), and check that (a) the chain replays
 * bit-exact against the temp operand and (b) tampering one byte makes both
 * replay_pre and a full chain replay fail with the surfaced error.
 *
 * imoper_data is file-local to command.c; we mirror its layout (the on-disk
 * keys are the contract) so the real serializer/hook read our heap struct. */
struct t_imoper_data { void (*destructor)(void *); int oper; char *filename; gboolean ftf; };

/* Mirror free_imoper_data (command.c, file-local): free filename then struct. */
static void t_imoper_destroy(void *p) {
	struct t_imoper_data *d = p;
	if (!d) return;
	free(d->filename);
	free(d);
}

/* Save @op_fit to a fresh temp .fit and return its heap path (caller g_free +
 * g_remove).  Uses a temp dir + a ".fit" name so set_right_extension keeps it
 * verbatim. */
static gchar *save_temp_operand_fits(fits *op_fit, gchar **dir_out) {
	GError *e = NULL;
	gchar *dir = g_dir_make_tmp("nde_replay_XXXXXX", &e);
	cr_assert_not_null(dir, "temp dir: %s", e ? e->message : "?");
	gchar *path = g_build_filename(dir, "operand.fit", NULL);
	cr_assert_eq(savefits(path, op_fit), 0, "savefits failed for %s", path);
	if (dir_out) *dir_out = dir; else g_free(dir);
	return path;
}

Test(nde_replay, golden_imoper_file_operand) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* operand: a constant 0.1 image of the same geometry */
	fits *op_fit = flis_test_make_mono_fits(16, 12, 0.1f);
	gchar *dir = NULL;
	gchar *path = save_temp_operand_fits(op_fit, &dir);
	clearfits(op_fit); free(op_fit);

	struct t_imoper_data *u = calloc(1, sizeof(*u));
	cr_assert_not_null(u);
	/* mirror process_imoper: OPER_ADD, strdup filename, force_to_float TRUE.
	 * The worker frees this heap struct via its destructor after the op. */
	u->destructor = t_imoper_destroy;
	u->oper = 0;                            /* OPER_ADD */
	u->filename = strdup(path);
	u->ftf = TRUE;

	cr_assert_eq(apply_op_real(&op_desc_imoper, u), 0);

	/* the captured record is Tier A and pins the operand path + hash */
	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert_eq(rec->tier, NDE_TIER_A, "imoper with operand keys is Tier A");
	cr_assert(strstr(rec->params, "operand_path=") != NULL);
	cr_assert(strstr(rec->params, "operand_sha256=") != NULL);
	g_ptr_array_unref(snap);

	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "imoper-file-operand");
	golden_teardown(result, f);

	g_remove(path); g_free(path);
	g_rmdir(dir); g_free(dir);
}

Test(nde_replay, imoper_replay_pre_verifies_and_tamper_fails) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	fits *op_fit = flis_test_make_mono_fits(16, 12, 0.2f);
	gchar *dir = NULL;
	gchar *path = save_temp_operand_fits(op_fit, &dir);
	clearfits(op_fit); free(op_fit);

	struct t_imoper_data *u = calloc(1, sizeof(*u));
	u->destructor = t_imoper_destroy;
	u->oper = 0;
	u->filename = strdup(path);
	u->ftf = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_imoper, u), 0);

	/* build the chain BEFORE tampering — it captures the record's params */
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable);

	/* matched file: replay_pre must accept (a clean replay succeeds) */
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *ok_res = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(ok_res, "clean replay failed: %s", err ? err : "?");
	clearfits(ok_res); free(ok_res); g_free(err); err = NULL;

	/* tamper ONE byte in the operand file: replay_pre must now fail, and the
	 * whole chain replay must fail with the surfaced error. */
	gchar *bytes = NULL; gsize len = 0;
	cr_assert(g_file_get_contents(path, &bytes, &len, NULL));
	cr_assert(len > 3000, "fits should be larger than its header");
	bytes[len - 1] ^= 0xFF;     /* flip the last data byte */
	cr_assert(g_file_set_contents(path, bytes, len, NULL));
	g_free(bytes);

	cr_assert(reserve_thread());
	fits *bad_res = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_null(bad_res, "replay must fail after the operand was tampered");
	cr_assert_not_null(err);
	g_free(err);
	nde_chain_free(chain);

	golden_teardown(NULL, f);
	g_remove(path); g_free(path);
	g_rmdir(dir); g_free(dir);
}

Test(nde_replay, truncation_drops_output_checkpoints) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.5f);
	gfit = f;
	nde_checkpoint_baseline_ensure(f, -1);
	/* barrier capture wires the checkpoint from the post fits (P4.3). */
	gint64 b = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "b", f);
	cr_assert(nde_checkpoint_output_exists(b));

	/* undo the record, then append: the dead tail truncation must drop
	 * the checkpoint with the record */
	nde_history_on_undo(b);
	nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "b2", NULL);
	cr_assert(!nde_checkpoint_output_exists(b),
	          "truncated records must lose their output checkpoints");

	golden_teardown(NULL, f);
}

/* ---------------- phase 4.5 Convention 2: synthstar golden replay ----------
 *
 * synthstar is detection-free when com.stars is populated: install a
 * manufactured star list, apply synthstar via the real capture path (which
 * stashes the effective star list into the record), then replay from the
 * baseline into a scratch fits.  replay_pre reinstalls the stashed stars, so
 * the hook renders the identical PSFs → bit-exact. */
static psf_star **make_replay_test_stars(int n, int rx, int ry) {
	psf_star **stars = malloc((size_t)(n + 1) * sizeof(psf_star *));
	for (int i = 0; i < n; i++) {
		psf_star *s = new_psf_star();
		s->xpos = 4.0 + 3.0 * i;
		s->ypos = 4.0 + 2.0 * i;
		if (s->xpos > rx - 2) s->xpos = rx - 2;
		if (s->ypos > ry - 2) s->ypos = ry - 2;
		s->A = 0.5 + 0.1 * i;
		s->fwhmx = 2.5;
		s->fwhmy = 2.5;
		s->beta = 4.0;
		s->has_saturated = FALSE;
		s->profile = PSF_GAUSSIAN;
		stars[i] = s;
	}
	stars[n] = NULL;
	return stars;
}

Test(nde_replay, golden_synthstar_from_stashed_stars) {
	fits *f = flis_test_make_mono_fits(24, 20, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* manufactured com.stars so synthstar skips detection entirely */
	replace_com_stars(make_replay_test_stars(3, f->rx, f->ry));

	struct synthstar_data *u = new_synthstar_data();
	cr_assert_eq(apply_op_real(&op_desc_synthstar, u), 0);

	/* the captured record is Tier A and carries the effective star list */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 1);
	nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert_eq(rec->tier, NDE_TIER_A, "synthstar with stars is Tier A");
	cr_assert(strstr(rec->params, "stars=") != NULL, "params: %s", rec->params);
	g_ptr_array_unref(snap);

	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "synthstar-stashed-stars");

	clear_stars_list(FALSE);
	golden_teardown(result, f);
}
