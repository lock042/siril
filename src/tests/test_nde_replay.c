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
#include "core/nde_snapstore.h"
#include "core/nde_replay.h"
#include "core/fits_region.h"
#include "core/masks.h"
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
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m1", NULL, FALSE) > 0);
	cr_assert(nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m2", NULL, FALSE) > 0);

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
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", NULL, FALSE);
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
	 * destroy_fn = free_curve_params, which owns and frees the per-channel
	 * points lists.  Mirrors src/gui-gtk4/curves.c
	 * build_curve_params_from_gui(). */
	struct curve_params *u = new_curve_params();
	cr_assert_not_null(u);
	GList *pts = NULL;
	point *p0 = g_new(point, 1); p0->x = 0.0; p0->y = 0.0; pts = g_list_append(pts, p0);
	point *p1 = g_new(point, 1); p1->x = 0.4; p1->y = 0.6; pts = g_list_append(pts, p1);
	point *p2 = g_new(point, 1); p2->x = 1.0; p2->y = 1.0; pts = g_list_append(pts, p2);
	u->channels[CHAN_RGB_K].points = pts;
	/* a range-masked L* curve on top, so the replay has to reproduce the
	 * mask parameters as well as the points */
	GList *lpts = NULL;
	point *l0 = g_new(point, 1); l0->x = 0.0; l0->y = 0.0; lpts = g_list_append(lpts, l0);
	point *l1 = g_new(point, 1); l1->x = 0.5; l1->y = 0.35; lpts = g_list_append(lpts, l1);
	point *l2 = g_new(point, 1); l2->x = 1.0; l2->y = 1.0; lpts = g_list_append(lpts, l2);
	u->channels[CHAN_L].points = lpts;
	u->channels[CHAN_L].range_enabled = TRUE;
	u->channels[CHAN_L].lum_min = 0.2f;
	u->channels[CHAN_L].lum_max = 0.8f;
	u->channels[CHAN_L].feather = 0.3f;
	u->algorithm = LINEAR;
	u->fit = gfit;
	u->verbose = FALSE;
	u->for_preview = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_curves, u), 0);
	fits *result = replay_current_chain(1);
	assert_pixels_bit_exact(result, gfit, "curves");
	golden_teardown(result, f);
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
	fake->mask_fd = -1;   /* snap NULL via g_new0 — no pixels */
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

	/* Deleting a property record takes the compositing fast path (graph step
	 * 1): it commits the log and re-folds the layer, without building a chain
	 * or replaying a single pixel.  It used to be refused outright. */
	gint64 prop_id;
	{
		GPtrArray *snap = nde_history_snapshot(NULL);
		prop_id = ((nde_record *)g_ptr_array_index(snap, 1))->record_id;
		g_ptr_array_unref(snap);
	}
	cr_assert(reserve_thread());
	gchar *err = NULL;
	cr_assert(nde_delete_execute(prop_id, &err), "%s", err ? err : "");
	unreserve_thread();
	cr_assert_eq(nde_history_live_count(), 1,
	             "the property record is gone; the pixel op survives");

	/* The pixel chain is untouched by the property record's removal. */
	chain = nde_chain_build(-1);
	cr_assert(chain->replayable);
	cr_assert_eq(chain->records->len, 1);
	nde_chain_free(chain);

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
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", NULL, FALSE);
	nde_capture_opaque("icc.convert", NDE_SCOPE_DOCUMENT, -1, "Converted profile", NULL);
	nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m2", NULL, FALSE);

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
 * a post fits is supplied.  missedasinh.png: a masked dialog capture used to
 * be a barrier outright — the flag was set from gfit but neither WHICH mask
 * nor its pixels were recorded, so a stretch applied from a dialog froze its
 * whole chain.  It now pins the mask and keeps its state, like the worker's
 * own capture, so it stays replayable; and a NON-mask-aware dialog (curves
 * before it learned to blend, median, background extraction) records an
 * unmasked step even with a mask sitting active, because its op never read
 * it. */
Test(nde_replay, masked_dialog_capture_pins_the_mask_instead_of_freezing) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.4f);
	gfit = f;
	nde_checkpoint_baseline_ensure(f, -1);

	/* Tier A, no mask → not a barrier → no checkpoint even with a post fits */
	struct mirror_args ma = { 0 };
	ma.x_axis = TRUE;
	gint64 a = nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m", f, FALSE);
	cr_assert(a > 0);
	cr_assert(!nde_checkpoint_output_exists(a),
	          "a Tier-A capture without a mask is not a barrier — no checkpoint");

	/* Tier A WITH a mask active, caller mask-aware: pinned + kept. */
	cr_assert_null(gfit->mask);
	gfit->mask = calloc(1, sizeof(mask_t));
	gfit->mask->bitpix = BYTE_IMG;
	gfit->mask->data = calloc((size_t)gfit->rx * gfit->ry, 1);
	gfit->mask_active = TRUE;

	gint64 am = nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m-masked", f, TRUE);
	cr_assert(am > 0);
	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *rec = g_ptr_array_index(snap, snap->len - 1);
	cr_assert(rec->mask_active);
	const nde_input_pin *pin = nde_record_input(rec, "mask");
	cr_assert_not_null(pin, "the capture must say WHICH mask");
	cr_assert_eq(pin->src_item_id, NDE_ITEM_PLAIN_MASK);
	cr_assert(nde_mask_pin_resolvable(rec), "and keep its pixels");
	cr_assert(!nde_checkpoint_output_exists(am),
	          "with the mask kept there is nothing to freeze — not a barrier");

	/* Same mask still active, but the caller's op ignored it. */
	gint64 nm = nde_capture_from_descriptor(&op_desc_mirrorx, &ma, "m-unaware", f, FALSE);
	cr_assert(nm > 0);
	g_ptr_array_unref(snap);
	snap = nde_history_snapshot(NULL);
	rec = g_ptr_array_index(snap, snap->len - 1);
	cr_assert(!rec->mask_active,
	          "an op that never read the mask must not claim it was masked");
	cr_assert_null(nde_record_input(rec, "mask"));
	g_ptr_array_unref(snap);

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

/* Stamp a round Gaussian star into a mono fits at (cx,cy). */
static void stamp_star(fits *f, float cx, float cy, float amp, float fwhm) {
	float sigma = fwhm / 2.3548f;
	float inv2s2 = 1.0f / (2.0f * sigma * sigma);
	int r = (int)(4.0f * sigma) + 1;
	for (int y = (int)cy - r; y <= (int)cy + r; y++) {
		for (int x = (int)cx - r; x <= (int)cx + r; x++) {
			if (x < 0 || y < 0 || x >= (int)f->rx || y >= (int)f->ry)
				continue;
			float dx = x - cx, dy = y - cy;
			float v = amp * expf(-(dx * dx + dy * dy) * inv2s2);
			float *p = &f->fdata[x + (size_t)y * f->rx];
			if (*p < v) *p = v;
		}
	}
}

/* DELEGATED provenance (star-detection refinement): unclip ALWAYS auto-detects
 * (it never consumes com.stars), so its record must be Tier A carrying the
 * detection PARAMETERS (stars_auto=1 + sf_*), NOT a pinned star list — and it
 * must replay by re-detecting.  com.stars is left empty here so capture takes
 * the auto-detect path; a detectable star field makes detection deterministic. */
Test(nde_replay, unclip_delegated_records_conf_and_replays) {
	fits *f = flis_test_make_mono_fits(128, 128, 0.02f);
	stamp_star(f, 32.f, 40.f, 0.9f, 3.2f);
	stamp_star(f, 80.f, 64.f, 0.8f, 3.0f);
	stamp_star(f, 50.f, 96.f, 0.85f, 3.4f);
	stamp_star(f, 100.f, 30.f, 0.75f, 3.1f);
	gfit = f;

	clear_stars_list(FALSE);   /* force the auto-detect (DELEGATED) path */

	struct synthstar_data *u = new_synthstar_data();
	cr_assert_eq(apply_op_real(&op_desc_unclip, u), 0);

	/* the record is Tier A, delegated: conf recorded, no pinned list */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 1);
	nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert_eq(rec->tier, NDE_TIER_A, "unclip is now Tier A (delegated conf)");
	cr_assert(strstr(rec->params, "stars_auto=1") != NULL, "params: %s", rec->params);
	cr_assert(strstr(rec->params, "sf_sigma=") != NULL, "params: %s", rec->params);
	cr_assert_null(strstr(rec->params, "stars="),
	               "delegated unclip must not pin a star list: %s", rec->params);
	g_ptr_array_unref(snap);

	/* replay re-detects with the recorded conf and rebuilds a valid image */
	fits *result = replay_current_chain(1);
	cr_assert_not_null(result);
	cr_assert_eq(result->rx, gfit->rx);
	cr_assert_eq(result->ry, gfit->ry);

	clear_stars_list(FALSE);
	golden_teardown(result, f);
}

/* ---------------- C3: cached restarts + deposits + invalidation -------- */

Test(nde_replay, second_amend_restarts_from_cached_deposit) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);
	gint64 id1 = 1, id2 = 2;

	/* first amend of record 2: no cache yet — replays from the baseline
	 * and deposits POST(1), POST(2) as it goes */
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(id2, "beta=30;offset=0;human=0;clip_mode=1", &err),
	          "first amend failed: %s", err ? err : "?");
	unreserve_thread();
	cr_assert(nde_snapstore_has(-1, id1, TRUE),
	          "the replay must deposit POST(record 1)");

	/* second amend of record 2: must restart from the POST(1) deposit —
	 * exactly one record replayed, and a cache hit in the stats */
	nde_snapstore_stats_reset();
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(id2, "beta=40;offset=0;human=0;clip_mode=1", &err),
	          "second amend failed: %s", err ? err : "?");
	unreserve_thread();
	nde_snapstore_stats_t st;
	nde_snapstore_stats(&st);
	cr_assert(st.hits >= 1, "the second amend must hit the cached restart");
	cr_assert_eq(st.deposits, 1,
	             "restarting at the edit means exactly ONE record replayed (one deposit)");

	/* cache-restart result must equal a from-baseline computation */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 10.0f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, expected);
		asinh_params *e2 = calloc(1, sizeof(*e2));
		e2->beta = 40.0f; e2->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e2, expected);
	}
	assert_pixels_bit_exact(gfit, expected, "cached-restart-amend");
	clearfits(expected); free(expected);

	golden_teardown(NULL, f);
}

Test(nde_replay, upstream_amend_invalidates_stale_deposits) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);

	gchar *err = NULL;
	/* seed the cache: amend record 2 (deposits POST(1) w/ beta 10 state) */
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(2, "beta=30;offset=0;human=0;clip_mode=1", &err));
	/* amend record 1: the old POST(1)/POST(2) are stale for the new chain
	 * and must be replaced, not reused */
	cr_assert(nde_amend_execute(1, "beta=15;offset=0;human=0;clip_mode=1", &err));
	/* now amend record 2 again: its cached restart POST(1) must reflect
	 * beta 15, or the result diverges */
	cr_assert(nde_amend_execute(2, "beta=40;offset=0;human=0;clip_mode=1", &err));
	unreserve_thread();

	fits *expected = nde_checkpoint_baseline_get(-1);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 15.0f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, expected);
		asinh_params *e2 = calloc(1, sizeof(*e2));
		e2->beta = 40.0f; e2->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e2, expected);
	}
	assert_pixels_bit_exact(gfit, expected, "post-invalidation-amend");
	clearfits(expected); free(expected);

	golden_teardown(NULL, f);
}

Test(nde_replay, truncation_evicts_pool_deposits) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);

	/* seed pool deposits via a verify replay */
	nde_chain *chain = nde_chain_build(-1);
	gchar *err = NULL;
	cr_assert(reserve_thread());
	fits *r = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(r);
	clearfits(r); free(r);
	nde_chain_free(chain);
	cr_assert(nde_snapstore_has(-1, 2, TRUE));

	/* undo record 2 then append: truncation must evict its deposit */
	nde_history_on_undo(2);
	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);
	cr_assert(!nde_snapstore_has(-1, 2, TRUE),
	          "a truncated record's pool deposit must be evicted");

	golden_teardown(NULL, f);
}

/* ---------------- C3.5: reordering ---------------- */

/* Build the standard reorder fixture: asinh(15) then mtf — deliberately
 * non-commuting so order changes the pixels. */
static void reorder_fixture_apply(void) {
	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 15.0f; u1->offset = 0.02f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	struct mtf_data *u2 = create_mtf_data();
	u2->linked = TRUE;
	u2->params.midtones = 0.35f;
	u2->params.shadows = 0.02f;
	u2->params.highlights = 0.98f;
	u2->params.do_red = u2->params.do_green = u2->params.do_blue = TRUE;
	u2->auto_display_compensation = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_mtf, u2), 0);
}

static void reorder_expected_mtf_then_asinh(fits *dst, float asinh_beta) {
	struct mtf_data *e2 = create_mtf_data();
	e2->linked = TRUE;
	e2->params.midtones = 0.35f;
	e2->params.shadows = 0.02f;
	e2->params.highlights = 0.98f;
	e2->params.do_red = e2->params.do_green = e2->params.do_blue = TRUE;
	e2->auto_display_compensation = FALSE;
	apply_direct(&op_desc_mtf, e2, dst);
	asinh_params *e1 = calloc(1, sizeof(*e1));
	e1->beta = asinh_beta; e1->offset = 0.02f; e1->clip_mode = RESCALE;
	apply_direct(&op_desc_asinh, e1, dst);
}

Test(nde_replay, reorder_replays_in_new_order) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	reorder_fixture_apply();   /* records: 1=asinh, 2=mtf */

	/* no-op move: after record 1 == current position */
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_reorder_execute(2, 1, TRUE, &err), "no-op move failed: %s", err ? err : "?");
	/* the real move: mtf before asinh */
	gboolean ok = nde_reorder_execute(2, 1, FALSE, &err);
	unreserve_thread();
	cr_assert(ok, "reorder failed: %s", err ? err : "?");

	/* log order is now [2, 1] with ids and timestamps preserved */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 0))->record_id, 2);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 1))->record_id, 1);
	g_ptr_array_unref(snap);

	/* pixels equal mtf-then-asinh applied from the baseline */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	reorder_expected_mtf_then_asinh(expected, 15.0f);
	assert_pixels_bit_exact(gfit, expected, "reorder");
	clearfits(expected); free(expected);

	/* the undo history was cleared (no meta-undo) */
	cr_assert_null(com.undo_stack);

	/* post-reorder amend with NON-monotonic member ids: amend record 1
	 * (now the LAST member) — the positional resolver must restart from
	 * the state after record 2 (first member), not misread ids */
	cr_assert(reserve_thread());
	ok = nde_amend_execute(1, "beta=40;offset=0.0199999996;human=0;clip_mode=1", &err);
	unreserve_thread();
	cr_assert(ok, "post-reorder amend failed: %s", err ? err : "?");
	expected = nde_checkpoint_baseline_get(-1);
	reorder_expected_mtf_then_asinh(expected, 40.0f);
	assert_pixels_bit_exact(gfit, expected, "post-reorder-amend");
	clearfits(expected); free(expected);

	golden_teardown(NULL, f);
}

Test(nde_replay, reorder_refusals) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* A1(asinh), B(opaque w/ checkpoint), A2(mirrorx): moving A2 before A1
	 * crosses the barrier — refused */
	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 15.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	gfit->fdata[3] = 0.99f;
	gint64 b_id = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "freehand", gfit);
	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(!nde_reorder_execute(3, 1, FALSE, &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* the barrier itself cannot be reordered */
	cr_assert(!nde_reorder_execute(b_id, 3, TRUE, &err));
	cr_assert(strstr(err, "cannot be reordered") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* structural anchor from another scope refuses via membership */
	nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER, -1,
	                       g_strdup("opacity=0.5"), "opacity");
	cr_assert(!nde_reorder_execute(3, 4, TRUE, &err));
	cr_assert(strstr(err, "not part of the same editable history") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	unreserve_thread();

	golden_teardown(NULL, f);
}

/* ---------------- C4: amend preview ---------------- */

Test(nde_replay, amend_preview_installs_pre_state_and_restores_bit_exact) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);

	fits *expected_pre_k = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected_pre_k);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 10.0f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, expected_pre_k);
	}
	fits expected_full = { 0 };
	copyfits(gfit, &expected_full, CP_DEEPCOPY | CP_ALLOC, -1);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_preview_begin_execute(2, &err),
	          "begin failed: %s", err ? err : "?");

	/* the display now shows the chain up to record 1 only */
	assert_pixels_bit_exact(gfit, expected_pre_k, "amend-preview-install");
	cr_assert(nde_amend_preview_active());
	cr_assert_eq(nde_amend_preview_record_id(), 2);
	cr_assert_str_eq(nde_amend_preview_op_id(), "stretch.asinh");
	cr_assert(strstr(nde_amend_preview_params(), "beta=20") != NULL,
	          "pre-fill params: %s", nde_amend_preview_params());

	/* every history edit is blocked while the preview is installed */
	cr_assert(!nde_amend_execute(1, "beta=15;offset=0;human=0;clip_mode=1", &err));
	cr_assert(strstr(err, "being edited") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	cr_assert(!nde_delete_execute(1, &err));
	g_clear_pointer(&err, g_free);
	cr_assert(!nde_reorder_execute(2, 1, FALSE, &err));
	g_clear_pointer(&err, g_free);

	/* cancel: the true pixels come back bit-exactly, metadata included */
	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err),
	          "cancel failed: %s", err ? err : "?");
	unreserve_thread();
	assert_pixels_bit_exact(gfit, &expected_full, "amend-preview-restore");
	cr_assert(!nde_amend_preview_active());

	clearfits(expected_pre_k); free(expected_pre_k);
	clearfits(&expected_full);
	golden_teardown(NULL, f);
}

Test(nde_replay, amend_preview_apply_uses_deposited_restart) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_preview_begin_execute(2, &err),
	          "begin failed: %s", err ? err : "?");
	/* synthesizing pre-K from the baseline deposited POST(1) */
	cr_assert(nde_snapstore_has(-1, 1, TRUE),
	          "pre-K synthesis must deposit intermediate states");

	/* apply: restore-first, then the amend — whose tail replay must
	 * restart from the POST(1) this preview session deposited (exactly
	 * one record replayed = one deposit) */
	nde_snapstore_stats_reset();
	cr_assert(nde_amend_preview_end_execute(TRUE,
			"beta=40;offset=0;human=0;clip_mode=1", &err),
	          "apply failed: %s", err ? err : "?");
	unreserve_thread();
	cr_assert(!nde_amend_preview_active());
	nde_snapstore_stats_t st;
	nde_snapstore_stats(&st);
	cr_assert(st.hits >= 1, "the amend must hit the preview's cached restart");
	cr_assert_eq(st.deposits, 1,
	             "restarting at the edit means exactly ONE record replayed");

	/* result equals a from-baseline computation with the new params */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 10.0f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, expected);
		asinh_params *e2 = calloc(1, sizeof(*e2));
		e2->beta = 40.0f; e2->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e2, expected);
	}
	assert_pixels_bit_exact(gfit, expected, "amend-preview-apply");
	clearfits(expected); free(expected);

	/* the log carries the new params */
	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *rec = g_ptr_array_index(snap, 1);
	cr_assert_eq(rec->record_id, 2);
	cr_assert(strstr(rec->params, "beta=40") != NULL, "params: %s", rec->params);
	g_ptr_array_unref(snap);

	golden_teardown(NULL, f);
}

/* ---------------- phase 9: region-scoped tail replay ---------------- */

/* The worker's region-preview algorithm, spelled out here so the claim can be
 * checked against the only thing it means: crop the pre-op image, grown by
 * the edited op's halo PLUS the tail's; run the edited op on the crop; replay
 * the tail over it; take the requested rectangle back out.  That rectangle
 * must equal a full-image computation of the whole amended chain, cropped to
 * the same place. */
Test(nde_replay, region_tail_replay_matches_a_full_recompute) {
	com.pref.nde_cache_mb = 256;   /* fixture memsets com — enable the pool */
	com.max_thread = 1;
	fits *f = flis_test_make_mono_fits(48, 40, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* asinh(10) — the record we amend — then a median with a real reach, then
	 * a second stretch, so the tail is neither trivially pixel-local nor a
	 * single step. */
	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	struct median_filter_data *u2 = calloc(1, sizeof(*u2));
	u2->ksize = 5; u2->iterations = 2; u2->amount = 1.0;
	cr_assert_eq(apply_op_real(&op_desc_median, u2), 0);
	asinh_params *u3 = calloc(1, sizeof(*u3));
	u3->beta = 20.0f; u3->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u3), 0);

	/* The truth: the whole chain, full-image, with record 1 amended. */
	fits *truth = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(truth);
	{
		asinh_params *e1 = calloc(1, sizeof(*e1));
		e1->beta = 25.0f; e1->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e1, truth);
		struct median_filter_data *e2 = calloc(1, sizeof(*e2));
		e2->ksize = 5; e2->iterations = 2; e2->amount = 1.0;
		apply_direct(&op_desc_median, e2, truth);
		asinh_params *e3 = calloc(1, sizeof(*e3));
		e3->beta = 20.0f; e3->clip_mode = RESCALE;
		apply_direct(&op_desc_asinh, e3, truth);
	}

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_preview_begin_execute(1, &err),
	          "begin failed: %s", err ? err : "?");

	int tail_halo = -1;
	nde_region_tail *plan = nde_region_tail_begin(&op_desc_asinh, &tail_halo);
	cr_assert_not_null(plan, "a median + asinh tail is region-replayable");
	cr_assert_eq(tail_halo, 4,
	             "the tail's halo is the median's radius x iterations (2x2); "
	             "the trailing asinh is pixel-local");

	/* Crop as the worker would: the requested rectangle grown by 0 (asinh)
	 * + tail_halo, clipped to the image. */
	const rectangle want = { 18, 14, 12, 10 };
	const rectangle grown = { want.x - tail_halo, want.y - tail_halo,
	                          want.w + 2 * tail_halo, want.h + 2 * tail_halo };
	fits region = { 0 };
	cr_assert_eq(crop_fits_region(gfit, &grown, &region), 0);

	asinh_params *k = calloc(1, sizeof(*k));
	k->beta = 25.0f; k->clip_mode = RESCALE;
	apply_direct(&op_desc_asinh, k, &region);

	/* Before replaying the tail, keep what the rectangle WOULD have shown
	 * without it — the discrimination control below. */
	fits without_tail = { 0 };
	copyfits(&region, &without_tail, CP_DEEPCOPY | CP_ALLOC, -1);

	cr_assert(nde_region_tail_apply(plan, &region, &grown, 1));
	nde_region_tail_free(plan);

	const rectangle inner = { tail_halo, tail_halo, want.w, want.h };
	fits got = { 0 }, expected = { 0 }, control = { 0 };
	cr_assert_eq(crop_fits_region(&region, &inner, &got), 0);
	cr_assert_eq(crop_fits_region(truth, &want, &expected), 0);
	cr_assert_eq(crop_fits_region(&without_tail, &inner, &control), 0);

	/* Not bit-exact: the median is OpenCV-backed and sums the same pixels in
	 * a different order at a different image width (op_descriptor.h).  1e-5 is
	 * two thirds of one ADU in 65535. */
	assert_pixels_close(&got, &expected, 1e-5f, "region tail replay");

	/* ...and the tail is what got it there.  Without this the test would pass
	 * just as happily on a tail that never ran, since a stretch of a stretch
	 * still looks like an image. */
	float ctl_dev = 0.f;
	for (size_t i = 0; i < (size_t)got.rx * got.ry; i++)
		ctl_dev = fmaxf(ctl_dev, fabsf(control.fdata[i] - expected.fdata[i]));
	cr_assert(ctl_dev > 1e-3f,
	          "the un-replayed rectangle must differ substantially (max %g)", ctl_dev);

	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err));
	unreserve_thread();

	clearfits(&got); clearfits(&expected); clearfits(&control);
	clearfits(&region); clearfits(&without_tail);
	clearfits(truth); free(truth);
	golden_teardown(NULL, f);
}

/* The regime the amend banner states, and the reason it gives when the answer
 * is no.  A geometry op in the tail is the clearest case: the rectangle it
 * would produce is not the rectangle that was asked for. */
Test(nde_replay, region_tail_reports_why_it_cannot_run) {
	com.pref.nde_cache_mb = 256;
	com.max_thread = 1;
	fits *f = flis_test_make_mono_fits(24, 20, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* Nothing installed: not a refusal, a different situation — and the
	 * banner must be able to tell them apart. */
	gchar *why = (gchar *)0x1;
	cr_assert(!nde_region_tail_available(NULL, &why));
	cr_assert_null(why, "no amend preview is not a reason to report");

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	struct mirror_args *u2 = calloc(1, sizeof(*u2));
	u2->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, u2), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_preview_begin_execute(1, &err),
	          "begin failed: %s", err ? err : "?");

	cr_assert(!nde_region_tail_available(NULL, &why));
	cr_assert_not_null(why, "a refusal must name the step responsible");
	cr_assert(strstr(why, "geometry") != NULL, "got: %s", why);
	g_clear_pointer(&why, g_free);
	cr_assert_null(nde_region_tail_begin(&op_desc_asinh, NULL));

	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err));
	unreserve_thread();
	golden_teardown(NULL, f);
}

/* Amending the LAST record: the tail is empty, so there is nothing hidden and
 * nothing to replay — which is a plan with a zero halo, not a refusal.  And a
 * preview of some OTHER op running while the amend preview is installed is
 * not this dialog's preview, so it gets no tail. */
Test(nde_replay, region_tail_empty_and_mismatched) {
	com.pref.nde_cache_mb = 256;
	com.max_thread = 1;
	fits *f = flis_test_make_mono_fits(24, 20, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 10.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	asinh_params *u2 = calloc(1, sizeof(*u2));
	u2->beta = 20.0f; u2->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u2), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_preview_begin_execute(2, &err),
	          "begin failed: %s", err ? err : "?");

	int halo = -1;
	cr_assert(nde_region_tail_available(&halo, NULL));
	cr_assert_eq(halo, 0, "an empty tail needs no context");

	nde_region_tail *plan = nde_region_tail_begin(&op_desc_asinh, NULL);
	cr_assert_not_null(plan);
	/* Applying an empty tail is a no-op, not an error. */
	fits region = { 0 };
	const rectangle r = { 4, 4, 8, 8 };
	cr_assert_eq(crop_fits_region(gfit, &r, &region), 0);
	fits before = { 0 };
	copyfits(&region, &before, CP_DEEPCOPY | CP_ALLOC, -1);
	cr_assert(nde_region_tail_apply(plan, &region, &r, 1));
	assert_pixels_bit_exact(&region, &before, "empty tail");
	nde_region_tail_free(plan);
	clearfits(&region); clearfits(&before);

	cr_assert_null(nde_region_tail_begin(&op_desc_median, NULL),
	               "a different op's preview is not the amend dialog's");

	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err));
	unreserve_thread();
	golden_teardown(NULL, f);
}

Test(nde_replay, amend_preview_refusals) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* asinh, opaque-with-checkpoint, mirrorx: record 1 is frozen */
	asinh_params *u1 = calloc(1, sizeof(*u1));
	u1->beta = 15.0f; u1->clip_mode = RESCALE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, u1), 0);
	gfit->fdata[3] = 0.99f;
	gint64 b_id = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1, "freehand", gfit);
	struct mirror_args *m = calloc(1, sizeof(*m));
	m->x_axis = TRUE;
	cr_assert_eq(apply_op_real(&op_desc_mirrorx, m), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	/* unknown id */
	cr_assert(!nde_amend_preview_begin_execute(99, &err));
	cr_assert(strstr(err, "no live record") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* opaque record */
	cr_assert(!nde_amend_preview_begin_execute(b_id, &err));
	cr_assert(strstr(err, "opaque") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* frozen prefix record */
	cr_assert(!nde_amend_preview_begin_execute(1, &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	/* apply-end without an active preview fails; cancel-end is tolerated */
	cr_assert(!nde_amend_preview_end_execute(TRUE, "beta=1", &err));
	g_clear_pointer(&err, g_free);
	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err));
	/* a second begin while one is installed refuses */
	cr_assert(nde_amend_preview_begin_execute(3, &err), "begin failed: %s", err ? err : "?");
	cr_assert(!nde_amend_preview_begin_execute(3, &err));
	cr_assert(strstr(err, "already being edited") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	cr_assert(nde_amend_preview_end_execute(FALSE, NULL, &err));
	unreserve_thread();

	golden_teardown(NULL, f);
}

/* ---------------- graph step 2: edit at / insert before K ---------------- */

/* Heap asinh params, the shape apply_op_real / apply_direct expect. */
static asinh_params *asinh_beta(float beta) {
	asinh_params *p = calloc(1, sizeof(*p));
	p->beta = beta;
	p->clip_mode = RESCALE;
	return p;
}

/* Record ids in log order. */
static void assert_log_order(const gint64 *want, guint n, const char *ctx) {
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap, "%s: no history", ctx);
	cr_assert_eq(snap->len, n, "%s: expected %u records, got %u", ctx, n, snap->len);
	for (guint i = 0; i < n; i++)
		cr_assert_eq(((nde_record *)g_ptr_array_index(snap, i))->record_id, want[i],
		             "%s: position %u", ctx, i);
	g_ptr_array_unref(snap);
}

Test(nde_replay, edit_at_inserts_a_step_and_replays_forward) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);

	/* what the image must become: the inserted step runs BETWEEN the two */
	fits *expected = nde_checkpoint_baseline_get(-1);
	cr_assert_not_null(expected);
	apply_direct(&op_desc_asinh, asinh_beta(10.f), expected);
	apply_direct(&op_desc_asinh, asinh_beta(30.f), expected);
	apply_direct(&op_desc_asinh, asinh_beta(20.f), expected);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_edit_at_begin_execute(2, &err), "begin failed: %s", err ? err : "?");
	cr_assert(nde_edit_at_active());
	cr_assert(nde_amend_preview_active(), "edit-at shares the single-instance lock");
	cr_assert_eq(nde_edit_at_record_id(), 2);

	/* every history edit is refused while the point is open */
	cr_assert(!nde_amend_execute(1, "beta=15;offset=0;human=0;clip_mode=1", &err));
	cr_assert(strstr(err, "insertion point") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	cr_assert(!nde_delete_execute(1, &err));
	g_clear_pointer(&err, g_free);
	cr_assert(!nde_reorder_execute(2, 1, FALSE, &err));
	g_clear_pointer(&err, g_free);
	unreserve_thread();

	/* the user runs an ordinary operation against the pre-anchor state */
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(30.f)), 0);
	gint64 want_open[] = { 1, 3, 2 };
	assert_log_order(want_open, 3, "while open");

	cr_assert(reserve_thread());
	cr_assert(nde_edit_at_end_execute(TRUE, &err), "finish failed: %s", err ? err : "?");
	unreserve_thread();
	cr_assert(!nde_edit_at_active());
	cr_assert(!nde_amend_preview_active());

	assert_log_order(want_open, 3, "after finish");
	assert_pixels_bit_exact(gfit, expected, "edit-at-apply");

	clearfits(expected); free(expected);
	golden_teardown(NULL, f);
}

Test(nde_replay, edit_at_cancel_discards_the_inserted_step) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);
	fits expected_full = { 0 };
	copyfits(gfit, &expected_full, CP_DEEPCOPY | CP_ALLOC, -1);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_edit_at_begin_execute(2, &err), "begin failed: %s", err ? err : "?");
	unreserve_thread();

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(30.f)), 0);
	cr_assert_eq(nde_history_live_count(), 3);

	cr_assert(reserve_thread());
	cr_assert(nde_edit_at_end_execute(FALSE, &err), "cancel failed: %s", err ? err : "?");
	unreserve_thread();

	gint64 want[] = { 1, 2 };
	assert_log_order(want, 2, "after cancel");
	assert_pixels_bit_exact(gfit, &expected_full, "edit-at-cancel-restore");

	clearfits(&expected_full);
	golden_teardown(NULL, f);
}

Test(nde_replay, edit_at_inserting_nothing_is_a_no_op) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);
	fits expected_full = { 0 };
	copyfits(gfit, &expected_full, CP_DEEPCOPY | CP_ALLOC, -1);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_edit_at_begin_execute(2, &err), "begin failed: %s", err ? err : "?");
	cr_assert(nde_edit_at_end_execute(TRUE, &err), "finish failed: %s", err ? err : "?");
	unreserve_thread();

	gint64 want[] = { 1, 2 };
	assert_log_order(want, 2, "after empty finish");
	assert_pixels_bit_exact(gfit, &expected_full, "edit-at-empty");

	clearfits(&expected_full);
	golden_teardown(NULL, f);
}

Test(nde_replay, edit_at_refusals) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	/* an opaque step with a checkpoint: everything before it is frozen */
	gint64 b_id = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1,
	                                 "freehand", gfit);
	cr_assert(b_id > 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);

	gchar *err = NULL;
	cr_assert(reserve_thread());

	/* unknown record */
	cr_assert(!nde_edit_at_begin_execute(99, &err));
	cr_assert_not_null(err);
	g_clear_pointer(&err, g_free);

	/* frozen: record 1 sits before the barrier */
	cr_assert(!nde_edit_at_begin_execute(1, &err));
	cr_assert(strstr(err, "locked") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);

	/* finishing when nothing is open */
	cr_assert(!nde_edit_at_end_execute(TRUE, &err));
	g_clear_pointer(&err, g_free);
	cr_assert(nde_edit_at_end_execute(FALSE, &err), "a stray cancel is tolerated");

	/* the opaque step itself is no anchor either: closing the point would
	 * have to re-run it over the inserted work, which is the one thing an
	 * opaque record cannot do */
	cr_assert(!nde_edit_at_begin_execute(b_id, &err));
	cr_assert(strstr(err, "opaque step") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);

	/* record 3 is in the tail: a valid anchor, and only one at a time */
	cr_assert(nde_edit_at_begin_execute(3, &err), "begin failed: %s", err ? err : "?");
	cr_assert(!nde_edit_at_begin_execute(3, &err));
	cr_assert(strstr(err, "already") != NULL, "got: %s", err);
	g_clear_pointer(&err, g_free);
	cr_assert(nde_edit_at_end_execute(FALSE, &err), "cancel failed: %s", err ? err : "?");
	unreserve_thread();

	golden_teardown(NULL, f);
}

/* ---------------- graph step 4: pinned mask inputs ---------------- */

/* Attach a non-uniform float mask to @f so a masked op differs visibly from
 * an unmasked one — a uniform mask would make the assertions vacuous. */
static void attach_gradient_mask(fits *f) {
	size_t n = (size_t)f->rx * f->ry;
	mask_t *m = calloc(1, sizeof(mask_t));
	m->bitpix = 32;
	m->data = malloc(n * sizeof(float));
	for (size_t i = 0; i < n; i++)
		((float *)m->data)[i] = (float)i / (float)n;
	f->mask = m;
	f->mask_active = TRUE;
}

/* Run an op through the real capture path with the mask honoured. */
static int apply_op_masked(const op_descriptor *op, gpointer user_heap) {
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->op = op;
	args->user = user_heap;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;
	args->mask_aware = TRUE;
	/* Match the replay driver's thread count: an op whose result depends on
	 * how the work is divided is only bit-comparable against a run divided
	 * the same way. */
	args->max_threads = com.max_thread;
	args->mem_ratio = -1.0f;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_image_worker(args));
	com.headless = prev;
	return rc;
}

/* A masked record used to freeze its whole chain: the flag said a mask was
 * active but nothing said which, so there was nothing to reproduce. */
Test(nde_replay, masked_record_pins_its_mask_and_stays_replayable) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	attach_gradient_mask(gfit);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	const nde_record *masked = g_ptr_array_index(snap, 1);
	cr_assert(masked->mask_active, "the capture must still note the mask");
	const nde_input_pin *pin = nde_record_input(masked, "mask");
	cr_assert_not_null(pin, "and must pin WHICH mask");
	cr_assert_eq(pin->src_item_id, NDE_ITEM_PLAIN_MASK, "plain image: the sentinel item");
	cr_assert(nde_mask_pin_resolvable(masked), "its pixels must have been kept");
	g_ptr_array_unref(snap);

	/* The chain is replayable end to end, and the masked record is not a
	 * barrier: nothing before it is frozen. */
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable, "masked record must not freeze the chain: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert_eq(chain->tail_start, 0, "no freeze cause anywhere");
	nde_chain_free(chain);
}

/* The point of restoring the mask is that the replay reproduces the ORIGINAL
 * pixels.  Replaying the same chain without the mask would not. */
Test(nde_replay, replay_restores_the_pinned_mask_bit_exactly) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	attach_gradient_mask(gfit);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	fits expected = { 0 };
	copyfits(gfit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	assert_pixels_bit_exact(result, &expected, "masked-replay");
	/* the mask does not leak out of the replay onto the result */
	cr_assert_null(result->mask, "the mask belongs to the record, not the chain");
	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	golden_teardown(NULL, f);
}

/* Losing the pinned pixels degrades honestly: the record becomes a barrier
 * again rather than replaying unmasked and claiming success. */
Test(nde_replay, a_masked_record_whose_mask_was_dropped_is_a_barrier) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	attach_gradient_mask(gfit);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *masked = g_ptr_array_index(snap, 1);
	const nde_input_pin *pin = nde_record_input(masked, "mask");
	cr_assert_not_null(pin);
	/* The mask had no records of its own here, so it was pinned at the
	 * item's baseline; drop it the way a storage limit eventually will. */
	cr_assert_eq(pin->src_record_id, 0);
	nde_checkpoint_drop(pin->src_item_id);
	cr_assert(!nde_mask_pin_resolvable(masked));
	g_ptr_array_unref(snap);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable, "an unresolvable mask must freeze the chain");
	cr_assert(chain->reasons->len > 0);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "mask") != NULL,
	          "and say why: %s", (char *)g_ptr_array_index(chain->reasons, 0));
	nde_chain_free(chain);
	golden_teardown(NULL, f);
}

/* An unmasked record next to a masked one must not inherit its mask. */
Test(nde_replay, an_unmasked_record_replays_without_a_mask) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	attach_gradient_mask(gfit);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(10.f)), 0);
	/* drop the mask, then run an ordinary op */
	free_mask(gfit->mask);
	gfit->mask = NULL;
	gfit->mask_active = FALSE;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	cr_assert_not_null(nde_record_input(g_ptr_array_index(snap, 0), "mask"));
	cr_assert_null(nde_record_input(g_ptr_array_index(snap, 1), "mask"),
	               "the second op ran unmasked and must pin nothing");
	g_ptr_array_unref(snap);

	fits expected = { 0 };
	copyfits(gfit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	assert_pixels_bit_exact(result, &expected, "mixed-mask-replay");
	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	golden_teardown(NULL, f);
}

/* ---------------- graph step 5: geometry on a layered document ---------- */

/* Geometry on a FLIS layer used to be a hard blocker: the op replays fine on
 * a scratch fits, but a layer's value is pixels AND its position, and the
 * position had nowhere to come from.  With one recorded against the baseline
 * the chain replays like any other. */
Test(nde_replay, flis_geometry_chain_replays_and_lands_the_layer) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *lay = flis_test_add_layer(
	    flis_test_make_mono_fits(64, 48, 0.f), "base");
	lay->position_x = 20;
	lay->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	fill_mono_gradient(gfit);
	gint item = lay->item_id;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	struct crop_args *ca = calloc(1, sizeof(*ca));
	ca->area = (rectangle){ .x = 8, .y = 6, .w = 32, .h = 24 };
	cr_assert_eq(apply_op_real(&op_desc_crop, ca), 0);

	/* the live op moved the layer to its offset plus the selection origin */
	cr_assert_eq(lay->position_x, 28);
	cr_assert_eq(lay->position_y, 36);
	/* and the baseline kept where it started, not where the crop left it */
	gint bx = 0, by = 0;
	cr_assert(nde_checkpoint_baseline_get_offset(item, &bx, &by));
	cr_assert_eq(bx, 20, "the baseline's position is the PRE-first-op one");
	cr_assert_eq(by, 30);

	nde_chain *chain = nde_chain_build(item);
	cr_assert(chain->replayable, "geometry must no longer block: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert(chain->has_geometry, "and the chain must know it carries a position");
	cr_assert_eq(chain->records->len, 2);

	fits expected = { 0 };
	copyfits(gfit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");
	assert_pixels_bit_exact(result, &expected, "flis-geometry-replay");
	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	nde_history_attach(NULL);
	gfit = NULL;
}

/* Amending the crop must move the layer to match the new pixels: the offset
 * is re-derived from the baseline position, not left where the old crop put
 * it. */
Test(nde_replay, amending_a_crop_moves_the_layer) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *lay = flis_test_add_layer(
	    flis_test_make_mono_fits(64, 48, 0.f), "base");
	lay->position_x = 20;
	lay->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	fill_mono_gradient(gfit);

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	struct crop_args *ca = calloc(1, sizeof(*ca));
	ca->area = (rectangle){ .x = 8, .y = 6, .w = 32, .h = 24 };
	cr_assert_eq(apply_op_real(&op_desc_crop, ca), 0);
	cr_assert_eq(lay->position_x, 28);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(2, "x=4;y=2;w=16;h=12", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();

	/* re-derived from the BASELINE position (20,30), not from (28,36) */
	cr_assert_eq(lay->position_x, 24, "layer must move to match the new crop");
	cr_assert_eq(lay->position_y, 32);
	cr_assert_eq(lay->fit->rx, 16u);
	cr_assert_eq(lay->fit->ry, 12u);

	nde_history_attach(NULL);
	gfit = NULL;
}

/* An amended geometry step on a CONSUMED input must move where that input
 * lands in the composite.  A chain with a geometry member re-derives its
 * position on replay — the same rule the live layer's own amend follows
 * (amending_a_crop_moves_the_layer above) — so the render must use the
 * replayed position, not the offset recorded at capture time.  Before the
 * fix the amended pixels were rendered at the stale recorded offset. */
Test(nde_replay, amending_a_crop_on_a_consumed_input_moves_it_in_the_composite) {
	com.pref.nde_cache_mb = 256;
	flis_test_add_layer(flis_test_make_mono_fits(64, 48, 0.f), "bottom");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(64, 48, 0.f), "top");
	uniq_set_active_layer(com.uniq, 1);
	gfit = flis_active_layer_fit();
	fill_mono_gradient(gfit);
	gint top_item = top->item_id;

	/* a reference copy of the gradient, indexed in canvas coordinates */
	fits *ref = flis_test_make_mono_fits(64, 48, 0.f);
	fill_mono_gradient(ref);

	struct crop_args *ca = calloc(1, sizeof(*ca));
	ca->area = (rectangle){ .x = 8, .y = 6, .w = 32, .h = 24 };
	cr_assert_eq(apply_op_real(&op_desc_crop, ca), 0);
	cr_assert_eq(top->position_x, 8);
	cr_assert_eq(top->position_y, 6);

	cr_assert_eq(flis_merge_down_layer(top), 0);
	cr_assert_eq(flis_layer_count(), 1);
	flis_layer_t *merged = (flis_layer_t *)com.uniq->layers->data;
	gfit = merged->fit;

	/* the crop record on the consumed input */
	gint64 crop_id = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (r->target_item_id == top_item && !g_strcmp0(r->op_id, "geometry.crop"))
			crop_id = r->record_id;
	}
	g_ptr_array_unref(snap);
	cr_assert_neq(crop_id, 0, "the consumed input's crop record must survive");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(crop_id, "x=4;y=2;w=32;h=24", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();

	/* The patch now sits at (4,2), so inside it the composite reproduces the
	 * gradient IN PLACE — content and position both follow the amended crop —
	 * and outside it the bottom layer's zeros show.  With the stale-offset
	 * bug the same content rendered at (8,6). */
	/* Positions are display coordinates; fdata rows run bottom-up, so the
	 * patch at position (4,2), 32x24 on the 64x48 canvas occupies fdata
	 * x[4..35], rows [48-2-24 .. 48-2) = [22..45].  Both the crop window and
	 * the placement flip the same way, so inside the patch the composite
	 * reproduces the gradient at IDENTICAL raw indices. */
	const fits *m = merged->fit;
	cr_assert_eq(m->rx, 64u);
	cr_assert_eq(m->ry, 48u);
	#define MPX(x, y) m->fdata[(x) + (size_t)(y) * m->rx]
	#define RPX(x, y) ref->fdata[(x) + (size_t)(y) * ref->rx]
	cr_assert_float_eq(MPX(4, 22), RPX(4, 22), 1e-6,
	                   "patch corner must land at the amended origin");
	cr_assert_float_eq(MPX(35, 45), RPX(35, 45), 1e-6,
	                   "opposite patch corner too");
	cr_assert_float_eq(MPX(2, 1), 0.f, 1e-6,
	                   "outside the patch the bottom layer shows");
	cr_assert_float_eq(MPX(38, 30), 0.f, 1e-6,
	                   "nothing may render at the stale pre-amend offset");
	#undef MPX
	#undef RPX

	clearfits(ref); free(ref);
	nde_history_attach(NULL);
	gfit = NULL;
}

/* No recorded starting position (an older file) degrades honestly rather than
 * replaying from a position it made up. */
Test(nde_replay, flis_geometry_without_a_start_position_is_a_blocker) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *lay = flis_test_add_layer(
	    flis_test_make_mono_fits(64, 48, 0.f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	fill_mono_gradient(gfit);
	gint item = lay->item_id;

	struct crop_args *ca = calloc(1, sizeof(*ca));
	ca->area = (rectangle){ .x = 8, .y = 6, .w = 32, .h = 24 };
	cr_assert_eq(apply_op_real(&op_desc_crop, ca), 0);
	cr_assert(nde_checkpoint_baseline_get_offset(item, NULL, NULL));

	/* Drop just the position, the way an older file arrives: pixels present,
	 * no FLIS_POSX/Y.  nde_checkpoint_drop would take the baseline too, so
	 * re-adopt the pixels afterwards. */
	fits *base = nde_checkpoint_baseline_get(item);
	cr_assert_not_null(base);
	nde_checkpoint_drop(item);
	nde_checkpoint_baseline_adopt(base, item);
	clearfits(base); free(base);
	cr_assert(nde_checkpoint_baseline_exists(item));
	cr_assert(!nde_checkpoint_baseline_get_offset(item, NULL, NULL));

	nde_chain *chain = nde_chain_build(item);
	cr_assert(!chain->replayable);
	cr_assert(!chain->has_geometry);
	cr_assert(chain->reasons->len > 0);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "position") != NULL,
	          "and say why: %s", (char *)g_ptr_array_index(chain->reasons, 0));
	nde_chain_free(chain);
	nde_history_attach(NULL);
	gfit = NULL;
}

/* ---------------- graph step 6a: mask chains replay ---------------- */

/* Run a mask op through the real capture path, so it lands in the log the
 * way a user's would. */
static int apply_mask_op(const op_descriptor *op, gpointer user_heap) {
	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit = gfit;
	args->op = op;
	args->user = user_heap;
	args->command = TRUE;
	args->mask_creation = (op->flags & OP_MASK_FROM_IMAGE) != 0;
	args->max_threads = 1;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rc = GPOINTER_TO_INT(generic_mask_worker(args));
	com.headless = prev;
	return rc;
}

static mask_from_channel_data *from_channel(int chan) {
	mask_from_channel_data *d = calloc(1, sizeof(*d));
	d->channel = chan;
	d->bitpix = 8;
	return d;
}

static mask_blur_data *blur_mask(float radius) {
	mask_blur_data *d = calloc(1, sizeof(*d));
	d->radius = radius;
	return d;
}

/* Copy the live processing mask so a later state can be compared against it. */
static guint8 *snapshot_mask(const fits *f, size_t *n_out) {
	cr_assert_not_null(f->mask);
	cr_assert_not_null(f->mask->data);
	size_t n = (size_t)f->rx * f->ry * (f->mask->bitpix / 8);
	guint8 *copy = g_malloc(n);
	memcpy(copy, f->mask->data, n);
	*n_out = n;
	return copy;
}

static void assert_mask_differs(const fits *f, const guint8 *before, size_t n,
                                const char *ctx) {
	cr_assert_not_null(f->mask);
	cr_assert_eq((size_t)f->rx * f->ry * (f->mask->bitpix / 8), n,
	             "%s: mask changed size", ctx);
	cr_assert_neq(memcmp(f->mask->data, before, n), 0,
	              "%s: the mask was not rebuilt", ctx);
}

/* A mask's chain is recognised by its ops, and it starts at the image its
 * first member read rather than at a baseline of its own. */
Test(nde_replay, mask_chain_is_recognised_and_replayable) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(1.f)), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	const nde_record *create = g_ptr_array_index(snap, 0);
	cr_assert_eq(create->target_item_id, NDE_ITEM_PLAIN_MASK);
	const nde_input_pin *img = nde_record_input(create, "image");
	cr_assert_not_null(img, "a mask built from the image must say so");
	cr_assert_eq(img->src_item_id, NDE_ITEM_IMAGE);
	cr_assert_null(nde_record_input(g_ptr_array_index(snap, 1), "image"),
	               "editing a mask reads the mask, not the image");
	g_ptr_array_unref(snap);

	nde_chain *chain = nde_chain_build(NDE_ITEM_PLAIN_MASK);
	cr_assert(chain->is_mask);
	cr_assert_eq(chain->records->len, 2);
	cr_assert(chain->replayable,
	          "a mask chain needs no baseline of its own: %s",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);
	golden_teardown(NULL, f);
}

/* The maintainer's 5.2 example: build a mask, then go back and change the
 * parameters of the operation that built it. */
Test(nde_replay, amending_a_mask_op_rebuilds_the_mask) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(1.f)), 0);
	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(2, "radius=6", &err), "amend failed: %s", err ? err : "?");
	unreserve_thread();
	assert_mask_differs(gfit, before, n, "amend-blur");

	/* and the log carries the new value */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 1))->params, "radius=6");
	g_ptr_array_unref(snap);
	g_free(before);
	golden_teardown(NULL, f);
}

/* Amending the CREATION op re-derives from the image, which is the edge the
 * pin exists for. */
Test(nde_replay, amending_the_creating_op_re_derives_from_the_image) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_rgb_fits(16, 12, 0.2f, 0.5f, 0.8f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	/* a different channel of the same image gives a different mask */
	cr_assert(nde_amend_execute(1, "channel=2;autostretch=0;invert=0;bitpix=8", &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();
	assert_mask_differs(gfit, before, n, "amend-from_channel");
	g_free(before);
	golden_teardown(NULL, f);
}

/* Deleting a step from a mask's chain rebuilds it without that step. */
Test(nde_replay, deleting_a_mask_step_rebuilds_the_mask) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	size_t n = 0;
	guint8 *unblurred = snapshot_mask(gfit, &n);
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(3.f)), 0);
	size_t n2 = 0;
	guint8 *blurred = snapshot_mask(gfit, &n2);
	cr_assert_eq(n, n2);
	cr_assert_neq(memcmp(unblurred, blurred, n), 0, "fixture: the blur did something");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_delete_execute(2, &err), "delete failed: %s", err ? err : "?");
	unreserve_thread();
	cr_assert_eq(memcmp(gfit->mask->data, unblurred, n), 0,
	             "removing the blur must give back the unblurred mask exactly");
	g_free(unblurred);
	g_free(blurred);
	golden_teardown(NULL, f);
}

/* A mask whose origin is not recorded cannot be rebuilt, and says so rather
 * than inventing a source. */
Test(nde_replay, a_mask_with_no_recorded_origin_is_a_blocker) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	/* A mask EDIT with no creation step before it: nothing says where the
	 * mask came from.  (Reached in practice by a file whose creation record
	 * was deleted, or one written before mask provenance existed.) */
	gfit->mask = calloc(1, sizeof(mask_t));
	gfit->mask->bitpix = 8;
	gfit->mask->data = calloc(1, (size_t)16 * 12);
	gfit->mask_active = TRUE;
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(1.f)), 0);

	nde_chain *chain = nde_chain_build(NDE_ITEM_PLAIN_MASK);
	cr_assert(chain->is_mask);
	cr_assert(!chain->replayable);
	cr_assert(chain->reasons->len > 0);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "built from") != NULL,
	          "and say why: %s", (char *)g_ptr_array_index(chain->reasons, 0));
	nde_chain_free(chain);
	golden_teardown(NULL, f);
}

/* ---------------- graph step 6b: reverse invalidation ---------------- */

/* Build: mask from channel 0, blur it by @radius, then a masked stretch.
 * Returns the resulting image pixels (caller g_free()s) and its size. */
static float *run_masked_sequence(float radius, size_t *n_out) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(radius)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	size_t n = (size_t)gfit->rx * gfit->ry;
	float *out = g_malloc(n * sizeof(float));
	memcpy(out, gfit->fdata, n * sizeof(float));
	*n_out = n;
	nde_history_attach(NULL);      /* also purges the checkpoint store */
	clearfits(f);
	free(f);
	gfit = NULL;
	return out;
}

/* The whole point of the edges: amending the mask must reach the image that
 * used it.  The result has to equal the same work done in that order from the
 * start — "it changed" would not prove it changed to the RIGHT thing. */
Test(nde_replay, amending_a_mask_recomputes_what_consumed_it) {
	com.pref.nde_cache_mb = 256;
	size_t n_native = 0;
	float *native = run_masked_sequence(6.f, &n_native);

	/* now the same thing built at radius 1 and amended to 6 */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(1.f)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	size_t n = (size_t)gfit->rx * gfit->ry;
	cr_assert_eq(n, n_native);
	float *before = g_malloc(n * sizeof(float));
	memcpy(before, gfit->fdata, n * sizeof(float));
	cr_assert_neq(memcmp(before, native, n * sizeof(float)), 0,
	              "fixture: radius 1 and radius 6 must differ");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(2, "radius=6", &err), "amend failed: %s", err ? err : "?");
	unreserve_thread();

	/* Converges on the native result rather than merely "changing": a wrong
	 * recomputation would land somewhere else entirely.  Not asserted
	 * bit-exact because replaying an op that ran under an 8-bit mask is not
	 * bit-exact TODAY, independently of any of this — a plain nde_chain_replay
	 * of one deviates by the same single float ulp (see the regression bound
	 * below).  The bound here is two ulps at 0.5, i.e. still exact to within
	 * the existing gap, so a real error could not hide under it. */
	double maxd = 0.0; size_t worst = 0;
	for (size_t i = 0; i < n; i++) {
		double d = fabs((double)gfit->fdata[i] - (double)native[i]);
		if (d > maxd) { maxd = d; worst = i; }
	}
	cr_assert(maxd <= 1.2e-7, "max deviation %.9g at pixel %zu (%.9g vs %.9g)",
	          maxd, worst, gfit->fdata[worst], native[worst]);
	/* and it really did move: the pre-amend state was much further away */
	double pre = 0.0;
	for (size_t i = 0; i < n; i++) {
		double d = fabs((double)before[i] - (double)native[i]);
		if (d > pre) pre = d;
	}
	cr_assert(pre > 100.0 * maxd,
	          "the amend must account for the whole difference, not a sliver "
	          "(pre-amend %.9g vs post %.9g)", pre, maxd);

	g_free(before);
	g_free(native);
	golden_teardown(NULL, f);
}

/* A consumer pinned to the mask as it stood BEFORE the edited step used a
 * state the edit did not change, so it must be left alone. */
Test(nde_replay, a_pin_before_the_edit_is_not_disturbed) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);   /* pins mask@:1 */
	cr_assert_eq(apply_mask_op(&op_desc_mask_blur, blur_mask(1.f)), 0);   /* record 3 */

	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_input_pin *early = nde_record_input(g_ptr_array_index(snap, 1), "mask");
	cr_assert_not_null(early);
	cr_assert_eq(early->src_record_id, 1, "the stretch used the unblurred mask");
	g_ptr_array_unref(snap);

	/* the mask state that consumer pinned, as stored */
	fits *pinned_before = nde_checkpoint_output_get(1);
	cr_assert_not_null(pinned_before);
	size_t mn = (size_t)pinned_before->rx * pinned_before->ry;
	guint8 *keep = g_malloc(mn * sizeof(WORD));
	memcpy(keep, pinned_before->data, mn * sizeof(WORD));
	clearfits(pinned_before); free(pinned_before);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(3, "radius=5", &err), "amend failed: %s", err ? err : "?");
	unreserve_thread();

	fits *pinned_after = nde_checkpoint_output_get(1);
	cr_assert_not_null(pinned_after, "the earlier pin must still resolve");
	cr_assert_eq(memcmp(pinned_after->data, keep, mn * sizeof(WORD)), 0,
	             "editing a LATER mask step must not touch an earlier pin");
	clearfits(pinned_after); free(pinned_after);
	g_free(keep);
	golden_teardown(NULL, f);
}

/* ---------------- the image's own origin ---------------- */

/* The origin record says where the BASELINE came from, and a replay starts
 * from the baseline itself — so it must not be a chain member.  It is
 * DOCUMENT-scope, and the chain builder FAILS CLOSED on a DOCUMENT-scope
 * record it does not recognise: get this wrong and every image that has ever
 * been opened becomes unreplayable. */
Test(nde_replay, the_image_origin_record_is_not_a_step_to_re_run) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_neq(nde_capture_image_origin("file", "x.fit", "Opened 'x.fit'"), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);

	nde_chain *c = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert(c->replayable, "chain not replayable: %s",
	          c->reasons->len ? (char *)g_ptr_array_index(c->reasons, 0) : "none");
	cr_assert_eq(c->records->len, 1, "the stretch is the only step to re-run");
	nde_chain_free(c);

	/* And it offers nothing to edit: where an image came from is not a step
	 * whose parameters you can change, nor one you can remove. */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	const nde_record *origin = g_ptr_array_index(snap, 0);
	cr_assert_str_eq(origin->op_id, NDE_OP_IMAGE_ORIGIN);
	cr_assert(!nde_record_amendable(origin));
	cr_assert(!nde_record_deletable(origin));
	g_ptr_array_unref(snap);
	golden_teardown(NULL, f);
}

/* ------------- the reverse edge: image -> the masks built from it -------- */

/* Asinh params as the record's own serializer writes them, with beta varied:
 * an amend has to hand back a complete blob, not just the field that moved. */
#define ASINH_PARAMS(b) "beta=" b ";offset=0;human=0;clip_mode=1"

/* Stretch a gradient, then derive a mask from the stretched pixels.  Returns
 * the mask (caller g_free()s) and its size, and leaves no history behind. */
static guint8 *run_derived_mask(float beta, size_t *n_out) {
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(beta)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	guint8 *out = snapshot_mask(gfit, n_out);
	nde_history_attach(NULL);      /* also purges the checkpoint store */
	clearfits(f);
	free(f);
	gfit = NULL;
	return out;
}

/* The edge the pin exists for, read the other way round.  A mask derived from
 * an image records WHICH point of that image's history it read; editing the
 * image at or before that point changes the pixels the mask was built from,
 * so the mask has to be re-derived.  Nothing propagated that: the mask went on
 * describing pixels that no longer existed anywhere.
 *
 * Asserted against the same work done natively in that order — "it changed"
 * would not prove it changed to the right thing. */
Test(nde_replay, amending_the_image_re_derives_a_mask_built_from_it) {
	com.pref.nde_cache_mb = 256;
	size_t n_native = 0;
	guint8 *native = run_derived_mask(200.f, &n_native);

	/* the same document built at beta 10, then amended to 200 */
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);
	cr_assert_eq(n, n_native);
	cr_assert_neq(memcmp(before, native, n), 0,
	              "fixture: the two stretches must give different masks");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(1, ASINH_PARAMS("200"), &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();

	cr_assert_not_null(gfit->mask, "the commit emptied the mask slot");
	cr_assert_eq(memcmp(gfit->mask->data, native, n), 0,
	             "the mask must be re-derived from the amended pixels");
	g_free(before);
	g_free(native);
	golden_teardown(NULL, f);
}

/* Deleting a step is the same edge: what the mask read is gone. */
Test(nde_replay, deleting_an_image_step_re_derives_a_mask_built_from_it) {
	com.pref.nde_cache_mb = 256;
	size_t n = 0;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	/* the unstretched mask, for comparison, before the stretch exists */
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	guint8 *unstretched = snapshot_mask(gfit, &n);
	nde_history_attach(NULL);
	clearfits(f);
	free(f);

	f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	size_t n2 = 0;
	guint8 *stretched = snapshot_mask(gfit, &n2);
	cr_assert_eq(n, n2);
	cr_assert_neq(memcmp(unstretched, stretched, n), 0, "fixture: the stretch matters");

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_delete_execute(1, &err), "delete failed: %s", err ? err : "?");
	unreserve_thread();

	cr_assert_not_null(gfit->mask, "the commit emptied the mask slot");
	cr_assert_eq(memcmp(gfit->mask->data, unstretched, n), 0,
	             "without the stretch the mask must be the unstretched one");
	g_free(unstretched);
	g_free(stretched);
	golden_teardown(NULL, f);
}

/* A mask is a DIFFERENT ITEM with its own history, so replaying an image's
 * pixels has no business touching the mask slot — and the replay produces no
 * mask at all (mask_pin_clear runs after every masked member), so the
 * wholesale commit swap moved that nothing in and freed the user's mask with
 * the discarded pixels.  This mask is pinned to the image's BASELINE, which
 * the reverse cascade deliberately skips: what has to keep it is the commit. */
Test(nde_replay, a_history_edit_leaves_the_mask_slot_alone) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);

	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(2, ASINH_PARAMS("200"), &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();

	cr_assert_not_null(gfit->mask, "the edit emptied the mask slot");
	cr_assert_eq(memcmp(gfit->mask->data, before, n), 0,
	             "and nothing this mask was built from moved, so it must be "
	             "the same mask");
	g_free(before);
	golden_teardown(NULL, f);
}

/* A mask can read the image MORE THAN ONCE: built from one channel, then
 * rebuilt from another after the image was stretched again.  Each read is
 * pinned to its own point, and both halves of the machinery have to honour
 * that — the replay, which used to resolve only the FIRST record's pin and so
 * fed every later read the first one's pixels, and the reverse cascade, which
 * used to stop looking at a mask after its first pin and so missed an edit
 * that moved only the second. */
Test(nde_replay, a_mask_that_reads_the_image_twice_reads_each_point_it_pinned) {
	com.pref.nde_cache_mb = 256;

	/* the mask that sequence gives natively at beta 50 */
	fits *nat = flis_test_make_rgb_fits(16, 12, 0.f, 0.f, 0.f);
	fill_rgb_gradient(nat);
	gfit = nat;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(50.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(2)), 0);
	size_t n_native = 0;
	guint8 *native = snapshot_mask(gfit, &n_native);
	nde_history_attach(NULL);
	clearfits(nat);
	free(nat);

	/* the same thing with the SECOND stretch at 200, then amended to 50 */
	fits *f = flis_test_make_rgb_fits(16, 12, 0.f, 0.f, 0.f);
	fill_rgb_gradient(f);
	gfit = f;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(200.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(2)), 0);
	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);
	cr_assert_eq(n, n_native);
	cr_assert_neq(memcmp(before, native, n), 0,
	              "fixture: the second stretch must reach the second read");

	/* record 3 is the second stretch: the FIRST read is pinned before it and
	 * must be left alone, the second is pinned at it and must be rebuilt. */
	gchar *err = NULL;
	cr_assert(reserve_thread());
	cr_assert(nde_amend_execute(3, ASINH_PARAMS("50"), &err),
	          "amend failed: %s", err ? err : "?");
	unreserve_thread();

	cr_assert_not_null(gfit->mask);
	cr_assert_eq(memcmp(gfit->mask->data, native, n), 0,
	             "the second read must re-derive from the amended pixels");
	g_free(before);
	g_free(native);
	golden_teardown(NULL, f);
}


/* procmasksnag.png: an opaque step AFTER everything the mask's pin needs made
 * resolve_item_state refuse the whole input — "the source of this input
 * cannot be rebuilt: " with NOTHING after the colon, because a checkpointed
 * barrier sets replayable FALSE without adding a reason.  Only the replayed
 * prefix matters: the pin predates the barrier entirely. */
Test(nde_replay, a_mask_pin_survives_a_barrier_recorded_after_it) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);
	/* An opaque freehand with its restart checkpoint: a barrier on the image
	 * chain, well after the state the mask derives from. */
	gfit->fdata[0] += 0.125f;
	cr_assert(nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1,
	                             "freehand", gfit) > 0);

	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);
	gchar *err = NULL;
	cr_assert(reserve_thread());
	gboolean ok = nde_amend_execute(1, "channel=0;autostretch=0;invert=1;bitpix=8",
	                                &err);
	unreserve_thread();
	cr_assert(ok, "amending the mask's origin failed: %s", err ? err : "?");
	assert_mask_differs(gfit, before, n, "amend-past-barrier");
	g_free(before);
	golden_teardown(NULL, f);
}

/* The pin can also land AFTER a barrier: a mask built once an opaque step had
 * already run derives from pixels no from-baseline replay can produce.  The
 * barrier's output checkpoint is that state — restart there. */
Test(nde_replay, a_mask_built_after_a_barrier_restarts_from_its_checkpoint) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	gfit->fdata[5] += 0.25f;
	cr_assert(nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, -1,
	                             "freehand", gfit) > 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_input_pin *img = nde_record_input(g_ptr_array_index(snap, 2), "image");
	cr_assert_not_null(img);
	cr_assert_eq(img->src_record_id, 2, "fixture: the mask derives from the barrier's output");
	g_ptr_array_unref(snap);

	size_t n = 0;
	guint8 *before = snapshot_mask(gfit, &n);
	gchar *err = NULL;
	cr_assert(reserve_thread());
	gboolean ok = nde_amend_execute(3, "channel=0;autostretch=0;invert=1;bitpix=8",
	                                &err);
	unreserve_thread();
	cr_assert(ok, "amending a mask derived from a barrier state failed: %s",
	          err ? err : "?");
	assert_mask_differs(gfit, before, n, "amend-from-barrier-state");
	g_free(before);
	golden_teardown(NULL, f);
}

/* procmasksnag: clearing a USED processing mask sends it dormant, not away.
 * Its records stay editable — the rebuilt value refreshes the stored copies
 * its consumers read and the consumers recompute — and the edit must not
 * resurrect the mask the user cleared. */
Test(nde_replay, editing_a_cleared_masks_step_still_reaches_its_consumers) {
	com.pref.nde_cache_mb = 256;
	flis_layer_t *lay = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 12, 0.f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	fill_mono_gradient(gfit);

	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	gint mask_item = flis_layer_pmask_id(lay);
	cr_assert_neq(mask_item, 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_clear, NULL), 0);
	cr_assert_null(lay->fit->mask, "fixture: the mask is cleared");

	size_t n = (size_t)gfit->rx * gfit->ry;
	float *before = g_malloc(n * sizeof(float));
	memcpy(before, gfit->fdata, n * sizeof(float));

	gchar *err = NULL;
	cr_assert(reserve_thread());
	gboolean ok = nde_amend_execute(1, "channel=0;autostretch=0;invert=1;bitpix=8",
	                                &err);
	unreserve_thread();
	cr_assert(ok, "editing a cleared mask's step failed: %s", err ? err : "?");
	cr_assert_neq(memcmp(gfit->fdata, before, n * sizeof(float)), 0,
	              "the step that ran under the mask must be recomputed");
	cr_assert_null(lay->fit->mask, "the edit must not resurrect a cleared mask");

	g_free(before);
	nde_history_attach(NULL);
	gfit = NULL;
}

/* An 8-bit mask survives its trip through the snapshot store as an 8-bit mask.
 * It used to come back 16-bit: mask_to_fits widens 8-bit mask bytes into WORDs
 * and fits_to_mask reads orig_bitpix to decide what to rebuild, but the
 * snapshot did not carry orig_bitpix.  The pixels were unharmed (the widening
 * is exactly invertible), so the only visible trace was the blend then scaling
 * by 1/65535 instead of 1/255 — one float ulp on a handful of pixels. */
Test(nde_replay, an_eight_bit_masked_op_replays_bit_exactly) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	fits expected = { 0 };
	copyfits(gfit, &expected, CP_DEEPCOPY | CP_ALLOC, -1);
	nde_chain *chain = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert(chain->replayable);
	cr_assert(reserve_thread());
	gchar *err = NULL;
	fits *result = nde_chain_replay(chain, &err);
	unreserve_thread();
	cr_assert_not_null(result, "replay failed: %s", err ? err : "?");

	size_t n = (size_t)result->rx * result->ry;
	cr_assert_eq(memcmp(result->fdata, expected.fdata, n * sizeof(float)), 0,
	             "replaying under an 8-bit mask must reproduce the live pixels exactly");
	clearfits(result); free(result);
	nde_chain_free(chain);
	clearfits(&expected);
	golden_teardown(NULL, f);
}

/* Retention's promise is that an evicted pin is a VISIBLE loss.  This is the
 * other end of it: once the pinned mask is gone the step says so and becomes a
 * barrier, freezing what precedes it, rather than silently replaying without
 * the mask and producing different pixels (design note §7 — a silent
 * capability downgrade is worse than a refusal). */
Test(nde_replay, losing_a_pinned_mask_freezes_the_step_rather_than_lying) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;

	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_mask_op(&op_desc_mask_from_channel, from_channel(0)), 0);
	cr_assert_eq(apply_op_masked(&op_desc_asinh, asinh_beta(20.f)), 0);

	nde_chain *chain = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert(chain->replayable, "sanity: replayable while the pin is kept");
	cr_assert_eq(chain->tail_start, 0, "and nothing frozen");
	nde_chain_free(chain);

	/* Evict exactly what retention would evict: the checkpoint the mask pin
	 * names.  (Retention's own eviction ORDER is covered in
	 * test_nde_retention; here it is the consequence that matters.) */
	GPtrArray *snap = nde_history_snapshot(NULL);
	const nde_record *masked = g_ptr_array_index(snap, snap->len - 1);
	const nde_input_pin *pin = nde_record_input(masked, "mask");
	cr_assert_not_null(pin);
	cr_assert_neq(pin->src_record_id, 0, "this mask has its own history");
	gint64 pinned = pin->src_record_id;
	g_ptr_array_unref(snap);

	nde_checkpoint_output_drop(pinned);

	chain = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert_gt(chain->tail_start, 0,
	             "the step must freeze what precedes it, not replay maskless");
	guint8 flags = g_array_index(chain->member_flags, guint8, chain->tail_start - 1);
	cr_assert(flags & NDE_CHAIN_MEMBER_BARRIER,
	          "and the step itself becomes the barrier");
	nde_chain_free(chain);

	golden_teardown(NULL, f);
}

/* Losing an item's BASELINE — which is what happens to a layer consumed by a
 * merge, via nde_checkpoint_drop — blocks the whole chain rather than any
 * individual member.  The History keys its "cannot be re-run" marks off
 * exactly this, because tail_start and the member flags say nothing here:
 * without the chain-level verdict every row read as freely editable and
 * offered an Edit that could only fail. */
Test(nde_replay, losing_the_baseline_blocks_the_chain_not_a_member) {
	com.pref.nde_cache_mb = 256;
	fits *f = flis_test_make_mono_fits(16, 12, 0.f);
	fill_mono_gradient(f);
	gfit = f;
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(10.f)), 0);
	cr_assert_eq(apply_op_real(&op_desc_asinh, asinh_beta(20.f)), 0);

	nde_chain *before = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert(before->replayable);
	cr_assert_eq(before->reasons->len, 0);
	nde_chain_free(before);

	nde_checkpoint_drop(NDE_ITEM_IMAGE);

	nde_chain *after = nde_chain_build(NDE_ITEM_IMAGE);
	cr_assert(!after->replayable, "no baseline, nothing to replay from");
	cr_assert_gt(after->reasons->len, 0, "and it must say why");
	/* The members themselves are unremarkable — no barrier, empty frozen
	 * prefix — so a reader that only consults them learns nothing. */
	cr_assert_eq(after->tail_start, 0);
	for (guint i = 0; i < after->records->len; i++)
		cr_assert_eq(g_array_index(after->member_flags, guint8, i), 0,
		             "member %u should carry no flag of its own", i);
	nde_chain_free(after);

	golden_teardown(NULL, f);
}

