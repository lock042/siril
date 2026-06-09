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
 * test_flis_mask_route — exercises the §5.2 FLIS layermask routing in
 * generic_mask_worker.  When args->target_layer_id is set and the
 * mask_hook produces a mask on args->fit->mask, the worker moves that
 * mask onto the named layer's lmask slot and clears args->fit->mask.
 *
 * The mask creation path itself (mask_from_channel_hook etc.) is
 * tested elsewhere; this file uses a trivial test hook that fills a
 * known mask so the assertions can focus on routing semantics.
 *
 * The worker is called directly on the test thread (com.script +
 * args->command run the synchronous-completion / no-undo path).
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/masks.h"

cominfo com;
fits *gfit;

TestSuite(flis_mask_route, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Test hook: creates a constant 8-bit mask on the fit. */
static int fill_mask_hook(struct generic_mask_args *args) {
	fits *f = args->fit;
	size_t n = (size_t)f->rx * f->ry;
	f->mask = calloc(1, sizeof(mask_t));
	f->mask->bitpix = 8;
	f->mask->data = malloc(n);
	memset(f->mask->data, 0xAA, n);
	return 0;
}

/* Hook that does nothing — used to verify "no routing without -layermask=". */
static int noop_hook(struct generic_mask_args *args) {
	(void)args;
	return 0;
}

/* Routing happens: mask hook produces fit->mask, worker moves it to the
 * target layer's lmask and clears the processing mask. */
Test(flis_mask_route, target_id_set_routes_to_layer_lmask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.5f), "patch");
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();
	cr_assert_not_null(gfit);

	cr_assert_null(top->lmask, "precondition: top has no lmask");

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test routing to top";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = top->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	/* Mask now lives on top->lmask; gfit->mask cleared. */
	cr_assert_null(gfit->mask, "processing mask should be cleared");
	cr_assert_not_null(top->lmask, "target layer lmask should be populated");
	cr_assert_eq(top->lmask->w, 16u);
	cr_assert_eq(top->lmask->h, 16u);
	cr_assert_eq(top->lmask->bitpix, 8);
	const uint8_t *p = (const uint8_t *)top->lmask->data;
	cr_assert_eq(p[0], 0xAA);
	cr_assert_eq(p[16 * 16 - 1], 0xAA);
}

/* Routing target a different (base) layer than the active layer:
 * routing key is item_id, not "active layer". */
Test(flis_mask_route, target_id_can_be_non_active_layer) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, 0.5f), "patch");
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test routing to base";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = base->item_id;  /* not the active one */
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask);
	cr_assert_not_null(base->lmask);
	cr_assert_null(top->lmask,  "non-target layer must remain untouched");
}

/* When target_layer_id == 0 (the default), the worker leaves the mask on
 * the processing slot and activates it (mask_creation path). */
Test(flis_mask_route, target_id_zero_keeps_processing_mask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.5f), "base");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "test default routing";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = 0;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_not_null(gfit->mask, "processing mask should remain");
	cr_assert_eq(gfit->mask->bitpix, 8);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_null(base->lmask, "no lmask routing should have happened");
}

/* Dimension mismatch: target layer is smaller than args->fit; routing
 * is refused (logged), mask stays on processing slot — a safer fallback
 * than dropping the mask on the floor. */
Test(flis_mask_route, dim_mismatch_falls_back_to_processing_mask) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.0f), "base");
	flis_layer_t *small = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "small");
	uniq_set_active_layer(com.uniq, 0);  /* base active (16x16) */
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = fill_mask_hook;
	args->description      = "mismatch test";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = small->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_not_null(gfit->mask, "mask should fall back to processing slot");
	cr_assert_null(small->lmask,   "mismatch layer must not receive mask");
}

/* No mask produced by the hook (hook returned success but didn't set
 * fit->mask): routing block is skipped, no crash. */
Test(flis_mask_route, no_mask_produced_no_op) {
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(8, 8, 0.5f), "x");
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct generic_mask_args *args = calloc(1, sizeof(*args));
	args->fit              = gfit;
	args->mask_hook        = noop_hook;
	args->description      = "no mask test";
	args->command          = TRUE;
	args->mask_creation    = TRUE;
	args->target_layer_id  = top->item_id;
	args->max_threads      = 1;

	generic_mask_worker(args);

	cr_assert_null(gfit->mask);
	cr_assert_null(top->lmask);
}
