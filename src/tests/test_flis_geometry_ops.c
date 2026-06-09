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
 * test_flis_geometry_ops — calls the geometry image hooks directly to
 * verify they invoke the correct flis_update_layer_offset_after_*
 * helper after the in-place fit modification.
 *
 * §7 canvas decoupling: the helpers update only the active layer's
 * canvas position; other layers are never touched.  Canvas-scoped
 * operations live in flis_canvas_* and are covered separately.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "algos/geometry.h"

cominfo com;
fits *gfit;

TestSuite(flis_geometry_ops, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

static struct generic_img_args make_hook_args(void *user) {
	struct generic_img_args args = { 0 };
	args.fit = gfit;
	args.user = user;
	args.max_threads = 1;
	args.geometry_changing = TRUE;
	return args;
}

/* Crop on the bottom layer leaves other layers in place; the bottom
 * layer itself moves to the selection top-left. */
Test(flis_geometry_ops, crop_hook_only_moves_active) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(20, 20, 1.0f), "patch");
	top->position_x = 40; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	cr_assert_not_null(gfit);

	struct crop_args params = { .destroy_fn = NULL,
	    .area = { .x = 10, .y = 15, .w = 50, .h = 50 } };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(crop_image_hook_single(&args, gfit, 1), 0);

	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 10, "cropped layer lands at selection top-left");
	cr_assert_eq(base->position_y, 15);
	cr_assert_eq(top->position_x, 40, "non-active layer must not move");
	cr_assert_eq(top->position_y, 30);
	cr_assert_eq(gfit->rx, 50u);
	cr_assert_eq(gfit->ry, 50u);
}

/* Crop on a top layer: only that layer's position changes (translated
 * from fit coords to canvas coords by adding the layer's current pos). */
Test(flis_geometry_ops, crop_hook_top_translates_to_canvas) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *patch1 = flis_test_add_layer(
	    flis_test_make_mono_fits(40, 40, 0.5f), "patch1");
	flis_layer_t *patch2 = flis_test_add_layer(
	    flis_test_make_mono_fits(40, 40, 0.75f), "patch2");
	patch1->position_x = 10; patch1->position_y = 10;
	patch2->position_x = 50; patch2->position_y = 50;
	uniq_set_active_layer(com.uniq, 2);  /* patch2 active */
	gfit = flis_active_layer_fit();

	/* Selection in patch2 fit coords: (5, 5) → canvas (55, 55). */
	struct crop_args params = { .destroy_fn = NULL,
	    .area = { .x = 5, .y = 5, .w = 20, .h = 20 } };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(crop_image_hook_single(&args, gfit, 1), 0);
	cr_assert_eq(patch2->position_x, 55);
	cr_assert_eq(patch2->position_y, 55);
	cr_assert_eq(patch1->position_x, 10);  /* unchanged */
	cr_assert_eq(patch1->position_y, 10);
}

/* Rotate 180° on the bottom layer: the layer's centre stays put.  Other
 * layers are not touched. */
Test(flis_geometry_ops, rotate_180_hook_preserves_active_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 20;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct rotation_args params = { .destroy_fn = NULL,
	    .area = { .x = 0, .y = 0, .w = gfit->rx, .h = gfit->ry },
	    .angle = 180.0,
	    .interpolation = -1,
	    .cropped = 0,
	    .clamp = FALSE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(rotation_image_hook(&args, gfit, 1), 0);
	/* base dims unchanged (100x100 rotated by 180 = same dims); centre
	 * preserved → base position stays 0.  Top layer untouched. */
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 0);
	cr_assert_eq(base->position_y, 0);
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 20);
}

/* Resize doubles the bottom layer's dimensions: centre preserved, so
 * position shifts by (-old/2) = (-50,-50) from origin. */
Test(flis_geometry_ops, resize_hook_preserves_active_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct resample_args params = { .destroy_fn = NULL,
	    .toX = 200, .toY = 200,
	    .interpolation = OPENCV_NEAREST,
	    .clamp = FALSE,
	    .update_wcs = FALSE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(resample_image_hook(&args, gfit, 1), 0);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	/* (100-200)/2 = -50 added to each axis → base centre stays at (50,50). */
	cr_assert_eq(base->position_x, -50);
	cr_assert_eq(base->position_y, -50);
	cr_assert_eq(top->position_x, 40);  /* untouched */
	cr_assert_eq(top->position_y, 30);
}

/* Mirror is position-invariant: active layer stays where it was, and
 * other layers stay where they were. */
Test(flis_geometry_ops, mirrorx_hook_position_invariant) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct mirror_args params = { .destroy_fn = NULL, .x_axis = TRUE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(mirrorx_image_hook(&args, gfit, 1), 0);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 0);
	cr_assert_eq(base->position_y, 0);
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 30);
}

Test(flis_geometry_ops, mirrory_hook_position_invariant) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(10, 10, 0.5f), "patch");
	top->position_x = 20; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 1);  /* top active */
	gfit = flis_active_layer_fit();

	struct mirror_args params = { .destroy_fn = NULL, .x_axis = FALSE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(mirrory_image_hook(&args, gfit, 1), 0);
	cr_assert_eq(top->position_x, 20);
	cr_assert_eq(top->position_y, 30);
}

/* Binning halves the bottom layer's dimensions: centre preserved, so
 * position shifts by (100-50)/2 = +25 on each axis. */
Test(flis_geometry_ops, binning_hook_preserves_active_centre) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();

	struct binning_args params = { .destroy_fn = NULL,
	    .factor = 2, .mean = TRUE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(binning_image_hook(&args, gfit, 1), 0);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_eq(base->position_x, 25);
	cr_assert_eq(base->position_y, 25);
	cr_assert_eq(top->position_x, 40);  /* untouched */
	cr_assert_eq(top->position_y, 30);
}

/* Sanity: when no FLIS is loaded the hooks must not crash and simply
 * forward the underlying fits op. */
Test(flis_geometry_ops, hook_is_safe_when_no_flis_loaded) {
	fits *f = flis_test_make_mono_fits(20, 20, 0.5f);
	gfit = f;

	struct mirror_args params = { .destroy_fn = NULL, .x_axis = TRUE };
	struct generic_img_args args = make_hook_args(&params);

	cr_assert_eq(mirrorx_image_hook(&args, gfit, 1), 0);
	cr_assert_eq(gfit->rx, 20u);
	cr_assert_eq(gfit->ry, 20u);

	clearfits(f);
	free(f);
	gfit = NULL;
}
