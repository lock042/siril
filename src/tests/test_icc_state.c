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
 * test_icc_state — ICC assign / remove / convert operate on com.uniq.
 *
 * Regression coverage for the swap-path ICC loss: the old
 * icc_assign/remove/convert hooks rode generic_image_worker, whose swap
 * path hands the hook a private deep copy of gfit.  With ICC state on
 * com.uniq (no per-fits fields), the copy had no ICC identity, so the
 * profile change was silently discarded while "assignment succeeded"
 * was logged.  Assign/Remove now run through icc_state_worker (com.uniq
 * direct, no pixel copy); Convert's hooks read the source profile from
 * com.uniq and land the new profile there explicitly.
 *
 * These tests drive the worker / hooks directly (the command functions
 * add only argument parsing plus start_in_new_thread, which needs the
 * live processing-thread machinery).  The workers run with
 * command=FALSE in script mode, where siril_add_idle is a no-op, so the
 * args survive the call for inspection.
 */

#include <criterion/criterion.h>
#include <lcms2.h>
#include "flis_test_helpers.h"
#include "core/icc_profile.h"
#include "core/processing.h"

cominfo com;
fits *gfit;

static void icc_test_init(void) {
	flis_test_init_com();
	/* Skip the display-transform rebuild inside refresh_icc_transforms
	 * (needs monitor profiles the test env doesn't have). */
	com.headless = TRUE;
}

TestSuite(icc_state, .init = icc_test_init, .fini = flis_test_cleanup_com);

/* Install a plain (non-FLIS) image as the current image. */
static fits *load_plain(int chans) {
	fits *f = chans == 1 ? flis_test_make_mono_fits(8, 8, 0.5f)
	                     : flis_test_make_rgb_fits(8, 8, 0.5f, 0.5f, 0.5f);
	cr_assert_not_null(f);
	gfit = f;
	com.uniq->fit = f;
	com.uniq->chans = chans;
	com.uniq->filename = g_strdup("icc_test");
	return f;
}

/* Run icc_state_worker synchronously; takes ownership of @profile
 * (NULL = remove).  Returns the worker's retval. */
static int run_icc_state(cmsHPROFILE profile) {
	struct icc_state_args *args = calloc(1, sizeof(struct icc_state_args));
	args->destroy_fn = free_icc_state_args;
	args->profile = profile;
	gpointer r = icc_state_worker(args);
	int rc = GPOINTER_TO_INT(r);
	free_icc_state_args(args);
	return rc;
}

static struct icc_data *make_icc_data(cmsHPROFILE profile) {
	struct icc_data *icc = calloc(1, sizeof(struct icc_data));
	icc->destroy_fn = free_icc_data;
	icc->profile = profile;
	icc->intent = com.pref.icc.export_intent;
	return icc;
}

/* sRGB TRC 0.5 → linear ≈ 0.2140 */
#define SRGB_LINEAR_OF_HALF 0.2140f

Test(icc_state, assign_lands_on_com_uniq) {
	load_plain(3);
	cr_assert_null(current_icc_profile());
	cr_assert_not(current_image_color_managed());
	cr_assert_eq(run_icc_state(srgb_trc()), 0);
	cr_assert_not_null(current_icc_profile());
	cr_assert(current_image_color_managed());
}

Test(icc_state, remove_clears_com_uniq) {
	load_plain(3);
	cr_assert_eq(run_icc_state(srgb_trc()), 0);
	cr_assert_eq(run_icc_state(NULL), 0);
	cr_assert_null(current_icc_profile());
	cr_assert_not(current_image_color_managed());
}

Test(icc_state, assign_channel_mismatch_refused) {
	load_plain(1);
	cr_assert_eq(run_icc_state(srgb_trc()), 1); /* RGB profile, mono image */
	cr_assert_null(current_icc_profile());
	cr_assert_not(current_image_color_managed());
}

Test(icc_state, flis_assign_is_composite_aware) {
	/* Mono active layer, but the FLIS composite is always RGB: an RGB
	 * profile must be accepted and a gray one refused. */
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "L");
	uniq_set_active_layer(com.uniq, 0);
	cr_assert_eq(run_icc_state(gray_srgbtrc()), 1);
	cr_assert_null(current_icc_profile());
	cr_assert_eq(run_icc_state(srgb_trc()), 0);
	cr_assert_not_null(current_icc_profile());
	cr_assert_eq(cmsChannelsOf(cmsGetColorSpace(current_icc_profile())), 3);
}

Test(icc_state, convert_hook_lands_profile_from_swap_copy) {
	/* Emulate generic_image_worker's swap path: the hook receives a
	 * private deep copy of gfit.  The profile must land on com.uniq and
	 * the copy's pixels must be transformed; gfit itself stays untouched
	 * (the worker's swap brings the copy in afterwards). */
	fits *f = load_plain(3);
	cr_assert_eq(run_icc_state(srgb_trc()), 0);
	gchar *before = siril_color_profile_get_description(current_icc_profile());

	fits *clone = calloc(1, sizeof(fits));
	cr_assert_eq(copyfits(f, clone, CP_DEEPCOPY | CP_ALLOC, -1), 0);

	struct icc_data *icc = make_icc_data(srgb_linear());
	struct generic_img_args gargs = { 0 };
	gargs.fit = gfit;
	gargs.user = icc;
	cr_assert_eq(icc_convert_to_hook(&gargs, clone, 1), 0);

	gchar *after = siril_color_profile_get_description(current_icc_profile());
	cr_assert(g_strcmp0(before, after) != 0,
	          "profile on com.uniq did not change (before='%s' after='%s')",
	          before, after);
	cr_assert(current_image_color_managed());
	cr_assert(flis_test_approx_eq(clone->fpdata[0][0], SRGB_LINEAR_OF_HALF, 0.02f),
	          "swap copy not transformed: %f", clone->fpdata[0][0]);
	cr_assert(flis_test_approx_eq(f->fpdata[0][0], 0.5f, 1e-6f),
	          "gfit must stay untouched until the swap");

	g_free(before);
	g_free(after);
	free_icc_data(icc);
	clearfits(clone);
	free(clone);
}

Test(icc_state, flis_convert_transforms_every_layer) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.5f, 0.5f, 0.5f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "top");
	uniq_set_active_layer(com.uniq, 0);
	cr_assert_eq(run_icc_state(srgb_trc()), 0);
	gchar *before = siril_color_profile_get_description(current_icc_profile());

	struct icc_data *icc = make_icc_data(srgb_linear());
	struct generic_layer_args largs = { 0 };
	largs.user = icc;

	/* The layer worker holds the stack writer lock across the hook. */
	flis_stack_writer_lock();
	int rc = icc_convert_flis_layer_hook(&largs);
	flis_stack_writer_unlock();
	cr_assert_eq(rc, 0);

	gchar *after = siril_color_profile_get_description(current_icc_profile());
	cr_assert(g_strcmp0(before, after) != 0);

	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *top = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert(flis_test_approx_eq(base->fit->fpdata[0][0], SRGB_LINEAR_OF_HALF, 0.02f),
	          "base layer not transformed: %f", base->fit->fpdata[0][0]);
	cr_assert(flis_test_approx_eq(top->fit->fdata[0], SRGB_LINEAR_OF_HALF, 0.02f),
	          "non-active mono layer not transformed: %f", top->fit->fdata[0]);

	g_free(before);
	g_free(after);
	free_icc_data(icc);
}

Test(icc_state, convert_refused_without_profile) {
	load_plain(3);
	struct icc_data *icc = make_icc_data(srgb_linear());
	struct generic_img_args gargs = { 0 };
	gargs.fit = gfit;
	gargs.user = icc;
	cr_assert_eq(icc_convert_to_hook(&gargs, gfit, 1), 1);
	free_icc_data(icc);
}
