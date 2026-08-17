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
 * test_icc_state — display-transform cache lifecycle (#1948).
 *
 * Master subset of the flis-gtk4-nde test_icc_state.c: the ICC
 * assign / remove / convert worker tests there depend on that branch's
 * com.uniq ICC storage and FLIS test helpers, so only the
 * display-transform cache test is carried here.  Keep it textually
 * identical to the branch copy to ease merges.
 */

#include <criterion/criterion.h>
#include <string.h>
#include <lcms2.h>
#include "core/siril.h"
#include "core/icc_profile.h"

cominfo com;
fits *gfit;

static fits icc_test_fit;

static void icc_test_init(void) {
	memset(&com, 0, sizeof(com));
	memset(&icc_test_fit, 0, sizeof(icc_test_fit));
	gfit = &icc_test_fit;
	/* Skip the display-transform rebuild inside refresh_icc_transforms
	 * (needs monitor profiles the test env doesn't have). */
	com.headless = TRUE;
}

TestSuite(icc_state, .init = icc_test_init);

/* #1948 display-transform caches: clear_proofing_transforms() must delete all
 * three cached transforms and re-arm the lazy gamut build (gamut_transform_
 * tried), and get_gamut_transform() — called under the display transform lock,
 * per its contract — must cache a failed build until the caches are cleared.
 * In this headless test env there is no monitor profile, so the gamut build
 * legitimately fails: exactly the case gamut_transform_tried exists for. */
Test(icc_state, display_transform_caches_clear_and_rearm) {
	/* Stand-in transforms; trivial but real, so cmsDeleteTransform is legal. */
	cmsHPROFILE p = cmsCreate_sRGBProfile();
	cr_assert_not_null(p);
	com.gui_icc.proofing_transform = cmsCreateTransform(p, TYPE_RGB_8,
			p, TYPE_RGB_8, INTENT_PERCEPTUAL, 0);
	com.gui_icc.proofing_lut_transform = cmsCreateTransform(p, TYPE_RGB_FLT,
			p, TYPE_RGB_8, INTENT_PERCEPTUAL, 0);
	com.gui_icc.gamut_transform = cmsCreateTransform(p, TYPE_RGB_8,
			p, TYPE_RGB_8, INTENT_PERCEPTUAL, 0);
	cr_assert_not_null(com.gui_icc.proofing_transform);
	cr_assert_not_null(com.gui_icc.proofing_lut_transform);
	cr_assert_not_null(com.gui_icc.gamut_transform);
	com.gui_icc.gamut_transform_tried = TRUE;

	lock_display_transform();
	clear_proofing_transforms();
	unlock_display_transform();
	cr_assert_null(com.gui_icc.proofing_transform);
	cr_assert_null(com.gui_icc.proofing_lut_transform);
	cr_assert_null(com.gui_icc.gamut_transform);
	cr_assert_not(com.gui_icc.gamut_transform_tried);

	/* Lazy build: no monitor profile → build fails, result (NULL) is cached
	 * via the tried flag so it is not retried every lookup... */
	lock_display_transform();
	cr_assert_null(get_gamut_transform());
	cr_assert(com.gui_icc.gamut_transform_tried);
	cr_assert_null(get_gamut_transform());
	/* ...until a cache clear re-arms it. */
	clear_proofing_transforms();
	cr_assert_not(com.gui_icc.gamut_transform_tried);
	unlock_display_transform();

	cmsCloseProfile(p);
}
