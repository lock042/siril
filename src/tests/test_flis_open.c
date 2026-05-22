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
 * test_flis_open — stage 2 integration tests.  Covers:
 *   • readfits() FLIS dispatch: a file with FLIS=T in HDU 0 is routed to
 *     load_flis, populating com.uniq->layers without going through the
 *     normal FITS read path.
 *   • flis_promote: turns a loaded plain FITS into a 1-layer FLIS.
 *   • flis_promote idempotent refusal: calling on an already-FLIS image
 *     succeeds (with a warning) and changes nothing.
 *   • Full promote → save → load round-trip through the readfits dispatch
 *     (the test_flis_roundtrip tests go through save_flis/load_flis
 *     directly; this one exercises the public open_single_image path
 *     that scripts and the GUI actually use).
 */

#include <criterion/criterion.h>
#include <unistd.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "io/image_format_fits.h"

cominfo com;
fits *gfit;
/* needed by process_flis_promote — arg parser reads word[1..nb-1] */
char *word[MAX_COMMAND_WORDS];

static char *tmpdir = NULL;
static char *tmppath = NULL;

static void setup(void) {
	flis_test_init_com();
	tmpdir = g_dir_make_tmp("flis-open-XXXXXX", NULL);
	tmppath = g_build_filename(tmpdir, "test.flis", NULL);
}
static void teardown(void) {
	if (tmppath) { g_unlink(tmppath); g_free(tmppath); tmppath = NULL; }
	if (tmpdir)  { g_rmdir(tmpdir);   g_free(tmpdir);   tmpdir = NULL; }
	flis_test_cleanup_com();
}

TestSuite(flis_open, .init = setup, .fini = teardown);

/* Build a 2-layer fixture in memory, save it, free, reload via readfits
 * (the public path scripts and the GUI take), and assert the layer stack
 * is restored.  This is the key §2.1 round-trip — exercises the FLIS=T
 * header check inside readfits, not just save_flis/load_flis directly. */
Test(flis_open, readfits_dispatches_flis_files_to_load_flis) {
	/* Build and save a 2-layer FLIS */
	flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.1f, 0.2f, 0.3f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.7f), "top");
	cr_assert_eq(save_flis(tmppath), 0);
	cr_assert_eq(flis_layer_count(), 2);

	/* Reload through readfits — the public path.  Allocate gfit so
	 * readfits has somewhere to write to (the FLIS dispatch in readfits
	 * requires fit == gfit). */
	flis_free_layers(com.uniq);
	free(com.uniq); com.uniq = NULL;
	gfit = calloc(1, sizeof(fits));
	cr_assert_eq(readfits(tmppath, gfit, NULL, FALSE), 0);

	/* load_flis should have populated com.uniq via the readfits delegation,
	 * not the normal FITS read path. */
	cr_assert_not_null(com.uniq, "com.uniq should be populated by load_flis");
	cr_assert_eq(flis_layer_count(), 2,
	             "two layers should round-trip through readfits dispatch");
	cr_assert(is_current_image_flis());
}

/* Plain-FITS readfits unchanged: writing a non-FLIS FITS and reading it
 * back via readfits should NOT populate com.uniq->layers. */
Test(flis_open, readfits_plain_fits_not_treated_as_flis) {
	gfit = flis_test_make_mono_fits(4, 4, 0.5f);
	gchar *plain_fits = g_build_filename(tmpdir, "plain.fit", NULL);
	cr_assert_eq(savefits(plain_fits, gfit), 0);
	clearfits(gfit);
	free(gfit);
	gfit = calloc(1, sizeof(fits));

	cr_assert_eq(readfits(plain_fits, gfit, NULL, FALSE), 0);
	cr_assert(!is_current_image_flis() || !com.uniq || !com.uniq->layers,
	          "plain FITS must not be auto-treated as FLIS");
	g_unlink(plain_fits);
	g_free(plain_fits);
}

/* ----- flis_promote command coverage ----- */

Test(flis_open, promote_turns_plain_fits_into_single_layer_flis) {
	/* Simulate a freshly-loaded plain FITS: gfit allocated, com.uniq
	 * present with chans set but no layers (the post-create_uniq_from_gfit
	 * state for a plain FITS). */
	gfit = flis_test_make_mono_fits(8, 8, 0.5f);
	com.uniq->fit = gfit;
	com.uniq->chans = 1;
	cr_assert(!is_current_image_flis());

	word[0] = "flis_promote";
	word[1] = NULL;
	cr_assert_eq(process_flis_promote(1), CMD_OK);
	cr_assert(is_current_image_flis());
	cr_assert_eq(flis_layer_count(), 1);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_str_eq(base->layer_name, "Background");
}

Test(flis_open, promote_with_name_uses_given_name) {
	gfit = flis_test_make_mono_fits(8, 8, 0.5f);
	com.uniq->fit = gfit;
	com.uniq->chans = 1;

	word[0] = "flis_promote";
	word[1] = "-name=Jupiter";
	word[2] = NULL;
	cr_assert_eq(process_flis_promote(2), CMD_OK);
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_str_eq(base->layer_name, "Jupiter");
}

Test(flis_open, promote_idempotent_when_already_flis) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "existing");
	cr_assert(is_current_image_flis());

	word[0] = "flis_promote";
	word[1] = NULL;
	/* Second invocation succeeds with a warning; layer stack unchanged. */
	cr_assert_eq(process_flis_promote(1), CMD_OK);
	cr_assert_eq(flis_layer_count(), 1, "no extra layer should be added");
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_str_eq(base->layer_name, "existing");
}

Test(flis_open, promote_rejects_unknown_argument) {
	gfit = flis_test_make_mono_fits(8, 8, 0.5f);
	com.uniq->fit = gfit;
	com.uniq->chans = 1;

	word[0] = "flis_promote";
	word[1] = "-bogus=42";
	word[2] = NULL;
	cr_assert_eq(process_flis_promote(2), CMD_ARG_ERROR);
}

/* Data-driven save: a multi-layer FLIS saved to a plain .fit name still
 * comes back as multi-layer.  Format follows the data state (which is
 * FLIS, with >1 layer), not the chosen extension. */
Test(flis_open, save_flis_preserves_layers_with_fit_extension) {
	/* Two-layer FLIS in memory */
	flis_test_add_layer(flis_test_make_rgb_fits(4, 4, 0.1f, 0.1f, 0.1f), "base");
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.7f), "top");
	cr_assert_eq(flis_layer_count(), 2);

	/* Write to a non-.flis extension — save_flis is still the right
	 * writer because the data state is FLIS.  The file is a valid FITS
	 * container by spec, so the .fit extension is correct. */
	gchar *fit_path = g_build_filename(tmpdir, "test.fit", NULL);
	cr_assert_eq(save_flis(fit_path), 0);

	/* Round-trip through readfits — the FLIS=T header is the signal,
	 * not the extension.  Both layers should come back. */
	flis_free_layers(com.uniq);
	free(com.uniq); com.uniq = NULL;
	gfit = calloc(1, sizeof(fits));
	cr_assert_eq(readfits(fit_path, gfit, NULL, FALSE), 0);
	cr_assert_eq(flis_layer_count(), 2,
	             ".fit-extensioned FLIS must round-trip with all layers");

	g_unlink(fit_path);
	g_free(fit_path);
}

/* Full integration: promote, save (through save_flis), close, reopen via
 * readfits dispatch.  Mirrors what the e2e siril-cli script does. */
Test(flis_open, promote_then_save_then_readfits_roundtrip) {
	gfit = flis_test_make_rgb_fits(8, 8, 0.2f, 0.3f, 0.4f);
	com.uniq->fit = gfit;
	com.uniq->chans = 3;

	word[0] = "flis_promote";
	word[1] = "-name=base";
	word[2] = NULL;
	cr_assert_eq(process_flis_promote(2), CMD_OK);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	free(com.uniq); com.uniq = NULL;
	gfit = calloc(1, sizeof(fits));

	cr_assert_eq(readfits(tmppath, gfit, NULL, FALSE), 0);
	cr_assert(is_current_image_flis());
	cr_assert_eq(flis_layer_count(), 1);
	flis_layer_t *reloaded = (flis_layer_t *)com.uniq->layers->data;
	cr_assert_str_eq(reloaded->layer_name, "base");
	cr_assert_eq(reloaded->fit->naxes[2], 3);
	cr_assert_float_eq(reloaded->fit->fpdata[0][0], 0.2f, 1e-5);
	cr_assert_float_eq(reloaded->fit->fpdata[1][0], 0.3f, 1e-5);
	cr_assert_float_eq(reloaded->fit->fpdata[2][0], 0.4f, 1e-5);
}
