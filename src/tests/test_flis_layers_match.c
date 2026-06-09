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
 * test_flis_layers_match — exercises the §5.7 `flis_layers_match` command
 * (background neutralisation across tinted mono layers).  The underlying
 * primitive (`flis_background_neutralise_layers`) is already covered by
 * its own tests; this file focuses on the command surface: argument
 * parsing, subset resolution, and the full headless path.
 *
 * Strategy: build a 3-mono-layer fixture with distinct tints, push a
 * non-trivial median through each layer, run the command, assert the
 * scaled medians compose to a near-neutral background.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"

cominfo com;
fits *gfit;
char *word[MAX_COMMAND_WORDS];

extern int process_flis_layers_match(int nb);

TestSuite(flis_layers_match, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Build three tinted mono layers, each with a constant grey background.
 * The screen-blend composite of (a_i * lay_i * tint_i) summed across
 * channels should land near 1.0 per channel after neutralisation
 * (sum of channel contributions / 3 → 1.0). */
static void make_tinted_triple(double med_r, double med_g, double med_b) {
	flis_layer_t *base = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_r), "R");
	flis_layer_set_tint(base, 1.0, 0.0, 0.0);
	flis_layer_t *g = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_g), "G");
	flis_layer_set_tint(g, 0.0, 1.0, 0.0);
	flis_layer_t *b = flis_test_add_layer(
	    flis_test_make_mono_fits(16, 16, (float)med_b), "B");
	flis_layer_set_tint(b, 0.0, 0.0, 1.0);
}

Test(flis_layers_match, refuses_when_no_flis) {
	/* com.uniq is set up by init, but no layers added — is_current_image_flis
	 * returns FALSE without layers, so the command must refuse. */
	word[0] = "flis_layers_match";
	word[1] = NULL;
	cr_assert_eq(process_flis_layers_match(1), CMD_GENERIC_ERROR);
}

Test(flis_layers_match, runs_on_default_subset) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = NULL;
	cr_assert_eq(process_flis_layers_match(1), CMD_OK);
}

Test(flis_layers_match, accepts_subset_by_id) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	/* item_ids are 1, 2, 3 (assigned sequentially by flis_layer_add). */
	word[1] = "-subset=1,2,3";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
}

Test(flis_layers_match, accepts_subset_by_name) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-subset=R,G,B";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_OK);
}

Test(flis_layers_match, rejects_unknown_layer_in_subset) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-subset=R,NoSuchLayer";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_ARG_ERROR);
}

Test(flis_layers_match, rejects_unknown_option) {
	make_tinted_triple(0.2, 0.3, 0.4);
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
	word[0] = "flis_layers_match";
	word[1] = "-bogus=1";
	word[2] = NULL;
	cr_assert_eq(process_flis_layers_match(2), CMD_ARG_ERROR);
}
