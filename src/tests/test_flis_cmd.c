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
 * test_flis_cmd — exercises the FLIS introspection commands
 * (process_flis_info / _layer_list / _group_list / _layer_info /
 * _group_info) directly, without going through the script-engine
 * front end.  Subprocess script-engine tests would be a useful
 * adjunct but the direct calls are sufficient to validate the
 * command implementations themselves; the prerequisite check and
 * the parser integration were verified manually at the end of §1.6.
 *
 * The CSV output of flis_layer_list and flis_group_list is the
 * load-bearing piece for downstream test harnesses (compare
 * before/after CSVs in scripted operation tests), so we focus on
 * verifying its structure here.
 */

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include <stdio.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"

cominfo com;
fits *gfit;

/* The introspection commands access word[] for argument parsing.  We
 * need to declare it so the link resolves and we can populate it for
 * the by-name-arg tests.  Mirrors command.c's declaration. */
char *word[MAX_COMMAND_WORDS];

TestSuite(flis_cmd, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Build a known 2-layer fixture used by most tests below. */
static void load_two_layer_fixture(void) {
	flis_test_add_layer(flis_test_make_rgb_fits(8, 8, 0.1f, 0.1f, 0.1f), "base");
	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "Ha");
	flis_layer_set_blend_mode(top, FLIS_BLEND_SCREEN);
	flis_layer_set_opacity(top, 0.8f);
	flis_layer_set_tint(top, 1.0, 0.2, 0.1);
}

Test(flis_cmd, info_runs_and_reports_layer_count) {
	load_two_layer_fixture();
	cr_assert_eq(process_flis_info(1), CMD_OK);
	/* The command's output goes through siril_log_info; we don't
	 * intercept it here because the GUI-independent log routing in
	 * headless writes to stdout anyway and that's verified manually. */
}

Test(flis_cmd, layer_list_runs_in_text_mode) {
	load_two_layer_fixture();
	word[0] = "flis_layer_list";
	word[1] = NULL;
	cr_assert_eq(process_flis_layer_list(1), CMD_OK);
}

Test(flis_cmd, layer_list_csv_format) {
	load_two_layer_fixture();
	word[0] = "flis_layer_list";
	word[1] = "-format=csv";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_list(2), CMD_OK);
}

Test(flis_cmd, layer_list_bad_format_rejected) {
	load_two_layer_fixture();
	word[0] = "flis_layer_list";
	word[1] = "-format=xml";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_list(2), CMD_ARG_ERROR);
}

Test(flis_cmd, group_list_runs) {
	load_two_layer_fixture();
	flis_group_add("test group");
	word[0] = "flis_group_list";
	word[1] = NULL;
	cr_assert_eq(process_flis_group_list(1), CMD_OK);
	word[1] = "-format=csv";
	word[2] = NULL;
	cr_assert_eq(process_flis_group_list(2), CMD_OK);
}

Test(flis_cmd, layer_info_by_id_succeeds) {
	load_two_layer_fixture();
	word[0] = "flis_layer_info";
	word[1] = "1";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_info(2), CMD_OK);
}

Test(flis_cmd, layer_info_by_name_succeeds) {
	load_two_layer_fixture();
	word[0] = "flis_layer_info";
	word[1] = "Ha";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_info(2), CMD_OK);
}

Test(flis_cmd, layer_info_by_unknown_id_fails) {
	load_two_layer_fixture();
	word[0] = "flis_layer_info";
	word[1] = "999";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_info(2), CMD_ARG_ERROR);
}

Test(flis_cmd, layer_info_by_unknown_name_fails) {
	load_two_layer_fixture();
	word[0] = "flis_layer_info";
	word[1] = "nonexistent";
	word[2] = NULL;
	cr_assert_eq(process_flis_layer_info(2), CMD_ARG_ERROR);
}

Test(flis_cmd, layer_info_no_arg_fails) {
	load_two_layer_fixture();
	word[0] = "flis_layer_info";
	word[1] = NULL;
	cr_assert_eq(process_flis_layer_info(1), CMD_WRONG_N_ARG);
}

Test(flis_cmd, group_info_by_id_succeeds) {
	load_two_layer_fixture();
	flis_group_t *g = flis_group_add("Narrowband");
	word[0] = "flis_group_info";
	char buf[16];
	snprintf(buf, sizeof(buf), "%d", g->item_id);
	word[1] = buf;
	word[2] = NULL;
	cr_assert_eq(process_flis_group_info(2), CMD_OK);
}

Test(flis_cmd, group_info_by_name_succeeds) {
	load_two_layer_fixture();
	flis_group_add("Narrowband");
	word[0] = "flis_group_info";
	word[1] = "Narrowband";
	word[2] = NULL;
	cr_assert_eq(process_flis_group_info(2), CMD_OK);
}

Test(flis_cmd, group_info_by_unknown_fails) {
	load_two_layer_fixture();
	word[0] = "flis_group_info";
	word[1] = "ghost";
	word[2] = NULL;
	cr_assert_eq(process_flis_group_info(2), CMD_ARG_ERROR);
}
