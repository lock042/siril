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
 * test_flis_canvas_cmd — exercises the §7 canvas-scoped commands
 * (process_flis_canvas_resize / _fit / _rotate / _mirrorx / _mirrory).
 *
 * Each command wraps a primitive that's tested in detail by
 * test_flis_canvas; this file covers the surface (argument parsing,
 * refusal paths, end-to-end mutation).
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"

cominfo com;
fits *gfit;
char *word[MAX_COMMAND_WORDS];

TestSuite(flis_canvas_cmd, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

static void make_fixture(void) {
	flis_test_add_layer(flis_test_make_mono_fits(100, 100, 0.0f), "base");
	flis_layer_t *top = flis_test_add_layer(
	    flis_test_make_mono_fits(20, 20, 0.5f), "patch");
	top->position_x = 40; top->position_y = 30;
	com.uniq->canvas_w = 100;
	com.uniq->canvas_h = 100;
	uniq_set_active_layer(com.uniq, 0);
	gfit = flis_active_layer_fit();
}

Test(flis_canvas_cmd, resize_runs) {
	make_fixture();
	word[0] = "flis_canvas_resize";
	word[1] = "-w=200";
	word[2] = "-h=150";
	word[3] = NULL;
	cr_assert_eq(process_flis_canvas_resize(3), CMD_OK);
	cr_assert_eq(com.uniq->canvas_w, 200u);
	cr_assert_eq(com.uniq->canvas_h, 150u);
}

Test(flis_canvas_cmd, resize_with_offsets) {
	make_fixture();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	word[0] = "flis_canvas_resize";
	word[1] = "-w=200";
	word[2] = "-h=150";
	word[3] = "-dx=5";
	word[4] = "-dy=10";
	word[5] = NULL;
	cr_assert_eq(process_flis_canvas_resize(5), CMD_OK);
	cr_assert_eq(base->position_x, 5);
	cr_assert_eq(base->position_y, 10);
}

Test(flis_canvas_cmd, resize_requires_dims) {
	make_fixture();
	word[0] = "flis_canvas_resize";
	word[1] = "-dx=5";
	word[2] = NULL;
	cr_assert_eq(process_flis_canvas_resize(2), CMD_ARG_ERROR);
}

Test(flis_canvas_cmd, resize_rejects_negative_dims) {
	make_fixture();
	word[0] = "flis_canvas_resize";
	word[1] = "-w=-100";
	word[2] = "-h=50";
	word[3] = NULL;
	cr_assert_eq(process_flis_canvas_resize(3), CMD_ARG_ERROR);
}

Test(flis_canvas_cmd, resize_rejects_unknown_option) {
	make_fixture();
	word[0] = "flis_canvas_resize";
	word[1] = "-bogus=1";
	word[2] = NULL;
	cr_assert_eq(process_flis_canvas_resize(2), CMD_ARG_ERROR);
}

Test(flis_canvas_cmd, fit_runs) {
	make_fixture();
	/* base 100x100 at (0,0), top 20x20 at (40,30) — bbox is 100x100. */
	word[0] = "flis_canvas_fit";
	word[1] = NULL;
	cr_assert_eq(process_flis_canvas_fit(1), CMD_OK);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 100u);
}

Test(flis_canvas_cmd, fit_all_flag) {
	make_fixture();
	flis_layer_t *top = (flis_layer_t *)com.uniq->layers->next->data;
	top->visible = FALSE;
	top->position_x = 200; top->position_y = 200;
	word[0] = "flis_canvas_fit";
	word[1] = "-all";
	word[2] = NULL;
	cr_assert_eq(process_flis_canvas_fit(2), CMD_OK);
	/* bbox includes the hidden top layer extending to (220, 220). */
	cr_assert_eq(com.uniq->canvas_w, 220u);
	cr_assert_eq(com.uniq->canvas_h, 220u);
}

Test(flis_canvas_cmd, rotate_runs) {
	make_fixture();
	word[0] = "flis_canvas_rotate";
	word[1] = "180";
	word[2] = NULL;
	cr_assert_eq(process_flis_canvas_rotate(2), CMD_OK);
	cr_assert_eq(com.uniq->canvas_w, 100u);
	cr_assert_eq(com.uniq->canvas_h, 100u);
}

Test(flis_canvas_cmd, rotate_requires_angle) {
	make_fixture();
	word[0] = "flis_canvas_rotate";
	word[1] = NULL;
	cr_assert_eq(process_flis_canvas_rotate(1), CMD_WRONG_N_ARG);
}

Test(flis_canvas_cmd, rotate_rejects_non_numeric_angle) {
	make_fixture();
	word[0] = "flis_canvas_rotate";
	word[1] = "abc";
	word[2] = NULL;
	cr_assert_eq(process_flis_canvas_rotate(2), CMD_ARG_ERROR);
}

Test(flis_canvas_cmd, mirrorx_runs) {
	make_fixture();
	flis_layer_t *top = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(process_flis_canvas_mirrorx(1), CMD_OK);
	cr_assert_eq(top->position_y, 100 - 30 - 20);  /* flipped about y-axis */
}

Test(flis_canvas_cmd, mirrory_runs) {
	make_fixture();
	flis_layer_t *top = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(process_flis_canvas_mirrory(1), CMD_OK);
	cr_assert_eq(top->position_x, 100 - 40 - 20);  /* flipped about x-axis */
}

Test(flis_canvas_cmd, all_refuse_when_no_flis) {
	word[0] = "flis_canvas_resize";
	word[1] = "-w=100"; word[2] = "-h=100"; word[3] = NULL;
	cr_assert_eq(process_flis_canvas_resize(3), CMD_GENERIC_ERROR);

	word[0] = "flis_canvas_fit"; word[1] = NULL;
	cr_assert_eq(process_flis_canvas_fit(1), CMD_GENERIC_ERROR);

	word[0] = "flis_canvas_rotate"; word[1] = "90"; word[2] = NULL;
	cr_assert_eq(process_flis_canvas_rotate(2), CMD_GENERIC_ERROR);

	word[0] = "flis_canvas_mirrorx"; word[1] = NULL;
	cr_assert_eq(process_flis_canvas_mirrorx(1), CMD_GENERIC_ERROR);

	word[0] = "flis_canvas_mirrory"; word[1] = NULL;
	cr_assert_eq(process_flis_canvas_mirrory(1), CMD_GENERIC_ERROR);
}
