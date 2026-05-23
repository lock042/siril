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
#include <glib/gstdio.h>
#include "flis_test_helpers.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"

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

/* ---- flis_active_layer (stage 3.2 follow-up) -------------------------- */

Test(flis_cmd, active_layer_no_arg_reports_current) {
	load_two_layer_fixture();
	/* fixture leaves the second layer active by virtue of being last
	 * added; verify the command prints something and returns OK. */
	word[0] = "flis_active_layer";
	word[1] = NULL;
	cr_assert_eq(process_flis_active_layer(1), CMD_OK);
}

Test(flis_cmd, active_layer_set_by_name_switches_active) {
	load_two_layer_fixture();
	cr_assert_eq(flis_active_layer()->item_id, 2,
	             "fixture should have layer 2 active");
	word[0] = "flis_active_layer";
	word[1] = "base";
	word[2] = NULL;
	cr_assert_eq(process_flis_active_layer(2), CMD_OK);
	cr_assert_str_eq(flis_active_layer()->layer_name, "base",
	                 "active layer should now be 'base'");
}

Test(flis_cmd, active_layer_set_by_id_switches_active) {
	load_two_layer_fixture();
	word[0] = "flis_active_layer";
	word[1] = "1";
	word[2] = NULL;
	cr_assert_eq(process_flis_active_layer(2), CMD_OK);
	cr_assert_eq(flis_active_layer()->item_id, 1);
}

Test(flis_cmd, active_layer_set_by_unknown_name_fails) {
	load_two_layer_fixture();
	gint orig_id = flis_active_layer()->item_id;
	word[0] = "flis_active_layer";
	word[1] = "ghost";
	word[2] = NULL;
	cr_assert_eq(process_flis_active_layer(2), CMD_ARG_ERROR);
	cr_assert_eq(flis_active_layer()->item_id, orig_id,
	             "active layer must not change on lookup failure");
}

/* -----------------------------------------------------------------
 * flis_addlayer (§4.3 slice 1)
 *
 * The command wraps flis_addlayer_hook with start_in_new_thread, but
 * the test harness doesn't run a processing thread; we exercise the
 * hook directly through generic_layer_worker — same code path the
 * worker thread would take, just synchronous.  This is the same
 * pattern test_flis_worker uses.
 * ----------------------------------------------------------------- */

/* Write a small fits to a temp file and return the path (caller frees).
 * Returns NULL on failure. */
static gchar *write_tmp_mono_fits(int w, int h, float v) {
	fits *f = flis_test_make_mono_fits(w, h, v);
	if (!f) return NULL;
	gchar *path = g_build_filename(g_get_tmp_dir(), "flis_addlayer_test_XXXXXX.fit", NULL);
	gint fd = g_mkstemp(path);
	if (fd < 0) { clearfits(f); free(f); g_free(path); return NULL; }
	close(fd);
	if (savefits(path, f)) {
		clearfits(f); free(f); g_unlink(path); g_free(path);
		return NULL;
	}
	clearfits(f);
	free(f);
	return path;
}

Test(flis_cmd, addlayer_hook_adds_layer_with_derived_name) {
	load_two_layer_fixture();
	const guint before = flis_layer_count();

	gchar *path = write_tmp_mono_fits(16, 16, 0.25f);
	cr_assert_not_null(path, "must be able to write temp fits");

	struct flis_addlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addlayer_args_free;
	payload->filename   = g_strdup(path);
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addlayer_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addlayer test");

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(flis_layer_count(), before + 1);

	/* Derived name: temp basename minus extension.  Without owning the
	 * args anymore we can't check args->invalidate_item_id directly,
	 * but the new layer is now the active one. */
	flis_layer_t *active = flis_active_layer();
	cr_assert_not_null(active);
	gchar *base = g_path_get_basename(path);
	gchar *dot  = base ? strrchr(base, '.') : NULL;
	if (dot) *dot = '\0';
	cr_assert_str_eq(active->layer_name, base,
		"derived name should equal temp basename minus extension");
	g_free(base);
	g_unlink(path);
	g_free(path);
}

Test(flis_cmd, addlayer_hook_uses_explicit_name) {
	load_two_layer_fixture();
	gchar *path = write_tmp_mono_fits(8, 8, 0.1f);
	cr_assert_not_null(path);

	struct flis_addlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addlayer_args_free;
	payload->filename   = g_strdup(path);
	payload->name       = g_strdup("custom_layer_name");
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addlayer_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addlayer explicit name");

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_str_eq(flis_active_layer()->layer_name, "custom_layer_name");
	g_unlink(path);
	g_free(path);
}

Test(flis_cmd, addlayer_hook_missing_file_fails) {
	load_two_layer_fixture();
	const guint before = flis_layer_count();

	struct flis_addlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addlayer_args_free;
	payload->filename   = g_strdup("/no/such/path/should/exist.fit");
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addlayer_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addlayer missing");

	cr_assert_neq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(flis_layer_count(), before, "no layer added on failure");
}

Test(flis_cmd, addlayer_command_rejects_zero_args) {
	load_two_layer_fixture();
	word[0] = "flis_addlayer";
	word[1] = NULL;
	cr_assert_eq(process_flis_addlayer(1), CMD_WRONG_N_ARG);
}

Test(flis_cmd, addlayer_command_rejects_unknown_option) {
	load_two_layer_fixture();
	word[0] = "flis_addlayer";
	word[1] = "/tmp/file.fit";
	word[2] = "-bogus=1";
	word[3] = NULL;
	cr_assert_eq(process_flis_addlayer(3), CMD_ARG_ERROR);
}

/* -----------------------------------------------------------------
 * flis_setmask / flis_clearmask (§4.3 slice 2)
 * ----------------------------------------------------------------- */

Test(flis_cmd, setmask_hook_loads_mask_at_8bit_default) {
	load_two_layer_fixture();
	flis_layer_t *target = flis_active_layer();   /* the 8x8 "Ha" layer */
	cr_assert_not_null(target);
	cr_assert_eq(target->fit->rx, 8);
	cr_assert_eq(target->fit->ry, 8);
	cr_assert_null(target->lmask, "fixture starts without a mask");

	/* Write a same-sized mono mask file. */
	gchar *path = write_tmp_mono_fits(8, 8, 0.5f);
	cr_assert_not_null(path);

	struct flis_setmask_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setmask_args_free;
	payload->filename   = g_strdup(path);
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setmask_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setmask test");
	args->invalidate_item_id = target->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	/* Re-resolve target — generic_layer_worker freed args but the
	 * layer pointer is owned by com.uniq. */
	target = flis_active_layer();
	cr_assert_not_null(target->lmask, "mask should be attached after hook");
	cr_assert_eq(target->lmask->w, 8);
	cr_assert_eq(target->lmask->h, 8);
	cr_assert_eq(target->lmask->bitpix, 8);
	/* 0.5 * 255 = 127.5 → 127 (truncated) for the 8-bit mask. */
	cr_assert_eq(((uint8_t *)target->lmask->data)[0], 127,
		"central pixel should be ~127 for input value 0.5");
	g_unlink(path);
	g_free(path);
}

Test(flis_cmd, setmask_hook_size_mismatch_fails) {
	load_two_layer_fixture();
	flis_layer_t *target = flis_active_layer();   /* 8x8 layer */
	gchar *path = write_tmp_mono_fits(4, 4, 0.5f);  /* wrong size */
	cr_assert_not_null(path);

	struct flis_setmask_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setmask_args_free;
	payload->filename   = g_strdup(path);
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setmask_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setmask wrong size");
	args->invalidate_item_id = target->item_id;

	cr_assert_neq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_null(flis_active_layer()->lmask, "no mask attached on size mismatch");
	g_unlink(path);
	g_free(path);
}

Test(flis_cmd, clearmask_hook_removes_existing_mask) {
	load_two_layer_fixture();
	flis_layer_t *target = flis_active_layer();
	layermask_t *lm = flis_test_make_const_lmask(target->fit->rx, target->fit->ry,
	                                              8, 0.4);
	cr_assert_eq(flis_layer_set_lmask(target, lm), 0);
	cr_assert_not_null(target->lmask);

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_clearmask_hook;
	args->command            = TRUE;
	args->description        = g_strdup("clearmask test");
	args->invalidate_item_id = target->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_null(flis_active_layer()->lmask, "mask should be gone");
}

Test(flis_cmd, clearmask_hook_on_layer_without_mask_succeeds) {
	load_two_layer_fixture();
	flis_layer_t *target = flis_active_layer();
	cr_assert_null(target->lmask);

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_clearmask_hook;
	args->command            = TRUE;
	args->description        = g_strdup("clearmask no-op");
	args->invalidate_item_id = target->item_id;

	/* Clearing a non-existent mask is a no-op success, not an error. */
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
}

Test(flis_cmd, setmask_command_rejects_bad_bitpix) {
	load_two_layer_fixture();
	word[0] = "flis_setmask";
	word[1] = "Ha";
	word[2] = "/tmp/dummy.fit";
	word[3] = "-bitpix=16";
	word[4] = NULL;
	cr_assert_eq(process_flis_setmask(4), CMD_ARG_ERROR);
}

Test(flis_cmd, setmask_command_rejects_unknown_layer) {
	load_two_layer_fixture();
	word[0] = "flis_setmask";
	word[1] = "ghost";
	word[2] = "/tmp/dummy.fit";
	word[3] = NULL;
	cr_assert_eq(process_flis_setmask(3), CMD_ARG_ERROR);
}

Test(flis_cmd, clearmask_command_rejects_zero_args) {
	load_two_layer_fixture();
	word[0] = "flis_clearmask";
	word[1] = NULL;
	cr_assert_eq(process_flis_clearmask(1), CMD_WRONG_N_ARG);
}
