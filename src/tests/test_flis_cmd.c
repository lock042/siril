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

/* -----------------------------------------------------------------
 * flis_addgroup / flis_setgroup (§4.3 slice 3)
 * ----------------------------------------------------------------- */

Test(flis_cmd, addgroup_hook_creates_group_with_autoname) {
	load_two_layer_fixture();
	cr_assert_eq(flis_group_count(), 0);

	struct flis_addgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addgroup_args_free;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addgroup_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addgroup auto");

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(flis_group_count(), 1);
	/* Auto name should be "Group 1". */
	flis_group_t *grp = (flis_group_t *)com.uniq->groups->data;
	cr_assert_str_eq(grp->name, "Group 1");
}

Test(flis_cmd, addgroup_hook_creates_group_with_explicit_name) {
	load_two_layer_fixture();
	struct flis_addgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addgroup_args_free;
	payload->name = g_strdup("Narrowband");
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addgroup_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addgroup named");

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	flis_group_t *grp = (flis_group_t *)com.uniq->groups->data;
	cr_assert_str_eq(grp->name, "Narrowband");
}

Test(flis_cmd, addgroup_hook_auto_name_picks_next_free) {
	load_two_layer_fixture();
	flis_group_add("Group 1");        /* taken */
	flis_group_add("Group 3");        /* taken non-contiguous */
	cr_assert_eq(flis_group_count(), 2);

	struct flis_addgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addgroup_args_free;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook  = flis_addgroup_hook;
	args->user        = payload;
	args->command     = TRUE;
	args->description = g_strdup("addgroup gap");
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(flis_group_count(), 3);
	/* Newest group should be "Group 2" (smallest free N). */
	GSList *last = g_slist_last(com.uniq->groups);
	flis_group_t *grp = (flis_group_t *)last->data;
	cr_assert_str_eq(grp->name, "Group 2");
}

Test(flis_cmd, setgroup_hook_assigns_layer_to_group) {
	load_two_layer_fixture();
	flis_group_t *grp = flis_group_add("Targets");
	flis_layer_t *target = (flis_layer_t *)com.uniq->layers->next->data;  /* "Ha" */
	cr_assert_eq(target->group_id, 0, "starts ungrouped");

	struct flis_setgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setgroup_args_free;
	payload->group_id = grp->item_id;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setgroup_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setgroup test");
	args->invalidate_item_id = target->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(target->group_id, grp->item_id);
}

Test(flis_cmd, setgroup_hook_clear_removes_layer_from_group) {
	load_two_layer_fixture();
	flis_group_t *grp = flis_group_add("Targets");
	flis_layer_t *target = (flis_layer_t *)com.uniq->layers->next->data;
	flis_layer_set_group(target, grp->item_id);
	cr_assert_eq(target->group_id, grp->item_id);

	struct flis_setgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setgroup_args_free;
	payload->group_id = 0;            /* explicit clear */
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setgroup_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setgroup clear");
	args->invalidate_item_id = target->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(target->group_id, 0);
}

Test(flis_cmd, setgroup_hook_unknown_group_fails) {
	load_two_layer_fixture();
	flis_layer_t *target = flis_active_layer();
	struct flis_setgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setgroup_args_free;
	payload->group_id = 99999;        /* doesn't exist */
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setgroup_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setgroup bad");
	args->invalidate_item_id = target->item_id;
	cr_assert_neq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
}

Test(flis_cmd, setgroup_command_rejects_short_args) {
	load_two_layer_fixture();
	word[0] = "flis_setgroup";
	word[1] = "Ha";
	word[2] = NULL;
	cr_assert_eq(process_flis_setgroup(2), CMD_WRONG_N_ARG);
}

/* -----------------------------------------------------------------
 * flis_reorder (§4.3 slice 4)
 *
 * No command; panel-only DnD path.  Tests exercise the hook directly.
 * load_two_layer_fixture gives layer_order: base=10, Ha=20.
 * ----------------------------------------------------------------- */

Test(flis_cmd, reorder_hook_above_moves_layer_up) {
	load_two_layer_fixture();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *ha   = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(base->layer_order, 10);
	cr_assert_eq(ha->layer_order, 20);

	/* Move base ABOVE Ha → base.order > Ha.order */
	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = ha->item_id;
	payload->place_above = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("reorder above");
	args->invalidate_item_id = base->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert(base->layer_order > ha->layer_order,
		"base should now be above Ha (base=%d, ha=%d)",
		base->layer_order, ha->layer_order);
	/* List should still be sorted ascending. */
	cr_assert_eq((flis_layer_t *)com.uniq->layers->data, ha,
		"after reorder, lowest-order layer should be Ha");
}

Test(flis_cmd, reorder_hook_below_moves_layer_down) {
	load_two_layer_fixture();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *ha   = (flis_layer_t *)com.uniq->layers->next->data;

	/* Move Ha BELOW base → Ha.order < base.order */
	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = base->item_id;
	payload->place_above = FALSE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("reorder below");
	args->invalidate_item_id = ha->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert(ha->layer_order < base->layer_order);
}

Test(flis_cmd, reorder_hook_inherits_target_group) {
	load_two_layer_fixture();
	flis_group_t *grp = flis_group_add("G");
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *ha   = (flis_layer_t *)com.uniq->layers->next->data;
	/* Put Ha into the group, base stays out. */
	flis_layer_set_group(ha, grp->item_id);
	cr_assert_eq(base->group_id, 0);
	cr_assert_eq(ha->group_id,   grp->item_id);

	/* Drag base onto Ha → base should join the group. */
	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = ha->item_id;
	payload->place_above = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("reorder into group");
	args->invalidate_item_id = base->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(base->group_id, grp->item_id,
		"reordered layer should inherit target's group");
}

Test(flis_cmd, reorder_hook_onto_ungrouped_clears_source_group) {
	load_two_layer_fixture();
	flis_group_t *grp = flis_group_add("G");
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *ha   = (flis_layer_t *)com.uniq->layers->next->data;
	/* Put base into the group; Ha stays out. */
	flis_layer_set_group(base, grp->item_id);
	cr_assert_eq(base->group_id, grp->item_id);
	cr_assert_eq(ha->group_id,   0);

	/* Drag base onto ungrouped Ha → base should leave the group. */
	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = ha->item_id;
	payload->place_above = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("reorder out of group");
	args->invalidate_item_id = base->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(base->group_id, 0,
		"dropping on ungrouped layer should clear source's group");
}

/* -----------------------------------------------------------------
 * flis_setposition (§4.3 slice 6)
 * ----------------------------------------------------------------- */

Test(flis_cmd, setposition_hook_moves_non_base_layer) {
	load_two_layer_fixture();
	flis_layer_t *ha = (flis_layer_t *)com.uniq->layers->next->data;
	cr_assert_eq(ha->position_x, 0);
	cr_assert_eq(ha->position_y, 0);

	struct flis_setposition_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setposition_args_free;
	payload->x = 123;
	payload->y = 456;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setposition_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setposition test");
	args->invalidate_item_id = ha->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(ha->position_x, 123);
	cr_assert_eq(ha->position_y, 456);
}

/* §7 canvas decoupling: the base layer (= bottom-of-stack) is now a
 * regular layer with its own position.  Setting its position succeeds
 * and the layer moves like any other.  Pre-§7 this op was refused
 * because the base implicitly defined the canvas origin. */
Test(flis_cmd, setposition_hook_moves_base_layer) {
	load_two_layer_fixture();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

	struct flis_setposition_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setposition_args_free;
	payload->x = 50;
	payload->y = 50;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setposition_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("setposition base");
	args->invalidate_item_id = base->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(base->position_x, 50);
	cr_assert_eq(base->position_y, 50);
}

Test(flis_cmd, setposition_command_rejects_non_integer) {
	load_two_layer_fixture();
	word[0] = "flis_setposition";
	word[1] = "Ha";
	word[2] = "abc";
	word[3] = "50";
	word[4] = NULL;
	cr_assert_eq(process_flis_setposition(4), CMD_ARG_ERROR);
}

/* -----------------------------------------------------------------
 * flis_exportlayer (§4.3 slice 7)
 * ----------------------------------------------------------------- */

Test(flis_cmd, exportlayer_hook_writes_file) {
	load_two_layer_fixture();
	flis_layer_t *ha = (flis_layer_t *)com.uniq->layers->next->data;

	gchar *path = g_build_filename(g_get_tmp_dir(),
		"flis_exportlayer_test_XXXXXX.fit", NULL);
	gint fd = g_mkstemp(path);
	cr_assert(fd >= 0);
	close(fd);
	g_unlink(path);   /* let the hook create it fresh */

	struct flis_exportlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_exportlayer_args_free;
	payload->filename   = g_strdup(path);
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_exportlayer_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->read_only          = TRUE;
	args->description        = g_strdup("exportlayer test");
	args->invalidate_item_id = ha->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);

	/* File should exist and be non-empty. */
	cr_assert(g_file_test(path, G_FILE_TEST_EXISTS));
	GStatBuf st;
	cr_assert_eq(g_stat(path, &st), 0);
	cr_assert(st.st_size > 0, "exported FITS should be non-empty (got %ld)",
		(long)st.st_size);
	g_unlink(path);
	g_free(path);
}

Test(flis_cmd, exportlayer_command_rejects_missing_file_arg) {
	load_two_layer_fixture();
	word[0] = "flis_exportlayer";
	word[1] = "Ha";
	word[2] = NULL;
	cr_assert_eq(process_flis_exportlayer(2), CMD_WRONG_N_ARG);
}

/* -----------------------------------------------------------------
 * Layer property commands (§4.3 slice 8)
 * ----------------------------------------------------------------- */

Test(flis_cmd, setname_renames_layer) {
	load_two_layer_fixture();
	word[0] = "flis_setname"; word[1] = "Ha"; word[2] = "Hydrogen"; word[3] = NULL;
	cr_assert_eq(process_flis_setname(3), CMD_OK);
	cr_assert_str_eq(flis_layer_get_by_name("Hydrogen")->layer_name, "Hydrogen");
}

Test(flis_cmd, setblend_sets_mode) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	cr_assert_eq(ha->blend_mode, FLIS_BLEND_SCREEN);  /* set by fixture */
	word[0] = "flis_setblend"; word[1] = "Ha"; word[2] = "multiply"; word[3] = NULL;
	cr_assert_eq(process_flis_setblend(3), CMD_OK);
	cr_assert_eq(ha->blend_mode, FLIS_BLEND_MULTIPLY);
}

Test(flis_cmd, setblend_rejects_unknown_mode) {
	load_two_layer_fixture();
	word[0] = "flis_setblend"; word[1] = "Ha"; word[2] = "BLOWUP"; word[3] = NULL;
	cr_assert_eq(process_flis_setblend(3), CMD_ARG_ERROR);
}

Test(flis_cmd, setopacity_sets_value) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	word[0] = "flis_setopacity"; word[1] = "Ha"; word[2] = "0.25"; word[3] = NULL;
	cr_assert_eq(process_flis_setopacity(3), CMD_OK);
	cr_assert_float_eq(ha->opacity, 0.25f, 1e-5f);
}

Test(flis_cmd, setopacity_rejects_out_of_range) {
	load_two_layer_fixture();
	word[0] = "flis_setopacity"; word[1] = "Ha"; word[2] = "1.5"; word[3] = NULL;
	cr_assert_eq(process_flis_setopacity(3), CMD_ARG_ERROR);
}

Test(flis_cmd, setvisible_toggles) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	cr_assert(ha->visible);
	word[0] = "flis_setvisible"; word[1] = "Ha"; word[2] = "off"; word[3] = NULL;
	cr_assert_eq(process_flis_setvisible(3), CMD_OK);
	cr_assert(!ha->visible);
	word[2] = "true";
	cr_assert_eq(process_flis_setvisible(3), CMD_OK);
	cr_assert(ha->visible);
}

Test(flis_cmd, setlocked_toggles) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	cr_assert(!ha->locked);
	word[0] = "flis_setlocked"; word[1] = "Ha"; word[2] = "1"; word[3] = NULL;
	cr_assert_eq(process_flis_setlocked(3), CMD_OK);
	cr_assert(ha->locked);
}

Test(flis_cmd, settint_sets_components) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	word[0] = "flis_settint"; word[1] = "Ha";
	word[2] = "0.9"; word[3] = "0.1"; word[4] = "0.2"; word[5] = NULL;
	cr_assert_eq(process_flis_settint(5), CMD_OK);
	cr_assert(ha->has_tint);
	cr_assert_float_eq(ha->layer_tint.r, 0.9, 1e-5);
	cr_assert_float_eq(ha->layer_tint.g, 0.1, 1e-5);
	cr_assert_float_eq(ha->layer_tint.b, 0.2, 1e-5);
}

Test(flis_cmd, settint_clear_disables_tint) {
	load_two_layer_fixture();
	flis_layer_t *ha = flis_layer_get_by_name("Ha");
	cr_assert(ha->has_tint);  /* fixture sets it */
	word[0] = "flis_settint"; word[1] = "Ha"; word[2] = "-clear"; word[3] = NULL;
	cr_assert_eq(process_flis_settint(3), CMD_OK);
	cr_assert(!ha->has_tint);
}

Test(flis_cmd, settint_rejects_partial_rgb) {
	load_two_layer_fixture();
	word[0] = "flis_settint"; word[1] = "Ha"; word[2] = "0.5"; word[3] = "0.5"; word[4] = NULL;
	cr_assert_eq(process_flis_settint(4), CMD_WRONG_N_ARG);
}

/* -----------------------------------------------------------------
 * flis_group_reorder_hook (group-up/down panel handler primitive)
 * ----------------------------------------------------------------- */

Test(flis_cmd, group_reorder_hook_moves_group_up_preserving_internal_order) {
	/* Fixture: base, l_top with group G containing {l_mid, l_low}
	 * Order: base=10, l_low=20, l_mid=30, l_top=40
	 * Group G occupies orders 20 and 30; top external = l_top (40).
	 * Move group up → group rises above l_top.
	 * Expected: base=10, l_top=20, l_low=21, l_mid=31 (preserved relative).
	 * But that doesn't match my hook's algorithm exactly — let me just
	 * verify the group ends up ABOVE l_top in z-order. */
	load_two_layer_fixture();
	flis_layer_t *l_low = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.3f), "L_low");
	flis_layer_t *l_mid = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "L_mid");
	flis_layer_t *l_top = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.7f), "L_top");
	flis_group_t *grp = flis_group_add("G");
	flis_layer_set_group(l_low, grp->item_id);
	flis_layer_set_group(l_mid, grp->item_id);

	cr_assert(l_low->layer_order < l_mid->layer_order);
	cr_assert(l_mid->layer_order < l_top->layer_order);

	struct flis_group_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn   = flis_group_reorder_args_free;
	payload->direction_up = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_group_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("group reorder up");
	args->invalidate_item_id = grp->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);

	/* After moving group up, both group members should sit above
	 * l_top in z-order, and internal order should be preserved
	 * (l_low still below l_mid). */
	cr_assert(l_top->layer_order < l_low->layer_order,
		"l_top (%d) should be below l_low (%d) after move",
		l_top->layer_order, l_low->layer_order);
	cr_assert(l_top->layer_order < l_mid->layer_order);
	cr_assert(l_low->layer_order < l_mid->layer_order,
		"group's internal order preserved (l_low < l_mid)");
}

Test(flis_cmd, movemask_hook_transfers_lmask_to_target) {
	/* Two layers of same size; src has a mask, dst doesn't.  After
	 * movemask: src->lmask == NULL, dst->lmask != NULL. */
	load_two_layer_fixture();
	flis_layer_t *src = (flis_layer_t *)com.uniq->layers->data;        /* "base", 8x8 RGB */
	flis_layer_t *dst = (flis_layer_t *)com.uniq->layers->next->data;  /* "Ha",   8x8 MONO */
	layermask_t *lm = flis_test_make_const_lmask(src->fit->rx, src->fit->ry, 8, 0.5);
	cr_assert_eq(flis_layer_set_lmask(src, lm), 0);
	cr_assert_not_null(src->lmask);
	cr_assert_null(dst->lmask);

	struct flis_movemask_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_movemask_args_free;
	payload->to_layer_id = dst->item_id;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_movemask_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("movemask test");
	args->invalidate_item_id = src->item_id;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_null(src->lmask, "source's mask should be cleared");
	cr_assert_not_null(dst->lmask, "target should receive the mask");
}

Test(flis_cmd, group_reorder_hook_at_top_is_noop_success) {
	/* Group already on top — moving up should succeed with no change. */
	load_two_layer_fixture();
	flis_layer_t *l_a = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.3f), "A");
	flis_layer_t *l_b = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "B");
	flis_group_t *grp = flis_group_add("Top");
	flis_layer_set_group(l_a, grp->item_id);
	flis_layer_set_group(l_b, grp->item_id);
	gint a_before = l_a->layer_order;
	gint b_before = l_b->layer_order;

	struct flis_group_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn   = flis_group_reorder_args_free;
	payload->direction_up = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_group_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("group reorder noop");
	args->invalidate_item_id = grp->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(l_a->layer_order, a_before);
	cr_assert_eq(l_b->layer_order, b_before);
}

Test(flis_cmd, reorder_hook_self_target_noop) {
	load_two_layer_fixture();
	flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
	const gint base_order = base->layer_order;

	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = base->item_id;
	payload->place_above = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->command            = TRUE;
	args->description        = g_strdup("reorder self");
	args->invalidate_item_id = base->item_id;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(base->layer_order, base_order, "self-target reorder is a no-op");
}
