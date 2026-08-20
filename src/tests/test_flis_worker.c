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
 * test_flis_worker — exercises the generic_layer_worker dispatch
 * machinery introduced at stage 1.5.  The worker is normally invoked
 * via start_in_new_thread (asynchronous, with an idle callback to do
 * the GUI refresh).  For the test we call the worker function directly
 * — runs synchronously on the test thread, makes assertions simple,
 * and skips the GUI/idle plumbing that has no analogue in a headless
 * test process.
 *
 * Two principal axes are tested:
 *   • the hook fires and its retval is propagated;
 *   • the headless undo gate per §C.1a — a hook that calls
 *     undo_save_flis_* saves nothing when com.script is TRUE.
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/op_descriptor.h"
#include "core/nde/nde_history.h"

cominfo com;
fits *gfit;

TestSuite(flis_worker, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* Trivial hook that increments a counter via args->user. */
struct counter_arg {
	destructor d;
	int count;
};
static void counter_arg_destroy(void *p) { free(p); }
static int counter_hook(struct generic_layer_args *args) {
	struct counter_arg *u = (struct counter_arg *)args->user;
	u->count++;
	return 0;
}
static int failing_hook(struct generic_layer_args *args) {
	(void)args;
	return 1;  /* simulate operation failure */
}

/* Hook that calls undo_save_flis_layer_props with the worker's target
 * layer — used to verify the §C.1a headless gate. */
static int undo_recording_hook(struct generic_layer_args *args) {
	undo_save_flis_layer_props(args->layer, "test op");
	return 0;
}

/* ---- colour-model requirements, refused before anything happens --------- */

/* A file-static counter, not args->user: the worker frees user through
 * destroy_any_args(), which reads its first field as a destructor. */
static int image_hook_runs;

static int counting_image_hook(struct generic_img_args *args, fits *fit, int nb) {
	(void)args; (void)fit; (void)nb;
	image_hook_runs++;
	return 0;
}

/* Not in op_descriptors.def, so the registry conformance tests do not see
 * these two; they exist only to isolate the flag as the single difference. */
static const op_descriptor op_desc_needs_rgb = {
	.id = "test.needs_rgb", .version = 1,
	.image_hook = counting_image_hook,
	.description = N_("Needs colour"),
	.flags = OP_REQ_RGB,
};
static const op_descriptor op_desc_needs_nothing = {
	.id = "test.needs_nothing", .version = 1,
	.image_hook = counting_image_hook,
	.description = N_("Needs nothing"),
	.flags = 0,
};

static int run_image_op(const op_descriptor *op) {
	image_hook_runs = 0;
	struct generic_img_args *args = calloc(1, sizeof(*args));
	args->fit = gfit;
	args->op = op;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;
	args->max_threads = 1;
	args->mem_ratio = -1.0f;
	gboolean prev = com.headless;
	com.headless = TRUE;
	int rv = GPOINTER_TO_INT(generic_image_worker(args));
	com.headless = prev;
	return rv;
}

/* The failure this replaces: PCC on a mono layer logged "will do nothing for
 * monochrome images", returned SUCCESS, and so earned an undo entry and a
 * history step for an operation that did nothing.  The refusal has to happen
 * before the hook, not inside it, or the worker has already committed. */
Test(flis_worker, an_rgb_only_op_is_refused_on_a_mono_image) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "mono");
	gfit = flis_active_layer_fit();

	cr_assert_neq(run_image_op(&op_desc_needs_rgb), 0,
	              "an RGB-only op must fail on a mono image");
	cr_assert_eq(image_hook_runs, 0, "the hook must not have run at all");

	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++)
		cr_assert(g_strcmp0(((nde_record *)g_ptr_array_index(snap, i))->op_id,
		                    "test.needs_rgb"),
		          "a refused op must leave no history step");
	if (snap)
		g_ptr_array_unref(snap);
}

/* The control: the same hook on the same image, with only the flag removed.
 * Without this the test above would pass for any reason the worker failed. */
Test(flis_worker, the_same_op_without_the_flag_runs_on_a_mono_image) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "mono");
	gfit = flis_active_layer_fit();

	cr_assert_eq(run_image_op(&op_desc_needs_nothing), 0);
	cr_assert_eq(image_hook_runs, 1,
	             "the flag is the only thing that stopped the other one");
}

Test(flis_worker, hook_fires_and_returns_success) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");

	struct counter_arg *u = calloc(1, sizeof(*u));
	u->d = counter_arg_destroy;

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer       = l;
	args->layer_hook  = counter_hook;
	args->user        = u;
	args->command     = TRUE;  /* synchronous-completion path; OK in headless */
	args->description = g_strdup("counter");

	gpointer rv = generic_layer_worker(args);
	cr_assert_eq(GPOINTER_TO_INT(rv), 0);
	/* args is freed by the worker in headless+command mode, so we can't
	 * dereference u afterwards.  But the hook was called and its retval
	 * (0) was passed through, which is what we care about here. */
}

Test(flis_worker, hook_failure_propagates_to_retval) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer       = l;
	args->layer_hook  = failing_hook;
	args->command     = TRUE;
	args->description = g_strdup("failing");

	gpointer rv = generic_layer_worker(args);
	cr_assert_eq(GPOINTER_TO_INT(rv), 1, "failing hook should yield retval=1");
}

/* §C.1a: with com.script = TRUE, a hook that calls undo_save_flis_*
 * saves nothing — the undo machinery short-circuits at the entry to
 * each save function.  The worker itself does not gate undo; the gate
 * lives in the undo save functions per the design. */
Test(flis_worker, undo_gated_off_when_script_TRUE) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	com.script = TRUE;  /* flis_test_init_com already sets this, but be explicit */

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer       = l;
	args->layer_hook  = undo_recording_hook;
	args->command     = TRUE;
	args->description = g_strdup("undo-attempting op");

	generic_layer_worker(args);
	cr_assert_null(com.undo_stack,
	               "headless mode: hook's undo_save should have no effect");
}

Test(flis_worker, undo_saved_when_script_FALSE) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	com.script = FALSE;  /* simulate GUI mode */

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer       = l;
	args->layer_hook  = undo_recording_hook;
	args->command     = TRUE;
	args->description = g_strdup("undo-attempting op");

	generic_layer_worker(args);
	cr_assert_eq(g_list_length(com.undo_stack), 1,
	             "GUI mode: hook's undo_save should push exactly one entry");
}

/* ---- NDE descriptor-primary capture (steps 6-8 Part C) ---- */

/* Descriptor without serialize → Tier-B capture; layer_hook succeeds. */
static int nde_test_layer_hook(struct generic_layer_args *args) {
	args->invalidate_item_id = 5;   /* pretend it targeted layer id 5 */
	return 0;
}
static gchar *nde_test_log_hook(gpointer u, log_hook_detail d) {
	(void)u; (void)d;
	return g_strdup("test summary");
}
static const op_descriptor op_desc_test_layer = {
	.id = "test.layer", .version = 3,
	.layer_hook = nde_test_layer_hook,
	.log_hook = nde_test_log_hook,
	.description = "test layer op",
	.mem_ratio = 0.0f,
	.flags = 0,
};

/* args->op set → fill_layer_args populates layer_hook/log_hook/description
 * and the worker appends a Tier-B provenance record on success. */
Test(flis_worker, descriptor_layer_op_captures_record) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->op      = &op_desc_test_layer;   /* hooks come from the descriptor */
	args->command = TRUE;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(nde_history_live_count(), 1,
	             "descriptor layer op must capture exactly one record");
	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *r = g_ptr_array_index(snap, 0);
	cr_assert_str_eq(r->op_id, "test.layer");
	cr_assert_eq(r->op_version, 3);
	cr_assert_eq(r->tier, NDE_TIER_B, "no serialize → Tier B");
	cr_assert_null(r->params);
	cr_assert_str_eq(r->summary, "test summary");
	cr_assert_eq(r->scope, NDE_SCOPE_DOCUMENT);
	cr_assert_eq(r->target_item_id, 5);
	g_ptr_array_unref(snap);
}

/* Legacy op (args->op == NULL) must NOT capture in the worker — those
 * self-record at their own sites (Part A). */
Test(flis_worker, legacy_layer_op_does_not_capture) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer      = l;
	args->layer_hook = counter_hook;
	struct counter_arg *u = calloc(1, sizeof(*u));
	u->d = counter_arg_destroy;
	args->user       = u;
	args->command    = TRUE;

	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 0);
	cr_assert_eq(nde_history_live_count(), 0,
	             "legacy op_hook path must not double-record in the worker");
}

/* Failed hook → no record. */
Test(flis_worker, descriptor_layer_op_failure_no_record) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "x");
	static const op_descriptor op_desc_fail = {
		.id = "test.fail", .version = 1,
		.layer_hook = failing_hook,
		.description = "failing op",
	};
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->op      = &op_desc_fail;
	args->command = TRUE;
	cr_assert_eq(GPOINTER_TO_INT(generic_layer_worker(args)), 1);
	cr_assert_eq(nde_history_live_count(), 0, "no record on failure");
}
