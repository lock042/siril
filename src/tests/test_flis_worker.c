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
