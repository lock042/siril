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
 * test_nde_conductor — the off-worker NDE replay conductor's reservation and
 * identity gate (phase 5, C-ops foundation).  These are MULTITHREADED tests:
 * they drive the real processing worker and the SLOT_REPLAY primitives with
 * concurrent threads to exercise the points where the gate could go wrong,
 * rather than confirming via the synchronous _execute cores (which never touch
 * the conductor thread or the reservation).
 *
 * The adversarial case is `admitted_script_not_blocked_by_stranger`: while the
 * conductor holds the slot, a resident helper script's cmd() (a "stranger")
 * blocks in the submit gate; the replayed script's own cmd() must still be
 * admitted and run.  Because start_in_new_thread claims job_active_flag before
 * the gate, a naive design lets the parked stranger's held flag refuse the
 * replay's own submission — the deadlock-class hole this test guards.
 */

#include <criterion/criterion.h>
#include <glib.h>
#include "flis_test_helpers.h"
#include "core/siril.h"
#include "core/processing.h"
#include "core/processing_thread.h"

cominfo com;
fits *gfit;

static void conductor_test_init(void) {
	flis_test_init_com();
	com.headless = TRUE;          /* clean worker join at shutdown */
	processing_system_init();
}

static void conductor_test_fini(void) {
	processing_system_shutdown();
	flis_test_cleanup_com();
}

TestSuite(nde_conductor, .init = conductor_test_init, .fini = conductor_test_fini);

/* A trivial worker job: mark that it ran.  Runs on the processing worker. */
static gpointer mark_job(gpointer p) {
	g_atomic_int_set((gint *)p, 1);
	return GINT_TO_POINTER(0);
}

struct submitter {
	gint ran;              /* set by the job when it executes */
	gboolean submit_ret;   /* return value of start_in_new_thread */
	gboolean register_self;/* register as the replay script before submitting */
};

static gpointer submitter_fn(gpointer p) {
	struct submitter *s = p;
	if (s->register_self)
		processing_register_replay_script(g_thread_self());
	/* com.script is TRUE in tests, so this blocks until the job completes
	 * (admitted) or returns FALSE immediately (refused). */
	s->submit_ret = start_in_new_thread(mark_job, &s->ran);
	return NULL;
}

/* Spin until the slot shows a job claimed, or give up after ~2s. */
static gboolean wait_job_active(void) {
	for (int i = 0; i < 10000; i++) {
		if (processing_is_job_active())
			return TRUE;
		g_usleep(200);
	}
	return FALSE;
}

/* ---- reservation basics ------------------------------------------------- */

Test(nde_conductor, reserve_clears_cancel_and_is_exclusive) {
	/* A stale cancel from a prior aborted op must not survive into a replay:
	 * the off-worker conductor never passes the worker's per-job reset. */
	processing_request_cancel();
	cr_assert(!processing_should_continue(), "precondition: cancel set");

	cr_assert(replay_reserve_slot(), "reserve should succeed when free");
	cr_assert(processing_should_continue(), "reserve must clear cancel_flag");
	cr_assert(processing_is_reserved_for_replay());
	cr_assert(!processing_is_reserved_for_python(),
	          "replay reservation must not read as a python reservation");

	cr_assert(!replay_reserve_slot(), "reserve must be exclusive");

	replay_release_slot();
	cr_assert(!processing_is_reserved_for_replay(), "release must free the slot");
}

Test(nde_conductor, reserve_refused_while_a_job_is_active) {
	cr_assert(reserve_thread(), "reserve_thread should claim the idle slot");
	cr_assert(!replay_reserve_slot(),
	          "replay reserve must refuse while a job/reservation is active");
	unreserve_thread();
	cr_assert(replay_reserve_slot(), "reserve should succeed once free again");
	replay_release_slot();
}

/* ---- claim/release no-op semantics for the replayed script -------------- */

Test(nde_conductor, script_claim_and_release_are_noops_under_replay) {
	cr_assert(replay_reserve_slot());
	replay_bind_conductor();
	processing_register_replay_script(g_thread_self());  /* we are "the script" */

	/* The conductor already holds exclusivity, so the script's claim is a
	 * no-op success and must NOT flip the owner to python. */
	cr_assert_eq(claim_thread_for_python(), 0, "script claim should no-op succeed");
	cr_assert(processing_is_reserved_for_replay(), "owner must stay REPLAY after claim");
	cr_assert(!processing_is_reserved_for_python());

	/* And its matching release must NOT tear down the replay reservation. */
	python_releases_thread();
	cr_assert(processing_is_reserved_for_replay(),
	          "script release must leave SLOT_REPLAY held by the conductor");

	processing_unregister_replay_script();
	replay_release_slot();
}

static gpointer stranger_claim_fn(gpointer p) {
	*(int *)p = claim_thread_for_python();
	return NULL;
}

Test(nde_conductor, stranger_claim_refused_under_replay) {
	cr_assert(replay_reserve_slot());
	replay_bind_conductor();
	/* No script registered: a claim from any other thread is a stranger. */
	int ret = -99;
	GThread *t = g_thread_new("stranger-claim", stranger_claim_fn, &ret);
	g_thread_join(t);
	cr_assert_eq(ret, 1, "a stranger's claim must be refused (busy) under replay");
	replay_release_slot();
}

/* ---- the adversarial gate case ----------------------------------------- */

Test(nde_conductor, admitted_script_not_blocked_by_stranger) {
	struct submitter stranger = { .register_self = FALSE };
	struct submitter script   = { .register_self = TRUE  };

	cr_assert(replay_reserve_slot());
	replay_bind_conductor();

	/* A resident helper's cmd() enters the gate and parks (it is not the
	 * replay's script), claiming job_active_flag on the way in. */
	GThread *xt = g_thread_new("stranger", submitter_fn, &stranger);
	cr_assert(wait_job_active(), "stranger did not reach the gate");
	g_usleep(30000);   /* let the stranger settle into the cond_wait */
	cr_assert(!g_atomic_int_get(&stranger.ran), "stranger must not run while parked");

	/* The replayed script's own cmd() must be admitted and run despite the
	 * stranger holding job_active_flag. */
	GThread *st = g_thread_new("script", submitter_fn, &script);
	g_thread_join(st);
	cr_assert(script.submit_ret,
	          "replay script submission was refused while a stranger was parked");
	cr_assert(g_atomic_int_get(&script.ran),
	          "replay script job did not run");

	/* The stranger is still parked — admitted only after the slot is freed. */
	cr_assert(!g_atomic_int_get(&stranger.ran), "stranger ran before release");

	replay_release_slot();
	g_thread_join(xt);
	cr_assert(g_atomic_int_get(&stranger.ran), "stranger did not run after release");
	cr_assert(stranger.submit_ret, "stranger submission ultimately failed");
}

/* A GUI submission (represented by the main thread submitting while reserved)
 * must be refused outright, never queued. */
Test(nde_conductor, background_stranger_waits_then_runs_after_release) {
	struct submitter stranger = { .register_self = FALSE };

	cr_assert(replay_reserve_slot());
	replay_bind_conductor();

	GThread *xt = g_thread_new("stranger", submitter_fn, &stranger);
	cr_assert(wait_job_active(), "stranger did not reach the gate");
	g_usleep(30000);
	cr_assert(!g_atomic_int_get(&stranger.ran), "stranger must wait while reserved");

	replay_release_slot();
	g_thread_join(xt);
	cr_assert(g_atomic_int_get(&stranger.ran), "stranger must run once released");
}

/* ---- job origin attribution --------------------------------------------- */

/* The scope-suppression gates in the generic workers must be able to tell a
 * script's own job from a user's: a resident script's open provenance scope
 * swallows only script-submitted ops.  The origin rides on the job, captured
 * from the submitting thread's python-comm mark. */

static gpointer origin_probe_job(gpointer p) {
	*(gboolean *)p = processing_current_job_from_python();
	return GINT_TO_POINTER(0);
}

struct origin_probe {
	gboolean mark_python;   /* mark the submitting thread as a comm thread */
	gboolean from_python;   /* what the job observed */
	gboolean submit_ret;
};

static gpointer origin_submitter_fn(gpointer p) {
	struct origin_probe *o = p;
	if (o->mark_python)
		processing_mark_python_comm_thread();
	o->from_python = TRUE;  /* poisoned: the job must overwrite it */
	o->submit_ret = start_in_new_thread(origin_probe_job, &o->from_python);
	return NULL;
}

Test(nde_conductor, job_origin_tracks_submitting_thread) {
	/* A job submitted from a python comm thread reads as from-python... */
	struct origin_probe script = { .mark_python = TRUE };
	GThread *st = g_thread_new("comm-thread", origin_submitter_fn, &script);
	g_thread_join(st);
	cr_assert(script.submit_ret, "script-side submission failed");
	cr_assert(script.from_python,
	          "a comm-thread submission must attribute its job to the script");

	/* ...and one submitted from any other thread does not. */
	struct origin_probe user = { .mark_python = FALSE };
	GThread *ut = g_thread_new("user-thread", origin_submitter_fn, &user);
	g_thread_join(ut);
	cr_assert(user.submit_ret, "user-side submission failed");
	cr_assert(!user.from_python,
	          "a non-comm-thread submission must not read as from-python");
}

Test(nde_conductor, job_origin_meaningless_outside_a_job) {
	/* Outside the worker there is no current job: always FALSE, even on a
	 * thread that marked itself as a comm thread. */
	cr_assert(!processing_current_job_from_python());
}
