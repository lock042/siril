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
 * test_nde_tier_c — chain-validity logic for Tier-C (replayable Python
 * script) records, nde-phase5 replay half.  Exercises the re-run gate through
 * nde_chain_build: a valid recipe (existing script, matching sha256, GUI
 * mode) makes the chain replayable; a failing gate (stale script, missing
 * file, headless) degrades the record to a barrier with a restart point when
 * its output checkpoint exists, and to a hard blocker when it does not.  The
 * actual script re-execution needs a python venv and is covered by the GUI
 * e2e script, not here.
 */

#include <criterion/criterion.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <string.h>
#include "flis_test_helpers.h"
#include "core/siril.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "io/image_format_fits.h"

cominfo com;
fits *gfit;

static gchar *tmpdir;

static void tier_c_test_init(void) {
	flis_test_init_com();
	com.pref.nde_cache_mb = 2048;   /* else the checkpoint pool is silently off */
	com.headless = FALSE;           /* the gate refuses re-runs when headless */
	gfit = flis_test_make_mono_fits(8, 8, 0.25f);
	cr_assert_not_null(gfit);
	tmpdir = g_dir_make_tmp("nde_tier_c_XXXXXX", NULL);
	cr_assert_not_null(tmpdir);
}

static void tier_c_test_fini(void) {
	if (gfit) {
		clearfits(gfit);
		free(gfit);
		gfit = NULL;
	}
	g_free(tmpdir);
	tmpdir = NULL;
	flis_test_cleanup_com();
}

TestSuite(nde_tier_c, .init = tier_c_test_init, .fini = tier_c_test_fini);

/* Write a small script file and return its path (caller frees). */
static gchar *make_script(const char *name) {
	gchar *path = g_build_filename(tmpdir, name, NULL);
	cr_assert(g_file_set_contents(path, "print('replay me')\n", -1, NULL));
	return path;
}

/* Capture a Tier-C record whose recipe carries @path/@sha and a canonical
 * argv blob.  @post non-NULL stores the output checkpoint (as the real scope
 * engine always does). */
static gint64 capture_tier_c(const char *path, const char *sha, const fits *post) {
	GString *kv = nde_kv_start();
	nde_kv_add_str(kv, "script", path);
	nde_kv_add_str(kv, "sha256", sha);
	nde_kv_add_str(kv, "args", "{\"version\":1,\"argv\":[\"-npoints\",\"25\"]}");
	gchar *params = nde_kv_end(kv);
	gint64 id = nde_capture_script("python.script", NDE_SCOPE_LAYER, -1,
	                               params, "myscript.py", post);
	cr_assert(id > 0);
	return id;
}

/* ---- valid recipe: the chain replays through the script ----------------- */

Test(nde_tier_c, valid_recipe_is_replayable) {
	gchar *path = make_script("ok.py");
	gchar *sha = nde_file_sha256(path, NULL);
	cr_assert_not_null(sha);

	nde_checkpoint_baseline_ensure(gfit, -1);
	capture_tier_c(path, sha, gfit);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(chain->replayable,
	          "valid Tier-C member must leave the chain replayable (first reason: %s)",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	cr_assert_eq(chain->records->len, 1);
	nde_chain_free(chain);

	nde_history_attach(NULL);
	g_free(sha);
	g_free(path);
}

/* ---- stale script: degrade to a barrier with a restart point ------------ */

Test(nde_tier_c, stale_script_degrades_to_restartable_barrier) {
	gchar *path = make_script("stale.py");

	nde_checkpoint_baseline_ensure(gfit, -1);
	gint64 id = capture_tier_c(path, "0000000000000000", gfit);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable,
	          "a stale script cannot be re-run — full replay must be off");
	cr_assert(chain->tail_replayable,
	          "the output checkpoint must leave a restart point");
	cr_assert_eq(chain->restart_ckpt_id, id);
	cr_assert_eq(chain->reasons->len, 0,
	             "a restartable barrier is informational, not a blocker");
	nde_chain_free(chain);

	nde_history_attach(NULL);
	g_free(path);
}

/* ---- headless: the gate refuses re-execution ---------------------------- */

Test(nde_tier_c, headless_refuses_rerun) {
	gchar *path = make_script("headless.py");
	gchar *sha = nde_file_sha256(path, NULL);

	nde_checkpoint_baseline_ensure(gfit, -1);
	capture_tier_c(path, sha, gfit);

	com.headless = TRUE;
	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable,
	          "headless mode must refuse script re-execution");
	cr_assert(chain->tail_replayable, "checkpoint fallback must survive");
	nde_chain_free(chain);
	com.headless = FALSE;

	nde_history_attach(NULL);
	g_free(sha);
	g_free(path);
}

/* ---- deleted script file: degrade like stale ---------------------------- */

Test(nde_tier_c, missing_script_degrades) {
	gchar *path = make_script("gone.py");
	gchar *sha = nde_file_sha256(path, NULL);

	nde_checkpoint_baseline_ensure(gfit, -1);
	capture_tier_c(path, sha, gfit);
	cr_assert_eq(g_unlink(path), 0);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(chain->tail_replayable);
	nde_chain_free(chain);

	nde_history_attach(NULL);
	g_free(sha);
	g_free(path);
}

/* ---- failing gate with NO checkpoint: hard blocker ---------------------- */

Test(nde_tier_c, stale_without_checkpoint_blocks) {
	gchar *path = make_script("blocked.py");

	nde_checkpoint_baseline_ensure(gfit, -1);
	capture_tier_c(path, "0000000000000000", NULL);   /* no output checkpoint */

	nde_chain *chain = nde_chain_build(-1);
	cr_assert(!chain->replayable);
	cr_assert(!chain->tail_replayable,
	          "no checkpoint and no re-run: everything after is frozen");
	cr_assert_eq(chain->reasons->len, 1);
	cr_assert(strstr(g_ptr_array_index(chain->reasons, 0), "changed") != NULL,
	          "the reason must say the script changed, got: %s",
	          (char *)g_ptr_array_index(chain->reasons, 0));
	nde_chain_free(chain);

	nde_history_attach(NULL);
	g_free(path);
}

/* ---- valid Tier-C between Tier-A members: whole chain replayable -------- */

Test(nde_tier_c, mixed_chain_membership) {
	gchar *path = make_script("mid.py");
	gchar *sha = nde_file_sha256(path, NULL);

	nde_checkpoint_baseline_ensure(gfit, -1);
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("geometry.mirrorx");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("x_axis=1");
	nde_history_append(rec);
	capture_tier_c(path, sha, gfit);
	rec = nde_record_new();
	rec->op_id = g_strdup("geometry.mirrorx");
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup("x_axis=1");
	nde_history_append(rec);

	nde_chain *chain = nde_chain_build(-1);
	cr_assert_eq(chain->records->len, 3);
	cr_assert(chain->replayable,
	          "Tier-A / Tier-C / Tier-A must be fully replayable (first reason: %s)",
	          chain->reasons->len ? (char *)g_ptr_array_index(chain->reasons, 0) : "none");
	nde_chain_free(chain);

	nde_history_attach(NULL);
	g_free(sha);
	g_free(path);
}
