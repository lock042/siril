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
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * NDE replay engine — see nde_replay.h for the contract and
 * nde-phase2-3-plan.md P2.D for the design.
 */

#include <math.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "core/op_descriptor.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_snapstore.h"
#include "core/nde_replay.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/siril_pythonmodule.h"
#include "yyjson.h"

/* Destructor-first convention shared by every op params struct (the same
 * contract free_generic_img_args relies on). */
typedef void (*user_destructor)(void *);
static void destroy_user(gpointer user) {
	if (!user)
		return;
	user_destructor d = *(user_destructor *)user;
	if (d)
		d(user);
	else
		free(user);
}

/* C4: TRUE while an amend preview holds pre-K pixels in a target fits —
 * every history edit must refuse then (a commit would swap pixels the
 * preview's restore path is about to overwrite).  Defined with the amend
 * preview machinery at the end of this file. */
static gboolean amend_preview_installed(void);

static void add_reason(nde_chain *chain, const char *fmt, ...) G_GNUC_PRINTF(2, 3);
static void add_reason(nde_chain *chain, const char *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	g_ptr_array_add(chain->reasons, g_strdup_vprintf(fmt, ap));
	va_end(ap);
}

/* ---- Tier-C: replayable Python scripts (nde-phase5) -------------------- */

/* Gate + validity check for re-running a Tier-C record's script.  Returns the
 * first reason the re-run is not possible (heap string), or NULL when it is.
 * The sha256 comparison is the safety core: a matching hash means we re-execute
 * exactly the bytes the user already ran interactively when the record was
 * captured.  A failing gate is NOT a hard blocker — Tier-C records are always
 * output-checkpointed, so the chain build degrades them to a barrier with a
 * restart point (stale, but editable around).
 * TODO(securescripts): once the script sandbox lands, route the re-run through
 * it and drop the headless refusal in favour of sandbox policy. */
static gchar *tier_c_invalid_reason(const nde_record *rec) {
	if (rec->mask_active)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) ran with an active mask — mask replay is not supported yet"),
		                       rec->record_id, rec->op_id);
	GHashTable *kv = rec->params ? nde_kv_parse(rec->params) : NULL;
	if (!kv)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT ": script re-run recipe failed to parse"),
		                       rec->record_id);
	const gchar *script = g_hash_table_lookup(kv, "script");
	const gchar *sha = g_hash_table_lookup(kv, "sha256");
	gchar *reason = NULL;
	if (!script || !*script || !sha || !*sha) {
		reason = g_strdup_printf(_("record %" G_GINT64_FORMAT ": script re-run recipe is incomplete"),
		                         rec->record_id);
	} else if (com.headless) {
		/* Fail closed: automatic script re-execution needs a consenting user
		 * (or, later, the securescripts sandbox) — refuse headless. */
		reason = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): script re-execution is disabled in headless mode"),
		                         rec->record_id, rec->summary ? rec->summary : "?");
	} else if (!g_file_test(script, G_FILE_TEST_IS_REGULAR)) {
		reason = g_strdup_printf(_("record %" G_GINT64_FORMAT ": script no longer exists: %s"),
		                         rec->record_id, script);
	} else {
		gchar *cur = nde_file_sha256(script, NULL);
		if (!cur || g_strcmp0(cur, sha) != 0)
			reason = g_strdup_printf(_("record %" G_GINT64_FORMAT ": script has changed since this step was recorded: %s"),
			                         rec->record_id, script);
		g_free(cur);
	}
	g_hash_table_unref(kv);
	return reason;
}

/* Parse the recorded {"version":1,"argv":[...]} JSON into a NULL-terminated
 * argv vector (g_strfreev).  Absent/empty args yield an empty vector, not an
 * error; a malformed blob yields NULL with @err set. */
static gchar **tier_c_parse_argv(const gchar *args_json, gchar **err) {
	if (!args_json || !*args_json)
		return g_new0(gchar *, 1);
	yyjson_doc *doc = yyjson_read(args_json, strlen(args_json), 0);
	yyjson_val *root = doc ? yyjson_doc_get_root(doc) : NULL;
	yyjson_val *arr = yyjson_is_obj(root) ? yyjson_obj_get(root, "argv") : NULL;
	if (!yyjson_is_arr(arr)) {
		if (doc)
			yyjson_doc_free(doc);
		*err = g_strdup(_("recorded script arguments failed to parse"));
		return NULL;
	}
	GPtrArray *out = g_ptr_array_new();
	yyjson_val *v;
	yyjson_arr_iter it = yyjson_arr_iter_with(arr);
	while ((v = yyjson_arr_iter_next(&it)) != NULL) {
		if (!yyjson_is_str(v)) {
			g_ptr_array_set_free_func(out, g_free);
			g_ptr_array_unref(out);
			yyjson_doc_free(doc);
			*err = g_strdup(_("recorded script arguments failed to parse"));
			return NULL;
		}
		g_ptr_array_add(out, g_strdup(yyjson_get_str(v)));
	}
	yyjson_doc_free(doc);
	g_ptr_array_add(out, NULL);
	return (gchar **)g_ptr_array_free(out, FALSE);
}

/* Re-run a Tier-C record's script with its recorded arguments, transforming
 * @scratch in place (approach B, nde-phase5-plan).  The script can only see
 * the document image, so the replay state is presented AS gfit for the run:
 * swap scratch's pixels into the gfit struct, run the script synchronously on
 * the conductor (for_replay — no capture scope, comm thread admitted under
 * SLOT_REPLAY, CLI mode so no GUI pops up), then swap back, leaving the
 * script's output in @scratch and the user's image back in gfit.  Capture,
 * undo and per-op provenance are all suppressed during the run by the
 * SLOT_REPLAY guards.  Returns 0 on success. */
static int tier_c_rerun(fits *scratch, const nde_record *rec, gchar **err) {
	gchar *invalid = tier_c_invalid_reason(rec);
	if (invalid) {
		*err = invalid;
		return 1;
	}
	GHashTable *kv = nde_kv_parse(rec->params);
	const gchar *script = g_hash_table_lookup(kv, "script");
	const gchar *args_json = g_hash_table_lookup(kv, "args");
	gchar **argv = tier_c_parse_argv(args_json, err);
	if (!argv) {
		g_hash_table_unref(kv);
		return 1;
	}
	siril_log_message(_("Replaying script step: %s\n"),
	                  rec->summary ? rec->summary : script);

	g_rw_lock_writer_lock(&gfit->rwlock);
	fits_swap_all_except_rwlock(gfit, scratch);
	g_rw_lock_writer_unlock(&gfit->rwlock);

	int rc = execute_python_script(g_strdup(script),
	                               TRUE,   /* from_file */
	                               TRUE,   /* sync */
	                               argv,
	                               FALSE,  /* is_temp_file */
	                               TRUE,   /* from_cli: headless script path */
	                               FALSE,  /* debug_mode */
	                               TRUE);  /* for_replay */

	g_rw_lock_writer_lock(&gfit->rwlock);
	fits_swap_all_except_rwlock(gfit, scratch);
	g_rw_lock_writer_unlock(&gfit->rwlock);

	g_strfreev(argv);
	g_hash_table_unref(kv);
	if (rc) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): script re-run failed"),
		                       rec->record_id,
		                       rec->summary ? rec->summary : "?");
		return 1;
	}
	return 0;
}

/* First reason @rec cannot be replayed (heap string, caller frees), or NULL
 * when it replays.  The chain build decides whether an invalid member is a
 * hard blocker (no output checkpoint) or a barrier restart point. */
static gchar *member_invalid_reason(const nde_record *rec) {
	if (rec->tier == NDE_TIER_C)
		return tier_c_invalid_reason(rec);
	if (rec->tier != NDE_TIER_A)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque — not replayable"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
	if (rec->mask_active)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) ran with an active mask — mask replay is not supported yet"),
		                       rec->record_id, rec->op_id);
	const op_descriptor *op = op_descriptor_by_id(rec->op_id);
	if (!op)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT ": unknown operation '%s'"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
	if (!op->deserialize)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) has no parameter deserializer"),
		                       rec->record_id, rec->op_id);
	if (rec->op_version > op->version)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) was written by a newer version (v%d > v%d)"),
		                       rec->record_id, rec->op_id, rec->op_version, op->version);
	gpointer user = op->deserialize(rec->params, rec->op_version);
	if (!user)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		                       rec->record_id, rec->op_id);
	destroy_user(user);
	return NULL;
}

/* Compositing-state records (opacity/blend/visibility commits) are not
 * pixel operations: the compositor applies them from live FLIS state, so
 * they are neither chain members nor blockers. */
static gboolean is_compositing_state_op(const char *op_id) {
	return op_id && (!g_strcmp0(op_id, "layer.set_opacity") ||
	                 !g_strcmp0(op_id, "layer.set_blend") ||
	                 !g_strcmp0(op_id, "layer.set_visible"));
}

/* POLICY predicates — see nde_replay.h.  Liveness and trial-chain
 * replayability are the execute path's job, not these. */
gboolean nde_record_amendable(const nde_record *rec) {
	if (!rec || rec->tier != NDE_TIER_A)
		return FALSE;
	const op_descriptor *op = op_descriptor_by_id(rec->op_id);
	return op && op->deserialize;
}

gboolean nde_record_deletable(const nde_record *rec) {
	if (!rec)
		return FALSE;
	if (rec->scope == NDE_SCOPE_DOCUMENT)
		return FALSE;
	if (is_compositing_state_op(rec->op_id))
		return FALSE;
	return TRUE;
}

/* @exclude_record_id != 0 builds the TRIAL chain for a pending delete: the
 * excluded record is neither a member nor a blocker — a delete never needs
 * to replay the deleted step, only the survivors (so deleting the single
 * opaque record in a chain regains editability of everything around it). */
static nde_chain *chain_build_excluding(gint item_id, gint64 exclude_record_id) {
	nde_chain *chain = g_new0(nde_chain, 1);
	chain->item_id = item_id;
	chain->records = g_ptr_array_new();
	chain->reasons = g_ptr_array_new_with_free_func(g_free);
	chain->snapshot = nde_history_snapshot(NULL);
	chain->member_flags = g_array_new(FALSE, TRUE, sizeof(guint8));
	gboolean is_flis = is_current_image_flis();
	/* Barrier tracking (phase 4): tail_possible follows the LAST freeze
	 * cause in document order — TRUE when it left a restart point. */
	gboolean tail_possible = TRUE;

	for (guint i = 0; chain->snapshot && i < chain->snapshot->len; i++) {
		nde_record *rec = g_ptr_array_index(chain->snapshot, i);
		gboolean member = FALSE;
		if (exclude_record_id && rec->record_id == exclude_record_id)
			continue;
		if (is_compositing_state_op(rec->op_id))
			continue;
		switch (rec->scope) {
		case NDE_SCOPE_LAYER:
			member = (rec->target_item_id == item_id);
			break;
		case NDE_SCOPE_CANVAS:
			if (!is_flis) {
				/* Plain image: the whole image is the "layer"; geometry
				 * records replay as ordinary pixel ops (no canvas). */
				member = TRUE;
			} else if (rec->target_item_id == item_id) {
				/* Geometry on this layer of a FLIS document: the pixel op
				 * itself is Tier A, but its layer-offset/canvas side
				 * effects cannot be reproduced or verified on a scratch
				 * fits — phase-2 blocker (plan P2.D). */
				add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) changes layer geometry on a layered image — not replayable yet"),
				           rec->record_id, rec->op_id ? rec->op_id : "?");
			}
			/* Canvas records targeting other layers (or the canvas itself,
			 * target -1) move positions, not this item's pixels: ignore. */
			break;
		case NDE_SCOPE_DOCUMENT: {
			gboolean destructive = !g_strcmp0(rec->op_id, "layer.merge_down") ||
			                       !g_strcmp0(rec->op_id, "document.flatten");
			gboolean structural = destructive ||
			                      !g_strcmp0(rec->op_id, "layer.add") ||
			                      !g_strcmp0(rec->op_id, "layer.duplicate") ||
			                      !g_strcmp0(rec->op_id, "layer.remove") ||
			                      !g_strcmp0(rec->op_id, "layer.reorder");
			if (destructive && rec->target_item_id == item_id) {
				/* The item's pixels were destructively replaced by a
				 * composite of other layers — nothing to replay from. */
				add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) replaced this layer's pixels with a composite — not replayable"),
				           rec->record_id, rec->op_id);
				chain->tail_start = chain->records->len;
				tail_possible = FALSE;
			} else if (!structural) {
				/* FAIL CLOSED: a non-structural DOCUMENT-scope record
				 * mutated pixels document-wide (icc.convert via the layer
				 * worker today; unknown ops from newer builds tomorrow).
				 * Every layer's chain spans it, so no layer chain that
				 * contains records on both sides can replay without it. */
				add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) applies to the whole document — not replayable"),
				           rec->record_id, rec->op_id ? rec->op_id : "?");
				chain->tail_start = chain->records->len;
				tail_possible = FALSE;
			}
			/* Known structural records not targeting this item: layer.add
			 * is embodied in the baseline; the rest don't touch this
			 * item's pixels — ignore. */
			break;
		}
		default:
			break;
		}
		if (member) {
			gchar *invalid = member_invalid_reason(rec);
			guint8 flags = 0;
			if (invalid) {
				flags = NDE_CHAIN_MEMBER_BARRIER;
				if (nde_checkpoint_output_exists(rec->record_id)) {
					/* Barrier WITH a restart point: everything after it
					 * stays editable — informational, not a blocker. */
					g_free(invalid);
					tail_possible = TRUE;
					chain->restart_ckpt_id = rec->record_id;
				} else {
					/* Checkpoint-less barrier (pre-phase-4 capture): hard
					 * blocker for itself and everything after it. */
					g_ptr_array_add(chain->reasons, invalid);
					tail_possible = FALSE;
				}
			}
			g_ptr_array_add(chain->records, rec);
			g_array_append_val(chain->member_flags, flags);
			if (flags)
				chain->tail_start = chain->records->len;
		}
	}

	gboolean baseline_missing = chain->records->len &&
	                            !nde_checkpoint_baseline_exists(item_id);
	if (baseline_missing)
		add_reason(chain, _("no baseline checkpoint — the file predates baselines, or the history began before this build"));

	/* Full-chain verdict: no freeze cause anywhere and the baseline exists. */
	chain->replayable = (chain->reasons->len == 0) && chain->tail_start == 0;
	/* Tail verdict: with no freeze cause the tail IS the whole chain; with
	 * one, the last freeze cause must have left a restart checkpoint (the
	 * tail members themselves are valid by construction — an invalid member
	 * always moves tail_start past itself). */
	if (chain->tail_start == 0)
		chain->tail_replayable = chain->replayable;
	else
		chain->tail_replayable = tail_possible;
	return chain;
}

nde_chain *nde_chain_build(gint item_id) {
	return chain_build_excluding(item_id, 0);
}

void nde_chain_free(nde_chain *chain) {
	if (!chain)
		return;
	g_ptr_array_unref(chain->records);
	if (chain->snapshot)
		g_ptr_array_unref(chain->snapshot);
	g_ptr_array_unref(chain->reasons);
	g_array_unref(chain->member_flags);
	g_free(chain);
}

/* Apply members [from..upto) to @scratch (consumed on failure).  Returns
 * @scratch on success, NULL + @err on failure. */
static fits *replay_apply_records(fits *scratch, const nde_chain *chain,
                                  guint from, guint upto, gchar **err) {
	for (guint i = from; i < upto; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (!processing_should_continue()) {
			*err = g_strdup(_("cancelled"));
			goto fail;
		}
		if (rec->tier == NDE_TIER_C) {
			/* Replayable Python script: re-run it on the accumulated state
			 * (nde-phase5).  Chain membership already vetted the recipe via
			 * member_invalid_reason, but the gate re-checks — the file can
			 * change between chain build and execution. */
			if (tier_c_rerun(scratch, rec, err))
				goto fail;
			nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
			continue;
		}
		const op_descriptor *op = op_descriptor_by_id(rec->op_id);
		gpointer user = op->deserialize(rec->params, rec->op_version);
		if (!user) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
			                       rec->record_id, rec->op_id);
			goto fail;
		}
		if (op->replay_pre) {
			GHashTable *kv = nde_kv_parse(rec->params);
			int rc = op->replay_pre(user, kv, scratch);
			g_hash_table_unref(kv);
			if (rc) {
				destroy_user(user);
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): replay preparation failed"),
				                       rec->record_id, rec->op_id);
				goto fail;
			}
		}
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_user(user);
			*err = g_strdup(_("out of memory"));
			goto fail;
		}
		args->fit = scratch;       /* PRIVATE fits — the whole point */
		args->op = op;
		args->user = user;
		args->nde_replay = TRUE;
		args->max_threads = com.max_thread;
		int rc = GPOINTER_TO_INT(generic_image_worker(args));
		free_generic_img_args(args);   /* frees user via its destructor too */
		if (rc) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) failed to apply"),
			                       rec->record_id, rec->op_id);
			goto fail;
		}
		/* Convergence C3: deposit the intermediate state so the NEXT edit
		 * restarts here instead of the baseline.  Pure cache — silent when
		 * the budget pref is 0 or the write fails. */
		nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
	}
	return scratch;

fail:
	clearfits(scratch);
	free(scratch);
	return NULL;
}

fits *nde_chain_replay(const nde_chain *chain, gchar **err) {
	g_return_val_if_fail(chain != NULL, NULL);
	g_return_val_if_fail(err != NULL, NULL);
	*err = NULL;
	if (!chain->replayable) {
		*err = g_strdup(_("chain is not replayable"));
		return NULL;
	}
	fits *scratch = nde_checkpoint_baseline_get(chain->item_id);
	if (!scratch) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		return NULL;
	}
	return replay_apply_records(scratch, chain, 0, chain->records->len, err);
}

fits *nde_chain_replay_tail(const nde_chain *chain, gchar **err) {
	g_return_val_if_fail(chain != NULL, NULL);
	g_return_val_if_fail(err != NULL, NULL);
	*err = NULL;
	if (!chain->tail_replayable) {
		*err = g_strdup(_("the editable tail is not replayable"));
		return NULL;
	}
	fits *start = chain->restart_ckpt_id > 0 ?
			nde_checkpoint_output_get(chain->restart_ckpt_id) :
			nde_checkpoint_baseline_get(chain->item_id);
	if (!start) {
		*err = g_strdup(chain->restart_ckpt_id > 0 ?
				_("failed to load the barrier checkpoint") :
				_("failed to load the baseline checkpoint"));
		return NULL;
	}
	return replay_apply_records(start, chain, chain->tail_start, chain->records->len, err);
}

/* ======================================================================= */
/* Amend-and-replay commit machinery (phase 3, P3.B)                       */
/* ======================================================================= */

/* Convergence C3: find the LATEST usable cached restart for an edit of
 * record @edit_id on the (trial) @chain, at-or-before the edit boundary.
 * Candidate states for restarting at member index j (state embodying
 * members [0..j-1]) are PRE(boundary id) and POST(members[j-1].id); for a
 * DELETE, PRE-states of members positioned after the deleted record embody
 * the deleted op and are invalid — bounding the j==e candidate by @edit_id
 * itself (its own PRE) covers both edit kinds.  Falls back to the chain's
 * phase-4 restart (barrier checkpoint or baseline) at tail_start.
 * Returns the restart state (caller owns) and *start_idx, or NULL+@err. */
static fits *resolve_edit_restart(const nde_chain *chain, guint e,
                                  gint64 boundary_pre_id,
                                  guint *start_idx, gchar **err) {
	/* @e is the POSITIONAL edit boundary in @chain (callers know it: the
	 * amended member's index, the deleted member's former index, or the
	 * earliest affected index of a reorder).  Position, not id — after a
	 * committed reorder chain member ids are no longer monotonic in array
	 * order.  @boundary_pre_id is the id whose PRE-state marks boundary e
	 * (the edit record itself for amend/delete; the first affected member
	 * of the ORIGINAL order for a reorder). */
	for (guint j = e; j > chain->tail_start; j--) {
		gint64 pre_id = (j == e) ? boundary_pre_id :
				((const nde_record *)g_ptr_array_index(chain->records, j))->record_id;
		fits *state = nde_snapstore_lookup(chain->item_id, pre_id, FALSE);
		if (!state) {
			const nde_record *prev = g_ptr_array_index(chain->records, j - 1);
			state = nde_snapstore_lookup(chain->item_id, prev->record_id, TRUE);
		}
		if (state) {
			*start_idx = j;
			return state;
		}
	}
	/* Fall back to the phase-4 restart point. */
	fits *start = chain->restart_ckpt_id > 0 ?
			nde_checkpoint_output_get(chain->restart_ckpt_id) :
			nde_checkpoint_baseline_get(chain->item_id);
	if (!start) {
		*err = g_strdup(chain->restart_ckpt_id > 0 ?
				_("failed to load the barrier checkpoint") :
				_("failed to load the baseline checkpoint"));
		return NULL;
	}
	*start_idx = chain->tail_start;
	return start;
}

/* Shared core of amend (new_params != NULL) and delete (new_params == NULL).
 * Conductor context: the caller holds SLOT_REPLAY, so capture, undo and
 * python cannot interleave with the replay or the commit. */
/* Checkpoint-built replay results carry pixels + WCS only (nde_snap stores
 * no header text, keywords or HISTORY).  A history edit rewrites pixel
 * parameters, not the document's identity: after a commit swap, move the
 * pre-edit metadata (keywords, header text, unknown keys, FITS HISTORY)
 * from the superseded fits back onto the target, keeping the REPLAYED WCS
 * and focal length — the chain may legitimately have changed those.
 * @old is the pre-edit fits the swap left holding the rich metadata; it is
 * about to be cleared, so ownership moves are plain pointer swaps. */
static void commit_restore_metadata(fits *target, fits *old) {
	g_rw_lock_writer_lock(&target->rwlock);
	/* Move the replayed WCS aside, swap the keyword structs wholesale, then
	 * put the replayed WCS back.  The pre-edit wcslib travels to @old inside
	 * its keywords struct, so clearfits(@old) frees it properly. */
	wcs_info replay_wcsdata      = target->keywords.wcsdata;
	struct wcsprm *replay_wcslib = target->keywords.wcslib;
	double replay_flen           = target->keywords.focal_length;
	fkeywords rich = old->keywords;
	old->keywords = target->keywords;         /* replay-minimal set... */
	old->keywords.wcslib = rich.wcslib;       /* ...plus the superseded WCS */
	target->keywords = rich;
	target->keywords.wcsdata      = replay_wcsdata;
	target->keywords.wcslib       = replay_wcslib;
	target->keywords.focal_length = replay_flen;

	char *h = target->header;
	target->header = old->header;
	old->header = h;
	gchar *u = target->unknown_keys;
	target->unknown_keys = old->unknown_keys;
	old->unknown_keys = u;
	GSList *hist = target->history;
	target->history = old->history;
	old->history = hist;
	g_rw_lock_writer_unlock(&target->rwlock);
}

static gboolean edit_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	if (amend_preview_installed()) {
		*err = g_strdup(_("another history step is being edited — close its dialog first"));
		return FALSE;
	}

	/* Locate the record among the LIVE records and pre-check it, for
	 * user-facing errors before any heavy work.  (nde_history_amend/delete
	 * re-validate under the mutex at commit time.) */
	gint item_id = 0;
	gboolean found = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->record_id == record_id) {
			found = TRUE;
			item_id = rec->target_item_id;
			if (new_params && !nde_record_amendable(rec)) {
				/* Amend needs editable params.  DELETE does not: removing
				 * an opaque record is well-defined — the trial chain never
				 * replays the deleted step, so deleting the one opaque
				 * record in a chain regains editability around it. */
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque and cannot be edited"),
				                       record_id, rec->op_id ? rec->op_id : "?");
			} else if (!new_params && !nde_record_deletable(rec)) {
				/* Structural (DOCUMENT scope) steps cannot resurrect what
				 * they destroyed or un-create a layer; compositing-state
				 * records (opacity/blend/visibility) describe live FLIS
				 * state, not pixels — deleting either would only make the
				 * log lie. */
				if (rec->scope == NDE_SCOPE_DOCUMENT)
					*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is a structural step and cannot be deleted"),
					                       record_id, rec->op_id ? rec->op_id : "?");
				else
					*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) records a layer-property change and cannot be deleted"),
					                       record_id, rec->op_id ? rec->op_id : "?");
			}
			break;
		}
	}
	if (live)
		g_ptr_array_unref(live);
	if (!found) {
		*err = g_strdup_printf(_("no live record with id %" G_GINT64_FORMAT), record_id);
		return FALSE;
	}
	if (*err)
		return FALSE;

	/* Position/freeze check on the CURRENT chain (phase-4 barrier model):
	 * edits are permitted only in the editable tail — records strictly
	 * after the last freeze cause — plus the delete of the last barrier
	 * itself.  Everything earlier is hard-frozen: recomputing it would
	 * require recomputing a later opaque step, which is impossible. */
	nde_chain *pos = nde_chain_build(item_id);
	gint pos_idx = -1;
	for (guint i = 0; i < pos->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(pos->records, i);
		if (rec->record_id == record_id) {
			pos_idx = (gint)i;
			break;
		}
	}
	if (pos_idx < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " does not affect this image's pixels"),
		                       record_id);
		nde_chain_free(pos);
		return FALSE;
	}
	gboolean in_tail = (guint)pos_idx >= pos->tail_start;
	gboolean is_last_barrier = !new_params && pos->tail_start > 0 &&
			(guint)pos_idx == pos->tail_start - 1 &&
			(g_array_index(pos->member_flags, guint8, pos_idx) & NDE_CHAIN_MEMBER_BARRIER);
	if (!in_tail && !is_last_barrier) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is locked by a later opaque step"),
		                       record_id);
		nde_chain_free(pos);
		return FALSE;
	}
	nde_chain_free(pos);

	/* For a delete, build the TRIAL chain: the deleted record is excluded
	 * from membership AND from the blockers, so its own opaqueness/mask
	 * state cannot veto its removal.  (Deleting the last barrier makes the
	 * trial's restart point fall back to the previous barrier/baseline.) */
	nde_chain *chain = chain_build_excluding(item_id, new_params ? 0 : record_id);
	if (!chain->tail_replayable) {
		GString *s = g_string_new(_("the history is not replayable: "));
		for (guint i = 0; i < chain->reasons->len; i++) {
			if (i)
				g_string_append(s, "; ");
			g_string_append(s, g_ptr_array_index(chain->reasons, i));
		}
		*err = g_string_free(s, FALSE);
		nde_chain_free(chain);
		return FALSE;
	}

	/* Amend: substitute the record's params in the trial chain.  (Delete
	 * needs no edit here — chain_build_excluding already omitted it.)
	 * Track the POSITIONAL boundary for the restart resolver: the amended
	 * member's index, or — for a delete — the index the removed member
	 * used to occupy, which equals its surviving-prefix length. */
	guint boundary = chain->records->len;
	if (new_params) {
		nde_record *target_rec = NULL;
		for (guint i = chain->tail_start; i < chain->records->len; i++) {
			nde_record *rec = g_ptr_array_index(chain->records, i);
			if (rec->record_id == record_id) {
				target_rec = rec;
				boundary = i;
				break;
			}
		}
		if (!target_rec) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " does not affect this image's pixels"),
			                       record_id);
			nde_chain_free(chain);
			return FALSE;
		}
		/* Validate the new params against the op before replaying. */
		const op_descriptor *op = op_descriptor_by_id(target_rec->op_id);
		gpointer trial = (op && op->deserialize) ?
				op->deserialize(new_params, target_rec->op_version) : NULL;
		if (!trial) {
			*err = g_strdup_printf(_("the new parameters for '%s' failed to parse"),
			                       target_rec->op_id ? target_rec->op_id : "?");
			nde_chain_free(chain);
			return FALSE;
		}
		destroy_user(trial);
		g_free(target_rec->params);
		target_rec->params = g_strdup(new_params);
	}

	gui_iface.set_progress(0.f, _("Recomputing edit history..."));
	/* Convergence C3 invalidation, BEFORE replay: pool states at-or-after
	 * the edit describe the OLD chain and must not survive (they are never
	 * consulted as restarts for THIS edit, but they would be stale for the
	 * next one).  The replay below re-deposits fresh states as it goes. */
	nde_snapstore_invalidate_from(item_id, record_id);
	/* Restart from the LATEST cached state at-or-before the edit (undo/redo
	 * entries, prior deposits, barrier checkpoints) rather than always the
	 * tail start — the frozen prefix's effect is embodied in the restart
	 * pixels either way. */
	if (!new_params)
		boundary = (guint)pos_idx;   /* the deleted member's former index */
	guint start_idx = 0;
	fits *start = resolve_edit_restart(chain, boundary, record_id, &start_idx, err);
	fits *result = start ? replay_apply_records(start, chain, start_idx, chain->records->len, err) : NULL;
	nde_chain_free(chain);
	if (!result) {
		/* Deposits made by a failed replay describe an uncommitted chain —
		 * drop them along with anything else at-or-after the edit. */
		nde_snapstore_invalidate_from(item_id, record_id);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}

	/* Resolve the target fits.  gfit for plain images; the layer's fit for
	 * FLIS (identical pointer when it is the active layer). */
	fits *target = gfit;
	if (item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(item_id);
		target = lay ? lay->fit : NULL;
	}
	if (!target) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		clearfits(result);
		free(result);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}

	/* Atomic commit: swap pixels, then the log.  `result` holds the OLD
	 * pixels after the swap, so a log-commit failure can restore them. */
	g_rw_lock_writer_lock(&target->rwlock);
	fits_swap_all_except_rwlock(target, result);
	g_rw_lock_writer_unlock(&target->rwlock);

	gboolean log_ok = new_params ? nde_history_amend(record_id, new_params, err)
	                             : nde_history_delete(record_id, err);
	if (!log_ok) {
		/* Should be unreachable (everything was validated, we own the
		 * slot); restore the old pixels so nothing is half-committed. */
		g_rw_lock_writer_lock(&target->rwlock);
		fits_swap_all_except_rwlock(target, result);
		g_rw_lock_writer_unlock(&target->rwlock);
		clearfits(result);
		free(result);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}
	commit_restore_metadata(target, result);
	clearfits(result);   /* the pre-edit pixels — superseded */
	free(result);

	/* No meta-undo (sketch §7): stale undo entries would restore pixels the
	 * log no longer describes. */
	undo_flush();

	invalidate_stats_from_fit(target);
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	if (target == gfit)
		notify_gfit_data_modified();
	gui_iface.set_progress(PROGRESS_DONE, _("Edit history updated"));
	return TRUE;
}

gboolean nde_amend_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(new_params != NULL, FALSE);
	return edit_execute(record_id, new_params, err);
}

gboolean nde_delete_execute(gint64 record_id, gchar **err) {
	return edit_execute(record_id, NULL, err);
}

gboolean nde_reorder_execute(gint64 record_id, gint64 anchor_id, gboolean after, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (amend_preview_installed()) {
		*err = g_strdup(_("another history step is being edited — close its dialog first"));
		return FALSE;
	}
	if (anchor_id <= 0 || record_id == anchor_id) {
		*err = g_strdup(_("invalid move target"));
		return FALSE;
	}

	/* Locate + policy-check the moved record among the LIVE records. */
	gint item_id = 0;
	gboolean found = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->record_id == record_id) {
			found = TRUE;
			item_id = rec->target_item_id;
			if (!nde_record_amendable(rec) || !nde_record_deletable(rec))
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) cannot be reordered"),
				                       record_id, rec->op_id ? rec->op_id : "?");
			break;
		}
	}
	if (live)
		g_ptr_array_unref(live);
	if (!found) {
		*err = g_strdup_printf(_("no live record with id %" G_GINT64_FORMAT), record_id);
		return FALSE;
	}
	if (*err)
		return FALSE;

	/* Member positions come from the chain — position, not id, is the
	 * order (ids stop being monotonic after the first reorder). */
	nde_chain *chain = nde_chain_build(item_id);
	gint i_old = -1, i_anchor_member = -1;
	for (guint i = 0; i < chain->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (rec->record_id == record_id)
			i_old = (gint)i;
		if (rec->record_id == anchor_id)
			i_anchor_member = (gint)i;
	}
	if (i_old < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " does not affect this image's pixels"),
		                       record_id);
		nde_chain_free(chain);
		return FALSE;
	}
	if (i_anchor_member < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is not part of the same editable history"),
		                       anchor_id);
		nde_chain_free(chain);
		return FALSE;
	}
	gint i_insert = i_anchor_member + (after ? 1 : 0);   /* pre-removal coords */
	gint dest = (i_insert > i_old) ? i_insert - 1 : i_insert;
	if (dest == i_old) {
		nde_chain_free(chain);
		return TRUE;   /* no-op move */
	}
	guint min_idx = (guint)MIN(i_old, dest);
	if (min_idx < chain->tail_start) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is locked by a later opaque step"),
		                       record_id);
		nde_chain_free(chain);
		return FALSE;
	}
	if (!chain->tail_replayable) {
		GString *rs = g_string_new(_("the history is not replayable: "));
		for (guint i = 0; i < chain->reasons->len; i++) {
			if (i)
				g_string_append(rs, "; ");
			g_string_append(rs, g_ptr_array_index(chain->reasons, i));
		}
		*err = g_string_free(rs, FALSE);
		nde_chain_free(chain);
		return FALSE;
	}

	/* Restart boundary: the first affected position of the ORIGINAL order;
	 * conservative invalidation floor: the smallest id among affected
	 * members (over-eviction is safe — the pool is cache). */
	gint64 boundary_pre_id =
			((const nde_record *)g_ptr_array_index(chain->records, min_idx))->record_id;
	gint64 inval_min = boundary_pre_id;
	for (guint i = min_idx; i < chain->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (rec->record_id < inval_min)
			inval_min = rec->record_id;
	}

	/* Log-side target: a GLOBAL before-id for nde_history_reorder — the
	 * live record right after the anchor for an 'after' move (0 = end). */
	gint64 log_before_id = anchor_id;
	if (after) {
		log_before_id = 0;
		GPtrArray *snap2 = nde_history_snapshot(NULL);
		for (guint i = 0; snap2 && i < snap2->len; i++) {
			const nde_record *rec = g_ptr_array_index(snap2, i);
			if (rec->record_id == anchor_id) {
				if (i + 1 < snap2->len)
					log_before_id = ((const nde_record *)g_ptr_array_index(snap2, i + 1))->record_id;
				break;
			}
		}
		if (snap2)
			g_ptr_array_unref(snap2);
	}

	/* Permute the trial chain in place (records + parallel flags). */
	nde_record *moved = g_ptr_array_steal_index(chain->records, (guint)i_old);
	g_ptr_array_insert(chain->records, dest, moved);
	guint8 fl = g_array_index(chain->member_flags, guint8, i_old);
	g_array_remove_index(chain->member_flags, (guint)i_old);
	g_array_insert_val(chain->member_flags, (guint)dest, fl);

	gui_iface.set_progress(0.f, _("Recomputing edit history..."));
	nde_snapstore_invalidate_from(item_id, inval_min);
	guint start_idx = 0;
	fits *start = resolve_edit_restart(chain, min_idx, boundary_pre_id, &start_idx, err);
	fits *result = start ? replay_apply_records(start, chain, start_idx, chain->records->len, err) : NULL;
	nde_chain_free(chain);
	if (!result) {
		nde_snapstore_invalidate_from(item_id, inval_min);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}

	/* Atomic commit — mirrors edit_execute's tail. */
	fits *target = gfit;
	if (item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(item_id);
		target = lay ? lay->fit : NULL;
	}
	if (!target) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		clearfits(result);
		free(result);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}
	g_rw_lock_writer_lock(&target->rwlock);
	fits_swap_all_except_rwlock(target, result);
	g_rw_lock_writer_unlock(&target->rwlock);

	gboolean log_ok = nde_history_reorder(record_id, log_before_id, err);
	if (!log_ok) {
		g_rw_lock_writer_lock(&target->rwlock);
		fits_swap_all_except_rwlock(target, result);
		g_rw_lock_writer_unlock(&target->rwlock);
		clearfits(result);
		free(result);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}
	commit_restore_metadata(target, result);
	clearfits(result);
	free(result);

	undo_flush();   /* no meta-undo (sketch §7) */
	invalidate_stats_from_fit(target);
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	if (target == gfit)
		notify_gfit_data_modified();
	gui_iface.set_progress(PROGRESS_DONE, _("Edit history updated"));
	return TRUE;
}

/* ---- off-worker conductor ---------------------------------------------- */
/*
 * Every history edit (amend / delete / reorder / preview) runs on a dedicated
 * conductor thread rather than the processing worker.  This is the single
 * replay path: a chain may contain a Tier-C step whose script must drive gfit
 * via cmd() on the (now free) worker while the conductor blocks on it, and a
 * thread cannot free the worker it is itself running on.  Pure Tier-A chains
 * take the same path — one code path, always exercised.
 *
 * The conductor holds SLOT_REPLAY for its whole run so GUI actions and other
 * scripts cannot steal the worker between the replayed script's commands; the
 * reservation gate admits the conductor's own submissions and those of the one
 * script it launches (see processing_thread.c).  Tier-A records still run
 * inline via generic_image_worker — its nde_replay path is thread-agnostic
 * (no idles, no worker-TLS use), so off-worker execution is equivalent.
 */
static gpointer conductor_trampoline(gpointer p) {
	GThreadFunc fn = ((gpointer *)p)[0];
	gpointer data = ((gpointer *)p)[1];
	g_free(p);
	replay_bind_conductor();
	fn(data);                 /* owns data; posts its own completion idle */
	replay_release_slot();
	return NULL;
}

/* Reserve the slot and run @fn on a fresh off-worker conductor thread.  Returns
 * FALSE (reserving nothing) if the slot is busy or reserved — the caller then
 * reports the failure and must free @data itself.
 *
 * In a command-script / CLI / python-command context the call BLOCKS until the
 * conductor finishes, preserving the "each command completes before the next
 * one runs" contract (a scripted flis_amend followed by save/flis_replay_check
 * must see the committed result, not race the replay).  A GUI panel action
 * (neither com.script nor com.python_command set) returns immediately and the
 * completion idle refreshes the display. */
static gboolean replay_conductor_start(GThreadFunc fn, gpointer data) {
	if (!replay_reserve_slot()) {
		siril_log_error(_("The processing thread is busy; try again when the "
		                  "current operation has finished\n"));
		return FALSE;
	}
	gpointer *ctx = g_new(gpointer, 2);
	ctx[0] = (gpointer)fn;
	ctx[1] = data;
	GThread *t = g_thread_new("nde-replay", conductor_trampoline, ctx);
	if (com.script || com.python_command)
		g_thread_join(t);     /* synchronous: caller waits for the commit */
	else
		g_thread_unref(t);    /* async: completion handled via idle */
	return TRUE;
}

/* ---- job wrappers ------------------------------------------------------ */

struct nde_edit_job {
	gint64 record_id;
	gchar *new_params;   /* NULL = delete or reorder */
	gint64 anchor_id;    /* != 0 = reorder */
	gboolean after;
};

static gboolean nde_edit_done_idle(gpointer p) {
	(void)p;
	gui_iface.redraw_image(REDRAW_ALL);
	gui_iface.flis_gui_update();
	/* The edit flushed the undo stacks (no meta-undo) — refresh the
	 * undo/redo buttons or they stay sensitive over an empty stack. */
	gui_iface.update_menu_item();
	/* No end_generic: the conductor is not a worker job and frees the slot
	 * itself via replay_release_slot (calling stop_processing_thread here
	 * would set cancel_flag and poison the next operation). */
	return G_SOURCE_REMOVE;
}

static gpointer nde_edit_worker(gpointer p) {
	struct nde_edit_job *job = p;
	gchar *errmsg = NULL;
	gboolean ok = job->anchor_id ?
			nde_reorder_execute(job->record_id, job->anchor_id, job->after, &errmsg) :
			edit_execute(job->record_id, job->new_params, &errmsg);
	if (!ok)
		siril_log_error(_("Edit history change failed: %s\n"), errmsg ? errmsg : "?");
	g_free(errmsg);
	g_free(job->new_params);
	g_free(job);
	siril_add_idle(nde_edit_done_idle, NULL);
	return GINT_TO_POINTER(ok ? 0 : 1);
}

/* flis_replay_check worker (P2.D, conductor-hosted since nde-phase5): builds
 * and replays the chain, logging validation reasons or the deviation report.
 * Runs under SLOT_REPLAY so Tier-C script re-runs behave exactly as they do
 * in the edit paths (capture/undo suppressed, script commands admitted). */
static gpointer replay_check_worker(gpointer p) {
	gint item_id = GPOINTER_TO_INT(p);
	nde_chain *chain = nde_chain_build(item_id);
	if (chain->records->len == 0) {
		siril_log_info(_("No replayable records — nothing to check\n"));
	} else if (!chain->replayable && !chain->tail_replayable) {
		siril_log_warning(_("History is not replayable:\n"));
		for (guint i = 0; i < chain->reasons->len; i++)
			siril_log_message("  - %s\n", (char *)g_ptr_array_index(chain->reasons, i));
	} else if (!chain->replayable && chain->records->len == chain->tail_start) {
		/* barrier-last: nothing beyond its checkpoint to verify */
		siril_log_info(_("%u step(s) are frozen behind an opaque barrier and the last step is the barrier itself — nothing to verify\n"),
		               chain->tail_start);
	} else {
		guint frozen = chain->replayable ? 0 : chain->tail_start;
		if (frozen)
			siril_log_info(_("%u step(s) are frozen behind an opaque barrier; verifying the last %u step(s) from its checkpoint\n"),
			               frozen, chain->records->len - frozen);
		gchar *errmsg = NULL;
		fits *result = chain->replayable ? nde_chain_replay(chain, &errmsg)
		                                 : nde_chain_replay_tail(chain, &errmsg);
		if (!result) {
			siril_log_error(_("Replay failed: %s\n"), errmsg ? errmsg : "?");
		} else {
			fits *current = gfit;
			if (item_id >= 0) {
				flis_layer_t *lay = flis_layer_get_by_id(item_id);
				current = lay ? lay->fit : NULL;
			}
			if (!current) {
				siril_log_error(_("Cannot locate the current pixels to compare against\n"));
			} else {
				g_rw_lock_reader_lock(&current->rwlock);
				if (result->rx != current->rx || result->ry != current->ry
				    || result->naxes[2] != current->naxes[2]
				    || result->type != current->type) {
					siril_log_warning(_("Replayed %u record(s), but the result geometry differs from the current image (%ux%ux%ld vs %ux%ux%ld)\n"),
					                  chain->records->len - (chain->replayable ? 0 : chain->tail_start),
					                  result->rx, result->ry, result->naxes[2],
					                  current->rx, current->ry, current->naxes[2]);
				} else {
					size_t n = (size_t)current->rx * current->ry
					           * (current->naxes[2] ? current->naxes[2] : 1);
					double max_dev = 0.0, sum_dev = 0.0;
					for (size_t i = 0; i < n; i++) {
						double a, b;
						if (current->type == DATA_FLOAT) {
							a = result->fdata[i];
							b = current->fdata[i];
						} else {
							a = result->data[i];
							b = current->data[i];
						}
						double d = fabs(a - b);
						if (d > max_dev) max_dev = d;
						sum_dev += d;
					}
					siril_log_info(_("Replayed %u record(s) successfully: max deviation %.3g, mean %.3g (small numerical drift is expected)\n"),
					               chain->records->len - (chain->replayable ? 0 : chain->tail_start),
					               max_dev, sum_dev / n);
				}
				g_rw_lock_reader_unlock(&current->rwlock);
			}
			clearfits(result);
			free(result);
		}
		g_free(errmsg);
	}
	nde_chain_free(chain);
	/* No end_generic / completion idle: nothing was modified, and the
	 * conductor frees the slot itself via replay_release_slot. */
	return GINT_TO_POINTER(0);
}

gboolean nde_replay_check_start(gint item_id) {
	return replay_conductor_start(replay_check_worker, GINT_TO_POINTER(item_id));
}

static gboolean nde_edit_start(gint64 record_id, const gchar *new_params) {
	struct nde_edit_job *job = g_new0(struct nde_edit_job, 1);
	job->record_id = record_id;
	job->new_params = g_strdup(new_params);
	if (!replay_conductor_start(nde_edit_worker, job)) {
		g_free(job->new_params);
		g_free(job);
		return FALSE;
	}
	return TRUE;
}

gboolean nde_amend_start(gint64 record_id, const gchar *new_params) {
	g_return_val_if_fail(new_params != NULL, FALSE);
	return nde_edit_start(record_id, new_params);
}

gboolean nde_delete_start(gint64 record_id) {
	return nde_edit_start(record_id, NULL);
}

gboolean nde_reorder_start(gint64 record_id, gint64 anchor_id, gboolean after) {
	struct nde_edit_job *job = g_new0(struct nde_edit_job, 1);
	job->record_id = record_id;
	job->anchor_id = anchor_id;
	job->after = after;
	if (!replay_conductor_start(nde_edit_worker, job)) {
		g_free(job);
		return FALSE;
	}
	return TRUE;
}

/* ======================================================================= */
/* Amend preview (convergence C4)                                          */
/* ======================================================================= */

/* Single instance.  While `installed`, the record's target fits holds the
 * synthesized pre-K state and `saved` holds the true pixels, swapped out
 * wholesale with fits_swap_all_except_rwlock — the restore is the reverse
 * swap, bit-exact including metadata.  The heavy transitions (begin/end
 * _execute) run on the replay conductor holding SLOT_REPLAY, so the
 * reservation serializes them; apv_mutex is a leaf guard for the flag reads
 * that happen on other threads (GUI enablement, the edit_execute guard). */
static GMutex apv_mutex;
static struct {
	gboolean active;     /* start accepted, until end / failed begin */
	gboolean installed;  /* pre-K currently swapped into the target */
	gint64   record_id;
	gint     item_id;
	gchar   *op_id;
	gint     op_version;
	gchar   *params;     /* the record's current kv blob, for pre-fill */
	fits    *saved;      /* true pixels while installed (owned) */
	nde_amend_preview_ready_fn on_ready;
	gpointer on_ready_user;
} apv;

static gboolean amend_preview_installed(void) {
	g_mutex_lock(&apv_mutex);
	gboolean r = apv.installed;
	g_mutex_unlock(&apv_mutex);
	return r;
}

gboolean nde_amend_preview_active(void) {
	g_mutex_lock(&apv_mutex);
	gboolean r = apv.active;
	g_mutex_unlock(&apv_mutex);
	return r;
}

gint64       nde_amend_preview_record_id(void)  { return apv.record_id; }
const gchar *nde_amend_preview_op_id(void)      { return apv.op_id; }
gint         nde_amend_preview_op_version(void) { return apv.op_version; }
const gchar *nde_amend_preview_params(void)     { return apv.params; }

/* The fits a chain edit commits into: gfit for plain images, the layer's
 * fit for FLIS items (same pointer when it is the active layer). */
static fits *edit_target_fits(gint item_id) {
	if (item_id < 0)
		return gfit;
	flis_layer_t *lay = flis_layer_get_by_id(item_id);
	return lay ? lay->fit : NULL;
}

/* Swap @incoming into @target under the display-quiesce + writer-lock
 * discipline of the commit path, then refresh the derived state.  After
 * the call @incoming holds what @target held before. */
static void apv_swap_into_target(fits *target, fits *incoming) {
	gboolean is_display = (target == gfit);
	if (is_display)
		gui_iface.set_suppress_redraws(TRUE);
	g_rw_lock_writer_lock(&target->rwlock);
	fits_swap_all_except_rwlock(target, incoming);
	g_rw_lock_writer_unlock(&target->rwlock);
	if (is_display)
		gui_iface.set_suppress_redraws(FALSE);
	invalidate_stats_from_fit(target);
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	if (is_display)
		notify_gfit_data_modified();
}

static void apv_clear_state_locked(void) {
	apv.active = FALSE;
	apv.installed = FALSE;
	apv.record_id = 0;
	apv.item_id = -1;
	g_free(apv.op_id);    apv.op_id = NULL;
	g_free(apv.params);   apv.params = NULL;
	apv.op_version = 0;
	apv.saved = NULL;
}

gboolean nde_amend_preview_begin_execute(gint64 record_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	if (amend_preview_installed()) {
		*err = g_strdup(_("another history step is already being edited"));
		return FALSE;
	}
	if (gui_iface.is_preview_active()) {
		/* Some dialog's preview backup is armed — its backup/restore
		 * cycle and this install would fight over the same pixels. */
		*err = g_strdup(_("close the open preview dialog first"));
		return FALSE;
	}

	/* Locate + policy-check the record among the LIVE records. */
	gint item_id = 0;
	const gchar *op_id = NULL, *params = NULL;
	gint op_version = 0;
	gboolean found = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->record_id == record_id) {
			found = TRUE;
			item_id = rec->target_item_id;
			op_id = rec->op_id;
			op_version = rec->op_version;
			params = rec->params;
			if (!nde_record_amendable(rec) || !rec->params)
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque and cannot be edited"),
				                       record_id, rec->op_id ? rec->op_id : "?");
			break;
		}
	}
	gchar *op_id_copy = op_id ? g_strdup(op_id) : NULL;
	gchar *params_copy = params ? g_strdup(params) : NULL;
	if (live)
		g_ptr_array_unref(live);
	if (!found) {
		*err = g_strdup_printf(_("no live record with id %" G_GINT64_FORMAT), record_id);
		goto fail_free;
	}
	if (*err)
		goto fail_free;

	/* The dialog's preview pipeline works on the DISPLAYED image; refuse a
	 * record targeting a non-active FLIS layer rather than previewing one
	 * image while showing another. */
	if (item_id != nde_checkpoint_active_item_id()) {
		*err = g_strdup(_("this step targets another layer — make that layer active first"));
		goto fail_free;
	}

	/* Position/freeze check + synthesis on the CURRENT chain. */
	nde_chain *chain = nde_chain_build(item_id);
	gint e = -1;
	for (guint i = 0; i < chain->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (rec->record_id == record_id) {
			e = (gint)i;
			break;
		}
	}
	if (e < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " does not affect this image's pixels"),
		                       record_id);
		nde_chain_free(chain);
		goto fail_free;
	}
	if ((guint)e < chain->tail_start) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is locked by a later opaque step"),
		                       record_id);
		nde_chain_free(chain);
		goto fail_free;
	}

	/* Synthesize pre-K: latest cached state at-or-before K, then apply the
	 * members between it and K.  The deposits made along the way are valid
	 * for the CURRENT chain (nothing has been edited yet) and make the
	 * eventual amend's tail replay restart adjacent to K. */
	guint start_idx = 0;
	fits *start = resolve_edit_restart(chain, (guint)e, record_id, &start_idx, err);
	fits *pre_k = start ? replay_apply_records(start, chain, start_idx, (guint)e, err) : NULL;
	nde_chain_free(chain);
	if (!pre_k)
		goto fail_free;

	fits *target = edit_target_fits(item_id);
	if (!target) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		clearfits(pre_k);
		free(pre_k);
		goto fail_free;
	}

	/* Install: after the swap pre_k holds the TRUE pixels — that is the
	 * stash the end path restores from. */
	apv_swap_into_target(target, pre_k);

	g_mutex_lock(&apv_mutex);
	apv.active = TRUE;
	apv.installed = TRUE;
	apv.record_id = record_id;
	apv.item_id = item_id;
	apv.op_id = op_id_copy;
	apv.op_version = op_version;
	apv.params = params_copy;
	apv.saved = pre_k;
	g_mutex_unlock(&apv_mutex);
	return TRUE;

fail_free:
	g_free(op_id_copy);
	g_free(params_copy);
	return FALSE;
}

gboolean nde_amend_preview_end_execute(gboolean apply, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	g_mutex_lock(&apv_mutex);
	if (!apv.active) {
		g_mutex_unlock(&apv_mutex);
		if (apply) {
			*err = g_strdup(_("no amend preview is active"));
			return FALSE;
		}
		return TRUE;   /* tolerated: defensive cancels from destroy handlers */
	}
	gint64 record_id = apv.record_id;
	gint item_id = apv.item_id;
	fits *saved = apv.saved;
	apv_clear_state_locked();
	g_mutex_unlock(&apv_mutex);

	/* Restore the true pixels FIRST (the plan's ordering contract): gfit
	 * must hold the real image again before any amend runs, so a failed
	 * amend changes nothing.  `saved` receives the pre-K/preview pixels
	 * from the swap — superseded either way. */
	if (saved) {
		fits *target = edit_target_fits(item_id);
		if (target) {
			apv_swap_into_target(target, saved);
		} else {
			/* The layer vanished under the preview (should be unreachable:
			 * edits are blocked while installed).  Nothing to restore into. */
			siril_log_error(_("Amend preview: the target layer no longer exists\n"));
		}
		clearfits(saved);
		free(saved);
	}

	if (apply)
		return edit_execute(record_id, new_params, err);
	return TRUE;
}

/* ---- GUI wrappers (job spawning + main-thread ready callback) ---------- */

static gboolean apv_ready_idle(gpointer p) {
	gboolean ok = (GPOINTER_TO_INT(p) == 0);
	nde_amend_preview_ready_fn cb = apv.on_ready;
	gpointer user = apv.on_ready_user;
	apv.on_ready = NULL;
	apv.on_ready_user = NULL;
	gui_iface.redraw_image(REDRAW_ALL);
	/* The conductor already released SLOT_REPLAY before this idle ran, so the
	 * slot is free: the dialog the callback opens can run its preview jobs
	 * immediately. */
	if (cb)
		cb(ok, user);
	return G_SOURCE_REMOVE;
}

static gpointer apv_begin_worker(gpointer p) {
	gint64 record_id = *(gint64 *)p;
	g_free(p);
	gchar *errmsg = NULL;
	gboolean ok = nde_amend_preview_begin_execute(record_id, &errmsg);
	if (!ok)
		siril_log_error(_("Cannot edit this history step: %s\n"), errmsg ? errmsg : "?");
	g_free(errmsg);
	siril_add_idle(apv_ready_idle, GINT_TO_POINTER(ok ? 0 : 1));
	return GINT_TO_POINTER(ok ? 0 : 1);
}

gboolean nde_amend_preview_start(gint64 record_id,
                                 nde_amend_preview_ready_fn on_ready,
                                 gpointer user) {
	if (nde_amend_preview_active()) {
		siril_log_error(_("Another history step is already being edited\n"));
		return FALSE;
	}
	apv.on_ready = on_ready;
	apv.on_ready_user = user;
	gint64 *id = g_new(gint64, 1);
	*id = record_id;
	if (!replay_conductor_start(apv_begin_worker, id)) {
		apv.on_ready = NULL;
		apv.on_ready_user = NULL;
		g_free(id);
		return FALSE;
	}
	return TRUE;
}

struct apv_end_job {
	gboolean apply;
	gchar *new_params;
};

static gpointer apv_end_worker(gpointer p) {
	struct apv_end_job *job = p;
	gchar *errmsg = NULL;
	gboolean ok = nde_amend_preview_end_execute(job->apply, job->new_params, &errmsg);
	if (!ok)
		siril_log_error(_("Edit history change failed: %s\n"), errmsg ? errmsg : "?");
	g_free(errmsg);
	g_free(job->new_params);
	g_free(job);
	siril_add_idle(nde_edit_done_idle, NULL);
	return GINT_TO_POINTER(ok ? 0 : 1);
}

gboolean nde_amend_preview_end(gboolean apply, const gchar *new_params) {
	if (!nde_amend_preview_active())
		return !apply;
	struct apv_end_job *job = g_new0(struct apv_end_job, 1);
	job->apply = apply;
	job->new_params = g_strdup(new_params);
	if (!replay_conductor_start(apv_end_worker, job)) {
		g_free(job->new_params);
		g_free(job);
		return FALSE;
	}
	return TRUE;
}
