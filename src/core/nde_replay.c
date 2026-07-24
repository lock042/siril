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

static void add_reason(nde_chain *chain, const char *fmt, ...) G_GNUC_PRINTF(2, 3);
static void add_reason(nde_chain *chain, const char *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	g_ptr_array_add(chain->reasons, g_strdup_vprintf(fmt, ap));
	va_end(ap);
}

/* First reason @rec cannot be replayed (heap string, caller frees), or NULL
 * when it replays.  The chain build decides whether an invalid member is a
 * hard blocker (no output checkpoint) or a barrier restart point. */
static gchar *member_invalid_reason(const nde_record *rec) {
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

/* Apply members [from..end) to @scratch (consumed on failure).  Returns
 * @scratch on success, NULL + @err on failure. */
static fits *replay_apply_records(fits *scratch, const nde_chain *chain,
                                  guint from, gchar **err) {
	for (guint i = from; i < chain->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (!processing_should_continue()) {
			*err = g_strdup(_("cancelled"));
			goto fail;
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
	return replay_apply_records(scratch, chain, 0, err);
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
	return replay_apply_records(start, chain, chain->tail_start, err);
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
 * Job context: the caller owns the processing slot, so capture, undo and
 * python cannot interleave with the replay or the commit. */
static gboolean edit_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

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
	fits *result = start ? replay_apply_records(start, chain, start_idx, err) : NULL;
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
	fits *result = start ? replay_apply_records(start, chain, start_idx, err) : NULL;
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
	return end_generic(NULL);
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

static gboolean nde_edit_start(gint64 record_id, const gchar *new_params) {
	if (processing_is_reserved_for_python()) {
		siril_log_error(_("The processing thread is reserved by a Python script; try again later\n"));
		return FALSE;
	}
	struct nde_edit_job *job = g_new0(struct nde_edit_job, 1);
	job->record_id = record_id;
	job->new_params = g_strdup(new_params);
	if (!start_in_new_thread(nde_edit_worker, job)) {
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
	if (processing_is_reserved_for_python()) {
		siril_log_error(_("The processing thread is reserved by a Python script; try again later\n"));
		return FALSE;
	}
	struct nde_edit_job *job = g_new0(struct nde_edit_job, 1);
	job->record_id = record_id;
	job->anchor_id = anchor_id;
	job->after = after;
	if (!start_in_new_thread(nde_edit_worker, job)) {
		g_free(job);
		return FALSE;
	}
	return TRUE;
}
