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

/* Validate one chain member; every problem becomes a reason. */
static void validate_member(nde_chain *chain, const nde_record *rec) {
	if (rec->tier != NDE_TIER_A) {
		add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) is opaque — not replayable"),
		           rec->record_id, rec->op_id ? rec->op_id : "?");
		return;
	}
	if (rec->mask_active)
		add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) ran with an active mask — mask replay is not supported yet"),
		           rec->record_id, rec->op_id);
	const op_descriptor *op = op_descriptor_by_id(rec->op_id);
	if (!op) {
		add_reason(chain, _("record %" G_GINT64_FORMAT ": unknown operation '%s'"),
		           rec->record_id, rec->op_id ? rec->op_id : "?");
		return;
	}
	if (!op->deserialize) {
		add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) has no parameter deserializer"),
		           rec->record_id, rec->op_id);
		return;
	}
	if (rec->op_version > op->version) {
		add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) was written by a newer version (v%d > v%d)"),
		           rec->record_id, rec->op_id, rec->op_version, op->version);
		return;
	}
	gpointer user = op->deserialize(rec->params, rec->op_version);
	if (!user) {
		add_reason(chain, _("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		           rec->record_id, rec->op_id);
		return;
	}
	destroy_user(user);
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
	gboolean is_flis = is_current_image_flis();

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
		case NDE_SCOPE_DOCUMENT:
			if (rec->target_item_id == item_id &&
			    (!g_strcmp0(rec->op_id, "layer.merge_down") ||
			     !g_strcmp0(rec->op_id, "document.flatten"))) {
				/* The item's pixels were destructively replaced by a
				 * composite of other layers — nothing to replay from. */
				add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) replaced this layer's pixels with a composite — not replayable"),
				           rec->record_id, rec->op_id);
			}
			/* layer.add is embodied in the baseline; property/structure
			 * records don't touch this item's pixels: ignore. */
			break;
		default:
			break;
		}
		if (member) {
			validate_member(chain, rec);
			g_ptr_array_add(chain->records, rec);
		}
	}

	if (chain->records->len && !nde_checkpoint_baseline_exists(item_id))
		add_reason(chain, _("no baseline checkpoint — the file predates baselines, or the history began before this build"));

	chain->replayable = (chain->reasons->len == 0);
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
	g_free(chain);
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

	for (guint i = 0; i < chain->records->len; i++) {
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
	}
	return scratch;

fail:
	clearfits(scratch);
	free(scratch);
	return NULL;
}

/* ======================================================================= */
/* Amend-and-replay commit machinery (phase 3, P3.B)                       */
/* ======================================================================= */

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

	/* For a delete, build the TRIAL chain directly: the deleted record is
	 * excluded from membership AND from the blockers, so its own
	 * opaqueness/mask state cannot veto its removal. */
	nde_chain *chain = chain_build_excluding(item_id, new_params ? 0 : record_id);
	if (!chain->replayable) {
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
	 * needs no edit here — chain_build_excluding already omitted it.) */
	if (new_params) {
		nde_record *target_rec = NULL;
		for (guint i = 0; i < chain->records->len; i++) {
			nde_record *rec = g_ptr_array_index(chain->records, i);
			if (rec->record_id == record_id) {
				target_rec = rec;
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
	fits *result = nde_chain_replay(chain, err);
	nde_chain_free(chain);
	if (!result) {
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

/* ---- job wrappers ------------------------------------------------------ */

struct nde_edit_job {
	gint64 record_id;
	gchar *new_params;   /* NULL = delete */
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
	gboolean ok = edit_execute(job->record_id, job->new_params, &errmsg);
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
