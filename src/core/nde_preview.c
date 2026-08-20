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
 * The asynchronous, GUI-facing half of NDE editing: the off-worker conductor
 * every history edit runs on, the amend preview and edit-at state machine it
 * drives, and the region-scoped tail replay that keeps a preview cheap.
 *
 * These share one piece of state — the `apv` block — and one rule: a preview
 * has the target's real pixels swapped out, so nothing else may commit while
 * one is installed.  Keeping them together is what makes that rule checkable.
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
#include "core/fits_region.h"
#include "core/nde_history.h"
#include "core/nde_op_class.h"
#include "core/nde_checkpoint.h"
#include "core/masks.h"
#include "core/nde_compositing.h"
#include "core/nde_composite.h"
#include "core/nde_joint.h"
#include "core/nde_snapstore.h"
#include "core/nde_replay.h"
#include "core/nde_replay_internal.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/siril_pythonmodule.h"

#include "core/nde_replay_internal.h"

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
	/* After the release, so the slot is free the moment the pointer says it
	 * is — and after fn queued its completion idle, so the two land in that
	 * order and the cursor clears onto a redrawn image rather than a stale
	 * one (both go through g_idle_add at the same priority). */
	gui_iface.set_busy(FALSE);
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
	/* Replaying a chain re-runs every step of an item's history, which on a
	 * long one is the slowest thing the panel can ask for and shows nothing
	 * while it happens: the conductor runs off the worker, so the window stays
	 * responsive and looks idle.  Set once here rather than at each of the
	 * eight dispatch sites, and cleared by the trampoline, so the two cannot
	 * drift apart — and after the reservation, so a refused start leaves the
	 * pointer alone. */
	gui_iface.set_busy(TRUE);
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
	/* The notify calls made DURING the edit ran on the conductor with the
	 * job slot held, and notify_gfit_data_modified skips its remap in
	 * that state (single_image.c's mid-job guard) — so the view buffers
	 * still hold the pre-edit pixels.  Remap here, on the main thread
	 * with the pixel work finished, or the recomputed result stays
	 * invisible until the next incidental remap (the "recompute runs but
	 * nothing changes on screen" report). */
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	gui_iface.remap_all_vports();
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
			nde_edit_execute(job->record_id, job->new_params, &errmsg);
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

/* Undoing a composite runs on the conductor like every other history edit:
 * it replays each consumed layer's chain, which needs the job slot. */
static gpointer nde_composite_undo_worker(gpointer p) {
	gint64 record_id = *(gint64 *)p;
	g_free(p);
	gchar *errmsg = NULL;
	if (!nde_composite_undo_execute(record_id, &errmsg))
		siril_log_error(_("The merge could not be undone: %s\n"),
		                errmsg ? errmsg : "?");
	g_free(errmsg);
	siril_add_idle(nde_edit_done_idle, NULL);
	return GINT_TO_POINTER(0);
}

gboolean nde_composite_undo_start(gint64 record_id) {
	gint64 *id = g_new(gint64, 1);
	*id = record_id;
	if (!replay_conductor_start(nde_composite_undo_worker, id)) {
		g_free(id);
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
 * wholesale with nde_commit_pixels — the restore is the reverse swap, bit-exact
 * including metadata.  The mask slot stays put through both, as it does for
 * every other commit: the mask is a different item, and a preview of these
 * pixels says nothing about it.  The heavy transitions (begin/end
 * _execute) run on the replay conductor holding SLOT_REPLAY, so the
 * reservation serializes them; apv_mutex is a leaf guard for the flag reads
 * that happen on other threads (GUI enablement, the nde_edit_execute guard). */
static GMutex apv_mutex;
static struct {
	gboolean active;     /* start accepted, until end / failed begin */
	gboolean installed;  /* pre-K currently swapped into the target */
	gboolean insert;     /* edit-at mode: new steps are inserted before K */
	/* The item has no layer of its own — a merge or flatten consumed it — so
	 * its pre-K state is shown on the DISPLAY instead, and captures made
	 * while the mode is open are stamped with the borrowed item rather than
	 * with whatever layer happens to be active.  Nothing is committed to a
	 * layer on the way out; the composites that consume the item are
	 * recomputed instead. */
	gboolean borrowed;
	/* A region tail replay failed once; stop attempting it for the rest of
	 * this preview (nde_region_tail_apply). */
	gboolean tail_failed;
	gint64   record_id;
	gint     item_id;
	gchar   *op_id;
	gint     op_version;
	gchar   *params;     /* the record's current kv blob, for pre-fill */
	fits    *saved;      /* true pixels while installed (owned) */
	nde_amend_preview_ready_fn on_ready;
	gpointer on_ready_user;
} apv;

gboolean nde_amend_preview_installed(void) {
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

/* The item a capture must be stamped with while an insertion is open on an
 * item that has no layer: without this the record would carry whatever layer
 * is active, which is how an insertion aimed at a merge input used to end up
 * appended to the merged result instead.  0 when not borrowing. */
gint nde_edit_at_borrowed_item(void) {
	g_mutex_lock(&apv_mutex);
	gint id = (apv.active && apv.insert && apv.borrowed) ? apv.item_id : 0;
	g_mutex_unlock(&apv_mutex);
	return id;
}

gint64       nde_amend_preview_record_id(void)  { return apv.record_id; }
const gchar *nde_amend_preview_op_id(void)      { return apv.op_id; }
gint         nde_amend_preview_op_version(void) { return apv.op_version; }
const gchar *nde_amend_preview_params(void)     { return apv.params; }

/* Swap @incoming into @target under the display-quiesce + writer-lock
 * discipline of the commit path, then refresh the derived state.  After
 * the call @incoming holds what @target held before. */
static void apv_swap_into_target(fits *target, fits *incoming) {
	gboolean is_display = (target == gfit);
	if (is_display)
		gui_iface.set_suppress_redraws(TRUE);
	g_rw_lock_writer_lock(&target->rwlock);
	nde_commit_pixels(target, incoming);
	g_rw_lock_writer_unlock(&target->rwlock);
	if (is_display)
		gui_iface.set_suppress_redraws(FALSE);
	invalidate_stats_from_fit(target);
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	if (is_display)
		notify_gfit_data_modified();
	gui_iface.redraw_image(REDRAW_ALL);
}

static void apv_clear_state_locked(void) {
	apv.active = FALSE;
	apv.installed = FALSE;
	apv.insert = FALSE;
	apv.borrowed = FALSE;
	apv.tail_failed = FALSE;
	apv.record_id = 0;
	apv.item_id = -1;
	g_free(apv.op_id);    apv.op_id = NULL;
	g_free(apv.params);   apv.params = NULL;
	apv.op_version = 0;
	apv.saved = NULL;
}

/* Install the state just BEFORE record @record_id into its target image.
 * Shared by the two modes that need it: amend preview (@insert FALSE — the
 * record's own dialog reopens against the pre-K state) and edit-at (@insert
 * TRUE — ordinary operations run there and are inserted before the record). */
static gboolean apv_begin_execute(gint64 record_id, gboolean insert, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	if (nde_amend_preview_installed()) {
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
			/* Edit-at only needs the record as a POSITION: nothing about it
			 * is re-run, so an opaque anchor is fine. */
			if (!insert && (!nde_record_amendable(rec) || !rec->params))
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
	 * image while showing another.  An item a merge or flatten consumed is the
	 * exception that proves it: there is no layer to make active, so the
	 * display is LENT to it for the duration — for an insertion and for an
	 * amend preview alike.  The amend used to be refused here, on the grounds
	 * that the History could edit the step without a preview; what that meant
	 * in practice was that flattening a document demoted every one of its
	 * stretches to the raw key/value grid. */
	gboolean borrow = FALSE;
	if (item_id != nde_checkpoint_active_item_id()) {
		if (nde_item_is_retained_input(item_id)) {
			borrow = TRUE;
		} else {
			*err = g_strdup(_("this step targets another layer — make that layer active first"));
			goto fail_free;
		}
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
		/* Two ways to be outside the tail, and for an insertion the second
		 * one is worth naming: the anchor is ITSELF the opaque step, so the
		 * refusal is not "something later locked you" but "you could not be
		 * re-run over what you inserted". */
		gboolean self_barrier = insert &&
				(g_array_index(chain->member_flags, guint8, e) & NDE_CHAIN_MEMBER_BARRIER);
		*err = self_barrier ?
			g_strdup_printf(_("record %" G_GINT64_FORMAT " is an opaque step: it cannot be recomputed, so nothing can be inserted before it"),
			                record_id) :
			g_strdup_printf(_("record %" G_GINT64_FORMAT " is locked by a later opaque step"),
			                record_id);
		nde_chain_free(chain);
		goto fail_free;
	}
	/* The forward replay would have to restart from the layer's position
	 * just before the anchor, and the live position already embodies every
	 * record after it.  Nothing records the intermediate positions, so refuse
	 * rather than replay from a position that is not the right one. */
	if (insert && chain->has_geometry) {
		*err = g_strdup(_("this image's history moves the layer on the canvas; "
		                  "steps cannot be inserted into it yet"));
		nde_chain_free(chain);
		goto fail_free;
	}
	/* Edit-at ends by replaying members [e..end] forward over the inserted
	 * work, so refuse up front if that replay could not run — the check is
	 * cheap here and there is no way back once pixels have changed. */
	if (insert && !chain->tail_replayable) {
		GString *s = g_string_new(_("the steps after this one cannot be recomputed: "));
		for (guint i = 0; i < chain->reasons->len; i++) {
			if (i)
				g_string_append(s, "; ");
			g_string_append(s, g_ptr_array_index(chain->reasons, i));
		}
		*err = g_string_free(s, FALSE);
		nde_chain_free(chain);
		goto fail_free;
	}

	/* Synthesize pre-K: latest cached state at-or-before K, then apply the
	 * members between it and K.  The deposits made along the way are valid
	 * for the CURRENT chain (nothing has been edited yet) and make the
	 * eventual amend's tail replay restart adjacent to K. */
	guint start_idx = 0;
	nde_state *start = nde_edit_restart_state(chain, (guint)e, record_id,
	                                          &start_idx, err);
	/* Preview only: the true pixels — and the layer's real position — come
	 * back untouched on exit, so nothing is committed and nothing is carried. */
	if (start)
		start->has_pos = FALSE;
	/* A NULL restart with no error means the item was born of a composite and
	 * starts from no state at all — the same convention nde_edit_execute and
	 * nde_reorder_execute follow.  At K == 0 that leaves nothing to synthesize:
	 * there IS no state before an item's own origin, and saying so beats
	 * failing with an empty reason, which is what this did. */
	if (!start && !*err && e == 0)
		*err = g_strdup(_("this step is what produced this image — there is no "
		                  "earlier state of it to edit against"));
	fits *pre_k = (!start && *err) ? NULL :
			nde_state_release(nde_replay_apply_records(start, chain, start_idx,
			                                           (guint)e, err));
	nde_chain_free(chain);
	if (!pre_k) {
		if (!*err)
			*err = g_strdup(_("the state before this step could not be rebuilt"));
		goto fail_free;
	}

	/* A borrowed item has no fits of its own; the display holds its state for
	 * the duration and gives it back on the way out. */
	fits *target = borrow ? gfit : nde_edit_target_fits(item_id);
	if (!target) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		clearfits(pre_k);
		free(pre_k);
		goto fail_free;
	}

	/* Arm the log BEFORE installing the pixels: once the pre-K state is on
	 * screen the user can start an operation, and a capture that arrived
	 * before the point was armed would append to the end instead. */
	if (insert) {
		if (!nde_history_insert_point_set(record_id, item_id, err)) {
			clearfits(pre_k);
			free(pre_k);
			goto fail_free;
		}
		/* The undo stack describes the true lineage, which is about to leave
		 * the screen; Ctrl-Z against the pre-K state would restore pixels
		 * nothing in the log describes.  (Every other history edit flushes
		 * undo too — sketch §7, there is no meta-undo.) */
		undo_flush();
	}

	/* Install: after the swap pre_k holds the TRUE pixels — that is the
	 * stash the end path restores from. */
	apv_swap_into_target(target, pre_k);

	g_mutex_lock(&apv_mutex);
	apv.active = TRUE;
	apv.installed = TRUE;
	apv.insert = insert;
	apv.borrowed = borrow;
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

gboolean nde_amend_preview_begin_execute(gint64 record_id, gchar **err) {
	return apv_begin_execute(record_id, FALSE, err);
}

gboolean nde_amend_preview_end_execute(gboolean apply, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	g_mutex_lock(&apv_mutex);
	if (!apv.active || apv.insert) {
		gboolean wrong_mode = apv.insert;
		g_mutex_unlock(&apv_mutex);
		if (apply) {
			*err = g_strdup(wrong_mode ? _("an insertion point is open — finish it instead")
			                           : _("no amend preview is active"));
			return FALSE;
		}
		return TRUE;   /* tolerated: defensive cancels from destroy handlers */
	}
	gint64 record_id = apv.record_id;
	gint item_id = apv.item_id;
	gboolean borrowed = apv.borrowed;
	fits *saved = apv.saved;
	apv_clear_state_locked();
	g_mutex_unlock(&apv_mutex);

	/* Restore the true pixels FIRST (the plan's ordering contract): gfit
	 * must hold the real image again before any amend runs, so a failed
	 * amend changes nothing.  `saved` receives the pre-K/preview pixels
	 * from the swap — superseded either way.  A borrowed item has no fits of
	 * its own: the display held its state and takes its own back here, and
	 * what carries the amend to the image is the composite recompute the
	 * edit does. */
	if (saved) {
		fits *target = borrowed ? gfit : nde_edit_target_fits(item_id);
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
		return nde_edit_execute(record_id, new_params, err);
	return TRUE;
}

/* ======================================================================= */
/* Region-scoped tail replay (roi-nde-plan.md phase 9)                     */
/* ======================================================================= */

/* Taken once per preview run.  The halo the crop is grown by and the records
 * replayed into that crop MUST come from the same list, or the tail runs with
 * less context than it was measured to need — hence a plan object rather than
 * two independent queries. */
struct nde_region_tail {
	nde_chain *chain;
	guint      from;   /* first tail member index (== len ⇒ empty tail) */
	int        halo;   /* Σ of the tail's declared halos */
};

/* Why @rec cannot be recomputed on a rectangle — a translated heap string —
 * or NULL when it can, in which case *halo_out receives its declared halo.
 *
 * The op's own description is preferred to its id for the message: this is
 * shown in the amend banner, and "Asinh transformation" is what the user
 * clicked, whereas "colors.asinh" is what we called it. */
static gchar *region_tail_member_reason(const nde_record *rec, guint8 flags,
                                        int *halo_out) {
	const nde_op_class *cls = nde_op_class_for(rec->op_id);
	const op_descriptor *op = cls->desc;
	const gchar *name = cls->label ? _(cls->label)
	                  : (rec->op_id ? rec->op_id : "?");

	if (flags & NDE_CHAIN_MEMBER_BARRIER)
		return g_strdup_printf(_("\"%s\" cannot be recomputed"), name);
	if (rec->tier != NDE_TIER_A)
		return g_strdup_printf(_("\"%s\" is not a replayable step"), name);
	/* Composites and joint records read OTHER items, at full size.  A
	 * composite would also be caught by the descriptor test below (it has
	 * none), but naming both separately is what makes the banner's reason
	 * useful. */
	if (cls->family == NDE_OPC_COMPOSITE || cls->family == NDE_OPC_JOINT)
		return g_strdup_printf(_("\"%s\" combines several images"), name);
	if (!op || !op->deserialize)
		return g_strdup_printf(_("\"%s\" is not a known operation"), name);
	if (cls->traits & NDE_OPT_GEOMETRIC)
		return g_strdup_printf(_("\"%s\" changes the image geometry"), name);
	if (!op_descriptor_is_roi_capable(op))
		return g_strdup_printf(_("\"%s\" cannot be computed on a region"), name);
	/* Parameter-completeness, in the one form the code can actually test for:
	 * replay_pre exists so a record can re-derive something from the image it
	 * is about to run on (background extraction refits from its recorded
	 * sample positions), and a crop is not that image. */
	if (op->replay_pre)
		return g_strdup_printf(_("\"%s\" re-derives its settings from the whole image"), name);

	gpointer user = op->deserialize(rec->params, rec->op_version);
	if (!user)
		return g_strdup_printf(_("\"%s\": its settings could not be read"), name);
	*halo_out = op_descriptor_roi_halo(op, user);
	destroy_any_args(user);
	return NULL;
}

/* Build the plan, or NULL.  @editing NULL asks only "what is the regime?" —
 * the banner's question — and skips the op match. */
static nde_region_tail *region_tail_plan(const op_descriptor *editing,
                                         gchar **why) {
	if (why)
		*why = NULL;

	g_mutex_lock(&apv_mutex);
	/* insert mode has no "record being edited" to sit in front of the tail:
	 * what is hidden there is [anchor..head], replayed forward over inserted
	 * work at the end, and the ops that produce that work are ordinary
	 * full-image runs.  A later step. */
	gboolean applies  = apv.installed && !apv.insert && !apv.tail_failed;
	gboolean borrowed = apv.borrowed;
	gint64   record_id = apv.record_id;
	gint     item_id   = apv.item_id;
	gboolean op_match = !editing || (editing->id && apv.op_id
	                                 && !g_strcmp0(editing->id, apv.op_id));
	g_mutex_unlock(&apv_mutex);

	if (!applies || !op_match)
		return NULL;
	if (borrowed) {
		/* A merge or a flatten consumed this item, and the refusal is about
		 * the CHAIN BOUNDARY, not about coordinates — the composite record
		 * does keep the input's canvas position (nde_composite_input), so an
		 * offset is not what is missing.
		 *
		 * What is missing is everything past the boundary.  Chain membership
		 * is target_item_id, and the composite targets the item it PRODUCED,
		 * so nde_chain_build(input) stops where the input was consumed: the
		 * composite and every step built on its result sit on another item's
		 * chain and this scan cannot see them.  It would therefore report a
		 * replayable tail — truthfully, about the input in isolation — while
		 * the user is looking at a display fed through a composite that is not
		 * being replayed at all, and the banner would promise "the whole
		 * chain".
		 *
		 * The composite could not be replayed on a rectangle anyway: it reads
		 * its other inputs at full size.  Making it able to is exactly the
		 * windowed composite (roi-nde-plan.md phase 9 items 3 and 5), which is
		 * why that is where this belongs and not in a coordinate patch. */
		if (why)
			*why = g_strdup(_("this image was merged into another, so only the "
			                  "step being edited can be previewed"));
		return NULL;
	}

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
		nde_chain_free(chain);
		return NULL;
	}

	int halo = 0;
	gchar *reason = NULL;
	for (guint i = (guint)e + 1; i < chain->records->len && !reason; i++) {
		int h = 0;
		reason = region_tail_member_reason(g_ptr_array_index(chain->records, i),
		                                   g_array_index(chain->member_flags, guint8, i),
		                                   &h);
		halo += h;
	}
	if (reason) {
		nde_chain_free(chain);
		if (why)
			*why = reason;
		else
			g_free(reason);
		return NULL;
	}

	nde_region_tail *plan = g_new0(nde_region_tail, 1);
	plan->chain = chain;
	plan->from  = (guint)e + 1;
	plan->halo  = halo;
	return plan;
}

nde_region_tail *nde_region_tail_begin(const op_descriptor *editing, int *halo) {
	nde_region_tail *plan = region_tail_plan(editing, NULL);
	if (plan && halo)
		*halo = plan->halo;
	return plan;
}

gboolean nde_region_tail_available(int *halo, gchar **why) {
	nde_region_tail *plan = region_tail_plan(NULL, why);
	if (!plan)
		return FALSE;
	if (halo)
		*halo = plan->halo;
	nde_region_tail_free(plan);
	return TRUE;
}

void nde_region_tail_free(nde_region_tail *plan) {
	if (!plan)
		return;
	nde_chain_free(plan->chain);
	g_free(plan);
}

/* nde_mask_pin_install for a crop: the pinned mask is item-sized, so it is cropped
 * to the same rectangle the pixels were.  Masking is per-pixel, so this is
 * exact — a masked record needs no halo of its own beyond the op's. */
static gboolean region_mask_pin_install(fits *region, const nde_record *rec,
                                        const rectangle *rect, gchar **err) {
	const nde_pin *pin = nde_record_input(rec, "mask");
	if (!pin)
		return TRUE;
	fits *mfit = nde_state_release(nde_checkpoint_get_at(pin->item_id,
	                                                     pin->record_id));
	if (!mfit) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): its mask is no longer stored"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}
	gboolean ok = FALSE;
	fits cropped = { 0 };
	if (rect->x >= 0 && rect->y >= 0
	    && rect->x + rect->w <= (int)mfit->rx
	    && rect->y + rect->h <= (int)mfit->ry
	    && !crop_fits_region(mfit, rect, &cropped)) {
		mask_t *m = fits_to_mask(&cropped);
		if (m) {
			region->mask = m;
			region->mask_active = TRUE;
			ok = TRUE;
		}
	}
	clearfits(&cropped);
	clearfits(mfit);
	free(mfit);
	if (!ok)
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): its mask does not fit the image"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
	return ok;
}

/* Has this run been superseded?
 *
 * Cancellation here is ROUTINE, not a fault: notify_update cancels the running
 * preview whenever a newer tick arrives (a slider drag, a moved rectangle), so
 * it fires constantly in ordinary use.  It must neither disarm the feature nor
 * say anything — the tick that superseded us is about to redraw the rectangle.
 *
 * The worker-thread test is what makes the question answerable at all.  A LIVE
 * preview runs generic_image_worker SYNCHRONOUSLY on the GTK main thread
 * (asinh.c and every other live-preview dialog), never through the job queue —
 * so nothing clears cancel_flag for it, and it still holds the 1 that the last
 * queued job's stop_processing_thread() left behind.  Polling it there reported
 * "cancelled" on every tick after the first, which silently reduced the
 * rectangle to the edited step alone.  Off the processing thread the flag says
 * nothing about us, and a tail bounded by the rectangle is short enough that
 * running it to completion is the right answer anyway. */
static gboolean region_tail_cancelled(void) {
	return processing_in_worker_thread() && !processing_should_continue();
}

gboolean nde_region_tail_apply(nde_region_tail *plan, fits *region,
                               const rectangle *rect, int max_threads) {
	g_return_val_if_fail(plan != NULL && region != NULL && rect != NULL, FALSE);

	/* The crop carries the IMAGE's own processing mask, copied by
	 * crop_fits_region.  The tail's records carry their own pinned ones, and
	 * nde_mask_pin_clear frees whatever it finds — so put the crop's aside and
	 * give it back untouched, leaving this function non-destructive of its
	 * input's mask.  (The hook for record K has already run, under that mask,
	 * before we get here.) */
	mask_t  *own_mask   = region->mask;
	gboolean own_active = region->mask_active;
	region->mask = NULL;
	region->mask_active = FALSE;

	gchar *err = NULL;
	gboolean cancelled = FALSE;
	for (guint i = plan->from; i < plan->chain->records->len && !err; i++) {
		const nde_record *rec = g_ptr_array_index(plan->chain->records, i);
		if (region_tail_cancelled()) {
			cancelled = TRUE;
			break;
		}
		/* Vetted by region_tail_member_reason when the plan was built. */
		const op_descriptor *op = op_descriptor_by_id(rec->op_id);
		gpointer user = op->deserialize(rec->params, rec->op_version);
		if (!user) {
			err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
			                      rec->record_id, rec->op_id);
			break;
		}
		if (!region_mask_pin_install(region, rec, rect, &err)) {
			destroy_any_args(user);
			break;
		}
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_any_args(user);
			nde_mask_pin_clear(region);
			err = g_strdup(_("out of memory"));
			break;
		}
		args->fit         = region;     /* PRIVATE fits, as every replay demands */
		args->op          = op;
		args->user        = user;
		args->nde_replay  = TRUE;
		args->mask_aware  = (region->mask != NULL);
		args->max_threads = max_threads;
		gint64 outer = nde_replay_current_record_id();
		nde_replay_set_current_record(rec->record_id);
		int rc = GPOINTER_TO_INT(generic_image_worker(args));
		nde_replay_set_current_record(outer);
		free_generic_img_args(args);    /* frees user through its destructor */
		nde_mask_pin_clear(region);
		/* NO nde_snapstore_deposit here, deliberately: these are region-sized
		 * intermediates and depositing one would let a later FULL replay
		 * restart from a rectangle. */
		if (rc) {
			/* A hook that polls the cancel flag reports failure when a newer
			 * tick supersedes this one; same routine cancellation, one layer
			 * down. */
			if (region_tail_cancelled())
				cancelled = TRUE;
			else
				err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) failed to apply"),
				                      rec->record_id, rec->op_id);
			break;
		}
	}

	region->mask = own_mask;
	region->mask_active = own_active;

	if (cancelled)
		return FALSE;
	if (!err)
		return TRUE;

	/* Sticky degrade: say it once, then stop trying.  The alternative is this
	 * message on every slider tick, which buries itself. */
	g_mutex_lock(&apv_mutex);
	gboolean first = !apv.tail_failed;
	apv.tail_failed = TRUE;
	g_mutex_unlock(&apv_mutex);
	if (first)
		siril_log_error(_("Region preview: the later steps could not be recomputed "
		                  "(%s); showing this step alone.\n"), err);
	g_free(err);
	return FALSE;
}

/* ======================================================================= */
/* Edit at / insert before K (graph step 2)                                */
/* ======================================================================= */

/* Same install as the amend preview, but the exit verb differs: instead of
 * re-running the anchor with new parameters, the operations the user ran
 * against the pre-anchor state are already IN the log (nde_history_append
 * inserted them before the anchor), and finishing replays the anchor and
 * everything after it forward over them.
 *
 * Note what is NOT replayed: the inserted steps themselves.  Their result is
 * the pixels sitting in the target, so their tier is irrelevant — inserting
 * an opaque python write is as legal as inserting a Tier-A op, it just
 * freezes everything before it for future edits, which is the ordinary
 * barrier rule.  What must be replayable is [anchor..end], and that is
 * checked when the point is armed, before anything changes. */

gboolean nde_edit_at_begin_execute(gint64 record_id, gchar **err) {
	return apv_begin_execute(record_id, TRUE, err);
}

gboolean nde_edit_at_end_execute(gboolean apply, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	g_mutex_lock(&apv_mutex);
	if (!apv.active || !apv.insert) {
		g_mutex_unlock(&apv_mutex);
		if (apply) {
			*err = g_strdup(_("no insertion point is open"));
			return FALSE;
		}
		return TRUE;   /* tolerated: defensive cancels from destroy handlers */
	}
	gint64 anchor_id = apv.record_id;
	gint item_id = apv.item_id;
	gboolean borrowed = apv.borrowed;
	fits *saved = apv.saved;
	apv_clear_state_locked();
	g_mutex_unlock(&apv_mutex);

	gboolean disturbed = nde_history_insert_point_disturbed();
	GArray *inserted = nde_history_insert_point_clear();

	/* Restore the true pixels and metadata first, exactly as the amend
	 * preview does.  `saved` comes back holding the pre-anchor state PLUS
	 * whatever was inserted — the starting point for the forward replay. */
	fits *target = borrowed ? gfit : nde_edit_target_fits(item_id);
	if (target)
		apv_swap_into_target(target, saved);
	else
		siril_log_error(_("Edit at: the target layer no longer exists\n"));

	gboolean commit = apply && target && inserted->len > 0 && !disturbed;
	if (!commit) {
		clearfits(saved);
		free(saved);
		/* Nothing is committed, so the inserted records describe pixels that
		 * no longer exist anywhere.  Drop them. */
		nde_history_drop_records(inserted);
		gboolean abandoned = apply && disturbed && inserted->len > 0;
		g_array_unref(inserted);
		undo_flush();
		gui_iface.set_progress(PROGRESS_RESET,
		                       abandoned ? _("Insertion abandoned")
		                                 : _("Insertion point closed"));
		if (abandoned) {
			*err = g_strdup(_("an operation changed the document while the "
			                  "insertion point was open — the inserted steps "
			                  "were discarded"));
			return FALSE;
		}
		return TRUE;
	}

	gint64 first_inserted = g_array_index(inserted, gint64, 0);
	/* The last step of this item the insertion left alone, for the reverse
	 * mask cascade below.  The inserted records are already in the log, so
	 * the predecessor of the first of them is exactly that. */
	const gint64 unchanged_upto = nde_log_predecessor(item_id, first_inserted);

	/* Replay the anchor and everything after it over the inserted work. */
	nde_chain *chain = nde_chain_build(item_id);
	gint k = -1;
	for (guint i = 0; i < chain->records->len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (rec->record_id == anchor_id) {
			k = (gint)i;
			break;
		}
	}
	nde_state *result = NULL;
	if (k < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is no longer part of this image's history"),
		                       anchor_id);
	} else if ((guint)k < chain->tail_start || !chain->tail_replayable) {
		/* An inserted opaque step became a barrier that swallowed the anchor
		 * (it has no output checkpoint, or one of the later steps cannot
		 * replay from it). */
		*err = g_strdup(_("the steps after the insertion point can no longer be recomputed"));
	} else if (borrowed) {
		/* The item has no layer to receive a result: what the insertion
		 * changed reaches the image only through the composites that consumed
		 * it, exactly as an AMEND on the same item does.  The display already
		 * holds its true pixels again, and the cascade overwrites them with
		 * the recomputed merge. */
		gui_iface.set_progress(0.f, _("Recomputing edit history..."));
		nde_snapstore_invalidate_from(item_id, first_inserted);
		clearfits(saved);
		free(saved);
		saved = NULL;
		/* Joint records after the insertion read this item's changed prefix —
		 * their other participants need recomputing too. */
		GArray *jt = nde_edit_joint_targets(item_id, first_inserted, FALSE);
		nde_cascade_from(item_id, unchanged_upto, jt, FALSE);
		g_array_unref(jt);
		g_array_unref(inserted);
		nde_edit_finish(NULL, _("Insertion applied"));
		return TRUE;
	} else {
		gui_iface.set_progress(0.f, _("Recomputing edit history..."));
		/* Cached states at or after the insertion describe the pre-insert
		 * lineage; the replay re-deposits fresh ones as it goes. */
		nde_snapstore_invalidate_from(item_id, first_inserted);
		/* A positionless state: a chain that moves the layer is refused an
		 * insertion point in the first place (apv_begin_execute), so there is
		 * no position to restart from — and the layer's current one is not it,
		 * embodying as it does every record after the anchor. */
		result = nde_replay_apply_records(nde_state_new(saved), chain, (guint)k,
		                                  chain->records->len, err);
		saved = NULL;   /* consumed either way */
	}
	nde_chain_free(chain);

	if (!result) {
		if (saved) {
			clearfits(saved);
			free(saved);
		}
		nde_snapstore_invalidate_from(item_id, first_inserted);
		/* The target already holds the true pixels; take the log back to
		 * match them so the two never disagree. */
		nde_history_drop_records(inserted);
		g_array_unref(inserted);
		undo_flush();
		gui_iface.set_progress(PROGRESS_RESET, _("Insertion failed — nothing was changed"));
		return FALSE;
	}
	g_array_unref(inserted);

	/* No log commit: the inserted records were written to the log as they ran,
	 * so it already says what these pixels are. */
	nde_commit_ctx cc = {
		.item_id = item_id,
		.target  = target,
		.result  = result,
	};
	if (!nde_commit_replayed(&cc, NULL, NULL, err)) {
		gui_iface.set_progress(PROGRESS_RESET, _("Insertion failed — nothing was changed"));
		return FALSE;
	}

	/* A mask built from these pixels read a prefix the inserted steps now sit
	 * inside; joint records after the insertion read the changed prefix. */
	GArray *jt = nde_edit_joint_targets(item_id, first_inserted, FALSE);
	nde_cascade_from(item_id, unchanged_upto, jt, FALSE);
	g_array_unref(jt);

	nde_edit_finish(target, _("Edit history updated"));
	return TRUE;
}

gboolean nde_edit_at_refuses_op(const char *what) {
	g_mutex_lock(&apv_mutex);
	gboolean open = apv.active && apv.insert;
	g_mutex_unlock(&apv_mutex);
	if (!open)
		return FALSE;
	siril_log_error(_("%s cannot be done while a step is being inserted into the "
	                  "edit history — finish or cancel the insertion first\n"),
	                what ? what : _("This operation"));
	return TRUE;
}

gboolean nde_edit_at_active(void) {
	g_mutex_lock(&apv_mutex);
	gboolean r = apv.active && apv.insert;
	g_mutex_unlock(&apv_mutex);
	return r;
}

gint64 nde_edit_at_record_id(void) {
	g_mutex_lock(&apv_mutex);
	gint64 id = (apv.active && apv.insert) ? apv.record_id : 0;
	g_mutex_unlock(&apv_mutex);
	return id;
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

static gpointer edit_at_begin_worker(gpointer p) {
	gint64 record_id = *(gint64 *)p;
	g_free(p);
	gchar *errmsg = NULL;
	gboolean ok = nde_edit_at_begin_execute(record_id, &errmsg);
	if (!ok)
		siril_log_error(_("Cannot insert before this history step: %s\n"),
		                errmsg ? errmsg : "?");
	g_free(errmsg);
	siril_add_idle(apv_ready_idle, GINT_TO_POINTER(ok ? 0 : 1));
	return GINT_TO_POINTER(ok ? 0 : 1);
}

gboolean nde_edit_at_start(gint64 record_id,
                           nde_amend_preview_ready_fn on_ready, gpointer user) {
	if (nde_amend_preview_active()) {
		siril_log_error(_("Another history step is already being edited\n"));
		return FALSE;
	}
	apv.on_ready = on_ready;
	apv.on_ready_user = user;
	gint64 *id = g_new(gint64, 1);
	*id = record_id;
	if (!replay_conductor_start(edit_at_begin_worker, id)) {
		apv.on_ready = NULL;
		apv.on_ready_user = NULL;
		g_free(id);
		return FALSE;
	}
	return TRUE;
}

static gpointer edit_at_end_worker(gpointer p) {
	gboolean apply = (GPOINTER_TO_INT(p) != 0);
	gchar *errmsg = NULL;
	gboolean ok = nde_edit_at_end_execute(apply, &errmsg);
	if (!ok)
		siril_log_error(_("Edit history change failed: %s\n"), errmsg ? errmsg : "?");
	g_free(errmsg);
	siril_add_idle(nde_edit_done_idle, NULL);
	return GINT_TO_POINTER(ok ? 0 : 1);
}

gboolean nde_edit_at_end(gboolean apply) {
	if (!nde_edit_at_active())
		return !apply;
	return replay_conductor_start(edit_at_end_worker, GINT_TO_POINTER(apply ? 1 : 0));
}
