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
 * NDE history editing: amend, delete, reorder, and undoing a composite.
 *
 * An edit changes what the log says, replays the item from a restart point to
 * get the pixels the new log describes, commits both, and then chases down
 * everything that was derived from the old ones.  The replay itself is
 * nde_replay.c's; what is here is deciding what to replay, committing the
 * result, and the cascade.
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

static void cascade_joint_targets(GArray *targets);
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/siril_pythonmodule.h"

#include "core/nde_replay_internal.h"

static void cascade_joint_targets(GArray *targets);

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
nde_state *nde_edit_restart_state(const nde_chain *chain, guint e,
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
		nde_state *state = nde_snapstore_lookup(chain->item_id, pre_id, FALSE);
		if (!state) {
			const nde_record *prev = g_ptr_array_index(chain->records, j - 1);
			state = nde_snapstore_lookup(chain->item_id, prev->record_id, TRUE);
		}
		if (state) {
			/* A cached entry is a whole value, position included, so a chain
			 * that moves the layer can restart from one — it could not while
			 * the cache held pixels alone and geometry had nothing to anchor
			 * to.  But only from one that DOES say where the layer was: undo's
			 * PRE-states carry no position, and neither does anything
			 * deposited before this chain gained a geometry member.  Restarting
			 * from those would anchor the moves that follow to the baseline
			 * position, which is where the layer was before the FIRST op, not
			 * before this one.  Keep looking further back. */
			if (chain->has_geometry && !state->has_pos) {
				nde_state_free(state);
				continue;
			}
			*start_idx = j;
			return state;
		}
	}
	/* Fall back to the phase-4 restart point. */
	nde_state *start = nde_replay_chain_restart_state(chain, err);
	if (!start && *err)
		return NULL;
	*start_idx = chain->tail_start;
	return start;
}

/* Install a rebuilt mask into whatever slot the mask item names.  Takes the
 * mask out of @built (which the caller frees).  FALSE + @err when the slot no
 * longer exists — a mask can be cleared or routed away between the edit
 * starting and finishing. */
static gboolean commit_mask_value(gint item_id, fits *built, gchar **err) {
	if (!built->mask || !built->mask->data) {
		*err = g_strdup(_("the rebuilt mask is empty"));
		return FALSE;
	}
	flis_layer_t *owner = NULL;
	flis_item_kind kind = (item_id == NDE_ITEM_PLAIN_MASK) ?
			FLIS_ITEM_MASK : flis_item_lookup(item_id, &owner);

	if (item_id == NDE_ITEM_PLAIN_MASK || kind == FLIS_ITEM_MASK) {
		fits *host = owner ? owner->fit : gfit;
		if (!host) {
			*err = g_strdup(_("the mask's image no longer exists"));
			return FALSE;
		}
		if (host->mask)
			free_mask(host->mask);
		host->mask = built->mask;    /* moved */
		built->mask = NULL;
		host->mask_active = TRUE;
		return TRUE;
	}
	if (kind == FLIS_ITEM_LMASK && owner) {
		/* A layer mask is the same pixels in a different wrapper. */
		layermask_t *lm = calloc(1, sizeof(layermask_t));
		if (!lm) {
			*err = g_strdup(_("out of memory"));
			return FALSE;
		}
		lm->w      = built->rx;
		lm->h      = built->ry;
		lm->bitpix = built->mask->bitpix;
		lm->data   = built->mask->data;
		built->mask->data = NULL;
		free_mask(built->mask);
		built->mask = NULL;
		if (flis_layer_set_lmask(owner, lm)) {
			layermask_free(lm);
			*err = g_strdup(_("the rebuilt mask does not fit its layer"));
			return FALSE;
		}
		gui_iface.flis_invalidate_composite();
		return TRUE;
	}
	*err = g_strdup(_("the edited mask no longer exists"));
	return FALSE;
}

/* ---- reverse invalidation (graph step 6) --------------------------------
 * An input pin points one way, from consumer to source, because that is the
 * direction provenance is captured in.  Propagation runs the other way: when
 * a source changes, everything pinned to it is stale.  There is no reverse
 * index to maintain — the pins are few and the log is short, so the consumers
 * are found by scanning it, which cannot go out of date.                    */

/* Recompute an item's pixels from its (unchanged) log and commit them.  For
 * when something the item DEPENDS ON changed: the log is already right, only
 * the pixels are stale. */
static gboolean recompute_item(gint item_id, gchar **err) {
	nde_chain *chain = nde_chain_build(item_id);
	if (!chain->replayable) {
		GString *m = g_string_new(NULL);
		for (guint i = 0; i < chain->reasons->len; i++) {
			if (i)
				g_string_append(m, "; ");
			g_string_append(m, g_ptr_array_index(chain->reasons, i));
		}
		*err = g_string_free(m, FALSE);
		nde_chain_free(chain);
		return FALSE;
	}
	/* An item born of a merge or flatten has no baseline: its first record IS
	 * its origin, and the replay starts by rendering that composite's inputs
	 * (nde_chain_replay does the same). */
	nde_state *start = chain->from_composite ?
			NULL : nde_checkpoint_baseline_get(item_id);
	if (!start && !chain->from_composite) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		nde_chain_free(chain);
		return FALSE;
	}
	nde_replay_arm_position(start, chain);
	nde_snapstore_invalidate_from(item_id, 0);
	nde_state *result = nde_replay_apply_records(start, chain, 0,
	                                             chain->records->len, err);
	nde_chain_free(chain);
	if (!result)
		return FALSE;

	fits *target = nde_edit_target_fits(item_id);
	if (!target) {
		*err = g_strdup(_("the target layer no longer exists"));
		nde_state_free(result);
		return FALSE;
	}
	gboolean quiesced = nde_commit_lock(target);
	nde_commit_pixels(target, result->pix);
	nde_commit_unlock(target, quiesced);
	nde_commit_restore_metadata(target, result->pix);
	nde_commit_layer_position(item_id, result);
	nde_state_free(result);
	invalidate_stats_from_fit(target);
	return TRUE;
}

/* TRUE when @item_id names nothing in the document any more but IS still an
 * input to a composite — a RETAINED INPUT (design note §5.1).  Its records are
 * live and editable; what an edit to them recomputes is not its own value,
 * which has nowhere to go, but every composite downstream of it.
 *
 * "Names nothing" rather than "has no layer": a layer consumed by a composite
 * takes its layer mask with it, and that mask item is retained on exactly the
 * same terms — the composite pins it too. */
gboolean nde_item_is_retained_input(gint item_id) {
	if (item_id < 0 || flis_item_lookup(item_id, NULL) != FLIS_ITEM_NONE)
		return FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	gboolean found = FALSE;
	for (guint i = 0; !found && live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (!nde_composite_record_replayable(rec))
			continue;
		/* The composite's own target still has a layer to commit into, so it
		 * is not what makes an input retained. */
		found = rec->target_item_id != item_id &&
		        nde_record_input_by_item(rec, item_id) != NULL;
	}
	if (live)
		g_ptr_array_unref(live);
	return found;
}

/* Propagate a change to @item_id across composite edges: recompute every item
 * whose chain composites it in.  Unlike the mask cascade there is no stored
 * copy to refresh first — a composite resolves its inputs by replay, so the
 * consumer's own recompute picks the new state up (nde_composite.h).  A
 * consumer that cannot be recomputed is named rather than left silently
 * stale. */
static void cascade_composite_consumers(gint item_id) {
	GPtrArray *live = nde_history_snapshot(NULL);
	GHashTable *seen = g_hash_table_new(g_direct_hash, g_direct_equal);
	GQueue *frontier = g_queue_new();
	g_queue_push_tail(frontier, GINT_TO_POINTER(item_id));
	g_hash_table_add(seen, GINT_TO_POINTER(item_id));

	guint redone = 0;
	while (!g_queue_is_empty(frontier)) {
		gint src = GPOINTER_TO_INT(g_queue_pop_head(frontier));
		for (guint i = 0; live && i < live->len; i++) {
			const nde_record *rec = g_ptr_array_index(live, i);
			if (rec->target_item_id == src ||
			    !nde_composite_record_replayable(rec) ||
			    !nde_record_input_by_item(rec, src))
				continue;
			gint consumer = rec->target_item_id;
			if (!g_hash_table_add(seen, GINT_TO_POINTER(consumer)))
				continue;
			/* The composite reads its inputs by replay, so anything cached for
			 * the consumer at or after this record describes the OLD input. */
			nde_snapstore_invalidate_from(consumer, 0);
			/* A consumer that was itself composited away has no layer to commit
			 * into: nothing shows its pixels, and what has to be recomputed is
			 * whatever consumed IT.  Merging a group is exactly this — a run of
			 * merge-downs, each one's survivor consumed by the next. */
			if (nde_item_is_retained_input(consumer)) {
				g_queue_push_tail(frontier, GINT_TO_POINTER(consumer));
				continue;
			}
			gchar *cerr = NULL;
			if (recompute_item(consumer, &cerr)) {
				redone++;
			} else {
				siril_log_warning(_("The layer this was merged into could not be recomputed, so it still shows the old result: %s\n"),
				                  cerr ? cerr : "?");
			}
			g_free(cerr);
		}
	}
	g_queue_free(frontier);
	g_hash_table_destroy(seen);
	if (live)
		g_ptr_array_unref(live);
	if (redone) {
		gui_iface.flis_invalidate_composite();
		siril_log_info(_("Change applied to %u merged image(s)\n"), redone);
	}
}

/* Delete @record_id, the one step that built the layer mask @mask_item, when a
 * composite is the only thing still holding that mask.  There is no mask left
 * to commit anywhere, so what changes is the composites: each stops applying
 * it, and each is recomputed.
 *
 * The stored copy needs no explicit drop — it lives at the deleted record's
 * own coordinate, so nde_history_delete takes it with the record. */
static gboolean delete_retained_mask(gint64 record_id, gint mask_item, gchar **err) {
	GArray *layers = g_array_new(FALSE, FALSE, sizeof(gint));   /* masked inputs */
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (!nde_composite_is_op(rec->op_id) ||
		    !nde_record_input_by_item(rec, mask_item))
			continue;
		nde_composite_state *st = nde_composite_state_parse(rec->params);
		if (!st)
			continue;
		for (guint k = 0; k < st->inputs->len; k++) {
			const nde_composite_input *in =
					&g_array_index(st->inputs, nde_composite_input, k);
			if (in->mask_item_id != mask_item)
				continue;
			if (nde_history_drop_mask_input(rec->record_id, mask_item)) {
				siril_log_message(_("The layer mask no longer applies to '%s' in "
				                    "step %" G_GINT64_FORMAT "\n"),
				                  in->name ? in->name : "?", rec->record_id);
				g_array_append_val(layers, in->item_id);
			}
		}
		nde_composite_state_free(st);
	}
	if (live)
		g_ptr_array_unref(live);

	if (!nde_history_delete(record_id, err)) {
		g_array_free(layers, TRUE);
		return FALSE;
	}
	/* Recompute through the LAYER the mask was applied to: that is the pin the
	 * composites still carry, and it reaches the retained and transitive cases
	 * the same way any other upstream change does. */
	for (guint i = 0; i < layers->len; i++)
		cascade_composite_consumers(g_array_index(layers, gint, i));
	g_array_free(layers, TRUE);
	nde_edit_finish(NULL, _("Edit history updated"));
	return TRUE;
}

/* Re-derive @mask_item's mask as it stood after chain member @upto-1 and
 * replace the stored copy its consumers read (nde_replay.h's pin store). */
static gboolean refresh_pinned_mask(const nde_chain *mask_chain, guint upto,
                                    gint64 coord_record, gint mask_item,
                                    gchar **err) {
	fits *built = nde_mask_chain_replay(mask_chain, upto, err);
	if (!built)
		return FALSE;
	gboolean ok = FALSE;
	if (built->mask && built->mask->data) {
		fits *mfit = mask_to_fits(built);
		if (mfit) {
			nde_checkpoint_output_store(mfit, coord_record, mask_item);
			clearfits(mfit);
			free(mfit);
			ok = TRUE;
		}
	}
	if (!ok)
		*err = g_strdup(_("the rebuilt mask could not be stored"));
	clearfits(built);
	free(built);
	return ok;
}

/* Propagate a change to @mask_item from chain position @from_pos onwards.
 * Refreshes every stored pin state that moved, then recomputes each item that
 * consumed one.  Reports what it did; a consumer that cannot be recomputed is
 * named rather than silently left stale. */
static void cascade_mask_consumers(gint mask_item, guint from_pos) {
	nde_chain *mask_chain = nde_chain_build(mask_item);
	if (!mask_chain->replayable || !mask_chain->records->len) {
		nde_chain_free(mask_chain);
		return;
	}
	/* Position of each of the mask's records, so "at or after the edit" is a
	 * comparison rather than a guess. */
	GHashTable *pos_of = g_hash_table_new_full(g_int64_hash, g_int64_equal, g_free, NULL);
	for (guint i = 0; i < mask_chain->records->len; i++) {
		const nde_record *r = g_ptr_array_index(mask_chain->records, i);
		gint64 *k = g_new(gint64, 1);
		*k = r->record_id;
		g_hash_table_insert(pos_of, k, GUINT_TO_POINTER(i + 1));
	}

	GHashTable *items = g_hash_table_new(g_direct_hash, g_direct_equal);
	GHashTable *done_coords = g_hash_table_new_full(g_int64_hash, g_int64_equal, g_free, NULL);
	guint refreshed = 0;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		/* Matched by SOURCE, not by role: an op pins its mask as "mask", a
		 * composite pins one per masked input as "mask0", "mask1", …  Both are
		 * stored copies at the same kind of coordinate, so both refresh here. */
		const nde_pin *pin = nde_record_input_by_item(rec, mask_item);
		if (!pin)
			continue;
		gpointer p = g_hash_table_lookup(pos_of, &pin->record_id);
		if (!p)
			continue;   /* pinned to a record that is no longer in the chain */
		guint upto = GPOINTER_TO_UINT(p);
		if (upto <= from_pos)
			continue;   /* the edit is after this pin: its mask is unchanged */
		if (!g_hash_table_contains(done_coords, &pin->record_id)) {
			gchar *err = NULL;
			if (refresh_pinned_mask(mask_chain, upto, pin->record_id,
			                        mask_item, &err)) {
				refreshed++;
				gint64 *k = g_new(gint64, 1);
				*k = pin->record_id;
				g_hash_table_insert(done_coords, k, GINT_TO_POINTER(1));
			} else {
				siril_log_warning(_("Could not rebuild the mask used by step %" G_GINT64_FORMAT ": %s\n"),
				                  rec->record_id, err ? err : "?");
				g_free(err);
				continue;
			}
		}
		g_hash_table_add(items, GINT_TO_POINTER(rec->target_item_id));
	}
	if (live)
		g_ptr_array_unref(live);
	g_hash_table_destroy(pos_of);
	g_hash_table_destroy(done_coords);

	guint redone = 0, stale = 0;
	GHashTableIter it;
	gpointer k, v;
	g_hash_table_iter_init(&it, items);
	while (g_hash_table_iter_next(&it, &k, &v)) {
		gint item = GPOINTER_TO_INT(k);
		gchar *err = NULL;
		/* A consumer that was itself composited away has no pixels of its own
		 * to write; what shows the change is whatever consumed IT. */
		if (nde_item_is_retained_input(item)) {
			cascade_composite_consumers(item);
			continue;
		}
		if (recompute_item(item, &err)) {
			redone++;
			/* The host's pixels changed, so a joint record it participates
			 * in derives new factors for its siblings (nde_joint.h). */
			GArray *jt = nde_edit_joint_targets(item, 0, FALSE);
			cascade_joint_targets(jt);
			g_array_unref(jt);
		} else {
			stale++;
			siril_log_warning(_("Item %d used this mask but could not be recomputed, so it still shows the old result: %s\n"),
			                  item, err ? err : "?");
		}
		g_free(err);
	}
	g_hash_table_destroy(items);
	nde_chain_free(mask_chain);
	if (redone || stale)
		siril_log_info(_("Mask change applied to %u earlier step(s) and %u image(s)\n"),
		               refreshed, redone);
}

/* Position of @record_id among @item_id's OWN records, in log order, or -1
 * when the log no longer holds it.  The log, not the pixel chain: a pin names
 * a record in the item's history, which is the wider set (nde_history.h).
 * Position, never id — insertion hands ids out of order. */
static gint log_position_of(const GPtrArray *live, gint item_id, gint64 record_id) {
	gint pos = 0;
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *r = g_ptr_array_index(live, i);
		if (r->target_item_id != item_id)
			continue;
		if (r->record_id == record_id)
			return pos;
		pos++;
	}
	return -1;
}

/* The record just BEFORE @record_id in @item_id's log, or 0 when it is the
 * first — "the last step of this item an edit at @record_id leaves alone".
 * Call it BEFORE the mutation for anything that moves records about. */
gint64 nde_log_predecessor(gint item_id, gint64 record_id) {
	GPtrArray *live = nde_history_snapshot(NULL);
	gint64 prev = 0;
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *r = g_ptr_array_index(live, i);
		if (r->target_item_id != item_id)
			continue;
		if (r->record_id == record_id)
			break;
		prev = r->record_id;
	}
	if (live)
		g_ptr_array_unref(live);
	return prev;
}

/* THE OTHER DIRECTION of cascade_mask_consumers.  A mask built from an image
 * — mask.from_stars, mask.from_channel and the rest — pins the exact point in
 * that image's history its pixels came from, and nde_mask_chain_replay resolves
 * the pin live.  So editing the MASK re-derives it from the right place; what
 * had no answer at all was editing the IMAGE.  The mask then went on
 * describing stars, or a channel, of pixels that no longer existed anywhere,
 * and so did every step that had used it.
 *
 * @unchanged_upto is the last record of @item_id's log the edit left alone
 * (0 = none of it — the change reaches back to the baseline).  A mask pinned
 * at or before that point read the same pixels as before and is not touched:
 * re-deriving it would give identical results, but only after recomputing
 * every image that consumed it, which is a visible cost for no change.
 *
 * Runs AFTER the log and pixel commits — it re-derives by replaying the image
 * from the amended log. */
static void cascade_derived_masks(gint item_id, gint64 unchanged_upto) {
	GPtrArray *live = nde_history_snapshot(NULL);
	if (!live)
		return;
	/* -1 = "nothing of this item is unchanged", which is what a missing
	 * anchor should also mean: rebuilding too much is a cost, rebuilding too
	 * little is a lie. */
	const gint safe_pos = unchanged_upto ?
			log_position_of(live, item_id, unchanged_upto) : -1;

	GHashTable *seen = g_hash_table_new(g_direct_hash, g_direct_equal);
	guint rebuilt = 0, stale = 0;
	for (guint i = 0; i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		/* An "image" pin is what a mask records its origin as, and nothing
		 * else uses that role (processing.c).  Only a mask chain's FIRST
		 * record carries one — an edit of an existing mask reads the mask. */
		const nde_pin *pin = nde_record_input(rec, "image");
		if (!pin || pin->item_id != item_id)
			continue;
		const gint mask_item = rec->target_item_id;
		if (!pin->record_id)
			continue;   /* pinned to the BASELINE, which no edit can move */
		const gint pinned = log_position_of(live, item_id, pin->record_id);
		if (pinned >= 0 && pinned <= safe_pos)
			continue;   /* the edit lands after the state this mask read */
		/* Claimed only once a record has actually asked for the rebuild.  A
		 * mask reads the image more than once when it was rebuilt partway
		 * through its own history, and marking it seen at the first pin let
		 * an unaffected early read suppress an affected later one. */
		if (!g_hash_table_add(seen, GINT_TO_POINTER(mask_item)))
			continue;

		gchar *err = NULL;
		nde_chain *mc = nde_chain_build(mask_item);
		fits *built = (mc->is_mask && mc->replayable && mc->records->len)
				? nde_mask_chain_replay(mc, mc->records->len, &err) : NULL;
		nde_chain_free(mc);
		if (!built) {
			stale++;
			siril_log_warning(_("The mask built from this image could not be "
			                    "rebuilt, so it still describes the old pixels: %s\n"),
			                  err ? err : "?");
			g_free(err);
			continue;
		}
		/* Same three homes as an edit of the mask's own history: a live slot,
		 * or none at all when the mask is dormant or held only by a composite
		 * — and then what shows the change is the cascade below. */
		const gboolean dormant = mask_item != NDE_ITEM_PLAIN_MASK &&
				flis_item_lookup(mask_item, NULL) == FLIS_ITEM_NONE;
		if (!nde_item_is_retained_input(mask_item) && !dormant &&
		    !commit_mask_value(mask_item, built, &err)) {
			stale++;
			siril_log_warning(_("The mask built from this image was rebuilt but "
			                    "could not be stored: %s\n"), err ? err : "?");
			g_free(err);
		} else {
			rebuilt++;
		}
		clearfits(built);
		free(built);
		/* Every step that ran under this mask holds a stored copy of it at a
		 * pinned coordinate; those copies are now the old mask.  from_pos 0:
		 * the mask changed from its first step, because its INPUT did. */
		cascade_mask_consumers(mask_item, 0);
	}
	g_hash_table_destroy(seen);
	g_ptr_array_unref(live);
	if (rebuilt || stale)
		siril_log_info(_("Image change applied to %u mask(s) built from it\n"),
		               rebuilt);
}

/* ---- joint-record cascade (nde_joint.h) ---------------------------------
 * A joint record scales EVERY participant from one shared analysis, so a
 * change to any participant — its pixels upstream of the record, or its tint
 * — changes what the record derives for all the OTHERS.  The cascade
 * recomputes them from their (unchanged) logs: the factors come from the
 * generation cache, so however many layers replay, the analysis runs once.
 *
 * Collection happens BEFORE the log commit (the edited record's position is
 * still measurable, and a deleted joint record's participants are still
 * readable), execution AFTER it (the recomputes must read the amended log).*/

/* Every OTHER participant of every live joint record naming @item_id at or
 * after @from_record_id's log position, closed transitively over joint edges
 * (a rescaled sibling may itself feed a later joint record).  A
 * @from_record_id of 0 — or one the log no longer holds — means "anywhere".
 * @include_self additionally lists @item_id itself when a qualifying record
 * names it: for the compositing-state (tint) edits, whose cheap path never
 * replays the item's own chain, though the joint member on it now derives a
 * different factor.  Caller g_array_unref()s. */
GArray *nde_edit_joint_targets(gint item_id, gint64 from_record_id,
                                     gboolean include_self) {
	GArray *out = g_array_new(FALSE, FALSE, sizeof(gint));
	GPtrArray *live = nde_history_snapshot(NULL);
	if (!live)
		return out;
	gint from_pos = -1;
	if (from_record_id) {
		for (guint i = 0; i < live->len; i++) {
			const nde_record *r = g_ptr_array_index(live, i);
			if (r->record_id == from_record_id) {
				from_pos = (gint)i;
				break;
			}
		}
	}
	GHashTable *seen = g_hash_table_new(g_direct_hash, g_direct_equal);
	GQueue *frontier = g_queue_new();
	g_hash_table_add(seen, GINT_TO_POINTER(item_id));
	g_queue_push_tail(frontier, GINT_TO_POINTER(item_id));
	gboolean self_named = FALSE;
	while (!g_queue_is_empty(frontier)) {
		gint src = GPOINTER_TO_INT(g_queue_pop_head(frontier));
		for (guint i = 0; i < live->len; i++) {
			const nde_record *rec = g_ptr_array_index(live, i);
			if (!nde_joint_is_op(rec->op_id))
				continue;
			/* The positional filter applies only to the edited item: a joint
			 * record BEFORE the edit read a prefix the edit does not touch.
			 * Propagated items carry no position — their whole value moved. */
			if (src == item_id && from_pos >= 0 && (gint)i < from_pos)
				continue;
			if (!nde_joint_record_names_item(rec, src))
				continue;
			if (src == item_id)
				self_named = TRUE;
			guint n = 0;
			gint *parts = nde_joint_record_participants(rec, &n);
			for (guint k = 0; parts && k < n; k++) {
				if (parts[k] == src)
					continue;
				if (!g_hash_table_add(seen, GINT_TO_POINTER(parts[k])))
					continue;
				g_array_append_val(out, parts[k]);
				g_queue_push_tail(frontier, GINT_TO_POINTER(parts[k]));
			}
			g_free(parts);
		}
	}
	if (include_self && self_named)
		g_array_append_val(out, item_id);
	g_queue_free(frontier);
	g_hash_table_destroy(seen);
	g_ptr_array_unref(live);
	return out;
}

/* A subset amend ADDS participants the pre-commit collection could not see
 * (they were not named by any live record then): merge a post-commit
 * collection into @targets, deduplicated. */
static void joint_merge_post_commit_targets(GArray *targets, gint item_id,
                                            gint64 record_id) {
	GArray *fresh = nde_edit_joint_targets(item_id, record_id, FALSE);
	for (guint i = 0; i < fresh->len; i++) {
		gint p = g_array_index(fresh, gint, i);
		gboolean seen = FALSE;
		for (guint j = 0; j < targets->len && !seen; j++)
			seen = g_array_index(targets, gint, j) == p;
		if (!seen)
			g_array_append_val(targets, p);
	}
	g_array_unref(fresh);
}

/* Guard against mutual recursion: a recomputed participant's derived-mask
 * and composite cascades can wind back into another joint cascade, and the
 * generation cache already guarantees a nested pass would derive the same
 * factors — so a nested cascade recomputes nothing new and is skipped. */
static gboolean joint_cascading = FALSE;

static void cascade_joint_targets(GArray *targets) {
	if (!targets || !targets->len || joint_cascading)
		return;
	joint_cascading = TRUE;
	guint redone = 0;
	for (guint i = 0; i < targets->len; i++) {
		gint p = g_array_index(targets, gint, i);
		/* Everything cached for the participant may embody the old factor. */
		nde_snapstore_invalidate_from(p, 0);
		if (nde_item_is_retained_input(p)) {
			cascade_composite_consumers(p);
			redone++;
			continue;
		}
		gchar *cerr = NULL;
		if (recompute_item(p, &cerr)) {
			redone++;
			cascade_derived_masks(p, 0);
			cascade_composite_consumers(p);
		} else {
			siril_log_warning(_("A jointly calibrated layer could not be recomputed, so it still shows the old scaling: %s\n"),
			                  cerr ? cerr : "?");
		}
		g_free(cerr);
	}
	joint_cascading = FALSE;
	if (redone) {
		gui_iface.flis_invalidate_composite();
		siril_log_info(_("Change applied to %u jointly calibrated layer(s)\n"), redone);
	}
}

/* ---- the shared tail of every history edit ------------------------------ */
/*
 * An edit that changed @item_id's pixels has three consequences, and they are
 * taken in this order because each reads the results of the one before:
 *
 *   1. masks DERIVED from those pixels are rebuilt (they pinned a state the
 *      edit moved),
 *   2. composites that CONSUMED the item are recomputed (they resolve their
 *      inputs by replay, so this picks the new pixels up),
 *   3. joint records the item participates in re-derive their siblings'
 *      factors (its contribution changed).
 *
 * @unchanged_upto is the last of the item's own records the edit left alone,
 * measured BEFORE the log was touched.  @joint_targets is collected before the
 * commit too — the participants of a deleted joint record are gone from the
 * log afterwards — and stays the caller's to free.
 *
 * Runs AFTER the pixel and log commits: every step of it recomputes from the
 * committed log.
 */
void nde_cascade_from(gint item_id, gint64 unchanged_upto,
                         GArray *joint_targets, gboolean joint_reanchor) {
	cascade_derived_masks(item_id, unchanged_upto);
	cascade_composite_consumers(item_id);
	nde_replay_set_joint_reanchor(joint_reanchor);
	cascade_joint_targets(joint_targets);
	nde_replay_set_joint_reanchor(FALSE);
}

/* Commit a replayed value: pixels first, then the log, and put the pixels back
 * if the log refuses.  @log_commit is what the caller's kind of edit does to
 * the history — amend, delete, reorder — or NULL when the log already says
 * what these pixels are (an insertion's records are written as they run).
 *
 * On failure nothing is changed and @err is set.  Success leaves the cascade
 * (nde_cascade_from) and the display (nde_edit_finish) to the caller, because what to
 * cascade was collected before the commit and what to say afterwards differs.
 */
gboolean nde_commit_replayed(nde_commit_ctx *c,
                                gboolean (*log_commit)(gpointer, gchar **),
                                gpointer log_user, gchar **err) {
	if (c->retained) {
		/* A retained input has no layer to swap into: the replay was run to
		 * prove the edited chain still applies, and its pixels are an
		 * intermediate value only the composite consumes.  Commit the log and
		 * let the cascade recompute the consumers — they resolve this item by
		 * replaying it again, now from the amended log. */
		g_clear_pointer(&c->result, nde_state_free);
		return !log_commit || log_commit(log_user, err);
	}
	if (!c->target) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		g_clear_pointer(&c->result, nde_state_free);
		return FALSE;
	}

	/* Atomic commit: swap pixels, then the log.  `result` holds the OLD
	 * pixels after the swap, so a log-commit failure can restore them. */
	gboolean quiesced = nde_commit_lock(c->target);
	nde_commit_pixels(c->target, c->result->pix);
	nde_commit_unlock(c->target, quiesced);

	if (log_commit && !log_commit(log_user, err)) {
		/* Should be unreachable (everything was validated, we own the
		 * slot); restore the old pixels so nothing is half-committed. */
		quiesced = nde_commit_lock(c->target);
		nde_commit_pixels(c->target, c->result->pix);
		nde_commit_unlock(c->target, quiesced);
		g_clear_pointer(&c->result, nde_state_free);
		return FALSE;
	}
	nde_commit_restore_metadata(c->target, c->result->pix);
	nde_commit_layer_position(c->item_id, c->result);
	g_clear_pointer(&c->result, nde_state_free);   /* pre-edit pixels — superseded */
	return TRUE;
}

/* Why a history edit is refused right now.  Both modes install a synthesized
 * state over the target's real pixels, so a second edit would compute against
 * something that is not the committed image. */
static const char *edit_in_progress_reason(void) {
	return nde_edit_at_active() ?
			_("an insertion point is open — finish or cancel it first") :
			_("another history step is being edited — close its dialog first");
}

/* The log change an amend or a delete makes: they differ only in whether there
 * are new params to write, so both edits carry one of these to nde_commit_replayed
 * rather than repeating the choice at every commit point. */
typedef struct {
	gint64       record_id;
	const gchar *new_params;   /* NULL = delete */
} edit_log_change;

static gboolean edit_log_commit(gpointer user, gchar **err) {
	const edit_log_change *l = user;
	return l->new_params ? nde_history_amend(l->record_id, l->new_params, err)
	                     : nde_history_delete(l->record_id, err);
}

gboolean nde_edit_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	edit_log_change log_change = { record_id, new_params };

	if (nde_amend_preview_installed()) {
		*err = g_strdup(edit_in_progress_reason());
		return FALSE;
	}

	/* Locate the record among the LIVE records and pre-check it, for
	 * user-facing errors before any heavy work.  (nde_history_amend/delete
	 * re-validate under the mutex at commit time.) */
	gint item_id = 0;
	gboolean found = FALSE, is_compositing = FALSE, is_joint = FALSE;
	gboolean is_joint_geometry = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->record_id == record_id) {
			found = TRUE;
			item_id = rec->target_item_id;
			is_compositing = nde_compositing_is_op(rec->op_id);
			is_joint = nde_joint_is_op(rec->op_id);
			is_joint_geometry = nde_joint_is_geometric_op(rec->op_id);
			if (new_params && is_compositing) {
				/* Validated here rather than by an op descriptor — these
				 * records have none (nde_compositing.h). */
				gchar *why = NULL;
				if (!nde_compositing_validate(rec->op_id, new_params, &why))
					*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): %s"),
					                       record_id, rec->op_id, why ? why : "?");
				g_free(why);
			} else if (new_params && !nde_record_amendable(rec)) {
				/* Amend needs editable params.  DELETE does not: removing
				 * an opaque record is well-defined — the trial chain never
				 * replays the deleted step, so deleting the one opaque
				 * record in a chain regains editability around it. */
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque and cannot be edited"),
				                       record_id, rec->op_id ? rec->op_id : "?");
			} else if (!new_params && !nde_record_deletable(rec)) {
				/* Structural (DOCUMENT scope) steps cannot resurrect what
				 * they destroyed or un-create a layer. */
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is a structural step and cannot be deleted"),
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

	/* Compositing-state edits take the cheap path: these records are inputs to
	 * the COMPOSITOR, not to any pixel chain, so there is nothing to replay —
	 * commit the log, re-fold the layer's properties, re-render.  No chain
	 * build (they are not members, so the freeze check below would reject them
	 * as "does not affect this image's pixels"), no pixel swap, no snapstore
	 * invalidation.  Undo is still flushed: the coupled props-only entries
	 * would redo the pre-edit value (nde_compositing.h). */
	if (is_compositing) {
		/* Collected while the edited record's position is still measurable:
		 * a tint edit changes what every joint record downstream of it
		 * derives, for ALL participants including this layer itself (whose
		 * own chain this cheap path never replays). */
		GArray *jt = nde_edit_joint_targets(item_id, record_id, TRUE);
		if (!edit_log_commit(&log_change, err)) {
			g_array_unref(jt);
			return FALSE;
		}
		if (item_id >= 0 && !flis_layer_get_by_id(item_id)
		    && nde_item_is_retained_input(item_id)) {
			/* The layer was consumed by a composite (flatten / merge), so
			 * there is nothing live to re-fold onto — but its recorded
			 * compositing state lives on inside the consuming composite's
			 * params, which is what the re-render reads.  Patch that from
			 * the amended log, then recompute every consumer: the edit
			 * flows through the composite boundary instead of failing
			 * with "the target layer no longer exists". */
			if (!nde_composite_refresh_input_state(item_id, err)) {
				g_array_unref(jt);
				return FALSE;
			}
			cascade_composite_consumers(item_id);
		} else if (!nde_compositing_recompute(item_id, err)) {
			g_array_unref(jt);
			return FALSE;
		}
		cascade_joint_targets(jt);
		g_array_unref(jt);
		undo_flush();
		gui_iface.set_progress(PROGRESS_DONE, _("Edit history updated"));
		return TRUE;
	}

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
	/* Deleting the ONE step a mask was built by would leave a mask item with no
	 * history to rebuild it from.  Caught here because the trial chain below is
	 * empty in that case, and an empty chain cannot be recognised as a mask's
	 * (that test reads the first member's op) — so it fell through to the image
	 * path and failed with "failed to load the baseline checkpoint", about an
	 * item the user never named. */
	if (!new_params && pos->is_mask && pos->records->len == 1) {
		nde_chain_free(pos);
		/* A mask a composite kept a copy of: deleting the step that built it
		 * means those composites composite WITHOUT that mask, which is what
		 * deleting it asks for.  Recorded rather than inferred — the consumers
		 * lose the mask from their params, so the log says plainly that none
		 * applies instead of leaving an input pointing at an item with no
		 * history and hoping the replay reads that as "no mask" rather than
		 * "lost". */
		if (nde_item_is_retained_input(item_id))
			return delete_retained_mask(record_id, item_id, err);
		*err = g_strdup(_("this is the step that built the mask — remove the mask "
		                  "itself from the Layers panel instead of deleting it here"));
		return FALSE;
	}
	nde_chain_free(pos);

	/* A JOINT record is a member of EVERY participant's chain (nde_joint.h):
	 * editing it recomputes them all, so the freeze rule must hold on each —
	 * a later opaque step on ANY participant locks the record, not just one
	 * on the anchor's chain. */
	if (is_joint) {
		GPtrArray *jsnap = nde_history_snapshot(NULL);
		guint jn = 0;
		gint *jparts = NULL;
		gchar *jold_params = NULL;
		for (guint i = 0; jsnap && i < jsnap->len; i++) {
			const nde_record *r = g_ptr_array_index(jsnap, i);
			if (r->record_id == record_id) {
				jparts = nde_joint_record_participants(r, &jn);
				jold_params = g_strdup(r->params);
				break;
			}
		}
		if (jsnap)
			g_ptr_array_unref(jsnap);
		/* A SUBSET amend (participant list changed): the commit re-derives
		 * the pins (nde_history_amend), removed participants revert through
		 * the pre-commit cascade collection, added ones through a second
		 * post-commit pass.  Two rules first: the anchor must stay a
		 * participant (its chain hosts the record), and a brand-new
		 * participant needs a replay starting point — while it has no
		 * records its current pixels ARE its original state, so a baseline
		 * is taken from them now. */
		if (new_params && jparts &&
		    !nde_joint_params_same_participants(jold_params, new_params)) {
			guint nn = 0;
			gint *nparts = nde_joint_params_participants(new_params, &nn);
			gboolean anchor_kept = FALSE;
			for (guint k = 0; nparts && k < nn; k++)
				anchor_kept = anchor_kept || nparts[k] == item_id;
			if (!anchor_kept) {
				*err = nparts ?
					g_strdup(_("the layer this step is recorded on must stay a participant — delete the step and re-run the operation instead")) :
					g_strdup(_("the new parameters name no readable participant list"));
				g_free(nparts);
				g_free(jold_params);
				g_free(jparts);
				return FALSE;
			}
			for (guint k = 0; k < nn && !*err; k++) {
				gboolean was = FALSE;
				for (guint j = 0; j < jn; j++)
					was = was || jparts[j] == nparts[k];
				if (was)
					continue;
				flis_layer_t *lay = flis_layer_get_by_id(nparts[k]);
				if (!lay || !lay->fit) {
					*err = g_strdup_printf(_("layer %d cannot join this step — it is not in the document"),
					                       nparts[k]);
				} else {
					/* Pixel-CHAIN membership, not any record: a layer whose
					 * only history is compositing state (a tint) still shows
					 * its original pixels, and those are its baseline. */
					nde_chain *pc = nde_chain_build(nparts[k]);
					gboolean no_pixel_history = pc->records->len == 0;
					nde_chain_free(pc);
					if (no_pixel_history)
						nde_checkpoint_baseline_ensure(lay->fit, nparts[k]);
				}
			}
			g_free(nparts);
			if (*err) {
				g_free(jold_params);
				g_free(jparts);
				return FALSE;
			}
		}
		g_free(jold_params);
		if (!jparts) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
			                       record_id, "flis joint");
			return FALSE;
		}
		for (guint k = 0; k < jn && !*err; k++) {
			if (jparts[k] == item_id)
				continue;   /* the anchor chain was checked above */
			nde_chain *pc = nde_chain_build(jparts[k]);
			gint pidx = -1;
			for (guint i = 0; i < pc->records->len; i++) {
				const nde_record *r = g_ptr_array_index(pc->records, i);
				if (r->record_id == record_id) {
					pidx = (gint)i;
					break;
				}
			}
			if (pidx >= 0 && (guint)pidx < pc->tail_start) {
				flis_layer_t *play = flis_layer_get_by_id(jparts[k]);
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is locked by a later opaque step on layer '%s'"),
				                       record_id,
				                       play && play->layer_name ? play->layer_name : "?");
			}
			nde_chain_free(pc);
		}
		g_free(jparts);
		if (*err)
			return FALSE;
	}

	/* For a delete, build the TRIAL chain: the deleted record is excluded
	 * from membership AND from the blockers, so its own opaqueness/mask
	 * state cannot veto its removal.  (Deleting the last barrier makes the
	 * trial's restart point fall back to the previous barrier/baseline.) */
	nde_chain *chain = nde_chain_build_excluding(item_id, new_params ? 0 : record_id);
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
	 * needs no edit here — nde_chain_build_excluding already omitted it.)
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
		/* Validate the new params against the op before replaying.  A composite
		 * node has none, and is checked against its own recorded blob instead:
		 * the compositing state may change, what it consumed may not. */
		if (nde_composite_is_op(target_rec->op_id)) {
			if (!nde_composite_validate(target_rec->params, new_params, err)) {
				nde_chain_free(chain);
				return FALSE;
			}
		} else {
			const op_descriptor *op = op_descriptor_by_id(target_rec->op_id);
			gpointer trial = (op && op->deserialize) ?
					op->deserialize(new_params, target_rec->op_version) : NULL;
			if (!trial) {
				*err = g_strdup_printf(_("the new parameters for '%s' failed to parse"),
				                       target_rec->op_id ? target_rec->op_id : "?");
				nde_chain_free(chain);
				return FALSE;
			}
			destroy_any_args(trial);
		}
		g_free(target_rec->params);
		target_rec->params = g_strdup(new_params);
	}

	/* A mask item's value is a mask, not an image: it is rebuilt from the
	 * image its first member read rather than from a baseline, and it lands
	 * in a mask slot rather than in a fits.  Everything before this point —
	 * policy, freeze position, the trial chain, the parameter substitution —
	 * is the same for both. */
	if (chain->is_mask) {
		/* Where in the mask's own chain the edit lands: pins BEFORE it see an
		 * unchanged mask and must not be disturbed. */
		guint mask_pos = 0;
		for (guint i = 0; i < chain->records->len; i++) {
			const nde_record *r = g_ptr_array_index(chain->records, i);
			if (r->record_id == record_id) {
				mask_pos = i;
				break;
			}
		}
		gui_iface.set_progress(0.f, _("Rebuilding the mask..."));
		fits *built = nde_mask_chain_replay(chain, chain->records->len, err);
		nde_chain_free(chain);
		if (!built) {
			gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
			return FALSE;
		}
		/* A mask whose layer was composited away has nowhere to be committed
		 * — the layer it belonged to is gone.  The replay above still had to
		 * run, to prove the edited chain applies; what shows the change is the
		 * composite that kept a copy of this mask, refreshed by the cascade
		 * below.
		 *
		 * A mask that was CLEARED after being used is the same shape: dormant,
		 * no live slot to commit into, but its consumers' pinned copies remain
		 * the states their replays read.  Editing its history refreshes those
		 * copies and recomputes the consumers; it does not resurrect the mask. */
		gboolean dormant = item_id != NDE_ITEM_PLAIN_MASK &&
		                   flis_item_lookup(item_id, NULL) == FLIS_ITEM_NONE;
		gboolean ok = nde_item_is_retained_input(item_id) || dormant ||
		              commit_mask_value(item_id, built, err);
		clearfits(built);
		free(built);
		if (!ok) {
			gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
			return FALSE;
		}
		if (!edit_log_commit(&log_change, err)) {
			gui_iface.set_progress(PROGRESS_RESET, _("The mask was rebuilt but the history could not be updated"));
			return FALSE;
		}
		/* Reverse invalidation: the mask changed, so anything that used it
		 * is now showing a result built from the old one.  Runs AFTER the
		 * log commit, because the consumers are recomputed from the log. */
		cascade_mask_consumers(item_id, mask_pos);
		nde_edit_finish(NULL, _("Edit history updated"));
		return TRUE;
	}

	gui_iface.set_progress(0.f, _("Recomputing edit history..."));
	/* Taken BEFORE the log is touched, because a delete removes the very
	 * record the position would be measured against (cascade_derived_masks). */
	const gint64 unchanged_upto = nde_log_predecessor(item_id, record_id);
	/* Convergence C3 invalidation, BEFORE replay: pool states at-or-after
	 * the edit describe the OLD chain and must not survive (they are never
	 * consulted as restarts for THIS edit, but they would be stale for the
	 * next one).  The replay below re-deposits fresh states as it goes.
	 * Joint factors likewise: the trial replays against a log that is not
	 * yet committed, so factors cached before it must not leak in — and the
	 * ones the trial caches must not survive a failed commit (below). */
	nde_snapstore_invalidate_from(item_id, record_id);
	nde_joint_cache_invalidate();
	/* Restart from the LATEST cached state at-or-before the edit (undo/redo
	 * entries, prior deposits, barrier checkpoints) rather than always the
	 * tail start — the frozen prefix's effect is embodied in the restart
	 * pixels either way. */
	if (!new_params)
		boundary = (guint)pos_idx;   /* the deleted member's former index */
	guint start_idx = 0;
	nde_state *start = nde_edit_restart_state(chain, boundary, record_id, &start_idx, err);
	nde_replay_arm_position(start, chain);
	/* A NULL restart state is not always a failure: an item born of a merge
	 * restarts from no state at all, its first member rendering its own
	 * inputs.  nde_edit_restart_state distinguishes the two by whether it set
	 * @err, and only returns NULL-without-error at start_idx 0. */
	nde_state *result = (!start && *err) ? NULL :
			nde_replay_apply_records(start, chain, start_idx,
			                         chain->records->len, err);
	nde_chain_free(chain);
	if (!result) {
		/* Deposits made by a failed replay describe an uncommitted chain —
		 * drop them along with anything else at-or-after the edit.  Factors
		 * the trial cached describe the uncommitted params — same fate. */
		nde_snapstore_invalidate_from(item_id, record_id);
		nde_joint_cache_invalidate();
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}
	/* The trial's factors were computed against the uncommitted log; the
	 * commit bumps the generation, so post-commit cascades recompute from
	 * the committed state.  Drop the trial's entries outright — a failed
	 * log commit below would otherwise leave them poisoning the still-
	 * current generation. */
	nde_joint_cache_invalidate();
	/* Collected while the edited record is still live: an amend or delete
	 * of anything on this chain changes what every joint record at-or-after
	 * it derives for its OTHER participants (a deleted joint record's own
	 * participants included — post-commit its list is gone). */
	GArray *joint_targets = nde_edit_joint_targets(item_id, record_id, FALSE);

	/* Resolve the target fits.  gfit for plain images; the layer's fit for
	 * FLIS (identical pointer when it is the active layer). */
	fits *target = gfit;
	gboolean retained = FALSE;
	if (item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(item_id);
		target = lay ? lay->fit : NULL;
		retained = nde_item_is_retained_input(item_id);   /* implies !lay */
	}
	nde_commit_ctx commit = {
		.item_id  = item_id,
		.target   = target,
		.result   = result,
		.retained = retained,
	};
	if (!nde_commit_replayed(&commit, edit_log_commit, &log_change, err)) {
		g_array_unref(joint_targets);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}

	/* The edit moved this participant's contribution, so every sibling's
	 * scaling is now derived from stale factors.  A subset amend may have
	 * ADDED participants the pre-commit collection could not see — merge a
	 * post-commit pass. */
	if (is_joint && new_params)
		joint_merge_post_commit_targets(joint_targets, item_id, record_id);
	/* Participants of a GEOMETRIC joint record are sitting where its warp put
	 * them, so their recompute has to re-anchor from the baseline — the chain
	 * it replays may no longer contain anything that moves them (deleting the
	 * registration is exactly that case). */
	nde_cascade_from(item_id, unchanged_upto, joint_targets, is_joint_geometry);
	g_array_unref(joint_targets);

	/* NULL for a retained input — nothing live to drop statistics for. */
	nde_edit_finish(target, _("Edit history updated"));
	return TRUE;
}

gboolean nde_amend_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(new_params != NULL, FALSE);
	return nde_edit_execute(record_id, new_params, err);
}

/* ======================================================================= */
/* Undoing a composite: bringing the consumed layers back                  */
/* ======================================================================= */

/* Everything the log holds against @item at or after @from_pos, plus anything
 * held against that item's masks — all of it describes an image that is about
 * to stop existing.  Positions, because ids do not order the log. */
static GArray *records_from(GPtrArray *snap, guint from_pos, gint item) {
	GArray *ids = g_array_new(FALSE, FALSE, sizeof(gint64));
	gint lmask = 0, pmask = 0;
	flis_layer_t *lay = flis_layer_get_by_id(item);
	if (lay) {
		lmask = lay->lmask_item_id;
		pmask = lay->pmask_item_id;
	}
	for (guint i = from_pos; i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (r->target_item_id == item ||
		    (lmask && r->target_item_id == lmask) ||
		    (pmask && r->target_item_id == pmask))
			g_array_append_val(ids, r->record_id);
	}
	return ids;
}

/* Locate the composite and check that undoing it is possible at all.  Every
 * refusal here names what is in the way: the whole point of offering this is
 * that the alternative was silence. */
static const nde_record *composite_undo_check(GPtrArray *snap, gint64 record_id,
                                              guint *pos, nde_composite_state **out_st,
                                              gchar **err) {
	const nde_record *rec = NULL;
	for (guint i = 0; i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (r->record_id == record_id) {
			rec = r;
			*pos = i;
			break;
		}
	}
	if (!rec) {
		*err = g_strdup(_("that step is no longer in the history"));
		return NULL;
	}
	if (!nde_composite_is_op(rec->op_id)) {
		*err = g_strdup(_("only a merge or a flatten can be undone this way"));
		return NULL;
	}
	nde_composite_state *st = nde_composite_state_parse(rec->params);
	if (!st) {
		*err = g_strdup(_("this step predates the format that records what it "
		                  "consumed, so the layers it merged cannot be "
		                  "reconstructed"));
		return NULL;
	}
	const gint produced = rec->target_item_id;
	if (!flis_layer_get_by_id(produced)) {
		nde_composite_state_free(st);
		*err = nde_item_is_retained_input(produced) ?
			g_strdup(_("a later merge or flatten consumed the image this step "
			           "produced — undo that one first")) :
			g_strdup(_("the image this step produced is no longer in the document"));
		return NULL;
	}
	/* Every consumed layer has to be rebuildable from its own history, or the
	 * undo would bring a layer back holding pixels that are not what was
	 * merged — worse than not offering it. */
	for (guint i = 0; i < st->inputs->len; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		nde_chain *c = nde_chain_build(in->item_id);
		gboolean ok = c->replayable ||
		              (c->records->len == 0 && nde_checkpoint_baseline_exists(in->item_id));
		if (!ok) {
			*err = g_strdup_printf(_("'%s' cannot be rebuilt, so undoing this step "
			                         "would bring it back with the wrong pixels: %s"),
			                       in->name ? in->name : _("a consumed layer"),
			                       c->reasons->len ?
			                           (char *)g_ptr_array_index(c->reasons, 0) :
			                           _("its original pixels were not kept"));
			nde_chain_free(c);
			nde_composite_state_free(st);
			return NULL;
		}
		nde_chain_free(c);
	}
	*out_st = st;
	return rec;
}

gchar *nde_composite_undo_describe(gint64 record_id, guint *n_layers,
                                   guint *n_discarded, gchar **err) {
	g_return_val_if_fail(err != NULL, NULL);
	*err = NULL;
	if (n_layers) *n_layers = 0;
	if (n_discarded) *n_discarded = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	if (!snap) {
		*err = g_strdup(_("no edit history"));
		return NULL;
	}
	guint pos = 0;
	nde_composite_state *st = NULL;
	const nde_record *rec = composite_undo_check(snap, record_id, &pos, &st, err);
	if (!rec) {
		g_ptr_array_unref(snap);
		return NULL;
	}
	GArray *doomed = records_from(snap, pos, rec->target_item_id);
	GString *s = g_string_new(NULL);
	for (guint i = 0; i < doomed->len; i++) {
		const gint64 id = g_array_index(doomed, gint64, i);
		for (guint j = 0; j < snap->len; j++) {
			const nde_record *r = g_ptr_array_index(snap, j);
			if (r->record_id != id)
				continue;
			g_string_append_printf(s, "  • %s\n",
			                       r->summary ? r->summary : r->op_id);
			break;
		}
	}
	if (n_layers) *n_layers = st->inputs->len;
	if (n_discarded) *n_discarded = doomed->len;
	g_array_unref(doomed);
	nde_composite_state_free(st);
	g_ptr_array_unref(snap);
	return g_string_free(s, FALSE);
}

gboolean nde_composite_undo_execute(gint64 record_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	if (!snap) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}
	guint pos = 0;
	nde_composite_state *st = NULL;
	const nde_record *rec = composite_undo_check(snap, record_id, &pos, &st, err);
	if (!rec) {
		g_ptr_array_unref(snap);
		return FALSE;
	}
	const gint produced = rec->target_item_id;
	GArray *doomed = records_from(snap, pos, produced);

	/* Rebuild every consumed layer's pixels — and fetch every recorded layer
	 * mask's stored copy — BEFORE touching the document, so a failure half
	 * way leaves the image exactly as it was.  The mask is a pinned copy, not
	 * a replay (nde_composite.h); a recorded copy that is no longer stored
	 * refuses the undo the same way an unbuildable layer does, because
	 * restoring the layer without its mask would not restore the document. */
	guint n = st->inputs->len;
	fits **pix = g_new0(fits *, n);
	fits **msk = g_new0(fits *, n);
	gboolean ok = TRUE;
	for (guint i = 0; i < n && ok; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		nde_chain *c = nde_chain_build(in->item_id);
		/* A layer that was never edited has no chain, only the baseline the
		 * merge took of it — which is exactly its pixels at that moment. */
		pix[i] = c->records->len ? nde_chain_replay(c, err)
		                         : nde_state_release(nde_checkpoint_baseline_get(in->item_id));
		if (!pix[i]) {
			if (!*err)
				*err = g_strdup_printf(_("'%s' could not be rebuilt"),
				                       in->name ? in->name : _("a consumed layer"));
			ok = FALSE;
		}
		nde_chain_free(c);
		if (ok && in->mask_item_id) {
			const nde_pin *mp = nde_record_input_by_item(rec, in->mask_item_id);
			msk[i] = mp ? nde_state_release(nde_checkpoint_get_at(mp->item_id,
			                                                      mp->record_id))
			            : NULL;
			if (!msk[i]) {
				*err = g_strdup_printf(_("the stored layer mask of '%s' is no "
				                         "longer available"),
				                       in->name ? in->name : _("a consumed layer"));
				ok = FALSE;
			}
		}
	}
	if (!ok) {
		for (guint i = 0; i < n; i++) {
			if (pix[i]) { clearfits(pix[i]); free(pix[i]); }
			if (msk[i]) { clearfits(msk[i]); free(msk[i]); }
		}
		g_free(pix);
		g_free(msk);
		g_array_unref(doomed);
		nde_composite_state_free(st);
		g_ptr_array_unref(snap);
		return FALSE;
	}

	/* The document: out with the flattened layer, back with the ones it ate. */
	flis_layer_t *flat = flis_layer_get_by_id(produced);
	if (com.uniq) {
		com.uniq->canvas_w = st->canvas_w ? st->canvas_w : com.uniq->canvas_w;
		com.uniq->canvas_h = st->canvas_h ? st->canvas_h : com.uniq->canvas_h;
		com.uniq->canvas_bg_r = st->bg_r;
		com.uniq->canvas_bg_g = st->bg_g;
		com.uniq->canvas_bg_b = st->bg_b;
	}
	GHashTable *group_map = g_hash_table_new(g_direct_hash, g_direct_equal);
	for (guint i = 0; i < st->groups->len; i++) {
		const nde_composite_group *g = &g_array_index(st->groups, nde_composite_group, i);
		/* A flatten/merge only frees the layers it consumed, never the group
		 * objects — those still live in com.uniq->groups.  Blindly re-adding
		 * every recorded group would therefore append a duplicate flis_group_t
		 * with the same item_id (writing duplicate rows on save, a doubled
		 * entry in "Move to group", and a leaked per-render sub-composite).
		 * So reuse the survivor when it exists and only reissue its recorded
		 * properties; add a fresh group only when none is present. */
		flis_group_t *grp = flis_group_get_by_id(g->item_id);
		if (grp) {
			if (g->name) {
				g_free(grp->name);
				grp->name = g_strdup(g->name);
			}
		} else {
			grp = flis_group_add(g->name ? g->name : _("Group"));
			if (!grp)
				continue;
			grp->item_id = g->item_id;   /* the id its members were recorded with */
		}
		grp->blend_mode = (flis_blend_mode_t)g->blend_mode;
		grp->opacity    = (gfloat)g->opacity;
		grp->visible    = g->visible;
		g_hash_table_insert(group_map, GINT_TO_POINTER(g->item_id), grp);
	}
	for (guint i = 0; i < n; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		flis_layer_t *lay = flis_layer_add(pix[i], in->name ? in->name : _("Layer"));
		if (!lay) {
			clearfits(pix[i]);
			free(pix[i]);
			if (msk[i]) { clearfits(msk[i]); free(msk[i]); msk[i] = NULL; }
			continue;
		}
		/* The identity is restored, not reissued: the layer's own history is
		 * keyed on it and would otherwise be orphaned by its own undo. */
		lay->item_id     = in->item_id;
		lay->blend_mode  = (flis_blend_mode_t)in->blend_mode;
		lay->opacity     = (gfloat)in->opacity;
		lay->position_x  = in->position_x;
		lay->position_y  = in->position_y;
		lay->visible     = in->visible;
		/* The tint was recorded with the rest of the input state but never
		 * restored — an undone flatten/merge came back untinted. */
		lay->has_tint    = in->has_tint;
		lay->layer_tint.r = in->tint_r;
		lay->layer_tint.g = in->tint_g;
		lay->layer_tint.b = in->tint_b;
		/* Back into its recorded stacking slot, not on top of layers the
		 * composite never consumed — a merge-down undo in a 3+ layer document
		 * would otherwise reorder the stack and change the render.  0 = an
		 * older record that carried no order: keep flis_layer_add's top slot,
		 * as before.  The list is re-sorted after the loop. */
		if (in->layer_order)
			lay->layer_order = in->layer_order;
		if (in->group_id && g_hash_table_contains(group_map, GINT_TO_POINTER(in->group_id)))
			flis_layer_set_group(lay, in->group_id);
		/* Reinstall the layer mask from its stored copy — same unwrapping as
		 * nde_composite_render — and re-link the mask ITEM, so records that
		 * target it name a mask a layer carries again. */
		if (msk[i]) {
			mask_t *m = fits_to_mask(msk[i]);
			layermask_t *lm = m ? calloc(1, sizeof(layermask_t)) : NULL;
			if (lm) {
				lm->w      = msk[i]->rx;
				lm->h      = msk[i]->ry;
				lm->bitpix = m->bitpix;
				lm->data   = m->data;   /* moved */
				free(m);
				if (flis_layer_set_lmask(lay, lm)) {
					/* Cannot happen for a copy stored with the layer's own
					 * dimensions; degrade rather than abort mid-restore. */
					layermask_free(lm);
					siril_log_warning(_("the layer mask of '%s' could not be "
					                    "restored\n"),
					                  in->name ? in->name : "?");
				} else {
					lay->lmask_item_id = in->mask_item_id;
					lay->lmask_active  = TRUE;
				}
			} else if (m) {
				free_mask(m);
			}
			clearfits(msk[i]);
			free(msk[i]);
		}
		/* The layer took the fits POINTER (flis_layer_new: layer->fit = fit),
		 * so there is nothing here left to free — freeing it left the layer
		 * pointing at released memory, which only showed up as a crash in the
		 * test teardown under MALLOC_PERTURB_. */
	}
	g_hash_table_destroy(group_map);
	g_free(pix);
	g_free(msk);
	if (flat)
		flis_layer_remove(flat);
	flis_sort_layer_stack();

	/* The log last, so a failure above leaves it describing what is there. */
	nde_history_drop_records(doomed);
	guint discarded = doomed->len;
	g_array_unref(doomed);

	/* Activate a layer the undo actually brought back — the first input (the
	 * merge's survivor) — not whatever happens to sit at the bottom of the
	 * stack. */
	flis_layer_t *back = flis_layer_get_by_id(
	    g_array_index(st->inputs, nde_composite_input, 0).item_id);
	gint back_idx = back ? g_slist_index(com.uniq->layers, back) : -1;
	nde_composite_state_free(st);
	g_ptr_array_unref(snap);

	uniq_set_active_layer(com.uniq, back_idx >= 0 ? back_idx : 0);
	gfit = flis_active_layer_fit();
	undo_flush();
	gui_iface.flis_invalidate_composite();
	gui_iface.invalidate_histogram();
	siril_log_message(_("The merge was undone: %u layer(s) came back and "
	                    "%u step(s) after it were discarded\n"),
	                  n, discarded);
	return TRUE;
}

gboolean nde_delete_execute(gint64 record_id, gchar **err) {
	return nde_edit_execute(record_id, NULL, err);
}

/* Where a reorder puts the moved record: before @before_id, or last when 0. */
typedef struct {
	gint64 record_id;
	gint64 before_id;
} reorder_log_change;

static gboolean reorder_log_commit(gpointer user, gchar **err) {
	const reorder_log_change *l = user;
	return nde_history_reorder(l->record_id, l->before_id, err);
}

gboolean nde_reorder_execute(gint64 record_id, gint64 anchor_id, gboolean after, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (nde_amend_preview_installed()) {
		*err = g_strdup(edit_in_progress_reason());
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
			else if (nde_joint_is_op(rec->op_id))
				/* One log position shared by every participant's chain: a
				 * sound multi-chain move semantics is not defined (v1). */
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) spans several layers and cannot be reordered"),
				                       record_id, rec->op_id);
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
	/* An item born of a merge has no state before that merge: member 0 IS its
	 * origin, rendering the composite's inputs where every other chain starts
	 * from a baseline (nde_replay.h).  Move something in front of it and the
	 * replay has an ordinary op to run against no state at all, which it duly
	 * dereferenced.  Only this direction needs stating: the composite itself
	 * can never be the record being moved, the reorderable-record check above
	 * having refused it already. */
	if (chain->from_composite && dest == 0) {
		*err = g_strdup(_("this image was produced by that merge, so nothing "
		                  "can come before it — move the step within one of "
		                  "the layers the merge consumed instead"));
		nde_chain_free(chain);
		return FALSE;
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

	/* Taken while the log still holds the pre-move order (cascade_derived_masks). */
	const gint64 unchanged_upto = nde_log_predecessor(item_id, boundary_pre_id);

	/* Permute the trial chain in place (records + parallel flags). */
	nde_record *moved = g_ptr_array_steal_index(chain->records, (guint)i_old);
	g_ptr_array_insert(chain->records, dest, moved);
	guint8 fl = g_array_index(chain->member_flags, guint8, i_old);
	g_array_remove_index(chain->member_flags, (guint)i_old);
	g_array_insert_val(chain->member_flags, (guint)dest, fl);

	gui_iface.set_progress(0.f, _("Recomputing edit history..."));
	nde_snapstore_invalidate_from(item_id, inval_min);
	nde_joint_cache_invalidate();
	guint start_idx = 0;
	nde_state *start = nde_edit_restart_state(chain, min_idx, boundary_pre_id,
	                                          &start_idx, err);
	nde_replay_arm_position(start, chain);
	/* A NULL restart state is not always a failure: an item born of a merge
	 * restarts from no state at all, its first member rendering its own
	 * inputs.  nde_edit_restart_state distinguishes the two by whether it set
	 * @err, and only returns NULL-without-error at start_idx 0. */
	nde_state *result = (!start && *err) ? NULL :
			nde_replay_apply_records(start, chain, start_idx,
			                         chain->records->len, err);
	nde_chain_free(chain);
	if (!result) {
		nde_snapstore_invalidate_from(item_id, inval_min);
		nde_joint_cache_invalidate();
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}
	/* Same trial-cache discipline as nde_edit_execute: the factors above were
	 * derived from the permuted-but-uncommitted order. */
	nde_joint_cache_invalidate();
	/* A move across a joint member changes the prefix it reads — collected
	 * pre-commit, cascaded post-commit like every other joint disturbance. */
	GArray *joint_targets = nde_edit_joint_targets(item_id, boundary_pre_id, FALSE);

	/* Atomic commit — the same tail as nde_edit_execute, including its retained
	 * branch: an item a merge consumed has no layer to swap into, and the
	 * replay above was run to prove the reordered chain still applies rather
	 * than to produce pixels anyone keeps.  Without this, reordering a step
	 * on a consumed layer failed with "the record's target layer no longer
	 * exists" — true, and beside the point. */
	fits *target = gfit;
	gboolean retained = FALSE;
	if (item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(item_id);
		target = lay ? lay->fit : NULL;
		retained = nde_item_is_retained_input(item_id);   /* implies !lay */
	}
	reorder_log_change log_change = { record_id, log_before_id };
	nde_commit_ctx commit = {
		.item_id  = item_id,
		.target   = target,
		.result   = result,
		.retained = retained,
	};
	if (!nde_commit_replayed(&commit, reorder_log_commit, &log_change, err)) {
		/* The replay deposited states for the permuted order, which is not
		 * what the log now says. */
		nde_snapstore_invalidate_from(item_id, inval_min);
		g_array_unref(joint_targets);
		gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
		return FALSE;
	}

	/* A mask built from these pixels read a prefix the move may have
	 * reordered; a joint record reads the participant's changed prefix. */
	nde_cascade_from(item_id, unchanged_upto, joint_targets, FALSE);
	g_array_unref(joint_targets);

	nde_edit_finish(target, _("Edit history updated"));
	return TRUE;
}
