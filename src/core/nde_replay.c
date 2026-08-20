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
#include "yyjson.h"

/* ---- Tier-C: replayable Python scripts (nde-phase5) -------------------- */

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
	gchar *invalid = nde_tier_c_invalid_reason(rec);
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

	gboolean quiesced = nde_commit_lock(gfit);
	fits_swap_all_except_rwlock(gfit, scratch);
	nde_commit_unlock(gfit, quiesced);

	int rc = execute_python_script(g_strdup(script),
	                               TRUE,   /* from_file */
	                               TRUE,   /* sync */
	                               argv,
	                               FALSE,  /* is_temp_file */
	                               TRUE,   /* from_cli: headless script path */
	                               FALSE,  /* debug_mode */
	                               TRUE);  /* for_replay */

	quiesced = nde_commit_lock(gfit);
	fits_swap_all_except_rwlock(gfit, scratch);
	nde_commit_unlock(gfit, quiesced);

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

/* ---- pinned mask inputs (see nde_replay.h) ----------------------------- */

/* A pin's state lives at the same coordinate the checkpoint store already
 * uses: record 0 means "that item's baseline", anything else is that record's
 * output.  Masks therefore need no store of their own — a mask is a mono
 * image, and its item and record ids come from the same sequences. */
void nde_mask_pin_store(const nde_record *rec, const fits *fit) {
	const nde_input_pin *pin = nde_record_input(rec, "mask");
	if (!pin || !fit || !fit->mask || !fit->mask->data)
		return;
	fits *mfit = mask_to_fits((fits *)fit);
	if (!mfit)
		return;
	nde_checkpoint_store_at(mfit, pin->src_item_id, pin->src_record_id);
	clearfits(mfit);
	free(mfit);
}

gboolean nde_mask_pin_resolvable(const nde_record *rec) {
	const nde_input_pin *pin = nde_record_input(rec, "mask");
	if (!pin)
		return FALSE;
	return nde_checkpoint_exists_at(pin->src_item_id, pin->src_record_id);
}

/* Install @rec's pinned mask on @scratch for the duration of one op.  Returns
 * FALSE + @err when the pin does not resolve; the caller must not proceed,
 * since running the op unmasked would produce different pixels while claiming
 * to have reproduced the record. */
gboolean nde_mask_pin_install(fits *scratch, const nde_record *rec, gchar **err) {
	const nde_input_pin *pin = nde_record_input(rec, "mask");
	if (!pin)
		return TRUE;
	fits *mfit = nde_checkpoint_get_at(pin->src_item_id, pin->src_record_id);
	if (!mfit) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): its mask is no longer stored"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}
	gboolean ok = FALSE;
	if (mfit->rx == scratch->rx && mfit->ry == scratch->ry) {
		mask_t *m = fits_to_mask(mfit);
		if (m) {
			if (scratch->mask)
				free_mask(scratch->mask);
			scratch->mask = m;
			scratch->mask_active = TRUE;
			ok = TRUE;
		}
	}
	if (!ok)
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): its mask does not fit the image"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
	clearfits(mfit);
	free(mfit);
	return ok;
}

/* The mask belongs to the record, not to the chain: leave nothing behind for
 * the next one, which has its own pin or none. */
void nde_mask_pin_clear(fits *scratch) {
	if (scratch->mask) {
		free_mask(scratch->mask);
		scratch->mask = NULL;
	}
	scratch->mask_active = FALSE;
}

/* Swap @result's pixels into @target — the commit every replay ends with —
 * LEAVING THE MASK SLOT WHERE IT IS.
 *
 * A mask is an item in its own right, with a history of its own; a replay of
 * the pixels neither produces one (nde_mask_pin_clear runs after every masked
 * member, so a replayed result carries no mask at all) nor has any business
 * replacing one.  The wholesale swap did both: every amend, delete, reorder or
 * insertion on an image silently emptied that image's processing mask slot,
 * because the replay's NULL went in and the user's mask came out and was freed
 * with the discarded pixels.  Pre-swapping the slots means the wholesale swap
 * below puts them back, and the value @result then holds is the replay's own
 * (usually none), freed with it.
 *
 * Caller holds @target's writer lock, as it did for the bare swap. */
/* Tail refresh for every history edit that changed pixels: notify rebuilds
 * the remap/tile buffers, but queues NO paint — without the explicit redraw
 * the screen keeps the pre-edit frame until the next mouseover (the same
 * trap end_generic_layer documents).  Unconditional notify: on a FLIS the
 * displayed COMPOSITE reflects non-active layers too, so gating it on
 * target == gfit left non-active-layer edits invisible. */
static void edit_refresh_display(void) {
	notify_gfit_data_modified();
	gui_iface.redraw_image(REDRAW_ALL);
}

void nde_commit_pixels(fits *target, fits *result) {
	mask_t *m = target->mask;
	target->mask = result->mask;
	result->mask = m;
	gboolean active = target->mask_active;
	target->mask_active = result->mask_active;
	result->mask_active = active;
	fits_swap_all_except_rwlock(target, result);
}

/* Writer-lock @fit for a short commit window, quiescing the lazy-tile
 * materialise pool first when @fit IS the displayed image: GRWLock has no
 * writer preference, so an unquiesced writer starves behind the pool's
 * back-to-back reader-locked tile fills on a large image.  Counting
 * suppression (gui_iface_impl.c), so nesting under an outer suppression
 * is safe.  Returns whether the pool was quiesced; pass it back to
 * nde_commit_unlock(). */
gboolean nde_commit_lock(fits *fit) {
	gboolean quiesce = (fit == gfit);
	if (quiesce)
		gui_iface.set_suppress_redraws(TRUE);
	g_rw_lock_writer_lock(&fit->rwlock);
	return quiesce;
}

void nde_commit_unlock(fits *fit, gboolean quiesced) {
	g_rw_lock_writer_unlock(&fit->rwlock);
	if (quiesced)
		gui_iface.set_suppress_redraws(FALSE);
}

/* Where a replay's layer starts.  @restart_id 0 means the baseline; anything
 * else is that record's output checkpoint.  FALSE when the chain carries no
 * geometry member (a plain image, or a layer nothing ever moved), in which
 * case the hooks are handed NULL and move nothing. */

/* Set while an edit to a GEOMETRIC joint record (flis.register) cascades to
 * its participants.  Those layers are sitting where the record's warp put
 * them, so the recompute has to re-anchor them from the baseline even when
 * the chain it replays now has no geometry member left to move them — which
 * is exactly what deleting the registration produces.  Deliberately NOT the
 * general rule: a layer the user dragged carries a position no record
 * describes, and forcing the baseline on every recompute would undo it. */
static gboolean joint_geometry_reanchor = FALSE;

void nde_replay_set_joint_reanchor(gboolean on) {
	joint_geometry_reanchor = on;
}

gboolean nde_replay_start_offset(const nde_chain *chain, gint64 restart_id,
                                    gint *pos_x, gint *pos_y) {
	if (!chain->has_geometry && !joint_geometry_reanchor)
		return FALSE;
	if (joint_geometry_reanchor && !chain->has_geometry)
		return nde_checkpoint_baseline_get_offset(chain->item_id, pos_x, pos_y);
	if (restart_id > 0 &&
	    nde_checkpoint_output_get_offset(restart_id, pos_x, pos_y))
		return TRUE;
	return nde_checkpoint_baseline_get_offset(chain->item_id, pos_x, pos_y);
}

/* Move the layer to where the replay left it.  The pixels are committed by
 * the caller's swap; without this the layer would keep the position its
 * pre-edit geometry produced, which the new pixels no longer match. */
void nde_commit_layer_offset(gint item_id, gint pos_x, gint pos_y) {
	if (item_id < 0)
		return;
	flis_layer_t *lay = flis_layer_get_by_id(item_id);
	if (!lay)
		return;
	lay->position_x = pos_x;
	lay->position_y = pos_y;
	gui_iface.flis_invalidate_composite();
}

/* What every successful history edit does once its own work is committed.
 * @target is the fits whose pixels changed, or NULL when the edit changed no
 * live pixels of its own (a mask value, or an item a composite consumed) —
 * there are then no cached statistics to drop.  @done_msg is the caller's,
 * already translated, because what was done differs and the user is told. */
void nde_edit_finish(fits *target, const char *done_msg) {
	/* No meta-undo (sketch §7): stale undo entries would restore pixels the
	 * log no longer describes. */
	undo_flush();
	if (target)
		invalidate_stats_from_fit(target);
	gui_iface.invalidate_histogram();
	if (is_current_image_flis())
		gui_iface.flis_invalidate_composite();
	edit_refresh_display();
	gui_iface.set_progress(PROGRESS_DONE, done_msg);
}

static fits *resolve_item_state(gint item_id, gint64 upto_record_id, gchar **err);
static fits *resolve_item_state_pos(gint item_id, gint64 upto_record_id,
                                    gint *pos_x, gint *pos_y,
                                    gboolean *pos_valid, gchar **err);
static fits *resolve_item_state_pos_bound(gint item_id, gint64 upto_record_id,
                                          gboolean exclusive,
                                          gint *pos_x, gint *pos_y,
                                          gboolean *pos_valid, gchar **err);

/* Run one composite member (graph step 7, nde_composite.h) — merge-down with
 * two inputs, flatten with all of them.
 *
 * @base is the accumulated replay state, which is one of the inputs already
 * resolved for free: the one whose item this chain belongs to.  That input's
 * pin is deliberately NOT followed — resolving it through the chain machinery
 * would recurse, since its chain contains this very record.  Every other input
 * is a live edge, re-derived by replaying its own chain.
 *
 * Returns a new fits; @base is left alone on both paths and stays the caller's
 * to free.  On success the layer's position becomes the canvas origin, because
 * the composite is canvas-sized and canvas-aligned — the same reset the live
 * merge and flatten perform.
 */
static fits *composite_apply(fits *base, const nde_record *rec, gint item_id,
                             gint *pos_x, gint *pos_y, gchar **err) {
	nde_composite_state *st = nde_composite_state_parse(rec->params);
	if (!st) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT ": the composite's inputs were not recorded"),
		                       rec->record_id);
		return NULL;
	}
	guint n = st->inputs->len;
	fits **pixels = g_new0(fits *, n);
	fits **masks = g_new0(fits *, n);
	gboolean *owned = g_new0(gboolean, n);
	gboolean ok = TRUE;
	for (guint i = 0; i < n && ok; i++) {
		nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		/* The mask is a stored copy, not a replay: it is read for every visible
		 * masked input, including the one whose pixels are @base. */
		if (in->visible && in->was_masked) {
			const nde_input_pin *mp = in->mask_item_id ?
					nde_record_input_by_item(rec, in->mask_item_id) : NULL;
			masks[i] = mp ? nde_checkpoint_get_at(mp->src_item_id, mp->src_record_id) : NULL;
			if (!masks[i]) {
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT ": the layer mask of '%s' is no longer stored"),
				                       rec->record_id, in->name ? in->name : "?");
				ok = FALSE;
				break;
			}
		}
		if (in->item_id == item_id) {
			pixels[i] = base;
			/* Reached only if a composite record ever sits mid-chain on its
			 * own item.  The current capture reissues the survivor's identity
			 * before committing (image_format_flis.c), so today no recorded
			 * input matches the target and base is always NULL — resolved
			 * inputs are re-anchored below instead. */
			if (pos_x && pos_y) {
				in->position_x = *pos_x;
				in->position_y = *pos_y;
			}
			continue;
		}
		if (!in->visible)
			continue;   /* contributes nothing: not worth a replay */
		const nde_input_pin *pin = nde_record_input_by_item(rec, in->item_id);
		if (!pin) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT ": the layer '%s' it consumed is not recorded as an input"),
			                       rec->record_id, in->name ? in->name : "?");
			ok = FALSE;
			break;
		}
		gint rpx = 0, rpy = 0;
		gboolean rpos = FALSE;
		pixels[i] = resolve_item_state_pos(pin->src_item_id, pin->src_record_id,
		                                   &rpx, &rpy, &rpos, err);
		owned[i] = TRUE;
		ok = pixels[i] != NULL;
		if (ok && rpos) {
			/* A geometry member owns this input's position — the same rule
			 * nde_commit_layer_offset applies to a live layer's replay — so an
			 * amended geometry step (a crop moved to a new origin) moves
			 * where the input lands in the composite.  The recorded position
			 * stays authoritative only for chains that never moved. */
			in->position_x = rpx;
			in->position_y = rpy;
		}
	}
	fits *out = ok ? nde_composite_render(st, pixels, masks, err) : NULL;
	for (guint i = 0; i < n; i++) {
		if (owned[i] && pixels[i]) {
			clearfits(pixels[i]);
			free(pixels[i]);
		}
		if (masks[i]) {
			clearfits(masks[i]);
			free(masks[i]);
		}
	}
	g_free(pixels);
	g_free(masks);
	g_free(owned);
	nde_composite_state_free(st);
	if (out && pos_x && pos_y)
		*pos_x = *pos_y = 0;
	return out;
}

/* Apply members [from..upto) to @scratch (consumed on failure).  Returns
 * @scratch on success, NULL + @err on failure.
 *
 * @pos_x / @pos_y carry the LAYER VALUE's other half (graph step 5): a
 * geometry member moves the layer as well as changing its pixels, and each
 * move is relative to the previous position, so the driver threads the
 * position through the run and hands it to the hooks.  NULL for a plain
 * image, which has no position, and for chains with no geometry member. */
/* The record currently being applied by the replay driver, for hooks whose
 * recompute needs record-scoped side data — the photometric pipeline reads
 * its embedded star catalogue (nde_cat.h) under this id.  0 outside a
 * replayed record.  Single conductor thread; a plain static suffices. */
static gint64 replay_current_record;

gint64 nde_replay_current_record_id(void) {
	return replay_current_record;
}

void nde_replay_set_current_record(gint64 record_id) {
	replay_current_record = record_id;
}

fits *nde_replay_apply_records(fits *scratch, const nde_chain *chain,
                                  guint from, guint upto,
                                  gint *pos_x, gint *pos_y, gchar **err) {
	for (guint i = from; i < upto; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		replay_current_record = rec->record_id;
		if (!processing_should_continue()) {
			*err = g_strdup(_("cancelled"));
			goto fail;
		}
		const nde_op_class *cls = nde_op_class_for(rec->op_id);
		/* A NULL scratch is legal for exactly one kind of record: the
		 * composite that gives an item born of a merge its origin, which
		 * renders its own inputs and wants no prior state.  Anything else
		 * needs pixels to work on, and used to take the absence of them as
		 * far as scratch->mask.  Callers guard their own restart points; this
		 * is the last one, so that a gap upstream is a failed replay and not
		 * a lost session. */
		if (!scratch && cls->family != NDE_OPC_COMPOSITE) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) has no "
			                         "image to run on — the history has no "
			                         "starting state before it"),
			                       rec->record_id, rec->op_id ? rec->op_id : "?");
			goto fail;
		}
		if (rec->tier == NDE_TIER_C) {
			/* Replayable Python script: re-run it on the accumulated state
			 * (nde-phase5).  Chain membership already vetted the recipe via
			 * member_invalid_reason, but the gate re-checks — the file can
			 * change between chain build and execution. */
			if (tier_c_rerun(scratch, rec, err))
				goto fail;
			/* The script's own commands went through start_in_new_thread,
			 * which sets the busy cursor on submission and clears it when the
			 * job ends — and the flag is a boolean, not a count, so the last
			 * of them cleared the conductor's.  Take it back: there may be
			 * many records still to replay after this one. */
			gui_iface.set_busy(TRUE);
			nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
			continue;
		}
		/* Every family is listed so that adding one is a compiler warning
		 * here rather than a record that silently takes the descriptor
		 * path below and does the wrong thing. */
		switch (cls->family) {
		case NDE_OPC_JOINT: {
			/* A joint record (nde_joint.h): recompute the multi-layer
			 * analysis — siblings resolved positionally, this chain's state
			 * supplied by the accumulated scratch — and apply only this
			 * participant's share.  The result is generation-cached, so the
			 * N participant replays of one edit share one analysis.
			 *
			 * Dispatch on the op: the SCALAR joint ops derive a per-layer
			 * affine, but registration derives a WARP and a new frame size,
			 * which the affine path cannot express — handing it to
			 * nde_joint_factor_for_item would multiply the layer by a
			 * meaningless number and leave it the wrong size. */
			if (nde_joint_is_geometric_op(rec->op_id)) {
				if (!nde_joint_register_apply(rec, scratch, chain->item_id,
				                              pos_x, pos_y, err))
					goto fail;
				nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
				continue;
			}
			double s = 1.0, o = 0.0;
			if (!nde_joint_factor_for_item(rec, scratch, chain->item_id,
			                               &s, &o, err))
				goto fail;
			flis_affine_layer_pixels(scratch, s, o);
			invalidate_stats_from_fit(scratch);
			nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
			continue;
		}
		case NDE_OPC_COMPOSITE: {
			fits *merged = composite_apply(scratch, rec, chain->item_id,
			                               pos_x, pos_y, err);
			if (!merged)
				goto fail;   /* scratch untouched — the fail path frees it */
			/* NULL when this composite is the item's origin: an item born of
			 * a merge has no prior state for the composite to replace. */
			if (scratch) {
				clearfits(scratch);
				free(scratch);
			}
			scratch = merged;
			nde_snapstore_deposit(scratch, chain->item_id, rec->record_id);
			continue;
		}
		/* Everything else runs its descriptor's hook below.  The families that
		 * never become chain members (compositing state, structural, analysis)
		 * were filtered long before this loop; they are listed to keep the
		 * switch exhaustive, not because they can arrive here. */
		case NDE_OPC_PIXEL:
		case NDE_OPC_MASK:
		case NDE_OPC_DOCUMENT:
		case NDE_OPC_COMPOSITING:
		case NDE_OPC_STRUCTURAL:
		case NDE_OPC_ANALYSIS:
		case NDE_OPC_UNKNOWN:
			break;
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
				destroy_any_args(user);
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): replay preparation failed"),
				                       rec->record_id, rec->op_id);
				goto fail;
			}
		}
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_any_args(user);
			*err = g_strdup(_("out of memory"));
			goto fail;
		}
		/* Restore the record's mask input before running it.  Without this
		 * the op would run unmasked and the replay would claim to have
		 * reproduced a record it did not. */
		if (!nde_mask_pin_install(scratch, rec, err)) {
			destroy_any_args(user);
			free(args);
			goto fail;
		}
		args->fit = scratch;       /* PRIVATE fits — the whole point */
		args->op = op;
		args->user = user;
		args->nde_replay = TRUE;
		args->mask_aware = (scratch->mask != NULL);
		args->layer_pos_x = pos_x;   /* NULL unless the chain has geometry */
		args->layer_pos_y = pos_y;
		args->max_threads = com.max_thread;
		int rc = GPOINTER_TO_INT(generic_image_worker(args));
		free_generic_img_args(args);   /* frees user via its destructor too */
		nde_mask_pin_clear(scratch);
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
	replay_current_record = 0;
	return scratch;

fail:
	replay_current_record = 0;
	clearfits(scratch);
	free(scratch);
	return NULL;
}

/* The state of @item_id just after @upto_record_id (0 = its whole chain).
 * A LIVE edge: the value is re-derived by replaying the item rather than read
 * from a stored copy, which is what lets an amend upstream of it take effect.
 * Caller owns the result. */
/* The fits a chain edit commits into: gfit for plain images, the layer's
 * fit for FLIS items (same pointer when it is the active layer). */
fits *nde_edit_target_fits(gint item_id) {
	if (item_id < 0)
		return gfit;
	flis_layer_t *lay = flis_layer_get_by_id(item_id);
	return lay ? lay->fit : NULL;
}

/* Shared core of amend (new_params != NULL) and delete (new_params == NULL).
 * Conductor context: the caller holds SLOT_REPLAY, so capture, undo and
 * python cannot interleave with the replay or the commit. */
/* Checkpoint-built replay results carry pixels + WCS only (nde_snap stores
 * no header text, keywords or HISTORY).  A history edit rewrites pixel
 * parameters, not the document's identity: after a commit swap, move the
 * pre-edit metadata (keywords, header text, unknown keys, FITS HISTORY)
 * from the superseded fits back onto the target, keeping the REPLAYED WCS
 * and focal length WHEN THE REPLAY PRODUCED ONE — the chain may legitimately
 * have changed those, but a replay carrying no WCS has not un-solved the
 * image and must not clear the solve it already had.
 * @old is the pre-edit fits the swap left holding the rich metadata; it is
 * about to be cleared, so ownership moves are plain pointer swaps. */
void nde_commit_restore_metadata(fits *target, fits *old) {
	gboolean quiesced = nde_commit_lock(target);
	/* Move the replayed WCS aside, swap the keyword structs wholesale, then
	 * put the replayed WCS back.  The superseded wcslib travels to @old inside
	 * its keywords struct, so clearfits(@old) frees it properly.
	 *
	 * ONLY when the replay actually produced a WCS, though.  A replay that
	 * carries none did not UN-SOLVE the image: composite results
	 * (composite_apply) are pixels only, because there is no plate solve to
	 * derive from a blend.  Keeping the replayed NULL silently dropped the
	 * merged image's solve on the FIRST amend, and the next edit needing a WCS
	 * donor — the group-calibration replay hunting for a solved layer to
	 * donate one to its composite — then failed outright with "no WCS data or
	 * it is not supported" and fell back to the stored calibration. */
	const gboolean replay_solved = target->keywords.wcslib != NULL;
	wcs_info replay_wcsdata      = target->keywords.wcsdata;
	struct wcsprm *replay_wcslib = target->keywords.wcslib;
	double replay_flen           = target->keywords.focal_length;
	fkeywords rich = old->keywords;
	old->keywords = target->keywords;         /* replay-minimal set... */
	/* ...plus whichever WCS is being superseded (none, when we keep @rich's). */
	old->keywords.wcslib = replay_solved ? rich.wcslib : NULL;
	target->keywords = rich;
	if (replay_solved) {
		target->keywords.wcsdata      = replay_wcsdata;
		target->keywords.wcslib       = replay_wcslib;
		target->keywords.focal_length = replay_flen;
	}

	char *h = target->header;
	target->header = old->header;
	old->header = h;
	gchar *u = target->unknown_keys;
	target->unknown_keys = old->unknown_keys;
	old->unknown_keys = u;
	GSList *hist = target->history;
	target->history = old->history;
	old->history = hist;
	nde_commit_unlock(target, quiesced);
}

/* As resolve_item_state, additionally reporting where the replayed chain left
 * the input's layer.  *pos_valid is set TRUE only when the chain carries an
 * offset (a geometry member with a stored start offset — the same condition
 * nde_replay_start_offset applies everywhere else); otherwise the recorded
 * capture-time position remains the only truth and the outs are untouched. */
static fits *resolve_item_state_pos(gint item_id, gint64 upto_record_id,
                                    gint *pos_x, gint *pos_y,
                                    gboolean *pos_valid, gchar **err) {
	return resolve_item_state_pos_bound(item_id, upto_record_id, FALSE,
	                                    pos_x, pos_y, pos_valid, err);
}

/* The bound-aware core.  @exclusive says which side of the anchor the prefix
 * ends on: FALSE replays up to AND INCLUDING the anchor's log position (a
 * pin's "state after this record"), TRUE stops JUST BEFORE it (a joint
 * record resolving a sibling's pre-record state, nde_joint.h — inclusive
 * would replay the joint record itself and recurse forever). */
static fits *resolve_item_state_pos_bound(gint item_id, gint64 upto_record_id,
                                          gboolean exclusive,
                                          gint *pos_x, gint *pos_y,
                                          gboolean *pos_valid, gchar **err) {
	nde_chain *c = nde_chain_build(item_id);
	/* A pin's src_record_id of 0 means the item's BASELINE, not "all of it"
	 * (nde_history.h) — the state before anything was recorded against it.
	 * Reading it as "the whole chain" made a mask pinned to the untouched
	 * image re-derive from a later state once the image gained records.
	 *
	 * The pin names a record in the item's HISTORY, and that is a wider set
	 * than its pixel chain: a compositing-state change (opacity, blend,
	 * visibility) or a mask step can be the last record against an item
	 * without being a member here, and a member can later be deleted.  The
	 * meaning is positional either way — "the state as of that point in the
	 * log" — so take the prefix of members recorded at or before it rather
	 * than demanding an exact hit and failing when there is none.
	 *
	 * POSITIONAL MEANS POSITION, not id.  Comparing ids assumes they increase
	 * down the log, and insertion breaks that: a step inserted BEFORE the
	 * pinned record has a HIGHER id, so an id comparison stopped the prefix
	 * dead at it and an operation inserted into a consumed input changed
	 * nothing at all.  The chain's members are views into a snapshot held in
	 * log order, so the snapshot is what says which came first. */
	guint upto = 0;
	if (upto_record_id) {
		gint anchor_pos = -1;
		for (guint i = 0; c->snapshot && i < c->snapshot->len; i++) {
			const nde_record *r = g_ptr_array_index(c->snapshot, i);
			if (r->record_id == upto_record_id) {
				anchor_pos = (gint)i;
				break;
			}
		}
		for (guint i = 0; i < c->records->len; i++) {
			const nde_record *m = g_ptr_array_index(c->records, i);
			if (anchor_pos < 0) {
				/* The pin names a record the log no longer holds; with no
				 * position to compare against, fall back to the id. */
				if (m->record_id > upto_record_id ||
				    (exclusive && m->record_id == upto_record_id))
					break;
			} else {
				gint member_pos = -1;
				for (guint j = 0; c->snapshot && j < c->snapshot->len; j++) {
					if (g_ptr_array_index(c->snapshot, j) == m) {
						member_pos = (gint)j;
						break;
					}
				}
				if (member_pos > anchor_pos ||
				    (exclusive && member_pos == anchor_pos))
					break;
			}
			upto = i + 1;
		}
	}
	/* Only the members actually replayed matter: [0..upto).  Demanding
	 * whole-chain replayability here refused every pin into a chain with a
	 * barrier anywhere — including barriers far past the pin — and did so
	 * with an EMPTY reason list when the freeze cause was a checkpointed
	 * barrier (those add no reasons, they only move tail_start). */
	if (upto > c->first_block) {
		GString *m = g_string_new(NULL);
		g_string_append_printf(m, _("the source of this input cannot be rebuilt: "));
		if (!c->reasons->len)
			g_string_append(m, _("an earlier step on it cannot be replayed"));
		for (guint i = 0; i < c->reasons->len; i++) {
			if (i)
				g_string_append(m, "; ");
			g_string_append(m, g_ptr_array_index(c->reasons, i));
		}
		*err = g_string_free(m, FALSE);
		nde_chain_free(c);
		return NULL;
	}
	/* A barrier member inside the prefix cannot itself be re-run, but it left
	 * its post-op pixels as a restart point (a checkpoint-less barrier is a
	 * hard blocker, refused above).  Restart from the last one before the
	 * pin and replay only what follows it. */
	gint64 prefix_restart = 0;
	guint prefix_start = 0;
	for (guint i = 0; i < upto; i++) {
		if (g_array_index(c->member_flags, guint8, i) & NDE_CHAIN_MEMBER_BARRIER) {
			const nde_record *b = g_ptr_array_index(c->records, i);
			prefix_restart = b->record_id;
			prefix_start = i + 1;
		}
	}
	if (prefix_restart) {
		fits *start = nde_checkpoint_output_get(prefix_restart);
		if (!start) {
			*err = g_strdup_printf(_("the source of this input cannot be rebuilt: "
			                         "the stored pixels of step %" G_GINT64_FORMAT
			                         " are no longer kept"), prefix_restart);
			nde_chain_free(c);
			return NULL;
		}
		gint px = 0, py = 0;
		gboolean carry = nde_replay_start_offset(c, prefix_restart, &px, &py);
		fits *out = nde_replay_apply_records(start, c, prefix_start, upto,
		                                 carry ? &px : NULL, carry ? &py : NULL,
		                                 err);
		if (out && carry && pos_x && pos_y) {
			*pos_x = px;
			*pos_y = py;
			if (pos_valid)
				*pos_valid = TRUE;
		}
		nde_chain_free(c);
		return out;
	}
	/* Born of a merge: there is no baseline and none is wanted — the first
	 * member renders the inputs it was given.  The same exemption its
	 * siblings have (nde_chain_replay, nde_replay_chain_restart_state, recompute_item);
	 * without it a composite-born input whose baseline had been evicted
	 * failed here while nde_chain_replay on the same item succeeded. */
	if (c->from_composite) {
		if (upto == 0) {
			/* The pin names a state before the item's origin — nothing can
			 * produce it. */
			*err = g_strdup(_("the input predates the merge that created it"));
			nde_chain_free(c);
			return NULL;
		}
		fits *out = nde_replay_apply_records(NULL, c, 0, upto, NULL, NULL, err);
		nde_chain_free(c);
		return out;
	}
	fits *start = nde_checkpoint_baseline_get(item_id);
	if (!start && c->records->len == 0) {
		/* Nothing was ever recorded against this item, so no baseline was
		 * ever taken — and none is needed: with no operations to undo, its
		 * current pixels ARE its original state.  This is the ordinary case
		 * for a mask built as the very first thing done to an image. */
		fits *live = nde_edit_target_fits(item_id);
		if (live) {
			start = calloc(1, sizeof(fits));
			if (start) {
				g_rw_lock_reader_lock(&live->rwlock);
				int rc = copyfits(live, start, CP_DEEPCOPY | CP_ALLOC, -1);
				g_rw_lock_reader_unlock(&live->rwlock);
				if (rc) {
					free(start);
					start = NULL;
				}
			}
		}
	}
	if (!start) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		nde_chain_free(c);
		return NULL;
	}
	gint px = 0, py = 0;
	gboolean carry = nde_replay_start_offset(c, 0, &px, &py);
	fits *out = nde_replay_apply_records(start, c, 0, upto,
	                                 carry ? &px : NULL, carry ? &py : NULL,
	                                 err);
	if (out && carry && pos_x && pos_y) {
		*pos_x = px;
		*pos_y = py;
		if (pos_valid)
			*pos_valid = TRUE;
	}
	nde_chain_free(c);
	return out;
}

static fits *resolve_item_state(gint item_id, gint64 upto_record_id, gchar **err) {
	return resolve_item_state_pos(item_id, upto_record_id, NULL, NULL, NULL, err);
}

fits *nde_replay_resolve_before(gint item_id, gint64 before_record_id, gchar **err) {
	g_return_val_if_fail(before_record_id != 0, NULL);
	return resolve_item_state_pos_bound(item_id, before_record_id, TRUE,
	                                    NULL, NULL, NULL, err);
}

/* Replay members [0..upto) of a MASK chain.  The result is a fits carrying
 * the rebuilt mask in ->mask; its pixels are whatever the chain started from
 * and are not themselves meaningful.  Mask ops go straight to their hook
 * rather than through generic_mask_worker: everything the worker adds —
 * undo entries, lmask routing, capture, idles — is precisely what a replay
 * must not do. */
fits *nde_mask_chain_replay(const nde_chain *chain, guint upto, gchar **err) {
	g_return_val_if_fail(chain->records->len > 0, NULL);
	fits *scratch = NULL;
	/* Which image state @scratch currently holds, so a run of records reading
	 * the same one costs a single replay. */
	gint     src_item = 0;
	gint64   src_rec  = 0;
	gboolean have_src = FALSE;

	for (guint i = 0; i < upto; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (!processing_should_continue()) {
			*err = g_strdup(_("cancelled"));
			goto fail;
		}
		/* EVERY record that reads the image reads it AT ITS OWN recorded
		 * point.  Only the first record's pin used to be resolved, so a mask
		 * rebuilt later in its own history — mask.from_stars over a stretch
		 * added after mask.from_channel — silently read the pixels the FIRST
		 * record had named, and the replayed mask disagreed with the one the
		 * user had actually built. */
		const nde_input_pin *img = nde_record_input(rec, "image");
		if (img && (!have_src || img->src_item_id != src_item ||
		            img->src_record_id != src_rec)) {
			fits *next = resolve_item_state(img->src_item_id, img->src_record_id, err);
			if (!next)
				goto fail;
			/* The chain rebuilds the mask from nothing: whatever the resolved
			 * image happened to carry is not part of this chain's value. */
			if (next->mask) {
				free_mask(next->mask);
				next->mask = NULL;
			}
			next->mask_active = FALSE;
			/* But the mask built SO FAR crosses the change of source when it
			 * still fits.  A from-image op replaces it in practice; nothing in
			 * the model says one has to, and losing it here would be silent. */
			if (scratch && scratch->mask &&
			    next->rx == scratch->rx && next->ry == scratch->ry) {
				next->mask = scratch->mask;
				scratch->mask = NULL;
				next->mask_active = scratch->mask_active;
			}
			if (scratch) {
				clearfits(scratch);
				free(scratch);
			}
			scratch  = next;
			src_item = img->src_item_id;
			src_rec  = img->src_record_id;
			have_src = TRUE;
		}
		if (!scratch) {
			/* Built from a file the hook loads itself; it still needs
			 * somewhere to put the mask, and the file decides the size. */
			scratch = calloc(1, sizeof(fits));
			if (!scratch) {
				*err = g_strdup(_("out of memory"));
				return NULL;
			}
		}
		const op_descriptor *op = op_descriptor_by_id(rec->op_id);
		gpointer user = op ? op->deserialize(rec->params, rec->op_version) : NULL;
		if (!user) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
			                       rec->record_id, rec->op_id ? rec->op_id : "?");
			goto fail;
		}
		struct generic_mask_args margs = { 0 };
		margs.fit = scratch;
		margs.op = op;
		margs.user = user;
		margs.command = TRUE;      /* no idles, no GUI completion */
		margs.max_threads = com.max_thread;
		op_descriptor_fill_mask_args(&margs);
		int rc = margs.mask_hook(&margs);
		destroy_any_args(user);
		if (rc) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) failed to apply"),
			                       rec->record_id, rec->op_id ? rec->op_id : "?");
			goto fail;
		}
	}
	return scratch;

fail:
	if (scratch) {
		clearfits(scratch);
		free(scratch);
	}
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
	if (chain->is_mask)
		return nde_mask_chain_replay(chain, chain->records->len, err);
	/* Born of a merge: there is no baseline to start from, and none is
	 * wanted — the first member renders the inputs it was given. */
	fits *scratch = chain->from_composite ?
			NULL : nde_checkpoint_baseline_get(chain->item_id);
	if (!scratch && !chain->from_composite) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		return NULL;
	}
	/* Verification path: the offset does not affect a single pixel (it is a
	 * pure side output of the geometry hooks), so there is nothing to carry. */
	return nde_replay_apply_records(scratch, chain, 0, chain->records->len,
	                            NULL, NULL, err);
}

/* The state a chain's tail restarts from: the barrier checkpoint when there is
 * one, otherwise the baseline.  An item born of a composite has neither when
 * its tail begins at that composite — its origin IS the first member, so the
 * restart state is legitimately NULL and *err is left unset.  Callers must
 * therefore test @err, not the return value. */
fits *nde_replay_chain_restart_state(const nde_chain *chain, gchar **err) {
	if (chain->restart_ckpt_id > 0) {
		fits *f = nde_checkpoint_output_get(chain->restart_ckpt_id);
		if (!f)
			*err = g_strdup(_("failed to load the barrier checkpoint"));
		return f;
	}
	if (chain->from_composite && chain->tail_start == 0)
		return NULL;
	fits *f = nde_checkpoint_baseline_get(chain->item_id);
	if (!f)
		*err = g_strdup(_("failed to load the baseline checkpoint"));
	return f;
}

fits *nde_chain_replay_tail(const nde_chain *chain, gchar **err) {
	g_return_val_if_fail(chain != NULL, NULL);
	g_return_val_if_fail(err != NULL, NULL);
	*err = NULL;
	if (!chain->tail_replayable) {
		*err = g_strdup(_("the editable tail is not replayable"));
		return NULL;
	}
	fits *start = nde_replay_chain_restart_state(chain, err);
	if (!start && *err)
		return NULL;
	return nde_replay_apply_records(start, chain, chain->tail_start, chain->records->len,
	                            NULL, NULL, err);
}
