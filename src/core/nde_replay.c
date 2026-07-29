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
#include "core/masks.h"
#include "core/nde_compositing.h"
#include "core/nde_composite.h"
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
static gboolean mask_pin_install(fits *scratch, const nde_record *rec, gchar **err) {
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
static void mask_pin_clear(fits *scratch) {
	if (scratch->mask) {
		free_mask(scratch->mask);
		scratch->mask = NULL;
	}
	scratch->mask_active = FALSE;
}

/* TRUE when the record's params pin an external image file (Convention 1) —
 * a mask built from one needs no image input, because the hook loads it. */
static gboolean record_names_a_file(const nde_record *rec) {
	if (!rec->params)
		return FALSE;
	GHashTable *kv = nde_kv_parse(rec->params);
	const char *path = nde_kv_get_str(kv, "operand_path");
	gboolean has = path && *path;
	g_hash_table_unref(kv);
	return has;
}

static gchar *member_invalid_reason(const nde_record *rec) {
	if (rec->tier == NDE_TIER_C)
		return tier_c_invalid_reason(rec);
	/* A composite node has no op descriptor: its "deserializer" is
	 * nde_composite_params_decode, and chain membership already required it to
	 * succeed (nde_composite_record_replayable). */
	if (nde_composite_is_op(rec->op_id))
		return NULL;
	if (rec->tier != NDE_TIER_A)
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque — not replayable"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
	/* A masked record used to be a barrier outright: the flag said a mask had
	 * been active but nothing said WHICH, so there was nothing to reproduce.
	 * A pin that resolves to a stored state gives the replay the mask back;
	 * only an unpinned or no-longer-stored one still freezes the chain.
	 * (Tier C keeps the blanket refusal above: a re-run script decides for
	 * itself what to do with a mask, so installing one proves nothing.) */
	if (rec->mask_active && !nde_mask_pin_resolvable(rec))
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) ran with a mask whose pixels were not kept — not replayable"),
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

/* POLICY predicates — see nde_replay.h.  Liveness and trial-chain
 * replayability are the execute path's job, not these. */
gboolean nde_record_amendable(const nde_record *rec) {
	if (!rec || rec->tier != NDE_TIER_A)
		return FALSE;
	/* Compositing-state records carry editable params but no op descriptor:
	 * nde_compositing_validate() is their deserializer (see nde_compositing.h). */
	if (nde_compositing_is_op(rec->op_id))
		return TRUE;
	const op_descriptor *op = op_descriptor_by_id(rec->op_id);
	return op && op->deserialize;
}

gboolean nde_record_deletable(const nde_record *rec) {
	if (!rec)
		return FALSE;
	if (rec->scope == NDE_SCOPE_DOCUMENT)
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
		/* Compositing-state records are inputs to the compositor, not pixel
		 * operations: neither chain members nor blockers (nde_compositing.h). */
		if (nde_compositing_is_op(rec->op_id))
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
				/* Geometry on a layer moves it on the canvas as well as
				 * changing its pixels, and each move is relative to where
				 * the layer already was.  So the replay has to start from a
				 * known position: with one recorded against the baseline
				 * this is an ordinary member, without one there is nothing
				 * to anchor to (graph step 5). */
				if (nde_checkpoint_baseline_get_offset(item_id, NULL, NULL)) {
					member = TRUE;
					chain->has_geometry = TRUE;
				} else {
					add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) moves this layer on the canvas, but no starting position was recorded — not replayable"),
					           rec->record_id, rec->op_id ? rec->op_id : "?");
				}
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
				if (nde_composite_record_replayable(rec)) {
					/* A composite node (graph step 7): its inputs are
					 * pinned and its per-input state recorded, so the
					 * merge re-runs like any other member and the steps
					 * before it stay editable. */
					member = TRUE;
					chain->has_composite = TRUE;
					break;
				}
				/* Recorded before step 7, or with an input this build
				 * cannot rebuild (a masked overlay): the item's pixels
				 * were replaced by a composite and nothing here can
				 * reproduce it. */
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

	/* A MASK item's chain is recognised by its ops, not by its id: a mask op
	 * is one with a mask_hook.  It has no baseline of its own — a mask is
	 * DERIVED, and its chain starts at the image its first member read
	 * (that record's "image" pin) or at the file that member names. */
	if (chain->records->len) {
		const nde_record *first = g_ptr_array_index(chain->records, 0);
		const op_descriptor *fop = op_descriptor_by_id(first->op_id);
		chain->is_mask = fop && fop->mask_hook;
	}
	if (chain->is_mask) {
		const nde_record *first = g_ptr_array_index(chain->records, 0);
		if (!nde_record_input(first, "image") && !record_names_a_file(first))
			add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) does not say what this mask was built from — not replayable"),
			           first->record_id, first->op_id ? first->op_id : "?");
	} else if (chain->records->len && !nde_checkpoint_baseline_exists(item_id)) {
		add_reason(chain, _("no baseline checkpoint — the file predates baselines, or the history began before this build"));
	}

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

/* Where a replay's layer starts.  @restart_id 0 means the baseline; anything
 * else is that record's output checkpoint.  FALSE when the chain carries no
 * geometry member (a plain image, or a layer nothing ever moved), in which
 * case the hooks are handed NULL and move nothing. */
static gboolean replay_start_offset(const nde_chain *chain, gint64 restart_id,
                                    gint *pos_x, gint *pos_y) {
	if (!chain->has_geometry)
		return FALSE;
	if (restart_id > 0 &&
	    nde_checkpoint_output_get_offset(restart_id, pos_x, pos_y))
		return TRUE;
	return nde_checkpoint_baseline_get_offset(chain->item_id, pos_x, pos_y);
}

/* Move the layer to where the replay left it.  The pixels are committed by
 * the caller's swap; without this the layer would keep the position its
 * pre-edit geometry produced, which the new pixels no longer match. */
static void commit_layer_offset(gint item_id, gint pos_x, gint pos_y) {
	if (item_id < 0)
		return;
	flis_layer_t *lay = flis_layer_get_by_id(item_id);
	if (!lay)
		return;
	lay->position_x = pos_x;
	lay->position_y = pos_y;
	gui_iface.flis_invalidate_composite();
}

static fits *resolve_item_state(gint item_id, gint64 upto_record_id, gchar **err);

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
			/* A geometry step before the composite moves this input, so an
			 * amend of that step must move where it lands here too. */
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
		pixels[i] = resolve_item_state(pin->src_item_id, pin->src_record_id, err);
		owned[i] = TRUE;
		ok = pixels[i] != NULL;
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
static fits *replay_apply_records(fits *scratch, const nde_chain *chain,
                                  guint from, guint upto,
                                  gint *pos_x, gint *pos_y, gchar **err) {
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
		if (nde_composite_is_op(rec->op_id)) {
			fits *merged = composite_apply(scratch, rec, chain->item_id,
			                               pos_x, pos_y, err);
			if (!merged)
				goto fail;   /* scratch untouched — the fail path frees it */
			clearfits(scratch);
			free(scratch);
			scratch = merged;
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
		/* Restore the record's mask input before running it.  Without this
		 * the op would run unmasked and the replay would claim to have
		 * reproduced a record it did not. */
		if (!mask_pin_install(scratch, rec, err)) {
			destroy_user(user);
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
		mask_pin_clear(scratch);
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

/* The state of @item_id just after @upto_record_id (0 = its whole chain).
 * A LIVE edge: the value is re-derived by replaying the item rather than read
 * from a stored copy, which is what lets an amend upstream of it take effect.
 * Caller owns the result. */
static fits *edit_target_fits(gint item_id);
static void commit_restore_metadata(fits *target, fits *old);

static fits *resolve_item_state(gint item_id, gint64 upto_record_id, gchar **err) {
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
	 * than demanding an exact hit and failing when there is none. */
	guint upto = 0;
	if (upto_record_id) {
		for (guint i = 0; i < c->records->len; i++) {
			const nde_record *r = g_ptr_array_index(c->records, i);
			if (r->record_id > upto_record_id)
				break;
			upto = i + 1;
		}
	}
	if (!c->replayable) {
		GString *m = g_string_new(NULL);
		g_string_append_printf(m, _("the source of this input cannot be rebuilt: "));
		for (guint i = 0; i < c->reasons->len; i++) {
			if (i)
				g_string_append(m, "; ");
			g_string_append(m, g_ptr_array_index(c->reasons, i));
		}
		*err = g_string_free(m, FALSE);
		nde_chain_free(c);
		return NULL;
	}
	fits *start = nde_checkpoint_baseline_get(item_id);
	if (!start && c->records->len == 0) {
		/* Nothing was ever recorded against this item, so no baseline was
		 * ever taken — and none is needed: with no operations to undo, its
		 * current pixels ARE its original state.  This is the ordinary case
		 * for a mask built as the very first thing done to an image. */
		fits *live = edit_target_fits(item_id);
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
	fits *out = replay_apply_records(start, c, 0, upto, NULL, NULL, err);
	nde_chain_free(c);
	return out;
}

/* Replay members [0..upto) of a MASK chain.  The result is a fits carrying
 * the rebuilt mask in ->mask; its pixels are whatever the chain started from
 * and are not themselves meaningful.  Mask ops go straight to their hook
 * rather than through generic_mask_worker: everything the worker adds —
 * undo entries, lmask routing, capture, idles — is precisely what a replay
 * must not do. */
static fits *mask_chain_replay(const nde_chain *chain, guint upto, gchar **err) {
	g_return_val_if_fail(chain->records->len > 0, NULL);
	const nde_record *first = g_ptr_array_index(chain->records, 0);
	const nde_input_pin *img = nde_record_input(first, "image");
	fits *scratch = NULL;
	if (img) {
		scratch = resolve_item_state(img->src_item_id, img->src_record_id, err);
		if (!scratch)
			return NULL;
	} else {
		/* Built from a file the hook loads itself; it still needs somewhere
		 * to put the mask, and the file decides the size. */
		scratch = calloc(1, sizeof(fits));
		if (!scratch) {
			*err = g_strdup(_("out of memory"));
			return NULL;
		}
	}
	/* The chain rebuilds the mask from nothing: drop anything the starting
	 * image happened to carry. */
	if (scratch->mask) {
		free_mask(scratch->mask);
		scratch->mask = NULL;
	}
	scratch->mask_active = FALSE;

	for (guint i = 0; i < upto; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		if (!processing_should_continue()) {
			*err = g_strdup(_("cancelled"));
			goto fail;
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
		destroy_user(user);
		if (rc) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) failed to apply"),
			                       rec->record_id, rec->op_id ? rec->op_id : "?");
			goto fail;
		}
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
	if (chain->is_mask)
		return mask_chain_replay(chain, chain->records->len, err);
	fits *scratch = nde_checkpoint_baseline_get(chain->item_id);
	if (!scratch) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		return NULL;
	}
	/* Verification path: the offset does not affect a single pixel (it is a
	 * pure side output of the geometry hooks), so there is nothing to carry. */
	return replay_apply_records(scratch, chain, 0, chain->records->len,
	                            NULL, NULL, err);
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
	return replay_apply_records(start, chain, chain->tail_start, chain->records->len,
	                            NULL, NULL, err);
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
	/* The pool is a PIXEL cache: its states carry no layer position, so a
	 * chain that moves the layer cannot restart from one — there would be
	 * nothing to anchor the geometry to.  Fall straight through to the
	 * checkpoint restart, which does record a position (graph step 5). */
	if (chain->has_geometry)
		e = chain->tail_start;
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
	gint pos_x = 0, pos_y = 0;
	gboolean carry = replay_start_offset(chain, 0, &pos_x, &pos_y);
	fits *start = nde_checkpoint_baseline_get(item_id);
	if (!start) {
		*err = g_strdup(_("failed to load the baseline checkpoint"));
		nde_chain_free(chain);
		return FALSE;
	}
	nde_snapstore_invalidate_from(item_id, 0);
	fits *result = replay_apply_records(start, chain, 0, chain->records->len,
	                                    carry ? &pos_x : NULL,
	                                    carry ? &pos_y : NULL, err);
	nde_chain_free(chain);
	if (!result)
		return FALSE;

	fits *target = edit_target_fits(item_id);
	if (!target) {
		*err = g_strdup(_("the target layer no longer exists"));
		clearfits(result);
		free(result);
		return FALSE;
	}
	g_rw_lock_writer_lock(&target->rwlock);
	fits_swap_all_except_rwlock(target, result);
	g_rw_lock_writer_unlock(&target->rwlock);
	commit_restore_metadata(target, result);
	clearfits(result);
	free(result);
	if (carry)
		commit_layer_offset(item_id, pos_x, pos_y);
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
static gboolean item_is_retained_input(gint item_id) {
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
			if (item_is_retained_input(consumer)) {
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

/* Re-derive @mask_item's mask as it stood after chain member @upto-1 and
 * replace the stored copy its consumers read (nde_replay.h's pin store). */
static gboolean refresh_pinned_mask(const nde_chain *mask_chain, guint upto,
                                    gint64 coord_record, gint mask_item,
                                    gchar **err) {
	fits *built = mask_chain_replay(mask_chain, upto, err);
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
		const nde_input_pin *pin = nde_record_input_by_item(rec, mask_item);
		if (!pin)
			continue;
		gpointer p = g_hash_table_lookup(pos_of, &pin->src_record_id);
		if (!p)
			continue;   /* pinned to a record that is no longer in the chain */
		guint upto = GPOINTER_TO_UINT(p);
		if (upto <= from_pos)
			continue;   /* the edit is after this pin: its mask is unchanged */
		if (!g_hash_table_contains(done_coords, &pin->src_record_id)) {
			gchar *err = NULL;
			if (refresh_pinned_mask(mask_chain, upto, pin->src_record_id,
			                        mask_item, &err)) {
				refreshed++;
				gint64 *k = g_new(gint64, 1);
				*k = pin->src_record_id;
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
		if (item_is_retained_input(item)) {
			cascade_composite_consumers(item);
			continue;
		}
		if (recompute_item(item, &err)) {
			redone++;
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

/* Why a history edit is refused right now.  Both modes install a synthesized
 * state over the target's real pixels, so a second edit would compute against
 * something that is not the committed image. */
static const char *edit_in_progress_reason(void) {
	return nde_edit_at_active() ?
			_("an insertion point is open — finish or cancel it first") :
			_("another history step is being edited — close its dialog first");
}

static gboolean edit_execute(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;

	if (amend_preview_installed()) {
		*err = g_strdup(edit_in_progress_reason());
		return FALSE;
	}

	/* Locate the record among the LIVE records and pre-check it, for
	 * user-facing errors before any heavy work.  (nde_history_amend/delete
	 * re-validate under the mutex at commit time.) */
	gint item_id = 0;
	gboolean found = FALSE, is_compositing = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->record_id == record_id) {
			found = TRUE;
			item_id = rec->target_item_id;
			is_compositing = nde_compositing_is_op(rec->op_id);
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
		gboolean log_ok = new_params ? nde_history_amend(record_id, new_params, err)
		                             : nde_history_delete(record_id, err);
		if (!log_ok)
			return FALSE;
		if (!nde_compositing_recompute(item_id, err))
			return FALSE;
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
		fits *built = mask_chain_replay(chain, chain->records->len, err);
		nde_chain_free(chain);
		if (!built) {
			gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
			return FALSE;
		}
		/* A mask whose layer was composited away has nowhere to be committed
		 * — the layer it belonged to is gone.  The replay above still had to
		 * run, to prove the edited chain applies; what shows the change is the
		 * composite that kept a copy of this mask, refreshed by the cascade
		 * below. */
		gboolean ok = item_is_retained_input(item_id) ||
		              commit_mask_value(item_id, built, err);
		clearfits(built);
		free(built);
		if (!ok) {
			gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
			return FALSE;
		}
		gboolean log_ok = new_params ? nde_history_amend(record_id, new_params, err)
		                             : nde_history_delete(record_id, err);
		if (!log_ok) {
			gui_iface.set_progress(PROGRESS_RESET, _("The mask was rebuilt but the history could not be updated"));
			return FALSE;
		}
		/* Reverse invalidation: the mask changed, so anything that used it
		 * is now showing a result built from the old one.  Runs AFTER the
		 * log commit, because the consumers are recomputed from the log. */
		cascade_mask_consumers(item_id, mask_pos);
		undo_flush();
		gui_iface.invalidate_histogram();
		if (is_current_image_flis())
			gui_iface.flis_invalidate_composite();
		notify_gfit_data_modified();
		gui_iface.set_progress(PROGRESS_DONE, _("Edit history updated"));
		return TRUE;
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
	gint pos_x = 0, pos_y = 0;
	gboolean carry = replay_start_offset(chain, chain->restart_ckpt_id, &pos_x, &pos_y);
	fits *result = start ? replay_apply_records(start, chain, start_idx, chain->records->len,
	                                            carry ? &pos_x : NULL,
	                                            carry ? &pos_y : NULL, err) : NULL;
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
	gboolean retained = FALSE;
	if (item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(item_id);
		target = lay ? lay->fit : NULL;
		retained = !lay && item_is_retained_input(item_id);
	}
	if (retained) {
		/* A retained input has no layer to swap into: the replay above was
		 * run to prove the edited chain still applies, and its pixels are an
		 * intermediate value that only the composite consumes.  Commit the
		 * log, then recompute the consumers — which resolve this item by
		 * replaying it again, now from the amended log. */
		clearfits(result);
		free(result);
		gboolean log_ok = new_params ? nde_history_amend(record_id, new_params, err)
		                             : nde_history_delete(record_id, err);
		if (!log_ok) {
			gui_iface.set_progress(PROGRESS_RESET, _("Edit failed — nothing was changed"));
			return FALSE;
		}
		cascade_composite_consumers(item_id);
		undo_flush();
		gui_iface.invalidate_histogram();
		notify_gfit_data_modified();
		gui_iface.set_progress(PROGRESS_DONE, _("Edit history updated"));
		return TRUE;
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
	if (carry)
		commit_layer_offset(item_id, pos_x, pos_y);

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
	gint pos_x = 0, pos_y = 0;
	gboolean carry = replay_start_offset(chain, chain->restart_ckpt_id, &pos_x, &pos_y);
	fits *result = start ? replay_apply_records(start, chain, start_idx, chain->records->len,
	                                            carry ? &pos_x : NULL,
	                                            carry ? &pos_y : NULL, err) : NULL;
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
	if (carry)
		commit_layer_offset(item_id, pos_x, pos_y);

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
	gboolean insert;     /* edit-at mode: new steps are inserted before K */
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
	apv.insert = FALSE;
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
	fits *start = resolve_edit_restart(chain, (guint)e, record_id, &start_idx, err);
	/* Preview only: the true pixels (and the layer's real position) come back
	 * untouched on exit, so nothing is committed and nothing is carried. */
	fits *pre_k = start ? replay_apply_records(start, chain, start_idx, (guint)e,
	                                           NULL, NULL, err) : NULL;
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
	fits *saved = apv.saved;
	apv_clear_state_locked();
	g_mutex_unlock(&apv_mutex);

	gboolean disturbed = nde_history_insert_point_disturbed();
	GArray *inserted = nde_history_insert_point_clear();

	/* Restore the true pixels and metadata first, exactly as the amend
	 * preview does.  `saved` comes back holding the pre-anchor state PLUS
	 * whatever was inserted — the starting point for the forward replay. */
	fits *target = edit_target_fits(item_id);
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
	fits *result = NULL;
	if (k < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " is no longer part of this image's history"),
		                       anchor_id);
	} else if ((guint)k < chain->tail_start || !chain->tail_replayable) {
		/* An inserted opaque step became a barrier that swallowed the anchor
		 * (it has no output checkpoint, or one of the later steps cannot
		 * replay from it). */
		*err = g_strdup(_("the steps after the insertion point can no longer be recomputed"));
	} else {
		gui_iface.set_progress(0.f, _("Recomputing edit history..."));
		/* Cached states at or after the insertion describe the pre-insert
		 * lineage; the replay re-deposits fresh ones as it goes. */
		nde_snapstore_invalidate_from(item_id, first_inserted);
		/* No offset to carry: a chain that moves the layer is refused an
		 * insertion point in the first place (apv_begin_execute) — the
		 * position to restart from is not the layer's current one, which
		 * already embodies every record after the anchor. */
		result = replay_apply_records(saved, chain, (guint)k, chain->records->len,
		                              NULL, NULL, err);
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

	g_rw_lock_writer_lock(&target->rwlock);
	fits_swap_all_except_rwlock(target, result);
	g_rw_lock_writer_unlock(&target->rwlock);
	commit_restore_metadata(target, result);
	clearfits(result);
	free(result);

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
