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
 * NDE chain building and replayability policy.
 *
 * An item's history is a list of records; its CHAIN is the subsequence of them
 * that describes its pixels, together with the verdict on whether that
 * subsequence can be replayed and where its editable tail starts.  Everything
 * here answers a question ABOUT records without running one: may this record
 * be re-run, may the user amend or delete it, what blocks the ones before it.
 * The engine that does re-run them is nde_replay.c.
 */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/op_descriptor.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_op_class.h"
#include "core/nde/nde_checkpoint.h"
#include "core/nde/nde_compositing.h"
#include "core/nde/nde_composite.h"
#include "core/nde/nde_joint.h"
#include "core/nde/nde_replay.h"
#include "core/nde/nde_replay_internal.h"
#include "io/image_format_flis.h"
#include "yyjson.h"

static void add_reason(nde_chain *chain, const char *fmt, ...) G_GNUC_PRINTF(2, 3);
static void add_reason(nde_chain *chain, const char *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	g_ptr_array_add(chain->reasons, g_strdup_vprintf(fmt, ap));
	va_end(ap);
	/* Every reason is added while the build stands at the blocked position:
	 * records->len is the index the offending member would have (or has not
	 * yet) taken, so everything before it is still clean. */
	if (chain->records->len < chain->first_block)
		chain->first_block = chain->records->len;
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
gchar *nde_tier_c_invalid_reason(const nde_record *rec) {
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

/* First reason @rec cannot be replayed (heap string, caller frees), or NULL
 * when it replays.  The chain build decides whether an invalid member is a
 * hard blocker (no output checkpoint) or a barrier restart point. */
static gchar *member_invalid_reason(const nde_record *rec) {
	if (rec->tier == NDE_TIER_C)
		return nde_tier_c_invalid_reason(rec);
	/* A composite node's params are read by nde_composite_state_parse rather
	 * than by a descriptor, and chain membership already required that to
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
	/* Several outputs, and no way to produce them.  An ordinary op's hook
	 * returns ONE image, so replaying a record that wrote two would reproduce
	 * one of them and claim to have reproduced the step.  Joint records are
	 * exempt because theirs is a solved case: the hook applies one
	 * participant's share per replay of that participant's chain, which is why
	 * a joint record can be a member of several chains at once.  Refuse the
	 * rest until the apply path can carry N results (nde-simplify-plan S5.7). */
	if (nde_record_output_count(rec) > 1 && !nde_joint_is_op(rec->op_id))
		return g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) produces more than one image — replaying that is not supported yet"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
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
	destroy_any_args(user);
	return NULL;
}

/* POLICY predicates — see nde_replay.h.  Liveness and trial-chain
 * replayability are the execute path's job, not these. */
gboolean nde_record_amendable(const nde_record *rec) {
	if (!rec || rec->tier != NDE_TIER_A)
		return FALSE;
	const nde_op_class *cls = nde_op_class_for(rec->op_id);
	switch (cls->family) {
	case NDE_OPC_COMPOSITING:
		/* Editable params, no descriptor: the registry validates them
		 * (nde_op_class_params_valid), and nothing else gates the edit. */
		return TRUE;
	case NDE_OPC_COMPOSITE:
		/* The same, but only while the node is replayable — amending the
		 * opacity of a merge nobody can re-run would change the log and not
		 * the image. */
		return nde_composite_record_replayable(rec);
	default:
		/* An op's params are editable exactly when they can be read back. */
		return cls->desc && cls->desc->deserialize;
	}
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
nde_chain *nde_chain_build_excluding(gint item_id, gint64 exclude_record_id) {
	nde_chain *chain = g_new0(nde_chain, 1);
	chain->item_id = item_id;
	chain->records = g_ptr_array_new();
	chain->reasons = g_ptr_array_new_with_free_func(g_free);
	chain->snapshot = nde_history_snapshot(NULL);
	chain->member_flags = g_array_new(FALSE, TRUE, sizeof(guint8));
	chain->first_block = G_MAXUINT;
	gboolean is_flis = is_current_image_flis();
	/* Barrier tracking (phase 4): tail_possible follows the LAST freeze
	 * cause in document order — TRUE when it left a restart point. */
	gboolean tail_possible = TRUE;

	for (guint i = 0; chain->snapshot && i < chain->snapshot->len; i++) {
		nde_record *rec = g_ptr_array_index(chain->snapshot, i);
		gboolean member = FALSE;
		if (exclude_record_id && rec->record_id == exclude_record_id) {
			/* A DELETED geometry step still has to be undone on the canvas,
			 * not merely skipped: the layer is currently sitting where that
			 * step put it, and dropping it must return the layer to the
			 * position the surviving members produce.  So the chain keeps
			 * carrying geometry even though the step itself is gone —
			 * without this the pixels revert and the offset does not, and
			 * the layer ends up correct-but-misplaced. */
			gboolean geometric = nde_op_class_for(rec->op_id)->traits &
			                     NDE_OPT_GEOMETRIC;
			gboolean mine = nde_record_writes_item(rec, item_id);
			if (geometric && mine && is_flis &&
			    nde_checkpoint_baseline_has_position(item_id))
				chain->has_geometry = TRUE;
			continue;
		}
		/* Compositing-state records are inputs to the compositor, not pixel
		 * operations: neither chain members nor blockers (nde_compositing.h). */
		if (nde_compositing_is_op(rec->op_id))
			continue;
		switch (rec->scope) {
		case NDE_SCOPE_LAYER:
			/* Membership is "this record writes that item".  Usually that
			 * is its target and nothing else; a JOINT record (nde_joint.h)
			 * targets its anchor participant but writes every one of them,
			 * so it is a member of each of their chains — one record at one
			 * log position, shared. */
			member = nde_record_writes_item(rec, item_id);
			/* Registration is a joint op that MOVES its participants as well
			 * as warping them.  It is scope LAYER (its target is an anchor
			 * participant, not the canvas), so the NDE_SCOPE_CANVAS branch
			 * below never sees it — but the replay still has to thread the
			 * layer position through it, and that needs a recorded starting
			 * position exactly as an ordinary geometry step does. */
			if (member && nde_joint_is_geometric_op(rec->op_id) && is_flis) {
				if (nde_checkpoint_baseline_has_position(item_id)) {
					chain->has_geometry = TRUE;
				} else {
					member = FALSE;
					add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) moves this layer on the canvas, but no starting position was recorded — not replayable"),
					           rec->record_id, rec->op_id ? rec->op_id : "?");
				}
			}
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
				if (nde_checkpoint_baseline_has_position(item_id)) {
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
			const nde_op_class *cls = nde_op_class_for(rec->op_id);
			gboolean destructive = cls->traits & NDE_OPT_DESTRUCTIVE;
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
			} else if (!(cls->traits & NDE_OPT_CHAIN_IGNORE) && !destructive) {
				/* FAIL CLOSED: a DOCUMENT-scope record that is not known to
				 * be ignorable mutated pixels document-wide (icc.convert via
				 * the layer worker today; unknown ops from newer builds
				 * tomorrow — NDE_OPC_UNKNOWN carries no traits precisely so
				 * that it lands here).  Every layer's chain spans it, so no
				 * layer chain with records on both sides can replay without
				 * it.
				 *
				 * !destructive is not redundant: a merge or flatten targeting
				 * ANOTHER item falls through to here, and it consumed someone
				 * else's pixels, not ours. */
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
					if (chain->records->len < chain->first_block)
						chain->first_block = chain->records->len;
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
		/* An item whose first record is the composite that produced it was
		 * born of a merge or a flatten.  Like a mask, it has no baseline of
		 * its own — its origin is the composite's inputs. */
		chain->from_composite = nde_composite_is_op(first->op_id);
	}
	if (chain->is_mask) {
		const nde_record *first = g_ptr_array_index(chain->records, 0);
		if (!nde_record_input(first, "image") && !record_names_a_file(first)) {
			add_reason(chain, _("record %" G_GINT64_FORMAT " (%s) does not say what this mask was built from — not replayable"),
			           first->record_id, first->op_id ? first->op_id : "?");
			chain->first_block = 0;   /* the origin itself: no prefix survives */
		}
	} else if (chain->records->len && !chain->from_composite &&
	           !nde_checkpoint_baseline_exists(item_id)) {
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
	return nde_chain_build_excluding(item_id, 0);
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
