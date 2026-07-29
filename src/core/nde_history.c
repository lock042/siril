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
 * Nondestructive-editing operation log — see nde_history.h for the design
 * contract and flis-nde-sketch.md §10–20 for the phase-1 plan.
 */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/op_descriptor.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_compositing.h"
#include "core/nde_composite.h"
#include "core/nde_snapstore.h"
#include "core/nde_replay.h"    /* nde_item_is_retained_input() */
#include "io/image_format_flis.h"

/* Leaf lock guarding the CURRENT document's log (com.uniq->nde_history and
 * the fields inside it).  Nothing else is ever acquired while this is held —
 * keep it that way (sketch §6a). */
static GMutex nde_mutex;

/* ======================================================================= */
/* Record lifecycle                                                        */
/* ======================================================================= */

nde_record *nde_record_new(void) {
	nde_record *rec = g_new0(nde_record, 1);
	rec->target_item_id = -1;
	return rec;
}

nde_input_pin *nde_input_pin_new(const char *role, gint src_item_id,
                                 gint64 src_record_id) {
	nde_input_pin *pin = g_new0(nde_input_pin, 1);
	pin->role = g_strdup(role);
	pin->src_item_id = src_item_id;
	pin->src_record_id = src_record_id;
	return pin;
}

void nde_input_pin_free(nde_input_pin *pin) {
	if (!pin)
		return;
	g_free(pin->role);
	g_free(pin);
}

static GPtrArray *pins_new(void) {
	return g_ptr_array_new_with_free_func((GDestroyNotify)nde_input_pin_free);
}

void nde_record_add_input(nde_record *rec, const char *role,
                          gint src_item_id, gint64 src_record_id) {
	g_return_if_fail(rec != NULL && role != NULL);
	if (!rec->inputs)
		rec->inputs = pins_new();
	for (guint i = 0; i < rec->inputs->len; i++) {
		nde_input_pin *p = g_ptr_array_index(rec->inputs, i);
		if (!g_strcmp0(p->role, role)) {
			p->src_item_id = src_item_id;
			p->src_record_id = src_record_id;
			return;
		}
	}
	g_ptr_array_add(rec->inputs, nde_input_pin_new(role, src_item_id, src_record_id));
}

const nde_input_pin *nde_record_input(const nde_record *rec, const char *role) {
	if (!rec || !rec->inputs || !role)
		return NULL;
	for (guint i = 0; i < rec->inputs->len; i++) {
		const nde_input_pin *p = g_ptr_array_index(rec->inputs, i);
		if (!g_strcmp0(p->role, role))
			return p;
	}
	return NULL;
}

const nde_input_pin *nde_record_input_by_item(const nde_record *rec, gint item_id) {
	if (!rec || !rec->inputs)
		return NULL;
	for (guint i = 0; i < rec->inputs->len; i++) {
		const nde_input_pin *p = g_ptr_array_index(rec->inputs, i);
		if (p->src_item_id == item_id)
			return p;
	}
	return NULL;
}

nde_record *nde_record_copy(const nde_record *rec) {
	if (!rec)
		return NULL;
	nde_record *copy = g_new0(nde_record, 1);
	*copy = *rec;
	copy->op_id     = g_strdup(rec->op_id);
	copy->params    = g_strdup(rec->params);
	copy->summary   = g_strdup(rec->summary);
	copy->timestamp = g_strdup(rec->timestamp);
	copy->impl      = g_strdup(rec->impl);
	copy->mask_ref  = g_strdup(rec->mask_ref);
	copy->inputs    = NULL;
	for (guint i = 0; rec->inputs && i < rec->inputs->len; i++) {
		const nde_input_pin *p = g_ptr_array_index(rec->inputs, i);
		nde_record_add_input(copy, p->role, p->src_item_id, p->src_record_id);
	}
	return copy;
}

void nde_record_free(nde_record *rec) {
	if (!rec)
		return;
	g_free(rec->op_id);
	g_free(rec->params);
	g_free(rec->summary);
	g_free(rec->timestamp);
	g_free(rec->impl);
	g_free(rec->mask_ref);
	if (rec->inputs)
		g_ptr_array_unref(rec->inputs);
	g_free(rec);
}

/* ---- pin list codec ---------------------------------------------------- */

gchar *nde_pins_serialize(GPtrArray *pins) {
	if (!pins || !pins->len)
		return NULL;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "n", pins->len);
	for (guint i = 0; i < pins->len; i++) {
		const nde_input_pin *p = g_ptr_array_index(pins, i);
		gchar key[32];
		g_snprintf(key, sizeof(key), "role%u", i);
		nde_kv_add_str(kv, key, p->role ? p->role : "");
		g_snprintf(key, sizeof(key), "item%u", i);
		nde_kv_add_int(kv, key, p->src_item_id);
		g_snprintf(key, sizeof(key), "rec%u", i);
		nde_kv_add_int(kv, key, p->src_record_id);
	}
	return nde_kv_end(kv);
}

GPtrArray *nde_pins_parse(const char *blob) {
	if (!blob || !*blob)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	gint64 n = 0;
	GPtrArray *out = NULL;
	if (nde_kv_get_int(kv, "n", &n) && n > 0) {
		for (gint64 i = 0; i < n; i++) {
			gchar key[32];
			g_snprintf(key, sizeof(key), "role%" G_GINT64_FORMAT, i);
			const char *role = nde_kv_get_str(kv, key);
			g_snprintf(key, sizeof(key), "item%" G_GINT64_FORMAT, i);
			gint64 item = 0;
			gboolean have_item = nde_kv_get_int(kv, key, &item);
			g_snprintf(key, sizeof(key), "rec%" G_GINT64_FORMAT, i);
			gint64 src = 0;
			nde_kv_get_int(kv, key, &src);   /* absent = baseline */
			/* Skip an incomplete pin rather than inventing one: a pin that
			 * names no role or no source item points nowhere. */
			if (!role || !*role || !have_item)
				continue;
			if (!out)
				out = pins_new();
			g_ptr_array_add(out, nde_input_pin_new(role, (gint)item, src));
		}
	}
	g_hash_table_unref(kv);
	return out;
}

gchar *nde_pins_to_string(GPtrArray *pins) {
	if (!pins || !pins->len)
		return NULL;
	GString *s = g_string_new(NULL);
	for (guint i = 0; i < pins->len; i++) {
		const nde_input_pin *p = g_ptr_array_index(pins, i);
		if (i)
			g_string_append(s, ", ");
		if (p->src_record_id)
			g_string_append_printf(s, "%s@%d:%" G_GINT64_FORMAT,
			                       p->role ? p->role : "?", p->src_item_id,
			                       p->src_record_id);
		else
			g_string_append_printf(s, "%s@%d:baseline",
			                       p->role ? p->role : "?", p->src_item_id);
	}
	return g_string_free(s, FALSE);
}

/* ======================================================================= */
/* Canonical log                                                           */
/* ======================================================================= */

static nde_history *nde_history_new(void) {
	nde_history *h = g_new0(nde_history, 1);
	h->records = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
	h->next_record_id = 1;
	h->ins_ids = g_array_new(FALSE, FALSE, sizeof(gint64));
	h->ins_stash = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
	return h;
}

void nde_history_free(nde_history *h) {
	if (!h)
		return;
	g_ptr_array_unref(h->records);
	if (h->ins_ids)
		g_array_unref(h->ins_ids);
	if (h->ins_stash)
		g_ptr_array_unref(h->ins_stash);
	g_free(h);
}

/* Drop records beyond live_count.  Mutex held by caller.  The ids of the
 * dropped records are collected into @dropped (when non-NULL) so the caller
 * can release their output checkpoints AFTER unlocking — nde_checkpoint's
 * lock is a separate leaf and must never nest inside this one. */
static void truncate_dead_locked(nde_history *h, GArray *dropped) {
	while (h->records->len > h->live_count) {
		if (dropped) {
			nde_record *rec = g_ptr_array_index(h->records, h->records->len - 1);
			g_array_append_val(dropped, rec->record_id);
		}
		g_ptr_array_remove_index(h->records, h->records->len - 1);
	}
}

/* Release the output checkpoints AND cache-pool entries of dropped
 * records; consumes @ids. */
static void drop_output_checkpoints(GArray *ids) {
	if (!ids)
		return;
	for (guint i = 0; i < ids->len; i++) {
		gint64 id = g_array_index(ids, gint64, i);
		nde_checkpoint_output_drop(id);
		nde_snapstore_evict_record(id);
	}
	g_array_unref(ids);
}

static gint find_index_locked(nde_history *h, gint64 record_id);

/* Does @rec belong in the armed insertion point's item chain?  Anything that
 * nde_chain_build() would call a member of ins_item's chain qualifies, because
 * membership is exactly "this record's pixels are part of that lineage":
 * LAYER-scope records targeting the item, plus — on a plain image, where the
 * image IS the layer — CANVAS-scope geometry.  Mutex held. */
static gboolean insert_qualifies_locked(const nde_history *h, const nde_record *rec) {
	if (rec->scope == NDE_SCOPE_LAYER)
		return rec->target_item_id == h->ins_item;
	if (rec->scope == NDE_SCOPE_CANVAS)
		return h->ins_item < 0;   /* plain image */
	return FALSE;
}

/* Records the insertion cannot survive.  The end path restores the true
 * pixels before replaying forward, which silently undoes any pixel change the
 * insertion did not account for — so a record that changed the target's
 * pixels WITHOUT qualifying for insertion leaves the log describing something
 * the image no longer shows.  Callers refuse these operations up front
 * (nde_edit_at_refuses_op); this flag is the backstop for a path that does
 * not, and it makes the insertion abandon rather than lie.  Mutex held. */
static gboolean insert_disturbs(const nde_history *h, const nde_record *rec) {
	return (rec->scope == NDE_SCOPE_CANVAS && h->ins_item >= 0) ||
	       !g_strcmp0(rec->op_id, "document.flatten") ||
	       !g_strcmp0(rec->op_id, "layer.merge_down");
}

gint64 nde_history_append(nde_record *rec) {
	g_return_val_if_fail(rec != NULL, 0);
	if (!com.uniq) {
		/* No single image (sequence frame or nothing loaded): provenance
		 * has no document to attach to. */
		nde_record_free(rec);
		return 0;
	}
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	if (!com.uniq->nde_history)
		com.uniq->nde_history = nde_history_new();
	nde_history *h = com.uniq->nde_history;
	gint insert_at = -1;
	if (h->ins_before) {
		if (insert_qualifies_locked(h, rec)) {
			insert_at = find_index_locked(h, h->ins_before);
			if (insert_at < 0)
				/* The anchor vanished under us (unreachable: history edits
				 * are refused while the point is armed).  Fall back to an
				 * ordinary append and let the pixel side abandon. */
				h->ins_disturbed = TRUE;
		} else if (insert_disturbs(h, rec)) {
			h->ins_disturbed = TRUE;
		}
	}
	rec->record_id = h->next_record_id++;
	gint64 id = rec->record_id;
	if (insert_at >= 0) {
		/* A new step supersedes anything undone at this point, exactly as an
		 * ordinary append discards the dead tail. */
		for (guint i = 0; i < h->ins_stash->len; i++) {
			nde_record *dead = g_ptr_array_index(h->ins_stash, i);
			g_array_append_val(dropped, dead->record_id);
		}
		g_ptr_array_set_size(h->ins_stash, 0);
		g_ptr_array_insert(h->records, insert_at, rec);
		h->live_count++;
		g_array_append_val(h->ins_ids, id);
	} else {
		/* Ordinary append.  While the point is armed there is no dead tail
		 * to truncate (it went when the point was armed) — the call is a
		 * cheap no-op then. */
		truncate_dead_locked(h, dropped);
		g_ptr_array_add(h->records, rec);
		h->live_count = h->records->len;
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	return id;
}

/* Index of the record with @record_id, or -1.  Mutex held by caller. */
static gint find_index_locked(nde_history *h, gint64 record_id) {
	for (guint i = 0; i < h->records->len; i++) {
		nde_record *rec = g_ptr_array_index(h->records, i);
		if (rec->record_id == record_id)
			return (gint)i;
	}
	return -1;
}

/* Undo/redo while an insertion point is armed.  live_count still counts the
 * WHOLE log, so the ordinary "move live_count back to idx" rule would declare
 * every record after the insertion dead — the lineage after the anchor is not
 * being undone, one inserted step is.  Instead the inserted record is lifted
 * out of the log and stashed, and redo puts it back at the anchor.  Undo is
 * LIFO, so the stash is a stack and only the last inserted id can be undone.
 * Returns TRUE when the event was handled here.  Mutex held. */
static gboolean insert_on_undo_locked(nde_history *h, gint64 record_id) {
	if (!h->ins_before)
		return FALSE;
	guint n = h->ins_ids->len;
	if (n && g_array_index(h->ins_ids, gint64, n - 1) == record_id) {
		gint idx = find_index_locked(h, record_id);
		if (idx >= 0) {
			g_ptr_array_add(h->ins_stash, g_ptr_array_steal_index(h->records, idx));
			h->live_count--;
			g_array_set_size(h->ins_ids, n - 1);
		}
		return TRUE;
	}
	/* Not one of ours (or not the newest): the pre-anchor state was swapped
	 * in wholesale and undo was flushed when the point was armed, so nothing
	 * else on the stack describes this log.  Ignore rather than corrupt
	 * live_count. */
	siril_log_debug("nde: undo of record %" G_GINT64_FORMAT " ignored while inserting\n",
	                record_id);
	return TRUE;
}

static gboolean insert_on_redo_locked(nde_history *h, gint64 record_id) {
	if (!h->ins_before)
		return FALSE;
	guint n = h->ins_stash->len;
	nde_record *top = n ? g_ptr_array_index(h->ins_stash, n - 1) : NULL;
	if (top && top->record_id == record_id) {
		gint at = find_index_locked(h, h->ins_before);
		if (at >= 0) {
			g_ptr_array_insert(h->records, at, g_ptr_array_steal_index(h->ins_stash, n - 1));
			h->live_count++;
			g_array_append_val(h->ins_ids, record_id);
		}
		return TRUE;
	}
	siril_log_debug("nde: redo of record %" G_GINT64_FORMAT " ignored while inserting\n",
	                record_id);
	return TRUE;
}

void nde_history_on_undo(gint64 record_id) {
	if (!record_id || !com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h && !insert_on_undo_locked(h, record_id)) {
		gint idx = find_index_locked(h, record_id);
		if (idx < 0)
			siril_log_debug("nde: undo of unknown record %" G_GINT64_FORMAT "\n", record_id);
		else if (h->live_count > (guint)idx)
			h->live_count = (guint)idx;
	}
	g_mutex_unlock(&nde_mutex);
	nde_history_notify_panel();
}

/* Rebind every record targeting @from_item to @to_item.  Used by the plain
 * image → FLIS promote: the whole flat-image history (target -1) now
 * describes the new single base layer, so the chain stays one editable
 * lineage across the promote instead of splitting at the save point. */
void nde_history_rebind_item(gint from_item, gint to_item) {
	if (!com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h) {
		for (guint i = 0; i < h->records->len; i++) {
			nde_record *rec = g_ptr_array_index(h->records, i);
			if (rec->target_item_id == from_item)
				rec->target_item_id = to_item;
		}
	}
	g_mutex_unlock(&nde_mutex);
	nde_history_notify_panel();
}

void nde_history_on_redo(gint64 record_id) {
	if (!record_id || !com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h && !insert_on_redo_locked(h, record_id)) {
		gint idx = find_index_locked(h, record_id);
		if (idx < 0)
			siril_log_debug("nde: redo of unknown record %" G_GINT64_FORMAT "\n", record_id);
		else if (h->live_count < (guint)idx + 1)
			h->live_count = (guint)idx + 1;
	}
	g_mutex_unlock(&nde_mutex);
	nde_history_notify_panel();
}

GPtrArray *nde_history_snapshot(gint64 *next_id_out) {
	if (next_id_out)
		*next_id_out = 1;
	if (!com.uniq)
		return NULL;
	GPtrArray *out = NULL;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h && h->live_count > 0) {
		out = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
		for (guint i = 0; i < h->live_count; i++)
			g_ptr_array_add(out, nde_record_copy(g_ptr_array_index(h->records, i)));
	}
	if (h && next_id_out)
		*next_id_out = h->next_record_id;
	g_mutex_unlock(&nde_mutex);
	return out;
}

GPtrArray *nde_history_snapshot_all(guint *live_count_out) {
	if (live_count_out)
		*live_count_out = 0;
	if (!com.uniq)
		return NULL;
	GPtrArray *out = NULL;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h && h->records->len > 0) {
		out = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
		for (guint i = 0; i < h->records->len; i++)
			g_ptr_array_add(out, nde_record_copy(g_ptr_array_index(h->records, i)));
		if (live_count_out)
			*live_count_out = h->live_count;
	}
	g_mutex_unlock(&nde_mutex);
	return out;
}

gint64 nde_history_last_record_for_item(gint item_id) {
	if (!com.uniq)
		return 0;
	gint64 id = 0;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	for (guint i = 0; h && i < h->live_count; i++) {
		const nde_record *rec = g_ptr_array_index(h->records, i);
		if (rec->target_item_id == item_id)
			id = rec->record_id;   /* last one wins: records are in order */
	}
	g_mutex_unlock(&nde_mutex);
	return id;
}

guint nde_history_live_count(void) {
	if (!com.uniq)
		return 0;
	g_mutex_lock(&nde_mutex);
	guint n = com.uniq->nde_history ? com.uniq->nde_history->live_count : 0;
	g_mutex_unlock(&nde_mutex);
	return n;
}

gint nde_history_max_item_id(const nde_history *h) {
	gint max = 0;
	if (!h || !h->records)
		return 0;
	for (guint i = 0; i < h->records->len; i++) {
		const nde_record *r = g_ptr_array_index(h->records, i);
		if (r->target_item_id > max)
			max = r->target_item_id;
		for (guint p = 0; r->inputs && p < r->inputs->len; p++) {
			const nde_input_pin *pin = g_ptr_array_index(r->inputs, p);
			if (pin->src_item_id > max)
				max = pin->src_item_id;
		}
	}
	return max;
}

void nde_history_attach(nde_history *h) {
	/* The document's log is being replaced (load / close / new document),
	 * so its in-session baselines and cached intermediate states no longer
	 * describe anything — drop them.  The FLIS loader adopts fresh
	 * NDE_BASE baselines AFTER this call.  Purges take separate LEAF
	 * locks, outside the nde history mutex. */
	nde_checkpoint_purge();
	nde_snapstore_pool_purge();
	if (!com.uniq) {
		nde_history_free(h);
		return;
	}
	nde_history *old;
	g_mutex_lock(&nde_mutex);
	old = com.uniq->nde_history;
	com.uniq->nde_history = h;
	g_mutex_unlock(&nde_mutex);
	nde_history_free(old);
	nde_history_notify_panel();
}

/* Shared head of amend/delete: find the record and enforce the LIVE
 * contract (plus Tier A when @require_tier_a — amend needs editable
 * params; delete may remove opaque records, whose removal never requires
 * replaying them).  Mutex held by the caller.  Returns the index or -1
 * with *err set. */
static gint find_mutable_locked(nde_history *h, gint64 record_id,
                                gboolean require_tier_a, gchar **err) {
	gint idx = h ? find_index_locked(h, record_id) : -1;
	if (idx < 0) {
		*err = g_strdup_printf(_("no record with id %" G_GINT64_FORMAT), record_id);
		return -1;
	}
	if ((guint)idx >= h->live_count) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " has been undone — redo it first"), record_id);
		return -1;
	}
	nde_record *rec = g_ptr_array_index(h->records, idx);
	if (require_tier_a && rec->tier != NDE_TIER_A) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) is opaque and cannot be edited"),
		                       record_id, rec->op_id ? rec->op_id : "?");
		return -1;
	}
	return idx;
}

/* Log-side commit shared by both amend flavours (descriptor-validated ops and
 * compositing-state records): swap in the new params, refresh timestamp/impl,
 * optionally replace the summary, and truncate the dead tail.  @new_summary is
 * owned by this function (freed or transferred).  Params are validated by the
 * CALLER — by here the edit is known good. */
guint nde_history_drop_mask_input(gint64 record_id, gint mask_item_id) {
	if (!com.uniq || !com.uniq->nde_history || mask_item_id <= 0)
		return 0;
	/* The params transform runs OUTSIDE the mutex (§6a: no allocation-heavy
	 * work under the leaf lock), so read the blob, rewrite it, then swap. */
	gchar *old_params = NULL;
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		for (guint i = 0; i < h->live_count; i++) {
			nde_record *rec = g_ptr_array_index(h->records, i);
			if (rec->record_id == record_id) {
				old_params = g_strdup(rec->params);
				break;
			}
		}
	}
	g_mutex_unlock(&nde_mutex);
	if (!old_params)
		return 0;
	gchar *new_params = nde_composite_params_drop_mask(old_params, mask_item_id);
	g_free(old_params);
	if (!new_params)
		return 0;

	guint changed = 0;
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		for (guint i = 0; i < h->live_count; i++) {
			nde_record *rec = g_ptr_array_index(h->records, i);
			if (rec->record_id != record_id)
				continue;
			g_free(rec->params);
			rec->params = new_params;   /* ownership transferred */
			new_params = NULL;
			for (guint p = 0; rec->inputs && p < rec->inputs->len; p++) {
				const nde_input_pin *pin = g_ptr_array_index(rec->inputs, p);
				if (pin->src_item_id == mask_item_id) {
					g_ptr_array_remove_index(rec->inputs, p);
					break;
				}
			}
			changed = 1;
			break;
		}
	}
	g_mutex_unlock(&nde_mutex);
	g_free(new_params);   /* no-op when it was transferred */
	if (changed)
		nde_history_notify_panel();
	return changed;
}

static gboolean amend_commit(gint64 record_id, const gchar *new_params,
                             gchar *new_summary, gchar **err) {
	/* Strings prepared before locking (§6a discipline). */
	gchar *params_copy = g_strdup(new_params);
	gchar *ts = nde_iso8601_now();
	gchar *impl = nde_impl_string();

	gboolean ok = FALSE;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		gint idx = find_mutable_locked(h, record_id, TRUE, err);
		if (idx >= 0) {
			nde_record *rec = g_ptr_array_index(h->records, idx);
			g_free(rec->params);
			rec->params = params_copy;
			g_free(rec->timestamp);
			rec->timestamp = ts;
			g_free(rec->impl);
			rec->impl = impl;
			if (new_summary) {
				g_free(rec->summary);
				rec->summary = new_summary;   /* ownership transferred */
				new_summary = NULL;
			}
			truncate_dead_locked(h, dropped);
			ok = TRUE;
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	if (!ok) {
		g_free(params_copy);
		g_free(ts);
		g_free(impl);
		g_free(new_summary);
		return FALSE;
	}
	nde_history_notify_panel();
	return TRUE;
}

gboolean nde_history_amend(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!com.uniq || !com.uniq->nde_history) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}

	/* Validate the new params OUTSIDE the leaf mutex (deserializers may
	 * allocate freely).  The op id is read under the mutex first. */
	gchar *op_id = NULL, *old_params = NULL;
	int op_version = 0;
	g_mutex_lock(&nde_mutex);
	{
		gint idx = find_mutable_locked(com.uniq->nde_history, record_id, TRUE, err);
		if (idx >= 0) {
			nde_record *rec = g_ptr_array_index(com.uniq->nde_history->records, idx);
			op_id = g_strdup(rec->op_id);
			old_params = g_strdup(rec->params);
			op_version = rec->op_version;
		}
	}
	g_mutex_unlock(&nde_mutex);
	if (!op_id) {
		g_free(old_params);
		return FALSE;
	}

	/* A composite node has no op descriptor either, and its params say what it
	 * consumed as well as how: validation is against the RECORDED blob, so an
	 * amend can change the compositing state and nothing else (nde_composite.h). */
	if (nde_composite_is_op(op_id)) {
		gboolean ok = nde_composite_validate(old_params, new_params, err);
		g_free(op_id);
		g_free(old_params);
		return ok ? amend_commit(record_id, new_params, NULL, err) : FALSE;
	}
	g_free(old_params);

	/* Compositing-state records have no op descriptor: their params are
	 * validated by range/enum instead, and their summary ("Set opacity", …)
	 * is param-independent so it never goes stale.  See nde_compositing.h. */
	if (nde_compositing_is_op(op_id)) {
		gboolean ok = nde_compositing_validate(op_id, new_params, err);
		g_free(op_id);
		if (!ok)
			return FALSE;
		return amend_commit(record_id, new_params, NULL, err);
	}

	const op_descriptor *op = op_descriptor_by_id(op_id);
	if (!op || !op->deserialize) {
		*err = g_strdup_printf(_("operation '%s' cannot be edited by this build"), op_id);
		g_free(op_id);
		return FALSE;
	}
	gpointer trial = op->deserialize(new_params, op_version);
	if (!trial) {
		*err = g_strdup_printf(_("the new parameters for '%s' failed to parse"), op_id);
		g_free(op_id);
		return FALSE;
	}
	/* Regenerate the human summary from the NEW params while the deserialized
	 * struct is still alive, so the panel row and log stop showing the
	 * pre-amend values (params and pixels are already updated by the amend —
	 * only the label lagged).  log_hook(SUMMARY) formats the params struct and
	 * is param-pure for the amendable ops.  Ops without a log_hook captured a
	 * param-independent description that does not go stale, so leave theirs
	 * untouched (new_summary == NULL). */
	gchar *new_summary = op->log_hook ? op->log_hook(trial, SUMMARY) : NULL;
	/* destructor-first convention */
	void (*destroy)(void *) = *(void (**)(void *))trial;
	if (destroy)
		destroy(trial);
	else
		free(trial);
	g_free(op_id);

	return amend_commit(record_id, new_params, new_summary, err);
}

gboolean nde_history_delete(gint64 record_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!com.uniq || !com.uniq->nde_history) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}
	gboolean ok = FALSE;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		gint idx = find_mutable_locked(h, record_id, FALSE, err);
		if (idx >= 0) {
			g_array_append_val(dropped, record_id);  /* its own checkpoint */
			g_ptr_array_remove_index(h->records, idx);  /* free func frees it */
			h->live_count--;
			truncate_dead_locked(h, dropped);
			ok = TRUE;
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	if (ok)
		nde_history_notify_panel();
	return ok;
}

gboolean nde_history_reorder(gint64 record_id, gint64 before_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!com.uniq || !com.uniq->nde_history) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}
	gboolean ok = FALSE;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		gint idx = find_mutable_locked(h, record_id, FALSE, err);
		gint dest = -1;
		if (idx >= 0) {
			if (before_id == 0) {
				dest = (gint)h->live_count;
			} else {
				dest = find_index_locked(h, before_id);
				if (dest < 0 || (guint)dest >= h->live_count) {
					*err = g_strdup_printf(_("no live record with id %" G_GINT64_FORMAT " to move before"),
					                       before_id);
					dest = -1;
				} else if (before_id == record_id) {
					*err = g_strdup(_("cannot move a record before itself"));
					dest = -1;
				}
			}
		}
		if (idx >= 0 && dest >= 0) {
			nde_record *rec = g_ptr_array_steal_index(h->records, idx);
			if (dest > idx)
				dest--;   /* removal shifted the target left */
			g_ptr_array_insert(h->records, dest, rec);
			truncate_dead_locked(h, dropped);
			ok = TRUE;
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	if (ok)
		nde_history_notify_panel();
	return ok;
}

/* ---- edit-at insertion point (graph step 2) ---------------------------- */

gboolean nde_history_insert_point_set(gint64 before_id, gint item_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!com.uniq || !com.uniq->nde_history) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	gboolean ok = FALSE;
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		if (!h->ins_ids)
			h->ins_ids = g_array_new(FALSE, FALSE, sizeof(gint64));
		if (!h->ins_stash)
			h->ins_stash = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
		if (before_id == 0) {
			h->ins_before = 0;
			ok = TRUE;
		} else if (find_mutable_locked(h, before_id, FALSE, err) >= 0) {
			/* The dead tail goes now, not at the first insert: the mode
			 * flushes undo, so its redo lineage is unreachable from here
			 * on and leaving it would confuse the local liveness model. */
			truncate_dead_locked(h, dropped);
			h->ins_before = before_id;
			h->ins_item = item_id;
			h->ins_disturbed = FALSE;
			g_array_set_size(h->ins_ids, 0);
			g_ptr_array_set_size(h->ins_stash, 0);
			ok = TRUE;
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	return ok;
}

gint64 nde_history_insert_point(void) {
	if (!com.uniq)
		return 0;
	g_mutex_lock(&nde_mutex);
	gint64 id = com.uniq->nde_history ? com.uniq->nde_history->ins_before : 0;
	g_mutex_unlock(&nde_mutex);
	return id;
}

gboolean nde_history_insert_point_disturbed(void) {
	if (!com.uniq)
		return FALSE;
	g_mutex_lock(&nde_mutex);
	gboolean d = com.uniq->nde_history && com.uniq->nde_history->ins_disturbed;
	g_mutex_unlock(&nde_mutex);
	return d;
}

GArray *nde_history_insert_point_clear(void) {
	GArray *out = g_array_new(FALSE, FALSE, sizeof(gint64));
	if (!com.uniq)
		return out;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		if (h) {
			h->ins_before = 0;
			h->ins_disturbed = FALSE;
			if (h->ins_ids) {
				for (guint i = 0; i < h->ins_ids->len; i++) {
					gint64 id = g_array_index(h->ins_ids, gint64, i);
					g_array_append_val(out, id);
				}
				g_array_set_size(h->ins_ids, 0);
			}
			/* Stashed records were undone and never redone: they are gone
			 * for good, so release what they own. */
			if (h->ins_stash) {
				for (guint i = 0; i < h->ins_stash->len; i++) {
					nde_record *r = g_ptr_array_index(h->ins_stash, i);
					g_array_append_val(dropped, r->record_id);
				}
				g_ptr_array_set_size(h->ins_stash, 0);
			}
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	return out;
}

void nde_history_drop_records(GArray *ids) {
	if (!ids || !ids->len || !com.uniq)
		return;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	gboolean any = FALSE;
	g_mutex_lock(&nde_mutex);
	{
		nde_history *h = com.uniq->nde_history;
		for (guint i = 0; h && i < ids->len; i++) {
			gint64 id = g_array_index(ids, gint64, i);
			gint idx = find_index_locked(h, id);
			if (idx < 0 || (guint)idx >= h->live_count)
				continue;
			g_array_append_val(dropped, id);
			g_ptr_array_remove_index(h->records, idx);
			h->live_count--;
			any = TRUE;
		}
	}
	g_mutex_unlock(&nde_mutex);
	drop_output_checkpoints(dropped);
	if (any)
		nde_history_notify_panel();
}

void nde_history_truncate_dead(void) {
	if (!com.uniq)
		return;
	GArray *dropped = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_mutex_lock(&nde_mutex);
	if (com.uniq->nde_history)
		truncate_dead_locked(com.uniq->nde_history, dropped);
	g_mutex_unlock(&nde_mutex);
	gboolean any = dropped->len > 0;
	drop_output_checkpoints(dropped);
	if (any)
		nde_history_notify_panel();
}

void nde_history_set_stale(gboolean stale) {
	if (!com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	if (com.uniq->nde_history)
		com.uniq->nde_history->stale = stale;
	g_mutex_unlock(&nde_mutex);
}

gboolean nde_history_is_stale(void) {
	if (!com.uniq)
		return FALSE;
	g_mutex_lock(&nde_mutex);
	gboolean stale = com.uniq->nde_history && com.uniq->nde_history->stale;
	g_mutex_unlock(&nde_mutex);
	return stale;
}

void nde_history_notify_panel(void) {
	/* Never called with the leaf mutex held.  The GUI implementation
	 * schedules ONE coalesced async idle that snapshots the log and
	 * rebuilds the panel's main-thread mirror; headless is a stub. */
	gui_iface.nde_history_changed();
}

/* ======================================================================= */
/* Capture-site helpers                                                    */
/* ======================================================================= */

/* TRUE when @rec is a barrier: a live chain-member record that cannot be
 * replayed — Tier B, or Tier A with a mask active at capture (nde-phase4
 * plan "barrier model").  Barriers store their POST-op pixels as an output
 * checkpoint so everything after them stays editable. */
static gboolean record_is_barrier(const nde_record *rec) {
	return rec->tier == NDE_TIER_B ||
	       rec->tier == NDE_TIER_C ||   /* checkpointed: restart point for the tail */
	       (rec->tier == NDE_TIER_A && rec->mask_active);
}

/* Shared tail for the one-call capture helpers: fill timestamp/impl, append
 * (which takes ownership and assigns the id) and notify the panel.  @rec is
 * fully owned; on no-document the append frees it and returns 0.  When @post
 * is non-NULL and the record is a barrier, store its POST-op pixels as the
 * record's output checkpoint (nde-phase4 P4.3) — done AFTER append so the id
 * is assigned, and OUTSIDE the history mutex (append has already unlocked). */
static gint64 capture_finish(nde_record *rec, const char *summary,
                             const fits *post) {
	rec->summary   = g_strdup(summary);
	rec->timestamp = nde_iso8601_now();
	rec->impl      = nde_impl_string();
	gboolean barrier   = record_is_barrier(rec);
	gint     target_id = rec->target_item_id;
	gint64 id = nde_history_append(rec);   /* takes ownership */
	if (id > 0 && barrier && post)
		nde_checkpoint_output_store(post, id, target_id);
	nde_history_notify_panel();
	return id;
}

/* Whose history a capture belongs to.  There were two answers to this in two
 * files — the generic worker's and this one's — and they disagreed the moment
 * an insertion borrowed the display for an item a merge had consumed: the
 * worker's captures went to the borrowed item and the dialogs' captures went
 * to the active layer, so a stretch aimed at a merge input landed after the
 * merge.  One answer now, in one place. */
gint nde_capture_target_item(void) {
	gint borrowed = nde_edit_at_borrowed_item();
	if (borrowed)
		return borrowed;
	if (!is_current_image_flis())
		return -1;
	flis_layer_t *lay = flis_active_layer();
	return lay ? lay->item_id : -1;
}

gint64 nde_capture_structural(const char *op_id, gint scope,
                              gint target_item_id, gchar *params,
                              const char *summary) {
	return nde_capture_structural_pinned(op_id, scope, target_item_id, params,
	                                     summary, NULL, 0);
}

gint64 nde_capture_structural_pinned(const char *op_id, gint scope,
                                     gint target_item_id, gchar *params,
                                     const char *summary,
                                     const nde_pin_spec *pins, guint n_pins) {
	nde_record *rec = nde_record_new();
	rec->op_id          = g_strdup(op_id);
	rec->op_version     = 1;
	rec->scope          = scope;
	rec->target_item_id = target_item_id;
	rec->tier           = NDE_TIER_A;
	rec->params         = params;   /* ownership transferred */
	for (guint i = 0; i < n_pins; i++)
		nde_record_add_input(rec, pins[i].role, pins[i].src_item_id,
		                     pins[i].src_record_id);
	/* Structural records are not barriers and take no output checkpoint
	 * (post == NULL): an unpinned one is not a chain member at all, and a
	 * pinned composite re-runs from its inputs rather than from a stored
	 * copy of its result (nde_composite.h). */
	return capture_finish(rec, summary, NULL);
}

gint64 nde_capture_opaque(const char *op_id, gint scope,
                          gint target_item_id, const char *summary,
                          const fits *post) {
	nde_record *rec = nde_record_new();
	rec->op_id          = g_strdup(op_id);
	rec->op_version     = 1;
	rec->scope          = scope;
	rec->target_item_id = target_item_id;
	rec->tier           = NDE_TIER_B;
	rec->params         = NULL;
	return capture_finish(rec, summary, post);
}

gint64 nde_capture_script(const char *op_id, gint scope,
                          gint target_item_id, gchar *params,
                          const char *summary, const fits *post) {
	nde_record *rec = nde_record_new();
	rec->op_id          = g_strdup(op_id);
	rec->op_version     = 1;
	rec->scope          = scope;
	rec->target_item_id = target_item_id;
	rec->tier           = NDE_TIER_C;
	rec->params         = params;   /* ownership transferred (may be NULL) */
	return capture_finish(rec, summary, post);
}

gint64 nde_capture_from_descriptor(const op_descriptor *op,
                                   gconstpointer params, const char *summary,
                                   const fits *post) {
	g_return_val_if_fail(op != NULL, 0);
	nde_record *rec = nde_record_new();
	gboolean tier_a = op->serialize != NULL;
	rec->op_id      = g_strdup(op->id);
	rec->op_version = op->version;
	rec->tier       = tier_a ? NDE_TIER_A : NDE_TIER_B;
	rec->params     = tier_a ? op->serialize(params) : NULL;
	rec->scope      = (op->flags & OP_GEOMETRY_CHANGING) ?
	                  NDE_SCOPE_CANVAS : NDE_SCOPE_LAYER;
	rec->target_item_id = nde_capture_target_item();
	rec->mask_active = gfit && gfit->mask && gfit->mask_active;
	return capture_finish(rec, summary, post);
}

gchar *nde_iso8601_now(void) {
	GDateTime *now = g_date_time_new_now_utc();
	gchar *s = g_date_time_format_iso8601(now);
	g_date_time_unref(now);
	return s;
}

gchar *nde_impl_string(void) {
	return g_strdup_printf("Siril-%s-flis-gtk4", VERSION);
}

/* ======================================================================= */
/* key=value params codec (sketch §12)                                     */
/* ======================================================================= */

static gboolean key_is_valid(const char *key) {
	if (!key || !*key)
		return FALSE;
	for (const char *p = key; *p; p++) {
		if (!g_ascii_isalnum(*p) && *p != '_')
			return FALSE;
	}
	return TRUE;
}

/* Escape ';', '=', '\' and newline in @val, appending to @out. */
static void append_escaped(GString *out, const char *val) {
	for (const char *p = val; *p; p++) {
		switch (*p) {
		case ';':  g_string_append(out, "\\;");  break;
		case '=':  g_string_append(out, "\\=");  break;
		case '\\': g_string_append(out, "\\\\"); break;
		case '\n': g_string_append(out, "\\n");  break;
		default:   g_string_append_c(out, *p);   break;
		}
	}
}

GString *nde_kv_start(void) {
	return g_string_new(NULL);
}

void nde_kv_add_str(GString *kv, const char *key, const char *val) {
	g_return_if_fail(kv != NULL);
	g_return_if_fail(key_is_valid(key));
	if (kv->len)
		g_string_append_c(kv, ';');
	g_string_append(kv, key);
	g_string_append_c(kv, '=');
	if (val)
		append_escaped(kv, val);
}

void nde_kv_add_int(GString *kv, const char *key, gint64 val) {
	gchar buf[32];
	g_snprintf(buf, sizeof(buf), "%" G_GINT64_FORMAT, val);
	nde_kv_add_str(kv, key, buf);
}

void nde_kv_add_bool(GString *kv, const char *key, gboolean val) {
	nde_kv_add_str(kv, key, val ? "1" : "0");
}

void nde_kv_add_float(GString *kv, const char *key, float val) {
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
	g_ascii_formatd(buf, sizeof(buf), "%.9g", (double)val);
	nde_kv_add_str(kv, key, buf);
}

void nde_kv_add_double(GString *kv, const char *key, double val) {
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
	g_ascii_formatd(buf, sizeof(buf), "%.17g", val);
	nde_kv_add_str(kv, key, buf);
}

gchar *nde_kv_end(GString *kv) {
	g_return_val_if_fail(kv != NULL, NULL);
	return g_string_free(kv, FALSE);
}

GHashTable *nde_kv_parse(const char *blob) {
	GHashTable *kv = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
	if (!blob || !*blob)
		return kv;
	GString *key = g_string_new(NULL);
	GString *val = g_string_new(NULL);
	GString *cur = key;   /* filling the key until an unescaped '=' */
	const char *p = blob;
	for (;;) {
		char c = *p;
		if (c == '\\' && p[1]) {
			/* Unescape: "\n" is newline, anything else is the literal
			 * next character ("\;", "\=", "\\"). */
			g_string_append_c(cur, p[1] == 'n' ? '\n' : p[1]);
			p += 2;
			continue;
		}
		if (c == '=' && cur == key) {
			cur = val;
			p++;
			continue;
		}
		if (c == ';' || c == '\0') {
			if (cur == val)
				g_hash_table_replace(kv, g_strdup(key->str), g_strdup(val->str));
			else if (key->len)
				siril_log_debug("nde: skipping malformed params pair '%s'\n", key->str);
			if (c == '\0')
				break;
			g_string_truncate(key, 0);
			g_string_truncate(val, 0);
			cur = key;
			p++;
			continue;
		}
		g_string_append_c(cur, c);
		p++;
	}
	g_string_free(key, TRUE);
	g_string_free(val, TRUE);
	return kv;
}

const char *nde_kv_get_str(GHashTable *kv, const char *key) {
	return kv ? g_hash_table_lookup(kv, key) : NULL;
}

gboolean nde_kv_get_int(GHashTable *kv, const char *key, gint64 *out) {
	const char *s = nde_kv_get_str(kv, key);
	if (!s || !*s)
		return FALSE;
	gchar *end = NULL;
	gint64 v = g_ascii_strtoll(s, &end, 10);
	if (!end || *end)
		return FALSE;
	*out = v;
	return TRUE;
}

gboolean nde_kv_get_bool(GHashTable *kv, const char *key, gboolean *out) {
	gint64 v;
	if (!nde_kv_get_int(kv, key, &v))
		return FALSE;
	*out = (v != 0);
	return TRUE;
}

gboolean nde_kv_get_double(GHashTable *kv, const char *key, double *out) {
	const char *s = nde_kv_get_str(kv, key);
	if (!s || !*s)
		return FALSE;
	gchar *end = NULL;
	double v = g_ascii_strtod(s, &end);
	if (!end || *end)
		return FALSE;
	*out = v;
	return TRUE;
}

gboolean nde_kv_get_float(GHashTable *kv, const char *key, float *out) {
	double v;
	if (!nde_kv_get_double(kv, key, &v))
		return FALSE;
	*out = (float)v;
	return TRUE;
}

/* ======================================================================= */
/* External-input operands (phase 4.5 Convention 1)                        */
/* ======================================================================= */

gchar *nde_file_sha256(const char *path, gint64 *size_out) {
	if (!path || !*path)
		return NULL;
	FILE *f = g_fopen(path, "rb");
	if (!f)
		return NULL;
	GChecksum *ck = g_checksum_new(G_CHECKSUM_SHA256);
	if (!ck) {
		fclose(f);
		return NULL;
	}
	guchar buf[64 * 1024];
	gint64 total = 0;
	size_t n;
	gboolean io_error = FALSE;
	while ((n = fread(buf, 1, sizeof buf, f)) > 0) {
		g_checksum_update(ck, buf, n);
		total += (gint64)n;
	}
	if (ferror(f))
		io_error = TRUE;
	fclose(f);
	if (io_error) {
		g_checksum_free(ck);
		return NULL;
	}
	gchar *hex = g_strdup(g_checksum_get_string(ck));
	g_checksum_free(ck);
	if (size_out)
		*size_out = total;
	return hex;
}

gboolean nde_operand_verify(const char *path, gint64 expect_size,
                            const char *expect_sha256) {
	if (!path || !*path) {
		siril_log_error(_("operand file missing or changed: (no path recorded)\n"));
		return FALSE;
	}
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	if (!sha) {
		siril_log_error(_("operand file missing or changed: %s\n"), path);
		return FALSE;
	}
	gboolean ok = TRUE;
	if (size != expect_size)
		ok = FALSE;
	if (expect_sha256 && *expect_sha256 && g_ascii_strcasecmp(sha, expect_sha256) != 0)
		ok = FALSE;
	g_free(sha);
	if (!ok)
		siril_log_error(_("operand file missing or changed: %s\n"), path);
	return ok;
}
