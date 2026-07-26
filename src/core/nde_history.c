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
#include "core/nde_snapstore.h"
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
	g_free(rec);
}

/* ======================================================================= */
/* Canonical log                                                           */
/* ======================================================================= */

static nde_history *nde_history_new(void) {
	nde_history *h = g_new0(nde_history, 1);
	h->records = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
	h->next_record_id = 1;
	return h;
}

void nde_history_free(nde_history *h) {
	if (!h)
		return;
	g_ptr_array_unref(h->records);
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
	truncate_dead_locked(h, dropped);
	rec->record_id = h->next_record_id++;
	g_ptr_array_add(h->records, rec);
	h->live_count = h->records->len;
	gint64 id = rec->record_id;
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

void nde_history_on_undo(gint64 record_id) {
	if (!record_id || !com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h) {
		gint idx = find_index_locked(h, record_id);
		if (idx < 0)
			siril_log_debug("nde: undo of unknown record %" G_GINT64_FORMAT "\n", record_id);
		else if (h->live_count > (guint)idx)
			h->live_count = (guint)idx;
	}
	g_mutex_unlock(&nde_mutex);
	nde_history_notify_panel();
}

void nde_history_on_redo(gint64 record_id) {
	if (!record_id || !com.uniq)
		return;
	g_mutex_lock(&nde_mutex);
	nde_history *h = com.uniq->nde_history;
	if (h) {
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

guint nde_history_live_count(void) {
	if (!com.uniq)
		return 0;
	g_mutex_lock(&nde_mutex);
	guint n = com.uniq->nde_history ? com.uniq->nde_history->live_count : 0;
	g_mutex_unlock(&nde_mutex);
	return n;
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

gboolean nde_history_amend(gint64 record_id, const gchar *new_params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!com.uniq || !com.uniq->nde_history) {
		*err = g_strdup(_("no edit history"));
		return FALSE;
	}

	/* Validate the new params OUTSIDE the leaf mutex (deserializers may
	 * allocate freely).  The op id is read under the mutex first. */
	gchar *op_id = NULL;
	int op_version = 0;
	g_mutex_lock(&nde_mutex);
	{
		gint idx = find_mutable_locked(com.uniq->nde_history, record_id, TRUE, err);
		if (idx >= 0) {
			nde_record *rec = g_ptr_array_index(com.uniq->nde_history->records, idx);
			op_id = g_strdup(rec->op_id);
			op_version = rec->op_version;
		}
	}
	g_mutex_unlock(&nde_mutex);
	if (!op_id)
		return FALSE;

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
	/* destructor-first convention */
	void (*destroy)(void *) = *(void (**)(void *))trial;
	if (destroy)
		destroy(trial);
	else
		free(trial);
	g_free(op_id);

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
		return FALSE;
	}
	nde_history_notify_panel();
	return TRUE;
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

gint64 nde_capture_structural(const char *op_id, gint scope,
                              gint target_item_id, gchar *params,
                              const char *summary) {
	nde_record *rec = nde_record_new();
	rec->op_id          = g_strdup(op_id);
	rec->op_version     = 1;
	rec->scope          = scope;
	rec->target_item_id = target_item_id;
	rec->tier           = NDE_TIER_A;
	rec->params         = params;   /* ownership transferred */
	/* Structural records are not chain members — never a barrier, no
	 * checkpoint (post == NULL). */
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
	if (is_current_image_flis()) {
		flis_layer_t *lay = flis_active_layer();
		rec->target_item_id = lay ? lay->item_id : -1;
	}
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
