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
 * NDE checkpoint tables — baselines per layer item and barrier output
 * checkpoints per record.  Since convergence phase C1
 * (nde-convergence-plan.md) the pixel payloads live in the shared
 * refcounted snapshot store (nde_snapstore.[ch]); this module only OWNS
 * references and keys them.  The public API is unchanged from phase 2.
 *
 * Locking: cp_mutex guards the two tables and is a strict LEAF — but
 * nde_snap_unref() may take the snapstore's own leaf mutex, so releases
 * follow the STEAL-THEN-UNREF pattern: detach the reference under
 * cp_mutex, unref after unlocking.  Snapshot creation/reading (swap I/O)
 * always runs outside cp_mutex.  Table lookups take a reference with the
 * atomic nde_snap_ref() (no snapstore lock) before unlocking.
 *
 * Tagging: baselines register in the snapstore tag registry as
 * POST(record 0) of their item; barrier outputs as POST(record_id).  The
 * registry is how C3's cached-restart resolution finds them alongside
 * undo-owned and pool snapshots — one lookup path for every source.
 */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/nde_snapstore.h"
#include "core/nde_checkpoint.h"
#include "io/image_format_flis.h"

/* ======================================================================= */
/* Tables (LEAF lock; steal-then-unref for releases)                       */
/* ======================================================================= */

static GMutex      cp_mutex;
static GHashTable *cp_table;    /* item_id (GINT_TO_POINTER) -> nde_snap* (ref held) */
static GHashTable *out_table;   /* gint64* (owned key) -> nde_snap* (ref held) */

static void cp_tables_ensure_locked(void) {
	if (!cp_table)
		cp_table = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		                                 NULL, NULL /* values unref'd manually */);
	if (!out_table)
		out_table = g_hash_table_new_full(g_int64_hash, g_int64_equal,
		                                  g_free, NULL /* values unref'd manually */);
}

/* ======================================================================= */
/* Baselines                                                               */
/* ======================================================================= */

void nde_checkpoint_baseline_ensure(const fits *pre, gint item_id) {
	if (!pre)
		return;
	/* Fast path: bail without any copy if a baseline already exists. */
	g_mutex_lock(&cp_mutex);
	gboolean exists = cp_table &&
	    g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id));
	g_mutex_unlock(&cp_mutex);
	if (exists)
		return;

	/* Prepare the (expensive) snapshot OUTSIDE the lock. */
	nde_snap *s = nde_snap_create(pre);
	if (!s)
		return;
	nde_snap_set_tag(s, item_id, 0, TRUE);   /* baseline = POST(0) */

	/* Re-check under the lock: first writer keeps its snapshot — the
	 * baseline must reflect the pre-FIRST-op pixels. */
	g_mutex_lock(&cp_mutex);
	cp_tables_ensure_locked();
	if (g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id))) {
		g_mutex_unlock(&cp_mutex);
		nde_snap_unref(s);
		return;
	}
	g_hash_table_insert(cp_table, GINT_TO_POINTER(item_id), s);
	g_mutex_unlock(&cp_mutex);
}

/* Rebind every checkpoint tagged for @from_item to @to_item — the plain →
 * FLIS promote moves the whole flat-image history onto the new base layer.
 * Retagging takes the snapstore's own leaf lock, so gather the snaps under
 * cp_mutex and retag outside it (same discipline as everywhere else in this
 * file).  Promote runs synchronously on the main thread with no jobs
 * active, so the gathered pointers cannot be dropped concurrently. */
void nde_checkpoint_rebind_item(gint from_item, gint to_item) {
	GPtrArray *retag = g_ptr_array_new();
	g_mutex_lock(&cp_mutex);
	if (cp_table) {
		nde_snap *s = g_hash_table_lookup(cp_table, GINT_TO_POINTER(from_item));
		if (s && !g_hash_table_contains(cp_table, GINT_TO_POINTER(to_item))) {
			g_hash_table_steal(cp_table, GINT_TO_POINTER(from_item));
			g_hash_table_insert(cp_table, GINT_TO_POINTER(to_item), s);
			g_ptr_array_add(retag, s);
		}
	}
	if (out_table) {
		GHashTableIter it;
		gpointer k, v;
		g_hash_table_iter_init(&it, out_table);
		while (g_hash_table_iter_next(&it, &k, &v))
			g_ptr_array_add(retag, v);
	}
	g_mutex_unlock(&cp_mutex);

	for (guint i = 0; i < retag->len; i++) {
		nde_snap *s = g_ptr_array_index(retag, i);
		gint item = 0;
		gint64 rid = 0;
		gboolean post = FALSE;
		if (nde_snap_tag_get(s, &item, &rid, &post) && item == from_item)
			nde_snap_set_tag(s, to_item, rid, post);
	}
	g_ptr_array_unref(retag);
}

void nde_checkpoint_baseline_adopt(const fits *src, gint item_id) {
	if (!src)
		return;
	nde_snap *s = nde_snap_create(src);   /* outside the lock */
	if (!s)
		return;
	nde_snap_set_tag(s, item_id, 0, TRUE);
	nde_snap *old = NULL;
	g_mutex_lock(&cp_mutex);
	cp_tables_ensure_locked();
	old = g_hash_table_lookup(cp_table, GINT_TO_POINTER(item_id));
	if (old)
		g_hash_table_steal(cp_table, GINT_TO_POINTER(item_id));
	g_hash_table_insert(cp_table, GINT_TO_POINTER(item_id), s);
	g_mutex_unlock(&cp_mutex);
	nde_snap_unref(old);   /* on-disk baseline is authoritative on load */
}

fits *nde_checkpoint_baseline_get(gint item_id) {
	g_mutex_lock(&cp_mutex);
	nde_snap *s = cp_table ?
	    g_hash_table_lookup(cp_table, GINT_TO_POINTER(item_id)) : NULL;
	if (s)
		nde_snap_ref(s);   /* atomic — no snapstore lock */
	g_mutex_unlock(&cp_mutex);
	if (!s)
		return NULL;
	fits *out = nde_snap_read(s);   /* I/O outside the lock */
	nde_snap_unref(s);
	return out;
}

gboolean nde_checkpoint_baseline_exists(gint item_id) {
	g_mutex_lock(&cp_mutex);
	gboolean e = cp_table &&
	    g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id));
	g_mutex_unlock(&cp_mutex);
	return e;
}

/* ======================================================================= */
/* Output checkpoints (phase-4 barriers)                                   */
/* ======================================================================= */

static void output_insert(const fits *src, gint64 record_id, gint item_id) {
	if (!src || record_id <= 0)
		return;
	nde_snap *s = nde_snap_create(src);   /* outside the lock */
	if (!s)
		return;
	nde_snap_set_tag(s, item_id, record_id, TRUE);
	gint64 *key = g_new(gint64, 1);
	*key = record_id;
	nde_snap *old = NULL;
	gpointer old_key = NULL;
	g_mutex_lock(&cp_mutex);
	cp_tables_ensure_locked();
	g_hash_table_steal_extended(out_table, &record_id, &old_key, (gpointer *)&old);
	g_hash_table_insert(out_table, key, s);
	g_mutex_unlock(&cp_mutex);
	g_free(old_key);
	nde_snap_unref(old);   /* NULL-safe */
}

void nde_checkpoint_output_store(const fits *post, gint64 record_id, gint item_id) {
	output_insert(post, record_id, item_id);
}

void nde_checkpoint_output_adopt(const fits *src, gint64 record_id, gint item_id) {
	output_insert(src, record_id, item_id);
}

fits *nde_checkpoint_output_get(gint64 record_id) {
	g_mutex_lock(&cp_mutex);
	nde_snap *s = out_table ? g_hash_table_lookup(out_table, &record_id) : NULL;
	if (s)
		nde_snap_ref(s);
	g_mutex_unlock(&cp_mutex);
	if (!s)
		return NULL;
	fits *out = nde_snap_read(s);
	nde_snap_unref(s);
	return out;
}

gboolean nde_checkpoint_output_exists(gint64 record_id) {
	g_mutex_lock(&cp_mutex);
	gboolean e = out_table && g_hash_table_contains(out_table, &record_id);
	g_mutex_unlock(&cp_mutex);
	return e;
}

void nde_checkpoint_output_drop(gint64 record_id) {
	nde_snap *old = NULL;
	gpointer old_key = NULL;
	g_mutex_lock(&cp_mutex);
	if (out_table)
		g_hash_table_steal_extended(out_table, &record_id, &old_key, (gpointer *)&old);
	g_mutex_unlock(&cp_mutex);
	g_free(old_key);
	nde_snap_unref(old);
}

/* ======================================================================= */
/* Drop / purge                                                            */
/* ======================================================================= */

void nde_checkpoint_drop(gint item_id) {
	GPtrArray *doomed = g_ptr_array_new();
	GPtrArray *doomed_keys = g_ptr_array_new_with_free_func(g_free);
	g_mutex_lock(&cp_mutex);
	if (cp_table) {
		nde_snap *s = g_hash_table_lookup(cp_table, GINT_TO_POINTER(item_id));
		if (s) {
			g_hash_table_steal(cp_table, GINT_TO_POINTER(item_id));
			g_ptr_array_add(doomed, s);
		}
	}
	if (out_table) {
		/* Output checkpoints of the dying layer die with it.  The owning
		 * item is carried in the snapstore tag, set before insertion and
		 * immutable afterwards — safe to read here. */
		GHashTableIter it;
		gpointer k, v;
		g_hash_table_iter_init(&it, out_table);
		while (g_hash_table_iter_next(&it, &k, &v)) {
			nde_snap *s = v;
			gint tag_item;
			if (nde_snap_tag_get(s, &tag_item, NULL, NULL) && tag_item == item_id) {
				g_ptr_array_add(doomed, s);
				g_ptr_array_add(doomed_keys, k);
				g_hash_table_iter_steal(&it);
			}
		}
	}
	g_mutex_unlock(&cp_mutex);
	for (guint i = 0; i < doomed->len; i++)
		nde_snap_unref(g_ptr_array_index(doomed, i));
	g_ptr_array_unref(doomed);
	g_ptr_array_unref(doomed_keys);   /* frees the stolen gint64 keys */
}

void nde_checkpoint_purge(void) {
	GHashTable *t1 = NULL, *t2 = NULL;
	g_mutex_lock(&cp_mutex);
	t1 = cp_table;
	cp_table = NULL;
	t2 = out_table;
	out_table = NULL;
	g_mutex_unlock(&cp_mutex);
	if (t1) {
		GHashTableIter it;
		gpointer k, v;
		g_hash_table_iter_init(&it, t1);
		while (g_hash_table_iter_next(&it, &k, &v))
			nde_snap_unref(v);
		g_hash_table_destroy(t1);
	}
	if (t2) {
		GHashTableIter it;
		gpointer k, v;
		g_hash_table_iter_init(&it, t2);
		while (g_hash_table_iter_next(&it, &k, &v))
			nde_snap_unref(v);
		g_hash_table_destroy(t2);
	}
}

gint nde_checkpoint_active_item_id(void) {
	if (is_current_image_flis()) {
		flis_layer_t *lay = flis_active_layer();
		if (lay)
			return lay->item_id;
	}
	return -1;
}
