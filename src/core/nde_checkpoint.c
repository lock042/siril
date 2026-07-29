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
#include "core/nde_retention.h"
#include "io/image_format_flis.h"

/* ======================================================================= */
/* Tables (LEAF lock; steal-then-unref for releases)                       */
/* ======================================================================= */

static GMutex      cp_mutex;
static GHashTable *cp_table;    /* item_id (GINT_TO_POINTER) -> nde_snap* (ref held) */
static GHashTable *out_table;   /* gint64* (owned key) -> nde_snap* (ref held) */

/* Layer offsets, kept BESIDE the snapshots rather than inside them so every
 * existing call site keeps its signature (nde_checkpoint.h).  Same keys as
 * the tables above; an entry is present only when an offset was recorded. */
typedef struct { gint x, y; } cp_offset;
static GHashTable *cp_pos;      /* item_id  -> cp_offset* (owned) */
static GHashTable *out_pos;     /* gint64*  -> cp_offset* (owned) */

static void cp_tables_ensure_locked(void) {
	if (!cp_table)
		cp_table = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		                                 NULL, NULL /* values unref'd manually */);
	if (!out_table)
		out_table = g_hash_table_new_full(g_int64_hash, g_int64_equal,
		                                  g_free, NULL /* values unref'd manually */);
	if (!cp_pos)
		cp_pos = g_hash_table_new_full(g_direct_hash, g_direct_equal, NULL, g_free);
	if (!out_pos)
		out_pos = g_hash_table_new_full(g_int64_hash, g_int64_equal, g_free, g_free);
}

/* ---- layer offsets ----------------------------------------------------- */

void nde_checkpoint_baseline_set_offset(gint item_id, gint pos_x, gint pos_y) {
	cp_offset *o = g_new(cp_offset, 1);
	o->x = pos_x;
	o->y = pos_y;
	g_mutex_lock(&cp_mutex);
	cp_tables_ensure_locked();
	/* First writer wins, exactly as the baseline itself does: the baseline is
	 * the pre-FIRST-op state, so the position that goes with it is the one
	 * the layer had then — not wherever a later capture found it.  And only
	 * meaningful alongside stored pixels: an offset for a baseline we do not
	 * have would outlive what it describes. */
	if (g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id)) &&
	    !g_hash_table_contains(cp_pos, GINT_TO_POINTER(item_id)))
		g_hash_table_insert(cp_pos, GINT_TO_POINTER(item_id), o);
	else
		g_clear_pointer(&o, g_free);
	g_mutex_unlock(&cp_mutex);
}

gboolean nde_checkpoint_baseline_get_offset(gint item_id, gint *pos_x, gint *pos_y) {
	g_mutex_lock(&cp_mutex);
	const cp_offset *o = cp_pos ? g_hash_table_lookup(cp_pos, GINT_TO_POINTER(item_id)) : NULL;
	if (o) {
		if (pos_x) *pos_x = o->x;
		if (pos_y) *pos_y = o->y;
	}
	g_mutex_unlock(&cp_mutex);
	return o != NULL;
}

void nde_checkpoint_output_set_offset(gint64 record_id, gint pos_x, gint pos_y) {
	cp_offset *o = g_new(cp_offset, 1);
	o->x = pos_x;
	o->y = pos_y;
	gint64 *key = g_new(gint64, 1);
	*key = record_id;
	g_mutex_lock(&cp_mutex);
	cp_tables_ensure_locked();
	if (g_hash_table_contains(out_table, &record_id)) {
		g_hash_table_insert(out_pos, key, o);
		key = NULL;
		o = NULL;
	}
	g_mutex_unlock(&cp_mutex);
	g_free(key);
	g_free(o);
}

gboolean nde_checkpoint_output_get_offset(gint64 record_id, gint *pos_x, gint *pos_y) {
	g_mutex_lock(&cp_mutex);
	const cp_offset *o = out_pos ? g_hash_table_lookup(out_pos, &record_id) : NULL;
	if (o) {
		if (pos_x) *pos_x = o->x;
		if (pos_y) *pos_y = o->y;
	}
	g_mutex_unlock(&cp_mutex);
	return o != NULL;
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

	nde_retention_enforce(0);   /* outside the lock; may drop OLDER pins */
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
			if (cp_pos) {
				gpointer o = g_hash_table_lookup(cp_pos, GINT_TO_POINTER(from_item));
				if (o) {
					g_hash_table_steal(cp_pos, GINT_TO_POINTER(from_item));
					g_hash_table_insert(cp_pos, GINT_TO_POINTER(to_item), o);
				}
			}
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

	/* Protect what was just stored: over a budget this small, keeping the
	 * newest checkpoint beats storing it and dropping it again. */
	nde_retention_enforce(record_id);
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

/* ---- pin coordinates (see the header) ----------------------------------- */

void nde_checkpoint_store_at(const fits *f, gint item_id, gint64 record_id) {
	if (record_id)
		nde_checkpoint_output_store(f, record_id, item_id);
	else
		nde_checkpoint_baseline_ensure(f, item_id);
}

fits *nde_checkpoint_get_at(gint item_id, gint64 record_id) {
	return record_id ? nde_checkpoint_output_get(record_id)
	                 : nde_checkpoint_baseline_get(item_id);
}

gboolean nde_checkpoint_exists_at(gint item_id, gint64 record_id) {
	return record_id ? nde_checkpoint_output_exists(record_id)
	                 : nde_checkpoint_baseline_exists(item_id);
}

void nde_checkpoint_output_drop(gint64 record_id) {
	nde_snap *old = NULL;
	gpointer old_key = NULL;
	g_mutex_lock(&cp_mutex);
	if (out_table)
		g_hash_table_steal_extended(out_table, &record_id, &old_key, (gpointer *)&old);
	if (out_pos)
		g_hash_table_remove(out_pos, &record_id);   /* the offset describes those pixels */
	g_mutex_unlock(&cp_mutex);
	g_free(old_key);
	nde_snap_unref(old);
}

/* ======================================================================= */
/* Retention (design note §7 — driven by nde_retention.c)                  */
/* ======================================================================= */

gint64 nde_checkpoint_bytes(void) {
	gint64 total = 0;
	g_mutex_lock(&cp_mutex);
	GHashTableIter it;
	gpointer k, v;
	if (cp_table) {
		g_hash_table_iter_init(&it, cp_table);
		while (g_hash_table_iter_next(&it, &k, &v))
			total += nde_snap_size(v);
	}
	if (out_table) {
		g_hash_table_iter_init(&it, out_table);
		while (g_hash_table_iter_next(&it, &k, &v))
			total += nde_snap_size(v);
	}
	g_mutex_unlock(&cp_mutex);
	return total;
}

gboolean nde_checkpoint_evict_oldest(nde_checkpoint_eviction *out, gint64 protect_record) {
	gint protect = nde_checkpoint_active_item_id();
	nde_snap *victim = NULL;
	gpointer victim_key = NULL;
	gboolean is_baseline = FALSE;
	gint64  victim_record = 0;
	gint    victim_item = 0;

	g_mutex_lock(&cp_mutex);
	GHashTableIter it;
	gpointer k, v;
	/* Output checkpoints first, lowest record_id (= oldest) first.  Losing one
	 * costs the steps that depended on it their editability; losing a baseline
	 * costs an item its WHOLE chain, so baselines only go once nothing else is
	 * left. */
	if (out_table) {
		gint64 best = 0;
		g_hash_table_iter_init(&it, out_table);
		while (g_hash_table_iter_next(&it, &k, &v)) {
			gint64 rid = *(gint64 *)k;
			if (rid == protect_record)
				continue;
			if (!victim || rid < best) { best = rid; victim = v; victim_key = k; }
		}
		if (victim) {
			victim_record = best;
			nde_snap_tag_get(victim, &victim_item, NULL, NULL);
			g_hash_table_steal(out_table, victim_key);
			if (out_pos)
				g_hash_table_remove(out_pos, &victim_record);
		}
	}
	if (!victim && cp_table) {
		gint best = 0;
		g_hash_table_iter_init(&it, cp_table);
		while (g_hash_table_iter_next(&it, &k, &v)) {
			gint item = GPOINTER_TO_INT(k);
			/* Never the item being worked on: evicting its baseline would make
			 * the very next operation unreplayable, which is not a degradation
			 * the user could act on. */
			if (item == protect)
				continue;
			if (!victim || item < best) { best = item; victim = v; victim_key = k; }
		}
		if (victim) {
			is_baseline = TRUE;
			victim_item = best;
			g_hash_table_steal(cp_table, victim_key);
			g_hash_table_remove(cp_pos, victim_key);
		}
	}
	g_mutex_unlock(&cp_mutex);

	if (!victim)
		return FALSE;
	if (out) {
		out->bytes       = nde_snap_size(victim);
		out->is_baseline = is_baseline;
		out->item_id     = victim_item;
		out->record_id   = victim_record;
	}
	if (!is_baseline)
		g_free(victim_key);        /* the stolen gint64 key */
	nde_snap_unref(victim);
	return TRUE;
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
				if (out_pos)
					g_hash_table_remove(out_pos, k);
				g_hash_table_iter_steal(&it);
			}
		}
	}
	if (cp_pos)
		g_hash_table_remove(cp_pos, GINT_TO_POINTER(item_id));
	g_mutex_unlock(&cp_mutex);
	for (guint i = 0; i < doomed->len; i++)
		nde_snap_unref(g_ptr_array_index(doomed, i));
	g_ptr_array_unref(doomed);
	g_ptr_array_unref(doomed_keys);   /* frees the stolen gint64 keys */
}

void nde_checkpoint_purge(void) {
	/* A new document explains itself once; it must not stay silent because a
	 * previous one already hit the limit. */
	nde_retention_notice_reset();
	GHashTable *t1 = NULL, *t2 = NULL;
	g_mutex_lock(&cp_mutex);
	t1 = cp_table;
	cp_table = NULL;
	t2 = out_table;
	out_table = NULL;
	GHashTable *p1 = cp_pos, *p2 = out_pos;
	cp_pos = NULL;
	out_pos = NULL;
	g_mutex_unlock(&cp_mutex);
	if (p1) g_hash_table_destroy(p1);
	if (p2) g_hash_table_destroy(p2);
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
