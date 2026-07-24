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
#ifndef _NDE_HISTORY_H_
#define _NDE_HISTORY_H_

/**
 * \file nde_history.h
 * \brief Nondestructive-editing operation log (provenance phase).
 *
 * The current document (com.uniq) carries an ordered log of the operations
 * that produced its pixels.  Phase 1 is provenance-only: records are captured
 * from generic_image_worker and structural FLIS commands, persisted in .flis
 * files (FLIS_HIST table), and displayed — nothing replays yet.
 *
 * Threading: the canonical log is guarded by an internal GMutex that is a
 * strict LEAF in the lock order — no other lock is acquired and no blocking
 * call is made while it is held (see flis-nde-sketch.md §6a).  All public
 * functions may be called from any thread.  Callers prepare strings BEFORE
 * calling so the critical sections stay a few pointer operations long.
 *
 * Records ≤ live_count reflect the current pixel state; records beyond it
 * were undone (Ctrl-Z) and are dropped ("truncated dead") when a new
 * operation is appended, mirroring the redo-stack discard semantics.
 * Saving persists only the live prefix so a file's history always matches
 * its pixels.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
	NDE_SCOPE_LAYER    = 0,  /* affects the target layer's pixels only */
	NDE_SCOPE_CANVAS   = 1,  /* changes canvas/layer geometry (crop, resample...) */
	NDE_SCOPE_DOCUMENT = 2,  /* structural: add/remove/reorder layer, flatten... */
} nde_scope;

typedef enum {
	NDE_TIER_A = 0,  /* replayable: op has serialize/deserialize */
	NDE_TIER_B = 1,  /* opaque: recorded for provenance, not replayable */
} nde_tier;

typedef struct nde_record {
	gint64    record_id;      /* monotonic per document, never reused; 0 = unassigned */
	gchar    *op_id;          /* owned copy — survives unknown ops on reload */
	gint      op_version;
	gint      target_item_id; /* layer item_id; -1 for canvas/document scope or plain image */
	gint      scope;          /* nde_scope */
	gchar    *params;         /* serialized key=value blob (NULL for Tier B) */
	gchar    *summary;        /* human label (log_hook SUMMARY / description) */
	gchar    *timestamp;      /* ISO 8601 UTC, matching layer CREATED/MODIFIED style */
	gchar    *impl;           /* producing implementation, FLISIMPL-style string */
	gint      tier;           /* nde_tier */
	gboolean  mask_active;    /* a processing/layer mask was in effect at capture */
	gchar    *mask_ref;       /* reserved for mask snapshots; always NULL in phase 1 */
} nde_record;

typedef struct nde_history {
	GPtrArray *records;      /* of nde_record*, position == order */
	guint      live_count;   /* records[0..live_count-1] reflect current pixels */
	gint64     next_record_id;
	gboolean   stale;        /* load-time PIXHASH mismatch — see §14.4 */
} nde_history;

/* ---- record lifecycle -------------------------------------------------- */

nde_record *nde_record_new(void);
nde_record *nde_record_copy(const nde_record *rec);
void        nde_record_free(nde_record *rec);

/* ---- canonical log (current document, com.uniq) ------------------------ */

/**
 * Append @rec to the current document's log, creating the log lazily.
 * Takes ownership of @rec (also on failure).  Any dead tail beyond
 * live_count is dropped first.  Assigns and returns the record id;
 * returns 0 (record freed, nothing recorded) when no single image is
 * loaded (com.uniq == NULL).
 */
gint64 nde_history_append(nde_record *rec);

/**
 * Undo/redo coherence (sketch §13.3): called by undo_display_data() after a
 * successful restore of an undo entry tagged with @record_id.  on_undo moves
 * live_count back to just before that record; on_redo forward to include it.
 * Unknown ids are ignored with a debug log (the record may have been
 * truncated by an intervening append — the tagged entry dies with the redo
 * stack in that case, so this is defensive only).
 */
void nde_history_on_undo(gint64 record_id);
void nde_history_on_redo(gint64 record_id);

/**
 * Deep-copy the live prefix for saving or display: returns a GPtrArray of
 * nde_record* with a free-func set (g_ptr_array_unref releases everything).
 * NULL when there is no history or it is empty.  When @next_id_out is
 * non-NULL it receives next_record_id (persisted as FLISHSEQ).
 * No I/O and no other lock is taken while the internal mutex is held.
 */
GPtrArray *nde_history_snapshot(gint64 *next_id_out);

/**
 * Like nde_history_snapshot() but INCLUDING the dead tail (records undone
 * but not yet truncated), for display: the History panel dims records
 * beyond the live count.  When @live_count_out is non-NULL it receives the
 * live count.  The save path must keep using nde_history_snapshot() — only
 * the live prefix is ever persisted.  NULL when there are no records.
 */
GPtrArray *nde_history_snapshot_all(guint *live_count_out);

/** Live-record count of the current document's log (0 if none). */
guint nde_history_live_count(void);

/**
 * Install @h as the current document's history, freeing any existing one.
 * Used by load_flis(); @h may be NULL (pure reset).  Ownership transfers.
 */
void nde_history_attach(nde_history *h);

/** Free a detached history object (also used by free_image_data via attach(NULL)). */
void nde_history_free(nde_history *h);

/* ---- amend-and-replay log mutations (phase 3, nde-phase2-3-plan.md P3.A) --
 * These are the LOG-side commits: the pixel side must already have been
 * replayed and swapped in successfully (nde_amend_execute / P3.B owns the
 * ordering).  Both operate on LIVE Tier-A records only.  Callers must own
 * the processing slot so capture/undo cannot interleave.                    */

/**
 * Replace @record_id's params with @new_params (validated by a deserialize
 * round-trip against the record's op), refresh its timestamp and impl, keep
 * its record_id, and truncate the dead tail (it described a pixel lineage
 * that no longer exists).  FALSE + heap message in @err on unknown/dead
 * record, Tier-B record, unknown op or unparsable params.
 */
gboolean nde_history_amend(gint64 record_id, const gchar *new_params, gchar **err);

/**
 * Remove LIVE record @record_id from the log and truncate the dead tail.
 * Unlike amend this accepts Tier-B records — deleting an opaque step is
 * well-defined (its removal never requires replaying it) and is how a
 * chain regains editability around a one-off opaque operation.  Policy
 * checks beyond liveness (structural / compositing-state records) live in
 * nde_delete_execute, which is the only sanctioned caller path.
 */
gboolean nde_history_delete(gint64 record_id, gchar **err);

/** Stale flag accessors (load-time PIXHASH mismatch). */
void     nde_history_set_stale(gboolean stale);
gboolean nde_history_is_stale(void);

/**
 * Schedule a panel refresh.  Phase-1 placeholder until the History panel
 * lands (work-order step 8): currently a no-op kept so capture sites are
 * already wired.
 */
void nde_history_notify_panel(void);

/* ---- capture-site helpers ---------------------------------------------- */

/**
 * Build+append a structural (Tier A) record in one call and return the new
 * record id (0 = nothing recorded, e.g. no document loaded).  Structural ops
 * carry their params inline (a handful of ints, sketch §13.2) — no serializer
 * machinery; phase-2 replays them from op_id + params.  Takes ownership of
 * @params (may be NULL); @summary and @op_id are copied.  timestamp/impl are
 * filled here.  Callers MUST invoke this only AFTER the mutation succeeded,
 * never on a failure path (records must reflect real pixel/document changes).
 */
gint64 nde_capture_structural(const char *op_id, gint scope,
                              gint target_item_id, gchar *params,
                              const char *summary);

/**
 * Build+append an opaque (Tier B, params == NULL) record and return the new
 * record id.  For mutations recorded for provenance only, with no replayable
 * parameters (python pixel/mask writes — sketch §13.2).  @summary and @op_id
 * are copied.  Same success-only capture discipline as nde_capture_structural.
 */
gint64 nde_capture_opaque(const char *op_id, gint scope,
                          gint target_item_id, const char *summary);

/**
 * Build+append a record for a descriptor-identified operation applied
 * OUTSIDE generic_image_worker — the commit-the-preview dialogs and other
 * GUI paths that manage their own undo (and therefore run the worker with
 * skip_generic_undo, or not at all).  Mirrors the worker's capture block:
 * Tier A with op->serialize(@params) when the descriptor has a serializer,
 * Tier B otherwise; target/scope/mask state resolved from the live
 * document.  Same success-only discipline: call it only once the pixels
 * have actually changed, right before saving the operation's undo entry,
 * and couple with undo_tag_top_nde_record() when that save succeeds.
 */
struct op_descriptor;
gint64 nde_capture_from_descriptor(const struct op_descriptor *op,
                                   gconstpointer params, const char *summary);

/** Heap ISO 8601 UTC timestamp, matching FLIS layer CREATED/MODIFIED style. */
gchar *nde_iso8601_now(void);

/** Heap producing-implementation string, matching the FLISIMPL keyword. */
gchar *nde_impl_string(void);

/* ---- key=value params codec (sketch §12) -------------------------------- *
 * Semicolon-delimited key=value pairs, same family as the layer METADATA
 * column.  Values are backslash-escaped for ';', '=', '\' and newline; keys
 * are restricted to [A-Za-z0-9_] (enforced by assertion — keys are
 * programmer-supplied constants).  Floats use %.9g (round-trips float32),
 * doubles %.17g.  The parser keeps every key it finds, known or not.       */

GString *nde_kv_start(void);
void     nde_kv_add_str(GString *kv, const char *key, const char *val);
void     nde_kv_add_int(GString *kv, const char *key, gint64 val);
void     nde_kv_add_bool(GString *kv, const char *key, gboolean val);
void     nde_kv_add_float(GString *kv, const char *key, float val);
void     nde_kv_add_double(GString *kv, const char *key, double val);
/** Finish building: frees @kv and returns the heap blob (caller g_free()s). */
gchar   *nde_kv_end(GString *kv);

/**
 * Parse a params blob into a table of key → unescaped value (both owned by
 * the table; free with g_hash_table_unref).  Never NULL; empty table for
 * NULL/empty blobs.  Malformed pairs (no '=') are skipped with a debug log.
 */
GHashTable *nde_kv_parse(const char *blob);

/** Lookup helpers; the typed variants return FALSE when absent/unparsable. */
const char *nde_kv_get_str(GHashTable *kv, const char *key);
gboolean    nde_kv_get_int(GHashTable *kv, const char *key, gint64 *out);
gboolean    nde_kv_get_bool(GHashTable *kv, const char *key, gboolean *out);
gboolean    nde_kv_get_float(GHashTable *kv, const char *key, float *out);
gboolean    nde_kv_get_double(GHashTable *kv, const char *key, double *out);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_HISTORY_H_ */
