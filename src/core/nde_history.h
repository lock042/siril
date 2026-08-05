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

/* The image type is `fits` (typedef of `struct ffit` in siril.h); forward-
 * declare the tag so this header stays light (includers pull in siril.h). */
struct ffit;
typedef struct ffit fits;

/* Item-id sentinels for the two things that are not FLIS items but still need
 * naming in the log.  Real ids are positive and come from next_item_id; masks
 * owned by a layer get theirs from flis_layer_lmask_id/pmask_id. */
#define NDE_ITEM_IMAGE      (-1)  /* the plain single image (no layer stack) */
#define NDE_ITEM_PLAIN_MASK (-2)  /* that image's processing mask.  Rebound onto
                                   * the base layer's mask id by promote, the
                                   * same way NDE_ITEM_IMAGE is rebound onto the
                                   * base layer. */

typedef enum {
	NDE_SCOPE_LAYER    = 0,  /* affects the target layer's pixels only */
	NDE_SCOPE_CANVAS   = 1,  /* changes canvas/layer geometry (crop, resample...) */
	NDE_SCOPE_DOCUMENT = 2,  /* structural: add/remove/reorder layer, flatten... */
} nde_scope;

typedef enum {
	NDE_TIER_A = 0,  /* replayable: op has serialize/deserialize */
	NDE_TIER_B = 1,  /* opaque: recorded for provenance, not replayable */
	NDE_TIER_C = 2,  /* replayable Python script (phase 5): re-run the script
	                  * with recorded args (params carry script/sha256/args).
	                  * Always checkpointed — until the replay half lands it
	                  * degrades safely to Tier-B-barrier semantics. */
} nde_tier;

/* An INPUT EDGE of a record (graph step 4).  A record has always named its
 * one output — target_item_id — and nothing else; that implicit model says
 * "this record depends on the previous record of the same item, and on
 * nothing else", and every gap in the history is a case where that is false.
 * A pin makes one of those dependencies explicit.
 *
 * The coordinate is the snapstore's: item @src_item_id as it stood just after
 * record @src_record_id.  Record ids are stable across reorder and never
 * reused, so a pin survives everything except deleting what it points at.
 *
 * The "image" input every pixel op consumes is implicit and never pinned —
 * it is the chain itself.  Roles name the EXTRA inputs: "mask" today. */
typedef struct nde_input_pin {
	gchar  *role;           /* owned; "mask", later "overlay", "base"… */
	gint    src_item_id;    /* the item supplying this input */
	gint64  src_record_id;  /* its state after this record; 0 = its baseline */
} nde_input_pin;

nde_input_pin *nde_input_pin_new(const char *role, gint src_item_id, gint64 src_record_id);
void           nde_input_pin_free(nde_input_pin *pin);

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
	gchar    *mask_ref;       /* dead: the pin list superseded it before it was
	                           * ever populated.  Still written (empty) to
	                           * FLIS_HIST so older builds find every column
	                           * they look up by name. */
	GPtrArray *inputs;        /* of nde_input_pin*, or NULL when the record
	                           * consumes only its chain's implicit image */
} nde_record;

/* Attach an input edge.  @role and the ids are copied; the record owns the
 * pin.  Adding the same role twice replaces the earlier pin — a record has at
 * most one input per role. */
void nde_record_add_input(nde_record *rec, const char *role,
                          gint src_item_id, gint64 src_record_id);

/** The pin for @role, or NULL.  Borrowed — valid while @rec is. */
const nde_input_pin *nde_record_input(const nde_record *rec, const char *role);

/** The pin naming @item_id as its source, or NULL.  A composite matches its
 *  pins to its recorded inputs this way: roles are labels for the graph view,
 *  and an item appears at most once among a record's inputs. */
const nde_input_pin *nde_record_input_by_item(const nde_record *rec, gint item_id);

/**
 * Stop @record_id from consuming the layer mask @mask_item_id: clears the flag
 * in its params and drops the pin.  Returns the number of records changed (0
 * or 1) — nothing happens if that record does not use that mask.
 *
 * Deliberately not reachable through nde_history_amend: this is not a change a
 * user makes to a step, it is the consequence of deleting the step that BUILT
 * the mask.  Leaving the pin behind would leave a phantom edge to an item with
 * no history and a record that says it used a mask nobody can produce.
 */
guint nde_history_drop_mask_input(gint64 record_id, gint mask_item_id);

/* ---- pin list codec ------------------------------------------------------
 * Encoded with the ordinary kv codec so the escaping rules are shared:
 * "n=2;role0=mask;item0=3;rec0=7;role1=base;item1=1;rec1=0".              */

/** Serialize @pins; NULL for a NULL/empty list (nothing to persist). */
gchar *nde_pins_serialize(GPtrArray *pins);

/** Parse a blob back into a GPtrArray of nde_input_pin* with a free func, or
 *  NULL when the blob is empty or names no complete pin. */
GPtrArray *nde_pins_parse(const char *blob);

/** Human-readable one-liner for the History and flis_history:
 *  "mask@3:7, base@1:baseline".  NULL when there are no pins. */
gchar *nde_pins_to_string(GPtrArray *pins);

/**
 * Record id of the last LIVE record targeting @item_id, or 0 when the item
 * has no history yet.  This is the "as it stood now" end of a fresh pin.
 */
gint64 nde_history_last_record_for_item(gint item_id);

/**
 * The coordinates of every stored mask state some record still reads: each
 * mask-role pin ("mask", composite "mask0", "mask1"…) in the log.  Fills
 * @record_ids (heap gint64 keys, table must own+free them) with the
 * src_record_id != 0 coordinates, and @baseline_items (direct gint keys) with
 * the src_record_id == 0 sources.  Either table may be NULL.  These are the
 * states retention must not evict: a mask survives being cleared only through
 * them, and losing one turns every step that used it into a permanent
 * barrier.
 */
void nde_history_mask_pin_coords(GHashTable *record_ids, GHashTable *baseline_items);

/** TRUE when any record in the log (live or undone) pins @item_id as an
 *  input.  "Has this mask ever been used?" is exactly this question. */
gboolean nde_history_item_is_pinned(gint item_id);

/**
 * Remove every record targeting @item_id from the log outright — live and
 * dead — and release their output checkpoints and pool states.  Returns the
 * number of records removed.
 *
 * NOT a history edit: there is no replay and no validation, so this is only
 * sound for an item whose records influenced no other pixels — a mask that
 * was cleared without ever having masked anything.  A used mask must keep its
 * records (and its pinned states) so its consumers stay replayable.
 */
guint nde_history_forget_item(gint item_id);

typedef struct nde_history {
	GPtrArray *records;      /* of nde_record*, position == order */
	guint      live_count;   /* records[0..live_count-1] reflect current pixels */
	gint64     next_record_id;
	gboolean   stale;        /* load-time PIXHASH mismatch — see §14.4 */

	/* Edit-at insertion point (graph step 2).  While ins_before is non-zero
	 * new qualifying records are INSERTED before that record instead of
	 * appended; see nde_history_insert_point_set(). */
	gint64     ins_before;    /* anchor record id; 0 = ordinary append mode */
	gint       ins_item;      /* item whose chain the insertion belongs to */
	gboolean   ins_disturbed; /* something ran that invalidates the forward replay */
	GArray    *ins_ids;       /* gint64: inserted record ids, in log order */
	GPtrArray *ins_stash;     /* undone inserted records, LIFO (owns them) */
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

/**
 * The highest item id @h mentions, as a target or through a pin — including
 * items the document no longer has: a merged-away layer, or an input a
 * composite consumed.  Seeding next_item_id from the LIVE layers alone lets a
 * consumed id be handed out again after a save and reload, and the new layer
 * then inherits a dead one's provenance.  Works on a detached history, since
 * load_flis() needs it before nde_history_attach().
 */
gint nde_history_max_item_id(const nde_history *h);

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

/**
 * Move LIVE record @record_id to sit immediately before LIVE record
 * @before_id (@before_id == 0 ⇒ move to the end of the live prefix), then
 * truncate the dead tail.  Positions are global log positions; pixel-
 * safety (chain membership, freeze rules) is nde_reorder_execute's job —
 * this is the log-side commit only.  Record ids and timestamps are
 * unchanged (ids are therefore no longer monotonic in log order after a
 * reorder — position, not id, is the order).
 */
gboolean nde_history_reorder(gint64 record_id, gint64 before_id, gchar **err);

/* ---- edit-at insertion point (graph step 2) -----------------------------
 * "Edit at K" installs the state just before record K into the image and
 * lets the user run ordinary operations there.  The log side of that is an
 * armed insertion point: while it is armed nde_history_append() inserts
 * before K rather than appending, so the new steps land in the middle of
 * the lineage.  nde_replay.c owns the pixel side (nde_edit_at_*).
 *
 * A record QUALIFIES for insertion when it would be a member of the item's
 * replay chain: LAYER-scope targeting the item, plus CANVAS-scope geometry
 * on a plain image (where the image is the layer).  Records for other
 * layers, and non-destructive structural records, are appended normally —
 * they are not part of this lineage and their position relative to K is
 * immaterial.
 *
 * A third class exists because the pixel side restores the true image before
 * replaying forward: a record that changed the target's pixels but did not
 * qualify would have that change silently reverted, leaving the log
 * describing something the image no longer shows.  Those (canvas geometry on
 * a LAYERED document, document.flatten, layer.merge_down) raise the
 * "disturbed" flag and the insertion is abandoned instead.  They are refused
 * before they run wherever possible (nde_edit_at_refuses_op) — the flag is
 * the backstop, not the mechanism.
 *
 * The dead tail is truncated when the point is armed, not when a record is
 * inserted: entering the mode flushes undo, so a redo lineage beyond the
 * live prefix is unreachable from that moment on.  While armed, undo and
 * redo of the INSERTED records are tracked locally (ins_stash) so live_count
 * keeps meaning "the whole log", and undo/redo of anything else is ignored.
 */

/** Arm the point before live record @before_id for item @item_id.  FALSE +
 *  heap @err when @before_id is not a live record. */
gboolean nde_history_insert_point_set(gint64 before_id, gint item_id, gchar **err);

/** Anchor record id, or 0 when no insertion point is armed. */
gint64   nde_history_insert_point(void);

/** TRUE while an insertion point is armed; *item_out gets the item whose
 *  chain it belongs to (may be -1: the plain image).  The Layers panel uses
 *  it to refuse removing the very layer an insertion is armed on. */
gboolean nde_history_insert_point_target(gint *item_out);

/** TRUE when an operation ran that invalidates a forward replay past the
 *  anchor.  Read it BEFORE nde_history_insert_point_clear(). */
gboolean nde_history_insert_point_disturbed(void);

/** Disarm, discard the undo stash (and those records' checkpoints), and
 *  return the ids inserted so far in log order.  Never NULL; the caller
 *  g_array_unref()s.  The inserted records themselves stay in the log —
 *  nde_history_drop_records() is how an abandoned insertion undoes them. */
GArray  *nde_history_insert_point_clear(void);

/** Remove the live records with these ids from the log, dropping their
 *  output checkpoints and pool entries.  Ids that are not live are skipped.
 *  Does not touch the dead tail. */
void     nde_history_drop_records(GArray *ids);

/** Drop the dead tail (records undone but not yet superseded by an append)
 *  and release their output checkpoints. */
void     nde_history_truncate_dead(void);

/** Stale flag accessors (load-time PIXHASH mismatch). */
void     nde_history_set_stale(gboolean stale);
gboolean nde_history_is_stale(void);

/**
 * Schedule a panel refresh.  Phase-1 placeholder until the History panel
 * lands (work-order step 8): currently a no-op kept so capture sites are
 * already wired.
 */
void nde_history_notify_panel(void);

/* Plain → FLIS promote: move every record targeting @from_item onto
 * @to_item (see flis_promote_from_gfit). */
void nde_history_rebind_item(gint from_item, gint to_item);

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
/**
 * The item a capture belongs to: the active layer normally, or the item an
 * open insertion has borrowed the display for.  Every capture site must use
 * this rather than reading the active layer itself — two sites answering the
 * question separately is what filed a dialog-driven op against the wrong
 * layer while an insertion was open.
 */
gint nde_capture_target_item(void);

gint64 nde_capture_structural(const char *op_id, gint scope,
                              gint target_item_id, gchar *params,
                              const char *summary);

/** The op_id of the record every image's history opens with. */
#define NDE_OP_IMAGE_ORIGIN "image.origin"

/**
 * nde_capture_image_origin:
 * @kind: how the image came to exist — "file", "stack" or "generated".
 * @src:  the file it was read from, or the label its producer gave it.
 * @summary: the history line, already translated.
 *
 * The FIRST step of a plain image's history: where its pixels came from.
 * Every other item already says this for itself — a layer's history opens
 * with the layer.add or layer.duplicate that made it, an item born of a merge
 * with the composite that produced it — and the image was the one item whose
 * origin went unrecorded.  So a freshly opened file had no history at all and
 * no node in the graph until something was done to it: add a second layer and
 * the panel showed the NEW layer and not the one that was already there.
 *
 * DOCUMENT-scope and structural, so it is never a chain member: what a replay
 * starts from is the baseline checkpoint, and this record only says where that
 * baseline came from.  It is therefore neither amendable nor deletable, which
 * is right — where an image came from is not a step you can edit.
 *
 * A no-op when the document already has a history: a FLIS file brings its own,
 * and promoting a plain image REBINDS this record onto the base layer
 * (flis_promote_from_gfit) rather than wanting a second one.  Also a no-op
 * during replay, which must not grow the log it is reproducing.
 *
 * Returns the new record id, or 0 when nothing was recorded.
 */
gint64 nde_capture_image_origin(const char *kind, const char *src,
                                const char *summary);

/** One input pin, as a capture site states it. */
typedef struct {
	const char *role;
	gint        src_item_id;
	gint64      src_record_id;   /* 0 = that item's baseline */
} nde_pin_spec;

/**
 * nde_capture_structural with input pins — for the composite node, whose
 * inputs are the whole point (nde_composite.h).  @pins is copied; n_pins may
 * be 0, in which case this is exactly nde_capture_structural.
 */
gint64 nde_capture_structural_pinned(const char *op_id, gint scope,
                                     gint target_item_id, gchar *params,
                                     const char *summary,
                                     const nde_pin_spec *pins, guint n_pins);

/**
 * Build+append an opaque (Tier B, params == NULL) record and return the new
 * record id.  For mutations recorded for provenance only, with no replayable
 * parameters (python pixel/mask writes — sketch §13.2).  @summary and @op_id
 * are copied.  Same success-only capture discipline as nde_capture_structural.
 *
 * @post (may be NULL) is the POST-op image at capture: an opaque record is a
 * barrier, so when @post is non-NULL its pixels are stored as the record's
 * output checkpoint (the restart point that keeps the tail editable —
 * nde-phase4 P4.3).  Pass NULL when no post-op image is available.
 */
gint64 nde_capture_opaque(const char *op_id, gint scope,
                          gint target_item_id, const char *summary,
                          const fits *post);

/**
 * Build+append a replayable-script (Tier C, phase 5) record and return the
 * new record id.  For a declaring Python script's committed effect: @params
 * is the kv blob carrying the script path, its sha256 and the recorded arg
 * values (ownership transferred, may be NULL — then it is a plain barrier).
 * @post is the committed POST-script image, stored as the output checkpoint
 * (Tier C is always a barrier / restart point — scripts are slow, so
 * downstream edits restart here).  Same success-only discipline as the other
 * capture helpers.
 */
gint64 nde_capture_script(const char *op_id, gint scope,
                          gint target_item_id, gchar *params,
                          const char *summary, const fits *post);

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
 *
 * @post (may be NULL) is the POST-op image at capture: when the resolved
 * record is a barrier (Tier B, or Tier A whose mask was not kept) and @post
 * is non-NULL, its pixels are stored as the record's output checkpoint
 * (nde-phase4 P4.3).  Non-barrier records ignore @post.
 *
 * @mask_aware says whether the caller's worker run blended through the
 * processing mask (args->mask_aware on the apply).  TRUE + a live active
 * mask records the mask as a pinned input with its pixels kept, so the
 * step stays replayable; FALSE captures an unmasked record even when a mask
 * is sitting active on the image — because the op never read it.
 */
struct op_descriptor;
gint64 nde_capture_from_descriptor(const struct op_descriptor *op,
                                   gconstpointer params, const char *summary,
                                   const fits *post, gboolean mask_aware);

/**
 * nde_capture_from_descriptor_for_item:
 *
 * Like nde_capture_from_descriptor() but records against an EXPLICIT layer
 * item (LAYER scope, no mask pinning).  For multi-layer operations that
 * factor into independent per-layer records — layers match and group colour
 * calibration capture one record per affected layer, each targeting that
 * layer's item_id rather than the active layer.
 */
gint64 nde_capture_from_descriptor_for_item(const struct op_descriptor *op,
                                            gconstpointer params,
                                            const char *summary,
                                            const fits *post, gint item_id);

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

/**
 * Streamed SHA-256 of the file at @path (lowercase hex of the raw FILE bytes,
 * NOT the decoded pixels — Convention 1 of the external-input serializers).
 * On success returns a newly-allocated hex string (caller g_free()s) and,
 * when @size_out is non-NULL, stores the file size in bytes.  Returns NULL
 * (and leaves @size_out untouched) when the file is missing or unreadable.
 * Reads in bounded chunks so arbitrarily large operands hash in constant RAM.
 */
gchar *nde_file_sha256(const char *path, gint64 *size_out);

/**
 * Verify a pinned file operand for replay (Convention 1): re-hash the file at
 * @path and check it against @expect_size / @expect_sha256 (lowercase hex).
 * Returns TRUE only when the file exists, its size matches, and its SHA-256
 * matches.  On any failure logs a clear siril_log_error naming @path first
 * ("operand file missing or changed") so the replay driver's generic
 * "replay preparation failed" is preceded by an actionable line.  A NULL or
 * empty @path fails.  @expect_sha256 NULL/empty skips the hash comparison
 * (size-only) but the streamed read still detects unreadable files.
 */
gboolean nde_operand_verify(const char *path, gint64 expect_size,
                            const char *expect_sha256);

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
