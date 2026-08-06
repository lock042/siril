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
#ifndef _NDE_JOINT_H_
#define _NDE_JOINT_H_

/**
 * \file nde_joint.h
 * \brief JOINT multi-layer records: one history entry, N participating chains.
 *
 * The multi-channel colour operations (layers match, group colour
 * calibration) analyse SEVERAL layers together and scale each of them.
 * Recording that as one flis.layer_scale per layer froze the RESULT — the
 * factors — while the operation's meaning is its INTENT: neutralise the
 * composite background given the layers and tints as they now stand.  A
 * joint record therefore stores the operation's parameters and recomputes
 * the analysis at replay, so amending anything upstream of any participant
 * (the canonical case: the tint assigned to one narrowband layer) re-derives
 * fresh factors instead of replaying stale ones.
 *
 * Encoding: an ordinary NDE_SCOPE_LAYER Tier-A record targeting an ANCHOR
 * participant (the first recorded one), with an input pin per participant
 * ("in0"…, anchor included) and the participant list in params.  Chain
 * membership is extended (nde_replay.c): the record is a member of EVERY
 * participant's chain — one record at one log position, shared by N chains,
 * which also makes undo/redo across it a single step by construction.
 *
 * Replay (the engine special-cases these ops the way it does composites):
 * when a participant's chain reaches the record, the OTHER participants are
 * resolved positionally — their chain prefixes strictly BEFORE the record's
 * log position, so an edit-at insertion between a sibling's last step and
 * the record is honoured — their compositing state is folded positionally
 * too (nde_compositing_fold_upto), the analysis is re-run, and only the
 * replaying participant's factor is applied to the accumulated scratch.
 * The factors are cached per (record, history generation) so one edit plus
 * its cascades costs a single analysis however many participants replay.
 */

#include "core/siril.h"   /* fits, destructor */
#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct nde_record;
struct op_descriptor;

extern const struct op_descriptor op_desc_flis_layers_match;

/* One participating layer, as recorded in the joint record's params.  The
 * tint fields are the CAPTURE-TIME fallback: replay folds the live history
 * for the layer's tint as of the record's position and only falls back to
 * these when the log says nothing about it (a document recorded before tint
 * capture existed).  diag_* are the last factors the analysis produced —
 * diagnostics for display, never replayed while the analysis can run. */
typedef struct {
	gint     item_id;
	gchar   *name;         /* owned; label for summaries and the graph */
	gboolean tinted;
	double   tint[3];
	double   diag_scale;
	double   diag_offset;
} nde_joint_participant;

/* Params struct for "flis.layers_match" (background neutralisation across
 * tinted mono layers).  destroy_fn first, like every op user struct. */
struct nde_joint_layers_match_data {
	destructor             destroy_fn;
	guint                  n;
	nde_joint_participant *parts;   /* owned array of n */
};

struct nde_joint_layers_match_data *nde_joint_layers_match_data_new(guint n);
void nde_joint_layers_match_data_free(void *p);

/* Params struct for "flis.group_calibration" (CC / PCC / SPCC run on a
 * group composite and distributed into per-layer factors).  @kind follows
 * flis_group_calib_kind, plus DIRECT for the raw-affine API where the K/O
 * pair itself IS the operation parameter (no analysis to re-run).  diag_K /
 * diag_O are the composite-level affine the analysis last produced — the L1
 * fallback that still adapts the DISTRIBUTION to current tints and medians
 * when the full analysis cannot re-run (offline, no composite). */
#define NDE_JOINT_GROUP_DIRECT 3

struct nde_joint_group_calib_data {
	destructor             destroy_fn;
	guint                  n;
	nde_joint_participant *parts;    /* owned array of n */
	gint      kind;                  /* flis_group_calib_kind or DIRECT */
	gint      group_id;
	gboolean  manual;                /* CC: sliders are the parameters */
	double    manual_kw[3];
	rectangle white_sel, black_sel;  /* CC: sampling rectangles */
	double    low, high;
	gchar    *pcc_blob;              /* owned nested color.photometric_cc
	                                  * params blob, or NULL (CC/DIRECT) */
	double    diag_K[3], diag_O[3];
};

struct nde_joint_group_calib_data *nde_joint_group_calib_data_new(guint n);
void nde_joint_group_calib_data_free(void *p);

extern const struct op_descriptor op_desc_flis_group_calibration;

/* One participating layer of "flis.register".  Registration is a GEOMETRIC
 * joint op, so a participant's recorded result is a warp and a framing, not
 * a scalar: @H is the transform that was actually handed to cvTransformImage
 * (post-framing, i.e. already carrying the FRAMING_MAX bounding-box shift),
 * @out_rx/@out_ry the size that warp produced, and @pos_x/@pos_y where the
 * layer landed on the canvas afterwards.  @geom_sig is the digest of every
 * geometry-changing step upstream of the record on this layer AS OF capture:
 * replay compares it against a freshly computed one to decide whether the
 * stored transform is still the right answer (L1) or the solve has to be
 * re-run (L2) — see nde_joint_geometry_signature. */
typedef struct {
	gint     item_id;
	gchar   *name;           /* owned; label for summaries and the graph */
	double   H[9];           /* solved homography, row-major */
	gint     pos_x, pos_y;   /* resulting canvas offset */
	gint     out_rx, out_ry; /* framed output size the warp produced */
	gchar   *geom_sig;       /* owned; signature at capture */
} nde_joint_register_participant;

/* Params struct for "flis.register" (multi-layer layer-to-layer
 * registration).  It stores BOTH halves of the operation — the solved
 * transforms and the settings that produced them — because neither alone
 * replays correctly: replaying the settings always would re-solve (and so
 * re-jitter the stars) even when nothing upstream moved, and replaying the
 * transforms always would silently apply a stale warp once a rotation or a
 * crop was inserted upstream of a participant. */
struct nde_joint_register_data {
	destructor  destroy_fn;
	guint       n;
	nde_joint_register_participant *parts;   /* owned array of n */
	gint        method;         /* flis_reg_method_id */
	gint        tx_type;        /* transformation_type */
	gint        interpolation;  /* opencv_interpolation */
	gboolean    clamp;
	gint        ref_item;       /* reference layer's item id */
	rectangle   selection;      /* com.selection at capture */
	gint        canvas_w, canvas_h;  /* canvas AFTER registration */
};

struct nde_joint_register_data *nde_joint_register_data_new(guint n);
void nde_joint_register_data_free(void *p);

extern const struct op_descriptor op_desc_flis_register;

/** TRUE for the joint multi-layer op ids ("flis.layers_match",
 *  "flis.group_calibration", "flis.register"). */
gboolean nde_joint_is_op(const char *op_id);

/** TRUE for the joint ops that WARP their participants rather than scale
 *  them, so the replay driver must thread the layer position through them
 *  and the chain builder must treat the chain as carrying geometry. */
gboolean nde_joint_is_geometric_op(const char *op_id);

/**
 * A digest of every geometry-changing step in @item_id's chain strictly
 * BEFORE @before_record_id's log position (0 = the whole live log).  Two
 * equal digests mean the layer reaches that point with the same geometry, so
 * a stored transform is still valid; any difference means the solve must be
 * re-run.  Never NULL — the empty string when the layer has no geometry
 * upstream at all.  Caller g_free()s.
 */
gchar *nde_joint_geometry_signature(gint item_id, gint64 before_record_id);

/**
 * Replay one participant of a "flis.register" record: warp @scratch in place
 * (L1 with the stored transform, L2 with a freshly solved one) and report
 * where the layer lands via @pos_x / @pos_y.  @item_id is the chain being
 * replayed; it must be one of the record's participants.  FALSE + heap @err
 * when the record does not parse or the warp fails.
 */
gboolean nde_joint_register_apply(const struct nde_record *rec,
                                  fits *scratch, gint item_id,
                                  gint *pos_x, gint *pos_y, gchar **err);

/** TRUE when @rec is a joint record and @item_id is one of its recorded
 *  participants (parsed from params).  The chain-membership extension. */
gboolean nde_joint_record_names_item(const struct nde_record *rec, gint item_id);

/** TRUE when two params blobs name the same participants in the same order.
 *  An amend that changes the list takes the subset-amend path (nde_replay.c):
 *  pins re-derived, removed participants reverted, added ones scaled. */
gboolean nde_joint_params_same_participants(const char *a, const char *b);

/** The participant item ids named by a params blob (the record-less variant
 *  of nde_joint_record_participants, for validating an amend's NEW params).
 *  Caller g_free()s; NULL when the blob does not parse. */
gint *nde_joint_params_participants(const char *params, guint *n_out);

/**
 * The participant item ids recorded in @rec's params, in recorded order, or
 * NULL when @rec is not a joint record / its params do not parse.  Caller
 * g_free()s; *n_out receives the count.
 */
gint *nde_joint_record_participants(const struct nde_record *rec, guint *n_out);

/**
 * Drop the factor cache.  Called by edit_execute before building a TRIAL
 * chain: the trial replays against a not-yet-committed log, so factors
 * computed for it must not survive into (or come from) the committed
 * generation.  The committed log's own mutations are covered by the history
 * generation counter without this.
 */
void nde_joint_cache_invalidate(void);

/**
 * The factors the joint record @rec derives for ALL its participants, in
 * recorded participant order: on success *a_out / *b_out are heap arrays of
 * n scale / offset values (caller g_free()s both).  @self_item is the chain
 * being replayed and @self_scratch its accumulated pre-record state, used
 * as that participant's pixels instead of resolving its chain (exactly as
 * composite_apply treats its base).  Other participants are resolved
 * positionally to just before the record.  Results are cached per (record,
 * history generation) — cache hits return copies and never touch pixels.
 * FALSE + heap @err when a participant cannot be rebuilt or the analysis
 * fails with no stored fallback.
 */
gboolean nde_joint_factors(const struct nde_record *rec, fits *self_scratch,
                           gint self_item, double **a_out, double **b_out,
                           gchar **err);

/** Lifetime count of actual analysis runs (factor-cache misses) — for the
 *  tests that prove one edit costs one analysis regardless of participant
 *  count.  Monotonic, never reset. */
guint nde_joint_analysis_runs(void);

/**
 * The single (scale, offset) the joint record @rec derives for @self_item —
 * nde_joint_factors plus the participant-index lookup, packaged for the
 * replay driver.  Same contract and caching.
 */
gboolean nde_joint_factor_for_item(const struct nde_record *rec,
                                   fits *self_scratch, gint self_item,
                                   double *scale_out, double *offset_out,
                                   gchar **err);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_JOINT_H_ */
