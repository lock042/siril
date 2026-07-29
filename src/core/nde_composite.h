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
#ifndef _NDE_COMPOSITE_H_
#define _NDE_COMPOSITE_H_

/**
 * \file nde_composite.h
 * \brief The composite node: N inputs blended into one (design note §3, §9.0).
 *
 * Merge-down and flatten are the same computation — the compositor — differing
 * in how many inputs they consume and in one contract detail (below).  This
 * module holds the part of that node the history needs: the document state the
 * composite consumed, serialized into the record, and a render that takes
 * resolved pixels rather than live layers.
 *
 * WHY THE STATE MUST BE SERIALIZED.  A consumed input's pixels are recoverable
 * — its baseline and its records both survive, so replaying its chain
 * reproduces them exactly.  What does NOT survive is the layer struct: blend
 * mode, opacity, canvas offset, tint, group membership all die with the layer,
 * and those are inputs to the composite just as much as the pixels are.  Nor
 * does the document state around it: the canvas size, its background colour,
 * and the groups whose members were pre-composited.  So the record carries all
 * of it, and the replay renders against that rather than against whatever the
 * document happens to look like now (flis_render_ctx).
 *
 * THE ONE CONTRACT DIFFERENCE.  Merge-down paints its bottom input RAW — tint
 * baked, but no blend mode, opacity or mask — because those are RETAINED on the
 * surviving layer and baking them in would apply them twice.  Flatten blends
 * every input, including the bottom one, over the canvas background, and then
 * resets the survivor's properties to the defaults.  That is @raw_first, and it
 * is the only thing separating the two.
 *
 * LIVE EDGES, NOT PINNED COPIES.  Inputs are resolved by replaying their chains
 * at replay time rather than read from a stored copy.  That costs a replay per
 * re-run and buys two things: no extra checkpoint storage per composite, and —
 * the reason it matters — an amend anywhere upstream is picked up
 * automatically, with no pin to refresh and no cascade to keep in step.  Mask
 * edges chose stored copies because a mask op's input is an image state that
 * may be many steps back; a composite's inputs are whole items, which the chain
 * machinery already knows how to rebuild.
 *
 * MASKS ARE THE EXCEPTION — A PINNED COPY.  An input's layer mask is stored at
 * the capture site and read back from the checkpoint store, not re-derived.
 * Not for want of a chain: a layer mask can be painted, or loaded from a file,
 * and then it has no records to replay at all.  A copy is the only thing that
 * always exists.  It costs one mono image per masked input, and the mask
 * cascade already refreshes such copies when the mask's own chain is edited,
 * so a mask built by ops still propagates.  Only the PARTICIPATING masks are
 * stored: merge-down paints its bottom input raw, so that one's mask never
 * reaches the composite.
 *
 * WHAT IS STILL REFUSED.  A record written before this format existed: it has
 * no pins and no state, and degrading honestly beats guessing at an opacity
 * nobody wrote down.  Likewise a masked input whose stored mask has been
 * evicted or was never taken.
 *
 * The layer's PROCESSING mask (fit->mask) is not part of this: the compositor
 * never reads it — it restricts operations, not compositing — so it neither
 * needs storing nor stands in the way.
 *
 * Threading: the render is pure.  Capture reads live layer and document state,
 * so it runs at the capture site; parsing and rendering run on the replay
 * conductor.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;
struct nde_record;

/** Role names for a merge-down's two input pins.  Flatten uses "in0", "in1", …
 *  and a masked input adds "mask0", "mask1", …  Roles are labels for the graph
 *  view: the replay matches pins to recorded inputs by item id (which the
 *  params carry for both the layer and its mask), so a role is never parsed. */
#define NDE_COMPOSITE_ROLE_BASE    "base"
#define NDE_COMPOSITE_ROLE_OVERLAY "overlay"

/**
 * One input's compositing state — the fields of a flis_layer_t that decide how
 * its pixels reach the composite.  Everything else about the layer (its name
 * aside, which the graph view shows) is either reproduced by its chain or
 * irrelevant here.
 */
typedef struct {
	gint     item_id;
	gchar   *name;          /* owned; may be NULL */
	gint     blend_mode;    /* flis_blend_mode_t */
	gdouble  opacity;
	gint     position_x;
	gint     position_y;
	gboolean visible;
	gboolean has_tint;
	gdouble  tint_r, tint_g, tint_b;
	gint     group_id;      /* 0 = ungrouped */
	gboolean was_masked;    /* carried an active layer mask into the composite */
	gint     mask_item_id;  /* the LMASK item it came from; 0 if unmasked */
} nde_composite_input;

/** A group the composite pre-composited, or whose visibility/opacity modified
 *  its members.  Flat: FLIS groups do not nest. */
typedef struct {
	gint     item_id;
	gint     blend_mode;
	gdouble  opacity;
	gboolean visible;
} nde_composite_group;

/** Everything a composite record consumed, decoded. */
typedef struct {
	gboolean raw_first;            /* the merge-down contract (see above) */
	guint    canvas_w, canvas_h;
	gdouble  bg_r, bg_g, bg_b;
	GArray  *inputs;               /* nde_composite_input, bottom-first */
	GArray  *groups;               /* nde_composite_group */
} nde_composite_state;

void nde_composite_state_free(nde_composite_state *st);

/** TRUE for the op ids that ARE composite nodes. */
gboolean nde_composite_is_op(const char *op_id);

/**
 * Decode a composite record's params, or NULL if they predate this format or
 * are missing a field — in which case the record is provenance only.  Caller
 * frees with nde_composite_state_free().
 *
 * Both the current indexed form and the first merge-down form (base_… / top_…)
 * are accepted; the latter decodes to the same two-input state.
 */
nde_composite_state *nde_composite_state_parse(const char *params);

/**
 * Validate @new_params as a replacement for @old_params on a composite record.
 * A composite has no op descriptor, so this stands in for a deserializer.
 *
 * What it accepts is the COMPOSITING state — each input's opacity, blend mode,
 * visibility, tint and offset, and the same for the groups.  Those are exactly
 * the parameters that died with the layers and are exactly what a user would
 * want to reconsider: "that layer was merged in at 50%, make it 70%".
 *
 * What it refuses is everything else.  Changing which items the node consumes
 * would be rewiring the graph, which no edit does (design note §5.4); changing
 * the canvas, the input count, or the raw-first contract would describe an
 * operation that never happened.
 */
gboolean nde_composite_validate(const char *old_params, const char *new_params,
                                gchar **err);

/**
 * @params with the input masked by @mask_item_id no longer masked, or NULL if
 * no input was.  Caller owns the result.
 *
 * Not something a user edit may do — nde_composite_validate() refuses it,
 * because the mask is part of what the step consumed.  It is what DELETING the
 * step that built the mask means: the mask no longer exists, so the composite
 * that used it composites without one.  Recorded rather than inferred, so the
 * log says what will happen instead of leaving a mask input pointing at
 * nothing and hoping the replay reads that as "no mask" rather than "lost".
 */
gchar *nde_composite_params_drop_mask(const char *params, gint mask_item_id);

/**
 * TRUE when @rec is a composite record carrying both the input pins and the
 * state its re-run needs — that is, when it can be a replayable chain member
 * rather than a blocker.
 */
gboolean nde_composite_record_replayable(const struct nde_record *rec);

/**
 * Render @st from @pixels — one entry per recorded input, in the same order.
 * An invisible input contributes nothing and its entry may be NULL; a visible
 * one may not.  @masks is parallel to it and may be NULL entirely: entry i is
 * the mono mask image for input i, as stored by the capture, and is required
 * wherever the recorded input says it was masked.  All entries are borrowed.
 * Returns a newly allocated canvas-sized composite, or NULL and a heap @err.
 *
 * To place the input that a chain is replaying, overwrite that input's
 * position_x/position_y in @st before calling: a geometry step before the
 * composite moves the layer, so amending that step must move where it lands
 * here too.
 */
struct ffit *nde_composite_render(const nde_composite_state *st,
                                  struct ffit *const *pixels,
                                  struct ffit *const *masks, gchar **err);

/**
 * Capture, in two halves.
 *
 * _begin() is called at the last moment every input still exists.  It records
 * their compositing state and the document's, pins each input to its latest
 * record, and ensures a BASELINE for each — a layer consumed without ever
 * having been edited has none, and afterwards its pre-composite pixels are
 * unreachable.  @layers is bottom-first; @raw_first selects the merge-down
 * contract.
 *
 * _commit() appends the record once the mutation is done, so that the undo
 * purge in between cannot leave it coupled to an entry that no longer exists.
 * It consumes @cap either way.  _free() is the abort path.
 */
typedef struct _nde_composite_capture nde_composite_capture;

nde_composite_capture *nde_composite_capture_begin(GSList *layers,
                                                   gboolean raw_first);
gint64 nde_composite_capture_commit(nde_composite_capture *cap,
                                    const char *op_id, gint target_item_id,
                                    const char *summary);
void nde_composite_capture_free(nde_composite_capture *cap);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_COMPOSITE_H_ */
