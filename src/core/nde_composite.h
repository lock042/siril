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
 * Merge-down, group merge and flatten are the same computation — the
 * compositor — differing only in how many inputs they consume and at what
 * point the result is committed.  This module holds the part of that node the
 * history needs: the per-input compositing state, serialized into the record,
 * and a render that takes resolved pixels rather than live layers.
 *
 * WHY THE STATE MUST BE SERIALIZED.  A merged-away input's pixels are
 * recoverable — its baseline and its records both survive the merge, so
 * replaying its chain reproduces them exactly.  What does NOT survive is the
 * layer struct: its blend mode, opacity, canvas offset and tint die with it,
 * and those are inputs to the composite just as much as the pixels are.  So
 * the record carries them.  The surviving (base) input's own parameters are
 * deliberately absent: merge-down RETAINS them on the survivor, so baking
 * them here would apply them twice (see flis_render_layers_merge).
 *
 * LIVE EDGES, NOT PINNED COPIES.  Both inputs are resolved by replaying their
 * chains at replay time rather than read from a stored copy.  That costs a
 * replay per re-run and buys two things: no extra checkpoint storage per
 * merge, and — the reason it matters — an amend anywhere upstream is picked up
 * automatically, with no pin to refresh and no cascade to keep in step.  Mask
 * edges chose stored copies because a mask op's input is an image state that
 * may be many steps back; a composite's inputs are whole items, which the
 * chain machinery already knows how to rebuild.
 *
 * BACKWARD COMPATIBILITY.  Records written before this existed carry no input
 * pins and no per-input state.  nde_composite_record_replayable() is false for
 * them, and the chain builder keeps treating them as it always did: a hard
 * blocker naming what it cannot rebuild.  Degrading honestly beats guessing at
 * an opacity nobody recorded.
 *
 * Threading: the render is pure.  Encoding reads live layer state, so it runs
 * at the capture site; decoding and rendering run on the replay conductor.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;
struct _flis_layer_t;
struct nde_record;

/** Role names for a composite's input pins.  "base" is the surviving item the
 *  record targets; "overlay" is the input consumed by the merge. */
#define NDE_COMPOSITE_ROLE_BASE    "base"
#define NDE_COMPOSITE_ROLE_OVERLAY "overlay"

/**
 * One input's compositing state — the fields of a flis_layer_t that decide how
 * its pixels reach the composite.  Everything else about the layer (its name
 * aside) is either reproduced by its chain or irrelevant here.
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
} nde_composite_input;

void nde_composite_input_clear(nde_composite_input *in);

/** TRUE for the op ids that ARE composite nodes. */
gboolean nde_composite_is_op(const char *op_id);

/**
 * Encode the compositing state of both inputs as the params blob for a
 * composite record.  Caller owns the result (it is handed straight to the
 * capture helper).
 *
 * @base contributes only its canvas offset and tint: merge-down RETAINS its
 * blend mode, opacity and visibility on the survivor, so recording them would
 * invite a second application.  @top contributes everything, because nothing
 * of it survives.
 */
gchar *nde_composite_params_encode(const struct _flis_layer_t *base,
                                   const struct _flis_layer_t *top);

/**
 * Decode a composite record's params.  Either output may be NULL.  Caller
 * clears what it asked for with nde_composite_input_clear().  FALSE when the
 * blob predates this format or is missing a field, in which case the record is
 * provenance only.
 */
gboolean nde_composite_params_decode(const char *params,
                                     nde_composite_input *base_out,
                                     nde_composite_input *top_out);

/**
 * TRUE when @rec is a composite record carrying BOTH the input pins and the
 * per-input state its re-run needs — that is, when it can be a replayable
 * chain member rather than a blocker.
 */
gboolean nde_composite_record_replayable(const struct nde_record *rec);

/**
 * Blend @overlay (state @ov) over @base (state @bs) and return a newly
 * allocated canvas-sized composite, or NULL + heap @err.  @base is painted
 * raw — tint baked, no blend mode, opacity or mask — because the surviving
 * layer keeps those; this is exactly flis_render_layers_merge's contract, and
 * it is that function which does the work.
 *
 * @base_x / @base_y override @bs's recorded offset when the caller has a live
 * one: a geometry step before the merge moves the base, so an amend of that
 * step must move where it lands here too.  Pass the recorded values (or
 * @bs->position_*) when there is no threaded position.
 */
struct ffit *nde_composite_render(const struct ffit *base,
                                  const nde_composite_input *bs,
                                  gint base_x, gint base_y,
                                  const struct ffit *overlay,
                                  const nde_composite_input *ov,
                                  gchar **err);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_COMPOSITE_H_ */
