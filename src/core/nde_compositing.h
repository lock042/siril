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
#ifndef _NDE_COMPOSITING_H_
#define _NDE_COMPOSITING_H_

/**
 * \file nde_compositing.h
 * \brief Live compositing-state history records (nde-graph-design-note.md §6.1).
 *
 * `layer.set_opacity`, `layer.set_blend` and `layer.set_visible` are NOT pixel
 * operations: they are inputs to the compositor, which reads them from live
 * FLIS layer state at render time.  They are therefore never members of a
 * replay chain (nde_chain_build skips them) — but they ARE editable, and until
 * now editing them did nothing at all.
 *
 * The model here is a FOLD, not a replay: a layer's compositing state is the
 * last value recorded for each property, and amending / deleting / reordering
 * such a record simply re-folds the log onto the layer and re-renders.  No
 * pixels move, no chain is built, nothing is replayed.
 *
 * THE LOG IS AUTHORITATIVE.  The fold starts from the creation defaults and
 * applies every record, so a layer's compositing state is exactly what its
 * history says it is — including when the last record for a property is
 * DELETED, which must return the property to its default rather than strand
 * it at the deleted record's value.  The alternative (leaving properties the
 * log does not mention at their live value) was tried and rejected: it lets
 * the layer sit in a state no part of the log describes, so a later save
 * writes FLIS_META that disagrees with FLIS_HIST.
 *
 * The cost of that choice is a narrow one: a layer carrying a non-default
 * property with NO record to explain it (a file from a build predating the
 * record, or one whose history was truncated) is folded back to the default
 * the first time any compositing record on that layer is edited.  Files
 * written by this build always record what they set, so the fold reproduces
 * their FLIS_META exactly.  A per-item compositing baseline would close the
 * gap properly; that is the step-4 input-port generalisation, where the
 * baseline is simply another input (nde-graph-design-note.md §6).
 *
 * RESET POINTS.  Two structural records restore an item to the compositing
 * defaults and so make every property "known" from that point on:
 * `document.flatten` (which explicitly resets blend/opacity/visible on the
 * surviving base) and `layer.add` (a fresh layer starts at the defaults).
 * `layer.merge_down` is deliberately NOT a reset — it preserves the surviving
 * bottom layer's blend and opacity — and `layer.duplicate` copies the source
 * layer's properties, so neither makes the defaults knowable.
 *
 * Threading: reads the history through nde_history_snapshot() and mutates
 * live FLIS layer state, so it must run on a thread that owns the processing
 * slot (the edit job) or the main thread.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** TRUE for the three compositing-state op ids. */
gboolean nde_compositing_is_op(const char *op_id);

/**
 * Validate a candidate params blob for @op_id (range- and enum-checked), the
 * compositing counterpart of an op descriptor's deserialize round-trip.
 * FALSE + heap message in @err when the blob lacks the key or is out of range.
 */
gboolean nde_compositing_validate(const char *op_id, const char *params, gchar **err);

/**
 * Re-fold the live history onto layer @item_id's compositing state and
 * invalidate the composite.  A no-op returning TRUE for @item_id < 0 (a plain
 * image has no compositing state).  FALSE + heap @err when the layer is gone.
 */
gboolean nde_compositing_recompute(gint item_id, gchar **err);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_COMPOSITING_H_ */
