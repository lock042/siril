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
#ifndef SRC_GUI_NDE_EDITORS_H_
#define SRC_GUI_NDE_EDITORS_H_

#include <glib.h>

/*
 * Native amend editors (convergence C4).
 *
 * Maps an op_descriptor id to the dialog that can edit a history record
 * of that op in place (amend mode: pre-filled widgets, live preview
 * against the pre-record state, Apply routed through the amend path).
 * This table lives GUI-side by design — core must not reference GUI
 * code, so it is NOT a field on op_descriptor (recorded deviation from
 * sketch §7 decision 12).
 *
 * The History panel's Edit button consults it; the generic kv-grid
 * editor remains the fallback for unregistered ops.
 */

/** Open the native editor for a history record of op @op_id, if one is
 *  registered.  Returns TRUE when the record is handled by a native
 *  editor (even if entering amend mode was refused — the opener logs the
 *  reason); FALSE when the caller should fall back to the kv grid: no
 *  editor is registered for the op, or the editor vetoed this particular
 *  record (e.g. an unlinked MTF stretch the histogram sliders cannot
 *  represent).
 *
 *  @target_is_retained_input says the record's item has been merged or
 *  flattened away and survives only as a composite's input.  That rules out
 *  the LIVE-PREVIEW editors — there is no on-screen image to preview against
 *  — but not the apply-on-OK ones, which read the record's parameters and
 *  commit.  Deciding it here rather than at the call site is what keeps the
 *  distinction next to the table that knows which editor is which: gating
 *  the whole registry on it sent a flattened document's SPCC step to the raw
 *  kv grid. */
gboolean nde_editor_open(const gchar *op_id, gint64 record_id,
                         gboolean target_is_retained_input);

#endif /* SRC_GUI_NDE_EDITORS_H_ */
