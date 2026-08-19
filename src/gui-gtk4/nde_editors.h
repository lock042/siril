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
#include <gtk/gtk.h>

struct op_descriptor;

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
 *  An item a merge or flatten consumed keeps its editors: it has no layer to
 *  be made active, so the live-preview ones BORROW the display for the
 *  duration (nde_replay.h) and the change reaches the image through the
 *  composites that consumed the item. */
gboolean nde_editor_open(const gchar *op_id, gint64 record_id);

/** Show or hide a dialog's amend banner, and while showing it, append what a
 *  region preview would do with the steps AFTER the one being edited —
 *  recompute them live, or not, and why (nde_replay.h).  The .ui files carry
 *  the base sentence; this owns the rest, so the regime is written once.
 *
 *  @op is the descriptor this dialog will region-preview with in amend mode,
 *  or NULL when it offers no region preview there (it clears the ROI on
 *  entry).  Passing it is what keeps the banner honest: the tail can be
 *  perfectly region-replayable while THIS dialog has no rectangle to show it
 *  in, and promising a live recompute that cannot happen is worse than saying
 *  nothing.  It is a parameter rather than a read of
 *  gui.roi.operation_supports_roi because the dialogs call this before their
 *  startup declares the op.
 *
 *  Call wherever the dialog used to do gtk_widget_set_visible() on it. */
void nde_amend_note_update(GtkWidget *note, gboolean amend_mode,
                           const struct op_descriptor *op);

#endif /* SRC_GUI_NDE_EDITORS_H_ */
