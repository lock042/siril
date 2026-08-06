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
#ifndef _NDE_JOINT_EDITOR_H_
#define _NDE_JOINT_EDITOR_H_

#include <glib.h>

/* The dedicated editor for the JOINT multi-layer records (nde_joint.h):
 * flis.layers_match and flis.group_calibration.  One window for both —
 * participants and their last-derived factors shown read-only, the
 * operation's own parameters (CC selections / limits / manual multipliers)
 * editable where the record has them, and a Re-run that amends with the
 * (possibly updated) parameters so the whole analysis recomputes.  No live
 * preview: a joint record writes to several layers and the amend-preview
 * machinery synthesizes one target's state, so changes apply on OK. */
gboolean nde_joint_open_amend(gint64 record_id);

/* The editor for flis.register, the GEOMETRIC joint record.  Its params carry
 * two very different halves — the settings the user chose and the transforms
 * those settings solved for — and only the settings belong in a dialog: the
 * homographies, framing sizes, canvas positions and geometry signatures are
 * machine-derived state that hand-editing corrupts.  So this window exposes
 * method, transformation, interpolation, clamp and reference layer, and
 * re-serializes everything else verbatim.  Changing a setting invalidates the
 * stored solve, which the amend records by poisoning the participants'
 * signatures so replay solves the alignment again. */
gboolean nde_register_open_amend(gint64 record_id);

#endif /* _NDE_JOINT_EDITOR_H_ */
