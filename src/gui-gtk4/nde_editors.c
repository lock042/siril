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

#include "gui-gtk4/nde_editors.h"
#include "core/siril.h"
#include "core/nde/nde_replay.h"
#include "core/op_descriptor.h"
#include "gui-gtk4/asinh.h"
#include "gui-gtk4/curves.h"
#include "gui-gtk4/histogram.h"
#include "gui-gtk4/scnr.h"
#include "gui-gtk4/median.h"
#include "gui-gtk4/saturation.h"
#include "gui-gtk4/epf.h"
#include "gui-gtk4/clahe.h"
#include "gui-gtk4/denoisegui.h"
#include "gui-gtk4/banding.h"
#include "gui-gtk4/colors.h"
#include "gui-gtk4/cosmetic_correction.h"
#include "gui-gtk4/rgradient.h"
#include "gui-gtk4/nde_joint_editor.h"
#include "gui-gtk4/photometric_cc.h"
#include "algos/background_extraction.h"

/* Returns TRUE when the editor takes the record (even if entering amend
 * mode is then refused — the core logs why); FALSE vetoes it, sending the
 * caller to the kv-grid fallback (e.g. an unlinked MTF record the linked
 * histogram sliders cannot represent). */
typedef gboolean (*nde_editor_open_fn)(gint64 record_id);

/* A merge or a flatten does not cost an op its editor.  The live-preview
 * editors here amend against the display, and an item those operations
 * consumed has no layer to be made active — so the core lends it the display
 * for the duration (nde_replay.h) and the composites that consumed the item
 * are recomputed on the way out.  This table used to record which editors
 * previewed, in order to send exactly those to the kv grid once the layer was
 * gone; flattening a document then demoted every stretch in it to raw
 * key/value text, which is what the borrow removes. */
static const struct {
	const char *op_id;
	nde_editor_open_fn open;
} editors[] = {
	{ "stretch.asinh",        asinh_open_amend         },
	{ "stretch.curves",       curves_open_amend        },
	{ "stretch.mtf",          histogram_mtf_open_amend },
	{ "stretch.ghs",          histogram_ghs_open_amend },
	{ "filters.scnr",         scnr_open_amend          },
	{ "filters.median",       median_open_amend        },
	{ "color.saturation",     satu_open_amend          },
	{ "bkg.remove_gradient",  bge_open_amend           },
	{ "filters.epf",          epf_open_amend           },
	{ "filters.clahe",        clahe_open_amend         },
	{ "filters.denoise",      denoise_open_amend       },
	{ "filters.banding",      banding_open_amend       },
	{ "color.ccm",            ccm_open_amend           },
	{ "filters.cosmetic",     cosmetic_open_amend      },
	{ "filters.rgradient",    rgradient_open_amend     },
	/* Photometric CC / SPCC re-open their real dialog pre-filled; the amend
	 * commits without a live preview — the pipeline runs at commit, against
	 * the record's embedded star catalogue. */
	{ "color.photometric_cc",    pcc_open_amend        },
	/* JOINT multi-layer records (nde_joint.h): group PCC/SPCC route to the
	 * same real dialog (nde_joint_open_amend dispatches); layers match has
	 * no parameters of its own and CC keeps the compact selections window.
	 * Apply-on-OK throughout — the amend-preview machinery synthesizes ONE
	 * target's state and a joint record writes to many. */
	{ "flis.layers_match",       nde_joint_open_amend  },
	{ "flis.group_calibration",  nde_joint_open_amend  },
	/* Registration is joint too, but its params are half machine-derived
	 * state (transforms, framing, signatures), so it gets a window that shows
	 * only the settings — the kv-grid fallback would have offered the nine
	 * homography coefficients per layer as editable text. */
	{ "flis.register",           nde_register_open_amend },
};

gboolean nde_editor_open(const gchar *op_id, gint64 record_id) {
	if (!op_id)
		return FALSE;
	for (guint i = 0; i < G_N_ELEMENTS(editors); i++) {
		if (g_strcmp0(editors[i].op_id, op_id))
			continue;
		return editors[i].open(record_id);
	}
	return FALSE;
}

/* ---- the amend banner ---------------------------------------------------
 * Every amend-capable dialog carries the same GtkLabel from its .ui file,
 * shown while the dialog is in amend mode.  Its base sentence — later steps
 * are hidden and recomputed on apply — stopped being the whole truth when
 * region previews learned to replay the tail (nde_replay.h), so the regime is
 * appended here rather than duplicated fourteen times in the .ui files.
 *
 * Stated once, at open, and phrased so it holds whether or not a region is
 * drawn: the dialogs do not re-run this from their ROI callbacks, and a
 * sentence that went stale the moment the user dragged a selection would be
 * worse than no sentence at all. */
void nde_amend_note_update(GtkWidget *note, gboolean amend_mode,
                           const op_descriptor *op) {
	if (!note)
		return;
	gtk_widget_set_visible(note, amend_mode);
	if (!amend_mode || !GTK_IS_LABEL(note))
		return;

	/* Keep the .ui text: this runs on every open, and appending to the
	 * previous result would grow the label without bound. */
	const gchar *base = g_object_get_data(G_OBJECT(note), "nde-amend-base");
	if (!base) {
		base = g_strdup(gtk_label_get_text(GTK_LABEL(note)));
		g_object_set_data_full(G_OBJECT(note), "nde-amend-base",
		                       (gpointer)base, g_free);
	}

	gchar *why = NULL;
	gchar *text;
	if (!op_descriptor_is_roi_capable(op)) {
		/* This dialog shows no rectangle in amend mode; the regime is not
		 * the user's to reach, so there is nothing to tell them about it. */
		text = g_strdup(base);
	} else if (nde_region_tail_available(NULL, &why)) {
		text = g_strconcat(base, " ",
		                   _("A region preview recomputes them live inside the "
		                     "region outline."), NULL);
	} else if (why) {
		gchar *tail = g_strdup_printf(_("A region preview cannot include them: %s."),
		                              why);
		text = g_strconcat(base, " ", tail, NULL);
		g_free(tail);
	} else {
		text = g_strdup(base);   /* no amend preview installed — nothing to add */
	}
	g_free(why);
	gtk_label_set_text(GTK_LABEL(note), text);
	g_free(text);
}
