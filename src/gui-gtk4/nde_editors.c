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

static const struct {
	const char *op_id;
	nde_editor_open_fn open;
} editors[] = {
	{ "stretch.asinh",        asinh_open_amend },
	{ "stretch.curves",       curves_open_amend },
	{ "stretch.mtf",          histogram_mtf_open_amend },
	{ "stretch.ghs",          histogram_ghs_open_amend },
	{ "filters.scnr",         scnr_open_amend },
	{ "filters.median",       median_open_amend },
	{ "color.saturation",     satu_open_amend },
	{ "bkg.remove_gradient",  bge_open_amend },
	{ "filters.epf",          epf_open_amend },
	{ "filters.clahe",        clahe_open_amend },
	{ "filters.denoise",      denoise_open_amend },
	{ "filters.banding",      banding_open_amend },
	{ "color.ccm",            ccm_open_amend },
	{ "filters.cosmetic",     cosmetic_open_amend },
	{ "filters.rgradient",    rgradient_open_amend },
	/* Photometric CC / SPCC re-open their real dialog pre-filled; the amend
	 * commits without a live preview — the pipeline runs at commit, against
	 * the record's embedded star catalogue. */
	{ "color.photometric_cc",    pcc_open_amend },
	/* JOINT multi-layer records (nde_joint.h): group PCC/SPCC route to the
	 * same real dialog (nde_joint_open_amend dispatches); layers match has
	 * no parameters of its own and CC keeps the compact selections window.
	 * Apply-on-OK throughout — the amend-preview machinery synthesizes ONE
	 * target's state and a joint record writes to many. */
	{ "flis.layers_match",       nde_joint_open_amend },
	{ "flis.group_calibration",  nde_joint_open_amend },
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
		if (!g_strcmp0(editors[i].op_id, op_id))
			return editors[i].open(record_id);
	}
	return FALSE;
}
