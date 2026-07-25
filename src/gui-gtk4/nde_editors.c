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
#include "algos/background_extraction.h"

typedef void (*nde_editor_open_fn)(gint64 record_id);

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
};

gboolean nde_editor_open(const gchar *op_id, gint64 record_id) {
	if (!op_id)
		return FALSE;
	for (guint i = 0; i < G_N_ELEMENTS(editors); i++) {
		if (!g_strcmp0(editors[i].op_id, op_id)) {
			editors[i].open(record_id);
			return TRUE;
		}
	}
	return FALSE;
}
