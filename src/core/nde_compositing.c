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

/* See nde_compositing.h for the model (a fold, not a replay). */

#include "core/siril.h"
#include "core/nde_compositing.h"
#include "core/nde_history.h"
#include "core/gui_iface.h"
#include "io/image_format_flis.h"

#define OP_SET_OPACITY "layer.set_opacity"
#define OP_SET_BLEND   "layer.set_blend"
#define OP_SET_VISIBLE "layer.set_visible"

/* Compositing defaults, matching flis_layer_new() and the reset performed by
 * flis_flatten_all() on the surviving base layer. */
#define COMP_DEFAULT_OPACITY 1.0f
#define COMP_DEFAULT_BLEND   FLIS_BLEND_NORMAL
#define COMP_DEFAULT_VISIBLE TRUE

gboolean nde_compositing_is_op(const char *op_id) {
	return op_id && (!g_strcmp0(op_id, OP_SET_OPACITY) ||
	                 !g_strcmp0(op_id, OP_SET_BLEND) ||
	                 !g_strcmp0(op_id, OP_SET_VISIBLE));
}

/* Every layer blend mode.  PASS_THROUGH is excluded on purpose: it is a group
 * mode (see the enum comment in image_format_flis.h) and would make a layer
 * composite in a way the renderer does not define. */
gboolean nde_compositing_blend_valid(gint64 v) {
	switch ((flis_blend_mode_t)v) {
	case FLIS_BLEND_NORMAL:
	case FLIS_BLEND_MULTIPLY:
	case FLIS_BLEND_SCREEN:
	case FLIS_BLEND_OVERLAY:
	case FLIS_BLEND_SOFT_LIGHT:
	case FLIS_BLEND_HARD_LIGHT:
	case FLIS_BLEND_COLOR_DODGE:
	case FLIS_BLEND_COLOR_BURN:
	case FLIS_BLEND_DARKEN:
	case FLIS_BLEND_LIGHTEN:
	case FLIS_BLEND_DIFFERENCE:
	case FLIS_BLEND_EXCLUSION:
	case FLIS_BLEND_HUE:
	case FLIS_BLEND_SATURATION:
	case FLIS_BLEND_COLOR:
	case FLIS_BLEND_LUMINOSITY:
	case FLIS_BLEND_CHROMA:
		return TRUE;
	default:
		return FALSE;
	}
}

gboolean nde_compositing_validate(const char *op_id, const char *params, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (!nde_compositing_is_op(op_id)) {
		*err = g_strdup_printf(_("'%s' is not a layer-property operation"),
		                       op_id ? op_id : "?");
		return FALSE;
	}
	GHashTable *kv = nde_kv_parse(params);
	gboolean ok = FALSE;
	if (!g_strcmp0(op_id, OP_SET_OPACITY)) {
		float v;
		if (!nde_kv_get_float(kv, "opacity", &v))
			*err = g_strdup(_("no opacity value"));
		else if (!(v >= 0.f && v <= 1.f))   /* rejects NaN too */
			*err = g_strdup_printf(_("opacity %g is outside the range 0 to 1"), v);
		else
			ok = TRUE;
	} else if (!g_strcmp0(op_id, OP_SET_BLEND)) {
		gint64 v;
		if (!nde_kv_get_int(kv, "blend", &v))
			*err = g_strdup(_("no blend mode value"));
		else if (!nde_compositing_blend_valid(v))
			*err = g_strdup_printf(_("%" G_GINT64_FORMAT " is not a valid layer blend mode"), v);
		else
			ok = TRUE;
	} else {
		gboolean v;
		if (!nde_kv_get_bool(kv, "visible", &v))
			*err = g_strdup(_("no visibility value"));
		else
			ok = TRUE;
	}
	g_hash_table_unref(kv);
	return ok;
}

gboolean nde_compositing_recompute(gint item_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	if (item_id < 0)
		return TRUE;   /* plain single image: nothing composites it */

	flis_layer_t *lay = flis_layer_get_by_id(item_id);
	if (!lay) {
		*err = g_strdup(_("the record's target layer no longer exists"));
		return FALSE;
	}

	/* The log is authoritative: fold from the creation defaults through every
	 * record, so the result is exactly what the history describes (header). */
	gfloat            opacity = COMP_DEFAULT_OPACITY;
	flis_blend_mode_t blend   = COMP_DEFAULT_BLEND;
	gboolean          visible = COMP_DEFAULT_VISIBLE;

	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->target_item_id != item_id)
			continue;

		/* Reset points: after these the layer provably sits at the defaults,
		 * so everything recorded earlier is superseded. */
		if (!g_strcmp0(rec->op_id, "document.flatten") ||
		    !g_strcmp0(rec->op_id, "layer.add")) {
			opacity = COMP_DEFAULT_OPACITY;
			blend   = COMP_DEFAULT_BLEND;
			visible = COMP_DEFAULT_VISIBLE;
			continue;
		}
		if (!nde_compositing_is_op(rec->op_id))
			continue;

		/* A record whose params are unreadable leaves the property at its
		 * running value — the alternative (bailing out) would strand the
		 * layer in a state no part of the log describes. */
		GHashTable *kv = nde_kv_parse(rec->params);
		if (!g_strcmp0(rec->op_id, OP_SET_OPACITY)) {
			float v;
			if (nde_kv_get_float(kv, "opacity", &v) && v >= 0.f && v <= 1.f)
				opacity = v;
		} else if (!g_strcmp0(rec->op_id, OP_SET_BLEND)) {
			gint64 v;
			if (nde_kv_get_int(kv, "blend", &v) && nde_compositing_blend_valid(v))
				blend = (flis_blend_mode_t)v;
		} else {
			gboolean v;
			if (nde_kv_get_bool(kv, "visible", &v))
				visible = v;
		}
		g_hash_table_unref(kv);
	}
	if (live)
		g_ptr_array_unref(live);

	flis_layer_set_opacity(lay, opacity);
	flis_layer_set_blend_mode(lay, blend);
	flis_layer_set_visible(lay, visible);

	gui_iface.flis_invalidate_composite();
	/* The Layers panel mirrors these three properties in its own widgets
	 * (opacity slider, blend combo, visibility toggle) — without this they
	 * would keep showing the pre-edit values. */
	gui_iface.flis_gui_update();
	return TRUE;
}
