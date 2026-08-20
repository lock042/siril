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
#include "core/nde/nde_compositing.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_op_class.h"
#include "core/gui_iface.h"
#include "io/image_format_flis.h"

#define OP_SET_OPACITY "layer.set_opacity"
#define OP_SET_BLEND   "layer.set_blend"
#define OP_SET_VISIBLE "layer.set_visible"
#define OP_SET_TINT    "layer.set_tint"

/* Compositing defaults, matching flis_layer_new() and the reset performed by
 * flis_flatten_all() on the surviving base layer. */
#define COMP_DEFAULT_OPACITY 1.0f
#define COMP_DEFAULT_BLEND   FLIS_BLEND_NORMAL
#define COMP_DEFAULT_VISIBLE TRUE

gboolean nde_compositing_is_op(const char *op_id) {
	return nde_op_class_for(op_id)->family == NDE_OPC_COMPOSITING;
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
	} else if (!g_strcmp0(op_id, OP_SET_TINT)) {
		gboolean tinted;
		if (!nde_kv_get_bool(kv, "tinted", &tinted)) {
			*err = g_strdup(_("no tinted flag"));
		} else if (!tinted) {
			ok = TRUE;   /* clearing needs no colour */
		} else {
			double r, g, b;
			if (!nde_kv_get_double(kv, "r", &r) ||
			    !nde_kv_get_double(kv, "g", &g) ||
			    !nde_kv_get_double(kv, "b", &b))
				*err = g_strdup(_("no tint colour components"));
			else if (!(r >= 0. && r <= 1. && g >= 0. && g <= 1. && b >= 0. && b <= 1.))
				*err = g_strdup(_("tint components must be in the range 0 to 1"));
			else
				ok = TRUE;
		}
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

/* TRUE when the live history contains any tint record for @item_id — the
 * gate that keeps the fold from clearing tints on documents recorded before
 * tint capture existed.  @upto_record_id bounds the search to records
 * strictly before that log position (0 = the whole live log). */
gboolean nde_compositing_has_tint_record_upto(gint item_id, gint64 upto_record_id) {
	gboolean found = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; !found && live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (upto_record_id && rec->record_id == upto_record_id)
			break;
		found = rec->target_item_id == item_id &&
		        !g_strcmp0(rec->op_id, OP_SET_TINT);
	}
	if (live)
		g_ptr_array_unref(live);
	return found;
}

gboolean nde_compositing_has_tint_record(gint item_id) {
	return nde_compositing_has_tint_record_upto(item_id, 0);
}

/* The fold itself, shared by the live-layer recompute, the retained-input
 * refresh (nde_composite_refresh_input_state) and the joint-record replay
 * (nde_joint.c, which needs the state AS OF a record's log position): walk
 * the live log and derive the compositing state the history describes for
 * @item_id.  @upto_record_id stops the fold just BEFORE that record
 * (exclusive; 0 folds the whole live log) — record ids are stable across
 * reorder, so the bound survives everything a position would not. */
void nde_compositing_fold_upto(gint item_id, gint64 upto_record_id,
                               gfloat *out_opacity,
                               gint *out_blend, gboolean *out_visible,
                               gboolean *out_tinted, double *out_tint /* [3] */) {
	gfloat            opacity = COMP_DEFAULT_OPACITY;
	flis_blend_mode_t blend   = COMP_DEFAULT_BLEND;
	gboolean          visible = COMP_DEFAULT_VISIBLE;
	gboolean          tinted  = FALSE;   /* creation default: no tint */
	double            tint[3] = { 1.0, 1.0, 1.0 };

	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (upto_record_id && rec->record_id == upto_record_id)
			break;
		if (rec->target_item_id != item_id)
			continue;

		/* Reset points: after these the layer provably sits at the defaults,
		 * so everything recorded earlier is superseded. */
		if (nde_op_class_for(rec->op_id)->traits & NDE_OPT_COMPOSITING_RESET) {
			opacity = COMP_DEFAULT_OPACITY;
			blend   = COMP_DEFAULT_BLEND;
			visible = COMP_DEFAULT_VISIBLE;
			tinted  = FALSE;
			tint[0] = tint[1] = tint[2] = 1.0;
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
		} else if (!g_strcmp0(rec->op_id, OP_SET_TINT)) {
			gboolean t;
			double r, g, b;
			if (nde_kv_get_bool(kv, "tinted", &t)) {
				if (!t) {
					tinted = FALSE;
				} else if (nde_kv_get_double(kv, "r", &r) &&
				           nde_kv_get_double(kv, "g", &g) &&
				           nde_kv_get_double(kv, "b", &b) &&
				           r >= 0. && r <= 1. && g >= 0. && g <= 1. &&
				           b >= 0. && b <= 1.) {
					tinted = TRUE;
					tint[0] = r; tint[1] = g; tint[2] = b;
				}
			}
		} else {
			gboolean v;
			if (nde_kv_get_bool(kv, "visible", &v))
				visible = v;
		}
		g_hash_table_unref(kv);
	}
	if (live)
		g_ptr_array_unref(live);
	if (out_opacity) *out_opacity = opacity;
	if (out_blend)   *out_blend   = (gint)blend;
	if (out_visible) *out_visible = visible;
	if (out_tinted)  *out_tinted  = tinted;
	if (out_tint) {
		out_tint[0] = tint[0];
		out_tint[1] = tint[1];
		out_tint[2] = tint[2];
	}
}

void nde_compositing_fold(gint item_id, gfloat *out_opacity,
                          gint *out_blend, gboolean *out_visible,
                          gboolean *out_tinted, double *out_tint /* [3] */) {
	nde_compositing_fold_upto(item_id, 0, out_opacity, out_blend, out_visible,
	                          out_tinted, out_tint);
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
	gfloat   opacity;
	gint     blend;
	gboolean visible;
	gboolean tinted;
	double   tint[3];
	nde_compositing_fold(item_id, &opacity, &blend, &visible, &tinted, tint);

	flis_layer_set_opacity(lay, opacity);
	flis_layer_set_blend_mode(lay, (flis_blend_mode_t)blend);
	flis_layer_set_visible(lay, visible);
	/* Tint joined the recorded set after blend/opacity/visibility: only let
	 * the fold drive it when the history actually says something about it,
	 * so an amend on an older document cannot clear a tint the log never
	 * described (nde_compositing_has_tint_record). */
	if (nde_compositing_has_tint_record(item_id)) {
		if (tinted)
			flis_layer_set_tint(lay, tint[0], tint[1], tint[2]);
		else
			flis_layer_clear_tint(lay);
	}

	gui_iface.flis_invalidate_composite();
	/* The Layers panel mirrors these three properties in its own widgets
	 * (opacity slider, blend combo, visibility toggle) — without this they
	 * would keep showing the pre-edit values. */
	gui_iface.flis_gui_update();
	return TRUE;
}
