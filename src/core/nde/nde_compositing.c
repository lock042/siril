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

#include "core/gui_iface.h"
#include "core/nde/nde_compositing.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_op_class.h"
#include "core/siril.h"
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

/* ---- where the fold starts ----------------------------------------------- */

typedef struct {
	gfloat            opacity;
	flis_blend_mode_t blend;
	gboolean          visible;
	gboolean          tinted;
	double            tint[3];
} comp_state;

static void comp_state_defaults(comp_state *st) {
	st->opacity = COMP_DEFAULT_OPACITY;
	st->blend   = COMP_DEFAULT_BLEND;
	st->visible = COMP_DEFAULT_VISIBLE;
	st->tinted  = FALSE;   /* creation default: no tint */
	st->tint[0] = st->tint[1] = st->tint[2] = 1.0;
}

/* The starting state a compositing record carries (base_* keys): what its
 * layer looked like before the first such record.  All or nothing. */
static gboolean comp_base_read(GHashTable *kv, comp_state *st) {
	float o;
	gint64 bl;
	gboolean vis, tinted;
	double r = 1.0, g = 1.0, b = 1.0;
	if (!nde_kv_get_float(kv, "base_opacity", &o) || o < 0.f || o > 1.f ||
	    !nde_kv_get_int(kv, "base_blend", &bl) || !nde_compositing_blend_valid(bl) ||
	    !nde_kv_get_bool(kv, "base_visible", &vis) ||
	    !nde_kv_get_bool(kv, "base_tinted", &tinted))
		return FALSE;
	if (tinted && (!nde_kv_get_double(kv, "base_r", &r) ||
	               !nde_kv_get_double(kv, "base_g", &g) ||
	               !nde_kv_get_double(kv, "base_b", &b) ||
	               r < 0. || r > 1. || g < 0. || g > 1. || b < 0. || b > 1.))
		return FALSE;
	st->opacity = o;
	st->blend   = (flis_blend_mode_t)bl;
	st->visible = vis;
	st->tinted  = tinted;
	st->tint[0] = r; st->tint[1] = g; st->tint[2] = b;
	return TRUE;
}

static void comp_base_write(GString *kv, const comp_state *st) {
	nde_kv_add_float(kv, "base_opacity", st->opacity);
	nde_kv_add_int(kv, "base_blend", (gint64)st->blend);
	nde_kv_add_bool(kv, "base_visible", st->visible);
	nde_kv_add_bool(kv, "base_tinted", st->tinted);
	if (st->tinted) {
		nde_kv_add_double(kv, "base_r", st->tint[0]);
		nde_kv_add_double(kv, "base_g", st->tint[1]);
		nde_kv_add_double(kv, "base_b", st->tint[2]);
	}
}

/* The state @item_id was in before the log's first word about it.  The
 * new-layer defaults are right only for a layer the log saw being created;
 * one that arrived with its own blend mode, opacity or tint — every layer of
 * a file opened without a history, every duplicate — starts from those, and
 * folding it from the defaults made deleting its one opacity step reset its
 * blend mode as well.  In order of authority:
 *   - the base a compositing record of the item carries (it travels with the
 *     file, so it is the answer in a later session too);
 *   - what the live layer was marked with when it joined the document;
 *   - the defaults.
 * TRUE when the answer is one of the first two, i.e. actually known. */
static gboolean comp_start(GPtrArray *live, gint item_id, comp_state *st) {
	gboolean known = FALSE;
	comp_state_defaults(st);
	const flis_layer_t *lay = item_id >= 0 ? flis_layer_get_by_id(item_id) : NULL;
	if (lay && lay->comp_base_valid) {
		st->opacity = lay->comp_base_opacity;
		st->blend   = lay->comp_base_blend;
		st->visible = lay->comp_base_visible;
		st->tinted  = lay->comp_base_tinted;
		st->tint[0] = lay->comp_base_tint.r;
		st->tint[1] = lay->comp_base_tint.g;
		st->tint[2] = lay->comp_base_tint.b;
		known = TRUE;
	}
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (rec->target_item_id != item_id)
			continue;
		/* Past a reset the fold starts over from the defaults anyway. */
		if (nde_op_class_for(rec->op_id)->traits & NDE_OPT_COMPOSITING_RESET)
			break;
		if (!nde_compositing_is_op(rec->op_id))
			continue;
		GHashTable *kv = nde_kv_parse(rec->params);
		gboolean has = comp_base_read(kv, st);
		g_hash_table_unref(kv);
		if (has) {
			known = TRUE;
			break;
		}
	}
	return known;
}

gchar *nde_compositing_params_with_base(gchar *params, gint item_id) {
	comp_state st;
	GPtrArray *live = nde_history_snapshot(NULL);
	gboolean known = comp_start(live, item_id, &st);
	if (live)
		g_ptr_array_unref(live);
	if (!known)
		return params;
	GString *kv = g_string_new(params ? params : "");
	g_free(params);
	comp_base_write(kv, &st);
	return g_string_free(kv, FALSE);
}

void nde_compositing_adopt_recorded_bases(void) {
	if (!com.uniq)
		return;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = l->data;
		comp_state st;
		if (!comp_start(live, lay->item_id, &st))
			continue;
		lay->comp_base_opacity = st.opacity;
		lay->comp_base_blend   = st.blend;
		lay->comp_base_visible = st.visible;
		lay->comp_base_tinted  = st.tinted;
		lay->comp_base_tint.r  = st.tint[0];
		lay->comp_base_tint.g  = st.tint[1];
		lay->comp_base_tint.b  = st.tint[2];
		lay->comp_base_valid   = TRUE;
	}
	if (live)
		g_ptr_array_unref(live);
}

gchar *nde_compositing_carry_base(const char *old_params, const char *new_params) {
	comp_state st;
	GHashTable *okv = nde_kv_parse(old_params);
	GHashTable *nkv = nde_kv_parse(new_params);
	gboolean carry = comp_base_read(okv, &st) &&
	                 !g_hash_table_contains(nkv, "base_opacity");
	g_hash_table_unref(okv);
	g_hash_table_unref(nkv);
	if (!carry)
		return g_strdup(new_params);
	GString *kv = g_string_new(new_params ? new_params : "");
	comp_base_write(kv, &st);
	return g_string_free(kv, FALSE);
}

/* The fold itself, shared by the live-layer recompute, the retained-input
 * refresh (nde_composite_refresh_input_state) and the joint-record replay
 * (nde_joint.c, which needs the state AS OF a record's log position): walk
 * the live log and derive the compositing state the history describes for
 * @item_id, from the state it started in (comp_start).  @upto_record_id stops
 * the fold just BEFORE that record (exclusive; 0 folds the whole live log) —
 * record ids are stable across reorder, so the bound survives everything a
 * position would not.  Returns whether the starting state was known. */
static gboolean comp_fold(gint item_id, gint64 upto_record_id, comp_state *st) {
	GPtrArray *live = nde_history_snapshot(NULL);
	gboolean known = comp_start(live, item_id, st);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (upto_record_id && rec->record_id == upto_record_id)
			break;
		if (rec->target_item_id != item_id)
			continue;

		/* Reset points: after these the layer provably sits at the defaults,
		 * so everything recorded earlier is superseded. */
		if (nde_op_class_for(rec->op_id)->traits & NDE_OPT_COMPOSITING_RESET) {
			comp_state_defaults(st);
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
				st->opacity = v;
		} else if (!g_strcmp0(rec->op_id, OP_SET_BLEND)) {
			gint64 v;
			if (nde_kv_get_int(kv, "blend", &v) && nde_compositing_blend_valid(v))
				st->blend = (flis_blend_mode_t)v;
		} else if (!g_strcmp0(rec->op_id, OP_SET_TINT)) {
			gboolean t;
			double r, g, b;
			if (nde_kv_get_bool(kv, "tinted", &t)) {
				if (!t) {
					st->tinted = FALSE;
				} else if (nde_kv_get_double(kv, "r", &r) &&
				           nde_kv_get_double(kv, "g", &g) &&
				           nde_kv_get_double(kv, "b", &b) &&
				           r >= 0. && r <= 1. && g >= 0. && g <= 1. &&
				           b >= 0. && b <= 1.) {
					st->tinted = TRUE;
					st->tint[0] = r; st->tint[1] = g; st->tint[2] = b;
				}
			}
		} else {
			gboolean v;
			if (nde_kv_get_bool(kv, "visible", &v))
				st->visible = v;
		}
		g_hash_table_unref(kv);
	}
	if (live)
		g_ptr_array_unref(live);
	return known;
}

void nde_compositing_fold_upto(gint item_id, gint64 upto_record_id,
                               gfloat *out_opacity,
                               gint *out_blend, gboolean *out_visible,
                               gboolean *out_tinted, double *out_tint /* [3] */) {
	comp_state st;
	comp_fold(item_id, upto_record_id, &st);
	if (out_opacity) *out_opacity = st.opacity;
	if (out_blend)   *out_blend   = (gint)st.blend;
	if (out_visible) *out_visible = st.visible;
	if (out_tinted)  *out_tinted  = st.tinted;
	if (out_tint) {
		out_tint[0] = st.tint[0];
		out_tint[1] = st.tint[1];
		out_tint[2] = st.tint[2];
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

	/* The log is authoritative: fold from the state the layer started in
	 * through every record, so the result is exactly what the history
	 * describes (header). */
	comp_state st;
	gboolean known = comp_fold(item_id, 0, &st);

	flis_layer_set_opacity(lay, st.opacity);
	flis_layer_set_blend_mode(lay, st.blend);
	flis_layer_set_visible(lay, st.visible);
	/* Tint joined the recorded set after blend/opacity/visibility: only let
	 * the fold drive it when the history actually says something about it
	 * or the starting state is known, so an amend on an older document
	 * cannot clear a tint the log never described
	 * (nde_compositing_has_tint_record). */
	if (known || nde_compositing_has_tint_record(item_id)) {
		if (st.tinted)
			flis_layer_set_tint(lay, st.tint[0], st.tint[1], st.tint[2]);
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
