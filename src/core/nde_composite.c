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

/* See nde_composite.h for what this node is and why its state is recorded. */

#include "core/siril.h"
#include "core/nde_history.h"
#include "core/nde_composite.h"
#include "io/image_format_flis.h"
#include "io/flis_compose.h"

void nde_composite_input_clear(nde_composite_input *in) {
	if (!in)
		return;
	g_free(in->name);
	in->name = NULL;
}

gboolean nde_composite_is_op(const char *op_id) {
	return !g_strcmp0(op_id, "layer.merge_down");
}

/* The keys the decoder demands.  Their presence is the format version: a
 * record written before this module existed carries top_item and top_name
 * only, and answers FALSE to every replayability question rather than being
 * re-run from parameters nobody stored. */
gchar *nde_composite_params_encode(const flis_layer_t *base,
                                   const flis_layer_t *top) {
	g_return_val_if_fail(base != NULL && top != NULL, NULL);
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "base_item", base->item_id);
	nde_kv_add_int(kv, "base_x", base->position_x);
	nde_kv_add_int(kv, "base_y", base->position_y);
	nde_kv_add_bool(kv, "base_tinted", base->has_tint);
	nde_kv_add_double(kv, "base_tint_r", base->layer_tint.r);
	nde_kv_add_double(kv, "base_tint_g", base->layer_tint.g);
	nde_kv_add_double(kv, "base_tint_b", base->layer_tint.b);
	nde_kv_add_int(kv, "top_item", top->item_id);
	nde_kv_add_str(kv, "top_name", top->layer_name ? top->layer_name : "");
	nde_kv_add_int(kv, "top_blend", (gint64)top->blend_mode);
	nde_kv_add_double(kv, "top_opacity", (double)top->opacity);
	nde_kv_add_int(kv, "top_x", top->position_x);
	nde_kv_add_int(kv, "top_y", top->position_y);
	nde_kv_add_bool(kv, "top_visible", top->visible);
	nde_kv_add_bool(kv, "top_tinted", top->has_tint);
	nde_kv_add_double(kv, "top_tint_r", top->layer_tint.r);
	nde_kv_add_double(kv, "top_tint_g", top->layer_tint.g);
	nde_kv_add_double(kv, "top_tint_b", top->layer_tint.b);
	/* A masked input is recorded but NOT re-runnable yet: the composite would
	 * have to rebuild a layer mask that died with the layer, which is a
	 * different resolver from the one image chains use.  Recording the fact
	 * lets the chain builder say so instead of quietly compositing without
	 * the mask and calling the result a reproduction. */
	nde_kv_add_bool(kv, "top_masked",
	                (top->lmask && top->lmask_active) || (top->fit && top->fit->mask));
	return nde_kv_end(kv);
}

gboolean nde_composite_params_decode(const char *params,
                                     nde_composite_input *base_out,
                                     nde_composite_input *top_out) {
	if (base_out)
		memset(base_out, 0, sizeof(*base_out));
	if (top_out)
		memset(top_out, 0, sizeof(*top_out));
	if (!params || !*params)
		return FALSE;
	GHashTable *kv = nde_kv_parse(params);
	nde_composite_input base = { 0 }, top = { 0 };
	gint64 bitem = 0, bx = 0, by = 0, titem = 0, blend = 0, tx = 0, ty = 0;
	gboolean ok = nde_kv_get_int(kv, "base_item", &bitem)
	           && nde_kv_get_int(kv, "base_x", &bx)
	           && nde_kv_get_int(kv, "base_y", &by)
	           && nde_kv_get_bool(kv, "base_tinted", &base.has_tint)
	           && nde_kv_get_double(kv, "base_tint_r", &base.tint_r)
	           && nde_kv_get_double(kv, "base_tint_g", &base.tint_g)
	           && nde_kv_get_double(kv, "base_tint_b", &base.tint_b)
	           && nde_kv_get_int(kv, "top_item", &titem)
	           && nde_kv_get_int(kv, "top_blend", &blend)
	           && nde_kv_get_double(kv, "top_opacity", &top.opacity)
	           && nde_kv_get_int(kv, "top_x", &tx)
	           && nde_kv_get_int(kv, "top_y", &ty)
	           && nde_kv_get_bool(kv, "top_visible", &top.visible)
	           && nde_kv_get_bool(kv, "top_tinted", &top.has_tint)
	           && nde_kv_get_double(kv, "top_tint_r", &top.tint_r)
	           && nde_kv_get_double(kv, "top_tint_g", &top.tint_g)
	           && nde_kv_get_double(kv, "top_tint_b", &top.tint_b);
	if (ok) {
		const char *name = nde_kv_get_str(kv, "top_name");
		base.item_id    = (gint)bitem;
		base.position_x = (gint)bx;
		base.position_y = (gint)by;
		base.visible    = TRUE;
		base.opacity    = 1.0;
		top.item_id     = (gint)titem;
		top.blend_mode  = (gint)blend;
		top.position_x  = (gint)tx;
		top.position_y  = (gint)ty;
		top.name        = (name && *name) ? g_strdup(name) : NULL;
		if (base_out)
			*base_out = base;
		if (top_out)
			*top_out = top;
		else
			nde_composite_input_clear(&top);
	}
	g_hash_table_unref(kv);
	return ok;
}

/* TRUE when the record says its overlay carried a mask.  Separate from the
 * decode so the caller can name that specific limitation. */
static gboolean overlay_was_masked(const char *params) {
	if (!params)
		return FALSE;
	GHashTable *kv = nde_kv_parse(params);
	gboolean masked = FALSE;
	nde_kv_get_bool(kv, "top_masked", &masked);
	g_hash_table_unref(kv);
	return masked;
}

gboolean nde_composite_record_replayable(const nde_record *rec) {
	if (!rec || !nde_composite_is_op(rec->op_id))
		return FALSE;
	if (!nde_record_input(rec, NDE_COMPOSITE_ROLE_BASE)
	    || !nde_record_input(rec, NDE_COMPOSITE_ROLE_OVERLAY))
		return FALSE;
	if (overlay_was_masked(rec->params))
		return FALSE;
	nde_composite_input in;
	gboolean ok = nde_composite_params_decode(rec->params, NULL, &in);
	nde_composite_input_clear(&in);
	return ok;
}

fits *nde_composite_render(const fits *base, const nde_composite_input *bs,
                           gint base_x, gint base_y,
                           const fits *overlay, const nde_composite_input *ov,
                           gchar **err) {
	if (!base || !overlay || !bs || !ov) {
		if (err)
			*err = g_strdup(_("the composite is missing an input"));
		return NULL;
	}
	/* Synthetic layers: the compositor's input is a flis_layer_t, and after a
	 * merge neither input is one any more.  Everything it reads is either
	 * resolved pixels or recorded state, so borrowing the struct is honest —
	 * these are stack copies that own nothing and are never linked into the
	 * document.
	 *
	 * The base is painted raw by the merge variant (first_raw), so its blend
	 * mode, opacity and visibility are not read: they are RETAINED on the
	 * surviving layer and baking them here would apply them twice.  Its tint
	 * IS read, because the merge bakes it. */
	flis_layer_t lb = { 0 }, lt = { 0 };
	lb.fit = (fits *)base;
	lb.blend_mode = FLIS_BLEND_NORMAL;
	lb.opacity = 1.f;
	lb.visible = TRUE;
	lb.has_tint = bs->has_tint;
	lb.layer_tint = (flis_tint_t){ bs->tint_r, bs->tint_g, bs->tint_b };
	lb.position_x = base_x;
	lb.position_y = base_y;

	lt.fit = (fits *)overlay;
	lt.blend_mode = (flis_blend_mode_t)ov->blend_mode;
	lt.opacity = (gfloat)ov->opacity;
	lt.visible = ov->visible;
	lt.has_tint = ov->has_tint;
	lt.layer_tint = (flis_tint_t){ ov->tint_r, ov->tint_g, ov->tint_b };
	lt.position_x = ov->position_x;
	lt.position_y = ov->position_y;

	GSList *pair = g_slist_append(NULL, &lb);
	pair = g_slist_append(pair, &lt);
	fits *out = flis_render_layers_merge(pair);
	g_slist_free(pair);
	if (!out && err)
		*err = g_strdup(_("the composite could not be rendered"));
	return out;
}
