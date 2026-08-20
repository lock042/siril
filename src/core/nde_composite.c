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
#include "core/proto.h"
#include "core/masks.h"
#include "core/nde_history.h"
#include "core/nde_op_class.h"
#include "core/nde_checkpoint.h"
#include "core/nde_compositing.h"
#include "core/nde_composite.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/flis_compose.h"

/* Indexed key names: "i3_opacity", "g0_blend".  One buffer per call site. */
#define KEYBUF 32
static const char *ikey(char *buf, const char *prefix, guint i, const char *field) {
	g_snprintf(buf, KEYBUF, "%s%u_%s", prefix, i, field);
	return buf;
}

void nde_composite_state_free(nde_composite_state *st) {
	if (!st)
		return;
	if (st->inputs) {
		for (guint i = 0; i < st->inputs->len; i++)
			g_free(g_array_index(st->inputs, nde_composite_input, i).name);
		g_array_free(st->inputs, TRUE);
	}
	if (st->groups) {
		for (guint i = 0; i < st->groups->len; i++)
			g_free(g_array_index(st->groups, nde_composite_group, i).name);
		g_array_free(st->groups, TRUE);
	}
	g_free(st);
}

static nde_composite_state *state_new(void) {
	nde_composite_state *st = g_new0(nde_composite_state, 1);
	st->inputs = g_array_new(FALSE, TRUE, sizeof(nde_composite_input));
	st->groups = g_array_new(FALSE, TRUE, sizeof(nde_composite_group));
	return st;
}

gboolean nde_composite_is_op(const char *op_id) {
	return nde_op_class_for(op_id)->family == NDE_OPC_COMPOSITE;
}

/* ---- encoding ----------------------------------------------------------- */

/* Only the LAYER mask: the compositor reads lay->lmask and never fit->mask,
 * which restricts operations rather than compositing. */
static gboolean layer_is_masked(const flis_layer_t *lay) {
	return lay->lmask && lay->lmask_active && lay->lmask->data;
}

/* A layer mask as the mono fits the checkpoint store keeps.  layermask_t and
 * mask_t are the same pixels in different wrappers (commit_mask_value does the
 * reverse); mask_to_fits needs a host only for its dimensions and format. */
static fits *lmask_to_fits(const flis_layer_t *lay) {
	fits host = { 0 };
	if (copyfits((fits *)lay->fit, &host, CP_FORMAT, -1))
		return NULL;
	mask_t m = { .bitpix = lay->lmask->bitpix, .data = lay->lmask->data };
	host.mask = &m;
	fits *out = mask_to_fits(&host);
	host.mask = NULL;   /* borrowed from the layer — not ours to free */
	clearfits(&host);
	return out;
}

/* The document state @layers were composited against.  Only the groups those
 * layers actually belong to are recorded — the rest did not participate. */
static void collect_groups(nde_composite_state *st, GSList *layers) {
	for (GSList *l = layers; l; l = l->next) {
		const flis_layer_t *lay = l->data;
		if (!lay || !lay->group_id)
			continue;
		gboolean seen = FALSE;
		for (guint i = 0; i < st->groups->len && !seen; i++)
			seen = g_array_index(st->groups, nde_composite_group, i).item_id == lay->group_id;
		if (seen)
			continue;
		flis_group_t *grp = flis_group_get_by_id(lay->group_id);
		if (!grp)
			continue;
		nde_composite_group g = {
			.item_id    = grp->item_id,
			.name       = (grp->name && *grp->name) ? g_strdup(grp->name) : NULL,
			.blend_mode = (gint)grp->blend_mode,
			.opacity    = (gdouble)grp->opacity,
			.visible    = grp->visible,
		};
		g_array_append_val(st->groups, g);
	}
}

static nde_composite_state *state_from_layers(GSList *layers, gboolean raw_first) {
	nde_composite_state *st = state_new();
	st->raw_first = raw_first;
	if (com.uniq) {
		st->canvas_w = com.uniq->canvas_w;
		st->canvas_h = com.uniq->canvas_h;
		st->bg_r = com.uniq->canvas_bg_r;
		st->bg_g = com.uniq->canvas_bg_g;
		st->bg_b = com.uniq->canvas_bg_b;
	}
	guint idx = 0;
	for (GSList *l = layers; l; l = l->next, idx++) {
		flis_layer_t *lay = l->data;
		if (!lay)
			continue;
		/* A mask on the raw-painted first input is not applied — that is what
		 * raw means — so it does not participate in the render (was_masked
		 * stays FALSE).  Its item id is still recorded and its copy stored so
		 * that UNDOING the composite can put the pre-merge mask back
		 * (nde_composite_undo_execute); records written before this carry no
		 * id for it and restore that layer maskless, as they always did. */
		gboolean has_mask = layer_is_masked(lay);
		gboolean masked = has_mask && !(idx == 0 && raw_first);
		nde_composite_input in = {
			.item_id    = lay->item_id,
			.name       = g_strdup(lay->layer_name ? lay->layer_name : ""),
			.blend_mode = (gint)lay->blend_mode,
			.opacity    = (gdouble)lay->opacity,
			.position_x = lay->position_x,
			.position_y = lay->position_y,
			.visible    = lay->visible,
			.has_tint   = lay->has_tint,
			.tint_r     = lay->layer_tint.r,
			.tint_g     = lay->layer_tint.g,
			.tint_b     = lay->layer_tint.b,
			.group_id   = lay->group_id,
			.layer_order = lay->layer_order,
			.was_masked = masked,
			/* flis_layer_lmask_id hands out an id on first use, so ask for one
			 * only when the mask is going to be stored under it. */
			.mask_item_id = has_mask ? flis_layer_lmask_id(lay) : 0,
		};
		g_array_append_val(st->inputs, in);
	}
	collect_groups(st, layers);
	return st;
}

static gchar *state_encode(const nde_composite_state *st) {
	char k[KEYBUF];
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "v", 2);
	nde_kv_add_bool(kv, "raw_first", st->raw_first);
	nde_kv_add_int(kv, "canvas_w", st->canvas_w);
	nde_kv_add_int(kv, "canvas_h", st->canvas_h);
	nde_kv_add_double(kv, "bg_r", st->bg_r);
	nde_kv_add_double(kv, "bg_g", st->bg_g);
	nde_kv_add_double(kv, "bg_b", st->bg_b);
	nde_kv_add_int(kv, "n", st->inputs->len);
	for (guint i = 0; i < st->inputs->len; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		nde_kv_add_int(kv, ikey(k, "i", i, "item"), in->item_id);
		nde_kv_add_str(kv, ikey(k, "i", i, "name"), in->name ? in->name : "");
		nde_kv_add_int(kv, ikey(k, "i", i, "blend"), in->blend_mode);
		nde_kv_add_double(kv, ikey(k, "i", i, "opacity"), in->opacity);
		nde_kv_add_int(kv, ikey(k, "i", i, "x"), in->position_x);
		nde_kv_add_int(kv, ikey(k, "i", i, "y"), in->position_y);
		nde_kv_add_bool(kv, ikey(k, "i", i, "visible"), in->visible);
		nde_kv_add_bool(kv, ikey(k, "i", i, "tinted"), in->has_tint);
		nde_kv_add_double(kv, ikey(k, "i", i, "tint_r"), in->tint_r);
		nde_kv_add_double(kv, ikey(k, "i", i, "tint_g"), in->tint_g);
		nde_kv_add_double(kv, ikey(k, "i", i, "tint_b"), in->tint_b);
		nde_kv_add_int(kv, ikey(k, "i", i, "group"), in->group_id);
		nde_kv_add_int(kv, ikey(k, "i", i, "order"), in->layer_order);
		nde_kv_add_bool(kv, ikey(k, "i", i, "masked"), in->was_masked);
		nde_kv_add_int(kv, ikey(k, "i", i, "maskitem"), in->mask_item_id);
	}
	nde_kv_add_int(kv, "ng", st->groups->len);
	for (guint i = 0; i < st->groups->len; i++) {
		const nde_composite_group *g = &g_array_index(st->groups, nde_composite_group, i);
		nde_kv_add_int(kv, ikey(k, "g", i, "item"), g->item_id);
		nde_kv_add_str(kv, ikey(k, "g", i, "name"), g->name ? g->name : "");
		nde_kv_add_int(kv, ikey(k, "g", i, "blend"), g->blend_mode);
		nde_kv_add_double(kv, ikey(k, "g", i, "opacity"), g->opacity);
		nde_kv_add_bool(kv, ikey(k, "g", i, "visible"), g->visible);
	}
	return nde_kv_end(kv);
}

/* ---- decoding ----------------------------------------------------------- */

/* The first merge-down format (shipped before flatten joined this node): two
 * inputs, the base contributing only its offset and tint because the merge
 * retains the rest on the survivor.  Kept so that documents saved with it stay
 * replayable rather than silently reverting to blockers. */
static nde_composite_state *parse_legacy_merge(GHashTable *kv) {
	nde_composite_state *st = state_new();
	st->raw_first = TRUE;
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
	if (!ok) {
		nde_composite_state_free(st);
		return NULL;
	}
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
	nde_kv_get_bool(kv, "top_masked", &top.was_masked);
	g_array_append_val(st->inputs, base);
	g_array_append_val(st->inputs, top);
	return st;
}

nde_composite_state *nde_composite_state_parse(const char *params) {
	if (!params || !*params)
		return NULL;
	GHashTable *kv = nde_kv_parse(params);
	gint64 version = 0;
	if (!nde_kv_get_int(kv, "v", &version)) {
		nde_composite_state *legacy = parse_legacy_merge(kv);
		g_hash_table_unref(kv);
		return legacy;
	}

	char k[KEYBUF];
	nde_composite_state *st = state_new();
	gint64 n = 0, ng = 0, cw = 0, ch = 0;
	gboolean ok = nde_kv_get_bool(kv, "raw_first", &st->raw_first)
	           && nde_kv_get_int(kv, "canvas_w", &cw)
	           && nde_kv_get_int(kv, "canvas_h", &ch)
	           && nde_kv_get_double(kv, "bg_r", &st->bg_r)
	           && nde_kv_get_double(kv, "bg_g", &st->bg_g)
	           && nde_kv_get_double(kv, "bg_b", &st->bg_b)
	           && nde_kv_get_int(kv, "n", &n)
	           && nde_kv_get_int(kv, "ng", &ng)
	           && n > 0;
	st->canvas_w = (guint)MAX(cw, 0);
	st->canvas_h = (guint)MAX(ch, 0);
	for (gint64 i = 0; ok && i < n; i++) {
		nde_composite_input in = { 0 };
		gint64 item = 0, blend = 0, x = 0, y = 0, group = 0;
		ok = nde_kv_get_int(kv, ikey(k, "i", (guint)i, "item"), &item)
		  && nde_kv_get_int(kv, ikey(k, "i", (guint)i, "blend"), &blend)
		  && nde_kv_get_double(kv, ikey(k, "i", (guint)i, "opacity"), &in.opacity)
		  && nde_kv_get_int(kv, ikey(k, "i", (guint)i, "x"), &x)
		  && nde_kv_get_int(kv, ikey(k, "i", (guint)i, "y"), &y)
		  && nde_kv_get_bool(kv, ikey(k, "i", (guint)i, "visible"), &in.visible)
		  && nde_kv_get_bool(kv, ikey(k, "i", (guint)i, "tinted"), &in.has_tint)
		  && nde_kv_get_double(kv, ikey(k, "i", (guint)i, "tint_r"), &in.tint_r)
		  && nde_kv_get_double(kv, ikey(k, "i", (guint)i, "tint_g"), &in.tint_g)
		  && nde_kv_get_double(kv, ikey(k, "i", (guint)i, "tint_b"), &in.tint_b)
		  && nde_kv_get_int(kv, ikey(k, "i", (guint)i, "group"), &group)
		  && nde_kv_get_bool(kv, ikey(k, "i", (guint)i, "masked"), &in.was_masked);
		if (!ok)
			break;
		gint64 maskitem = 0;
		nde_kv_get_int(kv, ikey(k, "i", (guint)i, "maskitem"), &maskitem);
		in.mask_item_id = (gint)maskitem;
		/* Optional, like maskitem: 0 in records written before it existed. */
		gint64 order = 0;
		nde_kv_get_int(kv, ikey(k, "i", (guint)i, "order"), &order);
		in.layer_order = (gint)order;
		const char *name = nde_kv_get_str(kv, ikey(k, "i", (guint)i, "name"));
		in.item_id    = (gint)item;
		in.blend_mode = (gint)blend;
		in.position_x = (gint)x;
		in.position_y = (gint)y;
		in.group_id   = (gint)group;
		in.name       = (name && *name) ? g_strdup(name) : NULL;
		g_array_append_val(st->inputs, in);
	}
	for (gint64 i = 0; ok && i < ng; i++) {
		nde_composite_group g = { 0 };
		gint64 item = 0, blend = 0;
		ok = nde_kv_get_int(kv, ikey(k, "g", (guint)i, "item"), &item)
		  && nde_kv_get_int(kv, ikey(k, "g", (guint)i, "blend"), &blend)
		  && nde_kv_get_double(kv, ikey(k, "g", (guint)i, "opacity"), &g.opacity)
		  && nde_kv_get_bool(kv, ikey(k, "g", (guint)i, "visible"), &g.visible);
		if (!ok)
			break;
		/* Optional: records written before groups carried a name still decode,
		 * they just have nothing to label the group with. */
		const char *gname = nde_kv_get_str(kv, ikey(k, "g", (guint)i, "name"));
		g.item_id    = (gint)item;
		g.name       = (gname && *gname) ? g_strdup(gname) : NULL;
		g.blend_mode = (gint)blend;
		g_array_append_val(st->groups, g);
	}
	g_hash_table_unref(kv);
	if (!ok) {
		nde_composite_state_free(st);
		return NULL;
	}
	return st;
}

/* ---- validating an amend ------------------------------------------------ */

static gboolean in_unit_range(gdouble v) {
	return v >= 0.0 && v <= 1.0;   /* rejects NaN too */
}

gboolean nde_composite_validate(const char *old_params, const char *new_params,
                                gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	nde_composite_state *was = nde_composite_state_parse(old_params);
	nde_composite_state *now = nde_composite_state_parse(new_params);
	if (!now) {
		*err = g_strdup(_("the new compositing parameters could not be read"));
	} else if (!was) {
		/* Nothing to compare against: the record predates this format, and it
		 * is not replayable either, so there is nothing an edit could do. */
		*err = g_strdup(_("this step did not record what it composited"));
	} else if (was->raw_first != now->raw_first ||
	           was->canvas_w != now->canvas_w || was->canvas_h != now->canvas_h ||
	           was->inputs->len != now->inputs->len ||
	           was->groups->len != now->groups->len) {
		*err = g_strdup(_("only the compositing parameters of this step can be changed, "
		                  "not what it composited"));
	}
	for (guint i = 0; !*err && i < now->inputs->len; i++) {
		const nde_composite_input *a = &g_array_index(was->inputs, nde_composite_input, i);
		const nde_composite_input *b = &g_array_index(now->inputs, nde_composite_input, i);
		if (a->item_id != b->item_id || a->mask_item_id != b->mask_item_id ||
		    a->group_id != b->group_id || a->was_masked != b->was_masked)
			*err = g_strdup_printf(_("input %u names a different layer — a step cannot be "
			                         "rewired, only its compositing parameters changed"), i);
		else if (!in_unit_range(b->opacity))
			*err = g_strdup_printf(_("opacity %g for '%s' is outside the range 0 to 1"),
			                       b->opacity, b->name ? b->name : "?");
		else if (!nde_compositing_blend_valid(b->blend_mode))
			*err = g_strdup_printf(_("%d is not a valid layer blend mode"), b->blend_mode);
		else if (b->has_tint && !(in_unit_range(b->tint_r) && in_unit_range(b->tint_g)
		                          && in_unit_range(b->tint_b)))
			*err = g_strdup_printf(_("the tint for '%s' is outside the range 0 to 1"),
			                       b->name ? b->name : "?");
	}
	for (guint i = 0; !*err && i < now->groups->len; i++) {
		const nde_composite_group *a = &g_array_index(was->groups, nde_composite_group, i);
		const nde_composite_group *b = &g_array_index(now->groups, nde_composite_group, i);
		if (a->item_id != b->item_id)
			*err = g_strdup_printf(_("group %u names a different group"), i);
		else if (!in_unit_range(b->opacity))
			*err = g_strdup_printf(_("group opacity %g is outside the range 0 to 1"),
			                       b->opacity);
		/* A group MAY be PASS_THROUGH, unlike a layer, so its blend mode is
		 * checked against the wider set. */
		else if (!nde_compositing_blend_valid(b->blend_mode) &&
		         b->blend_mode != FLIS_BLEND_PASS_THROUGH)
			*err = g_strdup_printf(_("%d is not a valid group blend mode"), b->blend_mode);
	}
	nde_composite_state_free(was);
	nde_composite_state_free(now);
	return *err == NULL;
}

gchar *nde_composite_params_drop_mask(const char *params, gint mask_item_id) {
	nde_composite_state *st = nde_composite_state_parse(params);
	if (!st)
		return NULL;
	gboolean hit = FALSE;
	for (guint i = 0; i < st->inputs->len; i++) {
		nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		if (!in->mask_item_id || in->mask_item_id != mask_item_id)
			continue;
		in->mask_item_id = 0;
		in->was_masked   = FALSE;
		hit = TRUE;
	}
	gchar *out = hit ? state_encode(st) : NULL;
	nde_composite_state_free(st);
	return out;
}

/* Re-derive a RETAINED input's recorded compositing state (blend / opacity /
 * visibility) from the item's amended history and patch it into every live
 * replayable composite record that pins the item.  This is what makes a
 * compositing amend on a flattened/merged-away layer effective: the layer no
 * longer exists to re-fold onto, but the composite's recorded copy of its
 * state is what the re-render reads — so the log stays authoritative through
 * the composite boundary.  Tint is left as captured (tint edits are pixel-
 * adjacent, not compositing records).  Returns FALSE with @err when no
 * consuming record could be found. */
gboolean nde_composite_refresh_input_state(gint item_id, gchar **err) {
	g_return_val_if_fail(err != NULL, FALSE);
	*err = NULL;
	gfloat   opacity;
	gint     blend;
	gboolean visible;
	gboolean tinted;
	double   tint[3];
	nde_compositing_fold(item_id, &opacity, &blend, &visible, &tinted, tint);
	const gboolean fold_tint = nde_compositing_has_tint_record(item_id);

	gboolean patched = FALSE;
	GPtrArray *live = nde_history_snapshot(NULL);
	for (guint i = 0; live && i < live->len; i++) {
		const nde_record *rec = g_ptr_array_index(live, i);
		if (!nde_composite_record_replayable(rec))
			continue;
		if (rec->target_item_id == item_id ||
		    !nde_record_input_by_item(rec, item_id))
			continue;
		nde_composite_state *st = nde_composite_state_parse(rec->params);
		if (!st)
			continue;
		for (guint j = 0; j < st->inputs->len; j++) {
			nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, j);
			if (in->item_id != item_id)
				continue;
			in->blend_mode = blend;
			in->opacity    = (gdouble)opacity;
			in->visible    = visible;
			/* Same late-arrival gate as the live recompute: only let the
			 * fold drive the tint when the history says something about
			 * it, so older records keep their captured tint. */
			if (fold_tint) {
				in->has_tint = tinted;
				in->tint_r = tint[0];
				in->tint_g = tint[1];
				in->tint_b = tint[2];
			}
		}
		gchar *blob = state_encode(st);
		nde_composite_state_free(st);
		gchar *aerr = NULL;
		if (!nde_history_amend(rec->record_id, blob, &aerr)) {
			g_free(blob);
			*err = aerr ? aerr : g_strdup(_("could not update the composite record"));
			if (live)
				g_ptr_array_unref(live);
			return FALSE;
		}
		g_free(blob);
		patched = TRUE;
	}
	if (live)
		g_ptr_array_unref(live);
	if (!patched)
		*err = g_strdup(_("no composite consumes this layer"));
	return patched;
}

gboolean nde_composite_record_replayable(const nde_record *rec) {
	if (!rec || !nde_composite_is_op(rec->op_id))
		return FALSE;
	nde_composite_state *st = nde_composite_state_parse(rec->params);
	if (!st)
		return FALSE;
	gboolean ok = TRUE;
	for (guint i = 0; i < st->inputs->len && ok; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		if (!in->visible)
			continue;   /* contributes nothing; never resolved */
		/* Every visible input needs a pin to resolve its pixels through. */
		if (!nde_record_input_by_item(rec, in->item_id))
			ok = FALSE;
		/* A masked one also needs its stored mask still to be there.  The
		 * store is a cache under a budget, so this is a live question, not a
		 * property of the record. */
		if (in->was_masked) {
			const nde_pin *mp = in->mask_item_id ?
					nde_record_input_by_item(rec, in->mask_item_id) : NULL;
			if (!mp || !nde_checkpoint_exists_at(mp->item_id, mp->record_id))
				ok = FALSE;
		}
	}
	nde_composite_state_free(st);
	return ok;
}

/* ---- rendering ---------------------------------------------------------- */

fits *nde_composite_render(const nde_composite_state *st,
                           fits *const *pixels, fits *const *masks,
                           gchar **err) {
	if (!st || !st->inputs->len || !pixels) {
		if (err)
			*err = g_strdup(_("the composite is missing an input"));
		return NULL;
	}
	/* Synthetic layers and groups: the compositor's inputs are flis_layer_t and
	 * flis_group_t, and after a composite none of them are one any more.
	 * Everything it reads is either resolved pixels or recorded state, so
	 * borrowing the structs is honest — these are stack copies that own nothing
	 * and are never linked into the document. */
	guint n = st->inputs->len;
	flis_layer_t *lays = g_new0(flis_layer_t, n);
	layermask_t  *lms  = g_new0(layermask_t, n);
	flis_group_t *grps = st->groups->len ? g_new0(flis_group_t, st->groups->len) : NULL;
	GSList *list = NULL, *glist = NULL;
	gchar *fail = NULL;
	for (guint i = 0; i < n; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		if (in->visible && !pixels[i]) {
			fail = g_strdup_printf(_("the pixels of input '%s' could not be resolved"),
			                       in->name ? in->name : "?");
			break;
		}
		if (in->visible && in->was_masked) {
			/* Same unwrapping as commit_mask_value, in the other direction: the
			 * stored mono image becomes the 8-bit alpha the compositor reads. */
			mask_t *m = masks && masks[i] ? fits_to_mask(masks[i]) : NULL;
			if (!m) {
				fail = g_strdup_printf(_("the layer mask of input '%s' is no longer stored"),
				                       in->name ? in->name : "?");
				break;
			}
			lms[i].w      = masks[i]->rx;
			lms[i].h      = masks[i]->ry;
			lms[i].bitpix = m->bitpix;
			lms[i].data   = m->data;   /* moved */
			free(m);
			lays[i].lmask        = &lms[i];
			lays[i].lmask_active = TRUE;
		}
		lays[i].fit        = pixels[i];
		lays[i].item_id    = in->item_id;
		lays[i].blend_mode = (flis_blend_mode_t)in->blend_mode;
		lays[i].opacity    = (gfloat)in->opacity;
		lays[i].visible    = in->visible;
		lays[i].has_tint   = in->has_tint;
		lays[i].layer_tint = (flis_tint_t){ in->tint_r, in->tint_g, in->tint_b };
		lays[i].position_x = in->position_x;
		lays[i].position_y = in->position_y;
		lays[i].group_id   = in->group_id;
		list = g_slist_append(list, &lays[i]);
	}
	for (guint i = 0; !fail && i < st->groups->len; i++) {
		const nde_composite_group *g = &g_array_index(st->groups, nde_composite_group, i);
		grps[i].item_id    = g->item_id;
		grps[i].blend_mode = (flis_blend_mode_t)g->blend_mode;
		grps[i].opacity    = (gfloat)g->opacity;
		grps[i].visible    = g->visible;
		glist = g_slist_append(glist, &grps[i]);
	}

	fits *out = NULL;
	if (!fail) {
		flis_render_ctx ctx = {
			.canvas_w = st->canvas_w,
			.canvas_h = st->canvas_h,
			.bg_r = st->bg_r, .bg_g = st->bg_g, .bg_b = st->bg_b,
			.groups = glist,
		};
		/* raw_first also selects the sub-composite contract, and for the same
		 * reason: merge-down composites two layers in isolation — no canvas
		 * background beneath them, no group pre-pass — while flatten renders the
		 * document exactly as the display does. */
		out = flis_render_layers_ctx(list, &ctx, st->raw_first, st->raw_first);
		if (!out)
			fail = g_strdup(_("the composite could not be rendered"));
	}
	g_slist_free(list);
	g_slist_free(glist);
	for (guint i = 0; i < n; i++)
		free(lms[i].data);
	g_free(lms);
	g_free(lays);
	g_free(grps);
	if (fail) {
		if (err)
			*err = fail;
		else
			g_free(fail);
	}
	return out;
}

/* ---- capture ------------------------------------------------------------ */

struct _nde_composite_capture {
	gchar        *params;
	nde_pin_spec *pins;
	gchar       **roles;   /* owns the role strings the pins point at */
	guint         n_pins;
};

nde_composite_capture *nde_composite_capture_begin(GSList *layers,
                                                   gboolean raw_first) {
	guint n = g_slist_length(layers);
	if (!n)
		return NULL;
	nde_composite_state *st = state_from_layers(layers, raw_first);
	guint n_inputs = st->inputs->len;
	guint n_masked = 0;
	for (guint i = 0; i < n_inputs; i++)
		if (g_array_index(st->inputs, nde_composite_input, i).mask_item_id)
			n_masked++;

	nde_composite_capture *cap = g_new0(nde_composite_capture, 1);
	cap->params = state_encode(st);
	cap->n_pins = n_inputs + n_masked;
	cap->pins   = g_new0(nde_pin_spec, cap->n_pins);
	cap->roles  = g_new0(gchar *, cap->n_pins);
	guint p = 0;
	for (guint i = 0; i < n_inputs; i++) {
		const nde_composite_input *in = &g_array_index(st->inputs, nde_composite_input, i);
		/* Merge-down keeps its two named roles; they read better in the graph
		 * view than in0/in1 and predate it.  Nothing parses them. */
		if (raw_first && n_inputs == 2)
			cap->roles[p] = g_strdup(i == 0 ? NDE_COMPOSITE_ROLE_BASE
			                               : NDE_COMPOSITE_ROLE_OVERLAY);
		else
			cap->roles[p] = g_strdup_printf("in%u", i);
		cap->pins[p].role      = cap->roles[p];
		cap->pins[p].item_id   = in->item_id;
		cap->pins[p].record_id = nde_history_last_record_for_item(in->item_id);
		p++;
		if (in->mask_item_id) {
			cap->roles[p] = g_strdup_printf("mask%u", i);
			cap->pins[p].role      = cap->roles[p];
			cap->pins[p].item_id   = in->mask_item_id;
			cap->pins[p].record_id = nde_history_last_record_for_item(in->mask_item_id);
			p++;
		}
	}

	guint i = 0;
	for (GSList *l = layers; l; l = l->next, i++) {
		const flis_layer_t *lay = l->data;
		if (!lay)
			continue;
		nde_checkpoint_baseline_ensure(lay->fit, lay->item_id);
		/* The mask is stored, not re-derived: a painted or loaded one has no
		 * chain to replay.  The coordinate is its pin's, so the mask cascade
		 * refreshes this copy when the mask IS built by ops and one of them is
		 * amended (nde_composite.h). */
		const nde_composite_input *in = (i < n_inputs) ?
				&g_array_index(st->inputs, nde_composite_input, i) : NULL;
		if (in && in->mask_item_id) {
			fits *mfit = lmask_to_fits(lay);
			if (mfit) {
				nde_checkpoint_store_at(mfit, in->mask_item_id,
				                        nde_history_last_record_for_item(in->mask_item_id));
				clearfits(mfit);
				free(mfit);
			}
		}
	}
	nde_composite_state_free(st);
	return cap;
}

void nde_composite_capture_free(nde_composite_capture *cap) {
	if (!cap)
		return;
	g_free(cap->params);
	for (guint i = 0; i < cap->n_pins; i++)
		g_free(cap->roles[i]);
	g_free(cap->roles);
	g_free(cap->pins);
	g_free(cap);
}

gint64 nde_composite_capture_commit(nde_composite_capture *cap, const char *op_id,
                                    gint target_item_id, const char *summary) {
	if (!cap)
		return 0;
	gint64 id = nde_capture_structural_pinned(op_id, NDE_SCOPE_DOCUMENT,
	                                          target_item_id, cap->params,
	                                          summary, cap->pins, cap->n_pins);
	cap->params = NULL;   /* consumed by the capture */
	nde_composite_capture_free(cap);
	return id;
}
