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

#include <glib.h>
#include <glib/gi18n.h>

#include "core/nde_op_class.h"
#include "core/nde_history.h"
#include "core/nde_composite.h"
#include "core/nde_compositing.h"
#include "core/processing.h"   /* destroy_any_args */

/* The ids with no op_descriptor.  THIS IS THE ONLY PLACE THEY APPEAR: every
 * other site asks nde_op_class_for() instead of comparing strings, which is
 * the whole point of the module.
 *
 * These thirteen are not an op_descriptor's business and are not waiting to
 * become one.  A descriptor is what the generic workers fill their args from —
 * it exists to be RUN — and none of these is ever run by any worker: a
 * structural record is never replayed at all, and the composite and
 * compositing-state families have their own strategies in nde_composite.c and
 * nde_compositing.c.  Giving them descriptors would mean thirteen entries with
 * every hook NULL, and the conformance tests that protect the other 85 (at
 * least one hook; serialize/deserialize paired or a recorded reason) would each
 * grow an exemption list.  That is more for an op author to know, not less.
 *
 * What they DO share with descriptor-backed ops is asked through this module:
 * see nde_op_class_params_valid(). */
static const nde_op_class descriptorless[] = {
	/* ---- composites: consume other items and replace this one's pixels --- */
	{ "layer.merge_down",  NDE_OPC_COMPOSITE, NDE_OPT_DESTRUCTIVE |
	                                          NDE_OPT_INSERT_DISTURBS,
	  NULL, N_("Merge down") },
	/* A flatten additionally resets the surviving base layer's compositing to
	 * the defaults (flis_flatten_all), which is what makes it a fold reset
	 * point as well as a composite. */
	{ "document.flatten",  NDE_OPC_COMPOSITE, NDE_OPT_DESTRUCTIVE |
	                                          NDE_OPT_INSERT_DISTURBS |
	                                          NDE_OPT_COMPOSITING_RESET,
	  NULL, N_("Flatten document") },

	/* ---- compositing state: inputs to the compositor, never chain members - */
	{ "layer.set_opacity", NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Set opacity") },
	{ "layer.set_blend",   NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Set blend mode") },
	{ "layer.set_visible", NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Set visibility") },
	{ "layer.set_tint",    NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Set tint") },

	/* ---- structural: change the document's shape, not any item's pixels --- */
	/* A new layer starts at the compositing defaults, so like a flatten it is
	 * a point past which earlier compositing records are superseded. */
	{ "layer.add",         NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE |
	                                           NDE_OPT_COMPOSITING_RESET,
	  NULL, N_("Add layer") },
	{ "layer.duplicate",   NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Duplicate layer") },
	{ "layer.remove",      NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE |
	                                           NDE_OPT_DELETES_ITEM,
	  NULL, N_("Remove layer") },
	{ "layer.reorder",     NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Reorder layers") },
	{ "canvas.resize",     NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Resize canvas") },
	{ "canvas.fit",        NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Fit canvas to layers") },
	/* Says where the BASELINE came from, which is where a replay already
	 * starts; running it would mean re-opening the file. */
	{ NDE_OP_IMAGE_ORIGIN, NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE,
	  NULL, N_("Image origin") },
};

/* Ops that HAVE a descriptor but whose family cannot be read off it.
 *
 * The joint ops look like ordinary image ops from the descriptor's side — they
 * have an image_hook and a serializer — but they read several items and are
 * members of every participant's chain.
 *
 * The analysis ops report a number or write a file and never modify the item.
 * That is a property of their CALL SITES (each one passes skip_generic_undo, or
 * for_preview, which is the same gate the NDE capture tests), so it cannot be
 * derived from the descriptor either; naming them here is what lets the worker
 * assert it once instead of trusting twenty call sites to remember.  The list
 * is cross-checked against opaque_by_design[] in test_op_descriptor.c. */
static const struct {
	const char   *id;
	nde_op_family family;
} classified_descriptors[] = {
	{ "flis.layers_match",       NDE_OPC_JOINT },
	{ "flis.group_calibration",  NDE_OPC_JOINT },
	{ "flis.register",           NDE_OPC_JOINT },

	{ "stats.bg",                NDE_OPC_ANALYSIS },
	{ "stats.bgnoise",           NDE_OPC_ANALYSIS },
	{ "stats.cdg",               NDE_OPC_ANALYSIS },
	{ "stats.entropy",           NDE_OPC_ANALYSIS },
	{ "stats.stat",              NDE_OPC_ANALYSIS },
	{ "psf.estimate",            NDE_OPC_ANALYSIS },
	{ "catalog.search",          NDE_OPC_ANALYSIS },
	{ "cfa.split",               NDE_OPC_ANALYSIS },
	{ "cfa.extract_green",       NDE_OPC_ANALYSIS },
	{ "cfa.extract_ha",          NDE_OPC_ANALYSIS },
	{ "cfa.extract_haoiii",      NDE_OPC_ANALYSIS },
	{ "cfa.findhot",             NDE_OPC_ANALYSIS },
};

/* An id from a newer build, or a corrupt one.  No traits: every "is this
 * ignorable / harmless" test is written as a trait CHECK, so absence of the
 * trait is refusal, and a class we have never heard of refuses everything. */
static const nde_op_class unknown_class = {
	NULL, NDE_OPC_UNKNOWN, 0, NULL, NULL
};

static GHashTable *class_table;   /* const char* -> const nde_op_class* */

static nde_op_family family_from_descriptor(const op_descriptor *op) {
	/* Order matters: icc.convert has BOTH an image_hook (plain images) and a
	 * layer_hook (FLIS documents), and it is the layer_hook that makes it
	 * document-wide. */
	if (op->layer_hook)
		return NDE_OPC_DOCUMENT;
	if (op->mask_hook)
		return NDE_OPC_MASK;
	return NDE_OPC_PIXEL;
}

static gpointer class_table_init(gpointer unused) {
	(void)unused;
	GHashTable *t = g_hash_table_new(g_str_hash, g_str_equal);

	for (guint i = 0; i < G_N_ELEMENTS(descriptorless); i++)
		g_hash_table_insert(t, (gpointer)descriptorless[i].id,
		                    (gpointer)&descriptorless[i]);

	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	for (size_t i = 0; i < n; i++) {
		const op_descriptor *op = all[i];
		if (!op || !op->id)
			continue;
		nde_op_class *cls = g_new0(nde_op_class, 1);
		cls->id     = op->id;
		cls->desc   = op;
		cls->label  = op->description;
		cls->family = family_from_descriptor(op);
		if (op->flags & OP_GEOMETRY_CHANGING)
			cls->traits |= NDE_OPT_GEOMETRIC;
		for (guint k = 0; k < G_N_ELEMENTS(classified_descriptors); k++) {
			if (!g_strcmp0(op->id, classified_descriptors[k].id)) {
				cls->family = classified_descriptors[k].family;
				break;
			}
		}
		/* The two tables must be disjoint: an id with a descriptor is
		 * classified from it, and one without is classified by hand.  An id in
		 * both means one of them is stale. */
		g_warn_if_fail(!g_hash_table_contains(t, op->id));
		g_hash_table_insert(t, (gpointer)cls->id, cls);
	}
	class_table = t;
	return NULL;
}

const nde_op_class *nde_op_class_for(const char *op_id) {
	static GOnce once = G_ONCE_INIT;
	g_once(&once, class_table_init, NULL);
	if (!op_id)
		return &unknown_class;
	const nde_op_class *cls = g_hash_table_lookup(class_table, op_id);
	return cls ? cls : &unknown_class;
}

gboolean nde_op_class_params_valid(const char *op_id, int op_version,
                                   const char *old_params,
                                   const char *new_params, gchar **err) {
	const nde_op_class *cls = nde_op_class_for(op_id);
	/* A switch, not a vtable.  There is a small closed set of ways to check a
	 * params blob and they differ in what they need to see, not merely in
	 * behaviour — the composite is the only one that compares against the
	 * recorded blob — so a function pointer would buy an indirection with one
	 * consumer to justify it. */
	switch (cls->family) {
	case NDE_OPC_COMPOSITE:
		/* Against the RECORDED blob: the compositing state may change, what
		 * the node consumed may not (nde_composite.h). */
		return nde_composite_validate(old_params, new_params, err);
	case NDE_OPC_COMPOSITING:
		return nde_compositing_validate(op_id, new_params, err);
	default:
		break;
	}
	if (!cls->desc || !cls->desc->deserialize) {
		/* Either a structural step, which describes something that happened to
		 * the document rather than a computation with settings, or an id from a
		 * build that knows more ops than this one. */
		*err = g_strdup_printf(_("operation '%s' cannot be edited by this build"),
		                       op_id ? op_id : "?");
		return FALSE;
	}
	/* The round trip IS the validation: a blob the deserializer accepts is one
	 * the replay can run. */
	gpointer trial = cls->desc->deserialize(new_params, op_version);
	if (!trial) {
		*err = g_strdup_printf(_("the new parameters for '%s' failed to parse"),
		                       op_id ? op_id : "?");
		return FALSE;
	}
	destroy_any_args(trial);
	return TRUE;
}
