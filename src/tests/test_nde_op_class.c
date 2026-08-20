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
 * along with Siril. If not, see <http://www.gnu.org/licenses/ >.
 */

#include <criterion/criterion.h>
#include <string.h>
#include "core/siril.h"
#include "core/op_descriptor.h"
#include "core/nde_op_class.h"

cominfo com;
fits *gfit;

/* The classification every other module now depends on.  These are pinned
 * expectations, not derivations: the point of the table is that changing a
 * family or a trait is a deliberate act, so a test that recomputed them from
 * the same table would assert nothing. */
static const struct {
	const char   *id;
	nde_op_family family;
	guint32       traits;
} expected[] = {
	/* ---- descriptor-less families ---- */
	{ "layer.merge_down",  NDE_OPC_COMPOSITE,
	  NDE_OPT_DESTRUCTIVE | NDE_OPT_INSERT_DISTURBS },
	{ "document.flatten",  NDE_OPC_COMPOSITE,
	  NDE_OPT_DESTRUCTIVE | NDE_OPT_INSERT_DISTURBS | NDE_OPT_COMPOSITING_RESET },
	{ "layer.set_opacity", NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE },
	{ "layer.set_blend",   NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE },
	{ "layer.set_visible", NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE },
	{ "layer.set_tint",    NDE_OPC_COMPOSITING, NDE_OPT_CHAIN_IGNORE },
	{ "layer.add",         NDE_OPC_STRUCTURAL,
	  NDE_OPT_CHAIN_IGNORE | NDE_OPT_COMPOSITING_RESET },
	{ "layer.duplicate",   NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE },
	{ "layer.remove",      NDE_OPC_STRUCTURAL,
	  NDE_OPT_CHAIN_IGNORE | NDE_OPT_DELETES_ITEM },
	{ "layer.reorder",     NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE },
	{ "canvas.resize",     NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE },
	{ "canvas.fit",        NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE },
	{ "image.origin",      NDE_OPC_STRUCTURAL, NDE_OPT_CHAIN_IGNORE },

	/* ---- descriptors whose family is not readable off the descriptor ---- */
	{ "flis.layers_match",      NDE_OPC_JOINT, 0 },
	{ "flis.group_calibration", NDE_OPC_JOINT, 0 },
	/* the only geometric joint op: it moves participants as well as warping */
	{ "flis.register",          NDE_OPC_JOINT, NDE_OPT_GEOMETRIC },
	{ "stats.bg",               NDE_OPC_ANALYSIS, 0 },
	{ "psf.estimate",           NDE_OPC_ANALYSIS, 0 },
	{ "cfa.split",              NDE_OPC_ANALYSIS, 0 },

	/* ---- families read straight off the descriptor ---- */
	{ "stretch.mtf",       NDE_OPC_PIXEL,    0 },
	{ "mask.blur",         NDE_OPC_MASK,     0 },
	{ "icc.convert",       NDE_OPC_DOCUMENT, 0 },
	{ "geometry.rotation", NDE_OPC_PIXEL,    NDE_OPT_GEOMETRIC },
	{ "geometry.crop",     NDE_OPC_PIXEL,    NDE_OPT_GEOMETRIC },
};

Test(nde_op_class, pinned_ids_classify_as_expected) {
	for (size_t i = 0; i < G_N_ELEMENTS(expected); i++) {
		const nde_op_class *cls = nde_op_class_for(expected[i].id);
		cr_assert_neq(cls->family, NDE_OPC_UNKNOWN,
		              "'%s' is not in the registry at all", expected[i].id);
		cr_assert_eq(cls->family, expected[i].family,
		             "'%s': family %d, expected %d",
		             expected[i].id, cls->family, expected[i].family);
		cr_assert_eq(cls->traits, expected[i].traits,
		             "'%s': traits 0x%x, expected 0x%x",
		             expected[i].id, cls->traits, expected[i].traits);
	}
}

/* Every descriptor must classify, and must carry its descriptor back — the
 * label lookup in region_tail_member_reason() reads cls->desc and cls->label
 * instead of doing its own op_descriptor_by_id(). */
Test(nde_op_class, every_descriptor_resolves) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	cr_assert(n > 0);
	for (size_t i = 0; i < n; i++) {
		const nde_op_class *cls = nde_op_class_for(all[i]->id);
		cr_assert_neq(cls->family, NDE_OPC_UNKNOWN,
		              "descriptor '%s' does not classify", all[i]->id);
		cr_assert_eq(cls->desc, all[i],
		             "'%s' does not carry its descriptor", all[i]->id);
		cr_assert_eq(cls->label, all[i]->description,
		             "'%s' does not carry its description as a label", all[i]->id);
	}
}

/* The descriptor-less table and the descriptor registry must be disjoint: an
 * id in both would be classified twice, and whichever won would be silent. */
Test(nde_op_class, descriptorless_ids_have_no_descriptor) {
	for (size_t i = 0; i < G_N_ELEMENTS(expected); i++) {
		const nde_op_class *cls = nde_op_class_for(expected[i].id);
		if (cls->family != NDE_OPC_COMPOSITE &&
		    cls->family != NDE_OPC_COMPOSITING &&
		    cls->family != NDE_OPC_STRUCTURAL)
			continue;
		cr_assert_null(cls->desc,
		               "'%s' is in the descriptor-less table but HAS a descriptor",
		               expected[i].id);
		cr_assert_null(op_descriptor_by_id(expected[i].id),
		               "'%s' is classified by hand but the registry knows it too",
		               expected[i].id);
	}
}

/* THE fail-closed property.  An id from a newer build must not pick up
 * NDE_OPT_CHAIN_IGNORE by accident: chain_build_excluding() spells "harmless"
 * as the PRESENCE of that trait, so an unknown op keeps hard-blocking a
 * document-scope chain exactly as it did when this was a string comparison. */
Test(nde_op_class, unknown_ids_fail_closed) {
	const char *strangers[] = {
		"something.new", "layer.teleport", "", "layer.set_", "nonsense",
	};
	for (size_t i = 0; i < G_N_ELEMENTS(strangers); i++) {
		const nde_op_class *cls = nde_op_class_for(strangers[i]);
		cr_assert_eq(cls->family, NDE_OPC_UNKNOWN,
		             "'%s' should be unknown", strangers[i]);
		cr_assert_eq(cls->traits, 0,
		             "'%s' must carry NO traits, got 0x%x", strangers[i], cls->traits);
		cr_assert_eq(cls->traits & NDE_OPT_CHAIN_IGNORE, 0,
		             "'%s' must not be treated as ignorable", strangers[i]);
	}
	const nde_op_class *nul = nde_op_class_for(NULL);
	cr_assert_eq(nul->family, NDE_OPC_UNKNOWN, "NULL classifies as unknown");
	cr_assert_eq(nul->traits, 0, "NULL carries no traits");
}

/* Only composites are destructive, and every destructive op disturbs an armed
 * insertion.  Both directions, because a new destructive family that forgot
 * INSERT_DISTURBS would let an insertion survive something that invalidates
 * it — a silent corruption of the log rather than a crash. */
Test(nde_op_class, destructive_implies_insert_disturbs) {
	for (size_t i = 0; i < G_N_ELEMENTS(expected); i++) {
		const nde_op_class *cls = nde_op_class_for(expected[i].id);
		if (cls->traits & NDE_OPT_DESTRUCTIVE) {
			cr_assert_eq(cls->family, NDE_OPC_COMPOSITE,
			             "'%s' is destructive but not a composite", expected[i].id);
			cr_assert(cls->traits & NDE_OPT_INSERT_DISTURBS,
			          "'%s' is destructive but an insertion could survive it",
			          expected[i].id);
		}
		/* And the converse for CHAIN_IGNORE: nothing may be both ignorable
		 * and destructive, or chain_build_excluding's two arms disagree. */
		if (cls->traits & NDE_OPT_CHAIN_IGNORE)
			cr_assert_eq(cls->traits & NDE_OPT_DESTRUCTIVE, 0,
			             "'%s' is both ignorable and destructive", expected[i].id);
	}
}
