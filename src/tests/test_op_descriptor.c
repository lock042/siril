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
#include "core/processing.h"
#include "core/op_descriptor.h"
#include "core/nde_op_class.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image (now a pointer)

/* ------------------------------------------------------------------ *
 *  Registry invariants — checked over every descriptor in the build  *
 * ------------------------------------------------------------------ */

/* matches ^[a-z0-9_]+\.[a-z0-9_]+$ */
static gboolean valid_id(const char *id) {
	if (!id)
		return FALSE;
	const char *dot = strchr(id, '.');
	if (!dot || dot == id || *(dot + 1) == '\0')
		return FALSE;
	if (strchr(dot + 1, '.'))
		return FALSE;   /* exactly one dot */
	for (const char *p = id; *p; p++) {
		if (p == dot)
			continue;
		if (!((*p >= 'a' && *p <= 'z') || (*p >= '0' && *p <= '9') || *p == '_'))
			return FALSE;
	}
	return TRUE;
}

Test(op_descriptor, registry_invariants) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	cr_assert(all != NULL);

	for (size_t i = 0; i < n; i++) {
		const op_descriptor *op = all[i];
		cr_assert(op != NULL, "descriptor %zu is NULL", i);
		cr_assert(valid_id(op->id), "descriptor '%s' has malformed id", op->id ? op->id : "(null)");
		cr_assert(op->version >= 1, "descriptor '%s' has version < 1", op->id);
		/* Hook invariant (relaxed for the NDE layer_hook, sketch §11):
		 * mask_hook is exclusive of the other two — a mask op is never also an
		 * image/layer op — but image_hook and layer_hook MAY coexist on one
		 * descriptor (the same logical op run on a plain image vs a FLIS
		 * document; icc.convert is exactly this).  At least one hook must be
		 * set. */
		gboolean has_img   = op->image_hook != NULL;
		gboolean has_mask  = op->mask_hook  != NULL;
		gboolean has_layer = op->layer_hook != NULL;
		cr_assert(has_img || has_mask || has_layer,
		          "descriptor '%s' must set at least one hook", op->id);
		if (has_mask)
			cr_assert(!has_img && !has_layer,
			          "descriptor '%s' mask_hook is exclusive of image/layer hooks", op->id);
		cr_assert(op->description != NULL, "descriptor '%s' has NULL description", op->id);
		/* NDE members: serialize/deserialize are paired — a descriptor that
		 * can serialize its params must also be able to deserialize them, and
		 * vice versa.  Either both set (Tier A) or both NULL (Tier B). */
		cr_assert_eq(!!op->serialize, !!op->deserialize,
		             "descriptor '%s' must set serialize and deserialize together", op->id);
	}
}

Test(op_descriptor, ids_unique) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			cr_assert(strcmp(all[i]->id, all[j]->id) != 0,
			          "duplicate descriptor id '%s'", all[i]->id);
		}
	}
}

/* The exact set of ops that carry no serializer.
 *
 * This is the counterpart of the registry_invariants pairing assertion above.
 * That one catches a HALF-written serializer (one slot set, the other not);
 * this catches a MISSING one — an op added, or a params struct changed, where
 * nobody wrote the serializer at all.  Nothing else catches that: the op runs,
 * the record is captured, and the only symptom is that the user's edit history
 * quietly stops being editable from that step backwards.
 *
 * "No serializer" is a supported outcome, not a bug — but it is supposed to be
 * a DECISION.  Naming the whole set here is what makes it one: adding an op
 * without a serializer fails this test, and the fix is either to write the
 * serializer or to add the id below with its reason.
 *
 * The reason is recorded per entry rather than in a prose block, because the
 * two kinds are NOT equivalent: only NOT_CAPTURED is free, and only TIER_B
 * costs the user anything. */
typedef enum {
	/* Never reaches the NDE capture block at all, so there is no record and
	 * therefore no barrier.  Every call site is measurement-only and sets
	 * skip_generic_undo — or, for the deconvolution GUI's PSF estimate,
	 * for_preview — and that is the SAME gate the capture uses (see "NDE
	 * provenance capture" in processing.c).  A serializer here would be dead
	 * code: there is nothing to attach it to.
	 *
	 * That property lives at the CALL SITE, not on the descriptor, so this
	 * test cannot check it directly.  What it can do is require every id here
	 * to be declared NDE_OPC_ANALYSIS in the op-class registry, which is what
	 * the capture block asserts against — so the two lists cannot drift, and
	 * a call site that drops skip_generic_undo trips the warning there. */
	OPAQUE_NOT_CAPTURED,
	/* Genuinely captured, as a Tier-B barrier.  These are the ones that cost
	 * the user editability, and each is a considered verdict. */
	OPAQUE_TIER_B,
} opaque_reason;

static const struct {
	const char    *id;
	opaque_reason  why;
} opaque_by_design[] = {
	/* Analysis only — report a number, never touch pixels. */
	{ "stats.bg",             OPAQUE_NOT_CAPTURED },
	{ "stats.bgnoise",        OPAQUE_NOT_CAPTURED },
	{ "stats.cdg",            OPAQUE_NOT_CAPTURED },
	{ "stats.entropy",        OPAQUE_NOT_CAPTURED },
	{ "stats.stat",           OPAQUE_NOT_CAPTURED },
	{ "psf.estimate",         OPAQUE_NOT_CAPTURED },
	/* Reaches the network, so its result is not reproducible — but it
	 * annotates, it does not modify, so the question never arises. */
	{ "catalog.search",       OPAQUE_NOT_CAPTURED },
	/* Write NEW images out; they leave the loaded one alone. */
	{ "cfa.split",            OPAQUE_NOT_CAPTURED },
	{ "cfa.extract_green",    OPAQUE_NOT_CAPTURED },
	{ "cfa.extract_ha",       OPAQUE_NOT_CAPTURED },
	{ "cfa.extract_haoiii",   OPAQUE_NOT_CAPTURED },
	{ "cfa.findhot",          OPAQUE_NOT_CAPTURED },

	/* Their operands are files a SEPARATE earlier command left in the temp
	 * directory, so a params blob cannot describe what they will consume.
	 * Permanent, by maintainer verdict recorded in the descriptors (NDE.rst).
	 * Note wavelets.wrecons is also OP_ROI_CAPABLE: it can compute a rectangle
	 * for a live preview and still be unreplayable afterwards.  The two
	 * properties are independent. */
	{ "filters.fft",          OPAQUE_TIER_B },
	{ "wavelets.wrecons",     OPAQUE_TIER_B },
	/* Mask op: derives a mask from the image it is handed.  The mask worker
	 * captures unconditionally — there is no skip_generic_undo on that path —
	 * and OP_MASK_FROM_IMAGE pins the image so the mask chain re-derives the
	 * mask instead of replaying it from parameters. */
	{ "mask.mtf_autostretch", OPAQUE_TIER_B },
};

#define N_OPAQUE (sizeof(opaque_by_design) / sizeof(opaque_by_design[0]))

static gboolean id_in_opaque_list(const char *id) {
	for (size_t i = 0; i < N_OPAQUE; i++)
		if (!strcmp(id, opaque_by_design[i].id))
			return TRUE;
	return FALSE;
}

Test(op_descriptor, every_descriptor_is_replayable_or_deliberately_opaque) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	size_t seen = 0;

	for (size_t i = 0; i < n; i++) {
		const gboolean replayable = all[i]->serialize != NULL;
		const gboolean expected_opaque = id_in_opaque_list(all[i]->id);
		cr_assert_eq(replayable, !expected_opaque,
		             "'%s': %s, but the pinned list says it should be %s.  "
		             "Write serialize/deserialize for it, or add it to "
		             "opaque_by_design[] with the reason.",
		             all[i]->id,
		             replayable ? "it has a serializer" : "it has no serializer",
		             expected_opaque ? "opaque" : "replayable");
		if (!replayable)
			seen++;
	}

	/* Catches an id removed from the registry, or given a serializer, but left
	 * in the list above — which the per-descriptor loop alone would not see. */
	cr_assert_eq(seen, N_OPAQUE,
	             "%zu descriptors carry no serializer, list names %zu",
	             seen, N_OPAQUE);
}

/* The Tier-B set is the part of the list that actually costs the user
 * editability, so it is pinned separately and kept small on purpose.  An op
 * moving from NOT_CAPTURED to TIER_B is a real regression in the edit history
 * — and, because the difference lives at the call site, one the sweep above
 * cannot see.  Asserting the count is what forces that move to be argued for
 * in review rather than discovered later. */
Test(op_descriptor, the_tier_b_barrier_set_stays_small) {
	size_t barriers = 0;
	for (size_t i = 0; i < N_OPAQUE; i++)
		if (opaque_by_design[i].why == OPAQUE_TIER_B)
			barriers++;
	cr_assert_eq(barriers, 3,
	             "%zu ops are opaque Tier-B barriers, expected 3.  Adding one "
	             "means a step the user can no longer edit past: say why in the "
	             "entry's comment and update this count.",
	             barriers);
}

/* The cross-check this file could not write on its own: every OPAQUE_NOT_
 * CAPTURED id must be declared NDE_OPC_ANALYSIS in the op-class registry, and
 * nothing else here may be.  The capture block warns when an NDE_OPC_ANALYSIS
 * op reaches it, so this is what connects "has no serializer because it never
 * captures" to a statement the running code actually enforces. */
Test(op_descriptor, analysis_ops_agree_with_the_op_class_registry) {
	for (size_t i = 0; i < N_OPAQUE; i++) {
		nde_op_family fam = nde_op_class_for(opaque_by_design[i].id)->family;
		if (opaque_by_design[i].why == OPAQUE_NOT_CAPTURED)
			cr_assert_eq(fam, NDE_OPC_ANALYSIS,
			             "'%s' is pinned as never-captured here, but the op-class "
			             "registry does not call it NDE_OPC_ANALYSIS — so nothing "
			             "checks that its call sites still skip the capture",
			             opaque_by_design[i].id);
		else
			cr_assert_neq(fam, NDE_OPC_ANALYSIS,
			              "'%s' is a Tier-B barrier, so it IS captured, but the "
			              "registry calls it an analysis op",
			              opaque_by_design[i].id);
	}
}

/* The exact set of ops that may be computed on a sub-rectangle.
 *
 * The flag is now what generic_image_worker consults to decide whether a
 * preview run is region-scoped, so an op silently gaining or losing it changes
 * behaviour.  Pinning the whole set is what makes such a change fail here
 * rather than in a user's viewport.  The dialogs no longer carry a second
 * copy of this truth: roi_declare_op() names the op and reads the flag.
 *
 * wavelets is in the set for a different reason from the rest: its hook cannot
 * compute from a crop, because its input is a decomposition of the whole
 * image — but it can produce just the window, reading only the transform rows
 * that window covers. */
static const char *const roi_capable_ids[] = {
	"color.saturation",
	"filters.deconvolve",
	"filters.denoise",
	"filters.epf",
	"filters.median",
	"filters.scnr",
	"filters.unpurple",
	"stretch.asinh",
	"stretch.curves",
	"stretch.ghs",
	"stretch.mtf",
	"wavelets.wrecons",
};

static gboolean id_in_roi_list(const char *id) {
	for (size_t i = 0; i < sizeof(roi_capable_ids) / sizeof(roi_capable_ids[0]); i++)
		if (!strcmp(id, roi_capable_ids[i]))
			return TRUE;
	return FALSE;
}

Test(op_descriptor, roi_capable_set_is_exactly_the_expected_ops) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	size_t seen = 0;

	for (size_t i = 0; i < n; i++) {
		const gboolean flagged = (all[i]->flags & OP_ROI_CAPABLE) != 0;
		const gboolean expected = id_in_roi_list(all[i]->id);
		cr_assert_eq(flagged, expected,
		             "'%s': OP_ROI_CAPABLE is %s but the pinned list says %s",
		             all[i]->id, flagged ? "set" : "unset",
		             expected ? "it should be set" : "it should not be");
		if (flagged)
			seen++;
	}

	/* Catches an id removed from the registry but left in the list above,
	 * which the per-descriptor loop alone would not see. */
	cr_assert_eq(seen, sizeof(roi_capable_ids) / sizeof(roi_capable_ids[0]),
	             "%zu descriptors carry OP_ROI_CAPABLE, list names %zu",
	             seen, sizeof(roi_capable_ids) / sizeof(roi_capable_ids[0]));
}

/* ------------------------------------------------------------------ *
 *  Fill semantics                                                    *
 * ------------------------------------------------------------------ */

static int dummy_image_hook(struct generic_img_args *args, fits *fit, int nb) {
	(void)args; (void)fit; (void)nb;
	return 0;
}

static gchar *dummy_log_hook(gpointer user, log_hook_detail detail) {
	(void)user; (void)detail;
	return g_strdup("dummy");
}

static const op_descriptor op_desc_test = {
	.id = "test.fill",
	.version = 1,
	.image_hook = dummy_image_hook,
	.log_hook = dummy_log_hook,
	.description = N_("Test op"),
	.mem_ratio = 2.5f,
	.flags = OP_MASK_CAPABLE,
};

Test(op_descriptor, fill_populates_from_descriptor) {
	struct generic_img_args args = { 0 };
	args.op = &op_desc_test;
	op_descriptor_fill_img_args(&args);
	cr_assert(args.image_hook == dummy_image_hook);
	cr_assert(args.log_hook == dummy_log_hook);
	cr_assert(args.description != NULL);
	cr_assert_str_eq(args.description, _("Test op"));
	cr_assert_float_eq(args.mem_ratio, 2.5f, 1e-6);
}

Test(op_descriptor, fill_respects_site_overrides) {
	struct generic_img_args args = { 0 };
	args.op = &op_desc_test;
	args.description = "override";   /* pre-set: descriptor must not clobber */
	args.mem_ratio = 1.0f;           /* non-zero: descriptor default not applied */
	op_descriptor_fill_img_args(&args);
	cr_assert_str_eq(args.description, "override");
	cr_assert_float_eq(args.mem_ratio, 1.0f, 1e-6);
}

Test(op_descriptor, fill_noop_when_no_op) {
	struct generic_img_args args = { 0 };
	args.image_hook = dummy_image_hook;
	op_descriptor_fill_img_args(&args);   /* op == NULL → no-op, no assert */
	cr_assert(args.image_hook == dummy_image_hook);
	cr_assert(args.log_hook == NULL);
}

/* Setting both a descriptor and a per-site image_hook is a programming error;
 * op_descriptor_fill_img_args g_assert()s against it. */
Test(op_descriptor, fill_asserts_on_double_hook, .signal = SIGABRT) {
	struct generic_img_args args = { 0 };
	args.op = &op_desc_test;
	args.image_hook = dummy_image_hook;   /* both set → abort */
	op_descriptor_fill_img_args(&args);
}

/* ------------------------------------------------------------------ *
 *  Layer-args fill (NDE, steps 6-8 Part C)                           *
 * ------------------------------------------------------------------ */

static int dummy_layer_hook(struct generic_layer_args *a) { (void)a; return 0; }

/* A descriptor that carries BOTH an image and a layer hook — the new
 * coexistence invariant (icc.convert is real-world proof). */
static const op_descriptor op_desc_layer_test = {
	.id = "test.layerfill",
	.version = 2,
	.image_hook = dummy_image_hook,
	.layer_hook = dummy_layer_hook,
	.log_hook = dummy_log_hook,
	.description = N_("Layer test op"),
	.mem_ratio = 3.0f,
};

Test(op_descriptor, fill_layer_populates_from_descriptor) {
	struct generic_layer_args args = { 0 };
	args.op = &op_desc_layer_test;
	op_descriptor_fill_layer_args(&args);
	cr_assert(args.layer_hook == dummy_layer_hook);
	cr_assert(args.log_hook == dummy_log_hook);
	cr_assert_str_eq(args.description, _("Layer test op"));
	cr_assert_float_eq(args.mem_ratio, 3.0f, 1e-6);
}

Test(op_descriptor, fill_layer_noop_when_no_op) {
	struct generic_layer_args args = { 0 };
	args.layer_hook = dummy_layer_hook;
	op_descriptor_fill_layer_args(&args);   /* op == NULL → no-op, no assert */
	cr_assert(args.layer_hook == dummy_layer_hook);
	cr_assert(args.log_hook == NULL);
}

Test(op_descriptor, fill_layer_asserts_on_double_hook, .signal = SIGABRT) {
	struct generic_layer_args args = { 0 };
	args.op = &op_desc_layer_test;
	args.layer_hook = dummy_layer_hook;   /* both set → abort */
	op_descriptor_fill_layer_args(&args);
}
