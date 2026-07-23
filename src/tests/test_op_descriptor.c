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
		/* exactly one of image_hook / mask_hook */
		cr_assert((op->image_hook != NULL) != (op->mask_hook != NULL),
		          "descriptor '%s' must have exactly one of image_hook/mask_hook", op->id);
		cr_assert(op->description != NULL, "descriptor '%s' has NULL description", op->id);
		/* NDE members reserved — always NULL in this MR */
		cr_assert(op->serialize == NULL, "descriptor '%s' sets serialize", op->id);
		cr_assert(op->deserialize == NULL, "descriptor '%s' sets deserialize", op->id);
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
