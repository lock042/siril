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

/*
 * test_nde_serializers — NDE phase-1 serialize/deserialize round-trips
 * (flis-nde-sketch.md §10-§15).  Every op with an op_descriptor serializer
 * builds a params struct with distinctive non-default values (negatives,
 * binary-inexact floats, enums > 0), serializes, deserializes, and asserts
 * bit-exact field equality plus a valid destructor.  Malformed and
 * newer-version blobs must deserialize to NULL.
 */

#include <criterion/criterion.h>
#include <string.h>
#include "core/siril.h"
#include "core/op_descriptor.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"

#include "algos/geometry.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image (a pointer)

/* Common malformed-blob checks: empty blob and a version one past the
 * descriptor's own version both yield NULL. */
#define CHECK_MALFORMED(op, goodblob) do {                                    \
	cr_assert_null((op)->deserialize("", 1),                              \
	               "%s: empty blob must fail", (op)->id);                 \
	cr_assert_null((op)->deserialize((goodblob), (op)->version + 1),      \
	               "%s: newer version must fail", (op)->id);              \
} while (0)

/* destructor-first free helper: assert non-NULL then invoke it. */
#define FREE_VIA_DESTRUCTOR(copy) do {                                        \
	destructor *df = (destructor *)(copy);                                \
	cr_assert_not_null(*df, "deserialized struct has NULL destroy_fn");   \
	(*df)((copy));                                                        \
} while (0)

/* ------------------------------------------------------------------ *
 *  Batch 1 — geometry (all POD)                                      *
 * ------------------------------------------------------------------ */

Test(nde_serializers, crop_roundtrip) {
	struct crop_args in = { 0 };
	in.area.x = -7; in.area.y = 13; in.area.w = 640; in.area.h = 481;

	gchar *blob = op_desc_crop.serialize(&in);
	cr_assert_not_null(blob);
	struct crop_args *out = op_desc_crop.deserialize(blob, op_desc_crop.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->area.x, in.area.x);
	cr_assert_eq(out->area.y, in.area.y);
	cr_assert_eq(out->area.w, in.area.w);
	cr_assert_eq(out->area.h, in.area.h);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_crop, blob);
	g_free(blob);
}

Test(nde_serializers, rotation_roundtrip) {
	struct rotation_args in = { 0 };
	in.area.x = 3; in.area.y = -5; in.area.w = 800; in.area.h = 600;
	in.angle = 37.125;
	in.interpolation = OPENCV_LANCZOS4;   /* 4 */
	in.cropped = 1;
	in.clamp = TRUE;

	gchar *blob = op_desc_rotation.serialize(&in);
	cr_assert_not_null(blob);
	struct rotation_args *out = op_desc_rotation.deserialize(blob, op_desc_rotation.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->area.x, in.area.x);
	cr_assert_eq(out->area.y, in.area.y);
	cr_assert_eq(out->area.w, in.area.w);
	cr_assert_eq(out->area.h, in.area.h);
	cr_assert(memcmp(&out->angle, &in.angle, sizeof(double)) == 0);
	cr_assert_eq(out->interpolation, in.interpolation);
	cr_assert_eq(out->cropped, in.cropped);
	cr_assert_eq(out->clamp, in.clamp);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_rotation, blob);
	g_free(blob);
}

Test(nde_serializers, mirror_roundtrip) {
	/* mirrorx and mirrory share the serializer pair; test both x_axis values. */
	struct mirror_args in = { 0 };
	in.x_axis = TRUE;
	gchar *blob = op_desc_mirrorx.serialize(&in);
	cr_assert_not_null(blob);
	struct mirror_args *out = op_desc_mirrorx.deserialize(blob, op_desc_mirrorx.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->x_axis, TRUE);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mirrorx, blob);
	g_free(blob);

	in.x_axis = FALSE;
	blob = op_desc_mirrory.serialize(&in);
	out = op_desc_mirrory.deserialize(blob, op_desc_mirrory.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->x_axis, FALSE);
	FREE_VIA_DESTRUCTOR(out);
	g_free(blob);
}

Test(nde_serializers, binning_roundtrip) {
	struct binning_args in = { 0 };
	in.factor = 3;
	in.mean = TRUE;
	gchar *blob = op_desc_binning.serialize(&in);
	cr_assert_not_null(blob);
	struct binning_args *out = op_desc_binning.deserialize(blob, op_desc_binning.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->factor, in.factor);
	cr_assert_eq(out->mean, in.mean);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_binning, blob);
	g_free(blob);
}

Test(nde_serializers, resample_roundtrip) {
	struct resample_args in = { 0 };
	in.toX = 1024; in.toY = 768;
	in.interpolation = OPENCV_CUBIC;   /* 2 */
	in.clamp = TRUE;
	in.update_wcs = FALSE;
	gchar *blob = op_desc_resample.serialize(&in);
	cr_assert_not_null(blob);
	struct resample_args *out = op_desc_resample.deserialize(blob, op_desc_resample.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->toX, in.toX);
	cr_assert_eq(out->toY, in.toY);
	cr_assert_eq(out->interpolation, in.interpolation);
	cr_assert_eq(out->clamp, in.clamp);
	cr_assert_eq(out->update_wcs, in.update_wcs);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_resample, blob);
	g_free(blob);
}

/* ------------------------------------------------------------------ *
 *  Registry: the set of ops with serializers is exactly phase 1.     *
 *  Keeps the set deliberate — a new serializer without an entry here *
 *  (or an accidental one) fails the build.                           *
 * ------------------------------------------------------------------ */

static const char *phase1_ids[] = {
	/* batch 1 — geometry */
	"geometry.crop", "geometry.rotation", "geometry.mirrorx",
	"geometry.mirrory", "geometry.binning", "geometry.resample",
};

Test(nde_serializers, serializer_set_is_phase1) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	cr_assert(all != NULL);
	for (size_t i = 0; i < n; i++) {
		const op_descriptor *op = all[i];
		if (!op->serialize)
			continue;
		gboolean found = FALSE;
		for (size_t j = 0; j < G_N_ELEMENTS(phase1_ids); j++) {
			if (!strcmp(op->id, phase1_ids[j])) { found = TRUE; break; }
		}
		cr_assert(found, "descriptor '%s' has a serializer but is not in the "
		                 "expected phase-1 set", op->id);
	}
}
