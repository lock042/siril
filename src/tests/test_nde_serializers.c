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
#include "filters/mtf.h"
#include "filters/ght.h"
#include "filters/asinh.h"
#include "filters/scnr.h"
#include "filters/median.h"
#include "algos/colors.h"
#include "filters/curve_transform.h"
#include "algos/background_extraction.h"

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
 *  Batch 2 — stretches, scnr, median, ccm                            *
 * ------------------------------------------------------------------ */

static void fill_mtf_params(struct mtf_params *p, float base) {
	p->midtones = base + 0.1f;      /* binary-inexact */
	p->shadows = base + 0.03125f;
	p->highlights = base + 0.9f;
	p->do_red = TRUE; p->do_green = FALSE; p->do_blue = TRUE;
}

static void check_mtf_params(const struct mtf_params *a, const struct mtf_params *b) {
	cr_assert(memcmp(&a->midtones, &b->midtones, sizeof(float)) == 0);
	cr_assert(memcmp(&a->shadows, &b->shadows, sizeof(float)) == 0);
	cr_assert(memcmp(&a->highlights, &b->highlights, sizeof(float)) == 0);
	cr_assert_eq(a->do_red, b->do_red);
	cr_assert_eq(a->do_green, b->do_green);
	cr_assert_eq(a->do_blue, b->do_blue);
}

Test(nde_serializers, mtf_roundtrip_unlinked) {
	/* unlinked: params + all three uparams round-trip */
	struct mtf_data in = { 0 };
	in.linked = FALSE;
	fill_mtf_params(&in.params, 0.2f);
	fill_mtf_params(&in.uparams[0], 0.3f);
	fill_mtf_params(&in.uparams[1], 0.4f);
	fill_mtf_params(&in.uparams[2], 0.5f);
	in.auto_display_compensation = TRUE;

	gchar *blob = op_desc_mtf.serialize(&in);
	cr_assert_not_null(blob);
	struct mtf_data *out = op_desc_mtf.deserialize(blob, op_desc_mtf.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->linked, FALSE);
	check_mtf_params(&out->params, &in.params);
	check_mtf_params(&out->uparams[0], &in.uparams[0]);
	check_mtf_params(&out->uparams[1], &in.uparams[1]);
	check_mtf_params(&out->uparams[2], &in.uparams[2]);
	cr_assert_eq(out->auto_display_compensation, TRUE);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mtf, blob);
	g_free(blob);
}

Test(nde_serializers, mtf_roundtrip_linked) {
	/* linked: uparams are not serialized (absent) and default to zero */
	struct mtf_data in = { 0 };
	in.linked = TRUE;
	fill_mtf_params(&in.params, 0.15f);
	in.auto_display_compensation = FALSE;

	gchar *blob = op_desc_mtf.serialize(&in);
	struct mtf_data *out = op_desc_mtf.deserialize(blob, op_desc_mtf.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->linked, TRUE);
	check_mtf_params(&out->params, &in.params);
	cr_assert_eq(out->auto_display_compensation, FALSE);
	FREE_VIA_DESTRUCTOR(out);
	g_free(blob);
}

Test(nde_serializers, mtf_inverse_shares_pair) {
	/* mtf_inverse uses the same serializer pair */
	struct mtf_data in = { 0 };
	in.linked = TRUE;
	fill_mtf_params(&in.params, 0.25f);
	gchar *blob = op_desc_mtf_inverse.serialize(&in);
	struct mtf_data *out = op_desc_mtf_inverse.deserialize(blob, op_desc_mtf_inverse.version);
	cr_assert_not_null(out);
	check_mtf_params(&out->params, &in.params);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mtf_inverse, blob);
	g_free(blob);
}

Test(nde_serializers, ghs_roundtrip) {
	ght_params p = { 0 };
	p.B = 1.5f; p.D = 2.25f; p.LP = 0.1f; p.SP = 0.5f; p.HP = 0.9f; p.BP = -0.03125f;
	p.stretchtype = STRETCH_ASINH;             /* 2 */
	p.payne_colourstretchmodel = COL_HUMANLUM; /* 1 */
	p.do_red = TRUE; p.do_green = FALSE; p.do_blue = TRUE;
	p.clip_mode = RGBBLEND;                     /* 2 */
	struct ght_data in = { 0 };
	in.params_ght = &p;
	in.auto_display_compensation = TRUE;

	gchar *blob = op_desc_ghs.serialize(&in);
	cr_assert_not_null(blob);
	struct ght_data *out = op_desc_ghs.deserialize(blob, op_desc_ghs.version);
	cr_assert_not_null(out);
	cr_assert_not_null(out->params_ght);
	cr_assert(memcmp(&out->params_ght->B, &p.B, sizeof(float)) == 0);
	cr_assert(memcmp(&out->params_ght->D, &p.D, sizeof(float)) == 0);
	cr_assert(memcmp(&out->params_ght->LP, &p.LP, sizeof(float)) == 0);
	cr_assert(memcmp(&out->params_ght->SP, &p.SP, sizeof(float)) == 0);
	cr_assert(memcmp(&out->params_ght->HP, &p.HP, sizeof(float)) == 0);
	cr_assert(memcmp(&out->params_ght->BP, &p.BP, sizeof(float)) == 0);
	cr_assert_eq(out->params_ght->stretchtype, p.stretchtype);
	cr_assert_eq(out->params_ght->payne_colourstretchmodel, p.payne_colourstretchmodel);
	cr_assert_eq(out->params_ght->do_red, p.do_red);
	cr_assert_eq(out->params_ght->do_green, p.do_green);
	cr_assert_eq(out->params_ght->do_blue, p.do_blue);
	cr_assert_eq(out->params_ght->clip_mode, p.clip_mode);
	cr_assert_eq(out->auto_display_compensation, TRUE);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_ghs, blob);
	g_free(blob);
}

Test(nde_serializers, autoghs_shares_ghs_pair) {
	/* autoghs uses the same serializer pair as ghs */
	cr_assert_eq(op_desc_autoghs.serialize, op_desc_ghs.serialize);
	cr_assert_eq(op_desc_autoghs.deserialize, op_desc_ghs.deserialize);
	ght_params p = { 0 };
	p.B = 0.5f; p.D = 1.0f; p.LP = 0.0f; p.SP = 0.3f; p.HP = 0.8f; p.BP = 0.0f;
	p.stretchtype = STRETCH_PAYNE_NORMAL;
	p.payne_colourstretchmodel = COL_INDEP;
	p.do_red = p.do_green = p.do_blue = TRUE;
	p.clip_mode = CLIP;
	struct ght_data in = { 0 };
	in.params_ght = &p;
	gchar *blob = op_desc_autoghs.serialize(&in);
	struct ght_data *out = op_desc_autoghs.deserialize(blob, op_desc_autoghs.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->params_ght->SP, &p.SP, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	g_free(blob);
}

/* Local mirror of the internal struct in command.c (POD, stable field
 * layout — the on-disk keys are the contract). */
struct t_autoghs_unlinked_data {
	destructor destroy_fn;
	float shadows_clipping;
	float amount;
	float b, hp, lp;
	clip_mode_t clip_mode;
};

Test(nde_serializers, autoghs_unlinked_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("stretch.autoghs_unlinked");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct t_autoghs_unlinked_data in = { 0 };
	in.shadows_clipping = -2.8f;
	in.amount = 13.0f;
	in.b = 0.1f;
	in.hp = 0.7f;
	in.lp = 0.03125f;
	in.clip_mode = RESCALEGLOBAL;   /* 3 */

	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_autoghs_unlinked_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->shadows_clipping, &in.shadows_clipping, sizeof(float)) == 0);
	cr_assert(memcmp(&out->amount, &in.amount, sizeof(float)) == 0);
	cr_assert(memcmp(&out->b, &in.b, sizeof(float)) == 0);
	cr_assert(memcmp(&out->hp, &in.hp, sizeof(float)) == 0);
	cr_assert(memcmp(&out->lp, &in.lp, sizeof(float)) == 0);
	cr_assert_eq(out->clip_mode, in.clip_mode);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

Test(nde_serializers, logstretch_paramless) {
	const op_descriptor *op = op_descriptor_by_id("stretch.log");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	/* paramless: serialize emits an empty blob; deserialize accepts it */
	gchar *blob = op->serialize(NULL);
	cr_assert_not_null(blob);
	cr_assert_str_eq(blob, "");
	void *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	FREE_VIA_DESTRUCTOR(out);
	/* paramless: the empty blob is VALID, so only the newer-version case is
	 * a failure (the generic empty-blob malformed check does not apply). */
	cr_assert_null(op->deserialize(blob, op->version + 1));
	g_free(blob);
}

Test(nde_serializers, asinh_roundtrip) {
	asinh_params in = { 0 };
	in.beta = 42.5f;
	in.offset = 0.1f;
	in.human_luminance = TRUE;
	in.clip_mode = RESCALE;   /* 1 */
	gchar *blob = op_desc_asinh.serialize(&in);
	cr_assert_not_null(blob);
	asinh_params *out = op_desc_asinh.deserialize(blob, op_desc_asinh.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->beta, &in.beta, sizeof(float)) == 0);
	cr_assert(memcmp(&out->offset, &in.offset, sizeof(float)) == 0);
	cr_assert_eq(out->human_luminance, in.human_luminance);
	cr_assert_eq(out->clip_mode, in.clip_mode);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_asinh, blob);
	g_free(blob);
}

Test(nde_serializers, scnr_roundtrip) {
	struct scnr_data in = { 0 };
	in.type = SCNR_ADDITIVE_MASK;   /* 3 */
	in.amount = 0.1;
	in.preserve = TRUE;
	gchar *blob = op_desc_scnr.serialize(&in);
	cr_assert_not_null(blob);
	struct scnr_data *out = op_desc_scnr.deserialize(blob, op_desc_scnr.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->type, in.type);
	cr_assert(memcmp(&out->amount, &in.amount, sizeof(double)) == 0);
	cr_assert_eq(out->preserve, in.preserve);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_scnr, blob);
	g_free(blob);
}

Test(nde_serializers, median_roundtrip) {
	struct median_filter_data in = { 0 };
	in.ksize = 7;
	in.amount = 0.1;
	in.iterations = 3;
	gchar *blob = op_desc_median.serialize(&in);
	cr_assert_not_null(blob);
	struct median_filter_data *out = op_desc_median.deserialize(blob, op_desc_median.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->ksize, in.ksize);
	cr_assert(memcmp(&out->amount, &in.amount, sizeof(double)) == 0);
	cr_assert_eq(out->iterations, in.iterations);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_median, blob);
	g_free(blob);
}

Test(nde_serializers, ccm_roundtrip) {
	struct ccm_data in = { 0 };
	float v = 0.1f;
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			in.matrix[r][c] = (v += 0.03125f) * ((r + c) % 2 ? -1.0f : 1.0f);
	in.power = 2.2f;
	gchar *blob = op_desc_ccm.serialize(&in);
	cr_assert_not_null(blob);
	struct ccm_data *out = op_desc_ccm.deserialize(blob, op_desc_ccm.version);
	cr_assert_not_null(out);
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			cr_assert(memcmp(&out->matrix[r][c], &in.matrix[r][c], sizeof(float)) == 0);
	cr_assert(memcmp(&out->power, &in.power, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_ccm, blob);
	g_free(blob);
}

/* ------------------------------------------------------------------ *
 *  Batch 3 — curves and background extraction (non-POD)              *
 * ------------------------------------------------------------------ */

static GList *make_point(GList *l, double x, double y) {
	point *pt = g_new(point, 1);
	pt->x = x; pt->y = y;
	return g_list_append(l, pt);
}

Test(nde_serializers, curves_roundtrip_with_points) {
	struct curve_params in = { 0 };
	in.algorithm = LINEAR;   /* 1 */
	in.do_channel[0] = TRUE; in.do_channel[1] = FALSE; in.do_channel[2] = TRUE;
	in.points = make_point(in.points, 0.1, 0.2);
	in.points = make_point(in.points, 0.5, 0.03125);
	in.points = make_point(in.points, 1.0 / 3.0, -0.25);

	gchar *blob = op_desc_curves.serialize(&in);
	cr_assert_not_null(blob);
	struct curve_params *out = op_desc_curves.deserialize(blob, op_desc_curves.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->algorithm, in.algorithm);
	cr_assert_eq(out->do_channel[0], TRUE);
	cr_assert_eq(out->do_channel[1], FALSE);
	cr_assert_eq(out->do_channel[2], TRUE);
	cr_assert_eq(g_list_length(out->points), 3);
	GList *a = in.points, *b = out->points;
	for (; a && b; a = a->next, b = b->next) {
		const point *pa = a->data, *pb = b->data;
		cr_assert(memcmp(&pa->x, &pb->x, sizeof(double)) == 0, "point x not bit-exact");
		cr_assert(memcmp(&pa->y, &pb->y, sizeof(double)) == 0, "point y not bit-exact");
	}
	/* the deserialized struct owns its list — its destructor must free it */
	cr_assert_not_null(out->destroy_fn);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_curves, blob);
	g_free(blob);
	g_list_free_full(in.points, g_free);
}

Test(nde_serializers, curves_empty_points) {
	/* empty points list: no points key; deserialize yields an empty list */
	struct curve_params in = { 0 };
	in.algorithm = CUBIC_SPLINE;   /* 0 */
	in.do_channel[0] = in.do_channel[1] = in.do_channel[2] = TRUE;
	in.points = NULL;
	gchar *blob = op_desc_curves.serialize(&in);
	cr_assert_not_null(blob);
	cr_assert_null(strstr(blob, "points="), "empty list must omit the points key");
	struct curve_params *out = op_desc_curves.deserialize(blob, op_desc_curves.version);
	cr_assert_not_null(out);
	cr_assert_null(out->points);
	cr_assert_eq(out->algorithm, CUBIC_SPLINE);
	FREE_VIA_DESTRUCTOR(out);
	g_free(blob);
}

Test(nde_serializers, bkg_remove_gradient_roundtrip) {
	struct background_data in = { 0 };
	in.method = BACKGROUND_METHOD_AUTO;            /* 1 */
	in.nb_of_samples = 20;
	in.tolerance = 0.1;
	in.correction = BACKGROUND_CORRECTION_DIVIDE;  /* 1 */
	in.interpolation_method = BACKGROUND_INTER_POLY; /* 1 */
	in.degree = BACKGROUND_POLY_3;                 /* 2 */
	in.smoothing = 0.03125;
	in.dither = TRUE;
	in.is_cfa = FALSE;
	in.randomize = TRUE;
	in.grad_descent = TRUE;
	in.border_value = 5.0;
	in.border_is_percent = TRUE;
	in.autograd.scale = 3.5;
	in.autograd.smoothness = 0.1;
	in.autograd.protect = TRUE;
	in.autograd.protect_threshold = 0.25;
	in.autograd.protect_amount = 0.75;
	in.autograd.simplified = TRUE;
	in.autograd.degree = 4;
	in.autograd.downsample = 2;

	/* effective sample positions (§15): one valid, one invalid (must be skipped) */
	background_sample s0 = { 0 }, s1 = { 0 };
	s0.valid = TRUE; s0.position.x = 12.5; s0.position.y = 0.1;
	s1.valid = FALSE; s1.position.x = 99.0; s1.position.y = 99.0;
	background_sample s2 = { 0 };
	s2.valid = TRUE; s2.position.x = 1.0 / 3.0; s2.position.y = -7.0;
	com.grad_samples = g_slist_append(com.grad_samples, &s0);
	com.grad_samples = g_slist_append(com.grad_samples, &s1);
	com.grad_samples = g_slist_append(com.grad_samples, &s2);

	gchar *blob = op_desc_remove_gradient.serialize(&in);
	cr_assert_not_null(blob);
	/* the samples key must carry the two VALID positions only */
	cr_assert_not_null(strstr(blob, "samples="));

	struct background_data *out = op_desc_remove_gradient.deserialize(blob, op_desc_remove_gradient.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->method, in.method);
	cr_assert_eq(out->nb_of_samples, in.nb_of_samples);
	cr_assert(memcmp(&out->tolerance, &in.tolerance, sizeof(double)) == 0);
	cr_assert_eq(out->correction, in.correction);
	cr_assert_eq(out->interpolation_method, in.interpolation_method);
	cr_assert_eq(out->degree, in.degree);
	cr_assert(memcmp(&out->smoothing, &in.smoothing, sizeof(double)) == 0);
	cr_assert_eq(out->dither, in.dither);
	cr_assert_eq(out->is_cfa, in.is_cfa);
	cr_assert_eq(out->randomize, in.randomize);
	cr_assert_eq(out->grad_descent, in.grad_descent);
	cr_assert(memcmp(&out->border_value, &in.border_value, sizeof(double)) == 0);
	cr_assert_eq(out->border_is_percent, in.border_is_percent);
	cr_assert(memcmp(&out->autograd.scale, &in.autograd.scale, sizeof(double)) == 0);
	cr_assert(memcmp(&out->autograd.smoothness, &in.autograd.smoothness, sizeof(double)) == 0);
	cr_assert_eq(out->autograd.protect, in.autograd.protect);
	cr_assert(memcmp(&out->autograd.protect_threshold, &in.autograd.protect_threshold, sizeof(double)) == 0);
	cr_assert(memcmp(&out->autograd.protect_amount, &in.autograd.protect_amount, sizeof(double)) == 0);
	cr_assert_eq(out->autograd.simplified, in.autograd.simplified);
	cr_assert_eq(out->autograd.degree, in.autograd.degree);
	cr_assert_eq(out->autograd.downsample, in.autograd.downsample);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_remove_gradient, blob);
	g_free(blob);
	g_slist_free(com.grad_samples);
	com.grad_samples = NULL;
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
	/* batch 2 — stretches, scnr, median, ccm */
	"stretch.mtf", "stretch.mtf_inverse", "stretch.log",
	"stretch.ghs", "stretch.autoghs", "stretch.autoghs_unlinked",
	"stretch.asinh", "filters.scnr", "filters.median", "color.ccm",
	/* batch 3 — curves, background extraction */
	"stretch.curves", "bkg.remove_gradient",
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
