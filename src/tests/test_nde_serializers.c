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
#include <unistd.h>
#include "core/siril.h"
#include "core/op_descriptor.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"
#include "core/nde_joint.h"

#include "algos/geometry.h"
#include "filters/mtf.h"
#include "filters/ght.h"
#include "filters/asinh.h"
#include "filters/scnr.h"
#include "filters/saturation.h"
#include "filters/median.h"
#include "algos/colors.h"
#include "filters/curve_transform.h"
#include "algos/background_extraction.h"
#include "filters/banding.h"
#include "filters/clahe.h"
#include "filters/rgradient.h"
#include "filters/wavelets.h"
#include "algos/wavelet_denoise.h"
#include "filters/cosmetic_correction.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "filters/linear_match.h"
#include "filters/epf.h"
#include "filters/synthstar.h"
#include "filters/unpurple.h"
#include "algos/PSF.h"
#include "algos/photometric_cc.h"
#include "compositing/align_rgb.h"
#include "io/image_format_fits.h"
#include "core/processing.h"
#include "filters/deconvolution/deconvolution.h"
#include "core/icc_profile.h"
#include "core/masks.h"

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
	/* post-phase-3 additions (maintainer-requested coverage) */
	"color.saturation",
	/* Opus batch — remaining mechanically-serializable ops */
	"arith.thresh", "arith.nozero", "arith.ffill", "arith.fill",
	"arith.offset", "arith.limit", "arith.neg", "arith.fmul",
	"filters.gauss", "filters.unsharp", "filters.ddp", "filters.denoise",
	"color.grey_flat", "cfa.fix_xtrans",
	"filters.banding", "filters.clahe", "filters.cosmetic",
	"filters.rgradient", "wavelets.atrous",
	/* phase 4.5 batch 1 — file operands (Convention 1) */
	"arith.imoper", "arith.addmax", "arith.fdiv",
	"color.linear_match", "filters.cosme", "filters.epf",
	/* phase 4.5 batch 2 — star lists (Conv. 2), PCC/SPCC (Conv. 3), ICC (Conv. 4).
	 * fft + wrecons stay permanently Tier B by design. */
	"star.synthstar", "filters.unpurple", "color.photometric_cc", "icc.convert",
	/* star provenance refinement: unclip promoted to Tier A (DELEGATED conf) */
	"star.unclip",
	/* graph step 4.2 — the mask.* ops (core/masks.c).  All 13 live mask ops,
	 * including the four param-less ones: a Tier-B mask record would be a
	 * barrier in its chain.  mask.mtf_autostretch is an image op with no live
	 * call sites and is deliberately excluded. */
	"mask.from_stars", "mask.from_channel", "mask.from_luminance",
	"mask.clear", "mask.threshold", "mask.blur", "mask.feather",
	"mask.multiply", "mask.invert", "mask.autostretch", "mask.bitpix",
	"mask.from_color", "mask.from_gradient",
	/* The last of the direct-apply colour tools.  Background neutralization
	 * and colour calibration were never routed through generic_image_worker,
	 * so they had no descriptor and nowhere to hang a serializer; they have
	 * one now, the calibration recording its computed coefficients per
	 * Convention 3.  RGB alignment records its method and the effective
	 * registration area and re-runs the detection at replay — the delegated
	 * provenance star.unclip settled on. */
	"color.background_neutralization", "color.calibration", "color.rgb_align",
	/* Convention 2 applied to a PSF: the kernel the deconvolution actually
	 * used is stashed and replayed, never re-estimated. */
	"filters.deconvolve",
	/* Per-layer factorisation of the multi-layer colour operations (layers
	 * match, group colour calibration): one record per affected layer. */
	"flis.layer_scale",
	/* JOINT records (nde_joint.h): ONE record per multi-layer operation,
	 * member of every participant's chain, recomputing its analysis at
	 * replay.  Supersedes the per-layer flis.layer_scale factorisation for
	 * new captures; layer_scale stays for legacy files and manual scaling. */
	"flis.layers_match", "flis.group_calibration",
	/* The GEOMETRIC joint record: multi-layer registration.  It stores the
	 * solved transforms AND the settings that produced them and picks between
	 * them at replay (nde_joint.h), which is why it is the one joint op that
	 * is also OP_GEOMETRY_CHANGING. */
	"flis.register",
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

Test(nde_serializers, layers_match_roundtrip) {
	struct nde_joint_layers_match_data *in = nde_joint_layers_match_data_new(3);
	int items[3] = { 4, 7, 11 };
	const char *names[3] = { "Ha", "OIII", "SII" };
	for (guint i = 0; i < 3; i++) {
		in->parts[i].item_id = items[i];
		in->parts[i].name = g_strdup(names[i]);
		in->parts[i].tinted = i != 1;
		in->parts[i].tint[0] = 0.83 + i * 0.01;
		in->parts[i].tint[1] = 0.25;
		in->parts[i].tint[2] = 0.0625;
		in->parts[i].diag_scale = 1.0 / 3.0 + i;
		in->parts[i].diag_offset = i == 2 ? -0.0043 : 0.0;
	}
	gchar *blob = op_desc_flis_layers_match.serialize(in);
	cr_assert_not_null(blob);
	struct nde_joint_layers_match_data *out =
		op_desc_flis_layers_match.deserialize(blob, op_desc_flis_layers_match.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->n, 3);
	for (guint i = 0; i < 3; i++) {
		cr_assert_eq(out->parts[i].item_id, items[i]);
		cr_assert_str_eq(out->parts[i].name, names[i]);
		cr_assert_eq(out->parts[i].tinted, in->parts[i].tinted);
		if (out->parts[i].tinted) {
			cr_assert(memcmp(out->parts[i].tint, in->parts[i].tint,
			                 3 * sizeof(double)) == 0);
		}
		cr_assert(memcmp(&out->parts[i].diag_scale, &in->parts[i].diag_scale,
		                 sizeof(double)) == 0);
		cr_assert(memcmp(&out->parts[i].diag_offset, &in->parts[i].diag_offset,
		                 sizeof(double)) == 0);
	}
	cr_assert_not_null(out->destroy_fn);
	out->destroy_fn(out);

	/* Newer-version blobs must refuse. */
	cr_assert_null(op_desc_flis_layers_match.deserialize(blob,
	               op_desc_flis_layers_match.version + 1));
	g_free(blob);
	nde_joint_layers_match_data_free(in);
}

Test(nde_serializers, layers_match_rejects_bad_blobs) {
	const op_descriptor *op = &op_desc_flis_layers_match;
	/* One participant is not a joint operation. */
	cr_assert_null(op->deserialize("n=1;i0_item=3;i0_name=x;i0_tint=0;i0_s=1", op->version));
	/* A layer listed twice would be scaled twice. */
	cr_assert_null(op->deserialize(
		"n=2;i0_item=3;i0_name=a;i0_tint=0;i0_s=1;"
		"i1_item=3;i1_name=b;i1_tint=0;i1_s=1", op->version));
	/* Tinted without colour components. */
	cr_assert_null(op->deserialize(
		"n=2;i0_item=3;i0_name=a;i0_tint=1;i0_s=1;"
		"i1_item=4;i1_name=b;i1_tint=0;i1_s=1", op->version));
	/* Negative diagnostic scale. */
	cr_assert_null(op->deserialize(
		"n=2;i0_item=3;i0_name=a;i0_tint=0;i0_s=-2;"
		"i1_item=4;i1_name=b;i1_tint=0;i1_s=1", op->version));
	/* Missing participant entry. */
	cr_assert_null(op->deserialize("n=2;i0_item=3;i0_name=a;i0_tint=0;i0_s=1", op->version));
}

Test(nde_serializers, saturation_roundtrip) {
	saturation_params in = { 0 };
	in.coeff = 0.35;
	in.background_factor = 1.1;
	in.h_min = 346.0;
	in.h_max = 20.0;
	gchar *blob = op_desc_saturation.serialize(&in);
	cr_assert_not_null(blob);
	saturation_params *out = op_desc_saturation.deserialize(blob, op_desc_saturation.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->coeff, &in.coeff, sizeof(double)) == 0);
	cr_assert(memcmp(&out->background_factor, &in.background_factor, sizeof(double)) == 0);
	cr_assert(memcmp(&out->h_min, &in.h_min, sizeof(double)) == 0);
	cr_assert(memcmp(&out->h_max, &in.h_max, sizeof(double)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_saturation, blob);
	g_free(blob);
}

/* ------------------------------------------------------------------ *
 *  Opus batch — remaining mechanically-serializable ops.             *
 *                                                                    *
 *  The arith/gauss/unsharp/ddp/denoise/fmul structs are file-local   *
 *  to command.c, so (as with autoghs_unlinked above) we mirror their *
 *  POD layout here — the on-disk keys are the contract, not the      *
 *  struct symbol.  The filter-module ops (banding/clahe/cosmetic/    *
 *  rgradient/atrous) have public headers and use the real structs.   *
 * ------------------------------------------------------------------ */

/* command.c: struct thresh_data { destructor; int type; int lo; int hi; } */
struct t_thresh_data { destructor destroy_fn; int type; int lo; int hi; };
Test(nde_serializers, thresh_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.thresh");
	cr_assert_not_null(op);
	struct t_thresh_data in = { 0 };
	in.type = 2; in.lo = -3; in.hi = 41234;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_thresh_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->type, in.type);
	cr_assert_eq(out->lo, in.lo);
	cr_assert_eq(out->hi, in.hi);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct nozero_data { destructor; WORD level; } */
struct t_nozero_data { destructor destroy_fn; guint16 level; };
Test(nde_serializers, nozero_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.nozero");
	cr_assert_not_null(op);
	struct t_nozero_data in = { 0 };
	in.level = 40000;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_nozero_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->level, in.level);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct ffill_data / fill_data { destructor; int level; rectangle area; } */
struct t_fill_data { destructor destroy_fn; int level; rectangle area; };
Test(nde_serializers, ffill_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.ffill");
	cr_assert_not_null(op);
	struct t_fill_data in = { 0 };
	in.level = 137;
	in.area.x = 5; in.area.y = -9; in.area.w = 320; in.area.h = 200;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_fill_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->level, in.level);
	cr_assert_eq(out->area.x, in.area.x);
	cr_assert_eq(out->area.y, in.area.y);
	cr_assert_eq(out->area.w, in.area.w);
	cr_assert_eq(out->area.h, in.area.h);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

Test(nde_serializers, fill_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.fill");
	cr_assert_not_null(op);
	struct t_fill_data in = { 0 };
	in.level = 65000;
	in.area.x = 1; in.area.y = 2; in.area.w = 3; in.area.h = 4;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_fill_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->level, in.level);
	cr_assert_eq(out->area.w, in.area.w);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct offset_data { destructor; float level; } */
struct t_offset_data { destructor destroy_fn; float level; };
Test(nde_serializers, offset_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.offset");
	cr_assert_not_null(op);
	struct t_offset_data in = { 0 };
	in.level = -12.5f + 0.1f;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_offset_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->level, &in.level, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct limit_data { destroy_fn; OverrangeResponse method; } */
struct t_limit_data { destructor destroy_fn; OverrangeResponse method; };
Test(nde_serializers, limit_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.limit");
	cr_assert_not_null(op);
	struct t_limit_data in = { 0 };
	in.method = RESPONSE_RESCALE_ALL;   /* 4 */
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_limit_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->method, in.method);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* Paramless ops (arith.neg, color.grey_flat, cfa.fix_xtrans, stretch.log, and
 * the four mask.* param-less ops): empty blob is VALID; only a newer version
 * fails.  CHECK_MALFORMED must NOT be used on these. */
static void check_paramless(const char *id) {
	const op_descriptor *op = op_descriptor_by_id(id);
	cr_assert_not_null(op, "%s missing", id);
	cr_assert_not_null(op->serialize, "%s has no serializer", id);
	gchar *blob = op->serialize(NULL);
	cr_assert_not_null(blob);
	cr_assert_str_eq(blob, "");
	void *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	FREE_VIA_DESTRUCTOR(out);
	cr_assert_null(op->deserialize(blob, op->version + 1));
	g_free(blob);
}

Test(nde_serializers, neg_paramless)        { check_paramless("arith.neg"); }
Test(nde_serializers, grey_flat_paramless)  { check_paramless("color.grey_flat"); }
Test(nde_serializers, fix_xtrans_paramless) { check_paramless("cfa.fix_xtrans"); }

/* command.c: struct fmul_data { destructor; float coeff; gboolean from8b; }
 * from8b is intentionally NOT serialized (runtime UI-only). */
struct t_fmul_data { destructor destroy_fn; float coeff; gboolean from8b; };
Test(nde_serializers, fmul_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.fmul");
	cr_assert_not_null(op);
	struct t_fmul_data in = { 0 };
	in.coeff = 2.5f + 0.1f;
	in.from8b = TRUE;   /* must not affect the blob */
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	cr_assert_null(strstr(blob, "from8b"), "from8b must not be serialized");
	struct t_fmul_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->coeff, &in.coeff, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct gauss_data { destructor; double sigma; } */
struct t_gauss_data { destructor destroy_fn; double sigma; };
Test(nde_serializers, gauss_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("filters.gauss");
	cr_assert_not_null(op);
	struct t_gauss_data in = { 0 };
	in.sigma = 3.0 + 0.1;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_gauss_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->sigma, &in.sigma, sizeof(double)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct unsharp_data { destructor; double sigma; double multi; } */
struct t_unsharp_data { destructor destroy_fn; double sigma; double multi; };
Test(nde_serializers, unsharp_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("filters.unsharp");
	cr_assert_not_null(op);
	struct t_unsharp_data in = { 0 };
	in.sigma = 5.0 + 0.1; in.multi = -0.5 + 0.03125;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_unsharp_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->sigma, &in.sigma, sizeof(double)) == 0);
	cr_assert(memcmp(&out->multi, &in.multi, sizeof(double)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

/* command.c: struct ddp_data { destructor; float level; float coeff; float sigma; } */
struct t_ddp_data { destructor destroy_fn; float level; float coeff; float sigma; };
Test(nde_serializers, ddp_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("filters.ddp");
	cr_assert_not_null(op);
	struct t_ddp_data in = { 0 };
	in.level = 100.0f + 0.1f; in.coeff = 1.5f + 0.03125f; in.sigma = 3.0f + 0.1f;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_ddp_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->level, &in.level, sizeof(float)) == 0);
	cr_assert(memcmp(&out->coeff, &in.coeff, sizeof(float)) == 0);
	cr_assert(memcmp(&out->sigma, &in.sigma, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

Test(nde_serializers, denoise_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("filters.denoise");
	cr_assert_not_null(op);
	denoise_args in = { 0 };
	in.modulation = 0.75f + 0.03125f;
	in.sos = 4;
	in.da3d = 1;
	in.rho = 0.2f + 0.1f;
	in.do_anscombe = TRUE;
	in.do_cosme = FALSE;
	in.suppress_artefacts = TRUE;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	denoise_args *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->modulation, &in.modulation, sizeof(float)) == 0);
	cr_assert_eq(out->sos, in.sos);
	cr_assert_eq(out->da3d, in.da3d);
	cr_assert(memcmp(&out->rho, &in.rho, sizeof(float)) == 0);
	cr_assert_eq(out->do_anscombe, in.do_anscombe);
	cr_assert_eq(out->do_cosme, in.do_cosme);
	cr_assert_eq(out->suppress_artefacts, in.suppress_artefacts);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	g_free(blob);
}

Test(nde_serializers, banding_roundtrip) {
	struct banding_data in = { 0 };
	in.sigma = 2.0 + 0.1;
	in.amount = 0.5 + 0.03125;
	in.protect_highlights = TRUE;
	in.applyRotation = FALSE;
	gchar *blob = op_desc_banding.serialize(&in);
	cr_assert_not_null(blob);
	struct banding_data *out = op_desc_banding.deserialize(blob, op_desc_banding.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->sigma, &in.sigma, sizeof(double)) == 0);
	cr_assert(memcmp(&out->amount, &in.amount, sizeof(double)) == 0);
	cr_assert_eq(out->protect_highlights, in.protect_highlights);
	cr_assert_eq(out->applyRotation, in.applyRotation);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_banding, blob);
	g_free(blob);
}

Test(nde_serializers, clahe_roundtrip) {
	clahe_params in = { 0 };
	in.clip = 2.0 + 0.1;
	in.tileSize = 16;
	gchar *blob = op_desc_clahe.serialize(&in);
	cr_assert_not_null(blob);
	clahe_params *out = op_desc_clahe.deserialize(blob, op_desc_clahe.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->clip, &in.clip, sizeof(double)) == 0);
	cr_assert_eq(out->tileSize, in.tileSize);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_clahe, blob);
	g_free(blob);
}

Test(nde_serializers, cosmetic_roundtrip) {
	struct cosmetic_data in = { 0 };
	in.sigma[0] = 3.0 + 0.1; in.sigma[1] = 2.5 + 0.03125;
	in.amount = 0.8 + 0.1;
	in.is_cfa = TRUE;
	gchar *blob = op_desc_cosmetic.serialize(&in);
	cr_assert_not_null(blob);
	struct cosmetic_data *out = op_desc_cosmetic.deserialize(blob, op_desc_cosmetic.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->sigma[0], &in.sigma[0], sizeof(double)) == 0);
	cr_assert(memcmp(&out->sigma[1], &in.sigma[1], sizeof(double)) == 0);
	cr_assert(memcmp(&out->amount, &in.amount, sizeof(double)) == 0);
	cr_assert_eq(out->is_cfa, in.is_cfa);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_cosmetic, blob);
	g_free(blob);
}

Test(nde_serializers, rgradient_roundtrip) {
	struct rgradient_data in = { 0 };
	in.xc = 512.5 + 0.1; in.yc = 384.0 + 0.03125;
	in.dR = 3.0 + 0.1; in.da = -15.0 + 0.03125;
	gchar *blob = op_desc_rgradient.serialize(&in);
	cr_assert_not_null(blob);
	struct rgradient_data *out = op_desc_rgradient.deserialize(blob, op_desc_rgradient.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->xc, &in.xc, sizeof(double)) == 0);
	cr_assert(memcmp(&out->yc, &in.yc, sizeof(double)) == 0);
	cr_assert(memcmp(&out->dR, &in.dR, sizeof(double)) == 0);
	cr_assert(memcmp(&out->da, &in.da, sizeof(double)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_rgradient, blob);
	g_free(blob);
}

Test(nde_serializers, atrous_roundtrip) {
	struct atrous_data in = { 0 };
	in.nbr_plan = 5;
	in.type = 1;
	in.anscombe = TRUE;
	for (int i = 0; i < 7; i++)
		in.coef[i] = 0.5f + (float)i * 0.03125f;
	in.denoise.enabled = TRUE;
	in.denoise.method = 2;
	in.denoise.k = 3.0f + 0.1f;
	for (int i = 0; i < WD_MAX_PLAN; i++)
		in.denoise.f[i] = 1.0f + (float)i * 0.03125f;
	in.denoise.sigma_source = 1;
	in.denoise.soft = TRUE;
	in.denoise.anscombe = FALSE;
	gchar *blob = op_desc_atrous.serialize(&in);
	cr_assert_not_null(blob);
	struct atrous_data *out = op_desc_atrous.deserialize(blob, op_desc_atrous.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->nbr_plan, in.nbr_plan);
	cr_assert_eq(out->type, in.type);
	cr_assert_eq(out->anscombe, in.anscombe);
	for (int i = 0; i < 7; i++)
		cr_assert(memcmp(&out->coef[i], &in.coef[i], sizeof(float)) == 0);
	cr_assert_eq(out->denoise.enabled, in.denoise.enabled);
	cr_assert_eq(out->denoise.method, in.denoise.method);
	cr_assert(memcmp(&out->denoise.k, &in.denoise.k, sizeof(float)) == 0);
	for (int i = 0; i < WD_MAX_PLAN; i++)
		cr_assert(memcmp(&out->denoise.f[i], &in.denoise.f[i], sizeof(float)) == 0);
	cr_assert_eq(out->denoise.sigma_source, in.denoise.sigma_source);
	cr_assert_eq(out->denoise.soft, in.denoise.soft);
	cr_assert_eq(out->denoise.anscombe, in.denoise.anscombe);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_atrous, blob);
	g_free(blob);
}

/* ------------------------------------------------------------------ *
 *  phase 4.5 batch 1 — file operands (Convention 1).                 *
 *                                                                    *
 *  These serializers hash the operand FILE at capture, so each test *
 *  writes a small temp file and asserts operand_path/_sha256/_size   *
 *  round-trip alongside the POD.  The arith structs are file-local   *
 *  to command.c, so (as elsewhere) we mirror their layout — the      *
 *  on-disk keys are the contract, not the struct symbol.             *
 * ------------------------------------------------------------------ */

/* Write @content to a fresh temp file; return the heap path (caller g_free). */
static gchar *write_temp_operand(const char *content) {
	GError *e = NULL;
	gchar *path = NULL;
	gint fd = g_file_open_tmp("nde_operand_XXXXXX", &path, &e);
	cr_assert(fd >= 0, "could not create temp operand: %s", e ? e->message : "?");
	close(fd);
	cr_assert(g_file_set_contents(path, content, -1, &e),
	          "could not write temp operand: %s", e ? e->message : "?");
	return path;
}

/* command.c: struct imoper_data { destructor; image_operator oper;
 *                                 char *filename; gboolean force_to_float; } */
struct t_imoper_data { destructor destroy_fn; int oper; char *filename; gboolean ftf; };
Test(nde_serializers, imoper_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.imoper");
	cr_assert_not_null(op);
	gchar *path = write_temp_operand("imoper operand bytes\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	struct t_imoper_data in = { 0 };
	in.oper = 2;            /* OPER_MUL */
	in.filename = path;
	in.ftf = TRUE;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_imoper_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->oper, in.oper);
	cr_assert_eq(out->ftf, in.ftf);
	cr_assert_str_eq(out->filename, path);
	/* the blob pins the file hash+size computed at capture */
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	gint64 blob_size = 0;
	cr_assert(nde_kv_get_int(kv, "operand_size", &blob_size));
	cr_assert_eq(blob_size, size);
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* a record without the operand keys must not deserialize */
	cr_assert_null(op->deserialize("oper=0;force_to_float=1", op->version),
	               "missing operand_path must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

/* command.c: struct addmax_data { destructor; fits *operand_fit;
 *                                 char *operand_path; gboolean force_to_float; } */
struct t_addmax_data { destructor destroy_fn; void *operand_fit; char *operand_path; gboolean ftf; };
Test(nde_serializers, addmax_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.addmax");
	cr_assert_not_null(op);
	gchar *path = write_temp_operand("addmax operand bytes 12345\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	struct t_addmax_data in = { 0 };
	in.operand_path = path;
	in.ftf = FALSE;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_addmax_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->ftf, in.ftf);
	cr_assert_str_eq(out->operand_path, path);
	cr_assert_null(out->operand_fit, "deserialize must not touch the filesystem");
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	cr_assert_null(op->deserialize("force_to_float=0", op->version),
	               "missing operand_path must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

/* command.c: struct fdiv_data { destructor; fits *operand_fit; float norm;
 *                               char *operand_path; gboolean force_to_float; } */
struct t_fdiv_data { destructor destroy_fn; void *operand_fit; float norm; char *operand_path; gboolean ftf; };
Test(nde_serializers, fdiv_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("arith.fdiv");
	cr_assert_not_null(op);
	gchar *path = write_temp_operand("fdiv operand bytes abcdef\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	struct t_fdiv_data in = { 0 };
	in.norm = 65535.0f + 0.03125f;
	in.operand_path = path;
	in.ftf = TRUE;
	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct t_fdiv_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->norm, &in.norm, sizeof(float)) == 0);
	cr_assert_eq(out->ftf, in.ftf);
	cr_assert_str_eq(out->operand_path, path);
	cr_assert_null(out->operand_fit, "deserialize must not touch the filesystem");
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	cr_assert_null(op->deserialize("norm=1;force_to_float=1", op->version),
	               "missing operand_path must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

Test(nde_serializers, linear_match_roundtrip) {
	gchar *path = write_temp_operand("linear match ref bytes\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	struct linear_match_data in = { 0 };
	in.low = 0.05 + 0.03125; in.high = 0.95 + 0.03125;
	in.force_to_float = TRUE;
	in.operand_path = path;
	gchar *blob = op_desc_linear_match.serialize(&in);
	cr_assert_not_null(blob);
	struct linear_match_data *out = op_desc_linear_match.deserialize(blob, op_desc_linear_match.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->low, &in.low, sizeof(double)) == 0);
	cr_assert(memcmp(&out->high, &in.high, sizeof(double)) == 0);
	cr_assert_eq(out->force_to_float, in.force_to_float);
	cr_assert_str_eq(out->operand_path, path);
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_linear_match, blob);
	cr_assert_null(op_desc_linear_match.deserialize("low=0;high=1;force_to_float=1", op_desc_linear_match.version),
	               "missing operand_path must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

Test(nde_serializers, cosme_roundtrip) {
	gchar *path = write_temp_operand("P 10 20 H\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	struct cosme_data in = { 0 };
	in.is_cfa = 1;
	in.file = g_file_new_for_path(path);
	gchar *blob = op_desc_cosme.serialize(&in);
	cr_assert_not_null(blob);
	g_object_unref(in.file);
	struct cosme_data *out = op_desc_cosme.deserialize(blob, op_desc_cosme.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->is_cfa, in.is_cfa);
	cr_assert_null(out->file, "deserialize must not touch the filesystem");
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_path"), path);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_cosme, blob);
	cr_assert_null(op_desc_cosme.deserialize("is_cfa=0", op_desc_cosme.version),
	               "missing operand_path must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

/* EPF: bilateral serializes plain POD (no operand keys); self-guided emits
 * guide=self; a file guide pins the guide image. */
Test(nde_serializers, epf_bilateral_roundtrip) {
	struct epfargs in = { 0 };
	in.d = 5.0 + 0.03125; in.sigma_col = 11.0 + 0.1;
	in.sigma_space = 8.0 + 0.03125; in.mod = 0.75 + 0.03125;
	in.filter = EP_BILATERAL;
	in.guidefit = NULL;
	gchar *blob = op_desc_epf.serialize(&in);
	cr_assert_not_null(blob);
	struct epfargs *out = op_desc_epf.deserialize(blob, op_desc_epf.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->d, &in.d, sizeof(double)) == 0);
	cr_assert(memcmp(&out->sigma_col, &in.sigma_col, sizeof(double)) == 0);
	cr_assert(memcmp(&out->sigma_space, &in.sigma_space, sizeof(double)) == 0);
	cr_assert(memcmp(&out->mod, &in.mod, sizeof(double)) == 0);
	cr_assert_eq(out->filter, in.filter);
	/* no operand keys for a non-guided filter */
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_null(nde_kv_get_str(kv, "guide"));
	cr_assert_null(nde_kv_get_str(kv, "operand_path"));
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_epf, blob);
	g_free(blob);
}

Test(nde_serializers, epf_self_guide_roundtrip) {
	fits target = { 0 };
	struct epfargs in = { 0 };
	in.d = 5.0; in.sigma_col = 11.0; in.sigma_space = 8.0; in.mod = 1.0;
	in.filter = EP_GUIDED;
	in.fit = &target;
	in.guidefit = &target;    /* self-guide: guidefit == fit */
	gchar *blob = op_desc_epf.serialize(&in);
	cr_assert_not_null(blob);
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "guide"), "self");
	cr_assert_null(nde_kv_get_str(kv, "operand_path"), "self-guide has no file operand");
	g_hash_table_unref(kv);
	struct epfargs *out = op_desc_epf.deserialize(blob, op_desc_epf.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->filter, EP_GUIDED);
	cr_assert_null(out->guide_path);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_epf, blob);
	g_free(blob);
}

Test(nde_serializers, epf_file_guide_roundtrip) {
	gchar *path = write_temp_operand("epf guide image bytes\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	fits target = { 0 }, guide = { 0 };
	struct epfargs in = { 0 };
	in.d = 5.0; in.sigma_col = 11.0; in.sigma_space = 8.0; in.mod = 1.0;
	in.filter = EP_GUIDED;
	in.fit = &target;
	in.guidefit = &guide;      /* separate guide → file operand */
	in.guide_path = path;
	gchar *blob = op_desc_epf.serialize(&in);
	cr_assert_not_null(blob);
	struct epfargs *out = op_desc_epf.deserialize(blob, op_desc_epf.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->filter, EP_GUIDED);
	cr_assert_str_eq(out->guide_path, path);
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "guide"), "file");
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_path"), path);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	g_hash_table_unref(kv);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_epf, blob);
	/* a guided record missing its guide token is not replayable */
	cr_assert_null(op_desc_epf.deserialize("d=5;sigma_col=11;sigma_space=8;mod=1;filter=1",
	               op_desc_epf.version), "guided record without guide token must yield NULL");
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

/* ==================================================================== *
 *  phase 4.5 batch 2 — Conventions 2 (stars), 3 (PCC/SPCC), 4 (ICC).   *
 * ==================================================================== */

/* ---- Convention 2: star-list blob codec + synthstar/unpurple ---- */

/* Build a small manufactured psf_star array with distinctive fields. */
static psf_star **make_test_stars(int n) {
	psf_star **stars = malloc((size_t)(n + 1) * sizeof(psf_star *));
	for (int i = 0; i < n; i++) {
		psf_star *s = new_psf_star();
		s->xpos = 10.5 + i;            /* exactly representable at %.9g */
		s->ypos = 20.25 + 2 * i;
		s->A = 0.75 + 0.03125 * i;
		s->fwhmx = 3.5 + 0.5 * i;
		s->fwhmy = 2.75 + 0.25 * i;
		s->beta = 4.0 + i;
		s->has_saturated = (i % 2) ? TRUE : FALSE;
		s->profile = (i % 2) ? PSF_MOFFAT_BFREE : PSF_GAUSSIAN;
		stars[i] = s;
	}
	stars[n] = NULL;
	return stars;
}

Test(nde_serializers, star_blob_roundtrip) {
	const int n = 4;
	psf_star **in = make_test_stars(n);
	gchar *blob = synthstar_stars_to_blob(in, n);
	cr_assert_not_null(blob);

	int nb = 0;
	psf_star **out = synthstar_stars_from_blob(blob, &nb);
	cr_assert_not_null(out);
	cr_assert_eq(nb, n, "star count mismatch: %d vs %d", nb, n);
	for (int i = 0; i < n; i++) {
		/* the stashed fields must survive with %.9g precision (exactly
		 * representable test values → bit-exact) */
		cr_assert_float_eq(out[i]->xpos,  in[i]->xpos,  0.0, "xpos %d", i);
		cr_assert_float_eq(out[i]->ypos,  in[i]->ypos,  0.0, "ypos %d", i);
		cr_assert_float_eq(out[i]->A,     in[i]->A,     0.0, "A %d", i);
		cr_assert_float_eq(out[i]->fwhmx, in[i]->fwhmx, 0.0, "fwhmx %d", i);
		cr_assert_float_eq(out[i]->fwhmy, in[i]->fwhmy, 0.0, "fwhmy %d", i);
		cr_assert_float_eq(out[i]->beta,  in[i]->beta,  0.0, "beta %d", i);
		cr_assert_eq(out[i]->has_saturated, in[i]->has_saturated, "sat %d", i);
		cr_assert_eq((int)out[i]->profile, (int)in[i]->profile, "profile %d", i);
	}
	/* re-serialization must reproduce the identical on-disk string (the
	 * on-disk representation is the contract) */
	gchar *blob2 = synthstar_stars_to_blob(out, nb);
	cr_assert_str_eq(blob2, blob);

	free_fitted_stars(in);
	free_fitted_stars(out);
	g_free(blob); g_free(blob2);

	/* empty / NULL cases */
	cr_assert_null(synthstar_stars_to_blob(NULL, 0));
	cr_assert_null(synthstar_stars_from_blob("", &nb));
	cr_assert_eq(nb, 0);
}

Test(nde_serializers, synthstar_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("star.synthstar");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct synthstar_data *in = new_synthstar_data();
	psf_star **stars = make_test_stars(3);
	in->stars_blob = synthstar_stars_to_blob(stars, 3);
	free_fitted_stars(stars);

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct synthstar_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_str_eq(out->stars_blob, in->stars_blob);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* a record with no stars= key predates the batch → non-replayable */
	cr_assert_null(op->deserialize("", op->version),
	               "synthstar without a stars key must yield NULL");
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, unpurple_roundtrip_pod) {
	const op_descriptor *op = op_descriptor_by_id("filters.unpurple");
	cr_assert_not_null(op);
	struct unpurpleargs *in = new_unpurple_args();
	in->mod_b = 0.375;
	in->thresh = 0.125;
	in->withstarmask = FALSE;   /* POD variant: no external input */
	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct unpurpleargs *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_float_eq(out->mod_b, in->mod_b, 0.0);
	cr_assert_float_eq(out->thresh, in->thresh, 0.0);
	cr_assert_eq(out->withstarmask, FALSE);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, unpurple_roundtrip_starmask) {
	const op_descriptor *op = op_descriptor_by_id("filters.unpurple");
	cr_assert_not_null(op);
	struct unpurpleargs *in = new_unpurple_args();
	in->mod_b = 0.5;
	in->thresh = 0.25;
	in->withstarmask = TRUE;
	psf_star **stars = make_test_stars(2);
	in->stars_blob = synthstar_stars_to_blob(stars, 2);
	free_fitted_stars(stars);

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct unpurpleargs *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->withstarmask, TRUE);
	cr_assert_str_eq(out->stars_blob, in->stars_blob);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* a withstarmask record missing its stars= key is not replayable */
	cr_assert_null(op->deserialize("mod_b=0.5;thresh=0.25;withstarmask=1",
	               op->version),
	               "withstarmask record without stars must yield NULL");
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

/* ---- DELEGATED (auto-detect) provenance: detection-parameter codec ---- */

/* A star_finder_params filled with distinct, recognisable values. */
static star_finder_params make_test_sf_conf(void) {
	star_finder_params sf = { 0 };
	sf.radius = 13;
	sf.sigma = 0.75;
	sf.roundness = 0.42;
	sf.focal_length = 530.0;
	sf.pixel_size_x = 3.76;
	sf.convergence = 2;
	sf.relax_checks = TRUE;
	sf.profile = PSF_MOFFAT_BFREE;
	sf.min_beta = 1.5;
	sf.min_A = 0.001;
	sf.max_A = 0.9;
	sf.max_r = 1.0;
	return sf;
}

static void assert_sf_conf_eq(const star_finder_params *a, const star_finder_params *b) {
	cr_assert_eq(a->radius, b->radius);
	cr_assert_float_eq(a->sigma, b->sigma, 0.0);
	cr_assert_float_eq(a->roundness, b->roundness, 0.0);
	cr_assert_float_eq(a->focal_length, b->focal_length, 0.0);
	cr_assert_float_eq(a->pixel_size_x, b->pixel_size_x, 0.0);
	cr_assert_eq(a->convergence, b->convergence);
	cr_assert_eq(a->relax_checks, b->relax_checks);
	cr_assert_eq(a->profile, b->profile);
	cr_assert_float_eq(a->min_beta, b->min_beta, 0.0);
	cr_assert_float_eq(a->min_A, b->min_A, 0.0);
	cr_assert_float_eq(a->max_A, b->max_A, 0.0);
	cr_assert_float_eq(a->max_r, b->max_r, 0.0);
}

Test(nde_serializers, starfinder_conf_codec_roundtrip) {
	star_finder_params in = make_test_sf_conf();
	GString *kv = nde_kv_start();
	synthstar_conf_to_kv(kv, &in);
	gchar *blob = nde_kv_end(kv);
	cr_assert_not_null(blob);
	GHashTable *parsed = nde_kv_parse(blob);
	star_finder_params out = { 0 };
	cr_assert(synthstar_conf_from_kv(parsed, &out), "conf must parse");
	assert_sf_conf_eq(&out, &in);
	g_hash_table_unref(parsed);
	/* a blob missing a required key must fail closed */
	GHashTable *partial = nde_kv_parse("sf_radius=5;sf_sigma=1");
	star_finder_params dummy = { 0 };
	cr_assert(!synthstar_conf_from_kv(partial, &dummy),
	          "an incomplete conf must not parse");
	g_hash_table_unref(partial);
	g_free(blob);
}

Test(nde_serializers, synthstar_roundtrip_delegated) {
	const op_descriptor *op = op_descriptor_by_id("star.synthstar");
	cr_assert_not_null(op);
	struct synthstar_data *in = new_synthstar_data();
	in->star_auto = TRUE;
	in->star_conf = make_test_sf_conf();

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	/* the delegated form records the conf, NOT a star list */
	cr_assert(strstr(blob, "stars_auto=1") != NULL, "blob: %s", blob);
	cr_assert(strstr(blob, "sf_sigma=") != NULL, "blob: %s", blob);
	cr_assert_null(strstr(blob, "stars="), "delegated must not pin a list: %s", blob);

	struct synthstar_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->star_auto, TRUE);
	cr_assert_null(out->stars_blob);
	assert_sf_conf_eq(&out->star_conf, &in->star_conf);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, unclip_roundtrip_delegated) {
	const op_descriptor *op = op_descriptor_by_id("star.unclip");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct synthstar_data *in = new_synthstar_data();
	in->star_auto = TRUE;
	in->star_conf = make_test_sf_conf();

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct synthstar_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->star_auto, TRUE);
	assert_sf_conf_eq(&out->star_conf, &in->star_conf);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* unclip is always delegated: a blob with no conf is non-replayable */
	cr_assert_null(op->deserialize("stars_auto=1", op->version),
	               "unclip without a conf must yield NULL");
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, unpurple_roundtrip_delegated) {
	const op_descriptor *op = op_descriptor_by_id("filters.unpurple");
	cr_assert_not_null(op);
	struct unpurpleargs *in = new_unpurple_args();
	in->mod_b = 0.5;
	in->thresh = 0.25;
	in->withstarmask = TRUE;
	in->star_auto = TRUE;
	in->star_conf = make_test_sf_conf();

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	cr_assert_null(strstr(blob, "stars="), "delegated must not pin a list: %s", blob);
	struct unpurpleargs *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->withstarmask, TRUE);
	cr_assert_eq(out->star_auto, TRUE);
	cr_assert_null(out->stars_blob);
	assert_sf_conf_eq(&out->star_conf, &in->star_conf);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

/* ---- PCC/SPCC: v2 recompute-from-inputs, v1 effective white balance ---- */

Test(nde_serializers, photometric_cc_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("color.photometric_cc");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct photometric_cc_data in = { 0 };
	in.spcc = TRUE;
	in.t0 = -2.5f;
	in.t1 = 2.0f;
	in.have_effective = TRUE;
	in.eff_kw[0] = 1.0f;    in.eff_kw[1] = 0.75f;   in.eff_kw[2] = 0.5f;
	in.eff_bg[0] = 0.125f;  in.eff_bg[1] = 0.0625f; in.eff_bg[2] = 0.25f;
	in.bg_auto = FALSE;
	in.bg_area = (rectangle){ 5, 6, 70, 80 };
	in.mag_mode = LIMIT_MAG_ABSOLUTE;
	in.magnitude_arg = 15.5;
	in.atmos_corr = TRUE;
	in.atmos_obs_height = 1250.0;
	in.atmos_pressure = 1013.25;
	in.spcc_mono_sensor = TRUE;
	in.selected_sensor_m = 7;
	in.selected_filter_r = 3;
	in.selected_white_ref = 2;
	in.nb_mode = TRUE;
	in.nb_center[0] = 656.28; in.nb_center[1] = 500.7; in.nb_center[2] = 672.4;
	in.nb_bandwidth[0] = 3.0; in.nb_bandwidth[1] = 3.5; in.nb_bandwidth[2] = 4.0;

	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct photometric_cc_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	/* v2 contract: the ANALYSIS re-runs at replay; the coefficients survive
	 * as the offline fallback. */
	cr_assert(!out->have_effective, "a v2 record must recompute at replay");
	cr_assert(out->have_fallback, "the recorded coefficients must survive as fallback");
	cr_assert_eq(out->spcc, TRUE);
	for (int c = 0; c < 3; c++) {
		cr_assert(memcmp(&out->eff_kw[c], &in.eff_kw[c], sizeof(float)) == 0,
		          "kw%d not bit-exact", c);
		cr_assert(memcmp(&out->eff_bg[c], &in.eff_bg[c], sizeof(float)) == 0,
		          "bg%d not bit-exact", c);
	}
	/* The analysis inputs round-trip. */
	cr_assert_eq(out->bg_auto, FALSE);
	cr_assert_eq(out->bg_area.x, 5);
	cr_assert_eq(out->bg_area.h, 80);
	cr_assert_eq(out->mag_mode, LIMIT_MAG_ABSOLUTE);
	cr_assert(memcmp(&out->magnitude_arg, &in.magnitude_arg, sizeof(double)) == 0);
	cr_assert_eq(out->atmos_corr, TRUE);
	cr_assert(memcmp(&out->atmos_obs_height, &in.atmos_obs_height, sizeof(double)) == 0);
	cr_assert_eq(out->spcc_mono_sensor, TRUE);
	cr_assert_eq(out->selected_sensor_m, 7);
	cr_assert_eq(out->selected_filter_r, 3);
	cr_assert_eq(out->selected_white_ref, 2);
	cr_assert_eq(out->nb_mode, TRUE);
	cr_assert(memcmp(out->nb_center, in.nb_center, sizeof(in.nb_center)) == 0);
	cr_assert(memcmp(out->nb_bandwidth, in.nb_bandwidth, sizeof(in.nb_bandwidth)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);

	/* LEGACY v1 blob (no auto_replay): applies its coefficients verbatim. */
	struct photometric_cc_data *lout = op->deserialize(
		"catalog=1;spcc=0;t0=-2;t1=2;kw0=1.5;kw1=1;kw2=0.5;"
		"bg0=0.1;bg1=0.1;bg2=0.1", 1);
	cr_assert_not_null(lout);
	cr_assert(lout->have_effective, "a v1 record must replay its stored coefficients");
	FREE_VIA_DESTRUCTOR(lout);

	/* a record without inputs OR coefficients (pre-batch) stays honestly
	 * non-replayable */
	cr_assert_null(op->deserialize("catalog=1;spcc=0;t0=-2;t1=2", op->version),
	               "PCC record without kw keys must yield NULL");
	g_free(blob);
}

/* Golden apply-the-branch test: construct params with have_effective and known
 * kw/bg on a small RGB fixture, run the hook's fast path, and compare against a
 * hand-computed expected using the exact formula in
 * apply_photometric_color_correction: bg_mean = (bg0+bg1+bg2)/3;
 * offset[c] = -bg[c]*kw[c] + bg_mean; out = in*kw[c] + offset[c]. */
Test(nde_serializers, photometric_cc_apply_branch_golden) {
	/* 2x2 RGB float fixture */
	fits f = { 0 };
	f.rx = f.naxes[0] = 2;
	f.ry = f.naxes[1] = 2;
	f.naxes[2] = 3;
	f.type = DATA_FLOAT;
	size_t n = 4;
	f.fdata = calloc(n * 3, sizeof(float));
	f.fpdata[0] = f.fdata;
	f.fpdata[1] = f.fdata + n;
	f.fpdata[2] = f.fdata + 2 * n;
	float inpix[3] = { 0.4f, 0.6f, 0.2f };
	for (int c = 0; c < 3; c++)
		for (size_t i = 0; i < n; i++)
			f.fpdata[c][i] = inpix[c];

	float kw[3] = { 1.0f, 0.75f, 0.5f };
	float bg[3] = { 0.1f, 0.2f, 0.05f };
	cr_assert_eq(apply_photometric_color_correction(&f, kw, bg), 0);

	float bg_mean = (bg[0] + bg[1] + bg[2]) / 3.f;
	for (int c = 0; c < 3; c++) {
		float offset = -bg[c] * kw[c] + bg_mean;
		float expect = inpix[c] * kw[c] + offset;
		for (size_t i = 0; i < n; i++)
			cr_assert_float_eq(f.fpdata[c][i], expect, 1e-6,
			                   "channel %d: got %g expect %g", c,
			                   f.fpdata[c][i], expect);
	}
	free(f.fdata);
}

/* ---- Convention 4: ICC profile source ---- */

Test(nde_serializers, icc_convert_builtin_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("icc.convert");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct icc_data in = { 0 };
	in.destroy_fn = free_icc_data;
	in.intent = 1;   /* INTENT_RELATIVE_COLORIMETRIC */
	in.profile_source = icc_source_for_builtin("srgb");
	in.profile = icc_profile_from_source(in.profile_source, NULL);
	cr_assert_not_null(in.profile);

	gchar *blob = op->serialize(&in);
	cr_assert_not_null(blob);
	struct icc_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_str_eq(out->profile_source, in.profile_source);
	cr_assert_eq(out->intent, in.intent);
	cr_assert_not_null(out->profile);
	/* same builtin → same colour space */
	cr_assert_eq(cmsGetColorSpace(out->profile), cmsGetColorSpace(in.profile));
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* pre-batch record (no profile_source) → non-replayable */
	cr_assert_null(op->deserialize("intent=1", op->version),
	               "icc record without profile_source must yield NULL");
	g_free(blob);
	cmsCloseProfile(in.profile);
	g_free(in.profile_source);
}

Test(nde_serializers, icc_convert_blob_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("icc.convert");
	cr_assert_not_null(op);
	cmsHPROFILE src = srgb_trc();
	cr_assert_not_null(src);
	struct icc_data in = { 0 };
	in.destroy_fn = free_icc_data;
	in.intent = 0;
	in.profile_source = icc_source_blob_from_profile(src);   /* embedded blob */
	cr_assert_not_null(in.profile_source);
	cr_assert(g_str_has_prefix(in.profile_source, "blob:"));
	in.profile = icc_profile_from_source(in.profile_source, NULL);
	cr_assert_not_null(in.profile);

	gchar *blob = op->serialize(&in);
	struct icc_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_not_null(out->profile);
	/* the reopened profile describes the same colour space as the original */
	cr_assert_eq(cmsGetColorSpace(out->profile), cmsGetColorSpace(src));
	FREE_VIA_DESTRUCTOR(out);
	g_free(blob);
	cmsCloseProfile(src);
	cmsCloseProfile(in.profile);
	g_free(in.profile_source);
}

Test(nde_serializers, icc_convert_file_hash_mismatch_refused) {
	const op_descriptor *op = op_descriptor_by_id("icc.convert");
	cr_assert_not_null(op);
	/* write a real profile file to a temp path */
	GError *e = NULL;
	gchar *dir = g_dir_make_tmp("siril_icc_XXXXXX", &e);
	cr_assert_not_null(dir, "temp dir: %s", e ? e->message : "?");
	gchar *path = g_build_filename(dir, "prof.icc", NULL);
	cmsHPROFILE src = srgb_trc();
	guint32 len = 0;
	unsigned char *bytes = get_icc_profile_data(src, &len);
	cr_assert(len > 0);
	cr_assert(g_file_set_contents(path, (char *)bytes, len, NULL));

	struct icc_data in = { 0 };
	in.destroy_fn = free_icc_data;
	in.intent = 0;
	in.profile_source = g_strconcat("file:", path, NULL);
	in.profile = icc_profile_from_source(in.profile_source, NULL);
	cr_assert_not_null(in.profile);

	gchar *blob = op->serialize(&in);
	cr_assert(strstr(blob, "profile_sha256=") != NULL,
	          "file source must pin the hash");
	/* clean deserialize succeeds */
	struct icc_data *ok = op->deserialize(blob, op->version);
	cr_assert_not_null(ok, "unchanged file must deserialize");
	FREE_VIA_DESTRUCTOR(ok);

	/* tamper one byte → hash mismatch → refusal (NULL) */
	bytes[0] ^= 0xFF;
	cr_assert(g_file_set_contents(path, (char *)bytes, len, NULL));
	cr_assert_null(op->deserialize(blob, op->version),
	               "a changed profile file must be refused at replay");

	g_free(blob);
	free(bytes);
	cmsCloseProfile(src);
	cmsCloseProfile(in.profile);
	g_free(in.profile_source);
	g_remove(path); g_free(path);
	g_rmdir(dir); g_free(dir);
}

/* ------------------------------------------------------------------ *
 *  Batch — mask ops (core/masks.c).  Serializers make these Tier A so *
 *  they no longer act as barriers in an NDE chain.                    *
 * ------------------------------------------------------------------ */

Test(nde_serializers, mask_from_stars_roundtrip) {
	mask_from_stars_data in = { 0 };
	in.r = 1.5f + 0.3f;
	in.feather = 0.3f;
	in.invert = TRUE;
	in.bitdepth = 32;

	gchar *blob = op_desc_mask_from_stars.serialize(&in);
	cr_assert_not_null(blob);
	mask_from_stars_data *out = op_desc_mask_from_stars.deserialize(blob, op_desc_mask_from_stars.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->r, &in.r, sizeof(float)) == 0);
	cr_assert(memcmp(&out->feather, &in.feather, sizeof(float)) == 0);
	cr_assert_eq(out->invert, in.invert);
	cr_assert_eq(out->bitdepth, in.bitdepth);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_from_stars, blob);
	g_free(blob);
}

Test(nde_serializers, mask_from_channel_roundtrip_with_file) {
	gchar *path = write_temp_operand("mask from channel operand bytes\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	mask_from_channel_data in = { 0 };
	in.channel = 2;
	in.autostretch = TRUE;
	in.invert = TRUE;
	in.bitpix = 32;
	in.filename = path;
	in.fit = (fits *)0xdeadbeef;   /* transient: must NOT be serialized */

	gchar *blob = op_desc_mask_from_channel.serialize(&in);
	cr_assert_not_null(blob);
	mask_from_channel_data *out = op_desc_mask_from_channel.deserialize(blob, op_desc_mask_from_channel.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->channel, in.channel);
	cr_assert_eq(out->autostretch, in.autostretch);
	cr_assert_eq(out->invert, in.invert);
	cr_assert_eq(out->bitpix, in.bitpix);
	cr_assert_str_eq(out->filename, path);
	cr_assert_null(out->fit, "fit is transient and must not round-trip");
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	gint64 stored_size = 0;
	cr_assert(nde_kv_get_int(kv, "operand_size", &stored_size));
	cr_assert_eq(stored_size, size);
	g_hash_table_unref(kv);
	/* the struct owns filename now — the destructor frees it */
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_from_channel, blob);
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

Test(nde_serializers, mask_from_channel_roundtrip_no_file) {
	mask_from_channel_data in = { 0 };
	in.channel = 1;
	in.autostretch = FALSE;
	in.invert = TRUE;
	in.bitpix = 16;
	in.filename = NULL;   /* build from the current image */

	gchar *blob = op_desc_mask_from_channel.serialize(&in);
	cr_assert_not_null(blob);
	cr_assert_null(strstr(blob, "operand_path"),
	               "no filename must not emit an operand_path key");
	mask_from_channel_data *out = op_desc_mask_from_channel.deserialize(blob, op_desc_mask_from_channel.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->channel, in.channel);
	cr_assert_eq(out->autostretch, in.autostretch);
	cr_assert_eq(out->invert, in.invert);
	cr_assert_eq(out->bitpix, in.bitpix);
	cr_assert_null(out->filename, "absent operand_path must round-trip as NULL");
	FREE_VIA_DESTRUCTOR(out);
	/* an explicitly empty operand_path must also mean "no file", not "" */
	mask_from_channel_data *empty = op_desc_mask_from_channel.deserialize(
	    "channel=1;autostretch=0;invert=1;bitpix=16;operand_path=",
	    op_desc_mask_from_channel.version);
	cr_assert_not_null(empty);
	cr_assert_null(empty->filename, "empty operand_path must round-trip as NULL");
	FREE_VIA_DESTRUCTOR(empty);
	CHECK_MALFORMED(&op_desc_mask_from_channel, blob);
	g_free(blob);
}

Test(nde_serializers, mask_from_luminance_roundtrip_with_file) {
	gchar *path = write_temp_operand("mask from luminance operand bytes\n");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(path, &size);
	cr_assert_not_null(sha);

	mask_from_lum_data in = { 0 };
	in.rw = 0.3f; in.gw = 0.6f; in.bw = 0.1f;
	in.autostretch = TRUE;
	in.invert = TRUE;
	in.use_human = TRUE;
	in.use_even = FALSE;
	in.bitpix = 32;
	in.filename = path;
	in.fit = (fits *)0xdeadbeef;   /* transient: must NOT be serialized */

	gchar *blob = op_desc_mask_from_luminance.serialize(&in);
	cr_assert_not_null(blob);
	mask_from_lum_data *out = op_desc_mask_from_luminance.deserialize(blob, op_desc_mask_from_luminance.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->rw, &in.rw, sizeof(float)) == 0);
	cr_assert(memcmp(&out->gw, &in.gw, sizeof(float)) == 0);
	cr_assert(memcmp(&out->bw, &in.bw, sizeof(float)) == 0);
	cr_assert_eq(out->autostretch, in.autostretch);
	cr_assert_eq(out->invert, in.invert);
	cr_assert_eq(out->use_human, in.use_human);
	cr_assert_eq(out->use_even, in.use_even);
	cr_assert_eq(out->bitpix, in.bitpix);
	cr_assert_str_eq(out->filename, path);
	cr_assert_null(out->fit, "fit is transient and must not round-trip");
	GHashTable *kv = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(kv, "operand_sha256"), sha);
	gint64 stored_size = 0;
	cr_assert(nde_kv_get_int(kv, "operand_size", &stored_size));
	cr_assert_eq(stored_size, size);
	g_hash_table_unref(kv);
	/* the struct owns filename now — the destructor frees it */
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_from_luminance, blob);
	g_free(blob); g_free(sha);
	g_remove(path); g_free(path);
}

Test(nde_serializers, mask_from_luminance_roundtrip_no_file) {
	mask_from_lum_data in = { 0 };
	in.rw = -0.25f; in.gw = 0.7f; in.bw = 0.3f;   /* rw < 0 is a legal sentinel */
	in.autostretch = FALSE;
	in.invert = FALSE;
	in.use_human = FALSE;
	in.use_even = TRUE;
	in.bitpix = 8;
	in.filename = NULL;

	gchar *blob = op_desc_mask_from_luminance.serialize(&in);
	cr_assert_not_null(blob);
	cr_assert_null(strstr(blob, "operand_path"),
	               "no filename must not emit an operand_path key");
	mask_from_lum_data *out = op_desc_mask_from_luminance.deserialize(blob, op_desc_mask_from_luminance.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->rw, &in.rw, sizeof(float)) == 0);
	cr_assert(memcmp(&out->gw, &in.gw, sizeof(float)) == 0);
	cr_assert(memcmp(&out->bw, &in.bw, sizeof(float)) == 0);
	cr_assert_eq(out->use_even, in.use_even);
	cr_assert_eq(out->bitpix, in.bitpix);
	cr_assert_null(out->filename, "absent operand_path must round-trip as NULL");
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_from_luminance, blob);
	g_free(blob);
}

Test(nde_serializers, mask_threshold_roundtrip) {
	mask_thresh_data in = { 0 };
	in.min_val = 0.3f;
	in.max_val = 0.7f + 0.03125f;
	in.range = 1.5f;

	gchar *blob = op_desc_mask_threshold.serialize(&in);
	cr_assert_not_null(blob);
	mask_thresh_data *out = op_desc_mask_threshold.deserialize(blob, op_desc_mask_threshold.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->min_val, &in.min_val, sizeof(float)) == 0);
	cr_assert(memcmp(&out->max_val, &in.max_val, sizeof(float)) == 0);
	cr_assert(memcmp(&out->range, &in.range, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_threshold, blob);
	g_free(blob);
}

Test(nde_serializers, mask_blur_roundtrip) {
	mask_blur_data in = { 0 };
	in.radius = 3.3f;

	gchar *blob = op_desc_mask_blur.serialize(&in);
	cr_assert_not_null(blob);
	mask_blur_data *out = op_desc_mask_blur.deserialize(blob, op_desc_mask_blur.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->radius, &in.radius, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_blur, blob);
	g_free(blob);
}

Test(nde_serializers, mask_feather_roundtrip) {
	mask_feather_data in = { 0 };
	in.distance = 7.3f;
	in.mode = FEATHER_EDGE;   /* 2 */

	gchar *blob = op_desc_mask_feather.serialize(&in);
	cr_assert_not_null(blob);
	mask_feather_data *out = op_desc_mask_feather.deserialize(blob, op_desc_mask_feather.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->distance, &in.distance, sizeof(float)) == 0);
	cr_assert_eq(out->mode, in.mode);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_feather, blob);
	g_free(blob);
}

Test(nde_serializers, mask_multiply_roundtrip) {
	mask_fmul_data in = { 0 };
	in.factor = 0.3f;

	gchar *blob = op_desc_mask_multiply.serialize(&in);
	cr_assert_not_null(blob);
	mask_fmul_data *out = op_desc_mask_multiply.deserialize(blob, op_desc_mask_multiply.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->factor, &in.factor, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_multiply, blob);
	g_free(blob);
}

Test(nde_serializers, mask_bitpix_roundtrip) {
	mask_bitpix_data in = { 0 };
	in.bitpix = 32;

	gchar *blob = op_desc_mask_bitpix.serialize(&in);
	cr_assert_not_null(blob);
	mask_bitpix_data *out = op_desc_mask_bitpix.deserialize(blob, op_desc_mask_bitpix.version);
	cr_assert_not_null(out);
	cr_assert_eq(out->bitpix, in.bitpix);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_bitpix, blob);
	g_free(blob);
}

Test(nde_serializers, mask_from_color_roundtrip) {
	mask_from_color_data in = { 0 };
	in.chrom_center_r = 0.3f;
	in.chrom_center_g = 0.6f;
	in.chrom_center_b = 0.1f;
	in.chrom_tolerance = 0.125f + 0.3f;
	in.lum_min = -0.25f;   /* negatives are representable here */
	in.lum_max = 0.9f;
	in.feather_radius = 5;
	in.invert = TRUE;
	in.bitpix = 16;
	in.cleanup = TRUE;

	gchar *blob = op_desc_mask_from_color.serialize(&in);
	cr_assert_not_null(blob);
	mask_from_color_data *out = op_desc_mask_from_color.deserialize(blob, op_desc_mask_from_color.version);
	cr_assert_not_null(out);
	cr_assert(memcmp(&out->chrom_center_r, &in.chrom_center_r, sizeof(float)) == 0);
	cr_assert(memcmp(&out->chrom_center_g, &in.chrom_center_g, sizeof(float)) == 0);
	cr_assert(memcmp(&out->chrom_center_b, &in.chrom_center_b, sizeof(float)) == 0);
	cr_assert(memcmp(&out->chrom_tolerance, &in.chrom_tolerance, sizeof(float)) == 0);
	cr_assert(memcmp(&out->lum_min, &in.lum_min, sizeof(float)) == 0);
	cr_assert(memcmp(&out->lum_max, &in.lum_max, sizeof(float)) == 0);
	cr_assert_eq(out->feather_radius, in.feather_radius);
	cr_assert_eq(out->invert, in.invert);
	cr_assert_eq(out->bitpix, in.bitpix);
	cr_assert_eq(out->cleanup, in.cleanup);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(&op_desc_mask_from_color, blob);
	g_free(blob);
}

/* The four param-less mask ops: an empty blob is VALID (their hooks ignore
 * args->user), so CHECK_MALFORMED must NOT be used — only a newer version
 * fails.  Same contract as check_paramless() above. */
Test(nde_serializers, mask_clear_paramless)        { check_paramless("mask.clear"); }
Test(nde_serializers, mask_invert_paramless)       { check_paramless("mask.invert"); }
Test(nde_serializers, mask_autostretch_paramless)  { check_paramless("mask.autostretch"); }
Test(nde_serializers, mask_from_gradient_paramless){ check_paramless("mask.from_gradient"); }

/* ------------------------------------------------------------------ *
 *  The direct-apply colour tools, and the PSF stash                   *
 * ------------------------------------------------------------------ */

Test(nde_serializers, bkg_neutral_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("color.background_neutralization");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize, "the tool must no longer be opaque");
	struct bkg_neutral_data *in = new_bkg_neutral_data();
	in->black_selection.x = 41;
	in->black_selection.y = 97;
	in->black_selection.w = 256;
	in->black_selection.h = 128;
	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct bkg_neutral_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->black_selection.x, in->black_selection.x);
	cr_assert_eq(out->black_selection.y, in->black_selection.y);
	cr_assert_eq(out->black_selection.w, in->black_selection.w);
	cr_assert_eq(out->black_selection.h, in->black_selection.h);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* An empty rectangle samples nothing, so it cannot be a recipe. */
	cr_assert_null(op->deserialize("bkg_x=0;bkg_y=0;bkg_w=0;bkg_h=0", op->version));
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

/* Since v2 the automatic mode's contract is INTENT (nde_joint.h lineage):
 * the record carries the selections and re-samples them over the replayed
 * upstream pixels, so an upstream edit re-derives the calibration instead
 * of replaying frozen coefficients.  The coefficients stay recorded as
 * diagnostics — and as the exact contract for manual and legacy records. */
Test(nde_serializers, color_calibration_records_the_computed_coefficients) {
	const op_descriptor *op = op_descriptor_by_id("color.calibration");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct color_calib_data *in = new_color_calib_data();
	in->is_manual = FALSE;
	in->white_selection = (rectangle){ 10, 20, 30, 40 };
	in->black_selection = (rectangle){ 50, 60, 70, 80 };
	in->low = 0.25;
	in->high = 0.875;
	in->kw[0] = 1.125; in->kw[1] = 1.0; in->kw[2] = 0.9375;
	in->bg[0] = 0.03125; in->bg[1] = 0.0625; in->bg[2] = 0.125;
	in->have_effective = TRUE;
	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct color_calib_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	for (int i = 0; i < 3; i++) {
		cr_assert(memcmp(&out->kw[i], &in->kw[i], sizeof(double)) == 0);
		cr_assert(memcmp(&out->bg[i], &in->bg[i], sizeof(double)) == 0);
	}
	/* The whole point: an automatic v2 record re-runs the analysis at
	 * replay (the hook's !have_effective branch). */
	cr_assert(!out->have_effective, "an automatic calibration must recompute at replay");
	cr_assert_eq(out->is_manual, FALSE);
	cr_assert_eq(out->white_selection.w, 30);
	cr_assert_eq(out->black_selection.h, 80);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);

	/* MANUAL records keep the stored coefficients: sliders ARE the params. */
	in->is_manual = TRUE;
	gchar *mblob = op->serialize(in);
	struct color_calib_data *mout = op->deserialize(mblob, op->version);
	cr_assert_not_null(mout);
	cr_assert(mout->have_effective, "a manual calibration must not re-sample");
	FREE_VIA_DESTRUCTOR(mout);
	g_free(mblob);

	/* LEGACY v1 blobs (no auto_replay key) keep the exact contract too. */
	struct color_calib_data *lout = op->deserialize(
		"kw0=1.5;bg0=0.1;kw1=1;bg1=0.1;kw2=0.5;bg2=0.1;manual=0;"
		"white_x=1;white_y=2;white_w=3;white_h=4", 1);
	cr_assert_not_null(lout);
	cr_assert(lout->have_effective, "a legacy record must replay its stored coefficients");
	FREE_VIA_DESTRUCTOR(lout);

	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, a_calibration_without_coefficients_is_not_replayable) {
	const op_descriptor *op = op_descriptor_by_id("color.calibration");
	cr_assert_not_null(op);
	/* Rectangles and limits alone: what a pre-feature record would look like.
	 * It must stay honestly non-replayable rather than re-sample. */
	cr_assert_null(op->deserialize("manual=0;white_x=1;white_y=2;white_w=3;"
	                               "white_h=4;low=0.1;high=0.9", op->version));
}

Test(nde_serializers, rgb_align_roundtrip) {
	const op_descriptor *op = op_descriptor_by_id("color.rgb_align");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	struct rgb_align_data *in = new_rgb_align_data();
	in->method = 2;
	in->area = (rectangle){ 128, 256, 512, 512 };
	in->have_area = TRUE;
	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	struct rgb_align_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->method, 2);
	cr_assert(out->have_area, "the registration area is the one thing replay "
	                          "cannot re-derive from the image");
	cr_assert_eq(out->area.x, 128);
	cr_assert_eq(out->area.w, 512);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	/* There are four methods; anything else is not one of ours. */
	cr_assert_null(op->deserialize("method=9", op->version));
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

/* Convention 2 applied to a PSF: the kernel travels with the record, because
 * it never lived in the params struct — it lives in com.kernel, put there by a
 * makepsf run this record knows nothing about. */
Test(nde_serializers, deconvolve_carries_the_kernel_it_used) {
	const op_descriptor *op = op_descriptor_by_id("filters.deconvolve");
	cr_assert_not_null(op);
	cr_assert_not_null(op->serialize);
	estk_data *in = calloc(1, sizeof(estk_data));
	cr_assert_not_null(in);
	in->destroy_fn = free_estk_data;
	const int ks = 5, kc = 1;
	in->eff_ks = ks;
	in->eff_kchans = kc;
	in->eff_kernel = malloc(ks * ks * kc * sizeof(float));
	cr_assert_not_null(in->eff_kernel);
	for (int i = 0; i < ks * ks * kc; i++)
		in->eff_kernel[i] = (float)i / 32.0f;   /* exact in binary */
	in->psftype = PSF_BLIND;
	in->alpha = 1.0f / 3000.0f;
	in->finaliters = 7;
	in->stepsize = 0.125f;

	gchar *blob = op->serialize(in);
	cr_assert_not_null(blob);
	estk_data *out = op->deserialize(blob, op->version);
	cr_assert_not_null(out);
	cr_assert_eq(out->eff_ks, ks);
	cr_assert_eq(out->eff_kchans, kc);
	cr_assert_not_null(out->eff_kernel);
	cr_assert(memcmp(out->eff_kernel, in->eff_kernel,
	                 ks * ks * kc * sizeof(float)) == 0,
	          "the kernel must survive base64 bit-for-bit");
	/* Replay uses the recorded PSF; it must never re-run the estimator. */
	cr_assert_eq(out->psftype, PSF_PREVIOUS);
	cr_assert_eq(out->finaliters, 7);
	cr_assert(memcmp(&out->stepsize, &in->stepsize, sizeof(float)) == 0);
	FREE_VIA_DESTRUCTOR(out);
	CHECK_MALFORMED(op, blob);
	FREE_VIA_DESTRUCTOR(in);
	g_free(blob);
}

Test(nde_serializers, a_deconvolution_without_a_kernel_is_not_replayable) {
	const op_descriptor *op = op_descriptor_by_id("filters.deconvolve");
	cr_assert_not_null(op);
	/* psftype alone is a way of looking for a PSF, not a PSF. */
	cr_assert_null(op->deserialize("psftype=0;alpha=0.000333;finaliters=1",
	                               op->version));
	/* A blob whose length disagrees with ks is corrupt, not "close enough". */
	cr_assert_null(op->deserialize("kernel_ks=5;kernel_kchans=1;kernel=AAAA",
	                               op->version));
}

/* The Convention-3 contract is a claim about the HOOK, not just the codec:
 * a replayed calibration must apply the recorded coefficients and must not
 * go back to the rectangles.  These two tests are the claim. */

static fits *rgb_fit(int rx, int ry, float v) {
	fits *f = NULL;
	if (new_fit_image(&f, rx, ry, 3, DATA_FLOAT))
		return NULL;
	for (int c = 0; c < 3; c++)
		for (size_t i = 0; i < (size_t)rx * ry; i++)
			f->fpdata[c][i] = v * (float)(c + 1);
	return f;
}

Test(nde_serializers, a_replayed_calibration_applies_what_was_recorded) {
	fits *f = rgb_fit(16, 16, 0.1f);
	cr_assert_not_null(f);
	struct color_calib_data *p = new_color_calib_data();
	p->have_effective = TRUE;              /* as deserialize leaves it */
	p->kw[0] = 2.0; p->kw[1] = 1.0; p->kw[2] = 0.5;
	p->bg[0] = 0.0; p->bg[1] = 0.0; p->bg[2] = 0.0;
	/* Rectangles that would sample somewhere quite different if consulted. */
	p->white_selection = (rectangle){ 0, 0, 4, 4 };
	p->black_selection = (rectangle){ 8, 8, 4, 4 };
	p->is_manual = FALSE;
	const double kw_before[3] = { p->kw[0], p->kw[1], p->kw[2] };

	struct generic_img_args args = { 0 };
	args.user = p;
	cr_assert_eq(color_calib_image_hook(&args, f, 1), 0);

	/* Untouched: sampling would have overwritten these. */
	for (int i = 0; i < 3; i++)
		cr_assert(memcmp(&p->kw[i], &kw_before[i], sizeof(double)) == 0,
		          "channel %d: replay re-sampled instead of using the record", i);
	/* kw == 1.0 is the skip case, so the green channel must be untouched. */
	cr_assert_float_eq(f->fpdata[1][0], 0.2f, 1e-6f);
	cr_assert_float_eq(f->fpdata[0][0], 0.2f, 1e-6f);   /* 0.1 * 2.0 */
	cr_assert_float_eq(f->fpdata[2][0], 0.15f, 1e-6f);  /* 0.3 * 0.5 */
	FREE_VIA_DESTRUCTOR(p);
	clearfits(f);
	free(f);
}

Test(nde_serializers, a_first_run_works_the_coefficients_out_and_keeps_them) {
	fits *f = rgb_fit(32, 32, 0.25f);
	cr_assert_not_null(f);
	struct color_calib_data *p = new_color_calib_data();
	p->is_manual = FALSE;
	p->have_effective = FALSE;             /* as the dialog leaves it */
	p->white_selection = (rectangle){ 0, 0, 16, 16 };
	p->black_selection = (rectangle){ 16, 16, 16, 16 };
	p->low = 0.0;
	p->high = 1.0;

	struct generic_img_args args = { 0 };
	args.user = p;
	cr_assert_eq(color_calib_image_hook(&args, f, 1), 0);
	/* The capture that follows this run is what the record is made of, so the
	 * hook has to leave the worked-out values behind for the serializer. */
	cr_assert(p->have_effective,
	          "the capturing run must leave its coefficients for the record");
	FREE_VIA_DESTRUCTOR(p);
	clearfits(f);
	free(f);
}

Test(nde_serializers, the_colour_tools_refuse_a_mono_image) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, 16, 16, 1, DATA_FLOAT), 0);
	struct generic_img_args args = { 0 };

	/* background_neutralize used to assert() here, which was survivable only
	 * while a dialog button was its one caller.  A replay can hand it anything
	 * the history claims was there, so it has to answer instead. */
	struct bkg_neutral_data *b = new_bkg_neutral_data();
	b->black_selection = (rectangle){ 0, 0, 8, 8 };
	args.user = b;
	cr_assert_neq(bkg_neutral_image_hook(&args, f, 1), 0);
	FREE_VIA_DESTRUCTOR(b);

	struct color_calib_data *c = new_color_calib_data();
	c->have_effective = TRUE;
	args.user = c;
	cr_assert_neq(color_calib_image_hook(&args, f, 1), 0);
	FREE_VIA_DESTRUCTOR(c);

	clearfits(f);
	free(f);
}
