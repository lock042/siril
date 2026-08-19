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

/* Does an op's declared halo actually make its region result right?
 *
 * A region preview crops W (grown by the op's roi_halo), runs the hook on the
 * crop, and writes back W.  The claim a halo makes is:
 *
 *     run_on_crop(W + h) restricted to W  ==  run_on_full(image) restricted to W
 *
 * This file measures the deviation between those two, directly against the
 * image hooks — no dialogs, no GUI, so it runs in CI.  For an op declaring
 * roi_halo_exact it ASSERTS zero, which is the point of the file: a halo that
 * is one pixel short produces a result that is plausible, stable and quietly
 * wrong at the edges of every ROI, and nothing else in the suite would notice.
 * For a deliberately approximate halo (deconvolution's iteration count, a
 * non-local search radius) it only reports, because there is no true bound to
 * assert against.
 *
 * The sweep tests report deviation as a function of h.  They are what turns
 * "this halo looks right" into evidence: the curve should fall to zero exactly
 * at the op's declared halo and stay there. */

#include <criterion/criterion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/op_descriptor.h"
#include "core/op_descriptors.h"
#include "core/fits_region.h"
#include "io/image_format_fits.h"
#include "filters/median.h"
#include "filters/epf.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

/* The hooks read com.max_thread directly in their OpenMP clauses; leaving it
 * at 0 is not merely slow, num_threads(0) is invalid.  One thread also makes
 * every run here bit-reproducible. */
static void setup(void) {
	com.max_thread = 1;
}

/* The image under test.  Deterministic but structured: a smooth gradient plus
 * a reproducible pseudo-random speckle, so a median has something to bite on
 * (a pure gradient is its own median and would hide a wrong halo) and the
 * result still does not depend on the platform's RNG. */
static fits *make_image(uint32_t rx, uint32_t ry, int nchans, data_type type) {
	fits *f = calloc(1, sizeof(fits));
	cr_assert_not_null(f);
	if (new_fit_image(&f, rx, ry, nchans, type))
		cr_assert(0, "new_fit_image failed");

	const size_t npix = (size_t)rx * ry;
	uint32_t state = 12345u;
	for (int c = 0; c < nchans; c++) {
		for (uint32_t y = 0; y < ry; y++) {
			for (uint32_t x = 0; x < rx; x++) {
				state = state * 1103515245u + 12345u;
				const uint32_t noise = (state >> 16) & 0x7fffu;
				const double smooth = 20000.0
					+ 12000.0 * sin(x * 0.07 + c) * cos(y * 0.05);
				double v = smooth + (double)(noise % 9000u) - 4500.0;
				if (v < 0.0) v = 0.0;
				if (v > 65535.0) v = 65535.0;
				const size_t i = (size_t)c * npix + (size_t)y * rx + x;
				if (type == DATA_FLOAT)
					f->fdata[i] = (float)(v / 65535.0);
				else
					f->data[i] = (WORD)v;
			}
		}
	}
	return f;
}

static void free_image(fits *f) {
	clearfits(f);
	free(f);
}

/* Largest absolute difference between two same-sized regions, in ADU for
 * ushort and in normalised units scaled to ADU for float, so the two are
 * comparable and a "0" means the same thing in both. */
static double max_abs_diff(const fits *a, const fits *b) {
	cr_assert_eq(a->rx, b->rx);
	cr_assert_eq(a->ry, b->ry);
	cr_assert_eq(a->naxes[2], b->naxes[2]);
	cr_assert_eq(a->type, b->type);
	const size_t n = (size_t)a->rx * a->ry * a->naxes[2];
	double worst = 0.0;
	for (size_t i = 0; i < n; i++) {
		double d;
		if (a->type == DATA_FLOAT)
			d = fabs((double)a->fdata[i] - (double)b->fdata[i]) * 65535.0;
		else
			d = fabs((double)a->data[i] - (double)b->data[i]);
		if (d > worst)
			worst = d;
	}
	return worst;
}

/* Run @op's image hook over a private copy of @src and hand the result back.
 * for_roi/roi_rect are set exactly as generic_image_worker sets them, so a
 * hook that reads them (wavelets windows its output, epf crops its guide)
 * behaves here as it does in the application. */
static fits *run_hook(const op_descriptor *op, gpointer user, fits *src,
                      gboolean for_roi, rectangle crop_rect) {
	fits *work = calloc(1, sizeof(fits));
	cr_assert_not_null(work);
	cr_assert_eq(copyfits(src, work, CP_ALLOC | CP_FORMAT | CP_COPYA, -1), 0);

	struct generic_img_args args = { 0 };
	args.fit = work;
	args.op = op;
	args.user = user;
	args.max_threads = 1;          /* determinism over speed */
	args.for_preview = for_roi;
	args.for_roi = for_roi;
	args.roi_rect = crop_rect;
	op_descriptor_fill_img_args(&args);

	cr_assert_eq(args.image_hook(&args, work, 1), 0, "hook failed");
	return work;
}

/* The measurement this whole file exists for.
 *
 * Returns max |run_on_crop(W (+) halo)|W - run_on_full(image)|W|. */
static double region_deviation(const op_descriptor *op, gpointer user,
                               fits *img, rectangle W, int halo) {
	/* reference: the whole image, then take W */
	fits *full = run_hook(op, user, img, FALSE, (rectangle){ 0, 0, img->rx, img->ry });
	fits ref = { 0 };
	cr_assert_eq(crop_fits_region(full, &W, &ref), 0);

	/* region: crop W (+) halo out of the ORIGINAL, run, take the inner W */
	const gint x0 = MAX(0, W.x - halo);
	const gint y0 = MAX(0, W.y - halo);
	const gint x1 = MIN((gint)img->rx, W.x + W.w + halo);
	const gint y1 = MIN((gint)img->ry, W.y + W.h + halo);
	const rectangle grown = { x0, y0, x1 - x0, y1 - y0 };

	fits crop = { 0 };
	cr_assert_eq(crop_fits_region(img, &grown, &crop), 0);
	fits *region = run_hook(op, user, &crop, TRUE, grown);

	const rectangle inner = { W.x - grown.x, W.y - grown.y, W.w, W.h };
	fits got = { 0 };
	cr_assert_eq(crop_fits_region(region, &inner, &got), 0);

	const double dev = max_abs_diff(&ref, &got);

	clearfits(&ref);
	clearfits(&got);
	clearfits(&crop);
	free_image(full);
	free_image(region);
	return dev;
}

/* Deliberately off-centre and away from every border, so the comparison is
 * about the halo and not about how the hook handles the image edge. */
static const rectangle test_window = { 37, 23, 40, 32 };

/* What "exact" is allowed to mean.  roi_halo_exact claims the halo is a true
 * bound on the op's SPATIAL SUPPORT — that no pixel further away can influence
 * the result.  It does not claim the arithmetic is bit-reproducible, and for
 * OpenCV-backed ops it is not: the same pixels summed in a different order
 * (different image width, different SIMD blocking) land a few units in the
 * last place apart.  Measured, that is ~0.01 ADU out of 65535 for the
 * edge-preserving filters, so anything at or below 1 ADU is noise and the
 * DISCRIMINATION assertion below is what actually catches a wrong halo.
 * Median is integer-exact and asserts a true zero. */
#define NEGLIGIBLE_ADU 1.0

/* ------------------------------------------------------------------ *
 *  Median                                                            *
 * ------------------------------------------------------------------ */

static struct median_filter_data median_params(int ksize, int iterations) {
	struct median_filter_data p = { 0 };
	p.ksize = ksize;
	p.iterations = iterations;
	p.amount = 1.0;
	return p;
}

/* The sweep that produces the evidence: deviation must fall to zero at the
 * declared halo, and must NOT be zero below it — otherwise the test would
 * pass for an op whose halo is irrelevant and prove nothing. */
Test(roi_halo, median_deviation_vanishes_at_its_declared_halo, .init = setup) {
	fits *img = make_image(128, 112, 1, DATA_USHORT);

	const struct { int ksize, iters; } cases[] = {
		{ 3, 1 }, { 5, 1 }, { 3, 3 }, { 9, 2 },
	};

	for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
		struct median_filter_data p = median_params(cases[c].ksize, cases[c].iters);
		const int declared = op_descriptor_roi_halo(&op_desc_median, &p);
		const int expect = ((cases[c].ksize - 1) / 2) * cases[c].iters;

		cr_assert_eq(declared, expect,
			     "median ksize=%d iters=%d: declared halo %d, expected %d",
			     cases[c].ksize, cases[c].iters, declared, expect);

		const double at = region_deviation(&op_desc_median, &p, img, test_window, declared);
		cr_assert_eq(at, 0.0,
			     "median ksize=%d iters=%d: deviation %.1f ADU at its own halo (%d px)",
			     cases[c].ksize, cases[c].iters, at, declared);

		/* one pixel short must be measurably wrong, or the halo is not
		 * what makes this work and the assertion above is vacuous */
		if (declared > 0) {
			const double short_by_one =
				region_deviation(&op_desc_median, &p, img, test_window, declared - 1);
			cr_assert_gt(short_by_one, 0.0,
				     "median ksize=%d iters=%d: halo %d gives the same answer as "
				     "%d, so the declared halo is not doing the work",
				     cases[c].ksize, cases[c].iters, declared, declared - 1);
		}
	}

	free_image(img);
}

/* Same claim on float data and on colour, since the hook takes a different
 * branch for each and a halo bug could hide in either. */
Test(roi_halo, median_is_exact_for_float_and_rgb, .init = setup) {
	struct median_filter_data p = median_params(5, 2);
	const int halo = op_descriptor_roi_halo(&op_desc_median, &p);

	fits *f = make_image(128, 112, 1, DATA_FLOAT);
	cr_assert_eq(region_deviation(&op_desc_median, &p, f, test_window, halo), 0.0,
		     "float median deviates at its declared halo");
	free_image(f);

	fits *rgb = make_image(128, 112, 3, DATA_USHORT);
	cr_assert_eq(region_deviation(&op_desc_median, &p, rgb, test_window, halo), 0.0,
		     "rgb median deviates at its declared halo");
	free_image(rgb);
}

/* ------------------------------------------------------------------ *
 *  Edge-preserving filter                                            *
 * ------------------------------------------------------------------ */

/* sigma_col matters more than it looks: edge_preserving_filter scales it down
 * by 100, then 25 for mono, then 5 for guided, and squares it for the OpenCV
 * call.  With the GUI's typical values the filter is very nearly the identity
 * on synthetic data, and a near-identity filter has no spatial support to
 * measure — every halo scores zero and the test proves nothing.  These values
 * are chosen so the filter actually mixes neighbours. */
static struct epfargs epf_params(ep_filter_t filter, double d, double sigma_space,
                                 double sigma_col) {
	struct epfargs p = { 0 };
	p.filter = filter;
	p.d = d;
	p.sigma_col = sigma_col;
	p.sigma_space = sigma_space;
	p.mod = 1.0;
	p.guidefit = NULL;   /* self-guided: the hook points it at its own fit */
	return p;
}

/* Float, not ushort: a 16-bit round trip quantises the comparison to +/-1 ADU
 * and hides exactly the sub-ADU behaviour this file is trying to measure. */
Test(roi_halo, epf_deviation_vanishes_at_its_declared_halo, .init = setup) {
	fits *img = make_image(128, 112, 1, DATA_FLOAT);

	const struct { const char *what; ep_filter_t f; double d, ss; int expect; } cases[] = {
		/* bilateral: OpenCV radius = d/2 ... */
		{ "bilateral d=9",              EP_BILATERAL, 9.0, 3.0, 4 },
		/* ... or cvRound(sigma_space * 1.5) when d <= 0 */
		{ "bilateral d=0 sigma=5",      EP_BILATERAL, 0.0, 5.0, 8 },
		/* guided: r = d/3, box-filtered twice, so the support is 2r */
		{ "guided d=12",                EP_GUIDED,   12.0, 4.0, 8 },
	};

	for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
		struct epfargs p = epf_params(cases[c].f, cases[c].d, cases[c].ss, 700.0);
		const int declared = op_descriptor_roi_halo(&op_desc_epf, &p);
		cr_assert_eq(declared, cases[c].expect,
			     "%s: declared halo %d, expected %d",
			     cases[c].what, declared, cases[c].expect);

		const int starved = declared / 2;
		const double at       = region_deviation(&op_desc_epf, &p, img, test_window, declared);
		const double starving = region_deviation(&op_desc_epf, &p, img, test_window, starved);

		cr_assert_leq(at, NEGLIGIBLE_ADU,
			      "%s: deviation %.4f ADU at its own halo (%d px)",
			      cases[c].what, at, declared);
		/* The load-bearing half.  Without it the assertion above would pass
		 * just as well for a halo ten times too large, or for parameters
		 * that make the filter a no-op.
		 *
		 * The comparison is against HALF the declared halo rather than one
		 * pixel less, because these filters weight distant neighbours very
		 * lightly: the guided filter is already within 4 ADU at one pixel
		 * short, which would make a one-pixel test a coin toss even though
		 * the halo is provably right (it is 248 ADU out at half).  Half is
		 * still a claim about THIS halo — nothing else in the op changes —
		 * and the measured margins are 2-4 orders of magnitude. */
		cr_assert_gt(starving, 10.0 * NEGLIGIBLE_ADU,
			     "%s: halo %d scores %.4f and halo %d scores %.4f — the "
			     "declared halo is not what makes this work",
			     cases[c].what, declared, at, starved, starving);
	}

	free_image(img);
}

/* ------------------------------------------------------------------ *
 *  Ops that have NOT declared a halo                                 *
 * ------------------------------------------------------------------ */

/* Deconvolution and NL-means denoise are ROI-capable but declare no halo, so
 * their region previews are computed from truncated context and their edges
 * are approximate.  That is a deliberate stopping point rather than an
 * oversight: both are iterative or non-local, so no exact bound exists (the
 * agreed policy is a flat one radius, documented as an approximation), and
 * enabling one is a visible change to what a preview shows — slower, and
 * different near the rectangle edge — which wants a release note rather than
 * a quiet commit.
 *
 * This pins the current state so that enabling a halo for either is a
 * deliberate act that updates this test, not something that happens by
 * accident. */
Test(roi_halo, ops_without_a_declared_halo_are_listed, .init = setup) {
	const struct { const char *id; const op_descriptor *op; } undeclared[] = {
		{ "filters.deconvolve", &op_desc_deconvolve },
		{ "filters.denoise",    &op_desc_denoise },
	};
	for (size_t i = 0; i < sizeof(undeclared) / sizeof(undeclared[0]); i++) {
		cr_assert(undeclared[i].op->flags & OP_ROI_CAPABLE,
			  "%s is expected to be ROI-capable", undeclared[i].id);
		cr_assert_null(undeclared[i].op->roi_halo,
			       "%s now declares a halo — add it to the measured set above",
			       undeclared[i].id);
	}
}
