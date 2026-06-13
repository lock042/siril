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

#include <criterion/criterion.h>
#include <math.h>
#include <stdlib.h>
#include <glib/gstdio.h>

#include "core/siril.h"
#include "algos/Def_Wavelet.h"
#include "algos/wavelet_denoise.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image (now a pointer)

static void setup(void) {
	com.headless = TRUE;
	com.max_thread = 1; /* pave.c OpenMP pragmas read com.max_thread */
}

TestSuite(wavelet_denoise, .init = setup);

/* Published Starck/Murtagh noise std factors for the 2D B3-spline starlet. */
static const double B3_FACTORS[] = {
	0.8908, 0.2007, 0.0856, 0.0413, 0.0205
};

Test(wavelet_denoise, bspline_factors_match_published) {
	double e[WD_MAX_PLAN];
	int ret = wavelet_noise_factors(TO_PAVE_BSPLINE, 6, e);
	cr_assert_eq(ret, 0, "wavelet_noise_factors failed");

	for (int j = 0; j < 5; j++) {
		/* 1% relative tolerance; the impulse method is exact up to float
		 * roundoff in the kernels, so this is comfortably met. */
		double tol = 0.01 * B3_FACTORS[j];
		cr_assert_float_eq(e[j], B3_FACTORS[j], tol,
				"B3 scale %d factor %.5f != expected %.5f", j, e[j],
				B3_FACTORS[j]);
	}
	/* Factors strictly decrease with scale. */
	for (int j = 1; j < 5; j++)
		cr_assert(e[j] < e[j - 1], "factors not decreasing at scale %d", j);
}

Test(wavelet_denoise, linear_scale0_factor) {
	double e[WD_MAX_PLAN];
	int ret = wavelet_noise_factors(TO_PAVE_LINEAR, 6, e);
	cr_assert_eq(ret, 0, "wavelet_noise_factors failed");

	/* Linear [1,2,1]/4 separable smooth: scale-0 detail = delta - smooth, so
	 * e0^2 = 1 - 2*0.25 + 0.140625 = 0.640625, e0 = 0.80039. */
	cr_assert_float_eq(e[0], 0.80039, 0.01, "linear scale-0 factor %.5f", e[0]);
	for (int j = 1; j < 5; j++)
		cr_assert(e[j] < e[j - 1], "factors not decreasing at scale %d", j);
}

Test(wavelet_denoise, factors_cached_stable) {
	double a[WD_MAX_PLAN], b[WD_MAX_PLAN];
	cr_assert_eq(wavelet_noise_factors(TO_PAVE_BSPLINE, 6, a), 0);
	cr_assert_eq(wavelet_noise_factors(TO_PAVE_BSPLINE, 6, b), 0);
	for (int j = 0; j < 5; j++)
		cr_assert_float_eq(a[j], b[j], 1e-12, "cache mismatch at %d", j);
}

Test(wavelet_denoise, bad_arguments_rejected) {
	double e[WD_MAX_PLAN];
	cr_assert_neq(wavelet_noise_factors(TO_PAVE_BSPLINE, 6, NULL), 0);
	cr_assert_neq(wavelet_noise_factors(99, 6, e), 0);
	cr_assert_neq(wavelet_noise_factors(TO_PAVE_BSPLINE, 1, e), 0);
	cr_assert_neq(wavelet_noise_factors(TO_PAVE_BSPLINE, WD_MAX_PLAN + 1, e), 0);
}

/* Deterministic standard normal sample (Box-Muller). */
static double gauss(void) {
	double u1 = (rand() + 1.0) / ((double) RAND_MAX + 2.0);
	double u2 = (rand() + 1.0) / ((double) RAND_MAX + 2.0);
	return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

Test(wavelet_denoise, mad_sigma_recovers_gaussian) {
	const size_t n = 1000000;
	const double sigma = 10.0;
	float *band = malloc(n * sizeof(float));
	cr_assert_not_null(band);
	srand(42);
	for (size_t i = 0; i < n; i++)
		band[i] = (float) (sigma * gauss());

	double est = wavelet_mad_sigma_float(band, n);
	/* MAD estimator is consistent; 2% is generous for n = 1e6. */
	cr_assert_float_eq(est, sigma, 0.02 * sigma, "MAD sigma %.4f != %.4f", est,
			sigma);
	free(band);
}

/* Standard deviation of (a - b) over a sub-rectangle of an N-wide image. */
static double region_resid_std(const float *a, const float *b, int N, int x0,
		int y0, int w, int h) {
	double s = 0.0, s2 = 0.0;
	int cnt = 0;
	for (int y = y0; y < y0 + h; y++)
		for (int x = x0; x < x0 + w; x++) {
			double d = (double) a[y * N + x] - (double) b[y * N + x];
			s += d;
			s2 += d * d;
			cnt++;
		}
	double mean = s / cnt;
	return sqrt(s2 / cnt - mean * mean);
}

Test(wavelet_denoise, threshold_reconstruction_denoises) {
	const int N = 256;
	const int nplan = 5;
	const double sigma = 0.02;
	const size_t n = (size_t) N * N;

	float *clean = malloc(n * sizeof(float));
	float *noisy = malloc(n * sizeof(float));
	float *out = malloc(n * sizeof(float));
	cr_assert(clean && noisy && out);

	/* flat background + a few Gaussian "stars" near the centre */
	const float bg = 0.1f, amp = 0.6f, psf = 2.5f;
	const int sx[3] = { 120, 140, 130 }, sy[3] = { 120, 130, 145 };
	for (int y = 0; y < N; y++)
		for (int x = 0; x < N; x++) {
			float v = bg;
			for (int s = 0; s < 3; s++) {
				float dx = x - sx[s], dy = y - sy[s];
				v += amp * expf(-(dx * dx + dy * dy) / (2.f * psf * psf));
			}
			clean[y * N + x] = v;
		}
	srand(123);
	for (size_t i = 0; i < n; i++)
		noisy[i] = clean[i] + (float) (sigma * gauss());

	const char *tmpdir = g_get_tmp_dir();
	gchar *fname = g_build_filename(tmpdir, "siril_wd_test.wave", NULL);
	cr_assert_eq(wavelet_transform_file_float(noisy, N, N, fname,
			TO_PAVE_BSPLINE, nplan, 0), 0, "transform failed");

	float coef[7];
	for (int i = 0; i < 7; i++)
		coef[i] = 1.0f;

	/* 1) disabled denoise reconstructs the (noisy) input near-losslessly */
	struct denoise_params dp;
	denoise_params_init(&dp); /* enabled == FALSE */
	cr_assert_eq(wavelet_reconstruct_file_float(fname, coef, &dp, out), 0);
	float maxdiff = 0.f;
	for (size_t i = 0; i < n; i++)
		maxdiff = fmaxf(maxdiff, fabsf(out[i] - noisy[i]));
	cr_assert_lt(maxdiff, 1e-3f, "lossless reconstruction off by %.2e", maxdiff);

	/* 2) enabled denoise (soft threshold, k=3) cuts background noise hard */
	dp.enabled = TRUE;
	dp.method = WD_THRESHOLD;
	cr_assert_eq(wavelet_reconstruct_file_float(fname, coef, &dp, out), 0);

	double noisy_bg = region_resid_std(noisy, clean, N, 8, 8, 48, 48);
	double den_bg = region_resid_std(out, clean, N, 8, 8, 48, 48);
	cr_assert_float_eq(noisy_bg, sigma, 0.15 * sigma,
			"input background noise %.4f != %.4f", noisy_bg, sigma);
	cr_assert_lt(den_bg, 0.6 * noisy_bg,
			"denoise did not reduce background noise (%.4f vs %.4f)", den_bg,
			noisy_bg);

	/* 3) bright signal is preserved (peak within a small margin) */
	float clean_peak = clean[sy[0] * N + sx[0]];
	float den_peak = out[sy[0] * N + sx[0]];
	cr_assert_float_eq(den_peak, clean_peak, 0.12f,
			"star peak not preserved: %.3f vs %.3f", den_peak, clean_peak);

	g_remove(fname);
	g_free(fname);
	free(clean);
	free(noisy);
	free(out);
}

/* Build the synthetic scene (flat bg + 3 Gaussian stars) into clean/noisy. */
static void make_scene(float *clean, float *noisy, int N, double sigma,
		int sx[3], int sy[3]) {
	const float bg = 0.1f, amp = 0.6f, psf = 2.5f;
	for (int y = 0; y < N; y++)
		for (int x = 0; x < N; x++) {
			float v = bg;
			for (int s = 0; s < 3; s++) {
				float dx = x - sx[s], dy = y - sy[s];
				v += amp * expf(-(dx * dx + dy * dy) / (2.f * psf * psf));
			}
			clean[y * N + x] = v;
		}
	srand(123);
	for (size_t i = 0; i < (size_t) N * N; i++)
		noisy[i] = clean[i] + (float) (sigma * gauss());
}

Test(wavelet_denoise, bishrink_reconstruction_denoises) {
	const int N = 256;
	const int nplan = 5;
	const double sigma = 0.02;
	const size_t n = (size_t) N * N;

	float *clean = malloc(n * sizeof(float));
	float *noisy = malloc(n * sizeof(float));
	float *out = malloc(n * sizeof(float));
	cr_assert(clean && noisy && out);

	int sx[3] = { 120, 140, 130 }, sy[3] = { 120, 130, 145 };
	make_scene(clean, noisy, N, sigma, sx, sy);

	const char *tmpdir = g_get_tmp_dir();
	gchar *fname = g_build_filename(tmpdir, "siril_wd_bishrink.wave", NULL);
	cr_assert_eq(wavelet_transform_file_float(noisy, N, N, fname,
			TO_PAVE_BSPLINE, nplan, 0), 0);

	float coef[7];
	for (int i = 0; i < 7; i++)
		coef[i] = 1.0f;

	/* default method is BiShrink */
	struct denoise_params dp;
	denoise_params_init(&dp);
	cr_assert_eq(dp.method, WD_BISHRINK, "BiShrink must be the default method");
	dp.enabled = TRUE;
	cr_assert_eq(wavelet_reconstruct_file_float(fname, coef, &dp, out), 0);

	double noisy_bg = region_resid_std(noisy, clean, N, 8, 8, 48, 48);
	double den_bg = region_resid_std(out, clean, N, 8, 8, 48, 48);
	cr_assert_lt(den_bg, 0.6 * noisy_bg,
			"BiShrink did not reduce background noise (%.4f vs %.4f)", den_bg,
			noisy_bg);

	/* BiShrink preserves bright signal well (parent + local variance keep
	 * strong coefficients). */
	float clean_peak = clean[sy[0] * N + sx[0]];
	float den_peak = out[sy[0] * N + sx[0]];
	cr_assert_float_eq(den_peak, clean_peak, 0.1f,
			"star peak not preserved: %.3f vs %.3f", den_peak, clean_peak);

	g_remove(fname);
	g_free(fname);
	free(clean);
	free(noisy);
	free(out);
}

Test(wavelet_denoise, anscombe_roundtrip_exact) {
	const size_t n = 2048;
	float *a = malloc(n * sizeof(float));
	float *b = malloc(n * sizeof(float));
	cr_assert(a && b);
	for (size_t i = 0; i < n; i++)
		a[i] = b[i] = (float) i / (float) n; /* [0,1) */
	anscombe_forward(b, n, ANSCOMBE_FLOAT_SCALE);
	anscombe_inverse(b, n, ANSCOMBE_FLOAT_SCALE);
	for (size_t i = 0; i < n; i++)
		cr_assert_float_eq(b[i], a[i], 1e-4, "VST round-trip off at %zu: %.6f vs %.6f",
				i, b[i], a[i]);
	free(a);
	free(b);
}

Test(wavelet_denoise, vst_decompose_reconstruct_identity) {
	const int N = 256;
	const int nplan = 5;
	const size_t n = (size_t) N * N;
	float *clean = malloc(n * sizeof(float));
	float *noisy = malloc(n * sizeof(float));
	float *out = malloc(n * sizeof(float));
	cr_assert(clean && noisy && out);

	int sx[3] = { 120, 140, 130 }, sy[3] = { 120, 130, 145 };
	make_scene(clean, noisy, N, 0.02, sx, sy);

	const char *tmpdir = g_get_tmp_dir();
	gchar *fname = g_build_filename(tmpdir, "siril_wd_vst.wave", NULL);
	/* decompose in the VST domain */
	cr_assert_eq(wavelet_transform_file_float(noisy, N, N, fname,
			TO_PAVE_BSPLINE, nplan, 1), 0);

	float coef[7];
	for (int i = 0; i < 7; i++)
		coef[i] = 1.0f;

	/* reconstruct with VST inverse but denoising disabled -> original image */
	struct denoise_params dp;
	denoise_params_init(&dp);
	dp.anscombe = TRUE; /* enabled stays FALSE */
	cr_assert_eq(wavelet_reconstruct_file_float(fname, coef, &dp, out), 0);

	float maxdiff = 0.f;
	for (size_t i = 0; i < n; i++)
		maxdiff = fmaxf(maxdiff, fabsf(out[i] - noisy[i]));
	cr_assert_lt(maxdiff, 2e-3f, "VST round-trip reconstruction off by %.2e", maxdiff);

	g_remove(fname);
	g_free(fname);
	free(clean);
	free(noisy);
	free(out);
}

Test(wavelet_denoise, roi_reconstruct_matches_full) {
	const int N = 256;
	const int nplan = 5;
	const size_t n = (size_t) N * N;
	float *clean = malloc(n * sizeof(float));
	float *noisy = malloc(n * sizeof(float));
	float *full = malloc(n * sizeof(float));
	cr_assert(clean && noisy && full);

	int sx[3] = { 120, 140, 130 }, sy[3] = { 120, 130, 145 };
	make_scene(clean, noisy, N, 0.02, sx, sy);

	const char *tmpdir = g_get_tmp_dir();
	gchar *fname = g_build_filename(tmpdir, "siril_wd_roi.wave", NULL);
	cr_assert_eq(wavelet_transform_file_float(noisy, N, N, fname,
			TO_PAVE_BSPLINE, nplan, 0), 0);

	float coef[7];
	for (int i = 0; i < 7; i++)
		coef[i] = 1.0f;

	struct denoise_params dp;
	denoise_params_init(&dp); /* disabled: pure reconstruction */

	/* full reconstruction (FITS-ordered buffer) */
	cr_assert_eq(wavelet_reconstruct_file_float(fname, coef, &dp, full), 0);

	/* ROI reconstruction into a minimal float fits */
	const int rx = 40, ry = 50, w = 64, h = 48;
	fits roifit = { 0 };
	roifit.type = DATA_FLOAT;
	roifit.rx = roifit.naxes[0] = w;
	roifit.ry = roifit.naxes[1] = h;
	roifit.naxes[2] = 1;
	roifit.naxis = 2;
	roifit.fdata = malloc((size_t) w * h * sizeof(float));
	cr_assert_not_null(roifit.fdata);
	roifit.fpdata[0] = roifit.fdata;

	cr_assert_eq(wavelet_reconstruct_file_roi(fname, coef, &dp, rx, ry, w, h, 0, &roifit), 0);

	/* ROI row y maps to FITS row N-1-ry-y (top-down -> bottom-up) */
	float maxdiff = 0.f;
	for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++) {
			float roiv = roifit.fdata[(size_t) y * w + x];
			float fullv = full[(size_t) (N - 1 - ry - y) * N + (rx + x)];
			maxdiff = fmaxf(maxdiff, fabsf(roiv - fullv));
		}
	cr_assert_lt(maxdiff, 1e-3f, "ROI reconstruction differs from full by %.2e", maxdiff);

	g_remove(fname);
	g_free(fname);
	free(roifit.fdata);
	free(clean);
	free(noisy);
	free(full);
}

Test(wavelet_denoise, estimate_noise_divides_by_factor) {
	const size_t n = 1000000;
	const double sigma_global = 5.0;
	const double e1 = 0.8908;
	float *band = malloc(n * sizeof(float));
	cr_assert_not_null(band);
	srand(7);
	/* finest detail band carries sigma_global * e1 of noise */
	for (size_t i = 0; i < n; i++)
		band[i] = (float) (sigma_global * e1 * gauss());

	double est = wavelet_estimate_noise_float(band, n, e1);
	cr_assert_float_eq(est, sigma_global, 0.02 * sigma_global,
			"recovered sigma %.4f != %.4f", est, sigma_global);
	cr_assert(wavelet_estimate_noise_float(band, n, 0.0) < 0.0,
			"e1=0 must return error");
	free(band);
}
