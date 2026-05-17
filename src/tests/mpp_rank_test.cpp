/*
 * Phase 1 — mpp_rank Criterion tests.
 *
 * Covers PSS-faithful Laplace-σ ranking against the private C++ entry points,
 * plus the C-level mpp_rank_sequence stub (NB: sequence integration ships in
 * Phase 1.3; the stub still returns ENOTIMPL).
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_rank_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_config.h"
#include "registration/mpp_rank.h"
}

cominfo com;
fits *gfit = nullptr;

namespace {

/* Synthesise a textured 16-bit mono frame with a bright disk + speckles, so
 * the Laplacian has plenty of structure to measure. Deterministic seed. */
cv::Mat make_sharp_frame(int rows = 200, int cols = 240, int seed = 1) {
	cv::Mat f(rows, cols, CV_16U, cv::Scalar(0));
	const double cy = rows / 2.0, cx = cols / 2.0;
	const double r0 = std::min(rows, cols) * 0.30;
	cv::RNG rng(seed);
	for (int y = 0; y < rows; ++y) {
		for (int x = 0; x < cols; ++x) {
			const double r = std::hypot(y - cy, x - cx);
			double v = 0.0;
			if (r < r0) {
				v = 0.55 + 0.25 * std::tanh((r0 - r) / 4.0);
				v += 0.03 * rng.gaussian(1.0);
			}
			v = std::clamp(v, 0.0, 1.0);
			f.at<uint16_t>(y, x) = static_cast<uint16_t>(v * 65535.0);
		}
	}
	/* Sprinkle bright speckles to feed the Laplacian. */
	rng.state = seed;
	for (int i = 0; i < 40; ++i) {
		const int sy = rng.uniform(int(cy - r0 * 0.6), int(cy + r0 * 0.6));
		const int sx = rng.uniform(int(cx - r0 * 0.6), int(cx + r0 * 0.6));
		cv::circle(f, {sx, sy}, 2, cv::Scalar(60000), -1);
	}
	return f;
}

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

}  // namespace

/* ------------------------------------------------------------------------- */

Test(mpp_rank, stub_sequence_returns_enotimpl) {
	double *q = nullptr;
	cr_assert_eq(mpp_rank_sequence(nullptr, nullptr, &q), MPP_ENOTIMPL);
	cr_assert_null(q);
}

Test(mpp_rank, default_config_has_pss_values) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.frames_gauss_width, 7);
	cr_assert_eq(cfg.align_frames_sampling_stride, 2);
	cr_assert_float_eq(cfg.rank_laplacian_alpha, 1.0 / 256.0, 1e-12);
	cr_assert(cfg.frames_normalization);
	cr_assert_eq(cfg.frames_normalization_threshold, 15);
	cr_assert_eq(cfg.bitdepth, 16);
	cr_assert_float_eq(mpp_cfg_threshold_scale(&cfg), 256.0, 1e-12);
}

Test(mpp_rank, fits_bitpix_maps_correctly) {
	cr_assert_eq(mpp_bitdepth_from_fits_bitpix(8),   8);   /* BYTE_IMG */
	cr_assert_eq(mpp_bitdepth_from_fits_bitpix(16),  16);  /* SHORT_IMG (signed) */
	cr_assert_eq(mpp_bitdepth_from_fits_bitpix(20),  16);  /* USHORT_IMG — Siril's normal 16-bit */
	cr_assert_eq(mpp_bitdepth_from_fits_bitpix(-32), 16);  /* FLOAT_IMG */
}

Test(mpp_rank, score_is_positive_on_sharp_frame) {
	const auto cfg = default_cfg();
	const cv::Mat f = make_sharp_frame();
	const double s = mpp::rank_score_mat(f, cfg);
	cr_assert_gt(s, 0.0, "σ should be > 0 on a textured frame, got %.6g", s);
}

/* PSS sanity: blurring a frame more must reduce its Laplace-σ. */
Test(mpp_rank, blur_decreases_quality) {
	const auto cfg = default_cfg();
	const cv::Mat sharp = make_sharp_frame();
	cv::Mat soft;
	cv::GaussianBlur(sharp, soft, cv::Size(13, 13), 0);  /* extra blur beyond the 7×7 PSS uses */

	const double q_sharp = mpp::rank_score_mat(sharp, cfg);
	const double q_soft  = mpp::rank_score_mat(soft,  cfg);

	cr_assert_gt(q_sharp, q_soft,
	             "sharp σ (%.4f) should exceed softened σ (%.4f)",
	             q_sharp, q_soft);
}

/* PSS sanity: the normalized rank should be approximately invariant to a
 * uniform scaling of the input. The raw σ does scale with brightness, but
 * `rank_score_normalized` divides it back out. We test with `THRESH_TOZERO`-
 * style background to mimic PSS's brightness measure. */
Test(mpp_rank, brightness_normalization_is_scale_invariant) {
	auto cfg = default_cfg();
	cfg.frames_normalization = true;

	const cv::Mat f = make_sharp_frame();
	cv::Mat dim;
	f.convertTo(dim, CV_16U, 0.5);  /* halve every pixel */

	const double q_full = mpp::rank_score_normalized(f,   cfg);
	const double q_dim  = mpp::rank_score_normalized(dim, cfg);

	/* Quantisation through convertScaleAbs+uint8 means the equivalence isn't
	 * exact; we just want the dynamic range to be tamed compared to the raw σ. */
	const double q_full_raw = mpp::rank_score_mat(f,   cfg);
	const double q_dim_raw  = mpp::rank_score_mat(dim, cfg);

	const double normalized_rel = std::abs(q_full - q_dim) / q_full;
	const double raw_rel        = std::abs(q_full_raw - q_dim_raw) / q_full_raw;

	cr_assert_lt(normalized_rel, raw_rel,
	             "normalization should reduce the brightness sensitivity "
	             "(normalized Δrel=%.4f, raw Δrel=%.4f)",
	             normalized_rel, raw_rel);
}

/* PSS uses 16-bit input by default; threshold scales by 256 in that case. */
Test(mpp_rank, brightness_threshold_scales_with_bitdepth) {
	const cv::Mat f = make_sharp_frame();

	auto cfg16 = default_cfg();
	cfg16.bitdepth = 16;
	cfg16.frames_normalization_threshold = 15;
	const double ab16 = mpp::rank_average_brightness(f, cfg16);

	auto cfg8 = default_cfg();
	cfg8.bitdepth = 8;
	cfg8.frames_normalization_threshold = 15;
	const double ab8 = mpp::rank_average_brightness(f, cfg8);

	/* For 16-bit input, thr16 = 15 × 256 = 3840 in the 16-bit range; thr8 = 15
	 * (essentially no rejection at 16-bit scale). The 8-bit-threshold reading
	 * should therefore be ≥ the 16-bit-threshold reading. */
	cr_assert_geq(ab8, ab16,
	              "ab(thr=15)=%.2f should be ≥ ab(thr=3840)=%.2f", ab8, ab16);
}

/* 8-bit Siril SER data is stored as WORD with values 0..255. align_average_frame
 * must upscale these by 256/N so mean_frame lands in the 16-bit-equivalent
 * range that downstream AP code assumes (matching PSS frames.average_frame). */
Test(mpp_rank, average_frame_upscales_8bit_input) {
	/* This test lives here because it exercises the same bitdepth knob; the
	 * function under test is mpp::align_average_frame. */
	auto cfg = default_cfg();
	/* Build a CV_16U-stored "8-bit-range" frame: values in 0..255. */
	cv::Mat truth(64, 80, CV_16U);
	for (int y = 0; y < truth.rows; ++y)
		for (int x = 0; x < truth.cols; ++x)
			truth.at<uint16_t>(y, x) = static_cast<uint16_t>((x + y) % 256);

	std::vector<cv::Mat> frames = {truth, truth, truth, truth};
	std::vector<double> q = {1.0, 0.9, 0.8, 0.7};
	std::vector<cv::Vec2i> shifts(4, cv::Vec2i(0, 0));

	cfg.bitdepth = 8;
	cfg.align_frames_fast_changing_object = false;
	const auto r8 = mpp::align_average_frame(frames, q, shifts, cfg);
	cfg.bitdepth = 16;
	const auto r16 = mpp::align_average_frame(frames, q, shifts, cfg);
	/* For 8-bit input, mean = 256 * sum/N = 256 * value (since all frames
	 * are the same). For 16-bit input, mean = sum/N = value. So r8 should
	 * be 256× r16 cell-for-cell. */
	cr_assert_eq(r8.mean_frame.size(), r16.mean_frame.size());
	int worst = 0;
	for (int y = 0; y < r8.mean_frame.rows; ++y)
		for (int x = 0; x < r8.mean_frame.cols; ++x) {
			const int v8  = r8.mean_frame.at<int>(y, x);
			const int v16 = r16.mean_frame.at<int>(y, x);
			worst = std::max(worst, std::abs(v8 - 256 * v16));
		}
	cr_assert_lt(worst, 256,
	             "8-bit mean_frame should be ~256× the 16-bit one; worst Δ = %d",
	             worst);
}
