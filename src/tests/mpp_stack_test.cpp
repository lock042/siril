/*
 * Phase 5a — mpp_stack Criterion tests.
 *
 * Stacking is broken into separately-testable pieces; this file builds out
 * coverage in lock-step with the implementation so divergences from PSS
 * surface at the smallest unit they appear at.
 */
#include <criterion/criterion.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_ap_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_shift_priv.hpp"
#include "registration/mpp/mpp_stack_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_stack.h"
}

cominfo com;
fits *gfit = nullptr;

namespace {

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

}  // namespace

/* ------------------------------------------------------------------------- */

Test(mpp_stack, stub_returns_enotimpl) {
	cr_assert_eq(mpp_stack(nullptr, nullptr, nullptr, nullptr, nullptr,
	                       MPP_RESAMPLE_BICUBIC, 1, nullptr),
	             MPP_ENOTIMPL);
}

Test(mpp_stack, default_config_phase5a) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.alignment_points_frame_percent, 100);
	cr_assert_eq(cfg.alignment_points_frame_number, -1);
	cr_assert_eq(cfg.alignment_points_rank_pixel_stride, 2);
	cr_assert(cfg.alignment_points_de_warp);
	cr_assert_float_eq(cfg.alignment_points_penalty_factor, 0.00025, 1e-12);
	cr_assert_float_eq(cfg.stack_frames_background_fraction, 0.3, 1e-12);
	cr_assert_float_eq(cfg.stack_frames_background_blend_threshold, 0.2, 1e-12);
	cr_assert_eq(cfg.stack_frames_background_patch_size, 100);
	cr_assert_float_eq(cfg.drizzle_scale, 1.0, 1e-12);
}

/* PSS one_dim_weight reference values for a symmetric ramp:
 *   patch=[100, 172], box_center=136, no extend
 *   center_offset = 36; high_len = patch_size - center_offset = 36
 *   weights[0] = 1/37 ≈ 0.027027
 *   weights[35] = 36/37 ≈ 0.972973
 *   weights[36] = 1.0   (the centre — first index of high ramp)
 *   weights[37] = 35/36 ≈ 0.972222
 *   weights[71] = 1/36  ≈ 0.027778
 */
Test(mpp_stack, one_dim_weight_symmetric_ramp) {
	const auto w = mpp::stack_one_dim_weight(100, 172, 136, false, false);
	cr_assert_eq(w.size(), 72u);
	cr_assert_float_eq(w[0],  1.0f / 37.0f, 1e-6);
	cr_assert_float_eq(w[35], 36.0f / 37.0f, 1e-6);
	cr_assert_float_eq(w[36], 1.0f, 1e-6);
	cr_assert_float_eq(w[37], 35.0f / 36.0f, 1e-6);
	cr_assert_float_eq(w[71], 1.0f / 36.0f, 1e-6);
}

Test(mpp_stack, one_dim_weight_extend_low) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 30, true, false);
	cr_assert_eq(w.size(), 50u);
	/* Lower ramp replaced by 1.0 from idx 0 to centre_offset (=30, exclusive). */
	for (int i = 0; i < 30; ++i)
		cr_assert_float_eq(w[i], 1.0f, 1e-9, "w[%d] should be 1.0", i);
	/* High ramp at idx 30 is 1.0; at idx 49 is 1/20. */
	cr_assert_float_eq(w[30], 1.0f, 1e-9);
	cr_assert_float_eq(w[49], 1.0f / 20.0f, 1e-6);
}

Test(mpp_stack, one_dim_weight_extend_high) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 20, false, true);
	cr_assert_eq(w.size(), 50u);
	/* Low ramp: [1/21, 2/21, …, 20/21]. */
	cr_assert_float_eq(w[0], 1.0f / 21.0f, 1e-6);
	cr_assert_float_eq(w[19], 20.0f / 21.0f, 1e-6);
	/* High ramp replaced by 1.0 from idx 20 to end. */
	for (int i = 20; i < 50; ++i)
		cr_assert_float_eq(w[i], 1.0f, 1e-9, "w[%d] should be 1.0", i);
}

/* weight_matrix_first_phase: centre is exactly 1.0; the four corners are
 * 1 − penalty × 2 since (x/sw1 − 1)² + (y/sw1 − 1)² = 1 + 1 = 2 at
 * (0,0). For default search_width=14: sw1 = 5, extent = 11. */
Test(mpp_stack, weight_matrix_centre_and_corners) {
	const auto cfg = default_cfg();
	const cv::Mat m = mpp::stack_build_first_phase_weight_matrix(cfg);
	cr_assert_eq(m.rows, 11);
	cr_assert_eq(m.cols, 11);
	cr_assert_float_eq(m.at<float>(5, 5), 1.0f, 1e-9);

	const double pen = cfg.alignment_points_penalty_factor;
	const float expected_corner = (float) (1.0 - pen * 2.0);
	cr_assert_float_eq(m.at<float>(0, 0),  expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(0, 10), expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(10, 0), expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(10, 10), expected_corner, 1e-7);

	/* (0, 5): x at centre (no x penalty), y at 0 (full y penalty of 1). */
	const float expected_edge = (float) (1.0 - pen * 1.0);
	cr_assert_float_eq(m.at<float>(0, 5), expected_edge, 1e-7);
}

/* remap_rigid: a 10×10 patch shifted (0, 0) lands exactly. */
Test(mpp_stack, remap_rigid_no_shift) {
	cv::Mat frame(50, 50, CV_32F, cv::Scalar(0));
	for (int y = 20; y < 30; ++y)
		for (int x = 20; x < 30; ++x)
			frame.at<float>(y, x) = (float) (y * 100 + x);
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 20, 30, 20, 30, b);

	for (int y = 0; y < 10; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x),
			                   (float) ((y + 20) * 100 + (x + 20)), 1e-6);
	cr_assert_eq(b.y_low,  0); cr_assert_eq(b.y_high, 0);
	cr_assert_eq(b.x_low,  0); cr_assert_eq(b.x_high, 0);
}

/* remap_rigid accumulates: calling twice doubles the contribution. */
Test(mpp_stack, remap_rigid_accumulates) {
	cv::Mat frame(20, 20, CV_32F, cv::Scalar(5.0f));
	cv::Mat buffer(8, 8, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 5, 13, 5, 13, b);
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 5, 13, 5, 13, b);
	cr_assert_float_eq(buffer.at<float>(4, 4), 10.0f, 1e-6);
}


/* remap_rigid clips when the shifted patch reaches off the frame and
 * records the maximum cropped extent. */
Test(mpp_stack, remap_rigid_clips_and_records_border) {
	cv::Mat frame(20, 20, CV_32F);
	for (int y = 0; y < 20; ++y)
		for (int x = 0; x < 20; ++x)
			frame.at<float>(y, x) = (float) (y * 100 + x);
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	/* Patch [0,10]×[0,10] with shift (-3, -3): src starts at (-3, -3), so
	 * the top-left 3 rows × 3 cols of the buffer remain zero (the clipped
	 * portion). Border counters record the clipped extents. */
	mpp::stack_remap_rigid(frame, buffer, -3, -3, 0, 10, 0, 10, b);
	cr_assert_eq(b.y_low, 3);
	cr_assert_eq(b.x_low, 3);
	cr_assert_eq(b.y_high, 0);
	cr_assert_eq(b.x_high, 0);
	for (int y = 0; y < 3; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x), 0.0f, 1e-9);
	/* The non-clipped region holds frame[y-3, x-3] for y ≥ 3, x ≥ 3. */
	cr_assert_float_eq(buffer.at<float>(3, 3), frame.at<float>(0, 0), 1e-6);
	cr_assert_float_eq(buffer.at<float>(9, 9), frame.at<float>(6, 6), 1e-6);
}
