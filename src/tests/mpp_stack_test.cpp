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
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_ap_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_shift_priv.hpp"
#include "registration/mpp/mpp_stack_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_shift.h"
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

cv::Mat make_textured_frame(int rows = 200, int cols = 240, int seed = 1) {
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
	rng.state = seed;
	for (int i = 0; i < 40; ++i) {
		const int sy = rng.uniform(int(cy - r0 * 0.6), int(cy + r0 * 0.6));
		const int sx = rng.uniform(int(cx - r0 * 0.6), int(cx + r0 * 0.6));
		cv::circle(f, {sx, sy}, 2, cv::Scalar(60000), -1);
	}
	return f;
}

cv::Mat blurred(const cv::Mat &mono, const mpp_config_t &cfg) {
	cv::Mat b;
	cv::GaussianBlur(mono, b, cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	return b;
}

cv::Mat shifted(const cv::Mat &mono, int dy, int dx) {
	const cv::Mat M = (cv::Mat_<double>(2, 3) << 1, 0, dx, 0, 1, dy);
	cv::Mat out;
	cv::warpAffine(mono, out, M, mono.size(), cv::INTER_NEAREST,
	               cv::BORDER_REPLICATE);
	return out;
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

/* The per-frame AP accumulation is parallelised over alignment points
 * (disjoint per-AP buffers + an order-independent max-merge of the border).
 * Running the full apply pass single-threaded and multi-threaded over the
 * same in-memory frames must produce bit-identical buffers, borders, and
 * final merged image. */
Test(mpp_stack, apply_shifts_parallel_matches_serial) {
	const auto cfg = default_cfg();
	const cv::Mat truth = make_textured_frame();
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -1}, {1, 2}, {-1, -2}, {2, 3}};
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q;
	for (auto p : jit) {
		const cv::Mat f = shifted(truth, p.first, p.second);
		frames_raw.push_back(f);
		frames_blurred.push_back(blurred(f, cfg));
		q.push_back(0.5);
	}
	q[0] = 1.0;  /* frame 0 is the reference */

	const auto align = mpp::align_global_from_frames(frames_blurred, q, cfg);
	auto cfg_no_fc = cfg;
	cfg_no_fc.align_frames_fast_changing_object = false;
	const auto avg = mpp::align_average_frame(frames_raw, q, align.shifts, cfg_no_fc);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	cr_assert_gt(aps->count, 1);

	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);

	std::vector<double> brightness;
	for (const auto &f : frames_raw)
		brightness.push_back(mpp::rank_average_brightness(f, cfg));

	const auto apq = mpp::ap_compute_frame_qualities(
	    frames_raw, brightness, *aps, offsets, truth.rows, truth.cols, cfg);
	/* Default frame_percent = 100 ⇒ every frame is in every AP's top-N, so
	 * the inner loop spans all M APs — the heavy case we're parallelising. */
	cr_assert_eq(apq.stack_size, (int) frames_raw.size());

	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	mpp_shifts_t *shifts = mpp::shift_compute_all(frames_blurred, *aps,
	                                              ref_boxes, offsets, cfg);
	cr_assert_not_null(shifts);

	std::vector<int> sorted(frames_raw.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(),
	          [&q](int a, int b) { return q[a] > q[b]; });

	const auto serial = mpp::stack_apply_shifts(frames_raw, *aps, apq, shifts,
	                                            offsets, brightness, sorted,
	                                            avg.intersection, cfg, nullptr,
	                                            /*max_threads=*/1);
	const auto par = mpp::stack_apply_shifts(frames_raw, *aps, apq, shifts,
	                                         offsets, brightness, sorted,
	                                         avg.intersection, cfg, nullptr,
	                                         /*max_threads=*/8);

	cr_assert_eq(serial.border.y_low,  par.border.y_low);
	cr_assert_eq(serial.border.y_high, par.border.y_high);
	cr_assert_eq(serial.border.x_low,  par.border.x_low);
	cr_assert_eq(serial.border.x_high, par.border.x_high);
	cr_assert_eq(serial.shift_failure_counter, par.shift_failure_counter);

	cr_assert_eq((int) serial.state.stacking_buffers.size(),
	             (int) par.state.stacking_buffers.size());
	for (int a = 0; a < (int) serial.state.stacking_buffers.size(); ++a) {
		const cv::Mat &s = serial.state.stacking_buffers[a];
		const cv::Mat &p = par.state.stacking_buffers[a];
		cr_assert_eq(s.rows, p.rows);
		cr_assert_eq(s.cols, p.cols);
		cr_assert_eq(s.type(), p.type());
		const double d = cv::norm(s, p, cv::NORM_INF);
		cr_assert_float_eq(d, 0.0, 0.0,
		    "AP %d stacking buffer differs (NORM_INF=%.6g)", a, d);
	}

	const cv::Mat ms = mpp::stack_merge_alignment_point_buffers(
	    serial.state, serial.border, *aps, cfg);
	const cv::Mat mp = mpp::stack_merge_alignment_point_buffers(
	    par.state, par.border, *aps, cfg);
	cr_assert_eq(ms.rows, mp.rows);
	cr_assert_eq(ms.cols, mp.cols);
	cr_assert_eq(ms.type(), mp.type());
	cr_assert_eq(0, std::memcmp(ms.data, mp.data, ms.total() * ms.elemSize()));

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}
