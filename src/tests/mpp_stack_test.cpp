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

/* one_dim_weight: PSS's linear ramp t is mapped through the raised-cosine
 * sin²(π/2·t) (a deliberate divergence — C1 blend window, see
 * stack_one_dim_weight), so the expected values are hann(t) of the PSS
 * reference ramp. For the symmetric case:
 *   patch=[100, 172], box_center=136, no extend
 *   center_offset = 36; high_len = patch_size - center_offset = 36
 *   weights[0]  = hann(1/37)
 *   weights[35] = hann(36/37)
 *   weights[36] = hann(1) = 1.0   (the centre — first index of high ramp)
 *   weights[37] = hann(35/36)
 *   weights[71] = hann(1/36)
 */
namespace {
float hann_of(double t) {
	const double s = std::sin(M_PI_2 * t);
	return (float) (s * s);
}
}  // namespace

Test(mpp_stack, one_dim_weight_symmetric_ramp) {
	const auto w = mpp::stack_one_dim_weight(100, 172, 136, false, false);
	cr_assert_eq(w.size(), 72u);
	cr_assert_float_eq(w[0],  hann_of(1.0 / 37.0), 1e-6);
	cr_assert_float_eq(w[35], hann_of(36.0 / 37.0), 1e-6);
	cr_assert_float_eq(w[36], 1.0f, 1e-6);
	cr_assert_float_eq(w[37], hann_of(35.0 / 36.0), 1e-6);
	cr_assert_float_eq(w[71], hann_of(1.0 / 36.0), 1e-6);
	/* Zero slope at the patch edge: the first step of the Hann taper is
	 * much smaller than the first step of the linear ramp it replaced. */
	cr_assert_lt(w[1] - w[0], 1.0f / 37.0f);
}

Test(mpp_stack, one_dim_weight_extend_low) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 30, true, false);
	cr_assert_eq(w.size(), 50u);
	/* Lower taper replaced by 1.0 from idx 0 to centre_offset (=30, exclusive). */
	for (int i = 0; i < 30; ++i)
		cr_assert_float_eq(w[i], 1.0f, 1e-9, "w[%d] should be 1.0", i);
	/* High taper at idx 30 is 1.0; at idx 49 is hann(1/20). */
	cr_assert_float_eq(w[30], 1.0f, 1e-9);
	cr_assert_float_eq(w[49], hann_of(1.0 / 20.0), 1e-6);
}

Test(mpp_stack, one_dim_weight_extend_high) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 20, false, true);
	cr_assert_eq(w.size(), 50u);
	/* Low taper: hann of [1/21, 2/21, …, 20/21]. */
	cr_assert_float_eq(w[0], hann_of(1.0 / 21.0), 1e-6);
	cr_assert_float_eq(w[19], hann_of(20.0 / 21.0), 1e-6);
	/* High taper replaced by 1.0 from idx 20 to end. */
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

/* Selection weight: hard cliff at taper 0; plateau / raised-cosine ramp /
 * zero zones at taper > 0; the half-sample-centred cosine sums to exactly
 * stack_size. */
Test(mpp_stack, selection_weight_profile) {
	/* Hard top-N. */
	cr_assert_float_eq(mpp::stack_selection_weight(0, 5, 0), 1.0f, 1e-9);
	cr_assert_float_eq(mpp::stack_selection_weight(4, 5, 0), 1.0f, 1e-9);
	cr_assert_float_eq(mpp::stack_selection_weight(5, 5, 0), 0.0f, 1e-9);

	/* N = 10, T = 3: plateau ends at rank 7, zero from rank 13. */
	const int N = 10, T = 3;
	for (int r = 0; r < 7; ++r)
		cr_assert_float_eq(mpp::stack_selection_weight(r, N, T), 1.0f, 1e-9,
		                   "rank %d should be on the plateau", r);
	cr_assert_float_eq(mpp::stack_selection_weight(13, N, T), 0.0f, 1e-9);
	/* Monotone decrease through the ramp, symmetric about rank N. */
	float prev = 1.0f;
	double sum = 7.0;
	for (int r = 7; r < 13; ++r) {
		const float w = mpp::stack_selection_weight(r, N, T);
		cr_assert_lt(w, prev);
		cr_assert_gt(w, 0.0f);
		prev = w;
		sum += w;
	}
	cr_assert_float_eq((float) sum, (float) N, 1e-5,
	                   "weights must sum to stack_size (%f)", sum);
	/* Complementary pairs across the centre sum to 1. */
	cr_assert_float_eq(mpp::stack_selection_weight(7, N, T)
	                   + mpp::stack_selection_weight(12, N, T), 1.0f, 1e-6);
}

/* Soft selection wiring: at 50% of 6 frames the target is 3 and the taper
 * clamps to 1 (≤ N/2), so each AP keeps 4 frames with weights
 * [1, 1, hann, hann]. Identical input frames make the per-AP ranking a
 * pure index tie-break, so the expected lists are deterministic. */
Test(mpp_stack, soft_selection_taper_and_weights) {
	auto cfg = default_cfg();
	cfg.alignment_points_frame_percent = 50;
	const cv::Mat truth = make_textured_frame();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q;
	for (int i = 0; i < 6; ++i) {
		frames_raw.push_back(truth.clone());
		frames_blurred.push_back(blurred(truth, cfg));
		q.push_back(0.5);
	}
	q[0] = 1.0;

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
	cr_assert_eq(apq.stack_size, 3);
	cr_assert_eq(apq.taper, 1);
	for (int a = 0; a < aps->count; ++a)
		cr_assert_eq((int) apq.best_frame_indices[a].size(), 4);

	/* Identical frames ⇒ ranking falls back to the index tie-break, so
	 * frame k sits at rank k in every AP. */
	const float w2 = mpp::stack_selection_weight(2, 3, 1);
	const float w3 = mpp::stack_selection_weight(3, 3, 1);
	cr_assert_float_eq(w2 + w3, 1.0f, 1e-6);
	for (int f = 0; f < 6; ++f) {
		const auto &uses = apq.used_alignment_points[f];
		if (f <= 3) {
			cr_assert_eq((int) uses.size(), aps->count,
			             "frame %d should be used by every AP", f);
			const float expect = f < 2 ? 1.0f : (f == 2 ? w2 : w3);
			for (const auto &u : uses)
				cr_assert_float_eq(u.weight, expect, 1e-6,
				                   "frame %d weight %f != %f",
				                   f, u.weight, expect);
		} else {
			cr_assert_eq((int) uses.size(), 0,
			             "frame %d is past the taper and must be unused", f);
		}
	}

	/* Normalisation invariance: identical frames mean ANY selection
	 * weighting must produce the same stacked image. Compare the 50%
	 * soft-selection stack against the 100% full stack. */
	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	mpp_shifts_t *shifts = mpp::shift_compute_all(frames_blurred, *aps,
	                                              ref_boxes, offsets, cfg);
	cr_assert_not_null(shifts);
	std::vector<int> sorted(frames_raw.size());
	std::iota(sorted.begin(), sorted.end(), 0);

	auto cfg_full = cfg;
	cfg_full.alignment_points_frame_percent = 100;
	const auto apq_full = mpp::ap_compute_frame_qualities(
	    frames_raw, brightness, *aps, offsets, truth.rows, truth.cols, cfg_full);
	cr_assert_eq(apq_full.taper, 0);

	const auto soft = mpp::stack_apply_shifts(frames_raw, *aps, apq, shifts,
	                                          offsets, brightness, sorted,
	                                          avg.intersection, cfg, nullptr, 1);
	const auto full = mpp::stack_apply_shifts(frames_raw, *aps, apq_full, shifts,
	                                          offsets, brightness, sorted,
	                                          avg.intersection, cfg_full, nullptr, 1);
	const cv::Mat ms = mpp::stack_merge_alignment_point_buffers(
	    soft.state, soft.border, *aps, cfg);
	const cv::Mat mf = mpp::stack_merge_alignment_point_buffers(
	    full.state, full.border, *aps, cfg_full);
	cr_assert_eq(ms.rows, mf.rows);
	cr_assert_eq(ms.cols, mf.cols);
	cv::Mat diff;
	cv::absdiff(ms, mf, diff);
	double dmin = 0.0, dmax = 0.0;
	cv::minMaxLoc(diff.reshape(1), &dmin, &dmax);
	cr_assert_leq(dmax, 2.0,
	              "weighted vs unweighted stack of identical frames diverges "
	              "by %.0f LSB", dmax);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}

/* remap_subpixel with an (effectively) integer shift must take the exact
 * integer blit — identical buffer and border to stack_remap_rigid. */
Test(mpp_stack, remap_subpixel_integer_falls_back_to_rigid) {
	cv::Mat frame(40, 40, CV_32F);
	for (int y = 0; y < 40; ++y)
		for (int x = 0; x < 40; ++x)
			frame.at<float>(y, x) = (float) (y * 100 + x);
	cv::Mat b_rigid(12, 12, CV_32F, cv::Scalar(0));
	cv::Mat b_sub(12, 12, CV_32F, cv::Scalar(0));
	mpp::RemapBorder br, bs;
	mpp::stack_remap_rigid(frame, b_rigid, -3, 2, 0, 12, 5, 17, br);
	mpp::stack_remap_subpixel(frame, b_sub, -3.0 + 1e-7, 2.0 - 1e-7,
	                          0, 12, 5, 17, bs);
	cv::Mat diff;
	cv::absdiff(b_rigid, b_sub, diff);
	double dmin = 0.0, dmax = 0.0;
	cv::minMaxLoc(diff, &dmin, &dmax);
	cr_assert_eq(dmax, 0.0, "integer-shift fast path must be exact (%g)", dmax);
	cr_assert_eq(bs.y_low, br.y_low);
	cr_assert_eq(bs.x_low, br.x_low);
}

/* A constant frame must stay constant under a fractional shift (the Lanczos
 * kernel is normalised), and the accumulate contract must hold. */
Test(mpp_stack, remap_subpixel_constant_preserved) {
	cv::Mat frame(40, 40, CV_32F, cv::Scalar(5.0f));
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_subpixel(frame, buffer, 0.5, 0.25, 10, 20, 10, 20, b);
	mpp::stack_remap_subpixel(frame, buffer, 0.5, 0.25, 10, 20, 10, 20, b);
	for (int y = 0; y < 10; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x), 10.0f, 1e-3,
			                   "buffer[%d][%d] = %f", y, x,
			                   buffer.at<float>(y, x));
	cr_assert_eq(b.y_low, 0); cr_assert_eq(b.y_high, 0);
	cr_assert_eq(b.x_low, 0); cr_assert_eq(b.x_high, 0);
}

/* A half-pixel shift of a linear ramp lands halfway between samples
 * (Lanczos reproduces affine signals to well under 1% in the interior). */
Test(mpp_stack, remap_subpixel_half_pixel_ramp) {
	cv::Mat frame(40, 40, CV_32F);
	for (int y = 0; y < 40; ++y)
		for (int x = 0; x < 40; ++x)
			frame.at<float>(y, x) = (float) (10 * x);
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_subpixel(frame, buffer, 0.0, 0.5, 15, 25, 15, 25, b);
	for (int y = 0; y < 10; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x),
			                   (float) (10.0 * (x + 15) + 5.0), 0.5,
			                   "buffer[%d][%d] = %f", y, x,
			                   buffer.at<float>(y, x));
}

/* Fractional shift past the frame edge: clipped rows are recorded in the
 * border counters so the final trim removes them. */
Test(mpp_stack, remap_subpixel_clips_and_records_border) {
	cv::Mat frame(20, 20, CV_32F, cv::Scalar(7.0f));
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	/* Sample rows start at -2.5 → first 3 destination rows clipped. */
	mpp::stack_remap_subpixel(frame, buffer, -2.5, 0.5, 0, 10, 5, 15, b);
	cr_assert_eq(b.y_low, 3);
	cr_assert_eq(b.y_high, 0);
	for (int x = 0; x < 10; ++x) {
		cr_assert_float_eq(buffer.at<float>(0, x), 0.0f, 1e-9);
		cr_assert_float_eq(buffer.at<float>(2, x), 0.0f, 1e-9);
		cr_assert_float_eq(buffer.at<float>(3, x), 7.0f, 1e-3);
	}
}

/* The apply pass is frame-parallel: each thread accumulates its frames into
 * private per-AP buffers that are then summed. The float-sum reduction
 * reorders, so single-threaded and multi-threaded runs over the same frames
 * match to within ≤1 16-bit LSB (not bit-identical) — verified on the merged
 * image. The border (a max-reduction) stays exact. */
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

	/* The frame-parallel per-thread-buffer reduction reorders the float sum,
	 * so the merged image is not bit-identical across thread counts. The
	 * reorder is ≤1 16-bit LSB on representative data; a broken reduction
	 * would diverge by thousands. Tolerance 2 leaves headroom without
	 * masking breakage. */
	const cv::Mat ms = mpp::stack_merge_alignment_point_buffers(
	    serial.state, serial.border, *aps, cfg);
	const cv::Mat mp = mpp::stack_merge_alignment_point_buffers(
	    par.state, par.border, *aps, cfg);
	cr_assert_eq(ms.rows, mp.rows);
	cr_assert_eq(ms.cols, mp.cols);
	cr_assert_eq(ms.type(), mp.type());
	cv::Mat diff;
	cv::absdiff(ms, mp, diff);
	double dmin = 0.0, dmax = 0.0;
	cv::minMaxLoc(diff.reshape(1), &dmin, &dmax);
	cr_assert_leq(dmax, 2.0,
	              "parallel vs serial merged image diverges by %.0f LSB", dmax);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}

/* Fractional output scale (e.g. 1.5x) is honoured exactly — the drizzled
 * canvas is round(intersection × scale), NOT the nearest integer factor —
 * and accumulation stays deterministic across thread counts. */
Test(mpp_stack, apply_shifts_fractional_scale) {
	auto cfg = default_cfg();
	cfg.drizzle_scale = 1.5;
	const cv::Mat truth = make_textured_frame();
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -1}, {1, 2}};
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q;
	for (auto p : jit) {
		const cv::Mat f = shifted(truth, p.first, p.second);
		frames_raw.push_back(f);
		frames_blurred.push_back(blurred(f, cfg));
		q.push_back(0.5);
	}
	q[0] = 1.0;

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
	                                            avg.intersection, cfg, nullptr, 1);
	const auto par = mpp::stack_apply_shifts(frames_raw, *aps, apq, shifts,
	                                         offsets, brightness, sorted,
	                                         avg.intersection, cfg, nullptr, 8);

	const int dim_y = avg.intersection[1] - avg.intersection[0];
	const int dim_x = avg.intersection[3] - avg.intersection[2];
	cr_assert_float_eq(serial.state.drizzle_scale, 1.5, 1e-12);
	/* Exact fractional geometry, not integer-rounded scale. */
	cr_assert_eq(serial.state.dim_y_drizzled, (int) std::lround(dim_y * 1.5));
	cr_assert_eq(serial.state.dim_x_drizzled, (int) std::lround(dim_x * 1.5));
	cr_assert_neq(serial.state.dim_y_drizzled, dim_y * 2);

	/* 1-thread vs 8-thread stays close at fractional scale (frame-parallel
	 * reduction reorders the float sum → ≤1 16-bit LSB, not bit-identical). */
	cr_assert_eq((int) serial.state.stacking_buffers.size(),
	             (int) par.state.stacking_buffers.size());
	const cv::Mat ms = mpp::stack_merge_alignment_point_buffers(
	    serial.state, serial.border, *aps, cfg);
	const cv::Mat mp = mpp::stack_merge_alignment_point_buffers(
	    par.state, par.border, *aps, cfg);
	cr_assert_eq(ms.rows, mp.rows);
	cr_assert_eq(ms.cols, mp.cols);
	cv::Mat diff;
	cv::absdiff(ms, mp, diff);
	double dmin = 0.0, dmax = 0.0;
	cv::minMaxLoc(diff.reshape(1), &dmin, &dmax);
	cr_assert_leq(dmax, 2.0,
	              "fractional-scale parallel vs serial diverges by %.0f LSB", dmax);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}

/* Per-AP DC equalisation: two flat AP buffers at different levels must merge
 * to a flat image at the common mean instead of a ramp between the levels.
 * With two APs the pairwise least-squares solve is exact: the overlap mean
 * difference is the full DC gap, the zero-mean constraint splits it evenly,
 * and both corrected buffers land on the same constant. */
Test(mpp_stack, merge_equalises_ap_dc_offsets) {
	const auto cfg = default_cfg();

	/* Hand-built two-AP geometry on a 72×126 canvas: both patches span the
	 * full height; AP 0 covers x [0, 72), AP 1 covers x [54, 126) — an
	 * 18-column overlap, like neighbouring grid APs. */
	mpp_ap_record_t recs[2];
	std::memset(recs, 0, sizeof(recs));
	recs[0].y = 36;          recs[0].x = 36;
	recs[0].box_y_low = 12;  recs[0].box_y_high = 60;
	recs[0].box_x_low = 12;  recs[0].box_x_high = 60;
	recs[0].patch_y_low = 0; recs[0].patch_y_high = 72;
	recs[0].patch_x_low = 0; recs[0].patch_x_high = 72;
	recs[1] = recs[0];
	recs[1].x = 90;
	recs[1].box_x_low = 66;  recs[1].box_x_high = 114;
	recs[1].patch_x_low = 54; recs[1].patch_x_high = 126;

	mpp_aps_t aps;
	std::memset(&aps, 0, sizeof(aps));
	aps.count = 2;
	aps.records = recs;

	const cv::Vec4i intersection(0, 72, 0, 126);
	auto state = mpp::stack_prepare_for_blending(aps, intersection,
	                                             /*stack_size=*/1,
	                                             /*drizzle_scale=*/1.0,
	                                             /*num_layers=*/1, cfg);
	cr_assert_eq(state.number_stacking_holes, 0);
	cr_assert_eq((int) state.ap_frame_counts.size(), 2);

	state.stacking_buffers[0].setTo(cv::Scalar(10000.0f));
	state.stacking_buffers[1].setTo(cv::Scalar(12000.0f));

	mpp::RemapBorder border;
	const cv::Mat img = mpp::stack_merge_alignment_point_buffers(
	    state, border, aps, cfg);

	double dmin = 0.0, dmax = 0.0;
	cv::minMaxLoc(img, &dmin, &dmax);
	/* Without equalisation this is a 10000 → 12000 blend ramp; with it,
	 * both buffers shift to the common mean and the output is flat. */
	cr_assert_leq(dmax - dmin, 1.0,
	              "merged image not flat: min=%.0f max=%.0f", dmin, dmax);
	cr_assert_geq(dmin, 10998.0);
	cr_assert_leq(dmax, 11002.0);
}

/* stack_skip_failed_aps: a (frame, AP) pair whose Stage B measurement failed
 * is dropped from the stack (with its per-AP weight reduced to match), so a
 * corrupted shift can no longer smear the patch.  An AP whose pairs ALL
 * failed keeps stacking everything (fallback: misaligned beats hole). */
Test(mpp_stack, skip_failed_aps_drops_corrupted_pairs) {
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
	q[0] = 1.0;

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
	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	mpp_shifts_t *shifts = mpp::shift_compute_all(frames_blurred, *aps,
	                                              ref_boxes, offsets, cfg);
	cr_assert_not_null(shifts);

	std::vector<int> sorted(frames_raw.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(),
	          [&q](int a, int b) { return q[a] > q[b]; });

	const int M = aps->count;
	/* The merged image is trimmed by the per-run RemapBorder, which a
	 * corrupted shift can change — so compare on a fixed canvas window,
	 * mapped into each run's image through its own border origin. */
	auto run_merge = [&](const mpp_config_t &c) {
		const auto out = mpp::stack_apply_shifts(frames_raw, *aps, apq, shifts,
		                                         offsets, brightness, sorted,
		                                         avg.intersection, c, nullptr,
		                                         /*max_threads=*/1);
		const cv::Mat img = mpp::stack_merge_alignment_point_buffers(
		    out.state, out.border, *aps, c);
		const int m = 24;   /* canvas margin large enough to cover any border */
		const int y0 = m - out.border.y_low, x0 = m - out.border.x_low;
		const int h = out.state.dim_y_drizzled - 2 * m;
		const int w = out.state.dim_x_drizzled - 2 * m;
		return img(cv::Range(y0, y0 + h), cv::Range(x0, x0 + w)).clone();
	};

	/* Clean baseline (all measurements succeeded). */
	const cv::Mat clean = run_merge(cfg);

	/* Sabotage one pair: frame 3 at AP 0 gets a wildly wrong shift and a
	 * failed flag — exactly what a Stage B failure looks like since
	 * failures now retain a (possibly bad) coarse estimate. */
	const int bad_f = 3, bad_a = 0;
	shifts->shifts[(size_t)(bad_f * M + bad_a) * 2 + 0] = 9.0;
	shifts->shifts[(size_t)(bad_f * M + bad_a) * 2 + 1] = -9.0;
	shifts->success[(size_t)bad_f * M + bad_a] = 0;

	cv::Mat diff;
	double dmin = 0.0, err_noskip = 0.0, err_skip = 0.0;

	const cv::Mat with_corruption = run_merge(cfg);
	cv::absdiff(with_corruption, clean, diff);
	cv::minMaxLoc(diff.reshape(1), &dmin, &err_noskip);

	auto cfg_skip = cfg;
	cfg_skip.stack_skip_failed_aps = true;
	const cv::Mat skipped = run_merge(cfg_skip);
	cv::absdiff(skipped, clean, diff);
	cv::minMaxLoc(diff.reshape(1), &dmin, &err_skip);

	/* The corrupted pair must visibly damage the no-skip stack, and the
	 * skip option must essentially restore the clean result (dropping one
	 * of six near-identical frames at one AP only moves the local average
	 * by interpolation noise). */
	cr_assert_gt(err_noskip, 100.0,
	             "sabotaged pair should disturb the stack (err=%.0f)", err_noskip);
	cr_assert_lt(err_skip, err_noskip / 4.0,
	             "skip should beat no-skip (skip=%.0f, noskip=%.0f)",
	             err_skip, err_noskip);

	/* Saturated-AP fallback: when every pair of an AP failed, skipping is
	 * disabled for that AP and the result matches the no-skip run. */
	for (int f = 0; f < (int) frames_raw.size(); ++f)
		shifts->success[(size_t)f * M + bad_a] = 0;
	const cv::Mat sat_skip   = run_merge(cfg_skip);
	const cv::Mat sat_noskip = run_merge(cfg);
	cv::absdiff(sat_skip, sat_noskip, diff);
	double sat_max = 0.0;
	cv::minMaxLoc(diff.reshape(1), &dmin, &sat_max);
	cr_assert_leq(sat_max, 0.0,
	              "saturated AP must stack unfiltered (diff=%.0f)", sat_max);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}
