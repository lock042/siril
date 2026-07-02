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
	const double s = std::sin(G_PI_2 * t);
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

	auto stack_with = [&](const mpp::APQualities &q_, const mpp_config_t &c_) {
		return mpp::stack_warp_apply_streamed(
		    [&frames_raw](int i) { return frames_raw[i]; },
		    (int) frames_raw.size(), 1, *aps, q_, shifts, offsets,
		    brightness, sorted, avg.intersection, c_, nullptr,
		    1, false, 0, nullptr);
	};
	const cv::Mat ms = stack_with(apq, cfg).image;
	const cv::Mat mf = stack_with(apq_full, cfg_full).image;
	cr_assert(!ms.empty());
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

/* Warp-field engine: it must (a) produce the scene from a set of identical
 * frames, and (b) apply a supplied derotation base map. With identical frames
 * and zero global/AP shifts the engine reduces to remap(frame, base_map), so a
 * pure -3px x translation in the derot provider shifts the output by exactly 3
 * columns — the load-bearing property for folding derotation into the stack's
 * single resample. mu = 1 everywhere, so weighting is unchanged. */
Test(mpp_stack, warp_engine_applies_derot_base) {
	const auto cfg = default_cfg();
	const cv::Mat truth = make_textured_frame();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q;
	for (int i = 0; i < 4; ++i) {
		frames_raw.push_back(truth.clone());           /* identical frames */
		frames_blurred.push_back(blurred(truth, cfg));
		q.push_back(i == 0 ? 1.0 : 0.5);
	}
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
	const cv::Vec4i isect = avg.intersection;
	mpp::FrameProvider provider = [&](int i) { return frames_raw[i]; };

	auto run = [&](const mpp::DerotMapProvider *derot) {
		return mpp::stack_warp_apply_streamed(
		    provider, (int) frames_raw.size(), 1, *aps, apq, shifts, offsets,
		    brightness, sorted, isect, cfg, nullptr, 1, false, 0, derot);
	};

	const auto ident = run(nullptr);
	cr_assert(!ident.image.empty());
	cr_assert(!ident.oom && !ident.cancelled);
	const int DY = isect[1] - isect[0], DX = isect[3] - isect[2];
	cr_assert_eq(ident.image.rows, DY);
	cr_assert_eq(ident.image.cols, DX);
	/* the engine reproduced the disk: centre much brighter than a corner */
	cr_assert_gt(ident.image.at<uint16_t>(DY / 2, DX / 2),
	             ident.image.at<uint16_t>(4, 4) + 5000);

	const int tx = 3;
	mpp::DerotMapProvider shift_provider =
	    [tx](int, int dy, int dx, cv::Mat &mx, cv::Mat &my, cv::Mat &mu) {
		mx.create(dy, dx, CV_32F); my.create(dy, dx, CV_32F); mu.create(dy, dx, CV_32F);
		for (int y = 0; y < dy; ++y)
			for (int x = 0; x < dx; ++x) {
				mx.at<float>(y, x) = (float) (x - tx);
				my.at<float>(y, x) = (float) y;
				mu.at<float>(y, x) = 1.0f;
			}
		return true;
	};
	const auto sh = run(&shift_provider);
	cr_assert(!sh.image.empty());

	/* interior: shifted output column x equals identity column x - tx */
	double max_d = 0.0;
	for (int y = 30; y < DY - 30; ++y)
		for (int x = 30; x < DX - 30; ++x) {
			const double d = std::abs((double) sh.image.at<uint16_t>(y, x)
			                          - (double) ident.image.at<uint16_t>(y, x - tx));
			max_d = std::max(max_d, d);
		}
	cr_assert_leq(max_d, 3.0,
	              "derot base map not applied: interior max diff %.0f LSB", max_d);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);
}
