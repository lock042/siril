/*
 * Phase 4 — mpp_shift Criterion tests.
 *
 * Per-AP per-frame shift via PSS's MultiLevelCorrelation, sharing the same
 * multilevel_correlation kernel as Phase 2's global aligner but with a
 * smaller search width (alignment_points_search_width = 14 vs 34).
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_ap_priv.hpp"
#include "registration/mpp_rank_priv.hpp"
#include "registration/mpp_shift_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_shift.h"
}

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

Test(mpp_shift, stub_returns_enotimpl) {
	mpp_shifts_t *s = nullptr;
	cr_assert_eq(mpp_shift_compute(nullptr, nullptr, nullptr, nullptr, nullptr, &s),
	             MPP_ENOTIMPL);
	cr_assert_null(s);
	mpp_shift_free(s);
}

Test(mpp_shift, default_config_phase4) {
	const auto cfg = default_cfg();
	cr_assert(!cfg.alignment_points_local_search_subpixel);
	cr_assert_eq(cfg.alignment_points_search_width, 14);
}

/* Build a small AP set on a static reference; the reference frame itself is
 * frame 0, plus we add three jittered copies. After Phase 2 global alignment,
 * per-AP shifts should all be near zero (the global aligner already corrected
 * the bulk translation). */
Test(mpp_shift, per_ap_shifts_near_zero_after_global_alignment) {
	const auto cfg = default_cfg();
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
	q[0] = 1.0;  /* frame 0 is the reference */

	const auto align = mpp::align_global_from_frames(frames_blurred, q, cfg);
	auto cfg_no_fc = cfg;
	cfg_no_fc.align_frames_fast_changing_object = false;
	const auto avg = mpp::align_average_frame(frames_raw, q, align.shifts, cfg_no_fc);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	cr_assert_gt(aps->count, 0);

	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	const auto offsets   = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	mpp_shifts_t *out = mpp::shift_compute_all(frames_blurred, *aps,
	                                           ref_boxes, offsets, cfg);
	cr_assert_not_null(out);
	cr_assert_eq(out->num_frames, (int) frames_blurred.size());
	cr_assert_eq(out->num_aps, aps->count);

	/* Best frame: every AP should resolve to exactly (0, 0). */
	for (int a = 0; a < aps->count; ++a) {
		const size_t off = (size_t) (align.best_frame_idx * aps->count + a) * 2;
		cr_assert_float_eq(out->shifts[off + 0], 0.0, 1e-9,
		                   "best-frame AP %d dy: expected 0, got %.4f",
		                   a, out->shifts[off + 0]);
		cr_assert_float_eq(out->shifts[off + 1], 0.0, 1e-9,
		                   "best-frame AP %d dx: expected 0, got %.4f",
		                   a, out->shifts[off + 1]);
	}

	/* Other frames: per-AP residuals should be 0 too within the integer
	 * resolution, because the global aligner already corrected the bulk
	 * translation and the simulated scene has no per-AP warp. */
	int worst_y = 0, worst_x = 0;
	int worst_count = 0;
	for (int f = 0; f < out->num_frames; ++f)
		for (int a = 0; a < aps->count; ++a) {
			const size_t off = (size_t) (f * aps->count + a) * 2;
			worst_y = std::max(worst_y, (int) std::abs(out->shifts[off + 0]));
			worst_x = std::max(worst_x, (int) std::abs(out->shifts[off + 1]));
			if (out->shifts[off + 0] || out->shifts[off + 1]) ++worst_count;
		}
	cr_assert_leq(worst_y, 1);
	cr_assert_leq(worst_x, 1);
	mpp_shift_free(out);
	mpp_ap_free(aps);
}

