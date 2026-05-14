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

std::vector<int> read_ints_csv(const std::filesystem::path &path) {
	std::vector<int> v;
	std::ifstream fin(path);
	int x;
	while (fin >> x) v.push_back(x);
	return v;
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

/* Oracle test: per-AP per-frame integer shifts match PSS exactly on the
 * bundled synthetic dataset. Acceptance bar from pss_port_plan.md §4.5. */
Test(mpp_shift, oracle_equivalence_synthetic) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir       = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path ap_shifts_csv    = fs::path(data_root) / "oracle_out" / "ap_shifts.csv";
	const fs::path ap_yx_csv        = fs::path(data_root) / "oracle_out" / "ap_yx.csv";

	if (!fs::exists(frames_dir) || !fs::exists(ap_shifts_csv) || !fs::exists(ap_yx_csv))
		cr_skip_test("oracle artifacts not present — regenerate via run_pss.py");

	std::vector<fs::path> paths;
	for (const auto &entry : fs::directory_iterator(frames_dir)) {
		const std::string n = entry.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && entry.path().extension() == ".png")
			paths.push_back(entry.path());
	}
	std::sort(paths.begin(), paths.end());
	cr_assert_gt(paths.size(), 0u);

	const auto cfg = default_cfg();

	/* Replicate the full Phase 1..3 chain to get the same APs and shifts as PSS. */
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		cr_assert_eq(m.depth(), CV_16U);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		q.push_back(mpp::rank_score_normalized(m, cfg));
	}
	const auto align = mpp::align_global_from_frames(frames_blurred, q, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q, align.shifts, cfg);
	/* mpp::align_average_frame is verified bit-equivalent to PSS's mean_frame
	 * by mpp_align::average_frame_oracle_equivalence — we use its output
	 * directly here. PSS's reference boxes for per-AP shift are sliced from
	 * the *unblurred* mean (set_reference_boxes_correlation, alignment_points.py:506),
	 * while AP placement uses the blurred mean. */
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);

	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	const auto offsets   = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	mpp_shifts_t *out = mpp::shift_compute_all(frames_blurred, *aps,
	                                           ref_boxes, offsets, cfg);
	cr_assert_not_null(out);

	/* Oracle ap_shifts.csv is shape (num_frames, 2 * num_aps), each row
	 * holding (y0, x0, y1, x1, ..., y_{m-1}, x_{m-1}). */
	const std::vector<int> oracle_flat = read_ints_csv(ap_shifts_csv);
	cr_assert_eq(oracle_flat.size(), (size_t) (out->num_frames * out->num_aps * 2),
	             "shift count mismatch: ours=%d oracle=%zu",
	             out->num_frames * out->num_aps * 2, oracle_flat.size());

	/* Compare per-frame, per-AP. PSS returns integer shifts in default mode;
	 * we should match exactly. */
	int worst_diff = 0;
	int worst_f = -1, worst_a = -1;
	for (int f = 0; f < out->num_frames; ++f) {
		for (int a = 0; a < out->num_aps; ++a) {
			const size_t off = (size_t) (f * out->num_aps + a) * 2;
			const int ours_y  = (int) std::lround(out->shifts[off + 0]);
			const int ours_x  = (int) std::lround(out->shifts[off + 1]);
			const int oracle_y = oracle_flat[off + 0];
			const int oracle_x = oracle_flat[off + 1];
			const int d = std::max(std::abs(ours_y - oracle_y),
			                       std::abs(ours_x - oracle_x));
			if (d > worst_diff) {
				worst_diff = d; worst_f = f; worst_a = a;
			}
		}
	}
	cr_assert_eq(worst_diff, 0,
	             "shift mismatch at frame %d AP %d: ours=(%.4f, %.4f) oracle=(%d, %d)",
	             worst_f, worst_a,
	             out->shifts[(size_t)(worst_f * out->num_aps + worst_a) * 2 + 0],
	             out->shifts[(size_t)(worst_f * out->num_aps + worst_a) * 2 + 1],
	             oracle_flat[(size_t)(worst_f * out->num_aps + worst_a) * 2 + 0],
	             oracle_flat[(size_t)(worst_f * out->num_aps + worst_a) * 2 + 1]);

	mpp_shift_free(out);
	mpp_ap_free(aps);
}
