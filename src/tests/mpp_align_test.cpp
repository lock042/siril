/*
 * Phase 2 — mpp_align Criterion tests.
 *
 * Covers PSS-faithful patch picker (compute_alignment_rect /
 * quality_measure_threshold_weighted) and two-phase MultiLevelCorrelation
 * via the private C++ entry points, plus the C-level mpp_align_global stub.
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_rank_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_align.h"
#include "registration/mpp_config.h"
}

namespace {

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

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
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

Test(mpp_align, stub_global_returns_enotimpl) {
	int shifts[4] = {0};
	int patch[4] = {0};
	cr_assert_eq(mpp_align_global(nullptr, nullptr, nullptr, shifts, patch, nullptr),
	             MPP_ENOTIMPL);
}

Test(mpp_align, default_config_has_pss_values) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.align_frames_search_width, 34);
	cr_assert_float_eq(cfg.align_frames_rectangle_scale_factor, 3.0, 1e-12);
	cr_assert_eq(cfg.align_frames_border_width, 10);
	cr_assert_eq(cfg.align_frames_rectangle_stride, 2);
	cr_assert_eq(cfg.align_frames_rectangle_black_threshold, 10240);
	cr_assert_float_eq(cfg.align_frames_rectangle_min_fraction, 0.7, 1e-12);
	cr_assert_eq(cfg.align_frames_average_frame_percent, 5);
}

Test(mpp_align, patch_picker_returns_bounds_within_frame) {
	const auto cfg = default_cfg();
	const cv::Mat f = blurred(make_textured_frame(), cfg);
	const cv::Vec4i p = mpp::align_pick_patch(f, cfg);
	cr_assert_lt(p[0], p[1]);
	cr_assert_lt(p[2], p[3]);
	cr_assert_geq(p[0], cfg.align_frames_border_width + cfg.align_frames_search_width);
	cr_assert_leq(p[1], f.rows - (cfg.align_frames_border_width + cfg.align_frames_search_width));
	cr_assert_geq(p[2], cfg.align_frames_border_width + cfg.align_frames_search_width);
	cr_assert_leq(p[3], f.cols - (cfg.align_frames_border_width + cfg.align_frames_search_width));
}

/* MultiLevelCorrelation must recover a known integer shift. */
Test(mpp_align, multilevel_correlation_recovers_known_shift) {
	const auto cfg = default_cfg();
	const cv::Mat truth = blurred(make_textured_frame(), cfg);

	const cv::Vec4i patch = mpp::align_pick_patch(truth, cfg);
	cv::Mat ref_f32;
	truth.convertTo(ref_f32, CV_32F);
	const cv::Mat ref_window = ref_f32(cv::Range(patch[0], patch[1]),
	                                   cv::Range(patch[2], patch[3]));
	cv::Mat ref_window_first;
	{
		const int s = 2;
		const int new_h = (ref_window.rows + s - 1) / s;
		const int new_w = (ref_window.cols + s - 1) / s;
		ref_window_first.create(new_h, new_w, CV_32F);
		for (int y = 0; y < new_h; ++y)
			for (int x = 0; x < new_w; ++x)
				ref_window_first.at<float>(y, x) = ref_window.at<float>(y * s, x * s);
	}

	for (int dy : {-3, -1, 0, 2, 5}) {
		for (int dx : {-4, -2, 0, 1, 3}) {
			const cv::Mat moved = shifted(truth, dy, dx);
			const auto r = mpp::align_shift_one_frame(
			    ref_window, ref_window_first, moved,
			    patch[0], patch[1], patch[2], patch[3], cfg);
			cr_assert(r.success, "shift=(%d,%d) should have succeeded", dy, dx);
			/* PSS sign convention: shift_y reported = shift to apply to *moved*
			 * to get back to truth.  truth(y,x) ≈ moved(y+dy, x+dx).  Aligning
			 * moved to truth means shifting moved by (-dy, -dx).  Wait — let's
			 * verify experimentally: warpAffine with translation (dx,dy) means
			 * the output pixel at (y,x) = input pixel at (y-dy, x-dx).  To
			 * realign, we need shift = (-dy, -dx). */
			cr_assert_eq(r.dy, -dy,
			             "dy mismatch at injected (%d,%d): expected %d, got %d",
			             dy, dx, -dy, r.dy);
			cr_assert_eq(r.dx, -dx,
			             "dx mismatch at injected (%d,%d): expected %d, got %d",
			             dy, dx, -dx, r.dx);
		}
	}
}

/* Multi-frame: each frame has known displacement, recovered shifts agree. */
Test(mpp_align, align_global_recovers_known_shifts_per_frame) {
	const auto cfg = default_cfg();
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -2}, {1, 2}, {-1, -1}, {0, 3}};

	std::vector<cv::Mat> frames;
	for (auto p : jit) frames.push_back(shifted(truth, p.first, p.second));

	/* Quality: synth makes them all "equal-ish"; bias so frame 0 wins (ref). */
	std::vector<double> q(frames.size(), 0.5);
	q[0] = 1.0;

	const auto r = mpp::align_global_from_frames(frames, q, cfg);
	cr_assert_eq(r.best_frame_idx, 0);
	for (size_t i = 0; i < jit.size(); ++i) {
		cr_assert_eq(r.shifts[i][0], -jit[i].first,
		             "frame %zu dy: expected %d, got %d", i, -jit[i].first, r.shifts[i][0]);
		cr_assert_eq(r.shifts[i][1], -jit[i].second,
		             "frame %zu dx: expected %d, got %d", i, -jit[i].second, r.shifts[i][1]);
	}
}

/* find_best_frames must pick the highest-quality subset constrained to a
 * sliding window of size region_size. */
Test(mpp_align, find_best_frames_sliding_window) {
	/* Global top-2 (0.95 at idx 3, 0.9 at idx 1) are 2 apart, so they don't
	 * fit in a window of size 2. Sliding 2-windows have sums:
	 *   [0,1]:0.1+0.9=1.0   [1,2]:0.9+0.4=1.3
	 *   [2,3]:0.4+0.95=1.35 [3,4]:0.95+0.1=1.05
	 * Best is [2,3] → indices {2, 3}. */
	const std::vector<double> q = {0.1, 0.9, 0.4, 0.95, 0.1};
	auto best = mpp::align_find_best_frames(q, 2, 2);
	cr_assert_eq(best.size(), 2u);
	std::sort(best.begin(), best.end());
	cr_assert(best == (std::vector<int>{2, 3}),
	          "expected {2,3}, got {%d, %d}", best[0], best[1]);

	/* Window covering whole sequence should return the global top-2 ({1,3}). */
	auto best_full = mpp::align_find_best_frames(q, 2, (int) q.size());
	std::sort(best_full.begin(), best_full.end());
	cr_assert(best_full == (std::vector<int>{1, 3}));
}

Test(mpp_align, average_frame_basic) {
	const auto cfg = default_cfg();
	const cv::Mat truth = make_textured_frame();
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -2}, {1, 2}, {-1, -1}, {0, 3}};
	std::vector<cv::Mat> frames;
	for (auto p : jit) frames.push_back(shifted(truth, p.first, p.second));
	const std::vector<double> q = {1.0, 0.7, 0.5, 0.85, 0.6, 0.55};
	std::vector<cv::Vec2i> shifts;
	for (auto p : jit) shifts.push_back(cv::Vec2i(-p.first, -p.second));

	auto cfg_no_fc = cfg;
	cfg_no_fc.align_frames_fast_changing_object = false;
	const auto r = mpp::align_average_frame(frames, q, shifts, cfg_no_fc);

	cr_assert_gt(r.mean_frame.rows, 0);
	cr_assert_gt(r.mean_frame.cols, 0);
	cr_assert_eq(r.mean_frame.type(), CV_32S);
	/* The intersection must cover the common region across all jittered
	 * frames; for {-2..3, -2..3} shifts and 200x240 frames, the intersection
	 * shape is (200 - 5) × (240 - 5) = 195 × 235. */
	cr_assert_eq(r.mean_frame.rows, 200 - 5);
	cr_assert_eq(r.mean_frame.cols, 240 - 5);
	const int n_target = std::max((int) std::ceil(6 * cfg.align_frames_average_frame_percent / 100.0), 1);
	cr_assert_eq((int) r.indices_used.size(), n_target);
}

/* Oracle test: patch_yxyx and integer shifts match PSS exactly on the bundled
 * synthetic dataset. Acceptance bar from pss_port_plan.md §2.4. */
Test(mpp_align, oracle_equivalence_synthetic) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir     = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path quality_csv    = fs::path(data_root) / "oracle_out" / "quality.csv";
	const fs::path shifts_csv     = fs::path(data_root) / "oracle_out" / "global_shifts.csv";
	const fs::path patch_yxyx_csv = fs::path(data_root) / "oracle_out" / "align_patch.csv";

	if (!fs::exists(frames_dir) || !fs::exists(quality_csv)
	    || !fs::exists(shifts_csv)) {
		cr_skip_test("oracle artifacts not present — regenerate via run_pss.py");
	}

	/* Load + sort the synthetic frames. */
	std::vector<fs::path> paths;
	for (const auto &entry : fs::directory_iterator(frames_dir)) {
		const std::string n = entry.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && entry.path().extension() == ".png")
			paths.push_back(entry.path());
	}
	std::sort(paths.begin(), paths.end());
	cr_assert_gt(paths.size(), 0u);

	const auto cfg = default_cfg();

	/* Build blurred-mono frames (PSS frames_mono_blurred). */
	std::vector<cv::Mat> frames_blurred;
	frames_blurred.reserve(paths.size());
	for (const auto &p : paths) {
		const cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		cr_assert_eq(m.depth(), CV_16U);
		frames_blurred.push_back(blurred(m, cfg));
	}

	/* Use σ (without normalization-by-max) ordering to find best frame; the
	 * argmax is what align_global uses. */
	std::vector<double> q(frames_blurred.size(), 0.0);
	for (size_t i = 0; i < paths.size(); ++i) {
		const cv::Mat m = cv::imread(paths[i].string(), cv::IMREAD_UNCHANGED);
		q[i] = mpp::rank_score_normalized(m, cfg);
	}

	const auto r = mpp::align_global_from_frames(frames_blurred, q, cfg);

	/* Compare per-frame integer shifts to oracle global_shifts.csv (which is
	 * (n, 2) — y first, x second, whitespace-separated). */
	const std::vector<int> oracle_shifts = read_ints_csv(shifts_csv);
	cr_assert_eq(oracle_shifts.size(), r.shifts.size() * 2);
	int worst_idx = -1;
	int worst_diff = 0;
	for (size_t i = 0; i < r.shifts.size(); ++i) {
		const int dy_diff = std::abs(r.shifts[i][0] - oracle_shifts[2 * i]);
		const int dx_diff = std::abs(r.shifts[i][1] - oracle_shifts[2 * i + 1]);
		const int diff = std::max(dy_diff, dx_diff);
		if (diff > worst_diff) { worst_diff = diff; worst_idx = (int) i; }
	}
	cr_assert_eq(worst_diff, 0,
	             "shift mismatch at frame %d: ours=(%d,%d) oracle=(%d,%d)",
	             worst_idx, r.shifts[worst_idx][0], r.shifts[worst_idx][1],
	             oracle_shifts[2 * worst_idx], oracle_shifts[2 * worst_idx + 1]);

	/* Patch coords: oracle.npz stored as `align_patch` (4-tuple); also dumped
	 * to align_patch.csv if available. */
	if (fs::exists(patch_yxyx_csv)) {
		const std::vector<int> ply = read_ints_csv(patch_yxyx_csv);
		cr_assert_eq(ply.size(), 4u);
		cr_assert_eq(r.patch_yxyx[0], ply[0]);
		cr_assert_eq(r.patch_yxyx[1], ply[1]);
		cr_assert_eq(r.patch_yxyx[2], ply[2]);
		cr_assert_eq(r.patch_yxyx[3], ply[3]);
	}
}
