/*
 * Phase 3 — mpp_ap Criterion tests.
 *
 * Covers PSS-faithful staggered AP grid placement + filtering via the
 * private C++ entry points, plus the C-level mpp_ap_place stub.
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

#include "registration/mpp_ap_priv.hpp"
#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_rank_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
}

namespace {

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

cv::Mat blurred(const cv::Mat &mono, const mpp_config_t &cfg) {
	cv::Mat b;
	cv::GaussianBlur(mono, b, cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	return b;
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

Test(mpp_ap, stub_place_returns_enotimpl) {
	mpp_aps_t *aps = nullptr;
	cr_assert_eq(mpp_ap_place(nullptr, nullptr, &aps), MPP_ENOTIMPL);
	cr_assert_null(aps);
	mpp_ap_free(aps);
}

Test(mpp_ap, default_config_has_pss_values) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.alignment_points_half_box_width, 24);
	cr_assert_eq(cfg.alignment_points_search_width, 14);
	cr_assert_eq(cfg.alignment_points_brightness_threshold, 10);
	cr_assert_eq(cfg.alignment_points_contrast_threshold, 0);
	cr_assert_float_eq(cfg.alignment_points_structure_threshold, 0.04, 1e-12);
	cr_assert_float_eq(cfg.alignment_points_dim_fraction_threshold, 0.6, 1e-12);
	/* Derived helpers reproduce PSS's values. */
	cr_assert_eq(mpp_cfg_half_patch_width(&cfg), 36);
	cr_assert_eq(mpp_cfg_step_size(&cfg), 54);
}

/* PSS ap_locations: distance_corrected = (num_pixels - 2*min_dist) / N_odd
 *   where N_odd = ceil((num_pixels - 2*min_dist) / step_size).
 * For 200 px frame, min_dist=38, step=54: span=124, N_odd=ceil(124/54)=3,
 * distance_corrected = 124/3 ≈ 41.333.
 * Even locations: [38 + 0, 38 + 41.33, 38 + 82.67, 38 + 124] = [38, 79, 120, 162]
 * Odd locations:  [38 + 20.67, 38 + 62, 38 + 103.33]      = [58, 100, 141] */
Test(mpp_ap, ap_locations_matches_pss_formula) {
	const auto even = mpp::ap_locations(200, 38, 54, true);
	cr_assert_eq(even.size(), 4u);
	cr_assert_eq(even[0], 38);
	cr_assert_eq(even[1], 79);
	cr_assert_eq(even[2], 120);
	cr_assert_eq(even[3], 162);

	const auto odd = mpp::ap_locations(200, 38, 54, false);
	cr_assert_eq(odd.size(), 3u);
	cr_assert_eq(odd[0], 58);
	cr_assert_eq(odd[1], 100);
	cr_assert_eq(odd[2], 141);
}

Test(mpp_ap, ap_locations_empty_when_frame_too_small) {
	cr_assert(mpp::ap_locations(50, 38, 54, true).empty());
}

/* Quality measure: min(avg|∂x|, avg|∂y|). A constant box has 0 quality. */
Test(mpp_ap, quality_measure_zero_on_constant_box) {
	cv::Mat box(48, 48, CV_32S, cv::Scalar(1000));
	cr_assert_float_eq(mpp::ap_quality_measure(box), 0.0, 1e-12);
}

/* A monotonically increasing ramp has equal ∂x and ∂y if symmetric. */
Test(mpp_ap, quality_measure_finds_structure_in_ramp) {
	cv::Mat box(48, 48, CV_32S);
	for (int y = 0; y < 48; ++y)
		for (int x = 0; x < 48; ++x)
			box.at<int32_t>(y, x) = 100 * (x + y);  /* shared gradient */
	const double q = mpp::ap_quality_measure(box);
	cr_assert_float_eq(q, 100.0, 1e-6,
	                   "avg |∂x| = avg |∂y| = 100 expected, got %.6f", q);
}

/* COM of a single bright pixel at (3, 5) in a 16x16 box should round to (3, 5). */
Test(mpp_ap, center_of_mass_single_bright_pixel) {
	cv::Mat box(16, 16, CV_32S, cv::Scalar(0));
	box.at<int32_t>(3, 5) = 1000;
	const auto com = mpp::ap_center_of_mass(box);
	cr_assert_eq(com.first, 3);
	cr_assert_eq(com.second, 5);
}

/* Build a grid on a sufficiently bright frame and verify the placement has
 * reasonable invariants. */
Test(mpp_ap, create_grid_basic_invariants) {
	const auto cfg = default_cfg();
	cv::Mat ref(200, 240, CV_32S, cv::Scalar(0));
	/* Put a bright textured disk in the centre. */
	for (int y = 0; y < ref.rows; ++y)
		for (int x = 0; x < ref.cols; ++x) {
			const double dy = y - ref.rows / 2.0, dx = x - ref.cols / 2.0;
			if (std::hypot(dy, dx) < 50) ref.at<int32_t>(y, x) = 8000 + (int) (200 * std::sin(0.4 * x + 0.3 * y));
		}

	mpp_aps_t *aps = mpp::ap_create_grid(ref, cfg);
	cr_assert_not_null(aps);
	for (int i = 0; i < aps->count; ++i) {
		const auto &ap = aps->records[i];
		cr_assert_lt(ap.box_y_low, ap.box_y_high);
		cr_assert_lt(ap.box_x_low, ap.box_x_high);
		cr_assert_eq(ap.box_y_high - ap.box_y_low, 2 * cfg.alignment_points_half_box_width);
		cr_assert_geq(ap.patch_y_low, 0);
		cr_assert_leq(ap.patch_y_high, ref.rows);
		cr_assert_geq(ap.patch_x_low, 0);
		cr_assert_leq(ap.patch_x_high, ref.cols);
		cr_assert_geq(ap.structure, cfg.alignment_points_structure_threshold);
		cr_assert_leq(ap.structure, 1.0 + 1e-9);
	}
	mpp_ap_free(aps);
}

/* Oracle test: AP centres and box coordinates must match PSS exactly on the
 * bundled synthetic dataset. Acceptance bar from pss_port_plan.md §3.3. */
Test(mpp_ap, oracle_equivalence_synthetic) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir  = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path ap_yx_csv   = fs::path(data_root) / "oracle_out" / "ap_yx.csv";
	const fs::path ap_box_csv  = fs::path(data_root) / "oracle_out" / "ap_box.csv";

	if (!fs::exists(frames_dir) || !fs::exists(ap_yx_csv) || !fs::exists(ap_box_csv))
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

	/* Load raw + blurred views. Raw is for the average frame; blurred is for
	 * the global aligner. */
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> quality;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		cr_assert_eq(m.depth(), CV_16U);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		quality.push_back(mpp::rank_score_normalized(m, cfg));
	}

	const auto align = mpp::align_global_from_frames(frames_blurred, quality, cfg);
	const auto avg = mpp::align_average_frame(frames_raw, quality, align.shifts, cfg);
	cr_assert_gt(avg.mean_frame.rows, 0);
	cr_assert_gt(avg.mean_frame.cols, 0);

	mpp_aps_t *aps = mpp::ap_create_grid(avg.mean_frame, cfg);
	cr_assert_not_null(aps);

	const std::vector<int> oracle_yx  = read_ints_csv(ap_yx_csv);   /* n*2 */
	const std::vector<int> oracle_box = read_ints_csv(ap_box_csv);  /* n*4 */
	cr_assert_eq(oracle_yx.size() % 2, 0u);
	cr_assert_eq(oracle_box.size() % 4, 0u);
	const int oracle_count = (int) (oracle_yx.size() / 2);

	cr_assert_eq(aps->count, oracle_count,
	             "AP count mismatch: ours=%d oracle=%d", aps->count, oracle_count);

	/* PSS iterates rows-then-columns in the same order we do; compare by
	 * index. */
	for (int i = 0; i < aps->count; ++i) {
		const auto &ap = aps->records[i];
		cr_assert_eq(ap.y, oracle_yx[2 * i],
		             "AP %d y mismatch: ours=%d oracle=%d", i, ap.y, oracle_yx[2 * i]);
		cr_assert_eq(ap.x, oracle_yx[2 * i + 1],
		             "AP %d x mismatch: ours=%d oracle=%d", i, ap.x, oracle_yx[2 * i + 1]);
		cr_assert_eq(ap.box_y_low,  oracle_box[4 * i + 0], "AP %d box_y_low",  i);
		cr_assert_eq(ap.box_y_high, oracle_box[4 * i + 1], "AP %d box_y_high", i);
		cr_assert_eq(ap.box_x_low,  oracle_box[4 * i + 2], "AP %d box_x_low",  i);
		cr_assert_eq(ap.box_x_high, oracle_box[4 * i + 3], "AP %d box_x_high", i);
	}
	mpp_ap_free(aps);
}
