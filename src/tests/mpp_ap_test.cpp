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

	const cv::Mat ref_blurred = mpp::blur_mean_frame_for_ap(ref, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(ref_blurred, cfg);
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

