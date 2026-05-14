/*
 * Phase 3: alignment-point grid placement.
 *
 * Ports PSS alignment_points.create_ap_grid (alignment_points.py:135) plus
 * its helpers ap_locations and new_alignment_point, and
 * Miscellaneous.quality_measure (miscellaneous.py:50).
 *
 * Pipeline:
 *   1. Compute staggered y/x grid locations via ap_locations.
 *   2. For each (y, x), build a reference box of half-box width centred on
 *      it and a patch (1.5× larger) that extends to the frame border for
 *      edge APs.
 *   3. Filter by brightness (max(box) > 10×256) and contrast (max−min > 0).
 *      Drop the AP if the box fails.
 *   4. If too many pixels are dim (fraction < bright > dim_fraction_threshold),
 *      recompute the AP centre as the brightness centre-of-mass within the
 *      box and rebuild the AP.
 *   5. Compute the per-AP structure score = min(avg|∂x|, avg|∂y|) on the
 *      reference box.
 *   6. After all APs are placed: normalise structure by the max across the
 *      kept APs, then drop APs whose normalised structure < 0.04.
 *
 * Coordinates use the mean (reference) frame's pixel space. Sign convention
 * matches Phase 2.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_ap_priv.hpp"

namespace mpp {

std::vector<int> ap_locations(int num_pixels, int min_boundary_distance,
                              int step_size, bool even) {
	std::vector<int> locations;
	const int span = num_pixels - 2 * min_boundary_distance;
	if (span <= 0 || step_size <= 0)
		return locations;
	const int num_interior_odd = (int) std::ceil((double) span / (double) step_size);
	if (num_interior_odd <= 0)
		return locations;
	const int num_interior_even = num_interior_odd + 1;
	const double distance_corrected = (double) span / (double) num_interior_odd;
	const int count = even ? num_interior_even : num_interior_odd;
	const double offset = even ? 0.0 : 0.5;
	for (int i = 0; i < count; ++i)
		locations.push_back((int) (min_boundary_distance + (offset + i) * distance_corrected));
	return locations;
}

double ap_quality_measure(const cv::Mat &box) {
	/* PSS: dx = box[:, 1:] - box[:, :-1]; dy = box[1:, :] - box[:-1, :];
	 *      return min(mean(|dx|), mean(|dy|))
	 * Input is the mean-frame patch, typically CV_32S; PSS treats negative
	 * differences correctly because the dtype is signed. */
	double sum_dx = 0.0, sum_dy = 0.0;
	const int H = box.rows, W = box.cols;
	if (H < 2 || W < 2) return 0.0;

	auto loop = [&](auto pix) {
		for (int y = 0; y < H; ++y)
			for (int x = 0; x < W - 1; ++x)
				sum_dx += std::abs((double) (pix(y, x + 1) - pix(y, x)));
		for (int y = 0; y < H - 1; ++y)
			for (int x = 0; x < W; ++x)
				sum_dy += std::abs((double) (pix(y + 1, x) - pix(y, x)));
	};
	if (box.type() == CV_32S) {
		loop([&](int y, int x) { return box.at<int32_t>(y, x); });
	} else if (box.type() == CV_32F) {
		loop([&](int y, int x) { return box.at<float>(y, x); });
	} else if (box.type() == CV_16U) {
		loop([&](int y, int x) { return (int) box.at<uint16_t>(y, x); });
	} else if (box.type() == CV_8U) {
		loop([&](int y, int x) { return (int) box.at<uint8_t>(y, x); });
	} else {
		return 0.0;  /* unsupported */
	}
	const double mean_dx = sum_dx / (double) (H * (W - 1));
	const double mean_dy = sum_dy / (double) ((H - 1) * W);
	return std::min(mean_dx, mean_dy);
}

std::pair<int, int> ap_center_of_mass(const cv::Mat &box) {
	/* ndimage.measurements.center_of_mass: weighted by pixel value.
	 *   y_com = sum(y * box) / sum(box); x_com analogously. */
	const int H = box.rows, W = box.cols;
	double total = 0.0, weighted_y = 0.0, weighted_x = 0.0;
	auto loop = [&](auto pix) {
		for (int y = 0; y < H; ++y)
			for (int x = 0; x < W; ++x) {
				const double v = (double) pix(y, x);
				total += v;
				weighted_y += y * v;
				weighted_x += x * v;
			}
	};
	if (box.type() == CV_32S)      loop([&](int y, int x) { return box.at<int32_t>(y, x); });
	else if (box.type() == CV_32F) loop([&](int y, int x) { return box.at<float>(y, x); });
	else if (box.type() == CV_16U) loop([&](int y, int x) { return (int) box.at<uint16_t>(y, x); });
	else if (box.type() == CV_8U)  loop([&](int y, int x) { return (int) box.at<uint8_t>(y, x); });
	if (total == 0.0) return {H / 2, W / 2};  /* defensive */
	return { (int) (weighted_y / total), (int) (weighted_x / total) };
}

namespace {

/* Populate a single AP record. Equivalent to PSS's new_alignment_point but
 * returns success/failure inline. */
bool fill_ap_record(mpp_ap_record_t &ap, int y, int x,
                    bool extend_x_low, bool extend_x_high,
                    bool extend_y_low, bool extend_y_high,
                    int num_pixels_y, int num_pixels_x,
                    const mpp_config_t &cfg) {
	const int hb = cfg.alignment_points_half_box_width;
	const int hp = mpp_cfg_half_patch_width(&cfg);
	const int sw = cfg.alignment_points_search_width;
	const int min_dist = std::max(hb + sw, hp);
	if (y < min_dist || y > num_pixels_y - min_dist
	 || x < min_dist || x > num_pixels_x - min_dist)
		return false;

	ap.y = y;
	ap.x = x;
	ap.box_y_low  = y - hb;
	ap.box_y_high = y + hb;
	ap.box_x_low  = x - hb;
	ap.box_x_high = x + hb;

	ap.patch_y_low  = extend_y_low  ? 0             : (y - hp);
	ap.patch_y_high = extend_y_high ? num_pixels_y  : (y + hp);
	ap.patch_x_low  = extend_x_low  ? 0             : (x - hp);
	ap.patch_x_high = extend_x_high ? num_pixels_x  : (x + hp);

	ap.structure = 0.0;
	return true;
}

}  // namespace

cv::Mat blur_mean_frame_for_ap(const cv::Mat &mean_frame_raw, const mpp_config_t &cfg) {
	/* PSS AlignmentPoints.__init__ (alignment_points.py:76):
	 *   mean_frame = GaussianBlur(align.mean_frame.astype(uint16),
	 *                             (frames_gauss_width,)*2, 0).astype(int32) */
	cv::Mat u16;
	mean_frame_raw.convertTo(u16, CV_16U);
	cv::Mat blurred;
	cv::GaussianBlur(u16, blurred,
	                 cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	cv::Mat out;
	blurred.convertTo(out, CV_32S);
	return out;
}

mpp_aps_t *ap_create_grid(const cv::Mat &mean_frame, const mpp_config_t &cfg) {
	const int H = mean_frame.rows, W = mean_frame.cols;
	const int hb = cfg.alignment_points_half_box_width;
	const int hp = mpp_cfg_half_patch_width(&cfg);
	const int sw = cfg.alignment_points_search_width;
	const int step = mpp_cfg_step_size(&cfg);
	const int min_dist = std::max(hb + sw, hp);
	const double brightness_thresh = (double) cfg.alignment_points_brightness_threshold * 256.0;
	const double contrast_thresh   = (double) cfg.alignment_points_contrast_threshold * 256.0;

	const auto ys     = ap_locations(H, min_dist, step, true);
	const auto xs_e   = ap_locations(W, min_dist, step, true);
	const auto xs_o   = ap_locations(W, min_dist, step, false);

	mpp_aps_t *out = (mpp_aps_t *) std::calloc(1, sizeof(*out));
	if (!out) return nullptr;
	if (ys.empty() || xs_e.empty() || xs_o.empty()) return out;

	std::vector<mpp_ap_record_t> kept;
	kept.reserve(ys.size() * xs_e.size());
	int dropped_dim = 0;
	bool even_row = true;

	for (size_t iy = 0; iy < ys.size(); ++iy) {
		const int y = ys[iy];
		const bool extend_y_low  = (iy == 0);
		const bool extend_y_high = (iy == ys.size() - 1);
		const auto &xs = even_row ? xs_e : xs_o;
		even_row = !even_row;

		for (size_t ix = 0; ix < xs.size(); ++ix) {
			const int x = xs[ix];
			const bool extend_x_low  = (ix == 0);
			const bool extend_x_high = (ix == xs.size() - 1);

			mpp_ap_record_t ap;
			if (!fill_ap_record(ap, y, x, extend_x_low, extend_x_high,
			                    extend_y_low, extend_y_high, H, W, cfg))
				continue;

			const cv::Mat box = mean_frame(cv::Range(ap.box_y_low, ap.box_y_high),
			                               cv::Range(ap.box_x_low, ap.box_x_high));

			double max_b = 0, min_b = 0;
			cv::minMaxLoc(box, &min_b, &max_b);

			if (!(max_b > brightness_thresh && (max_b - min_b) > contrast_thresh)) {
				++dropped_dim;
				continue;
			}

			/* PSS condition is `fraction > threshold` where fraction is the
			 * proportion of pixels BELOW the brightness threshold. */
			cv::Mat dim_mask;
			cv::compare(box, brightness_thresh, dim_mask, cv::CMP_LT);
			const double dim_fraction = (double) cv::countNonZero(dim_mask)
			                          / (double) (box.rows * box.cols);
			if (dim_fraction > cfg.alignment_points_dim_fraction_threshold) {
				const auto com = ap_center_of_mass(box);
				int y_adapt = y + com.first  - hb;
				int x_adapt = x + com.second - hb;
				y_adapt = std::max(std::min(y_adapt, H - min_dist), min_dist);
				x_adapt = std::max(std::min(x_adapt, W - min_dist), min_dist);
				if (!fill_ap_record(ap, y_adapt, x_adapt,
				                    extend_x_low, extend_x_high,
				                    extend_y_low, extend_y_high, H, W, cfg)) {
					++dropped_dim;  /* still no valid placement */
					continue;
				}
			}

			const cv::Mat box_final = mean_frame(cv::Range(ap.box_y_low, ap.box_y_high),
			                                     cv::Range(ap.box_x_low, ap.box_x_high));
			ap.structure = ap_quality_measure(box_final);
			kept.push_back(ap);
		}
	}

	/* Normalize structure by max and drop weak APs. */
	int dropped_structure = 0;
	if (!kept.empty()) {
		double smax = 0.0;
		for (const auto &ap : kept) smax = std::max(smax, ap.structure);
		if (smax > 0.0) {
			std::vector<mpp_ap_record_t> filtered;
			filtered.reserve(kept.size());
			for (auto ap : kept) {
				ap.structure /= smax;
				if (ap.structure < cfg.alignment_points_structure_threshold) {
					++dropped_structure;
					continue;
				}
				filtered.push_back(ap);
			}
			kept = std::move(filtered);
		}
	}

	out->count = (int) kept.size();
	out->dropped_dim = dropped_dim;
	out->dropped_structure = dropped_structure;
	if (out->count > 0) {
		out->records = (mpp_ap_record_t *) std::malloc(sizeof(*out->records) * out->count);
		if (!out->records) {
			std::free(out);
			return nullptr;
		}
		std::memcpy(out->records, kept.data(), sizeof(*out->records) * out->count);
	}
	return out;
}

}  // namespace mpp

/* ----------------- Public C interface ----------------- */

extern "C" mpp_status_t mpp_ap_place(const fits *ref,
                                     const mpp_config_t *cfg,
                                     mpp_aps_t **aps_out) {
	(void) ref; (void) cfg;
	/* Pending mpp_frames integration: tests drive ap_create_grid directly. */
	if (aps_out) *aps_out = nullptr;
	return MPP_ENOTIMPL;
}

extern "C" void mpp_ap_free(mpp_aps_t *aps) {
	if (!aps) return;
	if (aps->records) std::free(aps->records);
	std::free(aps);
}
