/*
 * Private C++ entry point for mpp_ap. Tests drive this directly with a
 * cv::Mat reference (mean) frame.
 */
#ifndef SRC_REGISTRATION_MPP_AP_PRIV_HPP_
#define SRC_REGISTRATION_MPP_AP_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"

namespace mpp {

/* Compute 1-D AP locations along one axis for either even or odd staggered
 * rows. Mirrors PSS alignment_points.ap_locations exactly. */
std::vector<int> ap_locations(int num_pixels, int min_boundary_distance,
                              int step_size, bool even);

/* PSS Miscellaneous.quality_measure: average |∂x|, average |∂y|, return the
 * smaller. Input may be CV_32S (typical: mean_frame from Phase 2) or CV_32F. */
double ap_quality_measure(const cv::Mat &box);

/* Centre-of-mass (y, x) of the brightness distribution inside `box`. Returns
 * floored integer coordinates relative to box origin, matching
 * `int(ndimage.measurements.center_of_mass(box)[i])`. */
std::pair<int, int> ap_center_of_mass(const cv::Mat &box);

/* PSS AlignmentPoints.__init__ re-blurs the mean frame from align_frames
 * before any AP placement (alignment_points.py:76):
 *   mean_frame = GaussianBlur(align.mean_frame.astype(uint16),
 *                             (frames_gauss_width,)*2, 0).astype(int32)
 * This helper applies that step. Both ap_create_grid and per-AP shift
 * computation consume the same blurred output. */
cv::Mat blur_mean_frame_for_ap(const cv::Mat &mean_frame_raw, const mpp_config_t &cfg);

/* Build the staggered AP grid on a *blurred* mean_frame and apply the PSS
 * filtering chain. Pass the output of blur_mean_frame_for_ap. */
mpp_aps_t *ap_create_grid(const cv::Mat &mean_frame_blurred, const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_AP_PRIV_HPP_ */
