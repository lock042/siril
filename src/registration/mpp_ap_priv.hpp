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

/* Build the staggered AP grid on the reference frame `mean_frame` and apply
 * the PSS filtering chain. Output indices_used / structure are populated. */
mpp_aps_t *ap_create_grid(const cv::Mat &mean_frame, const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_AP_PRIV_HPP_ */
