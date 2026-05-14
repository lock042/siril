/*
 * Private C++ entry points for mpp_shift. Tests drive these directly using
 * cv::Mat inputs (no Siril sequence required).
 */
#ifndef SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_
#define SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_shift.h"

namespace mpp {

struct APRefBoxes {
	cv::Mat second_phase;   /* CV_32F, full-resolution patch from blurred mean_frame */
	cv::Mat first_phase;    /* CV_32F, stride-2 view of second_phase */
};

/* Build per-AP reference boxes. PSS's set_reference_boxes_correlation
 * (alignment_points.py:506) reads from the *unblurred* `align_frames.mean_frame`
 * — NOT the blurred mean frame used for AP placement. Pass the raw output of
 * mpp::align_average_frame, not blur_mean_frame_for_ap's output. */
std::vector<APRefBoxes> shift_prepare_ref_boxes(const cv::Mat &mean_frame_raw,
                                                const mpp_aps_t &aps);

/* Per-frame offsets from PSS align_frames.dy / dx:
 *   dy = intersection_y_low - global_shifts[idx].dy
 * Used to translate AP box bounds (mean-frame coords) into each frame's
 * native coordinates. */
struct FrameOffset { int dy; int dx; };
std::vector<FrameOffset> shift_frame_offsets(const std::vector<cv::Vec2i> &global_shifts,
                                             const cv::Vec4i &intersection);

/* Compute per-AP per-frame shifts. Returns a heap-allocated mpp_shifts_t —
 * caller frees via mpp_shift_free. */
mpp_shifts_t *shift_compute_all(const std::vector<cv::Mat> &frames_mono_blurred,
                                const mpp_aps_t &aps,
                                const std::vector<APRefBoxes> &ref_boxes,
                                const std::vector<FrameOffset> &offsets,
                                const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_ */
