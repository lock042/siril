/*
 * Private C++ entry points for mpp_shift. Tests drive these directly using
 * cv::Mat inputs (no Siril sequence required).
 */
#ifndef SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_
#define SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_shift.h"

namespace mpp {

struct APRefBoxes {
	cv::Mat second_phase;   /* CV_32F, full-resolution patch from blurred mean_frame */
	cv::Mat first_phase;    /* CV_32F, stride-2 view of second_phase */
};

/* Build per-AP reference boxes. Reference boxes for correlation are
 * built from the *unblurred* mean frame — NOT the blurred mean frame
 * used for AP placement. Pass the raw output of mpp::align_average_frame,
 * not blur_mean_frame_for_ap's output. */
std::vector<APRefBoxes> shift_prepare_ref_boxes(const cv::Mat &mean_frame_raw,
                                                const mpp_aps_t &aps);

/* Per-frame offsets:
 *   dy = intersection_y_low - global_shifts[idx].dy
 * Used to translate AP box bounds (mean-frame coords) into each frame's
 * native coordinates. Sub-pixel since the global align pass switched to
 * parabolic refinement; consumers that need integer indices (Stage B
 * box positioning) round explicitly, the output resample uses the
 * doubles directly. */
struct FrameOffset { double dy; double dx; };
std::vector<FrameOffset> shift_frame_offsets(const std::vector<cv::Vec2d> &global_shifts,
                                             const cv::Vec4i &intersection);

/* Compute per-AP per-frame shifts. Returns a heap-allocated mpp_shifts_t —
 * caller frees via mpp_shift_free. */
mpp_shifts_t *shift_compute_all(const std::vector<cv::Mat> &frames_mono_blurred,
                                const mpp_aps_t &aps,
                                const std::vector<APRefBoxes> &ref_boxes,
                                const std::vector<FrameOffset> &offsets,
                                const mpp_config_t &cfg);

/* Per-frame robust smoothing of the shift field across APs (mpp_improve;
 * no PSS equivalent — cfg.alignment_points_smooth_radius, 0 disables).
 * For each frame, each AP's (dy, dx) is replaced by a robust local plane
 * fit (tricube distance weights × Tukey bisquare residual reweighting,
 * 3 IRLS rounds) over the successful measurements of the APs within
 * radius × grid-step of it. Rationale: the true warp field is smooth at
 * grid-step scale while the per-(frame, AP) measurement noise is
 * independent, so the fit averages the noise down by ~sqrt(n/3) without
 * biasing genuine warp; the bisquare loop rejects railed/aperture-
 * problem outliers. Failed pairs inside a neighbourhood with at least 5
 * successful measurements receive the fit prediction (success flags are
 * left untouched — the skip-failed stack option keeps its meaning).
 * APs with fewer successful neighbours keep their raw values. Excluded
 * frames (included[f] == 0) are skipped. In-place; thread-safe across
 * frames (max_threads > 1 parallelises the frame loop). */
void shift_field_smooth(mpp_shifts_t *shifts, const mpp_aps_t &aps,
                        const mpp_config_t &cfg, const int *included,
                        int max_threads = 1);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_SHIFT_PRIV_HPP_ */
