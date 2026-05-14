/*
 * Private C++ entry points for mpp_align. Tests drive these directly with
 * cv::Mat inputs (without going through a Siril sequence).
 */
#ifndef SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_
#define SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_config.h"

namespace mpp {

/* PSS quality_measure_threshold_weighted from miscellaneous.py:88.
 * `frame` is the candidate patch (best-frame mono_blurred cropped to (H, W)).
 * Returns the quality score (higher = more structure). Replicates PSS's
 * NumPy uint16 modular subtraction for u16 input (see oracle equivalence notes
 * in pss_port_plan.md / commit log). For CV_32F input, uses signed math. */
double align_quality_measure_threshold_weighted(const cv::Mat &frame,
                                                int stride,
                                                double black_threshold,
                                                double min_fraction);

/* Picks the (y_low, y_high, x_low, x_high) patch with the highest quality
 * over the coarse grid PSS searches in `compute_alignment_rect`. The input
 * must be the best frame's `frames_mono_blurred`. */
cv::Vec4i align_pick_patch(const cv::Mat &best_frame_mono_blurred,
                           const mpp_config_t &cfg);

struct AlignShiftResult {
	int dy = 0;          /* PSS sign convention: shift the frame by (dy,dx) to align with ref */
	int dx = 0;
	bool success = false;
};

/* PSS multilevel_correlation. Reference windows are pre-prepared float32
 * boxes from the best frame's `frames_mono_blurred`:
 *   reference_window_f32              — full-resolution patch
 *   reference_window_first_phase_f32  — stride-2 view of the above
 * The input `frame_mono_blurred` is whatever dtype the frame source produces
 * (typically uint16); we'll convert internally where matchTemplate requires
 * float32. */
AlignShiftResult align_shift_one_frame(const cv::Mat &reference_window_f32,
                                       const cv::Mat &reference_window_first_phase_f32,
                                       const cv::Mat &frame_mono_blurred,
                                       int patch_y_low, int patch_y_high,
                                       int patch_x_low, int patch_x_high,
                                       const mpp_config_t &cfg);

struct AlignGlobalResult {
	std::vector<cv::Vec2i> shifts;  /* (dy, dx) per frame, in PSS sign convention */
	cv::Vec4i patch_yxyx;           /* (y_low, y_high, x_low, x_high) on best frame */
	int best_frame_idx = -1;
};

/* Pick a patch on the best (max-quality) frame and compute per-frame integer
 * shifts via PSS's backward-then-forward cumulative-shift loop. `frames` are
 * the Gaussian-blurred mono frames (same as PSS's frames_mono_blurred). */
AlignGlobalResult align_global_from_frames(const std::vector<cv::Mat> &frames_mono_blurred,
                                           const std::vector<double> &quality,
                                           const mpp_config_t &cfg);

/* PSS rank_frames.find_best_frames: pick the `number_frames` indices that
 * maximise the summed quality within a sliding window of size `region_size`. */
std::vector<int> align_find_best_frames(const std::vector<double> &quality,
                                        int number_frames, int region_size);

struct AlignAverageResult {
	cv::Mat mean_frame;     /* CV_32S, scaled per PSS's uint8/uint16 convention */
	cv::Vec4i intersection; /* (y_low, y_high, x_low, x_high) in best-frame coords */
	std::vector<int> indices_used;
};

/* PSS align_frames.average_frame on `frames_mono_raw` (NOT blurred).
 * Uses the global `shifts` to align frames, picks the top-N by quality (with
 * the fast-changing-object sliding window if cfg.align_frames_fast_changing_object),
 * and averages within the inter-frame intersection. */
AlignAverageResult align_average_frame(const std::vector<cv::Mat> &frames_mono_raw,
                                       const std::vector<double> &quality,
                                       const std::vector<cv::Vec2i> &shifts,
                                       const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_ */
