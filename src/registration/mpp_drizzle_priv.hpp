/*
 * C++-API entry for the STScI stack path. Tests call this directly with
 * a vector<cv::Mat> instead of going through Siril's sequence I/O.
 */
#ifndef SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_
#define SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_config.h"

namespace mpp {

/* Inner STScI stack path. Same semantics as the extern "C"
 * mpp_stack_apply_stsci but operates on pre-loaded raw frames rather
 * than reading from a Siril sequence.
 *
 * `frames_raw[i]` is the input frame (single-channel mono or 3-channel
 * interleaved RGB, in OpenCV's R/G/B layout the way Phase 5a expects).
 * `included[i]` is 1 if frame i contributes to the stack, 0 otherwise
 * (must have length frames_raw.size()).
 * `frame_brightness[i]` is the per-frame average brightness used for
 * the brightness-equalisation pre-scale.
 *
 * `run` must have aps, shifts, global_shifts, intersection, and
 * frame_brightness populated. `cfg` drives drizzle_factor, pixfrac,
 * and kernel.
 *
 * Output dimensions: drizzle_factor * intersection_w × intersection_h,
 * num_layers channels.
 *
 * Returns MPP_OK on success; output uint16 stacked image filled into
 * `*out` (caller owns out->data — clearfits to free).
 */
mpp_status_t stack_apply_stsci(const std::vector<cv::Mat> &frames_raw,
                               const std::vector<int> &included,
                               const std::vector<double> &frame_brightness,
                               const mpp_run_t *run,
                               const mpp_config_t *cfg,
                               fits *out);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_ */
