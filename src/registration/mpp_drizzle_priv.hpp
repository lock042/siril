/*
 * C++-API entry for the STScI stack path. Tests call this directly with
 * a vector<cv::Mat> instead of going through Siril's sequence I/O.
 */
#ifndef SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_
#define SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_align_priv.hpp"   /* FrameProvider */
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
 * frame_brightness populated. `cfg` drives drizzle_scale, pixfrac,
 * and kernel.
 *
 * Output dimensions: drizzle_scale * intersection_w × intersection_h,
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

/* Streamed overload. provider(f) returns the raw frame (channel count
 * matching run->num_layers). num_frames must equal included.size().
 * One frame in flight at a time; sequential. */
mpp_status_t stack_apply_stsci_streamed(const FrameProvider &provider,
                                        int num_frames,
                                        const std::vector<int> &included,
                                        const std::vector<double> &frame_brightness,
                                        const mpp_run_t *run,
                                        const mpp_config_t *cfg,
                                        fits *out);

/* Inner Bayer-drizzle stack path. Inputs: per-frame single-channel raw
 * Bayer-mosaic frames (NO debayering applied), and the same run state
 * that Stage A produced on the corresponding debayered frames. dobox's
 * CFA branch routes each input pixel to the appropriate channel of a
 * 3-channel output canvas via `cfa[]` + `cfadim`, so the geometry from
 * the run's pixmap is reused unchanged.
 *
 * `cfa[]` is a 4-element BYTE array (2x2 pattern): RLAYER=0, GLAYER=1,
 * BLAYER=2 — same encoding Siril's get_compiled_pattern produces.
 * `cfadim` must be 2 (only 2x2 Bayer is supported; X-Trans uses cfadim=6
 * but isn't tested by this path yet).
 *
 * Output is always 3-channel uint16 (planar R/G/B with pdata[0..2] set). */
mpp_status_t stack_apply_bayer(const std::vector<cv::Mat> &frames_raw_bayer,
                               const std::vector<int> &included,
                               const std::vector<double> &frame_brightness,
                               const mpp_run_t *run,
                               const mpp_config_t *cfg,
                               const unsigned char *cfa,
                               int cfadim,
                               fits *out);

/* Streamed Bayer overload. provider(f) returns the single-channel raw
 * mosaic frame; num_frames must equal included.size(). One frame in
 * flight at a time; sequential. */
mpp_status_t stack_apply_bayer_streamed(const FrameProvider &provider,
                                        int num_frames,
                                        const std::vector<int> &included,
                                        const std::vector<double> &frame_brightness,
                                        const mpp_run_t *run,
                                        const mpp_config_t *cfg,
                                        const unsigned char *cfa,
                                        int cfadim,
                                        fits *out);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_DRIZZLE_PRIV_HPP_ */
