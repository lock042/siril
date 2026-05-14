/*
 * Private C++ entry points for mpp_rank, exposed so unit tests can drive the
 * ranker directly with cv::Mat inputs without going through a Siril sequence.
 */
#ifndef SRC_REGISTRATION_MPP_RANK_PRIV_HPP_
#define SRC_REGISTRATION_MPP_RANK_PRIV_HPP_

#include <opencv2/core.hpp>

#include "registration/mpp.h"
#include "registration/mpp_config.h"

namespace mpp {

/* Score a single mono frame with PSS's Laplace pipeline:
 *   GaussianBlur(7x7, sigma=0)
 *   → stride-`align_frames_sampling_stride` sub-sample
 *   → Laplacian → convertScaleAbs(alpha=1/256)
 *   → stddev of resulting uint8 image.
 * Mono input may be CV_8U, CV_16U, or CV_32F. Output is the raw σ (NOT yet
 * normalized across the sequence by the max). */
double rank_score_mat(const cv::Mat &mono, const mpp_config_t &cfg);

/* Frame brightness used for the optional normalization step:
 *   threshold(mono, thr * scale, 255, THRESH_TOZERO) → cv::mean[0] + 1e-10
 * where scale = 256 for 16-bit input, 1 for 8-bit. */
double rank_average_brightness(const cv::Mat &mono, const mpp_config_t &cfg);

/* σ / brightness if cfg.frames_normalization, otherwise just σ. */
double rank_score_normalized(const cv::Mat &mono, const mpp_config_t &cfg);

/* PSS frames.frames_mono_blurred_laplacian: GaussianBlur(7×7) → stride-2
 * sub-sample → cv::Laplacian(CV_32F) → cv::convertScaleAbs(α=1/256).
 * Returns a CV_8U image. Per-AP quality ranking in Phase 5a samples patches
 * from this image. */
cv::Mat rank_blurred_laplacian_u8(const cv::Mat &mono, const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_RANK_PRIV_HPP_ */
