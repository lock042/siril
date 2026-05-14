/*
 * Phase 1: frame quality ranking.
 *
 * Ports PSS's "Laplace" ranking path (rank_frames.frame_score + the helper
 * Frames.frames_mono_blurred_laplacian). The algorithm is exactly:
 *
 *   blurred = cv::GaussianBlur(mono, (gauss_width, gauss_width), 0)
 *   strided = blurred[::sampling_stride, ::sampling_stride]    // NumPy slicing
 *   lap     = cv::Laplacian(strided, CV_32F)                   // default ksize=1
 *   lap_u8  = cv::convertScaleAbs(lap, alpha=1/256)
 *   σ       = cv::meanStdDev(lap_u8)[1][0][0]
 *
 * Optional brightness normalization divides σ by:
 *   cv::mean(cv::threshold(mono, thr, 255, THRESH_TOZERO))[0] + 1e-10
 * where thr = frames_normalization_threshold * 256 for 16-bit input,
 *           = frames_normalization_threshold        for  8-bit input.
 */

#include "registration/mpp.h"
#include "registration/mpp_rank.h"
#include "registration/mpp_rank_priv.hpp"

#include <cstring>
#include <opencv2/imgproc.hpp>

namespace mpp {

namespace {

/* NumPy-style [::s, ::s] sub-sampling. Output shape is ceil(N/s) along each
 * axis; element (y, x) is src(y*s, x*s). */
cv::Mat stride_subsample(const cv::Mat &src, int s) {
	if (s <= 1)
		return src;
	const int new_rows = (src.rows + s - 1) / s;
	const int new_cols = (src.cols + s - 1) / s;
	cv::Mat out(new_rows, new_cols, src.type());
	const int elem_size = (int) src.elemSize();
	for (int y = 0; y < new_rows; ++y) {
		const uchar *src_row = src.ptr<uchar>(y * s);
		uchar *dst_row = out.ptr<uchar>(y);
		for (int x = 0; x < new_cols; ++x) {
			std::memcpy(dst_row + x * elem_size,
			            src_row + (x * s) * elem_size,
			            elem_size);
		}
	}
	return out;
}

double brightness_threshold(const mpp_config_t &cfg) {
	return (double) cfg.frames_normalization_threshold * mpp_cfg_threshold_scale(&cfg);
}

}  // namespace

cv::Mat blur_mono_for_align(const cv::Mat &mono, const mpp_config_t &cfg) {
	/* PSS frames.frames_mono_blurred (frames.py:1505-1511): upscales 8-bit
	 * inputs to 16-bit range before the blur so downstream operations have
	 * useful Laplacian magnitude after the α=1/256 convertScaleAbs. */
	cv::Mat src = mono;
	if (cfg.bitdepth == 8 && mono.depth() == CV_8U) {
		cv::Mat upscaled;
		mono.convertTo(upscaled, CV_16U, 256.0);
		src = upscaled;
	}
	cv::Mat blurred;
	cv::GaussianBlur(src, blurred,
	                 cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	return blurred;
}

cv::Mat rank_blurred_laplacian_u8(const cv::Mat &mono, const mpp_config_t &cfg) {
	const cv::Mat blurred = blur_mono_for_align(mono, cfg);
	const cv::Mat sub = stride_subsample(blurred, cfg.align_frames_sampling_stride);
	cv::Mat lap;
	cv::Laplacian(sub, lap, CV_32F);
	cv::Mat lap_u8;
	cv::convertScaleAbs(lap, lap_u8, cfg.rank_laplacian_alpha);
	return lap_u8;
}

double rank_score_mat(const cv::Mat &mono, const mpp_config_t &cfg) {
	const cv::Mat lap_u8 = rank_blurred_laplacian_u8(mono, cfg);
	cv::Scalar mean, stddev;
	cv::meanStdDev(lap_u8, mean, stddev);
	return stddev[0];
}

double rank_average_brightness(const cv::Mat &mono, const mpp_config_t &cfg) {
	cv::Mat thresholded;
	cv::threshold(mono, thresholded, brightness_threshold(cfg), 255.0,
	              cv::THRESH_TOZERO);
	const cv::Scalar m = cv::mean(thresholded);
	return m[0] + 1e-10;
}

double rank_score_normalized(const cv::Mat &mono, const mpp_config_t &cfg) {
	const double sigma = rank_score_mat(mono, cfg);
	if (!cfg.frames_normalization)
		return sigma;
	const double ab = rank_average_brightness(mono, cfg);
	return sigma / ab;
}

}  // namespace mpp

/* ----------------- Public C interface ----------------- */

extern "C" mpp_status_t mpp_rank_sequence(struct sequence *seq,
                                          const mpp_config_t *cfg,
                                          double **quality_out) {
	(void) seq; (void) cfg;
	/* Sequence integration lands once mpp_frames is wired up (Phase 1.3). */
	if (quality_out) *quality_out = nullptr;
	return MPP_ENOTIMPL;
}
