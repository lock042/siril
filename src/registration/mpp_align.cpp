/*
 * Phase 2: global frame alignment.
 *
 * Ports PSS's compute_alignment_rect (patch picker, align_frames.py:85),
 * multilevel_correlation (miscellaneous.py:202), and the
 * backward-then-forward cumulative-shift loop in align_frames.align_frames.
 *
 * MultiLevelCorrelation summary:
 *   reference_window               = best_frame.mono_blurred[patch].astype(f32)
 *   reference_window_first_phase   = reference_window[::2, ::2]
 *   Phase 1: stride-2 sub-sample of the frame window, additional 7x7
 *     Gaussian blur, TM_CCORR_NORMED matchTemplate. Search half-width is
 *     (search_width - 4) / 2 on the coarse grid; final shift is doubled.
 *   Phase 2: full-resolution refinement, ±4 search around phase-1 result,
 *     TM_CCORR_NORMED matchTemplate again.
 *   Either phase is regarded as failed if the optimum lies on the search
 *     boundary, signalling that the true optimum is likely outside the box.
 *
 * Sign convention (matches PSS frame_shifts):
 *   shift_y, shift_x are what we need to apply to the *frame* so that it
 *   aligns with the reference frame.
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp_align.h"
#include "registration/mpp_align_priv.hpp"

namespace mpp {

namespace {

/* NumPy-style [::s, ::s] sub-sampling (also used by mpp_rank). */
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

/* PSS uses NumPy's `abs(frame[:, s2:] - frame[:, :-s2])` on uint16 inputs,
 * which underflows on negative differences (becomes huge positive). To match
 * its quality_measure_threshold_weighted bit-exact we replicate the modular
 * arithmetic. For float inputs we use signed math (no underflow). */
template <typename T>
double quality_sum_one_axis(const cv::Mat &frame, int axis_h_diff, int stride,
                            const cv::Mat &mask) {
	/* axis_h_diff = 1 → horizontal differences (left vs right by 2*stride cols).
	 * axis_h_diff = 0 → vertical (top vs bottom by 2*stride rows). */
	const int s2 = 2 * stride;
	const int H = frame.rows, W = frame.cols;
	double sum = 0.0;
	if (axis_h_diff) {
		/* diff = frame(:, s2:) - frame(:, :-s2), mask = mask(:, s:-s) */
		for (int y = 0; y < H; ++y) {
			const T *row = frame.ptr<T>(y);
			const uchar *mrow = mask.ptr<uchar>(y);
			for (int x = 0; x < W - s2; ++x) {
				if (mrow[x + stride]) {
					const T d = (T) (row[x + s2] - row[x]);
					sum += (double) d;
				}
			}
		}
	} else {
		for (int y = 0; y < H - s2; ++y) {
			const T *row_top = frame.ptr<T>(y);
			const T *row_bot = frame.ptr<T>(y + s2);
			const uchar *mrow = mask.ptr<uchar>(y + stride);
			for (int x = 0; x < W; ++x) {
				if (mrow[x]) {
					const T d = (T) (row_bot[x] - row_top[x]);
					sum += (double) d;
				}
			}
		}
	}
	return sum;
}

}  // namespace

double align_quality_measure_threshold_weighted(const cv::Mat &frame, int stride,
                                                double black_threshold,
                                                double min_fraction) {
	const int H = frame.rows, W = frame.cols;
	const double frame_size = (double) H * (double) W;

	/* PSS: mask = frame > black_threshold; mask_fraction = mask.sum() / frame.size */
	cv::Mat mask;
	cv::compare(frame, black_threshold, mask, cv::CMP_GT);  /* CV_8U, 0 or 255 */
	const int mask_nonzero = cv::countNonZero(mask);
	const double mask_fraction = mask_nonzero / frame_size;

	/* `mask` from cv::compare is 0/255; downstream code treats it as truthy. */
	double sum_h = 0.0, sum_v = 0.0;
	if (frame.type() == CV_16U) {
		sum_h = quality_sum_one_axis<uint16_t>(frame, 1, stride, mask);
		sum_v = quality_sum_one_axis<uint16_t>(frame, 0, stride, mask);
	} else if (frame.type() == CV_8U) {
		sum_h = quality_sum_one_axis<uint8_t>(frame, 1, stride, mask);
		sum_v = quality_sum_one_axis<uint8_t>(frame, 0, stride, mask);
	} else if (frame.type() == CV_32F) {
		/* For float, PSS's NumPy abs() is correct; use signed math + abs. */
		const int s2 = 2 * stride;
		for (int y = 0; y < H; ++y) {
			const float *row = frame.ptr<float>(y);
			const uchar *mrow = mask.ptr<uchar>(y);
			for (int x = 0; x < W - s2; ++x)
				if (mrow[x + stride])
					sum_h += std::abs(row[x + s2] - row[x]);
		}
		for (int y = 0; y < H - s2; ++y) {
			const float *row_top = frame.ptr<float>(y);
			const float *row_bot = frame.ptr<float>(y + s2);
			const uchar *mrow = mask.ptr<uchar>(y + stride);
			for (int x = 0; x < W; ++x)
				if (mrow[x])
					sum_v += std::abs(row_bot[x] - row_top[x]);
		}
	} else {
		return 0.0;  /* Unsupported input type for this PSS pathway. */
	}

	if (mask_fraction > min_fraction) {
		sum_h /= mask_fraction;
		sum_v /= mask_fraction;
	}
	return std::min(sum_h, sum_v);
}

cv::Vec4i align_pick_patch(const cv::Mat &best, const mpp_config_t &cfg) {
	const int dim_y = best.rows, dim_x = best.cols;
	const int border = cfg.align_frames_border_width + cfg.align_frames_search_width;
	const double scale = cfg.align_frames_rectangle_scale_factor;

	const int rect_y = (int) ((dim_y - 2 * border) / scale);
	const int rect_x = (int) ((dim_x - 2 * border) / scale);
	const int step_y = rect_y / 2;
	const int step_x = rect_x / 2;
	if (rect_y <= 0 || rect_x <= 0 || step_y <= 0 || step_x <= 0)
		return cv::Vec4i(0, 0, 0, 0);

	double best_quality = -1.0;
	cv::Vec4i best_bounds(border, border + rect_y, border, border + rect_x);

	int x_low = border;
	while (x_low + rect_x <= dim_x - border) {
		int y_low = border;
		while (y_low + rect_y <= dim_y - border) {
			const cv::Mat patch = best(cv::Range(y_low, y_low + rect_y),
			                           cv::Range(x_low, x_low + rect_x));
			const double q = align_quality_measure_threshold_weighted(
			    patch, cfg.align_frames_rectangle_stride,
			    (double) cfg.align_frames_rectangle_black_threshold,
			    cfg.align_frames_rectangle_min_fraction);
			/* PSS uses a stable Python sorted(); std::stable_sort would give the
			 * same answer, but here we just need argmax with a strict-greater
			 * tiebreak that picks the *first* candidate seen (matching Python's
			 * "reverse=True" sort on a list built in the loop order below). */
			if (q > best_quality) {
				best_quality = q;
				best_bounds = cv::Vec4i(y_low, y_low + rect_y, x_low, x_low + rect_x);
			}
			y_low += step_y;
		}
		x_low += step_x;
	}
	return best_bounds;
}

/* PSS sub_pixel_solve_matrix from miscellaneous.py:464 — precomputed
 * pseudo-inverse `(AᵀA)⁻¹ Aᵀ` for fitting f(x,y) = a x² + b y² + c xy + d x + e y + g
 * to 9 correlation values on a 3×3 grid (indexed row-major y=0..2, x=0..2 with
 * the grid centred at (0,0)). */
static const double SUB_PIXEL_SOLVE[6][9] = {
	{ 0.16666667, -0.33333333,  0.16666667,  0.16666667, -0.33333333,  0.16666667,  0.16666667, -0.33333333,  0.16666667},
	{ 0.16666667,  0.16666667,  0.16666667, -0.33333333, -0.33333333, -0.33333333,  0.16666667,  0.16666667,  0.16666667},
	{ 0.25,        0.,         -0.25,        0.,          0.,          0.,         -0.25,        0.,          0.25      },
	{-0.16666667,  0.,          0.16666667, -0.16666667,  0.,          0.16666667, -0.16666667,  0.,          0.16666667},
	{-0.16666667, -0.16666667, -0.16666667,  0.,          0.,          0.,          0.16666667,  0.16666667,  0.16666667},
	{-0.11111111,  0.22222222, -0.11111111,  0.22222222,  0.55555556,  0.22222222, -0.11111111,  0.22222222, -0.11111111}
};

namespace {

/* Fit a quadratic surface to a 3×3 patch of correlation values and return
 * the (y_corr, x_corr) offset to the peak. Returns false on a near-singular
 * solve. Mirrors PSS Miscellaneous.sub_pixel_solve. */
bool sub_pixel_solve(const float vals[9], double *y_corr, double *x_corr) {
	double coeffs[6] = {0, 0, 0, 0, 0, 0};
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 9; ++j)
			coeffs[i] += SUB_PIXEL_SOLVE[i][j] * (double) vals[j];
	const double a_f = coeffs[0], b_f = coeffs[1], c_f = coeffs[2];
	const double d_f = coeffs[3], e_f = coeffs[4];
	const double denom_y = c_f * c_f - 4.0 * a_f * b_f;
	if (std::abs(denom_y) > 1e-10 && std::abs(a_f) > 1e-10) {
		const double yc = (2.0 * a_f * e_f - c_f * d_f) / denom_y;
		*y_corr = yc;
		*x_corr = (-c_f * yc - d_f) / (2.0 * a_f);
		return true;
	}
	if (std::abs(denom_y) > 1e-10 && std::abs(c_f) > 1e-10) {
		const double yc = (2.0 * a_f * e_f - c_f * d_f) / denom_y;
		*y_corr = yc;
		*x_corr = (-2.0 * b_f * yc - e_f) / c_f;
		return true;
	}
	return false;
}

}  // namespace

MultilevelShiftResult multilevel_correlation(const cv::Mat &ref_full_f32,
                                             const cv::Mat &ref_first_phase_f32,
                                             const cv::Mat &frame_mono_blurred,
                                             int y_low, int y_high,
                                             int x_low, int x_high,
                                             int gauss_width, int search_width,
                                             bool subpixel_solve,
                                             const cv::Mat &weight_matrix_first_phase) {
	MultilevelShiftResult result;
	const int sw2 = 4;
	const int sw1 = (search_width - sw2) / 2;
	const int index_ext = sw1 * 2;
	const int H = frame_mono_blurred.rows, W = frame_mono_blurred.cols;

	/* Phase 1: stride-2 frame window + GaussianBlur, full-frame reference at
	 * stride-2. */
	const int y_lo1 = y_low - index_ext, y_hi1 = y_high + index_ext;
	const int x_lo1 = x_low - index_ext, x_hi1 = x_high + index_ext;
	if (y_lo1 < 0 || x_lo1 < 0 || y_hi1 > H || x_hi1 > W)
		return result;
	const cv::Mat fwin = frame_mono_blurred(cv::Range(y_lo1, y_hi1),
	                                        cv::Range(x_lo1, x_hi1));
	const cv::Mat fwin_strided = stride_subsample(fwin, 2);
	cv::Mat fwin_blurred;
	cv::GaussianBlur(fwin_strided, fwin_blurred,
	                 cv::Size(gauss_width, gauss_width), 0);
	cv::Mat fwin_f32;
	fwin_blurred.convertTo(fwin_f32, CV_32F);

	cv::Mat ccr1;
	cv::matchTemplate(fwin_f32, ref_first_phase_f32, ccr1, cv::TM_CCORR_NORMED);
	cv::Point max1;
	if (!weight_matrix_first_phase.empty()) {
		cv::Mat weighted;
		cv::multiply(ccr1, weight_matrix_first_phase, weighted);
		cv::minMaxLoc(weighted, nullptr, nullptr, nullptr, &max1);
	} else {
		cv::minMaxLoc(ccr1, nullptr, nullptr, nullptr, &max1);
	}
	const int shift_y1 = (sw1 - max1.y) * 2;
	const int shift_x1 = (sw1 - max1.x) * 2;
	if (std::abs(shift_y1) == index_ext || std::abs(shift_x1) == index_ext)
		return result;

	/* Phase 2: ±4 around phase-1 result, full resolution. */
	const int y_lo2 = y_low - shift_y1 - sw2;
	const int y_hi2 = y_high - shift_y1 + sw2;
	const int x_lo2 = x_low - shift_x1 - sw2;
	const int x_hi2 = x_high - shift_x1 + sw2;
	if (y_lo2 < 0 || x_lo2 < 0 || y_hi2 > H || x_hi2 > W)
		return result;
	const cv::Mat fwin2 = frame_mono_blurred(cv::Range(y_lo2, y_hi2),
	                                         cv::Range(x_lo2, x_hi2));
	cv::Mat fwin2_f32;
	fwin2.convertTo(fwin2_f32, CV_32F);

	cv::Mat ccr2;
	cv::matchTemplate(fwin2_f32, ref_full_f32, ccr2, cv::TM_CCORR_NORMED);
	cv::Point max2;
	cv::minMaxLoc(ccr2, nullptr, nullptr, nullptr, &max2);
	const int shift_y2 = sw2 - max2.y;
	const int shift_x2 = sw2 - max2.x;
	if (std::abs(shift_y2) == sw2 || std::abs(shift_x2) == sw2)
		return result;

	double y_total = (double) (shift_y1 + shift_y2);
	double x_total = (double) (shift_x1 + shift_x2);

	if (subpixel_solve && max2.x > 0 && max2.x < ccr2.cols - 1
	                  && max2.y > 0 && max2.y < ccr2.rows - 1) {
		float vals[9];
		for (int j = 0; j < 3; ++j)
			for (int i = 0; i < 3; ++i)
				vals[j * 3 + i] = ccr2.at<float>(max2.y - 1 + j, max2.x - 1 + i);
		double y_corr = 0.0, x_corr = 0.0;
		if (sub_pixel_solve(vals, &y_corr, &x_corr)
		    && std::abs(y_corr) <= 1.0 && std::abs(x_corr) <= 1.0) {
			y_total -= y_corr;
			x_total -= x_corr;
		}
	}

	result.dy = y_total;
	result.dx = x_total;
	result.success = true;
	return result;
}

AlignShiftResult align_shift_one_frame(const cv::Mat &ref_window_f32,
                                       const cv::Mat &ref_window_first_phase_f32,
                                       const cv::Mat &frame_mono_blurred,
                                       int y_low, int y_high,
                                       int x_low, int x_high,
                                       const mpp_config_t &cfg) {
	const MultilevelShiftResult r = multilevel_correlation(
	    ref_window_f32, ref_window_first_phase_f32, frame_mono_blurred,
	    y_low, y_high, x_low, x_high,
	    cfg.frames_gauss_width, cfg.align_frames_search_width,
	    /*subpixel_solve=*/false);
	AlignShiftResult out;
	out.dy = (int) std::lround(r.dy);
	out.dx = (int) std::lround(r.dx);
	out.success = r.success;
	return out;
}

AlignGlobalResult align_global_from_provider(const FrameProvider &provider,
                                             int N,
                                             const std::vector<double> &quality,
                                             const mpp_config_t &cfg,
                                             progress_cb_fn progress,
                                             void *progress_user) {
	AlignGlobalResult out;
	if (N <= 0 || (int) quality.size() != N)
		return out;
	int done = 0;
	const int total = N;   /* every frame is visited once across the two passes */

	/* Best frame is argmax of quality. */
	int best = 0;
	for (int i = 1; i < N; ++i)
		if (quality[i] > quality[best]) best = i;
	out.best_frame_idx = best;

	/* Patch on best frame's mono_blurred. Hold the best frame for the
	 * full duration of both sweeps — we need to re-derive ref_window
	 * from it (a view, not a copy). */
	const cv::Mat best_blurred = provider(best);
	const cv::Vec4i patch = align_pick_patch(best_blurred, cfg);
	out.patch_yxyx = patch;
	const int py_lo = patch[0], py_hi = patch[1];
	const int px_lo = patch[2], px_hi = patch[3];

	/* Reference windows: best frame mono_blurred cropped & cast to float32. */
	cv::Mat best_f32;
	best_blurred.convertTo(best_f32, CV_32F);
	const cv::Mat ref_window = best_f32(cv::Range(py_lo, py_hi),
	                                    cv::Range(px_lo, px_hi));
	const cv::Mat ref_window_first = stride_subsample(ref_window, 2);

	out.shifts.assign(N, cv::Vec2i(0, 0));

	/* Backward pass: max_idx, max_idx-1, ..., 0. */
	{
		int dy_cum = 0, dx_cum = 0;
		for (int idx = best; idx >= 0; --idx) {
			if (idx == best) { out.shifts[idx] = {0, 0}; ++done; continue; }
			const cv::Mat frame = provider(idx);
			const AlignShiftResult r = align_shift_one_frame(
			    ref_window, ref_window_first, frame,
			    py_lo - dy_cum, py_hi - dy_cum,
			    px_lo - dx_cum, px_hi - dx_cum, cfg);
			if (!r.success) {
				/* PSS raises InternalError; tests will see a zero shift and
				 * notice the divergence. */
				out.shifts[idx] = {dy_cum, dx_cum};
				++done;
				if (progress) progress((double) done / (double) total, progress_user);
				continue;
			}
			dy_cum += r.dy;
			dx_cum += r.dx;
			out.shifts[idx] = {dy_cum, dx_cum};
			++done;
			if (progress) progress((double) done / (double) total, progress_user);
		}
	}

	/* Forward pass: max_idx+1, ..., N-1. PSS re-enters at max_idx but only to
	 * reset dy_cum; we skip the redundant visit. */
	{
		int dy_cum = 0, dx_cum = 0;
		for (int idx = best + 1; idx < N; ++idx) {
			const cv::Mat frame = provider(idx);
			const AlignShiftResult r = align_shift_one_frame(
			    ref_window, ref_window_first, frame,
			    py_lo - dy_cum, py_hi - dy_cum,
			    px_lo - dx_cum, px_hi - dx_cum, cfg);
			if (!r.success) {
				out.shifts[idx] = {dy_cum, dx_cum};
				++done;
				if (progress) progress((double) done / (double) total, progress_user);
				continue;
			}
			dy_cum += r.dy;
			dx_cum += r.dx;
			out.shifts[idx] = {dy_cum, dx_cum};
			++done;
			if (progress) progress((double) done / (double) total, progress_user);
		}
	}

	return out;
}

AlignGlobalResult align_global_from_frames(const std::vector<cv::Mat> &frames,
                                           const std::vector<double> &quality,
                                           const mpp_config_t &cfg,
                                           progress_cb_fn progress,
                                           void *progress_user) {
	return align_global_from_provider(
	    [&frames](int i) -> cv::Mat { return frames[i]; },
	    (int) frames.size(), quality, cfg, progress, progress_user);
}

std::vector<int> align_find_best_frames(const std::vector<double> &quality,
                                        int number_frames, int region_size) {
	const int N = (int) quality.size();
	if (number_frames > region_size || region_size > N || number_frames <= 0)
		return {};
	std::vector<int> best;
	double best_rank_sum = 0.0;
	std::vector<int> in_window(region_size);
	for (int start = 0; start <= N - region_size; ++start) {
		std::iota(in_window.begin(), in_window.end(), start);
		std::partial_sort(in_window.begin(), in_window.begin() + number_frames,
		                  in_window.end(),
		                  [&quality](int a, int b) { return quality[a] > quality[b]; });
		double s = 0.0;
		for (int i = 0; i < number_frames; ++i) s += quality[in_window[i]];
		if (s > best_rank_sum) {
			best_rank_sum = s;
			best.assign(in_window.begin(), in_window.begin() + number_frames);
		}
	}
	return best;
}

AlignAverageResult align_average_frame_streamed(const FrameProvider &provider,
                                                int N,
                                                int H, int W,
                                                const std::vector<double> &quality,
                                                const std::vector<cv::Vec2i> &shifts,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress,
                                                void *progress_user) {
	AlignAverageResult out;
	if (N == 0 || (int) shifts.size() != N || (int) quality.size() != N)
		return out;

	int max_dy = INT_MIN, min_dy = INT_MAX, max_dx = INT_MIN, min_dx = INT_MAX;
	for (const auto &s : shifts) {
		max_dy = std::max(max_dy, s[0]);
		min_dy = std::min(min_dy, s[0]);
		max_dx = std::max(max_dx, s[1]);
		min_dx = std::min(min_dx, s[1]);
	}
	out.intersection = cv::Vec4i(max_dy, min_dy + H, max_dx, min_dx + W);
	const int int_y_lo = out.intersection[0], int_y_hi = out.intersection[1];
	const int int_x_lo = out.intersection[2], int_x_hi = out.intersection[3];
	if (int_y_hi <= int_y_lo || int_x_hi <= int_x_lo)
		return out;

	const int n_target = std::max(
	    (int) std::ceil(N * cfg.align_frames_average_frame_percent / 100.0), 1);

	std::vector<int> indices;
	if (cfg.align_frames_fast_changing_object) {
		const int window = std::min(
		    n_target * cfg.align_frames_best_frames_window_extension, N);
		indices = align_find_best_frames(quality, n_target, window);
	} else {
		indices.resize(N);
		std::iota(indices.begin(), indices.end(), 0);
		std::partial_sort(indices.begin(), indices.begin() + n_target, indices.end(),
		                  [&quality](int a, int b) { return quality[a] > quality[b]; });
		indices.resize(n_target);
	}
	out.indices_used = indices;

	const int oh = int_y_hi - int_y_lo, ow = int_x_hi - int_x_lo;
	cv::Mat accum(oh, ow, CV_32F, cv::Scalar(0));
	int done = 0;
	const int total = (int) indices.size();
	for (int idx : indices) {
		const auto &s = shifts[idx];
		const cv::Mat frame = provider(idx);
		if (frame.empty()) { ++done; continue; }
		const cv::Mat patch = frame(
		    cv::Range(int_y_lo - s[0], int_y_hi - s[0]),
		    cv::Range(int_x_lo - s[1], int_x_hi - s[1]));
		cv::Mat patch_f32;
		patch.convertTo(patch_f32, CV_32F);
		accum += patch_f32;
		++done;
		if (progress) progress((double) done / (double) total, progress_user);
	}

	/* PSS scales by 256/N for uint8 inputs and 1/N for uint16. Siril stores
	 * 8-bit SER data as WORD without upscaling (values 0..255), so we look
	 * at cfg.bitdepth — NOT cv::Mat::depth() — to decide. Output mean_frame
	 * always lands in the 0..65535 range, matching PSS so downstream AP
	 * code can use a single `× 256` threshold scale. */
	const double scale = (cfg.bitdepth == 8 ? 256.0 : 1.0) / (double) indices.size();
	accum *= scale;
	/* PSS does `.astype(int32)`, which truncates toward zero (Python/NumPy
	 * convention). OpenCV's convertTo(CV_32S) uses cvRound (round half to
	 * even on x86), so we'd get off-by-one diffs at half-pixel values and
	 * those propagate into the per-AP correlation peaks downstream. Do the
	 * truncation manually to keep bit-equivalence. */
	out.mean_frame.create(accum.size(), CV_32S);
	for (int y = 0; y < accum.rows; ++y) {
		const float *src = accum.ptr<float>(y);
		int32_t *dst = out.mean_frame.ptr<int32_t>(y);
		for (int x = 0; x < accum.cols; ++x)
			dst[x] = (int32_t) src[x];  /* truncation toward zero */
	}
	return out;
}

AlignAverageResult align_average_frame(const std::vector<cv::Mat> &frames_raw,
                                       const std::vector<double> &quality,
                                       const std::vector<cv::Vec2i> &shifts,
                                       const mpp_config_t &cfg,
                                       progress_cb_fn progress,
                                       void *progress_user) {
	const int N = (int) frames_raw.size();
	const int H = N > 0 ? frames_raw[0].rows : 0;
	const int W = N > 0 ? frames_raw[0].cols : 0;
	return align_average_frame_streamed(
	    [&frames_raw](int i) -> cv::Mat { return frames_raw[i]; },
	    N, H, W, quality, shifts, cfg, progress, progress_user);
}

}  // namespace mpp

/* ----------------- Public C interface ----------------- */

extern "C" mpp_status_t mpp_align_global(sequence *seq,
                                         const mpp_config_t *cfg,
                                         const double *quality,
                                         int *shifts_out,
                                         int patch_yxyx_out[4],
                                         fits *avg_ref_out) {
	(void) seq; (void) cfg; (void) quality;
	(void) shifts_out; (void) patch_yxyx_out; (void) avg_ref_out;
	/* Sequence integration lands with mpp_frames (Phase 1.3 / 2). */
	return MPP_ENOTIMPL;
}
