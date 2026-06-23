/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 *
 * This implementation of multipoint registration & stacking is based on
 * PlanetarySystemStacker by Rolf Hempel:
 *     https://github.com/Rolf-Hempel/PlanetarySystemStacker
 */

/*
 * Phase 2: global frame alignment.
 *
 * Picks the best-structured patch on the best frame, then aligns every
 * other frame to it via two-phase cross-correlation and a backward-then-
 * forward cumulative-shift loop.
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
 * Sign convention:
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

#include "core/processing.h"
#include "registration/mpp.h"
#include "registration/mpp/mpp_align.h"
#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_image_priv.hpp"

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

/* The reference algorithm uses NumPy's
 * `abs(frame[:, s2:] - frame[:, :-s2])` on uint16 inputs, which underflows
 * on negative differences (becomes huge positive). To stay bit-exact
 * with the reference we replicate the modular arithmetic here. For
 * float inputs we use signed math (no underflow). */
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

	/* mask = frame > black_threshold; mask_fraction = mask.sum() / frame.size */
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
		/* For float, NumPy's abs() is well-defined; use signed math + abs. */
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
		return 0.0;  /* Unsupported input type for this pathway. */
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
			/* argmax with a strict-greater tiebreak picks the *first*
			 * candidate seen, matching a stable descending sort on a
			 * list built in the loop order below. */
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

cv::Vec2d center_of_gravity(const cv::Mat &frame_mono_blurred) {
	/* PSS AlignFrames.center_of_gravity (align_frames.py):
	 *   minVal, maxVal = minMaxLoc(frame)
	 *   threshold = int((minVal + maxVal) / 2)        # truncates toward zero
	 *   thresh    = clip(frame, threshold, None) - threshold   # = max(frame-thr, 0)
	 *   M = moments(thresh); cog = (M01/M00, M10/M00)
	 * PSS then round()s cog to integer; we keep it sub-pixel (see header). */
	double min_val = 0.0, max_val = 0.0;
	cv::minMaxLoc(frame_mono_blurred, &min_val, &max_val);
	const double threshold = std::trunc((min_val + max_val) / 2.0);

	cv::Mat weights;
	frame_mono_blurred.convertTo(weights, CV_64F);
	weights -= threshold;
	cv::max(weights, 0.0, weights);   /* clip(frame, thr, None) - thr */

	const cv::Moments m = cv::moments(weights, /*binaryImage=*/false);
	if (m.m00 <= 0.0) {
		/* Blank/flat frame: no centroid. Fall back to the geometric
		 * centre so a uniform sequence yields zero relative shift. */
		return cv::Vec2d((frame_mono_blurred.rows - 1) / 2.0,
		                 (frame_mono_blurred.cols - 1) / 2.0);
	}
	const double cog_x = m.m10 / m.m00;
	const double cog_y = m.m01 / m.m00;
	return cv::Vec2d(cog_y, cog_x);
}

/* Precomputed pseudo-inverse `(AᵀA)⁻¹ Aᵀ` for fitting
 * f(x,y) = a x² + b y² + c xy + d x + e y + g to 9 correlation values
 * on a 3×3 grid (indexed row-major y=0..2, x=0..2 with the grid centred
 * at (0,0)). */
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
 * solve. */
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

	/* From here on phase 1 has produced a usable coarse estimate: report
	 * it even when phase 2 fails, with success=false.  Upstream PSS keeps
	 * the phase-1 shift in exactly this way (miscellaneous.py
	 * multilevel_correlation zeroes only the SECOND-phase component on a
	 * phase-2 border hit) and its stacker then stacks the patch at the
	 * phase-1 estimate.  Returning {0,0} here instead — as this port used
	 * to — made every phase-2 failure stack its patch at global-only
	 * alignment, several pixels off, which shows up as ghosted AP boxes
	 * on difficult data.  Callers that must not consume a coarse-only
	 * estimate (the global-align sweeps) already gate on `success`. */
	result.dy = (double) shift_y1;
	result.dx = (double) shift_x1;

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
	/* Sub-pixel ON unconditionally: the parabolic-fit residual is the
	 * only thing that gives Bayer drizzle's cross-frame CFA-phase
	 * coverage. Cheap (4 extra reads per frame on the phase-2 correlation
	 * surface, gated on the same success path). The previous integer
	 * rounding wasted that information. */
	const MultilevelShiftResult r = multilevel_correlation(
	    ref_window_f32, ref_window_first_phase_f32, frame_mono_blurred,
	    y_low, y_high, x_low, x_high,
	    cfg.frames_gauss_width, cfg.align_frames_search_width,
	    /*subpixel_solve=*/true);
	AlignShiftResult out;
	out.dy = r.dy;
	out.dx = r.dx;
	out.success = r.success;
	return out;
}

AlignGlobalResult align_global_from_provider(const FrameProvider &provider,
                                             int N,
                                             const std::vector<double> &quality,
                                             const mpp_config_t &cfg,
                                             progress_cb_fn progress,
                                             void *progress_user,
                                             bool provider_is_thread_safe,
                                             const std::vector<cv::Vec2d> &seed) {
	AlignGlobalResult out;
	if (N <= 0 || (int) quality.size() != N)
		return out;
	gint done = 0;
	const int total = N;   /* every frame is visited once across the two passes */

	/* Best frame is argmax of quality. */
	int best = 0;
	for (int i = 1; i < N; ++i)
		if (quality[i] > quality[best]) best = i;
	out.best_frame_idx = best;

	/* Planet mode (PSS "Planet"): align on the brightness centroid rather
	 * than correlating a surface patch. The per-frame shift is the centroid
	 * displacement relative to the best frame — this "cannot fail" (no
	 * search window, no patch). Frames are independent, so a single
	 * parallel-for over all frames replaces the two cumulative sweeps. */
	if (cfg.align_frames_mode == MPP_ALIGN_PLANET) {
		const cv::Mat best_blurred = provider(best);
		if (best_blurred.empty()) {        /* best frame unreadable */
			out.best_frame_idx = -1;
			return out;
		}
		const cv::Vec2d cog_ref = center_of_gravity(best_blurred);
		out.patch_yxyx = cv::Vec4i(0, 0, 0, 0);  /* no patch in planet mode */
		out.shifts.assign(N, cv::Vec2d(0.0, 0.0));

		gint cancelled = 0;
		gint oom = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (provider_is_thread_safe)
#endif
		for (int idx = 0; idx < N; ++idx) {
			if (g_atomic_int_get(&cancelled) || g_atomic_int_get(&oom)) continue;
			if (!processing_should_continue()) {
				g_atomic_int_set(&cancelled, 1);
				continue;
			}
			/* Contain allocation failures — an exception escaping an
			 * OpenMP structured block is std::terminate(). */
			try {
				if (idx != best) {
					const cv::Mat frame = provider(idx);
					/* Empty frame → leave the zero-init shift (excluded). */
					if (!frame.empty()) {
						const cv::Vec2d cog = center_of_gravity(frame);
						out.shifts[idx] = cv::Vec2d(cog_ref[0] - cog[0],
						                            cog_ref[1] - cog[1]);
					}
				}
			} catch (const std::exception &) {
				g_atomic_int_set(&oom, 1);
				continue;
			}
			const gint d = g_atomic_int_add(&done, 1) + 1;
			if (progress) progress((double) d / (double) total, progress_user);
		}
		if (g_atomic_int_get(&oom)) {
			out.best_frame_idx = -1;
			out.oom = true;
			return out;
		}
		if (g_atomic_int_get(&cancelled)) {
			out.best_frame_idx = -1;
			return out;
		}
		return out;
	}

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

	out.shifts.assign(N, cv::Vec2d(0.0, 0.0));

	/* Optional gross-shift seed (see header). When sized to N, each frame's
	 * search window is centred at its own seed re-referenced to the best
	 * frame, so the cumulative frame-to-frame chain is bypassed entirely —
	 * frames become independent and arbitrarily large jumps are tolerated
	 * as long as the seed's residual error is within search_width. */
	const bool use_seed = ((int) seed.size() == N);
	const cv::Vec2d seed_base = use_seed ? seed[best] : cv::Vec2d(0.0, 0.0);

	/* Backward and forward sweeps are independent: they read disjoint
	 * frame indices, share only const inputs (ref_window*, py/px_lo/hi,
	 * cfg) and write to disjoint slices of out.shifts. They each carry
	 * their own dy_cum/dx_cum chain. With a thread-safe provider we run
	 * them as two OpenMP sections for ~2x speedup on this stage when
	 * `best` lands mid-sequence. SEQ_AVI providers don't get the
	 * upgrade (ffms2 isn't reentrant); the if-clause collapses the
	 * parallel region to a single thread when provider_is_thread_safe
	 * is false. Cancellation and progress go through atomics so neither
	 * thread blocks the other or trips a memory race.
	 *
	 * Cumulation is integer-only — the search window slice (py_lo -
	 * dy_int_cum, …) needs integer indices, and accumulating the
	 * parabolic-fit residual would drift over N frames (each step's
	 * residual is independent noise, not signal). The fractional part
	 * is added at store-time so each stored shift = (sum of integer
	 * peaks up to here) + (this frame's local refinement). */
	gint cancelled = 0;
	gint oom = 0;
	/* Each sweep contains allocation failures locally (try/catch around
	 * its whole loop): the sweeps run as OpenMP sections, and an
	 * exception escaping a section is std::terminate(). An aborted sweep
	 * leaves the remaining shifts at their zero init; the oom flag below
	 * discards the whole result, so partial data never leaks out. */
	auto run_backward_sweep = [&] {
		try {
		int dy_int_cum = 0, dx_int_cum = 0;
		for (int idx = best; idx >= 0; --idx) {
			if (g_atomic_int_get(&cancelled)) break;
			if (!processing_should_continue()) {
				g_atomic_int_set(&cancelled, 1);
				break;
			}
			if (idx == best) {
				out.shifts[idx] = cv::Vec2d(0.0, 0.0);
				const gint d = g_atomic_int_add(&done, 1) + 1;
				if (progress) progress((double) d / (double) total, progress_user);
				continue;
			}
			const cv::Mat frame = provider(idx);
			if (use_seed) {
				/* Independent per-frame: window centred at this frame's
				 * seed (re-referenced to best); total = seed + residual. */
				const int sdy = (int) std::lround(seed[idx][0] - seed_base[0]);
				const int sdx = (int) std::lround(seed[idx][1] - seed_base[1]);
				const AlignShiftResult rs = align_shift_one_frame(
				    ref_window, ref_window_first, frame,
				    py_lo - sdy, py_hi - sdy, px_lo - sdx, px_hi - sdx, cfg);
				double dy_s = (double) sdy, dx_s = (double) sdx;
				if (rs.success) { dy_s = (double) sdy + rs.dy; dx_s = (double) sdx + rs.dx; }
				out.shifts[idx] = cv::Vec2d(dy_s, dx_s);
				const gint d = g_atomic_int_add(&done, 1) + 1;
				if (progress) progress((double) d / (double) total, progress_user);
				continue;
			}
			const AlignShiftResult r = align_shift_one_frame(
			    ref_window, ref_window_first, frame,
			    py_lo - dy_int_cum, py_hi - dy_int_cum,
			    px_lo - dx_int_cum, px_hi - dx_int_cum, cfg);
			double dy_store = (double) dy_int_cum;
			double dx_store = (double) dx_int_cum;
			if (r.success) {
				const int dy_int = (int) std::lround(r.dy);
				const int dx_int = (int) std::lround(r.dx);
				dy_int_cum += dy_int;
				dx_int_cum += dx_int;
				dy_store = (double) dy_int_cum + (r.dy - (double) dy_int);
				dx_store = (double) dx_int_cum + (r.dx - (double) dx_int);
			}
			/* On failure we leave the cumulative shift chain frozen
			 * at the last successful step and store that. */
			out.shifts[idx] = cv::Vec2d(dy_store, dx_store);
			const gint d = g_atomic_int_add(&done, 1) + 1;
			if (progress) progress((double) d / (double) total, progress_user);
		}
		} catch (const std::exception &) {
			g_atomic_int_set(&oom, 1);
		}
	};
	auto run_forward_sweep = [&] {
		try {
		int dy_int_cum = 0, dx_int_cum = 0;
		for (int idx = best + 1; idx < N; ++idx) {
			if (g_atomic_int_get(&cancelled)) break;
			if (!processing_should_continue()) {
				g_atomic_int_set(&cancelled, 1);
				break;
			}
			const cv::Mat frame = provider(idx);
			if (use_seed) {
				/* Independent per-frame: window centred at this frame's
				 * seed (re-referenced to best); total = seed + residual. */
				const int sdy = (int) std::lround(seed[idx][0] - seed_base[0]);
				const int sdx = (int) std::lround(seed[idx][1] - seed_base[1]);
				const AlignShiftResult rs = align_shift_one_frame(
				    ref_window, ref_window_first, frame,
				    py_lo - sdy, py_hi - sdy, px_lo - sdx, px_hi - sdx, cfg);
				double dy_s = (double) sdy, dx_s = (double) sdx;
				if (rs.success) { dy_s = (double) sdy + rs.dy; dx_s = (double) sdx + rs.dx; }
				out.shifts[idx] = cv::Vec2d(dy_s, dx_s);
				const gint d = g_atomic_int_add(&done, 1) + 1;
				if (progress) progress((double) d / (double) total, progress_user);
				continue;
			}
			const AlignShiftResult r = align_shift_one_frame(
			    ref_window, ref_window_first, frame,
			    py_lo - dy_int_cum, py_hi - dy_int_cum,
			    px_lo - dx_int_cum, px_hi - dx_int_cum, cfg);
			double dy_store = (double) dy_int_cum;
			double dx_store = (double) dx_int_cum;
			if (r.success) {
				const int dy_int = (int) std::lround(r.dy);
				const int dx_int = (int) std::lround(r.dx);
				dy_int_cum += dy_int;
				dx_int_cum += dx_int;
				dy_store = (double) dy_int_cum + (r.dy - (double) dy_int);
				dx_store = (double) dx_int_cum + (r.dx - (double) dx_int);
			}
			out.shifts[idx] = cv::Vec2d(dy_store, dx_store);
			const gint d = g_atomic_int_add(&done, 1) + 1;
			if (progress) progress((double) d / (double) total, progress_user);
		}
		} catch (const std::exception &) {
			g_atomic_int_set(&oom, 1);
		}
	};

#ifdef _OPENMP
#pragma omp parallel sections num_threads(2) if (provider_is_thread_safe)
	{
#pragma omp section
		run_backward_sweep();
#pragma omp section
		run_forward_sweep();
	}
#else
	(void) provider_is_thread_safe;
	run_backward_sweep();
	run_forward_sweep();
#endif

	if (g_atomic_int_get(&oom)) {
		out.best_frame_idx = -1;
		out.oom = true;
		return out;
	}
	if (g_atomic_int_get(&cancelled)) {
		out.best_frame_idx = -1;   /* cancellation sentinel */
		return out;
	}
	return out;
}

AlignGlobalResult align_global_from_frames(const std::vector<cv::Mat> &frames,
                                           const std::vector<double> &quality,
                                           const mpp_config_t &cfg,
                                           progress_cb_fn progress,
                                           void *progress_user,
                                           const std::vector<cv::Vec2d> &seed) {
	/* Reading from a const cached vector is trivially reentrant — enable
	 * the two-section parallel sweep. */
	return align_global_from_provider(
	    [&frames](int i) -> cv::Mat { return frames[i]; },
	    (int) frames.size(), quality, cfg, progress, progress_user,
	    /*provider_is_thread_safe=*/true, seed);
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
                                                const std::vector<cv::Vec2d> &shifts,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress,
                                                void *progress_user,
                                                const std::vector<int> &included) {
	AlignAverageResult out;
	if (N == 0 || (int) shifts.size() != N || (int) quality.size() != N)
		return out;

	const bool have_mask = (int) included.size() == N;
	auto is_included = [&](int i) { return !have_mask || included[i]; };
	int n_included = 0;
	for (int i = 0; i < N; ++i) if (is_included(i)) ++n_included;
	if (n_included == 0) return out;   /* nothing to build a reference from */

	/* Intersection bounds are integer — they're used as cv::Range slice
	 * indices everywhere downstream. Round sub-pixel shifts to nearest
	 * integer for the min/max; the 0.5-pixel rounding doesn't affect
	 * the outer bounds materially. Excluded frames are skipped: a glitch
	 * frame with a wild global shift would otherwise collapse the common
	 * overlap area to nothing. */
	int max_dy = INT_MIN, min_dy = INT_MAX, max_dx = INT_MIN, min_dx = INT_MAX;
	for (int i = 0; i < N; ++i) {
		if (!is_included(i)) continue;
		const int sy = (int) std::lround(shifts[i][0]);
		const int sx = (int) std::lround(shifts[i][1]);
		max_dy = std::max(max_dy, sy);
		min_dy = std::min(min_dy, sy);
		max_dx = std::max(max_dx, sx);
		min_dx = std::min(min_dx, sx);
	}
	out.intersection = cv::Vec4i(max_dy, min_dy + H, max_dx, min_dx + W);
	const int int_y_lo = out.intersection[0], int_y_hi = out.intersection[1];
	const int int_x_lo = out.intersection[2], int_x_hi = out.intersection[3];
	if (int_y_hi <= int_y_lo || int_x_hi <= int_x_lo)
		return out;

	/* Clamp to the included count so the masked quality (excluded frames
	 * pre-set to a sentinel by the caller) can never pull an excluded
	 * frame into the averaged reference. */
	const int n_target = std::min(n_included, std::max(
	    (int) std::ceil(N * cfg.align_frames_average_frame_percent / 100.0), 1));

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
		if (!processing_should_continue()) {
			out.mean_frame.release();   /* cancellation sentinel */
			return out;
		}
		const auto &s = shifts[idx];
		const int sy = (int) std::lround(s[0]);
		const int sx = (int) std::lround(s[1]);
		const cv::Mat frame = provider(idx);
		if (frame.empty()) { ++done; continue; }
		const cv::Mat patch = frame(
		    cv::Range(int_y_lo - sy, int_y_hi - sy),
		    cv::Range(int_x_lo - sx, int_x_hi - sx));
		cv::Mat patch_f32;
		patch.convertTo(patch_f32, CV_32F);
		accum += patch_f32;
		++done;
		if (progress) progress((double) done / (double) total, progress_user);
	}

	/* Scale by 256/N for 8-bit inputs and 1/N for 16-bit. Siril stores
	 * 8-bit SER data as WORD without upscaling (values 0..255), so we
	 * look at cfg.bitdepth — NOT cv::Mat::depth() — to decide. Output
	 * mean_frame always lands in the 0..65535 range, so downstream AP
	 * code can use a single `× 256` threshold scale. */
	const double scale = (cfg.bitdepth == 8 ? 256.0 : 1.0) / (double) indices.size();
	accum *= scale;
	/* Use NumPy-style `.astype(int32)` truncation toward zero. OpenCV's
	 * convertTo(CV_32S) uses cvRound (round half to even on x86), so a
	 * naive cast would produce off-by-one diffs at half-pixel values and
	 * those propagate into the per-AP correlation peaks downstream. Do
	 * the truncation manually to keep bit-equivalence with the reference
	 * algorithm. */
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
                                       const std::vector<cv::Vec2d> &shifts,
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

extern "C" gboolean mpp_frame_has_disc(const fits *fit) {
	cv::Mat layer = mpp::wrap_fits_layer(fit, mpp::analysis_layer_for(fit));
	if (layer.empty() || layer.rows < 16 || layer.cols < 16)
		return FALSE;

	/* Bin down to <= 256 px on the long side with area averaging. This
	 * suppresses hot pixels and noise, and averages over CFA mosaics so
	 * the test also works on raw Bayer frames. */
	cv::Mat small_f;
	layer.convertTo(small_f, CV_32F);
	const int maxdim = std::max(small_f.rows, small_f.cols);
	if (maxdim > 256) {
		const double scale = 256.0 / maxdim;
		cv::resize(small_f, small_f, cv::Size(), scale, scale, cv::INTER_AREA);
	}
	cv::GaussianBlur(small_f, small_f, cv::Size(3, 3), 0.0);

	/* Mid-range threshold, same convention as center_of_gravity(). */
	double min_val = 0.0, max_val = 0.0;
	cv::minMaxLoc(small_f, &min_val, &max_val);
	if (max_val <= min_val)
		return FALSE;	/* flat frame, nothing to detect */
	const cv::Mat mask = small_f > (min_val + max_val) / 2.0;

	/* The object must have some extent: a handful of bright pixels is a
	 * star (or noise), not a disc. */
	const int npix = small_f.rows * small_f.cols;
	const int bright = cv::countNonZero(mask);
	if (bright < std::max(9, npix / 1000))
		return FALSE;

	/* A disc — round or elongated like Saturn — is a bright object
	 * completely surrounded by sky, so (almost) no bright pixel may lie
	 * in a thin band along the frame edges. A surface close-up always
	 * reaches at least one edge. The 0.5 % tolerance forgives residual
	 * noise in the band without letting a real limb through. */
	const int band = std::max(2, std::min(small_f.rows, small_f.cols) / 50);
	const cv::Rect inner_rect(band, band,
	                          small_f.cols - 2 * band, small_f.rows - 2 * band);
	const int border_bright = bright - cv::countNonZero(mask(inner_rect));
	const int border_pix = npix - inner_rect.width * inner_rect.height;
	return border_bright <= border_pix / 200;
}
