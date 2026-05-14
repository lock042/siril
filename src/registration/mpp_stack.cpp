/*
 * Phase 5a: PSS-faithful multipoint stacking.
 *
 * Builds the final stacked image by:
 *   1) Compute per-AP per-frame qualities; keep top-N per AP (PSS
 *      alignment_points.compute_frame_qualities).
 *   2) Pre-compute per-AP weights_yx and the global sum_single_frame_weights
 *      buffer (prepare_for_stack_blending).
 *   3) For each frame: brightness-equalise to the per-frame brightness
 *      median, resize for drizzle (INTER_LINEAR), recompute per-AP shift
 *      with the weight_matrix_first_phase penalty, then remap_rigid the
 *      patch into the AP's stacking buffer.
 *   4) Merge AP buffers into the global stacking buffer with the weights
 *      from (2), blend in the background where APs don't cover, trim
 *      shift-induced borders, scale + clip to uint16.
 *
 * Each piece lives behind a small C++ entry point in mpp_stack_priv.hpp
 * with its own oracle equivalence test. See pss_port_plan.md §5a.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp_align_priv.hpp"     /* multilevel_correlation */
#include "registration/mpp_rank_priv.hpp"
#include "registration/mpp_shift_priv.hpp"     /* shift_prepare_ref_boxes */
#include "registration/mpp_stack.h"
#include "registration/mpp_stack_priv.hpp"

#include <set>

namespace mpp {

std::vector<float> stack_one_dim_weight(int patch_low, int patch_high, int box_center,
                                        bool extend_low, bool extend_high) {
	const int patch_size    = patch_high - patch_low;
	const int center_offset = box_center - patch_low;
	std::vector<float> w(patch_size, 0.0f);

	if (extend_low) {
		for (int i = 0; i < center_offset; ++i) w[i] = 1.0f;
	} else {
		/* PSS: arange(1, c+1) / float32(c+1) → [1/(c+1), 2/(c+1), …, c/(c+1)]. */
		const float denom = (float) (center_offset + 1);
		for (int i = 0; i < center_offset; ++i)
			w[i] = (float) (i + 1) / denom;
	}

	const int high_len = patch_size - center_offset;  /* = patch_high - box_center */
	if (extend_high) {
		for (int i = 0; i < high_len; ++i) w[center_offset + i] = 1.0f;
	} else {
		/* PSS: arange(N, 0, -1) / float32(N) → [N/N, (N-1)/N, …, 1/N]. */
		const float denom = (float) high_len;
		for (int i = 0; i < high_len; ++i)
			w[center_offset + i] = (float) (high_len - i) / denom;
	}
	return w;
}

cv::Mat stack_build_first_phase_weight_matrix(const mpp_config_t &cfg) {
	const int sw2 = 4;
	const int sw1 = (cfg.alignment_points_search_width - sw2) / 2;
	const int extent = 2 * sw1 + 1;
	cv::Mat m(extent, extent, CV_32F);
	const double pen = cfg.alignment_points_penalty_factor;
	for (int y = 0; y < extent; ++y) {
		float *row = m.ptr<float>(y);
		const double ty = (double) y / (double) sw1 - 1.0;
		for (int x = 0; x < extent; ++x) {
			const double tx = (double) x / (double) sw1 - 1.0;
			row[x] = (float) (1.0 - pen * (ty * ty + tx * tx));
		}
	}
	return m;
}

void stack_remap_rigid(const cv::Mat &frame, cv::Mat &buffer,
                       int shift_y, int shift_x,
                       int y_low, int y_high, int x_low, int x_high,
                       RemapBorder &border) {
	/* Mirrors PSS stack_frames.remap_rigid. Tracks the maximum clipped
	 * extent on each border so the final trim step can remove artifacts. */
	const int frame_h = frame.rows, frame_w = frame.cols;
	int y_low_src  = y_low  + shift_y;
	int y_high_src = y_high + shift_y;
	int y_low_dst  = 0;
	if (y_low_src < 0) {
		y_low_dst = -y_low_src;
		y_low_src = 0;
		border.y_low = std::max(border.y_low, y_low_dst);
	}
	if (y_high_src > frame_h) {
		border.y_high = std::max(border.y_high, y_high_src - frame_h);
		y_high_src = frame_h;
	}
	const int y_high_dst = y_low_dst + (y_high_src - y_low_src);

	int x_low_src  = x_low  + shift_x;
	int x_high_src = x_high + shift_x;
	int x_low_dst  = 0;
	if (x_low_src < 0) {
		x_low_dst = -x_low_src;
		x_low_src = 0;
		border.x_low = std::max(border.x_low, x_low_dst);
	}
	if (x_high_src > frame_w) {
		border.x_high = std::max(border.x_high, x_high_src - frame_w);
		x_high_src = frame_w;
	}
	const int x_high_dst = x_low_dst + (x_high_src - x_low_src);

	if (y_high_dst <= y_low_dst || x_high_dst <= x_low_dst)
		return;
	const cv::Mat src = frame(cv::Range(y_low_src, y_high_src),
	                          cv::Range(x_low_src, x_high_src));
	cv::Mat dst = buffer(cv::Range(y_low_dst, y_high_dst),
	                     cv::Range(x_low_dst, x_high_dst));
	dst += src;
}

APQualities ap_compute_frame_qualities(const std::vector<cv::Mat> &frames,
                                       const std::vector<double> &frame_brightness,
                                       const mpp_aps_t &aps,
                                       const std::vector<FrameOffset> &offsets,
                                       int frame_rows, int frame_cols,
                                       const mpp_config_t &cfg) {
	APQualities out;
	const int N = (int) frames.size();
	const int M = aps.count;
	const int stride = cfg.align_frames_sampling_stride;
	const bool normalise = cfg.frames_normalization;
	if (N == 0 || M == 0 || (int) frame_brightness.size() != N
	 || (int) offsets.size() != N || stride <= 0)
		return out;

	/* PSS stack_size: explicit override, else %% of N (PSS rounds up). */
	if (cfg.alignment_points_frame_number > 0)
		out.stack_size = cfg.alignment_points_frame_number;
	else
		out.stack_size = std::max(
		    (int) std::ceil((double) N * cfg.alignment_points_frame_percent / 100.0),
		    1);
	out.stack_size = std::min(out.stack_size, N);

	out.qualities.assign(M, std::vector<double>(N, 0.0));

	for (int f = 0; f < N; ++f) {
		const cv::Mat lap = rank_blurred_laplacian_u8(frames[f], cfg);
		const int dy = offsets[f].dy, dx = offsets[f].dx;
		for (int a = 0; a < M; ++a) {
			const auto &ap = aps.records[a];
			/* PSS:
			 *   y_low  = int(max(0, patch_y_low  + dy) / stride);
			 *   y_high = int(min(frame_h, patch_y_high + dy) / stride);
			 * Python `int()` truncates toward zero. C integer division
			 * does the same for non-negatives, so / stride matches. */
			const int yl_raw = ap.patch_y_low  + dy;
			const int yh_raw = ap.patch_y_high + dy;
			const int xl_raw = ap.patch_x_low  + dx;
			const int xh_raw = ap.patch_x_high + dx;
			const int yl = std::max(0, yl_raw) / stride;
			const int yh = std::min(frame_rows, yh_raw) / stride;
			const int xl = std::max(0, xl_raw) / stride;
			const int xh = std::min(frame_cols, xh_raw) / stride;
			double sigma = 0.0;
			if (yh > yl && xh > xl
			 && yh <= lap.rows && xh <= lap.cols) {
				const cv::Mat box = lap(cv::Range(yl, yh), cv::Range(xl, xh));
				cv::Scalar mean, stddev;
				cv::meanStdDev(box, mean, stddev);
				sigma = stddev[0];
			}
			out.qualities[a][f] = normalise ? sigma / frame_brightness[f] : sigma;
		}
	}

	/* Per-AP: pick top-N frame indices (PSS sorts in descending quality
	 * order and slices [:stack_size]). */
	out.best_frame_indices.assign(M, {});
	out.used_alignment_points.assign(N, {});
	for (int a = 0; a < M; ++a) {
		std::vector<int> indices(N);
		std::iota(indices.begin(), indices.end(), 0);
		const auto &q = out.qualities[a];
		std::partial_sort(indices.begin(), indices.begin() + out.stack_size,
		                  indices.end(),
		                  [&q](int x, int y) { return q[x] > q[y]; });
		out.best_frame_indices[a].assign(indices.begin(),
		                                 indices.begin() + out.stack_size);
		for (int i = 0; i < out.stack_size; ++i)
			out.used_alignment_points[indices[i]].push_back(a);
	}
	return out;
}

StackState stack_prepare_for_blending(const mpp_aps_t &aps,
                                      const cv::Vec4i &intersection,
                                      int stack_size,
                                      int drizzle_factor,
                                      const mpp_config_t &cfg) {
	StackState s;
	s.drizzle_factor = drizzle_factor;
	s.stack_size     = stack_size;
	s.dim_y = intersection[1] - intersection[0];
	s.dim_x = intersection[3] - intersection[2];
	s.dim_y_drizzled = s.dim_y * drizzle_factor;
	s.dim_x_drizzled = s.dim_x * drizzle_factor;

	const int M = aps.count;
	s.patch_drizzled.reserve(M);
	s.ap_drizzled.reserve(M);
	s.weights_yx.reserve(M);
	s.stacking_buffers.reserve(M);
	s.sum_single_frame_weights.create(s.dim_y_drizzled, s.dim_x_drizzled, CV_32F);
	s.sum_single_frame_weights.setTo(cv::Scalar(1e-30f));

	const float stack_size_f = (float) stack_size;

	for (int a = 0; a < M; ++a) {
		const auto &ap = aps.records[a];
		const int py_lo = ap.patch_y_low  * drizzle_factor;
		const int py_hi = ap.patch_y_high * drizzle_factor;
		const int px_lo = ap.patch_x_low  * drizzle_factor;
		const int px_hi = ap.patch_x_high * drizzle_factor;
		const int yc    = ap.y * drizzle_factor;
		const int xc    = ap.x * drizzle_factor;

		s.patch_drizzled.push_back(cv::Vec4i(py_lo, py_hi, px_lo, px_hi));
		s.ap_drizzled.push_back(cv::Vec2i(yc, xc));

		const bool extend_y_low  = (py_lo == 0);
		const bool extend_y_high = (py_hi == s.dim_y_drizzled);
		const bool extend_x_low  = (px_lo == 0);
		const bool extend_x_high = (px_hi == s.dim_x_drizzled);
		const auto wy = stack_one_dim_weight(py_lo, py_hi, yc, extend_y_low,  extend_y_high);
		const auto wx = stack_one_dim_weight(px_lo, px_hi, xc, extend_x_low,  extend_x_high);

		cv::Mat wyx((int) wy.size(), (int) wx.size(), CV_32F);
		for (int y = 0; y < (int) wy.size(); ++y) {
			float *row = wyx.ptr<float>(y);
			for (int x = 0; x < (int) wx.size(); ++x)
				row[x] = std::min(wy[y], wx[x]);
		}
		s.weights_yx.push_back(wyx);
		s.stacking_buffers.push_back(cv::Mat(wyx.rows, wyx.cols, CV_32F,
		                                     cv::Scalar(0)));

		/* Accumulate stack_size × weights_yx into the global buffer. */
		cv::Mat dst = s.sum_single_frame_weights(
		    cv::Range(py_lo, py_hi), cv::Range(px_lo, px_hi));
		dst += stack_size_f * wyx;
	}

	cv::Mat hole_mask;
	cv::compare(s.sum_single_frame_weights, 1e-10, hole_mask, cv::CMP_LT);
	s.number_stacking_holes = cv::countNonZero(hole_mask);
	if (s.number_stacking_holes == 0)
		return s;

	/* Background buffer (pre-drizzle, mono). */
	s.averaged_background = cv::Mat::zeros(s.dim_y, s.dim_x, CV_32F);

	/* PSS stack_frames.py:233-277. Count pixels where sum_weights <
	 * background_blend_threshold * stack_size; if a small fraction, subdivide
	 * the intersection into background_patch_size patches and keep only those
	 * patches that contain at least one such pixel. */
	const float blend_t = (float) cfg.stack_frames_background_blend_threshold;
	cv::Mat needed_mask;
	cv::compare(s.sum_single_frame_weights, blend_t * stack_size_f, needed_mask, cv::CMP_LT);
	const int points_where_used = cv::countNonZero(needed_mask);
	const double total_drizzled = (double) s.dim_y_drizzled * (double) s.dim_x_drizzled;
	if ((double) points_where_used / total_drizzled
	    < cfg.stack_frames_background_fraction) {
		const int patch_size = cfg.stack_frames_background_patch_size;
		for (int py = 0; py < s.dim_y; py += patch_size) {
			/* PSS uses `dim_y - 1` as the cap — replicate. */
			const int py_hi = std::min(py + patch_size, s.dim_y - 1);
			if (py == py_hi) continue;
			for (int px = 0; px < s.dim_x; px += patch_size) {
				const int px_hi = std::min(px + patch_size, s.dim_x - 1);
				if (px == px_hi) continue;
				const cv::Mat slice = needed_mask(
				    cv::Range(py * drizzle_factor, py_hi * drizzle_factor),
				    cv::Range(px * drizzle_factor, px_hi * drizzle_factor));
				if (cv::countNonZero(slice) > 0)
					s.background_patches.push_back({py, py_hi, px, px_hi});
			}
		}
	}
	return s;
}

mpp_shifts_t *stack_compute_shifts(const std::vector<cv::Mat> &frames_mono_blurred,
                                   const cv::Mat &mean_frame_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const std::vector<FrameOffset> &offsets,
                                   const mpp_config_t &cfg) {
	const int N = (int) frames_mono_blurred.size();
	const int M = aps.count;
	const int K = cfg.drizzle_factor;
	const bool use_subpixel = (K != 1);

	mpp_shifts_t *out = (mpp_shifts_t *) std::calloc(1, sizeof(*out));
	if (!out) return nullptr;
	out->num_frames = N;
	out->num_aps = M;
	out->shifts  = (double *) std::calloc((size_t) N * (size_t) M * 2, sizeof(double));
	out->success = (uint8_t *) std::calloc((size_t) N * (size_t) M, sizeof(uint8_t));
	if (!out->shifts || !out->success) {
		mpp_shift_free(out);
		return nullptr;
	}

	const cv::Mat weight_matrix = stack_build_first_phase_weight_matrix(cfg);
	const auto ref_boxes = shift_prepare_ref_boxes(mean_frame_raw, aps);
	int failures = 0;

	for (int f = 0; f < N; ++f) {
		const int dy = offsets[f].dy, dx = offsets[f].dx;
		for (int a : apq.used_alignment_points[f]) {
			const auto &ap = aps.records[a];
			const MultilevelShiftResult r = multilevel_correlation(
			    ref_boxes[a].second_phase, ref_boxes[a].first_phase,
			    frames_mono_blurred[f],
			    ap.box_y_low + dy, ap.box_y_high + dy,
			    ap.box_x_low + dx, ap.box_x_high + dx,
			    cfg.frames_gauss_width, cfg.alignment_points_search_width,
			    use_subpixel, weight_matrix);
			const size_t off = (size_t) (f * M + a) * 2;
			out->shifts[off + 0] = r.dy;
			out->shifts[off + 1] = r.dx;
			out->success[f * M + a] = r.success ? 1 : 0;
			if (!r.success) ++failures;
		}
	}
	out->failure_counter = failures;
	return out;
}

StackLoopOutput stack_apply_shifts(const std::vector<cv::Mat> &frames_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const mpp_shifts_t *shifts,
                                   const std::vector<FrameOffset> &offsets,
                                   const std::vector<double> &frame_brightness,
                                   const std::vector<int> &quality_sorted_idx,
                                   const cv::Vec4i &intersection,
                                   const mpp_config_t &cfg) {
	StackLoopOutput out;
	out.state = stack_prepare_for_blending(aps, intersection, apq.stack_size,
	                                       cfg.drizzle_factor, cfg);
	out.shift_failure_counter = shifts ? shifts->failure_counter : 0;

	const int N = (int) frames_raw.size();
	const int M = aps.count;
	const int K = cfg.drizzle_factor;

	double median_brightness = 0.0;
	if (cfg.frames_normalization) {
		std::vector<double> sorted_b(frame_brightness);
		std::sort(sorted_b.begin(), sorted_b.end());
		const size_t n = sorted_b.size();
		median_brightness = (n % 2 == 1)
		    ? sorted_b[n / 2]
		    : 0.5 * (sorted_b[n / 2 - 1] + sorted_b[n / 2]);
	}

	std::set<int> top_for_bg;
	for (int i = 0; i < apq.stack_size && i < (int) quality_sorted_idx.size(); ++i)
		top_for_bg.insert(quality_sorted_idx[i]);

	for (int f = 0; f < N; ++f) {
		const cv::Mat &frame_raw = frames_raw[f];

		cv::Mat frame_f32;
		if (cfg.frames_normalization) {
			const double scale = median_brightness
			                   / (frame_brightness[f] + 1e-7);
			frame_raw.convertTo(frame_f32, CV_32F, scale);
		} else {
			frame_raw.convertTo(frame_f32, CV_32F);
		}

		cv::Mat frame_drizzled;
		if (K != 1)
			cv::resize(frame_f32, frame_drizzled,
			           cv::Size(frame_f32.cols * K, frame_f32.rows * K),
			           0, 0, cv::INTER_LINEAR);
		else
			frame_drizzled = frame_f32;

		const int dy = offsets[f].dy;
		const int dx = offsets[f].dx;
		const int dy_K = dy * K;
		const int dx_K = dx * K;

		for (int a : apq.used_alignment_points[f]) {
			const size_t soff = (size_t) (f * M + a) * 2;
			const double shift_y = shifts ? shifts->shifts[soff + 0] : 0.0;
			const double shift_x = shifts ? shifts->shifts[soff + 1] : 0.0;
			const int shift_y_K = (int) std::lround(shift_y * (double) K);
			const int shift_x_K = (int) std::lround(shift_x * (double) K);
			const int total_y_K = dy_K - shift_y_K;
			const int total_x_K = dx_K - shift_x_K;

			const auto &p = out.state.patch_drizzled[a];
			stack_remap_rigid(frame_drizzled, out.state.stacking_buffers[a],
			                  total_y_K, total_x_K, p[0], p[1], p[2], p[3],
			                  out.border);
		}

		if (out.state.number_stacking_holes > 0
		 && top_for_bg.find(f) != top_for_bg.end()) {
			if (!out.state.background_patches.empty()) {
				for (const auto &bp : out.state.background_patches) {
					const cv::Mat src = frame_f32(
					    cv::Range(bp.y_low + dy, bp.y_high + dy),
					    cv::Range(bp.x_low + dx, bp.x_high + dx));
					cv::Mat dst = out.state.averaged_background(
					    cv::Range(bp.y_low, bp.y_high),
					    cv::Range(bp.x_low, bp.x_high));
					dst += src;
				}
			} else {
				const cv::Mat src = frame_f32(
				    cv::Range(dy, out.state.dim_y + dy),
				    cv::Range(dx, out.state.dim_x + dx));
				out.state.averaged_background += src;
			}
		}
	}

	if (out.state.number_stacking_holes > 0) {
		if (K != 1) {
			cv::Mat r;
			cv::resize(out.state.averaged_background, r,
			           cv::Size(out.state.dim_x_drizzled, out.state.dim_y_drizzled),
			           0, 0, cv::INTER_CUBIC);
			out.state.averaged_background = r;
		}
		out.state.averaged_background /= (float) apq.stack_size;
	}

	return out;
}

StackLoopOutput stack_frames_loop(const std::vector<cv::Mat> &frames_raw,
                                  const std::vector<cv::Mat> &frames_mono_blurred,
                                  const cv::Mat &mean_frame_raw,
                                  const mpp_aps_t &aps,
                                  const APQualities &apq,
                                  const std::vector<FrameOffset> &offsets,
                                  const std::vector<double> &frame_brightness,
                                  const std::vector<int> &quality_sorted_idx,
                                  const cv::Vec4i &intersection,
                                  const mpp_config_t &cfg) {
	/* Thin wrapper preserving the Phase-5a API: Stage B then Stage C. */
	mpp_shifts_t *shifts = stack_compute_shifts(frames_mono_blurred,
	                                            mean_frame_raw, aps, apq,
	                                            offsets, cfg);
	StackLoopOutput out = stack_apply_shifts(frames_raw, aps, apq, shifts,
	                                         offsets, frame_brightness,
	                                         quality_sorted_idx, intersection, cfg);
	mpp_shift_free(shifts);
	return out;
}

cv::Mat stack_merge_alignment_point_buffers(const StackState &state,
                                            const RemapBorder &border,
                                            const mpp_aps_t &aps,
                                            const mpp_config_t &cfg) {
	/* Global drizzled image buffer (float32). */
	cv::Mat buf = cv::Mat::zeros(state.dim_y_drizzled, state.dim_x_drizzled, CV_32F);

	/* For each AP: buf[patch] += stacking_buffer * weights_yx. */
	for (int a = 0; a < aps.count; ++a) {
		const auto &p = state.patch_drizzled[a];
		cv::Mat dst = buf(cv::Range(p[0], p[1]), cv::Range(p[2], p[3]));
		cv::Mat weighted;
		cv::multiply(state.stacking_buffers[a], state.weights_yx[a], weighted);
		dst += weighted;
	}

	/* buf /= sum_single_frame_weights (cell-wise). PSS guards by initialising
	 * sum_weights to 1e-30, so we never divide by zero. */
	cv::divide(buf, state.sum_single_frame_weights, buf);

	/* Background blend. PSS computes a foreground_weight ramp and lerps. */
	if (state.number_stacking_holes > 0) {
		const float denom = (float) (cfg.stack_frames_background_blend_threshold
		                            * (double) state.stack_size);
		cv::Mat fg_weight;
		state.sum_single_frame_weights.convertTo(fg_weight, CV_32F, 1.0 / denom);
		cv::min(fg_weight, 1.0f, fg_weight);
		cv::max(fg_weight, 0.0f, fg_weight);
		/* buf = (buf - bg) * fg_weight + bg = bg + (buf - bg) * fg_weight. */
		cv::Mat diff;
		cv::subtract(buf, state.averaged_background, diff);
		cv::multiply(diff, fg_weight, diff);
		cv::add(state.averaged_background, diff, buf);
	}

	/* Border trim. */
	const int trim_y_lo = border.y_low;
	const int trim_y_hi = state.dim_y_drizzled - border.y_high;
	const int trim_x_lo = border.x_low;
	const int trim_x_hi = state.dim_x_drizzled - border.x_high;
	if (border.y_low || border.y_high || border.x_low || border.x_high)
		buf = buf(cv::Range(trim_y_lo, trim_y_hi),
		          cv::Range(trim_x_lo, trim_x_hi)).clone();

	/* Scale and cast to uint16.
	 * PSS final cast: `img_as_uint(clip(buf / (255 or 65535), 0., 1., ...))`.
	 * skimage.img_as_uint on a float in [0, 1] uses round + saturate_cast<u16>.
	 * For Siril values that landed in 0..65535 the scale is 1/65535; the
	 * round-trip is therefore approximately identity (modulo rounding). */
	const double in_max = (cfg.bitdepth == 8) ? 255.0 : 65535.0;
	cv::Mat normalised;
	buf.convertTo(normalised, CV_32F, 1.0 / in_max);
	cv::min(normalised, 1.0f, normalised);
	cv::max(normalised, 0.0f, normalised);

	cv::Mat out(normalised.rows, normalised.cols, CV_16U);
	/* skimage img_as_uint does floor(x * 65535 + 0.5).  Replicate exactly to
	 * dodge OpenCV's round-half-to-even at .5. */
	for (int y = 0; y < normalised.rows; ++y) {
		const float *src = normalised.ptr<float>(y);
		uint16_t *dst = out.ptr<uint16_t>(y);
		for (int x = 0; x < normalised.cols; ++x) {
			const double v = (double) src[x] * 65535.0 + 0.5;
			const int iv = (int) v;
			dst[x] = (uint16_t) std::clamp(iv, 0, 65535);
		}
	}
	return out;
}

}  // namespace mpp

/* ----------------- Public C interface (still stubbed) ----------------- */

extern "C" mpp_status_t mpp_stack(sequence *seq,
                                  const mpp_config_t *cfg,
                                  const mpp_aps_t *aps,
                                  const mpp_shifts_t *shifts,
                                  const int *global_shifts,
                                  mpp_resample_kind_t backend, int upscale,
                                  fits *stacked_out) {
	(void) seq; (void) cfg; (void) aps; (void) shifts; (void) global_shifts;
	(void) backend; (void) upscale; (void) stacked_out;
	return MPP_ENOTIMPL;
}
