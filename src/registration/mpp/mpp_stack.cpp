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
 * Phase 5a: multipoint stacking.
 *
 * Builds the final stacked image by:
 *   1) Compute per-AP per-frame qualities; keep top-N per AP.
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
 * with its own oracle equivalence test.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp/mpp_align_priv.hpp"     /* multilevel_correlation */
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_shift_priv.hpp"     /* shift_prepare_ref_boxes */
#include "registration/mpp/mpp_stack.h"
#include "registration/mpp/mpp_stack_priv.hpp"

#include "core/gui_iface.h"                    /* set_progress (no-op stub in headless) */
#include "core/processing.h"                   /* processing_should_continue */

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
		/* arange(1, c+1) / float32(c+1) → [1/(c+1), 2/(c+1), …, c/(c+1)]. */
		const float denom = (float) (center_offset + 1);
		for (int i = 0; i < center_offset; ++i)
			w[i] = (float) (i + 1) / denom;
	}

	const int high_len = patch_size - center_offset;  /* = patch_high - box_center */
	if (extend_high) {
		for (int i = 0; i < high_len; ++i) w[center_offset + i] = 1.0f;
	} else {
		/* arange(N, 0, -1) / float32(N) → [N/N, (N-1)/N, …, 1/N]. */
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
	/* Rigid remap of the AP's patch into its stacking buffer. Tracks the
	 * maximum clipped extent on each border so the final trim step can
	 * remove the resulting artefact ring. */
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

APQualities ap_compute_frame_qualities_streamed(const FrameProvider &provider,
                                                int N,
                                                const std::vector<double> &frame_brightness,
                                                const mpp_aps_t &aps,
                                                const std::vector<FrameOffset> &offsets,
                                                int frame_rows, int frame_cols,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress,
                                                void *progress_user) {
	APQualities out;
	const int M = aps.count;
	const int stride = cfg.align_frames_sampling_stride;
	const bool normalise = cfg.frames_normalization;
	if (N == 0 || M == 0 || (int) frame_brightness.size() != N
	 || (int) offsets.size() != N || stride <= 0)
		return out;

	/* stack_size: explicit override, else %% of N (rounded up). */
	if (cfg.alignment_points_frame_number > 0)
		out.stack_size = cfg.alignment_points_frame_number;
	else
		out.stack_size = std::max(
		    (int) std::ceil((double) N * cfg.alignment_points_frame_percent / 100.0),
		    1);
	out.stack_size = std::min(out.stack_size, N);

	out.qualities.assign(M, std::vector<double>(N, 0.0));

	for (int f = 0; f < N; ++f) {
		if (!processing_should_continue()) {
			out.stack_size = 0;   /* cancellation sentinel */
			return out;
		}
		const cv::Mat frame = provider(f);
		if (frame.empty()) {
			if (progress) progress((double)(f + 1) / (double) N, progress_user);
			continue;
		}
		const cv::Mat lap = rank_blurred_laplacian_u8(frame, cfg);
		/* Round sub-pixel offsets — Stage A's box-slicing math wants
		 * integer indices; the sub-pixel residual matters only at the
		 * drizzle pixmap stage. */
		const int dy = (int) std::lround(offsets[f].dy);
		const int dx = (int) std::lround(offsets[f].dx);
		if (progress) progress((double)(f + 1) / (double) N, progress_user);
		for (int a = 0; a < M; ++a) {
			const auto &ap = aps.records[a];
			/* y_low  = int(max(0, patch_y_low  + dy) / stride);
			 * y_high = int(min(frame_h, patch_y_high + dy) / stride);
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

	/* Per-AP: pick top-N frame indices (descending quality, sliced to
	 * [:stack_size]). */
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

APQualities ap_compute_frame_qualities(const std::vector<cv::Mat> &frames,
                                       const std::vector<double> &frame_brightness,
                                       const mpp_aps_t &aps,
                                       const std::vector<FrameOffset> &offsets,
                                       int frame_rows, int frame_cols,
                                       const mpp_config_t &cfg,
                                       progress_cb_fn progress,
                                       void *progress_user) {
	return ap_compute_frame_qualities_streamed(
	    [&frames](int i) { return frames[i]; },
	    (int) frames.size(), frame_brightness, aps, offsets,
	    frame_rows, frame_cols, cfg, progress, progress_user);
}

StackState stack_prepare_for_blending(const mpp_aps_t &aps,
                                      const cv::Vec4i &intersection,
                                      int stack_size,
                                      double drizzle_scale,
                                      int num_layers,
                                      const mpp_config_t &cfg) {
	StackState s;
	s.drizzle_scale  = drizzle_scale;
	s.stack_size     = stack_size;
	s.num_layers     = num_layers;
	s.dim_y = intersection[1] - intersection[0];
	s.dim_x = intersection[3] - intersection[2];
	/* Output (drizzled) dimensions. drizzle_scale may be fractional (e.g.
	 * 1.5); every coordinate below is the original × scale rounded to the
	 * nearest pixel. lround is monotonic so patch bounds never exceed the
	 * canvas. */
	auto scl = [drizzle_scale](int v) { return (int) std::lround((double) v * drizzle_scale); };
	s.dim_y_drizzled = scl(s.dim_y);
	s.dim_x_drizzled = scl(s.dim_x);
	const int buf_type = (num_layers == 3) ? CV_32FC3 : CV_32F;

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
		const int py_lo = scl(ap.patch_y_low);
		const int py_hi = scl(ap.patch_y_high);
		const int px_lo = scl(ap.patch_x_low);
		const int px_hi = scl(ap.patch_x_high);
		const int yc    = scl(ap.y);
		const int xc    = scl(ap.x);

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
		s.stacking_buffers.push_back(cv::Mat(wyx.rows, wyx.cols, buf_type,
		                                     cv::Scalar(0, 0, 0)));

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

	/* Background buffer (pre-drizzle). */
	s.averaged_background = cv::Mat::zeros(s.dim_y, s.dim_x, buf_type);

	/* Count pixels where sum_weights < background_blend_threshold *
	 * stack_size; if a small fraction, subdivide the intersection into
	 * background_patch_size patches and keep only those patches that
	 * contain at least one such pixel. */
	const float blend_t = (float) cfg.stack_frames_background_blend_threshold;
	cv::Mat needed_mask;
	cv::compare(s.sum_single_frame_weights, blend_t * stack_size_f, needed_mask, cv::CMP_LT);
	const int points_where_used = cv::countNonZero(needed_mask);
	const double total_drizzled = (double) s.dim_y_drizzled * (double) s.dim_x_drizzled;
	if ((double) points_where_used / total_drizzled
	    < cfg.stack_frames_background_fraction) {
		const int patch_size = cfg.stack_frames_background_patch_size;
		for (int py = 0; py < s.dim_y; py += patch_size) {
			/* Use `dim_y - 1` as the cap to match the reference algorithm. */
			const int py_hi = std::min(py + patch_size, s.dim_y - 1);
			if (py == py_hi) continue;
			for (int px = 0; px < s.dim_x; px += patch_size) {
				const int px_hi = std::min(px + patch_size, s.dim_x - 1);
				if (px == px_hi) continue;
				const cv::Mat slice = needed_mask(
				    cv::Range(scl(py), scl(py_hi)),
				    cv::Range(scl(px), scl(px_hi)));
				if (cv::countNonZero(slice) > 0)
					s.background_patches.push_back({py, py_hi, px, px_hi});
			}
		}
	}
	return s;
}

mpp_shifts_t *stack_compute_shifts_streamed(const FrameProvider &provider,
                                            int N,
                                            const cv::Mat &mean_frame_raw,
                                            const mpp_aps_t &aps,
                                            const APQualities &apq,
                                            const std::vector<FrameOffset> &offsets,
                                            const mpp_config_t &cfg,
                                            const int *included) {
	const int M = aps.count;
	/* Sub-pixel correlation is worthwhile when any drizzle upscale is
	 * being asked for — the dobox path uses the fractional shifts and
	 * the bicubic no-op path does not, but turning it on for the latter
	 * costs nothing relative to the integer fit. */
	const bool use_subpixel = (cfg.drizzle_scale > 1.001);

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

	int n_included = 0;
	if (included) {
		for (int f = 0; f < N; ++f) if (included[f]) ++n_included;
	} else {
		n_included = N;
	}
	int cur_nb = 0;
	for (int f = 0; f < N; ++f) {
		if (included && !included[f]) continue;
		if (!processing_should_continue()) {
			out->failure_counter = -1;   /* cancellation sentinel */
			return out;
		}
		/* Stage B correlation works at integer pixel grids — search-box
		 * cropping needs integer indices. Round the sub-pixel global
		 * offset; the rounded-away residual is the sub-pixel position
		 * where the AP's content actually sits relative to the integer
		 * search box. Sign trace for global_shifts = -0.3 (frame
		 * content drifted DOWN by 0.3):
		 *   offsets[f].dy = intersection[0] - (-0.3) = +0.3.
		 *   dy_rounded     = 0.
		 *   sub_y          = +0.3.
		 *   Search box at frame_y = ap.box_y_low + 0.
		 *   AP content sits at ap.box_y_low + 0.3 → correlation peak
		 *     at +0.3 from search centre → r.dy = -0.3 (matches the
		 *     algorithm's sign convention).
		 *   We want stored per-AP shift = 0 so that drizzle's
		 *     gdy + weighted_ap = -0.3 + 0 = -0.3 (matches the true
		 *     alignment correction). That means r.dy + sub_y = 0. ✓
		 * If we used r.dy − sub_y the sub-pixel global residual would
		 * be DOUBLED in the drizzle path and you get the high-frequency
		 * CFA-phase striping the user sees on real data. The classical
		 * stack path rounds shifts to whole drizzled pixels and is
		 * unaffected by this ≤0.5 px correction. */
		const double dy_full = offsets[f].dy;
		const double dx_full = offsets[f].dx;
		const int dy = (int) std::lround(dy_full);
		const int dx = (int) std::lround(dx_full);
		const double sub_y = dy_full - (double) dy;
		const double sub_x = dx_full - (double) dx;
		const cv::Mat frame = provider(f);   /* fresh blur in streaming mode;
		                                      * cached vector view in cached mode */
		for (int a : apq.used_alignment_points[f]) {
			const auto &ap = aps.records[a];
			const MultilevelShiftResult r = multilevel_correlation(
			    ref_boxes[a].second_phase, ref_boxes[a].first_phase,
			    frame,
			    ap.box_y_low + dy, ap.box_y_high + dy,
			    ap.box_x_low + dx, ap.box_x_high + dx,
			    cfg.frames_gauss_width, cfg.alignment_points_search_width,
			    use_subpixel, weight_matrix);
			const size_t off = (size_t) (f * M + a) * 2;
			out->shifts[off + 0] = r.dy + sub_y;
			out->shifts[off + 1] = r.dx + sub_x;
			out->success[f * M + a] = r.success ? 1 : 0;
			if (!r.success) ++failures;
		}
		g_atomic_int_inc(&cur_nb);
		/* Streaming mode pays full read+blur cost in this loop (no
		 * separate read pass), so map progress across the full bar.
		 * Cached mode's caller already moved through 0..0.5 during its
		 * own read pass and the 0.5 base here means we restart the
		 * second half — preserves the old visible behaviour. */
		gui_iface.set_progress(0.5 + 0.5 * (double) cur_nb / (double) n_included, NULL);
	}
	out->failure_counter = failures;
	return out;
}

mpp_shifts_t *stack_compute_shifts(const std::vector<cv::Mat> &frames_mono_blurred,
                                   const cv::Mat &mean_frame_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const std::vector<FrameOffset> &offsets,
                                   const mpp_config_t &cfg,
                                   const int *included) {
	return stack_compute_shifts_streamed(
	    [&frames_mono_blurred](int i) { return frames_mono_blurred[i]; },
	    (int) frames_mono_blurred.size(),
	    mean_frame_raw, aps, apq, offsets, cfg, included);
}

StackLoopOutput stack_apply_shifts_streamed(const FrameProvider &provider,
                                            int N,
                                            int num_layers,
                                            const mpp_aps_t &aps,
                                            const APQualities &apq,
                                            const mpp_shifts_t *shifts,
                                            const std::vector<FrameOffset> &offsets,
                                            const std::vector<double> &frame_brightness,
                                            const std::vector<int> &quality_sorted_idx,
                                            const cv::Vec4i &intersection,
                                            const mpp_config_t &cfg,
                                            const int *included,
                                            int max_threads,
                                            bool provider_thread_safe,
                                            size_t mem_budget_bytes) {
	StackLoopOutput out;
	/* Two parallelism axes, both bit-identical to the serial result:
	 *   - Decode (the dominant per-frame cost: read + debayer + float
	 *     convert + optional resize) is done a batch at a time across
	 *     threads, when the provider is reentrant.
	 *   - The per-frame AP accumulation parallelises over the disjoint
	 *     per-AP stacking buffers (each AP appears at most once in
	 *     used_alignment_points[f] and writes only its own buffer).
	 * Frames are accumulated serially in index order, so every buffer sees
	 * the same summation order ⇒ deterministic regardless of thread count. */
	const int nt = max_threads > 0 ? max_threads : 1;
	/* Output scale. Fractional values (e.g. 1.5) are supported: every
	 * coordinate is original × S rounded to the nearest pixel, and each
	 * frame is cv::resize'd to S× before its patches are accumulated. */
	const double S = std::max(1.0, cfg.drizzle_scale);
	out.state = stack_prepare_for_blending(aps, intersection, apq.stack_size,
	                                       S, num_layers, cfg);
	out.shift_failure_counter = shifts ? shifts->failure_counter : 0;

	const int M = aps.count;

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

	/* Frames that actually participate, in index order. */
	std::vector<int> work;
	work.reserve(N);
	for (int f = 0; f < N; ++f)
		if (!included || included[f]) work.push_back(f);
	const int W = (int) work.size();
	const int n_included = W > 0 ? W : 1;
	const bool holes = out.state.number_stacking_holes > 0;

	/* Frame-parallel accumulation. Each thread runs the full per-frame
	 * pipeline (read + brightness + resize + per-AP remap) for its share of
	 * frames into PRIVATE per-AP buffers, then the buffers are summed. This
	 * keeps every core busy on both the (cheap, I/O-bound for mono) decode
	 * and the heavier AP accumulation without a phase barrier or a fixed
	 * I/O-vs-AP pool split — a thread blocked on the SER read lock simply
	 * lets others get on with AP work. The float sum is reordered by the
	 * reduction, so the output is not bit-identical across thread counts —
	 * measured ≤1 16-bit LSB (~1 % of pixels off by one) on an 8000-frame
	 * stack, negligible for a statistical mean. Reads must be reentrant to
	 * parallelise; otherwise run single-threaded. */
	const bool par = provider_thread_safe && nt > 1;
	int n_threads = par ? nt : 1;
	if (par && mem_budget_bytes > 0) {
		/* Per thread: one private copy of all per-AP buffers (+ background)
		 * plus one in-flight frame (float + resized). */
		size_t set_bytes = 0;
		for (const auto &b : out.state.stacking_buffers)
			set_bytes += b.total() * b.elemSize();
		if (holes && !out.state.averaged_background.empty())
			set_bytes += out.state.averaged_background.total()
			           * out.state.averaged_background.elemSize();
		const size_t frame_bytes = (size_t) out.state.dim_x_drizzled
		                         * (size_t) out.state.dim_y_drizzled
		                         * (size_t) num_layers * sizeof(float) * 2;
		const size_t per_thread = set_bytes + frame_bytes;
		int fit = (int) (mem_budget_bytes / std::max<size_t>(1, per_thread));
		if (fit < 1) fit = 1;
		if (fit < n_threads) n_threads = fit;
	}

	auto clamped_blit = [](const cv::Mat &src, cv::Mat &dst,
	                       int src_y, int src_x, int len_y, int len_x) {
		const int src_y_lo = std::max(0, src_y);
		const int src_y_hi = std::min(src.rows, src_y + len_y);
		const int src_x_lo = std::max(0, src_x);
		const int src_x_hi = std::min(src.cols, src_x + len_x);
		if (src_y_hi <= src_y_lo || src_x_hi <= src_x_lo) return;
		const int dst_y_lo = src_y_lo - src_y;
		const int dst_x_lo = src_x_lo - src_x;
		cv::Mat dst_view = dst(
		    cv::Range(dst_y_lo, dst_y_lo + (src_y_hi - src_y_lo)),
		    cv::Range(dst_x_lo, dst_x_lo + (src_x_hi - src_x_lo)));
		dst_view += src(cv::Range(src_y_lo, src_y_hi),
		                cv::Range(src_x_lo, src_x_hi));
	};

	/* Per-frame progress maps into the second half of Stage C's bar — the
	 * outer mpp_stack_apply used 0..0.5 for the frame read pass. */
	gint cur_nb = 0;
	gint cancelled = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
#endif
	{
		/* Per-thread private accumulators (zero-init). */
		std::vector<cv::Mat> lbuf(M);
		for (int a = 0; a < M; ++a)
			lbuf[a] = cv::Mat::zeros(out.state.stacking_buffers[a].size(),
			                         out.state.stacking_buffers[a].type());
		cv::Mat lbg;
		if (holes)
			lbg = cv::Mat::zeros(out.state.averaged_background.size(),
			                     out.state.averaged_background.type());
		RemapBorder lb;

#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
		for (int wi = 0; wi < W; ++wi) {
			if (g_atomic_int_get(&cancelled)) continue;
			if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
			const int f = work[wi];

			const cv::Mat frame_raw = provider(f);
			if (frame_raw.empty()) {
				gui_iface.set_progress(0.5 + 0.5 * (double) g_atomic_int_add(&cur_nb, 1)
				                       / (double) n_included, NULL);
				continue;
			}
			cv::Mat frame_f32;
			if (cfg.frames_normalization) {
				const double sc = median_brightness / (frame_brightness[f] + 1e-7);
				frame_raw.convertTo(frame_f32, CV_32F, sc);
			} else {
				frame_raw.convertTo(frame_f32, CV_32F);
			}
			cv::Mat frame_drizzled;
			if (S != 1.0)
				cv::resize(frame_f32, frame_drizzled,
				           cv::Size((int) std::lround(frame_f32.cols * S),
				                    (int) std::lround(frame_f32.rows * S)),
				           0, 0, cv::INTER_LINEAR);
			else
				frame_drizzled = frame_f32;
			/* Integer rigid placement in drizzled coords; round the global
			 * offset (sub-pixel residual isn't carried by this path). */
			const int dy = (int) std::lround(offsets[f].dy);
			const int dx = (int) std::lround(offsets[f].dx);
			const int dy_K = (int) std::lround(dy * S);
			const int dx_K = (int) std::lround(dx * S);

			for (int a : apq.used_alignment_points[f]) {
				const size_t soff = (size_t) (f * M + a) * 2;
				const double shift_y = shifts ? shifts->shifts[soff + 0] : 0.0;
				const double shift_x = shifts ? shifts->shifts[soff + 1] : 0.0;
				const int shift_y_K = (int) std::lround(shift_y * S);
				const int shift_x_K = (int) std::lround(shift_x * S);
				const auto &p = out.state.patch_drizzled[a];
				stack_remap_rigid(frame_drizzled, lbuf[a],
				                  dy_K - shift_y_K, dx_K - shift_x_K,
				                  p[0], p[1], p[2], p[3], lb);
			}

			if (holes && top_for_bg.find(f) != top_for_bg.end()) {
				if (!out.state.background_patches.empty()) {
					for (const auto &bp : out.state.background_patches) {
						cv::Mat dst = lbg(cv::Range(bp.y_low, bp.y_high),
						                  cv::Range(bp.x_low, bp.x_high));
						clamped_blit(frame_f32, dst, bp.y_low + dy, bp.x_low + dx,
						             bp.y_high - bp.y_low, bp.x_high - bp.x_low);
					}
				} else {
					clamped_blit(frame_f32, lbg, dy, dx,
					             out.state.dim_y, out.state.dim_x);
				}
			}
			gui_iface.set_progress(0.5 + 0.5 * (double) g_atomic_int_add(&cur_nb, 1)
			                       / (double) n_included, NULL);
		}

		/* Sum the private buffers into the shared accumulators (border is
		 * order-independent max). */
#ifdef _OPENMP
#pragma omp critical(mpp_stack_reduce)
#endif
		{
			for (int a = 0; a < M; ++a)
				out.state.stacking_buffers[a] += lbuf[a];
			if (holes)
				out.state.averaged_background += lbg;
			out.border.y_low  = std::max(out.border.y_low,  lb.y_low);
			out.border.y_high = std::max(out.border.y_high, lb.y_high);
			out.border.x_low  = std::max(out.border.x_low,  lb.x_low);
			out.border.x_high = std::max(out.border.x_high, lb.x_high);
		}
	}

	if (cancelled) {
		out.state.dim_y = 0;   /* cancellation sentinel — outer caller checks */
		return out;
	}

	if (out.state.number_stacking_holes > 0) {
		if (S != 1.0) {
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

StackLoopOutput stack_apply_shifts(const std::vector<cv::Mat> &frames_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const mpp_shifts_t *shifts,
                                   const std::vector<FrameOffset> &offsets,
                                   const std::vector<double> &frame_brightness,
                                   const std::vector<int> &quality_sorted_idx,
                                   const cv::Vec4i &intersection,
                                   const mpp_config_t &cfg,
                                   const int *included,
                                   int max_threads,
                                   size_t mem_budget_bytes) {
	const int N = (int) frames_raw.size();
	const int num_layers = N == 0 ? 1 : frames_raw[0].channels();
	/* In-memory frames are trivially safe to read concurrently. */
	return stack_apply_shifts_streamed(
	    [&frames_raw](int i) -> cv::Mat { return frames_raw[i]; },
	    N, num_layers, aps, apq, shifts, offsets, frame_brightness,
	    quality_sorted_idx, intersection, cfg, included, max_threads,
	    /*provider_thread_safe=*/true, mem_budget_bytes);
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

namespace {
/* Replicate a single-channel matrix into a multi-channel one of the same
 * spatial size, so cv::multiply / cv::divide can be used against a
 * multi-channel partner. No-op for n=1. */
cv::Mat broadcast_channels(const cv::Mat &single, int n) {
	if (n <= 1) return single;
	std::vector<cv::Mat> chans(n, single);
	cv::Mat out;
	cv::merge(chans, out);
	return out;
}
}  // namespace

cv::Mat stack_merge_alignment_point_buffers(const StackState &state,
                                            const RemapBorder &border,
                                            const mpp_aps_t &aps,
                                            const mpp_config_t &cfg) {
	const int C = state.num_layers;
	const int buf_type = (C == 3) ? CV_32FC3 : CV_32F;

	/* Global drizzled image buffer. */
	cv::Mat buf = cv::Mat::zeros(state.dim_y_drizzled, state.dim_x_drizzled, buf_type);

	/* For each AP: buf[patch] += stacking_buffer * weights_yx. */
	for (int a = 0; a < aps.count; ++a) {
		const auto &p = state.patch_drizzled[a];
		cv::Mat dst = buf(cv::Range(p[0], p[1]), cv::Range(p[2], p[3]));
		cv::Mat weighted;
		const cv::Mat w = broadcast_channels(state.weights_yx[a], C);
		cv::multiply(state.stacking_buffers[a], w, weighted);
		dst += weighted;
	}

	/* buf /= sum_single_frame_weights. sum_weights was initialised to
	 * 1e-30, so we never divide by zero. */
	const cv::Mat sw = broadcast_channels(state.sum_single_frame_weights, C);
	cv::divide(buf, sw, buf);

	/* Background blend: compute a foreground_weight ramp and lerp. */
	if (state.number_stacking_holes > 0) {
		const float denom = (float) (cfg.stack_frames_background_blend_threshold
		                            * (double) state.stack_size);
		cv::Mat fg_weight;
		state.sum_single_frame_weights.convertTo(fg_weight, CV_32F, 1.0 / denom);
		cv::min(fg_weight, 1.0f, fg_weight);
		cv::max(fg_weight, 0.0f, fg_weight);
		const cv::Mat fgw = broadcast_channels(fg_weight, C);
		cv::Mat diff;
		cv::subtract(buf, state.averaged_background, diff);
		cv::multiply(diff, fgw, diff);
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

	const double in_max = (cfg.bitdepth == 8) ? 255.0 : 65535.0;
	cv::Mat normalised;
	buf.convertTo(normalised, buf_type, 1.0 / in_max);
	cv::min(normalised, 1.0f, normalised);
	cv::max(normalised, 0.0f, normalised);

	const int out_type = (C == 3) ? CV_16UC3 : CV_16U;
	cv::Mat out(normalised.rows, normalised.cols, out_type);
	/* skimage img_as_uint does floor(x * 65535 + 0.5). Multi-channel-aware
	 * flat iteration over `cols × channels()`. */
	for (int y = 0; y < normalised.rows; ++y) {
		const float *src = normalised.ptr<float>(y);
		uint16_t *dst = out.ptr<uint16_t>(y);
		const int n = normalised.cols * C;
		for (int x = 0; x < n; ++x) {
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
