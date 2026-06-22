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
#include "core/siril_log.h"

#include <set>

namespace mpp {

float stack_selection_weight(int rank, int stack_size, int taper) {
	if (rank < 0) return 0.0f;
	if (taper <= 0)
		return rank < stack_size ? 1.0f : 0.0f;
	const int plateau = stack_size - taper;
	if (rank < plateau) return 1.0f;
	if (rank >= stack_size + taper) return 0.0f;
	/* Half-sample-centred raised cosine across the 2·taper ramp: the
	 * cosine terms over the symmetric grid sum to zero, so Σw over all
	 * ranks is exactly plateau + taper = stack_size. */
	const double t = ((double) rank + 0.5 - (double) plateau)
	               / (2.0 * (double) taper);
	return (float) (0.5 * (1.0 + std::cos(G_PI * t)));
}

namespace {

/* Soft-selection taper half-width for this run.
 *
 * The taper exists to make adjacent APs give a borderline frame nearly
 * the same weight, so its natural width is the amount a frame's quality
 * RANK jitters between adjacent APs around the cut at rank N. Where the
 * quality-vs-rank curve is flat (many near-interchangeable frames) the
 * jitter is large and a wide taper costs nothing in sharpness; where it
 * is steep the jitter is small and the taper stays narrow, preserving
 * selectivity. Measured as: for a bounded sample of adjacent AP pairs,
 * the median |rank_b(f) − rank_a(f)| over frames ranked within ±K of N
 * in AP a; T = 2 × median, clamped to [max(2, N/10), N/2] and to the
 * hard limits T ≤ N − 1 (the plateau must be non-empty) and
 * T ≤ num_frames − N (no taper without room above the cut — this also
 * keeps the default 100% selection bit-identical to the hard path).
 *
 * Explicit frame-number mode (alignment_points_frame_number > 0) always
 * returns 0: that option's contract is an exact N-frame stack. */
int stack_compute_selection_taper(const std::vector<std::vector<double>> &qualities,
                                  const mpp_aps_t &aps, int num_frames,
                                  int stack_size, const mpp_config_t &cfg) {
	if (cfg.alignment_points_frame_number > 0)
		return 0;
	const int N = stack_size;
	const int room = num_frames - N;
	const int hard_cap = std::min({N / 2, N - 1, room});
	if (hard_cap <= 0 || aps.count < 2)
		return 0;

	/* Adjacent AP pairs (centre distance within 1.5 × grid step),
	 * deterministically subsampled to a bounded count. */
	const int step = mpp_cfg_step_size(&cfg);
	const int thr = (step * 3) / 2;
	std::vector<std::pair<int, int>> pairs;
	for (int a = 0; a < aps.count; ++a) {
		const auto &ra = aps.records[a];
		for (int b = a + 1; b < aps.count; ++b) {
			const auto &rb = aps.records[b];
			if (std::abs(ra.y - rb.y) <= thr && std::abs(ra.x - rb.x) <= thr)
				pairs.emplace_back(a, b);
		}
	}
	const int MAX_PAIRS = 200;
	const int stride = (int) pairs.size() > MAX_PAIRS
	                 ? (int) (pairs.size() / MAX_PAIRS) + 1 : 1;

	/* Per-AP rank-of-frame arrays, computed lazily for sampled APs only
	 * (full argsort; the partial sort used for selection doesn't give
	 * ranks past the cut). Ties broken by frame index for determinism. */
	std::vector<std::vector<int>> rank_of(aps.count);
	auto ranks = [&](int a) -> const std::vector<int> & {
		if (rank_of[a].empty()) {
			const auto &q = qualities[a];
			std::vector<int> order(num_frames);
			std::iota(order.begin(), order.end(), 0);
			std::sort(order.begin(), order.end(),
			          [&q](int x, int y) {
				return q[x] > q[y] || (q[x] == q[y] && x < y);
			});
			rank_of[a].resize(num_frames);
			for (int r = 0; r < num_frames; ++r)
				rank_of[a][order[r]] = r;
		}
		return rank_of[a];
	};

	const int K = std::max(3, N / 10);
	std::vector<int> jitter;
	for (size_t pi = 0; pi < pairs.size(); pi += stride) {
		const int a = pairs[pi].first, b = pairs[pi].second;
		const auto &rka = ranks(a);
		const auto &rkb = ranks(b);
		for (int f = 0; f < num_frames; ++f) {
			const int r = rka[f];
			if (r < N - K || r >= N + K) continue;
			jitter.push_back(std::abs(rkb[f] - r));
		}
	}

	int T;
	if (jitter.empty()) {
		T = std::max(2, N / 10);   /* no measurable pairs — floor */
	} else {
		std::nth_element(jitter.begin(), jitter.begin() + jitter.size() / 2,
		                 jitter.end());
		const int median = jitter[jitter.size() / 2];
		T = std::max(2 * median, std::max(2, N / 10));
	}
	return std::min(T, hard_cap);
}

}  // namespace

std::vector<float> stack_one_dim_weight(int patch_low, int patch_high, int box_center,
                                        bool extend_low, bool extend_high) {
	const int patch_size    = patch_high - patch_low;
	const int center_offset = box_center - patch_low;
	std::vector<float> w(patch_size, 0.0f);

	/* Raised-cosine (Hann) taper. PSS uses the linear ramp t directly;
	 * mapping it through sin²(π/2·t) keeps the endpoints (t=1 → 1, t→0
	 * → 0) but zeroes the derivative at both, so the 2D blend window has
	 * no gradient kink at AP centres or patch boundaries — small per-AP
	 * level differences no longer imprint the AP lattice as creases the
	 * way the C0 triangular window did. The per-pixel division by
	 * sum_single_frame_weights normalises whatever window is used, so no
	 * partition-of-unity property is required. */
	auto hann = [](double t) {
		const double s = std::sin(G_PI_2 * t);
		return (float) (s * s);
	};

	if (extend_low) {
		for (int i = 0; i < center_offset; ++i) w[i] = 1.0f;
	} else {
		/* ramp t = [1/(c+1), 2/(c+1), …, c/(c+1)]. */
		const double denom = (double) (center_offset + 1);
		for (int i = 0; i < center_offset; ++i)
			w[i] = hann((double) (i + 1) / denom);
	}

	const int high_len = patch_size - center_offset;  /* = patch_high - box_center */
	if (extend_high) {
		for (int i = 0; i < high_len; ++i) w[center_offset + i] = 1.0f;
	} else {
		/* ramp t = [N/N, (N-1)/N, …, 1/N]. */
		const double denom = (double) high_len;
		for (int i = 0; i < high_len; ++i)
			w[center_offset + i] = hann((double) (high_len - i) / denom);
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
                       RemapBorder &border, float weight) {
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
	if (weight == 1.0f)
		dst += src;
	else
		cv::scaleAdd(src, weight, dst, dst);
}

void stack_remap_subpixel(const cv::Mat &frame, cv::Mat &buffer,
                          double shift_y, double shift_x,
                          int y_low, int y_high, int x_low, int x_high,
                          RemapBorder &border, float weight) {
	/* Fractional-shift companion to stack_remap_rigid: adds
	 * frame[y_low+shift_y : y_high+shift_y, ...] into the buffer with the
	 * fractional part of the shift resolved by Lanczos interpolation
	 * instead of being rounded away. A shift within EPS of an integer
	 * takes the exact integer blit — no interpolation softening, and the
	 * historic path stays bit-reproducible for integer data. */
	const double EPS = 1e-3;
	const double fy = shift_y - std::floor(shift_y);
	const double fx = shift_x - std::floor(shift_x);
	if ((fy < EPS || fy > 1.0 - EPS) && (fx < EPS || fx > 1.0 - EPS)) {
		stack_remap_rigid(frame, buffer,
		                  (int) std::lround(shift_y), (int) std::lround(shift_x),
		                  y_low, y_high, x_low, x_high, border, weight);
		return;
	}

	/* Clip the destination so every sample position stays inside the
	 * frame; track the clipped extent for the final border trim, mirroring
	 * the integer path. */
	const int patch_h = y_high - y_low;
	const int patch_w = x_high - x_low;
	int y_low_dst = 0;
	const double y_first = (double) y_low + shift_y;
	if (y_first < 0.0) {
		y_low_dst = (int) std::ceil(-y_first);
		border.y_low = std::max(border.y_low, y_low_dst);
	}
	int y_high_dst = patch_h;
	const double y_last = y_first + (double) (patch_h - 1);
	if (y_last > (double) (frame.rows - 1)) {
		const int clip = (int) std::ceil(y_last - (double) (frame.rows - 1));
		y_high_dst = patch_h - clip;
		border.y_high = std::max(border.y_high, clip);
	}
	int x_low_dst = 0;
	const double x_first = (double) x_low + shift_x;
	if (x_first < 0.0) {
		x_low_dst = (int) std::ceil(-x_first);
		border.x_low = std::max(border.x_low, x_low_dst);
	}
	int x_high_dst = patch_w;
	const double x_last = x_first + (double) (patch_w - 1);
	if (x_last > (double) (frame.cols - 1)) {
		const int clip = (int) std::ceil(x_last - (double) (frame.cols - 1));
		x_high_dst = patch_w - clip;
		border.x_high = std::max(border.x_high, clip);
	}
	if (y_high_dst <= y_low_dst || x_high_dst <= x_low_dst)
		return;

	/* Inverse-map translation: dst(x', y') samples frame at
	 * (x' + x_low + x_low_dst + shift_x, y' + y_low + y_low_dst + shift_y).
	 * Lanczos support pixels that fall just outside the frame (≤3 px,
	 * only at frame edges) are replicated — the affected ring is inside
	 * the border trim. warpAffine's cost scales with the dst patch, not
	 * the frame, so passing the whole frame as src costs nothing. */
	const cv::Mat M = (cv::Mat_<double>(2, 3)
	    << 1.0, 0.0, (double) (x_low + x_low_dst) + shift_x,
	       0.0, 1.0, (double) (y_low + y_low_dst) + shift_y);
	cv::Mat patch;
	cv::warpAffine(frame, patch, M,
	               cv::Size(x_high_dst - x_low_dst, y_high_dst - y_low_dst),
	               cv::INTER_LANCZOS4 | cv::WARP_INVERSE_MAP,
	               cv::BORDER_REPLICATE);
	cv::Mat dst = buffer(cv::Range(y_low_dst, y_high_dst),
	                     cv::Range(x_low_dst, x_high_dst));
	if (weight == 1.0f)
		dst += patch;
	else
		cv::scaleAdd(patch, weight, dst, dst);
}

APQualities ap_compute_frame_qualities_streamed(const FrameProvider &provider,
                                                int N,
                                                const std::vector<double> &frame_brightness,
                                                const mpp_aps_t &aps,
                                                const std::vector<FrameOffset> &offsets,
                                                int frame_rows, int frame_cols,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress,
                                                void *progress_user,
                                                const std::vector<int> &included,
                                                int max_threads,
                                                bool provider_thread_safe) {
	APQualities out;
	const int M = aps.count;
	const int stride = cfg.align_frames_sampling_stride;
	const bool normalise = cfg.frames_normalization;
	if (N == 0 || M == 0 || (int) frame_brightness.size() != N
	 || (int) offsets.size() != N || stride <= 0)
		return out;

	/* Inclusion mask (empty = all). Excluded frames get a sentinel quality
	 * below the lowest real σ (≥ 0), so std::partial_sort never lifts them
	 * into any AP's best_frame_indices; the selection count is clamped to
	 * the included-frame count so a sentinel can't fill a slot either. */
	const bool have_mask = (int) included.size() == N;
	auto is_included = [&](int i) { return !have_mask || included[i]; };
	int n_included = 0;
	for (int i = 0; i < N; ++i) if (is_included(i)) ++n_included;
	if (n_included == 0) return out;

	/* stack_size: explicit override, else %% of N (rounded up). */
	if (cfg.alignment_points_frame_number > 0)
		out.stack_size = cfg.alignment_points_frame_number;
	else
		out.stack_size = std::max(
		    (int) std::ceil((double) N * cfg.alignment_points_frame_percent / 100.0),
		    1);
	out.stack_size = std::min(out.stack_size, n_included);

	out.qualities.assign(M, std::vector<double>(N, 0.0));

	/* Frame-parallel: each frame writes its own column out.qualities[*][f], so
	 * the result is independent of thread count. The heavy per-frame work
	 * (read + optional derotation remap + Laplacian + per-AP meanStdDev)
	 * dominates once derotation is active, so spreading frames across cores is
	 * what keeps them busy. Only engaged when the provider is reentrant.
	 * Cancellation / OOM leave the region via flags, never a return or an
	 * escaping exception. */
	const int nt = (provider_thread_safe && max_threads > 1) ? max_threads : 1;
	/* The frame loop is the parallel level — stop OpenCV from *also* threading
	 * each per-frame op (the derot remap, the Laplacian, every meanStdDev).
	 * With OpenCV's TBB backend, nt OpenMP workers each submitting to TBB's one
	 * global arena oversubscribe and stall on the arena, so utilisation
	 * collapses to a fraction of the cores. Pinning OpenCV to one thread makes
	 * each frame's ops run inline in their OpenMP worker; the cores stay busy
	 * because many frames are in flight. Restored right after the loop. When
	 * nt==1 we leave OpenCV alone so the lone frame still uses every core. */
	const int cv_saved_threads = cv::getNumThreads();
	if (nt > 1) cv::setNumThreads(1);
	gint cancelled = 0, cur_nb = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic) if (nt > 1)
#endif
	for (int f = 0; f < N; ++f) {
		if (g_atomic_int_get(&cancelled)) continue;
		if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
		if (!is_included(f)) {
			/* Sentinel below the lowest real σ — never selected. No read. */
			for (int a = 0; a < M; ++a) out.qualities[a][f] = -1.0;
			if (progress) progress((double) (g_atomic_int_add(&cur_nb, 1) + 1)
			                       / (double) N, progress_user);
			continue;
		}
		try {
			const cv::Mat frame = provider(f);
			if (!frame.empty()) {
				const cv::Mat lap = rank_blurred_laplacian_u8(frame, cfg);
				/* Round sub-pixel offsets — Stage A's box-slicing math wants
				 * integer indices; the sub-pixel residual matters only at the
				 * drizzle pixmap stage. */
				const int dy = (int) std::lround(offsets[f].dy);
				const int dx = (int) std::lround(offsets[f].dx);
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
		} catch (const std::exception &) {
			/* OOM in a worker: abandon the pass cleanly (sentinel below). */
			g_atomic_int_set(&cancelled, 1);
		}
		if (progress) progress((double) (g_atomic_int_add(&cur_nb, 1) + 1)
		                       / (double) N, progress_user);
	}
	if (nt > 1) cv::setNumThreads(cv_saved_threads);
	if (g_atomic_int_get(&cancelled)) {
		out.stack_size = 0;   /* cancellation / OOM sentinel */
		return out;
	}

	/* Soft-selection taper (0 in explicit frame-number mode and whenever
	 * there is no room above the cut — e.g. the default 100%). */
	out.taper = stack_compute_selection_taper(out.qualities, aps, N,
	                                          out.stack_size, cfg);
	const int selected = std::min(n_included, out.stack_size + out.taper);

	/* Per-AP: pick the top `selected` frame indices by descending quality
	 * (index tie-break for determinism); entry k carries selection weight
	 * stack_selection_weight(k, stack_size, taper). */
	out.best_frame_indices.assign(M, {});
	out.used_alignment_points.assign(N, {});
	for (int a = 0; a < M; ++a) {
		std::vector<int> indices(N);
		std::iota(indices.begin(), indices.end(), 0);
		const auto &q = out.qualities[a];
		std::partial_sort(indices.begin(), indices.begin() + selected,
		                  indices.end(),
		                  [&q](int x, int y) {
			return q[x] > q[y] || (q[x] == q[y] && x < y);
		});
		out.best_frame_indices[a].assign(indices.begin(),
		                                 indices.begin() + selected);
		for (int i = 0; i < selected; ++i) {
			const float w = stack_selection_weight(i, out.stack_size,
			                                       out.taper);
			if (w > 0.0f)
				out.used_alignment_points[indices[i]].push_back({a, w});
		}
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
                                       void *progress_user,
                                       const std::vector<int> &included) {
	return ap_compute_frame_qualities_streamed(
	    [&frames](int i) { return frames[i]; },
	    (int) frames.size(), frame_brightness, aps, offsets,
	    frame_rows, frame_cols, cfg, progress, progress_user, included);
}

StackState stack_prepare_for_blending(const mpp_aps_t &aps,
                                      const cv::Vec4i &intersection,
                                      int stack_size,
                                      double drizzle_scale,
                                      int num_layers,
                                      const mpp_config_t &cfg,
                                      const std::vector<float> *ap_effective_counts) {
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
	s.ap_frame_counts.reserve(M);
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

		/* Separable product, not PSS's min(wy, wx): the product of two
		 * C1 Hann tapers is C1 everywhere, while min() recreates
		 * diagonal derivative creases along |wy| = |wx|. */
		cv::Mat wyx((int) wy.size(), (int) wx.size(), CV_32F);
		for (int y = 0; y < (int) wy.size(); ++y) {
			float *row = wyx.ptr<float>(y);
			for (int x = 0; x < (int) wx.size(); ++x)
				row[x] = wy[y] * wx[x];
		}
		s.weights_yx.push_back(wyx);
		s.stacking_buffers.push_back(cv::Mat(wyx.rows, wyx.cols, buf_type,
		                                     cv::Scalar(0, 0, 0)));

		/* Accumulate (per-AP weight total) × weights_yx into the global
		 * buffer.  The caller passes the per-AP accumulated-weight totals
		 * (soft selection weights over included frames, less any dropped
		 * failed pairs) so the normalisation matches exactly what gets
		 * accumulated; stack_size is the fallback for direct callers. */
		const float ap_count_f = ap_effective_counts
		    ? (*ap_effective_counts)[a] : stack_size_f;
		s.ap_frame_counts.push_back(ap_count_f);
		cv::Mat dst = s.sum_single_frame_weights(
		    cv::Range(py_lo, py_hi), cv::Range(px_lo, px_hi));
		dst += ap_count_f * wyx;
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
                                            const int *included,
                                            int max_threads,
                                            bool provider_thread_safe) {
	const int M = aps.count;
	/* Sub-pixel correlation is always on: the accumulation path places
	 * each patch at its fractional shift (stack_remap_subpixel), so an
	 * integer-only fit would quantise every AP to the pixel grid and
	 * re-introduce the block-lattice phase noise the sub-pixel paste
	 * exists to remove. The parabolic refinement is a 3×3 fit on the
	 * already-computed correlation surface — negligible cost. */
	const bool use_subpixel = true;

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

	int n_included = 0;
	if (included) {
		for (int f = 0; f < N; ++f) if (included[f]) ++n_included;
	} else {
		n_included = N;
	}
	if (n_included <= 0) n_included = 1;

	/* Frame-parallel: every (frame, AP) shift is written to a disjoint slot, so
	 * the result is independent of thread count — bit-identical to the serial
	 * loop. The heavy per-frame work (read + optional derotation remap + blur)
	 * dominates once derotation is active, so spreading frames across cores is
	 * what keeps them busy. Only engaged when the provider is reentrant
	 * (cached frames or thread-safe reads); otherwise nt collapses to 1.
	 * Cancellation and OOM leave the OpenMP region via flags — never a bare
	 * return or an escaping exception, both of which are UB inside a parallel
	 * for. */
	const int nt = (provider_thread_safe && max_threads > 1) ? max_threads : 1;
	/* See ap_compute_frame_qualities_streamed: the frame loop is the parallel
	 * level, so pin OpenCV to one thread for the per-frame ops (derot remap +
	 * blur + the multilevel correlation's cv calls). Otherwise nt OpenMP
	 * workers contend on OpenCV's single TBB arena and utilisation collapses.
	 * Restored right after the loop; left untouched when nt==1. */
	const int cv_saved_threads = cv::getNumThreads();
	if (nt > 1) cv::setNumThreads(1);
	gint cancelled = 0, failures = 0, cur_nb = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic) if (nt > 1)
#endif
	for (int f = 0; f < N; ++f) {
		if (included && !included[f]) continue;
		if (g_atomic_int_get(&cancelled)) continue;
		if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
		try {
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
			 * be DOUBLED in the paste and you get high-frequency CFA-phase
			 * striping on real data. The stack path places each patch at
			 * offsets.dy − shift (sub-pixel, see stack_apply_shifts), so
			 * this ≤0.5 px correction lands directly in the output. */
			const double dy_full = offsets[f].dy;
			const double dx_full = offsets[f].dx;
			const int dy = (int) std::lround(dy_full);
			const int dx = (int) std::lround(dx_full);
			const double sub_y = dy_full - (double) dy;
			const double sub_x = dx_full - (double) dx;
			const cv::Mat frame = provider(f);   /* fresh blur in streaming mode;
			                                      * cached vector view in cached mode */
			if (!frame.empty()) {
				for (const auto &u : apq.used_alignment_points[f]) {
					const int a = u.ap;
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
					if (!r.success) g_atomic_int_inc(&failures);
				}
			}
		} catch (const std::exception &) {
			/* OOM in a worker: abandon the pass cleanly (sentinel below). */
			g_atomic_int_set(&cancelled, 1);
		}
		/* Single streaming pass (read + blur + correlate happen together),
		 * so map progress across Stage B's whole 0..1 bar. */
		gui_iface.set_progress((double) (g_atomic_int_add(&cur_nb, 1) + 1)
		                       / (double) n_included, NULL);
	}
	if (nt > 1) cv::setNumThreads(cv_saved_threads);
	if (g_atomic_int_get(&cancelled)) {
		out->failure_counter = -1;   /* cancellation / OOM sentinel */
		return out;
	}
	out->failure_counter = (int) g_atomic_int_get(&failures);
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
	const int M = aps.count;

	/* Per-AP accumulated-weight totals (selection weights over included
	 * frames). These feed stack_prepare_for_blending so the normalisation
	 * matches exactly what the loop below accumulates — both for soft
	 * selection (fractional weights) and for frames excluded after
	 * Analyze (whose missing contributions would otherwise render the
	 * patch dark).
	 *
	 * Optional drop of failed (frame, AP) contributions.  Decided up
	 * front from shifts->success so (a) the weight totals can be reduced
	 * to match, and (b) the frame-parallel accumulation below needs no
	 * extra synchronisation.  An AP whose used pairs ALL failed keeps
	 * stacking everything: a coarsely aligned patch beats a hole in the
	 * output. */
	const bool skip_failed = cfg.stack_skip_failed_aps && shifts != nullptr;
	std::vector<uint8_t> ap_skip_enabled;
	std::vector<float> ap_eff_counts(M, 0.0f);
	for (int f = 0; f < N; ++f) {
		if (included && !included[f]) continue;
		for (const auto &u : apq.used_alignment_points[f])
			ap_eff_counts[u.ap] += u.weight;
	}
	if (skip_failed) {
		std::vector<float> failed_w(M, 0.0f);
		std::vector<int> failed_n(M, 0);
		for (int f = 0; f < N; ++f) {
			if (included && !included[f]) continue;
			for (const auto &u : apq.used_alignment_points[f]) {
				if (!shifts->success[(size_t) f * M + u.ap]) {
					failed_w[u.ap] += u.weight;
					++failed_n[u.ap];
				}
			}
		}
		ap_skip_enabled.assign(M, 0);
		int dropped = 0, saturated_aps = 0, skipped_aps = 0;
		for (int a = 0; a < M; ++a) {
			if (failed_n[a] == 0) continue;
			if (failed_w[a] >= ap_eff_counts[a] - 1e-6f) {
				++saturated_aps;   /* nothing measurable at this AP */
				continue;
			}
			ap_skip_enabled[a] = 1;
			++skipped_aps;
			ap_eff_counts[a] -= failed_w[a];
			dropped += failed_n[a];
		}
		if (dropped > 0)
			siril_log_message(_("Stack: skipped failed shift measurements at "
			                    "%d of %d alignment points (%d frame/alignment-"
			                    "point pairs dropped)\n"),
			                  skipped_aps, M, dropped);
		if (saturated_aps > 0)
			siril_log_message(_("Stack: %d alignment point(s) had no "
			                    "successful shifts — stacking them "
			                    "unfiltered\n"), saturated_aps);
	}

	out.state = stack_prepare_for_blending(aps, intersection, apq.stack_size,
	                                       S, num_layers, cfg,
	                                       &ap_eff_counts);
	out.shift_failure_counter = shifts ? shifts->failure_counter : 0;

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

	/* Single streaming pass (read + drizzle + remap happen together), so
	 * per-frame progress maps across Stage C's whole 0..1 bar. */
	gint cur_nb = 0;
	gint cancelled = 0;
	gint oom = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
#endif
	{
		/* Per-thread private accumulators (zero-init). Allocation and the
		 * per-frame work below are try/catch-contained: an exception
		 * escaping an OpenMP region is std::terminate(), and on Windows
		 * (no overcommit) cv allocations genuinely fail under memory
		 * pressure. A thread that trips the flag keeps encountering the
		 * worksharing constructs (required by OpenMP) but skips the work;
		 * the caller discards the partial result. */
		bool local_ok = true;
		std::vector<cv::Mat> lbuf(M);
		cv::Mat lbg;
		try {
			for (int a = 0; a < M; ++a)
				lbuf[a] = cv::Mat::zeros(out.state.stacking_buffers[a].size(),
				                         out.state.stacking_buffers[a].type());
			if (holes)
				lbg = cv::Mat::zeros(out.state.averaged_background.size(),
				                     out.state.averaged_background.type());
		} catch (const std::exception &) {
			local_ok = false;
			g_atomic_int_set(&oom, 1);
		}
		RemapBorder lb;

#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
		for (int wi = 0; wi < W; ++wi) {
			if (g_atomic_int_get(&cancelled) || g_atomic_int_get(&oom)) continue;
			if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
			const int f = work[wi];
			try {

			const cv::Mat frame_raw = provider(f);
			if (frame_raw.empty()) {
				gui_iface.set_progress((double) g_atomic_int_add(&cur_nb, 1)
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
			/* Sub-pixel rigid placement in drizzled coords. The full
			 * source position for a patch is (offsets − per-AP shift) ×
			 * scale; stack_remap_subpixel resolves the fractional part by
			 * interpolation, so neighbouring APs no longer disagree by the
			 * up-to-±0.5 px rounding that imprints the AP lattice on the
			 * stack. The background blit below stays integer — it is a
			 * heavily blended fallback where sub-pixel accuracy is moot. */
			const int dy = (int) std::lround(offsets[f].dy);
			const int dx = (int) std::lround(offsets[f].dx);

			for (const auto &u : apq.used_alignment_points[f]) {
				const int a = u.ap;
				if (skip_failed && ap_skip_enabled[a]
				    && !shifts->success[(size_t) f * M + a])
					continue;
				const size_t soff = (size_t) (f * M + a) * 2;
				const double shift_y = shifts ? shifts->shifts[soff + 0] : 0.0;
				const double shift_x = shifts ? shifts->shifts[soff + 1] : 0.0;
				const auto &p = out.state.patch_drizzled[a];
				stack_remap_subpixel(frame_drizzled, lbuf[a],
				                     (offsets[f].dy - shift_y) * S,
				                     (offsets[f].dx - shift_x) * S,
				                     p[0], p[1], p[2], p[3], lb, u.weight);
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
			gui_iface.set_progress((double) g_atomic_int_add(&cur_nb, 1)
			                       / (double) n_included, NULL);

			} catch (const std::exception &) {
				g_atomic_int_set(&oom, 1);
				continue;
			}
		}

		/* Sum the private buffers into the shared accumulators (border is
		 * order-independent max). Skipped once oom trips — the result is
		 * discarded anyway, and a thread whose lbuf allocation failed
		 * holds empty Mats that cv refuses to add. */
#ifdef _OPENMP
#pragma omp critical(mpp_stack_reduce)
#endif
		{
			if (local_ok && !g_atomic_int_get(&oom)) {
				try {
					for (int a = 0; a < M; ++a)
						out.state.stacking_buffers[a] += lbuf[a];
					if (holes)
						out.state.averaged_background += lbg;
					out.border.y_low  = std::max(out.border.y_low,  lb.y_low);
					out.border.y_high = std::max(out.border.y_high, lb.y_high);
					out.border.x_low  = std::max(out.border.x_low,  lb.x_low);
					out.border.x_high = std::max(out.border.x_high, lb.x_high);
				} catch (const std::exception &) {
					g_atomic_int_set(&oom, 1);
				}
			}
		}
	}

	if (g_atomic_int_get(&oom)) {
		out.oom = true;
		return out;
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

/* Per-AP DC equalisation offsets (Siril enhancement; PSS has no
 * equivalent). Each AP buffer is a mean over its own top-N frame subset,
 * so neighbouring APs settle at slightly different local levels
 * (disjoint frame picks, transparency drift); the blend window then
 * turns each level step into block-pitch shading. Over the overlap of
 * two patches the scene content is identical, so the mean difference of
 * the two normalised buffers there is pure stacking inconsistency:
 *   d_ab = mean_overlap(B_a / c_a − B_b / c_b)   (per channel).
 * Solve the least-squares offsets δ minimising
 *   Σ_pairs n_ab (δ_a − δ_b + d_ab)²
 * (n_ab = overlap area) by Gauss-Seidel — the normal equations are
 *   δ_a = Σ_b n_ab (δ_b − d_ab) / Σ_b n_ab,
 * a strictly diagonally dominant relaxation that converges in a few
 * dozen sweeps for any AP graph. The system is invariant under a global
 * constant, so pin the mean δ of connected APs to zero (preserves the
 * overall brightness). APs with no overlapping neighbour keep δ = 0. */
std::vector<cv::Scalar> stack_solve_ap_dc_offsets(const mpp::StackState &state) {
	const int M = (int) state.stacking_buffers.size();
	std::vector<cv::Scalar> delta(M, cv::Scalar::all(0.0));
	if (M < 2 || (int) state.ap_frame_counts.size() != M)
		return delta;

	struct Edge { int b; double n; cv::Scalar d; };
	std::vector<std::vector<Edge>> nbrs(M);
	for (int a = 0; a < M; ++a) {
		const auto &pa = state.patch_drizzled[a];
		const double ca = state.ap_frame_counts[a];
		if (ca <= 0.0) continue;
		for (int b = a + 1; b < M; ++b) {
			const double cb = state.ap_frame_counts[b];
			if (cb <= 0.0) continue;
			const auto &pb = state.patch_drizzled[b];
			const int y0 = std::max(pa[0], pb[0]), y1 = std::min(pa[1], pb[1]);
			const int x0 = std::max(pa[2], pb[2]), x1 = std::min(pa[3], pb[3]);
			if (y1 <= y0 || x1 <= x0) continue;
			const cv::Mat sa = state.stacking_buffers[a](
			    cv::Range(y0 - pa[0], y1 - pa[0]),
			    cv::Range(x0 - pa[2], x1 - pa[2]));
			const cv::Mat sb = state.stacking_buffers[b](
			    cv::Range(y0 - pb[0], y1 - pb[0]),
			    cv::Range(x0 - pb[2], x1 - pb[2]));
			const double n = (double) (y1 - y0) * (double) (x1 - x0);
			const cv::Scalar d = cv::sum(sa) * (1.0 / (ca * n))
			                   - cv::sum(sb) * (1.0 / (cb * n));
			nbrs[a].push_back({b, n, d});
			nbrs[b].push_back({a, n, -d});
		}
	}

	for (int it = 0; it < 32; ++it) {
		for (int a = 0; a < M; ++a) {
			if (nbrs[a].empty()) continue;
			cv::Scalar acc = cv::Scalar::all(0.0);
			double nsum = 0.0;
			for (const auto &e : nbrs[a]) {
				acc += (delta[e.b] - e.d) * e.n;
				nsum += e.n;
			}
			delta[a] = acc * (1.0 / nsum);
		}
	}

	cv::Scalar mean = cv::Scalar::all(0.0);
	int connected = 0;
	for (int a = 0; a < M; ++a) {
		if (nbrs[a].empty()) continue;
		mean += delta[a];
		++connected;
	}
	if (connected > 0) {
		mean *= 1.0 / connected;
		for (int a = 0; a < M; ++a)
			if (!nbrs[a].empty()) delta[a] -= mean;
	}
	return delta;
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

	/* DC equalisation offsets — see stack_solve_ap_dc_offsets. */
	const std::vector<cv::Scalar> dc = stack_solve_ap_dc_offsets(state);
	{
		double max_dc = 0.0;
		for (const auto &d : dc)
			for (int c = 0; c < C; ++c)
				max_dc = std::max(max_dc, std::abs(d[c]));
		if (max_dc >= 0.005)
			siril_log_debug("Stack (mpp): AP DC equalisation, max "
			                "correction %.2f ADU (16-bit)\n", max_dc);
	}

	/* For each AP: buf[patch] += (stacking_buffer + c_a·δ_a) * weights_yx.
	 * The δ term rides under the same window, so after the division by
	 * sum_single_frame_weights it contributes exactly +δ_a wherever the
	 * AP carries weight. */
	for (int a = 0; a < aps.count; ++a) {
		const auto &p = state.patch_drizzled[a];
		cv::Mat dst = buf(cv::Range(p[0], p[1]), cv::Range(p[2], p[3]));
		cv::Mat weighted;
		const cv::Mat w = broadcast_channels(state.weights_yx[a], C);
		cv::multiply(state.stacking_buffers[a], w, weighted);
		dst += weighted;
		if (a < (int) state.ap_frame_counts.size()) {
			const cv::Scalar term = dc[a] * (double) state.ap_frame_counts[a];
			if (std::abs(term[0]) + std::abs(term[1]) + std::abs(term[2]) > 0.0) {
				cv::Mat corr(w.rows, w.cols, buf_type, term);
				cv::multiply(corr, w, corr);
				dst += corr;
			}
		}
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
