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
 * Multipoint stacking — the engine-independent pieces:
 *
 *   - per-AP per-frame quality ranking and soft selection (with spatial
 *     harmonisation of weak-structure APs),
 *   - Stage B per-(frame, AP) shift measurement (with escalation retry),
 *   - the blend-window / patch-bound geometry and level-ramp helpers the
 *     warp engine builds its dense fields from.
 *
 * Stage C accumulation itself lives in mpp_warp_stack.cpp: the per-AP
 * patch-blend engine (the PSS architecture this port started from) was
 * retired once the warp-field engine matched or beat it on solar, lunar
 * and planetary data. Each piece here sits behind a small C++ entry
 * point in mpp_stack_priv.hpp with its own unit test.
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

	/* Selection harmonisation for weak-structure APs. At APs whose box
	 * holds little rankable structure (solar limb, prominence fringes,
	 * near-sky patches) the per-frame Laplace σ is noise, so adjacent APs
	 * pick gratuitously different frame subsets; each subset carries its
	 * own seeing mix, the aureole each AP stacks differs in blur profile,
	 * and the blend imprints AP-pitch scallops along the limb — a shape
	 * difference no level matching can remove. Blend such APs' quality
	 * vectors toward a structure-weighted consensus of their neighbours:
	 * α = clamp(1 − 2·s/median(s), 0, 1) is zero for anything at or above
	 * half the grid's median structure (object selection untouched) and
	 * approaches full consensus for structure-free APs, whose ranking
	 * then follows the regional seeing measured by the nearest APs that
	 * can measure it — and converges across the ring. The published
	 * qualities stay raw (diagnostic contract); only the ranking below
	 * uses the harmonised copy. */
	const std::vector<std::vector<double>> *rank_q = &out.qualities;
	std::vector<std::vector<double>> harmonised;
	{
		std::vector<double> sv(M);
		for (int a = 0; a < M; ++a) sv[a] = aps.records[a].structure;
		std::vector<double> tmp(sv);
		std::nth_element(tmp.begin(), tmp.begin() + M / 2, tmp.end());
		const double med = tmp[M / 2];
		const int thr = (mpp_cfg_step_size(&cfg) * 5) / 2;
		bool any = false;
		if (med > 0.0) {
			harmonised = out.qualities;
			for (int a = 0; a < M; ++a) {
				const auto &ra = aps.records[a];
				/* Structure-weighted consensus of the neighbourhood
				 * (~2 grid rings). */
				std::vector<double> cons(N, 0.0);
				double wsum = 0.0;
				for (int b = 0; b < M; ++b) {
					if (b == a) continue;
					const auto &rb = aps.records[b];
					if (std::abs(ra.y - rb.y) > thr
					 || std::abs(ra.x - rb.x) > thr)
						continue;
					const double w = sv[b] + 1e-6;
					for (int f = 0; f < N; ++f)
						cons[f] += w * out.qualities[b][f];
					wsum += w;
				}
				if (wsum <= 0.0) continue;
				for (int f = 0; f < N; ++f) cons[f] /= wsum;
				/* Rank coherence: Pearson correlation of this AP's
				 * quality vector against the consensus over included
				 * frames. Seeing is regional at grid-step scale, so a
				 * measurable AP tracks its neighbours closely; a low
				 * correlation means the ranking is noise (limb-edge
				 * boxes score high structure yet rank incoherently). */
				double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
				int n = 0;
				for (int f = 0; f < N; ++f) {
					if (!is_included(f)) continue;
					const double x = out.qualities[a][f], y = cons[f];
					sx += x; sy += y; sxx += x * x; syy += y * y;
					sxy += x * y; ++n;
				}
				if (n < 8) continue;
				const double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
				const double cov = sxy - sx * sy / n;
				const double r = (vx > 0 && vy > 0)
				    ? cov / std::sqrt(vx * vy) : 1.0;
				const double alpha_corr =
				    std::min(1.0, std::max(0.0, (0.9 - r) / 0.5));
				const double alpha_struct =
				    std::min(1.0, std::max(0.0, 1.0 - 2.0 * sv[a] / med));
				const double alpha = std::max(alpha_corr, alpha_struct);
				if (alpha <= 0.0) continue;
				auto &dst = harmonised[a];
				for (int f = 0; f < N; ++f)
					dst[f] = (1.0 - alpha) * out.qualities[a][f]
					       + alpha * cons[f];
				any = true;
			}
		}
		if (any) rank_q = &harmonised;
	}

	/* Soft-selection taper (0 in explicit frame-number mode and whenever
	 * there is no room above the cut — e.g. the default 100%). */
	out.taper = stack_compute_selection_taper(*rank_q, aps, N,
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
		const auto &q = (*rank_q)[a];
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

/* Boundary-patch extension. The AP grid stops where the reference frame
 * has structure, so along an object's rim (e.g. the solar limb) the
 * coverage boundary — where Stage C hands over to the softer,
 * global-aligned background — runs right across the faint structure just
 * outside it (spicule layer, prominence fringes), and the hand-off band
 * follows the scalloped outer window contours, imprinting AP-pitch arcs
 * under sharpening. Push that boundary one grid step outward wherever no
 * other AP's patch covers the space beyond an edge: the extended area
 * stacks the AP's own frames at its own measured shifts (locally
 * aligned, unlike the background), and the hand-off moves into empty sky
 * where it is invisible. Interior edges (fully covered by a neighbour's
 * patch) are left untouched, so the mosaic's overlap topology on the
 * object is unchanged. Returns per-AP extended bounds in ORIGINAL frame
 * coordinates; shared by both Stage C engines. */
std::vector<cv::Vec4i> stack_extended_patch_bounds(const mpp_aps_t &aps,
                                                   int dim_y, int dim_x,
                                                   const mpp_config_t &cfg) {
	const int M = aps.count;
	const int ext_step = mpp_cfg_step_size(&cfg);
	std::vector<cv::Vec4i> patch(M), ext(M);
	for (int a = 0; a < M; ++a) {
		const auto &ap = aps.records[a];
		patch[a] = cv::Vec4i(ap.patch_y_low, ap.patch_y_high,
		                     ap.patch_x_low, ap.patch_x_high);
		ext[a] = patch[a];
	}
	/* Is the 1-px strip beyond `edge` (a row for horiz_strip, else a
	 * column), across [lo, hi) of the perpendicular axis, ≥98 % covered
	 * by the union of the OTHER patches? */
	auto covered = [&](int a, int lo, int hi, int edge, bool horiz_strip) {
		std::vector<std::pair<int, int>> iv;
		for (int b = 0; b < M; ++b) {
			if (b == a) continue;
			const auto &p = patch[b];
			const int b_lo  = horiz_strip ? p[2] : p[0];
			const int b_hi  = horiz_strip ? p[3] : p[1];
			const int b_plo = horiz_strip ? p[0] : p[2];
			const int b_phi = horiz_strip ? p[1] : p[3];
			if (b_plo > edge || b_phi < edge + 1) continue;
			const int s0 = std::max(lo, b_lo), s1 = std::min(hi, b_hi);
			if (s1 > s0) iv.emplace_back(s0, s1);
		}
		std::sort(iv.begin(), iv.end());
		long len = 0;
		int cur = lo;
		for (const auto &v : iv) {
			if (v.second <= cur) continue;
			len += v.second - std::max(cur, v.first);
			cur = std::max(cur, v.second);
		}
		return (double) len >= 0.98 * (double) (hi - lo);
	};
	for (int a = 0; a < M; ++a) {
		const auto &p = patch[a];
		if (p[0] > 0     && !covered(a, p[2], p[3], p[0] - 1, true))
			ext[a][0] = std::max(0, p[0] - ext_step);
		if (p[1] < dim_y && !covered(a, p[2], p[3], p[1], true))
			ext[a][1] = std::min(dim_y, p[1] + ext_step);
		if (p[2] > 0     && !covered(a, p[0], p[1], p[2] - 1, false))
			ext[a][2] = std::max(0, p[2] - ext_step);
		if (p[3] < dim_x && !covered(a, p[0], p[1], p[3], false))
			ext[a][3] = std::min(dim_x, p[3] + ext_step);
	}
	return ext;
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
	/* Sub-pixel correlation is always on: the warp engine samples each
	 * frame at the fractional displacement field, so an integer-only fit
	 * would quantise every AP to the pixel grid and re-introduce
	 * block-lattice phase noise. The parabolic refinement is a 3×3 fit on
	 * the already-computed correlation surface — negligible cost. */
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
	/* Escalation retry for failed local searches. A failure means the
	 * correlation peak railed at the search border — usually because the
	 * true local shift exceeds the window (large radial limb excursions
	 * in poor solar frames rail phase 1, whose {0,0} fallback then stacks
	 * the patch at global-only alignment and imprints AP-pitch beads
	 * along the limb). Re-measure just the failed pairs with a wider
	 * phase-1 window; only they pay the extra cost. +8 raises the
	 * phase-1 reach from ±(search_width−4) to ±(search_width+4) px. */
	mpp_config_t cfg_retry = cfg;
	cfg_retry.alignment_points_search_width = cfg.alignment_points_search_width + 8;
	const cv::Mat weight_matrix_retry = stack_build_first_phase_weight_matrix(cfg_retry);
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
	gint retried = 0, recovered = 0;
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
					MultilevelShiftResult r = multilevel_correlation(
					    ref_boxes[a].second_phase, ref_boxes[a].first_phase,
					    frame,
					    ap.box_y_low + dy, ap.box_y_high + dy,
					    ap.box_x_low + dx, ap.box_x_high + dx,
					    cfg.frames_gauss_width, cfg.alignment_points_search_width,
					    use_subpixel, weight_matrix);
					if (!r.success) {
						const MultilevelShiftResult r2 = multilevel_correlation(
						    ref_boxes[a].second_phase, ref_boxes[a].first_phase,
						    frame,
						    ap.box_y_low + dy, ap.box_y_high + dy,
						    ap.box_x_low + dx, ap.box_x_high + dx,
						    cfg.frames_gauss_width,
						    cfg_retry.alignment_points_search_width,
						    use_subpixel, weight_matrix_retry);
						/* Keep the original coarse estimate over a retry
						 * that produced none (railed again / window
						 * off-frame). */
						const bool r2_useless = !r2.success
						    && r2.dy == 0.0 && r2.dx == 0.0
						    && (r.dy != 0.0 || r.dx != 0.0);
						if (!r2_useless) r = r2;
						g_atomic_int_inc(&retried);
						if (r.success) g_atomic_int_inc(&recovered);
					}
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
	if (retried > 0)
		siril_log_message(_("Register: %d failed local searches retried at "
		                    "width %d, %d recovered\n"),
		                  g_atomic_int_get(&retried),
		                  cfg_retry.alignment_points_search_width,
		                  g_atomic_int_get(&recovered));
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

/* Signal ramp for the brightness-aware background blend:
 * s = clamp((v − lo) / lo, 0, 1), i.e. 0 at/below the background line
 * `lo`, full weight from 2·lo up. `v` is a single-channel CV_32F
 * luminance; returns CV_32F. */
cv::Mat dc_signal_ramp(const cv::Mat &v_gray, double lo) {
	cv::Mat s;
	v_gray.convertTo(s, CV_32F, 1.0 / lo, -1.0);
	cv::min(s, 1.0f, s);
	cv::max(s, 0.0f, s);
	return s;
}

/* Background line for the signal ramp, in pixel units. Reuses the
 * AP-placement brightness threshold — the existing "background vs
 * structure" knob — so ≤ 0 disables the ramped blend boost. */
double dc_signal_floor(const mpp_config_t &cfg) {
	return (double) cfg.alignment_points_brightness_threshold
	     * mpp_cfg_threshold_scale(&cfg);
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
