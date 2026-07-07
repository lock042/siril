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
 * Phase 4: per-AP local shifts via the default MultiLevelCorrelation
 * method. The two-phase cross-correlation kernel itself is shared with
 * Phase 2's global aligner via mpp::multilevel_correlation.
 *
 * Pipeline per (frame, AP):
 *   1. Pre-built reference boxes:
 *        second_phase = mean_frame[box_y_low:box_y_high, box_x_low:box_x_high]
 *                       .astype(float32)
 *        first_phase  = second_phase[::2, ::2]
 *      `mean_frame` here is the post-blur version used by Phase 3.
 *   2. Per-frame offset folds in both the intersection-vs-original-frame
 *      offset *and* the per-frame global shift from Phase 2. So the box
 *      bounds for the search become
 *        (y_low + dy, y_high + dy, x_low + dx, x_high + dx)
 *      in the frame's own coordinates.
 *   3. multilevel_correlation with alignment_points_search_width (14 by
 *      default — smaller than the 34 used for global alignment) and the
 *      shared frames_gauss_width (7). Optionally sub-pixel via cfg's
 *      alignment_points_local_search_subpixel.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_shift.h"
#include "registration/mpp/mpp_shift_priv.hpp"

namespace mpp {

namespace {

/* NumPy [::s, ::s] sub-sampling — duplicated from mpp_align.cpp to keep the
 * modules independent. */
cv::Mat stride_subsample_local(const cv::Mat &src, int s) {
	if (s <= 1) return src;
	const int new_rows = (src.rows + s - 1) / s;
	const int new_cols = (src.cols + s - 1) / s;
	cv::Mat out(new_rows, new_cols, src.type());
	const int elem_size = (int) src.elemSize();
	for (int y = 0; y < new_rows; ++y) {
		const uchar *src_row = src.ptr<uchar>(y * s);
		uchar *dst_row = out.ptr<uchar>(y);
		for (int x = 0; x < new_cols; ++x)
			std::memcpy(dst_row + x * elem_size,
			            src_row + (x * s) * elem_size, elem_size);
	}
	return out;
}

}  // namespace

std::vector<APRefBoxes> shift_prepare_ref_boxes(const cv::Mat &mean_frame,
                                                const mpp_aps_t &aps) {
	std::vector<APRefBoxes> out;
	out.reserve(aps.count);
	cv::Mat mean_f32;
	mean_frame.convertTo(mean_f32, CV_32F);
	for (int i = 0; i < aps.count; ++i) {
		const auto &ap = aps.records[i];
		APRefBoxes b;
		b.second_phase = mean_f32(cv::Range(ap.box_y_low, ap.box_y_high),
		                          cv::Range(ap.box_x_low, ap.box_x_high)).clone();
		b.first_phase = stride_subsample_local(b.second_phase, 2);
		out.push_back(std::move(b));
	}
	return out;
}

std::vector<FrameOffset> shift_frame_offsets(const std::vector<cv::Vec2d> &global_shifts,
                                             const cv::Vec4i &intersection) {
	std::vector<FrameOffset> out(global_shifts.size());
	for (size_t i = 0; i < global_shifts.size(); ++i) {
		out[i].dy = (double) intersection[0] - global_shifts[i][0];
		out[i].dx = (double) intersection[2] - global_shifts[i][1];
	}
	return out;
}

mpp_shifts_t *shift_compute_all(const std::vector<cv::Mat> &frames_blurred,
                                const mpp_aps_t &aps,
                                const std::vector<APRefBoxes> &ref_boxes,
                                const std::vector<FrameOffset> &offsets,
                                const mpp_config_t &cfg) {
	const int N = (int) frames_blurred.size();
	const int M = aps.count;
	if (N == 0 || M == 0 || (int) ref_boxes.size() != M
	 || (int) offsets.size() != N)
		return nullptr;

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

	const int gw = cfg.frames_gauss_width;
	const int sw = cfg.alignment_points_search_width;
	const bool subpix = cfg.alignment_points_local_search_subpixel;

	for (int f = 0; f < N; ++f) {
		/* See the long sign-trace comment in
		 * stack_compute_shifts_streamed — per-AP shift stored =
		 * r.dy + sub_y so that drizzle's gdy + weighted_ap equals
		 * the full sub-pixel alignment correction. */
		const double dy_full = offsets[f].dy;
		const double dx_full = offsets[f].dx;
		const int dy = (int) std::lround(dy_full);
		const int dx = (int) std::lround(dx_full);
		const double sub_y = dy_full - (double) dy;
		const double sub_x = dx_full - (double) dx;
		for (int a = 0; a < M; ++a) {
			const auto &ap = aps.records[a];
			const MultilevelShiftResult r = multilevel_correlation(
			    ref_boxes[a].second_phase, ref_boxes[a].first_phase,
			    frames_blurred[f],
			    ap.box_y_low  + dy, ap.box_y_high + dy,
			    ap.box_x_low  + dx, ap.box_x_high + dx,
			    gw, sw, subpix);
			const size_t s_off = (size_t) (f * M + a) * 2;
			out->shifts[s_off + 0] = r.dy + sub_y;
			out->shifts[s_off + 1] = r.dx + sub_x;
			out->success[f * M + a] = r.success ? 1 : 0;
		}
	}
	return out;
}

namespace {

/* Solve the symmetric 3×3 system A·p = rhs (A packed as a00 a01 a02 a11
 * a12 a22) by Gaussian elimination with partial pivoting. Returns false
 * on a (near-)singular system — collinear/degenerate neighbourhoods. */
bool solve_sym3(const double A[6], const double rhs[3], double p[3]) {
	double m[3][4] = {
		{A[0], A[1], A[2], rhs[0]},
		{A[1], A[3], A[4], rhs[1]},
		{A[2], A[4], A[5], rhs[2]},
	};
	for (int col = 0; col < 3; ++col) {
		int piv = col;
		for (int r = col + 1; r < 3; ++r)
			if (std::abs(m[r][col]) > std::abs(m[piv][col])) piv = r;
		if (std::abs(m[piv][col]) < 1e-9)
			return false;
		if (piv != col)
			for (int k = col; k < 4; ++k) std::swap(m[piv][k], m[col][k]);
		for (int r = col + 1; r < 3; ++r) {
			const double t = m[r][col] / m[col][col];
			for (int k = col; k < 4; ++k) m[r][k] -= t * m[col][k];
		}
	}
	for (int r = 2; r >= 0; --r) {
		double s = m[r][3];
		for (int k = r + 1; k < 3; ++k) s -= m[r][k] * p[k];
		p[r] = s / m[r][r];
	}
	return true;
}

}  // namespace

void shift_field_smooth(mpp_shifts_t *shifts, const mpp_aps_t &aps,
                        const mpp_config_t &cfg, const int *included,
                        int max_threads) {
	if (!shifts || !shifts->shifts || !shifts->success) return;
	if (cfg.alignment_points_smooth_radius <= 0.0) return;
	const int N = shifts->num_frames;
	const int M = shifts->num_aps;
	const int MIN_OK = 5;   /* successful neighbours needed to fit a plane */
	if (M != aps.count || M < MIN_OK || N <= 0) return;

	const double R = cfg.alignment_points_smooth_radius
	               * (double) mpp_cfg_step_size(&cfg);
	if (R <= 0.0) return;

	/* Static geometry, shared by every frame: per-AP neighbour lists
	 * (self included) with tricube distance weights and R-normalised
	 * position offsets (keeps the normal equations well-scaled). */
	struct Nbr { int b; double w_dist, dy, dx; };
	std::vector<std::vector<Nbr>> nbrs(M);
	for (int a = 0; a < M; ++a) {
		for (int b = 0; b < M; ++b) {
			const double dy = (double) (aps.records[b].y - aps.records[a].y);
			const double dx = (double) (aps.records[b].x - aps.records[a].x);
			const double d = std::hypot(dy, dx);
			if (d > R) continue;
			const double t = d / R;
			const double u = 1.0 - t * t * t;
			nbrs[a].push_back({b, u * u * u, dy / R, dx / R});
		}
	}

	const int nt = max_threads > 1 ? max_threads : 1;
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (int f = 0; f < N; ++f) {
		if (included && !included[f]) continue;
		const double *S = shifts->shifts + (size_t) f * M * 2;
		const uint8_t *ok = shifts->success + (size_t) f * M;

		/* Fit from the RAW values of the whole frame, then commit — an
		 * in-place update would feed already-smoothed neighbours into
		 * later APs' fits and bias the field toward the low-index side. */
		std::vector<double> out(2 * (size_t) M);
		std::vector<uint8_t> have(M, 0);
		std::vector<double> w, absr;

		for (int a = 0; a < M; ++a) {
			const auto &nb = nbrs[a];
			int n_ok = 0;
			for (const auto &v : nb) if (ok[v.b]) ++n_ok;
			if (n_ok < MIN_OK) continue;

			for (int c = 0; c < 2; ++c) {
				w.assign(nb.size(), 0.0);
				for (size_t i = 0; i < nb.size(); ++i)
					if (ok[nb[i].b]) w[i] = nb[i].w_dist;

				double p[3] = {0.0, 0.0, 0.0};
				for (int iter = 0; iter < 3; ++iter) {
					double A[6] = {0}, rhs[3] = {0};
					double wsum = 0.0, wval = 0.0;
					for (size_t i = 0; i < nb.size(); ++i) {
						const double wi = w[i];
						if (wi <= 0.0) continue;
						const double v  = S[(size_t) nb[i].b * 2 + c];
						const double b1 = nb[i].dy, b2 = nb[i].dx;
						A[0] += wi;           A[1] += wi * b1;      A[2] += wi * b2;
						A[3] += wi * b1 * b1; A[4] += wi * b1 * b2; A[5] += wi * b2 * b2;
						rhs[0] += wi * v; rhs[1] += wi * b1 * v; rhs[2] += wi * b2 * v;
						wsum += wi; wval += wi * v;
					}
					if (wsum <= 0.0) { p[0] = S[(size_t) a * 2 + c]; p[1] = p[2] = 0.0; break; }
					if (!solve_sym3(A, rhs, p)) {
						/* Degenerate neighbourhood (collinear APs):
						 * weighted-mean fallback still averages noise. */
						p[0] = wval / wsum; p[1] = p[2] = 0.0;
					}
					if (iter == 2) break;

					/* Tukey bisquare reweight on the residual scale
					 * (1.4826 × MAD, c = 4.685 — Siril's robust-fit
					 * convention, cf. gsl_multifit_robust_bisquare). */
					absr.clear();
					for (size_t i = 0; i < nb.size(); ++i) {
						if (!ok[nb[i].b]) continue;
						const double v = S[(size_t) nb[i].b * 2 + c];
						absr.push_back(std::abs(v - (p[0] + p[1] * nb[i].dy + p[2] * nb[i].dx)));
					}
					if (absr.empty()) break;
					const size_t h = absr.size() / 2;
					std::nth_element(absr.begin(), absr.begin() + h, absr.end());
					const double s = 1.4826 * absr[h];
					if (s < 1e-3) break;   /* already converged (px units) */
					const double cs = 4.685 * s;
					for (size_t i = 0; i < nb.size(); ++i) {
						if (!ok[nb[i].b]) { w[i] = 0.0; continue; }
						const double v = S[(size_t) nb[i].b * 2 + c];
						const double r = (v - (p[0] + p[1] * nb[i].dy + p[2] * nb[i].dx)) / cs;
						w[i] = (std::abs(r) < 1.0)
						     ? nb[i].w_dist * (1.0 - r * r) * (1.0 - r * r)
						     : 0.0;
					}
				}
				out[(size_t) a * 2 + c] = p[0];   /* fit at the AP centre */
			}
			have[a] = 1;
		}

		double *Sw = shifts->shifts + (size_t) f * M * 2;
		for (int a = 0; a < M; ++a) {
			if (!have[a]) continue;
			Sw[(size_t) a * 2 + 0] = out[(size_t) a * 2 + 0];
			Sw[(size_t) a * 2 + 1] = out[(size_t) a * 2 + 1];
		}
	}
}

}  // namespace mpp

/* ----------------- Public C interface ----------------- */

extern "C" mpp_status_t mpp_shift_compute(sequence *seq,
                                          const mpp_config_t *cfg,
                                          const mpp_aps_t *aps,
                                          const fits *ref,
                                          const int *global_shifts,
                                          mpp_shifts_t **shifts_out) {
	(void) seq; (void) cfg; (void) aps;
	(void) ref; (void) global_shifts;
	/* Sequence integration lands with mpp_frames; the algorithm is exercised
	 * via the C++ entry points (mpp_shift_priv.hpp). */
	if (shifts_out) *shifts_out = nullptr;
	return MPP_ENOTIMPL;
}

extern "C" void mpp_shift_free(mpp_shifts_t *shifts) {
	if (!shifts) return;
	if (shifts->shifts)  std::free(shifts->shifts);
	if (shifts->success) std::free(shifts->success);
	std::free(shifts);
}
