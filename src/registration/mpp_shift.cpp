/*
 * Phase 4: per-AP local shifts.
 *
 * Ports PSS's alignment_points.compute_shift_alignment_point for the default
 * MultiLevelCorrelation method (alignment_points.py:716). The two-phase
 * cross-correlation kernel itself is shared with Phase 2's global aligner
 * via mpp::multilevel_correlation (mpp_align_priv.hpp / mpp_align.cpp).
 *
 * Pipeline per (frame, AP):
 *   1. Pre-built reference boxes (PSS set_reference_boxes_correlation):
 *        second_phase = mean_frame[box_y_low:box_y_high, box_x_low:box_x_high]
 *                       .astype(float32)
 *        first_phase  = second_phase[::2, ::2]
 *      `mean_frame` here is the post-blur version used by Phase 3.
 *   2. Per-frame offset: align with PSS align_frames.dy / dx, which folds in
 *      both the intersection-vs-original-frame offset *and* the per-frame
 *      global shift from Phase 2. So the box bounds for the search become
 *        (y_low + dy, y_high + dy, x_low + dx, x_high + dx)
 *      in the frame's own coordinates.
 *   3. multilevel_correlation with alignment_points_search_width (14 by
 *      default — smaller than the 34 used for global alignment) and the
 *      shared frames_gauss_width (7). Optionally sub-pixel via cfg's
 *      alignment_points_local_search_subpixel.
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_shift.h"
#include "registration/mpp_shift_priv.hpp"

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
		 * stack_compute_shifts_streamed (mpp_stack.cpp) — per-AP shift
		 * stored = r.dy + sub_y so that drizzle's gdy + weighted_ap
		 * equals the full sub-pixel alignment correction. */
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
