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
 * mpp_drizzle — Phase 5b STScI drizzle backend for MPP stacking.
 *
 * Slice 5b.1 (the pixmap builder): encodes the full MPP warp as a per-
 * pixel output-coordinate map suitable for handing to dobox().
 *
 * Slice 5b.2 (this file's STScI stack path): drives dobox() one frame
 * at a time, accumulating into a shared float output canvas + counts
 * buffer, then normalises and casts to uint16.
 *
 * Slice 5b.3 (Bayer-drizzle) lands in a follow-up.
 *
 * NB: include ordering. driz_portability.h does `#define private` for C
 * compat, which wrecks subsequent C++ standard-library includes (their
 * class definitions use `private:` as a keyword). All C++ standard
 * headers and C++ helpers MUST be included before the drizzle headers.
 * The drizzle block is included last and bookended with `#undef private`.
 */

/* === C++ std + Siril's C++-side headers FIRST === */
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

/* registration/mpp.h transitively pulls in core/siril.h → <omp.h>, whose
 * libgomp 15 headers declare C++ templates that can't sit inside an
 * extern "C" block. Include core/siril.h first so its header guard
 * suppresses the re-include from inside the wrap. */
#include "core/siril.h"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_drizzle.h"
#include "registration/mpp/mpp_shift.h"

#include "core/proto.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "algos/demosaicing.h"
#include "io/image_format_fits.h"
#include "io/ser.h"
}

#include "registration/mpp/mpp_drizzle_priv.hpp"  /* mpp::stack_apply_stsci */
#include "registration/mpp/mpp_stack_priv.hpp"    /* mpp::stack_one_dim_weight */

/* read_full_frame is in mpp.cpp's mpp:: namespace (moved out of the
 * inner anonymous namespace in slice 5b.2). Forward-declared here so we
 * can read frames without exposing it through a public header. The
 * trailing avi_pattern arg (default MPP_AVI_BAYER_AUTO) lets callers
 * trigger an OpenCV debayer on AVI Bayer mosaics. avi_bayer_to_cfa is
 * declared here too for the AVI Bayer-drizzle wrapper below. */
namespace mpp {
extern cv::Mat read_full_frame(sequence *seq, int idx, int avi_pattern);
extern bool avi_bayer_to_cfa(int avi_pattern, int rows, unsigned char cfa[4]);
}

/* === Drizzle headers LAST (they redefine `private`) === */
/* cdrizzlemap.h is deliberately NOT included — it pulls in wcslib.h
 * which has a conflicting `int64` typedef vs OpenCV's. cdrizzlebox.h
 * (with dobox + driz_param_t + e_kernel_t via cdrizzleutil.h) is all
 * we need for the stack path. */
extern "C" {
#include "drizzle/cdrizzleutil.h"
#include "drizzle/cdrizzlebox.h"
}
#undef private


/* ============================================================
 * Pixmap (slice 5b.1) — unchanged in this commit.
 * ============================================================ */

extern "C" mpp_status_t mpp_imgmap_alloc(imgmap_t *out, int rx, int ry) {
	if (!out || rx <= 0 || ry <= 0) return MPP_EINVAL;
	const size_t n = (size_t) rx * (size_t) ry;
	out->xmap = (float *) std::malloc(n * sizeof(float));
	if (!out->xmap) { out->ymap = nullptr; out->rx = 0; out->ry = 0; return MPP_ENOMEM; }
	out->ymap = (float *) std::malloc(n * sizeof(float));
	if (!out->ymap) {
		std::free(out->xmap);
		out->xmap = nullptr;
		out->rx = 0; out->ry = 0;
		return MPP_ENOMEM;
	}
	out->rx = rx;
	out->ry = ry;
	return MPP_OK;
}

extern "C" void mpp_imgmap_free(imgmap_t *m) {
	if (!m) return;
	std::free(m->xmap);
	std::free(m->ymap);
	m->xmap = nullptr;
	m->ymap = nullptr;
	m->rx = 0;
	m->ry = 0;
}

extern "C" mpp_status_t mpp_pixmap_build_filtered(const mpp_run_t *run,
                                                  int frame_idx,
                                                  double drizzle_scale,
                                                  const int *ap_active,
                                                  int include_background,
                                                  imgmap_t *out) {
	if (!run || !out || !out->xmap || !out->ymap) return MPP_EINVAL;
	if (!run->aps || !run->global_shifts || !run->shifts || !run->shifts->shifts)
		return MPP_EINVAL;
	if (frame_idx < 0 || frame_idx >= run->num_frames) return MPP_EINVAL;
	if (frame_idx >= run->shifts->num_frames) return MPP_EINVAL;
	if (out->rx != run->frame_cols || out->ry != run->frame_rows) return MPP_EINVAL;
	if (drizzle_scale <= 0.0) return MPP_EINVAL;

	const int rx = out->rx;
	const int ry = out->ry;
	const int M  = run->aps->count;

	const double gdy = run->global_shifts[2 * frame_idx + 0];
	const double gdx = run->global_shifts[2 * frame_idx + 1];

	const double iy_lo = (double) run->intersection[0];
	const double ix_lo = (double) run->intersection[2];

	std::vector<std::vector<float>> wy(M), wx(M);
	for (int a = 0; a < M; ++a) {
		const auto &ap = run->aps->records[a];
		const bool extend_y_low  = (ap.patch_y_low  == 0);
		const bool extend_y_high = (ap.patch_y_high == ry);
		const bool extend_x_low  = (ap.patch_x_low  == 0);
		const bool extend_x_high = (ap.patch_x_high == rx);
		wy[a] = mpp::stack_one_dim_weight(ap.patch_y_low, ap.patch_y_high,
		                                  ap.y, extend_y_low, extend_y_high);
		wx[a] = mpp::stack_one_dim_weight(ap.patch_x_low, ap.patch_x_high,
		                                  ap.x, extend_x_low, extend_x_high);
	}

	const size_t frame_base = (size_t) frame_idx * (size_t) M;
	const double *ap_shifts = &run->shifts->shifts[2 * frame_base];
	/* NaN sentinel makes dobox's map_pixel (cdrizzlemap.c:381) return
	 * "out of bounds" so the input pixel is skipped — this is how per-AP
	 * top-K filtering selectively excludes a frame from regions whose
	 * covering APs don't have it in their top-stack_size. */
	const float skip = std::nanf("");

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 0; j < ry; ++j) {
		for (int i = 0; i < rx; ++i) {
			double wsum = 0.0, wsy = 0.0, wsx = 0.0;
			bool any_ap_covers = false;
			for (int a = 0; a < M; ++a) {
				const auto &ap = run->aps->records[a];
				if (j < ap.patch_y_low || j >= ap.patch_y_high) continue;
				if (i < ap.patch_x_low || i >= ap.patch_x_high) continue;
				any_ap_covers = true;
				if (ap_active && !ap_active[a]) continue;
				const float wy_v = wy[a][j - ap.patch_y_low];
				const float wx_v = wx[a][i - ap.patch_x_low];
				const double w   = (double) std::min(wy_v, wx_v);
				if (w <= 0.0) continue;
				wsum += w;
				wsy  += w * ap_shifts[2 * a + 0];
				wsx  += w * ap_shifts[2 * a + 1];
			}
			const size_t idx = (size_t) j * (size_t) rx + (size_t) i;
			if (wsum > 0.0) {
				/* Pixel falls inside ≥1 AP that has this frame in its
				 * top-K — accumulate with the (filtered) AP-weighted
				 * shift. */
				const double sy = gdy + wsy / wsum;
				const double sx = gdx + wsx / wsum;
				out->xmap[idx] = (float) (drizzle_scale * ((double) i + sx - ix_lo));
				out->ymap[idx] = (float) (drizzle_scale * ((double) j + sy - iy_lo));
			} else if (!any_ap_covers && include_background) {
				/* True background region (no AP covers) and this frame
				 * is allowed to contribute to background — use the
				 * global shift only. */
				out->xmap[idx] = (float) (drizzle_scale * ((double) i + gdx - ix_lo));
				out->ymap[idx] = (float) (drizzle_scale * ((double) j + gdy - iy_lo));
			} else {
				/* Covered only by APs that don't have this frame in
				 * their top-K, or background but frame is below the
				 * overall top-stack_size threshold — skip. */
				out->xmap[idx] = skip;
				out->ymap[idx] = skip;
			}
		}
	}
	return MPP_OK;
}

extern "C" mpp_status_t mpp_pixmap_build(const mpp_run_t *run,
                                         int frame_idx,
                                         double drizzle_scale,
                                         imgmap_t *out) {
	return mpp_pixmap_build_filtered(run, frame_idx, drizzle_scale,
	                                 /*ap_active=*/nullptr,
	                                 /*include_background=*/1, out);
}


/* ============================================================
 * STScI stack path (slice 5b.2).
 * ============================================================ */

namespace {

enum e_kernel_t kernel_from_cfg(int mpp_kernel) {
	switch (mpp_kernel) {
		case MPP_KERNEL_GAUSSIAN: return kernel_gaussian;
		case MPP_KERNEL_POINT:    return kernel_point;
		case MPP_KERNEL_TURBO:    return kernel_turbo;
		case MPP_KERNEL_LANCZOS2: return kernel_lanczos2;
		case MPP_KERNEL_LANCZOS3: return kernel_lanczos3;
		case MPP_KERNEL_SQUARE:
		default:                  return kernel_square;
	}
}

/* Allocate a float-typed `fits` with the given shape, layout planar
 * (R then G then B for multi-channel) and pdata pointers set up for
 * dobox's get_pixel accessors. Returns false on allocation failure. */
bool alloc_float_fits(fits *out, int rx, int ry, int channels) {
	std::memset(out, 0, sizeof(fits));
	out->rx = rx;
	out->ry = ry;
	out->naxes[0] = rx;
	out->naxes[1] = ry;
	out->naxes[2] = channels;
	out->naxis = (channels == 1) ? 2 : 3;
	out->bitpix = FLOAT_IMG;
	out->orig_bitpix = FLOAT_IMG;
	out->type = DATA_FLOAT;
	const size_t plane = (size_t) rx * (size_t) ry;
	out->fdata = (float *) std::calloc(plane * (size_t) channels, sizeof(float));
	if (!out->fdata) return false;
	out->fpdata[0] = out->fdata;
	out->fpdata[1] = (channels >= 2) ? out->fdata + plane     : nullptr;
	out->fpdata[2] = (channels >= 3) ? out->fdata + 2 * plane : nullptr;
	return true;
}

/* Build a planar single-frame `fits` from a cv::Mat. The Mat may be
 * single-channel or 3-channel (interleaved BGR per OpenCV convention);
 * the fits is always planar. Applies the brightness-equalisation scale
 * during conversion. Caller owns `out->fdata` (clearfits to free). */
bool cv_to_planar_fits(const cv::Mat &raw, double scale, fits *out) {
	const int rx = raw.cols, ry = raw.rows, channels = raw.channels();
	if (!alloc_float_fits(out, rx, ry, channels)) return false;
	const size_t plane = (size_t) rx * (size_t) ry;
	if (channels == 1) {
		cv::Mat as_float;
		raw.convertTo(as_float, CV_32F, scale);
		std::memcpy(out->fdata, as_float.data, plane * sizeof(float));
	} else {
		/* Split into R/G/B planes. OpenCV stores BGR by default, but
		 * mpp::read_full_frame returns frames in Siril's R/G/B order
		 * (see mpp.cpp read_full_frame's cv::merge of pdata[0..2]). */
		cv::Mat as_float;
		raw.convertTo(as_float, CV_32FC3, scale);
		std::vector<cv::Mat> planes;
		cv::split(as_float, planes);
		for (int c = 0; c < channels; ++c)
			std::memcpy(out->fdata + (size_t) c * plane,
			            planes[c].data, plane * sizeof(float));
	}
	return true;
}

/* Per-AP top-K filter state. Built once per Stage C call, queried per
 * frame to derive (ap_active mask, include_background flag) for
 * mpp_pixmap_build_filtered. Mirrors the bicubic stack path's
 * apq.used_alignment_points / top_for_bg construction (mpp_stack.cpp:
 * stack_apply_shifts_streamed). Frames that contribute nowhere
 * (no AP top-K AND not in background top-K) can be skipped entirely. */
struct FrameFilterSetup {
	bool apply;                              /* false → no filtering (pre-filter behaviour) */
	int M;
	std::vector<std::vector<int>> used_alignment_points;
	std::vector<int> top_for_bg_flag;        /* length N, 0/1 */
	std::vector<int> ap_active_buf;          /* length M, reused per frame */
};

FrameFilterSetup make_frame_filter(const mpp_run_t *run, int N,
                                    const std::vector<int> &included) {
	FrameFilterSetup s;
	s.M     = run->aps ? run->aps->count : 0;
	s.apply = (run->stack_size > 0
	        && run->best_frame_indices != nullptr
	        && run->quality != nullptr
	        && s.M > 0);
	if (!s.apply) return s;

	s.used_alignment_points.assign(N, {});
	for (int a = 0; a < s.M; ++a) {
		for (int k = 0; k < run->stack_size; ++k) {
			const int idx = run->best_frame_indices[a * run->stack_size + k];
			if (idx >= 0 && idx < N)
				s.used_alignment_points[idx].push_back(a);
		}
	}

	/* Top-stack_size frames globally for the background. Restricted to
	 * the included set so an excluded-but-high-quality frame doesn't
	 * displace an included one. */
	std::vector<int> sorted_idx;
	sorted_idx.reserve(N);
	for (int i = 0; i < N; ++i)
		if (included[i]) sorted_idx.push_back(i);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [run](int a, int b) { return run->quality[a] > run->quality[b]; });
	s.top_for_bg_flag.assign(N, 0);
	const int K_bg = std::min((int) sorted_idx.size(), run->stack_size);
	for (int i = 0; i < K_bg; ++i)
		s.top_for_bg_flag[sorted_idx[i]] = 1;

	s.ap_active_buf.assign(s.M, 0);
	return s;
}

/* Populate ap_active_buf for frame f. Returns true if this frame should
 * be drizzled (contributes to at least one AP or to the background);
 * false to skip the disk read + dobox call entirely. */
bool prepare_frame_filter(FrameFilterSetup &s, int f,
                          const int **ap_active_ptr_out,
                          int *include_bg_flag_out) {
	if (!s.apply) {
		*ap_active_ptr_out   = nullptr;
		*include_bg_flag_out = 1;
		return true;
	}
	std::fill(s.ap_active_buf.begin(), s.ap_active_buf.end(), 0);
	for (int a : s.used_alignment_points[f]) s.ap_active_buf[a] = 1;
	const bool any_ap = !s.used_alignment_points[f].empty();
	const bool any_bg = s.top_for_bg_flag[f] != 0;
	*ap_active_ptr_out   = s.ap_active_buf.data();
	*include_bg_flag_out = any_bg ? 1 : 0;
	return any_ap || any_bg;
}

/* Per-frame drizzle iteration plan. dobox uses [xmin,xmax]×[ymin,ymax]
 * inclusive bounds. */
struct DrizRect {
	int xmin, ymin, xmax, ymax;
};

/* Strategy picker for per-frame drizzle iteration. Two candidates:
 *
 *   • bbox: one dobox call covering the bounding box of active patches.
 *     NaN entries in the pixmap (already set by mpp_pixmap_build_filtered)
 *     skip inactive pixels inside the bbox. Constant per-frame overhead;
 *     pixel visits = bbox area.
 *
 *   • per-AP: one dobox call per active AP, restricted to that AP's
 *     patch rect. Between calls we NaN-out the pixmap inside the
 *     processed rect (claim mask) so overlap with subsequent APs is
 *     visited once. Pixel visits = union of active patches; per-AP
 *     overhead scales with K (number of active APs).
 *
 * Heuristic: per-AP wins on pixel visits when Σ(patch_area) < bbox_area.
 * When patches tile densely (typical of frames where most APs are
 * active) the sum saturates near the bbox and bbox is at least as
 * cheap — falling back guarantees worst-case ≤ bbox-only behaviour.
 *
 * Background contribution (include_bg_flag = 1) forces bbox = full
 * frame because true-background pixels can fall anywhere outside the
 * APs; per-AP iteration would skip them.
 *
 * out_rects/out_needs_claim are caller-owned to avoid per-frame vector
 * allocation across 10⁴-10⁵ frames. */
void compute_driz_plan(const mpp_run_t *run,
                        const int *ap_active,
                        int include_bg_flag,
                        int frame_rx, int frame_ry,
                        std::vector<DrizRect> &out_rects,
                        bool &out_needs_claim) {
	out_rects.clear();
	out_needs_claim = false;

	const DrizRect full = {0, 0, frame_rx - 1, frame_ry - 1};
	if (include_bg_flag || !ap_active || !run->aps || run->aps->count <= 0) {
		out_rects.push_back(full);
		return;
	}

	const int M = run->aps->count;
	int K = 0;
	int bbox_xmin = INT_MAX, bbox_ymin = INT_MAX;
	int bbox_xmax = INT_MIN, bbox_ymax = INT_MIN;
	long long sum_patch_area = 0;
	for (int a = 0; a < M; ++a) {
		if (!ap_active[a]) continue;
		const auto &ap = run->aps->records[a];
		const int px_lo = std::max(0, ap.patch_x_low);
		const int px_hi = std::min(frame_rx, ap.patch_x_high);
		const int py_lo = std::max(0, ap.patch_y_low);
		const int py_hi = std::min(frame_ry, ap.patch_y_high);
		if (px_hi <= px_lo || py_hi <= py_lo) continue;
		++K;
		bbox_xmin = std::min(bbox_xmin, px_lo);
		bbox_xmax = std::max(bbox_xmax, px_hi - 1);
		bbox_ymin = std::min(bbox_ymin, py_lo);
		bbox_ymax = std::max(bbox_ymax, py_hi - 1);
		sum_patch_area += (long long)(px_hi - px_lo) * (long long)(py_hi - py_lo);
	}
	if (K == 0) {
		/* Shouldn't reach here — prepare_frame_filter would have
		 * returned false. Defensive fallback. */
		out_rects.push_back(full);
		return;
	}

	const long long bbox_area = (long long)(bbox_xmax - bbox_xmin + 1)
	                          * (long long)(bbox_ymax - bbox_ymin + 1);
	if (sum_patch_area < bbox_area) {
		out_rects.reserve(K);
		for (int a = 0; a < M; ++a) {
			if (!ap_active[a]) continue;
			const auto &ap = run->aps->records[a];
			const int px_lo = std::max(0, ap.patch_x_low);
			const int px_hi = std::min(frame_rx, ap.patch_x_high);
			const int py_lo = std::max(0, ap.patch_y_low);
			const int py_hi = std::min(frame_ry, ap.patch_y_high);
			if (px_hi <= px_lo || py_hi <= py_lo) continue;
			out_rects.push_back({px_lo, py_lo, px_hi - 1, py_hi - 1});
		}
		out_needs_claim = true;
	} else {
		out_rects.push_back({bbox_xmin, bbox_ymin, bbox_xmax, bbox_ymax});
	}
}

/* Mark a rectangular region of the pixmap as already-drizzled by NaN-ing
 * both xmap and ymap entries. Subsequent dobox calls' map_pixel checks
 * (cdrizzlemap.c:381) will treat the pixel as out-of-bounds and skip
 * it, so input pixels in patch overlap regions are accumulated exactly
 * once with the pre-claim AP-weighted shift. */
void pixmap_claim_rect(imgmap_t *pixmap, const DrizRect &r) {
	const float skip = std::nanf("");
	for (int j = r.ymin; j <= r.ymax; ++j) {
		const size_t row = (size_t) j * (size_t) pixmap->rx;
		for (int i = r.xmin; i <= r.xmax; ++i) {
			pixmap->xmap[row + i] = skip;
			pixmap->ymap[row + i] = skip;
		}
	}
}

}  // namespace

/* Inner C++ entry point. Tests drive this directly with vector<cv::Mat>
 * so they can run on PNG-extracted frames without constructing a Siril
 * sequence. The extern "C" mpp_stack_apply_stsci wrapper just reads
 * frames via mpp::read_full_frame and forwards. */
mpp_status_t mpp::stack_apply_stsci_streamed(const FrameProvider &provider,
                                             int N,
                                             const std::vector<int> &included,
                                             const std::vector<double> &frame_brightness,
                                             const mpp_run_t *run,
                                             const mpp_config_t *cfg,
                                             fits *out,
                                             int max_threads,
                                             bool provider_thread_safe) {
	if (!cfg || !run || !out || !run->aps || !run->shifts) return MPP_EINVAL;
	if (N <= 0 || (int) included.size() != N
	 || (int) frame_brightness.size() != N) return MPP_EINVAL;
	if (N != run->num_frames) return MPP_EINVAL;

	const int frame_rx = run->frame_cols;
	const int frame_ry = run->frame_rows;
	if (frame_rx <= 0 || frame_ry <= 0) return MPP_EINVAL;

	const int intersection_rx = run->intersection[3] - run->intersection[2];
	const int intersection_ry = run->intersection[1] - run->intersection[0];
	if (intersection_rx <= 0 || intersection_ry <= 0) return MPP_EINVAL;

	const double scale = cfg->drizzle_scale < 1.0 ? 1.0 : cfg->drizzle_scale;
	const int out_rx   = (int) std::lround((double) intersection_rx * scale);
	const int out_ry   = (int) std::lround((double) intersection_ry * scale);
	const int channels = run->num_layers > 0 ? run->num_layers : 1;

	fits output_data{};
	fits output_counts{};
	if (!alloc_float_fits(&output_data, out_rx, out_ry, channels)) return MPP_ENOMEM;
	if (!alloc_float_fits(&output_counts, out_rx, out_ry, channels)) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}

	std::vector<double> bright_inc;
	bright_inc.reserve(N);
	for (int i = 0; i < N; ++i)
		if (included[i]) bright_inc.push_back(frame_brightness[i]);
	if (bright_inc.empty()) {
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENODATA;
	}
	std::nth_element(bright_inc.begin(),
	                 bright_inc.begin() + bright_inc.size() / 2,
	                 bright_inc.end());
	const double median_brightness = bright_inc[bright_inc.size() / 2];

	FrameFilterSetup filter = make_frame_filter(run, N, included);

	struct driz_args_t driz_args;
	std::memset(&driz_args, 0, sizeof(driz_args));
	driz_args.kernel         = kernel_from_cfg(cfg->drizzle_kernel);
	driz_args.scale          = (float) scale;
	driz_args.pixel_fraction = (float) cfg->drizzle_pixfrac;
	driz_args.weight_scale   = 1.0f;
	driz_args.is_bayer       = FALSE;
	driz_args.cfadim         = 0;
	driz_args.use_flats      = FALSE;
	driz_args.flat           = nullptr;

	struct driz_error_t error;
	std::memset(&error, 0, sizeof(error));

	/* Frames that participate, in index order. */
	std::vector<int> work;
	work.reserve(N);
	for (int f = 0; f < N; ++f)
		if (included[f]) work.push_back(f);
	const int W = (int) work.size();
	const int n_included = W;

	/* dobox's do_kernel_square hardcodes get_pixel(p->data, i, j, 0) — it
	 * always reads channel 0 of the input fits, so multi-channel input is
	 * looped per channel with single-channel view fits pointing at the
	 * matching planes. For 1-channel input this is a single iteration. */
	const size_t in_plane  = (size_t) frame_rx * (size_t) frame_ry;
	const size_t out_plane = (size_t) out_rx * (size_t) out_ry;

	/* Frame-parallel drizzle. dobox (the per-frame kernel splat) is ~95 % of
	 * the runtime and can't be parallelised in place — it accumulates a
	 * running weighted mean into a shared canvas. So each thread drizzles
	 * its frames into a PRIVATE output canvas, and the canvases are then
	 * reduced as a weighted mean:  combined = Σ(meanₜ·countsₜ) / Σ countsₜ.
	 * Mathematically identical to the single-thread mean; floating-point
	 * reordering makes it differ at the sub-quantisation level, so the
	 * output is NOT bit-reproducible across thread counts — an accepted
	 * trade for parallelising the dominant cost of a statistical stack.
	 *
	 * Per-thread memory = two float output canvases + one in-flight frame
	 * (planar float + pixmap + decoded raw). Clamp the thread count so the
	 * private canvases fit the available-memory budget (also bounded by the
	 * user's configured memory cap). */
	const int max_t = max_threads > 0 ? max_threads : 1;
	int n_threads = max_t;
#ifdef _OPENMP
	if (provider_thread_safe && max_t > 1) {
		const size_t per_thread_canvas = (size_t) 2 * out_plane * (size_t) channels * sizeof(float);
		const size_t per_thread_frame  = (size_t) channels * in_plane * sizeof(float)  /* planar float */
		                               + (size_t) 2 * in_plane * sizeof(float)         /* pixmap xmap+ymap */
		                               + (size_t) channels * in_plane * 2;             /* decoded raw (~16-bit) */
		const size_t per_thread = per_thread_canvas + per_thread_frame;
		const size_t shared     = (size_t) 2 * out_plane * (size_t) channels * sizeof(float); /* reduction targets */
		guint64 budget_b = get_available_memory();
		const long long cap_mb = get_max_memory_in_MB();
		if (cap_mb > 0) {
			const guint64 cap_b = (guint64) cap_mb * 1024ULL * 1024ULL;
			if (cap_b < budget_b) budget_b = cap_b;
		}
		long double budget = (long double) budget_b * 0.8L;   /* 20 % headroom */
		budget = (budget > (long double) shared) ? budget - (long double) shared : 0.0L;
		int fit = (per_thread > 0) ? (int) (budget / (long double) per_thread) : max_t;
		if (fit < 1) fit = 1;
		if (fit < n_threads) {
			siril_log_message(_("Stack (mpp): STScI drizzle limited to %d of %d "
			                    "thread(s) by available memory\n"), fit, max_t);
			n_threads = fit;
		}
	} else {
		n_threads = 1;   /* non-reentrant provider → serial reads */
	}
#else
	n_threads = 1;
#endif

	gui_iface.set_progress(PROGRESS_RESET, _("Stack (STScI): drizzling frames"));
	gint cur_nb = 0;
	gint cancelled = 0;
	mpp_status_t rc = MPP_OK;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
#endif
	{
		/* Per-thread private output canvas + counts (calloc-zeroed), a
		 * reusable pixmap, and dobox error scratch. */
		fits od_local{}, oc_local{};
		imgmap_t pixmap_local{};
		struct driz_error_t error;
		std::memset(&error, 0, sizeof(error));
		const bool alloc_ok =
		    alloc_float_fits(&od_local, out_rx, out_ry, channels)
		    && alloc_float_fits(&oc_local, out_rx, out_ry, channels)
		    && (mpp_imgmap_alloc(&pixmap_local, frame_rx, frame_ry) == MPP_OK);
		mpp_status_t local_rc = alloc_ok ? MPP_OK : MPP_ENOMEM;

#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
		for (int wi = 0; wi < W; ++wi) {
			if (local_rc != MPP_OK || g_atomic_int_get(&cancelled)) continue;
			if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
			const int f = work[wi];

			/* Per-frame AP-active mask (thread-local; prepare_frame_filter's
			 * shared scratch buffer is not reentrant). */
			const int *ap_active_ptr = nullptr;
			int include_bg_flag = 1;
			std::vector<int> ap_active_local;
			if (filter.apply) {
				const bool any_ap = !filter.used_alignment_points[f].empty();
				const bool any_bg = filter.top_for_bg_flag[f] != 0;
				if (!(any_ap || any_bg)) {     /* contributes nowhere */
					gui_iface.set_progress((double) g_atomic_int_add(&cur_nb, 1) / (double) n_included, NULL);
					continue;
				}
				ap_active_local.assign(filter.M, 0);
				for (int a : filter.used_alignment_points[f]) ap_active_local[a] = 1;
				ap_active_ptr   = ap_active_local.data();
				include_bg_flag = any_bg ? 1 : 0;
			}

			const cv::Mat frame_raw = provider(f);
			if (frame_raw.empty()) {
				gui_iface.set_progress((double) g_atomic_int_add(&cur_nb, 1) / (double) n_included, NULL);
				continue;
			}
			const double brightness_scale = median_brightness / (frame_brightness[f] + 1e-7);
			fits frame_fit{};
			if (!cv_to_planar_fits(frame_raw, brightness_scale, &frame_fit)) {
				local_rc = MPP_ENOMEM; continue;
			}
			if (mpp_pixmap_build_filtered(run, f, scale, ap_active_ptr,
			                              include_bg_flag, &pixmap_local) != MPP_OK) {
				clearfits(&frame_fit); local_rc = MPP_EINVAL; continue;
			}
			std::vector<DrizRect> rects;
			bool needs_claim = false;
			compute_driz_plan(run, ap_active_ptr, include_bg_flag,
			                  frame_rx, frame_ry, rects, needs_claim);

			bool dobox_failed = false;
			int failed_channel = 0;
			for (size_t r = 0; r < rects.size() && !dobox_failed; ++r) {
				const DrizRect &rect = rects[r];
				for (int c = 0; c < channels; ++c) {
					fits frame_view = frame_fit;
					frame_view.fdata = frame_fit.fdata + (size_t) c * in_plane;
					frame_view.fpdata[0] = frame_view.fdata;
					frame_view.fpdata[1] = nullptr;
					frame_view.fpdata[2] = nullptr;
					frame_view.naxes[2]  = 1;
					frame_view.naxis     = 2;

					fits od_view = od_local;
					od_view.fdata = od_local.fdata + (size_t) c * out_plane;
					od_view.fpdata[0] = od_view.fdata;
					od_view.fpdata[1] = nullptr;
					od_view.fpdata[2] = nullptr;
					od_view.naxes[2]  = 1;
					od_view.naxis     = 2;

					fits oc_view = oc_local;
					oc_view.fdata = oc_local.fdata + (size_t) c * out_plane;
					oc_view.fpdata[0] = oc_view.fdata;
					oc_view.fpdata[1] = nullptr;
					oc_view.fpdata[2] = nullptr;
					oc_view.naxes[2]  = 1;
					oc_view.naxis     = 2;

					struct driz_param_t p;
					std::memset(&p, 0, sizeof(p));
					driz_param_init(&p);
					p.driz           = &driz_args;
					p.kernel         = driz_args.kernel;
					/* p.scale = 1.0 (not driz_args.scale): dobox uses it only
					 * for a scale² flux multiplier; the output geometry is
					 * carried by the pixmap. */
					p.scale          = 1.0f;
					p.pixel_fraction = driz_args.pixel_fraction;
					p.weight_scale   = 1.0f;
					p.in_units       = unit_counts;
					p.out_units      = unit_counts;
					p.exposure_time  = 1.0f;
					p.data           = &frame_view;
					p.weights        = nullptr;
					p.pixmap         = &pixmap_local;
					p.output_data    = &od_view;
					p.output_counts  = &oc_view;
					p.xmin           = rect.xmin;
					p.ymin           = rect.ymin;
					p.xmax           = rect.xmax;
					p.ymax           = rect.ymax;
					p.threads        = 1;
					p.error          = &error;

					if (dobox(&p)) {
						dobox_failed = true;
						failed_channel = c;
						break;
					}
				}
				/* In per-AP mode, NaN-out the just-processed patch so the
				 * next AP's dobox call skips overlap pixels. */
				if (!dobox_failed && needs_claim)
					pixmap_claim_rect(&pixmap_local, rect);
			}
			clearfits(&frame_fit);
			if (dobox_failed) {
				siril_log_error(
				    _("Stack (Drizzle): dobox failed on frame %d channel %d: %s\n"),
				    f, failed_channel,
				    error.last_message[0] ? error.last_message : _("(no detail)"));
				local_rc = MPP_EIO;
				continue;
			}
			gui_iface.set_progress((double) g_atomic_int_add(&cur_nb, 1) / (double) n_included, NULL);
		}

		/* Weighted-mean reduction into the shared canvases:
		 *   output_data   += meanₜ · countsₜ
		 *   output_counts += countsₜ
		 * The final divide (after the region) recovers the global mean. */
#ifdef _OPENMP
#pragma omp critical(stsci_reduce)
#endif
		{
			if (local_rc != MPP_OK) {
				if (rc == MPP_OK) rc = local_rc;
			} else {
				const size_t total = out_plane * (size_t) channels;
				for (size_t i = 0; i < total; ++i) {
					output_data.fdata[i]   += od_local.fdata[i] * oc_local.fdata[i];
					output_counts.fdata[i] += oc_local.fdata[i];
				}
			}
		}
		clearfits(&od_local);
		clearfits(&oc_local);
		mpp_imgmap_free(&pixmap_local);
	}

	if (cancelled && rc == MPP_OK) rc = MPP_EINTR;
	if (rc != MPP_OK) {
		clearfits(&output_data);
		clearfits(&output_counts);
		return rc;
	}

	/* Combine: output_data holds Σ(meanₜ·countsₜ); divide by the total
	 * weight to recover the global weighted mean. Cells no frame covered
	 * have zero weight and stay zero. The downstream cast reads output_data
	 * directly — it is already the mean, do NOT divide again. */
	{
		const size_t total = out_plane * (size_t) channels;
		for (size_t i = 0; i < total; ++i) {
			const float w = output_counts.fdata[i];
			output_data.fdata[i] = (w > 0.0f) ? (output_data.fdata[i] / w) : 0.0f;
		}
	}
	clearfits(&output_counts);

	clearfits(out);
	out->rx = out_rx;
	out->ry = out_ry;
	out->naxes[0] = out_rx;
	out->naxes[1] = out_ry;
	out->naxes[2] = channels;
	out->naxis = (channels == 1) ? 2 : 3;
	out->bitpix = USHORT_IMG;
	out->orig_bitpix = USHORT_IMG;
	out->type = DATA_USHORT;
	const size_t plane = (size_t) out_rx * (size_t) out_ry;
	out->data = (WORD *) std::calloc(plane * (size_t) channels, sizeof(WORD));
	if (!out->data) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}
	out->pdata[0] = out->data;
	/* Siril's mono convention: pdata[1]/pdata[2] alias pdata[0] (not
	 * NULL). Various GUI paths — notably remap_all_vports's per-vport
	 * memcpy — iterate c=0..2 unconditionally, so a NULL plane pointer
	 * crashes during display refresh after the stack swap. */
	out->pdata[1] = (channels >= 2) ? out->data + plane     : out->data;
	out->pdata[2] = (channels >= 3) ? out->data + 2 * plane : out->data;

	/* Scale to 16-bit uint range. For 8-bit input, output_data values
	 * are still in [0, 255] (the raw input range we fed dobox) so we
	 * multiply by 256 here — mean_frame is always carried in the 16-bit
	 * equivalent regardless of input bitdepth. For 16-bit input the
	 * scale is 1.0 and values pass through unchanged.
	 * NB: mpp_cfg_threshold_scale is the INVERSE direction (scales
	 * thresholds down to match input units), so we compute the output
	 * scale inline rather than reusing it. */
	const double cast_scale = (cfg->bitdepth == 8) ? 256.0 : 1.0;
	for (size_t k = 0; k < plane * (size_t) channels; ++k) {
		const double v = (double) output_data.fdata[k] * cast_scale;
		out->data[k] = (WORD) (v < 0.0 ? 0
		                     : (v > 65535.0 ? 65535
		                                    : (int) (v + 0.5)));
	}
	clearfits(&output_data);
	return MPP_OK;
}

mpp_status_t mpp::stack_apply_stsci(const std::vector<cv::Mat> &frames_raw,
                                    const std::vector<int> &included,
                                    const std::vector<double> &frame_brightness,
                                    const mpp_run_t *run,
                                    const mpp_config_t *cfg,
                                    fits *out) {
	/* In-memory frames are trivially safe to read concurrently. */
	return mpp::stack_apply_stsci_streamed(
	    [&frames_raw](int i) -> cv::Mat { return frames_raw[i]; },
	    (int) frames_raw.size(), included, frame_brightness, run, cfg, out,
	    /*max_threads=*/std::max(1, com.max_thread), /*provider_thread_safe=*/true);
}

/* ============================================================
 * Bayer drizzle (slice 5b.3).
 *
 * Same as the STScI path but:
 *   - input frames are single-channel raw Bayer (no debayer);
 *   - driz_args.is_bayer = TRUE with cfa[] + cfadim populated;
 *   - output is always 3-channel planar (R, G, B), even if input is
 *     single-channel.
 *
 * dobox's do_kernel_square (cdrizzlebox.c:998) reads the CFA channel
 * from the *input* pixel position, so the per-pixel routing is
 * intrinsic to (i, j) and works fine with our non-affine pixmap
 * (every input pixel still has a fixed colour identity regardless of
 * where the pixmap sends it).
 * ============================================================ */

mpp_status_t mpp::stack_apply_bayer_streamed(const FrameProvider &provider,
                                             int N,
                                             const std::vector<int> &included,
                                             const std::vector<double> &frame_brightness,
                                             const mpp_run_t *run,
                                             const mpp_config_t *cfg,
                                             const unsigned char *cfa,
                                             int cfadim,
                                             fits *out) {
	if (!cfg || !run || !out || !run->aps || !run->shifts || !cfa) return MPP_EINVAL;
	if (cfadim != 2) return MPP_EINVAL;     /* 2x2 Bayer only for now */
	if (N <= 0 || (int) included.size() != N
	 || (int) frame_brightness.size() != N) return MPP_EINVAL;
	if (N != run->num_frames) return MPP_EINVAL;

	const int frame_rx = run->frame_cols;
	const int frame_ry = run->frame_rows;
	if (frame_rx <= 0 || frame_ry <= 0) return MPP_EINVAL;

	const int intersection_rx = run->intersection[3] - run->intersection[2];
	const int intersection_ry = run->intersection[1] - run->intersection[0];
	if (intersection_rx <= 0 || intersection_ry <= 0) return MPP_EINVAL;

	const double scale = cfg->drizzle_scale < 1.0 ? 1.0 : cfg->drizzle_scale;
	const int out_rx   = (int) std::lround((double) intersection_rx * scale);
	const int out_ry   = (int) std::lround((double) intersection_ry * scale);
	const int channels = 3;                  /* Bayer output is always RGB */

	fits output_data{};
	fits output_counts{};
	if (!alloc_float_fits(&output_data, out_rx, out_ry, channels)) return MPP_ENOMEM;
	if (!alloc_float_fits(&output_counts, out_rx, out_ry, channels)) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}

	imgmap_t pixmap{};
	if (mpp_imgmap_alloc(&pixmap, frame_rx, frame_ry) != MPP_OK) {
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENOMEM;
	}

	std::vector<double> bright_inc;
	bright_inc.reserve(N);
	for (int i = 0; i < N; ++i)
		if (included[i]) bright_inc.push_back(frame_brightness[i]);
	if (bright_inc.empty()) {
		mpp_imgmap_free(&pixmap);
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENODATA;
	}
	std::nth_element(bright_inc.begin(),
	                 bright_inc.begin() + bright_inc.size() / 2,
	                 bright_inc.end());
	const double median_brightness = bright_inc[bright_inc.size() / 2];

	FrameFilterSetup filter = make_frame_filter(run, N, included);
	std::vector<DrizRect> driz_rects;        /* re-used per frame */
	bool driz_needs_claim = false;

	struct driz_args_t driz_args;
	std::memset(&driz_args, 0, sizeof(driz_args));
	driz_args.kernel         = kernel_from_cfg(cfg->drizzle_kernel);
	driz_args.scale          = (float) scale;
	driz_args.pixel_fraction = (float) cfg->drizzle_pixfrac;
	driz_args.weight_scale   = 1.0f;
	driz_args.is_bayer       = TRUE;
	std::memcpy(driz_args.cfa, cfa, 4);
	driz_args.cfadim         = (size_t) cfadim;
	driz_args.use_flats      = FALSE;
	driz_args.flat           = nullptr;

	struct driz_error_t error;
	std::memset(&error, 0, sizeof(error));

	const int n_included = (int) bright_inc.size();
	int cur_nb = 0;
	gui_iface.set_progress(PROGRESS_RESET, _("Stack (Bayer drizzle): processing frames"));

	for (int f = 0; f < N; ++f) {
		if (!included[f]) continue;
		if (!processing_should_continue()) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EINTR;
		}

		const int *ap_active_ptr = nullptr;
		int include_bg_flag = 1;
		if (!prepare_frame_filter(filter, f, &ap_active_ptr, &include_bg_flag)) {
			g_atomic_int_inc(&cur_nb);
			gui_iface.set_progress((double) cur_nb / (double) n_included, NULL);
			continue;
		}

		const cv::Mat frame_raw = provider(f);
		if (frame_raw.empty()) continue;
		if (frame_raw.channels() != 1) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EINVAL;
		}

		const double brightness_scale = median_brightness
		                              / (frame_brightness[f] + 1e-7);

		fits frame_fit{};
		if (!cv_to_planar_fits(frame_raw, brightness_scale, &frame_fit)) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_ENOMEM;
		}

		if (mpp_pixmap_build_filtered(run, f, scale, ap_active_ptr,
		                              include_bg_flag, &pixmap) != MPP_OK) {
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EINVAL;
		}

		compute_driz_plan(run, ap_active_ptr, include_bg_flag,
		                  frame_rx, frame_ry,
		                  driz_rects, driz_needs_claim);

		bool dobox_failed = false;
		for (size_t r = 0; r < driz_rects.size() && !dobox_failed; ++r) {
			const DrizRect &rect = driz_rects[r];
			struct driz_param_t p;
			std::memset(&p, 0, sizeof(p));
			driz_param_init(&p);
			p.driz           = &driz_args;
			p.kernel         = driz_args.kernel;
			/* Bayer-drizzle drop-size correction: dobox's turbo and
			 * gaussian kernels compute drop half-width as pfo =
			 * pixfrac · p.scale / 2 (cdrizzlebox.c around lines 829
			 * and 504). The pre-existing code set p.scale = 1.0,
			 * leaving the drop always 1 output pixel wide regardless
			 * of drizzle_scale. At scale=2 the inter-input-pixel gap
			 * in output coords is 2, and same-CFA-colour inputs (R,
			 * B) are 2 input pixels apart = 4 output pixels — 1-wide
			 * drops leave 3-pixel gaps where that colour is never
			 * deposited from any frame without large sub-pixel
			 * diversity, producing the 2-pixel-period vertical-stripe
			 * artefact wavelets bring out. (Square's kernel sets
			 * drop bounds via the pixmap transform of input corners,
			 * so it is correct regardless; lanczos uses
			 * input_pixel_size_in_output too — both kernels' tests
			 * pass either way. Turbo is the GUI default.)
			 *
			 * Setting p.scale = drizzle_scale fixes pfo uniformly
			 * across all affected kernels. As a side effect it
			 * scales the per-input-pixel value by drizzle_scale²
			 * (the flux-conservation factor `d *= scale2` inside
			 * each kernel), inflating output brightness by the same
			 * factor — undone post-stack by dividing the
			 * accumulator by drizzle_scale². Brightness is then
			 * indistinguishable from the previous p.scale=1 path. */
			p.scale          = (float) scale;
			p.pixel_fraction = driz_args.pixel_fraction;
			p.weight_scale   = 1.0f;
			p.in_units       = unit_counts;
			p.out_units      = unit_counts;
			p.exposure_time  = 1.0f;
			p.data           = &frame_fit;
			p.weights        = nullptr;
			p.pixmap         = &pixmap;
			p.output_data    = &output_data;
			p.output_counts  = &output_counts;
			std::memcpy(p.cfa, cfa, 4);
			p.cfadim         = (size_t) cfadim;
			p.xmin           = rect.xmin;
			p.ymin           = rect.ymin;
			p.xmax           = rect.xmax;
			p.ymax           = rect.ymax;
			p.threads        = 1;
			p.error          = &error;

			if (dobox(&p)) {
				dobox_failed = true;
				break;
			}
			if (driz_needs_claim) {
				pixmap_claim_rect(&pixmap, rect);
			}
		}
		if (dobox_failed) {
			siril_log_error(
			    _("Stack (Bayer Drizzle): dobox failed on frame %d: %s\n"),
			    f, error.last_message[0] ? error.last_message : _("(no detail)"));
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EIO;
		}
		clearfits(&frame_fit);

		++cur_nb;
		gui_iface.set_progress((double) cur_nb / (double) n_included, NULL);
	}

	mpp_imgmap_free(&pixmap);
	clearfits(&output_counts);

	clearfits(out);
	out->rx = out_rx;
	out->ry = out_ry;
	out->naxes[0] = out_rx;
	out->naxes[1] = out_ry;
	out->naxes[2] = channels;
	out->naxis = 3;
	out->bitpix = USHORT_IMG;
	out->orig_bitpix = USHORT_IMG;
	out->type = DATA_USHORT;
	const size_t plane = (size_t) out_rx * (size_t) out_ry;
	out->data = (WORD *) std::calloc(plane * (size_t) channels, sizeof(WORD));
	if (!out->data) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}
	out->pdata[0] = out->data;
	out->pdata[1] = out->data + plane;
	out->pdata[2] = out->data + 2 * plane;

	/* Undo the drizzle_scale² flux-conservation multiplier that
	 * p.scale = drizzle_scale introduced inside dobox (see comment
	 * around p.scale in the per-frame loop). With p.scale=scale,
	 * each input pixel's value got multiplied by scale²; the
	 * weighted-mean accumulator preserves that, so output is
	 * scale²× too bright. Dividing here restores input-comparable
	 * brightness without affecting the spatial drop coverage that
	 * scale² in d was paid to obtain. */
	const double cast_scale = ((cfg->bitdepth == 8) ? 256.0 : 1.0) / (scale * scale);
	for (size_t k = 0; k < plane * (size_t) channels; ++k) {
		const double v = (double) output_data.fdata[k] * cast_scale;
		out->data[k] = (WORD) (v < 0.0 ? 0
		                     : (v > 65535.0 ? 65535
		                                    : (int) (v + 0.5)));
	}
	clearfits(&output_data);
	return MPP_OK;
}

mpp_status_t mpp::stack_apply_bayer(const std::vector<cv::Mat> &frames_raw_bayer,
                                    const std::vector<int> &included,
                                    const std::vector<double> &frame_brightness,
                                    const mpp_run_t *run,
                                    const mpp_config_t *cfg,
                                    const unsigned char *cfa,
                                    int cfadim,
                                    fits *out) {
	for (const auto &f : frames_raw_bayer) {
		if (!f.empty() && f.channels() != 1) return MPP_EINVAL;
	}
	return mpp::stack_apply_bayer_streamed(
	    [&frames_raw_bayer](int i) -> cv::Mat { return frames_raw_bayer[i]; },
	    (int) frames_raw_bayer.size(), included, frame_brightness, run, cfg,
	    cfa, cfadim, out);
}

/* Sequence-side wrapper: streams each included frame via Siril's
 * sequence I/O and forwards to the inner streamed C++ function. The
 * dobox loop pays the per-frame disk read inline rather than pre-
 * loading the whole sequence — memory peak collapses from N frames to
 * one, disk traffic is identical (one read per included frame either
 * way). */
extern "C" mpp_status_t mpp_stack_apply_stsci(sequence *seq,
                                              const mpp_config_t *cfg,
                                              const mpp_run_t *run,
                                              fits *out) {
	if (!seq || !cfg || !run || !run->included || !run->frame_brightness)
		return MPP_EINVAL;
	const int N = run->num_frames;
	std::vector<int> included(run->included, run->included + N);
	std::vector<double> brightness(run->frame_brightness,
	                                run->frame_brightness + N);
	auto provider = [seq, cfg](int i) -> cv::Mat {
		return mpp::read_full_frame(seq, i, cfg->avi_bayer_pattern);
	};
	/* Parallel decode needs a reentrant provider — same predicate as Stage
	 * A's parallel frame read (SER serialises only the raw fread; FITS only
	 * when cfitsio was built reentrant). */
	const bool provider_safe =
	    (seq->type == SEQ_SER
	     || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ
	          || seq->type == SEQ_INTERNAL) && fits_is_reentrant()));
	return mpp::stack_apply_stsci_streamed(provider, N, included, brightness,
	                                       run, cfg, out,
	                                       std::max(1, com.max_thread), provider_safe);
}

/* Sequence-side wrapper for Bayer drizzle. Reads each frame raw (no
 * debayer, bypassing com.pref.debayer.open_debayer) and forwards to
 * the inner mpp::stack_apply_bayer.
 *
 * Only SER sequences with a Bayer ColorID are supported. The CFA
 * pattern is derived from the SER header's color_id, matching the
 * 2x2 RLAYER/GLAYER/BLAYER convention Siril's get_compiled_pattern
 * produces (see demosaicing.c:333). */
namespace {
/* Returns FALSE if the sequence has no usable Bayer pattern. cfa is
 * filled with 4 bytes on success.
 *
 * ser_read_frame Y-flips the image buffer to convert SER's top-down
 * raster order to FITS' bottom-up convention (ser.c::ser_read_frame
 * line 970). The CFA pattern from the SER header describes the on-disk
 * orientation; after the flip, what arrives at our cv::Mat is rotated.
 * We use Siril's adjust_Bayer_pattern_orientation helper so the cfa[]
 * we feed to dobox matches the in-memory mosaic, not the file's
 * declared pattern. */
bool bayer_cfa_from_ser(struct sequ *seq, unsigned int ry, unsigned char cfa[4]) {
	if (!seq || seq->type != SEQ_SER || !seq->ser_file) return false;
	sensor_pattern pat;
	switch (seq->ser_file->color_id) {
		case SER_BAYER_RGGB: pat = BAYER_FILTER_RGGB; break;
		case SER_BAYER_GRBG: pat = BAYER_FILTER_GRBG; break;
		case SER_BAYER_GBRG: pat = BAYER_FILTER_GBRG; break;
		case SER_BAYER_BGGR: pat = BAYER_FILTER_BGGR; break;
		default: return false;
	}
	adjust_Bayer_pattern_orientation(&pat, ry, TRUE);
	const int R = 0, G = 1, B = 2;
	switch (pat) {
		case BAYER_FILTER_RGGB: cfa[0]=R; cfa[1]=G; cfa[2]=G; cfa[3]=B; return true;
		case BAYER_FILTER_GRBG: cfa[0]=G; cfa[1]=R; cfa[2]=B; cfa[3]=G; return true;
		case BAYER_FILTER_GBRG: cfa[0]=G; cfa[1]=B; cfa[2]=R; cfa[3]=G; return true;
		case BAYER_FILTER_BGGR: cfa[0]=B; cfa[1]=G; cfa[2]=G; cfa[3]=R; return true;
		default: return false;
	}
}
}  // namespace

extern "C" mpp_status_t mpp_stack_apply_bayer(sequence *seq,
                                              const mpp_config_t *cfg,
                                              const mpp_run_t *run,
                                              fits *out) {
	if (!seq || !cfg || !run || !run->included || !run->frame_brightness)
		return MPP_EINVAL;

	unsigned char cfa[4];
	if (seq->type == SEQ_SER && seq->ser_file) {
		if (!bayer_cfa_from_ser(seq, (unsigned int) seq->ry, cfa)) {
			siril_log_error(_("Bayer drizzle: SER ColorID is not RGGB / GRBG / "
			                          "GBRG / BGGR — cannot derive CFA pattern.\n"));
			return MPP_EINVAL;
		}
	} else if (seq->type == SEQ_AVI) {
		/* AVI Bayer: pattern comes from the MPP-tab combo (no container-
		 * level marker exists). film_read_frame's existing Y-flip means
		 * the user's top-down pick must be adjusted to bottom-up here. */
		if (!mpp::avi_bayer_to_cfa(cfg->avi_bayer_pattern,
		                            (int) seq->ry, cfa)) {
			siril_log_error(_("Bayer drizzle: AVI sequence has no Bayer pattern set "
			                          "(MPP tab → \"AVI Bayer pattern\").\n"));
			return MPP_EINVAL;
		}
	} else {
		siril_log_error(_("Bayer drizzle: requires a SER or AVI sequence.\n"));
		return MPP_EINVAL;
	}

	const int N = run->num_frames;
	std::vector<int> included(run->included, run->included + N);
	std::vector<double> brightness(run->frame_brightness,
	                                run->frame_brightness + N);

	/* For SER, ser_read_frame's `open_debayer` parameter is a latent
	 * NOP for SERs already opened with debayer-on-open: the actual
	 * debayer call inside ser_read_frame dispatches on
	 * ser_file->debayer_type_ser (set at SER open time from
	 * com.pref.debayer.open_debayer), not on the parameter. To
	 * actually get raw mosaic bytes back we have to flip that field
	 * to SER_MONO for the duration of the provider's reads and
	 * restore on the way out. Stage A/B already ran upstream on
	 * debayered frames; nothing else is reading from this SER
	 * concurrently because Stage C is sequential by construction.
	 * AVI doesn't need this dance — film_read_frame already returns
	 * a 1-layer mosaic, so we just call read_full_frame with no
	 * pattern hint (avi_pattern=AUTO suppresses the cv::cvtColor
	 * debayer step that the bicubic path uses).
	 *
	 * RAII so an exception escaping stack_apply_bayer_streamed (cv::
	 * Exception, std::bad_alloc, …) still restores the SER state. */
	struct SerRawGuard {
		sequence *seq = nullptr;
		ser_color saved = SER_MONO;
		bool active = false;
		~SerRawGuard() {
			if (active && seq && seq->ser_file)
				seq->ser_file->debayer_type_ser = saved;
		}
	} ser_raw_guard;
	const bool is_ser = (seq->type == SEQ_SER);
	if (is_ser) {
		ser_raw_guard.seq = seq;
		ser_raw_guard.saved = seq->ser_file->debayer_type_ser;
		ser_raw_guard.active = true;
		seq->ser_file->debayer_type_ser = SER_MONO;
	}

	auto provider = [seq, is_ser](int i) -> cv::Mat {
		if (is_ser) {
			fits f{};
			if (ser_read_frame(seq->ser_file, i, &f, /*force_float=*/FALSE,
			                   /*open_debayer=*/FALSE)) {
				clearfits(&f);
				return cv::Mat();
			}
			cv::Mat m(f.ry, f.rx, CV_16U, f.data);
			cv::Mat owned = m.clone();   /* own the data before clearfits frees f */
			clearfits(&f);
			return owned;
		}
		/* AVI: ask read_full_frame for the raw 1-layer mosaic (avi_pattern
		 * AUTO → no cv::cvtColor debayer). */
		return mpp::read_full_frame(seq, i, MPP_AVI_BAYER_AUTO);
	};
	return mpp::stack_apply_bayer_streamed(
	    provider, N, included, brightness, run, cfg, cfa, /*cfadim=*/2, out);
}
