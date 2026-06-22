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
 */

#include <algorithm>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp/mpp_multistack.hpp"
#include "registration/mpp/mpp_multisource.hpp"
#include "registration/mpp/mpp_derot.h"          /* mpp_derot_frame_map_ms */
#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_ap_priv.hpp"
#include "registration/mpp/mpp_stack_priv.hpp"

extern "C" {
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_shift.h"
}

namespace mpp {

MultiStackResult multistack_channel(const std::vector<MsSource> &srcs,
                                    const mpp_derot_t *out_derot, int num_layers,
                                    const mpp_config_t &cfg,
                                    int max_threads) {
	MultiStackResult res;
	if (srcs.empty() || !out_derot) {
		res.error = true;
		return res;
	}
	const mpp_derot_t *refd = out_derot;

	/* output (reference epoch) canvas */
	const int out_w = refd->frame_cols;
	const int out_h = refd->frame_rows;
	if (out_w <= 0 || out_h <= 0) { res.error = true; return res; }

	/* global index space over every source's frames */
	MultiSourceLayout layout;
	for (const auto &s : srcs)
		layout.add(s.num_frames);
	const int N = layout.total();
	if (N <= 0) { res.error = true; return res; }
	res.num_frames = N;

	const int nt = std::max(1, max_threads);

	/* Relocate+derotate a source frame into the reference epoch canvas: the C3
	 * ms map sends each reference-canvas pixel back to this source's frame
	 * (out side = reference disk, source side = this sequence's own disk). */
	auto relocate = [&](const cv::Mat &m, int s, int local) -> cv::Mat {
		if (m.empty()) return cv::Mat();
		cv::Mat mx, my, mu, out;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, local, out_w, out_h,
		                       0.0, 0.0, 1.0, 0.0, 0.0, 1.0, mx, my, mu,
		                       /*mask_outside=*/true);
		cv::remap(m, out, mx, my, cv::INTER_LANCZOS4, cv::BORDER_CONSTANT,
		          cv::Scalar(0));
		return out;
	};

	/* Analysis providers — raw mono and blurred mono, both in the reference
	 * canvas. Used by the global aligner, the average reference, the AP
	 * qualities and the per-AP shift correlation. */
	auto raw_provider = [&](int g) -> cv::Mat {
		int s, l;
		if (!layout.locate(g, &s, &l)) return cv::Mat();
		return relocate(srcs[s].analysis_read(l), s, l);
	};
	auto blurred_provider = [&](int g) -> cv::Mat {
		cv::Mat m = raw_provider(g);
		if (m.empty()) return cv::Mat();
		return mpp::blur_mono_for_align(m, cfg);
	};

	/* Pass 1 equivalent: per-frame quality + brightness on the relocated raw
	 * frames (no caching in this first multi-source implementation — every
	 * stage streams from the source readers). */
	std::vector<double> quality(N, 0.0), brightness(N, 0.0);
	std::vector<int> included(N, 1);
	for (int g = 0; g < N; ++g) {
		const cv::Mat m = raw_provider(g);
		if (m.empty()) { res.error = true; return res; }
		quality[g]    = mpp::rank_score_normalized(m, cfg);
		brightness[g] = mpp::rank_average_brightness(m, cfg);
	}

	/* Global alignment of the relocated frames (small residuals — derotation
	 * already co-registers the disks). */
	const auto align = mpp::align_global_from_provider(
	    blurred_provider, N, quality, cfg, nullptr, nullptr,
	    /*provider_is_thread_safe=*/false);
	if (align.oom) { res.oom = true; return res; }
	if (align.best_frame_idx < 0) { res.cancelled = true; return res; }

	/* Averaged reference + frame intersection, in the reference canvas. */
	const auto avg = mpp::align_average_frame_streamed(
	    raw_provider, N, out_h, out_w, quality, align.shifts, cfg,
	    nullptr, nullptr, included);
	if (avg.mean_frame.empty()) { res.cancelled = true; return res; }

	/* AP grid on the reference. */
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	if (!aps || aps->count <= 0) { mpp_ap_free(aps); res.error = true; return res; }
	res.num_aps = aps->count;

	std::vector<FrameOffset> offsets(N);
	for (int g = 0; g < N; ++g) {
		offsets[g].dy = align.shifts[g][0];
		offsets[g].dx = align.shifts[g][1];
	}

	/* Per-AP per-frame qualities over the union. */
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, brightness, *aps, offsets, out_h, out_w, cfg,
	    nullptr, nullptr, included, nt, /*provider_thread_safe=*/false);

	/* Per-AP per-frame shifts (epoch space, against the shared reference). */
	mpp_shifts_t *shifts = mpp::stack_compute_shifts_streamed(
	    blurred_provider, N, avg.mean_frame, *aps, apq, offsets, cfg,
	    included.data(), nt, /*provider_thread_safe=*/false);
	if (!shifts) { mpp_ap_free(aps); res.oom = true; return res; }
	if (shifts->failure_counter == -1) {
		mpp_shift_free(shifts);
		mpp_ap_free(aps);
		res.cancelled = true;
		return res;
	}

	/* Quality-sorted index (included frames only) for the background pick. */
	std::vector<int> sorted_idx;
	sorted_idx.reserve(N);
	for (int g = 0; g < N; ++g)
		if (included[g]) sorted_idx.push_back(g);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [&](int a, int b) { return quality[a] > quality[b]; });

	const cv::Vec4i intersection(avg.intersection[0], avg.intersection[1],
	                             avg.intersection[2], avg.intersection[3]);
	res.intersection = intersection;

	/* Full-frame stack provider: returns each source's NATIVE full frame; the
	 * derot map (below) does the relocation+derotation in the single resample. */
	auto full_provider = [&](int g) -> cv::Mat {
		int s, l;
		if (!layout.locate(g, &s, &l)) return cv::Mat();
		return srcs[s].full_read(l);
	};

	/* Warp-stack derot map: output = reference disk in the drizzled
	 * intersection canvas, source = this frame's own sequence disk. */
	const double S = std::max(1.0, cfg.drizzle_scale);
	DerotMapProvider derot_provider =
	    [&](int g, int DY, int DX, cv::Mat &mx, cv::Mat &my, cv::Mat &mu) -> bool {
		int s, l;
		if (!layout.locate(g, &s, &l)) return false;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, l, DX, DY,
		                       (double) intersection[2], (double) intersection[0], S,
		                       0.0, 0.0, 1.0, mx, my, mu, /*mask_outside=*/true);
		return true;
	};

	const auto wr = mpp::stack_warp_apply_streamed(
	    full_provider, N, num_layers, *aps, apq, shifts, offsets, brightness,
	    sorted_idx, intersection, cfg, included.data(), nt,
	    /*provider_thread_safe=*/false, /*mem_budget_bytes=*/0, &derot_provider);

	mpp_shift_free(shifts);
	mpp_ap_free(aps);

	if (wr.oom) { res.oom = true; return res; }
	if (wr.cancelled || wr.image.empty()) { res.cancelled = true; return res; }
	res.image = wr.image;
	return res;
}

}  // namespace mpp
