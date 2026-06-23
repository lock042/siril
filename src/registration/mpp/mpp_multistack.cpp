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
#include <atomic>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
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

/* The stack-family engines self-report into the GUI progress bar; the
 * align/average/AP-quality engines are callback-driven, so feed them the same
 * sink (a stage fills 0..1, the caller sets the descriptive text). Safe to call
 * headless — gui_iface.set_progress is a no-op stub there. */
namespace {
void ms_progress(double f, void *) { gui_iface.set_progress(f, NULL); }
}  // namespace

MultiStackResult multistack_channel(const std::vector<MsSource> &srcs,
                                    const mpp_derot_t *out_derot, int num_layers,
                                    const mpp_config_t &cfg,
                                    int max_threads, bool provider_thread_safe) {
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

	/* Mask everything off the globe only for sources whose disk differs from the
	 * reference (a displaced/other sequence) — there, identity passthrough would
	 * ghost their globe into the reference sky. The reference's OWN frames keep
	 * the single-sequence behaviour: rings and sky pass through and stack
	 * normally. (A single-sequence combine is then identical to the ordinary
	 * derotation stack.) */
	auto src_masks = [&](int s) -> bool { return srcs[s].derot != refd; };

	/* Relocate+derotate a source frame into the reference epoch canvas: the C3
	 * ms map sends each reference-canvas pixel back to this source's frame
	 * (out side = reference disk, source side = this sequence's own disk). */
	auto relocate = [&](const cv::Mat &m, int s, int local) -> cv::Mat {
		if (m.empty()) return cv::Mat();
		cv::Mat mx, my, mu, out;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, local, out_w, out_h,
		                       0.0, 0.0, 1.0, 0.0, 0.0, 1.0, mx, my, mu,
		                       src_masks(s));
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
	 * stage streams from the source readers). Parallelised over frames when the
	 * readers are reentrant (SER / reentrant FITS), like the single-sequence
	 * Stage A: this is the first of several streamed passes, so serialising it
	 * — and the engines below — is the main multi-source throughput limiter. */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: reading and ranking frames"));
	std::vector<double> quality(N, 0.0), brightness(N, 0.0);
	std::vector<int> included(N, 1);
	std::atomic<bool> p1_fail{false};
	std::atomic<int> p1_done{0};
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(guided) \
    if (provider_thread_safe && nt > 1)
#endif
	for (int g = 0; g < N; ++g) {
		if (p1_fail.load()) continue;
		const cv::Mat m = raw_provider(g);
		if (m.empty()) { p1_fail.store(true); continue; }
		quality[g]    = mpp::rank_score_normalized(m, cfg);
		brightness[g] = mpp::rank_average_brightness(m, cfg);
		const int d = ++p1_done;
		gui_iface.set_progress((double) d / (double) N, NULL);
	}
	if (p1_fail.load()) {
		siril_log_error(_("Derotation stack: could not read/relocate every frame "
		                  "for analysis.\n"));
		res.error = true; return res;
	}

	/* Global alignment of the relocated frames (small residuals — derotation
	 * already co-registers the disks). */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: aligning frames"));
	const auto align = mpp::align_global_from_provider(
	    blurred_provider, N, quality, cfg, ms_progress, nullptr,
	    provider_thread_safe);
	if (align.oom) { res.oom = true; return res; }
	if (align.best_frame_idx < 0) { res.cancelled = true; return res; }

	/* Averaged reference + frame intersection, in the reference canvas. */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: building reference"));
	const auto avg = mpp::align_average_frame_streamed(
	    raw_provider, N, out_h, out_w, quality, align.shifts, cfg,
	    ms_progress, nullptr, included);
	if (avg.mean_frame.empty()) { res.cancelled = true; return res; }

	/* AP grid on the reference. */
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	if (!aps || aps->count <= 0) {
		mpp_ap_free(aps);
		siril_log_error(_("Derotation stack: no alignment points could be placed "
		                  "on the reference — check the disk fit.\n"));
		res.error = true; return res;
	}
	res.num_aps = aps->count;

	std::vector<FrameOffset> offsets(N);
	for (int g = 0; g < N; ++g) {
		offsets[g].dy = align.shifts[g][0];
		offsets[g].dx = align.shifts[g][1];
	}

	/* AP grid built — log a count so a degenerate fit (e.g. an over-masked or
	 * empty reference) is diagnosable rather than a bare error code. */
	siril_log_message(_("Derotation stack: %d alignment point(s) on the reference\n"),
	                  aps->count);

	/* Per-AP per-frame qualities over the union. */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: ranking alignment points"));
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, brightness, *aps, offsets, out_h, out_w, cfg,
	    ms_progress, nullptr, included, nt, provider_thread_safe);

	/* Per-AP per-frame shifts (epoch space, against the shared reference). */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: measuring local shifts"));
	mpp_shifts_t *shifts = mpp::stack_compute_shifts_streamed(
	    blurred_provider, N, avg.mean_frame, *aps, apq, offsets, cfg,
	    included.data(), nt, provider_thread_safe);
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
		                       0.0, 0.0, 1.0, mx, my, mu, src_masks(s));
		return true;
	};

	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: stacking frames"));
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
