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
#include "core/processing.h"                   /* processing_should_continue */
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

/* Shared per-call context: layout + relocating providers over the sources.
 * Rebuilt cheaply in each phase (the lambdas close over `srcs`/`refd`, which
 * must outlive the phase call — they do: both phases run synchronously). */
namespace {
struct MsContext {
	const mpp_derot_t *refd = nullptr;
	int out_w = 0, out_h = 0;
	MultiSourceLayout layout;
	int N = 0;
	std::vector<int> included;
	int n_excluded = 0;
	bool ok = false;
};

MsContext ms_context(const std::vector<MsSource> &srcs,
                     const mpp_derot_t *out_derot) {
	MsContext c;
	if (srcs.empty() || !out_derot) return c;
	c.refd = out_derot;
	c.out_w = out_derot->frame_cols;
	c.out_h = out_derot->frame_rows;
	if (c.out_w <= 0 || c.out_h <= 0) return c;
	for (const auto &s : srcs)
		c.layout.add(s.num_frames);
	c.N = c.layout.total();
	if (c.N <= 0) return c;

	/* Frame-selector state, folded into the global index space. Excluded
	 * frames are never read: they stay out of the ranking, the reference,
	 * the AP selection and the stack — the single-sequence contract. */
	c.included.assign(c.N, 1);
	for (int s = 0; s < (int) srcs.size(); ++s) {
		if (srcs[s].included.empty()) continue;
		for (int l = 0; l < srcs[s].num_frames
		             && l < (int) srcs[s].included.size(); ++l)
			if (!srcs[s].included[l]) {
				c.included[c.layout.base_of(s) + l] = 0;
				++c.n_excluded;
			}
	}
	c.ok = true;
	return c;
}
}  // namespace

MultiStackAnalysis multistack_analyze(const std::vector<MsSource> &srcs,
                                      const mpp_derot_t *out_derot,
                                      const mpp_config_t &cfg,
                                      int max_threads, bool provider_thread_safe) {
	MultiStackAnalysis res;
	MsContext ctx = ms_context(srcs, out_derot);
	if (!ctx.ok) { res.error = true; return res; }
	const mpp_derot_t *refd = ctx.refd;
	const int out_w = ctx.out_w, out_h = ctx.out_h;
	const int N = ctx.N;
	const auto &layout = ctx.layout;
	res.num_frames = N;

	const int nt = std::max(1, max_threads);

	/* Relocate+derotate a source frame into the reference epoch canvas: the C3
	 * ms map sends each reference-canvas pixel back to this source's frame
	 * (out side = reference disk, source side = this sequence's own disk).
	 *
	 * The analysis relocation must NOT hard-mask off-globe content (the stack
	 * maps do, see multistack_stack): the mask is anchored to the FITTED disk
	 * while the content drifts frame to frame (and the fit itself carries an
	 * error), so masking would cookie-cut the real planet at a static ellipse
	 * — the global aligner then locks onto that artificial edge instead of
	 * the content, the mean reference grows a dark crescent, the AP
	 * thresholds refuse APs there, and the final stack inherits the bite
	 * (the badmulti regression). */
	auto relocate = [&](const cv::Mat &m, int s, int local) -> cv::Mat {
		if (m.empty()) return cv::Mat();
		cv::Mat mx, my, mu, out;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, local, out_w, out_h,
		                       0.0, 0.0, 1.0, 0.0, 0.0, 1.0, mx, my, mu,
		                       /*mask_outside=*/false);
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

	const std::vector<int> &included = ctx.included;
	if (ctx.n_excluded > 0)
		siril_log_message(_("Derotation stack: %d frame(s) excluded by the "
		                    "frame selector\n"), ctx.n_excluded);
	if (ctx.n_excluded >= N) {
		siril_log_error(_("Derotation stack: every frame is excluded.\n"));
		res.error = true; return res;
	}

	std::atomic<bool> p1_fail{false};
	std::atomic<bool> p1_cancel{false};
	std::atomic<int> p1_done{0};
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(guided) \
    if (provider_thread_safe && nt > 1)
#endif
	for (int g = 0; g < N; ++g) {
		if (p1_fail.load() || p1_cancel.load() || !included[g]) continue;
		if (!processing_should_continue()) { p1_cancel.store(true); continue; }
		const cv::Mat m = raw_provider(g);
		if (m.empty()) { p1_fail.store(true); continue; }
		quality[g]    = mpp::rank_score_normalized(m, cfg);
		brightness[g] = mpp::rank_average_brightness(m, cfg);
		const int d = ++p1_done;
		gui_iface.set_progress((double) d / (double) N, NULL);
	}
	if (p1_cancel.load()) { res.cancelled = true; return res; }
	if (p1_fail.load()) {
		siril_log_error(_("Derotation stack: could not read/relocate every frame "
		                  "for analysis.\n"));
		res.error = true; return res;
	}

	/* Masked quality for selection: excluded frames sit below the lowest
	 * real score (>= 0) so they can never be picked as the global reference
	 * or averaged into the mean frame. */
	std::vector<double> q_sel = quality;
	if (ctx.n_excluded > 0)
		for (int g = 0; g < N; ++g)
			if (!included[g]) q_sel[g] = -1.0;

	/* Global alignment of the relocated frames (small residuals — derotation
	 * already co-registers the disks). */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: aligning frames"));
	const auto align = mpp::align_global_from_provider(
	    blurred_provider, N, q_sel, cfg, ms_progress, nullptr,
	    provider_thread_safe);
	if (align.oom) { res.oom = true; return res; }
	if (align.best_frame_idx < 0) { res.cancelled = true; return res; }

	/* Averaged reference + frame intersection, in the reference canvas.
	 * The visibility provider (the derot map's mu, in the same reference
	 * canvas as the relocated frame) makes the mean per-pixel normalised
	 * over the frames that actually see each pixel: a surface region hidden
	 * behind the limb in one video keeps its full brightness from the other
	 * — so the AP thresholds still place APs there and the region stacks
	 * from the video that covers it. */
	auto vis_provider = [&](int g) -> cv::Mat {
		int s, l;
		if (!layout.locate(g, &s, &l)) return cv::Mat();
		cv::Mat mx, my, mu;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, l, out_w, out_h,
		                       0.0, 0.0, 1.0, 0.0, 0.0, 1.0, mx, my, mu,
		                       /*mask_outside=*/false);
		return mu;
	};
	const mpp::FrameProvider vis_fn = vis_provider;
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: building reference"));
	const auto avg = mpp::align_average_frame_streamed(
	    raw_provider, N, out_h, out_w, q_sel, align.shifts, cfg,
	    ms_progress, nullptr, included, &vis_fn);
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
	/* AP grid built — log a count so a degenerate fit (e.g. an over-masked or
	 * empty reference) is diagnosable rather than a bare error code. */
	siril_log_message(_("Derotation stack: %d alignment point(s) on the reference\n"),
	                  aps->count);

	res.best_frame_idx = align.best_frame_idx;
	res.quality = std::move(quality);
	res.brightness = std::move(brightness);
	res.shifts = align.shifts;
	res.intersection = cv::Vec4i(avg.intersection[0], avg.intersection[1],
	                             avg.intersection[2], avg.intersection[3]);
	res.mean_frame = avg.mean_frame;
	res.aps = aps;
	return res;
}

MultiStackResult multistack_stack(const std::vector<MsSource> &srcs,
                                  const mpp_derot_t *out_derot, int num_layers,
                                  const mpp_config_t &cfg,
                                  const MultiStackAnalysis &an,
                                  const mpp_aps_t &aps_in,
                                  int max_threads, bool provider_thread_safe) {
	MultiStackResult res;
	MsContext ctx = ms_context(srcs, out_derot);
	if (!ctx.ok || an.error || an.oom || an.cancelled
	    || an.num_frames != ctx.N || (int) an.shifts.size() != ctx.N
	    || an.mean_frame.empty() || aps_in.count <= 0) {
		res.error = true;
		return res;
	}
	const mpp_derot_t *refd = ctx.refd;
	const int out_w = ctx.out_w, out_h = ctx.out_h;
	const int N = ctx.N;
	const auto &layout = ctx.layout;
	res.num_frames = N;
	res.num_aps = aps_in.count;
	res.intersection = an.intersection;
	const mpp_aps_t *aps = &aps_in;

	const int nt = std::max(1, max_threads);

	/* Frame-selector state re-read live from the sources, mirroring the
	 * single-sequence Stage C inclusion refresh (the user may have excluded
	 * frames between analysis and stack). */
	const std::vector<int> &included = ctx.included;
	if (ctx.n_excluded >= N) {
		siril_log_error(_("Derotation stack: every frame is excluded.\n"));
		res.error = true;
		return res;
	}

	/* Analysis providers in the reference canvas — same relocation as the
	 * analysis phase (no off-globe masking; see multistack_analyze). */
	auto raw_provider = [&](int g) -> cv::Mat {
		int s, l;
		if (!layout.locate(g, &s, &l)) return cv::Mat();
		const cv::Mat m = srcs[s].analysis_read(l);
		if (m.empty()) return cv::Mat();
		cv::Mat mx, my, mu, out;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, l, out_w, out_h,
		                       0.0, 0.0, 1.0, 0.0, 0.0, 1.0, mx, my, mu,
		                       /*mask_outside=*/false);
		cv::remap(m, out, mx, my, cv::INTER_LANCZOS4, cv::BORDER_CONSTANT,
		          cv::Scalar(0));
		return out;
	};
	auto blurred_provider = [&](int g) -> cv::Mat {
		cv::Mat m = raw_provider(g);
		if (m.empty()) return cv::Mat();
		return mpp::blur_mono_for_align(m, cfg);
	};

	/* Mask everything off the globe only for sources that do not share the
	 * output canvas disk (MsSource.shares_output_disk, decided by the
	 * caller from plan CONTENT, not pointer identity) — there, identity
	 * passthrough would ghost their globe into the reference sky. Sources
	 * sharing the output disk keep the single-sequence behaviour: rings
	 * and sky pass through and stack normally. (A single-sequence combine
	 * is then identical to the ordinary derotation stack.) */
	auto src_masks = [&](int s) -> bool { return !srcs[s].shares_output_disk; };

	/* Frame offsets in the SAME convention the engines expect (matches
	 * offsets_from_run): the intersection origin minus the global shift, NOT the
	 * raw shift. The AP boxes and the warp map are positioned relative to the
	 * intersection canvas, so this term carries the crop origin and the
	 * alignment in one — passing the raw shift mis-places every frame by the
	 * intersection origin (and with the wrong sign), smearing the stack. */
	std::vector<FrameOffset> offsets(N);
	for (int g = 0; g < N; ++g) {
		offsets[g].dy = (double) an.intersection[0] - an.shifts[g][0];
		offsets[g].dx = (double) an.intersection[2] - an.shifts[g][1];
	}

	/* Per-AP per-frame qualities over the union. */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: ranking alignment points"));
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, an.brightness, *aps, offsets, out_h, out_w, cfg,
	    ms_progress, nullptr, included, nt, provider_thread_safe);
	if (apq.stack_size <= 0) {   /* cancellation / OOM sentinel */
		res.cancelled = true;
		return res;
	}

	/* Per-AP per-frame shifts (epoch space, against the shared reference). */
	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: measuring local shifts"));
	mpp_shifts_t *shifts = mpp::stack_compute_shifts_streamed(
	    blurred_provider, N, an.mean_frame, *aps, apq, offsets, cfg,
	    included.data(), nt, provider_thread_safe);
	if (!shifts) { res.oom = true; return res; }
	if (shifts->failure_counter == -1) {
		mpp_shift_free(shifts);
		res.cancelled = true;
		return res;
	}

	/* Quality-sorted index (included frames only) for the background pick. */
	std::vector<int> sorted_idx;
	sorted_idx.reserve(N);
	for (int g = 0; g < N; ++g)
		if (included[g]) sorted_idx.push_back(g);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [&](int a, int b) { return an.quality[a] > an.quality[b]; });

	const cv::Vec4i intersection = an.intersection;

	/* Full-frame stack provider: returns each source's NATIVE full frame; the
	 * derot map (below) does the relocation+derotation in the single resample. */
	auto full_provider = [&](int g) -> cv::Mat {
		int s, l;
		if (!layout.locate(g, &s, &l)) return cv::Mat();
		return srcs[s].full_read(l);
	};

	/* Warp-stack derot map: output = reference disk in the drizzled intersection
	 * canvas, source = this frame's own sequence disk. BOTH sides use the same
	 * intersection origin and drizzle scale: the warp engine's map is
	 *   sample = derot_base + global_offset*S - residual*S,
	 * and the global offset (measured in the reference canvas) carries the
	 * intersection origin, so the source disk must be placed at
	 * (cx - intersection_x)*S to cancel it — exactly as the single-sequence
	 * mpp_derot_frame_map does (which places one disk for both sides). Using
	 * (0,0,1) for the source instead samples every frame off by the intersection
	 * origin whenever it is cropped (sharp on a full-frame intersection, smeared
	 * otherwise). The source sequence still contributes its OWN disk centre via
	 * srcs[s].derot, so a displaced sequence is still relocated correctly. */
	const double S = std::max(1.0, cfg.drizzle_scale);
	DerotMapProvider derot_provider =
	    [&](int g, int DY, int DX, cv::Mat &mx, cv::Mat &my, cv::Mat &mu) -> bool {
		int s, l;
		if (!layout.locate(g, &s, &l)) return false;
		mpp_derot_frame_map_ms(refd, srcs[s].derot, l, DX, DY,
		                       (double) intersection[2], (double) intersection[0], S,
		                       (double) intersection[2], (double) intersection[0], S,
		                       mx, my, mu, src_masks(s));
		return true;
	};

	gui_iface.set_progress(PROGRESS_RESET, _("Derotation: stacking frames"));
	const auto wr = mpp::stack_warp_apply_streamed(
	    full_provider, N, num_layers, *aps, apq, shifts, offsets, an.brightness,
	    sorted_idx, intersection, cfg, included.data(), nt,
	    provider_thread_safe, /*mem_budget_bytes=*/0, &derot_provider);

	mpp_shift_free(shifts);

	if (wr.oom) { res.oom = true; return res; }
	if (wr.cancelled || wr.image.empty()) { res.cancelled = true; return res; }
	res.image = wr.image;
	return res;
}

MultiStackResult multistack_channel(const std::vector<MsSource> &srcs,
                                    const mpp_derot_t *out_derot, int num_layers,
                                    const mpp_config_t &cfg,
                                    int max_threads, bool provider_thread_safe) {
	MultiStackAnalysis an = multistack_analyze(srcs, out_derot, cfg,
	                                           max_threads, provider_thread_safe);
	if (an.error || an.oom || an.cancelled) {
		MultiStackResult res;
		res.error = an.error;
		res.oom = an.oom;
		res.cancelled = an.cancelled;
		res.num_frames = an.num_frames;
		return res;
	}
	MultiStackResult res = multistack_stack(srcs, out_derot, num_layers, cfg,
	                                        an, *an.aps,
	                                        max_threads, provider_thread_safe);
	mpp_ap_free(an.aps);
	return res;
}

}  // namespace mpp
