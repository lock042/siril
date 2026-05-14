/*
 * Multipoint registration & stacking — PlanetarySystemStacker port.
 * Orchestrator. See pss_port_plan.md Phase 6 for the full design.
 *
 * Public entry points (C-callable, declared in mpp.h):
 *   - mpp_run_alloc / mpp_run_free
 *   - mpp_analyze         (Stage A: rank + global align + avg ref + APs)
 *   - mpp_compute_shifts  (Stage B: per-AP per-frame multilevel correlation)
 *   - mpp_stack_apply     (Stage C: brightness equalise + stack + uint16)
 *   - register_mpp        (Siril `register -method=mpp`: A + B + sidecar)
 *
 * Stages re-read frames from the sequence each time. The trade-off is
 * simplicity over memory; a future optimisation could keep frames in
 * memory across A → B → C for GUI workflows where the user pauses
 * between stages.
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

extern "C" {
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_frames.h"
#include "registration/mpp_shift.h"
#include "registration/mpp_sidecar.h"
#include "registration/registration.h"
}

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_ap_priv.hpp"
#include "registration/mpp_image_priv.hpp"
#include "registration/mpp_rank_priv.hpp"
#include "registration/mpp_shift_priv.hpp"
#include "registration/mpp_stack_priv.hpp"

/* mpp_run_alloc / mpp_run_free live in mpp_run.c (no Siril runtime deps). */

namespace mpp {

namespace {

/* RAII-style fits storage so we don't leak Siril's WORD buffers on early
 * returns. Siril's `fits` is POD with malloc-owned `data` etc.; clearfits
 * frees them. */
struct FitsBuf {
	fits f{};
	~FitsBuf() { clearfits(&f); }
};

/* Helper: read a frame and return a mono analysis Mat (cloned, owns its
 * memory after the fits is cleared). For RGB input, mirrors PSS's
 * frames_mono default (`frames_mono_channel = 'panchromatic'` →
 * `color_index = 3` → `cv::cvtColor(RGB, GRAY)`). Mono input passes
 * through. */
cv::Mat read_analysis_frame(sequence *seq, int idx) {
	FitsBuf buf;
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		return view.empty() ? cv::Mat() : view.clone();
	}
	/* RGB → luminance. Siril's pdata[0/1/2] are R/G/B. */
	std::vector<cv::Mat> planes;
	planes.reserve(n);
	for (int l = 0; l < n; ++l) {
		cv::Mat v = wrap_fits_layer(&buf.f, l);
		if (v.empty()) return cv::Mat();
		planes.push_back(v);
	}
	cv::Mat merged;
	cv::merge(planes, merged);
	cv::Mat gray;
	cv::cvtColor(merged, gray, cv::COLOR_RGB2GRAY);
	return gray;
}

/* Read the frame as a multi-channel cv::Mat preserving all layers (1 or 3).
 * Used by mpp_stack_apply so colour frames carry through to the stacked
 * output. */
cv::Mat read_full_frame(sequence *seq, int idx) {
	FitsBuf buf;
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		if (view.empty()) return cv::Mat();
		return view.clone();
	}
	std::vector<cv::Mat> planes;
	planes.reserve(n);
	for (int l = 0; l < n; ++l) {
		cv::Mat view = wrap_fits_layer(&buf.f, l);
		if (view.empty()) return cv::Mat();
		planes.push_back(view.clone());
	}
	cv::Mat out;
	cv::merge(planes, out);
	return out;
}

/* Reconstruct APQualities from a populated mpp_run_t. Stage B and Stage C
 * both consume this. */
APQualities apq_from_run(const mpp_run_t *run) {
	APQualities apq;
	apq.stack_size = run->stack_size;
	apq.best_frame_indices.assign(run->aps->count, {});
	apq.used_alignment_points.assign(run->num_frames, {});
	for (int a = 0; a < run->aps->count; ++a) {
		apq.best_frame_indices[a].reserve(run->stack_size);
		for (int k = 0; k < run->stack_size; ++k) {
			const int frame_idx = run->best_frame_indices[a * run->stack_size + k];
			apq.best_frame_indices[a].push_back(frame_idx);
			if (frame_idx >= 0 && frame_idx < run->num_frames)
				apq.used_alignment_points[frame_idx].push_back(a);
		}
	}
	return apq;
}

std::vector<FrameOffset> offsets_from_run(const mpp_run_t *run) {
	std::vector<FrameOffset> out(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i) {
		out[i].dy = run->intersection[0] - run->global_shifts[2 * i + 0];
		out[i].dx = run->intersection[2] - run->global_shifts[2 * i + 1];
	}
	return out;
}

cv::Mat mean_frame_from_run(const mpp_run_t *run) {
	return cv::Mat(run->mean_frame_rows, run->mean_frame_cols, CV_32S,
	               run->mean_frame_data);
}

}  // namespace
}  // namespace mpp

/* -------------------- Stage A: mpp_analyze -------------------- */

extern "C" mpp_status_t mpp_analyze(sequence *seq, const mpp_config_t *cfg,
                                    mpp_run_t **run_out) {
	if (!seq || !cfg || !run_out)
		return MPP_EINVAL;
	if (seq->number <= 0)
		return MPP_ENODATA;

	const int N = seq->number;
	std::vector<cv::Mat> analysis_frames;
	std::vector<cv::Mat> blurred_frames;
	analysis_frames.reserve(N);
	blurred_frames.reserve(N);

	std::vector<double> q_rank;
	std::vector<double> frame_brightness;
	q_rank.reserve(N);
	frame_brightness.reserve(N);

	int frame_rows = 0, frame_cols = 0, num_layers = 1, bitpix = 16;
	for (int i = 0; i < N; ++i) {
		const cv::Mat mono = mpp::read_analysis_frame(seq, i);
		if (mono.empty()) return MPP_EIO;
		if (i == 0) { frame_rows = mono.rows; frame_cols = mono.cols; }
		analysis_frames.push_back(mono);
		blurred_frames.push_back(mpp::blur_mono_for_align(mono, *cfg));
		q_rank.push_back(mpp::rank_score_normalized(mono, *cfg));
		frame_brightness.push_back(mpp::rank_average_brightness(mono, *cfg));
	}

	num_layers = seq->nb_layers > 0 ? seq->nb_layers : 1;
	bitpix = seq->bitpix;

	const auto align = mpp::align_global_from_frames(blurred_frames, q_rank, *cfg);
	const auto avg = mpp::align_average_frame(analysis_frames, q_rank,
	                                          align.shifts, *cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, *cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, *cfg);
	if (!aps || aps->count <= 0) {
		mpp_ap_free(aps);
		return MPP_ENODATA;
	}
	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	const auto apq = mpp::ap_compute_frame_qualities(analysis_frames, frame_brightness,
	                                                 *aps, offsets,
	                                                 frame_rows, frame_cols, *cfg);

	mpp_run_t *run = mpp_run_alloc();
	if (!run) { mpp_ap_free(aps); return MPP_ENOMEM; }
	run->cfg = (mpp_config_t *) std::malloc(sizeof(mpp_config_t));
	if (!run->cfg) { mpp_ap_free(aps); mpp_run_free(run); return MPP_ENOMEM; }
	*run->cfg = *cfg;
	run->num_frames = N;
	run->frame_rows = frame_rows;
	run->frame_cols = frame_cols;
	run->num_layers = num_layers;
	run->bitdepth = mpp_bitdepth_from_fits_bitpix(bitpix);

	run->quality          = (double *) std::malloc(N * sizeof(double));
	run->frame_brightness = (double *) std::malloc(N * sizeof(double));
	run->included         = (int *)    std::malloc(N * sizeof(int));
	run->global_shifts    = (int *)    std::malloc(2 * N * sizeof(int));
	if (!run->quality || !run->frame_brightness
	 || !run->included || !run->global_shifts) {
		mpp_run_free(run);
		return MPP_ENOMEM;
	}
	for (int i = 0; i < N; ++i) {
		run->quality[i]          = q_rank[i];
		run->frame_brightness[i] = frame_brightness[i];
		run->included[i]         = 1;
		run->global_shifts[2*i+0] = align.shifts[i][0];
		run->global_shifts[2*i+1] = align.shifts[i][1];
	}
	run->best_frame_idx = align.best_frame_idx;
	for (int k = 0; k < 4; ++k) {
		run->patch_yxyx[k]   = align.patch_yxyx[k];
		run->intersection[k] = avg.intersection[k];
	}

	run->mean_frame_rows = avg.mean_frame.rows;
	run->mean_frame_cols = avg.mean_frame.cols;
	const size_t mean_bytes = (size_t) run->mean_frame_rows
	                        * run->mean_frame_cols * sizeof(int32_t);
	run->mean_frame_data = (int32_t *) std::malloc(mean_bytes);
	if (!run->mean_frame_data) { mpp_run_free(run); return MPP_ENOMEM; }
	std::memcpy(run->mean_frame_data, avg.mean_frame.data, mean_bytes);

	run->aps = aps;
	run->stack_size = apq.stack_size;
	run->best_frame_indices = (int *) std::malloc((size_t) aps->count
	                                              * run->stack_size * sizeof(int));
	if (!run->best_frame_indices) { mpp_run_free(run); return MPP_ENOMEM; }
	for (int a = 0; a < aps->count; ++a)
		for (int k = 0; k < run->stack_size; ++k)
			run->best_frame_indices[a * run->stack_size + k] = apq.best_frame_indices[a][k];

	*run_out = run;
	return MPP_OK;
}

/* -------------------- Stage B: mpp_compute_shifts -------------------- */

extern "C" mpp_status_t mpp_compute_shifts(sequence *seq, const mpp_config_t *cfg,
                                           mpp_run_t *run) {
	if (!seq || !cfg || !run || !run->aps) return MPP_EINVAL;

	std::vector<cv::Mat> blurred_frames;
	blurred_frames.reserve(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i) {
		const cv::Mat mono = mpp::read_analysis_frame(seq, i);
		if (mono.empty()) return MPP_EIO;
		blurred_frames.push_back(mpp::blur_mono_for_align(mono, *cfg));
	}

	const auto apq = mpp::apq_from_run(run);
	const auto offsets = mpp::offsets_from_run(run);
	const cv::Mat mean_raw = mpp::mean_frame_from_run(run);

	mpp_shifts_t *shifts = mpp::stack_compute_shifts(blurred_frames, mean_raw,
	                                                 *run->aps, apq, offsets, *cfg);
	if (!shifts) return MPP_ENOMEM;
	if (run->shifts) mpp_shift_free(run->shifts);
	run->shifts = shifts;
	return MPP_OK;
}

/* -------------------- Stage C: mpp_stack_apply -------------------- */

extern "C" mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                                        const mpp_run_t *run, fits *out) {
	if (!seq || !cfg || !run || !run->aps || !run->shifts || !out)
		return MPP_EINVAL;

	/* Read frames preserving all layers — colour data flows through Stage C
	 * even though Stages A and B use only the green analysis layer. */
	std::vector<cv::Mat> raw_frames;
	raw_frames.reserve(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i) {
		const cv::Mat frame = mpp::read_full_frame(seq, i);
		if (frame.empty()) return MPP_EIO;
		raw_frames.push_back(frame);
	}

	const auto apq = mpp::apq_from_run(run);
	const auto offsets = mpp::offsets_from_run(run);

	std::vector<int> sorted_idx(run->num_frames);
	std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [run](int a, int b) { return run->quality[a] > run->quality[b]; });

	std::vector<double> frame_brightness(run->frame_brightness,
	                                     run->frame_brightness + run->num_frames);
	const cv::Vec4i intersection(run->intersection[0], run->intersection[1],
	                             run->intersection[2], run->intersection[3]);

	const auto loop = mpp::stack_apply_shifts(raw_frames, *run->aps, apq,
	                                          run->shifts, offsets,
	                                          frame_brightness, sorted_idx,
	                                          intersection, *cfg);
	const cv::Mat stacked = mpp::stack_merge_alignment_point_buffers(
	    loop.state, loop.border, *run->aps, *cfg);

	/* Pack stacked into the caller-supplied fits. The caller allocates the
	 * shell (`fits out = {0};`); we fill its members and own `data`. For
	 * multi-channel stacked (CV_16UC3), Siril expects PLANAR layout in
	 * `data` (R plane then G plane then B plane), with pdata[0..2]
	 * pointing into the contiguous buffer. cv::Mat stores interleaved
	 * BGRBGR... — we split + repack. */
	clearfits(out);
	const int H = stacked.rows;
	const int W = stacked.cols;
	const int C = stacked.channels();
	const size_t plane = (size_t) H * W;
	out->data = (WORD *) std::calloc(plane * C, sizeof(WORD));
	if (!out->data) return MPP_ENOMEM;
	if (C == 1) {
		std::memcpy(out->data, stacked.data, plane * sizeof(WORD));
	} else {
		std::vector<cv::Mat> planes;
		cv::split(stacked, planes);
		for (int c = 0; c < C; ++c)
			std::memcpy(out->data + c * plane, planes[c].data, plane * sizeof(WORD));
	}
	out->rx = W;
	out->ry = H;
	out->naxes[0] = W;
	out->naxes[1] = H;
	out->naxes[2] = C;
	out->naxis = (C == 1) ? 2 : 3;
	out->bitpix = USHORT_IMG;
	out->orig_bitpix = USHORT_IMG;
	out->type = DATA_USHORT;
	out->pdata[0] = out->data;
	out->pdata[1] = (C >= 2) ? out->data + plane     : NULL;
	out->pdata[2] = (C >= 3) ? out->data + plane * 2 : NULL;
	return MPP_OK;
}

/* -------------------- Siril register-framework entry point -------------------- */

/* Populate seq->regparam[layer][i].quality from a completed Stage-A run so
 * Siril's existing frame selector / quality plot displays PSS Laplace-σ.
 * Public via mpp.h so the CLI (process_pss, process_register_mpp) and the
 * GUI register_mpp(args) entry can all wire it in. The layer choice is
 * caller-supplied; PSS's quality is panchromatic (same score on every
 * channel), but we only populate one layer to mirror Siril's per-layer
 * regdata model — the existing frame selector reads one layer at a time. */
extern "C" void mpp_write_quality_to_regdata(sequence *seq, int layer,
                                             const mpp_run_t *run) {
	if (!seq || !run || seq->nb_layers <= 0) return;
	if (layer < 0 || layer >= seq->nb_layers) layer = 0;

	if (!seq->regparam) {
		seq->regparam = (regdata **) std::calloc(seq->nb_layers, sizeof(regdata *));
		if (!seq->regparam) return;
	}
	if (!seq->regparam[layer]) {
		seq->regparam[layer] = (regdata *) std::calloc(seq->number, sizeof(regdata));
		if (!seq->regparam[layer]) return;
	}
	for (int i = 0; i < seq->number && i < run->num_frames; ++i) {
		seq->regparam[layer][i].quality = run->quality[i];
		/* Match other methods: leave H as the zero homography (identity is
		 * set by other paths via cvGetEye; mpp shifts are per-AP per-frame
		 * and not expressible as a single global homography). */
	}
	seq->needs_saving = TRUE;
}

extern "C" int register_mpp(struct registration_args *regargs) {
	if (!regargs || !regargs->seq)
		return MPP_EINVAL;

	mpp_config_t cfg;
	if (regargs->mpp_cfg) {
		cfg = *(const mpp_config_t *) regargs->mpp_cfg;
	} else {
		mpp_config_defaults(&cfg);
	}
	cfg.bitdepth = mpp_bitdepth_from_fits_bitpix(regargs->seq->bitpix);

	siril_log_message(_("pss: Stage A — analyse %d frames\n"), regargs->seq->number);
	mpp_run_t *run = NULL;
	const int rc_a = mpp_analyze(regargs->seq, &cfg, &run);
	if (rc_a != MPP_OK) {
		siril_log_message(_("pss: Stage A failed (code %d)\n"), rc_a);
		return rc_a;
	}
	siril_log_message(_("pss: Stage A done — best=%d, %d APs, stack_size=%d\n"),
	                  run->best_frame_idx, run->aps->count, run->stack_size);

	/* Phase 9.2: surface per-frame quality through Siril's regdata so the
	 * existing frame selector and quality plot pick it up. */
	mpp_write_quality_to_regdata(regargs->seq, regargs->layer, run);

	siril_log_message(_("pss: Stage B — per-AP per-frame shifts\n"));
	const int rc_b = mpp_compute_shifts(regargs->seq, &cfg, run);
	if (rc_b != MPP_OK) {
		siril_log_message(_("pss: Stage B failed (code %d)\n"), rc_b);
		mpp_run_free(run);
		return rc_b;
	}
	siril_log_message(_("pss: Stage B done — %d shifts failed\n"),
	                  run->shifts ? run->shifts->failure_counter : -1);

	char sidecar_path[2048];
	snprintf(sidecar_path, sizeof(sidecar_path), "%s.mpp", regargs->seq->seqname);
	if (mpp_sidecar_write(sidecar_path, run) != MPP_OK) {
		siril_log_message(_("pss: sidecar write to %s failed\n"), sidecar_path);
		mpp_run_free(run);
		return MPP_EIO;
	}
	siril_log_message(_("pss: sidecar written to %s\n"), sidecar_path);

	mpp_run_free(run);
	return MPP_OK;
}
