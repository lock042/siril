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
 * Multipoint registration & stacking — orchestrator.
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
#include <queue>
#include <utility>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

/* core/siril.h pulls in <omp.h>, which under GCC 15's libgomp declares
 * C++ templates; those can't sit inside an `extern "C"` block. Include
 * it first so the header guard fires before the wrap below. */
#include "core/siril.h"

extern "C" {
#include "core/siril_log.h"
#include "core/proto.h"
#include "core/gui_iface.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "algos/demosaicing.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/ser.h"

#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_frames.h"
#include "registration/mpp/mpp_shift.h"
#include "registration/mpp/mpp_sidecar.h"
#include "registration/mpp/mpp_derot_sidecar.h"
#include "registration/registration.h"
}

#include "registration/mpp/mpp_derot.h"   /* C++ (cv::Mat) — outside extern "C" */
#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_ap_priv.hpp"
#include "registration/mpp/mpp_image_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_shift_priv.hpp"
#include "registration/mpp/mpp_stack_priv.hpp"

/* mpp_run_alloc / mpp_run_free live in mpp_run.c (no Siril runtime deps). */

/* Analysis-frame cache shared across Stage A → Stage B → recompute.
 * Backed by a vector<cv::Mat> sized to N; populated slots are top-K-by-
 * quality cv::Mats (mono, CV_8U or CV_16U matching cfg.bitdepth);
 * empty slots are cache misses. Empty top-level struct ABI by design —
 * the C side never indexes into it, just stashes the pointer on
 * mpp_run_t and calls mpp_cache_free when the run ends. */
struct mpp_cache {
	std::vector<cv::Mat> raw_frames;  /* size N; empty Mat = miss */
	int bitdepth;                      /* CV_8U if 8, else CV_16U */
	int populated_count;               /* sum of non-empty slots, for logs */
};

extern "C" void mpp_cache_free(mpp_cache_t *cache) {
	delete cache;
}

namespace mpp {

namespace {

/* RAII-style fits storage so we don't leak Siril's WORD buffers on early
 * returns. Siril's `fits` is POD with malloc-owned `data` etc.; clearfits
 * frees them. */
struct FitsBuf {
	fits f{};
	~FitsBuf() { clearfits(&f); }
};

}  /* close inner anonymous namespace */

/* sensor_pattern → OpenCV cvtColor demosaic code, where the pattern
 * describes the buffer rows the Mat actually wraps. OpenCV's 2-letter
 * COLOR_BayerXX2RGB constants describe the 2×2 tile by its second-row
 * first-two pixels (the 4-letter RGGB/BGGR/GBRG/GRBG aliases only exist
 * on OpenCV ≥ 4.x and break the oldstable AppImage build). */
static int bayer_pattern_to_cv_code(sensor_pattern pat) {
	switch (pat) {
		case BAYER_FILTER_RGGB: return cv::COLOR_BayerBG2RGB;
		case BAYER_FILTER_BGGR: return cv::COLOR_BayerRG2RGB;
		case BAYER_FILTER_GBRG: return cv::COLOR_BayerGR2RGB;
		case BAYER_FILTER_GRBG: return cv::COLOR_BayerGB2RGB;
		default: return -1;
	}
}

/* AVI Bayer pattern → OpenCV cvtColor code, with Y-flip awareness so the
 * pattern the user specifies (top-down, as it appears in the AVI raster)
 * stays correct after film_read_frame's fits_flip_top_to_bottom. Returns
 * -1 when no debayer should happen (AUTO/NONE/non-Bayer).
 *
 * We feed in the Mat we got back from the pipeline (bottom-up after the
 * Y-flip), so we ask Siril's helper to adjust the user's top-down pick. */
int avi_bayer_to_cv_code(int avi_pattern, int rows) {
	if (avi_pattern == MPP_AVI_BAYER_AUTO || avi_pattern == MPP_AVI_BAYER_NONE)
		return -1;
	sensor_pattern pat;
	switch (avi_pattern) {
		case MPP_AVI_BAYER_RGGB: pat = BAYER_FILTER_RGGB; break;
		case MPP_AVI_BAYER_BGGR: pat = BAYER_FILTER_BGGR; break;
		case MPP_AVI_BAYER_GBRG: pat = BAYER_FILTER_GBRG; break;
		case MPP_AVI_BAYER_GRBG: pat = BAYER_FILTER_GRBG; break;
		default: return -1;
	}
	adjust_Bayer_pattern_orientation(&pat, (unsigned int) rows, TRUE);
	return bayer_pattern_to_cv_code(pat);
}

/* CFA SER → OpenCV cvtColor code for a raw-mosaic read, -1 for non-Bayer.
 * The header's color_id (already preference-adjusted at open time)
 * describes the buffer in stored order — Siril's own SER debayer hands
 * the same pattern to librtprocess without orientation adjustment, and
 * ser_read_frame doesn't flip rows — so it maps directly. */
static int ser_bayer_to_cv_code(const sequence *seq) {
	if (!seq || seq->type != SEQ_SER || !seq->ser_file)
		return -1;
	switch (seq->ser_file->color_id) {
		case SER_BAYER_RGGB: return bayer_pattern_to_cv_code(BAYER_FILTER_RGGB);
		case SER_BAYER_BGGR: return bayer_pattern_to_cv_code(BAYER_FILTER_BGGR);
		case SER_BAYER_GBRG: return bayer_pattern_to_cv_code(BAYER_FILTER_GBRG);
		case SER_BAYER_GRBG: return bayer_pattern_to_cv_code(BAYER_FILTER_GRBG);
		default: return -1;
	}
}

/* CFA FITS / FITSEQ → OpenCV cvtColor code for the raw mosaic that
 * readfits / fitseq_read_frame return (those readers never debayer — that
 * happens higher up in Siril's open path, which the mpp readers bypass).
 * Delegates pattern detection to Siril's get_validated_cfa_pattern, which
 * reads the BAYERPAT header (and honours the user's debayer preferences /
 * forced pattern), then resolves orientation and Bayer offsets to one of
 * the four canonical 2×2 patterns as it applies to the in-memory buffer —
 * the same buffer wrap_fits_layer wraps and cvtColor reads, so it maps
 * directly with no extra flip. Returns -1 for a mono frame (no BAYERPAT →
 * BAYER_FILTER_NONE), for X-Trans / CMYG patterns cvtColor can't handle,
 * and for already-debayered multi-channel frames (handled by the caller's
 * n>1 path). force_debayer=FALSE so a genuine mono frame is left alone. */
static int fits_bayer_to_cv_code(fits *f) {
	if (!f || f->naxes[2] != 1) return -1;
	const sensor_pattern pat = get_validated_cfa_pattern(f, FALSE, FALSE);
	return bayer_pattern_to_cv_code(pat);
}

/* Probe whether a FITS / FITSEQ sequence is a raw single-layer CFA mosaic
 * by reading frame 0 and running the same detection the per-frame readers
 * use. nb_layers separates mono/CFA (1) from already-RGB (>1) cheaply, but
 * tells mono from CFA only by reading a frame — one read at stack time,
 * negligible against the whole stack, and it keeps the stack's input
 * classification in lock-step with read_full_frame's per-frame debayer so
 * the accumulator's channel count always matches the frames it receives. */
static bool fits_seq_is_cfa(sequence *seq) {
	if (!seq || (seq->type != SEQ_REGULAR && seq->type != SEQ_FITSEQ))
		return false;
	if (seq->nb_layers > 1) return false;   /* already multi-channel */
	FitsBuf probe;
	if (mpp_seq_read_frame(seq, 0, &probe.f, false, 0,
	                       /*ser_force_debayer=*/false) != MPP_OK)
		return false;
	return fits_bayer_to_cv_code(&probe.f) >= 0;
}

namespace {

/* Helper: read a frame and return a mono analysis Mat (cloned, owns its
 * memory after the fits is cleared). RGB input is converted to a
 * panchromatic luminance via cv::cvtColor(RGB, GRAY); mono input passes
 * through.
 *
 * CFA input (a Bayer SER read raw — the analysis stages pin the SER
 * reader to raw via SerAnalysisRawGuard — or a SEQ_AVI with
 * `avi_pattern` set to a Bayer value) is OpenCV-debayered before
 * luminance extraction so the rank/align passes see a properly
 * demosaiced image rather than a CFA-modulated brightness map.
 * cv::cvtColor's linear demosaic is plenty for a luminance-only
 * consumer and costs a fraction of the transient memory and CPU of the
 * RCD path in ser_read_frame (which peaks at ~9× the mosaic size in
 * float scratch per in-flight frame — enough to break the memory
 * pre-flight's in-flight model at high thread counts).
 *
 * When `bitdepth == 8` the result is downconverted to CV_8U. Siril
 * stores 8-bit data as WORD with values still in 0..255, so the
 * convertTo is a pure storage-layout change with no information loss.
 * Halves the per-frame footprint of the analysis cache (decisive for
 * datasets where N · rx · ry · 2B exceeds the user's memory budget
 * but the CV_8U equivalent fits). The downstream pipeline already
 * handles CV_8U: blur_mono_for_align upscales 8-bit input internally
 * (cfg.bitdepth gate), align_average_frame's accumulator runs through
 * CV_32F regardless of input depth, and rank_average_brightness's
 * threshold scales by cfg.bitdepth too. read_full_frame keeps CV_16U
 * so Stage C's stack accumulator retains the full dynamic range. */
cv::Mat read_analysis_frame(sequence *seq, int idx, int avi_pattern, int bitdepth) {
	FitsBuf buf;
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0,
	                       /*ser_force_debayer=*/false) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	cv::Mat out;
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		if (view.empty()) return cv::Mat();
		int code = -1;
		if (seq && seq->type == SEQ_AVI)
			code = mpp::avi_bayer_to_cv_code(avi_pattern, view.rows);
		else if (seq && seq->type == SEQ_SER)
			code = mpp::ser_bayer_to_cv_code(seq);
		else if (seq && (seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ))
			code = mpp::fits_bayer_to_cv_code(&buf.f);
		if (code >= 0) {
			cv::Mat rgb;
			cv::cvtColor(view, rgb, code);
			cv::cvtColor(rgb, out, cv::COLOR_RGB2GRAY);
		} else {
			out = view.clone();
		}
	} else {
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
		cv::cvtColor(merged, out, cv::COLOR_RGB2GRAY);
	}
	if (bitdepth == 8 && out.type() == CV_16U) {
		cv::Mat narrow;
		out.convertTo(narrow, CV_8U);
		return narrow;
	}
	return out;
}

}  /* close inner anonymous namespace */

/* Read the frame as a multi-channel cv::Mat preserving all layers (1 or 3).
 * Used by mpp_stack_apply so colour frames carry through to the stacked
 * output.
 *
 * For SEQ_AVI with a Bayer pattern set, the mono mosaic is debayered to
 * 3-channel RGB before return so the classical stack sees proper colour
 * input. Bayer SERs are RCD-debayered inside ser_read_frame
 * (ser_force_debayer=true overrides the user's debayer-on-open pref) —
 * stacking output keeps the high-quality demosaic. */
cv::Mat read_full_frame(sequence *seq, int idx, int avi_pattern) {
	FitsBuf buf;
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0,
	                       /*ser_force_debayer=*/true) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		if (view.empty()) return cv::Mat();
		int code = -1;
		if (seq && seq->type == SEQ_AVI)
			code = mpp::avi_bayer_to_cv_code(avi_pattern, view.rows);
		else if (seq && (seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ))
			code = mpp::fits_bayer_to_cv_code(&buf.f);
		if (code >= 0) {
			cv::Mat rgb;
			cv::cvtColor(view, rgb, code);
			return rgb;
		}
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

/* See mpp_align_priv.hpp. translation_from_H yields exactly the to-align
 * content translation that shift_fit_from_reg applies, in the same buffer
 * orientation the mpp reader produces — so it maps to the port's (dy, dx)
 * with no extra flip. Kept free of logging so unit tests can call it
 * without a live Siril log/gui.
 *
 * A layer whose translations are all zero is skipped: it carries no
 * alignment information, only quality scores — exactly what
 * mpp_write_quality_to_regdata publishes after every Analyze/Register
 * (zero H on every frame). Feeding that back in as a seed would be
 * actively harmful, not neutral: a seed bypasses the cumulative
 * frame-to-frame tracking chain and pins every frame's search window at
 * its seed value, so an all-zero seed confines the whole search to
 * ±align_frames_search_width around zero shift. Frames drifting further
 * than that get wrong global shifts, and the stack shows doubled
 * features. */
std::vector<cv::Vec2d> seed_from_regdata(sequence *seq) {
	if (!seq || !seq->regparam || seq->number <= 0 || seq->nb_layers <= 0)
		return {};
	for (int l = 0; l < seq->nb_layers; ++l) {
		const regdata *rp = seq->regparam[l];
		if (!rp)
			continue;
		std::vector<cv::Vec2d> seed((size_t) seq->number, cv::Vec2d(0.0, 0.0));
		bool any_translation = false;
		for (int i = 0; i < seq->number; ++i) {
			double dx = 0.0, dy = 0.0;
			translation_from_H(rp[i].H, &dx, &dy);
			if (dx != 0.0 || dy != 0.0)
				any_translation = true;
			seed[(size_t) i] = cv::Vec2d(dy, dx);  /* port (dy, dx) to-align convention */
		}
		if (any_translation)
			return seed;
	}
	return {};
}

namespace {   /* reopen anonymous namespace for the remaining helpers */

/* REGISTER-time per-AP frame count (the baked upper bound): round-up % of N
 * from the register control, clamped to [1, N]. Mirrors the formula in
 * ap_compute_frame_qualities_streamed so the Stage-A provisional matches
 * what the per-AP ranking bakes. (number override wins, though the register
 * GUI/CLI leaves frame_number at -1.) */
int register_stack_size_from_cfg(const mpp_config_t &cfg, int N) {
	if (N <= 0) return 0;
	int s;
	if (cfg.alignment_points_frame_number > 0)
		s = cfg.alignment_points_frame_number;
	else
		s = std::max(1, (int) std::ceil(
		        (double) N * cfg.alignment_points_frame_percent / 100.0));
	return std::min(s, N);
}

/* STACK-time per-AP target from the stack control (number override wins,
 * else round-up % of N), clamped to [1, N]. Stage C clamps this against the
 * register-baked ceiling (run->stack_size). */
int stack_size_from_cfg(const mpp_config_t &cfg, int N) {
	if (N <= 0) return 0;
	int s;
	if (cfg.stack_frame_number > 0)
		s = cfg.stack_frame_number;
	else
		s = std::max(1, (int) std::ceil(
		        (double) N * cfg.stack_frame_percent / 100.0));
	return std::min(s, N);
}

/* Reconstruct APQualities from a populated mpp_run_t. Stage B and Stage C
 * both consume this.
 *
 * The row stride of best_frame_indices is the baked selected_per_ap
 * (stack_size + taper at Analyze time), NOT the current stack_size — a
 * Stack-tab override may have shrunk stack_size after bake, and using it
 * as the stride would scramble the per-AP rows. Selection weights are a
 * pure function of rank, so a shrunk stack_size just re-derives them; the
 * taper is re-clamped against the new target (it can't exceed the room
 * the baked selection actually stored). */
APQualities apq_from_run(const mpp_run_t *run) {
	APQualities apq;
	apq.stack_size = run->stack_size;
	const int stride = run->selected_per_ap > 0 ? run->selected_per_ap
	                                            : run->stack_size;
	int taper = std::min({run->taper, run->stack_size / 2,
	                      run->stack_size - 1,
	                      stride - run->stack_size});
	if (taper < 0) taper = 0;
	apq.taper = taper;
	const int use = std::min(stride, run->stack_size + taper);
	apq.best_frame_indices.assign(run->aps->count, {});
	apq.used_alignment_points.assign(run->num_frames, {});
	for (int a = 0; a < run->aps->count; ++a) {
		apq.best_frame_indices[a].reserve(use);
		for (int k = 0; k < use; ++k) {
			const int frame_idx = run->best_frame_indices[a * stride + k];
			apq.best_frame_indices[a].push_back(frame_idx);
			const float w = stack_selection_weight(k, run->stack_size, taper);
			if (frame_idx >= 0 && frame_idx < run->num_frames && w > 0.0f)
				apq.used_alignment_points[frame_idx].push_back({a, w});
		}
	}
	return apq;
}

std::vector<FrameOffset> offsets_from_run(const mpp_run_t *run) {
	std::vector<FrameOffset> out(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i) {
		out[i].dy = (double) run->intersection[0] - run->global_shifts[2 * i + 0];
		out[i].dx = (double) run->intersection[2] - run->global_shifts[2 * i + 1];
	}
	return out;
}

cv::Mat mean_frame_from_run(const mpp_run_t *run) {
	return cv::Mat(run->mean_frame_rows, run->mean_frame_cols, CV_32S,
	               run->mean_frame_data);
}

}  // namespace
}  // namespace mpp

/* -------------------- Memory pre-flight -------------------- */

/* Each MPP stage may hold all sequence frames in RAM simultaneously
 * (index-stable vector<cv::Mat>) — this is structurally different from
 * Siril's normal sequence ops, which stream frame-by-frame. On planetary
 * workloads (8000+ frames) the cached peak can run to tens of GB.
 *
 * mpp_pick_memory_mode tries four tiers in order of decreasing speed
 * and increasing safety, and returns the most performant mode whose
 * working set fits in get_max_memory_in_MB() (the user's memory-mode
 * preference — RATIO of free RAM, or AMOUNT cap).
 *
 *   CACHED         — hold every frame, often in multiple representations
 *                    (e.g. raw + blurred in Stage A). One disk read pass.
 *   HALF_CACHED    — hold every raw frame, drop the blurred cache and
 *                    re-blur on demand. One disk pass; ~2× the blur
 *                    compute (cheap vs disk).
 *   PARTIAL_CACHED — hold the first K raw frames (K < N) and stream the
 *                    rest from disk. K is sized to fill the available
 *                    budget. Saves K/N of the disk re-reads that pure
 *                    STREAMING would pay; cache misses fall through to
 *                    the same disk path as STREAMING.
 *   STREAMING      — keep only the in-flight frame(s). 2-3× disk reads
 *                    but ~constant memory. Not implemented for every
 *                    stage yet — callers pass half_cached_copies=0 to
 *                    disable PARTIAL too.
 *
 * Each tier's byte cost is N · rx · ry · channels · bytes_per_sample · copies
 * (for cached/half_cached) or K · ... (for partial), or
 * in_flight · rx · ry · channels · bytes_per_sample (for streaming),
 * plus output_bytes (mean_frame accumulator, shifts arrays, drizzle
 * output canvas).
 *
 * `bytes_per_sample` should match what the caller actually stores: 2
 * for CV_16U analysis cache, 1 for CV_8U (8-bit source downconverted
 * by read_analysis_frame).
 *
 * Pass cached_copies = 0 and half_cached_copies = 0 to gate caching off
 * entirely (Stages B and C, which stream by design). Returns MPP_ENOMEM
 * with a clear message if no enabled tier fits. */
typedef enum {
	MPP_MEM_CACHED         = 0,
	MPP_MEM_HALF_CACHED    = 1,
	MPP_MEM_PARTIAL_CACHED = 2,
	MPP_MEM_STREAMING      = 3,
} mpp_mem_mode;

typedef struct {
	int mode;          /* mpp_mem_mode */
	int threads;       /* 1..max_threads, clamped to fit the budget */
	int cached_count;  /* number of raw-frame slots to cache:
	                    * N for CACHED/HALF_CACHED, 1..N-1 for
	                    * PARTIAL_CACHED, 0 for STREAMING. */
	int blurred_count; /* number of blurred-frame slots to pre-compute,
	                    * 0..N. = N in CACHED; for HALF_CACHED, sized
	                    * to the leftover budget so headroom past the
	                    * raw cache pays into a partial blur cache too
	                    * (every blurred hit saves a GaussianBlur in
	                    * Pass 2). 0 for PARTIAL_CACHED and STREAMING. */
} mpp_mem_pick;

static const char *mpp_mem_mode_name(int m) {
	switch (m) {
		case MPP_MEM_CACHED:         return _("cached");
		case MPP_MEM_HALF_CACHED:    return _("half-cached");
		case MPP_MEM_PARTIAL_CACHED: return _("partial-cached");
		case MPP_MEM_STREAMING:      return _("streaming");
		default:                     return "?";
	}
}

/* Picks the most performant memory tier whose working set fits in
 * get_max_memory_in_MB(), and the largest thread count whose
 * in-flight working set also fits. Iterates tiers CACHED →
 * HALF_CACHED → PARTIAL_CACHED → STREAMING; within each tier, picks
 * T in [1, max_threads] inversely proportional to the per-thread cost.
 *
 * Total peak per tier:
 *   CACHED:         cached_copies      · N · F  +  T · per_thread_in_flight · F + output
 *   HALF_CACHED:    half_cached_copies · N · F  +  T · per_thread_in_flight · F + output
 *   PARTIAL_CACHED: half_cached_copies · K · F  +  T · per_thread_in_flight · F + output
 *                   (K = largest k ≤ N fitting at full max_threads)
 *   STREAMING:                                     T · per_thread_in_flight · F + output
 *
 * Per-thread cost reflects what each parallel iteration of the
 * caller's hottest loop holds in scratch (e.g. raw + transient blur
 * during the rank pass = 2). For purely sequential callers, pass
 * max_threads=1, per_thread_in_flight=1.
 *
 * Pass copies=0 to disable a tier the caller doesn't support. Setting
 * half_cached_copies=0 also disables PARTIAL_CACHED (it shares the
 * half-cache layout — raw frames only).
 *
 * bytes_per_sample matches the caller's actual storage: 2 for CV_16U,
 * 1 for CV_8U (8-bit source downconverted). The hot dataset shrinks
 * proportionally — for an 8-bit AVI-origin SER this can lift the
 * pick from PARTIAL to HALF_CACHED outright.
 *
 * Returns MPP_ENOMEM if no enabled tier fits even at T=1. */
static mpp_status_t mpp_pick_memory_mode(int N, int rx, int ry, int channels,
                                         int bytes_per_sample,
                                         int cached_copies,
                                         int half_cached_copies,
                                         int max_threads,
                                         int per_thread_in_flight,
                                         uint64_t output_bytes,
                                         const char *stage_name,
                                         mpp_mem_pick *out) {
	if (N <= 0 || rx <= 0 || ry <= 0 || channels < 1 || !out
	    || max_threads < 1 || per_thread_in_flight < 1
	    || bytes_per_sample < 1)
		return MPP_EINVAL;
	const uint64_t per_frame = (uint64_t) rx * (uint64_t) ry
	                         * (uint64_t) channels * (uint64_t) bytes_per_sample;
	/* Blurred slot size: blur_mono_for_align upscales 8-bit input to
	 * CV_16U before the GaussianBlur, so the cached blurred frame is
	 * always 2 bytes/sample regardless of source bitdepth. For 16-bit
	 * input this equals per_frame; for 8-bit it is 2× per_frame, which
	 * is why the blur top-up must size against blurred_per_frame rather
	 * than per_frame. */
	const uint64_t blurred_per_frame = (uint64_t) rx * (uint64_t) ry
	                                 * (uint64_t) channels * 2;
	const uint64_t budget    = (uint64_t) get_max_memory_in_MB()
	                         * (uint64_t) BYTES_IN_A_MB;
	const uint64_t per_thread_bytes = (uint64_t) per_thread_in_flight * per_frame;

	/* For each tier, the largest T in [1, max_threads] such that
	 *     tier_persistent + T · per_thread_bytes + output ≤ budget.
	 * Returns 0 if even T=1 doesn't fit. */
	auto best_threads_for = [&](uint64_t tier_persistent) -> int {
		const uint64_t fixed = tier_persistent + output_bytes;
		if (fixed >= budget) return 0;
		const uint64_t headroom = budget - fixed;
		const uint64_t max_T_fit = headroom / per_thread_bytes;
		if (max_T_fit < 1) return 0;
		return (int) std::min((uint64_t) max_threads, max_T_fit);
	};

	/* Cache layout (caller-defined via copies arg): one raw analysis
	 * frame plus (copies − 1) blurred frames. Sizes the raw and blurred
	 * components separately because the blurred frame is always CV_16U
	 * (see blurred_per_frame above) while the raw frame depth follows
	 * bytes_per_sample. */
	auto persistent_for = [&](int copies) -> uint64_t {
		if (copies <= 0) return 0;
		return (uint64_t) N * per_frame
		     + (uint64_t)(copies - 1) * (uint64_t) N * blurred_per_frame;
	};
	const uint64_t cached_persistent = persistent_for(cached_copies);
	const uint64_t half_persistent   = persistent_for(half_cached_copies);

	int picked_mode = -1;
	int picked_T = 0;
	int picked_K = 0;
	int picked_K_blur = 0;
	uint64_t picked_bytes = 0;

	if (cached_copies > 0) {
		const int T = best_threads_for(cached_persistent);
		if (T > 0) {
			picked_mode = MPP_MEM_CACHED;
			picked_T = T;
			picked_K = N;
			picked_K_blur = N;
			picked_bytes = cached_persistent + (uint64_t) T * per_thread_bytes + output_bytes;
		}
	}
	if (picked_mode < 0 && half_cached_copies > 0) {
		const int T = best_threads_for(half_persistent);
		if (T > 0) {
			picked_mode = MPP_MEM_HALF_CACHED;
			picked_T = T;
			picked_K = N;
			picked_K_blur = 0;
			picked_bytes = half_persistent + (uint64_t) T * per_thread_bytes + output_bytes;
			/* Blur top-up: HALF_CACHED's chosen persistent is N raw
			 * frames. If the budget has headroom past that and past
			 * the in-flight working set, spend it on a partial blur
			 * cache (each cached blurred slot saves one GaussianBlur
			 * in Pass 2). Sizes against blurred_per_frame, which for
			 * 8-bit input is 2× per_frame because the blur upscales
			 * to CV_16U — using per_frame here is what caused the
			 * blur cache to double-overshoot the budget. cached_copies
			 * must be ≥ 2 for this to make sense — i.e., the CACHED
			 * tier we just rejected was "raw + blurred", so the
			 * caller supports a blur cache. */
			if (cached_copies >= 2) {
				const uint64_t fixed = half_persistent + (uint64_t) T * per_thread_bytes
				                     + output_bytes;
				if (budget > fixed) {
					const uint64_t headroom = budget - fixed;
					const uint64_t K_blur_fit = headroom / blurred_per_frame;
					if (K_blur_fit >= 1) {
						picked_K_blur = (int) std::min((uint64_t) N, K_blur_fit);
						picked_bytes += (uint64_t) picked_K_blur * blurred_per_frame;
					}
				}
			}
		}
	}
	/* Partial cache: HALF_CACHED's layout but K < N frames. Sized to
	 * fill whatever budget remains after max_threads worth of in-flight
	 * frames and the output. Always runs at max_threads since the
	 * persistent cost is what's being capped — adding more threads at
	 * fixed cache size is free given the picker only got here because
	 * K is the binding constraint. K is rounded down to fit. */
	if (picked_mode < 0 && half_cached_copies > 0) {
		const uint64_t fixed = (uint64_t) max_threads * per_thread_bytes + output_bytes;
		if (budget > fixed) {
			const uint64_t headroom = budget - fixed;
			const uint64_t per_cached_frame = (uint64_t) half_cached_copies * per_frame;
			const uint64_t K_fit = headroom / per_cached_frame;
			if (K_fit >= 1) {
				const int K = (int) std::min((uint64_t) (N - 1), K_fit);
				if (K >= 1) {
					picked_mode = MPP_MEM_PARTIAL_CACHED;
					picked_T = max_threads;
					picked_K = K;
					picked_K_blur = 0;
					picked_bytes = (uint64_t) half_cached_copies * (uint64_t) K * per_frame
					             + (uint64_t) picked_T * per_thread_bytes + output_bytes;
				}
			}
		}
	}
	if (picked_mode < 0) {
		const int T = best_threads_for(0);
		if (T > 0) {
			picked_mode = MPP_MEM_STREAMING;
			picked_T = T;
			picked_K = 0;
			picked_K_blur = 0;
			picked_bytes = (uint64_t) T * per_thread_bytes + output_bytes;
		}
	}

	if (picked_mode < 0) {
		/* Nothing fits, not even single-threaded streaming. Report the
		 * absolute minimum (one in-flight frame + output) as the need
		 * figure so the user sees an achievable target. */
		const uint64_t need = per_thread_bytes + output_bytes;
		gchar *need_str  = g_format_size_full(need, G_FORMAT_SIZE_IEC_UNITS);
		gchar *avail_str = g_format_size_full(budget, G_FORMAT_SIZE_IEC_UNITS);
		siril_log_error(_("%s: not enough memory — even a single in-flight "
		                  "frame (%s for %dx%dx%d) exceeds the %s budget. "
		                  "Use a smaller crop or increase the memory "
		                  "ratio in Preferences.\n"),
		                stage_name, need_str, rx, ry, channels, avail_str);
		g_free(need_str);
		g_free(avail_str);
		return MPP_ENOMEM;
	}

	out->mode          = picked_mode;
	out->threads       = picked_T;
	out->cached_count  = picked_K;
	out->blurred_count = picked_K_blur;

	/* Log when we drop below CACHED (the fastest tier) or when threads
	 * are clamped — both have user-visible perf implications. */
	const bool dropped_tier = (cached_copies > 0 && picked_mode != MPP_MEM_CACHED);
	const bool clamped_threads = (picked_T < max_threads);
	if (dropped_tier || clamped_threads) {
		gchar *use_str   = g_format_size_full(picked_bytes, G_FORMAT_SIZE_IEC_UNITS);
		gchar *avail_str = g_format_size_full(budget,       G_FORMAT_SIZE_IEC_UNITS);
		if (dropped_tier) {
			gchar *cached_str = g_format_size_full(
			    cached_persistent + (uint64_t) max_threads * per_thread_bytes + output_bytes,
			    G_FORMAT_SIZE_IEC_UNITS);
			if (picked_mode == MPP_MEM_PARTIAL_CACHED) {
				siril_log_message(_("%s: partial-cached mode, %d/%d frames cached, "
				                    "%d thread(s) (cached would need %s, %s available, "
				                    "using %s)\n"),
				                  stage_name, picked_K, N, picked_T,
				                  cached_str, avail_str, use_str);
			} else if (picked_mode == MPP_MEM_HALF_CACHED && picked_K_blur > 0) {
				siril_log_message(_("%s: half-cached mode + %d/%d blurred top-up, "
				                    "%d thread(s) (cached would need %s, %s available, "
				                    "using %s)\n"),
				                  stage_name, picked_K_blur, N, picked_T,
				                  cached_str, avail_str, use_str);
			} else {
				siril_log_message(_("%s: %s mode, %d thread(s) "
				                    "(cached would need %s, %s available, using %s)\n"),
				                  stage_name, mpp_mem_mode_name(picked_mode),
				                  picked_T, cached_str, avail_str, use_str);
			}
			g_free(cached_str);
		} else {
			siril_log_message(_("%s: %d/%d threads (%s used, %s available)\n"),
			                  stage_name, picked_T, max_threads, use_str, avail_str);
		}
		g_free(use_str);
		g_free(avail_str);
	}
	return MPP_OK;
}

/* Worst-case transient cost of one analysis-read iteration
 * (read_analysis_frame + blur/rank scratch), expressed in per-frame
 * slots for mpp_pick_memory_mode's per_thread_in_flight argument
 * (slot = rx·ry·bytes_per_sample bytes).
 *
 * Peaks, in bytes/pixel (WORD fits storage; cv Mats CV_16U until the
 * final 8-bit narrow):
 *   mono: fits 2 + cloned mono 2 + blur/Laplace scratch ~4        ≈ 8
 *   CFA:  fits mosaic 2 + cvtColor RGB 6 + gray 2 + scratch ~4    ≈ 14
 *   RGB:  ser_read_frame tmp-copy peak 12, then fits 6 + merged 6
 *         + gray 2 + scratch ~4                                   ≈ 18
 * Deliberately a little fat: overestimating costs a few Pass-1 threads
 * when the budget is tight; underestimating means allocation failure
 * mid-run — on Windows (no overcommit) that's how the old flat
 * 2-slot model crashed Stage A on Bayer SERs. */
static int mpp_analysis_in_flight_slots(const sequence *seq, const mpp_config_t *cfg) {
	int peak_bpp = 8;   /* mono */
	const bool ser_cfa = (seq->type == SEQ_SER && seq->ser_file
	                      && seq->ser_file->color_id >= SER_BAYER_RGGB
	                      && seq->ser_file->color_id <= SER_BAYER_BGGR);
	const bool avi_cfa = (seq->type == SEQ_AVI
	                      && cfg->avi_bayer_pattern >= MPP_AVI_BAYER_RGGB
	                      && cfg->avi_bayer_pattern <= MPP_AVI_BAYER_GRBG);
	/* Single-layer FITS / FITSEQ may be a raw CFA mosaic that
	 * read_analysis_frame now cvtColor-debayers (fits_bayer_to_cv_code).
	 * We can't tell mono from CFA without reading a frame here, so assume
	 * CFA and budget for the RGB+gray transient — over-budgeting a mono
	 * FITS only costs a few Pass-1 threads under memory pressure, whereas
	 * under-budgeting a CFA one risks the Windows no-overcommit OOM the
	 * 8-bpp mono estimate used to cause on Bayer SERs. */
	const bool fits_maybe_cfa = (seq->type == SEQ_REGULAR
	                             || seq->type == SEQ_FITSEQ)
	                            && seq->nb_layers <= 1;
	const bool rgb = (seq->type == SEQ_SER && seq->ser_file
	                  && (seq->ser_file->color_id == SER_RGB
	                      || seq->ser_file->color_id == SER_BGR))
	              || (seq->type != SEQ_SER && seq->nb_layers > 1);
	if (ser_cfa || avi_cfa || fits_maybe_cfa)
		peak_bpp = 14;
	else if (rgb)
		peak_bpp = 18;
	const int bps = (cfg->bitdepth == 8) ? 1 : 2;
	return (peak_bpp + bps - 1) / bps;
}

/* RAII helper: temporarily pins a CFA SER to raw-mosaic reads by setting
 * debayer_type_ser = SER_MONO for the duration of an analysis stage,
 * restoring on destruction. Analysis (rank / global align / per-AP
 * quality) only consumes luminance, so read_analysis_frame demosaics the
 * mosaic with cv::cvtColor — a fraction of the transient memory and CPU
 * of ser_read_frame's RCD path, and identical regardless of the user's
 * debayer-on-open preference. The full-frame stacking reads are
 * unaffected: read_full_frame passes ser_force_debayer=true, which
 * overrides debayer_type_ser inside ser_read_frame. */
namespace {
struct SerAnalysisRawGuard {
	sequence *seq = nullptr;
	ser_color saved_debayer_type = SER_MONO;
	bool active = false;

	~SerAnalysisRawGuard() {
		if (active && seq && seq->ser_file)
			seq->ser_file->debayer_type_ser = saved_debayer_type;
	}
};

bool maybe_engage_ser_raw_for_analysis(SerAnalysisRawGuard &g, sequence *seq) {
	if (!seq || seq->type != SEQ_SER || !seq->ser_file) return false;
	const ser_color cid = seq->ser_file->color_id;
	if (cid < SER_BAYER_RGGB || cid > SER_BAYER_BGGR) return false;
	if (seq->ser_file->debayer_type_ser == SER_MONO) return false; /* already raw */
	g.seq = seq;
	g.saved_debayer_type = seq->ser_file->debayer_type_ser;
	g.active = true;
	seq->ser_file->debayer_type_ser = SER_MONO;
	return true;
}

/* Run a provider body, converting OpenCV/std allocation failures into an
 * empty-Mat miss — every streamed consumer already skips empty frames.
 * Without this, a cv::Exception / std::bad_alloc thrown inside an OpenMP
 * region calls std::terminate(): instant abort, no message, no backtrace.
 * Windows hits this for real — there is no memory overcommit, so
 * cv::Mat::create throws as soon as the commit charge is exhausted. */
gint provider_oom_warned = 0;
template <typename Fn>
cv::Mat guarded_frame(Fn &&fn) {
	try {
		return fn();
	} catch (const std::exception &e) {
		if (g_atomic_int_compare_and_exchange(&provider_oom_warned, 0, 1))
			siril_log_error(_("mpp: frame decode/processing failed (%s) — "
			                  "affected frames will be skipped. If this is "
			                  "an out-of-memory condition, reduce the memory "
			                  "ratio in Preferences.\n"), e.what());
		return cv::Mat();
	}
}
}  // namespace

/* Derotation-aware registration (5b): when a .derot plan is active, every
 * analysis frame is warped to the reference epoch before global alignment,
 * mean-frame building and per-AP shift measurement. The disk position is
 * rotation-invariant, so global alignment still recovers per-frame jitter on
 * the derotated frames; the mean is then sharp (not rotationally smeared) and
 * the per-AP shifts are small epoch-space residuals. Stage C composes the same
 * derotation with those residuals in one resample, so the rotation is applied
 * exactly once and is not double-counted. Disk fit is in full-frame pixels, so
 * the map is sized to the frame itself (org 0,0, scale 1). */
/* A .derot plan applies to this run only if its frame count AND the frame
 * dimensions its disk fit was measured in both match — otherwise the disk would
 * be mis-placed (e.g. a plan from a different-resolution capture with the same
 * frame count). */
static bool mpp_derot_applies(const struct mpp_derot *d, int num_frames,
                              int rows, int cols) {
	return d && d->num_frames == num_frames
	        && d->frame_rows == rows && d->frame_cols == cols;
}

static cv::Mat mpp_derot_warp_frame(const struct mpp_derot *d, int frame,
                                    const cv::Mat &m) {
	if (!d || m.empty()) return m;
	cv::Mat mapx, mapy, mu, out;
	mpp_derot_frame_map(d, frame, m.cols, m.rows, 0.0, 0.0, 1.0, mapx, mapy, mu);
	cv::remap(m, out, mapx, mapy, cv::INTER_LANCZOS4, cv::BORDER_CONSTANT,
	          cv::Scalar(0));
	return out;
}

/* -------------------- Stage A: mpp_analyze -------------------- */

/* Driver state for the sub-range scaled progress callback. The C++
 * helpers report [0, 1] within their own stage; the callback maps that
 * onto a slice of the full Analyze 0..1 bar so the bar moves continuously
 * across the whole flow. */
struct stage_progress_range {
	double base;    /* fraction at stage start */
	double span;    /* fraction width allotted to this stage */
};

static void stage_progress_cb(double fraction, void *user) {
	const stage_progress_range *r = (const stage_progress_range *) user;
	gui_iface.set_progress(r->base + r->span * fraction, NULL);
}

/* Allocation of the 0..1 bar across the Analyze stages. The percentages
 * are rough estimates of relative wall-clock time. AP placement is a
 * single-image op that finishes essentially instantly — it shares its
 * fraction with the average build (no separate slice, just a text
 * update). */
/* Stage A sub-phase boundaries on its own 0..1 progress bar. Per-AP
 * ranking was deferred to Stage B (mpp_recompute_qualities runs its own
 * 0..1 bar), so the three remaining passes — read, global-align, average —
 * fill the whole bar; AP placement that follows is a fast single-image step
 * left at 100%. */
#define MPP_PROG_READ_END         0.60   /* parallel frame read */
#define MPP_PROG_GLOBAL_END       0.90   /* per-frame correlation */
#define MPP_PROG_AVERAGE_END      1.00   /* per-frame accumulator */

static mpp_status_t mpp_analyze_impl(sequence *seq, const mpp_config_t *cfg,
                                     mpp_run_t **run_out,
                                     const struct mpp_derot *derot) {
	if (!seq || !cfg || !run_out)
		return MPP_EINVAL;
	if (seq->number <= 0)
		return MPP_ENODATA;
	const bool derot_active = mpp_derot_applies(derot, seq->number,
	                                            seq->ry, seq->rx);
	if (derot && !derot_active)
		siril_log_warning(_("Analyze: the .derot plan does not match this "
		                    "sequence (frame count or dimensions); derotation "
		                    "skipped.\n"));
	if (derot_active)
		siril_log_message(_("Analyze: derotation active — registering in the "
		                    "reference-epoch frame.\n"));

	/* Pin CFA SERs to raw-mosaic reads for the whole stage — all passes
	 * and the persistent analysis cache they populate go through
	 * read_analysis_frame's cv::cvtColor demosaic, regardless of the
	 * user's debayer-on-open preference. */
	SerAnalysisRawGuard raw_guard;
	maybe_engage_ser_raw_for_analysis(raw_guard, seq);

	const int N = seq->number;

	/* Honour the sequence's frame-inclusion state from the very start.
	 * Frames the user excluded in the frame selector are kept out of the
	 * global reference pick, the averaged reference frame, AP placement and
	 * per-AP ranking — not just the final stack — so a corrupt or
	 * out-of-frame frame can't poison the shared analysis artifacts. Real
	 * per-frame quality is still computed for every frame below (the Plot
	 * tab shows it); only the selection decisions consult the mask. */
	std::vector<int> included(N, 1);
	if (seq->imgparam)
		for (int i = 0; i < N; ++i)
			included[i] = seq->imgparam[i].incl ? 1 : 0;
	int n_included = 0;
	for (int i = 0; i < N; ++i) n_included += included[i];
	if (n_included == 0) {
		siril_log_error(_("Analyze: every frame is excluded — include at "
		                  "least one frame in the frame selector.\n"));
		return MPP_ENODATA;
	}
	if (n_included < N)
		siril_log_message(_("mpp: %d of %d frames excluded by the frame "
		                    "selector; analysing the remaining %d\n"),
		                  N - n_included, N, n_included);

	/* Analysis is always mono (RGB inputs are luminance-converted in
	 * read_analysis_frame). Memory tiers:
	 *   CACHED         — raw + blurred (2 copies); 1 disk pass.
	 *   HALF_CACHED    — raw only; Pass 2 re-blurs on demand. 1 disk pass.
	 *   PARTIAL_CACHED — first K of N raw frames cached; Passes 2/3/4 hit
	 *                    cache for [0..K) and re-read from disk for [K..N).
	 *                    1 + (1-K/N)·(passes) disk passes.
	 *   STREAMING      — neither cached; Passes 2/3/4 re-read from disk.
	 *                    ~3 disk passes total but only T frames in RAM.
	 * Pass 1 (rank + brightness) is parallel in all four tiers; the picker
	 * clamps T down to whatever fits. Passes 2/3/4 stay sequential in
	 * P3/P4 (single-frame accumulator); P2's two-section forward/backward
	 * parallelism only engages when the provider is thread-safe (SER or
	 * reentrant-FITS, including cached vectors), which still holds with
	 * PARTIAL_CACHED on SER because both cache hits and disk reads are
	 * thread-safe there.
	 * per_thread_in_flight comes from mpp_analysis_in_flight_slots: the
	 * peak transient of one read+rank iteration, which for CFA/RGB input
	 * is several times the bare mono frame (decode intermediates) — the
	 * old flat "2 mono frames" model under-budgeted Bayer SERs badly
	 * enough to fail allocations mid-run on Windows.
	 * bytes_per_sample tracks read_analysis_frame's storage choice: 8-bit
	 * sources are downconverted to CV_8U (1 B/sample), so the picker
	 * sees half the cache pressure for AVI / 8-bit-SER datasets and can
	 * land them in HALF_CACHED outright. */
#ifdef _OPENMP
	const int max_threads = std::max(1, com.max_thread);
#else
	const int max_threads = 1;
#endif
	const int analysis_bytes_per_sample = (cfg->bitdepth == 8) ? 1 : 2;
	mpp_mem_pick mem_pick{MPP_MEM_CACHED, max_threads, N, 0};
	{
		const mpp_status_t mem = mpp_pick_memory_mode(
		    N, (int) seq->rx, (int) seq->ry, /*channels=*/1,
		    analysis_bytes_per_sample,
		    /*cached_copies=*/2,
		    /*half_cached_copies=*/1,
		    max_threads,
		    mpp_analysis_in_flight_slots(seq, cfg),
		    /*output_bytes=*/0,
		    _("Analyze"), &mem_pick);
		if (mem != MPP_OK) return mem;
	}
	const int cached_count   = mem_pick.cached_count;  /* 0..N, see picker */
	const int blurred_count  = mem_pick.blurred_count; /* 0..N, see picker */
	const int pass1_threads  = mem_pick.threads;

	/* Index-stable storage so parallel threads can write to per-frame slots
	 * without contention. push_back would force serialisation.
	 *
	 * analysis_frames (raw cache): sized to N whenever cached_count > 0.
	 * Slots populated during Pass 1 hold top-K-by-quality frames (rolling
	 * eviction below); unpopulated slots stay empty (sentinel for
	 * "re-read from disk"). With K = N (full HALF/CACHED) every slot
	 * is filled; with K < N (PARTIAL) the chosen slots are sparse.
	 * Top-K-by-quality (over first-K-by-index) matters for the
	 * persistent cache reused by Stage B / recompute: those passes
	 * preferentially read top-quality frames, so a quality-aligned
	 * cache pushes their hit rate toward 100% even at modest K/N.
	 *
	 * blurred_frames (blurred cache): sized to N whenever blurred_count > 0.
	 * Slots [0..blurred_count) are pre-computed during Pass 1; the rest
	 * stay empty. CACHED tier has blurred_count = N (every blurred slot
	 * filled); HALF_CACHED may have 0 < blurred_count < N when the
	 * budget has headroom past the raw cache (each cached blurred slot
	 * saves one GaussianBlur in Pass 2). Pass 2's provider checks the
	 * blurred slot first, then the raw slot (re-blur), then disk. */
	const bool any_cache = (cached_count > 0);
	const bool any_blur  = (blurred_count > 0);
	std::vector<cv::Mat> analysis_frames;
	if (any_cache) analysis_frames.resize(N);
	std::vector<cv::Mat> blurred_frames;
	if (any_blur)  blurred_frames.resize(N);
	std::vector<double> q_rank(N, 0.0);
	std::vector<double> frame_brightness(N, 0.0);

	int frame_rows = (int) seq->ry, frame_cols = (int) seq->rx;
	int num_layers = 1, bitpix = 16;
	int cur_nb = 0;
	int read_failed = 0;
	int alloc_failed = 0;

	/* Rolling min-heap of (quality, frame_idx) for cached slots.
	 * `cached_count` < N triggers eviction: when a newly-ranked frame
	 * beats the worst cached one, evict the worst and insert the new.
	 * Mutex-guarded for parallel Pass 1; AVI is serial so the lock is
	 * uncontended there; SER under 16 threads sees ~µs contention
	 * which is negligible vs the 5-15 ms of per-frame work (read +
	 * blur + Laplace) outside the lock.
	 *
	 * When cached_count == N every frame fits — the heap fills and
	 * never evicts, equivalent to the prior unconditional fill. */
	GMutex cache_mtx;
	g_mutex_init(&cache_mtx);
	using QIdx = std::pair<double, int>;
	std::priority_queue<QIdx, std::vector<QIdx>, std::greater<QIdx>> low_q;
	auto offer_to_cache = [&](int idx, double q, const cv::Mat &mono) {
		if (!any_cache) return;
		g_mutex_lock(&cache_mtx);
		if ((int) low_q.size() < cached_count) {
			analysis_frames[idx] = mono;
			low_q.push({q, idx});
		} else if (q > low_q.top().first) {
			const int evict = low_q.top().second;
			low_q.pop();
			analysis_frames[evict].release();
			analysis_frames[idx] = mono;
			low_q.push({q, idx});
		}
		g_mutex_unlock(&cache_mtx);
	};

	/* Parallel frame read: SER's fd_lock serialises only the actual fread
	 * call inside ser_read_frame so the rest (decompress, debayer, cvtColor,
	 * GaussianBlur, Laplace) runs concurrently. Pattern + guard mirror
	 * shift_methods.c:432 — anything fits-backed needs cfitsio's
	 * fits_is_reentrant() flag. */
	int cancelled = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(pass1_threads) schedule(guided) \
    if (seq->type == SEQ_SER || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ \
                                  || seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
	for (int i = 0; i < N; ++i) {
		if (g_atomic_int_get(&read_failed) || g_atomic_int_get(&cancelled)
		    || g_atomic_int_get(&alloc_failed)) continue;
		if (!processing_should_continue()) {
			g_atomic_int_set(&cancelled, 1);
			continue;
		}
		/* Contain allocation failures (cv::Exception / std::bad_alloc):
		 * an exception escaping an OpenMP structured block is
		 * std::terminate(). On Windows allocations genuinely fail under
		 * memory pressure (no overcommit), so this is a survival path,
		 * not paranoia. */
		try {
			const cv::Mat mono = mpp::read_analysis_frame(seq, i, cfg->avi_bayer_pattern,
			                                              cfg->bitdepth);
			if (mono.empty()) {
				g_atomic_int_set(&read_failed, 1);
				continue;
			}
			if (any_blur && i < blurred_count)
				blurred_frames[i] = mpp::blur_mono_for_align(mono, *cfg);
			const double q = mpp::rank_score_normalized(mono, *cfg);
			q_rank[i]          = q;
			frame_brightness[i] = mpp::rank_average_brightness(mono, *cfg);
			offer_to_cache(i, q, mono);
		} catch (const std::exception &) {
			g_atomic_int_set(&alloc_failed, 1);
			continue;
		}
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(MPP_PROG_READ_END * (double) cur_nb / (double) N, NULL);
	}
	g_mutex_clear(&cache_mtx);
	if (cancelled) return MPP_EINTR;
	if (alloc_failed) {
		siril_log_error(_("Analyze: out of memory while reading frames — "
		                  "reduce the memory ratio in Preferences or use a "
		                  "smaller crop.\n"));
		return MPP_ENOMEM;
	}
	if (read_failed) return MPP_EIO;

	num_layers = seq->nb_layers > 0 ? seq->nb_layers : 1;
	bitpix = seq->bitpix;

	/* Unified providers for Passes 2/3/4. Hit the analysis_frames slot
	 * if populated; otherwise re-read from disk. Works identically for
	 * HALF_CACHED (every slot hits), PARTIAL_CACHED (first K hit), and
	 * STREAMING (every slot misses since analysis_frames is empty). */
	auto raw_provider = [&analysis_frames, seq, cfg, derot, derot_active](int i) -> cv::Mat {
		return guarded_frame([&]() -> cv::Mat {
			cv::Mat m;
			if (i >= 0 && i < (int) analysis_frames.size()
			 && !analysis_frames[i].empty())
				m = analysis_frames[i];
			else
				m = mpp::read_analysis_frame(seq, i, cfg->avi_bayer_pattern, cfg->bitdepth);
			return derot_active ? mpp_derot_warp_frame(derot, i, m) : m;
		});
	};
	/* When derotation is active the pre-blurred cache (blurred raw frames) no
	 * longer matches — the frame must be derotated before blurring — so that
	 * cache is bypassed and each blurred frame is built fresh. */
	auto blurred_provider = [&analysis_frames, &blurred_frames, seq, cfg, derot, derot_active](int i) -> cv::Mat {
		return guarded_frame([&]() -> cv::Mat {
			if (!derot_active && i >= 0 && i < (int) blurred_frames.size()
			 && !blurred_frames[i].empty())
				return blurred_frames[i];
			cv::Mat m;
			if (i >= 0 && i < (int) analysis_frames.size()
			 && !analysis_frames[i].empty())
				m = analysis_frames[i];
			else
				m = mpp::read_analysis_frame(seq, i, cfg->avi_bayer_pattern, cfg->bitdepth);
			if (m.empty()) return cv::Mat();
			if (derot_active) m = mpp_derot_warp_frame(derot, i, m);
			return mpp::blur_mono_for_align(m, *cfg);
		});
	};

	gui_iface.set_progress(MPP_PROG_READ_END, _("Analyze: aligning frames globally"));
	stage_progress_range gp_global{MPP_PROG_READ_END,
	                               MPP_PROG_GLOBAL_END - MPP_PROG_READ_END};
	/* Pass 2 — global align. The blurred provider's lookup order is:
	 *   1. blurred_frames cache (pre-computed in Pass 1; covers
	 *      [0..blurred_count). Best case — no blur, no I/O.
	 *   2. analysis_frames raw cache → blur on demand. HALF_CACHED's
	 *      common case past the blurred top-up; also PARTIAL's hit path.
	 *   3. Disk read → blur on demand. STREAMING and the PARTIAL miss tail.
	 * The forward+backward sweeps parallelise when the provider is
	 * reentrant; cache hits are trivially reentrant, and disk-read
	 * misses are reentrant for SER / reentrant-FITS. */
	const bool streaming_provider_safe =
	    (seq->type == SEQ_SER
	     || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ
	          || seq->type == SEQ_INTERNAL) && fits_is_reentrant()));
	/* Auto-seed the global aligner from any existing .seq shift registration
	 * (Siril enhancement; opt-out via cfg->align_frames_seed_from_regdata).
	 * Lets jumpy video whose gross shifts exceed align_frames_search_width
	 * still align — the bare MLC search only covers the seed's residual. */
	std::vector<cv::Vec2d> seed;
	if (cfg->align_frames_seed_from_regdata) {
		seed = mpp::seed_from_regdata(seq);
		if (!seed.empty())
			siril_log_message(_("mpp: seeding global alignment from existing "
			                    "registration data\n"));
	}
	/* Masked quality for selection: excluded frames sit below the lowest
	 * real score (≥ 0) so they can never be picked as the global reference
	 * or averaged into the mean frame. The unmasked q_rank is kept for
	 * display (run->quality / the Plot tab). */
	std::vector<double> q_rank_sel = q_rank;
	if (n_included < N)
		for (int i = 0; i < N; ++i)
			if (!included[i]) q_rank_sel[i] = -1.0;
	const auto align = mpp::align_global_from_provider(
	    blurred_provider, N, q_rank_sel, *cfg, stage_progress_cb, &gp_global,
	    streaming_provider_safe, seed);
	/* Drop the blurred cache now — Passes 3/4 don't use it, and the
	 * freed bytes give them more headroom (matters most when the blur
	 * top-up was large). */
	blurred_frames.clear();
	blurred_frames.shrink_to_fit();
	if (align.oom) {
		siril_log_error(_("Analyze: out of memory during global alignment — "
		                  "reduce the memory ratio in Preferences.\n"));
		return MPP_ENOMEM;
	}
	if (align.best_frame_idx < 0) return MPP_EINTR;   /* cancellation sentinel */
	gui_iface.set_progress(MPP_PROG_GLOBAL_END, _("Analyze: building reference frame"));
	stage_progress_range gp_avg{MPP_PROG_GLOBAL_END,
	                            MPP_PROG_AVERAGE_END - MPP_PROG_GLOBAL_END};
	/* Pass 3 — average frame. raw_provider serves cache hits or disk
	 * reads transparently; the streamed variant is sequential so thread
	 * safety doesn't apply. */
	const auto avg = mpp::align_average_frame_streamed(
	    raw_provider, N, frame_rows, frame_cols, q_rank_sel, align.shifts, *cfg,
	    stage_progress_cb, &gp_avg, included);
	if (avg.mean_frame.empty()) return MPP_EINTR;
	gui_iface.set_progress(PROGRESS_NONE, _("Analyze: placing alignment points"));
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, *cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, *cfg);
	if (!aps || aps->count <= 0) {
		mpp_ap_free(aps);
		return MPP_ENODATA;
	}
	/* Per-AP frame ranking (formerly "Pass 4") is deferred to Stage B entry
	 * via mpp_recompute_qualities. Analyze stops here — as soon as the AP
	 * grid exists — so the user can review/edit APs without first paying a
	 * full-sequence ranking pass that any edit would invalidate anyway. */

	mpp_run_t *run = mpp_run_alloc();
	if (!run) { mpp_ap_free(aps); return MPP_ENOMEM; }
	/* Attach aps to the run immediately so every subsequent OOM path
	 * below (which calls mpp_run_free) reclaims it too. */
	run->aps = aps;
	run->cfg = (mpp_config_t *) std::malloc(sizeof(mpp_config_t));
	if (!run->cfg) { mpp_run_free(run); return MPP_ENOMEM; }
	*run->cfg = *cfg;
	run->num_frames = N;
	run->frame_rows = frame_rows;
	run->frame_cols = frame_cols;
	run->num_layers = num_layers;
	run->bitdepth = mpp_bitdepth_from_fits_bitpix(bitpix);

	run->quality          = (double *) std::malloc(N * sizeof(double));
	run->frame_brightness = (double *) std::malloc(N * sizeof(double));
	run->included         = (int *)    std::malloc(N * sizeof(int));
	run->global_shifts    = (double *) std::malloc(2 * N * sizeof(double));
	if (!run->quality || !run->frame_brightness
	 || !run->included || !run->global_shifts) {
		mpp_run_free(run);
		return MPP_ENOMEM;
	}
	for (int i = 0; i < N; ++i) {
		run->quality[i]          = q_rank[i];
		run->frame_brightness[i] = frame_brightness[i];
		run->included[i]         = included[i];
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

	/* Per-AP ranking deferred to Stage B (see above). best_frame_indices
	 * stays null until mpp_recompute_qualities fills it; stack_size is a
	 * provisional value from the same cfg formula Stage B/C use (over the
	 * included frames), so the Analyze log and any pre-Register display show
	 * the right number. */
	run->stack_size = mpp::register_stack_size_from_cfg(*cfg, n_included);
	run->taper = 0;
	run->selected_per_ap = 0;
	run->best_frame_indices = nullptr;

	/* Hand the analysis-frame cache off to the run so Stage B / recompute
	 * can hit it. Zero-copy: cv::Mat headers move, pixel data refcounts
	 * stay put. Skipped for STREAMING (analysis_frames is empty) and
	 * harmless otherwise — the cache lives until mpp_run_free. */
	if (any_cache) {
		int populated = 0;
		for (const auto &m : analysis_frames) if (!m.empty()) ++populated;
		if (populated > 0) {
			auto *cache = new mpp_cache;
			cache->raw_frames     = std::move(analysis_frames);
			cache->bitdepth       = cfg->bitdepth;
			cache->populated_count = populated;
			run->cache = cache;
		}
	}

	/* Stamp the fingerprint of the derotation plan the AP geometry was
	 * registered in (0 when derotation was off) so Stage C can detect a
	 * mismatched/added/edited .derot. */
	run->derot_fingerprint = derot_active ? mpp_derot_fingerprint(derot) : 0;

	*run_out = run;
	return MPP_OK;
}

/* Exception boundary: mpp_analyze is called from C, and most cv:: work
 * not already contained inside the OpenMP loops above runs in serial
 * sections whose exceptions would otherwise unwind into the C caller
 * (undefined behaviour / terminate). Out-of-memory is the realistic
 * thrower — map it to MPP_ENOMEM. */
extern "C" mpp_status_t mpp_analyze(sequence *seq, const mpp_config_t *cfg,
                                    mpp_run_t **run_out,
                                    const struct mpp_derot *derot) {
	try {
		return mpp_analyze_impl(seq, cfg, run_out, derot);
	} catch (const std::exception &e) {
		siril_log_error(_("Analyze: failed (%s) — likely out of memory; "
		                  "reduce the memory ratio in Preferences.\n"), e.what());
		return MPP_ENOMEM;
	}
}

/* -------------------- Stage B: mpp_compute_shifts -------------------- */

static mpp_status_t mpp_compute_shifts_impl(sequence *seq, const mpp_config_t *cfg,
                                            mpp_run_t *run,
                                            const struct mpp_derot *derot) {
	if (!seq || !cfg || !run || !run->aps) return MPP_EINVAL;
	const bool derot_active = mpp_derot_applies(derot, run->num_frames,
	                                            run->frame_rows, run->frame_cols);

	/* Stage B re-reads analysis frames on cache misses — pin CFA SERs to
	 * raw-mosaic reads so those reads match the Stage A frames
	 * (cv::cvtColor demosaic) byte for byte. */
	SerAnalysisRawGuard raw_guard;
	maybe_engage_ser_raw_for_analysis(raw_guard, seq);

	/* Stage B streams: each frame is read + (optionally derotated) + blurred
	 * inside the correlation loop and discarded at iteration end. The frame
	 * loop is frame-parallel (each (frame, AP) shift is a disjoint slot), so
	 * the picker is asked for up to max_threads in-flight reads — the per-AP
	 * correlation plus the per-frame derotation remap then run across cores
	 * instead of serialising behind OpenCV's internal threading. Disk traffic
	 * is still 1 read per included frame. */
#ifdef _OPENMP
	const int stageb_max_threads = std::max(1, com.max_thread);
#else
	const int stageb_max_threads = 1;
#endif
	/* The provider is reentrant when its reads are: cache hits are trivially
	 * so, and disk-read misses are reentrant for SER / reentrant-FITS — same
	 * predicate Stage A uses for its parallel sweeps. */
	const bool stageb_provider_safe =
	    (seq->type == SEQ_SER
	     || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ
	          || seq->type == SEQ_INTERNAL) && fits_is_reentrant()));
	mpp_mem_pick mem_pick{MPP_MEM_STREAMING, 1, 0, 0};
	{
		const int bps = (cfg->bitdepth == 8) ? 1 : 2;
		const mpp_status_t mem = mpp_pick_memory_mode(
		    run->num_frames, (int) seq->rx, (int) seq->ry, /*channels=*/1,
		    bps,
		    /*cached_copies=*/0,
		    /*half_cached_copies=*/0,
		    stageb_provider_safe ? stageb_max_threads : 1,
		    mpp_analysis_in_flight_slots(seq, cfg),
		    /*output_bytes=*/0,
		    _("Register"), &mem_pick);
		if (mem != MPP_OK) return mem;
	}
	const int stageb_threads = std::max(1, mem_pick.threads);

	const auto apq = mpp::apq_from_run(run);
	const auto offsets = mpp::offsets_from_run(run);
	const cv::Mat mean_raw = mpp::mean_frame_from_run(run);

	gui_iface.set_progress(PROGRESS_RESET, _("Register: correlating per-AP boxes"));
	/* Provider hits the persistent analysis cache (populated by Stage A)
	 * for top-K-by-quality frames and falls through to disk on miss.
	 * For the typical planetary workflow where stack frame counts are
	 * substantial fractions of N and the user runs Stage B immediately
	 * after Stage A, this eliminates almost all disk reads here. */
	const mpp_cache_t *cache = run->cache;
	const int cache_N = cache ? (int) cache->raw_frames.size() : 0;
	if (cache && cache->populated_count > 0) {
		siril_log_message(_("Register: re-using Stage A analysis cache "
		                    "(%d/%d frames in memory)\n"),
		                  cache->populated_count, run->num_frames);
	}
	auto provider = [seq, cfg, cache, cache_N, derot, derot_active](int f) -> cv::Mat {
		return guarded_frame([&]() -> cv::Mat {
			cv::Mat mono;
			if (cache && f >= 0 && f < cache_N
			 && !cache->raw_frames[f].empty()) {
				mono = cache->raw_frames[f];
			} else {
				mono = mpp::read_analysis_frame(seq, f, cfg->avi_bayer_pattern,
				                                cfg->bitdepth);
			}
			if (mono.empty()) return cv::Mat();
			/* Match Stage A: shifts are measured against the (derotated) mean
			 * in epoch space, so each frame is derotated before correlation. */
			if (derot_active) mono = mpp_derot_warp_frame(derot, f, mono);
			return mpp::blur_mono_for_align(mono, *cfg);
		});
	};
	mpp_shifts_t *shifts = mpp::stack_compute_shifts_streamed(
	    provider, run->num_frames, mean_raw, *run->aps, apq, offsets,
	    *cfg, run->included, stageb_threads, stageb_provider_safe);
	if (!shifts) return MPP_ENOMEM;
	if (shifts->failure_counter == -1) {   /* cancellation sentinel */
		mpp_shift_free(shifts);
		return MPP_EINTR;
	}
	if (run->shifts) mpp_shift_free(run->shifts);
	run->shifts = shifts;
	return MPP_OK;
}

/* Exception boundary — see mpp_analyze. Stage B's correlation loop is
 * serial, so any cv:: allocation failure unwinds to here instead of
 * terminating inside an OpenMP region. */
extern "C" mpp_status_t mpp_compute_shifts(sequence *seq, const mpp_config_t *cfg,
                                           mpp_run_t *run,
                                           const struct mpp_derot *derot) {
	try {
		return mpp_compute_shifts_impl(seq, cfg, run, derot);
	} catch (const std::exception &e) {
		siril_log_error(_("Register: failed (%s) — likely out of memory; "
		                  "reduce the memory ratio in Preferences.\n"), e.what());
		return MPP_ENOMEM;
	}
}

/* -------------------- Stage C: mpp_stack_apply -------------------- */

extern "C" mpp_input_type mpp_classify_sequence_input(const sequence *seq) {
	if (!seq) return MPP_INPUT_MONO;
	if (seq->type == SEQ_SER && seq->ser_file) {
		const ser_color cid = seq->ser_file->color_id;
		if (cid >= SER_BAYER_RGGB && cid <= SER_BAYER_BGGR) return MPP_INPUT_CFA;
		if (cid == SER_RGB || cid == SER_BGR)               return MPP_INPUT_RGB;
		return MPP_INPUT_MONO;
	}
	/* FITS / FITSEQ: nb_layers separates already-RGB frames from
	 * single-layer ones, and a single-layer frame may still be a raw CFA
	 * mosaic (BAYERPAT in the header) that the mpp readers now debayer —
	 * fits_seq_is_cfa peeks frame 0 to tell mono from CFA. */
	if (seq->nb_layers > 1) return MPP_INPUT_RGB;
	if (mpp::fits_seq_is_cfa((sequence *) seq)) return MPP_INPUT_CFA;
	return MPP_INPUT_MONO;
}

static mpp_status_t mpp_stack_apply_impl(sequence *seq, const mpp_config_t *cfg,
                                         mpp_run_t *run,
                                         const struct mpp_derot *derot,
                                         fits *out) {
	if (!seq || !cfg || !run || !run->aps || !run->shifts || !out)
		return MPP_EINVAL;

	/* Refresh inclusion from the live frame selector: the user may have
	 * excluded more frames after registration (the sidecar / cached run
	 * carries the registration-time selection). AND the live incl state in
	 * so those frames are dropped from the stack. The per-AP weight
	 * normalisation rebuilds its denominator from whatever stays included,
	 * so dropped frames don't darken their patches, and an AP that loses
	 * all its frames falls back to the background composite. Dropped frames
	 * are NOT replaced by the next-best (Stage B computed no shifts for
	 * them) — a large post-registration selection change is best followed
	 * by re-running Register. */
	if (seq->imgparam && run->included && run->num_frames == seq->number) {
		int dropped = 0;
		for (int i = 0; i < run->num_frames; ++i) {
			const int live = seq->imgparam[i].incl ? 1 : 0;
			if (run->included[i] && !live) ++dropped;
			run->included[i] = run->included[i] && live;
		}
		if (dropped > 0)
			siril_log_message(_("Stack (mpp): %d frame(s) excluded in the "
			                    "frame selector since registration — dropping "
			                    "them from the stack\n"), dropped);
	}

	/* Honour Stack-tab per-AP top-K overrides. Stage A baked
	 * run->stack_size from the Analyze-time cfg; the GUI Stack tab and
	 * sidecar overrides land in cfg here without affecting run->stack_size
	 * automatically. Re-derive the target from cfg and clamp:
	 *   target ≤ baked: shrink in-place — best_frame_indices is sorted
	 *     descending by per-AP quality, so taking the first `target`
	 *     entries per AP is correct (apq_from_run keeps reading rows at
	 *     the baked selected_per_ap stride and re-derives the selection
	 *     weights from rank against the new stack_size).
	 *   target > baked: warn and keep baked. Analyze discarded the per-AP
	 *     quality matrix past the top selected_per_ap, so we can't grow
	 *     without re-running Stage A. */
	if (run->num_frames > 0) {
		const int target = mpp::stack_size_from_cfg(*cfg, run->num_frames);
		if (target > run->stack_size) {
			siril_log_warning(
			    _("Stack (mpp): requested %d top-quality frames per AP but "
			      "registration baked %d — keeping %d. Re-register with a "
			      "higher Register percent to raise the ceiling.\n"),
			    target, run->stack_size, run->stack_size);
		} else if (target < run->stack_size) {
			run->stack_size = target;
		}
		/* Always announce the effective per-AP stack size so the user
		 * can verify Stack-tab/Register-tab settings are taking effect
		 * (silence when stack_size == num_frames means no per-AP filter
		 * is active and the message would be noise). */
		if (run->stack_size < run->num_frames) {
			if (run->taper > 0)
				siril_log_message(
				    _("Stack (mpp): per-AP frame selection active — "
				      "%d/%d frames per AP, soft taper ±%d.\n"),
				    run->stack_size, run->num_frames, run->taper);
			else
				siril_log_message(
				    _("Stack (mpp): per-AP top-K filter active — %d/%d "
				      "frames per AP.\n"),
				    run->stack_size, run->num_frames);
		}
	}

	/* Stage C can't consume the analysis cache (single-channel mono at
	 * analysis bitdepth vs. read_full_frame's multi-channel CV_16U), so
	 * release it before the streaming-only memory picker fires.
	 * Otherwise the cache sits resident through Stage C, silently eating
	 * budget the picker doesn't account for and risking OOM. The GUI's
	 * com.mpp_run typically holds a separate copy of the cache (sidecar
	 * loads strip cache → `run` here is cache-less in that flow, but
	 * com.mpp_run isn't); drop it too. AP re-edits after stacking
	 * repopulate the cache via mpp_recompute_qualities at one extra
	 * disk pass. */
	mpp_run_drop_cache(run);
	mpp_run_t *cached = mpp_get_cached_run();
	if (cached && cached != run) mpp_run_drop_cache(cached);

	/* Dispatch: all stacking goes through the classical engine
	 * (stack_apply_shifts_streamed) — sub-pixel rigid-shift per-AP patch
	 * accumulation + Hann-weight mosaic, parallelised over frames, with
	 * cv::resize interpolation for scale > 1.
	 *
	 * The STScI / Bayer dobox drizzle backends were removed: both produced
	 * sampling / mosaic-phase artefacts on planetary data. (The code is
	 * preserved on the pss-drizzle branch in case a future CFA-aware
	 * revival is wanted.)
	 *
	 * CFA input is debayered first (read_full_frame forces ser_read_frame
	 * to debayer regardless of the user's debayer-on-open pref), then
	 * stacked classically like any RGB sequence. */
	/* AVI Bayer routing: an AVI marked as a raw mosaic (cfg->avi_bayer_pattern
	 * RGGB/BGGR/GBRG/GRBG) is the AVI analogue of MPP_INPUT_CFA. The
	 * classifier can't see this (sequence-only signature) so it's
	 * tracked here for the channel-count math below. */
	const bool avi_cfa = (seq->type == SEQ_AVI
	                      && cfg->avi_bayer_pattern >= MPP_AVI_BAYER_RGGB
	                      && cfg->avi_bayer_pattern <= MPP_AVI_BAYER_GRBG);
	const mpp_input_type input_type = mpp_classify_sequence_input(seq);
	const bool is_cfa = (input_type == MPP_INPUT_CFA) || avi_cfa;

	/* The classical (OpenCV-interpolation) path handles every input and
	 * scale; a sidecar/script that explicitly pinned a dobox backend is
	 * overridden — drizzle is no longer applied for scaling. */
	if (cfg->drizzle_mode != MPP_DRIZZLE_OFF)
		siril_log_message(_("Stack (mpp): dobox drizzle backend disabled — "
		                    "scaling uses OpenCV interpolation\n"));

	/* Pre-flight memory check for the classical path. */
	{
		/* CFA inputs (SER mosaic or AVI Bayer) decode to 3-channel RGB via
		 * forced RCD debayer / cv::cvtColor; run->num_layers may still read 1
		 * from the Stage-A capture, so force 3 for the estimate. */
		const int channels = is_cfa
		                   ? 3
		                   : (run->num_layers > 0 ? run->num_layers : 1);
		mpp_mem_pick mem_pick{MPP_MEM_STREAMING, 1, 0, 0};
		/* The classical path decodes a batch of ~2 frames per thread
		 * concurrently (CV_32F, hence bytes_per_sample 4) before
		 * accumulating, so its working set scales with the thread budget. */
#ifdef _OPENMP
		const int mem_threads = std::max(1, com.max_thread);
#else
		const int mem_threads = 1;
#endif
		const mpp_status_t mem = mpp_pick_memory_mode(
		    run->num_frames, run->frame_cols, run->frame_rows, channels,
		    /*bytes_per_sample=*/4,
		    /*cached_copies=*/0,
		    /*half_cached_copies=*/0,
		    /*max_threads=*/mem_threads,
		    /*per_thread_in_flight=*/2,
		    /*output_bytes=*/0,
		    _("Stack (mpp)"), &mem_pick);
		if (mem != MPP_OK) return mem;
		(void) mem_pick;
	}

	const auto apq = mpp::apq_from_run(run);
	const auto offsets = mpp::offsets_from_run(run);

	/* Quality-sorted index for top-N-for-background pick — only included
	 * frames participate. Otherwise excluded frames would contribute to
	 * the background blend even though they're skipped from the AP buffers. */
	std::vector<int> sorted_idx;
	sorted_idx.reserve(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i)
		if (run->included[i]) sorted_idx.push_back(i);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [run](int a, int b) { return run->quality[a] > run->quality[b]; });

	std::vector<double> frame_brightness(run->frame_brightness,
	                                     run->frame_brightness + run->num_frames);
	const cv::Vec4i intersection(run->intersection[0], run->intersection[1],
	                             run->intersection[2], run->intersection[3]);

	gui_iface.set_progress(PROGRESS_RESET, _("Stack: accumulating frames"));
	auto provider = [seq, cfg](int i) -> cv::Mat {
		return guarded_frame([&]() -> cv::Mat {
			return mpp::read_full_frame(seq, i, cfg->avi_bayer_pattern);
		});
	};
	/* read_full_frame returns 3-channel RGB for any CFA input (AVI Bayer
	 * mosaics via cv::cvtColor, CFA SERs via forced RCD debayer), even
	 * though run->num_layers (captured from seq->nb_layers at Stage A time)
	 * may read 1. The accumulator type depends on num_layers, so force 3
	 * for CFA. */
	int num_layers = run->num_layers > 0 ? run->num_layers : 1;
	if (is_cfa) num_layers = 3;
#ifdef _OPENMP
	const int stack_threads = std::max(1, com.max_thread);
#else
	const int stack_threads = 1;
#endif
	/* Parallel decode needs a reentrant provider — same predicate as Stage
	 * A's parallel frame read: SER serialises only the raw fread, and FITS
	 * is safe only when cfitsio was built reentrant. */
	const bool stack_provider_safe =
	    (seq->type == SEQ_SER
	     || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ
	          || seq->type == SEQ_INTERNAL) && fits_is_reentrant()));
	siril_log_message(_("Stack (mpp): classical accumulation, %d thread(s), "
	                    "%s frame reads\n"),
	                  stack_threads,
	                  stack_provider_safe ? _("parallel") : _("serial"));
	/* Memory budget for the per-thread private accumulators: 80% of the
	 * lesser of free RAM and the configured cap. stack_apply_shifts_streamed
	 * clamps its thread count to what fits. */
	size_t stack_mem_budget = 0;
	{
		guint64 budget_b = get_available_memory();
		const long long cap_mb = get_max_memory_in_MB();
		if (cap_mb > 0) {
			const guint64 cap_b = (guint64) cap_mb * 1024ULL * 1024ULL;
			if (cap_b < budget_b) budget_b = cap_b;
		}
		stack_mem_budget = (size_t) ((long double) budget_b * 0.8L);
	}
	/* Stage C engine: the warp-field engine (forced when derotation is
	 * active, see stack_mpp.c) warps each frame once through a dense
	 * displacement field; otherwise the PSS per-AP patch-blend engine. */
	cv::Mat stacked;
	if (cfg->stack_method == MPP_STACK_WARP) {
		/* Per-frame derotation maps in the drizzled intersection canvas: pixel
		 * (0,0) is original-frame coord (intersection x_low, y_low), scaled by
		 * the drizzle factor. Only built when a valid .derot plan is present. */
		const double S = std::max(1.0, cfg->drizzle_scale);
		mpp::DerotMapProvider derot_provider =
		    [derot, &intersection, S](int f, int DY, int DX,
		                              cv::Mat &mx, cv::Mat &my, cv::Mat &mu) {
			mpp_derot_frame_map(derot, f, DX, DY,
			                    (double) intersection[2], (double) intersection[0],
			                    S, mx, my, mu);
			return true;
		};
		const bool use_derot = mpp_derot_applies(derot, run->num_frames,
		                                         run->frame_rows, run->frame_cols);
		if (derot && !use_derot)
			siril_log_warning(_("Stack (mpp): the .derot plan does not match this "
			                    "sequence (frame count or dimensions); derotation "
			                    "skipped.\n"));
		const mpp::DerotMapProvider *derot_pp = use_derot ? &derot_provider : nullptr;
		const mpp::WarpStackResult wr = mpp::stack_warp_apply_streamed(
		    provider, run->num_frames, num_layers, *run->aps, apq,
		    run->shifts, offsets, frame_brightness, sorted_idx,
		    intersection, *cfg, run->included, stack_threads, stack_provider_safe,
		    stack_mem_budget, derot_pp);
		if (wr.oom) {
			siril_log_error(_("Stack (mpp): out of memory while accumulating "
			                  "frames — reduce the memory ratio in Preferences.\n"));
			return MPP_ENOMEM;
		}
		if (wr.cancelled || wr.image.empty()) return MPP_EINTR;
		stacked = wr.image;
	} else {
		const auto loop = mpp::stack_apply_shifts_streamed(
		    provider, run->num_frames, num_layers, *run->aps, apq,
		    run->shifts, offsets, frame_brightness, sorted_idx,
		    intersection, *cfg, run->included, stack_threads, stack_provider_safe,
		    stack_mem_budget);
		if (loop.oom) {
			siril_log_error(_("Stack (mpp): out of memory while accumulating "
			                  "frames — reduce the memory ratio in Preferences.\n"));
			return MPP_ENOMEM;
		}
		if (loop.state.dim_y == 0) return MPP_EINTR;   /* cancellation sentinel */
		stacked = mpp::stack_merge_alignment_point_buffers(
		    loop.state, loop.border, *run->aps, *cfg);
	}

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
	/* Siril's mono convention: pdata[1]/pdata[2] alias pdata[0] (not
	 * NULL). remap_all_vports's per-vport memcpy loops c=0..2
	 * unconditionally, so a NULL plane pointer crashes the GUI
	 * refresh that runs immediately after stack_function_handler
	 * swaps gfit. */
	out->pdata[1] = (C >= 2) ? out->data + plane     : out->data;
	out->pdata[2] = (C >= 3) ? out->data + plane * 2 : out->data;

	/* Optional 32-bit float output. The merge above runs in 16-bit; when
	 * requested (`-32b` / "Force 32b output") convert the packed USHORT
	 * buffer to normalised float in place. fit_replace_buffer frees the
	 * old WORD buffer and rewires pdata/fpdata. */
	if (cfg->output_32bit) {
		float *fbuf = ushort_buffer_to_float(out->data, plane * C);
		if (!fbuf) return MPP_ENOMEM;
		fit_replace_buffer(out, fbuf, DATA_FLOAT);
		siril_log_message(_("Stacking result will be stored as a 32-bit image\n"));
	}
	return MPP_OK;
}

/* Exception boundary — see mpp_analyze. */
extern "C" mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                                        mpp_run_t *run,
                                        const struct mpp_derot *derot,
                                        fits *out) {
	try {
		return mpp_stack_apply_impl(seq, cfg, run, derot, out);
	} catch (const std::exception &e) {
		siril_log_error(_("Stack (mpp): failed (%s) — likely out of memory; "
		                  "reduce the memory ratio in Preferences.\n"), e.what());
		return MPP_ENOMEM;
	}
}

/* -------------------- GUI cache for the current run --------------------
 *
 * The cached run is stored on com.mpp_run (declared as void* in siril.h to
 * avoid cycles). Lifetime: installed by Stage A on success, freed by
 * close_sequence (or replaced by a subsequent Stage A). The AP overlay in
 * gui/image_display.c reads it under the mutex; the AP editor mutates it
 * via mpp_get_cached_run + the mpp_ap_* helpers.
 */
static GMutex mpp_cache_mutex;

extern "C" void mpp_set_cached_run(mpp_run_t *run) {
	g_mutex_lock(&mpp_cache_mutex);
	mpp_run_t *prev = (mpp_run_t *) com.mpp_run;
	com.mpp_run = run;
	g_mutex_unlock(&mpp_cache_mutex);
	if (prev) mpp_run_free(prev);
}

extern "C" void mpp_clear_cached_run(void) {
	g_mutex_lock(&mpp_cache_mutex);
	mpp_run_t *prev = (mpp_run_t *) com.mpp_run;
	com.mpp_run = NULL;
	g_mutex_unlock(&mpp_cache_mutex);
	if (prev) mpp_run_free(prev);
}

static mpp_status_t mpp_recompute_qualities_impl(sequence *seq, mpp_run_t *run,
                                                 const struct mpp_derot *derot) {
	if (!seq || !run || !run->aps || run->aps->count <= 0 || !run->cfg)
		return MPP_EINVAL;
	if (run->best_frame_indices) return MPP_OK;   /* nothing to do */
	if (!run->frame_brightness || !run->global_shifts) return MPP_EINVAL;
	if (run->num_frames != seq->number) return MPP_EINVAL;
	const bool derot_active = mpp_derot_applies(derot, run->num_frames,
	                                            run->frame_rows, run->frame_cols);

	/* Analysis-frame reads — pin CFA SERs to raw-mosaic reads so they
	 * match Stage A's cv::cvtColor demosaic. */
	SerAnalysisRawGuard raw_guard;
	maybe_engage_ser_raw_for_analysis(raw_guard, seq);

	const int N = run->num_frames;
#ifdef _OPENMP
	const int max_threads = std::max(1, com.max_thread);
#else
	const int max_threads = 1;
#endif
	const int bps = (run->cfg->bitdepth == 8) ? 1 : 2;

	/* If Stage A left an analysis cache on the run, skip the pre-read
	 * entirely and feed the streamed pass directly from it (hits for
	 * cached slots, disk reads for the rest). Without a cache (sidecar-
	 * loaded run, or Stage A picked STREAMING), do the pre-read of the
	 * largest top-by-index slice that fits. */
	mpp_cache_t *fresh_cache_owner = nullptr;
	std::vector<cv::Mat> fresh_cache_vec;
	int fresh_cache_count = 0;
	/* When a pre-read pass runs it owns the first half of this stage's bar
	 * (0..0.5) and the quality compute owns the second (0.5..1.0). With the
	 * frames already cached there is no read pass, so the quality compute
	 * sweeps the whole bar. */
	const bool did_read_pass = (run->cache == nullptr);
	if (!run->cache) {
		mpp_mem_pick mem_pick{MPP_MEM_CACHED, max_threads, N, 0};
		/* cached_copies=1 and half_cached_copies=1 are identical here
		 * because there is no blur cache to drop — the distinction
		 * exists only so the picker also considers PARTIAL_CACHED. */
		const mpp_status_t mem = mpp_pick_memory_mode(
		    N, (int) seq->rx, (int) seq->ry, /*channels=*/1,
		    bps,
		    /*cached_copies=*/1,
		    /*half_cached_copies=*/1,
		    max_threads,
		    mpp_analysis_in_flight_slots(seq, run->cfg),
		    /*output_bytes=*/0,
		    _("Register (per-AP quality refresh)"), &mem_pick);
		if (mem != MPP_OK) return mem;
		fresh_cache_count = mem_pick.cached_count;
		const int preread_threads = mem_pick.threads;

		siril_log_message(_("mpp: computing per-AP frame qualities "
		                    "(reading %d frames)\n"), N);
		gui_iface.set_progress(PROGRESS_RESET, _("Register: reading frames for per-AP ranking"));

		if (fresh_cache_count > 0) {
			fresh_cache_vec.resize(N);
			int cur_nb = 0;
			int read_failed = 0;
			int alloc_failed = 0;
			int cancelled = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(preread_threads) schedule(guided) \
    if (seq->type == SEQ_SER || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ \
                                  || seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
			for (int i = 0; i < fresh_cache_count; ++i) {
				if (g_atomic_int_get(&read_failed) || g_atomic_int_get(&cancelled)
				    || g_atomic_int_get(&alloc_failed)) continue;
				if (!processing_should_continue()) {
					g_atomic_int_set(&cancelled, 1);
					continue;
				}
				/* Same OOM containment as Stage A Pass 1 — an exception
				 * escaping the OpenMP region is std::terminate(). */
				try {
					cv::Mat mono = mpp::read_analysis_frame(seq, i, run->cfg->avi_bayer_pattern,
					                                        run->cfg->bitdepth);
					if (mono.empty()) {
						g_atomic_int_set(&read_failed, 1);
						continue;
					}
					fresh_cache_vec[i] = mono;
				} catch (const std::exception &) {
					g_atomic_int_set(&alloc_failed, 1);
					continue;
				}
				g_atomic_int_inc(&cur_nb);
				gui_iface.set_progress(0.5 * (double) cur_nb / (double) fresh_cache_count, NULL);
			}
			if (cancelled)  return MPP_EINTR;
			if (alloc_failed) {
				siril_log_error(_("Register: out of memory while re-reading "
				                  "frames — reduce the memory ratio in "
				                  "Preferences.\n"));
				return MPP_ENOMEM;
			}
			if (read_failed) return MPP_EIO;
		}
		/* Promote the freshly populated cache onto the run so future
		 * recomputes (and Stage B / Stage C) also hit it. */
		if (fresh_cache_count > 0) {
			fresh_cache_owner = new mpp_cache;
			fresh_cache_owner->raw_frames     = std::move(fresh_cache_vec);
			fresh_cache_owner->bitdepth       = run->cfg->bitdepth;
			fresh_cache_owner->populated_count = fresh_cache_count;
			run->cache = fresh_cache_owner;
		}
	} else {
		siril_log_message(_("mpp: computing per-AP frame qualities "
		                    "(%d/%d frames already cached)\n"),
		                  run->cache->populated_count, N);
		gui_iface.set_progress(PROGRESS_RESET, _("Register: re-using cached analysis frames"));
	}

	const double qual_base = did_read_pass ? 0.5 : 0.0;
	gui_iface.set_progress(qual_base, _("Register: computing per-AP qualities"));
	const std::vector<double> frame_brightness(run->frame_brightness,
	                                           run->frame_brightness + N);
	const auto offsets = mpp::offsets_from_run(run);
	stage_progress_range rp{qual_base, 1.0 - qual_base};
	const mpp_cache_t *cache = run->cache;
	const int cache_N = cache ? (int) cache->raw_frames.size() : 0;
	auto raw_provider = [cache, cache_N, seq, run, derot, derot_active](int i) -> cv::Mat {
		return guarded_frame([&]() -> cv::Mat {
			cv::Mat m;
			if (cache && i >= 0 && i < cache_N
			 && !cache->raw_frames[i].empty())
				m = cache->raw_frames[i];
			else
				m = mpp::read_analysis_frame(seq, i, run->cfg->avi_bayer_pattern,
				                             run->cfg->bitdepth);
			/* Per-AP ranking samples the AP boxes, so it must see the same
			 * (derotated) frames as the mean and the shift stage. */
			return derot_active ? mpp_derot_warp_frame(derot, i, m) : m;
		});
	};
	/* Rank only included frames — excluded frames must not occupy any AP's
	 * top-N slot (they'd be dropped at Stage C anyway, wasting the slot). */
	const std::vector<int> included(run->included, run->included + N);
	/* Frame-parallel quality compute when the provider is reentrant — same
	 * predicate as Stage A / Stage B (cache hits trivially, SER / reentrant-
	 * FITS on disk-read misses). Heavy once derotation adds a per-frame
	 * remap, so this is what keeps the cores busy here. */
	const bool quality_provider_safe =
	    (seq->type == SEQ_SER
	     || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ
	          || seq->type == SEQ_INTERNAL) && fits_is_reentrant()));
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, frame_brightness, *run->aps, offsets,
	    run->frame_rows, run->frame_cols, *run->cfg,
	    stage_progress_cb, &rp, included,
	    quality_provider_safe ? max_threads : 1, quality_provider_safe);
	if (apq.stack_size <= 0) return MPP_EINTR;

	const int stack_size = apq.stack_size > 0 ? apq.stack_size : run->stack_size;
	const int selected = run->aps->count > 0 && !apq.best_frame_indices.empty()
	    ? (int) apq.best_frame_indices[0].size() : stack_size;
	const size_t total = (size_t) run->aps->count * (size_t) selected;
	int *bfi = (int *) std::malloc(total * sizeof(int));
	if (!bfi) return MPP_ENOMEM;
	for (int a = 0; a < run->aps->count; ++a) {
		for (int k = 0; k < selected; ++k) {
			bfi[a * selected + k] =
			    (k < (int) apq.best_frame_indices[a].size())
			        ? apq.best_frame_indices[a][k]
			        : -1;
		}
	}
	std::free(run->best_frame_indices);
	run->best_frame_indices = bfi;
	run->stack_size         = stack_size;
	run->taper              = apq.taper;
	run->selected_per_ap    = selected;
	siril_log_message(_("mpp: per-AP qualities computed for %d APs\n"),
	                  run->aps->count);
	return MPP_OK;
}

/* Exception boundary — see mpp_analyze. */
extern "C" mpp_status_t mpp_recompute_qualities(sequence *seq, mpp_run_t *run,
                                                const struct mpp_derot *derot) {
	try {
		return mpp_recompute_qualities_impl(seq, run, derot);
	} catch (const std::exception &e) {
		siril_log_error(_("Register: per-AP quality refresh failed (%s) — "
		                  "likely out of memory; reduce the memory ratio in "
		                  "Preferences.\n"), e.what());
		return MPP_ENOMEM;
	}
}

extern "C" mpp_status_t mpp_ap_replace(mpp_run_t *run, const mpp_config_t *cfg) {
	if (!run || !cfg || !run->mean_frame_data
	    || run->mean_frame_rows <= 0 || run->mean_frame_cols <= 0)
		return MPP_EINVAL;
	cv::Mat mean(run->mean_frame_rows, run->mean_frame_cols, CV_32S,
	             run->mean_frame_data);
	cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(mean, *cfg);
	mpp_aps_t *new_aps = mpp::ap_create_grid(mean_blurred, *cfg);
	if (!new_aps) return MPP_ENOMEM;
	if (new_aps->count <= 0) {
		mpp_ap_free(new_aps);
		return MPP_ENODATA;
	}
	mpp_ap_free(run->aps);
	run->aps = new_aps;
	if (run->cfg) *run->cfg = *cfg;
	/* Per-AP qualities are now stale (different AP set). Clear so
	 * register_mpp can detect at Stage B time and refuse / recompute. */
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
	return MPP_OK;
}

extern "C" void mpp_ap_clear_all(mpp_run_t *run) {
	if (!run) return;
	if (run->aps) {
		std::free(run->aps->records);
		run->aps->records = nullptr;
		run->aps->count = 0;
	}
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
}

/* Fill the box/patch bounds for an AP at (x, y) for an explicit half-box
 * (hb) and half-patch (hp), clamped to the given frame dimensions. */
static void fill_ap_bounds_hb(mpp_ap_record_t *r, int x, int y,
                              int hb, int hp,
                              int frame_rows, int frame_cols) {
	r->x = x;
	r->y = y;
	r->box_x_low  = std::max(0, x - hb);
	r->box_x_high = std::min(frame_cols, x + hb);
	r->box_y_low  = std::max(0, y - hb);
	r->box_y_high = std::min(frame_rows, y + hb);
	r->patch_x_low  = std::max(0, x - hp);
	r->patch_x_high = std::min(frame_cols, x + hp);
	r->patch_y_low  = std::max(0, y - hp);
	r->patch_y_high = std::min(frame_rows, y + hp);
}

/* Fill the box/patch bounds for an AP at (x, y) from cfg, extending the
 * patch out to the frame border for edge APs. */
static void fill_ap_bounds(mpp_ap_record_t *r, int x, int y,
                           const mpp_config_t *cfg,
                           int frame_rows, int frame_cols) {
	fill_ap_bounds_hb(r, x, y, cfg->alignment_points_half_box_width,
	                  mpp_cfg_half_patch_width(cfg), frame_rows, frame_cols);
}

extern "C" mpp_status_t mpp_ap_add(mpp_run_t *run, int x, int y) {
	if (!run || !run->cfg) return MPP_EINVAL;
	if (!run->aps) {
		run->aps = (mpp_aps_t *) std::calloc(1, sizeof(*run->aps));
		if (!run->aps) return MPP_ENOMEM;
	}
	const int new_count = run->aps->count + 1;
	mpp_ap_record_t *grown = (mpp_ap_record_t *) std::realloc(
		run->aps->records, (size_t) new_count * sizeof(mpp_ap_record_t));
	if (!grown) return MPP_ENOMEM;
	run->aps->records = grown;
	mpp_ap_record_t *r = &grown[run->aps->count];
	std::memset(r, 0, sizeof(*r));
	fill_ap_bounds(r, x, y, run->cfg, run->frame_rows, run->frame_cols);
	r->structure = 1.0;   /* user-added APs bypass the structure threshold */
	run->aps->count = new_count;
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
	return MPP_OK;
}

extern "C" mpp_status_t mpp_ap_remove(mpp_run_t *run, int i) {
	if (!run || !run->aps || i < 0 || i >= run->aps->count) return MPP_EINVAL;
	const int last = run->aps->count - 1;
	if (i < last) {
		std::memmove(&run->aps->records[i], &run->aps->records[i + 1],
		             (size_t)(last - i) * sizeof(mpp_ap_record_t));
	}
	run->aps->count = last;
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
	return MPP_OK;
}

extern "C" mpp_status_t mpp_ap_move(mpp_run_t *run, int i, int x, int y) {
	if (!run || !run->aps || !run->cfg || i < 0 || i >= run->aps->count)
		return MPP_EINVAL;
	mpp_ap_record_t *r = &run->aps->records[i];
	fill_ap_bounds(r, x, y, run->cfg, run->frame_rows, run->frame_cols);
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
	return MPP_OK;
}

/* Current half-box of an AP, recovered as the largest centre-to-edge
 * distance — the box is clamped at the frame border for edge APs, so the
 * unclamped side gives the true value. Boxes produced by mpp_ap_set_half_box
 * are symmetric (the new half-box is capped below the border), so this is
 * exact for any AP the editor has touched. */
static int ap_current_half_box(const mpp_ap_record_t *r) {
	return std::max(std::max(r->x - r->box_x_low, r->box_x_high - r->x),
	                std::max(r->y - r->box_y_low, r->box_y_high - r->y));
}

extern "C" int mpp_ap_get_half_box(const mpp_run_t *run, int i) {
	if (!run || !run->aps || i < 0 || i >= run->aps->count) return -1;
	return ap_current_half_box(&run->aps->records[i]);
}

/* Set AP `i`'s half-box to `hb` (and the patch to 1.5× hb, matching
 * mpp_cfg_half_patch_width), recomputing both around the AP centre. `hb` is
 * clamped to [6, border − search_width] so the symmetric box stays inside the
 * mean (reference) frame — which is what the boxes index — with room for the
 * per-AP search. GUI-only: the editor edits the cached run; resetting
 * best_frame_indices makes Register recompute the per-AP qualities. */
extern "C" mpp_status_t mpp_ap_set_half_box(mpp_run_t *run, int i, int hb) {
	if (!run || !run->aps || !run->cfg || i < 0 || i >= run->aps->count)
		return MPP_EINVAL;
	mpp_ap_record_t *r = &run->aps->records[i];
	const int x = r->x, y = r->y;
	const int rows = run->mean_frame_rows, cols = run->mean_frame_cols;
	if (rows <= 0 || cols <= 0) return MPP_EINVAL;
	const int sw = run->cfg->alignment_points_search_width;
	const int min_hb = 6;
	int max_hb = std::min(std::min(x, cols - 1 - x),
	                      std::min(y, rows - 1 - y)) - sw;
	if (max_hb < min_hb) max_hb = min_hb;
	int new_hb = hb;
	if (new_hb < min_hb) new_hb = min_hb;
	if (new_hb > max_hb) new_hb = max_hb;
	if (new_hb == ap_current_half_box(r)) return MPP_OK;   /* no change */
	fill_ap_bounds_hb(r, x, y, new_hb, (new_hb * 3) / 2, rows, cols);
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
	return MPP_OK;
}

/* Grow/shrink AP `i`'s box (and patch) by `step` px of half-box. */
extern "C" mpp_status_t mpp_ap_resize(mpp_run_t *run, int i, int step) {
	if (!run || !run->aps || i < 0 || i >= run->aps->count) return MPP_EINVAL;
	return mpp_ap_set_half_box(run, i, ap_current_half_box(&run->aps->records[i]) + step);
}

extern "C" int mpp_ap_hit_test(const mpp_run_t *run, int x, int y) {
	if (!run || !run->aps) return -1;
	/* Iterate from newest (highest index) backward so user-added APs win
	 * over older auto-placed ones when boxes overlap. */
	for (int i = run->aps->count - 1; i >= 0; --i) {
		const mpp_ap_record_t *r = &run->aps->records[i];
		if (x >= r->box_x_low && x < r->box_x_high &&
		    y >= r->box_y_low && y < r->box_y_high)
			return i;
	}
	return -1;
}

extern "C" void mpp_display_to_ap_coord(const mpp_run_t *run,
                                        int gfit_rx, int gfit_ry,
                                        int current_frame_idx,
                                        int click_x, int click_y,
                                        int *out_x, int *out_y) {
	if (!out_x || !out_y) return;
	if (!run || gfit_ry <= 0) {
		*out_x = click_x;
		*out_y = click_y;
		return;
	}
	int dx = 0, dy = 0;
	const bool showing_ref = (gfit_rx == run->mean_frame_cols &&
	                          gfit_ry == run->mean_frame_rows);
	if (!showing_ref && run->global_shifts &&
	    current_frame_idx >= 0 && current_frame_idx < run->num_frames) {
		/* Click → integer pixel coord. The sub-pixel residual in
		 * global_shifts is irrelevant for hit-testing AP boxes. */
		dy = (int) std::lround((double) run->intersection[0]
		                       - run->global_shifts[2 * current_frame_idx + 0]);
		dx = (int) std::lround((double) run->intersection[2]
		                       - run->global_shifts[2 * current_frame_idx + 1]);
	}
	/* Inverse of the draw-time mapping:
	 *   display_x = ap.x + dx
	 *   display_y = (gfit_ry - 1) - (ap.y + dy)
	 */
	*out_x = click_x - dx;
	*out_y = (gfit_ry - 1) - click_y - dy;
}

extern "C" mpp_aps_t *mpp_aps_snapshot(const mpp_aps_t *src) {
	if (!src) return nullptr;
	mpp_aps_t *out = (mpp_aps_t *) std::calloc(1, sizeof(*out));
	if (!out) return nullptr;
	out->count            = src->count;
	out->dropped_dim      = src->dropped_dim;
	out->dropped_structure = src->dropped_structure;
	if (src->count > 0 && src->records) {
		const size_t bytes = (size_t) src->count * sizeof(mpp_ap_record_t);
		out->records = (mpp_ap_record_t *) std::malloc(bytes);
		if (!out->records) { std::free(out); return nullptr; }
		std::memcpy(out->records, src->records, bytes);
	}
	return out;
}

extern "C" void mpp_aps_restore(mpp_run_t *run, mpp_aps_t *snapshot) {
	if (!run || !snapshot) return;
	mpp_ap_free(run->aps);
	run->aps = snapshot;
	std::free(run->best_frame_indices);
	run->best_frame_indices = nullptr;
}

extern "C" mpp_run_t *mpp_get_cached_run(void) {
	g_mutex_lock(&mpp_cache_mutex);
	mpp_run_t *r = (mpp_run_t *) com.mpp_run;
	g_mutex_unlock(&mpp_cache_mutex);
	return r;
}

/* -------------------- Siril register-framework entry point -------------------- */

/* Populate seq->regparam[layer][i].quality from a completed Stage-A run so
 * Siril's existing frame selector / quality plot displays the per-frame
 * Laplace-σ. Public via mpp.h so the CLI and the GUI register entry can
 * all wire it in. The layer choice is caller-supplied; mpp quality is
 * panchromatic (same score on every channel), but we only populate one
 * layer to mirror Siril's per-layer regdata model — the existing frame
 * selector reads one layer at a time. */
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

/* True if two cfgs differ in any field that feeds AP placement
 * (ap_create_grid / blur_mean_frame_for_ap). When Register reuses a cached
 * or sidecar-loaded run, a difference here means the baked AP grid no longer
 * matches the requested settings and must be regenerated. Fields with no GUI
 * widget (frames_gauss_width, dim_fraction_threshold) always carry their
 * default in both cfgs, so they never spuriously trigger a regeneration. */
static bool ap_placement_cfg_differs(const mpp_config_t *a, const mpp_config_t *b) {
	return a->alignment_points_half_box_width        != b->alignment_points_half_box_width
	    || a->alignment_points_search_width          != b->alignment_points_search_width
	    || a->alignment_points_brightness_threshold  != b->alignment_points_brightness_threshold
	    || a->alignment_points_contrast_threshold    != b->alignment_points_contrast_threshold
	    || a->alignment_points_structure_threshold   != b->alignment_points_structure_threshold
	    || a->alignment_points_dim_fraction_threshold != b->alignment_points_dim_fraction_threshold
	    || a->frames_gauss_width                     != b->frames_gauss_width;
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

	const gboolean stage_a_only = regargs->mpp_stage_a_only;
	siril_log_status(stage_a_only
	                        ? _("mpp: MODE = ANALYZE (Stage A only, %d frames)\n")
	                        : _("mpp: MODE = REGISTER (Stages A + B, %d frames)\n"), regargs->seq->number);
	/* Belt-and-braces: prove the flag value to the script log so a regression
	 * here is immediately visible. mpp_stage_a_only is the only thing that
	 * gates the early return after Stage A. */
	siril_log_debug("register_mpp: regargs=%p mpp_stage_a_only=%d filter_included=%d\n",
	                  (void *) regargs, (int) regargs->mpp_stage_a_only,
	                  (int) regargs->filters.filter_included);
	/* Reuse-cached path: if Stage A already ran (Analyze button) and the
	 * cached run matches this sequence, skip Stage A and go straight to
	 * Stage B. This is what makes the AP editor's edits actually take
	 * effect — Register uses the (possibly-edited) cached AP grid rather
	 * than rebuilding it from scratch. The cache is cleared on
	 * close_sequence so there's no risk of using a stale run from a
	 * different sequence. */
	mpp_run_t *cached = mpp_get_cached_run();
	const bool reuse_cache = (cached && cached->aps && cached->aps->count > 0
	                          && cached->num_frames == regargs->seq->number);
	mpp_run_t *run = nullptr;
	bool run_owned = true;   /* if false, run is borrowed from com.mpp_run; do not free */

	if (reuse_cache) {
		run = cached;
		run_owned = false;
		siril_log_status(stage_a_only
		                        ? _("mpp: re-Analyze — reusing cached Stage A run; "
		                            "AP grid will be replaced.\n")
		                        : _("mpp: reusing cached Stage A run from previous Analyze "
		                            "(%d APs).\n"),
		                        cached->aps->count);
		if (stage_a_only) {
			/* Re-Analyze with cached run: drop the cache here so the fresh
			 * mpp_analyze below builds a new run that supersedes it. */
			run = nullptr;
			run_owned = true;
		}
	}

	if (!run) {
		gui_iface.set_progress(PROGRESS_RESET,
		                       stage_a_only ? _("Analyze: ranking frames")
		                                    : _("Register: ranking frames"));
		mpp_derot_t *derot_a = NULL;
		{
			gchar *dp = g_strdup_printf("%s.derot", regargs->seq->seqname);
			mpp_derot_read(dp, &derot_a);
			g_free(dp);
		}
		const int rc_a = mpp_analyze(regargs->seq, &cfg, &run, derot_a);
		mpp_derot_free(derot_a);
		if (rc_a == MPP_EINTR) {
			siril_log_message(_("mpp: Stage A cancelled by user.\n"));
			gui_iface.set_progress(PROGRESS_DONE, _("Cancelled"));
			return rc_a;
		}
		if (rc_a != MPP_OK) {
			siril_log_error(_("mpp: Stage A failed (code %d)\n"), rc_a);
			gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
			return rc_a;
		}
	}
	siril_log_message(_("mpp: Stage A done — best=%d, %d APs, stack_size=%d\n"),
	                  run->best_frame_idx, run->aps->count, run->stack_size);

	/* Phase 9.2: surface per-frame quality through Siril's regdata so the
	 * existing frame selector and quality plot pick it up. */
	mpp_write_quality_to_regdata(regargs->seq, regargs->layer, run);

	if (stage_a_only) {
		/* Persist Stage A to the .mpp sidecar (mean reference frame + AP
		 * grid; per-AP shifts omitted, so it reads back as "analysed,
		 * registration pending"). This means closing and reopening the
		 * sequence restores the APs and reference frame via the auto-load
		 * in load_sequence — no need to re-Analyze. Non-fatal on failure:
		 * the in-memory cached run installed below still drives this
		 * session's overlay. */
		char sidecar_path[2048];
		snprintf(sidecar_path, sizeof(sidecar_path), "%s.mpp", regargs->seq->seqname);
		if (mpp_sidecar_write(sidecar_path, run) != MPP_OK)
			siril_log_message(_("mpp: Stage A sidecar write to %s failed — "
			                    "APs will be lost when the sequence is closed\n"),
			                  sidecar_path);
		else
			siril_log_message(_("mpp: Stage A sidecar written to %s\n"), sidecar_path);
		/* Install the run in the GUI cache so end_register_idle can
		 * paint the mean frame and the AP overlay can pick up the AP
		 * grid. mpp_set_cached_run takes ownership and frees any prior
		 * cached run from a previous Analyze in this session. */
		gui_iface.set_progress(PROGRESS_DONE, _("Analyze: done"));
		mpp_set_cached_run(run);
		return MPP_OK;
	}

	/* If the AP-placement settings changed since this run's grid was built,
	 * regenerate the grid at the new settings before Stage B. Without this,
	 * Register would silently reuse the cached/loaded grid (e.g. the previous
	 * AP size) and ignore the spin buttons. Re-placement runs on the existing
	 * mean reference frame, so it does not re-read frames; it does discard
	 * manual AP-editor edits, which is the expected outcome of a global
	 * placement change. A fresh Analyze run has run->cfg == cfg, so this is a
	 * no-op there; it only fires for a reused cache or a sidecar-loaded run
	 * whose baked settings differ from the current widgets. */
	if (run->cfg && ap_placement_cfg_differs(run->cfg, &cfg)) {
		siril_log_message(_("mpp: alignment-point settings changed since "
		                    "analysis — regenerating the AP grid (any manual "
		                    "AP edits are discarded)\n"));
		const mpp_status_t rc_r = mpp_ap_replace(run, &cfg);
		if (rc_r != MPP_OK) {
			siril_log_error(_("mpp: AP grid regeneration failed (code %d)\n"), rc_r);
			gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
			if (run_owned) mpp_run_free(run);
			return rc_r;
		}
		siril_log_message(_("mpp: regenerated %d APs at the new settings\n"),
		                  run->aps->count);
	}

	/* Refresh inclusion from the live frame selector before ranking, so the
	 * per-AP top-N is taken over the user's current selection. Stage A
	 * already honours incl, so for a fresh Analyze this matches what it just
	 * read; for a reused cached run it picks up selection edits made since.
	 * If the selection changed under a run whose ranking is already baked,
	 * invalidate best_frame_indices so the recompute below re-ranks over the
	 * new set rather than early-returning. */
	if (regargs->seq->imgparam && run->num_frames == regargs->seq->number) {
		int n_selected = 0, changed = 0;
		for (int i = 0; i < run->num_frames; ++i) {
			const int live = regargs->seq->imgparam[i].incl ? 1 : 0;
			if (live != run->included[i]) changed = 1;
			run->included[i] = live;
			if (live) ++n_selected;
		}
		if (n_selected == 0) {
			siril_log_error(_("mpp: every frame is excluded — aborting\n"));
			gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
			if (run_owned) mpp_run_free(run);
			return MPP_ENODATA;
		}
		if (changed && run->best_frame_indices) {
			std::free(run->best_frame_indices);
			run->best_frame_indices = nullptr;
		}
	}

	/* One derotation plan for the whole of Stage B (per-AP ranking + shifts). */
	mpp_derot_t *derot_b = NULL;
	{
		gchar *dp = g_strdup_printf("%s.derot", regargs->seq->seqname);
		mpp_derot_read(dp, &derot_b);
		g_free(dp);
	}

	/* Per-AP frame ranking — deferred from Stage A. Runs here for a fresh
	 * Analyze (best_frame_indices null) and for a reused cached run (edited
	 * or not); an idempotent no-op for a sidecar-loaded run that already
	 * carries the ranking. */
	const int rc_q = mpp_recompute_qualities(regargs->seq, run, derot_b);
	if (rc_q != MPP_OK) {
		siril_log_error(_("mpp: per-AP frame ranking failed (code %d)\n"), rc_q);
		gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
		mpp_derot_free(derot_b);
		if (run_owned) mpp_run_free(run);
		return rc_q;
	}

	siril_log_message(_("mpp: Stage B — per-AP per-frame shifts\n"));
	gui_iface.set_progress(PROGRESS_RESET, _("Register: computing per-AP shifts"));
	const int rc_b = mpp_compute_shifts(regargs->seq, &cfg, run, derot_b);
	mpp_derot_free(derot_b);
	if (rc_b == MPP_EINTR) {
		siril_log_message(_("mpp: Stage B cancelled by user.\n"));
		gui_iface.set_progress(PROGRESS_DONE, _("Cancelled"));
		if (run_owned) mpp_run_free(run);
		return rc_b;
	}
	if (rc_b != MPP_OK) {
		siril_log_error(_("mpp: Stage B failed (code %d)\n"), rc_b);
		gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
		if (run_owned) mpp_run_free(run);
		return rc_b;
	}
	siril_log_message(_("mpp: Stage B done — %d shifts failed\n"),
	                  run->shifts ? run->shifts->failure_counter : -1);

	char sidecar_path[2048];
	snprintf(sidecar_path, sizeof(sidecar_path), "%s.mpp", regargs->seq->seqname);
	gui_iface.set_progress(PROGRESS_NONE, _("Register: writing sidecar"));
	if (mpp_sidecar_write(sidecar_path, run) != MPP_OK) {
		siril_log_message(_("mpp: sidecar write to %s failed\n"), sidecar_path);
		gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
		if (run_owned) mpp_run_free(run);
		return MPP_EIO;
	}
	siril_log_message(_("mpp: sidecar written to %s\n"), sidecar_path);
	gui_iface.set_progress(PROGRESS_DONE, _("Register: done"));

	/* Install the run in the GUI cache so end_register_idle's
	 * View APs / View Shifts sensitivity refresh — and any subsequent
	 * AP-editor or shift-viewer click — can find it. Same rationale as
	 * the stage_a_only branch above. When run_owned is false the cache
	 * already holds run from a prior Analyze and Stage B updated it in
	 * place, so there's nothing to install. */
	if (run_owned) mpp_set_cached_run(run);
	return MPP_OK;
}
