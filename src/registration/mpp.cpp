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
#include "registration/mpp/mpp_drizzle.h"
#include "registration/mpp/mpp_frames.h"
#include "registration/mpp/mpp_shift.h"
#include "registration/mpp/mpp_sidecar.h"
#include "registration/registration.h"
}

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

/* AVI Bayer pattern → OpenCV cvtColor code, with Y-flip awareness so the
 * pattern the user specifies (top-down, as it appears in the AVI raster)
 * stays correct after film_read_frame's fits_flip_top_to_bottom. Returns
 * -1 when no debayer should happen (AUTO/NONE/non-Bayer).
 *
 * OpenCV's 2-letter COLOR_BayerXX2RGB constants describe the 2×2 tile
 * by its second-row first-two pixels (the 4-letter RGGB/BGGR/GBRG/GRBG
 * aliases only exist on OpenCV ≥ 4.x and break the oldstable AppImage
 * build). We feed in the Mat we got back from the pipeline (bottom-up
 * after the Y-flip), so we ask Siril's helper to adjust the user's
 * top-down pick. */
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
	switch (pat) {
		case BAYER_FILTER_RGGB: return cv::COLOR_BayerBG2RGB;
		case BAYER_FILTER_BGGR: return cv::COLOR_BayerRG2RGB;
		case BAYER_FILTER_GBRG: return cv::COLOR_BayerGR2RGB;
		case BAYER_FILTER_GRBG: return cv::COLOR_BayerGB2RGB;
		default: return -1;
	}
}

/* Map an mpp_avi_bayer value to the Y-flip-aware cfa[4] used by the
 * Bayer-drizzle dobox path. Output channels are R=0, G=1, B=2. Returns
 * false when the pattern is AUTO / NONE / out-of-range. */
bool avi_bayer_to_cfa(int avi_pattern, int rows, unsigned char cfa[4]) {
	if (avi_pattern == MPP_AVI_BAYER_AUTO || avi_pattern == MPP_AVI_BAYER_NONE)
		return false;
	sensor_pattern pat;
	switch (avi_pattern) {
		case MPP_AVI_BAYER_RGGB: pat = BAYER_FILTER_RGGB; break;
		case MPP_AVI_BAYER_BGGR: pat = BAYER_FILTER_BGGR; break;
		case MPP_AVI_BAYER_GBRG: pat = BAYER_FILTER_GBRG; break;
		case MPP_AVI_BAYER_GRBG: pat = BAYER_FILTER_GRBG; break;
		default: return false;
	}
	adjust_Bayer_pattern_orientation(&pat, (unsigned int) rows, TRUE);
	const int R = 0, G = 1, B = 2;
	switch (pat) {
		case BAYER_FILTER_RGGB: cfa[0]=R; cfa[1]=G; cfa[2]=G; cfa[3]=B; return true;
		case BAYER_FILTER_GRBG: cfa[0]=G; cfa[1]=R; cfa[2]=B; cfa[3]=G; return true;
		case BAYER_FILTER_GBRG: cfa[0]=G; cfa[1]=B; cfa[2]=R; cfa[3]=G; return true;
		case BAYER_FILTER_BGGR: cfa[0]=B; cfa[1]=G; cfa[2]=G; cfa[3]=R; return true;
		default: return false;
	}
}

namespace {

/* Helper: read a frame and return a mono analysis Mat (cloned, owns its
 * memory after the fits is cleared). RGB input is converted to a
 * panchromatic luminance via cv::cvtColor(RGB, GRAY); mono input passes
 * through.
 *
 * For SEQ_AVI sequences with `avi_pattern` set to a Bayer value, the
 * 1-layer mosaic is OpenCV-debayered before luminance extraction so
 * Stage A sees a properly demosaiced image rather than a CFA-modulated
 * brightness map.
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
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	cv::Mat out;
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		if (view.empty()) return cv::Mat();
		if (seq && seq->type == SEQ_AVI) {
			const int code = mpp::avi_bayer_to_cv_code(avi_pattern, view.rows);
			if (code >= 0) {
				cv::Mat rgb;
				cv::cvtColor(view, rgb, code);
				cv::cvtColor(rgb, out, cv::COLOR_RGB2GRAY);
			} else {
				out = view.clone();
			}
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

}  /* close inner anonymous namespace — read_full_frame needs external
    * linkage so mpp_drizzle.cpp (the STScI stack path) can call it. */

/* Read the frame as a multi-channel cv::Mat preserving all layers (1 or 3).
 * Used by mpp_stack_apply and by mpp_stack_apply_stsci so colour frames
 * carry through to the stacked output.
 *
 * For SEQ_AVI with a Bayer pattern set, the mono mosaic is debayered to
 * 3-channel RGB before return so the non-drizzle stack paths see proper
 * colour input. The Bayer-drizzle path bypasses this helper. */
cv::Mat read_full_frame(sequence *seq, int idx, int avi_pattern) {
	FitsBuf buf;
	if (mpp_seq_read_frame(seq, idx, &buf.f, false, 0) != MPP_OK)
		return cv::Mat();
	const int n = fits_layer_count(&buf.f);
	if (n == 1) {
		cv::Mat view = wrap_fits_layer(&buf.f, 0);
		if (view.empty()) return cv::Mat();
		if (seq && seq->type == SEQ_AVI) {
			const int code = mpp::avi_bayer_to_cv_code(avi_pattern, view.rows);
			if (code >= 0) {
				cv::Mat rgb;
				cv::cvtColor(view, rgb, code);
				return rgb;
			}
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

namespace {   /* reopen anonymous namespace for the remaining helpers */

/* Per-AP target stack size from cfg (number override wins, else
 * round-up % of N), clamped to [1, N]. Mirrors the formula in
 * ap_compute_frame_qualities_streamed so Stage C entry can honour
 * Stack-tab overrides against the value Stage A baked at Analyze
 * time. */
int stack_size_from_cfg(const mpp_config_t &cfg, int N) {
	if (N <= 0) return 0;
	int s;
	if (cfg.alignment_points_frame_number > 0)
		s = cfg.alignment_points_frame_number;
	else
		s = std::max(1, (int) std::ceil(
		        (double) N * cfg.alignment_points_frame_percent / 100.0));
	return std::min(s, N);
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
#define MPP_PROG_READ_END         0.40   /* parallel frame read */
#define MPP_PROG_GLOBAL_END       0.60   /* per-frame correlation */
#define MPP_PROG_AVERAGE_END      0.65   /* per-frame accumulator */
#define MPP_PROG_QUALITIES_END    1.00   /* per-frame Laplace + per-AP stddev */

extern "C" mpp_status_t mpp_analyze(sequence *seq, const mpp_config_t *cfg,
                                    mpp_run_t **run_out) {
	if (!seq || !cfg || !run_out)
		return MPP_EINVAL;
	if (seq->number <= 0)
		return MPP_ENODATA;

	const int N = seq->number;
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
	 * per_thread_in_flight = 2: each Pass 1 iteration holds the raw mono
	 * plus a transient blur (inside rank_score_normalized) at peak.
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
		    /*per_thread_in_flight=*/2,
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
		if (g_atomic_int_get(&read_failed) || g_atomic_int_get(&cancelled)) continue;
		if (!processing_should_continue()) {
			g_atomic_int_set(&cancelled, 1);
			continue;
		}
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
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(MPP_PROG_READ_END * (double) cur_nb / (double) N, NULL);
	}
	g_mutex_clear(&cache_mtx);
	if (cancelled) return MPP_EINTR;
	if (read_failed) return MPP_EIO;

	num_layers = seq->nb_layers > 0 ? seq->nb_layers : 1;
	bitpix = seq->bitpix;

	/* Unified providers for Passes 2/3/4. Hit the analysis_frames slot
	 * if populated; otherwise re-read from disk. Works identically for
	 * HALF_CACHED (every slot hits), PARTIAL_CACHED (first K hit), and
	 * STREAMING (every slot misses since analysis_frames is empty). */
	auto raw_provider = [&analysis_frames, seq, cfg](int i) -> cv::Mat {
		if (i >= 0 && i < (int) analysis_frames.size()
		 && !analysis_frames[i].empty())
			return analysis_frames[i];
		return mpp::read_analysis_frame(seq, i, cfg->avi_bayer_pattern, cfg->bitdepth);
	};
	auto blurred_provider = [&analysis_frames, &blurred_frames, seq, cfg](int i) -> cv::Mat {
		if (i >= 0 && i < (int) blurred_frames.size()
		 && !blurred_frames[i].empty())
			return blurred_frames[i];
		cv::Mat m;
		if (i >= 0 && i < (int) analysis_frames.size()
		 && !analysis_frames[i].empty())
			m = analysis_frames[i];
		else
			m = mpp::read_analysis_frame(seq, i, cfg->avi_bayer_pattern, cfg->bitdepth);
		return m.empty() ? cv::Mat() : mpp::blur_mono_for_align(m, *cfg);
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
	const auto align = mpp::align_global_from_provider(
	    blurred_provider, N, q_rank, *cfg, stage_progress_cb, &gp_global,
	    streaming_provider_safe);
	/* Drop the blurred cache now — Passes 3/4 don't use it, and the
	 * freed bytes give them more headroom (matters most when the blur
	 * top-up was large). */
	blurred_frames.clear();
	blurred_frames.shrink_to_fit();
	if (align.best_frame_idx < 0) return MPP_EINTR;   /* cancellation sentinel */
	gui_iface.set_progress(MPP_PROG_GLOBAL_END, _("Analyze: building reference frame"));
	stage_progress_range gp_avg{MPP_PROG_GLOBAL_END,
	                            MPP_PROG_AVERAGE_END - MPP_PROG_GLOBAL_END};
	/* Pass 3 — average frame. raw_provider serves cache hits or disk
	 * reads transparently; the streamed variant is sequential so thread
	 * safety doesn't apply. */
	const auto avg = mpp::align_average_frame_streamed(
	    raw_provider, N, frame_rows, frame_cols, q_rank, align.shifts, *cfg,
	    stage_progress_cb, &gp_avg);
	if (avg.mean_frame.empty()) return MPP_EINTR;
	gui_iface.set_progress(PROGRESS_NONE, _("Analyze: placing alignment points"));
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, *cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, *cfg);
	if (!aps || aps->count <= 0) {
		mpp_ap_free(aps);
		return MPP_ENODATA;
	}
	gui_iface.set_progress(MPP_PROG_AVERAGE_END, _("Analyze: per-AP frame qualities"));
	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	stage_progress_range gp_qual{MPP_PROG_AVERAGE_END,
	                             MPP_PROG_QUALITIES_END - MPP_PROG_AVERAGE_END};
	/* Pass 4 — per-AP qualities. raw_provider handles cache hits vs disk
	 * reads; full CACHED/HALF_CACHED see 100% hits, PARTIAL sees K/N. */
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, frame_brightness, *aps, offsets,
	    frame_rows, frame_cols, *cfg, stage_progress_cb, &gp_qual);
	if (apq.stack_size <= 0) { mpp_ap_free(aps); return MPP_EINTR; }

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
	run->global_shifts    = (double *) std::malloc(2 * N * sizeof(double));
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

	*run_out = run;
	return MPP_OK;
}

/* -------------------- Stage B: mpp_compute_shifts -------------------- */

extern "C" mpp_status_t mpp_compute_shifts(sequence *seq, const mpp_config_t *cfg,
                                           mpp_run_t *run) {
	if (!seq || !cfg || !run || !run->aps) return MPP_EINVAL;

	/* Stage B streams: each frame is read + blurred inside the
	 * correlation loop and discarded at iteration end. Sequential, so
	 * max_threads=1 and per_thread_in_flight=1 — the picker just
	 * verifies that a single frame fits. Disk traffic is identical to
	 * the previous cached two-pass design (1 read per included frame). */
	mpp_mem_pick mem_pick{MPP_MEM_STREAMING, 1, 0, 0};
	{
		const int bps = (cfg->bitdepth == 8) ? 1 : 2;
		const mpp_status_t mem = mpp_pick_memory_mode(
		    run->num_frames, (int) seq->rx, (int) seq->ry, /*channels=*/1,
		    bps,
		    /*cached_copies=*/0,
		    /*half_cached_copies=*/0,
		    /*max_threads=*/1,
		    /*per_thread_in_flight=*/1,
		    /*output_bytes=*/0,
		    _("Register"), &mem_pick);
		if (mem != MPP_OK) return mem;
	}
	(void) mem_pick;

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
	auto provider = [seq, cfg, cache, cache_N](int f) -> cv::Mat {
		cv::Mat mono;
		if (cache && f >= 0 && f < cache_N
		 && !cache->raw_frames[f].empty()) {
			mono = cache->raw_frames[f];
		} else {
			mono = mpp::read_analysis_frame(seq, f, cfg->avi_bayer_pattern,
			                                cfg->bitdepth);
		}
		if (mono.empty()) return cv::Mat();
		return mpp::blur_mono_for_align(mono, *cfg);
	};
	mpp_shifts_t *shifts = mpp::stack_compute_shifts_streamed(
	    provider, run->num_frames, mean_raw, *run->aps, apq, offsets,
	    *cfg, run->included);
	if (!shifts) return MPP_ENOMEM;
	if (shifts->failure_counter == -1) {   /* cancellation sentinel */
		mpp_shift_free(shifts);
		return MPP_EINTR;
	}
	if (run->shifts) mpp_shift_free(run->shifts);
	run->shifts = shifts;
	return MPP_OK;
}

/* -------------------- Stage C: mpp_stack_apply -------------------- */

/* RAII helper: temporarily flips a CFA SER's debayer_type_ser from
 * SER_MONO to its actual color_id (and bumps nb_layers from 1 to 3) for
 * the duration of an analysis pass, restoring on destruction. The rank
 * measure (GaussianBlur 7×7 → Laplace σ) cannot operate sensibly on a
 * raw Bayer mosaic — the CFA's checker pattern dominates the
 * measurement, masking the actual seeing-driven sharpness signal — so
 * analysis always needs debayered input. The bayer-* drizzle stack
 * path re-reads raw bytes via its own wrapper regardless of this
 * setting, so flipping for analysis doesn't disturb the raw-mosaic
 * path that bayer-drizzle relies on at stack time. */
namespace {
struct SerAnalysisDebayerGuard {
	sequence *seq = nullptr;
	ser_color saved_debayer_type = SER_MONO;
	int saved_nb_layers = -1;
	bool active = false;

	~SerAnalysisDebayerGuard() {
		if (active && seq && seq->ser_file) {
			seq->ser_file->debayer_type_ser = saved_debayer_type;
			seq->nb_layers = saved_nb_layers;
		}
	}
};

bool maybe_engage_ser_debayer_for_analysis(SerAnalysisDebayerGuard &g, sequence *seq) {
	if (!seq || seq->type != SEQ_SER || !seq->ser_file) return false;
	if (seq->nb_layers != 1) return false;
	const ser_color cid = seq->ser_file->color_id;
	if (cid < SER_BAYER_RGGB || cid > SER_BAYER_BGGR) return false;
	g.seq = seq;
	g.saved_debayer_type = seq->ser_file->debayer_type_ser;
	g.saved_nb_layers = seq->nb_layers;
	g.active = true;
	seq->ser_file->debayer_type_ser = cid;
	seq->nb_layers = 3;
	return true;
}
}  // namespace

extern "C" mpp_input_type mpp_classify_sequence_input(const sequence *seq) {
	if (!seq) return MPP_INPUT_MONO;
	if (seq->type == SEQ_SER && seq->ser_file) {
		const ser_color cid = seq->ser_file->color_id;
		if (cid >= SER_BAYER_RGGB && cid <= SER_BAYER_BGGR) return MPP_INPUT_CFA;
		if (cid == SER_RGB || cid == SER_BGR)               return MPP_INPUT_RGB;
		return MPP_INPUT_MONO;
	}
	/* FITS / other: nb_layers heuristic. A CFA-FITS sequence with
	 * nb_layers=1 + bayer_pattern keyword would be misclassified as
	 * MONO here — that's a corner case to address if it shows up in the
	 * field; the SER path covers the common case. */
	if (seq->nb_layers > 1) return MPP_INPUT_RGB;
	return MPP_INPUT_MONO;
}

extern "C" mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                                        mpp_run_t *run, fits *out) {
	if (!seq || !cfg || !run || !run->aps || !run->shifts || !out)
		return MPP_EINVAL;

	/* Honour Stack-tab per-AP top-K overrides. Stage A baked
	 * run->stack_size from the Analyze-time cfg; the GUI Stack tab and
	 * sidecar overrides land in cfg here without affecting run->stack_size
	 * automatically. Re-derive the target from cfg and clamp:
	 *   target ≤ baked: shrink in-place — best_frame_indices is sorted
	 *     descending by per-AP quality, so taking the first `target`
	 *     entries per AP is correct (used_alignment_points naturally
	 *     follows when apq_from_run reads the new stack_size).
	 *   target > baked: warn and keep baked. Analyze discarded the per-AP
	 *     quality matrix past the top run->stack_size, so we can't grow
	 *     without re-running Stage A. */
	if (run->num_frames > 0) {
		const int target = mpp::stack_size_from_cfg(*cfg, run->num_frames);
		if (target > run->stack_size) {
			siril_log_warning(
			    _("Stack (mpp): requested %d top-quality frames per AP but "
			      "Analyze baked %d — keeping %d. Re-run Analyze with a "
			      "higher frame %%/count to raise the ceiling.\n"),
			    target, run->stack_size, run->stack_size);
		} else if (target < run->stack_size) {
			run->stack_size = target;
		}
		/* Always announce the effective per-AP stack size so the user
		 * can verify Stack-tab/Register-tab settings are taking effect
		 * (silence when stack_size == num_frames means no per-AP filter
		 * is active and the message would be noise). */
		if (run->stack_size < run->num_frames) {
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

	/* Dispatch: the picked scaling method determines the backend.
	 *   - MPP_KERNEL_UPSCALE: cv::resize bicubic path. At scale=1 this
	 *     is a no-op upscale + classical-debayer-and-average for CFA
	 *     input. At scale > 1 it's a plain bicubic upscale + average,
	 *     the reliable fallback for cases where true drizzle resonates
	 *     on fine high-contrast structure.
	 *   - Any drizzle kernel: dobox STScI path. CFA inputs are fed in
	 *     pre-debayered (SerAnalysisDebayerGuard below forces the SER
	 *     reader to debayer regardless of the user's debayer-on-open
	 *     pref) — STScI's pixmap accumulator gets full-RGB input pixels,
	 *     no CFA-phase coupling.
	 *
	 * The MPP_DRIZZLE_BAYER backend (CFA-position-aware splat with
	 * dobox) is preserved in-tree (mpp_stack_apply_bayer +
	 * stack_apply_bayer_streamed in mpp_drizzle.cpp, tests in
	 * mpp_drizzle_test.cpp) but is no longer reachable through normal
	 * dispatch. It depends on sub-pixel diversity being faithfully
	 * captured by the alignment chain; with the current bilinear-
	 * demosaic-then-grayscale analysis pipeline, phase-1 stride-2
	 * subsampling (= CFA period) and the debayer-induced bias in the
	 * correlation surface together suppress that diversity for
	 * low-jitter planetary captures, producing the wavelet-amplified
	 * 2-pixel-period stripe artefacts users hit on real data. A future
	 * CFA-aware correlation pass (operating directly on the raw
	 * mosaic per CFA channel) could revive this path; until then,
	 * STScI on debayered input is the correct default.
	 *
	 * cfg->drizzle_mode is honoured if explicitly set (sidecar /
	 * scripts); the default MPP_DRIZZLE_OFF triggers the kernel-driven
	 * routing here. */
	int mode = cfg->drizzle_mode;
	/* AVI Bayer routing: an AVI marked as a raw mosaic (cfg->avi_bayer_pattern
	 * RGGB/BGGR/GBRG/GRBG) is the AVI analogue of MPP_INPUT_CFA. The
	 * classifier can't see this (sequence-only signature) so it's
	 * tracked here for the channel-count math below. */
	const bool avi_cfa = (seq->type == SEQ_AVI
	                      && cfg->avi_bayer_pattern >= MPP_AVI_BAYER_RGGB
	                      && cfg->avi_bayer_pattern <= MPP_AVI_BAYER_GRBG);
	if (mode == MPP_DRIZZLE_OFF) {
		if (cfg->drizzle_kernel == MPP_KERNEL_UPSCALE) {
			/* Stay on MPP_DRIZZLE_OFF — falls through to the cv::resize
			 * bicubic path below. At scale=1 cv::resize is skipped
			 * entirely (K=1) and frames are stacked at native res. */
		} else {
			/* Everything else → STScI dobox. CFA inputs are pre-debayered
			 * by the guard below; non-CFA inputs (mono / native RGB SER /
			 * FITS) flow straight through. */
			mode = MPP_DRIZZLE_STSCI;
		}
	} else if (mode == MPP_DRIZZLE_BAYER) {
		/* Sidecar pinned the retired Bayer backend (or a script set
		 * drizzle_mode=BAYER explicitly). Promote to STScI rather
		 * than reject — the user's intent ("drizzle this CFA
		 * sequence") is honoured, just via the more robust path. */
		siril_log_warning(_("Stack (mpp): Bayer drizzle backend retired "
		                    "(see commit notes); routing to STScI on "
		                    "auto-debayered input instead.\n"));
		mode = MPP_DRIZZLE_STSCI;
	}

	/* Auto-debayer for any CFA SER opened raw — STScI, Bayer, and the
	 * bicubic fallthrough all need RGB input now. The guard restores
	 * the prior debayer_type_ser on destruction. */
	SerAnalysisDebayerGuard ser_guard;
	maybe_engage_ser_debayer_for_analysis(ser_guard, seq);

	/* Pre-flight memory check. All Stack paths stream sequentially —
	 * the dobox wrappers and the bicubic loop below both hold one
	 * in-flight frame at a time, plus the drizzle output canvases
	 * (small relative to per-frame). */
	{
		/* AVI Bayer mosaics decode to 3-channel via cv::cvtColor in the
		 * bicubic path; the retired Bayer drizzle path is the only one
		 * that would have read 1-channel; everything else carries
		 * channels through. */
		const int channels = (mode == MPP_DRIZZLE_BAYER)
		                   ? 1
		                   : (avi_cfa ? 3
		                              : (run->num_layers > 0 ? run->num_layers : 1));
		mpp_mem_pick mem_pick{MPP_MEM_STREAMING, 1, 0, 0};
		/* Stack uses read_full_frame (CV_16U) regardless of source depth
		 * so the accumulator retains the full dynamic range. */
		const mpp_status_t mem = mpp_pick_memory_mode(
		    run->num_frames, run->frame_cols, run->frame_rows, channels,
		    /*bytes_per_sample=*/2,
		    /*cached_copies=*/0,
		    /*half_cached_copies=*/0,
		    /*max_threads=*/1,
		    /*per_thread_in_flight=*/1,
		    /*output_bytes=*/0,
		    _("Stack (mpp)"), &mem_pick);
		if (mem != MPP_OK) return mem;
		(void) mem_pick;
	}
	switch (mode) {
		case MPP_DRIZZLE_STSCI: return mpp_stack_apply_stsci(seq, cfg, run, out);
		case MPP_DRIZZLE_BAYER:
			/* Unreachable in normal use — the OFF and explicit-BAYER
			 * branches above both promote to STSCI. Kept as a switch
			 * arm so the enum value stays callable for tests and any
			 * future CFA-aware revival. */
			return mpp_stack_apply_bayer(seq, cfg, run, out);
		case MPP_DRIZZLE_OFF:
		default:
			break;   /* fall through to the cv::resize no-upscale path */
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

	gui_iface.set_progress(PROGRESS_RESET, _("Stack: drizzling frames"));
	auto provider = [seq, cfg](int i) -> cv::Mat {
		return mpp::read_full_frame(seq, i, cfg->avi_bayer_pattern);
	};
	/* read_full_frame's Bayer path returns 3-channel RGB for AVI mosaics,
	 * even though run->num_layers (captured from seq->nb_layers at Stage A
	 * time) is 1. The accumulator type depends on num_layers, so override
	 * to 3 for the AVI Bayer case. */
	int num_layers = run->num_layers > 0 ? run->num_layers : 1;
	if (avi_cfa) num_layers = 3;
	const auto loop = mpp::stack_apply_shifts_streamed(
	    provider, run->num_frames, num_layers, *run->aps, apq,
	    run->shifts, offsets, frame_brightness, sorted_idx,
	    intersection, *cfg, run->included);
	if (loop.state.dim_y == 0) return MPP_EINTR;   /* cancellation sentinel */
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
	/* Siril's mono convention: pdata[1]/pdata[2] alias pdata[0] (not
	 * NULL). remap_all_vports's per-vport memcpy loops c=0..2
	 * unconditionally, so a NULL plane pointer crashes the GUI
	 * refresh that runs immediately after stack_function_handler
	 * swaps gfit. */
	out->pdata[1] = (C >= 2) ? out->data + plane     : out->data;
	out->pdata[2] = (C >= 3) ? out->data + plane * 2 : out->data;
	return MPP_OK;
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

extern "C" mpp_status_t mpp_recompute_qualities(sequence *seq, mpp_run_t *run) {
	if (!seq || !run || !run->aps || run->aps->count <= 0 || !run->cfg)
		return MPP_EINVAL;
	if (run->best_frame_indices) return MPP_OK;   /* nothing to do */
	if (!run->frame_brightness || !run->global_shifts) return MPP_EINVAL;
	if (run->num_frames != seq->number) return MPP_EINVAL;

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
		    /*per_thread_in_flight=*/1,
		    /*output_bytes=*/0,
		    _("Register (per-AP quality refresh)"), &mem_pick);
		if (mem != MPP_OK) return mem;
		fresh_cache_count = mem_pick.cached_count;
		const int preread_threads = mem_pick.threads;

		siril_log_message(_("mpp: recomputing per-AP frame qualities for edited AP grid "
		                    "(re-reading %d frames)\n"), N);
		gui_iface.set_progress(PROGRESS_RESET, _("Register: re-reading frames for edited APs"));

		if (fresh_cache_count > 0) {
			fresh_cache_vec.resize(N);
			int cur_nb = 0;
			int read_failed = 0;
			int cancelled = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(preread_threads) schedule(guided) \
    if (seq->type == SEQ_SER || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ \
                                  || seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
			for (int i = 0; i < fresh_cache_count; ++i) {
				if (g_atomic_int_get(&read_failed) || g_atomic_int_get(&cancelled)) continue;
				if (!processing_should_continue()) {
					g_atomic_int_set(&cancelled, 1);
					continue;
				}
				cv::Mat mono = mpp::read_analysis_frame(seq, i, run->cfg->avi_bayer_pattern,
				                                        run->cfg->bitdepth);
				if (mono.empty()) {
					g_atomic_int_set(&read_failed, 1);
					continue;
				}
				fresh_cache_vec[i] = mono;
				g_atomic_int_inc(&cur_nb);
				gui_iface.set_progress(0.5 * (double) cur_nb / (double) fresh_cache_count, NULL);
			}
			if (cancelled)  return MPP_EINTR;
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
		siril_log_message(_("mpp: recomputing per-AP frame qualities for edited AP grid "
		                    "(%d/%d frames already cached)\n"),
		                  run->cache->populated_count, N);
		gui_iface.set_progress(PROGRESS_RESET, _("Register: re-using cached analysis frames"));
	}

	gui_iface.set_progress(0.5, _("Register: computing per-AP qualities for edited grid"));
	const std::vector<double> frame_brightness(run->frame_brightness,
	                                           run->frame_brightness + N);
	const auto offsets = mpp::offsets_from_run(run);
	stage_progress_range rp{0.5, 0.5};
	const mpp_cache_t *cache = run->cache;
	const int cache_N = cache ? (int) cache->raw_frames.size() : 0;
	auto raw_provider = [cache, cache_N, seq, run](int i) -> cv::Mat {
		if (cache && i >= 0 && i < cache_N
		 && !cache->raw_frames[i].empty())
			return cache->raw_frames[i];
		return mpp::read_analysis_frame(seq, i, run->cfg->avi_bayer_pattern,
		                                run->cfg->bitdepth);
	};
	const auto apq = mpp::ap_compute_frame_qualities_streamed(
	    raw_provider, N, frame_brightness, *run->aps, offsets,
	    run->frame_rows, run->frame_cols, *run->cfg,
	    stage_progress_cb, &rp);
	if (apq.stack_size <= 0) return MPP_EINTR;

	const int stack_size = apq.stack_size > 0 ? apq.stack_size : run->stack_size;
	const size_t total = (size_t) run->aps->count * (size_t) stack_size;
	int *bfi = (int *) std::malloc(total * sizeof(int));
	if (!bfi) return MPP_ENOMEM;
	for (int a = 0; a < run->aps->count; ++a) {
		for (int k = 0; k < stack_size; ++k) {
			bfi[a * stack_size + k] =
			    (k < (int) apq.best_frame_indices[a].size())
			        ? apq.best_frame_indices[a][k]
			        : -1;
		}
	}
	std::free(run->best_frame_indices);
	run->best_frame_indices = bfi;
	run->stack_size         = stack_size;
	siril_log_message(_("mpp: per-AP qualities refreshed for %d edited APs\n"),
	                  run->aps->count);
	return MPP_OK;
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

/* Fill the box/patch bounds for an AP at (x, y) from cfg, extending the
 * patch out to the frame border for edge APs. */
static void fill_ap_bounds(mpp_ap_record_t *r, int x, int y,
                           const mpp_config_t *cfg,
                           int frame_rows, int frame_cols) {
	const int hb = cfg->alignment_points_half_box_width;
	const int hp = mpp_cfg_half_patch_width(cfg);
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

extern "C" int register_mpp(struct registration_args *regargs) {
	if (!regargs || !regargs->seq)
		return MPP_EINVAL;

	/* Auto-debayer for analysis when a CFA SER was opened raw (i.e.
	 * Preferences → Debayering → "Debayer SER files at opening" is OFF).
	 * The guard flips the SER state for the duration of this function
	 * and restores on every return path. See SerAnalysisDebayerGuard
	 * above for the rationale. */
	SerAnalysisDebayerGuard ser_guard;
	if (maybe_engage_ser_debayer_for_analysis(ser_guard, regargs->seq)) {
		siril_log_warning(
		    _("mpp: auto-debayering CFA SER for analysis.\n"));
	}

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
		                            "(%d APs%s).\n"),
		                        cached->aps->count,
		                        cached->best_frame_indices ? "" : ", AP-edited");
		if (stage_a_only) {
			/* Re-Analyze with cached run: drop the cache here so the fresh
			 * mpp_analyze below builds a new run that supersedes it. */
			run = nullptr;
			run_owned = true;
		} else if (!cached->best_frame_indices) {
			/* APs were edited since Analyze — refresh per-AP qualities
			 * for the edited grid before Stage B. Re-reads analysis frames
			 * but skips rank/global-align/mean-frame build (cached). */
			const int rc_q = mpp_recompute_qualities(regargs->seq, cached);
			if (rc_q != MPP_OK) {
				siril_log_error(_("mpp: failed to recompute per-AP qualities "
				                          "for edited grid (code %d)\n"), rc_q);
				gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
				return rc_q;
			}
		}
	}

	if (!run) {
		gui_iface.set_progress(PROGRESS_RESET,
		                       stage_a_only ? _("Analyze: ranking frames")
		                                    : _("Register: ranking frames"));
		const int rc_a = mpp_analyze(regargs->seq, &cfg, &run);
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
		/* Install the run in the GUI cache so end_register_idle can
		 * paint the mean frame and the AP overlay can pick up the AP
		 * grid. mpp_set_cached_run takes ownership and frees any prior
		 * cached run from a previous Analyze in this session. */
		gui_iface.set_progress(PROGRESS_DONE, _("Analyze: done"));
		mpp_set_cached_run(run);
		return MPP_OK;
	}

	/* Honour the GUI's "Selected images only" combo (CLI -selected) for
	 * Stages B and C. Stage A always processes every frame so the rank
	 * data covers the whole sequence and the user can review it. */
	if (regargs->filters.filter_included) {
		int n_selected = 0;
		for (int i = 0; i < run->num_frames && i < regargs->seq->number; ++i) {
			run->included[i] = regargs->seq->imgparam[i].incl ? 1 : 0;
			if (run->included[i]) ++n_selected;
		}
		siril_log_message(_("mpp: filtering to %d selected frames of %d total\n"),
		                  n_selected, run->num_frames);
		if (n_selected == 0) {
			siril_log_error(_("mpp: no frames selected — aborting\n"));
			gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
			if (run_owned) mpp_run_free(run);
			return MPP_ENODATA;
		}
	}

	siril_log_message(_("mpp: Stage B — per-AP per-frame shifts\n"));
	gui_iface.set_progress(PROGRESS_RESET, _("Register: computing per-AP shifts"));
	const int rc_b = mpp_compute_shifts(regargs->seq, &cfg, run);
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
