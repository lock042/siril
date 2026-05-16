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

/* core/siril.h pulls in <omp.h>, which under GCC 15's libgomp declares
 * C++ templates; those can't sit inside an `extern "C"` block. Include
 * it first so the header guard fires before the wrap below. */
#include "core/siril.h"

extern "C" {
#include "core/siril_log.h"
#include "core/proto.h"
#include "core/gui_iface.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/ser.h"

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_drizzle.h"
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

}  /* close inner anonymous namespace — read_full_frame needs external
    * linkage so mpp_drizzle.cpp (the STScI stack path) can call it. */

/* Read the frame as a multi-channel cv::Mat preserving all layers (1 or 3).
 * Used by mpp_stack_apply and by mpp_stack_apply_stsci so colour frames
 * carry through to the stacked output. */
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

namespace {   /* reopen anonymous namespace for the remaining helpers */

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
	/* Index-stable storage so parallel threads can write to per-frame slots
	 * without contention. push_back would force serialisation. */
	std::vector<cv::Mat> analysis_frames(N);
	std::vector<cv::Mat> blurred_frames(N);
	std::vector<double> q_rank(N, 0.0);
	std::vector<double> frame_brightness(N, 0.0);

	int frame_rows = 0, frame_cols = 0, num_layers = 1, bitpix = 16;
	int cur_nb = 0;
	int read_failed = 0;
	/* Parallel frame read: SER's fd_lock serialises only the actual fread
	 * call inside ser_read_frame so the rest (decompress, debayer, cvtColor,
	 * GaussianBlur, Laplace) runs concurrently. Pattern + guard mirror
	 * shift_methods.c:432 — anything fits-backed needs cfitsio's
	 * fits_is_reentrant() flag. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
    if (seq->type == SEQ_SER || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ \
                                  || seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
	for (int i = 0; i < N; ++i) {
		if (g_atomic_int_get(&read_failed)) continue;
		const cv::Mat mono = mpp::read_analysis_frame(seq, i);
		if (mono.empty()) {
			g_atomic_int_set(&read_failed, 1);
			continue;
		}
		analysis_frames[i] = mono;
		blurred_frames[i]  = mpp::blur_mono_for_align(mono, *cfg);
		q_rank[i]          = mpp::rank_score_normalized(mono, *cfg);
		frame_brightness[i] = mpp::rank_average_brightness(mono, *cfg);
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(MPP_PROG_READ_END * (double) cur_nb / (double) N, NULL);
	}
	if (read_failed) return MPP_EIO;
	frame_rows = analysis_frames[0].rows;
	frame_cols = analysis_frames[0].cols;

	num_layers = seq->nb_layers > 0 ? seq->nb_layers : 1;
	bitpix = seq->bitpix;

	gui_iface.set_progress(MPP_PROG_READ_END, _("Analyze: aligning frames globally"));
	stage_progress_range gp_global{MPP_PROG_READ_END,
	                               MPP_PROG_GLOBAL_END - MPP_PROG_READ_END};
	const auto align = mpp::align_global_from_frames(blurred_frames, q_rank, *cfg,
	                                                 stage_progress_cb, &gp_global);
	gui_iface.set_progress(MPP_PROG_GLOBAL_END, _("Analyze: building reference frame"));
	stage_progress_range gp_avg{MPP_PROG_GLOBAL_END,
	                            MPP_PROG_AVERAGE_END - MPP_PROG_GLOBAL_END};
	const auto avg = mpp::align_average_frame(analysis_frames, q_rank,
	                                          align.shifts, *cfg,
	                                          stage_progress_cb, &gp_avg);
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
	const auto apq = mpp::ap_compute_frame_qualities(analysis_frames, frame_brightness,
	                                                 *aps, offsets,
	                                                 frame_rows, frame_cols, *cfg,
	                                                 stage_progress_cb, &gp_qual);

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
	blurred_frames.resize(run->num_frames);  /* index-stable; empty Mat for excluded */
	int n_included = 0;
	for (int i = 0; i < run->num_frames; ++i)
		if (run->included[i]) ++n_included;
	/* The frame-read pass is roughly half the Stage B cost; report it as
	 * the first half of the [0..1] progress range, leaving the second half
	 * for the inner stack_compute_shifts (which emits its own 0.5..1). */
	int cur_nb = 0;
	for (int i = 0; i < run->num_frames; ++i) {
		if (!run->included[i]) continue;  /* skip frame entirely — keep empty slot */
		const cv::Mat mono = mpp::read_analysis_frame(seq, i);
		if (mono.empty()) return MPP_EIO;
		blurred_frames[i] = mpp::blur_mono_for_align(mono, *cfg);
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(0.5 * (double) cur_nb / (double) n_included, NULL);
	}

	const auto apq = mpp::apq_from_run(run);
	const auto offsets = mpp::offsets_from_run(run);
	const cv::Mat mean_raw = mpp::mean_frame_from_run(run);

	gui_iface.set_progress(PROGRESS_NONE, _("Register: correlating per-AP boxes"));
	mpp_shifts_t *shifts = mpp::stack_compute_shifts(blurred_frames, mean_raw,
	                                                 *run->aps, apq, offsets,
	                                                 *cfg, run->included);
	if (!shifts) return MPP_ENOMEM;
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
                                        const mpp_run_t *run, fits *out) {
	if (!seq || !cfg || !run || !run->aps || !run->shifts || !out)
		return MPP_EINVAL;

	/* Phase 5b dispatch: STScI / Bayer paths route through dobox via the
	 * pixmap built in mpp_drizzle.cpp. Bicubic (and the bicubic-with-
	 * drizzle_factor=1 default) stays on the Phase 5a cv::resize loop. */
	switch (cfg->drizzle_mode) {
		case MPP_DRIZZLE_STSCI: return mpp_stack_apply_stsci(seq, cfg, run, out);
		case MPP_DRIZZLE_BAYER: return mpp_stack_apply_bayer(seq, cfg, run, out);
		case MPP_DRIZZLE_OFF:
		case MPP_DRIZZLE_BICUBIC:
		default:
			break;   /* fall through to the Phase 5a bicubic path below */
	}

	/* Auto-debayer for the bicubic path when a CFA SER was opened raw.
	 * Without this, mpp::read_full_frame would return a 1-channel mosaic
	 * and the bicubic stack would produce a single-channel "image" that
	 * is really the CFA pattern smoothed by INTER_LINEAR — visually
	 * grey-scale and structurally wrong as a colour stack. The bayer-*
	 * path above doesn't need the guard because its wrapper reads the
	 * raw mosaic directly and routes per-CFA-position to the matching
	 * output channel. */
	SerAnalysisDebayerGuard ser_guard;
	maybe_engage_ser_debayer_for_analysis(ser_guard, seq);

	/* Read frames preserving all layers — colour data flows through Stage C
	 * even though Stages A and B use only the green analysis layer. */
	std::vector<cv::Mat> raw_frames;
	raw_frames.resize(run->num_frames);
	gui_iface.set_progress(PROGRESS_RESET, _("Stack: reading frames"));
	int n_included = 0;
	for (int i = 0; i < run->num_frames; ++i)
		if (run->included[i]) ++n_included;
	int cur_nb = 0;
	for (int i = 0; i < run->num_frames; ++i) {
		if (!run->included[i]) continue;
		const cv::Mat frame = mpp::read_full_frame(seq, i);
		if (frame.empty()) return MPP_EIO;
		raw_frames[i] = frame;
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(0.5 * (double) cur_nb / (double) n_included, NULL);
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

	const auto loop = mpp::stack_apply_shifts(raw_frames, *run->aps, apq,
	                                          run->shifts, offsets,
	                                          frame_brightness, sorted_idx,
	                                          intersection, *cfg, run->included);
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

	siril_log_message(_("mpp: recomputing per-AP frame qualities for edited AP grid "
	                    "(re-reading %d frames)\n"), N);
	gui_iface.set_progress(PROGRESS_RESET, _("Register: re-reading frames for edited APs"));

	/* Re-read analysis frames using the same parallel pattern as Stage A. */
	std::vector<cv::Mat> analysis_frames(N);
	int cur_nb = 0;
	int read_failed = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
    if (seq->type == SEQ_SER || ((seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ \
                                  || seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
	for (int i = 0; i < N; ++i) {
		if (g_atomic_int_get(&read_failed)) continue;
		cv::Mat mono = mpp::read_analysis_frame(seq, i);
		if (mono.empty()) {
			g_atomic_int_set(&read_failed, 1);
			continue;
		}
		analysis_frames[i] = mono;
		g_atomic_int_inc(&cur_nb);
		gui_iface.set_progress(0.5 * (double) cur_nb / (double) N, NULL);
	}
	if (read_failed) return MPP_EIO;

	gui_iface.set_progress(0.5, _("Register: computing per-AP qualities for edited grid"));
	const std::vector<double> frame_brightness(run->frame_brightness,
	                                           run->frame_brightness + N);
	const auto offsets = mpp::offsets_from_run(run);
	stage_progress_range rp{0.5, 0.5};
	const auto apq = mpp::ap_compute_frame_qualities(analysis_frames, frame_brightness,
	                                                 *run->aps, offsets,
	                                                 run->frame_rows, run->frame_cols,
	                                                 *run->cfg,
	                                                 stage_progress_cb, &rp);

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

/* Fill the box/patch bounds for an AP at (x, y) from cfg. Mirrors PSS's
 * extend-to-border logic at the frame edges. */
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

	/* Auto-debayer for analysis when a CFA SER was opened raw (i.e.
	 * Preferences → Debayering → "Debayer SER files at opening" is OFF).
	 * The guard flips the SER state for the duration of this function
	 * and restores on every return path. See SerAnalysisDebayerGuard
	 * above for the rationale. */
	SerAnalysisDebayerGuard ser_guard;
	if (maybe_engage_ser_debayer_for_analysis(ser_guard, regargs->seq)) {
		siril_log_warning(
		    _("mpp: auto-debayering CFA SER for analysis "
		      "(Preferences → Debayering → \"Debayer SER files at opening\" is OFF). "
		      "The raw mosaic is still available to the bayer-* drizzle stack path.\n"));
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

	if (run_owned) mpp_run_free(run);
	return MPP_OK;
}
