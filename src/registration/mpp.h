#ifndef SRC_REGISTRATION_MPP_H_
#define SRC_REGISTRATION_MPP_H_

#include <stddef.h>
#include <stdint.h>

/* Pull in Siril's `sequence` and `fits` typedefs so the orchestrator entry
 * points below can refer to them by their proper typedef names. Siril
 * tags these as `struct sequ` / `struct ffit`; using `struct sequence`
 * would create a different (incompatible) type. */
#include "core/siril.h"

#ifdef __cplusplus
extern "C" {
#endif

struct registration_args;
struct mpp_derot;
typedef struct mpp_aps mpp_aps_t;
typedef struct mpp_shifts mpp_shifts_t;

typedef enum {
	MPP_OK = 0,
	MPP_ENOTIMPL = -1,
	MPP_EINVAL = -2,
	MPP_ENOMEM = -3,
	MPP_EIO = -4,
	MPP_ENODATA = -5,
	MPP_EINTR = -6      /* user-requested cancellation */
} mpp_status_t;

/* mpp_config.h includes us back, so we can't include it here — would
 * create a parse-order cycle. The `cfg` field of mpp_run is therefore a
 * pointer; mpp_run_alloc allocates it and mpp_run_free releases it. */
typedef struct mpp_config mpp_config_t;

/* Opaque holder for the analysis-frame cache populated by Stage A
 * (mono cv::Mat per frame, populated for the K top-by-quality frames
 * the memory picker found room for). When non-NULL Stage B and the
 * mpp_recompute_qualities helper hit it directly instead of re-reading
 * from disk; cache misses fall through to read_analysis_frame. Stage C
 * doesn't share layout with the analysis cache so it doesn't consult
 * it. Not serialised to the sidecar — purely an in-memory perf
 * accelerator owned by the current run. */
typedef struct mpp_cache mpp_cache_t;

void mpp_cache_free(mpp_cache_t *cache);

/* mpp_run_t — in-memory state spanning the orchestrator stages.
 * Stage A populates everything up to and including `aps`/`stack_size`/
 * `best_frame_indices`. Stage B fills `shifts`. Stage C consumes the whole
 * thing and produces an output fits. The sidecar serialises every field
 * here for the register/stack split workflow.
 *
 * All pointer fields are heap-owned; free with mpp_run_free(). */
typedef struct mpp_run {
	/* Cfg snapshot used by the run (heap-allocated; owned by run). */
	mpp_config_t *cfg;

	/* Sequence properties at analyze-time. */
	int num_frames;
	int frame_rows;
	int frame_cols;
	int num_layers;        /* 1 = mono, 3 = RGB */
	int bitdepth;          /* 8 or 16 (PSS convention; see mpp_bitdepth_from_fits_bitpix) */

	/* Stage A — per-frame outputs (each length num_frames, malloc'd). */
	double *quality;
	double *frame_brightness;
	int *included;             /* 0 or 1; 1 by default. CLI -selected filters here. */
	/* (dy, dx) per frame, contiguous, length 2 * num_frames. Sub-pixel
	 * (parabolic-fit refinement on the global correlation peak): with
	 * drizzle_scale > 1 the fractional component is the dominant source
	 * of cross-frame CFA-phase coverage in Bayer drizzle, so integer-only
	 * values would leave the same CFA phase visited by every frame and
	 * produce vertical-stripe / partial-Bayer-grid artefacts at the
	 * output. v3-and-earlier sidecars stored these as int32 — readers
	 * keyed on version promote on the fly. */
	double *global_shifts;
	int best_frame_idx;

	/* Stage A — alignment + intersection. */
	int patch_yxyx[4];         /* (y_low, y_high, x_low, x_high) on best frame */
	int intersection[4];       /* (y_low, y_high, x_low, x_high) — mean-frame bounds in orig coords */

	/* Stage A — average reference frame (CV_32S equivalent layout).
	 * Always 16-bit-equivalent range regardless of bitdepth. */
	int32_t *mean_frame_data;
	int mean_frame_rows;
	int mean_frame_cols;

	/* Stage A — AP grid (mpp_aps_t owns its records). */
	mpp_aps_t *aps;

	/* Per-AP frame ranking. Computed at Stage B entry by
	 * mpp_recompute_qualities, NOT by Stage A — Analyze leaves
	 * best_frame_indices null (selected_per_ap 0, taper 0) and sets only a
	 * provisional stack_size from cfg, so the AP grid is editable without a
	 * ranking pass that an edit would discard. `selected_per_ap` is the
	 * baked row stride of best_frame_indices (stack_size + taper at ranking
	 * time); it never changes after bake, while stack_size may later be
	 * shrunk by a Stack-tab override (entry weights are re-derived from
	 * rank — see apq_from_run). `taper` is the soft-selection taper
	 * half-width in ranks; 0 = hard top-N. */
	int stack_size;
	int taper;
	int selected_per_ap;
	int *best_frame_indices;   /* aps->count × selected_per_ap, flat row-major; null until ranking */

	/* Stage B — per-AP per-frame shifts. Null until Stage B runs. */
	mpp_shifts_t *shifts;

	/* Fingerprint (mpp_derot_fingerprint) of the .derot plan this run was
	 * registered against when derotation was active (0 = registered without
	 * derotation). Persisted so Stage C can refuse to apply a .derot that
	 * differs in ANY way (epoch, disk fit, body/system, per-frame geometry)
	 * from the one the AP shifts were measured in — applying derotation on top
	 * of mismatched shifts double-counts the planet's rotation. */
	uint64_t derot_fingerprint;

	/* Analysis-frame cache populated by Stage A and reused by Stage B
	 * and mpp_recompute_qualities. Owned by the run; freed via
	 * mpp_cache_free(). Null when caching was skipped (STREAMING) or
	 * when the run was loaded from a sidecar (sidecar doesn't carry
	 * the cache; subsequent stages re-read from disk). */
	mpp_cache_t *cache;
} mpp_run_t;

/* Allocate / free / shape helpers. */
mpp_run_t *mpp_run_alloc(void);
void mpp_run_free(mpp_run_t *run);

/* Release the run's analysis cache and null the field. No-op when the
 * run is NULL or already cache-less. Called before Stage C (which can't
 * consume the mono single-channel cache layout — read_full_frame returns
 * full-depth multi-channel) to free budget the streaming-only Stage C
 * memory picker doesn't otherwise see; subsequent AP re-edits repopulate
 * the cache via mpp_recompute_qualities. */
void mpp_run_drop_cache(mpp_run_t *run);

/* Stage A: read frames from `seq`, rank, globally align, build mean reference
 * frame, auto-place APs. Per-AP frame ranking is NOT done here — it is
 * deferred to Stage B entry (mpp_recompute_qualities) so the AP grid can be
 * reviewed/edited before paying that pass. Output `*run_out`
 * (caller frees via mpp_run_free). Honours seq->imgparam[i].incl if
 * `cfg->frames_normalization == true`? No — incl filtering is a separate
 * step; this Stage A processes every frame so the GUI's frame selector has
 * data for all of them. */
mpp_status_t mpp_analyze(sequence *seq, const mpp_config_t *cfg,
                         mpp_run_t **run_out, const struct mpp_derot *derot);

/* Stage B: compute per-AP per-frame shifts (with PSS's penalty weight
 * matrix). Reads `run->aps`, `run->best_frame_indices`, etc.; writes
 * `run->shifts`. Re-reads frames from `seq`. */
mpp_status_t mpp_compute_shifts(sequence *seq, const mpp_config_t *cfg,
                                mpp_run_t *run, const struct mpp_derot *derot);

/* Stage C: apply pre-computed Stage-B shifts and produce the final stacked
 * fits. Re-reads frames from `seq`. `run` is non-const because Stage C
 * drops the analysis cache at entry — see mpp_run_drop_cache.
 *
 * `derot` (may be NULL) is the optional derotation plan from a .derot sidecar;
 * when non-NULL and cfg->stack_method == MPP_STACK_WARP, each frame is
 * derotated to the reference epoch as part of the warp engine's single
 * resample. */
mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                             mpp_run_t *run, const struct mpp_derot *derot,
                             fits *out);

/* Multi-source derotation stack (Option B). Combine `n` sequences — each with a
 * derotation plan `derots[i]` (built/loaded from <seqname>.derot) referencing a
 * common epoch — into one stack in `seqs[ref]`'s epoch canvas, written to `out`.
 * Every frame of every sequence is relocated + derotated into the reference
 * canvas and stacked as one long sequence. `cfg->stack_method` is forced to the
 * warp engine internally. The sequences must stay open for the call. */
mpp_status_t mpp_multistack(sequence **seqs, struct mpp_derot **derots,
                            int n, int ref, const mpp_config_t *cfg, fits *out);

/* Top-level Siril `register -method=mpp` entry. Drives Stages A + B and
 * writes the sidecar next to the sequence. */
int register_mpp(struct registration_args *regargs);

/* Classify a sequence's input layout for drizzle-mode dispatch. Same
 * answer regardless of the user's debayer-on-open preference, because
 * the classifier reads the SER ColorID directly (or falls back to
 * nb_layers for non-SER inputs).
 *
 *   MONO: true single-channel (SER_MONO, mono FITS)        — stsci-* OK
 *   CFA:  Bayer pattern, raw or debayered                  — bayer-* OK
 *   RGB:  true 3-channel non-Bayer (RGB SER, debayered FITS without bayer_pattern keyword)
 *                                                            — only bicubic/Off
 *
 * Used by:
 *   - process_pss / process_stack_mpp to reject incompatible -drizzle modes early
 *   - the GUI drizzle-mode combo to populate only valid choices
 *   - register_mpp to decide whether to auto-debayer for analysis */
typedef enum {
	MPP_INPUT_MONO = 0,
	MPP_INPUT_CFA  = 1,
	MPP_INPUT_RGB  = 2,
} mpp_input_type;

mpp_input_type mpp_classify_sequence_input(const sequence *seq);

/* Publish per-frame Stage-A quality scores into seq->regparam[layer][i].quality
 * so Siril's existing frame selector and quality plot surface them. */
void mpp_write_quality_to_regdata(sequence *seq, int layer, const mpp_run_t *run);

/* AP editor support. Mutate the cached run's AP grid in place. Per-AP
 * frame qualities (run->best_frame_indices) are already null coming out of
 * Analyze (ranking is deferred to Stage B), and every edit re-clears them
 * defensively; Stage B (mpp_recompute_qualities) ranks the final grid.
 *
 * mpp_ap_replace: re-runs auto-placement on the cached mean frame with
 *   the supplied cfg (typically read from the editor's spinners).
 *   Replaces run->aps. Updates run->cfg. Returns MPP_ENODATA if zero APs
 *   would survive the threshold filters with the new cfg.
 * mpp_ap_clear_all: drops every AP. Useful pre-cursor to manual placement.
 *
 * Mouse-driven editing (slice 10):
 * mpp_ap_add: append an AP centred at (x, y). Box/patch bounds are
 *   derived from run->cfg->alignment_points_half_box_width. structure
 *   is set to 1.0 (user-added APs bypass the structure threshold).
 * mpp_ap_remove: drop AP at index i.
 * mpp_ap_move: relocate AP at index i to (x, y).
 * mpp_ap_hit_test: return the topmost (last) AP whose box contains
 *   (x, y), or -1 if none.
 *
 * AP-snapshot helpers (slice 10) for Cancel/revert in the editor:
 * mpp_aps_snapshot: deep-copy the AP records (returns a new mpp_aps_t).
 * mpp_aps_restore: replace run->aps with the snapshot (transfers
 *   ownership; the snapshot pointer must not be freed by the caller
 *   after this call).
 *
 * All are no-ops if run is NULL (return MPP_EINVAL for the int ones). */
mpp_status_t mpp_ap_replace(mpp_run_t *run, const mpp_config_t *cfg);
void         mpp_ap_clear_all(mpp_run_t *run);
mpp_status_t mpp_ap_add(mpp_run_t *run, int x, int y);
mpp_status_t mpp_ap_remove(mpp_run_t *run, int i);
mpp_status_t mpp_ap_move(mpp_run_t *run, int i, int x, int y);
/* Grow (step>0) or shrink (step<0) AP `i`'s box/patch by `step` px of
 * half-box, clamped so the box stays inside the reference frame. GUI-only. */
mpp_status_t mpp_ap_resize(mpp_run_t *run, int i, int step);
/* Current half-box of AP `i` (max centre-to-edge), or -1 on error. */
int          mpp_ap_get_half_box(const mpp_run_t *run, int i);
/* Set AP `i`'s half-box absolutely (clamped); used to revert a resize. */
mpp_status_t mpp_ap_set_half_box(mpp_run_t *run, int i, int hb);
int          mpp_ap_hit_test(const mpp_run_t *run, int x, int y);

mpp_aps_t *mpp_aps_snapshot(const mpp_aps_t *src);
void       mpp_aps_restore(mpp_run_t *run, mpp_aps_t *snapshot);

/* Compute per-AP frame qualities for the current run->aps, populating
 * run->best_frame_indices (left null by Analyze, and cleared by any AP
 * edit). Re-reads analysis frames from `seq`; uses cached
 * run->frame_brightness and run->global_shifts (no re-rank, no re-align).
 * Idempotent: returns MPP_OK immediately if best_frame_indices is already
 * populated (e.g. a sidecar-loaded run).
 *
 * Called at Stage B start on every Register path (CLI and GUI). This is
 * where per-AP ranking actually happens — deferred out of Stage A so the
 * final (possibly hand-edited) AP grid is ranked exactly once. */
mpp_status_t mpp_recompute_qualities(sequence *seq, mpp_run_t *run,
                                     const struct mpp_derot *derot);

/* GUI-side cache of the current run. Stage A installs the run via
 * mpp_set_cached_run; the AP overlay and the AP editor read it back via
 * mpp_get_cached_run; close_sequence calls mpp_clear_cached_run.
 *
 * Ownership transfers to the cache on set: the cache frees the prior run
 * (if any), and frees the current run on clear or replacement. Thread-safe
 * via an internal mutex; callers must not hold a returned pointer across
 * a call to mpp_set/clear_cached_run from another thread.
 *
 * The actual storage is com.mpp_run (declared as void* in core/siril.h to
 * avoid include cycles); these helpers cast for you. */
void mpp_set_cached_run(mpp_run_t *run);
void mpp_clear_cached_run(void);
mpp_run_t *mpp_get_cached_run(void);

/* Convert a click in display (overlay) coords — top-left origin, matching
 * what cairo_matrix_transform_point(image_matrix, ...) yields — to the
 * mean_frame pdata-row coord the AP records use. Inverse of the draw-time
 * transform in draw_mpp_aps: applies the Y flip, and when viewing a
 * sequence frame (current_frame_idx in [0, num_frames)) also subtracts
 * the per-frame global shift so a click on a feature on frame i is stored
 * as that feature's position in the mean_frame. */
void mpp_display_to_ap_coord(const mpp_run_t *run, int gfit_rx, int gfit_ry,
                             int current_frame_idx,
                             int click_x, int click_y,
                             int *out_x, int *out_y);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_H_ */
