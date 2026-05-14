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
typedef struct mpp_aps mpp_aps_t;
typedef struct mpp_shifts mpp_shifts_t;

typedef enum {
	MPP_OK = 0,
	MPP_ENOTIMPL = -1,
	MPP_EINVAL = -2,
	MPP_ENOMEM = -3,
	MPP_EIO = -4,
	MPP_ENODATA = -5
} mpp_status_t;

/* mpp_config.h includes us back, so we can't include it here — would
 * create a parse-order cycle. The `cfg` field of mpp_run is therefore a
 * pointer; mpp_run_alloc allocates it and mpp_run_free releases it. */
typedef struct mpp_config mpp_config_t;

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
	int *global_shifts;        /* (dy, dx) per frame, contiguous, length 2 * num_frames */
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

	/* Stage A — per-AP frame ranking. */
	int stack_size;
	int *best_frame_indices;   /* aps->count × stack_size, flat row-major */

	/* Stage B — per-AP per-frame shifts. Null until Stage B runs. */
	mpp_shifts_t *shifts;
} mpp_run_t;

/* Allocate / free / shape helpers. */
mpp_run_t *mpp_run_alloc(void);
void mpp_run_free(mpp_run_t *run);

/* Stage A: read frames from `seq`, rank, globally align, build mean reference
 * frame, auto-place APs, compute per-AP frame qualities. Output `*run_out`
 * (caller frees via mpp_run_free). Honours seq->imgparam[i].incl if
 * `cfg->frames_normalization == true`? No — incl filtering is a separate
 * step; this Stage A processes every frame so the GUI's frame selector has
 * data for all of them. */
mpp_status_t mpp_analyze(sequence *seq, const mpp_config_t *cfg,
                         mpp_run_t **run_out);

/* Stage B: compute per-AP per-frame shifts (with PSS's penalty weight
 * matrix). Reads `run->aps`, `run->best_frame_indices`, etc.; writes
 * `run->shifts`. Re-reads frames from `seq`. */
mpp_status_t mpp_compute_shifts(sequence *seq, const mpp_config_t *cfg,
                                mpp_run_t *run);

/* Stage C: apply pre-computed Stage-B shifts and produce the final stacked
 * fits. Re-reads frames from `seq`. */
mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                             const mpp_run_t *run, fits *out);

/* Top-level Siril `register -method=mpp` entry. Drives Stages A + B and
 * writes the sidecar next to the sequence. */
int register_mpp(struct registration_args *regargs);

/* Publish per-frame Stage-A quality scores into seq->regparam[layer][i].quality
 * so Siril's existing frame selector and quality plot surface them. */
void mpp_write_quality_to_regdata(sequence *seq, int layer, const mpp_run_t *run);

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

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_H_ */
