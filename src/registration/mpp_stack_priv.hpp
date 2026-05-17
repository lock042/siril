/*
 * Private C++ entry points for mpp_stack. Tests drive these directly using
 * cv::Mat inputs (no Siril sequence required).
 *
 * Phase 5a is broken into small, separately-testable pieces so divergences
 * from PSS can be caught early. Each piece has an oracle equivalence test
 * against artifacts dumped by tools/pss_reference/run_pss.py.
 */
#ifndef SRC_REGISTRATION_MPP_STACK_PRIV_HPP_
#define SRC_REGISTRATION_MPP_STACK_PRIV_HPP_

#include <opencv2/core.hpp>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_shift.h"
#include "registration/mpp_shift_priv.hpp"   /* for FrameOffset */
#include "registration/mpp_align_priv.hpp"   /* for progress_cb_fn */

namespace mpp {

/* PSS one_dim_weight (stack_frames.py:658). Length = patch_high - patch_low.
 * The peak (1.0) lands at index (box_center - patch_low). Ramp from
 * 1/(c+1) at the lower edge up to 1.0 at the centre, then 1.0 → 1/N at the
 * upper edge, where N = patch_high - box_center. extend_low / extend_high
 * replace the corresponding ramp with all-ones (edge-of-frame patches). */
std::vector<float> stack_one_dim_weight(int patch_low, int patch_high, int box_center,
                                        bool extend_low, bool extend_high);

/* PSS weight_matrix_first_phase (stack_frames.py:311-317). Shape =
 * (2*sw1+1, 2*sw1+1) where sw1 = (alignment_points_search_width - 4) / 2.
 * Centre is 1.0; off-centre values get a quadratic penalty. */
cv::Mat stack_build_first_phase_weight_matrix(const mpp_config_t &cfg);

/* PSS remap_rigid (stack_frames.py:508). Adds frame[y_low+shift_y:y_high+shift_y,
 * x_low+shift_x:x_high+shift_x] into buffer[..]. Out-of-frame portions are
 * clipped and the four border_* counters are updated to track the maximum
 * cropped extent in any frame. Both frame and buffer are CV_32F (single-
 * channel for mono). */
struct RemapBorder {
	int y_low = 0;
	int y_high = 0;
	int x_low = 0;
	int x_high = 0;
};
void stack_remap_rigid(const cv::Mat &frame_f32, cv::Mat &buffer_f32,
                       int shift_y, int shift_x,
                       int y_low, int y_high, int x_low, int x_high,
                       RemapBorder &border);

/* PSS alignment_points.compute_frame_qualities for the default Laplace
 * (frame-rank) + Laplace (AP-rank) path. The strided LoG of each frame is
 * computed once via rank_blurred_laplacian_u8 (shared with mpp_rank), then
 * each AP's patch is sliced from it and meanStdDev gives σ. With
 * frames_normalization on, σ is divided by frame avg_brightness. */
struct APQualities {
	int stack_size = 0;
	/* [num_aps][num_frames] — raw σ (or σ/brightness when normalised).
	 * Diagnostic; tests compare against PSS. */
	std::vector<std::vector<double>> qualities;
	/* [num_aps][stack_size] — top-N frame indices by descending quality. */
	std::vector<std::vector<int>> best_frame_indices;
	/* [num_frames][variable] — APs for which this frame is among the best. */
	std::vector<std::vector<int>> used_alignment_points;
};

/* `frames` are the raw mono frames (NOT blurred; we'll Gaussian-blur and
 * compute the strided LoG ourselves to match PSS). `frame_brightness[i]` is
 * the PSS frames.average_brightness for frame i (pre-computed via
 * mpp_rank::rank_average_brightness). `frame_rows`/`frame_cols` are the
 * common original-frame shape (used for the PSS clamp before stride
 * division). */
APQualities ap_compute_frame_qualities(const std::vector<cv::Mat> &frames,
                                       const std::vector<double> &frame_brightness,
                                       const mpp_aps_t &aps,
                                       const std::vector<FrameOffset> &offsets,
                                       int frame_rows, int frame_cols,
                                       const mpp_config_t &cfg,
                                       progress_cb_fn progress = nullptr,
                                       void *progress_user = nullptr);

/* Streamed overload — the provider returns the raw mono frame for the
 * requested index. Each iteration reads one frame, computes its strided
 * LoG once, scores every AP against it, then discards. Sequential by
 * construction (parallel reads would multiply the in-flight working
 * set by thread count). */
APQualities ap_compute_frame_qualities_streamed(const FrameProvider &provider,
                                                int num_frames,
                                                const std::vector<double> &frame_brightness,
                                                const mpp_aps_t &aps,
                                                const std::vector<FrameOffset> &offsets,
                                                int frame_rows, int frame_cols,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress = nullptr,
                                                void *progress_user = nullptr);

/* PSS prepare_for_stack_blending (stack_frames.py:154-213).
 *
 * Per-AP: drizzled patch bounds + drizzled centre + 2D weights_yx
 * (`min(weight_y[:, None], weight_x[None, :])`, ramping from the AP centre
 * to the patch boundaries, with the boundary ramp suppressed at frame
 * edges). Global: a sum_single_frame_weights buffer in drizzled coords,
 * accumulating `stack_size × weights_yx` over every AP. The count of pixels
 * with summed weight < 1e-10 is the "background hole" count.
 *
 * Initialised stacking_buffers (per-AP zeros, drizzled patch size) live
 * here so the main loop can accumulate into them. Mono only for Phase 5a;
 * colour comes later. */
/* A pre-drizzle rectangular region of the intersection that contains at
 * least one "hole" pixel needing background fill. See
 * stack_frames.py:240-277. Empty list ⇒ full-frame background. */
struct StackBackgroundPatch {
	int y_low, y_high, x_low, x_high;
};

struct StackState {
	int drizzle_factor = 1;
	int stack_size = 0;
	int num_layers = 1;       /* 1 = mono, 3 = RGB */
	int dim_y = 0;            /* pre-drizzle intersection size (for background) */
	int dim_x = 0;
	int dim_y_drizzled = 0;
	int dim_x_drizzled = 0;

	/* Per-AP — all length = aps->count. */
	std::vector<cv::Vec4i> patch_drizzled;  /* (y_low, y_high, x_low, x_high) */
	std::vector<cv::Vec2i> ap_drizzled;     /* (y, x) */
	std::vector<cv::Mat> weights_yx;        /* CV_32F, patch-sized */
	std::vector<cv::Mat> stacking_buffers;  /* CV_32F, patch-sized, zero-init */

	cv::Mat sum_single_frame_weights;       /* CV_32F, dim_y_drizzled × dim_x_drizzled */
	int number_stacking_holes = 0;

	/* Background plumbing (only meaningful when number_stacking_holes > 0). */
	std::vector<StackBackgroundPatch> background_patches;  /* empty ⇒ full frame */
	cv::Mat averaged_background;            /* CV_32F, dim_y × dim_x; empty if no holes */
};

StackState stack_prepare_for_blending(const mpp_aps_t &aps,
                                      const cv::Vec4i &intersection,
                                      int stack_size,
                                      int drizzle_factor,
                                      int num_layers,
                                      const mpp_config_t &cfg);

/* Stage B: per-AP per-frame shift compute (PSS stack_frames.py's per-(frame, AP)
 * call to compute_shift_alignment_point with weight_matrix_first_phase set).
 *
 * Stage B is split out from the main stacking loop so the GUI workflow can
 * pause after AP placement (Stage A) and frame selection but BEFORE the
 * heavy per-AP correlation. The resulting mpp_shifts_t is the canonical
 * Stage-B output: sparse where (frame, AP) is not in apq.used_alignment_points;
 * dense in the integer-mode default. Caller frees via mpp_shift_free. */
mpp_shifts_t *stack_compute_shifts(const std::vector<cv::Mat> &frames_mono_blurred,
                                   const cv::Mat &mean_frame_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const std::vector<FrameOffset> &offsets,
                                   const mpp_config_t &cfg,
                                   const int *included = nullptr);

/* Streamed overload: identical algorithm, but the caller supplies a
 * provider lambda that returns the blurred mono frame for the requested
 * index. Used by mpp_compute_shifts to read+blur each frame in the
 * correlation loop and immediately discard it, keeping memory to a
 * single frame instead of N. The thin overload above is implemented in
 * terms of this one. */
mpp_shifts_t *stack_compute_shifts_streamed(const FrameProvider &provider,
                                            int num_frames,
                                            const cv::Mat &mean_frame_raw,
                                            const mpp_aps_t &aps,
                                            const APQualities &apq,
                                            const std::vector<FrameOffset> &offsets,
                                            const mpp_config_t &cfg,
                                            const int *included = nullptr);

/* Stage C: apply pre-computed Stage-B shifts to produce per-AP buffers and
 * averaged_background. Takes the same shifts vector as Stage B; the rest of
 * the bookkeeping (brightness equalise, drizzle resize, remap_rigid,
 * background accumulate) is identical to stack_frames_loop. `included` is an
 * optional per-frame mask (length frames_raw.size()); when non-null, frames
 * with included[f]==0 are skipped entirely from the stack. NULL ⇒ all frames
 * participate (preserves the Phase-5a unit-test API). */
struct StackLoopOutput;  /* defined below */
StackLoopOutput stack_apply_shifts(const std::vector<cv::Mat> &frames_raw,
                                   const mpp_aps_t &aps,
                                   const APQualities &apq,
                                   const mpp_shifts_t *shifts,
                                   const std::vector<FrameOffset> &offsets,
                                   const std::vector<double> &frame_brightness,
                                   const std::vector<int> &quality_sorted_idx,
                                   const cv::Vec4i &intersection,
                                   const mpp_config_t &cfg,
                                   const int *included = nullptr);

/* Streamed overload — `num_frames` and `num_layers` are passed
 * explicitly because there's no upfront vector to introspect. The
 * provider returns the raw frame for each index; iteration is
 * sequential by construction (the per-AP and background accumulators
 * have shared state so parallelising would need per-thread buffers).
 * Cached overload above is a thin wrapper. */
StackLoopOutput stack_apply_shifts_streamed(const FrameProvider &provider,
                                            int num_frames,
                                            int num_layers,
                                            const mpp_aps_t &aps,
                                            const APQualities &apq,
                                            const mpp_shifts_t *shifts,
                                            const std::vector<FrameOffset> &offsets,
                                            const std::vector<double> &frame_brightness,
                                            const std::vector<int> &quality_sorted_idx,
                                            const cv::Vec4i &intersection,
                                            const mpp_config_t &cfg,
                                            const int *included = nullptr);

/* PSS stack_frames.stack_frames main loop (stack_frames.py:281-506).
 *
 * Inputs:
 *   frames_raw            — raw mono (uint16-stored in Siril). Brightness-
 *                           equalised and resampled to drizzled coords
 *                           inside.
 *   frames_mono_blurred   — Gaussian-blurred mono frames, fed into
 *                           multilevel_correlation for per-AP shift.
 *   mean_frame_raw        — unblurred mean_frame from align_average_frame
 *                           (NOT blurred — PSS's
 *                           set_reference_boxes_correlation samples from
 *                           align_frames.mean_frame, not from the
 *                           AlignmentPoints-blurred copy).
 *   apq                   — output of ap_compute_frame_qualities.
 *   offsets               — output of shift_frame_offsets (PSS
 *                           align_frames.dy/dx).
 *   frame_brightness      — per-frame avg brightness from
 *                           mpp::rank_average_brightness.
 *   quality_sorted_idx    — frames sorted by descending GLOBAL Phase-1
 *                           quality (the top stack_size of these contribute
 *                           to averaged_background).
 *   intersection          — intersection bounds in original-frame coords.
 *   cfg                   — Phase 5a config.
 *
 * Output: filled StackState (stacking_buffers, averaged_background) plus
 * border counts and shift_failure_counter. */
struct StackLoopOutput {
	StackState state;
	RemapBorder border;
	int shift_failure_counter = 0;
};

StackLoopOutput stack_frames_loop(const std::vector<cv::Mat> &frames_raw,
                                  const std::vector<cv::Mat> &frames_mono_blurred,
                                  const cv::Mat &mean_frame_raw,
                                  const mpp_aps_t &aps,
                                  const APQualities &apq,
                                  const std::vector<FrameOffset> &offsets,
                                  const std::vector<double> &frame_brightness,
                                  const std::vector<int> &quality_sorted_idx,
                                  const cv::Vec4i &intersection,
                                  const mpp_config_t &cfg);

/* PSS merge_alignment_point_buffers (stack_frames.py:561).
 *
 * 1. For each AP: add `stacking_buffer * weights_yx` into the global
 *    stacked_image_buffer at the AP's drizzled patch bounds.
 * 2. Divide stacked_image_buffer by sum_single_frame_weights (per-pixel).
 * 3. If there are background holes, blend in averaged_background using a
 *    clipped `sum_weights / (blend_threshold × stack_size)` foreground
 *    weight (PSS line 612-624).
 * 4. Trim the four border counters off the edges of the result.
 * 5. Scale by 1/255 (for 8-bit) or 1/65535 (for 16-bit), clip to [0, 1],
 *    cast to uint16 via skimage `img_as_uint` (which is equivalent to
 *    `(x * 65535).round().clip(0, 65535).astype(uint16)` for our case).
 *
 * Returns the final uint16 (CV_16U) stacked image. */
cv::Mat stack_merge_alignment_point_buffers(const StackState &state,
                                            const RemapBorder &border,
                                            const mpp_aps_t &aps,
                                            const mpp_config_t &cfg);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_STACK_PRIV_HPP_ */
