/*
 * Private C++ entry points for mpp_stack. Tests drive these directly using
 * cv::Mat inputs (no Siril sequence required).
 *
 */
#ifndef SRC_REGISTRATION_MPP_STACK_PRIV_HPP_
#define SRC_REGISTRATION_MPP_STACK_PRIV_HPP_

#include <opencv2/core.hpp>
#include <functional>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_shift.h"
#include "registration/mpp/mpp_shift_priv.hpp"   /* for FrameOffset */
#include "registration/mpp/mpp_align_priv.hpp"   /* for progress_cb_fn */

namespace mpp {

/* Length = patch_high - patch_low.
 * The peak (1.0) lands at index (box_center - patch_low). Raised-cosine
 * (Hann) taper: the linear ramp t (1/(c+1) … c/(c+1) below the centre,
 * N/N … 1/N above it, N = patch_high - box_center) is mapped through
 * sin²(π/2·t), giving zero slope at the centre and at the patch edges so
 * the blend window is C1 (PSS uses the C0 triangular ramp directly).
 * extend_low / extend_high replace the corresponding taper with all-ones
 * (edge-of-frame patches). */
std::vector<float> stack_one_dim_weight(int patch_low, int patch_high, int box_center,
                                        bool extend_low, bool extend_high);

/* Shape = (2*sw1+1, 2*sw1+1) where sw1 = (alignment_points_search_width - 4) / 2.
 * Centre is 1.0; off-centre values get a quadratic penalty. */
cv::Mat stack_build_first_phase_weight_matrix(const mpp_config_t &cfg);

/* Adds frame[y_low+shift_y:y_high+shift_y,
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
                       RemapBorder &border, float weight = 1.0f);

/* Fractional-shift variant: same accumulate-into-buffer contract, with the
 * fractional part of (shift_y, shift_x) resolved by Lanczos interpolation
 * (cv::warpAffine, translation-only). Shifts within 1e-3 of an integer take
 * the exact stack_remap_rigid blit. Destination rows/columns whose sample
 * position would fall outside the frame are clipped and recorded in
 * `border`, like the integer path. `weight` scales the contribution
 * (soft frame selection); 1.0 is a plain add. */
void stack_remap_subpixel(const cv::Mat &frame_f32, cv::Mat &buffer_f32,
                          double shift_y, double shift_x,
                          int y_low, int y_high, int x_low, int x_high,
                          RemapBorder &border, float weight = 1.0f);

/* Per-AP per-frame quality for the default Laplace (frame-rank) +
 * Laplace (AP-rank) path. The strided LoG of each frame is computed
 * once via rank_blurred_laplacian_u8 (shared with mpp_rank), then each
 * AP's patch is sliced from it and meanStdDev gives σ. With
 * frames_normalization on, σ is divided by frame avg_brightness. */
/* A (frame → AP) selection entry: the AP index plus the frame's selection
 * weight at that AP (1.0 on the plateau, raised-cosine in the taper zone —
 * see stack_selection_weight). */
struct APUse {
	int ap;
	float weight;
};

struct APQualities {
	int stack_size = 0;
	/* Soft-selection taper half-width in ranks. 0 = hard top-N (PSS
	 * behaviour; always the case in explicit frame-number mode). When
	 * > 0, each AP keeps stack_size + taper frames, with ranks
	 * stack_size − taper … stack_size + taper − 1 weighted by a raised
	 * cosine so adjacent APs that rank a borderline frame slightly
	 * differently give it nearly the same weight instead of a 1-vs-0
	 * cliff. The cosine is symmetric around rank stack_size, so the
	 * effective frame count Σw stays exactly stack_size. */
	int taper = 0;
	/* [num_aps][num_frames] — raw σ (or σ/brightness when normalised).
	 * Diagnostic; tests compare against the reference implementation. */
	std::vector<std::vector<double>> qualities;
	/* [num_aps][stack_size + taper] — top frame indices by descending
	 * quality. The selection weight of entry k is
	 * stack_selection_weight(k, stack_size, taper). */
	std::vector<std::vector<int>> best_frame_indices;
	/* [num_frames][variable] — (AP, weight) pairs for which this frame
	 * carries selection weight. */
	std::vector<std::vector<APUse>> used_alignment_points;
};

/* Selection weight of the frame at 0-based quality rank `rank` within an
 * AP, for a target stack size N and taper half-width T:
 *   1.0 for rank < N − T, 0.0 for rank ≥ N + T, raised-cosine in between
 *   (half-sample-centred so Σ over all ranks is exactly N).
 * T = 0 reproduces the hard top-N cliff. */
float stack_selection_weight(int rank, int stack_size, int taper);

/* `frames` are the raw mono frames (NOT blurred; we'll Gaussian-blur and
 * compute the strided LoG ourselves). `frame_brightness[i]` is the
 * average brightness for frame i (pre-computed via
 * mpp_rank::rank_average_brightness). `frame_rows`/`frame_cols` are the
 * common original-frame shape (used for the frame-bounds clamp before
 * stride division). */
/* `included` (empty = all): when sized to num_frames, excluded frames are
 * given a sentinel quality so they never enter any AP's best_frame_indices,
 * and the per-AP selection count is clamped to the included-frame count. */
APQualities ap_compute_frame_qualities(const std::vector<cv::Mat> &frames,
                                       const std::vector<double> &frame_brightness,
                                       const mpp_aps_t &aps,
                                       const std::vector<FrameOffset> &offsets,
                                       int frame_rows, int frame_cols,
                                       const mpp_config_t &cfg,
                                       progress_cb_fn progress = nullptr,
                                       void *progress_user = nullptr,
                                       const std::vector<int> &included = {});

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
                                                void *progress_user = nullptr,
                                                const std::vector<int> &included = {});

/*
 * Per-AP: drizzled patch bounds + drizzled centre + 2D weights_yx
 * (`weight_y[:, None] * weight_x[None, :]`, Hann-tapering from the AP
 * centre to the patch boundaries, with the boundary taper suppressed at
 * frame edges). Global: a sum_single_frame_weights buffer in drizzled coords,
 * accumulating `stack_size × weights_yx` over every AP. The count of pixels
 * with summed weight < 1e-10 is the "background hole" count.
 *
 * Initialised stacking_buffers (per-AP zeros, drizzled patch size) live
 * here so the main loop can accumulate into them.
 * A pre-drizzle rectangular region of the intersection that contains at
 * least one "hole" pixel needing background fill. Empty list ⇒ full-frame
 * background. */
struct StackBackgroundPatch {
	int y_low, y_high, x_low, x_high;
};

struct StackState {
	double drizzle_scale = 1.0;   /* output scale; fractional supported (cv::resize) */
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
	/* Frames accumulated per AP (stack_size, or the reduced effective
	 * count under stack_skip_failed_aps). Normalises a stacking buffer to
	 * a mean patch — used by the merge's per-AP DC equalisation. */
	std::vector<float> ap_frame_counts;

	cv::Mat sum_single_frame_weights;       /* CV_32F, dim_y_drizzled × dim_x_drizzled */
	int number_stacking_holes = 0;

	/* Background plumbing (only meaningful when number_stacking_holes > 0). */
	std::vector<StackBackgroundPatch> background_patches;  /* empty ⇒ full frame */
	cv::Mat averaged_background;            /* CV_32F, dim_y × dim_x; empty if no holes */
};

StackState stack_prepare_for_blending(const mpp_aps_t &aps,
                                      const cv::Vec4i &intersection,
                                      int stack_size,
                                      double drizzle_scale,
                                      int num_layers,
                                      const mpp_config_t &cfg,
                                      const std::vector<float> *ap_effective_counts = nullptr);

/* Stage B: per-AP per-frame shift compute.
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
                                   const int *included = nullptr,
                                   int max_threads = 0,
                                   size_t mem_budget_bytes = 0);

/* Streamed overload — `num_frames` and `num_layers` are passed
 * explicitly because there's no upfront vector to introspect. The
 * provider returns the raw frame for each index. Cached overload above
 * is a thin wrapper.
 *
 * Parallelism (bounded by `max_threads`): the frame loop is parallelised —
 * each thread runs the full per-frame pipeline (read + brightness +
 * resize + per-AP remap) for its share of frames into PRIVATE per-AP
 * buffers, which are then summed. This keeps every core busy on both the
 * (cheap, I/O-bound) decode and the heavier AP accumulation with no phase
 * barrier and no fixed I/O-vs-AP pool split. Only enabled when
 * `provider_thread_safe` (reentrant provider, e.g. SER / reentrant FITS);
 * otherwise the read — and hence the whole loop — stays single-threaded.
 *
 * The float sum is reordered by the per-thread reduction, so the output is
 * not bit-identical across thread counts (measured ≤1 16-bit LSB on an
 * 8000-frame stack) — the same trade the drizzle paths make, negligible for
 * a statistical stack.
 *
 * `mem_budget_bytes` (0 = unbounded) caps the thread count so the private
 * per-AP buffer sets fit in the available memory. */
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
                                            const int *included = nullptr,
                                            int max_threads = 0,
                                            bool provider_thread_safe = false,
                                            size_t mem_budget_bytes = 0);

/* stack_frames main loop.
 *
 * Inputs:
 *   frames_raw            — raw mono (uint16-stored in Siril). Brightness-
 *                           equalised and resampled to drizzled coords
 *                           inside.
 *   frames_mono_blurred   — Gaussian-blurred mono frames, fed into
 *                           multilevel_correlation for per-AP shift.
 *   mean_frame_raw        — unblurred mean_frame from align_average_frame
 *                           (NOT blurred — Stage B's per-AP reference
 *                           boxes are sampled from the unblurred mean,
 *                           not from the AP-placement blurred copy).
 *   apq                   — output of ap_compute_frame_qualities.
 *   offsets               — output of shift_frame_offsets.
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
	/* An allocation failure was contained inside the frame-parallel
	 * accumulation (the partial result must be discarded). Lets the
	 * caller report MPP_ENOMEM instead of crashing — on Windows there
	 * is no overcommit, so cv allocations genuinely fail under memory
	 * pressure. */
	bool oom = false;
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

/* merge_alignment_point_buffers.
 *
 * 1. For each AP: add `stacking_buffer * weights_yx` into the global
 *    stacked_image_buffer at the AP's drizzled patch bounds.
 * 2. Divide stacked_image_buffer by sum_single_frame_weights (per-pixel).
 * 3. If there are background holes, blend in averaged_background using a
 *    clipped `sum_weights / (blend_threshold × stack_size)` foreground
 *    weight.
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

/* Final-output conversion shared by both stacking engines: scale the float
 * canvas by 1/255 (8-bit) or 1/65535, clip to [0,1] and cast to uint16 via
 * skimage img_as_uint rounding. */
cv::Mat stack_float_to_uint16(const cv::Mat &buf, int num_layers, int bitdepth);

/* Optional per-frame derotation hook for the warp engine. Fills, in the
 * engine's drizzled canvas (rows=DY, cols=DX):
 *   mapx/mapy (CV_32F) — source position in the MPP-aligned frame canvas of the
 *                        surface point seen at each epoch-canvas pixel (replaces
 *                        the identity base of the remap);
 *   mu (CV_32F)        — cosine of the emission angle at frame time (0 off-disk),
 *                        folded into the accumulation weight so foreshortened
 *                        near-limb samples are down-weighted and off-disk points
 *                        drop out.
 * Returns false to treat the frame as un-derotated (identity base, mu = 1). The
 * provider must be thread-safe (called concurrently across frames). */
using DerotMapProvider =
    std::function<bool(int frame, int DY, int DX,
                       cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu)>;

/* Warp-field stacking engine (MPP_STACK_WARP) — Stage C alternative to the
 * patch-blend pipeline. Per frame: per-AP shifts -> smooth dense displacement
 * field D, per-(frame,AP) selection weights -> dense weight field W; the frame
 * is warped once (cv::remap, Lanczos4) through
 *   map = derot_base + global_offset*S - D*S
 * and accumulated acc += (W*mu) (.) warped, wsum += W*mu. Output acc/wsum plus
 * the PSS background blend where AP coverage thins. `derot` is optional; when
 * null the base is the identity and mu = 1 (pure de-warp, no rotation). Same
 * frame-parallel private-accumulator scheme and mem_budget clamp as
 * stack_apply_shifts_streamed. */
struct WarpStackResult {
	cv::Mat image;          /* final CV_16U / CV_16UC3; empty on cancel/error */
	bool oom = false;
	bool cancelled = false;
};

WarpStackResult stack_warp_apply_streamed(const FrameProvider &provider,
                                          int num_frames, int num_layers,
                                          const mpp_aps_t &aps,
                                          const APQualities &apq,
                                          const mpp_shifts_t *shifts,
                                          const std::vector<FrameOffset> &offsets,
                                          const std::vector<double> &frame_brightness,
                                          const std::vector<int> &quality_sorted_idx,
                                          const cv::Vec4i &intersection,
                                          const mpp_config_t &cfg,
                                          const int *included = nullptr,
                                          int max_threads = 0,
                                          bool provider_thread_safe = false,
                                          size_t mem_budget_bytes = 0,
                                          const DerotMapProvider *derot = nullptr);

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_STACK_PRIV_HPP_ */
