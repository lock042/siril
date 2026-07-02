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

/* Per-AP patch bounds extended one grid step beyond every edge that no
 * other AP's patch covers (ORIGINAL frame coordinates), so Stage C
 * coverage — and with it the hand-off to the background blend — moves
 * off the faint structure hugging an object's rim into empty sky. */
std::vector<cv::Vec4i> stack_extended_patch_bounds(const mpp_aps_t &aps,
                                                   int dim_y, int dim_x,
                                                   const mpp_config_t &cfg);

/* Signal ramp clamp((v − lo)/lo, 0, 1) on a CV_32F luminance, and the
 * background line `lo` it is anchored to (AP brightness threshold in
 * pixel units; ≤ 0 disables the ramped blend boost). */
cv::Mat dc_signal_ramp(const cv::Mat &v_gray, double lo);
double dc_signal_floor(const mpp_config_t &cfg);

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
                                                const std::vector<int> &included = {},
                                                int max_threads = 1,
                                                bool provider_thread_safe = false);

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
                                            const int *included = nullptr,
                                            int max_threads = 1,
                                            bool provider_thread_safe = false);

/* Final-output conversion: scale the float
 * canvas by 1/255 (8-bit) or 1/65535, clip to [0,1] and cast to uint16 via
 * skimage img_as_uint rounding. */
cv::Mat stack_float_to_uint16(const cv::Mat &buf, int num_layers, int bitdepth);

/* Optional per-frame derotation hook for the warp engine. Fills, in the
 * engine's drizzled canvas (rows=DY, cols=DX):
 *   mapx/mapy (CV_32F) — source position in the MPP-aligned frame canvas of the
 *                        surface point seen at each epoch-canvas pixel (replaces
 *                        the identity base of the remap);
 *   mu (CV_32F)        — 1 where the map is usable, 0 only where a globe point
 *                        rotated behind the limb; folded into the accumulation
 *                        weight to drop those masked points. Pixels outside the
 *                        globe pass through (identity, mu=1) so rings/sky stack
 *                        normally.
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
