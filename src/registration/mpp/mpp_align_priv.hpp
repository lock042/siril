/*
 * Private C++ entry points for mpp_align. Tests drive these directly with
 * cv::Mat inputs (without going through a Siril sequence).
 */
#ifndef SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_
#define SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_

#include <opencv2/core.hpp>
#include <functional>
#include <vector>

#include "registration/mpp.h"
#include "registration/mpp/mpp_config.h"

struct sequ;   /* Siril sequence (typedef struct sequ sequence) */

namespace mpp {

/* Generic progress callback for long stages. fraction is in [0, 1] over
 * the stage's total work. Called periodically (typically once per
 * frame). Passing nullptr disables. */
typedef void (*progress_cb_fn)(double fraction, void *user);

/* Threshold-weighted quality measure. `frame` is the candidate patch
 * (best-frame mono_blurred cropped to (H, W)). Returns the quality
 * score (higher = more structure). Replicates NumPy uint16 modular
 * subtraction for u16 input; CV_32F input uses signed math. */
double align_quality_measure_threshold_weighted(const cv::Mat &frame,
                                                int stride,
                                                double black_threshold,
                                                double min_fraction);

/* Picks the (y_low, y_high, x_low, x_high) patch with the highest quality
 * over a coarse grid stepped across the best frame's blurred mono image. */
cv::Vec4i align_pick_patch(const cv::Mat &best_frame_mono_blurred,
                           const mpp_config_t &cfg);

/* Brightness-weighted centroid of a mono frame — PSS Planet mode
 * (AlignFrames.center_of_gravity, align_frames.py). The threshold/clip
 * pixel-selection matches PSS exactly: threshold = trunc((min+max)/2),
 * weights = max(frame - threshold, 0). The returned (y, x) is the
 * moment centroid m01/m00, m10/m00 kept at SUB-PIXEL precision —
 * PSS rounds these to integer, but the port keeps the fractional part
 * to seed drizzle's cross-frame CFA-phase diversity (same rationale as
 * the always-on sub-pixel global correlation residual). Oracle tests
 * must lround() before comparing to PSS's integer cog. Returns the
 * frame centre if the frame has no pixels above threshold (m00 == 0). */
cv::Vec2d center_of_gravity(const cv::Mat &frame_mono_blurred);

struct AlignShiftResult {
	/* Sign convention: shift the frame by (dy,dx) to align with ref.
	 * Sub-pixel (multilevel_correlation's parabolic refinement) is
	 * always enabled for the global pass — feeds Bayer drizzle's
	 * cross-frame CFA-phase diversity, see mpp_pixmap_build_filtered. */
	double dy = 0.0;
	double dx = 0.0;
	bool   success = false;
};

/* Sub-pixel variant; used by per-AP shift compute when subpixel_solve is on. */
struct MultilevelShiftResult {
	double dy = 0.0;     /* same sign convention; fractional when subpixel_solve=true */
	double dx = 0.0;
	bool success = false;
};

/* Generic two-phase multilevel correlation. `ref_full_f32` is the
 * full-resolution float32 reference patch; `ref_first_phase_f32` is its
 * stride-2 view. The frame is whatever dtype the source produces;
 * converted to f32 internally where matchTemplate needs it. `gauss_width`
 * is the kernel size of the extra blur applied to the frame's phase-1
 * window (default 7). `search_width` is the search half-width — total
 * search range is split between phases (phase 2 fixed at ±4, phase 1 gets
 * (search_width − 4)/2 on the coarse grid). If `subpixel_solve`, the
 * phase-2 peak is refined by fitting a quadratic over its 3×3
 * neighbourhood. `weight_matrix_first_phase`, when non-empty, is
 * element-wise multiplied into the phase-1 correlation result before
 * argmax — used during stacking to apply an off-centre penalty bias
 * (alignment_points_penalty_factor). */
MultilevelShiftResult multilevel_correlation(const cv::Mat &ref_full_f32,
                                             const cv::Mat &ref_first_phase_f32,
                                             const cv::Mat &frame_mono_blurred,
                                             int patch_y_low, int patch_y_high,
                                             int patch_x_low, int patch_x_high,
                                             int gauss_width, int search_width,
                                             bool subpixel_solve,
                                             const cv::Mat &weight_matrix_first_phase = cv::Mat());

/* Phase 2 wrapper. */
AlignShiftResult align_shift_one_frame(const cv::Mat &reference_window_f32,
                                       const cv::Mat &reference_window_first_phase_f32,
                                       const cv::Mat &frame_mono_blurred,
                                       int patch_y_low, int patch_y_high,
                                       int patch_x_low, int patch_x_high,
                                       const mpp_config_t &cfg);

struct AlignGlobalResult {
	std::vector<cv::Vec2d> shifts;  /* (dy, dx) per frame, sub-pixel */
	cv::Vec4i patch_yxyx;           /* (y_low, y_high, x_low, x_high) on best frame */
	int best_frame_idx = -1;
	/* An allocation failure was contained inside one of the parallel
	 * sweeps (best_frame_idx is also -1). Distinguishes out-of-memory
	 * from user cancellation so the caller can report MPP_ENOMEM. */
	bool oom = false;
};

/* Returns the cv::Mat for frame index `i`. The pixel type / channel
 * count is whatever the consumer requires (blurred mono for the
 * align/shift paths, raw mono for the average-frame and per-AP-quality
 * paths, raw multi-channel for the stack paths) — the typedef is
 * intentionally semantic-free.
 *
 * Used by every "_streamed" / "_from_provider" overload so the caller
 * can decide whether each frame is fetched from a cached vector
 * (cheap view), lazily derived from a raw cache (one blur per call),
 * or read fresh from disk. The Mat only needs to stay valid for the
 * single iteration it's used in. Empty Mat signals read failure — the
 * inner loops treat this the same as an excluded frame. */
using FrameProvider = std::function<cv::Mat(int idx)>;

/* Optional per-frame gross-shift seed (Siril enhancement, no PSS
 * equivalent). When `seed.size() == num_frames`, each frame's MLC search
 * window is pre-centred at that frame's seed value (re-referenced to the
 * best frame: seed[idx] - seed[best]) instead of the previous frame's
 * cumulative shift. This lets a frame whose true displacement exceeds
 * `align_frames_search_width` still be aligned — the residual MLC search
 * only has to cover the seed's error, not the whole jump. The seed must
 * be expressed in the SAME (dy, dx) convention and orientation as the
 * returned `shifts` (the to-align correction); its absolute origin is
 * irrelevant since only seed[idx] - seed[best] is used. An empty seed (or
 * one whose size != num_frames) restores the exact cumulative behaviour.
 * Sourced from Siril's regdata via seed_from_regdata() in mpp.cpp. */

/* Pick a patch on the best (max-quality) frame and compute per-frame integer
 * shifts via a backward-then-forward cumulative-shift loop. `frames` are
 * the Gaussian-blurred mono frames. */
AlignGlobalResult align_global_from_frames(const std::vector<cv::Mat> &frames_mono_blurred,
                                           const std::vector<double> &quality,
                                           const mpp_config_t &cfg,
                                           progress_cb_fn progress = nullptr,
                                           void *progress_user = nullptr,
                                           const std::vector<cv::Vec2d> &seed = {});

/* Streamed overload — identical algorithm, but the provider supplies
 * one frame at a time so the caller controls how the blurred frames are
 * sourced (cached vector, lazy blur from raw cache, or fresh disk read).
 * The thin overload above is implemented in terms of this one.
 *
 * provider_is_thread_safe: when true, the backward and forward sweeps
 * around the best frame run concurrently (two OpenMP sections). Up to
 * ~2x speedup on this stage when `best` lands near the middle of the
 * sequence. Set false (default) for sources whose providers aren't
 * reentrant — notably SEQ_AVI's film_read_frame. Cached-vector providers
 * are always safe. */
AlignGlobalResult align_global_from_provider(const FrameProvider &provider,
                                             int num_frames,
                                             const std::vector<double> &quality,
                                             const mpp_config_t &cfg,
                                             progress_cb_fn progress = nullptr,
                                             void *progress_user = nullptr,
                                             bool provider_is_thread_safe = false,
                                             const std::vector<cv::Vec2d> &seed = {});

/* Build a per-frame gross-shift seed (see the seed param above) from the
 * sequence's existing shift registration data. Uses the first layer whose
 * regparam carries at least one nonzero translation; per frame, extracts the
 * translation via Siril's canonical translation_from_H (dx = H.h02,
 * dy = -H.h12) — the same to-align content translation shift_fit_from_reg
 * applies — and stores it as the port's (dy, dx). Returns an empty vector
 * when no layer carries usable regdata (so the aligner falls back to its
 * cumulative path). Layers whose translations are all zero are treated as
 * quality-only regdata — mpp itself publishes such a layer via
 * mpp_write_quality_to_regdata — and must NOT seed: an all-zero seed pins
 * every search window at zero shift instead of restoring cumulative
 * tracking, which breaks alignment for any frame drifting beyond
 * align_frames_search_width. */
std::vector<cv::Vec2d> seed_from_regdata(struct sequ *seq);

/* Pick the `number_frames` indices that maximise the summed quality
 * within a sliding window of size `region_size`. */
std::vector<int> align_find_best_frames(const std::vector<double> &quality,
                                        int number_frames, int region_size);

struct AlignAverageResult {
	cv::Mat mean_frame;     /* CV_32S, scaled for the 0..65535 convention */
	cv::Vec4i intersection; /* (y_low, y_high, x_low, x_high) in best-frame coords */
	std::vector<int> indices_used;
};

/* Average the raw mono frames (NOT blurred) to build the reference. Uses
 * the global `shifts` to align frames, picks the top-N by quality (with
 * the fast-changing-object sliding window if
 * cfg.align_frames_fast_changing_object), and averages within the inter-
 * frame intersection. */
AlignAverageResult align_average_frame(const std::vector<cv::Mat> &frames_mono_raw,
                                       const std::vector<double> &quality,
                                       const std::vector<cv::Vec2d> &shifts,
                                       const mpp_config_t &cfg,
                                       progress_cb_fn progress = nullptr,
                                       void *progress_user = nullptr);

/* Streamed overload — the provider returns the raw mono frame for
 * index `i`. Frame dims are also passed because no upfront frame is
 * available to read them from. Only the frames listed in the computed
 * `indices` are actually requested (top-N by quality), so on a 10% top-N
 * the streaming path reads ~N/10 frames. */
/* `included` (empty = all): when sized to num_frames, excluded frames are
 * skipped both from the intersection/overlap bounds (a glitch frame with a
 * wild global shift must not collapse the common area) and from the
 * averaged-reference top-N (the caller also passes `quality` pre-masked so
 * excluded frames sort last). */
AlignAverageResult align_average_frame_streamed(const FrameProvider &provider,
                                                int num_frames,
                                                int frame_rows, int frame_cols,
                                                const std::vector<double> &quality,
                                                const std::vector<cv::Vec2d> &shifts,
                                                const mpp_config_t &cfg,
                                                progress_cb_fn progress = nullptr,
                                                void *progress_user = nullptr,
                                                const std::vector<int> &included = {});

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_ALIGN_PRIV_HPP_ */
