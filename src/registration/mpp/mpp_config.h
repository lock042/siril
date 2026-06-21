#ifndef SRC_REGISTRATION_MPP_CONFIG_H_
#define SRC_REGISTRATION_MPP_CONFIG_H_

#include <stdbool.h>
#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* MPP configuration. Fields are added as the corresponding phase needs
 * them. Members not bracketed by a phase comment are Phase-1 (ranking)
 * and active now. */
struct mpp_config {
	/* Frame preprocessing (Phase 1+) */
	int frames_gauss_width;            /* default 7 */
	int align_frames_sampling_stride;  /* default 2 (used by Laplace rank) */
	double rank_laplacian_alpha;       /* default 1/256 (convertScaleAbs scale) */

	/* Brightness normalization (Phase 1) */
	bool frames_normalization;             /* default true */
	int frames_normalization_threshold;    /* default 15, on 0..255 scale */

	/* Source bit depth (the value range of pixels, NOT the CFITSIO bitpix
	 * code). Either 8 (values in 0..255) or 16 (values in 0..65535).
	 * Drives threshold scaling and average-frame upscaling so downstream
	 * code (post-mean_frame) can always treat values as 16-bit-range.
	 * Use mpp_bitdepth_from_fits_bitpix() to convert from Siril's
	 * fit->bitpix. */
	int bitdepth;

	/* Global frame alignment (Phase 2). */
	int align_frames_mode;                      /* enum mpp_align_mode; default MPP_ALIGN_SURFACE */
	int align_frames_search_width;              /* 34 */
	double align_frames_rectangle_scale_factor; /* 3.0 */
	int align_frames_border_width;              /* 10 */
	int align_frames_rectangle_stride;          /* 2 */
	int align_frames_rectangle_black_threshold; /* 10240 */
	double align_frames_rectangle_min_fraction; /* 0.7 */
	int align_frames_average_frame_percent;     /* 5 */
	bool align_frames_fast_changing_object;     /* false (default; PSS defaults true) — build the reference from a short interval */
	int align_frames_best_frames_window_extension; /* 2 */
	bool align_frames_seed_from_regdata;        /* true — seed global align from existing .seq shift regdata (Siril enhancement; no PSS equivalent) */

	/* Alignment-point grid (Phase 3). Per-AP search and stacking parameters
	 * live with their respective modules. */
	int alignment_points_half_box_width;        /* 24  → derives half_patch_width and step_size */
	int alignment_points_search_width;          /* 14 */
	int alignment_points_brightness_threshold;  /* 10  (gets multiplied by 256 for 16-bit) */
	int alignment_points_contrast_threshold;    /* 0   (gets multiplied by 256 for 16-bit) */
	double alignment_points_structure_threshold; /* 0.04 (post-normalisation) */
	double alignment_points_dim_fraction_threshold; /* 0.6 — triggers COM re-centring */
	bool alignment_points_local_search_subpixel; /* false — phase-2 parabolic fit (drizzle path) */

	/* Per-AP frame selection. Two independent controls:
	 *  - REGISTER (alignment_points_frame_*): the upper bound baked at
	 *    registration — how many top-quality frames per AP get a shift
	 *    computed (selected_per_ap). It caps what the stack can use.
	 *  - STACK (stack_frame_*): how many of those the final stack actually
	 *    blends; may be ≤ the register bound, never more (raising it needs
	 *    a re-register). */
	int alignment_points_frame_percent;          /* 100 — register %% of frames per AP */
	int alignment_points_frame_number;           /* -1  — register explicit override (>0 wins; register CLI/GUI leave it -1) */
	int stack_frame_percent;                     /* 100 — stack %% of frames per AP (used when stack_frame_number ≤ 0) */
	int stack_frame_number;                      /* -1  — stack explicit override; positive value overrides stack percent */
	int alignment_points_rank_pixel_stride;      /* 2   — used only for xy/Sobel per-AP rank methods */
	bool alignment_points_de_warp;               /* true */
	double alignment_points_penalty_factor;       /* 0.00025 — weight matrix off-centre penalty */

	double stack_frames_background_fraction;     /* 0.3 */
	double stack_frames_background_blend_threshold; /* 0.2 */
	int stack_frames_background_patch_size;      /* 100 */
	bool stack_skip_failed_aps;                  /* false — drop (frame, AP)
	    contributions whose Stage B shift measurement failed instead of
	    stacking them at the coarse phase-1 estimate.  Per-AP weight sums
	    are reduced accordingly so brightness is unaffected.  An AP whose
	    measurements ALL failed keeps stacking everything (a misaligned
	    patch beats a hole). */

	/* Output bit depth. The weighted merge always runs internally in
	 * 32-bit float; by default the result is packed down to a USHORT
	 * fits. Set true to keep the float result and emit a 32-bit float
	 * fits instead. Stack-time only — driven by the `-32b` CLI flag
	 * (pss / stack_mpp) and the "Force 32b output" stacking-tab
	 * checkbutton. */
	bool output_32bit;                           /* false */

	/* Output scale factor. 1.0 = no upscale (bicubic-no-op path); > 1.0
	 * routes to STScI dobox for mono / RGB input and Bayer dobox for raw
	 * CFA input. Non-integer values (e.g. 1.5) are supported natively —
	 * the dobox pixmap is built with a double scale and output dimensions
	 * round to nearest. */
	double drizzle_scale;

	/* Phase 5b — drizzle backend selection. Now derived from input type
	 * at stack time (see mpp_stack_apply); user-facing surface is just
	 * the scale factor. The enum still exists so the dispatcher can
	 * record which path it picked, and so older sidecar code keeps
	 * compiling. */
	int drizzle_mode;          /* enum mpp_drizzle_mode; default MPP_DRIZZLE_OFF */
	double drizzle_pixfrac;    /* (0, 1]; default 0.7 — drizzle-only */
	int drizzle_kernel;        /* enum mpp_drizzle_kernel; default MPP_KERNEL_TURBO */

	/* AVI Bayer-pattern hint. AVI / film containers carry no Bayer
	 * marker, so when an OSC capture is saved to AVI Siril has no way
	 * to know the mosaic layout — the heuristic in film_read_frame
	 * just sees R==G==B (the raw mosaic value broadcast to three
	 * channels) and classifies the frame as mono, dropping all colour
	 * information. The user must tell us. Default AUTO leaves the
	 * existing film_read_frame heuristic in charge (current behaviour);
	 * MONO is an explicit "yes it is mono, don't probe"; the four
	 * Bayer values route analysis through cv::cvtColor BayerXX2BGR
	 * and Stage C through the dobox Bayer path. Only consulted for
	 * SEQ_AVI sequences; ignored for SER / FITS. */
	int avi_bayer_pattern;     /* enum mpp_avi_bayer; default MPP_AVI_BAYER_AUTO */

	/* Stage C stacking engine (enum mpp_stack_method). PATCH is the PSS
	 * architecture: per-AP patches placed rigidly and mosaiced with blend
	 * windows. WARP builds a smooth dense displacement field and warps each
	 * frame once with cv::remap — required for derotation (a smooth rotation
	 * cannot be represented by rigid per-AP patches) and forced on at stack
	 * time when a .derot sidecar is present. Stack-time only. */
	int stack_method;          /* default MPP_STACK_PATCH */
};

/* Global frame alignment mode (PSS configuration.align_frames_mode).
 * SURFACE uses the auto-picked alignment rectangle + MultiLevelCorrelation
 * (the only Phase-2 path historically); PLANET aligns on the brightness
 * centroid (PSS "Planet" mode / AlignFrames.center_of_gravity), suitable
 * for discs on a dark background where there is no surface detail to
 * correlate. */
enum mpp_align_mode {
	MPP_ALIGN_SURFACE = 0,
	MPP_ALIGN_PLANET  = 1,
};

enum mpp_drizzle_mode {
	MPP_DRIZZLE_OFF     = 0,  /* scale = 1.0; cv::resize no-op path */
	/* slot 1 was MPP_DRIZZLE_BICUBIC (Phase 5a cv::resize upscale) —
	 * deprecated, the OFF path covers the scale=1 case and dobox covers
	 * upscale better. Numbering preserved so saved configs keep meaning. */
	MPP_DRIZZLE_STSCI   = 2,  /* dobox with debayered / mono / RGB input */
	MPP_DRIZZLE_BAYER   = 3,  /* dobox with raw Bayer input → 3-channel output */
};

/* Stage C stacking engine (see mpp_config.stack_method). */
enum mpp_stack_method {
	MPP_STACK_PATCH = 0,   /* per-AP patch blend (PSS architecture) */
	MPP_STACK_WARP  = 1,   /* dense warp-field accumulation (derotation) */
};

/* AVI Bayer-pattern hint values. Ordering matches the GUI combo's item
 * indices (see siril.ui's combo_mpp_avi_bayer); persisted to the sidecar
 * via cfg->avi_bayer_pattern. */
enum mpp_avi_bayer {
	MPP_AVI_BAYER_AUTO = 0,    /* default — film_read_frame heuristic */
	MPP_AVI_BAYER_NONE = 1,    /* explicit mono — skip any debayer */
	MPP_AVI_BAYER_RGGB = 2,
	MPP_AVI_BAYER_BGGR = 3,
	MPP_AVI_BAYER_GBRG = 4,
	MPP_AVI_BAYER_GRBG = 5,
};

enum mpp_drizzle_kernel {
	MPP_KERNEL_SQUARE   = 0,
	MPP_KERNEL_GAUSSIAN = 1,
	MPP_KERNEL_POINT    = 2,
	MPP_KERNEL_TURBO    = 3,
	MPP_KERNEL_LANCZOS2 = 4,
	MPP_KERNEL_LANCZOS3 = 5,
	/* Not a real drizzle kernel — when picked, mpp_stack_apply
	 * routes through the cv::resize bicubic path instead of dobox.
	 * Reliable fallback for rare cases where true drizzle produces
	 * resonance / ringing on fine high-contrast structure (Saturn's
	 * rings, planetary terminator edges). */
	MPP_KERNEL_UPSCALE  = 6,
};

typedef struct mpp_config mpp_config_t;

/* Half-patch width and step size derive from
 * alignment_points_half_box_width; helpers keep the arithmetic in one
 * place. For default half_box_width = 24 both formulas give exact
 * integers; for odd values they fall back to a floor of half_box * 1.5. */
static inline int mpp_cfg_half_patch_width(const struct mpp_config *c) {
	return (c->alignment_points_half_box_width * 3) / 2;  /* round(half_box * 1.5) */
}
static inline int mpp_cfg_step_size(const struct mpp_config *c) {
	/* round((half_patch_width * 4.5) / 3) == round(half_patch_width * 1.5) */
	return (mpp_cfg_half_patch_width(c) * 3) / 2;
}

/* Threshold scale: brightness thresholds are expressed in 0..255 units
 * but compared against pixel values in their native range. Returns 1.0
 * for 8-bit input (0..255) and 256.0 for 16-bit (0..65535). */
static inline double mpp_cfg_threshold_scale(const struct mpp_config *c) {
	return c->bitdepth == 8 ? 1.0 : 256.0;
}

/* Convert Siril's fit->bitpix (CFITSIO BITPIX codes: BYTE_IMG=8,
 * SHORT_IMG=16, USHORT_IMG=20, FLOAT_IMG=-32, ...) to our 8/16 bit-depth
 * convention. Only BYTE_IMG maps to 8; everything else is treated as
 * 16-bit-equivalent, which is how Siril stores both 16-bit integer and
 * float frames internally (WORD or fnormalised float). */
static inline int mpp_bitdepth_from_fits_bitpix(int fits_bitpix) {
	return fits_bitpix == 8 ? 8 : 16;
}

mpp_status_t mpp_config_defaults(mpp_config_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_CONFIG_H_ */
