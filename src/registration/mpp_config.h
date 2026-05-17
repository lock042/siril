#ifndef SRC_REGISTRATION_MPP_CONFIG_H_
#define SRC_REGISTRATION_MPP_CONFIG_H_

#include <stdbool.h>
#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* PSS configuration. Defaults match upstream configuration.py exactly. Fields
 * are added as the corresponding phase needs them; commented with PSS field
 * names where they map 1:1. Members not bracketed by a phase comment are
 * Phase-1 (ranking) and active now. */
struct mpp_config {
	/* Frame preprocessing (Phase 1+) */
	int frames_gauss_width;            /* PSS frames_gauss_width = 7 */
	int align_frames_sampling_stride;  /* PSS align_frames_sampling_stride = 2 (used by Laplace rank) */
	double rank_laplacian_alpha;       /* PSS Frames.alpha = 1/256 (convertScaleAbs scale) */

	/* Brightness normalization (Phase 1) */
	bool frames_normalization;             /* PSS frames_normalization = true */
	int frames_normalization_threshold;    /* PSS frames_normalization_threshold = 15, on 0..255 scale */

	/* Source bit depth (PSS-conceptual: the value range of pixels, NOT the
	 * CFITSIO bitpix code). Either 8 (values in 0..255) or 16 (values in
	 * 0..65535). Drives threshold scaling and average-frame upscaling so
	 * downstream code (post-mean_frame) can always treat values as 16-bit-
	 * range, matching PSS. Use mpp_bitdepth_from_fits_bitpix() to convert
	 * from Siril's fit->bitpix. */
	int bitdepth;

	/* Global frame alignment (Phase 2). All names mirror PSS configuration.py. */
	int align_frames_search_width;              /* 34 */
	double align_frames_rectangle_scale_factor; /* 3.0 */
	int align_frames_border_width;              /* 10 */
	int align_frames_rectangle_stride;          /* 2 */
	int align_frames_rectangle_black_threshold; /* 10240 */
	double align_frames_rectangle_min_fraction; /* 0.7 */
	int align_frames_average_frame_percent;     /* 5 */
	bool align_frames_fast_changing_object;     /* true */
	int align_frames_best_frames_window_extension; /* 2 */

	/* Alignment-point grid (Phase 3). Per-AP search and stacking parameters
	 * live with their respective modules. */
	int alignment_points_half_box_width;        /* 24  → derives half_patch_width and step_size */
	int alignment_points_search_width;          /* 14 */
	int alignment_points_brightness_threshold;  /* 10  (gets multiplied by 256 for 16-bit) */
	int alignment_points_contrast_threshold;    /* 0   (gets multiplied by 256 for 16-bit) */
	double alignment_points_structure_threshold; /* 0.04 (post-normalisation) */
	double alignment_points_dim_fraction_threshold; /* 0.6 — triggers COM re-centring */
	bool alignment_points_local_search_subpixel; /* false — phase-2 parabolic fit (drizzle path) */

	/* Stacking (Phase 5a). */
	int alignment_points_frame_percent;          /* 10  — %% of frames per AP (used when frame_number ≤ 0) */
	int alignment_points_frame_number;           /* -1  — explicit override; positive value overrides percent */
	int alignment_points_rank_pixel_stride;      /* 2   — used only for xy/Sobel per-AP rank methods */
	bool alignment_points_de_warp;               /* true */
	double alignment_points_penalty_factor;       /* 0.00025 — weight matrix off-centre penalty */

	double stack_frames_background_fraction;     /* 0.3 */
	double stack_frames_background_blend_threshold; /* 0.2 */
	int stack_frames_background_patch_size;      /* 100 */

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
};

enum mpp_drizzle_mode {
	MPP_DRIZZLE_OFF     = 0,  /* scale = 1.0; cv::resize no-op path */
	/* slot 1 was MPP_DRIZZLE_BICUBIC (Phase 5a cv::resize upscale) —
	 * deprecated, the OFF path covers the scale=1 case and dobox covers
	 * upscale better. Numbering preserved so saved configs keep meaning. */
	MPP_DRIZZLE_STSCI   = 2,  /* dobox with debayered / mono / RGB input */
	MPP_DRIZZLE_BAYER   = 3,  /* dobox with raw Bayer input → 3-channel output */
};

enum mpp_drizzle_kernel {
	MPP_KERNEL_SQUARE   = 0,
	MPP_KERNEL_GAUSSIAN = 1,
	MPP_KERNEL_POINT    = 2,
	MPP_KERNEL_TURBO    = 3,
	MPP_KERNEL_LANCZOS2 = 4,
	MPP_KERNEL_LANCZOS3 = 5,
};

typedef struct mpp_config mpp_config_t;

/* PSS derives these from alignment_points_half_box_width; helpers keep the
 * arithmetic in one place. For default half_box_width = 24 both formulas
 * give exact integers; for odd values they fall back to Python's
 * banker's-rounding floor behaviour. */
static inline int mpp_cfg_half_patch_width(const struct mpp_config *c) {
	return (c->alignment_points_half_box_width * 3) / 2;  /* PSS: round(half_box * 1.5) */
}
static inline int mpp_cfg_step_size(const struct mpp_config *c) {
	/* PSS: round((half_patch_width * 4.5) / 3) == round(half_patch_width * 1.5) */
	return (mpp_cfg_half_patch_width(c) * 3) / 2;
}

/* Threshold scale: PSS expresses brightness thresholds in 0..255 units but
 * compares against pixel values in their native range. Returns 1.0 for 8-bit
 * input (0..255) and 256.0 for 16-bit (0..65535). */
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
