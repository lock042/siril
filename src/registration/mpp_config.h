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
};

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
