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

	/* Input bit depth — drives threshold scaling. 8 or 16. */
	int bitpix;

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
};

mpp_status_t mpp_config_defaults(mpp_config_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_CONFIG_H_ */
