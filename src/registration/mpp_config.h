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
};

mpp_status_t mpp_config_defaults(mpp_config_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_CONFIG_H_ */
