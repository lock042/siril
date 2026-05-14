#include "registration/mpp_config.h"

mpp_status_t mpp_config_defaults(mpp_config_t *cfg) {
	if (!cfg)
		return MPP_EINVAL;
	cfg->frames_gauss_width = 7;
	cfg->align_frames_sampling_stride = 2;
	cfg->rank_laplacian_alpha = 1.0 / 256.0;
	cfg->frames_normalization = true;
	cfg->frames_normalization_threshold = 15;
	cfg->bitpix = 16;
	return MPP_OK;
}
