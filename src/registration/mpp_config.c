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

	/* Phase 2 — global alignment */
	cfg->align_frames_search_width = 34;
	cfg->align_frames_rectangle_scale_factor = 3.0;
	cfg->align_frames_border_width = 10;
	cfg->align_frames_rectangle_stride = 2;
	cfg->align_frames_rectangle_black_threshold = 10240;
	cfg->align_frames_rectangle_min_fraction = 0.7;
	cfg->align_frames_average_frame_percent = 5;
	cfg->align_frames_fast_changing_object = true;
	cfg->align_frames_best_frames_window_extension = 2;
	return MPP_OK;
}
