#include "registration/mpp_config.h"

mpp_status_t mpp_config_defaults(mpp_config_t *cfg) {
	if (!cfg)
		return MPP_EINVAL;
	cfg->frames_gauss_width = 7;
	cfg->align_frames_sampling_stride = 2;
	cfg->rank_laplacian_alpha = 1.0 / 256.0;
	cfg->frames_normalization = true;
	cfg->frames_normalization_threshold = 15;
	cfg->bitdepth = 16;

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

	/* Phase 3 — alignment-point grid */
	cfg->alignment_points_half_box_width = 24;
	cfg->alignment_points_search_width = 14;
	cfg->alignment_points_brightness_threshold = 10;
	cfg->alignment_points_contrast_threshold = 0;
	cfg->alignment_points_structure_threshold = 0.04;
	cfg->alignment_points_dim_fraction_threshold = 0.6;
	cfg->alignment_points_local_search_subpixel = false;
	return MPP_OK;
}
