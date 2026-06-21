/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 *
 * This implementation of multipoint registration & stacking is based on
 * PlanetarySystemStacker by Rolf Hempel:
 *     https://github.com/Rolf-Hempel/PlanetarySystemStacker
 */

#include "registration/mpp/mpp_config.h"

mpp_status_t mpp_config_defaults(mpp_config_t *cfg) {
	if (!cfg)
		return MPP_EINVAL;
	cfg->frames_gauss_width = 7;
	cfg->align_frames_sampling_stride = 2;
	cfg->rank_laplacian_alpha = 1.0 / 256.0;
	cfg->frames_normalization = true;
	cfg->frames_normalization_threshold = 15;
	cfg->bitdepth = 16;

	/* Phase 2 — global alignment. Default Planet (centroid, full-disc
	 * planetary) — a deliberate divergence from PSS, which defaults to
	 * Surface. Surface remains available via -align=surface / the combo. */
	cfg->align_frames_mode = MPP_ALIGN_PLANET;
	cfg->align_frames_search_width = 34;
	cfg->align_frames_rectangle_scale_factor = 3.0;
	cfg->align_frames_border_width = 10;
	cfg->align_frames_rectangle_stride = 2;
	cfg->align_frames_rectangle_black_threshold = 10240;
	cfg->align_frames_rectangle_min_fraction = 0.7;
	cfg->align_frames_average_frame_percent = 5;
	cfg->align_frames_fast_changing_object = false;   /* Siril default off; PSS defaults on */
	cfg->align_frames_best_frames_window_extension = 2;
	cfg->align_frames_seed_from_regdata = true;

	/* Phase 3 — alignment-point grid */
	cfg->alignment_points_half_box_width = 24;
	cfg->alignment_points_search_width = 14;
	cfg->alignment_points_brightness_threshold = 10;
	cfg->alignment_points_contrast_threshold = 0;
	cfg->alignment_points_structure_threshold = 0.04;
	cfg->alignment_points_dim_fraction_threshold = 0.6;
	cfg->alignment_points_local_search_subpixel = false;

	/* Phase 5a — per-AP frame selection. Both controls default to 100%:
	 * register every frame (the cap) and stack every registered frame.
	 * Lower the register percent only to save compute when you know you
	 * won't stack more than that fraction; lower the stack percent for
	 * auto top-N culling at stack time. */
	cfg->alignment_points_frame_percent = 100;   /* register cap */
	cfg->alignment_points_frame_number = -1;
	cfg->stack_frame_percent = 100;              /* stack selection */
	cfg->stack_frame_number = -1;
	cfg->alignment_points_rank_pixel_stride = 2;
	cfg->alignment_points_de_warp = true;
	cfg->alignment_points_penalty_factor = 0.00025;
	cfg->stack_frames_background_fraction = 0.3;
	cfg->stack_frames_background_blend_threshold = 0.2;
	cfg->stack_frames_background_patch_size = 100;
	cfg->stack_skip_failed_aps = false;
	cfg->output_32bit = false;
	cfg->drizzle_scale = 1.0;

	/* Drizzle backend defaults: off. drizzle_mode is auto-set at stack
	 * time from input type (mpp_classify_sequence_input); pixfrac /
	 * kernel only consulted when scale > 1.0.
	 *
	 * Turbo kernel selected as default after the Phase 7.4 / 7.5 kernel
	 * sweep: same wall-clock as bicubic-2x on test-big.ser (174.6 s vs
	 * 177.1 s) and 2.8x faster than square with no measurable quality
	 * difference (SSIM avg 0.9686 vs 0.9701, span across all four kernels
	 * < 0.15 %). */
	cfg->drizzle_mode    = MPP_DRIZZLE_OFF;
	cfg->drizzle_pixfrac = 0.7;
	cfg->drizzle_kernel  = MPP_KERNEL_TURBO;
	cfg->avi_bayer_pattern = MPP_AVI_BAYER_AUTO;
	cfg->stack_method = MPP_STACK_PATCH;
	return MPP_OK;
}
