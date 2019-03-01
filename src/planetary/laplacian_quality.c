/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
 */

/* This file provides the frame or area quality evaluation based on the
 * Laplacian of the Gaussian-filtered frames or areas. */

#include "core/siril.h"
#include "core/processing.h"
#include "caching.h"
#include "laplacian_quality.h"
#include "gui/progress_and_log.h"
#include "gui/zones.h"
#include "registration/registration.h"
#include "core/proto.h"	// writeseqfile
#include "stacking/stacking.h"	// stack_filter_all

static double variance(WORD *set, int size);
static double variance_area(WORD *image, int width, int height, rectangle *area);

int lapl_prepare_hook(struct generic_seq_args *args) {
	struct lapl_data *ldata = args->user;
	ldata->cache = calloc(1, sizeof(struct planetary_cache));
	ldata->cache->use_cached = ldata->use_caching;
	ldata->cache->cache_data = ldata->use_caching;
	init_caching(args->seq->seqname, ldata->cache, ldata->kernel_size);

	if (args->seq->regparam[0]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		ldata->current_regdata = args->seq->regparam[0];
		if (!ldata->for_zones) {
			/* we reset all global values as we may register different images */
			memset(ldata->current_regdata, 0, args->seq->number * sizeof(regdata));
		}
	} else {
		ldata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (ldata->current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
	}

	if (ldata->for_zones) {
		// allocate mpregparam if needed
		int zone_idx;
		ldata->nb_zones = get_number_of_zones();
		if (ldata->nb_zones < 1)
			return 1;

		for (zone_idx = 0; zone_idx < ldata->nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			if (!zone->mpregparam) {
				zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
				if (!zone->mpregparam) {
					fprintf(stderr, "memory allocation failed\n");
					return -2;
				}
			}
		}

		ldata->max = calloc(ldata->nb_zones, sizeof(double));
	}
	return 0;
}

int lapl_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct lapl_data *ldata = args->user;

	WORD *gaussian_data = get_gaussian_data_for_image(in_index, fit, ldata->cache);
	WORD *laplacian_data = get_laplacian_data_for_image(in_index,
			gaussian_data, fit->rx, fit->ry, ldata->cache);

	// optional: apply downsampling?

	if (!ldata->for_zones) {
		ldata->current_regdata[in_index].quality = variance(laplacian_data, fit->rx * fit->ry);
		//fprintf(stdout, "variance for %d: %g\n", in_index, ldata->current_regdata[in_index].quality);
	} else {
		int zone_idx;
		for (zone_idx = 0; zone_idx < ldata->nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			int side = get_side(zone);
			// TODO: make sure zone is in image
			rectangle shifted_zone = {
				.x = zone->centre.x - ldata->current_regdata[in_index].shiftx,
				.y = zone->centre.y + ldata->current_regdata[in_index].shifty,
				.w = side, .h = side };

			/* compute the quality for each zone */
			zone->mpregparam[in_index].quality =
				variance_area(laplacian_data, fit->rx, fit->ry, &shifted_zone);
			if (ldata->max[zone_idx] < zone->mpregparam[in_index].quality)
				ldata->max[zone_idx] = zone->mpregparam[in_index].quality;
		}
	}

	free(laplacian_data);
	free(gaussian_data);
	return 0;
}

int lapl_finalize_hook(struct generic_seq_args *args) {
	struct lapl_data *ldata = args->user;
	finalize_caching(ldata->cache);
	free(ldata->cache);
	if (!args->retval) {
		if (!ldata->for_zones) {
			args->seq->regparam[0] = ldata->current_regdata;
			/* find best frame and normalize quality */

			int q_index = normalize_quality_data(ldata->current_regdata, args->seq->number);
			siril_log_message(_("Laplacian quality finished.\n"));
			siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);
			args->seq->reference_image = q_index;
			writeseqfile(args->seq);
		} else {
			// normalize zone quality
			int frame, zone_idx;
			for (frame = 0; frame < args->seq->number; frame++) {
				for (zone_idx = 0; zone_idx < ldata->nb_zones; zone_idx++) {
					stacking_zone *zone = &com.stacking_zones[zone_idx];
					zone->mpregparam[frame].quality /= ldata->max[zone_idx];
				}
			}
		}
	}
	else {
		siril_log_message(_("Laplacian quality aborted.\n"));
	}
	if (ldata->for_zones) free(ldata->max);
	free(ldata);
	args->user = NULL;
	return 0;
}

static double variance(WORD *set, int size) {
	int i;
	unsigned long sum = 0ul, squared_sum = 0ul;
	// could be parallelized with sum reductions, is it worth it?
	for (i = 0; i < size; i++) {
		sum += set[i];
		squared_sum += set[i] * set[i];
	}

	double mean = (double)sum / (double)size;
	return (double)squared_sum / (mean * mean);
}

static double variance_area(WORD *image, int width, int height, rectangle *area) {
	int x, y, size = width * height;
	unsigned long sum = 0ul, squared_sum = 0ul;
	int i = area->y * width + area->x;
	for (y = 0; y < area->h; y++) {
		for (x = 0; x < area->w; x++) {
			sum += image[i];
			squared_sum += image[i] * image[i];
			i++;
		}
		i += width - area->w + 1;
	}

	double mean = (double)sum / (double)size;
	return (double)squared_sum / (mean * mean);
}

int laplace_quality_for_zones(sequence *seq, gboolean process_all_frames, gboolean use_caching, int kernel_size) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = seq;
	args->partial_image = FALSE;
	if (process_all_frames) {
		args->filtering_criterion = stack_filter_all;
		args->nb_filtered_images = seq->number;
	} else {
		args->filtering_criterion = stack_filter_included;
		args->nb_filtered_images = seq->selnum;
	}
	args->layer = 0;
	args->prepare_hook = lapl_prepare_hook;
	args->image_hook = lapl_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = lapl_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Image quality computation");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;	// ?
	args->parallel = TRUE;

	struct lapl_data *ldata = malloc(sizeof(struct lapl_data));
	if (!ldata) {
		free(args);
		return -1;
	}
	ldata->use_caching = use_caching;
	ldata->kernel_size = kernel_size;
	ldata->for_zones = TRUE;
	args->user = ldata;

	generic_sequence_worker(args);

	return args->retval;
}
