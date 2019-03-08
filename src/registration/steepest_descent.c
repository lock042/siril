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

/* This is the steepest descent registration for full images.
 * It is defined here as a generic registration method, but is used only 
 * in planetary mode.
 * It runs two passes on the sequence: the first computes image quality, the
 * second computes the shifts. If caching is enabled, during the first pass,
 * gaussian-filtered images are stored for later use.
 * Registration data is always stored on layer 0, because it uses a monochrome
 * version of images so it does not depend on the layer.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "core/siril.h"
#include "core/processing.h"
#include "registration.h"
#include "algos/quality.h"
#include "gui/progress_and_log.h"
#include "stacking/stacking.h"
#include "core/proto.h"	// writeseqfile
#include "planetary/caching.h"
#include "planetary/laplacian_quality.h"
#include "io/sequence.h"

struct regsd_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	struct planetary_cache *cache;
	WORD *reference_image;	// gaussian filtered data for the reference image
	int search_radius;
};

static int regsd_prepare_hook(struct generic_seq_args *args) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	rsdata->cache = calloc(1, sizeof(struct planetary_cache));
	rsdata->cache->use_cached = regargs->use_caching;
	rsdata->cache->cache_data = regargs->use_caching;
	init_caching(args->seq->seqname, rsdata->cache, regargs->kernel_size);

	int ref_idx = args->seq->reference_image;
	assert(ref_idx >= 0);

	if (args->seq->regparam[args->layer]) {
		siril_log_message(_("Recomputing already existing registration\n"));
		rsdata->current_regdata = args->seq->regparam[args->layer];
	} else {
		int i;
		rsdata->current_regdata = malloc(args->seq->number * sizeof(regdata));
		if (rsdata->current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
		for (i = 0; i < args->seq->number; i++) {
			rsdata->current_regdata[i].shiftx = NAN;
			rsdata->current_regdata[i].shifty = NAN;
		}
	}

	/* Get reference frame */
	rsdata->reference_image = get_gaussian_data_for_image(ref_idx, NULL, rsdata->cache);
	if (!rsdata->reference_image) {
		fits fit = { 0 };
		if (seq_read_frame(args->seq, ref_idx, &fit)) {
			siril_log_color_message(_("Could not read reference image for steepest descent registration\n"), "red");
			clearfits(&fit);
			return 1;
		}
		rsdata->reference_image = get_gaussian_data_for_image(ref_idx, &fit, rsdata->cache);
		clearfits(&fit);
		if (!rsdata->reference_image) {
			siril_log_color_message(_("Could not read filtered reference image for steepest descent registration\n"), "red");
			return 1;
		}
	}
	return 0;
}

static int regsd_image_hook(struct generic_seq_args *args,
		int out_index, int in_index, fits *fit, rectangle *_) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	
	WORD *gaussian_data = get_gaussian_data_for_image(in_index, fit, rsdata->cache);
	if (!gaussian_data) {
		siril_log_color_message(_("Could not get filtered image for steepest descent registration\n"), "red");
		return 1;
	}

	// compute shift using the steepest descent algorithm
	int max_radius = rsdata->search_radius;
	rectangle area = { .x = max_radius, .y = max_radius,
		.w = fit->rx - max_radius * 2, .h = fit->ry - max_radius * 2};
	int shiftx = 0, shifty, error;
	// initialize shifts with the previous image's
	if (in_index > 0 && !isnan(rsdata->current_regdata[in_index-1].shiftx)) {
		shiftx = rsdata->current_regdata[in_index-1].shiftx;
		shifty = rsdata->current_regdata[in_index-1].shifty;
	}
	error = search_local_match_gradient(rsdata->reference_image, gaussian_data,
			fit->rx, fit->ry, &area, &area, max_radius, 1, &shiftx, &shifty);
	if (error) {
		siril_log_message(_("Image %d alignment failed\n"), in_index);
		rsdata->current_regdata[in_index].quality = -1.0;
		rsdata->current_regdata[in_index].shiftx = NAN;
		rsdata->current_regdata[in_index].shifty = NAN;
		return 1;
	}

	rsdata->current_regdata[in_index].shiftx = -shiftx;
	rsdata->current_regdata[in_index].shifty = shifty;
	return 0;
}

static int regsd_finalize_hook(struct generic_seq_args *args) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	finalize_caching(rsdata->cache);
	free(rsdata->cache);
	if (!args->retval) {
		args->seq->regparam[args->layer] = rsdata->current_regdata;
		siril_log_message(_("Registration finished.\n"));
		writeseqfile(args->seq);
	} 
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	free(rsdata);
	args->user = NULL;
	return 0;
}

int register_sd(struct registration_args *regargs) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = regargs->seq;
	args->partial_image = FALSE;
	if (regargs->process_all_frames) {
		args->filtering_criterion = stack_filter_all;
		args->nb_filtered_images = regargs->seq->number;
	} else {
		args->filtering_criterion = stack_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	/* we assume layer 0 for this registration */
	args->layer = 0;

	/* check if it's already computed before doing it for nothing */
	if (args->seq->regparam[args->layer]) {
		int i;
		for (i = 0; i < args->seq->number; i++) {
			if (args->filtering_criterion(args->seq, args->layer, i, 0) &&
					args->seq->regparam[args->layer][i].quality < 0.0)
				break;
		}
		if (i == args->seq->number) {
			siril_log_message(_("Registration data already available, not recomputing.\n"));
			return 0;
		}
	}

	/* first pass: compute frame quality, using laplacian of gaussian
	 * filtered frames. This also builds the cache. */
	args->prepare_hook = lapl_prepare_hook;
	args->image_hook = lapl_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = lapl_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Steepest descent registration");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct lapl_data *ldata = malloc(sizeof(struct lapl_data));
	if (!ldata) {
		free(args);
		return -1;
	}
	ldata->use_caching = regargs->use_caching;
	ldata->kernel_size = regargs->kernel_size;
	ldata->for_zones = FALSE;
	args->user = ldata;
	generic_sequence_worker(args);
	if (args->retval)
		return args->retval;

	/* second pass: compute frame shifts relative to best frame */
	args->prepare_hook = regsd_prepare_hook;
	args->image_hook = regsd_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = regsd_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = FALSE;
	args->description = _("Steepest descent registration");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct regsd_data *rsdata = calloc(1, sizeof(struct regsd_data));
	if (!rsdata) {
		free(args);
		return -1;
	}
	rsdata->regargs = regargs;
	rsdata->search_radius = 50;
	args->user = rsdata;

	generic_sequence_worker(args);
	return args->retval;
}

/* makes the sum of differences between pixels of the reference image and the
 * tested image, on the given areas.
 * width and height are the size of both images, sampling occurs one every
 * 'stride' pixels, stride has to be lower than the area width
 * areas must have same width and height
 */
static unsigned long compute_deviation(WORD *ref_frame, WORD *frame,
		int width, int height, rectangle *ref_area, rectangle *area,
		int stride) {
	unsigned long sum = 0ul;
	int x, y;
	int ref_i = ref_area->y * width + ref_area->x;
	int i = area->y * width + area->x;
	for (y = 0; y < area->h; y++) {
		for (x = 0; x < area->w; x += stride) {
			// implementing abs() with unsigned data type
			if (ref_frame[ref_i] > frame[i])
				sum += ref_frame[ref_i] - frame[i];
			else sum += frame[i] - ref_frame[ref_i];
			ref_i += stride;
			i += stride;
		}
		ref_i += width - area->w + 1;
		i += width - area->w + 1;
	}
	return sum;
}

/* The steepest descent image alignment algorithm.
 * Compares areas of a reference frame and the tested frame and looks for the
 * best shift between the two.
 * Sampling occurs one every 'stride' pixels, stride has to be lower than the
 * area width.
 */
int search_local_match_gradient(WORD *ref_frame, WORD *frame, int width, int height,
		rectangle *ref_area, rectangle *area, int search_width,
		int sampling_stride, int *dx_result, int *dy_result)
{
	int iterations = 0;	// for stats

        // Initialize the global optimum with the value at dy=dx=0 or the values
	// of the previous frame
	int dx_min = *dx_result, dy_min = *dy_result;
	unsigned long deviation_min = compute_deviation(ref_frame, frame, width,
			height, area, area, sampling_stride);

	// Start with shift [0, 0]. Stop when a circle with radius 1 around the
	// current optimum reaches beyond the search area.
	while (max(abs(dy_min), abs(dx_min)) < search_width-1) {
		int dx, dy;
		// Go through the neighbours of radius 1 and compute the difference
		// (deviation) between the shifted frame and the corresponding box in
		// the mean frame. Find the minimum "deviation_min_1".
		int dx_min_1 = INT_MAX, dy_min_1 = INT_MAX;
		unsigned long deviation_min_1 = ULONG_MAX;
		for (dx = dx_min-1; dx <= dx_min+1; dx++) {
			for (dy = dy_min-1; dy <= dy_min+1; dy++) {
				unsigned long deviation;
				if (dy == 0 && dx == 0) continue;
				rectangle test_area = {
					.x = area->x - dx, .y = area->y - dy,
					.w = area->w, .h = area->h };

				deviation = compute_deviation(ref_frame, frame,
						width, height, area, &test_area,
						sampling_stride);

				if (deviation < deviation_min_1) {
					deviation_min_1 = deviation;
					dx_min_1 = dx;
					dy_min_1 = dy;
				}
			}
		}

		// If for the current center the match is better than for all
		// neighboring points, a local optimum is found.
		if (deviation_min_1 >= deviation_min) {
			*dx_result = -dx_min;
			*dy_result = dy_min;
			fprintf(stdout, "found shift after %d iterations (%d, %d)\n", iterations, dx_min, dy_min);
			return 0;
		}

		// Otherwise, update the current optimum and continue.
		deviation_min = deviation_min_1;
		dx_min = dx_min_1;
		dy_min = dy_min_1;
		iterations++;
	}
	// If within the maximum search radius no optimum could be found, return [0, 0].
	*dx_result = 0;
	*dy_result = 0;
	fprintf(stdout, "shift not found after %d iterations\n", iterations);
	return -1;
}
