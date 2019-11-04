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

/* This file provide the steepest descent registration algorithm for alignment
 * points then stacks the areas to produce the final result in one pass.
 *
 * The monochrome gaussian-filtered reference image is used for registration,
 * then stacking is done with the original images.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "planetary.h"
#include "caching.h"
#include "quality.h"
#include "stacking.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "gui/progress_and_log.h"
#include "gui/zones.h"
#include "registration/registration.h"
#include "algos/statistics.h"

#define DEBUG_MPP
#define USE_REF_FROM_SEQUENCE_INSTEAD_OF_THE_STACKED_ONE 0

static void zone_to_rectangle(stacking_zone *zone, rectangle *rectangle);

struct mppsd_data {
	BYTE **best_zones;	// index of best zones per image

	/* registration data */
	//struct mpr_args *mpr;
	struct planetary_cache *cache;
	int kernel_size;
	int search_radius;
	int nb_zones;		// caching the number of declared zones
	fits *refimage;

	WORD *reference_data;	// gaussian filtered data for the reference image
	float *normalized_refdata;	// the same normalized
	float *ref_gradient;	// gradients of the normalized data

	/* stacking data */
	unsigned long *sum[3];	// the new image's channels
	int *count[3];	// the new image's contributions count
};

static int mppsd_prepare_hook(struct generic_seq_args *args) {
	struct mppsd_data *mppdata = args->user;
	mppdata->cache = calloc(1, sizeof(struct planetary_cache));
	mppdata->cache->use_cached = TRUE;
	mppdata->cache->cache_data = FALSE;
	// ^ should have been done before, and there's only one pass here
	init_caching(args->seq->seqname, mppdata->cache, mppdata->kernel_size);
	unsigned int nbdata = args->seq->ry * args->seq->rx;

	/* Get reference frame */
	WORD ref_min, ref_max;
#if USE_REF_FROM_SEQUENCE_INSTEAD_OF_THE_STACKED_ONE
	int ref_idx = args->seq->reference_image;
	mppdata->reference_data = get_gaussian_data_for_image_in_seq(args->seq, ref_idx, mppdata->cache);
	if (!mppdata->reference_data) {
		siril_log_color_message(_("Could not read reference image for multipoint processing\n"), "red");
		return 1;
	}
	siril_stats_ushort_minmax(&ref_min, &ref_max, mppdata->reference_data, nbdata*args->seq->nb_layers);
#else
	// the stacked reference image is much less noisy, so for now we don't
	// take the gaussian-filtered version
	// TODO: use monochrome conversion if it's not a monochrome image
	mppdata->reference_data = mppdata->refimage->data;
	image_find_minmax(mppdata->refimage);
	ref_min = mppdata->refimage->mini;
	ref_max = mppdata->refimage->maxi;
#endif
	// Normalization
	mppdata->normalized_refdata = malloc(nbdata * args->seq->nb_layers * sizeof(float));
	mppdata->ref_gradient = malloc(nbdata * args->seq->nb_layers * sizeof(float));
	normalize_data(mppdata->reference_data, nbdata*args->seq->nb_layers,
			ref_min, ref_max, mppdata->normalized_refdata);
	compute_gradients_of_buffer(mppdata->normalized_refdata,
			args->seq->rx, args->seq->ry, mppdata->ref_gradient);

	mppdata->sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	mppdata->count[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!mppdata->sum[0] || !mppdata->count[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		mppdata->sum[1] = mppdata->sum[0] + nbdata;	// index of green layer in sum[0]
		mppdata->sum[2] = mppdata->sum[0] + nbdata*2;	// index of blue layer in sum[0]
		mppdata->count[1] = mppdata->count[0] + nbdata;
		mppdata->count[2] = mppdata->count[0] + nbdata*2;
	} else {
		mppdata->sum[1] = NULL;
		mppdata->sum[2] = NULL;
		mppdata->count[1] = NULL;
		mppdata->count[2] = NULL;
	}

	// allocate the mpregparams in a thread-safe environment
	int zone_idx;
	for (zone_idx = 0; zone_idx < mppdata->nb_zones; zone_idx++) {
		stacking_zone *zone = &com.stacking_zones[zone_idx];
		if (!zone->mpregparam) {
			zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
			if (!zone->mpregparam) {
				fprintf(stderr, "memory allocation failed\n");
				return 1;
			}
		}
	}
	return 0;
}


static int mppsd_image_hook(struct generic_seq_args *args,
		int out_index, int in_index, fits *fit, rectangle *_) {
	struct mppsd_data *mppdata = args->user;
	int zone_idx;
	WORD *gaussian_data = get_gaussian_data_for_image(in_index, fit, mppdata->cache);
	if (!gaussian_data) {
		siril_log_color_message(_("Could not get filtered image for multipoint processing\n"), "red");
		return 1;
	}
	WORD min, max;
	int nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
	siril_stats_ushort_minmax(&min, &max, gaussian_data, nbdata);
	float *normalized_gaussian_data = malloc(nbdata * sizeof(float));
	normalize_data(gaussian_data, nbdata*args->seq->nb_layers,
			min, max, normalized_gaussian_data);

	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;

	/* reading zones in the sequence images */
	for (zone_idx = 0; zone_idx < mppdata->nb_zones; zone_idx++) {
		if (!mppdata->best_zones[zone_idx][in_index])
			continue;

		stacking_zone *zone = &com.stacking_zones[zone_idx];
		stacking_zone shifted_zone = { .centre =
			{ .x = round_to_int(zone->centre.x - regparam[in_index].shiftx),
				.y = round_to_int(zone->centre.y + regparam[in_index].shifty) },
			.half_side = zone->half_side };

		/* AP registration */
		int max_radius = mppdata->search_radius;
		int shiftx = 0, shifty = 0, error;
		if (in_index > 0 && !isnan(zone->mpregparam[in_index-1].x)) {
			shiftx = zone->mpregparam[in_index-1].x;
			shiftx = zone->mpregparam[in_index-1].y;
		}
		rectangle ref_area, im_area;
		zone_to_rectangle(zone, &ref_area);
		zone_to_rectangle(&shifted_zone, &im_area);
		if (is_area_in_image(&im_area, args->seq)) {
			/* error = search_local_match_gradient(mppdata->reference_data,
					gaussian_data, fit->rx, fit->ry, &ref_area,
					&im_area, max_radius, 1, &shiftx, &shifty); */
			error = search_local_match_gradient_float(mppdata->normalized_refdata,
					mppdata->ref_gradient, normalized_gaussian_data,
					AP_DEVIATION, fit->rx, fit->ry, &ref_area,
					&im_area, max_radius, 1, &shiftx, &shifty);
		}
		else error = 1;

		if (error) {
			if (error == -1)
				fprintf(stderr, "could not match zone %d of image %d with steepest gradient\n", zone_idx, in_index);
			else if (error == 1)
				fprintf(stderr, "zone %d was out of image %d, not aligning (steepest gradient)\n", zone_idx, in_index);
			zone->mpregparam[in_index].x = NAN;
			zone->mpregparam[in_index].y = NAN;
			zone->mpregparam[in_index].quality = -1.0;
			// TODO: handle failed APs
			continue;
		}

		zone->mpregparam[in_index].x =  -shiftx + regparam[in_index].shiftx;
		zone->mpregparam[in_index].y =   shifty + regparam[in_index].shifty;
		fprintf(stdout, "frame %d, zone %d local shifts: %d,%d\n",
				in_index, zone_idx, -shiftx, shifty);

		/* AP stacking */
		add_image_zone_to_stacking_sum(fit, zone, in_index,
				mppdata->sum, mppdata->count);
#ifdef DEBUG_MPP
		// to see what's happening with the shifts, use this
		{
			// dump the gaussian data or the original data
			int dump_gaussian = 1;
			// dump shifted zone or the compared zone
			int dump_shifted = 1;

			int side = get_side(zone);
			WORD *buffer = malloc(side * side * sizeof(WORD));
			if (dump_shifted) {
				// same coordinates as in add_image_zone_to_stacking_sum()
				shifted_zone.centre.x = round_to_int(zone->centre.x - zone->mpregparam[in_index].x);
				shifted_zone.centre.y = round_to_int(zone->centre.y + zone->mpregparam[in_index].y);
			}
			if (dump_gaussian) {
				copy_image_buffer_zone_to_buffer(gaussian_data, fit->rx, fit->ry, &shifted_zone, buffer);
			} else {
				copy_image_zone_to_buffer(fit, &shifted_zone, buffer, args->layer);
			}
			save_buffer_tmp(in_index, zone_idx, buffer, side);
			free(buffer);
		}
#endif

	}
	free(gaussian_data);
	free(normalized_gaussian_data);

	return 0;
}

/* This is still the old normalization, to be rewritten */
static int mppsd_finalize_hook(struct generic_seq_args *args) {
	struct mppsd_data *mppdata = args->user;
	int zone_idx;
	for (zone_idx = 0; zone_idx < mppdata->nb_zones; zone_idx++)
		free(mppdata->best_zones[zone_idx]);
	free(mppdata->best_zones);

	fprintf(stdout, "multipoint stacking ended, creating final image\n");

	/* compute averages for the zones and merge with the reference image for
	 * areas outside zones, store the result in gfit.
	 */
	int nb_images = args->seq->number;
	int i;
        //minzones = nb_images / 3 + 1;	// at least 1
	//int from_mpp = 0, from_both = 0, from_ref = 0;
	int nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;

	// normalize the sum to the contribution count
	// find the max to normalize to 0..1
	double max = 0.0, invmax, *buf = malloc(nbdata * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
	for (i = 0; i < nbdata; i++) {
		buf[i] = (double)mppdata->sum[0][i] / (double)mppdata->count[0][i];
		if (buf[i] > max)
			max = buf[i];
	}

	invmax = 1.0 / max;
	for (i = 0; i < nbdata; i++)
		buf[i] *= invmax;

#if 0
	for (i = 0; i < nbdata; i++) {
		if (mppdata->count[0][i] > minzones) {
			// already in buffy
			from_mpp++;
		}
		/* TODO: put the reference image in global_image */
		else if (mppdata->count[0][i] > 0) {
			buf[i] = (args->global_image[i] + buf[i]) / 2.0;
			from_both++;
		} else {
			buf[i] = args->global_image[i];
			from_ref++;
		}
	}
	fprintf(stdout, " pixels from best local zones: %d\n", from_mpp);
	fprintf(stdout, " pixels from best local zones mixed with reference: %d\n", from_both);
	fprintf(stdout, " pixels from reference image (global): %d\n", from_ref);
#endif

	// make a copy of the reference image to gfit to initialize it
	clearfits(&gfit);
	fits *fit = &gfit;
	if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers, NULL))
		return -1;
	gfit.hi = USHRT_MAX;

	double scale = USHRT_MAX_DOUBLE;

	for (i = 0; i < nbdata; i++)
		fit->data[i] = round_to_WORD(buf[i] * scale);
	free(buf);

	free(mppdata->sum[0]);
	free(mppdata->count[0]);
	//free(args->global_image);

	return 0;
}


void the_multipoint_steepest_descent_registration_and_stacking(struct mpr_args *mprargs) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	if (!args) return;
	args->seq = mprargs->seq;
	args->partial_image = FALSE;
	// the criterion should be: has_image_some_good_zones
	args->filtering_criterion = mprargs->filtering_criterion;
	args->filtering_parameter = mprargs->filtering_parameter;
	args->nb_filtered_images = -1;
	args->layer = mprargs->layer;
	args->prepare_hook = mppsd_prepare_hook;
	args->image_hook = mppsd_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = mppsd_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Multi-point processing");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct mppsd_data *mppdata = calloc(1, sizeof(struct mppsd_data));
	if (!mppdata) {
		PRINT_ALLOC_ERR;
		free(args);
		return;
	}
	int nb_images = mprargs->seq->number;
	int nb_best = round_to_int(nb_images * mprargs->filtering_percent);
	create_index_of_best_zones(&mppdata->best_zones, nb_best, nb_images);
	mppdata->kernel_size = mprargs->kernel_size;
	mppdata->nb_zones = get_number_of_zones();
	mppdata->search_radius = 50;	// TODO: pass it
	mppdata->refimage = mprargs->refimage;
	args->user = mppdata;

	generic_sequence_worker(args);
}

static void zone_to_rectangle(stacking_zone *zone, rectangle *rectangle) {
	int side = get_side(zone);
	rectangle->x = zone->centre.x - zone->half_side;
	rectangle->y = zone->centre.y - zone->half_side;
	rectangle->w = side;
	rectangle->h = side;
}

void compute_gradients_of_buffer(float *ref_zone, int w, int h, float *gradient) {
	int x, y, i;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			i = y * w + x;
			if (x == 0)
				gradient[i] = fabs(ref_zone[i] - ref_zone[i+1]);
			else if (x == w-1)
				gradient[i] = fabs(ref_zone[i] - ref_zone[i-1]);
			else {
				gradient[i] = fabs(ref_zone[i] - ref_zone[i-1]);
				gradient[i] += fabs(ref_zone[i] - ref_zone[i+1]);
			}

			if (y == 0)
				gradient[i] += fabs(ref_zone[i] - ref_zone[i+w]);
			else if (y == h-1)
				gradient[i] += fabs(ref_zone[i] - ref_zone[i-w]);
			else {
				gradient[i] += fabs(ref_zone[i] - ref_zone[i-w]);
				gradient[i] += fabs(ref_zone[i] - ref_zone[i+w]);
			}
		}
	}
}

