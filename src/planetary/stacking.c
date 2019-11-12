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

/* This file contains the multi-point stacking for planetary processing.
 */

#include "planetary.h"
#include "caching.h"
#include "stacking.h"
#include "io/sequence.h"
#include "quality.h"
#include "core/proto.h"
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "gui/zones.h"

/* this function is a stacking sum that takes shifts from the multipoint registration instead
 * of taking them from the global registration.  Additionally, it takes only parts of images
 * defined as best for a zone instead of taking the best global images for the sequence.
 * If it works well, it will have to be exploded as a generic processing function. */
int the_old_local_multipoint_sum_stacking(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, nb_images, abort = 0;
	struct weighted_AP *closest_zones_map;	// list of nb_closest_AP AP (fixed) for each pixel

	nb_zones = get_number_of_zones();
	nb_images = round_to_int((double)args->seq->number * args->filtering_percent / 100.0);
	fprintf(stdout, "keeping %d best images for each zone\n", nb_images);

	/* 1. sort images indices from the list of best quality for zones *
	 * In com.stacking_zones[zone].mpregparam we have the normalized quality for each
	 * image for the concerned zone. To speed up the look-up, we create an index here.
	 * */
	BYTE **best_zones = malloc(nb_zones * sizeof(BYTE *));
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		best_zones[zone_idx] = sort_zones_quality(zone_idx, nb_images, args->seq->number);
	}

	/* 2. allocate stacking data *
	 * A sum and a contribution map (number of pixel added for each sum)
	 */
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	unsigned long *sum[3];	// the new image's channels
	int *count[3];	// the new image's contributions count
	sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	count[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!sum[0] || !count[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		sum[1] = sum[0] + nbdata;	// index of green layer in sum[0]
		sum[2] = sum[0] + nbdata*2;	// index of blue layer in sum[0]
		count[1] = count[0] + nbdata;
		count[2] = count[0] + nbdata*2;
	} else {
		sum[1] = NULL;
		sum[2] = NULL;
		count[1] = NULL;
		count[2] = NULL;
	}

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		int zone_idx = 0;
		if (abort) continue;

		while (!best_zones[zone_idx++][frame] && zone_idx < nb_zones);
		if (zone_idx == nb_zones) {
			fprintf(stdout, "frame %d had only bad zones, skipping\n", frame);
			continue;
		}

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "reading best zones for image %d\n", frame);

		/* 3. extract best zones and stack them *
		 * From each image we add the pixels of the best zones to the sum and keep track
		 * of how many times we add each pixel to keep track of the average
		 */
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			if (!best_zones[zone_idx][frame])
				continue;

			stacking_zone *zone = &com.stacking_zones[zone_idx];
			// TODO: iterate over channels
#if 0
			if (args->using_homography && zone->mpregparam[frame].transform) {
				// ECC registration with image transformation
				int side = get_side(zone);
				WORD *buffer = malloc(side * side * sizeof(WORD)); // TODO: prealloc
				if (!buffer) {
					fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
					abort = 1;
					break;
				}

				// read the zone data to buffer
				regdata *regparam = args->seq->regparam[args->layer];
				stacking_zone shifted_zone = { .centre =
					{ .x = zone->centre.x - regparam[frame].shiftx,
						.y = zone->centre.y + regparam[frame].shifty },
					.half_side = zone->half_side };
				copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
				cvTransformBuf(buffer, side, zone->mpregparam[frame].transform);
				add_buf_zone_to_stacking_sum(buffer, 0, zone, frame,
						sum, count, fit.rx, fit.ry);
#ifdef DEBUG_MPP
				save_buffer_tmp(frame, zone_idx, buffer, side);
#endif
				free(buffer);
			} else if (!args->using_homography) {
#endif
				// DFT registration or ECC with translation
				//add_image_zone_to_stacking_sum(&fit, zone, frame, sum, count);
#ifdef DEBUG_MPP
				// to see what's happening with the shifts, use this
				int side = get_side(zone);
				stacking_zone shifted_zone = { .centre =
					{ .x = round_to_int(zone->centre.x - zone->mpregparam[frame].x),
						.y = round_to_int(zone->centre.y + zone->mpregparam[frame].y) },
					.half_side = zone->half_side };

				WORD *buffer = malloc(side * side * sizeof(WORD));
				if (!buffer) {
					fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
					break;
				}

				copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
				save_buffer_tmp(frame, zone_idx, buffer, side);
				free(buffer);
#endif
			//}
			/* TODO someday: instead of copying data like this, create
			 * a mask from the zones, process the mask to unsharpen it
			 * and make it follow the path of the turbulence, and then
			 * use the mask to copy the data. */
		}

		clearfits(&fit);
		// this progress is not as bad as usual, but could be improved by counting
		// the number of zones instead of the number of images
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++)
		free(best_zones[zone_idx]);
	free(best_zones);

	fprintf(stdout, "multipoint stacking ended, creating final image\n");

	/* 4. compute averages for the zones and merge with the reference image for areas
	 * outside zones, store the result in gfit.
	 */
	if (!abort) {
		int i, minzones = nb_images / 3 + 1;	// at least 1
		int from_mpp = 0, from_both = 0, from_ref = 0;
		nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;

		// normalize the sum to the contribution count
		// find the max to normalize to 0..1
		double max = 0.0, invmax, *buf = malloc(nbdata * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i = 0; i < nbdata; i++) {
			buf[i] = (double)sum[0][i] / (double)count[0][i];
			if (buf[i] > max)
				max = buf[i];
		}

		invmax = 1.0 / max;
		for (i = 0; i < nbdata; i++)
			buf[i] *= invmax;

		for (i = 0; i < nbdata; i++) {
			if (count[0][i] > minzones) {
				// already in buffy
				from_mpp++;
			} else if (count[0][i] > 0) {
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

		// make a copy of the reference image to gfit to initialize it
		clearfits(&gfit);
		fits *fit = &gfit;
		if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers, NULL, FALSE))
			return -1;
		gfit.hi = USHRT_MAX;

		double ratio = 1.0;
		ratio = USHRT_MAX_DOUBLE;

		for (i = 0; i < nbdata; i++)
			fit->data[i] = round_to_WORD(buf[i] * ratio);
		free(buf);
	}

	free(sum[0]);
	free(count[0]);
	free(args->global_image);
	return 0;
}

#if 0
/* copy a buffer representing an area into the sum of pixels being stacked.
 * It's done only for one layer. rx and ry are the dimensions of the stacked
 * image and of sum and count matrices.
 */
static void add_buf_zone_to_stacking_sum(WORD *buf, int layer, const stacking_zone *zone,
		int frame, unsigned long *sum[3], int *count[3], unsigned int rx, unsigned int ry) {
	int side = get_side(zone);
	int dst_startx = zone->centre.x - zone->half_side;
	int dst_starty = zone->centre.y - zone->half_side;

	unsigned long *to = sum[layer];
	int *lcount = count[layer];
	int x, y, stride = rx- side;
	int i = 0;
	int o = (ry - dst_starty - side - 1) * rx + dst_startx;

	for (y = 0; y < side; ++y) {
		for (x = 0; x < side; ++x) {
			if (buf[i]) {
				// TODO: is this a good check for empty pixels?
				to[o] += buf[i];
				lcount[o]++;
			}
			i++; o++;
		}
		i += side;
		o += stride;
	}
}
#endif

/* copy the zone from an image to the same zone shifted in the stacked image */
void add_image_zone_to_stacking_sum(fits *fit, const stacking_zone *zone, int frame,
		regdata *regparam, unsigned long *sum[3], int *count[3]) {
	int layer;
	int side = get_side(zone);
	int src_startx = round_to_int(zone->centre.x - zone->half_side - regparam->shiftx - zone->mpregparam[frame].x);
	int src_starty = round_to_int(zone->centre.y - zone->half_side + regparam->shifty - zone->mpregparam[frame].y);
	int dst_startx = zone->centre.x - zone->half_side;
	int dst_starty = zone->centre.y - zone->half_side;

	if (src_startx < 0 || src_startx >= fit->rx - side ||
			src_starty < 0 || src_starty >= fit->ry - side) {
		/* this zone is partly outside the image, we could partially
		 * read it, but for now we just ignore it for stacking */
		return;
	}

	for (layer = 0; layer < fit->naxes[2]; layer++) {
		WORD *from = fit->pdata[layer];
		unsigned long *to = sum[layer];
		int *lcount = count[layer];
		int x, y, stride = fit->rx - side;
		int i = (fit->ry - src_starty - side - 1) * fit->rx + src_startx;
		int o = (fit->ry - dst_starty - side - 1) * fit->rx + dst_startx;

		for (y = 0; y < side; ++y) {
			for (x = 0; x < side; ++x) {
				to[o] += from[i];
				lcount[o]++;
				i++; o++;
			}
			i += stride;
			o += stride;
		}
	}
}

