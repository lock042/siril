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

/* This is the main file for quality evaluation of frames (global) or alignment
 * points (local) for planetary stacking.
 * The laplacian-based quality evaluation is defined in laplacian_quality.c for
 * more clarity.
 * The so-called siril quality uses the algorithm from src/algos/quantize.c to
 * compute quality and this sequence processing is implemented below. If it is
 * kept as an alternative, it should be converted to a generic processing
 * function like laplacian_quality.
 */

#include "core/proto.h"
#include "laplacian_quality.h"
#include "quality.h"
#include "gui.h"
#include "stacking.h"
#include "algos/quality.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"
#include "gui/zones.h"
#include <stdlib.h>
#include <stdio.h>

/*********************** ANALYSIS *************************/
gpointer the_multipoint_analysis(gpointer ptr) {
	struct mpr_args *args = (struct mpr_args*)ptr;
	int retval = the_laplace_multipoint_quality_analysis(args);
	//int retval = the_siril_multipoint_quality_analysis(args);
	siril_add_idle(end_mpp_analysis, NULL);
	free(args);
	return GINT_TO_POINTER(retval);
}

/* evaluates quality of each zone for each image and saves it in mpregparam */
int the_laplace_multipoint_quality_analysis(struct mpr_args *args) {
	return laplace_quality_for_zones(args->seq, FALSE, TRUE, 3);
}

int the_siril_multipoint_quality_analysis(struct mpr_args *args) {
	int retval = 0, frame, zone_idx, abort = 0;

	int nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}
	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;
	WORD **buffer = calloc(nb_zones, sizeof(WORD *));
	double *max = calloc(nb_zones, sizeof(double));

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "analysing zones for image %d\n", frame);

		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			stacking_zone shifted_zone = { .centre =
				{ .x = zone->centre.x - regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };
			int side = get_side(zone);

			if (!zone->mpregparam) {
				zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
				if (!zone->mpregparam) {
					fprintf(stderr, "memory allocation failed\n");
					abort = 1;
					break;
				}
			}

			// copy the zone to a buffer and evaluate quality,
			// cannot be done in-place
			if (!buffer[zone_idx]) {
				/* this is not thread-safe, a buffer is required per thread, it
				 * can be done with a single allocation with round-robin
				 * distribution of zones across threads */
				buffer[zone_idx] = malloc(side * side * sizeof(WORD));
			}

			if (!copy_image_zone_to_buffer(&fit, &shifted_zone,
						buffer[zone_idx], args->layer)) {
				zone->mpregparam[frame].quality =
					QualityEstimateBuf(buffer[zone_idx], side, side);
				if (max[zone_idx] < zone->mpregparam[frame].quality)
					max[zone_idx] = zone->mpregparam[frame].quality;
			}
		}

		clearfits(&fit);
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	// normalization
	for (frame = 0; frame < args->seq->number; frame++) {
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			zone->mpregparam[frame].quality /= max[zone_idx];
		}
	}	

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		if (buffer[zone_idx])
			free(buffer[zone_idx]);
	}
	free(buffer);
	free(max);

	return retval;
}


/* sort images indices from the list of best quality for zones *
 * In com.stacking_zones[zone].mpregparam we have the normalized quality for each
 * image for the concerned zone. To speed up the look-up, we create an index here.
 */
void create_index_of_best_zones(BYTE ***best_zones, int nb_best, int nb_images) {
	int zone_idx, nb_zones = get_number_of_zones();
	*best_zones = malloc(nb_zones * sizeof(BYTE *));
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		(*best_zones)[zone_idx] = sort_zones_quality(zone_idx, nb_best, nb_images);
	}
}

/* Sort images indices from the list of best quality for zones *
 * In com.stacking_zones[zone].mpregparam we have the normalized quality for each
 * image for the concerned zone. To speed up the look-up, we create an index here.
 * */
BYTE * sort_zones_quality(int zone_idx, int nb_best, int nb_seq_images) {
	int i;
	BYTE *index = malloc(nb_seq_images);
	// finding the nth element of an unsorted set: sort it, it's easier
	int *indices = apregdata_best(com.stacking_zones[zone_idx].mpregparam, nb_seq_images);

	// output: index with ones if image is among the best
	for (i = 0; i < nb_seq_images; i++)
		index[indices[i]] = i < nb_best;

	free(indices);
	return index;
}

