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

/* This file provide a planetary-tailored ECC registration.
 * It should not be use anymore because the steepest descent algorithm is more
 * stable and faster.
 */

#include "planetary.h"
#include "caching.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "gui/zones.h"

int the_multipoint_ecc_registration(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, abort = 0;
	WORD **references;	// an array of pointers because their size varies
	stacking_zone *zone;

	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;

	nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}

	/* read zones of the reference image */
	references = malloc(nb_zones*sizeof(WORD *));
	if (!references) {
		fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
		return -1;
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	fits reffits = { 0 };
	if (seq_read_frame(args->seq, args->seq->reference_image, &reffits))
		return -1;
	fprintf(stdout, "Using the reference image from sequence (%d) instead "
			"of the stacked reference\n", args->seq->reference_image);
#endif
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		zone = &com.stacking_zones[zone_idx];
		int side = get_side(zone);
		references[zone_idx] = malloc(side * side * sizeof(WORD));
		if (!references[zone_idx]) {
			fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
			return -1;
		}
#ifdef USE_SEQUENCE_REF_IMAGE
		copy_image_zone_to_buffer(&reffits, zone, references[zone_idx], args->layer);
#else
		copy_image_zone_to_buffer(args->refimage, zone, references[zone_idx], args->layer);
#endif
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	clearfits(&reffits);
#endif

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		/* TODO: this filtering is required for global mode, not for local */
		/*if (!args->filtering_criterion(args->seq, args->layer,
					frame, args->filtering_parameter) || abort)
			continue;*/

		if (seq_read_frame(args->seq, frame, &fit) || image_find_minmax(&fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "\naligning zones for image %d\n", frame);

		/* reading zones in the reference image */
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			reg_ecc reg_param = { 0 };
			if (abort) continue;
			zone = &com.stacking_zones[zone_idx];
			stacking_zone shifted_zone = { .centre =
				{ .x = round_to_int(zone->centre.x - regparam[frame].shiftx),
					.y = round_to_int(zone->centre.y + regparam[frame].shifty) },
				.half_side = zone->half_side };
			int side = get_side(zone);
			WORD *buffer = malloc(side * side * sizeof(WORD));	// TODO: prealloc
			if (!buffer) {
				fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
				abort = 1;
				continue;
			}

			// read the zone data to buffer
			copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
#ifdef DEBUG_MPP1
			save_buffer_tmp(frame, zone_idx, buffer, side);
#endif
			if (!zone->mpregparam) {
				zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
				if (!zone->mpregparam) {
					fprintf(stderr, "memory allocation failed\n");
					abort = 1;
					continue;
				}
			}

			int regretval; 
			/*if (args->using_homography) {
				if (!zone->mpregparam[frame].transform)
					zone->mpregparam[frame].transform = malloc(sizeof(Homography));
				regretval = ecc_find_transform_buf(references[zone_idx], buffer,
						side, fit.maxi, &reg_param,
						zone->mpregparam[frame].transform);
				if (regretval) {
					free(zone->mpregparam[frame].transform);
					zone->mpregparam[frame].transform = NULL;
					// TODO can we use the downsampled results? we need a backup
				}
			} else {*/
				regretval = ecc_find_translation_buf(references[zone_idx], buffer,
						side, fit.maxi, &reg_param);
				/*if (zone->mpregparam[frame].transform)
					free(zone->mpregparam[frame].transform);
				zone->mpregparam[frame].transform = NULL;
			}*/
			if (regretval) {
				fprintf(stdout, "ECC alignment failed for full def zone %d of frame %d\n", zone_idx, frame);
				// setting the quality to -1 will remove it from stacking, but it
				// may trigger the sequence quality reanalysis
				zone->mpregparam[frame].quality = -1.0;
			} else {
				// in ecc registration + sum stacking it's - and -
				zone->mpregparam[frame].x = reg_param.dx + regparam[frame].shiftx;
				zone->mpregparam[frame].y = -reg_param.dy + regparam[frame].shifty;
				fprintf(stdout, "frame %d, zone %d local shifts: %f,%f\n",
						frame, zone_idx, reg_param.dx, reg_param.dy);
			}

			free(buffer);
		}
		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}
	return abort;
}


