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

/* This stacking method is not working well and is not currently used.
 *
 * It is a multi-point method, uses the global quality sorting of frames for
 * point selection, and computes a weighted sum shift for pixels outside
 * stacking zones based on neighbours' shifts.
 */

struct weighted_AP {
	int zone_index;
	double distance;
	// double direction; someday
};

static void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP);
//static void filter_closest_list_owned(struct weighted_AP *list_for_this_point,
//		int max_AP, double owned_distance);


/* this function is a stacking sum that takes shifts from the multipoint registration
 * instead of taking them from the global registration. If it works well, it will have to
 * be exploded as a generic processing function. */
static int the_global_multipoint_barycentric_sum_stacking(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, abort = 0;
	regdata *regparam = args->seq->regparam[args->layer];
	struct weighted_AP *closest_zones_map;	// list of nb_closest_AP AP (fixed) for each pixel

	nb_zones = get_number_of_zones();
	if (nb_zones < args->nb_closest_AP)
		args->nb_closest_AP = nb_zones;
	closest_zones_map = malloc(args->seq->rx * args->seq->ry * args->nb_closest_AP * sizeof(struct weighted_AP));
	if (!closest_zones_map) {
		fprintf(stderr, "Stacking: memory allocation failure for zone mapping\n");
		return -1;
	}
	siril_log_message("Using %d closest zones for multipoint stacking refinement\n", args->nb_closest_AP);

	/* precompute the closest zones from each pixel */
	int x, y;
	for (y = 0; y < args->seq->ry; y++) {
		for (x = 0; x < args->seq->rx; x++) {
			int ap;
			struct weighted_AP *list_for_this_pixel = closest_zones_map + (x + y * args->seq->rx) * args->nb_closest_AP;
			for (ap = 0; ap < args->nb_closest_AP; ap++)	// init the struct
				list_for_this_pixel[ap].distance = -1.0;

			for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
				stacking_zone *zone = &com.stacking_zones[zone_idx];
				// compute the distance from the centre of the pixel
				double xdist = zone->centre.x - x + 0.5;
				double ydist = zone->centre.y - y + 0.5;
				double distance = sqrt(xdist * xdist + ydist * ydist);
				if (distance < args->max_distance)
					check_closest_list(list_for_this_pixel, distance, zone_idx, args->nb_closest_AP);
			}
			/* TODO: doesn't work, zone size has to be known, it has to be
			 * done in the main check */
			//filter_closest_list_owned(list_for_this_point, args->nb_closest_AP, args->own_distance_f * zone
		}
	}

	// init stacking data (copied from sum_stacking_prepare_hook)
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	unsigned long *sum[3];	// the new image's channels
	int *count[3];	// the new image's contributions count
	sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!sum[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		sum[1] = sum[0] + nbdata;	// index of green layer in sum[0]
		sum[2] = sum[0] + nbdata*2;	// index of blue layer in sum[0]
	} else {
		sum[1] = NULL;
		sum[2] = NULL;
	}

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		if (!args->filtering_criterion(args->seq, args->layer, frame, args->filtering_parameter))
			continue;

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "barycentre stacking for image %d\n", frame);

		/* for each image, we use the closest zones' shifts computed in the
		 * multi-point registration and weight them with their distance to get
		 * the coordinates of the pixel of the image to stack for each end-image
		 * pixel. The direction of each AP would also be important, the
		 * distribution of the AP having a role in the weighing. */
		int x, y;
		int pixel = 0;	// index in sum[0]
		for (y = 0; y < fit.ry; y++) {
			for (x = 0; x < fit.rx; x++) {
				struct weighted_AP *list_for_this_pixel = closest_zones_map + (x + y * args->seq->rx) * args->nb_closest_AP;
				int i;
				double total_weight = 0.0, sumx = 0.0, sumy = 0.0, shiftx, shifty;

				for (i = 0; i < args->nb_closest_AP; i++) {
					if (list_for_this_pixel[i].distance < 0.0)
						continue;
					double this_weight;
					if (list_for_this_pixel[i].distance <= 1.0)
						this_weight = args->max_distance;
					else this_weight = args->max_distance / list_for_this_pixel[i].distance;
					int zone_index = list_for_this_pixel[i].zone_index;
					struct ap_regdata *regparam_zone_image = 
						&com.stacking_zones[zone_index].mpregparam[frame];
					if (regparam_zone_image->quality < 0.0)
						continue;

					total_weight += this_weight;
					sumx += regparam_zone_image->x * this_weight;
					sumy += regparam_zone_image->y * this_weight;
				}
				// if zones are too far away, which will happen for pixels
				// away from the planet, just take the global shift
				if (total_weight == 0.0) {
					shiftx = regparam[frame].shiftx;
					shifty = regparam[frame].shifty;
				} else {
					shiftx = sumx / total_weight;
					shifty = sumy / total_weight;
				}

				// then do the regular sum stacking with shift
				int nx = round_to_int(x - shiftx);
				int ny = round_to_int(y - shifty);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					int ii = ny * fit.rx + nx;	// index in source image
					int layer;
					for (layer = 0; layer < args->seq->nb_layers; ++layer) {
#ifdef _OPENMP
//#pragma omp atomic
#endif
						sum[layer][pixel] += fit.pdata[layer][ii];
					}
				}
				++pixel;
			}
		}
		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	fprintf(stdout, "barycentre stacking ended, creating final image\n");

	if (!abort) {
		/* copying result into gfit and args->global_image
		 * code copied from sum_stacking_finalize_hook() */
		nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
		int i, layer;
		// find the max first
		unsigned long max = 0L;	// max value of the image's channels
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i=0; i < nbdata; ++i)
			if (sum[0][i] > max)
				max = sum[0][i];

		clearfits(&gfit);
		fits *fit = &gfit;
		if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers, NULL))
			return -1;
		gfit.hi = round_to_WORD(max);

		double ratio = 1.0;
		if (max > USHRT_MAX)
			ratio = USHRT_MAX_DOUBLE / (double)max;

		double norm_ratio = 1.0 / (double)max;
		args->global_image = malloc(sizeof(double) * nbdata);

		unsigned long* from = sum[0];
		WORD *to = gfit.data;
		for (i=0; i < nbdata; ++i) {
			if (ratio == 1.0)
				to[i] = round_to_WORD(from[i]);
			else	to[i] = round_to_WORD((double)(from[i]) * ratio);
			args->global_image[i] = from[i] * norm_ratio;
		}
	}

	free(sum[0]);

	return 0;
}

/* we build a list of max_AP of the closest alignment points for a pixel. The
 * list is passed, max_AP too, and we try a new candidate AP for the top max_AP.
 * The built list is not ordered.
 */
static void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP) {
	int i, max = 0, max_idx = 0, found_closer = 0;
	for (i = 0; i < max_AP; i++) {
		double ap_dist = list_for_this_point[i].distance;
		if (ap_dist >= 0) {
			// we search the max to replace it by the new
			if (ap_dist > max) {
				max = ap_dist;
				max_idx = i;
			}
			if (distance < ap_dist)
				found_closer = 1;
		}
		else {
			found_closer = 1;
			max_idx = i;
			break;
		}
	}
	if (found_closer) {
		list_for_this_point[max_idx].distance = distance;
		list_for_this_point[max_idx].zone_index = zone_idx;
	}
}

#if 0
/* when the list of closest AP is complete, we check for special case of pixel position
 * being very close to AP centres. Below this owned_distance threshold, other AP are not
 * considered relevant and are removed of the list. */
static void filter_closest_list_owned(struct weighted_AP *list_for_this_point,
		int max_AP, double owned_distance) {
	int i, unowned_idx = -1, owned = 0;
	for (i = 0; i < max_AP; i++) {
		double ap_dist = list_for_this_point[i].distance;
		if (ap_dist < 0) continue;
		if (ap_dist < owned_distance) {
			owned = 1;
			if (unowned_idx != -1) {
				list_for_this_point[unowned_idx].distance = ap_dist;
				list_for_this_point[unowned_idx].zone_index = list_for_this_point[i].zone_index;
				list_for_this_point[i].distance = -1.0;
				list_for_this_point[i].zone_index = -1;
				while (unowned_idx < i &&
						list_for_this_point[unowned_idx].distance < owned_distance);
				if (i == unowned_idx)
					unowned_idx = -1;
			}
		}
		else {
		       if (owned) {
			       list_for_this_point[i].distance = -1.0;
			       list_for_this_point[i].zone_index = -1;
		       }
		       if (unowned_idx == -1)
			       unowned_idx = i;
		}
	}
}
#endif


