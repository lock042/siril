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
#include "algos/planetary_caching.h"
#include "laplacian_quality.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"
#include "core/proto.h"	// writeseqfile

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
		/* we reset all values as we may register different images */
		memset(ldata->current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		ldata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (ldata->current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
	}
	return 0;
}

int lapl_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct lapl_data *ldata = args->user;
	
	WORD * gaussian_data = get_gaussian_data_for_image(in_index, fit, ldata->cache);
	WORD * laplacian_data = get_laplacian_data_for_image(in_index,
			gaussian_data, fit->rx, fit->ry, ldata->cache);

	// optional: apply downsampling?

	ldata->current_regdata[in_index].quality = variance(laplacian_data, fit->rx * fit->ry);
	//fprintf(stdout, "variance for %d: %g\n", in_index, ldata->current_regdata[in_index].quality);

	free(laplacian_data);
	free(gaussian_data);
	return 0;
}

int lapl_finalize_hook(struct generic_seq_args *args) {
	struct lapl_data *ldata = args->user;
	finalize_caching(ldata->cache);
	free(ldata->cache);
	if (!args->retval) {
		args->seq->regparam[0] = ldata->current_regdata;
		/* find best frame and normalize quality */

		int q_index = normalize_quality_data(ldata->current_regdata, args->seq->number);
		siril_log_message(_("Laplacian quality finished.\n"));
		siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);
		args->seq->reference_image = q_index;
		writeseqfile(args->seq);
	} 
	else {
		siril_log_message(_("Laplacian quality aborted.\n"));
	}
	free(ldata);
	args->user = NULL;
	return 0;
}

double variance(WORD *set, int size) {
	int i;
	unsigned long sum = 0ul, squared_sum = 0ul;
	// could be parallelized with sum reductions, is it worth it?
	for (i = 0; i < size; i++) {
		sum += set[i];
		squared_sum += set[i] * set[i];
	}

	double mean = sum / size;
	return (double)squared_sum / (mean * mean);
}

