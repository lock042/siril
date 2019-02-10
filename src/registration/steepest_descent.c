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

#include "core/siril.h"
#include "core/processing.h"
#include "registration.h"
#include "algos/quality.h"
#include "gui/progress_and_log.h"
#include "stacking/stacking.h"
#include "core/proto.h"	// writeseqfile
#include "algos/planetary_caching.h"
#include "algos/laplacian_quality.h"

struct regsd_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	struct planetary_cache *cache;
	// TODO: add reference frame data
};

static int regsd_prepare_hook(struct generic_seq_args *args) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	rsdata->cache = calloc(1, sizeof(struct planetary_cache));
	rsdata->cache->use_cached = regargs->use_caching;
	rsdata->cache->cache_data = regargs->use_caching;
	init_caching(args->seq->seqname, rsdata->cache, regargs->kernel_size);

	/* TODO: check that layer = 0 is alright here */
	if (args->seq->regparam[0]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		rsdata->current_regdata = args->seq->regparam[0];
		/* we reset all values as we may register different images */
		memset(rsdata->current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		rsdata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (rsdata->current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
	}
	return 0;
}

static int regsd_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	
	WORD * gaussian_data = get_gaussian_data_for_image(in_index, fit, rsdata->cache);

	// compute shift using the steepest descent algorithm

	return 0;
}

static int regsd_finalize_hook(struct generic_seq_args *args) {
	struct regsd_data *rsdata = args->user;
	struct registration_args *regargs = rsdata->regargs;
	finalize_caching(rsdata->cache);
	free(rsdata->cache);
	if (!args->retval) {
		args->seq->regparam[0] = rsdata->current_regdata;
		/* find best frame and normalize quality */

		/* get shifts relative to the reference image */

		siril_log_message(_("Registration finished.\n"));
		//siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);
		//args->seq->reference_image = q_index;
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
	//args->layer = regargs->layer;
	args->layer = 0;

	/* check if it's already computed before doing it for nothing */
	if (args->seq->regparam[args->layer]) {
		int i;
		for (i = 0; i < args->seq->number; i++) {
			if (args->filtering_criterion(args->seq, args->layer, i, 0) &&
					args->seq->regparam[args->layer][i].quality < 0.0)
				break;
			// TODO: check shifts too ^
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
	args->description = _("Image quality computation");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;	// ?
	args->parallel = TRUE;

	struct lapl_data *ldata = malloc(sizeof(struct lapl_data));
	if (!ldata) {
		free(args);
		return -1;
	}
	ldata->use_caching = regargs->use_caching;
	ldata->kernel_size = regargs->kernel_size;
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
	args->stop_on_error = TRUE;
	args->description = _("Steepest descent registration");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;	// ?
	args->parallel = TRUE;

	struct regsd_data *rsdata = calloc(1, sizeof(struct regsd_data));
	if (!rsdata) {
		free(args);
		return -1;
	}
	rsdata->regargs = regargs;
	args->user = rsdata;

	generic_sequence_worker(args);
	return args->retval;
}

