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

#include <float.h>

#include "core/siril.h"
#include "core/processing.h"
#include "registration.h"
#include "algos/quality.h"
#include "gui/progress_and_log.h"
#include "stacking/stacking.h"
#include "core/proto.h"	// writeseqfile

struct regcog_data {
	struct registration_args *regargs;
	regdata *current_regdata;
};

static int regcog_prepare_hook(struct generic_seq_args *args) {
	struct regcog_data *rcdata = args->user;
	struct registration_args *regargs = rcdata->regargs;

	if (args->seq->regparam[regargs->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		rcdata->current_regdata = args->seq->regparam[regargs->layer];
		/* we reset all values as we may register different images */
		memset(rcdata->current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		rcdata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (rcdata->current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
	}
	return 0;
}

static int regcog_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct regcog_data *rcdata = args->user;
	struct registration_args *regargs = rcdata->regargs;
	double shiftx, shifty, quality;
	int retval;

	// Compute barycentre for the image
	retval = FindCentre(fit, regargs->layer, &shiftx, &shifty);
	if (retval) return retval;
	rcdata->current_regdata[in_index].shiftx = (float)shiftx;
	rcdata->current_regdata[in_index].shifty = (float)shifty;

	// Compute quality for the image
	// Warning: QualityEstimate modifies fit on this layer
	quality = QualityEstimate(fit, regargs->layer);
	rcdata->current_regdata[in_index].quality = quality;
	return retval;
}

static int regcog_finalize_hook(struct generic_seq_args *args) {
	struct regcog_data *rcdata = args->user;
	struct registration_args *regargs = rcdata->regargs;
	if (!args->retval) {
		args->seq->regparam[regargs->layer] = rcdata->current_regdata;

		/* find best frame and normalize quality */
		int i, q_index;
		double q_max = 0, q_min = DBL_MAX;
		for (i = 0; i < args->seq->number; i++) {
			double qual = rcdata->current_regdata[i].quality;
			if (qual <= 0.0)
				continue;
			if (qual > q_max) {
				q_max = qual;
				q_index = i;
			}
			q_min = min(q_min, qual);
		}
		normalizeQualityData(regargs, q_min, q_max);

		/* get shifts relative to the reference image */
		double ref_x = rcdata->current_regdata[q_index].shiftx;
		double ref_y = rcdata->current_regdata[q_index].shifty;
		for (i = 0; i < args->seq->number; i++) {
			if (rcdata->current_regdata[i].quality < 0.0) continue;

			rcdata->current_regdata[i].shiftx = ref_x - rcdata->current_regdata[i].shiftx;
			/* for Y, FITS are upside-down */
			if (args->seq->type == SEQ_REGULAR)
				rcdata->current_regdata[i].shifty = ref_y + rcdata->current_regdata[i].shifty;
			else rcdata->current_regdata[i].shifty = ref_y - rcdata->current_regdata[i].shifty;
		}
		siril_log_message(_("Registration finished.\n"));
		siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);
		args->seq->reference_image = q_index;
		writeseqfile(args->seq);
	}
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	free(rcdata);
	args->user = NULL;
	return 0;
}

int register_cog(struct registration_args *regargs) {
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
	args->layer = regargs->layer;

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

	args->prepare_hook = regcog_prepare_hook;
	args->image_hook = regcog_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = regcog_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Centre of gravity registration");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;	// ?
	args->parallel = TRUE;

	struct regcog_data *rcdata = calloc(1, sizeof(struct regcog_data));
	if (!rcdata) {
		free(args);
		return -1;
	}
	rcdata->regargs = regargs;
	args->user = rcdata;

	generic_sequence_worker(args);
	return args->retval;
}

