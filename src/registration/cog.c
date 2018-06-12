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
	double shiftx, shifty;
	int retval;

	retval = FindCentre(fit, &shiftx, &shifty);
	rcdata->current_regdata[in_index].shiftx = (float)shiftx;
	rcdata->current_regdata[in_index].shifty = (float)shifty;
	return retval;
}

static int regcog_finalize_hook(struct generic_seq_args *args) {
	struct regcog_data *rcdata = args->user;
	struct registration_args *regargs = rcdata->regargs;
	if (!args->retval) {
		args->seq->regparam[regargs->layer] = rcdata->current_regdata;
		siril_log_message(_("Registration finished.\n"));
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
		args->filtering_criterion = seq_filter_all;
		args->nb_filtered_images = regargs->seq->number;
	} else {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
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

