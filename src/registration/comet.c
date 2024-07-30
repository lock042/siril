/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
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
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "registration.h"

pointf compute_velocity(GDateTime *t1, GDateTime *t2, point d1, point d2) {
	pointf delta_d, px_per_hour = { 0.f, 0.f };

	if (t1 && t2) {
		float delta_t = g_date_time_difference(t2, t1);
		delta_d.x = (float)(d2.x - d1.x);
		delta_d.y = (float)(d2.y - d1.y);

		px_per_hour.x = delta_d.x / delta_t * 3600000000.f;
		px_per_hour.y = delta_d.y / delta_t * 3600000000.f;
	}

	return px_per_hour;
}

int get_comet_shift(GDateTime *ref, GDateTime *img, pointf px_per_hour, pointf *reg) {
	if (img && ref) {
		float delta_t = (float) g_date_time_difference(img, ref);
		delta_t /= 3600000000.f;
		reg->x = delta_t * px_per_hour.x;
		reg->y = delta_t * px_per_hour.y;
	}
	return 0;
}

/***** generic moving object registration *****/

struct comet_align_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	GDateTime *reference_date;
};

static int comet_align_prepare_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;
	int ref_image;
	fits ref = { 0 };

	if (args->seq->regparam[regargs->layer]) {
		cadata->current_regdata = args->seq->regparam[regargs->layer];
	} else {
		cadata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (cadata->current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return -2;
		}
		args->seq->regparam[regargs->layer] = cadata->current_regdata;
	}

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	if (seq_read_frame_metadata(args->seq, ref_image, &ref)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(cadata->current_regdata);
		return 1;
	}
	cadata->reference_date = g_date_time_ref(ref.keywords.date_obs);
	clearfits(&ref);

	return 0;
}

static int comet_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;
	pointf reg = { 0.f, 0.f };

	get_comet_shift(cadata->reference_date, fit->keywords.date_obs, regargs->velocity, &reg);

	/* get_comet_shift does not care about orientation of image */
	cum_shifts(args->seq, in_index, regargs->layer, -reg.x, -reg.y);
	return 0;
}

static int comet_align_finalize_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;

	if (args->retval) {
		free(args->seq->regparam[regargs->layer]);
		args->seq->regparam[regargs->layer] = NULL;
	}

	if (cadata->reference_date)
		g_date_time_unref(cadata->reference_date);

	free(cadata);
	args->user = NULL;
	return 0;
}

int register_comet(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	/* we don't need to read image data, for simplicity we just read one
	 * pixel from it, making sure the header is read */
	args->partial_image = TRUE;
	args->area.x = 0; args->area.y = 0;
	args->area.w = 1; args->area.h = 1;
	args->layer_for_partial = 0;
	args->get_photometry_data_for_partial = TRUE;

	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->prepare_hook = comet_align_prepare_hook;
	args->image_hook = comet_align_image_hook;
	args->finalize_hook = comet_align_finalize_hook;
	args->description = _("Moving object registration");
	args->already_in_a_thread = TRUE;

	struct comet_align_data *cadata = calloc(1, sizeof(struct comet_align_data));
	if (!cadata) {
		free(args);
		return -1;
	}
	cadata->regargs = regargs;
	args->user = cadata;

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}
