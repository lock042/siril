/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "stacking.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"

struct sum_stacking_data {
	guint64 *sum[3];	// the new image's channels
	double *fsum[3];	// the new image's channels, for float input image
	GList *list_date; // list of dates of every FITS file
	double livetime;	// sum of the exposures
	int reglayer;		// layer used for registration data
	int ref_image;		// reference image index in the stacked sequence
	gboolean input_32bits;	// input is a sequence of 32-bit float images
	gboolean output_32bits;	// output a 32-bit float image instead of the default ushort
	fits result;
};

static int sum_stacking_prepare_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	size_t nbdata = args->seq->ry * args->seq->rx;

	if (ssdata->input_32bits) {
		ssdata->fsum[0] = calloc(nbdata, sizeof(double) * args->seq->nb_layers);
		if (ssdata->fsum[0] == NULL){
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		if(args->seq->nb_layers == 3){
			ssdata->fsum[1] = ssdata->fsum[0] + nbdata;
			ssdata->fsum[2] = ssdata->fsum[0] + nbdata*2;
		} else {
			ssdata->fsum[1] = NULL;
			ssdata->fsum[2] = NULL;
		}
		ssdata->sum[0] = NULL;
	} else {
		ssdata->sum[0] = calloc(nbdata, sizeof(guint64) * args->seq->nb_layers);
		if (ssdata->sum[0] == NULL){
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		if(args->seq->nb_layers == 3){
			ssdata->sum[1] = ssdata->sum[0] + nbdata;	// index of green layer in sum[0]
			ssdata->sum[2] = ssdata->sum[0] + nbdata*2;	// index of blue layer in sum[0]
		} else {
			ssdata->sum[1] = NULL;
			ssdata->sum[2] = NULL;
		}
		ssdata->fsum[0] = NULL;
	}

	ssdata->livetime = 0.0;
	ssdata->list_date = NULL;
	return ST_OK;
}

static int sum_stacking_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct sum_stacking_data *ssdata = args->user;
	int shiftx = 0, shifty = 0, nx, ny, x, y, layer;
	size_t ii, pixel = 0;	// index in sum[0]
	/* we get some metadata at the same time: date, exposure ... */

#ifdef _OPENMP
#pragma omp atomic
#endif
	ssdata->livetime += fit->exposure;

	if (fit->date_obs) {
		GDateTime *date = g_date_time_ref(fit->date_obs);
		ssdata->list_date = g_list_prepend(ssdata->list_date, new_date_item(date, fit->exposure));
	}

	if (ssdata->reglayer != -1 && args->seq->regparam[ssdata->reglayer]) {
		double scale = args->seq->upscale_at_stacking;
		double dx, dy;
		translation_from_H(args->seq->regparam[ssdata->reglayer][i].H, &dx, &dy);
		shiftx = round_to_int(dx * scale);
		shifty = round_to_int(dy * scale);
	}

	for (y = 0; y < fit->ry; ++y) {
		for (x = 0; x < fit->rx; ++x) {
			nx = x - shiftx;
			ny = y - shifty;
			if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
				// we have data for this pixel
				ii = ny * fit->rx + nx;		// index in source image
				if (ii < fit->rx * fit->ry) {
					for (layer = 0; layer < args->seq->nb_layers; ++layer) {
						if (ssdata->input_32bits) {
#ifdef _OPENMP
#pragma omp atomic
#endif
							ssdata->fsum[layer][pixel] += (double)fit->fpdata[layer][ii];
						}
						else
#ifdef _OPENMP
#pragma omp atomic
#endif
							ssdata->sum[layer][pixel] += fit->pdata[layer][ii];
					}
				}
			}
			++pixel;
		}
	}
	return ST_OK;
}

// convert the result and store it
static int sum_stacking_finalize_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	guint64 max = 0L;	// max value of the image's channels
	double fmax = 0.0;
	size_t i, nbdata;
	int layer;

	if (args->retval) {
		if (ssdata->sum[0]) free(ssdata->sum[0]);
		if (ssdata->fsum[0]) free(ssdata->fsum[0]);
		args->user = NULL;
		return args->retval;
	}

	nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
	// find the max first
	if (ssdata->input_32bits) {
#ifdef _OPENMP
#pragma omp parallel for reduction(max:fmax)
#endif
		for (i = 0; i < nbdata; ++i) {
			if (ssdata->fsum[0][i] > fmax)
				fmax = ssdata->fsum[0][i];
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i = 0; i < nbdata; ++i)
			if (ssdata->sum[0][i] > max)
				max = ssdata->sum[0][i];
	}

	fits *fit = &ssdata->result;
	if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers, ssdata->output_32bits ? DATA_FLOAT : DATA_USHORT))
		return ST_GENERIC_ERROR;

	/* We copy metadata from reference to the final fit */
	int ref = ssdata->ref_image;
	if (args->seq->type == SEQ_REGULAR) {
		if (!seq_open_image(args->seq, ref)) {
			import_metadata_from_fitsfile(args->seq->fptr[ref], fit);
			seq_close_image(args->seq, ref);
		}
		fit->livetime = ssdata->livetime;
	} else if (args->seq->type == SEQ_FITSEQ) {
		if (!fitseq_set_current_frame(args->seq->fitseq_file, ref))
			import_metadata_from_fitsfile(args->seq->fitseq_file->fptr, fit);
		fit->livetime = ssdata->livetime;
	} else if (args->seq->type == SEQ_SER) {
		import_metadata_from_serfile(args->seq->ser_file, fit);
		fit->livetime = fit->exposure * args->nb_filtered_images; // ssdata->livetime is null for ser as fit has no exposure data
	}

	fit->stackcnt = args->nb_filtered_images;
	nbdata = args->seq->ry * args->seq->rx;
	compute_date_time_keywords(ssdata->list_date, fit);
	g_list_free_full(ssdata->list_date, (GDestroyNotify) free_list_date);

	if (ssdata->output_32bits) {
		if (ssdata->input_32bits) {
			double ratio = 1.0 / fmax;
			for (layer=0; layer<args->seq->nb_layers; ++layer){
				double *from = ssdata->fsum[layer];
				float *to = fit->fpdata[layer];
				for (i=0; i < nbdata; ++i) {
					*to++ = (float)((double)(*from++) * ratio);
				}
			}
		} else {
			double ratio = 1.0 / (max == 0 ? 1 : (double)max);
			for (layer=0; layer<args->seq->nb_layers; ++layer){
				guint64 *from = ssdata->sum[layer];
				float *to = fit->fpdata[layer];
				for (i=0; i < nbdata; ++i) {
					*to++ = (float)((double)(*from++) * ratio);
				}
			}
		}
	} else {
		double ratio = 1.0;
		if (max > USHRT_MAX) {
			ratio = USHRT_MAX_DOUBLE / (double)max;
			siril_log_color_message(_("Reducing the stacking output to a 16-bit image will result in precision loss\n"), "salmon");
		}

		for (layer=0; layer<args->seq->nb_layers; ++layer){
			guint64 *from = ssdata->sum[layer];
			WORD *to = fit->pdata[layer];
			for (i=0; i < nbdata; ++i) {
				if (ratio == 1.0)
					*to++ = round_to_WORD(*from++);
				else *to++ = round_to_WORD((double)(*from++) * ratio);
			}
		}
	}

	if (ssdata->sum[0]) free(ssdata->sum[0]);
	if (ssdata->fsum[0]) free(ssdata->fsum[0]);
	args->user = NULL;

	return ST_OK;
}

int stack_summing_generic(struct stacking_args *stackargs) {
	struct generic_seq_args *args = create_default_seqargs(stackargs->seq);
	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = sum_stacking_prepare_hook;
	args->image_hook = sum_stacking_image_hook;
	args->finalize_hook = sum_stacking_finalize_hook;
	args->description = _("Sum stacking");
	args->already_in_a_thread = TRUE;

	struct sum_stacking_data *ssdata = calloc(1, sizeof(struct sum_stacking_data));
	ssdata->reglayer = stackargs->reglayer;
	ssdata->ref_image = stackargs->ref_image;
	assert(ssdata->ref_image >= 0 && ssdata->ref_image < args->seq->number);
	ssdata->input_32bits = get_data_type(args->seq->bitpix) == DATA_FLOAT;
	ssdata->output_32bits = stackargs->use_32bit_output;
	if (ssdata->input_32bits)
		assert(ssdata->output_32bits);
	args->user = ssdata;

	generic_sequence_worker(args);
	memcpy(&stackargs->result, &ssdata->result, sizeof(fits));
	free(ssdata);
	return args->retval;
}

