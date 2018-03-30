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
#include "core/proto.h"		// FITS functions
#include "stacking.h"
#include "algos/demosaicing.h"

/* TODO:
 * make sure nb_layers is 1 for the sequence on start
 * make sure pattern is not BAYER_FILTER_NONE and that all other types are handled
 * make sure there is a registration layer, i.e. not -1 and data present
 * pass pattern in the args
 * open the result as rgb while the input was b&w
 * warn about missing data
 * interpolate missing data
 * change the main atomic condition to a lock on the channel
 */

struct bayer_stacking_data {
	unsigned long *stack[3];// the new image's channels
	double exposure;	// sum of the exposures
	int reglayer;		// layer used for registration data
	sensor_pattern pattern;	// bayer matrix pattern
};

static int bayer_stacking_prepare_hook(struct generic_seq_args *args) {
	struct bayer_stacking_data *ssdata = args->user;
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	ssdata->stack[0] = calloc(3 * nbdata, sizeof(unsigned long));
	if (ssdata->stack[0] == NULL){
		printf("Stacking: memory allocation failure\n");
		return -1;
	}
	ssdata->stack[1] = ssdata->stack[0] + nbdata;	// index of green layer in stack[0]
	ssdata->stack[2] = ssdata->stack[0] + nbdata*2;	// index of blue layer in stack[0]

	ssdata->exposure = 0.0;
	ssdata->pattern = BAYER_FILTER_RGGB;	// TODO: pass the actual value
	return 0;
}

static int bayer_stacking_image_hook(struct generic_seq_args *args, int i, fits *fit, rectangle *_) {
	struct bayer_stacking_data *ssdata = args->user;
	int shiftx, shifty, nx, ny, x, y, ii;
	int pixel = 0;	// pixel index in result image
	int first = 1;

#ifdef _OPENMP
#pragma omp atomic
#endif
	ssdata->exposure += fit->exposure;
	
	if (ssdata->reglayer != -1 && args->seq->regparam[ssdata->reglayer]) {
		shiftx = roundf_to_int(args->seq->regparam[ssdata->reglayer][i].shiftx * args->seq->upscale_at_stacking);
		shifty = roundf_to_int(args->seq->regparam[ssdata->reglayer][i].shifty * args->seq->upscale_at_stacking);
	} else {
		shiftx = 0;
		shifty = 0;
	}

	for (y=0; y < fit->ry; ++y){
		for (x=0; x < fit->rx; ++x){
			nx = x - shiftx;
			ny = y - shifty;
			if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
				// we have data for this pixel
				int color = get_bayer_color(ssdata->pattern, nx, ny);
				ii = ny * fit->rx + nx;		// index in source image
				if (first) {
					fprintf(stdout, "bayer %d - shiftx %d, shifty %d, color %d, ii %d\n", i, shiftx, shifty, color, ii);
					first = 0;
				}
#ifdef _OPENMP
#pragma omp atomic
#endif
				ssdata->stack[color][pixel] += fit->pdata[0][ii];
			}
			++pixel;
		}
	}
	return 0;
}

// convert the result and store it into gfit
static int bayer_stacking_finalize_hook(struct generic_seq_args *args) {
	struct bayer_stacking_data *ssdata = args->user;
	unsigned long max = 0L;	// max value of the image's channels
	unsigned int i, nbdata;
	int layer;

	nbdata = args->seq->ry * args->seq->rx * 3;
	// find the max first
	#pragma omp parallel for reduction(max:max)
	for (i=0; i < nbdata; ++i)
		if (ssdata->stack[0][i] > max)
			max = ssdata->stack[0][i];

	clearfits(&gfit);
	fits *fit = &gfit;
	if (new_fit_image(&fit, args->seq->rx, args->seq->ry, 3))
		return -1;
	gfit.hi = round_to_WORD(max);
	gfit.exposure = ssdata->exposure;

	double ratio = 1.0;
	if (max > USHRT_MAX)
		ratio = USHRT_MAX_DOUBLE / (double)max;

	nbdata = args->seq->ry * args->seq->rx;
	for (layer=0; layer<3; ++layer){
		unsigned long* from = ssdata->stack[layer];
		WORD *to = gfit.pdata[layer];
		for (i=0; i < nbdata; ++i) {
			if (ratio == 1.0)
				*to++ = round_to_WORD(*from++);
			else	*to++ = round_to_WORD((double)(*from++) * ratio);
		}
	}

	free(ssdata->stack[0]);
	return 0;
}

int stack_bayer_generic(struct stacking_args *stackargs) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = stackargs->seq;
	args->partial_image = FALSE;
	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = bayer_stacking_prepare_hook;
	args->image_hook = bayer_stacking_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = bayer_stacking_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Bayer stacking");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct bayer_stacking_data *ssdata = malloc(sizeof(struct bayer_stacking_data));
	ssdata->reglayer = stackargs->reglayer;
	args->user = ssdata;

	//start_in_new_thread(generic_sequence_worker, args);
	generic_sequence_worker(args);
	return args->retval;
}


