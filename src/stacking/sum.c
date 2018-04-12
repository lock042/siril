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
#include "gui/progress_and_log.h"

struct sum_stacking_data {
	unsigned long *sum[3];	// the new image's channels
	unsigned int *contributions[3];	// number of actual values per pixel
	double exposure;	// sum of the exposures
	int reglayer;		// layer used for registration data
	int debayer;		// debayer drizzle
	sensor_pattern pattern;	// color filter array pattern
};

static int sum_stacking_prepare_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	int nb_layers;
	unsigned int nbdata = args->seq->ry * args->seq->rx;

	nb_layers = ssdata->debayer ? 3 : args->seq->nb_layers;
	ssdata->sum[0] = calloc(nb_layers, nbdata * sizeof(unsigned long));

	if (ssdata->sum[0] == NULL) {
		printf("Stacking: memory allocation failure\n");
		return -1;
	}
	if (ssdata->debayer || args->seq->nb_layers == 3) {
		ssdata->sum[1] = ssdata->sum[0] + nbdata; // index of green layer in sum[0]
		ssdata->sum[2] = ssdata->sum[1] + nbdata; // index of blue layer in sum[0]
	} else {
		ssdata->sum[1] = NULL;
		ssdata->sum[2] = NULL;
	}
	if (ssdata->debayer) {
		ssdata->contributions[0] = calloc(3, nbdata * sizeof(unsigned int));
		ssdata->contributions[1] = ssdata->contributions[0] + nbdata;
		ssdata->contributions[2] = ssdata->contributions[1] + nbdata;
	}

	ssdata->exposure = 0.0;
	return 0;
}

static int sum_stacking_image_hook(struct generic_seq_args *args, int i, fits *fit, rectangle *_) {
	struct sum_stacking_data *ssdata = args->user;
	int shiftx, shifty, nx, ny, x, y, ii, layer;
	int pixel = 0;	// index in sum[0]

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
				ii = ny * fit->rx + nx;		// index in source image
				if (ssdata->debayer) {
					int color = get_bayer_color(ssdata->pattern, nx, ny);
#ifdef _OPENMP
#pragma omp atomic
#endif
					ssdata->sum[color][pixel] += fit->pdata[0][ii];
#ifdef _OPENMP
#pragma omp atomic
#endif
					ssdata->contributions[color][pixel]++;
				} else {
					for(layer=0; layer<args->seq->nb_layers; ++layer) {
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
	return 0;
}

// convert the result and store it into gfit
static int sum_stacking_finalize_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	unsigned int i, nbdata;
	int nb_layers;

	nb_layers = ssdata->debayer ? 3 : args->seq->nb_layers;
	nbdata = args->seq->ry * args->seq->rx * nb_layers;

	clearfits(&gfit);
	fits *fit = &gfit;
	if (new_fit_image(&fit, args->seq->rx, args->seq->ry, nb_layers))
		return -1;
	gfit.exposure = ssdata->exposure;
	
	if (ssdata->debayer) {
		/* normalize each pixel with the number of contributions */
		double *normalized_image = malloc(nbdata * sizeof(double));
		if (!normalized_image) {
			siril_log_color_message(_("Stacking result normalization failed because it reached the memory limit, using a less accurate version\n"), "red");
			for (i=0; i < nbdata; ++i)
				if (ssdata->contributions[0][i] > 1)
					ssdata->sum[0][i] /= ssdata->contributions[0][i];
		} else {
			double dmax = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (i=0; i < nbdata; ++i)
				if (ssdata->contributions[0][i] > 1)
					normalized_image[i] = (double)ssdata->sum[0][i] / (double)ssdata->contributions[0][i];

			// find the max for normalization
#ifdef _OPENMP
#pragma omp parallel for reduction(max:dmax) schedule(static)
#endif
			for (i=0; i < nbdata; ++i)
				if (normalized_image[i] > dmax)
					dmax = normalized_image[i];

			double ratio = USHRT_MAX_DOUBLE / dmax;

			WORD *to = gfit.data;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (i=0; i < nbdata; ++i)
				to[i] = round_to_WORD(normalized_image[i] * ratio);

			free(normalized_image);
		}

		free(ssdata->contributions[0]);
		gfit.hi = 65535;
	} else {
		unsigned long max = 0L;	// max value of the image's channels
		// find the max for normalization
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max) schedule(static)
#endif
		for (i=0; i < nbdata; ++i)
			if (ssdata->sum[0][i] > max)
				max = ssdata->sum[0][i];
		double ratio = 1.0;
		if (max > USHRT_MAX)
			ratio = USHRT_MAX_DOUBLE / (double)max;

		unsigned long* from = ssdata->sum[0];
		WORD *to = gfit.data;
		if (ratio == 1.0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (i=0; i < nbdata; ++i)
				*to++ = round_to_WORD(*from++);
			gfit.hi = round_to_WORD(max);	// wrong!
		} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (i=0; i < nbdata; ++i)
				*to++ = round_to_WORD((double)(*from++) * ratio);
			gfit.hi = 65535;	// wrong!
		}
	}

	free(ssdata->sum[0]);
	return 0;
}

int stack_summing_generic(struct stacking_args *stackargs) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = stackargs->seq;
	args->partial_image = FALSE;
	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = sum_stacking_prepare_hook;
	args->image_hook = sum_stacking_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = sum_stacking_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Sum stacking");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct sum_stacking_data *ssdata = malloc(sizeof(struct sum_stacking_data));
	ssdata->reglayer = stackargs->reglayer;
	ssdata->debayer = stackargs->debayer;
	ssdata->pattern = stackargs->pattern;
	args->user = ssdata;

	//start_in_new_thread(generic_sequence_worker, args);
	generic_sequence_worker(args);
	return args->retval;
}

