/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "algos/siril_wcs.h"
#include "opencv/opencv.h"

struct sum_stacking_data {
	guint64 *sum[3];	// the new image's channels
	double *fsum[3];	// the new image's channels, for float input image
	double *fweight[3];	// the new image's weights for drizzle stacking
	GList *list_date; // list of dates of every FITS file
	double livetime;	// sum of the exposures
	int reglayer;		// layer used for registration data
	int ref_image;		// reference image index in the stacked sequence
	gboolean input_32bits;	// input is a sequence of 32-bit float images
	gboolean output_32bits;	// output a 32-bit float image instead of the default ushort
	gboolean maximize_framing; // outputs an image that encompasses all frames including registration
	gboolean upscale_at_stacking; // outputs an image twice the size of the original sequence
	int output_size[2]; // stacked image size
	int offset[2]; // reference offset
	fits result;
};

static int sum_stacking_prepare_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	size_t nbdata = ssdata->output_size[0] * ssdata->output_size[1];

	if (ssdata->input_32bits || args->seq->is_drizzle) {
		ssdata->fsum[0] = calloc(nbdata, sizeof(double) * args->seq->nb_layers);
		if (ssdata->fsum[0] == NULL){
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		if(args->seq->nb_layers == 3){
			ssdata->fsum[1] = ssdata->fsum[0] + nbdata;
			ssdata->fsum[2] = ssdata->fsum[0] + nbdata * 2;
		} else {
			ssdata->fsum[1] = NULL;
			ssdata->fsum[2] = NULL;
		}
		if (args->seq->is_drizzle) {
			ssdata->fweight[0] = calloc(nbdata, sizeof(double) * args->seq->nb_layers);
			if (ssdata->fweight[0] == NULL){
				PRINT_ALLOC_ERR;
				return ST_ALLOC_ERROR;
			}
			if(args->seq->nb_layers == 3){
				ssdata->fweight[1] = ssdata->fweight[0] + nbdata;	// index of green layer in fweight[0]
				ssdata->fweight[2] = ssdata->fweight[0] + nbdata * 2;	// index of blue layer in fweight[0]
			} else {
				ssdata->fweight[1] = NULL;
				ssdata->fweight[2] = NULL;
			}
		} else
			ssdata->fweight[0] = NULL;
		ssdata->sum[0] = NULL;
	} else {
		ssdata->sum[0] = calloc(nbdata, sizeof(guint64) * args->seq->nb_layers);
		if (ssdata->sum[0] == NULL){
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		if(args->seq->nb_layers == 3){
			ssdata->sum[1] = ssdata->sum[0] + nbdata;	// index of green layer in sum[0]
			ssdata->sum[2] = ssdata->sum[0] + nbdata * 2;	// index of blue layer in sum[0]
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
	size_t nbdata = ssdata->output_size[0] * ssdata->output_size[1];
	gboolean is_drizzle = args->seq->is_drizzle;
	float **weights = NULL;
	if (is_drizzle) {
		const gchar *drizzfile = get_sequence_cache_filename(args->seq, i, "drizztmp", "fit", NULL);
		int rx = (args->seq->is_variable) ? args->seq->imgparam[i].rx : args->seq->rx;
		int ry = (args->seq->is_variable) ? args->seq->imgparam[i].ry : args->seq->ry;
		rectangle drizz_area = { 0, 0, rx, ry};
		float *dweights = calloc(rx * ry, sizeof(float) * args->seq->nb_layers);
		int layer = (args->seq->nb_layers == 1) ? -1 : 4; // 4 means all layers, 1 means only the first layer
		if (read_drizz_fits_area(drizzfile, layer, &drizz_area, ry, dweights)) {
			siril_log_color_message(_("Error reading one of the drizzle weights areas (%d: %d %d %d %d)\n"), "red", i + 1,
			drizz_area.x, drizz_area.y, drizz_area.w, drizz_area.h);
			return ST_SEQUENCE_ERROR;
		}
		weights = malloc(sizeof(float *) * args->seq->nb_layers);
		for (int l = 0; l < args->seq->nb_layers; ++l) {
			weights[l] = dweights + (l * rx * ry);
		}
	}
	/* we get some metadata at the same time: date, exposure ... */

#ifdef _OPENMP
#pragma omp atomic
#endif
	ssdata->livetime += fit->keywords.exposure;

	if (fit->keywords.date_obs) {
		GDateTime *date = g_date_time_ref(fit->keywords.date_obs);
		ssdata->list_date = g_list_prepend(ssdata->list_date, new_date_item(date, fit->keywords.exposure));
	}

	if (ssdata->reglayer != -1 && args->seq->regparam[ssdata->reglayer]) {
		double scale = (ssdata->upscale_at_stacking) ? 2. : 1.;
		double dx, dy;
		translation_from_H(args->seq->regparam[ssdata->reglayer][i].H, &dx, &dy);
		dx *= scale;
		dy *= scale;
		dx -= ssdata->offset[0];
		dy -= ssdata->offset[1];
		if (ssdata->maximize_framing)
			dy -= (double)fit->ry;
		shiftx = round_to_int(dx);
		shifty = round_to_int(dy);
		siril_debug_print("img %d dx %d dy %d\n", o, shiftx, shifty);
	}

	for (y = 0; y < ssdata->output_size[1]; ++y) {
		ny = y - shifty;
		size_t rny = ny * fit->rx;
		for (x = 0; x < ssdata->output_size[0]; ++x) {
			nx = x - shiftx;
			if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
				// we have data for this pixel
				ii = rny + nx;		// index in source image
				if (ii < nbdata) {
					for (layer = 0; layer < args->seq->nb_layers; ++layer) {
						if (!is_drizzle) {
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
						} else {
							if (ssdata->input_32bits) {
#ifdef _OPENMP
#pragma omp atomic
#endif
								ssdata->fsum[layer][pixel] += (double)(fit->fpdata[layer][ii] * weights[layer][ii]);
							}
							else {
#ifdef _OPENMP
#pragma omp atomic
#endif
								ssdata->fsum[layer][pixel] += (double)(fit->pdata[layer][ii] * INV_USHRT_MAX_SINGLE * weights[layer][ii]);
							}
#ifdef _OPENMP
#pragma omp atomic
#endif
							ssdata->fweight[layer][pixel] += (double)weights[layer][ii]; // weight is always in the first layer
						}
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
		if (ssdata->fweight[0]) free(ssdata->fweight[0]);
		args->user = NULL;
		return args->retval;
	}

	nbdata = ssdata->output_size[0] * ssdata->output_size[1] * args->seq->nb_layers;

	if (!args->seq->is_drizzle) {
		// find the max first
		if (ssdata->input_32bits) {
#ifdef _OPENMP
#pragma omp parallel for reduction(max:fmax)
#endif
			for (i = 0; i < nbdata; ++i) {
				if (ssdata->fsum[0][i] > fmax)
					fmax = ssdata->fsum[0][i];
			}
			if (!fmax)
				return ST_GENERIC_ERROR;
		} else {
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
			for (i = 0; i < nbdata; ++i)
				if (ssdata->sum[0][i] > max)
					max = ssdata->sum[0][i];
			if (!max)
				return ST_GENERIC_ERROR;
		}
	}

	fits *fit = &ssdata->result;
	if (new_fit_image(&fit, ssdata->output_size[0], ssdata->output_size[1], args->seq->nb_layers, ssdata->output_32bits ? DATA_FLOAT : DATA_USHORT))
		return ST_GENERIC_ERROR;

	/* We copy metadata from reference to the final fit */
	int ref = ssdata->ref_image;
	if (args->seq->type == SEQ_REGULAR) {
		if (!seq_open_image(args->seq, ref)) {
			import_metadata_from_fitsfile(args->seq->fptr[ref], fit);
			seq_close_image(args->seq, ref);
		}
		fit->keywords.livetime = ssdata->livetime;
	} else if (args->seq->type == SEQ_FITSEQ) {
		if (!fitseq_set_current_frame(args->seq->fitseq_file, ref))
			import_metadata_from_fitsfile(args->seq->fitseq_file->fptr, fit);
		fit->keywords.livetime = ssdata->livetime;
	} else if (args->seq->type == SEQ_SER) {
		import_metadata_from_serfile(args->seq->ser_file, fit);
		fit->keywords.livetime = fit->keywords.exposure * args->nb_filtered_images; // ssdata->livetime is null for ser as fit has no exposure data
	}

	fit->keywords.stackcnt = args->nb_filtered_images;
	nbdata = ssdata->output_size[0] * ssdata->output_size[1];
	compute_date_time_keywords(ssdata->list_date, fit);
	g_list_free_full(ssdata->list_date, (GDestroyNotify) free_list_date);

	if (!args->seq->is_drizzle) {
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
	} else { // divide by drizzle weights
		for (layer=0; layer<args->seq->nb_layers; ++layer){
			double *from = ssdata->fsum[layer];
			double *fromw = ssdata->fweight[layer];
			float *to = fit->fpdata[layer];
			for (i = 0; i < nbdata; ++i) {
				if (ssdata->output_32bits) {
					if (*fromw > 0.0)
						*to++ = (float)((*from++) / (*fromw++));
					else {
						*to++ = 0.0f; // avoid division by zero
						++from;
						++fromw;
					}
				} else {
					if (*fromw > 0.0)
						*to++ = round_to_WORD((*from++) / (*fromw++));
					else {
						*to++ = 0; // avoid division by zero
						++from;
						++fromw;
					}
				}
			}
		}
	}

	if (ssdata->sum[0]) free(ssdata->sum[0]);
	if (ssdata->fsum[0]) free(ssdata->fsum[0]);
	if (ssdata->fweight[0]) free(ssdata->fweight[0]);
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
	gboolean update_wcs = FALSE;
	if (stackargs->maximize_framing) {
		compute_max_framing(stackargs, ssdata->output_size, ssdata->offset);
		ssdata->maximize_framing = TRUE;
		update_wcs = TRUE;
	} else {
		ssdata->output_size[0] = args->seq->rx;
		ssdata->output_size[1] = args->seq->ry;
		double dx = 0., dy = 0.;
		if (ssdata->reglayer >= 0) {
			translation_from_H(args->seq->regparam[ssdata->reglayer][ssdata->ref_image].H, &dx, &dy);
			update_wcs = TRUE;
		}
		ssdata->offset[0] = (int)dx;
		ssdata->offset[1] = (int)dy;
	}
	ssdata->upscale_at_stacking = stackargs->upscale_at_stacking;
	args->user = ssdata;

	generic_sequence_worker(args);
	memcpy(&stackargs->result, &ssdata->result, sizeof(fits));
	if (update_wcs && has_wcs(&ssdata->result)) {
		fits *result = &ssdata->result;
		Homography Hs = { 0 };
		cvGetEye(&Hs);
		double dx, dy;
		translation_from_H(args->seq->regparam[stackargs->reglayer][stackargs->ref_image].H, &dx, &dy);
		siril_debug_print("ref shift: %d %d\n", (int)dx, (int)dy);
		siril_debug_print("crpix: %.1f %.1f\n", result->keywords.wcslib->crpix[0], result->keywords.wcslib->crpix[1]);
		Hs.h02  = dx - ssdata->offset[0];
		Hs.h12 -= dy - ssdata->offset[1];
		int orig_rx = (args->seq->is_variable) ? args->seq->imgparam[args->seq->reference_image].rx : args->seq->rx;
		int orig_ry = (args->seq->is_variable) ? args->seq->imgparam[args->seq->reference_image].ry : args->seq->ry;
		siril_debug_print("size: %d %d\n", orig_rx, orig_ry);
		cvApplyFlips(&Hs, orig_ry, 0);
		reframe_wcs(result->keywords.wcslib, &Hs);
		update_wcsdata_from_wcs(result);
	}
	free(ssdata);
	return args->retval;
}

