/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "registration/registration.h"
#include "algos/siril_wcs.h"
#include "opencv/opencv.h"

#include "stacking/stacking.h"

/******************************* ADDMIN AND ADDMAX STACKING ******************************
 * These methods are very close to summing stacking except that the result
 * takes only the pixel if it is brighter (max) or dimmer (min) than the
 * previous one at the same coordinates.
 *
 * Images are processed sequentially (max_parallel_images=1) so that the
 * per-pixel comparison-update needs no locking. Parallelism is achieved by
 * splitting the pixel loop across threads within each image hook.
 */

struct minmax_stacking_data {
	WORD  *final_pixel[3];		// integer accumulation buffers
	float *ffinal_pixel[3];		// float accumulation buffers
	GList *list_date;
	double livetime;
	int reglayer;
	int ref_image;
	gboolean input_32bits;
	gboolean ismax;
	gboolean maximize_framing;
	gboolean upscale_at_stacking;
	int output_size[2];
	int offset[2];
	fits result;
};

static int minmax_stacking_prepare_hook(struct generic_seq_args *args) {
	struct minmax_stacking_data *mmdata = args->user;
	size_t nbdata = (size_t)mmdata->output_size[0] * mmdata->output_size[1];
	size_t nbpixels = nbdata * args->seq->nb_layers;

	if (mmdata->input_32bits) {
		mmdata->ffinal_pixel[0] = malloc(nbpixels * sizeof(float));
		if (!mmdata->ffinal_pixel[0]) {
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		float init = mmdata->ismax ? 0.0f : 1.0f;
		for (size_t k = 0; k < nbpixels; k++)
			mmdata->ffinal_pixel[0][k] = init;
		if (args->seq->nb_layers == 3) {
			mmdata->ffinal_pixel[1] = mmdata->ffinal_pixel[0] + nbdata;
			mmdata->ffinal_pixel[2] = mmdata->ffinal_pixel[0] + nbdata * 2;
		}
		mmdata->final_pixel[0] = NULL;
	} else {
		mmdata->final_pixel[0] = malloc(nbpixels * sizeof(WORD));
		if (!mmdata->final_pixel[0]) {
			PRINT_ALLOC_ERR;
			return ST_ALLOC_ERROR;
		}
		WORD init = mmdata->ismax ? 0 : USHRT_MAX;
		for (size_t k = 0; k < nbpixels; k++)
			mmdata->final_pixel[0][k] = init;
		if (args->seq->nb_layers == 3) {
			mmdata->final_pixel[1] = mmdata->final_pixel[0] + nbdata;
			mmdata->final_pixel[2] = mmdata->final_pixel[0] + nbdata * 2;
		}
		mmdata->ffinal_pixel[0] = NULL;
	}

	mmdata->livetime = 0.0;
	mmdata->list_date = NULL;
	return ST_OK;
}

/* Images are processed one at a time (max_parallel_images=1), so livetime and
 * list_date updates here are safe without locks. The pixel loop is parallelized
 * across threads using the 'threads' subthread count. */
static int minmax_stacking_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct minmax_stacking_data *mmdata = args->user;
	int shiftx = 0, shifty = 0;
	size_t nbdata = (size_t)mmdata->output_size[0] * mmdata->output_size[1];
	gboolean ismax = mmdata->ismax;

	mmdata->livetime += fit->keywords.exposure;

	if (fit->keywords.date_obs) {
		GDateTime *date = g_date_time_ref(fit->keywords.date_obs);
		mmdata->list_date = g_list_prepend(mmdata->list_date, new_date_item(date, fit->keywords.exposure));
	}

	if (mmdata->reglayer != -1 && args->seq->regparam[mmdata->reglayer]) {
		double scale = mmdata->upscale_at_stacking ? 2. : 1.;
		double dx, dy;
		translation_from_H(args->seq->regparam[mmdata->reglayer][i].H, &dx, &dy);
		dx *= scale;
		dy *= scale;
		dx -= mmdata->offset[0];
		dy -= mmdata->offset[1];
		if (mmdata->maximize_framing)
			dy -= (double)fit->ry;
		shiftx = round_to_int(dx);
		shifty = round_to_int(dy);
		siril_debug_print("img %d dx %d dy %d\n", o, shiftx, shifty);
	}
	if (shiftx == INT_MIN) {
		siril_debug_print("Error: image #%d has a wrong shiftx value\n", o + 1);
		shiftx += 1;
	}
	if (shifty == INT_MIN) {
		siril_debug_print("Error: image #%d has a wrong shifty value\n", o + 1);
		shifty += 1;
	}

	/* Each y row is independent: different y -> different pixel range in output
	 * buffer -> no write conflicts between threads. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for (int y = 0; y < mmdata->output_size[1]; ++y) {
		int ny = y - shifty;
		if (ny < 0 || ny >= fit->ry)
			continue;
		size_t rny = (size_t)ny * fit->rx;
		size_t pixel = (size_t)y * mmdata->output_size[0];
		for (int x = 0; x < mmdata->output_size[0]; ++x, ++pixel) {
			int nx = x - shiftx;
			if (nx >= 0 && nx < fit->rx) {
				size_t ii = rny + nx;
				if (ii < nbdata) {
					for (int layer = 0; layer < args->seq->nb_layers; ++layer) {
						if (mmdata->input_32bits) {
							float cur = fit->fpdata[layer][ii];
							if ((ismax && cur > mmdata->ffinal_pixel[layer][pixel]) ||
									(!ismax && cur < mmdata->ffinal_pixel[layer][pixel]))
								mmdata->ffinal_pixel[layer][pixel] = cur;
						} else {
							WORD cur = fit->pdata[layer][ii];
							if ((ismax && cur > mmdata->final_pixel[layer][pixel]) ||
									(!ismax && cur < mmdata->final_pixel[layer][pixel]))
								mmdata->final_pixel[layer][pixel] = cur;
						}
					}
				}
			}
		}
	}
	return ST_OK;
}

static int minmax_stacking_finalize_hook(struct generic_seq_args *args) {
	struct minmax_stacking_data *mmdata = args->user;

	if (args->retval) {
		if (mmdata->final_pixel[0])  free(mmdata->final_pixel[0]);
		if (mmdata->ffinal_pixel[0]) free(mmdata->ffinal_pixel[0]);
		args->user = NULL;
		return args->retval;
	}

	fits *result = &mmdata->result;
	if (mmdata->input_32bits) {
		if (new_fit_image_with_data(&result, mmdata->output_size[0], mmdata->output_size[1],
				args->seq->nb_layers, DATA_FLOAT, mmdata->ffinal_pixel[0]))
			return ST_GENERIC_ERROR;
	} else {
		if (new_fit_image_with_data(&result, mmdata->output_size[0], mmdata->output_size[1],
				args->seq->nb_layers, DATA_USHORT, mmdata->final_pixel[0]))
			return ST_GENERIC_ERROR;
	}
	/* buffers are now owned by result, do not free them */

	int ref = mmdata->ref_image;
	if (args->seq->type == SEQ_REGULAR) {
		if (!seq_open_image(args->seq, ref)) {
			import_metadata_from_fitsfile(args->seq->fptr[ref], result);
			result->orig_bitpix = result->bitpix = args->seq->bitpix;
			seq_close_image(args->seq, ref);
		}
		result->keywords.livetime = mmdata->livetime;
	} else if (args->seq->type == SEQ_FITSEQ) {
		if (!fitseq_set_current_frame(args->seq->fitseq_file, ref)) {
			import_metadata_from_fitsfile(args->seq->fitseq_file->fptr, result);
			result->orig_bitpix = result->bitpix = args->seq->fitseq_file->bitpix;
		}
		result->keywords.livetime = mmdata->livetime;
	} else if (args->seq->type == SEQ_SER) {
		import_metadata_from_serfile(args->seq->ser_file, result);
		result->orig_bitpix = result->bitpix = (args->seq->ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) ? BYTE_IMG : USHORT_IMG;
		result->keywords.livetime = result->keywords.exposure * args->nb_filtered_images;
	}

	result->keywords.stackcnt = args->nb_filtered_images;
	compute_date_time_keywords(mmdata->list_date, result);
	g_list_free_full(mmdata->list_date, (GDestroyNotify) free_list_date);

	args->user = NULL;
	return ST_OK;
}

static int stack_minmax_generic(struct stacking_args *stackargs, gboolean ismax) {
	struct generic_seq_args *args = create_default_seqargs(stackargs->seq);
	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = minmax_stacking_prepare_hook;
	args->image_hook = minmax_stacking_image_hook;
	args->finalize_hook = minmax_stacking_finalize_hook;
	args->description = ismax ? _("Max stacking") : _("Min stacking");
	args->already_in_a_thread = TRUE;
	args->max_parallel_images = 1; // sequential images; pixel loop is parallelized inside image_hook

	struct minmax_stacking_data *mmdata = calloc(1, sizeof(struct minmax_stacking_data));
	mmdata->reglayer = stackargs->reglayer;
	mmdata->ref_image = stackargs->ref_image;
	mmdata->input_32bits = get_data_type(args->seq->bitpix) == DATA_FLOAT;
	mmdata->ismax = ismax;
	mmdata->upscale_at_stacking = stackargs->upscale_at_stacking;

	gboolean update_wcs = FALSE;
	if (stackargs->maximize_framing) {
		compute_max_framing(stackargs, mmdata->output_size, mmdata->offset);
		mmdata->maximize_framing = TRUE;
		update_wcs = TRUE;
	} else {
		mmdata->output_size[0] = args->seq->rx;
		mmdata->output_size[1] = args->seq->ry;
		double dx = 0., dy = 0.;
		if (mmdata->reglayer >= 0) {
			translation_from_H(args->seq->regparam[mmdata->reglayer][mmdata->ref_image].H, &dx, &dy);
			update_wcs = TRUE;
		}
		mmdata->offset[0] = (int)dx;
		mmdata->offset[1] = (int)dy;
	}

	args->user = mmdata;
	generic_sequence_worker(args);
	memcpy(&stackargs->result, &mmdata->result, sizeof(fits));

	if (update_wcs && has_wcs(&mmdata->result)) {
		fits *result = &mmdata->result;
		Homography Hs = { 0 };
		cvGetEye(&Hs);
		double dx, dy;
		translation_from_H(args->seq->regparam[stackargs->reglayer][stackargs->ref_image].H, &dx, &dy);
		siril_debug_print("ref shift: %d %d\n", (int)dx, (int)dy);
		siril_debug_print("crpix: %.1f %.1f\n", result->keywords.wcslib->crpix[0], result->keywords.wcslib->crpix[1]);
		Hs.h02  = dx - mmdata->offset[0];
		Hs.h12 -= dy - mmdata->offset[1];
		int orig_rx = (args->seq->is_variable) ? args->seq->imgparam[args->seq->reference_image].rx : args->seq->rx;
		int orig_ry = (args->seq->is_variable) ? args->seq->imgparam[args->seq->reference_image].ry : args->seq->ry;
		siril_debug_print("size: %d %d\n", orig_rx, orig_ry);
		cvApplyFlips(&Hs, orig_ry, 0);
		reframe_wcs(result->keywords.wcslib, &Hs);
		update_wcsdata_from_wcs(result);
	}

	free(mmdata);
	return args->retval;
}

int stack_addmax(struct stacking_args *args) {
	return stack_minmax_generic(args, TRUE);
}

int stack_addmin(struct stacking_args *args) {
	return stack_minmax_generic(args, FALSE);
}
