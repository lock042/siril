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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "drizzle/cdrizzlebox.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzleutil.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"
#include "algos/demosaicing.h"

#include "opencv/opencv.h"

#define MIN_RATIO 0.1 // minimum fraction of ref image dimensions to validate output sequence is worth creating

static void create_output_sequence_for_apply_reg(struct registration_args *args);
static int new_ref_index = -1;

static Homography Htransf = {0};
static int rx_out, ry_out, ry_out_unscaled;

regdata *apply_reg_get_current_regdata(struct registration_args *regargs) {
	regdata *current_regdata;
	if (regargs->seq->regparam[regargs->layer]) {
		siril_log_message(
				_("Applying existing registration from layer #%d to transform the images\n"), regargs->layer);
		current_regdata = regargs->seq->regparam[regargs->layer];
	} else {
		siril_log_message(
				_("No registration data exists for this layer\n"));
		return NULL;
	}
	return current_regdata;
}

static gboolean compute_framing(struct registration_args *regargs) {
	// validity of matrices has already been checked before this call
	// and null matrices have been discarded
	Homography Href = regargs->seq->regparam[regargs->layer][regargs->reference_image].H;
	Homography Hshift = {0};
	cvGetEye(&Hshift);
	int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	int x0, y0, rx_0 = rx, ry_0 = ry, n;
	double xmin, xmax, ymin, ymax, cogx, cogy;

	regframe framing = { 0 };
	framing.pt[0].x = 0.;
	framing.pt[0].y = 0.;
	framing.pt[1].x = (double)rx;
	framing.pt[1].y = 0.;
	framing.pt[2].x = (double)rx;
	framing.pt[2].y = (double)ry;
	framing.pt[3].x = 0.;
	framing.pt[3].y = (double)ry;

	switch (regargs->framing) {
		case FRAMING_CURRENT:
			break;
		case FRAMING_MAX:
			xmin = DBL_MAX;
			xmax = -DBL_MAX;
			ymin = DBL_MAX;
			ymax = -DBL_MAX;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				if (regargs->seq->is_variable) {
					double rx2 = (double)regargs->seq->imgparam[i].rx;
					double ry2 = (double)regargs->seq->imgparam[i].ry;
					current_framing.pt[1].x = rx2;
					current_framing.pt[2].x = rx2;
					current_framing.pt[2].y = ry2;
					current_framing.pt[3].y = ry2;
				}
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y, regargs->seq->regparam[regargs->layer][i].H, Href);
					if (xmin > current_framing.pt[j].x) xmin = current_framing.pt[j].x;
					if (ymin > current_framing.pt[j].y) ymin = current_framing.pt[j].y;
					if (xmax < current_framing.pt[j].x) xmax = current_framing.pt[j].x;
					if (ymax < current_framing.pt[j].y) ymax = current_framing.pt[j].y;
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
				}
			}
			rx_0 = (int)(ceil(xmax) - floor(xmin));
			ry_0 = (int)(ceil(ymax) - floor(ymin));
			x0 = floor(xmin);
			y0 = floor(ymin);
			siril_debug_print("new size: %d %d\n", rx_0, ry_0);
			siril_debug_print("new origin: %d %d\n", x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		case FRAMING_MIN:
			xmin = -DBL_MAX;
			xmax = DBL_MAX;
			ymin = -DBL_MAX;
			ymax = DBL_MAX;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				if (regargs->seq->is_variable) {
					double rx2 = (double)regargs->seq->imgparam[i].rx;
					double ry2 = (double)regargs->seq->imgparam[i].ry;
					current_framing.pt[1].x = rx2;
					current_framing.pt[2].x = rx2;
					current_framing.pt[2].y = ry2;
					current_framing.pt[3].y = ry2;
				}
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				if (xmin < xs[1]) xmin = xs[1];
				if (ymin < ys[1]) ymin = ys[1];
				if (xmax > xs[2]) xmax = xs[2];
				if (ymax > ys[2]) ymax = ys[2];
			}
			rx_0 = (int)(floor(xmax) - ceil(xmin));
			ry_0 = (int)(floor(ymax) - ceil(ymin));
			x0 = ceil(xmin);
			y0 = ceil(ymin);
			siril_debug_print("new size: %d %d\n", rx_0, ry_0);
			siril_debug_print("new origin: %d %d\n", x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		case FRAMING_COG:
			cogx = 0.;
			cogy = 0;
			n = 0;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				double currcogx = 0., currcogy = 0.;
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					currcogx += current_framing.pt[j].x;
					currcogy += current_framing.pt[j].y;
				}
				cogx += currcogx * 0.25;
				cogy += currcogy * 0.25;
				n++;
			}
			cogx /= (double)n;
			cogy /= (double)n;
			x0 = (int)(cogx - (double)rx * 0.5);
			y0 = (int)(cogy - (double)ry * 0.5);
			siril_log_message(_("Framing: Shift from reference origin: %d, %d\n"), x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		default:
			return FALSE;
	}
	cvMultH(Href, Hshift, &Htransf);
	ry_out_unscaled = ry_0;
	if (regargs->driz) {
		rx_out = rx_0 * regargs->driz->scale;
		ry_out = ry_0 * regargs->driz->scale;
	} else {
		rx_out = rx_0 * ((regargs->x2upscale) ? 2. : 1.);
		ry_out = ry_0 * ((regargs->x2upscale) ? 2. : 1.);
	}
	return TRUE;
}

int apply_reg_prepare_results(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	// allocate destination sequence data
	regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
	regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
	if (!regargs->imgparam  || !regargs->regparam) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	if (seq_prepare_hook(args))
		return 1;

	sadata->success = calloc(args->nb_filtered_images, sizeof(BYTE));
	if (!sadata->success) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}

struct _drizzle_pair {
	int index;
	fits *out;
	fits *output_counts;
};

static int apply_drz_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;

	/* preparing reference data from reference fit and making sanity checks*/
	sadata->current_regdata = apply_reg_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;

	int number_of_outputs = 1;
	// we call the generic prepare twice with different prefixes
	args->new_seq_prefix = regargs->prefix;
	if (apply_reg_prepare_results(args))
		return 1;
	// but we copy the result between each call
	driz->new_ser_drz = args->new_ser;
	driz->new_fitseq_drz = args->new_fitseq;

	if (driz->keep_counts) {
		args->new_seq_prefix = "oc_"; // This is OK here, it does not get freed in the
		// end_generic_sequence later
		if (apply_reg_prepare_results(args))
			return 1;
		driz->new_ser_pxcnt = args->new_ser;
		driz->new_fitseq_pxcnt = args->new_fitseq;

		args->new_seq_prefix = regargs->prefix; // Put it back so it gets loaded on completion
		args->new_ser = NULL;
		args->new_fitseq = NULL;
		number_of_outputs = 2;
	}

	seqwriter_set_number_of_outputs(number_of_outputs);

	return 0;
}

struct _double_driz {
	int index;
	fits *output;
	fits *pixel_count;
};

int apply_reg_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	fits fit = { 0 };

	/* preparing reference data from reference fit and making sanity checks*/
	sadata->current_regdata = apply_reg_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;


	if (seq_read_frame_metadata(args->seq, regargs->reference_image, &fit)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		return 1;
	}
	if (fit.naxes[2] == 1 && fit.keywords.bayer_pattern[0] != '\0')
		siril_log_color_message(_("Applying transformation on a sequence opened as CFA is a bad idea.\n"), "red");
	free_wcs(&fit);
	reset_wcsdata(&fit);
	return apply_reg_prepare_results(args);
}

/* reads the image and apply existing transformation */
int apply_reg_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;
	float scale = 2.f;

	Homography H = { 0 };
	Homography Himg = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
		siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
	}

	// Composing transformation wrt reference image
	Himg = regargs->seq->regparam[regargs->layer][in_index].H;
	if (guess_transform_from_H(Himg) == NULL_TRANSFORMATION)
		return 1; // in case H is null and -selected was not passed
	cvTransfH(Himg, Htransf, &H);

	struct driz_param_t *p = NULL;
	if (regargs->driz) {
		scale = driz->scale;
		p = calloc(1, sizeof(struct driz_param_t));
		driz_param_init(p);
		p->kernel = driz->kernel;
		p->driz = driz;
		p->error = malloc(sizeof(struct driz_error_t));
		p->scale = driz->scale;
		p->pixel_fraction = driz->pixel_fraction;
		p->cfa = driz->cfa;
		// Set bounds equal to whole image
		p->xmin = p->ymin = 0;
		p->xmax = fit->rx - 1;
		p->ymax = fit->ry - 1;
		p->pixmap = calloc(1, sizeof(imgmap_t));
		p->pixmap->rx = fit->rx;
		p->pixmap->ry = fit->ry;

		map_image_coordinates_h(fit, H, p->pixmap, ry_out_unscaled, driz->scale);
		if (!p->pixmap->xmap) {
			siril_log_color_message(_("Error generating mapping array.\n"), "red");
			free(p->error);
			free(p->pixmap);
			free(p);
			return 1;
		}
		p->data = fit;
		// Convert fit to 32-bit float if required
		float *newbuf = NULL;
		if (fit->type == DATA_USHORT) {
			siril_debug_print("Replacing ushort buffer for drizzling\n");
			size_t ndata = fit->rx * fit->ry * fit->naxes[2];
			newbuf = malloc(ndata * sizeof(float));
			float invnorm = 1.f / USHRT_MAX_SINGLE;
			for (size_t i = 0 ; i < ndata ; i++) {
				newbuf[i] = fit->data[i] * invnorm;
			}
			fit_replace_buffer(fit, newbuf, DATA_FLOAT);
		}
		/* Set up output fits */
		fits out;
		copyfits(fit, &out, CP_FORMAT, -1);
		out.rx = out.naxes[0] = rx_out;
		out.ry = out.naxes[1] = ry_out;
		out.naxes[2] = driz->is_bayer ? 3 : 1;
		siril_debug_print("Output image %d x %d x %ld\n", out.rx, out.ry, out.naxes[2]);
		size_t chansize = out.rx * out.ry * sizeof(float);
		out.fdata = calloc(out.naxes[2] * chansize, 1);
		out.fpdata[0] = out.fdata;
		out.fpdata[1] = out.naxes[2] == 1 ? out.fdata : out.fdata + chansize;
		out.fpdata[2] = out.naxes[2] == 1 ? out.fdata : out.fdata + 2 * chansize;
		p->output_data = &out;

		// Set up the output_counts fits to store pixel hit counts
		fits *output_counts = calloc(1, sizeof(fits));
		copyfits(&out, output_counts, CP_FORMAT, -1);
		siril_debug_print("Output counts image %d x %d x %ld\n", out.rx, out.ry, out.naxes[2]);
		output_counts->fdata = calloc(output_counts->rx * output_counts->ry * output_counts->naxes[2], sizeof(float));
		p->output_counts = output_counts;

		p->weights = driz->flat;

		if (dobox(p)) { // Do the drizzle
			siril_log_color_message("s\n", p->error->last_message);
			return 1;
		}
		clearfits(fit);
		copyfits(&out, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		clearfits(&out);

		if (driz->is_bayer) {
			/* we need to do something special here because it's a 1-channel sequence and
			* image that will become 3-channel after this call, so stats caching will
			* not be correct. Preprocessing does not compute new stats so we don't need
			* to save them into cache now.
			* Destroying the stats from the fit will also prevent saving them in the
			* sequence, which may be done automatically by the caller.
			*/
			full_stats_invalidation_from_fit(fit);
			fit->keywords.lo = 0;
		}

		free(p->pixmap->xmap);
		free(p->pixmap);

		// Get rid of the output_counts if not required
		if (!driz->keep_counts) {
			clearfits(output_counts);
			free(output_counts);
			output_counts = NULL;
		}

		struct _double_driz *double_data = calloc(1, sizeof(struct _double_driz));
		double_data->index = out_index;
		double_data->output = fit;
		double_data->pixel_count = output_counts;
#ifdef _OPENMP
		omp_set_lock(&args->lock);
#endif
		driz->processed_images = g_list_append(driz->processed_images, double_data);
#ifdef _OPENMP
		omp_unset_lock(&args->lock);
#endif
		siril_debug_print("drizzle: processed image and output_counts (if specified) added to the save list (index %d)\n", out_index);
	}
	else {
		if (regargs->interpolation <= OPENCV_LANCZOS4) {
			if (cvTransformImage(fit, rx_out, ry_out, H, regargs->x2upscale, regargs->interpolation, regargs->clamp)) {
				return 1;
			}
		} else {
			if (shift_fit_from_reg(fit, H)) {
				return 1;
			}
		}
	}
	if (in_index == regargs->reference_image)
		new_ref_index = out_index; // keeping track of the new ref index in output sequence

	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	regargs->imgparam[out_index].rx = rx_out;
	regargs->imgparam[out_index].ry = ry_out;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
	regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
	regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
	cvGetEye(&regargs->regparam[out_index].H);

	if (regargs->x2upscale) {
		fit->keywords.pixel_size_x /= scale;
		fit->keywords.pixel_size_y /= scale;
		regargs->regparam[out_index].fwhm *= scale;
		regargs->regparam[out_index].weighted_fwhm *= scale;
		regargs->imgparam[out_index].rx *= scale;
		regargs->imgparam[out_index].ry *= scale;
	}

	sadata->success[out_index] = 1;
	return 0;
}

static int apply_drz_save_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;
	struct _double_driz *double_data = NULL;
	// images are passed from the image_hook to the save in a list, because
	// there are two, which is unsupported by the generic arguments
#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	GList *list = driz->processed_images;
	while (list) {
		if (((struct _double_driz *)list->data)->index == out_index) {
			double_data = list->data;
			break;
		}
		list = g_list_next(driz->processed_images);
	}
	if (double_data)
		driz->processed_images = g_list_remove(driz->processed_images, double_data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	if (!double_data) {
		siril_log_color_message(_("Image %d not found for writing\n"), "red", in_index);
		return 1;
	}

	siril_debug_print("Drizzle: images to be saved (index %d)\n", out_index);
	if (double_data->output->naxes[0] == 0) {
		siril_debug_print("empty drizzle data\n");
		return 1;
	}

	int retval1 = 0, retval2 = 0;
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		if (double_data->output) {
			if (double_data->output->type != DATA_USHORT) {
				fit_replace_buffer(	double_data->output,
									float_buffer_to_ushort( double_data->output->fdata,
															(double_data->output->rx *
															 double_data->output->ry *
															 double_data->output->naxes[2])),
									DATA_USHORT);
			}
			retval1 = ser_write_frame_from_fit(driz->new_ser_drz, double_data->output, out_index);
		}
		if (double_data->pixel_count) {
			float factor = max(1.f, (1.f / (driz->scale * driz->scale)));
			soper(double_data->pixel_count, factor, OPER_MUL, FALSE);
			if (double_data->pixel_count->type != DATA_USHORT) {
				fit_replace_buffer(	double_data->pixel_count,
									float_buffer_to_ushort( double_data->pixel_count->fdata,
															(double_data->pixel_count->rx *
															 double_data->pixel_count->ry *
															 double_data->pixel_count->naxes[2])),
									DATA_USHORT);
			}
			retval2 = ser_write_frame_from_fit(driz->new_ser_pxcnt, double_data->pixel_count, out_index);
		}
		// the two fits are freed by the writing thread
		if (retval1 || retval2) {
			// special case because it's not done in the generic
			clearfits(fit);
			free(fit);
		}
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		if (double_data->output) {
			if (com.pref.force_16bit && double_data->output->type != DATA_USHORT) {
				fit_replace_buffer(	double_data->output,
									float_buffer_to_ushort( double_data->output->fdata,
															(double_data->output->rx *
															 double_data->output->ry *
															 double_data->output->naxes[2])),
									DATA_USHORT);
			}
			retval1 = fitseq_write_image(driz->new_fitseq_drz, double_data->output, out_index);
		}
		if (double_data->pixel_count) {
			if (double_data->output->type == DATA_USHORT) {
				float factor = max(1.f, (1.f / (driz->scale * driz->scale)));
				soper(double_data->pixel_count, factor, OPER_MUL, FALSE);
				fit_replace_buffer(	double_data->pixel_count,
									float_buffer_to_ushort( double_data->pixel_count->fdata,
															(double_data->pixel_count->rx *
															 double_data->pixel_count->ry *
															 double_data->pixel_count->naxes[2])),
									DATA_USHORT);
			}
			retval2 = fitseq_write_image(driz->new_fitseq_pxcnt, double_data->pixel_count, out_index);
		}
		// the two fits are freed by the writing thread
		if (retval1 || retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, regargs->prefix, in_index);
		if (double_data->output) {
			if (com.pref.force_16bit) {
				fit_replace_buffer(	double_data->output,
									float_buffer_to_ushort( double_data->output->fdata,
															(double_data->output->rx *
															 double_data->output->ry *
															 double_data->output->naxes[2])),
									DATA_USHORT);
			}
			retval1 = savefits(dest, double_data->output);
		}
		free(dest);
		if (double_data->pixel_count) {
			driz->pixcnt_prefix = g_strdup_printf("oc_%s", regargs->prefix);
			dest = fit_sequence_get_image_filename_prefixed(args->seq, driz->pixcnt_prefix, in_index);
			if (com.pref.force_16bit) {
				float factor = max(1.f, (1.f / (driz->scale * driz->scale)));
				if (fabs(factor - 1.f) > FLT_EPSILON)
					soper(double_data->pixel_count, factor, OPER_MUL, FALSE);
				fit_replace_buffer(	double_data->pixel_count,
									float_buffer_to_ushort( double_data->pixel_count->fdata,
															(double_data->pixel_count->rx *
															 double_data->pixel_count->ry *
															 double_data->pixel_count->naxes[2])),
									DATA_USHORT);
			}
			retval2 = savefits(dest, double_data->pixel_count);
			free(dest);
		}
	}
	if (!driz->new_fitseq_drz && !driz->new_ser_drz) { // detect if there is a seqwriter
		clearfits(double_data->output);
		// Don't free double_data->output as this == fit and needs to be freed by the sequence worker
		clearfits(double_data->pixel_count);
		free(double_data->pixel_count);
	}
	free(double_data);
	return retval1 || retval2;
}

int apply_drz_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;
	int retval = 0;

	if (driz->keep_counts) {
		// We deal with the pixel_count sequence first
		args->new_ser = driz->new_ser_pxcnt;
		args->new_fitseq = driz->new_fitseq_pxcnt;
		if (!args->retval) {
			retval = seq_finalize_hook(args);
		} else {
			// same as seq_finalize_hook but with file deletion
			if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
				ser_close_and_delete_file(args->new_ser);
				free(args->new_ser);
			} else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
				fitseq_close_and_delete_file(args->new_fitseq);
				free(args->new_fitseq);
			} else if (args->seq->type == SEQ_REGULAR) {
				remove_prefixed_sequence_files(args->seq, driz->pixcnt_prefix);
			}
		}
		driz->new_ser_pxcnt = NULL;
		driz->new_fitseq_pxcnt = NULL;
	}

	// We deal with the drizzled sequence next
	args->new_ser = driz->new_ser_drz;
	args->new_fitseq = driz->new_fitseq_drz;

	int failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;
		if (failed) {
			// driz->imgparam and driz->regparam may have holes caused by images
			// that failed to be registered - compact them
			for (int i = 0, j = 0; i < regargs->new_total; i++, j++) {
				while (!sadata->success[j] && j < args->nb_filtered_images) j++;
				g_assert(sadata->success[j]);
				if (i != j) {
					regargs->imgparam[i] = regargs->imgparam[j];
					regargs->regparam[i] = regargs->regparam[j];
				}
			}
		}
		seq_finalize_hook(args);
	} else {
		regargs->new_total = 0;

		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		} else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		} else if (args->seq->type == SEQ_REGULAR) {
			remove_prefixed_sequence_files(args->seq, regargs->prefix);
		}
	}

	if (sadata->success) free(sadata->success);
	//	Do not free driz here as it will be needed for blot and second drizzle
	args->user = NULL;

	if (!args->retval) {
		siril_log_message(_("Applying drizzle completed.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d exported.\n"), "green", failed, regargs->new_total);
		regargs->load_new_sequence = TRUE;
		g_free(str);
		if (!(args->seq->type == SEQ_INTERNAL)) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_apply_reg(regargs);
			// will be loaded in the idle function if (load_new_sequence)
			args->load_new_sequence = regargs->load_new_sequence; // always want to load drizzled sequence
		}
	}
	else {
		siril_log_message(_("Drizzle aborted.\n"));
	}

	if (regargs->new_total == 0)
		return 1;
	driz->new_ser_drz = NULL;
	driz->new_fitseq_drz = NULL;

	return retval;
}

int apply_reg_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;
		if (failed) {
			// regargs->imgparam and regargs->regparam may have holes caused by images
			// that failed to be registered - compact them
			for (int i = 0, j = 0; i < regargs->new_total; i++, j++) {
				while (!sadata->success[j] && j < args->nb_filtered_images) j++;
				g_assert(sadata->success[j]);
				if (i != j) {
					regargs->imgparam[i] = regargs->imgparam[j];
					regargs->regparam[i] = regargs->regparam[j];
				}
			}
		}
		seq_finalize_hook(args);
	} else {
		regargs->new_total = 0;

		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		} else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		} else if (args->seq->type == SEQ_REGULAR) {
			remove_prefixed_sequence_files(regargs->seq, regargs->prefix);
		}
	}

	if (sadata->success) free(sadata->success);
	free(sadata);
	args->user = NULL;

	if (!args->retval) {
		siril_log_message(_("Applying registration completed.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d exported.\n"), "green", failed, regargs->new_total);

		g_free(str);
		if (!(args->seq->type == SEQ_INTERNAL)) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_apply_reg(regargs);
			// will be loaded in the idle function if (load_new_sequence)
			regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Transformation aborted.\n"));
	}
	return regargs->new_total == 0;
}

int apply_drz_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, driz->scale, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
	int float_multiplier = (is_float) ? 1 : 2;
	int is_bayer = driz->is_bayer;
	MB_per_scaled_image *= is_float ? 1 : 2; // Output is always float
	if (is_bayer) {
		MB_per_scaled_image *= 3;
	}
	unsigned int MB_per_float_image = MB_per_orig_image * float_multiplier;
	unsigned int MB_per_float_channel = MB_per_float_image;

	/* The drizzle memory consumption is:
		* the original image
		* Two 2 * rx * ry * float arrays for computing mapping
		* (note this could become double arrays with WCS?)
		  (the mapping file reallocs one of these so never exceeds that amount)
		* the weights file (1 * scaled image float data)
		* the transformed image, including scaling factor if required
		* the output counts image, the same size as the scaled image
	*/
	unsigned int required = MB_per_orig_image + 4 * MB_per_float_channel + 2 * MB_per_scaled_image;

	if (limit > 0) {

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int apply_reg_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;

	/* The transformation memory consumption is:
		* the original image
		* the transformed image, including upscale if required (4x)
		*/
	unsigned int required = MB_per_orig_image + MB_per_scaled_image;
	// If interpolation clamping is set, 2x additional Mats of the same format
	// as the original image are required
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	if (regargs->clamp && (regargs->interpolation == OPENCV_CUBIC ||
			regargs->interpolation == OPENCV_LANCZOS4)) {
		float factor = (is_float) ? 0.25 : 0.5;
		required += (1 + factor) * MB_per_scaled_image;
	}
	regargs = NULL;
	sadata = NULL;

	if (limit > 0) {

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int register_apply_reg(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	struct driz_args_t *driz = regargs->driz;
	control_window_switch_to_tab(OUTPUT_LOGS);

	if (regargs->driz) {
		set_progress_bar_data(_("Initializing drizzle data..."), PROGRESS_PULSATE);
	}

	if (!check_before_applyreg(regargs)) {
		free(args);
		return -1;
	}

	if (driz) {
		args->prepare_hook = apply_drz_prepare_hook;
		args->save_hook = apply_drz_save_hook;
		args->finalize_hook = apply_drz_finalize_hook;
		args->upscale_ratio = 1.0; // Drizzle scale dealt with separately
		driz_param_dump(driz); // Print some info to the log
		/* preparing reference data from reference fit and making sanity checks*/
		fits fit = { 0 };

		/* fit will now hold the reference frame */
		if (seq_read_frame_metadata(regargs->seq, regargs->reference_image, &fit)) {
			siril_log_message(_("Could not load reference image\n"));
			args->seq->regparam[0] = NULL;
			free(args);
			return 1;
		}
		sensor_pattern pattern;
		if (args->seq->type == SEQ_SER && args->seq->ser_file ) {
			driz->is_bayer = TRUE;
			switch (args->seq->ser_file->color_id) {
				case SER_MONO:
					pattern = get_cfa_pattern_index_from_string("");
					driz->is_bayer = FALSE;
					break;
				case SER_BAYER_RGGB:
					pattern = get_cfa_pattern_index_from_string("RGGB");
					break;
				case SER_BAYER_GRBG:
					pattern = get_cfa_pattern_index_from_string("GRBG");
					break;
				case SER_BAYER_GBRG:
					pattern = get_cfa_pattern_index_from_string("GBRG");
					break;
				case SER_BAYER_BGGR:
					pattern = get_cfa_pattern_index_from_string("BGGR");
					break;
				default:
					siril_log_message(_("Unsupported SER CFA pattern detected. Treating as mono.\n"));
					driz->is_bayer = FALSE;
			}
		} else {
			driz->is_bayer = (fit.keywords.bayer_pattern[0] != '\0'); // If there is a CFA pattern we need to CFA drizzle
			pattern = com.pref.debayer.use_bayer_header ? get_cfa_pattern_index_from_string(fit.keywords.bayer_pattern) : com.pref.debayer.bayer_pattern;
		}
		if (driz->is_bayer) {
			adjust_Bayer_pattern(&fit, &pattern);
			driz->cfa = get_cfa_from_pattern(pattern);
			if (!driz->cfa) // if fit.bayer_pattern exists and get_cfa_from_pattern returns NULL then there is a problem.
				return 1;
		}
	} else {
		args->prepare_hook = apply_reg_prepare_hook;
		args->finalize_hook = apply_reg_finalize_hook;
		args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	}
	args->filtering_criterion = regargs->filtering_criterion;
	args->filtering_parameter = regargs->filtering_parameter;
	args->nb_filtered_images = regargs->new_total;
	args->compute_mem_limits_hook = apply_reg_compute_mem_limits;
	args->image_hook = apply_reg_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Apply registration");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = regargs->prefix;
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free(args);
		return -1;
	}
	sadata->regargs = regargs;
	args->user = sadata;

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}

static void create_output_sequence_for_apply_reg(struct registration_args *args) {
	sequence seq = { 0 };
	initialize_sequence(&seq, TRUE);

	/* we are not interested in the whole path */
	gchar *seqname = g_path_get_basename(args->seq->seqname);
	char *rseqname = malloc(
			strlen(args->prefix) + strlen(seqname) + 5);
	sprintf(rseqname, "%s%s.seq", args->prefix, seqname);
	g_free(seqname);
	g_unlink(rseqname);	// remove previous to overwrite
	args->new_seq_name = remove_ext_from_filename(rseqname);
	free(rseqname);
	seq.seqname = strdup(args->new_seq_name);
	if (args->driz) {
		gchar *tmpbuf = g_strdup_printf("oc_%s", seq.seqname);
		siril_debug_print("seq name: %s\npix_cnt seq name: %s\n", args->seq->seqname, tmpbuf);
		seq.pixcnt_seqname = strdup(tmpbuf);
		g_free(tmpbuf);
	}
	seq.number = args->new_total;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = (args->driz && args->driz->is_bayer) ? 3 : args->seq->nb_layers;
	seq.rx = rx_out;
	seq.ry = ry_out;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number-1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	seq.is_variable = FALSE;
	seq.fz = com.pref.comp.fits_enabled;
	// update with the new numbering
	seq.reference_image = new_ref_index;
	seq.needs_saving = TRUE;
	writeseqfile(&seq);
	free_sequence(&seq, FALSE);
	new_ref_index = -1; // resetting
}

transformation_type guess_transform_from_H(Homography H) {
	if (fabs(H.h00 + H.h01 + H.h02 + H.h10 + H.h11 + H.h12 + H.h20 + H.h21 + H.h22) < __DBL_EPSILON__)
		return NULL_TRANSFORMATION;
	if (fabs(H.h20) > __DBL_EPSILON__ || fabs(H.h21) > __DBL_EPSILON__)
		return HOMOGRAPHY_TRANSFORMATION;
	if (fabs(H.h00 - 1.) < __DBL_EPSILON__ && fabs(H.h11 - 1.) < __DBL_EPSILON__ &&
			fabs(H.h10) < __DBL_EPSILON__ && fabs(H.h01) < __DBL_EPSILON__) {
		if (fabs(H.h02) > __DBL_EPSILON__ || fabs(H.h12) > __DBL_EPSILON__)
			return SHIFT_TRANSFORMATION;
		return IDENTITY_TRANSFORMATION;
	}
	if (fabs(H.h10 - H.h00  + H.h01 + H.h11) < __DBL_EPSILON__)
		return SIMILARITY_TRANSFORMATION;
	return AFFINE_TRANSFORMATION;
}

void guess_transform_from_seq(sequence *seq, int layer,
		transformation_type *min, transformation_type *max, gboolean excludenull) {
	*min = HOMOGRAPHY_TRANSFORMATION; // highest value
	*max = UNDEFINED_TRANSFORMATION;  // lowest value
	gboolean needs_sel_update = FALSE;

	if (!layer_has_registration(seq, layer)) {
		siril_debug_print("No registration data found in sequence or layer\n");
		return;
	}
	for (int i = 0; i < seq->number; i++){
		transformation_type val = guess_transform_from_H(seq->regparam[layer][i].H);
		//siril_debug_print("Image #%d - transf = %d\n", i+1, val);
		if (*max < val) *max = val;
		if (*min > val) *min = val;
		if (val == NULL_TRANSFORMATION && excludenull) {
			seq->imgparam[i].incl = FALSE;
			needs_sel_update = TRUE;
		}
	}
	if (excludenull && needs_sel_update)
		fix_selnum(seq, FALSE);
}

gboolean check_before_applyreg(struct registration_args *regargs) {
	struct driz_args_t *driz = regargs->driz;
		// check the reference image matrix is not null
	transformation_type checkH = guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H);
	if (checkH == NULL_TRANSFORMATION) {
		siril_log_color_message(_("The reference image has a null matrix and was not previously aligned, choose another one, aborting\n"), "red");
		return FALSE;
	}
	// check the number of dof if -interp=none
	transformation_type min, max;
	guess_transform_from_seq(regargs->seq, regargs->layer, &min, &max, TRUE);
	if (!driz && max > SHIFT_TRANSFORMATION && regargs->interpolation == OPENCV_NONE) {
		siril_log_color_message(_("Applying registration computed with higher degree of freedom (%d) than shift is not allowed when interpolation is set to none, aborting\n"), "red", ((int)max + 1) * 2);
		return FALSE;
	}

	// check the consistency of output images size if -interp=none
	if (!driz && regargs->interpolation == OPENCV_NONE && (regargs->x2upscale || regargs->framing == FRAMING_MAX || regargs->framing == FRAMING_MIN)) {
		siril_log_color_message(_("Applying registration with changes in output image sizeis not allowed when interpolation is set to none , aborting\n"), "red");
		return FALSE;
	}

	// check the consistency of images size if -interp=none
	if (!driz && regargs->interpolation == OPENCV_NONE && regargs->seq->is_variable) {
		siril_log_color_message(_("Applying registration on images with different sizes when interpolation is set to none is not allowed, aborting\n"), "red");
		return FALSE;
	}

	// check that we are not trying to apply identity transform to all the images
	if (max == IDENTITY_TRANSFORMATION) {
		siril_log_color_message(_("Existing registration data is a set of identity matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// check that we are not trying to apply null transform to all the images
	if (max == NULL_TRANSFORMATION || (regargs->seq->selnum <= 1) ) {
		siril_log_color_message(_("Existing registration data is a set of null matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// force -selected if some matrices were null
	if (min == NULL_TRANSFORMATION) {
		siril_log_color_message(_("Some images were not registered, excluding them\n"), "salmon");
		regargs->filters.filter_included = TRUE;
	}

	// cog framing method requires all images to be of same size
	if (regargs->framing == FRAMING_COG && regargs->seq->is_variable) {
		siril_log_color_message(_("Framing method \"cog\" requires all images to be of same size, aborting\n"), "red");
		return FALSE;
	}

	/* compute_framing uses the filtered list of images, so we compute the filter here */
	if (!regargs->filtering_criterion &&
			convert_parsed_filter_to_filter(&regargs->filters,
				regargs->seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter)) {
		return FALSE;
	}
	int nb_frames = compute_nb_filtered_images(regargs->seq,
			regargs->filtering_criterion, regargs->filtering_parameter);
	regargs->new_total = nb_frames;	// to avoid recomputing it later
	gchar *str = describe_filter(regargs->seq, regargs->filtering_criterion,
			regargs->filtering_parameter);
	siril_log_message(str);
	g_free(str);

	// determines the reference homography (including framing shift) and output size
	if (!compute_framing(regargs)) {
		siril_log_color_message(_("Unknown framing method, aborting\n"), "red");
		return FALSE;
	}

	// make sure we apply registration only if the output sequence has a meaningful size
	int rx0 = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry0 = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	if (rx_out < rx0 * MIN_RATIO || ry_out < ry0 * MIN_RATIO) {
		siril_log_color_message(_("The output sequence is too small compared to reference image (too much rotation or little overlap?)\n"), "red");
		siril_log_color_message(_("You should change framing method, aborting\n"), "red");
		return FALSE;
	}

	// cannot use seq_compute_size as rx_out/ry_out are not necessarily consistent with seq->rx/ry
	// rx_out/ry_out already account for 2x upscale if any
	int64_t size = (int64_t) rx_out * ry_out * regargs->seq->nb_layers;
	if (regargs->seq->type == SEQ_SER) {
		size *= regargs->seq->ser_file->byte_pixel_depth;
		size *= nb_frames;
		size += SER_HEADER_LEN;
	} else {
		size *= (get_data_type(regargs->seq->bitpix) == DATA_USHORT) ? sizeof(WORD) : sizeof(float);
		size += FITS_DOUBLE_BLOC_SIZE; // FITS double HDU size
		size *= nb_frames;
	}
	gchar* size_msg = g_format_size_full(size, G_FORMAT_SIZE_IEC_UNITS);
	siril_debug_print("Apply Registration: sequence out size: %s\n", size_msg);
	g_free(size_msg);
	if (test_available_space(size)) {
		siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
		return FALSE;
	}
	return TRUE;
}

