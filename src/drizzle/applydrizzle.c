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
#include "algos/siril_wcs.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/preprocess.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"
#include "algos/demosaicing.h"

#include "drizzle/driz_portability.h"
#include "drizzle/cdrizzleutil.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzlebox.h"

#include "opencv/opencv.h"

#define MIN_RATIO 0.1 // minimum fraction of ref image dimensions to validate output sequence is worth creating

static void create_output_sequence_for_apply_driz(struct driz_args_t *args);
static int new_ref_index = -1;

static Homography Htransf = {0};
static int rx_out, ry_out, ry_out_unscaled;

static regdata *apply_driz_get_current_regdata(struct driz_args_t *driz) {
	regdata *current_regdata;
	if (driz->seq->regparam[0]) {
		siril_log_message(
				_("Applying existing registration from layer 0 to transform the images\n"));
		current_regdata = driz->seq->regparam[0];
	} else {
		siril_log_message(
				_("No registration data exists for this layer\n"));
		return NULL;
	}
	return current_regdata;
}

// TODO: once WCS drizzle mapping is reinstated, the framing calculation needs to be written properly
/*
static gboolean driz_compute_wcs_framing(struct driz_args_t *driz) {
	int rx = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].rx : driz->seq->rx;
	int ry = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].ry : driz->seq->ry;
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
	rx_out = rx_0 * driz->scale;
	ry_out = ry_0 * driz->scale;
	return TRUE;
}
*/

// confirmed generates the correct Homography
static gboolean driz_compute_framing(struct driz_args_t *driz) {
	// validity of matrices has already been checked before this call
	// and null matrices have been discarded
	Homography Href = driz->seq->regparam[0][driz->reference_image].H;
	Homography Hshift = {0};
	cvGetEye(&Hshift);
	int rx = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].rx : driz->seq->rx;
	int ry = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].ry : driz->seq->ry;
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

	switch (driz->framing) {
		case FRAMING_CURRENT:
			break;
		case FRAMING_MAX:
			xmin = DBL_MAX;
			xmax = -DBL_MAX;
			ymin = DBL_MAX;
			ymax = -DBL_MAX;
			for (int i = 0; i < driz->seq->number; i++) {
				if (!driz->filtering_criterion(driz->seq, i, driz->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				if (driz->seq->is_variable) {
					double rx2 = (double)driz->seq->imgparam[i].rx;
					double ry2 = (double)driz->seq->imgparam[i].ry;
					current_framing.pt[1].x = rx2;
					current_framing.pt[2].x = rx2;
					current_framing.pt[2].y = ry2;
					current_framing.pt[3].y = ry2;
				}
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y, driz->seq->regparam[0][i].H, Href);
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
			for (int i = 0; i < driz->seq->number; i++) {
				if (!driz->filtering_criterion(driz->seq, i, driz->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				if (driz->seq->is_variable) {
					double rx2 = (double)driz->seq->imgparam[i].rx;
					double ry2 = (double)driz->seq->imgparam[i].ry;
					current_framing.pt[1].x = rx2;
					current_framing.pt[2].x = rx2;
					current_framing.pt[2].y = ry2;
					current_framing.pt[3].y = ry2;
				}
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,driz->seq->regparam[0][i].H, Href);
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
			for (int i = 0; i < driz->seq->number; i++) {
				if (!driz->filtering_criterion(driz->seq, i, driz->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = {0};
				memcpy(&current_framing, &framing, sizeof(regframe));
				double currcogx = 0., currcogy = 0.;
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,driz->seq->regparam[0][i].H, Href);
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
			siril_log_message(_("Framing for drizzle: shift from reference origin: %d, %d\n"), x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		default:
			return FALSE;
	}
	cvMultH(Href, Hshift, &Htransf);
	rx_out = rx_0 * driz->scale;
	ry_out = ry_0 * driz->scale;
	ry_out_unscaled = ry_0;
	return TRUE;
}

int apply_drz_prepare_results(struct generic_seq_args *args) {
	struct driz_args_t *driz = args->user;

	// allocate destination sequence data
	driz->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
	driz->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
	if (!driz->imgparam  || !driz->regparam) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	if (seq_prepare_hook(args))
		return 1;

	driz->success = calloc(args->nb_filtered_images, sizeof(BYTE));
	if (!driz->success) {
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
	struct driz_args_t *driz = args->user;
	int number_of_outputs = 1;
	// we call the generic prepare twice with different prefixes
	args->new_seq_prefix = driz->prefix;
	if (apply_drz_prepare_results(args))
		return 1;
	// but we copy the result between each call
	driz->new_ser_drz = args->new_ser;
	driz->new_fitseq_drz = args->new_fitseq;

	if (driz->keep_counts) {
		args->new_seq_prefix = "oc_"; // This is OK here, it does not get freed in the
		// end_generic_sequence later
		if (apply_drz_prepare_results(args))
			return 1;
		driz->new_ser_pxcnt = args->new_ser;
		driz->new_fitseq_pxcnt = args->new_fitseq;

		args->new_seq_prefix = driz->prefix; // Put it back so it gets loaded on completion
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

/* reads the image and apply existing transformation */
int apply_drz_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct driz_args_t *driz = args->user;
//	struct wcsprm *refwcs = driz->refwcs;
	/* Set up the per-image drizzle parameters */
	struct driz_param_t *p = calloc(1, sizeof(struct driz_param_t));

	driz_param_init(p);
	p->kernel = driz->kernel;
	p->driz = driz;
	p->error = malloc(sizeof(struct driz_error_t));
	p->scale = driz->scale;
	p->pixel_fraction = driz->pixel_fraction;
	p->cfa = driz->cfa;

	// Set bounds. TODO: make this account for a selection
	p->xmin = p->ymin = 0;
	p->xmax = fit->rx - 1;
	p->ymax = fit->ry - 1;

	Homography H = { 0 };
	Homography Himg = { 0 };
	if (!driz->use_wcs) {
		p->current_regdata = apply_driz_get_current_regdata(driz);
		if (!p->current_regdata) {
			free(p);
			return -2;
		}
	}

	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes
	if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
		siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
	}

	// Composing transformation wrt reference image
	if (driz->use_wcs) {
		if (!fit->wcslib) {
			siril_log_color_message(_("Error: drizzle configured to use WCS transforms but this frame "
						"is not plate solved. Ensure all frames are plate solved, e.g. by running "
						"seqplatesolve.\n"), "red");
			return 1;
		}
	} else {
		Himg = driz->seq->regparam[0][in_index].H;
		if (guess_transform_from_H(Himg) == NULL_TRANSFORMATION)
			return 1; // in case H is null and -selected was not passed
		cvTransfH(Himg, Htransf, &H);
	}
	/* Populate the mapping array. This maps pixels from the current frame to
	 * the reference frame. Either a Homography mapping can be used based on
	 * image registration or a WCS mapping can be used based on plate solving */
	p->pixmap = calloc(1, sizeof(imgmap_t));
	p->pixmap->rx = fit->rx;
	p->pixmap->ry = fit->ry;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
//	if (p->driz->use_wcs) {
//		map_image_coordinates_wcs(fit->rx, fit->ry, fit->wcslib, refwcs, p->pixmap, driz->scale);
//	} else {
		map_image_coordinates_h(fit, H, p->pixmap, ry_out_unscaled, driz->scale);
//	}
	if (!p->pixmap->pixmap) {
		siril_log_color_message(_("Error generating mapping array.\n"), "red");
		return 1;
	}
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, _("Remapping"));
	/* Populate the data fits to be drizzled */
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

	gettimeofday(&t_start, NULL);
	if (dobox(p)) // Do the drizzle
		return 1;
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, _("Drizzle"));

	// Copy driz->output_data to fit so that it is saved as the output sequence frame
	clearfits(fit);
	copyfits(&out, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	clearfits(&out);

	free(p->pixmap->pixmap);
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

	if (in_index == driz->reference_image)
		new_ref_index = out_index; // keeping track of the new ref index in output sequence

	driz->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	driz->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	driz->imgparam[out_index].rx = rx_out;
	driz->imgparam[out_index].ry = ry_out;
	if (!driz->use_wcs) {
		driz->regparam[out_index].fwhm = p->current_regdata[in_index].fwhm;
		driz->regparam[out_index].weighted_fwhm = p->current_regdata[in_index].weighted_fwhm;
		driz->regparam[out_index].roundness = p->current_regdata[in_index].roundness;
		driz->regparam[out_index].background_lvl = p->current_regdata[in_index].background_lvl;
		driz->regparam[out_index].number_of_stars = p->current_regdata[in_index].number_of_stars;
		cvGetEye(&driz->regparam[out_index].H);
	}

	// Compensate metadata for any change in scale
	fit->pixel_size_x /= p->scale;
	fit->pixel_size_y /= p->scale;
	driz->regparam[out_index].fwhm *= p->scale;
	driz->regparam[out_index].weighted_fwhm *= p->scale;
	driz->imgparam[out_index].rx *= p->scale;
	driz->imgparam[out_index].ry *= p->scale;
	driz->success[out_index] = 1;
	driz->new_total++; // is this thread safe?
	return 0;
}

static int apply_drz_save_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct driz_args_t *driz = (struct driz_args_t*) args->user;
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
		if (double_data->output)
			retval1 = ser_write_frame_from_fit(driz->new_ser_drz, double_data->output, out_index);
		if (double_data->pixel_count)
			retval2 = ser_write_frame_from_fit(driz->new_ser_pxcnt, double_data->pixel_count, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		if (double_data->output)
			retval1 = fitseq_write_image(driz->new_fitseq_drz, double_data->output, out_index);
		if (double_data->pixel_count)
			retval2 = fitseq_write_image(driz->new_fitseq_pxcnt, double_data->pixel_count, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, driz->prefix, in_index);
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
		gchar *ocprefix = g_strdup_printf("oc_%s", driz->prefix);
		dest = fit_sequence_get_image_filename_prefixed(args->seq, ocprefix, in_index);
		g_free(ocprefix);
		if (double_data->pixel_count) {
			// No 16-bit conversion here, the output_counts is always float.
			retval2 = savefits(dest, double_data->pixel_count);
		}
		free(dest);
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
	struct driz_args_t *driz = (struct driz_args_t*) args->user;

	// We deal with the pixel_count sequence first
	args->new_ser = driz->new_ser_pxcnt;
	args->new_fitseq = driz->new_fitseq_pxcnt;
	int retval = 0;
	if (args->new_ser || args->new_fitseq) {
		retval = seq_finalize_hook(args);
	}
	driz->new_ser_pxcnt = NULL;
	driz->new_fitseq_pxcnt = NULL;

	// We deal with the drizzled sequence next
	args->new_ser = driz->new_ser_drz;
	args->new_fitseq = driz->new_fitseq_drz;

	int failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!driz->success[i])
				failed++;
		driz->new_total = args->nb_filtered_images - failed;
		if (failed) {
			// driz->imgparam and driz->regparam may have holes caused by images
			// that failed to be registered - compact them
			for (int i = 0, j = 0; i < driz->new_total; i++, j++) {
				while (!driz->success[j] && j < args->nb_filtered_images) j++;
				g_assert(driz->success[j]);
				if (i != j) {
					driz->imgparam[i] = driz->imgparam[j];
					driz->regparam[i] = driz->regparam[j];
				}
			}
		}
		seq_finalize_hook(args);
	} else {
		driz->new_total = 0;

		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		} else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		} else if (args->seq->type == SEQ_REGULAR) {
			remove_prefixed_sequence_files(driz->seq, driz->prefix);
		}
	}

	if (driz->success) free(driz->success);
	//	Do not free driz here as it will be needed for blot and second drizzle
	args->user = NULL;

	if (!args->retval) {
		siril_log_message(_("Applying drizzle completed.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d exported.\n"), "green", failed, driz->new_total);

		g_free(str);
		if (!(args->seq->type == SEQ_INTERNAL)) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_apply_driz(driz);
			// will be loaded in the idle function if (load_new_sequence)
			args->load_new_sequence = driz->load_new_sequence; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Drizzle aborted.\n"));
	}

	if (driz->new_total == 0)
		return 1;
	driz->new_ser_drz = NULL;
	driz->new_fitseq_drz = NULL;

	return retval;
}

int apply_drz_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	struct driz_args_t *driz = args->user;
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

gboolean check_before_applydrizzle(struct driz_args_t *driz);

int apply_drizzle(struct driz_args_t *driz) {
	struct generic_seq_args *args = create_default_seqargs(driz->seq);
	set_progress_bar_data(_("Initializing drizzle data..."), PROGRESS_PULSATE);

	if (!check_before_applydrizzle(driz)) {
		free(args);
		return -1;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);

	args->seq = driz->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = driz->seq->selnum;
	args->compute_mem_limits_hook = apply_drz_compute_mem_limits;
	args->prepare_hook = apply_drz_prepare_hook;
	args->image_hook = apply_drz_image_hook;
	args->save_hook = apply_drz_save_hook;
	args->finalize_hook = apply_drz_finalize_hook;
	args->description = _("Apply drizzle");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(driz->prefix);
	args->force_float = FALSE;
	args->user = driz;

	driz_param_dump(driz); // Print some info to the log
	fits fit = { 0 };

	/* preparing reference data from reference fit and making sanity checks*/
	/* fit will now hold the reference frame */
	if (seq_read_frame_metadata(args->seq, driz->reference_image, &fit)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[0] = NULL;
		free(args);
		return 1;
	}

	if (driz->use_flats) {
		fits reffit = { 0 };
		GtkEntry *entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		const gchar *flat_filename = gtk_entry_get_text(entry);
		gchar *error = NULL;
		int status;
		gchar *expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
		if (status) {
			error = _("NOT USING FLAT: could not parse the expression");
			driz->use_flats = FALSE;
		} else {
			free(expression);
			if (flat_filename[0] == '\0') {
				siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
				return 1;
			} else {
				set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
				driz->flat = calloc(1, sizeof(fits));
				if (!readfits(flat_filename, driz->flat, NULL, TRUE)) {
					if (driz->flat->naxes[2] != fit.naxes[2]) {
						error = _("NOT USING FLAT: number of channels is different");
					} else if (driz->flat->naxes[0] != fit.naxes[0] ||
							driz->flat->naxes[1] != fit.naxes[1]) {
						error = _("NOT USING FLAT: image dimensions are different");
					} else {
						// no need to deal with bitdepth conversion as readfits has already forced conversion to float
						siril_log_message(_("Master flat read for use as initial pixel weight\n"));
					}

				} else error = _("NOT USING FLAT: cannot open the file");
				if (error) {
					siril_log_color_message("%s\n", "red", error);
					set_progress_bar_data(error, PROGRESS_DONE);
					if (driz->flat) {
						clearfits(driz->flat);
						free(driz->flat);
					}
					return 1;
				}
				if (driz->use_bias) {
					entry = GTK_ENTRY(lookup_widget("offsetname_entry"));
					const gchar *offset_filename = gtk_entry_get_text(entry);
					if (offset_filename[0] == '\0') {
						siril_log_message(_("Error: no master bias specified in the preprocessing tab.\n"));
						return 1;
					} else if (offset_filename[0] == '=') { // offset is specified as a level not a file
						set_progress_bar_data(_("Checking offset level..."), PROGRESS_NONE);
						int offsetlevel = evaluateoffsetlevel(offset_filename + 1, &fit);
						if (!offsetlevel) {
							error = _("NOT USING OFFSET: the offset value could not be parsed");
							driz->use_bias = FALSE;
						} else {
							siril_log_message(_("Synthetic offset: Level = %d\n"),offsetlevel);
							int maxlevel = (fit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
							if ((offsetlevel > maxlevel) || (offsetlevel < -maxlevel) ) {   // not excluding all neg values here to allow defining a pedestal
								error = _("NOT USING OFFSET: the offset value is not consistent with image bitdepth");
								driz->use_bias = FALSE;
							} else {
								float bias_level = (float)offsetlevel;
								bias_level *= (fit.orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE; //converting to [0 1] to use with soper
								soper(driz->flat, bias_level, OPER_SUB, TRUE);
								siril_log_message(_("Bias level subtracted from master flat for use as initial pixel weighting\n"));
							}
						}
					} else {
						int status;
						expression = path_parse(&reffit, offset_filename, PATHPARSE_MODE_READ, &status);
						if (status) {
							error = _("NOT USING OFFSET: could not parse the expression\n");
						} else {
							set_progress_bar_data(_("Opening master bias image..."), PROGRESS_NONE);
							fits bias = { 0 };
							if (!readfits(offset_filename, &bias, NULL, TRUE)) {
								if (bias.naxes[2] != fit.naxes[2]) {
									error = _("NOT USING BIAS: number of channels is different");
								} else if (bias.naxes[0] != fit.naxes[0] ||
										bias.naxes[1] != fit.naxes[1]) {
									error = _("NOT USING BIAS: image dimensions are different");
								} else {
									imoper(driz->flat, &bias, OPER_SUB, TRUE);
									siril_log_message(_("Master bias subtracted from master flat for use as initial pixel weighting\n"));
									// no need to deal with bitdepth conversion as flat is just a division (unlike darks which need to be on same scale)
								}
							} else error = _("NOT USING BIAS: cannot open the file");
							clearfits(&bias);
						}
					}
					if (error) {
						siril_log_color_message("%s\n", "red", error);
						set_progress_bar_data(error, PROGRESS_DONE);
						if (driz->flat) {
							clearfits(driz->flat);
							free(driz->flat);
						}
						return 1;
					}
				}
			}
		}
	}

/*	if (driz->use_wcs) {
		if (!fit.wcslib) {
			siril_log_color_message(_("Error: reference image is not plate solved. Unable to drizzle using WCS data.\n"), "red");
			return 1;
		}
		// Set the reference WCS data
		int copy_status;
		driz->refwcs = wcs_deepcopy(fit.wcslib, &copy_status);
		if (copy_status) {
			siril_log_color_message(_("Error: failed to set the reference WCS.\n"), "red");
			return 1;
		}
		if (driz->refwcs->lin.dispre) {
			// Disabling distortion in the reference frame so that the distorted
			// input images are mapped to a flat output reference
			free(driz->refwcs->lin.dispre);
			driz->refwcs->lin.dispre = NULL;
		}
	}
*/

	driz->is_bayer = (fit.bayer_pattern[0] != '\0'); // If there is a CFA pattern we need to CFA drizzle
	sensor_pattern pattern = com.pref.debayer.use_bayer_header ? get_cfa_pattern_index_from_string(fit.bayer_pattern) : com.pref.debayer.bayer_pattern;
	if (driz->is_bayer) {
		driz->cfa = get_cfa_from_pattern(pattern);
		if (!driz->cfa) // if fit.bayer_pattern exists and get_cfa_from_pattern returns NULL then there is a problem.
			return 1;
	}

	start_in_new_thread(generic_sequence_worker, args);
	return 0;
}

static void create_output_sequence_for_apply_driz(struct driz_args_t *args) {
	sequence seq = { 0 };
	initialize_sequence(&seq, TRUE);

	/* we are not interested in the whole path */
	gchar *seqname = g_path_get_basename(args->seq->seqname);
	char *rseqname = malloc(strlen(args->prefix) + strlen(seqname) + 5);
	sprintf(rseqname, "%s%s.seq", args->prefix, seqname);
	g_free(seqname);
	g_unlink(rseqname);	// remove previous to overwrite
	args->new_seq_name = remove_ext_from_filename(rseqname);
	free(rseqname);
	seq.seqname = strdup(args->new_seq_name);
	gchar *tmpbuf = g_strdup_printf("oc_%s", seq.seqname);
	siril_debug_print("seq name: %s\npix_cnt seq name: %s\n", args->seq->seqname, tmpbuf);
	seq.pixcnt_seqname = strdup(tmpbuf);
	g_free(tmpbuf);
	seq.number = args->new_total;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = args->is_bayer ? 3 : args->seq->nb_layers;
	seq.rx = rx_out;
	seq.ry = ry_out;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[0] = args->regparam;
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

gboolean check_before_applydrizzle(struct driz_args_t *driz) {
	// check the reference image matrix is not null
	if (!driz->use_wcs) {
		if (!(driz->seq && (driz->seq->regparam[0]))) {
		siril_log_color_message(_("Error: registration parameters not found, aborting\n"), "red");
		return FALSE;
		}

		transformation_type checkH = guess_transform_from_H(driz->seq->regparam[0][driz->seq->reference_image].H);
		if (checkH == NULL_TRANSFORMATION) {
			siril_log_color_message(_("The reference image has a null matrix and was not previously aligned, choose another one, aborting\n"), "red");
			return FALSE;
		}

		transformation_type min, max;
		guess_transform_from_seq(driz->seq, 0, &min, &max, TRUE);

		// check that we are not trying to apply identity transform to all the images
		if (max == IDENTITY_TRANSFORMATION) {
			siril_log_color_message(_("Existing registration data is a set of identity matrices, no transformation would be applied, aborting\n"), "red");
			return FALSE;
		}

		// check that we are not trying to apply null transform to all the images
		if (max == NULL_TRANSFORMATION || (driz->seq->selnum <= 1) ) {
			siril_log_color_message(_("Existing registration data is a set of null matrices, no transformation would be applied, aborting\n"), "red");
			return FALSE;
		}

		// force -selected if some matrices were null
		if (min == NULL_TRANSFORMATION) {
			siril_log_color_message(_("Some images were not registered, excluding them\n"), "salmon");
			driz->filters.filter_included = TRUE;
		}
	}

	// cog framing method requires all images to be of same size
	if (driz->framing == FRAMING_COG && driz->seq->is_variable) {
		siril_log_color_message(_("Framing method \"cog\" requires all images to be of same size, aborting\n"), "red");
		return FALSE;
	}

	// compute_framing uses the filtered list of images, so we compute the filter here
	if (!driz->filtering_criterion &&
			convert_parsed_filter_to_filter(&driz->filters,
				driz->seq, &driz->filtering_criterion,
				&driz->filtering_parameter)) {
		return FALSE;
	}
	int nb_frames = compute_nb_filtered_images(driz->seq,
			driz->filtering_criterion, driz->filtering_parameter);
	driz->new_total = nb_frames;	// to avoid recomputing it later
	gchar *str = describe_filter(driz->seq, driz->filtering_criterion,
			driz->filtering_parameter);
	siril_log_message(str);
	g_free(str);

	// determines the reference homography (including framing shift) and output size
	// This line can be restored once the WCS drizzle mapping is reintroduced following the mosaic work
	//	int ret = (driz->use_wcs) ? (driz_compute_wcs_framing(driz)) : (driz_compute_framing(driz));

	int ret = driz_compute_framing(driz);
	if (!ret) {
		siril_log_color_message(_("Unknown framing method, aborting\n"), "red");
		return FALSE;
	}

	// make sure we apply registration only if the output sequence has a meaningful size
	int rx0 = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].rx : driz->seq->rx;
	int ry0 = (driz->seq->is_variable) ? driz->seq->imgparam[driz->reference_image].ry : driz->seq->ry;
	if (rx_out < rx0 * MIN_RATIO || ry_out < ry0 * MIN_RATIO) {
		siril_log_color_message(_("The output sequence is too small compared to reference image (too much rotation or little overlap?)\n"), "red");
		siril_log_color_message(_("You should change framing method, aborting\n"), "red");
		return FALSE;
	}

	// cannot use seq_compute_size as rx_out/ry_out are not necessarily consistent with seq->rx/ry
	// rx_out/ry_out already account for scale
	int64_t size = (int64_t) rx_out * ry_out * driz->seq->nb_layers;
	if (driz->seq->type == SEQ_SER) {
		size *= driz->seq->ser_file->byte_pixel_depth;
		size *= nb_frames;
		size += SER_HEADER_LEN;
	} else {
		size *= (get_data_type(driz->seq->bitpix) == DATA_USHORT) ? sizeof(WORD) : sizeof(float);
		size += 5760; // FITS double HDU size
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

