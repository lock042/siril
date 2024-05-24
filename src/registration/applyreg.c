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
static int rx_out, ry_out;

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

static void update_framing(regframe *framing, sequence *seq, int index) {
	if (seq->is_variable) {
		double rx2 = (double)seq->imgparam[index].rx;
		double ry2 = (double)seq->imgparam[index].ry;
		framing->pt[1].x = rx2;
		framing->pt[2].x = rx2;
		framing->pt[2].y = ry2;
		framing->pt[3].y = ry2;
	}
}

static gboolean compute_framing(struct registration_args *regargs) {
	// validity of matrices has already been checked before this call
	// and null matrices have been discarded
	Homography Href = regargs->seq->regparam[regargs->layer][regargs->reference_image].H;
	Homography Hshift = { 0 };
	cvGetEye(&Hshift);
	int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	int x0, y0, rx_0 = rx, ry_0 = ry, n;
	double xmin, xmax, ymin, ymax, cogx, cogy;

	regframe framing = (regframe){(point){0., 0.}, (point){(double)rx, 0.}, (point){(double)rx, (double)ry}, (point){0., (double)ry}};
	gboolean retval = TRUE;

	switch (regargs->framing) {
		case FRAMING_CURRENT:
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				// check we have overlap with the reference
				if (xs[3] < 0. || xs[0] > (double)rx || ys[3] < 0. || ys[0] > (double)ry) {
					siril_log_color_message(_("Image %d has no overlap with the reference\n"), "red", i + 1);
					retval = FALSE;
				}
			}
			break;
		case FRAMING_MAX:
			break;
		case FRAMING_MIN:
			xmin = -DBL_MAX;
			xmax = DBL_MAX;
			ymin = -DBL_MAX;
			ymax = DBL_MAX;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
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
			if (rx_0 < 0 || ry_0 < 0) {
				siril_log_color_message(_("The intersection of all images is null or negative\n"), "red");
				retval = FALSE;
			}
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
				int rx2 = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx;
				int ry2 = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry;
				double currcogx = (double)rx2 * 0.5;
				double currcogy = (double)ry2 * 0.5;
				cvTransfPoint(&currcogx, &currcogy,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
				cogx += currcogx;
				cogy += currcogy;
				n++;
			}
			cogx /= (double)n;
			cogy /= (double)n;
			x0 = (int)(cogx - (double)rx * 0.5);
			y0 = (int)(cogy - (double)ry * 0.5);
			siril_log_message(_("Framing: Shift from reference origin: %d, %d\n"), x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			Homography newHref = { 0 };
			cvMultH(Href, Hshift, &newHref);
		// we now check overlaps
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, newHref, 1.);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				// check we have overlap with the reference
				if (xs[3] < 0. || xs[0] > (double)rx || ys[3] < 0. || ys[0] > (double)ry) {
					siril_log_color_message(_("Image %d has no overlap with the reference\n"), "red", i + 1);
					retval = FALSE;
				}
			}
			break;
		default:
			return FALSE;
	}
	cvMultH(Href, Hshift, &Htransf);
	if (regargs->driz) {
		rx_out = rx_0 * regargs->driz->scale;
		ry_out = ry_0 * regargs->driz->scale;
	} else {
		rx_out = rx_0 * ((regargs->x2upscale) ? 2. : 1.);
		ry_out = ry_0 * ((regargs->x2upscale) ? 2. : 1.);
	}
	return retval;
}

// For framing max, we don't want to export the image with black borders
// From the homographies and the input image size
// we compute:
// - H transf, the warping transformation in place
// - H shift, the residual shift wrt to ref image
// - the size of the transformed image after applying the transformation
static void compute_Hmax(Homography *Himg, Homography *Href, int src_rx_in, int src_ry_in, double scale, Homography *H, Homography *Hshift, int *dst_rx_out, int *dst_ry_out) {
	*dst_rx_out = 0;
	*dst_ry_out = 0;
	regframe framing = { 0 };
	framing = (regframe){(point){0., 0.}, (point){(double)src_rx_in, 0.}, (point){(double)src_rx_in, (double)src_ry_in}, (point){0., (double)src_ry_in}};
	double xmin, xmax, ymin, ymax;
	xmin = DBL_MAX;
	xmax = -DBL_MAX;
	ymin = DBL_MAX;
	ymax = -DBL_MAX;
	for (int j = 0; j < 4; j++) {
		cvTransfPoint(&framing.pt[j].x, &framing.pt[j].y, *Himg, *Href, scale);
		if (xmin > framing.pt[j].x) xmin = framing.pt[j].x;
		if (ymin > framing.pt[j].y) ymin = framing.pt[j].y;
		if (xmax < framing.pt[j].x) xmax = framing.pt[j].x;
		if (ymax < framing.pt[j].y) ymax = framing.pt[j].y;
		siril_debug_print("Point #%d: %3.2f %3.2f\n", j, framing.pt[j].x, framing.pt[j].y);
	}
	*dst_rx_out = (int)(ceil(xmax) - floor(xmin)) + 1;
	*dst_ry_out = (int)(ceil(ymax) - floor(ymin)) + 1;
	Hshift->h02 = (double)floor(xmin);
	Hshift->h12 = (double)floor(ymin);
	cvTransfH(*Himg, *Href, H);
	// the shift matrix is at the final scale, while the transformation matrix
	// is still at the orginal scale (will be upscaled in cvTransformImage)
	H->h02 -= Hshift->h02 / scale;
	H->h12 -= Hshift->h12 / scale;
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
	if (!regargs->driz && fit.naxes[2] == 1 && fit.keywords.bayer_pattern[0] != '\0')
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
	float scale = 1.f;
	if (regargs->driz)
		scale = driz->scale;
	else if (regargs->x2upscale)
		scale = 2.f;

	Homography H = { 0 };
	Homography Himg = { 0 };
	Homography Hs = { 0 };
	cvGetEye(&Hs);
	int dst_rx = rx_out;
	int dst_ry = ry_out;
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
		siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
	}

	// Composing transformation wrt reference image
	Himg = regargs->seq->regparam[regargs->layer][in_index].H;
	if (guess_transform_from_H(Himg) == NULL_TRANSFORMATION)
		return 1; // in case H is null and -selected was not passed

	if (regargs->framing != FRAMING_MAX)
		cvTransfH(Himg, Htransf, &H);
	else {
		compute_Hmax(&Himg, &Htransf, fit->rx, fit->ry, scale, &H, &Hs, &dst_rx, &dst_ry);
	}

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
		p->threads = threads;

		map_image_coordinates_h(fit, H, p->pixmap, dst_ry, driz->scale, threads);
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
		fits out = { 0 };
		copyfits(fit, &out, CP_FORMAT, -1);
		// copy the DATE_OBS
		out.keywords.date_obs = g_date_time_ref(fit->keywords.date_obs);
		out.rx = out.naxes[0] = dst_rx;
		out.ry = out.naxes[1] = dst_ry;
		out.naxes[2] = driz->is_bayer ? 3 : 1;
		size_t chansize = out.rx * out.ry * sizeof(float);
		out.fdata = calloc(out.naxes[2] * chansize, 1);
		out.fpdata[0] = out.fdata;
		out.fpdata[1] = out.naxes[2] == 1 ? out.fdata : out.fdata + chansize;
		out.fpdata[2] = out.naxes[2] == 1 ? out.fdata : out.fdata + 2 * chansize;
		p->output_data = &out;

		// Set up the output_counts fits to store pixel hit counts
		fits *output_counts = calloc(1, sizeof(fits));
		copyfits(&out, output_counts, CP_FORMAT, -1);
		output_counts->fdata = calloc(output_counts->rx * output_counts->ry * output_counts->naxes[2], sizeof(float));
		p->output_counts = output_counts;

		p->weights = driz->flat;

		if (dobox(p)) { // Do the drizzle
			siril_log_color_message("s\n", p->error->last_message);
			return 1;
		}
		clearfits(fit);
		// copy the DATE_OBS
		copyfits(&out, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		fit->keywords.date_obs = g_date_time_ref(out.keywords.date_obs);
		clearfits(&out);
		if (args->seq->type == SEQ_SER || com.pref.force_16bit) {
			fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, fit->rx * fit->ry * fit->naxes[2]), DATA_USHORT);
		}
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

		// Get rid of the output_counts image, no longer required
		clearfits(output_counts);
		free(output_counts);
		output_counts = NULL;
	}
	else {
		if (regargs->interpolation <= OPENCV_LANCZOS4) {
			if (cvTransformImage(fit, dst_rx, dst_ry, H, regargs->x2upscale, regargs->interpolation, regargs->clamp)) {
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
	regargs->imgparam[out_index].rx = dst_rx;
	regargs->imgparam[out_index].ry = dst_ry;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
	regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
	regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
	regargs->regparam[out_index].H = Hs;

	if (regargs->driz || regargs->x2upscale) {
		fit->keywords.pixel_size_x /= scale;
		fit->keywords.pixel_size_y /= scale;
		regargs->regparam[out_index].fwhm *= scale;
		regargs->regparam[out_index].weighted_fwhm *= scale;
	}

	sadata->success[out_index] = 1;
	return 0;
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
		args->compute_mem_limits_hook = apply_drz_compute_mem_limits;
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
		args->compute_mem_limits_hook = apply_reg_compute_mem_limits;
		args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	}
	args->prepare_hook = apply_reg_prepare_hook;
	args->finalize_hook = apply_reg_finalize_hook;
	args->filtering_criterion = regargs->filtering_criterion;
	args->filtering_parameter = regargs->filtering_parameter;
	args->nb_filtered_images = regargs->new_total;
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
	seq.number = args->new_total;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = (args->driz && args->driz->is_bayer) ? 3 : args->seq->nb_layers;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number-1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	seq.is_variable = check_seq_is_variable(&seq);
	if (!seq.is_variable) {
		seq.rx = args->seq->rx;
		seq.ry = args->seq->ry;
	}
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
		siril_log_color_message(_("Unselect the images generating the error or change framing method to max\n"), "red");
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

