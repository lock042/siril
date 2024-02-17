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
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"

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
	rx_out = rx_0 * ((regargs->x2upscale) ? 2. : 1.);
	ry_out = ry_0 * ((regargs->x2upscale) ? 2. : 1.);
	return TRUE;
}

int apply_drz_prepare_results(struct generic_seq_args *args) {
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

int apply_drz_prepare_hook(struct generic_seq_args *args) {
	struct driz_param_t *driz = args->user;
	struct registration_args *regargs = driz->regargs;

	fits fit = { 0 };

	/* preparing reference data from reference fit and making sanity checks*/
	driz->current_regdata = apply_reg_get_current_regdata(regargs);
	if (!driz->current_regdata) return -2;


	if (seq_read_frame_metadata(args->seq, regargs->reference_image, &fit)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(driz->current_regdata);
		return 1;
	}
	if (fit.naxes[2] == 1 && fit.bayer_pattern[0] != '\0' && !driz->is_bayer)
		siril_log_color_message(_("Applying drizzle on a sequence opened as CFA without setting the Bayer flag is a bad idea.\n"), "red");
	if (driz->use_wcs && !fit->wcslib) {
		// TODO: Attempt to platesolve the reference image
		// Code goes here, then try again to see if we have a viable struct wcslib...
		if (driz->use_wcs && !fit->wcslib) {
			siril_log_color_message(_("Error: platesolver failed. Unable to drizzle using WCS data.\n"), "red");
			return 1;
		}
	}
	driz->ref_wcs = fit->wcslib;
}

/* reads the image and apply existing transformation */
int apply_drz_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct driz_param_t *driz = args->user;
	struct registration_args *regargs = driz->regargs;
	struct wcsprm *refwcs = driz->refwcs;

	Homography H = { 0 };
	Homography Himg = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
		siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
	}

	// Composing transformation wrt reference image
	if (driz->use_wcs) {
		// TODO: check the image is plate solved, solve it if not
	} else {
		Himg = regargs->seq->regparam[regargs->layer][in_index].H;
		if (guess_transform_from_H(Himg) == NULL_TRANSFORMATION)
			return 1; // in case H is null and -selected was not passed
		cvTransfH(Himg, Htransf, &H);
	}

	/* Populate the mapping array. This maps pixels from the current frame to
	 * the reference frame. Either a Homography mapping can be used based on
	 * image registration or a WCS mapping can be used based on plate solving */

	driz->pixmap = malloc(sizeof(imgmap_t));
	driz->pixmap->rx = fit->rx;
	driz->pixmap->ry = fit->ry;
	if (driz->use_wcs) {
		map_image_coordinates_wcs(&fit, wcs_ref, driz->pixmap);
	} else {
		map_image_coordinates_h(&fit, H, pixmap);
	}

	/* NOTE: on the first pass there is no weights file, everything is evenly
	 *       weighted. This will be used to remove outliers after drz_medstack.
	 *       So we don't need to initialize driz->weights here */

	if (dobox(driz)) // Do the drizzle
		return 1;

	// TODO: copy driz->output_data to fit so that it is saved as the output sequence frame
	// TODO: do we still need the mapping file for the blot or the second drizzle or is it
	//       cheaper to recalculate it (if needed) than to save / load it?
	// TODO: do we need to save driz->output_counts? Or just discard it once dobox() is done?

	// Set up the output_counts fits to store pixel hit counts
	fits output_counts = { 0 };
	output_counts->rx = fit->rx;
	output_counts->ry = fit->ry;
	output_counts->type = DATA_FLOAT;
	output_counts->fdata = calloc(output_counts->rx * output_counts->ry * sizeof(float));
	driz->output_counts = &output_counts;

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

	// Compensate metadata for any change in scale
	fit->pixel_size_x /= driz->scale;
	fit->pixel_size_y /= driz->scale;
	regargs->regparam[out_index].fwhm *= driz->scale;
	regargs->regparam[out_index].weighted_fwhm *= driz->scale;
	regargs->imgparam[out_index].rx *= driz->scale;
	regargs->imgparam[out_index].ry *= driz->scale;

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
	struct driz_param_t *driz = args->user;
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, driz->scale, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
	int float_multiplier = (is_float) ? 1 : 2;
	int is_color = args->seq->nb_layers == 3;
	int is_bayer = driz->is_bayer;
	MB_per_scaled_image *= is_float ? 1 : 2; // Output is always float
	if (is_bayer)
		MB_per_scaled_image *= 3;
	unsigned int MB_per_float_image = MB_per_orig_image * float_multiplier;
	unsigned int MB_per_float_channel = is_color ? MB_per_float_image / 3 : MB_per_float_image;
	unsigned int MB_per_double_channel = MB_per_float_channel * 2;
	/* The drizzle memory consumption is:
		* the original image
		* Two 2 * rx * ry * double arrays for computing mapping
		  (the mapping file reallocs one of these so never exceeds that amount)
		* the weights file (1 * input channel float data)
		* the transformed image, including scaling factor if required
		* the output counts image, the same size as the scaled image
		*/
	unsigned int required = MB_per_orig_image + 2 * MB_per_double_channel + MB_per_float_channel + 2 * MB_per_scaled_image;
	// If interpolation clamping is set, 2x additional Mats of the same format
	// as the original image are required
	struct registration_args *regargs = sadata->regargs;

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

int apply_drizzle(/*TODO: what parameters here? */) {
	struct generic_seq_args *args = create_default_seqargs(p->seq);

	if (!check_before_applydrizzle(regargs)) {
		free(args);
		return -1;
	}

	args->filtering_criterion = regargs->filtering_criterion;
	args->filtering_parameter = regargs->filtering_parameter;
	args->nb_filtered_images = regargs->new_total;
	args->compute_mem_limits_hook = apply_drz_compute_mem_limits;
	args->prepare_hook = apply_drz_prepare_hook;
	args->image_hook = apply_drz_image_hook;
	args->finalize_hook = apply_drz_finalize_hook;
	args->idle_function = apply_drz_medstack; /* When this sequence finishes,
		it kicks off another one to perform a median stack of the drizzled
		images. */
	args->stop_on_error = FALSE;
	args->description = _("Apply drizzle");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
//	args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	args->new_seq_prefix = g_strdup("drz");
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct driz_param_t drizdata = calloc(1, sizeof(struct driz_param_t));
	if (!drizdata) {
		free(args);
		return -1;
	}
	driz_param_init(drizdata);
	args->user = drizdata;

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
	seq.nb_layers = args->seq->nb_layers;
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
		// check the reference image matrix is not null
	transformation_type checkH = guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H);
	if (checkH == NULL_TRANSFORMATION) {
		siril_log_color_message(_("The reference image has a null matrix and was not previously aligned, choose another one, aborting\n"), "red");
		return FALSE;
	}
	// check the number of dof if -interp=none
	transformation_type min, max;
	guess_transform_from_seq(regargs->seq, regargs->layer, &min, &max, TRUE);
	if (max > SHIFT_TRANSFORMATION && regargs->interpolation == OPENCV_NONE) {
		siril_log_color_message(_("Applying registration computed with higher degree of freedom (%d) than shift is not allowed when interpolation is set to none, aborting\n"), "red", ((int)max + 1) * 2);
		return FALSE;
	}

	// check the consistency of output images size if -interp=none
	if (regargs->interpolation == OPENCV_NONE && (regargs->x2upscale || regargs->framing == FRAMING_MAX || regargs->framing == FRAMING_MIN)) {
		siril_log_color_message(_("Applying registration with changes in output image sizeis not allowed when interpolation is set to none , aborting\n"), "red");
		return FALSE;
	}

	// check the consistency of images size if -interp=none
	if (regargs->interpolation == OPENCV_NONE && regargs->seq->is_variable) {
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

	// cog frmaing method requires all images to be of same size
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

