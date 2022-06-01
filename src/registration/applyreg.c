/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"

#include "opencv/opencv.h"


static void create_output_sequence_for_apply_reg(struct registration_args *args);

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
	if (fit.naxes[2] == 1 && fit.bayer_pattern[0] != '\0')
		siril_log_color_message(_("Applying transformation on a sequence opened as CFA is a bad idea.\n"), "red");


	sadata->ref.x = fit.rx;
	sadata->ref.y = fit.ry;

	clearfits(&fit);

	if (regargs->x2upscale) {
			sadata->ref.x *= 2.0;
			sadata->ref.y *= 2.0;
	}
	return apply_reg_prepare_results(args);
}

/* reads the image and apply existing transformation */
int apply_reg_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	Homography H = { 0 };
	Homography Href = { 0 };
	Homography Himg = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (in_index != regargs->reference_image) {
		if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
			siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
		}

		// Composing transformation wrt reference image
		Himg = regargs->seq->regparam[regargs->layer][in_index].H;
		if (guess_transform_from_H(Himg) == -2)  return 1; // in case H is null and -selected was not passed
		Href = regargs->seq->regparam[regargs->layer][regargs->reference_image].H;
		cvTransfH(Himg, Href, &H);

		if (regargs->interpolation <= OPENCV_LANCZOS4) {
			if (cvTransformImage(fit, sadata->ref.x, sadata->ref.y, H, regargs->x2upscale, regargs->interpolation)) {
				return 1;
			}
		} else {
			if (shift_fit_from_reg(fit, regargs, H)) {
				return 1;
			}
		}
	} else {
		// reference image
		if (regargs->x2upscale) {
			if (cvResizeGaussian(fit, fit->rx * 2, fit->ry * 2, OPENCV_NEAREST))
				return 1;
		}
	}

	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	regargs->imgparam[out_index].rx = args->seq->imgparam[in_index].rx;
	regargs->imgparam[out_index].ry = args->seq->imgparam[in_index].ry;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
	regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
	regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
	cvGetEye(&regargs->regparam[out_index].H);

	if (regargs->x2upscale) {
		fit->pixel_size_x /= 2;
		fit->pixel_size_y /= 2;
		regargs->regparam[out_index].fwhm *= 2.0;
		regargs->regparam[out_index].weighted_fwhm *= 2.0;
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
		free(args->seq->regparam[regargs->layer]);
		args->seq->regparam[regargs->layer] = NULL;

		// args->new_ser can be null if stars were not detected in the reference image
		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		}
		else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
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
		
		// explicit sequence creation to copy imgparam and regparam
		create_output_sequence_for_apply_reg(regargs);
		// will be loaded in the idle function if (load_new_sequence)
		regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded

	}
	else {
		siril_log_message(_("Transformation aborted.\n"));
	}
	return regargs->new_total == 0;
	// TODO: args is never freed because we don't call an end function for
	// this generic processing function. The register idle is called for
	// everything else, but does not know this pointer, and we cannot free
	// it here because it's still used in the generic processing function.
}

int apply_reg_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) { 
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	
	/* The transformation memory consumption is:
		* the original image
		* the transformed image, including upscale if required (4x)
		*/
	unsigned int required = MB_per_orig_image + MB_per_scaled_image;
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

	// need to do the checks before prepare_hook as they may alter regargs->process_all_frames
	if (!check_before_applyreg(regargs)) {
		free(args);
		return -1;
	}

	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = apply_reg_compute_mem_limits;
	args->prepare_hook = apply_reg_prepare_hook;
	args->image_hook = apply_reg_image_hook;
	args->finalize_hook = apply_reg_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Apply registration");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
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
	seq.nb_layers = args->seq->nb_layers;
	seq.rx = args->seq->rx;
	seq.ry = args->seq->ry;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number-1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	seq.is_variable = args->seq->is_variable;
	// Copy from old sequence, we want to keep them consistent
	seq.reference_image = args->seq->reference_image;
	seq.needs_saving = TRUE;
	writeseqfile(&seq);
	free_sequence(&seq, FALSE);
}

int guess_transform_from_H(Homography H) {
	if (fabs(H.h00 + H.h01 + H.h02 + H.h10 + H.h11 + H.h12 + H.h20 + H.h21 + H.h22) < __DBL_EPSILON__) return -2; //null matrix
	if (fabs(H.h20) > __DBL_EPSILON__ || fabs(H.h21) > __DBL_EPSILON__) return HOMOGRAPHY_TRANSFORMATION;
	if (fabs(H.h00 - 1.) < __DBL_EPSILON__ && fabs(H.h11 - 1.) < __DBL_EPSILON__ && fabs(H.h10) < __DBL_EPSILON__ && fabs(H.h01) < __DBL_EPSILON__) {
		if (fabs(H.h02) > __DBL_EPSILON__ || fabs(H.h12) > __DBL_EPSILON__) return SHIFT_TRANSFORMATION;
		return -1; //identity matrix
	}
	if (fabs(H.h10 - H.h00  + H.h01 + H.h11) < __DBL_EPSILON__) return SIMILARITY_TRANSFORMATION;
	return AFFINE_TRANSFORMATION;
}

void guess_transform_from_seq(sequence *seq, int layer, int *min, int *max, gboolean excludenull) {
	int val;
	*min = HOMOGRAPHY_TRANSFORMATION;
	*max = -3;
	gboolean needs_sel_update = FALSE;

	if (!layer_has_registration(seq, layer)) {
		siril_debug_print("No registration data found in sequence or layer\n");
		return;
	}
	for (int i = 0; i < seq->number; i++){
		if (!(&seq->regparam[layer][i])) continue;
		val = guess_transform_from_H(seq->regparam[layer][i].H);
		siril_debug_print("Image #%d - transf = %d\n", i+1, val);
		if (*max < val) *max = val;
		if (*min > val) *min = val;
		if ((val == -2) && excludenull) {
			seq->imgparam[i].incl = FALSE;
			needs_sel_update = TRUE;
		}
	}
	if (excludenull && needs_sel_update)
		fix_selnum(seq, FALSE);
}

gboolean check_before_applyreg(struct registration_args *regargs) {
		// check the reference image matrix is not null
	int checkH = guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H);
	if (checkH == -2) {
		siril_log_color_message(_("The reference image has a null matrix and was not previously aligned, choose another one, aborting\n"), "red");
		return FALSE;
	}
	// check the number of dof if -interp=none
	int min, max; 
	guess_transform_from_seq(regargs->seq, regargs->layer, &min, &max, TRUE);
	if (max > SHIFT_TRANSFORMATION && regargs->interpolation == OPENCV_NONE) {
		siril_log_color_message(_("Applying registration computed with higher degree of freedom (%d) than shift is not allowed when interpolation is set to none, aborting\n"), "red", (max + 1) * 2);
		return FALSE;
	}
	// check that we are not trying to apply identity transform to all the images
	if (max == -1) {
		siril_log_color_message(_("Existing registration data is a set of identity matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// check that we are not trying to apply null transform to all the images
	if (max == -2 || (regargs->seq->selnum <= 1) ) {
		siril_log_color_message(_("Existing registration data is a set of null matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// force -selected if some matrices were null
	if (min == -2) {
		siril_log_color_message(_("Some images were not registered, excluding them\n"), "salmon");
		regargs->process_all_frames = FALSE;
	}

	// Remove the files that we are about to create
	remove_prefixed_sequence_files(regargs->seq, regargs->prefix);

	// testing free space
	int nb_frames = regargs->process_all_frames ? regargs->seq->number : regargs->seq->selnum;
	int64_t size = seq_compute_size(regargs->seq, nb_frames, get_data_type(regargs->seq->bitpix));
	if (regargs->x2upscale)
		size *= 4;
	if (test_available_space(size)) {
		siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
		return FALSE;
	}
	return TRUE;
}

