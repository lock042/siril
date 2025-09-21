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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/siril_wcs.h"
#include "algos/star_finder.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "drizzle/cdrizzleutil.h"
#include "registration/registration.h"
#include "registration/distorsion.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/match.h"
#include "registration/matching/misc.h"
#include "opencv/opencv.h"

# define MIN_RATIO_INLIERS 30 //percentage of inliers after transformation fitting (shift, affine or homography)
# define AMPLITUDE_CUT 0.05 // percentile to clip the lower tail of amplitudes distribtion when filtering out stars for 2pass reg
# define MAX_SHIFT_RATIO 0.25f // max ratio of image rx for ref image offset from sequence cog
# define MAX_TRIALS_2PASS 5 // max number of trialsto find the best ref

static void print_alignment_results(Homography H, int filenum, float FWHMx, float FWHMy, char *units);

static int get_min_requires_stars(transformation_type type) {
	switch(type) {
	case SHIFT_TRANSFORMATION:
	case SIMILARITY_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		return 3;
	default:
	case HOMOGRAPHY_TRANSFORMATION:
		return 4;
	}
}

const char *describe_transformation_type(transformation_type type) {
	switch(type) {
		case SHIFT_TRANSFORMATION:
			return _("shift-only (2 degrees)");
		case SIMILARITY_TRANSFORMATION:
			return _("similarity (4 degrees)");
		case AFFINE_TRANSFORMATION:
			return _("affine (6 degrees)");
		case HOMOGRAPHY_TRANSFORMATION:
			return _("homography (8 degrees)");
		default:
			return _("INVALID");
	}
}

int registration_prepare_results(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	if (!regargs->no_output) {
		// allocate destination sequence data
		regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
		regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
		if (!regargs->imgparam  || !regargs->regparam) {
			PRINT_ALLOC_ERR;
			return 1;
		}

		if (seq_prepare_hook(args))
			return 1;
	}

	sadata->success = calloc(args->nb_filtered_images, sizeof(BYTE));
	if (!sadata->success) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}

int star_align_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	float FWHMx, FWHMy, B;
	char *units;
	fits fit = { 0 };
	int nb_stars = 0;

	sadata->current_regdata = registration_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;

	/* first we're looking for stars in reference image */
	if (seq_read_frame(args->seq, regargs->reference_image, &fit, FALSE, -1)) {
		siril_log_color_message(_("Could not load reference image\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		return 1;
	}

	siril_log_color_message(_("Reference Image:\n"), "green");

	// peaking or using cached list in lst(if no selection is made)
	// For now selection is using com.selection
	struct starfinder_data *sf_data = findstar_image_worker(regargs->sfargs, -1, regargs->reference_image, &fit, NULL, com.max_thread);
	if (sf_data) {
		sadata->refstars = *sf_data->stars;
		nb_stars = *sf_data->nb_stars;
		free(sf_data);
	}

	siril_log_message(_("Found %d stars in reference, channel #%d\n"), nb_stars, regargs->layer);

	if (!sadata->refstars || nb_stars < get_min_requires_stars(regargs->type)) {
		siril_log_color_message(
				_("There are not enough stars in reference image to perform alignment\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}

	/* copying refstars to com.stars for display
	/ * we make this cop before distortion correction */
	if (!com.script && &com.seq == args->seq && com.seq.current == regargs->reference_image) {
		psf_star **stars = new_fitted_stars(MAX_STARS);
		if (stars) {
			for (int i = 0; i < nb_stars && sadata->refstars[i]; i++) {
				psf_star *tmp = new_psf_star();
				if (!tmp) {
					PRINT_ALLOC_ERR;
					stars[i] = NULL;
					break;
				}
				memcpy(tmp, sadata->refstars[i], sizeof(psf_star));
				stars[i] = tmp;
				stars[i+1] = NULL;
			}
		}
		update_star_list(stars, FALSE, FALSE);
	}

	// We prepare the distortion structure maps if required
	if (regargs->undistort && init_disto_map(fit.rx, fit.ry, regargs->disto)) {
		siril_log_color_message(
				_("Could not init distortion mapping\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}
	
	if (!regargs->no_output && regargs->driz && initialize_drizzle_params(args, regargs)) {
		siril_log_color_message(
				_("Could not init drizzle\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}

	// we correct the reference stars for distortions
	// we'll do the same for each imagebefore matching the 2 lists
	// The star lists are not saved with distortions as they may be re-used
	// later on with a different distortion setting
	// Instead, when registration has been computed accounting for SIP,
	// the distortion source is logged together with the registration info in the .seq file
	// Note: for global reg, disto source can only by DISTO_IMAGE or DISTO_FILE
	// i.e regargs->disto is of size 1
	if (regargs->undistort && disto_correct_stars(sadata->refstars, regargs->disto)) {
		siril_log_color_message(
			_("Could not correct the stars position with SIP coeffients\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}

	sadata->ref.x = fit.rx;
	sadata->ref.y = fit.ry;

	// For internal sequences the data / fdata pointer still
	// points to the original memory in seq->internal_fits.
	// It must not be freed by clearfits here so we set the
	// pointers in fit to NULL
	if (args->seq->type == SEQ_INTERNAL) {
		fit.data = NULL;
		fit.fdata = NULL;
	}
	clearfits(&fit);

	if (!regargs->no_output) {
		sadata->ref.x *= regargs->output_scale;
		sadata->ref.y *= regargs->output_scale;
		regargs->framing = FRAMING_CURRENT;
		cvGetEye(&regargs->framingd.Htransf);
		regargs->framingd.roi_out.w = sadata->ref.x;
		regargs->framingd.roi_out.h = sadata->ref.y;
	}

	if (nb_stars >= MAX_STARS_FITTED) {
		sadata->fitted_stars = MAX_STARS_FITTED;
		siril_log_color_message(_("Reference Image: Limiting to %d brightest stars\n"), "green", MAX_STARS_FITTED);
	} else {
		sadata->fitted_stars = nb_stars;
	}
	FWHM_stats(sadata->refstars, sadata->fitted_stars, args->seq->bitpix, &FWHMx, &FWHMy, &units, &B, NULL, 0.);
	siril_log_message(_("FWHMx:%*.2f %s\n"), 12, FWHMx, units);
	siril_log_message(_("FWHMy:%*.2f %s\n"), 12, FWHMy, units);
	sadata->current_regdata[regargs->reference_image].roundness = FWHMy/FWHMx;
	sadata->current_regdata[regargs->reference_image].fwhm = FWHMx;
	sadata->current_regdata[regargs->reference_image].weighted_fwhm = FWHMx;
	sadata->current_regdata[regargs->reference_image].background_lvl = B;
	sadata->current_regdata[regargs->reference_image].number_of_stars = sadata->fitted_stars;

	return registration_prepare_results(args);
}

static int star_match_and_checks(psf_star **ref_stars, psf_star **stars, int nb_ref_stars, int nb_stars, struct registration_args *regargs, int filenum, Homography *H) {
	double scale_min = 0.9;
	double scale_max = 1.1;
	int attempt = 1;
	int nobj = 0;
	int failure = 1;
	/* make a loop with different tries in order to align the two sets of data */
	while (failure && attempt < NB_OF_MATCHING_TRY) {
		failure = new_star_match(stars, ref_stars, nb_stars, nb_ref_stars, nobj,
				scale_min, scale_max, H, NULL, FALSE, regargs->type, AT_TRANS_UNDEFINED,
				NULL, NULL);
		if (attempt == 1) {
			scale_min = -1.0;
			scale_max = -1.0;
		} else {
			nobj += 50;
		}
		attempt++;
	}
	if (failure) {
		siril_log_color_message(_("Cannot perform star matching: try #%d. Image %d skipped\n"),
				"red", attempt, filenum);
	}
	else if (H->Inliers < regargs->min_pairs) {
		siril_log_color_message(_("Not enough star pairs (%d): Image %d skipped\n"),
				"red", H->Inliers, filenum);
		failure = 1;
	}
	else if (((double)H->Inliers / (double)H->pair_matched) < ((double)MIN_RATIO_INLIERS / 100.)) {
		switch (regargs->type) {
			case SHIFT_TRANSFORMATION:
				siril_log_color_message(_("Less than %d%% star pairs kept by shift model, it may be too rigid for your data: Image %d skipped\n"),
						"red", MIN_RATIO_INLIERS, filenum);
				failure = 1;
				break;
			case SIMILARITY_TRANSFORMATION:
				siril_log_color_message(_("Less than %d%% star pairs kept by similarity model, it may be too rigid for your data: Image %d skipped\n"),
						"red", MIN_RATIO_INLIERS, filenum);
				failure = 1;
				break;
			case AFFINE_TRANSFORMATION:
				siril_log_color_message(_("Less than %d%% star pairs kept by affine model, it may be too rigid for your data: Image %d skipped\n"),
						"red", MIN_RATIO_INLIERS, filenum);
				failure = 1;
				break;
			case HOMOGRAPHY_TRANSFORMATION:
				siril_log_color_message(_("Less than %d%% star pairs kept by homography model, Image %d may show important distortion\n"),
						"salmon", MIN_RATIO_INLIERS, filenum);
				break;
			default:
				printf("Should not happen\n");
		}
	}
	return failure;
}

/* reads the image, searches for stars in it, tries to match them with
 * reference stars, computes the homography matrix, applies it on the image,
 * possibly up-scales the image and stores registration data */
int star_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int nb_stars = 0;
	float FWHMx, FWHMy, B;
	char *units;
	Homography H = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes
	siril_debug_print("registration of image %d using %d threads\n", in_index, threads);

	if (regargs->no_output) {  // used by livestacking but not by register which always outputs sequence
		args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
	}

	if (in_index != regargs->reference_image) {
		psf_star **stars = NULL;
		if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
			siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
		}

		/* sometimes sequence are not consistent.... They shouldn't but ..... */
		// This now normally checked in processing... but ok, livestacking may use this check
		// Wouldn't it be safer to compare fit->naxes[2] to regargs->seq->nb_layers?
		if (regargs->layer > RLAYER && !isrgb(fit)) {
			siril_log_color_message(_("It looks like your sequence contains a mix of monochrome and RGB images.\n"), "salmon");
		}

		// peaking or using cached list in lst(if no selection is made)
		// For now selection is using com.selection
		struct starfinder_data *sf_data = findstar_image_worker(regargs->sfargs, -1, in_index, fit, NULL, com.max_thread);
		if (sf_data) {
			stars = *sf_data->stars;
			nb_stars = *sf_data->nb_stars;
			free(sf_data);
		}

		if (!stars || nb_stars < get_min_requires_stars(regargs->type)) {
			siril_log_message(
					_("Not enough stars. Image %d skipped\n"), filenum);
			if (stars)
				free_fitted_stars(stars);
			args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}
		// if distortion correction is included, we correct the current starlist
		// before performing the match
		if (regargs->undistort && disto_correct_stars(stars, regargs->disto)) {
			siril_log_color_message(
				_("Could not correct the stars position with SIP coeffients\n"), "red");
			if (stars)
				free_fitted_stars(stars);
			args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}

		int not_matched = star_match_and_checks(sadata->refstars, stars, sadata->fitted_stars, nb_stars, regargs, filenum, &H);
		if (!not_matched)
			FWHM_stats(stars, nb_stars, args->seq->bitpix, &FWHMx, &FWHMy, &units, &B, NULL, 0.);
		free_fitted_stars(stars);
		if (not_matched) {
			args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}

#ifdef _OPENMP
#pragma omp critical
#endif
		print_alignment_results(H, filenum, FWHMx, FWHMy/FWHMx, "px");

		sadata->current_regdata[in_index].roundness = FWHMy/FWHMx;
		sadata->current_regdata[in_index].fwhm = FWHMx;
		sadata->current_regdata[in_index].weighted_fwhm = 2 * FWHMx
				* (((double) sadata->fitted_stars) - (double) H.Inliers)
				/ (double) sadata->fitted_stars + FWHMx;
		sadata->current_regdata[in_index].background_lvl = B;
		sadata->current_regdata[in_index].number_of_stars = nb_stars;
		sadata->current_regdata[in_index].H = H;

		if (!regargs->no_output) {
			if (regargs->interpolation <= OPENCV_LANCZOS4) {
				if (apply_reg_image_hook(args, out_index, in_index, fit, _, threads)) {
					args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
					return 1;
				}
			} else {
				if (shift_fit_from_reg(fit, H)) {
					args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
					return 1;
				}
			}
		}
	} else {
		// reference image
		cvGetEye(&H);
		sadata->current_regdata[in_index].H = H;
		if (!regargs->no_output && (regargs->output_scale != 1.f || regargs->driz || regargs->undistort)) {
			if (apply_reg_image_hook(args, out_index, in_index, fit, _, threads)) {
				args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
				return 1;
			}
		}
	}

	// updating regparam for the output sequence if any and updating pixel size 
	// is handled in apply_reg_image_hook
	args->seq->imgparam[in_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	return 0;
}

static int star_align_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	free_fitted_stars(sadata->refstars);

	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;
		if (regargs->new_total <= 1) {
			siril_log_color_message(_("No image was registered to the reference\n"), "red");
			args->retval = 1;
		}
	}

	if (!args->retval) {
		regargs->seq->distoparam[regargs->layer] = regargs->distoparam;

		if (!regargs->no_output) {
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
		}
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
		} else if (args->seq->type == SEQ_REGULAR) {
			remove_prefixed_sequence_files(regargs->seq, regargs->prefix);
			remove_prefixed_drizzle_files(regargs->seq, regargs->prefix);
		}
	}

	if (!args->retval) {
		siril_log_message(_("Registration finished.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->new_total);

		g_free(str);
		if (!regargs->no_output && (args->seq->type != SEQ_INTERNAL)) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_registration(regargs, -1);
			// will be loaded in the idle function if (load_new_sequence)
			regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	if (sadata->success)
		free(sadata->success);
	free(sadata);
	args->user = NULL;
	clear_stars_list(FALSE);
	return regargs->new_total == 0;
}

int star_align_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
		/* The registration memory consumption, n is image size and m channel size.
		 * First, a threshold is computed for star pixel value, using statistics:
		 *	O(m), data is duplicated for median computation if
		 *	there are nil values, O(1) otherwise
		 * Then, still in peaker(), image is filtered using Gaussian blur, duplicating
		 * the reference channel to act as input and output of the filter as float O(m
		 * as float for RT, current), O(2m as float for openCV).
		 * All this is in addition to the image being already loaded.
		 * In the special case of CFA image as input, we also have a copy of the
		 * orig image so O(n=m) to interpolate non green pixels

		 * Then, we use the same function as apply_reg_compute_mem_limits is used to compute
		 * consumption for producing the registered images (if not 2-pass)

		 * Since these three operations are in sequence, we need room only for the
		 * largest.
		 *
		 * Step 1 - statistics:
		 * 1 original image O(n)
		 * 1 copy of original image O(n) if CFA
		 * 1 original image layer O(m)

		 * Step 2 - detection (always larger than step 1):
		 * 1 original image O(n)
		 * 1 copy of original image O(n) if CFA
		 * 1 float original image layer O(m as float)

		 * Step 3 - transformation:
		 * See apply_reg_compute_mem_consumption()
		 * Note: Because scale can be between 0.1 and 3, we can't know how this
		 * compares to step 2 before actually doing the calc

// TODO: do we keep this?
		 * rotated color scaled float	mem needed
		 *       0     0      0     0	O(m as float)
		 *       1     0      0     0	O(m as float)
		 *       1     0      0     1	O(m as float, same as n)
		 *       1     0      1     0	O(4m, same as 2m as float)
		 *       1     0      1     1	O(4m)
		 *       1     1      0     0	O(n)
		 *       1     1      0     1	O(n)
		 *       1     1      1     0	O(8n or 2nscaled)
		 *       1     1      1     1	O(8n or 2nscaled)
		 */

	struct registration_args *regargs = ((struct star_align_data *)args->user)->regargs;

	unsigned int MB_per_orig_image, MB_per_orig_channel_float, MB_per_scaled_image, MB_avail, required_step2 = 0, required_step3 = 0, required = 0;
	uint64_t memory_per_orig_channel_float, memory_per_orig_image;
	int limit_step2 = INT16_MAX, limit_step3 = INT16_MAX, limit;
	gboolean is_float = (get_data_type(args->seq->bitpix) == DATA_FLOAT || args->force_float);

	memory_per_orig_channel_float = (uint64_t) args->seq->rx * args->seq->ry * sizeof(float);
	memory_per_orig_image = memory_per_orig_channel_float * args->seq->nb_layers * ((is_float) ? 1 : 0.5);

	MB_per_orig_image = max(1, memory_per_orig_image / BYTES_IN_A_MB); // i.e n
	MB_per_orig_channel_float = max(1, memory_per_orig_channel_float / BYTES_IN_A_MB); // i.e m as float

	required_step2 = MB_per_orig_image + MB_per_orig_channel_float;

	if (regargs->driz && regargs->driz->is_bayer)
		required_step2 += MB_per_orig_image; // the copy for interpolating nongreen pixels

	limit_step3 = apply_reg_compute_mem_consumption(args, &required_step2, &MB_per_scaled_image, &MB_avail);
	limit_step2 = (int)(MB_avail / required_step2);

	limit = min(limit_step2, limit_step3);
	required = max(required_step2, required_step3);

	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (!regargs->no_output && for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else
			limit = thread_limit;
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

int register_star_alignment(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	args->force_float = !com.pref.force_16bit && regargs->seq->type != SEQ_SER;
	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}

	args->compute_mem_limits_hook = star_align_compute_mem_limits;
	args->upscale_ratio = regargs->output_scale;

	if (regargs->undistort) {
		int status = 1;
		regargs->disto = init_disto_data(&regargs->distoparam, regargs->seq, NULL, regargs->driz != NULL, &status);
		if (status) {
			free(args);
			siril_log_color_message(_("Could not initialize distortion data, aborting\n"), "red");
			return -1;
		}
		if (!regargs->disto) {
			regargs->undistort = DISTO_UNDEF;
		}
	}

	// preparing detection params
	regargs->sfargs = calloc(1, sizeof(struct starfinder_data));
	regargs->sfargs->im.from_seq = regargs->seq;
	regargs->sfargs->layer = regargs->layer;
	regargs->sfargs->keep_stars = TRUE;
	regargs->sfargs->save_to_file = !regargs->matchSelection && !regargs->no_starlist;
	regargs->sfargs->max_stars_fitted = regargs->max_stars_candidates;
	if (regargs->matchSelection && com.selection.w > 0 && com.selection.h > 0)
		regargs->sfargs->selection = com.selection;
	else {
		regargs->sfargs->selection = (rectangle) {0, 0, 0, 0};
		regargs->matchSelection = FALSE;
	}

	args->compute_size_hook = compute_registration_size_hook;
	args->prepare_hook = star_align_prepare_hook;
	args->image_hook = star_align_image_hook;
	args->finalize_hook = star_align_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Global star registration");
	args->has_output = !regargs->no_output;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(regargs->prefix);
	args->load_new_sequence = !regargs->no_output;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free_generic_seq_args(args, FALSE);
		return -1;
	}
	sadata->regargs = regargs;
	args->user = sadata;

	if (!regargs->no_output) {
		regargs->wcsref = get_wcs_ref(regargs->seq);
		if (regargs->wcsref && regargs->undistort && regargs->wcsref->lin.dispre) {
			remove_dis_from_wcs(regargs->wcsref); // we remove distortions as the output is undistorted
		}
	}

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free_generic_seq_args(args, FALSE);
	return regargs->retval;
}

static void print_alignment_results(Homography H, int filenum, float fwhm, float roundness, char *units) {
	double rotation, scale, scaleX, scaleY;
	point shift;
	double inliers;

	/* Matching information */
	siril_log_color_message(_("Matching stars in image %d: done\n"), "green", filenum);
	siril_log_message(_("Initial pair matches: %d\n"), H.pair_matched);
	siril_log_message(_("Pair matches after fitting: %d\n"), H.Inliers);
	inliers = 1.0 - ((((double) H.pair_matched - (double) H.Inliers)) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 11, inliers);

	/* Scale */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	scale = (scaleX + scaleY) * 0.5;
	siril_log_message(_("scaleX:%*.3f\n"), 12, scaleX);
	siril_log_message(_("scaleY:%*.3f\n"), 12, scaleY);
	siril_log_message(_("scale:%*.3f\n"), 13, scale);

	/* Rotation */
	rotation = atan2(H.h01, H.h00) * 180 / M_PI;
	siril_log_message(_("rotation:%+*.3f deg\n"), 10, rotation);

	/* Translation */
	shift.x = -H.h02;
	shift.y = -H.h12;
	siril_log_message(_("dx:%+*.2f px\n"), 15, shift.x);
	siril_log_message(_("dy:%+*.2f px\n"), 15, shift.y);
	siril_log_message(_("FWHM:%*.2f %s\n"), 13, fwhm, units);
	siril_log_message(_("roundness:%*.2f\n"), 8, roundness);
}

static int compute_transform(struct registration_args *regargs, struct starfinder_data *sfargs, gboolean *included, int *failed, const float *fwhm, const float *roundness, const float *B, gboolean verbose) {
	regdata *current_regdata = registration_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata) return -1;
	int nb_ref_stars = sfargs->nb_stars[regargs->seq->reference_image];
	int nb_aligned = 0;
	int nbfail = *failed;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) shared(nbfail, nb_aligned)
#endif
	for (int i = 0; i < regargs->seq->number; i++) {
		if (!included[i])
			continue;
		Homography H = { 0 };
		if (i == regargs->seq->reference_image) {
			cvGetEye(&H);
		} else {
			int filenum = regargs->seq->imgparam[i].filenum;	// for display purposes
			int not_matched = star_match_and_checks(sfargs->stars[regargs->seq->reference_image], sfargs->stars[i],
					sfargs->nb_stars[regargs->seq->reference_image], sfargs->nb_stars[i], regargs, filenum, &H);
			if (not_matched) {
				g_atomic_int_inc(&nbfail);
				included[i] = FALSE;
				continue;
			}
#ifdef _OPENMP
#pragma omp critical
#endif
			if (verbose) print_alignment_results(H, filenum, fwhm[i], roundness[i], "px");

		}
		g_atomic_int_inc(&nb_aligned);
		current_regdata[i].roundness = roundness[i];
		current_regdata[i].fwhm = fwhm[i];
		current_regdata[i].weighted_fwhm = 2 * fwhm[i]
			* ((double)nb_ref_stars - sfargs->nb_stars[i])
			/ (double)nb_ref_stars + fwhm[i];
		current_regdata[i].background_lvl = B[i];
		current_regdata[i].number_of_stars = sfargs->nb_stars[i];
		current_regdata[i].H = H;
	}
	*failed = nbfail;
	return nb_aligned;
}

static void compute_dist(struct registration_args *regargs, float *dist, const gboolean *included) {
	Homography Href = regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H;
	Homography Hshift = {0};
	Homography Htransf = {0};
	cvGetEye(&Hshift);
	int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	int x0, y0, n;
	double cogx, cogy;

	point center = { 0 };
	center.x =  (double)rx * 0.5;
	center.y =  (double)ry * 0.5;

	cogx = 0.;
	cogy = 0;
	n = 0;
	for (int i = 0; i < regargs->seq->number; i++) {
		if (!included[i]) continue;
		point currcenter;
		memcpy(&currcenter, &center, sizeof(point));
		cvTransfPoint(&currcenter.x, &currcenter.y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
		cogx += currcenter.x;
		cogy += currcenter.y;
		n++;
	}
	cogx /= (n == 0 ? 1. : (double)n);
	cogy /= (n == 0 ? 1. : (double)n);
	x0 = (int)(cogx - (double)rx * 0.5);
	y0 = (int)(cogy - (double)ry * 0.5);
	Hshift.h02 = (double)x0;
	Hshift.h12 = (double)y0;
	siril_debug_print("cog at (%3.2f, %3.2f)\n", Hshift.h02, Hshift.h12);
	cvMultH(Href, Hshift, &Htransf);
	for (int i = 0; i < regargs->seq->number; i++) {
		if (!included[i]) continue;
		point currcenter;
		double dx, dy;
		memcpy(&currcenter, &center, sizeof(point));
		cvTransfPoint(&currcenter.x, &currcenter.y, regargs->seq->regparam[regargs->layer][i].H, Htransf, 1.);
		dx = currcenter.x - center.x;
		dy = currcenter.y - center.y;
		dist[i] = (float)sqrt(dx * dx + dy * dy);
	}

}

// returns the index of the minimum element of float array arr of size nb
// if gboolean mask array is passed, it only includes elements where mask is TRUE
// if *vval is passed, it contains the best value
int minidx(const float *arr, const gboolean *mask, int nb, float *val) {
	if (!arr) return -1;
	if (nb < 1) return -1;
	int idx = -1;
	float best = FLT_MAX;
	gboolean check_included = TRUE;
	if (!mask) check_included = FALSE;
	for (int i = 0; i < nb; i++) {
		if ((!check_included || mask[i]) && (arr[i] < best)) {
			best = arr[i];
			idx = i;
		}
	}
	if (val) *val = best;
	return idx;
}

/********************** the new multi-step registration ***********************/
/* This is a modification of the global registration without output, separating
 * the two main steps of finding stars and computing the transforms to add an
 * extra step in the middle:
* 1. finding stars
* 2. searching for the best image and set it as reference
* 3. compute the transforms and store them in regparams
*/

int register_multi_step_global(struct registration_args *regargs) {
	// 1. finding stars
	int retval = 0;
	float *fwhm = NULL, *roundness = NULL, *A = NULL, *B = NULL, *Acut = NULL, *scores = NULL;
	float *dist = NULL;
	// local flag (and its copy) accounting both for process_all_frames flag and collecting failures along the process
	gboolean *included = NULL, *tmp_included = NULL;
	// local flag to make checks only on frames that matter
	gboolean *meaningful = NULL;

	if (regargs->undistort) {
		int status = 1;
		regargs->disto = init_disto_data(&regargs->distoparam, regargs->seq, NULL, regargs->driz != NULL, &status);
		if (status) {
			siril_log_color_message(_("Could not initialize distortion data, aborting\n"), "red");
		}
		if (!regargs->disto) {
			regargs->undistort = DISTO_UNDEF;
		}
	}

	struct starfinder_data *sfargs = calloc(1, sizeof(struct starfinder_data));
	sfargs->im.from_seq = regargs->seq;
	sfargs->layer = regargs->layer;
	sfargs->max_stars_fitted = regargs->max_stars_candidates;
	sfargs->stars = calloc(regargs->seq->number, sizeof(psf_star **));
	if (!sfargs->stars) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	sfargs->nb_stars = calloc(regargs->seq->number, sizeof(int));
	if (!sfargs->nb_stars) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	sfargs->update_GUI = FALSE;
	sfargs->already_in_thread = TRUE;
	sfargs->process_all_images = !regargs->filters.filter_included;
	sfargs->save_to_file = !regargs->no_starlist;
	sfargs->save_eqcoords = FALSE;
	struct timeval t_start, t_end;

	gettimeofday(&t_start, NULL);
	if (apply_findstar_to_sequence(sfargs)) {
		siril_debug_print("finding stars failed\n");	// aborted probably
		retval = 1;
		goto free_all;
	}

	// 2. searching for the best image and set it as reference
	fwhm = calloc(regargs->seq->number, sizeof(float));
	roundness = calloc(regargs->seq->number, sizeof(float));
	A = calloc(regargs->seq->number, sizeof(float));
	B = calloc(regargs->seq->number, sizeof(float));
	Acut = calloc(regargs->seq->number, sizeof(float));
	included = calloc(regargs->seq->number, sizeof(gboolean));
	tmp_included = calloc(regargs->seq->number, sizeof(gboolean));
	meaningful = calloc(regargs->seq->number, sizeof(gboolean));
	scores = calloc(regargs->seq->number, sizeof(float));
	dist = calloc(regargs->seq->number, sizeof(float));
	if (!fwhm || !roundness || !B || !A || !included || !tmp_included || !meaningful || !scores || !dist) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	int maxstars = 0;
	int failed = 0;		// number of images to fail registration at any step
	int best_index = -1;
	for (int i = 0; i < regargs->seq->number; i++) {
		if (regargs->filters.filter_included && !regargs->seq->imgparam[i].incl)
			continue;
		if (!sfargs->stars[i]) {
			// star finder failed, we exclude the frame from the sequence
			// need this check to be after filters.filter_included otherwise, we would report false failures
			regargs->seq->imgparam[i].incl = FALSE;
			failed++;
			continue;
		}
		// we apply distortion (if any) before matching
		if (regargs->undistort && disto_correct_stars(sfargs->stars[i], regargs->disto)) {
			siril_log_color_message(_("Could not correct the stars position with SIP coeffients\n"), "red");
			retval = 1;
			goto free_all;
		}
		included[i] = TRUE;
		float FWHMx, FWHMy;
		char *units;
		FWHM_stats(sfargs->stars[i], sfargs->nb_stars[i], regargs->seq->bitpix, &FWHMx, &FWHMy, &units, B + i, Acut + i, AMPLITUDE_CUT);
		fwhm[i] = FWHMx;
		roundness[i] = FWHMy/FWHMx;
		if (sfargs->nb_stars[i] > maxstars) maxstars = sfargs->nb_stars[i];
	}


	if (maxstars == sfargs->max_stars_fitted) {
		siril_log_message(_("The number of stars has capped, readapting threshold and filtering\n"));
		float Athreshold = 0.0f;
		// Determine the updated A threshold
		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			if (sfargs->nb_stars[i] == sfargs->max_stars_fitted && Acut[i] > Athreshold) Athreshold = Acut[i];
		}
		// Filter the whole series against new amplitude threshold
		maxstars = maxstars * (1 - AMPLITUDE_CUT);

		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			sfargs->stars[i] = filter_stars_by_amplitude(sfargs->stars[i], Athreshold, &sfargs->nb_stars[i]);
			if (!sfargs->stars[i]) {
				regargs->seq->imgparam[i].incl = FALSE;
				included[i] = FALSE;
				failed++;
				continue;
			}
			float FWHMx, FWHMy;
			char *units;
			FWHM_stats(sfargs->stars[i], sfargs->nb_stars[i], regargs->seq->bitpix, &FWHMx, &FWHMy, &units, B + i, NULL, 0.);
			// we only rank images with at least half the maximum number of stars (and we count them as meaningful)
			scores[i]  = (sfargs->nb_stars[i] >= maxstars / 2) ? 2. * FWHMx * (maxstars - sfargs->nb_stars[i]) / (maxstars == 0 ? 1 : maxstars) + FWHMx : FLT_MAX;
			if (sfargs->nb_stars[i] >= maxstars / 2) meaningful[i] = TRUE;
			fwhm[i] = FWHMx;
			roundness[i] = FWHMy/FWHMx;
		}
		best_index = minidx(scores, included, regargs->seq->number, NULL);
	} else {
		float FWHMx;
		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			FWHMx = fwhm[i];
			// we only rank images with at least half the maximum number of stars (and we count them as meaningful)
			scores[i]  = (sfargs->nb_stars[i] >= maxstars / 2) ? 2. * FWHMx * (maxstars - sfargs->nb_stars[i]) / (maxstars == 0 ? 1 : maxstars) + FWHMx : FLT_MAX;
			if (sfargs->nb_stars[i] >= maxstars / 2) meaningful[i] = TRUE;
		}
		best_index = minidx(scores, included, regargs->seq->number, NULL);
	}

	regargs->seq->reference_image = best_index;
	int reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
	siril_log_message(_("Trial #%d: After sequence analysis, we are choosing image %d as new reference for registration\n"), 1, reffilenum);

	// 3. compute the transforms and store them in regparams
	// we will check that we have registered enough meaningful frames to the new ref before going to distance checking (step 4.)

	//intializing some useful info
	int trials = 0;
	int nb_meaningful = 0;
	for (int i = 0; i < regargs->seq->number; i++) {
		if (meaningful[i]) nb_meaningful++;
	}
	int best_indexes[MAX_TRIALS_2PASS];
	int nb_aligned[MAX_TRIALS_2PASS];
	best_indexes[trials] = best_index;
	float allowable_dist = (float)regargs->seq->imgparam[regargs->reference_image].rx * MAX_SHIFT_RATIO;
	int tmp_failed = 0;

	int max_trials = min(MAX_TRIALS_2PASS, regargs->seq->number);
	while (trials < max_trials) {
		tmp_failed = failed;
		for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
		nb_aligned[trials] = compute_transform(regargs, sfargs, tmp_included, &tmp_failed, fwhm, roundness, B, FALSE);
		// if number of aligned frames is less than half the number of meaningful frames (those with enough stars)
		// we have chosen a reference which is not framed well enough to align the sequence (the computed cog is probably meaningless as well)
		// we set its score to FLT_MAX and start again with the next best frame
		// we do not remove it from the meaningful frames count or included ones, as it still may be worth aligning with the new ref
		if (nb_aligned[trials] < nb_meaningful / 2) {
			scores[best_index] = FLT_MAX;
			siril_log_message(_("Trial #%d: After sequence alignment, image #%d could not align more than half of the frames, recomputing\n"), trials + 1, reffilenum);
			best_index = minidx(scores, included, regargs->seq->number, NULL);
			regargs->seq->reference_image = best_index;
			reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
			trials++;
			if (trials < max_trials) {
				siril_log_message(_("Trial #%d: After sequence analysis, we are choosing image %d as new reference for registration\n"), trials + 1, reffilenum);
				best_indexes[trials] = best_index;
			}
		} else { // not necessary but a simple to have print_alignment_results
			tmp_failed = failed;
			for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
			compute_transform(regargs, sfargs, tmp_included, &tmp_failed, fwhm, roundness, B, TRUE);
			break;
		}
	}

	if (trials == max_trials) {
	// We have tried many times over, we need to select the "best" frame
	// even though it cannot align more than half of the good images
		int best_try = -1;
		int maxreg = 0;
		siril_log_message(_("After %d trials, no reference image could align more than half of the frames, selecting best candidate\n"), max_trials);
		for (int i = 0; i < max_trials; i++) {
			siril_log_message(_("Trial #%d: Reference image: #%d - Frames aligned: %d\n"), i + 1, regargs->seq->imgparam[best_indexes[i]].filenum, nb_aligned[i]);
			if (nb_aligned[i] > maxreg) {
				best_try = i;
				maxreg = nb_aligned[i];
			}
		}
		if (maxreg <= 1) {
			siril_log_color_message(_("Could not find an image that aligns more than itself, aborting\n"), "red");
			retval = 1;
			goto free_all;
		}
		best_index = best_indexes[best_try];
		regargs->seq->reference_image = best_index;
		reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
		siril_log_message(_("After sequence analysis, we are choosing image %d as new reference for registration\n"), reffilenum);
		tmp_failed = failed;
		for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
		compute_transform(regargs, sfargs, tmp_included, &tmp_failed, fwhm, roundness, B, TRUE);
	}
	// and we copy back to the initial arrays
	for (int i = 0; i < regargs->seq->number; i++) included[i] = tmp_included[i];
	failed = tmp_failed;

	// 4. check distance to cog
	compute_dist(regargs, dist, included);
	//if larger than cog, we recompute a score accoutning for the distance
	if (dist[best_index] > allowable_dist) {
		siril_log_message(_("After sequence alignment, image %d is too far from the sequence cog, recomputing\n"), reffilenum);
		siril_debug_print(_("Distance to cog is  %2.1fpx while threshold is set at %2.1fpx\n"), dist[best_indexes[trials]], allowable_dist);
		float FWHMx;
		int new_best_index = -1;
		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			FWHMx = fwhm[i];
			// images are now scored ONLY if they are within allowable distance from cog
			scores[i]  = (sfargs->nb_stars[i] >= maxstars / 2 && dist[i] < allowable_dist) ? 2. * FWHMx * (maxstars - sfargs->nb_stars[i]) / maxstars + FWHMx: FLT_MAX;
		}
		new_best_index = minidx(scores, included, regargs->seq->number, NULL);
		if (new_best_index != best_index && new_best_index > -1) { // do not recompute if none or same is found (same should not happen)
			regargs->seq->reference_image = new_best_index;
			reffilenum = regargs->seq->imgparam[new_best_index].filenum;	// for display purposes
			siril_log_message(_("After sequence analysis, we are choosing image %d as new reference for registration\n"), reffilenum);
			// back to 3. compute the transforms and store them in regparams
			compute_transform(regargs, sfargs, included, &failed, fwhm, roundness, B, TRUE);
		} else {
			siril_log_message(_("Could not find a better frame, keeping image %d as the reference for the sequence\n"), reffilenum);
		}
	} else {
		siril_log_message(_("After sequence alignment, we find that image %d is %2.1fpx away from sequence cog, i.e. within allowable bounds\n"), reffilenum, dist[best_index]);
	}

	// we finally copy the local included flag to seq->imgparam
	for (int i = 0; i < regargs->seq->number; i++) {
		regargs->seq->imgparam[i].incl = included[i];
	}
	regargs->seq->distoparam[regargs->layer] = regargs->distoparam;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(regargs->seq, FALSE);
	siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->seq->selnum);
	clear_stars_list(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

free_all:
	if (sfargs->stars) {
		for (int i = 0; i < regargs->seq->number; i++)
			free_fitted_stars(sfargs->stars[i]);
	}
	free(sfargs->stars);
	free(sfargs->nb_stars);
	free(sfargs);
	free(fwhm);
	free(roundness);
	free(B);
	free(A);
	free(Acut);
	free(included);
	free(tmp_included);
	free(meaningful);
	free(scores);
	free(dist);
	return retval;
}

