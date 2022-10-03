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
#include "core/siril_log.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/match.h"
#include "registration/matching/misc.h"
#include "opencv/opencv.h"

# define MIN_RATIO_INLIERS 30 //percentage of inliers after transformation fitting (shift, affine or homography)
# define AMPLITUDE_CUT 0.05 // percentile to clip the lower tail of amplitudes distribtion when filtering out stars for 2pass reg
# define MAX_SHIFT_RATIO 0.25f // max ratio of image rx for ref image offset from sequence cog
# define MAX_TRIALS_2PASS 5 // max number of trialsto find the best ref

static void create_output_sequence_for_global_star(struct registration_args *args);
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

regdata *star_align_get_current_regdata(struct registration_args *regargs) {
	regdata *current_regdata;
	if (regargs->seq->regparam[regargs->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = regargs->seq->regparam[regargs->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, regargs->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(regargs->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		regargs->seq->regparam[regargs->layer] = current_regdata;
	}
	return current_regdata;
}

int star_align_prepare_results(struct generic_seq_args *args) {
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
	int i, nb_stars = 0;

	sadata->current_regdata = star_align_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;

	/* first we're looking for stars in reference image */
	if (seq_read_frame(args->seq, regargs->reference_image, &fit, FALSE, -1)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		return 1;
	}
	if (fit.naxes[2] == 1 && fit.bayer_pattern[0] != '\0')
		siril_log_color_message(_("Registering a sequence opened as CFA is a bad idea.\n"), "red");

	siril_log_color_message(_("Reference Image:\n"), "green");
	image refimage = { .fit = &fit, .from_seq = args->seq, .index_in_seq = regargs->reference_image };

	if (regargs->matchSelection && regargs->selection.w > 0 && regargs->selection.h > 0) {
		sadata->refstars = peaker(&refimage, regargs->layer, &com.pref.starfinder_conf, &nb_stars, &regargs->selection, FALSE, TRUE, regargs->max_stars_candidates, com.max_thread);
	}
	else {
		sadata->refstars = peaker(&refimage, regargs->layer, &com.pref.starfinder_conf, &nb_stars, NULL, FALSE, TRUE, regargs->max_stars_candidates, com.max_thread);
	}

	siril_log_message(_("Found %d stars in reference, channel #%d\n"), nb_stars, regargs->layer);


	if (!sadata->refstars || nb_stars < get_min_requires_stars(regargs->type)) {
		siril_log_message(
				_("There are not enough stars in reference image to perform alignment\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}
	if (!com.script && &com.seq == args->seq && com.seq.current == regargs->reference_image)
		queue_redraw(REDRAW_OVERLAY); // draw stars

	sadata->ref.x = fit.rx;
	sadata->ref.y = fit.ry;

	clearfits(&fit);

	if (regargs->x2upscale) {
		if (regargs->no_output) {
			args->seq->upscale_at_stacking = 2.0;
		} else {
			sadata->ref.x *= 2.0;
			sadata->ref.y *= 2.0;
		}
	}
	else {
		if (regargs->no_output) {
			args->seq->upscale_at_stacking = 1.0;
		}
	}

	/* copying refstars to com.stars for display */
	if (sequence_is_loaded()) {
		com.stars = new_fitted_stars(MAX_STARS);
		if (com.stars) {
			i = 0;
			while (i < MAX_STARS && sadata->refstars[i]) {
				psf_star *tmp = new_psf_star();
				if (!tmp) {
					PRINT_ALLOC_ERR;
					com.stars[i] = NULL;
					break;
				}
				memcpy(tmp, sadata->refstars[i], sizeof(psf_star));
				com.stars[i] = tmp;
				com.stars[i+1] = NULL;
				i++;
			}
		}
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

	return star_align_prepare_results(args);
}

static int star_match_and_checks(psf_star **ref_stars, psf_star **stars, int nb_stars, struct registration_args *regargs, int filenum, Homography *H) {
	double scale_min = 0.9;
	double scale_max = 1.1;
	int attempt = 1;
	int nobj = 0;
	int failure = 1;
	/* make a loop with different tries in order to align the two sets of data */
	while (failure && attempt < NB_OF_MATCHING_TRY) {
		failure = new_star_match(stars, ref_stars, nb_stars, nobj,
				scale_min, scale_max, H, FALSE, NULL, NULL, regargs->type,
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
	int nbpoints, nb_stars = 0;
	float FWHMx, FWHMy, B;
	char *units;
	Homography H = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes
	siril_debug_print("registration of image %d using %d threads\n", in_index, threads);

	if (regargs->no_output) {
		/* if "save transformation only", we choose to initialize all frames
		 * to exclude status. If registration is ok, the status is
		 * set to include */
		args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
	}

	if (in_index != regargs->reference_image) {
		psf_star **stars;
		if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
			siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
		}

		/* sometimes sequence are not consistent.... They shouldn't but ..... */
		int layer;
		if (regargs->layer > RLAYER && !isrgb(fit)) {
			layer = RLAYER;
			siril_log_color_message(_("It looks like your sequence contains a mix of monochrome and RGB images.\n"), "salmon");
		} else {
			layer = regargs->layer;
		}

		image im = { .fit = fit, .from_seq = args->seq, .index_in_seq = in_index };
		if (regargs->matchSelection && regargs->selection.w > 0 && regargs->selection.h > 0) {
			stars = peaker(&im, layer, &com.pref.starfinder_conf, &nb_stars, &regargs->selection, FALSE, TRUE, regargs->max_stars_candidates, threads);
		}
		else {
			stars = peaker(&im, layer, &com.pref.starfinder_conf, &nb_stars, NULL, FALSE, TRUE, regargs->max_stars_candidates, threads);
		}

		siril_log_message(_("Found %d stars in image %d, channel #%d\n"), nb_stars, filenum, regargs->layer);

		if (!stars || nb_stars < get_min_requires_stars(regargs->type)) {
			siril_log_message(
					_("Not enough stars. Image %d skipped\n"), filenum);
			if (stars) free_fitted_stars(stars);
			args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}

		if (nb_stars >= sadata->fitted_stars) {
			if (nb_stars >= MAX_STARS_FITTED) {
				siril_log_color_message(_("Target Image: Limiting to %d brightest stars\n"), "green", MAX_STARS_FITTED);
			}
			nbpoints = sadata->fitted_stars;
		}
		else {
			nbpoints = nb_stars;
		}

		int not_matched = star_match_and_checks(sadata->refstars, stars, nbpoints, regargs, filenum, &H);
		if (!not_matched)
			FWHM_stats(stars, nbpoints, args->seq->bitpix, &FWHMx, &FWHMy, &units, &B, NULL, 0.);
		free_fitted_stars(stars);
		if (not_matched) {
			args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		print_alignment_results(H, filenum, FWHMx, FWHMy/FWHMx, units);

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
				gboolean originally_WORD;
				fitlog(fit, &originally_WORD); // Experiment: log transform to try to reduce ringing
				if (cvTransformImage(fit, sadata->ref.x, sadata->ref.y, H, regargs->x2upscale, regargs->interpolation)) {
					args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
					return 1;
				}
				invfitlog(fit, originally_WORD); // Experiment: invert the log transform on completion of LANCZOS4
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
		if (regargs->x2upscale && !regargs->no_output) {
			if (cvResizeGaussian(fit, fit->rx * 2, fit->ry * 2, OPENCV_NEAREST)) {
				args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
				return 1;
			}
		}
	}

	if (!regargs->no_output) {
		regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
		regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
		regargs->imgparam[out_index].rx = sadata->ref.x;
		regargs->imgparam[out_index].ry = sadata->ref.y;
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
	} else {
		// TODO: check if H matrix needs to include a flip or not based on fit->top_down
		// seems like not but this could backfire at some point
		args->seq->imgparam[in_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	}
	sadata->success[out_index] = 1;
	return 0;
}

int star_align_finalize_hook(struct generic_seq_args *args) {
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
		}
	}

	if (sadata->success) free(sadata->success);
	free(sadata);
	args->user = NULL;
	clear_stars_list(FALSE);

	if (!args->retval) {
		siril_log_message(_("Registration finished.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->new_total);

		g_free(str);
		if (!regargs->no_output) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_global_star(regargs);
			// will be loaded in the idle function if (load_new_sequence)
			regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	return regargs->new_total == 0;
	// TODO: args is never freed because we don't call an end function for
	// this generic processing function. The register idle is called for
	// everything else, but does not know this pointer, and we cannot free
	// it here because it's still used in the generic processing function.
}

int star_align_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	unsigned int required = MB_per_scaled_image;
	if (limit > 0) {
		/* The registration memory consumption, n is image size and m channel size.
		 * First, a threshold is computed for star pixel value, using statistics:
		 *	O(m), data is duplicated for median computation if
		 *	there are nil values, O(1) otherwise
		 * Then, still in peaker(), image is filtered using Gaussian blur, duplicating
		 * the reference channel to act as input and output of the filter as float O(m
		 * as float for RT, current), O(2m as float for openCV).
		 * Then, the image is rotated and upscaled by the generic function if enabled:
		 * cvTransformImage is O(n) in mem for unscaled, O(nscaled)=O(4m) for
		 * monochrome scaled and O(2nscaled)=O(21m) for color scaled
		 * All this is in addition to the image being already loaded, except for the
		 * color scaled image.
		 *
		 * Since these three operations are in sequence, we need room only for the
		 * largest.
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
		int is_color = args->seq->nb_layers == 3;
		int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
		int is_scaled = args->upscale_ratio == 2.0;
		unsigned int float_multiplier = is_float ? 1 : 2;
		unsigned int MB_per_float_image = MB_per_orig_image * float_multiplier;
		unsigned int MB_per_float_channel = is_color ? MB_per_float_image / 3 : MB_per_float_image;
		unsigned int MB_per_orig_channel = is_color ? MB_per_orig_image / 3 : MB_per_float_image;
		MB_per_float_channel = min(1, MB_per_float_channel);
		MB_per_orig_channel = min(1, MB_per_orig_channel);
		if (!args->has_output || (!is_scaled && !is_color)) {
			required = MB_per_orig_image + MB_per_float_channel;
		}
		else if (args->has_output && !is_color && is_scaled) {
			required = MB_per_orig_image + 4 * MB_per_orig_channel;
		}
		else if (args->has_output && is_color && !is_scaled) {
			required = 2 * MB_per_orig_image;
		}
		else {
			required = 2 * MB_per_scaled_image;
		}
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

int register_star_alignment(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = star_align_compute_mem_limits;
	args->prepare_hook = star_align_prepare_hook;
	args->image_hook = star_align_image_hook;
	args->finalize_hook = star_align_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Global star registration");
	args->has_output = !regargs->no_output;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	args->new_seq_prefix = regargs->prefix;
	args->load_new_sequence = !regargs->no_output;
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

static void create_output_sequence_for_global_star(struct registration_args *args) {
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
	seq.is_variable = FALSE;
	// don't copy from old sequence, it may not be the same image
	seq.reference_image = sequence_find_refimage(&seq);
	seq.needs_saving = TRUE;
	writeseqfile(&seq);
	free_sequence(&seq, FALSE);
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

static int compute_transform(struct registration_args *regargs, struct starfinder_data *sf_args, gboolean *included, int *failed, float *fwhm, float *roundness, const float *B, gboolean verbose) {
	regdata *current_regdata = star_align_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata) return -1;
	int nb_ref_stars = sf_args->nb_stars[regargs->seq->reference_image];
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
			int not_matched = star_match_and_checks(sf_args->stars[regargs->seq->reference_image], sf_args->stars[i],
					sf_args->nb_stars[i], regargs, filenum, &H);
			if (not_matched) {
#ifdef _OPENMP
#pragma omp atomic
#endif
				nbfail++;
				included[i] = FALSE;
				continue;
			}
#ifdef _OPENMP
#pragma omp critical
#endif
			if (verbose) print_alignment_results(H, filenum, fwhm[i], roundness[i], "px");

		}
#ifdef _OPENMP
#pragma omp atomic
#endif
		nb_aligned++;
		current_regdata[i].roundness = roundness[i];
		current_regdata[i].fwhm = fwhm[i];
		current_regdata[i].weighted_fwhm = 2 * fwhm[i]
			* ((double)nb_ref_stars - sf_args->nb_stars[i])
			/ (double)nb_ref_stars + fwhm[i];
		current_regdata[i].background_lvl = B[i];
		current_regdata[i].number_of_stars = sf_args->nb_stars[i];
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
		cvTransfPoint(&currcenter.x, &currcenter.y,regargs->seq->regparam[regargs->layer][i].H, Href);
		cogx += currcenter.x;
		cogy += currcenter.y;
		n++;
	}
	cogx /= (double)n;
	cogy /= (double)n;
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
		cvTransfPoint(&currcenter.x, &currcenter.y, regargs->seq->regparam[regargs->layer][i].H, Htransf);
		dx = currcenter.x - center.x;
		dy = currcenter.y - center.y;
		dist[i] = (float)sqrt(dx * dx + dy * dy);
	}

}

// returns the index of the minimum element of float array arr of size nb
// if gboolean mask array is passed, it only includes elements where mask is TRUE
// if *vval is passed, it contains the best value
static int minidx(const float *arr, const gboolean *mask, int nb, float *val) {
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
	struct starfinder_data *sf_args = calloc(1, sizeof(struct starfinder_data));
	sf_args->im.from_seq = regargs->seq;
	sf_args->layer = regargs->layer;
	sf_args->max_stars_fitted = regargs->max_stars_candidates;
	sf_args->stars = calloc(regargs->seq->number, sizeof(psf_star **));
	sf_args->nb_stars = calloc(regargs->seq->number, sizeof(int));
	sf_args->update_GUI = FALSE;
	sf_args->already_in_thread = TRUE;
	sf_args->process_all_images = !regargs->filters.filter_included;
	sf_args->save_to_file = !regargs->no_starlist;
	float *fwhm = NULL, *roundness = NULL, *A = NULL, *B = NULL, *Acut = NULL, *scores = NULL;
	float *dist = NULL;
	// local flag (and its copy)accounting both for process_all_frames flag and collecting failures along the process
	gboolean *included = NULL, *tmp_included = NULL;
	// local flag to make checks only on frames that matter
	gboolean *meaningful = NULL;
	int retval = 0;
	struct timeval t_start, t_end;

	if (!sf_args->stars || !sf_args->nb_stars) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	gettimeofday(&t_start, NULL);
	if (apply_findstar_to_sequence(sf_args)) {
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
		if (!sf_args->stars[i]) {
			// star finder failed, we exclude the frame from the sequence
			// need this check to be after filters.filter_included otherwise, we would report false failures
			regargs->seq->imgparam[i].incl = FALSE;
			failed++;
			continue;
		}
		included[i] = TRUE;
		float FWHMx, FWHMy;
		char *units;
		FWHM_stats(sf_args->stars[i], sf_args->nb_stars[i], regargs->seq->bitpix, &FWHMx, &FWHMy, &units, B + i, Acut + i, AMPLITUDE_CUT);
		fwhm[i] = FWHMx;
		roundness[i] = FWHMy/FWHMx;
		if (sf_args->nb_stars[i] > maxstars) maxstars = sf_args->nb_stars[i];
	}


	if (maxstars == sf_args->max_stars_fitted) {
		siril_log_message(_("The number of stars has capped, readapting threshold and filtering\n"));
		float Athreshold = 0.0f;
		// Determine the updated A threshold
		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			if (sf_args->nb_stars[i] == sf_args->max_stars_fitted && Acut[i] > Athreshold) Athreshold = Acut[i];
		}
		// Filter the whole series against new amplitude threshold
		maxstars = maxstars * (1 - AMPLITUDE_CUT);

		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			sf_args->stars[i] = filter_stars_by_amplitude(sf_args->stars[i], Athreshold, &sf_args->nb_stars[i]);
			if (!sf_args->stars[i]) {
				regargs->seq->imgparam[i].incl = FALSE;
				included[i] = FALSE;
				failed++;
				continue;
			}
			float FWHMx, FWHMy;
			char *units;
			FWHM_stats(sf_args->stars[i], sf_args->nb_stars[i], regargs->seq->bitpix, &FWHMx, &FWHMy, &units, B + i, NULL, 0.);
			// we only rank images with at least half the maximum number of stars (and we count them as meaningful)
			scores[i]  = (sf_args->nb_stars[i] >= maxstars / 2) ? 2. * FWHMx * (maxstars - sf_args->nb_stars[i]) / maxstars + FWHMx : FLT_MAX;
			if (sf_args->nb_stars[i] >= maxstars / 2) meaningful[i] = TRUE;
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
			scores[i]  = (sf_args->nb_stars[i] >= maxstars / 2) ? 2. * FWHMx * (maxstars - sf_args->nb_stars[i]) / maxstars + FWHMx : FLT_MAX;
			if (sf_args->nb_stars[i] >= maxstars / 2) meaningful[i] = TRUE;
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
	int nb_valid_frames = 0, nb_meaningful = 0;
	for (int i = 0; i < regargs->seq->number; i++) {
		if (included[i]) nb_valid_frames++;
		if (meaningful[i]) nb_meaningful++;
	}
	int best_indexes[MAX_TRIALS_2PASS];
	int nb_aligned[MAX_TRIALS_2PASS];
	best_indexes[trials] = best_index;
	float allowable_dist = (float)regargs->seq->imgparam[regargs->reference_image].rx * MAX_SHIFT_RATIO;
	int tmp_failed;

	while (trials < MAX_TRIALS_2PASS) {
		tmp_failed = failed;
		for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
		nb_aligned[trials] = compute_transform(regargs, sf_args, tmp_included, &tmp_failed, fwhm, roundness, B, FALSE);
		// if number of aligned frames is less than half the number of meaningful frames (those with enough stars)
		// we have chosen a reference which is not framed well enough to align the sequence (the computed cog is probably meaningless as well)
		// we set its score to FLT_MAX and start again with the next best frame
		// we do not remove it from the meaningful frames count or included ones, as it still may be worth aligning with the new ref
		if (nb_aligned[trials] < nb_meaningful / 2) {
			scores[best_index] = FLT_MAX;
			siril_log_message(_("Trial #%d: After sequence alignement, image #%d could not align more than half of the frames, recomputing\n"), trials + 1, reffilenum);
			best_index = minidx(scores, included, regargs->seq->number, NULL);
			regargs->seq->reference_image = best_index;
			reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
			trials++;
			if (trials < MAX_TRIALS_2PASS) {
				siril_log_message(_("Trial #%d: After sequence analysis, we are choosing image %d as new reference for registration\n"), trials + 1, reffilenum);
				best_indexes[trials] = best_index;
			}
		} else { // not necessary but a simple to have print_alignment_results
			tmp_failed = failed;
			for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
			compute_transform(regargs, sf_args, tmp_included, &tmp_failed, fwhm, roundness, B, TRUE);
			break;
		}
	}

	if (trials == MAX_TRIALS_2PASS) {
	// We have tried many times over, we need to select the "best" frame
	// even though it cannot align more than half of the good images
		int best_try = -1;
		int maxreg = 0;
		siril_log_message(_("After %d trials, no reference image could align more than half of the frames, selecting best candidate\n"), MAX_TRIALS_2PASS);
		for (int i = 0; i < MAX_TRIALS_2PASS; i++) {
			siril_log_message(_("Trial #%d: Reference image: #%d - Frames aligned: %d\n"), i + 1, regargs->seq->imgparam[best_indexes[i]].filenum, nb_aligned[i]);
			if (nb_aligned[i] > maxreg) {
				best_try = i;
				maxreg = nb_aligned[i];
			}
		}
		if (maxreg <= 1) goto free_all;
		best_index = best_indexes[best_try];
		regargs->seq->reference_image = best_index;
		reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
		siril_log_message(_("After sequence analysis, we are choosing image %d as new reference for registration\n"), reffilenum);
		tmp_failed = failed;
		for (int i = 0; i < regargs->seq->number; i++) tmp_included[i] = included[i];
		compute_transform(regargs, sf_args, tmp_included, &tmp_failed, fwhm, roundness, B, TRUE);
	}
	// and we copy back to the initial arrays
	for (int i = 0; i < regargs->seq->number; i++) included[i] = tmp_included[i];
	failed = tmp_failed;

	// 4. check distance to cog
	compute_dist(regargs, dist, included);
	//if larger than cog, we recompute a score accoutning for the distance
	if (dist[best_index] > allowable_dist) {
		siril_log_message(_("After sequence alignement, image %d is too far from the sequence cog, recomputing\n"), reffilenum);
		siril_debug_print(_("Distance to cog is  %2.1fpx while threshold is set at %2.1fpx\n"), dist[best_indexes[trials]], allowable_dist);
		float FWHMx;
		int new_best_index = -1;
		for (int i = 0; i < regargs->seq->number; i++) {
			if (!included[i]) continue;
			FWHMx = fwhm[i];
			// images are now scored ONLY if they are within allowable distance from cog
			scores[i]  = (sf_args->nb_stars[i] >= maxstars / 2 && dist[i] < allowable_dist) ? 2. * FWHMx * (maxstars - sf_args->nb_stars[i]) / maxstars + FWHMx: FLT_MAX;
		}
		new_best_index = minidx(scores, included, regargs->seq->number, NULL);
		if (new_best_index != best_index && new_best_index > -1) { // do not recompute if none or same is found (same should not happen)
			regargs->seq->reference_image = new_best_index;
			reffilenum = regargs->seq->imgparam[new_best_index].filenum;	// for display purposes
			siril_log_message(_("After sequence analysis, we are choosing image %d as new reference for registration\n"), reffilenum);
			// back to 3. compute the transforms and store them in regparams
			compute_transform(regargs, sf_args, included, &failed, fwhm, roundness, B, TRUE);
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

	// images may have been excluded but selnum wasn't updated
	fix_selnum(regargs->seq, FALSE);
	siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->seq->selnum);
	clear_stars_list(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

free_all:
	for (int i = 0; i < regargs->seq->number; i++)
		free_fitted_stars(sf_args->stars[i]);
	free(sf_args->stars);
	free(sf_args->nb_stars);
	free(sf_args);
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

