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

/* This is the main file for the multi-point planetary processing */

/**************************** DOCUMENTATION ****************************
 * I. THE ALGORITHM
 * 1. Analysis of the sequence
 * -> provides global registration and quality value for all frames
 * Uses the register_cog (centre of gravity) registration function that
 * computes both values (see regcog_image_hook() in registration/cog.c)
 * 	image quality -> seq->regparam[layer][frame].quality
 * 	image shift -> seq->regparam[layer][frame].shiftx|y
 * 	selects best frame -> seq->reference_image
 * - image shifts are given relative to image of highest quality
 * - quality is normalized between 0 and 1, defaults to -1
 * 
 * 2. Creation of the reference frame used for display and placing AP
 * -> creates the reference image, saves it and displays it immediately
 * Makes a sum stacking of the images with quality above 0.75 (variable number
 * of images), using regparam shifts computed above
 * - image is also copied refimage, making it available for the next processing
 *   steps even if not displayed
 *
 * Operations 1. and 2. are done in the analysis thread that runs the function
 * sequence_analysis_thread_func()
 * 
 * 3. Adding the alignment points (AP), a.k.a. stacking zones, a.k.a. boxes
 * -> User adds boxes manually or uses the automatic placement algorithm
 *
 * 4. Evaluating stacking zones' quality
 * -> evaluates the quality of each zone for each image and keep it in memory
 * Zones are extracted from frames with coordinates centred on
 * (zone->centre.x - seq->regparam[layer][frame].shiftx,
 *  zone->centre.y + seq->regparam[layer][frame].shifty) and with a size
 * depending on the evaluated zone. Zones are always square.
 * 	zone quality -> zone->mpregparam[frame].quality
 * Zone quality is then normalized for all zones of all frames.
 * This is done in the_multipoint_quality_analysis() below.
 * From this step, quality for a zone over the sequence can be displayed in the
 * GUI, in addition to the global image quality. This step is done in
 * the_multipoint_processing(), the main mpp processing function, but can also
 * be triggered by the user before launching the main processing.
 * 
 * 5. Aligning stacking zones
 * -> computes shifts for all zones of all frames relative to the same zone
 *    from the reference image, using global registration data as a starting
 *    point.
 * Zones are extracted from frames with coordinates centred on
 * (zone->centre.x - seq->regparam[layer][frame].shiftx,
 *  zone->centre.y + seq->regparam[layer][frame].shifty) and with a size
 * depending on the evaluated zone.
 * The resulting shifts are added to the regparam shifts and saved like that:
 * 	zone->mpregparam[frame].x =  shiftx - regparam[frame].shiftx;
 * 	zone->mpregparam[frame].y = -shifty + regparam[frame].shifty;
 * We could do it only for the best zones to save time, since at this stage we
 * already have filtering information, but if the registration of a zone fails,
 * the next best one would be used in stacking.
 * This is currently done with the ECC algorithm, in
 * the_multipoint_ecc_registration(). Before that it was done with DFT phase
 * correlation in the_multipoint_dft_registration().
 * 
 * 6. Preparing the background image
 * -> an optional step that creates the image that will be used where no
 *    stacking zone has been defined
 * This could be the reference image, but we tried something different with a
 * global multi-point stacked image. It is computed by taking the best frames
 * as given by the global quality in seq->regparam and stacking a normalized
 * reference image that takes into account the displacement of zones with a
 * weighted distance (called barycentric) approach. Results are in general
 * slightly better than the reference image but may show some artifacts.
 * This is defined in the_global_multipoint_barycentric_sum_stacking().
 * The 5 nearest zones are computed for each pixel. Then for each pixel, the
 * shift for each zone (mpregparam[frame].x) is weighted with the distance of
 * the zone, giving a shift that follows the context.
 * Stacked image is built by taking pixels of each image with + shiftx and -
 * shifty shifts.
 * - image data is saved in gfit(ushort) and global_image(double)
 * 
 * 7. Stacking the final image
 * -> for each zone, stacks the best image of the sequence and blends them with
 *    the background image
 * For each zone, the best x% of images are selected. If an image is part of
 * the best of a zone, its zone is extracted and added to the result.
 * Zones are extracted from frames with coordinates centred on
 * (zone->centre.x - seq->mpregparam[frame].x,
 *  zone->centre.y + seq->mpregparam[frame].y) and with a size
 * depending on the evaluated zone.
 * Stacking result is normalized with the number of contributions per pixel,
 * since the number of zones stacked on it can vary. If there were no
 * contribution of zones on a pixel, the background image is used. If there are
 * less than a threshold (nb_images / 3 + 1), the result is blended with the
 * background image.
 * This is defined in the_local_multipoint_sum_stacking().
 * - resulting image is multiplied by 65535 and saved in gfit (the image that
 *   is displayed)
 */

 
/* Determining zone placement parameters for automation:
 *  - AP position: given by some feature detection algorithm.
 *  - zone size: can it be guessed from the quality values? We could count the
 *  number of pixels above the background level to estimate the size of the
 *  planet in pixels. Zone size should depend on image quality and noise. When
 *  an AP is placed, we could iteratively grow the zone around it until a
 *  minimum quality or contrast is reached.
 *  - overlap: the number of pixels that are common between two zones and of
 *  the number of pixel in the zones, as a percentage. A global setting of
 *  turbulence strength could drive this value, as more turbulence causes shear
 *  to appear at zones' borders. Overalap blends them and appear better. The
 *  value should not be higher than 50%, which would mean that the centre of
 *  two zones are both in the two zones. It's probably wiser to have the AP
 *  only in one zone.
 *  - extra read for zones: when trying to align zones, even if we read them in
 *  images taking into account the shifts from the global registration, they
 *  might not perfectly fit. This value indicates how many pixels extra in the
 *  half-side of the zones read from the reference image we need to read. This
 *  depends probably on the same data as the overlap, with emphasis on the
 *  standard deviation of the set of shifts. This depends on the registration
 *  algorithm used, some may not allow different sample sizes.
 *
 *  In automatic mode, if the user changes the zone size, we could update the
 *  overlap amount.
 */

/* TODO:
 * handle coloured images
 * extra read for zones for the registration, does it make sens DFT-wise?
 *   probably not
 * sub-pixel phase correlation is possible and explained in the method section here
 *   https://en.wikipedia.org/wiki/Phase_correlation but "Subpixel methods are also
 *   particularly sensitive to noise in the images"
 * allow the ref/first image not to be displayed when loading a sequence if we
 *   display the planetary reference image
 * demosaicing setting: manage prepro in on_prepro_button_clicked
 * fix the new arg->layer and filters in processing and elsewhere
 * GUI mode switch: maybe next start, maybe fork+exec?
 * remove menuitemcolor and menuitemgray in planetary mode
 * chain hi/lo/mode when RGB vport is modified
 * FFTW-related, if we keep DFT registration it for mpp:
 *    make a background thread that precomputes FFTW wisdom, it's too slow for big zones
 *    restrict zone size (see FFTW comment below) to some values
 * restrict zones inside the image
 * allow the sequence list to be ellipsized because it can be huge -> feature not found!
 * displays 'No sequence selected' when selecting a sequence that has no global regdata
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <math.h>
#include <complex.h>	// to include before fftw3.h
#include <fftw3.h>
#include <time.h>

#include "planetary.h"
#include "core/siril.h"
#include "core/processing.h"
#include "core/proto.h"
#include "algos/quality.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/progress_and_log.h"
#include "gui/planetary_callbacks.h"
#include "gui/callbacks.h"
#include "gui/plot.h"
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"

#define FFTW_WISDOM_FILE "fftw_wisdom"
#define DEBUG_MPP
//#define USE_SEQUENCE_REF_IMAGE

static fits refimage;
static char *refimage_filename;

static gpointer sequence_analysis_thread_func(gpointer p);
static gboolean end_reference_image_stacking(gpointer p);

static int the_multipoint_quality_analysis(struct mpr_args *args);

static int the_multipoint_ecc_registration(struct mpr_args *args);
static int the_multipoint_dft_registration(struct mpr_args *args);

static int the_global_multipoint_barycentric_sum_stacking(struct mpr_args *args);
static int the_local_multipoint_sum_stacking(struct mpr_args *args);

static int get_number_of_zones();
static int copy_image_zone_to_fftw(fits *fit, const stacking_zone *zone, fftw_complex *dest,
		int layer);
static int copy_image_zone_to_buffer(fits *fit, const stacking_zone *zone, WORD *dest, int layer);
static void add_image_zone_to_stacking_sum(fits *fit, const stacking_zone *zone, int frame,
		unsigned long *sum[3], int *count[3]);
static void add_buf_zone_to_stacking_sum(WORD *buf, int layer, const stacking_zone *zone,
		int frame, unsigned long *sum[3], int *count[3], unsigned int rx, unsigned int ry);
static BYTE * sort_zones_quality(int zone_idx, int nb_best, int nb_seq_images);

#if defined DEBUG_MPP || defined DEBUG_MPP1
static void save_buffer_tmp(int frame_index, int zone_idx, WORD *buffer, int square_size);
#endif

/* First step: running the global registration and building the reference image */
void on_planetary_analysis_clicked(GtkButton *button, gpointer user_data) {
	if (!sequence_is_loaded()) return;
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));

	struct registration_args *reg_args;
	reg_args = calloc(1, sizeof(struct registration_args));
	reg_args->func = register_cog;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->process_all_frames = FALSE;
	reg_args->follow_star = FALSE;
	reg_args->matchSelection = FALSE;
	reg_args->translation_only = TRUE;
	reg_args->x2upscale = FALSE;
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;

	char *msg = siril_log_color_message(_("Starting sequence analysis\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	set_cursor_waiting(TRUE);

	start_in_new_thread(sequence_analysis_thread_func, reg_args);
}

static gpointer sequence_analysis_thread_func(gpointer p) {
	struct registration_args *reg_args = (struct registration_args *) p;

	// COG registration
	reg_args->retval = reg_args->func(reg_args);

	if (reg_args->retval)
		return GINT_TO_POINTER(-1);

	/* sequence quality evaluation is over, now we generate the reference
	 * image by taking a few of the best images */

	struct stacking_args *stack_args = calloc(1, sizeof(struct stacking_args));
	stack_args->method = stack_summing_generic;
	stack_args->seq = reg_args->seq;
	stack_args->reglayer = reg_args->layer;
	stack_args->filtering_criterion = stack_filter_quality;
	stack_args->filtering_parameter = 0.75;	// not the right way do to it, sorting is
	stack_args->nb_images_to_stack = compute_nb_filtered_images(stack_args);
	stack_args->image_indices = malloc(stack_args->nb_images_to_stack * sizeof(int));
	stack_fill_list_of_unfiltered_images(stack_args);
	snprintf(stack_args->description, sizeof stack_args->description,
			_("Creating the reference image for the multi-point planetary processing"));
	stack_args->output_filename = get_reference_image_name(reg_args->seq, reg_args->layer);
	stack_args->output_overwrite = TRUE;
	gettimeofday(&stack_args->t_start, NULL);
	free(reg_args);

	stack_args->retval = stack_summing_generic(stack_args);
	if (stack_args->retval) {
		free(stack_args);
		return GINT_TO_POINTER(-1);
	}

	/* the reference image is now in gfit */
	copyfits(&gfit, &refimage, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	siril_add_idle(end_stacking, stack_args);
	siril_add_idle(end_reference_image_stacking, stack_args);

	return NULL;
}

static gboolean end_reference_image_stacking(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;

	/* adapting automatic settings based on registration data */
	regdata *regparam = args->seq->regparam[args->reglayer];
	WORD *shifts = malloc(args->seq->number * sizeof(WORD));
	int i, nb_images;
	for (i = 0, nb_images = 0; i < args->seq->number; i++) {
		if (regparam[i].quality < 0.0) continue;
		shifts[nb_images++] = (WORD)roundf_to_int(
				fabs(regparam[i].shiftx) + fabs(regparam[i].shifty));
	}

	double mean, sigma;
	int status = 0;
	FnMeanSigma_ushort(shifts, nb_images, 0, 0, NULL, &mean, &sigma, &status);
	// the set of shifts is converted to WORD which will lose some precision but it's
	// probably not significant for this use.

	fprintf(stdout, "Average shifts for the sequence: %f (sigma: %f)\n", mean, sigma);

	/* things from end_register_idle */
	set_layers_for_registration();	// update display of available reg data
	drawPlot();

	free(args->output_filename);
	free(args);
	return FALSE;
}

// 2be3d
char *get_reference_image_name(sequence *seq, int layer) {
	char *seqname = seq->seqname;
	char *output_filename = malloc(strlen(seqname) + 20);
	sprintf(output_filename, "%s%sreference%d%s", seqname,
			ends_with(seqname, "_") ?  "" :
			(ends_with(seqname, "-") ? "" : "_"),
			layer, com.ext);
	return output_filename;
}

int refimage_is_set() {
	return refimage.naxes[0] > 0 && refimage_filename;
}

const fits *get_refimage() {
	return &refimage;
}

const char *get_refimage_filename() {
	return (const char *)refimage_filename;
}

void update_refimage_on_layer_change(sequence *seq, int layer) {
	if (refimage_is_set())
		clearfits(&refimage);
	if (layer == -1) {
		if (refimage_filename) free(refimage_filename);
		refimage_filename = NULL;
		refimage.naxes[0] = 0;
		return;
	}
	char *ref_name = get_reference_image_name(seq, layer);
	if (is_readable_file(ref_name) && !readfits(ref_name, &refimage, NULL)) {
		siril_log_message(_("A previously computed reference image was found"
				       " for layer %d, using it\n"), layer);
		if (refimage_filename) free(refimage_filename);
		refimage_filename = ref_name;
		display_refimage_if_needed();
	} else {
		free(ref_name);
		refimage.naxes[0] = 0;
	}
}

/*********************** ANALYSIS *************************/
gpointer the_multipoint_analysis(gpointer ptr) {
	struct mpr_args *args = (struct mpr_args*)ptr;
	int retval = the_multipoint_quality_analysis(args);
	siril_add_idle(end_generic, NULL);
	free(args);
	return GINT_TO_POINTER(retval);
}

/* evaluates quality of each zone for each image and saves it in mpregparam */
static int the_multipoint_quality_analysis(struct mpr_args *args) {
	int retval = 0, frame, zone_idx, abort = 0;

	int nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}
	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;
	WORD **buffer = calloc(nb_zones, sizeof(WORD *));
	double *max = calloc(nb_zones, sizeof(double));

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "analysing zones for image %d\n", frame);

		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			stacking_zone shifted_zone = { .centre =
				{ .x = zone->centre.x - regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };
			int side = round_to_int(zone->half_side * 2.0);

			if (!zone->mpregparam) {
				zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
				if (!zone->mpregparam) {
					fprintf(stderr, "memory allocation failed\n");
					abort = 1;
					break;
				}
			}

			// copy the zone to a buffer and evaluate quality,
			// cannot be done in-place
			if (!buffer[zone_idx]) {
				/* this is not thread-safe, a buffer is required per thread, it
				 * can be done with a single allocation with round-robin
				 * distribution of zones across threads */
				buffer[zone_idx] = malloc(side * side * sizeof(WORD));
			}

			if (!copy_image_zone_to_buffer(&fit, &shifted_zone,
						buffer[zone_idx], args->layer)) {
				zone->mpregparam[frame].quality =
					QualityEstimateBuf(buffer[zone_idx], side, side);
				if (max[zone_idx] < zone->mpregparam[frame].quality)
					max[zone_idx] = zone->mpregparam[frame].quality;
			}
		}

		clearfits(&fit);
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	// normalization
	for (frame = 0; frame < args->seq->number; frame++) {
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			stacking_zone *zone = &com.stacking_zones[zone_idx];
			zone->mpregparam[frame].quality /= max[zone_idx];
		}
	}	

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		if (buffer[zone_idx])
			free(buffer[zone_idx]);
	}
	free(buffer);
	free(max);

	return retval;
}

/*********************** PROCESSING *************************/
gpointer the_multipoint_processing(gpointer ptr) {
	struct mpr_args *args = (struct mpr_args*)ptr;
	int retval;
	if (!com.stacking_zones || com.stacking_zones[0].centre.x == -1.0) {
		siril_log_message("Place some zones on the image before starting the processing\n");
		return GINT_TO_POINTER(-1);
	}
	/* make sure the zones analysis has been done above */
	if (!com.stacking_zones[0].mpregparam ||
			com.stacking_zones[0].mpregparam[0].quality < 0.0 ||
			com.stacking_zones[0].mpregparam[0].quality > 1.0) {
		// TODO: do this in the registration instead to improve performance
		retval = the_multipoint_quality_analysis(args);
		if (retval) return GINT_TO_POINTER(retval);
	}
	siril_log_message("zones will be stacked using %s\n",
			args->using_homography ? "homography" : "translation only");

	/* multi-point registration: compute the local shifts */
	//retval = the_multipoint_dft_registration(args);
	retval = the_multipoint_ecc_registration(args);
	if (retval) return GINT_TO_POINTER(retval);

	/* reference image stacking: stack global images with local shifts.
	 * This image is used for areas of the result where there is no zone defined */
	retval = the_global_multipoint_barycentric_sum_stacking(args);
	if (retval) return GINT_TO_POINTER(retval);

	/* multi-point stacking: stack zones with local shifts and creates the result in gfit */
	retval = the_local_multipoint_sum_stacking(args);

	/* registering the generic stacking idle with the required args */
	struct stacking_args *stackargs = malloc(sizeof(struct stacking_args));;
	stackargs->retval = retval;
	stackargs->output_filename = args->output_filename;
	stackargs->output_overwrite = args->output_overwrite;
	stackargs->seq = args->seq;
	// TODO: nb_images_to_stack for summary output
	siril_add_idle(end_stacking, stackargs);

	return GINT_TO_POINTER(retval);
}

static int the_multipoint_ecc_registration(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, abort = 0;
	WORD **references;	// an array of pointers because their size varies
	stacking_zone *zone;

	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;

	nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}

	/* read zones of the reference image */
	references = malloc(nb_zones*sizeof(WORD *));
	if (!references) {
		fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
		return -1;
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	fits reffits = { 0 };
	if (seq_read_frame(args->seq, args->seq->reference_image, &reffits))
		return -1;
	fprintf(stdout, "Using the reference image from sequence (%d) instead "
			"of the stacked reference\n", args->seq->reference_image);
#endif
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		zone = &com.stacking_zones[zone_idx];
		int side = round_to_int(zone->half_side * 2.0);
		references[zone_idx] = malloc(side * side * sizeof(WORD));
		if (!references[zone_idx]) {
			fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
			return -1;
		}
#ifdef USE_SEQUENCE_REF_IMAGE
		copy_image_zone_to_buffer(&reffits, zone, references[zone_idx], args->layer);
#else
		copy_image_zone_to_buffer(&refimage, zone, references[zone_idx], args->layer);
#endif
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	clearfits(&reffits);
#endif

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		/* TODO: this filtering is required for global mode, not for local */
		/*if (!args->filtering_criterion(args->seq, args->layer,
					frame, args->filtering_parameter) || abort)
			continue;*/

		if (seq_read_frame(args->seq, frame, &fit) || image_find_minmax(&fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "\naligning zones for image %d\n", frame);

		/* reading zones in the reference image */
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			reg_ecc reg_param = { 0 };
			if (abort) continue;
			zone = &com.stacking_zones[zone_idx];
			stacking_zone shifted_zone = { .centre =
				{ .x = zone->centre.x - regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };
			int side = round_to_int(zone->half_side * 2.0);
			WORD *buffer = malloc(side * side * sizeof(WORD));	// TODO: prealloc
			if (!buffer) {
				fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
				abort = 1;
				continue;
			}

			// read the zone data to buffer
			copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
#ifdef DEBUG_MPP1
			save_buffer_tmp(frame, zone_idx, buffer, side);
#endif
			if (!zone->mpregparam) {
				zone->mpregparam = calloc(args->seq->number, sizeof(struct ap_regdata));
				if (!zone->mpregparam) {
					fprintf(stderr, "memory allocation failed\n");
					abort = 1;
					continue;
				}
			}

			int regretval; 
			if (args->using_homography) {
				if (!zone->mpregparam[frame].transform)
					zone->mpregparam[frame].transform = malloc(sizeof(Homography));
				regretval = ecc_find_transform_buf(references[zone_idx], buffer,
						side, fit.maxi, &reg_param,
						zone->mpregparam[frame].transform);
				if (regretval) {
					free(zone->mpregparam[frame].transform);
					zone->mpregparam[frame].transform = NULL;
					// TODO can we use the downsampled results? we need a backup
				}
			} else {
				regretval = ecc_find_translation_buf(references[zone_idx], buffer,
						side, fit.maxi, &reg_param);
				if (zone->mpregparam[frame].transform)
					free(zone->mpregparam[frame].transform);
				zone->mpregparam[frame].transform = NULL;
			}
			if (regretval) {
				fprintf(stdout, "ECC alignment failed for full def zone %d of frame %d\n", zone_idx, frame);
				// setting the quality to -1 will remove it from stacking, but it
				// may trigger the sequence quality reanalysis
				zone->mpregparam[frame].quality = -1.0;
			} else {
				// in ecc registration + sum stacking it's - and -
				zone->mpregparam[frame].x = reg_param.dx + regparam[frame].shiftx;
				zone->mpregparam[frame].y = -reg_param.dy + regparam[frame].shifty;
				fprintf(stdout, "frame %d, zone %d local shifts: %f,%f\n",
						frame, zone_idx, reg_param.dx, reg_param.dy);
			}

			free(buffer);
		}
		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}
	return abort;
}

/* about FFTW: http://www.fftw.org/fftw3_doc/Introduction.html
 * The standard FFTW distribution works most efficiently for arrays whose size
 * can be factored into small primes (2, 3, 5, and 7), and otherwise it uses a
 * slower general-purpose routine.
 * Maybe we need to constrain the zone size to these numbers, that's not too
 * many that do not register (11 13 17 19 23 31 37 41 43 47 53 59 61 67 71 ...)
 * wait, do they need to decompose down to only 2 3 5 and 7? If so the list is
 * longer.
 *
 * Applying the phase correlation method to a pair of images produces a third
 * image which contains a single peak. The location of this peak corresponds to
 * the relative translation between the images.
 * The Fourier-Mellin transform extends phase correlation to handle images
 * transformed by both translation and rotation.
 */
static int the_multipoint_dft_registration(struct mpr_args *args) {
	int zone_idx, nb_zones, frame;
	int abort = 0;
	int retval = 0;
	stacking_zone *zone;
	fftw_complex **ref, **in, **out, **convol;
	fftw_plan *fplan, *bplan;	// forward and backward plans
	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;

	nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}
	ref = malloc(nb_zones * sizeof(fftw_complex*));
	in = malloc(nb_zones * sizeof(fftw_complex*));
	out = malloc(nb_zones * sizeof(fftw_complex*));
	convol = malloc(nb_zones * sizeof(fftw_complex*));
	fplan = malloc(nb_zones * sizeof(fftw_plan));
	bplan = malloc(nb_zones * sizeof(fftw_plan));

	/* reading zones in the reference image, single threaded init */
	fprintf(stdout, "loading reference zones\n");
	gchar *wisdom_file = get_configdir_file_path(FFTW_WISDOM_FILE);
	if (wisdom_file) {
		if (fftw_import_wisdom_from_filename(wisdom_file))
			fprintf(stdout, "FFTW wisdom restored\n");
	}

#ifdef USE_SEQUENCE_REF_IMAGE
	fits reffits = { 0 };
	if (seq_read_frame(args->seq, args->seq->reference_image, &reffits))
		return -1;
	fprintf(stdout, "Using the reference image from sequence (%d) instead "
			"of the stacked reference\n", args->seq->reference_image);
#endif

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		zone = &com.stacking_zones[zone_idx];
		int side = round_to_int(zone->half_side * 2.0);
		int nb_pixels = side * side;
		// allocate aligned DFT buffers with fftw_malloc
		ref[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		in[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		out[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		convol[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);

		start_timer();
		fplan[zone_idx] = fftw_plan_dft_2d(side, side, ref[zone_idx], out[zone_idx], FFTW_FORWARD, FFTW_PATIENT);
		bplan[zone_idx] = fftw_plan_dft_2d(side, side, convol[zone_idx], out[zone_idx], FFTW_BACKWARD, FFTW_PATIENT);
		fprintf(stdout, "plan %d creation time: %ld microsec\n", zone_idx,
				stop_timer_elapsed_mus());
		// the plan creation time should be long only the first time a zone size is used

#ifdef USE_SEQUENCE_REF_IMAGE
		copy_image_zone_to_fftw(&reffits, zone, ref[zone_idx], args->layer);
#else
		// use the refimage stacked from the best of sequence
		copy_image_zone_to_fftw(&refimage, zone, ref[zone_idx], args->layer);
#endif
		fftw_execute_dft(fplan[zone_idx], ref[zone_idx], in[zone_idx]);
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	clearfits(&reffits);
#endif

	// in out and convol can probably be freed here, or even allocated once and reused
	// for the above loop
	if (wisdom_file) {
		if (fftw_export_wisdom_to_filename(wisdom_file))
			fprintf(stdout, "FFTW wisdom saved\n");
		g_free(wisdom_file);
	}

	fftw_complex **zones = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **out2 = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **convol2 = calloc(nb_zones, sizeof(fftw_complex*));

	/* for each image, we read the zones with global shift and register them */
	/* for sequences that require demosaicing, seq_read_frame is usually the longest
	 * operation in this loop, so parallelizing the zone alignment does not help
	 * much, either with fftw3_threads or with OpenMP in the zone loop below. The
	 * best way to speed things up is probably to execute this loop in parallel.
	 * The fftw plan execution is thread-safe so it should not be a problem. */
	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		int i;
		if (abort) continue;
		/* TODO: this filtering is required for global mode, not for local */
		/*if (!args->filtering_criterion(args->seq, args->layer,
					frame, args->filtering_parameter) || abort)
			continue;*/

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "aligning zones for image %d\n", frame);

		/* reading zones in the reference image */
#ifdef _OPENMP
//#pragma omp parallel for num_threads(com.max_thread) private(zone_idx, zone) schedule(static) if((args->seq->type == SEQ_REGULAR && fits_is_reentrant()) || args->seq->type == SEQ_SER)
#endif
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			if (abort) continue;
			zone = &com.stacking_zones[zone_idx];
			int side = round_to_int(zone->half_side * 2.0);
			int nb_pixels = side * side;
			stacking_zone shifted_zone = { .centre =
				{ .x = zone->centre.x - regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };

			if (!zones[zone_idx]) {	// keep across images
			       zones[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       out2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       convol2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			}

			copy_image_zone_to_fftw(&fit, &shifted_zone, zones[zone_idx],
					args->layer);

			// forward transformation zone -> out2
			fftw_execute_dft(fplan[zone_idx], zones[zone_idx], out2[zone_idx]);

			// compute a new fourier domain image defined as the fourier
			// representation of the reference image on this zone * conj(out2)
			for (i = 0; i < nb_pixels; i++) {
				convol2[zone_idx][i] = in[zone_idx][i] * conj(out2[zone_idx][i]);
			}

			// backward transformation of this new image to out2
			fftw_execute_dft(bplan[zone_idx], convol2[zone_idx], out2[zone_idx]);

			// searching for the real part peak in out2, which is the shift
			// between the reference image and this image in this zone
			int shift = 0;
			for (i = 1; i < nb_pixels; ++i) {
				if (creal(out2[zone_idx][i]) > creal(out2[zone_idx][shift])) {
					shift = i;
				}
			}
			int shifty = shift / side;
			int shiftx = shift % side;
			if (shifty > zone->half_side) {
				shifty -= side;
			}
			if (shiftx > zone->half_side) {
				shiftx -= side;
			}

			/* shitfs for this image and this zone is (shiftx, shifty) + the global shifts */
			fprintf(stdout, "frame %d, zone %d adjustment shifts: %d,%d\n", frame, zone_idx, shiftx, shifty);
			zone->mpregparam[frame].x = (double)shiftx + regparam[frame].shiftx;
			zone->mpregparam[frame].y = (double)shifty - regparam[frame].shifty;
		}

		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

cleaning_all:
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		fftw_free(zones[zone_idx]);
		fftw_free(out2[zone_idx]);
		fftw_free(convol2[zone_idx]);

		fftw_destroy_plan(fplan[zone_idx]);
		fftw_destroy_plan(bplan[zone_idx]);
		fftw_free(in[zone_idx]);
		fftw_free(out[zone_idx]);
		fftw_free(ref[zone_idx]);
		fftw_free(convol[zone_idx]);
	}
	free(zones); free(out2); free(convol2);
	free(fplan); free(bplan);
	free(in); free(out); free(ref); free(convol);

	return retval;
}


struct weighted_AP {
	int zone_index;
	double distance;
	// double direction; someday
};

static void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP);
//static void filter_closest_list_owned(struct weighted_AP *list_for_this_point,
//		int max_AP, double owned_distance);

/* this function is a stacking sum that takes shifts from the multipoint registration
 * instead of taking them from the global registration. If it works well, it will have to
 * be exploded as a generic processing function. */
static int the_global_multipoint_barycentric_sum_stacking(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, abort = 0;
	regdata *regparam = args->seq->regparam[args->layer];
	struct weighted_AP *closest_zones_map;	// list of nb_closest_AP AP (fixed) for each pixel

	nb_zones = get_number_of_zones();
	if (nb_zones < args->nb_closest_AP)
		args->nb_closest_AP = nb_zones;
	closest_zones_map = malloc(args->seq->rx * args->seq->ry * args->nb_closest_AP * sizeof(struct weighted_AP));
	if (!closest_zones_map) {
		fprintf(stderr, "Stacking: memory allocation failure for zone mapping\n");
		return -1;
	}
	siril_log_message("Using %d closest zones for multipoint stacking refinement\n", args->nb_closest_AP);

	/* precompute the closest zones from each pixel */
	int x, y;
	for (y = 0; y < args->seq->ry; y++) {
		for (x = 0; x < args->seq->rx; x++) {
			int ap;
			struct weighted_AP *list_for_this_pixel = closest_zones_map + (x + y * args->seq->rx) * args->nb_closest_AP;
			for (ap = 0; ap < args->nb_closest_AP; ap++)	// init the struct
				list_for_this_pixel[ap].distance = -1.0;

			for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
				stacking_zone *zone = &com.stacking_zones[zone_idx];
				// compute the distance from the centre of the pixel
				double xdist = zone->centre.x - x + 0.5;
				double ydist = zone->centre.y - y + 0.5;
				double distance = sqrt(xdist * xdist + ydist * ydist);
				if (distance < args->max_distance)
					check_closest_list(list_for_this_pixel, distance, zone_idx, args->nb_closest_AP);
			}
			/* TODO: doesn't work, zone size has to be known, it has to be
			 * done in the main check */
			//filter_closest_list_owned(list_for_this_point, args->nb_closest_AP, args->own_distance_f * zone
		}
	}

	// init stacking data (copied from sum_stacking_prepare_hook)
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	unsigned long *sum[3];	// the new image's channels
	int *count[3];	// the new image's contributions count
	sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!sum[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		sum[1] = sum[0] + nbdata;	// index of green layer in sum[0]
		sum[2] = sum[0] + nbdata*2;	// index of blue layer in sum[0]
	} else {
		sum[1] = NULL;
		sum[2] = NULL;
	}

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		if (abort) continue;
		if (!args->filtering_criterion(args->seq, args->layer, frame, args->filtering_parameter))
			continue;

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "barycentre stacking for image %d\n", frame);

		/* for each image, we use the closest zones' shifts computed in the
		 * multi-point registration and weight them with their distance to get
		 * the coordinates of the pixel of the image to stack for each end-image
		 * pixel. The direction of each AP would also be important, the
		 * distribution of the AP having a role in the weighing. */
		int x, y;
		int pixel = 0;	// index in sum[0]
		for (y = 0; y < fit.ry; y++) {
			for (x = 0; x < fit.rx; x++) {
				struct weighted_AP *list_for_this_pixel = closest_zones_map + (x + y * args->seq->rx) * args->nb_closest_AP;
				int i;
				double total_weight = 0.0, sumx = 0.0, sumy = 0.0, shiftx, shifty;

				for (i = 0; i < args->nb_closest_AP; i++) {
					if (list_for_this_pixel[i].distance < 0.0)
						continue;
					double this_weight;
					if (list_for_this_pixel[i].distance <= 1.0)
						this_weight = args->max_distance;
					else this_weight = args->max_distance / list_for_this_pixel[i].distance;
					int zone_index = list_for_this_pixel[i].zone_index;
					struct ap_regdata *regparam_zone_image = 
						&com.stacking_zones[zone_index].mpregparam[frame];
					if (regparam_zone_image->quality < 0.0)
						continue;

					total_weight += this_weight;
					sumx += regparam_zone_image->x * this_weight;
					sumy += regparam_zone_image->y * this_weight;
				}
				// if zones are too far away, which will happen for pixels
				// away from the planet, just take the global shift
				if (total_weight == 0.0) {
					shiftx = regparam[frame].shiftx;
					shifty = regparam[frame].shifty;
				} else {
					shiftx = sumx / total_weight;
					shifty = sumy / total_weight;
				}

				// then do the regular sum stacking with shift
				int nx = round_to_int(x - shiftx);
				int ny = round_to_int(y - shifty);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					int ii = ny * fit.rx + nx;	// index in source image
					int layer;
					for (layer = 0; layer < args->seq->nb_layers; ++layer) {
#ifdef _OPENMP
//#pragma omp atomic
#endif
						sum[layer][pixel] += fit.pdata[layer][ii];
					}
				}
				++pixel;
			}
		}
		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	fprintf(stdout, "barycentre stacking ended, creating final image\n");

	if (!abort) {
		/* copying result into gfit and args->global_image
		 * code copied from sum_stacking_finalize_hook() */
		nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
		int i, layer;
		// find the max first
		unsigned long max = 0L;	// max value of the image's channels
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i=0; i < nbdata; ++i)
			if (sum[0][i] > max)
				max = sum[0][i];

		clearfits(&gfit);
		fits *fit = &gfit;
		if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers))
			return -1;
		gfit.hi = round_to_WORD(max);

		double ratio = 1.0;
		if (max > USHRT_MAX)
			ratio = USHRT_MAX_DOUBLE / (double)max;

		double norm_ratio = 1.0 / (double)max;
		args->global_image = malloc(sizeof(double) * nbdata);

		unsigned long* from = sum[0];
		WORD *to = gfit.data;
		for (i=0; i < nbdata; ++i) {
			if (ratio == 1.0)
				to[i] = round_to_WORD(from[i]);
			else	to[i] = round_to_WORD((double)(from[i]) * ratio);
			args->global_image[i] = from[i] * norm_ratio;
		}
	}

	free(sum[0]);

	return 0;
}

/* this function is a stacking sum that takes shifts from the multipoint registration instead
 * of taking them from the global registration.  Additionally, it takes only parts of images
 * defined as best for a zone instead of taking the best global images for the seuqence.
 * If it works well, it will have to be exploded as a generic processing function. */
static int the_local_multipoint_sum_stacking(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, nb_images, abort = 0;
	struct weighted_AP *closest_zones_map;	// list of nb_closest_AP AP (fixed) for each pixel

	nb_zones = get_number_of_zones();
	nb_images = round_to_int((double)args->seq->number * args->filtering_percent / 100.0);
	fprintf(stdout, "keeping %d best images for each zone\n", nb_images);
	
	/* 1. sort images indices from the list of best quality for zones *
	 * In com.stacking_zones[zone].mpregparam we have the normalized quality for each
	 * image for the concerned zone. To speed up the look-up, we create an index here.
	 * */
	BYTE **best_zones = malloc(nb_zones * sizeof(BYTE *));
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		best_zones[zone_idx] = sort_zones_quality(zone_idx, nb_images, args->seq->number);
	}

	/* 2. allocate stacking data *
	 * A sum and a contribution map (number of pixel added for each sum)
	 */
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	unsigned long *sum[3];	// the new image's channels
	int *count[3];	// the new image's contributions count
	sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	count[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!sum[0] || !count[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		sum[1] = sum[0] + nbdata;	// index of green layer in sum[0]
		sum[2] = sum[0] + nbdata*2;	// index of blue layer in sum[0]
		count[1] = count[0] + nbdata;
		count[2] = count[0] + nbdata*2;
	} else {
		sum[1] = NULL;
		sum[2] = NULL;
		count[1] = NULL;
		count[2] = NULL;
	}

	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		int zone_idx = 0;
		if (abort) continue;

		while (!best_zones[zone_idx++][frame] && zone_idx < nb_zones);
		if (zone_idx == nb_zones) {
			fprintf(stdout, "frame %d had only bad zones, skipping\n", frame);
			continue;
		}

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "reading best zones for image %d\n", frame);

		/* 3. extract best zones and stack them *
		 * From each image we add the pixels of the best zones to the sum and keep track
		 * of how many times we add each pixel to keep track of the average
		 */
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			if (!best_zones[zone_idx][frame])
				continue;

			stacking_zone *zone = &com.stacking_zones[zone_idx];
			// TODO: iterate over channels
			if (args->using_homography && zone->mpregparam[frame].transform) {
				// ECC registration with image transformation
				int side = round_to_int(zone->half_side * 2.0);
				WORD *buffer = malloc(side * side * sizeof(WORD)); // TODO: prealloc
				if (!buffer) {
					fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
					abort = 1;
					break;
				}

				// read the zone data to buffer
				regdata *regparam = args->seq->regparam[args->layer];
				stacking_zone shifted_zone = { .centre =
					{ .x = zone->centre.x - regparam[frame].shiftx,
						.y = zone->centre.y + regparam[frame].shifty },
					.half_side = zone->half_side };
				copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
				cvTransformBuf(buffer, side, zone->mpregparam[frame].transform);
				add_buf_zone_to_stacking_sum(buffer, 0, zone, frame,
						sum, count, fit.rx, fit.ry);
#ifdef DEBUG_MPP
				save_buffer_tmp(frame, zone_idx, buffer, side);
#endif
				free(buffer);
			} else if (!args->using_homography) {
				// DFT registration or ECC with translation
				add_image_zone_to_stacking_sum(&fit, zone, frame, sum, count);
#ifdef DEBUG_MPP
				// to see what's happening with the shifts, use this
				int side = round_to_int(zone->half_side * 2.0);
				stacking_zone shifted_zone = { .centre =
					{ .x = zone->centre.x - zone->mpregparam[frame].x,
						.y = zone->centre.y + zone->mpregparam[frame].y },
					.half_side = zone->half_side };

				WORD *buffer = malloc(side * side * sizeof(WORD));
				if (!buffer) {
					fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
					break;
				}

				copy_image_zone_to_buffer(&fit, &shifted_zone, buffer, args->layer);
				save_buffer_tmp(frame, zone_idx, buffer, side);
				free(buffer);
#endif
			}
			/* TODO someday: instead of copying data like this, create
			 * a mask from the zones, process the mask to unsharpen it
			 * and make it follow the path of the turbulence, and then
			 * use the mask to copy the data. */
		}

		clearfits(&fit);
		// this progress is not as bad as usual, but could be improved by counting
		// the number of zones instead of the number of images
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++)
		free(best_zones[zone_idx]);
	free(best_zones);

	fprintf(stdout, "multipoint stacking ended, creating final image\n");

	/* 4. compute averages for the zones and merge with the reference image for areas
	 * outside zones, store the result in gfit.
	 */
	if (!abort) {
		int i, minzones = nb_images / 3 + 1;	// at least 1
		int from_mpp = 0, from_both = 0, from_ref = 0;
		nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;

		// normalize the sum to the contribution count
		// find the max to normalize to 0..1
		double max = 0.0, invmax, *buf = malloc(nbdata * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i = 0; i < nbdata; i++) {
			buf[i] = (double)sum[0][i] / (double)count[0][i];
			if (buf[i] > max)
				max = buf[i];
		}

		invmax = 1.0 / max;
		for (i = 0; i < nbdata; i++)
			buf[i] *= invmax;

		for (i = 0; i < nbdata; i++) {
			if (count[0][i] > minzones) {
				// already in buffy
				from_mpp++;
			} else if (count[0][i] > 0) {
				buf[i] = (args->global_image[i] + buf[i]) / 2.0;
				from_both++;
			} else {
				buf[i] = args->global_image[i];
				from_ref++;
			}
		}
		fprintf(stdout, " pixels from best local zones: %d\n", from_mpp);
		fprintf(stdout, " pixels from best local zones mixed with reference: %d\n", from_both);
		fprintf(stdout, " pixels from reference image (global): %d\n", from_ref);

		// make a copy of the reference image to gfit to initialize it
		clearfits(&gfit);
		fits *fit = &gfit;
		if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers))
			return -1;
		gfit.hi = USHRT_MAX;

		double ratio = 1.0;
		ratio = USHRT_MAX_DOUBLE;

		for (i = 0; i < nbdata; i++)
			fit->data[i] = round_to_WORD(buf[i] * ratio);
		free(buf);
	}

	free(sum[0]);
	free(count[0]);
	free(args->global_image);
	return 0;
}


static int get_number_of_zones() {
	int i = 0;
	if (!com.stacking_zones)
		return 0;
	while (com.stacking_zones[i].centre.x >= 0.0) i++;
	return i;
}

// copy the image zone into a double buffer, it remains upside-down
static int copy_image_zone_to_fftw(fits *fit, const stacking_zone *zone, fftw_complex *dest,
		int layer) {
	int side = round_to_int(zone->half_side * 2.0);
	// start coordinates on the displayed image, but images are read upside-down
	int startx = round_to_int(zone->centre.x - zone->half_side);
	int starty = round_to_int(zone->centre.y - zone->half_side);
	if (startx < 0 || startx >= fit->rx - side || starty < 0 || starty >= fit->ry - side) {
		/* this zone is partly outside the image, I don't think there's
		 * much we can do for it, it just has to be ignored for this
		 * image for the stacking. */
		return -1;
	}

	WORD *from = fit->pdata[layer] + (fit->ry - starty - side - 1) * fit->rx + startx;
	int stridefrom = fit->rx - side;
	int i, j;

	for (i = 0; i < side; ++i) {
		for (j = 0; j < side; ++j) {
			*dest++ = (double)*from++;
		}
		from += stridefrom;
	}
	return 0;
}

// the same as above with a WORD buffer instead of double
static int copy_image_zone_to_buffer(fits *fit, const stacking_zone *zone, WORD *dest, int layer) {
	int side = round_to_int(zone->half_side * 2.0);
	// start coordinates on the displayed image, but images are read upside-down
	int startx = round_to_int(zone->centre.x - zone->half_side);
	int starty = round_to_int(zone->centre.y - zone->half_side);

	if (startx < 0 || startx >= fit->rx - side || starty < 0 || starty >= fit->ry - side) {
		/* this zone is partly outside the image, I don't think there's
		 * much we can do for it, it just has to be ignored for this
		 * image for the stacking. */
		return -1;
	}

	WORD *from = fit->pdata[layer] + (fit->ry - starty - side - 1) * fit->rx + startx;
	int stride = fit->rx - side;
	int i, j;

	for (i = 0; i < side; ++i) {
		for (j = 0; j < side; ++j) {
			*dest++ = *from++;
		}
		from += stride;
	}
	return 0;
}

/* copy the zone from an image to the same zone shifted in the stacked image */
static void add_image_zone_to_stacking_sum(fits *fit, const stacking_zone *zone, int frame,
		unsigned long *sum[3], int *count[3]) {
	int layer;
	int side = round_to_int(zone->half_side * 2.0);
	int src_startx = round_to_int(zone->centre.x - zone->half_side - zone->mpregparam[frame].x);
	int src_starty = round_to_int(zone->centre.y - zone->half_side + zone->mpregparam[frame].y);
	int dst_startx = round_to_int(zone->centre.x - zone->half_side);
	int dst_starty = round_to_int(zone->centre.y - zone->half_side);

	if (src_startx < 0 || src_startx >= fit->rx - side ||
			src_starty < 0 || src_starty >= fit->ry - side) {
		/* this zone is partly outside the image, we could partially
		 * read it, but for now we just ignore it for stacking */
		return;
	}

	for (layer = 0; layer < fit->naxes[2]; layer++) {
		WORD *from = fit->pdata[layer];
		unsigned long *to = sum[layer];
		int *lcount = count[layer];
		int x, y, stride = fit->rx - side;
		int i = (fit->ry - src_starty - side - 1) * fit->rx + src_startx;
		int o = (fit->ry - dst_starty - side - 1) * fit->rx + dst_startx;

		for (y = 0; y < side; ++y) {
			for (x = 0; x < side; ++x) {
				to[o] += from[i];
				lcount[o]++;
				i++; o++;
			}
			i += stride;
			o += stride;
		}
	}
}

/* copy a buffer representing an area into the sum of pixels being stacked.
 * It's done only for one layer. rx and ry are the dimensions of the stacked
 * image and of sum and count matrices.
 */
static void add_buf_zone_to_stacking_sum(WORD *buf, int layer, const stacking_zone *zone,
		int frame, unsigned long *sum[3], int *count[3], unsigned int rx, unsigned int ry) {
	int side = round_to_int(zone->half_side * 2.0);
	int dst_startx = round_to_int(zone->centre.x - zone->half_side);
	int dst_starty = round_to_int(zone->centre.y - zone->half_side);

	unsigned long *to = sum[layer];
	int *lcount = count[layer];
	int x, y, stride = rx- side;
	int i = 0;
	int o = (ry - dst_starty - side - 1) * rx + dst_startx;

	for (y = 0; y < side; ++y) {
		for (x = 0; x < side; ++x) {
			if (buf[i]) {
				// TODO: is this a good check for empty pixels?
				to[o] += buf[i];
				lcount[o]++;
			}
			i++; o++;
		}
		i += side;
		o += stride;
	}
}

int point_is_inside_zone(int px, int py, stacking_zone *zone) {
	int side = round_to_int(zone->half_side * 2.0);
	int startx = round_to_int(zone->centre.x - zone->half_side);
	int starty = round_to_int(zone->centre.y - zone->half_side);
	return px > startx && px < startx + side && py > starty && py < starty + side;
}

/* we build a list of max_AP of the closest alignment points for a pixel. The
 * list is passed, max_AP too, and we try a new candidate AP for the top max_AP.
 * The built list is not ordered.
 */
static void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP) {
	int i, max = 0, max_idx = 0, found_closer = 0;
	for (i = 0; i < max_AP; i++) {
		double ap_dist = list_for_this_point[i].distance;
		if (ap_dist >= 0) {
			// we search the max to replace it by the new
			if (ap_dist > max) {
				max = ap_dist;
				max_idx = i;
			}
			if (distance < ap_dist)
				found_closer = 1;
		}
		else {
			found_closer = 1;
			max_idx = i;
			break;
		}
	}
	if (found_closer) {
		list_for_this_point[max_idx].distance = distance;
		list_for_this_point[max_idx].zone_index = zone_idx;
	}
}

#if 0
/* when the list of closest AP is complete, we check for special case of pixel position
 * being very close to AP centres. Below this owned_distance threshold, other AP are not
 * considered relevant and are removed of the list. */
static void filter_closest_list_owned(struct weighted_AP *list_for_this_point,
		int max_AP, double owned_distance) {
	int i, unowned_idx = -1, owned = 0;
	for (i = 0; i < max_AP; i++) {
		double ap_dist = list_for_this_point[i].distance;
		if (ap_dist < 0) continue;
		if (ap_dist < owned_distance) {
			owned = 1;
			if (unowned_idx != -1) {
				list_for_this_point[unowned_idx].distance = ap_dist;
				list_for_this_point[unowned_idx].zone_index = list_for_this_point[i].zone_index;
				list_for_this_point[i].distance = -1.0;
				list_for_this_point[i].zone_index = -1;
				while (unowned_idx < i &&
						list_for_this_point[unowned_idx].distance < owned_distance);
				if (i == unowned_idx)
					unowned_idx = -1;
			}
		}
		else {
		       if (owned) {
			       list_for_this_point[i].distance = -1.0;
			       list_for_this_point[i].zone_index = -1;
		       }
		       if (unowned_idx == -1)
			       unowned_idx = i;
		}
	}
}
#endif

/* Sort images indices from the list of best quality for zones *
 * In com.stacking_zones[zone].mpregparam we have the normalized quality for each
 * image for the concerned zone. To speed up the look-up, we create an index here.
 * */
BYTE * sort_zones_quality(int zone_idx, int nb_best, int nb_seq_images) {
	int i;
	BYTE *index = malloc(nb_seq_images);
	// finding the nth element of an unsorted set: sort it, it's easier
	int *indices = apregdata_best(com.stacking_zones[zone_idx].mpregparam, nb_seq_images);

	// output: index with ones if image is among the best
	for (i = 0; i < nb_seq_images; i++)
		index[indices[i]] = i < nb_best;

	free(indices);
	return index;
}

#if defined DEBUG_MPP || defined DEBUG_MPP1
static void save_buffer_tmp(int frame_index, int zone_idx, WORD *buffer, int square_size) {
	char tmpfn[100];	// this is for debug purposes
	sprintf(tmpfn, "/tmp/zone_%d_image_%d.fit", zone_idx, frame_index);
	fits *tmp = NULL;
	new_fit_image(&tmp, square_size, square_size, 1);
	tmp->data = buffer;
	tmp->pdata[0] = tmp->data;
	tmp->pdata[1] = tmp->data;
	tmp->pdata[2] = tmp->data;
	savefits(tmpfn, tmp);
	tmp->data = NULL; // don't free the original buffer
	clearfits(tmp);
}
#endif
