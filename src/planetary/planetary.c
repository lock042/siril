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
 * This is defined in the_old_local_multipoint_sum_stacking().
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
#include "laplacian_quality.h"
#include "quality.h"
#include "registration.h"
#include "stacking.h"
#include "gui.h"
#include "core/siril.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/sequence_filtering.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/plot.h"
#include "gui/zones.h"
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"

#define FFTW_WISDOM_FILE "fftw_wisdom"
//#define USE_SEQUENCE_REF_IMAGE

static fits refimage;
static char *refimage_filename;

static gboolean end_reference_image_stacking(gpointer p);

gpointer sequence_analysis_thread_func(gpointer p) {
	struct registration_args *reg_args = (struct registration_args *) p;

	/* registration, as defined in on_planetary_analysis_clicked() */
	reg_args->retval = reg_args->func(reg_args);

	if (reg_args->retval)
		return GINT_TO_POINTER(-1);

	/* generating the reference image by taking a few of the best images */
	struct stacking_args *stack_args = calloc(1, sizeof(struct stacking_args));
	stack_args->method = stack_summing_generic;
	stack_args->seq = reg_args->seq;
	stack_args->reglayer = reg_args->layer;
	stack_args->filtering_criterion = seq_filter_quality;
	stack_args->filtering_parameter = 0.75;	// not the right way do to it, sorting is
	stack_args->nb_images_to_stack = compute_nb_filtered_images(reg_args->seq,
			seq_filter_quality, 0.75, reg_args->layer);
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
		retval = the_laplace_multipoint_quality_analysis(args);
		if (retval) return GINT_TO_POINTER(retval);
	}

	args->refimage = &refimage;

	/* multi-point registration: compute the local shifts */
	//retval = the_multipoint_dft_registration(args);
	//retval = the_multipoint_ecc_registration(args);
	retval = the_multipoint_steepest_descent_registration_and_stacking(args);
	if (retval) return GINT_TO_POINTER(retval);

	/* reference image stacking: stack global images with local shifts.
	 * This image is used for areas of the result where there is no zone defined */
	//retval = the_global_multipoint_barycentric_sum_stacking(args);
	//if (retval) return GINT_TO_POINTER(retval);

	/* multi-point stacking: stack zones with local shifts and creates the result in gfit */
	//retval = the_old_local_multipoint_sum_stacking(args);

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

int copy_image_zone_to_buffer(fits *fit, const stacking_zone *zone, WORD *dest, int layer) {
	int side = get_side(zone);
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

// for debug purposes
void save_buffer_tmp(int frame_index, int zone_idx, WORD *buffer, int square_size) {
	char tmpfn[100];
	sprintf(tmpfn, "/tmp/zone_%d_image_%d.fit", zone_idx, frame_index);
	fits *tmp = NULL;
	new_fit_image(&tmp, square_size, square_size, 1, buffer);
	savefits(tmpfn, tmp);
	tmp->data = NULL; // don't free the original buffer
	clearfits(tmp);
}

