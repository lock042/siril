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

/* stacking process:
 * 1. compute parameters:
 *  - zone size: can it be guessed from the quality values? Are they a good
 *  indicator of SNR? Otherwise, we may count the number of pixels above a
 *  threshold to estimate the size of the planet in pixels and the samping of the
 *  acquisition. Zone size should depend on sampling, SNR and sky quality.
 *  - overlap, the number of pixels that are common between two zones and of the
 *  number of pixel in the zones, as a percentage. It depends on how much the
 *  images move, how much spaced the zones are (which depends on their size)
 *  because the more zones are far apart, the more shear might appear if they
 *  don't align. The value should not be higher than 50%, which would mean that
 *  the centre of two zones are both in the two zones. It's probably wiser to have
 *  the centres only in one zone.
 *  - extra read for zones, when trying to align zones, even if we read them in
 *  images taking into account the global shifts from the first registration, they
 *  might not perfectly fit. This value indicates how many pixels extra half side
 *  we need to read for each zone. This depends probably on the same data as the
 *  overlap, with emphasis on the standard deviation of the set of shifts.
 *  - sigma low and high of the stacking may be guesses from quality standard
 *  deviation?
 *
 *  Then show these values in the GUI and wait for user input. When the user
 *  changes the zone size, it may update the overlap amount in automatic mode.
 *
 * 2. place the stacking zones
 *  for the automatic mode, we will have to run a filter on the reference image to
 *  find points of interest. The points should not be on the planet's edge, but
 *  since it's a high contrast zone, it needs to be filtered out. The points
 *  should be spaced in a way that somewhat respect the overlap parameter
 *
 * 3. do the mpp stacking
 *  - each zone has to be read, with the extra read borders, on each image of the
 *  X best percent selected by user, to run a DFT registration on them. This gives
 *  local shifts, or zone shifts. A texture of the size of the reference image is
 *  used to store the cumulated number of zones that contain each pixel. This
 *  helps managing the variable pixel count in the stack.
 *  - average stacking with rejection is then run almost normally, with the blocks
 *  computed for multi-threading almost as usual, but taking into account the
 *  overlap that will add some pixels to each pixel stacks.
 *    For each pixel, we get for each image the pixels of each zone slided with
 *    the local shift (there may be several zones that contain the pixel in
 *    overlapping zones) or globally shifted pixel if the pixel is not in a zone.
 *
 * There may be a way do to the zones registration and the stacking at the same
 * time, but this gets complicated with the stacking blocks for the
 * multi-threading. Separating in two as explained above should not give an
 * unreasonable speed penalty.
 *    
 */

/* TODO:
 * extra read for zones for the registration, does it make sens DFT-wise?
 * is it possible to have a confidence rating for a DFT match? To discard
 *   failures and adjust barycentre weights
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
 * make a background thread that precomputes FFTW wisdom, it's too slow for big zones
 * restrict zone size (see FFTW comment below) to some values
 * restrict zones inside the image
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
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"
#include "gui/planetary_callbacks.h"

#define FFTW_WISDOM_FILE "fftw_wisdom"

static fits refimage;
static char *refimage_filename;

static gpointer sequence_analysis_thread_func(gpointer p);
static gboolean end_reference_image_stacking(gpointer p);
static int the_multipoint_registration(struct mpr_args *args);
static int the_multipoint_barycentric_sum_stacking(struct mpr_args *args);
static int get_number_of_zones();
static int copy_image_zone_to_fftw(fits *fit, stacking_zone *zone, fftw_complex *dest,
		int layer, double shiftx, double shifty);
static void compute_zones_confidence(struct mpregdata *regparam, int nb_zones);

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
	if (stack_args->retval)
		return GINT_TO_POINTER(-1);

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

gpointer the_multipoint_processing(gpointer ptr) {
	struct mpr_args *args = (struct mpr_args*)ptr;
	int retval;
	/* multi-point registration: compute the local shifts */
	retval = the_multipoint_registration(args);
	if (!retval) {
		/* multi-point stacking: stack with local shifts */
		retval = the_multipoint_barycentric_sum_stacking(args);

		/* registering the generic stacking idle with the required args */
		struct stacking_args *stackargs = malloc(sizeof(struct stacking_args));;
		stackargs->retval = retval;
		stackargs->output_filename = args->output_filename;
		stackargs->output_overwrite = args->output_overwrite;
		stackargs->seq = args->seq;
		// TODO: nb_images_to_stack for summary output
		siril_add_idle(end_stacking, stackargs);
	}
	return GINT_TO_POINTER(retval);
}	

/* about FFTW: http://www.fftw.org/fftw3_doc/Introduction.html
 * The standard FFTW distribution works most efficiently for arrays whose size
 * can be factored into small primes (2, 3, 5, and 7), and otherwise it uses a
 * slower general-purpose routine.
 * Maybe we need to constrain the zone size to these numbers, that's not too
 * many that do not register (11 13 17 19 23 31 37 41 43 47 53 59 61 67 71 ...)
 *
 * Applying the phase correlation method to a pair of images
 * produces a third image which contains a single peak. The
 * location of this peak corresponds to the relative translation
 * between the images.
 * The Fourier-Mellin transform extends phase correlation to
 * handle images transformed by both translation and rotation.
 */
static int the_multipoint_registration(struct mpr_args *args) {
	int zone_idx, nb_zones, frame;
	int abort = 0;
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

		copy_image_zone_to_fftw(&refimage, zone, ref[zone_idx],
				args->layer, 0.0, 0.0);

		fftw_execute_dft(fplan[zone_idx], ref[zone_idx], in[zone_idx]);

	}
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
	args->regdata = malloc(args->seq->number * sizeof(struct mpregdata*));
	if (!args->regdata) {
		fprintf(stderr, "Stacking: memory allocation failure for registration data\n");
		return -1;
	}

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
		if (!args->filtering_criterion(args->seq, args->layer, frame, args->filtering_parameter))
			continue;

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "aligning zones for image %d\n", frame);

		args->regdata[frame] = malloc(nb_zones * sizeof(struct mpregdata));
		if (!args->regdata[frame]) {
			fprintf(stderr, "Stacking: memory allocation failure for registration data per frame\n");
			return -1;
		}

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
				{ .x = zone->centre.x + regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };

			if (!zones[zone_idx]) {	// keep across images
			       zones[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       out2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       convol2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			}

			copy_image_zone_to_fftw(&fit, &shifted_zone, zones[zone_idx],
					args->layer, 0.0, 0.0);

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
			double peak;
			for (i = 1; i < nb_pixels; ++i) {
				if (creal(out2[zone_idx][i]) > creal(out2[zone_idx][shift])) {
					shift = i;
				}
			}
			peak = creal(out2[zone_idx][shift]);
			int shifty = shift / side;
			int shiftx = shift % side;
			if (shifty > zone->half_side) {
				shifty -= side;
			}
			if (shiftx > zone->half_side) {
				shiftx -= side;
			}

			/* for Y, it's a bit special because FITS are upside-down */
			int sign = args->seq->type == SEQ_SER ? -1 : 1;

			/* shitfs for this image and this zone is shiftx, sign*shifty */
			fprintf(stdout, "frame %d, zone %d shifts: %d,%d,\tpeak: %f\n", frame, zone_idx, shiftx, shifty*sign, peak);
			args->regdata[frame][zone_idx].x = (double)shiftx;
			args->regdata[frame][zone_idx].y = (double)shifty * sign;
			args->regdata[frame][zone_idx].peak = peak;
		}

		compute_zones_confidence(args->regdata[frame], nb_zones);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

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

	return 0;
}


struct weighted_AP {
	int zone_index;
	double distance;
	float confidence;
};

static void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP);
//static void filter_closest_list_owned(struct weighted_AP *list_for_this_point,
//		int max_AP, double owned_distance);

/* this function is a stacking sum that takes shifts from the multipoint registration
 * instead of taking them from the global registration. If it works well, it will have to
 * be exploded as a generic processing function. */
static int the_multipoint_barycentric_sum_stacking(struct mpr_args *args) {
	int frame, zone_idx, nb_zones, abort = 0;
	regdata *regparam = args->seq->regparam[args->layer];
	struct weighted_AP *closest_zones_map = malloc(args->seq->rx * args->seq->ry * args->nb_closest_AP * sizeof(struct weighted_AP));
	if (!closest_zones_map) {
		fprintf(stderr, "Stacking: memory allocation failure for zone mapping\n");
		return -1;
	}
	nb_zones = get_number_of_zones();

	/* precompute the closest zones from each pixel */
	int x, y;
	for (y = 0; y < args->seq->ry; y++) {
		for (x = 0; x < args->seq->rx; x++) {
			int ap;
			struct weighted_AP *list_for_this_pixel = closest_zones_map + (x + y * args->seq->ry) * args->nb_closest_AP;
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
	args->sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (!args->sum[0]){
		fprintf(stderr, "Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		args->sum[1] = args->sum[0] + nbdata;	// index of green layer in sum[0]
		args->sum[2] = args->sum[0] + nbdata*2;	// index of blue layer in sum[0]
	} else {
		args->sum[1] = NULL;
		args->sum[2] = NULL;
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
		 * pixel. THe direction of each AP would also be important, the
		 * distribution of the AP having a role in the weighing. */
		int x, y;
		int pixel = 0;	// index in sum[0]
		for (y = 0; y < fit.ry; y++) {
			for (x = 0; x < fit.rx; x++) {
				struct weighted_AP *list_for_this_AP = closest_zones_map + (x + y * args->seq->ry) * args->nb_closest_AP;
				int i;
				double weight = 0.0, sumx = 0.0, sumy = 0.0, shiftx, shifty;

				for (i = 0; i < args->nb_closest_AP; i++) {
					double this_weight;
					if (list_for_this_AP[i].distance < 0.0)
						continue;
					if (list_for_this_AP[i].distance <= 1.0)
						this_weight = args->max_distance;
					else this_weight = args->max_distance / list_for_this_AP[i].distance;
					this_weight *= args->regdata[frame][list_for_this_AP[i].zone_index].peak;
					weight += this_weight;
					sumx += args->regdata[frame][list_for_this_AP[i].zone_index].x * this_weight;
					sumy += args->regdata[frame][list_for_this_AP[i].zone_index].y * this_weight;
				}
				// if zones are too far away, which will happen for pixels
				// away from the planet, just take the global shift
				if (weight == 0.0) {
					shiftx = regparam[frame].shiftx;
					shifty = regparam[frame].shifty;
				} else {
					shiftx = sumx / weight;
					shifty = sumy / weight;
				}

				// then do the regular sum stacking with shift
				int nx = round_to_int(x - shiftx);
				int ny = round_to_int(y - shifty);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					int ii = ny * fit.rx + nx;		// index in source image
					int layer;
					for (layer=0; layer<args->seq->nb_layers; ++layer) {
#ifdef _OPENMP
#pragma omp atomic
#endif
						args->sum[layer][pixel] += fit.pdata[layer][ii];
					}
				}
				++pixel;
			}
		}

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

	fprintf(stdout, "barycentre stacking ended, creating final image\n");

	if (!abort) {
		/* copying result into gfit, code copied from sum_stacking_finalize_hook() */
		nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
		int i, layer;
		unsigned long max = 0L;	// max value of the image's channels
		// find the max first
#ifdef _OPENMP
#pragma omp parallel for reduction(max:max)
#endif
		for (i=0; i < nbdata; ++i)
			if (args->sum[0][i] > max)
				max = args->sum[0][i];

		clearfits(&gfit);
		fits *fit = &gfit;
		if (new_fit_image(&fit, args->seq->rx, args->seq->ry, args->seq->nb_layers))
			return -1;
		gfit.hi = round_to_WORD(max);

		double ratio = 1.0;
		if (max > USHRT_MAX)
			ratio = USHRT_MAX_DOUBLE / (double)max;

		nbdata = args->seq->ry * args->seq->rx;
		for (layer=0; layer<args->seq->nb_layers; ++layer){
			unsigned long* from = args->sum[layer];
			WORD *to = gfit.pdata[layer];
			for (i=0; i < nbdata; ++i) {
				if (ratio == 1.0)
					*to++ = round_to_WORD(*from++);
				else	*to++ = round_to_WORD((double)(*from++) * ratio);
			}
		}
	}

	free(args->sum[0]);

	return 0;
}

static int get_number_of_zones() {
	int i = 0;
	while (com.stacking_zones[i].centre.x >= 0.0) i++;
	return i;
}

static int copy_image_zone_to_fftw(fits *fit, stacking_zone *zone, fftw_complex *dest,
		int layer, double shiftx, double shifty) {

	int side = round_to_int(zone->half_side * 2.0);
	int startx = round_to_int(zone->centre.x - zone->half_side + shiftx);
	int starty = round_to_int(zone->centre.y - zone->half_side + shiftx);
	if (startx < 0 || startx >= fit->rx - side || starty < 0 || starty >= fit->ry - side) {
		/* this zone is partly outside the image, I don't think there's
		 * much we can do for it, it just has to be ignored for this
		 * image for the stacking. */
		return -1;
	}

	WORD *from = fit->pdata[layer] + (fit->ry - starty) * fit->rx + starty;
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
			if (ap_dist > max) {
				max = ap_dist;
				max_idx = i;
			}
			if (distance > ap_dist)
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
 * being very closed to AP centres. Below this owned_distance threshold, other AP are not
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

/* normalize the peaks to have a confidence ratio [0; 1].
 * Typical values range from 0.01 to 1:
 * 0.7 to 1 for good confidence, 0.2 to 0.7 for medium, 0.01 to 0.2 for poor.
 * Would have been better to do it on all images at the same time, but the data structure
 * does not currently permit it.
 */
static void compute_zones_confidence(struct mpregdata *regparam, int nb_zones) {
	int zone_idx;
	double max = 0;
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		if (regparam[zone_idx].peak > max)
			max = regparam[zone_idx].peak;
	}

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++)
		regparam[zone_idx].peak /= max;
}
