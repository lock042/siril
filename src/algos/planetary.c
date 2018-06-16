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

/* TODO:
 * allow the ref image not to be displayed when loading a sequence if we
 *   display the planetary reference image 
 * demosaicing setting: manage prepro in on_prepro_button_clicked
 * fix the new arg->layer and filters in processing and elsewhere
 * GUI mode switch: maybe not
 * remove menuitemcolor and menuitemgray in planetary mode
 * chain hi/lo/mode when RGB vport is modified
 * add a 'clear zones' button to remove all existing
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

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

static fits refimage;
static char *refimage_filename;

static gpointer sequence_analysis_thread_func(gpointer p);
static gboolean end_reference_image_stacking(gpointer p);
static int get_number_of_zones();
static int copy_image_zone_to_fftw(fits *fit, stacking_zone *zone, fftw_complex *dest,
		int layer, double shiftx, double shifty);

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
	return refimage.naxes[0] > 0;
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
				       " for layer %d, using it"), layer);
		refimage_filename = ref_name;
		display_refimage_if_needed();
	} else {
		free(ref_name);
	}
}

gpointer the_multipoint_registration(gpointer ptr) {
	struct mpr_args *args = (struct mpr_args*)ptr;
	int zone_idx, nb_zones, frame;
	int abort = 0;
	stacking_zone *zone;
	fftw_complex **ref, **in, **out, **convol;
	fftw_plan *p, *q;
	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return GINT_TO_POINTER(-1);

	nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		// TODO: add the idle instead
		return GINT_TO_POINTER(-1);
	}
	ref = malloc(nb_zones * sizeof(fftw_complex*));
	in = malloc(nb_zones * sizeof(fftw_complex*));
	out = malloc(nb_zones * sizeof(fftw_complex*));
	convol = malloc(nb_zones * sizeof(fftw_complex*));
	p = malloc(nb_zones * sizeof(fftw_plan));
	q = malloc(nb_zones * sizeof(fftw_plan));
	
	/* reading zones in the reference image */
	fprintf(stdout, "loading reference zones\n");
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		zone = &com.stacking_zones[zone_idx];
		int side = round_to_int(zone->half_side * 2.0);
		int nb_pixels = side * side;
		ref[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		in[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		out[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		convol[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);

		/* TODO: what's that?
		if (nb_frames > 200.f)
			plan = FFTW_MEASURE;
		else plan = FFTW_ESTIMATE; */
		int plan = FFTW_MEASURE;

		// CRASH HERE v
		p[zone_idx] = fftw_plan_dft_2d(side, side, ref[zone_idx], out[zone_idx], FFTW_FORWARD, plan);
		q[zone_idx] = fftw_plan_dft_2d(side, side, convol[zone_idx], out[zone_idx], FFTW_BACKWARD, plan);

		copy_image_zone_to_fftw(&refimage, zone, ref[zone_idx],
				args->layer, 0.0, 0.0);

		fftw_execute_dft(p[zone_idx], ref[zone_idx], in[zone_idx]);
	}

	fftw_complex **zones = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **out2 = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **convol2 = calloc(nb_zones, sizeof(fftw_complex*));

	/* for each image, we read the zones with global shift and register them */
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

		/* reading zones in the reference image */
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

			fftw_execute_dft(p[zone_idx], zones[zone_idx], out2[zone_idx]);

			for (i = 0; i < nb_pixels; i++) {
				convol2[zone_idx][i] = in[zone_idx][i] * conj(out2[zone_idx][i]);
			}
			fftw_execute_dft(q[zone_idx], convol2[zone_idx], out2[zone_idx]);

			int shift = 0;
			for (i = 1; i < nb_pixels; ++i) {
				if (creal(out2[zone_idx][i]) > creal(out2[zone_idx][shift])) {
					shift = i;
					// TODO: break or get last or max value?
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

			/* for Y, it's a bit special because FITS are upside-down */
			int sign = args->seq->type == SEQ_SER ? -1 : 1;

			/* shitfs for this image and this zone is shiftx, sign*shifty */
			fprintf(stdout, "frame %d, zone %d shifts: %d,%d\n", frame, zone_idx, shiftx, shifty*sign);
		}
	}

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		fftw_free(zones[zone_idx]);
		fftw_free(out2[zone_idx]);
		fftw_free(convol2[zone_idx]);

		fftw_destroy_plan(p[zone_idx]);
		fftw_destroy_plan(q[zone_idx]);
		fftw_free(in[zone_idx]);
		fftw_free(out[zone_idx]);
		fftw_free(ref[zone_idx]);
		fftw_free(convol[zone_idx]);
	}
	free(zones); free(out2); free(convol2);
	free(p); free(q);
	free(in); free(out); free(ref); free(convol);

	//siril_add_idle();
	return NULL;
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

