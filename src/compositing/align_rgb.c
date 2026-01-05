/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

/* This file is currently not used by compositing, only by the RGB align menu
 * entry in the RGB image popup. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "registration/registration.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#define REGLAYER 0

static sequence *seq = NULL;		// the sequence of channels
static struct registration_method *reg_methods[5];


static void initialize_methods() {
	reg_methods[0] = new_reg_method(_("One star registration (deep-sky)"),
			&register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[1] = new_reg_method(_("Image pattern alignment (planetary/deep-sky)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[2] = new_reg_method(_("Global star registration (deep-sky)"),
			&register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[3] = new_reg_method(_("KOMBAT registration (planetary / deep-sky)"),
			&register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[4] = NULL;
}

// We cannot currently do this in free_sequence() because compositing still
// uses the references, so we have to do it here as a special case
static void free_internal_sequence(sequence *seq) {
	if (seq) {
		for (int i = 0; i < seq->number; i++)
			clearfits(internal_sequence_get(seq, i));
		free_sequence(seq, TRUE);
	}
	clear_stars_list(TRUE);
}

static int initialize_internal_rgb_sequence() {
	if (seq) free_internal_sequence(seq);

	seq = create_internal_sequence(3);
	for (int i = 0; i < 3; i++) {
		fits *fit = calloc(1, sizeof(fits));
		if (extract_fits(gfit, fit, i, FALSE)) {
			free(fit);
			free_sequence(seq, TRUE);
			return -1;
		}
		internal_sequence_set(seq, i, fit);
	}
	seq->rx = gfit->rx;
	seq->ry = gfit->ry;
	seq->bitpix = gfit->bitpix;

	return 0;
}

static void compose() {
	size_t npixels = gfit->rx * gfit->ry;
	fits *fit[3];
	for (int i = 0 ; i < 3 ; i++) {
		fit[i] = internal_sequence_get(seq, i);
	}
	if (gfit->type == DATA_FLOAT) {
		for (int i = 0 ; i < 3 ; i++) {
			memcpy(gfit->fpdata[i], fit[i]->fdata, sizeof(float) * npixels);
		}
	} else {
		for (int i = 0 ; i < 3 ; i++) {
			memcpy(gfit->pdata[i], fit[i]->data, sizeof(WORD) * npixels);
		}
	}
}

int rgb_align(int m) {
	struct registration_args regargs = { 0 };
	struct registration_method *method;
	framing_type framing = FRAMING_COG;
	int retval1 = 0, retval2 = 0;

	initialize_methods();
	initialize_internal_rgb_sequence();
	set_cursor_waiting(TRUE);
	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* align it */
	method = reg_methods[m];
	regargs.seq = seq;
	regargs.no_output = FALSE;
	get_the_registration_area(&regargs, method);
	regargs.layer = REGLAYER;
	seq->reference_image = 0;
	regargs.seq->nb_layers = 1;
	regargs.max_stars_candidates = MAX_STARS_FITTED;
	regargs.run_in_thread = FALSE;
	regargs.interpolation = OPENCV_LANCZOS4;
	regargs.clamp = TRUE;
	regargs.framing = framing;
	regargs.output_scale = 1.f;
	regargs.percent_moved = 0.50f; // Only needed for KOMBAT
	regargs.two_pass = TRUE;
	if (method->method_ptr == register_shift_fwhm || method->method_ptr == register_shift_dft)
		regargs.type = SHIFT_TRANSFORMATION;
	else
		regargs.type = HOMOGRAPHY_TRANSFORMATION;
	com.run_thread = TRUE;	// fix for the canceling check in processing

	retval1 = method->method_ptr(&regargs);
	free(regargs.imgparam);
	regargs.imgparam = NULL;
	free(regargs.regparam);
	regargs.regparam = NULL;
	if (retval1) {
		set_progress_bar_data(_("Error in channels alignment."), PROGRESS_DONE);
		set_cursor_waiting(FALSE);
		com.run_thread = FALSE;	// fix for the cancelling check in processing
		return retval1;
	}
	retval2 = register_apply_reg(&regargs);
	compose(); // Register_apply_reg has already done the alignment for us

	com.run_thread = FALSE;	// fix for the canceling check in processing

	if (retval2) {
		set_progress_bar_data(_("Error in layers alignment."), PROGRESS_DONE);
	} else {
		set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
		notify_gfit_modified();
		redraw(REMAP_ALL);
	}
	siril_log_message(_("Aligned RGB channels\n"));
	set_cursor_waiting(FALSE);
	free_internal_sequence(seq);
	seq =  NULL;
	return retval2;
}
