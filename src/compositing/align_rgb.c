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
 * entry in the RGB image popup and by the rgbalign command. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/gui_iface.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "registration/registration.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "compositing/align_rgb.h"

#define REGLAYER 0

/* the largest selection accepted by the one-star method: it must contain a
 * single star for the PSF fit to be meaningful */
#define MAX_PSF_SELECTION 300

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

const char *rgb_align_method_name(rgb_align_method m) {
	switch (m) {
		case RGBALIGN_PSF:
			return _("one star");
		case RGBALIGN_DFT:
			return _("image pattern");
		case RGBALIGN_KOMBAT:
			return _("KOMBAT");
		default:
			return _("global star");
	}
}

/* Checks that the loaded image and the current selection are suitable for the
 * requested alignment method, logging the reason if they are not. Shared by the
 * RGB align menu entries and the rgbalign command. */
gboolean rgb_align_prerequisites_met(rgb_align_method m) {
	if (!single_image_is_loaded() || !isrgb(gfit)) {
		siril_log_message(_("RGB alignment requires a loaded colour image.\n"));
		return FALSE;
	}
	/* the global star method registers on the whole image, the others all
	 * need an area to work on */
	if (m != RGBALIGN_GLOBAL && (com.selection.w <= 0 || com.selection.h <= 0)) {
		siril_log_message(_("The %s alignment method requires a selection. Make one in the "
					"image or use the boxselect command.\n"), rgb_align_method_name(m));
		return FALSE;
	}
	if (m == RGBALIGN_PSF && (com.selection.w > MAX_PSF_SELECTION || com.selection.h > MAX_PSF_SELECTION)) {
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return FALSE;
	}
	return TRUE;
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
	/* Channels the registration failed on are excluded from the sequence and
	 * so keep their original data: the result is only partially aligned, which
	 * is easy to miss in a script if we don't say so. */
	for (int i = 0; i < 3; i++) {
		if (!seq->imgparam[i].incl)
			siril_log_warning(_("The %s channel could not be aligned and is left unchanged.\n"),
					channel_number_to_name(i));
	}
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

int rgb_align(rgb_align_method m) {
	struct registration_args regargs = { 0 };
	struct registration_method *method;
	framing_type framing = FRAMING_COG;
	int retval1 = 0, retval2 = 0;

	initialize_methods();
	if (initialize_internal_rgb_sequence()) {
		siril_log_message(_("Could not extract the channels of the loaded image.\n"));
		return 1;
	}
	gui_iface.set_busy(TRUE);
	gui_iface.set_progress(PROGRESS_RESET, NULL);

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

	if (!reserve_thread()) {
		siril_log_message(_("A processing operation is already running.\n"));
		gui_iface.set_busy(FALSE);
		return 1;
	}
	retval1 = method->method_ptr(&regargs);
	free(regargs.imgparam);
	regargs.imgparam = NULL;
	free(regargs.regparam);
	regargs.regparam = NULL;
	if (retval1) {
		gui_iface.set_progress(PROGRESS_DONE, _("Error in channels alignment."));
		gui_iface.set_busy(FALSE);
		unreserve_thread();
		return retval1;
	}
	retval2 = register_apply_reg(&regargs);
	compose();
	unreserve_thread();

	if (retval2) {
		gui_iface.set_progress(PROGRESS_DONE, _("Error in layers alignment."));
	} else {
		gui_iface.set_progress(PROGRESS_DONE, _("Registration complete."));
		notify_gfit_data_modified();
		gfit_modified_update_gui();
		siril_log_message(_("Aligned RGB channels using the %s method\n"), rgb_align_method_name(m));
	}
	gui_iface.set_busy(FALSE);
	free_internal_sequence(seq);
	seq =  NULL;
	return retval2;
}
