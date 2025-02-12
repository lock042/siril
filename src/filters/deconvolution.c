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
#include <fftw3.h>
#include <math.h>
#include <locale.h>
#include <gdk/gdk.h>
#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_date.h"
#include "core/command.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/sequence_list.h"
#include "gui/message_dialog.h"
#include "gui/PSF_list.h"
#include "gui/registration_preview.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/progress_and_log.h"
#include "gui/newdeconv.h"
#include "filters/deconvolution/deconvolution.h"
#include "filters/synthstar.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "io/sequence.h"
#include "io/ser.h"

gboolean aperture_warning_given = FALSE;
gboolean bad_load = FALSE;
orientation_t imageorientation;
gboolean next_psf_is_previous = FALSE;
fits* the_fit = NULL;
int sequence_is_running;

estk_data args = { 0 };

orientation_t get_imageorientation() {
	orientation_t result;
	if (sequence_is_loaded() && com.seq.type == SEQ_SER) {
		result = TOP_DOWN;
	} else if (sequence_is_loaded()) {
		result = BOTTOM_UP; // All other sequences should be BOTTOM_UP
	} else if (single_image_is_loaded()) {
		if (!g_strcmp0(the_fit->keywords.row_order, "TOP-DOWN")) {
			result = TOP_DOWN;
		} else if (!g_strcmp0(the_fit->keywords.row_order, "BOTTOM-UP")) {
			result = BOTTOM_UP;
	} else {
			result = UNDEFINED;
		}
	} else {
		result = UNDEFINED;
	}
	return result;
}

void reset_conv_args(estk_data* args) {
	siril_debug_print("Resetting deconvolution args\n");

	// Basic image and kernel parameters
	args->savepsf_filename = NULL;
	args->save_after = FALSE;
	args->stars_need_clearing = FALSE;
	args->recalc_ks = FALSE;
	args->psftype = PSF_BLIND;
	the_fit = (!com.headless && gui.roi.active) ? &gui.roi.fit : &gfit;
	imageorientation = get_imageorientation();
	args->fdata = NULL;
	args->rx = 0;
	args->ry = 0;
	args->ks = 15;
	if (!com.kernel)
		args->kchans = 1;

	// Process parameters
	args->blindtype = BLIND_L0;

	// l0 Descent Kernel Estimation parameters
	args->lambda = 1.f / 3000.f;
	args->lambda_ratio = 1/1.1f;
	args->lambda_min = 1e-3f;
	args->gamma = 20.f;
	args->iterations = 3;
	args->multiscale = FALSE;
	args->scalefactor = 0.5f;
	args->kernel_threshold_max = 0.f;
	args->remove_isolated = FALSE;
	args->better_kernel = TRUE;
	args->upscaleblur = 0.f;
	args->downscaleblur = 1.6f;
	args->k_l1 = 0.5f;

	// Spectral irregularity kernel estimation parameters
	args->ninner = 300;
	args->ntries = 30;
	args->nouter = 3;
	args->compensationfactor = 2.1f;
	args->medianfilter = 1.f; // check
	args->intermediatedeconvolutionweight = 3000.f;
	args->finaldeconvolutionweight = 3000.f;

	// Synthetic kernel parameters
	args->profile = PROFILE_GAUSSIAN;
	args->psf_fwhm = 3.0f;
	args->psf_beta = 4.5f;
	args->psf_angle = 0.f;
	args->psf_ratio = 1.f;
	args->airy_wl = 525.f;
	args->airy_pixelsize = 0.1f;
	args->airy_fl = 1.f;
	args->airy_diameter = 1.f;
	args->airy_obstruction = 0.f;

	// Non-blind deconvolution parameters
	args->nonblindtype = DECONV_RL;
	args->finaliters = 10;
	args->alpha = 1.f / 3000.f;
	args->stopcriterion = 0.002f;
	args->rl_method = RL_GD;
	args->stepsize = 0.0003f;
	args->regtype = REG_TV_GRAD;
}

void reset_conv_kernel() {
	if (com.kernel != NULL) {
		free(com.kernel);
		com.kernel = NULL;
		com.kernelchannels = 1;
		com.kernelsize = 0;
	}
}

void reset_conv_controls_and_args() {
	if (!get_thread_run() || (the_fit == NULL))
		reset_conv_args(&args);
	gui_function(reset_conv_controls, NULL);
}

void check_orientation() {
	int npixels = com.kernelsize * com.kernelsize;
	int ndata = npixels * com.kernelchannels;
	if (get_imageorientation() != args.kernelorientation) {
		float *flip_the_kernel = (float*) malloc(ndata * sizeof(float));
		for (int c = 0 ; c < com.kernelchannels; c++) {
			for (int i = 0 ; i < com.kernelsize ; i++) {
				for (int j = 0 ; j < com.kernelsize ; j++) {
					flip_the_kernel[i + (com.kernelsize - j - 1) * com.kernelsize + npixels * c] = com.kernel[i + com.kernelsize * j + npixels *c];
				}
			}
		}
		free(com.kernel);
		com.kernel = flip_the_kernel;
		args.kernelorientation = (args.kernelorientation == BOTTOM_UP) ? TOP_DOWN : BOTTOM_UP;
	}
}

int load_kernel(gchar* filename) {
	int retval = 0;
	int orig_size;
	args.kernelorientation = BOTTOM_UP; // PSFs are always BOTTOM_UP when saved
	gboolean original_debayer_setting = com.pref.debayer.open_debayer;
	com.pref.debayer.open_debayer = FALSE;
	bad_load = FALSE;
	args.nchans = the_fit->naxes[2];
	fits load_fit = { 0 };
	if ((retval = read_single_image(filename, &load_fit, NULL, FALSE, NULL, FALSE, TRUE))) {
		bad_load = TRUE;
		goto ENDSAVE;
	}
	if (load_fit.rx != load_fit.ry){
		retval = 1;
		char *msg = siril_log_color_message(_("Error: PSF file does not contain a square PSF. Cannot load this file.\n"), "red");
		// no need to check com.script as it is built into the siril_message_dialog function
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong PSF size"), msg);
		bad_load = TRUE;
		goto ENDSAVE;
	}
	if (com.kernel)
		reset_conv_kernel();
	if (!(load_fit.rx % 2)) {
		com.kernelsize = load_fit.rx - 1;
		orig_size = load_fit.rx;
		siril_log_color_message(_("Warning: PSF file is even (%d x %d). PSFs should always be odd. Cropping by 1 pixel in each direction. "
				"This may not produce optimum results.\n"), "salmon", load_fit.rx, load_fit.rx);
	} else {
		com.kernelsize = load_fit.rx;
		orig_size = com.kernelsize;
	}
	com.kernelchannels = load_fit.naxes[2];
	args.kchans = com.kernelchannels;
	gui_function(set_kernel_size_in_gui, NULL);

	int npixels = com.kernelsize * com.kernelsize;
	int orig_pixels = orig_size * orig_size;
	int ndata = npixels * com.kernelchannels;
	com.kernel = (float*) malloc(ndata * sizeof(float));
	if (load_fit.type == DATA_FLOAT) {
		for (int c = 0 ; c < com.kernelchannels ; c++) {
			for (int i = 0 ; i < com.kernelsize ; i++) {
				for (int j = 0 ; j < com.kernelsize ; j++) {
				com.kernel[i + j * com.kernelsize + c * npixels] = load_fit.fdata[i + j * orig_size + c * orig_pixels];
				}
			}
		}
	} else {
		for (int c = 0 ; c < com.kernelchannels ; c++) {
			for (int i = 0 ; i < com.kernelsize ; i++) {
				for (int j = 0 ; j < com.kernelsize ; j++) {
					com.kernel[i + j * com.kernelsize + c * npixels] = (float) load_fit.data[i + j * orig_size + c * orig_pixels] / USHRT_MAX_SINGLE;
				}
			}
		}
	}
	// Handle SER orientation issues, if the image orientation is TOP_DOWN we need to match it
	if (get_imageorientation() == TOP_DOWN) {
		float *flip_the_kernel = (float*) malloc(ndata * sizeof(float));
		for (int c = 0 ; c < com.kernelchannels ; c++) {
			for (int i = 0 ; i < com.kernelsize ; i++) {
				for (int j = 0 ; j < com.kernelsize ; j++) {
					flip_the_kernel[i + (com.kernelsize - j - 1) * com.kernelsize + c * npixels] = com.kernel[i + com.kernelsize * j + c * npixels];
				}
			}
		}
		free(com.kernel);
		com.kernel = flip_the_kernel;
		args.kernelorientation = (args.kernelorientation == BOTTOM_UP) ? TOP_DOWN : BOTTOM_UP;
	}
	if (com.kernelchannels > args.nchans) { // If we have a color kernel but the open image is mono, log a warning and desaturate the kernel
		siril_log_message(_("The selected PSF is RGB but the loaded image is monochrome. The PSF will be converted to monochrome (luminance).\n"));
		float* desatkernel = malloc(com.kernelsize * com.kernelsize * sizeof(float));
		for (int i = 0 ; i < com.kernelsize ; i++) {
			for (int j = 0 ; j < com.kernelsize ; j++) {
				float val = 0.f;
				for (int c = 0 ; c < com.kernelchannels ; c++) {
					val += com.kernel[i + j * com.kernelsize + c * npixels];
				}
				desatkernel[i + j * com.kernelsize] = val / com.kernelchannels;
			}
		}
		free(com.kernel);
		com.kernel = desatkernel;
		com.kernelchannels = 1;
		args.kchans = 1;
	}
	com.pref.debayer.open_debayer = original_debayer_setting;
	clearfits(&load_fit);
	ENDSAVE:
	gui_function(DrawPSF, NULL);
	return retval;
}

int save_kernel(gchar* filename) {
	int retval = 0;
	fits *save_fit = NULL;
	//Check there is a PSF to save
	if (com.kernel == NULL) {
		retval = 1;
		siril_log_color_message(_("Error: no PSF has been computed, nothing to save.\n"), "red");
		return retval;
	}
	int npixels = com.kernelsize * com.kernelsize;
	int ndata = npixels * com.kernelchannels;
	// Need to make a sacrificial copy of com.kernel as the save_fit data will be freed when we call clearfits
	float* copy_kernel = malloc(ndata * com.kernelchannels * sizeof(float));
	memcpy(copy_kernel, com.kernel, ndata * com.kernelchannels * sizeof(float));

	// Handle SER orientation issues, if the kernel orientation is TOP_DOWN then we need to reverse it: we always save the
	// kernel in BOTTOM_UP orientation. Note that we don't actually change args->kernelorientation here as we are only flipping
	// the sacrificial copy.
	if (get_imageorientation() == TOP_DOWN) {
		float *flip_the_kernel = (float*) malloc(ndata * sizeof(float));
		for (int c = 0 ; c < com.kernelchannels ; c++) {
			for (int i = 0 ; i < com.kernelsize ; i++) {
				for (int j = 0 ; j < com.kernelsize ; j++) {
					flip_the_kernel[i + (com.kernelsize - j - 1) * com.kernelsize + c * npixels] = copy_kernel[i + com.kernelsize * j + c * npixels];
				}
			}
		}
		free(copy_kernel);
		copy_kernel = flip_the_kernel;
	}

	if ((retval = new_fit_image_with_data(&save_fit, com.kernelsize, com.kernelsize, com.kernelchannels, DATA_FLOAT, copy_kernel))) {
		siril_log_color_message(_("Error preparing PSF for save.\n"), "red");
		return retval;
	}

	if (g_str_has_suffix(filename, ".fit") || g_str_has_suffix(filename, ".fits") || g_str_has_suffix(filename, ".fts")) {
		retval = savefits(filename, save_fit);
	} else {
#ifdef HAVE_LIBTIFF
		retval = savetif(filename, save_fit, 32, "Saved Siril deconvolution PSF", NULL, FALSE, FALSE, TRUE);
#else
		// This needs to catch the case where a colour kernel is loaded, PGM does't support RGB.
		siril_log_color_message(_("This copy of Siril was compiled without libtiff support: saving PSF in FITS format.\n"), "salmon");
		retval = savefits(filename, save_fit);
#endif
	}
	clearfits(save_fit); // also frees copy_kernel
	free(save_fit);
	return retval;
}

int get_kernel() {
	int retval = 0;
	args.rx = the_fit->rx;
	args.ry = the_fit->ry;
	args.nchans = the_fit->naxes[2];
	args.kchans = 1;
	com.kernelchannels = 1;
	com.fftw_max_thread = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;
	switch (args.psftype) {
		case PSF_BLIND: // Blind deconvolution
			reset_conv_kernel();
			switch(args.blindtype) {
				case BLIND_SI:
					com.kernel = gf_estimate_kernel(&args, com.fftw_max_thread);
					break;
				case BLIND_L0:
					com.kernel = estimate_kernel(&args, com.fftw_max_thread);
					break;
			}
			if (get_thread_run())
				siril_log_message(_("Kernel estimation complete.\n"));
			break;
		case PSF_SELECTION: // Kernel from selection
			// Not implemented yet (maybe never)
			break;
		case PSF_STARS: // Kernel from com.stars
			if (com.kernel && sequence_is_running == 1) {
// We don't want to recompute stars on a running sequence, if we are using a kernel from stars
// detected in the first image then the only sensible approach for subsequent images is to reuse
// the kernel, as there is no way of automatically re-detecting the same stars in subsequent images.
				retval = 0;
				goto END;
			}
			if (!(com.stars && com.stars[0])) {
				siril_log_color_message(_("Error: no stars detected. Run findstar or select stars using Dynamic PSF dialog first.\n"), "red");
				retval = 1;
				goto END;
			}
			reset_conv_kernel();
			calculate_parameters();
			com.kernel = (float*) calloc(args.ks * args.ks * args.kchans, sizeof(float));
			if (com.stars[0]->profile == PSF_GAUSSIAN)
				makegaussian(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f,args.psf_ratio, -args.psf_angle);
			else
				makemoffat(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_beta, args.psf_ratio, -args.psf_angle);

			break;
		case PSF_MANUAL: // Kernel from provided parameters
			reset_conv_kernel();
			com.kernel = (float*) calloc(args.ks * args.ks * args.kchans, sizeof(float));
			switch (args.profile) {
				case PROFILE_GAUSSIAN:
					makegaussian(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_ratio, -args.psf_angle);
					break;
				case PROFILE_MOFFAT:
					makemoffat(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_beta, args.psf_ratio, -args.psf_angle);
					break;
				case PROFILE_DISK:
					makedisc(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f);
					break;
				case PROFILE_AIRY:
					if (args.airy_fl == 1.f)
						siril_log_message(_("Warning: focal length appears likely to be incorrect. Continuing anyway...\n"));
					if (args.airy_diameter == 1.f)
						siril_log_message(_("Warning: diameter appears likely to be incorrect. Continuing anyway...\n"));
					if (args.airy_fl == 0.1f)
						siril_log_message(_("Warning: sensor pixel size appears likely to be incorrect. Continuing anyway...\n"));

					makeairy(com.kernel, args.ks, 1.f, +0.5, -0.5, args.airy_wl, args.airy_diameter, args.airy_fl, args.airy_pixelsize, args.airy_obstruction);
					break;
			}
			break;
		case PSF_PREVIOUS:
			if (com.kernel == NULL) {
				siril_log_color_message(_("Error: no previous PSF found. A blind PSF estimator or calculation of PSF from stars or manual parameters must already have been done to use the Previous PSF option.\n"), "red");
				retval = 1;
				goto END;
			}
			break;
	}
	if (com.kernel == NULL) {
		com.kernelsize = 0;
		com.kernelchannels = 1;
		siril_log_color_message(_("Error: no PSF defined. Select blind deconvolution or define a PSF from selection or psf parameters.\n"), "red");
		retval = 1;
		goto END;
	}
#if DEBUG_PSF
// Ignoring the possibility of multichannel kernels here for readability: only the first channel will be shown
	for (int i = 0; i < args.ks; i++) {
		for (int j = 0; j < args.ks; j++) {
			siril_debug_print("%0.2f\t", com.kernel[i * args.ks + j]);
		}
		siril_debug_print("\n");
	}
#endif

	args.kernelorientation = get_imageorientation();
#ifdef DEBUG
	if (args.kernelorientation == BOTTOM_UP)
		siril_log_message(_("PSF made in bottom up orientation.\n"));
	else
		siril_log_message(_("PSF made in top down orientation.\n"));
#endif
	com.kernelsize = (!com.kernel) ? 0 : args.ks;
	com.kernelchannels = (!com.kernel) ? 0 : args.kchans;
	if (args.psftype != PSF_PREVIOUS) {
		gui_function(DrawPSF, NULL);
	}
END:
	return retval;
}

gpointer estimate_only(gpointer p) {
	lock_roi_mutex();
	int retval = 0;
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
	}
	if (args.psftype == PSF_STARS && (!(com.stars && com.stars[0]))) {
		// User wants PSF from stars but has not selected any stars
		siril_log_message(_("No stars detected. Finding suitable non-saturated stars...\n"));
		star_finder_params sfpar;
		memcpy(&sfpar, &com.pref.starfinder_conf, sizeof(star_finder_params));
		sfpar.min_A = 0.07;
		sfpar.max_A = 0.7;
		sfpar.profile = PSF_MOFFAT_BFREE;
		image *input_image = calloc(1, sizeof(image));
		input_image->fit = the_fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;

		int nb_stars;
		int chan = the_fit->naxes[2] > 1 ? 1 : 0; // G channel for color, mono channel for mono
		com.stars = peaker(input_image, chan, &sfpar, &nb_stars, NULL, FALSE, FALSE, MAX_STARS, com.pref.starfinder_conf.profile, com.max_thread);
		free(input_image);
		if (!com.stars || nb_stars == 0) {
			siril_log_color_message(_("No suitable stars detectable in this image. Aborting..."), "red");
			retval = 1;
			goto ENDEST;
		} else {
			args.stars_need_clearing = TRUE;
		}
		siril_log_message(_("Found %d suitably bright, non-saturated stars.\n"), nb_stars);

		calculate_parameters(); // Calculates some parameters based on detected stars
		if (args.recalc_ks && args.ks < (int)(args.psf_fwhm * 4.f)) {
			int newks = (int)(args.psf_fwhm * 4.f);
			if (!(newks%2))
				newks++;
			args.ks = newks;
			siril_log_message(_("Kernel size not large enough to fit PSF generated from detected stars. Increasing kernel size to %d. To override, specify kernel size using the -ks= option.\n"), newks);
		} else if (args.ks < (int)(args.psf_fwhm * 4.f)) {
			int recc_ks = (int)(args.psf_fwhm * 4.f);
				if (!(recc_ks%2))
					recc_ks++;
			siril_log_message(_("Warning: PSF generated from the stars detected in this image appears to be too big for the specified kernel size. Recommend increasing kernel size to %d.\n"), recc_ks);
		}
	}

	args.nchans = the_fit->naxes[2];
	com.fftw_max_thread = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;
	set_wisdom_file();
	if (args.psftype == PSF_BLIND) {
		if (fftwf_import_wisdom_from_filename(com.pref.fftw_conf.wisdom_file) == 1) {
			siril_log_message(_("Siril FFT wisdom imported successfully...\n"));
		} else if (fftwf_import_system_wisdom() == 1) {
			siril_log_message(_("System FFT wisdom imported successfully...\n"));
		} else {
			siril_log_message(_("No FFT wisdom found to import...\n"));
		}
	}
	args.ndata = the_fit->rx * the_fit->ry * the_fit->naxes[2];
	args.fdata = malloc(args.ndata * sizeof(float));
	if (the_fit->type == DATA_FLOAT)
		memcpy(args.fdata, the_fit->fdata, args.ndata * sizeof(float));
	else {
		float invnorm = 1.f / USHRT_MAX_SINGLE;
		for (size_t i = 0 ; i < args.ndata ; i++) {
			args.fdata[i] = (float) the_fit->data[i] * invnorm;
		}
	}

	set_progress_bar_data(_("Starting PSF computation..."), PROGRESS_PULSATE);
	START_TIMER;
	get_kernel();
	END_TIMER;
	if (args.psftype == PSF_BLIND) {
		if (fftwf_export_wisdom_to_filename(com.pref.fftw_conf.wisdom_file) == 1) {
			siril_log_message(_("Siril FFT wisdom updated successfully...\n"));
		} else {
			siril_log_message(_("Siril FFT wisdom update failed...\n"));
		}
	}
	siril_log_color_message(_("Deconvolution PSF generated.\n"), "green");
ENDEST:
	if (args.stars_need_clearing) {
		clear_stars_list(FALSE);
		args.stars_need_clearing = FALSE;
	}
	if(!retval && args.save_after) {
		save_kernel(args.savepsf_filename);
		free(args.savepsf_filename);
		args.savepsf_filename = NULL;
		args.save_after = FALSE;
	}
	unlock_roi_mutex();
	siril_add_idle(estimate_idle, NULL);
	if (com.script) {
		free(args.fdata);
		args.fdata = NULL;
	}
	return GINT_TO_POINTER(retval);
}

gpointer deconvolve(gpointer p) {
	lock_roi_mutex(); // Prevent the ROI change / clear callbacks running until we are
	// done, in order not to have the_fit changed under us
	struct timeval t_start, t_end;
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
		the_fit = &gfit;
	}
	gboolean stars_need_clearing = FALSE;
	check_orientation();
	if (sequence_is_running == 0)
		gui_function(DrawPSF, NULL);
	args.nchans = the_fit->naxes[2];
	args.rx = the_fit->rx;
	args.ry = the_fit->ry;
	int retval = 0;
	if (args.psftype == PSF_PREVIOUS && ((!com.kernel) || com.kernelsize == 0)) {
	// Refuse to process the image using previous PSF if there is no previous PSF defined
		siril_log_color_message(_("Error: trying to use previous PSF but no PSF has been generated. Aborting...\n"),"red");
		retval = 1;
		goto ENDDECONV;
	}
	if (args.psftype == PSF_STARS && (!(com.stars && com.stars[0]))) {
		// User wants PSF from stars but has not selected any stars
		siril_log_message(_("No stars detected. Finding suitable non-saturated stars...\n"));
		star_finder_params sfpar;
		memcpy(&sfpar, &com.pref.starfinder_conf, sizeof(star_finder_params));
		sfpar.min_A = 0.07;
		sfpar.max_A = 0.7;
		sfpar.profile = PSF_MOFFAT_BFREE;
		image *input_image = calloc(1, sizeof(image));
		input_image->fit = the_fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;

		int nb_stars;
		int chan = the_fit->naxes[2] > 1 ? 1 : 0; // G channel for color, mono channel for mono
		com.stars = peaker(input_image, chan, &sfpar, &nb_stars, NULL, FALSE, FALSE, MAX_STARS, com.pref.starfinder_conf.profile, com.max_thread);
		free(input_image);
		if (retval || nb_stars == 0) {
			siril_log_color_message(_("No suitable stars detectable in this image. Aborting..."), "red");
			goto ENDDECONV;
		} else
			stars_need_clearing = TRUE;
	}
	com.fftw_max_thread = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;

	set_cursor_waiting(TRUE);
	set_wisdom_file();
	if (fftwf_import_wisdom_from_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom imported successfully...\n"));
	} else if (fftwf_import_system_wisdom() == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("System FFT wisdom imported successfully...\n"));
	} else {
		if (sequence_is_running == 0)
			siril_log_message(_("No FFT wisdom found to import...\n"));
	}
	if (the_fit == &gfit || the_fit == &gui.roi.fit)
		if (!com.script && !com.headless && !args.previewing)
			undo_save_state(&gfit, _("Deconvolution"));
	args.ndata = the_fit->rx * the_fit->ry * the_fit->naxes[2];
	args.fdata = malloc(args.ndata * sizeof(float));
	if (the_fit->type == DATA_FLOAT)
		memcpy(args.fdata, the_fit->fdata, args.ndata * sizeof(float));
	else {
		float invnorm = 1.f / USHRT_MAX_SINGLE;
		for (size_t i = 0 ; i < args.ndata ; i++) {
			args.fdata[i] = (float) the_fit->data[i] * invnorm;
		}
	}

	// Get the kernel
	if (args.psftype != PSF_PREVIOUS) {
		if (sequence_is_running == 0)
			set_progress_bar_data(_("Starting PSF estimation..."), PROGRESS_PULSATE);
		get_kernel();
}
	if (!com.kernel) {
		siril_debug_print("Kernel missing!\n");
		retval = 1;
		goto ENDDECONV;
	}

	next_psf_is_previous = (args.psftype == PSF_BLIND || args.psftype == PSF_STARS) ? TRUE : FALSE;

	float *yuvdata = NULL;
	int threads = 1;
	if (the_fit->naxes[2] == 3 && com.kernelchannels == 1) {
		// Convert the fit to XYZ and only deconvolve Y
		int npixels = the_fit->rx * the_fit->ry;
		yuvdata = malloc(npixels * the_fit->naxes[2] * sizeof(float));
#ifdef _OPENMP
		threads = sequence_is_running ? 1 : com.max_thread;
#pragma omp parallel for simd num_threads(threads) schedule(static)
#endif
		for (int i = 0 ; i < npixels ; i++) {
			rgb_to_yuvf(args.fdata[i], args.fdata[i + npixels], args.fdata[i + 2 * npixels], &yuvdata[i], &yuvdata[i + npixels], &yuvdata[i + 2 * npixels]);
		}
		args.nchans = 1;
		free(args.fdata);
		args.fdata = yuvdata; // fdata now points to the Y part of yuvdata
	}

	if (get_thread_run() || sequence_is_running == 1) {
		if (sequence_is_running == 0)
			set_progress_bar_data(_("Starting non-blind deconvolution..."), 0);
		gettimeofday(&t_start, NULL);
		// Non-blind deconvolution stage
		switch (args.nonblindtype) {
			case DECONV_SB:
				split_bregman(args.fdata, args.rx, args.ry, args.nchans, com.kernel, args.ks, args.kchans, args.alpha, args.finaliters, com.fftw_max_thread);
				break;
			case DECONV_RL:
				if (args.rl_method == RL_MULT) {
					if (args.regtype == REG_TV_GRAD)
						args.regtype = REG_TV_MULT;
					else if (args.regtype == REG_FH_GRAD)
						args.regtype = REG_FH_MULT;
					else args.regtype = REG_NONE_MULT;
				}
				if (args.ks < com.pref.fftw_conf.fft_cutoff)
					naive_richardson_lucy(args.fdata, args.rx,args.ry, args.nchans, com.kernel, args.ks, args.kchans, args.alpha, args.finaliters, args.stopcriterion, com.fftw_max_thread, args.regtype, args.stepsize, args.stopcriterion_active);
				else
					fft_richardson_lucy(args.fdata, args.rx,args.ry, args.nchans, com.kernel, args.ks, args.kchans, args.alpha, args.finaliters, args.stopcriterion, com.fftw_max_thread, args.regtype, args.stepsize, args.stopcriterion_active);

				break;
			case DECONV_WIENER:
				wienerdec(args.fdata, args.rx, args.ry, args.nchans, com.kernel, args.ks, args.kchans, args.alpha, com.fftw_max_thread);
				break;
		}
		gettimeofday(&t_end, NULL);
		if (sequence_is_running == 0) {
			show_time(t_start, t_end);
		}
	}

	if (the_fit->naxes[2] == 3 && com.kernelchannels == 1) {
		// Put things back as they were
		int npixels = the_fit->rx * the_fit->ry;
		args.nchans = 3;
		args.fdata = malloc(npixels * args.nchans * sizeof(float));
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(threads) schedule(static)
#endif
		for (int i = 0 ; i < npixels ; i++) {
			yuv_to_rgbf(yuvdata[i], yuvdata[i + npixels], yuvdata[i + 2 * npixels], &args.fdata[i], &args.fdata[i + npixels], &args.fdata[i + 2 * npixels]);
		}
		free(yuvdata);
	}

	// Update the_fit with the result
	if (get_thread_run()) {
		if (the_fit->type == DATA_FLOAT) {
			memcpy(the_fit->fdata, args.fdata, args.ndata * sizeof(float));
			if (com.pref.force_16bit)
				fit_replace_buffer(the_fit, float_buffer_to_ushort(the_fit->fdata, args.ndata), DATA_USHORT);
			else invalidate_stats_from_fit(the_fit);
		} else {
			for (size_t i = 0 ; i < args.ndata ; i++) {
				the_fit->data[i] = roundf_to_WORD(args.fdata[i] * USHRT_MAX_SINGLE);
			}
			invalidate_stats_from_fit(the_fit);
		}
	}
ENDDECONV:
	// Do not free the PSF here as it is populated into com.kernel
	if (fftwf_export_wisdom_to_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom updated successfully...\n"));
	} else {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom update failed...\n"));
	}
	if (stars_need_clearing) {
		clear_stars_list(TRUE);
		com.stars = NULL;
		stars_need_clearing = FALSE;
	}
	if (sequence_is_running == 0)
		siril_add_idle(deconvolve_idle, NULL);
	else {
		free(args.fdata);
		args.fdata = NULL;
	}
	unlock_roi_mutex();
	return GINT_TO_POINTER(retval);
}

///////// ****** SEQUENCE PROCESSING ****** //////////

static int deconvolution_compute_mem_limits(struct generic_seq_args *seqargs, gboolean for_writer) {
// The deconvolution routines use FFTW which optimises its planning on the basis of the virtual
// memory model and available threads. It is therefore pragmatic to allow FFTW to optimise itself
// and run the sequence sequentially rather than in parallel. The memory hook therefore only checks
// that at least one image can be processed.

// Allocations table
// Assume kernel size is small c/w image size
// All: allocate float whc (fdata)
// Color images - allocate float whc (xyzdata) but then free fdata
// c == 1 when passed to deconvolve.cpp
// Allocate float wh (img_t<float> f)
// Wiener: 8 x wh (3 x complex float buffers) + 1 x complex buffer used by fftw as the transform is not in-place
// RL (FFT): 13 x wh (4 x complex float buffers, 3 x float buffers) + 1 x complex buffer used by fftw as the transform is not in-place
// RL (naive): 6 x wh (6 x float buffers)
// SB: 16 x wh (!) (6 x complex float buffers + 3 x float buffers) + 1 x complex buffer used by fftw as the transform is not in-place
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(seqargs->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	if (limit > 0) {
		int is_float = get_data_type(seqargs->seq->bitpix) == DATA_FLOAT;
		int float_multiplier = (is_float) ? 1 : 2;
		int MB_per_float_image = MB_per_image * float_multiplier;
		int required = MB_per_float_image * 2;
		int MB_per_float_chan = float_multiplier * ceil((double)MB_per_image / seqargs->seq->nb_layers);
		if (args.nonblindtype == DECONV_SB) {
			required += MB_per_float_chan * 16;
		} else if (args.nonblindtype == DECONV_RL) {
			if (args.ks < com.pref.fftw_conf.fft_cutoff) {
				required += MB_per_float_chan * 6; // Naive deconvolution
			} else {
				required += MB_per_float_chan * 13; // FFT deconvolution
			}
		} else if (args.nonblindtype == DECONV_WIENER) {
			required += MB_per_float_chan * 8;
		}
		limit = MB_avail / required;
	}
	// Limit to 1 image at a time anyway, to support FFTW planning optimizations
	limit = (limit >= 1 ? 1 : 0);
	return limit;
}

int deconvolution_finalize_hook(struct generic_seq_args *seqargs) {
	deconvolution_sequence_data *data = (deconvolution_sequence_data *) seqargs->user;
	int retval = seq_finalize_hook(seqargs);
	if (data && data->from_command && data->deconv_data)
		free(data->deconv_data);
	if (data)
		free(data);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	args.psftype = args.oldpsftype; // Restore consistency
	sequence_is_running = 0;
	return retval;
}


int deconvolution_image_hook(struct generic_seq_args *seqargs, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 0;
	the_fit = fit;
	ret = GPOINTER_TO_INT(deconvolve(NULL));
	the_fit = &gfit; // Prevent bad things happening if fit is freed and we then try to do things with the_fit
	args.oldpsftype = args.psftype; // Need to store the previous psf type so we can restore
	// it later and avoid inconsistency between the GTK widget and the parameter.
	args.psftype = PSF_PREVIOUS; // For all but the first image in the sequence we will reuse the kernel calculated for the first image.

	return ret;
}

int deconvolution_prepare_hook(struct generic_seq_args *seqargs) {
	set_progress_bar_data(_("Deconvolution. Processing sequence..."), 0.);
	remove_prefixed_sequence_files(seqargs->seq, seqargs->new_seq_prefix);
	remove_prefixed_star_files(seqargs->seq, seqargs->new_seq_prefix);
	if (args.psftype == 4 && ((!com.kernel) || com.kernelsize == 0)) {
	// Refuse to process the sequence using previous PSF if there is no previous PSF defined
		siril_log_color_message(_("Error: trying to use previous PSF but no PSF has been generated. Aborting...\n"),"red");
		return 1;
	}
	int retval = 0;
	g_assert(seqargs->has_output); // don't call this hook otherwise
	if (seqargs->force_ser_output || (seqargs->seq->type == SEQ_SER && !seqargs->force_fitseq_output)) {
		gchar *dest;
		const char *ptr = strrchr(seqargs->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s.ser", seqargs->new_seq_prefix, ptr + 1);
		else dest = g_strdup_printf("%s%s.ser", seqargs->new_seq_prefix, seqargs->seq->seqname);

		seqargs->new_ser = calloc(1, sizeof(struct ser_struct));
		if (ser_create_file(dest, seqargs->new_ser, TRUE, seqargs->seq->ser_file)) {
			free(seqargs->new_ser);
			seqargs->new_ser = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else if (seqargs->force_fitseq_output || (seqargs->seq->type == SEQ_FITSEQ && !seqargs->force_ser_output)) {
		gchar *dest;
		const char *ptr = strrchr(seqargs->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s%s", seqargs->new_seq_prefix, ptr + 1, com.pref.ext);
		else dest = g_strdup_printf("%s%s%s", seqargs->new_seq_prefix, seqargs->seq->seqname, com.pref.ext);

		seqargs->new_fitseq = calloc(1, sizeof(fitseq));
		if (fitseq_create_file(dest, seqargs->new_fitseq, seqargs->nb_filtered_images)) {
			free(seqargs->new_fitseq);
			seqargs->new_fitseq = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else return 0;
	if (!retval)
		retval = seq_prepare_writer(seqargs);
	return retval;
}

void apply_deconvolve_to_sequence(struct deconvolution_sequence_data *seqdata) {
	sequence_is_running = 1;
	struct generic_seq_args *seqargs = create_default_seqargs(seqdata->seq);
	seqargs->seq = seqdata->seq;
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seqdata->seq->selnum;
	seqargs->compute_mem_limits_hook = deconvolution_compute_mem_limits;
	seqargs->prepare_hook = deconvolution_prepare_hook;
	seqargs->finalize_hook = deconvolution_finalize_hook;
	seqargs->image_hook = deconvolution_image_hook;
	seqargs->description = _("Deconvolution");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	seqargs->new_seq_prefix = strdup(seqdata->seqEntry);
	seqargs->load_new_sequence = TRUE;
	seqargs->user = seqdata;

	if (!start_in_new_thread(generic_sequence_worker, seqargs)) {
		free(seqdata->seqEntry);
		free(seqdata);
		free_generic_seq_args(seqargs);
	}
}

gpointer deconvolve_sequence_command(gpointer p, sequence* seqname) {
	int retval = 0;
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
	}
	deconvolution_sequence_data* seqargs = calloc(1, sizeof(deconvolution_sequence_data));
	seqargs->seqEntry = strdup("dec_");
	seqargs->seq = seqname;
	seqargs->from_command = TRUE;
	apply_deconvolve_to_sequence(seqargs);
	return GINT_TO_POINTER(retval);
}
