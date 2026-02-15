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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_histogram.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <glib.h>
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#endif
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "yyjson.h"

#include "git-version.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/arithm.h"
#include "core/initfile.h"
#include "core/siril_app_dirs.h"
#include "core/masks.h"
#include "core/preprocess.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "core/siril_networking.h"
#include "core/siril_update.h"
#include "core/undo.h"
#include "io/Astro-TIFF.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/siril_catalogues.h"
#include "io/local_catalogues.h"
#include "io/FITS_symlink.h"
#include "io/fits_keywords.h"
#include "io/siril_pythonmodule.h"
#include "drizzle/cdrizzleutil.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "gui/cut.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/keywords_tree.h"
#include "gui/newdeconv.h"
#include "gui/siril_preview.h"
#include "gui/stacking.h"
#include "gui/registration.h"
#include "gui/registration_preview.h"
#include "gui/script_menu.h"
#include "gui/unpurple.h"
#include "filters/asinh.h"
#include "filters/banding.h"
#include "filters/epf.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "filters/clahe.h"
#include "filters/cosmetic_correction.h"
#include "filters/deconvolution/deconvolution.h"
#include "filters/linear_match.h"
#include "filters/median.h"
#include "filters/mtf.h"
#include "filters/fft.h"
#include "filters/rgradient.h"
#include "filters/saturation.h"
#include "filters/scnr.h"
#include "filters/starnet.h"
#include "filters/synthstar.h"
#include "filters/unpurple.h"
#include "filters/wavelets.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/search_objects.h"
#include "algos/star_finder.h"
#include "algos/Def_Wavelet.h"
#include "algos/background_extraction.h"
#include "algos/ccd-inspector.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"
#include "algos/fitting.h"
#include "algos/colors.h"
#include "algos/quality.h"
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "algos/geometry.h"
#include "algos/photometric_cc.h"
#include "algos/fix_xtrans_af.h"
#include "algos/comparison_stars.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "registration/registration.h"
#include "livestacking/livestacking.h"
#include "pixelMath/pixel_math_runner.h"
#include "io/healpix/healpix_cat.h"

#include "command.h"
#include "command_list.h"
#include "command_line_processor.h"

#define PRINT_DEPRECATED_WARNING(__new_function__) siril_log_color_message(_("This command is deprecated: %s should be used instead.\n"), "red", __new_function__)
#define PRINT_DEPRECATED_OPTION_WARNING(__option__, __new_function__) siril_log_color_message(_("The %s option is deprecated: %s should be used instead.\n"), "red", __option__, __new_function__)

/* process_command functions take the number of arguments (like argc) as
 * parameter and will be able to access the equivalent of argv with `word'
 * they return CMD_OK on success
 */

char *word[MAX_COMMAND_WORDS];	// NULL terminated

// Returns TRUE if the sequence does not contain CFA images
// Otherwise, returns FALSE and prints a warning
static gboolean sequence_cfa_warning_check(sequence* seq) {
	gboolean retval;
	fits tmpfit = { 0 };
	if (seq_read_frame_metadata(seq, sequence_find_refimage(seq), &tmpfit)) {
		return TRUE; // something has gone wrong but we will not highlight it here, the error will be detected and dealt with later
	}
	gboolean mono = (tmpfit.naxes[2] == 1);
	gboolean cfa = (tmpfit.keywords.bayer_pattern[0] != '\0');
	clearfits(&tmpfit);
	if (mono && cfa) {
		control_window_switch_to_tab(OUTPUT_LOGS);
		siril_log_color_message(_("Warning: sequence contains undebayered CFA images. Applying a sequence function not intended for this kind of image. This is likely to give poor results: check your intended workflow.\n"), "salmon");
		retval = FALSE;
	} else {
		retval = TRUE;
	}
	return retval;
}

// Returns TRUE if gfit does not contain a CFA image
// Otherwise, returns FALSE and prints a warning
gboolean image_cfa_warning_check() {
	gboolean retval;
	if (gfit->naxes[2] == 1 && gfit->keywords.bayer_pattern[0] != '\0') {
		control_window_switch_to_tab(OUTPUT_LOGS);
		siril_log_color_message(_("Warning: an undebayered CFA image is loaded. Applying an image function not intended for this kind of image. This is likely to give poor results: check your intended workflow.\n"), "salmon");
		retval = FALSE;
	} else {
		retval = TRUE;
	}
	return retval;
}

int process_load(int nb){
	long maxpath = get_pathmax();
	char filename[maxpath];
	size_t len = strlen(word[1]);
	strncpy(filename, word[1], maxpath - 1);
	filename[maxpath - 1] = '\0';
	for (int i = 1; i < nb - 1; ++i) {
		strncat(filename, " ", maxpath - 1 - len);
		len += 1;
		strncat(filename, word[i + 1], len);
	}
	expand_home_in_filename(filename, maxpath);

	int retval = open_single_image(filename);
	gui_function(launch_clipboard_survey, NULL);
	return (retval == 0) ? CMD_OK : CMD_FILE_NOT_FOUND;
}

int process_load_seq(int nb) {
	close_sequence(FALSE);
	close_single_image();

	// Load the sequence into com.seq
	set_seq(word[1]);
	if (com.seq.seqname) {
		siril_debug_print("Sequence loaded ok\n");
		return CMD_OK;
	} else {
		return CMD_GENERIC_ERROR;
	}
}

int process_dumpheader(int nb) {
	if (!gfit->header || gfit->header[0] == '\0') {
		siril_log_message(_("Currently loaded image has no FITS header\n"));
	} else {
		siril_log_message(_("=FITS header for currently loaded image=\n"));
		char *header = strdup(gfit->header);
		log_several_lines(header);
		free(header);
	}
	return CMD_OK;
}

int process_seq_clean(int nb) {
	gboolean cleanreg = FALSE, cleanstat = FALSE, cleansel = FALSE;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	if (nb > 2) {
		for (int i = 2; i < nb; i++) {
			if (word[i]) {
				if (!strcmp(word[i], "-reg")) {
					cleanreg = TRUE;
				}
				else if (!strcmp(word[i], "-stat")) {
					cleanstat = TRUE;
				}
				else if (!strcmp(word[i], "-sel")) {
					cleansel = TRUE;
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					return CMD_ARG_ERROR;
				}
			}
		}
	} else {
		cleanreg = TRUE;
		cleanstat = TRUE;
		cleansel = TRUE;
	}

	clean_sequence(seq, cleanreg, cleanstat, cleansel);
	if (check_seq_is_comseq(seq)) {
		fix_selnum(&com.seq, FALSE);
		update_stack_interface(TRUE);
		update_reg_interface(FALSE);
		adjust_sellabel();
		set_layers_for_registration();
		drawPlot();
	} else {
		free_sequence(seq, TRUE);
	}
	return CMD_OK;
}

int process_satu(int nb) {
	gboolean mask_aware = FALSE;
	const char *num_args[3] = { NULL, NULL, NULL };
	int num_index = 0;
	gchar *end;

	/* Scan arguments (order-independent for -mask) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (num_index < 3) {
			num_args[num_index++] = word[i];
		} else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	/* coeff (required) */
	if (num_index < 1) {
		siril_log_message(_("Missing saturation coefficient\n"));
		return CMD_ARG_ERROR;
	}

	double coeff = g_ascii_strtod(num_args[0], &end);
	if (end == num_args[0]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), num_args[0]);
		return CMD_ARG_ERROR;
	}

	/* background_factor (optional) */
	double background_factor = 1.0;
	if (num_index > 1) {
		background_factor = g_ascii_strtod(num_args[1], &end);
		if (end == num_args[1] || background_factor < 0.0) {
			siril_log_message(_("Background factor must be positive\n"));
			return CMD_ARG_ERROR;
		}
	}

	/* satu_hue_type (optional) */
	int satu_hue_type = 6;
	if (num_index > 2) {
		satu_hue_type = g_ascii_strtoull(num_args[2], &end, 10);
		if (end == num_args[2] || satu_hue_type > 6) {
			siril_log_message(_("Hue range must be [0, 6]\n"));
			return CMD_ARG_ERROR;
		}
	}

	/* Allocate parameters */
	saturation_params *params = calloc(1, sizeof(saturation_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->free = free; // Destructor
	params->coeff = coeff;
	params->background_factor = background_factor;
	satu_set_hues_from_types(params, satu_hue_type);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = saturation_image_hook;
	args->idle_function = NULL; // Synchronous execution
	args->description = _("Saturation");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = FALSE;
	args->user = params;
	args->mask_aware = mask_aware;
	args->log_hook = satu_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	/* Run synchronously */
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	int retval = args->retval;

	if (retval != 0) {
		siril_log_color_message(_("Saturation enhancement failed\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	char log[90];
	sprintf(log, "Color saturation %d%%, threshold %.2f",
		round_to_int(coeff * 100.0), background_factor);
	gfit->history = g_slist_append(gfit->history, strdup(log));
	siril_log_message(_("%s\n"), log);

	return CMD_OK;
}

int process_save(int nb){
	gchar *filename = g_strdup(word[1]);
	int status, retval;

	gfit->checksum = FALSE;
	for (int i = 2; i < nb; i++) {
		if (word[i] && !g_strcmp0(word[i], "-chksum")) {
			gfit->checksum = TRUE;
		}
	}

	if (!com.script) {
		gfit->keywords.lo = gui.lo;
		gfit->keywords.hi = gui.hi;
	}

	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = savefits(savename, gfit) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	if (com.uniq) {
		if (!g_str_has_suffix(savename, com.pref.ext)) {
			gchar* tempfilename = g_strdup_printf("%s%s", savename, com.pref.ext);
			com.uniq->filename = strdup(tempfilename);
			g_free(tempfilename);
		} else {
			com.uniq->filename = strdup(savename);
		}
		com.uniq->fileexist = TRUE;
		if (!com.headless) {
			display_filename();
		}
	}
	gui_function(set_precision_switch, NULL);

	g_free(filename);
	g_free(savename);
	return retval;
}

int process_savebmp(int nb){
	gchar *filename = g_strdup_printf("%s.bmp", word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = savebmp(savename, gfit) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}

int process_synthstar(int nb) {
	set_cursor_waiting(TRUE);
	image_cfa_warning_check();
	start_in_new_thread(do_synthstar, NULL);
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

int process_unclip(int nb) {
	set_cursor_waiting(TRUE);
	image_cfa_warning_check();
	start_in_new_thread(fix_saturated_stars, NULL);
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

// TODO: check this against the new denoise idles
/*
static gboolean end_denoise(gpointer p) {
	stop_processing_thread();
	struct denoise_args *args = (struct denoise_args *) p;
	if (!args->previewing) {
		copy_gfit_to_backup();
		populate_roi();
	}
	notify_gfit_modified();
	queue_redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
	free(args);
	return FALSE;
}*/

gchar *denoise_log_hook(gpointer p, log_hook_detail detail) {
	denoise_args *args = (denoise_args *) p;

	GString *gs = g_string_new(NULL);

	/* Base message */
	g_string_append_printf(gs, _("NL-Bayes denoise (mod=%.3f"), args->modulation);

	/* Optional DA3D / SOS / VST info */
	if (args->da3d) {
		g_string_append(gs, _(", DA3D enabled"));
	} else if (args->sos > 1) {
		g_string_append_printf(gs, _(", SOS enabled (iters=%d, rho=%.3f)"),
							args->sos, args->rho);
	} else if (args->do_anscombe) {
		g_string_append(gs, _(", VST enabled"));
	}

	/* Optional CC info */
	if (args->do_cosme) {
		g_string_append(gs, _(", CC enabled)"));
	} else {
		g_string_append(gs, _(")"));
	}

	/* Convert to gchar* */
	gchar *result = g_strdup(gs->str);
	g_string_free(gs, TRUE);

	return result;  /* caller must g_free() */
}

gpointer run_nlbayes_on_fit(gpointer p) {
	lock_roi_mutex();
	copy_backup_to_gfit();
	denoise_args *args = (denoise_args *) p;

	if (args->suppress_artefacts)
		siril_log_message(_("Colour artefact suppression active.\n"));

	int retval = 0;

	// Carry out cosmetic correction at the start, if selected
	if (args->do_cosme)
		denoise_hook_cosmetic(args->fit);

	// Apply NR to each channel independently
	if (args->fit == gfit && args->fit->naxes[2] == 3 && args->suppress_artefacts) {
		fits *loop = NULL;
		if (new_fit_image(&loop, args->fit->rx, args->fit->ry, 1, args->fit->type)) {
			retval = 1;
		}
		loop->naxis = 2;
		loop->naxes[2] = 1;
		size_t npixels = args->fit->naxes[0] * args->fit->naxes[1];
		if (retval == 0) {
			if (args->fit->type == DATA_FLOAT) {
				for (size_t i = 0; i < 3; i++) {
					float *loop_fdata = (float*) calloc(npixels, sizeof(float));
					free(loop->fdata);
					loop->fdata = loop_fdata;
					memcpy(loop_fdata, args->fit->fpdata[i], npixels * sizeof(float));
					retval = do_nlbayes(loop, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
					memcpy(args->fit->fpdata[i], loop->fdata, npixels * sizeof(float));
				}
				free(loop->fdata);
				loop->fdata = NULL;
			} else {
				for (size_t i = 0; i < 3; i++) {
					WORD *loop_data = (WORD*) calloc(npixels, sizeof(WORD));
					free(loop->data);
					loop->data = loop_data;
					memcpy(loop_data, args->fit->pdata[i], npixels * sizeof(WORD));
					retval = do_nlbayes(loop, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
					memcpy(args->fit->pdata[i], loop->data, npixels * sizeof(WORD));
				}
				free(loop->data);
				loop->data = NULL;
			}
		clearfits(loop);
		}
	} else {
		retval = do_nlbayes(args->fit, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
	}
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	unlock_roi_mutex();
	return GINT_TO_POINTER(retval | CMD_NOTIFY_GFIT_MODIFIED);
}

/*****************************************************************************
 *      D E N O I S E   A L L O C A T O R   A N D   D E S T R U C T O R     *
 ****************************************************************************/

/* Allocator for denoise_args */
struct denoise_args *new_denoise_args() {
	struct denoise_args *args = calloc(1, sizeof(struct denoise_args));
	if (args) {
		args->destroy_fn = free_denoise_args;
		// Set default values
		args->modulation = 1.0f;
		args->sos = 1;
		args->rho = 0.2f;
		args->da3d = 0;
		args->do_anscombe = FALSE;
		args->do_cosme = TRUE;
		args->suppress_artefacts = FALSE;
		args->previewing = FALSE;
	}
	return args;
}

/* Destructor for denoise_args */
void free_denoise_args(void *ptr) {
	struct denoise_args *args = (struct denoise_args *)ptr;
	if (!args)
		return;
	free(ptr);
}

/*****************************************************************************
 *      D E N O I S E   P R O C E S S I N G   H O O K                       *
 ****************************************************************************/

/* The actual denoising processing hook */
int denoise_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct denoise_args *params = (struct denoise_args *)args->user;
	if (!params) {
		siril_log_color_message(_("Denoising failed: invalid parameters\n"), "red");
		return 1;
	}

	// Update fit pointer to the one being processed
	params->fit = fit;

	// Log what we're doing (only if verbose)
	if (args->verbose) {
		if (params->do_anscombe)
			siril_log_message(_("Applying generalised Anscombe variance stabilising transform.\n"));
		if (params->da3d)
			siril_log_message(_("Carrying out final stage DA3D denoising.\n"));
		if (params->sos > 1)
			siril_log_message(_("Performing %d SOS iterations with rho=%.2f.\n"), params->sos, params->rho);
		siril_log_message(_("Cosmetic correction: %s, Suppress artefacts: %s\n"),
		                  params->do_cosme ? _("enabled") : _("disabled"),
		                  params->suppress_artefacts ? _("enabled") : _("disabled"));
	}

	// Call the actual denoising function
	// Note: do_nlbayes doesn't take cosmetic or suppress_artefacts parameters
	// These might need to be applied separately or the function signature needs updating
	int retval = do_nlbayes(fit, params->modulation, params->sos, params->da3d,
	                        params->rho, params->do_anscombe);

	if (retval != 0 && args->verbose) {
		siril_log_color_message(_("NL-Bayes denoising failed.\n"), "red");
	}

	return retval;
}

/*****************************************************************************
 *      C O M M A N D   P R O C E S S I N G                                 *
 ****************************************************************************/

int process_denoise(int nb) {
	set_cursor_waiting(TRUE);

	gboolean mask_aware = FALSE;

	// Allocate parameters
	struct denoise_args *params = new_denoise_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	// Scan arguments (order-independent)
	for (int i = 1; i < nb; i++) {
		char *arg = word[i], *end;
		if (!arg)
			continue;

		if (!g_strcmp0(arg, "-mask")) {
			mask_aware = TRUE;
		}
		else if (g_str_has_prefix(arg, "-vst")) {
			params->do_anscombe = TRUE;
		}
		else if (g_str_has_prefix(arg, "-da3d")) {
			params->da3d = 1;
		}
		else if (g_str_has_prefix(arg, "-indep")) {
			params->suppress_artefacts = TRUE;
		}
		else if (g_str_has_prefix(arg, "-nocosmetic")) {
			params->do_cosme = FALSE;
		}
		else if (g_str_has_prefix(arg, "-mod=")) {
			arg += 5;
			float mod = (float) g_ascii_strtod(arg, &end);
			if (end == arg || mod <= 0.f || mod > 1.f) {
				siril_log_message(_("Error: modulation must be > 0.0 and <= 1.0.\n"));
				free_denoise_args(params);
				set_cursor_waiting(FALSE);
				return CMD_ARG_ERROR;
			}
			params->modulation = mod;
		}
		else if (g_str_has_prefix(arg, "-rho=")) {
			arg += 5;
			float rho = g_ascii_strtod(arg, &end);
			if (end == arg || rho <= 0.f || rho >= 1.f) {
				siril_log_message(_("Error in rho parameter: must be strictly > 0 and < 1, aborting.\n"));
				free_denoise_args(params);
				set_cursor_waiting(FALSE);
				return CMD_ARG_ERROR;
			}
			params->rho = rho;
		}
		else if (g_str_has_prefix(arg, "-sos=")) {
			arg += 5;
			unsigned sos = (unsigned) g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Error parsing SOS iterations.\n"));
				free_denoise_args(params);
				free(params);
				set_cursor_waiting(FALSE);
				return CMD_ARG_ERROR;
			} else if (sos < 1) {
				siril_log_message(_("Note: SOS iterations < 1. Defaulting to 1.\n"));
				sos = 1;
			} else if (sos > 10) {
				siril_log_message(_("Note: high number of SOS iterations. Processing time may be lengthy...\n"));
			}
			params->sos = sos;
		}
		else {
			siril_log_message(_("Unknown argument %s\n"), arg);
			free_denoise_args(params);
			set_cursor_waiting(FALSE);
			return CMD_ARG_ERROR;
		}
	}

	// Validate parameter combinations
	if (params->do_anscombe && (params->sos != 1 || params->da3d)) {
		siril_log_color_message(
			_("Error: will not carry out DA3D or SOS iterations with Anscombe transform VST selected. Aborting.\n"),
			"red"
		);
		free_denoise_args(params);
		set_cursor_waiting(FALSE);
		return CMD_ARG_ERROR;
	}

	if (params->do_anscombe)
		siril_log_message(_("Will apply generalised Anscombe variance stabilising transform.\n"));

	if (params->da3d) {
		siril_log_message(_("Will carry out final stage DA3D denoising.\n"));
		if (params->sos != 1) {
			siril_log_message(_("Will not carry out both DA3D and SOS. SOS iterations set to 1.\n"));
			params->sos = 1;
		}
	}

	image_cfa_warning_check();

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_denoise_args(params);
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	// Set up generic_img_args
	args->fit = gfit;
	args->mem_ratio = 3.0f; // Denoising needs extra memory
	args->image_hook = denoise_image_hook;
	args->idle_function = denoise_apply_idle;
	args->description = _("NL-Bayes Denoising");
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = denoise_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->mask_aware = mask_aware;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	set_cursor_waiting(FALSE);

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_starnet(int nb) {
#ifdef HAVE_LIBTIFF
	if (!com.pref.starnet_exe || com.pref.starnet_exe[0] == '\0') {
		siril_log_color_message(_("Error: no StarNet executable set.\n"), "red");
		return CMD_FILE_NOT_FOUND;
	}
	if (starnet_executablecheck(com.pref.starnet_exe) == NIL) {
		siril_log_color_message(_("Error: StarNet executable (%s) is not valid.\n"),
			"red", com.pref.starnet_exe);
		return CMD_GENERIC_ERROR;
	}

	// Allocate StarNet parameters
	starnet_data *starnet_params = new_starnet_args();
	if (!starnet_params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	// Set defaults
	starnet_params->linear = FALSE;
	starnet_params->customstride = FALSE;
	starnet_params->upscale = FALSE;
	starnet_params->starmask = TRUE;
	starnet_params->follow_on = FALSE;
	starnet_params->starnet_fit = gfit;

	/* Scan arguments (order-independent for -mask) */
	gboolean mask_aware = FALSE;
	for (int i = 1; i < nb; i++) {
		char *arg = word[i];
		char *end;

		if (!arg) continue;

		if (!g_strcmp0(arg, "-mask")) {
			mask_aware = TRUE;
		} else if (g_str_has_prefix(arg, "-stretch")) {
			starnet_params->linear = TRUE;
		} else if (g_str_has_prefix(arg, "-upscale")) {
			starnet_params->upscale = TRUE;
		} else if (g_str_has_prefix(arg, "-nostarmask")) {
			starnet_params->starmask = FALSE;
		} else if (g_str_has_prefix(arg, "-stride=")) {
			arg += 8;
			int stride = (int) g_ascii_strtod(arg, &end);
			if (arg == end) {
				siril_log_message(_("Error parsing stride argument, aborting.\n"));
				free_starnet_args(starnet_params);
				return CMD_ARG_ERROR;
			} else if ((stride < 2) || (stride > 512) || (stride % 2)) {
				siril_log_message(_("Stride must be an even integer in [2, 512], aborting.\n"));
				free_starnet_args(starnet_params);
				return CMD_ARG_ERROR;
			}
			starnet_params->stride = g_strdup_printf("%d", stride);
			starnet_params->customstride = TRUE;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			free_starnet_args(starnet_params);
			return CMD_ARG_ERROR;
		}
	}

	image_cfa_warning_check();

	// Save backup for undo before processing
	copy_gfit_to_backup();

	// Allocate generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_starnet_args(starnet_params);
		return CMD_ALLOC_ERROR;
	}

	// Set up generic_img_args
	args->fit = gfit;
	args->mem_ratio = 3.0f;
	args->image_hook = starnet_single_image_hook;
	args->idle_function = starnet_single_image_idle;
	args->description = _("StarNet");
	args->verbose = TRUE;
	args->user = starnet_params;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->mask_aware = mask_aware;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
#else
	siril_log_message(_("starnet command unavailable as Siril has not been compiled with libtiff.\n"));
	return CMD_NOT_FOR_THIS_OS;
#endif
}

int process_seq_starnet(int nb){
#ifdef HAVE_LIBTIFF
	if (!com.pref.starnet_exe || com.pref.starnet_exe[0] == '\0') {
		siril_log_color_message(_("Error: no StarNet executable set.\n"), "red");
		return CMD_FILE_NOT_FOUND;
	}
	if (starnet_executablecheck(com.pref.starnet_exe) == NIL) {
		siril_log_color_message(_("Error: StarNet executable (%s) is not valid.\n"), "red", com.pref.starnet_exe);
		return CMD_GENERIC_ERROR;
	}
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}
	struct multi_output_data *multi_args = calloc(1, sizeof(struct multi_output_data));
	if (!multi_args)
		return CMD_ALLOC_ERROR;
	starnet_data *starnet_args = new_starnet_args();
	if (!starnet_args)
		return CMD_ALLOC_ERROR;
	starnet_args->linear = FALSE;
	starnet_args->customstride = FALSE;
	starnet_args->upscale = FALSE;
	starnet_args->starmask = TRUE;
	starnet_args->follow_on = FALSE;
	starnet_args->multi_args = multi_args;
	gboolean error = FALSE;
	multi_args->seq = seq;
	if (!multi_args->seq) {
		siril_log_message(_("Error: cannot open sequence\n"));
		free_starnet_args(starnet_args);
		free(multi_args);
		return CMD_SEQUENCE_NOT_FOUND;
	}

	for (int i = 2; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(arg, "-stretch")) {
			starnet_args->linear = TRUE;
		}
		else if (g_str_has_prefix(arg, "-upscale")) {
			starnet_args->upscale = TRUE;
		}
		else if (g_str_has_prefix(arg, "-nostarmask")) {
			starnet_args->starmask = FALSE;
		}
		else if (g_str_has_prefix(arg, "-stride=")) {
			arg += 8;
			int stride = (int) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((stride < 2.0) || (stride > 512) || (stride % 2)) {
				siril_log_message(_("Error in stride parameter: must be a positive even integer, max 512, aborting.\n"));
				if (!check_seq_is_comseq(multi_args->seq))
					free_sequence(multi_args->seq, TRUE);
				free_starnet_args(starnet_args);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				starnet_args->stride = g_strdup_printf("%d", stride);
				starnet_args->customstride = TRUE;
			}
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				if (!check_seq_is_comseq(multi_args->seq))
					free_sequence(multi_args->seq, TRUE);
				free_starnet_args(starnet_args);
			return CMD_ARG_ERROR;
		}
		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
			if (!check_seq_is_comseq(multi_args->seq))
				free_sequence(multi_args->seq, TRUE);
			free_starnet_args(starnet_args);
			return CMD_ARG_ERROR;
		}
	}
	multi_args->user_data = (gpointer) starnet_args;
	multi_args->n = starnet_args->starmask ? 2 : 1;
	multi_args->prefixes = calloc(multi_args->n, sizeof(gchar*));
	multi_args->prefixes[0] = g_strdup("starless_");
	if (starnet_args->starmask) {
		multi_args->prefixes[1] = g_strdup("starmask_");
	}
	multi_args->seqEntry = strdup(multi_args->prefixes[0]);
	sequence_cfa_warning_check(multi_args->seq);
	set_cursor_waiting(TRUE);
	apply_starnet_to_sequence(multi_args);

#else
	siril_log_message(_("starnet command unavailable as Siril has not been compiled with libtiff.\n"));
#endif
	return CMD_OK;
}

#ifdef HAVE_LIBJPEG
int process_savejpg(int nb){
	int quality = 100;

	if (nb == 3) {
		gchar *end;
		quality = g_ascii_strtoull(word[2], &end, 10);
		if (end == word[2] || quality < 10 || quality > 100) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[2]);
			return CMD_ARG_ERROR;
		}
	}

	gchar *filename = g_strdup_printf("%s.jpg", word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = savejpg(savename, gfit, quality, TRUE) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}
#endif

#ifdef HAVE_LIBJXL
int process_savejxl(int nb){
	int effort = 7;
	double quality = 94.0;
	gboolean force_8bit = FALSE;
	for (int i = 2; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(arg, "-8bit")) {
			force_8bit = TRUE;
		}
		else if (g_str_has_prefix(arg, "-quality=")) {
			arg += 9;
			double quality = g_ascii_strtod(arg, &end);
			if (quality <= 0.0 || quality > 100.0) {
				siril_log_message(_("Error: quality must be >= 0.0 and <= 100.0.\n"));
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-effort=")) {
			arg += 8;
			effort = (int) g_ascii_strtod(arg, &end);
			if (effort < 1.0 || effort > 9.0) {
				siril_log_message(_("Error: effort must be an integer between 1 and 9.\n"));
				return CMD_ARG_ERROR;
			}
		}
	}

	gchar *filename = g_strdup_printf("%s.jxl", word[1]);
	int status, retval = CMD_OK;
	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = savejxl(savename, gfit, effort, quality, force_8bit) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}
#endif

#ifdef HAVE_LIBPNG
int process_savepng(int nb){
	gchar *filename = g_strdup_printf("%s.png", word[1]);
	int status, retval;
	gchar *savename =update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		uint32_t bytes_per_sample = gfit->orig_bitpix != BYTE_IMG ? 2 : 1;
		retval = savepng(savename, gfit, bytes_per_sample, gfit->naxes[2] == 3) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	uint16_t bitspersample = 16;
	gboolean tiff_compression = FALSE;
	gchar *astro_tiff = NULL;

	if (strcasecmp(word[0], "savetif8") == 0)
		bitspersample = 8;
	else if (strcasecmp(word[0], "savetif32") == 0)
		bitspersample = 32;

	for (int i = 2; i < nb; i++) {
		if (word[i] && !g_strcmp0(word[i], "-astro")) {
			if (astro_tiff)
				g_free(astro_tiff);
			astro_tiff = AstroTiff_build_header(gfit);
		} else if (word[i] && !g_strcmp0(word[i], "-deflate")) {
			tiff_compression = TRUE;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			if (astro_tiff)
				g_free(astro_tiff);
			return CMD_ARG_ERROR;
		}
	}

	gchar *filename = g_strdup_printf("%s.tif", word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = savetif(savename, gfit, bitspersample, astro_tiff, com.pref.copyright, tiff_compression, TRUE, TRUE) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	free(astro_tiff);
	g_free(filename);
	g_free(savename);
	return retval;
}
#endif

int process_savepnm(int nb){
	gchar *filename = g_strdup(word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = CMD_GENERIC_ERROR;
	} else {
		set_cursor_waiting(TRUE);
		retval = saveNetPBM(savename, gfit) ? CMD_GENERIC_ERROR : CMD_OK;
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}

static char* normalize_rebayerfilename(char *filename_buffer, const char *input, long int maxpath) {
	strncpy(filename_buffer, input, maxpath);
	filename_buffer[maxpath-1] = '\0';
	expand_home_in_filename(filename_buffer, maxpath);
	return filename_buffer;
}

int process_rebayer(int nb){
	long maxpath = get_pathmax();
	char filename[maxpath];
	fits cfa0 = { 0 };
	fits cfa1 = { 0 };
	fits cfa2 = { 0 };
	fits cfa3 = { 0 };
	fits *out = NULL;
	sensor_pattern pattern = -1;
	if (nb < 5) {
		siril_log_color_message(_("Error, requires at least 4 arguments to specify the 4 files!\n"), "red");
		return CMD_WRONG_N_ARG;
	}

	if (!strcmp(word[5], "RGGB")) {
		pattern = BAYER_FILTER_RGGB;
		siril_log_message(_("Reconstructing RGGB Bayer matrix.\n"));
	} else if (!strcmp(word[5], "BGGR")) {
		pattern = BAYER_FILTER_BGGR;
		siril_log_message(_("Reconstructing BGGR Bayer matrix.\n"));
	} else if (!strcmp(word[5], "GBRG")) {
		pattern = BAYER_FILTER_GBRG;
		siril_log_message(_("Reconstructing GBRG Bayer matrix.\n"));
	} else if (!strcmp(word[5], "GRBG")) {
		pattern = BAYER_FILTER_GRBG;
		siril_log_message(_("Reconstructing GRBG Bayer matrix.\n"));
	} else {
		siril_log_color_message(_("Invalid Bayer matrix specified!\n"), "red");
		return CMD_ARG_ERROR;
	}

	set_cursor_waiting(TRUE);
	int retval = readfits(normalize_rebayerfilename(filename, word[1], maxpath), &cfa0, NULL, FALSE);
	retval += readfits(normalize_rebayerfilename(filename, word[2], maxpath), &cfa1, NULL, FALSE);
	retval += readfits(normalize_rebayerfilename(filename, word[3], maxpath), &cfa2, NULL, FALSE);
	retval += readfits(normalize_rebayerfilename(filename, word[4], maxpath), &cfa3, NULL, FALSE);
	if (retval) {
		siril_log_color_message(_("Error loading files!\n"), "red");
		set_cursor_waiting(FALSE);
		return CMD_FILE_NOT_FOUND;
	}

	out = merge_cfa(&cfa0, &cfa1, &cfa2, &cfa3, pattern);
	siril_log_message("Bayer pattern produced: 1 layer, %dx%d pixels\n", out->rx, out->ry);
	close_single_image();
	copyfits(out, gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copy_fits_metadata(out, gfit);
	update_sampling_information(gfit, 0.5f);
	update_bayer_pattern_information(gfit, pattern);
	free_wcs(gfit);
	update_fits_header(gfit);
	clearfits(out);
	free(out);

	clear_stars_list(TRUE);
	com.seq.current = UNRELATED_IMAGE;
	if (!create_uniq_from_gfit(strdup(_("Unsaved Bayer pattern merge")), FALSE))
		com.uniq->comment = strdup(_("Bayer pattern merge"));

	if (!com.script) {
		gui_function(open_single_image_from_gfit, NULL);
	}
	set_cursor_waiting(FALSE);
	return CMD_OK | CMD_NOTIFY_GFIT_MODIFIED;
}

// Structure to hold imoper-specific data
struct imoper_data {
	void (*destructor)(void *);  // Required as first member
	image_operator oper;
	char *filename;
	gboolean force_to_float;
};

// Destructor for imoper_data
static void free_imoper_data(void *data) {
	struct imoper_data *imoper = (struct imoper_data *)data;
	if (!imoper)
		return;
	free(imoper->filename);
	free(imoper);
}

// Image processing hook for imoper
static int imoper_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct imoper_data *data = (struct imoper_data *)args->user;

	// Read the second image
	fits operand = { 0 };
	if (readfits(data->filename, &operand, NULL, data->force_to_float)) {
		siril_log_color_message(_("Failed to read operand image\n"), "red");
		return 1;
	}

	// Perform the operation
	int retval = imoper(fit, &operand, data->oper, data->force_to_float);

	clearfits(&operand);
	return retval;
}

static gchar* imoper_log_hook(gpointer p, log_hook_detail detail) {
	struct imoper_data *data = (struct imoper_data*) p;
	gchar *message = NULL;
	if (detail == SUMMARY) {
	switch (data->oper) {
		case OPER_ADD:
			message = g_strdup_printf(_("Image addition"));
			break;
		case OPER_SUB:
			message = g_strdup_printf(_("Image subtraction"));
			break;
		case OPER_MUL:
			message = g_strdup_printf(_("Image multiplication"));
			break;
		default:
			message = g_strdup_printf(_("Image division"));
	}	} else {
		switch (data->oper) {
			case OPER_ADD:
				message = g_strdup_printf(_("Image addition of %s"), data->filename);
				break;
			case OPER_SUB:
				message = g_strdup_printf(_("Image subtraction of %s"), data->filename);
				break;
			case OPER_MUL:
				message = g_strdup_printf(_("Image multiplication by %s"), data->filename);
				break;
			default:
				message = g_strdup_printf(_("Image division by %s"), data->filename);
		}
	}
	return message;
}

// Main command function
int process_imoper(int nb) {
	gboolean mask_aware = FALSE;
	const char *filename = NULL;

	/* Scan arguments (order-independent) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!filename) {
			filename = word[i];
		} else {
			siril_log_color_message(_("Invalid argument %s\n"), "red", word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!filename) {
		siril_log_color_message(_("Missing operand image filename\n"), "red");
		return CMD_ARG_ERROR;
	}

	image_operator oper;
	switch (word[0][1]) {
		case 'a':
		case 'A':
			oper = OPER_ADD;
			break;
		case 's':
		case 'S':
			oper = OPER_SUB;
			break;
		case 'm':
		case 'M':
			oper = OPER_MUL;
			break;
		case 'd':
		case 'D':
			oper = OPER_DIV;
			break;
		default:
			siril_log_color_message(_("Could not understand the requested operator\n"), "red");
			return CMD_ARG_ERROR;
	}

	/* Allocate and initialize user data */
	struct imoper_data *imoper = calloc(1, sizeof(struct imoper_data));
	if (!imoper) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	imoper->destructor = free_imoper_data;
	imoper->oper = oper;
	imoper->filename = strdup(filename);
	imoper->force_to_float = !com.pref.force_16bit;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_imoper_data(imoper);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;  // Need memory for original + operand image
	args->image_hook = imoper_image_hook;
	args->idle_function = NULL;  // Use default
	switch (oper) {
		case OPER_ADD:
			args->description = _("Image addition");
			break;
		case OPER_SUB:
			args->description = _("Image subtraction");
			break;
		case OPER_MUL:
			args->description = _("Image multiplication");
			break;
		default:
			args->description = _("Image division");
	}
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;  // This command modifies gfit
	args->command = TRUE; // calling as command, not from GUI
	args->user = imoper;
	args->log_hook = imoper_log_hook;
	args->max_threads = 1;  // imoper likely doesn't need multi-threading
	args->for_preview = FALSE;
	args->mask_aware = mask_aware;
	args->for_roi = FALSE;

	/* Start the worker thread */
	start_in_new_thread(generic_image_worker, args);

	return CMD_OK;
}

struct addmax_data {
	void (*destructor)(void *);
	fits *operand_fit;
};

static void addmax_destructor(void *p) {
	struct addmax_data *data = (struct addmax_data *)p;
	if (data->operand_fit) {
		clearfits(data->operand_fit);
		free(data->operand_fit);
	}
	free(data);
}

static int addmax_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct addmax_data *data = (struct addmax_data *)args->user;
	addmax(fit, data->operand_fit);
	return 0;
}

static gchar *addmax_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup_printf(_("Add max operation"));
}

int process_addmax(int nb) {
	if (nb < 2) {
		siril_log_message(_("Missing operand filename\n"));
		return CMD_ARG_ERROR;
	}

	gboolean mask_aware = FALSE;
	char *filename = NULL;

	/* Scan arguments (order-independent) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!filename) {
			filename = word[i];
		} else {
			siril_log_message(_("Unknown argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!filename) {
		siril_log_message(_("Operand filename not specified\n"));
		return CMD_ARG_ERROR;
	}

	/* Allocate and read operand fits */
	fits *operand_fit = calloc(1, sizeof(fits));
	if (!operand_fit) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	if (readfits(filename, operand_fit, NULL, gfit->type == DATA_FLOAT)) {
		free(operand_fit);
		return CMD_INVALID_IMAGE;
	}

	/* Allocate and initialize user data */
	struct addmax_data *data = calloc(1, sizeof(struct addmax_data));
	if (!data) {
		clearfits(operand_fit);
		free(operand_fit);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	data->destructor = addmax_destructor;
	data->operand_fit = operand_fit;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		addmax_destructor(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = addmax_image_hook;
	args->idle_function = NULL;
	args->description = _("Add max");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = addmax_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

struct fdiv_data {
	void (*destructor)(void *);
	fits *operand_fit;
	float norm;
};

static void fdiv_destructor(void *p) {
	struct fdiv_data *data = (struct fdiv_data *)p;
	if (data->operand_fit) {
		clearfits(data->operand_fit);
		free(data->operand_fit);
	}
	free(data);
}

static int fdiv_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct fdiv_data *data = (struct fdiv_data *)args->user;
	siril_fdiv(fit, data->operand_fit, data->norm, TRUE);
	return 0;
}

static gchar *fdiv_log_hook(gpointer p, log_hook_detail detail) {
	struct fdiv_data *args = (struct fdiv_data*) p;
	return g_strdup_printf(_("Image division with normalization: %.6f"), args->norm);
}

int process_fdiv(int nb) {
	if (nb < 3) {
		siril_log_message(_("Missing arguments: <operand> <norm>\n"));
		return CMD_ARG_ERROR;
	}

	gboolean mask_aware = FALSE;
	char *filename = NULL;
	char *norm_arg = NULL;

	/* Scan arguments (order-independent) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!filename) {
			filename = word[i];
		} else if (!norm_arg) {
			norm_arg = word[i];
		} else {
			siril_log_message(_("Unknown argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!filename || !norm_arg) {
		siril_log_message(_("Missing operand filename or normalization factor.\n"));
		return CMD_ARG_ERROR;
	}

	/* Parse norm */
	gchar *end;
	float norm = g_ascii_strtod(norm_arg, &end);
	if (end == norm_arg) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), norm_arg);
		return CMD_ARG_ERROR;
	}

	/* Allocate and read operand fits */
	fits *operand_fit = calloc(1, sizeof(fits));
	if (!operand_fit) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	if (readfits(filename, operand_fit, NULL, !com.pref.force_16bit)) {
		free(operand_fit);
		return CMD_INVALID_IMAGE;
	}

	/* Allocate and initialize user data */
	struct fdiv_data *data = calloc(1, sizeof(struct fdiv_data));
	if (!data) {
		clearfits(operand_fit);
		free(operand_fit);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	data->destructor = fdiv_destructor;
	data->operand_fit = operand_fit;
	data->norm = norm;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		fdiv_destructor(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = fdiv_image_hook;
	args->idle_function = NULL;
	args->description = _("Image division");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = fdiv_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

// Structure to hold fmul-specific data
struct fmul_data {
	void (*destructor)(void *);  // Required as first member
	float coeff;
	gboolean from8b;
};

// Image processing hook for fmul
static int fmul_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct fmul_data *data = (struct fmul_data *)args->user;

	// Perform the scalar multiplication
	int retval = soper(fit, data->coeff, OPER_MUL, TRUE);

	if (data->from8b) {
		// Image is now 32b, need to reset slider max and update hi/lo
		invalidate_stats_from_fit(fit);
		image_find_minmax(fit);
		fit->keywords.hi = (WORD)(fit->maxi * USHRT_MAX_SINGLE);
		fit->keywords.lo = (WORD)(fit->mini * USHRT_MAX_SINGLE);
		set_cutoff_sliders_max_values();
	}

	return retval;
}

static gchar *fmul_log_hook(gpointer p, log_hook_detail detail) {
	struct fmul_data *args = (struct fmul_data*) p;
	return g_strdup_printf(_("Scalar multiplication: factor %.6f"), args->coeff);
}

// Main command function for fmul
int process_fmul(int nb) {
	gchar *end;
	gboolean mask_aware = FALSE;

	float coeff = g_ascii_strtod(word[1], &end);
	if (end == word[1] || coeff <= 0.f) {
		siril_log_message(_("Multiplying by a coefficient less than or equal to 0 is not possible.\n"));
		return CMD_ARG_ERROR;
	}

	if (word[2] && !g_strcmp0(word[2], "-mask")) {
		mask_aware = TRUE;
	}
	// Allocate and initialize user data
	struct fmul_data *data = calloc(1, sizeof(struct fmul_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->coeff = coeff;
	data->from8b = (gfit->orig_bitpix == BYTE_IMG);

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;  // soper works in-place
	args->image_hook = fmul_image_hook;
	args->idle_function = NULL;  // Use default
	args->description = _("Scalar multiplication");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;  // This command modifies gfit
	args->command = TRUE; // calling as command, not from GUI
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = fmul_log_hook;
	args->max_threads = 1;  // soper doesn't need multi-threading
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Start the worker thread
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

// Structure to hold gauss-specific data
struct gauss_data {
	void (*destructor)(void *);  // Required as first member
	double sigma;
};

// Image processing hook for gauss
static int gauss_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct gauss_data *data = (struct gauss_data *)args->user;

	// Perform the Gaussian blur
	int retval = unsharp(fit, data->sigma, 0.0, TRUE);

	// Add to history
	char log[90];
	sprintf(log, "Gaussian filtering, sigma: %.2f", data->sigma);
	fit->history = g_slist_append(fit->history, strdup(log));

	return retval;
}

static gchar *gauss_log_hook(gpointer p, log_hook_detail detail) {
	struct gauss_data *data = (struct gauss_data*) p;
	return g_strdup_printf(_("Gaussian blur: sigma = %f"), data->sigma);
}

// Main command function for gauss
int process_gauss(int nb) {
	gchar *end;
	gboolean mask_aware = FALSE;
	double sigma = g_ascii_strtod(word[1], &end);
	if (end == word[1] || sigma <= 0.0) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	if (word[2] && !g_strcmp0(word[2], "-mask")) {
		mask_aware = TRUE;
	}

	image_cfa_warning_check();

	// Allocate and initialize user data
	struct gauss_data *data = calloc(1, sizeof(struct gauss_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->sigma = sigma;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;  // Gaussian blur needs temporary buffers
	args->image_hook = gauss_image_hook;
	args->idle_function = NULL;  // Use default
	args->description = _("Gaussian blur");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;  // This command modifies gfit
	args->command = TRUE; // calling as command, not from GUI
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = gauss_log_hook;
	args->max_threads = com.max_thread;  // Gaussian blur can benefit from multi-threading
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Start the worker thread
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

struct entropy_data {
	void (*destructor)(void *);
	rectangle area;
	gboolean has_selection;
	float total_entropy;
};

static int entropy_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct entropy_data *data = (struct entropy_data *)args->user;
	float e = 0.f;

	image_cfa_warning_check();

	if (data->has_selection) {
		rectangle area = data->area;
		for (int c = 0; c < fit->naxes[2]; c++)
			e += entropy(fit, c, &area, NULL);
	} else {
		for (int c = 0; c < fit->naxes[2]; c++)
			e += entropy(fit, c, NULL, NULL);
	}

	data->total_entropy = e;
	return 0;
}

static gchar *entropy_log_hook(gpointer p, log_hook_detail detail) {
	struct entropy_data *args = (struct entropy_data*) p;
	if (detail == DETAILED) {
		return g_strdup_printf(_("Entropy: %.3f"), args->total_entropy);
	}
	return g_strdup(_("Entropy computation"));
}

int process_entropy(int nb) {
	// Allocate and initialize user data
	struct entropy_data *data = calloc(1, sizeof(struct entropy_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->has_selection = (com.selection.w > 0 && com.selection.h > 0);
	if (data->has_selection) {
		data->area = com.selection;
	}
	data->total_entropy = 0.f;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = entropy_image_hook;
	args->idle_function = NULL;
	args->description = _("Entropy");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;  // This doesn't modify gfit
	args->command = TRUE;
	args->user = data;
	args->log_hook = entropy_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_unpurple(int nb) {
	gboolean mask_aware = FALSE;
	gboolean withstarmask = FALSE;
	double mod = 1.0, thresh = 0.0;

	/* Scan arguments (order-independent) */
	for (int i = 1; i < nb; i++) {
		gchar *arg = word[i];

		if (!g_strcmp0(arg, "-mask")) {
			mask_aware = TRUE;
		}
		else if (!g_strcmp0(arg, "-starmask")) {
			withstarmask = TRUE;
		}
		else if (g_str_has_prefix(arg, "-mod=")) {
			gchar *val = arg + 5;
			gchar *end;
			mod = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-thresh=")) {
			gchar *val = arg + 8;
			gchar *end;
			thresh = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_color_message(_("Unknown argument %s\n"), "red", arg);
			return CMD_ARG_ERROR;
		}
	}

	/* Allocate parameters */
	struct unpurpleargs *params = new_unpurple_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_GENERIC_ERROR;
	}

	params->mod_b = mod;
	params->thresh = thresh;
	params->withstarmask = withstarmask;
	params->fit = gfit;
	params->verbose = TRUE;
	params->applying = TRUE;

	/* Set up starmask if needed */
	params->starmask_needs_freeing = FALSE;
	if (withstarmask) {
		params->starmask = calloc(1, sizeof(fits));
		if (!params->starmask) {
			PRINT_ALLOC_ERR;
			free_unpurple_args(params);
			return CMD_GENERIC_ERROR;
		}

		if (generate_binary_starmask(gfit, &params->starmask, thresh)) {
			free_unpurple_args(params);
			return CMD_GENERIC_ERROR;
		}
		params->starmask_needs_freeing = TRUE;
	} else {
		params->starmask = NULL;
	}

	/* Allocate worker args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_unpurple_args(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = unpurple_image_hook;
	args->idle_function = NULL;
	args->description = _("Unpurple Filter");
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = unpurple_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->mask_aware = mask_aware;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	siril_log_message(_("Unpurple mod: %.2f, threshold: %.2f, withstarmask: %d, mask_aware: %d\n"),
		mod, thresh, withstarmask, mask_aware);

	return CMD_OK;
}

int process_epf(int nb) {
	gchar *end;
	double	d = 0.0, mod = 1.0,
			sigma_col = 11.0,
			sigma_space = 11.0;
	ep_filter_t filter = EP_BILATERAL;
	gchar *filename = NULL;
	fits *guidefit = NULL;
	gboolean guide_needs_freeing = FALSE;
	gboolean mask_aware = FALSE;

	for (int i = 1 ; i < nb ; i++) {
		gchar *arg = word[i];
		if (!g_strcmp0(arg, "-guided")) {
			filter = EP_GUIDED;
		}
		else if (!g_strcmp0(arg, "-mask")) {
			mask_aware = TRUE;
		}
		else if (g_str_has_prefix(arg, "-guideimage=")) {
			gchar *val = arg + 12;
			if (filename) {
				g_free(filename);
			}
			filename = g_strdup(val);
		}
		else if (g_str_has_prefix(arg, "-d=")) {
			gchar *val = arg + 3;
			d = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				g_free(filename);
				return CMD_ARG_ERROR;
			}
			if (d > 25)
				siril_log_color_message(_("Warning: d > approx. 25 may result in lengthy execution times.\n"), "salmon");
		}
		else if (g_str_has_prefix(arg, "-mod=")) {
			gchar *val = arg + 5;
			mod = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				g_free(filename);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-si=")) {
			gchar *val = arg + 4;
			sigma_col = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				g_free(filename);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-ss=")) {
			gchar *val = arg + 4;
			sigma_space = g_ascii_strtod(val, &end);
			if (end == val) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", arg);
				g_free(filename);
				return CMD_ARG_ERROR;
			}
		}
	}
	if (sigma_col <= 0.0 || sigma_col > 65535) {
		 siril_log_color_message(_("-si= value must be > 0.0 and <= 65535\n"), "red");
		 g_free(filename);
		 return CMD_ARG_ERROR;
	}
	if (d < 0.0 || d > 20.0) {
		 siril_log_color_message(_("-d= value must be > 0.0 and <= 20.0\n"), "red");
		 g_free(filename);
		 return CMD_ARG_ERROR;
	}
	if (sigma_space <= 0.0 || sigma_space > 32.0) {
		 siril_log_color_message(_("-ss= value must be > 0.0 and <= 32.0\n"), "red");
		 g_free(filename);
		 return CMD_ARG_ERROR;
	}
	if (sigma_space > 20.0 && d == 0)
			siril_log_color_message(_("Warning: spatial sigma > approx. 20 with auto diameter may result in lengthy execution times.\n"), "salmon");

	if (mod <= 0.0 || mod > 1.0) {
		 siril_log_color_message(_("-mod= value must be > 0.0 and <= 1.0\n"), "red");
		 g_free(filename);
		 return CMD_ARG_ERROR;
	}
	if (filename != NULL)
		filter = EP_GUIDED; // passing guideimage name is enough to set to guided
	if (filter == EP_GUIDED) {
		if (filename) {
			guidefit = calloc(1, sizeof(fits));
			if (readfits(filename, guidefit, NULL, FALSE)) {
				siril_log_color_message(_("Error: guide image could not be loaded\n"), "red");
				clearfits(guidefit);
				free(guidefit);
				g_free(filename);
				return CMD_ARG_ERROR;
			}
			g_free(filename);
			guide_needs_freeing = TRUE;
		} else {
			guidefit = gfit;
		}
		if (guidefit->rx != gfit->rx || guidefit->ry != gfit->ry) {
			siril_log_color_message(_("Error: guide image dimensions do not match\n"), "red");
			if (guide_needs_freeing) {
				clearfits(guidefit);
				free(guidefit);
			}
			return CMD_ARG_ERROR;
		}
	}
	if (filter == EP_GUIDED && d == 0.0) {
		siril_log_color_message(_("Warning: d = 0.0 cannot be used to specify automatic diameter when using a guided filter. Setting d to default value of 5.\n"), "salmon");
		d = 5.0;
	}

	// Allocate parameters using the allocator
	struct epfargs *params = new_epf_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		if (guide_needs_freeing) {
			clearfits(guidefit);
			free(guidefit);
		}
		return CMD_GENERIC_ERROR;
	}

	params->fit = gfit;
	params->guidefit = guidefit;
	params->guide_needs_freeing = guide_needs_freeing;
	params->d = d;
	params->sigma_col = sigma_col;
	params->sigma_space = sigma_space;
	params->mod = mod;
	params->filter = filter;
	params->verbose = TRUE;
	params->applying = TRUE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_epf_args(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 3.0f;
	args->image_hook = epf_image_hook;
	args->idle_function = NULL; // Use default idle function for command-line
	args->description = _("Edge Preserving Filter");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->user = params;
	args->mask_aware = mask_aware;
	args->log_hook = epf_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	char log[90];
	if (filter == EP_BILATERAL) {
		sprintf(log, "Bilateral filtering, d: %.2f, sigma(color): %.2f, sigma(spatial): %.2f, modulation: %.2f", d, sigma_col, sigma_space, mod);
	} else {
		sprintf(log, "Guided filtering, d: %.2f, sigma: %.2f, modulation: %.2f", d, sigma_col, mod);
	}
	gfit->history = g_slist_append(gfit->history, strdup(log));

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_getref(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		siril_log_message(_("Error: cannot open sequence\n"));
		return CMD_SEQUENCE_NOT_FOUND;
	}

	int ref_image = seq->reference_image;
	if (seq->reference_image < 0) {
		siril_log_message(_("Reference image is undefined, the following would be used:\n"));
		ref_image = sequence_find_refimage(seq);
	}

	if (seq->type == SEQ_REGULAR) {
		long maxpath = get_pathmax();
		char filename[maxpath];
		fit_sequence_get_image_filename_checkext(seq, ref_image, filename);
		siril_log_message(_("Image %d: '%s'\n"), ref_image, filename);
	}
	else siril_log_message(_("Image %d\n"), ref_image);

	if (!seq->imgparam[ref_image].incl)
		siril_log_message(_("Warning: this image is excluded from the sequence main processing list\n"));
	notify_gfit_modified();
	return CMD_OK;
}

static int grey_flat_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	compute_grey_flat(fit);
	return 0;
}

int process_grey_flat(int nb) {
	if (isrgb(gfit)) {
		return CMD_FOR_CFA_IMAGE;
	}

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = grey_flat_image_hook;
	args->description = _("Grey flat");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

/* Wrapper hook for command-line PSF estimation with dynamic estk_data */
int estimate_only_cmd_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct estk_data *data = (struct estk_data *)args->user;
	if (!data)
		return 1;
	return GPOINTER_TO_INT(estimate_only(data));
}

/* Idle function for command-line PSF estimation */
gboolean estimate_only_cmd_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	// PSF estimation doesn't modify gfit directly, but may update the kernel
	// No need to notify_gfit_modified() here

	// Free using the generic cleanup which will call the destructor
	free_generic_img_args(args);

	stop_processing_thread();
	return FALSE;
}

int process_makepsf(int nb) {
	gboolean error = FALSE;
	// Use the new allocator to ensure function pointers (destructor) are set
	estk_data* data = alloc_estk_data();
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_GENERIC_ERROR;
	}

	reset_conv_args(data);
	data->fit = gfit; // Explicitly enforce gfit per instructions

	cmd_errors status = CMD_OK;

	char *arg_1 = word[1];
	if (!g_strcmp0(arg_1, "clear")) {
		if (get_thread_run()) {
			siril_log_message(_("Error: will not clear the PSF while a sequence is running.\n"));
			status = CMD_GENERIC_ERROR;
			goto terminate_makepsf;
		}
		reset_conv_kernel();
		siril_log_color_message(_("Deconvolution kernel cleared.\n"), "green");
		goto terminate_makepsf;
	} else {
		if (!g_strcmp0(arg_1, "save")) {
			siril_log_message(_("Save PSF to file:\n"));
			if (!word[2] || word[2][0] == '\0') {
				// We cannot call on_bdeconv_savekernel_clicked here as it now reads from GUI widgets.
				// Enforce filename argument in command line mode.
				siril_log_color_message(_("Error: a filename must be specified when saving a PSF from the command line.\n"), "red");
				status = CMD_ARG_ERROR;
			} else {
				if (!(g_str_has_suffix(word[2], ".fit") || g_str_has_suffix(word[2], ".fits") || g_str_has_suffix(word[2], ".fts") || g_str_has_suffix(word[2], ".tif"))) {
					siril_log_color_message(_("Error: filename must have the extension \".fit\", \".fits\", \".fts\" or \".tif\"\n"), "red");
					status = CMD_ARG_ERROR;
					goto terminate_makepsf;
				}
				// Update: pass data to save_kernel
				save_kernel(word[2], data);
			}
			goto terminate_makepsf;
		}

		// If we aren't loading, we are generating, so reset previous
		if (g_strcmp0(arg_1, "load") != 0) {
			reset_conv_kernel();
		}

		status = CMD_ARG_ERROR; // setting to this value as it will be the most likely error from now on

		if (!g_strcmp0(arg_1, "blind")) {
			if (!(single_image_is_loaded() || sequence_is_loaded())) {
				siril_log_message(_("Error: image or sequence must be loaded to carry out %s PSF estimation, aborting...\n"), "blind");
				status = CMD_GENERIC_ERROR;
				goto terminate_makepsf;
			}
			data->psftype = PSF_BLIND;
			siril_log_message(_("Blind kernel estimation:\n"));
			for (int i = 2; i < nb; i++) {
				char *arg = word[i], *end;
				if (!word[i])
					break;
				if (!g_strcmp0(arg, "-l0")) {
					siril_log_message(_("ℓ0 descent prior method\n"));
					data->blindtype = 1;
				}
				else if (!g_strcmp0(arg, "-si")) {
					siril_log_message(_("spectral irregularity method\n"));
					data->blindtype = 0;
				}
				else if (!g_strcmp0(arg, "-multiscale")) {
					siril_log_message(_("multiscale estimation\n"));
					data->multiscale = TRUE;
				}
				else if (g_str_has_prefix(arg, "-lambda=")) {
					arg += 8;
					float lambda = (float) g_ascii_strtod(arg, &end);
					if (arg == end)
						goto terminate_makepsf;
					if (error || lambda < 0.f || lambda > 100000.f) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "lambda", "]0,10000[");
						goto terminate_makepsf;
					}
					data->lambda = lambda;
					data->finaldeconvolutionweight = lambda;
					data->intermediatedeconvolutionweight = lambda;
				}
				else if (g_str_has_prefix(arg, "-comp=")) {
					arg += 6;
					float comp = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if (comp < 1.f || comp > 100000.f) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "compensation factor", "]1,10000[");
						goto terminate_makepsf;
					}
					data->compensationfactor = comp;
				}
				else if (g_str_has_prefix(arg, "-ks=")) {
					arg += 4;
					int ks = (int) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if (ks < 3 || !(ks %2) || ks > min(gfit->rx, gfit->ry)) {
						siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width), aborting.\n"));
						goto terminate_makepsf;
					}
					data->ks = ks;
				}
				else if (g_str_has_prefix(arg, "-savepsf=")) {
					if (data->savepsf_filename) {
						g_free(data->savepsf_filename);
						data->savepsf_filename = NULL;
					}
					arg += 9;
					if (arg[0] == '\0') {
						siril_log_message(_("Error: no filename specified, aborting.\n"));
						goto terminate_makepsf;
					} else {
						if (!(g_str_has_suffix(arg, ".fit") || g_str_has_suffix(arg, ".fits") || g_str_has_suffix(arg, ".fts") || g_str_has_suffix(arg, ".tif"))) {
							siril_log_color_message(_("Error: filename must have the extension \".fit\", \".fits\", \".fts\" or \".tif\"\n"), "red");
							goto terminate_makepsf;
						}
						data->savepsf_filename = g_strdup(arg);
						data->save_after = TRUE;
					}
				} else {
					siril_log_color_message(_("Unknown parameter %s, aborting.\n"), "red", arg);
					goto terminate_makepsf;
				}
			}
			image_cfa_warning_check();

			// Allocate generic worker args
			struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				PRINT_ALLOC_ERR;
				if (data->destroy_fn) data->destroy_fn(data);
				return CMD_GENERIC_ERROR;
			}

			args->fit = gfit;
			args->mem_ratio = 3.0f; // PSF estimation memory requirement
			args->image_hook = estimate_only_image_hook;
			args->idle_function = NULL;
			args->description = _("PSF Estimation (Blind)");
			args->verbose = TRUE;
			args->user = data; // Dynamic estk_data managed by generic worker
			args->max_threads = com.max_thread;
			args->for_preview = FALSE;
			args->for_roi = FALSE;

			if (!start_in_new_thread(generic_image_worker, args)) {
				// Free both args and data on thread start failure
				free_generic_img_args(args);
				return CMD_GENERIC_ERROR;
			}
			return CMD_OK;
		} else if (!g_strcmp0(arg_1, "stars")) {
			data->recalc_ks = FALSE;
			gboolean force_ks = FALSE;
			if (!(single_image_is_loaded() || sequence_is_loaded())) {
				siril_log_message(_("Error: image or sequence must be loaded to carry out %s PSF estimation, aborting...\n"), "star-based");
				status = CMD_GENERIC_ERROR;
				goto terminate_makepsf;
			}
			data->psftype = PSF_STARS;
			for (int i = 2; i < nb; i++) {
				char *arg = word[i], *end;
				if (!word[i])
					break;
				if (!g_strcmp0(arg, "-sym")) {
					siril_log_message(_("symmetric kernel\n"));
					data->symkern = TRUE;
				}
				else if (g_str_has_prefix(arg, "-ks=")) {
					arg += 4;
					int ks = (int) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if (ks < 3 || !(ks %2) || ks > min(gfit->rx, gfit->ry)) {
						siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width), aborting.\n"));
						goto terminate_makepsf;
					}
					data->ks = ks;
					force_ks = TRUE;
				}
				else if (g_str_has_prefix(arg, "-savepsf=")) {
					if (data->savepsf_filename) {
						g_free(data->savepsf_filename);
						data->savepsf_filename = NULL;
					}
					arg += 9;
					if (arg[0] == '\0') {
						siril_log_message(_("Error: no filename specified, aborting.\n"));
						goto terminate_makepsf;
					} else {
						if (!(g_str_has_suffix(arg, ".fit") || g_str_has_suffix(arg, ".fits") || g_str_has_suffix(arg, ".fts") || g_str_has_suffix(arg, ".tif"))) {
							siril_log_color_message(_("Error: filename must have the extension \".fit\", \".fits\", \".fts\" or \".tif\"\n"), "red");
							goto terminate_makepsf;
						}
						data->savepsf_filename = g_strdup(arg);
						data->save_after = TRUE;
					}
				} else {
					siril_log_color_message(_("Unknown parameter %s, aborting.\n"), "red", arg);
					goto terminate_makepsf;
				}
			}
			if (!force_ks) {
				data->recalc_ks = TRUE;
			}
			image_cfa_warning_check();

			// Allocate generic worker args
			struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				PRINT_ALLOC_ERR;
				if (data->destroy_fn) data->destroy_fn(data);
				return CMD_GENERIC_ERROR;
			}

			args->fit = gfit;
			args->mem_ratio = 3.0f; // PSF estimation memory requirement
			args->image_hook = estimate_only_image_hook;
			args->idle_function = NULL;
			args->description = _("PSF Estimation (Stars)");
			args->verbose = TRUE;
			args->user = data; // Dynamic estk_data managed by generic worker
			args->max_threads = com.max_thread;
			args->for_preview = FALSE;
			args->for_roi = FALSE;

			if (!start_in_new_thread(generic_image_worker, args)) {
				// Free both args and data on thread start failure
				free_generic_img_args(args);
				return CMD_GENERIC_ERROR;
			}
			return CMD_OK;
		} else if (!g_strcmp0(arg_1, "manual")) {
			siril_log_message(_("Manual PSF generation:\n"));
			data->psftype = PSF_MANUAL;
			for (int i = 2; i < nb; i++) {
				char *arg = word[i], *end;
				if (!word[i])
					break;
				if (!g_strcmp0(arg, "-gaussian")) {
					siril_log_message(_("Gaussian PSF\n"));
					data->profile = 0;
				}
				else if (!g_strcmp0(arg, "-moffat")) {
					siril_log_message(_("Moffat PSF\n"));
					data->profile = 1;
				}
				else if (!g_strcmp0(arg, "-disc")) {
					siril_log_message(_("Disc PSF\n"));
					data->profile = 2;
				}
				else if (!g_strcmp0(arg, "-airy")) {
					siril_log_message(_("Airy disc PSF\n"));
					data->profile = 3;
				}
				else if (g_str_has_prefix(arg, "-fwhm=")) {
					arg += 6;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val <= 0.f) || (val > 100.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "fwhm", "]0,100]");
						goto terminate_makepsf;
					}
					data->psf_fwhm = val;
				}
				else if (g_str_has_prefix(arg, "-angle=")) {
					arg += 7;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val <= -360.f) || (val > 360.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "angle", "]-360,360]");
						goto terminate_makepsf;
					}
					data->psf_angle = val;
				}
				else if (g_str_has_prefix(arg, "-ratio=")) {
					arg += 7;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					else if ((val < 1.f) || (val > 5.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "ratio", "[1,5]");
						goto terminate_makepsf;
					}
					data->psf_ratio = val;
				}
				else if (g_str_has_prefix(arg, "-beta=")) {
					arg += 6;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val <= 0.f) || (val > 10.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "beta", "]0,10]");
						goto terminate_makepsf;;
					}
					data->psf_beta = val;
				}
				else if (g_str_has_prefix(arg, "-dia=")) {
					arg += 5;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val <= 0.f) || (val > 5000.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "dia", "]0,5000]");
						goto terminate_makepsf;
					}
					data->airy_diameter = val;
				}
				else if (g_str_has_prefix(arg, "-fl=")) {
					arg += 4;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val <= 0.f) || (val > 60000.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "fl", "]0,60000]");
						goto terminate_makepsf;
					}
					data->airy_fl = val;
				}
				else if (g_str_has_prefix(arg, "-wl=")) {
					arg += 4;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val < 100.f) || (val > 30000.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "wl", "[100,30000]");
						goto terminate_makepsf;
					}
					data->airy_wl = val;
				}
				else if (g_str_has_prefix(arg, "-pixelsize=")) {
					arg += 11;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if ((val < 1.f) || (val > 30.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "pixelsize", "[1,30]");
						goto terminate_makepsf;
					}
					data->airy_pixelsize = val;
				}
				else if (g_str_has_prefix(arg, "-obstruct=")) {
					arg += 10;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val < 0.f) || (val >= 100.f)) {
						siril_log_message(_("Error in %s parameter: must be in %s, aborting.\n"), "obstruct", "[0,100[");
						goto terminate_makepsf;
					}
					data->airy_obstruction = val;
				}
				else if (g_str_has_prefix(arg, "-ks=")) {
					arg += 4;
					int ks = (int) g_ascii_strtod(arg, &end);
					if (arg == end) goto terminate_makepsf;
					if (ks < 3 || !(ks %2) || ks > min(gfit->rx, gfit->ry)) {
						siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width), aborting.\n"));
						goto terminate_makepsf;
					}
					data->ks = ks;
				}
				else if (g_str_has_prefix(arg, "-savepsf=")) {
					if (data->savepsf_filename) {
						g_free(data->savepsf_filename);
						data->savepsf_filename = NULL;
					}
					arg += 9;
					if (arg[0] == '\0') {
						siril_log_message(_("Error: no filename specified, aborting.\n"));
						goto terminate_makepsf;
					} else {
						if (!(g_str_has_suffix(arg, ".fit") || g_str_has_suffix(arg, ".fits") || g_str_has_suffix(arg, ".fts") || g_str_has_suffix(arg, ".tif"))) {
							siril_log_color_message(_("Error: filename must have the extension \".fit\", \".fits\", \".fts\" or \".tif\"\n"), "red");
							goto terminate_makepsf;
						}
						data->savepsf_filename = g_strdup(arg);
						data->save_after = TRUE;
					}
				} else {
					siril_log_color_message(_("Unknown parameter %s, aborting.\n"), "red", arg);
					goto terminate_makepsf;
				}
			}

			// Allocate generic worker args
			struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				PRINT_ALLOC_ERR;
				if (data->destroy_fn) data->destroy_fn(data);
				return CMD_GENERIC_ERROR;
			}

			args->fit = gfit;
			args->mem_ratio = 3.0f; // PSF estimation memory requirement
			args->image_hook = estimate_only_image_hook;
			args->idle_function = NULL;
			args->description = _("PSF Estimation (Manual)");
			args->verbose = TRUE;
			args->user = data; // Dynamic estk_data managed by generic worker
			args->max_threads = com.max_thread;
			args->for_preview = FALSE;
			args->for_roi = FALSE;

			if (!start_in_new_thread(generic_image_worker, args)) {
				// Free both args and data on thread start failure
				free_generic_img_args(args);
				return CMD_GENERIC_ERROR;
			}
			return CMD_OK;
		} else if (!g_strcmp0(arg_1, "load")) {
			siril_log_message(_("Load PSF from file:\n"));
			if (word[2] && word[2][0] != '\0') {
				// Update: pass data to load_kernel
				if (load_kernel(word[2], data)) {
					siril_log_color_message(_("Error loading PSF.\n"), "red");
					status = CMD_FILE_NOT_FOUND;
					goto terminate_makepsf;
				}
			}
			data->psftype = PSF_PREVIOUS;
			status = CMD_OK;
			goto terminate_makepsf;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg_1);
		}
	}
terminate_makepsf:
	// If we reached here, the data struct was not handed off to a worker thread.
	// Clean it up.
	if (data->destroy_fn)
		data->destroy_fn(data);
	else {
		// Fallback if destructor missing (shouldn't happen with alloc_estk_data)
		g_free(data->savepsf_filename);
		free(data);
	}
	return status;
}

int parse_deconvolve(int first_arg, int nb, estk_data* data, nonblind_t type) {
	gboolean error = FALSE;
	gboolean kernel_loaded = FALSE;

	reset_conv_args(data);
	data->fit = gfit; // Ensure we operate on the global fit for commands

	data->regtype = REG_NONE_GRAD;

	// If a kernel already exists globally, sync the kernel size to the struct
	if (com.kernel && com.kernelsize > 0)
		data->ks = com.kernelsize;

	if (type == DECONV_SB)
		data->finaliters = 1;
	if (type == DECONV_WIENER)
		data->alpha = 1.f / 500.f;

	for (int i = first_arg; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;

		if (g_str_has_prefix(arg, "-loadpsf=")) {
			arg += 9;
			// UPDATE: pass 'data' to load_kernel
			if (arg[0] != '\0' && load_kernel(arg, data)) {
				siril_log_message(_("Error loading PSF.\n"));
				return CMD_FILE_NOT_FOUND;
			}
			kernel_loaded = TRUE;
		}
		else if (!g_strcmp0(arg, "-mask")) {
			data->mask_aware = TRUE;
		}
		else if (g_str_has_prefix(arg, "-alpha=")) {
			arg += 7;
			float alpha = (float) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((alpha < 0.f) || (alpha > 100000.f)) {
				siril_log_message(_("Error in alpha parameter: must be between 0 and 1e5, aborting.\n"));
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->alpha = 1.f / alpha;
			}
		}
		else if (g_str_has_prefix(arg, "-iters=")) {
			arg += 7;
			float iters = (int) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((iters < 1) || (iters > 100000)) {
				siril_log_message(_("Error in iterations parameter: must be between 1 and 1e5, aborting.\n"));
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->finaliters = iters;
			}
		}
		else if (g_str_has_prefix(arg, "-stop=")) {
			arg += 6;
			float stopcriterion = (float) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((stopcriterion < 0.f) || (stopcriterion > 1.f)) {
				siril_log_message(_("Error in stop parameter: must be between 0 and 1, aborting.\n"));
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->stopcriterion = stopcriterion;
				data->stopcriterion_active = TRUE;
			}
		}
		else if (g_str_has_prefix(arg, "-gdstep=")) {
			arg += 8;
			float stepsize = (float) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((stepsize < 0.f) || (stepsize > 1.f)) {
				siril_log_message(_("Error in step size parameter: must be between 0 and 1, aborting.\n"));
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->stepsize = stepsize;
			}
		}
		else if (!g_strcmp0(arg, "-tv")) {
			data->regtype = REG_TV_GRAD;
		}
		else if (!g_strcmp0(arg, "-mul")) {
			data->rl_method = RL_MULT;
		}
		else if (!g_strcmp0(arg, "-fh")) {
			data->regtype = REG_FH_GRAD;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			return CMD_ARG_ERROR;
		}
	}

	if (error) {
		return CMD_ARG_ERROR;
	}

	// Guess the user's intentions for a kernel:
	// Order of preference is: if a PSF has specifically been loaded, use that, else
	// if com.stars is populated, assume the user wants to make a PSF from stars and
	// use that. If not, if com.kernel is already populated then use it. Otherwise,
	// do a blind deconvolution using the default parameters (and default l0 method).
	// A manual PSF, if required, must be created using the makepsf command and then rl
	// will detect it as an existing kernel and use it.

	if (kernel_loaded)
		data->psftype = PSF_PREVIOUS; // use loaded (existing) kernel
	else if (com.stars && com.stars[0])
		data->psftype = PSF_STARS; // PSF from stars
	else if (com.kernel && (com.kernelsize != 0))
		data->psftype = PSF_PREVIOUS; // use existing kernel
	else
		data->psftype = PSF_BLIND; // blind deconvolve

	data->nonblindtype = type;
	return CMD_OK;
}

int process_deconvolve(int nb, nonblind_t type) {
	// Use the new allocator
	estk_data* data = alloc_estk_data();
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	// Explicitly set gfit
	data->fit = gfit;

	int ret = parse_deconvolve(1, nb, data, type);
	if (ret == CMD_OK) {
		image_cfa_warning_check();

		// Allocate generic worker args
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			PRINT_ALLOC_ERR;
			if (data->destroy_fn)
				data->destroy_fn(data);
			return CMD_ALLOC_ERROR;
		}

		args->fit = gfit;
		args->mem_ratio = 4.0f; // Deconvolution needs significant memory
		args->image_hook = deconvolve_image_hook;
		args->idle_function = NULL;
		args->description = _("Deconvolution");
		args->command_updates_gfit = TRUE;
		args->command = TRUE; // calling as command, not from GUI
		args->verbose = TRUE;
		args->user = data; // Dynamic estk_data managed by generic worker
		args->mask_aware = data->mask_aware;
		args->max_threads = com.max_thread;
		args->for_preview = FALSE;
		args->for_roi = FALSE;

		if (!start_in_new_thread(generic_image_worker, args)) {
			// Free both args and data on thread start failure
			free_generic_img_args(args);
			return CMD_GENERIC_ERROR;
		}
		return CMD_OK;
	}

	// Call destructor on data if parse failed or thread didn't start
	if (data->destroy_fn)
		data->destroy_fn(data);
	return ret;
}

int process_seqdeconvolve(int nb, nonblind_t type) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	estk_data* data = calloc(1, sizeof(estk_data));
	reset_conv_args(data);

	int ret = parse_deconvolve(2, nb, data, type);
	if (ret == CMD_OK){
		sequence_cfa_warning_check(seq);
		deconvolve_sequence_command(data, seq);
		return CMD_OK;
	}

	if (!check_seq_is_comseq(seq))
		free_sequence(seq, TRUE);
	free(data);
	return ret;
}

int process_sb(int nb) {
	return process_deconvolve(nb, DECONV_SB);
}

int process_rl(int nb) {
	return process_deconvolve(nb, DECONV_RL);
}

int process_wiener(int nb) {
	return process_deconvolve(nb, DECONV_WIENER);
}

int process_seq_sb(int nb) {
	return process_seqdeconvolve(nb, DECONV_SB);
}

int process_seq_rl(int nb) {
	return process_seqdeconvolve(nb, DECONV_RL);
}

int process_seq_wiener(int nb) {
	return process_seqdeconvolve(nb, DECONV_WIENER);
}

// Structure to hold unsharp-specific data
struct unsharp_data {
	void (*destructor)(void *);  // Required as first member
	double sigma;
	double multi;
};

// Image processing hook for unsharp
static int unsharp_cmd_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct unsharp_data *data = (struct unsharp_data *)args->user;

	// Perform the unsharp mask operation
	int retval = unsharp(fit, data->sigma, data->multi, TRUE);

	// Add to history
	char log[90];
	sprintf(log, "Unsharp filtering, sigma: %.2f, coefficient: %.2f", data->sigma, data->multi);
	fit->history = g_slist_append(fit->history, strdup(log));

	return retval;
}

static gchar *unsharp_log_hook(gpointer p, log_hook_detail detail) {
	struct unsharp_data *data = (struct unsharp_data *) p;
	return g_strdup_printf(_("Unsharp mask: sigma=%f, amount=%f"), data->sigma, data->multi);
}

// Main command function for unsharp
int process_unsharp(int nb) {
	gchar *end;
	double sigma = g_ascii_strtod(word[1], &end);
	if (end == word[1] || sigma <= 0.0) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	double multi = g_ascii_strtod(word[2], &end);
	if (end == word[2]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[2]);
		return CMD_ARG_ERROR;
	}
	gboolean mask_aware = FALSE;
	if (word[3] && !g_strcmp0(word[3], "-mask")) {
		mask_aware = TRUE;
	}
	// Allocate and initialize user data
	struct unsharp_data *data = calloc(1, sizeof(struct unsharp_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}


	data->sigma = sigma;
	data->multi = multi;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;  // Unsharp mask needs temporary buffers
	args->image_hook = unsharp_cmd_image_hook;
	args->idle_function = NULL;
	args->description = _("Unsharp mask");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = unsharp_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Start the worker thread
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

#define CHECK_KEY_LENGTH(__key__) \
	do { \
		if (strlen(__key__) > 8) { \
			siril_log_color_message(_("The size of the key can't exceed 8 characters.\n"), "red"); \
			g_free(__key__); \
			return CMD_ARG_ERROR; \
		} \
	} while(0)

/**
 * update_key key value [comment]
 * update_key -delete key
 * update_key -modify key newkey
 * update_key -comment comment
 */
int process_update_key(int nb) {
	gchar *key = NULL, *value = NULL, *comment = NULL;

	/* manage options */
	if (word[1][0] == '-') {
		if (!g_strcmp0(word[1], "-delete") && word[2]) {
			key = replace_wide_char(word[2]);
			CHECK_KEY_LENGTH(key);
			updateFITSKeyword(gfit, key, NULL, NULL, NULL, TRUE, FALSE);
		} else if (!g_strcmp0(word[1], "-modify") && word[2] && word[3]) {
			key = replace_wide_char(word[2]);
			CHECK_KEY_LENGTH(key);
			value = replace_wide_char(word[3]);
			updateFITSKeyword(gfit, key, value, NULL, NULL, TRUE, FALSE);
		} else if (!g_strcmp0(word[1], "-comment") && word[2]) {
			comment = replace_wide_char(word[2]);
			updateFITSKeyword(gfit, NULL, NULL, NULL, comment, TRUE, FALSE);
		} else {
			return CMD_ARG_ERROR;
		}

	/* without options */
	} else {
		char valstring[FLEN_VALUE];

		key = replace_wide_char(word[1]);
		CHECK_KEY_LENGTH(key);
		value = replace_wide_char(word[2]);

		process_keyword_string_value(value, valstring, string_has_space(value));

		if (nb == 4)
			comment = replace_wide_char(word[3]);

		updateFITSKeyword(gfit, key, NULL, valstring, comment, TRUE, FALSE);
	}
	gui_function(refresh_keywords_dialog, NULL);

	g_free(key);
	g_free(value);
	g_free(comment);

	return CMD_OK;
}

#define CHECK_KEY_LENGTH_SEQ(__key__, __seq__, __args__) \
	do { \
		if (strlen(__key__) > 8) { \
			siril_log_color_message(_("The size of the key can't exceed 8 characters.\n"), "red"); \
			g_free(__key__); \
			if (!check_seq_is_comseq(__seq__)) \
				free_sequence(__seq__, TRUE); \
			free(__args__); \
			return CMD_ARG_ERROR; \
		} \
	} while(0)

int process_seq_update_key(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	struct keywords_data *args = calloc(1, sizeof(struct keywords_data));
	if (!args) {
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	siril_log_color_message(_("Updating keywords...\n"), "green");

	/* manage options */
	if (word[2][0] == '-') {
		if (!g_strcmp0(word[2], "-delete") && word[3]) {
			args->FITS_key = replace_wide_char(word[3]);
			CHECK_KEY_LENGTH_SEQ(args->FITS_key, seq, args);
		} else if (!g_strcmp0(word[2], "-modify") && word[3] && word[4]) {
			args->FITS_key = replace_wide_char(word[3]);
			CHECK_KEY_LENGTH_SEQ(args->FITS_key, seq, args);
			args->newkey = replace_wide_char(word[4]);
		} else if (!g_strcmp0(word[2], "-comment") && word[3]) {
			args->comment = replace_wide_char(word[3]);
		} else {
			free(args);
			return CMD_ARG_ERROR;
		}
	/* without options */
	} else {
		char valstring[FLEN_VALUE];

		args->FITS_key = replace_wide_char(word[2]);
		CHECK_KEY_LENGTH_SEQ(args->FITS_key, seq, args);

		args->value = replace_wide_char(word[3]);
		process_keyword_string_value(args->value, valstring, string_has_space(args->value));

		g_free(args->value);
		args->value = g_strdup(valstring);

		if (nb == 5)
			args->comment = replace_wide_char(word[4]);
	}
	start_sequence_keywords(seq, args);

	return CMD_OK;
}

/* Command interpreter function for CCM using generic_image_worker */
int process_ccm(int nb) {
	sequence *seq = NULL;
	char *prefix = NULL;
	int arg_index = 1, offset;
	gboolean is_sequence = (word[0][2] == 'q');

	if (is_sequence) {
		arg_index = 2;
		offset = 1;
		seq = load_sequence(word[1], NULL);
		if (!seq) {
			return CMD_SEQUENCE_NOT_FOUND;
		}
	} else {
		offset = 0;
		if (!single_image_is_loaded())
			return CMD_IMAGE_NOT_FOUND;
	}

	arg_index++;

	// Parse prefix argument for sequences
	while (arg_index < nb && word[arg_index]) {
		char *arg = word[arg_index];
		if (is_sequence && g_str_has_prefix(arg, "-prefix=")) {
			if (prefix) {
				siril_log_message(_("There can be only one prefix argument"));
				free(prefix);
				return CMD_ARG_ERROR;
			}
			char *value = arg + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
			prefix = strdup(value);
		}
		arg_index++;
	}

	// Allocate ccm_data structure
	struct ccm_data *args = new_ccm_data();
	if (!args) {
		PRINT_ALLOC_ERR;
		free(prefix);
		return CMD_GENERIC_ERROR;
	}

	// Parse matrix elements
	gchar *end;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int word_index = 3 * i + j + 1 + offset;
			args->matrix[i][j] = g_ascii_strtod(word[word_index], &end);
			if (end == word[word_index]) {
				siril_log_message(_("Invalid matrix element (%d, %d) %s, aborting.\n"), i, j, word[word_index]);
				free(prefix);
				free_ccm_data(args);
				return CMD_ARG_ERROR;
			}
		}
	}

	// Parse power parameter (optional)
	if (word[10 + offset]) {
		args->power = g_ascii_strtod(word[10 + offset], &end);
		if (end == word[10 + offset] || args->power < 0.f || args->power > 10.f) {
			siril_log_message(_("Invalid power %s, must be between 0.0 and 10.0: aborting.\n"), word[10 + offset]);
			free(prefix);
			free_ccm_data(args);
			return CMD_ARG_ERROR;
		}
	}

	if (is_sequence) {
		// Sequence processing
		args->seq = seq;
		args->seqEntry = prefix ? prefix : strdup("ccm_");
		apply_ccm_to_sequence(args);
		return CMD_OK;
	} else {
		// Single image processing
		if (!isrgb(gfit)) {
			siril_log_color_message(_("Color Conversion Matrices can only be applied to 3-channel images.\n"), "red");
			free(prefix);
			free_ccm_data(args);
			return CMD_INVALID_IMAGE;
		}

		// Allocate worker args
		struct generic_img_args *worker_args = calloc(1, sizeof(struct generic_img_args));
		if (!worker_args) {
			PRINT_ALLOC_ERR;
			free(prefix);
			free_ccm_data(args);
			return CMD_GENERIC_ERROR;
		}

		worker_args->fit = gfit;
		worker_args->mem_ratio = 1.5f;
		worker_args->image_hook = ccm_single_image_hook;
		worker_args->idle_function = NULL; // Use default idle for commands
		worker_args->description = _("Color Conversion Matrix");
		worker_args->verbose = TRUE;
		worker_args->user = args;
		worker_args->log_hook = ccm_log_hook;
		worker_args->max_threads = com.max_thread;
		worker_args->command = TRUE;  // This is being called from a command
		worker_args->command_updates_gfit = TRUE;  // We need gfit to be updated

		// Build history log
		char log[256];
		snprintf(log, 255, "Color correction matrix applied:");
		gfit->history = g_slist_append(gfit->history, strdup(log));
		snprintf(log, 255, "[ [%.4f %.4f %.4f ] [%.4f %.4f %.4f] [%.4f %.4f %.4f ] ]",
					args->matrix[0][0], args->matrix[0][1], args->matrix[0][2],
					args->matrix[1][0], args->matrix[1][1], args->matrix[1][2],
					args->matrix[2][0], args->matrix[2][1], args->matrix[2][2]);
		gfit->history = g_slist_append(gfit->history, strdup(log));
		snprintf(log, 255, "Power: %.4f", args->power);
		gfit->history = g_slist_append(gfit->history, strdup(log));

		free(prefix);

		if (!start_in_new_thread(generic_image_worker, worker_args)) {
			free_generic_img_args(worker_args);
			return CMD_GENERIC_ERROR;
		}

		// Note: We do NOT return CMD_NOTIFY_GFIT_MODIFIED here because
		// command_updates_gfit is set to TRUE, which handles the notification
		return CMD_OK;
	}
}

int process_cd(int nb) {
	long maxpath = get_pathmax();
	char filename[maxpath];
	int retval;

	g_strlcpy(filename, word[1], maxpath - 1);

	expand_home_in_filename(filename, maxpath);
	retval = siril_change_dir(filename, NULL);
	if (!retval){                   /* chdir ok */
		if (!com.script) {          /* if not a script, change GUI WD*/
			if (com.pref.wd)
				g_free(com.pref.wd);
			com.pref.wd = g_strdup(com.wd);
			writeinitfile();
		}
		gui_function(set_GUI_CWD, NULL);
	}
	else {   /* chdir failed */
	/*
	 *  if filename does not exist, siril_change_dir(filename, NULL) returns 2,
	 *  ie CMD_NO_WAIT which is "no error"; in a script such an error should be
	 *  considered fatal and end the script.
	 *  So we return CMD_DIR_NOT_FOUND instead:
	 */
		retval = CMD_DIR_NOT_FOUND;
	}
	return retval;
}

int process_wrecons(int nb) {

	float coef[7];
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];

	int nb_chan = gfit->naxes[2];

	g_assert(nb_chan == 1 || nb_chan == 3);

	const char *tmpdir = g_get_tmp_dir();

	for (int i = 0; i < nb - 1; ++i) {
		gchar *end;
		coef[i] = g_ascii_strtod(word[i + 1], &end);
		if (end == word[i + 1]) {
			siril_log_message(_("Wrong parameters.\n"));
			return CMD_ARG_ERROR;
		}
	}

	for (int i = 0; i < nb_chan; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		if (gfit->type == DATA_USHORT) {
			wavelet_reconstruct_file(dir[i], coef, gfit->pdata[i]);
		} else if (gfit->type == DATA_FLOAT) {
			wavelet_reconstruct_file_float(dir[i], coef, gfit->fpdata[i]);
		}
		else return CMD_GENERIC_ERROR;
		g_free(dir[i]);
	}
	siril_log_message(_("Wavelet reconstruction\n"));
	notify_gfit_modified();
	return CMD_OK | CMD_NOTIFY_GFIT_MODIFIED;
}

int process_ght_args(int nb, gboolean ght_seq, gboolean *mask_aware, int stretchtype, ght_params *params, struct ght_data *seqdata) {
	int stretch_colourmodel = COL_INDEP;
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
	clip_mode_t clip_mode = RGBBLEND;
	*mask_aware = FALSE;
	double D = NAN, B = 0.0 , LP = 0.0, SP = 0.0, HP = 1.0, BP = NAN;

	for (int i = (ght_seq ? 2 : 1) ; i < nb ; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (!strcmp(arg, "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		else if (!strcmp(arg, "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		else if (!strcmp(arg, "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		else if (!strcmp(arg, "RG")) {
			do_blue = FALSE;
		}
		else if (!strcmp(arg, "RB")) {
			do_green = FALSE;
		}
		else if (!strcmp(arg, "GB")) {
			do_red = FALSE;
		}
		else if (!g_strcmp0(arg, "-mask")) {
			*mask_aware = TRUE;
		}
		else if (g_str_has_prefix(word[i], "-clipmode=")) {
			char *argument = word[i] + 10;
			if (!g_ascii_strncasecmp(argument, "clip", 4))
				clip_mode = CLIP;
			else if (!g_ascii_strncasecmp(argument, "rescale", 7))
				clip_mode = RESCALE;
			else if (!g_ascii_strncasecmp(argument, "globalrescale", 13))
				clip_mode = RESCALEGLOBAL;
			else if (!g_ascii_strncasecmp(argument, "rgbblend", 8))
				clip_mode = RGBBLEND;
			else {
				siril_log_color_message(_("Error, unknown clip mode %s specified\n"), "red", argument);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[i], "-prefix=")) {
			if (ght_seq) {
				char *current = word[i], *value;
				value = current + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
				seqdata->seqEntry = strdup(value);
			} else {
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg,"-BP=")) {
			if (stretchtype == STRETCH_LINEAR) {
				arg += 4;
				BP = g_ascii_strtod(arg, &end);
			} else {
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-sat")) {
			stretch_colourmodel = COL_SAT;
		}
		else {
			if (stretchtype == STRETCH_LINEAR) {
				return CMD_ARG_ERROR;
			}
			else if (g_str_has_prefix(arg, "-human")) {
				stretch_colourmodel = COL_HUMANLUM;
			}
			else if (g_str_has_prefix(arg, "-even")) {
				stretch_colourmodel = COL_EVENLUM;
			}
			else if (g_str_has_prefix(arg, "-indep")) {
				stretch_colourmodel = COL_INDEP;
			}
			else if (g_str_has_prefix(arg,"-D=")) {
				arg += 3;
				D = g_ascii_strtod(arg, &end);
			}
			else if (g_str_has_prefix(arg,"-B=")) {
				if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
					arg += 3;
					B = g_ascii_strtod(arg, &end);
					if (fabs(B) < 1.e-3f)
						B = 0.f;
				} else {
					return CMD_ARG_ERROR;
				}
			}
			else if (g_str_has_prefix(arg,"-LP=")) {
				arg += 4;
				LP = g_ascii_strtod(arg, &end);
			}
			else if (g_str_has_prefix(arg,"-SP=")) {
				arg += 4;
				SP = g_ascii_strtod(arg, &end);
			}
			else if (g_str_has_prefix(arg,"-HP=")) {
				arg += 4;
				HP = g_ascii_strtod(arg, &end);
			}
		}
	}
	if (stretchtype == STRETCH_LINEAR) {
		if (BP != BP) {
			siril_log_message(_("Error: BP must be specified between 0.0 and 1.0 using -BP=\n"));
			return CMD_WRONG_N_ARG;
		}
		if (BP < 0. || BP > 1.) {
			siril_log_message(_("Error: BP must be >= 0.0 and <= 1.0.\n"));
			return CMD_ARG_ERROR;
		}
		D = 0.0;
		B = 0.0;
		LP = 0.0;
		SP = 0.0;
		HP = 0.0;
	} else {
		if (D != D) {
			siril_log_message(_("Error: D must be specified between 0.0 and 1.0 using -D=\n"));
			return CMD_WRONG_N_ARG;
		}

		if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
			if (B < 0. || B > 15.) {
				siril_log_message(_("Error: B must be >= 0.0 and <= 15.0.\n"));
				return CMD_ARG_ERROR;
			}
		} else {
			B = 0.0;
		}
		if (D <= 0. || D > 10.) {
			siril_log_message(_("Error: D must be > 0.0 and <= 10.0.\n"));
			return CMD_ARG_ERROR;
		}
		if (SP < 0. || SP > 1.) {
			siril_log_message(_("Error: SP must be >= 0.0 and <= 1.0.\n"));
			return CMD_ARG_ERROR;
		}
		if (LP < 0. || LP > SP) {
			siril_log_message(_("Error: LP must be >= 0.0 and <= SP.\n"));
			return CMD_ARG_ERROR;
		}
		if (HP < SP || HP > 1.) {
			siril_log_message(_("Error: HP must be >= SP and <= 1.0.\n"));
			return CMD_ARG_ERROR;
		}
		BP = 0.0;
		if (stretch_colourmodel == COL_SAT && (!(do_red && do_green && do_blue))) {
			siril_log_message(_("Error: saturation stretch requires that all channels must be selected.\n"));
			return CMD_ARG_ERROR;
		}
	}

	set_cursor_waiting(TRUE);
	params->B = (float) B;
	params->D = (float) expm1(D);
	params->LP = (float) LP;
	params->SP = (float) SP;
	params->HP = (float) HP;
	params->BP = (float) BP;
	params->stretchtype = stretchtype;
	params->payne_colourstretchmodel = stretch_colourmodel;
	params->do_red = do_red;
	params->do_green = do_green;
	params->do_blue = do_blue;
	params->clip_mode = clip_mode;
	return CMD_OK;
}

int process_seq_ghs(int nb, int stretchtype) {

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		siril_log_message(_("Error: cannot open sequence\n"));
		return CMD_SEQUENCE_NOT_FOUND;
	}

	ght_params *params = calloc(1, sizeof(ght_params));
	if (!params) {
		if (seq != &com.seq)
			free_sequence(seq, TRUE);
		return CMD_ALLOC_ERROR;
	}

	struct ght_data *seqdata = calloc(1, sizeof(struct ght_data));
	if (!seqdata) {
		free(params);
		if (seq != &com.seq)
			free_sequence(seq, TRUE);
		return CMD_ALLOC_ERROR;
	} else {
		seqdata->seq = seq;
		seqdata->params_ght = params;
	}

	gboolean mask_aware; // we don't use it for the sequence processing but need to provide it to process_ght_args
	int retval = process_ght_args(nb, TRUE, &mask_aware, stretchtype, params, seqdata);
	if (retval) {
		free(params);
		free(seqdata->seqEntry);
		free(seqdata);
		if (seq != &com.seq)
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}
	if (params->payne_colourstretchmodel == COL_SAT && seq->nb_layers != 3) {
		siril_log_message(_("Error: cannot apply saturation stretch to mono images.\n"));
		free(params);
		free(seqdata->seqEntry);
		free(seqdata);
		return CMD_ARG_ERROR;
	}
	if (!seqdata->seqEntry)
		seqdata->seqEntry = strdup("stretch_");
	sequence_cfa_warning_check(seq);
	apply_ght_to_sequence(seqdata);
	return CMD_OK;
}

int process_seq_ght(int nb) {
	return process_seq_ghs(nb, STRETCH_PAYNE_NORMAL);
}

int process_seq_invght(int nb) {
	return process_seq_ghs(nb, STRETCH_PAYNE_INVERSE);
}

int process_seq_modasinh(int nb) {
	return process_seq_ghs(nb, STRETCH_ASINH);
}

int process_seq_invmodasinh(int nb) {
	return process_seq_ghs(nb, STRETCH_INVASINH);
}

int process_seq_linstretch(int nb) {
	return process_seq_ghs(nb, STRETCH_LINEAR);
}

int process_mtf(int nb) {
	struct mtf_params params;
	gboolean inverse = word[0][0] == 'i' || word[0][0] == 'I';
	gboolean mask_aware = FALSE;
	gchar *end1, *end2, *end3;
	params.shadows = g_ascii_strtod(word[1], &end1);
	params.midtones = g_ascii_strtod(word[2], &end2);
	params.highlights = g_ascii_strtod(word[3], &end3);
	params.do_red = TRUE;
	params.do_green = TRUE;
	params.do_blue = TRUE;

	if (end1 == word[1] || end2 == word[2] || end3 == word[3] ||
			params.shadows < 0.0 || params.midtones <= 0.0 || params.highlights <= 0.0 ||
			params.shadows >= 1.0 || params.midtones >= 1.0 || params.highlights > 1.0) {
		siril_log_message(_("Invalid argument to %s, aborting.\n"), word[0]);
		return CMD_ARG_ERROR;
	}

	// Parse optional arguments (channel and/or -mask flag)
	for (int i = 4; word[i]; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!strcmp(word[i], "R")) {
			params.do_green = FALSE;
			params.do_blue = FALSE;
		} else if (!strcmp(word[i], "G")) {
			params.do_red = FALSE;
			params.do_blue = FALSE;
		} else if (!strcmp(word[i], "B")) {
			params.do_green = FALSE;
			params.do_red = FALSE;
		} else if (!strcmp(word[i], "RG")) {
			params.do_blue = FALSE;
		} else if (!strcmp(word[i], "RB")) {
			params.do_green = FALSE;
		} else if (!strcmp(word[i], "GB")) {
			params.do_red = FALSE;
		}
	}

	image_cfa_warning_check();

	// Create data structure
	struct mtf_data *data = create_mtf_data();
	if (!data) {
		return CMD_ALLOC_ERROR;
	}

	data->fit = gfit;
	data->params = params;
	data->auto_display_compensation = FALSE;
	data->is_preview = FALSE;

	// Create generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		destroy_mtf_data(data);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = inverse ? invmtf_single_image_hook : mtf_single_image_hook;
	args->log_hook = inverse ? invmtf_log_hook : mtf_log_hook;
	args->idle_function = NULL;  // No idle in command mode
	args->description = inverse ? _("Inverse Midtones Transfer Function") : _("Midtones Transfer Function");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Run worker synchronously - cleanup happens via destructor
	gpointer result = generic_image_worker(args);
	int retval = GPOINTER_TO_INT(result);

	if (retval != 0) {
		return CMD_GENERIC_ERROR;
	}

	siril_log_message(_("Applying MTF with values %f, %f, %f\n"),
		params.shadows, params.midtones, params.highlights);

	// Add to history
	char log[90];
	sprintf(log, "%s transfer (%.3f, %.4f, %.3f)",
			inverse ? "Inverse midtones" : "Midtones",
			params.shadows, params.midtones, params.highlights);
	gfit->history = g_slist_append(gfit->history, strdup(log));

	return CMD_OK;
}

int process_ghs(int nb, int stretchtype) {
	ght_params *params = calloc(1, sizeof(ght_params));
	if (!params)
		return CMD_ALLOC_ERROR;

	struct ght_data *seqdata = NULL;
	gboolean mask_aware;
	int retval = process_ght_args(nb, FALSE, &mask_aware, stretchtype, params, seqdata);
	if (retval) {
		free(params);
		return CMD_ARG_ERROR;
	}

	if (params->payne_colourstretchmodel == COL_SAT && gfit->naxes[2] != 3) {
		siril_log_message(_("Error: cannot apply saturation stretch to a mono image.\n"));
		free(params);
		return CMD_ARG_ERROR;
	}

	image_cfa_warning_check();

	// Create data structure
	struct ght_data *data = create_ght_data();
	if (!data) {
		free(params);
		return CMD_ALLOC_ERROR;
	}

	data->fit = gfit;
	data->params_ght = params;  // Take ownership of params
	data->auto_display_compensation = FALSE;
	data->is_preview = FALSE;

	// Create generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		destroy_ght_data(data);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = (params->payne_colourstretchmodel == COL_SAT) ? 2.0f : 1.0f;
	args->image_hook = ght_single_image_hook;
	args->log_hook = ght_log_hook;
	args->idle_function = NULL;  // No idle in command mode
	args->description = _("Generalised Hyperbolic Stretch");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Run worker synchronously - cleanup happens via destructor
	gpointer result = generic_image_worker(args);
	retval = GPOINTER_TO_INT(result);

	if (retval != 0) {
		return CMD_GENERIC_ERROR;
	}

	// Add to history
	char log[100];
	switch (stretchtype) {
		case STRETCH_PAYNE_NORMAL:
			sprintf(log, "GHS (pivot: %.3f, amount: %.2f, local: %.1f [%.2f, %.2f])",
				params->SP, params->D, params->B, params->LP, params->HP);
			break;
		case STRETCH_PAYNE_INVERSE:
			sprintf(log, "Inverse GHS (pivot: %.3f, amount: %.2f, local: %.1f [%.2f, %.2f])",
				params->SP, params->D, params->B, params->LP, params->HP);
			break;
		case STRETCH_ASINH:
			sprintf(log, "GHS asinh (pivot: %.3f, amount: %.2f [%.2f, %.2f])",
				params->SP, params->D, params->LP, params->HP);
			break;
		case STRETCH_INVASINH:
			sprintf(log, "GHS inverse asinh (pivot: %.3f, amount: %.2f [%.2f, %.2f])",
				params->SP, params->D, params->LP, params->HP);
			break;
		case STRETCH_LINEAR:
			sprintf(log, "GHS BP shift (new BP: %.3f)", params->BP);
			break;
	}
	gfit->history = g_slist_append(gfit->history, strdup(log));

	if (gui.roi.active)
		populate_roi();

	return CMD_OK;
}

int process_autoghs(int nb) {
	int argidx = 1;
	gboolean linked = FALSE, mask_aware = FALSE;
	clip_mode_t clip_mode = RGBBLEND;
	float shadows_clipping, b = 13.0f, hp = 0.7f, lp = 0.0f;

	if (!g_strcmp0(word[1], "-linked")) {
		linked = TRUE;
		argidx++;
	}

	gchar *end = NULL;
	shadows_clipping = g_ascii_strtod(word[argidx], &end);
	if (end == word[argidx]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
		return CMD_ARG_ERROR;
	}
	argidx++;

	float amount = g_ascii_strtod(word[argidx], &end);
	if (end == word[argidx]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
		return CMD_ARG_ERROR;
	}
	argidx++;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-b=")) {
			char *arg = word[argidx] + 3;
			b = g_ascii_strtod(arg, &end);
			if (fabsf(b) < 1.e-3f)
				b = 0.f;
			if (arg == end || b < -5.0f || b > 15.0f) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-hp=")) {
			char *arg = word[argidx] + 4;
			hp = g_ascii_strtod(arg, &end);
			if (arg == end || hp < 0.0f || hp > 1.0f) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-lp=")) {
			char *arg = word[argidx] + 4;
			lp = g_ascii_strtod(arg, &end);
			if (arg == end || lp < 0.0f || lp > 1.0f) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-clipmode=")) {
			char *argument = word[argidx] + 10;
			if (!g_ascii_strncasecmp(argument, "clip", 4))
				clip_mode = CLIP;
			else if (!g_ascii_strncasecmp(argument, "rescale", 7))
				clip_mode = RESCALE;
			else if (!g_ascii_strncasecmp(argument, "globalrescale", 13))
				clip_mode = RESCALEGLOBAL;
			else if (!g_ascii_strncasecmp(argument, "rgbblend", 8))
				clip_mode = RGBBLEND;
			else {
				siril_log_color_message(_("Error, unknown clip mode %s specified\n"), "red", argument);
				return CMD_ARG_ERROR;
			}
		}
		else if (!g_strcmp0(word[argidx], "-mask")) {
			mask_aware = TRUE;
		}
		argidx++;
	}

	int nb_channels = (int)gfit->naxes[2];
	imstats *stats[3] = { NULL };
	int ret = compute_all_channels_statistics_single_image(gfit, STATS_BASIC, MULTI_THREADED, stats);
	if (ret) {
		for (int i = 0; i < nb_channels; ++i) {
			free_stats(stats[i]);
		}
		return CMD_GENERIC_ERROR;
	}

	image_cfa_warning_check();

	if (linked) {
		double median = 0.0, sigma = 0.0;
		for (int i = 0; i < nb_channels; ++i) {
			median += stats[i]->median;
			sigma += stats[i]->sigma;
			free_stats(stats[i]);
		}
		median /= nb_channels;
		sigma /= nb_channels;
		float SP = median + shadows_clipping * sigma;
		if (gfit->type == DATA_USHORT)
			SP *= (gfit->orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE;
		siril_log_message(_("Symmetry point SP=%f\n"), SP);

		// Create data structure for linked processing
		struct ght_data *data = create_ght_data();
		if (!data) {
			return CMD_ALLOC_ERROR;
		}

		data->fit = gfit;
		data->params_ght = malloc(sizeof(ght_params));
		if (!data->params_ght) {
			destroy_ght_data(data);
			return CMD_ALLOC_ERROR;
		}

		data->params_ght->B = b;
		data->params_ght->D = amount;
		data->params_ght->LP = lp;
		data->params_ght->SP = SP;
		data->params_ght->HP = hp;
		data->params_ght->BP = 0.0;
		data->params_ght->stretchtype = STRETCH_PAYNE_NORMAL;
		data->params_ght->payne_colourstretchmodel = COL_INDEP;
		data->params_ght->do_red = TRUE;
		data->params_ght->do_green = TRUE;
		data->params_ght->do_blue = TRUE;
		data->params_ght->clip_mode = clip_mode;
		data->auto_display_compensation = FALSE;
		data->is_preview = FALSE;

		// Create generic_img_args
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_ght_data(data);
			return CMD_ALLOC_ERROR;
		}

		args->fit = gfit;
		args->mem_ratio = 1.0f;
		args->image_hook = ght_single_image_hook;
		args->log_hook = ght_log_hook;
		args->idle_function = NULL;  // No idle in command mode
		args->description = _("AutoGHS");
		args->command_updates_gfit = TRUE;
		args->command = TRUE; // calling as command, not from GUI
		args->verbose = TRUE;
		args->user = data;
		args->mask_aware = mask_aware;
		args->max_threads = com.max_thread;
		args->for_preview = FALSE;
		args->for_roi = FALSE;

		// Run worker synchronously - cleanup happens via destructor
		if (!start_in_new_thread(generic_image_worker, args)) {
			free_generic_img_args(args);
			return CMD_GENERIC_ERROR;
		}

	} else {
		// Unlinked mode - process each channel independently
		// For now, keep direct processing as it uses apply_ght_to_fits_channel
		// TODO: Implement per-channel worker support for mask compatibility
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(nb_channels > 1)
#endif
		for (int i = 0; i < nb_channels; ++i) {
			if (stats[i]) {
				float SP = stats[i]->median + shadows_clipping * stats[i]->sigma;
				if (gfit->type == DATA_USHORT)
					SP *= (gfit->orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE;
				siril_log_message(_("Symmetry point for channel %d: SP=%f\n"), i, SP);

				ght_params params = { .B = b, .D = amount, .LP = lp, .SP = SP, .HP = hp,
					.BP = 0.0, STRETCH_PAYNE_NORMAL, COL_INDEP, TRUE, TRUE, TRUE};
				apply_ght_to_fits_channel(gfit, gfit, i, &params, TRUE);

				free_stats(stats[i]);
			}
		}
	}

	char log[100];
	sprintf(log, "AutoGHS (%sk.sigma: %.2f, amount: %.2f, local: %.1f [%.2f, %.2f])",
			linked ? "linked, " : "", shadows_clipping, amount, b, lp, hp);
	gfit->history = g_slist_append(gfit->history, strdup(log));

	return CMD_OK;
}

int process_autostretch(int nb) {
	int arg_index = 1;
	gboolean linked = FALSE;
	gboolean mask_aware = FALSE;

	/* Parse optional flags (order-independent) */
	while (arg_index < nb) {
		if (!strcmp(word[arg_index], "-linked")) {
			linked = TRUE;
			arg_index++;
		} else if (!strcmp(word[arg_index], "-mask")) {
			mask_aware = TRUE;
			arg_index++;
		} else {
			break; /* first non-flag argument */
		}
	}

	gchar *end = NULL;
	float shadows_clipping = AS_DEFAULT_SHADOWS_CLIPPING;
	if (nb > arg_index) {
		shadows_clipping = g_ascii_strtod(word[arg_index], &end);
		if (end == word[arg_index]) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[arg_index]);
			return CMD_ARG_ERROR;
		}
		arg_index++;
	}

	float target_bg = AS_DEFAULT_TARGET_BACKGROUND;
	if (nb > arg_index) {
		target_bg = g_ascii_strtod(word[arg_index], &end);
		if (end == word[arg_index] || target_bg < 0.0f || target_bg > 1.0f) {
			siril_log_message(_("The target background value must be in the [0, 1] range\n"));
			return CMD_ARG_ERROR;
		}
	}

	siril_log_message(_("Computing %s auto-stretch with values %f and %f\n"),
			linked ? _("linked") : _("unlinked"), shadows_clipping, target_bg);

	image_cfa_warning_check();

	/* Create data structure */
	struct mtf_data *data = create_mtf_data();
	if (!data) {
		return CMD_ALLOC_ERROR;
	}
	data->fit = gfit;
	data->auto_display_compensation = FALSE;
	data->is_preview = FALSE;
	data->linked = linked;

	/* Create generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		destroy_mtf_data(data);
		return CMD_ALLOC_ERROR;
	}

	/* Compute the autostretch parameters */
	if (linked) {
		find_linked_midtones_balance(gfit, shadows_clipping, target_bg, &data->params);
		data->params.do_red = data->params.do_green = data->params.do_blue = TRUE;
	} else {
		find_unlinked_midtones_balance(gfit, shadows_clipping, target_bg, data->uparams);
		data->uparams[0].do_red = data->uparams[0].do_green = data->uparams[0].do_blue = TRUE;
		data->uparams[1].do_red = data->uparams[1].do_green = data->uparams[1].do_blue = TRUE;
		data->uparams[2].do_red = data->uparams[2].do_green = data->uparams[2].do_blue = TRUE;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = mtf_single_image_hook;
	args->log_hook = mtf_log_hook;
	args->idle_function = NULL;  // No idle in command mode
	args->description = _("Autostretch");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	/* Run worker synchronously - cleanup happens via destructor */
	gpointer result = generic_image_worker(args);
	int retval = GPOINTER_TO_INT(result);
	if (retval != 0) {
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_ght(int nb) {
	return process_ghs(nb, STRETCH_PAYNE_NORMAL);
}

int process_invght(int nb) {
	return process_ghs(nb, STRETCH_PAYNE_INVERSE);
}

int process_modasinh(int nb) {
	return process_ghs(nb, STRETCH_ASINH);
}

int process_invmodasinh(int nb) {
	return process_ghs(nb, STRETCH_INVASINH);
}

int process_linstretch(int nb) {
	return process_ghs(nb, STRETCH_LINEAR);
}

int process_wavelet(int nb) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];

	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;
	const char* tmpdir = g_get_tmp_dir();
	gchar *end1, *end2;

	Nbr_Plan = g_ascii_strtoull(word[1], &end1, 10);
	Type_Transform = g_ascii_strtoull(word[2], &end2, 10);

	nb_chan = gfit->naxes[2];
	g_assert(nb_chan <= 3);

	mins = min (gfit->rx, gfit->ry);
	maxplan = log(mins) / log(2) - 2;

	if (end1 == word[1] || Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return CMD_ARG_ERROR;
	}

	if(end2 == word[2] || (Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE)){
		siril_log_message(_("Wavelet: type must be %d or %d\n"), TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		return CMD_ARG_ERROR;
	}
	image_cfa_warning_check();
	if (gfit->type == DATA_USHORT) {
		float *Imag = f_vector_alloc(gfit->rx * gfit->ry);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return CMD_ALLOC_ERROR;
		}

		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan], NULL);
			wavelet_transform_file(Imag, gfit->ry, gfit->rx, dir[chan],
					Type_Transform, Nbr_Plan, gfit->pdata[chan]);
			g_free(dir[chan]);
		}

		free(Imag);
	} else if (gfit->type == DATA_FLOAT) {
		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan], NULL);
			wavelet_transform_file_float(gfit->fpdata[chan], gfit->ry, gfit->rx, dir[chan],
					Type_Transform, Nbr_Plan);
			g_free(dir[chan]);
		}
	}
	else return CMD_INVALID_IMAGE;
	return CMD_OK;
}

int process_log(int nb){
	loglut(gfit);
	notify_gfit_modified();
	return CMD_OK | CMD_NOTIFY_GFIT_MODIFIED;
}

int process_linear_match(int nb) {
	fits ref = { 0 };
	gchar *end1, *end2;
	double a[3] = { 0.0 }, b[3] = { 0.0 };
	double low = g_ascii_strtod(word[2], &end1);
	double high = g_ascii_strtod(word[3], &end2);

	if (end1 == word[2] || low < 0 || low > 1) {
		siril_log_message(_("Low value must be in the [0, 1] range.\n"));
		return CMD_ARG_ERROR;
	}

	if (end1 == word[3] || high < 0 || high > 1) {
		siril_log_message(_("High value must be in the [0, 1] range.\n"));
		return CMD_ARG_ERROR;
	}

	if (readfits(word[1], &ref, NULL, gfit->type == DATA_FLOAT))
		return CMD_INVALID_IMAGE;
	cmd_errors retval = CMD_OK;
	if (!find_linear_coeff(gfit, &ref, low, high, a, b, NULL)) {
		image_cfa_warning_check();
		set_cursor_waiting(TRUE);
		apply_linear_to_fits(gfit, a, b);
		notify_gfit_modified();
		retval |= CMD_NOTIFY_GFIT_MODIFIED;
	}
	clearfits(&ref);
	return retval;
}

int process_asinh(int nb) {
	gboolean human_luminance = FALSE;
	clip_mode_t clip_mode = RGBBLEND;
	gchar *end;
	int arg_offset = 1;
	if (!strcmp(word[1], "-human")) {
		human_luminance = TRUE;
		arg_offset++;
	}
	if (nb <= arg_offset)
		return CMD_WRONG_N_ARG;
	double beta = g_ascii_strtod(word[arg_offset], &end);
	if (end == word[arg_offset] || beta < 1.0) {
		siril_log_color_message(_("Stretch must be greater than or equal to 1\n"), "red");
		return CMD_ARG_ERROR;
	}
	arg_offset++;

	double offset = 0.0;
	gboolean use_mask = FALSE;
	while (arg_offset < nb) {
		if (g_str_has_prefix(word[arg_offset], "-clipmode=")) {
			char *argument = word[arg_offset] + 10;
			if (!g_ascii_strncasecmp(argument, "clip", 4))
				clip_mode = CLIP;
			else if (!g_ascii_strncasecmp(argument, "rescale", 7))
				clip_mode = RESCALE;
			else if (!g_ascii_strncasecmp(argument, "globalrescale", 13))
				clip_mode = RESCALEGLOBAL;
			else if (!g_ascii_strncasecmp(argument, "rgbblend", 8))
				clip_mode = RGBBLEND;
			else {
				siril_log_color_message(_("Invalid clip mode %s, aborting.\n"), "red", argument);
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(word[arg_offset], "-mask")) {
			use_mask = TRUE;
		} else {
			offset = g_ascii_strtod(word[arg_offset], &end);
			if (end == word[arg_offset]) {
				siril_log_color_message(_("Invalid argument %s, aborting.\n"), "red", word[arg_offset]);
				return CMD_ARG_ERROR;
			}
		}
		arg_offset++;
	}

	set_cursor_waiting(TRUE);
	image_cfa_warning_check();

	// Allocate parameters
	asinh_params *params = calloc(1, sizeof(asinh_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->beta = beta;
	params->offset = offset;
	params->human_luminance = human_luminance;
	params->clip_mode = clip_mode;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = asinh_image_hook;
	args->idle_function = NULL; // No idle function for synchronous execution
	args->description = _("Asinh stretch");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = FALSE;
	args->user = params;
	args->log_hook = asinh_log_hook;
	args->max_threads = com.max_thread;
	args->mask_aware = use_mask;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Run synchronously by calling the worker directly
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	// Check return value
	int retval = args->retval;

	// Resources are freed at the end of the generic image worker

	if (retval != 0) {
		siril_log_color_message(_("Asinh stretch failed\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	char log[90];
	sprintf(log, "Asinh stretch (amount: %.1f, offset: %.1f, human: %s)", beta, offset, human_luminance ? "yes" : "no");
	gfit->history = g_slist_append(gfit->history, strdup(log));
	siril_log_message(log);

	return CMD_OK;
}

int process_clahe(int nb) {
	gchar *end;
	double clip_limit = g_ascii_strtod(word[1], &end);
	if (end == word[1] || clip_limit <= 0.0) {
		siril_log_message(_("Clip limit must be > 0.\n"));
		return CMD_ARG_ERROR;
	}
	int size = g_ascii_strtoull(word[2], &end, 10);
	if (end == word[2] || size <= 0) {
		siril_log_message(_("Tile size must be > 0.\n"));
		return CMD_ARG_ERROR;
	}

	image_cfa_warning_check();

	// Allocate parameters
	clahe_params *params = calloc(1, sizeof(clahe_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_GENERIC_ERROR;
	}

	params->clip = clip_limit;
	params->tileSize = size;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = clahe_image_hook;
	args->idle_function = NULL; // Use default idle function for command-line
	args->description = _("CLAHE");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = TRUE;
	args->user = params;
	args->mask_aware = word[3] && g_strcmp0(word[3], "-mask") == 0; // handle the -mask flag to set mask_aware state
	args->log_hook = clahe_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free(params);
		free(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_ls(int nb){
	gchar *path = NULL;

	/* If a path is given in argument */
	if (nb > 1) {
		if (word[1][0] != '\0') {
			/* Absolute path */
			if (word[1][0] == G_DIR_SEPARATOR || word[1][0] == '~') {
				long maxpath = get_pathmax();
				char filename[maxpath];

				g_strlcpy(filename, word[1], 250);
				expand_home_in_filename(filename, maxpath);
				path = g_build_filename(filename, NULL);
			}
			/* Relative path */
			else {
				path = g_build_filename(com.wd, word[1], NULL);
			}
		}
		/* Should not happen */
		else {
			printf("Cannot list files in %s\n", word[1]);
			return CMD_GENERIC_ERROR;
		}
	}
	/* No paths are given in argument */
	else {
		if (!com.wd) {
			siril_log_message(_("Cannot list files, set working directory first.\n"));
			return CMD_GENERIC_ERROR;
		}
		path = g_strdup(com.wd);
	}
	if (path == NULL) {
		siril_log_message(_("Siril cannot open the directory.\n"));
		return CMD_NO_CWD;
	}

#ifndef _WIN32
	struct dirent **list = NULL;

	int n = scandir(path, &list, 0, alphasort);
	if (n < 0) {
		perror("scandir");
		siril_log_message(_("Siril cannot open the directory.\n"));
		g_free(path);
		return CMD_NO_CWD;
	}

	/* List the entries */
	for (int i = 0; i < n; ++i) {
		GStatBuf entrystat;
		gchar *filename;
		const char *ext;
		if (list[i]->d_name[0] == '.')
			continue; /* no hidden files */

		filename = g_build_filename(path, list[i]->d_name, NULL);

		if (g_lstat(filename, &entrystat)) {
			perror("stat");
			g_free(filename);
			break;
		}
		g_free(filename);
		if (S_ISLNK(entrystat.st_mode)) {
			siril_log_color_message(_("Link: %s\n"), "bold", list[i]->d_name);
			continue;
		}
		if (S_ISDIR(entrystat.st_mode)) {
			siril_log_color_message(_("Directory: %s\n"), "green",
					list[i]->d_name);
			continue;
		}
		ext = get_filename_ext(list[i]->d_name);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type != TYPEUNDEF) {
			if (type == TYPEAVI || type == TYPESER)
				siril_log_color_message(_("Sequence: %s\n"), "salmon",
						list[i]->d_name);
			else if (type == TYPEFITS)
				siril_log_color_message(_("Image: %s\n"), "plum", list[i]->d_name);
			else
				siril_log_color_message(_("Image: %s\n"), "red", list[i]->d_name);
		} else if (!strncmp(ext, "seq", 3))
			siril_log_color_message(_("Sequence: %s\n"), "blue", list[i]->d_name);
	}
	for (int i = 0; i < n; i++)
		free(list[i]);
	free(list);
#else
	WIN32_FIND_DATAW fdFile;
	HANDLE hFind = NULL;
	char sPath[2048];

	//Specify a file mask. *.* = We want everything
	sprintf(sPath, "%s\\*.*", path);

	wchar_t *wpath = g_utf8_to_utf16(sPath, -1, NULL, NULL, NULL);
	if (wpath == NULL)
		return CMD_ALLOC_ERROR;

	if ((hFind = FindFirstFileW(wpath, &fdFile)) == INVALID_HANDLE_VALUE) {
		siril_log_message(_("Siril cannot open the directory.\n"));
		g_free(wpath);
		return CMD_NO_CWD;
	}

	g_free(wpath);
	do {
		//Find first file will always return "."
		//    and ".." as the first two directories.
		if (wcscmp(fdFile.cFileName, L".") != 0
				&& wcscmp(fdFile.cFileName, L"..") != 0) {

			gchar *filename = g_utf16_to_utf8(fdFile.cFileName, -1, NULL, NULL, NULL);
			//Is the entity a File or Folder?
			if (fdFile.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
				siril_log_color_message(_("Directory: %s\n"), "green", filename);
			} else {
				const char *ext = get_filename_ext(filename);
				if (!ext)
					continue;
				image_type type = get_type_for_extension(ext);
				if (type != TYPEUNDEF) {
					if (type == TYPEAVI || type == TYPESER)
						siril_log_color_message(_("Sequence: %s\n"), "salmon", filename);
					else if (type == TYPEFITS)
						siril_log_color_message(_("Image: %s\n"), "plum", filename);
					else
						siril_log_color_message(_("Image: %s\n"), "red", filename);
				} else if (!strncmp(ext, "seq", 3))
					siril_log_color_message(_("Sequence: %s\n"), "blue", filename);

			}
			g_free(filename);
		}
	} while(FindNextFileW(hFind, &fdFile)); //Find the next file.

	FindClose(hFind); // Clean things up!
#endif
	siril_log_message(_("********* END OF THE LIST *********\n"));
	g_free(path);

	return CMD_OK;
}

int process_merge(int nb) {
	int retval = 0, nb_seq = nb-2;
	if (!com.wd) {
		siril_log_message(_("Merge: no working directory set.\n"));
		set_cursor_waiting(FALSE);
		return CMD_NO_CWD;
	}
	if (file_name_has_invalid_chars(word[nb - 1])) {
		siril_log_color_message(_("Specified output name %s contains forbidden characters, aborting\n"), "red", word[nb - 1]);
		return CMD_ARG_ERROR;
	}
	char *dest_dir = strdup(com.wd);
	char *outseq_name = NULL;
	sequence **seqs = calloc(nb_seq, sizeof(sequence *));
	GList *list = NULL;
	for (int i = 0; i < nb_seq; i++) {
		char *seqpath1 = strdup(word[i + 1]), *seqpath2 = strdup(word[i + 1]);
		char *dir = g_path_get_dirname(seqpath1);
		char *seqname = g_path_get_basename(seqpath2);
#ifdef _WIN32
		gchar **token = g_strsplit(dir, "/", -1);
		g_free(dir);
		dir = g_strjoinv(G_DIR_SEPARATOR_S, token);
		g_strfreev(token);
#endif
		if (dir[0] != '\0' && !(dir[0] == '.' && dir[1] == '\0'))
			siril_change_dir(dir, NULL);
		if (!(seqs[i] = load_sequence(seqname, NULL))) {
			siril_log_message(_("Could not open sequence `%s' for merging\n"), word[i + 1]);
			retval = CMD_SEQUENCE_NOT_FOUND;
			free(seqpath1); free(seqpath2);	g_free(seqname); g_free(dir);
			goto merge_clean_up;
		}
		g_free(seqname);
		if (seq_check_basic_data(seqs[i], FALSE) < 0) {
			siril_log_message(_("Sequence `%s' is invalid, could not merge\n"), word[i + 1]);
			retval = CMD_GENERIC_ERROR;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (i != 0 && (seqs[i]->rx != seqs[0]->rx ||
					seqs[i]->ry != seqs[0]->ry ||
					seqs[i]->nb_layers != seqs[0]->nb_layers ||
					seqs[i]->bitpix != seqs[0]->bitpix ||
					seqs[i]->type != seqs[0]->type)) {
			siril_log_message(_("All sequences must be the same format for merging. Sequence `%s' is different\n"), word[i + 1]);
			retval = CMD_GENERIC_ERROR;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (seqs[i]->type == SEQ_REGULAR) {
			// we need to build the list of files
			long maxpath = get_pathmax();
			char filename[maxpath];
			for (int image = 0; image < seqs[i]->number; image++) {
				fit_sequence_get_image_filename_checkext(seqs[i], image, filename);
				list = g_list_append(list, g_build_filename(dir, filename, NULL));
			}
		}
		free(seqpath1); free(seqpath2); g_free(dir);
		siril_change_dir(dest_dir, NULL);	// they're all relative to this one
	}

	struct ser_struct out_ser;
	struct _convert_data *args = NULL;
	fitseq out_fitseq;
	char *destroot;
	switch (seqs[0]->type) {
		case SEQ_REGULAR:
			// use the conversion, it makes symbolic links or copies as a fallback
			destroot = strdup(word[nb - 1]);
			args = calloc(1, sizeof(struct _convert_data));
			args->start = 0;
			args->total = 0; // init to get it from glist_to_array()
			args->list = glist_to_array(list, &args->total);
			args->destroot = destroot;
			args->input_has_a_seq = FALSE;
			args->input_has_a_film = FALSE;
			args->debayer = FALSE;
			args->multiple_output = FALSE;
			args->output_type = SEQ_REGULAR;
			args->make_link = TRUE;
			gettimeofday(&(args->t_start), NULL);
			if (!start_in_new_thread(convert_thread_worker, args)) {
				free(args->destroot);
				g_strfreev(args->list);
				free(args);
			}
			break;

		case SEQ_SER:
			if (g_str_has_suffix(word[nb - 1], ".ser"))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s.ser", word[nb - 1]);
			if (ser_create_file(outseq_name, &out_ser, TRUE, seqs[0]->ser_file)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = CMD_GENERIC_ERROR;
				goto merge_clean_up;
			}
			seqwriter_set_max_active_blocks(2);
			int written_frames = 0;
			for (int i = 0; i < nb_seq; i++) {
				for (unsigned int frame = 0; frame < seqs[i]->number; frame++) {
					seqwriter_wait_for_memory();
					fits *fit = calloc(1, sizeof(fits));
					if (ser_read_frame(seqs[i]->ser_file, frame, fit, FALSE, com.pref.debayer.open_debayer)) {
						siril_log_message(_("Failed to read frame %d from input sequence `%s'\n"), frame, word[i + 1]);
						retval = CMD_INVALID_IMAGE;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}

					if (ser_write_frame_from_fit(&out_ser, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = CMD_GENERIC_ERROR;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (ser_write_and_close(&out_ser)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = CMD_GENERIC_ERROR;
			}
			break;

		case SEQ_FITSEQ:
			if (g_str_has_suffix(word[nb - 1], com.pref.ext))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s%s", word[nb - 1], com.pref.ext);
			if (fitseq_create_file(outseq_name, &out_fitseq, -1)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = CMD_GENERIC_ERROR;
				goto merge_clean_up;
			}
			g_free(outseq_name);
			outseq_name = NULL;
			seqwriter_set_max_active_blocks(2);
			written_frames = 0;
			for (int i = 0; i < nb_seq; i++) {
				for (unsigned int frame = 0; frame < seqs[i]->number; frame++) {
					seqwriter_wait_for_memory();
					fits *fit = calloc(1, sizeof(fits));
					if (fitseq_read_frame(seqs[i]->fitseq_file, frame, fit, FALSE, -1)) {
						siril_log_message(_("Failed to read frame %d from input sequence `%s'\n"), frame, word[i + 1]);
						retval = CMD_INVALID_IMAGE;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}

					if (fitseq_write_image(&out_fitseq, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = CMD_GENERIC_ERROR;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (fitseq_close_file(&out_fitseq)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = CMD_GENERIC_ERROR;
			}
			break;
		default:
			siril_log_message(_("This type of sequence cannot be created by Siril, aborting the merge\n"));
			retval = CMD_GENERIC_ERROR;
	}

merge_clean_up:
	for (int i = 0; i < nb_seq; i++) {
		if (seqs[i])
			free_sequence(seqs[i], TRUE);
	}
	free(seqs);
	g_free(outseq_name);
	siril_change_dir(dest_dir, NULL);
	free(dest_dir);
	return retval;
}

/* Geometry operation image commands (updated to use generic_image_worker) */

int process_mirrorx_single(int nb){
	image_type imagetype;
	char *realname = NULL;
	if (stat_file(word[1], &imagetype, &realname)) {
		siril_log_color_message(_("Error opening image %s: file not found or not supported.\n"), "red", word[1]);
		free(realname);
		return CMD_FILE_NOT_FOUND;
	}
	if (imagetype != TYPEFITS && imagetype != TYPETIFF) {
		siril_log_color_message(_("This command is only supported with FITS and TIFF, able to contain orientation information\n"), "red");
		free(realname);
		return CMD_INVALID_IMAGE;
	}
	if (imagetype == TYPEFITS && fitseq_is_fitseq(realname, NULL)) {
		siril_log_color_message(_("This command is only supported with single FITS images, for the first HDU, not a FITS cube.\n"), "red");
		free(realname);
		return CMD_INVALID_IMAGE;
	}

	fits fit = { 0 };
	if (read_fits_metadata_from_path(realname, &fit)) {
		siril_log_color_message(_("Could not open file: %s\n"), "red", realname);
		clearfits(&fit);
		free(realname);
		return CMD_ARG_ERROR;
	}
	if (!strcmp(fit.keywords.row_order, "BOTTOM-UP")) {
		siril_log_message(_("Image data is already bottom-up\n"));
		clearfits(&fit);
		free(realname);
		return CMD_OK;
	}
	clearfits(&fit);
	siril_log_message(_("Mirroring image to convert to bottom-up data\n"));
	if (readfits(realname, &fit, NULL, FALSE)) {
		siril_log_color_message(_("Could not open file: %s\n"), "red", realname);
		clearfits(&fit);
		free(realname);
		return CMD_ARG_ERROR;
	}

	mirrorx(&fit, TRUE);

	int retval = CMD_OK;
	if (savefits(realname, &fit)) {
		siril_log_color_message(_("Could not save mirrored image: %s\n"), "red", realname);
		retval = CMD_ARG_ERROR;
	}
	clearfits(&fit);
	free(realname);
	return retval;
}

int process_mirrorx(int nb){
	if (nb == 2 && !strcmp(word[1], "-bottomup")) {
		if (!strcmp(gfit->keywords.row_order, "BOTTOM-UP")) {
			siril_log_message(_("Image data is already bottom-up\n"));
			return CMD_OK;
		}
		siril_log_message(_("Mirroring image to convert to bottom-up data\n"));
	}

	// Allocate parameters
	struct mirror_args *params = new_mirror_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->x_axis = TRUE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_mirror_args(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = mirrorx_image_hook;
	args->idle_function = NULL;  // Use default
	args->description = _("Mirror X");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;  // Don't add CMD_NOTIFY_GFIT_MODIFIED - handled by worker
}

int process_mirrory(int nb){
	// Allocate parameters
	struct mirror_args *params = new_mirror_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->x_axis = FALSE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_mirror_args(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = mirrory_image_hook;
	args->idle_function = NULL;
	args->description = _("Mirror Y");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_binxy(int nb) {
	gchar *end;
	int factor = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || factor <= 0) {
		siril_log_message(_("Factor must be a number greater than 0.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	gboolean mean = TRUE;

	if (nb > 2 && !g_ascii_strncasecmp(word[2], "-sum", 4)) {
		mean = FALSE;
	}
	image_cfa_warning_check();

	// Allocate parameters
	struct binning_args *params = new_binning_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->factor = factor;
	params->mean = mean;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_binning_args(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = binning_image_hook;
	args->log_hook = binning_log_hook;
	args->idle_function = NULL;
	args->description = _("Binning");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_resample(int nb) {
	gchar *end;
	gboolean clamp = TRUE;
	opencv_interpolation interpolation = OPENCV_LANCZOS4;
	int toX, toY;

	if (word[1][0] == '-') {
		if (g_str_has_prefix(word[1], "-height=")) {
			toY = g_ascii_strtoull(word[1] + 8, &end, 10);
			if (end == word[1] + 9) {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				return CMD_ARG_ERROR;
			}
			toX = round_to_int(gfit->rx * (double)toY / gfit->ry);
		} else if (g_str_has_prefix(word[1], "-width=")) {
			toX = g_ascii_strtoull(word[1] + 7, &end, 10);
			if (end == word[1]+8) {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				return CMD_ARG_ERROR;
			}
			toY = round_to_int(gfit->ry * (double)toX / gfit->rx);
		} else if (g_str_has_prefix(word[1], "-maxdim=")) {
			double maxdim = g_ascii_strtoull(word[1] + 8, &end, 10);
			if (end == word[1]+9) {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				return CMD_ARG_ERROR;
			}
			if (gfit->rx > gfit->ry) {
				toX = maxdim;
				toY = round_to_int(gfit->ry * (double)toX / gfit->rx);
			} else {
				toY = maxdim;
				toX = round_to_int(gfit->rx * (double)toY / gfit->ry);
			}
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[1]);
			return CMD_ARG_ERROR;
		}
	} else {
		double factor = g_ascii_strtod(word[1], &end);
		if (end == word[1] || factor <= 0.2 || factor > 5.0) {
			siril_log_message(_("Scale %lf not allowed. Should be between 0.2 and 5.0.\n"), factor);
			return CMD_ARG_ERROR;
		}
		if (factor == 1.0) {
			siril_log_message(_("Scale is 1.0. Not doing anything.\n"));
			return CMD_ARG_ERROR;
		}
		toX = round_to_int(factor * gfit->rx);
		toY = round_to_int(factor * gfit->ry);
	}

	for (int i = 2; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				interpolation = OPENCV_AREA;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			return CMD_ARG_ERROR;
		} else if (!strcmp(word[i], "-noclamp")) {
			clamp = FALSE;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}
	image_cfa_warning_check();
	set_cursor_waiting(TRUE);

	// Allocate parameters
	struct resample_args *params = new_resample_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	params->toX = toX;
	params->toY = toY;
	params->interpolation = interpolation;
	params->clamp = clamp;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_resample_args(params);
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f + ((toX / gfit->rx) * (toY / gfit->ry));
	args->image_hook = resample_image_hook;
	args->log_hook = resample_log_hook;
	args->idle_function = NULL;
	args->description = _("Resample");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_crop(int nb) {
	if (is_preview_active()) {
		siril_log_message(_("It is impossible to crop the image when a filter with preview session is active. "
						"Please consider to close the filter dialog first.\n"));
		return CMD_GENERIC_ERROR;
	}

	rectangle area;
	if (com.selection.h && com.selection.w) {
		area = com.selection;
	} else {
		if (nb == 5) {
			gchar *end1, *end2;
			area.x = g_ascii_strtoull(word[1], &end1, 10);
			area.y = g_ascii_strtoull(word[2], &end2, 10);
			if (end1 == word[1] || area.x < 0 || end2 == word[2] || area.y < 0) {
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return CMD_ARG_ERROR;
			}
			area.w = g_ascii_strtoull(word[3], &end1, 10);
			area.h = g_ascii_strtoull(word[4], &end2, 10);
			if (end1 == word[3] || area.w < 0 || end2 == word[4] || area.h < 0) {
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return CMD_ARG_ERROR;
			}
			if (area.x + area.w > gfit->rx || area.y + area.h > gfit->ry) {
				siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"), gfit->rx, gfit->ry);
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
			return CMD_ARG_ERROR;
		}
	}

	// Allocate parameters
	struct crop_args *params = new_crop_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->area = area;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_crop_args(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = crop_image_hook_single;
	args->idle_function = NULL;
	args->description = _("Crop");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = 1;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_rotate(int nb) {
	gchar *end;
	set_cursor_waiting(TRUE);
	int crop = 1;
	gboolean has_selection = FALSE;
	rectangle area = { 0, 0, gfit->rx, gfit->ry };
	if (com.selection.w > 0 && com.selection.h > 0) {
		siril_log_color_message(_("Rotation will apply only to current selection, the resulting image will be cropped.\n"), "salmon");
		area = com.selection;
		has_selection = TRUE;
	}

	gboolean clamp = TRUE;
	int interpolation = OPENCV_LANCZOS4;
	double angle=0.0;

	angle = g_ascii_strtod(word[1], &end);
	if (end == word[1] || angle < 0.0 || angle > 360.0) {
		siril_log_message(_("Angle %lf not allowed. Should be between 0.0 and 360.0.\n"), angle);
		return CMD_ARG_ERROR;
	}

	for (int i = 2; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				interpolation = OPENCV_AREA;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			return CMD_ARG_ERROR;
		} else if (!strcmp(word[i], "-noclamp")) {
			clamp = FALSE;
		} else if (!strcmp(word[i], "-nocrop")) {
			if (has_selection) {
				siril_log_color_message(_("-nocrop option is not valid if a selection is active. Ignoring\n"), "red");
			} else crop = 0;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (angle == 0.0 || angle == 360.0) {
		siril_log_message(_("Angle is 0.0 or 360.0 degrees. Doing nothing...\n"));
		return CMD_ARG_ERROR;
	}
	if (angle != 90.0 && angle != 180.0 && angle != 270.0)
		image_cfa_warning_check();

	// Allocate parameters
	struct rotation_args *params = new_rotation_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	params->area = area;
	params->angle = angle;
	params->interpolation = interpolation;
	params->cropped = crop;
	params->clamp = clamp;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rotation_args(params);
		set_cursor_waiting(FALSE);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = rotation_image_hook;
	args->log_hook = rotation_log_hook;
	args->idle_function = NULL;
	args->description = _("Rotation");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_rotatepi(int nb){
	// Allocate parameters
	struct rotation_args *params = new_rotation_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->area = (rectangle){0, 0, gfit->rx, gfit->ry};
	params->angle = 180.0;
	params->interpolation = -1;  // Fast rotation
	params->cropped = 0;
	params->clamp = FALSE;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rotation_args(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = rotation_image_hook;
	args->log_hook = rotation_log_hook;
	args->idle_function = NULL;
	args->description = _("Rotation 180°");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = 1;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_rgradient(int nb) {
	if (gfit->orig_bitpix == BYTE_IMG) {
		siril_log_color_message(_("This process cannot be applied to 8b images\n"), "red");
		return CMD_INVALID_IMAGE;
	}
	// Parse arguments
	gchar *endxc, *endyc, *enddR, *endda;
	double xc = g_ascii_strtod(word[1], &endxc);
	double yc = g_ascii_strtod(word[2], &endyc);
	double dR = g_ascii_strtod(word[3], &enddR);
	double da = g_ascii_strtod(word[4], &endda);
	gboolean mask_aware = FALSE;

	// Validate arguments
	if (endxc == word[1] || endyc == word[2] || enddR == word[3] || endda == word[4]) {
		siril_log_message(_("Invalid numeric arguments. Usage: rgradient xc yc dR da\n"));
		return CMD_ARG_ERROR;
	}
	if (xc >= gfit->rx || yc >= gfit->ry) {
		siril_log_message(_("The coordinates cannot be greater than the size of the image. "
			"Please change their values and retry.\n"));
		return CMD_ARG_ERROR;
	}

	// Check for optional -mask flag
	if (word[5] && !strcmp(word[5], "-mask")) {
		mask_aware = TRUE;
	}

	image_cfa_warning_check();
	// Allocate parameters
	struct rgradient_data *params = new_rgradient_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	// Set parameters
	params->xc = xc;
	params->yc = yc;
	params->dR = dR;
	params->da = da;
	params->fit = gfit;
	params->verbose = TRUE;
	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rgradient_data(params);
		return CMD_ALLOC_ERROR;
	}
	// Set up generic_img_args
	args->fit = gfit;
	args->mem_ratio = 3.0f; // Need memory for two temporary images
	args->image_hook = rgradient_image_hook;
	args->log_hook = rgradient_log_hook;
	args->description = _("Rotational Gradient");
	args->verbose = TRUE;
	args->user = params;
	args->mask_aware = mask_aware;
	args->max_threads = com.max_thread;
	args->command_updates_gfit = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_set(int nb) {
	gboolean is_get = word[0][0] == 'g' || word[0][0] == 'G';
	char *input = word[1];
	if (input[0] == '-') {
		if (is_get) {
			if (!strcmp(input, "-a"))
				return print_all_settings(FALSE);
			if (!strcmp(input, "-A"))
				return print_all_settings(TRUE);
			else return CMD_ARG_ERROR;
		}
		else {
			if (g_str_has_prefix(input, "-import=")) {
				char *filename = g_shell_unquote(input + 8, NULL);
				if (!filename || filename[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), input);
					return CMD_ARG_ERROR;
				}
				return readinitfile(filename);
			}
			else return CMD_ARG_ERROR;
		}
	}

	/* parsing a single variable command */
	int sep, len = strlen(input);
	for (sep = 1; sep < len; sep++)
		if (input[sep] == '.')
			break;
	if (sep == len) {
		siril_log_message("syntax: group.key=value\n");
		return CMD_ARG_ERROR;
	}
	input[sep] = '\0';
	if (is_get) {
		print_settings_key(input, input+sep+1, FALSE);
	} else {
		/* set */
		char fakefile[1024];
		int filelen = snprintf(fakefile, 1024, "[%s]\n%s\n", input, input+sep+1);
		GKeyFile *kf = g_key_file_new();
		g_key_file_load_from_data(kf, fakefile, filelen, G_KEY_FILE_NONE, NULL);
		return read_keyfile(kf) == 0;
	}
	return 0;
}

int process_set_mag(int nb) {
	gchar *end;
	double mag_reference = g_ascii_strtod(word[1], &end);
	if (end == word[1]) {
		siril_log_message(_("Wrong magnitude value.\n"));
		return CMD_ARG_ERROR;
	}

	gboolean found = FALSE;
	double mag = 0.0;
	if (gui.qphot) {
		mag = gui.qphot->phot->mag;
		found = TRUE;
	} else {
		if (com.selection.w > 300 || com.selection.h > 300){
			siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
			return CMD_SELECTION_ERROR;
		}
		if (com.selection.w <= 0 || com.selection.h <= 0){
			siril_log_message(_("Select an area first\n"));
			return CMD_SELECTION_ERROR;
		}
		psf_error error;
		struct phot_config *ps = phot_set_adjusted_for_image(gfit);
		psf_star *result = psf_get_minimisation(gfit, select_vport(gui.cvport), &com.selection, TRUE, FALSE, ps, TRUE, com.pref.starfinder_conf.profile, &error);
		free(ps);
		if (result && result->phot_is_valid && error == PSF_NO_ERR) {
			found = TRUE;
			mag = result->phot->mag;
		}
		else
			siril_log_message(_("PSF minimisation failed with error %d\n"), error);
		free_psf(result);
	}
	if (found) {
		com.magOffset = mag_reference - mag;
		siril_log_message(_(
					"Relative magnitude: %.3lf, "
					"True reduced magnitude: %.3lf, "
					"Offset: %.3lf\n"
				   ), mag, mag_reference, com.magOffset);
	}
	return CMD_OK;
}

int process_set_photometry(int nb) {
	if (nb > 0) {
		double inner = -1.0, outer = -1.0, aperture = -1.0, gain = -1.0, force = -1.0;
		int min = -65536, max = -1;	//, force = -1;
		gboolean error = FALSE;
		for (int i = 1; i < nb; i++) {
			char *arg = word[i], *end;
			if (!word[i])
				break;
			if (g_str_has_prefix(arg, "-inner=")) {
				arg += 7;
				inner = g_ascii_strtod(arg, &end);
				if (arg == end) error = TRUE;
				else if (inner < 1.0) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-outer=")) {
				arg += 7;
				outer = g_ascii_strtod(arg, &end);
				if (arg == end) error = TRUE;
				else if (outer < 2.0) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-aperture=")) {
				arg += 10;
				aperture = g_ascii_strtod(arg, &end);
				if (arg == end) error = TRUE;
				else if (aperture < 1.0) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-gain=")) {
				arg += 6;
				gain = g_ascii_strtod(arg, &end);
				if (arg == end) error = TRUE;
				else if (gain <= 0.0) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-min_val=")) {
				arg += 9;
				min = g_ascii_strtoll(arg, &end, 10);
				if (arg == end) error = TRUE;
				else if (min >= 65535) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-max_val=")) {
				arg += 9;
				max = g_ascii_strtoull(arg, &end, 10);
				if (arg == end) error = TRUE;
				else if (max == 0 || max > 65535) error = TRUE;
			}
			else if (g_str_has_prefix(arg, "-dyn_ratio=")) {
				arg += 11;
				force = g_ascii_strtod(arg, &end);
				if (arg == end) error = TRUE;
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
		}

		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
			return CMD_ARG_ERROR;
		}
		if (inner > 0) {
			if (outer > 0) {
				if (outer > inner) {
					com.pref.phot_set.inner = inner;
					com.pref.phot_set.outer = outer;
				} else {
					siril_log_message(_("Inner radius value must be less than outer. Please change the value."));
					error = TRUE;
				}
			} else {
				if (com.pref.phot_set.outer > inner) {
					com.pref.phot_set.inner = inner;
				} else {
					siril_log_message(_("Inner radius value must be less than outer. Please change the value."));
					error = TRUE;
				}
			}
		}
		else if (outer > 0) {
			if (outer > com.pref.phot_set.inner) {
				com.pref.phot_set.outer = outer;
			} else {
					siril_log_message(_("Inner radius value must be less than outer. Please change the value."));
					error = TRUE;
				}
		}
		if (aperture > 0.0)
			com.pref.phot_set.aperture = aperture;
		if (force && force >= 1.0 && force <= 5.0) {
			com.pref.phot_set.force_radius = FALSE;
			com.pref.phot_set.auto_aperture_factor = (double)force;
		} else {
			com.pref.phot_set.force_radius = TRUE;
		}

		if (gain > 0.0)
			com.pref.phot_set.gain = gain;
		if (min >= -65536) {
			if (max > 0) {
				if (min < max) {
					com.pref.phot_set.minval = (double)min;
					com.pref.phot_set.maxval = (double)max;
				} else {
					siril_log_message(_("minimum value must be smaller than the maximum\n"));
					error = TRUE;
				}
			} else {
				if ((double)min < com.pref.phot_set.maxval) {
					com.pref.phot_set.minval = (double)min;
				} else {
					siril_log_message(_("minimum value must be smaller than the maximum\n"));
					error = TRUE;
				}
			}
		} else if (max > 0) {
			if (com.pref.phot_set.minval < max) {
				com.pref.phot_set.maxval = max;
			} else {
				siril_log_message(_("minimum value must be smaller than the maximum\n"));
				error = TRUE;
			}
		}

		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
			return CMD_ARG_ERROR;
		}
	}
	siril_log_message(_("Local background annulus radius: %.1f to %.1f, %s: %.1f (%s), camera conversion gain: %f e-/ADU, using pixels with values ]%f, %f[\n"),
			com.pref.phot_set.inner,
			com.pref.phot_set.outer,
			com.pref.phot_set.force_radius ? _("aperture") : _("dynamic aperture"),
			com.pref.phot_set.force_radius ? com.pref.phot_set.aperture : com.pref.phot_set.auto_aperture_factor,
			com.pref.phot_set.force_radius ? _("forced") : _("times the half-FWHM"),
			com.pref.phot_set.gain,
			com.pref.phot_set.minval,
			com.pref.phot_set.maxval);
	return CMD_OK;
}

int process_set_ref(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	gchar *end;

	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}
	int n = g_ascii_strtoull(word[2], &end, 10) - 1;
	if (end == word[2] || n < 0 || n >= seq->number) {
		siril_log_message(_("The reference image must be set between 1 and %d\n"), seq->number);
		return CMD_ARG_ERROR;
	}

	seq->reference_image = n;
	// a reference image should not be excluded to avoid confusion
	if (!seq->imgparam[seq->reference_image].incl) {
		seq->imgparam[seq->reference_image].incl = TRUE;
	}

	writeseqfile(seq);
	if (check_seq_is_comseq(seq)) {
		fix_selnum(&com.seq, FALSE);
		seq_load_image(&com.seq, n, TRUE);
		update_stack_interface(TRUE);
		update_reg_interface(FALSE);
		adjust_sellabel();
		drawPlot();
	} else {
		free_sequence(seq, FALSE);
	}

	return CMD_OK;
}

int process_unset_mag(int nb) {
	com.magOffset = 0.0;
	return CMD_OK;
}

int process_set_mag_seq(int nb) {
	gchar *end;
	double mag = g_ascii_strtod(word[1], &end);
	if (end == word[1]) {
		siril_log_message(_("Wrong magnitude value.\n"));
		return CMD_ARG_ERROR;
	}

	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++);
	com.seq.reference_star = i - 1;
	if (i == 0) {
		siril_log_message(_("Run a PSF for the sequence first (see seqpsf)\n"));
		return CMD_NEED_INIT_FIRST;
	}
	com.seq.reference_mag = mag;
	siril_log_message(_("Reference magnitude has been set for star %d to %f and will be computed for each image\n"), i - 1, mag);
	drawPlot();
	return CMD_OK;
}

int process_set_ext(int nb) {
	if (word[1]) {
		if ((g_ascii_strncasecmp(word[1], "fit", 3))
				&& (g_ascii_strncasecmp(word[1], "fts", 3))
				&& (g_ascii_strncasecmp(word[1], "fits", 4))) {
			siril_log_message(_("FITS extension unknown: %s\n"), word[1]);
			return CMD_ARG_ERROR;
		}

		g_free(com.pref.ext);
		gchar *lower = g_ascii_strdown(word[1], strlen(word[1]));
		com.pref.ext = g_strdup_printf(".%s", lower);
		g_free(lower);
		writeinitfile();
	}

	return CMD_OK;
}

int process_set_findstar(int nb) {
	int startoptargs = 1;
	double sigma = com.pref.starfinder_conf.sigma;
	double minbeta = com.pref.starfinder_conf.min_beta;
	double roundness = com.pref.starfinder_conf.roundness;
	int radius = com.pref.starfinder_conf.radius;
	double focal_length = com.pref.starfinder_conf.focal_length;
	double pixel_size_x = com.pref.starfinder_conf.pixel_size_x;
	gboolean relax_checks = com.pref.starfinder_conf.relax_checks;
	int convergence = com.pref.starfinder_conf.convergence;
	starprofile profile = com.pref.starfinder_conf.profile;
	double minA = com.pref.starfinder_conf.min_A;
	double maxA = com.pref.starfinder_conf.max_A;
	double maxr = com.pref.starfinder_conf.max_r;
	gchar *end;

	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			char *current = word[i], *value;
			if (current) {
				if (g_str_has_prefix(current, "-radius=")) {
					value = current + 8;
					radius = g_ascii_strtoull(value, &end, 10);
					if (end == value || radius < 3 || radius > 50) {
						siril_log_message(_("Wrong parameter values. Radius must be between 3 and 50, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-sigma=")) {
					value = current + 7;
					sigma = g_ascii_strtod(value, &end);
					if (end == value || sigma < 0.05) {
						siril_log_message(_("Wrong parameter values. Sigma must be greater than 0.05, aborting\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-minbeta=")) {
					value = current + 9;
					minbeta = g_ascii_strtod(value, &end);
					if (end == value || minbeta < 0.0 || minbeta >= MOFFAT_BETA_UBOUND) {
						siril_log_message(_("Wrong parameter values. Minimum beta must be greater than or equal to 0.0 and less than %.0f, aborting\n"), MOFFAT_BETA_UBOUND);
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-gaussian")) {
					profile = PSF_GAUSSIAN;
				} else if (g_str_has_prefix(current, "-moffat")) {
					profile = PSF_MOFFAT_BFREE;
				} else if (g_str_has_prefix(current, "-roundness=")) {
					value = current + 11;
					roundness = g_ascii_strtod(value, &end);
					if (end == value || roundness < 0.0 || roundness > 0.95) {
						siril_log_message(_("Wrong parameter values. Roundness must be between 0 and 0.95, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-focal=")) {
					value = current + 7;
					focal_length = g_ascii_strtod(value, &end);
					if (end == value || focal_length < 0) {
						siril_log_message(_("Wrong parameter values. Focal length must be greater than 0, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-pixelsize=")) {
					value = current + 11;
					pixel_size_x = g_ascii_strtod(value, &end);
					if (end == value || pixel_size_x < 0) {
						siril_log_message(_("Wrong parameter values. Pixel size must be greater than 0, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-relax=")) {
					value = current + 7;
					if (!(g_ascii_strcasecmp(value, "on"))) relax_checks = TRUE;
					else if (!(g_ascii_strcasecmp(value, "off"))) relax_checks = FALSE;
					else {
						siril_log_message(_("Wrong parameter values. Relax must be set to on or off, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-convergence=")) {
					value = current + 13;
					convergence = g_ascii_strtoull(value, &end, 10);
					if (end == value || convergence < 1 || convergence > 3) {
						siril_log_message(_("Wrong parameter values. Convergence must be between 1 and 3, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-minA=")) {
					value = current + 6;
					minA = g_ascii_strtod(value, &end);
					if (end == value) {
						siril_log_message(_("Wrong parameter value %s.\n"), current);
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-maxA=")) {
					value = current + 6;
					maxA = g_ascii_strtod(value, &end);
					if (end == value) {
						siril_log_message(_("Wrong parameter value %s.\n"), current);
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(current, "-maxR=")) {
					value = current + 6;
					maxr = g_ascii_strtod(value, &end);
					if (end == value || maxr <= 0.0 || maxr > 1.0) {
						siril_log_message(_("Wrong parameter value %s.\n"), current);
						return CMD_ARG_ERROR;
					}
				} else if (!g_ascii_strcasecmp(current, "reset")) {
					siril_log_message(_("Resetting findstar parameters to default values.\n"));
					sigma = 1.;
					roundness = 0.5;
					radius = DEF_BOX_RADIUS;
					focal_length = 0.;
					pixel_size_x = 0.;
					relax_checks = FALSE;
					convergence = 1;
					profile = PSF_GAUSSIAN;
					minbeta = 1.5;
					minA = 0.0;
					maxA = 0.0;
					maxr = 1.0;
				} else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
			}
		}
	}

	siril_log_message(_("sigma = %3.2f\n"), sigma);
	com.pref.starfinder_conf.sigma = sigma;
	if (maxr != 1.0 && maxr > roundness)
		siril_log_message(_("roundness range = [%f, %f]\n"), roundness, maxr);
	else {
		siril_log_message(_("roundness = %3.2f\n"), roundness);
		maxr = 1.0;
	}
	com.pref.starfinder_conf.roundness = roundness;
	com.pref.starfinder_conf.max_r = maxr;
	siril_log_message(_("radius = %d\n"), radius);
	com.pref.starfinder_conf.radius = radius;
	if ((minA > 0.0 || maxA > 0.0) && minA < maxA)
		siril_log_message(_("amplitude range = [%f, %f]\n"), minA, maxA);
	else {
		siril_log_message(_("amplitude range unlimited\n"));
		minA = 0.0; maxA = 0.0;
	}
	com.pref.starfinder_conf.min_A = minA;
	com.pref.starfinder_conf.max_A = maxA;
	siril_log_message(_("convergence = %d\n"), convergence);
	com.pref.starfinder_conf.convergence = convergence;
	siril_log_message(_("focal = %3.1f\n"), focal_length);
	com.pref.starfinder_conf.focal_length = focal_length;
	siril_log_message(_("pixelsize = %3.2f\n"), pixel_size_x);
	com.pref.starfinder_conf.pixel_size_x = pixel_size_x;
	siril_log_message(_("relax = %s\n"), (relax_checks) ? "on" : "off");
	com.pref.starfinder_conf.relax_checks = relax_checks;
	siril_log_message(_("profile = %s\n"), (profile == PSF_GAUSSIAN) ? "Gaussian" : "Moffat");
	com.pref.starfinder_conf.profile = profile;
	if (profile != PSF_GAUSSIAN) {
		siril_log_message(_("minimum beta = %3.2f\n"), minbeta);
		com.pref.starfinder_conf.min_beta = minbeta;
	}
	return CMD_OK;
}

int process_unset_mag_seq(int nb) {
	com.seq.reference_star = -1;
	com.seq.reference_mag = -1001.0;
	siril_log_message(_("Reference magnitude unset for sequence\n"));
	drawPlot();
	return CMD_OK;
}

int process_pm(int nb) {
	/* First we want to replace all variable by filename if exist. Return error if not
	 * Variables start and end by $ token.
	 */
	gchar *expression = g_shell_unquote(word[1], NULL);
	gchar *next, *cur;
	int count = 0;
	float min = -1.f;
	float max = -1.f;
	gboolean do_sum = TRUE;
	gboolean has_gfit = FALSE;

	// first replace $T with "gfit"
	GRegex *regex = g_regex_new("\\$T(?![A-Za-z0-9_])", 0, 0, NULL);
	gchar *cleaned_expression = g_regex_replace(regex, expression, -1, 0, "gfit", 0, NULL);
	g_regex_unref(regex);

	// Check if a replacement was made to set some flags and check if an image is really loaded
	if (g_strcmp0(expression, cleaned_expression) != 0) {
		if (!single_image_is_loaded()) {
			g_free(cleaned_expression);
			return CMD_ARG_ERROR;
		}
		has_gfit = TRUE;
	}

    g_free(expression);

    /* cleaned_expression is the expression we need to keep to work on it
     * Let's free it at the end
     */
    expression = g_strdup(cleaned_expression);

    /* as $T has already be changed, the number of tokens '$' must be even
     * Let's count them
     */
    cur = expression;
	while ((next = strchr(cur, '$')) != NULL) {
		cur = next + 1;
		count++;
	}

	if (count == 0 && !has_gfit) {
		siril_log_message(_("You need to add at least a loaded image or one image as variable. Use $ tokens to surround the file names .\n"));
		g_free(expression);
		g_free(cleaned_expression);
		return CMD_ARG_ERROR;
	} else if (count % 2 != 0) {
		siril_log_message(_("There is an unmatched $. Please check the expression.\n"));
		g_free(expression);
		g_free(cleaned_expression);
		return CMD_ARG_ERROR;
	}

	/* parse rescale and nosum options if they exist */
	if (nb > 1) {
		for (int i = 2; i < 5; i++) {
			if (!g_strcmp0(word[i], "-rescale")) {
				if (nb == 5) {
					gchar *end;
					min = g_ascii_strtod(word[i + 1], &end);
					if (end == word[i + 1] || min < 0 || min > 1) {
						siril_log_message(_("Rescale can only be done in the [0, 1] range.\n"));
						g_free(expression);
						g_free(cleaned_expression);
						return CMD_ARG_ERROR;
					}
					max = g_ascii_strtod(word[i + 2], &end);
					if (end == word[i + 2] || max < 0 || max > 1) {
						siril_log_message(_("Rescale can only be done in the [0, 1] range.\n"));
						g_free(expression);
						g_free(cleaned_expression);
						return CMD_ARG_ERROR;
					}
				} else {
					min = 0.f;
					max = 1.f;
				}
			} else if (!g_strcmp0(word[i], "-nosum")) {
				do_sum = FALSE;
			}
		}
	}

	struct pixel_math_data *args = calloc(1, sizeof(struct pixel_math_data));
	args->nb_rows = count / 2; // this is the number of variable
	args->varname = calloc(args->nb_rows, sizeof(gchar *));

	cur = expression;

	/* List all variables */
	char *start = cur;
	char *end = cur;
	gboolean first = TRUE;
	int i = 0;
	while (*cur) {
		if (*cur == '$') {
			if (first) {
				start = cur;
				first = FALSE;
			} else {
				end = cur;
				first = TRUE;
			}
		}
		if (start < end && *start) {
			*end = 0;
			gboolean found = FALSE;
			for (int j = 0; j < i; j++) {
/*				gchar *test = g_strrstr(args->varname[j], start + 1);
				if (test)*/
				if (!g_strcmp0(args->varname[j], start + 1))
					found = TRUE;
			}
			if (!found)
				args->varname[i++] = g_strdup(start + 1);
			start = cur = end;
		}
		cur++;
	}
	args->nb_rows = i; // this is the final number of variables
	args->varname = realloc(args->varname, i * sizeof(gchar *));

	int width = -1;
	int height = -1;
	int channel = -1;

	for (int j = 0; j < args->nb_rows; j++) {
		int w, h, c;
		if (args->varname && load_pm_var(args->varname[j], j, &w, &h, &c)) {
			if (j > 0)
				free_pm_var(j - 1);
			g_free(expression);
			g_free(cleaned_expression);
			free(args->varname);
			free(args);
			return CMD_INVALID_IMAGE;
		}

		if (channel == -1) {
			width = w;
			height = h;
			channel = c;
		} else {
			if (w != width || height != h || channel != c) {
				siril_log_message(_("Image must have same dimension\n"));
				free_pm_var(args->nb_rows);
				g_free(expression);
				g_free(cleaned_expression);
				free(args->varname);
				free(args);
				return CMD_INVALID_IMAGE;
			}
		}
	}

	/* gfit image MUST have same size of the others */
	if (has_gfit) {
		if (width == -1) {
			width = gfit->rx;
			height = gfit->ry;
			channel = gfit->naxes[2];
		}
		if (gfit->rx != width || height != gfit->ry || channel != gfit->naxes[2]) {
			siril_log_message(_("Image must have same dimension\n"));
			free_pm_var(args->nb_rows);
			g_free(expression);
			g_free(cleaned_expression);
			free(args->varname);
			free(args);
			return CMD_INVALID_IMAGE;
		}
	}

	/* OK. Now all prerequisites are ok. Let's remove the tokens
	 * cleaned_expression can be freed right after its last used
	 */
	g_free(expression);
	expression = g_shell_unquote(cleaned_expression, NULL);
	g_free(cleaned_expression);

	/* We must now replace the original variable names between the $ signs in the expression with
	 * the new generic variable names.
	 * This ensures that the variable names in the expression passed to pm match the variable names
	 * stored in args->varname
	 */
	gchar **chunks = g_strsplit(expression, "$", count + 1);
	for (int i = 0, j = 1; i < count / 2; i++) {
		int idx = 0;
		int k;
		for (k = 0; k < args->nb_rows; k++) {
			if (!g_strcmp0(chunks[2 * i + 1], args->varname[k])) {
				idx = k + 1;
				break;
			}
		}
		if (idx != k + 1) idx = j;
		g_free(chunks[2 * i + 1]);
		chunks[2 * i + 1] = g_strdup_printf("var_%d", idx);
		j++;
	}
	g_free(expression);
	expression = g_strjoinv(NULL, chunks);
	g_strfreev(chunks);

	/* Rewrite the variable names to var_1, var_2 etc. now the files are loaded.
	 * This avoids conflicts where characters are permitted in filenames but cannot
	 * be used in pixelmath variable names.
	 * We will amend the expression to match below.
	 */
	for (int j = 0; j < args->nb_rows; j++) {
		g_free(args->varname[j]);
		args->varname[j] = g_strdup_printf("var_%d", j + 1);
	}

	remove_spaces_from_str(expression);

	fits *fit = NULL;
	if (new_fit_image(&fit, width, height, channel, com.pref.force_16bit ? DATA_USHORT : DATA_FLOAT)) {
		free_pm_var(args->nb_rows);
		free(args->varname);
		free(args);
		return CMD_GENERIC_ERROR;
	}

	args->expression1 = expression;
	args->expression2 = NULL;
	args->expression3 = NULL;
	args->single_rgb = TRUE;
	args->fit = fit;
	args->ret = 0;
	args->from_ui = FALSE;
	args->do_sum = do_sum;
	args->has_gfit = has_gfit;
	if (min >= 0.f) {
		args->rescale = TRUE;
		args->min = min;
		args->max = max;
	} else {
		args->rescale = FALSE;
	}

	if (!start_in_new_thread(apply_pixel_math_operation, args)) {
		g_free(args->expression1);
		clearfits(fit);
		g_strfreev(args->varname);
		free(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_psf(int nb){
	if (com.selection.w > 300 || com.selection.h > 300){
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return CMD_SELECTION_ERROR;
	}
	if (com.selection.w <= 0 || com.selection.h <= 0){
		siril_log_message(_("Select an area first\n"));
		return CMD_SELECTION_ERROR;
	}

	if (gfit->naxes[2] > 1 && nb == 1 && com.headless) {
		siril_log_color_message(_("Please display the channel on which you want to compute the PSF or use -channel argument\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	int channel = com.headless ? 0 : select_vport(gui.cvport);
	if (nb == 2) {
		gchar *next;
		channel = g_ascii_strtoull(word[1], &next, 10);
		if (word[1] == next || channel > (int)gfit->naxes[2]) {
			siril_log_message(_("Please provide the channel number starting from 0 for red\n"));
			return CMD_ARG_ERROR;
		}
	}

	starprofile profile = com.pref.starfinder_conf.profile;
	psf_error error = PSF_NO_ERR;
	struct phot_config *ps = phot_set_adjusted_for_image(gfit);
	psf_star *result = psf_get_minimisation(gfit, channel, &com.selection, TRUE, FALSE, ps, TRUE, profile, &error);
	free(ps);
	if (result) {
		gchar *str = format_psf_result(result, &com.selection, gfit, NULL);
		siril_log_message("%s\n", str);
		g_free(str);
	}
	free_psf(result);
	return CMD_OK;
}

int process_seq_tilt(int nb) {
	gboolean draw_polygon = FALSE;

	sequence *seq;
	seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
		draw_polygon = TRUE;
	}
	// through GUI, in case the specified sequence is not the loaded sequence
	// load it before running
	if (!com.script && !com.python_script && !check_seq_is_comseq(seq)) {
		execute_idle_and_wait_for_it(set_seq, seq->seqname);
		free_sequence(seq, TRUE);
		seq = &com.seq;
		draw_polygon = TRUE;
	}

	struct tilt_data *args = calloc(sizeof(struct tilt_data), 1);
	args->seq = seq;
	args->draw_polygon = draw_polygon;

	apply_tilt_to_sequence(args);
	return CMD_OK;
}

/* com.pref.phot_set.outer needs to be set at this point */
static int parse_star_position_arg(char *arg, sequence *seq, fits *first, rectangle *seqpsf_area, gchar **target_descr) {
	rectangle area;
	if (g_str_has_prefix(arg, "-at=") || g_str_has_prefix(arg, "-refat=")) {
		gchar *value, *end;
		if (arg[1] == 'a')
			value = arg + 4;
		else value = arg + 7;
		int x = g_ascii_strtoull(value, &end, 10);
		if (end == value) return CMD_ARG_ERROR;
		if (*end != ',') return CMD_ARG_ERROR;
		end++;
		value = end;
		int y = g_ascii_strtoull(value, &end, 10);
		if (end == value) return CMD_ARG_ERROR;

		double start = 1.5 * com.pref.phot_set.outer;
		double size = 3 * com.pref.phot_set.outer;

		area.x = x - start;
		area.y = y - start;
		area.w = size;
		area.h = size;
		if (area.x < 0 || area.y < 0 ||
				area.h <= 0 || area.w <= 0 ||
				area.x + area.w >= first->rx ||
				area.y + area.h >= first->ry) {
			siril_log_message(_("The given coordinates are not in the image, aborting\n"));
			return CMD_ARG_ERROR;
		}

		if (enforce_area_in_image(&area, seq, 0)) {
			siril_log_message(_("Selection was not completely inside the first image of the sequence, aborting.\n"));
			return CMD_SELECTION_ERROR;
		}
		if (target_descr && arg[1] != 'a')
			*target_descr = g_strdup_printf("at %d, %d", x, y);
	}
	else if (g_str_has_prefix(arg, "-wcs=") || g_str_has_prefix(arg, "-refwcs=")) {
		char *value;
		if (arg[1] == 'w')
			value = arg + 5;
		else value = arg + 8;
		char *sep = strchr(value, ',');
		if (!sep) {
			siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
			return CMD_ARG_ERROR;
		}
		*sep++ = '\0';
		char *end;
		double ra = g_ascii_strtod(value, &end);
		if (end == value) {
			siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
			return CMD_ARG_ERROR;
		}
		double dec = g_ascii_strtod(sep, &end);
		if (end == sep) {
			siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
			return CMD_ARG_ERROR;
		}

		if (!has_wcs(first)) {
			siril_log_message(_("The reference image of the sequence does not have the WCS information required for star selection by equatorial coordinates\n"));
			return CMD_FOR_PLATE_SOLVED;
		}
		double fx, fy;
		if (wcs2pix(first, ra, dec, &fx, &fy)) {
			siril_log_message(_("The given coordinates are not in the image, aborting\n"));
			return CMD_ARG_ERROR;
		}
		double dx, dy;
		siril_to_display(fx, fy, &dx, &dy, first->ry);
		double start = 1.5 * com.pref.phot_set.outer;
		double size = 3 * com.pref.phot_set.outer;
		area.x = dx - start;
		area.y = dy - start;
		area.w = size;
		area.h = size;
		if (area.x < 0 || area.y < 0 ||
				area.h <= 0 || area.w <= 0 ||
				area.x + area.w >= first->rx ||
				area.y + area.h >= first->ry) {
			siril_log_message(_("The given coordinates are not in the image, aborting\n"));
			return CMD_ARG_ERROR;
		}
		siril_log_message(_("Coordinates of the star: %.1f, %.1f\n"), dx, dy);
		if (target_descr && arg[1] != 'a')
			*target_descr = g_strdup_printf("at %f, %f", ra, dec);
	}
	else {
		siril_log_message(_("Invalid argument %s, aborting.\n"), arg);
		return CMD_ARG_ERROR;
	}
	*seqpsf_area = area;
	return CMD_OK;
}

gboolean get_followstar_idle(gpointer user_data) {
	framing_mode *framing = (framing_mode*) user_data;
	GtkToggleButton *follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
	if (gtk_toggle_button_get_active(follow))
		*framing = FOLLOW_STAR_FRAME;
	// no need to have an else as framing is already initiated by the caller
	return FALSE;
}

/* seqpsf sequencename [channel] [{ -at=x,y | -wcs=ra,dec }] [-followstar] */
int process_seq_psf(int nb) {
	// Must have at least sequence name argument
	if (nb < 2) {
		siril_log_message(_("Sequence name is required as first argument (use '.' for current sequence)\n"));
		return CMD_ARG_ERROR;
	}

	sequence *seq = NULL;
	int layer = -1;
	gboolean has_area = FALSE;
	rectangle area = { 0 };
	if (com.selection.w > 0 && com.selection.h > 0) {
		area = com.selection;
		has_area = TRUE;
	}
	gboolean followstar_set = FALSE;

	// First argument is always sequence name
	seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	gboolean use_current_seq = check_seq_is_comseq(seq);
	if (!com.script && !com.python_script && !use_current_seq) {
		execute_idle_and_wait_for_it(set_seq, seq->seqname);
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	// Parse remaining arguments in any order
	for (int i = 2; i < nb; i++) {
		if (strcmp(word[i], "-followstar") == 0) {
			followstar_set = TRUE;
		} else if (layer == -1) {
			// Try to parse as layer number
			gchar *end;
			int potential_layer = g_ascii_strtoull(word[i], &end, 10);
			if (end != word[i] && potential_layer < seq->nb_layers) {
				layer = potential_layer;
			} else {
				// Not a valid layer number, try as star position
				if (!has_area) {
					fits first = { 0 };
					if (use_current_seq) {
						// For current sequence, we may not need to read metadata
						// if we're using current selection
						if (area.w > 0 && area.h > 0) {
							has_area = TRUE;
							continue;
						}
					}

					if (seq_read_frame_metadata(seq, 0, &first)) {
						if (seq != &com.seq) free_sequence(seq, TRUE);
						return CMD_GENERIC_ERROR;
					}
					seq->current = 0;

					if (parse_star_position_arg(word[i], seq, &first, &area, NULL) == 0) {
						has_area = TRUE;
						clearfits(&first);
					} else {
						clearfits(&first);
						siril_log_message(_("Invalid argument: %s\n"), word[i]);
						if (seq != &com.seq) free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
				} else {
					siril_log_message(_("Invalid argument: %s\n"), word[i]);
					if (seq != &com.seq) free_sequence(seq, TRUE);
					return CMD_ARG_ERROR;
				}
			}
		} else if (!has_area) {
			// Try to parse as star position
			fits first = { 0 };
			if (use_current_seq) {
				// For current sequence, we may not need to read metadata
				// if we're using current selection
				if (area.w > 0 && area.h > 0) {
					has_area = TRUE;
					continue;
				}
			}
			if (seq_read_frame_metadata(seq, 0, &first)) {
				if (seq != &com.seq)
					free_sequence(seq, TRUE);
				return CMD_GENERIC_ERROR;
			}
			seq->current = 0;

			if (parse_star_position_arg(word[i], seq, &first, &area, NULL) == 0) {
				has_area = TRUE;
				clearfits(&first);
			} else {
				clearfits(&first);
				siril_log_message(_("Invalid argument: %s\n"), word[i]);
				if (seq != &com.seq) free_sequence(seq, TRUE);
				return CMD_ARG_ERROR;
			}
		} else {
			siril_log_message(_("Invalid argument: %s\n"), word[i]);
			if (seq != &com.seq) free_sequence(seq, TRUE);
			return CMD_ARG_ERROR;
		}
	}

	// Set defaults for missing arguments
	if (layer == -1) {
		if (use_current_seq) {
			layer = select_vport(gui.cvport);
		} else {
			layer = 0; // Default to first layer
		}
	}

	if (!has_area) {
		if (use_current_seq && com.selection.w > 0 && com.selection.h > 0) {
			area = com.selection;
			has_area = TRUE;
		} else {
			if (com.script) {
				siril_log_message(_("Layer and star position arguments are required when called from a script\n"));
				if (seq != &com.seq) free_sequence(seq, TRUE);
				return CMD_ARG_ERROR;
			}
			if (!com.headless) {
				siril_log_message(_("Layer and star position arguments are required, or make a selection around a star\n"));
				if (seq != &com.seq) free_sequence(seq, TRUE);
				return CMD_ARG_ERROR;
			}
		}
	}

	// Validate layer
	if (layer >= seq->nb_layers) {
		siril_log_message(_("PSF cannot be computed on channel %d for this sequence of %d channels\n"), layer, seq->nb_layers);
		if (seq != &com.seq) free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	// Validate selection area
	if (has_area) {
		if (area.w > 300 || area.h > 300) {
			siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
			if (seq != &com.seq) free_sequence(seq, TRUE);
			return CMD_SELECTION_ERROR;
		}
		if (area.w <= 0 || area.h <= 0) {
			siril_log_message(_("Select an area first\n"));
			if (seq != &com.seq) free_sequence(seq, TRUE);
			return CMD_SELECTION_ERROR;
		}
		com.selection = area;
	}

	// Set framing mode
	framing_mode framing = ORIGINAL_FRAME;
	if (followstar_set) {
		framing = FOLLOW_STAR_FRAME;
	} else {
		// Try registered frame first if available
		if (seq->regparam[layer]) {
			framing = REGISTERED_FRAME;
		}
	}

	siril_log_message(_("Running the PSF on the sequence, layer %d\n"), layer);
	int retval = seqpsf(seq, layer, FALSE, FALSE, FALSE, framing, TRUE, (com.script || com.python_script)) ? CMD_GENERIC_ERROR : CMD_OK;
	return retval;
}

// light_curve sequencename channel [-autoring] { -ninalist=file | [-auto] { -at=x,y | -wcs=ra,dec [-refat=x,y] [-refwcs=ra,dec] ... } }
int process_light_curve(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	gboolean has_GUI = FALSE;
	gboolean seq_has_wcs = FALSE;
	int refimage = sequence_find_refimage(seq);
	if (check_seq_is_comseq(seq)) { // we are in GUI
		free_sequence(seq, TRUE);
		seq = &com.seq;
		clear_all_photometry_and_plot(); // calls GTK+ code, but we're not in a script here
		init_plot_colors();
		has_GUI = TRUE;
		seq_has_wcs = has_wcs(gfit);
	} else { // we are in script or headless, loading the seq has loaded the ref image, we check if it has wcs info
		siril_debug_print("reference image in seqfile is %d\n", refimage);
		seq_has_wcs = sequence_ref_has_wcs(seq);
		seq->current = refimage;
	}

	gchar *end;
	int layer = g_ascii_strtoull(word[2], &end, 10);
	if (end == word[2] || layer >= seq->nb_layers) {
		siril_log_message(_("PSF cannot be computed on channel %d for this sequence of %d channels\n"), layer, seq->nb_layers);
		if (seq != &com.seq)
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	if (sequence_drifts(seq, layer, seq->rx / 4))
		return CMD_GENERIC_ERROR;

	int argidx = 3;
	gboolean autoring = FALSE;
	if (!g_strcmp0(word[argidx], "-autoring")) {
		autoring = TRUE;
		if (nb == 3) {
			siril_log_message(_("Missing argument to command\n"));
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			return CMD_ARG_ERROR;
		}
		argidx++;
	}

	fits reffit = { 0 };
	if (autoring) { // we detect the stars in the reference frame
		if (seq_read_frame(seq, refimage, &reffit, FALSE, -1)) {
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			return CMD_GENERIC_ERROR;
		}
		float fwhm = measure_image_FWHM(&reffit, layer);
		if (fwhm <= 0.0f) {
			siril_log_color_message(_("Could not find stars in the reference image, aborting.\n"), "red");
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			clearfits(&reffit);
			return CMD_GENERIC_ERROR;
		}
		/* automatic aperture is 2 times each star's FWHM, which can be
		 * significantly larger than the mean FWHM that we get here */
		com.pref.phot_set.inner = com.pref.phot_set.auto_inner_factor * fwhm;
		com.pref.phot_set.outer = com.pref.phot_set.auto_outer_factor * fwhm;
		siril_log_message(_("Setting inner and outer photometry ring radii to %.1f and %.1f (FWHM is %f)\n"),
				com.pref.phot_set.inner, com.pref.phot_set.outer, fwhm);
		clearfits(&reffit);
	} else {
		if (seq_read_frame_metadata(seq, refimage, &reffit)) {
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			return CMD_GENERIC_ERROR;
		}
		siril_log_message(_("Using preconfigured inner and outer photometry ring radii of %.1f and %.1f\n"),
				com.pref.phot_set.inner, com.pref.phot_set.outer);
	}
	g_assert(com.pref.phot_set.inner < com.pref.phot_set.outer);

	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		if (seq != &com.seq)
			free_sequence(seq, TRUE);
		return CMD_ALLOC_ERROR;
	}

	/* we have a sequence and channel, now we iterate with seqpsf on the list of stars */
	if (g_str_has_prefix(word[argidx], "-ninastars=")) {
		char *file = word[argidx] + 11;
		if (file[0] == '\0') {
			siril_log_message(_("Missing argument to %s, aborting.\n"), word[argidx]);
			free(args);
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			return CMD_GENERIC_ERROR;
		}
		if (!seq_has_wcs) {
			siril_log_message(_("No image in the sequence was found with the WCS information required for star selection by equatorial coordinates, plate solve the reference or the current\n"));
			free(args);
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			return CMD_FOR_PLATE_SOLVED;
		}
		if (parse_nina_stars_file_using_WCS(args, file, TRUE, TRUE, has_GUI ? gfit : &reffit)) {
			if (seq != &com.seq)
				free_sequence(seq, TRUE);
			free_light_curve_args(args);
			return CMD_GENERIC_ERROR;
		}
	} else {
		args->areas = malloc((nb - argidx) * sizeof(rectangle));
		for (int arg_index = argidx; arg_index < nb; arg_index++) {
			if (parse_star_position_arg(word[arg_index], seq, has_GUI ? gfit : &reffit,
						&args->areas[arg_index - argidx], &args->target_descr)) {
				free_light_curve_args(args);
				return CMD_ARG_ERROR;;
			}
		}
		args->nb = nb - argidx;
	}
	clearfits(&reffit);
	args->seq = seq;
	args->layer = layer;
	args->display_graph = has_GUI;
	siril_debug_print("starting PSF analysis of %d stars\n", args->nb);

	if (!start_in_new_thread(light_curve_worker, args)) {
		g_free(args->target_descr);
		free(args->areas);
		free(args->target_descr);
		free(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_seq_crop(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	rectangle area;

	int startoptargs = 6;

	if (nb >= startoptargs) {
		gchar *endx, *endy, *endw, *endh;

		area.x = g_ascii_strtoull(word[2], &endx, 10);
		area.y = g_ascii_strtoull(word[3], &endy, 10);
		area.w = g_ascii_strtoull(word[4], &endw, 10);
		area.h = g_ascii_strtoull(word[5], &endh, 10);

		if (endx == word[2] || endy == word[3] || area.x < 0 || area.y < 0) {
			siril_log_message(_("Crop: x and y must be positive values.\n"));
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
			return CMD_ARG_ERROR;
		}
		if (endw == word[4] || endh == word[5] || area.w <= 0 || area.h <= 0) {
			siril_log_message(_("Crop: width and height must be greater than 0.\n"));
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
			return CMD_ARG_ERROR;
		}
	} else {
		siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}
	gchar *end1, *end2;
	int rx = g_ascii_strtoull(word[4], &end1, 10);
	int ry = g_ascii_strtoull(word[5], &end2, 10);

	if (end1 == word[4] || end2 == word[5] || rx > seq->rx || ry > seq->ry) {
		siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"),
				seq->rx, seq->ry);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	struct crop_sequence_data *args = calloc(1, sizeof(struct crop_sequence_data));

	args->seq = seq;
	args->area = area;
	args->prefix = strdup("cropped_");

	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), current);
						free(args->prefix);
						free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					free(args->prefix);
					args->prefix = strdup(value);
				}
			}
		}
	}

	crop_sequence(args);
	return CMD_OK;
}

int process_seq_resample(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	struct scale_sequence_data *args = calloc(1, sizeof(struct scale_sequence_data));

	args->seq = seq;
	args->prefix = strdup("scaled_");
	args->scale = 1.0;
	args->interpolation = OPENCV_AREA;
	args->clamp = FALSE;

	for (int i = 2; i < nb; i++) {
		if (word[i]) {
			if (g_str_has_prefix(word[i], "-prefix=")) {
				char *current = word[i], *value;
				value = current + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					free(args->prefix);
					free(args);
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					return CMD_ARG_ERROR;
				}
				free(args->prefix);
				args->prefix = strdup(value);
			}
			else if (g_str_has_prefix(word[1], "-height=")) {
				char *current = word[i] + 8, *end;
				double toY = g_ascii_strtod(current, &end);
				if (toY < 10) {
					siril_log_message(_("Error: height cannot be less than 10, aborting.\n"));
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					free(args->prefix);
					return CMD_ARG_ERROR;
				}
				args->scale = toY / (double) seq->ry;
			}
			else if (g_str_has_prefix(word[1], "-width=")) {
				char *current = word[i] + 7, *end;
				double toX = g_ascii_strtod(current, &end);
				if (toX < 10) {
					siril_log_message(_("Error: width cannot be less than 10, aborting.\n"));
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					free(args->prefix);
					return CMD_ARG_ERROR;
				}
				args->scale = toX / (double) seq->rx;
			} else if (g_str_has_prefix(word[1], "-maxdim=")) {
				char *current = word[i] + 8, *end;
				double maxdim = g_ascii_strtoull(current, &end, 10);
				if (maxdim < 10) {
					siril_log_message(_("Error: max dimension cannot be less than 10, aborting.\n"));
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						free(args->prefix);
					return CMD_ARG_ERROR;
				}
				if (gfit->rx > gfit->ry) {
					args->scale = maxdim / (double) seq->rx;
				} else {
					args->scale = maxdim / (double) seq->ry;
				}
			} else if (g_str_has_prefix(word[i], "-scale=")) {
				char *current = word[i] + 7, *end;
				args->scale = g_ascii_strtod(current, &end);
				if (args->scale < 0.01) {
					siril_log_message(_("Error: scale cannot be less than 0.01, aborting.\n"));
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					free(args->prefix);
					return CMD_ARG_ERROR;
				}
				if (current == end) {
					siril_log_message(_("Error parsing argument to -scale= (%s), aborting.\n"), current);
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					free(args->prefix);
					return CMD_ARG_ERROR;
				}
			}
			else if (g_str_has_prefix(word[i], "-interp=")) {
				char *current = word[i], *value;
				value = current + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
				if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
					args->interpolation = OPENCV_NEAREST;
					args->clamp = FALSE;
					continue;
				}
				if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
					args->interpolation = OPENCV_CUBIC;
					args->clamp = TRUE;
					continue;
				}
				if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
					args->interpolation = OPENCV_LANCZOS4;
					args->clamp = TRUE;
					continue;
				}
				if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
					args->interpolation = OPENCV_LINEAR;
					args->clamp = FALSE;
					continue;
				}
				if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
					args->interpolation = OPENCV_AREA;
					args->clamp = FALSE;
					continue;
				}
				siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
				return CMD_ARG_ERROR;
			}
		}
	}
	if (args->scale < 0.1) {
		siril_log_message(_("Error: scale cannot be less than 0.01, aborting.\n"));
		free(args->prefix);
		free(args);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	if (fabs(args->scale - 1.0) < 1.e-10) {
		siril_log_message(_("Scale is 1.0, nothing to do.\n"));
		free(args->prefix);
		free(args);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}

	scale_sequence(args);
	return CMD_OK;
}

// ============================================================================
// BG (Background)
// ============================================================================
struct bg_data {
	void (*destructor)(void *);
	rectangle selection;
	gboolean has_selection;
	double bg_values[3];  // Store background for up to 3 channels
	WORD us_bg_values[3];
	int nchannels;
};

static int bg_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct bg_data *data = (struct bg_data *)args->user;

	image_cfa_warning_check();

	rectangle *sel_ptr = data->has_selection ? &data->selection : NULL;
	data->nchannels = fit->naxes[2];

	for (int layer = 0; layer < fit->naxes[2]; layer++) {
		double bg = background(fit, layer, sel_ptr, MULTI_THREADED);
		data->bg_values[layer] = bg;

		if (fit->type == DATA_USHORT) {
			data->us_bg_values[layer] = round_to_WORD(bg);
			data->bg_values[layer] = bg / get_normalized_value(fit);
		} else if (fit->type == DATA_FLOAT) {
			data->us_bg_values[layer] = float_to_ushort_range(bg);
		}
	}

	return 0;
}
/*
static void bg_idle_function(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	struct bg_data *data = (struct bg_data *)args->user;

	for (int layer = 0; layer < data->nchannels; layer++) {
		siril_log_message(_("Background value (channel: #%d): %d (%.3e)\n"),
						  layer, data->us_bg_values[layer], data->bg_values[layer]);
	}
}
*/
static gchar *bg_log_hook(gpointer p, log_hook_detail detail) {
	struct bg_data *args = (struct bg_data*) p;
	if (detail == DETAILED && args->nchannels > 0) {
		if (args->nchannels == 1) {
			return g_strdup_printf(_("Background: %d (%.3e)"),
								   args->us_bg_values[0], args->bg_values[0]);
		} else {
			return g_strdup_printf(_("Background (RGB): %d, %d, %d"),
								   args->us_bg_values[0],
								   args->us_bg_values[1],
								   args->us_bg_values[2]);
		}
	}
	return g_strdup(_("Background computation"));
}

int process_bg(int nb) {
	// Allocate and initialize user data
	struct bg_data *data = calloc(1, sizeof(struct bg_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->has_selection = (com.selection.w > 0 && com.selection.h > 0);
	if (data->has_selection) {
		data->selection = com.selection;
	}
	data->nchannels = 0;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = bg_image_hook;
	args->idle_function = NULL;
	args->description = _("Background");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;  // This doesn't modify gfit
	args->command = TRUE;
	args->user = data;
	args->log_hook = bg_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

static int bgnoise_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct noise_data *data = (struct noise_data *)args->user;

	image_cfa_warning_check();

	// Call the noise worker directly
	noise_worker(data);

	return 0;
}

static gchar *bgnoise_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup(_("Background noise estimation"));
}

int process_bgnoise(int nb) {
	// Allocate noise_data
	struct noise_data *noise_args = calloc(1, sizeof(struct noise_data));
	if (!noise_args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	noise_args->fit = gfit;
	noise_args->use_idle = FALSE;
	noise_args->display_results = TRUE;
	noise_args->display_start_end = TRUE;
	memset(noise_args->bgnoise, 0.0, sizeof(double[3]));

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(noise_args);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = bgnoise_image_hook;
	args->idle_function = NULL;  // noise_worker handles its own idle function
	args->description = _("Background noise");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;  // This doesn't modify gfit
	args->command = TRUE;
	args->user = noise_args;
	args->log_hook = bgnoise_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_histo(int nb) {
	GError *error = NULL;
	gchar *end;
	int nlayer = g_ascii_strtoull(word[1], &end, 10);
	const gchar* clayer;

	if (end == word[1] || nlayer > 3 || nlayer < 0)
		return CMD_INVALID_IMAGE;
	image_cfa_warning_check();
	gsl_histogram *histo = computeHisto(gfit, nlayer);
	if (!isrgb(gfit))
		clayer = "bw";		//if B&W
	else
		clayer = channel_number_to_name(nlayer);
	gchar *filename = g_strdup_printf("histo_%s.dat", clayer);

	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "Cannot save histo\n");
		}
		g_object_unref(file);
		return CMD_FILE_NOT_FOUND;
	}
	for (size_t i = 0; i < USHRT_MAX + 1; i++) {
		gchar *buffer = g_strdup_printf("%zu %d\n", i, (int) gsl_histogram_get (histo, i));

		if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_free(buffer);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			return CMD_GENERIC_ERROR;
		}
		g_free(buffer);
	}

	siril_log_message(_("The file %s has been created for the %s layer.\n"), g_file_peek_path(file), clayer);

	g_object_unref(output_stream);
	g_object_unref(file);
	gsl_histogram_free(histo);
	return CMD_OK;
}

int process_tilt(int nb) {
	if (word[1] && !g_ascii_strcasecmp(word[1], "clear")) {
		clear_sensor_tilt();
		siril_log_message(_("Clearing tilt information\n"));
		queue_redraw(REDRAW_OVERLAY);
	} else {
		set_cursor_waiting(TRUE);
		draw_sensor_tilt(gfit);
		set_cursor_waiting(FALSE);
	}

	return CMD_OK;
}

int process_inspector(int nb) {
	if (com.headless) {
		siril_log_color_message(_("Error: cannot run inspector while headless"), "red");
		return CMD_ARG_ERROR;
	}
	compute_aberration_inspector();
	return CMD_OK;
}

// Enum to distinguish threshold operations
typedef enum {
	THRESH_BOTH,
	THRESH_LO,
	THRESH_HI
} thresh_type;

// Structure to hold threshold-specific data
struct thresh_data {
	void (*destructor)(void *);  // Required as first member
	thresh_type type;
	int lo;
	int hi;
};

// Image processing hook for threshold operations
static int thresh_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct thresh_data *data = (struct thresh_data *)args->user;

	// Perform the threshold operation(s)
	if (data->type == THRESH_BOTH || data->type == THRESH_LO) {
		threshlo(fit, data->lo);
	}
	if (data->type == THRESH_BOTH || data->type == THRESH_HI) {
		threshhi(fit, data->hi);
	}

	// Add to history
	char log[90];
	if (data->type == THRESH_BOTH) {
		sprintf(log, "Image clamped to [%d, %d]", data->lo, data->hi);
	} else if (data->type == THRESH_LO) {
		sprintf(log, "Image clamped to [%d, max]", data->lo);
	} else {
		sprintf(log, "Image clamped to [min, %d]", data->hi);
	}
	fit->history = g_slist_append(fit->history, strdup(log));

	return 0;
}

static gchar *thresh_log_hook(gpointer p, log_hook_detail detail) {
	struct thresh_data *args = (struct thresh_data *) p;
	if (args->type == THRESH_LO) {
		return g_strdup_printf(_("Low threshold: %d"), args->lo);
	} else if (args->type == THRESH_HI) {
		return g_strdup_printf(_("High threshold: %d"), args->hi);
	}
	return g_strdup_printf(_("Threshold: lo = %d, hi = %d"), args->lo, args->hi);
}

// Unified threshold processing function
int process_thresh(int nb) {
	gboolean mask_aware = FALSE;
	const char *num_args[2] = { NULL, NULL };
	int num_index = 0;
	gchar *end1, *end2;

	int maxlevel = (gfit->orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	thresh_type type;
	int lo = 0, hi = maxlevel;

	/* Scan arguments (order-independent for -mask) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (num_index < 2) {
			num_args[num_index++] = word[i];
		} else {
			siril_log_message(_("Invalid argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	/* Determine operation type by command name */
	if (strcmp(word[0], "thresh") == 0) {
		if (num_index != 2)
			return CMD_ARG_ERROR;

		type = THRESH_BOTH;

		lo = g_ascii_strtoull(num_args[0], &end1, 10);
		if (end1 == num_args[0] || lo < 0 || lo > maxlevel) {
			siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
			return CMD_ARG_ERROR;
		}

		hi = g_ascii_strtoull(num_args[1], &end2, 10);
		if (end2 == num_args[1] || hi < 0 || hi > maxlevel) {
			siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
			return CMD_ARG_ERROR;
		}

		if (lo >= hi) {
			siril_log_message(_("lo must be strictly smaller than hi\n"));
			return CMD_ARG_ERROR;
		}

	} else if (strcmp(word[0], "threshlo") == 0) {
		if (num_index != 1)
			return CMD_ARG_ERROR;

		type = THRESH_LO;

		lo = g_ascii_strtoull(num_args[0], &end1, 10);
		if (end1 == num_args[0] || lo < 0 || lo > maxlevel) {
			siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
			return CMD_ARG_ERROR;
		}

	} else if (strcmp(word[0], "threshhi") == 0) {
		if (num_index != 1)
			return CMD_ARG_ERROR;

		type = THRESH_HI;

		hi = g_ascii_strtoull(num_args[0], &end1, 10);
		if (end1 == num_args[0] || hi < 0 || hi > maxlevel) {
			siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
			return CMD_ARG_ERROR;
		}

	} else {
		return CMD_ARG_ERROR;
	}

	/* Allocate and initialize user data */
	struct thresh_data *data = calloc(1, sizeof(struct thresh_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}


	data->type = type;
	data->lo = lo;
	data->hi = hi;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;  // Threshold operations work in-place
	args->image_hook = thresh_image_hook;
	args->idle_function = NULL;
	args->description = _("Threshold operation");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = thresh_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	/* Start the worker thread */
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}


// Structure for neg operation (simple, no extra data needed)
struct neg_data {
	void (*destructor)(void *);  // Required as first member
};

// Image processing hook for neg
static int neg_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	pos_to_neg(fit);
	return 0;
}

// Main command function for neg
int process_neg(int nb) {
	// Allocate and initialize user data
	struct neg_data *data = calloc(1, sizeof(struct neg_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	gboolean use_mask = FALSE;
	if (!g_strcmp0(word[1], "-mask")) {
		use_mask = TRUE;
	}


	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;  // Neg operation works in-place
	args->image_hook = neg_image_hook;
	args->idle_function = NULL;
	args->description = _("Negative");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->user = data;
	args->mask_aware = use_mask;
	// No log_hook required here as the basic description is enough
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	// Start the worker thread
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

struct nozero_data {
	void (*destructor)(void *);
	WORD level;
};

static int nozero_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct nozero_data *data = (struct nozero_data *)args->user;
	nozero(fit, data->level);
	return 0;
}

static gchar *nozero_log_hook(gpointer p, log_hook_detail detail) {
	struct nozero_data *args = (struct nozero_data*) p;
	return g_strdup_printf(_("Replace zeros with %d"), args->level);
}

int process_nozero(int nb) {
	gchar *end;
	int level = g_ascii_strtoull(word[1], &end, 10);
	int maxlevel = (gfit->orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	if (end == word[1] || level < 0 || level > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}

	// Allocate and initialize user data
	struct nozero_data *data = calloc(1, sizeof(struct nozero_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->level = (WORD)level;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = nozero_image_hook;
	args->idle_function = NULL;
	args->description = _("Replace zeros");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->log_hook = nozero_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

// Structure to hold DDP-specific data
struct ddp_data {
	void (*destructor)(void *);  // Required as first member, can be NULL
	float level;
	float coeff;
	float sigma;
};

// Image processing hook for DDP
static int ddp_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct ddp_data *data = (struct ddp_data *)args->user;

	// Perform the DDP operation
	int retval = ddp(fit, data->level, data->coeff, data->sigma);

	return retval;
}

static gchar *ddp_log_hook(gpointer p, log_hook_detail detail) {
	struct ddp_data *args = (struct ddp_data*) p;
	if (detail == SUMMARY) {
		return g_strdup_printf(_("DDP: level %.3f, coeff %.3f, Sigma %.3f"), args->level, args->coeff, args->sigma);
	}
	return g_strdup_printf(_("Digital Development Processing: level %.3f, coeff %.3f, Sigma %.3f"), args->level, args->coeff, args->sigma);
}

// Main command function
int process_ddp(int nb) {
	gboolean mask_aware = FALSE;
	const char *num_args[3] = { NULL, NULL, NULL };
	int num_index = 0;
	gchar *end;

	/* Scan arguments (order-independent for -mask) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (num_index < 3) {
			num_args[num_index++] = word[i];
		} else {
			siril_log_message(_("Invalid argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (num_index != 3) {
		siril_log_message(_("DDP requires 3 numeric arguments: level coeff sigma\n"));
		return CMD_ARG_ERROR;
	}

	float level = (float) g_ascii_strtod(num_args[0], &end);
	if (end == num_args[0] || level < 0 || level > USHRT_MAX_SINGLE) {
		siril_log_message(_("Level value is incorrect\n"));
		return CMD_ARG_ERROR;
	}

	float coeff = g_ascii_strtod(num_args[1], &end);
	if (end == num_args[1] || coeff < 0) {
		siril_log_message(_("Coeff value is incorrect\n"));
		return CMD_ARG_ERROR;
	}

	float sigma = g_ascii_strtod(num_args[2], &end);
	if (end == num_args[2] || sigma < 0) {
		siril_log_message(_("Sigma value is incorrect\n"));
		return CMD_ARG_ERROR;
	}

	image_cfa_warning_check();

	/* Allocate and initialize user data */
	struct ddp_data *ddp_args = calloc(1, sizeof(struct ddp_data));
	if (!ddp_args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	ddp_args->destructor = NULL;  // Use default free()
	ddp_args->level = level;
	ddp_args->coeff = coeff;
	ddp_args->sigma = sigma;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(ddp_args);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.f;
	args->image_hook = ddp_image_hook;
	args->idle_function = NULL;  // Use default
	args->description = _("Digital Development Processing");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;  // This command modifies gfit
	args->command = TRUE; // calling as command, not from GUI
	args->user = ddp_args;
	args->mask_aware = mask_aware;
	args->log_hook = ddp_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	/* Start the worker thread */
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}


int process_new(int nb){
	int width, height, layers;
	gchar *end, *endw, *endh;
	char *filename = NULL;

	width = g_ascii_strtod(word[1], &endw);
	height = g_ascii_strtod(word[2], &endh);
	layers = g_ascii_strtoull(word[3], &end, 10);
	if (end == word[3] || (layers != 1 && layers != 3)) {
		siril_log_message(_("Number of layers MUST be 1 or 3\n"));
		return CMD_ARG_ERROR;
	}
	if (endw == word[1] || endh == word[2] || !height || !width) {
		siril_log_message(_("Width and Height must be positive values.\n"));
		return CMD_ARG_ERROR;
	}

	/* If a filename has been set */
	if (nb == 5 && (word[4][0] != '\0')) {
		filename = strdup(word[4]);
	} else {
		filename = strdup(_("new empty image"));
	}

	close_single_image();
	close_sequence(FALSE);

	fits *fit = gfit;
	if (new_fit_image(&fit, width, height, layers, DATA_FLOAT))
		return CMD_GENERIC_ERROR;
	memset(gfit->fdata, 0, width * height * layers * sizeof(float));

	com.seq.current = UNRELATED_IMAGE;
	create_uniq_from_gfit(filename, FALSE);
	gui_function(open_single_image_from_gfit, NULL);
	if (!com.headless) {
		if (g_main_context_is_owner(g_main_context_default())) {
			// it is safe to call the function directly
			open_single_image_from_gfit(NULL);
		} else {
			// we aren't in the GTK main thread or a script, so we run the idle and wait for it
			execute_idle_and_wait_for_it(open_single_image_from_gfit, NULL);
		}
	}
	return CMD_OK;
}

int process_visu(int nb) {
	gchar *end1, *end2;
	int low = g_ascii_strtoull(word[1], &end1, 10);
	int high = g_ascii_strtoull(word[2], &end2, 10);
	if ((end2 == word[2] || high > USHRT_MAX) || (end1 == word[1] || low < 0)) {
		siril_log_message(_("Values must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}
	visu(gfit, low, high);
	return CMD_OK;
}

struct ffill_data {
	void (*destructor)(void *);
	int level;
	rectangle area;
};

static int ffill_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct ffill_data *data = (struct ffill_data *)args->user;
	rectangle area = data->area;

	int retval = fill(fit, data->level, &area);
	if (retval) {
		return retval;
	}

	// Mirror operation
	area.x = fit->rx - area.x - area.w;
	area.y = fit->ry - area.y - area.h;
	return fill(fit, data->level, &area);
}

static gchar *ffill_log_hook(gpointer p, log_hook_detail detail) {
	struct ffill_data *args = (struct ffill_data*) p;
	return g_strdup_printf(_("Fill mirrored region with level %d"), args->level);
}

int process_ffill(int nb) {
	if (nb < 2) {
		siril_log_message(_("Missing fill level argument.\n"));
		return CMD_ARG_ERROR;
	}

	gboolean mask_aware = FALSE;
	char *level_arg = NULL;
	rectangle area;
	gchar *end;

	/* Scan arguments to find level and optional -mask */
	for (int i = 1; i < nb; i++) {
		if (!g_strcmp0(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!level_arg) {
			level_arg = word[i];
		} else if (i >= 2 && nb >= 6) {
			// Assume rectangle coordinates if provided
			continue;
		} else {
			siril_log_message(_("Unknown argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!level_arg) {
		siril_log_message(_("Missing fill level argument.\n"));
		return CMD_ARG_ERROR;
	}

	/* Parse fill level */
	int level = g_ascii_strtoull(level_arg, &end, 10);
	if (end == level_arg || level < 0 || level > USHRT_MAX) {
		siril_log_message(_("Value must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}

	/* Determine area */
	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			gchar *endx, *endy, *endw, *endh;
			area.x = g_ascii_strtoull(word[2], &endx, 10);
			area.y = g_ascii_strtoull(word[3], &endy, 10);
			area.w = g_ascii_strtoull(word[4], &endw, 10);
			area.h = g_ascii_strtoull(word[5], &endh, 10);
			if (endx == word[2] || endy == word[3] || endw == word[4]
				|| endh == word[5] || (area.w + area.x > gfit->rx)
				|| (area.h + area.y > gfit->ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return CMD_ARG_ERROR;
			}
		} else {
			siril_log_message(_("Fill2: select a region or provide x, y, width, height\n"));
			return CMD_ARG_ERROR;
		}
	} else {
		area = com.selection;
	}

	/* Allocate and initialize user data */
	struct ffill_data *data = calloc(1, sizeof(struct ffill_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->level = level;
	data->area = area;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = ffill_image_hook;
	args->idle_function = NULL;
	args->description = _("Fill mirrored region");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = ffill_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

cmd_errors parse_findstar(struct starfinder_data *args, int start, int nb) {
	for (int i = start; i < nb; i++) {
		char *current = word[i], *value;
		if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			/* Make sure path exists */
			gchar *dirname = g_path_get_dirname(value);
			if (siril_mkdir_with_parents(dirname, 0755) < 0) {
				g_free(dirname);
				return CMD_GENERIC_ERROR;
			}
			g_free(dirname);
			args->starfile = g_strdup(value);
			siril_debug_print("Findstar: saving at %s\n", args->starfile);
		} else if (g_str_has_prefix(word[i], "-layer=")) {
			int nb_layers = (start == 2) ? args->im.from_seq->nb_layers : args->im.fit->naxes[2];
			if (nb_layers == 1) {  // handling mono case
				siril_log_message(_("This sequence is mono, ignoring layer number.\n"));
				continue;
			}
			value = current + 7;
			gchar *end;
			int layer = g_ascii_strtoull(value, &end, 10);
			if (end == value || layer < 0 || layer > 2) {
				siril_log_message(_("Unknown layer number %s, must be between 0 and 2, will use green layer.\n"), value);
				layer = GLAYER;
				continue;
			}
			args->layer = layer;
		} else if (g_str_has_prefix(word[i], "-maxstars=")) {
			char *current = word[i], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			gchar *end;
			int max_stars = g_ascii_strtoull(value, &end, 10);
			if (end == value || max_stars > MAX_STARS_FITTED || max_stars < MIN_STARS_FITTED) {
				// limiting values to avoid too long computation or too low number of candidates
				siril_log_message(_("Max number of stars %s not allowed. Should be between %d and %d.\n"), value, MIN_STARS_FITTED, MAX_STARS_FITTED);
				return CMD_ARG_ERROR;
			}
			args->max_stars_fitted = max_stars;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	}
	return CMD_OK;
}

int process_findstar(int nb) {
	int layer;
	if (!(com.script || (com.python_script && !com.headless))) {
		layer = select_vport(gui.cvport);
	} else {
		layer = (gfit->naxes[2] > 1) ? GLAYER : RLAYER;
	}

	struct starfinder_data *args = calloc(1, sizeof(struct starfinder_data));
	args->layer = layer;
	args->im.fit = gfit;
	if (sequence_is_loaded() && com.seq.current >= 0) {
		args->im.from_seq = &com.seq;
		args->im.index_in_seq = com.seq.current;
	} else {
		args->im.from_seq = NULL;
		args->im.index_in_seq = -1;
	}

	// initializing args
	args->starfile = NULL;
	args->max_stars_fitted = 0;
	args->threading = MULTI_THREADED;
	args->update_GUI = !com.script || com.python_script;
	siril_debug_print("findstar profiling %s stars\n", (com.pref.starfinder_conf.profile == PSF_GAUSSIAN) ? "Gaussian" : "Moffat");

	cmd_errors argparsing = parse_findstar(args, 1, nb);

	if (argparsing) {
		if (args->starfile) g_free(args->starfile);
		free(args);
		return argparsing;
	}
	if (!com.script && com.selection.w != 0 && com.selection.h != 0) {
		args->selection = com.selection;
	}
	if (has_wcs(args->im.fit) && (args->starfile || args->update_GUI)) {
		args->save_eqcoords = TRUE;
		args->ref_wcs = args->im.fit->keywords.wcslib;
	} else {
		args->save_eqcoords = FALSE;
	}

	if (!start_in_new_thread(findstar_worker, args)) {
		g_free(args->starfile);
		free(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_seq_findstar(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	struct starfinder_data *args = calloc(1, sizeof(struct starfinder_data));
	int layer;
	if (!(com.script || (com.python_script && !com.headless)) && check_seq_is_comseq(seq)) { // we use vport only if seq is com.seq
		layer = select_vport(gui.cvport);
	} else {
		layer = (seq->nb_layers > 1) ? GLAYER : RLAYER;
	}
	// initializing findstar args
	args->layer = layer;
	args->im.fit = NULL;
	args->im.from_seq = seq;
	args->im.index_in_seq = -1;
	args->max_stars_fitted = 0;
	args->update_GUI = FALSE;
	args->save_to_file = TRUE;
	args->save_eqcoords = TRUE;	// managed in apply_findstar_to_sequence()

	args->starfile = NULL;
	cmd_errors argparsing = parse_findstar(args, 2, nb);

	if (argparsing) {
		if (args->starfile) g_free(args->starfile);
		free(args);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return argparsing;
	}
	if (args->starfile) {
		siril_log_message(_("Option -out= is not available for sequences, ignoring\n"));
		g_free(args->starfile);
		args->starfile = NULL;
	}

	fits tmpfit = { 0 };
	seq_read_frame_metadata(seq, sequence_find_refimage(seq), &tmpfit);
	gboolean mono = (tmpfit.naxes[2] == 1);
	gboolean cfa = (tmpfit.keywords.bayer_pattern[0] != '\0');
	clearfits(&tmpfit);
	if (mono && cfa) {
		control_window_switch_to_tab(OUTPUT_LOGS);
		siril_log_color_message(_("Warning: sequence contains undebayered CFA images. Star detection may produce results for this sequence but will not perform optimally and star parameters may be inaccurate.\n"), "salmon");
	}

	return apply_findstar_to_sequence(args);
}

int process_findhot(int nb){
	if (gfit->naxes[2] != 1) {
		siril_log_message(_("find_hot must be applied on an one-channel master-dark frame"));
		return CMD_NEED_INIT_FIRST;
	}
	GError *error = NULL;
	long icold, ihot;
	gchar type, *end;
	double sig[2];

	sig[0] = g_ascii_strtod(word[2], &end);
	if (end == word[2]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[2]);
		return CMD_ARG_ERROR;
	}
	sig[1] = g_ascii_strtod(word[3], &end);
	if (end == word[3]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[3]);
		return CMD_ARG_ERROR;
	}

	deviant_pixel *dev = find_deviant_pixels(gfit, sig, &icold, &ihot, FALSE);
	siril_log_message(_("%ld cold and %ld hot pixels\n"), icold, ihot);

	gchar *filename = g_strdup_printf("%s.lst", word[1]);
	GFile *file = g_file_new_for_path(filename);

	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			siril_log_message(_("Could not open file: %s\n"), filename);
		}
		g_object_unref(file);
		g_free(filename);
		return CMD_FILE_NOT_FOUND;
	}
	g_free(filename);

	for (int i = 0; i < icold + ihot; i++) {
		int y = gfit->ry - (int) dev[i].p.y - 1;  /* FITS is stored bottom to top */
		if (dev[i].type == HOT_PIXEL)
			type = 'H';
		else
			type = 'C';
		gchar *buffer = g_strdup_printf("P %d %d %c\n", (int) dev[i].p.x, y, type);
		if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_free(buffer);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			return CMD_GENERIC_ERROR;
		}
		g_free(buffer);
	}

	free(dev);
	g_object_unref(output_stream);
	g_object_unref(file);

	return CMD_OK;
}

static int fix_xtrans_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	return fix_xtrans_ac(fit);
}

int process_fix_xtrans(int nb) {
	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = fix_xtrans_image_hook;
	args->description = _("Fix X-Trans artefacts");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_cosme(int nb) {
	gchar *filename;
	if (!g_str_has_suffix(word[1], ".lst")) {
		filename = g_strdup_printf("%s.lst", word[1]);
	} else {
		filename = g_strdup(word[1]);
	}

	if (!g_file_test(filename, G_FILE_TEST_EXISTS)) {
		siril_log_color_message(_("File [%s] does not exist.\n"), "red", filename);
		g_free(filename);
		return CMD_FILE_NOT_FOUND;
	}

	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	int is_cfa = (word[0][5] == '_') ? 1 : 0;

	// Allocate parameters using the new allocator
	struct cosme_data *params = new_cosme_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		g_object_unref(file);
		return CMD_GENERIC_ERROR;
	}

	params->file = file;  // Transfer ownership
	params->is_cfa = is_cfa;
	params->fit = NULL;  // Not used in this context
	params->seq = NULL;  // Not used in this context
	params->prefix = NULL;  // Not used in this context

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_cosme_data(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;  // Cosmetic correction works in-place with minimal overhead
	args->image_hook = cosme_image_hook_generic;
	args->idle_function = NULL;
	args->description = _("Cosmetic Correction");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_seq_cosme(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	gchar *filename;

	if (!g_str_has_suffix(word[2], ".lst")) {
		filename = g_strdup_printf("%s.lst", word[2]);
	} else {
		filename = g_strdup(word[2]);
	}

	if (!g_file_test(filename, G_FILE_TEST_EXISTS)) {
		siril_log_color_message(_("File [%s] does not exist.\n"), "red", filename);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		g_free(filename);
		return CMD_FILE_NOT_FOUND;
	}

	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	struct cosme_data *args = calloc(1, sizeof(struct cosme_data));
	args->prefix = strdup("cosme_");

	if (g_str_has_prefix(word[3], "-prefix=")) {
		char *current = word[3], *value;
		value = current + 8;
		if (value[0] == '\0') {
			g_object_unref(file);
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
			free(args->prefix);
			free(args);
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
		free(args->prefix);
		args->prefix = strdup(value);
	}

	args->seq = seq;
	args->is_cfa = (word[0][8] == '_') ? 1 : 0;
	args->file = file;
	args->fit = gfit;

	apply_cosme_to_sequence(args);

	return CMD_OK;
}

int process_fmedian(int nb) {
	gboolean mask_aware = FALSE;
	const char *num_args[2] = { NULL, NULL };
	int num_index = 0;
	gchar *end1, *end2;

	/* Scan arguments (order-independent for -mask) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (num_index < 2) {
			num_args[num_index++] = word[i];
		} else {
			siril_log_message(_("Invalid argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (num_index != 2) {
		siril_log_message(_("Median filter requires 2 arguments: ksize amount\n"));
		return CMD_ARG_ERROR;
	}

	/* Allocate and initialize parameters */
	struct median_filter_data *params = calloc(1, sizeof(struct median_filter_data));
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	params->ksize = g_ascii_strtoull(num_args[0], &end1, 10);
	if (end1 == num_args[0] || !(params->ksize & 1) ||
	    params->ksize < 3 || params->ksize > 15) {
		siril_log_message(_("The size of the kernel MUST be odd and in the range [3, 15].\n"));
		free(params);
		return CMD_ARG_ERROR;
	}

	params->amount = g_ascii_strtod(num_args[1], &end2);
	if (end2 == num_args[1] || params->amount < 0.0 || params->amount > 1.0) {
		siril_log_message(_("Modulation value MUST be between 0 and 1\n"));
		free(params);
		return CMD_ARG_ERROR;
	}

	params->iterations = 1;
	params->fit = gfit;

	/* Allocate worker args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = median_image_hook;
	args->idle_function = NULL;
	args->description = _("Median filter");
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->verbose = FALSE;
	args->user = params;
	args->mask_aware = mask_aware;
	args->log_hook = median_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	image_cfa_warning_check();

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

/* The name of this command should be COG in english but this choice
 * was done to be consistent with IRIS
 */
struct cdg_data {
	void (*destructor)(void *);
	float x_avg;
	float y_avg;
	gboolean success;
};

static int cdg_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct cdg_data *data = (struct cdg_data *)args->user;
	data->success = !FindCentre(fit, &data->x_avg, &data->y_avg);
	if (data->success) {
		siril_log_message(_("Center of gravity coordinates are (%.3lf, %.3lf)\n"),
						  data->x_avg, data->y_avg);
	}
	return 0;
}

static gchar *cdg_log_hook(gpointer p, log_hook_detail detail) {
	struct cdg_data *args = (struct cdg_data*) p;
	if (detail == DETAILED && args->success) {
		return g_strdup_printf(_("Center of gravity: (%.3f, %.3f)"),
							   args->x_avg, args->y_avg);
	}
	return g_strdup(_("Center of gravity computation"));
}

int process_cdg(int nb) {
	// Allocate and initialize user data
	struct cdg_data *data = calloc(1, sizeof(struct cdg_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->success = FALSE;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = cdg_image_hook;
	args->idle_function = NULL;
	args->description = _("Center of gravity");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;  // This doesn't modify gfit
	args->command = TRUE;
	args->user = data;
	args->log_hook = cdg_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

static gboolean clear_log_buffer(gpointer user_data) {
	GtkTextView *text = GTK_TEXT_VIEW(lookup_widget("output"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
	return FALSE;
}

int process_clear(int nb) {
	gui_function(clear_log_buffer, NULL);
	return CMD_OK;
}

int process_clearstar(int nb){
	execute_idle_and_wait_for_it(clear_stars_list_as_idle, GINT_TO_POINTER(TRUE));
	notify_gfit_modified();
	queue_redraw(REDRAW_OVERLAY);
	gui_function(redraw_previews, NULL);
	return CMD_OK;
}

int process_close(int nb) {
	close_sequence(FALSE);
	close_single_image();
	return CMD_OK;
}

struct fill_data {
	void (*destructor)(void *);
	int level;
	rectangle area;
};

static int fill_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct fill_data *data = (struct fill_data *)args->user;
	rectangle area = data->area;
	return fill(fit, data->level, &area);
}

static gchar *fill_log_hook(gpointer p, log_hook_detail detail) {
	struct fill_data *args = (struct fill_data*) p;
	return g_strdup_printf(_("Fill region with level %d"), args->level);
}

int process_fill(int nb) {
	rectangle area;
	gchar *end;
	gboolean mask_aware = FALSE;

	/* Determine area */
	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb >= 6) {
			gchar *endx, *endy, *endw, *endh;
			area.x = g_ascii_strtoull(word[2], &endx, 10);
			area.y = g_ascii_strtoull(word[3], &endy, 10);
			area.w = g_ascii_strtoull(word[4], &endw, 10);
			area.h = g_ascii_strtoull(word[5], &endh, 10);
			if (endx == word[2] || endy == word[3] || endw == word[4]
				|| endh == word[5] || (area.w + area.x > gfit->rx)
				|| (area.h + area.y > gfit->ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return CMD_ARG_ERROR;
			}
		} else {
			area.w = gfit->rx;
			area.h = gfit->ry;
			area.x = 0;
			area.y = 0;
		}
	} else {
		area = com.selection;
	}

	/* Scan arguments to find level and optional -mask */
	char *level_arg = NULL;
	for (int i = 1; i < nb; i++) {
		if (!g_strcmp0(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!level_arg) {
			level_arg = word[i];
		} else if (i >= 2 && nb >= 6) {
			// Skip rectangle arguments if present
			continue;
		} else {
			siril_log_message(_("Unknown argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!level_arg) {
		siril_log_message(_("Missing fill value argument.\n"));
		return CMD_ARG_ERROR;
	}

	/* Parse fill level */
	int level = g_ascii_strtoull(level_arg, &end, 10);
	if (end == level_arg || level < 0 || level > USHRT_MAX) {
		siril_log_message(_("Value must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}

	/* Allocate and initialize user data */
	struct fill_data *data = calloc(1, sizeof(struct fill_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->level = level;
	data->area = area;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = fill_image_hook;
	args->idle_function = NULL;
	args->description = _("Fill region");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = fill_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

struct offset_data {
	void (*destructor)(void *);
	float level;
};

static int offset_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct offset_data *data = (struct offset_data *)args->user;
	off(fit, data->level);
	return 0;
}

static gchar *offset_log_hook(gpointer p, log_hook_detail detail) {
	struct offset_data *args = (struct offset_data*) p;
	return g_strdup_printf(_("Offset: %.0f"), args->level);
}

int process_offset(int nb) {
	if (nb < 2) {
		siril_log_message(_("Missing offset level argument.\n"));
		return CMD_ARG_ERROR;
	}

	gboolean mask_aware = FALSE;
	char *level_arg = NULL;

	/* Scan arguments (order-independent) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!level_arg) {
			level_arg = word[i];
		} else {
			siril_log_message(_("Unknown argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	if (!level_arg) {
		siril_log_message(_("Offset level argument not specified.\n"));
		return CMD_ARG_ERROR;
	}

	/* Parse level */
	gchar *end;
	int level = g_ascii_strtoull(level_arg, &end, 10);
	if (end == level_arg) {
		siril_log_message(_("Invalid offset level argument.\n"));
		return CMD_ARG_ERROR;
	}

	/* Allocate and initialize user data */
	struct offset_data *data = calloc(1, sizeof(struct offset_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	data->level = (float)level;

	/* Allocate and initialize generic_img_args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = offset_image_hook;
	args->idle_function = NULL;
	args->description = _("Offset");
	args->verbose = TRUE;
	args->command_updates_gfit = TRUE;
	args->command = TRUE;
	args->user = data;
	args->mask_aware = mask_aware;
	args->log_hook = offset_log_hook;
	args->max_threads = 1;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

/* Command interpreter function for SCNR using generic_image_worker */
int process_scnr(int nb) {
	gboolean mask_aware = FALSE;
	gboolean preserve = TRUE;

	const char *num_args[2] = { NULL, NULL };
	int num_index = 0;

	/* Scan arguments (order-independent flags) */
	for (int i = 1; i < nb; i++) {
		if (!strcmp(word[i], "-mask")) {
			mask_aware = TRUE;
		} else if (!strcmp(word[i], "-nopreserve")) {
			preserve = FALSE;
		} else if (num_index < 2) {
			num_args[num_index++] = word[i];
		} else {
			siril_log_message(_("Invalid argument %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	scnr_type type = SCNR_AVERAGE_NEUTRAL;
	double amount = 0.0;

	/* Parse type (optional) */
	if (num_index >= 1) {
		gchar *end;
		type = g_ascii_strtoull(num_args[0], &end, 10);
		if (end == num_args[0] || type > 3) {
			siril_log_message(
				_("Type can either be 0 (average neutral), "
				  "1 (maximum neutral), 2 (maximum mask) or 3 (additive mask)\n"));
			return CMD_ARG_ERROR;
		}
	}

	/* Parse amount (required only for mask types) */
	if (type == SCNR_MAXIMUM_MASK || type == SCNR_ADDITIVE_MASK) {
		if (num_index < 2) {
			siril_log_message(_("Amount argument is required for mask modes\n"));
			return CMD_ARG_ERROR;
		}

		gchar *end;
		amount = g_ascii_strtod(num_args[1], &end);
		if (end == num_args[1] || amount < 0.0 || amount > 1.0) {
			siril_log_message(_("Amount can only be [0, 1]\n"));
			return CMD_ARG_ERROR;
		}
	} else if (num_index > 1) {
		/* Too many numeric arguments */
		siril_log_message(_("Too many arguments\n"));
		return CMD_ARG_ERROR;
	}

	/* Allocate parameters */
	struct scnr_data *params = new_scnr_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_GENERIC_ERROR;
	}

	params->type = type;
	params->amount = amount;
	params->preserve = preserve;
	params->verbose = TRUE;
	params->applying = TRUE;

	/* Allocate worker args */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_scnr_data(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 1.5f;
	args->image_hook = scnr_image_hook;
	args->idle_function = NULL;
	args->description = _("SCNR");
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = scnr_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->mask_aware = mask_aware;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_fft(int nb){
	struct fft_data *args = calloc(1, sizeof(struct fft_data));

	args->fit = gfit;
	args->type = g_strdup(word[0]);
	args->modulus = g_strdup(word[1]);
	args->phase = g_strdup(word[2]);
	args->type_order = 0;
	image_cfa_warning_check();
	if (!start_in_new_thread(fourier_transform, args)) {
		free(args->type);
		free(args->modulus);
		free(args->phase);
		free(args);
	}

	return CMD_OK;
}

int process_fixbanding(int nb) {
	gchar *end1, *end2;
	double amount = g_ascii_strtod(word[1], &end1);
	if (end1 == word[1] || amount < 0 || amount > 4) {
		siril_log_message(_("Amount value must be in the [0, 4] range.\n"));
		return CMD_ARG_ERROR;
	}
	double sigma = g_ascii_strtod(word[2], &end2);
	if (end2 == word[2] || sigma < 0 || sigma > 5) {
		siril_log_message(_("1/sigma value must be in the [0, 5] range.\n"));
		return CMD_ARG_ERROR;
	}

	gboolean applyRotation = FALSE;
	if (nb > 3) {
		int arg_index = 3;
		while (arg_index < nb && word[arg_index]) {
			char *arg = word[arg_index];
			if (!g_strcmp0(arg, "-vertical")) {
				applyRotation = TRUE;
			} else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
			arg_index++;
		}
	}

	image_cfa_warning_check();

	// Allocate parameters using the allocator
	struct banding_data *params = new_banding_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return CMD_GENERIC_ERROR;
	}

	params->protect_highlights = TRUE;
	params->amount = amount;
	params->sigma = sigma;
	params->applyRotation = applyRotation;
	params->seqEntry = NULL;
	params->seq = NULL;
	params->fit = NULL;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_banding_data(params);
		free(params);
		return CMD_GENERIC_ERROR;
	}

	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = banding_single_image_hook;
	args->idle_function = NULL; // Use default idle function for command-line
	args->description = _("Canon Banding Reduction");
	args->command_updates_gfit = TRUE;
	args->command = TRUE; // calling as command, not from GUI
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = banding_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_banding_data(params);
		free(params);
		free(args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_seq_fixbanding(int nb) {
	struct banding_data *args = calloc(1, sizeof(struct banding_data));
	gchar *end1 = NULL, *end2 = NULL;
	args->seq = NULL;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		free (args);
		return CMD_SEQUENCE_NOT_FOUND;
	}
	args->seq = seq;
	args->amount = g_ascii_strtod(word[2], &end1);
	if (end1 == word[2] || args->amount < 0 || args->amount > 4) {
		siril_log_message(_("Amount value must be in the [0, 4] range.\n"));
		if (!check_seq_is_comseq(args->seq))
			free_sequence(args->seq, TRUE);
		free(args);
		return CMD_ARG_ERROR;
	}
	args->sigma = g_ascii_strtod(word[3], &end2);
	if (end2 == word[3] || args->sigma < 0 || args->sigma > 5) {
		siril_log_message(_("1/sigma value must be in the [0, 5] range.\n"));
		if (!check_seq_is_comseq(args->seq))
			free_sequence(args->seq, TRUE);
		free(args);
		return CMD_ARG_ERROR;
	}
	// settings default optional values
	args->protect_highlights = TRUE;
	args->applyRotation = FALSE;
	args->fit = NULL;

	if (nb > 4) {
		int arg_index = 4;
		while (arg_index < nb && word[arg_index]) {
			char *arg = word[arg_index];
			if (g_str_has_prefix(arg, "-prefix=")) {
				char *value = arg + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), arg);
					if (!check_seq_is_comseq(args->seq))
						free_sequence(args->seq, TRUE);
					free((char*) args->seqEntry);
					free(args);
					return CMD_ARG_ERROR;
				}
				args->seqEntry = strdup(value);
			} else if (!g_strcmp0(arg, "-vertical")) {
				args->applyRotation = TRUE;
			} else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				free((char*) args->seqEntry);
				printf("seqEntry freed");
				if (!check_seq_is_comseq(args->seq))
					free_sequence(args->seq, TRUE);
				free(args);
				return CMD_ARG_ERROR;
			}
			arg_index++;
		}
	}
	if (!args->seqEntry)
		args->seqEntry = strdup("unband_");
	sequence_cfa_warning_check(seq);
	apply_banding_to_sequence(args);
	return CMD_OK;
}

int process_subsky(int nb) {
	sequence *seq = NULL;
	int degree = 0, samples = 20;
	double tolerance = 1.0, smooth = 0.5;
	gboolean dithering;
	background_interpolation interp;
	char *prefix = NULL;
	gboolean use_existing = FALSE;

	int arg_index = 1;
	gboolean is_sequence = (word[0][2] == 'q');

	if (is_sequence) {
		arg_index = 2;
		dithering = TRUE;
		seq = load_sequence(word[1], NULL);
		if (!seq) {
			return CMD_SEQUENCE_NOT_FOUND;
		}
		if (check_seq_is_comseq(seq)) {
			free_sequence(seq, TRUE);
			seq = &com.seq;
		}
	} else {
		if (!single_image_is_loaded()) return CMD_IMAGE_NOT_FOUND;
		dithering = FALSE;
	}

	if (!strcmp(word[arg_index], "-rbf"))
		interp = BACKGROUND_INTER_RBF;
	else {
		gchar *end;
		interp = BACKGROUND_INTER_POLY;
		degree = g_ascii_strtoull(word[arg_index], &end, 10);
		if (end == word[arg_index] || degree < 1 || degree > 4) {
			siril_log_message(_("Polynomial degree order must be within the [1, 4] range.\n"));
			return CMD_ARG_ERROR;
		}
	}

	arg_index++;
	while (arg_index < nb && word[arg_index]) {
		char *arg = word[arg_index];
		if (is_sequence && g_str_has_prefix(arg, "-prefix=")) {
			if (prefix) {
				siril_log_message(_("There can be only one prefix argument"));
				free(prefix);
				return CMD_ARG_ERROR;
			}
			char *value = arg + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
			prefix = strdup(value);
		}
		else if (g_str_has_prefix(arg, "-samples=")) {
			gchar *end;
			char *value = arg + 9;
			samples = g_ascii_strtoull(value, &end, 10);
			if (end == value || samples <= 1) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				free(prefix);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-tolerance=")) {
			char *next, *value = arg + 11;
			tolerance = g_ascii_strtod(value, &next);
			if (next == value || tolerance < 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				free(prefix);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-smooth=")) {
			char *next, *value = arg + 8;
			smooth = g_ascii_strtod(value, &next);
			if (next == value || smooth < 0.0 || smooth > 1.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				free(prefix);
				return CMD_ARG_ERROR;
			}
			if (interp != BACKGROUND_INTER_RBF)
				siril_log_color_message(_("smooth parameter is unused with the polynomial model, ignoring.\n"), "salmon");
		}
		else if (is_sequence && !g_strcmp0(arg, "-nodither")) {
			dithering = FALSE;
		}
		else if (!is_sequence && !g_strcmp0(arg, "-dither")) {
			dithering = TRUE;
		}
		else if (!is_sequence && !g_strcmp0(arg, "-existing")) {
			use_existing = TRUE;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			free(prefix);
			return CMD_ARG_ERROR;
		}
		arg_index++;
	}

	struct background_data *args = calloc(1, sizeof(struct background_data));
	args->nb_of_samples = samples;
	args->tolerance = tolerance;
	args->correction = BACKGROUND_CORRECTION_SUBTRACT;
	args->interpolation_method = interp;
	args->degree = (poly_order) (degree - 1);
	args->smoothing = smooth;
	args->threads = com.max_thread;
	args->dither = dithering;
	args->from_ui = FALSE;
	siril_debug_print("dithering: %s\n", dithering ? "enabled" : "disabled");

	if (is_sequence) {
		args->seq = seq;
		args->seqEntry = prefix ? prefix : strdup("bkg_");
		apply_background_extraction_to_sequence(args);
	} else {
		args->seq = NULL;
		args->seqEntry = NULL;
		args->fit = gfit;

		// Check if the image has a Bayer CFA pattern
		sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit->keywords.bayer_pattern);
		gboolean is_cfa = gfit->naxes[2] == 1 && pattern >= BAYER_FILTER_MIN && pattern <= BAYER_FILTER_MAX;

		int retval = 1;
		if (use_existing) {
			if (!com.grad_samples) {
				siril_log_color_message(_("Error, no existing samples available\n"), "red");
				free(args);
				return CMD_GENERIC_ERROR;
			} else {
				retval = 0;
			}
		} else {
			retval = generate_background_samples(samples, tolerance);
		}
		if (!retval) {
			start_in_new_thread(is_cfa ? remove_gradient_from_cfa_image :
								remove_gradient_from_image, args);
		} else {
			siril_log_color_message(_("Error generating background samples\n"), "red");
			free(args);
			return CMD_GENERIC_ERROR;
		}
	}

	return CMD_OK;
}

int process_findcosme(int nb) {
	sequence *seq = NULL;
	int i = 0;

	gboolean is_sequence = (word[0][0] == 's');

	if (is_sequence) {
		seq = load_sequence(word[1], NULL);
		if (!seq) {
			return CMD_SEQUENCE_NOT_FOUND;
		}
		i++;
	}

	struct cosmetic_data *args = calloc(1, sizeof(struct cosmetic_data));
	gchar *end1, *end2;

	args->seq = seq;
	args->sigma[0] = g_ascii_strtod(word[1 + i], &end1);
	if (end1 == word[1 + i] || args->sigma[0] < 0) {
		siril_log_message(_("Sigma low must be positive.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}
	args->sigma[1] = g_ascii_strtod(word[2 + i], &end2);
	if (end2 == word[2 + i] || args->sigma[1] < 0) {
		siril_log_message(_("Sigma high must be positive.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}

	args->is_cfa = (word[0][10] == '_' || word[0][13] == '_');	// find_cosme_cfa or seqfind_cosme_cfa
	args->amount = 1.0;
	args->fit = gfit;

	if (is_sequence) {
		args->seqEntry = strdup("cc_");
		args->threading = SINGLE_THREADED;

		int startoptargs = i + 3;
		int nb_command_max = i + 4;
		if (nb > startoptargs) {
			for (int j = startoptargs; j < nb_command_max; j++) {
				if (word[j]) {
					if (g_str_has_prefix(word[j], "-prefix=")) {
						char *current = word[j], *value;
						value = current + 8;
						if (value[0] == '\0') {
							siril_log_message(_("Missing argument to %s, aborting.\n"), current);
							return CMD_ARG_ERROR;
						}
						args->seqEntry = strdup(value);
					}
				}
			}
		}
		apply_cosmetic_to_sequence(args);
	} else {
		args->threading = MULTI_THREADED;
		if (!start_in_new_thread(autoDetectThreaded, args)) {
			free(args->seqEntry);
			free(args);
			return CMD_GENERIC_ERROR;
		}
	}

	return CMD_OK;
}

static gboolean select_update_gui(gpointer user_data) {
	update_stack_interface(TRUE);
	update_reg_interface(FALSE);
	adjust_sellabel();
	drawPlot();
	return FALSE;
}

int select_unselect(gboolean select) {
	char *end1, *end2;
	int from = g_ascii_strtoull(word[2], &end1, 10);
	int to = g_ascii_strtoull(word[3], &end2, 10);
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}
	if (end1 == word[2] || from < 1 || from > seq->number) {
		siril_log_message(_("The second argument must be between 1 and the number of images.\n"));
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}
	if (end2 == word[3] || to < from) {
		siril_log_message(_("The third argument must be larger or equal than the \"from\" argument.\n"));
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}
	if (to > seq->number) {
		siril_log_color_message(_("The third argument is larger than the number of images.\n"), "salmon");
		siril_log_message(_("Re-adjusting to %d.\n"), "salmon", seq->number);
		to = seq->number;
	}
	for (int i = from - 1; i <= to - 1; i++) { // use real index
		if (i >= seq->number) break;
		if (seq->imgparam[i].incl != select) {
			seq->imgparam[i].incl = select;
			if (select)
				seq->selnum++;
			else
				seq->selnum--;
		}
		if (!select && seq->reference_image == i) {
			seq->reference_image = -1;
		}
	}
	writeseqfile(seq);
	if (check_seq_is_comseq(seq)) {
		fix_selnum(&com.seq, FALSE);
		if (!com.headless)
			execute_idle_and_wait_for_it(select_update_gui, NULL);
	}
	siril_log_message(_("Selection update finished, %d images are selected in the sequence\n"), seq->selnum);

	if (!check_seq_is_comseq(seq))
		free_sequence(seq, TRUE);

	return CMD_OK;
}

int process_select(int nb){
	return select_unselect(TRUE);
}

int process_unselect(int nb){
	return select_unselect(FALSE);
}

int process_split(int nb){
	struct extract_channels_data *args = calloc(1, sizeof(struct extract_channels_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->channel[0] = g_strdup_printf("%s%s", word[1], com.pref.ext);
	args->channel[1] = g_strdup_printf("%s%s", word[2], com.pref.ext);
	args->channel[2] = g_strdup_printf("%s%s", word[3], com.pref.ext);

	args->fit = calloc(1, sizeof(fits));
	if (copyfits(gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Could not copy the input image, aborting.\n"));
		clearfits(args->fit);
		free(args->fit);
		free(args->channel[0]);
		free(args->channel[1]);
		free(args->channel[2]);
		free(args);
		return CMD_ALLOC_ERROR;
	}

	if (nb > 4) {
		if (!g_ascii_strcasecmp(word[4], "-hsl")) {
			args->type = EXTRACT_HSL;
			args->str_type = _("HSL");
		} else if (!g_ascii_strcasecmp(word[4], "-hsv")) {
			args->type = EXTRACT_HSV;
			args->str_type = _("HSV");
		} else if (!g_ascii_strcasecmp(word[4], "-lab")) {
			args->type = EXTRACT_CIELAB;
			args->str_type = _("CieLAB");
		} else {
			args->type = EXTRACT_RGB;
			args->str_type = _("RGB");
		}
	} else {
		args->type = EXTRACT_RGB;
		args->str_type = _("RGB");
	}

	copy_fits_metadata(gfit, args->fit);

	args->fit->keywords.bayer_pattern[0] = '\0'; // Mark this as no longer having a Bayer pattern

	if (!start_in_new_thread(extract_channels, args)) {
		clearfits(args->fit);
		free(args->fit);
		free(args->channel[0]);
		free(args->channel[1]);
		free(args->channel[2]);
		free(args);
	}
	return CMD_OK;
}

int process_split_cfa(int nb) {
	char *filename = NULL;
	int ret = 1;

	fits f_cfa0 = { 0 }, f_cfa1 = { 0 }, f_cfa2 = { 0 }, f_cfa3 = { 0 };

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}

	gchar *cfa0 = g_strdup_printf("CFA0_%s%s", filename, com.pref.ext);
	gchar *cfa1 = g_strdup_printf("CFA1_%s%s", filename, com.pref.ext);
	gchar *cfa2 = g_strdup_printf("CFA2_%s%s", filename, com.pref.ext);
	gchar *cfa3 = g_strdup_printf("CFA3_%s%s", filename, com.pref.ext);

	if (gfit->type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, 0) ||
				save1fits16(cfa1, &f_cfa1, 0) ||
				save1fits16(cfa2, &f_cfa2, 0) ||
				save1fits16(cfa3, &f_cfa3, 0);
		}
	}
	else if (gfit->type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits32(cfa0, &f_cfa0, 0) ||
				save1fits32(cfa1, &f_cfa1, 0) ||
				save1fits32(cfa2, &f_cfa2, 0) ||
				save1fits32(cfa3, &f_cfa3, 0);
		}
	}

	g_free(cfa0); g_free(cfa1);
	g_free(cfa2); g_free(cfa3);
	clearfits(&f_cfa0); clearfits(&f_cfa1);
	clearfits(&f_cfa2); clearfits(&f_cfa3);
	free(filename);
	return ret;
}

int process_extractGreen(int nb) {
	char *filename = NULL;
	int ret = 1;

	fits f_green = { 0 };

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}

	sensor_pattern pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
		siril_log_color_message(_("This image does not have a Bayer CFA pattern, cannot extract green channel.\n"), "red");
		g_free(filename);
		return CMD_INVALID_IMAGE;
	}

	gchar *green = g_strdup_printf("Green_%s%s", filename, com.pref.ext);
	if (gfit->type == DATA_USHORT) {
		if (!(ret = extractGreen_ushort(gfit, &f_green, pattern))) {
			ret = save1fits16(green, &f_green, 0);
		}
	}
	else if (gfit->type == DATA_FLOAT) {
		if (!(ret = extractGreen_float(gfit, &f_green, pattern))) {
			ret = save1fits32(green, &f_green, 0);
		}
	} else {
		g_free(filename);
		g_free(green);
		clearfits(&f_green);
		return CMD_INVALID_IMAGE;
	}

	g_free(green);
	clearfits(&f_green);
	free(filename);
	return ret;

}

int extract_Ha(extraction_scaling scaling) {
	int ret = 1;
	char *filename = NULL;
	fits f_Ha = { 0 };
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}
	sensor_pattern pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
		siril_log_color_message(_("This image does not have a Bayer CFA pattern, cannot extract Ha.\n"), "red");
		g_free(filename);
		return CMD_INVALID_IMAGE;
	}
	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	if (gfit->type == DATA_USHORT) {
		if (!(ret = extractHa_ushort(gfit, &f_Ha, pattern, scaling))) {
			ret = save1fits16(Ha, &f_Ha, 0);
		}
	}
	else if (gfit->type == DATA_FLOAT) {
		if (!(ret = extractHa_float(gfit, &f_Ha, pattern, scaling))) {
			ret = save1fits32(Ha, &f_Ha, 0);
		}
	} else ret = CMD_INVALID_IMAGE;
	g_free(Ha);
	clearfits(&f_Ha);
	free(filename);
	return ret;
}

int process_extractHa(int nb) {
	int ret = 1;
	extraction_scaling scaling = SCALING_NONE;
	if (g_str_has_prefix(word[1], "-upscale")) {
		scaling = SCALING_HA_UP;
		siril_log_message(_("Upscaling x2\n"));
	}
	ret = extract_Ha(scaling);
	return ret;
}

int extract_HaOIII(extraction_scaling scaling) {
	gchar *filename = NULL;
	int ret = 1;
	fits f_Ha = { 0 }, f_OIII = { 0 };
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}
	sensor_pattern pattern = get_validated_cfa_pattern(gfit,FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
		siril_log_color_message(_("This image does not have a Bayer CFA pattern, cannot extract Ha/OIII channels.\n"), "red");
		g_free(filename);
		return CMD_INVALID_IMAGE;
	}
	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	gchar *OIII = g_strdup_printf("OIII_%s%s", filename, com.pref.ext);
	if (gfit->type == DATA_USHORT) {
		if (!(ret = extractHaOIII_ushort(gfit, &f_Ha, &f_OIII, pattern, scaling, com.max_thread))) {
			ret = save1fits16(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	}
	else if (gfit->type == DATA_FLOAT) {
		if (!(ret = extractHaOIII_float(gfit, &f_Ha, &f_OIII, pattern, scaling, com.max_thread))) {
			ret = save1fits32(Ha, &f_Ha, 0) ||
					save1fits32(OIII, &f_OIII, 0);
		}
	} else {
		g_free(Ha);
		g_free(OIII);
		g_free(filename);
		return CMD_INVALID_IMAGE;
	}
	g_free(Ha);
	g_free(OIII);
	clearfits(&f_Ha);
	clearfits(&f_OIII);
	g_free(filename);
	return ret;
}

int process_extractHaOIII(int nb) {
	extraction_scaling scaling = SCALING_NONE;
	if (word[1]) {
		if (g_str_has_prefix(word[1], "-resample=")) {
			char *current = word[1], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				return CMD_ARG_ERROR;
			} else if (!strcasecmp(value, "ha")) {
				scaling = SCALING_HA_UP;
			} else if (!strcasecmp(value, "oiii")) {
				scaling = SCALING_OIII_DOWN;
			}
		}
	}
	return extract_HaOIII(scaling);
}

int process_seq_mtf(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	struct mtf_data *args = calloc(1, sizeof(struct mtf_data));

	args->seq = seq;
	args->fit = gfit;
	args->seqEntry = strdup("mtf_");
	gchar *end1, *end2, *end3;
	args->params.shadows = g_ascii_strtod(word[2], &end1);
	args->params.midtones = g_ascii_strtod(word[3], &end2);
	args->params.highlights = g_ascii_strtod(word[4], &end3);
	args->params.do_red = TRUE;
	args->params.do_green = TRUE;
	args->params.do_blue = TRUE;
	if (end1 == word[2] || end2 == word[3] || end3 == word[4] ||
			args->params.shadows < 0.0 || args->params.midtones <= 0.0 || args->params.highlights <= 0.0 ||
			args->params.shadows >= 1.0 || args->params.midtones >= 1.0 || args->params.highlights > 1.0) {
		siril_log_message(_("Invalid argument to %s, aborting.\n"), word[0]);
		free(args->seqEntry);
		free(args);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	int startoptargs = 5;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (!strcmp(word[i], "R")) {
					args->params.do_green = FALSE;
					args->params.do_blue = FALSE;
				}
				if (!strcmp(word[i], "G")) {
					args->params.do_red = FALSE;
					args->params.do_blue = FALSE;
				}
				if (!strcmp(word[i], "B")) {
					args->params.do_green = FALSE;
					args->params.do_red = FALSE;
				}
				if (!strcmp(word[i], "RG")) {
					args->params.do_blue = FALSE;
				}
				if (!strcmp(word[i], "RB")) {
					args->params.do_green = FALSE;
				}
				if (!strcmp(word[i], "GB")) {
					args->params.do_red = FALSE;
				}
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), current);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
			}
		}
	}
	sequence_cfa_warning_check(seq);
	apply_mtf_to_sequence(args);

	return CMD_OK;
}

int process_seq_split_cfa(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
		return CMD_FOR_CFA_IMAGE;
	}

	struct multi_output_data *args = calloc(1, sizeof(struct multi_output_data));

	args->seq = seq;
	args->seqEntry = ("CFA"); // propose to default to "CFA" for consistency of output names with single image split_cfa
	args->n = 4;
	args->prefixes = calloc(5, sizeof(const char*));

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						free(args->seqEntry);
						free(args->prefixes);
						free(args);
						return CMD_ARG_ERROR;
					}
					free(args->seqEntry);
					args->seqEntry = strdup(value);
				}
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				free(args->seqEntry);
				free(args->prefixes);
				free(args);
				return CMD_ARG_ERROR;
			}
		}
	}
	if (args->seqEntry && args->seqEntry[0] == '\0') {
		free(args->seqEntry);
		args->seqEntry = strdup("CFA");
	}
	for (int i = 0 ; i < 4 ; i++) {
		args->prefixes[i] = g_strdup_printf("%s%d_", args->seqEntry, i);
	}

	apply_split_cfa_to_sequence(args);

	return CMD_OK;
}

int process_seq_merge_cfa(int nb) {
	sequence *seq0 = load_sequence(word[1], NULL);
	if (!seq0)
		return CMD_SEQUENCE_NOT_FOUND;
	if (seq0->nb_layers > 1) {
		if (!check_seq_is_comseq(seq0))
			free_sequence(seq0, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	sequence *seq1 = load_sequence(word[2], NULL);
	if (!seq1)
		return CMD_SEQUENCE_NOT_FOUND;
	if (seq1->nb_layers > 1) {
		if (!check_seq_is_comseq(seq0)) {
			free_sequence(seq0, TRUE);
		}
		free_sequence(seq1, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	sequence *seq2 = load_sequence(word[3], NULL);
	if (!seq2)
		return CMD_SEQUENCE_NOT_FOUND;
	if (seq2->nb_layers > 1) {
		if (!check_seq_is_comseq(seq0)) {
			free_sequence(seq0, TRUE);
		}
		free_sequence(seq1, TRUE);
		free_sequence(seq2, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	sequence *seq3 = load_sequence(word[4], NULL);
	if (!seq3)
		return CMD_SEQUENCE_NOT_FOUND;

	if (seq3->nb_layers > 1) {
		if (!check_seq_is_comseq(seq0)) {
			free_sequence(seq0, TRUE);
		}
		free_sequence(seq1, TRUE);
		free_sequence(seq2, TRUE);
		free_sequence(seq3, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}
	if (seq3->nb_layers > 1) {
		if (!check_seq_is_comseq(seq3))
			free_sequence(seq3, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	if ((seq0->rx != seq1->rx || seq0->rx != seq2->rx || seq0->rx != seq3->rx) ||
	(seq0->ry != seq1->ry || seq0->ry != seq2->ry || seq0->ry != seq3->ry) ||
	(seq0->nb_layers != seq1->nb_layers || seq0->nb_layers != seq2->nb_layers || seq0->nb_layers != seq3->nb_layers) ||
	(seq0->bitpix != seq1->bitpix || seq0->bitpix != seq2->bitpix || seq0->bitpix != seq3->bitpix) ||
	(seq0->number != seq1->number || seq0->number != seq2->number || seq0->number != seq3->number)) {
		siril_log_color_message(_("Error: sequences don't match (dimensions, bitdepth, number of images must all be the same)\n"), "red");
		if (!check_seq_is_comseq(seq0))
			free_sequence(seq0, TRUE);
		free_sequence(seq0, TRUE);
		free_sequence(seq1, TRUE);
		free_sequence(seq2, TRUE);
		return CMD_GENERIC_ERROR;
	}

	struct merge_cfa_data *args = calloc(1, sizeof(struct merge_cfa_data));
	if(!args) {
		if (!check_seq_is_comseq(seq0))
			free_sequence(seq0, TRUE);
		free_sequence(seq0, TRUE);
		free_sequence(seq1, TRUE);
		free_sequence(seq2, TRUE);
		return CMD_ALLOC_ERROR;
	}
	args->seq0 = seq0;
	args->seq1 = seq1;
	args->seq2 = seq2;
	args->seq3 = seq3;

	if (!strcmp(word[5], "RGGB")) {
		args->pattern = BAYER_FILTER_RGGB;
	} else if (!strcmp(word[5], "BGGR")) {
		args->pattern = BAYER_FILTER_BGGR;
	} else if (!strcmp(word[5], "GBRG")) {
		args->pattern = BAYER_FILTER_GBRG;
	} else if (!strcmp(word[5], "GRBG")) {
		args->pattern = BAYER_FILTER_GRBG;
	} else {
		siril_log_color_message(_("Invalid Bayer matrix specified!\n"), "red");
		if (!check_seq_is_comseq(seq0))
			free_sequence(seq0, TRUE);
		free_sequence(seq1, TRUE);
		free_sequence(seq2, TRUE);
		free_sequence(seq3, TRUE);
		free(args);
		return CMD_ARG_ERROR;
	}
	siril_log_message(_("Reconstructing %s Bayer matrix.\n"), word[2]);
	args->seqEntryOut = strdup("mCFA_");

	int startoptargs = 6;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-outprefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						if (!check_seq_is_comseq(seq0))
							free_sequence(seq0, TRUE);
						free_sequence(seq1, TRUE);
						free_sequence(seq2, TRUE);
						free_sequence(seq3, TRUE);
						free(args->seqEntryOut);
						free(args);
						return CMD_ARG_ERROR;
					}
					if (args->seqEntryOut)
						free(args->seqEntryOut);
					args->seqEntryOut = strdup(value);
				}
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
				if (!check_seq_is_comseq(seq0))
					free_sequence(seq0, TRUE);
				free_sequence(seq1, TRUE);
				free_sequence(seq2, TRUE);
				free_sequence(seq3, TRUE);
				free(args->seqEntryOut);
				free(args);
				return CMD_ARG_ERROR;
			}
		}
	}

	apply_mergecfa_to_sequence(args);

	return CMD_OK;
}

int process_seq_extractHa(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
		free_sequence(seq, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	struct simple_extract_data *args = calloc(1, sizeof(struct simple_extract_data));

	args->seq = seq;
	args->seqEntry = strdup("Ha_");

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						free(args->seqEntry);
						free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					free(args->seqEntry);
					args->seqEntry = strdup(value);
				}
				else if (g_str_has_prefix(word[i], "-upscale")) {
					args->scaling = SCALING_HA_UP;
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
					free(args->seqEntry);
					free(args);
					return CMD_ARG_ERROR;
				}
			}
		}
	}

	apply_extractHa_to_sequence(args);

	return CMD_OK;
}

int process_seq_extractGreen(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
		free_sequence(seq, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	struct simple_extract_data *args = calloc(1, sizeof(struct simple_extract_data));

	args->seq = seq;
	args->seqEntry = strdup("Green_");

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						free(args->seqEntry);
						free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					if (args->seqEntry) {
						free(args->seqEntry);
						args->seqEntry = strdup(value);
					}
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					free(args->seqEntry);
					free(args);
					return CMD_ARG_ERROR;
				}
			}
		}
	}

	apply_extractGreen_to_sequence(args);

	return CMD_OK;
}

int process_seq_extractHaOIII(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
		free_sequence(seq, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	struct multi_output_data *args = calloc(1, sizeof(struct multi_output_data));
	args->user_data = calloc(1, sizeof(extraction_scaling));
	*(extraction_scaling *) args->user_data = (extraction_scaling) SCALING_NONE;
	args->seq = seq;
	args->seqEntry = strdup("CFA"); // propose to default to "CFA" for consistency of output names with single image split_cfa
	args->n = 2;
	if (word[2]) {
		if (g_str_has_prefix(word[2], "-resample=")) {
			char *current = word[2], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[2]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				free(args->seqEntry);
				free(args->user_data);
				free(args);
				return CMD_ARG_ERROR;
			} else if (!strcmp(value, "ha")) {
				*(extraction_scaling*) args->user_data = SCALING_HA_UP;
			} else if (!strcmp(value, "oiii")) {
				*(extraction_scaling*) args->user_data = SCALING_OIII_DOWN;
			}
		}
	}
	args->prefixes = calloc(3, sizeof(char*));
	args->prefixes[0] = g_strdup("Ha_");
	args->prefixes[1] = g_strdup("OIII_");

	apply_extractHaOIII_to_sequence(args);

	return CMD_OK;
}

static int stat_cmd_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct stat_data *data = (struct stat_data *)args->user;
	rectangle *sel_ptr = data->has_selection ? &data->selection : NULL;

	for (int layer = 0; layer < data->nplane; layer++) {
		int super_layer = layer;
		if (data->cfa)
			super_layer = -layer - 1;

		data->stats[layer] = statistics(NULL, -1, fit, super_layer, sel_ptr, STATS_MAIN, MULTI_THREADED);
		if (!data->stats[layer]) {
			siril_log_message(_("Statistics computation failed for channel %d (all nil?).\n"), layer);
		}
	}

	return 0;
}

static gchar *stat_log_hook(gpointer p, log_hook_detail detail) {
	struct stat_data *data = (struct stat_data *)p;

	if (detail == SUMMARY) {
		if (data->cfa)
			return g_strdup(_("Statistics (CFA)"));
		else if (data->option == STATS_MAIN)
			return g_strdup(_("Statistics (main)"));
		else
			return g_strdup(_("Statistics"));
	}

	// DETAILED - recompute stats for logging
	int nplane = data->cfa ? 3 : data->fit->naxes[2];
	rectangle *sel_ptr = (data->selection.w > 0 && data->selection.h > 0) ? &data->selection : NULL;

	GString *result = g_string_new("");

	for (int layer = 0; layer < nplane; layer++) {
		int super_layer = layer;
		if (data->cfa)
			super_layer = -layer - 1;

		imstats *stat = statistics(NULL, -1, data->fit, super_layer, sel_ptr, STATS_MAIN, MULTI_THREADED);
		if (!stat)
			continue;

		const char *layername;
		switch (layer) {
			case 0:
				layername = (nplane == 1) ? "B&W" : "Red";
				break;
			case 1:
				layername = "Green";
				break;
			case 2:
				layername = "Blue";
				break;
			default:
				layername = "Unknown";
		}

		if (layer > 0)
			g_string_append(result, "\n");

		if (data->option == STATS_BASIC) {
			if (data->fit->type == DATA_USHORT) {
				g_string_append_printf(result,
					_("%s layer: Mean: %0.6f, Median: %0.1f, Sigma: %0.6f, "
					  "Min: %0.1f, Max: %0.1f, bgnoise: %0.6f"),
					layername, stat->mean, stat->median, stat->sigma,
					stat->min, stat->max, stat->bgnoise);
			} else {
				g_string_append_printf(result,
					_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
					  "Min: %0.1f, Max: %0.1f, bgnoise: %0.1f"),
					layername, stat->mean * USHRT_MAX_DOUBLE,
					stat->median * USHRT_MAX_DOUBLE,
					stat->sigma * USHRT_MAX_DOUBLE,
					stat->min * USHRT_MAX_DOUBLE,
					stat->max * USHRT_MAX_DOUBLE,
					stat->bgnoise * USHRT_MAX_DOUBLE);
			}
		} else if (data->option == STATS_MAIN) {
			if (data->fit->type == DATA_USHORT) {
				g_string_append_printf(result,
					_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
					  "Min: %0.1f, Max: %0.1f, bgnoise: %0.1f, "
					  "avgDev: %0.1f, MAD: %0.1f, sqrt(BWMV): %0.1f"),
					layername, stat->mean, stat->median, stat->sigma,
					stat->min, stat->max, stat->bgnoise, stat->avgDev,
					stat->mad, stat->sqrtbwmv);
			} else {
				g_string_append_printf(result,
					_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
					  "Min: %0.1f, Max: %0.1f, bgnoise: %0.1f, "
					  "avgDev: %0.1f, MAD: %0.1f, sqrt(BWMV): %0.1f"),
					layername, stat->mean * USHRT_MAX_DOUBLE,
					stat->median * USHRT_MAX_DOUBLE,
					stat->sigma * USHRT_MAX_DOUBLE,
					stat->min * USHRT_MAX_DOUBLE,
					stat->max * USHRT_MAX_DOUBLE,
					stat->bgnoise * USHRT_MAX_DOUBLE,
					stat->avgDev * USHRT_MAX_DOUBLE,
					stat->mad * USHRT_MAX_DOUBLE,
					stat->sqrtbwmv * USHRT_MAX_DOUBLE);
			}
		}

		free_stats(stat);
	}

	return g_string_free(result, FALSE);
}

int process_stat(int nb) {
	int nplane = gfit->naxes[2];
	gboolean cfa = FALSE;
	int option = STATS_BASIC;
	int argidx = 1;

	if (nb == 2 && !g_strcmp0(word[1], "-cfa") && nplane == 1 && gfit->keywords.bayer_pattern[0] != '\0') {
		siril_debug_print("Running stats on CFA\n");
		nplane = 3;
		cfa = TRUE;
		if ((com.selection.w && com.selection.w < 2) || (com.selection.h && com.selection.h < 2)) {
			siril_log_color_message(_("Statistics cannot be made on CFA images with a selection smaller than a 2x2 square, aborting\n"), "red");
			return CMD_GENERIC_ERROR;
		}
		argidx++;
	}

	if (nb == argidx + 1 && !g_strcmp0(word[argidx], "main"))
		option = STATS_MAIN;

	// Allocate and initialize user data
	struct stat_data *data = alloc_stat_data();
	if (!data) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	data->destroy_fn = free_stat_data;
	data->has_selection = (com.selection.w > 0 && com.selection.h > 0);
	if (data->has_selection) {
		data->selection = com.selection;
	}
	data->fit = gfit;
	data->option = option;
	data->cfa = cfa;
	data->nplane = nplane;
	memset(data->stats, 0, sizeof(data->stats));

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_stat_data(data);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = stat_cmd_image_hook;
	args->description = _("Statistics");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->log_hook = stat_log_hook;
	args->max_threads = 1;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_seq_stat(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	struct stat_data *args = calloc(1, sizeof(struct stat_data));
	args->seq = seq;
	args->csv_name = g_strdup(word[2]);
	args->selection = com.selection;
	args->option = STATS_MAIN;
	args->cfa = FALSE;

	if (nb > 3) {
		if (!g_strcmp0(word[3], "main")) {
			args->option = STATS_MAIN;
		} else if (!g_strcmp0(word[3], "full")) {
			args->option = STATS_NORM | STATS_MAIN; // adding STATS_MAIN to include also AVGDEV and SQRTBWMV
		} else if (!g_strcmp0(word[3], "basic")) {
			args->option = STATS_BASIC;
		} else if (!g_strcmp0(word[3], "-cfa")) {
			args->cfa = TRUE;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[3]);
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
			g_free(args->csv_name);
			free(args);
			return CMD_ARG_ERROR;
		}
		if (nb > 4) {
			if (!g_strcmp0(word[4], "-cfa")) {
				args->cfa = TRUE;
			} else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[4]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				g_free(args->csv_name);
				free(args);
				return CMD_ARG_ERROR;
			}
		}
	}

	apply_stats_to_sequence(args);

	return CMD_OK;
}

// Only for FITS images
int process_jsonmetadata(int nb) {
	char *input_filename = word[1];
	gchar *output_filename = NULL;
	gboolean use_gfit = FALSE, compute_stats = TRUE;

	// Process command arguments
	for (int i = 2; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-out=") && word[i][5] != '\0') {
			if (output_filename) g_free(output_filename);
			output_filename = g_strdup(word[i] + 5);
		} else if (!strcmp(word[i], "-stats_from_loaded")) {
			use_gfit = TRUE;
			if (!gfit->rx || !gfit->ry) {
				siril_log_color_message(_("No image appears to be loaded, reloading from '%s'\n"), "salmon", input_filename);
				use_gfit = FALSE;
			}
		} else if (!strcmp(word[i], "-nostats"))
			compute_stats = FALSE;
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			g_free(output_filename);
			return CMD_ARG_ERROR;
		}
	}

	if (!output_filename)
		output_filename = replace_ext(input_filename, ".json");

	int status = 0;
	fitsfile *fptr;
	if (siril_fits_open_diskfile_img(&fptr, input_filename, READONLY, &status)) {
		report_fits_error(status);
		g_free(output_filename);
		return CMD_GENERIC_ERROR;
	}

	GSList *header = read_header_keyvals_strings(fptr);

	imstats *stats[3] = { NULL };
	int nb_channels = 0;
	if (compute_stats) {
		if (use_gfit) {
			compute_all_channels_statistics_single_image(gfit, STATS_BASIC | STATS_FOR_CFA, MULTI_THREADED, stats);
			gboolean cfa = gfit->keywords.bayer_pattern[0] != '\0';
			nb_channels = cfa ? 3 : (int)gfit->naxes[2];
		} else {
			fits fit = { 0 };
			fit.fptr = fptr;
			if (read_fits_metadata(&fit) || read_fits_with_convert(&fit, input_filename, FALSE)) {
				fits_close_file(fptr, &status);
				return CMD_GENERIC_ERROR;
			}
			compute_all_channels_statistics_single_image(&fit, STATS_BASIC | STATS_FOR_CFA, MULTI_THREADED, stats);
			gboolean cfa = fit.keywords.bayer_pattern[0] != '\0';
			nb_channels = cfa ? 3 : (int)fit.naxes[2];
			clearfits(&fit);
		}
	}
	fits_close_file(fptr, &status);

	// Create the root object
	yyjson_mut_doc *doc = yyjson_mut_doc_new(NULL);
	yyjson_mut_val *root = yyjson_mut_obj(doc);
	yyjson_mut_doc_set_root(doc, root);

	// Add headers array
	yyjson_mut_val *headers_arr = yyjson_mut_arr(doc);
	yyjson_mut_obj_add_val(doc, root, "headers", headers_arr);

	// Add header records
	GSList *ptr = header;
	while (ptr) {
		header_record *r = ptr->data;
		yyjson_mut_val *header_obj = yyjson_mut_obj(doc);
		yyjson_mut_val *key_str = yyjson_mut_str(doc, r->key);
		yyjson_mut_val *value_str = yyjson_mut_str(doc, r->value);
		yyjson_mut_obj_add_val(doc, header_obj, "key", key_str);
		yyjson_mut_obj_add_val(doc, header_obj, "value", value_str);
		yyjson_mut_arr_append(headers_arr, header_obj);
		ptr = ptr->next;
	}

	// Add statistics if computed
	if (compute_stats) {
		yyjson_mut_val *stats_obj = yyjson_mut_obj(doc);
		yyjson_mut_obj_add_val(doc, root, "statistics", stats_obj);

		for (int i = 0; i < nb_channels; ++i) {
			if (stats[i]) {
				char channame[20];
				sprintf(channame, "channel%d", i);

				yyjson_mut_val *chan_obj = yyjson_mut_obj(doc);
				yyjson_mut_obj_add_val(doc, stats_obj, channame, chan_obj);

				yyjson_mut_val *mean_val = yyjson_mut_real(doc, stats[i]->mean);
				yyjson_mut_val *median_val = yyjson_mut_real(doc, stats[i]->median);
				yyjson_mut_val *sigma_val = yyjson_mut_real(doc, stats[i]->sigma);
				yyjson_mut_val *noise_val = yyjson_mut_real(doc, stats[i]->bgnoise);
				yyjson_mut_val *min_val = yyjson_mut_real(doc, stats[i]->min);
				yyjson_mut_val *max_val = yyjson_mut_real(doc, stats[i]->max);
				yyjson_mut_val *total_val = yyjson_mut_real(doc, stats[i]->total);
				yyjson_mut_val *good_val = yyjson_mut_real(doc, stats[i]->ngoodpix);

				yyjson_mut_obj_add_val(doc, chan_obj, "mean", mean_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "median", median_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "sigma", sigma_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "noise", noise_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "min", min_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "max", max_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "total_pix_count", total_val);
				yyjson_mut_obj_add_val(doc, chan_obj, "good_pix_count", good_val);

				free_stats(stats[i]);
			}
		}
	}

	// Write the JSON to file with pretty formatting
	yyjson_write_flag flg = YYJSON_WRITE_PRETTY;
	yyjson_write_err err;
	int retval = CMD_OK;

	if (!yyjson_mut_write_file(output_filename, doc, flg, NULL, &err)) {
		siril_log_message(_("Failed to save the JSON file %s: %s\n"),
						  output_filename, err.msg);
		retval = CMD_GENERIC_ERROR;
	} else {
		siril_log_message(_("Save metadata to the JSON file '%s'\n"), output_filename);
	}

	#ifdef DEBUG_TEST
	size_t len;
	char *json_str = yyjson_mut_write(doc, flg, &len);
	if (json_str) {
		printf("JSON:\n%s\n", json_str);
		free(json_str);
	}
	#endif

	yyjson_mut_doc_free(doc);
	g_free(output_filename);
	return retval;
}

int header_hook(struct generic_seq_metadata_args *args, fitsfile *fptr, int index) {
	GString *string = g_string_new(NULL);
	int status = 0;
	GSList *list = args->keys;

	for (int i = 0; i < g_slist_length(args->keys); i++) {
		char str[FLEN_VALUE] = { 0 };

		gchar *key = (gchar *)list->data;
		fits_read_keyword(fptr, key, str, NULL, &status);
		if (status) {
			strcpy(str, "N/A");
			status = 0;
		} else {
			if (i > 0) string = g_string_append(string, ",");
			string = g_string_append(string, str);
		}
		list = list->next;
	}

	if (args->output_stream) {
		GError *error = NULL;
		if (!g_output_stream_printf(args->output_stream, NULL, NULL, &error, "%d,%s\n", index + 1, string->str)) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			g_string_free(string, TRUE);
			return 1;
		}
	}
	else {
		gchar **token_keys = g_strsplit (args->header, ",", -1);
		gchar **token_values = g_strsplit (string->str, ",", -1);
		g_string_free(string, TRUE);
		gchar *output = NULL;
		for (int i = 0; i < g_strv_length(token_keys) && token_keys[i] && token_values[i]; i++) {
			gchar *tmp = g_strdup_printf("%s = %s,", token_keys[i], token_values[i]);
			if (output == NULL) {
				output = g_strdup(tmp);
			} else {
				gchar *old_output = output;
				output = g_strconcat(old_output, tmp, NULL);
				g_free(old_output);
			}
			g_free(tmp);
		}
		if (output) {
			siril_log_message(_("Image %d, %s\n"), index + 1, output);
			g_free(output);
		}
		g_strfreev(token_keys);
		g_strfreev(token_values);
	}
	return status;
}

int process_seq_header(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	gboolean filter = FALSE;
	GSList *list = NULL;
	int key;
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (seq->type != SEQ_REGULAR && seq->type != SEQ_FITSEQ) {
		siril_log_message(_("This command can only run for FITS images\n"));
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}

	GString *list_of_keys = g_string_new(NULL);

	for (key = 2; key < nb; key++) {
		if (!word[key] || word[key][0] == '-') {
			break;
		}
		list = g_slist_prepend(list, g_strdup(word[key]));
		if (key > 2) list_of_keys = g_string_append(list_of_keys, ", ");
		list_of_keys = g_string_append(list_of_keys, word[key]);
	}

	list = g_slist_reverse(list);
	gchar *header = g_string_free(list_of_keys, FALSE);

	if (!header) {
		return CMD_ARG_ERROR;
	}

	GOutputStream* output_stream = NULL;
	for (int i = key; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-out=")) {
			const char *arg = word[i] + 5;
			if (arg[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
				return CMD_ARG_ERROR;
			}
			GFile *file = g_file_new_for_path(arg);
			GError *error = NULL;
			output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
			g_object_unref(file);
			if (!output_stream) {
				if (error) {
					siril_log_color_message(_("Failed to create output file: %s\n"), "red", error->message);
					g_clear_error(&error);
				}
				return CMD_ARG_ERROR;
			}
			if (!g_output_stream_printf(output_stream, NULL, NULL, NULL, "# image number,%s\n", header)) {
				siril_log_color_message(_("Failed to write to output file\n"), "red");
				return CMD_ARG_ERROR;
			}
			siril_log_message(_("The file %s has been created.\n"), arg);
		} else if (!g_strcmp0(word[i], "-sel")) {
			filter = TRUE;
		}
	}

	struct generic_seq_metadata_args *args = calloc(1, sizeof(struct generic_seq_metadata_args));
	args->seq = seq;
	args->keys = list;
	args->image_hook = header_hook;
	args->output_stream = output_stream;
	args->header = header;
	args->filtering_criterion = filter ? seq_filter_included : seq_filter_all;
	if (!start_in_new_thread(generic_sequence_metadata_worker, args)) {
		g_slist_free_full(args->keys, g_free);
		g_free(header);
		g_object_unref(args->output_stream);
		free(args);
	}
	return 0;
}

struct file_time {
	gchar *file;
	GDateTime *time;
};
static gint sort_date(gconstpointer a, gconstpointer b) {
	return g_date_time_compare(((struct file_time*) a)->time,
			((struct file_time*) b)->time);
}
static gpointer extract(gconstpointer src, gpointer data) {
	g_date_time_unref(((struct file_time*) src)->time);
	return ((struct file_time*) src)->file;
}

int process_link(int nb) {
	if (word[1][0] == '-') {
		siril_log_message(_("First argument is the converted sequence name and shall not start with a -\n"));
		return CMD_ARG_ERROR;
	}
	if (file_name_has_invalid_chars(word[1])) {
		siril_log_color_message(_("Specified output name %s contains forbidden characters, aborting\n"), "red", word[1]);
		return CMD_ARG_ERROR;
	}
	char *destroot = strdup(word[1]);
	int idx = 1;
	gboolean sort_dateobs = FALSE;

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			gchar *end;
			idx = g_ascii_strtoull(value, &end, 10);
			if (end == value || idx <= 0 || idx >= INDEX_MAX) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), current);
				free(destroot);
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				free(destroot);
				return CMD_ARG_ERROR;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (siril_mkdir_with_parents(value, 0755) < 0) {
					free(destroot);
					return CMD_GENERIC_ERROR;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			free(destroot);
			destroot = strdup(filename);
			g_free(filename);
		}
		else if (!g_strcmp0(current, "-date")) {
			sort_dateobs = TRUE;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
			free(destroot);
			return CMD_ARG_ERROR;
		}
	}

	GDir *dir;
	GError *error = NULL;
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Link: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Link: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		free(destroot);
		return CMD_GENERIC_ERROR;
	}

	int count = 0;
	const gchar *file;
	GList *list = NULL;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext || file[0] == '.')
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPEFITS) {
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No FITS files were found for link\n"));
		free(destroot);
		return CMD_GENERIC_ERROR;
	}
	/* sort list */
	if (sort_dateobs) {
		GList *cur = list;
		GList *timed = NULL;
		gboolean failed = FALSE;
		siril_log_message("Sorting FITS files using DATE-OBS\n");
		while (cur) {
			gchar *filename = (gchar*) cur->data;
			GDateTime *date = get_date_from_fits(filename);
			if (!date) {
				siril_log_color_message(_("Sorting by date could not be done, file %s doesn't have DATE-OBS, not sorting\n"), "red", filename);
				failed = TRUE;
				break;
			}
			struct file_time *ft = calloc(1, sizeof(struct file_time));
			ft->file = filename;
			ft->time = date;
			timed = g_list_append(timed, ft);
			cur = cur->next;
		}
		if (!failed) {
			timed = g_list_sort(timed, (GCompareFunc) sort_date);
			list = g_list_copy_deep(timed, extract, NULL);
			g_list_free_full(timed, free);
		}
	} else {
		list = g_list_sort(list, (GCompareFunc) strcompare);
	}
	/* convert the list to an array for parallel processing */
	gchar **files_to_link = glist_to_array(list, &count);

	int nb_allowed;
	if (!allow_to_open_files(count, &nb_allowed)) {
		siril_log_message(_("Too many files: you should use the 'convert' command and pass the argument -fitseq to convert your sequence to fitseq format.\n"));
		g_strfreev(files_to_link);
		free(destroot);
		return CMD_GENERIC_ERROR;
	}

	gchar *str = ngettext("Link: processing %d FITS file...\n", "Link: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);

	if (!com.wd) {
		siril_log_message(_("Link: no working directory set.\n"));
		g_strfreev(files_to_link);
		free(destroot);
		return CMD_GENERIC_ERROR;
	}

	struct _convert_data *args = calloc(1, sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_link;
	args->total = count;
	args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = FALSE;
	args->multiple_output = FALSE;
	args->output_type = SEQ_REGULAR; // fallback if symlink does not work
	args->make_link = TRUE;
	gettimeofday(&(args->t_start), NULL);
	if (!start_in_new_thread(convert_thread_worker, args)) {
		free(args->destroot);
		g_strfreev(args->list);
		free(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_convert(int nb) {
	gboolean raw_only = word[0][7] != '\0';
	if (!com.wd) {
		siril_log_message(_("Conversion: no working directory set.\n"));
		return CMD_NO_CWD;
	}

	if (word[1][0] == '-') {
		siril_log_message(_("First argument is the converted sequence name and shall not start with a -\n"));
		return CMD_ARG_ERROR;
	}
	if (file_name_has_invalid_chars(word[1])) {
		siril_log_color_message(_("Specified output name %s contains forbidden characters, aborting\n"), "red", word[1]);
		return CMD_ARG_ERROR;
	}
	gchar *destroot = g_strdup(word[1]);
	int idx = 1;
	gboolean debayer = FALSE;
	gboolean make_link = !raw_only;
	sequence_type output = SEQ_REGULAR;
	gboolean in_cwd = TRUE;

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (!strcmp(current, "-debayer")) {
			debayer = TRUE;
			make_link = FALSE;
		} else if (!strcmp(current, "-fitseq")) {
			output = SEQ_FITSEQ;
			if (!g_str_has_suffix(destroot, com.pref.ext))
				str_append(&destroot, com.pref.ext);
		} else if (!strcmp(current, "-ser")) {
			output = SEQ_SER;
			if (!g_str_has_suffix(destroot, ".ser"))
				str_append(&destroot, ".ser");
		} else if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			gchar *end;
			idx = g_ascii_strtoull(value, &end, 10);
			if (end == value || idx <= 0 || idx >= INDEX_MAX) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), current);
				free(destroot);
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				free(destroot);
				return CMD_ARG_ERROR;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (siril_mkdir_with_parents(value, 0755) < 0) {
					free(destroot);
					return CMD_GENERIC_ERROR;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = strdup(filename);
			g_free(filename);
			in_cwd = FALSE;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
			free(destroot);
			return CMD_ARG_ERROR;
		}
	}

	GDir *dir;
	GError *error = NULL;
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		free(destroot);
		return CMD_NO_CWD;
	}

	if (in_cwd && (output == SEQ_SER || output == SEQ_FITSEQ) && g_file_test(destroot, G_FILE_TEST_EXISTS)) {
		siril_log_color_message(_("Destination sequence %s already exists in the current folder, cannot proceed.\n"), "red", destroot);
		free(destroot);
		g_dir_close(dir);
		return CMD_GENERIC_ERROR;
	}

	int count = 0;
	const gchar *file;
	GList *list = NULL;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext || file[0] == '.')
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPERAW && output == SEQ_SER && !g_ascii_strcasecmp(ext, "raf") && !debayer) {
			siril_log_message(_("FujiFilm XTRANS sensors are not supported by SER v2 (CFA-style) standard. You may use FITS sequences instead."));
			g_list_free_full(list, g_free);
			free(destroot);
			return CMD_GENERIC_ERROR;
		}
		if ((raw_only && type == TYPERAW) ||
				(!raw_only && type != TYPEUNDEF && type != TYPEAVI && type != TYPESER)) {
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		if (raw_only)
			siril_log_message(_("No RAW files were found for conversion\n"));
		else siril_log_message(_("No files were found for conversion\n"));
		free(destroot);
		return CMD_GENERIC_ERROR;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	gchar **files_to_convert = glist_to_array(list, &count);

	int nb_allowed;
	if (!allow_to_open_files(count, &nb_allowed) && output == SEQ_REGULAR) {
		siril_log_message(_("You should pass an extra argument -fitseq to convert your sequence to fitseq format.\n"));
		g_strfreev(files_to_convert);
		free(destroot);
		return CMD_GENERIC_ERROR;
	}

	gchar *str = ngettext("Convert: processing %d file...\n", "Convert: processing %d files...\n", count);
	siril_log_color_message(str, "green", count);

	struct _convert_data *args = calloc(1, sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_convert;
	args->total = count;
	args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = debayer;
	args->output_type = output;
	args->multiple_output = FALSE;
	args->make_link = make_link;
	gettimeofday(&(args->t_start), NULL);
	if (!start_in_new_thread(convert_thread_worker, args)) {
		free(args->destroot);
		g_strfreev(args->list);
		free(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_register(int nb) {
	cmd_errors errval = CMD_ARG_ERROR;
	struct registration_args *regargs = NULL;
	struct registration_method *method = NULL;
	char *msg = NULL;
	gboolean drizzle = FALSE;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	regargs = calloc(1, sizeof(struct registration_args));

	/* filling the arguments for registration */
	regargs->seq = seq;
	regargs->reference_image = sequence_find_refimage(seq);
	regargs->follow_star = FALSE;
	regargs->matchSelection = FALSE;
	regargs->no_output = FALSE;
	regargs->output_scale = 1.f;
	regargs->prefix = strdup("r_");
	regargs->min_pairs = 10; // 10 is good enough to ensure good matching
	regargs->max_stars_candidates = MAX_STARS_FITTED;
	regargs->type = HOMOGRAPHY_TRANSFORMATION;
	regargs->layer = (regargs->seq->nb_layers == 3) ? 1 : 0;
	regargs->interpolation = OPENCV_LANCZOS4;
	regargs->clamp = TRUE;
	regargs->undistort = DISTO_UNDEF;

	struct driz_args_t *driz = calloc(1, sizeof(struct driz_args_t));
	// Default values for the driz_args_t
	driz->use_flats = FALSE;
	driz->scale = 1.f;
	driz->kernel = kernel_square;
	driz->weight_scale = 1.f;
	driz->pixel_fraction = 1.f;

	/* check for options */
	for (int i = 2; i < nb; i++) {
		if (!strcmp(word[i], "-2pass")) {
			regargs->two_pass = TRUE;
			regargs->no_output = TRUE;
		} else if (!strcmp(word[i], "-noclamp")) {
			regargs->clamp = FALSE;
		} else if (!strcmp(word[i], "-nostarlist")) {
			regargs->no_starlist = TRUE;
		} else if (!strcmp(word[i], "-selected")) {
			regargs->filters.filter_included = TRUE;
		} else if (g_str_has_prefix(word[i], "-transf=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "shift", 5)) {
#ifdef HAVE_CV44
				regargs->type = SHIFT_TRANSFORMATION;
				continue;
#else
				siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
				goto terminate_register_on_error;
#endif
			}
			if(!g_ascii_strncasecmp(value, "similarity", 10)) {
				regargs->type = SIMILARITY_TRANSFORMATION;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "affine", 6)) {
				regargs->type = AFFINE_TRANSFORMATION;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "homography", 10)) {
				regargs->type = HOMOGRAPHY_TRANSFORMATION;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-layer=")) {
			if (regargs->seq->nb_layers == 1) {  // handling mono case
				siril_log_message(_("This sequence is mono, ignoring layer number.\n"));
				continue;
			}
			char *current = word[i], *value;
			value = current + 7;
			gchar *end;
			int layer = g_ascii_strtoull(value, &end, 10);
			if (end == value || layer < 0 || layer > 2) {
				siril_log_message(_("Unknown layer number %s, must be between 0 and 2, will use green layer.\n"), value);
				if (end == value) break;
				else continue;
			}
			regargs->layer = layer;
		} else if (g_str_has_prefix(word[i], "-prefix=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if (regargs->prefix)
				free(regargs->prefix);
			regargs->prefix = strdup(value);
		} else if (g_str_has_prefix(word[i], "-minpairs=")) {
			char *current = word[i], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			gchar *end;
			int min_pairs = g_ascii_strtoull(value, &end, 10);
			if (end == value || min_pairs < 4) { // using absolute min_pairs required by homography
				gchar *str = g_strdup_printf(_("%d smaller than minimum allowable star pairs: %d, aborting.\n"), min_pairs, regargs->min_pairs);
				siril_log_message(str);
				g_free(str);
				goto terminate_register_on_error;
			}
			regargs->min_pairs = min_pairs;
		} else if (g_str_has_prefix(word[i], "-maxstars=")) {
			char *current = word[i], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			gchar *end;
			int max_stars = g_ascii_strtoull(value, &end, 10);
			if (end == value || max_stars > MAX_STARS_FITTED || max_stars < MIN_STARS_FITTED) {
				// limiting values to avoid too long computation or too low number of candidates
				siril_log_message(_("Max number of stars %s not allowed. Should be between %d and %d.\n"), value, MIN_STARS_FITTED, MAX_STARS_FITTED);
				goto terminate_register_on_error;
			}
			regargs->max_stars_candidates = max_stars;
		} else if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				regargs->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				regargs->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				regargs->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				regargs->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "none", 4) || !g_ascii_strncasecmp(value, "no", 2)) {
				regargs->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				regargs->interpolation = OPENCV_AREA;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-disto=")) {
			char *current = word[i], *value;
			value = current + 7;
			gchar *filename = NULL;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if (!g_ascii_strncasecmp(value, "image", 5)) {
				regargs->undistort = DISTO_IMAGE;
			} else if (!g_ascii_strncasecmp(value, "file", 4)) {
				regargs->undistort = DISTO_FILE;
				if (i + 1 >= nb) {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					goto terminate_register_on_error;
				}
				filename = word[i + 1];
			} else if (!g_ascii_strncasecmp(value, "master", 6)) {
				regargs->undistort = DISTO_MASTER;
			} else {
				siril_log_message(_("Unknown distortion type %s, aborting.\n"), value);
				goto terminate_register_on_error;
			}
			// we now check the distortion params are ok
			fits reffit = { 0 };
			fits *preffit = &reffit;
			if (!check_seq_is_comseq(seq)) { // processing an image from the current sequence
				int image_to_load = sequence_find_refimage(seq);
				if (seq_read_frame_metadata(seq, image_to_load, preffit)) {
					siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
					goto terminate_register_on_error;
				}
			} else
				preffit = gfit;
			gchar *msgdisto = NULL;
			gboolean disto_valid = validate_disto_params(preffit, filename, regargs->undistort, &msgdisto, NULL);
			if (preffit != gfit)
				clearfits(preffit);
			if (!disto_valid) {
				siril_log_color_message("%s\n", "red", msgdisto);
				g_free(msgdisto);
				goto terminate_register_on_error;
			} else if(regargs->undistort == DISTO_FILE) {
				i++; // we have consumed one more argument for the filename
			}
			regargs->distoparam.index = regargs->undistort;
			if (filename) {
				regargs->distoparam.filename = g_strdup(filename);
			}
		} else if (!strcmp(word[i], "-drizzle")) {
			if (regargs->seq->nb_layers != 1) {  // handling mono case
				siril_log_message(_("This sequence is not mono / CFA, cannot drizzle.\n"));
				goto terminate_register_on_error;
			}
			drizzle = TRUE;
		// Drizzle options
		} else if (g_str_has_prefix(word[i], "-scale=")) {
			char *arg = word[i] + 7;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg || value < 0.1 || value > 3.) {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
			regargs->output_scale = (float)value;
			driz->scale = (float)value;
		} else if (g_str_has_prefix(word[i], "-pixfrac=")) {
			char *arg = word[i] + 9;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
			driz->pixel_fraction = (float) value;
		} else if (g_str_has_prefix(word[i], "-kernel=")) {
			char *arg = word[i] + 8;
			if (!g_ascii_strncasecmp(arg, "point", 5))
				driz->kernel = kernel_point;
			else if (!g_ascii_strncasecmp(arg, "turbo", 5))
				driz->kernel = kernel_turbo;
			else if (!g_ascii_strncasecmp(arg, "square", 6))
				driz->kernel = kernel_square;
			else if (!g_ascii_strncasecmp(arg, "gaussian", 8))
				driz->kernel = kernel_gaussian;
			else if (!g_ascii_strncasecmp(arg, "lanczos2", 8))
				driz->kernel = kernel_lanczos2;
			else if (!g_ascii_strncasecmp(arg, "lanczos3", 8))
				driz->kernel = kernel_lanczos3;
			else {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
		} else if (g_str_has_prefix(word[i], "-flat=")) {
			if (driz->flat) {
				siril_log_color_message(_("Error: flat image already set. Aborting.\n"), "red");
				goto terminate_register_on_error;
			}
			if (seq->is_variable) {
				siril_log_color_message(_("Error: flat image cannot work with variable sized sequence.\n"), "red");
				goto terminate_register_on_error;
			}
			char *flat_filename = word[i] + 6;
			fits reffit = { 0 };
			gchar *error = NULL;
			int status;
			if (seq_read_frame_metadata(seq, regargs->reference_image, &reffit)) {
				siril_log_color_message(_("NOT USING FLAT: Could not load reference image\n"), "red");
				clearfits(&reffit);
				goto terminate_register_on_error;
			}
			gchar *expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
			clearfits(&reffit);
			if (status) {
				error = _("NOT USING FLAT: could not parse the expression");
				// no need to free expression as we do not call path_parse with the NOFAIL mode
				goto terminate_register_on_error;
			} else {
				if (expression[0] == '\0') {
					siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
					g_free(expression);
					goto terminate_register_on_error;
				} else {
					driz->flat = calloc(1, sizeof(fits));
					if (!readfits(expression, driz->flat, NULL, TRUE)) {
						if (driz->flat->naxes[2] != seq->nb_layers) {
							error = _("NOT USING FLAT: number of channels is different");
						} else if (driz->flat->naxes[0] != seq->rx ||
								driz->flat->naxes[1] != seq->ry) {
							error = _("NOT USING FLAT: image dimensions are different");
						} else {
							// no need to deal with bitdepth conversion as readfits has already forced conversion to float
							siril_log_message(_("Master flat read for use as initial pixel weight\n"));
						}

					} else error = _("NOT USING FLAT: cannot open the file");
					g_free(expression);
					if (error) {
						goto terminate_register_on_error;
					}
				}
			}
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			goto terminate_register_on_error;
		}
	}

	if (drizzle) {
		if (seq->is_variable) {
			siril_log_color_message(_("Cannot drizzle sequences with images of different sizes, aborting.\n"), "red");
			goto terminate_register_on_error;
		}
		regargs->driz = driz;
		// we now check the drizzle params are ok
		fits reffit = { 0 };
		fits *preffit = &reffit;
		if (!check_seq_is_comseq(seq)) { // processing an image from the current sequence
			int image_to_load = sequence_find_refimage(seq);
			if (seq_read_frame_metadata(seq, image_to_load, preffit)) {
				siril_log_color_message(_("Could not load the reference image of the sequence, aborting.\n"), "red");
				goto terminate_register_on_error;
			}
		} else
			preffit = gfit;
		if (preffit->naxes[2] == 1 && preffit->keywords.bayer_pattern[0] != '\0') {
			// TODO: same remark as in gui/registration.c
			sensor_pattern pattern = get_cfa_pattern_index_from_string(preffit->keywords.bayer_pattern);
			if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
				siril_log_color_message(_("Cannot use drizzle on non-bayer sensors, aborting.\n"), "red");
				clearfits(preffit);
				goto terminate_register_on_error;
			}
			driz->is_bayer = TRUE;
		}
		if (preffit != gfit)
			clearfits(preffit);
	} else {
		free(driz);
		driz = NULL;
	}

	/* getting the selected registration method */
	method = calloc(1, sizeof(struct registration_method));
	if (regargs->two_pass) {
		method->name = _("Two-Pass Global Star Alignment (deep-sky)");
		method->method_ptr = &register_multi_step_global;
	} else {
		method->name = _("Global Star Alignment (deep-sky)");
		method->method_ptr = &register_star_alignment;
	}
	method->sel = REQUIRES_NO_SELECTION;
	method->type = REGTYPE_DEEPSKY;
	regargs->func = method->method_ptr;

	if (regargs->no_output && regargs->output_scale != 1.f) {
		siril_log_color_message(_("Scaling a sequence with -2pass has no effect, ignoring\n"), "salmon");
	}

	if (regargs->interpolation == OPENCV_NONE && !(regargs->type == SHIFT_TRANSFORMATION)) {
#ifdef HAVE_CV44
		regargs->type = SHIFT_TRANSFORMATION;
		siril_log_color_message(_("Forcing the registration transformation to shift, which is the only transformation compatible with no interpolation\n"), "salmon");

#else
		siril_log_color_message(_("Forcing the registration transformation to shift, which is the only transformation compatible with no interpolation, is not compatible with OpenCV below 4.4. Aborting\n"), "red");
		goto terminate_register_on_error;
#endif
	}

	if (regargs->interpolation == OPENCV_NONE && (regargs->output_scale != 1.f || regargs->seq->is_variable)) {
		siril_log_color_message(_("When interpolation is set to None, the images must be of same size and no scaling can be applied. Aborting\n"), "red");
		goto terminate_register_on_error;
	}

	get_the_registration_area(regargs, method);	// sets selection
	regargs->run_in_thread = TRUE;
	regargs->load_new_sequence = FALSE;	// don't load it for command line execution

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"), "green", method->name);
	free(method);
	msg[strlen(msg) - 1] = '\0';

	if (regargs->interpolation == OPENCV_AREA ||
			regargs->interpolation == OPENCV_LINEAR ||
			regargs->interpolation == OPENCV_NEAREST ||
			regargs->interpolation == OPENCV_NONE ||
			regargs->no_output || drizzle)
		regargs->clamp = FALSE;
	if (regargs->clamp)
		siril_log_message(_("Interpolation clamping active\n"));

	set_progress_bar_data(msg, PROGRESS_RESET);

	if (!start_in_new_thread(register_thread_func, regargs)) {
		errval = CMD_GENERIC_ERROR;
		goto terminate_register_on_error;
	}
	return CMD_OK;

terminate_register_on_error:
	free(regargs->prefix);
	free(regargs);
	if (!check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
	}
	free(driz);
	free(method);
	return errval;
}

// returns 1 if arg was parsed
static int parse_filter_args(char *current, struct seq_filter_config *arg) {
	char *value;
	if (g_str_has_prefix(current, "-filter-fwhm=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_fwhm_p = val;
				arg->f_fwhm_k = (*end == 'k');
			}
			else arg->f_fwhm = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-wfwhm=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_wfwhm_p = val;
				arg->f_wfwhm_k = (*end == 'k');
			}
			else arg->f_wfwhm = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-round=") ||
			g_str_has_prefix(current, "-filter-roundness=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_round_p = val;
				arg->f_round_k = (*end == 'k');
			}
			else arg->f_round = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-qual=") ||
			g_str_has_prefix(current, "-filter-quality=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_quality_p = val;
				arg->f_quality_k = (*end == 'k');
			}
			else arg->f_quality = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-bkg=") ||
			g_str_has_prefix(current, "-filter-background=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_bkg_p = val;
				arg->f_bkg_k = (*end == 'k');
			}
			else arg->f_bkg = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-nbstars=")) {
		value = strchr(current, '=') + 1;
		if (value[0] != '\0') {
			char *end;
			float val = strtof(value, &end);
			if (end == value) {
				siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
				return CMD_ARG_ERROR;
			}
			if (*end == '%' || *end == 'k') {
				arg->f_nbstars_p = val;
				arg->f_nbstars_k = (*end == 'k');
			}
			else arg->f_nbstars = val;
		} else {
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	} else if (g_str_has_prefix(current, "-filter-incl") ||
			g_str_has_prefix(current, "-filter-included")) {
		arg->filter_included = TRUE;
	}
	else return 0;
	return 1;
}

int process_seq_applyreg(int nb) {
	cmd_errors errval = CMD_ARG_ERROR;
	struct registration_args *regargs = NULL;
	gboolean drizzle = FALSE;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;

	regargs = calloc(1, sizeof(struct registration_args));
	struct driz_args_t *driz = calloc(1, sizeof(struct driz_args_t));
	// Default values for the driz_args_t
	driz->use_flats = FALSE;
	driz->scale = 1.f;
	driz->kernel = kernel_square;
	driz->weight_scale = 1.f;
	driz->pixel_fraction = 1.f;

	// check that registration exists for one layer at least
	int layer = -1;
	if (seq->regparam) {
		for (int i = 0; i < seq->nb_layers; i++) {
			if (seq->regparam[i]) {
				layer = i;
				break;
			}
		}
	}
	if (layer == -1) {
		siril_log_color_message(_("No registration data exists for this sequence, aborting\n"), "red");
		goto terminate_register_on_error;
	}

	/* filling the arguments for registration */
	regargs->func = &register_apply_reg;
	regargs->seq = seq;
	regargs->reference_image = sequence_find_refimage(seq);
	regargs->no_output = FALSE;
	regargs->prefix = strdup("r_");
	regargs->layer = layer;
	regargs->interpolation = OPENCV_LANCZOS4;
	regargs->clamp = TRUE;
	regargs->framing = FRAMING_CURRENT;
	regargs->output_scale = 1.f;

	/* check for options */
	for (int i = 2; i < nb; i++) {
		if (!strcmp(word[i], "-drizzle")) {
			if (regargs->seq->nb_layers != 1) {  // handling mono case
				siril_log_message(_("This sequence is not mono / CFA, cannot drizzle.\n"));
				goto terminate_register_on_error;
			}
			drizzle = TRUE;
		// Drizzle options
		} else if (g_str_has_prefix(word[i], "-scale=")) {
			char *arg = word[i] + 7;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg || value < 0.1 || value > 3.) {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
			regargs->output_scale = (float)value;
			driz->scale = (float)value;
		} else if (g_str_has_prefix(word[i], "-pixfrac=")) {
			char *arg = word[i] + 9;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
			driz->pixel_fraction = (float) value;
		} else if (g_str_has_prefix(word[i], "-kernel=")) {
			char *arg = word[i] + 8;
			if (!g_ascii_strncasecmp(arg, "point", 5))
				driz->kernel = kernel_point;
			else if (!g_ascii_strncasecmp(arg, "turbo", 5))
				driz->kernel = kernel_turbo;
			else if (!g_ascii_strncasecmp(arg, "square", 6))
				driz->kernel = kernel_square;
			else if (!g_ascii_strncasecmp(arg, "gaussian", 8))
				driz->kernel = kernel_gaussian;
			else if (!g_ascii_strncasecmp(arg, "lanczos2", 8))
				driz->kernel = kernel_lanczos2;
			else if (!g_ascii_strncasecmp(arg, "lanczos3", 8))
				driz->kernel = kernel_lanczos3;
			else {
				siril_log_color_message(_("Invalid argument to %s, aborting.\n"), "red", word[i]);
				goto terminate_register_on_error;
			}
		} else if (g_str_has_prefix(word[i], "-flat=")) {
			if (driz->flat) {
				siril_log_color_message(_("Error: flat image already set. Aborting.\n"), "red");
				goto terminate_register_on_error;
			}
			if (seq->is_variable) {
				siril_log_color_message(_("Error: flat image cannot work with variable sized sequence.\n"), "red");
				goto terminate_register_on_error;
			}
			char *flat_filename = word[i] + 6;
			fits reffit = { 0 };
			gchar *error = NULL;
			int status;
			if (seq_read_frame_metadata(seq, seq->reference_image, &reffit)) {
				siril_log_color_message(_("NOT USING FLAT: Could not load reference image\n"), "red");
				clearfits(&reffit);
				goto terminate_register_on_error;
			}
			gchar *expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
			clearfits(&reffit);
			if (status) {
				g_free(expression);
				error = _("NOT USING FLAT: could not parse the expression");
				// no need to free expression as we don't call path_parse with the NOFAIL mode
				goto terminate_register_on_error;
			} else {
				if (expression[0] == '\0') {
					siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
					g_free(expression);
					goto terminate_register_on_error;
				} else {
					driz->flat = calloc(1, sizeof(fits));
					if (!readfits(expression, driz->flat, NULL, TRUE)) {
						if (driz->flat->naxes[2] != seq->nb_layers) {
							error = _("NOT USING FLAT: number of channels is different");
						} else if (driz->flat->naxes[0] != seq->rx ||
								driz->flat->naxes[1] != seq->ry) {
							error = _("NOT USING FLAT: image dimensions are different");
						} else {
							// no need to deal with bitdepth conversion as readfits has already forced conversion to float
							siril_log_message(_("Master flat read for use as initial pixel weight\n"));
						}

					} else error = _("NOT USING FLAT: cannot open the file");
					g_free(expression);
					if (error) {
						goto terminate_register_on_error;
					}
				}
			}
		// Other registration options
		} else if (g_str_has_prefix(word[i], "-prefix=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			regargs->prefix = strdup(value);
		} else if (!strcmp(word[i], "-noclamp")) {
			regargs->clamp = FALSE;
		} else if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				regargs->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				regargs->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				regargs->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				regargs->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "none", 4) || !g_ascii_strncasecmp(value, "no", 2)) {
				regargs->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				regargs->interpolation = OPENCV_AREA;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-framing=")) {
			char *current = word[i], *value;
			value = current + 9;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "current", 7)) {
				regargs->framing = FRAMING_CURRENT;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "min", 3)) {
				regargs->framing = FRAMING_MIN;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "max", 3)) {
				regargs->framing = FRAMING_MAX;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cog", 3)) {
				regargs->framing = FRAMING_COG;
				continue;
			}
			siril_log_message(_("Unknown framing type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-layer=")) {
			if (regargs->seq->nb_layers == 1) {  // handling mono case
				siril_log_message(_("This sequence is mono, ignoring layer number.\n"));
				continue;
			}
			char *current = word[i], *value;
			value = current + 7;
			gchar *end;
			int layer2 = g_ascii_strtoull(value, &end, 10);
			if (end == value) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), value);
				continue;
			}
			if (!(seq->regparam[layer2])) {
				siril_log_color_message(_("Registration data does not exist for layer #%d, will use layer #%d instead.\n"), "red", layer2, layer);
				continue;
			}
			regargs->layer = layer2;
		} else if (parse_filter_args(word[i], &regargs->filters)) {
			;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			goto terminate_register_on_error;
		}
	}

	if (drizzle) {
		regargs->driz = driz;
		// we now check the distortion params are ok
		fits reffit = { 0 };
		fits *preffit = &reffit;
		if (!check_seq_is_comseq(seq)) { // processing an image from the current sequence
			int image_to_load = sequence_find_refimage(seq);
			if (seq_read_frame_metadata(seq, image_to_load, preffit)) {
				siril_log_color_message(_("Could not load the reference image of the sequence, aborting.\n"), "red");
				goto terminate_register_on_error;
			}
		} else
			preffit = gfit;
		if (preffit->naxes[2] == 1 && preffit->keywords.bayer_pattern[0] != '\0') {
			sensor_pattern pattern = get_cfa_pattern_index_from_string(preffit->keywords.bayer_pattern);
			if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
				siril_log_color_message(_("Cannot use drizzle on non-bayer sensors, aborting.\n"), "red");
				clearfits(preffit);
				goto terminate_register_on_error;
			}
			driz->is_bayer = TRUE;
		}
		clearfits(preffit);
	} else {
		free(driz);
		driz = NULL;
	}
	// sanity checks are done in register_apply_reg

	regargs->run_in_thread = TRUE;
	regargs->load_new_sequence = FALSE;	// don't load it for command line execution

	if (regargs->interpolation == OPENCV_AREA || regargs->interpolation == OPENCV_LINEAR || regargs->interpolation == OPENCV_NEAREST || regargs->interpolation == OPENCV_NONE)
		regargs->clamp = FALSE;
	if (regargs->clamp && !drizzle)
		siril_log_message(_("Interpolation clamping active\n"));

	set_progress_bar_data(_("Registration: Applying existing data"), PROGRESS_RESET);

	if (!start_in_new_thread(register_thread_func, regargs)) {
		errval = CMD_GENERIC_ERROR;
		goto terminate_register_on_error;
	}
	return CMD_OK;

terminate_register_on_error:
	free(driz);
	if (!check_seq_is_comseq(seq))
		free_sequence(seq, TRUE);
	free(regargs->prefix);
	free(regargs->new_seq_name);
	free(regargs);
	return errval;
}

// parse normalization and filters from the stack command line, starting at word `first'
static int parse_stack_command_line(struct stacking_configuration *arg, int first,
		gboolean med_options_allowed, gboolean rej_options_allowed, gboolean out_allowed) {
	while (word[first]) {
		char *current = word[first], *value;
		if (!strcmp(current, "-nonorm") || !strcmp(current, "-no_norm"))
			arg->force_no_norm = TRUE;
		else if (!strcmp(current, "-output_norm")) {
			if (!med_options_allowed) {
				siril_log_message(_("Output normalization is allowed only with median or mean stacking, ignoring.\n"));
			} else {
				arg->output_norm = TRUE;
			}
		} else if (!strcmp(current, "-32b")) {
			arg->force32b = TRUE;
		} else if (!strcmp(current, "-overlap_norm")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Overlap normalization is allowed only with mean stacking, ignoring.\n"));
			} else {
				arg->overlap_norm = TRUE;
			}
		} else if (g_str_has_prefix(current, "-weight=")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Weighting is allowed only with mean stacking, ignoring.\n"));
			} else {
				value = current + 8;
				if (!strcmp(value, "noise"))
					arg->weighting_type = NOISE_WEIGHT;
				else if (!strcmp(value, "nbstars"))
					arg->weighting_type = NBSTARS_WEIGHT;
				else if (!strcmp(value, "nbstack"))
					arg->weighting_type = NBSTACK_WEIGHT;
				else if (!strcmp(value, "wfwhm"))
					arg->weighting_type = WFWHM_WEIGHT;
				else {
					siril_log_message(_("Unknown argument to %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
			}
		} else if (!strcmp(current, "-fastnorm")) {
			if (!med_options_allowed) {
				siril_log_message(_("Fast normalization is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("Fast normalization is allowed only if normalization has been activated, ignoring.\n"));
			} else {
				arg->lite_norm = TRUE;
			}
		} else if (g_str_has_prefix(current, "-norm=")) {
			if (!med_options_allowed) {
				siril_log_message(_("Normalization options are not allowed in this context, ignoring.\n"));
			} else {
				value = current + 6;
				if (!strcmp(value, "add"))
					arg->norm = ADDITIVE;
				else if (!strcmp(value, "addscale"))
					arg->norm = ADDITIVE_SCALING;
				else if (!strcmp(value, "mul"))
					arg->norm = MULTIPLICATIVE;
				else if (!strcmp(value, "mulscale"))
					arg->norm = MULTIPLICATIVE_SCALING;
			}
		} else if (g_str_has_prefix(current, "-feather=")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Blending option is not allowed in this context, ignoring.\n"));
			} else {
				gchar *end;
				value = current + 9;
				int dist = g_ascii_strtoull(value, &end, 10);
				if (end == value || dist < 0) {
					siril_log_message(_("Unknown argument to %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
				if (dist > 2000) {
					siril_log_message(_("Blending distance must be between 0 and 2000 pixels, got %d, forcing to 2000.\n"), value);
					dist = 2000;
				}
				arg->feather_dist = dist;
			}
		} else if (!strcmp(current, "-rgb_equal")) {
			if (!med_options_allowed) {
				siril_log_message(_("RGB equalization is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("RGB equalization is allowed only if normalization has been activated, ignoring.\n"));
			} else {
				arg->equalizeRGB = TRUE;
			}
		} else if (parse_filter_args(current, &arg->filters)) {
			;
		} else if (g_str_has_prefix(current, "-out=")) {
			if (out_allowed) {
				value = current + 5;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
				arg->result_file = g_strdup(value);
			}
			else {
				siril_log_message(_("Output filename option is not allowed in this context, ignoring.\n"));
			}
		} else if (g_str_has_prefix(current, "-rejmap")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Rejection maps can only be created with rejection stacking, ignoring.\n"));
			} else if (arg->type_of_rejection == NO_REJEC) {
				siril_log_message(_("Rejection maps can only be created if rejection has been activated, ignoring.\n"));
			} else {
				arg->create_rejmaps = TRUE;
				arg->merge_lowhigh_rejmaps = TRUE;
				if (current[strlen(current)-1] == 's')
					arg->merge_lowhigh_rejmaps = FALSE;
			}
		} else if (!strcmp(current, "-maximize")) {
			arg->maximize_framing = TRUE;
		} else if (!strcmp(current, "-upscale")) {
			arg->upscale_at_stacking = TRUE;
		} else {
			siril_log_message(_("Unexpected argument to stacking `%s', aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
		first++;
	}
	return CMD_OK;
}

static gboolean stack_stop_thread(gpointer user_data) {
	stop_processing_thread();
	return FALSE;
}

static int stack_one_seq(struct stacking_configuration *arg) {
	sequence *seq = readseqfile(arg->seqfile);
	if (!seq) {
		siril_log_message(_("No sequence `%s' found.\n"), arg->seqfile);
		return -1;
	}
	gchar *error = NULL;
	if (seq_check_basic_data(seq, FALSE) == -1) {
		free(seq);
		return CMD_GENERIC_ERROR;
	}
	siril_log_message(_("Stacking sequence %s\n"), seq->seqname);

	struct stacking_args args = { 0 };
	init_stacking_args(&args);
	args.seq = seq;
	args.ref_image = sequence_find_refimage(seq);
	// the three below: used only if method is average w/ rejection
	if (arg->method == stack_mean_with_rejection && (arg->sig[0] != 0.0 || arg->sig[1] != 0.0)) {
		args.sig[0] = arg->sig[0];
		args.sig[1] = arg->sig[1];
		args.type_of_rejection = arg->type_of_rejection;
		args.create_rejmaps = arg->create_rejmaps;
		args.merge_lowhigh_rejmaps = arg->merge_lowhigh_rejmaps;
	} else {
		args.type_of_rejection = NO_REJEC;
		siril_log_message(_("Not using rejection for stacking\n"));
	}
	if (!arg->force_no_norm &&
			(arg->method == stack_median || arg->method == stack_mean_with_rejection))
		args.normalize = arg->norm;
	else args.normalize = NO_NORM;
	args.method = arg->method;
	args.output_norm = arg->output_norm;
	args.equalizeRGB = arg->equalizeRGB;
	args.lite_norm = arg->lite_norm;
	args.reglayer = get_registration_layer(args.seq);
	args.feather_dist = arg->feather_dist;
	args.overlap_norm = arg->overlap_norm;
	args.weighting_type = (arg->method == stack_mean_with_rejection) ? arg->weighting_type : NO_WEIGHT;

	// manage registration data
	if (!test_regdata_is_valid_and_shift(args.seq, args.reglayer)) {
		siril_log_color_message(_("Stacking has detected registration data on layer %d with more than simple shifts. You should apply existing registration before stacking\n"), "red", args.reglayer);
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	// manage reframing and upscale
	gboolean can_reframe = layer_has_usable_registration(seq, args.reglayer);
	gboolean can_upscale = can_reframe && !seq->is_variable && !seq->is_drizzle;
	gboolean must_reframe = can_reframe && seq->is_variable;
	args.maximize_framing = arg->maximize_framing;
	args.upscale_at_stacking = arg->upscale_at_stacking;
	if (args.maximize_framing && !can_reframe) {
		siril_log_color_message(_("No registration data in the sequence. Maximize framing will be ignored\n"), "red");
		args.maximize_framing = FALSE;
	}
	if (!args.maximize_framing && must_reframe) {
		siril_log_color_message(_("The sequence has different image sizes and registration data. Forcing to maximize framing\n"), "red");
		args.maximize_framing = TRUE;
	}
	if (!can_reframe && seq->is_variable) {
		siril_log_color_message(_("The sequence has different image sizes but no registration data, cannot stack. Aborting\n"), "red");
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	if (args.upscale_at_stacking && !can_upscale) {
		siril_log_color_message(_("No registration data in the sequence or images with different sizes or drizzled. Upscale at stacking will be ignored\n"), "red");
		args.upscale_at_stacking = FALSE;
	}
	if ((args.upscale_at_stacking || args.maximize_framing) && arg->method == stack_median) {
		siril_log_color_message(_("Cannot upscale or maximize framing with median stacking. Disabling\n"), "red");
		args.maximize_framing = FALSE;
		args.upscale_at_stacking = FALSE;
	}
	if (!args.maximize_framing && args.overlap_norm) {
		siril_log_color_message(_("Cannot compute overlap statistics if -maximize is not enabled. Disabling\n"), "red");
		args.overlap_norm = FALSE;
	}
	if (args.normalize == NO_NORM && (args.weighting_type == NOISE_WEIGHT)) {
		siril_log_color_message(_("Weighting by noise is allowed only if normalization has been activated, ignoring weights.\n"), "red");
		args.weighting_type = NO_WEIGHT;
	}
	if (args.overlap_norm && (args.weighting_type == NOISE_WEIGHT)) {
		siril_log_color_message(_("Weighting by noise cannot be used with overlap normalization, ignoring weights.\n"), "red");
		args.weighting_type = NO_WEIGHT;
	}


	// manage filters
	if (convert_parsed_filter_to_filter(&arg->filters, seq,
				&args.filtering_criterion, &args.filtering_parameter) ||
			setup_filtered_data(&args)) {
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	args.description = describe_filter(seq, args.filtering_criterion, args.filtering_parameter);
	args.use_32bit_output = arg->force32b || evaluate_stacking_should_output_32bits(args.method,
			args.seq, args.nb_images_to_stack, &error);
	if (error) {
		siril_log_color_message(error, "red");
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	if (args.overlap_norm && args.nb_images_to_stack > MAX_IMAGES_FOR_OVERLAP) {
		siril_log_color_message(_("Normalizing on overlaps for more than %d images can be slow\n"), "salmon", MAX_IMAGES_FOR_OVERLAP);
	}

	main_stack(&args);

	int retval = args.retval ? CMD_GENERIC_ERROR : CMD_OK;
	clean_end_stacking(&args);
	free(args.image_indices);
	free(args.description);
	free(args.critical_value);

	if (retval == CMD_OK) {
		execute_idle_and_wait_for_it(stack_stop_thread, NULL);
		bgnoise_async(&args.result, TRUE);
		// preparing the output filename
		// needs to be done after stack is completed to have
		// stack-specific keywords available for parsing if needed
		long maxpath = get_pathmax();
		if (!arg->result_file) {
			char filename[maxpath];
			char *suffix = g_str_has_suffix(seq->seqname, "_") ||
				g_str_has_suffix(seq->seqname, "-") ? "" : "_";
			snprintf(filename, maxpath - 1, "%s%sstacked%s", seq->seqname, suffix, com.pref.ext);
			arg->result_file = g_strdup(filename);
		} else { // the name is to be parsed (including folder creation if required)
			int status = PATHPARSE_ERR_OK;
			gchar *expression = g_strdup(arg->result_file);
			gchar *parsedname = update_header_and_parse(&args.result, expression, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
			char filename[maxpath];
			if (!parsedname || parsedname[0] == '\0') { // we cannot handout a NULL filename
				snprintf(filename, maxpath - 1, "unknown");
			} else {
				snprintf(filename, maxpath - 1, "%s", parsedname);
			}
			g_free(arg->result_file);
			arg->result_file = g_strdup(filename);
			g_free(parsedname);
			g_free(expression);
		}
		/* Make sure path exists */
		gchar *dirname = g_path_get_dirname(arg->result_file);
		if (siril_mkdir_with_parents(dirname, 0755) < 0) {
			g_free(dirname);
			retval = CMD_GENERIC_ERROR;
		}
		g_free(dirname);
		if (savefits(arg->result_file, &args.result)) {
			siril_log_color_message(_("Could not save the stacking result %s\n"),
					"red", arg->result_file);
			retval = CMD_GENERIC_ERROR;
		}
		else ++arg->number_of_loaded_sequences;

		if (args.create_rejmaps) {
			siril_log_message(_("Saving rejection maps\n"));
			if (args.merge_lowhigh_rejmaps) {
				char new_ext[30];
				sprintf(new_ext, "_low+high_rejmap%s", com.pref.ext);
				gchar *low_filename = replace_ext(arg->result_file, new_ext);
				soper_unscaled_div_ushort_to_float(args.rejmap_low, args.nb_images_to_stack);
				describe_stack_for_history(&args, &args.rejmap_low->history, TRUE, FALSE);
				savefits(low_filename, args.rejmap_low);
				g_free(low_filename);
			} else {
				char new_ext[30];
				sprintf(new_ext, "_low_rejmap%s", com.pref.ext);
				gchar *low_filename = replace_ext(arg->result_file, new_ext);
				soper_unscaled_div_ushort_to_float(args.rejmap_low, args.nb_images_to_stack);
				describe_stack_for_history(&args, &args.rejmap_low->history, TRUE, TRUE);
				savefits(low_filename, args.rejmap_low);
				g_free(low_filename);

				sprintf(new_ext, "_high_rejmap%s", com.pref.ext);
				gchar *high_filename = replace_ext(arg->result_file, new_ext);
				soper_unscaled_div_ushort_to_float(args.rejmap_high, args.nb_images_to_stack);
				describe_stack_for_history(&args, &args.rejmap_low->history, TRUE, FALSE);
				savefits(high_filename, args.rejmap_high);
				g_free(high_filename);
			}
		}

		bgnoise_await();
	} else {
		siril_log_color_message(_("Stacking failed, please check the log to fix your issue.\n"), "red");
	}

	free_sequence(seq, TRUE);
	clearfits(&args.result);
	if (args.create_rejmaps) {
		clearfits(args.rejmap_low);
		free(args.rejmap_low);
		if (!args.merge_lowhigh_rejmaps) {
			clearfits(args.rejmap_high);
			free(args.rejmap_high);
		}
	}
	return retval;
}

static gpointer stackall_worker(gpointer garg) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	struct timeval t_end;
	struct stacking_configuration *arg = (struct stacking_configuration *)garg;
	gboolean was_in_script = com.script;
	com.script = TRUE;

	siril_log_message(_("Looking for sequences in current working directory...\n"));
	if (check_seq() || (dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		siril_log_message(_("Error while searching sequences or opening the directory.\n"));
		if (error) {
			fprintf(stderr, "stackall: %s\n", error->message);
			g_clear_error(&error);
		}
		siril_add_idle(end_generic, NULL);
		return NULL;
	}

	siril_log_message(_("Starting stacking of found sequences...\n"));
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_str_has_suffix(file, ".seq")) {
			arg->seqfile = strdup(file);
			stack_one_seq(arg);

			g_free(arg->result_file);
			arg->result_file = NULL;
			g_free(arg->seqfile);
		}
	}

	siril_log_message(_("Stacked %d sequences successfully.\n"), arg->number_of_loaded_sequences);
	gettimeofday(&t_end, NULL);
	show_time(arg->t_start, t_end);
	g_dir_close(dir);
	free(arg);
	com.script = was_in_script;
	siril_add_idle(end_generic, NULL);
	return NULL;
}

int process_stackall(int nb) {
	struct stacking_configuration *arg;

	arg = calloc(1, sizeof(struct stacking_configuration));
	arg->norm = NO_NORM;

	// stackall { sum | min | max } [-filter-fwhm=value[%|k]] [-filter-wfwhm=value[%|k]] [-filter-round=value[%|k]] [-filter-quality=value[%|k]] [-filter-bkg=value[%|k]] [-filter-nbstars=value[%|k]] [-filter-incl[uded]]
	// stackall { med | median } [-nonorm, norm=] [-filter-incl[uded]]
	// stackall { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%|k]] [-filter-round=value[%|k]] [-filter-bkg=value[%|k]] [-filter-nbstars=value[%|k]] [-filter-quality=value[%|k]] [-filter-incl[uded]] [-weighted]
	if (!word[1]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 2;
		gboolean allow_rej_options = FALSE, allow_med_options = FALSE;
		if (!strcmp(word[1], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[1], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[1], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[1], "med") || !strcmp(word[1], "median")) {
			arg->method = stack_median;
			allow_med_options = TRUE;
		} else if (!strcmp(word[1], "rej") || !strcmp(word[1], "mean")) {
			int shift = 1;
			gchar *end1, *end2;

			if (!strcmp(word[3], "p") || !strcmp(word[3], "percentile")) {
				arg->type_of_rejection = PERCENTILE;
			} else if (!strcmp(word[3], "s") || !strcmp(word[3], "sigma")) {
				arg->type_of_rejection = SIGMA;
			} else if (!strcmp(word[3], "a") || !strcmp(word[3], "mad")) {
				arg->type_of_rejection = MAD;
			} else if (!strcmp(word[3], "m") || !strcmp(word[3], "median")) {
				arg->type_of_rejection = SIGMEDIAN;
			} else if (!strcmp(word[3], "l") || !strcmp(word[3], "linear")) {
				arg->type_of_rejection = LINEARFIT;
			} else if (!strcmp(word[3], "w") || !strcmp(word[3], "winsorized")) {
				arg->type_of_rejection = WINSORIZED;
			} else if (!strcmp(word[3], "g") || !strcmp(word[3], "generalized")) {
				arg->type_of_rejection = GESDT;
			} else {
				arg->type_of_rejection = WINSORIZED;
				shift = 0;
			}
			if (!word[2 + shift] || !word[3 + shift] || (arg->sig[0] = g_ascii_strtod(word[2 + shift], &end1)) < 0.0
					|| (arg->sig[1] = g_ascii_strtod(word[3 + shift], &end2)) < 0.0 || end1 == word[2 + shift] || end2 == word[3 + shift]) {
				siril_log_color_message(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"), "red");
				goto failure;
			}
			if (((arg->type_of_rejection == GESDT))
					&& (arg->sig[0] > 1.0 || (arg->sig[1] > 1.0))) {
				siril_log_color_message(_("Extra parameters of GESDT rejection algorithm must be between 0 and 1, default is 0.3 and 0.05.\n"), "red");
				goto failure;
			}
			if (((arg->type_of_rejection == PERCENTILE))
					&& (arg->sig[0] > 1.0 || (arg->sig[1] > 1.0))) {
				siril_log_color_message(_("Extra parameters of percentile rejection algorithm must be between 0 and 1, default is 0.2 and 0.1.\n"), "red");
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = 4 + shift;
			allow_med_options = TRUE;
			allow_rej_options = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_med_options, allow_rej_options, FALSE))
			goto failure;
	}

	gettimeofday(&arg->t_start, NULL);

	if (!start_in_new_thread(stackall_worker, arg)) {
		g_free(arg->seqfile);
		g_free(arg->result_file);
		free(arg);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;

failure:
	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	return CMD_ARG_ERROR;
}

static gpointer stackone_worker(gpointer garg) {
	int retval = 0;
	struct timeval t_end;
	struct stacking_configuration *arg = (struct stacking_configuration *)garg;
	gboolean was_in_script = com.script;
	com.script = TRUE;

	retval = stack_one_seq(arg);

	if (retval) {
		if (retval == ST_ALLOC_ERROR) {
			siril_log_message(_("It looks like there is a memory allocation error, change memory settings and try to fix it.\n"));
		}
	} else {
		siril_log_message(_("Stacked sequence successfully.\n"));
	}

	gettimeofday(&t_end, NULL);
	show_time(arg->t_start, t_end);

	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	com.script = was_in_script;
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

int process_stackone(int nb) {
	struct stacking_configuration *arg = calloc(1, sizeof(struct stacking_configuration));
	arg->norm = NO_NORM;

	sequence *seq = load_sequence(word[1], &arg->seqfile);
	if (!seq)
		goto failure;
	free_sequence(seq, TRUE);

	// stack seqfilename { sum | min | max } [-filter-fwhm=value[%|k]] [-filter-wfwhm=value[%|k]] [-filter-round=value[%|k]] [-filter-quality=value[%|k]] [-filter-bkg=value[%|k]] [-filter-nbstars=value[%|k]] [-filter-incl[uded]] [-out=result_filename]
	// stack seqfilename { med | median } [-nonorm, norm=] [-filter-incl[uded]] [-out=result_filename]
	// stack seqfilename { rej | mean } [type_of_rejection] sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%|k]] [-filter-round=value[%|k]] [-filter-quality=value[%|k]] [-filter-bkg=value[%|k]] [-filter-nbstars=value[%|k]] [-filter-incl[uded]] [-weighted] [-out=result_filename]
	if (!word[2]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 3;
		gchar *end1, *end2;
		gboolean allow_rej_options = FALSE, allow_med_options = FALSE;
		if (!strcmp(word[2], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[2], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[2], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[2], "med") || !strcmp(word[2], "median")) {
			arg->method = stack_median;
			allow_med_options = TRUE;
		} else if (!strcmp(word[2], "rej") || !strcmp(word[2], "mean")) {
			int shift = 1, base_shift = 5;
			if (nb < 4) {
				siril_log_color_message(_("Missing arguments for rejection stacking.\n"), "red");
				goto failure;
			}
			if (!strcmp(word[3], "p") || !strcmp(word[3], "percentile")) {
				arg->type_of_rejection = PERCENTILE;
			} else if (!strcmp(word[3], "s") || !strcmp(word[3], "sigma")) {
				arg->type_of_rejection = SIGMA;
			} else if (!strcmp(word[3], "a") || !strcmp(word[3], "mad")) {
				arg->type_of_rejection = MAD;
			} else if (!strcmp(word[3], "m") || !strcmp(word[3], "median")) {
				arg->type_of_rejection = SIGMEDIAN;
			} else if (!strcmp(word[3], "l") || !strcmp(word[3], "linear")) {
				arg->type_of_rejection = LINEARFIT;
			} else if (!strcmp(word[3], "w") || !strcmp(word[3], "winsorized")) {
				arg->type_of_rejection = WINSORIZED;
			} else if (!strcmp(word[3], "g") || !strcmp(word[3], "generalized")) {
				arg->type_of_rejection = GESDT;
			} else if (!strcmp(word[3], "n") || !strcmp(word[3], "none")) {
				arg->type_of_rejection = NO_REJEC;
			} else {
				arg->type_of_rejection = WINSORIZED;
				shift = 0;
			}
			if ((nb < 4 + shift + 1) || !word[3 + shift] || !word[4 + shift] ||
					!string_is_a_number(word[3 + shift]) ||
					!string_is_a_number(word[4 + shift]) ||
					(arg->sig[0] = g_ascii_strtod(word[3 + shift], &end1)) < 0.0 ||
					(arg->sig[1] = g_ascii_strtod(word[4 + shift], &end2)) < 0.0 ||
					end1 == word[3 + shift] || end2 == word[4 + shift]) {
				if (arg->type_of_rejection != NO_REJEC) {
					siril_log_color_message(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"), "red");
					goto failure;
				} else {
					base_shift = 3;
				}
			}
			if (((arg->type_of_rejection == GESDT))
					&& (arg->sig[0] > 1.0 || (arg->sig[1] > 1.0))) {
				siril_log_color_message(_("Extra parameters of GESDT rejection algorithm must be between 0 and 1, default is 0.3 and 0.05.\n"), "red");
				goto failure;
			}
			if (((arg->type_of_rejection == PERCENTILE))
					&& (arg->sig[0] > 1.0 || (arg->sig[1] > 1.0))) {
				siril_log_color_message(_("Extra parameters of percentile rejection algorithm must be between 0 and 1, default is 0.2 and 0.1.\n"), "red");
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = base_shift + shift;
			allow_med_options = TRUE;
			allow_rej_options = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_med_options, allow_rej_options, TRUE))
			goto failure;
	}

	gettimeofday(&arg->t_start, NULL);

	if (!start_in_new_thread(stackone_worker, arg)) {
		g_free(arg->result_file);
		g_free(arg->seqfile);
		free(arg);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;

failure:
	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	return CMD_ARG_ERROR;
}

/* calibrate sequencename [-bias=filename|value] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt[=exp]] [-prefix=] [-fitseq]
 * calibrate_single filename [-bias=filename|value] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=]
 */
struct preprocessing_data *parse_calibrate_args(int nb, sequence *seq) {
	int retvalue = 0;
	struct preprocessing_data *args = calloc(1, sizeof(struct preprocessing_data));
	fits reffit = { 0 };
	int bitpix;
	char *realname = NULL;
	image_type imagetype;
	if (seq) {
		if (seq->type == SEQ_SER) {
			// to be able to check allow_32bit_output. Overridden by -fitseq if required
			args->output_seqtype = SEQ_SER;
		}
		args->seq = seq;
		args->is_sequence = TRUE;
		bitpix = seq->bitpix;
		// loading the sequence reference image's metadata in case it's needed
		int image_to_load = sequence_find_refimage(seq);
		if (seq_read_frame_metadata(seq, image_to_load, &reffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			retvalue = CMD_INVALID_IMAGE;
			goto prepro_parse_end;
		}
	}
	else {
		if (stat_file(word[1], &imagetype, &realname)) {
			siril_log_color_message(_("Error opening image %s: file not found or not supported.\n"), "red", word[1]);
			retvalue = CMD_FILE_NOT_FOUND;
			goto prepro_parse_end;
		}
		if (read_fits_metadata_from_path(realname, &reffit)) {
			siril_log_message(_("Could not load the image, aborting.\n"));
			retvalue = CMD_INVALID_IMAGE;
			goto prepro_parse_end;
		}
		bitpix = reffit.bitpix;
	}
	args->ppprefix = "pp_";
	args->bias_level = FLT_MAX;

	//cosmetic_correction init
	args->sigma[0] = -1.00; /* cold pixels */
	args->sigma[1] =  3.00; /* hot pixels - used if -cc=dark but not sigmas specified */
	args->use_cosmetic_correction = FALSE;// dy default, CC is not activated
	args->cc_from_dark = FALSE;
	args->bad_pixel_map_file = NULL;
	args->ignore_exclusion = FALSE;

	/* checking for options */
	for (int i = 2; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-bias=")) {
			gchar *expression = g_shell_unquote(word[i] + 6, NULL);
			if (expression && expression[0] == '=') {
				// parsing offset level
				int offsetlevel = evaluateoffsetlevel(expression + 1, &reffit);
				if (!offsetlevel) {
					siril_log_message(_("The offset value could not be parsed from expression: %s, aborting.\n"), expression +1);
					retvalue = 1;
					g_free(expression);
					break;
				} else {
					g_free(expression);
					siril_log_message(_("Synthetic offset: Level = %d\n"), offsetlevel);
					int maxlevel = (bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
					if ((offsetlevel > maxlevel) || (offsetlevel < -maxlevel) ) {   // not excluding all neg values here to allow defining a pedestal
						siril_log_message(_("The offset value is out of allowable bounds [-%d,%d], aborting.\n"), maxlevel, maxlevel);
						retvalue = CMD_ARG_ERROR;
						break;
					} else {
						args->bias_level = (float)offsetlevel;
						args->bias_level *= (bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE; //converting to [0 1] to use with soper
						args->use_bias = TRUE;
					}
				}
			} else {
				g_free(expression);
				expression = NULL;
				int status;
				expression = path_parse(&reffit, word[i] + 6, PATHPARSE_MODE_READ, &status);
				if (status) {
					g_free(expression);
					retvalue = CMD_GENERIC_ERROR;
					break;
				}
				args->bias = calloc(1, sizeof(fits));
				if (!readfits(expression, args->bias, NULL, !com.pref.force_16bit)) {
					args->use_bias = TRUE;
					// if input is 8b, we assume 32b master needs to be rescaled
					if ((args->bias->type == DATA_FLOAT) && (bitpix == BYTE_IMG)) {
						soper(args->bias, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
					}
				} else {
					retvalue = CMD_INVALID_IMAGE;
					free(args->bias);
					g_free(expression);
					break;
				}
				g_free(expression);
			}
		} else if (g_str_has_prefix(word[i], "-dark=")) {
			args->dark = calloc(1, sizeof(fits));
			int status;
			gchar *expression = path_parse(&reffit, word[i] + 6, PATHPARSE_MODE_READ, &status);
			if (status > 0) { // negative status are warnings
				retvalue = CMD_GENERIC_ERROR;
				g_free(expression);
				free(args->dark);
				break;
			}
			if (!readfits(expression, args->dark, NULL, !com.pref.force_16bit)) {
				args->use_dark = TRUE;
				// if input is 8b, we assume 32b master needs to be rescaled
				if ((args->dark->type == DATA_FLOAT) && (bitpix == BYTE_IMG)) {
					soper(args->dark, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
				}
			} else {
				retvalue = CMD_INVALID_IMAGE;
				free(args->dark);
				g_free(expression);
				break;
			}
			g_free(expression);
		} else if (g_str_has_prefix(word[i], "-flat=")) {
			args->flat = calloc(1, sizeof(fits));
			int status;
			gchar *expression = path_parse(&reffit, word[i] + 6, PATHPARSE_MODE_READ, &status);
			if (status) {
				retvalue = CMD_GENERIC_ERROR;
				free(args->flat);
				g_free(expression);
				break;
			}
			if (!readfits(expression, args->flat, NULL, !com.pref.force_16bit)) {
				args->use_flat = TRUE;
				// no need to deal with bitdepth conversion as flat is just a division (unlike darks which need to be on same scale)
			} else {
				retvalue = CMD_INVALID_IMAGE;
				free(args->flat);
				g_free(expression);
				break;
			}
			g_free(expression);
		} else if (g_str_has_prefix(word[i], "-prefix=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				retvalue = CMD_ARG_ERROR;
				break;
			}
			args->ppprefix = strdup(value);
		} else if (!strcmp(word[i], "-opt") || g_str_has_prefix(word[i], "-opt=")) {
			if (bitpix == BYTE_IMG) {
				siril_log_color_message(_("Dark optimization: This process cannot be applied to 8b images\n"), "red");
				retvalue = CMD_INVALID_IMAGE;
				break;
			}
			args->use_exposure = FALSE;
			if (word[i][4] == '=') {
				char *current = word[i], *value;
				value = current + 5;
				if (value[0] == '\0' || strcmp(value, "exp")) {
					siril_log_message(_("Missing or wrong argument to %s, aborting.\n"), current);
					retvalue = CMD_ARG_ERROR;
					break;
				}
				args->use_exposure = TRUE;
			}
			args->use_dark_optim = TRUE;
		} else if (!strcmp(word[i], "-fix_xtrans")) {
			args->fix_xtrans = TRUE;
		} else if (!strcmp(word[i], "-all")) {
			args->ignore_exclusion = TRUE;
		} else if (!strcmp(word[i], "-cfa")) {
			args->is_cfa = TRUE;
		} else if (!strcmp(word[i], "-debayer")) {
			args->debayer = TRUE;
		} else if (!strcmp(word[i], "-equalize_cfa")) {
			args->equalize_cfa = TRUE;
		} else if (seq && !strcmp(word[i], "-fitseq")) {
			args->output_seqtype = SEQ_FITSEQ;
		} else if (g_str_has_prefix(word[i], "-cc=")) {
			char *current = word[i], *value;
			value = current + 4;
			args->use_cosmetic_correction = TRUE;// dy default, CC is not activated
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				retvalue = CMD_ARG_ERROR;
				break;
			}
			if (!strcmp(value, "dark")) {
				if (!args->use_dark){
					siril_log_message(_("You must specify a masterdark with -dark= before activating this option, aborting.\n"));
					retvalue = CMD_ARG_ERROR;
					break;
				}
				args->cc_from_dark = TRUE;
				//either we use the default sigmas or we try to read the next two checking for sigma low and high
				if (i + 2 < nb) {
					gchar *end1 = NULL, *end2 = NULL;
					double sigma0 = g_ascii_strtod(word[i + 1], &end1);
					double sigma1 = g_ascii_strtod(word[i + 2], &end2);
					if (word[i + 1] && word[i + 2] && word[i + 1] != end1
							&& word[i + 2] != end2) {
						i += 2;
						args->sigma[0] = (sigma0 == 0.0) ? -1.00 : sigma0;
						args->sigma[1] = (sigma1 == 0.0) ? -1.00 : sigma1;
					}
				}
				if (args->sigma[0] > 0)
					siril_log_message(_("Cosmetic correction from masterdark: using sigma %.2lf for cold pixels.\n"), args->sigma[0]);
				else
					siril_log_message(_("Cosmetic correction from masterdark: deactivated for cold pixels.\n"));
				if (args->sigma[1] > 0)
					siril_log_message(_("Cosmetic correction from masterdark: using sigma %.2lf for hot pixels.\n"), args->sigma[1]);
				else
					siril_log_message(_("Cosmetic correction from masterdark: deactivated for hot pixels.\n"));
			} else if (!strcmp(value, "bpm")) {
				if (word[i + 1] && word[i + 1][0] != '\0') {
					args->bad_pixel_map_file = g_file_new_for_path(word[i + 1]);
					if (!check_for_cosme_file_sanity(args->bad_pixel_map_file)) {
						//g_object_unref(args->bad_pixel_map_file); // This is unreferenced in check_for_cosme_file_sanity
						args->bad_pixel_map_file = NULL;
						siril_log_message(_("Could not open file %s, aborting.\n"), word[i + 1]);
						retvalue = 1;
						break;
					} else {
						siril_log_message(_("Cosmetic correction from Bad Pixel Map: %s.\n"), word[i + 1]);
						i++;
					}
				} else {
					siril_log_message(_("You must specify a bad pixel map file with -cc=bpm option, aborting.\n"));
					retvalue = CMD_ARG_ERROR;
					break;
				}
			} else {
				siril_log_message(_("Unknown argument %s, aborting.\n"), word[i]);
				retvalue = CMD_ARG_ERROR;
				break;
			}
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			retvalue = CMD_ARG_ERROR;
			break;
		}
	}

prepro_parse_end:
	clearfits(&reffit);
	free(realname);
	if (retvalue) {
		clear_preprocessing_data(args);
		free(args);
		return NULL;
	}
	return args;
}

int process_calibrate(int nb) {
	if (word[1][0] == '\0')
		return CMD_ARG_ERROR;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	struct preprocessing_data *args = parse_calibrate_args(nb, seq);
	if (!args) {
		free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	siril_log_color_message(_("Preprocessing...\n"), "green");
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway
	args->allow_32bit_output = (args->output_seqtype == SEQ_REGULAR
			|| args->output_seqtype == SEQ_FITSEQ) && !com.pref.force_16bit;

	start_sequence_preprocessing(args);
	return CMD_OK;
}

int process_calibrate_single(int nb) {
	if (word[1][0] == '\0')
		return CMD_ARG_ERROR;

	struct preprocessing_data *args = parse_calibrate_args(nb, NULL);
	if (!args)
		return CMD_ARG_ERROR;

	siril_log_color_message(_("Preprocessing...\n"), "green");
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway
	args->allow_32bit_output = !com.pref.force_16bit;

	return preprocess_given_image(word[1], args);
}

int process_set_32bits(int nb) {
	com.pref.force_16bit = word[0][3] == '1';
	if (com.pref.force_16bit)
		siril_log_message(_("16-bit per channel in processed images mode is active\n"));
	else siril_log_message(_("32-bit per channel in processed images mode is active\n"));
	writeinitfile();
	return CMD_OK;
}

int process_set_compress(int nb) {
	if (word[1][0] != '0' && word[1][0] != '1' && word[1][0] != 'y' && word[1][0] != 'n' && word[1][1] != '\0') {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	gboolean compress = word[1][0] == '1' || word[1][0] == 'y';
	int method = 0;
	double q = 16.0, hscale= 4.0;

	if (compress) {
		gchar *comp = NULL;
		if (!word[2] || !word[3] || (!g_str_has_prefix(word[2], "-type="))) {
			// take values from pref if arguments are missing
			method = com.pref.comp.fits_method;
			q = com.pref.comp.fits_quantization;
			comp = method == RICE_COMP ? g_strdup("rice") : (method == GZIP1_COMP ? g_strdup("GZIP1") : g_strdup("GZIP2"));
		} else {
			if (!g_ascii_strncasecmp(word[2] + 6, "rice", 4)) {
				method = RICE_COMP;
				comp = g_strdup("rice");
			} else if (!g_ascii_strncasecmp(word[2] + 6, "gzip1", 5)) {
				method = GZIP1_COMP;
				comp = g_strdup("GZIP1");
			} else if (!g_ascii_strncasecmp(word[2] + 6, "gzip2", 5)) {
				method = GZIP2_COMP;
				comp = g_strdup("GZIP2");
			}
			// hcompress with mode 2 is not working in cfitsio
			// see https://gitlab.com/free-astro/siril/-/issues/696#note_932398268
			/*else if (!g_ascii_strncasecmp(word[2] + 6, "hcompress", 9)) {
				method = HCOMPRESS_COMP;
				if (!word[4]) {
					siril_log_message(_("Please specify the value of hcompress scale factor.\n"));
					g_free(comp);
					return CMD_GENERIC_ERROR;
				}
				hscale = g_ascii_strtod(word[4], NULL);
				comp = g_strdup_printf("hcompress (scale factor = %.2lf) ", hscale);
			}*/ else {
	//			siril_log_message(_("Wrong type of compression. Choices are rice, gzip1, gzip2 or hcompress\n"));
				siril_log_message(_("Wrong type of compression. Choices are rice, gzip1, gzip2\n"));
				return CMD_ARG_ERROR;
			}

			gchar *end;
			q = g_ascii_strtod(word[3], &end);
			if (end == word[3] || (q == 0.0 && (method == RICE_COMP || method == HCOMPRESS_COMP))) {
				siril_log_message(_("Quantization can only be equal to 0 for GZIP1 and GZIP2 algorithms.\n"));
				g_free(comp);

				return CMD_ARG_ERROR;
			}
		}
		siril_log_message(_("Compression enabled with the %s algorithm and a quantization value of %.2lf\n"), comp, q);
		g_free(comp);
	} else {
		siril_log_message(_("No compression enabled.\n"));
	}
	com.pref.comp.fits_enabled = compress;
	com.pref.comp.fits_method = method;
	com.pref.comp.fits_quantization = q;
	com.pref.comp.fits_hcompress_scale = hscale;
	writeinitfile();
	return CMD_OK;
}

#ifdef _OPENMP
int process_set_cpu(int nb){
	int proc_in, proc_out, proc_max;
	gchar *end;

	proc_in = g_ascii_strtoull(word[1], &end, 10);
	proc_max = omp_get_num_procs();
	if (end == word[1] || proc_in > proc_max || proc_in < 1) {
		siril_log_message(_("Number of logical processors MUST be greater "
				"than 0 and lower or equal to %d.\n"), proc_max);
		return CMD_ARG_ERROR;
	}
	omp_set_num_threads(proc_in);

#pragma omp parallel
	{
		proc_out = omp_get_num_threads();
	}

	gchar *str = ngettext("Using now %d logical processor\n", "Using now %d logical processors\n", proc_out);
	str = g_strdup_printf(str, proc_out);
	siril_log_message(str);
	g_free(str);

	com.max_thread = proc_out;
	gui_function(update_spinCPU, GINT_TO_POINTER(0));

	return CMD_OK;
}
#endif

int process_set_mem(int nb) {
	gchar *end;
	double ratio = g_ascii_strtod(word[1], &end);
	if (end == word[1] || ratio < 0.05 || ratio > 4.0) {
		siril_log_message(_("The accepted range for the ratio of memory used for stacking is [0.05, 4], with values below the available memory recommended\n"));
		return CMD_ARG_ERROR;
	}
	if (ratio > 1.0) {
		siril_log_message(_("Setting the ratio of memory used for stacking above 1 will require the use of on-disk memory, which can be very slow and is unrecommended (%g requested)\n"), ratio);
	}
	com.pref.memory_ratio = ratio;
	com.pref.mem_mode = RATIO;
	writeinitfile();
	siril_log_message(_("Usable memory for stacking changed to %g\n"), ratio);
	return CMD_OK;
}

int process_help(int nb) {
	command *current = commands;
	if (nb == 1)
		siril_log_message(_("********* LIST OF AVAILABLE COMMANDS *********\n"));
	while (current->process) {
		if (nb == 2) {
			if (!g_ascii_strcasecmp(current->name, word[1])) {
				siril_log_message("%s\n", current->usage);
				char *def = strdup(current->definition);
				log_several_lines(def);
				free(def);
				siril_log_message(_("Can be used in a script: %s\n"), current->scriptable ? _("YES") : _("NO"));
				break;
			}
		}
		else siril_log_message("%s\n", current->usage);
		current++;
	}
	if (nb == 1)
		siril_log_message(_("********* END OF THE LIST *********\n"));
	if (nb == 2 && !current->process)
		siril_log_message(_("Error: command %s is not available\n"), word[1]);
	return CMD_OK;
}

int process_capabilities(int nb) {
	// don't translate these strings, they must be easy to parse
	gchar *version_string = get_siril_version_string();
	siril_log_message("%s\n", version_string);
	g_free(version_string);
#ifdef _OPENMP
	siril_log_message("OpenMP available (%d %s)\n", com.max_thread,
			ngettext("processor", "processors", com.max_thread));
#else
	siril_log_message("OpenMP unavailable\n");
#endif
	siril_log_message("Detected system available memory: %d MB\n", (int)(get_available_memory() / BYTES_IN_A_MB));
	siril_log_message("Can%s create symbolic links\n", test_if_symlink_is_ok(FALSE) ? "" : "not");
#ifdef _WIN32
	int open_max = _getmaxstdio();
	if (open_max < 8192) {
		/* extend the limit to 8192 if possible
		 * 8192 is the maximum on WINDOWS */
		int ret = _setmaxstdio(8192);
		if (ret == -1) {
			/* fallback to 2048 */
			_setmaxstdio(2048);
		}
		open_max = _getmaxstdio();
	}
	siril_log_message("Can open %d files simultaneously\n", open_max);
#endif
#ifndef HAVE_CV44
	siril_log_message("OpenCV 4.2 used, shift-only registration transformation unavailable\n");
#endif
#ifdef HAVE_LIBCURL
	siril_log_message("Built with libcurl\n");
#endif
#ifdef HAVE_EXIV2
	siril_log_message("Built with exiv2\n");
#endif

	siril_log_message("Can read and write FITS files\n");
	siril_log_message("Can read and write SER files\n");
	siril_log_message("Can read and write BMP files\n");
	siril_log_message("Can read and write NetPBM files\n");
#ifdef HAVE_LIBRAW
	siril_log_message("Can read DSLR RAW files\n");
#endif
#ifdef HAVE_LIBJPEG
	siril_log_message("Can read and write JPEG files\n");
#endif
#ifdef HAVE_LIBJXL
	siril_log_message("Can read and write JPEG XL files\n");
#endif
#ifdef HAVE_LIBPNG
	siril_log_message("Can read and write PNG files\n");
#endif
#ifdef HAVE_LIBTIFF
	siril_log_message("Can read and write TIFF and Astro-TIFF files\n");
#endif
#ifdef HAVE_LIBXISF
	siril_log_message("Can read XISF files\n");
#endif
#ifdef HAVE_LIBHEIF
	siril_log_message("Can read and write HEIF and AVIF files\n");
#endif
	siril_log_message("Can read IRIS PIC files\n");
#ifdef HAVE_FFMS2
	siril_log_message("Can read films\n");
#endif
#ifdef HAVE_FFMPEG
	siril_log_message("Can export films\n");
#endif
	siril_log_message("Can export uncompressed AVI\n");
	return CMD_OK;
}

int process_exit(int nb) {
	// This GTK function is OK to call regardless of thread, as it will quit the GTK main loop
	// and hand control back to main.c and the application will terminate.
	gtk_main_quit();
	return CMD_OK;
}

int process_extract(int nb) {
	int Nbr_Plan, maxplan, mins;
	gchar *end;

	Nbr_Plan = g_ascii_strtoull(word[1], &end, 10);

	mins = min (gfit->rx, gfit->ry);
	maxplan = log(mins) / log(2) - 2;

	if (end == word[1] || Nbr_Plan > maxplan){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return CMD_GENERIC_ERROR;
	}

	struct wavelets_filter_data *args = calloc(1, sizeof(struct wavelets_filter_data));

	args->Type = TO_PAVE_BSPLINE;
	args->Nbr_Plan = Nbr_Plan;
	args->fit = gfit;
	siril_log_message(_("Extracting wavelets (%d plans)\n"), Nbr_Plan);
	if (!start_in_new_thread(extract_plans, args)) {
		free(args);
		return CMD_ARG_ERROR;
	}

	return CMD_OK;
}

int process_reloadscripts(int nb){
	if (com.headless) {
		siril_log_color_message(_("Error: cannot reload script menu when running headless\n"), "red");
		return CMD_GENERIC_ERROR;
	} else {
		g_thread_unref(g_thread_new("refresh_scripts", refresh_scripts_in_thread, NULL));
	}
	return CMD_OK;
}

int process_requires(int nb) {
	version_number current_version = get_current_version_number();
	version_number required_version = get_version_number_from_string(word[1]);
	version_number expired_version = { 0 };
	int new_enough = compare_version(current_version, required_version);
	int too_new = -1;
	if (word[2]) {
		expired_version = get_version_number_from_string(word[2]);
		too_new = compare_version(current_version, expired_version);
	}
	if (new_enough < 0) {
		siril_log_color_message(_("A newer version (>= %lu.%lu.%lu-%s%lu) of Siril is required in order to run this script: please update your installation.\n"), "red", required_version.major_version, required_version.minor_version, required_version.micro_version, required_version.beta_version ? _("beta"): required_version.rc_version ?_("RC"): "", required_version.patched_version);
		return CMD_GENERIC_ERROR;
	} else if (too_new >= 0) {
		siril_log_color_message(_("This script has been marked as obsolete for this version of Siril, please check for an update to the script.\n"), "red");
		return CMD_GENERIC_ERROR;
	}
	siril_log_message(_("OK! This script is compatible with this version of Siril.\n"));
	return CMD_OK;
}

int process_boxselect(int nb){
	/* first case: no argument, printing current selection */
	if (nb == 1) {
		if (com.selection.w > 0 && com.selection.h > 0)
			siril_log_message(_("Current selection [x, y, w, h]: %d %d %d %d\n"),
					com.selection.x, com.selection.y,
					com.selection.w, com.selection.h);
		else siril_log_message(_("No current selection in image\n"));
		return CMD_OK;
	}

	/* second case: clearing the selection */
	if (nb > 1 && !strcmp(word[1], "-clear")) {
		if (nb > 2) {
			siril_log_message(_("Too many arguments to boxselect, use either -clear or the coordinates, not both.\n"));
			return CMD_ARG_ERROR;
		}
		delete_selected_area();
		siril_log_message(_("Selected area in image was cleared\n"));
		return CMD_OK;
	}

	/* third case: setting a new selection */
	if (nb != 5) {
		siril_log_message(_("Please specify x, y, w and h, aborting\n"));
		return CMD_ARG_ERROR;
	}

	gboolean parse_error = FALSE;
	int x, y, w, h;
	char *end;
	x = g_ascii_strtoull(word[1], &end, 10);
	if (word[1] == end) parse_error = TRUE;
	y = g_ascii_strtoull(word[2], &end, 10);
	if (word[2] == end) parse_error = TRUE;
	w = g_ascii_strtoull(word[3], &end, 10);
	if (word[3] == end) parse_error = TRUE;
	h = g_ascii_strtoull(word[4], &end, 10);
	if (word[4] == end) parse_error = TRUE;
	if (parse_error || w == 0 || h == 0) {
		siril_log_message(_("Please specify x, y, w and h, aborting\n"));
		return CMD_ARG_ERROR;
	}
	if (x+w > gfit->rx || y+h > gfit->ry) {
		siril_log_message(_("The provided coordinates are outside the dimension of the currently loaded image (%d x %d).\n"), gfit->rx, gfit->ry);
		return CMD_ARG_ERROR;
	}

	if (com.selection.w > 0 && com.selection.h > 0)
		siril_log_message(_("Overriding previous selection\n"));

	com.selection.x = x;
	com.selection.y = y;
	com.selection.w = w;
	com.selection.h = h;
	siril_log_message(_("Current selection [x, y, w, h]: %d %d %d %d\n"), x, y, w, h);
	if (!com.script)
		gui_function(new_selection_zone, NULL);
	return CMD_OK;
}

static void rgb_extract_last_options(int next_arg, gchar **result_filename,
		const gchar *default_filename, gboolean *do_sum) {
	gchar *filename = NULL;
	*do_sum = TRUE;

	for (int i = next_arg; word[i]; i++) {
		if (g_str_has_prefix(word[i], "-out=") && word[i][5] != '\0') {
			gchar* val = word[i] + 5;
			if (filename) g_free(filename);
			if (g_str_has_suffix(val, com.pref.ext)) {
				filename = g_strdup(val);
			} else {
				filename = g_strdup_printf("%s%s", val, com.pref.ext);
			}
		} else if (!g_strcmp0(word[i], "-nosum")) {
			*do_sum = FALSE;
		}
	}
	if (filename == NULL) {
		filename = g_strdup_printf("%s%s", default_filename, com.pref.ext);
	}
	*result_filename = filename;
}

int process_rgbcomp(int nb) {
	fits r = { 0 }, g = { 0 }, b = { 0 };
	fits rgb = { 0 }, *rgbptr = &rgb;
	int retval = 0, next_arg;
	gboolean do_sum;
	gchar *result_filename;

	if (g_str_has_prefix(word[1], "-lum=")) {
		char *lum_file = word[1] + 5;
		if (nb < 3) {
			return CMD_WRONG_N_ARG;
		}

		gboolean had_an_rgb_image = FALSE;
		fits l = { 0 };
		if (readfits(lum_file, &l, NULL, TRUE)) return CMD_INVALID_IMAGE;
		if (readfits(word[2], &r, NULL, TRUE)) { clearfits(&l); return CMD_INVALID_IMAGE; }
		/* colors can be a single color image or 3 mono */
		if (r.naxis == 3) {
			extract_fits(&r, &g, 1, FALSE);
			extract_fits(&r, &b, 2, FALSE);
			keep_only_first_channel(&r);
			had_an_rgb_image = TRUE;
			next_arg = 3;
		} else {
			if (nb < 5) {
				return CMD_WRONG_N_ARG;
			}
			if (readfits(word[3], &g, NULL, TRUE)) {
				clearfits(&l); clearfits(&r);
				return CMD_INVALID_IMAGE;
			}
			if (readfits(word[4], &b, NULL, TRUE)) {
				clearfits(&l); clearfits(&r); clearfits(&g);
				return CMD_INVALID_IMAGE;
			}
			next_arg = 5;
		}

		if (l.naxes[2] != 1 || check_loaded_fits_params(&l, &r, &g, &b, NULL)) {
			clearfits(&l); clearfits(&r); clearfits(&g); clearfits(&b);
			siril_log_message(_("Image must all have the same dimensions and be monochrome\n"));
			return CMD_ARG_ERROR;
		}

		if (new_fit_image(&rgbptr, l.rx, l.ry, 3, DATA_FLOAT)) {
			clearfits(&l); clearfits(&r); clearfits(&g); clearfits(&b);
			PRINT_ALLOC_ERR;
			return CMD_ALLOC_ERROR;
		}

		/* we need to parse last parameters before merge_fits_headers_to_result */
		rgb_extract_last_options(next_arg, &result_filename, "composed_lrgb", &do_sum);

		if (had_an_rgb_image)
			merge_fits_headers_to_result(rgbptr, do_sum, &l, &r, NULL);
		else merge_fits_headers_to_result(rgbptr, do_sum, &l, &r, &g, &b, NULL);
		rgbptr->history = g_slist_append(rgbptr->history, strdup("LRGB composition"));

		size_t nbpix = l.naxes[0] * l.naxes[1];
		for (size_t i = 0; i < nbpix; i++) {
			gdouble h, s, el, rd, gd, bd;
			rgb_to_hsl(r.fdata[i], g.fdata[i], b.fdata[i], &h, &s, &el);
			hsl_to_rgb(h, s, l.fdata[i], &rd, &gd, &bd);
			rgb.fpdata[RLAYER][i] = (float)rd;
			rgb.fpdata[GLAYER][i] = (float)gd;
			rgb.fpdata[BLAYER][i] = (float)bd;
		}
		clearfits(&l);
	} else {
		if (readfits(word[1], &r, NULL, TRUE)) return CMD_INVALID_IMAGE;
		if (readfits(word[2], &g, NULL, TRUE)) { clearfits(&r); return CMD_INVALID_IMAGE; }
		if (readfits(word[3], &b, NULL, TRUE)) { clearfits(&r); clearfits(&g); return CMD_INVALID_IMAGE; }
		if (r.naxes[2] != 1 || check_loaded_fits_params(&r, &g, &b, NULL)) {
			clearfits(&r); clearfits(&g); clearfits(&b);
			siril_log_message(_("Image must all have the same dimensions and be monochrome\n"));
			return CMD_ARG_ERROR;
		}

		if (new_fit_image(&rgbptr, r.rx, r.ry, 3, DATA_FLOAT)) {
			clearfits(&r); clearfits(&g); clearfits(&b);
			PRINT_ALLOC_ERR;
			return CMD_ALLOC_ERROR;
		}
		next_arg = 4;
		rgb_extract_last_options(next_arg, &result_filename, "composed_rgb", &do_sum);

		merge_fits_headers_to_result(rgbptr, do_sum, &r, &g, &b, NULL);
		rgbptr->history = g_slist_append(rgbptr->history, strdup("RGB composition"));
		size_t nbpix = r.naxes[0] * r.naxes[1];
		for (size_t i = 0; i < nbpix; i++) {
			rgb.fpdata[RLAYER][i] = r.fdata[i];
			rgb.fpdata[GLAYER][i] = g.fdata[i];
			rgb.fpdata[BLAYER][i] = b.fdata[i];
		}
		next_arg = 4;
	}

	clearfits(&r); clearfits(&g); clearfits(&b);

	retval = savefits(result_filename, rgbptr);
	siril_log_message(_("Successful RGB composition with metadata accumulate option %s.\n"), do_sum ? "enable" : "disabled");
	g_free(result_filename);
	clearfits(rgbptr);

	return retval ? CMD_GENERIC_ERROR : CMD_OK;
}

// used for PCC and SPCC commands
static int do_pcc(int nb, gboolean spectro) {
	if (!has_wcs(gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}
	siril_cat_index cat = CAT_AUTO;
	double mag_offset = 0.0, target_mag = -1.0;
	gchar *monosensor = NULL, *oscsensor = NULL, *rfilter = NULL, *gfilter = NULL, *bfilter = NULL, *oscfilter = NULL, *osclpf = NULL, *whiteref = NULL;
	gchar* spcc_strings_to_free[8] = { oscsensor, oscfilter, osclpf, monosensor, rfilter, gfilter, bfilter, whiteref };
	int mono_or_osc = com.pref.spcc.is_mono == 1 ? -1 : 1;
	int dslr; // for SPCC
	gboolean nb_mode = FALSE; // for SPCC
	double wl[3] = { -1.0 , -1.0, -1.0}; // for SPCC
	double bw[3] = { -1.0 , -1.0, -1.0}; // for SPCC
	double t0 = -2.8, t1 = 2.0; // background correction tolerance

	gboolean atmos = FALSE, slp = TRUE;
	double pressure = 1013.25; // standard atmosphere
	double obsheight = gfit->keywords.siteelev != DEFAULT_DOUBLE_VALUE ? gfit->keywords.siteelev : 10.0;
	gboolean local_kstars = local_kstars_available();
	gboolean local_gaia = local_gaia_available();
	int next_arg = 1;

	while (nb > next_arg && word[next_arg]) {
		if (g_str_has_prefix(word[next_arg], "-bgtol=")) {
			char *arg = word[next_arg] + 7;
			gchar *end;
			t0 = g_ascii_strtod(arg, &end);
			gchar *arg2 = end + 1;
			t1 = g_ascii_strtod(arg2, &end);
			if (end == arg || end == arg2 || t0 < 0.1 || t0 > 10.0 || t1 < 0.1 || t1 > 10.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				return CMD_ARG_ERROR;
			}
			t0 = 0.0 - t0;
		} else if (g_str_has_prefix(word[next_arg], "-limitmag=")) {
			char *arg = word[next_arg] + 10;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				return CMD_ARG_ERROR;
			}
			if (arg[0] == '-' || arg[0] == '+')
				mag_offset = value;
			else target_mag = value;
		} else if (g_str_has_prefix(word[next_arg], "-catalog=")) {
			if (!spectro) {
				char *arg = word[next_arg] + 9;
				if (!g_strcmp0(arg, "nomad"))
					cat = local_kstars ? CAT_LOCAL_KSTARS : CAT_NOMAD;
				else if (!g_strcmp0(arg, "gaia"))
					cat = CAT_GAIADR3;
				else if (!g_strcmp0(arg, "localgaia")) {
					cat = local_gaia ? CAT_LOCAL_GAIA_ASTRO : CAT_GAIADR3;
					if (cat == CAT_GAIADR3) {
						siril_log_color_message(_("Local Gaia catalog is unavailable, reverting to online Gaia catalog via Vizier\n"), "salmon");
					}
				}
				else if (!g_strcmp0(arg, "apass"))
					cat = CAT_APASS;
				else {
					siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
					for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
					return CMD_ARG_ERROR;
				}
			} else {
				char *arg = word[next_arg] + 9;
				if (!g_strcmp0(arg, "gaia"))
					cat = CAT_REMOTE_GAIA_XPSAMP;
				else if (!g_strcmp0(arg, "localgaia")) {
					cat = CAT_LOCAL_GAIA_XPSAMP;
					if (!local_gaia_xpsamp_available()) {
						siril_log_color_message(_("Local Gaia catalog is unavailable, reverting to online Gaia catalog via ESA\n"), "salmon");
						cat = CAT_REMOTE_GAIA_XPSAMP;
				
					}
				} else {
					siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
					for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
					return CMD_ARG_ERROR;
				}
			}
		} else if (spectro && !g_strcmp0(word[next_arg], "-narrowband")) {
			nb_mode = TRUE;
		} else if (spectro && g_str_has_prefix(word[next_arg], "-rwl=")) {
			char *arg = word[next_arg] + 5;
			wl[RLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-gwl=")) {
			char *arg = word[next_arg] + 5;
			wl[GLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-bwl=")) {
			char *arg = word[next_arg] + 5;
			wl[BLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-rbw=")) {
			char *arg = word[next_arg] + 5;
			bw[RLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-gbw=")) {
			char *arg = word[next_arg] + 5;
			bw[GLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-bbw=")) {
			char *arg = word[next_arg] + 5;
			bw[BLAYER] = g_ascii_strtod(arg, NULL);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-monosensor=")) {
			char *arg = word[next_arg] + 12;
			if (monosensor) g_free(monosensor);
			monosensor = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-oscsensor=")) {
			char *arg = word[next_arg] + 11;
			if (oscsensor) g_free(oscsensor);
			oscsensor = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-rfilter=")) {
			char *arg = word[next_arg] + 9;
			if (rfilter) g_free(rfilter);
			rfilter = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-gfilter=")) {
			char *arg = word[next_arg] + 9;
			if (gfilter) g_free(gfilter);
			gfilter = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-bfilter=")) {
			char *arg = word[next_arg] + 9;
			if (bfilter) g_free(bfilter);
			bfilter = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-oscfilter=")) {
			char *arg = word[next_arg] + 11;
			if (oscfilter) g_free(oscfilter);
			oscfilter = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-osclpf=")) {
			char *arg = word[next_arg] + 8;
			if (osclpf) g_free(osclpf);
			osclpf = g_strdup(arg);
		} else if (spectro && g_str_has_prefix(word[next_arg], "-whiteref=")) {
			char *arg = word[next_arg] + 10;
			if (whiteref) g_free(whiteref);
			whiteref = g_strdup(arg);
		} else if (spectro && !g_strcmp0(word[next_arg], "-atmos")) {
			atmos = TRUE;
		} else if (spectro && !g_strcmp0(word[next_arg], "-slp")) {
			slp = TRUE;
		} else if (spectro && g_str_has_prefix(word[next_arg], "-obsheight=")) {
			char *arg = word[next_arg] + 11, *end;
			obsheight = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				return CMD_ARG_ERROR;
			}
		} else if (spectro && g_str_has_prefix(word[next_arg], "-pressure=")) {
			char *arg = word[next_arg] + 10, *end;
			pressure = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				return CMD_ARG_ERROR;
			}
		} else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
			for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
			return CMD_ARG_ERROR;
		}
		next_arg++;
	}

	if (gfit->keywords.wcslib->lin.dispre == NULL) {
		siril_log_message(_("Found linear plate solve data, you may need to solve your image with distortions to ensure correct calibration of stars near image corners.\n"));
	}

	if (spectro && nb_mode) {
		wl[RLAYER] = wl[RLAYER] < 0.0 ? com.pref.spcc.red_wl : wl[RLAYER];
		wl[GLAYER] = wl[GLAYER] < 0.0 ? com.pref.spcc.green_wl : wl[GLAYER];
		wl[BLAYER] = wl[BLAYER] < 0.0 ? com.pref.spcc.blue_wl : wl[BLAYER];
		bw[RLAYER] = bw[RLAYER] < 0.0 ? com.pref.spcc.red_bw : bw[RLAYER];
		bw[GLAYER] = bw[GLAYER] < 0.0 ? com.pref.spcc.green_bw : bw[GLAYER];
		bw[BLAYER] = bw[BLAYER] < 0.0 ? com.pref.spcc.blue_bw : bw[BLAYER];
		if (wl[RLAYER] < 380.0 || wl[RLAYER] > 700.0 || wl[GLAYER] < 380.0 || wl[GLAYER] > 700.0 || wl[BLAYER] < 380.0 || wl[BLAYER] > 700.0) {
		siril_log_message(_("NB wavelength out of range (must be 380 <= wl <= 700)\n"));
			for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
			return CMD_ARG_ERROR;
		}
		if (bw[RLAYER] < 1.0 || bw[RLAYER] > 40.0 || bw[GLAYER] < 1.0 || bw[GLAYER] > 40.0 || bw[BLAYER] < 1.0 || bw[BLAYER] > 40.0) {
			siril_log_message(_("NB bandwidth out of range (must be 1.0 <= wl <= 40.0)\n"));
			for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
			return CMD_ARG_ERROR;
		}
	}

	if (!spectro && cat == CAT_AUTO) {
		cat = local_kstars ? CAT_LOCAL_KSTARS : CAT_NOMAD;
	} else if (spectro && cat == CAT_AUTO) {
		cat = local_gaia_xpsamp_available() ? CAT_LOCAL_GAIA_XPSAMP : CAT_REMOTE_GAIA_XPSAMP;
	}
	if (!spectro && local_kstars && cat != CAT_LOCAL_KSTARS) {
		siril_log_color_message(_("Using remote %s instead of local NOMAD catalogue\n"),
				"salmon", catalog_to_str(cat));
		local_kstars = FALSE;
	}

	struct photometric_cc_data *pcc_args = NULL;	// filled only if pcc_command

	pcc_args = calloc(1, sizeof(struct photometric_cc_data));
	pcc_args->fit = gfit;
	pcc_args->bg_auto = TRUE;
	pcc_args->spcc = spectro;
	pcc_args->t0 = t0;
	pcc_args->t1 = t1;
	if (spectro) {
		pcc_args->nb_mode = nb_mode;
		if (nb_mode) {
			memcpy(&pcc_args->nb_center, wl, sizeof(double[3]));
			memcpy(&pcc_args->nb_bandwidth, bw, sizeof(double[3]));
		}
		pcc_args->atmos_corr = atmos;
		pcc_args->atmos_obs_height = obsheight;
		pcc_args->atmos_pressure = pressure;
		pcc_args->atmos_pressure_is_slp = slp;
		if (oscsensor || mono_or_osc == 1) {
			pcc_args->selected_sensor_osc = get_favourite_oscsensor(com.spcc_data.osc_sensors, oscsensor ? oscsensor : com.pref.spcc.oscsensorpref);
			GList *osc = g_list_nth(com.spcc_data.osc_sensors, pcc_args->selected_sensor_osc);
			if (osc) {
				osc_sensor *oscsensor = (osc_sensor*) osc->data;
				dslr = oscsensor->channel[0].is_dslr;
			} else {
				dslr = com.pref.spcc.is_dslr;
			}
			pcc_args->selected_filter_osc = get_favourite_spccobject(com.spcc_data.osc_filters, oscfilter ? oscfilter : com.pref.spcc.oscfilterpref ? com.pref.spcc.oscfilterpref : "No filter");
			pcc_args->selected_filter_lpf = get_favourite_spccobject(com.spcc_data.osc_lpf, osclpf ? osclpf : com.pref.spcc.lpfpref ? com.pref.spcc.lpfpref : "Full spectrum");
			if (pcc_args->selected_sensor_osc == -1 || pcc_args->selected_filter_osc == -1 || (pcc_args->is_dslr && pcc_args->selected_filter_lpf == -1)) {
				siril_log_message(_("Either the sensor or a filter / LPF was not specified as argument or guessable from previous use. Ensure all necessary data is set.\n"));
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				free(pcc_args);
				return CMD_ARG_ERROR;
			}
			pcc_args->spcc_mono_sensor = FALSE;
			if (dslr)
				siril_log_message(_("SPCC will use OSC sensor \"%s\", filter \"%s\" and DSLR LPF \"%s\"\n"), oscsensor ? oscsensor : com.pref.spcc.oscsensorpref, oscfilter ? oscfilter : com.pref.spcc.oscfilterpref, osclpf ? osclpf : com.pref.spcc.lpfpref);
			else
				siril_log_message(_("SPCC will use OSC sensor \"%s\" and filter \"%s\"\n"), oscsensor ? oscsensor : com.pref.spcc.oscsensorpref, oscfilter ? oscfilter : com.pref.spcc.oscfilterpref);
			pcc_args->is_dslr = dslr == 1 ? TRUE : dslr == -1 ? FALSE : com.pref.spcc.is_dslr;
		} else {
			pcc_args->selected_sensor_m = get_favourite_spccobject(com.spcc_data.mono_sensors, monosensor ? monosensor : com.pref.spcc.monosensorpref);
			pcc_args->selected_filter_r = get_favourite_spccobject(com.spcc_data.mono_filters[0], rfilter ? rfilter : com.pref.spcc.redpref);
			pcc_args->selected_filter_g = get_favourite_spccobject(com.spcc_data.mono_filters[1], gfilter ? gfilter : com.pref.spcc.greenpref);
			pcc_args->selected_filter_b = get_favourite_spccobject(com.spcc_data.mono_filters[2], bfilter ? bfilter : com.pref.spcc.bluepref);
			if (pcc_args->selected_sensor_m == -1 || pcc_args->selected_filter_r == -1 || pcc_args->selected_filter_g == -1 || pcc_args->selected_filter_b == -1 ) {
				siril_log_message(_("Either the sensor or a filter was not specified as argument or guessable from previous use. Ensure all necessary data is set.\n"));
				for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }
				free(pcc_args);
				return CMD_ARG_ERROR;
			}
			pcc_args->spcc_mono_sensor = TRUE;
			siril_log_message(_("SPCC will use mono senor \"%s\" and filters \"%s\", \"%s\" and \"%s\n"), monosensor ? monosensor : com.pref.spcc.monosensorpref, rfilter ? rfilter : com.pref.spcc.redpref, gfilter ? gfilter : com.pref.spcc.greenpref, bfilter ? bfilter : com.pref.spcc.bluepref);
		}
		pcc_args->selected_white_ref = get_favourite_spccobject(com.spcc_data.wb_ref, whiteref ? whiteref : "Average Spiral Galaxy");
		pcc_args->do_plot = FALSE; // TODO: do we want to enable if from GUI?
	}
	for (int z = 0 ; z < 8 ; z++) { g_free(spcc_strings_to_free[z]); }

	if (target_mag > -1.0) {
		pcc_args->mag_mode = LIMIT_MAG_ABSOLUTE;
		pcc_args->magnitude_arg = target_mag;
	} else if (mag_offset != 0.0) {
		pcc_args->mag_mode = LIMIT_MAG_AUTO_WITH_OFFSET;
		pcc_args->magnitude_arg = mag_offset;
	}else {
		pcc_args->mag_mode = LIMIT_MAG_AUTO;
	}

	pcc_args->catalog = cat;
	if (spectro)
		load_spcc_metadata_if_needed();

	if (!start_in_new_thread(photometric_cc_standalone, pcc_args)) {
		g_free(pcc_args->datalink_path);
		free(pcc_args);
		return CMD_GENERIC_ERROR;
	}

	return CMD_OK;
}

int process_pcc(int nb) {
	return do_pcc(nb, FALSE);
}

int process_spcc(int nb) {
	siril_check_spcc_mirrors(TRUE, TRUE);
	return do_pcc(nb, TRUE);
}

// used for platesolve and seqplatesolve commands
int process_platesolve(int nb) {
	gboolean noflip = FALSE, force = FALSE, downsample = FALSE, autocrop = TRUE, nocache = FALSE,
		asnet_blind_pos = FALSE, asnet_blind_res = FALSE, noreg = FALSE;
	gboolean forced_metadata[3] = { 0 }; // used for sequences in the image hook, for center, pixel and focal
	SirilWorldCS *target_coords = NULL;
	double forced_focal = -1.0, forced_pixsize = -1.0;
	double mag_offset = 0.0, target_mag = -1.0, searchradius = com.pref.astrometry.radius_degrees;
	int order = com.pref.astrometry.sip_correction_order; // we default to the pref value
	siril_cat_index cat = CAT_AUTO;
	gboolean seqps = word[0][0] == 's';
	sequence *seq = NULL;
	platesolve_solver solver = SOLVER_SIRIL;
	gchar *distofilename = NULL;
	cmd_errors retval = CMD_OK;
	struct astrometry_data *args = NULL;

	int next_arg = 1;
	if (seqps) {
		if (!(seq = load_sequence(word[1], NULL)))
			return CMD_SEQUENCE_NOT_FOUND;
		if (seq->type == SEQ_SER) {
			siril_log_color_message(_("SER cannot contain WCS info, plate solving will export to FITS cube format\n"), "salmon");
		}
		next_arg++;
		if (check_seq_is_comseq(seq)) {
			free_sequence(seq, TRUE);
			seq = &com.seq;
		}
	} else {
		force = !has_wcs(gfit);
	}

	// check if we have target_coords
	if (nb > next_arg && (word[next_arg][0] != '-' || (word[next_arg][1] >= '0' && word[next_arg][1] <= '9'))) {
		char *sep = strchr(word[next_arg], ',');
		if (!sep) {
			if (nb == 2) {
				siril_log_message(_("Could not parse target coordinates\n"));
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			target_coords = siril_world_cs_new_from_objct_ra_dec(word[next_arg], word[next_arg+1]);
			next_arg += 2;
		}
		else {
			*sep++ = '\0';
			target_coords = siril_world_cs_new_from_objct_ra_dec(word[next_arg], sep);
			next_arg++;
		}
		if (!target_coords) {
			siril_log_message(_("Could not parse target coordinates\n"));
			retval = CMD_ARG_ERROR;
			goto clean_and_exit_platesolve;
		}
		forced_metadata[FORCED_CENTER] = TRUE;
	}

	while (nb > next_arg && word[next_arg]) {
		if (!strcmp(word[next_arg], "-noflip"))
			noflip = TRUE;
		else if (!strcmp(word[next_arg], "-nocrop"))
			autocrop = FALSE;
		else if (!strcmp(word[next_arg], "-force"))
			force = TRUE;
		else if (!strcmp(word[next_arg], "-downscale"))
			downsample = TRUE;
		else if (!strcmp(word[next_arg], "-nocache"))
			nocache = TRUE;
		else if (!strcmp(word[next_arg], "-noreg"))
			noreg = TRUE;
		else if (!strcmp(word[next_arg], "-blindpos"))
			asnet_blind_pos = TRUE;
		else if (!strcmp(word[next_arg], "-blindres"))
			asnet_blind_res = TRUE;
		else if (g_str_has_prefix(word[next_arg], "-focal=")) {
			char *arg = word[next_arg] + 7;
			gchar *end;
			forced_focal = g_ascii_strtod(arg, &end);
			if (end == arg || forced_focal <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			forced_metadata[FORCED_FOCAL] = TRUE;
		}
		else if (g_str_has_prefix(word[next_arg], "-pixelsize=")) {
			char *arg = word[next_arg] + 11;
			gchar *end;
			forced_pixsize = g_ascii_strtod(arg, &end);
			if (end == arg || forced_pixsize <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			forced_metadata[FORCED_PIXEL] = TRUE;
		}
		else if (g_str_has_prefix(word[next_arg], "-limitmag=")) {
			char *arg = word[next_arg] + 10;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			if (arg[0] == '-' || arg[0] == '+')
				mag_offset = value;
			else target_mag = value;
		}
		else if (g_str_has_prefix(word[next_arg], "-order=")) {
			char *arg = word[next_arg] + 7;
			gchar *end;
			int value = g_ascii_strtoull(arg, &end, 10);
			if (end == arg || value < 1 || value > 5) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			order = value;
		}
		else if (g_str_has_prefix(word[next_arg], "-radius=")) {
			char *arg = word[next_arg] + 8;
			gchar *end;
			searchradius = g_ascii_strtod(arg, &end);
			if (end == arg || searchradius < 0.0 || searchradius > 30.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
		}
		else if (g_str_has_prefix(word[next_arg], "-catalog=")) {
			char *arg = word[next_arg] + 9;
			if (!g_strcmp0(arg, "tycho2"))
				cat = CAT_TYCHO2;
			else if (!g_strcmp0(arg, "nomad"))
				cat = CAT_NOMAD;
			else if (!g_strcmp0(arg, "gaia"))
				cat = CAT_GAIADR3;
			else if (!g_strcmp0(arg, "localgaia")) {
				if (local_gaia_available())
					cat = CAT_LOCAL_GAIA_ASTRO;
				else {
					cat = CAT_GAIADR3;
					siril_log_color_message(_("Local Gaia catalog is unavailable, reverting to online Gaia catalog via Vizier\n"), "salmon");
				}
			}
			else if (!g_strcmp0(arg, "ppmxl"))
				cat = CAT_PPMXL;
			else if (!g_strcmp0(arg, "bsc"))
				cat = CAT_BSC;
			else if (!g_strcmp0(arg, "apass"))
				cat = CAT_APASS;
			else {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
		}
		else if (g_str_has_prefix(word[next_arg], "-disto=")) {
			gchar *arg  = word[next_arg] + 7;
			if (arg[0] == '\0') {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				retval = CMD_ARG_ERROR;
				goto clean_and_exit_platesolve;
			}
			if (distofilename)
				g_free(distofilename); // in case -disto= is passed twice, don't leak distofilename
			distofilename = g_strdup(arg);
		}
		else if (!g_ascii_strcasecmp(word[next_arg], "-localasnet")) {
			if (cat != CAT_AUTO)
				siril_log_message(_("Specifying a catalog has no effect for astrometry.net solving\n"));
			solver = SOLVER_LOCALASNET;
		} else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
			retval = CMD_ARG_ERROR;
			goto clean_and_exit_platesolve;
		}
		next_arg++;
	}

	if (!seqps && !force) {
		siril_log_message(_("Image is already plate solved. Nothing will be done.\n"));
		goto clean_and_exit_platesolve; // not an arror, retval is CMD_OK
	}

	if (solver == SOLVER_LOCALASNET && !asnet_is_available()) {
		siril_log_color_message(_("The local astrometry.net solver was not found, aborting. Please check the settings.\n"), "red");
		retval = CMD_GENERIC_ERROR;
		goto clean_and_exit_platesolve;
	}

	fits reffit = { 0 };
	fits *preffit = &reffit;
	if (seqps) {
		int image_to_load = sequence_find_refimage(seq);
		seq->reference_image = sequence_find_refimage(seq);
		if (seq_read_frame_metadata(seq, image_to_load, preffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			retval = CMD_SEQUENCE_NOT_FOUND;
			goto clean_and_exit_platesolve;
		}
	} else
		preffit = gfit;

	if (solver == SOLVER_SIRIL) {
		if (asnet_blind_pos) {
			siril_log_color_message(_("Siril internal solver cannot be set blind in %s\n"), "salmon", _("position"));
			asnet_blind_pos = FALSE;
		}
		if (asnet_blind_res) {
			siril_log_color_message(_("Siril internal solver cannot be set blind in %s\n"), "salmon", _("resolution"));
			asnet_blind_res = FALSE;
		}
	}

	if (!target_coords) {
		target_coords = get_eqs_from_header(preffit);
		if (solver != SOLVER_LOCALASNET && !target_coords) {
			siril_log_color_message(_("Cannot plate solve, no target coordinates passed and image header doesn't contain any either\n"), "red");
			if (seqps)
				clearfits(preffit);
			retval = CMD_INVALID_IMAGE;
			goto clean_and_exit_platesolve;
		}
		if (target_coords) {
			siril_log_message(_("Using target coordinate from image header: %f, %f\n"),
			siril_world_cs_get_alpha(target_coords),
			siril_world_cs_get_delta(target_coords));
		} else {
			asnet_blind_pos = TRUE;
		}
	}
	if (target_coords && asnet_blind_pos) {
		if (forced_metadata[FORCED_CENTER]) {
			siril_log_color_message(_("%s is ignored when using local astrometry.net blind in %s\n"), "salmon", _("Center position"), _("position"));
			forced_metadata[FORCED_CENTER] = FALSE;
		}
		siril_world_cs_unref(target_coords);
		target_coords = NULL;
	}
	gboolean iscfa = preffit->keywords.bayer_pattern[0] != '\0';

	// we are now ready to fill the structure
	args = calloc(1, sizeof(struct astrometry_data));

	if (asnet_blind_res) {
		if (forced_metadata[FORCED_PIXEL]) {
			siril_log_color_message(_("%s is ignored when using local astrometry.net blind in %s\n"), "salmon", _("Pixel size"), _("resolution"));
			forced_metadata[FORCED_PIXEL] = FALSE;
		}
		args->pixel_size = 0.;
	} else if (forced_pixsize > 0.0) {
		args->pixel_size = forced_pixsize;
		siril_log_message(_("Using provided pixel size: %.2f\n"), args->pixel_size);
	} else {
		args->pixel_size = max(preffit->keywords.pixel_size_x, preffit->keywords.pixel_size_y);
		if (args->pixel_size <= 0.0) {
			args->pixel_size = com.pref.starfinder_conf.pixel_size_x;
			if (args->pixel_size <= 0.0) {
				siril_log_color_message(_("Pixel size not found in image or in settings, cannot proceed\n"), "red");
				if (seqps)
					clearfits(preffit);
				retval = CMD_INVALID_IMAGE;
				goto clean_and_exit_platesolve;
			}
			siril_log_message(_("Using pixel size from preferences: %.2f\n"), args->pixel_size);
		}
		else siril_log_message(_("Using pixel size from image: %.2f\n"), args->pixel_size);
	}

	if (asnet_blind_res) {
		if (forced_metadata[FORCED_FOCAL]) {
			siril_log_color_message(_("%s is ignored when using local astrometry.net blind in %s\n"), "salmon", _("Focal length"), _("resolution"));
			forced_metadata[FORCED_FOCAL] = FALSE;
		}
		args->focal_length = 0.;
	} else if (forced_focal > 0.0) {
		args->focal_length = forced_focal;
		siril_log_message(_("Using provided focal length: %.2f\n"), args->focal_length);
	} else {
		args->focal_length = preffit->keywords.focal_length;
		if (args->focal_length <= 0.0) {
			args->focal_length = com.pref.starfinder_conf.focal_length;
			if (args->focal_length <= 0.0) {
				siril_log_color_message(_("Focal length not found in image or in settings, cannot proceed\n"), "red");
				if (seqps)
					clearfits(preffit);
				retval = CMD_INVALID_IMAGE;
				goto clean_and_exit_platesolve;
			}
			siril_log_message(_("Using focal length from preferences: %.2f\n"), args->focal_length);
		}
		else siril_log_message(_("Using focal length from image: %.2f\n"), args->focal_length);
	}

	if (target_mag > -1.0) {
		if (solver != SOLVER_SIRIL)
			siril_log_message(_("Magnitude alteration arguments are useless for astrometry.net plate solving\n"));
		else {
			args->mag_mode = LIMIT_MAG_ABSOLUTE;
			args->magnitude_arg = target_mag;
		}
	} else if (mag_offset != 0.0) {
		if (solver != SOLVER_SIRIL)
			siril_log_message(_("Magnitude alteration arguments are useless for astrometry.net plate solving\n"));
		else {
			args->mag_mode = LIMIT_MAG_AUTO_WITH_OFFSET;
			args->magnitude_arg = mag_offset;
		}
	} else {
		args->mag_mode = LIMIT_MAG_AUTO;
	}

	if (seqps)
		clearfits(preffit);
	args->solver = solver;
	args->downsample = downsample;
	args->autocrop = autocrop && solver == SOLVER_SIRIL; // we don't crop fov when using asnet
	args->nocache = nocache || 
					solver == SOLVER_LOCALASNET || 
					((cat == CAT_AUTO || cat == CAT_GAIADR3 || cat == CAT_NOMAD || cat == CAT_TYCHO2) && local_catalogues_available());
	if (!searchradius && solver == SOLVER_LOCALASNET && !asnet_blind_pos) {
		args->searchradius = com.pref.astrometry.radius_degrees;
		siril_log_color_message(_("Cannot force null radius for localasnet if not blind solving, using default instead\n"), "red");
	} else {
		args->searchradius = searchradius;
	}

	if (!com.script && com.selection.w > 0 && com.selection.h > 0 && (!seqps || !seq->is_variable)) {
		args->solvearea = com.selection;
		args->autocrop = FALSE;
	}

	if (distofilename) {
		args->distofilename = distofilename;
	}
	args->force = force;
	memcpy(&args->forced_metadata, forced_metadata, 3 * sizeof(gboolean));
	if (seqps || sequence_is_loaded()) { // we are platesolving an image from a sequence or a sequence, we can't allow to flip (may be registered)
		noflip = TRUE;
		siril_debug_print("forced no flip for solving an image from a sequence\n");
	}
	if (iscfa) { // prevent flipping for bayered images
		noflip = TRUE;
		siril_debug_print("forced no flip for CFA image\n");
	}
	args->flip_image = !noflip;
	args->manual = FALSE;
	args->trans_order = order;
	if (target_coords) {
		args->cat_center = target_coords;
	}
	// catalog query parameters
	if (solver == SOLVER_SIRIL) {
		args->ref_stars = siril_catalog_new(cat);
		args->ref_stars->phot = FALSE;
		args->ref_stars->center_ra = siril_world_cs_get_alpha(target_coords);
		args->ref_stars->center_dec = siril_world_cs_get_delta(target_coords);
	} else {
		args->asnet_blind_pos = asnet_blind_pos;
		args->asnet_blind_res = asnet_blind_res;
	}

	args->for_sequence = seqps;
	args->verbose = !seqps;
	args->asnet_checked = TRUE;

	// sequence
	if (seqps) {
		args->update_reg = !noreg && !(seq->type == SEQ_FITSEQ || seq->type == SEQ_SER);
		args->sfargs = calloc(1, sizeof(struct starfinder_data));
		args->sfargs->im.from_seq = seq;
		args->sfargs->layer = -1;
		args->sfargs->keep_stars = TRUE;
		args->sfargs->save_to_file = TRUE;
		args->sfargs->max_stars_fitted = BRIGHTEST_STARS;
		start_sequence_astrometry(seq, args);
		return CMD_OK;
	}
	// single-image
	if (solver == SOLVER_LOCALASNET) { // see !432
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			args->filename = g_strdup(com.uniq->filename);
		} else if (sequence_is_loaded()) {
			args->filename = g_strdup_printf("%s%.5d", com.seq.seqname, com.seq.current + 1);
		} else {
			args->filename = g_strdup_printf("image");
		}
	}
	args->fit = gfit;
	args->numthreads = com.max_thread;
	process_plate_solver_input(args);
	if (!start_in_new_thread(plate_solver, args)) {
		retval = CMD_GENERIC_ERROR;
		goto clean_and_exit_platesolve;
	}
	return CMD_OK;
clean_and_exit_platesolve:
	if (seq && !check_seq_is_comseq(seq))
		free_sequence(seq, TRUE);
	if (target_coords)
		siril_world_cs_unref(target_coords);
	if (args)
		free_astrometry_data(args);
	if (distofilename)
		g_free(distofilename);
	return retval;

}

static conesearch_params* parse_conesearch_args(int nb) {
	conesearch_params *params = init_conesearch_params();
	gboolean local_kstars = local_kstars_available();

	if (!has_wcs(gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		g_free(params->obscode);
		g_free(params);
		return NULL;
	}

	int arg_idx = 1;
	while (arg_idx < nb) {
		if (g_str_has_prefix(word[arg_idx], "-cat=")) {
			char *arg = word[arg_idx] + 5;
			if (!g_strcmp0(arg, "tycho2"))
				params->cat = CAT_TYCHO2;
			else if (!g_strcmp0(arg, "nomad"))
				params->cat = CAT_NOMAD;
			else if (!g_strcmp0(arg, "gaia"))
				params->cat = CAT_GAIADR3;
			else if (!g_strcmp0(arg, "localgaia"))
				params->cat = CAT_LOCAL_GAIA_ASTRO;
			else if (!g_strcmp0(arg, "ppmxl"))
				params->cat = CAT_PPMXL;
			else if (!g_strcmp0(arg, "bsc"))
				params->cat = CAT_BSC;
			else if (!g_strcmp0(arg, "apass"))
				params->cat = CAT_APASS;
			else if (!g_strcmp0(arg, "gcvs"))
				params->cat = CAT_GCVS;
			else if (!g_strcmp0(arg, "vsx"))
				params->cat = CAT_VSX;
			else if (!g_strcmp0(arg, "varisum"))
				params->cat = CAT_VARISUM;
			else if (!g_strcmp0(arg, "simbad"))
				params->cat = CAT_SIMBAD;
			else if (!g_strcmp0(arg, "exo"))
				params->cat = CAT_EXOPLANETARCHIVE;
			else if (!g_strcmp0(arg, "pgc"))
				params->cat = CAT_PGC;
			else if (!g_strcmp0(arg, "aavso_chart"))
				params->cat = CAT_AAVSO_CHART;
			else if (!g_strcmp0(arg, "solsys")) {
				params->cat = CAT_IMCCE;
				if (!gfit->keywords.date_obs) {
					siril_log_color_message(_("This option only works on images that have observation date information\n"), "red");
					g_free(params->obscode);
					g_free(params);
					return NULL;
				}
			} else {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
		} else if (g_str_has_prefix(word[arg_idx], "-obscode=")) {
			char *arg = word[arg_idx] + 9;
			if (strlen(arg) != 3) {
				siril_log_color_message(_("The observatory should be coded as a 3-letter word\n"), "red");
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
			if (params->obscode)
				g_free(params->obscode);
			params->default_obscode_used = FALSE;
			params->obscode = g_strdup(arg);
		} else if (g_str_has_prefix(word[arg_idx], "-trix=")) {
			if (!local_kstars) {
				siril_log_color_message(_("No local catalogues found, ignoring -trix option\n"), "red");
				continue;
			}
			gchar *end;
			int trix = (int) g_ascii_strtoull(word[arg_idx] + 6, &end, 10);
			if (trix < 0 || trix > 511) {
				siril_log_color_message(_("Trixel number must be between 0 and 511\n"), "red");
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
			params->trixel = trix;
		} else if (!strcmp(word[arg_idx], "-phot")) {
			params->photometric = TRUE;
		} else if (g_str_has_prefix(word[arg_idx], "-log=")) {
			char *arg = word[arg_idx] + 5;
			if (!(g_ascii_strcasecmp(arg, "on")))
				params->display_log = BOOL_TRUE;
			else if (!(g_ascii_strcasecmp(arg, "off")))
				params->display_log = BOOL_FALSE;
			else {
				siril_log_message(_("Wrong parameter values. Log must be set to on or off, aborting.\n"));
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
		} else if (g_str_has_prefix(word[arg_idx], "-tag=")) {
			char *arg = word[arg_idx] + 5;
			if (!(g_ascii_strcasecmp(arg, "on")))
				params->display_tag = BOOL_TRUE;
			else if (!(g_ascii_strcasecmp(arg, "off")))
				params->display_tag = BOOL_FALSE;
			else {
				siril_log_message(_("Wrong parameter values. Tag must be set to on or off, aborting.\n"));
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
		} else if (g_str_has_prefix(word[arg_idx], "-out=")) {
			char *arg = word[arg_idx] + 5;
			if (arg[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[arg_idx]);
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
			params->outfilename = g_strdup(arg);
		} else if (!(g_ascii_strcasecmp(word[arg_idx], "-compare"))) {
			params->compare = TRUE;
		} else {
			gchar *end;
			params->limit_mag = g_ascii_strtod(word[arg_idx], &end);
			if (end == word[arg_idx]) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[arg_idx]);
				g_free(params->obscode);
				g_free(params);
				return NULL;
			}
		}
		arg_idx++;
	}

	if (params->cat == CAT_AUTO) {
		params->cat = (local_kstars) ? CAT_LOCAL_KSTARS : CAT_NOMAD;
		if (params->trixel >= 0 && params->cat == CAT_LOCAL_KSTARS)
			params->cat = CAT_LOCAL_TRIX;
	}

	if (params->compare && !is_star_catalogue(params->cat)) {
		siril_log_message(_("Cannot use -compare argument with non-star catalogues, ignoring.\n"));
		params->compare = FALSE;
	}

	return params;
}

int process_conesearch(int nb) {
	conesearch_params *params = parse_conesearch_args(nb);
	if (!params) {
		return CMD_ARG_ERROR;
	}

	// Set destructor (assuming free_conesearch_params exists)
	params->destroy_fn = (destructor)free_conesearch_params;

	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_conesearch_params(params);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = conesearch_image_hook;
	args->description = _("Cone search");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;
	args->command = TRUE;
	args->user = params;
	args->log_hook = NULL;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_catsearch(int nb) {
	if (!has_wcs(gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}

	sky_object_query_args *query_args = init_sky_object_query();
	query_args->fit = gfit;
	if (nb > 2) {
		query_args->name = build_string_from_words(word + 1);
	} else {
		query_args->name = g_strdup(word[1]);
	}
	// Allocate and initialize generic_img_args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_sky_object_query(query_args);
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = catsearch_image_hook;
	args->description = _("Catalog search");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = query_args;
	args->log_hook = catsearch_log_hook;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_findcompstars(int nb) {
	// findcompstars star_name [-narrow] [-catalog={nomad|apass}] [-dvmag=3] [-dbv=0.5] [-emag=0.03] [-out=nina_file.csv]
	if (!has_wcs(gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}
	const char *target = word[1];
	gboolean narrow = FALSE;
	siril_cat_index used_cat = CAT_APASS;
	double delta_Vmag = 3.0, delta_BV = 0.5, emag = 0.03;
	const char *nina_file = NULL;

	int arg_idx = 2;
	while (arg_idx < nb) {
		if (!g_strcmp0(word[arg_idx], "-narrow"))
			narrow = TRUE;
		else if (g_str_has_prefix(word[arg_idx], "-catalog=")) {
			const char *cat = word[arg_idx] + 9;
			if (!g_ascii_strcasecmp(cat, "nomad"))
				used_cat = CAT_NOMAD;
			else if (!g_ascii_strcasecmp(cat, "apass"))
				used_cat = CAT_APASS;
			else {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (!g_ascii_strncasecmp(word[arg_idx], "-dvmag=", 7)) {
			const char *val = word[arg_idx] + 7;
			gchar *end;
			delta_Vmag = g_ascii_strtod(val, &end);
			if (end == val || delta_Vmag < 0.0 || delta_Vmag > 6.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (!g_ascii_strncasecmp(word[arg_idx], "-dbv=", 5)) {
			const char *val = word[arg_idx] + 5;
			gchar *end;
			delta_BV = g_ascii_strtod(val, &end);
			if (end == val || delta_BV < 0.0 || delta_BV > 0.7) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (!g_ascii_strncasecmp(word[arg_idx], "-emag=", 6)) {
			const char *val = word[arg_idx] + 6;
			gchar *end;
			emag = g_ascii_strtod(val, &end);
			if (end == val || emag < 0.0 || emag > 0.1) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[arg_idx], "-out=")) {
			nina_file = word[arg_idx] + 5;
			if (nina_file[0] == '\0') {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_message(_("Unknown option provided: %s\n"), word[arg_idx]);
			return CMD_ARG_ERROR;
		}
		arg_idx++;
	}

	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));
	args->fit = gfit;
	args->target_name = g_strdup(target);
	args->narrow_fov = narrow;
	args->cat = used_cat;
	args->delta_Vmag = delta_Vmag;
	args->delta_BV = delta_BV;
	args->max_emag = emag;
	args->nina_file = g_strdup(nina_file);

	if (!start_in_new_thread(compstars_worker, args)) {
		g_free(args->target_name);
		g_free(args->nina_file);
		free(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int process_start_ls(int nb) {
	// start_ls [-dark=filename] [-flat=filename] [-rotate] [-32bits] [-gradient_removal] [-watch_files]"
	gchar *dark_file = NULL, *flat_file = NULL;
	gboolean use_file_watcher = FALSE/*, remove_gradient = FALSE*/, shift_reg = TRUE, use_32b = FALSE;
	for (int i = 1; i < nb; i++) {
		if (!word[i])
			continue;
		if (g_str_has_prefix(word[i], "-dark="))
			dark_file = g_shell_unquote(word[i] + 6, NULL);
		else if (g_str_has_prefix(word[i], "-flat="))
			flat_file = g_shell_unquote(word[i] + 6, NULL);
		else if (!strcmp(word[i], "-gradient_removal")) {
			//remove_gradient = TRUE;
			siril_log_message("gradient removal in live stacking is not yet implemented\n");
			return CMD_GENERIC_ERROR;
		} else if (!strcmp(word[i], "-watch_files")) {
			//use_file_watcher = TRUE;
			siril_log_message("file watcher in headless live stacking is not yet implemented\n");
			return CMD_GENERIC_ERROR;
		} else if (!strcmp(word[i], "-rotate")) {
			shift_reg = FALSE;
		} else if (!strcmp(word[i], "-32bits")) {
			use_32b = TRUE;
		} else {
			siril_log_message(_("Unknown option provided: %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	return start_livestack_from_command(dark_file, flat_file, use_file_watcher/*, remove_gradient*/, shift_reg, use_32b);
}

int process_livestack(int nb) {
	// livestack filename [-out=result]
	if (!livestacking_is_started()) {
		siril_log_message(_("Live stacking needs to be initialized with the START_LS command first\n"));
		return CMD_NEED_INIT_FIRST;
	}
	if (livestacking_uses_filewatcher()) {
		siril_log_message(_("Live stacking was configured to use file watching, not this command\n"));
		return CMD_ARG_ERROR;
	}

	image_type type;
	char *filename;
	if (stat_file(word[1], &type, &filename)) {
		siril_log_message(_("Could not open file: %s\n"), word[1]);
		return CMD_INVALID_IMAGE;
	}

	livestacking_queue_file(filename);
	return CMD_NO_WAIT;
}

int process_stop_ls(int nb) {
	if (!livestacking_is_started()) {
		siril_log_message(_("Live stacking needs to be initialized with the START_LS command first\n"));
		return CMD_NEED_INIT_FIRST;
	}
	stop_live_stacking_engine();
	return CMD_OK;
}

int process_parse(int nb) {
	pathparse_mode mode = PATHPARSE_MODE_WRITE;
	if (nb == 3) {
		if (g_strcmp0(word[2], "-r")) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[2]);
			return CMD_ARG_ERROR;
		}
		mode = PATHPARSE_MODE_READ;
	} else if (nb > 3) {
		return CMD_WRONG_N_ARG;
	}
	int status;
	gchar *expression = NULL;
	if (gfit->header) { // fits or astrotiff - do not update the header
		expression = path_parse(gfit, word[1], mode, &status);
	} else {
		expression = update_header_and_parse(gfit, word[1], mode, FALSE, &status);
	}
	siril_log_message(_("String in: %s\n"), word[1]);
	siril_log_message(_("String out: %s\n"), expression);
	g_free(expression);
	return CMD_OK;
}

static show_params* parse_show_args(int nb) {
	show_params *params = g_new0(show_params, 1);
	params->display_log = BOOL_NOT_SET;
	params->display_tag = BOOL_NOT_SET;
	int next_arg = 1;

	if (!g_strcmp0(word[next_arg], "-clear")) {
		params->clear = TRUE;
		next_arg++;
		if (nb == 2) {
			return params;
		}
	}

	//passing a list
	if (g_str_has_prefix(word[next_arg], "-list=")) {
		params->file = g_strdup(word[next_arg] + 6);
		next_arg++;
		while (next_arg < nb) {
			if (!g_ascii_strcasecmp(word[next_arg], "-nolog")) {
				params->display_log = BOOL_FALSE;
			} else if (!g_ascii_strcasecmp(word[next_arg], "-notag")) {
				params->display_tag = BOOL_FALSE;
			} else {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
				g_free(params->file);
				g_free(params);
				return NULL;
			}
			next_arg++;
		}
		return params;
	}

	if (nb > next_arg && !isalpha(word[next_arg][0])
			&& (isdigit(word[next_arg][0]) || isdigit(word[next_arg][1]))) {
		goto parse_coords;
	} else {
		if (nb > next_arg + 1) {
			params->name = g_strdup(word[next_arg]);
			next_arg++;
			goto parse_coords;
		} else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
			g_free(params);
			return NULL;
		}
	}

parse_coords:
	if (nb > next_arg && !isalpha(word[next_arg][0]) && (isdigit(word[next_arg][0]) || isdigit(word[next_arg][1]))) {
		char *sep = strchr(word[next_arg], ',');
		if (!sep) {
			params->coords = siril_world_cs_new_from_objct_ra_dec(word[next_arg], word[next_arg + 1]);
			next_arg += 2;
		} else {
			*sep++ = '\0';
			params->coords = siril_world_cs_new_from_objct_ra_dec(word[next_arg], sep);
			next_arg++;
		}
		if (!params->coords) {
			siril_log_message(_("Could not parse target coordinates\n"));
			g_free(params->name);
			g_free(params);
			return NULL;
		}
		if (nb > next_arg) {
			g_free(params->name);
			params->name = g_strdup(word[next_arg]);
		} else {
			g_free(params->name);
			params->name = g_strdup("object");
		}
	}

	return params;
}

int process_show(int nb) {
	show_params *params = parse_show_args(nb);
	if (!params) {
		return CMD_ARG_ERROR;
	}
	int result = execute_show_command(params);
	g_free(params->name);
	g_free(params->file);
	g_free(params);
	return result;
}

int read_cut_pair(char *value, point *pair) {
	char *end;
	pair->x = (double) g_ascii_strtoull(value, &end, 10);
	if (end == value)
		return CMD_ARG_ERROR;
	if (*end != ',')
		return CMD_ARG_ERROR;
	end++;
	value = end;
	pair->y = (double) g_ascii_strtoull(value, &end, 10);
	if (end == value)
		return CMD_ARG_ERROR;
	return CMD_OK;
}

cut_struct *parse_cut_args(int nb, sequence *seq, cmd_errors *err) {
	cut_struct *cut_args = calloc(1, sizeof(cut_struct));
	initialize_cut_struct(cut_args);
	int start = (seq) ? 2 : 1;
	*err = CMD_OK;
	if (seq)
		cut_args->seq = seq;
	else
		cut_args->fit = gfit;
	int nb_layers = (cut_args->seq) ? cut_args->seq->nb_layers : cut_args->fit->naxes[2];
	for (int i = start; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(word[i], "-tri")/* || g_str_has_prefix(word[i], "-bgremove")*/) {
			cut_args->tri = TRUE;
		}
		else if (g_str_has_prefix(word[i], "-cfa")) {
			fits reffit = { 0 };
			if (seq) {
				// loading the sequence reference image's metadata to read its bayer pattern
				int image_to_load = sequence_find_refimage(seq);
				if (seq_read_frame_metadata(seq, image_to_load, &reffit)) {
					siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
					*err = CMD_SEQUENCE_NOT_FOUND;
					break;
				}
			} else {
				reffit = *gfit;
			}
			sensor_pattern pattern = get_cfa_pattern_index_from_string(reffit.keywords.bayer_pattern); // we don't need the validated value here because we just want to know if it's CFA or not
			if (reffit.naxes[2] > 1 || pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
				siril_log_color_message(_("Error: CFA mode cannot be used with color images or mono images with no Bayer pattern.\n"), "red");
				*err = CMD_ARG_ERROR;
				break;
			}
			cut_args->cfa = TRUE;
		}
		else if (g_str_has_prefix(word[i], "-savedat")) {
			cut_args->save_dat = TRUE;
		}
		else if (g_str_has_prefix(word[i], "-arcsec")) {
			cut_args->pref_as = TRUE;
		}
		else if (g_str_has_prefix(word[i], "-width=")) {
			arg += 7;
			cut_args->width = g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-spacing=")) {
			arg += 9;
			cut_args->step = g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-layer=")) {
			arg += 7;
			if (g_str_has_prefix(arg, "red"))
				cut_args->vport = 0;
			else if (g_str_has_prefix(arg, "green"))
				cut_args->vport = 1;
			else if (g_str_has_prefix(arg, "blue"))
				cut_args->vport = 2;
			else if (g_str_has_prefix(arg, "lum")) {
				if (nb_layers == 1) {
					cut_args->vport = 0;
					cut_args->mode = CUT_MONO;
				} else {
					cut_args->vport = 3;
					cut_args->mode = CUT_MONO;
				}
			} else if (g_str_has_prefix(arg, "col")) {
				cut_args->vport = 0;
				cut_args->mode = CUT_COLOR;
			} else {
				siril_log_message(_("Incorrect option follows -layer=, aborting.\n"));
				*err = CMD_ARG_ERROR;
				break;
			}
		}
/*		else if (g_str_has_prefix(arg, "-xaxis=wavenum")) {
			cut_args->plot_as_wavenumber = TRUE;
		}
		else if (g_str_has_prefix(arg, "-xaxis=wavelen")) {
			cut_args->plot_as_wavenumber = FALSE;
		}
		else if (g_str_has_prefix(arg, "-wavenumber1=")) {
			arg += 13;
			cut_args->wavenumber1 = g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-wavenumber2=")) {
			arg += 13;
			cut_args->wavenumber2 = g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-wavelength1=")) {
			arg += 13;
			cut_args->wavenumber1 = 10000000. / g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-wavelength2=")) {
			arg += 13;
			cut_args->wavenumber2 = 10000000. / g_ascii_strtod(arg, &end);
		}
		else if (g_str_has_prefix(arg, "-bgpoly=")) {
			arg += 8;
			cut_args->bg_poly_order = (int) g_ascii_strtod(arg, &end);
		}
*/		else if (g_str_has_prefix(arg, "-from=")) {
			gchar *value;
			value = arg + 6;
			if ((*err = read_cut_pair(value, &cut_args->cut_start))) {
				siril_log_color_message(_("Error: Could not parse -from values.\n"), "red");
				break;
			}
		}
		else if (g_str_has_prefix(arg, "-to=")) {
			gchar *value;
			value = arg + 4;
			if ((*err = read_cut_pair(value, &cut_args->cut_end))) {
				siril_log_color_message(_("Error: Could not parse -to values.\n"), "red");
				break;
			}
		}
/*		else if (g_str_has_prefix(arg, "-wn1at=")) {
			gchar *value;
			value = arg + 7;
			if ((*err = read_cut_pair(value, &cut_args->cut_wn1))) {
				siril_log_color_message(_("Error: Could not parse -wn1at values.\n"), "red");
				break;
			}
		}
		else if (g_str_has_prefix(arg, "-wn2at=")) {
			gchar *value;
			value = arg + 7;
			if ((*err = read_cut_pair(value, &cut_args->cut_wn2))) {
				siril_log_color_message(_("Error: Could not parse -wn2at values.\n"), "red");
				break;
			}
		}
*/		else if (g_str_has_prefix(arg, "-filename=")) {
			if (seq) {
				siril_log_color_message(_("Error: this option cannot be used for sequences.\n"), "red");
				*err = CMD_ARG_ERROR;
				break;
			}
			arg += 10;
			if (cut_args->filename)
				g_free(cut_args->filename);
			cut_args->filename = g_strdup(arg);
			cut_args->save_dat = TRUE;
		}
		else if (g_str_has_prefix(arg, "-title=")) {
			arg += 7;
			if (arg && arg[0] != '\0') {
				if (cut_args->user_title)
					g_free(cut_args->user_title);
				cut_args->user_title = g_strdup(arg);
				siril_debug_print("title: %s\n", arg);
			}
		}
	}
	if (cut_args->vport == -1) {
		if (nb_layers == 1) {
			cut_args->vport = 0;
			cut_args->mode = CUT_MONO;
		} else {
			cut_args->vport = 3;
			cut_args->mode = CUT_COLOR;
		}
	}
	if (seq && seq->is_variable) {
		siril_log_message(_("Error: sequence has variable sized images.\n"));
		*err = CMD_GENERIC_ERROR;
	} else if (!cut_struct_is_valid(cut_args)) {
		*err = CMD_ARG_ERROR;
	}
	if (*err) {
		free_cut_args(cut_args);
		cut_args = NULL;
	}
	return cut_args;
}

int process_profile(int nb) {
	cmd_errors err;
	cut_struct *cut_args = parse_cut_args(nb, NULL, &err);
	if (err)
		return err;

	cut_args->display_graph = (!com.script); // we can display the plot if not in a script
	cut_args->save_png_too = TRUE;

	if (cut_args->cfa) {
		if (!start_in_new_thread(cfa_cut, cut_args)) {
			free_cut_args(cut_args);
			return CMD_ARG_ERROR;
		}
	} else if (cut_args->tri) {
		if (!start_in_new_thread(tri_cut, cut_args)) {
			free_cut_args(cut_args);
			return CMD_ARG_ERROR;
		}
	} else {
		if (!start_in_new_thread(cut_profile, cut_args)) {
			free_cut_args(cut_args);
			return CMD_ARG_ERROR;
		}
	}
	return CMD_OK;
}

int process_seq_profile(int nb) {
	cmd_errors err;
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	cut_struct *cut_args = parse_cut_args(nb, seq, &err);
	if (err)
		return err;

	cut_args->display_graph = FALSE;
	cut_args->save_png_too = FALSE;

	apply_cut_to_sequence(cut_args);

	return CMD_OK;
}

int process_icc_assign(int nb) {
	if(!com.headless) on_clear_roi();
	char *arg = word[1];
	cmsHPROFILE profile = NULL;
	if (!g_ascii_strncasecmp(arg, "srgblinear", 10)) {
		profile = gfit->naxes[2] == 1 ? gray_linear() : srgb_linear();
	} else if (!g_ascii_strncasecmp(arg, "srgb", 4)) {
		profile = gfit->naxes[2] == 1 ? gray_srgbtrc() : srgb_trc();
	} else if (!g_ascii_strncasecmp(arg, "rec2020linear", 13)) {
		profile = gfit->naxes[2] == 1 ? gray_linear() : rec2020_linear();
	} else if (!g_ascii_strncasecmp(arg, "rec2020", 7)) {
		profile = gfit->naxes[2] == 1 ? gray_rec709trc() : rec2020_trc();
	} else if (!g_ascii_strncasecmp(arg, "linear", 6)) {
		profile = gfit->naxes[2] == 1 ? gray_linear() : siril_color_profile_linear_from_color_profile (com.icc.working_standard);
	} else if (!g_ascii_strncasecmp(arg, "working", 7)) {
		profile = copyICCProfile(gfit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard);
	} else if (g_file_test(arg, G_FILE_TEST_EXISTS) && g_file_test(arg, G_FILE_TEST_IS_REGULAR)) {
		profile = cmsOpenProfileFromFile(arg, "r");
	}
	if (profile) {
		if (gfit->icc_profile)
			cmsCloseProfile(gfit->icc_profile);
		gfit->icc_profile = NULL;
		siril_colorspace_transform(gfit, profile);
		cmsCloseProfile(profile);
	} else {
		siril_log_color_message(_("Error opening target ICC profile.\n"), "red");
		return CMD_GENERIC_ERROR;
	}
	if (gfit->icc_profile) {
		siril_log_color_message(_("Color profile assignment complete.\n"), "green");
		color_manage(gfit, TRUE);
	} else {
		siril_log_color_message(_("Error opening ICC profile.\n"), "red");
		color_manage(gfit, FALSE);
		return CMD_GENERIC_ERROR;
	}
	refresh_icc_transforms();
	if (!com.headless)
		notify_gfit_modified();

	return CMD_OK;
}

int process_icc_convert_to(int nb) {
	if (!com.headless) on_clear_roi();
	char *arg = word[1];
	cmsUInt32Number temp_intent = com.pref.icc.processing_intent;
	com.pref.icc.processing_intent = com.pref.icc.export_intent;
	if (!gfit->icc_profile) {
		siril_log_color_message(_("Image has no color profile assigned to convert from. Assign a profile first.\n"), "red");
		com.pref.icc.processing_intent = temp_intent;
		return CMD_GENERIC_ERROR;
	}
	if (word[2]) {
		if (!g_ascii_strcasecmp(word[2], "perceptual"))
			com.pref.icc.processing_intent = INTENT_PERCEPTUAL;
		else if (!g_ascii_strcasecmp(word[2], "relative"))
			com.pref.icc.processing_intent = INTENT_RELATIVE_COLORIMETRIC;
		else if (!g_ascii_strcasecmp(word[2], "saturation"))
			com.pref.icc.processing_intent = INTENT_SATURATION;
		else if (!g_ascii_strcasecmp(word[2], "absolute"))
			com.pref.icc.processing_intent = INTENT_ABSOLUTE_COLORIMETRIC;
		else {
			siril_log_color_message(_("Specified intent not recognized. Specify one of \"perceptual\", \"relative\", \"saturation\" or \"absolute\"\n"), "red");
			com.pref.icc.processing_intent = temp_intent;
			return CMD_GENERIC_ERROR;
		}
	}
	cmsHPROFILE profile = NULL;
	if (!g_ascii_strncasecmp(arg, "srgblinear", 10)) {
		profile = srgb_linear();
	} else if (!g_ascii_strncasecmp(arg, "srgb", 4)) {
		profile = srgb_trc();
	} else if (!g_ascii_strncasecmp(arg, "rec2020linear", 13)) {
		profile = rec2020_linear();
	} else if (!g_ascii_strncasecmp(arg, "rec2020", 7)) {
		profile = rec2020_trc();
	} else if (!g_ascii_strncasecmp(arg, "graysrgb", 8)) {
		profile = gray_srgbtrc();
	} else if (!g_ascii_strncasecmp(arg, "grayrec2020", 11)) {
		profile = gray_rec709trc();
	} else if (!g_ascii_strncasecmp(arg, "graylinear", 10)) {
		profile = gray_linear();
	} else if (!g_ascii_strncasecmp(arg, "working", 7)) {
		profile = copyICCProfile(gfit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard);
	} else if (g_file_test(arg, G_FILE_TEST_EXISTS) && g_file_test(arg, G_FILE_TEST_IS_REGULAR)) {
		profile = cmsOpenProfileFromFile(arg, "r");
	}
	if (profile) {
		siril_colorspace_transform(gfit, profile);
		gfit->icc_profile = copyICCProfile(profile);
		color_manage(gfit, TRUE);
		com.pref.icc.processing_intent = temp_intent;
		cmsCloseProfile(profile);
		siril_log_color_message(_("Color space conversion complete.\n"), "green");
	} else {
		siril_log_color_message(_("Error opening ICC profile.\n"), "red");
		com.pref.icc.processing_intent = temp_intent;
		return CMD_GENERIC_ERROR;
		// Don't call color_manage(gfit, FALSE) here: no change is made to
		// the pre-existing state of gfit color management
	}
	refresh_icc_transforms();
	if (!com.headless) {
		gui_function(close_tab, NULL);
		gui_function(init_right_tab, NULL);
	}
	if (!com.headless)
		notify_gfit_modified();
	return CMD_OK;
}

int process_icc_remove(int nb) {
	if (!com.headless) on_clear_roi();
	siril_colorspace_transform(gfit, NULL);
	refresh_icc_transforms();
	if (!com.headless)
		notify_gfit_modified();

	return CMD_OK;
}

int process_spcc_list(int nb) {
	GList *list = NULL;
	gboolean is_osc_sensor = FALSE;
	const gchar *list_name = NULL;
	load_spcc_metadata_if_needed();
	if (g_str_has_prefix(word[1], "oscsensor")) {
		list = com.spcc_data.osc_sensors;
		is_osc_sensor = TRUE;
		list_name = _("OSC Sensors");
	} else if (g_str_has_prefix(word[1], "monosensor")) {
		list = com.spcc_data.mono_sensors;
		list_name = _("Mono Sensors");
	} else if (g_str_has_prefix(word[1], "redfilter")) {
		list = com.spcc_data.mono_filters[RLAYER];
		list_name = _("Red Filters");
	} else if (g_str_has_prefix(word[1], "greenfilter")) {
		list = com.spcc_data.mono_filters[GLAYER];
		list_name = _("Green Filters");
	} else if (g_str_has_prefix(word[1], "bluefilter")) {
		list = com.spcc_data.mono_filters[BLAYER];
		list_name = _("Blue Filters");
	} else if (g_str_has_prefix(word[1], "oscfilter")) {
		list = com.spcc_data.osc_filters;
		list_name = _("OSC Filters");
	} else if (g_str_has_prefix(word[1], "osclpf")) {
		list = com.spcc_data.osc_lpf;
		list_name = _("OSC Low-Pass Filters");
	} else if (g_str_has_prefix(word[1], "whiteref")) {
		list = com.spcc_data.wb_ref;
		list_name = _("White References");
	} else {
		siril_log_message(_("Unknown SPCC list\n"));
		return CMD_ARG_ERROR;
	}
	siril_log_color_message("%s\n", "green", list_name);
	while (list) {
		const spcc_object *object = (const spcc_object*) list->data;
		siril_log_message("%s\n", is_osc_sensor ? object->model : object->name);
		list = list->next;
	}
	return CMD_OK;
}

int process_disto(int nb) {
	if (!has_wcs(gfit) || !gfit->keywords.wcslib->lin.dispre) {
		siril_log_color_message(_("This command only works on plate solved images with distortions included\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}
	if (nb > 2)
		return CMD_WRONG_N_ARG;
	if (nb == 1) {
		gui.show_wcs_disto	= TRUE;
		queue_redraw(REDRAW_OVERLAY);
		return CMD_OK;
	}
	if (!strcmp(word[1], "clear")) {
		gui.show_wcs_disto	= FALSE;
		queue_redraw(REDRAW_OVERLAY);
		return CMD_OK;
	} else {
		siril_log_message(_("Unknown parameter %s, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	return CMD_OK;
}

int process_trixel(int nb) {
	if (nb > 3)
		return CMD_WRONG_N_ARG;
	if (nb == 1) {
		if (!single_image_is_loaded())
			return CMD_LOAD_IMAGE_FIRST;
		if (!has_wcs(gfit))
			return CMD_FOR_PLATE_SOLVED;
		start_in_new_thread(list_trixels, NULL);
	} else if (!strcmp(word[2], "-p"))
		start_in_new_thread(write_trixels, NULL);
	else {
		siril_log_message(_("Unknown parameter %s, aborting.\n"), word[2]);
		return CMD_ARG_ERROR;
	}
	return CMD_OK;
}

int process_limit(int nb) {
	if (nb != 2)
		return CMD_WRONG_N_ARG;
	if (gfit->type == DATA_USHORT) {
		siril_log_message(_("16-bit images cannot have out-of-range pixels: nothing to do.\n"));
		return CMD_OK;
	}

	double maxval, minval;
	int retval = quick_minmax(gfit, &minval, &maxval);
	if (retval)
		return CMD_GENERIC_ERROR;

	if (maxval <= 1.0 && minval >= 0.0) {
		siril_log_message(_("No pixels require clipping. Nothing to do...\n"));
		return CMD_OK;
	}

	OverrangeResponse method;
	if (!g_ascii_strncasecmp(word[1], "-clip", 5)) {
		method = RESPONSE_CLIP;
	} else if (!g_ascii_strncasecmp(word[1], "-posrescale", 11)) {
		method = RESPONSE_RESCALE_CLIPNEG;
	} else if (!g_ascii_strncasecmp(word[1], "-rescale", 8)) {
		method = RESPONSE_RESCALE_ALL;
	} else {
		siril_log_color_message(_("Error: unknown argument!\n"), "red");
		return CMD_ARG_ERROR;
	}

	apply_limits(gfit, minval, maxval, method);
	siril_log_message(_("Pixel limits applied successfully.\n"));
	return CMD_OK;
}

int process_online(int nb) {
	set_online_status(TRUE);
	return CMD_OK;
}

int process_offline(int nb) {
	set_online_status(FALSE);
	return CMD_OK;
}

int process_pwd(int nb) {
	siril_log_message(_("Current working directory: '%s'\n"), com.wd);
	return CMD_OK;
}

typedef struct _pyscript_data {
	gchar *script_name;
	gchar **argv_script;
	gboolean from_cli;
} pyscript_data;

gpointer execute_python_script_wrapper(gpointer user_data) {
	pyscript_data *data = (pyscript_data*) user_data;
	execute_python_script(data->script_name, TRUE, TRUE, data->argv_script, FALSE, data->from_cli, FALSE);
	// execute_python_script() frees data->script_name
	g_strfreev(data->argv_script);
	free(data);
	return GINT_TO_POINTER(0);
}

int process_pyscript(int nb) {
	gchar *script_name = NULL;
	GStatBuf statbuf;
	if (g_stat(word[1], &statbuf) == 0) {
		// full path provided (or at least we can find it directly)
		script_name = g_strdup(word[1]);
	} else {
		// Search for the file in the user's set script directories and the scripts repository
		// We search the user's path first so any local modifications are used in preference to
		// the repository script with the same name.
		GSList *path = com.pref.gui.script_path;
		while (path) {
			siril_debug_print("Searching script path: %s\n", (gchar*) path->data);
			script_name = find_file_in_directory(word[1], (gchar*) path->data);
			if (script_name) // found it!
				break;
			path = path->next;
		}
		if (!script_name) {
			// If we still haven't found it, iterate over the siril-scripts local repository
			script_name = find_file_recursively(word[1], siril_get_scripts_repo_path());
		}
	}

	if (script_name) {
		gchar** argv_script = NULL;
		if (nb > 1) {
			// Treat additional arguments as arguments to be passed to the script
			argv_script = calloc(nb, sizeof(gchar*));
			for (int i = 2 ; i <= nb ; i++) {
				argv_script[i-2] = g_strdup(word[i]);
			}
		}
		pyscript_data *data = calloc(1, sizeof(pyscript_data));
		data->script_name = script_name;
		data->argv_script = argv_script;
		data->from_cli = TRUE;
		// Cannot use start_in_new_thread here because of the possibility of the python script
		// calling siril.cmd() and running commands that themselves require the processing thread
		// so we use a generic GThread
		gboolean already_in_a_python_script = com.python_script;
		GThread *thread = g_thread_new("pyscript_thread", execute_python_script_wrapper, data);
		// data->script_name freed by execute_python_script_wrapper
		if (!thread) {
			free(data);
			g_free(script_name);
			g_strfreev(argv_script);
			return CMD_GENERIC_ERROR;
		} else {
			if (com.script || already_in_a_python_script) {
				g_thread_join(thread);
			} else {
				g_thread_unref(thread);
			}
		}
		return CMD_OK;
	} else {
		return CMD_FILE_NOT_FOUND;
	}
}

// Process functions refactored
int process_mask_from_stars(int nb) {
	int argidx = 1;
	float r = 0.f, feather = 0.f;
	gboolean invert = FALSE;
	int bitdepth = com.pref.default_mask_bitpix;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-invert")) {
			invert = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-r=")) {
			char *arg = word[argidx] + 3;
			r = g_ascii_strtod(arg, &end);
			if (arg == end || r <= 0.0f || r > 100.f) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-feather=")) {
			char *arg = word[argidx] + 9;
			feather = g_ascii_strtod(arg, &end);
			if (arg == end || feather < 0.0f || feather > 2000.f) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-bitdepth=")) {
			char *arg = word[argidx] + 10;
			bitdepth = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || (bitdepth != 8 && bitdepth != 16 && bitdepth != 32)) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	mask_from_stars_data *data = calloc(1, sizeof(mask_from_stars_data));

	data->r = r;
	data->feather = feather;
	data->invert = invert;
	data->bitdepth = bitdepth;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_from_stars_hook;
	args->log_hook = mask_from_stars_log;
	args->idle_function = NULL;
	args->description = _("Mask from stars");
	args->verbose = TRUE;
	args->command = TRUE;
	args->mask_creation = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_mask_from_channel(int nb) {
	int argidx = 1;
	int channel = -1;
	gboolean autostretch = FALSE;
	gboolean invert = FALSE;
	int bitdepth = com.pref.default_mask_bitpix;
	gchar *filename = NULL;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-autostretch")) {
			autostretch = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-invert")) {
			invert = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-channel=")) {
			char *arg = word[argidx] + 9;
			channel = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || channel < 0 || channel > 2) {
				siril_log_message(_("Invalid argument %s, channel must be 0-2, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-bitdepth=")) {
			char *arg = word[argidx] + 10;
			bitdepth = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || (bitdepth != 8 && bitdepth != 16 && bitdepth != 32)) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-filename=")) {
			filename = word[argidx] + 10;
		}
		argidx++;
	}

	if (channel == -1) {
		siril_log_message(_("Channel parameter (-channel=) is required, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	mask_from_channel_data *data = calloc(1, sizeof(mask_from_channel_data));

	data->channel = channel;
	data->autostretch = autostretch;
	data->invert = invert;
	data->bitpix = (uint8_t)bitdepth;
	data->filename = filename;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_from_channel_hook;
	args->log_hook = mask_from_channel_log;
	args->idle_function = NULL;
	args->description = _("Mask from channel");
	args->verbose = TRUE;
	args->command = TRUE;
	args->mask_creation = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_mask_from_lum(int nb) {
	int argidx = 1;
	float rw = -1.f, gw = -1.f, bw = -1.f;
	gboolean autostretch = FALSE;
	gboolean invert = FALSE;
	gboolean use_human = FALSE;
	gboolean use_even = FALSE;
	int bitdepth = com.pref.default_mask_bitpix;
	gchar *filename = NULL;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-autostretch")) {
			autostretch = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-invert")) {
			invert = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-human")) {
			use_human = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-even")) {
			use_even = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-rw=")) {
			char *arg = word[argidx] + 4;
			rw = g_ascii_strtod(arg, &end);
			if (arg == end || rw < 0.0f || rw > 1.0f) {
				siril_log_message(_("Invalid argument %s, red weight must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-gw=")) {
			char *arg = word[argidx] + 4;
			gw = g_ascii_strtod(arg, &end);
			if (arg == end || gw < 0.0f || gw > 1.0f) {
				siril_log_message(_("Invalid argument %s, green weight must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-bw=")) {
			char *arg = word[argidx] + 4;
			bw = g_ascii_strtod(arg, &end);
			if (arg == end || bw < 0.0f || bw > 1.0f) {
				siril_log_message(_("Invalid argument %s, blue weight must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-bitdepth=")) {
			char *arg = word[argidx] + 10;
			bitdepth = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || (bitdepth != 8 && bitdepth != 16 && bitdepth != 32)) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-filename=")) {
			filename = word[argidx] + 10;
		}
		argidx++;
	}

	if (use_human && use_even) {
		siril_log_message(_("Cannot specify both -human and -even options, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if ((use_human || use_even) && (rw >= 0.f || gw >= 0.f || bw >= 0.f)) {
		siril_log_message(_("Cannot specify custom weights with -human or -even options, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	gboolean has_custom_weights = (rw >= 0.f || gw >= 0.f || bw >= 0.f);
	if (has_custom_weights && !(rw >= 0.f && gw >= 0.f && bw >= 0.f)) {
		siril_log_message(_("All three weights (-rw, -gw, -bw) must be specified, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (use_human) {
		rw = 0.2126f;
		gw = 0.7152f;
		bw = 0.0722f;
	}
	else if (use_even) {
		rw = 0.3333f;
		gw = 0.3333f;
		bw = 0.3334f;
	}
	else if (!has_custom_weights) {
		rw = 0.2126f;
		gw = 0.7152f;
		bw = 0.0722f;
	}

	mask_from_lum_data *data = calloc(1, sizeof(mask_from_lum_data));

	data->rw = rw;
	data->gw = gw;
	data->bw = bw;
	data->autostretch = autostretch;
	data->invert = invert;
	data->use_human = use_human;
	data->use_even = use_even;
	data->bitpix = (uint8_t)bitdepth;
	data->filename = filename;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_from_lum_hook;
	args->log_hook = mask_from_lum_log;
	args->idle_function = NULL;
	args->description = _("Mask from luminance");
	args->verbose = TRUE;
	args->command = TRUE;
	args->mask_creation = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_clear_mask(int nb) {
	if (gfit && gfit->mask) {
		set_mask_active(gfit, FALSE);
		free_mask(gfit->mask);
		gfit->mask = NULL;
		show_or_hide_mask_tab();
		siril_log_message(_("Mask cleared\n"));
	}
	if (!com.script) {
		execute_idle_and_wait_for_it(redraw_mask_idle, NULL);
	}
	return CMD_OK;
}

int process_binarize_mask(int nb) {
	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	int argidx = 1;
	float lo = 0.f, hi;
	switch (gfit->mask->bitpix) {
		case 8:
			hi = 255.f;
			break;
		case 16:
			hi = 65535.f;
			break;
		case 32:
			hi = 1.f;
			break;
		default:
			siril_log_message(_("Unknown mask bit depth, aborting.\n"));
			return CMD_GENERIC_ERROR;
	}
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-lo=")) {
			char *arg = word[argidx] + 4;
			lo = g_ascii_strtod(arg, &end);
			if (arg == end) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-hi=")) {
			char *arg = word[argidx] + 4;
			hi = g_ascii_strtod(arg, &end);
			if (arg == end) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	if (lo >= hi) {
		siril_log_message(_("Low value must be less than high value, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	mask_binarize_data *data = calloc(1, sizeof(mask_binarize_data));

	data->lo = lo;
	data->hi = hi;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_binarize_hook;
	args->log_hook = mask_binarize_log;
	args->idle_function = NULL;
	args->description = _("Binarize mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_blur_mask(int nb) {
	int argidx = 1;
	float radius = -1.f;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-r=")) {
			char *arg = word[argidx] + 3;
			radius = g_ascii_strtod(arg, &end);
			if (arg == end || radius <= 0.0f || radius > 1000.f) {
				siril_log_message(_("Invalid argument %s, radius must be > 0 and <= 1000, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	if (radius < 0.f) {
		siril_log_message(_("Radius parameter (-radius=) is required, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	mask_blur_data *data = calloc(1, sizeof(mask_blur_data));

	data->radius = radius;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 2.0f;
	args->mask_hook = mask_blur_hook;
	args->log_hook = mask_blur_log;
	args->idle_function = NULL;
	args->description = _("Blur mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_feather_mask(int nb) {
	int argidx = 1;
	float distance = -1.f;
	feather_mode mode = FEATHER_OUTER;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-dist=")) {
			char *arg = word[argidx] + 6;
			distance = g_ascii_strtod(arg, &end);
			if (arg == end || distance <= 0.0f || distance > 2000.f) {
				siril_log_message(_("Invalid argument %s, distance must be > 0 and <= 2000, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-mode=")) {
			char *arg = word[argidx] + 6;
			if (g_ascii_strcasecmp(arg, "inner") == 0) {
				mode = FEATHER_INNER;
			}
			else if (g_ascii_strcasecmp(arg, "outer") == 0) {
				mode = FEATHER_OUTER;
			}
			else if (g_ascii_strcasecmp(arg, "edge") == 0) {
				mode = FEATHER_EDGE;
			}
			else {
				siril_log_message(_("Invalid argument %s, mode must be 'inner', 'outer', or 'edge', aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	if (distance < 0.f) {
		siril_log_message(_("Distance parameter (-distance=) is required, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	mask_feather_data *data = calloc(1, sizeof(mask_feather_data));

	data->distance = distance;
	data->mode = mode;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 2.0f;
	args->mask_hook = mask_feather_hook;
	args->log_hook = mask_feather_log;
	args->idle_function = NULL;
	args->description = _("Feather mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_mask_fmul(int nb) {
	int argidx = 1;
	float factor = -1.f;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-factor=")) {
			char *arg = word[argidx] + 8;
			factor = g_ascii_strtod(arg, &end);
			if (arg == end || factor < 0.0f) {
				siril_log_message(_("Invalid argument %s, factor must be >= 0, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	if (factor < 0.f) {
		siril_log_message(_("Factor parameter (-factor=) is required, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	mask_fmul_data *data = calloc(1, sizeof(mask_fmul_data));

	data->factor = factor;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 0.0f;
	args->mask_hook = mask_fmul_hook;
	args->log_hook = mask_fmul_log;
	args->idle_function = NULL;
	args->description = _("Multiply mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_invert_mask(int nb) {
	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 0.0f;
	args->mask_hook = mask_invert_hook;
	args->log_hook = NULL;
	args->idle_function = NULL;
	args->description = _("Invert mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = NULL;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_autostretch_mask(int nb) {
	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 0.0f;
	args->mask_hook = mask_autostretch_hook;
	args->log_hook = NULL;
	args->idle_function = NULL;
	args->description = _("Autostretch mask");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = NULL;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_mask_bitpix(int nb) {
	uint8_t bitpix;
	char* end;
	char* arg = word[1];

	size_t bitdepth = g_ascii_strtoull(arg, &end, 10);
	if (arg == end || (bitdepth != 8 && bitdepth != 16 && bitdepth != 32)) {
		siril_log_message(_("Invalid argument %s, bitdepth must be 8, 16 or 32, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	bitpix = (uint8_t) bitdepth;

	if (!gfit->mask || !gfit->mask->data) {
		siril_log_message(_("No mask present, aborting.\n"));
		return CMD_GENERIC_ERROR;
	}

	mask_bitpix_data *data = calloc(1, sizeof(mask_bitpix_data));

	data->bitpix = bitpix;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_bitpix_hook;
	args->log_hook = mask_bitpix_log;
	args->idle_function = NULL;
	args->description = _("Change mask bitdepth");
	args->verbose = TRUE;
	args->command = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}

int process_mask_from_color(int nb) {
	int argidx = 1;
	float chrom_center_r = -1.f, chrom_center_g = -1.f, chrom_center_b = -1.f;
	float chrom_tolerance = -1.f;
	float lum_min = 0.f, lum_max = 1.f;
	int feather_radius = 0;
	gboolean invert = FALSE;
	int bitdepth = com.pref.default_mask_bitpix;
	char *end;

	while (argidx < nb) {
		if (g_str_has_prefix(word[argidx], "-invert")) {
			invert = TRUE;
		}
		else if (g_str_has_prefix(word[argidx], "-cr=")) {
			char *arg = word[argidx] + 4;
			chrom_center_r = g_ascii_strtod(arg, &end);
			if (arg == end || chrom_center_r < 0.0f || chrom_center_r > 1.0f) {
				siril_log_message(_("Invalid argument %s, chroma center red must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-cg=")) {
			char *arg = word[argidx] + 4;
			chrom_center_g = g_ascii_strtod(arg, &end);
			if (arg == end || chrom_center_g < 0.0f || chrom_center_g > 1.0f) {
				siril_log_message(_("Invalid argument %s, chroma center green must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-cb=")) {
			char *arg = word[argidx] + 4;
			chrom_center_b = g_ascii_strtod(arg, &end);
			if (arg == end || chrom_center_b < 0.0f || chrom_center_b > 1.0f) {
				siril_log_message(_("Invalid argument %s, chroma center blue must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-tol=")) {
			char *arg = word[argidx] + 5;
			chrom_tolerance = g_ascii_strtod(arg, &end);
			if (arg == end || chrom_tolerance < 0.0f || chrom_tolerance > 1.0f) {
				siril_log_message(_("Invalid argument %s, chroma tolerance must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-lum_min=")) {
			char *arg = word[argidx] + 9;
			lum_min = g_ascii_strtod(arg, &end);
			if (arg == end || lum_min < 0.0f || lum_min > 1.0f) {
				siril_log_message(_("Invalid argument %s, luminance minimum must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-lum_max=")) {
			char *arg = word[argidx] + 9;
			lum_max = g_ascii_strtod(arg, &end);
			if (arg == end || lum_max < 0.0f || lum_max > 1.0f) {
				siril_log_message(_("Invalid argument %s, luminance maximum must be 0-1, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-fr=")) {
			char *arg = word[argidx] + 4;
			feather_radius = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || feather_radius < 0) {
				siril_log_message(_("Invalid argument %s, feather radius must be >= 0, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[argidx], "-bitdepth=")) {
			char *arg = word[argidx] + 10;
			bitdepth = (int) g_ascii_strtoull(arg, &end, 10);
			if (arg == end || (bitdepth != 8 && bitdepth != 16 && bitdepth != 32)) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[argidx]);
				return CMD_ARG_ERROR;
			}
		}
		argidx++;
	}

	if (chrom_center_r < 0.f || chrom_center_g < 0.f || chrom_center_b < 0.f) {
		siril_log_message(_("All three chroma center values (-cr, -cg, -cb) must be specified, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (chrom_tolerance < 0.f) {
		siril_log_message(_("Chroma tolerance (-tolerance) must be specified, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	if (lum_min > lum_max) {
		siril_log_message(_("Luminance minimum cannot be greater than luminance maximum, aborting.\n"));
		return CMD_ARG_ERROR;
	}

	mask_from_color_data *data = calloc(1, sizeof(mask_from_color_data));

	data->chrom_center_r = chrom_center_r;
	data->chrom_center_g = chrom_center_g;
	data->chrom_center_b = chrom_center_b;
	data->chrom_tolerance = chrom_tolerance;
	data->lum_min = lum_min;
	data->lum_max = lum_max;
	data->feather_radius = feather_radius;
	data->invert = invert;
	data->bitpix = (uint8_t)bitdepth;

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit =  gfit;
	args->mem_ratio = 1.5f;
	args->mask_hook = mask_from_color_hook;
	args->log_hook = mask_from_color_log;
	args->idle_function = NULL;
	args->description = _("Mask from color");
	args->verbose = TRUE;
	args->command = TRUE;
	args->mask_creation = TRUE;
	args->user = data;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
	return CMD_OK;
}
