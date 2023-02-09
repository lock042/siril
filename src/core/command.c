/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_histogram.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <opencv2/core/version.hpp>
#include <glib.h>
#include <libgen.h>
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#endif
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_JSON_GLIB
#include <json-glib/json-glib.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/initfile.h"
#include "core/preprocess.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "io/Astro-TIFF.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/catalogues.h"
#include "io/FITS_symlink.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/linear_match.h"
#include "gui/newdeconv.h"
#include "gui/sequence_list.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"
#include "gui/script_menu.h"
#include "gui/preferences.h"
#include "filters/asinh.h"
#include "filters/banding.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "filters/clahe.h"
#include "filters/cosmetic_correction.h"
#include "filters/deconvolution/deconvolution.h"
#include "filters/median.h"
#include "filters/mtf.h"
#include "filters/fft.h"
#include "filters/rgradient.h"
#include "filters/saturation.h"
#include "filters/scnr.h"
#include "filters/starnet.h"
#include "filters/synthstar.h"
#include "filters/wavelets.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/search_objects.h"
#include "algos/star_finder.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/background_extraction.h"
#include "algos/ccd-inspector.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"
#include "algos/colors.h"
#include "algos/quality.h"
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "algos/geometry.h"
#include "algos/photometric_cc.h"
#include "opencv/opencv.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "algos/fix_xtrans_af.h"
#include "algos/annotate.h"
#include "livestacking/livestacking.h"
#include "pixelMath/pixel_math_runner.h"
#include "git-version.h"

#include "command.h"
#include "command_def.h"
#include "command_list.h"
#include "command_line_processor.h"

/* process_command functions take the number of arguments (like argc) as
 * parameter and will be able to access the equivalent of argv with `word'
 * they return CMD_OK on success
 */

char *word[MAX_COMMAND_WORDS];	// NULL terminated

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
	launch_clipboard_survey();
	return retval;
}

int process_dumpheader(int nb) {
	if (!gfit.header || gfit.header[0] == '\0') {
		siril_log_message(_("Currently loaded image has no FITS header\n"));
	} else {
		siril_log_message(_("=FITS header for currently loaded image=\n"));
		char *header = strdup(gfit.header);
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
		free_sequence(seq, FALSE);
	}
	return CMD_OK;
}

int process_satu(int nb){
	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	args->background_factor = 1.0;
	gchar *end;
	args->coeff = g_ascii_strtod(word[1], &end);
	if (end == word[1]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
		free(args);
		return CMD_ARG_ERROR;
	}
	if (nb > 2) {
		args->background_factor = g_ascii_strtod(word[2], &end);
		if (end == word[2] || args->background_factor < 0.0) {
			siril_log_message(_("Background factor must be positive\n"));
			free(args);
			return CMD_ARG_ERROR;
		}
	}
	int satu_hue_type = 6;
	if (nb > 3) {
		satu_hue_type = g_ascii_strtoull(word[3], &end, 10);
		if (end == word[3] || satu_hue_type > 6) {
			siril_log_message(_("Hue range must be [0, 6]\n"));
			free(args);
			return CMD_ARG_ERROR;
		}
	}

	satu_set_hues_from_types(args, satu_hue_type);
	args->input = &gfit;
	args->output = &gfit;
	args->for_preview = FALSE;
	siril_log_message(_("Applying saturation enhancement of %d%%, from level %g times (median + sigma).\n"), round_to_int(args->coeff * 100.0), args->background_factor);

	set_cursor_waiting(TRUE);
	return GPOINTER_TO_INT(enhance_saturation(args));
}

int process_save(int nb){
	gchar *filename = g_strdup(word[1]);
	if (!com.script) {
		gfit.lo = gui.lo;
		gfit.hi = gui.hi;
	}
	int status, retval;
	gchar *savename = update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		retval = savefits(savename, &gfit);
		set_cursor_waiting(FALSE);
	}
	set_precision_switch();
	g_free(filename);
	g_free(savename);
	return retval;
}

int process_savebmp(int nb){
	gchar *filename = g_strdup_printf("%s.bmp", word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		retval = savebmp(savename, &gfit);
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}

int process_synthstar(int nb) {
	set_cursor_waiting(TRUE);
	start_in_new_thread(do_synthstar, NULL);
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

int process_unclip(int nb) {
	set_cursor_waiting(TRUE);
	start_in_new_thread(fix_saturated_stars, NULL);
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

static gboolean end_denoise(gpointer p) {
	struct denoise_args *args = (struct denoise_args *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	free(args);
	return FALSE;
}

gpointer run_nlbayes_on_fit(gpointer p) {
	denoise_args *args = (denoise_args *) p;
	struct timeval t_start, t_end;
	char *msg1 = NULL, *msg2 = NULL, *msg3 = NULL, *log_msg = NULL;
	int n = 0, m = 0, q = 0;
	n = snprintf(NULL, 0, _("NL-Bayes denoise (mod=%.3f"), args->modulation);
	msg1 = malloc(n + 1);
	n = snprintf(msg1, n + 1, _("NL-Bayes denoise (mod=%.3f"), args->modulation);
	if(args->da3d) {
		m = snprintf(NULL, 0, _(", DA3D enabled"));
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", DA3D enabled"));
	} else if (args->sos > 1) {
		m = snprintf(NULL, 0, _(", SOS enabled (iters=%d, rho=%.3f)"), args->sos, args->rho);
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", SOS enabled (iters=%d, rho=%.3f)"), args->sos, args->rho);
	} else if (args->do_anscombe) {
		m = snprintf(NULL, 0, _(", VST enabled"));
		msg2 = malloc(m + 1);
		m = snprintf(msg2, m + 1, _(", VST enabled"));
	}
	if (args->do_cosme) {
		q = snprintf(NULL, 0, _(", CC enabled)"));
		msg3 = malloc(q + 1);
		q = snprintf(msg3, q + 1, _(", CC enabled)"));
	} else {
		q = 1;
		msg3 = malloc(q + 1);
		q = snprintf(msg3, q + 1, _(")"));
	}
	log_msg = malloc(n + m + q + 1);
	if (m == 0 && q == 0)
		snprintf(log_msg, n + 1, "%s", msg1);
	else if (m > 0 && q == 0)
		snprintf(log_msg, n + m + 1, "%s%s", msg1, msg2);
	else if(m == 0 && q > 0)
		snprintf(log_msg, n + q + 1, "%s%s", msg1, msg3);
	else if (m > 0 && q > 0)
		snprintf(log_msg, n + m + q + 1, "%s%s%s", msg1, msg2, msg3);
	else
		snprintf(log_msg, 26, _("Error, this can't happen!"));

	if (msg1) free(msg1);
	if (msg2) free(msg2);
	if (msg3) free(msg3);

	siril_log_message("%s\n", log_msg);
	if (args->suppress_artefacts)
		siril_log_message(_("Colour artefact suppression active.\n"));
	free(log_msg);
	gettimeofday(&t_start, NULL);
	set_progress_bar_data(_("Starting NL-Bayes denoising..."), 0.0);

	int retval = 0;

	// Carry out cosmetic correction at the start, if selected
	if (args->do_cosme)
		denoise_hook_cosmetic(args->fit);

	if (args->fit == &gfit && args->fit->naxes[2] == 3 && args->suppress_artefacts) {
		fits *loop = NULL;
		if (new_fit_image(&loop, args->fit->rx, args->fit->ry, 1, args->fit->type)) {
			retval = 1;
		}
		loop->naxis = 1;
		loop->naxes[2] = 1;
		size_t npixels = args->fit->naxes[0] * args->fit->naxes[1];
		if (retval == 0) {
			if (args->fit->type == DATA_FLOAT) {
				for (size_t i = 0; i < 3; i++) {
					float *loop_fdata = (float*) calloc(npixels, sizeof(float));
					loop->fdata = loop_fdata;
					for (size_t j = 0 ; j < npixels ; j++) {
						loop_fdata[j] = args->fit->fpdata[i][j];
					}
					retval = do_nlbayes(loop, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
					for (size_t j = 0 ;j < npixels ; j++) {
						args->fit->fpdata[i][j] = loop->fdata[j];
					}
				}
				free(loop->fdata);
				loop->fdata = NULL;
			} else {
				for (size_t i = 0; i < 3; i++) {
					WORD *loop_data = (WORD*) calloc(npixels, sizeof(WORD));
					loop->data = loop_data;
					for (size_t j = 0 ; j < npixels ; j++) {
						loop_data[j] = args->fit->pdata[i][j];
					}
					retval = do_nlbayes(loop, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
					for (size_t j = 0 ;j < npixels ; j++) {
						args->fit->pdata[i][j] = loop->data[j];
					}
				}
				free(loop->data);
				loop->fdata = NULL;
			}
		clearfits(loop);
		}
	} else {
		retval = do_nlbayes(args->fit, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe);
	}
	if (args->fit == &gfit)
		notify_gfit_modified();
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, _("NL-Bayes execution time"));
	set_progress_bar_data(_("Ready."), 0.0);
	siril_add_idle(end_denoise, args);
	return GINT_TO_POINTER(retval);
}

int process_denoise(int nb){
	gboolean error = FALSE;
	set_cursor_waiting(TRUE);
	denoise_args *args = calloc(1, sizeof(denoise_args));
	args->sos = 1;
	args->rho = 0.2f;
	args->modulation = 1.f;
	args->da3d = 0;
	args->do_anscombe = FALSE;
	args->do_cosme = TRUE;
	args->suppress_artefacts = FALSE;
	args->fit = &gfit;
	for (int i = 1; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(arg, "-vst")) {
			args->do_anscombe = TRUE;
		}
		else if (g_str_has_prefix(arg, "-da3d")) {
			args->da3d = 1;
		}
		else if (g_str_has_prefix(arg, "-indep")) {
			args->suppress_artefacts = TRUE;
		}
		else if (g_str_has_prefix(arg, "-nocosmetic")) {
			args->do_cosme = FALSE;
		}
		else if (g_str_has_prefix(arg, "-mod=")) {
			if (args->modulation == 0.f) {
				siril_log_message(_("Modulation is zero: doing nothing.\n"));
				g_free(args);
				return CMD_OK;
			}
		}
		else if (g_str_has_prefix(arg, "-rho=")) {
			arg += 5;
			float rho = g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((rho <= 0.f) || (rho >= 1.f)) {
				siril_log_message(_("Error in rho parameter: must be strictly > 0 and < 1, aborting.\n"));
				g_free(args);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				args->rho = rho;
			}
		}
		else if (g_str_has_prefix(arg, "-sos=")) {
			arg += 5;
			unsigned sos = (int) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if (sos < 1) {
				siril_log_message(_("Error: SOS iterations must be >= 1. Defaulting to 1.\n"));
				sos = 1;
			} else if (sos > 10)
				siril_log_message(_("Note: high number of SOS iterations. Processing time may be lengthy...\n"));
			if (!error) {
				args->sos = sos;
			}
		}
	}
	if (args->do_anscombe && (args->sos != 1 || args->da3d)) {
		siril_log_color_message(_("Error: will not carry out DA3D or SOS iterations with Anscombe transform VST selected. aborting.\n"), "red");
		g_free(args);
		return CMD_ARG_ERROR;
	}
	if (args->do_anscombe)
		siril_log_message(_("Will apply generalised Anscombe variance stabilising transform.\n"));
	if (args->da3d) {
		siril_log_message(_("Will carry out final stage DA3D denoising.\n"));
		if (args->sos != 1) {
			siril_log_message(_("Will not carry out both DA3D and SOS. SOS iterations set to 1.\n"));
			args->sos = 1;
		}
	}
	start_in_new_thread(run_nlbayes_on_fit, args);

	return CMD_OK;
}

int process_starnet(int nb){
#ifdef HAVE_LIBTIFF
	starnet_data *starnet_args = malloc(sizeof(starnet_data));
	memset(starnet_args->stride, 0, sizeof(starnet_args->stride));
	starnet_args->linear = FALSE;
	starnet_args->seq = NULL;
	starnet_args->customstride = FALSE;
	starnet_args->upscale = FALSE;
	starnet_args->starmask = TRUE;
	starnet_args->follow_on = FALSE;
	starnet_args->starnet_fit = &gfit;
	gboolean error = FALSE;
	for (int i = 1; i < nb; i++) {
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
			double stride = g_ascii_strtod(arg, &end);
			int intstride = stride;
			if (arg == end) error = TRUE;
			else if ((intstride < 2.0) || (intstride > 256) || (intstride % 2)) {
				siril_log_message(_("Error in stride parameter: must be a positive even integer, max 256, aborting.\n"));
				g_free(starnet_args);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				sprintf(starnet_args->stride, "%d", intstride);
				starnet_args->customstride = TRUE;
			}
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			g_free(starnet_args);
			return CMD_ARG_ERROR;
		}
		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
			g_free(starnet_args);
			return CMD_ARG_ERROR;
		}
	}

	set_cursor_waiting(TRUE);
	start_in_new_thread(do_starnet, starnet_args);

#else
	siril_log_message(_("starnet command unavailable as Siril has not been compiled with libtiff.\n"));
#endif
	return CMD_OK;
}

int process_seq_starnet(int nb){
#ifdef HAVE_LIBTIFF
	starnet_data *starnet_args = calloc(1, sizeof(starnet_data));
	if (!starnet_args)
		return CMD_ALLOC_ERROR;
	memset(starnet_args->stride, 0, sizeof(starnet_args->stride));
	starnet_args->linear = FALSE;
	starnet_args->customstride = FALSE;
	starnet_args->upscale = FALSE;
	starnet_args->starmask = TRUE;
	starnet_args->follow_on = FALSE;
	starnet_args->starnet_fit = &gfit;
	gboolean error = FALSE;
	starnet_args->seq = load_sequence(word[1], NULL);
	if (!starnet_args->seq) {
		siril_log_message(_("Error: cannot open sequence\n"));
		g_free(starnet_args);
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
			double stride = g_ascii_strtod(arg, &end);
			int intstride = stride;
			if (arg == end) error = TRUE;
			else if ((intstride < 2.0) || (intstride > 256) || (intstride % 2)) {
				siril_log_message(_("Error in stride parameter: must be a positive even integer, max 256, aborting.\n"));
				if (!check_seq_is_comseq(starnet_args->seq))
					free_sequence(starnet_args->seq, TRUE);
				g_free(starnet_args);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				sprintf(starnet_args->stride, "%d", intstride);
				starnet_args->customstride = TRUE;
			}
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				if (!check_seq_is_comseq(starnet_args->seq))
					free_sequence(starnet_args->seq, TRUE);
				g_free(starnet_args);
			return CMD_ARG_ERROR;
		}
		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
			if (!check_seq_is_comseq(starnet_args->seq))
				free_sequence(starnet_args->seq, TRUE);
			g_free(starnet_args);
			return CMD_ARG_ERROR;
		}
	}
	set_cursor_waiting(TRUE);
	apply_starnet_to_sequence(starnet_args);

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
	gchar *savename = update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		retval = savejpg(savename, &gfit, quality);
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
	gchar *savename =update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		uint32_t bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
		retval = savepng(savename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
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
	gchar *astro_tiff = NULL;

	if (strcasecmp(word[0], "savetif8") == 0)
		bitspersample = 8;
	else if (strcasecmp(word[0], "savetif32") == 0)
		bitspersample = 32;
	if (word[2] && !g_strcmp0(word[2], "-astro")) {
		astro_tiff = AstroTiff_build_header(&gfit);
	}

	gchar *filename = g_strdup_printf("%s.tif", word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		retval = savetif(filename, &gfit, bitspersample, astro_tiff, com.pref.copyright, TRUE);
		set_cursor_waiting(FALSE);
	}
	g_free(astro_tiff);
	g_free(filename);
	g_free(savename);
	return retval;
}
#endif

int process_savepnm(int nb){
	gchar *filename = g_strdup(word[1]);
	int status, retval;
	gchar *savename = update_header_and_parse(&gfit, filename, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	if (status > 0) {
		retval = 1;
	} else {
		set_cursor_waiting(TRUE);
		retval = saveNetPBM(savename, &gfit);
		set_cursor_waiting(FALSE);
	}
	g_free(filename);
	g_free(savename);
	return retval;
}

static gboolean merge_cfa_idle(gpointer arg) {
	initialize_display_mode();
	update_zoom_label();
	display_filename();
	set_precision_switch();
	sliders_mode_set_state(gui.sliders);
	init_layers_hi_and_lo_values(MIPSLOHI);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();
	set_display_mode();
	redraw(REMAP_ALL);
	sequence_list_change_current();
	return FALSE;
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

	set_cursor_waiting(TRUE);

	strncpy(filename, word[1], 250);
	filename[250] = '\0';
	expand_home_in_filename(filename, maxpath);
	int retval = readfits(filename, &cfa0, NULL, FALSE);
	strncpy(filename, word[2], 250);
	filename[250] = '\0';
	expand_home_in_filename(filename, maxpath);
	retval += readfits(filename, &cfa1, NULL, FALSE);
	strncpy(filename, word[3], 250);
	filename[250] = '\0';
	expand_home_in_filename(filename, maxpath);
	retval += readfits(filename, &cfa2, NULL, FALSE);
	strncpy(filename, word[4], 250);
	filename[250] = '\0';
	expand_home_in_filename(filename, maxpath);
	retval += readfits(filename, &cfa3, NULL, FALSE);
	if (retval) {
		siril_log_color_message(_("Error loading files!\n"), "red");
		return CMD_FILE_NOT_FOUND;
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

	out = merge_cfa(&cfa0, &cfa1, &cfa2, &cfa3, pattern);
	siril_log_message("Bayer pattern produced: 1 layer, %dx%d pixels\n", out->rx, out->ry);
	close_single_image();
	copyfits(out, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	clear_stars_list(TRUE);
	com.seq.current = UNRELATED_IMAGE;
	if (!create_uniq_from_gfit(strdup(_("Unsaved Bayer pattern merge")), FALSE))
		com.uniq->comment = strdup(_("Bayer pattern merge"));
	open_single_image_from_gfit();
	set_cursor_waiting(FALSE);
	siril_add_idle(merge_cfa_idle, NULL);
	return CMD_OK;
}

int process_imoper(int nb){
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

	fits fit = { 0 };
	if (readfits(word[1], &fit, NULL, !com.pref.force_16bit)) return CMD_INVALID_IMAGE;

	int retval = imoper(&gfit, &fit, oper, !com.pref.force_16bit);

	clearfits(&fit);
	notify_gfit_modified();
	return retval;
}

int process_addmax(int nb){
	fits fit = { 0 };

	if (readfits(word[1], &fit, NULL, gfit.type == DATA_FLOAT))
		return CMD_INVALID_IMAGE;
	if (addmax(&gfit, &fit) == 0) {
		adjust_cutoff_from_updated_gfit();
		redraw(REMAP_ALL);
		redraw_previews();
	}
	clearfits(&fit);
	return CMD_OK;
}

int process_fdiv(int nb){
	// combines an image division and a scalar multiplication.
	gchar *end;
	float norm = g_ascii_strtod(word[2], &end);
	if (end == word[2]) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[2]);
		return CMD_ARG_ERROR;
	}

	fits fit = { 0 };
	if (readfits(word[1], &fit, NULL, !com.pref.force_16bit)) return CMD_INVALID_IMAGE;
	siril_fdiv(&gfit, &fit, norm, TRUE);

	clearfits(&fit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_fmul(int nb){
	gboolean from8b = (gfit.orig_bitpix == BYTE_IMG); // get orig bitdepth of 8b before it gets converted to 32b by soper
	gchar *end;
	float coeff = g_ascii_strtod(word[1], &end);
	if (end == word[1] || coeff <= 0.f) {
		siril_log_message(_("Multiplying by a coefficient less than or equal to 0 is not possible.\n"));
		return CMD_ARG_ERROR;
	}
	soper(&gfit, coeff, OPER_MUL, TRUE);
	if (from8b) { // image is now 32b, need to reset slider max and update hi/lo
		invalidate_stats_from_fit(&gfit);
		image_find_minmax(&gfit);
		gfit.hi = (WORD)(gfit.maxi * USHRT_MAX_SINGLE);
		gfit.lo = (WORD)(gfit.mini * USHRT_MAX_SINGLE);
		set_cutoff_sliders_max_values();
	}

	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_entropy(int nb){
	rectangle area;
	float e = 0.f;

	if (com.selection.w > 0 && com.selection.h > 0) {
		area = com.selection;
		for (int c = 0; c < gfit.naxes[2]; c++)
			e += entropy(&gfit, c, &area, NULL);
	}
	else {
		for (int c = 0; c < gfit.naxes[2]; c++)
			e += entropy(&gfit, c, NULL, NULL);
	}
	siril_log_message(_("Entropy: %.3f\n"), e);
	return CMD_OK;
}

int process_gauss(int nb){
	gchar *end;
	double sigma = g_ascii_strtod(word[1], &end);
	if (end == word[1] || sigma <= 0.0) {
		siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
		return CMD_ARG_ERROR;
	}
	unsharp(&gfit, sigma, 0.0, TRUE);
	//gaussian_blur_RT(&gfit, sigma, com.max_thread);

	char log[90];
	sprintf(log, "Gaussian filtering, sigma: %.2f", sigma);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
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
		fit_sequence_get_image_filename(seq, ref_image, filename, TRUE);
		siril_log_message(_("Image %d: '%s'\n"), ref_image, filename);
	}
	else siril_log_message(_("Image %d\n"), ref_image);

	if (!seq->imgparam[ref_image].incl)
		siril_log_message(_("Warning: this image is excluded from the sequence main processing list\n"));
	return CMD_OK;
}

int process_grey_flat(int nb) {
	if (isrgb(&gfit)) {
		return CMD_FOR_CFA_IMAGE;
	}

	compute_grey_flat(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();

	return CMD_OK;
}

int process_makepsf(int nb) {
	gboolean error = FALSE;
	estk_data* data = malloc(sizeof(estk_data));
	reset_conv_args(data);

	char *arg = word[1];
	if (!word[1]) {
		g_free(data);
		return CMD_WRONG_N_ARG;
	}
	if (g_str_has_prefix(arg, "clear")) {
		if (get_thread_run()) {
			siril_log_message(_("Error: will not clear the PSF while a sequence is running.\n"));
			g_free(data);
			return CMD_GENERIC_ERROR;
		}
		if (com.kernel) {
			free(com.kernel);
			com.kernel = NULL;
			com.kernelsize = 0;
		}
		siril_log_color_message(_("Deconvolution kernel cleared.\n"), "green");
		g_free(data);
		return CMD_OK;
	} else {
		if (g_str_has_prefix(arg, "save")) {
			siril_log_message(_("Save PSF to file:\n"));
			on_bdeconv_savekernel_clicked(NULL, NULL);
			g_free(data);
			return CMD_OK;
		}
		if (com.kernel) {
			free(com.kernel);
			com.kernel = NULL;
			com.kernelsize = 0;
		}
		if (g_str_has_prefix(arg, "blind")) {
			if (!(single_image_is_loaded() || sequence_is_loaded())) {
				siril_log_message(_("Error: image or sequence must be loaded to carry out blind PSF estimation. Aborting...\n"));
				g_free(data);
				return CMD_GENERIC_ERROR;
			}
			data->psftype = PSF_BLIND;
			siril_log_message(_("Blind kernel estimation:\n"));
			for (int i = 2; i < nb; i++) {
				char *arg = word[i], *end;
				if (!word[i])
					break;
				if (g_str_has_prefix(arg, "-l0")) {
					siril_log_message(_("ℓ0 descent prior method\n"));
					data->blindtype = 1;
				}
				else if (g_str_has_prefix(arg, "-si")) {
					siril_log_message(_("spectral irregularity method\n"));
					data->blindtype = 0;
				}
				else if (g_str_has_prefix(arg, "-multiscale")) {
					siril_log_message(_("multiscale estimation\n"));
					data->multiscale = TRUE;
				}
				else if (g_str_has_prefix(arg, "-lambda=")) {
					arg += 8;
					float lambda = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((lambda < 0.f) || (lambda > 100000.f)) {
						siril_log_message(_("Error in alpha parameter: must be between 0 and 1e5, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->lambda = lambda;
						data->finaldeconvolutionweight = lambda;
						data->intermediatedeconvolutionweight = lambda;
					}
				}
				else if (g_str_has_prefix(arg, "-comp=")) {
					arg += 6;
					float comp = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((comp < 1.f) || (comp > 100000.f)) {
						siril_log_message(_("Error in compensation factor parameter: must be between 1 and 1e5, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->compensationfactor = comp;
					}
				}
				else if (g_str_has_prefix(arg, "-ks=")) {
					arg += 4;
					int ks = (int) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((ks < 3) || !(ks %2) || (ks > min(gfit.rx, gfit.ry))) {
						siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width): aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->ks = ks;
					}
				}
			}
			if (error) {
				g_free(data);
				return CMD_ARG_ERROR;
			}
			start_in_new_thread(estimate_only, data);
			return CMD_OK;
		} else if (g_str_has_prefix(arg, "stars")) {
			if (!(single_image_is_loaded() || sequence_is_loaded())) {
				siril_log_message(_("Error: image or sequence must be loaded to carry out blind PSF estimation. Aborting...\n"));
				g_free(data);
				return CMD_GENERIC_ERROR;
			}
			if (!(com.stars && com.stars[0])) {
				siril_log_message(_("Error: requested to generate PSF from stars but no stars have been selected. Run findstar first.\n"));
				g_free(data);
				return CMD_ARG_ERROR;
			} else {
				data->psftype = PSF_STARS;
				for (int i = 2; i < nb; i++) {
					char *arg = word[i], *end;
					if (!word[i])
						break;
					if (g_str_has_prefix(arg, "-sym")) {
					siril_log_message(_("symmetric kernel\n"));
						data->symkern = TRUE;
					}
					else if (g_str_has_prefix(arg, "-ks=")) {
						arg += 4;
						int ks = (int) g_ascii_strtod(arg, &end);
						if (arg == end) error = TRUE;
						else if ((ks < 3) || !(ks %2) || (ks > min(gfit.rx, gfit.ry))) {
							siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width): aborting.\n"));
							g_free(data);
							return CMD_ARG_ERROR;
						}
						if (!error) {
							data->ks = ks;
						}
					}
				}
				if (error) {
					g_free(data);
					return CMD_ARG_ERROR;
				}
				start_in_new_thread(estimate_only,data);
				return CMD_OK;
			}
		} else if (g_str_has_prefix(arg, "manual")) {
			siril_log_message(_("Manual PSF generation:\n"));
			data->psftype = PSF_MANUAL;
			for (int i = 2; i < nb; i++) {
				char *arg = word[i], *end;
				if (!word[i])
					break;
				if (g_str_has_prefix(arg, "-gaussian")) {
					siril_log_message(_("Gaussian PSF\n"));
					data->profile = 0;
				}
				else if (g_str_has_prefix(arg, "-moffat")) {
					siril_log_message(_("Moffat PSF\n"));
					data->profile = 1;
				}
				else if (g_str_has_prefix(arg, "-disc")) {
					siril_log_message(_("Disc PSF\n"));
					data->profile = 2;
				}
				else if (g_str_has_prefix(arg, "-airy")) {
					siril_log_message(_("Airy disc PSF\n"));
					data->profile = 3;
				}
				else if (g_str_has_prefix(arg, "-ks=")) {
					arg += 4;
					int ks = (int) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((ks < 3) || !(ks %2) || (ks > min(gfit.rx, gfit.ry))) {
						siril_log_message(_("Error in ks parameter: must be odd and between 3 and minimum of (image height, image width): aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->ks = ks;
					}
				}
				else if (g_str_has_prefix(arg, "-fwhm=")) {
					arg += 6;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= 0.f) || (val > 100.f)) {
						siril_log_message(_("Error in fwhm parameter: must be between 0 and 100, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->psf_fwhm = val;
					}
				}
				else if (g_str_has_prefix(arg, "-angle=")) {
					arg += 7;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= -360.f) || (val > 360.f)) {
						siril_log_message(_("Error in ratio parameter: must be between -360 and +360, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->psf_angle = val;
					}
				}
				else if (g_str_has_prefix(arg, "-ratio=")) {
					arg += 7;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val < 1.f) || (val > 5.f)) {
						siril_log_message(_("Error in ratio parameter: must be between 0 and 5, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->psf_ratio = val;
					}
				}
				else if (g_str_has_prefix(arg, "-beta=")) {
					arg += 6;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= 0.f) || (val > 10.f)) {
						siril_log_message(_("Error in beta parameter: must be between 0 and 10, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->psf_beta = val;
					}
				}
				else if (g_str_has_prefix(arg, "-dia=")) {
					arg += 5;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= 0.f) || (val > 5000.f)) {
						siril_log_message(_("Error in dia parameter: must be between 0 and 5000, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->airy_diameter = val;
					}
				}
				else if (g_str_has_prefix(arg, "-fl=")) {
					arg += 4;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= 0.f) || (val > 60000.f)) {
						siril_log_message(_("Error in fl parameter: must be between 0 and 60000, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->airy_fl = val;
					}
				}
				else if (g_str_has_prefix(arg, "-wl=")) {
					arg += 4;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val <= 100.f) || (val > 30000.f)) {
						siril_log_message(_("Error in wl parameter: must be between 100 and 30000, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->airy_wl = val;
					}
				}
				else if (g_str_has_prefix(arg, "-pixelsize=")) {
					arg += 11;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val < 1.f) || (val > 30.f)) {
						siril_log_message(_("Error in fwhm parameter: must be between 1 and 30, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->airy_pixelsize = val;
					}
				}
				else if (g_str_has_prefix(arg, "-obstruct=")) {
					arg += 10;
					float val = (float) g_ascii_strtod(arg, &end);
					if (arg == end) error = TRUE;
					else if ((val < 0.f) || (val >= 100.f)) {
						siril_log_message(_("Error in fwhm parameter: must be between 0 and 100, aborting.\n"));
						g_free(data);
						return CMD_ARG_ERROR;
					}
					if (!error) {
						data->airy_obstruction = val;
					}
				}
			}
			if (error) {
				g_free(data);
				return CMD_ARG_ERROR;
			}
			start_in_new_thread(estimate_only,data);
			return CMD_OK;
		} else if (g_str_has_prefix(arg, "load")) {
			siril_log_message(_("Load PSF from file:\n"));
			if (word[2] && word[2][0] != '\0') {
				if (load_kernel(word[2])) {
					siril_log_message(_("Error loading PSF.\n"));
					g_free(data);
					return CMD_FILE_NOT_FOUND;
				}
			}
			data->psftype = PSF_PREVIOUS;
			g_free(data);
			return CMD_OK;
		}
	}
	g_free(data);
	return CMD_ARG_ERROR;
}

int process_deconvolve(int nb, nonblind_t type) {
	gboolean error = FALSE;
	estk_data* data = malloc(sizeof(estk_data));
	reset_conv_args(data);
	data->regtype = REG_NONE_GRAD;
	if (type == 0)
		data->finaliters = 1;
	if (type == 2)
		data->alpha = 1.f / 500.f;
	for (int i = 1; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(arg, "-alpha=")) {
			arg += 7;
			float alpha = (float) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((alpha < 0.f) || (alpha > 100000.f)) {
				siril_log_message(_("Error in alpha parameter: must be between 0 and 1e5, aborting.\n"));
				g_free(data);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->alpha = alpha;
			}
		}
		else if (g_str_has_prefix(arg, "-iters=")) {
			arg += 7;
			float iters = (int) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((iters < 1) || (iters > 100000)) {
				siril_log_message(_("Error in iterations parameter: must be between 1 and 1e5, aborting.\n"));
				g_free(data);
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
				g_free(data);
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
				g_free(data);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->stepsize = stepsize;
			}
		}
		else if (g_str_has_prefix(arg, "-tv")) {
			 data->regtype = REG_TV_GRAD;
		}
		else if (g_str_has_prefix(arg, "-mul")) {
			data->rl_method = RL_MULT;
		}
		else if (g_str_has_prefix(arg, "-fh")) {
			data->regtype = REG_FH_GRAD;
		}
	}

	if (error) {
		g_free(data);
		return CMD_ARG_ERROR;
	}
	// Guess the user's intentions for a kernel:
	// Order of preference is: if com.stars is populated, assume the user wants to make
	// PSF from stars. If not, if com.kernel is already populated then use it. If not,
	// do a blind deconvolution using the default parameters (and default l0 method).
	// A manual PSF, if required, must be created using the makepsf command and then rl
	// will detect it as an existing kernel and use it.

	if (com.stars && com.stars[0])
		data->psftype = PSF_STARS; // PSF from stars
	else if (com.kernel && (com.kernelsize != 0))
		data->psftype = PSF_PREVIOUS; // use existing kernel
	else data->psftype = PSF_BLIND; // blind deconvolve

	data->nonblindtype = type;

	start_in_new_thread(deconvolve, data);

	return CMD_OK;
}

int process_seqdeconvolve(int nb, nonblind_t type) {
	gboolean error = FALSE;
	estk_data* data = malloc(sizeof(estk_data));
	sequence* seq = NULL;
	reset_conv_args(data);
	data->regtype = REG_NONE_GRAD;
	if (type == 0)
		data->finaliters = 1;
	if (type == 2)
		data->alpha = 1.f / 500.f;
	if (word[1] && word[1][0] != '\0') {
		seq = load_sequence(word[1], NULL);
	}
	if (seq == NULL) {
		siril_log_message(_("Error: cannot open sequence\n"));
		g_free(data);
		return CMD_SEQUENCE_NOT_FOUND;
	}
	for (int i = 2; i < nb; i++) {
		char *arg = word[i], *end;
		if (!word[i])
			break;
		if (g_str_has_prefix(arg, "-alpha=")) {
			arg += 7;
			float alpha = (float) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((alpha < 0.f) || (alpha > 100000.f)) {
				siril_log_message(_("Error in alpha parameter: must be between 0 and 1e5, aborting.\n"));
				g_free(data);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->alpha = alpha;
			}
		}
		else if (g_str_has_prefix(arg, "-iters=")) {
			arg += 7;
			float iters = (int) g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((iters < 1) || (iters > 100000)) {
				siril_log_message(_("Error in iterations parameter: must be between 1 and 1e5, aborting.\n"));
				g_free(data);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
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
				g_free(data);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
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
				g_free(data);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				return CMD_ARG_ERROR;
			}
			if (!error) {
				data->stepsize = stepsize;
			}
		}
		else if (g_str_has_prefix(arg, "-tv")) {
			 data->regtype = REG_TV_GRAD;
		}
		else if (g_str_has_prefix(arg, "-mul")) {
			data->rl_method = RL_MULT;
		}
		else if (g_str_has_prefix(arg, "-fh")) {
			data->regtype = REG_FH_GRAD;
		}
	}

	if (error) {
		g_free(data);
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}
	// Guess the user's intentions for a kernel:
	// Order of preference is: if com.stars is populated, assume the user wants to make
	// PSF from stars. If not, if com.kernel is already populated then use it. If not,
	// do a blind deconvolution using the default parameters (and default l0 method).
	// A manual PSF, if required, must be created using the makepsf command and then rl
	// will detect it as an existing kernel and use it.

	if (com.stars && com.stars[0])
		data->psftype = PSF_STARS; // PSF from stars
	else if (com.kernel && (com.kernelsize != 0))
		data->psftype = PSF_PREVIOUS; // use existing kernel
	else data->psftype = PSF_BLIND; // blind deconvolve

	data->nonblindtype = type;

	deconvolve_sequence_command(data, seq);

	return CMD_OK;
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
	unsharp(&(gfit), sigma, multi, TRUE);

	char log[90];
	sprintf(log, "Unsharp filtering, sigma: %.2f, coefficient: %.2f", sigma, multi);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
	return CMD_OK;
}

static gboolean crop_command_idle(gpointer arg) {
	// operations that are not in the generic idle of notify_gfit_modified()
	clear_stars_list(TRUE);
	delete_selected_area();
	reset_display_offset();
	update_zoom_label();
	return FALSE;
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
			if (area.x + area.w > gfit.rx || area.y + area.h > gfit.ry) {
				siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"), gfit.rx, gfit.ry);
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
			return CMD_ARG_ERROR;
		}
	}

	crop(&gfit, &area);

	char log[90];
	sprintf(log, _("Crop (x=%d, y=%d, w=%d, h=%d)"),
					area.x, area.y, area.w, area.h);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	siril_add_idle(crop_command_idle, NULL);

	return CMD_OK;
}

int process_cd(int nb) {
	long maxpath = get_pathmax();
	char filename[maxpath];
	int retval;

	g_strlcpy(filename, word[1], maxpath - 1);

	expand_home_in_filename(filename, maxpath);
	retval = siril_change_dir(filename, NULL);
	if (!retval && !com.script) {
		if (com.pref.wd)
			g_free(com.pref.wd);
		com.pref.wd = g_strdup(com.wd);
		writeinitfile();
		set_GUI_CWD();
	}
	return retval;
}

int process_wrecons(int nb) {

	float coef[7];
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];

	int nb_chan = gfit.naxes[2];

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
		if (gfit.type == DATA_USHORT) {
			wavelet_reconstruct_file(dir[i], coef, gfit.pdata[i]);
		} else if (gfit.type == DATA_FLOAT) {
			wavelet_reconstruct_file_float(dir[i], coef, gfit.fpdata[i]);
		}
		else return CMD_GENERIC_ERROR;
		g_free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_linstretch(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;

	if (nb <= 0) return CMD_ARG_ERROR;

	double BP = 0.0;
	gchar *end;

	BP = g_ascii_strtod(word[1], &end);
	if (end == word[1] || (BP < 0.0) || (BP > 1.0)) {
		siril_log_message(_("Black Point BP must be between 0 and 1\n"));
	}
	if (word[2]) {
		if (!strcmp(word[2], "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[2], "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[2], "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		if (!strcmp(word[2], "RG")) {
			do_blue = FALSE;
		}
		if (!strcmp(word[2], "RB")) {
			do_green = FALSE;
		}
		if (!strcmp(word[2], "GB")) {
			do_red = FALSE;
		}
	}
	set_cursor_waiting(TRUE);
	ght_params params = {0.0, 0.0, 0.0, 0.0, 0.0, BP, STRETCH_LINEAR, COL_INDEP, do_red, do_green, do_blue};
	apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

	char log[90];
	sprintf(log, "Linear stretch to black point %.3f", BP);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_ght(int nb) {
	int stretch_colourmodel = COL_INDEP;
	int arg_offset = 0;
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
	gchar *end;

	if (!strcmp(word[1], "-human")) {
		stretch_colourmodel = COL_HUMANLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-even")) {
		stretch_colourmodel = COL_EVENLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-independent")) {
		stretch_colourmodel = COL_INDEP;
		arg_offset = 1;
	}
	if (nb <= 4 + arg_offset)
	return CMD_WRONG_N_ARG;

	double D = g_ascii_strtod(word[1 + arg_offset], &end);
	if (end == word[1 + arg_offset] || (D < 0.0) || (D > 10.0)) {
		siril_log_message(_("Stretch factor D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double B = g_ascii_strtod(word[2 + arg_offset], &end);
	if (end == word[2 + arg_offset] || (B < -5.0) || (B > 15.0)) {
		siril_log_message(_("Stretch intensity B must be between -5 and +15\n"));
		return CMD_ARG_ERROR;
	}

	double SP = g_ascii_strtod(word[4 + arg_offset], &end);
	if (end == word[4 + arg_offset] || (SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("Stretch focal point SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[3 + arg_offset], &end);
	if (end == word[3 + arg_offset] || (LP < 0.0) || (LP > SP)) {
		siril_log_message(_("Shadow preservation point LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[5 + arg_offset], &end);
	if (end == word[5 + arg_offset] || (HP < SP) || (HP > 1.0)) {
		siril_log_message(_("Headroom preservation point HP must be between stretch focal point and 1\n"));
		return CMD_ARG_ERROR;
	}

	if (word[6 + arg_offset]) {
		if (!strcmp(word[6 + arg_offset], "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "RG")) {
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "RB")) {
			do_green = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "GB")) {
			do_red = FALSE;
		}
	}
	set_cursor_waiting(TRUE);
	ght_params params = {B, D, LP, SP, HP, 0.0, STRETCH_PAYNE_NORMAL, stretch_colourmodel, do_red, do_green, do_blue};
	apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

	char log[100];
	sprintf(log, "GHS (pivot: %.3f, amount: %.2f, local: %.1f [%.2f, %.2f])", SP, D, B, LP, HP);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_invght(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
	gchar *end;

	int stretch_colourmodel = COL_INDEP;
	int arg_offset = 0;
	if (!strcmp(word[1], "-human")) {
		stretch_colourmodel = COL_HUMANLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-even")) {
		stretch_colourmodel = COL_EVENLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-independent")) {
		stretch_colourmodel = COL_INDEP;
		arg_offset = 1;
	}
	if (nb <= 4 + arg_offset)

		return CMD_WRONG_N_ARG;

	double D = g_ascii_strtod(word[1 + arg_offset], &end);
	if (end == word[1 + arg_offset] || (D < 0.0) || (D > 10.0)) {
		siril_log_message(_("Stretch factor D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double B = g_ascii_strtod(word[2 + arg_offset], &end);
	if (end == word[2 + arg_offset] || (B < -5.0) || (B > 15.0)) {
		siril_log_message(_("Stretch intensity B must be between -5 and +15\n"));
		return CMD_ARG_ERROR;
	}

	double SP = g_ascii_strtod(word[4 + arg_offset], &end);
	if (end == word[4 + arg_offset] || (SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("Stretch focal point SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[3 + arg_offset], &end);
	if (end == word[3 + arg_offset] || (LP < 0.0) || (LP > SP)) {
		siril_log_message(_("Shadow preservation point LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[5 + arg_offset], &end);
	if (end == word[5 + arg_offset] || (HP < SP) || (HP > 1.0)) {
		siril_log_message(_("Headroom preservation point HP must be between stretch focal point and 1\n"));
		return CMD_ARG_ERROR;
	}
	if (word[6 + arg_offset]) {
		if (!strcmp(word[6 + arg_offset], "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "RG")) {
			do_blue = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "RB")) {
			do_green = FALSE;
		}
		if (!strcmp(word[6 + arg_offset], "GB")) {
			do_red = FALSE;
		}
	}

	set_cursor_waiting(TRUE);
	ght_params params = {B, D, LP, SP, HP, 0.0, STRETCH_PAYNE_INVERSE,
		stretch_colourmodel, do_red, do_green, do_blue};
	apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

	char log[100];
	sprintf(log, "Inverse GHS (pivot: %.3f, amount: %.2f, local: %.1f [%.2f, %.2f])",
			SP, D, B, LP, HP);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_modasinh(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
	gchar *end;

	int stretch_colourmodel = COL_INDEP;
	int arg_offset = 0;
	if (!strcmp(word[1], "-human")) {
		stretch_colourmodel = COL_HUMANLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-even")) {
		stretch_colourmodel = COL_EVENLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-independent")) {
		stretch_colourmodel = COL_INDEP;
		arg_offset = 1;
	}
	if (nb <= arg_offset + 3)
	return CMD_WRONG_N_ARG;

	double D = g_ascii_strtod(word[arg_offset + 1], &end);
	if (end == word[1 + arg_offset] || (D < 0.0) || (D > 10.0)) {
		siril_log_message(_("D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double SP = g_ascii_strtod(word[arg_offset + 3], &end);
	if (end == word[3 + arg_offset] || (SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[arg_offset + 2], &end);
	if (end == word[2 + arg_offset] || (LP < 0.0) || (LP > SP)) {
		siril_log_message(_("LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[arg_offset + 4], &end);
	if (end == word[4 + arg_offset] || (HP < SP) || (HP > 1.0)) {
		siril_log_message(_("HP must be between stretch focal point and 1\n"));
		return CMD_ARG_ERROR;
	}
	if (word[5 + arg_offset]) {
		if (!strcmp(word[5 + arg_offset], "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "RG")) {
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "RB")) {
			do_green = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "GB")) {
			do_red = FALSE;
		}
	}

	set_cursor_waiting(TRUE);
	ght_params params = {0.0, D, LP, SP, HP, 0.0, STRETCH_ASINH, stretch_colourmodel, do_red, do_green, do_blue};
	apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

	char log[100];
	sprintf(log, "Modified asinh (pivot: %.3f, amount: %.2f, [%.2f, %.2f])",
			SP, D, LP, HP);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_invmodasinh(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
	gchar *end;

	int stretch_colourmodel = COL_INDEP;
	int arg_offset = 0;
	if (!strcmp(word[1], "-human")) {
		stretch_colourmodel = COL_HUMANLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-even")) {
		stretch_colourmodel = COL_EVENLUM;
		arg_offset = 1;
	}
	else if (!strcmp(word[1], "-independent")) {
		stretch_colourmodel = COL_INDEP;
		arg_offset = 1;
	}
	if (nb <= arg_offset + 3)
	return CMD_WRONG_N_ARG;

	double D = g_ascii_strtod(word[arg_offset + 1], &end);
	if (end == word[1 + arg_offset] || (D < 0.0) || (D > 10.0)) {
		siril_log_message(_("D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float)D);

	double SP = g_ascii_strtod(word[arg_offset + 3], &end);
	if (end == word[3 + arg_offset] || (SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[arg_offset + 2], &end);
	if (end == word[2 + arg_offset] || (LP < 0.0) || (LP > SP)) {
		siril_log_message(_("LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[arg_offset + 4], &end);
	if (end == word[4 + arg_offset] || (HP < SP) || (HP > 1.0)) {
		siril_log_message(_("HP must be between stretch focal point and 1\n"));
		return CMD_ARG_ERROR;
	}
	if (word[5 + arg_offset]) {
		if (!strcmp(word[5 + arg_offset], "R")) {
			do_green = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "G")) {
			do_red = FALSE;
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "B")) {
			do_green = FALSE;
			do_red = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "RG")) {
			do_blue = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "RB")) {
			do_green = FALSE;
		}
		if (!strcmp(word[5 + arg_offset], "GB")) {
			do_red = FALSE;
		}
	}

	set_cursor_waiting(TRUE);
	ght_params params = {0.0, D, LP, SP, HP, 0.0, STRETCH_INVASINH, stretch_colourmodel, do_red, do_green, do_blue};
	apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

	char log[100];
	sprintf(log, "Inverse modified asinh (pivot: %.3f, amount: %.2f, [%.2f, %.2f])",
			SP, D, LP, HP);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_wavelet(int nb) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];

	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;
	const char* tmpdir = g_get_tmp_dir();
	gchar *end1, *end2;

	Nbr_Plan = g_ascii_strtoull(word[1], &end1, 10);
	Type_Transform = g_ascii_strtoull(word[2], &end2, 10);

	nb_chan = gfit.naxes[2];
	g_assert(nb_chan <= 3);

	mins = min (gfit.rx, gfit.ry);
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

	if (gfit.type == DATA_USHORT) {
		float *Imag = f_vector_alloc(gfit.rx * gfit.ry);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return CMD_ALLOC_ERROR;
		}

		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan], NULL);
			wavelet_transform_file(Imag, gfit.ry, gfit.rx, dir[chan],
					Type_Transform, Nbr_Plan, gfit.pdata[chan]);
			g_free(dir[chan]);
		}

		free(Imag);
	} else if (gfit.type == DATA_FLOAT) {
		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan], NULL);
			wavelet_transform_file_float(gfit.fpdata[chan], gfit.ry, gfit.rx, dir[chan],
					Type_Transform, Nbr_Plan);
			g_free(dir[chan]);
		}
	}
	else return CMD_INVALID_IMAGE;
	return CMD_OK;
}

int process_log(int nb){
	loglut(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
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

	if (readfits(word[1], &ref, NULL, gfit.type == DATA_FLOAT))
		return CMD_INVALID_IMAGE;
	if (!find_linear_coeff(&gfit, &ref, low, high, a, b, NULL)) {
		set_cursor_waiting(TRUE);
		apply_linear_to_fits(&gfit, a, b);

		adjust_cutoff_from_updated_gfit();
		redraw(REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
	clearfits(&ref);
	return CMD_OK;
}

int process_asinh(int nb) {
	gboolean human_luminance = FALSE;
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
		siril_log_message(_("Stretch must be greater than or equal to 1\n"));
		return CMD_ARG_ERROR;
	}
	arg_offset++;

	double offset = 0.0;
	if (nb > arg_offset) {
		offset = g_ascii_strtod(word[arg_offset], &end);
		if (end == word[arg_offset]) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[arg_offset]);
			return CMD_ARG_ERROR;
		}
	}

	set_cursor_waiting(TRUE);
	asinhlut(&gfit, beta, offset, human_luminance);

	char log[90];
	sprintf(log, "Asinh stretch (amount: %.1f, offset: %.1f, human: %s)", beta, offset, human_luminance ? "yes" : "no");
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
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

	if (end == word[2] || size <= 0.0) {
		siril_log_message(_("Tile size must be > 0.\n"));
		return CMD_ARG_ERROR;
	}

	struct CLAHE_data *args = malloc(sizeof(struct CLAHE_data));

	args->fit = &gfit;
	args->clip = clip_limit;
	args->tileSize = size;

	start_in_new_thread(clahe, args);

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
				filename[maxpath - 1] = '\0';
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
	struct dirent **list;

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
			retval = 1;
			free(seqpath1); free(seqpath2);	g_free(seqname); g_free(dir);
			goto merge_clean_up;
		}
		g_free(seqname);
		if (seq_check_basic_data(seqs[i], FALSE) < 0) {
			siril_log_message(_("Sequence `%s' is invalid, could not merge\n"), word[i + 1]);
			retval = 1;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (i != 0 && (seqs[i]->rx != seqs[0]->rx ||
					seqs[i]->ry != seqs[0]->ry ||
					seqs[i]->nb_layers != seqs[0]->nb_layers ||
					seqs[i]->bitpix != seqs[0]->bitpix ||
					seqs[i]->type != seqs[0]->type)) {
			siril_log_message(_("All sequences must be the same format for merging. Sequence `%s' is different\n"), word[i + 1]);
			retval = 1;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (seqs[i]->type == SEQ_REGULAR) {
			// we need to build the list of files
			long maxpath = get_pathmax();
			char filename[maxpath];
			for (int image = 0; image < seqs[i]->number; image++) {
				fit_sequence_get_image_filename(seqs[i], image, filename, TRUE);
				list = g_list_append(list, g_build_filename(dir, filename, NULL));
			}
		}
		free(seqpath1); free(seqpath2); g_free(dir);
		siril_change_dir(dest_dir, NULL);	// they're all relative to this one
	}

	struct ser_struct out_ser;
	struct _convert_data *args = NULL;
	fitseq out_fitseq;
	switch (seqs[0]->type) {
		case SEQ_REGULAR:
			// use the conversion, it makes symbolic links or copies as a fallback
			args = malloc(sizeof(struct _convert_data));
			args->start = 0;
			args->total = 0; // init to get it from glist_to_array()
			args->list = glist_to_array(list, &args->total);
			args->destroot = format_basename(word[nb - 1], FALSE);
			args->input_has_a_seq = FALSE;
			args->input_has_a_film = FALSE;
			args->debayer = FALSE;
			args->multiple_output = FALSE;
			args->output_type = SEQ_REGULAR;
			args->make_link = TRUE;
			gettimeofday(&(args->t_start), NULL);
			start_in_new_thread(convert_thread_worker, args);
			break;

		case SEQ_SER:
			if (g_str_has_suffix(word[nb - 1], ".ser"))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s.ser", word[nb - 1]);
			if (ser_create_file(outseq_name, &out_ser, TRUE, seqs[0]->ser_file)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = 1;
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
						retval = 1;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}

					if (ser_write_frame_from_fit(&out_ser, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = 1;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (ser_write_and_close(&out_ser)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = 1;
			}
			break;

		case SEQ_FITSEQ:
			if (g_str_has_suffix(word[nb - 1], com.pref.ext))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s%s", word[nb - 1], com.pref.ext);
			if (fitseq_create_file(outseq_name, &out_fitseq, -1)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = 1;
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
						retval = 1;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}

					if (fitseq_write_image(&out_fitseq, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = 1;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (fitseq_close_file(&out_fitseq)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = 1;
			}
			break;
		default:
			siril_log_message(_("This type of sequence cannot be created by Siril, aborting the merge\n"));
			retval = 1;
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
	if (!strcmp(fit.row_order, "BOTTOM-UP")) {
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
		if (!strcmp(gfit.row_order, "BOTTOM-UP")) {
			siril_log_message(_("Image data is already bottom-up\n"));
			return CMD_OK;
		}
		siril_log_message(_("Mirroring image to convert to bottom-up data\n"));
	}
	mirrorx(&gfit, TRUE);
	notify_gfit_modified();
	return CMD_OK;
}

int process_mirrory(int nb){
	mirrory(&gfit, TRUE);
	notify_gfit_modified();
	return CMD_OK;
}

int process_mtf(int nb) {
	struct mtf_params params;
	gboolean inverse = word[0][0] == 'i' || word[0][0] == 'I';
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
	if (word[4]) {
		if (!strcmp(word[4], "R")) {
			params.do_green = FALSE;
			params.do_blue = FALSE;
		}
		if (!strcmp(word[4], "G")) {
			params.do_red = FALSE;
			params.do_blue = FALSE;
		}
		if (!strcmp(word[4], "B")) {
			params.do_green = FALSE;
			params.do_red = FALSE;
		}
		if (!strcmp(word[4], "RG")) {
			params.do_blue = FALSE;
		}
		if (!strcmp(word[4], "RB")) {
			params.do_green = FALSE;
		}
		if (!strcmp(word[4], "GB")) {
			params.do_red = FALSE;
		}
	}

	if (inverse)
		apply_linked_pseudoinverse_mtf_to_fits(&gfit, &gfit, params, TRUE);
	else apply_linked_mtf_to_fits(&gfit, &gfit, params, TRUE);

	char log[90];
	sprintf(log, "%s transfer (%.3f, %.4f, %.3f)",
			inverse ? "Inverse midtones" : "Midtones",
			params.shadows, params.midtones, params.highlights);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_autoghs(int nb) {
	int argidx = 1;
	gboolean linked = FALSE;
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
		argidx++;
	}

	int nb_channels = (int)gfit.naxes[2];
	imstats *stats[3] = { NULL };
	int ret = compute_all_channels_statistics_single_image(&gfit, STATS_BASIC, MULTI_THREADED, stats);
	if (ret) {
		for (int i = 0; i < nb_channels; ++i) {
			free_stats(stats[i]);
		}
		return CMD_GENERIC_ERROR;
	}
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
		if (gfit.type == DATA_USHORT)
			SP *= (gfit.orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE;
		siril_log_message(_("Symmetry point SP=%f\n"), SP);

		ght_params params = { .B = b, .D = amount, .LP = lp, .SP = SP, .HP = hp,
			.BP = 0.0, STRETCH_PAYNE_NORMAL, COL_INDEP, TRUE, TRUE, TRUE };
		apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(nb_channels > 1)
#endif
		for (int i = 0; i < nb_channels; ++i) {
			gboolean do_red = i == 0, do_green = i == 1, do_blue = i == 2;
			if (stats[i]) {
				float SP = stats[i]->median + shadows_clipping * stats[i]->sigma;
				if (gfit.type == DATA_USHORT)
					SP *= (gfit.orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE;
				siril_log_message(_("Symmetry point for channel %d: SP=%f\n"), i, SP);

				ght_params params = { .B = b, .D = amount, .LP = lp, .SP = SP, .HP = hp,
					.BP = 0.0, STRETCH_PAYNE_NORMAL, COL_INDEP, do_red, do_green, do_blue};
				apply_linked_ght_to_fits(&gfit, &gfit, params, TRUE);

				free_stats(stats[i]);
			}
		}
	}

	char log[100];
	sprintf(log, "AutoGHS (%sk.sigma: %.2f, amount: %.2f, local: %.1f [%.2f, %.2f])",
			linked ? "linked, " : "", shadows_clipping, amount, b, lp, hp);
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	return CMD_OK;
}

int process_autostretch(int nb) {
	int arg_index = 1;
	gboolean linked = FALSE;
	if (nb > 1 && !strcmp(word[1], "-linked")) {
		linked = TRUE;
		arg_index++;
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
	if (linked) {
		struct mtf_params params;
		find_linked_midtones_balance(&gfit, shadows_clipping, target_bg, &params);
		params.do_red = TRUE;
		params.do_green = TRUE;
		params.do_blue = TRUE;
		apply_linked_mtf_to_fits(&gfit, &gfit, params, TRUE);
	} else {
		struct mtf_params params[3];
		find_unlinked_midtones_balance(&gfit, shadows_clipping, target_bg, params);
		apply_unlinked_mtf_to_fits(&gfit, &gfit, params);
	}

	char log[90];
	sprintf(log, "Autostretch (shadows: %.2f, target bg: %.2f, %s)",
			shadows_clipping, target_bg, linked ? "linked" : "unlinked");
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
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

	fits_binning(&gfit, factor, mean);

	notify_gfit_modified();
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
			toX = round_to_int(gfit.rx * (double)toY / gfit.ry);
		} else if (g_str_has_prefix(word[1], "-width=")) {
			toX = g_ascii_strtoull(word[1] + 7, &end, 10);
			if (end == word[1]+8) {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				return CMD_ARG_ERROR;
			}
			toY = round_to_int(gfit.ry * (double)toX / gfit.rx);
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[1]);
			return CMD_ARG_ERROR;
		}
	} else {
		double factor = g_ascii_strtod(word[1], &end);
		if (end == word[1] || factor <= 0.0 || factor > 5.0) {
			siril_log_message(_("Scale %lf not allowed. Should be between 0.0 and 5.0.\n"), factor);
			return CMD_ARG_ERROR;
		}
		if (factor == 1.0) {
			siril_log_message(_("Scale is 1.0. Not doing anything.\n"));
			return CMD_ARG_ERROR;
		}
		toX = round_to_int(factor * gfit.rx);
		toY = round_to_int(factor * gfit.ry);
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
	siril_log_message(_("Resampling to %d x %d pixels with %s interpolation\n"),
			toX, toY, interp_to_str(interpolation));
	int fromX = gfit.rx, fromY = gfit.ry;

	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, interpolation, clamp);

	char log[90];
	sprintf(log, "Resampled from %d x %d, %s interp%s", fromX, fromY, interp_to_str(interpolation),
			((interpolation == OPENCV_LANCZOS4 || interpolation == OPENCV_CUBIC) && clamp) ?
			", clamped" : "");
	gfit.history = g_slist_append(gfit.history, strdup(log));

	notify_gfit_modified();
	if (!com.script) update_MenuItem();
	return CMD_OK;
}

int process_rgradient(int nb) {
	if (gfit.orig_bitpix == BYTE_IMG) {
		siril_log_color_message(_("This process cannot be applied to 8b images\n"), "red");
		return CMD_INVALID_IMAGE;
	}

	gchar *endxc, *endyc, *enddR, *endda;

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = g_ascii_strtod(word[1], &endxc);
	args->yc = g_ascii_strtod(word[2], &endyc);
	args->dR = g_ascii_strtod(word[3], &enddR);
	args->da = g_ascii_strtod(word[4], &endda);
	args->fit = &gfit;

	if (endxc == word[1] || endyc == word[2] || enddR == word[3]
			|| endda == word[4] || (args->xc >= args->fit->rx)
			|| (args->yc >= args->fit->ry)) {
		siril_log_message(_("The coordinates cannot be greater than the size of the image. "
				"Please change their values and retry.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}

	start_in_new_thread(rgradient_filter, args);
	return CMD_OK;
}

int process_rotate(int nb) {
	gchar *end;
	set_cursor_waiting(TRUE);
	int crop = 1;
	gboolean has_selection = FALSE;
	rectangle area = { 0, 0, gfit.rx, gfit.ry };
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
	siril_debug_print("%f\n",angle);

	if (angle == 0.0 || angle == 360.0) {
		siril_log_message(_("Angle is 0.0 or 360.0 degrees. Doing nothing...\n"));
		return CMD_ARG_ERROR;
	}
	verbose_rotate_image(&gfit, area, angle, interpolation, crop, clamp);

	// the new selection will match the current image
	if (has_selection) {
		com.selection = (rectangle){ 0, 0, gfit.rx, gfit.ry };
		new_selection_zone();
	}
	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

int process_rotatepi(int nb){
	if (verbose_rotate_fast(&gfit, 180))
		return CMD_GENERIC_ERROR;

	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
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
		return 1;
	}
	input[sep] = '\0';
	if (is_get) {
		print_settings_key(input, input+sep+1, FALSE);
	} else {
		/* set */
		long pathmax = get_pathmax();
		char fakefile[pathmax];
		int filelen = snprintf(fakefile, pathmax - 1, "[%s]\n%s\n", input, input+sep+1);
		GKeyFile *kf = g_key_file_new();
		g_key_file_load_from_data(kf, fakefile, filelen, G_KEY_FILE_NONE, NULL);
		return read_keyfile(kf);
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
		mag = gui.qphot->mag;
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
		struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
		psf_star *result = psf_get_minimisation(&gfit, select_vport(gui.cvport), &com.selection, TRUE, ps, TRUE, com.pref.starfinder_conf.profile, &error);
		free(ps);
		if (result) {
			found = TRUE;
			mag = result->mag;
			free_psf(result);
		}
		else siril_log_message(_("PSF minimisation failed with error %d\n"), error);
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
		double inner = -1.0, outer = -1.0, aperture = -1.0, gain = -1.0;
		int min = -65536, max = -1, force = -1;
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
			else if (g_str_has_prefix(arg, "-force_radius=")) {
				arg += 14;
				if (*arg == 'y')
					force = 1;
				else if (*arg == 'n')
					force = 0;
				else error = TRUE;
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
		if (force == 0)
			com.pref.phot_set.force_radius = FALSE;
		else if (force == 1)
			com.pref.phot_set.force_radius = TRUE;
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
	siril_log_message(_("Local background annulus radius: %.1f to %.1f, aperture: %.1f (%s), camera conversion gain: %f e-/ADU, using pixels with values ]%f, %f[\n"),
			com.pref.phot_set.inner,
			com.pref.phot_set.outer,
			com.pref.phot_set.aperture,
			com.pref.phot_set.force_radius ? _("forced") : _("unused, dynamic"),
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
						siril_log_message(_("Wrong parameter values. Auto must be set to on or off, aborting.\n"));
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
				} else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
					return CMD_ARG_ERROR;
				}
			}
		}
	}

	siril_log_message(_("sigma = %3.2f\n"), sigma);
	com.pref.starfinder_conf.sigma = sigma;
	siril_log_message(_("roundness = %3.2f\n"), roundness);
	com.pref.starfinder_conf.roundness = roundness;
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

static void remove_char_from_str(char *str, const char toRemove) {
	int i, j;
	int len = strlen(str);

	for (i = 0; i < len; i++) {
		/*
		 * If the character to remove is found then shift all characters to one
		 * place left and decrement the length of string by 1.
		 */
		if (str[i] == toRemove) {
			for (j = i; j < len; j++) {
				str[j] = str[j + 1];
			}

			len--;
			// If a character is removed then make sure i doesn't increments
			i--;
		}
	}
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

	cur = expression;
	while ((next = strchr(cur, '$')) != NULL) {
		cur = next + 1;
		count++;
	}

	if (count == 0) {
		siril_log_message(_("You need to add at least one image as variable. Use $ tokens to surround the file names .\n"));
		return CMD_ARG_ERROR;
	} else if (count % 2 != 0) {
		siril_log_message(_("There is an unmatched $. Please check the expression.\n"));
		return CMD_ARG_ERROR;
	}

	/* parse rescale option if exist */
	if (nb > 1) {
		if (!g_strcmp0(word[2], "-rescale")) {
			if (nb == 5) {
				gchar *end;
				min = g_ascii_strtod(word[3], &end);
				if (end == word[3] || min < 0 || min > 1) {
					siril_log_message(_("Rescale can only be done in the [0, 1] range.\n"));
					return CMD_ARG_ERROR;
				}
				max = g_ascii_strtod(word[4], &end);
				if (end == word[4] || max < 0 || max > 1) {
					siril_log_message(_("Rescale can only be done in the [0, 1] range.\n"));
					return CMD_ARG_ERROR;
				}
			} else {
				min = 0.f;
				max = 1.f;
			}
		}
	}

	struct pixel_math_data *args = malloc(sizeof(struct pixel_math_data));
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
			args->varname[i] = g_strdup(start + 1);
			start = cur = end;
			i++;
		}
		cur++;
	}

	int width = -1;
	int height = -1;
	int channel = -1;

	for (int j = 0; j < args->nb_rows; j++) {
		int w, h, c;
		if (load_pm_var(args->varname[j], j, &w, &h, &c)) {
			if (j > 0)
				free_pm_var(j - 1);
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
				free(args->varname);
				free(args);
				return CMD_INVALID_IMAGE;
			}
		}
	}

	/* remove tokens */
	g_free(expression);
	expression = g_shell_unquote(word[1], NULL);
	remove_char_from_str(expression, '$');
	remove_spaces_from_str(expression);

	fits *fit = NULL;
	if (new_fit_image(&fit, width, height, channel, DATA_FLOAT)) {
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
	if (min >= 0.f) {
		args->rescale = TRUE;
		args->min = min;
		args->max = max;
	} else {
		args->rescale = FALSE;
	}

	start_in_new_thread(apply_pixel_math_operation, args);

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

	if (gfit.naxes[2] > 1 && nb == 1 && com.headless) {
		siril_log_color_message(_("Please display the channel on which you want to compute the PSF or use -channel argument\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	int channel = com.headless ? 0 : select_vport(gui.cvport);
	if (nb == 2) {
		gchar *next;
		channel = g_ascii_strtoull(word[1], &next, 10);
		if (word[1] == next || channel > (int)gfit.naxes[2]) {
			siril_log_message(_("Please provide the channel number starting from 0 for red\n"));
			return CMD_ARG_ERROR;
		}
	}

	starprofile profile = com.pref.starfinder_conf.profile;
	psf_error error;
	struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
	psf_star *result = psf_get_minimisation(&gfit, channel, &com.selection, TRUE, ps, TRUE, profile, &error);
	free(ps);
	if (result) {
		psf_display_result(result, &com.selection);
		free_psf(result);
	}
	else siril_log_message(_("PSF minimisation failed with error %d\n"), error);
	return CMD_OK;
}

int process_seq_tilt(int nb) {
	gboolean draw_polygon = FALSE;

	sequence *seq;

	if (word[1] && word[1][0] != '\0') {
		seq = load_sequence(word[1], NULL);
	} else {
		if (!sequence_is_loaded()) {
			return CMD_NOT_FOR_SINGLE;
		}
		seq = &com.seq;
		draw_polygon = TRUE;
	}
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;

	struct tilt_data *args = calloc(sizeof(struct tilt_data), 1);
	args->seq = seq;
	args->draw_polygon = draw_polygon;

	apply_tilt_to_sequence(args);
	return CMD_OK;
}

int parse_star_position_arg(char *arg, sequence *seq, fits *first, rectangle *seqpsf_area, gchar **target_descr) {
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
		double x, y;
		if (wcs2pix(first, ra, dec, &x, &y)) {
			siril_log_message(_("The given coordinates are not in the image, aborting\n"));
			return CMD_ARG_ERROR;
		}
		y = first->ry - y - 1;
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
		siril_log_message(_("Coordinates of the star: %.1f, %.1f\n"), x, y);
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

/* seqpsf [sequencename channel { -at=x,y | -wcs=ra,dec }] */
int process_seq_psf(int nb) {
	if (com.script && nb < 4) {
		siril_log_message(_("Arguments are not optional when called from a script\n"));
		return CMD_ARG_ERROR;
	}
	if (!com.headless && !sequence_is_loaded() && nb < 4) {
		siril_log_message(_("Arguments are optional only if a sequence is already loaded and selection around a star made\n"));
		return CMD_ARG_ERROR;
	}

	sequence *seq = NULL;
	int layer;
	if (nb < 4) {
		seq = &com.seq;
		layer = select_vport(gui.cvport);

		if (com.selection.w > 300 || com.selection.h > 300) {
			siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
			return CMD_SELECTION_ERROR;
		}
		if (com.selection.w <= 0 || com.selection.h <= 0) {
			siril_log_message(_("Select an area first\n"));
			return CMD_SELECTION_ERROR;
		}
	} else {
		seq = load_sequence(word[1], NULL);
		if (!seq) {
			return CMD_SEQUENCE_NOT_FOUND;
		}
		fits first = { 0 };
		if (seq_read_frame_metadata(seq, 0, &first)) {
			free_sequence(seq, TRUE);
			return CMD_GENERIC_ERROR;
		}
		seq->current = 0;

		gchar *end;
		layer = g_ascii_strtoull(word[2], &end, 10);
		if (end == word[2] || layer >= seq->nb_layers) {
			siril_log_message(_("PSF cannot be computed on channel %d for this sequence of %d channels\n"), layer, seq->nb_layers);
			free_sequence(seq, TRUE);
			return CMD_ARG_ERROR;
		}

		rectangle area;
		if (parse_star_position_arg(word[3], seq, &first, &area, NULL)) {
			free_sequence(seq, TRUE);
			return 1;
		}
		com.selection = area;
	}

	framing_mode framing = REGISTERED_FRAME;
	if (framing == REGISTERED_FRAME && !seq->regparam[layer])
		framing = ORIGINAL_FRAME;
	if (framing == ORIGINAL_FRAME) {
		if (com.headless)
			framing = FOLLOW_STAR_FRAME;
		else {
			GtkToggleButton *follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
			if (gtk_toggle_button_get_active(follow))
				framing = FOLLOW_STAR_FRAME;
		}
	}
	siril_log_message(_("Running the PSF on the sequence, layer %d\n"), layer);
	return seqpsf(seq, layer, FALSE, FALSE, framing, TRUE, com.script);
}

// light_curve sequencename channel { -ninalist=file | [-auto] { -at=x,y | -wcs=ra,dec [-refat=x,y] [-refwcs=ra,dec] ... } }
int process_light_curve(int nb) {
	sequence *seq;
	int layer;

	seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	siril_debug_print("reference image in seqfile is %d\n", seq->reference_image);
	int refimage = 0;
	gboolean seq_has_wcs;
	if ((seq_has_wcs = sequence_has_wcs(seq, &refimage)))
		seq->reference_image = refimage;
	seq->current = refimage;	// seqpsf computes transformations from current

	gchar *end;
	layer = g_ascii_strtoull(word[2], &end, 10);
	if (end == word[2] || layer >= seq->nb_layers) {
		siril_log_message(_("PSF cannot be computed on channel %d for this sequence of %d channels\n"), layer, seq->nb_layers);
		free_sequence(seq, TRUE);
		return CMD_ARG_ERROR;
	}

	if (sequence_drifts(seq, layer, seq->rx / 4))
		return CMD_GENERIC_ERROR;

	struct light_curve_args *args = malloc(sizeof(struct light_curve_args));

	/* we have a sequence and channel, now we iterate with seqpsf on the list of stars */
	if (g_str_has_prefix(word[3], "-ninastars=")) {
		char *file = word[3] + 11;
		if (file[0] == '\0') {
			siril_log_message(_("Missing argument to %s, aborting.\n"), word[3]);
			free(args);
			return 1;
		}
		if (!seq_has_wcs) {
			siril_log_message(_("No image in the sequence was found with the WCS information required for star selection by equatorial coordinates, plate solve the reference or the first\n"));
			free(args);
			return CMD_FOR_PLATE_SOLVED;
		}

		fits first = { 0 };
		if (seq_read_frame_metadata(seq, refimage, &first)) {
			free_sequence(seq, TRUE);
			free(args);
			return CMD_GENERIC_ERROR;
		}
		if (parse_nina_stars_file_using_WCS(args, word[3]+11, TRUE, TRUE, &first)) {
			free_sequence(seq, TRUE);
			free(args);
			return CMD_GENERIC_ERROR;
		}
	} else {
		fits first = { 0 };
		if (seq_read_frame_metadata(seq, refimage, &first)) {
			free_sequence(seq, TRUE);
			free(args);
			return CMD_GENERIC_ERROR;
		}
		args->areas = malloc((nb - 3) * sizeof(rectangle));
		for (int arg_index = 3; arg_index < nb; arg_index++) {
			if (parse_star_position_arg(word[arg_index], seq, &first, &args->areas[arg_index - 3], &args->target_descr)) {
				free_sequence(seq, TRUE);
				free(args->areas);
				free(args);
				return CMD_ARG_ERROR;;
			}
		}
		args->nb = nb - 3;
	}
	args->seq = seq;
	args->layer = layer;
	args->display_graph = FALSE;
	siril_debug_print("starting PSF analysis of %d stars\n", args->nb);

	// TODO: display stars if sequence is current and not headless

	start_in_new_thread(light_curve_worker, args);
	return 0;
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

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	args->seq = seq;
	args->area = area;
	args->prefix = "cropped_";

	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), current);
						g_free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					args->prefix = strdup(value);
				}
			}
		}
	}

	crop_sequence(args);
	return CMD_OK;
}

int process_bg(int nb) {
	WORD us_bg;

	for (int layer = 0; layer < gfit.naxes[2]; layer++) {
		double bg = background(&gfit, layer, &com.selection, MULTI_THREADED);
		if (gfit.type == DATA_USHORT) {
			us_bg = round_to_WORD(bg);
			bg = bg / get_normalized_value(&gfit);
		} else if (gfit.type == DATA_FLOAT) {
			us_bg = float_to_ushort_range(bg);
		} else return CMD_INVALID_IMAGE;
		siril_log_message(_("Background value (channel: #%d): %d (%.3e)\n"), layer, us_bg, bg);
	}
	return CMD_OK;
}

int process_bgnoise(int nb) {
	evaluate_noise_in_image();
	return CMD_OK;
}

int process_histo(int nb) {
	GError *error = NULL;
	gchar *end;
	int nlayer = g_ascii_strtoull(word[1], &end, 10);
	const gchar* clayer;

	if (end == word[1] || nlayer > 3 || nlayer < 0)
		return CMD_INVALID_IMAGE;
	gsl_histogram *histo = computeHisto(&gfit, nlayer);
	if (!isrgb(&gfit))
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
		redraw(REDRAW_OVERLAY);
	} else {
		set_cursor_waiting(TRUE);
		draw_sensor_tilt(&gfit);
		set_cursor_waiting(FALSE);
	}

	return CMD_OK;
}

int process_inspector(int nb) {
	compute_aberration_inspector();
	return CMD_OK;
}

int process_thresh(int nb){
	gchar *end1, *end2;
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int lo = g_ascii_strtoull(word[1], &end1, 10);
	if (end1 == word[1] || lo < 0 || lo > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	int hi = g_ascii_strtoull(word[2], &end2, 10);
	if (end2 == word[2] || hi < 0 || hi > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	if (lo >= hi) {
		siril_log_message(_("lo must be strictly smaller than hi\n"));
		return CMD_ARG_ERROR;
	}
	threshlo(&gfit, lo);
	threshhi(&gfit, hi);

	char log[90];
	sprintf(log, "Image clamped to [%d, %d]", lo, hi);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
	return CMD_OK;
}

int process_threshlo(int nb) {
	gchar *end;
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int lo = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || lo < 0 || lo > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	threshlo(&gfit, lo);

	char log[90];
	sprintf(log, "Image clamped to [%d, max]", lo);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
	return CMD_OK;
}

int process_threshhi(int nb) {
	gchar *end;
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int hi = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || hi < 0 || hi > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	threshhi(&gfit, hi);

	char log[90];
	sprintf(log, "Image clamped to [min, %d]", hi);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
	return CMD_OK;
}

int process_neg(int nb) {
	set_cursor_waiting(TRUE);
	pos_to_neg(&gfit);
	gfit.history = g_slist_append(gfit.history, strdup("Image made negative"));
	notify_gfit_modified();
	return CMD_OK;
}

int process_nozero(int nb){
	gchar *end;
	int level = g_ascii_strtoull(word[1], &end, 10);
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	if (end == word[1] || level < 0 || level > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	nozero(&gfit, (WORD)level);

	char log[90];
	sprintf(log, "Replaced zeros with %d", level);
	gfit.history = g_slist_append(gfit.history, strdup(log));
	notify_gfit_modified();
	return CMD_OK;
}

int process_ddp(int nb) {
	gchar *end;
	int level = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || level < 0) {
		siril_log_message(_("Level value is incorrect\n"));
		return CMD_ARG_ERROR;
	}
	float coeff = g_ascii_strtod(word[2], &end);
	if (end == word[2] || coeff < 0) {
		siril_log_message(_("Coeff value is incorrect\n"));
		return CMD_ARG_ERROR;
	}
	float sigma = g_ascii_strtod(word[3], &end);
	if (end == word[3] || sigma < 0) {
		siril_log_message(_("Sigma value is incorrect\n"));
		return CMD_ARG_ERROR;
	}
	ddp(&gfit, level, coeff, sigma);
	notify_gfit_modified();
	return CMD_OK;
}

int process_new(int nb){
	int width, height, layers;
	gchar *end, *endw, *endh;

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

	close_single_image();
	close_sequence(FALSE);

	fits *fit = &gfit;
	if (new_fit_image(&fit, width, height, layers, DATA_FLOAT))
		return CMD_GENERIC_ERROR;
	memset(gfit.fdata, 0, width * height * layers * sizeof(float));

	com.seq.current = UNRELATED_IMAGE;
	create_uniq_from_gfit(strdup(_("new empty image")), FALSE);
	open_single_image_from_gfit();
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
	visu(&gfit, low, high);
	return CMD_OK;
}

int process_ffill(int nb) {
	gchar *end;
	rectangle area;
	int level = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || level < 0 || level > USHRT_MAX) {
		siril_log_message(_("Value must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			gchar *endx, *endy, *endw, *endh;
			area.x = g_ascii_strtoull(word[2], &endx, 10);
			area.y = g_ascii_strtoull(word[3], &endy, 10);
			area.w = g_ascii_strtoull(word[4], &endw, 10);
			area.h = g_ascii_strtoull(word[5], &endh, 10);
			if (endx == word[2] || endy == word[3] || endw == word[4]
					|| endh == word[5] || (area.w + area.x > gfit.rx)
					|| (area.h + area.y > gfit.ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_message(_("Fill2: select a region or provide x, y, width, height\n"));
			return CMD_ARG_ERROR;
		}
	} else {
		area = com.selection;
	}
	int retval = fill(&gfit, level, &area);
	if (retval) {
		siril_log_message(_("Wrong parameters.\n"));
		return CMD_ARG_ERROR;
	}
	area.x = gfit.rx - area.x - area.w;
	area.y = gfit.ry - area.y - area.h;
	fill(&gfit, level, &area);
	redraw(REMAP_ALL);
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
			if (g_mkdir_with_parents(dirname, 0755) < 0) {
				siril_log_color_message(_("Cannot create output folder: %s\n"), "red", dirname);
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
	if (!com.script) {
		layer = select_vport(gui.cvport);
	} else {
		layer = (gfit.naxes[2] > 1) ? GLAYER : RLAYER;
	}

	struct starfinder_data *args = calloc(1, sizeof(struct starfinder_data));
	args->layer = layer;
	args->im.fit = &gfit;
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
	args->update_GUI = TRUE;
	siril_debug_print("findstar profiling %s stars\n", (com.pref.starfinder_conf.profile == PSF_GAUSSIAN) ? "Gaussian" : "Moffat");

	cmd_errors argparsing = parse_findstar(args, 1, nb);

	if (argparsing) {
		if (args->starfile) g_free(args->starfile);
		free(args);
		return argparsing;
	}
	start_in_new_thread(findstar_worker, args);

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
	if (!com.script && check_seq_is_comseq(seq)) { // we use vport only if seq is com.seq
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
#ifdef HAVE_WCSLIB
	args->save_eqcoords = TRUE;	// managed in apply_findstar_to_sequence()
#else
	args->save_eqcoords = FALSE;
#endif
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

	return apply_findstar_to_sequence(args);
}

int process_findhot(int nb){
	if (gfit.naxes[2] != 1) {
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

	deviant_pixel *dev = find_deviant_pixels(&gfit, sig, &icold, &ihot, FALSE);
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
		int y = gfit.ry - (int) dev[i].p.y - 1;  /* FITS is stored bottom to top */
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

int process_fix_xtrans(int nb) {
	fix_xtrans_ac(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	return CMD_OK;
}

int process_cosme(int nb) {
	gchar *filename;
	int retval = 0;

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

	int is_cfa = (word[0][5] == '_') ? 1 : 0;

	retval = apply_cosme_to_image(&gfit, file, is_cfa);

	g_free(filename);
	if (file)
		g_object_unref(file);
	if (retval)
		siril_log_color_message(_("There were some errors, please check your input file.\n"), "salmon");

	invalidate_stats_from_fit(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
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

	struct cosme_data *args = malloc(sizeof(struct cosme_data));

	if (g_str_has_prefix(word[3], "-prefix=")) {
		char *current = word[3], *value;
		value = current + 8;
		if (value[0] == '\0') {
			free_sequence(seq, TRUE);
			g_object_unref(file);
			if (!check_seq_is_comseq(seq))
				free_sequence(seq, TRUE);
			free(args);
			siril_log_message(_("Missing argument to %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
		args->prefix = strdup(value);
	} else {
		args->prefix = "cosme_";
	}

	args->seq = seq;
	args->is_cfa = (word[0][8] == '_') ? 1 : 0;
	args->file = file;
	args->fit = &gfit;

	apply_cosme_to_sequence(args);

	return CMD_OK;
}

int process_fmedian(int nb){
	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));
	gchar *end1, *end2;
	args->ksize = g_ascii_strtoull(word[1], &end1, 10);
	args->amount = g_ascii_strtod(word[2], &end2);
	args->iterations = 1;

	if (end1 == word[1] || !(args->ksize & 1) || args->ksize < 2 || args->ksize > 15) {
		siril_log_message(_("The size of the kernel MUST be odd and in the range [3, 15].\n"));
		free(args);
		return CMD_ARG_ERROR;
	}
	if (end2 == word[2] || args->amount < 0.0 || args->amount > 1.0) {
		siril_log_message(_("Modulation value MUST be between 0 and 1\n"));
		free(args);
		return CMD_ARG_ERROR;
	}
	args->fit = &gfit;

	start_in_new_thread(median_filter, args);

	return CMD_OK;
}

/* The name of this command should be COG in english but this choice
 * was done to be consistent with IRIS
 */
int process_cdg(int nb) {
	float x_avg, y_avg;

	if (!FindCentre(&gfit, &x_avg, &y_avg)) {
		siril_log_message(_("Center of gravity coordinates are (%.3lf, %.3lf)\n"), x_avg, y_avg);
		return CMD_OK;
	}
	return CMD_GENERIC_ERROR;
}

int process_clear(int nb) {
	if (com.script) return CMD_OK;
	GtkTextView *text = GTK_TEXT_VIEW(lookup_widget("output"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
	return CMD_OK;
}

int process_clearstar(int nb){
	clear_stars_list(TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(REDRAW_OVERLAY);
	redraw_previews();
	return CMD_OK;
}

int process_close(int nb) {
	close_sequence(FALSE);
	close_single_image();
	return CMD_OK;
}

int process_fill(int nb){
	rectangle area;
	gchar *end;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			gchar *endx, *endy, *endw, *endh;
			area.x = g_ascii_strtoull(word[2], &endx, 10);
			area.y = g_ascii_strtoull(word[3], &endy, 10);
			area.w = g_ascii_strtoull(word[4], &endw, 10);
			area.h = g_ascii_strtoull(word[5], &endh, 10);
			if (endx == word[2] || endy == word[3] || endw == word[4]
					|| endh == word[5] || (area.w + area.x > gfit.rx)
					|| (area.h + area.y > gfit.ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return CMD_ARG_ERROR;
			}
		}
		else {
			area.w = gfit.rx; area.h = gfit.ry;
			area.x = 0; area.y = 0;
		}
	} else {
		area = com.selection;
	}
	int level = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1] || level < 0 || level > USHRT_MAX) {
		siril_log_message(_("Value must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}
	int retval = fill(&gfit, level, &area);
	if (retval) {
		siril_log_message(_("Wrong parameters.\n"));
		return CMD_ARG_ERROR;
	}
	redraw(REMAP_ALL);
	return CMD_OK;
}

int process_offset(int nb) {
	gchar *end;

	int level = g_ascii_strtoull(word[1], &end, 10);
	if (end == word[1]) {
		siril_log_message(_("Wrong parameters.\n"));
		return CMD_ARG_ERROR;
	}
	off(&gfit, (float)level);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_scnr(int nb){
	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	args->type = SCNR_AVERAGE_NEUTRAL;
	args->amount = 0.0;
	args->fit = &gfit;
	args->preserve = TRUE;

	int argidx = 1;
	if (argidx < nb && !g_strcmp0(word[1], "-nopreserve")) {
		args->preserve = FALSE;
		argidx++;
	}

	if (argidx < nb) {
		gchar *end;
		args->type = g_ascii_strtoull(word[argidx], &end, 10);
		if (end == word[argidx] || args->type > 3) {
			siril_log_message(_("Type can either be 0 (average neutral) or 1 (maximum neutral), 2 (maximum mask) or 3 (additive mask)\n"));
			free(args);
			return CMD_ARG_ERROR;
		}
		argidx++;
	}

	if (argidx < nb && (args->type == SCNR_MAXIMUM_MASK || args->type == SCNR_ADDITIVE_MASK)) {
		gchar *end;
		args->amount = g_ascii_strtod(word[argidx], &end);
		if (end == word[argidx] || args->amount < 0.0 || args->amount > 1.0) {
			siril_log_message(_("Amount can only be [0, 1]\n"));
			free(args);
			return CMD_ARG_ERROR;
		}
	}

	char log[90];
	char amountstr[30] = "";
	if (args->type == SCNR_MAXIMUM_MASK || args->type == SCNR_ADDITIVE_MASK)
		sprintf(amountstr, "amount %.2f, ", args->amount);
	sprintf(log, "SCNR green removal (%s, %s%spreserving lightness)",
			scnr_type_to_string(args->type), amountstr,
			args->preserve ? "" : "not ");
	gfit.history = g_slist_append(gfit.history, strdup(log));

	start_in_new_thread(scnr, args);
	return CMD_OK;
}

int process_fft(int nb){
	struct fft_data *args = malloc(sizeof(struct fft_data));

	args->fit = &gfit;
	args->type = strdup(word[0]);
	args->modulus = strdup(word[1]);
	args->phase = strdup(word[2]);
	args->type_order = 0;

	start_in_new_thread(fourier_transform, args);

	return CMD_OK;
}

int process_fixbanding(int nb) {
	struct banding_data *args = malloc(sizeof(struct banding_data));
	args->seq = NULL;
	gchar *end1, *end2;

	args->amount = g_ascii_strtod(word[1], &end1);
	if (end1 == word[1] || args->amount < 0 || args->amount > 4) {
		siril_log_message(_("Amount value must be in the [0, 4] range.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}
	args->sigma = g_ascii_strtod(word[2], &end2);
	if (end2 == word[2] || args->sigma < 0 || args->sigma > 5) {
		siril_log_message(_("1/sigma value must be in the [0, 5] range.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}

	args->protect_highlights = TRUE;
	args->applyRotation = FALSE;
	args->fit = &gfit;
	if (nb > 3) {
		int arg_index = 3;
		while (arg_index < nb && word[arg_index]) {
			char *arg = word[arg_index];
			if (!g_strcmp0(arg, "-vertical")) {
				args->applyRotation = TRUE;
			} else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
			arg_index++;
		}
	}
	start_in_new_thread(BandingEngineThreaded, args);
	return CMD_OK;
}

int process_seq_fixbanding(int nb) {
	struct banding_data *args = malloc(sizeof(struct banding_data));
	gchar *end1 = NULL, *end2 = NULL;
	args->seq = NULL;

	if (word[1] && word[1][0] != '\0') {
		if (!(args->seq = load_sequence(word[1], NULL))) {
			if (!args->seq) {
				free(args);
				return CMD_SEQUENCE_NOT_FOUND;
			}
		}
	}
	if (!args->seq) {
		free(args);
		return CMD_SEQUENCE_NOT_FOUND;
	}
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
	args->seqEntry = g_strdup("unband_");
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
					g_free(args);
					return CMD_ARG_ERROR;
				}
				args->seqEntry = g_strdup(value);
			} else if (!g_strcmp0(arg, "-vertical")) {
				args->applyRotation = TRUE;
			} else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
				if (!check_seq_is_comseq(args->seq))
					free_sequence(args->seq, TRUE);
				g_free(args);
				return CMD_ARG_ERROR;
			}
			arg_index++;
		}
	}
	apply_banding_to_sequence(args);
	return CMD_OK;
}

int process_subsky(int nb) {
	sequence *seq = NULL;
	int degree = 0, samples = 20;
	double tolerance = 1.0, smooth = 0.5;
	background_interpolation interp;
	char *prefix = NULL;

	int arg_index = 1;
	gboolean is_sequence = (word[0][2] == 'q');

	if (is_sequence) {
		arg_index = 2;
		seq = load_sequence(word[1], NULL);
		if (!seq) {
			return CMD_SEQUENCE_NOT_FOUND;
		}
	} else {
		if (!single_image_is_loaded()) return CMD_IMAGE_NOT_FOUND;
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
		if (g_str_has_prefix(arg, "-prefix=")) {
			char *value = arg + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), arg);
				g_free(prefix);
				return CMD_ARG_ERROR;
			}
			g_free(prefix);
			prefix = strdup(value);
		}
		else if (g_str_has_prefix(arg, "-samples=")) {
			gchar *end;
			char *value = arg + 9;
			samples = g_ascii_strtoull(value, &end, 10);
			if (end == value || samples <= 1) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				g_free(prefix);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-tolerance=")) {
			char *next, *value = arg + 11;
			tolerance = g_ascii_strtod(value, &next);
			if (next == value || tolerance < 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				g_free(prefix);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-smooth=")) {
			char *next, *value = arg + 8;
			smooth = g_ascii_strtod(value, &next);
			if (next == value || smooth < 0.0 || smooth > 1.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				g_free(prefix);
				return CMD_ARG_ERROR;
			}
			if (interp != BACKGROUND_INTER_RBF)
				siril_log_color_message(_("smooth parameter is unused with the polynomial model, ignoring.\n"), "salmon");
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			g_free(prefix);
			return CMD_ARG_ERROR;
		}
		arg_index++;
	}

	struct background_data *args = malloc(sizeof(struct background_data));
	args->nb_of_samples = samples;
	args->tolerance = tolerance;
	args->correction = BACKGROUND_CORRECTION_SUBTRACT;
	args->interpolation_method = interp;
	args->degree = (poly_order) (degree - 1);
	args->smoothing = smooth;
	args->threads = com.max_thread;
	args->from_ui = FALSE;

	if (is_sequence) {
		args->seq = seq;
		args->dither = TRUE;
		args->seqEntry = prefix ? prefix : "bkg_";

		apply_background_extraction_to_sequence(args);
	} else {
		args->seq = NULL;
		args->dither = FALSE;
		args->seqEntry = NULL;
		args->fit = &gfit;
		g_free(prefix);

		if (!generate_background_samples(samples, tolerance))
			start_in_new_thread(remove_gradient_from_image, args);
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

	struct cosmetic_data *args = malloc(sizeof(struct cosmetic_data));
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
	args->fit = &gfit;

	if (is_sequence) {
		args->seqEntry = "cc_";
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
		start_in_new_thread(autoDetectThreaded, args);
	}

	return CMD_OK;
}

int select_unselect(gboolean select) {
	char *end1, *end2;
	int from = g_ascii_strtoull(word[1], &end1, 10);
	int to = g_ascii_strtoull(word[2], &end2, 10);
	if (end1 == word[1] || from < 1 || from >= com.seq.number) {
		siril_log_message(_("The first argument must be between 1 and the number of images.\n"));
		return CMD_ARG_ERROR;
	}

	if (end2 == word[2] || to <= from || to >= com.seq.number) {
		siril_log_message(_("The second argument must be between from and the number of images.\n"));
		return CMD_ARG_ERROR;
	}
	gboolean current_updated = FALSE;
	for (int i = from - 1; i <= to - 1; i++) { // use real index
		if (i >= com.seq.number) break;
		if (com.seq.imgparam[i].incl != select) {
			com.seq.imgparam[i].incl = select;
			if (!com.headless)
				sequence_list_change_selection_index(i, i);
			if (select)
				com.seq.selnum++;
			else	com.seq.selnum--;
			if (i + 1 == com.seq.current)
				current_updated = TRUE;
		}
		if (!select && com.seq.reference_image == i) {
			com.seq.reference_image = -1;
			if (!com.headless) {
				sequence_list_change_reference();
				adjust_refimage(com.seq.current);
			}
		}
	}

	if (!com.headless) {
		if (current_updated) {
			redraw(REDRAW_OVERLAY);
			adjust_sellabel();
		}
		drawPlot();
		update_reg_interface(FALSE);
		adjust_sellabel();
	}
	writeseqfile(&com.seq);
	siril_log_message(_("Selection update finished, %d images are selected in the sequence\n"), com.seq.selnum);

	return CMD_OK;
}

int process_select(int nb){
	return select_unselect(TRUE);
}

int process_unselect(int nb){
	return select_unselect(FALSE);
}

int process_split(int nb){
	struct extract_channels_data *args = malloc(sizeof(struct extract_channels_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return CMD_ALLOC_ERROR;
	}

	args->channel[0] = g_strdup_printf("%s%s", word[1], com.pref.ext);
	args->channel[1] = g_strdup_printf("%s%s", word[2], com.pref.ext);
	args->channel[2] = g_strdup_printf("%s%s", word[3], com.pref.ext);

	args->fit = calloc(1, sizeof(fits));
	if (copyfits(&gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Could not copy the input image, aborting.\n"));
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

	copy_fits_metadata(&gfit, args->fit);
	start_in_new_thread(extract_channels, args);
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

	if (gfit.type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(&gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, 0) ||
				save1fits16(cfa1, &f_cfa1, 0) ||
				save1fits16(cfa2, &f_cfa2, 0) ||
				save1fits16(cfa3, &f_cfa3, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(&gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
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

	sensor_pattern pattern = get_bayer_pattern(&gfit);

	gchar *green = g_strdup_printf("Green_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractGreen_ushort(&gfit, &f_green, pattern))) {
			ret = save1fits16(green, &f_green, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractGreen_float(&gfit, &f_green, pattern))) {
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

int process_extractHa(int nb) {
	char *filename = NULL;
	int ret = 1;

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

	sensor_pattern pattern = get_bayer_pattern(&gfit);

	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractHa_ushort(&gfit, &f_Ha, pattern))) {
			ret = save1fits16(Ha, &f_Ha, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractHa_float(&gfit, &f_Ha, pattern))) {
			ret = save1fits32(Ha, &f_Ha, 0);
		}
	} else ret = CMD_INVALID_IMAGE;

	g_free(Ha);
	clearfits(&f_Ha);
	free(filename);
	return ret;
}

int process_extractHaOIII(int nb) {
	char *filename = NULL;
	int ret = 1;
	extraction_scaling scaling = SCALING_NONE;

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
	if (word[1]) {
		if (g_str_has_prefix(word[1], "-resample=")) {
			char *current = word[1], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[1]);
				g_free(filename);
				return CMD_ARG_ERROR;
			} else if (!strcasecmp(value, "ha")) {
				scaling = SCALING_HA_UP;
			} else if (!strcasecmp(value, "oiii")) {
				scaling = SCALING_OIII_DOWN;
			}
		}
	}

	sensor_pattern pattern = get_bayer_pattern(&gfit);

	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	gchar *OIII = g_strdup_printf("OIII_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractHaOIII_ushort(&gfit, &f_Ha, &f_OIII, pattern, scaling))) {
			ret = save1fits16(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractHaOIII_float(&gfit, &f_Ha, &f_OIII, pattern, scaling))) {
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

int process_seq_mtf(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	struct mtf_data *args = malloc(sizeof(struct mtf_data));

	args->seq = seq;
	args->fit = &gfit;
	args->seqEntry = "mtf_";
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

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->fit = &gfit;
	args->seqEntry = "CFA_"; // propose to default to "CFA" for consistency of output names with single image split_cfa

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
						g_free(args);
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				g_free(args);
				return CMD_ARG_ERROR;
			}
		}
	}

	apply_split_cfa_to_sequence(args);

	return CMD_OK;
}

int process_seq_merge_cfa(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_FOR_CFA_IMAGE;
	}

	struct merge_cfa_data *args = calloc(1, sizeof(struct merge_cfa_data));

	if (!strcmp(word[2], "RGGB")) {
		args->pattern = BAYER_FILTER_RGGB;
	} else if (!strcmp(word[2], "BGGR")) {
		args->pattern = BAYER_FILTER_BGGR;
	} else if (!strcmp(word[2], "GBRG")) {
		args->pattern = BAYER_FILTER_GBRG;
	} else if (!strcmp(word[2], "GRBG")) {
		args->pattern = BAYER_FILTER_GRBG;
	} else {
		siril_log_color_message(_("Invalid Bayer matrix specified!\n"), "red");
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		g_free(args);
		return CMD_ARG_ERROR;
	}
	siril_log_message(_("Reconstructing %s Bayer matrix.\n"), word[2]);
	args->seq = seq;
	args->seqEntryIn = "CFA_"; // propose to default to "CFA" for consistency of output names with single image split_cfa
	args->seqEntryOut = "mCFA_"; // propose to default to "CFA" for consistency of output names with single image split_cfa

	int startoptargs = 3;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-inprefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						g_free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					args->seqEntryIn = strdup(value);
				} else if (g_str_has_prefix(word[i], "-outprefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						g_free(args);
						return CMD_ARG_ERROR;
					}
					args->seqEntryOut = strdup(value);
				}
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				g_free(args);
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

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->seqEntry = "Ha_";

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						g_free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
					g_free(args);
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

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->seqEntry = "Green_";

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"), word[i]);
						g_free(args);
						if (!check_seq_is_comseq(seq))
							free_sequence(seq, TRUE);
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
					if (!check_seq_is_comseq(seq))
						free_sequence(seq, TRUE);
					g_free(args);
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

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));
	args->scaling = SCALING_NONE;

	if (word[2]) {
		if (g_str_has_prefix(word[2], "-resample=")) {
			char *current = word[2], *value;
			value = current + 10;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), word[2]);
				if (!check_seq_is_comseq(seq))
					free_sequence(seq, TRUE);
				free(args);
				return CMD_ARG_ERROR;
			} else if (!strcmp(value, "ha")) {
				args->scaling = SCALING_HA_UP;
			} else if (!strcmp(value, "oiii")) {
				args->scaling = SCALING_OIII_DOWN;
			}
		}
	}

	args->seq = seq;
	args->seqEntry = ""; // not used

	apply_extractHaOIII_to_sequence(args);

	return CMD_OK;
}

int process_stat(int nb){
	int nplane, layer, option = STATS_BASIC;
	char layername[6];

	nplane = gfit.naxes[2];
	gboolean cfa = FALSE;
	int argidx = 1;
	if (nb == 2 && !g_strcmp0(word[1], "-cfa") && nplane == 1 && gfit.bayer_pattern[0] != '\0') {
		siril_debug_print("Running stats on CFA\n");
		nplane = 3;
		cfa = TRUE;
		argidx++;
	}
	if (nb == argidx + 1 && !g_strcmp0(word[argidx], "main"))
		option = STATS_MAIN;

	for (layer = 0; layer < nplane; layer++) {
		int super_layer = layer;
		if (cfa)
			super_layer = -layer - 1;
		imstats* stat = statistics(NULL, -1, &gfit, super_layer, &com.selection, STATS_MAIN, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Statistics computation failed for channel %d (all nil?).\n"), layer);
			continue;
		}

		switch (layer) {
			case 0:
				if (nplane == 1)
					strcpy(layername, "B&W");
				else
					strcpy(layername, "Red");
				break;
			case 1:
				strcpy(layername, "Green");
				break;
			case 2:
				strcpy(layername, "Blue");
				break;
		}

		if (option == STATS_BASIC) {
			if (gfit.type == DATA_USHORT) {
				siril_log_message(_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
							"Min: %0.1f, Max: %0.1f, bgnoise: %0.1f\n"),
						layername, stat->mean, stat->median, stat->sigma,
						stat->min, stat->max, stat->bgnoise);
			} else {
				siril_log_message(_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
							"Min: %0.1f, Max: %0.1f, bgnoise: %0.1f\n"),
						layername, stat->mean * USHRT_MAX_DOUBLE,
						stat->median * USHRT_MAX_DOUBLE,
						stat->sigma * USHRT_MAX_DOUBLE,
						stat->min * USHRT_MAX_DOUBLE,
						stat->max * USHRT_MAX_DOUBLE,
						stat->bgnoise * USHRT_MAX_DOUBLE);
			}

		} else if (option == STATS_MAIN) {
			if (gfit.type == DATA_USHORT) {
				siril_log_message(_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
							"Min: %0.1f, Max: %0.1f, bgnoise: %0.1f, "
							"avgDev: %0.1f, MAD: %0.1f, sqrt(BWMV): %0.1f\n"),
						layername, stat->mean, stat->median, stat->sigma,
						stat->min, stat->max, stat->bgnoise, stat->avgDev,
						stat->mad, stat->sqrtbwmv);
			} else {
				siril_log_message(_("%s layer: Mean: %0.1f, Median: %0.1f, Sigma: %0.1f, "
							"Min: %0.1f, Max: %0.1f, bgnoise: %0.1f, "
							"avgDev: %0.1f, MAD: %0.1f, sqrt(BWMV): %0.1f\n"),
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
	args->seqEntry = ""; // not used
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
#ifdef HAVE_JSON_GLIB
	/* we need both the file descriptor to read the header and the data to
	 * compute stats, so we need the file name even in the case of a loaded
	 * image, but we may use the loaded data if the option is provided, to
	 * avoid reading it for nothing */
	char *input_filename = word[1];
	gchar *output_filename = NULL;
	gboolean use_gfit = FALSE, compute_stats = TRUE;
	for (int i = 2; i < nb; i++) {
		if (g_str_has_prefix(word[i], "-out=") && word[i][5] != '\0')
			output_filename = g_strdup(word[i] + 5);
		else if (!strcmp(word[i], "-stats_from_loaded")) {
			use_gfit = TRUE;
			if (!gfit.rx || !gfit.ry) {
				siril_log_color_message(_("No image appears to be loaded, reloading from '%s'\n"), "salmon", input_filename);
				use_gfit = FALSE;
			}
		} else if (!strcmp(word[i], "-nostats"))
			compute_stats = FALSE;
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}
	if (!output_filename)
		output_filename = replace_ext(input_filename, ".json");

	int status = 0;
	fitsfile *fptr;
	if (siril_fits_open_diskfile_img(&fptr, input_filename, READONLY, &status)) {
		report_fits_error(status);
		return CMD_GENERIC_ERROR;
	}

	GSList *header = read_header_keyvals_strings(fptr);

	imstats *stats[3] = { NULL };
	int nb_channels;
	if (compute_stats) {
		if (use_gfit) {
			compute_all_channels_statistics_single_image(&gfit, STATS_BASIC | STATS_FOR_CFA, MULTI_THREADED, stats);
			gboolean cfa = gfit.bayer_pattern[0] != '\0';
			nb_channels = cfa ? 3 : (int)gfit.naxes[2];
		} else {
			fits fit = { 0 };
			fit.fptr = fptr;
			if (read_fits_metadata(&fit) || read_fits_with_convert(&fit, input_filename, FALSE)) {
				fits_close_file(fptr, &status);
				return CMD_GENERIC_ERROR;
			}
			compute_all_channels_statistics_single_image(&fit, STATS_BASIC | STATS_FOR_CFA, MULTI_THREADED, stats);
			gboolean cfa = fit.bayer_pattern[0] != '\0';
			nb_channels = cfa ? 3 : (int)fit.naxes[2];
			clearfits(&fit);
		}
	}
	fits_close_file(fptr, &status);

	// https://gnome.pages.gitlab.gnome.org/json-glib/class.Builder.html
	JsonBuilder *builder = json_builder_new();
	json_builder_begin_object(builder);
	json_builder_set_member_name(builder, "headers");
	json_builder_begin_array(builder);
	GSList *ptr = header;
	while (ptr) {
		json_builder_begin_object(builder);
		header_record *r = ptr->data;
		json_builder_set_member_name(builder, "key");
		json_builder_add_string_value(builder, r->key);
		json_builder_set_member_name(builder, "value");
		json_builder_add_string_value(builder, r->value);
		json_builder_end_object(builder);
		ptr = ptr->next;
	}
	json_builder_end_array(builder);

	if (compute_stats) {
		json_builder_set_member_name(builder, "statistics");
		json_builder_begin_object(builder);
		for (int i = 0; i < nb_channels; ++i) {
			if (stats[i]) {
				char channame[20];
				sprintf(channame, "channel%d", i);
				json_builder_set_member_name(builder, channame);
				json_builder_begin_object(builder);
				json_builder_set_member_name(builder, "mean");
				json_builder_add_double_value(builder, stats[i]->mean);
				json_builder_set_member_name(builder, "median");
				json_builder_add_double_value(builder, stats[i]->median);
				json_builder_set_member_name(builder, "sigma");
				json_builder_add_double_value(builder, stats[i]->sigma);
				json_builder_set_member_name(builder, "noise");
				json_builder_add_double_value(builder, stats[i]->bgnoise);
				json_builder_set_member_name(builder, "min");
				json_builder_add_double_value(builder, stats[i]->min);
				json_builder_set_member_name(builder, "max");
				json_builder_add_double_value(builder, stats[i]->max);
				json_builder_set_member_name(builder, "total_pix_count");
				json_builder_add_double_value(builder, stats[i]->total);
				json_builder_set_member_name(builder, "good_pix_count");
				json_builder_add_double_value(builder, stats[i]->ngoodpix);
				json_builder_end_object(builder);
				// BASIC: median, mean, sigma, noise, min, max
				free_stats(stats[i]);
			}
		}
		json_builder_end_object(builder);
	}
	json_builder_end_object(builder);

	JsonGenerator *gen = json_generator_new();
	JsonNode *root = json_builder_get_root(builder);
	json_generator_set_root(gen, root);
	json_generator_set_pretty(gen, TRUE);
	GError *err = NULL;
	int retval = CMD_OK;
	if (!json_generator_to_file(gen, output_filename, &err)) {
		siril_log_message(_("Failed to save the JSON file %s: %s\n"), output_filename, err->message);
		retval = CMD_GENERIC_ERROR;
	}
	else siril_log_message(_("Save metadata to the JSON file '%s'\n"), output_filename);
	g_free(output_filename);

#ifdef DEBUG_TEST
	gchar *str = json_generator_to_data(gen, NULL);
	printf("JSON:\n%s\n", str);
#endif

	json_node_free(root);
	g_object_unref(gen);
	g_object_unref(builder);
	return retval;
#else
	siril_log_message(_("json-glib was not found at build time, cannot proceed. Install and rebuild.\n"));
	return CMD_GENERIC_ERROR;
#endif
}

int header_hook(struct generic_seq_metadata_args *args, fitsfile *fptr, int index) {
	char str[FLEN_VALUE] = { 0 };
	int status = 0;
	fits_read_keyword(fptr, args->key, str, NULL, &status);
	if (status)
		strcpy(str, "not found");
	siril_log_message(_("Image %d, %s = %s\n"), index, args->key, str);
	return status;
}

int process_seq_header(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;
	if (seq->type != SEQ_REGULAR && seq->type != SEQ_FITSEQ) {
		siril_log_message(_("This command can only run for FITS images\n"));
		if (!check_seq_is_comseq(seq))
			free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}

	struct generic_seq_metadata_args *args = malloc(sizeof(struct generic_seq_metadata_args));
	args->seq = seq;
	args->key = g_strdup(word[2]);
	args->image_hook = header_hook;
	start_in_new_thread(generic_sequence_metadata_worker, args);
	return 0;
}

int process_convertraw(int nb) {
	if (!com.wd) {
		siril_log_message(_("Conversion: no working directory set.\n"));
		return CMD_NO_CWD;
	}

	if (word[1][0] == '-') {
		siril_log_message(_("First argument is the converted sequence name and shall not start with a -\n"));
		return CMD_ARG_ERROR;
	}
	gchar *destroot = g_strdup(word[1]);
	int idx = 1;
	gboolean debayer = FALSE;
	sequence_type output = SEQ_REGULAR;

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (!strcmp(current, "-debayer")) {
			debayer = TRUE;
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
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return CMD_GENERIC_ERROR;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	}

	GDir *dir;
	GError *error = NULL;
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_clear_error(&error);
		return CMD_NO_CWD;
	}

	int count = 0;
	const gchar *file;
	GList *list = NULL;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext || file[0] == '.')
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPERAW) {
			if (output == SEQ_SER && !g_ascii_strcasecmp(ext, "raf") && !debayer) {
				siril_log_message(_("FujiFilm XTRANS sensors are not supported by SER v2 (CFA-style) standard. You may use FITS sequences instead."));
				g_list_free_full(list, g_free);
				return CMD_GENERIC_ERROR;
			}
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No RAW files were found for conversion\n"));
		g_list_free_full(list, g_free);
		return CMD_GENERIC_ERROR;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	char **files_to_convert = glist_to_array(list, &count);

	siril_log_color_message(_("Conversion: processing %d RAW files...\n"), "green", count);

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_convert;
	args->total = count;
	if (output == SEQ_REGULAR)
		args->destroot = format_basename(destroot, TRUE);
	else
		args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = debayer;
	args->output_type = output;
	args->multiple_output = FALSE;
	args->make_link = FALSE;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);
	return CMD_OK;
}

int process_link(int nb) {
	if (word[1][0] == '-') {
		siril_log_message(_("First argument is the converted sequence name and shall not start with a -\n"));
		return CMD_ARG_ERROR;
	}
	gchar *destroot = g_strdup(word[1]);
	int idx = 1;

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			gchar *end;
			idx = g_ascii_strtoull(value, &end, 10);
			if (end == value || idx <= 0 || idx >= INDEX_MAX) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return CMD_GENERIC_ERROR;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
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
		return CMD_GENERIC_ERROR;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	char **files_to_link = glist_to_array(list, &count);

	gchar *str = ngettext("Link: processing %d FITS file...\n", "Link: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);

	if (!com.wd) {
		siril_log_message(_("Link: no working directory set.\n"));
		g_strfreev(files_to_link);
		return CMD_GENERIC_ERROR;
	}

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_link;
	args->total = count;
	args->destroot = format_basename(destroot, TRUE);
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = FALSE;
	args->multiple_output = FALSE;
	args->output_type = SEQ_REGULAR; // fallback if symlink does not work
	args->make_link = TRUE;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);

	return CMD_OK;
}

int process_convert(int nb) {
	if (word[1][0] == '-') {
		siril_log_message(_("First argument is the converted sequence name and shall not start with a -\n"));
		return CMD_ARG_ERROR;
	}
	gchar *destroot = g_strdup(word[1]);
	int idx = 1;
	gboolean debayer = FALSE;
	gboolean make_link = TRUE;
	sequence_type output = SEQ_REGULAR;

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
				return CMD_ARG_ERROR;
			}
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return CMD_ARG_ERROR;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return CMD_GENERIC_ERROR;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
	}

	GDir *dir;
	GError *error = NULL;
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Convert: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Convert: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
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
		if (type != TYPEUNDEF && type != TYPEAVI && type != TYPESER) {
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No files were found for convert\n"));
		return CMD_GENERIC_ERROR;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	gchar **files_to_link = glist_to_array(list, &count);

	gchar *str = ngettext("Convert: processing %d FITS file...\n", "Convert: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);

	if (!com.wd) {
		siril_log_message(_("Convert: no working directory set.\n"));
		g_strfreev(files_to_link);
		return CMD_NO_CWD;
	}

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_link;
	args->total = count;
	if (output == SEQ_REGULAR)
		args->destroot = format_basename(destroot, TRUE);
	else
		args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = debayer;
	args->multiple_output = FALSE;
	args->output_type = output;
	args->make_link = make_link;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);

	return CMD_OK;
}

int process_register(int nb) {
	struct registration_args *reg_args = NULL;
	struct registration_method *method = NULL;
	char *msg;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}
	if (check_seq_is_comseq(seq)) {
		free_sequence(seq, TRUE);
		seq = &com.seq;
	}

	reg_args = calloc(1, sizeof(struct registration_args));

	/* filling the arguments for registration */
	reg_args->seq = seq;
	reg_args->reference_image = sequence_find_refimage(seq);
	reg_args->follow_star = FALSE;
	reg_args->matchSelection = FALSE;
	reg_args->no_output = FALSE;
	reg_args->x2upscale = FALSE;
	reg_args->prefix = "r_";
	reg_args->min_pairs = 10; // 10 is good enough to ensure good matching
	reg_args->max_stars_candidates = MAX_STARS_FITTED;
	reg_args->type = HOMOGRAPHY_TRANSFORMATION;
	reg_args->layer = (reg_args->seq->nb_layers == 3) ? 1 : 0;
	reg_args->interpolation = OPENCV_LANCZOS4;
	reg_args->clamp = TRUE;

	/* check for options */
	for (int i = 2; i < nb; i++) {
		if (!strcmp(word[i], "-drizzle")) {
			reg_args->x2upscale = TRUE;
		} else if (!strcmp(word[i], "-noout")) {
			reg_args->no_output = TRUE;
		} else if (!strcmp(word[i], "-noclamp")) {
			reg_args->clamp = FALSE;
		} else if (!strcmp(word[i], "-2pass")) {
			reg_args->two_pass = TRUE;
			reg_args->no_output = TRUE;
		} else if (!strcmp(word[i], "-nostarlist")) {
			reg_args->no_starlist = TRUE;
		} else if (!strcmp(word[i], "-selected")) {
			reg_args->filters.filter_included = TRUE;
		} else if (g_str_has_prefix(word[i], "-transf=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "shift", 5)) {
#ifdef HAVE_CV44
				reg_args->type = SHIFT_TRANSFORMATION;
				continue;
#else
				siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
				goto terminate_register_on_error;
#endif
			}
			if(!g_ascii_strncasecmp(value, "similarity", 10)) {
				reg_args->type = SIMILARITY_TRANSFORMATION;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "affine", 6)) {
				reg_args->type = AFFINE_TRANSFORMATION;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "homography", 10)) {
				reg_args->type = HOMOGRAPHY_TRANSFORMATION;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-layer=")) {
			if (reg_args->seq->nb_layers == 1) {  // handling mono case
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
			reg_args->layer = layer;
		} else if (g_str_has_prefix(word[i], "-prefix=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			reg_args->prefix = strdup(value);
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
				gchar *str = g_strdup_printf(_("%d smaller than minimum allowable star pairs: %d, aborting.\n"), min_pairs, reg_args->min_pairs);
				siril_log_message(str);
				g_free(str);
				goto terminate_register_on_error;
			}
			reg_args->min_pairs = min_pairs;
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
			reg_args->max_stars_candidates = max_stars;
		} else if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				reg_args->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				reg_args->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				reg_args->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				reg_args->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "none", 4) || !g_ascii_strncasecmp(value, "no", 2)) {
				reg_args->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				reg_args->interpolation = OPENCV_AREA;
				continue;
			}
			siril_log_message(_("Unknown transformation type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			goto terminate_register_on_error;
		}
	}

	/* getting the selected registration method */
	method = malloc(sizeof(struct registration_method));
	if (reg_args->two_pass) {
		method->name = _("Two-Pass Global Star Alignment (deep-sky)");
		method->method_ptr = &register_multi_step_global;
	} else {
		method->name = _("Global Star Alignment (deep-sky)");
		method->method_ptr = &register_star_alignment;
	}
	method->sel = REQUIRES_NO_SELECTION;
	method->type = REGTYPE_DEEPSKY;
	reg_args->func = method->method_ptr;

	if (reg_args->no_starlist && !reg_args->two_pass)
		siril_log_message(_("The -nostarlist option has an effect only when -2pass is used, ignoring\n"));

	// testing free space
	if (!reg_args->no_output) {
		// first, remove the files that we are about to create
		remove_prefixed_sequence_files(reg_args->seq, reg_args->prefix);

		int nb_frames = reg_args->filters.filter_included ? reg_args->seq->selnum : reg_args->seq->number;
		int64_t size = seq_compute_size(reg_args->seq, nb_frames, get_data_type(seq->bitpix));
		if (reg_args->x2upscale)
			size *= 4;
		if (test_available_space(size)) {
			siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
			goto terminate_register_on_error;
		}
	}
	if (reg_args->interpolation == OPENCV_NONE && !(reg_args->type == SHIFT_TRANSFORMATION)) {
#ifdef HAVE_CV44
		reg_args->type = SHIFT_TRANSFORMATION;
		siril_log_color_message(_("Forcing the registration transformation to shift, which is the only transformation compatible with no interpolation\n"), "salmon");

#else
		siril_log_color_message(_("Forcing the registration transformation to shift, which is the only transformation compatible with no interpolation, is not compatible with OpenCV below 4.4. Aborting\n"), "red");
		goto terminate_register_on_error;
#endif
	}

	if (reg_args->interpolation == OPENCV_NONE && (reg_args->x2upscale || reg_args->seq->is_variable)) {
		siril_log_color_message(_("When interpolation is set to None, the images must be of same size and no upscaling can be applied. Aborting\n"), "red");
		goto terminate_register_on_error;
	}

	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;	// don't load it for command line execution

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"), "green", method->name);
	free(method);
	msg[strlen(msg) - 1] = '\0';

	if (reg_args->interpolation == OPENCV_AREA ||
			reg_args->interpolation == OPENCV_LINEAR ||
			reg_args->interpolation == OPENCV_NEAREST ||
			reg_args->interpolation == OPENCV_NONE ||
			reg_args->no_output)
		reg_args->clamp = FALSE;
	if (reg_args->clamp)
		siril_log_message(_("Interpolation clamping active\n"));

	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
	return CMD_OK;

terminate_register_on_error:
	g_free(reg_args);
	free_sequence(seq, TRUE);
	free(method);
	return CMD_ARG_ERROR;
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
	struct registration_args *reg_args = NULL;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;

	reg_args = calloc(1, sizeof(struct registration_args));

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
	reg_args->func = &register_apply_reg;
	reg_args->seq = seq;
	reg_args->reference_image = sequence_find_refimage(seq);
	reg_args->no_output = FALSE;
	reg_args->x2upscale = FALSE;
	reg_args->prefix = "r_";
	reg_args->layer = layer;
	reg_args->interpolation = OPENCV_LANCZOS4;
	reg_args->clamp = TRUE;
	reg_args->framing = FRAMING_CURRENT;

	/* check for options */
	for (int i = 2; i < nb; i++) {
		if (!strcmp(word[i], "-drizzle")) {
			reg_args->x2upscale = TRUE;
		} else if (g_str_has_prefix(word[i], "-prefix=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			reg_args->prefix = strdup(value);
		} else if (!strcmp(word[i], "-noclamp")) {
			reg_args->clamp = FALSE;
		} else if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_ascii_strncasecmp(value, "nearest", 7) || !g_ascii_strncasecmp(value, "ne", 2)) {
				reg_args->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cubic", 5) || !g_ascii_strncasecmp(value, "cu", 2)) {
				reg_args->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "lanczos4", 8) || !g_ascii_strncasecmp(value, "la", 2)) {
				reg_args->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "linear", 6) || !g_ascii_strncasecmp(value, "li", 2)) {
				reg_args->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "none", 4) || !g_ascii_strncasecmp(value, "no", 2)) {
				reg_args->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "area", 4) || !g_ascii_strncasecmp(value, "ar", 2)) {
				reg_args->interpolation = OPENCV_AREA;
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
				reg_args->framing = FRAMING_CURRENT;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "min", 3)) {
				reg_args->framing = FRAMING_MIN;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "max", 3)) {
				reg_args->framing = FRAMING_MAX;
				continue;
			}
			if(!g_ascii_strncasecmp(value, "cog", 3)) {
				reg_args->framing = FRAMING_COG;
				continue;
			}
			siril_log_message(_("Unknown framing type %s, aborting.\n"), value);
			goto terminate_register_on_error;
		} else if (g_str_has_prefix(word[i], "-layer=")) {
			if (reg_args->seq->nb_layers == 1) {  // handling mono case
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
			reg_args->layer = layer2;
		} else if (parse_filter_args(word[i], &reg_args->filters)) {
			;
		} else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
			goto terminate_register_on_error;
		}
	}

	// sanity checks are done in register_apply_reg

	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;	// don't load it for command line execution

	if (reg_args->interpolation == OPENCV_AREA || reg_args->interpolation == OPENCV_LINEAR || reg_args->interpolation == OPENCV_NEAREST || reg_args->interpolation == OPENCV_NONE)
		reg_args->clamp = FALSE;
	if (reg_args->clamp)
		siril_log_message(_("Interpolation clamping active\n"));

	set_progress_bar_data(_("Registration: Applying existing data"), PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
	return CMD_OK;

terminate_register_on_error:
	if (!check_seq_is_comseq(seq))
		free_sequence(seq, TRUE);
	g_free(reg_args);
	return CMD_ARG_ERROR;
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
		} else if (!strcmp(current, "-weight_from_noise")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Weighting is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("Weighting is allowed only if normalization has been activated, ignoring.\n"));
			} else if (arg->apply_nbstack_weights || arg->apply_wfwhm_weights || arg->apply_nbstars_weights) {
				siril_log_message(_("Only one weighting method can be used\n"));
			} else {
				arg->apply_noise_weights = TRUE;
			}
		} else if (!strcmp(current, "-weight_from_wfwhm")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Weighting is allowed only with average stacking, ignoring.\n"));
			} else if (arg->apply_nbstack_weights || arg->apply_noise_weights || arg->apply_nbstars_weights) {
				siril_log_message(_("Only one weighting method can be used\n"));
			} else {
				arg->apply_wfwhm_weights = TRUE;
			}
		} else if (!strcmp(current, "-weight_from_nbstars")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Weighting is allowed only with average stacking, ignoring.\n"));
			} else if (arg->apply_nbstack_weights || arg->apply_noise_weights || arg->apply_wfwhm_weights) {
				siril_log_message(_("Only one weighting method can be used\n"));
			} else {
				arg->apply_nbstars_weights = TRUE;
			}
		} else if (!strcmp(current, "-weight_from_nbstack")) {
			if (!rej_options_allowed) {
				siril_log_message(_("Weighting is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("Weighting is allowed only if normalization has been activated, ignoring.\n"));
			} else if (arg->apply_noise_weights || arg->apply_wfwhm_weights || arg->apply_nbstars_weights) {
				siril_log_message(_("Only one weighting method can be used\n"));
			} else {
				arg->apply_nbstack_weights = TRUE;
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
		} else {
			siril_log_message(_("Unexpected argument to stacking `%s', aborting.\n"), current);
			return CMD_ARG_ERROR;
		}
		first++;
	}
	return CMD_OK;
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
	args.coeff.offset = NULL;
	args.coeff.mul = NULL;
	args.coeff.scale = NULL;
	if (!arg->force_no_norm &&
			(arg->method == stack_median || arg->method == stack_mean_with_rejection))
		args.normalize = arg->norm;
	else args.normalize = NO_NORM;
	args.method = arg->method;
	args.force_norm = FALSE;
	args.output_norm = arg->output_norm;
	args.equalizeRGB = arg->equalizeRGB;
	args.lite_norm = arg->lite_norm;
	args.reglayer = get_registration_layer(args.seq);
	args.apply_noise_weights = arg->apply_noise_weights;
	args.apply_nbstack_weights = arg->apply_nbstack_weights;
	args.apply_wfwhm_weights = arg->apply_wfwhm_weights;
	args.apply_nbstars_weights = arg->apply_nbstars_weights;

	// manage registration data
	if (!stack_regdata_is_valid(&args)) {
		siril_log_color_message(_("Stacking has detected registration data on layer %d with more than simple shifts. You should apply existing registration before stacking\n"), "red", args.reglayer);
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}

	// manage filters
	if (convert_parsed_filter_to_filter(&arg->filters, seq,
				&args.filtering_criterion, &args.filtering_parameter) ||
			setup_filtered_data(&args)) {
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}
	args.description = describe_filter(seq, args.filtering_criterion, args.filtering_parameter);
	args.use_32bit_output = evaluate_stacking_should_output_32bits(args.method,
			args.seq, args.nb_images_to_stack, &error);
	if (error) {
		siril_log_color_message(error, "red");
		free_sequence(seq, TRUE);
		return CMD_GENERIC_ERROR;
	}

	main_stack(&args);

	int retval = args.retval;
	clean_end_stacking(&args);
	free(args.image_indices);
	free(args.description);
	free(args.critical_value);

	if (!retval) {
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
		if (g_mkdir_with_parents(dirname, 0755) < 0) {
			siril_log_color_message(_("Cannot create output folder: %s\n"), "red", dirname);
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

	start_in_new_thread(stackall_worker, arg);
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

	start_in_new_thread(stackone_worker, arg);
	return CMD_OK;

failure:
	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	return CMD_ARG_ERROR;
}

/* preprocess sequencename [-bias=filename|value] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=] [-fitseq]
 * preprocess_single filename [-bias=filename|value] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=]
 */
struct preprocessing_data *parse_preprocess_args(int nb, sequence *seq) {
	int retvalue = 0;
	struct preprocessing_data *args = calloc(1, sizeof(struct preprocessing_data));
	fits reffit = { 0 };
	int bitpix;
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
		if (read_fits_metadata_from_path(word[1], &reffit)) {
			siril_log_message(_("Could not load the image, aborting.\n"));
			clearfits(&reffit);
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
				int status;
				expression = path_parse(&reffit, word[i] + 6, PATHPARSE_MODE_READ, &status);
				if (status) {
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
		} else if (!strcmp(word[i], "-opt")) {
			if (bitpix == BYTE_IMG) {
				siril_log_color_message(_("Dark optimization: This process cannot be applied to 8b images\n"), "red");
				retvalue = CMD_INVALID_IMAGE;
				break;
			}
			args->use_dark_optim = TRUE;
		} else if (!strcmp(word[i], "-fix_xtrans")) {
			args->fix_xtrans = TRUE;
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
				gchar *end1, *end2;
				args->sigma[0] = g_ascii_strtod(word[i + 1], &end1);
				args->sigma[1] = g_ascii_strtod(word[i + 2], &end2);

				if (word[i + 1] && word[i + 2] && word[i + 1] != end1
						&& word[i + 2] != end2) {
					i+= 2;
				}

				if (args->sigma[0] == 0) args->sigma[0] = -1.00;
				if (args->sigma[1] == 0) args->sigma[1] = -1.00;
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
//						g_object_unref(args->bad_pixel_map_file); // This is unreferenced in check_for_cosme_file_sanity
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
				siril_log_message(_("Unknown argument %s, aborting.\n"), word[i + 1]);
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
	if (retvalue) {
		free(args);
		return NULL;
	}
	return args;
}

int process_preprocess(int nb) {
	if (word[1][0] == '\0')
		return CMD_ARG_ERROR;

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	struct preprocessing_data *args = parse_preprocess_args(nb, seq);
	if (!args)
		return CMD_ARG_ERROR;

	siril_log_color_message(_("Preprocessing...\n"), "green");
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway
	args->allow_32bit_output = (args->output_seqtype == SEQ_REGULAR
			|| args->output_seqtype == SEQ_FITSEQ) && !com.pref.force_16bit;

	start_sequence_preprocessing(args);
	return CMD_OK;
}

int process_preprocess_single(int nb) {
	if (word[1][0] == '\0')
		return CMD_ARG_ERROR;

	struct preprocessing_data *args = parse_preprocess_args(nb, NULL);
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
		if (!word[2] || !word[3] || (!g_str_has_prefix(word[2], "-type="))) {
			siril_log_message(_("Please specify the type of compression and quantization value.\n"));
			return CMD_ARG_ERROR;
		}
		gchar *comp = NULL;
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
		if (!word[3]) {
			siril_log_message(_("Please specify the value of quantization.\n"));
			g_free(comp);
			return CMD_ARG_ERROR;
		}
		gchar *end;
		q = g_ascii_strtod(word[3], &end);
		if (end == word[3] || (q == 0.0 && (method == RICE_COMP || method == HCOMPRESS_COMP))) {
			siril_log_message(_("Quantization can only be equal to 0 for GZIP1 and GZIP2 algorithms.\n"));
			return CMD_ARG_ERROR;
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
	if (!com.script)
		update_spinCPU(0);

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
#ifdef SIRIL_UNSTABLE
	siril_log_message("unreleased %s %s-%s for %s (%s)\n", PACKAGE, VERSION, SIRIL_GIT_VERSION_ABBREV,
			SIRIL_BUILD_PLATFORM_FAMILY, CPU_ARCH);
#else
	siril_log_message("%s %s for %s (%s)\n", PACKAGE, VERSION, SIRIL_BUILD_PLATFORM_FAMILY, CPU_ARCH);
#endif
#ifdef _OPENMP
	siril_log_message("OpenMP available (%d %s)\n", com.max_thread,
			ngettext("processor", "processors", com.max_thread));
#else
	siril_log_message("OpenMP unavailable\n");
#endif
	siril_log_message("Detected system available memory: %d MB\n", (int)(get_available_memory() / BYTES_IN_A_MB));
	siril_log_message("Can%s create symbolic links\n", test_if_symlink_is_ok(FALSE) ? "" : "not");
#ifndef HAVE_CV44
	siril_log_message("OpenCV 4.2 used, shift-only registration transformation unavailable\n");
#endif
#ifdef HAVE_LIBCURL
	siril_log_message("Built with libcurl\n");
#endif
#ifdef HAVE_JSON_GLIB
	siril_log_message("Built with json-glib\n");
#endif
//#ifdef HAVE_GLIB_NETWORKING
#ifdef HAVE_WCSLIB
	siril_log_message("Built with wcslib\n");
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
#ifdef HAVE_LIBPNG
	siril_log_message("Can read and write PNG files\n");
#endif
#ifdef HAVE_LIBTIFF
	siril_log_message("Can read and write TIFF and Astro-TIFF files\n");
#endif
#ifdef HAVE_LIBHEIF
	siril_log_message("Can read HEIF files\n");
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
	gtk_main_quit();
	return CMD_OK;
}

int process_extract(int nb) {
	int Nbr_Plan, maxplan, mins;
	gchar *end;

	Nbr_Plan = g_ascii_strtoull(word[1], &end, 10);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (end == word[1] || Nbr_Plan > maxplan){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return CMD_GENERIC_ERROR;
	}

	struct wavelets_filter_data *args = malloc(sizeof(struct wavelets_filter_data));

	args->Type = TO_PAVE_BSPLINE;
	args->Nbr_Plan = Nbr_Plan;
	args->fit = &gfit;
	start_in_new_thread(extract_plans, args);

	return CMD_OK;
}

int process_reloadscripts(int nb){
	return refresh_scripts(FALSE, NULL);
}

int process_requires(int nb) {
	gchar **version, **required;
	gint major, minor, micro;
	gint req_major, req_minor, req_micro;
	gchar *endmaj, *endmin, *endmicro, *endreqmaj, *endreqmin, *endreqmicro;

	version = g_strsplit(PACKAGE_VERSION, ".", 3);
	required = g_strsplit(word[1], ".", 3);

	if (g_strv_length(required) != 3) {
		siril_log_color_message(_("Required version is not correct.\n"), "red");

		g_strfreev(version);
		g_strfreev(required);
		return CMD_GENERIC_ERROR;
	}

	major = g_ascii_strtoull(version[0], &endmaj, 10);
	minor = g_ascii_strtoull(version[1], &endmin, 10);
	micro = g_ascii_strtoull(version[2], &endmicro, 10);

	req_major = g_ascii_strtoull(required[0], &endreqmaj, 10);
	req_minor = g_ascii_strtoull(required[1], &endreqmin, 10);
	req_micro = g_ascii_strtoull(required[2], &endreqmicro, 10);

	if (endmaj == version[0] || endmin == version[1] || endmicro == version[2]
			|| endreqmaj == required[0] || endreqmin == required[1]
			|| endreqmicro == required[2]) {
		siril_log_message(_("Wrong parameters.\n"));
		g_strfreev(version);
		g_strfreev(required);
		return CMD_ARG_ERROR;
	}

	g_strfreev(version);
	g_strfreev(required);

	if ((major > req_major || (major == req_major && minor > req_minor)
			|| (major == req_major && minor == req_minor && micro >= req_micro))) {
		// no need to output something in script conditions
		if (!com.script) {
			siril_log_message(_("The required version of Siril is ok.\n"));
		}
		return CMD_OK;
	} else {
		if (!com.script) {
			siril_log_color_message(_("A newer version of Siril is required, please update your version.\n"), "red");
		} else {
			siril_log_color_message(_("The script you are executing requires a newer version of Siril to run (%s), aborting.\n"), "red", word[1]);
		}
		return CMD_GENERIC_ERROR;
	}
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
	if (x+w > gfit.rx || y+h > gfit.ry) {
		siril_log_message(_("The provided coordinates are outside the dimension of the currently loaded image (%d x %d).\n"), gfit.rx, gfit.ry);
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
		new_selection_zone();
	return CMD_OK;
}

int process_rgbcomp(int nb) {
	fits r = { 0 }, g = { 0 }, b = { 0 };
	fits rgb = { 0 }, *rgbptr = &rgb;
	int retval = 0, next_arg;
	char *default_result_name;

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
		if (had_an_rgb_image)
			merge_fits_headers_to_result(rgbptr, &l, &r, NULL);
		else merge_fits_headers_to_result(rgbptr, &l, &r, &g, &b, NULL);
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
		default_result_name = "composed_lrgb";
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
		merge_fits_headers_to_result(rgbptr, &r, &g, &b, NULL);
		rgbptr->history = g_slist_append(rgbptr->history, strdup("RGB composition"));
		size_t nbpix = r.naxes[0] * r.naxes[1];
		for (size_t i = 0; i < nbpix; i++) {
			rgb.fpdata[RLAYER][i] = r.fdata[i];
			rgb.fpdata[GLAYER][i] = g.fdata[i];
			rgb.fpdata[BLAYER][i] = b.fdata[i];
		}
		next_arg = 4;
		default_result_name = "composed_rgb";
	}

	clearfits(&r); clearfits(&g); clearfits(&b);
	gchar *result_filename;
	if (nb == next_arg + 1 && g_str_has_prefix(word[next_arg], "-out=") &&
			word[next_arg][5] != '\0') {
		result_filename = word[next_arg] + 5;
		if (g_str_has_suffix(result_filename, com.pref.ext))
			result_filename = g_strdup(result_filename);
		else result_filename = g_strdup_printf("%s%s", result_filename, com.pref.ext);
	} else result_filename = g_strdup_printf("%s%s", default_result_name, com.pref.ext);

	retval = savefits(result_filename, rgbptr);
	g_free(result_filename);
	clearfits(rgbptr);

	return retval;
}

// used for plate solver and PCC commands
int process_pcc(int nb) {
	gboolean noflip = FALSE, plate_solve = !has_wcs(&gfit), downsample = FALSE;
	SirilWorldCS *target_coords = NULL;
	double forced_focal = -1.0, forced_pixsize = -1.0;
	double mag_offset = 0.0, target_mag = -1.0;
	online_catalog cat = CAT_AUTO;
	gboolean pcc_command = word[0][1] == 'c'; // not 'platesolve' or 'seqplatesolve'
	gboolean seqps = word[0][0] == 's';
	sequence *seq = NULL;

	gboolean local_cat = local_catalogues_available();
	int next_arg = 1;
	if (seqps) {
		if (!(seq = load_sequence(word[1], NULL)))
			return CMD_SEQUENCE_NOT_FOUND;
		next_arg++;
	}

	// check if we have target_coords
	if (nb > next_arg && (word[next_arg][0] != '-' || (word[next_arg][1] >= '0' && word[next_arg][1] <= '9'))) {
		char *sep = strchr(word[next_arg], ',');
		if (!sep) {
			if (nb == 2) {
				siril_log_message(_("Could not parse target coordinates\n"));
				return CMD_ARG_ERROR;
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
			return CMD_ARG_ERROR;
		}
	}

	while (nb > next_arg && word[next_arg]) {
		if (!strcmp(word[next_arg], "-noflip"))
			noflip = TRUE;
		else if (!strcmp(word[next_arg], "-platesolve"))
			plate_solve = TRUE;
		else if (!strcmp(word[next_arg], "-downscale"))
			downsample = TRUE;
		else if (g_str_has_prefix(word[next_arg], "-focal=")) {
			char *arg = word[next_arg] + 7;
			gchar *end;
			forced_focal = g_ascii_strtod(arg, &end);
			if (end == arg || forced_focal <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				if (target_coords)
					siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[next_arg], "-pixelsize=")) {
			char *arg = word[next_arg] + 11;
			gchar *end;
			forced_pixsize = g_ascii_strtod(arg, &end);
			if (end == arg || forced_pixsize <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				if (target_coords)
					siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[next_arg], "-limitmag=")) {
			char *arg = word[next_arg] + 10;
			gchar *end;
			double value;
			value = g_ascii_strtod(arg, &end);
			if (end == arg) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				if (target_coords)
					siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
			if (arg[0] == '-' || arg[0] == '+')
				mag_offset = value;
			else target_mag = value;
		}
		else if (g_str_has_prefix(word[next_arg], "-catalog=")) {
			char *arg = word[next_arg] + 9;
			if (!g_strcmp0(arg, "tycho2"))
				cat = CAT_TYCHO2;
			else if (!g_strcmp0(arg, "nomad"))
				cat = CAT_NOMAD;
			else if (!g_strcmp0(arg, "gaia"))
				cat = CAT_GAIADR3;
			else if (!g_strcmp0(arg, "ppmxl"))
				cat = CAT_PPMXL;
			else if (!g_strcmp0(arg, "brightstars"))
				cat = CAT_BRIGHT_STARS;
			else if (!g_strcmp0(arg, "apass"))
				cat = CAT_APASS;
			else {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				if (target_coords)
					siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
			if (pcc_command && cat != CAT_NOMAD && cat != CAT_APASS) {
				siril_log_message(_("Catalog can only be NOMAD or APASS for photometric usage\n"));
				if (target_coords)
					siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
		}
		else if (!pcc_command && !g_ascii_strcasecmp(word[next_arg], "-localasnet")) {
			if (cat != CAT_AUTO)
				siril_log_message(_("Specifying a catalog has no effect for astrometry.net solving\n"));
#ifndef HAVE_WCSLIB
			siril_log_color_message(_("Astrometry.net interaction relies on the missing WCSLIB software, cannot continue.\n"), "red");
			return CMD_ARG_ERROR;
#else
			cat = CAT_ASNET;
			local_cat = TRUE;
#endif
		} else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
			if (target_coords)
				siril_world_cs_unref(target_coords);
			return CMD_ARG_ERROR;
		}
		next_arg++;
	}

	if (local_cat && cat != CAT_AUTO && cat != CAT_ASNET) {
		siril_log_color_message(_("Using remote %s instead of local NOMAD catalogue\n"),
				"salmon", catalog_to_str(cat));
		local_cat = FALSE;
	}

	if (seqps) {
		struct astrometry_data *args = calloc(1, sizeof(struct astrometry_data));
		args->pixel_size = forced_pixsize;
		args->focal_length = forced_focal;
		args->use_local_cat = local_cat;
		args->onlineCatalog = cat;
		args->cat_center = target_coords;
		args->downsample = downsample;
		args->autocrop = TRUE;
		args->flip_image = !noflip;
		args->manual = FALSE;
		if (target_mag > -1.0) {
			args->mag_mode = LIMIT_MAG_ABSOLUTE;
			args->magnitude_arg = target_mag;
		}
		else if (mag_offset != 0.0) {
			args->mag_mode = LIMIT_MAG_AUTO_WITH_OFFSET;
			args->magnitude_arg = mag_offset;
		}
		else args->mag_mode = LIMIT_MAG_AUTO;
		args->for_sequence = TRUE;
		args->verbose = FALSE;
		args->for_photometry_cc = FALSE;

		start_sequence_astrometry(seq, args);
		return CMD_OK;
	}

	if (plate_solve && !target_coords) {
		target_coords = get_eqs_from_header(&gfit);
		if (cat != CAT_ASNET) {
			if (!target_coords) {
				siril_log_color_message(_("Cannot plate solve, no target coordinates passed and image header doesn't contain any either\n"), "red");
				return CMD_INVALID_IMAGE;
			}
		}
		if (target_coords) {
			siril_log_message(_("Using target coordinate from image header: %f, %f\n"),
					siril_world_cs_get_alpha(target_coords),
					siril_world_cs_get_delta(target_coords));
		}
	}

	if (pcc_command) {
		if (plate_solve)
			siril_log_message(_("Photometric color correction will use plate solving first\n"));
		else siril_log_message(_("Photometric color correction will use WCS information and bypass internal plate solver\n"));
	} else if (!plate_solve) {
		siril_log_message(_("Image is already plate solved. Nothing will be done.\n"));
		return CMD_OK;
	}

	struct astrometry_data *args = NULL;	// filled only if running a plate solve
	struct photometric_cc_data *pcc_args = NULL;	// filled only if pcc_command

	if (plate_solve) {
		args = calloc(1, sizeof(struct astrometry_data));
		args->fit = &gfit;
		args->verbose = TRUE;
	}
	if (pcc_command) {
		pcc_args = calloc(1, sizeof(struct photometric_cc_data));
		pcc_args->fit = &gfit;
		pcc_args->bg_auto = TRUE;
	}

	/* populating plate solve arguments */
	if (forced_pixsize > 0.0) {
		if (!plate_solve)
			siril_log_message(_("Focal length and pixel size are only used for plate solving, ignored now\n"));
		else {
			args->pixel_size = forced_pixsize;
			siril_log_message(_("Using provided pixel size: %.2f\n"), args->pixel_size);
		}
	} else if (plate_solve) {
		args->pixel_size = max(gfit.pixel_size_x, gfit.pixel_size_y);
		if (args->pixel_size <= 0.0) {
			args->pixel_size = com.pref.starfinder_conf.pixel_size_x;
			if (args->pixel_size <= 0.0) {
				siril_log_color_message(_("Pixel size not found in image or in settings, cannot proceed\n"), "red");
				if (target_coords)
					siril_world_cs_unref(target_coords);
				g_free(args);
				g_free(pcc_args);
				return CMD_INVALID_IMAGE;
			}
			siril_log_message(_("Using pixel size from preferences: %.2f\n"), args->pixel_size);
		}
		else siril_log_message(_("Using pixel size from image: %.2f\n"), args->pixel_size);
	}
	if (forced_focal > 0.0) {
		if (!plate_solve) {
			if (forced_pixsize <= 0.0)
				siril_log_message(_("Focal length and pixel size are only used for plate solving, ignored now\n"));
		} else {
			args->focal_length = forced_focal;
			siril_log_message(_("Using provided focal length: %.2f\n"), args->focal_length);
		}
	} else if (plate_solve) {
		args->focal_length = gfit.focal_length;
		if (args->focal_length <= 0.0) {
			args->focal_length = com.pref.starfinder_conf.focal_length;
			if (args->focal_length <= 0.0) {
				siril_log_color_message(_("Focal length not found in image or in settings, cannot proceed\n"), "red");
				if (target_coords)
					siril_world_cs_unref(target_coords);
				g_free(args);
				g_free(pcc_args);
				return CMD_INVALID_IMAGE;
			}
			siril_log_message(_("Using focal length from preferences: %.2f\n"), args->focal_length);
		}
		else siril_log_message(_("Using focal length from image: %.2f\n"), args->focal_length);
	}

	if (target_mag > -1.0) {
		if (cat == CAT_ASNET)
			siril_log_message(_("Magnitude alteration arguments are useless for astrometry.net plate solving\n"));
		if (args) {
			args->mag_mode = LIMIT_MAG_ABSOLUTE;
			args->magnitude_arg = target_mag;
		}
		if (pcc_args) {
			pcc_args->mag_mode = LIMIT_MAG_ABSOLUTE;
			pcc_args->magnitude_arg = target_mag;
		}
	}
	else if (mag_offset != 0.0) {
		if (cat == CAT_ASNET)
			siril_log_message(_("Magnitude alteration arguments are useless for astrometry.net plate solving\n"));
		if (args) {
			args->mag_mode = LIMIT_MAG_AUTO_WITH_OFFSET;
			args->magnitude_arg = mag_offset;
		}
		if (pcc_args) {
			pcc_args->mag_mode = LIMIT_MAG_AUTO_WITH_OFFSET;
			pcc_args->magnitude_arg = mag_offset;
		}
	}
	else {
		if (args)
			args->mag_mode = LIMIT_MAG_AUTO;
		if (pcc_args)
			pcc_args->mag_mode = LIMIT_MAG_AUTO;
	}

	if (plate_solve) {
		args->use_local_cat = local_cat;
		if (cat == CAT_ASNET) {
			args->onlineCatalog = CAT_ASNET;
			args->filename = g_strdup(com.uniq->filename);
		}
		else args->onlineCatalog = local_cat ? CAT_NOMAD : cat;
		args->cat_center = target_coords;
		args->downsample = downsample;
		args->autocrop = TRUE;
		args->flip_image = !noflip;
		args->manual = FALSE;
		process_plate_solver_input(args);
	}

	if (pcc_command) {
		pcc_args->use_local_cat = local_cat;
		pcc_args->catalog = local_cat ? CAT_NOMAD : cat;
		if (plate_solve) {
			args->for_photometry_cc = TRUE;
			args->pcc = pcc_args;
		}
	}

	if (plate_solve)
		start_in_new_thread(plate_solver, args);
	else {
		start_in_new_thread(photometric_cc_standalone, pcc_args);
	}

	return CMD_OK;
}

int process_nomad(int nb) {
	float limit_mag = 13.0f;
	gboolean photometric = FALSE;
	online_catalog cat = CAT_AUTO;
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}
	int arg_idx = 1;
	while (arg_idx < nb) {
		if (g_str_has_prefix(word[arg_idx], "-catalog=")) {
			char *arg = word[arg_idx] + 9;
			if (!g_strcmp0(arg, "tycho2"))
				cat = CAT_TYCHO2;
			else if (!g_strcmp0(arg, "nomad"))
				cat = CAT_NOMAD;
			else if (!g_strcmp0(arg, "gaia"))
				cat = CAT_GAIADR3;
			else if (!g_strcmp0(arg, "ppmxl"))
				cat = CAT_PPMXL;
			else if (!g_strcmp0(arg, "brightstars"))
				cat = CAT_BRIGHT_STARS;
			else if (!g_strcmp0(arg, "apass"))
				cat = CAT_APASS;
			else {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[arg_idx], "-phot")) {
			photometric = TRUE;
		} else {
			gchar *end;
			limit_mag = g_ascii_strtod(word[arg_idx], &end);
			if (end == word[arg_idx]) {
				siril_log_message(_("Invalid argument %s, aborting.\n"), word[arg_idx]);
				return CMD_ARG_ERROR;
			}
		}
		arg_idx++;
	}

	pcc_star *stars;
	int nb_stars;
	double ra, dec;
	center2wcs(&gfit, &ra, &dec);
	double resolution = get_wcs_image_resolution(&gfit);
	uint64_t sqr_radius = ((uint64_t) gfit.rx * gfit.rx + (uint64_t) gfit.ry * gfit.ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	siril_debug_print("centre coords: %f, %f, radius: %f\n", ra, dec, radius);

	if (cat == CAT_AUTO) {
		if (get_photo_stars_from_local_catalogues(ra, dec, radius, &gfit, limit_mag, &stars, &nb_stars)) {
			siril_log_color_message(_("Failed to get data from the local catalogue, is it installed?\n"), "red");
			return CMD_GENERIC_ERROR;
		}
		siril_debug_print("Got %d stars from the trixels of this target (mag limit %.2f)\n", nb_stars, limit_mag);
	} else {
		SirilWorldCS *center = siril_world_cs_new_from_a_d(ra, dec);
		GFile *catalog_file = download_catalog(cat, center, radius * 60.0, limit_mag);
		siril_world_cs_unref(center);
		if (!catalog_file) {
			siril_log_message(_("Could not download the online star catalog.\n"));
			return CMD_GENERIC_ERROR;
		}
		siril_log_message(_("The %s catalog has been successfully downloaded.\n"), catalog_to_str(cat));

		/* project using WCS */
		if (project_catalog_with_WCS(catalog_file, &gfit, photometric, &stars, &nb_stars))
			return CMD_GENERIC_ERROR;
	}

	clear_stars_list(FALSE);
	int j = 0;
	for (int i = 0; i < nb_stars && j < MAX_STARS; i++) {
		if (stars[i].x < 0.0 || stars[i].x >= gfit.rx ||
				stars[i].y < 0.0 || stars[i].y >= gfit.ry)
			continue;
		if (!com.stars)
			com.stars = new_fitted_stars(MAX_STARS);
		com.stars[j] = new_psf_star();
		com.stars[j]->xpos = stars[i].x;
		com.stars[j]->ypos = stars[i].y;
		com.stars[j]->fwhmx = 5.0f;
		com.stars[j]->fwhmy = 5.0f;
		com.stars[j]->layer = 0;
		com.stars[j]->angle = 0.0f;
		j++;
	}
	if (j > 0)
		com.stars[j] = NULL;
	free(stars);
	siril_log_message("%d stars found%s in the image (mag limit %.2f)\n", j,
			photometric ? " with valid photometry data" : "", limit_mag);
	redraw(REDRAW_OVERLAY);
	return CMD_OK;
}

int process_sso() {
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}

	if (!gfit.date_obs) {
		siril_log_color_message(_("This command only works on images that have observation date information\n"), "red");
		return CMD_INVALID_IMAGE;
	}

	purge_temp_user_catalogue();
	force_to_refresh_catalogue_list();

	double lim_mag = 20.0;
	struct astrometry_data *args = calloc(1, sizeof(struct astrometry_data));
	args->fit = &gfit;

	char *arg = word[1];
	if (word[1] && g_str_has_prefix(arg, "-mag=")) {
		char *next, *value = arg + 5;
		lim_mag = g_ascii_strtod(value, &next);
		if (next == value) {
			siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
			free(args);
			return CMD_ARG_ERROR;
		}
	}

	args->limit_mag = lim_mag;
	args->focal_length = gfit.focal_length;
	args->pixel_size = gfit.pixel_size_x;
	args->scale = get_resolution(args->focal_length, args->pixel_size);

	start_in_new_thread(search_in_online_conesearch, args);

	return CMD_OK;
}

static gboolean end_process_catsearch(gpointer p) {
	GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
	force_to_refresh_catalogue_list();
	if (!gtk_toggle_tool_button_get_active(button)) {
		gtk_toggle_tool_button_set_active(button, TRUE);
	} else {
		redraw(REDRAW_OVERLAY);
	}
	return end_generic(NULL);
}

int process_catsearch(int nb){
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}
	set_cursor_waiting(TRUE);

	gchar *result;
	if (nb > 2) {
		// replace with build_string_from_words(word+1);
		GString *str = g_string_new(word[1]);
		int i = 2;
		while (i < nb && word[i]) {
			g_string_append_printf(str, " %s", word[i]);
			i++;
		}
		gchar *name = g_string_free(str, FALSE);
		result = search_object(name);
		g_free(name);
	} else {
		result = search_object(word[1]);
	}

	if (result && !parse_buffer(result, 20.0)) {
		siril_add_idle(end_process_catsearch, NULL);
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
	if (gfit.header) { // fits or astrotiff - do not update the header
		expression = path_parse(&gfit, word[1], mode, &status);
	} else {
		expression = update_header_and_parse(&gfit, word[1], mode, FALSE, &status);
	}
	siril_log_message(_("String in: %s\n"), word[1]);
	siril_log_message(_("String out: %s\n"), expression);
	g_free(expression);
	return CMD_OK;
}
