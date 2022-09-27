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

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/initfile.h"
#include "core/preprocess.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"
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
#include "gui/sequence_list.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"
#include "gui/photometric_cc.h"
#include "gui/script_menu.h"
#include "gui/preferences.h"
#include "filters/asinh.h"
#include "filters/banding.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "filters/clahe.h"
#include "filters/cosmetic_correction.h"
#include "filters/deconv.h"
#include "filters/median.h"
#include "filters/mtf.h"
#include "filters/fft.h"
#include "filters/rgradient.h"
#include "filters/saturation.h"
#include "filters/scnr.h"
#include "filters/starnet.h"
#include "filters/wavelets.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
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
	char filename[256];

	strncpy(filename, word[1], 250);
	filename[250] = '\0';

	for (int i = 1; i < nb - 1; ++i) {
		strcat(filename, " ");
		strcat(filename, word[i + 1]);
	}
	expand_home_in_filename(filename, 256);

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
	// TODO: if sequence is loaded in the UI, it needs to be closed first
	// to avoid rewriting again the .seq upon closing

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return CMD_SEQUENCE_NOT_FOUND;

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
	free_sequence(seq, FALSE);
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
	siril_log_message(_("Applying saturation enhancement of %d%%, from level %g * sigma.\n"), round_to_int(args->coeff * 100.0), args->background_factor);

	set_cursor_waiting(TRUE);
	return GPOINTER_TO_INT(enhance_saturation(args));
}

int process_save(int nb){
	gchar *filename = g_strdup(word[1]);
	set_cursor_waiting(TRUE);
	if (!com.script) {
		gfit.lo = gui.lo;
		gfit.hi = gui.hi;
	}
	int retval = savefits(filename, &gfit);
	set_precision_switch();
	set_cursor_waiting(FALSE);
	g_free(filename);
	return retval;
}

int process_savebmp(int nb){
	gchar *filename = g_strdup_printf("%s.bmp", word[1]);

	set_cursor_waiting(TRUE);
	savebmp(filename, &gfit);
	set_cursor_waiting(FALSE);
	g_free(filename);
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

	siril_log_message("%s\n", log_msg); // This is the standard non-translated message to make things easy for log parsers.
	free(log_msg);
	gettimeofday(&t_start, NULL);
	set_progress_bar_data(_("Starting NL-Bayes denoising..."), 0.0);

	int retval = do_nlbayes(args->fit, args->modulation, args->sos, args->da3d, args->rho, args->do_anscombe, args->do_cosme);

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
		else if (g_str_has_prefix(arg, "-nocosmetic")) {
			args->do_cosme = FALSE;
		}
		else if (g_str_has_prefix(arg, "-mod=")) {
			arg += 5;
			float mod = g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((mod < 0.f) || (mod > 1.f)) {
				siril_log_message(_("Error in mod parameter: must be between 0 and 1, aborting.\n"));
				return CMD_ARG_ERROR;
			}
			if (!error) {
				args->modulation = mod;
			}
		}
		if (args->modulation == 0.f) {
			siril_log_message(_("Modulation is zero: doing nothing.\n"));
			return CMD_OK;
		}
		else if (g_str_has_prefix(arg, "-rho=")) {
			arg += 5;
			float rho = g_ascii_strtod(arg, &end);
			if (arg == end) error = TRUE;
			else if ((rho <= 0.f) || (rho >= 1.f)) {
				siril_log_message(_("Error in rho parameter: must be strictly > 0 and < 1, aborting.\n"));
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
	starnet_args->customstride = FALSE;
	starnet_args->upscale = FALSE;
	starnet_args->starmask = TRUE;
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
				return CMD_ARG_ERROR;
			}
			if (!error) {
				sprintf(starnet_args->stride, "%d", intstride);
				starnet_args->customstride = TRUE;
			}
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
			return CMD_ARG_ERROR;
		}
		if (error) {
			siril_log_message(_("Error parsing arguments, aborting.\n"));
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

	set_cursor_waiting(TRUE);
	savejpg(filename, &gfit, quality);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return CMD_OK;
}
#endif

#ifdef HAVE_LIBPNG
int process_savepng(int nb){

	gchar *filename = g_strdup_printf("%s.png", word[1]);

	set_cursor_waiting(TRUE);
	uint32_t bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
	savepng(filename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return CMD_OK;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	uint16_t bitspersample = 16;

	if (strcasecmp(word[0], "savetif8") == 0)
		bitspersample = 8;
	else if (strcasecmp(word[0], "savetif32") == 0)
		bitspersample = 32;
	gchar *filename = g_strdup_printf("%s.tif", word[1]);
	set_cursor_waiting(TRUE);
	savetif(filename, &gfit, bitspersample, NULL, com.pref.copyright, TRUE);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return CMD_OK;
}
#endif

int process_savepnm(int nb){
	saveNetPBM(word[1], &gfit);
	return CMD_OK;
}

int process_imoper(int nb){
	fits fit = { 0 };

	if (readfits(word[1], &fit, NULL, !com.pref.force_16bit)) return CMD_INVALID_IMAGE;

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
			clearfits(&fit);
			return CMD_ARG_ERROR;
	}
	int retval = imoper(&gfit, &fit, oper, !com.pref.force_16bit);

	clearfits(&fit);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
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

	float coeff = g_ascii_strtod(word[1], NULL);
	if (coeff <= 0.f) {
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
	notify_gfit_modified();
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

int process_rl(int nb) {
	double sigma, corner;
	int iter;
	gchar *end;

	sigma = g_ascii_strtod(word[1], &end);
	if (end == word[1] || sigma < 0.4 || sigma > 2.0) {
		siril_log_message(_("Sigma must be between [0.4, 2.0]\n"));
		return CMD_ARG_ERROR;
	}

	corner = g_ascii_strtod(word[2], &end);
	if (end == word[2] || corner < -0.5 || corner > 0.5) {
		siril_log_message(_("Corner radius boost must be between [-0.5, 0.5]\n"));
		return CMD_ARG_ERROR;
	}

	iter = g_ascii_strtoull(word[3], &end, 10);
	if (end == word[3] || iter <= 0) {
		siril_log_message(_("Number of iterations must be > 0.\n"));
		return CMD_ARG_ERROR;
	}

	struct deconv_data *args = malloc(sizeof(struct deconv_data));

	args->fit = &gfit;
	if (args->fit->type == DATA_USHORT) {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi;
	} else {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi * USHRT_MAX_DOUBLE;
	}
	args->auto_contrast_threshold = TRUE;
	args->sigma = sigma;
	args->corner_radius = corner;
	args->iterations = (size_t)iter;
	args->auto_limit = TRUE;

	start_in_new_thread(RTdeconv, args);

	return CMD_OK;
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
	notify_gfit_modified();
	return CMD_OK;
}

int process_crop(int nb) {
	if (is_preview_active()) {
		siril_log_message(_("It is impossible to crop the image when a filter with preview session is active. "
						"Please consider to close the filter dialog first.\n"));
		return CMD_GENERIC_ERROR;
	}

	rectangle area;
	if (!com.selection.h || !com.selection.w) {
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
	} else {
		area = com.selection;
	}

	crop(&gfit, &area);
	delete_selected_area();
	reset_display_offset();
	adjust_cutoff_from_updated_gfit();
	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();

	return CMD_OK;
}

int process_cd(int nb) {
	char filename[256];
	int retval;

	g_strlcpy(filename, word[1], 250);

	expand_home_in_filename(filename, 256);
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
		coef[i] = g_ascii_strtod(word[i + 1], NULL);
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
	if (nb <= 0)
	return CMD_ARG_ERROR;

	double BP = 0.0;
	BP = g_ascii_strtod(word[1], NULL);
	if ((BP < 0.0) || (BP > 1.0)) {
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
	ght_compute_params compute_params;
	GHTsetup(&compute_params, 0.0, 0.0, 0.0, 0.0, 0.0, STRETCH_LINEAR);
	apply_linked_ght_to_fits(&gfit, &gfit, params, compute_params, TRUE);

	notify_gfit_modified();
	return CMD_OK;
}

int process_ght(int nb) {
	int stretch_colourmodel = COL_INDEP;
	int arg_offset = 0;
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
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

	double D = g_ascii_strtod(word[1 + arg_offset], NULL);
	if ((D < 0.0) || (D > 10.0)) {
		siril_log_message(_("Stretch factor D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double B = g_ascii_strtod(word[2 + arg_offset],NULL);
	if ((B < -5.0) || (B > 15.0)) {
		siril_log_message(_("Stretch intensity B must be between -5 and +15\n"));
		return CMD_ARG_ERROR;
	}

	double SP = g_ascii_strtod(word[4 + arg_offset],NULL);
	if ((SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("Stretch focal point SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[3 + arg_offset],NULL);
	if ((LP < 0.0) || (LP > SP)) {
		siril_log_message(_("Shadow preservation point LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[5 + arg_offset],NULL);
	if ((HP < SP) || (HP > 1.0)) {
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
	ght_compute_params compute_params;
	GHTsetup(&compute_params, B, D, LP, SP, HP, STRETCH_PAYNE_NORMAL);
	apply_linked_ght_to_fits(&gfit, &gfit, params, compute_params, TRUE);

	notify_gfit_modified();
	return CMD_OK;
}

int process_invght(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
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

	double D = g_ascii_strtod(word[1 + arg_offset], NULL);
	if ((D < 0.0) || (D > 10.0)) {
		siril_log_message(_("Stretch factor D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double B = g_ascii_strtod(word[2 + arg_offset],NULL);
	if ((B < -5.0) || (B > 15.0)) {
		siril_log_message(_("Stretch intensity B must be between -5 and +15\n"));
		return CMD_ARG_ERROR;
	}

	double SP = g_ascii_strtod(word[4 + arg_offset],NULL);
	if ((SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("Stretch focal point SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[3 + arg_offset],NULL);
	if ((LP < 0.0) || (LP > SP)) {
		siril_log_message(_("Shadow preservation point LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[5 + arg_offset],NULL);
	if ((HP < SP) || (HP > 1.0)) {
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
	ght_params params = {B, D, LP, SP, HP, 0.0, STRETCH_PAYNE_INVERSE, stretch_colourmodel, do_red, do_green, do_blue};
	ght_compute_params compute_params;
	GHTsetup(&compute_params, B, D, LP, SP, HP, STRETCH_PAYNE_INVERSE);
	apply_linked_ght_to_fits(&gfit, &gfit, params, compute_params, TRUE);

	notify_gfit_modified();
	return CMD_OK;
}

int process_modasinh(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
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

	double D = g_ascii_strtod(word[arg_offset+1], NULL);
	if ((D < 0.0) || (D > 10.0)) {
		siril_log_message(_("D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float) D);

	double SP = g_ascii_strtod(word[arg_offset+3],NULL);
	if ((SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[arg_offset+2],NULL);
	if ((LP < 0.0) || (LP > SP)) {
		siril_log_message(_("LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[arg_offset+4],NULL);
	if ((HP < SP) || (HP > 1.0)) {
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
	ght_compute_params compute_params;
	GHTsetup(&compute_params, 0.0, D, LP, SP, HP, STRETCH_ASINH);
	apply_linked_ght_to_fits(&gfit, &gfit, params, compute_params, TRUE);

	notify_gfit_modified();
	return CMD_OK;
}

int process_invmodasinh(int nb) {
	gboolean do_red = TRUE;
	gboolean do_green = TRUE;
	gboolean do_blue = TRUE;
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

	double D = g_ascii_strtod(word[arg_offset+1], NULL);
	if ((D < 0.0) || (D > 10.0)) {
		siril_log_message(_("D must be between 0 and 10\n"));
		return CMD_ARG_ERROR;
	}
	D = expm1f((float)D);

	double SP = g_ascii_strtod(word[arg_offset+3],NULL);
	if ((SP < 0.0) || (SP > 1.0)) {
		siril_log_message(_("SP must be between 0 and 1\n"));
		return CMD_ARG_ERROR;
	}

	double LP = g_ascii_strtod(word[arg_offset+2],NULL);
	if ((LP < 0.0) || (LP > SP)) {
		siril_log_message(_("LP must be between 0 and stretch focal point\n"));
		return CMD_ARG_ERROR;
	}

	double HP = g_ascii_strtod(word[arg_offset+4],NULL);
	if ((HP < SP) || (HP > 1.0)) {
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
	ght_compute_params compute_params;
	GHTsetup(&compute_params, 0.0, D, LP, SP, HP, STRETCH_INVASINH);
	apply_linked_ght_to_fits(&gfit, &gfit, params, compute_params, TRUE);

	notify_gfit_modified();
	return CMD_OK;
}

int process_wavelet(int nb) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];

	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;

	const char* tmpdir = g_get_tmp_dir();

	Nbr_Plan = g_ascii_strtoull(word[1], NULL, 10);
	Type_Transform = g_ascii_strtoull(word[2], NULL, 10);

	nb_chan = gfit.naxes[2];
	g_assert(nb_chan <= 3);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return CMD_ARG_ERROR;
	}

	if(Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE){
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
	double a[3] = { 0.0 }, b[3] = { 0.0 };
	double low = g_ascii_strtod(word[2], NULL);
	double high = g_ascii_strtod(word[3], NULL);
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
	int arg_offset = 0;
	if (!strcmp(word[1], "-human")) {
		human_luminance = TRUE;
		arg_offset = 1;
	}
	if (nb <= arg_offset + 1)
		return CMD_WRONG_N_ARG;
	double beta = g_ascii_strtod(word[arg_offset + 1], NULL);
	if (beta < 1.0) {
		siril_log_message(_("Stretch must be greater than or equal to 1\n"));
		return CMD_ARG_ERROR;
	}

	double offset = 0.0;
	gchar *end;
	if (nb > 2 + arg_offset) {
		offset = g_ascii_strtod(word[2 + arg_offset], &end);
		if (end == word[2+arg_offset]) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[2 + arg_offset]);
			return CMD_ARG_ERROR;
		}
	}

	set_cursor_waiting(TRUE);
	asinhlut(&gfit, beta, offset, human_luminance);

	notify_gfit_modified();
	return CMD_OK;
}

int process_clahe(int nb) {
	double clip_limit = g_ascii_strtod(word[1], NULL);

	if (clip_limit <= 0.0) {
		siril_log_message(_("Clip limit must be > 0.\n"));
		return CMD_ARG_ERROR;
	}

	int size = g_ascii_strtoull(word[2], NULL, 10);

	if (size <= 0.0) {
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
				char filename[256];

				g_strlcpy(filename, word[1], 250);
				filename[250] = '\0';
				expand_home_in_filename(filename, 256);
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
			char filename[256];
			for (int image = 0; image < seqs[i]->number; image++) {
				fit_sequence_get_image_filename(seqs[i], image, filename, TRUE);
				list = g_list_append(list, g_build_filename(dir, filename, NULL));
			}
		}
		free(seqpath1); free(seqpath2); g_free(dir);
		siril_change_dir(dest_dir, NULL);	// they're all relative to this one
	}

	char *outseq_name;
	struct ser_struct out_ser;
	struct _convert_data *args;
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
			free(outseq_name);
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
			free(outseq_name);
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
	siril_change_dir(dest_dir, NULL);
	free(dest_dir);
	return retval;
}

int	process_mirrorx(int nb){
	mirrorx(&gfit, TRUE);
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int	process_mirrory(int nb){
	mirrory(&gfit, TRUE);
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_mtf(int nb) {
	struct mtf_params params;
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

	apply_linked_mtf_to_fits(&gfit, &gfit, params, TRUE);

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

	notify_gfit_modified();
	return CMD_OK;
}

int process_resample(int nb) {
	gchar *end;
	double factor = g_ascii_strtod(word[1], &end);
	if (end == word[1] || factor <= 0.0 || factor > 5.0) {
		siril_log_message(_("The scaling factor must be less than 5.0\n"));
		return CMD_ARG_ERROR;
	}
	int toX = round_to_int(factor * gfit.rx);
	int toY = round_to_int(factor * gfit.ry);

	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, OPENCV_AREA);

	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	if (!com.script) update_MenuItem();
	return CMD_OK;
}

int process_rgradient(int nb) {
	if (gfit.orig_bitpix == BYTE_IMG) {
		siril_log_color_message(_("This process cannot be applied to 8b images\n"), "red");
		return CMD_INVALID_IMAGE;
	}

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = g_ascii_strtod(word[1], NULL);
	args->yc = g_ascii_strtod(word[2], NULL);
	args->dR = g_ascii_strtod(word[3], NULL);
	args->da = g_ascii_strtod(word[4], NULL);
	args->fit = &gfit;

	if ((args->xc >= args->fit->rx) || (args->yc >= args->fit->ry)) {
		siril_log_message(_("The coordinates cannot be greater than the size of the image. "
				"Please change their values and retry.\n"));
		free(args);
		return CMD_ARG_ERROR;
	}

	start_in_new_thread(rgradient_filter, args);
	return CMD_OK;
}

int process_rotate(int nb) {
	set_cursor_waiting(TRUE);
	int crop = 1;
	gboolean has_selection = FALSE;
	rectangle area = { 0, 0, gfit.rx, gfit.ry };
	if (com.selection.w > 0 && com.selection.h > 0) {
		siril_log_color_message(_("Rotation will apply only to current selection, the resulting image will be cropped.\n"), "salmon");
		area = com.selection;
		has_selection = TRUE;
	}

	double degree = g_ascii_strtod(word[1], NULL);

	/* check for options */
	if (word[2] && (!strcmp(word[2], "-nocrop"))) {
		if (has_selection) {
			siril_log_color_message(_("-nocrop option is not valid if a selection is active. Ignoring\n"), "red");
		} else crop = 0;
	}

	verbose_rotate_image(&gfit, area, degree, OPENCV_AREA, crop);

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
	char *input = word[1];
	if (input[0] == '-') {
		if (word[0][0] == 'g') {
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
	if (word[0][0] == 'g') {
		print_settings_key(input, input+sep+1, FALSE);
	} else {
		/* set */
		char fakefile[1024];
		int filelen = snprintf(fakefile, 1024, "[%s]\n%s\n", input, input+sep+1);
		GKeyFile *kf = g_key_file_new();
		g_key_file_load_from_data(kf, fakefile, filelen, G_KEY_FILE_NONE, NULL);
		return read_keyfile(kf);
	}
	return 0;
}

int process_set_mag(int nb) {
	if (gui.cvport >= MAXGRAYVPORT) {
		siril_log_color_message(_("Please display the channel on which you set the reference magnitude\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	double mag_reference = g_ascii_strtod(word[1], NULL);

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
		psf_star *result = psf_get_minimisation(&gfit, gui.cvport, &com.selection, FALSE, TRUE, ps, TRUE, &error);
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
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	int n = g_ascii_strtoull(word[2], NULL, 10) - 1;
	if (n < 0 || n > seq->number) {
		siril_log_message(_("The reference image must be set between 1 and %d\n"), seq->number);
		return CMD_ARG_ERROR;
	}

	seq->reference_image = n;
	// a reference image should not be excluded to avoid confusion
	if (!seq->imgparam[seq->current].incl) {
		seq->imgparam[seq->current].incl = TRUE;
	}

	writeseqfile(seq);

	return CMD_OK;
}

int process_unset_mag(int nb) {
	com.magOffset = 0.0;
	return CMD_OK;
}

int process_set_mag_seq(int nb) {
	double mag = g_ascii_strtod(word[1], NULL);
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
	double roundness = com.pref.starfinder_conf.roundness;
	int radius = com.pref.starfinder_conf.radius;
	gboolean adjust = com.pref.starfinder_conf.adjust;
	double focal_length = com.pref.starfinder_conf.focal_length;
	double pixel_size_x = com.pref.starfinder_conf.pixel_size_x;
	gboolean relax_checks = com.pref.starfinder_conf.relax_checks;
	int convergence = com.pref.starfinder_conf.convergence;
	gchar *end;

	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-radius=")) {
					char *current = word[i], *value;
					value = current + 8;
					radius = g_ascii_strtoull(value, &end, 10);
					if (end == value || radius < 3 || radius > 50) {
						siril_log_message(_("Wrong parameter values. Radius must be between 3 and 50, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-sigma=")) {
					char *current = word[i], *value;
					value = current + 7;
					sigma = g_ascii_strtod(value, &end);
					if (end == value || sigma < 0.05) {
						siril_log_message(_("Wrong parameter values. Sigma must be greater than 0.05, aborting\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-roundness=")) {
					char *current = word[i], *value;
					value = current + 11;
					roundness = g_ascii_strtod(value, &end);
					if (end == value || roundness < 0.0 || roundness > 0.95) {
						siril_log_message(_("Wrong parameter values. Roundness must be between 0 and 0.95, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-focal=")) {
					char *current = word[i], *value;
					value = current + 7;
					focal_length = g_ascii_strtod(value, &end);
					if (end == value || focal_length < 0) {
						siril_log_message(_("Wrong parameter values. Focal length must be greater than 0, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-pixelsize=")) {
					char *current = word[i], *value;
					value = current + 11;
					pixel_size_x = g_ascii_strtod(value, &end);
					if (end == value || pixel_size_x < 0) {
						siril_log_message(_("Wrong parameter values. Pixel size must be greater than 0, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-auto=")) {
					char *current = word[i], *value;
					value = current + 6;
					if (!(g_ascii_strcasecmp(value, "on"))) adjust = TRUE;
					else if (!(g_ascii_strcasecmp(value, "off"))) adjust = FALSE;
					else {
						siril_log_message(_("Wrong parameter values. Auto must be set to on or off, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-relax=")) {
					char *current = word[i], *value;
					value = current + 7;
					if (!(g_ascii_strcasecmp(value, "on"))) relax_checks = TRUE;
					else if (!(g_ascii_strcasecmp(value, "off"))) relax_checks = FALSE;
					else {
						siril_log_message(_("Wrong parameter values. Auto must be set to on or off, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (g_str_has_prefix(word[i], "-convergence=")) {
					char *current = word[i], *value;
					value = current + 13;
					convergence = g_ascii_strtoull(value, &end, 10);
					if (end == value || convergence < 1 || convergence > 3) {
						siril_log_message(_("Wrong parameter values. Convergence must be between 1 and 3, aborting.\n"));
						return CMD_ARG_ERROR;
					}
				} else if (!g_ascii_strcasecmp(word[i], "reset")) {
					siril_log_message(_("Resetting findstar parameters to default values.\n"));
					sigma = 1.;
					roundness = 0.5;
					radius = 10;
					adjust = TRUE;
					focal_length = 0.;
					pixel_size_x = 0.;
					relax_checks = FALSE;
					convergence = 1;
				} else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
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
	siril_log_message(_("convergence = %d\n"), convergence);
	com.pref.starfinder_conf.convergence = convergence;
	siril_log_message(_("focal = %3.1f\n"), focal_length);
	com.pref.starfinder_conf.focal_length = focal_length;
	siril_log_message(_("pixelsize = %3.2f\n"), pixel_size_x);
	com.pref.starfinder_conf.pixel_size_x = pixel_size_x;
	siril_log_message(_("auto = %s\n"), (adjust) ? "on" : "off");
	com.pref.starfinder_conf.adjust = adjust;
	siril_log_message(_("relax = %s\n"), (relax_checks) ? "on" : "off");
	com.pref.starfinder_conf.relax_checks = relax_checks;

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
	args->varname = malloc(args->nb_rows * sizeof(gchar *));

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

	if (gfit.naxes[2] > 1 && nb == 1 && (com.headless || gui.cvport >= MAXGRAYVPORT)) {
		siril_log_color_message(_("Please display the channel on which you want to compute the PSF or use -channel argument\n"), "red");
		return CMD_GENERIC_ERROR;
	}

	int channel = com.headless ? 0 : gui.cvport;
	if (nb == 2) {
		char *next;
		channel = g_ascii_strtoull(word[1], &next, 10);
		if (word[1] == next || channel > (int)gfit.naxes[2]) {
			siril_log_message(_("Please provide the channel number starting from 0 for red\n"));
			return CMD_ARG_ERROR;
		}
	}

	psf_error error;
	struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
	psf_star *result = psf_get_minimisation(&gfit, channel, &com.selection, TRUE, TRUE, ps, TRUE, &error);
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

	sequence *seq;
	int layer;
	if (nb < 4) {
		seq = &com.seq;
		layer = gui.cvport;

		if (gui.cvport >= MAXGRAYVPORT) {
			siril_log_color_message(_("Please display the channel on which you want to compute the PSF\n"), "red");
			return CMD_GENERIC_ERROR;
		}
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
	return seqpsf(seq, layer, FALSE, FALSE, framing, TRUE, FALSE);
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
		if (g_ascii_strtoull(word[2], NULL, 10) < 0 || g_ascii_strtoull(word[3], NULL, 10) < 0) {
			siril_log_message(_("Crop: x and y must be positive values.\n"));
			return CMD_ARG_ERROR;
		}
		if (g_ascii_strtoull(word[4], NULL, 10) <= 0 || g_ascii_strtoull(word[5], NULL, 10) <= 0) {
			siril_log_message(_("Crop: width and height must be greater than 0.\n"));
			return CMD_ARG_ERROR;
		}
		area.x = g_ascii_strtoull(word[2], NULL, 10);
		area.y = g_ascii_strtoull(word[3], NULL, 10);
		area.w = g_ascii_strtoull(word[4], NULL, 10);
		area.h = g_ascii_strtoull(word[5], NULL, 10);
	} else {
		siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
		return CMD_ARG_ERROR;
	}

	if (g_ascii_strtoull(word[4], NULL, 10) > seq->rx || g_ascii_strtoull(word[5], NULL, 10) > seq->ry) {
		siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"),
				seq->rx, seq->ry);
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
	int nlayer = g_ascii_strtoull(word[1], NULL, 10);
	const gchar* clayer;

	if (nlayer > 3 || nlayer < 0)
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
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int lo = g_ascii_strtoull(word[1], NULL, 10);
	if (lo < 0 || lo > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	int hi = g_ascii_strtoull(word[2], NULL, 10);
	if (hi < 0 || hi > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	if (lo >= hi) {
		siril_log_message(_("lo must be strictly smaller than hi\n"));
		return CMD_ARG_ERROR;
	}
	threshlo(&gfit, lo);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_threshlo(int nb){
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int lo = g_ascii_strtoull(word[1], NULL, 10);
	if (lo < 0 || lo > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	threshlo(&gfit, lo);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_threshhi(int nb){
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	int hi = g_ascii_strtoull(word[1], NULL, 10);
	if (hi < 0 || hi > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_neg(int nb) {
	set_cursor_waiting(TRUE);
	pos_to_neg(&gfit);
	update_gfit_histogram_if_needed();
	invalidate_stats_from_fit(&gfit);
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return CMD_OK;
}

int process_nozero(int nb){
	int level = g_ascii_strtoull(word[1], NULL, 10);
	int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
	if (level < 0 || level > maxlevel) {
		siril_log_message(_("Replacement value is out of range (0 - %d)\n"), maxlevel);
		return CMD_ARG_ERROR;
	}
	nozero(&gfit, (WORD)level);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_ddp(int nb) {
	unsigned level = g_ascii_strtoull(word[1], NULL, 10);
	float coeff = g_ascii_strtod(word[2], NULL);
	float sigma = g_ascii_strtod(word[3], NULL);
	ddp(&gfit, level, coeff, sigma);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

int process_new(int nb){
	int width, height, layers;

	width = g_ascii_strtod(word[1], NULL);
	height = g_ascii_strtod(word[2], NULL);
	layers = g_ascii_strtoull(word[3], NULL, 10);
	if (layers != 1 && layers != 3) {
		siril_log_message(_("Number of layers MUST be 1 or 3\n"));
		return CMD_ARG_ERROR;
	}
	if (!height || !width) return CMD_ARG_ERROR;

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
	int low = g_ascii_strtoull(word[1], NULL, 10);
	int high = g_ascii_strtoull(word[2], NULL, 10);
	if ((high > USHRT_MAX) || (low < 0)) {
		siril_log_message(_("Values must be positive and less than %d.\n"), USHRT_MAX);
		return CMD_ARG_ERROR;
	}
	visu(&gfit, low, high);
	return CMD_OK;
}

int process_fill2(int nb){
	int level = g_ascii_strtoull(word[1], NULL, 10);
	rectangle area;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = g_ascii_strtoull(word[2], NULL, 10);
			area.y = g_ascii_strtoull(word[3], NULL, 10);
			area.w = g_ascii_strtoull(word[4], NULL, 10);
			area.h = g_ascii_strtoull(word[5], NULL, 10);
			if ((area.w + area.x > gfit.rx) || (area.h + area.y > gfit.ry)) {
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
			if (args->im.fit->naxes[2] == 1) {  // handling mono case
				siril_log_message(_("This sequence is mono, ignoring layer number.\n"));
				continue;
			}
			value = current + 7;
			gchar *end;
			int layer = g_ascii_strtoull(value, &end, 10);
			if (end == value || layer < 0 || layer > 2) {
				siril_log_message(_("Unknown layer number %s, must be between 0 and 2, will use green layer.\n"), value);
				if (end == value) break;
				else continue;
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
		layer = gui.cvport == RGB_VPORT ? GLAYER : gui.cvport;
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

	struct starfinder_data *args = calloc(1, sizeof(struct starfinder_data));
	int layer;
	if (!com.script) {
		layer = gui.cvport == RGB_VPORT ? GLAYER : gui.cvport;
	} else {
		layer = (gfit.naxes[2] > 1) ? GLAYER : RLAYER;
	}
	// initializing findstar args
	args->layer = layer;
	args->im.fit = NULL;
	args->im.from_seq = seq;
	args->im.index_in_seq = -1;
	args->max_stars_fitted = 0;
	args->update_GUI = FALSE;
	args->save_to_file = TRUE;

	cmd_errors argparsing = parse_findstar(args, 2, nb);

	if (argparsing) {
		if (args->starfile) g_free(args->starfile);
		free(args);
		return argparsing;
	}

	apply_findstar_to_sequence(args);
	return 0;
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
	g_free(filename);

	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			siril_log_message(_("Could not open file: %s\n"), filename);
		}
		g_object_unref(file);
		return CMD_FILE_NOT_FOUND;
	}

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
		free_sequence(seq, TRUE);
		g_free(filename);
		return CMD_FILE_NOT_FOUND;
	}

	GFile *file = g_file_new_for_path(filename);

	struct cosme_data *args = malloc(sizeof(struct cosme_data));

	if (g_str_has_prefix(word[3], "-prefix=")) {
		char *current = word[3], *value;
		value = current + 8;
		if (value[0] == '\0') {
			free_sequence(seq, TRUE);
			g_free(filename);
			g_object_unref(file);
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
	args->ksize = g_ascii_strtoull(word[1], NULL, 10);
	args->amount = g_ascii_strtod(word[2], NULL);
	args->iterations = 1;

	if (!(args->ksize & 1) || args->ksize < 2 || args->ksize > 15) {
		siril_log_message(_("The size of the kernel MUST be odd and in the range [3, 15].\n"));
		free(args);
		return CMD_ARG_ERROR;
	}
	if (args->amount < 0.0 || args->amount > 1.0) {
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

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = g_ascii_strtoull(word[2], NULL, 10);
			area.y = g_ascii_strtoull(word[3], NULL, 10);
			area.w = g_ascii_strtoull(word[4], NULL, 10);
			area.h = g_ascii_strtoull(word[5], NULL, 10);
			if ((area.w + area.x > gfit.rx) || (area.h + area.y > gfit.ry)) {
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
	int level = g_ascii_strtoull(word[1], NULL, 10);
	int retval = fill(&gfit, level, &area);
	if (retval) {
		siril_log_message(_("Wrong parameters.\n"));
		return CMD_ARG_ERROR;
	}
	redraw(REMAP_ALL);
	return CMD_OK;
}

int process_offset(int nb){
	int level = g_ascii_strtod(word[1], NULL);
	off(&gfit, (float)level);
	adjust_cutoff_from_updated_gfit();
	redraw(REMAP_ALL);
	redraw_previews();
	return CMD_OK;
}

/* The version in command line is a minimal version
 * Only neutral type are available (no amount needed),
 * then we always preserve the lightness */
int process_scnr(int nb){
	struct scnr_data *args = malloc(sizeof(struct scnr_data));

	args->type = 0;
	if (nb > 1) {
		gchar *end;
		args->type = g_ascii_strtoull(word[1], &end, 10);
		if (end == word[1] || args->type > 1) {
			siril_log_message(_("Type can either be 0 (average) or 1 (maximum) neutral protection\n"));
			free(args);
			return CMD_ARG_ERROR;
		}
	}
	args->fit = &gfit;
	args->amount = 0.0;
	args->preserve = TRUE;

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

	args->amount = g_ascii_strtod(word[1], NULL);
	args->sigma = g_ascii_strtod(word[2], NULL);
	args->protect_highlights = TRUE;
	args->fit = &gfit;

	start_in_new_thread(BandingEngineThreaded, args);

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
		interp = BACKGROUND_INTER_POLY;
		degree = g_ascii_strtoull(word[arg_index], NULL, 10);
		if (degree < 1 || degree > 4) {
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
				return CMD_ARG_ERROR;
			}
			prefix = strdup(value);
		}
		else if (g_str_has_prefix(arg, "-samples=")) {
			char *value = arg + 9;
			samples = g_ascii_strtoull(value, NULL, 10);
			if (samples <= 1) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-tolerance=")) {
			char *next, *value = arg + 11;
			tolerance = g_ascii_strtod(value, &next);
			if (next == value || tolerance < 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(arg, "-smooth=")) {
			char *next, *value = arg + 8;
			smooth = g_ascii_strtod(value, &next);
			if (next == value || smooth < 0.0 || smooth > 1.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), arg);
				return CMD_ARG_ERROR;
			}
			if (interp != BACKGROUND_INTER_RBF)
				siril_log_color_message(_("smooth parameter is unused with the polynomial model, ignoring.\n"), "salmon");
		}
		else {
			siril_log_message(_("Unknown parameter %s, aborting.\n"), arg);
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

	args->seq = seq;
	args->sigma[0] = g_ascii_strtod(word[1 + i], NULL);
	args->sigma[1] = g_ascii_strtod(word[2 + i], NULL);
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
	int from = g_ascii_strtoull(word[1], NULL, 10);
	int to = g_ascii_strtoull(word[2], NULL, 10);
	if (from < 1 || from >= com.seq.number) {
		siril_log_message(_("The first argument must be between 1 and the number of images.\n"));
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

	args->type = 0;
	args->str_type = _("RGB");

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
	} else return CMD_INVALID_IMAGE;

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
	} else return CMD_INVALID_IMAGE;

	g_free(Ha);
	clearfits(&f_Ha);
	free(filename);
	return ret;
}

int process_extractHaOIII(int nb) {
	char *filename = NULL;
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

	sensor_pattern pattern = get_bayer_pattern(&gfit);

	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	gchar *OIII = g_strdup_printf("OIII_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractHaOIII_ushort(&gfit, &f_Ha, &f_OIII, pattern))) {
			ret = save1fits16(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractHaOIII_float(&gfit, &f_Ha, &f_OIII, pattern))) {
			ret = save1fits32(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	} else return CMD_INVALID_IMAGE;

	g_free(Ha);
	g_free(OIII);
	clearfits(&f_Ha);
	clearfits(&f_OIII);
	free(filename);
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
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
			}
			else {
				siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
				return CMD_ARG_ERROR;
			}
		}
	}

	apply_split_cfa_to_sequence(args);

	return CMD_OK;
}

int process_seq_extractHa(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq) {
		return CMD_SEQUENCE_NOT_FOUND;
	}

	if (seq->nb_layers > 1) {
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
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
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
						return CMD_ARG_ERROR;
					}
					args->seqEntry = strdup(value);
				}
				else {
					siril_log_message(_("Unknown parameter %s, aborting.\n"), word[i]);
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
		return CMD_FOR_CFA_IMAGE;
	}

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->seqEntry = ""; // not used

	apply_extractHaOIII_to_sequence(args);

	return CMD_OK;
}

int process_stat(int nb){
	int nplane;
	int layer;
	char layername[6];

	nplane = gfit.naxes[2];

	for (layer = 0; layer < nplane; layer++) {
		imstats* stat = statistics(NULL, -1, &gfit, layer, &com.selection, STATS_MAIN, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return CMD_GENERIC_ERROR;
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

		if (gfit.type == DATA_USHORT) {
			siril_log_message(
					_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
							"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
					layername, stat->mean, stat->median, stat->sigma,
					stat->avgDev, stat->min, stat->max);
		} else {
			siril_log_message(
					_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
							"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
					layername, stat->mean * USHRT_MAX_DOUBLE,
					stat->median * USHRT_MAX_DOUBLE,
					stat->sigma * USHRT_MAX_DOUBLE,
					stat->avgDev * USHRT_MAX_DOUBLE,
					stat->min * USHRT_MAX_DOUBLE, stat->max * USHRT_MAX_DOUBLE);
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

	struct stat_data *args = calloc(1, sizeof(struct stat_data));

	args->seq = seq;
	args->seqEntry = ""; // not used
	args->csv_name = g_strdup(word[2]);

	if (word[3] && !g_strcmp0(word[3], "main")) {
		args->option = STATS_MAIN;
	} else if (word[3] && !g_strcmp0(word[3], "full")) {
		args->option = STATS_NORM | STATS_MAIN; // adding STATS_MAIN to include also AVGDEV and SQRTBWMV
	} else {
		args->option = STATS_BASIC;
	}
	com.selection = args->selection;

	apply_stats_to_sequence(args);

	return CMD_OK;
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
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gchar *destroot = g_strdup(word[1]);
	sequence_type output = SEQ_REGULAR;
	gboolean debayer = FALSE;

	if (!com.wd) {
		siril_log_message(_("Conversion: no working directory set.\n"));
		return CMD_NO_CWD;
	}

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

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_clear_error(&error);
		return CMD_NO_CWD;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
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
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gchar *destroot = g_strdup(word[1]);

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

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Link: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Link: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		return CMD_GENERIC_ERROR;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
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
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gboolean debayer = FALSE;
	gboolean make_link = TRUE;
	sequence_type output = SEQ_REGULAR;
	gchar *destroot = g_strdup(word[1]);

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

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Convert: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Convert: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		return CMD_GENERIC_ERROR;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
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
	char **files_to_link = glist_to_array(list, &count);

	gchar *str = ngettext("Convert: processing %d FITS file...\n", "Convert: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);


	if (!com.wd) {
		siril_log_message(_("Convert: no working directory set.\n"));
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
	reg_args->interpolation = OPENCV_AREA;

	/* check for options */
	for (int i = 2; i < nb; i++) {
		if (!strcmp(word[i], "-drizzle")) {
			reg_args->x2upscale = TRUE;
		} else if (!strcmp(word[i], "-noout")) {
			reg_args->no_output = TRUE;
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
			if(!g_strcmp0(g_ascii_strdown(value, -1),"shift")) {
#ifdef HAVE_CV44
				reg_args->type = SHIFT_TRANSFORMATION;
				continue;
#else
				siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
				goto terminate_register_on_error;
#endif
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"similarity")) {
				reg_args->type = SIMILARITY_TRANSFORMATION;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"affine")) {
				reg_args->type = AFFINE_TRANSFORMATION;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"homography")) {
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
			int min_pairs = g_ascii_strtoull(value, NULL, 10);
			if (min_pairs < 4) { // using absolute min_pairs required by homography
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
			if(!g_strcmp0(g_ascii_strdown(value, -1),"nearest") || !g_strcmp0(g_ascii_strdown(value, -1),"ne")) {
				reg_args->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"cubic") || !g_strcmp0(g_ascii_strdown(value, -1),"cu")) {
				reg_args->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"lanczos4") || !g_strcmp0(g_ascii_strdown(value, -1),"la")) {
				reg_args->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"linear") || !g_strcmp0(g_ascii_strdown(value, -1),"li")) {
				reg_args->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"none") || !g_strcmp0(g_ascii_strdown(value, -1),"no")) {
				reg_args->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"area") || !g_strcmp0(g_ascii_strdown(value, -1),"ar")) {
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
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
	return CMD_OK;

terminate_register_on_error:
	g_free(reg_args);
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
			if (*end == '%')
				arg->f_fwhm_p = val;
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
			if (*end == '%')
				arg->f_wfwhm_p = val;
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
			if (*end == '%')
				arg->f_round_p = val;
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
			if (*end == '%')
				arg->f_quality_p = val;
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
			if (*end == '%')
				arg->f_bkg_p = val;
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
			if (*end == '%')
				arg->f_nbstars_p = val;
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
	reg_args->interpolation = OPENCV_AREA;
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
		} else if (g_str_has_prefix(word[i], "-interp=")) {
			char *current = word[i], *value;
			value = current + 8;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				goto terminate_register_on_error;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"nearest") || !g_strcmp0(g_ascii_strdown(value, -1),"ne")) {
				reg_args->interpolation = OPENCV_NEAREST;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"cubic") || !g_strcmp0(g_ascii_strdown(value, -1),"cu")) {
				reg_args->interpolation = OPENCV_CUBIC;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"lanczos4") || !g_strcmp0(g_ascii_strdown(value, -1),"la")) {
				reg_args->interpolation = OPENCV_LANCZOS4;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"linear") || !g_strcmp0(g_ascii_strdown(value, -1),"li")) {
				reg_args->interpolation = OPENCV_LINEAR;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"none") || !g_strcmp0(g_ascii_strdown(value, -1),"no")) {
				reg_args->interpolation = OPENCV_NONE;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"area") || !g_strcmp0(g_ascii_strdown(value, -1),"ar")) {
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
			if(!g_strcmp0(g_ascii_strdown(value, -1),"current")) {
				reg_args->framing = FRAMING_CURRENT;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"min")) {
				reg_args->framing = FRAMING_MIN;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"max")) {
				reg_args->framing = FRAMING_MAX;
				continue;
			}
			if(!g_strcmp0(g_ascii_strdown(value, -1),"cog")) {
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

	set_progress_bar_data(_("Registration: Applying existing data"), PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
	return CMD_OK;

terminate_register_on_error:
	g_free(reg_args);
	return CMD_ARG_ERROR;
}

// parse normalization and filters from the stack command line, starting at word `first'
static int parse_stack_command_line(struct stacking_configuration *arg, int first, gboolean rej_options_allowed, gboolean out_allowed) {
	while (word[first]) {
		char *current = word[first], *value;
		if (!strcmp(current, "-nonorm") || !strcmp(current, "-no_norm"))
			arg->force_no_norm = TRUE;
		else if (!strcmp(current, "-output_norm")) {
			arg->output_norm = TRUE;
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
		} else if (!strcmp(current, "-rgb_equal")) {
			if (!rej_options_allowed) {
				siril_log_message(_("RGB equalization is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("RGB equalization is allowed only if normalization has been activated, ignoring.\n"));
			} else {
				arg->equalizeRGB = TRUE;
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
			if (!rej_options_allowed) {
				siril_log_message(_("Fast normalization is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("Fast normalization is allowed only if normalization has been activated, ignoring.\n"));
			} else {
				arg->lite_norm = TRUE;
			}
		} else if (g_str_has_prefix(current, "-norm=")) {
			if (!rej_options_allowed) {
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
		} else if (parse_filter_args(current, &arg->filters)) {
			;
		} else if (g_str_has_prefix(current, "-out=")) {
			if (out_allowed) {
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

				arg->result_file = g_strdup(value);
			}
			else {
				siril_log_message(_("Output filename option is not allowed in this context, ignoring.\n"));
			}
		}
		else {
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
	struct stacking_args args = { 0 };
	gchar *error = NULL;
	if (seq_check_basic_data(seq, FALSE) == -1) {
		free(seq);
		return CMD_GENERIC_ERROR;
	}
	siril_log_message(_("Stacking sequence %s\n"), seq->seqname);
	args.seq = seq;
	args.ref_image = sequence_find_refimage(seq);
	// the three below: used only if method is average w/ rejection
	if (arg->method == stack_mean_with_rejection && (arg->sig[0] != 0.0 || arg->sig[1] != 0.0)) {
		args.sig[0] = arg->sig[0];
		args.sig[1] = arg->sig[1];
		args.type_of_rejection = arg->type_of_rejection;
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
	args.reglayer = args.seq->nb_layers == 1 ? 0 : 1;
	args.apply_noise_weights = arg->apply_noise_weights;
	args.apply_nbstack_weights = arg->apply_nbstack_weights;
	args.apply_wfwhm_weights = arg->apply_wfwhm_weights;
	args.apply_nbstars_weights = arg->apply_nbstars_weights;
	args.equalizeRGB = arg->equalizeRGB;
	args.lite_norm = arg->lite_norm;

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

	if (!arg->result_file) {
		char filename[256];
		char *suffix = g_str_has_suffix(seq->seqname, "_") ||
			g_str_has_suffix(seq->seqname, "-") ? "" : "_";
		snprintf(filename, 256, "%s%sstacked%s",
				seq->seqname, suffix, com.pref.ext);
		arg->result_file = g_strdup(filename);
	}

	main_stack(&args);

	int retval = args.retval;
	clean_end_stacking(&args);
	free_sequence(seq, TRUE);
	free(args.image_indices);
	g_free(args.description);

	if (!retval) {
		bgnoise_async(&args.result, TRUE);
		if (savefits(arg->result_file, &args.result)) {
			siril_log_color_message(_("Could not save the stacking result %s\n"),
					"red", arg->result_file);
			retval = 1;
		}
		else ++arg->number_of_loaded_sequences;
		bgnoise_await();
	}
	clearfits(&args.result);
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

	// stackall { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-incl[uded]]
	// stackall { med | median } [-nonorm, norm=] [-filter-incl[uded]]
	// stackall { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weighted]
	if (!word[1]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 2;
		gboolean allow_rej_options = FALSE;
		if (!strcmp(word[1], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[1], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[1], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[1], "med") || !strcmp(word[1], "median")) {
			arg->method = stack_median;
			allow_rej_options = TRUE;
		} else if (!strcmp(word[1], "rej") || !strcmp(word[1], "mean")) {
			int shift = 1;
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
			if (!word[2 + shift] || !word[3 + shift] || (arg->sig[0] = g_ascii_strtod(word[2 + shift], NULL)) < 0.0
					|| (arg->sig[1] = g_ascii_strtod(word[3 + shift], NULL)) < 0.0) {
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
			allow_rej_options = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_rej_options, FALSE))
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

	// stack seqfilename { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-incl[uded]] [-out=result_filename]
	// stack seqfilename { med | median } [-nonorm, norm=] [-filter-incl[uded]] [-out=result_filename]
	// stack seqfilename { rej | mean } [type_of_rejection] sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-incl[uded]] [-weighted] [-out=result_filename]
	if (!word[2]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 3;
		gboolean allow_rej_options = FALSE;
		if (!strcmp(word[2], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[2], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[2], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[2], "med") || !strcmp(word[2], "median")) {
			arg->method = stack_median;
			allow_rej_options = TRUE;
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
					(arg->sig[0] = g_ascii_strtod(word[3 + shift], NULL)) < 0.0 ||
					(arg->sig[1] = g_ascii_strtod(word[4 + shift], NULL)) < 0.0) {
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
			allow_rej_options = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_rej_options, TRUE))
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
	}
	else {
		if (read_fits_metadata_from_path(word[1], &reffit)) {
			siril_log_message(_("Could not load the image, aborting.\n"));
			clearfits(&reffit);
			retvalue = 1;
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
				if (seq) {
					// loading the sequence reference image's metadata in case $OFFSET is passed in the expression
					int image_to_load = sequence_find_refimage(seq);
					if (seq_read_frame_metadata(seq, image_to_load, &reffit)) {
						siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
						g_free(expression);
						retvalue = CMD_INVALID_IMAGE;
						break;
					}
				}
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
				args->bias = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->bias, NULL, !com.pref.force_16bit)) {
					args->use_bias = TRUE;
					// if input is 8b, we assume 32b master needs to be rescaled
					if ((args->bias->type == DATA_FLOAT) && (bitpix == BYTE_IMG)) {
						soper(args->bias, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
					}
				} else {
					retvalue = CMD_INVALID_IMAGE;
					free(args->bias);
					break;
				}
			}
		} else if (g_str_has_prefix(word[i], "-dark=")) {
			args->dark = calloc(1, sizeof(fits));
			if (!readfits(word[i] + 6, args->dark, NULL, !com.pref.force_16bit)) {
				args->use_dark = TRUE;
				// if input is 8b, we assume 32b master needs to be rescaled
				if ((args->dark->type == DATA_FLOAT) && (bitpix == BYTE_IMG)) {
					soper(args->dark, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
				}
			} else {
				retvalue = CMD_INVALID_IMAGE;
				free(args->dark);
				break;
			}
		} else if (g_str_has_prefix(word[i], "-flat=")) {
			args->flat = calloc(1, sizeof(fits));
			if (!readfits(word[i] + 6, args->flat, NULL, !com.pref.force_16bit)) {
				args->use_flat = TRUE;
				// no need to deal with bitdepth conversion as flat is just a division (unlike darks which need to be on same scale)
			} else {
				retvalue = CMD_INVALID_IMAGE;
				free(args->flat);
				break;
			}
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
				if (word[i + 1] && word[i + 2] && (args->sigma[0] = g_ascii_strtod(word[i + 1], NULL)) < 0.0
						&& (args->sigma[1] = g_ascii_strtod(word[i + 2], NULL)) < 0.0) {
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
						g_object_unref(args->bad_pixel_map_file);
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
		} else if (!g_ascii_strncasecmp(word[2] + 6, "gzip2", 5))  {
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
		q = g_ascii_strtod(word[3], NULL);
		if (q == 0.0 && (method == RICE_COMP || method == HCOMPRESS_COMP)) {
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

	proc_in = g_ascii_strtoull(word[1], NULL, 10);
	proc_max = omp_get_num_procs();
	if (proc_in > proc_max || proc_in < 1) {
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

int process_set_mem(int nb){
	double ratio = g_ascii_strtod(word[1], NULL);
	if (ratio < 0.05 || ratio > 4.0) {
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
	siril_log_message("unreleased %s %s-%s for %s\n", PACKAGE, VERSION, SIRIL_GIT_VERSION_ABBREV,
			SIRIL_BUILD_PLATFORM_FAMILY);
#else
	siril_log_message("%s %s for %s\n", PACKAGE, VERSION, SIRIL_BUILD_PLATFORM_FAMILY);
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
//#ifdef HAVE_GLIB_NETWORKING
#ifdef HAVE_WCSLIB
	siril_log_message("Built with WCSLIB\n");
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

	Nbr_Plan = g_ascii_strtoull(word[1], NULL, 10);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
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

	version = g_strsplit(PACKAGE_VERSION, ".", 3);
	required = g_strsplit(word[1], ".", 3);

	if (g_strv_length(required) != 3) {
		siril_log_color_message(_("Required version is not correct.\n"), "red");

		g_strfreev(version);
		g_strfreev(required);
		return CMD_GENERIC_ERROR;
	}

	major = g_ascii_strtoull(version[0], NULL, 10);
	minor = g_ascii_strtoull(version[1], NULL, 10);
	micro = g_ascii_strtoull(version[2], NULL, 10);

	req_major = g_ascii_strtoull(required[0], NULL, 10);
	req_minor = g_ascii_strtoull(required[1], NULL, 10);
	req_micro = g_ascii_strtoull(required[2], NULL, 10);

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

int process_detect_trail(int nb) {
	//detect_trail [sigma layer minlen]
	//seq_detect_trail seqname [sigma layer minlen]
	gboolean is_sequence;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	is_sequence = (word[0][2] == 'q');

	float ksigma = 1.f;
	int layer = -1, minlen = 100, defminlen = 100; //default min length of 20px to be qualified as a trail
	int nblines;
	int startnb = (is_sequence) ? 1 : 0;

	for (int i = 1; i < nb - startnb; i++) {
		if (i == 1) ksigma = g_ascii_strtod(word[i + startnb], NULL);
		if (i == 2) layer = g_ascii_strtoull(word[i + startnb], NULL, 10);
		if (i == 3) minlen = g_ascii_strtoull(word[i + startnb], NULL, 10);
	}

	if (minlen < defminlen) {
		siril_log_color_message(_("Minimum length of %d pixels may detect many false positives\n"), "salmon", minlen);
	}
	if (minlen < 0) {
		siril_log_color_message(_("Minimum length cannot be negative, ignoring\n"), "salmon");
		minlen = defminlen;
	}
	if (ksigma < 0) {
		//siril_log_color_message(_("ksigma cannot be negative, ignoring\n"), "salmon");
		//ksigma = 1.f;
	}

	if (!is_sequence) {
		clear_stars_list(FALSE);
		if (layer > gfit.naxes[2]) {
			siril_log_color_message(_("Layer %d. does not exist.\n"), "red", layer);
			return 1;
		}
		if (layer == -1) layer = (gfit.naxes[2] == 3) ? 1 : 0;
		imstats* stat = statistics(NULL, -1, &gfit, layer, NULL, STATS_BASIC, TRUE);
		if (!stat) {
			siril_log_color_message(_("Error: statistics computation failed.\n"), "red");
			return 1;
		}

		float threshold;
		if (ksigma >= 0)
			threshold = stat->median + ksigma * stat->bgnoise;
		else threshold = stat->median - ksigma;
		// sanity checks
		if (threshold > stat->max) {
			siril_log_color_message(_("Detection threshold %f is larger than max value %f.\n"), "red", threshold, stat->max);
			free_stats(stat);
			return 1;
		}
		if (threshold < stat->median) {
			siril_log_color_message(_("Detection threshold is lower than median value.\n"), "salmon");
		}

		struct track *tracks;
		nblines = cvHoughLines(&gfit, layer, threshold, minlen, &tracks);

		if (nblines) {
			if (nblines > 2000) nblines = 2000;
			siril_log_message(_("Found %d trail(s) in current frame, displaying start points\n"), nblines);
			com.stars = malloc((2 * nblines + 1) * sizeof(psf_star *));
			for (int i = 0; i < nblines; i++) {
				com.stars[2*i] = new_psf_star();
				com.stars[2*i]->xpos = tracks[i].start.x;
				com.stars[2*i]->ypos = tracks[i].start.y;
				com.stars[2*i]->fwhmx = 8.0;
				com.stars[2*i]->has_saturated = FALSE;
				com.stars[2*i+1] = new_psf_star();
				com.stars[2*i+1]->xpos = tracks[i].end.x;
				com.stars[2*i+1]->ypos = tracks[i].end.y;
				com.stars[2*i+1]->fwhmx = 8.0;
				com.stars[2*i+1]->has_saturated = TRUE;
			}
			com.stars[nblines] = NULL;
			redraw(REDRAW_OVERLAY);
			free(tracks);
		} else {
			siril_log_message(_("No trails found\n"));
		}

		free_stats(stat);
	} else {
		/*
		sequence *seq = load_sequence(word[1], NULL);
		if (!seq) return 1;
		gboolean is_cfa = FALSE;
		switch (seq->type) {
			case SEQ_SER:
				is_cfa = ser_is_cfa(seq->ser_file);
				break;
			case SEQ_REGULAR:
			case SEQ_FITSEQ:
				break;
			default:
				return 1;
		}
		*/

		/* TODO */

	}
	return 0;
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

int process_pcc(int nb) {
	gboolean noflip = FALSE, plate_solve = !has_wcs(&gfit);
	SirilWorldCS *target_coords = NULL;
	double forced_focal = -1.0, forced_pixsize = -1.0;
	// parse args: [target_coords] [-noflip] [-platesolve] [-focal=len] [-pixelsize=ps]

	// check if we have target_coords
	int next_arg = 1;
	if (nb > 1 && (word[1][0] != '-' || (word[1][1] >= '0' && word[1][1] <= '9'))) {
		char *sep = strchr(word[1], ',');
		next_arg = 2;
		if (!sep) {
			target_coords = siril_world_cs_new_from_objct_ra_dec(word[1], word[2]);
			next_arg = 3;
		}
		else {
			*sep++ = '\0';
			target_coords = siril_world_cs_new_from_objct_ra_dec(word[1], sep);
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
		else if (g_str_has_prefix(word[next_arg], "-focal=")) {
			char *arg = word[next_arg] + 7;
			forced_focal = g_ascii_strtod(arg, NULL);
			if (forced_focal <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
		}
		else if (g_str_has_prefix(word[next_arg], "-pixelsize=")) {
			char *arg = word[next_arg] + 11;
			forced_pixsize = g_ascii_strtod(arg, NULL);
			if (forced_pixsize <= 0.0) {
				siril_log_message(_("Invalid argument to %s, aborting.\n"), word[next_arg]);
				siril_world_cs_unref(target_coords);
				return CMD_ARG_ERROR;
			}
		}
		else {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[next_arg]);
				siril_world_cs_unref(target_coords);
			return CMD_ARG_ERROR;
		}
		next_arg++;
	}

	if (plate_solve && !target_coords) {
		target_coords = get_eqs_from_header(&gfit);
		if (!target_coords) {
			siril_log_color_message(_("Cannot plate solve, no target coordinates passed and image header doesn't contain any either\n"), "red");
			return CMD_INVALID_IMAGE;
		}
		siril_log_message(_("Using target coordinate from image header: %f, %f\n"),
				siril_world_cs_get_alpha(target_coords), siril_world_cs_get_delta(target_coords));
	}

	if (plate_solve)
		siril_log_message(_("Photometric color correction will use plate solving first\n"));
	else  siril_log_message(_("Photometric color correction will use WCS information and bypass internal plate solver\n"));

	struct astrometry_data *args = NULL;
	struct photometric_cc_data *pcc_args = calloc(1, sizeof(struct photometric_cc_data));
	if (plate_solve) {
		args = calloc(1, sizeof(struct astrometry_data));
		args->fit = &gfit;
	}
	if (forced_pixsize > 0.0) {
		if (!plate_solve)
			siril_log_message(_("focal length and pixel size are only used for plate solving, ignored now\n"));
		else {
			args->pixel_size = forced_pixsize;
			siril_log_message(_("Using provided pixel size: %f\n"), args->pixel_size);
		}
	} else if (plate_solve) {
		args->pixel_size = max(gfit.pixel_size_x, gfit.pixel_size_y);
		if (args->pixel_size <= 0.0) {
			args->pixel_size = com.pref.starfinder_conf.pixel_size_x;
			//args->pixel_size = com.pref.pitch;
			if (args->pixel_size <= 0.0) {
				siril_log_color_message(_("Pixel size not found in image or in settings, cannot proceed\n"), "red");
				siril_world_cs_unref(target_coords);
				free(args);
				return CMD_INVALID_IMAGE;
			}
			siril_log_message(_("Using pixel size from preferences: %f\n"), args->pixel_size);
		}
		else siril_log_message(_("Using pixel size from image: %f\n"), args->pixel_size);
	}
	if (forced_focal > 0.0) {
		if (!plate_solve) {
			if (forced_pixsize <= 0.0)
				siril_log_message(_("focal length and pixel size are only used for plate solving, ignored now\n"));
		} else {
			args->focal_length = forced_focal;
			siril_log_message(_("Using provided focal length: %f\n"), args->focal_length);
		}
	} else if (plate_solve) {
		args->focal_length = gfit.focal_length;
		if (args->focal_length <= 0.0) {
			// TODO: which one should we use here?
			args->focal_length = com.pref.starfinder_conf.focal_length;
			//args->focal_length = com.pref.focal;
			if (args->focal_length <= 0.0) {
				siril_log_color_message(_("Focal length not found in image or in settings, cannot proceed\n"), "red");
				siril_world_cs_unref(target_coords);
				free(args);
				return CMD_INVALID_IMAGE;
			}
			siril_log_message(_("Using focal length from preferences: %f\n"), args->focal_length);
		}
		else siril_log_message(_("Using focal length from image: %f\n"), args->focal_length);
	}

	if (local_catalogues_available()) {
		siril_debug_print("using local star catalogues\n");
		pcc_args->use_local_cat = TRUE;
	}
	if (plate_solve) {
		args->use_local_cat = pcc_args->use_local_cat;
		args->onlineCatalog = NOMAD;	// could also be APASS if !use_local_cat
		args->for_photometry_cc = TRUE;
		args->cat_center = target_coords;
		args->downsample = gfit.rx > 6000;
		args->autocrop = TRUE;
		args->flip_image = !noflip;
		args->manual = FALSE;
		args->auto_magnitude = TRUE;
		args->pcc = pcc_args;
		process_plate_solver_input(args);
	}

	pcc_args->fit = &gfit;
	pcc_args->bg_auto = TRUE;
	pcc_args->n_channel = CHANNEL_MIDDLE;
	pcc_args->catalog = NOMAD;

	if (plate_solve)
		start_in_new_thread(match_catalog, args);
	else {
		start_in_new_thread(photometric_cc_standalone, pcc_args);
	}
	return CMD_OK;
}

int process_nomad(int nb) {
	pcc_star *stars;
	int nb_stars;
	double ra, dec;
	float limit_mag = 13.0f;
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return 1;
	}
	if (nb == 2) {
		gchar *end;
		limit_mag = g_ascii_strtod(word[1], &end);
		if (end == word[1]) {
			siril_log_message(_("Invalid argument %s, aborting.\n"), word[1]);
			return 1;
		}
	}
	center2wcs(&gfit, &ra, &dec);
	double resolution = get_wcs_image_resolution(&gfit);
	uint64_t sqr_radius = (gfit.rx * gfit.rx + gfit.ry * gfit.ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	siril_debug_print("centre coords: %f, %f, radius: %f\n", ra, dec, radius);
	if (get_stars_from_local_catalogues(ra, dec, radius, &gfit, limit_mag, &stars, &nb_stars)) {
		siril_log_color_message(_("Failed to get data from the local catalogue, is it installed?\n"), "red");
		return 1;
	}

	siril_debug_print("Got %d stars from the trixels of this target (mag limit %.2f)\n", nb_stars, limit_mag);
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
		com.stars[j]->layer = 0;
		j++;
	}
	if (j > 0)
		com.stars[j] = NULL;
	siril_log_message("%d stars from local catalogues found with valid photometry data in the image (mag limit %.2f)\n", j, limit_mag);
	redraw(REDRAW_OVERLAY);
	return 0;
}

int process_start_ls(int nb) {
	// start_ls [-dark=filename] [-flat=filename] [-gradient_removal] [-watch_files]"
	gchar *dark_file = NULL, *flat_file = NULL;
	gboolean use_file_watcher = FALSE, remove_gradient = FALSE;
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
		} else {
			siril_log_message(_("Unknown option provided: %s\n"), word[i]);
			return CMD_ARG_ERROR;
		}
	}

	return start_livestack_from_command(dark_file, flat_file, use_file_watcher, remove_gradient);
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
