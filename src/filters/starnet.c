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
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/colors.h"
#include "algos/extraction.h"
#include "algos/geometry.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "filters/mtf.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/remixer.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "opencv/opencv.h"

#include <unistd.h>
#include <sys/types.h>
#ifndef _WIN32
#include <signal.h>
#include <sys/wait.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "starnet.h"

#ifdef HAVE_LIBTIFF

fits *current_fit = NULL;
gboolean verbose = TRUE;

// Check maximum path length - OSes except for Windows
#ifndef _WIN32
long get_pathmax(void) {
	long pathmax = -1;

	errno = 0;
	pathmax = pathconf("/", _PC_PATH_MAX);
	if (-1 == pathmax) {
		if (0 == errno) {
#define PATHMAX_INFINITE_GUESS 4096
			pathmax = PATHMAX_INFINITE_GUESS;
		} else {
			fprintf(stderr, "pathconf() FAILED, %d, %s\n", errno,
					strerror(errno));
		}
	}
	return pathmax;
}
#endif

// Wrapper for execve
const char *my_argv[64];
static int forkerrors = 0;

#ifndef _WIN32
static pid_t child_pid;

static int exec_prog(const char **argv) {
	pid_t my_pid;
	int status;
	forkerrors = 0;
	int pipe_fds[2];
	int n;
	char buf[0x100] = {0};

	if (pipe(pipe_fds)) {
		perror("pipe creation error");
		forkerrors = 1;
		return 0;
	}

	if (0 == (my_pid = fork())) {
		child_pid = getpid();
		close(pipe_fds[0]);
		dup2(pipe_fds[1], 1);
		fprintf(stdout, "pid:%d\n",child_pid);
		if (-1 == execve(argv[0], (char **)argv , NULL)) {
			perror("child process execve failed");
			forkerrors = 1;
			return 0;
		}
		fflush(stdout);
		perror(argv[0]);
		exit(0);
	} else {
		close(pipe_fds[1]);
		while (0 == waitpid(my_pid , &status , WNOHANG)) {
			usleep(100000);	// Wait for starnet++ to finish before attempting to process the output
			if ((n = read(pipe_fds[0], buf, 0x100)) > 0) {
				buf[n - 1] = 0;
				char *m = strstr(buf, "pid:");
				if (m != 0) {
					m += 4;
					child_pid = atoi(m);
					com.child_is_running = TRUE;
					com.childpid = child_pid;
				} else {
					double value = g_ascii_strtod(buf, NULL);
					if (value != 0.0 && value == value && verbose) { //
						set_progress_bar_data(_("Running Starnet++"), (value / 100));
					}
				}
			}
		}
		com.child_is_running = FALSE;
		com.childpid = 0;
		if (1 != WIFEXITED(status) || 0 != WEXITSTATUS(status)) {
			siril_log_color_message(_("Error: external command %s failed...\n"), "red", argv[0]);
			forkerrors = 1;
			return 0;
		}
	}
	return 0;
}
#else
static gboolean
utf8_charv_to_wcharv(const gchar *const*utf8_charv,
		wchar_t ***wcharv, int *error_index, GError **error) {
	wchar_t **retval = NULL;

	*wcharv = NULL;
	if (utf8_charv != NULL) {
		int n = 0, i;

		while (utf8_charv[n])
			n++;
		retval = g_new(wchar_t*, n + 1);

		for (i = 0; i < n; i++) {
			retval[i] = g_utf8_to_utf16(utf8_charv[i], -1, NULL, NULL, error);
			if (retval[i] == NULL) {
				if (error_index)
					*error_index = i;
				while (i)
					g_free(retval[--i]);
				g_free(retval);
				return FALSE;
			}
		}

		retval[n] = NULL;
	}
	*wcharv = retval;
	return TRUE;
}

static int exec_prog_win32(const char **argv) {
	forkerrors = 0;
	int pipe_fds[2];
	int n;
	char buf[256] = { 0 };
	HANDLE hProcess;
	int nExitCode = STILL_ACTIVE;
	int fdStdOut;
	gboolean has_started = FALSE;
	wchar_t *wargv0, **wargv;
	GError *conv_error = NULL;
	GError *error = NULL;
	gint conv_error_index;

	if (_pipe(pipe_fds, 0x100, O_TEXT)) {
		perror("pipe creation error");
		forkerrors = 1;
		return 0;
	}

	wargv0 = g_utf8_to_utf16 (argv[0], -1, NULL, NULL, &conv_error);
	if (wargv0 == NULL) {
		g_set_error (&error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED, _("Invalid program name: %s"), conv_error->message);
		g_error_free (conv_error);
		return 0;
	}

	if (!utf8_charv_to_wcharv(argv, &wargv, &conv_error_index, &conv_error)) {
		g_set_error(&error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
				_("Invalid string in argument vector at %d: %s"),
				conv_error_index, conv_error->message);
		g_error_free(conv_error);
		g_free(wargv0);

		return 0;
	}

	//from https://docs.microsoft.com/fr-fr/cpp/c-runtime-library/reference/pipe?view=msvc-170
	fdStdOut = _dup(_fileno(stdout));
	if(_dup2(pipe_fds[1], _fileno(stdout)) != 0) return -1;
	_close(pipe_fds[1]);
	hProcess = (HANDLE)_wspawnve(P_NOWAIT, wargv0, (const wchar_t **) wargv, NULL);
	if(_dup2(fdStdOut, _fileno(stdout)) != 0) return -1;
	_close(fdStdOut);

	g_free(wargv0);
	g_strfreev((gchar**) wargv);

	if (hProcess) {
		com.child_is_running = TRUE;
		com.childhandle = (void*) hProcess;
		while (nExitCode == STILL_ACTIVE) {
			usleep(100000);	// Wait for starnet++ to finish before attempting to process the output
			if ((n = _read(pipe_fds[0], buf, 256)) > 0) {
				buf[n - 1] = '\0';
				double value = g_ascii_strtod(buf, NULL);
				if (value != 0.0 && value == value && verbose) { //
					set_progress_bar_data("Running Starnet++", (value / 100));
					if (!has_started) has_started = TRUE;
				}
			}
			if(!GetExitCodeProcess(hProcess,(unsigned long*)&nExitCode))
				return -1;
		}
		com.child_is_running = FALSE;
		com.childhandle = NULL;
	} else {
		siril_log_color_message(_("Error: external command %s failed...\n"), "red", argv[0]);
 		return -1;
 	}
	if (!has_started) return -1;
 	return 0;
}
#endif

gboolean starnet_executablecheck() {
	// Check for starnet executables (pre-v2.0.2 or v2.0.2+)
	gchar* fullpath = g_build_filename(com.pref.starnet_dir, STARNET_BIN, NULL);
	gchar* rgbpath = g_build_filename(com.pref.starnet_dir, STARNET_RGB, NULL);
	gchar* monopath = g_build_filename(com.pref.starnet_dir, STARNET_MONO, NULL);
	if (g_file_test(fullpath, G_FILE_TEST_IS_EXECUTABLE)) {
		g_free(fullpath);
		g_free(rgbpath);
		g_free(monopath);
		return TRUE;
	} else
		g_free(fullpath);
	if ((current_fit->naxes[2] == 3) && (g_file_test(rgbpath, G_FILE_TEST_IS_EXECUTABLE))) {
		g_free(rgbpath);
		g_free(monopath);
		return TRUE;
	} else
		g_free(rgbpath);
	if ((current_fit->naxes[2] == 1 ) && (g_file_test(monopath, G_FILE_TEST_IS_EXECUTABLE))) {
		g_free(monopath);
		return TRUE;
	}
	else {
		g_free(monopath);
		return FALSE;
	}
}

gboolean end_and_call_remixer(gpointer p)
{
	struct remixargs *blendargs = (remixargs *) p;
	toggle_remixer_window_visibility(CALL_FROM_STARNET, blendargs->fit1, blendargs->fit2);
	return end_generic(NULL);
}

/* Starnet++v2 star removal routine */

gpointer do_starnet(gpointer p) {
	verbose = single_image_is_loaded(); // To suppress log messages during seq working
	int retval = 0;
	fits workingfit, fit;
	starnet_data *args = (starnet_data *) p;
	args->follow_on = single_image_is_loaded() ? args->follow_on : FALSE;
	current_fit = args->starnet_fit;
	int orig_x = current_fit->rx, orig_y = current_fit->ry;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	// Only allocate as much space for filenames as required - we determine the max pathlength
#ifndef _WIN32
	long pathmax = get_pathmax();
#else
	long pathmax = MAX_PATH;	// On Windows use of MAX_PATH is fine as it is not
								// a configurable item
#endif
	gchar *currentdir = NULL;
	gchar starnetcommand[19] = "starnet++";
	gchar *temptif = NULL;
	gchar *starlesstif = NULL;
	gchar *starlessfit = NULL;
	gchar *starmaskfit = NULL;
#ifdef _WIN32
	gchar qtemptif[pathmax];
	gchar qstarlesstif[pathmax];
#endif
	char *imagenoext = NULL;
	char *starlessnoext = NULL;
	char *starmasknoext = NULL;
	char *imagenoextorig = NULL;
	char *temp = NULL;
	gchar starnetprefix[10] = "starnet_";
	gchar starlessprefix[10] = "starless_";
	gchar starmaskprefix[10] = "starmask_";
	memset(starnetcommand, 0, sizeof(starnetcommand));
	// Set up paths and filenames
	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		temp = g_path_get_basename(com.uniq->filename);
		imagenoextorig = g_strdup_printf("%s", temp);
		g_free(temp);
	} else if (sequence_is_loaded()) {
		imagenoextorig = g_strdup_printf("%s%.5d", args->seq->seqname, args->imgnumber + 1);
	} else {
		imagenoextorig = g_strdup_printf("image");
	}
	imagenoext = g_strdup(imagenoextorig);
	for (char *c = imagenoextorig, *q = imagenoext;  *c;  ++c, ++q)
        *q = *c == ' ' ? '_' : *c;
	if (g_strcmp0(imagenoext, imagenoextorig) && verbose)
		siril_log_color_message(_("Starnet++: spaces detected in filename. Starnet++ can't handle these so they have been replaced by underscores.\n"), "salmon");
	starlessnoext = g_strdup_printf("%s%s", starlessprefix, imagenoext);
	starmasknoext = g_strdup_printf("%s%s", starmaskprefix, imagenoext);
	imagenoext = g_strdup_printf("%s%s", starnetprefix, imagenoext);
	imagenoext = g_build_filename(com.wd, imagenoext, NULL);
	imagenoext = remove_ext_from_filename(imagenoext);
	temptif = g_strdup_printf("%s.tif", imagenoext);
	if (strlen(temptif) > pathmax) {
		retval = 1;
		siril_log_color_message(_("Error: file path too long!\n"), "red");
		goto CLEANUP3;
	}
	starlesstif = g_build_filename(com.wd, starlessnoext, NULL);
	temp = remove_ext_from_filename(starlesstif);
	g_free(starlesstif);
	starlesstif = g_strdup(temp);
	g_free(temp);
	starlessfit = g_strdup(starlesstif);
	temp = g_strdup_printf("%s.tif", starlesstif);
	g_free(starlesstif);
	starlesstif = g_strdup(temp);
	g_free(temp);
		if (strlen(starlesstif) > pathmax) {
		retval = 1;
		siril_log_color_message(_("Error: file path too long!\n"), "red");
		goto CLEANUP3;
	}

	temp = g_strdup_printf("%s%s", starlessfit, com.pref.ext);
	g_free(starlessfit);
	starlessfit = g_strdup(temp);
	g_free(temp);
	if (strlen(starlessfit) > pathmax) {
		retval = 1;
		siril_log_color_message(_("Error: file path too long!\n"), "red");
		goto CLEANUP3;
	}
	starmaskfit = g_build_filename(com.wd, starmasknoext, NULL);
	if (strlen(starmaskfit) > pathmax) {
		retval = 1;
		siril_log_color_message(_("Error: file path too long!\n"), "red");
		goto CLEANUP3;
	}
	temp = remove_ext_from_filename(starmaskfit);
	g_free(starmaskfit);
	starmaskfit = g_strdup(temp);
	g_free(temp);
	temp = g_strdup_printf("%s%s", starmaskfit, com.pref.ext);
	g_free(starmaskfit);
	starmaskfit = g_strdup(temp);
	g_free(temp);
	// ok, let's start
	if (verbose)
		set_progress_bar_data(_("Starting Starnet++"), PROGRESS_NONE);

	// Store current working directory
	currentdir = g_get_current_dir();

	if (verbose)
		siril_log_color_message(_("Starnet++: running. Please wait...\n"), "green");
	if (args->customstride && verbose) siril_log_message(_("Starnet++: stride = %s...\n"), args->stride);
	if (!args->starmask && verbose) siril_log_message(_("Starnet++: -nostarmask invoked, star mask will not be generated...\n"));

	// Check starnet directory is not NULL - can happen first time the new preference file is loaded
	if (!com.pref.starnet_dir) {
		retval = 1;
		siril_log_color_message(_("Incorrect permissions on the Starnet++ directory: %s\nEnsure it is correctly set in Preferences / Miscellaneous.\n"), "red", com.pref.starnet_dir);
		goto CLEANUP3;
	}

	// Check starnet directory is defined
	retval = g_access(com.pref.starnet_dir, R_OK);
	if (retval) {
		siril_log_color_message(_("Incorrect permissions on the Starnet++ directory: %s\nEnsure it is correctly set in Preferences / Miscellaneous.\n"), "red", com.pref.starnet_dir);
		goto CLEANUP3;
		// Dijkstra might be spinning in his grave but one of the few legitimate uses
		// of gotos is this type of error cleanup and return. Much more readable than a
		// mess of indentation.
	}

	// Change to starnet directory
	retval = g_chdir(com.pref.starnet_dir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to Starnet++ directory.\nEnsure it is set in Preferences / Miscellaneous...\n"), "red");
		goto CLEANUP3;
	}

	// Create a working copy of the image.
	retval = copyfits(current_fit, &workingfit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
	if (retval) {
		siril_log_color_message(_("Error: image copy failed...\n"), "red");
		goto CLEANUP3;
	}

	// If workingfit is type DATA_FLOAT, check maximum value of gfit, if > 1.0 renormalize
	// in order to avoid cropping brightness peaks on saving to 16 bit TIFF
	if (current_fit->type == DATA_FLOAT) {
		int nplane = current_fit->naxes[2];
		double max = 0.0;
		for (int layer = 0; layer < nplane; layer++) {
			imstats* stat = statistics(NULL, -1, current_fit, layer, &com.selection, STATS_MAIN, MULTI_THREADED);
			if (stat->max > max)
				max = stat->max;
			free_stats(stat);
		}
		if (max > 1.0 && verbose) {
			siril_log_message(_("Starnet++: Pixel values exceed 1.0. Rescaling to avoid clipping peaks.\n"));
			soper(&workingfit, max, OPER_DIV, FALSE);
		}
	}

	// If needed, make a second copy for later use in making the star mask.
	if (args->starmask) {
		retval = copyfits(&workingfit, &fit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
		if (retval) {
			siril_log_color_message(_("Error: image copy failed...\n"), "red");
			goto CLEANUP2;
		}
		// If the force_16bit preference is set and we have a 32bit image loaded, replace the buffer
		if (com.pref.force_16bit && fit.type == DATA_FLOAT) {
			const size_t ndata = fit.naxes[0] * fit.naxes[1] * fit.naxes[2];
			fit_replace_buffer(&fit, float_buffer_to_ushort(fit.fdata, ndata), DATA_USHORT);
		}
	}

	// If we need to pre-stretch a linear image, we use the auto histogram stretch with
	// default shadow clipping and target background. This does marginally clip the
	// shadows but generally by less than 0.001% of pixels. The result of starnet using
	// this stretch is much better than either asinh or GHT stretches.
	struct mtf_params params;
	params.do_red = TRUE;
	params.do_green = TRUE;
	params.do_blue = TRUE;
	find_linked_midtones_balance_default(&workingfit, &params);
	if (args->linear) {
		if (verbose)
			siril_log_message(_("Starnet++: linear mode. Applying Midtone Transfer Function (MTF) pre-stretch to image.\n"));
		apply_linked_mtf_to_fits(&workingfit, &workingfit, params, TRUE);
	}

	// Upscale if needed
	if (args->upscale) {
		if (verbose)
			siril_log_message(_("Starnet++: 2x upscaling selected. Upscaling image...\n"));
		retval = cvResizeGaussian(&workingfit, round_to_int(2*orig_x), round_to_int(2*orig_y), OPENCV_AREA, FALSE);
		if (retval) {
			siril_log_color_message(_("Error: image resize failed...\n"), "red");
			goto CLEANUP;
		}
	}

	// Save current stretched image as working 16-bit TIFF (post initial stretch if the image was linear)
	if (savetif(temptif, &workingfit, 16, NULL, com.pref.copyright, TRUE)) {
		siril_log_color_message(_("Error: unable to save working TIFF of original image...\n"), "red");
		goto CLEANUP;
	}

	// Check for starnet executables (pre-v2.0.2 or v2.0.2+)
	if (g_file_test(STARNET_BIN, G_FILE_TEST_IS_EXECUTABLE)) {
		snprintf(starnetcommand, 19, STARNET_BIN);
	} else if ((current_fit->naxes[2] == 3) && (g_file_test(STARNET_RGB, G_FILE_TEST_IS_EXECUTABLE))) {
		snprintf(starnetcommand, 19, STARNET_RGB);
	} else if ((current_fit->naxes[2] == 1 ) && (g_file_test(STARNET_MONO, G_FILE_TEST_IS_EXECUTABLE))) {
		snprintf(starnetcommand, 19, STARNET_MONO);
	}
	else {
		siril_log_color_message(_("No valid executable found in the Starnet++ directory\n"), "red");
		goto CLEANUP;
	}
	my_argv[0] = starnetcommand;
	my_argv[1] = temptif;
	my_argv[2] = starlesstif;
	if (args->customstride) my_argv[3] = args->stride;
	// *** Call starnet++ *** //
#ifdef _WIN32
	// add quotes around full path in case there are spaces
	memset(qtemptif, 0, sizeof(qtemptif));
	strcat(qtemptif, "\"");
	strcat(qtemptif, temptif);
	strcat(qtemptif, "\"");
	my_argv[1] = qtemptif;
	memset(qstarlesstif, 0, sizeof(qstarlesstif));
	strcat(qstarlesstif, "\"");
	strcat(qstarlesstif,  starlesstif);
	strcat(qstarlesstif,"\"");
	my_argv[2] = qstarlesstif;
	retval = exec_prog_win32(my_argv);
#else
	retval = exec_prog(my_argv);
#endif
	if (retval || forkerrors) {
		if (!retval && forkerrors)
		siril_log_color_message(_("Error: Starnet++ did not execute correctly...\n"), "red");
		goto CLEANUP;
	}

	// Read the starless stretched tiff. Successful return value of readtif() is nsamples
	retval = readtif(starlesstif, &workingfit, FALSE);
	if (retval < 1 || retval > 3) {
		siril_log_color_message(_("Error: unable to read starless image from TIFF...\n"), "red");
		goto CLEANUP;
	}
	/* we need to copy metadata as they have been removed with readtif     */
	copy_fits_metadata(current_fit, &workingfit);

	// Increase bit depth of starless image to 32 bit to improve precision
	// for subsequent processing. Only if !force_16bit otherwise there is an error on subtraction
	gboolean force_16bit = com.pref.force_16bit || (sequence_is_loaded() && args->seq && (args->seq->type == SEQ_SER || args->force_ser));
	if (!force_16bit) {
		const size_t ndata = workingfit.naxes[0] * workingfit.naxes[1] * workingfit.naxes[2];
		fit_replace_buffer(&workingfit, ushort_buffer_to_float(workingfit.data, ndata), DATA_FLOAT);
	}

	// Downscale again if needed
	if (args->upscale) {
		if (verbose)
			siril_log_message(_("Starnet++: 2x upscaling selected. Re-scaling starless image to original size...\n"));
		retval = cvResizeGaussian(&workingfit, orig_x, orig_y, OPENCV_AREA, FALSE);
		if (retval) {
			siril_log_color_message(_("Error: image resize failed...\n"), "red");
			goto CLEANUP;
		}
	}

	// If we are doing a pseudo-linear stretch we need to apply the inverse
	// stretch to the starless version and re-save the final result
	if (args->linear) {
		if (verbose)
			siril_log_message(_("Starnet++: linear mode. Applying inverse MTF stretch to starless image.\n"));
		apply_linked_pseudoinverse_mtf_to_fits(&workingfit, &workingfit, params, TRUE);
	}

	// Chdir back to the Siril working directory, we don't need to be in the starnet
	// directory any more
	retval = g_chdir(currentdir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
		goto CLEANUP1;
	}

	// Save workingfit as starless stretched image fits
	update_filter_information(&workingfit, "Starless", TRUE);
	if (single_image_is_loaded() && get_thread_run()) { // sequence worker will handle saving this in the sequence
		retval = savefits(starlessfit, &workingfit);
		if (retval) {
			siril_log_color_message(_("Error: unable to save starless image as FITS...\n"), "red");
			goto CLEANUP;
		}
		if (verbose)
			siril_log_color_message(_("Starnet++: starless image generated\n"), "green");
	}

	if (args->starmask) {
		// Subtract starless stretched from original stretched
		retval = imoper(&fit, &workingfit, OPER_SUB, !force_16bit);
		if (retval) {
			siril_log_color_message(_("Error: image subtraction failed...\n"), "red");
			goto CLEANUP;
		}
		update_filter_information(&fit, "StarMask", TRUE);

		// Save fit as starmask fits
		if (get_thread_run()) {
			if (single_image_is_loaded()) {
				retval = savefits(starmaskfit, &fit);
				if (retval) {
					siril_log_color_message(_("Error: unable to save starmask image as FITS...\n"), "red");
					goto CLEANUP;
				}
			} else {
				copyfits(&fit, args->starmask_fit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
			}
			if (verbose)
				siril_log_color_message(_("Starnet++: star mask generated\n"), "green");
		}
	}

	// Remove working files, they are no longer required
	retval = g_remove(starlesstif);
	if (retval) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
		// No goto here as even if it fails we want to try to remove the other TIFF
	}

	retval |= g_remove(temptif);
	if (retval) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
		goto CLEANUP;
	}

	// All done, now copy the working image back into gfit
	retval = copyfits(&workingfit, current_fit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
	if (retval) {
		siril_log_color_message(_("Error: image copy failed...\n"), "red");
		goto CLEANUP;
	}

	// Before CLEANUP so that this doesn't print on failure.
	if (verbose)
		siril_log_color_message(_("Starnet++: job completed.\n"), "green");

	if (single_image_is_loaded()) {
		free(com.uniq->filename);
		com.uniq->filename = strdup(starlessfit);
		if (args->follow_on) {
			struct remixargs *blendargs;
			blendargs = calloc(1, sizeof(struct remixargs));
			blendargs->fit1 = calloc(1, sizeof(fits));
			blendargs->fit2 = calloc(1, sizeof(fits));
			copyfits(&workingfit, blendargs->fit1, (CP_ALLOC | CP_COPYA |CP_FORMAT), -1);
			copyfits(&fit, blendargs->fit2, (CP_ALLOC | CP_COPYA |CP_FORMAT), -1);
			siril_add_idle(end_and_call_remixer, blendargs);
		}
	}

	CLEANUP:
	retval = g_chdir(currentdir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
	}
	CLEANUP1:
	if (args->starmask)
		clearfits(&fit);
	CLEANUP2:
	clearfits(&workingfit);
	CLEANUP3:
	if (com.child_is_running)
		com.child_is_running = FALSE;
	if (verbose)
		set_progress_bar_data("Ready.", PROGRESS_RESET);
	g_free(currentdir);
	g_free(starlesstif);
	g_free(starlessfit);
	g_free(starmaskfit);
	g_free(starlessnoext);
	g_free(starmasknoext);
	g_free(imagenoext);
	g_free(imagenoextorig);
	g_free(temptif);
	gettimeofday(&t_end, NULL);
	if (verbose)
		show_time(t_start, t_end);
	if (single_image_is_loaded()) {
		if (!args->follow_on)
			notify_gfit_modified();
		siril_add_idle(end_generic, args);
	}
	return GINT_TO_POINTER(retval);
}

static int starnet_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
// Starnet cannot run in parallel as it fully utilizes the GPU. This function therefore
// returns 1 and all images will be processed in series.
	return 1;
}

static int starnet_finalize_hook(struct generic_seq_args *args) {
	struct starnet_data *starnet_args = (struct starnet_data *) args->user;
	args->new_ser = starnet_args->new_ser_starless;
	args->new_fitseq = starnet_args->new_fitseq_starless;
	int retval = seq_finalize_hook(args);
	starnet_args->new_ser_starless = NULL;
	starnet_args->new_fitseq_starless = NULL;

	args->new_ser = starnet_args->new_ser_starmask;
	args->new_fitseq = starnet_args->new_fitseq_starmask;
	retval = seq_finalize_hook(args) || retval;
	starnet_args->new_ser_starmask = NULL;
	starnet_args->new_fitseq_starmask = NULL;
	seqwriter_set_number_of_outputs(1);
	return retval;
}

static int starnet_save_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	starnet_data *seqdata = (starnet_data *) args->user;
	if (!get_thread_run()) {
		return 1;
	}
	int retval1, retval2 = 0;
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		retval1 = ser_write_frame_from_fit(seqdata->new_ser_starless, seqdata->starnet_fit, out_index);
		if (seqdata->starmask)
			retval2 = ser_write_frame_from_fit(seqdata->new_ser_starmask, seqdata->starmask_fit, out_index);
		// the two fits are freed by the writing thread
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		retval1 = fitseq_write_image(seqdata->new_fitseq_starless, seqdata->starnet_fit, out_index);
		if (seqdata->starmask)
			retval2 = fitseq_write_image(seqdata->new_fitseq_starmask, seqdata->starmask_fit, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			g_free(fit);
			fit = NULL;
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, "starless_", in_index);
		if (fit->type == DATA_USHORT) {
			retval1 = save1fits16(dest, seqdata->starnet_fit, RLAYER);
		} else {
			retval1 = save1fits32(dest, seqdata->starnet_fit, RLAYER);
		}
		g_free(dest);
		clearfits(seqdata->starnet_fit);
		if (seqdata->starmask) {
			dest = fit_sequence_get_image_filename_prefixed(args->seq, "starmask_", in_index);
			if (fit->type == DATA_USHORT) {
				retval2 = save1fits16(dest, seqdata->starmask_fit, RLAYER);
			} else {
				retval2 = save1fits32(dest, seqdata->starmask_fit, RLAYER);
			}
			g_free(dest);
			clearfits(seqdata->starmask_fit);
			seqdata->starmask_fit = NULL;
		}
	}
	return retval1 || retval2;
}

int starnet_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 0;
	starnet_data *seqdata = (starnet_data *) args->user;
	seqdata->force_ser = args->force_ser_output;
	seqdata->starnet_fit = fit;
	seqdata->starmask_fit = calloc(1, sizeof(fits));
	seqdata->imgnumber = o;

	do_starnet(seqdata);
	printf("1");
	return ret;
}

static int starnet_basic_prepare_hook(struct generic_seq_args *args) {
	int retval = seq_prepare_hook(args);
	if (!retval && args->new_ser) {
		retval = ser_reset_to_monochrome(args->new_ser);
	}
	return retval;
}

static int starnet_prepare_hook(struct generic_seq_args *args) {
	struct starnet_data *starnet_args = (struct starnet_data *) args->user;
	// we call the generic prepare twice with different prefixes
	args->new_seq_prefix = "starless_";
	if (starnet_basic_prepare_hook(args))
		return 1;
	// but we copy the result between each call
	starnet_args->new_ser_starless = args->new_ser;
	starnet_args->new_fitseq_starless = args->new_fitseq;

	args->new_seq_prefix = "starmask_";
	if (starnet_basic_prepare_hook(args))
		return 1;
	starnet_args->new_ser_starmask = args->new_ser;
	starnet_args->new_fitseq_starmask = args->new_fitseq;

	args->new_seq_prefix = "starless_"; // Set the prefix for the sequence we want
										// loaded afterwards
	args->new_ser = NULL;
	args->new_fitseq = NULL;

	seqwriter_set_number_of_outputs(2);
	return 0;
}
void apply_starnet_to_sequence(struct starnet_data *seqdata) {
	seqdata->starnet_fit = NULL;
	struct generic_seq_args *seqargs = create_default_seqargs(seqdata->seq);
	seqargs->seq = seqdata->seq;
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seqdata->seq->selnum;
	seqargs->compute_mem_limits_hook = starnet_compute_mem_limits;
	seqargs->finalize_hook = starnet_finalize_hook;
	seqargs->save_hook = starnet_save_hook;
	seqargs->image_hook = starnet_image_hook;
	seqargs->prepare_hook = starnet_prepare_hook;
	seqargs->description = _("Starnet++");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	seqdata->seqEntry = g_strdup_printf("starless_");
	seqargs->new_seq_prefix = seqdata->seqEntry;
	seqargs->load_new_sequence = TRUE;
	seqargs->user = seqdata;
	const char *ptr = strrchr(seqdata->seq->seqname, G_DIR_SEPARATOR);
	if (ptr)
		seqdata->seqname = g_strdup_printf("%s%s%s", seqdata->seqEntry, ptr + 1, com.pref.ext);
	else seqdata->seqname = g_strdup_printf("%s%s%s", seqargs->new_seq_prefix, seqargs->seq->seqname, com.pref.ext);

	start_in_new_thread(generic_sequence_worker, seqargs);
}
#endif

