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
#include "algos/colors.h"
#include "algos/extraction.h"
#include "algos/geometry.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/mtf.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
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

// Check maximum path length - OSes except for Windows
#ifndef _WIN32
long get_pathmax(void)
{
	long pathmax = -1;

	errno = 0;
	pathmax = pathconf("/", _PC_PATH_MAX);
	if (-1 == pathmax) {
		if (0 == errno) {
#define PATHMAX_INFINITE_GUESS 4096
			pathmax = PATHMAX_INFINITE_GUESS;
		} else {
			fprintf (stderr, "pathconf() FAILED, %d, %s\n", errno, strerror(errno));
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

static int exec_prog(const char **argv)
{
	pid_t my_pid;
	int status;
	forkerrors = 0;
	int pipe_fds[2];
	int n;
	char buf[0x100] = {0};

	if (pipe (pipe_fds)) {
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
			if ((n = read(pipe_fds[0], buf, 0x100)) >= 0) {
				buf[n] = 0;
				char *m = strstr(buf, "pid:");
				if (m != 0) {
					m += 4;
					child_pid = atoi(m);
					com.child_is_running = TRUE;
					com.childpid = child_pid;
				} else {
					double value = g_ascii_strtod(buf, NULL);
					if (value != 0.0 && value == value) { //
						set_progress_bar_data("Running Starnet++", (value / 100));
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
 static int exec_prog_win32(const char **argv) {
	forkerrors = 0;
	int pipe_fds[2];
	int n;
	char buf[256] = {0};
	HANDLE hProcess;
	int nExitCode = STILL_ACTIVE;
	int fdStdOut;

	if (_pipe(pipe_fds, 0x100, O_TEXT)) {
		perror("pipe creation error");
		forkerrors = 1;
		return 0;
	}
	//from https://docs.microsoft.com/fr-fr/cpp/c-runtime-library/reference/pipe?view=msvc-170
	fdStdOut = _dup(_fileno(stdout));
	if(_dup2(pipe_fds[1], _fileno(stdout)) != 0) return -1;
	_close(pipe_fds[1]);
	hProcess = (HANDLE)_spawnve(P_NOWAIT, argv[0], argv, NULL);
	if(_dup2(fdStdOut, _fileno(stdout)) != 0) return -1;
	_close(fdStdOut);
	if (hProcess) {
		com.child_is_running = TRUE;
		com.childhandle = (void*) hProcess;
		while (nExitCode == STILL_ACTIVE) {
			usleep(100000);	// Wait for starnet++ to finish before attempting to process the output
			if ((n = _read(pipe_fds[0], buf, 256)) >= 0) {
				buf[n] = '\0';
				double value = g_ascii_strtod(buf, NULL);
				if (value != 0.0 && value == value) { //
					set_progress_bar_data("Running Starnet++", (value / 100));
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
	if ((gfit.naxes[2] == 3) && (g_file_test(rgbpath, G_FILE_TEST_IS_EXECUTABLE))) {
		g_free(rgbpath);
		g_free(monopath);
		return TRUE;
	} else
		g_free(rgbpath);
	if ((gfit.naxes[2] == 1 ) && (g_file_test(monopath, G_FILE_TEST_IS_EXECUTABLE))) {
		g_free(monopath);
		return TRUE;
	}
	else {
		g_free(monopath);
		return FALSE;
	}
}

/* Starnet++v2 star removal routine */

gpointer do_starnet(gpointer p) {
	int retval;
	int orig_x = gfit.rx, orig_y = gfit.ry;
	fits workingfit, fit;
	starnet_data *args = (starnet_data *) p;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);

	// Only allocate as much space for filenames as required - we determine the max pathlength
#ifndef _WIN32
	long pathmax = get_pathmax();
#else
	long pathmax = MAX_PATH;	// On Windows use of MAX_PATH is fine as it is not
								// a configurable item
#endif
	gchar *currentdir;
	gchar starnetcommand[16] = "starnet++";
	gchar temptif[pathmax];
	gchar starlesstif[pathmax];
	gchar starlessfit[pathmax];
	gchar starmaskfit[pathmax];
	char *imagenoext;
	gchar starnetsuffix[10] = "_starnet";
	gchar starlesssuffix[10] = "_starless";
	gchar starmasksuffix[10] = "_starmask";
	// Initialise the filename strings as empty strings
	memset(temptif, 0, sizeof(temptif));
	memset(starlesstif, 0, sizeof(starlesstif));
	memset(starlessfit, 0, sizeof(starlessfit));
	memset(starmaskfit, 0, sizeof(starmaskfit));
	memset(starnetcommand, 0, sizeof(starnetcommand));
	// Set up paths and filenames
	imagenoext = g_path_get_basename(com.uniq->filename);
	fprintf(stdout,"%s\n",imagenoext);
	imagenoext = g_build_filename(com.wd, imagenoext, NULL);
	fprintf(stdout,"%s\n",imagenoext);
	imagenoext = remove_ext_from_filename(imagenoext);
	fprintf(stdout,"%s\n",imagenoext);
	strncat(temptif,imagenoext,sizeof(temptif) - strlen(imagenoext));
	strncat(temptif,starnetsuffix, 10);
	strncat(temptif,".tif", 5);
	strncat(starlesstif,imagenoext,sizeof(starlesstif) - strlen(imagenoext));
	strncat(starlesstif,starlesssuffix, 10);
	strncat(starlesstif,".tif",5);
	strncat(starlessfit,imagenoext,sizeof(starlessfit) - strlen(imagenoext));
	strncat(starlessfit,starlesssuffix, 10);
	strncat(starlessfit,com.pref.ext,5);
	strncat(starmaskfit,imagenoext,sizeof(starmaskfit) - strlen(imagenoext));
	strncat(starmaskfit,starmasksuffix, 10);
	strncat(starmaskfit,com.pref.ext,5);

	// ok, let's start
	set_progress_bar_data(_("Starting Starnet++"), PROGRESS_NONE);

	// Store current working directory
	currentdir = g_get_current_dir();

	siril_log_color_message(_("Starnet++: running. Please wait...\n"), "green");
	if (args->customstride) siril_log_message(_("Starnet++: stride = %s...\n"), args->stride);
	if (!args->starmask) siril_log_message(_("Starnet++: -nostarmask invoked, star mask will not be generated...\n"));

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
	retval = copyfits(&gfit, &workingfit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
	if (retval) {
		siril_log_color_message(_("Error: image copy failed...\n"), "red");
		goto CLEANUP3;
	}

	// If workingfit is type DATA_FLOAT, check maximum value of gfit, if > 1.0 renormalize
	// in order to avoid cropping brightness peaks on saving to 16 bit TIFF
	if (gfit.type == DATA_FLOAT) {
		int nplane = gfit.naxes[2];
		double max = 0.0;
		for (int layer = 0; layer < nplane; layer++) {
			imstats* stat = statistics(NULL, -1, &gfit, layer, &com.selection, STATS_MAIN, MULTI_THREADED);
			if (stat->max> max)
				max = stat->max;
		}
		if (max > 1.0) {
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
	}

	// If we need to pre-stretch a linear image, we use the auto histogram stretch with
	// default shadow clipping and target background. This does marginally clip the
	// shadows but generally by less than 0.001% of pixels. The result of starnet using
	// this stretch is much better than either asinh or GHT stretches.
	struct mtf_params params;
	find_linked_midtones_balance_default(&workingfit, &params);
	if (args->linear) {
		siril_log_message(_("Starnet++: linear mode. Applying Midtone Transfer Function (MTF) pre-stretch to image.\n"));
		apply_linked_mtf_to_fits(&workingfit, &workingfit, params);
	}

	// Upscale if needed
	if (args->upscale) {
		siril_log_message(_("Starnet++: 2x upscaling selected. Upscaling image...\n"));
		retval = cvResizeGaussian(&workingfit, round_to_int(2*orig_x), round_to_int(2*orig_y), OPENCV_AREA);
		if (retval) {
			siril_log_color_message(_("Error: image resize failed...\n"), "red");
			goto CLEANUP;
		}
	}

	// Save current stretched image as working 16-bit TIFF (post initial stretch if the image was linear)
	if (savetif(temptif, &workingfit, 16)) {
		siril_log_color_message(_("Error: unable to save working TIFF of original image...\n"), "red");
		goto CLEANUP;
	}

	// Check for starnet executables (pre-v2.0.2 or v2.0.2+)
	if (g_file_test(STARNET_BIN, G_FILE_TEST_IS_EXECUTABLE)) {
		snprintf(starnetcommand,15, STARNET_BIN);
	} else if ((gfit.naxes[2] == 3) && (g_file_test(STARNET_RGB, G_FILE_TEST_IS_EXECUTABLE))) {
		snprintf(starnetcommand,15, STARNET_RGB);
	} else if ((gfit.naxes[2] == 1 ) && (g_file_test(STARNET_MONO, G_FILE_TEST_IS_EXECUTABLE))) {
		snprintf(starnetcommand, 15, STARNET_MONO);
	}
	else {
		retval = 1;
		siril_log_color_message(_("No valid executable found in the Starnet++ directory\n"), "red");
		goto CLEANUP;
	}
	my_argv[0] = starnetcommand;
	my_argv[1] = temptif;
	my_argv[2] = starlesstif;
	if (args->customstride) my_argv[3] = args->stride;
	// *** Call starnet++ *** //
#ifdef _WIN32
	retval = exec_prog_win32(my_argv);
#else
	retval = exec_prog(my_argv);
#endif
	if (retval || forkerrors) {
		if (!retval && forkerrors)
			retval = forkerrors;
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
	copy_fits_metadata(&gfit, &workingfit);

	// Increase bit depth of starless image to 32 bit to improve precision
	// for subsequent processing
	const size_t ndata = workingfit.naxes[0] * workingfit.naxes[1] * workingfit.naxes[2];
	fit_replace_buffer(&workingfit, ushort_buffer_to_float(workingfit.data, ndata), DATA_FLOAT);

	// Downscale again if needed
	if (args->upscale) {
		siril_log_message(_("Starnet++: 2x upscaling selected. Re-scaling starless image to original size...\n"));
		retval = cvResizeGaussian(&workingfit, orig_x, orig_y, OPENCV_AREA);
		if (retval) {
			siril_log_color_message(_("Error: image resize failed...\n"), "red");
			goto CLEANUP;
		}
	}

	// If we are doing a pseudo-linear stretch we need to apply the inverse
	// stretch to the starless version and re-save the final result
	if (args->linear) {
		siril_log_message(_("Starnet++: linear mode. Applying inverse MTF stretch to starless image.\n"));
		apply_linked_pseudoinverse_mtf_to_fits(&workingfit, &workingfit, params);
	}

	// Chdir back to the Siril working directory, we don't need to be in the starnet
	// directory any more
	retval = g_chdir(currentdir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
		goto CLEANUP1;
	}

	// Save starless stretched image as fits
	update_filter_information(&workingfit, "Starless", TRUE);
	retval = savefits(starlessfit, &workingfit);
	if (retval) {
		siril_log_color_message(_("Error: unable to save starless image as FITS...\n"), "red");
		goto CLEANUP;
	}
	siril_log_color_message(_("Starnet++: starless image generated\n"), "green");

	if (args->starmask) {
		// Subtract starless stretched from original stretched
		retval = imoper(&fit, &workingfit, OPER_SUB, !com.pref.force_16bit);
		if (retval) {
			siril_log_color_message(_("Error: image subtraction failed...\n"), "red");
			goto CLEANUP;
		}
		update_filter_information(&fit, "StarMask", TRUE);

		// Save starmask as fits
		retval = savefits(starmaskfit, &fit);
		if (retval) {
			siril_log_color_message(_("Error: unable to save starmask image as FITS...\n"), "red");
			goto CLEANUP;
		}
		siril_log_color_message(_("Starnet++: star mask generated\n"), "green");
	}

	// Remove working files, they are no longer required
	retval = remove(starlesstif);
	if (retval) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
		// No goto here as even if it fails we want to try to remove the other TIFF
	}

	retval |= remove(temptif);
	if (retval) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
		goto CLEANUP;
	}

	// All done, now copy the working image back into gfit
	retval = copyfits(&workingfit, &gfit, (CP_ALLOC | CP_INIT | CP_FORMAT | CP_COPYA), 0);
	if (retval) {
		siril_log_color_message(_("Error: image copy failed...\n"), "red");
		goto CLEANUP;
	}

	// Before CLEANUP so that this doesn't print on failure.
	siril_log_color_message(_("Starnet++: job completed.\n"), "green");

	free(com.uniq->filename);
	com.uniq->filename = strdup(_(starlessfit));


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
	set_progress_bar_data("Ready.", PROGRESS_RESET);
	free(imagenoext);
	free(currentdir);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	notify_gfit_modified();

	return GINT_TO_POINTER(retval);
}
