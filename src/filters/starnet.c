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

// Comment out this #define before public release
//#define STARNET_DEBUG

#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <gio/gwin32inputstream.h>
#else
#include <sys/types.h> // for waitpid(2)
#include <sys/wait.h> // for waitpid(2)
#include <gio/gunixinputstream.h>
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "starnet.h"

#ifdef HAVE_LIBTIFF
fits *current_fit = NULL;
gboolean verbose = TRUE;

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	g_spawn_close_pid(pid);
}

static int exec_prog_starnet(char **argv, starnet_version version) {
	gint child_stdout;
	GPid child_pid;
	g_autoptr(GError) error = NULL;
	int retval = -1;

	// g_spawn handles wchar so not need to convert
	g_spawn_async_with_pipes(NULL, argv, NULL,
			G_SPAWN_SEARCH_PATH |
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_STDERR_TO_DEV_NULL | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, &child_pid, NULL, &child_stdout,
			NULL, &error);

	if (error != NULL) {
		siril_log_color_message(_("Spawning starnet failed: %s\n"), "red", error->message);
		return retval;
	}
	g_child_watch_add(child_pid, child_watch_cb, NULL);
	com.child_is_running = EXT_STARNET;
#ifdef _WIN32
	com.childhandle = child_pid;		// For Windows, handle of a child process
#else
	com.child_pid = child_pid;			// For other OSes, PID of a child process
#endif

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stdout), FALSE);
#else
	stream = g_unix_input_stream_new(child_stdout, FALSE);
#endif

	gboolean doprint = TRUE;
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
#ifdef STARNET_DEBUG
	gchar *lastbuffer = NULL;
#endif
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
		if (doprint && verbose)
			siril_log_message("StarNet: %s\n", buffer);
		// TODO: the Mac M1 version definitely contains the substring "backend" in
		// the line before the progress percentages start. Need to check this is
		// still true when the Torch version comes out for other targets.
		if (g_str_has_prefix(buffer, "Total number of tiles") || g_strrstr(buffer, "backend")) {
			// need to set EOL to CR otherwise can't read the progress percentage
			// which just sends xx% progress/r to avoid flushing
			g_data_input_stream_set_newline_type(data_input, G_DATA_STREAM_NEWLINE_TYPE_CR);
			doprint = FALSE;
		}
		gchar *arg = buffer;
		if (g_str_has_prefix(buffer, "Working:"))
			arg += 9;
		double value = g_ascii_strtod(arg, NULL);
		if (value != 0.0 && value == value && verbose) {
			set_progress_bar_data(_("Running StarNet"), value / 100.0);
		}
		if (g_str_has_prefix(buffer, "100% finished") || g_strrstr(buffer, "Writing mask")) {
			retval = 0;
#ifdef STARNET_DEBUG
			if (verbose)
				siril_log_message(_("StarNet caught buffer: %s : exit successful\n"), buffer);
#endif
		}
#ifdef STARNET_DEBUG
		lastbuffer = g_strdup(buffer);
#endif
		g_free(buffer);
	}
#ifdef STARNET_DEBUG
	if (retval)
		siril_log_message(_("Starnet exit not caught, last buffer: %s\n"), lastbuffer);
	g_free(lastbuffer);
#endif
	g_object_unref(data_input);
	g_object_unref(stream);
	if (!g_close(child_stdout, &error))
		siril_debug_print("%s\n", error->message);
	return retval;
}

starnet_version starnet_executablecheck(gchar* executable) {
	char *test_argv[3] = {0};
	int retval = NIL;
	gint child_stdout;
	g_autoptr(GError) error = NULL;
	gchar *v1dir = NULL;
	if (!executable || executable[0] == '\0') {
		return NIL;
	}
	if (!g_file_test(executable, G_FILE_TEST_IS_EXECUTABLE)) {
		return NIL; // It's not executable so return NIL
	}
	// V1 tests
	if (g_str_has_suffix(executable, "rgb_starnet++") || g_str_has_suffix(executable, "rgb_starnet++.exe") || g_str_has_suffix(executable, "mono_starnet++") || g_str_has_suffix(executable, "mono_starnet++.exe")) {
		v1dir = g_path_get_dirname(executable);
		if (g_str_has_suffix(executable, "rgb_starnet++")) {
			retval = V1RGB;
			goto END;
		} else if (g_str_has_suffix(executable, "rgb_starnet++.exe")) {
			retval = V1RGB;
			goto END;
		} else if (g_str_has_suffix(executable, "mono_starnet++")) {
			retval = V1MONO;
			goto END;
		} else if (g_str_has_suffix(executable, "mono_starnet++.exe")) {
			retval = V1MONO;
			goto END;
		}
	}
	// Not V1: execute file and test output to determine version
	// Get the starnet installation directory and chdir to it
	gchar *dir = g_path_get_dirname(executable);
	gchar *currentdir = g_get_current_dir();
	int retval2 = g_chdir(dir);
	if (retval2) {
		retval = NIL;
		g_free(dir);
		g_free(currentdir);
		return NIL;
	}
	g_free(dir);

	int nb = 0;
	test_argv[nb++] = executable;
	gchar *versionarg = g_strdup("--version");
	test_argv[nb++] = versionarg;
	// g_spawn handles wchar so not need to convert
	g_spawn_async_with_pipes(NULL, test_argv, NULL,
			G_SPAWN_SEARCH_PATH |
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_STDERR_TO_DEV_NULL,
			NULL, NULL, NULL, NULL, &child_stdout,
			NULL, &error);

	if (error != NULL) {
		siril_log_color_message(_("Spawning starnet failed during version check: %s\n"), "red", error->message);
		g_free(versionarg);
		return NIL;
	}

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stdout), FALSE);
#else
	stream = g_unix_input_stream_new(child_stdout, FALSE);
#endif
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
	gboolean done = FALSE;
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL)) && !done) {
#ifdef STARNET_DEBUG
		siril_log_message("%s\n", buffer);
#endif
		if (g_strrstr(buffer, "StarNet++ v2")) {
#ifdef STARNET_DEBUG
			siril_log_message(_("V2 detected\n"));
#endif
			retval = V2;
			done = TRUE;
		} else if (g_strrstr(buffer, " version:")) {
#ifdef STARNET_DEBUG
			siril_log_message(_("V2-torch detected\n"));
#endif
			retval = TORCH; // We ignore the actual version number for now, it's enough that this substring is
							 // a unique identifier for the Torch-based StarNet versions.
			done = TRUE;
		}
		g_free(buffer);
	}
	g_object_unref(data_input);
	g_object_unref(stream);
	g_free(versionarg);
	if (!g_close(child_stdout, &error))
		siril_debug_print("%s\n", error->message);

	retval2 = g_chdir(currentdir);
	if (retval2) {
		retval2 = g_chdir(com.wd);
		if (retval2) {
			siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
		}
	}
	g_free(currentdir);

END:
	g_free(v1dir);

	return retval;
}

void free_starnet_args(starnet_data *args) {
	g_free(args->seqEntry);
	g_free(args->stride);
	free(args);
	siril_debug_print("starnet_args freed\n");
}

gboolean end_starnet(gpointer p) {
	starnet_data *args = (starnet_data *) p;
	free_starnet_args(args);
	return end_generic(NULL);
}

gboolean end_and_call_remixer(gpointer p)
{
	struct remixargs *blendargs = (remixargs *) p;
	toggle_remixer_window_visibility(CALL_FROM_STARNET, blendargs->fit1, blendargs->fit2);
	free(blendargs);
	return end_generic(NULL);
}

/* StarNet star removal routine */

gpointer do_starnet(gpointer p) {
	char *my_argv[64] = { 0 };
	starnet_version version = NIL;
	int retval = 0;
	int retval2 = 0;
	fits workingfit = { 0 }, fit = { 0 };
	starnet_data *args = (starnet_data *) p;
	verbose = (args->seq == NULL); // To suppress log messages during seq working
	args->follow_on = args->seq ? FALSE : args->follow_on;
	current_fit = args->starnet_fit;
	int orig_x = current_fit->rx, orig_y = current_fit->ry;
	struct remixargs *blendargs = NULL;
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
	gchar *starnetcommand = NULL;
	gchar *temptif = NULL;
	gchar *starlesstif = NULL;
	gchar *starmasktif = NULL;
	gchar *starlessfit = NULL;
	gchar *starmaskfit = NULL;
	gchar *imagenoext = NULL;
	gchar *starlessnoext = NULL;
	gchar *starmasknoext = NULL;
	gchar *imagenoextorig = NULL;
	gchar *torcharg_in = NULL;
	gchar *torcharg_out = NULL;
	gchar *torcharg_stride = NULL;
	gchar *torcharg_up = NULL;
	gchar *torcharg_weights = NULL;
	gchar *torcharg_mask = NULL; // We don't want the mask output, but the new starnet generates it regardless
			// so we may as well put it somewhere known so we can delete it so as not to leave behind clutter.
	gchar *temp = NULL;
	gchar starnetprefix[10] = "starnet_";
	gchar starlessprefix[10] = "starless_";
	gchar starmaskprefix[10] = "starmask_";

	// Check StarNet version
	version = starnet_executablecheck(com.pref.starnet_exe);

	// Check image size is sufficient
	unsigned int minsize = 512;
	if (args->upscale && version != TORCH)
		minsize = 256;
	if (min(current_fit->rx, current_fit->ry) < minsize) {
		args->too_small = TRUE;
		if (verbose) {
			siril_log_color_message(_("Warning: some versions of StarNet fail to process images with dimensions smaller than 512x512. This image may fail. Attempting to continue...\n"), "salmon");
		}
	}

	// If v1, check we have the right executable for the number of channels in current_fit
	if (version == NIL) {
		siril_log_color_message(_("No suitable StarNet executable found.\n"), "red");
		goto CLEANUP;
	} else if (version == V1RGB) {
		if (current_fit->naxes[2] == 3) {
			starnetcommand = g_strdup(com.pref.starnet_exe);
		} else {
			gchar *temp = g_path_get_dirname(com.pref.starnet_exe);
			gchar *winext = g_str_has_suffix(com.pref.starnet_exe, ".exe") ? g_strdup(".exe") : g_strdup("\0");
			starnetcommand = g_strdup_printf("%s/mono_starnet++%s", temp, winext);
			g_free(temp);
			g_free(winext);
			if (starnet_executablecheck(starnetcommand) != V1MONO) {
				siril_log_color_message(_("No suitable StarNet executable found for a mono image.\n"), "red");
				retval = 1;
				goto CLEANUP3;
			}
		}
	} else if (version == V1MONO) {
		if (current_fit->naxes[2] == 1) {
			starnetcommand = g_strdup(com.pref.starnet_exe);
		} else {
			gchar *temp = g_path_get_dirname(com.pref.starnet_exe);
			gchar *winext = g_str_has_suffix(com.pref.starnet_exe, ".exe") ? g_strdup(".exe") : g_strdup("\0");
			starnetcommand = g_strdup_printf("%s/rgb_starnet++%s", temp, winext);
			g_free(temp);
			g_free(winext);
			if (starnet_executablecheck(starnetcommand) != V1RGB) {
				siril_log_color_message(_("No suitable StarNet executable found for a RGB image.\n"), "red");
				retval = 1;
				goto CLEANUP3;
			}
		}
	} else if (version == V2 || version & TORCH) {
		starnetcommand = g_strdup(com.pref.starnet_exe);
	}

	siril_debug_print("Using %s\n", starnetcommand);

	// Set up paths and filenames
	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		temp = g_path_get_basename(com.uniq->filename);
		imagenoextorig = g_strdup_printf("%s", temp);
		g_free(temp);
	} else if (args->seq) {
		imagenoextorig = g_strdup_printf("%s%.5d", args->seq->seqname, args->imgnumber + 1);
	} else {
		imagenoextorig = g_strdup_printf("image");
	}
	imagenoext = g_strdup(imagenoextorig);
	starlessnoext = g_strdup_printf("%s%s", starlessprefix, imagenoext);
	starmasknoext = g_strdup_printf("%s%s", starmaskprefix, imagenoext);
	imagenoext = g_strdup_printf("%s%s", starnetprefix, imagenoext);
	temp = g_build_filename(com.wd, imagenoext, NULL);
	g_free(imagenoext);
	imagenoext = remove_ext_from_filename(temp);
	g_free(temp);
	temp = NULL;
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
	starlessfit = g_strdup(temp);
	g_free(temp);
	temp = g_strdup_printf("%s.tif", starlesstif);
	g_free(starlesstif);
	starlesstif = g_strdup(temp);
	g_free(temp);

	starmasktif = g_build_filename(com.wd, starmasknoext, NULL);
	temp = remove_ext_from_filename(starmasktif);
	g_free(starmasktif);
	starmasktif = g_strdup(temp);
	g_free(temp);
	temp = g_strdup_printf("%s.tif", starmasktif);
	g_free(starmasktif);
	starmasktif = g_strdup(temp);
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

//	printf("%s\n%s\n%s\n%s\n", starlesstif, starlessfit, starmasktif, starmaskfit);
	// ok, let's start
	if (verbose)
		set_progress_bar_data(_("Starting StarNet"), PROGRESS_NONE);

	// Store current working directory
	currentdir = g_get_current_dir();

	if (verbose)
		siril_log_color_message(_("StarNet: running. Please wait...\n"), "green");
	if (args->customstride && verbose) siril_log_message(_("StarNet: stride = %s...\n"), args->stride);
	if (!args->starmask && verbose) siril_log_message(_("StarNet: -nostarmask invoked, star mask will not be generated...\n"));

	// Check starnet directory is not NULL - can happen first time the new preference file is loaded
	if (!com.pref.starnet_exe) {
		retval = 1;
		siril_log_color_message(_("Incorrect permissions on the StarNet executable: %s\nEnsure it is correctly set in Preferences / Miscellaneous.\n"), "red", com.pref.starnet_exe);
		goto CLEANUP3;
	}

	// Check starnet directory is defined
	gchar* dir = g_path_get_dirname(com.pref.starnet_exe);
	retval = g_access(dir, R_OK);
	if (retval) {
		siril_log_color_message(_("Incorrect permissions on the StarNet directory: %s\nEnsure it is correctly set in Preferences / Miscellaneous.\n"), "red", dir);
		g_free(dir);
		goto CLEANUP3;
		// Dijkstra might be spinning in his grave but one of the few legitimate uses
		// of gotos is this type of error cleanup and return. Much more readable than a
		// mess of indentation.
	}

	// Change to starnet directory
	retval = g_chdir(dir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to StarNet directory.\nEnsure it is set in Preferences / Miscellaneous...\n"), "red");
		g_free(dir);
		goto CLEANUP3;
	}
	g_free(dir);

	// Create a working copy of the image.
	// copyfits(from, to, oper, layer)
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
			siril_log_message(_("StarNet: Pixel values exceed 1.0. Rescaling to avoid clipping peaks.\n"));
		}
		if (max > 1.0)
			soper(&workingfit, max, OPER_DIV, FALSE);
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
	retval = find_linked_midtones_balance_default(&workingfit, &params);
	if (retval && args->linear) {
		siril_log_color_message(_("Error: unable to find the MTF stretch factors...\n"), "red");
		goto CLEANUP;
	}
	if (args->linear) {
		if (verbose)
			siril_log_message(_("StarNet: linear mode. Applying Midtone Transfer Function (MTF) pre-stretch to image.\n"));
		apply_linked_mtf_to_fits(&workingfit, &workingfit, params, TRUE);
	}

	// Upscale if needed
	if (args->upscale && !(version & TORCH)) {
		if (verbose)
			siril_log_message(_("StarNet: 2x upscaling selected. Upscaling image...\n"));
		retval = cvResizeGaussian(&workingfit, round_to_int(2*orig_x), round_to_int(2*orig_y), OPENCV_AREA, FALSE);
		if (retval) {
			siril_log_color_message(_("Error: image resize failed...\n"), "red");
			goto CLEANUP;
		}
	}

	// Save current stretched image as working 16-bit TIFF (post initial stretch if the image was linear)
	retval = savetif(temptif, &workingfit, 16, NULL, com.pref.copyright, FALSE, TRUE, FALSE);
	if (retval) {
		siril_log_color_message(_("Error: unable to save working StarNet input file...\n"), "red");
		goto CLEANUP;
	}

	// Process StarNet arguments
	int nb = 0;
	my_argv[nb++] = starnetcommand;
	if (version & TORCH) {
		torcharg_in = g_strdup_printf("-i %s", temptif);
		my_argv[nb++] = torcharg_in;

		torcharg_out = g_strdup_printf("-o %s", starlesstif);
		my_argv[nb++] = torcharg_out;
	} else {
		my_argv[nb++] = temptif;
		my_argv[nb++] = starlesstif;
	}
	if (args->customstride) {
		if (version & TORCH) {
			torcharg_stride = g_strdup_printf("-s %s", args->stride);
			my_argv[nb++] = torcharg_stride;
		} else {
			my_argv[nb++] = args->stride;
		}
	}
	if (version & TORCH) {
		if (args->upscale) {
			torcharg_up = g_strdup("-u");
			my_argv[nb++] = torcharg_up;
		}
		if (com.pref.starnet_weights && com.pref.starnet_weights[0] != '\0') {
			if (g_access(com.pref.starnet_weights, R_OK)) {
				siril_log_color_message(_("Error: cannot read the neural net weights file.\n"), "red");
				goto CLEANUP;
			}
			torcharg_weights = g_strdup_printf("-w %s", com.pref.starnet_weights);
			my_argv[nb++] = torcharg_weights;
		}
		torcharg_mask = g_strdup_printf("-m %s", starmasktif);
		my_argv[nb++] = torcharg_mask;
	}

	// *** Call starnet++ *** //
	retval = exec_prog_starnet(my_argv, version);
	g_free(starnetcommand);
	starnetcommand = NULL;
	if (retval) {
		siril_log_color_message(_("Error: StarNet did not execute correctly...\n"), "red");
		goto CLEANUP;
	}

	// Read the starless stretched tiff. Successful return value of readtif() is nsamples
	clearfits(&workingfit); // Clear it first to free the data
	retval = readtif(starlesstif, &workingfit, FALSE, FALSE);
	if (retval < 1 || retval > 3) {
		siril_log_color_message(_("Error: unable to read StarNet output file...\n"), "red");
		goto CLEANUP;
	}

	// Remove working TIFF files, they are no longer required
	retval = g_remove(starlesstif);
	retval |= (g_remove(starmasktif) && (version & TORCH));
	retval |= g_remove(temptif);
	if (retval) {
		siril_log_color_message(_("Error: unable to remove temporary working file...\n"), "red");
		siril_log_message(_("Attempting to continue. You will need to clean up the working files manually.\n"));
	}

	/* we need to copy metadata as they have been removed with readtif     */
	copy_fits_metadata(current_fit, &workingfit);
	copy_fits_metadata(current_fit, &fit);

	// Increase bit depth of starless image to 32 bit to improve precision
	// for subsequent processing. Only if !force_16bit otherwise there is an error on subtraction

	//force_16_bit needs to be generated carefully because of the stacking result corner case
	gboolean force_16bit = com.pref.force_16bit;
	if (sequence_is_loaded()) {
		if ((!(com.seq.current == RESULT_IMAGE || com.seq.current == UNRELATED_IMAGE)) && args->seq && (args->seq->type == SEQ_SER || args->force_ser)) {
			force_16bit = TRUE;
		}
	} else if (args->seq) {
		if (args->seq->type == SEQ_SER || args->force_ser) {
			force_16bit = TRUE;
		}
	}

	if (!force_16bit) {
		const size_t ndata = workingfit.naxes[0] * workingfit.naxes[1] * workingfit.naxes[2];
		fit_replace_buffer(&workingfit, ushort_buffer_to_float(workingfit.data, ndata), DATA_FLOAT);
	}

	// Downscale again if needed
	if (args->upscale && !(version & TORCH)) {
		if (verbose)
			siril_log_message(_("StarNet: 2x upscaling selected. Re-scaling starless image to original size...\n"));
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
			siril_log_message(_("StarNet: linear mode. Applying inverse MTF stretch to starless image.\n"));
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
	if ((!args->seq) && get_thread_run()) { // sequence worker will handle saving this in the sequence
		retval = savefits(starlessfit, &workingfit);
		if (retval) {
			siril_log_color_message(_("Error: unable to save starless image as FITS...\n"), "red");
			goto CLEANUP;
		}
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
			if ((!args->seq)) {
				retval = savefits(starmaskfit, &fit);
				if (retval) {
					siril_log_color_message(_("Error: unable to save starmask image as FITS...\n"), "red");
					goto CLEANUP;
				}
			} else {
				if (args->starmask_fit)
					clearfits(args->starmask_fit);
				retval = copyfits(&fit, args->starmask_fit, (CP_ALLOC | CP_FORMAT | CP_COPYA), 0);
				if (retval) {
					siril_log_color_message(_("Error: image copy failed...\n"), "red");
					goto CLEANUP;
				}
				copy_fits_metadata(&fit, args->starmask_fit);
			}
		}
	}

	// All done, now copy the working image back into gfit
	clearfits(current_fit);
	retval = copyfits(&workingfit, current_fit, (CP_ALLOC | CP_FORMAT | CP_COPYA), 0);
	if (retval) {
		siril_log_color_message(_("Error: image copy failed...\n"), "red");
		goto CLEANUP;
	}
	copy_fits_metadata(&workingfit, current_fit);
	// Before CLEANUP so that this doesn't print on failure.
	if (verbose)
		siril_log_color_message(_("StarNet: job completed.\n"), "green");

	if ((!args->seq)) {
		free(com.uniq->filename);
		com.uniq->filename = strdup(starlessfit);
		if (args->follow_on) {
			blendargs = calloc(1, sizeof(struct remixargs));
			blendargs->fit1 = calloc(1, sizeof(fits));
			blendargs->fit2 = calloc(1, sizeof(fits));
			retval = copyfits(&workingfit, blendargs->fit1, (CP_ALLOC | CP_COPYA |CP_FORMAT), -1);
			if (retval) {
				siril_log_color_message(_("Error: image copy failed...\n"), "red");
				goto CLEANUP;
			}
			copy_fits_metadata(&workingfit, blendargs->fit1);
			retval = copyfits(&fit, blendargs->fit2, (CP_ALLOC | CP_COPYA |CP_FORMAT), -1);
			if (retval) {
				siril_log_color_message(_("Error: image copy failed...\n"), "red");
				goto CLEANUP;
			}
			copy_fits_metadata(&fit, blendargs->fit2);
		}
	}

	CLEANUP:
	if (currentdir) {
		retval2 = g_chdir(currentdir);
		if (retval2) {
			siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
			retval = retval2;
		}
	}
	CLEANUP1:
	if (args->starmask) {
		clearfits(&fit);
	}
	CLEANUP2:
	clearfits(&workingfit);
	CLEANUP3:
	if (com.child_is_running == EXT_STARNET)
		com.child_is_running = EXT_NONE;
	if (verbose)
		set_progress_bar_data("Ready.", PROGRESS_RESET);
	g_free(currentdir);
	g_free(starlesstif); // filename
	g_free(starlessfit); // filename
	g_free(starmaskfit); // filename
	g_free(starlessnoext); // part filename
	g_free(starmasknoext); // part filename
	g_free(imagenoext); // part filename
	g_free(imagenoextorig); // part filename
	g_free(temptif); // part filename
	g_free(torcharg_in);
	g_free(torcharg_out);
	g_free(torcharg_stride);
	g_free(torcharg_weights);
	g_free(torcharg_mask);
	g_free(torcharg_up);
	gettimeofday(&t_end, NULL);
	if (verbose)
		show_time(t_start, t_end);
	if ((!args->seq)) {
		if (args->follow_on) {
			free_starnet_args(args);
			if (!(retval)) {
				siril_add_idle(end_and_call_remixer, blendargs);
				return GINT_TO_POINTER(retval);
			} else {
				if (blendargs && blendargs->fit1) {
					clearfits(blendargs->fit1);
					free(blendargs->fit1);
				}
				if (blendargs && blendargs->fit2) {
					clearfits(blendargs->fit2);
					free(blendargs->fit2);
				}
				if (blendargs)
					free(blendargs);
				return GINT_TO_POINTER(retval);
			}
		} else {
			notify_gfit_modified();
			siril_add_idle(end_starnet, args);
			return GINT_TO_POINTER(retval);
		}
	}
	return GINT_TO_POINTER(retval);
}

static int starnet_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
// StarNet cannot run in parallel as it fully utilizes the GPU. This function therefore
// returns 1 and all images will be processed in series.
	return 1;
}

static int starnet_finalize_hook(struct generic_seq_args *args) {
	struct starnet_data *starnet_args = (struct starnet_data *) args->user;
	if (starnet_args->too_small)
		siril_log_color_message(_("Warning: some images in the sequence were smaller than 512x512. Some versions of StarNet may fail to process these images.\n"), "salmon");
	args->new_ser = starnet_args->new_ser_starless;
	args->new_fitseq = starnet_args->new_fitseq_starless;
	int retval = seq_finalize_hook(args);
	starnet_args->new_ser_starless = NULL;
	starnet_args->new_fitseq_starless = NULL;

	if (starnet_args->starmask) {
		args->new_ser = starnet_args->new_ser_starmask;
		args->new_fitseq = starnet_args->new_fitseq_starmask;
		retval |= seq_finalize_hook(args);
		starnet_args->new_ser_starmask = NULL;
		starnet_args->new_fitseq_starmask = NULL;
		seqwriter_set_number_of_outputs(1);
	}
	free_starnet_args(starnet_args);
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
		if (seqdata->starmask) {
			retval2 = ser_write_frame_from_fit(seqdata->new_ser_starmask, seqdata->starmask_fit, out_index);
		}
		// the two fits are freed by the writing thread
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		retval1 = fitseq_write_image(seqdata->new_fitseq_starless, seqdata->starnet_fit, out_index);
		if (seqdata->starmask) {
			retval2 = fitseq_write_image(seqdata->new_fitseq_starmask, seqdata->starmask_fit, out_index);
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, "starless_", in_index);
		retval1 = savefits(dest, seqdata->starnet_fit);
		free(dest);
		clearfits(seqdata->starnet_fit);
		seqdata->starnet_fit = NULL;
		if (seqdata->starmask) {
			dest = fit_sequence_get_image_filename_prefixed(args->seq, "starmask_", in_index);
			retval2 = savefits(dest, seqdata->starmask_fit);
			free(dest);
			clearfits(seqdata->starmask_fit);
			free(seqdata->starmask_fit);
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
	if (seqdata->starmask)
		seqdata->starmask_fit = calloc(1, sizeof(fits));
	seqdata->imgnumber = o;
	siril_log_color_message(_("Starnet: Processing image %d\n"), "green", o + 1);
	do_starnet(seqdata);
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
	args->new_seq_prefix = strdup("starless_");
	if (starnet_basic_prepare_hook(args))
		return 1;
	// but we copy the result between each call
	starnet_args->new_ser_starless = args->new_ser;
	starnet_args->new_fitseq_starless = args->new_fitseq;
	free(args->new_seq_prefix);

	if (starnet_args->starmask) {
		args->new_seq_prefix = strdup("starmask_");
		if (starnet_basic_prepare_hook(args))
			return 1;
		starnet_args->new_ser_starmask = args->new_ser;
		starnet_args->new_fitseq_starmask = args->new_fitseq;
		free(args->new_seq_prefix);
	}
	// Set the prefix for the sequence we want loaded afterwards
	args->new_seq_prefix = strdup("starless_");
	args->new_ser = NULL;
	args->new_fitseq = NULL;

	if (starnet_args->starmask)
		seqwriter_set_number_of_outputs(2);
	else
		seqwriter_set_number_of_outputs(1);

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
	seqargs->description = _("StarNet");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	seqargs->new_seq_prefix = seqdata->seqEntry;
	seqargs->load_new_sequence = TRUE;
	seqargs->user = seqdata;
	const char *ptr = strrchr(seqdata->seq->seqname, G_DIR_SEPARATOR);
	if (ptr)
		seqdata->seqname = g_strdup_printf("%s%s%s", seqdata->seqEntry, ptr + 1, com.pref.ext);
	else seqdata->seqname = g_strdup_printf("%s%s%s", seqargs->new_seq_prefix, seqargs->seq->seqname, com.pref.ext);
	set_progress_bar_data(_("StarNet: Processing..."), 0.);
	start_in_new_thread(generic_sequence_worker, seqargs);
}

#endif
