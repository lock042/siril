/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "core/icc_profile.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_update.h"
#include "core/siril_log.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "filters/graxpert.h"

//#define GRAXPERT_DEBUG

static gboolean verbose = TRUE;

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	siril_debug_print("GraXpert is being closed\n");
	g_spawn_close_pid(pid);
}

static int exec_prog_graxpert(char **argv, gboolean is_gui) {
	const gchar *progress_key = "Progress: ";
	gint child_stderr;
	GPid child_pid;
	g_autoptr(GError) error = NULL;
	int retval = -1;

	int index = 0;
	while (argv[index]) {
		fprintf(stdout, "%s ", argv[index++]);
	}
	fprintf(stdout, "\n");
	// g_spawn handles wchar so not need to convert
	set_progress_bar_data(_("Starting GraXpert..."), 0.0);
	g_spawn_async_with_pipes(NULL, argv, NULL,
			G_SPAWN_SEARCH_PATH |
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN | /*G_SPAWN_STDERR_TO_DEV_NULL |*/ G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, &child_pid, NULL, NULL,
			&child_stderr, &error);

	if (error != NULL) {
		siril_log_color_message(_("Spawning GraXpert failed: %s\n"), "red", error->message);
		return retval;
	}
	int i = 0;
	while (argv[i])
		g_free(argv[i++]);
	g_child_watch_add(child_pid, child_watch_cb, NULL);
	com.child_is_running = EXT_GRAXPERT;
#ifdef _WIN32
	com.childhandle = child_pid;		// For Windows, handle of a child process
#else
	com.childpid = child_pid;			// For other OSes, PID of a child process
#endif

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stderr), FALSE);
#else
	stream = g_unix_input_stream_new(child_stderr, FALSE);
#endif

	gboolean doprint = TRUE;
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
#ifdef GRAXPERT_DEBUG
	gchar *lastbuffer = NULL;
#endif
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
		gchar *arg = g_strstr_len(buffer, -1, progress_key);
		double value = -1.0;
		if (arg)
			value = g_ascii_strtod(arg + strlen(progress_key), NULL);
		if (value > 0.0 && value == value && verbose) {
			set_progress_bar_data(_("Running GraXpert"), value / 100.0);
		} else if (g_strrstr(buffer, "Finished")) {
			set_progress_bar_data(_("Done."), 1.0);
			retval = 0;
#ifdef GRAXPERT_DEBUG
			if (verbose)
				siril_log_message(_("GraXpert caught buffer: %s : exit successful\n"), buffer);
#endif
		} else if (doprint && verbose) {
			if (!g_strstr_len(buffer, -1, "ForkProcess")) { // These are not useful messages
				gchar *print_from = g_strstr_len(buffer, -1, "INFO");
				if (print_from) {
					print_from += 9;
					if (print_from) {
						siril_log_message("GraXpert: %s\n", print_from);
					}
				}
			}
		}
#ifdef GRAXPERT_DEBUG
		lastbuffer = g_strdup(buffer);
#endif
		g_free(buffer);
	}
	if (is_gui) {
		siril_log_message(_("GraXpert GUI finished.\n"));
		set_progress_bar_data(_("Done."), 1.0);
		retval = 0; // No "Finished" log message is printed when closing the GUI, so we assume success
	}
#ifdef GRAXPERT_DEBUG
	if (retval)
		siril_log_message(_("GraXpert exit not caught, last buffer: %s\n"), lastbuffer);
	g_free(lastbuffer);
#endif
	g_object_unref(data_input);
	g_object_unref(stream);
	if (!g_close(child_stderr, &error))
		siril_debug_print("%s\n", error->message);
	return retval;
}

static version_number graxpert_executablecheck(gchar* executable) {
	const gchar *version_key = "version: ";
	char *test_argv[3] = { NULL };
	gint child_stderr;
	g_autoptr(GError) error = NULL;
	version_number version = { 0 };
	if (!executable || executable[0] == '\0') {
		return version;
	}
	if (!g_file_test(executable, G_FILE_TEST_IS_EXECUTABLE)) {
		return version; // It's not executable so return a zero version number
	}

	// Change to the GraXpert installation directory
	gchar *dir = g_path_get_dirname(executable);
	gchar *currentdir = g_get_current_dir();
	int retval2 = g_chdir(dir);
	if (retval2) {
		g_free(dir);
		g_free(currentdir);
		return version;
	}
	g_free(dir);

	int nb = 0;
	test_argv[nb++] = executable;
	gchar *versionarg = g_strdup("-v");
	test_argv[nb++] = versionarg;
	// g_spawn handles wchar so not need to convert
	g_spawn_async_with_pipes(NULL, test_argv, NULL,
			G_SPAWN_SEARCH_PATH |
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN,
			NULL, NULL, NULL, NULL, NULL,
			&child_stderr, &error);

	if (error != NULL) {
		siril_log_color_message(_("Spawning GraXpert failed during version check: %s\n"), "red", error->message);
		g_free(versionarg);
		return version;
	}

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stderr), FALSE);
#else
	stream = g_unix_input_stream_new(child_stderr, FALSE);
#endif
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
	gboolean done = FALSE;
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL)) && !done) {
		// Find the start of the version substring
		gchar *version_start = g_strstr_len(buffer, -1, version_key);
		if (version_start) {
			version_start += strlen(version_key); // Move past "version: "

			// Find the end of the version substring
			gchar *version_end = g_strstr_len(version_start, -1, " ");
			*version_end = '\0';

			// Extract the version substring
			version = get_version_number_from_string(version_start);

			g_free(buffer);
			break;
		}
		g_free(buffer);
	}
	g_object_unref(data_input);
	g_object_unref(stream);
	g_free(versionarg);
	if (!g_close(child_stderr, &error))
		siril_debug_print("%s\n", error->message);

	// Change back to the prior working directory
	retval2 = g_chdir(currentdir);
	if (retval2) {
		retval2 = g_chdir(com.wd);
		if (retval2) {
			siril_log_color_message(_("Error: unable to change back to Siril working directory...\n"), "red");
		}
	}
	g_free(currentdir);

	return version;
}

void free_graxpert_data(graxpert_data *p) {
	g_free(p->path);
	free(p);
}

static gboolean end_graxpert(gpointer p) {
	stop_processing_thread();
	graxpert_data *args = (graxpert_data *) p;
	if (args->path) {
		open_single_image(args->path);
		if (args->backup_icc) {
			if (gfit.icc_profile)
				cmsCloseProfile(gfit.icc_profile);
			gfit.icc_profile = copyICCProfile(args->backup_icc);
			color_manage(&gfit, TRUE);
			cmsCloseProfile(args->backup_icc);
		}
	}
	free_graxpert_data(args);
	notify_gfit_modified();
	launch_clipboard_survey();
	return end_generic(NULL);
}

gpointer do_graxpert (gpointer p) {
	graxpert_data *args = (graxpert_data *) p;
	char *my_argv[64] = { 0 };
	int retval = 1;

	version_number graxpert_version = graxpert_executablecheck(com.pref.graxpert_path);
	if (compare_version(graxpert_version, (version_number) {.major_version = 3, .minor_version = 0, .micro_version = 2}) < 0) {
		siril_log_color_message(_("Error: GraXpert version is too old. You have version %d.%d.%d; at least version 3.0.2 is required.\n"), "red", graxpert_version.major_version, graxpert_version.minor_version, graxpert_version.micro_version);
		goto ERROR_OR_FINISHED;
	}

	// Configure input filename
	gchar *filename = NULL, *path = NULL;
	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		filename = g_path_get_basename(com.uniq->filename);
		gchar *temp = remove_ext_from_filename(filename);
		g_free(filename);
		filename = g_strdup_printf("graxpert_%s", temp);
		g_free(temp);
	} else {
		filename = g_strdup_printf("graxpert");
	}

	// We have to forcibly disable FITS compression as GraXpert cannot open compressed FITS (true as of GraXpert 3.0.2)
	gboolean pref_fitscomp = com.pref.comp.fits_enabled;
	com.pref.comp.fits_enabled = FALSE;
	gchar *temp = set_right_extension(filename);
	g_free(filename);
	filename = temp;

	path = g_build_filename(com.wd, filename, NULL);
	g_free(filename);
	// Save current image as input filename
	savefits(path, &gfit);
	com.pref.comp.fits_enabled = pref_fitscomp;

	// Configure GraXpert commandline
	int nb = 0;
	gboolean is_gui = FALSE;
	my_argv[nb++] = g_strdup(com.pref.graxpert_path);
	if (args->operation == GRAXPERT_BG) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("background-extraction");
		my_argv[nb++] = g_strdup("-correction");
		my_argv[nb++] = g_strdup_printf("%s", args->bg_mode == GRAXPERT_SUBTRACTION ? "Subtraction" : "Division");
		my_argv[nb++] = g_strdup("-smoothing");
		my_argv[nb++] = g_strdup_printf("%f", args->bg_smoothing);
		if (args->use_gpu) {
			my_argv[nb++] = g_strdup("-gpu");
			my_argv[nb++] = g_strdup("true");
		}
		if (args->keep_bg)
			my_argv[nb++] = g_strdup("-bg");
		my_argv[nb++] = g_strdup(" -output");
		my_argv[nb++] = g_strdup_printf("%s", path);
	} else if (args->operation == GRAXPERT_DENOISE) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("denoising");
		my_argv[nb++] = g_strdup_printf("-output");
		my_argv[nb++] = g_strdup_printf("%s", path);
		if (args->use_gpu) {
			my_argv[nb++] = g_strdup("-gpu");
			my_argv[nb++] = g_strdup("true");
		}
		my_argv[nb++] = g_strdup("-strength");
		my_argv[nb++] = g_strdup_printf("%.2f", args->denoise_strength);
	} else if (args->operation == GRAXPERT_GUI) {
		siril_log_message(_("GraXpert GUI will open with the current image. When you have finished, save the "
							"image and return to Siril. You will need to check and re-assign the ICC profile "
							"as GraXpert does not preserve ICC profiles embedded in FITS files.\n"));
		is_gui = TRUE;
	} else {
		siril_log_message(_("Error: unknown GraXpert operation\n"));
		g_free(my_argv[0]);
		goto ERROR_OR_FINISHED;
	}
	my_argv[nb++] = g_strdup_printf("%s", path);

	// Save a copy of the current ICC profile, as GraXpert does not preserve these
	if (gfit.icc_profile)
		args->backup_icc = copyICCProfile(gfit.icc_profile);

	// Execute GraXpert
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	retval = exec_prog_graxpert(my_argv, is_gui);
	gettimeofday(&t_end, NULL);
	if (verbose)
		show_time(t_start, t_end);
	if (retval)
		goto ERROR_OR_FINISHED;

	if (!is_gui) // In the GUI the user may change the saved file name, so we won't assume a filename
		args->path = g_strdup(path);
	g_free(path);

ERROR_OR_FINISHED:
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	siril_add_idle(end_graxpert, args); // this loads the result
	return GINT_TO_POINTER(retval);
}
