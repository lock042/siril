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
#include "algos/background_extraction.h"
#include "core/siril.h"
#include "core/icc_profile.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_update.h"
#include "core/siril_log.h"
#include "core/siril_spawn.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/siril_preview.h"
#include "gui/graxpert.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "yyjson.h"
#include "filters/graxpert.h"

// Uncomment the following line for highly verbose debugging messages
// #define GRAXPERT_DEBUG
// The following line keeps the config file
// #define GRAXPERT_CONFIG_DEBUG

// Define the minumum version numbers that support different operations
static const version_number min_bg_ver = { 3, 0, 0, 0 , FALSE, FALSE};
static const version_number min_denoise_ver = { 3, 0, 0, 0, FALSE, FALSE };
static const version_number min_deconv_ver = { 3, 1, 0, 0, FALSE, FALSE };

static gboolean verbose = TRUE;
static version_number graxpert_version = { 0 };
static gchar **background_ai_models = NULL;
static gchar **denoise_ai_models = NULL;
static gchar **deconv_ai_models = NULL;
static gchar **deconv_stellar_ai_models = NULL;
static gboolean graxpert_aborted = FALSE;
static GPid running_pid = (GPid) -1;

GPid get_running_graxpert_pid() {
	return running_pid;
}

void set_graxpert_aborted(gboolean state) {
	graxpert_aborted = state;
}

const gchar** get_ai_models(graxpert_operation operation) {
    return (const gchar**) (operation == GRAXPERT_DENOISE ? denoise_ai_models :
    					   (operation == GRAXPERT_DECONV ? deconv_ai_models :
    					   (operation == GRAXPERT_DECONV_STELLAR ? deconv_stellar_ai_models : background_ai_models)));
}

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	siril_debug_print("GraXpert exited with status %d\n", status);
	g_spawn_close_pid(pid);
}

// This ensures GraXpert is always called with a wide enough environment variable
// terminal width to prevent munging of the output we need to parse to get version
// information.

static GError *spawn_graxpert(gchar **argv, gint columns,
												GPid *child_pid, gint *stdin_fd,
												gint *stdout_fd, gint *stderr_fd) {
	GError *error = NULL;
	gchar **env = g_get_environ();
	gchar *columns_str = g_strdup_printf("%d", columns);

	env = g_environ_setenv(env, "COLUMNS", columns_str, TRUE);

	// On Windows, also set ANSICON_COLUMNS and ANSICON_LINES
	#ifdef G_OS_WIN32
	env = g_environ_setenv(env, "ANSICON_COLUMNS", columns_str, TRUE);
	#endif

	child_info *child = g_malloc(sizeof(child_info));
	gboolean spawn_result = siril_spawn_host_async_with_pipes(
		NULL,           // working directory
		argv,           // argument vector
		env,            // environment
		G_SPAWN_SEARCH_PATH | G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL,           // child setup function
		NULL,           // user data for child setup
		child_pid,      // child process id
		stdin_fd,       // stdin file descriptor
		stdout_fd,      // stdout file descriptor
		stderr_fd,      // stderr file descriptor
		&error
	);

	// At this point, remove the processing thread from the list of children and replace it
	// with the GraXpert process. This avoids tracking two children for the same task.
	if (get_thread_run())
		remove_child_from_children((GPid)-2);
	child->childpid = *child_pid;
	child->program = EXT_GRAXPERT;
	child->name = g_strdup("GraXpert");
	child->datetime = g_date_time_new_now_local();
	com.children = g_slist_prepend(com.children, child);

	// Set a static variable so we can recover the pid info from gui
	running_pid = *child_pid;
	g_child_watch_add(*child_pid, child_watch_cb, NULL);

	g_strfreev(env);
	g_free(columns_str);

	if (!spawn_result) {
		return error;
	}

	return NULL;
}

static int exec_prog_graxpert(char **argv, gboolean graxpert_no_exit_report, gboolean is_sequence) {
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

	if (!get_thread_run()) {
		return retval;
	}

	// g_spawn handles wchar so not need to convert
	if (!is_sequence) set_progress_bar_data(_("Starting GraXpert..."), 0.0);
	error = spawn_graxpert(argv, 200, &child_pid, NULL, NULL, &child_stderr);

	int i = 0;
	while (argv[i])
		g_free(argv[i++]);

	if (error != NULL) {
		siril_log_color_message(_("Spawning GraXpert failed: %s\n"), "red", error->message);
		return retval;
	}

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
#ifdef GRAXPERT_DEBUG
		siril_debug_print("%s\n", buffer);
#endif
		gchar *arg = g_strstr_len(buffer, -1, progress_key);
		double value = -1.0;
		gchar *errmsg = NULL;
		if (arg)
			value = g_ascii_strtod(arg + strlen(progress_key), NULL);
		if (value > 0.0 && value == value && verbose) {
			if (!is_sequence) set_progress_bar_data(_("Running GraXpert"), value / 100.0);
		} else if ( ((errmsg = g_strstr_len(buffer, -1, "ERROR")) && !graxpert_aborted) ) {
			if (!is_sequence) set_progress_bar_data(_("GraXpert reported an error"), max(value, 0.0) / 100);
			if (strlen(errmsg) > 9) {
				errmsg += 9;
				const gchar* color = g_strstr_len(buffer, -1, "Warning") ? "salmon" : "red";
				siril_log_color_message("GraXpert: %s\n", color, errmsg);
			}
			retval = 1;
		} else if (g_strrstr(buffer, "Finished") || g_strrstr(buffer, "finished")) {
			if (!is_sequence) set_progress_bar_data(_("Done."), 1.0);
			retval = 0;
#ifdef GRAXPERT_DEBUG
			if (verbose)
				siril_log_message(_("GraXpert caught buffer: %s : exit successful\n"), buffer);
#endif
		} else if (doprint && verbose) {
			if (!g_strstr_len(buffer, -1, "ForkProcess") && !g_strstr_len(buffer, -1, "ai_version")) { // These are not useful messages
				gchar *print_from = g_strstr_len(buffer, -1, "INFO");
				if (print_from) {
					if (strlen(print_from) > 9) {
						print_from += 9;
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
	// GraXpert has exited, remove from child list and reset the stored pid
	remove_child_from_children(child_pid);
	running_pid = (GPid) -1;
	if (graxpert_no_exit_report && retval == -1) {
		if (!is_sequence) siril_log_message(_("GraXpert GUI finished.\n"));
		if (!is_sequence) set_progress_bar_data(_("Done."), 1.0);
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

static gchar** parse_ai_versions(const char* version_line) {
	// Extract the versions string
	const char* versions_start = strchr(version_line, '[');
	if (!versions_start)
		return NULL;
	const char* versions_end = strchr(versions_start, ']');
	if (!versions_end || versions_start >= versions_end) {
		return NULL;
	}

	// Copy the versions string (excluding brackets)
	gchar* versions_str = g_strndup(versions_start + 1, versions_end - versions_start - 1);

	// Split the versions string into an array
	gchar** versions_array = g_strsplit(versions_str, ", ", -1);
	g_free(versions_str);

	// Count the number of versions
	int count = g_strv_length(versions_array);

	// Allocate the final array (including space for NULL terminator)
	gchar** result = g_malloc((count + 1) * sizeof(char*));

	// Copy each version string
	for (int i = 0; i < count; i++) {
		result[i] = g_strdup(g_strstrip(versions_array[i]));
		siril_debug_print("%s ", result[i]);
	}
	siril_debug_print("\n");
	result[count] = NULL;  // NULL-terminate the array

	g_strfreev(versions_array);

	return result;
}

static GMutex ai_version_check_mutex = { 0 };

static GError *spawn_graxpert_sync(gchar **argv, gint columns,
				GPid *child_pid, gint *exit_status, gchar **output) {
	GError *error = NULL;
	gchar **env = g_get_environ();
	gchar *columns_str = g_strdup_printf("%d", columns);

	env = g_environ_setenv(env, "COLUMNS", columns_str, TRUE);

	#ifdef G_OS_WIN32
	env = g_environ_setenv(env, "ANSICON_COLUMNS", columns_str, TRUE);
	#endif

	gboolean spawn_result = siril_spawn_host_sync(
		NULL,           // working directory
		argv,           // argument vector
		env,            // environment
		G_SPAWN_SEARCH_PATH, // g_spawn flags
		NULL,           // child setup function
		NULL,           // user data for child setup
		NULL,           // standard output
		output,         // standard error
		exit_status,    // exit status
		&error
	);

	g_strfreev(env);
	g_free(columns_str);

	if (!spawn_result) {
		return error;
	}

	return NULL;
}

gchar** ai_version_check(gchar* executable, graxpert_operation operation) {
	siril_debug_print("AI version check\n");
	g_mutex_lock(&ai_version_check_mutex);
	const gchar *key = "available remotely: [";
	gchar **result = NULL;
	g_autoptr(GError) error = NULL;
	version_number null_version = { 0 };

	if (memcmp(&graxpert_version, &null_version, sizeof(version_number))) {
		if (!executable || executable[0] == '\0') {
			siril_debug_print("No executable defined\n");
			g_mutex_unlock(&ai_version_check_mutex);
			return FALSE;
		}
		if (!g_file_test(executable, G_FILE_TEST_IS_EXECUTABLE)) {
			siril_debug_print("Executable indicates it is not executable\n");
			g_mutex_unlock(&ai_version_check_mutex);
			return FALSE;
		}

		gchar **test_argv = g_malloc_n(5, sizeof(gchar *));
		int nb = 0;
		test_argv[nb++] = executable;
		test_argv[nb++] = "-cmd";
		gchar *versionarg = g_strdup(operation == GRAXPERT_DENOISE ? "denoising" :
		(operation == GRAXPERT_DECONV ? "deconv-obj" :
		(operation == GRAXPERT_DECONV_STELLAR ? "deconv-stellar" : "background-extraction")));
		test_argv[nb++] = versionarg;
		test_argv[nb++] = "--help";
		test_argv[nb] = NULL;

		gint exit_status;
		gchar *output = NULL;
		error = spawn_graxpert_sync(test_argv, 500, NULL, &exit_status, &output);

		if (error != NULL) {
			siril_log_color_message(_("Spawning GraXpert failed during available AI model versions check: %s\n"), "red", error->message);
			g_free(versionarg);
			g_free(test_argv);
			g_mutex_unlock(&ai_version_check_mutex);
			return FALSE;
		}

		if (output) {
			// Find the start of the version substring
			gchar *start = g_strstr_len(output, -1, key);
			if (start) {
				siril_debug_print("Version string found for %d\n", (int) operation);
				result = parse_ai_versions(start);
			}
			g_free(output);
		}

		g_free(versionarg);
		g_free(test_argv);
	}
	g_mutex_unlock(&ai_version_check_mutex);
	return result;
}

typedef struct {
	graxpert_operation operation;
	version_number min_version;
	gchar **result;
} AIVersionCheckArgs;

static gpointer run_ai_version_check(gpointer data) {
	AIVersionCheckArgs *args = (AIVersionCheckArgs *)data;

	if (compare_version(args->min_version, graxpert_version) <= 0) {
		args->result = ai_version_check(com.pref.graxpert_path, args->operation);
	}

	return NULL;
}

void fill_graxpert_version_arrays() {
	GThread *bg_thread, *denoise_thread, *deconv_thread, *deconv_stellar_thread;

	AIVersionCheckArgs bg_args = {
		.operation = GRAXPERT_BG,
		.min_version = min_bg_ver,
		.result = NULL
	};

	AIVersionCheckArgs denoise_args = {
		.operation = GRAXPERT_DENOISE,
		.min_version = min_denoise_ver,
		.result = NULL
	};

	AIVersionCheckArgs deconv_args = {
		.operation = GRAXPERT_DECONV,
		.min_version = min_deconv_ver,
		.result = NULL
	};

	AIVersionCheckArgs deconv_stellar_args = {
		.operation = GRAXPERT_DECONV_STELLAR,
		.min_version = min_deconv_ver,
		.result = NULL
	};

	GError *error[4] = { NULL };
	bg_thread = g_thread_try_new("bg_check", run_ai_version_check, &bg_args, &error[0]);
	denoise_thread = g_thread_try_new("denoise_check", run_ai_version_check, &denoise_args, &error[1]);
	deconv_thread = g_thread_try_new("deconv_check", run_ai_version_check, &deconv_args, &error[2]);
	deconv_stellar_thread = g_thread_try_new("deconv_stellar_check", run_ai_version_check, &deconv_stellar_args, &error[3]);

	int errors = 0;
	for (int i = 0 ; i < 4 ; i++) {
		if (error[i]) {
			siril_log_color_message(_("Thread creation failed: %s\n"), "red", error[i]->message);
			errors++;
		}
	}
	if (errors) {
		return; // FALSE
	}

	g_thread_join(bg_thread);
	g_thread_join(denoise_thread);
	g_thread_join(deconv_thread);
	g_thread_join(deconv_stellar_thread);

	// Assign results to global variables
	background_ai_models = bg_args.result;
	denoise_ai_models = denoise_args.result;
	deconv_ai_models = deconv_args.result;
	deconv_stellar_ai_models = deconv_stellar_args.result;

	return; // TRUE
}

gboolean check_graxpert_version(const gchar *version, graxpert_operation operation) {
	if (version == NULL)
		return FALSE;
	if (operation == GRAXPERT_DENOISE) {
		if (denoise_ai_models == NULL)
			return FALSE;
		for (int i = 0 ; denoise_ai_models[i] != NULL ; i++) {
			if (!strcmp(version, denoise_ai_models[i]))
				return TRUE;
		}
	} else if (operation == GRAXPERT_DECONV) {
		if (deconv_ai_models == NULL)
			return FALSE;
		for (int i = 0 ; deconv_ai_models[i] != NULL ; i++) {
			if (!strcmp(version, deconv_ai_models[i]))
				return TRUE;
		}
	} else if (operation == GRAXPERT_DECONV_STELLAR) {
		if (deconv_stellar_ai_models == NULL)
			return FALSE;
		for (int i = 0 ; deconv_stellar_ai_models[i] != NULL ; i++) {
			if (!strcmp(version, deconv_stellar_ai_models[i]))
				return TRUE;
		}
	} else {
		if (background_ai_models == NULL)
			return FALSE;
		for (int i = 0 ; background_ai_models[i] != NULL ; i++) {
			if (!strcmp(version, background_ai_models[i]))
				return TRUE;
		}
	}
	return FALSE;
}

static GMutex graxpert_version_mutex = { 0 };

static gboolean graxpert_fetchversion(gchar* executable) {
	const gchar *version_key = "version: ";
	g_autoptr(GError) error = NULL;

	if (!executable || executable[0] == '\0') {
		return FALSE;
	}
	if (!g_file_test(executable, G_FILE_TEST_IS_EXECUTABLE)) {
		return FALSE;
	}

	g_mutex_lock(&graxpert_version_mutex);

	gchar **test_argv = g_malloc_n(3, sizeof(gchar *));
	int nb = 0;
	test_argv[nb++] = executable;
	gchar *versionarg = g_strdup("-v");
	test_argv[nb++] = versionarg;
	test_argv[nb] = NULL;

	gint exit_status;
	gchar *output = NULL;
	error = spawn_graxpert_sync(test_argv, 200, NULL, &exit_status, &output);

	if (error != NULL) {
		siril_log_color_message(_("Spawning GraXpert failed during version check: %s\n"), "red", error->message);
		g_free(versionarg);
		g_free(test_argv);
		g_mutex_unlock(&graxpert_version_mutex);
		return FALSE;
	}

	if (output) {
		// Find the start of the version substring
		gchar *version_start = g_strstr_len(output, -1, version_key);
		if (version_start) {
			version_start += strlen(version_key); // Move past "version: "
			// Find the end of the version substring
			gchar *version_end = g_strstr_len(version_start, -1, " ");
			if (version_end) {
				*version_end = '\0';
				// Extract the version substring
				graxpert_version = get_version_number_from_string(version_start);
			}
		}
		g_free(output);
	}

	g_free(versionarg);
	g_free(test_argv);
	g_mutex_unlock(&graxpert_version_mutex);
	return TRUE;
}

gboolean graxpert_executablecheck(gchar* executable, graxpert_operation operation) {
	version_number null_version = { 0 };
	if (!memcmp(&graxpert_version, &null_version, sizeof(version_number))) {
		if (!graxpert_fetchversion(executable)) {
			return FALSE;
		}
	}

	if (compare_version(graxpert_version, (version_number) {.major_version = 3, .minor_version = 0, .micro_version = 0}) < 0) {
		siril_log_color_message(_("Error: GraXpert version is too old. You have version %d.%d.%d; at least version 3.0.0 is required.\n"), "red",
				graxpert_version.major_version, graxpert_version.minor_version, graxpert_version.micro_version);
		return FALSE;
	} else {
		if (compare_version(graxpert_version, (version_number ) {.major_version = 3, .minor_version = 1, .micro_version = 0 }) < 0
				&& (operation == GRAXPERT_DECONV || operation == GRAXPERT_DECONV_STELLAR)) {
			return FALSE;
		}
		return TRUE;
	}
}

gpointer graxpert_setup_async(gpointer user_data) {
	if (graxpert_fetchversion(com.pref.graxpert_path)) {
		siril_debug_print("GraXpert version %d.%d.%d found\n", graxpert_version.major_version, graxpert_version.minor_version, graxpert_version.micro_version);
		fill_graxpert_version_arrays();
		siril_debug_print("GraXpert AI model arrays populated\n");
		version_number null_version = { 0 };
		if (!com.headless && memcmp(&graxpert_version, &null_version, sizeof(version_number))) {
			// initialize widgets in the GTK thread
			siril_add_idle(initialize_graxpert_widgets_if_needed, GINT_TO_POINTER(1));
		}
	} else {
		g_strfreev(background_ai_models);
		background_ai_models = NULL;
		g_strfreev(denoise_ai_models);
		denoise_ai_models = NULL;
		g_strfreev(deconv_ai_models);
		deconv_ai_models = NULL;
		g_strfreev(deconv_stellar_ai_models);
		deconv_stellar_ai_models = NULL;
	}
	return GINT_TO_POINTER(0);
}

void ai_versions_to_log(graxpert_operation operation) {
	gchar** array = operation == GRAXPERT_DENOISE ? denoise_ai_models :
			(operation == GRAXPERT_DECONV ? deconv_ai_models :
			(operation == GRAXPERT_DECONV_STELLAR ? deconv_stellar_ai_models : background_ai_models));
	if (!array) {
		siril_log_message(_("None!\n"));
	} else {
		for (int i = 0 ; array[i] ; i++) {
			siril_log_message("%s\n", array[i]);
		}
	}
}

gboolean save_graxpert_config(graxpert_data *args) {
	const gchar *filename = args->configfile;
	gboolean success = FALSE;

	// Create JSON document
	yyjson_mut_doc *doc = yyjson_mut_doc_new(NULL);
	yyjson_mut_val *root = yyjson_mut_obj(doc);
	yyjson_mut_doc_set_root(doc, root);

	// Add working directory
	yyjson_mut_obj_add_str(doc, root, "working_dir", com.wd);

	// Add width and height
	yyjson_mut_obj_add_int(doc, root, "width", args->fit->rx);
	yyjson_mut_obj_add_int(doc, root, "height", args->fit->ry);

	// Add background points
	if (args->bg_samples) {
		yyjson_mut_val *bg_points = yyjson_mut_arr(doc);
		yyjson_mut_obj_add_val(doc, root, "background_points", bg_points);

		for (GSList *l = args->bg_samples; l != NULL; l = l->next) {
			background_sample *s = (background_sample *)l->data;
			yyjson_mut_val *point = yyjson_mut_arr(doc);

			yyjson_mut_arr_add_int(doc, point, min(args->fit->rx - 1, round_to_int(s->position.x)));
			yyjson_mut_arr_add_int(doc, point, min(args->fit->ry - 1, round_to_int(s->position.y)));
			yyjson_mut_arr_add_int(doc, point, 1);

			yyjson_mut_arr_append(bg_points, point);
		}
	}

	// Add other options
	yyjson_mut_obj_add_int(doc, root, "bg_pts_option", args->bg_pts_option);

	// Add stretch option
	const char *stretch_str;
	switch (args->stretch_option) {
		case STRETCH_OPTION_10_BG_3_SIGMA: stretch_str = "10% Bg, 3 sigma"; break;
		case STRETCH_OPTION_15_BG_3_SIGMA: stretch_str = "15% Bg, 3 sigma"; break;
		case STRETCH_OPTION_20_BG_3_SIGMA: stretch_str = "20% Bg, 3 sigma"; break;
		case STRETCH_OPTION_30_BG_2_SIGMA: stretch_str = "30% Bg, 2 sigma"; break;
		default: stretch_str = "No Stretch";
	}
	yyjson_mut_obj_add_str(doc, root, "stretch_option", stretch_str);

	// Add bg tolerance option (with clamping)
	double bg_tol = args->bg_tol_option;
	bg_tol = fmin(fmax(bg_tol, -2.0), 10.0);
	yyjson_mut_obj_add_real(doc, root, "bg_tol_option", bg_tol);

	// Add interpolation type
	const char *interpol_str;
	switch (args->bg_algo) {
		case GRAXPERT_BG_KRIGING: interpol_str = "Kriging"; break;
		case GRAXPERT_BG_RBF: interpol_str = "RBF"; break;
		case GRAXPERT_BG_SPLINE: interpol_str = "Splines"; break;
		default: interpol_str = "AI";
	}
	yyjson_mut_obj_add_str(doc, root, "interpol_type_option", interpol_str);

	yyjson_mut_obj_add_real(doc, root, "smoothing_option", args->bg_smoothing);

	if (args->ai_batch_size != -1) {
		yyjson_mut_obj_add_int(doc, root, "ai_batch_size", args->ai_batch_size);
	}

	if (args->sample_size != -1) {
		yyjson_mut_obj_add_int(doc, root, "sample_size", args->sample_size);
	}

	if (args->spline_order != -1) {
		yyjson_mut_obj_add_int(doc, root, "spline_order", args->spline_order);
	}

	yyjson_mut_obj_add_str(doc, root, "corr_type",
						   args->bg_mode == GRAXPERT_SUBTRACTION ? "Subtraction" : "Division");

	const char *kernel_str;
	switch (args->kernel) {
		case GRAXPERT_THIN_PLATE: kernel_str = "thin_plate"; break;
		case GRAXPERT_QUINTIC: kernel_str = "quintic"; break;
		case GRAXPERT_CUBIC: kernel_str = "cubic"; break;
		default: kernel_str = "linear";
	}
	yyjson_mut_obj_add_str(doc, root, "RBF_kernel", kernel_str);

	yyjson_mut_obj_add_real(doc, root, "denoise_strength", args->denoise_strength);
	yyjson_mut_obj_add_real(doc, root, "deconvolution_strength", args->deconv_strength);
	yyjson_mut_obj_add_real(doc, root, "deconvolution_psfsize", args->deconv_blur_psf_size);

	yyjson_mut_obj_add_str(doc, root, "saveas_option",
						   args->fit->type == DATA_FLOAT ? "32 bit Fits" : "16 bit Fits");

	// Write to file
	success = yyjson_mut_write_file(filename, doc, YYJSON_WRITE_PRETTY, NULL, NULL);

	// Cleanup
	yyjson_mut_doc_free(doc);

	return success;
}

graxpert_data *new_graxpert_data() {
	graxpert_data *p = calloc(1, sizeof(graxpert_data));
	p->fit = &gfit;
	p->operation = GRAXPERT_BG;
	p->bg_smoothing = 0.5;
	p->bg_algo = GRAXPERT_BG_AI;
	p->bg_mode = GRAXPERT_SUBTRACTION;
	p->kernel = GRAXPERT_THIN_PLATE;
	p->stretch_option = STRETCH_OPTION_10_BG_3_SIGMA;
	p->sample_size = 25;
	p->spline_order = 3;
	p->bg_tol_option = 2;
	p->denoise_strength = 0.8;
	p->deconv_strength = 0.5;
	p->deconv_blur_psf_size = 0.3;
	p->use_gpu = TRUE;
	p->ai_batch_size = 4;
	p->bg_pts_option = 15;
	return p;
}

void free_graxpert_data(graxpert_data *args) {
	g_free(args->path);
	g_free(args->configfile);
	g_free(args->ai_version);
	free_background_sample_list(args->bg_samples);
	if (args->backup_icc)
		cmsCloseProfile(args->backup_icc);
	free(args);
}

static void open_graxpert_result(graxpert_data *args) {
	// Clean up config file if one was used
#ifndef GRAXPERT_CONFIG_DEBUG
	if (args->configfile && g_unlink(args->configfile))
		siril_debug_print("Failed to remove GraXpert config file\n");
#endif
	// If successful, open the result image
	/* Note: we do this even when sequence working: it's a bit inefficient to read it in,
	 * delete it and save it again but it works for all sequence types (incl. ser and
	 * FITSEQ) whereas simply moving the file saved by GraXpert would only work for
	 * sequences of individual FITS files. */
	if (args->path) {
		disable_profile_check_verbose();
		if (args->previewing) {
			siril_log_message(_("Reading result from the temporary working file...\n"));
			if (readfits(args->path, &gui.roi.fit, NULL, (gui.roi.fit.type == DATA_FLOAT))) {
				siril_log_color_message(_("Error opening GraXpert result. Check the file %s\n"), "red", args->path);
				enable_profile_check_verbose();
				goto END_AND_RETURN;
			}
		} else {
			fits *result = calloc(1, sizeof(fits));
			siril_log_message(_("Reading result from the temporary working file...\n"));
			if (readfits(args->path, result, NULL, !com.pref.force_16bit)) {
				siril_log_color_message(_("Error opening GraXpert result. Check the file %s\n"), "red", args->path);
				enable_profile_check_verbose();
				goto END_AND_RETURN;
			}
			// Check the result dimensions and bit depth match what we expect
			// This should never fail, but we shouldn't trust external software too
			// much and a dimension mismatch would result in a crash
			if (args->fit->rx != result->rx || args->fit->ry != result->ry || args->fit->naxes[2] != result->naxes[2]) {
				siril_log_color_message(_("Error: the GraXpert image dimensions and bit depth "
						"do not match those of the original image. Please report this as a bug.\n"), "red");
				goto END_AND_RETURN;
			}
			// Swap the image data pointers. This way the data originally in args->fit
			// gets freed by the call to clearfits(result) and the result data becomes
			// owned by args->fit
			if (fits_swap_image_data(args->fit, result)) {
				siril_debug_print("Error, NULL pointer passed to fits_swap_image_data\n");
				goto END_AND_RETURN;
			}
			if (args->fit->type == DATA_FLOAT && com.pref.force_16bit) {
				size_t npixels = args->fit->rx * args->fit->ry * args->fit->naxes[2];
				WORD *newbuf = malloc(npixels * sizeof(WORD));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (npixels > 50000)
#endif
				for (size_t i = 0 ; i < npixels; i++) {
					newbuf[i] = roundf_to_WORD(args->fit->fdata[i] * USHRT_MAX_SINGLE);
				}
				fit_replace_buffer(args->fit, newbuf, DATA_USHORT);
			}
			clearfits(result);
			free(result);
			if (args->fit == &gfit) {
				copy_gfit_to_backup();
				populate_roi();
			}
		}
		siril_log_message(_("Removing the temporary working file...\n"));
		if (g_unlink(args->path))
			siril_debug_print("Failed to unlink GraXpert working file\n");
		enable_profile_check_verbose();
	}
END_AND_RETURN:
	unlock_roi_mutex();
	return;
}

static gboolean end_graxpert(gpointer p) {
	stop_processing_thread();
	open_graxpert_result((graxpert_data *) p);
	free_graxpert_data((graxpert_data *) p);
	notify_gfit_modified();
	gui_function(launch_clipboard_survey, NULL);
	return end_generic(NULL);
}

gpointer do_graxpert (gpointer p) {
	lock_roi_mutex();
	set_graxpert_aborted(FALSE);
	gchar *text = NULL;
	graxpert_data *args = (graxpert_data *) p;
	if (args->fit == &gfit)
		copy_backup_to_gfit();
	if (!args->previewing && !com.script) {
		switch (args->operation) {
			case GRAXPERT_BG:
				text = g_strdup_printf(_("GraXpert BG extraction, smoothness %.3f"), args->bg_smoothing);
				break;
			case GRAXPERT_DENOISE:
				text = g_strdup_printf(_("GraXpert denoising, strength %.3f"), args->denoise_strength);
				break;
			case GRAXPERT_DECONV:
			case GRAXPERT_DECONV_STELLAR:
				text = g_strdup_printf(_("GraXpert deconv, strength %.3f, psf size %.3f"), args->deconv_strength, args->deconv_blur_psf_size);
				break;
			default:
				text = g_strdup(_("GraXpert operations using GUI"));
		}
	}
	char *my_argv[64] = { 0 };
	gchar *filename = NULL, *path = NULL, *outpath = NULL;
	int retval = 1;
	if (!graxpert_executablecheck(com.pref.graxpert_path, args->operation))
		goto ERROR_OR_FINISHED;

	// Configure input filename
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
	// We also have to force the extension ".fits" to avoid GraXpert adding it onto other extensions
	gboolean pref_fitscomp = com.pref.comp.fits_enabled;
	gchar *backup_ext = g_strdup(com.pref.ext);
	g_free(com.pref.ext);
	com.pref.ext = g_strdup(".fits");
	com.pref.comp.fits_enabled = FALSE;
	gchar *temp = set_right_extension(filename);
	g_free(filename);
	filename = temp;
	g_free(com.pref.ext);
	com.pref.ext = backup_ext;

	path = g_build_filename(com.wd, filename, NULL);
	temp = g_strdup(path);
	outpath = remove_ext_from_filename(temp);
	g_free(temp);
	g_free(filename);
	// Save current image as input filename
	siril_log_message(_("Saving temporary working file...\n"));
	if (savefits(path, args->fit)) {
		siril_log_color_message(_("Error: failed to save temporary FITS\n"), "red");
		goto ERROR_OR_FINISHED;
	}

	com.pref.comp.fits_enabled = pref_fitscomp;

	// Configure GraXpert commandline
	int nb = 0;
	gboolean graxpert_no_exit_report = FALSE;
	gboolean is_gui = FALSE;
	my_argv[nb++] = g_strdup(com.pref.graxpert_path);
	if (args->operation == GRAXPERT_BG) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("background-extraction");
		if (args->bg_algo == GRAXPERT_BG_AI) {
			my_argv[nb++] = g_strdup("-correction");
			my_argv[nb++] = g_strdup_printf("%s", args->bg_mode == GRAXPERT_SUBTRACTION ? "Subtraction" : "Division");
			my_argv[nb++] = g_strdup("-smoothing");
			my_argv[nb++] = g_strdup_printf("%f", args->bg_smoothing);
			my_argv[nb++] = g_strdup("-gpu");
			my_argv[nb++] = g_strdup(args->use_gpu ? "true" : "false");
			if (args->ai_version != NULL) {
				my_argv[nb++] = g_strdup("-ai_version");
				my_argv[nb++] = g_strdup(args->ai_version);
			}
		} else {
			if (!args->bg_samples) {
				siril_log_color_message(_("Background samples must be computed for GraXpert RBF, Spline and Kriging methods\n"), "red");
				goto ERROR_OR_FINISHED;
			}
			my_argv[nb++] = g_strdup("-preferences_file");
			args->configfile = g_build_filename(com.wd, "siril-graxpert.pref", NULL);
			my_argv[nb++] = g_strdup(args->configfile);
			save_graxpert_config(args);
		}
		if (args->keep_bg)
			my_argv[nb++] = g_strdup("-bg");
		my_argv[nb++] = g_strdup("-output");
		my_argv[nb++] = g_strdup_printf("%s", outpath);
		graxpert_no_exit_report = TRUE;
	} else if (args->operation == GRAXPERT_DENOISE) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("denoising");
		my_argv[nb++] = g_strdup_printf("-output");
		my_argv[nb++] = g_strdup_printf("%s", outpath);
		my_argv[nb++] = g_strdup("-gpu");
		my_argv[nb++] = g_strdup(args->use_gpu ? "true" : "false");
		my_argv[nb++] = g_strdup("-strength");
		my_argv[nb++] = g_strdup_printf("%.2f", args->denoise_strength);
		if (args->ai_version != NULL) {
			my_argv[nb++] = g_strdup("-ai_version");
			my_argv[nb++] = g_strdup(args->ai_version);
		}
	} else if (args->operation == GRAXPERT_DECONV || args->operation == GRAXPERT_DECONV_STELLAR) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = args->operation == GRAXPERT_DECONV ? g_strdup("deconv-obj") : g_strdup("deconv-stellar");
		my_argv[nb++] = g_strdup_printf("-output");
		my_argv[nb++] = g_strdup_printf("%s", outpath);
		my_argv[nb++] = g_strdup("-gpu");
		my_argv[nb++] = g_strdup(args->use_gpu ? "true" : "false");
		my_argv[nb++] = g_strdup("-strength");
		my_argv[nb++] = g_strdup_printf("%.2f", args->deconv_strength);
		my_argv[nb++] = g_strdup("-psfsize");
		my_argv[nb++] = g_strdup_printf("%.2f", args->deconv_blur_psf_size);
		if (args->ai_version != NULL) {
			my_argv[nb++] = g_strdup("-ai_version");
			my_argv[nb++] = g_strdup(args->ai_version);
		}
	} else if (args->operation == GRAXPERT_GUI) {
		siril_log_message(_("GraXpert GUI will open with the current image. When you have finished, save the "
							"image and return to Siril.\nYou will need to check and re-assign the ICC profile "
							"as GraXpert does not preserve ICC profiles embedded in FITS files.\n"));
		graxpert_no_exit_report = TRUE;
		is_gui = TRUE;
	} else {
		siril_log_message(_("Error: unknown GraXpert operation\n"));
		g_free(my_argv[0]);
		goto ERROR_OR_FINISHED;
	}
	my_argv[nb++] = g_strdup_printf("%s", path);

	// Save a copy of the current ICC profile, as GraXpert does not preserve these
	if (args->fit->icc_profile)
		args->backup_icc = copyICCProfile(args->fit->icc_profile);

	// Execute GraXpert
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	retval = exec_prog_graxpert(my_argv, graxpert_no_exit_report, args->seq != NULL);
	gettimeofday(&t_end, NULL);
	if (verbose)
		show_time(t_start, t_end);
	if (retval) {
		if (g_unlink(path))
			siril_debug_print("Error removing temporary file %s\n", path);
		goto ERROR_OR_FINISHED;
	}

	if (!is_gui) {
		if (args->path)
			g_free(args->path);
		args->path = g_strdup(path);
	}
	g_free(path);

ERROR_OR_FINISHED:
	g_free(outpath);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	if (!retval && text) {
		undo_save_state(&gfit, text);
	}
	g_free(text);
	if (!args->seq && !com.script)
		siril_add_idle(end_graxpert, args); // this loads the result
	else
		open_graxpert_result(args);
	return GINT_TO_POINTER(retval);
}

/******** SEQUENCE FUNCTIONS *********/

static int graxpert_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
//	GraXpert cannot run in parallel as it fully utilizes the GPU / CPU. This function therefore
//	returns a maximum of 1 and all images will be processed in series.
	unsigned int MB_per_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	if (limit > 0) {
		required = MB_per_image;
		limit = MB_avail / required;
	}
	limit = (limit >= 1 ? 1 : 0);
	return limit;
}

static int graxpert_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 0;
	graxpert_data *data = (graxpert_data *) args->user;
	if (data->operation == GRAXPERT_BG && data->bg_algo != GRAXPERT_BG_AI) {
		const char *err;
		data->bg_samples = generate_samples(fit, data->bg_pts_option, data->bg_tol_option, data->sample_size, &err, MULTI_THREADED);
		if (!data->bg_samples) {
			siril_log_color_message(_("Failed to generate background samples for image: %s\n"), "red", _(err));
			return 1;
		}
	}
	data->fit = fit;
	siril_log_color_message(_("GraXpert: Processing image %d\n"), "green", o + 1);
	verbose = FALSE;
	do_graxpert(data);
	verbose = TRUE;
	free_background_sample_list(data->bg_samples);
	return ret;
}

void apply_graxpert_to_sequence(graxpert_data *args) {
	if (!graxpert_executablecheck(com.pref.graxpert_path, args->operation)) {
		return;
	}
	args->fit = NULL;
	struct generic_seq_args *seqargs = create_default_seqargs(args->seq);
	seqargs->seq = args->seq;
	seqargs->parallel = FALSE;
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = args->seq->selnum;
	seqargs->compute_mem_limits_hook = graxpert_compute_mem_limits;
	seqargs->prepare_hook = seq_prepare_hook;
	seqargs->image_hook = graxpert_image_hook;
	seqargs->description = _("GraXpert");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	if (args->operation == GRAXPERT_BG)
		seqargs->new_seq_prefix = strdup("gxbg_");
	else if (args->operation == GRAXPERT_DENOISE)
		seqargs->new_seq_prefix = strdup("gxnr_");
	else if (args->operation == GRAXPERT_DECONV || args->operation == GRAXPERT_DECONV_STELLAR)
		seqargs->new_seq_prefix = strdup("gxdec_");
	else  {
		free_graxpert_data(args);
		free_generic_seq_args(seqargs);
		return;
	}
	seqargs->load_new_sequence = TRUE;
	seqargs->user = args;
	set_progress_bar_data(_("GraXpert: Processing..."), 0.);
	if (!start_in_new_thread(generic_sequence_worker, seqargs)) {
		free_graxpert_data(args);
		free_generic_seq_args(seqargs);
	}
}
