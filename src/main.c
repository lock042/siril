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

#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

// #define DEBUG_MAIN

#include <gsl/gsl_errno.h>
#include <gtk/gtk.h>
#include <gtksourceview/gtksource.h>
#include <stdio.h>
#include <string.h>
#include <locale.h>
#include <unistd.h>
#include <fftw3.h>
#ifdef OS_OSX
#import <AppKit/AppKit.h>
#if defined(ENABLE_RELOCATABLE_RESOURCES)
#include <sys/param.h> /* PATH_MAX */
#include <libgen.h> /* dirname */
#include <sys/stat.h>
#endif /* ENABLE_RELOCATABLE_RESOURCES */
#elif _WIN32
#include <windows.h>
#endif

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "siril_resource.h"
#include "git-version.h"
#include "core/siril.h"
#include "core/icc_profile.h"
#include "core/proto.h"
#include "algos/siril_random.h"
#include "core/siril_actions.h"
#include "core/initfile.h"
#include "core/command_line_processor.h"
#include "core/command.h"
#include "core/pipe.h"
#include "core/signals.h"
#include "core/siril_app_dirs.h"
#include "core/siril_language.h"
#include "core/siril_networking.h"
#include "io/siril_pythonmodule.h"
#include "core/siril_update.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "io/siril_git.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "gui/ui_files.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/siril_css.h"
#include "registration/registration.h"


/* the global variables of the whole project */
cominfo com = { 0 };	// the core data struct
guiinfo gui = { 0 };	// the gui data struct
fits gfit;	// currently loaded image

static gchar *main_option_directory = NULL;
static gchar *main_option_script = NULL;
static gchar *main_option_initfile = NULL;
static gchar *main_option_rpipe_path = NULL;
static gchar *main_option_wpipe_path = NULL;
static gboolean main_option_pipe = FALSE;

static gboolean _print_version_and_exit(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
#ifdef SIRIL_UNSTABLE
	g_print("%s %s-%s\n", PACKAGE, VERSION, SIRIL_GIT_VERSION_ABBREV);
#else
	g_print("%s %s\n", PACKAGE, VERSION);
#endif
	exit(EXIT_SUCCESS);
	return TRUE;
}

static gboolean _print_copyright_and_exit(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
	g_print("Copyright © 2012-%s team free-astro\n", SIRIL_GIT_LAST_COMMIT_YEAR);
	exit(EXIT_SUCCESS);
	return TRUE;
}

static gboolean _print_list_of_formats_and_exit(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
	list_format_available();
	exit(EXIT_SUCCESS);
	return TRUE;
}

static gboolean _set_offline(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
	set_online_status(FALSE);
	return TRUE;
}

static GOptionEntry main_option[] = {
	{ "directory", 'd', 0, G_OPTION_ARG_FILENAME, &main_option_directory, N_("changing the current working directory as the argument"), NULL },
	{ "script", 's', 0, G_OPTION_ARG_FILENAME, &main_option_script, N_("run the siril commands script in console mode. If argument is equal to \"-\", then siril will read stdin input"), NULL },
	{ "initfile", 'i', 0, G_OPTION_ARG_FILENAME, &main_option_initfile, N_("load configuration from file name instead of the default configuration file"), NULL },
	{ "pipe", 'p', 0, G_OPTION_ARG_NONE, &main_option_pipe, N_("run in console mode with command and log stream through named pipes"), NULL },
	{ "inpipe", 'r', 0, G_OPTION_ARG_FILENAME, &main_option_rpipe_path, N_("specify the path for the read pipe, the one receiving commands"), NULL },
	{ "outpipe", 'w', 0, G_OPTION_ARG_FILENAME, &main_option_wpipe_path, N_("specify the path for the write pipe, the one outputing messages"), NULL },
	{ "format", 'f', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_list_of_formats_and_exit, N_("print all supported image file formats (depending on installed libraries)" ), NULL },
	{ "offline", 'o', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _set_offline, N_("start in offline mode"), NULL },
	{ "version", 'v', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_version_and_exit, N_("print the application’s version"), NULL},
	{ "copyright", 'c', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_copyright_and_exit, N_("print the copyright"), NULL},
	{ NULL },
};

static GActionEntry app_entries[] = {
	{ "quit", quit_action_activate },
	{ "preferences", preferences_action_activate },
	{ "open",  open_action_activate },
	{ "save", save_action_activate },
	{ "save-as", save_as_action_activate },
	{ "about", about_action_activate }
};

gboolean builder_add_from_resource_with_replace(GtkBuilder *builder, const gchar* resource_path, GError** error) {
	// Load the resource data
	GBytes* resource_data = g_resources_lookup_data(resource_path, G_RESOURCE_LOOKUP_FLAGS_NONE, error);
	if (!resource_data) {
		return FALSE;
	}

	// Get the data as a string
	gsize data_size;
	const gchar* data = g_bytes_get_data(resource_data, &data_size);
#ifdef __APPLE__
	// Create a GString from the data
	GString *str = g_string_new_len(data, data_size);

	// Perform the replacement - no limit on number of replacements
	g_string_replace(str, "GDK_CONTROL_MASK", "GDK_META_MASK", 0);

	// Get the modified result
	const gchar *result = str->str;
#else
	const gchar* result = data;
#endif

	// Create and load the builder
	gboolean success = gtk_builder_add_from_string(builder, result, -1, error);

	// Clean up
#ifdef __APPLE__
	g_string_free(str, TRUE);  // TRUE to free the string data as well
#endif
	g_bytes_unref(resource_data);

	if (!success) {
		g_object_unref(builder);
		return FALSE;
	}

	return TRUE;
}

void load_ui_files() {
	GError *err = NULL;
	gboolean retval;

	/* try to load the first UI file, from the sources defined above */
	gui.builder = gtk_builder_new_from_resource(ui_files[0]);
	if (!gui.builder) {
		g_error(_("%s was not found or contains errors, "
					"cannot render GUI.\n Exiting.\n"), ui_files[0]);
		exit(EXIT_FAILURE);
	}
#ifdef DEBUG_MAIN
	siril_debug_print("Successfully loaded '%s'\n", ui_files[0]);
#endif

	uint32_t i = 1;
	while (*ui_files[i]) {

		/* try to load each successive UI file from ui_files */
		retval = gtk_builder_add_from_resource(gui.builder, ui_files[i], &err);
		if (!retval) {
			g_error(_("%s was not found or contains errors, "
						"cannot render GUI:\n%s\n Exiting.\n"), ui_files[i], err->message);
			g_clear_error(&err);
			exit(EXIT_FAILURE);
		}
#ifdef DEBUG_MAIN
		siril_debug_print("Successfully loaded '%s'\n", ui_files[i]);
#endif
		i++;
	}
	/* Now we load any UI files that use the primary accelerator.
	   These get processed specially so that on MacOS GDK_CONTROL_MASK is replaced
	   with GDK_META_MASK so that shortcuts work with Cmd as expected on that OS. */
	i = 0;
	while (*ui_files_with_primary_accelerator[i]) {
		retval = builder_add_from_resource_with_replace(gui.builder, ui_files_with_primary_accelerator[i], &err);
		retval = gtk_builder_add_from_resource(gui.builder, ui_files[i], &err);
		if (!retval) {
			g_error(_("%s was not found or contains errors, "
			"cannot render GUI:\n%s\n Exiting.\n"), ui_files[i], err->message);
			g_clear_error(&err);
			exit(EXIT_FAILURE);
		}
		#ifdef DEBUG_MAIN
		siril_debug_print("Successfully loaded '%s'\n", ui_files[i]);
		#endif
		i++;
	}
}

#if defined (HAVE_FFTW3F_THREADS) && (_OPENMP)
void parallel_loop(void *(*work)(char *), char *jobdata, size_t elsize, int njobs, void *data)
{
#pragma omp parallel for
	for (int i = 0; i < njobs; ++i)
		work(jobdata + elsize * i);
}
#endif

static void global_initialization() {
	gsl_set_error_handler_off();
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.tilt = NULL;
	com.uniq = NULL;
	com.kernel = NULL;
	com.kernelsize = 0;
	com.kernelchannels = 0;
	memset(&com.spcc_data, 0, sizeof(struct spcc_data_store));
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));
	gui.selected_star = -1;
	gui.repo_scripts = NULL;
	gui.qphot = NULL;
	gui.draw_extra = NULL;
	gui.cvport = RED_VPORT;
	gui.show_excluded = TRUE;
	gui.sliders = MINMAX;
	gui.zoom_value = ZOOM_DEFAULT;
	gui.ratio = 0.0;
	gui.use_hd_remap = FALSE;
	for (int i = 0; i < 3 ; i++)
		gui.hd_remap_index[i] = NULL;

	initialize_default_settings();	// com.pref
#ifdef HAVE_FFTW3F_MULTITHREAD
	fprintf(stdout, _("Initializing FFTW multithreading support...\n"));
	fftwf_init_threads(); // Should really only be called once so do it at startup
#endif
#if defined (HAVE_FFTW3F_THREADS) && (_OPENMP)
	// If we are using FFTW built against pthreads but are using OpenMP for other aspects of the
	// program, replace the parallel loop to avoid competing threads.
	// See https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
	void fftw_threads_set_callback(
		void (*parallel_loop)(void *(*work)(char *), char *jobdata, size_t elsize, int njobs,
		void *data), void *data);
#endif

#ifdef _OPENMP
	omp_set_max_active_levels(2);
#endif

}

static void siril_app_startup(GApplication *application) {
	signals_init();
	/*
	 * Force C locale for numbers to avoid "," being used as decimal separator.
	 * Called here and not in main() because setlocale(LC_ALL, "") is called as
	 * part of g_application_run().
	 */
	setlocale(LC_NUMERIC, "C");

	g_set_application_name(PACKAGE_NAME);
	gtk_window_set_default_icon_name("siril");
	g_application_set_resource_base_path(application, "/org/siril/Siril/pixmaps/");

	g_action_map_add_action_entries(G_ACTION_MAP(application), app_entries,
			G_N_ELEMENTS(app_entries), application);
	// Initialize GtkSource
	gtk_source_init();
}

static void siril_app_activate(GApplication *application) {
	/* the first thing we need to do is to know if we are headless or not */
	if (main_option_script || main_option_pipe) {
		com.script = TRUE;
		com.headless = TRUE;
	}

#if defined(_WIN32) && !defined(SIRIL_UNSTABLE)
	ShowWindow(GetConsoleWindow(), SW_MINIMIZE); //hiding the console
#endif
	siril_initialize_rng();
	global_initialization();
	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);

	siril_log_color_message(_("Welcome to %s v%s\n"), "bold", PACKAGE, VERSION);

	/* initialize converters (utilities used for different image types importing) */
	gchar *supported_files = initialize_converters();

	if (main_option_initfile) {
		com.initfile = g_strdup(main_option_initfile);
	}

	if (checkinitfile()) {
		fprintf(stderr,	_("Could not load or create settings file, exiting.\n"));
		exit(EXIT_FAILURE);
	}

	// After this point com.pref is populated
	siril_language_parser_init();
	if (com.pref.lang)
		language_init(com.pref.lang);

	if (main_option_directory) {
		gchar *cwd_forced;
		if (!g_path_is_absolute(main_option_directory))
			cwd_forced = g_build_filename(g_get_current_dir(), main_option_directory, NULL);
		else cwd_forced = g_strdup(main_option_directory);

		siril_change_dir(cwd_forced, NULL);
		if (com.pref.wd)
			g_free(com.pref.wd);
		com.pref.wd = g_strdup(cwd_forced);
		// if provided to the command line, make it persistent
		g_free(cwd_forced);
	}
	else {
		if (com.pref.wd && com.pref.wd[0] != '\0')
			siril_change_dir(com.pref.wd, NULL);
		else {
			// no other option
			siril_change_dir(siril_get_startup_dir(), NULL);
		}
	}

	init_num_procs();
	initialize_python_venv_in_thread();
	initialize_profiles_and_transforms(); // color management

#ifdef HAVE_LIBGIT2
	if (is_online()) {
		async_update_git_repositories();
	} else {
		siril_log_message(_("Siril started in offline mode. Will not attempt to update siril-scripts or siril-spcc-database...\n"));
	}
#else
	siril_log_message(_("Siril was compiled without libgit2 support. Remote repositories cannot be automatically fetched...\n"));
#endif

	if (com.headless) {
		if (main_option_script) {
			GInputStream *input_stream = NULL;

			if (g_strcmp0(main_option_script, "-") == 0) {
				input_stream = siril_input_stream_from_stdin();
			} else {
				GError *error = NULL;
				GFile *file = g_file_new_for_path(main_option_script);
				if (file)
					input_stream = (GInputStream *)g_file_read(file, NULL, &error);
				if (!input_stream) {
					if (error != NULL) {
						g_clear_error(&error);
						siril_log_message(_("File [%s] does not exist (from CWD, use absolute path?)\n"), main_option_script);
					}
					g_object_unref(file);
					exit(EXIT_FAILURE);
				}
				g_object_unref(file);
			}
#ifdef _WIN32
			ReconnectIO(1);
#endif
			if (execute_script(input_stream)) {
				exit(EXIT_FAILURE);
			}
		} else {
			pipe_start(main_option_rpipe_path, main_option_wpipe_path);
			read_pipe(main_option_rpipe_path);
		}
	}

	if (!com.headless) {
		/* Load GResource */
		com.resource = siril_resource_get_resource();
		/* Load preferred theme */
		load_prefered_theme(com.pref.gui.combo_theme);
		/* Load the css sheet for general style */
		load_css_style_sheet();
		/* Load UI files */
		load_ui_files();
		/* Passing GApplication to the control center */
		gtk_window_set_application(GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window"))), GTK_APPLICATION(application));
		/* Load state of the main windows (position and maximized) */
		gui_function(load_main_window_state, NULL);
#if defined(HAVE_LIBCURL)
		curl_global_init(CURL_GLOBAL_ALL);
		/* Check for update */
		if (is_online()) {
			if (com.pref.check_update) {
				siril_check_updates(FALSE);
			}
			siril_check_notifications(FALSE);
		}

#else
		gtk_widget_set_visible(lookup_widget("main_menu_updates"), FALSE);
		gtk_widget_set_visible(lookup_widget("frame24"), FALSE);
#endif
	}

	if (!com.headless) {
		gtk_builder_connect_signals(gui.builder, NULL);
		initialize_all_GUI(supported_files);
	}

	g_free(supported_files);
}

static void siril_app_open(GApplication *application, GFile **files, gint n_files, const gchar *hint) {
#if !defined(OS_OSX)
	g_application_activate(application);
#endif

	if (n_files > 0) {
		gchar *path = g_file_get_path(files[0]);
		const char *ext = get_filename_ext(path);
		if (ext && !strncmp(ext, "seq", 4)) {
			gchar *sequence_dir = g_path_get_dirname(path);
			if (!siril_change_dir(sequence_dir, NULL)) {
				if (check_seq()) {
					siril_log_message(_("No sequence `%s' found.\n"), path);
				} else {
					set_seq(path);
					if (!com.script) {
						populate_seqcombo(path);
						gui_function(set_GUI_CWD, NULL);
					}
				}
				g_free(sequence_dir);
			}
		} else {
			image_type type = get_type_from_filename(path);
			if (!main_option_directory && type != TYPEAVI && type != TYPESER
					&& type != TYPEUNDEF) {
				gchar *image_dir = g_path_get_dirname(path);
				siril_change_dir(image_dir, NULL);
				g_free(image_dir);
			}
			if (!com.script)
				gui_function(set_GUI_CWD, NULL);
			open_single_image(path);
		}
		g_free(path);
	}
}

#if defined(ENABLE_RELOCATABLE_RESOURCES) && defined(OS_OSX)
static void siril_macos_setenv(const char *progname) {
	/* helper to set environment variables for Siril to be relocatable.
	 * Due to the latest changes in Catalina it is not recommended
	 * to set it in the shell wrapper anymore.
	 */
	gchar resolved_path[PATH_MAX];

	if (realpath(progname, resolved_path)) {
		gchar tmp[PATH_MAX];
		gchar *exe_dir;           /* executable directory */
		gchar res_dir[PATH_MAX];  /* resources directory  */

		exe_dir = g_path_get_dirname(resolved_path);

		/* get canonical path to Foo.app/Contents/Resources directory */
		g_snprintf(tmp, sizeof(tmp), "%s/../Resources", exe_dir);
		struct stat sb;
		if (realpath(tmp, res_dir) && !stat(res_dir, &sb) && S_ISDIR(sb.st_mode)) {
			g_print("Siril is started as macOS application\n");
		}
		else {
			g_free(exe_dir);
			return;
		}

		/* store canonical path to resources directory in environment variable. */
		g_setenv("SIRIL_RELOCATED_RES_DIR", res_dir, TRUE);

    /* prepend PATH with our Contents/MacOS directory */
		gchar *path = g_try_malloc(PATH_MAX);
		if (path == NULL) {
			g_warning("Failed to allocate memory");
			exit(EXIT_FAILURE);
		}
		if (g_getenv("PATH"))
			g_snprintf(path, PATH_MAX, "%s:%s", exe_dir, g_getenv("PATH"));
		else
			g_snprintf(path, PATH_MAX, "%s", exe_dir);
  	g_setenv("PATH", path, TRUE);
		g_free(path);
		g_free(exe_dir);

		/* set XDG base directory specification variables */
		g_snprintf(tmp, sizeof(tmp), "%s/share", res_dir);
		g_setenv("XDG_DATA_DIRS", tmp, TRUE);
		if (g_getenv("HOME") != NULL) {
			g_snprintf (tmp, sizeof(tmp), "%s/Library/Caches/org.siril.Siril", g_getenv("HOME"));
			g_setenv ("XDG_CACHE_HOME", tmp, TRUE);
			g_snprintf(tmp, sizeof(tmp), "%s/Library/Application Support/org.siril.Siril", g_getenv("HOME"));
			g_setenv("XDG_CONFIG_HOME", tmp, TRUE);
			g_setenv("XDG_DATA_HOME", tmp, TRUE);
		}

		/* set GTK related environment variables */
		g_snprintf(tmp, sizeof(tmp), "%s/share/schemas", res_dir);
		g_setenv("GTK_PATH", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache", res_dir);
		g_setenv("GDK_PIXBUF_MODULE_FILE", tmp, TRUE);

		/* set fontconfig related variables */
		g_snprintf(tmp, sizeof(tmp), "%s/etc/fonts", res_dir);
		g_setenv("FONTCONFIG_PATH", tmp, TRUE);

		/* set curl related variables */
		g_snprintf(tmp, sizeof(tmp), "%s/etc/ca-certificates/cacert.pem", res_dir);
		g_setenv("CURL_CA_BUNDLE", tmp, TRUE);

		/* astropy does not create its director itself */
		g_snprintf(tmp, sizeof(tmp), "%s/astropy", g_getenv("XDG_CONFIG_HOME"));
		g_mkdir_with_parents(tmp, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH); // perm 755
	}
}
#endif


int main(int argc, char *argv[]) {
	GtkApplication *app;
	const gchar *dir;
	gint status;

#if defined(ENABLE_RELOCATABLE_RESOURCES) && defined(OS_OSX)
	// Remove macOS session identifier from command line arguments.
	// Code adopted from GIMP's app/main.c

	int new_argc = 0;
	for (int i = 0; i < argc; i++) {
		// Rewrite argv[] without "-psn_..." argument.
		if (!g_str_has_prefix(argv[i], "-psn_")) {
			argv[new_argc] = argv[i];
			new_argc++;
		}
	}
	if (argc > new_argc) {
		argv[new_argc] = NULL; // glib expects null-terminated array
		argc = new_argc;
	}

	siril_macos_setenv(argv[0]);
#elif _WIN32
	// suppression of annoying error boxes, hack from RawTherapee
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
#endif

#ifdef _WIN32
	// Turn off buffering for output to console on MSYS terminal on Windows
	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);
#endif

	initialize_siril_directories();

	dir = siril_get_locale_dir();
	setlocale(LC_ALL, "");
	bindtextdomain(PACKAGE, dir);
	bind_textdomain_codeset(PACKAGE, "UTF-8");
	textdomain(PACKAGE);

#if GLIB_CHECK_VERSION(2,74,0)
	app = gtk_application_new("org.siril.Siril", G_APPLICATION_DEFAULT_FLAGS | G_APPLICATION_HANDLES_OPEN | G_APPLICATION_NON_UNIQUE);
#else
	app = gtk_application_new("org.siril.Siril", G_APPLICATION_FLAGS_NONE | G_APPLICATION_HANDLES_OPEN | G_APPLICATION_NON_UNIQUE);
#endif

	g_signal_connect(app, "startup", G_CALLBACK(siril_app_startup), NULL);
	g_signal_connect(app, "activate", G_CALLBACK(siril_app_activate), NULL);
	g_signal_connect(app, "open", G_CALLBACK(siril_app_open), NULL);

	g_application_set_option_context_summary(G_APPLICATION(app), _("Siril - A free astronomical image processing software."));
	g_application_add_main_option_entries(G_APPLICATION(app), main_option);

	{	/* Setting the 'register-session' property is for macOS
		 * - to enable DnD via dock icon
		 * - to enable system menu "Quit"
		 */
		GValue value = G_VALUE_INIT;
		g_value_init(&value, G_TYPE_BOOLEAN);
		g_value_set_boolean(&value, TRUE);
		g_object_set_property(G_OBJECT(app), "register-session", &value);
		g_value_unset(&value);
	}

	status = g_application_run(G_APPLICATION(app), argc, argv);
	if (status) {
		gchar *help_msg;

		help_msg = g_strdup_printf(_("Run “%s --help” to see a full "
					"list of available command line "
					"options."), argv[0]);
		g_printerr("%s\n", help_msg);
		g_free(help_msg);
	}
	pipe_stop();		// close the pipes and their threads
	g_object_unref(app);
	return status;
}
