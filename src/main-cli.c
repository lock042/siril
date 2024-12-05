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

#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gsl/gsl_errno.h>
#include <gtk/gtk.h>
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

#include "git-version.h"
#include "core/siril.h"
#include "core/icc_profile.h"
#include "core/proto.h"
#include "core/siril_actions.h"
#include "core/initfile.h"
#include "core/command_line_processor.h"
#include "core/command.h"
#include "core/pipe.h"
#include "core/signals.h"
#include "core/siril_app_dirs.h"
#include "core/siril_language.h"
#include "core/siril_log.h"
#include "core/siril_networking.h"
#include "core/OS_utils.h"
#include "algos/siril_random.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "io/siril_pythonmodule.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"

/* the global variables of the whole project */
cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
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

static void global_initialization() {
	gsl_set_error_handler_off();
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.tilt = NULL;
	com.uniq = NULL;
	com.child_is_running = EXT_NONE;
	com.kernel = NULL;
	com.kernelsize = 0;
	com.kernelchannels = 0;
	com.childpid = 0;
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));
	initialize_default_settings();	// com.pref
#ifdef HAVE_FFTW3F_OMP
	fftwf_init_threads(); // Should really only be called once so do it at startup
#endif
#ifdef _OPENMP
	omp_set_max_active_levels(2);
#endif

}

static void siril_app_activate(GApplication *application) {
	/*
	 * Force C locale for numbers to avoid "," being used as decimal separator.
	 * Called here and not in main() because setlocale(LC_ALL, "") is called as
	 * part of g_application_run().
	 */
	setlocale(LC_NUMERIC, "C");

	com.script = TRUE;
	com.headless = TRUE;
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

#if defined(HAVE_LIBCURL)
	curl_global_init(CURL_GLOBAL_ALL);
#endif

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

	g_free(supported_files);
}

#if defined(ENABLE_RELOCATABLE_RESOURCES) && defined(OS_OSX)
static void siril_macos_setenv(const char *progname) {
	/* helper to set environment variables for Siril to be relocatable.
	 * Due to the latest changes in Catalina it is not recommended
	 * to set it in the shell wrapper anymore.
	 */
	gchar resolved_path[PATH_MAX];

	if (realpath(progname, resolved_path)) {
		gchar *path;
		gchar tmp[PATH_MAX];
		gchar *app_dir;
		gchar lib_dir[PATH_MAX];
		size_t path_len;
		struct stat sb;
		app_dir = g_path_get_dirname(resolved_path);

		g_snprintf(tmp, sizeof(tmp), "%s/../Resources", app_dir);
		if (realpath(tmp, lib_dir) && !stat(lib_dir, &sb) && S_ISDIR(sb.st_mode))
			g_print("Siril is started as MacOS application\n");
		else
			return;

		/* we define the relocated resources path */
		g_setenv("SIRIL_RELOCATED_RES_DIR", tmp, TRUE);

		g_snprintf(tmp, sizeof(tmp), "%s/../Resources/bin", app_dir);
		path_len = strlen(g_getenv("PATH") ? g_getenv("PATH") : "")
			+ strlen(tmp) + 2;
		path = g_try_malloc(path_len);
		if (path == NULL) {
			g_warning("Failed to allocate memory");
			exit(EXIT_FAILURE);
		}
		if (g_getenv("PATH"))
			g_snprintf(path, path_len, "%s:%s", tmp, g_getenv("PATH"));
		else
			g_snprintf(path, path_len, "%s", tmp);
		/* the relocated path is storred in this env. variable in order to be reused if needed */
		g_free(app_dir);
		g_setenv("PATH", path, TRUE);
		g_free(path);
		g_snprintf(tmp, sizeof(tmp), "%s/share", lib_dir);
		g_setenv("XDG_DATA_DIRS", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/share/schemas", lib_dir);
		g_setenv("GSETTINGS_SCHEMA_DIR", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/lib/gtk-3.0/3.0.0", lib_dir);
		g_setenv("GTK_PATH", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache", lib_dir);
		g_setenv("GDK_PIXBUF_MODULE_FILE", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/lib/gdk-pixbuf-2.0/2.10.0/loaders", lib_dir);
		g_setenv("GDK_PIXBUF_MODULE_DIR", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/etc/fonts", lib_dir);
		g_setenv("FONTCONFIG_PATH", tmp, TRUE);
		g_snprintf(tmp, sizeof(tmp), "%s/etc/ca-certificates/cacert.pem", lib_dir);
		g_setenv("CURL_CA_BUNDLE", tmp, TRUE);
		if (g_getenv("HOME") != NULL) {
			g_snprintf(tmp, sizeof(tmp), "%s/Library/Application Support/org.siril.Siril", g_getenv("HOME"));
			g_setenv("XDG_CONFIG_HOME", tmp, TRUE);
			g_snprintf (tmp, sizeof(tmp), "%s/Library/Caches/org.siril.Siril",
					g_getenv("HOME"));
			g_setenv ("XDG_CACHE_HOME", tmp, TRUE);
		}
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

	g_signal_connect(app, "activate", G_CALLBACK(siril_app_activate), NULL);
	//g_signal_connect(app, "open", G_CALLBACK(siril_app_open), NULL);

	g_application_set_option_context_summary(G_APPLICATION(app), _("Siril - A free astronomical image processing software."));
	g_application_add_main_option_entries(G_APPLICATION(app), main_option);

	status = g_application_run(G_APPLICATION(app), argc, argv);
	if (status) {
		gchar *help_msg;

		help_msg = g_strdup_printf(_("Run “%s --help” to see a full "
					"list of available command line "
					"options."), argv[0]);
		g_printerr("%s\n", help_msg);
		g_free(help_msg);
	}

	cmsUnregisterPlugins(); // unregister any lcms2 plugins

	pipe_stop();		// close the pipes and their threads
	g_object_unref(app);
	return status;
}
