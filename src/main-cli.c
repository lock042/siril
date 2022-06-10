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

#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <locale.h>
#include <unistd.h>
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

#include "git-version.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_actions.h"
#include "core/initfile.h"
#include "core/command_line_processor.h"
#include "core/command.h"
#include "core/pipe.h"
#include "core/signals.h"
#include "core/siril_app_dirs.h"
#include "core/OS_utils.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"

/* the global variables of the whole project */
cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image
const gchar *startup_cwd = NULL;
gboolean forcecwd = FALSE;

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

static GOptionEntry main_option[] = {
	{ "directory", 'd', 0, G_OPTION_ARG_FILENAME, &main_option_directory, N_("changing the current working directory as the argument"), NULL },
	{ "script", 's', 0, G_OPTION_ARG_FILENAME, &main_option_script, N_("run the siril commands script in console mode. If argument is equal to \"-\", then siril will read stdin input"), NULL },
	{ "initfile", 'i', 0, G_OPTION_ARG_FILENAME, &main_option_initfile, N_("load configuration from file name instead of the default configuration file"), NULL },
	{ "pipe", 'p', 0, G_OPTION_ARG_NONE, &main_option_pipe, N_("run in console mode with command and log stream through named pipes"), NULL },
	{ "inpipe", 'r', 0, G_OPTION_ARG_FILENAME, &main_option_rpipe_path, N_("specify the path for the read pipe"), NULL },
	{ "outpipe", 'w', 0, G_OPTION_ARG_FILENAME, &main_option_wpipe_path, N_("specify the path for the write pipe"), NULL },
	{ "format", 'f', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_list_of_formats_and_exit, N_("print all supported image file formats (depending on installed libraries)" ), NULL },
	{ "version", 'v', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_version_and_exit, N_("print the application’s version"), NULL},
	{ "copyright", 'c', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_copyright_and_exit, N_("print the copyright"), NULL},
	{ NULL },
};

static void global_initialization() {
	com.selected_star = -1;
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.qphot = NULL;
	com.tilt = NULL;
	com.uniq = NULL;
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));

	gui.cvport = RED_VPORT;
	gui.show_excluded = TRUE;
	//gui.color = NORMAL_COLOR;
	gui.sliders = MINMAX;
	gui.zoom_value = ZOOM_DEFAULT;

	initialize_default_preferences();
}

static void init_num_procs() {
	/* Get CPU number and set the number of threads */
#ifdef _OPENMP
	int num_proc = (int) g_get_num_processors();
	int omp_num_proc = omp_get_num_procs();
	if (num_proc != omp_num_proc) {
		siril_log_message(_("Questionable parallel processing efficiency - openmp reports %d %s. "
					"Possibly broken opencv/openblas installation.\n"),	omp_num_proc,
				ngettext("processor", "processors", omp_num_proc));
	}
	omp_set_nested(1);
	int supports_nesting = omp_get_nested() && omp_get_max_active_levels() > 1;
	siril_log_message(
			_("Parallel processing %s: using %d logical %s%s.\n"),
			_("enabled"), com.max_thread = num_proc,
			ngettext("processor", "processors", num_proc),
			supports_nesting ? "" : _(", nesting not supported"));
#else
	com.max_thread = 1;
	siril_log_message(_("Parallel processing disabled: using 1 logical processor.\n"));
#endif
}

static void siril_app_activate(GApplication *application) {
	gchar *cwd_forced = NULL;

	/*
	 * Force C locale for numbers to avoid "," being used as decimal separator.
	 * Called here and not in main() because setlocale(LC_ALL, "") is called as
	 * part of g_application_run().
	 */
	setlocale(LC_NUMERIC, "C");

	memset(&com, 0, sizeof(struct cominf));	// needed? doesn't hurt
	com.initfile = NULL;

	com.script = TRUE;
	com.headless = TRUE;
	/* need to force cwd to the current dir if no option -d */
	if (!forcecwd) {
		cwd_forced = g_strdup(g_get_current_dir());
		forcecwd = TRUE;
	}

	global_initialization();

	/* initialize peaker variables */
	init_peaker_default();
	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);

	siril_log_color_message(_("Welcome to %s v%s\n"), "bold", PACKAGE, VERSION);

	/* initialize converters (utilities used for different image types importing) */
	gchar *supported_files = initialize_converters();
	startup_cwd = g_get_current_dir();

	if (main_option_initfile) {
		com.initfile = g_strdup(main_option_initfile);
	}

	if (checkinitfile()) {
		fprintf(stderr,	_("Could not load or create settings file, exiting.\n"));
		exit(EXIT_FAILURE);
	}

	if (main_option_directory) {
		if (!g_path_is_absolute(main_option_directory)) {
			cwd_forced = g_build_filename(g_get_current_dir(), main_option_directory, NULL);
		} else {
			cwd_forced = g_strdup(main_option_directory);
		}
		forcecwd = TRUE;
	}

	if (forcecwd && cwd_forced) {
		siril_change_dir(cwd_forced, NULL);
		g_free(cwd_forced);
	}

	init_num_procs();

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

	if (siril_change_dir(com.wd, NULL)) {
		com.wd = g_strdup(siril_get_startup_dir());
		siril_change_dir(com.wd, NULL);
	}

	g_free(supported_files);
}

static void siril_app_open(GApplication *application, GFile **files, gint n_files,
		const gchar *hint) {

	g_application_activate(application);

	if (n_files > 0) {
		gchar *path = g_file_get_path (files[0]);
		const char *ext = get_filename_ext(path);
		if (ext && !strncmp(ext, "seq", 4)) {
			gchar *sequence_dir = g_path_get_dirname(path);
			if (!siril_change_dir(sequence_dir, NULL)) {
				if (check_seq(FALSE)) {
					siril_log_message(_("No sequence `%s' found.\n"), path);
				} else {
					set_seq(path);
					if (!com.script)
						set_GUI_CWD();
				}
				g_free(sequence_dir);
			}
		} else {
			image_type type = get_type_from_filename(path);
			if (!forcecwd && type != TYPEAVI && type != TYPESER && type != TYPEUNDEF) {
				gchar *image_dir = g_path_get_dirname(path);
				siril_change_dir(image_dir, NULL);
				g_free(image_dir);
			} else if (startup_cwd) {
				siril_change_dir(startup_cwd, NULL);
			}
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
		gchar *path;
		gchar tmp[PATH_MAX];
		gchar *app_dir;
		gchar lib_dir[PATH_MAX];
		size_t path_len;
		struct stat sb;
		app_dir = g_path_get_dirname(resolved_path);

		g_snprintf(tmp, sizeof(tmp), "%s/../Resources", app_dir);
		if (realpath(tmp, lib_dir) && !stat(lib_dir, &sb) && S_ISDIR(sb.st_mode))
			g_print("SiriL is started as MacOS application\n");
		else
			return;

		/* we define the relocated resources path */
		g_setenv("SIRIL_RELOCATED_RES_DIR", tmp, TRUE);

		path_len = strlen(g_getenv("PATH") ? g_getenv("PATH") : "")
			+ strlen(app_dir) + 2;
		path = g_try_malloc(path_len);
		if (path == NULL) {
			g_warning("Failed to allocate memory");
			exit(EXIT_FAILURE);
		}
		if (g_getenv("PATH"))
			g_snprintf(path, path_len, "%s:%s", app_dir, g_getenv("PATH"));
		else
			g_snprintf(path, path_len, "%s", app_dir);
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
		if (g_getenv("HOME") != NULL) {
			g_snprintf(tmp, sizeof(tmp), "%s/Library/Application Support", g_getenv("HOME"));
			g_setenv("XDG_CONFIG_HOME", tmp, TRUE);
			g_snprintf (tmp, sizeof(tmp), "%s/Library/Application Support/SiriL/1.00/cache",
					g_getenv("HOME"));
			g_setenv ("XDG_CACHE_HOME", tmp, TRUE);

		}
	}
}
#endif


int main(int argc, char *argv[]) {
	GApplication *app;
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

	app = g_application_new("org.free_astro.siril", G_APPLICATION_HANDLES_OPEN | G_APPLICATION_NON_UNIQUE);

	g_signal_connect(app, "activate", G_CALLBACK(siril_app_activate), NULL);
	g_signal_connect(app, "open", G_CALLBACK(siril_app_open), NULL);

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

	pipe_stop();		// close the pipes and their threads
	g_object_unref(app);
	return status;
}
