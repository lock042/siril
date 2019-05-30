/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#ifdef MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#include <io.h>
#include <fcntl.h>
#endif
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <locale.h>
#if (defined(__APPLE__) && defined(__MACH__))
#include <stdlib.h>
#include <libproc.h>
#endif
#include <getopt.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command.h"
#include "core/pipe.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "gui/callbacks.h"
#include "gui/gui.h"
#include "gui/script_menu.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "core/undo.h"
#include "io/single_image.h"
#include "algos/star_finder.h"
#include "algos/photometry.h"

#define GLADE_FILE "siril3.glade"

/* the global variables of the whole project */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
GtkBuilder *builder = NULL;	// get widget references anywhere

#ifdef _WIN32
/* origine du source: https://stackoverflow.com/questions/24171017/win32-console-application-that-can-open-windows */
int ReconnectIO(int OpenNewConsole)
{
    int    hConHandle;
    HANDLE lStdHandle;
    FILE  *fp;
    int    MadeConsole;

    MadeConsole=0;
    if(!AttachConsole(ATTACH_PARENT_PROCESS))
    {
        if(!OpenNewConsole)
            return 0;

        MadeConsole=1;
        if(!AllocConsole())
            return 0;  
    }

    // STDOUT to the console
    lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
    hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
    fp = _fdopen( hConHandle, "w" );
    *stdout = *fp;
    setvbuf( stdout, NULL, _IONBF, 0 );

     // STDIN to the console
    lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
    hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
    fp = _fdopen( hConHandle, "r" );
    *stdin = *fp;
    setvbuf( stdin, NULL, _IONBF, 0 );

    // STDERR to the console
    lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
    hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
    fp = _fdopen( hConHandle, "w" );
    *stderr = *fp;
    setvbuf( stderr, NULL, _IONBF, 0 );

    return MadeConsole;
}	
#endif

static void usage(const char *command) {
	printf("\nUsage:  %s [OPTIONS] [IMAGE_FILE_TO_OPEN]\n\n", command);
	puts("    -d, --directory CWD        changing the current working directory as the argument");
	puts("    -s, --script    SCRIPTFILE run the siril commands script in console mode");
	puts("    -i              INITFILE   load configuration from file name instead of the default configuration file");
	puts("    -p                         run in console mode with command and log stream through named pipes");
	puts("    -f, --format               print all supported image file formats (depending on installed libraries)");
	puts("    -v, --version              print program name and version and exit");
	puts("    -h, --help                 show this message");
}

static void signal_handled(int s) {
	// printf("Caught signal %d\n", s);
	gtk_main_quit();
}

struct option long_opts[] = {
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'},
	{"format", no_argument, 0, 'f'},
	{"directory", required_argument, 0, 'd'},
	{"script",    required_argument, 0, 's'},
	{0, 0, 0, 0}
};

int main(int argc, char *argv[]) {
	int i, c;
	extern char *optarg;
	extern int opterr;
	gchar *startup_cwd = NULL;
	gboolean forcecwd = FALSE;
	char *cwd_forced = NULL, *start_script = NULL;

	g_setenv ("LC_NUMERIC", "C", TRUE); // avoid possible bugs using french separator ","

	/* for translation */
#ifdef _WIN32
	setlocale(LC_ALL, "");

	gchar *localedir = g_build_filename(_getcwd(0, 0), "\\..\\share\\locale", NULL);
	gchar *localefilename = g_win32_locale_filename_from_utf8(localedir);
	bindtextdomain(PACKAGE, localefilename);
	bind_textdomain_codeset(PACKAGE, "UTF-8");
	g_free(localefilename);
	g_free(localedir);
#else
	bindtextdomain(PACKAGE, LOCALEDIR);
#endif
	textdomain(PACKAGE);

	opterr = 0;

	/* Caught signals */
	signal(SIGINT, signal_handled);

	while ((c = getopt_long(argc, argv, "i:phfvd:s:", long_opts, NULL)) != -1) {
		switch (c) {
			case 'i':
				com.initfile = g_strdup(optarg);
				break;
			case 'v':
				fprintf(stdout, "%s %s\n", PACKAGE, VERSION);
				exit(EXIT_SUCCESS);
				break;
			case 'f':
				list_format_available();
				exit(EXIT_SUCCESS);
				break;
			case 'd':
				cwd_forced = optarg;
				forcecwd = TRUE;
				break;
			case 's':
			case 'p':
				com.script = TRUE;
				com.headless = TRUE;
				/* need to force cwd to the current dir if no option -d */
				if (!forcecwd) {
					cwd_forced = g_get_current_dir();
					forcecwd = TRUE;
				}
				if (c == 's')
					start_script = optarg;
				break;
			default:
				fprintf(stderr, _("unknown command line parameter '%c'\n"), argv[argc - 1][1]);
				/* no break */
			case 'h':
				usage(argv[0]);
				exit(EXIT_SUCCESS);
		}
	}

	com.cvport = RED_VPORT;
	com.show_excluded = TRUE;
	com.selected_star = -1;
	com.stacking_zone_focus = -1;
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.uniq = NULL;
	com.color = NORMAL_COLOR;
	for (i = 0; i < MAXVPORT; i++)
		com.buf_is_dirty[i] = TRUE;
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));

	/* initialize the com struct and zoom level */
	com.sliders = MINMAX;
	com.zoom_value = ZOOM_DEFAULT;
	com.stack.memory_percent = 0.9;
	com.app_path = NULL;

	siril_log_color_message(_("Welcome to %s v%s\n"), "bold", PACKAGE, VERSION);

	/***************
	 *  initialization of some parameters that need to be done before
	 * checkinitfile
	 ***************/
	/* initialize converters (utilities used for different image types importing) */
	gchar *supported_files = initialize_converters();
	/* initialize photometric variables */
	initialize_photometric_param();
	/* initialize peaker variables */
	init_peaker_default();
	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);

	/* set default CWD, and load init file
	 * checkinitfile will load all saved parameters
	 * */
	com.wd = siril_get_startup_dir();
	/* load init file */
	startup_cwd = g_get_current_dir();
	if (checkinitfile()) {
		fprintf(stderr,	_("Could not load or create settings file, exiting.\n"));
		exit(1);
	}

	if (!com.headless) {
		gtk_init(&argc, &argv);
		com.siril_mode = MODE_NO_GUI;
		init_gui(MODE_PLANETARY, startup_cwd);

		show_supported_files(supported_files);
	} else {
		com.siril_mode = MODE_NO_GUI;	// TODO: transition headless mode to that
	}

	changedir(com.wd, NULL);

	if (!com.headless) {
		gtk_builder_connect_signals (builder, NULL);
		initialize_all_GUI(supported_files);
	}
	g_free(supported_files);

	changedir(com.wd, NULL);

	/* Get CPU number and set the number of threads */
	siril_log_message(_("Parallel processing %s: using %d logical processor(s).\n"),
#ifdef _OPENMP
			_("enabled"), com.max_thread = omp_get_num_procs()
#else
			_("disabled"), com.max_thread = 1
#endif
			);

	/* start Siril */
	if (argv[optind] != NULL) {
		if (startup_cwd) {
			changedir(startup_cwd, NULL);
		}
		open_single_image(argv[optind]);
		if (!forcecwd) {
			gchar *newpath = g_path_get_dirname(argv[optind]);
			changedir(newpath, NULL);
			g_free(newpath);
		}
	}
	g_free(startup_cwd);

	if (forcecwd && cwd_forced) {
		changedir(cwd_forced, NULL);
	}

	if (!com.script) {
		set_GUI_CWD();
	}

	if (com.headless) {
		if (start_script) {
			FILE* fp = g_fopen(start_script, "r");
			if (fp == NULL) {
				siril_log_message(_("File [%s] does not exist\n"), start_script);
				exit(1);
			}
#ifdef _WIN32			
			ReconnectIO(1);
#endif
			execute_script(fp);
		}
		else {
			pipe_start();
			read_pipe(NULL);
		}
	}
	else {
		do {
			gtk_main();

			process_close(0);
			// wait for clean-up, apparently it never gets there
			while (gtk_events_pending())
				gtk_main_iteration_do(FALSE);
			//uninit_gui();

			init_gui(com.siril_mode, startup_cwd);	// sets new mode
		} while (com.requested_mode != MODE_NO_GUI);
	}

	/* quit Siril */
	close_sequence(FALSE);	// closing a sequence if loaded
	close_single_image();	// close the previous image and free resources
	pipe_stop();		// close the pipes and their threads
	g_free(com.app_path);
	return 0;
}
