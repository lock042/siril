/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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
#include "gtkmacintegration/gtkosxapplication.h"
#endif
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command.h"
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

#define GLADE_FILE "siril3.glade"
#define PLANETARY_GLADE_FILE "siril_planetary3.glade"

/* the global variables of the whole project */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
GtkBuilder *builder; // get widget references anywhere

#ifdef MAC_INTEGRATION

static gboolean osx_open_file(GtkosxApplication *osx_app, gchar *path, gpointer data){
	if (path != NULL) {
		open_single_image(path);
		return FALSE;
	}
	return TRUE;
}

static void gui_add_osx_to_app_menu(GtkosxApplication *osx_app, const gchar *item_name, gint index) {
	GtkWidget *item;

	item = lookup_widget(item_name);
	if (GTK_IS_MENU_ITEM(item))
		gtkosx_application_insert_app_menu_item(osx_app, GTK_WIDGET(item), index);
}

static void set_osx_integration(GtkosxApplication *osx_app, gchar *siril_path) {
	GtkWidget *menubar = lookup_widget("menubar1");
	GtkWidget *file_quit_menu_item = lookup_widget("exit");
	GtkWidget *help_menu = lookup_widget("help1");
	GtkWidget *window_menu = lookup_widget("menuitemWindows");
	GtkWidget *sep;
	GdkPixbuf *icon;
	gchar *icon_path;
	
	g_signal_connect(osx_app, "NSApplicationOpenFile", G_CALLBACK(osx_open_file), NULL);

	gtk_widget_hide(menubar);

	gtkosx_application_set_menu_bar(osx_app, GTK_MENU_SHELL(menubar));
	gtkosx_application_set_window_menu(osx_app, GTK_MENU_ITEM(window_menu));

	gui_add_osx_to_app_menu(osx_app, "help_item1", 0);
	gui_add_osx_to_app_menu(osx_app, "help_update", 1);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 2);
	gui_add_osx_to_app_menu(osx_app, "settings", 3);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 4);

	gtk_widget_hide(file_quit_menu_item);
	gtk_widget_hide(help_menu);
	
	icon_path = g_build_filename(siril_path, "pixmaps/siril_1.svg", NULL);
	icon = gdk_pixbuf_new_from_file(icon_path, NULL);
	gtkosx_application_set_dock_icon_pixbuf(osx_app, icon);
		
	gtkosx_application_ready(osx_app);
	g_free(icon_path);
}

#endif

void usage(const char *command) {
	printf("\nUsage:  %s [OPTIONS] [IMAGE_FILE_TO_OPEN]\n\n", command);
	puts("    -d, --directory CWD        changing the current working directory as the argument");
	puts("    -s, --script    SCRIPTFILE run the siril commands script in console mode");
	puts("    -i                         load configuration from file name instead of the default configuration file");
	puts("    -f, --format               print all supported image file formats (depending on installed libraries)");
	puts("    -v, --version              print program name and version and exit");
	puts("    -h, --help                 show this message");
}

void signal_handled(int s) {
	// printf("Caught signal %d\n", s);
	undo_flush();
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
	int i;
	extern char *optarg;
	extern int opterr;
	gboolean forcecwd = FALSE;
	char *cwd_forced = NULL, *start_script = NULL;

	g_setenv ("LC_NUMERIC", "C", TRUE); // avoid possible bugs using french separator ","

	/* for translation */
#ifdef _WIN32
	setlocale(LC_ALL, "");

	gchar *localedir = g_build_filename(_getcwd(0, 0), "\\..\\share\\locale", NULL);
	bindtextdomain(PACKAGE, g_win32_locale_filename_from_utf8(localedir));
	bind_textdomain_codeset(PACKAGE, "UTF-8");
	g_free(localedir);
#else
	bindtextdomain(PACKAGE, LOCALEDIR);
#endif
	textdomain(PACKAGE);

	opterr = 0;

	/* Caught signals */
	signal(SIGINT, signal_handled);

	while (1) {
		signed char c = getopt(argc, argv, "i:hfvd:s:");
		if (c == '?') {
			for (i = 1; i < argc; i++) {
				if (argv[i][1] == '-') {
					if (!strcmp(argv[i], "--version"))
						c = 'v';
					else if (!strcmp(argv[i], "--help"))
						c = 'h';
					else if (!strcmp(argv[i], "--format"))
						c = 'f';
					else if (!strcmp(argv[i], "--directory"))
						c = 'd';
					else if (!strcmp(argv[i], "--script"))
						c = 's';
					else {
						usage(argv[0]);
						exit(EXIT_FAILURE);
					}
				}
			}
		}

		if (c == -1)
			break;
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
				start_script = optarg;
				com.script = TRUE;
				com.headless = TRUE;
				/* need to force cwd to the current dir if no option -d */
				if (!forcecwd) {
					cwd_forced = g_get_current_dir();
					forcecwd = TRUE;
				}
				break;
			default:
				fprintf(stderr, _("unknown command line parameter '%c'\n"), argv[argc - 1][1]);
				/* no break */
			case 'h':
				usage(argv[0]);
				exit(EXIT_SUCCESS);
		}
	}

	/* initializing internal structures with widgets (drawing areas) */
	com.cvport = RED_VPORT;
	com.show_excluded = TRUE;
	com.selected_star = -1;
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.uniq = NULL;
	com.grad = NULL;
	com.grad_boxes_drawn = TRUE;
	com.color = NORMAL_COLOR;
	for (i=0; i<MAXVPORT; i++)
		com.buf_is_dirty[i] = TRUE;
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));

	/* initialize the com struct and zoom level */
	com.sliders = MINMAX;
	com.zoom_value = ZOOM_DEFAULT;

	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);

	/* set default CWD */
	com.wd = siril_get_startup_dir();
	com.startup_dir = g_get_current_dir();

	/* load init file */
	if (checkinitfile()) {
		siril_log_message(_("Could not load or create settings file, exiting.\n"));
		exit(1);
	}

	if (!com.headless) {
		gtk_init(&argc, &argv);
		enum _siril_mode mode = com.siril_mode;
		com.siril_mode = MODE_NO_GUI;
		init_gui(mode);	// sets new mode
	} else {
		com.siril_mode = MODE_NO_GUI;	// TODO: transition headless mode to that
	}

	siril_log_color_message(_("Welcome to %s v%s\n"), "bold", PACKAGE, VERSION);

	/* initialize converters (supported file type) */
	initialize_converters();

	/* Get CPU number and set the number of threads */
	siril_log_message(_("Parallel processing %s: Using %d logical processor(s).\n"),
#ifdef _OPENMP
			_("enabled"), com.max_thread = omp_get_num_procs()
#else
			_("disabled"), com.max_thread = 1
#endif
			);

	if (com.headless) {
		init_peaker_default();
	}

	/* handling OS-X integration */
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	if (!com.headless) {
		set_osx_integration(osx_app, siril_path);
	}
#endif //MAC_INTEGRATION

	/* start Siril */
	if (argv[optind] != NULL) {
		changedir(com.startup_dir, NULL);
		open_single_image(argv[optind]);
		if (!forcecwd) {
			gchar *newpath = g_path_get_dirname(argv[optind]);
			changedir(newpath, NULL);
			g_free(newpath);
		}
	}

	if (forcecwd && cwd_forced) {
		changedir(cwd_forced, NULL);
	}

	if (com.headless) {
		FILE* fp = g_fopen(start_script, "r");
		if (fp == NULL) {
			siril_log_message(_("File [%s] does not exist\n"), start_script);
			exit(1);
		}
		execute_script(fp);
	}
	else {
		do {
			gtk_main();

			process_close(0);
			// wait for clean-up, apparently it never gets there
			while (gtk_events_pending())
				gtk_main_iteration_do(FALSE);
			//uninit_gui();

			init_gui(com.requested_mode);
		} while (com.requested_mode != MODE_NO_GUI);
	}

	/* quit Siril */
	close_sequence(FALSE);
	undo_flush();
#ifdef MAC_INTEGRATION
	g_object_unref(osx_app);
#endif //MAC_INTEGRATION
	return 0;
}
