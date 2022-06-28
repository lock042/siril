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
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"

#include <unistd.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#include <stdio.h>
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

#ifndef _WIN32
static int exec_prog(const char **argv)
{
	pid_t my_pid;
	int status;

	if (0 == (my_pid = fork())) {
		if (-1 == execve(argv[0], (char **)argv , NULL)) {
			perror("child process execve failed [%m]");
			return 0;
		}
	} else {
		while (0 == waitpid(my_pid , &status , WNOHANG)) {
			sleep(1);	// Wait for starnet++ to finish before attempting to process the output
		}

		if (1 != WIFEXITED(status) || 0 != WEXITSTATUS(status)) {
			siril_log_color_message(_("Error: external command %s failed...\n"), "red", argv[0]);
			return 0;
		}
	}

	return 0;
}
#else
static int exec_prog_win32(const char **argv) {

	if (-1 == _spawnve(_P_WAIT, argv[0], argv , NULL)) {
		perror("child process _spawnve failed [%m]");
		return -1;
	}
	return 0;
}
#endif

/* Starnet++v2 star removal routine */

gpointer do_starnet() {
	int retval;

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

	siril_log_color_message(_("Starnet++: running. Please wait...\n"), "green");

	// Check starnet directory is defined
	if (g_access(com.pref.starnet_dir, R_OK)) {
		siril_log_color_message(_("Incorrect permissions on the Starnet++ directory: %s\nEnsure it is correctly set in Preferences / Miscellaneous.\n"), "red", com.pref.starnet_dir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Set up paths and filenames
	imagenoext = remove_ext_from_filename(com.uniq->filename);

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
	my_argv[0] = starnetcommand;
	my_argv[1] = temptif;
	my_argv[2] = starlesstif;
	fits fit = { 0 };

	// Store current working directory
	currentdir = g_get_current_dir();

	// Change to starnet directory

	retval = g_chdir(com.pref.starnet_dir);
	if (retval) {
		siril_log_color_message(_("Error: unable to change to Starnet++ directory.\nEnsure it is set in Preferences / Miscellaneous...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Save current image as working 16-bit TIFF
	if (savetif(temptif, &gfit, 16)) {
		siril_log_color_message(_("Error: unable to save working TIFF of original image...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Check for starnet executables (pre-v2.0.2 or v2.0.2+)
	if (g_file_test("starnet++", G_FILE_TEST_EXISTS)) {
		snprintf(starnetcommand,15, "starnet++");
		snprintf(starnetcommand,15, "starnet++");
	} else if ((gfit.naxes[2] == 3) && (g_file_test("rgb_starnet++", G_FILE_TEST_EXISTS))) {
		snprintf(starnetcommand,15, "rgb_starnet++");
	} else if ((gfit.naxes[2] == 1 ) && (g_file_test("mono_starnet++", G_FILE_TEST_EXISTS))) {
		snprintf(starnetcommand, 15, "mono_starnet++");
	} else {
		siril_log_color_message(_("No valid executable found in the Starnet++ directory\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// *** Call starnet++ *** //
#ifdef _WIN32
	retval = exec_prog_win32(my_argv);
#else
	retval = exec_prog(my_argv);
#endif
	if (retval) {
		siril_log_color_message(_("Error: Starnet++ did not execute correctly...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Read the starless tiff
	retval = readtif(starlesstif, &gfit, FALSE);

	// Chdir back to the Siril working directory, we don't need to be in the starnet
	// directory any more
	if (g_chdir(currentdir)) {
		siril_log_color_message(_("Error: unable to change to Siril working directory...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Save starless image as fits
	if (savefits(starlessfit, &gfit)) {
		siril_log_color_message(_("Error: unable to save starless image as FITS...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Subtract starless from original to create starmask
	// Load the original image TIFF back into gfit
	retval = readtif(temptif, &gfit, FALSE);
	// Load the starless version into fit
	if (readfits(starlessfit, &fit, NULL, !com.pref.force_to_16bit)) {
		siril_log_color_message(_("Error: unable to load starless image for starmask generation...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}
	// Subtract starless from original
	retval = imoper(&gfit, &fit, OPER_SUB, !com.pref.force_to_16bit);
	clearfits(&fit);

	// Save starmask as fits
	if (savefits(starmaskfit, &gfit)) {
		siril_log_color_message(_("Error: unable to save starmask image as FITS...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Load starless - the final result we want to show is the starless version
	if (readfits(starlessfit, &gfit, NULL, !com.pref.force_to_16bit)) {
		siril_log_color_message(_("Error: unable to save starless image as FITS...\n"), "red");
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	// Remove working files, they are no longer required
	if (remove(starlesstif)) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
		siril_add_idle(end_generic, NULL);
		free(imagenoext);
		free(currentdir);
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
		notify_gfit_modified();
		return NULL;
	}

	if (remove(temptif)) {
		siril_log_color_message(_("Error: unable to remove working file...\n"), "red");
	siril_add_idle(end_generic, NULL);
	free(imagenoext);
	free(currentdir);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	notify_gfit_modified();
	return NULL;
	}

	siril_add_idle(end_generic, NULL);
	free(imagenoext);
	free(currentdir);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	notify_gfit_modified();

	return NULL;
}
/*
void on_starnet_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(gui.builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_starnet_apply_clicked(GtkButton *button, gpointer user_data) {

	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_scnr")));
	GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(gui.builder, "preserve_light"));
	gboolean preserve = gtk_toggle_button_get_active(light_button);
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_scnr")));

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	undo_save_state(&gfit, _("SCNR (type=%d, amount=%0.2lf, preserve=%s)"),
			type, amount, preserve ? "true" : "false");

	args->fit = &gfit;
	args->type = type;
	args->amount = amount;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	start_in_new_thread(scnr, args);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("starnet_dialog");
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));
	GtkSpinButton *spinButton = GTK_SPIN_BUTTON(lookup_widget("spin_scnr"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(spinButton), type > 1);
}

*/
