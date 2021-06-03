/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "io/FITS_symlink.h"

#define EXIT_TOKEN ":EXIT:"

static GFileMonitor *dirmon = NULL;
static gboolean do_links = FALSE;
static GThread *live_stacker_thread = NULL;
static GAsyncQueue *new_files_queue = NULL;

static gpointer live_stacker(gpointer arg);

void display(const char *str) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	gtk_label_set_text(label, str);
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide((GtkWidget *)user_data);
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gboolean is_fullscreen;

	window = GTK_WINDOW(lookup_widget("control_window"));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);

	if (dirmon) {
		g_object_unref(dirmon);
		dirmon = NULL;
	}
	if (new_files_queue) {
		g_async_queue_push(new_files_queue, EXIT_TOKEN);
		g_async_queue_unref(new_files_queue);
		new_files_queue = NULL;
	}
	if (live_stacker_thread) {
		g_thread_join(live_stacker_thread);
		live_stacker_thread = NULL;
	}
}

static void file_changed(GFileMonitor *monitor, GFile *file, GFile *other,
		GFileMonitorEvent evtype, gpointer user_data) {
	if (evtype == G_FILE_MONITOR_EVENT_CREATED) {
		gchar *filename = g_file_get_basename(file);
		siril_debug_print("File %s created\n", filename);

		image_type type;
		stat_file(filename, &type, NULL);
		if (type != TYPEFITS) {
			gchar *str = g_strdup_printf(_("File not supported for live stacking: %s"), filename);
			display(str);
			g_free(str);
		} else {
			if (strncmp(filename, "live_stack", 10) &&
					strncmp(filename, "r_live_stack", 12) &&
					strncmp(filename, "result_live_stack", 17)) {
				/* TODO: start stacking */
				g_async_queue_push(new_files_queue, filename);
			}
		}
		//g_free(filename);
	}
}

void on_livestacking_start() {
	display(_("Starting live stacking"));
	/* start monitoring CWD */
	GFile *cwd = g_file_new_for_path(com.wd);
	GError *err = NULL;
	dirmon = g_file_monitor_directory(cwd, G_FILE_MONITOR_WATCH_MOVES, NULL, &err);
	g_object_unref(cwd);
	if (err) {
		fprintf(stderr, "Unable to monitor CWD (%s): %s\n", com.wd, err->message);
		g_error_free(err);
		return;
	}

	g_signal_connect(G_OBJECT(dirmon), "changed", G_CALLBACK(file_changed), NULL);
	display(_("Live stacking waiting for files"));

	do_links = test_if_symlink_is_ok();

	new_files_queue = g_async_queue_new();
	live_stacker_thread = g_thread_new("live stacker", live_stacker, NULL);
}


static gpointer live_stacker(gpointer arg) {
	g_async_queue_ref(new_files_queue);
	int index = 1;
	do {
		gchar *filename = g_async_queue_pop(new_files_queue); // blocking
		if (!strcmp(filename, EXIT_TOKEN)) {
			siril_debug_print("Exiting thread\n");
			break;
		}
		siril_debug_print("Adding file to input sequence\n");
		gchar *target = g_strdup_printf("live_stack_%05d%s", index, com.pref.ext);
		if (symlink_uniq_file(filename, target, do_links)) {
			break;
		}
		g_free(filename);

		/* TODO: start processing */
		/* Create the sequence */
		/* Register the sequence */
		/* Stack the sequence */
		/* Update display */

		index++;
	} while (1);

	g_async_queue_unref(new_files_queue);
	return NULL;
}

