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

#include <gtk/gtk.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "core/proto.h"
#include "core/pipe.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/command.h"

#include "progress_and_log.h"

/*************************** P R O G R E S S    B A R ***************************/

static void progress_bar_set_text(const char *text) {
	static GtkProgressBar *pbar = NULL;
	if (pbar == NULL)
		pbar = GTK_PROGRESS_BAR(gtk_builder_get_object(gui.builder, "progressbar1"));
	/* It will not happen that text is NULL here, because it's
	 * catched by set_progress_bar_data() */
	if (!text || text[0] == '\0')
		text = _("Ready.");
	gtk_progress_bar_set_text(pbar, text);
}

/* http://developer.gnome.org/gtk3/stable/GtkProgressBar.html */
static void progress_bar_set_percent(double percent) {
	static GtkProgressBar *pbar = NULL;
	if (pbar == NULL)
		pbar = GTK_PROGRESS_BAR(gtk_builder_get_object(gui.builder, "progressbar1"));
	if (percent == PROGRESS_PULSATE) {
		gtk_progress_bar_pulse(pbar);
	}
	else {
		gtk_progress_bar_set_fraction(pbar, percent);
	}
}

struct progress_bar_idle_data {
	char *progress_bar_text;
	double progress_bar_percent;
};

static gboolean progress_bar_idle_callback(gpointer p) {
	struct progress_bar_idle_data *data = (struct progress_bar_idle_data *) p;

	if (data->progress_bar_text) {
		progress_bar_set_text(data->progress_bar_text);
		free(data->progress_bar_text);
	}
	if (data->progress_bar_percent != PROGRESS_NONE)
		progress_bar_set_percent(data->progress_bar_percent);
	free(data);
	return FALSE;	// only run once
}

// Thread-safe progress bar update.
// text can be NULL, percent can be -1 for pulsating, -2 for nothing, or between 0 and 1 for percent
void set_progress_bar_data(const char *text, double percent) {
	if (com.headless) {
		if (percent < 0.0) percent = 1.0;
		if (text)
			fprintf(stdout, "progress: %s, %4.2lf%%\n", text, percent*100.0);
		else fprintf(stdout, "\033[A\33[2KT\rprogress: %4.2lf%%\n", percent*100.0);
		// Warning: I don't know how to do that in other OS than GNU
		// On OS-X it works. On Windows ... well, I doubt it will be used

		/* progress update to the named pipe */
		char buf[30];
		snprintf(buf, 30, "progress: %4.2lf%%\n", percent*100.0);
		pipe_send_message(PIPE_PROGRESS, PIPE_NA, buf);
	} else {
		struct progress_bar_idle_data *data;
		data = malloc(sizeof(struct progress_bar_idle_data));
		data->progress_bar_text = text ? strdup(text) : NULL;
		data->progress_bar_percent = percent;
		assert(percent == PROGRESS_PULSATE || percent == PROGRESS_NONE ||
				(percent >= 0.0 && percent <= 1.0));
		gdk_threads_add_idle(progress_bar_idle_callback, data);
	}
}

/************************ M E S S A G E    L O G G I N G  ************************/

struct log_message {
	char *timestamp;
	char *message;
	const char* color;
};

static gboolean scroll_to_end(gpointer data) {
	GtkTextIter iter;
	GtkTextView *text = (GtkTextView *) data;
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	GtkTextMark *insert_mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	gtk_text_view_scroll_to_mark(text, insert_mark, 0.0, TRUE, 0.0, 1.0);

	return FALSE;
}

// The main thread internal function that does the printing.
static gboolean idle_messaging(gpointer p) {
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter iter;
	struct log_message *log = (struct log_message *) p;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(gtk_builder_get_object(gui.builder, "output"));
		tbuf = gtk_text_view_get_buffer(text);
	}

	if (log->message[0] == '\n' && log->message[1] == '\0') {
		gtk_text_buffer_get_start_iter(tbuf, &iter);
		gtk_text_buffer_insert(tbuf, &iter, log->message, strlen(log->message));
		free(log->timestamp);
		free(log->message);
		free(log);
		return FALSE;
	}

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->timestamp,
			strlen(log->timestamp), "bold", NULL);

	if (!log->color)
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->message,
				strlen(log->message), "normal", NULL);
	else
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->message,
				strlen(log->message), log->color, NULL);

	/* scroll to the end with a timeout
	 * just to be sure that everything is displayed
	 */
	g_timeout_add(50, scroll_to_end, (gpointer) text);

	free(log->timestamp);
	free(log->message);
	free(log);
	return FALSE;
}

// Send a log message to the console in the UI
void gui_log_message(const char* msg, const char* color) {
	if (com.headless)
		return;
	time_t now_sec = time(NULL);
	struct tm *now = localtime(&now_sec);
	char timestamp[30];
	g_snprintf(timestamp, sizeof(timestamp), "%.2d:%.2d:%.2d: ", now->tm_hour,
			now->tm_min, now->tm_sec);

	struct log_message *new_msg = malloc(sizeof(struct log_message));
	new_msg->timestamp = strdup(timestamp);
	new_msg->message = strdup(msg);
	new_msg->color = color;
	gdk_threads_add_idle(idle_messaging, new_msg);
}

void initialize_log_tags() {
	/* Create tags associated with the buffer for the output text. */
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(lookup_widget("output")));
	/* Tag with weight bold and tag name "bold" . */
	gtk_text_buffer_create_tag (tbuf, "bold", "weight", PANGO_WEIGHT_BOLD, NULL);
	/* Tag with style normal */
	gtk_text_buffer_create_tag (tbuf, "normal", "weight", PANGO_WEIGHT_NORMAL, NULL);
	/* Couleur Tags */
	gtk_text_buffer_create_tag (tbuf, "red", "foreground", "#e72828", NULL);
	gtk_text_buffer_create_tag (tbuf, "salmon", "foreground", "#ff9898", NULL);
	gtk_text_buffer_create_tag (tbuf, "green", "foreground", "#01b301", NULL);
	gtk_text_buffer_create_tag (tbuf, "blue", "foreground", "#7a7af8", NULL);
	gtk_text_buffer_create_tag (tbuf, "plum", "foreground", "#8e4585", NULL);
}

/********************** S A V I N G    T H E    L O G **********************/

static void save_log_file(gchar *filename) {
	GtkTextBuffer *log;
	GtkTextView *tv;
	GtkTextIter start, end;
	gchar *str;
	GError *error = NULL;

	tv = GTK_TEXT_VIEW(lookup_widget("output"));
	log = gtk_text_view_get_buffer(tv);
	gtk_text_buffer_get_bounds(log, &start, &end);
	str = gtk_text_buffer_get_text(log, &start, &end, FALSE);

	GFile *file = g_file_new_for_path(filename);
	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			siril_log_message(_("Cannot create logfile [%s]\n"), filename);
		}
		g_object_unref(file);
		return;
	}

	gsize bytes_written = 0;
	if (!g_output_stream_write_all(output_stream, str, strlen(str),
				&bytes_written, NULL, &error)) {
		g_warning("%s\n", error->message);
		g_clear_error(&error);
	}

	g_object_unref(output_stream);
	g_object_unref(file);
	g_free(str);
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("Log files (*.log)"));
	gtk_file_filter_add_pattern(f, "*.log");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

static void save_log_dialog() {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	gint res;
	gchar *filename;

	filename = build_timestamp_filename();
	filename = str_append(&filename, ".log");

	widgetdialog = siril_file_chooser_save(control_window, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, filename);
	gtk_file_chooser_set_local_only(dialog, FALSE);
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_file_chooser_get_filename(dialog);
		save_log_file(file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
	g_free(filename);
}

void on_export_log_button_clicked(GtkButton *button, gpointer user_data) {
	save_log_dialog();
}

void on_clear_log_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean ret = siril_confirm_dialog(_("Clear the log"),
			_("Are you sure you want to clear the log? There is no possible undo."), _("Clear the Log"));
	if (ret) {
		process_clear(0);
	}
}

/**************************** C U R S O R S ****************************/

struct _cursor_data {
	gboolean change;
	//	GdkCursorType cursor_type;
	const gchar* cursor_name;
};

/* thread-safe cursor change */
static gboolean idle_set_cursor(gpointer garg) {
	struct _cursor_data *arg = (struct _cursor_data*) garg;
	GdkCursor *cursor = NULL;
	static const gchar *current_name = NULL;

	GdkDisplay *display = gdk_display_get_default();
	GdkScreen *screen = gdk_screen_get_default();

	if (arg->change) {
		if (!current_name ||
				(arg->cursor_name && g_strcmp0(current_name, arg->cursor_name))) {
			cursor = gdk_cursor_new_from_name(display, arg->cursor_name);
			current_name = arg->cursor_name;
		}
		else {
			// it's the current cursor
			free(arg);
			return FALSE;
		}
	} else {
		if (!current_name) {
			// it's the current default
			free(arg);
			return FALSE;
		}
		current_name = NULL;
	}

	GList *list = gdk_screen_get_toplevel_windows(screen);
	for (GList *l = list; l; l = l->next) {
		GdkWindow *window = GDK_WINDOW(l->data);
		/* Passing NULL for the cursor argument to gdk_window_set_cursor()
		 * means that window will use the cursor of its parent window. Most
		 * windows should use this default.
		 */
		gdk_window_set_cursor(window, cursor);
		gdk_display_sync(gdk_window_get_display(window));
	}

	gdk_display_flush(display);
	g_list_free(list);

	free(arg);
	return FALSE;
}

void set_cursor_waiting(gboolean waiting) {
	if (com.headless) return;
	struct _cursor_data *arg = malloc(sizeof (struct _cursor_data));

	arg->change = waiting;
	arg->cursor_name = "progress";

	gui_function(idle_set_cursor, arg);
}

void set_cursor(const gchar* cursor_name) {
	struct _cursor_data *arg = malloc(sizeof (struct _cursor_data));

	arg->change = TRUE;
	arg->cursor_name = cursor_name;

	gui_function(idle_set_cursor, arg);
}
