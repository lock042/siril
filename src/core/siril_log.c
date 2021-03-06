/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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
#include "core/siril_date.h"
#include "core/command.h" // for process_clear()
#include "core/OS_utils.h"
#include "core/pipe.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"


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
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = gtk_file_chooser_get_filename(dialog);
		save_log_file(file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
	g_free(filename);
}

/* This function writes a message on Siril's console/log. It is not thread safe.
 * There is a limit in number of characters that it is able to write in one call: 1023.
 * Return value is the string printed from arguments, or NULL if argument was empty or
 * only newline. It is an allocated string and must not be freed. It can be
 * reused until next call to this function.
 */
char* siril_log_internal(const char* format, const char* color, va_list arglist) {
	static char *msg = NULL;

	if (msg == NULL) {
		msg = malloc(1024);
		msg[1023] = '\0';
	}

	vsnprintf(msg, 1023, format, arglist);

	if (msg == NULL || msg[0] == '\0')
		return NULL;

	if (msg[0] == '\n' && msg[1] == '\0') {
		fputc('\n', stdout);
		gui_log_message("\n", NULL);
		return NULL;
	}

	fprintf(stdout, "log: %s", msg);
	pipe_send_message(PIPE_LOG, PIPE_NA, msg);
	gui_log_message(msg, color);

	return msg;
}

char* siril_log_message(const char* format, ...) {
	va_list args;
	va_start(args, format);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, NULL, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

char* siril_log_color_message(const char* format, const char* color, ...) {
	va_list args;
	va_start(args, color);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, color, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

/************** Callbacks function ***********/

void on_clear_log_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean ret = siril_confirm_dialog(_("Clear the log"),
			_("Are you sure you want to clear the log? There is no possible undo."), _("Clear the Log"));
	if (ret) {
		process_clear(0);
	}
}

void on_export_log_button_clicked(GtkButton *button, gpointer user_data) {
	save_log_dialog();
}
