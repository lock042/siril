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

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "python/siril_python.h"
#include "gui/dialogs.h"
#include "gui/utils.h"

int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("python_dialog");
	return 0;
}


void on_python_pad_close_clicked(GtkWidget *widget, gpointer user_data) {
	siril_close_dialog("python_dialog");
}

void on_button_python_pad_save_clicked(GtkWidget *widget, gpointer user_data) {
	// Get the text buffer from the widget
	GtkTextBuffer *buffer = GTK_TEXT_BUFFER(lookup_widget("python_textbuffer"));

	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(buffer, &start, &end);

	// Get the text between start and end iterators
	char *text = gtk_text_buffer_get_text(buffer, &start, &end, FALSE);

	// Create a file chooser dialog for saving
	GtkWidget *dialog = gtk_file_chooser_dialog_new(_("Save Python Script"),
			GTK_WINDOW(lookup_widget("control_window")),
			GTK_FILE_CHOOSER_ACTION_SAVE,
			_("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);

	// Set up a filter for .py files
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.py");
	gtk_file_filter_set_name(filter, _("Python Files (*.py)"));
	gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

	// Show the dialog and capture the response
	gint res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		// Get the filename from the file chooser
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		char *filename = gtk_file_chooser_get_filename(chooser);

		// Expand the home directory in the filename (if applicable)
		expand_home_in_filename(filename, 256);

		// Switch to console tab (assuming this logs the save)
		control_window_switch_to_tab(OUTPUT_LOGS);

		// Create a GFile for the chosen path
		GFile *file = g_file_new_for_path(filename);
		GError *error = NULL;

		// Open an output stream for writing to the file
		GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, TRUE, G_FILE_CREATE_NONE, NULL, &error);
		if (output_stream == NULL) {
			// Handle file opening errors
			if (error != NULL) {
				siril_log_message(_("File [%s] could not be opened for writing: %s\n"), filename, error->message);
				g_clear_error(&error);
			}
			g_object_unref(file);
			g_free(filename);
			gtk_widget_destroy(dialog);
			return;
		}

		// Write the text to the output stream
		gsize bytes_written;
		gboolean success = g_output_stream_write_all(output_stream, text, strlen(text), &bytes_written, NULL, &error);
		if (!success) {
			// Handle writing errors
			siril_log_message(_("Error writing to file [%s]: %s\n"), filename, error->message);
			g_clear_error(&error);
		} else {
			siril_log_message(_("Successfully saved file [%s]\n"), filename);
		}

		// Clean up and release resources
		g_output_stream_close(output_stream, NULL, NULL);
		g_object_unref(output_stream);
		g_object_unref(file);
		g_free(filename);
	}

	// Free the text and destroy the dialog
	g_free(text);
	gtk_widget_destroy(dialog);
}


void on_button_python_pad_execute_clicked(GtkWidget *widget, gpointer user_data) {
	GtkTextBuffer *buffer = GTK_TEXT_BUFFER(lookup_gobject("python_textbuffer"));
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(buffer, &start, &end);
	// Get the text
	char *text = gtk_text_buffer_get_text(buffer, &start, &end, FALSE);
	run_python_script_in_python_thread(text, FALSE);
}

void on_button_python_pad_clear_clicked(GtkWidget *widget, gpointer user_data) {
	GtkTextBuffer *buffer = GTK_TEXT_BUFFER(lookup_widget("python_textbuffer"));
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(buffer, &start, &end);
	gtk_text_buffer_delete(buffer, &start, &end);
}
