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

#include <gtksourceview/gtksource.h>

// Statics declarations

GtkButton *button_python_pad_close = NULL, *button_python_pad_clear = NULL, *button_python_pad_open = NULL, *button_python_pad_save = NULL, *button_python_pad_execute = NULL;
GtkDialog *python_dialog = NULL;
GtkLabel *script_label = NULL;
GtkWidget *code_view = NULL;
GtkSourceBuffer *sourcebuffer = NULL;

void add_code_view(GtkBuilder *builder) {
	GtkSourceLanguageManager *language_manager = NULL;
	GtkSourceLanguage *language = NULL;
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(builder, "python_scrolled_window"));

	// Create a new GtkSourceView
	code_view = gtk_source_view_new();
	gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_widget_set_vexpand(code_view, TRUE);
	gtk_widget_set_hexpand(code_view, TRUE);

	// Create a new GtkSourceBuffer
	sourcebuffer = gtk_source_buffer_new(NULL);
	gtk_text_view_set_buffer(GTK_TEXT_VIEW(code_view), GTK_TEXT_BUFFER(sourcebuffer));

	// Set additional properties for GtkSourceView
	gtk_text_view_set_monospace(GTK_TEXT_VIEW(code_view), TRUE);
	gtk_source_view_set_auto_indent(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_insert_spaces_instead_of_tabs(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_indent_on_tab(GTK_SOURCE_VIEW(code_view), TRUE);

	// Get the GtkSourceLanguageManager and set the Python language
	language_manager = gtk_source_language_manager_get_default();
	language = gtk_source_language_manager_get_language(language_manager, "python");
	if (language == NULL) {
		siril_debug_print("Could not find Python language definition\n");
	} else {
		gtk_source_buffer_set_language(sourcebuffer, language);
		siril_debug_print("Set buffer language to Python\n");
	}

	// Enable syntax highlighting
	gtk_source_buffer_set_highlight_syntax(sourcebuffer, TRUE);

	// Add the GtkSourceView to the GtkScrolledWindow
	gtk_container_add(GTK_CONTAINER(scrolled_window), code_view);

	// Ensure the GtkSourceView is visible and expanded
	gtk_widget_set_visible(GTK_WIDGET(code_view), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(scrolled_window), TRUE);

	// Show all widgets in the dialog
	gtk_widget_show_all(GTK_WIDGET(python_dialog));
}

// Statics init
void python_scratchpad_init_statics() {
	if (python_dialog == NULL) {
		// GtkButton
		button_python_pad_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_close"));
		button_python_pad_clear = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_clear"));
		button_python_pad_open = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_open"));
		button_python_pad_save = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_save"));
		button_python_pad_execute = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_execute"));
		// GtkDialog
		python_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "python_dialog"));
		// GtkLabel
		script_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "script_label"));
		add_code_view(gui.builder);
	}
}

int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("python_dialog");
	return 0;
}


void on_python_pad_close_clicked(GtkWidget *widget, gpointer user_data) {
	siril_close_dialog("python_dialog");
}

static gchar* read_stream_into_gchar(GInputStream* stream, gsize* length, GError** error) {
	// TODO: is there already a Siril function that does this?
	gsize buffer_size = 4096;  // Initial buffer size, you can adjust this.
	gsize total_bytes_read = 0;
	gssize bytes_read;
	gchar* buffer = g_malloc(buffer_size);  // Allocate an initial buffer

	if (!buffer) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	while ((bytes_read = g_input_stream_read(stream, buffer + total_bytes_read, buffer_size - total_bytes_read, NULL, error)) > 0) {
		total_bytes_read += bytes_read;

		if (total_bytes_read == buffer_size) {
			// Expand buffer size if necessary
			buffer_size *= 2;
			buffer = g_realloc(buffer, buffer_size);

			if (!buffer) {
				PRINT_ALLOC_ERR;
				return NULL;
			}
		}
	}

	if (bytes_read == -1) {
		// Handle read error
		g_free(buffer);
		return NULL;
	}

	// Null-terminate the result and resize buffer to exact length
	buffer = g_realloc(buffer, total_bytes_read + 1);
	buffer[total_bytes_read] = '\0';

	if (length) {
		*length = total_bytes_read;  // Set the length if requested
	}

	return buffer;
}

void on_button_python_pad_open_clicked(GtkWidget *widget, gpointer user_data) {
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
		char *filename;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		filename = gtk_file_chooser_get_filename(chooser);

		/* Expand home directory in filename */
		expand_home_in_filename(filename, 256);

		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

		GFile *file = g_file_new_for_path(filename);
		GError *error = NULL;
		GInputStream *input_stream = (GInputStream*) g_file_read(file, NULL, &error);

		if (input_stream == NULL) {
			if (error != NULL) {
				g_clear_error(&error);
				siril_log_message(_("File [%s] does not exist\n"), filename);
			}
			g_object_unref(file);
			g_free(filename);
			gtk_widget_destroy(dialog);
			return;
		}
		gsize length;
		gchar *text = read_stream_into_gchar(input_stream, &length, &error);
		gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sourcebuffer), text, -1);
		g_input_stream_close(input_stream, NULL, NULL);
		g_object_unref(input_stream);
		g_free(text);
		gtk_widget_destroy(dialog);
	}
}

void on_python_dialog_show(GtkWidget *widget, gpointer user_data) {
	if (!python_dialog)
		python_scratchpad_init_statics();
}

void on_button_python_pad_save_clicked(GtkWidget *widget, gpointer user_data) {
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);

	// Get the text between start and end iterators
	char *text = gtk_text_buffer_get_text(GTK_TEXT_BUFFER(sourcebuffer), &start, &end, FALSE);

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
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	// Get the text
	char *text = gtk_text_buffer_get_text(GTK_TEXT_BUFFER(sourcebuffer), &start, &end, FALSE);
	run_python_script_in_python_thread(text, FALSE);
}

void on_button_python_pad_clear_clicked(GtkWidget *widget, gpointer user_data) {
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
}
