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
#include "core/processing.h"
#include "core/command_line_processor.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/script_menu.h"
#include "gui/utils.h"
#include "io/siril_pythonmodule.h"

#include <gtksourceview/gtksource.h>

// Statics declarations

GtkButton *button_python_pad_close = NULL, *button_python_pad_clear = NULL, *button_python_pad_open = NULL, *button_python_pad_save = NULL, *button_python_pad_execute = NULL;
GtkDialog *python_dialog = NULL;
GtkLabel *script_label = NULL;
GtkWidget *code_view = NULL;
GtkSourceBuffer *sourcebuffer = NULL;
GtkScrolledWindow *scrolled_window = NULL;
GtkComboBox *combo_language = NULL;
GtkSourceLanguageManager *language_manager = NULL;
GtkSourceLanguage *language = NULL;
GtkSourceStyleSchemeManager *stylemanager = NULL;
GtkSourceStyleScheme *scheme = NULL;

enum {
	LANG_PYTHON,
	LANG_SSF
};

static gint active_language = LANG_PYTHON;

void set_code_view_theme() {
	// Set the style scheme based on whether the Siril theme is light or dark
	// The core "Classic" and "Oblivion" themes are used: these should always be available
	stylemanager = gtk_source_style_scheme_manager_get_default();
	scheme = gtk_source_style_scheme_manager_get_scheme(stylemanager,
										com.pref.gui.combo_theme == 0 ? "oblivion" : "classic");
	if (scheme)
		gtk_source_buffer_set_style_scheme(sourcebuffer, scheme);
}

gboolean code_view_exists() {
	return (code_view != NULL);
}

void add_code_view(GtkBuilder *builder) {
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
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(code_view), GTK_WRAP_WORD_CHAR);
	gtk_source_view_set_auto_indent(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_insert_spaces_instead_of_tabs(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_indent_on_tab(GTK_SOURCE_VIEW(code_view), TRUE);

	gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_show_right_margin(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_view_set_highlight_current_line(GTK_SOURCE_VIEW(code_view), TRUE);
	gtk_source_buffer_set_highlight_matching_brackets(sourcebuffer, TRUE);


	///////

	// Set the GtkSourceView style depending on whether the Siril light or dark
	// theme is set.
	set_code_view_theme();

	// Get the GtkSourceLanguageManager and set the Python language
	language_manager = gtk_source_language_manager_get_default();
	language = gtk_source_language_manager_get_language(language_manager, "python3");
	if (language == NULL) {
		siril_debug_print("Could not find Python language definition\n");
	} else {
		gtk_source_buffer_set_language(sourcebuffer, language);
		gtk_source_view_set_show_line_marks(GTK_SOURCE_VIEW(code_view), TRUE);
		// Create mark attributes for folding
		GtkSourceMarkAttributes *fold_attributes = gtk_source_mark_attributes_new();
		// Customize fold mark appearance
		GdkRGBA color;
		gdk_rgba_parse(&color, "#0000FF");  // Blue color for fold markers
		gtk_source_mark_attributes_set_background(fold_attributes, &color);
		gtk_source_view_set_mark_attributes(GTK_SOURCE_VIEW(code_view),
											"fold",      // Category name
											fold_attributes,  // Attributes object
											1);          // Priority
		// Release references
		g_object_unref(fold_attributes);
		siril_debug_print("Set buffer language to python3\n");
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
		// GtkComboBox
		combo_language = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "python_pad_language"));
		// GtkButton
		button_python_pad_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_close"));
		button_python_pad_clear = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_clear"));
		button_python_pad_open = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_open"));
		button_python_pad_save = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_save"));
		button_python_pad_execute = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_execute"));
		// GtkDialog
		python_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "python_dialog"));
		// GtkSCrolledWindow
		scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "python_scrolled_window"));
		// GtkLabel
		script_label = GTK_LABEL(gtk_label_new("unsaved"));
		gtk_widget_show(GTK_WIDGET(script_label));
		GtkWidget *parent = lookup_widget("python_pad_box");
		if (parent) {
			gtk_box_pack_start(GTK_BOX(parent), GTK_WIDGET(script_label), TRUE, TRUE, 0);
			gtk_box_reorder_child(GTK_BOX(parent), GTK_WIDGET(scrolled_window), 1);
		}
		add_code_view(gui.builder);
	}
}

int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("python_dialog");
	gtk_combo_box_set_active(combo_language, active_language);
	gtk_widget_grab_focus(code_view);
	return 0;
}

void on_python_pad_close_clicked(GtkWidget *widget, gpointer user_data) {
	siril_close_dialog("python_dialog");
}

void on_python_pad_language_changed(GtkComboBox *combo, gpointer user_data) {
	active_language = gtk_combo_box_get_active(combo);
	if (active_language == LANG_PYTHON) {
		language = gtk_source_language_manager_get_language(language_manager, "python");
	} else {
		language = gtk_source_language_manager_get_language(language_manager, "bash"); // in the absence of a proper SSF definition, this will do
	}
	if (language == NULL) {
		siril_debug_print("Could not find  language definition\n");
	} else {
		gtk_source_buffer_set_language(sourcebuffer, language);
		gtk_source_view_set_show_line_marks(GTK_SOURCE_VIEW(code_view), TRUE);
		// Create mark attributes for folding
		GtkSourceMarkAttributes *fold_attributes = gtk_source_mark_attributes_new();
		// Customize fold mark appearance
		GdkRGBA color;
		gdk_rgba_parse(&color, "#0000FF");  // Blue color for fold markers
		gtk_source_mark_attributes_set_background(fold_attributes, &color);
		gtk_source_view_set_mark_attributes(GTK_SOURCE_VIEW(code_view),
											"fold",      // Category name
											fold_attributes,  // Attributes object
											1);          // Priority
		// Release references
		g_object_unref(fold_attributes);
		siril_debug_print("Set buffer language to python\n");
		siril_debug_print("Set buffer language\n");
	}
}

static gchar* read_stream_into_gchar(GInputStream* stream, gsize* length, GError** error) {
	// TODO: is there already a Siril function that does this?
	gsize buffer_size = 4096;  // Initial buffer size
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
	GtkWidget *dialog = gtk_file_chooser_dialog_new(_("Open Script"),
			GTK_WINDOW(lookup_widget("control_window")),
			GTK_FILE_CHOOSER_ACTION_OPEN,
			_("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);

	// Set up a filter for .py files
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.py");
	gtk_file_filter_add_pattern(filter, "*.ssf");
	gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
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
		gsize length = 0;
		gchar *text = read_stream_into_gchar(input_stream, &length, &error);
		gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sourcebuffer), text, -1);
		g_input_stream_close(input_stream, NULL, NULL);
		if (length > 0) {
			gtk_label_set_text(script_label, filename);
		} else {
			gtk_label_set_text(script_label, _("unsaved"));
		}
		if (!strcmp(get_filename_ext(filename), "py")) {
			gtk_combo_box_set_active(combo_language, LANG_PYTHON);
		} else if (strcmp(get_filename_ext(filename), "py")) {
			gtk_combo_box_set_active(combo_language, LANG_SSF);
		}
		gtk_widget_queue_draw(GTK_WIDGET(python_dialog));
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
	gtk_file_filter_add_pattern(filter, "*.ssf");
	gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
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
			gtk_label_set_text(script_label, filename);
			gtk_widget_queue_draw(GTK_WIDGET(python_dialog));
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
	switch (active_language) {
		case LANG_PYTHON:;
			execute_python_script_async(text, FALSE);
			break;
		case LANG_SSF:;
			GInputStream *input_stream = g_memory_input_stream_new_from_data(text, strlen(text), NULL);
			if (get_thread_run()) {
				PRINT_ANOTHER_THREAD_RUNNING;
				g_object_unref(input_stream);
				return;
			}

			if (com.script_thread)
				g_thread_join(com.script_thread);

			/* Switch to console tab */
			control_window_switch_to_tab(OUTPUT_LOGS);
			/* Last thing before running the script, disable widgets except for Stop */
			script_widgets_enable(FALSE);

			/* Run the script */
			siril_log_message(_("Starting script\n"));
			com.script_thread = g_thread_new("script", execute_script, input_stream);
			break;
		default:
			siril_debug_print("Error: unknown script language\n");
			break;
	}
	// TODO: Neither case properly cleans up text yet
}

void on_button_python_pad_clear_clicked(GtkWidget *widget, gpointer user_data) {
    g_print("script_label pointer: %p\n", (void*)script_label);  // Check if it's NULL
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	if (siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
		gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		g_print("Setting label to 'unsaved'\n");
		gtk_label_set_text(script_label, _("unsaved"));
		g_print("Label text is now: %s\n", gtk_label_get_text(script_label));
		gtk_widget_queue_draw(GTK_WIDGET(python_dialog));
	}
}
