// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/command_line_processor.h"
#include "gui/documentation.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/script_menu.h"
#include "gui/utils.h"
#include "io/siril_pythonmodule.h"

#include <gtksourceview/gtksource.h>

// Statics declarations
enum {
	LANG_PYTHON,
	LANG_SSF
};

static GtkButton *button_python_pad_close = NULL, *button_python_pad_clear = NULL, *button_python_pad_open = NULL, *button_python_pad_save = NULL, *button_python_pad_execute = NULL;
static GtkLabel *language_label = NULL;
static GtkLabel *script_label = NULL;
static GtkSourceView *code_view = NULL;
static GtkSourceBuffer *sourcebuffer = NULL;
static GtkScrolledWindow *scrolled_window = NULL;
static GtkComboBox *combo_language = NULL;
static GtkSourceLanguageManager *language_manager = NULL;
static GtkSourceLanguage *language = NULL;
static GtkSourceStyleSchemeManager *stylemanager = NULL;
static GtkSourceStyleScheme *scheme = NULL;
static GFile *current_file = NULL;
static GtkWindow *editor_window = NULL;
static GtkWindow *main_window = NULL;
static GtkCheckMenuItem *radio_py = NULL;
static GtkCheckMenuItem *radio_ssf = NULL;
static gint active_language = LANG_PYTHON;
static gboolean buffer_modified = FALSE;

// Forward declarations
void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_execute(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_new(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_select_language(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_python_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_command_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data);
static void on_buffer_modified_changed(GtkTextBuffer *buffer, gpointer user_data);
void set_language();

static GActionEntry editor_actions[] = {
	{ "open", on_action_file_open, NULL, NULL, NULL },
	{ "save", on_action_file_save, NULL, NULL, NULL },
	{ "save_as", on_action_file_save_as, NULL, NULL, NULL },
	{ "new", on_action_file_new, NULL, NULL, NULL },
	{ "close", on_action_file_close, NULL, NULL, NULL },
	{ "execute", on_action_file_execute, NULL, NULL, NULL },
	{ "python_doc", on_action_python_doc, NULL, NULL, NULL },
	{ "command_doc", on_action_command_doc, NULL, NULL, NULL }
};

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

// Add this function to initialize the actions
static void setup_editor_actions(GtkWindow *window) {
	GSimpleActionGroup *action_group;
	action_group = g_simple_action_group_new();

	g_action_map_add_action_entries(G_ACTION_MAP(action_group),
								editor_actions,
								G_N_ELEMENTS(editor_actions),
								window);

	// Insert the action group into the window's action map
	gtk_widget_insert_action_group(GTK_WIDGET(window),
								"editor",  // prefix for all actions
								G_ACTION_GROUP(action_group));

	// The action group can be unreferenced here as the window now owns it
	g_object_unref(action_group);
}

void add_code_view(GtkBuilder *builder) {
	// Create a new GtkSourceBuffer (glade doesn't handle this properly)
	sourcebuffer = gtk_source_buffer_new(NULL);
	gtk_text_view_set_buffer(GTK_TEXT_VIEW(code_view), GTK_TEXT_BUFFER(sourcebuffer));

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
}

// Statics init
void python_scratchpad_init_statics() {
	if (editor_window == NULL) {
		// GtkWindow
		editor_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "python_window"));
		main_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
		// GtkSourceView
		code_view = GTK_SOURCE_VIEW(gtk_builder_get_object(gui.builder, "code_view"));
		sourcebuffer = GTK_SOURCE_BUFFER(gtk_builder_get_object(gui.builder, "sourcebuffer"));
		// GtkCheckMenuItem
		radio_py = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "radio_py"));
		radio_ssf = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "radio_ssf"));
		// GtkComboBox
		combo_language = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "python_pad_language"));
		// GtkButton
		button_python_pad_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_close"));
		button_python_pad_clear = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_clear"));
		button_python_pad_open = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_open"));
		button_python_pad_save = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_save"));
		button_python_pad_execute = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_execute"));
		// GtkSCrolledWindow
		scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "python_scrolled_window"));
		// GtkLabel
		language_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "script_language_label"));
		script_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "script_editor_label"));
		gtk_widget_show(GTK_WIDGET(script_label));
		add_code_view(gui.builder);
		gtk_window_set_transient_for(GTK_WINDOW(editor_window), GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
		g_signal_connect(sourcebuffer, "modified-changed",
						G_CALLBACK(on_buffer_modified_changed), NULL);

		// Initialize with "unsaved" if no file is loaded
		if (!current_file) {
			gtk_label_set_text(script_label, "unsaved");
		}
	}
}

static void update_title_with_modification() {
	const gchar *current_text = gtk_label_get_text(script_label);
	if (!current_text || !*current_text) return;

	// If the title already ends with *, don't add another
	if (g_str_has_suffix(current_text, "*")) {
		if (!buffer_modified) {
			gchar *base_name = g_strndup(current_text, strlen(current_text) - 1);
			gtk_label_set_text(script_label, base_name);
			g_free(base_name);
		}
	} else if (buffer_modified) {
		gchar *new_title = g_strdup_printf("%s*", current_text);
		gtk_label_set_text(script_label, new_title);
		g_free(new_title);
	}
}

static void update_title(GFile *file) {
	if (file) {
		char *basename = g_file_get_basename(file);
		gtk_label_set_text(script_label, basename);

		char *suffix = strrchr(basename, '.');
		if (suffix != NULL) {
			if (g_ascii_strcasecmp(suffix, ".py") == 0) {
				gtk_check_menu_item_set_active(radio_py, TRUE);
			} else if (g_ascii_strcasecmp(suffix, ".ssf") == 0) {
				gtk_check_menu_item_set_active(radio_ssf, TRUE);
			}
		}

		g_free(basename);
	} else {
		gtk_label_set_text(script_label, "unsaved");
	}
	buffer_modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
}

// File handling

void load_file_complete(GObject *loader, GAsyncResult *result, gpointer user_data) {
	GError *error = NULL;
	if (!gtk_source_file_loader_load_finish(GTK_SOURCE_FILE_LOADER(loader), result, &error)) {
		g_printerr("Error loading file: %s\n", error->message);
		gtk_label_set_text(script_label, _(""));
		g_error_free(error);
	}
	g_object_unref(loader);
}

void load_file(GFile *file) {
	GtkSourceBuffer *buffer = sourcebuffer;
	GtkSourceFile *source_file = gtk_source_file_new();
	gtk_source_file_set_location(source_file, file);

	GtkSourceFileLoader *loader = gtk_source_file_loader_new(buffer, source_file);

	// Start the async load operation
	gtk_source_file_loader_load_async(loader,
									G_PRIORITY_DEFAULT,
									NULL,  // No cancellable
									NULL,  // No progress callback
									NULL,  // No progress data
									NULL,  // No progress notify
									load_file_complete,
									loader); // No user data
	g_object_unref(source_file);
	// loader will be unreferenced in the callback
}

// If you need to save the file:
void save_file_complete(GObject *saver, GAsyncResult *result, gpointer user_data) {
	GError *error = NULL;
	if (!gtk_source_file_saver_save_finish(GTK_SOURCE_FILE_SAVER(saver), result, &error)) {
		g_printerr("Error saving file: %s\n", error->message);
		g_error_free(error);
	}
	g_object_unref(saver);
}

void save_file(GFile *file) {
	GtkSourceBuffer *buffer = sourcebuffer;
	GtkSourceFile *source_file = gtk_source_file_new();
	gtk_source_file_set_location(source_file, file);

	GtkSourceFileSaver *saver = gtk_source_file_saver_new(buffer, source_file);

	// Start the async save operation
	gtk_source_file_saver_save_async(saver,
								G_PRIORITY_DEFAULT,
								NULL,  // No cancellable
								NULL,  // No progress callback
								NULL,  // No progress data
								NULL,  // No progress notify
								save_file_complete,
								NULL); // No user data
	buffer_modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
	update_title(file);
	g_object_unref(source_file);
}

void on_action_file_new(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	if (siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;
		gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		buffer_modified = FALSE;
		gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
		gtk_label_set_text(script_label, "unsaved");
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}
}

void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	gint char_count = gtk_text_buffer_get_char_count(GTK_TEXT_BUFFER(sourcebuffer));
	gboolean is_empty = (char_count == 0);
	if (!is_empty) {
		on_action_file_new(action, parameter, user_data);
	} else {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;
		gtk_label_set_text(script_label, "");
	}
	gtk_widget_hide(GTK_WIDGET(editor_window));
}

void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (G_IS_OBJECT(current_file))
		g_object_unref(current_file);
	current_file = NULL;
	GtkWidget *dialog = gtk_file_chooser_dialog_new(_("Open Script"),
			GTK_WINDOW(lookup_widget("control_window")),
			GTK_FILE_CHOOSER_ACTION_OPEN,
			_("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);

	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.py");
	gtk_file_filter_add_pattern(filter, "*.ssf");
	gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
	gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

	gint res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
		control_window_switch_to_tab(OUTPUT_LOGS);
		load_file(file);
		current_file = g_object_ref(file);
		update_title(current_file);
		g_object_unref(file);
	}

	gtk_widget_destroy(dialog);
	gtk_window_present(editor_window);
}

void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *dialog = gtk_file_chooser_dialog_new(_("Save Script As"),
			GTK_WINDOW(lookup_widget("control_window")),
			GTK_FILE_CHOOSER_ACTION_SAVE,
			_("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);

	// Enable overwrite confirmation
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);

	// Set up file filter for script files
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.py");
	gtk_file_filter_add_pattern(filter, "*.ssf");
	gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
	gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

	// If there's a current file, set it as the default
	if (G_IS_OBJECT(current_file)) {
		gtk_file_chooser_set_file(GTK_FILE_CHOOSER(dialog), current_file, NULL);
	}

	gint res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
		control_window_switch_to_tab(OUTPUT_LOGS);

		// Update current_file with the newly selected file
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = g_object_ref(file);
		update_title(current_file);

		save_file(file);  // This now uses the async version
		g_object_unref(file);
	}

	gtk_widget_destroy(dialog);
	gtk_window_present(editor_window);
}

void on_action_file_save(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_file) {
		// No file has been saved yet, so call Save As
		on_action_file_save_as(action, parameter, user_data);
		return;
	}

	// We have a current file, just save directly to it
	control_window_switch_to_tab(OUTPUT_LOGS);
	save_file(current_file);
}

int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data) {
	if (!editor_window) {
		python_scratchpad_init_statics();
	}
	setup_editor_actions(editor_window);

	// If the window isn't visible, let's position it relative to the main window
	if (!gtk_widget_get_visible(GTK_WIDGET(editor_window))) {
		// Get the main window's position and size
		gint main_x, main_y, main_width, main_height;
		gtk_window_get_position(main_window, &main_x, &main_y);
		gtk_window_get_size(main_window, &main_width, &main_height);

		// Position the editor window slightly offset from the main window
		// This creates a cascading effect common in MDI applications
		gtk_window_move(editor_window, main_x + 50, main_y + 50);
	}

	gtk_label_set_text(language_label, active_language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));

	// Show the window and bring it to front

	gtk_window_present_with_time(editor_window, GDK_CURRENT_TIME);

	// Set the correct check menu item active
	gtk_check_menu_item_set_active(radio_py, active_language == LANG_PYTHON);
	gtk_check_menu_item_set_active(radio_ssf, active_language == LANG_SSF);

	// Focus the SourceView
	gtk_widget_grab_focus(GTK_WIDGET(code_view));

	return 0;
}

// Handler for main window state changes
gboolean on_main_window_state_changed(GtkWidget *widget,
									GdkEventWindowState *event,
									gpointer user_data) {
	GtkWindow *editor_window = GTK_WINDOW(user_data);

	if (event->changed_mask & GDK_WINDOW_STATE_ICONIFIED) {
		if (event->new_window_state & GDK_WINDOW_STATE_ICONIFIED) {
			// Main window was minimized, minimize editor too
			gtk_window_iconify(editor_window);
		} else {
			// Main window was restored, restore editor too
			gtk_window_deiconify(editor_window);
		}
	}

	return FALSE;
}

void set_language() {
	if (active_language == LANG_PYTHON) {
		language = gtk_source_language_manager_get_language(language_manager, "python");
		gtk_label_set_text(language_label, _("Python Script"));
	} else if (active_language == LANG_SSF) {
		language = gtk_source_language_manager_get_language(language_manager, "sh");
		gtk_label_set_text(language_label, _("Siril Script File"));
	}
	if (language == NULL) {
		siril_debug_print("Could not find language definition\n");
	} else {
		gtk_source_buffer_set_language(sourcebuffer, language);
	}
}

void on_language_set(GtkRadioMenuItem *item, gpointer user_data) {
	if (gtk_check_menu_item_get_active(radio_py)) {
		active_language = LANG_PYTHON;
	} else {
		active_language = LANG_SSF;
	}
	set_language();
}

void setup_python_editor_window() {
	GtkWindow *editor_window = GTK_WINDOW(lookup_widget("python_window"));
	GtkWindow *main_window = GTK_WINDOW(lookup_widget("control_window"));

	// Make the editor window hide instead of destroy when closed
	g_signal_connect(editor_window, "delete-event",
					G_CALLBACK(gtk_widget_hide_on_delete), NULL);

	// Optional: Make editor window minimize when main window minimizes
	g_signal_connect(main_window, "window-state-event",
					G_CALLBACK(on_main_window_state_changed), editor_window);
}

void on_action_file_execute(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
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

static void on_buffer_modified_changed(GtkTextBuffer *buffer, gpointer user_data) {
	buffer_modified = gtk_text_buffer_get_modified(buffer);
	update_title_with_modification();
}

void on_action_python_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("scripts/api.html");
}

void on_action_command_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("Commands.html");
}

void on_editor_syntax_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_buffer_set_highlight_syntax(sourcebuffer, gtk_check_menu_item_get_active(item));
}

void on_editor_bracketmatch_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_buffer_set_highlight_matching_brackets(sourcebuffer, gtk_check_menu_item_get_active(item));
}

void on_editor_rmargin_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_show_right_margin(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_linenums_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_show_line_numbers(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_linemarks_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_show_line_marks(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_highlightcurrentline_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_highlight_current_line(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_autoindent_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_auto_indent(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_indentontab_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_indent_on_tab(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_smartbs_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_smart_backspace(code_view, gtk_check_menu_item_get_active(item));
}

void on_editor_smarthomeend_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_source_view_set_smart_home_end(code_view, gtk_check_menu_item_get_active(item));
}
