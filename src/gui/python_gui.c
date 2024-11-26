// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include <gtksourceview/gtksource.h>

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

#include "python_gui.h"

// Statics declarations
enum {
	LANG_PYTHON,
	LANG_SSF
};

static GtkButton *button_python_pad_close = NULL, *button_python_pad_clear = NULL, *button_python_pad_open = NULL, *button_python_pad_save = NULL, *button_python_pad_execute = NULL;
static GtkLabel *language_label = NULL;
static GtkLabel *find_label = NULL;
static GtkSourceView *code_view = NULL;
static GtkSourceBuffer *sourcebuffer = NULL;
static GtkSourceMap *map = NULL;
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
static GtkCheckMenuItem *syncheck = NULL;
static GtkCheckMenuItem *bracketcheck = NULL;
static GtkCheckMenuItem *rmargincheck = NULL;
static GtkCheckMenuItem *linenumcheck = NULL;
static GtkCheckMenuItem *linemarkcheck = NULL;
static GtkCheckMenuItem *currentlinecheck = NULL;
static GtkCheckMenuItem *autoindentcheck = NULL;
static GtkCheckMenuItem *indenttabcheck = NULL;
static GtkCheckMenuItem *smartbscheck = NULL;
static GtkCheckMenuItem *homeendcheck = NULL;
static GtkCheckMenuItem *showspaces = NULL;
static GtkCheckMenuItem *shownewlines = NULL;
static GtkCheckMenuItem *minimap = NULL;
static GtkSourceSpaceDrawer* space_drawer = NULL;
static GtkBox *codeviewbox = NULL;

GtkRevealer *find_revealer = NULL;
GtkEntry *find_entry = NULL;
GtkWidget *find_overlay = NULL;
static GtkButton *go_up_button = NULL;
static GtkButton *go_down_button = NULL;

static gint active_language = LANG_PYTHON;
static gboolean buffer_modified = FALSE;

// Forward declarations
void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_execute(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_new(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_python_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_command_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_set_rmarginpos(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_undo(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_redo(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_cut(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_copy(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_paste(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_find(GSimpleAction *action, GVariant *parameter, gpointer user_data);
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
	{ "command_doc", on_action_command_doc, NULL, NULL, NULL },
	{ "set_rmarginpos", on_set_rmarginpos, NULL, NULL, NULL },
	{ "undo", on_undo, NULL, NULL, NULL },
	{ "redo", on_redo, NULL, NULL, NULL },
	{ "cut", on_cut, NULL, NULL, NULL },
	{ "copy", on_copy, NULL, NULL, NULL },
	{ "paste", on_paste, NULL, NULL, NULL },
	{ "find", on_find, NULL, NULL, NULL }
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
	/* initialize minimap */
	map = GTK_SOURCE_MAP(gtk_source_map_new());

 // Apply CSS to force font size to 1

	GtkCssProvider* provider = gtk_css_provider_new();
	gtk_css_provider_load_from_data(provider,
		"* { font-size: 1px; }",
		-1, NULL);
	gtk_style_context_add_provider(
		gtk_widget_get_style_context(GTK_WIDGET(map)),
		GTK_STYLE_PROVIDER(provider),
		GTK_STYLE_PROVIDER_PRIORITY_USER
	);
	g_object_unref(provider);

	gtk_widget_show(GTK_WIDGET(map));
	gtk_source_map_set_view(map, code_view);
	gtk_box_pack_start(GTK_BOX(codeviewbox), GTK_WIDGET(map), FALSE, FALSE, 0);

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
	}

	// Enable syntax highlighting
	gtk_source_buffer_set_highlight_syntax(sourcebuffer, TRUE);

	// Get the space drawer
	space_drawer = gtk_source_view_get_space_drawer(code_view);

	// Configure which types of spaces to show
	GtkSourceSpaceTypeFlags space_types = 0;

	GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

	// Enable space drawing with marks
	gtk_source_space_drawer_set_types_for_locations(space_drawer,
												space_locations,
												space_types);

	gtk_source_space_drawer_set_enable_matrix(space_drawer, TRUE);
}

/** Code for the find feature ***/

static void setup_find_overlay() {
	// Create a new overlay
	GtkWidget *new_overlay = gtk_overlay_new();
	gtk_widget_show(new_overlay);

	// Replace scrolled window with overlay in its parent
	GtkWidget *scrolled_parent = gtk_widget_get_parent(GTK_WIDGET(scrolled_window));
	g_object_ref(scrolled_window);
	gtk_container_remove(GTK_CONTAINER(scrolled_parent), GTK_WIDGET(scrolled_window));
	gtk_container_add(GTK_CONTAINER(scrolled_parent), new_overlay);

	// Add scrolled window to overlay
	gtk_container_add(GTK_CONTAINER(new_overlay), GTK_WIDGET(scrolled_window));
	gtk_widget_set_vexpand(GTK_WIDGET(scrolled_window), TRUE);
	g_object_unref(scrolled_window);

	// Move find_overlay to new overlay
	g_object_ref(find_overlay);
	gtk_overlay_add_overlay(GTK_OVERLAY(new_overlay), GTK_WIDGET(find_overlay));
	g_object_unref(find_overlay);

	// Configure find_overlay position
	gtk_widget_set_halign(GTK_WIDGET(find_overlay), GTK_ALIGN_END);
	gtk_widget_set_valign(GTK_WIDGET(find_overlay), GTK_ALIGN_START);
	gtk_widget_set_margin_top(GTK_WIDGET(find_overlay), 6);
	gtk_widget_set_margin_end(GTK_WIDGET(find_overlay), 6);

	gtk_revealer_set_transition_type(find_revealer, GTK_REVEALER_TRANSITION_TYPE_SLIDE_DOWN);
	gtk_revealer_set_reveal_child(find_revealer, FALSE);
}

void toggle_find_overlay(gboolean show) {
	gtk_revealer_set_reveal_child(find_revealer, show);

	if (show) {
		gtk_widget_grab_focus(GTK_WIDGET(find_entry));
	}
}

void hide_find_box() {
	gtk_revealer_set_transition_type(find_revealer, GTK_REVEALER_TRANSITION_TYPE_SLIDE_DOWN);
	gtk_revealer_set_reveal_child(find_revealer, FALSE);
}

gboolean on_find_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	hide_find_box();
	return FALSE;
}


void on_find_entry_activate(GtkEntry *entry, gpointer user_data) {
	hide_find_box();
}

gboolean on_find_entry_key_press_event(GtkWidget *widget, GdkEventKey *event, gpointer data) {
	if (event->keyval == GDK_KEY_Escape) {
		hide_find_box();
		return TRUE;
	}
	return FALSE;
}

/**
 * Clear any existing search highlighting
 */
static void clear_search_highlighting(SearchData *search_data) {
    if (!search_data || !search_data->buffer || !search_data->search_tag)
        return;

    GtkTextBuffer *buffer = GTK_TEXT_BUFFER(search_data->buffer);
    GtkTextIter start, end;

    gtk_text_buffer_get_bounds(buffer, &start, &end);
    gtk_text_buffer_remove_tag(buffer, search_data->search_tag, &start, &end);
}

static void update_search_info(SearchData *search_data) {
	if (!search_data || !search_data->info_label)
		return;

	if (search_data->total_matches > 0) {
		gchar *info_text = g_strdup_printf("%d/%d", search_data->current_match, search_data->total_matches);
		gtk_label_set_text(search_data->info_label, info_text);
		g_free(info_text);
	} else {
		gtk_label_set_text(search_data->info_label, "0/0");
	}
}

static void free_match_position(MatchPosition *pos) {
	if (pos) {
		if (pos->start_mark)
			gtk_text_buffer_delete_mark(gtk_text_mark_get_buffer(pos->start_mark), pos->start_mark);
		if (pos->end_mark)
			gtk_text_buffer_delete_mark(gtk_text_mark_get_buffer(pos->end_mark), pos->end_mark);
		g_free(pos);
	}
}

static void clear_match_positions(SearchData *search_data) {
	if (search_data->match_positions) {
		g_slist_free_full(search_data->match_positions, (GDestroyNotify) free_match_position);
		search_data->match_positions = NULL;
	}
}

static void goto_match(SearchData *search_data, gint match_number) {
	if (!search_data || match_number < 1 || match_number > search_data->total_matches)
		return;

	GSList *link = g_slist_nth(search_data->match_positions, match_number - 1);
	if (!link) return;

	MatchPosition *pos = link->data;
	GtkTextIter start_iter, end_iter;

	// Get iters from marks
	gtk_text_buffer_get_iter_at_mark(search_data->buffer, &start_iter, pos->start_mark);
	gtk_text_buffer_get_iter_at_mark(search_data->buffer, &end_iter, pos->end_mark);

	// Select the text
	gtk_text_buffer_select_range(search_data->buffer, &start_iter, &end_iter);

	// Scroll to the match
	gtk_text_view_scroll_to_iter(GTK_TEXT_VIEW(search_data->source_view), &start_iter, 0.0, TRUE, 0.5, 0.5);

	// Update current match number
	search_data->current_match = match_number;
	update_search_info(search_data);
}

static gboolean perform_search(SearchData *search_data) {
	if (!search_data || !search_data->search_entry || !search_data->buffer)
		return FALSE;

	const gchar *search_text = gtk_entry_get_text(search_data->search_entry);
	GtkTextIter iter;
	gboolean found = FALSE;

	// Clear previous highlighting and positions
	clear_search_highlighting(search_data);
	clear_match_positions(search_data);

	search_data->total_matches = 0;
	search_data->current_match = 0;

	if (!search_text || strlen(search_text) == 0) {
		update_search_info(search_data);
		return FALSE;
	}

	// Get buffer bounds
	gtk_text_buffer_get_start_iter(search_data->buffer, &iter);

	// Search and highlight all occurrences
	while (gtk_text_iter_forward_search(&iter, search_text,
			GTK_TEXT_SEARCH_CASE_INSENSITIVE, &search_data->current_match_start,
			&search_data->current_match_end,
			NULL)) {
		search_data->total_matches++;

		// Store match position
		MatchPosition *pos = g_new0(MatchPosition, 1);
		pos->start_mark = gtk_text_buffer_create_mark(search_data->buffer, NULL, &search_data->current_match_start, TRUE);
		pos->end_mark = gtk_text_buffer_create_mark(search_data->buffer, NULL, &search_data->current_match_end, FALSE);
		search_data->match_positions = g_slist_append(search_data->match_positions, pos);

		// Apply highlighting
		gtk_text_buffer_apply_tag(search_data->buffer, search_data->search_tag,
				&search_data->current_match_start,
				&search_data->current_match_end);

		// Move iterator to end of match
		iter = search_data->current_match_end;
		found = TRUE;
	}

	if (found) {
		// Select first match
		search_data->current_match = 1;
		goto_match(search_data, 1);
	}

	update_search_info(search_data);
	return found;
}

void on_next_match(GtkButton *button, SearchData *search_data) {
	if (!search_data || search_data->total_matches == 0)
		return;

	gint next_match = search_data->current_match + 1;
	if (next_match > search_data->total_matches)
		next_match = 1;

	goto_match(search_data, next_match);
}

void on_previous_match(GtkButton *button, SearchData *search_data) {
	if (!search_data || search_data->total_matches == 0)
		return;

	gint prev_match = search_data->current_match - 1;
	if (prev_match < 1)
		prev_match = search_data->total_matches;

	goto_match(search_data, prev_match);
}

void on_find_entry_changed(GtkEntry *entry, SearchData *search_data) {
    perform_search(search_data);
}

void setup_search(GtkSourceView *source_view, GtkEntry *search_entry) {
	if (!source_view || !search_entry)
		return;

	// First, check if we already have search data
	SearchData *search_data = g_object_get_data(G_OBJECT(source_view), "search-data");

	if (search_data) {
		// Update existing search data
		search_data->search_entry = search_entry;
		return;
	}

	// Create new search data only if it doesn't exist
	search_data = g_new0(SearchData, 1);
	search_data->source_view = source_view;
	search_data->buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(source_view));
	search_data->search_entry = search_entry;
	search_data->info_label = find_label;

	// Get the tag table
	GtkTextTagTable *tag_table = gtk_text_buffer_get_tag_table(search_data->buffer);

	// Check if the tag already exists
	search_data->search_tag = gtk_text_tag_table_lookup(tag_table, "search_match");

	// Create the tag only if it doesn't exist
	if (search_data->search_tag == NULL) {
		search_data->search_tag = gtk_text_buffer_create_tag(
				search_data->buffer, "search_match", "background", "yellow",
				"foreground", "black",
				NULL);
	}

	g_signal_connect(search_entry, "changed", G_CALLBACK(on_find_entry_changed), search_data);
	g_signal_connect(go_down_button, "clicked", G_CALLBACK(on_next_match), search_data);
	g_signal_connect(go_up_button, "clicked", G_CALLBACK(on_previous_match), search_data);

	// Store the search data
	g_object_set_data_full(G_OBJECT(source_view), "search-data", search_data, g_free);
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
		syncheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_syntax"));
		bracketcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_bracketmatch"));
		rmargincheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_rmargin"));
		linenumcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_linenums"));
		linemarkcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_linemarks"));
		currentlinecheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_highlightcurrentline"));
		autoindentcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_autoindent"));
		indenttabcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_indentontab"));
		smartbscheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_smartbs"));
		homeendcheck = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_smarthomeend"));
		showspaces = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_showspaces"));
		shownewlines = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_shownewlines"));
		minimap = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_minimap"));
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
		// Findbox
		find_revealer = GTK_REVEALER(gtk_builder_get_object(gui.builder, "find_revealer"));
		find_entry = GTK_ENTRY(lookup_widget("find_entry"));
		find_overlay = GTK_WIDGET(gtk_builder_get_object(gui.builder, "find_overlay"));
		find_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "find_label"));
		go_up_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "go_up_button"));
		go_down_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "go_down_button"));
		// Box
		codeviewbox = GTK_BOX(gtk_builder_get_object(gui.builder, "codeviewbox"));
		add_code_view(gui.builder);
		gtk_window_set_transient_for(GTK_WINDOW(editor_window), GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
		g_signal_connect(sourcebuffer, "modified-changed", G_CALLBACK(on_buffer_modified_changed), NULL);

		/* initialize find box */
		setup_find_overlay();
		setup_search(code_view, find_entry);

		// Initialize with "unsaved" if no file is loaded
		if (!current_file) {
			gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
		}
	}
}

static void update_title_with_modification() {
	const gchar *current_text = gtk_window_get_title(GTK_WINDOW(editor_window));
	if (!current_text || !*current_text) return;

	// If the title already ends with *, don't add another
	if (g_str_has_suffix(current_text, "*")) {
		if (!buffer_modified) {
			gchar *base_name = g_strndup(current_text, strlen(current_text) - 1);
			gtk_window_set_title(GTK_WINDOW(editor_window), base_name);
			g_free(base_name);
		}
	} else if (buffer_modified) {
		gchar *new_title = g_strdup_printf("%s*", current_text);
		gtk_window_set_title(GTK_WINDOW(editor_window), new_title);
		g_free(new_title);
	}
}

static void update_title(GFile *file) {
	if (file) {
		char *basename = g_file_get_basename(file);
		gtk_window_set_title(GTK_WINDOW(editor_window), basename);


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
		gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
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
		gtk_window_set_title(GTK_WINDOW(editor_window), "");
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
	if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;
		gtk_source_buffer_begin_not_undoable_action(sourcebuffer);
		gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		buffer_modified = FALSE;
		gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
		gtk_source_buffer_end_not_undoable_action(sourcebuffer);
		gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
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
		gtk_window_set_title(GTK_WINDOW(editor_window), "");
	}
	gtk_widget_hide(GTK_WIDGET(editor_window));
}

void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
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
}

void on_scratchpad_recent_menu_activated(GtkRecentChooser *chooser, gpointer user_data) {
	gchar *uri;
	GFile *file;
	GError *error = NULL;

	uri = gtk_recent_chooser_get_current_uri(chooser);
	if (!uri) {
		g_warning("Failed to get URI from recent chooser");
		return;
	}

	file = g_file_new_for_uri(uri);
	if (!file) {
		g_warning("Failed to create GFile from URI: %s", uri);
		g_free(uri);
		return;
	}

	if (!g_file_query_exists(file, NULL)) {
		GtkWidget *dialog;
		GtkRecentManager *manager;

		manager = gtk_recent_manager_get_default();

		gtk_recent_manager_remove_item(manager, uri, &error);
		if (error != NULL) {
			g_warning("Failed to remove recent item: %s", error->message);
			g_error_free(error);
		}

		dialog = gtk_message_dialog_new(NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK, "The file '%s' no longer exists.", uri);
		gtk_dialog_run(GTK_DIALOG(dialog));
		gtk_widget_destroy(dialog);

		g_object_unref(file);
		g_free(uri);
		return;
	}

	load_file(file);

	current_file = g_object_ref(file);
	update_title(current_file);

	g_object_unref(file);
	g_free(uri);
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

	// Set the right margin position
	gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);

	// Focus the SourceView
	gtk_widget_grab_focus(GTK_WIDGET(code_view));

	// Initialize the View menu based on com.pref.gui.editor_cfg
	gtk_check_menu_item_set_active(syncheck, com.pref.gui.editor_cfg.highlight_syntax);
	gtk_check_menu_item_set_active(bracketcheck, com.pref.gui.editor_cfg.highlight_bracketmatch);
	gtk_check_menu_item_set_active(rmargincheck, com.pref.gui.editor_cfg.rmargin);
	gtk_check_menu_item_set_active(linenumcheck, com.pref.gui.editor_cfg.show_linenums);
	gtk_check_menu_item_set_active(linemarkcheck, com.pref.gui.editor_cfg.show_linemarks);
	gtk_check_menu_item_set_active(currentlinecheck, com.pref.gui.editor_cfg.highlight_currentline);
	gtk_check_menu_item_set_active(autoindentcheck, com.pref.gui.editor_cfg.autoindent);
	gtk_check_menu_item_set_active(indenttabcheck, com.pref.gui.editor_cfg.indentontab);
	gtk_check_menu_item_set_active(smartbscheck, com.pref.gui.editor_cfg.smartbs);
	gtk_check_menu_item_set_active(homeendcheck, com.pref.gui.editor_cfg.smarthomeend);
	gtk_check_menu_item_set_active(showspaces, com.pref.gui.editor_cfg.showspaces);
	gtk_check_menu_item_set_active(shownewlines, com.pref.gui.editor_cfg.shownewlines);
	gtk_check_menu_item_set_active(minimap, com.pref.gui.editor_cfg.minimap);
	gtk_widget_set_visible(GTK_WIDGET(map), com.pref.gui.editor_cfg.minimap);

	return 0;
}

void on_undo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_source_buffer_can_undo(sourcebuffer)) {
		gtk_source_buffer_undo(sourcebuffer);
	}
}

void on_redo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_source_buffer_can_redo(sourcebuffer)) {
		gtk_source_buffer_redo(sourcebuffer);
	}
}

void on_cut(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view), GDK_SELECTION_CLIPBOARD);
	gtk_text_buffer_cut_clipboard(GTK_TEXT_BUFFER(sourcebuffer),
								clipboard, gtk_text_view_get_editable(GTK_TEXT_VIEW(code_view)));
}

void on_copy(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view), GDK_SELECTION_CLIPBOARD);
	gtk_text_buffer_copy_clipboard(GTK_TEXT_BUFFER(sourcebuffer), clipboard);
}

void on_paste(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view), GDK_SELECTION_CLIPBOARD);
	gtk_text_buffer_paste_clipboard(GTK_TEXT_BUFFER(sourcebuffer),
								clipboard, NULL, gtk_text_view_get_editable(GTK_TEXT_VIEW(code_view)));
}

void on_find(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_find_overlay(TRUE);
}

void on_set_rmarginpos(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *dialog = gtk_dialog_new_with_buttons(
		_("Right Margin Position"),
		NULL,  // parent window could be passed via user_data if needed
		GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
		_("Close"),
		GTK_RESPONSE_CLOSE,
		_("Apply"),
		GTK_RESPONSE_APPLY,
		NULL
	);

	// Set up keyboard shortcuts
	GtkWidget *button_close = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_CLOSE);
	GtkWidget *button_apply = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);

	// Add suggested-action style to Apply button
	gtk_style_context_add_class(gtk_widget_get_style_context(button_apply), "suggested-action");

	gtk_widget_add_accelerator(button_close, "clicked", gtk_accel_group_new(), GDK_KEY_Escape, 0, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(button_apply, "clicked", gtk_accel_group_new(), GDK_KEY_Return, 0, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(button_apply, "clicked", gtk_accel_group_new(), GDK_KEY_KP_Enter, 0, GTK_ACCEL_VISIBLE);

	// Create a horizontal box with spacing
	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), 12);

	// Create and add the label
	GtkWidget *label = gtk_label_new_with_mnemonic(_("Right _margin position"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

	// Create a spin button
	gint current_pos = gtk_source_view_get_right_margin_position(code_view);
	GtkWidget *spin = gtk_spin_button_new_with_range(20, 200, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), current_pos);

	// Connect the label to the spin button for keyboard navigation
	gtk_label_set_mnemonic_widget(GTK_LABEL(label), spin);

	// Add spin button to the box
	gtk_box_pack_start(GTK_BOX(hbox), spin, TRUE, TRUE, 0);

	// Connect to the spin button's activate signal (Enter key)
	g_signal_connect_swapped(spin, "activate", G_CALLBACK(gtk_button_clicked), button_apply);

	// Add the box to the dialog's content area
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), hbox);
	gtk_widget_show_all(hbox);

	// Run the dialog and handle response
	gint response = gtk_dialog_run(GTK_DIALOG(dialog));
	if (response == GTK_RESPONSE_APPLY) {
		com.pref.gui.editor_cfg.rmargin_pos = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin));
		gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}

	gtk_widget_destroy(dialog);
}

// Handler for main window state changes
gboolean on_main_window_state_changed(GtkWidget *widget, GdkEventWindowState *event, gpointer user_data) {
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
	g_signal_connect(editor_window, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

	// Optional: Make editor window minimize when main window minimizes
	g_signal_connect(main_window, "window-state-event", G_CALLBACK(on_main_window_state_changed), editor_window);
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
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_buffer_set_highlight_syntax(sourcebuffer, status);
	com.pref.gui.editor_cfg.highlight_syntax = status;
}

void on_editor_bracketmatch_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_buffer_set_highlight_matching_brackets(sourcebuffer, status);
	com.pref.gui.editor_cfg.highlight_bracketmatch = status;
}

void on_editor_rmargin_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_show_right_margin(code_view, status);
	com.pref.gui.editor_cfg.rmargin = status;
}

void on_editor_linenums_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_show_line_numbers(code_view, status);
	com.pref.gui.editor_cfg.show_linenums = status;
}

void on_editor_linemarks_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_show_line_marks(code_view, status);
	com.pref.gui.editor_cfg.show_linemarks = status;
}

void on_editor_highlightcurrentline_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_highlight_current_line(code_view, status);
	com.pref.gui.editor_cfg.highlight_currentline = status;
}

void on_editor_autoindent_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_auto_indent(code_view, status);
	com.pref.gui.editor_cfg.autoindent = status;
}

void on_editor_indentontab_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_indent_on_tab(code_view, status);
	com.pref.gui.editor_cfg.indentontab = status;
}

void on_editor_smartbs_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_smart_backspace(code_view, status);
	com.pref.gui.editor_cfg.smartbs = status;
}

void on_editor_smarthomeend_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_source_view_set_smart_home_end(code_view, status);
	com.pref.gui.editor_cfg.smarthomeend = status;
}

void on_editor_showspaces_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;
	GtkSourceSpaceTypeFlags spaces_and_tabs = GTK_SOURCE_SPACE_TYPE_SPACE | GTK_SOURCE_SPACE_TYPE_TAB;
	GtkSourceSpaceLocationFlags space_types = gtk_source_space_drawer_get_types_for_locations(space_drawer, space_locations);
	gboolean status = gtk_check_menu_item_get_active(item);
	if (status)
		space_types = space_types | spaces_and_tabs;
	else
		space_types = space_types & ~spaces_and_tabs;

	gtk_source_space_drawer_set_types_for_locations(space_drawer,
												space_locations,
												space_types);
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	com.pref.gui.editor_cfg.showspaces = status;
}

void on_editor_shownewlines_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;
	GtkSourceSpaceLocationFlags space_types = gtk_source_space_drawer_get_types_for_locations(space_drawer, space_locations);
	gboolean status = gtk_check_menu_item_get_active(item);
	if (status)
		space_types = space_types | GTK_SOURCE_SPACE_TYPE_NEWLINE;
	else
		space_types = space_types & ~GTK_SOURCE_SPACE_TYPE_NEWLINE;

	gtk_source_space_drawer_set_types_for_locations(space_drawer,
												space_locations,
												space_types);
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	com.pref.gui.editor_cfg.shownewlines = status;
}

void on_editor_minimap_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	gtk_widget_set_visible(GTK_WIDGET(map), status);
	com.pref.gui.editor_cfg.minimap = status;
}
