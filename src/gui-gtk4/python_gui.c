// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <gtksourceview/gtksource.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/command_line_processor.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/documentation.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/script_menu.h"
#include "gui-gtk4/utils.h"
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
static GtkDropDown *combo_language = NULL;
static GtkSourceLanguageManager *language_manager = NULL;
static GtkSourceLanguage *language = NULL;
static GtkSourceStyleSchemeManager *stylemanager = NULL;
static GtkSourceStyleScheme *scheme = NULL;
static GFile *current_file = NULL;
static GtkWindow *editor_window = NULL;
static GtkWindow *main_window = NULL;
static GtkCheckButton *radio_py = NULL;
static GtkCheckButton *radio_ssf = NULL;
static GtkCheckButton *syncheck = NULL;
static GtkCheckButton *bracketcheck = NULL;
static GtkCheckButton *rmargincheck = NULL;
static GtkCheckButton *linenumcheck = NULL;
static GtkCheckButton *linemarkcheck = NULL;
static GtkCheckButton *currentlinecheck = NULL;
static GtkCheckButton *autoindentcheck = NULL;
static GtkCheckButton *indenttabcheck = NULL;
static GtkCheckButton *smartbscheck = NULL;
static GtkCheckButton *homeendcheck = NULL;
static GtkCheckButton *showspaces = NULL;
static GtkCheckButton *shownewlines = NULL;
static GtkCheckButton *minimap = NULL;
static GtkCheckButton *useargs = NULL;
static GtkSourceSpaceDrawer* space_drawer = NULL;
static GtkBox *codeviewbox = NULL;
static GtkBox *args_box = NULL;
static GtkEntry *args_entry = NULL;

GtkRevealer *find_revealer = NULL;
/* GTK4: GtkSearchEntry no longer derives from GtkEntry — both implement
 * GtkEditable instead.  Store the typed pointer to match the .ui class. */
GtkSearchEntry *find_entry = NULL;
GtkWidget *find_overlay = NULL;
static GtkButton *go_up_button = NULL;
static GtkButton *go_down_button = NULL;

static gint active_language = LANG_PYTHON;
static gboolean buffer_modified = FALSE;
static gboolean from_cli = FALSE;
static gboolean python_debug = FALSE;

// Forward declarations
void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_execute(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_new(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void on_action_file_reload(GSimpleAction *action, GVariant *parameter, gpointer user_data);
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
gboolean on_editor_key_press(GtkEventControllerKey *ctrl, guint keyval, guint keycode,
                              GdkModifierType state, gpointer user_data);

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
	{ "find", on_find, NULL, NULL, NULL },
	{ "reload", on_action_file_reload, NULL, NULL, NULL }
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

gboolean script_editor_has_unsaved_changes() {
	return buffer_modified;
}

void add_code_view(GtkBuilder *builder) {
	// Create a new GtkSourceBuffer (glade doesn't handle this properly)
	sourcebuffer = gtk_source_buffer_new(NULL);
	gtk_text_view_set_buffer(GTK_TEXT_VIEW(code_view), GTK_TEXT_BUFFER(sourcebuffer));
	/* initialize minimap */
	map = GTK_SOURCE_MAP(gtk_source_map_new());

	/* Apply a tiny font to the map's text view so the minimap looks like
	 * a thumbnail rather than a second editor.  The previous CSS used
	 * `* { font-size: 1px }` which also matched GtkSourceMapSlider and
	 * gave it font-driven sizing of zero — every snapshot during scroll
	 * then warned "Trying to snapshot GtkSourceMapSlider without a
	 * current allocation".  Targeting `textview` instead leaves the
	 * slider node at its natural size. */
	siril_register_css_for_display("siril-pyeditor-map",
		".pyeditor-map textview { font-size: 1px; }");
	gtk_widget_add_css_class(GTK_WIDGET(map), "pyeditor-map");

	gtk_widget_set_visible(GTK_WIDGET(map), TRUE);
	gtk_source_map_set_view(map, code_view);
	/* Map sits to the right of the scrolled view in codeviewbox; let it
	 * fill vertically so the slider gets a real allocation and only take
	 * its natural width horizontally. */
	gtk_widget_set_valign(GTK_WIDGET(map), GTK_ALIGN_FILL);
	gtk_box_append(GTK_BOX(codeviewbox), GTK_WIDGET(map));

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


	// Force monospace on the source view via a class-scoped CSS rule.
	siril_register_css_for_display("siril-pyeditor-code",
		".pyeditor-code, .pyeditor-code * { font-family: monospace; }");
	gtk_widget_add_css_class(GTK_WIDGET(code_view), "pyeditor-code");

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
	gtk_widget_set_visible(new_overlay, TRUE);
	gtk_widget_set_vexpand(new_overlay, TRUE);
	gtk_widget_set_hexpand(new_overlay, TRUE);

	// Replace scrolled window with overlay in its parent.  The parent is
	// codeviewbox (GtkBox); the previous code mistakenly cast it to a
	// GtkScrolledWindow and called gtk_scrolled_window_set_child on it,
	// which is a no-op for a GtkBox and left the new overlay unparented
	// (so the GtkSourceView never appeared on screen).
	GtkWidget *scrolled_parent = gtk_widget_get_parent(GTK_WIDGET(scrolled_window));
	g_object_ref(scrolled_window);
	gtk_widget_unparent(GTK_WIDGET(scrolled_window));
	if (GTK_IS_BOX(scrolled_parent)) {
		gtk_box_append(GTK_BOX(scrolled_parent), new_overlay);
	} else if (GTK_IS_SCROLLED_WINDOW(scrolled_parent)) {
		gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_parent), new_overlay);
	} else {
		g_warning("setup_find_overlay: unexpected scrolled-window parent type %s",
		          G_OBJECT_TYPE_NAME(scrolled_parent));
	}

	// Add scrolled window to overlay
	gtk_overlay_set_child(GTK_OVERLAY(new_overlay), GTK_WIDGET(scrolled_window));
	gtk_widget_set_vexpand(GTK_WIDGET(scrolled_window), TRUE);
	gtk_widget_set_hexpand(GTK_WIDGET(scrolled_window), TRUE);
	g_object_unref(scrolled_window);

	// Move find_overlay to new overlay
	g_object_ref(find_overlay);
	gtk_overlay_add_overlay(GTK_OVERLAY(new_overlay), GTK_WIDGET(find_overlay));

	// Configure find_overlay position
	gtk_widget_set_halign(GTK_WIDGET(find_overlay), GTK_ALIGN_END);
	gtk_widget_set_valign(GTK_WIDGET(find_overlay), GTK_ALIGN_START);
	gtk_widget_set_margin_top(GTK_WIDGET(find_overlay), 6);
	gtk_widget_set_margin_end(GTK_WIDGET(find_overlay), 6);

	gtk_revealer_set_transition_type(find_revealer, GTK_REVEALER_TRANSITION_TYPE_SLIDE_DOWN);
	gtk_revealer_set_reveal_child(find_revealer, FALSE);

	g_object_unref(find_overlay);
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


void on_find_entry_activate(GtkSearchEntry *entry, gpointer user_data) {
	hide_find_box();
}

static void update_search_info(SearchData *search_data) {
	if (!search_data || !search_data->info_label)
		return;

	if (search_data->total_matches > 0) {
		gchar *info_text = g_strdup_printf("%d/%d", search_data->current_match, search_data->total_matches);
		gtk_label_set_text(search_data->info_label, info_text);
		g_free(info_text);
	} else {
		const char *text = gtk_editable_get_text(GTK_EDITABLE(search_data->search_entry));
		if (text != NULL && text[0] != '\0') {
			gtk_widget_add_css_class(GTK_WIDGET(find_entry), "search_empty");
		}
		gtk_label_set_text(search_data->info_label, "0/0");
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

static void forward_search(SearchData *search_data) {
	if (!search_data || search_data->total_matches == 0)
		return;

	gint next_match = search_data->current_match + 1;
	if (next_match > search_data->total_matches)
		next_match = 1;

	goto_match(search_data, next_match);
}

static void backward_search(SearchData *search_data) {
	if (!search_data || search_data->total_matches == 0)
		return;

	gint prev_match = search_data->current_match - 1;
	if (prev_match < 1)
		prev_match = search_data->total_matches;

	goto_match(search_data, prev_match);
}

gboolean on_find_entry_key_press_event(GtkEventControllerKey *ctrl, guint keyval,
                                         guint keycode, GdkModifierType state, gpointer data) {
	SearchData *search_data = (SearchData *)data;
	/* Close window */
	if (keyval == GDK_KEY_Tab || keyval == GDK_KEY_Escape) {
		hide_find_box();
		gtk_widget_grab_focus(GTK_WIDGET(code_view));

		return TRUE;
	}

	/* select previous matching iter */
	if (keyval == GDK_KEY_Up || keyval == GDK_KEY_KP_Up) {
		backward_search(search_data);
		return TRUE;
	}

	/* select next matching iter */
	if (keyval == GDK_KEY_Down || keyval == GDK_KEY_KP_Down) {
		forward_search(search_data);
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

static gboolean perform_search(SearchData *search_data) {
	if (!search_data || !search_data->search_entry || !search_data->buffer)
		return FALSE;

	const gchar *search_text = gtk_editable_get_text(GTK_EDITABLE(search_data->search_entry));
	GtkTextIter iter;
	gboolean found = FALSE;
	gtk_widget_remove_css_class(GTK_WIDGET(find_entry), "search_empty");

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
	forward_search(search_data);
}

void on_previous_match(GtkButton *button, SearchData *search_data) {
	backward_search(search_data);
}

void on_find_entry_changed(GtkSearchEntry *entry, SearchData *search_data) {
    (void) entry;
    perform_search(search_data);
}

void setup_search(GtkSourceView *source_view, GtkSearchEntry *search_entry) {
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
	{
		GtkEventController *kctrl = gtk_event_controller_key_new();
		g_signal_connect(kctrl, "key-pressed", G_CALLBACK(on_find_entry_key_press_event), search_data);
		gtk_widget_add_controller(GTK_WIDGET(search_entry), kctrl);
	}
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
		{
			GtkEventController *kctrl = gtk_event_controller_key_new();
			g_signal_connect(kctrl, "key-pressed", G_CALLBACK(on_editor_key_press), NULL);
			gtk_widget_add_controller(GTK_WIDGET(editor_window), kctrl);
		}
		main_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
		// GtkSourceView
		code_view = GTK_SOURCE_VIEW(gtk_builder_get_object(gui.builder, "code_view"));
		sourcebuffer = GTK_SOURCE_BUFFER(gtk_builder_get_object(gui.builder, "sourcebuffer"));
		// GtkCheckMenuItem
		radio_py = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radio_py"));
		radio_ssf = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radio_ssf"));
		syncheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_syntax"));
		bracketcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_bracketmatch"));
		rmargincheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_rmargin"));
		linenumcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_linenums"));
		linemarkcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_linemarks"));
		currentlinecheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_highlightcurrentline"));
		autoindentcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_autoindent"));
		indenttabcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_indentontab"));
		smartbscheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_smartbs"));
		homeendcheck = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_smarthomeend"));
		showspaces = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_showspaces"));
		shownewlines = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_shownewlines"));
		minimap = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_minimap"));
		useargs = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_useargs"));
		// GtkDropDown
		combo_language = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "python_pad_language"));
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
		find_entry = GTK_SEARCH_ENTRY(gtk_builder_get_object(gui.builder, "find_entry"));
		find_overlay = GTK_WIDGET(gtk_builder_get_object(gui.builder, "find_overlay"));
		find_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "find_label"));
		go_up_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "go_up_button"));
		go_down_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "go_down_button"));
		// Box
		codeviewbox = GTK_BOX(gtk_builder_get_object(gui.builder, "codeviewbox"));
		args_box = GTK_BOX(gtk_builder_get_object(gui.builder, "editor_args_box"));
		// GtkEntry
		args_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "editor_args_entry"));

		// GtkSourceView setup
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
				siril_toggle_set_active(GTK_WIDGET(radio_py), TRUE);
			} else if (g_ascii_strcasecmp(suffix, ".ssf") == 0) {
				siril_toggle_set_active(GTK_WIDGET(radio_ssf), TRUE);
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
		gtk_text_buffer_begin_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
		gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		buffer_modified = FALSE;
		gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
		gtk_text_buffer_end_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
		gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}
}

gboolean on_editor_key_press(GtkEventControllerKey *ctrl, guint keyval, guint keycode,
                              GdkModifierType state, gpointer user_data) {
	if ((state & GDK_CONTROL_MASK) && (state & GDK_SHIFT_MASK) && keyval == GDK_KEY_R) {
		on_action_file_reload(NULL, NULL, NULL);
		return TRUE;
	}
	return FALSE;
}

void on_action_file_reload(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_file) {
		siril_log_message(_("No file is currently loaded to reload.\n"));
		return;
	}

	if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will replace the entry buffer with the last saved version of the file. You will not be able to recover any contents."), _("Proceed"))) {
		load_file(current_file);
		update_title(current_file);
	}
}

void new_script(const gchar *contents, gint length, const char *ext) {
	on_open_pythonpad(NULL, NULL);
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;
		gtk_text_buffer_begin_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
		gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sourcebuffer), contents, (gint) length);
		buffer_modified = FALSE;
		gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
		gtk_text_buffer_end_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
		gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
		if (!g_ascii_strcasecmp(ext, "ssf")) {
			active_language = LANG_SSF;
			siril_toggle_set_active(GTK_WIDGET(radio_ssf), TRUE);
		} else {
			active_language = LANG_PYTHON;
			siril_toggle_set_active(GTK_WIDGET(radio_py), TRUE);
		}
		gtk_window_present(editor_window);
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}
}

void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	gint char_count = gtk_text_buffer_get_char_count(GTK_TEXT_BUFFER(sourcebuffer));
	gboolean is_empty = (char_count == 0);
	if (!is_empty) {
		GtkTextIter start, end;
		gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
		if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will clear the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
			if (G_IS_OBJECT(current_file))
				g_object_unref(current_file);
			current_file = NULL;
			gtk_text_buffer_begin_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
			gtk_text_buffer_delete(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
			buffer_modified = FALSE;
			gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
			gtk_text_buffer_end_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
			gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
			gtk_widget_queue_draw(GTK_WIDGET(editor_window));
			gtk_widget_set_visible(GTK_WIDGET(editor_window), FALSE);
			return;
		} else {
			return;
		}
	} else {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;
		gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
	}
	gtk_widget_set_visible(GTK_WIDGET(editor_window), FALSE);
}

/* Phase 14G.2: GtkFileChooserDialog → GtkFileDialog. */
static void on_file_open_done(GObject *src, GAsyncResult *res, gpointer ud) {
	(void)ud;
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *error = NULL;
	GFile *file = gtk_file_dialog_open_finish(fd, res, &error);
	if (file) {
		control_window_switch_to_tab(OUTPUT_LOGS);
		load_file(file);
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = g_object_ref(file);
		update_title(current_file);
		g_object_unref(file);
	}
	if (error)
		g_clear_error(&error);
	gtk_window_present(editor_window);
}

void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	(void)action; (void)parameter; (void)user_data;
	if (!buffer_modified || siril_confirm_dialog(_("Are you sure?"), _("This will replace the entry buffer. You will not be able to recover any contents."), _("Proceed"))) {
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = NULL;

		GtkFileDialog *fd = gtk_file_dialog_new();
		gtk_file_dialog_set_title(fd, _("Open Script"));
		GtkFileFilter *filter = gtk_file_filter_new();
		gtk_file_filter_add_pattern(filter, "*.py");
		gtk_file_filter_add_pattern(filter, "*.ssf");
		gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
		GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
		g_list_store_append(filters, filter);
		gtk_file_dialog_set_default_filter(fd, filter);
		gtk_file_dialog_set_filters(fd, G_LIST_MODEL(filters));
		g_object_unref(filter);
		g_object_unref(filters);

		gtk_file_dialog_open(fd,
			GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")),
			NULL, on_file_open_done, NULL);
		g_object_unref(fd);
	}
}

/* Phase 14: GtkRecentChooser / GtkRecentManager removed in GTK4.
 * Stubbed out — full recent-files reimplementation lands in Phase 15
 * via GtkApplication recent-files. */
void on_scratchpad_recent_menu_activated(GObject *chooser, gpointer user_data) {
	(void)chooser; (void)user_data;
	/* TODO Phase 15: route through GtkApplication / GMenuModel recent. */
}

/* Phase 14G.2: GtkFileChooserDialog → GtkFileDialog. */
static void on_file_save_done(GObject *src, GAsyncResult *res, gpointer ud) {
	(void)ud;
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *error = NULL;
	GFile *file = gtk_file_dialog_save_finish(fd, res, &error);
	if (file) {
		control_window_switch_to_tab(OUTPUT_LOGS);
		if (G_IS_OBJECT(current_file))
			g_object_unref(current_file);
		current_file = g_object_ref(file);
		update_title(current_file);
		save_file(file);
		g_object_unref(file);
	}
	if (error)
		g_clear_error(&error);
	gtk_window_present(editor_window);
}

void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	(void)action; (void)parameter; (void)user_data;
	GtkFileDialog *fd = gtk_file_dialog_new();
	gtk_file_dialog_set_title(fd, _("Save Script As"));

	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.py");
	gtk_file_filter_add_pattern(filter, "*.ssf");
	gtk_file_filter_set_name(filter, _("Script Files (*.py, *.ssf)"));
	GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
	g_list_store_append(filters, filter);
	gtk_file_dialog_set_default_filter(fd, filter);
	gtk_file_dialog_set_filters(fd, G_LIST_MODEL(filters));
	g_object_unref(filter);
	g_object_unref(filters);

	if (G_IS_OBJECT(current_file))
		gtk_file_dialog_set_initial_file(fd, current_file);

	gtk_file_dialog_save(fd,
		GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")),
		NULL, on_file_save_done, NULL);
	g_object_unref(fd);
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

int on_open_pythonpad(GtkWidget *menuitem, gpointer user_data) {
	if (!editor_window) {
		python_scratchpad_init_statics();
	}

	// This is a normal window with normal decorations
	/* Phase 14: gtk_window_set_type_hint removed in GTK4; window-type
	 * hints were window-manager-specific X11 metadata that GTK4 dropped. */

	// Decouple window behaviour from main window
	gtk_window_set_transient_for(editor_window, NULL);
	// (Note, if the editor window is minimized and there isn't an icon in the system tray
	// it can be restored just by clicking on the script editor menu item again - no work
	// will be lost...)

	// Hide on close
	g_signal_connect(editor_window, "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	setup_editor_actions(editor_window);

	/* Phase 17.5: gtk_window_get_position / gtk_window_get_size /
	 * gtk_window_move are all gone in GTK4 — placement is
	 * compositor-managed.  We drop the cascading-offset logic and
	 * rely on the WM to place the editor window. */
	(void)main_window;

	gtk_label_set_text(language_label, active_language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));

	// Show the window and bring it to front
	gtk_window_present(editor_window);

	// Set the correct check menu item active
	siril_toggle_set_active(GTK_WIDGET(radio_py), active_language == LANG_PYTHON);
	siril_toggle_set_active(GTK_WIDGET(radio_ssf), active_language == LANG_SSF);

	// Set the right margin position
	gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);

	// Focus the SourceView
	gtk_widget_grab_focus(GTK_WIDGET(code_view));

	// Initialize the View menu based on com.pref.gui.editor_cfg
	siril_toggle_set_active(GTK_WIDGET(syncheck), com.pref.gui.editor_cfg.highlight_syntax);
	siril_toggle_set_active(GTK_WIDGET(bracketcheck), com.pref.gui.editor_cfg.highlight_bracketmatch);
	siril_toggle_set_active(GTK_WIDGET(rmargincheck), com.pref.gui.editor_cfg.rmargin);
	siril_toggle_set_active(GTK_WIDGET(linenumcheck), com.pref.gui.editor_cfg.show_linenums);
	siril_toggle_set_active(GTK_WIDGET(linemarkcheck), com.pref.gui.editor_cfg.show_linemarks);
	siril_toggle_set_active(GTK_WIDGET(currentlinecheck), com.pref.gui.editor_cfg.highlight_currentline);
	siril_toggle_set_active(GTK_WIDGET(autoindentcheck), com.pref.gui.editor_cfg.autoindent);
	siril_toggle_set_active(GTK_WIDGET(indenttabcheck), com.pref.gui.editor_cfg.indentontab);
	siril_toggle_set_active(GTK_WIDGET(smartbscheck), com.pref.gui.editor_cfg.smartbs);
	siril_toggle_set_active(GTK_WIDGET(homeendcheck), com.pref.gui.editor_cfg.smarthomeend);
	siril_toggle_set_active(GTK_WIDGET(showspaces), com.pref.gui.editor_cfg.showspaces);
	siril_toggle_set_active(GTK_WIDGET(shownewlines), com.pref.gui.editor_cfg.shownewlines);
	siril_toggle_set_active(GTK_WIDGET(minimap), com.pref.gui.editor_cfg.minimap);
	gtk_widget_set_visible(GTK_WIDGET(map), com.pref.gui.editor_cfg.minimap);
	gtk_widget_set_visible(GTK_WIDGET(args_box), siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(useargs))));
	return 0;
}

void on_undo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_text_buffer_get_can_undo(GTK_TEXT_BUFFER(sourcebuffer))) {
		gtk_text_buffer_undo(GTK_TEXT_BUFFER(sourcebuffer));
	}
}

void on_redo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_text_buffer_get_can_redo(GTK_TEXT_BUFFER(sourcebuffer))) {
		gtk_text_buffer_redo(GTK_TEXT_BUFFER(sourcebuffer));
	}
}

void on_cut(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkTextBuffer *buffer = GTK_TEXT_BUFFER(sourcebuffer);
	GdkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view));
	GtkTextIter start, end;

	// Check if there's a selection
	if (!gtk_text_buffer_get_selection_bounds(buffer, &start, &end)) {
		// No selection, so select the entire line
		gtk_text_buffer_get_iter_at_mark(buffer, &start, gtk_text_buffer_get_insert(buffer));
		gtk_text_iter_set_line_offset(&start, 0);

		gtk_text_buffer_get_iter_at_mark(buffer, &end, gtk_text_buffer_get_insert(buffer));
		gtk_text_iter_forward_to_line_end(&end);
		// Move end to the start of the next line to include the newline
		if (!gtk_text_iter_is_end(&end)) {
			gtk_text_iter_forward_line(&end);
		}
		gtk_text_buffer_select_range(buffer, &start, &end);
	}

	// Cut the selected text (either the original selection or the entire line)
	gtk_text_buffer_cut_clipboard(buffer, clipboard, gtk_text_view_get_editable(GTK_TEXT_VIEW(code_view)));
}

void on_copy(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GdkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view));
	gtk_text_buffer_copy_clipboard(GTK_TEXT_BUFFER(sourcebuffer), clipboard);
}

void on_paste(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_revealer_get_reveal_child(GTK_REVEALER(find_revealer))) {
		gtk_widget_activate_action(GTK_WIDGET(find_entry), "clipboard.paste", NULL);
	} else {
		GdkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(code_view));
		gtk_text_buffer_paste_clipboard(GTK_TEXT_BUFFER(sourcebuffer), clipboard, NULL, gtk_text_view_get_editable(GTK_TEXT_VIEW(code_view)));
	}
}

void on_find(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_find_overlay(TRUE);
}

/* Phase 14G.4: helper for the rmargin dialog — Close/Apply buttons feed
 * a small state struct that shuttles the response code out of the
 * nested GMainLoop. */
struct rmargin_dlg_ctx { GMainLoop *loop; gint result; };
static void rmargin_close_clicked(GtkButton *b, gpointer ud) { (void)b; struct rmargin_dlg_ctx *c = ud; c->result = GTK_RESPONSE_CLOSE; if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); }
static void rmargin_apply_clicked(GtkButton *b, gpointer ud) { (void)b; struct rmargin_dlg_ctx *c = ud; c->result = GTK_RESPONSE_APPLY; if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); }
static gboolean rmargin_close_request(GtkWindow *w, gpointer ud) { (void)w; struct rmargin_dlg_ctx *c = ud; c->result = GTK_RESPONSE_CLOSE; if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); return FALSE; }

static gboolean rmargin_key_pressed(GtkEventControllerKey *ctrl, guint keyval,
		guint keycode, GdkModifierType state, gpointer ud) {
	(void)ctrl; (void)keycode; (void)state;
	struct rmargin_dlg_ctx *c = ud;
	if (keyval == GDK_KEY_Escape) {
		c->result = GTK_RESPONSE_CLOSE;
		if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop);
		return TRUE;
	}
	if (keyval == GDK_KEY_Return || keyval == GDK_KEY_KP_Enter) {
		c->result = GTK_RESPONSE_APPLY;
		if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop);
		return TRUE;
	}
	return FALSE;
}

void on_set_rmarginpos(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	(void)action; (void)parameter; (void)user_data;
	/* Phase 14G.4: GtkDialog → GtkWindow.  Escape/Enter shortcuts now
	 * arrive via a GtkEventControllerKey because GTK4 removed
	 * gtk_widget_add_accelerator / gtk_accel_group_new. */
	GtkWidget *dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Right Margin Position"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(dialog), TRUE);

	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);

	// Create a horizontal box with the label + spinner
	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	GtkWidget *label = gtk_label_new_with_mnemonic(_("Right _margin position"));
	gtk_widget_set_halign(label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(label, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(hbox), label);

	gint current_pos = gtk_source_view_get_right_margin_position(code_view);
	GtkWidget *spin = gtk_spin_button_new_with_range(20, 200, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), current_pos);
	gtk_label_set_mnemonic_widget(GTK_LABEL(label), spin);
	gtk_widget_set_hexpand(spin, TRUE);
	gtk_widget_set_vexpand(spin, TRUE);
	gtk_box_append(GTK_BOX(hbox), spin);

	gtk_box_append(GTK_BOX(content_area), hbox);

	// Action area
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *button_close = gtk_button_new_with_label(_("Close"));
	GtkWidget *button_apply = gtk_button_new_with_label(_("Apply"));
	gtk_widget_add_css_class(button_apply, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), button_close);
	gtk_box_append(GTK_BOX(bbox), button_apply);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), button_apply);
	gtk_window_set_child(GTK_WINDOW(dialog), content_area);

	struct rmargin_dlg_ctx ctx = { g_main_loop_new(NULL, FALSE), GTK_RESPONSE_CLOSE };
	g_signal_connect(button_close, "clicked",       G_CALLBACK(rmargin_close_clicked), &ctx);
	g_signal_connect(button_apply, "clicked",       G_CALLBACK(rmargin_apply_clicked), &ctx);
	g_signal_connect(dialog,       "close-request", G_CALLBACK(rmargin_close_request), &ctx);

	GtkEventController *kctrl = gtk_event_controller_key_new();
	g_signal_connect(kctrl, "key-pressed", G_CALLBACK(rmargin_key_pressed), &ctx);
	gtk_widget_add_controller(dialog, kctrl);

	gtk_window_present(GTK_WINDOW(dialog));
	g_main_loop_run(ctx.loop);
	g_main_loop_unref(ctx.loop);

	if (ctx.result == GTK_RESPONSE_APPLY) {
		com.pref.gui.editor_cfg.rmargin_pos = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin));
		gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}

	gtk_window_destroy(GTK_WINDOW(dialog));
}

/* TODO Phase8: window-state-event is gone in GTK4.  Use GtkWindow's
 * "notify::minimized" / "notify::maximized" / "notify::fullscreened"
 * properties to track the equivalent state.  Stubbed for now. */
gboolean on_main_window_state_changed(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	(void)widget; (void)event; (void)user_data;
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

void on_language_set(GtkCheckButton *item, gpointer user_data) {
	if (siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(radio_py)))) {
		active_language = LANG_PYTHON;
	} else {
		active_language = LANG_SSF;
	}
	set_language();
}

void setup_python_editor_window() {
	GtkWindow *editor_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "python_window"));
	GtkWindow *main_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));

	// Make the editor window hide instead of destroy when closed
	g_signal_connect(editor_window, "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	// Optional: Make editor window minimize when main window minimizes
	/* GTK4: window-state-event removed; use the "minimized" property */
	g_signal_connect(main_window, "notify::minimized",
	                 G_CALLBACK(on_main_window_state_changed), editor_window);
}

void on_action_file_execute(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!accept_script_warning_dialog())
		return;
	gchar** script_args = NULL;
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(sourcebuffer), &start, &end);
	// Get the text
	char *text = gtk_text_buffer_get_text(GTK_TEXT_BUFFER(sourcebuffer), &start, &end, FALSE);
	switch (active_language) {
		case LANG_PYTHON:;
			// Create a temporary file for the script
			GError *error = NULL;
			gchar *temp_filename = NULL;
			int fd = g_file_open_tmp("siril-script-XXXXXX.py", &temp_filename, &error);
			if (fd == -1) {
				siril_log_message(_("Error creating temporary script file: %s\n"), error->message);
				g_error_free(error);
				g_free(text);
				return;
			}
			// Write script content to the temporary file
			if (write(fd, text, strlen(text)) == -1) {
				siril_log_message(_("Error writing to temporary script file\n"));
				close(fd);
				if (g_unlink(temp_filename))
					siril_debug_print("g_unlink() failed in on_action_file_execute()\n");
				g_free(temp_filename);
				g_free(text);
				return;
			}
			close(fd);
			// Add args if we need to
			if (siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(useargs)))) {
				const gchar *args_string = gtk_editable_get_text(GTK_EDITABLE(args_entry));
				script_args = g_strsplit(args_string, " ", -1);
				from_cli = TRUE;
			}

			// Execute the script with the path to the temp file instead of the text content
			// Passing TRUE as the last parameter to indicate this is a temporary file
			execute_python_script(temp_filename, TRUE, FALSE, script_args, TRUE, from_cli, python_debug);
			g_strfreev(script_args);
			g_free(text);
			break;
		case LANG_SSF:;
			GInputStream *input_stream = g_memory_input_stream_new_from_data(text, strlen(text), NULL);
			g_free(text);
			if (processing_is_job_active()) {
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
			g_free(text);
			break;
	}
}

static void on_buffer_modified_changed(GtkTextBuffer *buffer, gpointer user_data) {
	buffer_modified = gtk_text_buffer_get_modified(buffer);
	update_title_with_modification();
}

void on_action_python_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("Python-API.html");
}

void on_action_command_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("Commands.html");
}

void on_editor_syntax_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_buffer_set_highlight_syntax(sourcebuffer, status);
	com.pref.gui.editor_cfg.highlight_syntax = status;
}

void on_editor_bracketmatch_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_buffer_set_highlight_matching_brackets(sourcebuffer, status);
	com.pref.gui.editor_cfg.highlight_bracketmatch = status;
}

void on_editor_rmargin_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_show_right_margin(code_view, status);
	com.pref.gui.editor_cfg.rmargin = status;
}

void on_editor_linenums_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_show_line_numbers(code_view, status);
	com.pref.gui.editor_cfg.show_linenums = status;
}

void on_editor_linemarks_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_show_line_marks(code_view, status);
	com.pref.gui.editor_cfg.show_linemarks = status;
}

void on_editor_highlightcurrentline_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_highlight_current_line(code_view, status);
	com.pref.gui.editor_cfg.highlight_currentline = status;
}

void on_editor_autoindent_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_auto_indent(code_view, status);
	com.pref.gui.editor_cfg.autoindent = status;
}

void on_editor_indentontab_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_indent_on_tab(code_view, status);
	com.pref.gui.editor_cfg.indentontab = status;
}

void on_editor_smartbs_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_smart_backspace(code_view, status);
	com.pref.gui.editor_cfg.smartbs = status;
}

void on_editor_smarthomeend_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_source_view_set_smart_home_end(code_view, status);
	com.pref.gui.editor_cfg.smarthomeend = status;
}

void on_editor_showspaces_toggled(GtkCheckButton *item, gpointer user_data) {
	// Define where spaces should be displayed (e.g., everywhere)
	GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

	// Get the current space types for the specified locations
	GtkSourceSpaceTypeFlags space_types = gtk_source_space_drawer_get_types_for_locations(space_drawer, space_locations);

	// Determine the new status (enabled or disabled)
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));

	// Define the space types to toggle (spaces and tabs)
	GtkSourceSpaceTypeFlags spaces_and_tabs = GTK_SOURCE_SPACE_TYPE_SPACE | GTK_SOURCE_SPACE_TYPE_TAB;

	// Update the space types based on the status
	if (status) {
		space_types |= spaces_and_tabs; // Enable spaces and tabs
	} else {
		space_types &= ~spaces_and_tabs; // Disable spaces and tabs
	}

	// Apply the updated space types to the drawer
	gtk_source_space_drawer_set_types_for_locations(space_drawer, space_locations, space_types);

	// Redraw the editor window to reflect the changes
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));

	// Update the configuration
	com.pref.gui.editor_cfg.showspaces = status;
}

void on_editor_shownewlines_toggled(GtkCheckButton *item, gpointer user_data) {
	// Define where spaces should be displayed (e.g., everywhere)
	GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

	// Get the current space types for the specified locations
	GtkSourceSpaceTypeFlags space_types = gtk_source_space_drawer_get_types_for_locations(space_drawer, space_locations);

	// Determine the new status (enabled or disabled)
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));

	// Update the space types based on the status
	if (status) {
		space_types |= GTK_SOURCE_SPACE_TYPE_NEWLINE; // Enable newlines
	} else {
		space_types &= ~GTK_SOURCE_SPACE_TYPE_NEWLINE; // Disable newlines
	}

	// Apply the updated space types to the drawer
	gtk_source_space_drawer_set_types_for_locations(space_drawer, space_locations, space_types);

	// Redraw the editor window to reflect the changes
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));

	// Update the configuration
	com.pref.gui.editor_cfg.shownewlines = status;
}

void on_editor_minimap_toggled(GtkCheckButton *item, gpointer user_data) {
	gboolean status = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	gtk_widget_set_visible(GTK_WIDGET(map), status);
	com.pref.gui.editor_cfg.minimap = status;
}

void on_editor_useargs_toggled(GtkCheckButton *item, gpointer user_data) {
	gtk_widget_set_visible(GTK_WIDGET(args_box), siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item))));
	if (siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)))) {
		from_cli = TRUE;
	} else {
		from_cli = FALSE;
	}
}

void on_editor_args_clear_clicked(GtkButton *button, gpointer user_data) {
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(args_entry);
	gtk_entry_buffer_set_text(buffer, "", 0);
	gtk_widget_grab_focus(GTK_WIDGET(args_entry));
}

void on_pythondebug_toggled(GtkCheckButton *item, gpointer user_data) {
    gboolean state = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(item)));
	static GtkCheckButton *editorwidget = NULL, *scriptmenuwidget = NULL;
	if (!editorwidget) {
		editorwidget = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "editor_toggledebug"));
		scriptmenuwidget = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "pythondebugtoggle"));
	}

	// Synchronize the two checkbuttons
	g_signal_handlers_block_by_func(editorwidget, on_pythondebug_toggled, NULL);
	g_signal_handlers_block_by_func(scriptmenuwidget, on_pythondebug_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(editorwidget), state);
	siril_toggle_set_active(GTK_WIDGET(scriptmenuwidget), state);
	g_signal_handlers_unblock_by_func(editorwidget, on_pythondebug_toggled, NULL);
	g_signal_handlers_unblock_by_func(scriptmenuwidget, on_pythondebug_toggled, NULL);

	if (state) {
		python_debug = TRUE;  // TRUE to overwrite if it exists
	} else {
		python_debug = FALSE;
	}
}

gboolean get_python_debug_mode() {
	return python_debug;
}
