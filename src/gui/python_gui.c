// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "gui/documentation.h"
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
static GtkNotebook *notebook = NULL;
static GtkComboBox *combo_language = NULL;
static GtkSourceLanguageManager *language_manager = NULL;
static GtkSourceLanguage *language = NULL;
static GtkSourceStyleSchemeManager *stylemanager = NULL;
static GtkSourceStyleScheme *scheme = NULL;
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
static GtkCheckMenuItem *useargs = NULL;
static GtkBox *codeviewbox = NULL;
static GtkBox *args_box = NULL;
static GtkEntry *args_entry = NULL;

GtkRevealer *find_revealer = NULL;
GtkEntry *find_entry = NULL;
GtkWidget *find_overlay = NULL;
static GtkButton *go_up_button = NULL;
static GtkButton *go_down_button = NULL;

static GList *tabs = NULL;
static TabInfo *current_tab = NULL;
static gboolean from_cli = FALSE;
static gboolean python_debug = FALSE;

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
static void on_tab_buffer_modified_changed(GtkTextBuffer *buffer, TabInfo *tab);
void set_language();
void setup_search_for_tab(TabInfo *tab);
static void apply_editor_settings_to_tab(TabInfo *tab);

// Tab management functions
static TabInfo* find_tab_by_widget(GtkWidget *widget);
static void update_current_tab(gint current_page);
static void update_ui_for_current_tab(void);
static void update_tab_title(TabInfo *tab);
static TabInfo* create_new_tab(const gchar *title, const gchar *content);
static void close_tab(TabInfo *tab);
static void on_tab_close_clicked(GtkButton *button, TabInfo *tab);
static void on_notebook_switch_page(GtkNotebook *notebook, GtkWidget *page, guint page_num, gpointer user_data);
static void apply_editor_settings_to_tab(TabInfo *tab);

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
	if (scheme && current_tab) {
		gtk_source_buffer_set_style_scheme(current_tab->source_buffer, scheme);
	}
}

gboolean code_view_exists() {
	return (current_tab != NULL && current_tab->source_view != NULL);
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
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		if (tab->modified) {
			return TRUE;
		}
	}
	return FALSE;
}

static void set_source_view_properties(GtkSourceView *source_view, GtkSourceBuffer *source_buffer) {
	// Set the GtkSourceView style depending on whether the Siril light or dark
	// theme is set.
	if (scheme) {
		gtk_source_buffer_set_style_scheme(source_buffer, scheme);
	}

	// Get the GtkSourceLanguageManager and set the Python language
	if (!language_manager) {
		language_manager = gtk_source_language_manager_get_default();
	}
	language = gtk_source_language_manager_get_language(language_manager, "python3");
	if (language == NULL) {
		siril_debug_print("Could not find Python language definition\n");
	} else {
		gtk_source_buffer_set_language(source_buffer, language);
	}

	// Force monospace
	GtkCssProvider *css = gtk_css_provider_new();
	gtk_css_provider_load_from_data(css, "* { font-family: monospace;}",-1,NULL);
	GtkStyleContext *context = gtk_widget_get_style_context(GTK_WIDGET(source_view));
	gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
	g_object_unref(css);

	// Enable syntax highlighting
	gtk_source_buffer_set_highlight_syntax(source_buffer, TRUE);
}

static TabInfo* create_new_tab(const gchar *title, const gchar *content) {
    TabInfo *tab = g_new0(TabInfo, 1);

    // Create tab content container
    GtkWidget *tab_content = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);

    // Create source view and buffer
    tab->source_buffer = gtk_source_buffer_new(NULL);
    tab->source_view = GTK_SOURCE_VIEW(gtk_source_view_new_with_buffer(tab->source_buffer));

    // Apply styling and properties
    set_source_view_properties(tab->source_view, tab->source_buffer);

    // Create scrolled window for source view FIRST
    GtkWidget *scrolled = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled),
                                  GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(scrolled), GTK_WIDGET(tab->source_view));

    // Create minimap
    tab->minimap = GTK_SOURCE_MAP(gtk_source_map_new());

    // Apply CSS to force font size to 1 for minimap
    GtkCssProvider* provider = gtk_css_provider_new();
    gtk_css_provider_load_from_data(provider,
        "* { font-size: 1px; }",
        -1, NULL);
    gtk_style_context_add_provider(
        gtk_widget_get_style_context(GTK_WIDGET(tab->minimap)),
        GTK_STYLE_PROVIDER(provider),
        GTK_STYLE_PROVIDER_PRIORITY_USER
    );
    g_object_unref(provider);

    gtk_widget_show(GTK_WIDGET(tab->minimap));
    gtk_source_map_set_view(tab->minimap, tab->source_view);

    // Get the space drawer for this source view
    tab->space_drawer = gtk_source_view_get_space_drawer(tab->source_view);

    // Configure which types of spaces to show
    GtkSourceSpaceTypeFlags space_types = 0;
    GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

    // Enable space drawing with marks
    gtk_source_space_drawer_set_types_for_locations(tab->space_drawer,
                                                space_locations,
                                                space_types);
    gtk_source_space_drawer_set_enable_matrix(tab->space_drawer, TRUE);

    // Pack into tab content - FIX: Ensure proper packing order and properties
    gtk_box_pack_start(GTK_BOX(tab_content), scrolled, TRUE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(tab_content), GTK_WIDGET(tab->minimap), FALSE, FALSE, 0);

    // Create tab label with close button
    GtkWidget *tab_label_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
    tab->tab_label = GTK_LABEL(gtk_label_new(title ? title : "untitled"));
    GtkWidget *close_button = gtk_button_new_from_icon_name("window-close-symbolic", GTK_ICON_SIZE_MENU);
    gtk_button_set_relief(GTK_BUTTON(close_button), GTK_RELIEF_NONE);

    gtk_box_pack_start(GTK_BOX(tab_label_box), GTK_WIDGET(tab->tab_label), FALSE, FALSE, 0);
    gtk_box_pack_end(GTK_BOX(tab_label_box), close_button, FALSE, FALSE, 0);
    gtk_widget_show_all(tab_label_box);

    // Add tab to notebook
    gint page_num = gtk_notebook_append_page(notebook, tab_content, tab_label_box);
    tab->tab_widget = tab_content;

    // Set initial content
    if (content) {
        gtk_text_buffer_set_text(GTK_TEXT_BUFFER(tab->source_buffer), content, -1);
    }

    // Initialize other properties
    tab->file = NULL;
    tab->modified = FALSE;
    tab->language = LANG_PYTHON;
    tab->title = g_strdup(title ? title : "untitled");

    // Setup search for this tab
    setup_search_for_tab(tab);

    // Connect signals
    g_signal_connect(tab->source_buffer, "modified-changed",
                     G_CALLBACK(on_tab_buffer_modified_changed), tab);
    g_signal_connect(close_button, "clicked",
                     G_CALLBACK(on_tab_close_clicked), tab);

    // Add to tabs list
    tabs = g_list_append(tabs, tab);

    gtk_widget_show_all(tab_content);

    gtk_notebook_set_current_page(notebook, page_num);

    // Apply editor settings to the new tab
    apply_editor_settings_to_tab(tab);

    return tab;
}

static void close_tab(TabInfo *tab) {
    if (!tab) return;

    // If this is the last tab, don't close it - just clear it
    if (g_list_length(tabs) <= 1) {
        // Clear the content but keep the tab
        GtkTextIter start, end;
        gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(tab->source_buffer), &start, &end);
        if (tab->modified && !siril_confirm_dialog(_("Are you sure?"),
                                                   _("This will clear the entry buffer. You will not be able to recover any contents."),
                                                   _("Proceed"))) {
            return;
        }

        if (tab->file) {
            g_object_unref(tab->file);
            tab->file = NULL;
        }

        gtk_source_buffer_begin_not_undoable_action(tab->source_buffer);
        gtk_text_buffer_delete(GTK_TEXT_BUFFER(tab->source_buffer), &start, &end);
        tab->modified = FALSE;
        gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(tab->source_buffer), FALSE);
        gtk_source_buffer_end_not_undoable_action(tab->source_buffer);

        g_free(tab->title);
        tab->title = g_strdup("untitled");
        update_tab_title(tab);
        return;
    }

    // Check for unsaved changes
    if (tab->modified) {
        if (!siril_confirm_dialog(_("Are you sure?"),
                                  _("This tab has unsaved changes. Close anyway?"),
                                  _("Close"))) {
            return;
        }
    }

    // Remove from notebook
    gint page_num = gtk_notebook_page_num(notebook, tab->tab_widget);
    gtk_notebook_remove_page(notebook, page_num);

    // Clean up tab data
    if (tab->file) g_object_unref(tab->file);
    if (tab->search_data) g_free(tab->search_data);
    g_free(tab->title);

    // Remove from list
    tabs = g_list_remove(tabs, tab);

    // Update current tab
    if (current_tab == tab) {
        current_tab = NULL;
        gint current_page = gtk_notebook_get_current_page(notebook);
        if (current_page >= 0) {
            update_current_tab(current_page);
        }
    }

    g_free(tab);
}

static TabInfo* find_tab_by_widget(GtkWidget *widget) {
    for (GList *l = tabs; l; l = l->next) {
        TabInfo *tab = l->data;
        if (tab->tab_widget == widget) {
            return tab;
        }
    }
    return NULL;
}

static void update_current_tab(gint current_page) {
    if (current_page < 0) {
        current_tab = NULL;
        return;
    }

    GtkWidget *page_widget = gtk_notebook_get_nth_page(notebook, current_page);
    current_tab = find_tab_by_widget(page_widget);

    if (current_tab) {
        // Update UI to reflect current tab state
        update_ui_for_current_tab();
    }
}

static void update_ui_for_current_tab() {
    if (!current_tab) return;

    // Update language label
    gtk_label_set_text(language_label,
                      current_tab->language == LANG_PYTHON ?
                      _("Python Script") : _("Siril Script File"));

    // Update menu items
    gtk_check_menu_item_set_active(radio_py, current_tab->language == LANG_PYTHON);
    gtk_check_menu_item_set_active(radio_ssf, current_tab->language == LANG_SSF);

    // Update style scheme for current tab
    set_code_view_theme();

    // Update minimap visibility based on preferences
    gtk_widget_set_visible(GTK_WIDGET(current_tab->minimap), com.pref.gui.editor_cfg.minimap);

    // Apply all editor settings to ensure consistency
    apply_editor_settings_to_tab(current_tab);

    // Focus the current source view
    gtk_widget_grab_focus(GTK_WIDGET(current_tab->source_view));
}

static void update_tab_title(TabInfo *tab) {
    if (!tab) return;

    gchar *display_title;
    if (tab->file) {
        char *basename = g_file_get_basename(tab->file);
        if (tab->modified) {
            display_title = g_strdup_printf("%s*", basename);
        } else {
            display_title = g_strdup(basename);
        }
        g_free(basename);
    } else {
        if (tab->modified) {
            display_title = g_strdup("untitled*");
        } else {
            display_title = g_strdup("untitled");
        }
    }

    gtk_label_set_text(tab->tab_label, display_title);
    g_free(display_title);
}

static void on_tab_close_clicked(GtkButton *button, TabInfo *tab) {
    close_tab(tab);
}

static void apply_editor_settings_to_tab(TabInfo *tab) {
    if (!tab) return;

    // Apply all the editor preferences to this tab
    gtk_source_buffer_set_highlight_syntax(tab->source_buffer, com.pref.gui.editor_cfg.highlight_syntax);
    gtk_source_buffer_set_highlight_matching_brackets(tab->source_buffer, com.pref.gui.editor_cfg.highlight_bracketmatch);
    gtk_source_view_set_show_right_margin(tab->source_view, com.pref.gui.editor_cfg.rmargin);
    gtk_source_view_set_right_margin_position(tab->source_view, com.pref.gui.editor_cfg.rmargin_pos);
    gtk_source_view_set_show_line_numbers(tab->source_view, com.pref.gui.editor_cfg.show_linenums);
    gtk_source_view_set_show_line_marks(tab->source_view, com.pref.gui.editor_cfg.show_linemarks);
    gtk_source_view_set_highlight_current_line(tab->source_view, com.pref.gui.editor_cfg.highlight_currentline);
    gtk_source_view_set_auto_indent(tab->source_view, com.pref.gui.editor_cfg.autoindent);
    gtk_source_view_set_indent_on_tab(tab->source_view, com.pref.gui.editor_cfg.indentontab);
    gtk_source_view_set_smart_backspace(tab->source_view, com.pref.gui.editor_cfg.smartbs);
    gtk_source_view_set_smart_home_end(tab->source_view, com.pref.gui.editor_cfg.smarthomeend);

    // Apply space drawer settings
    GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;
    GtkSourceSpaceTypeFlags space_types = 0;

    if (com.pref.gui.editor_cfg.showspaces) {
        space_types |= GTK_SOURCE_SPACE_TYPE_SPACE | GTK_SOURCE_SPACE_TYPE_TAB;
    }
    if (com.pref.gui.editor_cfg.shownewlines) {
        space_types |= GTK_SOURCE_SPACE_TYPE_NEWLINE;
    }

    gtk_source_space_drawer_set_types_for_locations(tab->space_drawer, space_locations, space_types);

    // Update minimap visibility
    gtk_widget_set_visible(GTK_WIDGET(tab->minimap), com.pref.gui.editor_cfg.minimap);
}

static void on_notebook_switch_page(GtkNotebook *notebook, GtkWidget *page,
                                   guint page_num, gpointer user_data) {
    update_current_tab(page_num);

    // Additional UI updates when switching tabs
    if (current_tab) {
        // Apply editor settings to the newly active tab
        apply_editor_settings_to_tab(current_tab);
    }
}

void add_code_view(GtkBuilder *builder) {
	// Create notebook for tabs
	notebook = GTK_NOTEBOOK(gtk_notebook_new());
	gtk_notebook_set_tab_pos(notebook, GTK_POS_TOP);
	gtk_notebook_set_scrollable(notebook, TRUE);

	// Add notebook to codeviewbox
	gtk_box_pack_start(GTK_BOX(codeviewbox), GTK_WIDGET(notebook), TRUE, TRUE, 0);

	// Connect notebook signals
	g_signal_connect(notebook, "switch-page", G_CALLBACK(on_notebook_switch_page), NULL);

	// Set up style scheme manager and scheme
	stylemanager = gtk_source_style_scheme_manager_get_default();
	scheme = gtk_source_style_scheme_manager_get_scheme(stylemanager,
		com.pref.gui.combo_theme == 0 ? "oblivion" : "classic");

	// Create initial tab
	current_tab = create_new_tab("untitled", NULL);

	gtk_widget_show_all(GTK_WIDGET(notebook));
}

/** Code for the find feature ***/

static void setup_find_overlay(void) {
    GtkWidget *new_overlay = gtk_overlay_new();
    GtkWidget *target = NULL;
    GtkWidget *parent = NULL;
    GList *children = NULL;

    gtk_widget_show(new_overlay);

    /* Prefer wrapping the notebook (tabbed UI). Otherwise use the
	 * first child of codeviewbox. */
    if (GTK_IS_NOTEBOOK(notebook)) {
        target = GTK_WIDGET(notebook);
    } else if (codeviewbox && GTK_IS_WIDGET(GTK_WIDGET(codeviewbox))) {
        children = gtk_container_get_children(GTK_CONTAINER(codeviewbox));
        if (children) {
            target = GTK_WIDGET(children->data); /* first child */
        }
    }

    if (!target) {
        if (children) g_list_free(children);
        g_warning("setup_find_overlay: no suitable widget found to attach overlay to");
        return;
    }

    parent = gtk_widget_get_parent(target);

    if (parent) {
        /* Keep a ref while we reparent the widget */
        g_object_ref(target);
        gtk_container_remove(GTK_CONTAINER(parent), target);

        /* Put overlay into the parent's place */
        gtk_container_add(GTK_CONTAINER(parent), new_overlay);

        /* Make the original target the main child of the overlay */
        gtk_container_add(GTK_CONTAINER(new_overlay), target);
        gtk_widget_set_vexpand(target, TRUE);
        gtk_widget_set_hexpand(target, TRUE);

        /* release the ref we took earlier */
        g_object_unref(target);
    } else {
        /* No parent: pack overlay into codeviewbox if possible */
        if (codeviewbox && GTK_IS_BOX(codeviewbox)) {
            gtk_container_add(GTK_CONTAINER(new_overlay), target);
            gtk_box_pack_start(GTK_BOX(codeviewbox), new_overlay, TRUE, TRUE, 0);
            gtk_widget_show_all(new_overlay);
        } else {
            if (children) g_list_free(children);
            g_warning("setup_find_overlay: cannot attach overlay (no parent and cannot pack)");
            return;
        }
    }

    if (children) g_list_free(children);

    /* Move the find overlay widget from the builder into our new overlay */
    if (find_overlay) {
        g_object_ref(find_overlay);
        gtk_overlay_add_overlay(GTK_OVERLAY(new_overlay), find_overlay);

        /* position the find overlay (same as original behaviour) */
        gtk_widget_set_halign(GTK_WIDGET(find_overlay), GTK_ALIGN_END);
        gtk_widget_set_valign(GTK_WIDGET(find_overlay), GTK_ALIGN_START);
        gtk_widget_set_margin_top(GTK_WIDGET(find_overlay), 6);
        gtk_widget_set_margin_end(GTK_WIDGET(find_overlay), 6);

        gtk_revealer_set_transition_type(find_revealer, GTK_REVEALER_TRANSITION_TYPE_SLIDE_DOWN);
        gtk_revealer_set_reveal_child(find_revealer, FALSE);

        g_object_unref(find_overlay);
    } else {
        g_warning("setup_find_overlay: find_overlay is NULL");
    }
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

static void update_search_info(SearchData *search_data) {
	if (!search_data || !search_data->info_label)
		return;

	if (search_data->total_matches > 0) {
		gchar *info_text = g_strdup_printf("%d/%d", search_data->current_match, search_data->total_matches);
		gtk_label_set_text(search_data->info_label, info_text);
		g_free(info_text);
	} else {
		const char *text = gtk_entry_get_text(GTK_ENTRY(search_data->search_entry));
		if (text != NULL && text[0] != '\0') {
			GtkStyleContext *context;
			context = gtk_widget_get_style_context(GTK_WIDGET(find_entry));
			gtk_style_context_add_class(context, "search_empty");
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

gboolean on_find_entry_key_press_event(GtkWidget *widget, GdkEventKey *event, gpointer data) {
	SearchData *search_data = (SearchData *)data;
	/* Close window */
	if (event->keyval == GDK_KEY_Tab || event->keyval == GDK_KEY_Escape) {
		hide_find_box();
		if (current_tab) {
			gtk_widget_grab_focus(GTK_WIDGET(current_tab->source_view));
		}
		return TRUE;
	}

	/* select previous matching iter */
	if (event->keyval == GDK_KEY_Up || event->keyval == GDK_KEY_KP_Up) {
		backward_search(search_data);
		return TRUE;
	}

	/* select next matching iter */
	if (event->keyval == GDK_KEY_Down || event->keyval == GDK_KEY_KP_Down) {
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

	const gchar *search_text = gtk_entry_get_text(search_data->search_entry);
	GtkTextIter iter;
	gboolean found = FALSE;
	GtkStyleContext *context;
	context = gtk_widget_get_style_context(GTK_WIDGET(find_entry));
	gtk_style_context_remove_class(context, "search_empty");

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

void on_find_entry_changed(GtkEntry *entry, SearchData *search_data) {
    perform_search(search_data);
}

void setup_search_for_tab(TabInfo *tab) {
	if (!tab) return;

	SearchData *search_data = g_new0(SearchData, 1);
	search_data->source_view = tab->source_view;
	search_data->buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tab->source_view));
	search_data->search_entry = find_entry; // Shared search entry
	search_data->info_label = find_label;   // Shared info label

	// Create search tag for this buffer
	search_data->search_tag = gtk_text_buffer_create_tag(
			search_data->buffer, "search_match",
			"background", "yellow", "foreground", "black", NULL);

	tab->search_data = search_data;
}

void setup_search(GtkSourceView *source_view, GtkEntry *search_entry) {
	// This function is now handled by setup_search_for_tab for each individual tab
	// Keep for compatibility but don't use
}

// Statics init
void python_scratchpad_init_statics() {
	if (editor_window == NULL) {
		// GtkWindow
		editor_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "python_window"));
		main_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
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
		useargs = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "editor_useargs"));
		// GtkComboBox
		combo_language = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "python_pad_language"));
		// GtkButton
		button_python_pad_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_close"));
		button_python_pad_clear = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_clear"));
		button_python_pad_open = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_open"));
		button_python_pad_save = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_save"));
		button_python_pad_execute = GTK_BUTTON(gtk_builder_get_object(gui.builder, "button_python_pad_execute"));
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
		args_box = GTK_BOX(gtk_builder_get_object(gui.builder, "editor_args_box"));
		// GtkEntry
		args_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "editor_args_entry"));

		// GtkSourceView setup - now creates tabbed interface
		add_code_view(gui.builder);
		gtk_window_set_transient_for(GTK_WINDOW(editor_window), GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));

		/* initialize find box */
		setup_find_overlay();
	}
}

static void update_title(GFile *file) {
	if (!current_tab) return;

	if (file) {
		char *basename = g_file_get_basename(file);
		g_free(current_tab->title);
		current_tab->title = g_strdup(basename);

		char *suffix = strrchr(basename, '.');
		if (suffix != NULL) {
			if (g_ascii_strcasecmp(suffix, ".py") == 0) {
				current_tab->language = LANG_PYTHON;
				gtk_check_menu_item_set_active(radio_py, TRUE);
			} else if (g_ascii_strcasecmp(suffix, ".ssf") == 0) {
				current_tab->language = LANG_SSF;
				gtk_check_menu_item_set_active(radio_ssf, TRUE);
			}
		}

		g_free(basename);
	} else {
		g_free(current_tab->title);
		current_tab->title = g_strdup("untitled");
	}

	current_tab->modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(current_tab->source_buffer), FALSE);
	update_tab_title(current_tab);
	set_language();
}

// File handling
void load_file_complete(GObject *loader, GAsyncResult *result, gpointer user_data) {
	GError *error = NULL;
	if (!gtk_source_file_loader_load_finish(GTK_SOURCE_FILE_LOADER(loader), result, &error)) {
		g_printerr("Error loading file: %s\n", error->message);
		g_error_free(error);
	}
	g_object_unref(loader);
}

void load_file_into_tab(GFile *file, TabInfo *tab) {
	if (!tab) return;

	GtkSourceFile *source_file = gtk_source_file_new();
	gtk_source_file_set_location(source_file, file);

	GtkSourceFileLoader *loader = gtk_source_file_loader_new(tab->source_buffer, source_file);

	// Start the async load operation
	gtk_source_file_loader_load_async(loader,
									G_PRIORITY_DEFAULT,
									NULL,  // No cancellable
									NULL,  // No progress callback
									NULL,  // No progress data
									NULL,  // No progress notify
									load_file_complete,
									tab); // Pass tab as user data

	// Update tab file reference
	if (tab->file) g_object_unref(tab->file);
	tab->file = g_object_ref(file);
	update_title(file);

	g_object_unref(source_file);
	// loader will be unreferenced in the callback
}

void load_file(GFile *file) {
	if (current_tab) {
		load_file_into_tab(file, current_tab);
	}
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

void save_tab_file(TabInfo *tab) {
	if (!tab || !tab->file) return;

	GtkSourceFile *source_file = gtk_source_file_new();
	gtk_source_file_set_location(source_file, tab->file);

	GtkSourceFileSaver *saver = gtk_source_file_saver_new(tab->source_buffer, source_file);

	// Start the async save operation
	gtk_source_file_saver_save_async(saver,
								G_PRIORITY_DEFAULT,
								NULL,  // No cancellable
								NULL,  // No progress callback
								NULL,  // No progress data
								NULL,  // No progress notify
								save_file_complete,
								tab); // Pass tab as user data
	tab->modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(tab->source_buffer), FALSE);
	update_tab_title(tab);
	g_object_unref(source_file);
}

void save_file(GFile *file) {
	// Legacy function - now saves current tab
	if (current_tab) {
		if (current_tab->file) g_object_unref(current_tab->file);
		current_tab->file = g_object_ref(file);
		save_tab_file(current_tab);
	}
}

void on_action_file_new(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	TabInfo *tab = create_new_tab("untitled", NULL);
	current_tab = tab;
	update_ui_for_current_tab();
}

void new_script(const gchar *contents, gint length, const char *ext) {
	on_open_pythonpad(NULL, NULL);

	TabInfo *tab = create_new_tab("untitled", contents);
	current_tab = tab;

	if (!g_ascii_strcasecmp(ext, "ssf")) {
		current_tab->language = LANG_SSF;
		gtk_check_menu_item_set_active(radio_ssf, TRUE);
	} else {
		current_tab->language = LANG_PYTHON;
		gtk_check_menu_item_set_active(radio_py, TRUE);
	}

	set_language();
	gtk_window_present_with_time(editor_window, GDK_CURRENT_TIME);
	update_ui_for_current_tab();
}

void on_action_file_close(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	if (g_list_length(tabs) <= 1) {
		// This is the last tab, just hide the window instead of closing the tab
		gtk_widget_hide(GTK_WIDGET(editor_window));
	} else {
		close_tab(current_tab);
	}
}

void on_action_file_open(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
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
		char *basename = g_file_get_basename(file);

		// Create new tab for the file, unless in the "untitled" tab without unsaved changes
		TabInfo *tab = NULL;
		if (!current_tab ||
				g_strcmp0(gtk_label_get_text(GTK_LABEL(current_tab->tab_label)), "untitled") != 0 ||
				current_tab->modified) {
			tab = create_new_tab(basename, NULL);
			current_tab = tab;
		} else {
			tab = current_tab;
		}
		control_window_switch_to_tab(OUTPUT_LOGS);
		load_file_into_tab(file, tab);

		g_free(basename);
		g_object_unref(file);

		update_ui_for_current_tab();
	}

	gtk_widget_destroy(dialog);
	gtk_window_present(editor_window);
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

	// Create new tab for recent file
	char *basename = g_file_get_basename(file);
	// Create new tab for the file, unless in the "untitled" tab without unsaved changes
	TabInfo *tab = NULL;
	if (!current_tab ||
			g_strcmp0(gtk_label_get_text(GTK_LABEL(current_tab->tab_label)), "untitled") != 0 ||
			current_tab->modified) {
		tab = create_new_tab(basename, NULL);
		current_tab = tab;
	} else {
		tab = current_tab;
	}
	load_file_into_tab(file, tab);
	update_ui_for_current_tab();

	g_free(basename);
	g_object_unref(file);
	g_free(uri);
}

void on_action_file_save_as(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

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
	if (current_tab->file) {
		gtk_file_chooser_set_file(GTK_FILE_CHOOSER(dialog), current_tab->file, NULL);
	}

	gint res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
		control_window_switch_to_tab(OUTPUT_LOGS);

		// Update current_tab with the newly selected file
		if (current_tab->file)
			g_object_unref(current_tab->file);
		current_tab->file = g_object_ref(file);
		update_title(current_tab->file);

		save_tab_file(current_tab);
		g_object_unref(file);
	}

	gtk_widget_destroy(dialog);
	gtk_window_present(editor_window);
}

void on_action_file_save(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	if (!current_tab->file) {
		// No file has been saved yet, so call Save As
		on_action_file_save_as(action, parameter, user_data);
		return;
	}

	// We have a current file, just save directly to it
	control_window_switch_to_tab(OUTPUT_LOGS);
	save_tab_file(current_tab);
}

int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data) {
	if (!editor_window) {
		python_scratchpad_init_statics();
	}

	// This is a normal window with normal decorations
	gtk_window_set_type_hint(GTK_WINDOW(editor_window), GDK_WINDOW_TYPE_HINT_NORMAL);

	// Decouple window behaviour from main window
	gtk_window_set_transient_for(editor_window, NULL);
	// (Note, if the editor window is minimized and there isn't an icon in the system tray
	// it can be restored just by clicking on the script editor menu item again - no work
	// will be lost...)

	// Hide on delete
	g_signal_connect(editor_window, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

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

	if (current_tab) {
		gtk_label_set_text(language_label, current_tab->language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));
	}

	// Show the window and bring it to front
	gtk_window_present_with_time(editor_window, GDK_CURRENT_TIME);

	if (current_tab) {
		// Set the correct check menu item active
		gtk_check_menu_item_set_active(radio_py, current_tab->language == LANG_PYTHON);
		gtk_check_menu_item_set_active(radio_ssf, current_tab->language == LANG_SSF);

		// Set the right margin position
		gtk_source_view_set_right_margin_position(current_tab->source_view, com.pref.gui.editor_cfg.rmargin_pos);

		// Focus the SourceView
		gtk_widget_grab_focus(GTK_WIDGET(current_tab->source_view));
	}

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
	if (current_tab) {
		gtk_widget_set_visible(GTK_WIDGET(current_tab->minimap), com.pref.gui.editor_cfg.minimap);
	}
	gtk_widget_set_visible(GTK_WIDGET(args_box), gtk_check_menu_item_get_active(useargs));
	return 0;
}

void on_undo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	if (gtk_source_buffer_can_undo(current_tab->source_buffer)) {
		gtk_source_buffer_undo(current_tab->source_buffer);
	}
}

void on_redo(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	if (gtk_source_buffer_can_redo(current_tab->source_buffer)) {
		gtk_source_buffer_redo(current_tab->source_buffer);
	}
}

void on_cut(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	GtkTextBuffer *buffer = GTK_TEXT_BUFFER(current_tab->source_buffer);
	GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(current_tab->source_view), GDK_SELECTION_CLIPBOARD);
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
	gtk_text_buffer_cut_clipboard(buffer, clipboard, gtk_text_view_get_editable(GTK_TEXT_VIEW(current_tab->source_view)));
}

void on_copy(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!current_tab) return;

	GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(current_tab->source_view), GDK_SELECTION_CLIPBOARD);
	gtk_text_buffer_copy_clipboard(GTK_TEXT_BUFFER(current_tab->source_buffer), clipboard);
}

void on_paste(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_revealer_get_reveal_child(GTK_REVEALER(find_revealer))) {
		gtk_editable_paste_clipboard(GTK_EDITABLE(find_entry));
	} else if (current_tab) {
		GtkClipboard *clipboard = gtk_widget_get_clipboard(GTK_WIDGET(current_tab->source_view), GDK_SELECTION_CLIPBOARD);
		gtk_text_buffer_paste_clipboard(GTK_TEXT_BUFFER(current_tab->source_buffer), clipboard, NULL, gtk_text_view_get_editable(GTK_TEXT_VIEW(current_tab->source_view)));
	}
}

void on_find(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
//	update_current_tab(); // TODO: Shouldn't be needed if the switch-page handler works properly
	if (!current_tab) return;

	// Update search to work with current tab - disconnect any existing signals first
	g_signal_handlers_disconnect_by_func(find_entry, on_find_entry_changed, NULL);
	g_signal_handlers_disconnect_by_func(find_entry, on_find_entry_key_press_event, NULL);
	g_signal_handlers_disconnect_by_func(go_down_button, on_next_match, NULL);
	g_signal_handlers_disconnect_by_func(go_up_button, on_previous_match, NULL);

	SearchData *search_data = current_tab->search_data;
	if (search_data) {
		search_data->search_entry = find_entry;
		search_data->info_label = find_label;

		// Connect signals for current tab
		g_signal_connect(find_entry, "changed", G_CALLBACK(on_find_entry_changed), search_data);
		g_signal_connect(find_entry, "key-press-event", G_CALLBACK(on_find_entry_key_press_event), search_data);
		g_signal_connect(go_down_button, "clicked", G_CALLBACK(on_next_match), search_data);
		g_signal_connect(go_up_button, "clicked", G_CALLBACK(on_previous_match), search_data);
	}

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
	gint current_pos = current_tab ? gtk_source_view_get_right_margin_position(current_tab->source_view) : com.pref.gui.editor_cfg.rmargin_pos;
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
		// Apply to all tabs
		for (GList *l = tabs; l; l = l->next) {
			TabInfo *tab = l->data;
			gtk_source_view_set_right_margin_position(tab->source_view, com.pref.gui.editor_cfg.rmargin_pos);
		}
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
	}

	gtk_widget_destroy(dialog);
}

gboolean on_main_window_state_changed(GtkWidget *widget, GdkEventWindowState *event, gpointer user_data) {
	if (event->changed_mask & GDK_WINDOW_STATE_ICONIFIED) {
		GtkWindow *editor_window = GTK_WINDOW(user_data);
		if (gtk_widget_get_visible(GTK_WIDGET(editor_window))) {
			if (event->new_window_state & GDK_WINDOW_STATE_ICONIFIED) {
				gtk_window_iconify(editor_window);
			} else {
				gtk_window_present(editor_window);
			}
		}
	}
	return FALSE;
}

void set_language() {
	if (!current_tab) return;

	if (current_tab->language == LANG_PYTHON) {
		language = gtk_source_language_manager_get_language(language_manager, "python");
		gtk_label_set_text(language_label, _("Python Script"));
	} else if (current_tab->language == LANG_SSF) {
		language = gtk_source_language_manager_get_language(language_manager, "sh");
		gtk_label_set_text(language_label, _("Siril Script File"));
	}
	if (language == NULL) {
		siril_debug_print("Could not find language definition\n");
	} else {
		gtk_source_buffer_set_language(current_tab->source_buffer, language);
	}
}

void on_language_set(GtkRadioMenuItem *item, gpointer user_data) {
	if (!current_tab) return;

	if (gtk_check_menu_item_get_active(radio_py)) {
		current_tab->language = LANG_PYTHON;
	} else {
		current_tab->language = LANG_SSF;
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
	if (!current_tab) return;

	if (!accept_script_warning_dialog())
		return;
	gchar** script_args = NULL;
	// Get the start and end iterators
	GtkTextIter start, end;
	gtk_text_buffer_get_bounds(GTK_TEXT_BUFFER(current_tab->source_buffer), &start, &end);
	// Get the text
	char *text = gtk_text_buffer_get_text(GTK_TEXT_BUFFER(current_tab->source_buffer), &start, &end, FALSE);
	switch (current_tab->language) {
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
			if (gtk_check_menu_item_get_active(useargs)) {
				const gchar *args_string = gtk_entry_get_text(args_entry);
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
			g_free(text);
			break;
	}
}

static void on_tab_buffer_modified_changed(GtkTextBuffer *buffer, TabInfo *tab) {
	tab->modified = gtk_text_buffer_get_modified(buffer);
	update_tab_title(tab);
}

void on_action_python_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("Python-API.html");
}

void on_action_command_doc(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation("Commands.html");
}

void on_editor_syntax_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.highlight_syntax = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_buffer_set_highlight_syntax(tab->source_buffer, status);
	}
}

void on_editor_bracketmatch_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.highlight_bracketmatch = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_buffer_set_highlight_matching_brackets(tab->source_buffer, status);
	}
}

void on_editor_rmargin_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.rmargin = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_show_right_margin(tab->source_view, status);
	}
}

void on_editor_linenums_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.show_linenums = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_show_line_numbers(tab->source_view, status);
	}
}

void on_editor_linemarks_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.show_linemarks = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_show_line_marks(tab->source_view, status);
	}
}

void on_editor_highlightcurrentline_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.highlight_currentline = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_highlight_current_line(tab->source_view, status);
	}
}

void on_editor_autoindent_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.autoindent = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_auto_indent(tab->source_view, status);
	}
}

void on_editor_indentontab_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.indentontab = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_indent_on_tab(tab->source_view, status);
	}
}

void on_editor_smartbs_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.smartbs = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_smart_backspace(tab->source_view, status);
	}
}

void on_editor_smarthomeend_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.smarthomeend = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_source_view_set_smart_home_end(tab->source_view, status);
	}
}

void on_editor_showspaces_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.showspaces = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;

		// Define where spaces should be displayed (e.g., everywhere)
		GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

		// Get the current space types for the specified locations
		GtkSourceSpaceTypeFlags space_types = gtk_source_space_drawer_get_types_for_locations(tab->space_drawer, space_locations);

		// Define the space types to toggle (spaces and tabs)
		GtkSourceSpaceTypeFlags spaces_and_tabs = GTK_SOURCE_SPACE_TYPE_SPACE | GTK_SOURCE_SPACE_TYPE_TAB;

		// Update the space types based on the status
		if (status) {
			space_types |= spaces_and_tabs; // Enable spaces and tabs
		} else {
			space_types &= ~spaces_and_tabs; // Disable spaces and tabs
		}

		// Apply the updated space types to the drawer
		gtk_source_space_drawer_set_types_for_locations(tab->space_drawer, space_locations, space_types);
	}

	// Redraw the editor window to reflect the changes
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
}

void on_editor_shownewlines_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.shownewlines = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;

		// Define where spaces should be displayed (e.g., everywhere)
		GtkSourceSpaceLocationFlags space_locations = GTK_SOURCE_SPACE_LOCATION_ALL;

		// Get the current space types for the specified locations
		GtkSourceSpaceTypeFlags space_types = gtk_source_space_drawer_get_types_for_locations(tab->space_drawer, space_locations);

		// Update the space types based on the status
		if (status) {
			space_types |= GTK_SOURCE_SPACE_TYPE_NEWLINE; // Enable newlines
		} else {
			space_types &= ~GTK_SOURCE_SPACE_TYPE_NEWLINE; // Disable newlines
		}

		// Apply the updated space types to the drawer
		gtk_source_space_drawer_set_types_for_locations(tab->space_drawer, space_locations, space_types);
	}

	// Redraw the editor window to reflect the changes
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
}

void on_editor_minimap_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gboolean status = gtk_check_menu_item_get_active(item);
	com.pref.gui.editor_cfg.minimap = status;

	// Apply to all tabs
	for (GList *l = tabs; l; l = l->next) {
		TabInfo *tab = l->data;
		gtk_widget_set_visible(GTK_WIDGET(tab->minimap), status);
	}
}

void on_editor_useargs_toggled(GtkCheckMenuItem *item, gpointer user_data) {
	gtk_widget_set_visible(GTK_WIDGET(args_box), gtk_check_menu_item_get_active(item));
	if (gtk_check_menu_item_get_active(item)) {
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

void on_pythondebug_toggled(GtkCheckMenuItem *item, gpointer user_data) {
    gboolean state = gtk_check_menu_item_get_active(item);
	GtkCheckMenuItem *editorwidget = (GTK_CHECK_MENU_ITEM(lookup_widget("editor_toggledebug")));
	// This is created programatically and added to the builder, but we know it has been
	// done by now as the script menu must exist in order to be able to access either the script
	// menu item or the script editor menu item.
	GtkCheckMenuItem *scriptmenuwidget = (GTK_CHECK_MENU_ITEM(lookup_widget("pythondebugtoggle")));

	// Synchronize the two checkbuttons
	g_signal_handlers_block_by_func(editorwidget, on_pythondebug_toggled, NULL);
	g_signal_handlers_block_by_func(scriptmenuwidget, on_pythondebug_toggled, NULL);
	gtk_check_menu_item_set_active(editorwidget, state);
	gtk_check_menu_item_set_active(scriptmenuwidget, state);
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
