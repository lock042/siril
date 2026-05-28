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
static GtkWidget    *map_wrapper = NULL;  /* clip container holding `map` */
static GtkScrolledWindow *scrolled_window = NULL;
static GtkDropDown *combo_language = NULL;
static GtkSourceLanguageManager *language_manager = NULL;
static GtkSourceLanguage *language = NULL;
static GtkSourceStyleSchemeManager *stylemanager = NULL;
static GtkSourceStyleScheme *scheme = NULL;
static GFile *current_file = NULL;
static GtkWindow *editor_window = NULL;
static GtkWindow *main_window = NULL;
static GtkPopoverMenuBar *pyeditor_menubar = NULL;
/* Recent Files: rebuilt each time the editor opens; owned by the menubar
 * model and referenced from inside the File submenu. */
static GMenu *pyeditor_recent_menu_model = NULL;
/* Editor action group; GTK4 has no gtk_widget_get_action_group(), so we
 * retain a reference to set states from outside on_open_pythonpad. */
static GSimpleActionGroup *editor_action_group = NULL;
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

static void on_action_open_recent(GSimpleAction *action, GVariant *parameter, gpointer user_data);
static void on_toggle_change_state(GSimpleAction *action, GVariant *value, gpointer user_data);
static void on_language_change_state(GSimpleAction *action, GVariant *value, gpointer user_data);

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
	{ "reload", on_action_file_reload, NULL, NULL, NULL },
	{ "open_recent", on_action_open_recent, "s", NULL, NULL },

	/* Boolean stateful actions: rendered as check items in the menu.
	 * Activation with no parameter flips the state via GIO's default,
	 * which delegates to change_state — our handler updates the source
	 * view, the saved pref, and writes the new state back. */
	{ "syntax",               NULL, NULL, "true",  on_toggle_change_state },
	{ "rmargin",              NULL, NULL, "true",  on_toggle_change_state },
	{ "bracketmatch",         NULL, NULL, "true",  on_toggle_change_state },
	{ "linenums",             NULL, NULL, "true",  on_toggle_change_state },
	{ "linemarks",            NULL, NULL, "true",  on_toggle_change_state },
	{ "highlightcurrentline", NULL, NULL, "true",  on_toggle_change_state },
	{ "autoindent",           NULL, NULL, "true",  on_toggle_change_state },
	{ "indentontab",          NULL, NULL, "true",  on_toggle_change_state },
	{ "smartbs",              NULL, NULL, "true",  on_toggle_change_state },
	{ "smarthomeend",         NULL, NULL, "true",  on_toggle_change_state },
	{ "showspaces",           NULL, NULL, "false", on_toggle_change_state },
	{ "shownewlines",         NULL, NULL, "false", on_toggle_change_state },
	{ "minimap",              NULL, NULL, "false", on_toggle_change_state },
	{ "useargs",              NULL, NULL, "false", on_toggle_change_state },
	{ "pythondebug",          NULL, NULL, "false", on_toggle_change_state },

	/* String stateful action: rendered as a radio group when multiple
	 * menu items share the action name with different `target` values. */
	{ "language",             NULL, "s",  "'py'",  on_language_change_state },
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

/* Append a labelled action item to `menu` with an optional GTK accelerator
 * string (e.g. "<Control>c").  When `accel` is non-NULL, the popover
 * renders it on the right of the item, matching the GTK3 menubar look.
 * The accelerator only triggers the action because we install matching
 * shortcuts on the editor window via install_editor_shortcuts(). */
static void pyeditor_menu_append(GMenu *menu, const gchar *label,
                                 const gchar *action, const gchar *accel) {
	GMenuItem *item = g_menu_item_new(label, action);
	if (accel)
		g_menu_item_set_attribute_value(item, "accel", g_variant_new_string(accel));
	g_menu_append_item(menu, item);
	g_object_unref(item);
}

/* Build the GMenuModel that backs the script editor's PopoverMenuBar.
 * The Recent Files submenu links to `pyeditor_recent_menu_model`, which
 * is repopulated each time the editor is opened.  Labels carry GTK-style
 * mnemonics (underscore-prefixed letter) so Alt+<key> navigation matches
 * the GTK3 menu bar. */
static GMenuModel *build_pyeditor_menubar_model(void) {
	if (!pyeditor_recent_menu_model)
		pyeditor_recent_menu_model = g_menu_new();

	GMenu *bar = g_menu_new();

	/* File */
	GMenu *file = g_menu_new();
	GMenu *file_top = g_menu_new();
	pyeditor_menu_append(file_top, _("_New"),  "editor.new",  "<Control>n");
	pyeditor_menu_append(file_top, _("_Open"), "editor.open", "<Control>o");
	GMenuItem *recent_item = g_menu_item_new(_("_Recent Files"), NULL);
	g_menu_item_set_submenu(recent_item, G_MENU_MODEL(pyeditor_recent_menu_model));
	g_menu_append_item(file_top, recent_item);
	g_object_unref(recent_item);
	g_menu_append_section(file, NULL, G_MENU_MODEL(file_top));
	g_object_unref(file_top);
	GMenu *file_save = g_menu_new();
	pyeditor_menu_append(file_save, _("_Save"),    "editor.save",    "<Control>s");
	pyeditor_menu_append(file_save, _("Save _As"), "editor.save_as", "<Control><Shift>s");
	g_menu_append_section(file, NULL, G_MENU_MODEL(file_save));
	g_object_unref(file_save);
	GMenu *file_close = g_menu_new();
	pyeditor_menu_append(file_close, _("_Close"), "editor.close", "<Control>q");
	g_menu_append_section(file, NULL, G_MENU_MODEL(file_close));
	g_object_unref(file_close);
	g_menu_append_submenu(bar, _("_File"), G_MENU_MODEL(file));
	g_object_unref(file);

	/* Edit */
	GMenu *edit = g_menu_new();
	GMenu *edit_undo = g_menu_new();
	pyeditor_menu_append(edit_undo, _("_Undo"), "editor.undo", "<Control>z");
	pyeditor_menu_append(edit_undo, _("_Redo"), "editor.redo", "<Control>y");
	g_menu_append_section(edit, NULL, G_MENU_MODEL(edit_undo));
	g_object_unref(edit_undo);
	GMenu *edit_cp = g_menu_new();
	pyeditor_menu_append(edit_cp, _("Cu_t"),   "editor.cut",   "<Control>x");
	pyeditor_menu_append(edit_cp, _("_Copy"),  "editor.copy",  "<Control>c");
	pyeditor_menu_append(edit_cp, _("_Paste"), "editor.paste", "<Control>v");
	g_menu_append_section(edit, NULL, G_MENU_MODEL(edit_cp));
	g_object_unref(edit_cp);
	GMenu *edit_find = g_menu_new();
	pyeditor_menu_append(edit_find, _("_Find"), "editor.find", "<Control>f");
	g_menu_append_section(edit, NULL, G_MENU_MODEL(edit_find));
	g_object_unref(edit_find);
	g_menu_append_submenu(bar, _("_Edit"), G_MENU_MODEL(edit));
	g_object_unref(edit);

	/* Script */
	GMenu *script = g_menu_new();
	GMenu *script_run = g_menu_new();
	pyeditor_menu_append(script_run, _("_Run"), "editor.execute", "F5");
	g_menu_append_section(script, NULL, G_MENU_MODEL(script_run));
	g_object_unref(script_run);
	GMenu *script_lang = g_menu_new();
	GMenuItem *lang_py = g_menu_item_new(_("_Python scripts"), NULL);
	g_menu_item_set_action_and_target_value(lang_py, "editor.language",
		g_variant_new_string("py"));
	g_menu_append_item(script_lang, lang_py);
	g_object_unref(lang_py);
	GMenuItem *lang_ssf = g_menu_item_new(_("_Siril Script Files"), NULL);
	g_menu_item_set_action_and_target_value(lang_ssf, "editor.language",
		g_variant_new_string("ssf"));
	g_menu_append_item(script_lang, lang_ssf);
	g_object_unref(lang_ssf);
	g_menu_append_section(script, NULL, G_MENU_MODEL(script_lang));
	g_object_unref(script_lang);
	GMenu *script_opts = g_menu_new();
	g_menu_append(script_opts, _("Enable test _arguments"),  "editor.useargs");
	g_menu_append(script_opts, _("Enable python _debug mode"), "editor.pythondebug");
	g_menu_append_section(script, NULL, G_MENU_MODEL(script_opts));
	g_object_unref(script_opts);
	g_menu_append_submenu(bar, _("_Script"), G_MENU_MODEL(script));
	g_object_unref(script);

	/* Preferences */
	GMenu *prefs = g_menu_new();
	g_menu_append(prefs, _("_Highlight syntax"),             "editor.syntax");
	g_menu_append(prefs, _("Enable right _margin indicator"), "editor.rmargin");
	g_menu_append(prefs, _("Set right margin _position"),    "editor.set_rmarginpos");
	g_menu_append(prefs, _("Enable _bracket matching"),      "editor.bracketmatch");
	g_menu_append(prefs, _("Show _line numbers"),            "editor.linenums");
	g_menu_append(prefs, _("Show line ma_rks"),              "editor.linemarks");
	g_menu_append(prefs, _("Highlight c_urrent line"),       "editor.highlightcurrentline");
	g_menu_append(prefs, _("Enable _auto-indentation"),      "editor.autoindent");
	g_menu_append(prefs, _("_Indent on tab"),                "editor.indentontab");
	g_menu_append(prefs, _("Enable smart _backspace"),       "editor.smartbs");
	g_menu_append(prefs, _("Smart H_ome / End"),             "editor.smarthomeend");
	g_menu_append(prefs, _("Show spaces and _tabs"),         "editor.showspaces");
	g_menu_append(prefs, _("Show _newlines"),                "editor.shownewlines");
	g_menu_append(prefs, _("Show mi_nimap"),                 "editor.minimap");
	g_menu_append_submenu(bar, _("_Preferences"), G_MENU_MODEL(prefs));
	g_object_unref(prefs);

	/* Help */
	GMenu *help = g_menu_new();
	g_menu_append(help, _("Python _API reference"),     "editor.python_doc");
	g_menu_append(help, _("Siril _Command Reference"),  "editor.command_doc");
	g_menu_append_submenu(bar, _("_Help"), G_MENU_MODEL(help));
	g_object_unref(help);

	return G_MENU_MODEL(bar);
}

/* Set a stateful editor action's state directly, bypassing the
 * change_state handler.  Used during initialization so the menu's check
 * marks reflect the saved prefs without re-running the side-effects that
 * we apply explicitly in apply_editor_setting (avoids double-application
 * and any timing edge case where the handler doesn't fire). */
static void editor_action_set_state(const gchar *name, GVariant *value) {
	if (!editor_action_group) {
		g_variant_unref(g_variant_ref_sink(value));
		return;
	}
	GAction *act = g_action_map_lookup_action(G_ACTION_MAP(editor_action_group), name);
	if (act)
		g_simple_action_set_state(G_SIMPLE_ACTION(act), value);
	else
		g_variant_unref(g_variant_ref_sink(value));
}

/* Request a state change on an editor action (fires the change_state
 * handler).  Used when external code wants the handler's side-effect to
 * run — e.g. switching the language radio when loading a .py / .ssf
 * file needs set_language() to rebind syntax highlighting. */
static void editor_action_change_state(const gchar *name, GVariant *value) {
	if (!editor_action_group) {
		g_variant_unref(g_variant_ref_sink(value));
		return;
	}
	g_action_group_change_action_state(G_ACTION_GROUP(editor_action_group), name, value);
}

/* Apply a boolean editor setting to the source view / wrapper widgets and
 * persist it in com.pref.gui.editor_cfg.  Centralises the side-effect of
 * each toggle so both the change_state handler (when the user clicks a
 * menu item) and the init path (when the editor opens) reach the same
 * code with no duplication. */
static void apply_editor_setting(const gchar *name, gboolean v) {
	if (g_strcmp0(name, "syntax") == 0) {
		gtk_source_buffer_set_highlight_syntax(sourcebuffer, v);
		com.pref.gui.editor_cfg.highlight_syntax = v;
	} else if (g_strcmp0(name, "bracketmatch") == 0) {
		gtk_source_buffer_set_highlight_matching_brackets(sourcebuffer, v);
		com.pref.gui.editor_cfg.highlight_bracketmatch = v;
	} else if (g_strcmp0(name, "rmargin") == 0) {
		gtk_source_view_set_show_right_margin(code_view, v);
		com.pref.gui.editor_cfg.rmargin = v;
		gtk_widget_queue_draw(GTK_WIDGET(code_view));
	} else if (g_strcmp0(name, "linenums") == 0) {
		gtk_source_view_set_show_line_numbers(code_view, v);
		com.pref.gui.editor_cfg.show_linenums = v;
	} else if (g_strcmp0(name, "linemarks") == 0) {
		gtk_source_view_set_show_line_marks(code_view, v);
		com.pref.gui.editor_cfg.show_linemarks = v;
	} else if (g_strcmp0(name, "highlightcurrentline") == 0) {
		gtk_source_view_set_highlight_current_line(code_view, v);
		com.pref.gui.editor_cfg.highlight_currentline = v;
	} else if (g_strcmp0(name, "autoindent") == 0) {
		gtk_source_view_set_auto_indent(code_view, v);
		com.pref.gui.editor_cfg.autoindent = v;
	} else if (g_strcmp0(name, "indentontab") == 0) {
		gtk_source_view_set_indent_on_tab(code_view, v);
		com.pref.gui.editor_cfg.indentontab = v;
	} else if (g_strcmp0(name, "smartbs") == 0) {
		gtk_source_view_set_smart_backspace(code_view, v);
		com.pref.gui.editor_cfg.smartbs = v;
	} else if (g_strcmp0(name, "smarthomeend") == 0) {
		gtk_source_view_set_smart_home_end(code_view, v ?
			GTK_SOURCE_SMART_HOME_END_AFTER : GTK_SOURCE_SMART_HOME_END_DISABLED);
		com.pref.gui.editor_cfg.smarthomeend = v;
	} else if (g_strcmp0(name, "showspaces") == 0 ||
	           g_strcmp0(name, "shownewlines") == 0) {
		GtkSourceSpaceLocationFlags loc = GTK_SOURCE_SPACE_LOCATION_ALL;
		GtkSourceSpaceTypeFlags types = gtk_source_space_drawer_get_types_for_locations(space_drawer, loc);
		GtkSourceSpaceTypeFlags mask = g_strcmp0(name, "showspaces") == 0
			? (GTK_SOURCE_SPACE_TYPE_SPACE | GTK_SOURCE_SPACE_TYPE_TAB)
			: GTK_SOURCE_SPACE_TYPE_NEWLINE;
		types = v ? (types | mask) : (types & ~mask);
		gtk_source_space_drawer_set_types_for_locations(space_drawer, loc, types);
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
		if (g_strcmp0(name, "showspaces") == 0)
			com.pref.gui.editor_cfg.showspaces = v;
		else
			com.pref.gui.editor_cfg.shownewlines = v;
	} else if (g_strcmp0(name, "minimap") == 0) {
		if (map_wrapper)
			gtk_widget_set_visible(map_wrapper, v);
		com.pref.gui.editor_cfg.minimap = v;
	} else if (g_strcmp0(name, "useargs") == 0) {
		if (args_box)
			gtk_widget_set_visible(GTK_WIDGET(args_box), v);
		from_cli = v;
	} else if (g_strcmp0(name, "pythondebug") == 0) {
		python_debug = v;
	}
}

/* Install keyboard shortcuts that target the editor window's GAction
 * group.  Run once per window: a GTK_PHASE_CAPTURE controller would
 * intercept GtkSourceView's built-in editing keys, so we attach to the
 * window (default bubble phase) — that way the source view still handles
 * Ctrl+C / Ctrl+X / Ctrl+V / Ctrl+Z internally, while shortcuts like
 * Ctrl+S, Ctrl+O, F5 reach the editor's GAction handler. */
static void install_editor_shortcuts(GtkWindow *window) {
	static const struct { const gchar *trigger; const gchar *action; } accels[] = {
		{ "<Control>n",          "editor.new" },
		{ "<Control>o",          "editor.open" },
		{ "<Control>s",          "editor.save" },
		{ "<Control><Shift>s",   "editor.save_as" },
		{ "<Control>q",          "editor.close" },
		{ "<Control>z",          "editor.undo" },
		{ "<Control>y",          "editor.redo" },
		{ "<Control>x",          "editor.cut" },
		{ "<Control>c",          "editor.copy" },
		{ "<Control>v",          "editor.paste" },
		{ "<Control>f",          "editor.find" },
		{ "F5",                  "editor.execute" },
		{ "<Control>r",          "editor.execute" },
	};
	/* MANAGED scope: the shortcut fires when any widget in the editor
	 * window's focus chain is focused (typically the GtkSourceView),
	 * not just when the window itself has focus.  Default LOCAL scope
	 * would never trigger since focus is on a child widget. */
	GtkShortcutController *ctrl = GTK_SHORTCUT_CONTROLLER(gtk_shortcut_controller_new());
	gtk_shortcut_controller_set_scope(ctrl, GTK_SHORTCUT_SCOPE_MANAGED);
	for (size_t i = 0; i < G_N_ELEMENTS(accels); i++) {
		GtkShortcutTrigger *trigger = gtk_shortcut_trigger_parse_string(accels[i].trigger);
		GtkShortcutAction *act = gtk_named_action_new(accels[i].action);
		GtkShortcut *sc = gtk_shortcut_new(trigger, act);
		gtk_shortcut_controller_add_shortcut(ctrl, sc);
	}
	gtk_widget_add_controller(GTK_WIDGET(window), GTK_EVENT_CONTROLLER(ctrl));
}

// Add this function to initialize the actions
static void setup_editor_actions(GtkWindow *window) {
	/* The action group is built once and reused: rebuilding it would
	 * reset every check item to its declared initial state on each
	 * re-open, undoing in-session edits. */
	if (!editor_action_group) {
		editor_action_group = g_simple_action_group_new();
		g_action_map_add_action_entries(G_ACTION_MAP(editor_action_group),
		                                editor_actions,
		                                G_N_ELEMENTS(editor_actions),
		                                window);
	}

	gtk_widget_insert_action_group(GTK_WIDGET(window), "editor",
	                               G_ACTION_GROUP(editor_action_group));

	/* Attach the menu model and shortcut controller on first setup;
	 * both reference actions in the group we just installed.  The
	 * menu model itself is static — only the Recent Files submenu's
	 * contents change. */
	if (pyeditor_menubar && !gtk_popover_menu_bar_get_menu_model(pyeditor_menubar)) {
		GMenuModel *model = build_pyeditor_menubar_model();
		gtk_popover_menu_bar_set_menu_model(pyeditor_menubar, model);
		g_object_unref(model);
		install_editor_shortcuts(window);
	}
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
	 * a thumbnail rather than a second editor.  GtkSourceMap IS itself
	 * the textview CSS node (it inherits from GtkSourceView → GtkTextView)
	 * — the previous `.pyeditor-map textview` selector searched for a
	 * child node that doesn't exist, so the rule never applied and the
	 * map kept default-font sizing.  Match the map widget directly and
	 * scope to `text`/`textview` parts to leave the slider node alone. */
	siril_register_css_for_display("siril-pyeditor-map",
		".pyeditor-map, .pyeditor-map text { font-size: 1px; }");
	gtk_widget_add_css_class(GTK_WIDGET(map), "pyeditor-map");

	gtk_widget_set_visible(GTK_WIDGET(map), TRUE);
	gtk_source_map_set_view(map, code_view);
	gtk_widget_set_valign(GTK_WIDGET(map), GTK_ALIGN_FILL);
	gtk_widget_set_hexpand(GTK_WIDGET(map), FALSE);

	/* GTK4 has no max-width on widgets — `gtk_widget_set_size_request`
	 * is a MINIMUM, not a maximum, and the natural-width measurement of
	 * GtkSourceMap (font-driven) propagates upward through the box and
	 * makes the map grow huge if the font CSS hasn't applied yet.  Wrap
	 * the map in a GtkScrolledWindow with `propagate-natural-width =
	 * FALSE` so the scroller's own size_request decides the allocation;
	 * the inner map is then clipped (no scrollbars enabled). */
	map_wrapper = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(map_wrapper),
	                               GTK_POLICY_NEVER, GTK_POLICY_NEVER);
	gtk_scrolled_window_set_propagate_natural_width(
		GTK_SCROLLED_WINDOW(map_wrapper), FALSE);
	gtk_scrolled_window_set_propagate_natural_height(
		GTK_SCROLLED_WINDOW(map_wrapper), TRUE);
	gtk_widget_set_size_request(map_wrapper, 120, -1);
	gtk_widget_set_hexpand(map_wrapper, FALSE);
	gtk_widget_set_valign(map_wrapper, GTK_ALIGN_FILL);
	gtk_widget_set_halign(map_wrapper, GTK_ALIGN_END);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(map_wrapper), GTK_WIDGET(map));
	gtk_box_append(GTK_BOX(codeviewbox), map_wrapper);

	// Set the GtkSourceView style depending on whether the Siril light or dark
	// theme is set.
	set_code_view_theme();

	// Get the GtkSourceLanguageManager and set the Python language
	language_manager = gtk_source_language_manager_get_default();
	language = gtk_source_language_manager_get_language(language_manager, "python3");
	if (language == NULL) {
		siril_log_debug("Could not find Python language definition\n");
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
		// GtkPopoverMenuBar (toggle/radio state lives in the editor action
		// group; menu items are built from a GMenuModel in setup_editor_actions).
		pyeditor_menubar = GTK_POPOVER_MENU_BAR(gtk_builder_get_object(gui.builder, "pyeditor_menubar"));
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
				editor_action_change_state("language", g_variant_new_string("py"));
			} else if (g_ascii_strcasecmp(suffix, ".ssf") == 0) {
				editor_action_change_state("language", g_variant_new_string("ssf"));
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
			editor_action_change_state("language", g_variant_new_string("ssf"));
		} else {
			active_language = LANG_PYTHON;
			editor_action_change_state("language", g_variant_new_string("py"));
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

/* "editor.open_recent" action — takes a string GVariant holding the
 * recent file path.  Loads the file into the editor buffer. */
static void on_action_open_recent(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	(void)action; (void)user_data;
	if (!parameter) return;
	const gchar *path = g_variant_get_string(parameter, NULL);
	if (!path || !*path) return;
	if (!g_file_test(path, G_FILE_TEST_EXISTS)) {
		siril_log_warning(_("Recent script no longer exists: %s\n"), path);
		gchar *uri = g_filename_to_uri(path, NULL, NULL);
		if (uri) {
			gtk_recent_manager_remove_item(
				gtk_recent_manager_get_default(), uri, NULL);
			g_free(uri);
		}
		return;
	}
	if (buffer_modified &&
	    !siril_confirm_dialog(_("Are you sure?"),
		_("This will replace the entry buffer. You will not be able to recover any contents."),
		_("Proceed")))
		return;
	GFile *file = g_file_new_for_path(path);
	control_window_switch_to_tab(OUTPUT_LOGS);
	load_file(file);
	if (G_IS_OBJECT(current_file))
		g_object_unref(current_file);
	current_file = g_object_ref(file);
	update_title(current_file);
	g_object_unref(file);
}

static int pyeditor_recent_cmp_visited_desc(gconstpointer a, gconstpointer b) {
	GDateTime *ta = gtk_recent_info_get_visited((GtkRecentInfo *) a);
	GDateTime *tb = gtk_recent_info_get_visited((GtkRecentInfo *) b);
	if (!ta && !tb) return 0;
	if (!ta) return 1;
	if (!tb) return -1;
	return g_date_time_compare(tb, ta);
}

/* Build (or rebuild) the Recent submenu inside the script editor's File
 * menu.  Filters to script extensions (.py, .ssf).  Called each time the
 * editor opens so it picks up new recent entries.  The menubar
 * references this GMenu by pointer; clearing/appending here updates the
 * popover contents automatically. */
static void populate_pyeditor_recent_menu(void) {
	if (!pyeditor_recent_menu_model)
		pyeditor_recent_menu_model = g_menu_new();
	g_menu_remove_all(pyeditor_recent_menu_model);

	GtkRecentManager *mgr = gtk_recent_manager_get_default();
	GList *items = mgr ? g_list_sort(gtk_recent_manager_get_items(mgr),
	                                 pyeditor_recent_cmp_visited_desc)
	                   : NULL;
	int count = 0;
	const int max_items = 15;
	for (GList *l = items; l && count < max_items; l = l->next) {
		GtkRecentInfo *info = l->data;
		const gchar *uri = gtk_recent_info_get_uri(info);
		if (!uri) continue;
		gchar *path = g_filename_from_uri(uri, NULL, NULL);
		if (!path || !g_file_test(path, G_FILE_TEST_EXISTS)) {
			g_free(path);
			continue;
		}
		gchar *lower = g_ascii_strdown(path, -1);
		gboolean is_script = g_str_has_suffix(lower, ".py") ||
		                     g_str_has_suffix(lower, ".ssf");
		g_free(lower);
		if (!is_script) {
			g_free(path);
			continue;
		}
		gchar *basename = g_path_get_basename(path);
		GMenuItem *mi = g_menu_item_new(basename, NULL);
		g_menu_item_set_action_and_target_value(mi, "editor.open_recent",
			g_variant_new_string(path));
		g_menu_append_item(pyeditor_recent_menu_model, mi);
		g_object_unref(mi);
		g_free(basename);
		g_free(path);
		count++;
	}
	if (items) g_list_free_full(items, (GDestroyNotify) gtk_recent_info_unref);

	if (count == 0) {
		/* Empty placeholder so the submenu still opens but conveys
		 * there's nothing to pick.  Disabled by targeting a no-op
		 * action name that isn't registered. */
		GMenuItem *empty = g_menu_item_new(_("(no recent scripts)"), "editor.noop_disabled");
		g_menu_append_item(pyeditor_recent_menu_model, empty);
		g_object_unref(empty);
	}
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
	/* The menu bar references this submenu by pointer, so populate it
	 * before showing the window — first open builds the bar, subsequent
	 * opens just refresh the existing one. */
	populate_pyeditor_recent_menu();

	/* Phase 17.5: gtk_window_get_position / gtk_window_get_size /
	 * gtk_window_move are all gone in GTK4 — placement is
	 * compositor-managed.  We drop the cascading-offset logic and
	 * rely on the WM to place the editor window. */
	(void)main_window;

	gtk_label_set_text(language_label, active_language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));

	/* Right-margin position is a numeric pref (not a toggle), so it
	 * has no action; apply directly. */
	gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);

	/* Apply every editor toggle from com.pref.gui.editor_cfg directly
	 * to the source view and helper widgets BEFORE presenting the
	 * window, so the first paint is correct.  We then sync the GAction
	 * states so the menu's check marks match — without re-running the
	 * change_state handler, which would just redo the same work. */
	struct { const gchar *name; gboolean v; } setup[] = {
		{ "syntax",               com.pref.gui.editor_cfg.highlight_syntax },
		{ "bracketmatch",         com.pref.gui.editor_cfg.highlight_bracketmatch },
		{ "rmargin",              com.pref.gui.editor_cfg.rmargin },
		{ "linenums",             com.pref.gui.editor_cfg.show_linenums },
		{ "linemarks",            com.pref.gui.editor_cfg.show_linemarks },
		{ "highlightcurrentline", com.pref.gui.editor_cfg.highlight_currentline },
		{ "autoindent",           com.pref.gui.editor_cfg.autoindent },
		{ "indentontab",          com.pref.gui.editor_cfg.indentontab },
		{ "smartbs",              com.pref.gui.editor_cfg.smartbs },
		{ "smarthomeend",         com.pref.gui.editor_cfg.smarthomeend },
		{ "showspaces",           com.pref.gui.editor_cfg.showspaces },
		{ "shownewlines",         com.pref.gui.editor_cfg.shownewlines },
		{ "minimap",              com.pref.gui.editor_cfg.minimap },
		/* useargs / pythondebug keep their session values rather than
		 * being read from the pref (they're not persisted). */
		{ "useargs",              from_cli },
		{ "pythondebug",          python_debug },
	};
	for (size_t i = 0; i < G_N_ELEMENTS(setup); i++) {
		apply_editor_setting(setup[i].name, setup[i].v);
		editor_action_set_state(setup[i].name, g_variant_new_boolean(setup[i].v));
	}
	editor_action_set_state("language",
		g_variant_new_string(active_language == LANG_PYTHON ? "py" : "ssf"));

	// Show the window and bring it to front
	gtk_window_present(editor_window);

	// Focus the SourceView
	gtk_widget_grab_focus(GTK_WIDGET(code_view));

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

/* GTK4: replaces the GTK3 window-state-event handler.  Wired to the main
 * window's "notify::minimized" so that minimising the main window also
 * minimises the editor (and vice-versa on un-minimise).  GTK4 has no
 * gtk_window_is_minimized() accessor — read the property directly. */
static void on_main_window_state_changed(GObject *obj, GParamSpec *pspec,
                                         gpointer user_data) {
	(void)pspec;
	GtkWindow *editor = GTK_WINDOW(user_data);
	if (!editor || !gtk_widget_get_visible(GTK_WIDGET(editor)))
		return;
	gboolean minimized = FALSE;
	g_object_get(obj, "minimized", &minimized, NULL);
	if (minimized)
		gtk_window_minimize(editor);
	else
		gtk_window_unminimize(editor);
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
		siril_log_debug("Could not find language definition\n");
	} else {
		gtk_source_buffer_set_language(sourcebuffer, language);
	}
}

/* String-state action for the language radio.  Activating an item with a
 * different target string changes the state through here, which records
 * the language and rebinds the SourceView's syntax definition. */
static void on_language_change_state(GSimpleAction *action, GVariant *value, gpointer user_data) {
	(void)user_data;
	g_simple_action_set_state(action, value);
	const gchar *v = g_variant_get_string(value, NULL);
	if (g_strcmp0(v, "ssf") == 0) {
		active_language = LANG_SSF;
	} else {
		active_language = LANG_PYTHON;
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
					siril_log_debug("g_unlink() failed in on_action_file_execute()\n");
				g_free(temp_filename);
				g_free(text);
				return;
			}
			close(fd);
			// Add args if we need to
			if (from_cli) {
				const gchar *args_string = gtk_editable_get_text(GTK_EDITABLE(args_entry));
				script_args = g_strsplit(args_string, " ", -1);
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
			siril_log_debug("Error: unknown script language\n");
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

/* change_state dispatcher invoked when the user clicks a check item.
 * The default activate behavior for a stateful boolean action inverts
 * the state via change_state, so this is the single funnel for in-menu
 * toggles — it just records the new state and delegates the widget
 * side-effects to apply_editor_setting. */
static void on_toggle_change_state(GSimpleAction *action, GVariant *value, gpointer user_data) {
	(void)user_data;
	const gboolean v = g_variant_get_boolean(value);
	const gchar *name = g_action_get_name(G_ACTION(action));
	g_simple_action_set_state(action, value);
	apply_editor_setting(name, v);
}

void on_editor_args_clear_clicked(GtkButton *button, gpointer user_data) {
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(args_entry);
	gtk_entry_buffer_set_text(buffer, "", 0);
	gtk_widget_grab_focus(GTK_WIDGET(args_entry));
}

gboolean get_python_debug_mode() {
	return python_debug;
}
