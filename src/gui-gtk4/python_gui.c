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
#include "gui-gtk4/callbacks.h"
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
/* EXPERIMENTAL: transparent overlay that paints grey indent/tab guides in the
 * virtual hanging-indent area of wrapped continuation rows (Kate style). */
static GtkWidget    *wrap_guides_area = NULL;
/* Multi-tab: notebook holding one page per open document, and the GtkOverlay
 * (created by setup_find_overlay) that wraps the active editor scroller. */
static GtkNotebook  *editor_notebook = NULL;
static GtkWidget    *content_overlay = NULL;

/* One open document per notebook page.  The module-level statics above
 * (code_view, sourcebuffer, scrolled_window, map, map_wrapper, wrap_guides_area,
 * space_drawer, content_overlay, current_file, active_language, buffer_modified)
 * always mirror the ACTIVE document; on tab switch we save them into the
 * outgoing doc and load them from the incoming one, so the existing
 * statics-based code keeps operating on the active tab. */
typedef struct {
	GtkWidget            *page;          /* notebook page child (hbox) */
	GtkWidget            *overlay;       /* GtkOverlay wrapping the scroller */
	GtkScrolledWindow    *scroll;
	GtkSourceView        *view;
	GtkSourceBuffer      *buffer;
	GtkSourceMap         *map;
	GtkWidget            *map_wrapper;
	GtkWidget            *guides;        /* wrap guides drawing area */
	GtkSourceSpaceDrawer *space_drawer;
	GFile                *current_file;
	gint                  active_language;
	gboolean              buffer_modified;
	GtkWidget            *tab_label;     /* label widget shown in the tab */
} EditorDoc;
static GPtrArray  *editor_docs = NULL;   /* EditorDoc* */
static EditorDoc  *active_doc  = NULL;
static SearchData *active_search = NULL; /* search data of the active view */
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
static void wrap_set_enabled(gboolean on);
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
	{ "dynamicwrap",          NULL, NULL, "false", on_toggle_change_state },
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
										siril_current_theme_is_dark() ? "oblivion" : "classic");
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
	g_menu_append(prefs, _("Dynamic line _wrap"),            "editor.dynamicwrap");
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
	} else if (g_strcmp0(name, "dynamicwrap") == 0) {
		wrap_set_enabled(v);
		com.pref.gui.editor_cfg.dynamic_wrap = v;
	} else if (g_strcmp0(name, "useargs") == 0) {
		if (args_box)
			gtk_widget_set_visible(GTK_WIDGET(args_box), v);
		from_cli = v;
	} else if (g_strcmp0(name, "pythondebug") == 0) {
		set_python_debug_mode(v);
	}
}

/* Set the shared Python debug flag and keep both toggle UIs in sync: the
 * editor's Script ▸ debug check item (editor.pythondebug) and the headerbar
 * Scripts ▸ "Enable Python debug mode" check item (win.script-pythondebug).
 * g_simple_action_set_state updates a check mark WITHOUT re-invoking the
 * change-state handler, so reflecting the value onto the sibling action can't
 * loop back here.  Either action map may be absent (the editor group isn't
 * built until the editor is first opened), so each update is guarded. */
void set_python_debug_mode(gboolean enabled) {
	python_debug = enabled;

	if (editor_action_group) {
		GAction *a = g_action_map_lookup_action(G_ACTION_MAP(editor_action_group), "pythondebug");
		if (a)
			g_simple_action_set_state(G_SIMPLE_ACTION(a), g_variant_new_boolean(enabled));
	}

	GtkWidget *control_window = lookup_widget("control_window");
	if (control_window && GTK_IS_APPLICATION_WINDOW(control_window)) {
		GAction *a = g_action_map_lookup_action(G_ACTION_MAP(control_window), "script-pythondebug");
		if (a)
			g_simple_action_set_state(G_SIMPLE_ACTION(a), g_variant_new_boolean(enabled));
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

	/* The debug flag is shared with the headerbar Scripts menu's toggle and
	 * may have been flipped there before the editor was first opened, so
	 * reflect the current value onto this group's check item every time the
	 * editor is set up.  (Other toggles keep their in-session state.) */
	{
		GAction *dbg = g_action_map_lookup_action(G_ACTION_MAP(editor_action_group), "pythondebug");
		if (dbg)
			g_simple_action_set_state(G_SIMPLE_ACTION(dbg), g_variant_new_boolean(python_debug));
	}

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

/* ---------------------------------------------------------------------------
 * Changed-line marks ("change bar")
 *
 * A thin gutter bar marks every line that differs from the last-saved version
 * (amber) or that was changed earlier this session but is now saved (green),
 * plus a small caret where lines were deleted (red while unsaved, green once
 * the deletion has been saved).
 *
 * The state hangs off the GtkSourceBuffer via g_object_set_data_full (same
 * pattern as SearchData on the view) so it is naturally per-document and will
 * survive the move to one-buffer-per-tab unchanged.
 * ------------------------------------------------------------------------- */

typedef enum { CHG_NONE = 0, CHG_UNSAVED, CHG_SAVED } ChangeState;
typedef enum { DEL_NONE = 0, DEL_UNSAVED, DEL_SAVED } DelState;

typedef struct {
	GtkTextBuffer *buffer;          /* not owned */
	gchar **open_lines; guint open_n;   /* content at open / new */
	gchar **save_lines; guint save_n;   /* content at last save */
	GArray *line_state;             /* guint8 ChangeState, indexed by line */
	GArray *del_state;              /* guint8 DelState, deletion before line i (size n+1) */
	guint recompute_id;             /* debounce timeout source */
	GtkSourceGutterRenderer *renderer;  /* not owned (held by the gutter) */
} ChangeTracker;

/* Custom gutter renderer that paints the change bar by reading the
 * ChangeTracker attached to its buffer. */
#define SIRIL_TYPE_CHANGE_GUTTER (siril_change_gutter_get_type())
G_DECLARE_FINAL_TYPE(SirilChangeGutter, siril_change_gutter, SIRIL, CHANGE_GUTTER, GtkSourceGutterRenderer)
struct _SirilChangeGutter { GtkSourceGutterRenderer parent_instance; };
G_DEFINE_TYPE(SirilChangeGutter, siril_change_gutter, GTK_SOURCE_TYPE_GUTTER_RENDERER)

static void siril_change_gutter_snapshot_line(GtkSourceGutterRenderer *renderer,
		GtkSnapshot *snapshot, GtkSourceGutterLines *lines, guint line) {
	GtkTextBuffer *buffer = gtk_source_gutter_lines_get_buffer(lines);
	if (!buffer)
		return;
	ChangeTracker *ct = g_object_get_data(G_OBJECT(buffer), "change-tracker");
	if (!ct)
		return;

	int width = gtk_widget_get_width(GTK_WIDGET(renderer));
	int y = 0, height = 0;
	/* Line y/height in the gutter's window coordinates, via the text view.
	 * (gtk_source_gutter_lines_get_line_yrange is deprecated in newer
	 * gtksourceview; this path is version-independent.) */
	GtkTextView *view = gtk_source_gutter_lines_get_view(lines);
	GtkTextIter it;
	gtk_text_buffer_get_iter_at_line(buffer, &it, (int) line);
	int by = 0;
	gtk_text_view_get_line_yrange(view, &it, &by, &height);
	gtk_text_view_buffer_to_window_coords(view, GTK_TEXT_WINDOW_LEFT, 0, by, NULL, &y);

	static const GdkRGBA amber = { 0.95f, 0.61f, 0.07f, 1.0f };
	static const GdkRGBA green = { 0.18f, 0.70f, 0.32f, 1.0f };
	static const GdkRGBA red   = { 0.86f, 0.20f, 0.18f, 1.0f };

	guint8 st = (line < ct->line_state->len)
		? g_array_index(ct->line_state, guint8, line) : CHG_NONE;
	if (st == CHG_UNSAVED || st == CHG_SAVED) {
		graphene_rect_t bar = GRAPHENE_RECT_INIT((float) width - 4.0f, (float) y,
		                                         3.0f, (float) height);
		gtk_snapshot_append_color(snapshot,
			st == CHG_UNSAVED ? &amber : &green, &bar);
	}

	guint8 d = (line < ct->del_state->len)
		? g_array_index(ct->del_state, guint8, line) : DEL_NONE;
	if (d != DEL_NONE) {
		const GdkRGBA *c = (d == DEL_UNSAVED) ? &red : &green;
		/* Deletion caret: a small triangle at the top edge of the line,
		 * pointing down to indicate lines were removed above this one. */
		graphene_rect_t area = GRAPHENE_RECT_INIT(0.0f, (float) y,
		                                          (float) width, 7.0f);
		cairo_t *cr = gtk_snapshot_append_cairo(snapshot, &area);
		cairo_set_source_rgba(cr, c->red, c->green, c->blue, c->alpha);
		cairo_move_to(cr, width - 7.0, y);
		cairo_line_to(cr, width, y);
		cairo_line_to(cr, width - 3.5, y + 6.0);
		cairo_close_path(cr);
		cairo_fill(cr);
		cairo_destroy(cr);
	}
}

static void siril_change_gutter_class_init(SirilChangeGutterClass *klass) {
	GTK_SOURCE_GUTTER_RENDERER_CLASS(klass)->snapshot_line =
		siril_change_gutter_snapshot_line;
}
static void siril_change_gutter_init(SirilChangeGutter *self) { (void) self; }

/* Split the buffer into lines; *n receives the line count.  The returned
 * array is NULL-terminated (g_strfreev to release).  Index i corresponds to
 * GtkTextBuffer line i. */
static gchar **change_split_buffer(GtkTextBuffer *buffer, guint *n) {
	GtkTextIter s, e;
	gtk_text_buffer_get_bounds(buffer, &s, &e);
	gchar *text = gtk_text_buffer_get_text(buffer, &s, &e, FALSE);
	gchar **lines = g_strsplit(text, "\n", -1);
	g_free(text);
	*n = g_strv_length(lines);
	return lines;
}

/* Line-based LCS diff of current[] against base[].
 * cur_changed[i] (size n) is set TRUE when current line i is not part of the
 * longest common subsequence (i.e. added or modified relative to base).
 * del_before[i] (size n+1) counts base lines deleted at the gap before
 * current line i (del_before[n] = trailing deletions). */
static void change_diff(gchar **cur, guint n, gchar **base, guint m,
                        gboolean *cur_changed, guint *del_before) {
	for (guint i = 0; i < n; i++) cur_changed[i] = TRUE;
	for (guint i = 0; i <= n; i++) del_before[i] = 0;
	if (m == 0)
		return;                  /* nothing to match: all current lines added */
	if (n == 0) {
		del_before[0] = m;       /* whole baseline deleted */
		return;
	}

	/* Guard the O(n*m) table against pathologically large buffers. */
	guint64 cells = (guint64)(n + 1) * (guint64)(m + 1);
	if (cells > 6000000ULL) {
		guint k = MIN(n, m);
		for (guint i = 0; i < k; i++)
			cur_changed[i] = (g_strcmp0(cur[i], base[i]) != 0);
		if (m > n)
			del_before[n] = m - n;
		return;
	}

	guint *dp = g_new0(guint, (gsize) cells);
#define DP(i,j) dp[(gsize)(i) * (m + 1) + (j)]
	for (gint i = (gint) n - 1; i >= 0; i--)
		for (gint j = (gint) m - 1; j >= 0; j--)
			DP(i, j) = (g_strcmp0(cur[i], base[j]) == 0)
				? DP(i + 1, j + 1) + 1
				: MAX(DP(i + 1, j), DP(i, j + 1));

	guint i = 0, j = 0;
	while (i < n && j < m) {
		if (g_strcmp0(cur[i], base[j]) == 0) {
			cur_changed[i] = FALSE; i++; j++;
		} else if (DP(i + 1, j) >= DP(i, j + 1)) {
			i++;                 /* current line changed/added */
		} else {
			del_before[i]++; j++; /* base line deleted here */
		}
	}
	while (j < m) { del_before[n]++; j++; }
#undef DP
	g_free(dp);
}

static gboolean change_recompute_idle(gpointer data) {
	ChangeTracker *ct = data;
	ct->recompute_id = 0;

	guint n = 0;
	gchar **cur = change_split_buffer(ct->buffer, &n);

	g_array_set_size(ct->line_state, n);
	g_array_set_size(ct->del_state, n + 1);
	if (n) memset(ct->line_state->data, 0, n);
	memset(ct->del_state->data, 0, n + 1);

	gboolean *chg_s = g_new0(gboolean, n ? n : 1);
	guint    *del_s = g_new0(guint, n + 1);
	change_diff(cur, n, ct->save_lines, ct->save_n, chg_s, del_s);

	gboolean *chg_o = g_new0(gboolean, n ? n : 1);
	guint    *del_o = g_new0(guint, n + 1);
	change_diff(cur, n, ct->open_lines, ct->open_n, chg_o, del_o);

	for (guint i = 0; i < n; i++) {
		guint8 st = CHG_NONE;
		if (chg_s[i])           st = CHG_UNSAVED;   /* differs from last save */
		else if (chg_o[i])      st = CHG_SAVED;     /* saved, but changed this session */
		g_array_index(ct->line_state, guint8, i) = st;
	}
	/* A deletion gap that is adjacent to a changed/added current line is part
	 * of a modification (the amber/green bar already conveys it) — only show a
	 * caret for a "pure" deletion, where no current line at the gap changed. */
	for (guint i = 0; i <= n; i++) {
		guint8 d = DEL_NONE;
		gboolean pure_s = (del_s[i] > 0) &&
			(i == 0 || !chg_s[i - 1]) && (i >= n || !chg_s[i]);
		gboolean pure_o = (del_o[i] > 0) &&
			(i == 0 || !chg_o[i - 1]) && (i >= n || !chg_o[i]);
		if (pure_s)             d = DEL_UNSAVED;
		else if (pure_o)        d = DEL_SAVED;
		g_array_index(ct->del_state, guint8, i) = d;
	}

	g_free(chg_s); g_free(del_s); g_free(chg_o); g_free(del_o);
	g_strfreev(cur);

	if (ct->renderer)
		gtk_widget_queue_draw(GTK_WIDGET(ct->renderer));
	return G_SOURCE_REMOVE;
}

static void change_tracker_schedule(ChangeTracker *ct) {
	if (ct && ct->recompute_id == 0)
		ct->recompute_id = g_timeout_add(150, change_recompute_idle, ct);
}

/* Snapshot the current buffer as the new save baseline (and optionally the
 * open baseline too).  Reset open+save on load/new; advance only save on a
 * successful save so amber marks flip to green. */
static void change_tracker_set_baselines(ChangeTracker *ct, gboolean reset_open) {
	if (!ct) return;
	guint n = 0;
	gchar **cur = change_split_buffer(ct->buffer, &n);
	g_strfreev(ct->save_lines);
	ct->save_lines = g_strdupv(cur);
	ct->save_n = n;
	if (reset_open) {
		g_strfreev(ct->open_lines);
		ct->open_lines = g_strdupv(cur);
		ct->open_n = n;
	}
	g_strfreev(cur);
	change_tracker_schedule(ct);
}

static void change_tracker_free(gpointer data) {
	ChangeTracker *ct = data;
	if (ct->recompute_id)
		g_source_remove(ct->recompute_id);
	g_strfreev(ct->open_lines);
	g_strfreev(ct->save_lines);
	g_array_unref(ct->line_state);
	g_array_unref(ct->del_state);
	g_free(ct);
}

static void on_buffer_changed_for_tracking(GtkTextBuffer *buffer, gpointer ud) {
	(void) ud;
	change_tracker_schedule(g_object_get_data(G_OBJECT(buffer), "change-tracker"));
}

/* Build the change tracker for a view/buffer pair and attach it to the
 * buffer.  Inserts the change-bar gutter renderer at the far left, ahead of
 * the line-number column. */
static ChangeTracker *change_tracker_setup(GtkSourceView *view, GtkSourceBuffer *buffer) {
	ChangeTracker *ct = g_new0(ChangeTracker, 1);
	ct->buffer = GTK_TEXT_BUFFER(buffer);
	ct->line_state = g_array_new(FALSE, TRUE, sizeof(guint8));
	ct->del_state  = g_array_new(FALSE, TRUE, sizeof(guint8));

	ct->renderer = GTK_SOURCE_GUTTER_RENDERER(g_object_new(SIRIL_TYPE_CHANGE_GUTTER, NULL));
	gtk_widget_set_size_request(GTK_WIDGET(ct->renderer), 8, -1);
	GtkSourceGutter *gutter = gtk_source_view_get_gutter(view, GTK_TEXT_WINDOW_LEFT);
	gtk_source_gutter_insert(gutter, ct->renderer, 0);

	g_object_set_data_full(G_OBJECT(buffer), "change-tracker", ct, change_tracker_free);
	g_signal_connect(buffer, "changed", G_CALLBACK(on_buffer_changed_for_tracking), NULL);

	change_tracker_set_baselines(ct, TRUE);
	return ct;
}

/* Reset both baselines to the current buffer content (load / new / clear). */
static void change_tracker_reset_current(GtkSourceBuffer *buffer) {
	change_tracker_set_baselines(
		g_object_get_data(G_OBJECT(buffer), "change-tracker"), TRUE);
}

/* Advance only the save baseline (successful save): amber → green. */
static void change_tracker_mark_saved(GtkSourceBuffer *buffer) {
	change_tracker_set_baselines(
		g_object_get_data(G_OBJECT(buffer), "change-tracker"), FALSE);
}

/* Build the minimap (GtkSourceMap), wrapped in a clipping scrolled window.
 * Reverted from the experimental drawing-area minimap.  Appended to
 * codeviewbox later (after setup_find_overlay re-parents the editor scroller)
 * so the minimap stays on the right-hand side. */
static void add_minimap(void) {
	map = GTK_SOURCE_MAP(gtk_source_map_new());

	/* Tiny font so the minimap reads as a thumbnail rather than a second
	 * editor.  GtkSourceMap IS the textview CSS node (inherits GtkSourceView
	 * -> GtkTextView), so target the widget itself and its text part. */
	siril_register_css_for_display("siril-pyeditor-map",
		".pyeditor-map, .pyeditor-map text { font-size: 1px; }");
	gtk_widget_add_css_class(GTK_WIDGET(map), "pyeditor-map");

	gtk_widget_set_visible(GTK_WIDGET(map), TRUE);
	gtk_source_map_set_view(map, code_view);
	gtk_widget_set_valign(GTK_WIDGET(map), GTK_ALIGN_FILL);
	gtk_widget_set_hexpand(GTK_WIDGET(map), FALSE);

	/* Wrap in a scrolled window with propagate-natural-width = FALSE so the
	 * map's font-driven natural width can't blow up the layout; the inner map
	 * is clipped (no scrollbars). */
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
}

/* ---------------------------------------------------------------------------
 * Dynamic line wrap (Kate-style indent-aligned wrapping)
 *
 * When enabled, long lines wrap (GTK_WRAP_WORD_CHAR) and each wrapped row is
 * aligned to its line's indentation via a per-line GtkTextTag with a negative
 * `indent` (a hanging indent: first row at x=0, wrapped rows at the indent
 * depth).  A gutter renderer draws a continuation marker (↪) on wrapped rows.
 *
 * State hangs off the GtkSourceBuffer (like the change tracker) so it stays
 * per-document.
 * ------------------------------------------------------------------------- */

typedef struct {
	GtkSourceView *view;            /* not owned */
	gboolean enabled;
	int char_width;                 /* pixels per monospace column */
	GHashTable *tags;               /* int cols -> GtkTextTag* (owned by tag table) */
	GtkSourceGutterRenderer *renderer;  /* continuation marker (not owned) */
} WrapData;

static WrapData *wrap_data(GtkTextBuffer *buffer) {
	return buffer ? g_object_get_data(G_OBJECT(buffer), "wrap-data") : NULL;
}

/* Pixel advance of one monospace column in the view's current font. */
static int wrap_measure_char_width(GtkSourceView *view) {
	PangoLayout *layout = gtk_widget_create_pango_layout(GTK_WIDGET(view), "0000000000");
	int w = 0, h = 0;
	pango_layout_get_pixel_size(layout, &w, &h);
	g_object_unref(layout);
	int cw = w / 10;
	return cw > 0 ? cw : 8;
}

/* Visual columns of leading whitespace on a line (tabs expand to tab stops). */
static int wrap_leading_cols(GtkSourceView *view, const char *line) {
	int tabw = gtk_source_view_get_tab_width(view);
	if (tabw <= 0) tabw = 4;
	int col = 0;
	for (const char *c = line; *c; c++) {
		if (*c == ' ') col++;
		else if (*c == '\t') col += tabw - (col % tabw);
		else break;
	}
	return col;
}

static GtkTextTag *wrap_tag_for_cols(WrapData *wd, GtkTextBuffer *buf, int cols) {
	if (cols <= 0)
		return NULL;
	GtkTextTag *t = g_hash_table_lookup(wd->tags, GINT_TO_POINTER(cols));
	if (!t) {
		int px = cols * wd->char_width;
		/* negative indent + zero left-margin == hanging indent */
		t = gtk_text_buffer_create_tag(buf, NULL, "indent", -px, "left-margin", 0, NULL);
		g_hash_table_insert(wd->tags, GINT_TO_POINTER(cols), t);
	}
	return t;
}

/* Re-apply the hanging-indent tag to buffer lines [first..last]. */
static void wrap_apply_lines(WrapData *wd, GtkTextBuffer *buf, int first, int last) {
	if (!wd->enabled)
		return;
	int n = gtk_text_buffer_get_line_count(buf);
	if (first < 0) first = 0;
	if (last > n - 1) last = n - 1;
	for (int ln = first; ln <= last; ln++) {
		GtkTextIter s, e;
		gtk_text_buffer_get_iter_at_line(buf, &s, ln);
		e = s;
		if (!gtk_text_iter_forward_line(&e))
			gtk_text_buffer_get_end_iter(buf, &e);
		/* clear any of our existing wrap tags on the line first */
		GHashTableIter it; gpointer k, v;
		g_hash_table_iter_init(&it, wd->tags);
		while (g_hash_table_iter_next(&it, &k, &v))
			gtk_text_buffer_remove_tag(buf, GTK_TEXT_TAG(v), &s, &e);
		GtkTextIter le = s;
		gtk_text_iter_forward_to_line_end(&le);
		gchar *txt = gtk_text_buffer_get_text(buf, &s, &le, FALSE);
		int cols = wrap_leading_cols(wd->view, txt);
		g_free(txt);
		GtkTextTag *t = wrap_tag_for_cols(wd, buf, cols);
		if (t)
			gtk_text_buffer_apply_tag(buf, t, &s, &e);
	}
}

static void wrap_apply_all(WrapData *wd, GtkTextBuffer *buf) {
	wrap_apply_lines(wd, buf, 0, gtk_text_buffer_get_line_count(buf) - 1);
}

static void wrap_on_insert(GtkTextBuffer *b, GtkTextIter *loc, char *text, int len, gpointer u) {
	(void) len;
	WrapData *wd = u;
	if (!wd->enabled)
		return;
	int end_line = gtk_text_iter_get_line(loc);  /* loc is end of inserted text */
	int start_line = end_line;
	for (char *c = text; *c; c++)
		if (*c == '\n') start_line--;
	wrap_apply_lines(wd, b, start_line, end_line);
	if (wrap_guides_area)
		gtk_widget_queue_draw(wrap_guides_area);
}

static void wrap_on_delete(GtkTextBuffer *b, GtkTextIter *s, GtkTextIter *e, gpointer u) {
	(void) e;
	WrapData *wd = u;
	if (!wd->enabled)
		return;
	wrap_apply_lines(wd, b, gtk_text_iter_get_line(s), gtk_text_iter_get_line(s));
	if (wrap_guides_area)
		gtk_widget_queue_draw(wrap_guides_area);
}

/* Continuation-marker gutter renderer: draws ↪ on wrapped rows. */
#define SIRIL_TYPE_CONT_GUTTER (siril_cont_gutter_get_type())
G_DECLARE_FINAL_TYPE(SirilContGutter, siril_cont_gutter, SIRIL, CONT_GUTTER, GtkSourceGutterRenderer)
struct _SirilContGutter { GtkSourceGutterRenderer parent_instance; };
G_DEFINE_TYPE(SirilContGutter, siril_cont_gutter, GTK_SOURCE_TYPE_GUTTER_RENDERER)

static void siril_cont_gutter_snapshot_line(GtkSourceGutterRenderer *renderer,
		GtkSnapshot *snapshot, GtkSourceGutterLines *lines, guint line) {
	GtkTextBuffer *buf = gtk_source_gutter_lines_get_buffer(lines);
	WrapData *wd = wrap_data(buf);
	if (!wd || !wd->enabled)
		return;
	GtkTextView *view = gtk_source_gutter_lines_get_view(lines);
	GtkTextIter it;
	gtk_text_buffer_get_iter_at_line(buf, &it, (int) line);
	GdkRectangle loc;
	gtk_text_view_get_iter_location(view, &it, &loc);  /* loc.height = one row */
	int by = 0, full = 0;
	gtk_text_view_get_line_yrange(view, &it, &by, &full);
	int rowh = loc.height > 0 ? loc.height : full;
	int nrows = (rowh > 0) ? (full / rowh) : 1;
	if (nrows <= 1)
		return;
	int wy = 0;
	gtk_text_view_buffer_to_window_coords(view, GTK_TEXT_WINDOW_LEFT, 0, by, NULL, &wy);
	int width = gtk_widget_get_width(GTK_WIDGET(renderer));
	/* Draw a hooked continuation arrow (↪) with cairo so it doesn't depend on
	 * a font glyph being present. */
	graphene_rect_t cell = GRAPHENE_RECT_INIT(0.0f, (float) wy,
	                                          (float) width, (float) full);
	cairo_t *cr = gtk_snapshot_append_cairo(snapshot, &cell);
	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.9);
	cairo_set_line_width(cr, 1.2);
	for (int k = 1; k < nrows; k++) {
		double cx = width - 11.0;
		double yt = wy + k * rowh + rowh * 0.30;
		double yb = wy + k * rowh + rowh * 0.62;
		cairo_move_to(cr, cx, yt);          /* down ... */
		cairo_line_to(cr, cx, yb);          /* ... then right (the hook) */
		cairo_line_to(cr, cx + 6.0, yb);
		cairo_stroke(cr);
		cairo_move_to(cr, cx + 3.0, yb - 3.0);  /* arrowhead */
		cairo_line_to(cr, cx + 6.5, yb);
		cairo_line_to(cr, cx + 3.0, yb + 3.0);
		cairo_stroke(cr);
	}
	cairo_destroy(cr);
}
static void siril_cont_gutter_class_init(SirilContGutterClass *klass) {
	GTK_SOURCE_GUTTER_RENDERER_CLASS(klass)->snapshot_line = siril_cont_gutter_snapshot_line;
}
static void siril_cont_gutter_init(SirilContGutter *self) { (void) self; }

static void wrap_data_free(gpointer p) {
	WrapData *wd = p;
	g_hash_table_destroy(wd->tags);
	g_free(wd);
}

static void wrap_setup(GtkSourceView *view, GtkSourceBuffer *buffer) {
	WrapData *wd = g_new0(WrapData, 1);
	wd->view = view;
	wd->char_width = 8;
	wd->tags = g_hash_table_new(g_direct_hash, g_direct_equal);
	g_object_set_data_full(G_OBJECT(buffer), "wrap-data", wd, wrap_data_free);
	g_signal_connect_after(buffer, "insert-text", G_CALLBACK(wrap_on_insert), wd);
	g_signal_connect_after(buffer, "delete-range", G_CALLBACK(wrap_on_delete), wd);

	wd->renderer = GTK_SOURCE_GUTTER_RENDERER(g_object_new(SIRIL_TYPE_CONT_GUTTER, NULL));
	gtk_widget_set_size_request(GTK_WIDGET(wd->renderer), 18, -1);
	GtkSourceGutter *gutter = gtk_source_view_get_gutter(view, GTK_TEXT_WINDOW_LEFT);
	/* High position so the marker sits to the right of the line-number column
	 * (nearest the text), where a wrapped row's blank line-number would be. */
	gtk_source_gutter_insert(gutter, wd->renderer, 100);
}

/* Toggle dynamic wrap: set the wrap mode and (re)apply or clear the
 * hanging-indent tags. */
static void wrap_set_enabled(gboolean on) {
	WrapData *wd = wrap_data(GTK_TEXT_BUFFER(sourcebuffer));
	if (!wd)
		return;
	wd->enabled = on;
	GtkTextBuffer *buf = GTK_TEXT_BUFFER(sourcebuffer);
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(code_view),
		on ? GTK_WRAP_WORD_CHAR : GTK_WRAP_NONE);
	if (on) {
		wd->char_width = wrap_measure_char_width(code_view);
		wrap_apply_all(wd, buf);
	} else {
		GtkTextIter s, e;
		gtk_text_buffer_get_bounds(buf, &s, &e);
		GHashTableIter it; gpointer k, v;
		g_hash_table_iter_init(&it, wd->tags);
		while (g_hash_table_iter_next(&it, &k, &v))
			gtk_text_buffer_remove_tag(buf, GTK_TEXT_TAG(v), &s, &e);
	}
	if (wd->renderer)
		gtk_widget_queue_draw(GTK_WIDGET(wd->renderer));
	if (wrap_guides_area)
		gtk_widget_queue_draw(wrap_guides_area);
}

/* EXPERIMENTAL (Kate-style): paint grey tab/indent guides in the virtual
 * hanging-indent area of wrapped continuation rows.  Drawn on a transparent
 * overlay above the editor; only active while dynamic wrap is enabled. */
static void wrap_guides_draw(GtkDrawingArea *area, cairo_t *cr, int w, int h, gpointer u) {
	(void) area; (void) w; (void) h; (void) u;
	if (!code_view || !sourcebuffer)
		return;
	WrapData *wd = wrap_data(GTK_TEXT_BUFFER(sourcebuffer));
	if (!wd || !wd->enabled)
		return;
	GtkTextView *view = GTK_TEXT_VIEW(code_view);
	GtkTextBuffer *buf = GTK_TEXT_BUFFER(sourcebuffer);
	int cw = wd->char_width > 0 ? wd->char_width : 8;

	GdkRectangle vis;
	gtk_text_view_get_visible_rect(view, &vis);
	GtkTextIter top, bot;
	gtk_text_view_get_line_at_y(view, &top, vis.y, NULL);
	gtk_text_view_get_line_at_y(view, &bot, vis.y + vis.height, NULL);
	int first = gtk_text_iter_get_line(&top);
	int last  = gtk_text_iter_get_line(&bot);

	/* widget-space x of buffer x=0 (the text's left edge, after the gutter) */
	int basex = 0, basetmp = 0;
	gtk_text_view_buffer_to_window_coords(view, GTK_TEXT_WINDOW_WIDGET, 0, 0, &basex, &basetmp);

	cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.6);
	cairo_set_line_width(cr, 1.0);

	for (int ln = first; ln <= last; ln++) {
		GtkTextIter s, e;
		gtk_text_buffer_get_iter_at_line(buf, &s, ln);
		e = s;
		gtk_text_iter_forward_to_line_end(&e);
		gchar *txt = gtk_text_buffer_get_text(buf, &s, &e, FALSE);
		int cols = wrap_leading_cols(GTK_SOURCE_VIEW(code_view), txt);
		g_free(txt);
		if (cols <= 0)
			continue;
		GdkRectangle loc;
		gtk_text_view_get_iter_location(view, &s, &loc);
		int rowh = loc.height;
		int by = 0, full = 0;
		gtk_text_view_get_line_yrange(view, &s, &by, &full);
		if (rowh <= 0)
			continue;
		int nrows = full / rowh;
		if (nrows <= 1)
			continue;
		int lwx = 0, lwy = 0;
		gtk_text_view_buffer_to_window_coords(view, GTK_TEXT_WINDOW_WIDGET, 0, by, &lwx, &lwy);
		for (int k = 1; k < nrows; k++) {
			int rowtop = lwy + k * rowh;
			/* one solid grey block spanning the whole indent depth */
			double x0 = basex;
			double bw = cols * cw;
			double pad_y = rowh * 0.18;
			cairo_rectangle(cr, x0, rowtop + pad_y, bw, rowh - 2.0 * pad_y);
			cairo_fill(cr);
		}
	}
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

	// Attach the changed-line tracker / change-bar gutter renderer
	change_tracker_setup(code_view, sourcebuffer);

	// Dynamic line-wrap manager (indent-aligned continuations + marker)
	wrap_setup(code_view, sourcebuffer);

	// Build the custom minimap (decoupled from the editor buffer)
	add_minimap();
}

/** Code for the find feature ***/

static void setup_find_overlay() {
	// Create a new overlay
	GtkWidget *new_overlay = gtk_overlay_new();
	content_overlay = new_overlay;   /* the active editor's find/guides host */
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

	/* EXPERIMENTAL: transparent overlay for Kate-style wrap indent guides.
	 * Non-interactive (can-target FALSE) so it never steals editor input; it
	 * only paints when dynamic wrap is enabled. */
	wrap_guides_area = gtk_drawing_area_new();
	gtk_widget_set_can_target(wrap_guides_area, FALSE);
	gtk_widget_set_hexpand(wrap_guides_area, TRUE);
	gtk_widget_set_vexpand(wrap_guides_area, TRUE);
	gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(wrap_guides_area),
	                               wrap_guides_draw, NULL, NULL);
	gtk_overlay_add_overlay(GTK_OVERLAY(new_overlay), wrap_guides_area);
	{
		GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(scrolled_window);
		if (vadj)
			g_signal_connect_swapped(vadj, "value-changed",
				G_CALLBACK(gtk_widget_queue_draw), wrap_guides_area);
	}

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
	(void) data;
	SearchData *search_data = active_search;   /* operate on the active tab */
	if (!search_data)
		return FALSE;
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
	(void) button; (void) search_data;
	if (active_search) forward_search(active_search);
}

void on_previous_match(GtkButton *button, SearchData *search_data) {
	(void) button; (void) search_data;
	if (active_search) backward_search(active_search);
}

void on_find_entry_changed(GtkSearchEntry *entry, SearchData *search_data) {
    (void) entry; (void) search_data;
    if (active_search) perform_search(active_search);
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

	/* The shared find widgets (find_entry, go_up/down) are connected once in
	 * setup_find_signals(); those handlers act on `active_search`.  We only
	 * create and attach the per-view SearchData here. */
	g_object_set_data_full(G_OBJECT(source_view), "search-data", search_data, g_free);
}

/* Connect the shared find widgets once.  The handlers operate on the active
 * tab's SearchData (active_search), which is updated on tab switch. */
static void setup_find_signals(void) {
	static gboolean done = FALSE;
	if (done || !find_entry)
		return;
	done = TRUE;
	g_signal_connect(find_entry, "changed", G_CALLBACK(on_find_entry_changed), NULL);
	GtkEventController *kctrl = gtk_event_controller_key_new();
	g_signal_connect(kctrl, "key-pressed", G_CALLBACK(on_find_entry_key_press_event), NULL);
	gtk_widget_add_controller(GTK_WIDGET(find_entry), kctrl);
	if (go_down_button)
		g_signal_connect(go_down_button, "clicked", G_CALLBACK(on_next_match), NULL);
	if (go_up_button)
		g_signal_connect(go_up_button, "clicked", G_CALLBACK(on_previous_match), NULL);
}

/* Apply every editor preference to the active document's source view (and sync
 * the menu check states).  Shared by the per-show path and new-tab creation so
 * new tabs match the current settings. */
static void apply_editor_settings_to_active(void) {
	if (!code_view)
		return;
	gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);
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
		{ "dynamicwrap",          com.pref.gui.editor_cfg.dynamic_wrap },
		{ "useargs",              from_cli },
		{ "pythondebug",          python_debug },
	};
	for (size_t i = 0; i < G_N_ELEMENTS(setup); i++) {
		apply_editor_setting(setup[i].name, setup[i].v);
		editor_action_set_state(setup[i].name, g_variant_new_boolean(setup[i].v));
	}
	editor_action_set_state("language",
		g_variant_new_string(active_language == LANG_PYTHON ? "py" : "ssf"));
}

/* ---- Multi-tab document management ----
 * The module statics mirror the active document.  doc_save_statics copies them
 * into a doc; doc_load_statics copies a doc back into the statics. */
static void doc_save_statics(EditorDoc *d) {
	if (!d) return;
	d->view = code_view; d->buffer = sourcebuffer; d->scroll = scrolled_window;
	d->map = map; d->map_wrapper = map_wrapper; d->guides = wrap_guides_area;
	d->space_drawer = space_drawer; d->overlay = content_overlay;
	d->current_file = current_file; d->active_language = active_language;
	d->buffer_modified = buffer_modified;
}
static void doc_load_statics(EditorDoc *d) {
	code_view = d->view; sourcebuffer = d->buffer; scrolled_window = d->scroll;
	map = d->map; map_wrapper = d->map_wrapper; wrap_guides_area = d->guides;
	space_drawer = d->space_drawer; content_overlay = d->overlay;
	current_file = d->current_file; active_language = d->active_language;
	buffer_modified = d->buffer_modified;
}

/* Apply an editor preference to EVERY open tab (view prefs are global to the
 * script editor), then restore the active document's statics. */
static void apply_editor_setting_all(const gchar *name, gboolean v) {
	if (!editor_docs || editor_docs->len <= 1) {
		apply_editor_setting(name, v);
		return;
	}
	EditorDoc *cur = active_doc;
	if (cur) doc_save_statics(cur);
	for (guint i = 0; i < editor_docs->len; i++) {
		EditorDoc *d = g_ptr_array_index(editor_docs, i);
		doc_load_statics(d);
		apply_editor_setting(name, v);
		doc_save_statics(d);
	}
	if (cur) doc_load_statics(cur);
}

/* Refresh both the window title and the active tab's label from the active
 * document's file + modified state. */
static void editor_update_window_title(void) {
	char *base = current_file ? g_file_get_basename(current_file) : g_strdup(_("unsaved"));
	char *t = g_strdup_printf("%s%s", base, buffer_modified ? "*" : "");
	if (editor_window)
		gtk_window_set_title(editor_window, t);
	if (active_doc && active_doc->tab_label)
		gtk_label_set_text(GTK_LABEL(active_doc->tab_label), t);
	g_free(t);
	g_free(base);
}

/* Make `d` the active document: load its statics, move the shared find box into
 * its overlay, and point active_search at its view's search data. */
static void editor_make_active(EditorDoc *d) {
	if (!d) return;
	active_doc = d;
	doc_load_statics(d);
	if (find_overlay && content_overlay) {
		GtkWidget *par = gtk_widget_get_parent(GTK_WIDGET(find_overlay));
		if (par != content_overlay) {
			g_object_ref(find_overlay);
			if (GTK_IS_OVERLAY(par))
				gtk_overlay_remove_overlay(GTK_OVERLAY(par), GTK_WIDGET(find_overlay));
			else if (par)
				gtk_widget_unparent(GTK_WIDGET(find_overlay));
			gtk_overlay_add_overlay(GTK_OVERLAY(content_overlay), GTK_WIDGET(find_overlay));
			g_object_unref(find_overlay);
		}
	}
	active_search = code_view ? g_object_get_data(G_OBJECT(code_view), "search-data") : NULL;
	if (active_search) active_search->search_entry = find_entry;
	gtk_label_set_text(language_label,
		active_language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));
	editor_update_window_title();
}

static EditorDoc *editor_doc_for_page(GtkWidget *page) {
	if (!editor_docs) return NULL;
	for (guint i = 0; i < editor_docs->len; i++) {
		EditorDoc *d = g_ptr_array_index(editor_docs, i);
		if (d->page == page) return d;
	}
	return NULL;
}

static void on_notebook_switch_page(GtkNotebook *nb, GtkWidget *page, guint idx, gpointer u) {
	(void) nb; (void) idx; (void) u;
	EditorDoc *d = editor_doc_for_page(page);
	if (!d || d == active_doc)
		return;
	if (active_doc) doc_save_statics(active_doc);
	editor_make_active(d);
}

static GtkWidget *editor_make_tab_label(EditorDoc *d) {
	GtkWidget *lbl = gtk_label_new(_("unsaved"));
	d->tab_label = lbl;
	return lbl;
}

/* Build a fresh empty document in a new notebook tab and make it active. */
static EditorDoc *editor_new_tab(void) {
	if (!editor_docs) editor_docs = g_ptr_array_new();
	if (active_doc) doc_save_statics(active_doc);   /* preserve outgoing doc */

	GtkSourceView *view = GTK_SOURCE_VIEW(gtk_source_view_new());
	gtk_widget_set_focusable(GTK_WIDGET(view), TRUE);
	GtkWidget *scroll = gtk_scrolled_window_new();
	gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(scroll), 600);
	gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(scroll), 400);
	gtk_widget_set_hexpand(scroll, TRUE);
	gtk_widget_set_vexpand(scroll, TRUE);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroll), GTK_WIDGET(view));

	/* repoint statics at the new widgets, then build per-doc machinery */
	code_view = view;
	scrolled_window = GTK_SCROLLED_WINDOW(scroll);
	current_file = NULL;
	active_language = LANG_PYTHON;
	buffer_modified = FALSE;

	add_code_view(NULL);   /* sourcebuffer + map + change tracker + wrap */
	g_signal_connect(sourcebuffer, "modified-changed",
	                 G_CALLBACK(on_buffer_modified_changed), NULL);

	GtkWidget *overlay = gtk_overlay_new();
	content_overlay = overlay;
	gtk_widget_set_hexpand(overlay, TRUE);
	gtk_widget_set_vexpand(overlay, TRUE);
	gtk_overlay_set_child(GTK_OVERLAY(overlay), scroll);
	wrap_guides_area = gtk_drawing_area_new();
	gtk_widget_set_can_target(wrap_guides_area, FALSE);
	gtk_widget_set_hexpand(wrap_guides_area, TRUE);
	gtk_widget_set_vexpand(wrap_guides_area, TRUE);
	gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(wrap_guides_area), wrap_guides_draw, NULL, NULL);
	gtk_overlay_add_overlay(GTK_OVERLAY(overlay), wrap_guides_area);
	{
		GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scroll));
		if (vadj)
			g_signal_connect_swapped(vadj, "value-changed",
				G_CALLBACK(gtk_widget_queue_draw), wrap_guides_area);
	}

	GtkWidget *page = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_set_hexpand(page, TRUE);
	gtk_widget_set_vexpand(page, TRUE);
	gtk_box_append(GTK_BOX(page), overlay);
	if (map_wrapper)
		gtk_box_append(GTK_BOX(page), map_wrapper);

	setup_search(view, find_entry);

	EditorDoc *d = g_new0(EditorDoc, 1);
	doc_save_statics(d);
	d->page = page;
	g_ptr_array_add(editor_docs, d);

	int idx = gtk_notebook_append_page(editor_notebook, page, editor_make_tab_label(d));
	gtk_widget_set_visible(page, TRUE);
	/* set active_doc before switching so the switch-page handler no-ops */
	active_doc = d;
	gtk_notebook_set_current_page(editor_notebook, idx);
	editor_make_active(d);
	apply_editor_settings_to_active();
	return d;
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

		/* Append the minimap now, AFTER setup_find_overlay() re-parents the
		 * editor scroller into codeviewbox, so the minimap is the last child
		 * (right-hand side) rather than being pushed to the left. */
		if (map_wrapper && !gtk_widget_get_parent(map_wrapper))
			gtk_box_append(GTK_BOX(codeviewbox), map_wrapper);

		/* Multi-tab (stage 1): wrap the assembled editor content (find-overlay
		 * host + minimap) in a GtkNotebook page.  Subsequent stages add the
		 * per-document factory, tab switching and close handling. */
		editor_notebook = GTK_NOTEBOOK(gtk_notebook_new());
		gtk_widget_set_vexpand(GTK_WIDGET(editor_notebook), TRUE);
		gtk_widget_set_hexpand(GTK_WIDGET(editor_notebook), TRUE);
		gtk_notebook_set_scrollable(editor_notebook, TRUE);
		GtkWidget *page0 = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
		gtk_widget_set_hexpand(page0, TRUE);
		gtk_widget_set_vexpand(page0, TRUE);
		/* Move the existing content (overlay + minimap) out of codeviewbox and
		 * into the first notebook page. */
		if (content_overlay && gtk_widget_get_parent(content_overlay) == GTK_WIDGET(codeviewbox)) {
			g_object_ref(content_overlay);
			gtk_box_remove(GTK_BOX(codeviewbox), content_overlay);
			gtk_box_append(GTK_BOX(page0), content_overlay);
			g_object_unref(content_overlay);
		}
		if (map_wrapper && gtk_widget_get_parent(map_wrapper) == GTK_WIDGET(codeviewbox)) {
			g_object_ref(map_wrapper);
			gtk_box_remove(GTK_BOX(codeviewbox), map_wrapper);
			gtk_box_append(GTK_BOX(page0), map_wrapper);
			g_object_unref(map_wrapper);
		}
		GtkWidget *tab0_label = gtk_label_new(_("unsaved"));
		gtk_notebook_append_page(editor_notebook, page0, tab0_label);
		gtk_box_append(GTK_BOX(codeviewbox), GTK_WIDGET(editor_notebook));

		/* Capture the first document and wire tab switching. */
		setup_find_signals();
		editor_docs = g_ptr_array_new();
		EditorDoc *d0 = g_new0(EditorDoc, 1);
		doc_save_statics(d0);
		d0->page = page0;
		d0->tab_label = tab0_label;
		g_ptr_array_add(editor_docs, d0);
		active_doc = d0;
		active_search = g_object_get_data(G_OBJECT(code_view), "search-data");
		if (active_search) active_search->search_entry = find_entry;
		g_signal_connect(editor_notebook, "switch-page",
		                 G_CALLBACK(on_notebook_switch_page), NULL);

		// Initialize with "unsaved" if no file is loaded
		if (!current_file) {
			gtk_window_set_title(GTK_WINDOW(editor_window), "unsaved");
		}
	}
}

static void update_title_with_modification() {
	editor_update_window_title();
}

static void update_title(GFile *file) {
	if (file) {
		char *basename = g_file_get_basename(file);
		char *suffix = strrchr(basename, '.');
		if (suffix != NULL) {
			if (g_ascii_strcasecmp(suffix, ".py") == 0) {
				editor_action_change_state("language", g_variant_new_string("py"));
			} else if (g_ascii_strcasecmp(suffix, ".ssf") == 0) {
				editor_action_change_state("language", g_variant_new_string("ssf"));
			}
		}
		g_free(basename);
	}
	buffer_modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
	editor_update_window_title();
	if (editor_window)
		gtk_widget_queue_draw(GTK_WIDGET(editor_window));
}

// File handling
void load_file_complete(GObject *loader, GAsyncResult *result, gpointer user_data) {
	GError *error = NULL;
	if (!gtk_source_file_loader_load_finish(GTK_SOURCE_FILE_LOADER(loader), result, &error)) {
		g_printerr("Error loading file: %s\n", error->message);
		gtk_window_set_title(GTK_WINDOW(editor_window), "");
		g_error_free(error);
	} else {
		// Loaded content becomes the baseline: no lines are "changed" yet.
		change_tracker_reset_current(sourcebuffer);
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
	} else {
		// Save succeeded: advance the save baseline so amber marks turn green.
		change_tracker_mark_saved(sourcebuffer);
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
	(void) action; (void) parameter; (void) user_data;
	/* New always opens a fresh, independent tab. */
	editor_new_tab();
	if (code_view)
		gtk_widget_grab_focus(GTK_WIDGET(code_view));
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
	/* Always edit in a fresh, independent tab. */
	editor_new_tab();
	current_file = NULL;
	gtk_text_buffer_begin_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
	gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sourcebuffer), contents, (gint) length);
	buffer_modified = FALSE;
	gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sourcebuffer), FALSE);
	gtk_text_buffer_end_irreversible_action(GTK_TEXT_BUFFER(sourcebuffer));
	change_tracker_reset_current(sourcebuffer);
	if (!g_ascii_strcasecmp(ext, "ssf")) {
		active_language = LANG_SSF;
		editor_action_change_state("language", g_variant_new_string("ssf"));
	} else {
		active_language = LANG_PYTHON;
		editor_action_change_state("language", g_variant_new_string("py"));
	}
	editor_update_window_title();
	gtk_window_present(editor_window);
	gtk_widget_queue_draw(GTK_WIDGET(editor_window));
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
		editor_new_tab();          /* open the file in a fresh tab */
		control_window_switch_to_tab(OUTPUT_LOGS);
		current_file = g_object_ref(file);
		load_file(file);
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
	GFile *file = g_file_new_for_path(path);
	editor_new_tab();              /* open the recent file in a fresh tab */
	control_window_switch_to_tab(OUTPUT_LOGS);
	current_file = g_object_ref(file);
	load_file(file);
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

		if (com.wd) {
			GFile *initial = g_file_new_for_path(com.wd);
			gtk_file_dialog_set_initial_folder(fd, initial);
			g_object_unref(initial);
		}

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

	if (G_IS_OBJECT(current_file)) {
		gtk_file_dialog_set_initial_file(fd, current_file);
	} else if (com.wd) {
		GFile *initial = g_file_new_for_path(com.wd);
		gtk_file_dialog_set_initial_folder(fd, initial);
		g_object_unref(initial);
	}

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

/* Per-show preparation for the script editor.  Connected to
 * python_window's "show" signal in the .ui so it runs no matter HOW
 * the window is made visible — whether via on_open_pythonpad (the
 * normal menu-driven path) or via the new-user Introduction
 * (siril_intro.c) which makes the widget visible directly to anchor
 * a tip popover.  Without this signal-driven path, an intro run
 * before the user has opened the editor manually displayed the
 * window with none of the cfg applied — missing source-view
 * settings, no menu bar, no language label, etc.
 *
 * Idempotent: python_scratchpad_init_statics, setup_editor_actions
 * and populate_pyeditor_recent_menu all guard themselves; the editor
 * toggles are direct pref→widget applies and safe to repeat.  Runs
 * every time the window is shown, which is what we want — closes
 * hide-on-delete, so the next open's pref values are picked up. */
void on_python_window_show(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	python_scratchpad_init_statics();

	/* Decouple window behaviour from main window. */
	gtk_window_set_transient_for(editor_window, NULL);

	/* close-request → hide.  Reconnecting on every show is harmless;
	 * GTK4 ignores duplicate identical g_signal_connect calls? — no,
	 * actually it doesn't.  But python_scratchpad_init_statics is
	 * idempotent because of its NULL guard, and the signal is
	 * connected exactly once there at the bottom of init.  Leaving
	 * it out of the per-show path. */

	setup_editor_actions(editor_window);
	/* The menu bar references this submenu by pointer, so populate it
	 * before showing the window — first open builds the bar, subsequent
	 * opens just refresh the existing one. */
	populate_pyeditor_recent_menu();

	gtk_label_set_text(language_label, active_language == LANG_PYTHON ? _("Python Script") : _("Siril Script File"));

	/* Apply all editor prefs to the active document's view and sync the menu
	 * check states (also used when creating new tabs). */
	apply_editor_settings_to_active();
}

int on_open_pythonpad(GtkWidget *menuitem, gpointer user_data) {
	if (!editor_window) {
		python_scratchpad_init_statics();
	}
	(void) main_window;

	/* gtk_window_present triggers the window's "show" signal when
	 * the window was hidden — on_python_window_show runs there and
	 * does the editor configuration. */
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
		/* Right-margin position is global to the editor: apply to every tab. */
		if (editor_docs) {
			for (guint i = 0; i < editor_docs->len; i++) {
				EditorDoc *d = g_ptr_array_index(editor_docs, i);
				if (d->view)
					gtk_source_view_set_right_margin_position(d->view,
						com.pref.gui.editor_cfg.rmargin_pos);
			}
		} else {
			gtk_source_view_set_right_margin_position(code_view, com.pref.gui.editor_cfg.rmargin_pos);
		}
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
	/* Editor view preferences are global: apply to all open tabs. */
	apply_editor_setting_all(name, v);
}

void on_editor_args_clear_clicked(GtkButton *button, gpointer user_data) {
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(args_entry);
	gtk_entry_buffer_set_text(buffer, "", 0);
	gtk_widget_grab_focus(GTK_WIDGET(args_entry));
}

gboolean get_python_debug_mode() {
	return python_debug;
}
