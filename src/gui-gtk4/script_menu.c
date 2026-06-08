/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifdef _WIN32
#include <windows.h>
#include "core/OS_utils.h"
/* Constant available since Shell32.dll 4.72 */
#ifndef CSIDL_APPDATA
#define CSIDL_APPDATA 0x001a
#endif
#endif
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "io/siril_pythonmodule.h"
#include "io/siril_git.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/python_gui.h"
#include "algos/sorting.h"
#include "script_menu.h"

#define CONFIRM_RUN_SCRIPTS _("You are about to use scripts. Note that scripts execute code with your current user privileges. While Siril Script Files can only execute Siril commands and a very small number of specific external programs, Python scripts are considerably more powerful and execute code not written by the Siril team. Ensure you obtain scripts from a reputable source.")

void on_get_scripts_clicked(gpointer user_data);
static GtkWidget *menuscript = NULL;
static GtkTextView *script_path_textview = NULL;
static GtkWidget *script_tab_conversion = NULL;
static GtkWidget *script_tab_sequence = NULL;
static GtkWidget *script_tab_calibration = NULL;
static GtkWidget *script_tab_registration = NULL;
static GtkWidget *script_tab_plot = NULL;
static GtkWidget *script_tab_stacking = NULL;
static GtkWidget *script_command_entry = NULL;
static GtkWidget *script_notebook1 = NULL;
static GtkWidget *script_headerbar = NULL;
static GtkWidget *script_toolbarbox = NULL;
static GtkStack *script_stack_pref = NULL;
static GtkWidget *script_scripts_page = NULL;

typedef struct {
	gchar *display_name;
	gchar *full_path;
} ScriptEntry;

static GSList *script_search_entries = NULL;
static GtkWidget *script_search_entry_widget = NULL;
static GtkWidget *script_search_listbox = NULL;
static GtkWidget *script_search_scroll = NULL;

static void script_dispatch(const gchar *user_data);

static void free_script_entry(gpointer data) {
	ScriptEntry *e = data;
	if (!e) return;
	g_free(e->display_name);
	g_free(e->full_path);
	g_free(e);
}

static void add_to_script_search_list(const gchar *display_name, const gchar *full_path) {
	ScriptEntry *e = g_new(ScriptEntry, 1);
	e->display_name = g_strdup(display_name);
	e->full_path = g_strdup(full_path);
	script_search_entries = g_slist_prepend(script_search_entries, e);
}

static void on_script_search_row_activated(GtkListBox *box, GtkListBoxRow *row,
                                            gpointer user_data) {
	(void)user_data;
	if (!row) return;
	const gchar *path = g_object_get_data(G_OBJECT(row), "script-path");
	if (!path) return;
	GtkWidget *popover = gtk_widget_get_ancestor(GTK_WIDGET(box), GTK_TYPE_POPOVER);
	if (popover)
		gtk_popover_popdown(GTK_POPOVER(popover));
	script_dispatch(path);
}

static void on_script_search_changed(GtkSearchEntry *entry, gpointer user_data) {
	(void)user_data;
	if (!script_search_listbox || !script_search_scroll)
		return;

	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));

	GtkWidget *child;
	while ((child = gtk_widget_get_first_child(script_search_listbox)) != NULL)
		gtk_list_box_remove(GTK_LIST_BOX(script_search_listbox), child);

	if (!text || !*text) {
		gtk_widget_set_visible(script_search_scroll, FALSE);
		return;
	}

	gchar *needle = g_utf8_casefold(text, -1);
	gint count = 0;

	for (GSList *l = script_search_entries; l; l = l->next) {
		ScriptEntry *e = l->data;
		if (!e) continue;
		gchar *haystack = g_utf8_casefold(e->display_name, -1);
		if (g_strstr_len(haystack, -1, needle)) {
			GtkWidget *label = gtk_label_new(e->display_name);
			gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
			gtk_widget_set_margin_start(label, 6);
			gtk_widget_set_margin_end(label, 6);
			gtk_widget_set_margin_top(label, 3);
			gtk_widget_set_margin_bottom(label, 3);
			GtkWidget *lrow = gtk_list_box_row_new();
			gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(lrow), label);
			gtk_widget_set_tooltip_text(lrow, e->full_path);
			g_object_set_data_full(G_OBJECT(lrow), "script-path",
			                       g_strdup(e->full_path), g_free);
			gtk_list_box_append(GTK_LIST_BOX(script_search_listbox), lrow);
			count++;
		}
		g_free(haystack);
	}

	g_free(needle);
	gtk_widget_set_visible(script_search_scroll, count > 0);
}

static GtkWidget *create_script_search_widget(void) {
	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_widget_set_margin_start(box, 6);
	gtk_widget_set_margin_end(box, 6);
	gtk_widget_set_margin_top(box, 4);
	gtk_widget_set_margin_bottom(box, 6);

	script_search_entry_widget = gtk_search_entry_new();
	gtk_search_entry_set_placeholder_text(GTK_SEARCH_ENTRY(script_search_entry_widget),
	                                      _("Search scripts\xe2\x80\xa6"));
	gtk_box_append(GTK_BOX(box), script_search_entry_widget);

	script_search_listbox = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(script_search_listbox),
	                                GTK_SELECTION_NONE);
	g_signal_connect(script_search_listbox, "row-activated",
	                 G_CALLBACK(on_script_search_row_activated), NULL);

	script_search_scroll = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(script_search_scroll),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_max_content_height(
	    GTK_SCROLLED_WINDOW(script_search_scroll), 250);
	gtk_scrolled_window_set_propagate_natural_height(
	    GTK_SCROLLED_WINDOW(script_search_scroll), TRUE);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(script_search_scroll),
	                              script_search_listbox);
	gtk_widget_set_visible(script_search_scroll, FALSE);
	gtk_box_append(GTK_BOX(box), script_search_scroll);

	g_signal_connect(script_search_entry_widget, "search-changed",
	                 G_CALLBACK(on_script_search_changed), NULL);

	return box;
}

static void script_menu_init_statics(void) {
	if (menuscript) return;
	menuscript = GTK_WIDGET(gtk_builder_get_object(gui.builder, "header_scripts_button"));
	script_path_textview = GTK_TEXT_VIEW(gtk_builder_get_object(gui.builder, "GtkTxtScriptPath"));
	script_tab_conversion = GTK_WIDGET(gtk_builder_get_object(gui.builder, "conversion_tab"));
	script_tab_sequence = GTK_WIDGET(gtk_builder_get_object(gui.builder, "sequence_tab"));
	script_tab_calibration = GTK_WIDGET(gtk_builder_get_object(gui.builder, "calibration_tab"));
	script_tab_registration = GTK_WIDGET(gtk_builder_get_object(gui.builder, "registration_tab"));
	script_tab_plot = GTK_WIDGET(gtk_builder_get_object(gui.builder, "plot_tab"));
	script_tab_stacking = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stacking_tab"));
	script_command_entry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "command"));
	script_notebook1 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "notebook1"));
	script_headerbar = GTK_WIDGET(gtk_builder_get_object(gui.builder, "headerbar"));
	script_toolbarbox = GTK_WIDGET(gtk_builder_get_object(gui.builder, "toolbarbox"));
	script_stack_pref = GTK_STACK(gtk_builder_get_object(gui.builder, "stack_pref"));
	script_scripts_page = GTK_WIDGET(gtk_builder_get_object(gui.builder, "scripts_page"));
}

static GSList *initialize_script_paths(){
	GSList *list = NULL;
#ifdef _WIN32
	list = g_slist_prepend(list, g_build_filename(get_special_folder(CSIDL_APPDATA), "siril", "scripts", NULL));

	gchar *execpath = g_win32_get_package_installation_directory_of_module(NULL);

	list = g_slist_prepend(list, g_build_filename(execpath, "scripts", NULL));
	g_free(execpath);
#else
	list = g_slist_prepend(list, g_build_filename(siril_get_system_data_dir(), "scripts", NULL));
	if (g_getenv("XDG_CONFIG_HOME") != NULL) {
		list = g_slist_prepend(list, g_build_filename(getenv("XDG_CONFIG_HOME"), "scripts", NULL));
	}
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), ".siril", "scripts", NULL));
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), "siril", "scripts", NULL));
#endif
	list = g_slist_reverse(list);
	return list;
}

static void add_path_to_gtkText(gchar *path) {
	script_menu_init_statics();
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(script_path_textview);
	GtkTextIter iter;

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	gtk_text_buffer_insert(tbuf, &iter, path, strlen(path));
	gtk_text_buffer_insert(tbuf, &iter, "\n", strlen("\n"));

	/* scroll to end */
	gtk_text_buffer_get_end_iter(tbuf, &iter);
	GtkTextMark *insert_mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	gtk_text_view_scroll_to_mark(script_path_textview, insert_mark, 0.0, TRUE, 0.0, 1.0);
}

static void clear_gtk_list() {
	script_menu_init_statics();
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(script_path_textview);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
}

void script_widgets_enable(gboolean status) {
	script_menu_init_statics();
	gtk_widget_set_sensitive(script_tab_conversion, status);
	gtk_widget_set_sensitive(script_tab_sequence, status);
	gtk_widget_set_sensitive(script_tab_calibration, status);
	gtk_widget_set_sensitive(script_tab_registration, status);
	gtk_widget_set_sensitive(script_tab_plot, status);
	gtk_widget_set_sensitive(script_tab_stacking, status);
	gtk_widget_set_sensitive(script_command_entry, status);
	gtk_widget_set_sensitive(script_notebook1, status);
	gtk_widget_set_sensitive(script_headerbar, status);
	gtk_widget_set_sensitive(script_toolbarbox, status);
}

gboolean script_widgets_idle(gpointer user_data) {
	script_widgets_enable(TRUE);
	return FALSE;
}

static GSList *search_script(const char *path) {
	GSList *list = NULL;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;

	dir = g_dir_open(path, 0, &error);
	if (!dir) {
		g_print("scripts: %s\n", error->message);
		g_clear_error(&error);
		return NULL;
	}
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_str_has_suffix(file, SCRIPT_EXT) || g_str_has_suffix(file, PYSCRIPT_EXT) || g_str_has_suffix(file, PYCSCRIPT_EXT)) {
			list = g_slist_prepend(list, g_strdup(file));  // Keep the full filename with extension
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	g_dir_close(dir);

	return list;
}

gboolean accept_script_warning_dialog() {
	if (com.pref.gui.warn_scripts_run) {
		gboolean dont_show_again;
		gboolean confirm = siril_confirm_dialog_and_remember(_("Please read me before using scripts"),
				CONFIRM_RUN_SCRIPTS, _("Run Script"), &dont_show_again);
		if (!confirm)
			return FALSE;

		com.pref.gui.warn_scripts_run = !dont_show_again;
		writeinitfile();
	}
	return TRUE;
}

/* Forward decls: defined further down. */
static void script_edit_action_cb(GSimpleAction *action, GVariant *parameter,
                                  gpointer user_data);
static void install_script_popover_rmb(void);

/* Set by the popover's capture-phase press handler; consulted by
 * script_execute_action_cb to decide whether to redirect to script-edit. */
static gint64 script_secondary_press_time = 0;

/* Phase 15: dispatch helper called from the GAction "script-execute"
 * activation handler.  The legacy `on_script_execution` GtkMenuItem
 * trampoline that GTK3 used was deleted — its callers no longer exist
 * (GMenuModel items invoke the GAction directly with the path packed
 * into the action target). */
static void script_dispatch(const gchar *user_data) {
	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (!accept_script_warning_dialog())
		return;

	if (get_script_thread_run())
		wait_for_script_thread();

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	gchar *script_file = g_strdup((gchar *)user_data);

	/* Run the script */
	siril_log_message(_("Starting script %s\n"), script_file);

	if (g_str_has_suffix(script_file, PYSCRIPT_EXT) || g_str_has_suffix(script_file, PYCSCRIPT_EXT)) {
		// Run Python script. Use the async variant so the per-script venv
		// build (which can take many seconds on first run) does not block
		// the GTK main thread.
		execute_python_script_async(script_file, TRUE, NULL, FALSE, FALSE, get_python_debug_mode(),
				script_file /* venv_identity_path: canonical script file */,
				NULL /* pep723_source: read from file */);
		// script_file ownership transferred to execute_python_script_async
	} else if (g_str_has_suffix(script_file, SCRIPT_EXT)) {
		/* Last thing before running the script, disable widgets except for Stop */
		script_widgets_enable(FALSE);
		// Run regular script
		GFile *file = g_file_new_for_path(script_file);
		GError *error = NULL;
		const GFileInfo *info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_SIZE,
				G_FILE_QUERY_INFO_NONE, NULL, &error);
		if (info) {
			GInputStream *input_stream = (GInputStream*) g_file_read(file, NULL, &error);

			if (input_stream == NULL) {
				if (error != NULL) {
					g_clear_error(&error);
					siril_log_message(_("File [%s] does not exist\n"), script_file);
				}
				g_object_unref(file);
				script_widgets_enable(TRUE);
				g_free(script_file);
				return;
			}
			com.script_thread = g_thread_new("script", execute_script, input_stream);
		}
		g_free(script_file);
		g_object_unref(file);
	} else {
		siril_log_message(_("Unknown script type: %s\n"), script_file);
		script_widgets_enable(TRUE);
		g_free(script_file);
	}
}

/* GAction backer for the "win.script-execute" menu items.  Activation
 * parameter is a GVariant string holding the absolute script path.
 * If the most recent button press on the script popover was secondary
 * (within 500 ms), redirect to the edit action — see the comment on
 * script_secondary_press_time. */
static void script_execute_action_cb(GSimpleAction *action, GVariant *parameter,
                                     gpointer user_data) {
	if (!parameter)
		return;
	if (script_secondary_press_time != 0
	    && g_get_monotonic_time() - script_secondary_press_time < 500000) {
		script_secondary_press_time = 0;
		script_edit_action_cb(action, parameter, user_data);
		return;
	}
	script_secondary_press_time = 0;
	const gchar *path = g_variant_get_string(parameter, NULL);
	if (path && *path)
		script_dispatch(path);
}

/* GtkModelButton (the internal widget GtkPopoverMenu uses for each
 * menu entry) does NOT expose its bound action via the GtkActionable
 * interface — gtk_actionable_get_action_name returns NULL on it, so we
 * can't identify which script a click corresponds to by walking the
 * widget tree.  Worse, model-button activates on release of any mouse
 * button, and returning TRUE from a capture-phase legacy event handler
 * on the popover does not stop that release-time activation.
 *
 * Workaround: record the most recent button press from the script
 * popover and have script_execute_action_cb consult that flag.  If the
 * triggering press was secondary, redirect to script_edit_action_cb;
 * otherwise run the script as normal.  The 500 ms window is large
 * enough to outlast any plausible press-release latency but small
 * enough that a stale flag from an earlier interaction can't redirect
 * an unrelated future invocation. */
static gboolean on_script_popover_event(GtkEventControllerLegacy *ctl,
                                        GdkEvent *event,
                                        gpointer user_data) {
	(void) ctl; (void) user_data;
	if (gdk_event_get_event_type(event) != GDK_BUTTON_PRESS)
		return FALSE;
	guint button = gdk_button_event_get_button(event);
	if (button == GDK_BUTTON_SECONDARY)
		script_secondary_press_time = g_get_monotonic_time();
	else
		script_secondary_press_time = 0;  /* primary cancels any pending flag */
	return FALSE;  /* let the model-button activate; we redirect in the action cb */
}

/* Install the popover-level capture handler once per popover.
 * gtk_menu_button_get_popover may return NULL until the menu model is
 * set and (in some GTK builds) the menu has been opened at least once.
 * We retry from notify::active in case the popover wasn't there yet. */
static void install_capture_on_popover(GtkPopover *popover) {
	if (!popover) return;
	if (g_object_get_data(G_OBJECT(popover), "siril-rmb-installed"))
		return;
	GtkEventController *ctl = gtk_event_controller_legacy_new();
	gtk_event_controller_set_propagation_phase(ctl, GTK_PHASE_CAPTURE);
	g_signal_connect(ctl, "event", G_CALLBACK(on_script_popover_event), NULL);
	gtk_widget_add_controller(GTK_WIDGET(popover), ctl);
	g_object_set_data(G_OBJECT(popover), "siril-rmb-installed",
	                  GINT_TO_POINTER(1));
}

/* When the menu button becomes active, ensure the popover has its
 * capture-phase right-click handler attached.  Cheap on the second
 * and later opens (idempotent via the "siril-rmb-installed" flag). */
static void on_menuscript_active_notify(GObject *obj, GParamSpec *pspec,
                                        gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkMenuButton *mb = GTK_MENU_BUTTON(obj);
	if (!gtk_menu_button_get_active(mb)) return;
	install_capture_on_popover(gtk_menu_button_get_popover(mb));
}

/* Called after every gtk_menu_button_set_menu_model.  Tries to install
 * the popover-level capture handler immediately; if the popover hasn't
 * been built yet, connects to notify::active to install on first open. */
static void install_script_popover_rmb(void) {
	if (!menuscript) return;
	install_capture_on_popover(gtk_menu_button_get_popover(GTK_MENU_BUTTON(menuscript)));
	if (g_object_get_data(G_OBJECT(menuscript), "siril-rmb-active-watch"))
		return;
	g_signal_connect(menuscript, "notify::active",
	                 G_CALLBACK(on_menuscript_active_notify), NULL);
	g_object_set_data(G_OBJECT(menuscript), "siril-rmb-active-watch",
	                  GINT_TO_POINTER(1));
}

void install_script_actions_once(GActionMap *map) {
	if (!map)
		return;
	if (!g_action_map_lookup_action(map, "script-execute")) {
		GSimpleAction *act = g_simple_action_new("script-execute",
		                                         G_VARIANT_TYPE_STRING);
		g_signal_connect(act, "activate", G_CALLBACK(script_execute_action_cb), NULL);
		g_action_map_add_action(map, G_ACTION(act));
		g_object_unref(act);
	}
	if (!g_action_map_lookup_action(map, "script-edit")) {
		GSimpleAction *act = g_simple_action_new("script-edit",
		                                         G_VARIANT_TYPE_STRING);
		g_signal_connect(act, "activate", G_CALLBACK(script_edit_action_cb), NULL);
		g_action_map_add_action(map, G_ACTION(act));
		g_object_unref(act);
	}
}

gboolean test_last_subdir(const gchar *path, const gchar *expected_subdir) {
	g_return_val_if_fail(path != NULL, FALSE);
	g_return_val_if_fail(expected_subdir != NULL, FALSE);

	gchar *dir = g_path_get_dirname(path);
	gchar *last_dir_component = g_path_get_basename(dir);
	gboolean result = (g_strcmp0(last_dir_component, expected_subdir) == 0);

	g_free(last_dir_component);
	g_free(dir);

	return result;
}

static gint compare_basenames(gconstpointer a, gconstpointer b) {
	// Handle NULL inputs
	if (a == NULL && b == NULL) return 0;
	if (a == NULL) return -1;
	if (b == NULL) return 1;

	const gchar *path_a = (const gchar*) a;
	const gchar *path_b = (const gchar*) b;

	gchar *basename_a = g_path_get_basename(path_a);
	gchar *basename_b = g_path_get_basename(path_b);
	gint result;
	if (!path_a || !*path_a) {
		result = -1;
	} else if (!path_b || !*path_b) {
		result = 1;
	} else {
	// Use g_utf8_collate for proper Unicode comparison
	// This handles accented characters correctly
	result = g_utf8_collate(g_utf8_casefold(basename_a, -1),
								g_utf8_casefold(basename_b, -1));
	}
	g_free(basename_a);
	g_free(basename_b);

	return result;
}

/* GAction backer for the "win.script-edit" menu items.  Activation
 * parameter is a GVariant string holding the absolute script path; we
 * read its contents and hand them to new_script() to pop the editor.
 * GTK4 GtkPopoverMenu doesn't expose a per-item right-click hook, so we
 * recreate the GTK3 secondary-button shortcut by attaching a
 * CAPTURE-phase secondary-click gesture to the popover (see
 * install_script_popover_rmb).  The gesture picks the button under the
 * cursor, reads its action target, and dispatches `win.script-edit`. */
static void script_edit_action_cb(GSimpleAction *action, GVariant *parameter,
                                  gpointer user_data) {
	(void)action; (void)user_data;
	if (!parameter) return;
	const gchar *path = g_variant_get_string(parameter, NULL);
	if (!path || !*path) return;
	gchar *contents = NULL;
	gsize length = 0;
	GError *error = NULL;
	if (g_file_get_contents(path, &contents, &length, &error) && length > 0) {
		const char *ext = get_filename_ext(path);
		new_script(contents, length, ext);
		g_free(contents);
	} else {
		gchar *msg = g_strdup_printf(_("Error loading script contents: %s\n"),
		                             error ? error->message : "");
		siril_log_error(msg);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
		g_free(msg);
		if (error) g_error_free(error);
	}
}

/* Phase 15: GMenu-based dirname helper.  Returns or creates the GMenu
 * subsection corresponding to the script's directory (capitalised).
 *
 * `label_suffix` is appended to the submenu's display label.  GtkPopoverMenu
 * uses a single internal GtkStack across the whole popover; each submenu
 * label becomes a stack page name, so labels must be globally unique.
 * The edit-tree (parallel "Edit Scripts" submenu) passes " (edit)" so its
 * directory submenus don't collide with the run-tree's. */
static GMenu *get_py_submenu(const gchar *script_path, GMenu *menu_py,
                             GHashTable *py_submenus, const gchar *label_suffix) {
	gchar *dir_path = g_path_get_dirname(script_path);
	gchar *dir_name = g_path_get_basename(dir_path);
	gchar *capitalized = g_strdup(dir_name);
	if (capitalized && capitalized[0])
		capitalized[0] = g_ascii_toupper(capitalized[0]);
	gchar *display_label = label_suffix ? g_strconcat(capitalized, label_suffix, NULL)
	                                    : g_strdup(capitalized);

	GMenu *submenu = (GMenu *)g_hash_table_lookup(py_submenus, display_label);
	if (!submenu) {
		submenu = g_menu_new();
		g_menu_append_submenu(menu_py, display_label, G_MENU_MODEL(submenu));
		/* hash table takes ownership of the display_label key string and
		 * a borrowed reference to the GMenu (which is owned by menu_py). */
		g_hash_table_insert(py_submenus, g_strdup(display_label), submenu);
	}

	g_free(dir_path);
	g_free(dir_name);
	g_free(capitalized);
	g_free(display_label);
	return submenu;
}

/* Append a script entry to the appropriate run-tree submenu.  Each item
 * targets `win.script-execute` with the script path packed into the
 * action target.  Edit-on-right-click is wired via a CAPTURE-phase
 * gesture installed on the popover (see install_script_popover_rmb). */
static void append_script_entries(const gchar *display_name, const gchar *full_path,
                                  GMenu *menu_ssf, GMenu *menu_py, GHashTable *py_submenus) {
	const gchar *extension = get_filename_ext(display_name);
	if (!extension) return;
	gboolean match_ssf = !g_strcmp0(extension, SCRIPT_EXT);
	gboolean match_py  = !g_strcmp0(extension, PYSCRIPT_EXT)
	                  || !g_strcmp0(extension, PYCSCRIPT_EXT);
	if (!match_ssf && !match_py) return;

	add_to_script_search_list(display_name, full_path);

	GMenuItem *run = g_menu_item_new(display_name, NULL);
	g_menu_item_set_action_and_target_value(run, "win.script-execute",
	                                        g_variant_new_string(full_path));
	if (match_ssf) {
		g_menu_append_item(menu_ssf, run);
	} else {
		GMenu *psub = get_py_submenu(full_path, menu_py, py_submenus, NULL);
		g_menu_append_item(psub, run);
	}
	g_object_unref(run);
}

static int initialize_script_menu(gboolean verbose, gboolean first_run) {
	GSList *list, *script_paths, *s;

	script_menu_init_statics();

	if (script_search_entries) {
		g_slist_free_full(script_search_entries, free_script_entry);
		script_search_entries = NULL;
	}
	script_search_entry_widget = NULL;
	script_search_listbox = NULL;
	script_search_scroll = NULL;

	/* Phase 15: install GAction backers once on the control window. */
	GtkWidget *control_window = lookup_widget("control_window");
	if (control_window && GTK_IS_APPLICATION_WINDOW(control_window))
		install_script_actions_once(G_ACTION_MAP(control_window));

	script_paths = set_list_to_preferences_dialog(com.pref.gui.script_path);

	/* Phase 15: build a GMenu instead of a GtkMenu — each script item
	 * targets the "win.script-execute" action with the full path as the
	 * activation parameter.  The right-click-to-edit shortcut on per
	 * menu items is dropped (no GtkPopoverMenu equivalent). */
	GMenu *menu = g_menu_new();
	GMenu *menu_ssf = g_menu_new();
	GMenu *menu_py  = g_menu_new();

	GHashTable *py_submenus = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

	/* Wrap the submenu rows in their own (label-less) section so they sit at
	 * the same GtkPopoverMenu nesting depth as the "Get Scripts" tail
	 * section below.  A submenu appended directly to the top-level menu and
	 * an item appended inside a <section> render at different left insets on
	 * themes that indent sections, so the submenu labels and the tail item
	 * labels wouldn't line up.  Keeping both groups in sections makes every
	 * row align. */
	GMenu *head = g_menu_new();
	g_menu_append_submenu(head, _("Siril Script Files"), G_MENU_MODEL(menu_ssf));
	g_menu_append_submenu(head, _("Python Scripts"),     G_MENU_MODEL(menu_py));
	g_menu_append_section(menu, NULL, G_MENU_MODEL(head));
	g_object_unref(head);

	GMenu *tail = g_menu_new();
	g_menu_append(tail, _("Get Scripts"),               "win.script-getscripts");
	g_menu_append(tail, _("Script Editor…"),           "win.script-pythonpad");
	g_menu_append(tail, _("Enable Python debug mode"), "win.script-pythondebug");
	g_menu_append_section(menu, NULL, G_MENU_MODEL(tail));
	g_object_unref(tail);

	GMenu *search_section = g_menu_new();
	GMenuItem *search_placeholder = g_menu_item_new(NULL, NULL);
	g_menu_item_set_attribute(search_placeholder, "custom", "s", "scripts-search");
	g_menu_append_item(search_section, search_placeholder);
	g_object_unref(search_placeholder);
	g_menu_append_section(menu, NULL, G_MENU_MODEL(search_section));
	g_object_unref(search_section);

	for (s = script_paths; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			if (verbose)
				siril_log_info(_("Searching for scripts in: \"%s\"...\n"), s->data);

			for (GSList *l = list; l; l = l->next) {
				if (l->data == NULL)
					continue;
				gchar *display_name = l->data;
				gchar *full_path = g_build_filename(s->data, l->data, NULL);
				append_script_entries(display_name, full_path,
				                      menu_ssf, menu_py, py_submenus);
				if (verbose)
					siril_log_message(_("Loading script: %s\n"), l->data);
				g_free(full_path);
			}
			g_slist_free_full(list, g_free);
		}
	}

	// Add scripts from the selections made in preferences
	if (com.pref.use_scripts_repository && com.pref.selected_scripts) {
		// Remove NULL entries from selected_scripts in-place
		// Sort selected_scripts in-place
		com.pref.selected_scripts = g_slist_sort(com.pref.selected_scripts, compare_basenames);
		// Iterate and prune any items not found in repo (if purge_removed)
		GSList *l = com.pref.selected_scripts;
		l = com.pref.selected_scripts;
		while (l != NULL) {
			GSList *next = l->next;
			gchar *path = l->data;
			// Remove any scripts with a NULL path
			if (!path) {
				com.pref.selected_scripts = g_slist_delete_link(com.pref.selected_scripts, l);
				l = next;
				continue;
			}
			gboolean exists = g_file_test(path, G_FILE_TEST_EXISTS);
			gboolean included = !com.repo_scripts;

			if (com.repo_scripts != NULL && exists) {
				for (GSList *it = com.repo_scripts; it; it = it->next) {
					if (g_strrstr(path, it->data)) {
						included = TRUE;
						break;
					}
				}
			}

			if (!first_run && (!exists && included)) {
				siril_log_warning(_("Script %s no longer exists in repository, removing from Scripts menu...\n"), path);
				// Remove the list element and free it as well as its data
				g_free(path);
				com.pref.selected_scripts = g_slist_delete_link(com.pref.selected_scripts, l);
			} else if (included) {
				gchar *basename = g_path_get_basename(path);
				append_script_entries(basename, path, menu_ssf, menu_py, py_submenus);
				if (verbose)
					siril_log_message(_("Loading script from repository: %s\n"), basename);
				g_free(basename);
			}

			l = next;
		}
	}

	// Add core scripts if they're not already in the menu
	for (GSList *core_iter = com.repo_scripts; core_iter; core_iter = core_iter->next) {
		const gchar *script_path = (gchar*)core_iter->data;
		if (test_last_subdir(script_path, "core")) {
			gboolean already_added = FALSE;
			for (GSList *selected = com.pref.selected_scripts; selected; selected = selected->next) {
				if (!selected->data) continue;
				if (g_strrstr((gchar*)selected->data, script_path)) {
					already_added = TRUE;
					break;
				}
			}

			if (!already_added) {
				gchar *basename = g_path_get_basename(script_path);
				const char *extension = get_filename_ext(basename);
				if (!extension) {
					g_free(basename);
					continue;
				}
				gchar *full_path = g_build_filename(siril_get_scripts_repo_path(), script_path, NULL);
				append_script_entries(basename, full_path, menu_ssf, menu_py, py_submenus);
				if (verbose)
					siril_log_message(_("Adding core script to menu: %s\n"), basename);
				g_free(basename);
				g_free(full_path);
			}
		}
	}

	// Phase 15: install the populated GMenuModel on the menu button.
	gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(menuscript), G_MENU_MODEL(menu));

	GtkPopover *search_popover = gtk_menu_button_get_popover(GTK_MENU_BUTTON(menuscript));
	if (GTK_IS_POPOVER_MENU(search_popover)) {
		GtkWidget *search_widget = create_script_search_widget();
		gtk_popover_menu_add_child(GTK_POPOVER_MENU(search_popover),
		                           search_widget, "scripts-search");
	}

	install_script_popover_rmb();

	g_object_unref(menu_ssf);
	g_object_unref(menu_py);
	g_object_unref(menu);

	g_hash_table_destroy(py_submenus);

	return 0;
}

// Called when the specified scripts directories are initialized. Just a wrapper so that the function
// has the right signature to be called in an idle in the GTK thread.
// You must call this from a secondary thread and call gui_mutex_lock() / unlock() around it
gboolean initialize_script_menu_idle(gpointer data) {
	gboolean state = (gboolean) GPOINTER_TO_INT(data);
	initialize_script_menu(state, TRUE);
	return FALSE;
}

// This function updates the scripts menu, first removing the old one, it is called at startup
// after refreshing the repository
// You must call this from a secondary thread and call gui_mutex_lock() / unlock() around it
gboolean refresh_script_menu_idle(gpointer user_data) {
	gboolean verbose = (gboolean) GPOINTER_TO_INT(user_data);
	if (menuscript) {
		// Remove the popup while we refresh the menu
		gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(menuscript), NULL);
	}
	initialize_script_menu(verbose, FALSE);
#ifdef HAVE_LIBGIT2
	fill_script_repo_tree(FALSE);
#endif
	return FALSE;
}

// This is called from preferences or the reloadscripts command to refresh the
// script menu.
gpointer refresh_scripts_in_thread(gpointer user_data) {
	GSList *list = get_list_from_preferences_dialog();
	// TODO: is there anything to stop refreshscripts being called from a script run by siril-cli?
	// if not, we need to prevent it as get_list_from_preferences_dialog() uses GTK code and will fail,
	// probably badly.

	if (list == NULL) {
		gchar *err = siril_log_error(_("Cannot refresh the scripts if the list is empty.\n"));
		queue_warning_message_dialog(_("Warning"), err);
		g_free(err);
	} else {
		g_slist_free_full(com.pref.gui.script_path, g_free);
		com.pref.gui.script_path = list;
		gui_mutex_lock();
		execute_idle_and_wait_for_it(initialize_script_menu_idle, GINT_TO_POINTER(1));
		gui_mutex_unlock();
	}
	return FALSE;
}

// Called from preferences to refresh the script menu
gpointer refresh_script_menu_in_thread(gpointer user_data) {
	gui_mutex_lock();
	execute_idle_and_wait_for_it(refresh_script_menu_idle, user_data);
	gui_mutex_unlock();
	return GINT_TO_POINTER(0);
}

GSList *get_list_from_preferences_dialog() {
	GSList *list = NULL;
	script_menu_init_statics();
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(script_path_textview);
	GtkTextIter start, end;
	gchar *txt;
	gint i = 0;
	gtk_text_buffer_get_bounds(tbuf, &start, &end);
	txt = gtk_text_buffer_get_text(tbuf, &start, &end, TRUE);
	if (txt) {
		gchar **token = g_strsplit(txt, "\n", -1);
		while (token[i]) {
			if (*token[i] != '\0')
				list = g_slist_prepend(list, g_strdup(token[i]));
			i++;
		}
		list = g_slist_reverse(list);
		g_strfreev(token);
	}

	return list;
}

GSList *set_list_to_preferences_dialog(GSList *list) {
	clear_gtk_list();
	if (list == NULL) {
		list = initialize_script_paths();
	}
	for (GSList *l = list; l; l = l->next) {
		add_path_to_gtkText((gchar *) l->data);
	}
	return list;
}

void on_get_scripts_clicked(gpointer user_data) {
	siril_open_dialog("settings_window");
	script_menu_init_statics();
	gtk_stack_set_visible_child(script_stack_pref, script_scripts_page);
}
