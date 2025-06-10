/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/python_gui.h"
#include "algos/sorting.h"
#include "script_menu.h"

#define CONFIRM_RUN_SCRIPTS _("You are about to use scripts. Note that scripts execute code with your current user privileges. While Siril Script Files can only execute Siril commands and a very small number of specific external programs, Python scripts are considerably more powerful and execute code not written by the Siril team. Ensure you obtain scripts from a reputable source.")

static GtkWidget *menuscript = NULL;

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
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter iter;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
		tbuf = gtk_text_view_get_buffer(text);
	}

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	gtk_text_buffer_insert(tbuf, &iter, path, strlen(path));
	gtk_text_buffer_insert(tbuf, &iter, "\n", strlen("\n"));

	/* scroll to end */
	gtk_text_buffer_get_end_iter(tbuf, &iter);
	GtkTextMark *insert_mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	gtk_text_view_scroll_to_mark(text, insert_mark, 0.0, TRUE, 0.0, 1.0);
}

static void clear_gtk_list() {
	GtkTextView *text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
}

void script_widgets_enable(gboolean status) {
	GtkWidget *tab1 = lookup_widget("conversion_tab");
	GtkWidget *tab2 = lookup_widget("sequence_tab");
	GtkWidget *tab3 = lookup_widget("calibration_tab");
	GtkWidget *tab4 = lookup_widget("registration_tab");
	GtkWidget *tab5 = lookup_widget("plot_tab");
	GtkWidget *tab6 = lookup_widget("stacking_tab");
	GtkWidget *command = lookup_widget("command");
	GtkWidget *notebook1 = lookup_widget("notebook1");
	GtkWidget *headerbar = lookup_widget("headerbar");
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gtk_widget_set_sensitive(tab1, status);
	gtk_widget_set_sensitive(tab2, status);
	gtk_widget_set_sensitive(tab3, status);
	gtk_widget_set_sensitive(tab4, status);
	gtk_widget_set_sensitive(tab5, status);
	gtk_widget_set_sensitive(tab6, status);
	gtk_widget_set_sensitive(command, status);
	gtk_widget_set_sensitive(notebook1, status);
	gtk_widget_set_sensitive(headerbar, status);
	gtk_widget_set_sensitive(toolbarbox, status);
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

static void on_script_execution(GtkMenuItem *menuitem, gpointer user_data) {
	if (get_thread_run()) {
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
		// Run Python script
		execute_python_script(script_file, TRUE, FALSE, NULL, FALSE, FALSE, get_python_debug_mode());
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
				return;
			}
			com.script_thread = g_thread_new("script", execute_script, input_stream);
		}
		g_object_unref(file);
	} else {
		siril_log_message(_("Unknown script type: %s\n"), script_file);
		script_widgets_enable(TRUE);
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
	if (!path_a || !*path_a) return -1;
	if (!path_b || !*path_b) return 1;
	// Use g_utf8_collate for proper Unicode comparison
	// This handles accented characters correctly
	gint result = g_utf8_collate(g_utf8_casefold(basename_a, -1),
								g_utf8_casefold(basename_b, -1));

	g_free(basename_a);
	g_free(basename_b);

	return result;
}

// Helper function to get or create Python submenu based on script path
static GtkWidget* get_py_submenu(const gchar *script_path, GtkWidget *menu_py, GHashTable *py_submenus) {
	gchar *dir_path = g_path_get_dirname(script_path);
	gchar *dir_name = g_path_get_basename(dir_path);

	// Capitalize the directory name
	gchar *capitalized = g_strdup(dir_name);
	if (capitalized && capitalized[0]) {
		capitalized[0] = g_ascii_toupper(capitalized[0]);
	}

	GtkWidget *submenu = (GtkWidget*)g_hash_table_lookup(py_submenus, capitalized);
	if (!submenu) {
		submenu = gtk_menu_new();
		GtkWidget *submenu_item = gtk_menu_item_new_with_label(capitalized);
		gtk_menu_item_set_submenu(GTK_MENU_ITEM(submenu_item), submenu);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu_py), submenu_item);
		gtk_widget_show(submenu_item);

		g_hash_table_insert(py_submenus, g_strdup(capitalized), submenu);
	}

	g_free(dir_path);
	g_free(dir_name);
	g_free(capitalized);
	return submenu;
}

static int initialize_script_menu(gboolean verbose, gboolean first_run) {
	GSList *list, *script_paths, *s;

	if (!menuscript)
		menuscript = lookup_widget("header_scripts_button");

	script_paths = set_list_to_preferences_dialog(com.pref.gui.script_path);

	GtkWidget *menu = gtk_menu_new();
	GtkWidget *menu_ssf = gtk_menu_new();
	GtkWidget *menu_py = gtk_menu_new();

	// Hash table to store Python script submenus by directory name
	GHashTable *py_submenus = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

	GtkWidget *menu_item_ssf = gtk_menu_item_new_with_label(_("Siril Script Files"));
	GtkWidget *menu_item_py = gtk_menu_item_new_with_label(_("Python Scripts"));
	gtk_widget_set_tooltip_markup(menu_item_py,
			_("<b>EXPERIMENTAL</b>: python scripts are currently an experimental feature. "
			"Please read the documentation for details..."));
	GtkWidget *sep = gtk_separator_menu_item_new();
	GtkWidget *menu_item_pythonpad = gtk_menu_item_new_with_label(_("Script Editor..."));
	g_signal_connect(G_OBJECT(menu_item_pythonpad), "activate", G_CALLBACK(on_open_pythonpad), NULL);
	GtkWidget *menu_item_pythondebug = gtk_check_menu_item_new_with_label(_("Enable Python debug mode"));
	GObject *existing = gtk_builder_get_object(gui.builder, "pythondebugtoggle");
	if (!existing) {
		gtk_builder_expose_object(gui.builder, "pythondebugtoggle", G_OBJECT(menu_item_pythondebug));
	}
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item_pythondebug), FALSE);
	g_signal_connect(G_OBJECT(menu_item_pythondebug), "toggled", G_CALLBACK(on_pythondebug_toggled), NULL);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_ssf);
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item_ssf), menu_ssf);
	gtk_widget_show(menu_item_ssf);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_py);
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item_py), menu_py);
	gtk_widget_show(menu_item_py);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), sep);
	gtk_widget_show(sep);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_pythonpad);
	gtk_widget_show(menu_item_pythonpad);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_pythondebug);
	gtk_widget_show(menu_item_pythondebug);

	gchar *previous_directory_ssf = NULL;
	gchar *previous_directory_py = NULL;
	gboolean first_item_ssf = TRUE;

	for (s = script_paths; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			if (verbose)
				siril_log_color_message(_("Searching for scripts in: \"%s\"...\n"), "green", s->data);

			for (GSList *l = list; l; l = l->next) {
				if (l->data == NULL)
					continue;
				GtkWidget *menu_item;
				gchar *display_name = l->data;
				const gchar *extension = get_filename_ext(display_name);
				gchar *full_path = g_build_filename(s->data, l->data, NULL);
				gchar *current_directory = g_path_get_dirname(full_path);
				if (!current_directory || *current_directory == '\0') {
					siril_debug_print("Directory name not found for script %s, skipping...\n", display_name);
					g_free(current_directory);
					continue;
				}
				if (extension) {
					gboolean match_ssf = !g_strcmp0(extension, SCRIPT_EXT);
					gboolean match_py = !g_strcmp0(extension, PYSCRIPT_EXT) || !g_strcmp0(extension, PYCSCRIPT_EXT);
					if (match_ssf) {
						if (!first_item_ssf && (!previous_directory_ssf || g_strcmp0(current_directory, previous_directory_ssf) != 0)) {
							GtkWidget *separator = gtk_separator_menu_item_new();
							gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), separator);
							gtk_widget_show(separator);
						}
						first_item_ssf = FALSE;
						g_free(previous_directory_ssf);
						previous_directory_ssf = g_strdup(current_directory);
					} else {
						if ( match_py) {
							g_free(previous_directory_py);
							previous_directory_py = g_strdup(current_directory);
						}
					}

					menu_item = gtk_menu_item_new_with_label(display_name);

					if (match_ssf) {
						gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
					} else if (match_py) {
						GtkWidget *py_submenu = get_py_submenu(full_path, menu_py, py_submenus);
						gtk_menu_shell_append(GTK_MENU_SHELL(py_submenu), menu_item);
					}

					g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(on_script_execution), full_path);
					if (verbose)
						siril_log_message(_("Loading script: %s\n"), l->data);

					gtk_widget_show(menu_item);
				}
				g_free(current_directory);
			}
			g_slist_free_full(list, g_free);
		}
	}
	g_free(previous_directory_ssf);
	g_free(previous_directory_py);

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
			gboolean included = !gui.repo_scripts;

			if (gui.repo_scripts != NULL && exists) {
				for (GSList *it = gui.repo_scripts; it; it = it->next) {
					if (g_strrstr(path, it->data)) {
						included = TRUE;
						break;
					}
				}
			}

			if (!first_run && (!exists && included)) {
				siril_log_color_message(_("Script %s no longer exists in repository, removing from Scripts menu...\n"), "salmon", path);
				// Remove the list element and free it as well as its data
				g_free(path);
				com.pref.selected_scripts = g_slist_delete_link(com.pref.selected_scripts, l);
			} else if (included) {
				// Build menu entry
				GtkWidget *menu_item = NULL;
				gchar *basename = g_path_get_basename(path);
				const char *extension = get_filename_ext(basename);

				menu_item = gtk_menu_item_new_with_label(basename);

				if (extension && g_strcmp0(extension, SCRIPT_EXT) == 0) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
				} else if (extension && ((g_strcmp0(extension, PYSCRIPT_EXT) == 0) || (g_strcmp0(extension, PYCSCRIPT_EXT) == 0))) {
					GtkWidget *py_submenu = get_py_submenu(path, menu_py, py_submenus);
					gtk_menu_shell_append(GTK_MENU_SHELL(py_submenu), menu_item);
				}

				g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(on_script_execution), g_strdup(path));

				if (verbose) {
					siril_log_message(_("Loading script from repository: %s\n"), basename);
				}

				gtk_widget_show(menu_item);
				g_free(basename);
			}

			l = next;
		}
	}

	// Add core scripts if they're not already in the menu
	for (GSList *core_iter = gui.repo_scripts; core_iter; core_iter = core_iter->next) {
		const gchar *script_path = (gchar*)core_iter->data;
		if (test_last_subdir(script_path, "core")) {
			// Check if this core script is already in selected_scripts
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
				GtkWidget *menu_item = gtk_menu_item_new_with_label(basename);
				gchar *full_path = g_build_filename(siril_get_scripts_repo_path(), script_path, NULL);

				if (extension && g_strcmp0(extension, SCRIPT_EXT) == 0) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
				} else if (extension && ((g_strcmp0(extension, PYSCRIPT_EXT) == 0) ||
								(g_strcmp0(extension, PYCSCRIPT_EXT) == 0))) {
					GtkWidget *py_submenu = get_py_submenu(full_path, menu_py, py_submenus);
					gtk_menu_shell_append(GTK_MENU_SHELL(py_submenu), menu_item);
				} else {
					g_free(basename);
					g_free(full_path);
					continue;
				}

				g_signal_connect(G_OBJECT(menu_item), "activate",
								 G_CALLBACK(on_script_execution), full_path);

				if (verbose)
					siril_log_message(_("Adding core script to menu: %s\n"), basename);
				g_free(basename);
				gtk_widget_show(menu_item);

			}
		}
	}

	// Now we have finished populating it, set the menu_button popup
	gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), menu);

	// Clean up hash table
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
		gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), NULL);
	}
	initialize_script_menu(verbose, FALSE);
	fill_script_repo_tree(FALSE);
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
		gchar *err = siril_log_color_message(_("Cannot refresh the scripts if the list is empty.\n"), "red");
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
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter start, end;
	gchar *txt;
	gint i = 0;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
		tbuf = gtk_text_view_get_buffer(text);
	}
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
