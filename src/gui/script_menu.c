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
/* Constant available since Shell32.dll 4.72 */
#ifndef CSIDL_APPDATA
#define CSIDL_APPDATA 0x001a
#endif
#endif
#include <string.h>
#include <locale.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command.h" // for process_close()
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "io/siril_pythonmodule.h"
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
	list = g_slist_prepend(list, g_build_filename(get_special_folder(CSIDL_APPDATA), "siril",
					"scripts", NULL));

	gchar *execpath = g_win32_get_package_installation_directory_of_module(NULL);

	list = g_slist_prepend(list, g_build_filename(execpath, "scripts", NULL));
	g_free(execpath);
#else
	list = g_slist_prepend(list, g_build_filename(siril_get_system_data_dir(), "scripts", NULL));
	if (g_getenv("XDG_CONFIG_HOME") != NULL) {
		list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), getenv("XDG_CONFIG_HOME"), "scripts", NULL));
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
//	GtkWidget *notebook_center_box = lookup_widget("notebook_center_box");
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
//	gtk_widget_set_sensitive(notebook_center_box, status);
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

	if (com.script_thread)
		g_thread_join(com.script_thread);

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	gchar *script_file = g_strdup((gchar *)user_data);

	/* Run the script */
	siril_log_message(_("Starting script %s\n"), script_file);

	if (g_str_has_suffix(script_file, PYSCRIPT_EXT) || g_str_has_suffix(script_file, PYCSCRIPT_EXT)) {
		// Run Python script
		execute_python_script(script_file, TRUE, FALSE, NULL);
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

int initialize_script_menu(gboolean verbose) {
	GSList *list, *script_paths, *s;
#ifdef HAVE_LIBGIT2
	GList *ss;
#endif

	if (!menuscript)
		menuscript = lookup_widget("header_scripts_button");

	script_paths = set_list_to_preferences_dialog(com.pref.gui.script_path);

	GtkWidget *menu = gtk_menu_new();
	GtkWidget *menu_ssf = gtk_menu_new();
	GtkWidget *menu_py = gtk_menu_new();

	GtkWidget *menu_item_ssf = gtk_menu_item_new_with_label(_("Siril Script Files"));
	GtkWidget *menu_item_py = gtk_menu_item_new_with_label(_("Python Scripts"));
	gtk_widget_set_tooltip_markup(menu_item_py,
			"<b>EXPERIMENTAL</b>: python scripts are currently an experimental feature. "
			"Please read the documentation for details...");
	GtkWidget *sep = gtk_separator_menu_item_new();
	GtkWidget *menu_item_pythonpad = gtk_menu_item_new_with_label(_("Script Editor..."));
	g_signal_connect(G_OBJECT(menu_item_pythonpad), "activate", G_CALLBACK(on_open_pythonpad), NULL);

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

	gchar *previous_directory_ssf = NULL;
	gchar *previous_directory_py = NULL;
	gboolean first_item_ssf = TRUE;
	gboolean first_item_py = TRUE;

	for (s = script_paths; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			if (verbose)
				siril_log_color_message(_("Searching for scripts in: \"%s\"...\n"), "green", s->data);

			for (GSList *l = list; l; l = l->next) {
				GtkWidget *menu_item;

				gchar *display_name = g_strdup(l->data);
				const gchar *extension = get_filename_ext(display_name);
				gchar *current_directory = g_path_get_dirname(s->data);

				if (extension) {
					gboolean match_ssf = !g_strcmp0(extension, SCRIPT_EXT);
					gboolean match_py = !g_strcmp0(extension, PYSCRIPT_EXT);
					gboolean match_pyc = !g_strcmp0(extension, PYCSCRIPT_EXT);
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
						if ( match_py || match_pyc) {
							if (!first_item_py && (!previous_directory_py || g_strcmp0(current_directory, previous_directory_py) != 0)) {
								GtkWidget *separator = gtk_separator_menu_item_new();
								gtk_menu_shell_append(GTK_MENU_SHELL(menu_py), separator);
								gtk_widget_show(separator);
							}
							first_item_py = FALSE;
							g_free(previous_directory_py);
							previous_directory_py = g_strdup(current_directory);
						}
					}

					menu_item = gtk_menu_item_new_with_label(display_name);
					gchar *full_path = g_build_filename(s->data, l->data, NULL);

					if (match_ssf) {
						gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
					} else if (match_py || match_pyc) {
						gtk_menu_shell_append(GTK_MENU_SHELL(menu_py), menu_item);
					}

					g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(on_script_execution), full_path);
					if (verbose)
						siril_log_message(_("Loading script: %s\n"), l->data);

					gtk_widget_show(menu_item);
				}
				g_free(current_directory);
				g_free(display_name);
			}
			g_slist_free_full(list, g_free);
		}
	}
	g_free(previous_directory_ssf);
	g_free(previous_directory_py);

	#ifdef HAVE_LIBGIT2
	// Wait for git repository update to complete during startup
	if (com.update_scripts_thread) {
		g_thread_join(com.update_scripts_thread);
		com.update_scripts_thread = NULL;
	}
	// Add scripts from the selections made in preferences
	if (com.pref.use_scripts_repository && g_list_length(com.pref.selected_scripts) > 0) {
		GList *new_list = NULL;
		for (ss = com.pref.selected_scripts; ss; ss = ss->next) {
			gchar *full_path = g_strdup(ss->data);
			if (!g_file_test(full_path, G_FILE_TEST_EXISTS)) {
				siril_log_color_message(_("Script %s no longer exists in repository, removing from Scripts menu...\n"), "salmon", ss->data);
				g_free(full_path);
				continue;
			}

			gboolean included = FALSE;
			GList *iterator;
			for (iterator = gui.repo_scripts; iterator; iterator = iterator->next) {
				if (g_strrstr((gchar*) ss->data, (gchar*) iterator->data)) {
					included = TRUE;
					break;
				}
			}
			if (included) {
				GtkWidget *menu_item;
				gchar *basename = g_path_get_basename(ss->data);
				const char *extension = get_filename_ext(basename);

				menu_item = gtk_menu_item_new_with_label(basename);

				if (extension && g_strcmp0(extension, SCRIPT_EXT) == 0) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
				} else if (extension && ((g_strcmp0(extension, PYSCRIPT_EXT) == 0) || (g_strcmp0(extension, PYCSCRIPT_EXT) == 0))) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_py), menu_item);
				}

				g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(on_script_execution), full_path);
				if (verbose)
					siril_log_message(_("Loading script from repository: %s\n"), basename);
				gtk_widget_show(menu_item);
				new_list = g_list_prepend(new_list, g_strdup(ss->data));

				g_free(basename);
			} else {
				siril_log_color_message(_("Script %s no longer exists in repository, removing from Scripts menu...\n"), "salmon", ss->data);
				g_free(full_path);
			}
		}
		GList *tmp = com.pref.selected_scripts;
		com.pref.selected_scripts = new_list;
		g_list_free_full(tmp, g_free);
	}

	// Add core scripts if they're not already in the menu
	for (GList *core_iter = gui.repo_scripts; core_iter; core_iter = core_iter->next) {
		const gchar *script_path = (gchar*)core_iter->data;
		if (test_last_subdir(script_path, "core")) {
			// Check if this core script is already in selected_scripts
			gboolean already_added = FALSE;
			for (GList *selected = com.pref.selected_scripts; selected; selected = selected->next) {
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
				if (extension && g_strcmp0(extension, SCRIPT_EXT) == 0) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_ssf), menu_item);
				} else if (extension && ((g_strcmp0(extension, PYSCRIPT_EXT) == 0) ||
								(g_strcmp0(extension, PYCSCRIPT_EXT) == 0))) {
					gtk_menu_shell_append(GTK_MENU_SHELL(menu_py), menu_item);
				} else {
					g_free(basename);
					continue;
				}

				gchar *full_path = g_build_filename(siril_get_scripts_repo_path(), script_path, NULL);
				g_signal_connect(G_OBJECT(menu_item), "activate",
								 G_CALLBACK(on_script_execution), full_path);

				if (verbose)
					siril_log_message(_("Adding core script to menu: %s\n"), basename);
				g_free(basename);
				if(!menu_item)
					continue;
				gtk_widget_show(menu_item);

			}
		}
	}

	#endif
	gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), menu);
	return 0;
}

int refresh_script_menu(gboolean verbose) {
	if (menuscript) {
		gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), NULL);
	}
	initialize_script_menu(verbose);
	return 0;
}

int refresh_scripts(gboolean update_list, gchar **error) {
	gchar *err = NULL;
	int retval;
	GSList *list = get_list_from_preferences_dialog();
	if (list == NULL) {
		err = siril_log_color_message(_("Cannot refresh the scripts if the list is empty.\n"), "red");
		retval = 1;
	} else {
		g_slist_free_full(com.pref.gui.script_path, g_free);
		com.pref.gui.script_path = list;
		retval = initialize_script_menu(TRUE);
	}
	if (error) {
		*error = err;
	}
	return retval;
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
