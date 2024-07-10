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
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"
#include "script_menu.h"

#define CONFIRM_RUN_SCRIPTS _("You are about to use scripts. Running automatic scripts is something that is easy and generally it provides a nice image. However you have to keep in mind that scripts are not magic; automatic choices are made where human decision would probably be better. Also, every command used in a script is available on the interface with a better parameter control.")

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
	GtkWidget *notebook_center_box = lookup_widget("notebook_center_box");
	GtkWidget *command = lookup_widget("command");
	GtkWidget *notebook1 = lookup_widget("notebook1");
	GtkWidget *headerbar = lookup_widget("headerbar");
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gtk_widget_set_sensitive(notebook_center_box, status);
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
		if (g_str_has_suffix(file, SCRIPT_EXT)) {
			gchar *str = (gchar*) remove_ext_from_filename(file);

			list = g_slist_prepend(list, str);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	g_dir_close(dir);

	return list;
}

static void on_script_execution(GtkMenuItem *menuitem, gpointer user_data) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (com.pref.gui.warn_script_run) {
		gboolean dont_show_again;
		gboolean confirm = siril_confirm_dialog_and_remember(
				_("Please read me before using scripts"), CONFIRM_RUN_SCRIPTS, _("Run Script"), &dont_show_again);
		if (!confirm)
			return;

		com.pref.gui.warn_script_run = !dont_show_again;
		writeinitfile();
	}

	if (com.script_thread)
		g_thread_join(com.script_thread);

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	gchar *script_file;
	if (g_str_has_suffix((gchar *) user_data, SCRIPT_EXT)) {
		script_file= g_strdup_printf("%s", (gchar *) user_data); // remote scripts

	} else {
		script_file= g_strdup_printf("%s%s", (gchar *) user_data, SCRIPT_EXT); // local scripts
	}

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

			g_free(script_file);
			g_object_unref(file);
			return;
		}
		/* Last thing before running the script, disable widgets except for Stop */
		script_widgets_enable(FALSE);

		/* Run the script */
		siril_log_message(_("Starting script %s\n"), script_file);
		com.script_thread = g_thread_new("script", execute_script, input_stream);
	}

	g_free(script_file);
	g_object_unref(file);
}

int initialize_script_menu(gboolean verbose) {
	GSList *list, *script_paths, *s;
#ifdef HAVE_LIBGIT2
	GList *ss;
#endif
	gint nb_item = 0;

	if (!menuscript)
		menuscript = lookup_widget("header_scripts_button");

	script_paths = set_list_to_preferences_dialog(com.pref.gui.script_path);

	GtkWidget *menu = gtk_menu_new();

	for (s = script_paths; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			/* write separator but not for the first one */
			if (nb_item != 0) {
				GtkWidget *separator = gtk_separator_menu_item_new();
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), separator);
				gtk_widget_show(separator);
			}
			if (verbose)
				siril_log_color_message(_("Searching scripts in: \"%s\"...\n"), "green", s->data);

			for (GSList *l = list; l; l = l->next) {
				nb_item++;
				/* write an item per script file */
				GtkWidget *menu_item;

				menu_item = gtk_menu_item_new_with_label(l->data);
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
				gchar *full_path = g_build_filename(s->data, l->data,
						NULL);
				g_signal_connect(G_OBJECT(menu_item), "activate",
						G_CALLBACK(on_script_execution), (gchar * ) full_path);
				if (verbose)
					siril_log_message(_("Loading script: %s\n"), l->data);
				gtk_widget_show(menu_item);
			}
			g_slist_free_full(list, g_free);
		}
	}
#ifdef HAVE_LIBGIT2
	if (com.pref.use_scripts_repository && g_list_length(com.pref.selected_scripts) > 0 ) {
		GtkWidget *separator = gtk_separator_menu_item_new();
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), separator);
		gtk_widget_show(separator);
		GList *new_list = NULL;
		for (ss = com.pref.selected_scripts ; ss ; ss = ss->next) {
			nb_item++;
			/* check that the selected script is still in the repository */

			gboolean included = FALSE;
			GList *iterator;
			for (iterator = gui.repo_scripts ; iterator ; iterator = iterator->next) {
				if (g_strrstr((gchar*) ss->data, (gchar*)iterator->data)) {
					included = TRUE;
				}
			}
			if (included) {
				/* write an item per script file */
				GtkWidget *menu_item;
				gchar* basename = g_path_get_basename(ss->data);
				char* basename_no_ext = remove_ext_from_filename(basename);
				g_free(basename);
				menu_item = gtk_menu_item_new_with_label(basename_no_ext);
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
				g_signal_connect(G_OBJECT(menu_item), "activate",
						G_CALLBACK(on_script_execution), (gchar * ) ss->data);
				if (verbose)
					siril_log_message(_("Loading script: %s\n"), basename_no_ext);
				free(basename_no_ext);
				gtk_widget_show(menu_item);
				new_list = g_list_prepend(new_list, ss->data);
			} else {
				// Print this whether verbose or not, as it is likely to be unexpected
				siril_log_color_message(_("Script %s no longer exists in repository, removing from Scripts menu...\n"), "salmon", ss->data);
			}
		}
		GList *tmp = com.pref.selected_scripts;
		com.pref.selected_scripts = new_list;
		g_list_free (g_steal_pointer (&tmp));
	}
#endif
	if (!nb_item) {
		gtk_widget_hide(menuscript);
		return 0;
	}
	gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), menu);
	if (!gtk_widget_get_visible(menuscript))
		gtk_widget_show(menuscript);
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
