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
 *
 *
 * FITS sequences are not a sequence of FITS files but a FITS file containing a
 * sequence. It simply has as many elements in the third dimension as the
 * number of images in the sequence multiplied by the number of channels per
 * image. Given its use of the third dimension, it's sometimes called FITS cube.
 */
#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/photometric_cc.h" // for reset_spcc_filters() (this is not a GTK function)
#include "gui/preferences.h"
#include "gui/progress_and_log.h"
#include "gui/script_menu.h"
#include "gui/utils.h"
#include "io/siril_git.h"

#ifdef HAVE_LIBGIT2

static GtkListStore *list_store = NULL;

static const char *bg_color[] = {"WhiteSmoke", "#1B1B1B"};

enum {
	COLUMN_CATEGORY = 0, // string
	COLUMN_SCRIPTNAME,   // string
	COLUMN_SELECTED,     // gboolean
	COLUMN_SCRIPTPATH,   // full path to populate into the scripts menu
	COLUMN_BGCOLOR,      // background color
	COLUMN_TYPE,         // string, type of script
	N_COLUMNS
};

static int reset_scripts_repository() {
	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();
	int retval = reset_repository(local_path);
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"),
			_("Success! The local repository is up-to-date with the remote."));
	}
	return retval;
}

static int reset_spcc_repository() {
	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_spcc_repo_path();
	int retval = reset_repository(local_path);
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"),
				_("Success! The local repository is up-to-date with the remote."));
	}
	return retval;
}

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(
			gtk_builder_get_object(gui.builder, "liststore_script_repo"));
	}
}

static gboolean fill_script_repo_list_idle(gpointer p) {
	GtkTreeView *tview = (GtkTreeView *)p;
	GtkTreeIter iter;
	if (!tview)
		return FALSE;
	if (list_store)
		gtk_list_store_clear(list_store);
	get_list_store();
	gint sort_column_id;
	GtkSortType order;
	// store sorted state of list_store, disable sorting, disconnect from the
	// view, fill, reconnect and re-apply sort
	gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(list_store), &sort_column_id, &order);
	gtk_tree_sortable_set_sort_column_id(
		GTK_TREE_SORTABLE(list_store), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID,
		GTK_SORT_ASCENDING);
	gtk_tree_view_set_model(tview, NULL);
	if (gui.repo_scripts) {
		int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
		GList *iterator;
		for (iterator = gui.repo_scripts; iterator; iterator = iterator->next) {
			// here we populate the GtkTreeView from GList gui.repo_scripts
			const gchar *category;
			gboolean included = FALSE;
			gboolean core = FALSE;
			if (test_last_subdir((gchar *)iterator->data, "preprocessing")) {
				category = _("Preprocessing");
			} else if (test_last_subdir((gchar *)iterator->data, "processing")) {
				category = _("Processing");
			} else if (test_last_subdir((gchar *)iterator->data, "utility")) {
				category = _("Utility");
			} else if (test_last_subdir((gchar *)iterator->data, "core")) {
				category = _("Core");
				core = TRUE;
			} else {
				category = _("Other");
			}
			gchar *scriptname = g_path_get_basename((gchar *)iterator->data);
			gchar *scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), (gchar *)iterator->data, NULL);
			const gchar *scripttype;
			if (g_str_has_suffix(scriptname, SCRIPT_EXT))
				scripttype = _("Siril Script File");
			else if (g_str_has_suffix(scriptname, PYSCRIPT_EXT) || g_str_has_suffix(scriptname, PYCSCRIPT_EXT))
				scripttype = _("Python script");
			else scripttype = NULL;

#ifdef DEBUG_GITSCRIPTS
			printf("%s\n", scriptpath);
#endif
			// Check whether the script appears in the list
			GList *iterator2 = NULL;
			if (!included && !core) {
				for (iterator2 = com.pref.selected_scripts; iterator2;
					iterator2 = iterator2->next) {
					if (g_strrstr((gchar *)iterator2->data, (gchar *)iterator->data)) {
					included = TRUE;
					}
				}
			}
			if (!core) {
				gtk_list_store_append(list_store, &iter);
				gtk_list_store_set(list_store, &iter, COLUMN_CATEGORY, category,
								COLUMN_SCRIPTNAME, scriptname, COLUMN_TYPE, scripttype, COLUMN_SELECTED,
								included, COLUMN_SCRIPTPATH, scriptpath,
								COLUMN_BGCOLOR, bg_color[color], -1);
			}
			/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
			g_free(scriptname);
			g_free(scriptpath); // it's ok to free this as the list_store keeps a copy internally
		}
	}
	gtk_tree_view_set_model(tview, GTK_TREE_MODEL(list_store));
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), sort_column_id, order);
	return FALSE;
}

/* called on preference window loading.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_script_repo_list(gboolean as_idle) {

	GtkTreeView *tview = GTK_TREE_VIEW(lookup_widget("treeview2"));
	if (as_idle)
		gdk_threads_add_idle(fill_script_repo_list_idle, tview);
	else
		fill_script_repo_list_idle(tview);
}

void on_treeview2_row_activated(GtkTreeView *treeview, GtkTreePath *path,
                                GtkTreeViewColumn *column, gpointer user_data) {
	gchar *scriptname = NULL, *scriptpath = NULL;
	gchar *contents = NULL;
	gsize length;
	GError *error = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model =
		gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("treeview2")));

	if (gtk_tree_model_get_iter(model, &iter, path)) {
		gtk_tree_model_get(model, &iter, 1, &scriptname, 3, &scriptpath, -1);
		if (g_file_get_contents(scriptpath, &contents, &length, &error) &&
			length > 0) {
			GtkTextBuffer *script_textbuffer = gtk_text_view_get_buffer(
				GTK_TEXT_VIEW(lookup_widget("script_contents")));
			GtkLabel *script_label = (GtkLabel *)lookup_widget("script_label");
			gtk_label_set_text(script_label, scriptname);
			gtk_text_buffer_set_text(script_textbuffer, contents, (gint)length);
			g_free(contents);
			siril_open_dialog("script_contents_dialog");
		} else {
			gchar *msg = g_strdup_printf(_("Error loading script contents: %s\n"), error->message);
			siril_log_color_message(msg, "red");
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			g_free(msg);
			g_error_free(error);
		}
	}
	g_free(scriptname);
	g_free(scriptpath);

}

void on_script_text_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("script_contents_dialog");
}

void on_manual_script_sync_button_clicked(GtkButton *button,
                                          gpointer user_data) {
	GString *git_pending_commit_buffer = NULL;
	set_cursor_waiting(TRUE);

	switch (preview_scripts_update(&git_pending_commit_buffer)) {
		case 1:
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Error getting the list of unmerged changes"));
			return;
		case 2:
			// Merge cannot be fast forwarded
			if (!siril_confirm_dialog(
					_("Warning!"),
					_("Merge analysis shows that "
					"the merge cannot be fast-forwarded. This indicates you have "
					"made changes to the local repository. Siril does not "
					"provide full git functionality and cannot be used to merge "
					"upstream updates into an altered local repository.\n\nIf you "
					"accept the update, the local repository will be hard reset "
					"to match the remote repository and any local changes will "
					"be lost.\n\nIf you have made local changes that you wish to "
					"keep, you should cancel this update and copy your modified "
					"scripts to another location, and add this location to the "
					"list of script directories to be searched."),
					_("Accept"))) {
				g_string_free(git_pending_commit_buffer, TRUE);
				return;
			} else {
				reset_scripts_repository();
				g_string_free(git_pending_commit_buffer, TRUE);
				fill_script_repo_list(FALSE);
				return;
			}
		default:
			break;
	}
	if (git_pending_commit_buffer != NULL) {
		if (siril_confirm_data_dialog(
				GTK_MESSAGE_QUESTION, _("Manual Update"),
				_("Read and confirm the pending changes to be synced"),
				_("Confirm"), git_pending_commit_buffer->str)) {
			if (reset_scripts_repository()) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Manual Update"), _("Error! Script database failed to update."));
			}
			fill_script_repo_list(FALSE);
		} else {
		siril_message_dialog(
			GTK_MESSAGE_INFO, _("Manual Update"),
			_("Update cancelled. Updates have not been applied."));
		}
		g_string_free(git_pending_commit_buffer, TRUE);
	} else {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("The script repository is up to date."));
	}
	fill_script_repo_list(FALSE);
	set_cursor_waiting(FALSE);
}

void on_manual_spcc_sync_button_clicked(GtkButton *button, gpointer user_data) {
	GString *git_pending_commit_buffer = NULL;
	set_cursor_waiting(TRUE);

	switch (preview_spcc_update(&git_pending_commit_buffer)) {
	case 1:
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Error getting the list of unmerged changes"));
		return;
	case 2:
		// Merge cannot be fast forwarded
		if (!siril_confirm_dialog(
				_("Warning!"),
				_("Merge analysis shows that "
				"the merge cannot be fast-forwarded. This indicates you have "
				"made changes to the local scripts repository. Siril does not "
				"provide full git functionality and cannot be used to merge "
				"upstream updates into an altered local repository.\n\nIf you "
				"accept the update, the local repository will be hard reset "
				"to match the remote repository and any local changes will "
				"be lost.\n\nIf you have made local changes that you wish to "
				"keep, you should cancel this update and copy your modified "
				"scripts to another location, and add this location to the "
				"list of script directories to be searched."),
				_("Accept"))) {
		g_string_free(git_pending_commit_buffer, TRUE);
		return;
		} else {
		reset_spcc_repository();
		g_string_free(git_pending_commit_buffer, TRUE);
		return;
		}
	default:
		break;
	}
	if (git_pending_commit_buffer != NULL) {
		if (siril_confirm_data_dialog(GTK_MESSAGE_QUESTION, _("Manual Update"),
				_("Read and confirm the pending changes to be synced"),
				_("Confirm"), git_pending_commit_buffer->str)) {
		if (reset_spcc_repository()) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Manual Update"), _("Error! SPCC database failed to update."));
		}
		} else {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("Update cancelled. Updates have not been applied."));
		}
		g_string_free(git_pending_commit_buffer, TRUE);
	} else {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("The SPCC database repository is up to date."));
	}
	if (!com.headless) {
		reset_spcc_filters();
		// Check if the SPCC window is open, if so refresh the combo boxes
		GtkWidget *spcc_dialog = lookup_widget("s_pcc_dialog");
		if (gtk_widget_get_visible(spcc_dialog)) {
		siril_debug_print("Reloading SPCC comboboxes\n");
		/* populate SPCC combos in a thread */
		g_thread_unref(
			g_thread_new("spcc_combos", populate_spcc_combos_async, NULL));
		}
	}
	set_cursor_waiting(FALSE);
}

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer, gchar *char_path, gpointer user_data) {
	gboolean val;
	GtkTreeIter iter;
	GtkTreePath *path;
	GtkTreeModel *model;
	gchar *script_path = NULL;
	path = gtk_tree_path_new_from_string(char_path);
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("treeview2")));
	if (gtk_tree_model_get_iter(model, &iter, path) == FALSE)
		return;
	gtk_tree_model_get(model, &iter, 3, &script_path, -1);
	gtk_tree_model_get(model, &iter, 2, &val, -1);
	gtk_list_store_set(GTK_LIST_STORE(model), &iter, 2, !val, -1);

	if (!val) {
		if (!(g_list_find(com.pref.selected_scripts, script_path))) {
#ifdef DEBUG_GITSCRIPTS
		printf("%s\n", script_path);
#endif
		com.pref.selected_scripts =
			g_list_prepend(com.pref.selected_scripts, script_path);
		}
	} else {
		GList *iterator = com.pref.selected_scripts;
		while (iterator) {
		if (g_strrstr((gchar *)iterator->data, script_path)) {
			iterator = g_list_remove_all(iterator, iterator->data);
			break;
		}
		iterator = iterator->next;
		}
		com.pref.selected_scripts = g_list_first(iterator);
	}
	notify_script_update();
}

void on_disable_gitscripts() {
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("treeview2")));
	GtkListStore *liststore = GTK_LIST_STORE(model);
	com.pref.use_scripts_repository = FALSE;
	gtk_list_store_clear(liststore);
	liststore = NULL;
	g_list_free_full(gui.repo_scripts, g_free);
	gui.repo_scripts = NULL;
	if (com.pref.selected_scripts)
		g_list_free_full(com.pref.selected_scripts, g_free);
	com.pref.selected_scripts = NULL;
	refresh_script_menu(TRUE);
}

void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(button)) {
		com.pref.use_scripts_repository = TRUE;
		auto_update_gitscripts(FALSE);
		fill_script_repo_list(FALSE);
	}
	gtk_widget_set_sensitive(lookup_widget("pref_script_automatic_updates"), com.pref.use_scripts_repository);
	gtk_widget_set_sensitive(lookup_widget("manual_script_sync_button"), (com.pref.use_scripts_repository && gui.script_repo_available));
	gtk_widget_set_sensitive(lookup_widget("treeview2"), (com.pref.use_scripts_repository && gui.script_repo_available));
}

void on_spcc_repo_enable_toggled(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(button)) {
		com.pref.spcc.use_spcc_repository = TRUE;
		auto_update_gitspcc(FALSE);
	}
	gtk_widget_set_sensitive(lookup_widget("pref_script_automatic_updates"), com.pref.spcc.use_spcc_repository);
	gtk_widget_set_sensitive(lookup_widget("spcc_repo_manual_sync"), (com.pref.spcc.use_spcc_repository && gui.spcc_repo_available));
}
#else

void hide_git_widgets() {
	gtk_widget_set_visible(lookup_widget("frame_gitscripts"), FALSE);
}

// We still need to provide placeholder callbacks to prevent GTK critical
// warnings, even though the widgets are hidden with libgit2 disabled

void on_pref_use_gitscripts_toggled(GtkToggleButton *button,
                                    gpointer user_data) {
	return;
}

void on_spcc_repo_enable_toggled(GtkToggleButton *button, gpointer user_data) {
	return;
}

void on_treeview2_row_activated(GtkTreeView *treeview, GtkTreePath *path,
                                GtkTreeViewColumn *column, gpointer user_data) {
	return;
}

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer,
                                   gchar *char_path, gpointer user_data) {
	return;
}

void on_manual_script_sync_button_clicked(GtkButton *button,
                                          gpointer user_data) {
	return;
}

void on_script_text_close_clicked(GtkButton *button, gpointer user_data) {
	return;
}

#endif
