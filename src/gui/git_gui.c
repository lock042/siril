/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#include "gui/python_gui.h"
#include "io/siril_git.h"

#ifdef HAVE_LIBGIT2

static GtkListStore *list_store = NULL;
static gchar *current_search_text = NULL;
static gboolean filter_enabled = FALSE;
static GtkTreeModelFilter *filter_model = NULL;
static GtkTreeModelSort *sort_model = NULL;


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

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(
			gtk_builder_get_object(gui.builder, "liststore_script_repo"));
	}
}

void on_find_script_entry_changed(GtkEntry *entry, gpointer user_data) {
	const gchar *text = gtk_entry_get_text(entry);

	g_free(current_search_text);
	current_search_text = g_strdup(text);
	filter_enabled = (text && *text != '\0');

	if (filter_model) {
		gtk_tree_model_filter_refilter(filter_model);
	}
}

static gboolean tree_filter_visible_func(GtkTreeModel *model, GtkTreeIter *iter, gpointer data) {
	gchar *script_name = NULL;
	gboolean visible = TRUE;

	if (!filter_enabled || !current_search_text || strlen(current_search_text) == 0) {
		return TRUE;
	}

	gtk_tree_model_get(model, iter, COLUMN_SCRIPTNAME, &script_name, -1);

	if (script_name) {
		gchar *key_lower = g_ascii_strdown(current_search_text, -1);
		gchar *name_lower = g_ascii_strdown(script_name, -1);
		visible = (strstr(name_lower, key_lower) != NULL);
		g_free(key_lower);
		g_free(name_lower);
		g_free(script_name);
	}

	return visible;
}

static gboolean fill_script_repo_tree_idle(gpointer p) {
	GtkTreeView *tview = (GtkTreeView *)p;
	GtkTreeIter iter;
	if (!tview)
		return FALSE;

	if (sort_model) {
		gtk_tree_view_set_model(tview, NULL);
		g_object_unref(sort_model);
		sort_model = NULL;
	}
	if (filter_model) {
		g_object_unref(filter_model);
		filter_model = NULL;
	}

	if (list_store)
		gtk_list_store_clear(list_store);
	get_list_store();

	gint sort_column_id = GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID;
	GtkSortType order = GTK_SORT_ASCENDING;
	if (GTK_IS_TREE_SORTABLE(list_store)) {
		gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(list_store), &sort_column_id, &order);
		gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store),
			GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID,
			GTK_SORT_ASCENDING);
	}

	gui_repo_scripts_mutex_lock();
	if (gui.repo_scripts) {
		int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
		GSList *iterator;
		for (iterator = gui.repo_scripts; iterator; iterator = iterator->next) {
			// here we populate the GtkTreeView from GSList gui.repo_scripts
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
				// Extract last subdirectory and convert to title case
				gchar *path = (gchar *)iterator->data;
				gchar *last_slash = strrchr(path, G_DIR_SEPARATOR);

				if (last_slash && last_slash > path) {
					gchar *prev_slash = g_strrstr_len(path, last_slash - path, G_DIR_SEPARATOR_S);
					gchar *subdir_start = prev_slash ? prev_slash + 1 : path;

					gchar *subdir = g_strndup(subdir_start, last_slash - subdir_start);

					// Convert to title case
					if (subdir && subdir[0]) {
						subdir[0] = g_ascii_toupper(subdir[0]);
						for (gsize i = 1; subdir[i]; i++) {
							subdir[i] = g_ascii_tolower(subdir[i]);
						}
					}

					category = subdir;
				} else {
					category = _("Other");
				}
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
			GSList *iterator2 = NULL;
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
				gtk_list_store_set(list_store, &iter,
					COLUMN_CATEGORY, category,
					COLUMN_SCRIPTNAME, scriptname,
					COLUMN_TYPE, scripttype,
					COLUMN_SELECTED, included,
					COLUMN_SCRIPTPATH, scriptpath,
					COLUMN_BGCOLOR, bg_color[color],
					-1);
				// Free dynamically allocated category if it's not one of the static strings
				if (category != _("Preprocessing") &&
				    category != _("Processing") &&
				    category != _("Utility") &&
				    category != _("Core") &&
				    category != _("Other")) {
					g_free((gchar *)category);
				}
			}
			g_free(scriptname);
			g_free(scriptpath);
		}
	}
	gui_repo_scripts_mutex_unlock();

	filter_model = GTK_TREE_MODEL_FILTER(gtk_tree_model_filter_new(GTK_TREE_MODEL(list_store), NULL));
	gtk_tree_model_filter_set_visible_func(filter_model, tree_filter_visible_func, NULL, NULL);

	sort_model = GTK_TREE_MODEL_SORT(gtk_tree_model_sort_new_with_model(GTK_TREE_MODEL(filter_model)));
	gtk_tree_view_set_model(tview, GTK_TREE_MODEL(sort_model));

	if (GTK_IS_TREE_SORTABLE(list_store)) {
		gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), sort_column_id, order);
	}

	return FALSE;
}


/* called on preference window loading.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_script_repo_tree(gboolean as_idle) {

	GtkTreeView *tview = GTK_TREE_VIEW(lookup_widget("treeview_scripts"));
	if (as_idle)
		gdk_threads_add_idle(fill_script_repo_tree_idle, tview);
	else
		fill_script_repo_tree_idle(tview);
}

typedef struct {
	const gchar *scriptpath;
	GtkTextView *message_textview;
} CommitMessageUpdateData;

static void on_script_revision_spin_value_changed(GtkSpinButton *spin, gpointer user_data) {
	CommitMessageUpdateData *data = (CommitMessageUpdateData *)user_data;
	int revisions_back = gtk_spin_button_get_value_as_int(spin);

	GtkTextBuffer *buffer = gtk_text_view_get_buffer(data->message_textview);

	gchar *commit_message = NULL;
	size_t message_size = 0;
	gchar *content = get_script_content_string_from_file_revision(
		data->scriptpath, revisions_back, &(size_t){0}, &commit_message, &message_size);
	if (content != NULL) {
		g_free(content); // We only need the message
	}

	if (commit_message != NULL && message_size > 0) {
		gtk_text_buffer_set_text(buffer, commit_message, -1);
	} else {
		gtk_text_buffer_set_text(buffer, _("(No commit message found for this revision)"), -1);
	}
	g_free(commit_message);
}

void on_treeview_scripts_row_activated(GtkTreeView *treeview, GtkTreePath *path,
                                       GtkTreeViewColumn *column, gpointer user_data) {
	gchar *scriptname = NULL, *scriptpath = NULL;
	gchar *contents = NULL;
	gsize length;
	GtkTreeIter iter;
	GtkTreeModel *model =
		gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("treeview_scripts")));

	if (!gtk_tree_model_get_iter(model, &iter, path))
		return;

	gtk_tree_model_get(model, &iter, 1, &scriptname, 3, &scriptpath, -1);

	gchar *tmpscriptpath = g_canonicalize_filename(scriptpath, NULL);
	contents = get_script_content_string_from_file_revision(tmpscriptpath, 0, &length, NULL, NULL);
	if (length > 0 && contents != NULL) {
		const char *ext = get_filename_ext(tmpscriptpath);
		new_script(contents, length, ext);
		g_free(contents);
	} else {
		gchar *msg = g_strdup_printf(_("Error loading script contents: %s\n"), tmpscriptpath);
		siril_log_color_message(msg, "red");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
		g_free(msg);
	}

	g_free(scriptname);
	g_free(scriptpath);
}

// Callback function for right-click events
gboolean on_treeview_scripts_button_press(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	// Check if it's a right-click (button 3)
	if (event->type == GDK_BUTTON_PRESS && event->button == 3) {
		GtkTreeView *treeview = GTK_TREE_VIEW(widget);
		GtkTreePath *path = NULL;
		GtkTreeViewColumn *column = NULL;
		gint cell_x, cell_y;

		// Get the path at the click coordinates
		if (gtk_tree_view_get_path_at_pos(treeview,
										(gint)event->x,
										(gint)event->y,
										&path,
										&column,
										&cell_x,
										&cell_y)) {

			// We have a valid path - the click was on a row
			GtkTreeIter iter;
			GtkTreeModel *model = gtk_tree_view_get_model(treeview);

			if (gtk_tree_model_get_iter(model, &iter, path)) {
				gchar *scriptname = NULL, *scriptpath = NULL;
				gchar *contents = NULL;  // Declare contents here
				gsize length;            // Declare length here

				// Select the right-clicked row and make it active
				GtkTreeSelection *selection = gtk_tree_view_get_selection(treeview);
				gtk_tree_selection_select_path(selection, path);

				// Queue a redraw to show the selection immediately
				gtk_widget_queue_draw(GTK_WIDGET(treeview));

				gtk_tree_model_get(model, &iter, 1, &scriptname, 3, &scriptpath, -1);

				// Create dialog to ask for revision count
				GtkWidget *dialog = gtk_dialog_new_with_buttons(
					_("Select Revision"),
					GTK_WINDOW(lookup_widget("control_window")),
					GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
					_("_Cancel"), GTK_RESPONSE_CANCEL,
					_("_OK"), GTK_RESPONSE_OK,
					NULL);

				const gchar *tooltip = _("Leave at 0 to open the current version of the script. If you have a problem with "
					"a script update you can use this to go back to earlier revisions. Start by going back 1 revison "
					"and increase the number until you find the last version that worked for you. Note that this will "
					"only open the previous version in the script editor, it will not revert the script in the local "
					"repository, but it allows you to save a local copy in one of your script folders for use until the "
					"upstream script is fixed.");

				GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
				GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
				GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
				GtkWidget *label = gtk_label_new(_("Revisions to go back (leave at 0 to get current revision):"));
				gtk_widget_set_tooltip_text(label, tooltip);
				GtkWidget *spin_button = gtk_spin_button_new_with_range(0, 999, 1);
				gtk_widget_set_tooltip_text(spin_button, tooltip);
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_button), 0);

				gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
				gtk_box_pack_start(GTK_BOX(hbox), spin_button, FALSE, FALSE, 5);
				gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

				GtkWidget *message_heading = gtk_label_new(NULL);
				gtk_label_set_markup(GTK_LABEL(message_heading), "<b>Commit message</b>");
				gtk_label_set_xalign(GTK_LABEL(message_heading), 0.0);
				gtk_widget_set_margin_start(message_heading, 5);
				gtk_widget_set_margin_top(message_heading, 10);
				gtk_widget_set_margin_bottom(message_heading, 0);
				gtk_box_pack_start(GTK_BOX(vbox), message_heading, FALSE, FALSE, 5);

				// Create TextView with ScrolledWindow for commit message
				GtkWidget *scrolled_window = gtk_scrolled_window_new(NULL, NULL);
				gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window),
											  GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
				gtk_widget_set_size_request(scrolled_window, -1, 100); // ~5 lines height

				GtkWidget *message_textview = gtk_text_view_new();
				gtk_text_view_set_editable(GTK_TEXT_VIEW(message_textview), FALSE);
				gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(message_textview), FALSE);
				gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(message_textview), GTK_WRAP_WORD);
				gtk_widget_set_margin_start(scrolled_window, 5);
				gtk_widget_set_margin_end(scrolled_window, 5);
				gtk_widget_set_margin_top(scrolled_window, 2);
				gtk_widget_set_margin_bottom(scrolled_window, 5);

				// Set initial text
				GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(message_textview));
				gtk_text_buffer_set_text(buffer, _("Loading commit message..."), -1);

				gtk_container_add(GTK_CONTAINER(scrolled_window), message_textview);
				gtk_box_pack_start(GTK_BOX(vbox), scrolled_window, TRUE, TRUE, 5);

				// Add to dialog
				gtk_container_add(GTK_CONTAINER(content_area), vbox);
				gtk_container_set_border_width(GTK_CONTAINER(content_area), 10);
				gtk_widget_show_all(dialog);

				// Set up live updating of commit message
				CommitMessageUpdateData update_data = {
					.scriptpath = scriptpath,
					.message_textview = GTK_TEXT_VIEW(message_textview), // Changed from message_label
				};
				on_script_revision_spin_value_changed(GTK_SPIN_BUTTON(spin_button), &update_data); // Initial update
				g_signal_connect(spin_button, "value-changed", G_CALLBACK(on_script_revision_spin_value_changed), &update_data);

				gint response = gtk_dialog_run(GTK_DIALOG(dialog));

				// Check if user clicked OK before processing
				if (response == GTK_RESPONSE_OK) {
					// Get the spin button value BEFORE destroying the dialog
					int revisions_back = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin_button));
					gtk_widget_destroy(dialog);

					// Now use the revisions_back value, not the response
					contents = get_script_content_string_from_file_revision(scriptpath, revisions_back, &length, NULL, NULL);

					if (length > 0 && contents != NULL) {
						const char *ext = get_filename_ext(scriptpath);
						new_script(contents, length, ext);
						g_free(contents);
					} else {
						gchar *msg = g_strdup_printf(_("Error loading script contents from %d revisions back: %s\n"),
							revisions_back, scriptpath);
						siril_log_color_message(msg, "red");
						siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
						g_free(msg);
					}
				} else {
					// User cancelled - just destroy the dialog
					gtk_widget_destroy(dialog);
				}

				g_free(scriptname);
				g_free(scriptpath);
			}

			gtk_tree_path_free(path);
			return TRUE; // Event handled
		}

		// Right-click on empty area - you might want to show a different menu
		// or just ignore it
		return TRUE;
	}

	return FALSE; // Let other handlers process the event
}

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer, gchar *char_path, gpointer user_data) {
	gboolean val;
	GtkTreeIter iter, filter_iter, child_iter;
	GtkTreePath *path;
	GtkTreeModel *model;
	gchar *script_path = NULL;

	path = gtk_tree_path_new_from_string(char_path);
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("treeview_scripts")));

	if (gtk_tree_model_get_iter(model, &iter, path) == FALSE) {
		gtk_tree_path_free(path);
		return;
	}

	// Free the path as soon as we don't need it anymore
	gtk_tree_path_free(path);

	// Get the script path and current toggle value
	gtk_tree_model_get(model, &iter, 3, &script_path, -1);
	gtk_tree_model_get(model, &iter, 2, &val, -1);

	if (GTK_IS_TREE_MODEL_SORT(model)) {
		gtk_tree_model_sort_convert_iter_to_child_iter(GTK_TREE_MODEL_SORT(model), &filter_iter, &iter);
		if (GTK_IS_TREE_MODEL_FILTER(gtk_tree_model_sort_get_model(GTK_TREE_MODEL_SORT(model)))) {
			GtkTreeModelFilter *filter = GTK_TREE_MODEL_FILTER(gtk_tree_model_sort_get_model(GTK_TREE_MODEL_SORT(model)));
			gtk_tree_model_filter_convert_iter_to_child_iter(filter, &child_iter, &filter_iter);
			// Toggle the value in the actual list_store
			gtk_list_store_set(list_store, &child_iter, 2, !val, -1);
		}
	} else {
		gtk_list_store_set(GTK_LIST_STORE(model), &iter, 2, !val, -1);
	}

	if (!val) {
		// Checkbox is now checked - add to list if not already present
		if (!g_slist_find_custom(com.pref.selected_scripts, script_path, (GCompareFunc)g_strcmp0)) {
#ifdef DEBUG_GITSCRIPTS
			printf("Adding script: %s\n", script_path);
#endif
			// g_slist_prepend takes ownership of script_path
			com.pref.selected_scripts = g_slist_prepend(com.pref.selected_scripts, script_path);
		} else {
			// Already in list, siril_free our copy
			g_free(script_path);
		}
	} else {
		// Checkbox is now unchecked - remove from list
		GSList *found = g_slist_find_custom(com.pref.selected_scripts, script_path, (GCompareFunc)g_strcmp0);
		if (found) {
#ifdef DEBUG_GITSCRIPTS
			printf("Removing script: %s\n", script_path);
#endif
			// Free the data stored in the list
			g_free(found->data);
			// Remove the element and update the list
			com.pref.selected_scripts = g_slist_delete_link(com.pref.selected_scripts, found);
		}
		// Free our copy of script_path
		g_free(script_path);
	}

	notify_script_update();
}

void on_disable_gitscripts() {
	// Clear the actual list_store, not the filter
	com.pref.use_scripts_repository = FALSE;
	if (list_store) {
		gtk_list_store_clear(list_store);
	}

	gui_repo_scripts_mutex_lock();
	g_slist_free_full(gui.repo_scripts, g_free);
	gui_repo_scripts_mutex_unlock();

	gui.repo_scripts = NULL;
	if (com.pref.selected_scripts)
		g_slist_free_full(com.pref.selected_scripts, g_free);
	com.pref.selected_scripts = NULL;
	g_thread_unref(g_thread_new("refresh_script_menu", refresh_script_menu_in_thread, GINT_TO_POINTER(1)));
}

void on_manual_script_sync_button_clicked(GtkButton *button, gpointer user_data) {
	g_thread_unref(g_thread_new("update_scripts", update_scripts, NULL));
}

void on_manual_spcc_sync_button_clicked(GtkButton *button, gpointer user_data) {
	g_thread_unref(g_thread_new("update_spcc", update_spcc, NULL));
}

void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.use_scripts_repository = gtk_toggle_button_get_active(button);
	if (com.pref.use_scripts_repository) {
		g_thread_unref(g_thread_new("update_scripts", initialize_scripts, NULL));
	}
	gtk_widget_set_sensitive(lookup_widget("pref_script_automatic_updates"), com.pref.use_scripts_repository);
	gtk_widget_set_sensitive(lookup_widget("manual_script_sync_button"), (com.pref.use_scripts_repository && gui.script_repo_available));
	gtk_widget_set_sensitive(lookup_widget("treeview_scripts"), (com.pref.use_scripts_repository && gui.script_repo_available));
}

void on_spcc_repo_enable_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.spcc.use_spcc_repository = gtk_toggle_button_get_active(button);
	if (com.pref.spcc.use_spcc_repository) {
		g_thread_unref(g_thread_new("update_spcc", initialize_spcc, NULL));
	}
	gtk_widget_set_sensitive(lookup_widget("spcc_repo_sync_at_startup"), com.pref.spcc.use_spcc_repository);
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

void on_treeview_scripts_row_activated(GtkTreeView *treeview, GtkTreePath *path,
                                GtkTreeViewColumn *column, gpointer user_data) {
	return;
}

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer,
                                   gchar *char_path, gpointer user_data) {
	return;
}

void on_script_text_close_clicked(GtkButton *button, gpointer user_data) {
	return;
}

void on_manual_script_sync_button_clicked(GtkButton *button, gpointer user_data) {
	return;
}

void on_manual_spcc_sync_button_clicked(GtkButton *button, gpointer user_data) {
	return;
}

#endif
