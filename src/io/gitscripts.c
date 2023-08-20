#include <git2.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/script_menu.h"

/* Uncomment the next line for additional debug messages showing scripts and
 * script paths
 */

//#define DEBUG_SCRIPTS

static GtkListStore *list_store = NULL;
static gchar *git_pending_commit_buffer = NULL;

int update_gitscripts(gboolean sync) {
	int retval = 0;
	// Initialize libgit2
    git_libgit2_init();

    // URL of the remote repository
    const char *url = "https://gitlab.com/free-astro/siril-scripts.git";

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();

	// Clone options
    git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;

    git_repository *repo = NULL;

	// See if the repository already exists
	int error = git_repository_open(&repo, local_path);

    if (error != 0) {
        const git_error *e = giterr_last();
        siril_log_color_message(_("Cannot open repository: %s\nAttempting to clone from remote source...\n"), "salmon", e->message);
		// Perform the clone operation
		error = git_clone(&repo, url, local_path, &clone_opts);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error cloning repository: %s\n"), "red", e->message);
		} else {
			siril_log_message(_("Repository cloned successfully!\n"));
		}
	} else {
		siril_log_message(_("Local scripts repository opened successfully!\n"));
	}
	// Synchronise the repository
	if (error == 0 && sync) {
		// Fetch options
		git_fetch_options fetch_opts = GIT_FETCH_OPTIONS_INIT;

		// Set up credentials for authentication if needed
		// fetch_opts.callbacks.credentials = your_credentials_callback_function;

		// Reference to the remote
		const char *remote_name = "origin";

		// Reference to the branch to pull from
		const char *branch_name = "main";

		// Remote
		git_remote *remote = NULL;
		error = git_remote_lookup(&remote, repo, remote_name);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Cannot look up remote: %s (this is normal if you are offline)"), "salmon", e->message);
			retval = 1;
			goto FINISH;
		}

		// Fetch
		error = git_remote_fetch(remote, NULL, &fetch_opts, NULL);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Cannot fetch remote: %s (this is normal if you are offline)\n"), "salmon", e->message);
			git_remote_free(remote);
			retval = 1;
			goto FINISH;
		}

		// Merge fetched changes into the local branch
		git_reference *branch_ref = NULL;
		error = git_branch_lookup(&branch_ref, repo, branch_name, GIT_BRANCH_LOCAL);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error looking up local branch: %s\n"), "red", e->message);
			git_remote_free(remote);
			retval = 1;
			goto FINISH;
		}

		git_reference *head_ref = NULL;
		error = git_repository_head(&head_ref, repo);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error getting HEAD reference: %s\n"), "red", e->message);
			git_reference_free(branch_ref);
			git_remote_free(remote);
			retval = 1;
			goto FINISH;
		}

		git_annotated_commit *head_commit = NULL;
		error = git_annotated_commit_lookup(&head_commit, repo, git_reference_target(head_ref));

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error looking up HEAD commit: %s\n"), "red", e->message);
			git_reference_free(head_ref);
			git_reference_free(branch_ref);
			git_remote_free(remote);
			retval = 1;
			goto FINISH;
		}

		git_reference *update_ref = NULL;
		error = git_reference_set_target(&update_ref, branch_ref, git_annotated_commit_id(head_commit), NULL);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error updating branch reference: %s\n"), "red", e->message);
			git_annotated_commit_free(head_commit);
			git_reference_free(head_ref);
			git_reference_free(branch_ref);
			git_remote_free(remote);
			retval = 1;
			goto FINISH;
		}

		siril_log_message(_("Changes pulled and merged successfully!\n"));
	}

FINISH: ;

	/*** Populate the list of available repository scripts ***/
	size_t i;
	const git_index_entry *entry;
	git_index *index = NULL;
	if ((error = git_repository_index(&index, repo)) < 0)
		retval = 1;

	/* populate com.all_scripts with all the scripts in the index.
		* We ignore anything not ending in .ssf */
	size_t entry_count = git_index_entrycount(index);
	g_autoptr(GStrvBuilder) builder = g_strv_builder_new();
	for (i = 0; i < entry_count; i++) {
		entry = git_index_get_byindex(index, i);
		if (g_str_has_suffix(entry->path, ".ssf")) {
			g_strv_builder_add(builder, entry->path);
#ifdef DEBUG_SCRIPTS
			printf("%s\n", entry->path);
#endif
		}
	}
	gui.repo_scripts = g_strv_builder_end(builder);

    // Cleanup
    if (repo)
		git_repository_free(repo);
    git_libgit2_shutdown();

    return retval;
}

int preview_update() {
	// Initialize libgit2
	git_libgit2_init();

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();

	git_repository *repo = NULL;
	int error = git_repository_open(&repo, local_path);
	if (error < 0) {
		siril_log_color_message(_("Error opening repository: %s\n"), "red", giterr_last()->message);
		return 1;
	}
	git_remote *remote = NULL;
	error = git_remote_lookup(&remote, repo, "origin");
	if (error < 0) {
		siril_log_color_message(_("Error looking up remote: %s\n"), "red", giterr_last()->message);
		git_repository_free(repo);
		return 1;
	}

	error = git_remote_fetch(remote, NULL, NULL, NULL);
	if (error < 0) {
		siril_log_color_message(_("Error fetching remote: %s\n"), "red", giterr_last()->message);
		git_remote_free(remote);
		git_repository_free(repo);
		return 1;
	}
	git_reference *remote_head_ref = NULL;
	error = git_reference_dwim(&remote_head_ref, repo, "origin/HEAD");
	if (error < 0) {
		siril_log_color_message(_("Error getting remote HEAD: %s\n"), "red", giterr_last()->message);
		git_remote_free(remote);
		git_repository_free(repo);
		return 1;
	}
	git_reference *local_head_ref = NULL;
	error = git_repository_head(&local_head_ref, repo);
	if (error < 0) {
		siril_log_color_message(_("Error getting local HEAD: %s\n"), "red", giterr_last()->message);
		git_reference_free(remote_head_ref);
		git_remote_free(remote);
		git_repository_free(repo);
		return 1;
	}
	git_oid remote_commit_id, local_commit_id;
	git_commit *remote_commit, *local_commit;

	git_oid_cpy(&remote_commit_id, git_reference_target(remote_head_ref));
	git_oid_cpy(&local_commit_id, git_reference_target(local_head_ref));

	if (git_pending_commit_buffer) {
		g_free(git_pending_commit_buffer);
		git_pending_commit_buffer = NULL;
	}
	if (git_commit_lookup(&remote_commit, repo, &remote_commit_id) == 0 &&
		git_commit_lookup(&local_commit, repo, &local_commit_id) == 0) {

		if (git_commit_time(remote_commit) > git_commit_time(local_commit)) {
			git_revwalk *walker = NULL;
			git_revwalk_new(&walker, repo);
			git_revwalk_push(walker, git_object_id((git_object *)remote_commit));

			git_oid commit_id;
			while (!git_revwalk_next(&commit_id, walker)) {
				git_commit *commit = NULL;
				if (git_commit_lookup(&commit, repo, &commit_id) == 0) {
					printf("Commit Message: %s\n", git_commit_summary(commit));
					gchar *temp = g_strdup_printf("%s\n", git_commit_summary(commit));
					gchar* temp2 = NULL;
					if (git_pending_commit_buffer)
						temp2 = g_strdup_printf("%s%s",git_pending_commit_buffer, temp);
					else
						temp2 = g_strdup_printf("%s", temp);
					g_free(temp);
					g_free(git_pending_commit_buffer);
					git_pending_commit_buffer = temp2;
					git_commit_free(commit);
				}
			}

			git_revwalk_free(walker);
		}

		git_commit_free(remote_commit);
		git_commit_free(local_commit);
	}

	siril_log_message(_("Completed checking pending commits to scripts repository...\n"));
	if (!git_pending_commit_buffer)
		siril_log_color_message(_("Local repository is up to date.\n"), "green");
	git_reference_free(local_head_ref);
	git_reference_free(remote_head_ref);
	git_remote_free(remote);
	git_repository_free(repo);
	git_libgit2_shutdown();
	return 0;
}

/************* GUI code for the Preferences->Scripts TreeView ****************/
static const char *bg_color[] = { "WhiteSmoke", "#1B1B1B" };

enum {
	COLUMN_CATEGORY = 0,		// string
	COLUMN_SCRIPTNAME,		// string
	COLUMN_SELECTED,	// gboolean
	COLUMN_SCRIPTPATH,	// full path to populate into the scripts menu
	COLUMN_BGCOLOR,		// background color
	N_COLUMNS
};

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststore_script_repo"));
	}
}

static gboolean fill_script_repo_list_idle(gpointer p) {
	int i;

	GtkTreeView* tview = (GtkTreeView*) p;
	GtkTreeIter iter;
	if (!tview)
		return FALSE;
	if (list_store) gtk_list_store_clear(list_store);
	get_list_store();
	gint sort_column_id;
	GtkSortType order;
	// store sorted state of list_store, disable sorting, disconnect from the view, fill, reconnect and re-apply sort
	gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(list_store), &sort_column_id, &order);
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID, GTK_SORT_ASCENDING);
	gtk_tree_view_set_model(tview, NULL);
	if (gui.repo_scripts) {
		i = 0;
		int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
		while (gui.repo_scripts[i]) {
			// here we populate the GtkTreeView from GStrv gui.repo_scripts
			gchar* category = g_strrstr(gui.repo_scripts[i], "preprocessing") ? "Preprocessing" : "Processing";
			gchar* scriptname = g_path_get_basename(gui.repo_scripts[i]);
#ifdef _WIN32
			gchar* scriptpath = g_build_path("\\", siril_get_scripts_repo_path(), gui.repo_scripts[i], NULL);
#else
			gchar* scriptpath = g_build_path("/", siril_get_scripts_repo_path(), gui.repo_scripts[i], NULL);
#endif
#ifdef DEBUG_SCRIPTS
			printf("%s\n", scriptpath);
#endif
			// Check whether the script appears in the list
			GList* iterator;
			gboolean included = FALSE;
			for (iterator = com.pref.selected_scripts ; iterator ; iterator = iterator->next) {
				if (g_strrstr((gchar*) iterator->data, gui.repo_scripts[i])) {
					included = TRUE;
				}
			}
			gtk_list_store_append (list_store, &iter);
			gtk_list_store_set (list_store, &iter,
					COLUMN_CATEGORY, category,
					COLUMN_SCRIPTNAME, scriptname,
					COLUMN_SELECTED, included,
					COLUMN_SCRIPTPATH, scriptpath,
					COLUMN_BGCOLOR, bg_color[color],
					-1);
			/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
			g_free(scriptpath);
			i++;
		}
	}
	gtk_tree_view_set_model(tview, GTK_TREE_MODEL(list_store));
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), sort_column_id, order);

	return FALSE;
}

/* called on preference window loading.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_script_repo_list(gboolean as_idle) {

	GtkTreeView* tview = GTK_TREE_VIEW(lookup_widget("treeview2"));
	if (as_idle)
		gdk_threads_add_idle(fill_script_repo_list_idle, tview);
	else fill_script_repo_list_idle(tview);
}

void on_treeview2_row_activated(GtkTreeView *treeview, GtkTreePath *path,
                    GtkTreeViewColumn *column, gpointer user_data) {
	gchar *scriptname = NULL, *scriptpath = NULL;
	gchar* contents = NULL;
	gsize length;
	GError* error = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model = gtk_tree_view_get_model (GTK_TREE_VIEW(lookup_widget("treeview2")));

	if (gtk_tree_model_get_iter(model, &iter, path)) {
		gtk_tree_model_get (model, &iter, 1, &scriptname, 3, &scriptpath, -1);
		if (g_file_get_contents(scriptpath, &contents, &length, &error) && length > 0) {
			GtkTextBuffer *script_textbuffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(lookup_widget("script_contents")));;
			GtkLabel* script_label = (GtkLabel*) lookup_widget("script_label");
			gtk_label_set_text(script_label, scriptname);
			gtk_text_buffer_set_text(script_textbuffer, contents, (gint) length);
			siril_open_dialog("script_contents_dialog");
		}
	}
}

void on_script_text_close_clicked(GtkButton* button, gpointer user_data) {
	siril_close_dialog("script_contents_dialog");
}

void on_manual_script_sync_button_clicked(GtkButton* button, gpointer user_data) {
	preview_update();
	if (git_pending_commit_buffer != NULL && siril_confirm_dialog(_("Read and accept the pending changes to be synced"), git_pending_commit_buffer, _("Confirm"))) {
		if (!update_gitscripts(TRUE)) {
			siril_message_dialog(GTK_MESSAGE_INFO, _("Update complete"), _("Scripts updated successfully."));
		} else {
			siril_message_dialog(GTK_MESSAGE_INFO, _("Error"), _("Scripts failed to update."));
		}
		fill_script_repo_list(FALSE);
	}
}

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *char_path, gpointer user_data) {
   gboolean val;
   GtkTreeIter iter;
   GtkTreePath *path;
   GtkTreeModel *model;
   gchar* script_path = NULL;
   path = gtk_tree_path_new_from_string(char_path);
   model = gtk_tree_view_get_model (GTK_TREE_VIEW(lookup_widget("treeview2")));
   if (gtk_tree_model_get_iter (model, &iter, path) == FALSE) return;
   gtk_tree_model_get(model, &iter, 3, &script_path, -1);
   gtk_tree_model_get(model, &iter, 2, &val, -1);
   gtk_list_store_set(GTK_LIST_STORE(model), &iter, 2, !val, -1);

	if (!val) {
		if (!(g_list_find(com.pref.selected_scripts, script_path))) {
#ifdef DEBUG_SCRIPTS
			printf("%s\n", script_path);
#endif
			com.pref.selected_scripts = g_list_prepend(com.pref.selected_scripts, script_path);
		}
	} else {
			GList* iterator = com.pref.selected_scripts;
			while (iterator) {
				if (g_strrstr((gchar*) iterator->data, script_path)) {
					iterator = g_list_remove_all(iterator, iterator->data);
					break;
				}
				iterator = iterator->next;
			}
			com.pref.selected_scripts = g_list_first(iterator);
	}
}

void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(button)) {
		com.pref.use_scripts_repository = TRUE;
		update_gitscripts(FALSE);
		fill_script_repo_list(FALSE);
	} else {
		com.pref.use_scripts_repository = FALSE;
		GtkTreeModel *model = gtk_tree_view_get_model (GTK_TREE_VIEW(lookup_widget("treeview2")));
		GtkListStore *liststore = GTK_LIST_STORE(model);
		gtk_list_store_clear(liststore);
		liststore = NULL;
//		g_list_free_full(com.pref.selected_scripts, g_free);
		g_strfreev(gui.repo_scripts);
		gui.repo_scripts = NULL;
		com.pref.selected_scripts = NULL;
		refresh_script_menu();
	}
}
