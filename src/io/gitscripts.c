#include <git2.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui/utils.h"

static GtkListStore *list_store = NULL;

int update_gitscripts(void) {
    // Initialize libgit2
    git_libgit2_init();

    // URL of the remote repository
    const char *url = "https://gitlab.com/free-astro/siril-scripts.git";

    // Local directory where the repository will be cloned
    char *local_path = g_build_filename(g_get_user_data_dir(), "siril-scripts", NULL);

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
	}
	// Further operations can be performed on the 'repo' pointer here
	if (error == 0) {
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
			siril_log_color_message(_("Error looking up remote: %s\n"), "red", e->message);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
		}

		// Fetch
		error = git_remote_fetch(remote, NULL, &fetch_opts, NULL);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error fetching remote: %s\n"), "red", e->message);
			git_remote_free(remote);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
		}

		// Merge fetched changes into the local branch
		git_reference *branch_ref = NULL;
		error = git_branch_lookup(&branch_ref, repo, branch_name, GIT_BRANCH_LOCAL);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error looking up local branch: %s\n"), "red", e->message);
			git_remote_free(remote);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
		}

		git_reference *head_ref = NULL;
		error = git_repository_head(&head_ref, repo);

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error getting HEAD reference: %s\n"), "red", e->message);
			git_reference_free(branch_ref);
			git_remote_free(remote);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
		}

		git_annotated_commit *head_commit = NULL;
		error = git_annotated_commit_lookup(&head_commit, repo, git_reference_target(head_ref));

		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error looking up HEAD commit: %s\n"), "red", e->message);
			git_reference_free(head_ref);
			git_reference_free(branch_ref);
			git_remote_free(remote);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
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
			git_repository_free(repo);
			git_libgit2_shutdown();
			return 1;
		}

		siril_log_message(_("Changes pulled and merged successfully!\n"));

		/*** Populate the list of available repository scripts ***/
		size_t i;
		const git_index_entry *entry;
		git_index *index = NULL;
		if ((error = git_repository_index(&index, repo)) < 0)
			return 1;

		/* populate com.all_scripts with all the scripts in the index.
		 * We ignore anything not ending in .ssf */
		size_t entry_count = git_index_entrycount(index);
		g_autoptr(GStrvBuilder) builder = g_strv_builder_new();
		for (i = 0; i < entry_count; i++) {
			entry = git_index_get_byindex(index, i);
			if (g_str_has_suffix(entry->path, ".ssf")) {
				g_strv_builder_add(builder, entry->path);
				printf("%s\n", entry->path);
			}
		}
		gui.repo_scripts = g_strv_builder_end(builder);

	}

    // Cleanup
    git_repository_free(repo);
    git_libgit2_shutdown();

    return 0;
}

/************* GUI code for the Preferences->Scripts TreeView ****************/

enum {
	COLUMN_CATEGORY = 0,		// string
	COLUMN_SCRIPTPATH,		// string
	COLUMN_SELECTED,	// gboolean
	N_COLUMNS
};

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststore_script_repo"));
	}
}

gchar* get_script_filepath_from_path(GtkTreePath *path) {
	// TODO: Should actually return the script file path from the Path column of the GtkTreeView
	return NULL;
}

static gboolean fill_script_repo_list_idle(gpointer p) {
	int i;

	GtkTreeView* tview = (GtkTreeView*) p;
	GtkTreeIter iter;
	if (!tview)
		return FALSE;
	if (list_store) gtk_list_store_clear(list_store);
	get_list_store();
/*	gint sort_column_id;
	GtkSortType order;
	// store sorted state of list_store, disable sorting, disconnect from the view, fill, reconnect and re-apply sort
	gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(list_store), &sort_column_id, &order);
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID, GTK_SORT_ASCENDING);
	gtk_tree_view_set_model(tview, NULL);*/
	if (gui.repo_scripts) {
		i = 0;
		while (gui.repo_scripts[i]) {
			// TODO: here we should populate the GtkTreeView from GStrv gui.repo_scripts
			gchar* category = g_strrstr(gui.repo_scripts[i], "preprocessing") ? "Preprocessing" : "Processing";
			gchar* scriptpath = g_path_get_basename(gui.repo_scripts[i]);
			gboolean included = TRUE;
			gtk_list_store_append (list_store, &iter);
			gtk_list_store_set (list_store, &iter,
					COLUMN_CATEGORY, category,
					COLUMN_SCRIPTPATH, scriptpath,
					COLUMN_SELECTED, included,
					-1);
			/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
			g_free(scriptpath);
			i++;
		}
	}
/*	gtk_tree_view_set_model(tview, GTK_TREE_MODEL(list_store));
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), sort_column_id, order);*/

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

void on_script_list_active_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *char_path, gpointer user_data) {
	GtkTreePath *path = gtk_tree_path_new_from_string(char_path);
	gchar* scriptpath = get_script_filepath_from_path(path);
	gboolean active = gtk_cell_renderer_toggle_get_active(cell_renderer);
	gtk_tree_path_free(path);
	if (active) {
		// CHECK SCRIPT PATH IS IN LIST OF SCRIPTS IN SCRIPT MENU, IF NOT ADD IT
	} else {
		// CHECK SCRIPT PATH IS NOT IN LIST OF SCRIPTS IN SCRIPT MENU, IF IT IS REMOVE IT
	}
	// SAVE LIST OF ACTIVE SCRIPTS IN A FILE SOMEWHERE
}
