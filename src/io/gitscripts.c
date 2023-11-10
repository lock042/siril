#include <assert.h>
#include <inttypes.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/preferences.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/script_menu.h"
#include "core/siril_update.h" // for the version_number struct

//#define DEBUG_GITSCRIPTS

#ifdef HAVE_LIBGIT2
#include <git2.h>

const gchar *REPOSITORY_URL = "https://gitlab.com/free-astro/siril-scripts";

static GtkListStore *list_store = NULL;
static GString *git_pending_commit_buffer = NULL, *git_conflict_buffer = NULL;
static gboolean can_fastforward;

static void *xrealloc(void *oldp, size_t newsz)
{
	void *p = realloc(oldp, newsz);
	if (p == NULL) {
		PRINT_ALLOC_ERR;
		//exit(1);
	}
	return p;
}

enum {
	FORMAT_DEFAULT   = 0,
	FORMAT_LONG      = 1,
	FORMAT_SHORT     = 2,
	FORMAT_PORCELAIN = 3
};

struct merge_options {
	const char **heads;
	size_t heads_count;

	git_annotated_commit **annotated;
	size_t annotated_count;
};

static void merge_options_init(struct merge_options *opts)
{
	memset(opts, 0, sizeof(*opts));

	opts->heads = NULL;
	opts->heads_count = 0;
	opts->annotated = NULL;
	opts->annotated_count = 0;
}

static void opts_add_refish(struct merge_options *opts, const char *refish)
{
	size_t sz;

	assert(opts != NULL);

	sz = ++opts->heads_count * sizeof(opts->heads[0]);
	opts->heads = xrealloc((void *) opts->heads, sz);
	opts->heads[opts->heads_count - 1] = refish;
}

static int resolve_refish(git_annotated_commit **commit, git_repository *repo, const char *refish)
{
	git_reference *ref;
	git_object *obj;
	int err = 0;

	assert(commit != NULL);

	err = git_reference_dwim(&ref, repo, refish);
	if (err == GIT_OK) {
		git_annotated_commit_from_ref(commit, repo, ref);
		git_reference_free(ref);
		return 0;
	}

	err = git_revparse_single(&obj, repo, refish);
	if (err == GIT_OK) {
		err = git_annotated_commit_lookup(commit, repo, git_object_id(obj));
		git_object_free(obj);
	}

	return err;
}

static int resolve_heads(git_repository *repo, struct merge_options *opts)
{
	git_annotated_commit **annotated = calloc(opts->heads_count, sizeof(git_annotated_commit *));
	size_t annotated_count = 0, i;

	for (i = 0; i < opts->heads_count; i++) {
		int err = resolve_refish(&annotated[annotated_count++], repo, opts->heads[i]);
		if (err != 0) {
			siril_debug_print("libgit2: failed to resolve refish %s: %s\n", opts->heads[i], git_error_last()->message);
			annotated_count--;
			continue;
		}
	}

	if (annotated_count != opts->heads_count) {
		siril_log_color_message(_("libgit2: unable to parse some refish\n"), "red");
		free(annotated);
		return -1;
	}

	opts->annotated = annotated;
	opts->annotated_count = annotated_count;
	return 0;
}

static char *get_commit_from_oid(git_repository *repo, git_oid* oid_to_find) {
	char *retval = NULL;
	// Convert the OID to a commit
	git_commit *commit = NULL;
	if (git_commit_lookup(&commit, repo, oid_to_find) != 0) {
		siril_log_color_message(_("Failed to find commit with given OID.\n"), "red");
		return NULL;
	}

	// Iterate through references to find the ref name
	git_reference_iterator *ref_iterator = NULL;
	if (git_reference_iterator_new(&ref_iterator, repo) == 0) {
		git_reference *ref = NULL;
		while (git_reference_next(&ref, ref_iterator) == 0) {
			if (git_reference_type(ref) == GIT_REFERENCE_DIRECT) {
				const git_oid *target_oid = git_reference_target(ref);
				if (git_oid_equal(target_oid, git_commit_id(commit))) {
					retval = strdup(git_reference_name(ref));
					break;
				}
			}
			git_reference_free(ref);
		}
		git_reference_iterator_free(ref_iterator);
	}
	return retval;
}

static int fetchhead_cb(const char *ref_name, const char *remote_url, const git_oid *oid, unsigned int is_merge, void *payload)
{
    if (is_merge)
    {
        siril_debug_print("reference: '%s' is the reference we should merge\n", ref_name);
        git_oid_cpy((git_oid *)payload, oid);
    }
    return 0;
}

static int reset_repository() {
	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();

	// Initialisation
	git_libgit2_init();

	// Open the repository
	git_repository *repo = NULL;
	int error = git_repository_open(&repo, local_path);
	if (error != 0) {
		siril_log_color_message(_("Error performing hard reset. You may need to delete the local git repository and allow Siril to re-clone it.\n"), "red");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	// Get the FETCH_HEAD reference
	git_object *target_commit = NULL;
	error = 0;
	error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
	if (error != 0) {
		siril_log_color_message(_("Error performing hard reset. You may need to delete the local git repository and allow Siril to re-clone it.\n"), "red");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	// Perform the reset
	error = git_reset(repo, target_commit, GIT_RESET_HARD, NULL);
	if (error != 0) {
		siril_log_color_message(_("Error performing hard reset. You may need to delete the local git repository and allow Siril to re-clone it.\n"), "red");
		git_object_free(target_commit);
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("Success! The local repository is up-to-date with the remote."));

	git_repository_free(repo);
	git_libgit2_shutdown();
	return 0;
}

static int lg2_fetch(git_repository *repo)
{
	git_remote *remote = NULL;
	git_fetch_options fetch_opts = GIT_FETCH_OPTIONS_INIT;

	// Reference to the remote
	const char *remote_name = "origin";

	/* Figure out whether it's a named remote or a URL */
	siril_debug_print("Fetching %s for repo %p\n", remote_name, repo);

	if (git_remote_lookup(&remote, repo, remote_name)) {
		if (git_remote_create_anonymous(&remote, repo, remote_name)) {
			siril_log_message(_("Unable to create anonymous remote for the repository. Check your connectivity and try again.\n"));
			goto on_error;
		}
	}

	if (git_remote_fetch(remote, NULL, &fetch_opts, "fetch") < 0)
		siril_log_message(_("Error fetching remote. This may be a temporary error, check your connectivity and try again.\n"));

on_error:
	git_remote_free(remote);
	return -1;
}

static char* find_str_before_comment(const char* str1, const char* str2, const char* str3) {
	char* strpos = strstr(str1, str2);
	char* chrpos = strstr(str1, str3);
	return !strpos ? NULL : (chrpos && chrpos < strpos) ? NULL : strpos;
}

static gboolean script_version_check(const gchar* filename) {
	// Get the current version number
	gchar **fullVersionNumber = NULL;
	gchar **fullRequiresVersion = NULL;
	version_number version;
	fullVersionNumber = g_strsplit_set(PACKAGE_VERSION, ".-", -1);
	version.major_version = g_ascii_strtoull(fullVersionNumber[0], NULL, 10);
	version.minor_version = g_ascii_strtoull(fullVersionNumber[1], NULL, 10);
	version.micro_version = g_ascii_strtoull(fullVersionNumber[2], NULL, 10);

	// Open the script and look for the required version number
	GFile *file = NULL;
	GInputStream *stream = NULL;
	GDataInputStream *data_input = NULL;
	GError *error = NULL;
	gchar *buffer = NULL;
	gsize length = 0;
	gchar* scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), filename, NULL);
	gboolean retval = FALSE;
#ifdef DEBUG_GITSCRIPTS
	printf("checking script version requirements: %s\n", scriptpath);
#endif
	file = g_file_new_for_path(scriptpath);
	stream = (GInputStream*) g_file_read(file, NULL, &error);
	if (error)
		goto ERROR_OR_COMPLETE;
	data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, &error)) && !error) {
		gchar *ver = find_str_before_comment(buffer, "requires", "#");
		if (ver) {
			ver += 9;
			version_number requires;
			if (fullRequiresVersion)
				g_strfreev(fullRequiresVersion);
			fullRequiresVersion = g_strsplit_set(ver, ".-", -1);
			requires.major_version = g_ascii_strtoull(fullRequiresVersion[0], NULL, 10);
			requires.minor_version = g_ascii_strtoull(fullRequiresVersion[1], NULL, 10);
			requires.micro_version = g_ascii_strtoull(fullRequiresVersion[2], NULL, 10);
			// Detect badly formed requires command (bad input to g_ascii_strtoull returns 0) and ignore it
			if (requires.major_version == 0 && requires.minor_version == 0 && requires.micro_version == 0)
				continue;
#ifdef DEBUG_GITSCRIPTS
			printf("requires: %d.%d.%d; has %d.%d.%d\n", requires.major_version, requires.minor_version, requires.micro_version, version.major_version, version.minor_version, version.micro_version);
#endif
			if (requires.major_version < version.major_version) {
#ifdef DEBUG_GITSCRIPTS
				printf("requirement met\n");
#endif
				retval = TRUE;
			} else if (requires.major_version == version.major_version && requires.minor_version < version.minor_version) {
#ifdef DEBUG_GITSCRIPTS
				printf("requirement met\n");
#endif
				retval = TRUE;
			} else if (requires.major_version == version.major_version && requires.minor_version == version.minor_version &&
					 requires.micro_version <= version.micro_version) {
#ifdef DEBUG_GITSCRIPTS
				printf("requirement met\n");
#endif
				retval = TRUE;
			}
			g_free(buffer);
			buffer = NULL;
			if (retval)
				break;
		}
	}
ERROR_OR_COMPLETE:
	g_input_stream_close(stream, NULL, &error);
	if (error)
		siril_debug_print("Error closing data input stream from file\n");
	g_free(scriptpath);
	g_strfreev(fullVersionNumber);
	g_strfreev(fullRequiresVersion);
	g_object_unref(data_input);
	g_object_unref(stream);
	g_object_unref(file);
	return retval;
}

static int analyse(git_repository *repo) {
	git_remote *remote = NULL;

	// Carry out merge analysis
	struct merge_options opts;
	git_merge_analysis_t analysis;
	git_merge_preference_t preference;
	git_oid id_to_merge;
	git_repository_state_t state;
	int error = 0;

	// Figure out which branch to feed to git_merge()
	merge_options_init(&opts);
	git_repository_fetchhead_foreach(repo, fetchhead_cb, &id_to_merge);
	char *head_to_merge = get_commit_from_oid(repo, &id_to_merge);

	opts_add_refish(&opts, head_to_merge);
	state = git_repository_state(repo);
	if (state != GIT_REPOSITORY_STATE_NONE) {
		siril_log_color_message(_("libgit2: repository is in unexpected state %d. Cleaning up...\n"), "salmon", state);
		git_repository_state_cleanup(repo);
		state = git_repository_state(repo);
		if (state != GIT_REPOSITORY_STATE_NONE) {
			siril_log_color_message(_("libgit2: repository failed to clean up properly. Cannot continue.\n"), "red");
			free((char **)opts.heads);
			free(opts.annotated);
			return 1;
		}
	}

	error = resolve_heads(repo, &opts);
	if (error != 0) {
		free((char **)opts.heads);
		free(opts.annotated);
		return 1;
	}

	error = git_merge_analysis(&analysis, &preference,
	                         repo,
	                         (const git_annotated_commit **)opts.annotated,
	                         opts.annotated_count);

	// If the merge cannot be fast-forwarded, warn the user that local changes will
	// be lost if they proceed.
	if (error < 0) {
		siril_debug_print("Error carrying out merge analysis: %s\n",giterr_last()->message);
		return -1;
	}

	if ((analysis & GIT_MERGE_ANALYSIS_FASTFORWARD) || (analysis & GIT_MERGE_ANALYSIS_UP_TO_DATE)) {
		can_fastforward = TRUE;
	} else {
		can_fastforward = FALSE;
	}

	// If we already know we can't fast forward we can skip the rest of the function, we just return 2
	if (!can_fastforward) {
		siril_debug_print("Cannot be fast forwarded\n");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return 2;
	}

	//Prepare for looping through unmerged commit messages

	// Get the HEAD reference
	git_reference *head_ref = NULL;
	error = git_repository_head(&head_ref, repo);
	if (error != 0) {
		siril_debug_print("Error getting HEAD reference: %s\n", git_error_last()->message);
		goto on_error;
	}

	// Get the HEAD commit
	git_commit *head_commit = NULL;
	error = git_commit_lookup(&head_commit, repo, git_reference_target(head_ref));
	if (error != 0) {
		siril_debug_print("Error looking up HEAD commit: %s\n", git_error_last()->message);
		git_reference_free(head_ref);
		goto on_error;
	}

	// Get the reference to the FETCH_HEAD
	git_oid fetch_head_oid;
	if (git_reference_name_to_id(&fetch_head_oid, repo, "FETCH_HEAD") != 0) {
		siril_debug_print("Error getting FETCH_HEAD\n");
		goto on_error;
	}

	// Iterate through fetched commits and display commit messages
	git_commit *commit = NULL;
	git_oid parent_oid;
	git_commit *parent_commit = NULL;
	const char *commit_msg = NULL;

	// Start with the FETCH_HEAD
	if (git_commit_lookup(&commit, repo, &fetch_head_oid) != 0) {
		siril_debug_print("Error looking up commit\n");
		goto on_error;
	}

	gboolean found_head_ancestor = FALSE;

	while (1) {
		// Check if the current commit is the HEAD. If not, we print the commit message.
		// If so, we break and don't show any further messages.
		if (git_oid_equal(git_commit_id(head_commit), git_commit_id(commit))) {
			found_head_ancestor = TRUE;
			break;
		}

		// We have not yet reached the HEAD so we print the commit message.
		commit_msg = git_commit_message(commit);
		if (!git_pending_commit_buffer) {
			gchar* buf = g_strdup_printf(_("Commit message: %s\n"), commit_msg);
			git_pending_commit_buffer = g_string_new(buf);
			g_free(buf);
		} else {
			g_string_append_printf(git_pending_commit_buffer, _("Commit message: %s\n"), commit_msg);
		}

		if (git_commit_parentcount(commit) > 0) {
			parent_oid = *git_commit_parent_id(commit, 0);
			if (git_commit_lookup(&parent_commit, repo, &parent_oid) != 0) {
			break;
			}
		} else {
			break;
		}

		git_commit_free(commit);
		commit = parent_commit;
		parent_commit = NULL;
	}
	git_commit_free(commit);

	// If there is no ancestor commit found it indicates the merge cannot be fast-forwarded
	if (!found_head_ancestor) {
		can_fastforward = FALSE;
		git_remote_free(remote);
		return 2;
	}

	git_remote_free(remote);

	return 0;

 on_error:
	git_remote_free(remote);
	return -1;
}

int auto_update_gitscripts(gboolean sync) {
	int retval = 0;
	// Initialize libgit2
	git_libgit2_init();

    // URL of the remote repository
	siril_debug_print("Repository URL: %s\n", REPOSITORY_URL);

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
		error = git_clone(&repo, REPOSITORY_URL, local_path, &clone_opts);

		if (error != 0) {
			e = giterr_last();
			siril_log_color_message(_("Error cloning repository: %s\n"), "red", e->message);
			gui.script_repo_available = FALSE;
			git_libgit2_shutdown();
			return 1;
		} else {
			siril_log_message(_("Repository cloned successfully!\n"));
		}
	} else {
		siril_debug_print("Local scripts repository opened successfully!\n");
	}
	gui.script_repo_available = TRUE;

	// Check we are using the correct repository
	git_remote *remote = NULL;
	const char *remote_name = "origin";
	error = git_remote_lookup(&remote, repo, remote_name);
	if (error != 0) {
		siril_log_color_message(_("Failed to lookup remote.\n"), "red");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return 1;
	}

	const char *remote_url = git_remote_url(remote);
	if (remote_url != NULL) {
		siril_debug_print("Remote URL: %s\n", remote_url);
	} else {
		siril_log_color_message(_("Error: cannot identify local repository's configured remote.\n"), "red");
		return 1;
	}
	if (strcmp(remote_url, REPOSITORY_URL)) {
		gchar *msg = g_strdup_printf(_("Error: local siril-scripts repository folder is not "
				"configured with %s as its remote. You should remove the folder %s "
				"and restart Siril to re-clone the correct repository.\n"),
				REPOSITORY_URL, local_path);
		siril_log_color_message(msg, "red");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Repository Error"), msg);
		g_free(msg);
		// Make scripts unavailable, as the contents of a random git repository could be complete rubbish
		gui.script_repo_available = FALSE;
		goto cleanup;
	}

	// Synchronise the repository
	if (error == 0 && sync) {
		// fetch, analyse and merge changes from the remote
		lg2_fetch(repo);

		// Get the FETCH_HEAD reference
		git_object *target_commit = NULL;
		error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem persists you may need to delete the local git repository and allow Siril to re-clone it.\n"), "red");
			gui.script_repo_available = FALSE;
			git_repository_free(repo);
			git_libgit2_shutdown();
			return -1;
		}

		// Perform the reset
		error = git_reset(repo, target_commit, GIT_RESET_HARD, NULL);
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem persists you may need to delete the local git repository and allow Siril to re-clone it.\n"), "red");
			git_object_free(target_commit);
			git_repository_free(repo);
			git_libgit2_shutdown();
			return -1;
		}

		siril_log_color_message(_("Local scripts repository is up-to-date!\n"), "green");
	}

	/*** Populate the list of available repository scripts ***/
	size_t i;
	const git_index_entry *entry;
	git_index *index = NULL;
	error = git_repository_index(&index, repo);
	if (error < 0)
		retval = 1;

	/* populate gui.repo_scripts with all the scripts in the index.
		* We ignore anything not ending in SCRIPT_EXT */
	size_t entry_count = git_index_entrycount(index);
	if (gui.repo_scripts) {
		g_list_free_full(gui.repo_scripts, g_free);
		gui.repo_scripts = NULL;
	}
	for (i = 0; i < entry_count; i++) {
		entry = git_index_get_byindex(index, i);
		if (g_str_has_suffix(entry->path, SCRIPT_EXT) && script_version_check(entry->path)) {
			gui.repo_scripts = g_list_prepend(gui.repo_scripts, g_strdup(entry->path));
#ifdef DEBUG_GITSCRIPTS
			printf("%s\n", entry->path);
#endif
		}
	}

    // Cleanup
cleanup:
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
		siril_debug_print("Error opening repository: %s\n",giterr_last()->message);
		siril_log_color_message(_("Error: unable to open local scripts repository.\n"), "red");
		gui.script_repo_available = FALSE;
		return 1;
	}

	// Fetch changes
	lg2_fetch(repo);
	// Analyse the repository against the remote
	return analyse(repo);
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
		int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
		GList *iterator;
		for (iterator = gui.repo_scripts ; iterator ; iterator = iterator->next) {
			// here we populate the GtkTreeView from GList gui.repo_scripts
			gchar* category = g_strrstr((gchar*)iterator->data, "preprocessing") ? "Preprocessing" : "Processing";
			gchar* scriptname = g_path_get_basename((gchar*)iterator->data);
			gchar* scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), (gchar*)iterator->data, NULL);
#ifdef DEBUG_GITSCRIPTS
			printf("%s\n", scriptpath);
#endif
			// Check whether the script appears in the list
			GList* iterator2;
			gboolean included = FALSE;
			for (iterator2 = com.pref.selected_scripts ; iterator2 ; iterator2 = iterator2->next) {
				if (g_strrstr((gchar*) iterator2->data, (gchar*)iterator->data)) {
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
			g_free(contents);
			g_error_free(error);
			siril_open_dialog("script_contents_dialog");
		} else {
			gchar* msg = g_strdup_printf(_("Error loading script contents: %s\n"), error->message);
			siril_log_color_message(msg, "red");
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			g_free(msg);
			g_error_free(error);
		}
	}
}

void on_script_text_close_clicked(GtkButton* button, gpointer user_data) {
	siril_close_dialog("script_contents_dialog");
}

void on_manual_script_sync_button_clicked(GtkButton* button, gpointer user_data) {
	can_fastforward = FALSE;
	set_cursor_waiting(TRUE);
	if (git_pending_commit_buffer) {
		g_string_free(git_pending_commit_buffer, TRUE);
		git_pending_commit_buffer = NULL;
	}
	if (git_conflict_buffer) {
		g_string_free(git_conflict_buffer, TRUE);
		git_conflict_buffer = NULL;
	}

	switch (preview_update()) {
		case 1:
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Error getting the list of unmerged changes"));
			return;
		case 2:
			// Merge cannot be fast forwarded
			if (!siril_confirm_dialog(_("Warning!"), _("Merge analysis shows that "
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
				return;
			} else {
				reset_repository();
				return;
			}
		default:
			break;
	}
	if (git_pending_commit_buffer != NULL) {
		if (siril_confirm_data_dialog(GTK_MESSAGE_QUESTION, _("Manual Update"), _("Read and confirm the pending changes to be synced"), _("Confirm"), git_pending_commit_buffer->str)) {
			if (reset_repository()) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Manual Update"), _("Error! Scripts failed to update."));
			}
			fill_script_repo_list(FALSE);
		} else {
			siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("Update cancelled. Updates have not been applied."));
		}
	} else {
		siril_message_dialog(GTK_MESSAGE_INFO, _("Manual Update"), _("The scripts repository is up to date."));
	}
	fill_script_repo_list(TRUE);
	set_cursor_waiting(FALSE);
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
#ifdef DEBUG_GITSCRIPTS
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
	notify_script_update();
}

void on_disable_gitscripts() {
	GtkTreeModel *model = gtk_tree_view_get_model (GTK_TREE_VIEW(lookup_widget("treeview2")));
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

#else

void hide_git_widgets() {
	gtk_widget_set_visible(lookup_widget("frame_gitscripts"), FALSE);
}

// We still need to provide placeholder callbacks to prevent GTK critical warnings,
// even though the widgets are hidden with libgit2 disabled

void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data) {
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

void on_manual_script_sync_button_clicked(GtkButton* button, gpointer user_data) {
	return;
}

void on_script_text_close_clicked(GtkButton* button, gpointer user_data) {
	return;
}

#endif
