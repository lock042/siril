#include <git2.h>
#include <assert.h>
#include <inttypes.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/script_menu.h"

/* Uncomment the next line for additional debug messages showing scripts and
 * script paths
 */

//#define DEBUG_SCRIPTS

static GtkListStore *list_store = NULL;
static GString *git_pending_commit_buffer = NULL;

static void *xrealloc(void *oldp, size_t newsz)
{
	void *p = realloc(oldp, newsz);
	if (p == NULL) {
		PRINT_ALLOC_ERR;
		//exit(1);
	}
	return p;
}

#define MAX_PATHSPEC 8

enum {
	FORMAT_DEFAULT   = 0,
	FORMAT_LONG      = 1,
	FORMAT_SHORT     = 2,
	FORMAT_PORCELAIN = 3
};

struct status_opts {
	git_status_options statusopt;
	char *repodir;
	char *pathspec[MAX_PATHSPEC];
	int npaths;
	int format;
	int zterm;
	int showbranch;
	int showsubmod;
	int repeat;
};

struct merge_options {
	const char **heads;
	size_t heads_count;

	git_annotated_commit **annotated;
	size_t annotated_count;

	int no_commit : 1;
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

void check_lg2(int error, const char *message, const char *extra)
{
	const git_error *lg2err;
	const char *lg2msg = "", *lg2spacer = "";

	if (!error)
		return;

	if ((lg2err = git_error_last()) != NULL && lg2err->message != NULL) {
		lg2msg = lg2err->message;
		lg2spacer = " - ";
	}

	if (extra) {
		siril_debug_print("%s '%s' [%d]%s%s\n",
			message, extra, error, lg2spacer, lg2msg);
		siril_log_color_message(_("%s [%d]%s%s\n"), "red",
			message, error, lg2spacer, lg2msg);
	} else {
		siril_log_color_message(_("%s [%d]%s%s\n"), "red",
			message, error, lg2spacer, lg2msg);
	}
//	exit(1);
}

int resolve_refish(git_annotated_commit **commit, git_repository *repo, const char *refish)
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
	int err = 0;

	for (i = 0; i < opts->heads_count; i++) {
		err = resolve_refish(&annotated[annotated_count++], repo, opts->heads[i]);
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

static int perform_fastforward(git_repository *repo, const git_oid *target_oid, int is_unborn)
{
	git_checkout_options ff_checkout_options = GIT_CHECKOUT_OPTIONS_INIT;
	git_reference *target_ref;
	git_reference *new_target_ref;
	git_object *target = NULL;
	int err = 0;

	if (is_unborn) {
		const char *symbolic_ref;
		git_reference *head_ref;

		/* HEAD reference is unborn, lookup manually so we don't try to resolve it */
		err = git_reference_lookup(&head_ref, repo, "HEAD");
		if (err != 0) {
			siril_log_color_message(_("libgit2: failed to lookup HEAD ref\n"), "red");
			return -1;
		}

		/* Grab the reference HEAD should be pointing to */
		symbolic_ref = git_reference_symbolic_target(head_ref);

		/* Create our master reference on the target OID */
		err = git_reference_create(&target_ref, repo, symbolic_ref, target_oid, 0, NULL);
		if (err != 0) {
			siril_log_color_message(_("libgit2: failed to create master reference\n"), "red");
			return -1;
		}

		git_reference_free(head_ref);
	} else {
		/* HEAD exists, just lookup and resolve */
		err = git_repository_head(&target_ref, repo);
		if (err != 0) {
			siril_log_color_message(_("libgit2: failed to get HEAD reference\n"), "red");
			return -1;
		}
	}

	/* Lookup the target object */
	err = git_object_lookup(&target, repo, target_oid, GIT_OBJECT_COMMIT);
	if (err != 0) {
		siril_log_color_message(_("libgit2: failed to lookup OID %s\n"), "red", git_oid_tostr_s(target_oid));
		return -1;
	}

	/* Checkout the result so the workdir is in the expected state */
	ff_checkout_options.checkout_strategy = GIT_CHECKOUT_SAFE;
	err = git_checkout_tree(repo, target, &ff_checkout_options);
	if (err != 0) {
		siril_log_color_message(_("libgit2: failed to checkout HEAD reference\n"), "red");
		return -1;
	}

	/* Move the target reference to the target OID */
	err = git_reference_set_target(&new_target_ref, target_ref, target_oid, NULL);
	if (err != 0) {
		siril_log_color_message(_("libgit2: failed to move HEAD reference\n"), "red");
		return -1;
	}

	git_reference_free(target_ref);
	git_reference_free(new_target_ref);
	git_object_free(target);

	return 0;
}

static void output_conflicts(git_index *index)
{
	git_index_conflict_iterator *conflicts;
	const git_index_entry *ancestor;
	const git_index_entry *our;
	const git_index_entry *their;
	int err = 0;

	check_lg2(git_index_conflict_iterator_new(&conflicts, index), _("libgit2: failed to create conflict iterator"), NULL);

	while ((err = git_index_conflict_next(&ancestor, &our, &their, conflicts)) == 0) {
		siril_log_color_message(_("libgit2: conflict: a:%s o:%s t:%s\n"), "red",
		        ancestor ? ancestor->path : "NULL",
		        our->path ? our->path : "NULL",
		        their->path ? their->path : "NULL");
	}

	if (err != GIT_ITEROVER) {
		siril_log_color_message(_("error iterating conflicts\n"), "red");
	}

	git_index_conflict_iterator_free(conflicts);
}

static int create_merge_commit(git_repository *repo, git_index *index, struct merge_options *opts)
{
	git_oid tree_oid, commit_oid;
	git_tree *tree;
	git_signature *sign;
	git_reference *merge_ref = NULL;
	git_annotated_commit *merge_commit;
	git_reference *head_ref;
	git_commit **parents = calloc(opts->annotated_count + 1, sizeof(git_commit *));
	const char *msg_target = NULL;
	size_t msglen = 0;
	char *msg;
	size_t i;
	int err;

	/* Grab our needed references */
	check_lg2(git_repository_head(&head_ref, repo), _("libgit2: failed to get repo HEAD"), NULL);
	if (resolve_refish(&merge_commit, repo, opts->heads[0])) {
		siril_log_color_message(_("libgit2: failed to resolve refish %s"), "red", opts->heads[0]);
		free(parents);
		return -1;
	}

	/* Maybe that's a ref, so DWIM it */
	err = git_reference_dwim(&merge_ref, repo, opts->heads[0]);
	check_lg2(err, _("libgit2: failed to DWIM reference"), git_error_last()->message);

	/* Grab a signature */
	// TODO: I doubt this works as-is from the example.
	check_lg2(git_signature_now(&sign, "Me", "me@example.com"), _("libgit2: failed to create signature"), NULL);

#define MERGE_COMMIT_MSG "Merge %s '%s'"
	/* Prepare a standard merge commit message */
	if (merge_ref != NULL) {
		check_lg2(git_branch_name(&msg_target, merge_ref), _("libgit2: failed to get branch name of merged ref"), NULL);
	} else {
		msg_target = git_oid_tostr_s(git_annotated_commit_id(merge_commit));
	}

	msglen = snprintf(NULL, 0, MERGE_COMMIT_MSG, (merge_ref ? "branch" : "commit"), msg_target);
	if (msglen > 0) msglen++;
	msg = malloc(msglen);
	err = snprintf(msg, msglen, MERGE_COMMIT_MSG, (merge_ref ? "branch" : "commit"), msg_target);

	/* This is only to silence the compiler */
	if (err < 0) goto cleanup;

	/* Setup our parent commits */
	err = git_reference_peel((git_object **)&parents[0], head_ref, GIT_OBJECT_COMMIT);
	check_lg2(err, _("failed to peel head reference"), NULL);
	for (i = 0; i < opts->annotated_count; i++) {
		git_commit_lookup(&parents[i + 1], repo, git_annotated_commit_id(opts->annotated[i]));
	}

	/* Prepare our commit tree */
	check_lg2(git_index_write_tree(&tree_oid, index), _("failed to write merged tree"), NULL);
	check_lg2(git_tree_lookup(&tree, repo, &tree_oid), _("failed to lookup tree"), NULL);

	/* Commit time ! */
	err = git_commit_create(&commit_oid,
	                        repo, git_reference_name(head_ref),
	                        sign, sign,
	                        NULL, msg,
	                        tree,
	                        opts->annotated_count + 1, (const git_commit **)parents);
	check_lg2(err, _("failed to create commit"), NULL);

	/* We're done merging, cleanup the repository state */
	git_repository_state_cleanup(repo);

cleanup:
	free(parents);
	return err;
}

char *get_commit_from_oid(git_repository *repo, git_oid* oid_to_find) {
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

int fetchhead_cb(const char *ref_name, const char *remote_url, const git_oid *oid, unsigned int is_merge, void *payload)
{
    if (is_merge)
    {
        siril_debug_print("reference: '%s' is the reference we should merge\n", ref_name);
        git_oid_cpy((git_oid *)payload, oid);
    }
    return 0;
}

int lg2_merge(git_repository *repo) {
	struct merge_options opts;
	git_index *index;
	git_repository_state_t state;
	git_merge_analysis_t analysis;
	git_merge_preference_t preference;
	git_oid id_to_merge;
	int err = 0;

	merge_options_init(&opts);

	// Figure out which branch to feed to git_merge()
	git_repository_fetchhead_foreach(repo, fetchhead_cb, &id_to_merge);
	char *head_to_merge = get_commit_from_oid(repo, &id_to_merge);

	opts_add_refish(&opts, head_to_merge);
	state = git_repository_state(repo);
	if (state != GIT_REPOSITORY_STATE_NONE) {
		siril_log_color_message(_("libgit2: repository is in unexpected state %d\n"), "red", state);
		goto cleanup;
	}

	err = resolve_heads(repo, &opts);
	if (err != 0)
		goto cleanup;

	err = git_merge_analysis(&analysis, &preference,
	                         repo,
	                         (const git_annotated_commit **)opts.annotated,
	                         opts.annotated_count);
	if (err != 0) {
		check_lg2(err, _("libgit2: merge analysis failed"), NULL);
		goto cleanup;
	}

	if (analysis & GIT_MERGE_ANALYSIS_UP_TO_DATE) {
		siril_log_color_message(_("Already up-to-date\n"), "green");
		return 0;
	} else if (analysis & GIT_MERGE_ANALYSIS_UNBORN ||
	          (analysis & GIT_MERGE_ANALYSIS_FASTFORWARD &&
	          !(preference & GIT_MERGE_PREFERENCE_NO_FASTFORWARD))) {
		const git_oid *target_oid;
		if (analysis & GIT_MERGE_ANALYSIS_UNBORN) {
			siril_debug_print("Unborn\n");
		} else {
			siril_debug_print("Fast-forward\n");
		}

		/* Since this is a fast-forward, there can be only one merge head */
		target_oid = git_annotated_commit_id(opts.annotated[0]);
		assert(opts.annotated_count == 1);

		return perform_fastforward(repo, target_oid, (analysis & GIT_MERGE_ANALYSIS_UNBORN));
	} else if (analysis & GIT_MERGE_ANALYSIS_NORMAL) {
		git_merge_options merge_opts = GIT_MERGE_OPTIONS_INIT;
		git_checkout_options checkout_opts = GIT_CHECKOUT_OPTIONS_INIT;

		merge_opts.flags = 0;
		merge_opts.file_flags = GIT_MERGE_FILE_STYLE_DIFF3;

		checkout_opts.checkout_strategy = GIT_CHECKOUT_FORCE|GIT_CHECKOUT_ALLOW_CONFLICTS;

		if (preference & GIT_MERGE_PREFERENCE_FASTFORWARD_ONLY) {
			siril_log_color_message(_("libgit2: fast-forward is preferred, but only a merge is possible\n"), "red");
			return -1;
		}

		err = git_merge(repo,
		                (const git_annotated_commit **)opts.annotated, opts.annotated_count,
		                &merge_opts, &checkout_opts);
		check_lg2(err, _("libgit2: merge failed"), NULL);
	}

	/* If we get here, we actually performed the merge above */

	check_lg2(git_repository_index(&index, repo), _("libgit2: failed to get repository index"), NULL);

	if (git_index_has_conflicts(index)) {
		/* Handle conflicts */
		output_conflicts(index);
	} else if (!opts.no_commit) {
		create_merge_commit(repo, index, &opts);
		siril_log_color_message(_("Merge completed\n"), "green");
	}

cleanup:
	free((char **)opts.heads);
	free(opts.annotated);

	return 0;
}

int lg2_fetch(git_repository *repo)
{
	git_remote *remote = NULL;
//	const git_indexer_progress *stats;
	git_fetch_options fetch_opts = GIT_FETCH_OPTIONS_INIT;

	// Reference to the remote
	const char *remote_name = "origin";

	/* Figure out whether it's a named remote or a URL */
	siril_debug_print("Fetching %s for repo %p\n", remote_name, repo);

	if (git_remote_lookup(&remote, repo, remote_name))
//		if (git_remote_create_anonymous(&remote, repo, remote_name))
			goto on_error;

	if (git_remote_fetch(remote, NULL, &fetch_opts, "fetch") < 0)
		goto on_error;

	git_remote_free(remote);

	return 0;

 on_error:
	git_remote_free(remote);
	return -1;
}

int update_gitscripts(gboolean sync) {
	int retval = 0;
	// Initialize libgit2
    git_libgit2_init();

    // URL of the remote repository
//    const char *url = "https://gitlab.com/free-astro/siril-scripts.git";
    const char *url = "https://gitlab.com/aje.baugh/siril-scripts.git"; // For testing purposes

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
			git_libgit2_shutdown();
			return 1;
		} else {
			siril_log_message(_("Repository cloned successfully!\n"));
		}
	} else {
		siril_log_message(_("Local scripts repository opened successfully!\n"));
	}
	// Synchronise the repository
	if (error == 0 && sync) {
		// fetch and merge changed from the remote
		lg2_fetch(repo);

		lg2_merge(repo);

		siril_log_color_message(_("Changes fetched and merged successfully!\n"), "green");
	}

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
		siril_debug_print("Error opening repository: %s\n",giterr_last()->message);
		return 1;
	}

	// Fetch changes
	lg2_fetch(repo);

	//Prepare for looping through unmerged commit messages

	const char *remote_name = "origin";
	git_remote *remote = NULL;

	// Open the remote
	if (git_remote_lookup(&remote, repo, remote_name) != 0) {
		siril_debug_print("Error looking up remote\n");
		git_repository_free(repo);
		return 1;
	}

	// Get the HEAD reference
	git_reference *head_ref = NULL;
	error = git_repository_head(&head_ref, repo);
	if (error != 0) {
		siril_debug_print("Error getting HEAD reference: %s\n", git_error_last()->message);
		git_repository_free(repo);
		return 1;
	}

	// Get the HEAD commit
	git_commit *head_commit = NULL;
	error = git_commit_lookup(&head_commit, repo, git_reference_target(head_ref));
	if (error != 0) {
		siril_debug_print("Error looking up HEAD commit: %s\n", git_error_last()->message);
		git_reference_free(head_ref);
		git_repository_free(repo);
		return 1;
	}

	// Get the reference to the FETCH_HEAD
	git_oid fetch_head_oid;
	if (git_reference_name_to_id(&fetch_head_oid, repo, "FETCH_HEAD") != 0) {
		siril_debug_print("Error getting FETCH_HEAD\n");
		git_remote_free(remote);
		git_repository_free(repo);
		return 1;
	}

	// Iterate through fetched commits and display commit messages
	git_commit *commit = NULL;
	git_oid parent_oid;
	git_commit *parent_commit = NULL;
	const char *commit_msg = NULL;

	if (git_commit_lookup(&commit, repo, &fetch_head_oid) != 0) {
		siril_debug_print("Error looking up commit\n");
		git_remote_free(remote);
		git_repository_free(repo);
		git_libgit2_shutdown();
		return 1;
	}

	while (1) {
		// Check if the current commit is the HEAD. If not, we print the commit message.
		// If so, we break and don't show any further messages.
		if (git_oid_equal(git_commit_id(head_commit), git_commit_id(commit)))
			break;

		// We have not yet reached the HEAD so we print the commit message.
		commit_msg = git_commit_message(commit);
		if (!git_pending_commit_buffer) {
			gchar *tmp = g_strdup_printf(_("Commit message: %s\n"), commit_msg);
			git_pending_commit_buffer = g_string_new(tmp);
			g_free(tmp);
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

	if (!git_pending_commit_buffer)
		siril_log_color_message(_("Local repository is up to date.\n"), "green");
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
	if (git_pending_commit_buffer) {
		g_string_free(git_pending_commit_buffer, TRUE);
		git_pending_commit_buffer = NULL;
	}
	if (preview_update()) {
		siril_log_color_message(_("Error getting the list of unmerged changes.\n"), "red");
	}
	if (git_pending_commit_buffer != NULL && siril_confirm_dialog(_("Read and accept the pending changes to be synced"), git_pending_commit_buffer->str, _("Confirm"))) {
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
		g_strfreev(gui.repo_scripts);
		gui.repo_scripts = NULL;
		com.pref.selected_scripts = NULL;
		refresh_script_menu();
	}
	gtk_widget_set_sensitive(lookup_widget("pref_script_automatic_updates"), com.pref.use_scripts_repository);
	gtk_widget_set_sensitive(lookup_widget("manual_script_sync_button"), com.pref.use_scripts_repository);
	gtk_widget_set_sensitive(lookup_widget("treeview2"), com.pref.use_scripts_repository);

}
