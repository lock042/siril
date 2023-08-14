#include <git2.h>
#include "core/siril.h"
#include "core/siril_log.h"

int update_gitscripts(void) {
    // Initialize libgit2
    git_libgit2_init();

    // URL of the remote repository
    const char *url = "https://gitlab.com/free-astro/siril-scripts.git";

    // Local directory where the repository will be cloned
    char *local_path = g_build_filename(g_get_user_cache_dir(), "siril-scripts", NULL);

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

	}

    // Cleanup
    git_repository_free(repo);
    git_libgit2_shutdown();

    return 0;
}
