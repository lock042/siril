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

#include "algos/spcc.h"
#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_networking.h"
#include "core/siril_update.h" // for the version_number struct
#include "gui/message_dialog.h"
#include "gui/photometric_cc.h"
#include "gui/script_menu.h" // for SCRIPT_EXT TODO: after python3 is merged, move this out of src/gui
#include "io/siril_git.h"
#include <assert.h>
#include <inttypes.h>

#define REPO_REPAIRED 999

// #define DEBUG_GITSCRIPTS

#ifdef HAVE_LIBGIT2
#include <git2.h>

static GMutex gui_repo_scripts_mutex = { 0 };

void gui_repo_scripts_mutex_lock() {
	g_mutex_lock(&gui_repo_scripts_mutex);
}

void gui_repo_scripts_mutex_unlock() {
	g_mutex_unlock(&gui_repo_scripts_mutex);
}

const gchar *SCRIPT_REPOSITORY_URL = "https://gitlab.com/free-astro/siril-scripts";
const gchar *SPCC_REPOSITORY_URL = "https://gitlab.com/free-astro/siril-spcc-database";

static void *xrealloc(void *oldp, size_t newsz) {
	void *p = realloc(oldp, newsz);
	if (p == NULL) {
		PRINT_ALLOC_ERR;
		// exit(1);
	}
	return p;
}

enum {
	FORMAT_DEFAULT = 0,
	FORMAT_LONG = 1,
	FORMAT_SHORT = 2,
	FORMAT_PORCELAIN = 3
};

struct merge_options {
	const char **heads;
	size_t heads_count;

	git_annotated_commit **annotated;
	size_t annotated_count;
};

static void merge_options_init(struct merge_options *opts) {
	memset(opts, 0, sizeof(*opts));

	opts->heads = NULL;
	opts->heads_count = 0;
	opts->annotated = NULL;
	opts->annotated_count = 0;
}

static void opts_add_refish(struct merge_options *opts, const char *refish) {
	size_t sz;
	assert(opts != NULL);

	sz = ++opts->heads_count * sizeof(opts->heads[0]);
	opts->heads = xrealloc((void *)opts->heads, sz);
	opts->heads[opts->heads_count - 1] = refish;
}

static int resolve_refish(git_annotated_commit **commit, git_repository *repo,
                          const char *refish) {
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

static int resolve_heads(git_repository *repo, struct merge_options *opts) {
	git_annotated_commit **annotated =
		calloc(opts->heads_count, sizeof(git_annotated_commit *));
	size_t annotated_count = 0, i;

	for (i = 0; i < opts->heads_count; i++) {
		if (!opts->heads[i]) {
			annotated_count--;
			continue;
		}
		int err =
			resolve_refish(&annotated[annotated_count++], repo, opts->heads[i]);
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

static char *get_commit_from_oid(git_repository *repo, git_oid *oid_to_find) {
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

static int fetchhead_cb(const char *ref_name, const char *remote_url,
                        const git_oid *oid, unsigned int is_merge,
                        void *payload) {
	if (is_merge) {
		siril_debug_print("reference: '%s' is the reference we should merge\n", ref_name);
		git_oid_cpy((git_oid *)payload, oid);
	}
	return 0;
}

/**
 * Removes Git lock files from a repository directory path.
 * This is a path-based version that doesn't require an open repository.
 *
 * @param git_dir Path to the Git repository directory (typically the .git directory)
 * @return 0 on success, or the error code on failure
 */
static int remove_git_locks_by_path(const char *git_dir) {
    GError *error = NULL;
    gint ret = 0;

    // Remove index.lock
    gchar *index_lock_path = g_build_filename(git_dir, "index.lock", NULL);
    if (g_unlink(index_lock_path) != 0 && errno != ENOENT) {
        g_warning("Error removing index.lock: %s", g_strerror(errno));
        ret = errno;
    }
    g_free(index_lock_path);

    // Remove branch lock files
    GDir *dir = g_dir_open(git_dir, 0, &error);
    if (dir) {
        const gchar *entry;
        while ((entry = g_dir_read_name(dir)) != NULL) {
            if (g_str_has_suffix(entry, ".lock")) {
                gchar *lock_path = g_build_filename(git_dir, entry, NULL);
                if (g_unlink(lock_path) != 0 && errno != ENOENT) {
                    g_warning("Error removing lock file '%s': %s", entry, g_strerror(errno));
                    ret = errno;
                }
                g_free(lock_path);
            }
        }
        g_dir_close(dir);
    } else {
        g_warning("Error opening Git directory '%s': %s", git_dir, error->message);
        ret = error->code;
        g_clear_error(&error);
    }

    return ret;
}

/**
 * Attempts to open a Git repository. If that fails due to locks,
 * removes lock files and tries again.
 *
 * @param out Pointer to store the opened repository
 * @param path Path to the Git repository
 * @return 0 on success, error code on failure
 * WARNING: a failure here may be due to the path not existing so you cannot
 * rely on git_error_last being set. Check it is non-NULL before printing e->message!!
 */
static int siril_repository_open(git_repository **out, const gchar *path) {
	// Check if directory exists first
	if (!g_file_test(path, G_FILE_TEST_IS_DIR)) {
		// Directory doesn't exist, return with error
		return -1;
	}
	int retval = git_repository_open(out, path);

	if (retval != 0) {  // Opening failed, try to fix by removing locks
		// Construct the path to the .git directory
		gchar *git_dir = NULL;

		// Check if path itself is a .git directory or a working directory
		if (g_str_has_suffix(path, ".git") || g_file_test(g_build_filename(path, ".git", NULL), G_FILE_TEST_IS_DIR)) {
			// Path is either the .git directory or contains a .git directory
			if (g_str_has_suffix(path, ".git")) {
				// Path is the .git directory itself
				git_dir = g_strdup(path);
			} else {
				// Path is the working directory, need to append .git
				git_dir = g_build_filename(path, ".git", NULL);
			}

			// Remove any existing lock files
			int error = remove_git_locks_by_path(git_dir);
			g_free(git_dir);

			if (error != 0) {
				siril_log_color_message(_("Error removing Git lock files. You may need to delete the local "
										"git repository and allow Siril to re-clone it.\n"), "red");
				return -1;
			}

			// Try opening again after removing locks
			retval = git_repository_open(out, path);
		} else {
			// Not a valid git path
			siril_log_color_message(_("Invalid Git repository path.\n"), "red");
		}
	}

	return retval;
}

int reset_repository(const gchar *local_path) {
	// Initialisation
	git_libgit2_init();

	// Open the repository
	git_repository *repo = NULL;
	int error = siril_repository_open(&repo, local_path);
	if (error != 0) {
		siril_log_color_message(
			_("Error performing hard reset. You may need to delete the local git "
			"repository and allow Siril to re-clone it.\n"), "red");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	// Get the FETCH_HEAD reference
	git_object *target_commit = NULL;
	error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
	if (error != 0) {
		siril_log_color_message(
			_("Error performing hard reset. You may need to delete the local git "
			"repository and allow Siril to re-clone it.\n"),
			"red");
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	// Perform the reset
	error = git_reset(repo, target_commit, GIT_RESET_HARD, NULL);
	if (error != 0) {
		siril_log_color_message(
			_("Error performing hard reset. You may need to delete the local git "
			"repository and allow Siril to re-clone it.\n"), "red");
		git_object_free(target_commit);
		git_repository_free(repo);
		git_libgit2_shutdown();
		return -1;
	}

	git_repository_free(repo);
	git_libgit2_shutdown();
	return 0;
}

static int lg2_fetch(git_repository *repo) {
	git_remote *remote = NULL;
	git_fetch_options fetch_opts = GIT_FETCH_OPTIONS_INIT;
	const char *remote_name = "origin";

	siril_debug_print("Fetching %s for repo %p\n", remote_name, repo);

	int val;
	if (git_remote_lookup(&remote, repo, remote_name)) {
		if (git_remote_create_anonymous(&remote, repo, remote_name)) {
			siril_log_message(_("Unable to create anonymous remote for the repository. Check your connectivity and try again.\n"));
			val = 1;
			goto on_error;
		}
	}
	val = git_remote_fetch(remote, NULL, &fetch_opts, "fetch");
	if (val < 0) {
		if (val == GIT_ENOTFOUND) {
			// Test if we can connect to the remote to rule out network issues
			git_remote_callbacks callbacks = GIT_REMOTE_CALLBACKS_INIT;
			int can_connect = 0;

			// Try to connect to the remote
			if (git_remote_connect(remote, GIT_DIRECTION_FETCH, &callbacks, NULL, NULL) == 0) {
				can_connect = 1;
				git_remote_disconnect(remote);
			}
			// If we can connect, assume this is probably the 1.4.0-beta2 bug
			if (can_connect) {
				siril_log_message(_("Detected 1.4.0-beta2 bug condition. Removing and re-cloning repository...\n"));
				const char *repo_root = git_repository_workdir(repo);
				GError *gerror = NULL;
				delete_directory(repo_root, &gerror);
				if (gerror) {
					siril_log_color_message(_("Error removing repository directory: %s\n"), "red", gerror->message);
					g_clear_error(&gerror);
					goto on_error;
				}
				// Clone options
				git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;
				int value = git_clone(&repo, SCRIPT_REPOSITORY_URL, repo_root, &clone_opts);
				if (value != 0) {
					const git_error *e = giterr_last();
					siril_log_color_message(_("Error cloning repository into %s: %s\n"), "red", repo_root, e->message);
					gui.script_repo_available = FALSE;
					goto on_error;
				} else {
					siril_log_message(_("Repository re-cloned successfully!\n"));
					val = REPO_REPAIRED;
				}
			} else {
				siril_log_message(_("Error fetching remote: Remote repository or references not found. "
								"Check the repository URL and your connectivity.\n"));
			}
		} else {
			// Generic fallback
			siril_log_message(_("Error fetching remote. This may be a temporary error, check your connectivity and try again.\n"));
		}
	}

on_error:
	if (remote) {
		git_remote_free(remote);
	}
	return val;
}

static char *find_str_before_comment(const char *str1, const char *str2, const char *str3) {
	char *strpos = strstr(str1, str2);
	const char *chrpos = strstr(str1, str3);
	return !strpos ? NULL : (chrpos && chrpos < strpos) ? NULL : strpos;
}

static gboolean script_version_check(const gchar *filename) {
	version_number current_version = get_current_version_number();
	// Open the script and look for the required version number
	GFile *file = NULL;
	GInputStream *stream = NULL;
	GDataInputStream *data_input = NULL;
	GError *error = NULL;
	gchar *buffer = NULL;
	gsize length = 0;
	gchar *scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), filename, NULL);
	gboolean retval = FALSE;
	#ifdef DEBUG_GITSCRIPTS
	printf("checking script version requirements: %s\n", scriptpath);
	#endif
	file = g_file_new_for_path(scriptpath);
	stream = (GInputStream *)g_file_read(file, NULL, &error);
	if (error)
		goto ERROR_OR_COMPLETE;
	data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, &error)) && !error) {
		gchar *ver = find_str_before_comment(buffer, "requires", "#");
		if (ver) {
			gchar **versions = g_strsplit(ver, " ", 2);
			version_number required_version =
				get_version_number_from_string(versions[1]);
			version_number expired_version = {0};
			int new_enough = compare_version(current_version, required_version);
			int too_new = -1;
			if (versions[2]) {
				expired_version = get_version_number_from_string(versions[2]);
				too_new = compare_version(current_version, expired_version);
			}
			retval = ((new_enough >= 0) && (too_new < 0));
			g_strfreev(versions);
			g_free(buffer);
			buffer = NULL;
			break;
		}
		g_free(buffer);
	}
	ERROR_OR_COMPLETE:
	g_input_stream_close(stream, NULL, &error);
	if (error)
		siril_debug_print("Error closing data input stream from file\n");
	g_free(scriptpath);
	g_object_unref(data_input);
	g_object_unref(stream);
	g_object_unref(file);
	return retval;
}

static gboolean version_meets_constraint(version_number *current, version_number *constraint, const gchar *op) {
	if (strcmp(op, "==") == 0)
		return compare_version(*current, *constraint);
	if (strcmp(op, "!=") == 0)
		return !compare_version(*current, *constraint);
	if (strcmp(op, "<") == 0) {
		if (current->major_version < constraint->major_version) return TRUE;
		if (current->major_version > constraint->major_version) return FALSE;
		if (current->minor_version < constraint->minor_version) return TRUE;
		if (current->minor_version > constraint->minor_version) return FALSE;
		if (current->micro_version < constraint->micro_version) return TRUE;
		if (current->micro_version > constraint->micro_version) return FALSE;
		if (current->patched_version < constraint->patched_version) return TRUE;
		return FALSE;
	}
	if (strcmp(op, ">") == 0) {
		if (current->major_version > constraint->major_version) return TRUE;
		if (current->major_version < constraint->major_version) return FALSE;
		if (current->minor_version > constraint->minor_version) return TRUE;
		if (current->minor_version < constraint->minor_version) return FALSE;
		if (current->micro_version > constraint->micro_version) return TRUE;
		if (current->micro_version < constraint->micro_version) return FALSE;
		if (current->patched_version > constraint->patched_version) return TRUE;
		return FALSE;
	}
	if (strcmp(op, "<=") == 0) {
		return version_meets_constraint(current, constraint, "<") ||
			version_meets_constraint(current, constraint, "==");
	}
	if (strcmp(op, ">=") == 0) {
		return version_meets_constraint(current, constraint, ">") ||
			version_meets_constraint(current, constraint, "==");
	}
	return FALSE;
}

static gboolean parse_version_string(const gchar *version_str, version_number *ver) {
	gchar **parts = g_strsplit(version_str, ".", 4);
	memset(ver, 0, sizeof(version_number));

	if (parts[0]) ver->major_version = atoi(parts[0]);
	if (parts[1]) ver->minor_version = atoi(parts[1]);
	if (parts[2]) ver->micro_version = atoi(parts[2]);
	if (parts[3]) ver->patched_version = atoi(parts[3]);

	g_strfreev(parts);
	return TRUE;
}

gboolean check_module_version_constraint(const gchar *line, GMatchInfo *match_info) {
	// Regular expression to match check_module_version pattern
	gboolean result = TRUE;

	gchar *version_spec = g_match_info_fetch(match_info, 1);

	// If no version specified, return TRUE
	if (!version_spec || strlen(version_spec) == 0) {
		result = TRUE;
		return result;
	}

	// Split the version specifier into constraint parts
	gchar **constraints = g_strsplit_set(version_spec, ",", -1);
	version_number minver = {0}, maxver = {0};
	gboolean min_set = FALSE, max_set = FALSE;

	for (int i = 0; constraints[i] != NULL; i++) {
		gchar *constraint = g_strstrip(constraints[i]);

		if (strncmp(constraint, "==", 2) == 0) {
			parse_version_string(constraint + 2, &minver);
			parse_version_string(constraint + 2, &maxver);
			min_set = max_set = TRUE;
		} else if (strncmp(constraint, "!=", 2) == 0) {
			parse_version_string(constraint + 2, &minver);
			result = !version_meets_constraint(&com.python_version, &minver, "==");
		} else if (strncmp(constraint, ">=", 2) == 0) {
			parse_version_string(constraint + 2, &minver);
			result = version_meets_constraint(&com.python_version, &minver, ">=");
			min_set = TRUE;
		} else if (strncmp(constraint, "<=", 2) == 0) {
			parse_version_string(constraint + 2, &maxver);
			result = version_meets_constraint(&com.python_version, &maxver, "<=");
			max_set = TRUE;
		} else if (strncmp(constraint, ">", 1) == 0) {
			parse_version_string(constraint + 1, &minver);
			result = version_meets_constraint(&com.python_version, &minver, ">");
			min_set = TRUE;
		} else if (strncmp(constraint, "<", 1) == 0) {
			parse_version_string(constraint + 1, &maxver);
			result = version_meets_constraint(&com.python_version, &maxver, "<");
			max_set = TRUE;
		}
	}

	// If multiple constraints are specified
	if (min_set && max_set) {
		result = version_meets_constraint(&com.python_version, &minver, ">=") &&
					version_meets_constraint(&com.python_version, &maxver, "<");
	}

	g_strfreev(constraints);
	g_free(version_spec);

	return result;
}

static gboolean pyscript_version_check(const gchar *filename) {
	// Open the script and look for the required version number
	GFile *file = NULL;
	GInputStream *stream = NULL;
	GDataInputStream *data_input = NULL;
	GError *error = NULL;
	gchar *buffer = NULL;
	gsize length = 0;
	gchar *scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), filename, NULL);
	gboolean retval = TRUE; // default to TRUE - if check_module_version() is not called there are no version requirements
	#ifdef DEBUG_GITSCRIPTS
	printf("checking python script version requirements: %s\n", scriptpath);
	#endif
	file = g_file_new_for_path(scriptpath);
	stream = (GInputStream *)g_file_read(file, NULL, &error);
	if (error)
		goto ERROR_OR_COMPLETE;
	data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, &error)) && !error) {
		GRegex *regex;
		GMatchInfo *match_info;
		gboolean matched = FALSE;
		// Regex to capture version checks with optional version specifier
		regex = g_regex_new("check_module_version\\(\"([^\"]+)?\"\\)", 0, 0, NULL);
		if (g_regex_match(regex, buffer, 0, &match_info)) {
			retval = check_module_version_constraint(buffer, match_info);
			matched = TRUE;
		}
		g_match_info_free(match_info);
		g_regex_unref(regex);
		g_free(buffer);
		if (matched)
			break;
	}
	ERROR_OR_COMPLETE:
	g_input_stream_close(stream, NULL, &error);
	if (error)
		siril_debug_print("Error closing data input stream from file\n");
	g_free(scriptpath);
	g_object_unref(data_input);
	g_object_unref(stream);
	g_object_unref(file);
	return retval;
}

static int analyse(git_repository *repo, GString **git_pending_commit_buffer) {
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
			return -1;
		}
	}

	error = resolve_heads(repo, &opts);
	if (error != 0) {
		free((char **)opts.heads);
		free(opts.annotated);

		return 2;
	}

	error = git_merge_analysis(&analysis, &preference, repo, (const git_annotated_commit **)opts.annotated, opts.annotated_count);

	// If the merge cannot be fast-forwarded, warn the user that local changes
	// will be lost if they proceed.
	if (error < 0) {
		siril_debug_print("Error carrying out merge analysis: %s\n", giterr_last()->message);
		return 2;
	}
	gboolean can_fastforward = FALSE;
	if ((analysis & GIT_MERGE_ANALYSIS_FASTFORWARD) || (analysis & GIT_MERGE_ANALYSIS_UP_TO_DATE)) {
		can_fastforward = TRUE;
	}

	// If we already know we can't fast forward we can skip the rest of the
	// function, we just return 2
	if (!can_fastforward) {
		siril_debug_print("Cannot be fast forwarded\n");
		git_remote_free(remote);
		return 2;
	}

	// Prepare for looping through unmerged commit messages

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
		// Check if the current commit is the HEAD. If not, we print the commit
		// message. If so, we break and don't show any further messages.
		if (git_oid_equal(git_commit_id(head_commit), git_commit_id(commit))) {
			found_head_ancestor = TRUE;
			break;
		}

		// We have not yet reached the HEAD so we print the commit message.
		commit_msg = git_commit_message(commit);
		if (!*git_pending_commit_buffer) {
			gchar *buf = g_strdup_printf(_("Commit message: %s\n"), commit_msg);
			*git_pending_commit_buffer = g_string_new(buf);
			g_free(buf);
		} else {
		g_string_append_printf(*git_pending_commit_buffer, _("Commit message: %s\n"), commit_msg);
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

	// If there is no ancestor commit found it indicates the merge cannot be
	// fast-forwarded
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
	git_repository *repo = NULL;
	git_remote *remote = NULL;
	git_index *index = NULL;

	// Initialize libgit2
	git_libgit2_init();

	// URL of the remote repository
	siril_debug_print("Repository URL: %s\n", SCRIPT_REPOSITORY_URL);

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();

	// Clone options
	git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;

	// See if the repository already exists
	int error = siril_repository_open(&repo, local_path);

	if (error != 0) {
		const git_error *e = giterr_last();
		siril_debug_print("Cannot open repository: %s\n", e ? e->message : "");
		if (is_online()) {
			siril_log_message(_("Attempting to clone from remote source...\n"));
			// Perform the clone operation
			error = git_clone(&repo, SCRIPT_REPOSITORY_URL, local_path, &clone_opts);

			if (error != 0) {
				e = giterr_last();
				siril_log_color_message(_("Error cloning repository into %s: %s\n"), "red", local_path, e->message);
				gui.script_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			} else {
				siril_log_message(_("Repository cloned successfully!\n"));
			}
		} else {
			siril_log_message(_("Siril is in offline mode. Will not attempt to clone remote repository.\n"));
			gui.script_repo_available = FALSE;
			retval = 1;
			goto cleanup;
		}
	} else {
		siril_debug_print("Local scripts repository opened successfully!\n");
	}
	gui.script_repo_available = TRUE;

	// Check we are using the correct repository
	const char *remote_name = "origin";
	error = git_remote_lookup(&remote, repo, remote_name);
	if (error != 0) {
		siril_log_color_message(_("Failed to lookup remote.\n"), "red");
		retval = 1;
		goto cleanup;
	}

	const char *remote_url = git_remote_url(remote);
	if (remote_url == NULL) {
		siril_log_color_message(
			_("Error: cannot identify local repository's configured remote.\n"), "red");
		retval = 1;
		goto cleanup;
	}

	siril_debug_print("Remote URL: %s\n", remote_url);

	if (strcmp(remote_url, SCRIPT_REPOSITORY_URL)) {
		gchar *msg = g_strdup_printf(
			_("Error: local siril-scripts repository folder is not "
			"configured with %s as its remote. You should remove the folder %s "
			"and restart Siril to re-clone the correct repository.\n"),
			SCRIPT_REPOSITORY_URL, local_path);
		siril_log_color_message(msg, "red");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Repository Error"), msg);
		g_free(msg);
		// Make scripts unavailable, as the contents of a random git repository
		// could be complete rubbish
		gui.script_repo_available = FALSE;
		retval = 1;
		goto cleanup;
	}

	// Synchronise the repository
	if (sync) {
		// fetch, analyse and merge changes from the remote
		int fetch_val = lg2_fetch(repo);

		// Get the FETCH_HEAD reference
		git_object *target_commit = NULL;
		error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem "
					"persists you may need to delete the local git repository and "
					"allow Siril to re-clone it.\n"), "red");
			gui.script_repo_available = FALSE;
			retval = -1;
			goto cleanup;
		}

		// Perform the reset
		error = git_reset(repo, target_commit, GIT_RESET_HARD, NULL);
		git_object_free(target_commit);
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem persists "
					"you may need to delete the local git repository and allow Siril to "
					"re-clone it.\n"), "red");
			retval = -1;
			goto cleanup;
		}

		// Print the "up-to-date" status message
		if (!fetch_val || fetch_val == REPO_REPAIRED) {
			siril_log_color_message(_("Local scripts repository is up-to-date!\n"), "green");
		}
	}

	/*** Populate the list of available repository scripts ***/
	error = git_repository_index(&index, repo);
	if (error < 0) {
		siril_log_color_message(_("Error accessing repository index.\n"), "red");
		retval = 1;
		goto cleanup;
	}

	gui_repo_scripts_mutex_lock();

	/* populate gui.repo_scripts with all the scripts in the index.
	* We ignore anything not ending in SCRIPT_EXT */
	size_t entry_count = git_index_entrycount(index);
	if (gui.repo_scripts) {
		g_slist_free_full(gui.repo_scripts, g_free);
		gui.repo_scripts = NULL;
	}

	for (size_t i = 0; i < entry_count; i++) {
		const git_index_entry *entry = git_index_get_byindex(index, i);
		if (entry == NULL) {
			siril_log_color_message(_("Warning: failed to get repository index entry.\n"), "red");
			continue;  // Skip this entry but continue processing
		}

		// Validate that the path is valid UTF-8
		if (!g_utf8_validate(entry->path, -1, NULL)) {
			siril_log_color_message(_("Warning: skipping script with invalid UTF-8 path.\n"), "red");
			continue;
		}

		if ((g_str_has_suffix(entry->path, SCRIPT_EXT) && script_version_check(entry->path)) ||
			(g_str_has_suffix(entry->path, PYSCRIPT_EXT) && pyscript_version_check(entry->path)) ||
			(g_str_has_suffix(entry->path, PYCSCRIPT_EXT))) {

			gchar *path_copy = g_strdup(entry->path);
			if (path_copy == NULL) {
				siril_log_color_message(_("Warning: memory allocation failed for script path.\n"), "red");
				continue;  // Skip this entry but continue processing
			}

			gui.repo_scripts = g_slist_prepend(gui.repo_scripts, path_copy);
#ifdef DEBUG_GITSCRIPTS
			printf("%s\n", entry->path);
#endif
		}
	}
	gui_repo_scripts_mutex_unlock();

	// Cleanup
cleanup:
	if (index) {
		git_index_free(index);
	}
	if (remote) {
		git_remote_free(remote);
	}
	if (repo) {
		git_repository_free(repo);
	}
	git_libgit2_shutdown();

	return retval;
}

int auto_update_gitspcc(gboolean sync) {
	int retval = 0;
	git_repository *repo = NULL;
	git_remote *remote = NULL;

	// Initialize libgit2
	git_libgit2_init();

	// URL of the remote repository
	siril_debug_print("Repository URL: %s\n", SPCC_REPOSITORY_URL);

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_spcc_repo_path();

	// Clone options
	git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;

	// See if the repository already exists
	int error = siril_repository_open(&repo, local_path);

	if (error != 0) {
		const git_error *e = giterr_last();
		siril_debug_print("Cannot open repository: %s\nAttempting to clone from remote source...\n", e ? e->message : "");
		// Perform the clone operation
		error = git_clone(&repo, SPCC_REPOSITORY_URL, local_path, &clone_opts);
		if (error != 0) {
			e = giterr_last();
			siril_log_color_message(_("Error cloning repository: %s\n"), "red", e->message);
			gui.spcc_repo_available = FALSE;
			retval = 1;
			goto cleanup;
		} else {
			siril_log_message(_("Repository cloned successfully!\n"));
		}
	} else {
		siril_debug_print("Local SPCC database repository opened successfully!\n");
	}
	gui.spcc_repo_available = TRUE;

	// Check we are using the correct repository
	const char *remote_name = "origin";
	error = git_remote_lookup(&remote, repo, remote_name);
	if (error != 0) {
		siril_log_color_message(_("Failed to lookup remote.\n"), "red");
		retval = 1;
		goto cleanup;
	}

	const char *remote_url = git_remote_url(remote);
	if (remote_url == NULL) {
		siril_log_color_message(
			_("Error: cannot identify local repository's configured remote.\n"), "red");
		retval = 1;
		goto cleanup;
	}

	siril_debug_print("Remote URL: %s\n", remote_url);

	if (strcmp(remote_url, SPCC_REPOSITORY_URL)) {
		gchar *msg = g_strdup_printf(
			_("Error: local siril-spcc-database repository folder is not "
			"configured with %s as its remote. You should remove the folder %s "
			"and restart Siril to re-clone the correct repository.\n"),
			SPCC_REPOSITORY_URL, local_path);
		siril_log_color_message(msg, "red");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Repository Error"), msg);
		g_free(msg);
		// Make scripts unavailable, as the contents of a random git repository
		// could be complete rubbish
		gui.spcc_repo_available = FALSE;
		retval = 1;
		goto cleanup;
	}

	// Synchronise the repository
	if (sync) {
		// fetch, analyse and merge changes from the remote
		lg2_fetch(repo);

		// Get the FETCH_HEAD reference
		git_object *target_commit = NULL;
		error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem "
					"persists you may need to delete the local git repository and "
					"allow Siril to re-clone it.\n"), "red");
			gui.spcc_repo_available = FALSE;
			retval = -1;
			goto cleanup;
		}

		// Perform the reset
		error = git_reset(repo, target_commit, GIT_RESET_HARD, NULL);
		git_object_free(target_commit);
		if (error != 0) {
			siril_log_color_message(_("Error performing hard reset. If the problem "
					"persists you may need to delete the local git repository and "
					"allow Siril to re-clone it.\n"), "red");
			retval = -1;
			goto cleanup;
		}
		siril_log_color_message(_("Local SPCC database repository is up-to-date!\n"), "green");
	}

	// Cleanup
cleanup:
	if (remote) {
		git_remote_free(remote);
	}
	if (repo) {
		git_repository_free(repo);
	}
	git_libgit2_shutdown();

	return retval;
}

int preview_scripts_update(GString** git_pending_commit_buffer) {
	// Initialize libgit2
	git_libgit2_init();

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();

	git_repository *repo = NULL;
	int error = siril_repository_open(&repo, local_path);
	if (error < 0) {
		siril_debug_print("Error opening repository: %s\n", giterr_last()->message);
		siril_log_color_message(_("Error: unable to open local scripts repository.\n"), "red");
		gui.script_repo_available = FALSE;
		return 1;
	}

	// Fetch changes
	lg2_fetch(repo);
	// Analyse the repository against the remote
	int retval = analyse(repo, git_pending_commit_buffer);
	git_repository_free(repo);
	git_libgit2_shutdown();
	return retval;
}

int preview_spcc_update(GString **git_pending_commit_buffer) {
	// Initialize libgit2
	git_libgit2_init();

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_spcc_repo_path();

	git_repository *repo = NULL;
	int error = siril_repository_open(&repo, local_path);
	if (error < 0) {
		siril_debug_print("Error opening repository: %s\n", giterr_last()->message);
		siril_log_color_message(_("Error: unable to open local spcc-database repository.\n"), "red");
		gui.spcc_repo_available = FALSE;
		return 1;
	}

	// Fetch changes
	lg2_fetch(repo);
	// Analyse the repository against the remote
	int retval = analyse(repo, git_pending_commit_buffer);
	git_repository_free(repo);
	git_libgit2_shutdown();
	return retval;
}

#endif // HAVE_LIBGIT2
