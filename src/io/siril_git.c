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
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_networking.h"
#include "core/siril_update.h" // for the version_number struct
#include "gui/message_dialog.h"
#include "gui/photometric_cc.h"
#include "gui/script_menu.h" // for SCRIPT_EXT TODO: after python3 is merged, move this out of src/gui
#include "gui/utils.h"
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
	// Open the repository
	git_repository *repo = NULL;
	int error = siril_repository_open(&repo, local_path);
	if (error != 0) {
		siril_log_color_message(
			_("Error performing hard reset. You may need to delete the local git "
			"repository and allow Siril to re-clone it.\n"), "red");
		git_repository_free(repo);
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
		return -1;
	}

	git_repository_free(repo);
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

static gboolean
version_meets_constraint(const version_number *current, const version_number *required, const gchar *op)
{
	int cmp = compare_version(*current, *required);

	if (g_strcmp0(op, "==") == 0) return cmp == 0;
	if (g_strcmp0(op, "!=") == 0) return cmp != 0;
	if (g_strcmp0(op, ">=") == 0) return cmp >= 0;
	if (g_strcmp0(op, "<=") == 0) return cmp <= 0;
	if (g_strcmp0(op, ">")  == 0) return cmp > 0;
	if (g_strcmp0(op, "<")  == 0) return cmp < 0;

	// Unknown operator
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

// gitscripts Repository synchronization and management
int sync_gitscripts_repository(gboolean sync) {
	int retval = 0;
	git_repository *repo = NULL;
	git_remote *remote = NULL;

	// URL of the remote repository
	siril_debug_print("Repository URL: %s\n", SCRIPT_REPOSITORY_URL);

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_scripts_repo_path();
	const char *remote_name = "origin";

	// Clone options
	git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;

	// Check if directory exists
	if (g_file_test(local_path, G_FILE_TEST_IS_DIR)) {
		siril_debug_print("Directory exists, attempting to open repository...\n");

		// Try to open existing repository
		int error = siril_repository_open(&repo, local_path);
		if (error != 0) {
			gui.script_repo_available = FALSE;
			siril_log_color_message(_("Failed to open repository.\n"), "red");
		} else {
			// Check we are using the correct repository
			error = git_remote_lookup(&remote, repo, remote_name);
			if (error != 0) {
				gui.script_repo_available = FALSE;
				siril_log_color_message(_("Failed to lookup remote.\n"), "red");
			} else {
				const char *remote_url = git_remote_url(remote);
				if (remote_url == NULL) {
					gui.script_repo_available = FALSE;
					error = 1;
					siril_log_color_message(
						_("Cannot identify local repository's configured remote.\n"), "red");
				} else {
					siril_debug_print("Remote URL: %s\n", remote_url);
					if (strcmp(remote_url, SCRIPT_REPOSITORY_URL)) {
						error = 1;
						gui.script_repo_available = FALSE;
						gchar *msg = g_strdup_printf(
							_("Local siril-scripts repository folder %s is not configured with %s as its remote.\n"),
							local_path, SCRIPT_REPOSITORY_URL);
						siril_log_color_message(msg, "red");
						g_free(msg);
					}
				}
			}
		}

		if (error == 0) {
			// Repository opened successfully
			siril_debug_print("Local scripts repository opened successfully...\n");
			gui.script_repo_available = TRUE;
		} else {
			// Failed to open repository - directory exists but isn't a valid repo
			const git_error *e = giterr_last();
			siril_debug_print("Cannot open repository: %s\n", e ? e->message : "");
			siril_log_message(_("Existing directory is not a valid git repository or is "
					"misconfigured. Cleaning up...\n"));

			// Delete the existing directory
			GError *delete_error = NULL;
			if (!delete_directory(local_path, &delete_error)) {
				siril_log_color_message(_("Error deleting existing directory %s: %s\n"),
					"red", local_path, delete_error->message);
				g_error_free(delete_error);
				gui.script_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			}

			// Now attempt to clone into the clean directory
			if (is_online()) {
				siril_log_message(_("Attempting to clone from remote source...\n"));
				error = git_clone(&repo, SCRIPT_REPOSITORY_URL, local_path, &clone_opts);
				if (error != 0) {
					e = giterr_last();
					siril_log_color_message(_("Error cloning repository into %s: %s\n"),
						"red", local_path, e->message);
					gui.script_repo_available = FALSE;
					retval = 1;
					goto cleanup;
				} else {
					siril_log_message(_("Repository cloned successfully!\n"));
					gui.script_repo_available = TRUE;
				}
			} else {
				siril_log_message(_("Siril is in offline mode. Will not attempt to clone remote repository.\n"));
				gui.script_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			}
		}
	} else {
		// Directory doesn't exist - clone it
		siril_debug_print("Directory doesn't exist, attempting to clone...\n");

		if (is_online()) {
			siril_log_message(_("Attempting to clone from remote source...\n"));
			int error = git_clone(&repo, SCRIPT_REPOSITORY_URL, local_path, &clone_opts);
			if (error != 0) {
				const git_error *e = giterr_last();
				siril_log_color_message(_("Error cloning repository into %s: %s\n"),
					"red", local_path, e->message);
				gui.script_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			} else {
				siril_log_message(_("Repository cloned successfully!\n"));
				gui.script_repo_available = TRUE;
			}
		} else {
			siril_log_message(_("Siril is in offline mode. Will not attempt to clone remote repository.\n"));
			gui.script_repo_available = FALSE;
			retval = 1;
			goto cleanup;
		}
	}

	// Final verification of repository URL
	int error = git_remote_lookup(&remote, repo, remote_name);
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
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Repository Error"), msg);
		g_free(msg);
		// Make scripts unavailable, as the contents of a random git repository
		// could be complete rubbish
		gui.script_repo_available = FALSE;
		retval = 1;
		goto cleanup;
	}

	// Synchronise the repository if requested
	if (sync) {
		gui_repo_scripts_mutex_lock();

		// fetch, analyse and merge changes from the remote
		int fetch_val = lg2_fetch(repo);

		// Get the FETCH_HEAD reference
		git_object *target_commit = NULL;
		error = git_revparse_single(&target_commit, repo, "FETCH_HEAD");
		if (error != 0) {
			gui_repo_scripts_mutex_unlock();
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
			gui_repo_scripts_mutex_unlock();
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

		gui_repo_scripts_mutex_unlock();
	}

cleanup:
	if (remote) {
		git_remote_free(remote);
	}
	if (repo) {
		git_repository_free(repo);
	}

	return retval;
}

// Update GUI script list from repository
int update_repo_scripts_list() {
	int retval = 0;
	git_repository *repo = NULL;
	git_index *index = NULL;

	const gchar *local_path = siril_get_scripts_repo_path();

	// Check if repository directory exists
	if (!g_file_test(local_path, G_FILE_TEST_IS_DIR)) {
		siril_log_color_message(_("Scripts repository directory does not exist.\n"), "red");
		gui.script_repo_available = FALSE;
		retval = 1;
		goto cleanup;
	}

	// Try to open existing repository
	int error = siril_repository_open(&repo, local_path);
	if (error != 0) {
		siril_log_color_message(_("Failed to open scripts repository.\n"), "red");
		gui.script_repo_available = FALSE;
		retval = 1;
		goto cleanup;
	}

	// Get repository index
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

cleanup:
	if (index) {
		git_index_free(index);
	}
	if (repo) {
		git_repository_free(repo);
	}

	return retval;
}

// Called at the end of setting up the venv
gpointer update_repo_scripts_list_and_menu_in_thread() {
	update_repo_scripts_list();
	gui_mutex_lock();
	execute_idle_and_wait_for_it(refresh_script_menu_idle, NULL);
	gui_mutex_unlock();

	return GINT_TO_POINTER(0);
}

int auto_update_gitscripts(gboolean sync) {
	int retval = 0;

	// First, sync the repository
	retval = sync_gitscripts_repository(sync);
	if (retval != 0) {
		return retval;
	}

	// Only update script list if repository sync was successful
	if (gui.script_repo_available) {
		retval = update_repo_scripts_list();
	}

	return retval;
}

int auto_update_gitspcc(gboolean sync) {
	int retval = 0;
	git_repository *repo = NULL;
	git_remote *remote = NULL;

	// URL of the remote repository
	siril_debug_print("Repository URL: %s\n", SPCC_REPOSITORY_URL);

	// Local directory where the repository will be cloned
	const gchar *local_path = siril_get_spcc_repo_path();
	const char *remote_name = "origin";

	// Clone options
	git_clone_options clone_opts = GIT_CLONE_OPTIONS_INIT;
	// Check if directory exists
	if (g_file_test(local_path, G_FILE_TEST_IS_DIR)) {
		siril_debug_print("Directory exists, attempting to open repository...\n");

		// Try to open existing repository
		int error = siril_repository_open(&repo, local_path);
		if (error != 0) {
			gui.spcc_repo_available = FALSE;
			siril_log_color_message(_("Failed to open repository.\n"), "red");
		} else {
			// Check we are using the correct repository
			error = git_remote_lookup(&remote, repo, remote_name);
			if (error != 0) {
				gui.spcc_repo_available = FALSE;
				siril_log_color_message(_("Failed to lookup remote.\n"), "red");
			} else {

				const char *remote_url = git_remote_url(remote);
				if (remote_url == NULL) {
					gui.spcc_repo_available = FALSE;
					error = 1;
					siril_log_color_message(
						_("Cannot identify local repository's configured remote.\n"), "red");
				} else {
					siril_debug_print("Remote URL: %s\n", remote_url);
					if (strcmp(remote_url, SPCC_REPOSITORY_URL)) {
						error = 1;
						gui.spcc_repo_available = FALSE;
						gchar *msg = g_strdup_printf(
							_("Local siril-spcc-database repository folder %s is not configured with %s as its remote.\n"),
							local_path, SPCC_REPOSITORY_URL);
						siril_log_color_message(msg, "red");
						g_free(msg);
					}
				}
			}
		}

		if (error == 0) {
			// Repository opened successfully
			siril_debug_print("Local SPCC database repository opened successfully!\n");
			gui.spcc_repo_available = TRUE;
		} else {
			// Failed to open repository - directory exists but isn't a valid repo
			const git_error *e = giterr_last();
			siril_debug_print("Cannot open repository: %s\n", e ? e->message : "");
			siril_log_message(_("Existing directory is not a valid git repository. Cleaning up...\n"));

			// Delete the existing directory
			GError *delete_error = NULL;
			if (!delete_directory(local_path, &delete_error)) {
				siril_log_color_message(_("Error deleting existing directory %s: %s\n"),
					"red", local_path, delete_error->message);
				g_error_free(delete_error);
				gui.spcc_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			}

			// Now attempt to clone into the clean directory
			siril_log_message(_("Attempting to clone from remote source...\n"));
			error = git_clone(&repo, SPCC_REPOSITORY_URL, local_path, &clone_opts);
			if (error != 0) {
				e = giterr_last();
				siril_log_color_message(_("Error cloning repository: %s\n"), "red", e->message);
				gui.spcc_repo_available = FALSE;
				retval = 1;
				goto cleanup;
			} else {
				siril_log_message(_("Repository cloned successfully!\n"));
				gui.spcc_repo_available = TRUE;
			}
		}
	} else {
		// Directory doesn't exist - clone it
		siril_debug_print("Directory doesn't exist, attempting to clone...\n");
		siril_log_message(_("Attempting to clone from remote source...\n"));

		int error = git_clone(&repo, SPCC_REPOSITORY_URL, local_path, &clone_opts);
		if (error != 0) {
			const git_error *e = giterr_last();
			siril_log_color_message(_("Error cloning repository: %s\n"), "red", e->message);
			gui.spcc_repo_available = FALSE;
			retval = 1;
			goto cleanup;
		} else {
			siril_log_message(_("Repository cloned successfully!\n"));
			gui.spcc_repo_available = TRUE;
		}
	}

	// Check we are using the correct repository
	int error = git_remote_lookup(&remote, repo, remote_name);
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
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Repository Error"), msg);
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
	return retval;
}

/**
* Checks if a file was modified between two commits
*
* @param repo Pointer to the git repository
* @param filepath Path to the file relative to repository root
* @param commit1 First commit to compare
* @param commit2 Second commit to compare
* @return 1 if file was modified, 0 if not modified, -1 on error
*/
static int file_modified_between_commits(git_repository *repo,
									const char *filepath,
									git_commit *commit1,
									git_commit *commit2) {
	git_tree *tree1 = NULL, *tree2 = NULL;
	git_tree_entry *entry1 = NULL, *entry2 = NULL;
	int error = 0;
	int result = 0;

	// Get trees from both commits
	error = git_commit_tree(&tree1, commit1);
	if (error != 0) goto cleanup;

	error = git_commit_tree(&tree2, commit2);
	if (error != 0) goto cleanup;

	// Try to find the file in both trees
	int found1 = git_tree_entry_bypath(&entry1, tree1, filepath) == 0;
	int found2 = git_tree_entry_bypath(&entry2, tree2, filepath) == 0;

	// If file exists in one but not the other, it was modified (added/deleted)
	if (found1 != found2) {
		result = 1;
		goto cleanup;
	}

	// If file doesn't exist in either, no modification
	if (!found1 && !found2) {
		result = 0;
		goto cleanup;
	}

	// Both exist, compare their OIDs
	const git_oid *oid1 = git_tree_entry_id(entry1);
	const git_oid *oid2 = git_tree_entry_id(entry2);

	result = git_oid_equal(oid1, oid2) ? 0 : 1;

cleanup:
	if (entry1) git_tree_entry_free(entry1);
	if (entry2) git_tree_entry_free(entry2);
	if (tree1) git_tree_free(tree1);
	if (tree2) git_tree_free(tree2);

	if (error != 0) return -1;
	return result;
}

static int find_file_commit_by_modifications(git_repository *repo,
											const char *filepath,
											int file_revisions_to_backtrack,
											git_commit **out_commit,
											gchar **out_relative_path,
											gchar **commit_msg) {
	git_object *head_commit_obj = NULL;
	git_commit *current_commit = NULL;
	git_commit *parent_commit = NULL;
	git_commit *prev_commit = NULL;
	int error = 0;
	int modifications_found = 0;
	gchar *relative_path = NULL;

	if (!repo || !filepath || !out_commit || !out_relative_path || file_revisions_to_backtrack < 0)
		return -1;

	const char *workdir = git_repository_workdir(repo);
	const char *tmpworkdir = g_canonicalize_filename(workdir, NULL);
	if (tmpworkdir && g_path_is_absolute(filepath)) {
		size_t workdir_len = strlen(tmpworkdir);
		if (strncmp(filepath, tmpworkdir, workdir_len) == 0) {
			const char *rel_start = filepath + workdir_len;
			while (*rel_start == '/' || *rel_start == '\\') rel_start++;
			const gchar *tmprelpath =  posix_path_separators(rel_start);
			relative_path = g_strdup(tmprelpath);
		} else {
			return -1;
		}
	} else {
		relative_path = g_strdup(filepath);
	}

	error = git_revparse_single(&head_commit_obj, repo, "HEAD");
	if (error != 0) {
		g_free(relative_path);
		return error;
	}

	current_commit = (git_commit *)head_commit_obj;

	if (file_revisions_to_backtrack == 0) {
		git_commit_dup(out_commit, current_commit);
		*out_relative_path = relative_path;
		return 0;
	}

	prev_commit = current_commit;

	while (modifications_found < file_revisions_to_backtrack) {
		error = git_commit_parent(&parent_commit, prev_commit, 0);
		if (error != 0) {
			git_object_free(head_commit_obj);
			g_free(relative_path);
			return -1;
		}

		int modified = file_modified_between_commits(repo, relative_path, parent_commit, prev_commit);
		if (modified < 0) {
			git_commit_free(parent_commit);
			git_object_free(head_commit_obj);
			g_free(relative_path);
			return -1;
		}

		if (modified) {
			modifications_found++;
			if (modifications_found == file_revisions_to_backtrack) {
				*out_commit = parent_commit;
				*out_relative_path = relative_path;
				if (commit_msg) {
					const gchar *msg = git_commit_message(prev_commit);
					*commit_msg = g_strdup(msg);
				}
				git_object_free(head_commit_obj);
				if (prev_commit != current_commit)
					git_commit_free(prev_commit);
				return 0;
			}
		}

		if (prev_commit != current_commit)
			git_commit_free(prev_commit);
		prev_commit = parent_commit;
		parent_commit = NULL;
	}

	git_commit_free(prev_commit);
	git_object_free(head_commit_obj);
	g_free(relative_path);
	return -1;
}

static int get_file_content_from_file_revision(git_repository *repo,
											const char *filepath,
											int file_revisions_to_backtrack,
											gchar **content,
											size_t *content_size) {
	git_commit *target_commit = NULL;
	git_tree *tree = NULL;
	git_tree_entry *entry = NULL;
	git_blob *blob = NULL;
	int error = 0;
	gchar *relative_path = NULL;

	if (!repo || !filepath || !content || !content_size || file_revisions_to_backtrack < 0)
		return -1;

	*content = NULL;
	*content_size = 0;

	error = find_file_commit_by_modifications(repo, filepath, file_revisions_to_backtrack,
											&target_commit, &relative_path, NULL);
	if (error != 0)
		return -1;

	error = git_commit_tree(&tree, target_commit);
	if (error != 0)
		goto cleanup;

	error = git_tree_entry_bypath(&entry, tree, relative_path);
	if (error != 0)
		goto cleanup;

	if (git_tree_entry_type(entry) != GIT_OBJECT_BLOB) {
		error = -1;
		goto cleanup;
	}

	error = git_blob_lookup(&blob, repo, git_tree_entry_id(entry));
	if (error != 0)
		goto cleanup;

	const void *blob_content = git_blob_rawcontent(blob);
	git_off_t blob_size = git_blob_rawsize(blob);

	if (blob_size < 0) {
		error = -1;
		goto cleanup;
	}

	*content = g_malloc(blob_size + 1);
	if (!*content) {
		error = -1;
		goto cleanup;
	}

	memcpy(*content, blob_content, blob_size);
	(*content)[blob_size] = '\0';
	*content_size = blob_size;

cleanup:
	if (relative_path)
		g_free(relative_path);
	if (blob)
		git_blob_free(blob);
	if (entry)
		git_tree_entry_free(entry);
	if (tree)
		git_tree_free(tree);
	if (target_commit)
		git_commit_free(target_commit);

	if (error != 0 && *content) {
		g_free(*content);
		*content = NULL;
		*content_size = 0;
	}

	return error;
}

static int get_commit_from_file_revision(git_repository *repo,
										const char *filepath,
										int file_revisions_to_backtrack,
										gchar **message,
										size_t *message_size) {
	git_commit *target_commit = NULL;
	int error = 0;
	gchar *relative_path = NULL;

	if (!repo || !filepath || !message || !message_size || file_revisions_to_backtrack < 0)
		return -1;

	*message = NULL;
	*message_size = 0;

	gchar *msg = NULL;
	error = find_file_commit_by_modifications(repo, filepath, file_revisions_to_backtrack + 1,
											&target_commit, &relative_path, &msg);
	if (error != 0)
		return -1;

	if (!msg) {
		error = -1;
		goto cleanup;
	}

	size_t len = strlen(msg);
	*message = g_malloc(len + 1);
	if (!*message) {
		error = -1;
		goto cleanup;
	}

	memcpy(*message, msg, len);
	(*message)[len] = '\0';
	*message_size = len;

cleanup:
	if (relative_path)
		g_free(relative_path);
	if (target_commit)
		git_commit_free(target_commit);

	if (error != 0 && *message) {
		g_free(*message);
		*message = NULL;
		*message_size = 0;
	}

	return error;
}

/**
* Convenience wrapper that returns the file content as a null-terminated string,
* counting only commits that modified the specific file.
* Optionally returns the associated commit message.
* The caller is responsible for freeing the returned strings using g_free().
*
* @param filepath Path to the file relative to repository root
* @param file_revisions_to_backtrack Number of file modifications to go back
* @param content_size (out) Size of the returned content
* @param commit_message (out, optional) Pointer to receive commit message string
* @param message_size (out, optional) Size of the returned commit message
* @return Allocated string containing file content, or NULL on error.
*/
gchar *get_script_content_string_from_file_revision(const char *filepath,
													int file_revisions_to_backtrack,
													size_t *content_size,
													gchar **commit_message,
													size_t *message_size) {
	*content_size = 0;
	if (message_size) *message_size = 0;
	if (commit_message) *commit_message = NULL;

	git_repository *repo = NULL;
	gchar *content = NULL;
	gchar *message = NULL;

	// URL of the remote repository
	siril_debug_print("Repository URL: %s\n", SCRIPT_REPOSITORY_URL);

	// Local repository directory
	const gchar *local_path = siril_get_scripts_repo_path();

	// Check if directory exists
	gboolean success = FALSE;
	if (g_file_test(local_path, G_FILE_TEST_IS_DIR)) {
		siril_debug_print("Directory exists, attempting to open repository...\n");

		// Try to open existing repository
		int error = siril_repository_open(&repo, local_path);
		if (error == 0) {
			siril_debug_print("Scripts repository opened successfully!\n");
			success = TRUE;
		}
	}
	if (!success) {
		siril_log_color_message(_("Error opening scripts repository\n"), "red");
		goto cleanup;
	}

	int error = get_file_content_from_file_revision(repo, filepath, file_revisions_to_backtrack,
													&content, content_size);
	if (error != 0) {
		siril_debug_print("Error retrieving file content in get_script_content_string_from_file_revision()\n");
		goto cleanup;
	}

	if (commit_message) {
		error = get_commit_from_file_revision(repo, filepath, file_revisions_to_backtrack,
											&message, message_size);
		if (error != 0) {
			siril_debug_print("Error retrieving commit message in get_script_content_string_from_file_revision()\n");
			g_free(content);
			content = NULL;
			*content_size = 0;
		} else {
			*commit_message = message;
		}
	}

cleanup:
	if (repo)
		git_repository_free(repo);
	return content;
}


#endif // HAVE_LIBGIT2
