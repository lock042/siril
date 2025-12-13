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
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "core/siril_update.h"

#if defined(HAVE_LIBCURL)
#include "yyjson.h"
#include "core/siril_networking.h"
#include "core/processing.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"

#define SIRIL_DOMAIN "https://siril.org/"
#define SIRIL_VERSIONS SIRIL_DOMAIN"siril_versions.json"
#define SIRIL_DOWNLOAD SIRIL_DOMAIN"download/"
#define GITLAB_URL "https://gitlab.com/free-astro/siril/raw"
#define BRANCH "master"
#define SIRIL_NOTIFICATIONS "notifications/siril_notifications.json"
#endif

// ============================================================================
// UTILITY FUNCTIONS - Independent of libcurl
// These functions work with version numbers and strings, no network required
// ============================================================================

static void remove_alpha(gchar *str, gboolean *is_rc, gboolean *is_beta) {
	unsigned long i = 0;
	unsigned long j = 0;
	char c;

	if (g_str_has_prefix(str, "beta")) {
		*is_rc = FALSE;
		*is_beta = TRUE;
	} else if (g_str_has_prefix(str, "rc")) {
		*is_rc = TRUE;
		*is_beta = FALSE;
	} else {
		*is_rc = FALSE;
		*is_beta = FALSE;
	}

	while ((c = str[i++]) != '\0') {
		if (g_ascii_isdigit(c)) {
			str[j++] = c;
		}
	}
	str[j] = '\0';
}

/**
 * Check if the version is a patched version.
 * patched version are named like that x.y.z.patch where patch only contains digits.
 * if patch contains alpha char it is because that's a RC or beta version. Not a patched one.
 * @param version version to be tested
 * @return 0 if the version is not patched. The version of the patch is returned otherwise.
 */
static guint check_for_patch(gchar *version, gboolean *is_rc, gboolean *is_beta) {
	remove_alpha(version, is_rc, is_beta);
	return (g_ascii_strtoull(version, NULL, 10));
}

version_number get_version_number_from_string(const gchar *input) {
	version_number version = { 0 };
	gchar **version_string = NULL;
	const gchar *string = find_first_numeric(input);
	if (!string)
		goto the_end;
	version_string = g_strsplit_set(string, ".-", -1);
	version.major_version = g_ascii_strtoull(version_string[0], NULL, 10);
	if (version_string[1])
		version.minor_version = g_ascii_strtoull(version_string[1], NULL, 10);
	else
		goto the_end;
	if (version_string[2])
		version.micro_version = g_ascii_strtoull(version_string[2], NULL, 10);
	else
		goto the_end;
	if (version_string[3] == NULL) {
		version.patched_version = 0;
		version.rc_version = FALSE;
		version.beta_version = FALSE;
	} else {
		version.patched_version = check_for_patch(version_string[3], &version.rc_version, &version.beta_version);
	}
the_end:
	g_strfreev(version_string);
	return version;
}

version_number get_current_version_number() {
	return get_version_number_from_string(PACKAGE_VERSION);
}

/**
 * This function compares two version numbers following the pattern x.y.z[-type#]
 * where type can be beta, rc, or a stable patch number.
 * Version ordering: beta < rc < stable < stable-patch
 * Examples: 1.4.0-beta1 < 1.4.0-rc2 < 1.4.0 < 1.4.0-1
 * @param v1 First version number to be tested
 * @param v2 Second version number to be tested
 * @return -1 if v1 < v2, 1 if v1 > v2 and 0 if v1 is equal to v2
 */
int compare_version(version_number v1, version_number v2) {
	if (v1.major_version < v2.major_version)
		return -1;
	else if (v1.major_version > v2.major_version)
		return 1;
	else {
		if (v1.minor_version < v2.minor_version)
			return -1;
		else if (v1.minor_version > v2.minor_version)
			return 1;
		else {
			if (v1.micro_version < v2.micro_version)
				return -1;
			else if (v1.micro_version > v2.micro_version)
				return 1;
			else {
				// Determine version type
				int v1_is_stable = !v1.rc_version && !v1.beta_version;
				int v2_is_stable = !v2.rc_version && !v2.beta_version;

				// Order: beta < rc < stable
				if (v1.beta_version && !v2.beta_version) return -1;  // beta < (rc or stable)
				if (v2.beta_version && !v1.beta_version) return 1;   // (rc or stable) > beta
				if (v1.rc_version && v2_is_stable) return -1;        // rc < stable
				if (v2.rc_version && v1_is_stable) return 1;         // stable > rc

				// Same type: compare patched_version
				if (v1.beta_version && v2.beta_version) {
					// Both beta versions
					if (v1.patched_version < v2.patched_version)
						return -1;
					else if (v1.patched_version > v2.patched_version)
						return 1;
				}
				else if (v1.rc_version && v2.rc_version) {
					// Both rc versions
					if (v1.patched_version < v2.patched_version)
						return -1;
					else if (v1.patched_version > v2.patched_version)
						return 1;
				}
				else if (v1_is_stable && v2_is_stable) {
					// Both stable versions: compare patches
					// 1.4.0 (patch=0) < 1.4.0-1 (patch=1)
					if (v1.patched_version < v2.patched_version)
						return -1;
					else if (v1.patched_version > v2.patched_version)
						return 1;
				}
			}
		}
	}
	return 0;
}

// ============================================================================
// NETWORK-DEPENDENT FUNCTIONS - Require libcurl
// ============================================================================

#if defined(HAVE_LIBCURL)

static version_number get_last_version_number(gchar *version_str) {
	gchar **v;
	version_number version = { 0 };

	v = g_strsplit_set(version_str, ".-", -1);

	if (v[0])
		version.major_version = g_ascii_strtoull(v[0], NULL, 10);
	if (v[0] && v[1])
		version.minor_version = g_ascii_strtoull(v[1], NULL, 10);
	if (v[0] && v[1] && v[2])
		version.micro_version = g_ascii_strtoull(v[2], NULL, 10);
	if (v[0] && v[1] && v[2] && v[3]) {
		remove_alpha(v[3], &version.rc_version, &version.beta_version);
		version.patched_version = g_ascii_strtoull(v[3], NULL, 10);
	}

	g_strfreev(v);
	return version;
}

static gboolean siril_update_get_highest(yyjson_doc *doc,
		gchar **highest_version, gint64 *release_timestamp,
		gint *build_revision, gchar **build_comment) {
	g_return_val_if_fail(highest_version != NULL, FALSE);
	g_return_val_if_fail(release_timestamp != NULL, FALSE);
	g_return_val_if_fail(build_revision != NULL, FALSE);
	g_return_val_if_fail(build_comment != NULL, FALSE);

	*highest_version = NULL;
	*release_timestamp = 0;
	*build_revision = 0;
	*build_comment = NULL;

	// Get root object and verify it's an object
	yyjson_val *root = yyjson_doc_get_root(doc);
	if (!root || !yyjson_is_obj(root)) {
		g_warning("Root is not an object");
		return FALSE;
	}

	// Get RELEASE array and verify it's an array
	yyjson_val *releases = yyjson_obj_get(root, "RELEASE");
	if (!releases || !yyjson_is_arr(releases)) {
		g_warning("RELEASE is not an array");
		return FALSE;
	}

	// Determine platform
	const char *platform;
	if (g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "windows") == 0 || g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "macos") == 0) {
		platform = SIRIL_BUILD_PLATFORM_FAMILY;
	} else {
		platform = "source";
	}

	// Iterate through releases array (versions are ordered newest first)
	yyjson_val *version_obj;
	size_t idx, max;
	yyjson_arr_foreach(releases, idx, max, version_obj) {
		if (!yyjson_is_obj(version_obj)) {
			g_warning("Version entry is not an object");
			continue;
		}

		// Get the platform-specific builds array
		yyjson_val *platform_builds = yyjson_obj_get(version_obj, platform);
		if (!platform_builds) {
			g_debug("No builds found for platform %s", platform);
			continue;
		}
		if (!yyjson_is_arr(platform_builds)) {
			g_warning("Platform builds is not an array");
			continue;
		}

		// Get version string early - we'll need it for any valid build
		yyjson_val *version = yyjson_obj_get(version_obj, "version");
		if (!version || !yyjson_is_str(version)) {
			g_warning("Version is missing or not a string");
			continue;
		}

		// Iterate through builds for this platform
		yyjson_val *build;
		size_t build_idx, build_max;
		yyjson_arr_foreach(platform_builds, build_idx, build_max, build) {
			if (!yyjson_is_obj(build)) {
				g_warning("Build entry is not an object");
				continue;
			}

			// Check build ID
			yyjson_val *build_id = yyjson_obj_get(build, "build-id");
			if (!build_id || !yyjson_is_str(build_id)) {
				if (strcmp(platform, "source") != 0) {
					g_debug("Build ID missing or not a string");
					continue;
				}
			}

			const char *build_id_str = build_id ? yyjson_get_str(build_id) : NULL;
			gboolean valid_build = FALSE;

			// For source platform, we don't need a build ID
			if (strcmp(platform, "source") == 0) {
				valid_build = TRUE;
			}
			// For other platforms, we need specific build IDs
			else if (build_id_str && (
				strcmp(build_id_str, "org.siril.Siril") == 0 ||
				strcmp(build_id_str, "org.free_astro.siril") == 0)) {
				valid_build = TRUE;
			}

			if (!valid_build) {
				continue;
			}

			// Get release date (from build or version)
			const char *release_date = NULL;
			yyjson_val *date = yyjson_obj_get(build, "date");
			if (date && yyjson_is_str(date)) {
				release_date = yyjson_get_str(date);
			} else {
				date = yyjson_obj_get(version_obj, "date");
				if (date && yyjson_is_str(date)) {
					release_date = yyjson_get_str(date);
				}
			}

			if (!release_date) {
				g_warning("No valid release date found");
				continue;
			}

			// We found a valid build - get the additional data
			*highest_version = g_strdup(yyjson_get_str(version));

			// Get optional build data
			yyjson_val *revision = yyjson_obj_get(build, "revision");
			if (revision && yyjson_is_int(revision)) {
				*build_revision = (gint) yyjson_get_int(revision);
			}

			yyjson_val *comment = yyjson_obj_get(build, "comment");
			if (comment && yyjson_is_str(comment)) {
				*build_comment = g_strdup(yyjson_get_str(comment));
			}

			// Parse release date
			gchar *str = g_strdup_printf("%s 00:00:00Z", release_date);
			GDateTime *datetime = g_date_time_new_from_iso8601(str, NULL);
			g_free(str);

			if (datetime) {
				*release_timestamp = g_date_time_to_unix(datetime);
				g_date_time_unref(datetime);
				return TRUE;
			} else {
				g_warning("Failed to parse release date: %s", release_date);
				g_clear_pointer(highest_version, g_free);
				g_clear_pointer(build_comment, g_free);
				*build_revision = 0;
				return FALSE;
			}
		}
	}

	return FALSE;
}

static gchar *parse_changelog(gchar *changelog) {
	gchar **token;
	GString *strResult;
	guint nargs, i;

	token = g_strsplit(changelog, "\n", -1);
	nargs = g_strv_length(token);

	strResult = g_string_new(token[0]);
	strResult = g_string_append(strResult, "\n\n");
	/* we start at line 3 */
	i = 3;
	while (i < nargs && token[i][0] != '\0') {
		strResult = g_string_append(strResult, token[i]);
		strResult = g_string_append(strResult, "\n");
		i++;
	}
	g_strfreev(token);
	return g_string_free(strResult, FALSE);
}

static gchar *check_version(gchar *version, gboolean *verbose, gchar **data) {
	gchar *changelog = NULL;
	gchar *msg = NULL;

	version_number last_version_available = get_last_version_number(version);
	version_number current_version = get_current_version_number();
	guint x = last_version_available.major_version;
	guint y = last_version_available.minor_version;
	guint z = last_version_available.micro_version;
	if (x == 0 && y == 0 && z == 0) {
		if (*verbose)
			msg = siril_log_message(_("No update check: cannot fetch version file\n"));
	} else {
		if (compare_version(current_version, last_version_available) < 0) {
			gchar *url = NULL;

			if (last_version_available.patched_version != 0) {
				const gchar *str;
				if (last_version_available.beta_version || last_version_available.rc_version) {
					str = "%s%d.%d.%d-%d";
				} else {
					str = "%s%d.%d.%d.%d";
				}
				url = g_strdup_printf(str, SIRIL_DOWNLOAD, last_version_available.major_version, last_version_available.minor_version, last_version_available.micro_version, last_version_available.patched_version);
			} else {
				url = g_strdup_printf("%s%d.%d.%d",
						SIRIL_DOWNLOAD, last_version_available.major_version, last_version_available.minor_version, last_version_available.micro_version);
			}

			msg = siril_log_message(_("New version is available. You can download it at <a href=\"%s\">%s</a>\n"), url, url);
			g_free(url);
			GString *urlstring = g_string_new(GITLAB_URL);
			g_string_append_printf(urlstring, "/%s/ChangeLog", version);
			gchar *changelog_url = g_string_free(urlstring, FALSE);
			gsize length;
			int error;
			changelog = fetch_url(changelog_url, &length, &error, FALSE);
			g_free(changelog_url);
			if (error)
				return NULL;
			if (changelog) {
				*data = parse_changelog(changelog);
				/* force the verbose variable */
				*verbose = TRUE;
			}
		} else if (compare_version(current_version, last_version_available) > 0) {
			if (*verbose)
				msg = siril_log_message(_("No update check: this is a development version\n"));
		} else {
			if (*verbose)
				msg = siril_log_message(_("Siril is up to date\n"));
		}
		g_free(changelog);
	}
	return msg;
}

static gchar *check_update_version(fetch_url_async_data *args) {
	gchar *last_version = NULL;
	gchar *build_comment = NULL;
	gint64 release_timestamp = 0;
	gint build_revision = 0;
	gchar *msg = NULL;
	gchar *data = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;

	// Parse JSON
	yyjson_read_err err = { 0 };
	yyjson_doc *doc = yyjson_read(args->content, strlen(args->content), 0);
	if (!doc) {
		g_printerr("%s: parsing of %s failed: %s\n", G_STRFUNC,
				   args->url, err.msg);
		return NULL;
	}

	if (siril_update_get_highest(doc, &last_version, &release_timestamp, &build_revision, &build_comment)) {
		g_fprintf(stdout, "Last available version: %s\n", last_version);

		msg = check_version(last_version, &(args->verbose), &data);
		message_type = GTK_MESSAGE_INFO;
	} else {
		msg = siril_log_message(_("Cannot fetch version file\n"));
	}

	if (args->verbose) {
		set_cursor_waiting(FALSE);
		if (msg) {
			siril_data_dialog(message_type, _("Software Update"), msg, data);
		}
	}

	g_clear_pointer(&last_version, g_free);
	g_clear_pointer(&build_comment, g_free);
	yyjson_doc_free(doc);

	return msg;
}

static gboolean end_update_idle(gpointer p) {
	fetch_url_async_data *args = (fetch_url_async_data *) p;
	if (args->content)
		check_update_version(args);

	/* free data */
	set_cursor_waiting(FALSE);
	free(args->content);
	free(args);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

// Define the notification struct
typedef struct _notification {
	GString *messageString;
	int status;
} notification;

static int parseJsonNotificationsString(const char *jsonString, GSList **validNotifications) {
	// Parse JSON from string using yyjson
	yyjson_doc *doc = yyjson_read(jsonString, strlen(jsonString), YYJSON_READ_NOFLAG);
	if (!doc) {
		siril_log_color_message(_("Error parsing JSON from URL: Failed to parse JSON\n"), "red");
		return 1;
	}

	yyjson_val *root = yyjson_doc_get_root(doc);
	if (!yyjson_is_obj(root) && !yyjson_is_arr(root)) {
		siril_log_color_message(_("Error parsing JSON from URL: JSON root is not an object or array\n"), "red");
		yyjson_doc_free(doc);
		return 1;
	}

	// Get the current time
	GDateTime *currentTime = g_date_time_new_now_local();

	size_t length = 0;
	yyjson_val *messages = NULL;

	// Handle messages as an array or single object
	if (yyjson_is_arr(root)) {
		messages = root;
		length = yyjson_arr_size(messages);
	} else if (yyjson_is_obj(root)) {
		length = 1;  // Single object as message
	} else {
		siril_log_color_message(_("Error parsing JSON from URL: Invalid root JSON structure\n"), "red");
		yyjson_doc_free(doc);
		g_date_time_unref(currentTime);
		return 1;
	}

	// Iterate over messages
	for (size_t i = 0; i < length; i++) {
		yyjson_val *message = (!yyjson_is_arr(root)) ? root : yyjson_arr_get(messages, i);
		if (!yyjson_is_obj(message)) {
			siril_log_color_message(_("Error parsing JSON from URL: Message is not a valid object\n"), "red");
			continue;
		}

		// Check for required fields
		const char *validFromStr = yyjson_get_str(yyjson_obj_get(message, "valid-from"));
		const char *validToStr = yyjson_get_str(yyjson_obj_get(message, "valid-to"));
		const char *messageStr = yyjson_get_str(yyjson_obj_get(message, "message"));
		yyjson_val *priority_val = yyjson_obj_get(message, "priority");

		if (!validFromStr || !validToStr || !messageStr || !priority_val || !yyjson_is_int(priority_val)) {
			siril_log_color_message(_("Error parsing JSON from URL: Required fields missing or invalid\n"), "red");
			continue;
		}

		// Parse valid-from and valid-to fields
		GDateTime *validFrom = g_date_time_new_from_iso8601(validFromStr, NULL);
		GDateTime *validTo = g_date_time_new_from_iso8601(validToStr, NULL);
		if (!validFrom || !validTo) {
			siril_log_color_message(_("Error parsing JSON from URL: Invalid ISO8601 date format\n"), "red");
			if (validFrom) g_date_time_unref(validFrom);
			if (validTo) g_date_time_unref(validTo);
			continue;
		}

		// Check for optional version-from and version-to fields
		version_number empty_version = {0};
		version_number current_version = get_current_version_number();
		version_number valid_from_version = {0};
		version_number valid_to_version = {0};

		const char *versionFromStr = yyjson_get_str(yyjson_obj_get(message, "version-from"));
		if (versionFromStr) {
			valid_from_version = get_version_number_from_string(versionFromStr);
		}

		const char *versionToStr = yyjson_get_str(yyjson_obj_get(message, "version-to"));
		if (versionToStr) {
			valid_to_version = get_version_number_from_string(versionToStr);
		}

		// Parse priority
		int status = (int)yyjson_get_int(priority_val);
		gboolean valid = TRUE;

		// Check if current time is within the validity period
		if (!(g_date_time_compare(currentTime, validFrom) >= 0 && g_date_time_compare(currentTime, validTo) <= 0)) {
			valid = FALSE;
		}

		// Check if current version matches version constraints
		if (memcmp(&valid_from_version, &empty_version, sizeof(version_number)) &&
			compare_version(current_version, valid_from_version) < 0) {
			valid = FALSE;
			}
			if (memcmp(&valid_to_version, &empty_version, sizeof(version_number)) &&
				compare_version(current_version, valid_to_version) > 0) {
				valid = FALSE;
				}

				if (valid) {
					// Create and populate notification
					notification *notif = g_new(notification, 1);
					notif->messageString = g_string_new(messageStr);
					notif->status = status;

					// Append notification to the list
					*validNotifications = g_slist_append(*validNotifications, notif);
				}

				// Free allocated memory
				g_date_time_unref(validFrom);
				g_date_time_unref(validTo);
	}

	// Clean up
	g_date_time_unref(currentTime);
	yyjson_doc_free(doc);

	return 0;
}

static gboolean end_notifier_idle(gpointer p) {
	fetch_url_async_data *args = (fetch_url_async_data *) p;
	if (!args->content)
		goto end_notifier_idle_error;
	GSList *validNotifications = NULL;

	control_window_switch_to_tab(OUTPUT_LOGS);

	// Fetch and parse JSON file from URL and populate validNotifications list
	if (parseJsonNotificationsString(args->content, &validNotifications) != 0) {
		siril_log_message(_("Error fetching or parsing Siril notifications JSON file from URL\n"));
		goto end_notifier_idle_error;
	}

	// Print and then free valid notifications
	for (GSList *iter = validNotifications; iter; iter = iter->next) {
		notification *notif = (notification *) iter->data;
		char *color = notif->status == 1 ? "green" : notif->status == 2 ? "salmon" : "red";
		siril_log_color_message(_("*** SIRIL NOTIFICATION ***\n%s\n"), color, notif->messageString->str);

		// Free allocated memory for notification
		g_string_free(notif->messageString, TRUE);
		g_free(notif);
	}
	g_slist_free(validNotifications);

end_notifier_idle_error:

	set_cursor_waiting(FALSE);
	/* free data */
	free(args->content);
	free(args);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

void siril_check_updates(gboolean verbose) {
	if (!is_online()) {
		siril_log_color_message(_("Error: Siril is in offline mode, cannot check updates.\n"), "red");
		return;
	}
	fetch_url_async_data *args = calloc(1, sizeof(fetch_url_async_data));
	args->url = g_strdup(SIRIL_VERSIONS);
	args->content = NULL;
	args->verbose = verbose;
	args->idle_function = end_update_idle;

	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);

	// this is a graphical operation, we don't use the main processing thread for it, it could block file opening
	g_thread_new("siril-update", fetch_url_async, args);
}

void siril_check_notifications(gboolean verbose) {
	fetch_url_async_data *args;

	args = calloc(1, sizeof(fetch_url_async_data));
	GString *url = g_string_new(GITLAB_URL);
	g_string_append_printf(url, "/%s/%s", BRANCH, SIRIL_NOTIFICATIONS);
	args->url = g_string_free(url, FALSE);
	siril_debug_print("Notification URL: %s\n", args->url);
	args->content = NULL;
	args->verbose = verbose;
	args->idle_function = end_notifier_idle;
	siril_debug_print("Checking notifications...\n");
	set_progress_bar_data(_("Looking for notifications..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);

	// this is a graphical operation, we don't use the main processing thread for it, it could block file opening
	g_thread_new("siril-notifications", fetch_url_async, args);
}

#endif // HAVE_LIBCURL