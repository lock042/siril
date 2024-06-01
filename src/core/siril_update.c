/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#if defined(HAVE_LIBCURL)
#include <json-glib/json-glib.h>
#include <curl/curl.h>

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "core/siril_update.h"


#define DOMAIN "https://siril.org/"
#define SIRIL_VERSIONS DOMAIN"siril_versions.json"
#define SIRIL_DOWNLOAD DOMAIN"download/"
#define GITLAB_URL "https://gitlab.com/free-astro/siril/raw"
#define BRANCH "notifier"
//#define BRANCH "master"
#define SIRIL_NOTIFICATIONS "notifications/siril_notifications.json"

static gchar *get_changelog(gchar *str);

// taken from gimp
static gboolean siril_update_get_highest(JsonParser *parser,
		gchar **highest_version, gint64 *release_timestamp,
		gint *build_revision, gchar **build_comment) {
	JsonPath *path;
	JsonNode *result;
	JsonArray *versions;
	const gchar *platform;
	const gchar *path_str;
	const gchar *release_date = NULL;
	GError *error = NULL;
	gint i;

	g_return_val_if_fail(highest_version != NULL, FALSE);
	g_return_val_if_fail(release_timestamp != NULL, FALSE);
	g_return_val_if_fail(build_revision != NULL, FALSE);
	g_return_val_if_fail(build_comment != NULL, FALSE);

	*highest_version = NULL;
	*release_timestamp = 0;
	*build_revision = 0;
	*build_comment = NULL;

	path_str = "$['RELEASE'][*]";

	/* For Windows and macOS, let's look if installers are available.
	 * For other platforms, let's just look for source release.
	 */
	if (g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "windows") == 0
			|| g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "macos") == 0)
		platform = SIRIL_BUILD_PLATFORM_FAMILY;
	else
		platform = "source";

	path = json_path_new();
	/* Ideally we could just use Json path filters like this to
	 * retrieve only released binaries for a given platform:
	 * g_strdup_printf ("$['STABLE'][?(@.%s)]['version']", platform);
	 * json_array_get_string_element (result, 0);
	 * And that would be it! We'd have our last release for given
	 * platform.
	 * Unfortunately json-glib does not support filter syntax, so we
	 * end up looping through releases.
	 */
	if (!json_path_compile(path, path_str, &error)) {
		g_warning("%s: path compilation failed: %s\n", G_STRFUNC,
				error->message);
		g_clear_error(&error);
		g_object_unref(path);

		return FALSE;
	}
	result = json_path_match(path, json_parser_get_root(parser));
	if (!JSON_NODE_HOLDS_ARRAY(result)) {
		g_printerr("%s: match for \"%s\" is not a JSON array.\n",
		G_STRFUNC, path_str);
		g_object_unref(path);

		return FALSE;
	}

	versions = json_node_get_array(result);
	for (i = 0; i < (gint) json_array_get_length(versions); i++) {
		JsonObject *version;

		/* Note that we don't actually look for the highest version,
		 * but for the highest version for which a build for your
		 * platform (and optional build-id) is available.
		 *
		 * So we loop through the version list then the build array
		 * and break at first compatible release, since JSON arrays
		 * are ordered.
		 */
		version = json_array_get_object_element(versions, i);
		if (json_object_has_member(version, platform)) {
			JsonArray *builds;
			gint j;

			builds = json_object_get_array_member(version, platform);

			for (j = 0; j < (gint) json_array_get_length(builds); j++) {
				const gchar *build_id = NULL;
				JsonObject *build;

				build = json_array_get_object_element(builds, j);
				if (json_object_has_member(build, "build-id"))
					build_id = json_object_get_string_member (build, "build-id");
				if (g_strcmp0(build_id, "org.free_astro.siril") == 0
						|| g_strcmp0(platform, "source") == 0) {
					/* Release date is the build date if any set,
					 * otherwise the main version release date.
					 */
					if (json_object_has_member(build, "date"))
						release_date = json_object_get_string_member(build, "date");
					else
						release_date = json_object_get_string_member(version, "date");

					/* These are optional data. */
					if (json_object_has_member(build, "revision"))
						*build_revision = json_object_get_int_member(build, "revision");
					if (json_object_has_member(build, "comment"))
						*build_comment = g_strdup(json_object_get_string_member(build, "comment"));
					break;
				}
			}

			if (release_date) {
				*highest_version = g_strdup(json_object_get_string_member(version, "version"));
				break;
			}
		}
	}

	if (*highest_version && *release_date) {
		GDateTime *datetime;
		gchar *str;

		str = g_strdup_printf("%s 00:00:00Z", release_date);
		datetime = g_date_time_new_from_iso8601(str, NULL);
		g_free(str);

		if (datetime) {
			*release_timestamp = g_date_time_to_unix(datetime);
			g_date_time_unref(datetime);
		} else {
			/* JSON file data bug. */
			g_printerr("%s: release date for version %s not properly formatted: %s\n",
					G_STRFUNC, *highest_version, release_date);

			g_clear_pointer(highest_version, g_free);
			g_clear_pointer(build_comment, g_free);
			*build_revision = 0;
		}
	}

	json_node_unref(result);
	g_object_unref(path);

	return (*highest_version != NULL);
}

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

static version_number get_current_version_number() {
	gchar **fullVersionNumber;
	version_number version;

	fullVersionNumber = g_strsplit_set(PACKAGE_VERSION, ".-", -1);
	version.major_version = g_ascii_strtoull(fullVersionNumber[0], NULL, 10);
	version.minor_version = g_ascii_strtoull(fullVersionNumber[1], NULL, 10);
	version.micro_version = g_ascii_strtoull(fullVersionNumber[2], NULL, 10);
	if (fullVersionNumber[3] == NULL) {
		version.patched_version = 0;
		version.rc_version = FALSE;
		version.beta_version = FALSE;
	} else {
		version.patched_version = check_for_patch(fullVersionNumber[3], &version.rc_version, &version.beta_version);
	}

	g_strfreev(fullVersionNumber);

	return version;
}

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

/**
 * This function compare x1.y1.z1.patch1 vs x2.y2.z2.patch2
 * @param v1 First version number to be tested
 * @param v2 Second version number to be tested
 * @return -1 if v1 < v2, 1 if v1 > v2 and 0 if v1 is equal to v2
 */
static int compare_version(version_number v1, version_number v2) {
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
				if (v1.beta_version && v2.rc_version) return -1;
				if (v2.beta_version && v1.rc_version) return 1;
				if (v1.beta_version && !v2.rc_version && !v2.beta_version) return -1;
				if (v1.rc_version && !v2.rc_version && !v2.beta_version) return -1;
				if (v2.rc_version && !v1.rc_version && !v1.beta_version) return 1;

				/* check for patched version */
				if ((!v1.rc_version && !v2.rc_version) || (!v1.beta_version && !v2.beta_version) ||
						(v1.rc_version && v2.rc_version) || (v1.beta_version && v2.beta_version)) {
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
			changelog = get_changelog(version);
			if (changelog) {
				*data = parse_changelog(changelog);
				/* force the verbose variable */
				*verbose = TRUE;
			}
			g_free(url);
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

struct _update_data {
	gchar *url;
	long code;
	gchar *content;
	gboolean verbose;
	gboolean (*idle_function)(gpointer args);
};


static const int DEFAULT_FETCH_RETRIES = 5;
static CURL *curl;

struct ucontent {
	gchar *data;
	size_t len;
};

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;

	mem->data = g_try_realloc(mem->data, mem->len + realsize + 1);

	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;

	return realsize;
}

static void init() {
	if (!curl) {
		g_fprintf(stdout, "initializing CURL\n");
		curl_global_init(CURL_GLOBAL_ALL);
		curl = curl_easy_init();
		if (g_getenv("CURL_CA_BUNDLE"))
			if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE")))
				siril_debug_print("Error in curl_easy_setopt()\n");
	}

	if (!curl)
		exit(EXIT_FAILURE);
}

static void http_cleanup() {
	curl_easy_cleanup(curl);
	curl_global_cleanup();
	curl = NULL;
}

static gchar *get_changelog(gchar *str) {
	struct ucontent *changelog;
	gchar *result = NULL;
	gchar *changelog_url;
	long code;

	changelog = g_try_malloc(sizeof(struct ucontent));
	if (changelog == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	changelog->data = g_malloc(1);
	changelog->data[0] = '\0';
	changelog->len = 0;

	GString *url = g_string_new(GITLAB_URL);
	url = g_string_append(url, "/");
	url = g_string_append(url, str);
	url = g_string_append(url, "/ChangeLog");

	changelog_url = g_string_free(url, FALSE);
	CURLcode ret;
	ret = curl_easy_setopt(curl, CURLOPT_URL, changelog_url);
	ret |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	ret |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, changelog);
	if (ret)
		siril_debug_print("Error in curl_easy_setopt()\n");
	if (curl_easy_perform(curl) == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
		if (code == 200) {
			result = g_strdup(changelog->data);
		}
	}
	if (!result)
		g_free(changelog->data);
	g_free(changelog);
	g_free(changelog_url);

	return result;
}

static gchar *check_update_version(struct _update_data *args) {
	JsonParser *parser;
	gchar *last_version = NULL;
	gchar *build_comment = NULL;
	gint64 release_timestamp = 0;
	gint build_revision = 0;
	GError *error = NULL;
	gchar *msg = NULL;
	gchar *data = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;

	parser = json_parser_new();
	if (!json_parser_load_from_data(parser, args->content, -1, &error)) {
		g_printerr("%s: parsing of %s failed: %s\n", G_STRFUNC,
				args->url, error->message);
		g_clear_object(&parser);
		g_clear_error(&error);

		return NULL;
	}

	siril_update_get_highest(parser, &last_version, &release_timestamp,	&build_revision, &build_comment);

	if (last_version) {
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
	g_object_unref(parser);

	return msg;
}

static gboolean end_update_idle(gpointer p) {
	struct _update_data *args = (struct _update_data *) p;

	if (args->content == NULL) {
		switch(args->code) {
		case 0:
			siril_log_message(_("Unable to check updates! "
					"Please Check your network connection\n"));
			break;
		default:
			siril_log_message(_("Unable to check updates! Error: %ld\n"),
					args->code);
		}
	} else {
		check_update_version(args);
	}

	/* free data */
	set_cursor_waiting(FALSE);
	g_free(args->content);
	free(args);
	http_cleanup();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

static void array_string_element_clear(GString *element) {
	g_free(element->str); // free the gchar* inside the GString but not the string itself
}

static int parseJsonNotificationsString(const gchar *jsonString, GArray *validMessages, GArray *validStatus) {
	// Load JSON from string and check for errors
	GError *error = NULL;
	JsonParser *parser = json_parser_new();
	if (!(json_parser_load_from_data(parser, jsonString, -1, &error))) {
		siril_log_color_message(_("Error parsing JSON from URL: %s\n"), "red", error->message);
		g_object_unref(parser);
		return 1;
	}

	// Get root node
	JsonNode *node = json_parser_get_root(parser);
	if (!node) {
		siril_log_color_message(_("Error parsing JSON from URL: unable to get root node\n"), "red");
		g_object_unref(parser);
		return 1;
	}

	// Get current time
	GDateTime *currentTime = g_date_time_new_now_local();

	// The message array / object and number of messages (length is overwritten if an array is found)
	JsonArray *messages = NULL;
	JsonObject *single_message = NULL;
	guint length = 1;

	// Get array of messages
	if (JSON_NODE_HOLDS_ARRAY(node)) {
		messages = json_node_get_array(node);
		length = json_array_get_length(messages);
	} else {
		if (JSON_NODE_HOLDS_OBJECT(node)) {
			single_message = json_node_get_object(node);
		} else {
			siril_log_color_message(_("Error parsing JSON from URL: unable to find a valid JSON object\n"), "red");
			g_object_unref(parser);
			return 1;
		}
	}
	if (!messages && !single_message) {
		siril_log_color_message(_("Error parsing JSON from URL: unable to find any valid JSON objects\n"), "red");
		g_object_unref(parser);
		return 1;
	}

	// Iterate over messages
	for (guint i = 0; i < length; i++) {
		if (!single_message)
			single_message = json_array_get_object_element(messages, i);
		// Check the necessary json objects are there
		if (!(json_object_has_member(single_message, "valid-from")) &&
				(json_object_has_member(single_message, "valid-to")) &&
				(json_object_has_member(single_message, "message")) &&
				(json_object_has_member(single_message, "priority"))) {
			siril_log_color_message(_("Error parsing JSON from URL: required JSON members not found\n"), "red");
			// No need to clean up any partially-filled GArrays as they will be freed by the caller
		g_object_unref(parser);
			return 1;
		}

		// Parse valid-from and valid-to fields
		GDateTime *validFrom = g_date_time_new_from_iso8601(json_object_get_string_member(single_message, "valid-from"), NULL);
		GDateTime *validTo = g_date_time_new_from_iso8601(json_object_get_string_member(single_message, "valid-to"), NULL);

		// Parse status and append to GArray
		int *status = g_malloc(sizeof(int));
		*status = json_object_get_int_member(single_message, "priority");
		g_array_append_val(validStatus, status);

		// Check if current time is within validity period
		if (g_date_time_compare(currentTime, validFrom) >= 0 && g_date_time_compare(currentTime, validTo) <= 0) {
			// Store the message in a variable first
			GString *messageString = g_string_new(json_object_get_string_member(single_message, "message"));
			// Append the message to the array
			g_array_append_val(validMessages, messageString);
		}

		// Free allocated memory
		g_date_time_unref(validFrom);
		g_date_time_unref(validTo);

		// Reset the JsonObject
		single_message = NULL;
	}

	// Free allocated memory
	g_date_time_unref(currentTime);
	g_object_unref(parser);

	return 0;
}

static gboolean end_notifier_idle(gpointer p) {
	struct _update_data *args = (struct _update_data *) p;
	GArray *validMessages = NULL, *validStatus = NULL;
	if (args->content == NULL) {
		switch(args->code) {
		case 0:
			siril_log_message(_("Unable to check notifications! "
					"Please Check your network connection\n"));
			break;
		default:
			siril_log_message(_("Unable to check notifications! Error: %ld\n"), args->code);
		}
	} else {
		control_window_switch_to_tab(OUTPUT_LOGS);

		// Create an array to store valid messages
		validMessages = g_array_new(FALSE, FALSE, sizeof(GString*));
		g_array_set_clear_func(validMessages, (GDestroyNotify) array_string_element_clear);
		validStatus = g_array_new(FALSE, FALSE, sizeof(int*));

		// Fetch and parse JSON file from URL and populate validMessages array
		if (parseJsonNotificationsString(args->content, validMessages, validStatus) != 0) {
			siril_log_message(_("Error fetching or parsing Siril notifications JSON file from URL\n"));
			goto end_notifier_idle_error;
		}

		// Print and then free valid messages
		for (guint i = 0; i < validMessages->len; i++) {
			int *status = g_array_index(validStatus, int*, i);
			GString *msg = g_array_index(validMessages, GString*, i);
			char *color = *status == 1 ? "green" : *status == 2 ? "salmon" : "red";
			siril_log_color_message(_("*** SIRIL NOTIFICATION ***\n%s\n"), color, msg->str);
		}
	}

end_notifier_idle_error:

	set_cursor_waiting(FALSE);

	/* free data */
	if (validMessages)
		g_array_free(validMessages, TRUE);
	if (validStatus)
		g_array_free(validStatus, TRUE);
	g_free(args->content);
	g_free(args->url);
	free(args);
	http_cleanup();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

static gpointer fetch_url(gpointer p) {
	struct ucontent *content;
	gchar *result;
	long code = -1L;
	int retries;
	unsigned int s;
	struct _update_data *args = (struct _update_data *) p;

	content = g_try_malloc(sizeof(struct ucontent));
	if (content == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	g_fprintf(stdout, "fetch_url(): %s\n", args->url);

	init();
	set_progress_bar_data(NULL, 0.1);

	result = NULL;

	retries = DEFAULT_FETCH_RETRIES;

	retrieve: content->data = g_malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, args->url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, "siril/0.0");
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	if (retval)
		siril_debug_print("Error in curl_easy_setopt()\n");

	retval = curl_easy_perform(curl);
	if (retval == CURLE_OK) {
		if (retries == DEFAULT_FETCH_RETRIES) set_progress_bar_data(NULL, 0.4);
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
		if (retries == DEFAULT_FETCH_RETRIES) set_progress_bar_data(NULL, 0.6);

		switch (code) {
		case 200:
			result = content->data;
			break;
		case 500:
		case 502:
		case 503:
		case 504:
			g_fprintf(stderr, "Fetch failed with code %ld for URL %s\n", code,
					args->url);

			if (retries && get_thread_run()) {
				double progress = (DEFAULT_FETCH_RETRIES - retries) / (double) DEFAULT_FETCH_RETRIES;
				progress *= 0.4;
				progress += 0.6;
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				char *msg = siril_log_message(_("Error: %ld. Wait %us before retry\n"), code, s);
				msg[strlen(msg) - 1] = 0; /* we remove '\n' at the end */
				set_progress_bar_data(msg, progress);
				g_usleep(s * 1E6);

				g_free(content->data);
				retries--;
				goto retrieve;
			}

			break;
		default:
			g_fprintf(stderr, "Fetch failed with code %ld for URL %s\n", code,
					args->url);
		}
		g_fprintf(stderr, "Fetch succeeded with code %ld for URL %s\n", code,
				args->url);
		args->code = code;
	} else {
		siril_log_color_message(_("Cannot retrieve information from the update URL. Error: [%ld]\n"), "red", retval);
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);

	if (!result)
		g_free(content->data);
	g_free(content);

	args->content = result;

	gdk_threads_add_idle(args->idle_function, args);
	return NULL;
}

void siril_check_updates(gboolean verbose) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	args->url = g_strdup(SIRIL_VERSIONS);
	args->code = 0L;
	args->content = NULL;
	args->verbose = verbose;
	args->idle_function = end_update_idle;

	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);

	// this is a graphical operation, we don't use the main processing thread for it, it could block file opening
	g_thread_new("siril-update", fetch_url, args);
}

void siril_check_notifications(gboolean verbose) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	GString *url = g_string_new(GITLAB_URL);
	g_string_append_printf(url, "/%s/%s", BRANCH, SIRIL_NOTIFICATIONS);
	args->url = g_string_free(url, FALSE);
	siril_debug_print("Notification URL: %s\n", args->url);
	args->code = 0L;
	args->content = NULL;
	args->verbose = verbose;
	args->idle_function = end_notifier_idle;
	siril_debug_print("Checking notifications...\n");
	set_progress_bar_data(_("Looking for notifications..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);

	// this is a graphical operation, we don't use the main processing thread for it, it could block file opening
	g_thread_new("siril-notifications", fetch_url, args);
}

#endif
