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

#include <string.h>

#include "core/siril.h"
#include "core/siril_networking.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"

static gboolean online_status = TRUE;

static const int DEFAULT_FETCH_RETRIES = 5;

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;

	mem->data = realloc(mem->data, mem->len + realsize + 1);

	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;

	return realsize;
}

/* Fetches a URL asynchronously. The caller calls this function and sets a callback
 * (args->idle_function) to handle the result.
 * */

gpointer fetch_url_async(gpointer p) {
	CURL *curl = NULL;
	struct ucontent *content;
	gchar *result = NULL;
	long code = -1L;
	int retries;
	unsigned int s;
	fetch_url_async_data *args = (fetch_url_async_data *) p;
	g_fprintf(stdout, "fetch_url(): %s\n", args->url);

	curl = curl_easy_init();
	if (g_getenv("CURL_CA_BUNDLE"))
		if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE")))
			siril_log_color_message(_("Error configuring CURL with CA bundle. https functionality unavailable.\n"), "red");

	if (!curl) {
		siril_log_color_message(_("Error initialising CURL handle, URL functionality unavailable.\n"), "red");
		g_free(args->url);
		free(args);
		return NULL;
	}
	args->curl = curl;

	content = calloc(1, sizeof(struct ucontent));
	if (content == NULL) {
		PRINT_ALLOC_ERR;
		g_free(args->url);
		free(args);
		return NULL;
	}

	set_progress_bar_data(NULL, 0.1);

	retries = DEFAULT_FETCH_RETRIES;

	retrieve:

	content->data = calloc(1, 1);
	if (content->data == NULL) {
		PRINT_ALLOC_ERR;
		g_free(args->url);
		free(args);
		free(content);
		return NULL;
	}

	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, args->url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, "siril/0.0");
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	if (retval) {
		siril_debug_print("Error in curl_easy_setopt()\n");
		goto failed_curl_easy_setopt;
	}

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

				free(content->data);
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
	} else {
		siril_log_color_message(_("Cannot retrieve information from the update URL. Error: [%ld]\n"), "red", retval);
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);
	args->content = result;

failed_curl_easy_setopt:

	if (result) {
		gdk_threads_add_idle(args->idle_function, args);
	} else {
		g_free(args->url);
		free(args);
		free(content->data);
	}
	free(content);

	curl_easy_cleanup(curl);
	return NULL;
}

char *fetch_url(const gchar *url, gsize *length) {
	CURL *curl = NULL;
	struct ucontent *content = malloc(sizeof(struct ucontent));
	char *result = NULL;
	long code;
	int retries;
	unsigned int s;

	curl = curl_easy_init();
	if (g_getenv("CURL_CA_BUNDLE"))
		if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE")))
			siril_debug_print("Error in curl_easy_setopt()\n");

	if (!curl) {
		siril_log_color_message(_("Error initialising CURL handle, URL functionality unavailable.\n"), "red");
		free(content);
		return NULL;
	}
	retries = DEFAULT_FETCH_RETRIES;

retrieve:
	content->data = malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	CURLcode ret = curl_easy_setopt(curl, CURLOPT_URL, url);
	ret |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	ret |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	ret |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	ret |= curl_easy_setopt(curl, CURLOPT_USERAGENT, PACKAGE_STRING);
	if (ret)
		siril_debug_print("Error in curl_easy_setopt()\n");

	siril_debug_print("fetch_url(): %s\n", url);
	if (curl_easy_perform(curl) == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

		switch (code) {
		case 200:
			result = content->data;
			siril_debug_print("downloaded %zd bytes\n", content->len);
			break;
		case 500:
		case 502:
		case 503:
		case 504:
			siril_debug_print("Failed to download page %s (error %ld)\n", url, code);
			if (retries) {
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				siril_debug_print("Wait %us before retry\n", s);
				g_usleep(s * 1E6);
				free(content->data);
				retries--;
				goto retrieve;
			} else {
				siril_log_color_message(_("After %ld tries, Server unreachable or unresponsive. (%s)\n"), "salmon", DEFAULT_FETCH_RETRIES, content->data);
			}
			break;
		default:
			siril_log_color_message(_("Server unreachable or unresponsive (code %ld).\nServer returned: %s\n"), "red", code, content->data);
			break;
		}
	}
	else siril_log_color_message(_("Internet Connection failure.\n"), "red");

	curl_easy_cleanup(curl);
	curl = NULL;
	if (!result || content->len == 0) {
		free(content->data);
		result = NULL;
	}
	*length = content->len;
	free(content);
	return result;
}

int submit_post_request(const char *url, const char *post_data, char **post_response) {
	CURL *curl;
	CURLcode res = CURLE_OK;
	struct ucontent chunk;

	chunk.data = malloc(1);  /* will be grown as needed by realloc above */
	chunk.len = 0;    /* no data at this point */

	curl_global_init(CURL_GLOBAL_ALL);
	curl = curl_easy_init();
	if (curl) {
		res = curl_easy_setopt(curl, CURLOPT_URL, url);
		// send all data to this function
		if (!res) res = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
		// we pass our 'chunk' struct to the callback function
		if (!res) res = curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&chunk);
		// some servers do not like requests that are made without a user-agent
		// field, so we provide one
		if (!res) res = curl_easy_setopt(curl, CURLOPT_USERAGENT, "libcurl-agent/1.0");
		if (!res) res = curl_easy_setopt(curl, CURLOPT_POSTFIELDS, post_data);
		// if we do not provide POSTFIELDSIZE, libcurl will strlen() by itself
		if (!res) res = curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)strlen(post_data));
		// Perform the request, res will get the return code
		if (!res) res = curl_easy_perform(curl);
		// Check for errors
		if(res != CURLE_OK) {
			siril_log_color_message(_("Error fetching URL: %s\n"), "red", curl_easy_strerror(res));
		} else {
			// Set the response data
			*post_response = g_strdup(chunk.data);
		}
		// always cleanup
		curl_easy_cleanup(curl);
	} else {
		res = CURLE_FAILED_INIT;
	}

	free(chunk.data);
	curl_global_cleanup();
    return (res != CURLE_OK ? 1 : 0);
}

gboolean siril_compiled_with_networking() {
	return TRUE;
}

#else

gboolean siril_compiled_with_networking() {
	return FALSE;
}

#endif

gboolean is_online() {
	return (siril_compiled_with_networking() && online_status);
}

gboolean set_online_status(gboolean status) {
	online_status = status;
	if (!siril_compiled_with_networking())
		return FALSE;
	return online_status;
}
