/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#endif

#include "core/siril.h"
#include "core/siril_networking.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"

#define STR_INDIR(x) #x 
#define STR(x) STR_INDIR(x)
#define SIRIL_USER_AGENT "siril/" STR(SIRIL_MAJOR_VERSION) "." STR(SIRIL_MINOR_VERSION) " (https://gitlab.com/free-astro/siril/)"

// Uncomment the next line for some additional debug printing
// #define NETWORKING_DEBUG

static gboolean online_status = TRUE;

#if defined(HAVE_LIBCURL)

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;
	mem->data = realloc(mem->data, mem->len + realsize + 1);
	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;
	return realsize;
}

typedef enum {
	HTTP_GET,
	HTTP_POST
} HttpRequestType;

static CURL* initialize_curl(const gchar *url, struct ucontent *content, HttpRequestType request_type, const gchar *post_data) {
	CURL *curl = curl_easy_init();
	if (!curl) {
		siril_log_color_message(_("Error initialising CURL handle, URL functionality unavailable.\n"), "red");
		return NULL;
	}
	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, SIRIL_USER_AGENT);
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	if (request_type == HTTP_POST) {
		retval |= curl_easy_setopt(curl, CURLOPT_POSTFIELDS, post_data);
		retval |= curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)strlen(post_data));
	}
	if (retval) {
		siril_debug_print("Error in curl_easy_setopt()\n");
		curl_easy_cleanup(curl);
		return NULL;
	}
	if (g_getenv("CURL_CA_BUNDLE")) {
		if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE"))) {
			siril_log_color_message(_("Error configuring CURL with CA bundle. https functionality unavailable.\n"), "red");
		}
	}
	return curl;
}

static char* handle_curl_response(CURL *curl, struct ucontent *content, const gchar *url, long *code, gboolean verbose) {
	char *result = NULL;
	content->data = calloc(1, 1);
	if (content->data == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	CURLcode retval = curl_easy_perform(curl);
	if (retval == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, code);
		switch (*code) {
			case 200: // 200: OK.
				result = content->data;
				break;
			// 201: only relevant to POST requests which don't use this function.
			// 202: Accepted. Typically used with HTTP PUT, but not what we are looking for here.
			// 203: indicates a transforming proxy has modified the message (typically used to notify of suspected malware sites). Not the result we are looking for.
			// 204: no content. Technically success, but no data to report. May be used with HTTP DELETE
			// 205: Indicates the browser should reset itself. Not an error, but not what we want here.
			// 206: Partial content. Not an error, but not what we want here so we treat it as one.
			// 207: Only relevant to WebDAV. Treat as an error here.
			// 208: Only relevant to WebDAV in conjunction with 207. Treat as an error here.
			default:
			// No need to handle 3xx status codes as CURLOPT_FOLLOWLOCATION is set TRUE and these should be dealt with internally
			// Codes >= 400 are error codes
				if (verbose) {
					siril_debug_print("Fetch failed with code %ld for URL %s\n", *code, url);
					siril_log_color_message(_("Server unreachable or unresponsive (HTTP code %ld - for details see https://developer.mozilla.org/en-US/docs/Web/HTTP/Status)\n"), "red", *code);
					if (content->data) {
						gchar **lines = g_strsplit(content->data, "\n", 4);
						for (int i = 0; i < 3 && lines[i] != NULL; i++) {
							siril_log_message("%s\n", lines[i]);
						}
						g_strfreev(lines);
					}
					set_progress_bar_data(_("Server unreachable or unresponsive"), 1.0);
				}
		}
	} else {
		siril_log_color_message(_("URL retrieval failed. libcurl error: [%ld]\n"), "red", retval);
	}
	return result;
}

gpointer fetch_url_async(gpointer p) {
	fetch_url_async_data *args = (fetch_url_async_data *) p;
	g_assert(args->idle_function != NULL);
	struct ucontent content = {NULL, 0};
	set_progress_bar_data(NULL, 0.1);
	CURL *curl = initialize_curl(args->url, &content, HTTP_GET, NULL);
	if (!curl) {
		g_free(args->url);
		free(args);
		return NULL;
	}
	char *result = handle_curl_response(curl, &content, args->url, &args->code, args->verbose);
	curl_easy_cleanup(curl);
	g_free(args->url);
	args->url = NULL;
	args->length = content.len;
	args->content = result;
	set_progress_bar_data(NULL, PROGRESS_DONE);
	siril_add_idle(args->idle_function, args);
	return NULL;
}

char* fetch_url(const gchar *url, gsize *length, int *error, gboolean quiet) {
	*error = 0;
	struct ucontent content = {NULL, 0};
	set_progress_bar_data(NULL, 0.1);
	CURL *curl = initialize_curl(url, &content, HTTP_GET, NULL);
	if (!curl) {
		*error = 1;
		set_progress_bar_data(NULL, PROGRESS_DONE);
		return NULL;
	}
	long code;
	char *result = handle_curl_response(curl, &content, url, &code, (!quiet));
	curl_easy_cleanup(curl);
	set_progress_bar_data(NULL, PROGRESS_DONE);
	*length = content.len;
	if (!result || content.len == 0 || code != 200) {
		free(content.data);
		result = NULL;
		*error = 1;
	}
	return result;
}

char* fetch_url_range_with_curl(void* curlp, const gchar *url, size_t start, size_t length,
                                gsize *response_length, int *error, gboolean quiet) {
	*error = 0;
	*response_length = 0;
	struct ucontent content = {NULL, 0};

	if (!curlp) {
		if (!quiet) {
			siril_log_color_message(_("Error: NULL CURL handle provided.\n"), "red");
		}
		*error = 1;
		return NULL;
	}
	CURL *curl = (CURL *) curlp;

	// Construct the range header
	gchar *range_header = g_strdup_printf("%zu-%zu", start, start + length - 1);

	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, &content);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, SIRIL_USER_AGENT);
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	retval |= curl_easy_setopt(curl, CURLOPT_RANGE, range_header);

	g_free(range_header);

	if (retval) {
		if (!quiet) {
			siril_log_color_message(_("Error in curl_easy_setopt() for range request\n"), "red");
		}
		*error = 1;
		return NULL;
	}

	if (g_getenv("CURL_CA_BUNDLE")) {
		if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE"))) {
			if (!quiet) {
				siril_log_color_message(_("Error configuring CURL with CA bundle.\n"), "red");
			}
		}
	}

	content.data = calloc(1, 1);
	if (content.data == NULL) {
		PRINT_ALLOC_ERR;
		*error = 1;
		return NULL;
	}

	CURLcode res = curl_easy_perform(curl);
	char *result = NULL;

	if (res == CURLE_OK) {
		long code;
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

		// Accept both 206 (Partial Content) and 200 (OK)
		if (code == 206 || code == 200) {
			result = content.data;
			*response_length = content.len;
#ifdef NETWORKING_DEBUG
			siril_debug_print("Retrieved result from %s, length %lu\n", url, content.len);
#endif
		} else {
			if (!quiet) {
				siril_log_color_message(_("HTTP range request failed with code %ld for URL %s\n"),
				                       "red", code, url);
			}
			free(content.data);
			*error = 1;
		}
	} else {
		if (!quiet) {
			siril_log_color_message(_("URL range request failed. libcurl error: %s\n"),
			                       "red", curl_easy_strerror(res));
		}
		free(content.data);
		*error = 1;
	}

	return result;
}

char* fetch_url_range(const gchar *url, size_t start, size_t length,
                      gsize *response_length, int *error, gboolean quiet) {
	*error = 0;
	*response_length = 0;
	struct ucontent content = {NULL, 0};

	CURL *curl = curl_easy_init();
	if (!curl) {
		if (!quiet) {
			siril_log_color_message(_("Error initialising CURL handle for range request.\n"), "red");
		}
		*error = 1;
		return NULL;
	}

	// Construct the range header
	gchar *range_header = g_strdup_printf("%zu-%zu", start, start + length - 1);

	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, &content);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, SIRIL_USER_AGENT);
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	retval |= curl_easy_setopt(curl, CURLOPT_RANGE, range_header);

	g_free(range_header);

	if (retval) {
		if (!quiet) {
			siril_log_color_message(_("Error in curl_easy_setopt() for range request\n"), "red");
		}
		curl_easy_cleanup(curl);
		*error = 1;
		return NULL;
	}

	if (g_getenv("CURL_CA_BUNDLE")) {
		if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE"))) {
			if (!quiet) {
				siril_log_color_message(_("Error configuring CURL with CA bundle.\n"), "red");
			}
		}
	}

	content.data = calloc(1, 1);
	if (content.data == NULL) {
		PRINT_ALLOC_ERR;
		curl_easy_cleanup(curl);
		*error = 1;
		return NULL;
	}

	CURLcode res = curl_easy_perform(curl);
	char *result = NULL;

	if (res == CURLE_OK) {
		long code;
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

		// Accept both 206 (Partial Content) and 200 (OK)
		if (code == 206 || code == 200) {
			result = content.data;
			*response_length = content.len;
#ifdef NETWORKING_DEBUG
			siril_debug_print("Retrieved result from %s, length %lu\n", url, content.len);
#endif
		} else {
			if (!quiet) {
				siril_log_color_message(_("HTTP range request failed with code %ld for URL %s\n"),
				                       "red", code, url);
			}
			free(content.data);
			*error = 1;
		}
	} else {
		if (!quiet) {
			siril_log_color_message(_("URL range request failed. libcurl error: %s\n"),
			                       "red", curl_easy_strerror(res));
		}
		free(content.data);
		*error = 1;
	}

	curl_easy_cleanup(curl);
	return result;
}

int submit_post_request(const char *url, const char *post_data, char **post_response) {
	struct ucontent chunk = {NULL, 0};  // will be grown as needed by realloc above
	CURL *curl = initialize_curl(url, &chunk, HTTP_POST, post_data);
	if (!curl) {
		free(chunk.data);
		return 1;
	}
	CURLcode res = curl_easy_perform(curl);
	if (res != CURLE_OK) {
		siril_log_color_message(_("Error fetching URL: %s\n"), "red", curl_easy_strerror(res));
	} else {
		*post_response = g_strdup(chunk.data);
	}
	free(chunk.data);
	curl_easy_cleanup(curl);
	return (res != CURLE_OK ? 1 : 0);
}

// abort callback
static size_t cbk_abort(void *ptr, size_t size, size_t nmemb, void *data) {
	size_t total = size * nmemb;
	
	// If it's exactly 1 byte, we assume it's our requested range.
	// We "consume" it to let the connection close gracefully (CURLE_OK).
	if (total == 1) {
		return total;
	}

	// If it's anything else, it's likely the entire request body so we
	// return zero to kill the connection immediately (CURLE_WRITE_ERROR).
	return 0;
}

/**
 * Checks if a URL is up and supports Range requests.
 * Returns:
 * - >0 (ms) if UP and supports Range (206)
 * - -2      if UP but NO Range support (200)
 * - -1      if DOWN or Error
 */
int http_check(const gchar *url) {
	CURL *curl;
	CURLcode res;
	int result_val = -1;

	if (!url) return -1;

	curl = curl_easy_init();
	if (!curl) return -1;

	CURLcode retval;
	retval = curl_easy_setopt(curl, CURLOPT_URL, url);
	retval |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	retval |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_abort);
	retval |= curl_easy_setopt(curl, CURLOPT_USERAGENT, SIRIL_USER_AGENT);
	retval |= curl_easy_setopt(curl, CURLOPT_RANGE, "0-0");
	retval |= curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	retval |= curl_easy_setopt(curl, CURLOPT_TIMEOUT, 5L);
	retval |= curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);

	if (retval) {
		siril_debug_print("Error in curl_easy_setopt()\n");
		curl_easy_cleanup(curl);
		return -1;
	}

	res = curl_easy_perform(curl);

	if (res == CURLE_OK || res == CURLE_WRITE_ERROR) {
		long response_code;
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);

		if (response_code == 206) {
			// Success: Server is up and handled the range request
			double total_time = 0;
			curl_easy_getinfo(curl, CURLINFO_TOTAL_TIME, &total_time);
			result_val = (int)(total_time * 1000.0);
		}
		else if (response_code == 200) {
			// Server is up, but ignored the Range header
			result_val = -2;
		}
	}

	curl_easy_cleanup(curl);
	return result_val;
}

gboolean siril_compiled_with_networking() {
	return TRUE;
}

#else

// Stub functions when libcurl is not available
// These return errors but allow the code to compile and link

gpointer fetch_url_async(gpointer p) {
	fetch_url_async_data *args = (fetch_url_async_data *) p;
	g_assert(args->idle_function != NULL);

	siril_log_color_message(_("Error: Siril was compiled without libcurl support. Network features are unavailable.\n"), "red");

	// Clean up and call the idle function with NULL content to signal failure
	g_free(args->url);
	args->url = NULL;
	args->length = 0;
	args->content = NULL;
	args->code = 0;

	siril_add_idle(args->idle_function, args);
	return NULL;
}

char* fetch_url(const gchar *url, gsize *length, int *error, gboolean quiet) {
	if (!quiet) {
		siril_log_color_message(_("Error: Siril was compiled without libcurl support. Cannot fetch URL: %s\n"), "red", url);
	}
	*error = 1;
	*length = 0;
	return NULL;
}

char* fetch_url_range_with_curl(void* curlp, const gchar *url, size_t start, size_t length,
                                gsize *response_length, int *error, gboolean quiet) {
	if (!quiet) {
		siril_log_color_message(_("Error: Siril was compiled without libcurl support. Cannot fetch URL range: %s\n"),
		                       "red", url);
	}
	*error = 1;
	*response_length = 0;
	return NULL;
}


char* fetch_url_range(const gchar *url, size_t start, size_t length,
                      gsize *response_length, int *error, gboolean quiet) {
	if (!quiet) {
		siril_log_color_message(_("Error: Siril was compiled without libcurl support. Cannot fetch URL range: %s\n"),
		                       "red", url);
	}
	*error = 1;
	*response_length = 0;
	return NULL;
}

int submit_post_request(const char *url, const char *post_data, char **post_response) {
	siril_log_color_message(_("Error: Siril was compiled without libcurl support. Cannot submit POST request to: %s\n"), "red", url);
	*post_response = NULL;
	return 1;
}

gboolean siril_compiled_with_networking() {
	return FALSE;
}

int http_check(const gchar *url) {
	return -1;
}

#endif

gboolean is_online() {
	return (siril_compiled_with_networking() && online_status);
}

gboolean set_online_status(gboolean status) {
	online_status = status;
	siril_log_message(online_status ? N_("Siril is in online mode.\n") : N_("Siril is in offline mode.\n"));
	if (!siril_compiled_with_networking())
		return FALSE;
	return online_status;
}
