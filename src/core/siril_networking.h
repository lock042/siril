#ifndef SIRIL_NETWORKING_H
#define SIRIL_NETWORKING_H

#if defined(HAVE_LIBCURL)
#include <curl/curl.h>

struct ucontent {
	char *data;
	size_t len;
};

typedef struct _fetch_url_async_data {
	CURL *curl;
	gchar *url;
	gchar *content;
	gboolean verbose;
	gboolean (*idle_function)(gpointer args);
} fetch_url_async_data;

gpointer fetch_url_async(gpointer p);
char *fetch_url(const gchar *url, gsize *length);
int submit_post_request(const char *url, const char *post_data, char **post_response);

// HAVE_LIBCURL
#endif

gboolean siril_compiled_with_networking();
gboolean is_online();
gboolean set_online_status(gboolean status);

// SIRIL_NETWORKING_H
#endif
