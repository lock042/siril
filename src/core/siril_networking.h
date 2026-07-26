#ifndef SIRIL_NETWORKING_H
#define SIRIL_NETWORKING_H

#include <glib.h>
struct ucontent {
	char *data;
	size_t len;
};

typedef struct _fetch_url_async_data {
	gchar *url;
	gchar *content;
	gsize length;
	gboolean verbose;
	long code;
	gboolean (*idle_function)(gpointer args);
} fetch_url_async_data;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

gpointer fetch_url_async(gpointer p);
char *fetch_url(const gchar *url, gsize *length, int *error, gboolean quiet);
char* fetch_url_range(const gchar *url, size_t start, size_t length,
                      gsize *response_length, int *error, gboolean quiet);
char* fetch_url_range_with_curl(void* curlp, const gchar *url, size_t start, size_t length,
                                gsize *response_length, int *error, gboolean quiet);

int submit_post_request(const char *url, const char *post_data, char **post_response);
int http_check(const gchar *url);

gboolean siril_compiled_with_networking();
gboolean is_online();
gboolean set_online_status(gboolean status);
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // SIRIL_NETWORKING_H
