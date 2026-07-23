#include <stdio.h>
#include <sys/file.h>
#include <unistd.h>
#include "core/siril.h"
#include "core/proto.h"
#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#endif
#include "download_files.h"

#define TEST_FILES_BUCKET_URL "https://siril-share-public.s3.rbx.io.cloud.ovh.net"
#define TEST_FILES_DOWNLOAD_DIR "siril/tests" // in user cache dir
#define LOCK_FILE_FOR_DOWNLOAD "/tmp/siril-test-download.lock"

#define LOCK(retval) \
	int fd_lock = open(LOCK_FILE_FOR_DOWNLOAD, O_RDWR | O_CREAT, 0640); \
	if (fd_lock < 0) { \
		siril_log_debug("could not create lock file\n"); \
		return retval; \
	} \
	flock(fd_lock, LOCK_EX)

#define UNLOCK \
	flock(fd_lock, LOCK_UN); \
	close(fd_lock);

static size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
	size_t written = fwrite(ptr, size, nmemb, stream);
	return written;
}

static gboolean directory_created = FALSE;

/* each test runs in a separate process, we cannot use a mutex to protect check_or_download_test_file(),
 * we use a file lock instead, which may not be portable. The other solution would be to make a function
 * with all the test function and call this one instead in the single Test() */
void init_download() {
	LOCK();
	const gchar *cache_dir_path = g_get_user_cache_dir();
	gchar *test_files_dir_path = g_strdup_printf("%s/%s", cache_dir_path, TEST_FILES_DOWNLOAD_DIR);
	if (g_mkdir_with_parents(test_files_dir_path, 0700)) {
		siril_log_debug("failed to create dir %s\n", test_files_dir_path);
	} else {
		directory_created = TRUE;
	}
	g_free(test_files_dir_path);
	UNLOCK;
}

gchar *check_or_download_test_file(const char *filename) {
	LOCK(NULL);
	const gchar *cache_dir_path = g_get_user_cache_dir();
	gchar *file_path = g_strdup_printf("%s/%s/%s", cache_dir_path, TEST_FILES_DOWNLOAD_DIR, filename);
	if (is_readable_file(file_path)) {
		UNLOCK;
		return file_path;
	}
	if (!directory_created) {
		siril_log_debug("cannot download to expected directory which creation failed\n");
		g_free(file_path);
		UNLOCK;
		return NULL;
	}
	if (string_is_a_path(filename)) {
		gchar *dir = g_path_get_dirname(file_path);
		if (g_mkdir_with_parents(dir, 0700)) {
			siril_log_debug("failed to create dir %s\n", dir);
			g_free(dir);
			g_free(file_path);
			UNLOCK;
			return NULL;
		}
		g_free(dir);
	}
#ifdef HAVE_LIBCURL
	siril_log_debug("test file %s is missing, downloading...\n", file_path);
	FILE *dest_file = g_fopen(file_path, "wb");
	if (!dest_file) {
		perror("fopen");
		siril_log_debug("failed to create file %s for downloads\n", file_path);
		UNLOCK;
		return NULL;
	}
	gchar *file_url = g_strdup_printf("%s/%s", TEST_FILES_BUCKET_URL, filename);

	CURL *curl;
	curl = curl_easy_init();
	if (curl) {
		curl_easy_setopt(curl, CURLOPT_URL, file_url);
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, dest_file);
		curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L);
		CURLcode res = curl_easy_perform(curl);
		curl_easy_cleanup(curl);
		if (res != CURLE_OK) {
			siril_log_debug("failed to download from %s\n", file_url);
			g_unlink(file_path);
			g_free(file_path);
			file_path = NULL;
		}
	}
	fclose(dest_file);
	g_free(file_url);
	//siril_log_debug("test file %s was downloaded\n", file_path);
#else
	siril_log_debug("test file %s is missing and cannot be doawnloaded (libcurl missing)\n", file_path);
	g_free(file_path);
	file_path = NULL;
#endif
	UNLOCK;
	return file_path;
}

