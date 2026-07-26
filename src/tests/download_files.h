#ifndef _TEST_DOWNLOADS_H
#define  _TEST_DOWNLOADS_H

#include <glib.h>

void init_download();
gchar *check_or_download_test_file(const char *filename);

#endif
