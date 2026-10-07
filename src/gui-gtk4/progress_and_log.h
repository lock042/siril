#ifndef _PROGRESSLOG_H
#define _PROGRESSLOG_H

#include <sys/time.h>
#include <glib.h>
/* PROGRESS_* constants are defined in core/gui_iface.h and re-exported here
 * so existing callers of this header continue to compile unchanged. */
#include "core/gui_iface.h"

#ifdef __cplusplus
extern "C" {
#endif

void initialize_log_tags();
void gui_log_message(const char* msg, const char* color);

void set_progress_bar_data(const char *text, double percent);
void set_cursor_waiting(gboolean waiting);
void set_cursor(const gchar* cursor_name);

gchar *get_log_as_string();

#ifdef __cplusplus
}
#endif

#endif
