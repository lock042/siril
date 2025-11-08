#ifndef _PROGRESSLOG_H
#define _PROGRESSLOG_H

#include <sys/time.h>
#include <gtk/gtk.h>

#define PROGRESS_NONE -2.0		// don't update the progress bar value
#define PROGRESS_PULSATE -1.0		// pulsate the progress bar
#define PROGRESS_RESET 0.0		// reset the progress bar
#define PROGRESS_DONE 1.0		// fill the progress bar
#define PROGRESS_TEXT_RESET ""		// reset the progress bar's text

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
