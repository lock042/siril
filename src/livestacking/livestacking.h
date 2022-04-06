#ifndef _LIVESTACK_H
#define _LIVESTACK_H

#include <glib.h>

void on_livestacking_start();
void stop_live_stacking_engine();
void pause_live_stacking_engine();

int get_paused_status();

int start_livestack_from_command(gchar *dark, gchar *flat, gboolean use_file_watcher, gboolean remove_gradient);
void init_preprocessing_finalize();
void start_livestacking(gboolean with_filewatcher);
gboolean livestacking_is_started();
gboolean livestacking_uses_filewatcher();

void livestacking_queue_file(char *file);

#endif
