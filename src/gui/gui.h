#ifndef _GUI_H
#define _GUI_H

#include <gtk/gtk.h>
#include "core/siril.h"

GtkWidget* lookup_widget (const gchar *widget_name);

void init_gui(enum _siril_mode new_mode, char *start_cwd);

void load_prefered_theme();

void show_supported_files(gchar *supported_files);

#endif
