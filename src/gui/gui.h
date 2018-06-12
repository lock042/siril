#ifndef _GUI_H
#define _GUI_H

#include <gtk/gtk.h>
#include "core/siril.h"

GtkWidget* lookup_widget (const gchar *widget_name);

void init_gui(gchar *startup_dir, enum _siril_mode new_mode);

#endif
