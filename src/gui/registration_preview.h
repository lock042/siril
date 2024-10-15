#ifndef _REGISTRATION_PREVIEW_H
#define _REGISTRATION_PREVIEW_H

#include <gtk/gtk.h>

gboolean redraw_previews(gpointer user_data);
void clear_previews();
void set_preview_area(int preview_area, int centerX, int centerY);
void adjust_reginfo();
void on_spinbut_shift_value_change(GtkSpinButton *spinbutton,
		gpointer user_data);
void test_and_allocate_reference_image(int vport);
void enable_view_reference_checkbox(gboolean status);

#endif
