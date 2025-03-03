#ifndef _LIVESTACK_GUI_
#define _LIVESTACK_GUI_

#include <glib.h>
#include "registration/registration.h"

void show_hide_toolbox();
void force_unlinked_channels();
void livestacking_display(gchar *str, gboolean free_after_display);
void livestacking_display_config(gboolean use_dark, gboolean use_flat, transformation_type regtype);
void livestacking_update_number_of_images(int nb, double total_exposure, double noise, const char *process_time);

void update_debayer_button_status(gboolean new_state);
gboolean livestacking_first_result_idle(gpointer p);

void enable_debayer(gboolean arg);
gboolean end_image_loading(gpointer arg);

#endif
