#ifndef _LIVESTACK_GUI_
#define _LIVESTACK_GUI_

#include <glib.h>

void livestacking_display(const char *str);
void update_debayer_button_status(gboolean new_state);
gboolean livestacking_first_result_idle(gpointer p);

#endif
