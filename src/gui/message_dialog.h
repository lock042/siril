#ifndef SRC_MESSAGE_DIALOG_H_
#define SRC_MESSAGE_DIALOG_H_

#include "core/siril.h"

struct siril_dialog_data {
	GtkWindow *parent;
	GtkMessageType type;
	const char *primary_text;
	const char *secondary_text;
};

void siril_message_dialog(GtkMessageType type, char *title, char *text);
gboolean siril_confirm_dialog(gchar *title, gchar *msg, gboolean check_button);

#endif
