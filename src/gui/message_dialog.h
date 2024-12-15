#ifndef SRC_MESSAGE_DIALOG_H_
#define SRC_MESSAGE_DIALOG_H_

#include "core/siril.h"

struct siril_dialog_data {
	GtkWindow *parent;
	GtkMessageType type;
	gchar *data;
	const char *primary_text;
	const char *secondary_text;
	gpointer user_data;
};

struct message_data {
	GtkMessageType type;
	char *title;
	char *text;
};

gboolean siril_message_dialog_idle(gpointer p);
void siril_message_dialog(GtkMessageType type, char *title, char *text);
void queue_error_message_dialog(const char *title, const char *text);
void queue_warning_message_dialog(const char *title, const char *text);
void siril_data_dialog(GtkMessageType type, char *title, char *text, gchar *data);
gboolean siril_confirm_dialog(gchar *title, gchar *msg, gchar *button_accept);
gboolean siril_confirm_dialog_and_remember(gchar *title, gchar *msg, gchar *button_accept, gboolean *user_data);
gboolean siril_confirm_data_dialog(GtkMessageType type, char *title, char *text, gchar *button_accept, gchar *data);
#endif
