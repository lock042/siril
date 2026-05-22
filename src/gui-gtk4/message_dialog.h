#ifndef SRC_MESSAGE_DIALOG_H_
#include <gtk/gtk.h>
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
void queue_message_dialog(GtkMessageType type, const char *title, const char *text);
void queue_error_message_dialog(const char *title, const char *text);
void queue_warning_message_dialog(const char *title, const char *text);
void siril_data_dialog(GtkMessageType type, char *title, char *text, gchar *data);
gboolean siril_confirm_dialog(gchar *title, gchar *msg, gchar *button_accept);
gboolean siril_confirm_dialog_and_remember(gchar *title, gchar *msg, gchar *button_accept, gboolean *user_data);
gboolean siril_confirm_data_dialog(GtkMessageType type, char *title, char *text, gchar *button_accept, gchar *data);
gboolean siril_confirm_dialog_async(gchar *title, gchar *msg, gchar *button_accept);

/* Three-button alert: 0 = Cancel/dismissed, 1 = button_a, 2 = button_b.
 * Default button is B (preserve-data convention); Escape / close → Cancel. */
int siril_three_button_dialog(gchar *title, gchar *msg,
                              gchar *button_a, gchar *button_b);
#endif
