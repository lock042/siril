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
/* Same as siril_confirm_dialog, but also embeds an AVI Bayer-pattern
 * combo in the dialog content area. On Accept, fills *avi_bayer_pattern
 * with an `enum mpp_avi_bayer` value (0..5; 0 = Auto). Used by the
 * deprecated-AVI dialog so the user can pick the mosaic pattern before
 * convert_single_film_to_ser runs. */
gboolean siril_confirm_dialog_with_avi_bayer(gchar *title, gchar *msg,
		gchar *button_accept, int *avi_bayer_pattern);
#endif
