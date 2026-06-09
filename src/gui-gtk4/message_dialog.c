/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* Phase 14G.1: GtkMessageDialog removed in GTK 4.10.
 *
 * This module presents a synchronous-style API to the rest of Siril:
 *   siril_message_dialog (), siril_confirm_dialog (), …
 * but the GTK4 replacements (GtkAlertDialog and our custom GtkWindow)
 * are inherently asynchronous.  We bridge the gap by spinning a nested
 * GMainLoop while the dialog is shown — the loop quits when the user
 * picks a button or closes the window.
 *
 * Routing rules:
 *   plain message  / question with no extras → GtkAlertDialog
 *   message + scrollable text body          → custom GtkWindow
 *   confirm + "do not show again" checkbox  → custom GtkWindow
 */

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"

#include "message_dialog.h"

static GtkWindow *msg_control_window = NULL;

static void message_dialog_init_statics(void) {
	if (msg_control_window) return;
	msg_control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
}

gchar *strip_last_ret_char(gchar *str) {
	char *pch;
	int len;

	pch = strrchr(str, '\0');
	len = pch - str;

	if (str[len - 1] == '\n') {
		str[len - 1] = '\0';
	}
	return str;
}

/* ── Sync wrappers around the async GtkAlertDialog API ───────────────────── */

struct alert_sync_state {
	GMainLoop *loop;
	int        response;     /* button index, or -1 on error */
};

static void alert_sync_choose_cb(GObject *src, GAsyncResult *res, gpointer ud) {
	struct alert_sync_state *st = (struct alert_sync_state *)ud;
	GError *err = NULL;
	st->response = gtk_alert_dialog_choose_finish(GTK_ALERT_DIALOG(src), res, &err);
	if (err) {
		st->response = -1;
		g_clear_error(&err);
	}
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
}

static int alert_choose_sync(GtkAlertDialog *ad, GtkWindow *parent) {
	struct alert_sync_state st = { g_main_loop_new(NULL, FALSE), -1 };
	gtk_alert_dialog_choose(ad, parent, NULL, alert_sync_choose_cb, &st);
	g_main_loop_run(st.loop);
	g_main_loop_unref(st.loop);
	return st.response;
}

/* ── Custom GtkWindow used for dialogs that have extra content (text view
 *    or "do not show again" checkbox) ─────────────────────────────────── */

struct custom_dlg_state {
	GMainLoop *loop;
	int        result;       /* 0 = cancel, 1 = accept */
	GtkWidget *checkbtn;     /* may be NULL */
	gboolean   check_value;
};

static void custom_dlg_accept_cb(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct custom_dlg_state *st = (struct custom_dlg_state *)ud;
	st->result = 1;
	if (st->checkbtn)
		st->check_value = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(st->checkbtn)));
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
}

static void custom_dlg_cancel_cb(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct custom_dlg_state *st = (struct custom_dlg_state *)ud;
	st->result = 0;
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
}

static gboolean custom_dlg_close_request_cb(GtkWindow *win, gpointer ud) {
	(void)win;
	struct custom_dlg_state *st = (struct custom_dlg_state *)ud;
	st->result = 0;
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
	return FALSE;  /* allow default close handling */
}

/* Build and run a custom modal GtkWindow.  Returns 1 if user accepted,
 * 0 otherwise.  If checkbutton is TRUE *check_out_opt receives the
 * checkbox value on accept. */
static int run_custom_message_window(GtkWindow *parent,
		const char *title_text,
		const char *secondary_markup,
		const char *data_text_opt,
		const char *button_accept,
		gboolean    show_checkbutton,
		gboolean   *check_out_opt) {
	struct custom_dlg_state st = { g_main_loop_new(NULL, FALSE), 0, NULL, FALSE };

	GtkWidget *win = gtk_window_new();
	gtk_window_set_modal(GTK_WINDOW(win), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(win), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(win), parent);
	gtk_window_set_title(GTK_WINDOW(win), title_text ? title_text : "");
	gtk_window_set_default_size(GTK_WINDOW(win), 480, -1);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12);
	gtk_widget_set_margin_start(vbox, 18);
	gtk_widget_set_margin_end(vbox, 18);
	gtk_widget_set_margin_top(vbox, 18);
	gtk_widget_set_margin_bottom(vbox, 18);
	gtk_window_set_child(GTK_WINDOW(win), vbox);

	GtkWidget *primary = gtk_label_new(NULL);
	gchar *primary_markup = g_markup_printf_escaped("<big><b>%s</b></big>",
			title_text ? title_text : "");
	gtk_label_set_markup(GTK_LABEL(primary), primary_markup);
	g_free(primary_markup);
	gtk_label_set_xalign(GTK_LABEL(primary), 0);
	gtk_label_set_wrap(GTK_LABEL(primary), TRUE);
	gtk_box_append(GTK_BOX(vbox), primary);

	if (secondary_markup) {
		GtkWidget *sec = gtk_label_new(NULL);
		gtk_label_set_markup(GTK_LABEL(sec), secondary_markup);
		gtk_label_set_xalign(GTK_LABEL(sec), 0);
		gtk_label_set_wrap(GTK_LABEL(sec), TRUE);
		gtk_box_append(GTK_BOX(vbox), sec);
	}

	if (data_text_opt) {
		GtkWidget *tview = gtk_text_view_new();
		GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));
		gtk_text_buffer_set_text(buf, data_text_opt, -1);
		gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(tview), FALSE);
		gtk_text_view_set_editable(GTK_TEXT_VIEW(tview), FALSE);
		gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(tview), GTK_WRAP_WORD);
		gtk_widget_set_margin_start(tview, 6);

		GtkWidget *swin = gtk_scrolled_window_new();
		gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(swin), tview);
		gtk_widget_set_size_request(swin, -1, 200);
		gtk_widget_set_hexpand(swin, TRUE);
		gtk_widget_set_vexpand(swin, TRUE);
		gtk_box_append(GTK_BOX(vbox), swin);
	}

	if (show_checkbutton) {
		st.checkbtn = gtk_check_button_new_with_mnemonic(_("_Do not show this dialog again"));
		gtk_box_append(GTK_BOX(vbox), st.checkbtn);
	}

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	gtk_box_append(GTK_BOX(vbox), bbox);

	GtkWidget *btn_cancel = gtk_button_new_with_mnemonic(_("_Cancel"));
	gtk_box_append(GTK_BOX(bbox), btn_cancel);
	g_signal_connect(btn_cancel, "clicked", G_CALLBACK(custom_dlg_cancel_cb), &st);

	GtkWidget *btn_accept = gtk_button_new_with_mnemonic(button_accept ? button_accept : _("_OK"));
	gtk_widget_add_css_class(btn_accept, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), btn_accept);
	g_signal_connect(btn_accept, "clicked", G_CALLBACK(custom_dlg_accept_cb), &st);
	gtk_window_set_default_widget(GTK_WINDOW(win), btn_accept);

	g_signal_connect(win, "close-request", G_CALLBACK(custom_dlg_close_request_cb), &st);

	gtk_window_present(GTK_WINDOW(win));
	g_main_loop_run(st.loop);
	g_main_loop_unref(st.loop);

	if (check_out_opt)
		*check_out_opt = st.check_value;

	gtk_window_destroy(GTK_WINDOW(win));
	return st.result;
}

/* ── Internal sync runners that pick GtkAlertDialog vs custom GtkWindow ──── */

static gboolean siril_confirm_data_dialog_internal(gpointer p, gchar *button_accept) {
	struct siril_dialog_data *args = (struct siril_dialog_data *) p;
	GtkWindow *parent = args->parent;
	if (!GTK_IS_WINDOW(parent)) {
		message_dialog_init_statics();
		parent = msg_control_window;
	}
	int r = run_custom_message_window(parent, args->primary_text,
			args->secondary_text, args->data, button_accept,
			FALSE, NULL);
	return r == 1;
}

gboolean siril_confirm_data_dialog(GtkMessageType type, char *title, char *text, gchar *button_accept, gchar *data) {
	(void)type;  /* GtkAlertDialog/GtkWindow surrogate has no message-type icon */
	if (com.headless || com.script)
		return FALSE;
	struct siril_dialog_data *args = calloc(1, sizeof(struct siril_dialog_data));

	args->parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(args->parent)) {
		message_dialog_init_statics();
		args->parent = msg_control_window;
	}
	strip_last_ret_char(title);
	strip_last_ret_char(text);

	args->primary_text = title;
	args->secondary_text = text;
	args->data = data;
	args->type = type;
	gboolean retval = siril_confirm_data_dialog_internal(args, button_accept);
	free(args);
	return retval;
}

static gboolean show_modal_dialog(gpointer p) {
	struct siril_dialog_data *args = (struct siril_dialog_data *) p;
	GtkWindow *parent = args->parent;
	if (!GTK_IS_WINDOW(parent)) {
		message_dialog_init_statics();
		parent = msg_control_window;
	}

	if (args->data) {
		/* Information-style dialog with extra text body — custom window. */
		run_custom_message_window(parent, args->primary_text,
				args->secondary_text, args->data, _("_OK"),
				FALSE, NULL);
	} else {
		/* Plain message — use GtkAlertDialog with a single OK button. */
		GtkAlertDialog *ad = gtk_alert_dialog_new("%s",
				args->primary_text ? args->primary_text : "");
		if (args->secondary_text)
			gtk_alert_dialog_set_detail(ad, args->secondary_text);
		const char *btns[] = { _("_OK"), NULL };
		gtk_alert_dialog_set_buttons(ad, btns);
		gtk_alert_dialog_set_default_button(ad, 0);
		gtk_alert_dialog_set_cancel_button (ad, 0);
		alert_choose_sync(ad, parent);
		g_object_unref(ad);
	}
	g_free(args);
	return FALSE;
}

/******* Public functions *****************/

void siril_message_dialog(GtkMessageType type, char *title, char *text) {
	(void)type;
	if (com.headless || com.script)
		return;
	struct siril_dialog_data *args = calloc(1, sizeof(struct siril_dialog_data));

	args->parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(args->parent)) {
		message_dialog_init_statics();
		args->parent = msg_control_window;
	}
	strip_last_ret_char(title);
	strip_last_ret_char(text);

	args->primary_text = title;
	args->secondary_text = text;
	args->data = NULL;
	args->type = type;
	show_modal_dialog(args);
}

gboolean siril_message_dialog_idle(gpointer p) {
	struct message_data *data = (struct message_data *) p;
	siril_message_dialog(data->type, data->title, data->text);
	free(data->title);
	free(data->text);
	free(data);
	return FALSE;
}

void queue_message_dialog(GtkMessageType type, const char *title, const char *text) {
	if (com.headless || com.script)
		return;
	struct message_data *data = calloc(1, sizeof(struct message_data));
	data->type = type;
	data->title = strdup(title);
	data->text = strdup(text);
	g_idle_add(siril_message_dialog_idle, data);
}

void queue_error_message_dialog(const char *title, const char *text) {
	queue_message_dialog(GTK_MESSAGE_ERROR, title, text);
}

void queue_warning_message_dialog(const char *title, const char *text) {
	queue_message_dialog(GTK_MESSAGE_WARNING, title, text);
}

void siril_data_dialog(GtkMessageType type, char *title, char *text, gchar *data) {
	(void)type;
	if (com.headless || com.script)
		return;
	struct siril_dialog_data *args = calloc(1, sizeof(struct siril_dialog_data));

	args->parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(args->parent)) {
		message_dialog_init_statics();
		args->parent = msg_control_window;
	}
	strip_last_ret_char(title);
	strip_last_ret_char(text);

	args->primary_text = title;
	args->secondary_text = text;
	args->data = data;
	args->type = type;
	show_modal_dialog(args);
}

static gboolean siril_confirm_dialog_internal(gchar *title, gchar *msg, gchar *button_accept, gboolean checkbutton, gboolean *user_data) {
	GtkWindow *parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(parent)) {
		message_dialog_init_statics();
		parent = msg_control_window;
	}

	strip_last_ret_char(title);
	strip_last_ret_char(msg);

	if (checkbutton) {
		/* GtkAlertDialog has no extra-widget API, so the
		 * "do not show again" variant uses a custom window. */
		gboolean check_value = FALSE;
		int r = run_custom_message_window(parent, title, msg, NULL,
				button_accept, TRUE, &check_value);
		if (user_data)
			*user_data = (r == 1) ? check_value : FALSE;
		return r == 1;
	}

	GtkAlertDialog *ad = gtk_alert_dialog_new("%s", title ? title : "");
	if (msg)
		gtk_alert_dialog_set_detail(ad, msg);
	const char *btns[] = { _("_Cancel"), button_accept ? button_accept : _("_OK"), NULL };
	gtk_alert_dialog_set_buttons(ad, btns);
	gtk_alert_dialog_set_default_button(ad, 1);
	gtk_alert_dialog_set_cancel_button (ad, 0);

	int r = alert_choose_sync(ad, parent);
	g_object_unref(ad);

	if (user_data)
		*user_data = FALSE;
	return r == 1;
}

gboolean siril_confirm_dialog(gchar *title, gchar *msg, gchar *button_accept) {
	return siril_confirm_dialog_internal(title, msg, button_accept, FALSE, NULL);
}

struct confirm_dialog_data {
    gchar *title;
    gchar *msg;
    gchar *button_accept;
    gboolean result;
};

static gboolean confirm_dialog_idle(gpointer arg) {
    struct confirm_dialog_data *data = (struct confirm_dialog_data *)arg;
    data->result = siril_confirm_dialog(data->title, data->msg, data->button_accept);
    return FALSE;
}

gboolean siril_confirm_dialog_async(gchar *title, gchar *msg, gchar *button_accept) {
    struct confirm_dialog_data data;
    data.title = title;
    data.msg = msg;
    data.button_accept = button_accept;
    data.result = FALSE;

    execute_idle_and_wait_for_it(confirm_dialog_idle, &data);

    return data.result;
}

gboolean siril_confirm_dialog_and_remember(gchar *title, gchar *msg, gchar *button_accept, gboolean *user_data) {
	return siril_confirm_dialog_internal(title, msg, button_accept, TRUE, user_data);
}
