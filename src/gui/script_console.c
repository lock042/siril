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

/*
 * Script console GUI — status bar logging, command history, end-of-script
 * cleanup, and the interactive command-entry widget handlers.
 *
 * Functions moved here from core/command_line_processor.c so all GTK code
 * is confined to GUI-build translation units.
 */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/command_line_processor.h"
#include "core/command_list.h"
#include "gui/gui_state.h"
#include "gui/image_interactions.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "script_console.h"

/* ── Status-bar logging ──────────────────────────────────────────────────── */

static void update_log_icon(gboolean is_running) {
	GtkImage *image = GTK_IMAGE(GTK_WIDGET(gtk_builder_get_object(gui.builder, "image_log")));
	if (is_running)
		gtk_image_set_from_icon_name(image, "gtk-yes", GTK_ICON_SIZE_LARGE_TOOLBAR);
	else
		gtk_image_set_from_icon_name(image, "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);
}

struct log_status_bar_idle_data {
	gchar *myline;
	int line;
};

static gboolean log_status_bar_idle_callback(gpointer p) {
	struct log_status_bar_idle_data *data = (struct log_status_bar_idle_data *) p;
	GtkStatusbar *statusbar = GTK_STATUSBAR(GTK_WIDGET(gtk_builder_get_object(gui.builder, "statusbar_script")));
	gchar *newline = g_strdup(data->myline);
	gchar *status = g_strdup_printf(_("Processing line %d: %s"), data->line, newline);
	update_log_icon(TRUE);
	gtk_statusbar_push(statusbar, 0, status);
	g_free(newline);
	g_free(status);
	g_free(data->myline);
	free(data);
	return FALSE;
}

/* Called from gui_iface.console_set_status() on any thread. */
void console_log_status(const gchar *msg, int line) {
	if (com.headless)
		return;
	struct log_status_bar_idle_data *data = malloc(sizeof(*data));
	data->line = line;
	data->myline = msg ? g_strdup(msg) : NULL;
	gdk_threads_add_idle(log_status_bar_idle_callback, data);
}

void console_clear_status_bar(void) {
	GtkStatusbar *bar = GTK_STATUSBAR(GTK_WIDGET(gtk_builder_get_object(gui.builder, "statusbar_script")));
	gtk_statusbar_remove_all(bar, 0);
	update_log_icon(FALSE);
}

/* ── End-of-script cleanup ───────────────────────────────────────────────── */

gboolean end_script(gpointer p) {
	console_clear_status_bar();
	gui_iface.set_gui_cwd();
	gui_iface.update_menu_state();
	gfit_modified_update_gui();
	gui_iface.redraw_previews();
	update_zoom_label();
	update_display_fwhm();
	display_filename();
	gui_iface.new_selection_zone();
	gui_iface.update_spin_cpu();
	gui_iface.set_busy(FALSE);
	return FALSE;
}

/* ── Command history ─────────────────────────────────────────────────────── */

static void history_add_line(char *line) {
	if (!gui.cmd_history) {
		gui.cmd_hist_size = CMD_HISTORY_SIZE;
		gui.cmd_history = calloc(gui.cmd_hist_size, sizeof(const char*));
		gui.cmd_hist_current = 0;
		gui.cmd_hist_display = 0;
	}
	gui.cmd_history[gui.cmd_hist_current] = line;
	gui.cmd_hist_current++;
	if (gui.cmd_hist_current == gui.cmd_hist_size)
		gui.cmd_hist_current = 0;
	if (gui.cmd_history[gui.cmd_hist_current]) {
		free(gui.cmd_history[gui.cmd_hist_current]);
		gui.cmd_history[gui.cmd_hist_current] = NULL;
	}
	gui.cmd_hist_display = gui.cmd_hist_current;
}

/* ── Interactive console callbacks ───────────────────────────────────────── */

void on_command_activate(GtkEntry *entry, gpointer user_data) {
	const gchar *text = gtk_entry_get_text(entry);
	history_add_line(strdup(text));
	if (!(processcommand(text, FALSE))) {
		gtk_entry_set_text(entry, "");
		gui_iface.on_precision_changed();
	}
}

#if GTK_CHECK_VERSION(3, 24, 24)
static gboolean on_command_key_press_event(GtkEventController *controller,
		guint keyval, guint keycode, GdkModifierType modifiers,
		GtkWidget *widget) {
#else
static gboolean on_command_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	guint keyval = event->keyval;
#endif
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (keyval) {
	case GDK_KEY_Up:
		handled = 1;
		if (!gui.cmd_history)
			break;
		if (gui.cmd_hist_display > 0) {
			if (gui.cmd_history[gui.cmd_hist_display - 1])
				--gui.cmd_hist_display;
			gtk_entry_set_text(entry, gui.cmd_history[gui.cmd_hist_display]);
		} else if (gui.cmd_history[gui.cmd_hist_size - 1]) {
			gui.cmd_hist_display = gui.cmd_hist_size - 1;
			gtk_entry_set_text(entry, gui.cmd_history[gui.cmd_hist_display]);
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Down:
		handled = 1;
		if (!gui.cmd_history)
			break;
		if (gui.cmd_hist_display == gui.cmd_hist_current)
			break;
		if (gui.cmd_hist_display == gui.cmd_hist_size - 1) {
			if (gui.cmd_hist_current == 0) {
				gui.cmd_hist_display = 0;
			} else if (gui.cmd_history[0]) {
				gui.cmd_hist_display = 0;
			}
		} else {
			gui.cmd_hist_display++;
		}
		if (gui.cmd_history[gui.cmd_hist_display]) {
			gtk_entry_set_text(entry, gui.cmd_history[gui.cmd_hist_display]);
		} else {
			gtk_entry_set_text(entry, "");
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	default:
		break;
	}
	return (handled != 0);
}

/* ── Command-line completion and help popup ──────────────────────────────── */

#define COMPLETION_COLUMN 0

static gboolean on_match_selected(GtkEntryCompletion *completion, GtkTreeModel *model,
		GtkTreeIter *iter, gpointer user_data) {
	const gchar *cmd;
	GtkEntry *entry = GTK_ENTRY(gtk_entry_completion_get_entry(completion));
	gtk_tree_model_get(model, iter, COMPLETION_COLUMN, &cmd, -1);
	/* Block the completion's own "changed" handler while we update the text to
	 * prevent it from queuing a popup-show idle (which would reopen the popup
	 * immediately after GtkEntryCompletion hides it post-selection). */
	g_signal_handlers_block_matched(entry, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, completion);
	gtk_entry_set_text(entry, cmd);
	gtk_editable_set_position(GTK_EDITABLE(entry), -1);
	g_signal_handlers_unblock_matched(entry, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, completion);
	return TRUE;
}

static gboolean completion_match_func(GtkEntryCompletion *completion,
		const gchar *key, GtkTreeIter *iter, gpointer user_data) {
	if (*key == '\0') return FALSE;
	gchar *item = NULL;
	GtkTreeModel *model = gtk_entry_completion_get_model(completion);
	gtk_tree_model_get(model, iter, COMPLETION_COLUMN, &item, -1);
	if (item == NULL) return FALSE;
	gboolean ret = g_str_has_prefix(item, key);
	g_free(item);
	return ret;
}

gboolean show_command_help_popup(gpointer user_data) {
	GtkEntry *entry = (GtkEntry *) user_data;
	gchar *helper = NULL;
	const gchar *text = gtk_entry_get_text(entry);
	if (*text == '\0') {
		helper = g_strdup(_("Please enter an existing command before hitting this button"));
	} else {
		command *current = commands;
		gchar **command_line = g_strsplit_set(text, " ", -1);
		while (current->process) {
			if (!g_ascii_strcasecmp(current->name, command_line[0])) {
				gchar **token = g_strsplit_set(current->usage, " \n", -1);
				GString *str = g_string_new(token[0]);
				str = g_string_prepend(str, "<span foreground=\"red\" size=\"larger\"><b>");
				str = g_string_append(str, "</b>");
				if (token[1]) {
					for (int i = 1; token[i]; i++) {
						str = g_string_append(str, " ");
						if (!g_ascii_strcasecmp(current->name, token[i]))
							str = g_string_append(str, "\n<b>");
						str = g_string_append(str, token[i]);
						if (!g_ascii_strcasecmp(current->name, token[i]))
							str = g_string_append(str, "</b>");
					}
				}
				str = g_string_append(str, "</span>\n\n\t");
				str = g_string_append(str, _(current->definition));
				str = g_string_append(str, "\n\n<b>");
				str = g_string_append(str, _("Can be used in a script: "));
				if (current->scriptable) {
					str = g_string_append(str, "<span foreground=\"green\">");
					str = g_string_append(str, _("YES"));
				} else {
					str = g_string_append(str, "<span foreground=\"red\">");
					str = g_string_append(str, _("NO"));
				}
				str = g_string_append(str, "</span></b>");
				helper = g_string_free(str, FALSE);
				g_strfreev(token);
				break;
			}
			current++;
		}
		g_strfreev(command_line);
	}
	if (!helper)
		helper = g_strdup(_("No help for this command"));
	GtkWidget *popover = popover_new(GTK_WIDGET(gtk_builder_get_object(gui.builder, "command")), helper);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	g_free(helper);
	return FALSE;
}

static void init_completion_command(void) {
	GtkEntryCompletion *completion = gtk_entry_completion_new();
	GtkListStore *model = gtk_list_store_new(1, G_TYPE_STRING);
	GtkTreeIter iter;
	GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "command")));
	gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(model));
	gtk_entry_completion_set_text_column(completion, COMPLETION_COLUMN);
	gtk_entry_completion_set_minimum_key_length(completion, 2);
	gtk_entry_completion_set_popup_completion(completion, TRUE);
	gtk_entry_completion_set_inline_completion(completion, TRUE);
	gtk_entry_completion_set_popup_single_match(completion, FALSE);
	gtk_entry_completion_set_match_func(completion, completion_match_func, NULL, NULL);
	gtk_entry_set_completion(entry, completion);
	g_signal_connect(G_OBJECT(completion), "match-selected", G_CALLBACK(on_match_selected), NULL);
	command *current = commands;
	while (current->process) {
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter, COMPLETION_COLUMN, current->name, -1);
		current++;
	}
	g_object_unref(model);
}

static void init_controller_command(void) {
	GtkWidget *widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "command"));
#if GTK_CHECK_VERSION(3, 24, 24)
	GtkEventController *controller = gtk_event_controller_key_new(widget);
	g_signal_connect(controller, "key-pressed", G_CALLBACK(on_command_key_press_event), widget);
#else
	g_signal_connect(widget, "key-press-event", G_CALLBACK(on_command_key_press_event), NULL);
#endif
}

void init_command(void) {
	init_completion_command();
	init_controller_command();
}

void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data) {
	show_command_help_popup((GtkEntry *)user_data);
}
