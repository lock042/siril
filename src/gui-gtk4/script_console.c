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
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"
#include "script_console.h"

/* ── Status-bar logging ──────────────────────────────────────────────────── */

/* Traffic-light status indicator: the GtkImage is always
 * media-record-symbolic (a filled circle), tinted red while idle and
 * green while a script is running.  Colours live in siril.css under
 * #image_log / #image_log.script-running. */
static void update_log_icon(gboolean is_running) {
	GtkWidget *image = GTK_WIDGET(gtk_builder_get_object(gui.builder, "image_log"));
	if (is_running)
		gtk_widget_add_css_class(image, "script-running");
	else
		gtk_widget_remove_css_class(image, "script-running");
}

struct log_status_bar_idle_data {
	gchar *myline;
	int line;
};

static gboolean log_status_bar_idle_callback(gpointer p) {
	struct log_status_bar_idle_data *data = (struct log_status_bar_idle_data *) p;
	GtkLabel *statusbar = GTK_LABEL(gtk_builder_get_object(gui.builder, "statusbar_script"));
	gchar *newline = g_strdup(data->myline);
	/* line <= 0 is an interactive command: show it alone, without the
	 * "Processing line N:" prefix that only makes sense inside a script. */
	gchar *status = (data->line > 0)
			? g_strdup_printf(_("Processing line %d: %s"), data->line, newline)
			: g_strdup_printf(_("Running command: %s"), newline);
	update_log_icon(TRUE);
	gtk_label_set_text(statusbar, status);
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
	g_idle_add(log_status_bar_idle_callback, data);
}

void console_clear_status_bar(void) {
	GtkLabel *bar = GTK_LABEL(gtk_builder_get_object(gui.builder, "statusbar_script"));
	gtk_label_set_text(bar, "");
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

/* The completion popover and its associated model live here; their full
 * setup is below in init_completion_command().  Forward-declared so the
 * key handler can Tab into the popover. */
typedef struct {
	GtkEntry          *entry;
	GtkPopover        *popover;
	GtkStringList     *all_commands;
	GtkStringFilter   *filter;
	GtkSingleSelection *selection;
	GtkListView       *list_view;
} CommandCompletion;

static CommandCompletion g_cmd_complete = { 0 };

void on_command_activate(GtkEntry *entry, gpointer user_data) {
	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
	history_add_line(strdup(text));
	if (!(processcommand(text, FALSE))) {
		gtk_editable_set_text(GTK_EDITABLE(entry), "");
		gui_iface.on_precision_changed();
	}
}

static gboolean on_command_key_press_event(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType modifiers,
		GtkWidget *widget) {
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (keyval) {
	case GDK_KEY_Tab:
	case GDK_KEY_ISO_Left_Tab:
		/* Tab into the completion popover when it's open.  GTK's default
		 * forward-focus is racy here: gtk_popover_popup() defers actual
		 * realisation/focus-chain rebuild by a frame or two, so a Tab
		 * pressed quickly after the popover appears would miss it and
		 * move focus to whatever the next widget in the window is.
		 * Doing it explicitly here is deterministic.  Seed a selection
		 * so arrow keys immediately navigate. */
		if (gtk_widget_get_visible(GTK_WIDGET(g_cmd_complete.popover))) {
			GtkSingleSelection *sel = g_cmd_complete.selection;
			if (gtk_single_selection_get_selected(sel) == GTK_INVALID_LIST_POSITION
			    && g_list_model_get_n_items(G_LIST_MODEL(sel)) > 0) {
				gtk_single_selection_set_selected(sel, 0);
			}
			gtk_widget_grab_focus(GTK_WIDGET(g_cmd_complete.list_view));
			handled = 1;
		}
		break;
	case GDK_KEY_Up:
		handled = 1;
		if (!gui.cmd_history)
			break;
		if (gui.cmd_hist_display > 0) {
			if (gui.cmd_history[gui.cmd_hist_display - 1])
				--gui.cmd_hist_display;
			gtk_editable_set_text(GTK_EDITABLE(entry), gui.cmd_history[gui.cmd_hist_display]);
		} else if (gui.cmd_history[gui.cmd_hist_size - 1]) {
			gui.cmd_hist_display = gui.cmd_hist_size - 1;
			gtk_editable_set_text(GTK_EDITABLE(entry), gui.cmd_history[gui.cmd_hist_display]);
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
			gtk_editable_set_text(GTK_EDITABLE(entry), gui.cmd_history[gui.cmd_hist_display]);
		} else {
			gtk_editable_set_text(GTK_EDITABLE(entry), "");
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	default:
		break;
	}
	return (handled != 0);
}

/* ── Command-line completion (custom GtkPopover + GtkListView) ────────────
 * GtkEntryCompletion is deprecated in GTK4 (and removed in GTK5) and is
 * tightly coupled to GtkTreeModel.  Phase 11 replaces both with a self-
 * contained popup: a GtkStringList of all command names, wrapped in a
 * GtkStringFilter (prefix match) bound to the GtkEntry's text, displayed
 * in a GtkListView inside a GtkPopover anchored at the entry.            */

/* When focus returns to the entry after popover dismissal, GtkText's
 * select-on-focus-in highlights all the just-inserted text — so the user's
 * next keystroke wipes the completion they just chose.  Defer a
 * select_region(-1,-1) to the next main-loop turn, after focus has settled,
 * to collapse the selection to a zero-width cursor at the end.  This is
 * essential for the Tab-to-completion path. */
static gboolean clear_entry_selection_idle(gpointer data) {
	GtkEditable *e = GTK_EDITABLE(data);
	gtk_editable_select_region(e, -1, -1);
	return G_SOURCE_REMOVE;
}

static void cmd_complete_apply(CommandCompletion *cc, guint pos) {
	if (!cc || pos == GTK_INVALID_LIST_POSITION) return;
	GListModel *m = G_LIST_MODEL(cc->selection);
	GObject *o = g_list_model_get_item(m, pos);
	if (!o) return;
	const gchar *cmd = NULL;
	if (GTK_IS_STRING_OBJECT(o))
		cmd = gtk_string_object_get_string(GTK_STRING_OBJECT(o));
	if (cmd) {
		GtkEditable *e = GTK_EDITABLE(cc->entry);
		/* Block on_cmd_entry_changed while we replace the text — the
		 * delete + insert each fire "changed" synchronously, which would
		 * otherwise re-popup the completion list for the just-selected
		 * prefix (e.g. selecting "load" would refilter to load, loadfits,
		 * …).  The explicit popdown below can lose to GTK's queued show. */
		g_signal_handlers_block_matched(cc->entry, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, cc);
		gtk_editable_delete_text(e, 0, -1);
		gtk_editable_insert_text(e, cmd, -1, &(int){0});
		gtk_editable_set_position(e, -1);
		g_signal_handlers_unblock_matched(cc->entry, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, cc);
		g_idle_add(clear_entry_selection_idle, cc->entry);
	}
	g_object_unref(o);
	gtk_widget_set_visible(GTK_WIDGET(cc->popover), FALSE);
}

static void on_cmd_complete_activate(GtkListView *self, guint position, gpointer user_data) {
	(void)self;
	cmd_complete_apply((CommandCompletion *)user_data, position);
}

static void on_cmd_entry_changed(GtkEditable *editable, gpointer user_data) {
	CommandCompletion *cc = (CommandCompletion *)user_data;
	const gchar *text = gtk_editable_get_text(editable);
	if (!text || strlen(text) < 2) {
		gtk_widget_set_visible(GTK_WIDGET(cc->popover), FALSE);
		return;
	}
	gtk_string_filter_set_search(cc->filter, text);
	guint n = g_list_model_get_n_items(G_LIST_MODEL(cc->selection));
	if (n == 0 || (n == 1 && g_strcmp0(text, gtk_string_object_get_string(
			GTK_STRING_OBJECT(g_list_model_get_item(G_LIST_MODEL(cc->selection), 0)))) == 0)) {
		gtk_widget_set_visible(GTK_WIDGET(cc->popover), FALSE);
		return;
	}
	if (!gtk_widget_get_visible(GTK_WIDGET(cc->popover)))
		gtk_popover_popup(cc->popover);
}

static void on_cmd_entry_activate(GtkEntry *entry, gpointer user_data) {
	(void)entry;
	CommandCompletion *cc = (CommandCompletion *)user_data;
	if (!gtk_widget_get_visible(GTK_WIDGET(cc->popover))) return;
	guint pos = gtk_single_selection_get_selected(cc->selection);
	cmd_complete_apply(cc, pos);
}

static void cmd_completion_setup_factory(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	gtk_list_item_set_child(li, gtk_label_new(NULL));
}
static void cmd_completion_bind_factory(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	GObject *o = gtk_list_item_get_item(li);
	const gchar *s = (o && GTK_IS_STRING_OBJECT(o))
			? gtk_string_object_get_string(GTK_STRING_OBJECT(o)) : "";
	gtk_label_set_xalign(lbl, 0.0);
	gtk_label_set_text(lbl, s);
}

gboolean show_command_help_popup(gpointer user_data) {
	GtkEntry *entry = (GtkEntry *) user_data;
	gchar *helper = NULL;
	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
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
	gtk_widget_set_visible(popover, TRUE);
#endif
	g_free(helper);
	return FALSE;
}

static void init_completion_command(void) {
	GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "command")));
	g_cmd_complete.entry = entry;

	/* Build the command-name model */
	GtkStringList *all = gtk_string_list_new(NULL);
	for (command *current = commands; current->process; current++) {
		gtk_string_list_append(all, current->name);
	}
	g_cmd_complete.all_commands = all;

	/* Prefix-match filter using the GtkStringObject's "string" property */
	GtkExpression *expr = gtk_property_expression_new(GTK_TYPE_STRING_OBJECT, NULL, "string");
	GtkStringFilter *filter = gtk_string_filter_new(expr);
	gtk_string_filter_set_match_mode(filter, GTK_STRING_FILTER_MATCH_MODE_PREFIX);
	gtk_string_filter_set_ignore_case(filter, TRUE);
	g_cmd_complete.filter = filter;

	GtkFilterListModel *fm = gtk_filter_list_model_new(G_LIST_MODEL(all), GTK_FILTER(filter));
	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(fm));
	gtk_single_selection_set_autoselect(sel, FALSE);
	gtk_single_selection_set_can_unselect(sel, TRUE);
	g_cmd_complete.selection = sel;

	GtkSignalListItemFactory *fac = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fac, "setup", G_CALLBACK(cmd_completion_setup_factory), NULL);
	g_signal_connect(fac, "bind",  G_CALLBACK(cmd_completion_bind_factory),  NULL);

	GtkListView *lv = GTK_LIST_VIEW(gtk_list_view_new(GTK_SELECTION_MODEL(g_object_ref(sel)),
			GTK_LIST_ITEM_FACTORY(fac)));
	gtk_list_view_set_single_click_activate(lv, TRUE);
	g_signal_connect(lv, "activate", G_CALLBACK(on_cmd_complete_activate), &g_cmd_complete);
	g_cmd_complete.list_view = lv;

	GtkWidget *sw = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(sw), TRUE);
	gtk_widget_set_size_request(sw, 240, -1);
	gtk_scrolled_window_set_max_content_height(GTK_SCROLLED_WINDOW(sw), 200);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sw), GTK_WIDGET(lv));

	GtkPopover *pop = GTK_POPOVER(gtk_popover_new());
	gtk_popover_set_autohide(pop, FALSE);
	gtk_popover_set_has_arrow(pop, FALSE);
	gtk_popover_set_position(pop, GTK_POS_BOTTOM);
	gtk_widget_set_parent(GTK_WIDGET(pop), GTK_WIDGET(entry));
	gtk_popover_set_child(pop, sw);
	g_cmd_complete.popover = pop;

	g_signal_connect(entry, "changed",  G_CALLBACK(on_cmd_entry_changed),  &g_cmd_complete);
	g_signal_connect(entry, "activate", G_CALLBACK(on_cmd_entry_activate), &g_cmd_complete);
}

static void init_controller_command(void) {
	GtkWidget *widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "command"));
	GtkEventController *controller = gtk_event_controller_key_new();
	g_signal_connect(controller, "key-pressed", G_CALLBACK(on_command_key_press_event), widget);
	gtk_widget_add_controller(widget, controller);
}

void init_command(void) {
	init_completion_command();
	init_controller_command();
}

void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data) {
	show_command_help_popup((GtkEntry *)user_data);
}
