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

/* Search entry for the stacked menu popovers of siril.ui. Items are read
 * from the popover's widgets at each search, so any entry added to the
 * .ui file is found without further registration. */

#include <string.h>

#include <gtk/gtk.h>

#include "gui-gtk4/menu_search.h"
#include "gui-gtk4/utils.h"

typedef struct {
	GtkPopover *popover;
	GtkStack *stack;
	gchar *nav_action;
	GtkWidget *entry;
	GtkWidget *listbox;
	GtkWidget *scroll;
	GtkWidget *separator;
	int width;
} MenuSearch;

static void menu_search_free(gpointer data) {
	MenuSearch *ms = data;
	g_free(ms->nav_action);
	g_free(ms);
}

/* casefold + strip accents, so that "debruit" matches "Débruitage" */
static gchar *fold(const gchar *s) {
	gchar *cf = g_utf8_casefold(s, -1);
	gchar *nfd = g_utf8_normalize(cf, -1, G_NORMALIZE_NFD);
	g_free(cf);
	if (!nfd)
		return g_strdup("");
	GString *out = g_string_sized_new(strlen(nfd));
	for (const gchar *p = nfd; *p; p = g_utf8_next_char(p)) {
		gunichar c = g_utf8_get_char(p);
		if (!g_unichar_ismark(c))
			g_string_append_unichar(out, c);
	}
	g_free(nfd);
	return g_string_free(out, FALSE);
}

static gboolean contains_folded(const gchar *text, const gchar *needle) {
	if (!text)
		return FALSE;
	gchar *hay = fold(text);
	gboolean found = strstr(hay, needle) != NULL;
	g_free(hay);
	return found;
}

static const gchar *find_label_text(GtkWidget *w) {
	if (GTK_IS_LABEL(w))
		return gtk_label_get_text(GTK_LABEL(w));
	for (GtkWidget *c = gtk_widget_get_first_child(w); c; c = gtk_widget_get_next_sibling(c)) {
		const gchar *t = find_label_text(c);
		if (t && *t)
			return t;
	}
	return NULL;
}

static void collect_buttons(GtkWidget *w, GPtrArray *out) {
	if (GTK_IS_BUTTON(w)) {
		g_ptr_array_add(out, w);
		return;
	}
	for (GtkWidget *c = gtk_widget_get_first_child(w); c; c = gtk_widget_get_next_sibling(c))
		collect_buttons(c, out);
}

static const gchar *nav_target(MenuSearch *ms, GtkWidget *btn) {
	if (g_strcmp0(gtk_actionable_get_action_name(GTK_ACTIONABLE(btn)), ms->nav_action))
		return NULL;
	GVariant *t = gtk_actionable_get_action_target_value(GTK_ACTIONABLE(btn));
	if (!t || !g_variant_is_of_type(t, G_VARIANT_TYPE_STRING))
		return NULL;
	return g_variant_get_string(t, NULL);
}

/* "Photometry › Output" for nested pages, NULL for root */
static gchar *page_path(const gchar *page, GHashTable *titles, GHashTable *parents) {
	gchar *path = NULL;
	for (int depth = 0; page && depth < 8; depth++) {
		const gchar *title = g_hash_table_lookup(titles, page);
		if (!title)
			break;
		gchar *tmp = path ? g_strdup_printf("%s \xe2\x80\xba %s", title, path) : g_strdup(title);
		g_free(path);
		path = tmp;
		page = g_hash_table_lookup(parents, page);
	}
	return path;
}

static void activate_row(MenuSearch *ms, GtkWidget *row) {
	const gchar *action = g_object_get_data(G_OBJECT(row), "menu-action");
	GVariant *target = g_object_get_data(G_OBJECT(row), "menu-target");
	if (action)
		gtk_widget_activate_action_variant(GTK_WIDGET(ms->popover), action, target);
	else /* buttons wired with a "clicked" handler instead of an action */
		g_signal_emit_by_name(g_object_get_data(G_OBJECT(row), "menu-button"), "clicked");
	gtk_popover_popdown(ms->popover);
}

static void on_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer user_data) {
	(void) box;
	activate_row((MenuSearch *) user_data, GTK_WIDGET(row));
}

static void on_entry_activate(GtkSearchEntry *entry, gpointer user_data) {
	(void) entry;
	MenuSearch *ms = user_data;
	for (GtkWidget *r = gtk_widget_get_first_child(ms->listbox); r; r = gtk_widget_get_next_sibling(r)) {
		if (gtk_widget_get_sensitive(r)) {
			activate_row(ms, r);
			return;
		}
	}
}

static GtkWidget *make_row(GtkWidget *btn, const gchar *label, const gchar *path) {
	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
	gtk_widget_set_margin_start(hbox, 6);
	gtk_widget_set_margin_end(hbox, 6);
	gtk_widget_set_margin_top(hbox, 3);
	gtk_widget_set_margin_bottom(hbox, 3);

	GtkWidget *name = gtk_label_new(label);
	gtk_label_set_xalign(GTK_LABEL(name), 0.0f);
	gtk_label_set_ellipsize(GTK_LABEL(name), PANGO_ELLIPSIZE_END);
	gtk_widget_set_hexpand(name, TRUE);
	gtk_box_append(GTK_BOX(hbox), name);
	if (path) {
		GtkWidget *where = gtk_label_new(path);
		gtk_widget_add_css_class(where, "dim-label");
		gtk_box_append(GTK_BOX(hbox), where);
	}

	GtkWidget *row = gtk_list_box_row_new();
	gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), hbox);
	gtk_widget_set_tooltip_text(row, label);
	gtk_widget_set_sensitive(row, gtk_widget_is_sensitive(btn));
	g_object_set_data(G_OBJECT(row), "menu-button", btn);
	g_object_set_data_full(G_OBJECT(row), "menu-action",
			g_strdup(gtk_actionable_get_action_name(GTK_ACTIONABLE(btn))), g_free);
	GVariant *target = gtk_actionable_get_action_target_value(GTK_ACTIONABLE(btn));
	if (target)
		g_object_set_data_full(G_OBJECT(row), "menu-target", g_variant_ref(target),
				(GDestroyNotify) g_variant_unref);
	return row;
}

static void on_search_changed(GtkSearchEntry *entry, gpointer user_data) {
	MenuSearch *ms = user_data;
	GtkWidget *child;
	while ((child = gtk_widget_get_first_child(ms->listbox)) != NULL)
		gtk_list_box_remove(GTK_LIST_BOX(ms->listbox), child);

	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
	if (!text || !*text) {
		gtk_widget_set_visible(ms->scroll, FALSE);
		gtk_widget_set_size_request(ms->scroll, -1, -1);
		ms->width = 0;
		gtk_widget_set_visible(ms->separator, TRUE);
		gtk_widget_set_visible(GTK_WIDGET(ms->stack), TRUE);
		return;
	}
	/* results replace the menu, at least as wide as it */
	if (gtk_widget_get_visible(GTK_WIDGET(ms->stack))) {
		ms->width = gtk_widget_get_width(GTK_WIDGET(ms->stack));
		gtk_widget_set_visible(GTK_WIDGET(ms->stack), FALSE);
		gtk_widget_set_visible(ms->separator, FALSE);
	}

	GHashTable *titles = g_hash_table_new(g_str_hash, g_str_equal);
	GHashTable *parents = g_hash_table_new(g_str_hash, g_str_equal);
	GPtrArray *pages = g_ptr_array_new_with_free_func((GDestroyNotify) g_ptr_array_unref);
	GPtrArray *names = g_ptr_array_new();
	GHashTable *by_name = g_hash_table_new(g_str_hash, g_str_equal);

	for (GtkWidget *p = gtk_widget_get_first_child(GTK_WIDGET(ms->stack)); p; p = gtk_widget_get_next_sibling(p)) {
		const gchar *page = gtk_stack_page_get_name(gtk_stack_get_page(ms->stack, p));
		GPtrArray *buttons = g_ptr_array_new();
		collect_buttons(p, buttons);
		g_ptr_array_add(pages, buttons);
		g_ptr_array_add(names, (gpointer) page);
		if (page)
			g_hash_table_insert(by_name, (gpointer) page, buttons);
	}

	/* breadth-first from root: the first link to a page names it, back
	 * buttons point to pages already seen */
	GPtrArray *queue = g_ptr_array_new();
	g_ptr_array_add(queue, "root");
	for (guint q = 0; q < queue->len; q++) {
		GPtrArray *buttons = g_hash_table_lookup(by_name, queue->pdata[q]);
		for (guint i = 0; buttons && i < buttons->len; i++) {
			const gchar *target = nav_target(ms, buttons->pdata[i]);
			if (target && strcmp(target, "root") && !g_hash_table_contains(parents, target)) {
				g_hash_table_insert(titles, (gpointer) target, (gpointer) find_label_text(buttons->pdata[i]));
				g_hash_table_insert(parents, (gpointer) target, queue->pdata[q]);
				g_ptr_array_add(queue, (gpointer) target);
			}
		}
	}
	g_ptr_array_unref(queue);
	g_hash_table_destroy(by_name);

	gchar *needle = fold(text);
	gint count = 0;
	for (guint i = 0; i < pages->len; i++) {
		GPtrArray *buttons = pages->pdata[i];
		gchar *path = page_path(names->pdata[i], titles, parents);
		gboolean path_match = contains_folded(path, needle);
		for (guint j = 0; j < buttons->len; j++) {
			GtkWidget *btn = buttons->pdata[j];
			if (nav_target(ms, btn) || !gtk_widget_get_visible(btn))
				continue;
			const gchar *label = find_label_text(btn);
			if (!label || !*label)
				continue;
			if (path_match || contains_folded(label, needle)) {
				gtk_list_box_append(GTK_LIST_BOX(ms->listbox), make_row(btn, label, path));
				count++;
			}
		}
		g_free(path);
	}
	g_free(needle);
	g_ptr_array_unref(names);
	g_ptr_array_unref(pages);
	g_hash_table_destroy(parents);
	g_hash_table_destroy(titles);

	/* grow to fit the labels, never shrink while typing so the popover doesn't jump */
	int natural;
	gtk_widget_measure(ms->listbox, GTK_ORIENTATION_HORIZONTAL, -1, NULL, &natural, NULL, NULL);
	ms->width = MAX(ms->width, natural);
	gtk_widget_set_size_request(ms->scroll, ms->width, -1);
	gtk_widget_set_visible(ms->scroll, count > 0);
}

static void on_popover_closed(GtkPopover *popover, gpointer user_data) {
	(void) popover;
	MenuSearch *ms = user_data;
	gtk_editable_set_text(GTK_EDITABLE(ms->entry), "");
}

/* Puts a search entry above the popover's "<popover_id>_stack". Typing while
 * the menu is open goes to the entry; results replace the menu until it is
 * cleared, so the popover never grows taller than the menu. */
void attach_menu_search(const char *popover_id, const char *action_namespace, const char *placeholder) {
	GtkPopover *popover = GTK_POPOVER(lookup_widget(popover_id));
	gchar *stack_id = g_strdup_printf("%s_stack", popover_id);
	GtkStack *stack = GTK_STACK(lookup_widget(stack_id));
	g_free(stack_id);
	if (!popover || !stack)
		return;

	MenuSearch *ms = g_new0(MenuSearch, 1);
	ms->popover = popover;
	ms->stack = stack;
	ms->nav_action = g_strdup_printf("%s.page", action_namespace);
	g_object_set_data_full(G_OBJECT(popover), "menu-search", ms, menu_search_free);

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);

	ms->entry = gtk_search_entry_new();
	gtk_search_entry_set_placeholder_text(GTK_SEARCH_ENTRY(ms->entry), placeholder);
	gtk_search_entry_set_key_capture_widget(GTK_SEARCH_ENTRY(ms->entry), GTK_WIDGET(popover));
	g_signal_connect(ms->entry, "search-changed", G_CALLBACK(on_search_changed), ms);
	g_signal_connect(ms->entry, "activate", G_CALLBACK(on_entry_activate), ms);
	gtk_box_append(GTK_BOX(box), ms->entry);

	ms->listbox = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(ms->listbox), GTK_SELECTION_NONE);
	g_signal_connect(ms->listbox, "row-activated", G_CALLBACK(on_row_activated), ms);

	ms->scroll = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(ms->scroll), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_max_content_height(GTK_SCROLLED_WINDOW(ms->scroll), 350);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(ms->scroll), TRUE);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(ms->scroll), ms->listbox);
	gtk_widget_set_visible(ms->scroll, FALSE);
	gtk_box_append(GTK_BOX(box), ms->scroll);

	ms->separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_box_append(GTK_BOX(box), ms->separator);

	g_object_ref(stack);
	gtk_popover_set_child(popover, NULL);
	gtk_box_append(GTK_BOX(box), GTK_WIDGET(stack));
	g_object_unref(stack);
	gtk_popover_set_child(popover, box);

	g_signal_connect(popover, "closed", G_CALLBACK(on_popover_closed), ms);
}
