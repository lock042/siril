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
#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/photometric_cc.h" // for reset_spcc_filters() (this is not a GTK function)
#include "gui-gtk4/preferences.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/script_menu.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/python_gui.h"
#include "io/siril_git.h"

#ifdef HAVE_LIBGIT2

/* ---- Per-row GObject for the script repo list ------------------- */
#define SIRIL_TYPE_SCRIPT_ROW (siril_script_row_get_type())
G_DECLARE_FINAL_TYPE(SirilScriptRow, siril_script_row, SIRIL, SCRIPT_ROW, GObject)
struct _SirilScriptRow {
	GObject parent_instance;
	gchar   *category;
	gchar   *scriptname;
	gboolean selected;
	gchar   *scriptpath;
	gchar   *type;
	gboolean startup;
	gboolean is_python;        /* startup-capable */
};
G_DEFINE_TYPE(SirilScriptRow, siril_script_row, G_TYPE_OBJECT)
static void siril_script_row_finalize(GObject *obj) {
	SirilScriptRow *self = SIRIL_SCRIPT_ROW(obj);
	g_clear_pointer(&self->category,   g_free);
	g_clear_pointer(&self->scriptname, g_free);
	g_clear_pointer(&self->scriptpath, g_free);
	g_clear_pointer(&self->type,       g_free);
	G_OBJECT_CLASS(siril_script_row_parent_class)->finalize(obj);
}
static void siril_script_row_class_init(SirilScriptRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_script_row_finalize;
}
static void siril_script_row_init(SirilScriptRow *self) { (void)self; }

static GListStore       *list_store         = NULL;
/* Forward decl: defined later, referenced from ensure_git_view(). */
static void on_git_columnview_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static GtkCustomFilter  *script_filter      = NULL;
static GtkColumnView    *git_columnview     = NULL;
static GtkScrolledWindow *git_scrolled      = NULL;
static gchar            *current_search_text = NULL;
static gboolean          filter_enabled      = FALSE;
static GtkWindow        *git_control_window  = NULL;
static GtkWidget        *git_pref_auto_updates = NULL;
static GtkWidget        *git_manual_sync_btn = NULL;
static GtkWidget        *git_spcc_sync_startup = NULL;
static GtkWidget        *git_spcc_manual_sync = NULL;

static void git_gui_init_statics(void) {
	if (git_scrolled) return;
	git_scrolled        = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolled_script_repo"));
	git_control_window  = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	git_pref_auto_updates = GTK_WIDGET(gtk_builder_get_object(gui.builder, "pref_script_automatic_updates"));
	git_manual_sync_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "manual_script_sync_button"));
	git_spcc_sync_startup = GTK_WIDGET(gtk_builder_get_object(gui.builder, "spcc_repo_sync_at_startup"));
	git_spcc_manual_sync = GTK_WIDGET(gtk_builder_get_object(gui.builder, "spcc_repo_manual_sync"));
}

static gboolean script_is_startup_capable(const gchar *path) {
	FILE *f = g_fopen(path, "r");
	if (!f)
		return FALSE;

	gchar line[512];
	gboolean result = FALSE;

	while (fgets(line, sizeof(line), f)) {
		gsize len = strlen(line);
		if (len > 0 && line[len - 1] == '\n')
			line[--len] = '\0';

		if (line[0] != '#') {
			break;
		}
		if (strstr(line, "STARTUP_CAPABLE")) {
			result = TRUE;
			break;
		}
	}

	fclose(f);
	return result;
}

/* ---- GtkCustomFilter callback ----------------------------------- */
static gboolean script_filter_match(gpointer item, gpointer user_data) {
	(void)user_data;
	if (!filter_enabled || !current_search_text || !*current_search_text)
		return TRUE;
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(item);
	if (!row->scriptname) return FALSE;
	gchar *key  = g_ascii_strdown(current_search_text, -1);
	gchar *name = g_ascii_strdown(row->scriptname, -1);
	gboolean v  = (strstr(name, key) != NULL);
	g_free(key); g_free(name);
	return v;
}

void on_find_script_entry_changed(GtkEntry *entry, gpointer user_data) {
	(void)user_data;
	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));

	g_free(current_search_text);
	current_search_text = g_strdup(text);
	filter_enabled = (text && *text != '\0');

	if (script_filter)
		gtk_filter_changed(GTK_FILTER(script_filter), GTK_FILTER_CHANGE_DIFFERENT);
}

/* ---- Cell factories --------------------------------------------- */

typedef enum {
	SCRIPT_COL_CATEGORY, SCRIPT_COL_SCRIPTNAME, SCRIPT_COL_TYPE
} ScriptStringCol;

static void script_string_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void script_string_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	ScriptStringCol kind = (ScriptStringCol)GPOINTER_TO_INT(u);
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(gtk_list_item_get_item(li));
	const gchar *s = "";
	switch (kind) {
		case SCRIPT_COL_CATEGORY:   s = row->category   ? row->category   : ""; break;
		case SCRIPT_COL_SCRIPTNAME: s = row->scriptname ? row->scriptname : ""; break;
		case SCRIPT_COL_TYPE:       s = row->type       ? row->type       : ""; break;
	}
	gtk_label_set_text(lbl, s);
	/* Attach the row to the cell widget so the secondary-click handler
	 * on the column view can recover the SirilScriptRow from a
	 * gtk_widget_pick walk — the GtkColumnView's selection is unreliable
	 * for right-click (the user may right-click an unselected row). */
	g_object_set_data(G_OBJECT(lbl), "siril-script-row", row);
}
static GtkListItemFactory *make_string_factory(ScriptStringCol c) {
	GtkListItemFactory *f = gtk_signal_list_item_factory_new();
	g_signal_connect(f, "setup", G_CALLBACK(script_string_setup_cb), NULL);
	g_signal_connect(f, "bind",  G_CALLBACK(script_string_bind_cb),  GINT_TO_POINTER(c));
	return f;
}

/* "Sel" toggle column: simple inline GtkCheckButton wired to row->selected. */
static void on_selected_toggled(GtkCheckButton *btn, gpointer user_data) {
	(void)user_data;
	GtkListItem *li = g_object_get_data(G_OBJECT(btn), "siril-listitem");
	if (!li) return;
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gboolean new_val = siril_toggle_get_active(GTK_WIDGET(btn));
	row->selected = new_val;

	if (new_val) {
		if (!g_slist_find_custom(com.pref.selected_scripts, row->scriptpath, (GCompareFunc)g_strcmp0)) {
#ifdef DEBUG_GITSCRIPTS
			printf("Adding script: %s\n", row->scriptpath);
#endif
			com.pref.selected_scripts = g_slist_prepend(com.pref.selected_scripts, g_strdup(row->scriptpath));
		}
	} else {
		GSList *found = g_slist_find_custom(com.pref.selected_scripts, row->scriptpath, (GCompareFunc)g_strcmp0);
		if (found) {
#ifdef DEBUG_GITSCRIPTS
			printf("Removing script: %s\n", row->scriptpath);
#endif
			g_free(found->data);
			com.pref.selected_scripts = g_slist_delete_link(com.pref.selected_scripts, found);
		}
	}
	notify_script_update();
}

static void selected_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *btn = gtk_check_button_new();
	gtk_widget_set_halign(btn, GTK_ALIGN_CENTER);
	gtk_list_item_set_child(li, btn);
	g_signal_connect(btn, "toggled", G_CALLBACK(on_selected_toggled), NULL);
}
static void selected_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", li);
	g_object_set_data(G_OBJECT(btn), "siril-script-row", row);
	g_signal_handlers_block_by_func(btn, on_selected_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(btn), row->selected);
	g_signal_handlers_unblock_by_func(btn, on_selected_toggled, NULL);
}
static void selected_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", NULL);
}

/* "Startup" toggle: only shown when row->is_python is TRUE. */
static void on_startup_toggled(GtkCheckButton *btn, gpointer user_data) {
	(void)user_data;
	GtkListItem *li = g_object_get_data(G_OBJECT(btn), "siril-listitem");
	if (!li) return;
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	if (!row->is_python) return;

	gboolean new_val = siril_toggle_get_active(GTK_WIDGET(btn));
	row->startup = new_val;

	gchar *canonical_path = g_canonicalize_filename(row->scriptpath, NULL);

	if (new_val) {
		if (!g_slist_find_custom(com.pref.startup_scripts, canonical_path, (GCompareFunc)g_strcmp0)) {
#ifdef DEBUG_GITSCRIPTS
			printf("Adding startup script: %s\n", canonical_path);
#endif
			com.pref.startup_scripts = g_slist_prepend(com.pref.startup_scripts, canonical_path);
		} else {
			g_free(canonical_path);
		}
	} else {
		GSList *found = g_slist_find_custom(com.pref.startup_scripts, canonical_path, (GCompareFunc)g_strcmp0);
		if (found) {
#ifdef DEBUG_GITSCRIPTS
			printf("Removing startup script: %s\n", canonical_path);
#endif
			g_free(found->data);
			com.pref.startup_scripts = g_slist_delete_link(com.pref.startup_scripts, found);
		}
		g_free(canonical_path);
	}
}

static void startup_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *btn = gtk_check_button_new();
	gtk_widget_set_halign(btn, GTK_ALIGN_CENTER);
	gtk_list_item_set_child(li, btn);
	g_signal_connect(btn, "toggled", G_CALLBACK(on_startup_toggled), NULL);
}
static void startup_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", li);
	g_object_set_data(G_OBJECT(btn), "siril-script-row", row);
	gtk_widget_set_visible(GTK_WIDGET(btn), row->is_python);
	gtk_widget_set_sensitive(GTK_WIDGET(btn), row->is_python);
	g_signal_handlers_block_by_func(btn, on_startup_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(btn), row->startup);
	g_signal_handlers_unblock_by_func(btn, on_startup_toggled, NULL);
}
static void startup_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", NULL);
}

/* ---- Sort comparators on the row GObject ------------------------ */

static gint cmp_strings(const gchar *a, const gchar *b) {
	return g_strcmp0(a ? a : "", b ? b : "");
}
static gint cmp_category(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return cmp_strings(((const SirilScriptRow*)a)->category, ((const SirilScriptRow*)b)->category);
}
static gint cmp_scriptname(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return cmp_strings(((const SirilScriptRow*)a)->scriptname, ((const SirilScriptRow*)b)->scriptname);
}
static gint cmp_type(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return cmp_strings(((const SirilScriptRow*)a)->type, ((const SirilScriptRow*)b)->type);
}
static gint cmp_selected(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return ((const SirilScriptRow*)a)->selected - ((const SirilScriptRow*)b)->selected;
}
static gint cmp_startup(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return ((const SirilScriptRow*)a)->startup - ((const SirilScriptRow*)b)->startup;
}

/* ---- Row-activate handler (double-click → open in editor) ------- */

static void on_row_activated(GtkColumnView *cv, guint position, gpointer user_data) {
	(void)cv; (void)user_data;
	GtkSelectionModel *sm = gtk_column_view_get_model(cv);
	if (!sm) return;
	SirilScriptRow *row = SIRIL_SCRIPT_ROW(g_list_model_get_item(G_LIST_MODEL(sm), position));
	if (!row) return;

	gchar *tmpscriptpath = g_canonicalize_filename(row->scriptpath, NULL);
	gsize length = 0;
	gchar *contents = get_script_content_string_from_file_revision(tmpscriptpath, 0, &length, NULL, NULL);
	if (length > 0 && contents != NULL) {
		const char *ext = get_filename_ext(tmpscriptpath);
		new_script(contents, length, ext);
		g_free(contents);
	} else {
		gchar *msg = g_strdup_printf(_("Error loading script contents: %s\n"), tmpscriptpath);
		siril_log_error(msg);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
		g_free(msg);
	}
	g_free(tmpscriptpath);
	g_object_unref(row);
}

/* ---- Build the column view (idempotent) ------------------------- */

static void ensure_git_view(void) {
	if (git_columnview) return;
	git_gui_init_statics();
	if (!git_scrolled) return;

	if (!list_store)
		list_store = g_list_store_new(SIRIL_TYPE_SCRIPT_ROW);
	script_filter = gtk_custom_filter_new(script_filter_match, NULL, NULL);
	GtkFilterListModel *fm = gtk_filter_list_model_new(G_LIST_MODEL(g_object_ref(list_store)),
			GTK_FILTER(g_object_ref(script_filter)));

	/* SortListModel wraps the filtered model; per-column sorters set later. */
	GtkSortListModel *sm = gtk_sort_list_model_new(G_LIST_MODEL(fm), NULL);
	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(sm));
	gtk_single_selection_set_can_unselect(sel, TRUE);
	gtk_single_selection_set_autoselect(sel, FALSE);

	git_columnview = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));
	gtk_column_view_set_show_column_separators(git_columnview, TRUE);
	g_signal_connect(git_columnview, "activate", G_CALLBACK(on_row_activated), NULL);

	GtkColumnViewColumn *c;
	c = gtk_column_view_column_new(N_("Category"), make_string_factory(SCRIPT_COL_CATEGORY));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_expand(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_category, NULL, NULL)));
	gtk_column_view_append_column(git_columnview, c); g_object_unref(c);

	c = gtk_column_view_column_new(N_("Script"), make_string_factory(SCRIPT_COL_SCRIPTNAME));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_expand(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_scriptname, NULL, NULL)));
	gtk_column_view_append_column(git_columnview, c); g_object_unref(c);

	c = gtk_column_view_column_new(N_("Type"), make_string_factory(SCRIPT_COL_TYPE));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_type, NULL, NULL)));
	gtk_column_view_append_column(git_columnview, c); g_object_unref(c);

	GtkSignalListItemFactory *fsel = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fsel, "setup",  G_CALLBACK(selected_setup_cb),  NULL);
	g_signal_connect(fsel, "bind",   G_CALLBACK(selected_bind_cb),   NULL);
	g_signal_connect(fsel, "unbind", G_CALLBACK(selected_unbind_cb), NULL);
	c = gtk_column_view_column_new(N_("Sel"), GTK_LIST_ITEM_FACTORY(fsel));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_selected, NULL, NULL)));
	gtk_column_view_append_column(git_columnview, c); g_object_unref(c);

	GtkSignalListItemFactory *fst = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fst, "setup",  G_CALLBACK(startup_setup_cb),  NULL);
	g_signal_connect(fst, "bind",   G_CALLBACK(startup_bind_cb),   NULL);
	g_signal_connect(fst, "unbind", G_CALLBACK(startup_unbind_cb), NULL);
	c = gtk_column_view_column_new(N_("Startup"), GTK_LIST_ITEM_FACTORY(fst));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_startup, NULL, NULL)));
	gtk_column_view_append_column(git_columnview, c); g_object_unref(c);

	gtk_scrolled_window_set_child(git_scrolled, GTK_WIDGET(git_columnview));

	/* GTK4: secondary-click on a script row opens the revision dialog.
	 * Replaces the GTK3 "button-press-event" on the legacy GtkTreeView.
	 * CAPTURE phase ensures we run before the GtkColumnView row's own
	 * gestures consume the secondary click. */
	GtkGesture *click = gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), GDK_BUTTON_SECONDARY);
	gtk_event_controller_set_propagation_phase(GTK_EVENT_CONTROLLER(click),
	                                           GTK_PHASE_CAPTURE);
	g_signal_connect(click, "pressed", G_CALLBACK(on_git_columnview_pressed), NULL);
	gtk_widget_add_controller(GTK_WIDGET(git_columnview), GTK_EVENT_CONTROLLER(click));
}

/* ---- Populate the model ---------------------------------------- */

static gboolean fill_script_repo_tree_idle(gpointer p) {
	(void)p;
	ensure_git_view();
	if (!list_store) return G_SOURCE_REMOVE;

	g_list_store_remove_all(list_store);

	gui_repo_scripts_mutex_lock();
	if (com.repo_scripts) {
		GSList *iterator;
		for (iterator = com.repo_scripts; iterator; iterator = iterator->next) {
			const gchar *category;
			gboolean included = FALSE;
			gboolean startup = FALSE;
			gboolean core = FALSE;
			gchar *cat_owned = NULL;
			if (test_last_subdir((gchar *)iterator->data, "preprocessing")) {
				category = _("Preprocessing");
			} else if (test_last_subdir((gchar *)iterator->data, "processing")) {
				category = _("Processing");
			} else if (test_last_subdir((gchar *)iterator->data, "utility")) {
				category = _("Utility");
			} else if (test_last_subdir((gchar *)iterator->data, "core")) {
				category = _("Core");
				core = TRUE;
			} else {
				gchar *path = (gchar *)iterator->data;
				gchar *last_slash = strrchr(path, G_DIR_SEPARATOR);
				if (last_slash && last_slash > path) {
					gchar *prev_slash = g_strrstr_len(path, last_slash - path, G_DIR_SEPARATOR_S);
					gchar *subdir_start = prev_slash ? prev_slash + 1 : path;
					gchar *subdir = g_strndup(subdir_start, last_slash - subdir_start);
					if (subdir && subdir[0]) {
						subdir[0] = g_ascii_toupper(subdir[0]);
						for (gsize i = 1; subdir[i]; i++) {
							subdir[i] = g_ascii_tolower(subdir[i]);
						}
					}
					cat_owned = subdir;
					category = subdir;
				} else {
					category = _("Other");
				}
			}
			gchar *scriptname = g_path_get_basename((gchar *)iterator->data);
			gchar *scriptpath = g_build_path(G_DIR_SEPARATOR_S, siril_get_scripts_repo_path(), (gchar *)iterator->data, NULL);
			const gchar *scripttype;
			gboolean startup_capable = FALSE;
			if (g_str_has_suffix(scriptname, SCRIPT_EXT))
				scripttype = _("Siril Script File");
			else if (g_str_has_suffix(scriptname, PYSCRIPT_EXT) || g_str_has_suffix(scriptname, PYCSCRIPT_EXT)) {
				scripttype = _("Python script");
				startup_capable = script_is_startup_capable(scriptpath);
			} else scripttype = NULL;

			GSList *iterator2 = NULL;
			if (!included && !core) {
				for (iterator2 = com.pref.selected_scripts; iterator2; iterator2 = iterator2->next) {
					if (g_strrstr((gchar *)iterator2->data, (gchar *)iterator->data)) {
						included = TRUE;
					}
				}
			}

			if (!startup_capable && !core) {
				startup = FALSE;
			} else if (!startup && !core) {
				gchar *canonical_scriptpath = g_canonicalize_filename(scriptpath, NULL);
				for (iterator2 = com.pref.startup_scripts; iterator2; iterator2 = iterator2->next) {
					if (g_strcmp0((gchar *)iterator2->data, canonical_scriptpath) == 0) {
						startup = TRUE;
						break;
					}
				}
				g_free(canonical_scriptpath);
			}

			if (!core) {
				SirilScriptRow *row = g_object_new(SIRIL_TYPE_SCRIPT_ROW, NULL);
				row->category   = g_strdup(category);
				row->scriptname = scriptname; /* takes ownership */
				row->scriptpath = scriptpath; /* takes ownership */
				row->type       = g_strdup(scripttype ? scripttype : "");
				row->selected   = included;
				row->startup    = startup;
				row->is_python  = startup_capable;
				g_list_store_append(list_store, row);
				g_object_unref(row);
			} else {
				g_free(scriptname);
				g_free(scriptpath);
			}
			g_free(cat_owned);
		}
	}
	gui_repo_scripts_mutex_unlock();

	return G_SOURCE_REMOVE;
}


/* called on preference window loading.
* It is executed safely in the GTK thread if as_idle is true. */
void fill_script_repo_tree(gboolean as_idle) {
	git_gui_init_statics();
	if (as_idle)
		g_idle_add(fill_script_repo_tree_idle, NULL);
	else
		fill_script_repo_tree_idle(NULL);
}

typedef struct {
	const gchar *scriptpath;
	GtkTextView *message_textview;
} CommitMessageUpdateData;

/* Live commit-message updater for the script-revision spin button.
 * `user_data` is a CommitMessageUpdateData* describing which script the
 * spin is browsing. */
static void on_script_revision_spin_value_changed(GtkSpinButton *spin, gpointer user_data) {
	CommitMessageUpdateData *data = (CommitMessageUpdateData *)user_data;
	int revisions_back = gtk_spin_button_get_value_as_int(spin);
	GtkTextBuffer *buffer = gtk_text_view_get_buffer(data->message_textview);
	gchar *commit_message = NULL;
	size_t message_size = 0;
	gchar *content = get_script_content_string_from_file_revision(
		data->scriptpath, revisions_back, &(size_t){0}, &commit_message, &message_size);
	if (content) g_free(content); /* we only want the message here */
	if (commit_message && message_size > 0)
		gtk_text_buffer_set_text(buffer, commit_message, -1);
	else
		gtk_text_buffer_set_text(buffer, _("(No commit message found for this revision)"), -1);
	g_free(commit_message);
}

/* Result-bridge for the synchronous "Select Revision" dialog.  GtkDialog
 * with run() is gone; we use a GtkWindow + nested GMainLoop, with the
 * Cancel/OK buttons writing into this struct before quitting the loop. */
typedef struct {
	GMainLoop *loop;
	gint       result;
} rev_dlg_ctx;

static void rev_dlg_cancel(GtkButton *b, gpointer p) {
	(void)b;
	rev_dlg_ctx *ctx = (rev_dlg_ctx *)p;
	ctx->result = GTK_RESPONSE_CANCEL;
	if (g_main_loop_is_running(ctx->loop)) g_main_loop_quit(ctx->loop);
}
static void rev_dlg_ok(GtkButton *b, gpointer p) {
	(void)b;
	rev_dlg_ctx *ctx = (rev_dlg_ctx *)p;
	ctx->result = GTK_RESPONSE_OK;
	if (g_main_loop_is_running(ctx->loop)) g_main_loop_quit(ctx->loop);
}
static gboolean rev_dlg_close_request(GtkWindow *w, gpointer p) {
	(void)w;
	rev_dlg_ctx *ctx = (rev_dlg_ctx *)p;
	ctx->result = GTK_RESPONSE_CANCEL;
	if (g_main_loop_is_running(ctx->loop)) g_main_loop_quit(ctx->loop);
	return TRUE;
}

/* Pop the "Select Revision" dialog for the given script.  Modal, synchronous
 * (nested GMainLoop), and on OK loads the chosen revision into the editor. */
static void show_script_revision_dialog(const gchar *scriptpath) {
	if (!scriptpath) return;
	GtkWidget *dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Select Revision"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(dialog), TRUE);
	if (git_control_window)
		gtk_window_set_transient_for(GTK_WINDOW(dialog), git_control_window);
	gtk_window_set_resizable(GTK_WINDOW(dialog), TRUE);

	const gchar *tooltip = _("Leave at 0 to open the current version of the script. If you have a problem with "
		"a script update you can use this to go back to earlier revisions. Start by going back 1 revison "
		"and increase the number until you find the last version that worked for you. Note that this will "
		"only open the previous version in the script editor, it will not revert the script in the local "
		"repository, but it allows you to save a local copy in one of your script folders for use until the "
		"upstream script is fixed.");

	GtkWidget *content = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
	gtk_widget_set_margin_start(content, 10);
	gtk_widget_set_margin_end(content, 10);
	gtk_widget_set_margin_top(content, 10);
	gtk_widget_set_margin_bottom(content, 10);
	gtk_window_set_child(GTK_WINDOW(dialog), content);

	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
	GtkWidget *label = gtk_label_new(_("Revisions to go back (leave at 0 to get current revision):"));
	gtk_widget_set_tooltip_text(label, tooltip);
	GtkWidget *spin = gtk_spin_button_new_with_range(0, 999, 1);
	gtk_widget_set_tooltip_text(spin, tooltip);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), 0);
	gtk_box_append(GTK_BOX(hbox), label);
	gtk_box_append(GTK_BOX(hbox), spin);
	gtk_box_append(GTK_BOX(content), hbox);

	GtkWidget *heading = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(heading), "<b>Commit message</b>");
	gtk_label_set_xalign(GTK_LABEL(heading), 0.0);
	gtk_box_append(GTK_BOX(content), heading);

	GtkWidget *sw = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
	                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
	gtk_widget_set_size_request(sw, -1, 100);
	gtk_widget_set_hexpand(sw, TRUE);
	gtk_widget_set_vexpand(sw, TRUE);
	GtkWidget *textview = gtk_text_view_new();
	gtk_text_view_set_editable(GTK_TEXT_VIEW(textview), FALSE);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(textview), FALSE);
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(textview), GTK_WRAP_WORD);
	gtk_text_buffer_set_text(gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview)),
	                         _("Loading commit message..."), -1);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sw), textview);
	gtk_box_append(GTK_BOX(content), sw);

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	gtk_widget_set_margin_top(bbox, 6);
	GtkWidget *btn_cancel = gtk_button_new_with_mnemonic(_("_Cancel"));
	GtkWidget *btn_ok     = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(btn_ok, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), btn_cancel);
	gtk_box_append(GTK_BOX(bbox), btn_ok);
	gtk_box_append(GTK_BOX(content), bbox);

	CommitMessageUpdateData update_data = {
		.scriptpath       = scriptpath,
		.message_textview = GTK_TEXT_VIEW(textview),
	};
	on_script_revision_spin_value_changed(GTK_SPIN_BUTTON(spin), &update_data);
	g_signal_connect(spin, "value-changed",
	                 G_CALLBACK(on_script_revision_spin_value_changed), &update_data);

	rev_dlg_ctx ctx = { g_main_loop_new(NULL, FALSE), GTK_RESPONSE_CANCEL };
	g_signal_connect(btn_cancel, "clicked",       G_CALLBACK(rev_dlg_cancel),         &ctx);
	g_signal_connect(btn_ok,     "clicked",       G_CALLBACK(rev_dlg_ok),             &ctx);
	g_signal_connect(dialog,     "close-request", G_CALLBACK(rev_dlg_close_request),  &ctx);

	gtk_widget_set_visible(dialog, TRUE);
	g_main_loop_run(ctx.loop);

	int revisions_back = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin));
	gint response = ctx.result;
	gtk_window_destroy(GTK_WINDOW(dialog));
	g_main_loop_unref(ctx.loop);

	if (response == GTK_RESPONSE_OK) {
		size_t length = 0;
		gchar *contents = get_script_content_string_from_file_revision(
		    scriptpath, revisions_back, &length, NULL, NULL);
		if (length > 0 && contents) {
			const char *ext = get_filename_ext(scriptpath);
			new_script(contents, length, ext);
			g_free(contents);
		}
	}
}

/* Secondary-click handler on the script GtkColumnView: pulls the
 * scriptpath of the *clicked* row (not the selected one — they aren't
 * the same on right-click) and opens the revision dialog.  The cell
 * factories' bind callbacks attach the SirilScriptRow as
 * "siril-script-row" data on every cell widget; we walk up from the
 * picked widget until we find one. */
static void on_git_columnview_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)n_press; (void)user_data;
	if (!git_columnview) return;
	GtkWidget *picked = gtk_widget_pick(GTK_WIDGET(git_columnview), x, y, GTK_PICK_DEFAULT);
	SirilScriptRow *row = NULL;
	for (GtkWidget *w = picked; w; w = gtk_widget_get_parent(w)) {
		gpointer d = g_object_get_data(G_OBJECT(w), "siril-script-row");
		if (d) { row = SIRIL_SCRIPT_ROW(d); break; }
		if (w == GTK_WIDGET(git_columnview)) break;
	}
	if (!row || !row->scriptpath) return;
	show_script_revision_dialog(row->scriptpath);
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
}

void on_disable_gitscripts() {
	com.pref.use_scripts_repository = FALSE;
	if (list_store) {
		g_list_store_remove_all(list_store);
	}
	gui_repo_scripts_mutex_lock();
	g_slist_free_full(com.repo_scripts, g_free);
	gui_repo_scripts_mutex_unlock();
	com.repo_scripts = NULL;
	if (com.pref.selected_scripts)
		g_slist_free_full(com.pref.selected_scripts, g_free);
	com.pref.selected_scripts = NULL;
	if (com.pref.startup_scripts)
		g_slist_free_full(com.pref.startup_scripts, g_free);
	com.pref.startup_scripts = NULL;
	g_thread_unref(g_thread_new("refresh_script_menu", refresh_script_menu_in_thread, GINT_TO_POINTER(1)));
}

void on_manual_script_sync_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	g_thread_unref(g_thread_new("update_scripts", update_scripts, NULL));
}

void on_manual_spcc_sync_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	g_thread_unref(g_thread_new("update_spcc", update_spcc, NULL));
}

void on_pref_use_gitscripts_toggled(GtkCheckButton *button, gpointer user_data) {
	(void)user_data;
	com.pref.use_scripts_repository = siril_toggle_get_active(GTK_WIDGET(button));
	if (com.pref.use_scripts_repository) {
		g_thread_unref(g_thread_new("update_scripts", initialize_scripts, NULL));
	}
	git_gui_init_statics();
	gtk_widget_set_sensitive(git_pref_auto_updates, com.pref.use_scripts_repository);
	gtk_widget_set_sensitive(git_manual_sync_btn, (com.pref.use_scripts_repository && com.script_repo_available));
	if (git_columnview)
		gtk_widget_set_sensitive(GTK_WIDGET(git_columnview),
				(com.pref.use_scripts_repository && com.script_repo_available));
}

void on_spcc_repo_enable_toggled(GtkCheckButton *button, gpointer user_data) {
	(void)user_data;
	com.pref.spcc.use_spcc_repository = siril_toggle_get_active(GTK_WIDGET(button));
	if (com.pref.spcc.use_spcc_repository) {
		g_thread_unref(g_thread_new("update_spcc", initialize_spcc, NULL));
	}
	git_gui_init_statics();
	gtk_widget_set_sensitive(git_spcc_sync_startup, com.pref.spcc.use_spcc_repository);
	gtk_widget_set_sensitive(git_spcc_manual_sync, (com.pref.spcc.use_spcc_repository && com.spcc_repo_available));
}

/* The GtkCellRendererToggle / GtkTreeView signal handlers wired in the old
 * .ui are gone; keep stub implementations of the public symbols for any
 * lingering reference (linker safety). */
void on_treeview_scripts_row_activated(void) { }
void on_script_list_active_toggled(void) { }
void on_script_list_startup_toggled(void) { }

#else

void hide_git_widgets() {
	GtkWidget *frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "frame_gitscripts"));
	gtk_widget_set_visible(frame, FALSE);
}

void on_pref_use_gitscripts_toggled(GtkCheckButton *button,
									gpointer user_data) {
	(void)button; (void)user_data;
	return;
}

void on_spcc_repo_enable_toggled(GtkCheckButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	return;
}

void on_treeview_scripts_row_activated(void) { }
void on_script_list_active_toggled(void) { }
void on_script_list_startup_toggled(void) { }

void on_script_text_close_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	return;
}

void on_manual_script_sync_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	return;
}

void on_manual_spcc_sync_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	return;
}

#endif
