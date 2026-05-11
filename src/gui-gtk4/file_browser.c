/*
 * SirilFileBrowser: custom modal file picker with preview pane.
 *
 * GTK4's GtkFileDialog dropped the preview slot that GtkFileChooserDialog
 * had — the GTK team's recommendation for apps that need previews is to
 * build their own browser using GtkDirectoryList / GtkListView /
 * GtkSingleSelection / GtkPicture / GdkTexture and a custom thumbnail
 * pipeline.  This file implements that, with back/forward navigation and
 * a sidebar of standard locations.
 */

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui-gtk4/file_browser.h"
#include "gui-gtk4/utils.h"

#include <string.h>

/* From io/image_format_fits.c — produces a small RGB GdkPixbuf for any
 * FITS file plus a textual description string.  Used by the default
 * preview handler so the browser can show FITS thumbnails. */
extern GdkPixbuf *get_thumbnail_from_fits(char *filename, gchar **descr);

typedef struct {
	gchar *title;
	GPtrArray *specs;   /* of GPatternSpec* */
} BrowserFilter;

static void browser_filter_free(BrowserFilter *bf) {
	if (!bf) return;
	g_free(bf->title);
	if (bf->specs) {
		for (guint i = 0; i < bf->specs->len; i++)
			g_pattern_spec_free(g_ptr_array_index(bf->specs, i));
		g_ptr_array_unref(bf->specs);
	}
	g_free(bf);
}

struct _SirilFileBrowser {
	GtkWindow              *parent;
	GtkWindow              *window;
	GtkLabel               *path_label;
	GtkListView            *listview;
	GtkPicture             *preview;
	GtkLabel               *metadata_label;
	GtkDropDown            *filter_combo;
	GtkWidget              *open_button;
	GtkWidget              *back_button;
	GtkWidget              *forward_button;
	GtkWidget              *up_button;

	/* Owned model objects.  Kept in fields so callbacks (filter, set_file,
	 * selection-changed) can reach them; references are taken via
	 * g_object_ref so my fields stay valid even after the widget tree
	 * teardown drops its own refs. */
	GtkDirectoryList       *dir_list;        /* +1 ref */
	GListStore             *recent_store;    /* +1 ref; GFileInfo items, full paths stashed via g_object_set_data */
	GtkFilterListModel     *filter_model;    /* +1 ref */
	GtkCustomFilter        *file_filter;     /* +1 ref */
	GtkSelectionModel      *selection;       /* +1 ref; GtkSingleSelection or GtkMultiSelection */
	gboolean                select_multiple;
	gboolean                in_recent_mode;

	GFile                  *current_folder;     /* owned */
	GFile                  *initial_file;       /* owned, set before run */

	/* History stacks for back / forward.  Hold owned GFile* refs. */
	GPtrArray              *back_history;
	GPtrArray              *forward_history;

	GPtrArray              *filters;            /* of BrowserFilter*, owned */
	gint                    active_filter_idx;  /* index into filters, -1 = none */

	SirilFileBrowserPreview preview_cb;
	gpointer                preview_cb_data;

	/* Result of run */
	gchar                  *result_path;
	GSList                 *result_paths;     /* full set in multi mode */
	gint                    response;
	GMainLoop              *loop;
};

/* ── filter helpers ────────────────────────────────────────────────── */

static gboolean filter_accepts(const BrowserFilter *bf, const gchar *name) {
	if (!bf || !bf->specs || bf->specs->len == 0) return TRUE;
	gsize nlen = strlen(name);
	gchar *reverse = g_utf8_strreverse(name, -1);
	for (guint i = 0; i < bf->specs->len; i++) {
		GPatternSpec *spec = g_ptr_array_index(bf->specs, i);
		if (g_pattern_spec_match(spec, nlen, name, reverse)) {
			g_free(reverse);
			return TRUE;
		}
	}
	g_free(reverse);
	return FALSE;
}

static gboolean fileinfo_filter_func(gpointer item, gpointer user_data) {
	GFileInfo *info = G_FILE_INFO(item);
	SirilFileBrowser *fb = user_data;
	/* Always show directories (so user can navigate into them). */
	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY)
		return TRUE;
	/* Hide hidden files. */
	if (g_file_info_get_is_hidden(info))
		return FALSE;
	if (fb->active_filter_idx < 0 || fb->active_filter_idx >= (gint)fb->filters->len)
		return TRUE;
	BrowserFilter *bf = g_ptr_array_index(fb->filters, fb->active_filter_idx);
	return filter_accepts(bf, g_file_info_get_name(info));
}

/* ── navigation / history ──────────────────────────────────────────── */

static int sidebar_recent_cmp_visited_desc(gconstpointer a, gconstpointer b);
static void update_nav_sensitivity(SirilFileBrowser *fb);

static void build_recent_store(SirilFileBrowser *fb) {
	if (!fb->recent_store)
		fb->recent_store = g_list_store_new(G_TYPE_FILE_INFO);
	else
		g_list_store_remove_all(fb->recent_store);

	GtkRecentManager *mgr = gtk_recent_manager_get_default();
	if (!mgr) return;
	GList *items = g_list_sort(gtk_recent_manager_get_items(mgr),
	                           sidebar_recent_cmp_visited_desc);
	int count = 0;
	const int max_items = 50;
	for (GList *l = items; l && count < max_items; l = l->next) {
		GtkRecentInfo *ri = l->data;
		const gchar *uri = gtk_recent_info_get_uri(ri);
		if (!uri) continue;
		gchar *path = g_filename_from_uri(uri, NULL, NULL);
		if (!path || !g_file_test(path, G_FILE_TEST_EXISTS)) {
			g_free(path);
			continue;
		}
		if (get_type_from_filename(path) == TYPEUNDEF) {
			g_free(path);
			continue;
		}
		gchar *basename = g_path_get_basename(path);
		GFileInfo *fi = g_file_info_new();
		g_file_info_set_name(fi, basename);
		g_file_info_set_display_name(fi, basename);
		g_file_info_set_file_type(fi, G_FILE_TYPE_REGULAR);
		/* Without this attribute, g_file_info_get_is_hidden() in
		 * fileinfo_filter_func trips a GLib-GIO critical (GTK 4.14+). */
		g_file_info_set_attribute_boolean(fi,
			G_FILE_ATTRIBUTE_STANDARD_IS_HIDDEN, FALSE);
		g_object_set_data_full(G_OBJECT(fi), "siril-full-path",
		                       g_strdup(path), g_free);
		g_list_store_append(fb->recent_store, fi);
		g_object_unref(fi);
		g_free(basename);
		g_free(path);
		count++;
	}
	g_list_free_full(items, (GDestroyNotify) gtk_recent_info_unref);
}

static void exit_recent_mode(SirilFileBrowser *fb) {
	if (!fb->in_recent_mode) return;
	fb->in_recent_mode = FALSE;
	gtk_filter_list_model_set_model(fb->filter_model,
		G_LIST_MODEL(fb->dir_list));
}

static void enter_recent_mode(SirilFileBrowser *fb) {
	if (fb->in_recent_mode) return;
	/* Push current folder onto back history; Back from Recent mode then
	 * pops it and returns to the previous directory (which exits Recent
	 * mode via apply_current_folder). */
	if (fb->current_folder)
		g_ptr_array_add(fb->back_history, g_object_ref(fb->current_folder));
	g_ptr_array_set_size(fb->forward_history, 0);
	build_recent_store(fb);
	fb->in_recent_mode = TRUE;
	gtk_filter_list_model_set_model(fb->filter_model,
		G_LIST_MODEL(fb->recent_store));
	gtk_label_set_text(fb->path_label, _("Recent files"));
	update_nav_sensitivity(fb);
}

static void update_path_label(SirilFileBrowser *fb) {
	if (!fb->current_folder) return;
	gchar *path = g_file_get_path(fb->current_folder);
	gtk_label_set_text(fb->path_label, path ? path : "");
	g_free(path);
}

static void update_nav_sensitivity(SirilFileBrowser *fb) {
	if (fb->back_button)
		gtk_widget_set_sensitive(fb->back_button, fb->back_history->len > 0);
	if (fb->forward_button)
		gtk_widget_set_sensitive(fb->forward_button, fb->forward_history->len > 0);
	if (fb->up_button) {
		gboolean can_up = FALSE;
		if (!fb->in_recent_mode && fb->current_folder) {
			GFile *p = g_file_get_parent(fb->current_folder);
			can_up = p != NULL;
			if (p) g_object_unref(p);
		}
		gtk_widget_set_sensitive(fb->up_button, can_up);
	}
}

/* Internal helper: navigate to `folder` without touching history.  All
 * external callers go through navigate_to(), which updates history. */
static void apply_current_folder(SirilFileBrowser *fb, GFile *folder) {
	if (!folder) return;
	exit_recent_mode(fb);  /* any folder navigation cancels Recent view */
	if (fb->current_folder) g_object_unref(fb->current_folder);
	fb->current_folder = g_object_ref(folder);
	gtk_directory_list_set_file(fb->dir_list, folder);
	update_path_label(fb);
	update_nav_sensitivity(fb);
}

/* User-facing navigation: push current onto back, clear forward. */
static void navigate_to(SirilFileBrowser *fb, GFile *folder) {
	if (!folder) return;
	if (fb->current_folder && g_file_equal(folder, fb->current_folder))
		return;
	if (fb->current_folder)
		g_ptr_array_add(fb->back_history, g_object_ref(fb->current_folder));
	g_ptr_array_set_size(fb->forward_history, 0);
	apply_current_folder(fb, folder);
}

static void on_up_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (!fb->current_folder) return;
	GFile *parent = g_file_get_parent(fb->current_folder);
	if (parent) {
		navigate_to(fb, parent);
		g_object_unref(parent);
	}
}

static void on_back_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (fb->back_history->len == 0) return;
	GFile *prev = g_ptr_array_steal_index(fb->back_history, fb->back_history->len - 1);
	if (fb->current_folder)
		g_ptr_array_add(fb->forward_history, g_object_ref(fb->current_folder));
	apply_current_folder(fb, prev);
	g_object_unref(prev);
}

static void on_forward_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (fb->forward_history->len == 0) return;
	GFile *next = g_ptr_array_steal_index(fb->forward_history, fb->forward_history->len - 1);
	if (fb->current_folder)
		g_ptr_array_add(fb->back_history, g_object_ref(fb->current_folder));
	apply_current_folder(fb, next);
	g_object_unref(next);
}

/* ── list factory ──────────────────────────────────────────────────── */

static void on_setup_item(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	GtkWidget *icon = gtk_image_new();
	GtkWidget *label = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_widget_set_hexpand(label, TRUE);
	gtk_box_append(GTK_BOX(box), icon);
	gtk_box_append(GTK_BOX(box), label);
	gtk_list_item_set_child(item, box);
}

static void on_bind_item(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GFileInfo *info = G_FILE_INFO(gtk_list_item_get_item(item));
	GtkWidget *box = gtk_list_item_get_child(item);
	GtkWidget *icon = gtk_widget_get_first_child(box);
	GtkWidget *label = gtk_widget_get_next_sibling(icon);

	const char *icon_name =
		(g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY)
			? "folder-symbolic" : "text-x-generic-symbolic";
	gtk_image_set_from_icon_name(GTK_IMAGE(icon), icon_name);
	gtk_label_set_text(GTK_LABEL(label), g_file_info_get_display_name(info));
}

/* ── selection / activation ────────────────────────────────────────── */

/* Returns the first selected item (or NULL).  Works for both
 * GtkSingleSelection and GtkMultiSelection. */
static GFileInfo *first_selected_info(SirilFileBrowser *fb) {
	if (!fb->selection) return NULL;
	if (!fb->select_multiple) {
		GObject *item = gtk_single_selection_get_selected_item(
			GTK_SINGLE_SELECTION(fb->selection));
		return (item && G_IS_FILE_INFO(item)) ? G_FILE_INFO(item) : NULL;
	}
	GtkBitset *sel = gtk_selection_model_get_selection(fb->selection);
	if (!sel) return NULL;
	GFileInfo *out = NULL;
	if (!gtk_bitset_is_empty(sel)) {
		guint idx = gtk_bitset_get_minimum(sel);
		GObject *item = g_list_model_get_item(G_LIST_MODEL(fb->selection), idx);
		if (item && G_IS_FILE_INFO(item)) out = G_FILE_INFO(item);
		if (item) g_object_unref(item);
	}
	gtk_bitset_unref(sel);
	return out;
}

static gchar *path_from_info(SirilFileBrowser *fb, GFileInfo *info) {
	if (!info) return NULL;
	/* Items injected by Recent mode carry the full absolute path. */
	const char *abs = g_object_get_data(G_OBJECT(info), "siril-full-path");
	if (abs) return g_strdup(abs);
	if (!fb->current_folder) return NULL;
	const char *name = g_file_info_get_name(info);
	if (!name) return NULL;
	GFile *gf = g_file_get_child(fb->current_folder, name);
	gchar *p = g_file_get_path(gf);
	g_object_unref(gf);
	return p;
}

static gchar *current_selected_path(SirilFileBrowser *fb) {
	return path_from_info(fb, first_selected_info(fb));
}

static gboolean current_is_directory(SirilFileBrowser *fb) {
	GFileInfo *info = first_selected_info(fb);
	if (!info) return FALSE;
	return g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY;
}

static gboolean any_selected(SirilFileBrowser *fb) {
	if (!fb->selection) return FALSE;
	if (!fb->select_multiple)
		return gtk_single_selection_get_selected_item(
			GTK_SINGLE_SELECTION(fb->selection)) != NULL;
	GtkBitset *sel = gtk_selection_model_get_selection(fb->selection);
	gboolean has = sel && !gtk_bitset_is_empty(sel);
	if (sel) gtk_bitset_unref(sel);
	return has;
}

static void update_preview_for_selection(SirilFileBrowser *fb) {
	gchar *path = (current_is_directory(fb)) ? NULL : current_selected_path(fb);
	SirilFileBrowserPreview cb = fb->preview_cb ? fb->preview_cb
	                                            : siril_file_browser_default_preview;
	cb(path, fb->preview, fb->metadata_label,
	   fb->preview_cb ? fb->preview_cb_data : fb);
	g_free(path);
}

static void on_selection_changed(GtkSelectionModel *m, guint pos, guint n, gpointer ud) {
	(void)m; (void)pos; (void)n;
	SirilFileBrowser *fb = ud;
	update_preview_for_selection(fb);
	gtk_widget_set_sensitive(fb->open_button, any_selected(fb));
}

static void on_row_activated(GtkListView *lv, guint position, gpointer ud) {
	(void)lv; (void)position;
	SirilFileBrowser *fb = ud;
	GFileInfo *info = first_selected_info(fb);
	if (!info) return;
	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY) {
		const char *name = g_file_info_get_name(info);
		if (!name || !fb->current_folder) return;
		GFile *target = g_file_get_child(fb->current_folder, name);
		navigate_to(fb, target);
		g_object_unref(target);
	} else {
		gchar *p = path_from_info(fb, info);
		if (!p) return;
		g_free(fb->result_path);
		fb->result_path = p;
		fb->response = GTK_RESPONSE_ACCEPT;
		if (fb->loop) g_main_loop_quit(fb->loop);
	}
}

/* ── filter dropdown ───────────────────────────────────────────────── */

static void on_filter_changed(GObject *obj, GParamSpec *pspec, gpointer ud) {
	(void)pspec;
	SirilFileBrowser *fb = ud;
	fb->active_filter_idx = (gint)gtk_drop_down_get_selected(GTK_DROP_DOWN(obj));
	gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);
}

/* ── buttons ───────────────────────────────────────────────────────── */

static void on_cancel_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	fb->response = GTK_RESPONSE_CANCEL;
	if (fb->loop) g_main_loop_quit(fb->loop);
}

static void on_open_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (current_is_directory(fb)) {
		on_row_activated(fb->listview, 0, fb);
		return;
	}
	/* Stash result(s). */
	g_slist_free_full(fb->result_paths, g_free);
	fb->result_paths = NULL;
	g_free(fb->result_path);
	fb->result_path = NULL;
	if (fb->select_multiple) {
		GtkBitset *sel = gtk_selection_model_get_selection(fb->selection);
		if (sel) {
			GtkBitsetIter it;
			guint idx;
			if (gtk_bitset_iter_init_first(&it, sel, &idx)) {
				do {
					GObject *item = g_list_model_get_item(
						G_LIST_MODEL(fb->selection), idx);
					if (item && G_IS_FILE_INFO(item)) {
						gchar *p = path_from_info(fb, G_FILE_INFO(item));
						if (p) {
							if (!fb->result_path) fb->result_path = g_strdup(p);
							fb->result_paths = g_slist_append(fb->result_paths, p);
						}
					}
					if (item) g_object_unref(item);
				} while (gtk_bitset_iter_next(&it, &idx));
			}
			gtk_bitset_unref(sel);
		}
		if (!fb->result_path) return;  /* nothing selected */
	} else {
		gchar *p = current_selected_path(fb);
		if (!p) return;
		fb->result_path = p;
	}
	fb->response = GTK_RESPONSE_ACCEPT;
	if (fb->loop) g_main_loop_quit(fb->loop);
}

static gboolean on_close_request(GtkWindow *w, gpointer ud) {
	(void)w;
	SirilFileBrowser *fb = ud;
	fb->response = GTK_RESPONSE_CANCEL;
	if (fb->loop) g_main_loop_quit(fb->loop);
	return TRUE;
}

/* ── sidebar (standard locations + recent files + mounted volumes) ─── */

static void on_sidebar_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer ud) {
	(void)box;
	SirilFileBrowser *fb = ud;
	/* The single "Recent" entry: swap the main list to recent-files mode. */
	if (g_object_get_data(G_OBJECT(row), "siril-is-recent")) {
		enter_recent_mode(fb);
		return;
	}
	/* Folder row: navigate into it. */
	const char *path = g_object_get_data(G_OBJECT(row), "siril-path");
	if (!path) return;
	GFile *gf = g_file_new_for_path(path);
	navigate_to(fb, gf);
	g_object_unref(gf);
}

static GtkWidget *sidebar_make_row(const char *icon_name, const char *label,
                                   const char *secondary) {
	GtkWidget *row = gtk_list_box_row_new();
	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_widget_set_margin_start(hbox, 6);
	gtk_widget_set_margin_end(hbox, 6);
	gtk_widget_set_margin_top(hbox, 4);
	gtk_widget_set_margin_bottom(hbox, 4);
	GtkWidget *icon = gtk_image_new_from_icon_name(icon_name);
	GtkWidget *labels = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	GtkWidget *lbl = gtk_label_new(label);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_label_set_ellipsize(GTK_LABEL(lbl), PANGO_ELLIPSIZE_END);
	gtk_widget_set_hexpand(labels, TRUE);
	gtk_box_append(GTK_BOX(labels), lbl);
	if (secondary && *secondary) {
		GtkWidget *sub = gtk_label_new(secondary);
		gtk_label_set_xalign(GTK_LABEL(sub), 0.0);
		gtk_label_set_ellipsize(GTK_LABEL(sub), PANGO_ELLIPSIZE_START);
		gtk_widget_add_css_class(sub, "dim-label");
		PangoAttrList *attrs = pango_attr_list_new();
		pango_attr_list_insert(attrs, pango_attr_scale_new(PANGO_SCALE_SMALL));
		gtk_label_set_attributes(GTK_LABEL(sub), attrs);
		pango_attr_list_unref(attrs);
		gtk_box_append(GTK_BOX(labels), sub);
	}
	gtk_box_append(GTK_BOX(hbox), icon);
	gtk_box_append(GTK_BOX(hbox), labels);
	gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), hbox);
	return row;
}

static void sidebar_add_location(GtkListBox *box, const char *icon_name,
                                 const char *label, const char *path) {
	if (!path || !*path) return;
	if (!g_file_test(path, G_FILE_TEST_IS_DIR)) return;
	GtkWidget *row = sidebar_make_row(icon_name, label, NULL);
	g_object_set_data_full(G_OBJECT(row), "siril-path", g_strdup(path), g_free);
	gtk_list_box_append(box, row);
}

static void sidebar_add_recent_entry(GtkListBox *box) {
	GtkWidget *row = sidebar_make_row("document-open-recent-symbolic",
	                                  _("Recent"), NULL);
	g_object_set_data(G_OBJECT(row), "siril-is-recent", GINT_TO_POINTER(1));
	gtk_list_box_append(box, row);
}

static void sidebar_add_section_header(GtkListBox *box, const char *text) {
	GtkWidget *row = gtk_list_box_row_new();
	gtk_list_box_row_set_selectable(GTK_LIST_BOX_ROW(row), FALSE);
	gtk_list_box_row_set_activatable(GTK_LIST_BOX_ROW(row), FALSE);
	GtkWidget *lbl = gtk_label_new(text);
	gtk_widget_add_css_class(lbl, "dim-label");
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_widget_set_margin_top(lbl, 8);
	gtk_widget_set_margin_bottom(lbl, 2);
	gtk_widget_set_margin_start(lbl, 6);
	PangoAttrList *attrs = pango_attr_list_new();
	pango_attr_list_insert(attrs, pango_attr_weight_new(PANGO_WEIGHT_BOLD));
	pango_attr_list_insert(attrs, pango_attr_scale_new(PANGO_SCALE_SMALL));
	gtk_label_set_attributes(GTK_LABEL(lbl), attrs);
	pango_attr_list_unref(attrs);
	gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), lbl);
	gtk_list_box_append(box, row);
}

static int sidebar_recent_cmp_visited_desc(gconstpointer a, gconstpointer b) {
	GDateTime *ta = gtk_recent_info_get_visited((GtkRecentInfo *) a);
	GDateTime *tb = gtk_recent_info_get_visited((GtkRecentInfo *) b);
	if (!ta && !tb) return 0;
	if (!ta) return 1;
	if (!tb) return -1;
	return g_date_time_compare(tb, ta);  /* descending */
}

static int sidebar_populate_volumes(GtkListBox *box) {
	GVolumeMonitor *monitor = g_volume_monitor_get();
	if (!monitor) return 0;
	int count = 0;
	GList *mounts = g_volume_monitor_get_mounts(monitor);
	for (GList *l = mounts; l; l = l->next) {
		GMount *mount = G_MOUNT(l->data);
		if (g_mount_is_shadowed(mount)) continue;
		gchar *name = g_mount_get_name(mount);
		GFile *root = g_mount_get_root(mount);
		gchar *path = root ? g_file_get_path(root) : NULL;
		if (path && g_file_test(path, G_FILE_TEST_IS_DIR)) {
			GIcon *gicon = g_mount_get_icon(mount);
			const char *icon_name = "drive-harddisk-symbolic";
			if (G_IS_THEMED_ICON(gicon)) {
				const gchar * const *names = g_themed_icon_get_names(G_THEMED_ICON(gicon));
				if (names && names[0]) icon_name = names[0];
			}
			sidebar_add_location(box, icon_name, name, path);
			count++;
			if (gicon) g_object_unref(gicon);
		}
		g_free(path);
		g_free(name);
		if (root) g_object_unref(root);
	}
	g_list_free_full(mounts, g_object_unref);
	g_object_unref(monitor);
	return count;
}

static GtkWidget *build_sidebar(SirilFileBrowser *fb) {
	GtkWidget *list = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(list), GTK_SELECTION_NONE);
	gtk_widget_add_css_class(list, "navigation-sidebar");

	/* Section: Standard places.  "Recent" is the first entry — clicking it
	 * swaps the main list to show recent files (instead of navigating to
	 * a directory). */
	sidebar_add_section_header(GTK_LIST_BOX(list), _("Places"));
	sidebar_add_recent_entry(GTK_LIST_BOX(list));
	sidebar_add_location(GTK_LIST_BOX(list), "user-home-symbolic",
	                     _("Home"), g_get_home_dir());
	sidebar_add_location(GTK_LIST_BOX(list), "user-desktop-symbolic",
	                     _("Desktop"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_DESKTOP));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-documents-symbolic",
	                     _("Documents"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_DOCUMENTS));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-download-symbolic",
	                     _("Downloads"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_DOWNLOAD));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-pictures-symbolic",
	                     _("Pictures"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_PICTURES));
	if (com.wd && *com.wd && g_file_test(com.wd, G_FILE_TEST_IS_DIR))
		sidebar_add_location(GTK_LIST_BOX(list), "folder-symbolic",
		                     _("Working Directory"), com.wd);

	/* Section: Mounted volumes (only if any are present). */
	{
		gint header_pos = 0;
		GtkWidget *r;
		while ((r = GTK_WIDGET(gtk_list_box_get_row_at_index(GTK_LIST_BOX(list), header_pos))))
			header_pos++;
		sidebar_add_section_header(GTK_LIST_BOX(list), _("Devices"));
		if (sidebar_populate_volumes(GTK_LIST_BOX(list)) == 0) {
			GtkListBoxRow *hdr = gtk_list_box_get_row_at_index(GTK_LIST_BOX(list), header_pos);
			if (hdr) gtk_list_box_remove(GTK_LIST_BOX(list), GTK_WIDGET(hdr));
		}
	}

	g_signal_connect(list, "row-activated",
	                 G_CALLBACK(on_sidebar_row_activated), fb);

	GtkWidget *scrolled = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled), list);
	gtk_widget_set_size_request(scrolled, 180, -1);
	return scrolled;
}

/* ── default preview handler ───────────────────────────────────────── */

void siril_file_browser_default_preview(const gchar *path,
                                        GtkPicture  *picture,
                                        GtkLabel    *metadata_label,
                                        gpointer     user_data) {
	(void)user_data;
	gtk_picture_set_paintable(picture, NULL);
	if (metadata_label) gtk_label_set_text(metadata_label, "");
	if (!path) return;

	/* FITS first — recognise common extensions and use Siril's loader. */
	gchar *low = g_ascii_strdown(path, -1);
	gboolean is_fits = g_str_has_suffix(low, ".fit") ||
	                   g_str_has_suffix(low, ".fits") ||
	                   g_str_has_suffix(low, ".fts") ||
	                   g_str_has_suffix(low, ".fit.fz") ||
	                   g_str_has_suffix(low, ".fits.fz") ||
	                   g_str_has_suffix(low, ".fts.fz");
	g_free(low);

	if (is_fits) {
		gchar *descr = NULL;
		GdkPixbuf *pb = get_thumbnail_from_fits((char *)path, &descr);
		if (pb) {
			GdkTexture *tex = siril_texture_from_pixbuf(pb);
			if (tex) {
				gtk_picture_set_paintable(picture, GDK_PAINTABLE(tex));
				g_object_unref(tex);
			}
			g_object_unref(pb);
		}
		if (metadata_label && descr)
			gtk_label_set_text(metadata_label, descr);
		g_free(descr);
		return;
	}

	/* Try GdkTexture for raster formats (PNG/JPEG/TIFF/WebP/...). */
	GError *err = NULL;
	GFile *gf = g_file_new_for_path(path);
	GdkTexture *tex = gdk_texture_new_from_file(gf, &err);
	g_object_unref(gf);
	if (tex) {
		gtk_picture_set_paintable(picture, GDK_PAINTABLE(tex));
		if (metadata_label) {
			gchar *meta = g_strdup_printf("%d x %d",
				gdk_texture_get_width(tex), gdk_texture_get_height(tex));
			gtk_label_set_text(metadata_label, meta);
			g_free(meta);
		}
		g_object_unref(tex);
	} else {
		g_clear_error(&err);
	}
}

/* ── construction / API ────────────────────────────────────────────── */

SirilFileBrowser *siril_file_browser_new(GtkWindow *parent, const gchar *title) {
	SirilFileBrowser *fb = g_new0(SirilFileBrowser, 1);
	fb->parent  = parent;
	fb->filters = g_ptr_array_new_with_free_func((GDestroyNotify)browser_filter_free);
	fb->back_history    = g_ptr_array_new_with_free_func(g_object_unref);
	fb->forward_history = g_ptr_array_new_with_free_func(g_object_unref);
	fb->active_filter_idx = -1;
	fb->response = GTK_RESPONSE_NONE;

	fb->window = GTK_WINDOW(gtk_window_new());
	gtk_window_set_title(fb->window, title ? title : _("Open File"));
	gtk_window_set_modal(fb->window, TRUE);
	if (parent) gtk_window_set_transient_for(fb->window, parent);
	gtk_window_set_default_size(fb->window, 1000, 620);

	/* Toolbar: Back / Forward / Up + path label. */
	GtkWidget *toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	gtk_widget_set_margin_start(toolbar, 6);
	gtk_widget_set_margin_end(toolbar, 6);
	gtk_widget_set_margin_top(toolbar, 6);
	gtk_widget_set_margin_bottom(toolbar, 3);
	fb->back_button = gtk_button_new_from_icon_name("go-previous-symbolic");
	gtk_widget_set_tooltip_text(fb->back_button, _("Back"));
	g_signal_connect(fb->back_button, "clicked", G_CALLBACK(on_back_clicked), fb);
	fb->forward_button = gtk_button_new_from_icon_name("go-next-symbolic");
	gtk_widget_set_tooltip_text(fb->forward_button, _("Forward"));
	g_signal_connect(fb->forward_button, "clicked", G_CALLBACK(on_forward_clicked), fb);
	fb->up_button = gtk_button_new_from_icon_name("go-up-symbolic");
	gtk_widget_set_tooltip_text(fb->up_button, _("Up a level"));
	g_signal_connect(fb->up_button, "clicked", G_CALLBACK(on_up_clicked), fb);
	gtk_box_append(GTK_BOX(toolbar), fb->back_button);
	gtk_box_append(GTK_BOX(toolbar), fb->forward_button);
	gtk_box_append(GTK_BOX(toolbar), fb->up_button);
	fb->path_label = GTK_LABEL(gtk_label_new(""));
	gtk_label_set_xalign(fb->path_label, 0.0);
	gtk_label_set_ellipsize(fb->path_label, PANGO_ELLIPSIZE_START);
	gtk_widget_set_hexpand(GTK_WIDGET(fb->path_label), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(fb->path_label), 6);
	gtk_box_append(GTK_BOX(toolbar), GTK_WIDGET(fb->path_label));

	/* Models: DirectoryList → FilterListModel → SingleSelection.
	 *
	 * Reference ownership: every chain step takes a transfer-full model
	 * arg, and we want to keep our own ref too — so we g_object_ref each
	 * field that we hand off, and we g_object_ref the selection again
	 * when handing it to gtk_list_view_new (which also takes transfer-full).
	 */
	fb->dir_list     = gtk_directory_list_new("standard::*", NULL);
	fb->file_filter  = gtk_custom_filter_new(fileinfo_filter_func, fb, NULL);
	fb->filter_model = gtk_filter_list_model_new(
		G_LIST_MODEL(g_object_ref(fb->dir_list)),
		GTK_FILTER(g_object_ref(fb->file_filter)));
	{
		GtkSingleSelection *ss = gtk_single_selection_new(
			G_LIST_MODEL(g_object_ref(fb->filter_model)));
		gtk_single_selection_set_autoselect(ss, FALSE);
		gtk_single_selection_set_can_unselect(ss, TRUE);
		fb->selection = GTK_SELECTION_MODEL(ss);
	}

	/* ListView. */
	GtkListItemFactory *factory = gtk_signal_list_item_factory_new();
	g_signal_connect(factory, "setup", G_CALLBACK(on_setup_item), NULL);
	g_signal_connect(factory, "bind",  G_CALLBACK(on_bind_item),  NULL);
	fb->listview = GTK_LIST_VIEW(
		gtk_list_view_new(GTK_SELECTION_MODEL(g_object_ref(fb->selection)),
		                  factory));
	gtk_list_view_set_single_click_activate(fb->listview, FALSE);
	g_signal_connect(fb->listview, "activate", G_CALLBACK(on_row_activated), fb);
	g_signal_connect(fb->selection, "selection-changed",
	                 G_CALLBACK(on_selection_changed), fb);

	GtkWidget *list_scrolled = gtk_scrolled_window_new();
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(list_scrolled),
	                              GTK_WIDGET(fb->listview));
	gtk_widget_set_hexpand(list_scrolled, TRUE);
	gtk_widget_set_vexpand(list_scrolled, TRUE);

	/* Preview pane. */
	GtkWidget *preview_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start(preview_box, 6);
	gtk_widget_set_margin_end(preview_box, 6);
	fb->preview = GTK_PICTURE(gtk_picture_new());
	gtk_widget_set_size_request(GTK_WIDGET(fb->preview), 280, 280);
	gtk_picture_set_can_shrink(fb->preview, TRUE);
	gtk_picture_set_content_fit(fb->preview, GTK_CONTENT_FIT_CONTAIN);
	fb->metadata_label = GTK_LABEL(gtk_label_new(""));
	gtk_label_set_wrap(fb->metadata_label, TRUE);
	gtk_label_set_xalign(fb->metadata_label, 0.0);
	gtk_label_set_selectable(fb->metadata_label, TRUE);
	gtk_box_append(GTK_BOX(preview_box), GTK_WIDGET(fb->preview));
	gtk_box_append(GTK_BOX(preview_box), GTK_WIDGET(fb->metadata_label));

	/* Inner paned: list | preview. */
	GtkWidget *list_preview_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_paned_set_start_child(GTK_PANED(list_preview_paned), list_scrolled);
	gtk_paned_set_end_child(GTK_PANED(list_preview_paned), preview_box);
	gtk_paned_set_position(GTK_PANED(list_preview_paned), 480);

	/* Outer paned: sidebar | (list | preview). */
	GtkWidget *outer_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_paned_set_start_child(GTK_PANED(outer_paned), build_sidebar(fb));
	gtk_paned_set_end_child(GTK_PANED(outer_paned), list_preview_paned);
	gtk_paned_set_position(GTK_PANED(outer_paned), 180);
	gtk_paned_set_resize_start_child(GTK_PANED(outer_paned), FALSE);
	gtk_widget_set_vexpand(outer_paned, TRUE);

	/* Filter dropdown (initially empty; populated as filters are added). */
	fb->filter_combo = GTK_DROP_DOWN(gtk_drop_down_new_from_strings((const char *[]){NULL}));
	gtk_widget_set_visible(GTK_WIDGET(fb->filter_combo), FALSE);
	g_signal_connect(fb->filter_combo, "notify::selected",
	                 G_CALLBACK(on_filter_changed), fb);

	/* Action row. */
	GtkWidget *action_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start(action_row, 6);
	gtk_widget_set_margin_end(action_row, 6);
	gtk_widget_set_margin_bottom(action_row, 6);
	gtk_widget_set_margin_top(action_row, 3);
	gtk_box_append(GTK_BOX(action_row), GTK_WIDGET(fb->filter_combo));
	GtkWidget *spacer = gtk_label_new(NULL);
	gtk_widget_set_hexpand(spacer, TRUE);
	gtk_box_append(GTK_BOX(action_row), spacer);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	g_signal_connect(cancel, "clicked", G_CALLBACK(on_cancel_clicked), fb);
	gtk_box_append(GTK_BOX(action_row), cancel);
	fb->open_button = gtk_button_new_with_label(_("Open"));
	gtk_widget_add_css_class(fb->open_button, "suggested-action");
	gtk_widget_set_sensitive(fb->open_button, FALSE);
	g_signal_connect(fb->open_button, "clicked", G_CALLBACK(on_open_clicked), fb);
	gtk_box_append(GTK_BOX(action_row), fb->open_button);

	/* Assemble. */
	GtkWidget *root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_box_append(GTK_BOX(root), toolbar);
	gtk_box_append(GTK_BOX(root), outer_paned);
	gtk_box_append(GTK_BOX(root), action_row);
	gtk_window_set_child(fb->window, root);

	g_signal_connect(fb->window, "close-request",
	                 G_CALLBACK(on_close_request), fb);

	update_nav_sensitivity(fb);
	return fb;
}

void siril_file_browser_set_initial_folder(SirilFileBrowser *fb, const gchar *path) {
	if (!fb || !path) return;
	GFile *gf = g_file_new_for_path(path);
	apply_current_folder(fb, gf);
	g_object_unref(gf);
}

void siril_file_browser_set_initial_file(SirilFileBrowser *fb, const gchar *path) {
	if (!fb || !path) return;
	g_clear_object(&fb->initial_file);
	fb->initial_file = g_file_new_for_path(path);
	GFile *parent = g_file_get_parent(fb->initial_file);
	if (parent) {
		apply_current_folder(fb, parent);
		g_object_unref(parent);
	}
}

void siril_file_browser_add_filter_pattern(SirilFileBrowser *fb,
                                           const gchar *title,
                                           const gchar *pattern,
                                           gboolean set_default) {
	if (!fb || !title || !pattern) return;
	BrowserFilter *bf = g_new0(BrowserFilter, 1);
	bf->title = g_strdup(title);
	bf->specs = g_ptr_array_new();
	gchar **parts = g_strsplit(pattern, ";", -1);
	for (gint i = 0; parts[i]; i++) {
		if (parts[i][0])
			g_ptr_array_add(bf->specs, g_pattern_spec_new(parts[i]));
	}
	g_strfreev(parts);
	g_ptr_array_add(fb->filters, bf);

	/* Rebuild the dropdown's string list. */
	const char **strs = g_new0(const char *, fb->filters->len + 1);
	for (guint i = 0; i < fb->filters->len; i++) {
		BrowserFilter *f = g_ptr_array_index(fb->filters, i);
		strs[i] = f->title;
	}
	GtkStringList *sl = gtk_string_list_new(strs);
	g_free(strs);
	gtk_drop_down_set_model(fb->filter_combo, G_LIST_MODEL(sl));
	g_object_unref(sl);
	gtk_widget_set_visible(GTK_WIDGET(fb->filter_combo), TRUE);

	if (set_default || fb->active_filter_idx < 0) {
		fb->active_filter_idx = (gint)fb->filters->len - 1;
		gtk_drop_down_set_selected(fb->filter_combo, (guint)fb->active_filter_idx);
	}
	gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);
}

void siril_file_browser_set_preview_callback(SirilFileBrowser       *fb,
                                             SirilFileBrowserPreview cb,
                                             gpointer                user_data) {
	if (!fb) return;
	fb->preview_cb      = cb;
	fb->preview_cb_data = user_data;
}

void siril_file_browser_set_select_multiple(SirilFileBrowser *fb, gboolean multi) {
	if (!fb || multi == fb->select_multiple) return;
	fb->select_multiple = multi;

	/* Disconnect the existing selection-changed handler and drop our ref
	 * on the old selection model.  The listview still holds a ref, which
	 * we override below by handing it the new model. */
	g_signal_handlers_disconnect_by_func(fb->selection,
		G_CALLBACK(on_selection_changed), fb);
	GtkSelectionModel *fresh;
	if (multi) {
		fresh = GTK_SELECTION_MODEL(gtk_multi_selection_new(
			G_LIST_MODEL(g_object_ref(fb->filter_model))));
	} else {
		GtkSingleSelection *ss = gtk_single_selection_new(
			G_LIST_MODEL(g_object_ref(fb->filter_model)));
		gtk_single_selection_set_autoselect(ss, FALSE);
		gtk_single_selection_set_can_unselect(ss, TRUE);
		fresh = GTK_SELECTION_MODEL(ss);
	}
	g_clear_object(&fb->selection);
	fb->selection = fresh;  /* takes the +1 ref we got from _new */
	gtk_list_view_set_model(fb->listview, fb->selection);
	g_signal_connect(fb->selection, "selection-changed",
	                 G_CALLBACK(on_selection_changed), fb);
}

gint siril_file_browser_run(SirilFileBrowser *fb) {
	if (!fb) return GTK_RESPONSE_NONE;
	if (!fb->current_folder) {
		const gchar *home = g_get_home_dir();
		if (home) siril_file_browser_set_initial_folder(fb, home);
	}
	fb->response = GTK_RESPONSE_NONE;
	g_clear_pointer(&fb->result_path, g_free);
	g_slist_free_full(fb->result_paths, g_free);
	fb->result_paths = NULL;
	fb->loop = g_main_loop_new(NULL, FALSE);
	gtk_window_present(fb->window);
	g_main_loop_run(fb->loop);
	g_main_loop_unref(fb->loop);
	fb->loop = NULL;
	gtk_widget_set_visible(GTK_WIDGET(fb->window), FALSE);
	return fb->response;
}

gchar *siril_file_browser_get_path(SirilFileBrowser *fb) {
	if (!fb || !fb->result_path) return NULL;
	return g_strdup(fb->result_path);
}

GSList *siril_file_browser_get_paths(SirilFileBrowser *fb) {
	if (!fb || fb->response != GTK_RESPONSE_ACCEPT) return NULL;
	GSList *out = NULL;
	for (GSList *l = fb->result_paths; l; l = l->next)
		out = g_slist_append(out, g_strdup((const gchar *)l->data));
	if (!out && fb->result_path)
		out = g_slist_append(out, g_strdup(fb->result_path));
	return out;
}

/* ── inline image-picker button helper ─────────────────────────────── */

typedef struct {
	gchar                     *title;
	gchar                     *filter_title;
	gchar                     *filter_pattern;
	SirilImageButtonCallback   on_picked;
	gpointer                   user_data;
} SirilImageButtonSpec;

static void siril_image_button_spec_free(gpointer p) {
	SirilImageButtonSpec *s = p;
	if (!s) return;
	g_free(s->title);
	g_free(s->filter_title);
	g_free(s->filter_pattern);
	g_free(s);
}

static void on_image_button_clicked(GtkButton *btn, gpointer user_data) {
	SirilImageButtonSpec *spec = user_data;
	GtkRoot *root = gtk_widget_get_root(GTK_WIDGET(btn));
	GtkWindow *parent = GTK_IS_WINDOW(root) ? GTK_WINDOW(root) : NULL;
	SirilFileBrowser *fb = siril_file_browser_new(parent,
		spec->title ? spec->title : _("Open Image"));
	if (spec->filter_title && spec->filter_pattern)
		siril_file_browser_add_filter_pattern(fb,
			spec->filter_title, spec->filter_pattern, TRUE);
	const gchar *prev = g_object_get_data(G_OBJECT(btn), "siril-path");
	const gchar *hint = g_object_get_data(G_OBJECT(btn), "siril-initial-folder");
	if (prev && *prev)
		siril_file_browser_set_initial_file(fb, prev);
	else if (hint && *hint)
		siril_file_browser_set_initial_folder(fb, hint);
	else if (com.wd && g_file_test(com.wd, G_FILE_TEST_IS_DIR))
		siril_file_browser_set_initial_folder(fb, com.wd);

	if (siril_file_browser_run(fb) == GTK_RESPONSE_ACCEPT) {
		gchar *path = siril_file_browser_get_path(fb);
		if (path) {
			g_object_set_data_full(G_OBJECT(btn), "siril-path",
			                       g_strdup(path), g_free);
			gchar *base = g_path_get_basename(path);
			gtk_button_set_label(btn, base && *base ? base : "(None)");
			g_free(base);
			if (spec->on_picked)
				spec->on_picked(GTK_WIDGET(btn), path, spec->user_data);
			g_free(path);
		}
	}
	siril_file_browser_destroy(fb);
}

void siril_image_button_init(GtkWidget               *button,
                             const gchar             *title,
                             const gchar             *filter_title,
                             const gchar             *filter_pattern,
                             SirilImageButtonCallback on_picked,
                             gpointer                 user_data) {
	if (!button || !GTK_IS_BUTTON(button)) return;
	if (g_object_get_data(G_OBJECT(button), "siril-image-spec")) return;
	SirilImageButtonSpec *spec = g_new0(SirilImageButtonSpec, 1);
	spec->title          = title          ? g_strdup(title)          : NULL;
	spec->filter_title   = filter_title   ? g_strdup(filter_title)   : NULL;
	spec->filter_pattern = filter_pattern ? g_strdup(filter_pattern) : NULL;
	spec->on_picked      = on_picked;
	spec->user_data      = user_data;
	g_object_set_data_full(G_OBJECT(button), "siril-image-spec", spec,
	                       siril_image_button_spec_free);
	g_signal_connect(button, "clicked",
	                 G_CALLBACK(on_image_button_clicked), spec);
}

void siril_file_browser_destroy(SirilFileBrowser *fb) {
	if (!fb) return;
	/* Drop our refs on model objects FIRST, while the widget tree still
	 * holds its own refs.  Then destroy the window — which finalizes the
	 * widget tree, which drops the widget-owned refs. */
	g_clear_object(&fb->selection);
	g_clear_object(&fb->filter_model);
	g_clear_object(&fb->file_filter);
	g_clear_object(&fb->dir_list);
	g_clear_object(&fb->recent_store);
	if (fb->window) {
		gtk_window_destroy(fb->window);
		fb->window = NULL;
	}
	g_clear_object(&fb->current_folder);
	g_clear_object(&fb->initial_file);
	if (fb->back_history)    g_ptr_array_unref(fb->back_history);
	if (fb->forward_history) g_ptr_array_unref(fb->forward_history);
	if (fb->filters)         g_ptr_array_unref(fb->filters);
	g_slist_free_full(fb->result_paths, g_free);
	g_free(fb->result_path);
	g_free(fb);
}
