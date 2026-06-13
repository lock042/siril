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
#include "core/initfile.h"
#include "io/avi_preview.h"
#include "gui-gtk4/file_browser.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"

#include <string.h>

/* From io/image_format_fits.c — produces a small RGB byte buffer for any
 * FITS file plus a textual description string.  Used by the default
 * preview handler so the browser can show FITS thumbnails. */
extern guchar *extract_thumbnail_from_fits(const char *filename, gchar **descr,
                                            int *width_out, int *height_out);
/* From gui-gtk4/dialog_preview.c — first-frame thumbnail for SER files,
 * with MTF stretch applied for higher bit depths.  Same return contract
 * as the FITS extractor. */
extern guchar *extract_thumbnail_from_ser(const char *filename, gchar **descr,
                                           int *width_out, int *height_out);

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
	GtkEntry               *path_entry;       /* editable; Enter navigates */
	GtkStack               *path_stack;       /* "breadcrumb" / "entry" pages */
	GtkWidget              *breadcrumb_box;   /* horizontal box of segment buttons */
	GtkWidget              *breadcrumb_scroller; /* GtkScrolledWindow wrapping breadcrumb_box */
	GtkWidget              *breadcrumb_left_btn;  /* scroll-left chevron */
	GtkWidget              *breadcrumb_right_btn; /* scroll-right chevron */
	GtkWidget              *edit_path_btn;     /* toggle between breadcrumb / entry */
	GtkWidget              *new_folder_btn;    /* create a subfolder; shown in directory-only mode */
	GtkPopover             *new_folder_popover;/* name-entry popover for new_folder_btn */
	GtkEntry               *new_folder_entry;  /* folder name input inside the popover */
	guint                   pending_breadcrumb_idle; /* g_idle_add id, 0 = none */
	GtkColumnView          *columnview;
	GtkWidget              *outer_paned;        /* sidebar | content divider */
	GtkWidget              *list_preview_paned; /* list | preview divider */
	GtkPicture             *preview;
	GtkLabel               *metadata_label;
	GtkDropDown            *filter_combo;
	GtkCheckButton         *debayer_check;     /* shown when show_debayer_toggle is TRUE */
	GtkWidget              *open_button;
	GtkWidget              *back_button;
	GtkWidget              *forward_button;
	GtkWidget              *up_button;
	GtkSearchBar           *search_bar;
	GtkSearchEntry         *search_entry;
	gchar                  *search_text;       /* lower-cased substring, NULL/"" disables */
	gboolean                show_hidden_files; /* honors org.gtk.Settings.FileChooser show-hidden; toggle with Ctrl+H */

	/* Owned model objects.  Kept in fields so callbacks (filter, set_file,
	 * selection-changed) can reach them; references are taken via
	 * g_object_ref so my fields stay valid even after the widget tree
	 * teardown drops its own refs. */
	GtkDirectoryList       *dir_list;        /* +1 ref */
	GListStore             *recent_store;    /* +1 ref; GFileInfo items, full paths stashed via g_object_set_data */
	GtkFilterListModel     *filter_model;    /* +1 ref */
	GtkSortListModel       *sort_model;      /* +1 ref */
	GtkSorter              *combined_sorter; /* +1 ref; multi-sorter wrapping dirs-first + column-view's sorter */
	GtkCustomFilter        *file_filter;     /* +1 ref */
	GtkSelectionModel      *selection;       /* +1 ref; GtkSingleSelection or GtkMultiSelection */
	gboolean                select_multiple;
	gboolean                directory_only;  /* folder picker: files greyed + unselectable */
	gboolean                in_recent_mode;

	/* Sibling demosaic toggle (Convert tab) — kept in sync with our local
	 * `debayer_check` so the user sees consistent state in both places. */
	GtkCheckButton         *external_demosaic_btn;
	gulong                  external_demosaic_handler;
	gulong                  debayer_check_handler;
	gboolean                show_debayer_toggle;

	GFile                  *current_folder;     /* owned */
	GFile                  *initial_file;       /* owned, set before run */
	GFile                  *breadcrumb_deepest; /* owned; deepest native path
	                                             * visited along the current
	                                             * trail.  Lets the bar keep
	                                             * showing /a/b/c/d after the
	                                             * user goes up to /a/b/c so
	                                             * they can click back into d
	                                             * (or use the right arrow). */

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
	gboolean is_dir = g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY;
	/* Hide hidden entries (dotfiles, dotdirs) when show_hidden_files is off.
	 * Honors fb->show_hidden_files (initialised from GSettings
	 * org.gtk.Settings.FileChooser show-hidden, toggled with Ctrl+H). */
	if (g_file_info_get_is_hidden(info) && !fb->show_hidden_files)
		return FALSE;
	/* Search substring (case-insensitive).  Applies to both dirs and files
	 * when active — when the user is searching, an unmatched dir would be
	 * just noise, so we filter it too.  When search is inactive we keep the
	 * historical "dirs always visible" behaviour. */
	if (fb->search_text && *fb->search_text) {
		const char *display = g_file_info_get_display_name(info);
		if (!display) return FALSE;
		gchar *lower = g_utf8_strdown(display, -1);
		gboolean match = lower && strstr(lower, fb->search_text) != NULL;
		g_free(lower);
		if (!match) return FALSE;
	} else if (is_dir) {
		return TRUE;
	}
	if (is_dir) return TRUE;
	if (fb->active_filter_idx < 0 || fb->active_filter_idx >= (gint)fb->filters->len)
		return TRUE;
	BrowserFilter *bf = g_ptr_array_index(fb->filters, fb->active_filter_idx);
	return filter_accepts(bf, g_file_info_get_name(info));
}

/* ── sort comparators ─────────────────────────────────────────────── */

/* First-pass sorter that keeps directories pinned to the top irrespective
 * of which column is the active sort key.  Returns 0 within the dir / file
 * groups so the column-view's own sorter (multi-sorter fallthrough) takes
 * over to order siblings. */
static gint dirs_first_compare(gconstpointer a, gconstpointer b, gpointer ud) {
	(void)ud;
	GFileInfo *ia = (GFileInfo *) a;
	GFileInfo *ib = (GFileInfo *) b;
	gboolean da = g_file_info_get_file_type(ia) == G_FILE_TYPE_DIRECTORY;
	gboolean db = g_file_info_get_file_type(ib) == G_FILE_TYPE_DIRECTORY;
	if (da != db) return da ? -1 : 1;
	return 0;
}

static gint name_column_compare(gconstpointer a, gconstpointer b, gpointer ud) {
	(void)ud;
	const gchar *na = g_file_info_get_display_name((GFileInfo *)a);
	const gchar *nb = g_file_info_get_display_name((GFileInfo *)b);
	return g_utf8_collate(na ? na : "", nb ? nb : "");
}

static gint modified_column_compare(gconstpointer a, gconstpointer b, gpointer ud) {
	(void)ud;
	GDateTime *ta = g_file_info_get_modification_date_time((GFileInfo *)a);
	GDateTime *tb = g_file_info_get_modification_date_time((GFileInfo *)b);
	gint cmp;
	if (!ta && !tb)      cmp = 0;
	else if (!ta)        cmp = -1;
	else if (!tb)        cmp = 1;
	else                 cmp = g_date_time_compare(ta, tb);
	if (ta) g_date_time_unref(ta);
	if (tb) g_date_time_unref(tb);
	return cmp;
}

/* ── navigation / history ──────────────────────────────────────────── */

static int sidebar_recent_cmp_visited_desc(gconstpointer a, gconstpointer b);
static void update_nav_sensitivity(SirilFileBrowser *fb);
static void rebuild_breadcrumb(SirilFileBrowser *fb);
static void return_to_breadcrumb_view(SirilFileBrowser *fb);
static void update_breadcrumb_arrow_sensitivity(SirilFileBrowser *fb);

/* Initial position of the list | preview divider.  The divider is a soft
 * one now (the user may drag it); the preview's min size_request and the
 * paned's shrink-*-child=FALSE flags keep the preview from ever being
 * clipped — see the inner-paned setup in build_browser_widgets. */
#define LIST_PREVIEW_POSITION 480

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
		/* Query the real file so the GFileInfo carries the full standard::*
		 * attribute set (size, icon, type, display-name, mtime, is-hidden).
		 * Hand-building a partial GFileInfo means every column bind callback
		 * that reads an unset attribute raises a GLib-GIO CRITICAL (icon,
		 * size, …) — querying once gets them all and matches exactly what
		 * GtkDirectoryList produces for the normal listing. */
		GFile *gf = g_file_new_for_path(path);
		GFileInfo *fi = g_file_query_info(gf,
			"standard::*,time::modified",
			G_FILE_QUERY_INFO_NONE, NULL, NULL);
		g_object_unref(gf);
		if (!fi) {  /* vanished or unreadable between the test and the query */
			g_free(path);
			continue;
		}
		g_object_set_data_full(G_OBJECT(fi), "siril-full-path",
		                       g_strdup(path), g_free);
		g_list_store_append(fb->recent_store, fi);
		g_object_unref(fi);
		g_free(path);
		count++;
	}
	g_list_free_full(items, (GDestroyNotify) gtk_recent_info_unref);
}

static void exit_recent_mode(SirilFileBrowser *fb) {
	if (!fb->in_recent_mode) return;
	fb->in_recent_mode = FALSE;
	/* Same stale-focus precaution as apply_current_folder — the filter
	 * model swap fires items-changed which causes the column view to
	 * dispose its current rows. */
	if (fb->window) gtk_window_set_focus(fb->window, NULL);
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
	if (fb->window) gtk_window_set_focus(fb->window, NULL);
	gtk_filter_list_model_set_model(fb->filter_model,
		G_LIST_MODEL(fb->recent_store));
	gtk_editable_set_text(GTK_EDITABLE(fb->path_entry), _("Recent files"));
	rebuild_breadcrumb(fb);
	update_nav_sensitivity(fb);
}

static void update_path_label(SirilFileBrowser *fb) {
	if (!fb->current_folder) return;
	/* Non-native locations (trash:/// and friends) have no filesystem path;
	 * fall back to the URI so the bar shows *something* meaningful and the
	 * user knows where they are. */
	gchar *display = g_file_get_path(fb->current_folder);
	if (!display) display = g_file_get_uri(fb->current_folder);
	gtk_editable_set_text(GTK_EDITABLE(fb->path_entry), display ? display : "");
	/* Park the cursor at the end so the deepest part is visible (the
	 * GtkEntry scrolls to keep the cursor in view). */
	if (display)
		gtk_editable_set_position(GTK_EDITABLE(fb->path_entry), -1);
	g_free(display);

	/* Mirror the same destination into the breadcrumb so both views stay
	 * in sync — clicking a sidebar entry or back/forward updates both. */
	rebuild_breadcrumb(fb);
	/* Refresh the trail-navigation arrow sensitivity now that
	 * current_folder is up-to-date. */
	update_breadcrumb_arrow_sensitivity(fb);
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
 * external callers go through navigate_to(), which updates history.
 *
 * Live updates: GtkDirectoryList monitors `folder` itself once
 * gtk_directory_list_set_monitored is TRUE (the default; we set it
 * explicitly at construction).  When a file is added / modified /
 * deleted inside the current folder, the model emits items-changed,
 * the filter/sort/selection models pass that through to the view, and
 * the listing updates in place — selection of unaffected items is
 * preserved.  We rely on that built-in monitor rather than spinning
 * up a parallel GFileMonitor that would force a full re-enumeration. */
/* Maintain `fb->breadcrumb_deepest` so the breadcrumb shows a persistent
 * trail of the deepest folder visited along the current chain.  Going UP
 * or clicking an ancestor segment keeps the deeper levels visible (the
 * user can click forward into them); only navigating into a SIBLING or
 * UNRELATED branch truncates the trail to the new location.
 *
 *   ancestor:    keep existing deepest, just update which segment is current
 *   descendant:  extend deepest to the new (deeper) location
 *   same:        keep deepest as-is
 *   unrelated:   reset deepest to the new location
 *
 * Non-native paths (trash:///, network:///) can't participate in a chain
 * comparison, so we treat them as a reset to NULL — the breadcrumb will
 * render the URI scheme as a single segment instead. */
static void set_breadcrumb_deepest(SirilFileBrowser *fb, GFile *newpath) {
	if (!newpath || !g_file_is_native(newpath)) {
		g_clear_object(&fb->breadcrumb_deepest);
		return;
	}
	/* If newpath is *above* $HOME (e.g. user went up past /home/user to
	 * /home or /), reset the trail to newpath.  The chain walk in
	 * populate_breadcrumb_into stops at $HOME for display, so a stale
	 * trail below $HOME wouldn't contain newpath and the "current"
	 * highlight would land on no segment. */
	const gchar *home = g_get_home_dir();
	if (home) {
		GFile *home_gf = g_file_new_for_path(home);
		gboolean above_home = !g_file_equal(newpath, home_gf) &&
		                      g_file_has_prefix(home_gf, newpath);
		g_object_unref(home_gf);
		if (above_home) {
			g_clear_object(&fb->breadcrumb_deepest);
			fb->breadcrumb_deepest = g_object_ref(newpath);
			return;
		}
	}
	if (!fb->breadcrumb_deepest) {
		fb->breadcrumb_deepest = g_object_ref(newpath);
		return;
	}
	if (g_file_equal(newpath, fb->breadcrumb_deepest) ||
	    g_file_has_prefix(fb->breadcrumb_deepest, newpath)) {
		/* newpath is ancestor-of (or equal-to) deepest → keep deepest. */
		return;
	}
	if (g_file_has_prefix(newpath, fb->breadcrumb_deepest)) {
		/* newpath is strictly deeper than deepest → extend to newpath. */
		g_object_unref(fb->breadcrumb_deepest);
		fb->breadcrumb_deepest = g_object_ref(newpath);
		return;
	}
	/* Unrelated branch — start a fresh trail at newpath. */
	g_object_unref(fb->breadcrumb_deepest);
	fb->breadcrumb_deepest = g_object_ref(newpath);
}

static void apply_current_folder(SirilFileBrowser *fb, GFile *folder) {
	if (!folder) return;
	exit_recent_mode(fb);  /* any folder navigation cancels Recent view */
	if (fb->current_folder) g_object_unref(fb->current_folder);
	fb->current_folder = g_object_ref(folder);
	set_breadcrumb_deepest(fb, folder);
	/* Clear window focus BEFORE swapping the directory list's file.
	 * gtk_directory_list_set_file synchronously fires items-changed to
	 * remove every current item, and the column view recycles / disposes
	 * its row widgets in response.  If a row still holds focus when its
	 * widget is disposed, GTK4's recursive focus / state propagation
	 * walker hits a stale widget pointer and trips
	 *   `gtk_widget_is_ancestor: assertion 'GTK_IS_WIDGET (widget)' failed`.
	 * Parking focus on the window itself before the model swap means no
	 * row widget is in the focus chain when it gets disposed. */
	if (fb->window)
		gtk_window_set_focus(fb->window, NULL);
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

/* Up / Back / Forward button handlers retired — the breadcrumb's
 * left/right arrows cover step-up and step-down-along-trail, and the
 * trail itself replaces history.  The fb->back_history /
 * forward_history arrays are still populated by navigate_to() for now
 * (cheap insurance in case we want Alt+Left / Alt+Right shortcuts
 * later) but they no longer drive any UI. */

static void on_path_entry_activate(GtkEntry *entry, gpointer ud) {
	SirilFileBrowser *fb = ud;
	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
	if (!text || !*text) return;

	/* URI form (e.g. trash:///, file:///, sftp://host/path) — navigate
	 * straight through GIO without touching the filesystem.  This lets the
	 * user re-trigger navigation to non-native locations from the address
	 * bar after we routed them there via the sidebar. */
	if (strstr(text, "://")) {
		gtk_widget_remove_css_class(GTK_WIDGET(entry), "error");
		GFile *gf = g_file_new_for_uri(text);
		navigate_to(fb, gf);
		g_object_unref(gf);
		return_to_breadcrumb_view(fb);
		return;
	}

	/* Accept ~ for $HOME — convenient when typing by hand. */
	gchar *expanded = NULL;
	if (text[0] == '~' && (text[1] == '\0' || text[1] == G_DIR_SEPARATOR)) {
		const gchar *home = g_get_home_dir();
		if (home)
			expanded = g_strconcat(home, text + 1, NULL);
	}
	const gchar *path = expanded ? expanded : text;

	if (!g_file_test(path, G_FILE_TEST_EXISTS)) {
		/* Don't disrupt the user's typing — flag the entry with the
		 * "error" CSS class so it's obvious nothing happened.  Cleared
		 * the next time update_path_label runs. */
		gtk_widget_add_css_class(GTK_WIDGET(entry), "error");
		g_free(expanded);
		return;
	}
	gtk_widget_remove_css_class(GTK_WIDGET(entry), "error");

	GFile *gf = g_file_new_for_path(path);
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		navigate_to(fb, gf);
	} else {
		/* Plain file: navigate to its parent and stop (selection is
		 * cheap to do by eye in a small directory). */
		GFile *parent = g_file_get_parent(gf);
		if (parent) {
			navigate_to(fb, parent);
			g_object_unref(parent);
		}
	}
	g_object_unref(gf);
	g_free(expanded);
	return_to_breadcrumb_view(fb);
}

/* ── breadcrumb path bar ──────────────────────────────────────────── */

static void on_breadcrumb_clicked(GtkButton *btn, gpointer ud) {
	SirilFileBrowser *fb = ud;
	GFile *gf = g_object_get_data(G_OBJECT(btn), "siril-gfile");
	if (!gf) return;
	if (fb->current_folder && g_file_equal(gf, fb->current_folder))
		return;
	navigate_to(fb, gf);
}

/* Build one segment button.  `label`, `icon_name`, or both may be NULL/empty
 * — at least one should be non-NULL.  The button holds an owned ref to `gf`
 * so it can be reused after the GFile chain is freed. */
static GtkWidget *make_breadcrumb_button(SirilFileBrowser *fb,
                                         const char *label,
                                         const char *icon_name,
                                         GFile *gf) {
	GtkWidget *btn = gtk_button_new();
	/* No `.flat` — keep the default button chrome so segments read as
	 * clickable buttons rather than mere labels.  The fine-grained look
	 * (subtle background, hover, current-segment bold) is in siril.css
	 * under .siril-breadcrumb-segment. */
	gtk_widget_add_css_class(btn, "siril-breadcrumb-segment");
	/* Keep these out of the focus chain.  If a focusable segment is
	 * removed during rebuild_breadcrumb (which fires whenever the user
	 * navigates), GTK4 walks the now-disposed widget's ancestor chain
	 * to update focus and trips a `gtk_widget_is_ancestor` assertion.
	 * Breadcrumbs are mouse-driven anyway — keyboard users have the
	 * entry mode via the edit-path toggle. */
	gtk_widget_set_focusable(btn, FALSE);
	gtk_widget_set_can_focus(btn, FALSE);
	if (icon_name && label && *label) {
		GtkWidget *box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
		gtk_box_append(GTK_BOX(box), gtk_image_new_from_icon_name(icon_name));
		gtk_box_append(GTK_BOX(box), gtk_label_new(label));
		gtk_button_set_child(GTK_BUTTON(btn), box);
	} else if (icon_name) {
		gtk_button_set_icon_name(GTK_BUTTON(btn), icon_name);
	} else {
		gtk_button_set_label(GTK_BUTTON(btn), label ? label : "");
	}
	g_object_set_data_full(G_OBJECT(btn), "siril-gfile",
	                       g_object_ref(gf), g_object_unref);
	g_signal_connect(btn, "clicked",
	                 G_CALLBACK(on_breadcrumb_clicked), fb);
	return btn;
}

/* Walk the trail from `breadcrumb_deepest` up until we find the GFile
 * whose direct parent is `current_folder`.  That's the "next deeper" step
 * along the trail — clicking the right arrow navigates to it.  Returns
 * a transfer-full GFile, or NULL if the trail doesn't extend below
 * current. */
static GFile *trail_next_below_current(SirilFileBrowser *fb) {
	if (!fb->breadcrumb_deepest || !fb->current_folder) return NULL;
	if (g_file_equal(fb->current_folder, fb->breadcrumb_deepest)) return NULL;
	if (!g_file_has_prefix(fb->breadcrumb_deepest, fb->current_folder))
		return NULL;
	GFile *p = g_object_ref(fb->breadcrumb_deepest);
	while (p) {
		GFile *parent = g_file_get_parent(p);
		if (!parent) {
			g_object_unref(p);
			return NULL;
		}
		if (g_file_equal(parent, fb->current_folder)) {
			g_object_unref(parent);
			return p;
		}
		g_object_unref(p);
		p = parent;
	}
	return NULL;
}

/* Sensitivity of the trail-navigation arrows:
 *   left  (step-up):   enabled if current_folder has a parent
 *   right (step-down): enabled if the trail extends below current_folder */
static void update_breadcrumb_arrow_sensitivity(SirilFileBrowser *fb) {
	gboolean can_left = FALSE, can_right = FALSE;
	if (fb->current_folder) {
		GFile *p = g_file_get_parent(fb->current_folder);
		if (p) {
			can_left = TRUE;
			g_object_unref(p);
		}
	}
	if (fb->breadcrumb_deepest && fb->current_folder &&
	    !g_file_equal(fb->current_folder, fb->breadcrumb_deepest) &&
	    g_file_has_prefix(fb->breadcrumb_deepest, fb->current_folder)) {
		can_right = TRUE;
	}
	if (fb->breadcrumb_left_btn)
		gtk_widget_set_sensitive(fb->breadcrumb_left_btn, can_left);
	if (fb->breadcrumb_right_btn)
		gtk_widget_set_sensitive(fb->breadcrumb_right_btn, can_right);
}

static void on_breadcrumb_step_up(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (!fb->current_folder) return;
	GFile *parent = g_file_get_parent(fb->current_folder);
	if (parent) {
		navigate_to(fb, parent);
		g_object_unref(parent);
	}
}

static void on_breadcrumb_step_down(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	GFile *next = trail_next_below_current(fb);
	if (next) {
		navigate_to(fb, next);
		g_object_unref(next);
	}
}

/* Populate `box` with segment buttons for the given folder.  No chevron
 * separators — segments butt up against one another and rely on
 * background contrast / hover effects for visual separation.  The
 * topmost segment (current folder) gets the `.suggested-action` class
 * so it reads as the "selected" / active cell in the bar.  Used by
 * rebuild_breadcrumb, which builds a fresh box rather than mutating
 * the existing one. */
static void populate_breadcrumb_into(SirilFileBrowser *fb, GtkWidget *box) {
	/* Recent mode: single dimmed label, no clickable segments. */
	if (fb->in_recent_mode) {
		GtkWidget *lbl = gtk_label_new(_("Recent files"));
		gtk_widget_add_css_class(lbl, "dim-label");
		gtk_widget_set_margin_start(lbl, 8);
		gtk_box_append(GTK_BOX(box), lbl);
		return;
	}

	if (!fb->current_folder) return;

	/* Non-native location (trash://, network://, sftp://…): we can't walk
	 * a parent chain meaningfully, so render a single segment named after
	 * the URI scheme.  Trash gets the system trash icon for recognisability. */
	if (!g_file_is_native(fb->current_folder)) {
		gchar *uri = g_file_get_uri(fb->current_folder);
		const char *scheme = uri ? g_uri_peek_scheme(uri) : NULL;
		const char *label  = scheme ? scheme : (uri ? uri : "?");
		const char *icon   = "folder-symbolic";
		const char *display = label;
		if (scheme && g_str_equal(scheme, "trash")) {
			icon = "user-trash-symbolic";
			display = _("Trash");
		}
		GtkWidget *seg = make_breadcrumb_button(fb, display, icon,
		                                        fb->current_folder);
		gtk_box_append(GTK_BOX(box), seg);
		g_free(uri);
		return;
	}

	/* Native filesystem: walk parents up to root (or up to $HOME) from
	 * the *deepest* path we've visited along this trail, not from
	 * fb->current_folder — so going up a level keeps the deeper segments
	 * visible / clickable.  fb->breadcrumb_deepest holds the deepest
	 * point; if it's NULL (first navigation, or just reset) the trail
	 * starts at the current folder. */
	GFile *trail_tip = fb->breadcrumb_deepest ? fb->breadcrumb_deepest
	                                          : fb->current_folder;
	const gchar *home = g_get_home_dir();
	GFile *home_gfile = home ? g_file_new_for_path(home) : NULL;
	GPtrArray *chain = g_ptr_array_new_with_free_func(g_object_unref);
	g_ptr_array_add(chain, g_object_ref(trail_tip));
	gboolean stopped_at_home = home_gfile &&
	                           g_file_equal(trail_tip, home_gfile);
	if (!stopped_at_home) {
		GFile *cur = trail_tip;
		GFile *p;
		while ((p = g_file_get_parent(cur))) {
			g_ptr_array_add(chain, p);    /* owns */
			if (home_gfile && g_file_equal(p, home_gfile)) {
				stopped_at_home = TRUE;
				break;
			}
			cur = p;
		}
	}

	/* chain[0] = trail tip (deepest), chain[len-1] = top (root or $HOME).
	 * The "current" segment is whichever GFile in the chain equals
	 * fb->current_folder — that's the one we render bold. */
	for (gint i = chain->len - 1; i >= 0; i--) {
		GFile *seg = g_ptr_array_index(chain, i);
		gboolean is_top = (i == (gint)chain->len - 1);
		gboolean is_current = g_file_equal(seg, fb->current_folder);
		const char *segname = NULL;
		const char *segicon = NULL;
		gchar *base = NULL;
		if (is_top && stopped_at_home) {
			/* Home: house icon + username label (e.g. "🏠 alice"). */
			segicon = "user-home-symbolic";
			segname = g_get_user_name();
		} else if (is_top) {
			segname = "/";
			segicon = "drive-harddisk-symbolic";
		} else {
			base = g_file_get_basename(seg);
			segname = base ? base : "?";
		}
		GtkWidget *btn = make_breadcrumb_button(fb, segname, segicon, seg);
		/* Highlight the segment that matches the current folder (which
		 * may be anywhere along the trail, not just the rightmost). */
		if (is_current)
			gtk_widget_add_css_class(btn, "siril-current");
		gtk_box_append(GTK_BOX(box), btn);
		g_free(base);
	}
	g_ptr_array_unref(chain);
	g_clear_object(&home_gfile);
}

/* Build a brand-new GtkBox and swap it into the scrolled window
 * atomically rather than incrementally removing children from the
 * existing box.  Incremental removal made GTK4 walk a partially-
 * disposed focus chain during each gtk_box_remove and trip a
 *   `gtk_widget_is_ancestor: assertion 'GTK_IS_WIDGET (widget)' failed`
 * critical. */
static void rebuild_breadcrumb_now(SirilFileBrowser *fb) {
	if (!fb->breadcrumb_scroller) return;

	GtkWidget *new_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_add_css_class(new_box, "siril-breadcrumb");
	populate_breadcrumb_into(fb, new_box);

	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(fb->breadcrumb_scroller),
	                              new_box);
	fb->breadcrumb_box = new_box;

	/* Sensitivity of the scroll-arrow buttons follows the hadjustment's
	 * notify signals once layout has settled.  The notify::upper handler
	 * we wired at construction time fires after the new box is allocated,
	 * so the arrows update themselves without us having to chase the
	 * upper bound from here. */
}

static gboolean rebuild_breadcrumb_idle(gpointer ud) {
	SirilFileBrowser *fb = ud;
	fb->pending_breadcrumb_idle = 0;
	rebuild_breadcrumb_now(fb);
	return G_SOURCE_REMOVE;
}

/* Rebuild the breadcrumb from fb->current_folder.  Called whenever the
 * folder changes (via update_path_label) and on Recent-mode entry.
 *
 * Deferred to an idle so the rebuild doesn't run inside the column-view's
 * "activate" signal dispatch — incremental focus / accessibility updates
 * during that dispatch were triggering `gtk_widget_is_ancestor: assertion
 * 'GTK_IS_WIDGET (widget)' failed` even with the atomic box swap.
 * Running on idle lets GTK drain its existing dispatch chain first, so
 * the rebuild operates on a settled widget tree.
 *
 * Coalesces multiple rapid calls: if an idle is already pending, the
 * existing one will pick up the latest fb->current_folder when it fires.
 * The pending source id is cancelled in reset_browser_state so a dialog
 * reopened mid-navigation doesn't fire the idle against a wiped state. */
static void rebuild_breadcrumb(SirilFileBrowser *fb) {
	if (!fb) return;
	if (fb->pending_breadcrumb_idle) return;
	fb->pending_breadcrumb_idle = g_idle_add(rebuild_breadcrumb_idle, fb);
}

static void on_edit_path_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	SirilFileBrowser *fb = ud;
	if (!fb->path_stack) return;
	/* Proper toggle: clicking the button when already in entry mode
	 * collapses back to breadcrumb mode (and discards any unsubmitted
	 * typing in the entry, matching Escape behaviour).  Previously
	 * clicking again was a no-op so the user had to press Enter on
	 * something valid or hit Escape to get back. */
	const char *cur = gtk_stack_get_visible_child_name(fb->path_stack);
	if (g_strcmp0(cur, "entry") == 0) {
		return_to_breadcrumb_view(fb);
	} else {
		gtk_stack_set_visible_child_name(fb->path_stack, "entry");
		gtk_widget_grab_focus(GTK_WIDGET(fb->path_entry));
		gtk_editable_select_region(GTK_EDITABLE(fb->path_entry), 0, -1);
	}
}

static void return_to_breadcrumb_view(SirilFileBrowser *fb) {
	if (!fb->path_stack) return;
	/* Clear the error indicator so a stale red border doesn't carry into
	 * the next edit session. */
	gtk_widget_remove_css_class(GTK_WIDGET(fb->path_entry), "error");
	gtk_stack_set_visible_child_name(fb->path_stack, "breadcrumb");
}

static gboolean on_path_entry_key_pressed(GtkEventControllerKey *c,
                                          guint keyval, guint keycode,
                                          GdkModifierType state,
                                          gpointer ud) {
	(void)c; (void)keycode; (void)state;
	SirilFileBrowser *fb = ud;
	if (keyval == GDK_KEY_Escape) {
		return_to_breadcrumb_view(fb);
		return TRUE;
	}
	return FALSE;
}

static void on_new_folder_create(SirilFileBrowser *fb) {
	if (!fb || !fb->current_folder) return;
	const gchar *raw = gtk_editable_get_text(GTK_EDITABLE(fb->new_folder_entry));
	gchar *name = raw ? g_strstrip(g_strdup(raw)) : NULL;
	if (!name || !*name) {
		g_free(name);
		return;  /* nothing typed yet — keep the popover open */
	}
	/* Reject path separators: this creates a single child, not a tree. */
	if (strchr(name, G_DIR_SEPARATOR) ||
	    !g_strcmp0(name, ".") || !g_strcmp0(name, "..")) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid folder name"),
			_("The folder name must not be empty or contain a path separator."));
		g_free(name);
		return;
	}

	GFile *child = g_file_get_child(fb->current_folder, name);
	GError *err = NULL;
	gboolean ok = g_file_make_directory(child, NULL, &err);
	if (ok) {
		gtk_editable_set_text(GTK_EDITABLE(fb->new_folder_entry), "");
		if (fb->new_folder_popover)
			gtk_popover_popdown(fb->new_folder_popover);
		navigate_to(fb, child);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not create folder"),
			err && err->message ? err->message : _("Unknown error"));
	}
	g_clear_error(&err);
	g_object_unref(child);
	g_free(name);
}

static void on_new_folder_entry_activate(GtkEntry *entry, gpointer ud) {
	(void)entry;
	on_new_folder_create(ud);
}

static void on_new_folder_confirm_clicked(GtkButton *b, gpointer ud) {
	(void)b;
	on_new_folder_create(ud);
}

/* Clear any stale text and focus the entry each time the popover opens. */
static void on_new_folder_popover_show(GtkPopover *p, gpointer ud) {
	(void)p;
	SirilFileBrowser *fb = ud;
	gtk_editable_set_text(GTK_EDITABLE(fb->new_folder_entry), "");
	gtk_widget_grab_focus(GTK_WIDGET(fb->new_folder_entry));
}

/* ── column factories ─────────────────────────────────────────────── */

/* Forward decl — defined further down (after the picture / preview helpers)
 * but needed by the Type column factory. */
static const char *type_short_label(image_type t, const char *path);

/* Helper: zero vertical margins on a cell child so its row hugs the text
 * baseline.  Row density is enforced by the `siril-dense-rows` CSS class
 * on the columnview; this just ensures the factory widget doesn't itself
 * contribute extra padding. */
static void tighten_cell_widget(GtkWidget *w) {
	gtk_widget_set_margin_top(w, 0);
	gtk_widget_set_margin_bottom(w, 0);
}

/* Name column: icon + label. */
/* In directory-only mode, files (non-folders) are shown but greyed out
 * (insensitive cell) and made non-selectable / non-activatable so only
 * folders can be picked.  Called from every column's bind callback so the
 * whole row reflects the state consistently; ud is the SirilFileBrowser
 * (passed as the bind handlers' user_data). */
static void apply_dir_only_row(gpointer ud, GtkListItem *item,
                               GFileInfo *info, GtkWidget *child) {
	SirilFileBrowser *fb = ud;
	gboolean is_dir = g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY;
	gboolean blocked = fb && fb->directory_only && !is_dir;
	if (child)
		gtk_widget_set_sensitive(child, !blocked);
	gtk_list_item_set_selectable(item, !blocked);
	gtk_list_item_set_activatable(item, !blocked);
}

static void name_setup_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	tighten_cell_widget(box);
	GtkWidget *icon = gtk_image_new();
	GtkWidget *label = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_widget_set_hexpand(label, TRUE);
	gtk_box_append(GTK_BOX(box), icon);
	gtk_box_append(GTK_BOX(box), label);
	gtk_list_item_set_child(item, box);
}

static void name_bind_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f;
	GFileInfo *info = G_FILE_INFO(gtk_list_item_get_item(item));
	GtkWidget *box = gtk_list_item_get_child(item);
	GtkWidget *icon = gtk_widget_get_first_child(box);
	GtkWidget *label = gtk_widget_get_next_sibling(icon);
	apply_dir_only_row(ud, item, info, box);

	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY) {
		gtk_image_set_from_icon_name(GTK_IMAGE(icon), "folder-symbolic");
	} else {
		/* Use the file's mime-type-derived GIcon if available so different
		 * file kinds get distinguishable icons (image-jpeg, text-plain,
		 * application-pdf, etc.).  Falls back to the generic glyph if GIO
		 * couldn't resolve one.  Guard with has_attribute: GFileInfos we
		 * build by hand (the Recent store) never set standard::icon, and
		 * g_file_info_get_icon() raises a GLib-GIO CRITICAL when the
		 * attribute was never set rather than just returning NULL. */
		GIcon *gicon = g_file_info_has_attribute(info, G_FILE_ATTRIBUTE_STANDARD_ICON)
			? g_file_info_get_icon(info) : NULL;
		if (gicon)
			gtk_image_set_from_gicon(GTK_IMAGE(icon), gicon);
		else
			gtk_image_set_from_icon_name(GTK_IMAGE(icon), "text-x-generic-symbolic");
	}
	gtk_label_set_text(GTK_LABEL(label), g_file_info_get_display_name(info));
}

/* Modified column: a single right-aligned label with the mtime. */
static void modified_setup_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GtkWidget *label = gtk_label_new(NULL);
	tighten_cell_widget(label);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_widget_add_css_class(label, "dim-label");
	gtk_list_item_set_child(item, label);
}

static void modified_bind_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f;
	GFileInfo *info = G_FILE_INFO(gtk_list_item_get_item(item));
	GtkLabel *label = GTK_LABEL(gtk_list_item_get_child(item));
	apply_dir_only_row(ud, item, info, GTK_WIDGET(label));
	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY) {
		/* Folders: empty mtime column to keep the eye on the name. */
		gtk_label_set_text(label, "");
		return;
	}
	GDateTime *mt = g_file_info_get_modification_date_time(info);
	if (!mt) {
		gtk_label_set_text(label, "");
		return;
	}
	gchar *s = g_date_time_format(mt, "%Y-%m-%d %H:%M");
	gtk_label_set_text(label, s ? s : "");
	g_free(s);
	g_date_time_unref(mt);
}

/* Size column: right-aligned human-readable byte count. */
static void size_setup_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GtkWidget *label = gtk_label_new(NULL);
	tighten_cell_widget(label);
	gtk_label_set_xalign(GTK_LABEL(label), 1.0);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_widget_add_css_class(label, "dim-label");
	gtk_list_item_set_child(item, label);
}

static void size_bind_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f;
	GFileInfo *info = G_FILE_INFO(gtk_list_item_get_item(item));
	GtkLabel *label = GTK_LABEL(gtk_list_item_get_child(item));
	apply_dir_only_row(ud, item, info, GTK_WIDGET(label));
	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY) {
		/* Don't recurse to total folder contents — too expensive in a
		 * file picker.  An em-dash makes the empty column self-evidently
		 * intentional rather than missing data. */
		gtk_label_set_text(label, "—");
		return;
	}
	goffset bytes = g_file_info_get_size(info);
	gchar *s = g_format_size(bytes);
	gtk_label_set_text(label, s ? s : "");
	g_free(s);
}

/* Type column: short uppercase format label (FITS, JPEG, ...) or "Folder". */
static void type_setup_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f; (void)ud;
	GtkWidget *label = gtk_label_new(NULL);
	tighten_cell_widget(label);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_widget_add_css_class(label, "dim-label");
	gtk_list_item_set_child(item, label);
}

static void type_bind_cb(GtkSignalListItemFactory *f, GtkListItem *item, gpointer ud) {
	(void)f;
	GFileInfo *info = G_FILE_INFO(gtk_list_item_get_item(item));
	GtkLabel *label = GTK_LABEL(gtk_list_item_get_child(item));
	apply_dir_only_row(ud, item, info, GTK_WIDGET(label));
	if (g_file_info_get_file_type(info) == G_FILE_TYPE_DIRECTORY) {
		gtk_label_set_text(label, _("Folder"));
		return;
	}
	/* get_type_from_filename only inspects the extension, so passing
	 * just the basename is enough — no need to resolve a full path. */
	const char *name = g_file_info_get_name(info);
	image_type t = name ? get_type_from_filename(name) : TYPEUNDEF;
	gtk_label_set_text(label, type_short_label(t, name));
}

/* ── column sort comparators ───────────────────────────────────────── */

static gint size_column_compare(gconstpointer a, gconstpointer b, gpointer ud) {
	(void)ud;
	GFileInfo *ia = (GFileInfo *)a;
	GFileInfo *ib = (GFileInfo *)b;
	gboolean da = g_file_info_get_file_type(ia) == G_FILE_TYPE_DIRECTORY;
	gboolean db = g_file_info_get_file_type(ib) == G_FILE_TYPE_DIRECTORY;
	/* Folders carry no meaningful size; the dirs-first primary sorter
	 * already groups them, but compare equal among themselves so the
	 * column-view's secondary key (name) can break ties cleanly. */
	if (da && db) return 0;
	goffset sa = g_file_info_get_size(ia);
	goffset sb = g_file_info_get_size(ib);
	return (sa > sb) - (sa < sb);
}

static gint type_column_compare(gconstpointer a, gconstpointer b, gpointer ud) {
	(void)ud;
	GFileInfo *ia = (GFileInfo *)a;
	GFileInfo *ib = (GFileInfo *)b;
	const char *na = g_file_info_get_name(ia);
	const char *nb = g_file_info_get_name(ib);
	image_type ta = (g_file_info_get_file_type(ia) == G_FILE_TYPE_DIRECTORY)
	                ? TYPEUNDEF : (na ? get_type_from_filename(na) : TYPEUNDEF);
	image_type tb = (g_file_info_get_file_type(ib) == G_FILE_TYPE_DIRECTORY)
	                ? TYPEUNDEF : (nb ? get_type_from_filename(nb) : TYPEUNDEF);
	const char *la = type_short_label(ta, na);
	const char *lb = type_short_label(tb, nb);
	return g_ascii_strcasecmp(la ? la : "", lb ? lb : "");
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

	/* For items inside non-native locations (notably trash://), the
	 * GIO backend stashes the real on-disk location in
	 * standard::target-uri.  The trashed file actually lives in
	 * $XDG_DATA_HOME/Trash/files/<name> (or the per-mount trash dir),
	 * and that path *is* native and openable.  Prefer the target-uri
	 * when present so Open / preview work just like on a normal file. */
	const char *target_uri = g_file_info_get_attribute_string(info,
		G_FILE_ATTRIBUTE_STANDARD_TARGET_URI);
	if (target_uri && *target_uri) {
		GFile *t = g_file_new_for_uri(target_uri);
		gchar *p = g_file_get_path(t);
		g_object_unref(t);
		if (p) return p;
	}

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
	/* The path is passed through for directories too — the default
	 * preview handler renders a folder icon for them, and custom
	 * callbacks can choose to do likewise (or short-circuit). */
	gchar *path = current_selected_path(fb);
	SirilFileBrowserPreview cb = fb->preview_cb ? fb->preview_cb
	                                            : siril_file_browser_default_preview;
	cb(path, fb->preview, fb->metadata_label,
	   fb->preview_cb ? fb->preview_cb_data : fb);
	g_free(path);
}

static void on_selection_changed(GtkSelectionModel *m, guint pos, guint n, gpointer ud) {
	(void)m; (void)pos; (void)n;
	SirilFileBrowser *fb = ud;
	/* Directory-only mode: a file should never end up selected.  The bind
	 * callbacks already mark file rows non-selectable, but guard here too
	 * (re-selecting NULL is harmless and re-enters with nothing selected). */
	if (fb->directory_only && !fb->select_multiple && fb->selection) {
		GFileInfo *info = first_selected_info(fb);
		if (info && g_file_info_get_file_type(info) != G_FILE_TYPE_DIRECTORY) {
			gtk_single_selection_set_selected(
				GTK_SINGLE_SELECTION(fb->selection), GTK_INVALID_LIST_POSITION);
			return;
		}
	}
	update_preview_for_selection(fb);
	/* In folder mode the Open button targets the current (or highlighted)
	 * directory, so keep it usable even with nothing highlighted. */
	gtk_widget_set_sensitive(fb->open_button,
	                         fb->directory_only || any_selected(fb));
}

static void on_row_activated(GtkColumnView *cv, guint position, gpointer ud) {
	(void)cv; (void)position;
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

/* ── search bar ───────────────────────────────────────────────────── */

static void on_search_changed(GtkSearchEntry *entry, gpointer ud) {
	SirilFileBrowser *fb = ud;
	const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
	g_clear_pointer(&fb->search_text, g_free);
	if (text && *text) {
		/* Stash lower-cased so the filter func can do a single
		 * strstr() per item instead of an allocation per item. */
		fb->search_text = g_utf8_strdown(text, -1);
	}
	gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);
}

/* Track the search bar's mode so closing it (Escape, or untoggling the
 * button) wipes the substring — otherwise reopening the search shows a
 * stale entry value but no longer filters, which is confusing. */
static void on_search_mode_changed(GObject *obj, GParamSpec *p, gpointer ud) {
	(void)p;
	SirilFileBrowser *fb = ud;
	if (gtk_search_bar_get_search_mode(GTK_SEARCH_BAR(obj))) {
		gtk_widget_grab_focus(GTK_WIDGET(fb->search_entry));
	} else {
		gtk_editable_set_text(GTK_EDITABLE(fb->search_entry), "");
		/* `set_text("")` already triggers on_search_changed → filter
		 * refresh, so no manual refresh needed here. */
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
	/* Folder picker: commit the highlighted directory, or fall back to the
	 * folder we're currently browsing when nothing is highlighted. */
	if (fb->directory_only) {
		gchar *p = NULL;
		if (any_selected(fb) && current_is_directory(fb))
			p = current_selected_path(fb);
		else if (fb->current_folder)
			p = g_file_get_path(fb->current_folder);
		if (!p) return;
		g_free(fb->result_path);
		fb->result_path = p;
		fb->response = GTK_RESPONSE_ACCEPT;
		if (fb->loop) g_main_loop_quit(fb->loop);
		return;
	}
	if (current_is_directory(fb)) {
		on_row_activated(fb->columnview, 0, fb);
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

/* Forward decl — used by the row-activated handler before its definition. */
static void sidebar_mount_volume_async(SirilFileBrowser *fb, GVolume *vol,
                                       GtkListBoxRow *row);

/* Sidebar row data keys, all attached via g_object_set_data{,_full}:
 *   "siril-is-recent"  GINT_TO_POINTER(1)     — Recent virtual entry
 *   "siril-path"       char* (owned)          — filesystem path
 *   "siril-gfile"      GFile* (owned)         — URI-backed entry (trash://, ...)
 *   "siril-volume"     GVolume* (owned)       — unmounted volume; activate mounts
 * Activation precedence: recent > volume > gfile > path. */
static void on_sidebar_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer ud) {
	(void)box;
	SirilFileBrowser *fb = ud;
	if (g_object_get_data(G_OBJECT(row), "siril-is-recent")) {
		enter_recent_mode(fb);
		return;
	}
	GVolume *vol = g_object_get_data(G_OBJECT(row), "siril-volume");
	if (vol) {
		sidebar_mount_volume_async(fb, vol, row);
		return;
	}
	GFile *gf = g_object_get_data(G_OBJECT(row), "siril-gfile");
	if (gf) {
		navigate_to(fb, gf);
		return;
	}
	const char *path = g_object_get_data(G_OBJECT(row), "siril-path");
	if (!path) return;
	GFile *gpf = g_file_new_for_path(path);
	navigate_to(fb, gpf);
	g_object_unref(gpf);
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
	gtk_widget_set_valign(icon, GTK_ALIGN_CENTER);
	GtkWidget *labels = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	/* Vertically centre the label stack against the icon — default
	 * valign=fill stretches the VBox to the row height, leaving the
	 * single label flush with the top of the row instead of in line
	 * with the icon. */
	gtk_widget_set_valign(labels, GTK_ALIGN_CENTER);
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

/* Attach a row that navigates to an arbitrary GFile (URI or path).  Used
 * for entries like trash:/// where g_file_test()-based filtering is
 * impossible and the activation handler must keep the GFile alive. */
static void sidebar_add_gfile(GtkListBox *box, const char *icon_name,
                              const char *label, GFile *gf) {
	if (!gf) return;
	GtkWidget *row = sidebar_make_row(icon_name, label, NULL);
	g_object_set_data_full(G_OBJECT(row), "siril-gfile",
	                       g_object_ref(gf), g_object_unref);
	gtk_list_box_append(box, row);
}

/* Pick the first themed icon name from a GIcon, or `fallback` if none. */
static const char *first_themed_icon_name(GIcon *gicon, const char *fallback) {
	if (G_IS_THEMED_ICON(gicon)) {
		const gchar * const *names = g_themed_icon_get_names(G_THEMED_ICON(gicon));
		if (names && names[0]) return names[0];
	}
	return fallback;
}

/* Add a row for an unmounted GVolume.  Activation triggers an async mount
 * (handled by sidebar_mount_volume_async) and on success navigates to the
 * resulting mount root.  Carries a `siril-unmounted` CSS class so the row
 * can be visually distinguished from already-mounted entries. */
static void sidebar_add_volume_unmounted(GtkListBox *box, GVolume *vol) {
	gchar *name = g_volume_get_name(vol);
	GIcon *gicon = g_volume_get_symbolic_icon(vol);
	if (!gicon) gicon = g_volume_get_icon(vol);
	const char *icon_name = first_themed_icon_name(gicon, "drive-removable-media-symbolic");
	GtkWidget *row = sidebar_make_row(icon_name, name ? name : _("Volume"),
	                                  _("Not mounted"));
	gtk_widget_add_css_class(row, "siril-unmounted");
	gtk_widget_set_tooltip_text(row, _("Click to mount this device"));
	g_object_set_data_full(G_OBJECT(row), "siril-volume",
	                       g_object_ref(vol), g_object_unref);
	gtk_list_box_append(box, row);
	if (gicon) g_object_unref(gicon);
	g_free(name);
}

/* Walk g_volume_monitor_get_volumes() and append a row for each volume
 * that is mountable but not currently mounted (mounted ones are already
 * shown by sidebar_populate_volumes via get_mounts).  Returns the number
 * of rows added.  Drives with no media are filtered out — a "Mount X"
 * row that can never succeed is just noise. */
static int sidebar_populate_unmounted_volumes(GtkListBox *box) {
	GVolumeMonitor *monitor = g_volume_monitor_get();
	if (!monitor) return 0;
	int count = 0;
	GList *vols = g_volume_monitor_get_volumes(monitor);
	for (GList *l = vols; l; l = l->next) {
		GVolume *vol = G_VOLUME(l->data);
		GMount *mount = g_volume_get_mount(vol);
		if (mount) {
			/* Already mounted — sidebar_populate_volumes covered it. */
			g_object_unref(mount);
			continue;
		}
		if (!g_volume_can_mount(vol)) continue;
		sidebar_add_volume_unmounted(box, vol);
		count++;
	}
	g_list_free_full(vols, g_object_unref);
	g_object_unref(monitor);
	return count;
}

typedef struct {
	SirilFileBrowser *fb;
	GVolume          *volume;  /* owned ref */
	GtkListBoxRow    *row;     /* not owned (weak ref via GWeakRef pattern below) */
	GWeakRef          row_ref;
} MountClosure;

/* In-place conversion of a "Not mounted" sidebar row into a normal
 * filesystem-path row.  Called after a successful mount (whether ours or
 * one that happened externally — e.g. another file manager auto-mounted
 * the same device while our dialog is still open).
 *
 * Mutations:
 *   - Drop the ".siril-unmounted" CSS class (kills the dimming).
 *   - Clear the "Click to mount" tooltip.
 *   - Remove the "Not mounted" subtitle.  sidebar_make_row builds the row
 *     as listrow → hbox → [icon, labels(vbox)], with the subtitle as the
 *     second child of `labels`; walk to it and detach.
 *   - Swap "siril-volume" for "siril-path" so subsequent clicks navigate
 *     instead of trying to mount the volume again. */
static void sidebar_row_demote_to_path(GtkListBoxRow *row, const char *path) {
	if (!row) return;
	gtk_widget_remove_css_class(GTK_WIDGET(row), "siril-unmounted");
	gtk_widget_set_tooltip_text(GTK_WIDGET(row), NULL);
	GtkWidget *hbox = gtk_list_box_row_get_child(row);
	if (hbox) {
		/* Find the labels vbox.  In sidebar_make_row that's the second
		 * child of the hbox (after the icon image). */
		GtkWidget *icon = gtk_widget_get_first_child(hbox);
		GtkWidget *labels = icon ? gtk_widget_get_next_sibling(icon) : NULL;
		if (labels && GTK_IS_BOX(labels)) {
			GtkWidget *main_lbl = gtk_widget_get_first_child(labels);
			GtkWidget *sub_lbl = main_lbl ? gtk_widget_get_next_sibling(main_lbl) : NULL;
			if (sub_lbl) gtk_box_remove(GTK_BOX(labels), sub_lbl);
		}
	}
	/* g_object_set_data(NULL) invokes the previous destroy notifier, so
	 * the held GVolume is unref'd. */
	g_object_set_data(G_OBJECT(row), "siril-volume", NULL);
	if (path && *path)
		g_object_set_data_full(G_OBJECT(row), "siril-path",
		                       g_strdup(path), g_free);
}

static void on_volume_mount_done(GObject *src, GAsyncResult *res, gpointer ud) {
	MountClosure *mc = ud;
	GError *err = NULL;
	gboolean ok = g_volume_mount_finish(G_VOLUME(src), res, &err);

	/* Some "errors" really mean success-by-side-channel: the volume might
	 * already be mounted (because something outside our dialog mounted it
	 * a moment earlier, or because the same row was clicked twice while
	 * the first mount was still running).  If the volume now has a mount,
	 * treat the operation as successful regardless of `ok` — and reuse
	 * the existing mount's root.  This silences the spurious "Device is
	 * already mounted" warning users were seeing on a second click. */
	GMount *mount = g_volume_get_mount(G_VOLUME(src));
	if (mount) {
		GFile *root = g_mount_get_root(mount);
		if (root) {
			navigate_to(mc->fb, root);
			gchar *path = g_file_get_path(root);
			GtkListBoxRow *row = g_weak_ref_get(&mc->row_ref);
			if (row) {
				sidebar_row_demote_to_path(row, path);
				g_object_unref(row);
			}
			g_free(path);
			g_object_unref(root);
		}
		g_object_unref(mount);
		g_clear_error(&err);
	} else if (!ok) {
		/* G_IO_ERROR_FAILED_HANDLED is the user dismissing the auth
		 * dialog — don't spam the log for that. */
		if (!err || !g_error_matches(err, G_IO_ERROR, G_IO_ERROR_FAILED_HANDLED)) {
			gchar *name = g_volume_get_name(G_VOLUME(src));
			siril_log_warning(_("Could not mount %s: %s\n"),
				name ? name : "?", err ? err->message : "unknown");
			g_free(name);
		}
		g_clear_error(&err);
	}
	g_weak_ref_clear(&mc->row_ref);
	g_object_unref(mc->volume);
	g_free(mc);
}

/* Kick off an async mount on `vol`, with a GtkMountOperation so any
 * password / authentication is parented to the file-browser window.  When
 * the mount completes (success), on_volume_mount_done navigates to the
 * mount's root and demotes the sidebar row.  `row` may be NULL when
 * called from a context other than the sidebar — in which case the row
 * mutation is simply skipped. */
static void sidebar_mount_volume_async(SirilFileBrowser *fb, GVolume *vol,
                                       GtkListBoxRow *row) {
	if (!fb || !vol) return;
	GMountOperation *op = gtk_mount_operation_new(GTK_WINDOW(fb->window));
	MountClosure *mc = g_new0(MountClosure, 1);
	mc->fb = fb;
	mc->volume = g_object_ref(vol);
	/* Weak ref so that if the file browser is torn down before the mount
	 * finishes (race on quick Cancel), the row pointer comes back NULL
	 * instead of dangling. */
	g_weak_ref_init(&mc->row_ref, row);
	g_volume_mount(vol, G_MOUNT_MOUNT_NONE, op, NULL,
	               on_volume_mount_done, mc);
	g_object_unref(op);
}

/* Read GTK bookmarks.  Format: one bookmark per line, "<uri> [label]".
 * Tries gtk-4.0 first, falls back to gtk-3.0 then ~/.gtk-bookmarks.
 * Returns the number of rows added. */
static int sidebar_populate_bookmarks(GtkListBox *box) {
	const gchar *config_home = g_get_user_config_dir();
	gchar *candidates[] = {
		g_build_filename(config_home, "gtk-4.0", "bookmarks", NULL),
		g_build_filename(config_home, "gtk-3.0", "bookmarks", NULL),
		g_build_filename(g_get_home_dir(), ".gtk-bookmarks", NULL),
		NULL,
	};
	gchar *contents = NULL;
	gsize len = 0;
	gchar *picked = NULL;
	for (int i = 0; candidates[i]; i++) {
		if (g_file_get_contents(candidates[i], &contents, &len, NULL)) {
			picked = candidates[i];
			break;
		}
	}
	int count = 0;
	if (contents) {
		gchar **lines = g_strsplit(contents, "\n", -1);
		for (int i = 0; lines[i]; i++) {
			gchar *line = g_strstrip(lines[i]);
			if (!*line) continue;
			/* Split on the first space: "<uri> [label...]". */
			gchar *space = strchr(line, ' ');
			gchar *uri, *label;
			if (space) {
				*space = '\0';
				uri = line;
				label = g_strstrip(space + 1);
			} else {
				uri = line;
				label = NULL;
			}
			GFile *gf = g_file_new_for_uri(uri);
			gchar *display = NULL;
			if (label && *label) {
				display = g_strdup(label);
			} else {
				/* No label: use the basename (or a friendly tail of
				 * the URI).  For local paths g_file_get_basename
				 * gives "foo" for /a/b/foo. */
				display = g_file_get_basename(gf);
				if (!display) display = g_file_get_uri(gf);
			}
			/* Skip bookmarks pointing at filesystems we can't see —
			 * deleted local paths produce noise. */
			gboolean ok = TRUE;
			if (g_file_is_native(gf)) {
				gchar *p = g_file_get_path(gf);
				if (!p || !g_file_test(p, G_FILE_TEST_IS_DIR)) ok = FALSE;
				g_free(p);
			}
			if (ok) {
				sidebar_add_gfile(box, "user-bookmarks-symbolic",
				                  display, gf);
				count++;
			}
			g_free(display);
			g_object_unref(gf);
		}
		g_strfreev(lines);
		g_free(contents);
	}
	for (int i = 0; candidates[i]; i++)
		if (candidates[i] != picked) g_free(candidates[i]);
	g_free(picked);
	return count;
}

/* Helper: append a section header and a populator's rows; if the populator
 * adds nothing, drop the header so the sidebar doesn't show an empty
 * heading.  Mirrors the GtkPlacesSidebar behaviour: sections appear only
 * when they have content. */
static void sidebar_section(GtkListBox *list, const char *title,
                            int (*populate)(GtkListBox *)) {
	gint header_pos = 0;
	GtkWidget *r;
	while ((r = GTK_WIDGET(gtk_list_box_get_row_at_index(list, header_pos))))
		header_pos++;
	sidebar_add_section_header(list, title);
	if (populate(list) == 0) {
		GtkListBoxRow *hdr = gtk_list_box_get_row_at_index(list, header_pos);
		if (hdr) gtk_list_box_remove(list, GTK_WIDGET(hdr));
	}
}

/* Mounted + unmounted volumes share the "Devices" section.  Mounted ones
 * render as plain navigation rows; unmounted ones get a "Not mounted"
 * subtitle and the .siril-unmounted class for a dimmed look. */
static int sidebar_populate_all_devices(GtkListBox *box) {
	int n = sidebar_populate_volumes(box);
	n += sidebar_populate_unmounted_volumes(box);
	return n;
}

/* Build a heading label styled like the preview pane's "Preview" header.
 * Used only for the preview pane now — Places is back inside the listbox
 * as a section header alongside Bookmarks / Devices. */
static GtkWidget *make_pane_heading(const char *text) {
	GtkWidget *heading = gtk_label_new(NULL);
	gchar *m = g_markup_printf_escaped("<b>%s</b>", text);
	gtk_label_set_markup(GTK_LABEL(heading), m);
	g_free(m);
	gtk_label_set_xalign(GTK_LABEL(heading), 0.0);
	return heading;
}

static GtkWidget *build_sidebar(SirilFileBrowser *fb) {
	GtkWidget *list = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(list), GTK_SELECTION_NONE);
	gtk_widget_add_css_class(list, "navigation-sidebar");

	/* Sidebar entries.  Order mirrors the GTK3 GtkFileChooser sidebar that
	 * Siril shipped previously, plus the Working Directory shortcut that
	 * jumps to com.wd (the folder Siril is operating on).  Templates /
	 * Public / Computer are intentionally omitted — see the brief. */
	sidebar_add_section_header(GTK_LIST_BOX(list), _("Places"));
	sidebar_add_recent_entry(GTK_LIST_BOX(list));
	sidebar_add_location(GTK_LIST_BOX(list), "user-home-symbolic",
	                     _("Home"), g_get_home_dir());
	sidebar_add_location(GTK_LIST_BOX(list), "folder-documents-symbolic",
	                     _("Documents"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_DOCUMENTS));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-pictures-symbolic",
	                     _("Pictures"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_PICTURES));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-music-symbolic",
	                     _("Music"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_MUSIC));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-download-symbolic",
	                     _("Downloads"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_DOWNLOAD));
	sidebar_add_location(GTK_LIST_BOX(list), "folder-videos-symbolic",
	                     _("Videos"),
	                     g_get_user_special_dir(G_USER_DIRECTORY_VIDEOS));
	{
		GFile *trash = g_file_new_for_uri("trash:///");
		sidebar_add_gfile(GTK_LIST_BOX(list), "user-trash-symbolic",
		                  _("Trash"), trash);
		g_object_unref(trash);
	}
	if (com.wd && *com.wd && g_file_test(com.wd, G_FILE_TEST_IS_DIR))
		sidebar_add_location(GTK_LIST_BOX(list), "folder-symbolic",
		                     _("Working Directory"), com.wd);

	/* Section: Bookmarks (only if the user has any). */
	sidebar_section(GTK_LIST_BOX(list), _("Bookmarks"), sidebar_populate_bookmarks);

	/* Section: Devices (mounted first, then unmounted with a Mount badge). */
	sidebar_section(GTK_LIST_BOX(list), _("Devices"), sidebar_populate_all_devices);

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

/* Render a themed icon into the preview's GtkPicture at native size,
 * centred (SCALE_DOWN), so it doesn't blow up to fill the 280x280 area
 * the way a thumbnail-sized texture would.  Falls back silently if the
 * icon lookup returns NULL (then the picture is just cleared). */
static void picture_set_icon(GtkPicture *picture, const char *icon_name) {
	if (!picture || !icon_name) return;
	GdkDisplay *disp = gtk_widget_get_display(GTK_WIDGET(picture));
	if (!disp) { gtk_picture_set_paintable(picture, NULL); return; }
	GtkIconTheme *theme = gtk_icon_theme_get_for_display(disp);
	/* Request a large icon so themed SVGs render crisply when scaled to
	 * fit the preview area; raster fallbacks stay close to native size. */
	GtkIconPaintable *icon = gtk_icon_theme_lookup_icon(theme, icon_name,
		NULL, 256, 1, GTK_TEXT_DIR_NONE, 0);
	if (!icon) { gtk_picture_set_paintable(picture, NULL); return; }
	gtk_picture_set_paintable(picture, GDK_PAINTABLE(icon));
	/* SCALE_DOWN, like the thumbnail path: show the icon at (up to) its
	 * looked-up native size (256 px) centred in the pane, rather than
	 * blowing it up to fill the whole preview area — keeps the icon within
	 * the same size constraints as image thumbnails. */
	gtk_picture_set_content_fit(picture, GTK_CONTENT_FIT_SCALE_DOWN);
	g_object_unref(icon);
}

/* Types that GdkTexture / gdk-pixbuf can never load.  Hitting
 * gdk_texture_new_from_file on these is wasteful at best and
 * pathologically slow at worst — pixbuf walks every registered loader
 * probing magic bytes, and some loaders read large chunks before
 * deciding the file isn't theirs.  For multi-GB SER / AVI files that
 * means a multi-second stall on every selection change. */
static gboolean type_known_non_raster(image_type t) {
	switch (t) {
	case TYPEAVI:
	case TYPESER:
	case TYPERAW:
	case TYPEXISF:
	case TYPEPIC:
		return TRUE;
	default:
		return FALSE;
	}
}

/* Pick a fallback icon name for files we couldn't generate a thumbnail
 * for.  We split by image_type so AVI/SER get the video glyph, raster
 * and FITS get the image glyph, and everything else falls back to a
 * generic file icon. */
static const char *fallback_icon_for_type(image_type t) {
	switch (t) {
	case TYPEAVI:
	case TYPESER:
		return "video-x-generic";
	case TYPEFITS:
	case TYPETIFF:
	case TYPEBMP:
	case TYPEPNG:
	case TYPEJPG:
	case TYPEHEIF:
	case TYPEPNM:
	case TYPEPIC:
	case TYPERAW:
	case TYPEXISF:
	case TYPEJXL:
	case TYPEAVIF:
		return "image-x-generic";
	default:
		return "text-x-generic";
	}
}

/* Map a Siril image_type to a short uppercase format label. */
static const char *type_short_label(image_type t, const char *path) {
	switch (t) {
	case TYPEFITS: return "FITS";
	case TYPETIFF: return "TIFF";
	case TYPEBMP:  return "BMP";
	case TYPEPNG:  return "PNG";
	case TYPEJPG:  return "JPEG";
	case TYPEHEIF: return "HEIF";
	case TYPEPNM:  return "PNM";
	case TYPEPIC:  return "PIC";
	case TYPERAW:  return "RAW";
	case TYPEXISF: return "XISF";
	case TYPEJXL:  return "JXL";
	case TYPEAVIF: return "AVIF";
	case TYPEAVI:  return "AVI";
	case TYPESER:  return "SER";
	default: break;
	}
	/* TYPEUNDEF — fall back to the extension if any. */
	if (path) {
		const char *ext = strrchr(path, '.');
		if (ext && *(ext + 1)) {
			static char buf[16];
			g_snprintf(buf, sizeof(buf), "%s", ext + 1);
			for (char *c = buf; *c; c++) *c = g_ascii_toupper(*c);
			return buf;
		}
	}
	return _("Unknown");
}

/* Middle-truncate a (UTF-8) basename with an internal ellipsis so a long
 * filename stays on a single line under the preview instead of wrapping
 * across several: "verylongprefix…_suffix.fits".  Both the leading text and
 * the trailing text (including the extension) are kept.  The budget is sized
 * from the preview pane width (≈ thumbnail_size px) so it roughly matches how
 * many average glyphs actually fit.  Caller frees the result. */
static gchar *ellipsize_basename_middle(const gchar *base) {
	glong max_chars = CLAMP(com.pref.gui.thumbnail_size / 8, 20, 60);
	glong len = g_utf8_strlen(base, -1);
	if (len <= max_chars)
		return g_strdup(base);
	glong keep = max_chars - 1;          /* reserve one glyph for the ellipsis */
	glong head = (keep + 1) / 2;         /* favour the head when keep is odd */
	glong tail = keep - head;
	const gchar *head_end   = g_utf8_offset_to_pointer(base, head);
	const gchar *tail_start = g_utf8_offset_to_pointer(base, len - tail);
	gchar *head_str = g_strndup(base, head_end - base);
	gchar *out = g_strconcat(head_str, "\xE2\x80\xA6" /* U+2026 … */, tail_start, NULL);
	g_free(head_str);
	return out;
}

void siril_file_browser_default_preview(const gchar *path,
                                        GtkPicture  *picture,
                                        GtkLabel    *metadata_label,
                                        gpointer     user_data) {
	(void)user_data;
	gtk_picture_set_paintable(picture, NULL);
	if (metadata_label) gtk_label_set_text(metadata_label, "");
	if (!path) return;

	/* Directory: show a folder glyph and the basename — no further info
	 * lookup needed (and we don't want to do an expensive recursive size
	 * count). */
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		picture_set_icon(picture, "folder");
		if (metadata_label) {
			gchar *base = g_path_get_basename(path);
			gchar *shrt = ellipsize_basename_middle(base ? base : "");
			gchar *esc  = g_markup_escape_text(shrt, -1);
			g_free(shrt);
			/* Name above type: bold filename on top, plain "Folder" beneath. */
			gchar *m = g_strdup_printf("<b>%s</b>\n%s", esc, _("Folder"));
			gtk_label_set_markup(metadata_label, m);
			g_free(m);
			g_free(esc);
			g_free(base);
		}
		return;
	}

	image_type itype = get_type_from_filename(path);
	const char *type_label = type_short_label(itype, path);

	/* On-disk size — always available via stat. */
	gchar *size_str = NULL;
	{
		GFile *gf = g_file_new_for_path(path);
		GFileInfo *fi = g_file_query_info(gf, G_FILE_ATTRIBUTE_STANDARD_SIZE,
		                                  G_FILE_QUERY_INFO_NONE, NULL, NULL);
		if (fi) {
			size_str = g_format_size(g_file_info_get_size(fi));
			g_object_unref(fi);
		}
		g_object_unref(gf);
	}

	/* Try to fill in dimensions / extra info via a preview load. */
	gchar *dims_str = NULL;        /* "1024 x 1024 pixels" */
	gchar *extra_str = NULL;       /* e.g. "3 channels (16 bits)" */
	gboolean got_texture = FALSE;

	/* FITS / SER / AVI thumbnails go through Siril-side extractors.  All
	 * three follow the same contract: return a malloc'd RGB888 buffer
	 * plus dimensions plus a multi-line description.  The FITS extractor
	 * recognises .fit / .fits / .fits.fz / etc; the SER one reads frame 0
	 * with an MTF stretch for higher bit depths; the AVI one decodes the
	 * first video frame via libavformat without a full FFMS2 index. */
	if (itype == TYPEFITS || itype == TYPESER || itype == TYPEAVI) {
		gchar *descr = NULL;
		int w = 0, h = 0;
		guchar *data = NULL;
		if (itype == TYPEFITS)
			data = extract_thumbnail_from_fits(path, &descr, &w, &h);
		else if (itype == TYPESER)
			data = extract_thumbnail_from_ser(path, &descr, &w, &h);
		else /* TYPEAVI */
			data = extract_thumbnail_from_avi(path, &descr, &w, &h);
		if (data) {
			GdkTexture *tex = siril_texture_from_rgb_bytes(data,
				(gsize)w * h * 3, w, h, w * 3, FALSE,
				(GDestroyNotify)free, data);
			if (tex) {
				gtk_picture_set_paintable(picture, GDK_PAINTABLE(tex));
				/* SCALE_DOWN, not CONTAIN: show the thumbnail at its native
				 * rendered size (thumbnail_size) rather than upscaling it
				 * blurry to fill the pane. */
				gtk_picture_set_content_fit(picture, GTK_CONTENT_FIT_SCALE_DOWN);
				g_object_unref(tex);
				got_texture = TRUE;
			} else {
				free(data);
			}
		}
		if (descr && *descr) {
			/* The descriptor's first line is "WxH pixels"; remaining
			 * lines hold channel/bit/frame info — split so we can lay
			 * it out predictably (type, dims, extra, size). */
			gchar **lines = g_strsplit(descr, "\n", 2);
			if (lines[0] && *lines[0]) dims_str  = g_strdup(lines[0]);
			if (lines[0] && lines[1]) extra_str = g_strdup(lines[1]);
			g_strfreev(lines);
		}
		g_free(descr);
	} else if (!type_known_non_raster(itype)) {
		int max_size = com.pref.gui.thumbnail_size;  /* 128, 256 or 512 */
		/* Original dimensions for the metadata line, read from the header
		 * without decoding the pixels. */
		int orig_w = 0, orig_h = 0;
		gdk_pixbuf_get_file_info(path, &orig_w, &orig_h);
		GError *err = NULL;
		GdkPixbuf *pix = gdk_pixbuf_new_from_file_at_scale(path, max_size,
			max_size, TRUE, &err);
		if (pix) {
			/* Wrap the pixbuf's pixels in a GdkMemoryTexture (the
			 * non-deprecated path; gdk_texture_new_for_pixbuf is deprecated).
			 * The GBytes keeps the pixbuf alive for as long as the texture
			 * references its pixel data — ownership of our pixbuf ref is
			 * handed to the bytes' free func, so we don't unref pix here. */
			int pw = gdk_pixbuf_get_width(pix);
			int ph = gdk_pixbuf_get_height(pix);
			int stride = gdk_pixbuf_get_rowstride(pix);
			GdkMemoryFormat fmt = gdk_pixbuf_get_has_alpha(pix)
				? GDK_MEMORY_R8G8B8A8 : GDK_MEMORY_R8G8B8;
			GBytes *bytes = g_bytes_new_with_free_func(
				gdk_pixbuf_get_pixels(pix), (gsize)stride * ph,
				(GDestroyNotify)g_object_unref, pix);
			GdkTexture *tex = gdk_memory_texture_new(pw, ph, fmt, bytes, stride);
			g_bytes_unref(bytes);  /* texture holds its own ref */
			if (tex) {
				gtk_picture_set_paintable(picture, GDK_PAINTABLE(tex));
				gtk_picture_set_content_fit(picture, GTK_CONTENT_FIT_SCALE_DOWN);
				if (orig_w > 0 && orig_h > 0)
					dims_str = g_strdup_printf("%d x %d %s",
						orig_w, orig_h, _("pixels"));
				else
					dims_str = g_strdup_printf("%d x %d %s",
						pw, ph, _("pixels"));
				g_object_unref(tex);
				got_texture = TRUE;
			}
		} else {
			g_clear_error(&err);
		}
	}

	/* No thumbnail available — show a type-appropriate icon so the
	 * preview area isn't a blank rectangle. */
	if (!got_texture)
		picture_set_icon(picture, fallback_icon_for_type(itype));

	if (metadata_label) {
		GString *meta = g_string_new(NULL);
		/* Name line — bold so it reads as the dominant identifier.  Name
		 * goes above type (per UX brief) so the user's eye lands on the
		 * file they've selected before reading format / dimensions. */
		gchar *base = g_path_get_basename(path);
		if (base && *base) {
			gchar *shrt = ellipsize_basename_middle(base);
			gchar *e = g_markup_escape_text(shrt, -1);
			g_string_append_printf(meta, "<b>%s</b>", e);
			g_free(e);
			g_free(shrt);
		}
		g_free(base);
		gchar *type_esc = g_markup_escape_text(type_label, -1);
		g_string_append_printf(meta, "%s%s", meta->len ? "\n" : "", type_esc);
		g_free(type_esc);
		if (dims_str && *dims_str) {
			gchar *e = g_markup_escape_text(dims_str, -1);
			g_string_append_printf(meta, "\n%s", e);
			g_free(e);
		}
		if (extra_str && *extra_str) {
			gchar *e = g_markup_escape_text(extra_str, -1);
			g_string_append_printf(meta, "\n%s", e);
			g_free(e);
		}
		if (size_str && *size_str) {
			gchar *e = g_markup_escape_text(size_str, -1);
			g_string_append_printf(meta, "\n%s", e);
			g_free(e);
		}
		gtk_label_set_markup(metadata_label, meta->str);
		g_string_free(meta, TRUE);
	}

	g_free(dims_str);
	g_free(extra_str);
	g_free(size_str);
}

/* ── debayer toggle (cross-tab sync) ──────────────────────────────── */

/* When the user clicks the browser's Debayer check, write the pref so
 * future opens honour it, and mirror the state onto the Convert tab's
 * demosaicingButton.  Mirroring fires the latter's own "toggled" handler
 * (on_demosaicing_toggled in conversion.c) but that handler just re-
 * writes the same pref, so there's no real loop. */
static void on_browser_debayer_toggled(GtkCheckButton *cb, gpointer ud) {
	SirilFileBrowser *fb = ud;
	gboolean active = gtk_check_button_get_active(cb);
	com.pref.debayer.open_debayer = active;
	if (fb->external_demosaic_btn &&
	    gtk_check_button_get_active(fb->external_demosaic_btn) != active) {
		gtk_check_button_set_active(fb->external_demosaic_btn, active);
	}
}

/* Reverse direction: if the Convert-tab toggle changes while the browser
 * is open (rare, but possible via scripting / hotkeys), mirror it back. */
static void on_external_demosaic_toggled(GtkCheckButton *cb, gpointer ud) {
	SirilFileBrowser *fb = ud;
	gboolean active = gtk_check_button_get_active(cb);
	if (fb->debayer_check &&
	    gtk_check_button_get_active(fb->debayer_check) != active) {
		/* Block our own handler while we mirror, so we don't write the
		 * pref twice. */
		if (fb->debayer_check_handler)
			g_signal_handler_block(fb->debayer_check, fb->debayer_check_handler);
		gtk_check_button_set_active(fb->debayer_check, active);
		if (fb->debayer_check_handler)
			g_signal_handler_unblock(fb->debayer_check, fb->debayer_check_handler);
	}
}

void siril_file_browser_set_show_debayer_toggle(SirilFileBrowser *fb, gboolean show) {
	if (!fb || !fb->debayer_check) return;
	fb->show_debayer_toggle = show;
	gtk_widget_set_visible(GTK_WIDGET(fb->debayer_check), show);
	if (!show) return;

	/* Initial state from the pref. */
	gtk_check_button_set_active(fb->debayer_check, com.pref.debayer.open_debayer);

	/* Connect our own toggle handler once. */
	if (!fb->debayer_check_handler) {
		fb->debayer_check_handler = g_signal_connect(fb->debayer_check,
			"toggled", G_CALLBACK(on_browser_debayer_toggled), fb);
	}

	/* Bind to the Convert tab's demosaicingButton (if the main UI builder
	 * is available — it is for any GUI session, but guard for unit tests
	 * that construct a browser without a full UI tree). */
	if (!fb->external_demosaic_btn && gui.builder) {
		GObject *o = gtk_builder_get_object(gui.builder, "demosaicingButton");
		if (o && GTK_IS_CHECK_BUTTON(o)) {
			fb->external_demosaic_btn = GTK_CHECK_BUTTON(o);
			fb->external_demosaic_handler = g_signal_connect(
				fb->external_demosaic_btn, "toggled",
				G_CALLBACK(on_external_demosaic_toggled), fb);
		}
	}
}

/* ── construction / API ────────────────────────────────────────────── */

/* Read the system-wide "show hidden files" setting from
 * org.gtk.Settings.FileChooser (the same key GTK's own GtkFileDialog
 * honors).  Returns FALSE when the schema is unavailable (typical on
 * macOS / Windows where there's no dconf backend), preserving the
 * historical hide-dotfiles behavior. */
static gboolean browser_initial_show_hidden(void) {
	GSettingsSchemaSource *src = g_settings_schema_source_get_default();
	if (!src) return FALSE;
	GSettingsSchema *schema = g_settings_schema_source_lookup(src,
		"org.gtk.Settings.FileChooser", TRUE);
	if (!schema) return FALSE;
	gboolean has_key = g_settings_schema_has_key(schema, "show-hidden");
	g_settings_schema_unref(schema);
	if (!has_key) return FALSE;
	GSettings *gs = g_settings_new("org.gtk.Settings.FileChooser");
	gboolean v = g_settings_get_boolean(gs, "show-hidden");
	g_object_unref(gs);
	return v;
}

/* Ctrl+H toggles hidden-file visibility, matching every other GTK file
 * chooser.  We refilter the model so the change is immediate. */
static gboolean browser_window_key_pressed(GtkEventControllerKey *kc, guint keyval,
                                           guint keycode, GdkModifierType state,
                                           gpointer ud) {
	(void)kc; (void)keycode;
	SirilFileBrowser *fb = ud;
	/* Primary+H toggles hidden files (Cmd+H on macOS, Ctrl+H elsewhere). */
	if ((state & get_primary()) && (keyval == GDK_KEY_h || keyval == GDK_KEY_H)) {
		fb->show_hidden_files = !fb->show_hidden_files;
		if (fb->file_filter)
			gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);
		return TRUE;
	}
	/* Primary+L toggles the path into text-entry mode, matching the GTK3
	 * GtkFileChooser shortcut (and the edit-path toolbar button).  Uses the
	 * platform primary modifier so it is Cmd+L on macOS, Ctrl+L elsewhere. */
	if ((state & get_primary()) && (keyval == GDK_KEY_l || keyval == GDK_KEY_L)) {
		on_edit_path_clicked(NULL, fb);
		return TRUE;
	}
	/* Escape dismisses the dialog */
	if (keyval == GDK_KEY_Escape) {
		if (fb->path_stack &&
		    g_strcmp0(gtk_stack_get_visible_child_name(fb->path_stack), "entry") == 0)
			return FALSE;
		if (fb->search_bar && gtk_search_bar_get_search_mode(fb->search_bar))
			return FALSE;
		fb->response = GTK_RESPONSE_CANCEL;
		if (fb->loop) g_main_loop_quit(fb->loop);
		return TRUE;
	}
	return FALSE;
}

/* The file browser is a process-wide singleton.  We never destroy it —
 * the window is built once on first use and hidden on close (see
 * gtk_window_set_hide_on_close + siril_file_browser_run's set_visible).
 * Each siril_file_browser_new() call resets per-open state on the same
 * instance and reconfigures title + transient-for.  This sidesteps the
 * macOS AppKit teardown problems (multiple earlier "fix the destroy
 * lifecycle" attempts all failed: drain, modal/transient toggling,
 * defer-load, drop-modal) by removing the destroy step entirely. */
static SirilFileBrowser *g_singleton = NULL;

static void build_browser_widgets(SirilFileBrowser *fb);
static void reset_browser_state(SirilFileBrowser *fb);

SirilFileBrowser *siril_file_browser_new(GtkWindow *parent, const gchar *title) {
	if (!g_singleton) {
		SirilFileBrowser *fb = g_new0(SirilFileBrowser, 1);
		fb->filters = g_ptr_array_new_with_free_func((GDestroyNotify)browser_filter_free);
		fb->back_history    = g_ptr_array_new_with_free_func(g_object_unref);
		fb->forward_history = g_ptr_array_new_with_free_func(g_object_unref);
		fb->active_filter_idx = -1;
		fb->response = GTK_RESPONSE_NONE;
		fb->show_hidden_files = browser_initial_show_hidden();

		fb->window = GTK_WINDOW(gtk_window_new());
		/* Hide on close, never destroy — see comment on g_singleton above. */
		gtk_window_set_hide_on_close(fb->window, TRUE);
		/* Restore the size the user last left the dialog at (persisted in
		 * the init file so it survives across sessions and is portable —
		 * unlike the GTK3 GtkFileChooser which leaned on the desktop's own
		 * geometry store).  0 means "never saved", so fall back to the
		 * built-in default. */
		int init_w = (com.pref.gui.remember_windows && com.pref.gui.open_dialog_w > 0)
			? com.pref.gui.open_dialog_w : 1000;
		int init_h = (com.pref.gui.remember_windows && com.pref.gui.open_dialog_h > 0)
			? com.pref.gui.open_dialog_h : 620;
		gtk_window_set_default_size(fb->window, init_w, init_h);
		build_browser_widgets(fb);
		g_singleton = fb;
	} else {
		reset_browser_state(g_singleton);
	}

	SirilFileBrowser *fb = g_singleton;
	fb->parent = parent;
	gtk_window_set_title(fb->window, title ? title : _("Open File"));
	if (parent) {
		gtk_window_set_transient_for(fb->window, parent);
		/* Plan C lives on: block parent input via insensitive instead
		 * of gtk_window_set_modal(TRUE).  See g_singleton comment above
		 * for the long history of attempts at the modal path on macOS. */
		gtk_widget_set_sensitive(GTK_WIDGET(parent), FALSE);
	} else {
		gtk_window_set_transient_for(fb->window, NULL);
	}
	return fb;
}

static void build_browser_widgets(SirilFileBrowser *fb) {
	{
		GtkEventController *kc = gtk_event_controller_key_new();
		/* CAPTURE phase: intercept Ctrl+H before GtkColumnView consumes it.
		 * In GTK4 the default BUBBLE phase never reaches the window because
		 * the focused column view handles (and drops) the key event first. */
		gtk_event_controller_set_propagation_phase(kc, GTK_PHASE_CAPTURE);
		g_signal_connect(kc, "key-pressed",
		                 G_CALLBACK(browser_window_key_pressed), fb);
		gtk_widget_add_controller(GTK_WIDGET(fb->window), kc);
	}
	gtk_widget_set_size_request(GTK_WIDGET(fb->window), 880, 520);

	/* Title-only header bar: shows the dialog title centred and the window
	 * controls (close / minimise, native traffic lights on macOS).  The
	 * action buttons themselves (Cancel / Open) live in the toolbar row
	 * just below, on the same line as the path breadcrumb. */
	GtkWidget *header = gtk_header_bar_new();
	gtk_header_bar_set_show_title_buttons(GTK_HEADER_BAR(header), TRUE);
#if defined(OS_OSX) && GTK_CHECK_VERSION(4, 18, 0)
	gtk_header_bar_set_use_native_controls(GTK_HEADER_BAR(header), TRUE);
#endif
	gtk_window_set_titlebar(fb->window, header);

	/* Toolbar (single row): Cancel | breadcrumb path | edit | search | Open.
	 * Cancel anchors the left and Open the right, with the navigation
	 * cluster between them — every control the user acts on lives on one
	 * line, directly under the title bar. */
	GtkWidget *toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	gtk_widget_set_margin_start(toolbar, 6);
	gtk_widget_set_margin_end(toolbar, 6);
	gtk_widget_set_margin_top(toolbar, 6);
	gtk_widget_set_margin_bottom(toolbar, 3);

	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	g_signal_connect(cancel, "clicked", G_CALLBACK(on_cancel_clicked), fb);
	gtk_widget_set_margin_end(cancel, 8);
	gtk_box_append(GTK_BOX(toolbar), cancel);

	/* Back / Forward / Up buttons removed — the breadcrumb's own
	 * left/right arrows (created below as fb->breadcrumb_left_btn /
	 * _right_btn) cover the same navigation actions, and the trail
	 * itself visualises history far better than a stack-of-frames
	 * back button could.  The fb->back_button / forward_button /
	 * up_button fields stay NULL; update_nav_sensitivity is NULL-safe. */
	/* Path display: a two-page GtkStack.
	 *
	 * "breadcrumb" page is a single styled bar (`.siril-breadcrumb-bar`)
	 * laid out as:
	 *     [ ← ][ scroll-clipped segment buttons ][ → ]
	 * The arrow buttons drive the inner scrolled window's hadjustment
	 * when the breadcrumb overflows; the scrolled window itself uses
	 * GTK_POLICY_EXTERNAL so no native scrollbar appears.  Segments
	 * have no chevron separators — they butt up against each other and
	 * rely on the bar background / hover contrast for separation.  The
	 * rightmost (current) segment carries `.suggested-action` for the
	 * blue active-cell highlight.
	 *
	 * "entry" page is the editable GtkEntry — surfaced via the edit
	 * toggle to the right or by pressing Escape (toggles back).  Enter
	 * in entry mode navigates and flips back to breadcrumb. */
	fb->breadcrumb_scroller = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(fb->breadcrumb_scroller),
	                               GTK_POLICY_EXTERNAL, GTK_POLICY_NEVER);
	gtk_widget_set_hexpand(fb->breadcrumb_scroller, TRUE);
	fb->breadcrumb_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_add_css_class(fb->breadcrumb_box, "siril-breadcrumb");
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(fb->breadcrumb_scroller),
	                              fb->breadcrumb_box);

	/* Trail-navigation arrows — left + right ends of the bar.
	 *   left  = step up one level (navigate to current_folder's parent)
	 *   right = step down one level along the persistent trail (the
	 *           next deeper folder we previously visited but have since
	 *           backed away from).
	 * These replace the old back/forward/up buttons that used to live in
	 * the toolbar — the breadcrumb trail already shows the navigation
	 * context, so dedicated history buttons are redundant. */
	fb->breadcrumb_left_btn = gtk_button_new_from_icon_name("go-previous-symbolic");
	gtk_widget_add_css_class(fb->breadcrumb_left_btn, "flat");
	gtk_widget_set_focusable(fb->breadcrumb_left_btn, FALSE);
	gtk_widget_set_can_focus(fb->breadcrumb_left_btn, FALSE);
	gtk_widget_set_tooltip_text(fb->breadcrumb_left_btn, _("Up a level"));
	g_signal_connect(fb->breadcrumb_left_btn, "clicked",
	                 G_CALLBACK(on_breadcrumb_step_up), fb);

	fb->breadcrumb_right_btn = gtk_button_new_from_icon_name("go-next-symbolic");
	gtk_widget_add_css_class(fb->breadcrumb_right_btn, "flat");
	gtk_widget_set_focusable(fb->breadcrumb_right_btn, FALSE);
	gtk_widget_set_can_focus(fb->breadcrumb_right_btn, FALSE);
	gtk_widget_set_tooltip_text(fb->breadcrumb_right_btn,
		_("Forward along the trail"));
	g_signal_connect(fb->breadcrumb_right_btn, "clicked",
	                 G_CALLBACK(on_breadcrumb_step_down), fb);

	GtkWidget *breadcrumb_bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_widget_add_css_class(breadcrumb_bar, "siril-breadcrumb-bar");
	gtk_box_append(GTK_BOX(breadcrumb_bar), fb->breadcrumb_left_btn);
	gtk_box_append(GTK_BOX(breadcrumb_bar), fb->breadcrumb_scroller);
	gtk_box_append(GTK_BOX(breadcrumb_bar), fb->breadcrumb_right_btn);

	fb->path_entry = GTK_ENTRY(gtk_entry_new());
	gtk_entry_set_placeholder_text(fb->path_entry,
		_("Type a path and press Enter to navigate"));
	g_signal_connect(fb->path_entry, "activate",
	                 G_CALLBACK(on_path_entry_activate), fb);
	{
		GtkEventController *kc = gtk_event_controller_key_new();
		g_signal_connect(kc, "key-pressed",
		                 G_CALLBACK(on_path_entry_key_pressed), fb);
		gtk_widget_add_controller(GTK_WIDGET(fb->path_entry), kc);
	}

	fb->path_stack = GTK_STACK(gtk_stack_new());
	gtk_stack_set_transition_type(fb->path_stack, GTK_STACK_TRANSITION_TYPE_NONE);
	gtk_widget_set_hexpand(GTK_WIDGET(fb->path_stack), TRUE);
	gtk_widget_set_margin_start(GTK_WIDGET(fb->path_stack), 6);
	gtk_stack_add_named(fb->path_stack, breadcrumb_bar, "breadcrumb");
	gtk_stack_add_named(fb->path_stack, GTK_WIDGET(fb->path_entry), "entry");
	gtk_stack_set_visible_child_name(fb->path_stack, "breadcrumb");
	gtk_box_append(GTK_BOX(toolbar), GTK_WIDGET(fb->path_stack));

	/* Edit-path toggle: clicking swaps breadcrumb ↔ entry mode. */
	fb->edit_path_btn = gtk_button_new_from_icon_name("document-edit-symbolic");
	gtk_widget_set_tooltip_text(fb->edit_path_btn, _("Edit path as text"));
	gtk_widget_add_css_class(fb->edit_path_btn, "flat");
	gtk_widget_set_margin_start(fb->edit_path_btn, 6);
	g_signal_connect(fb->edit_path_btn, "clicked",
	                 G_CALLBACK(on_edit_path_clicked), fb);
	gtk_box_append(GTK_BOX(toolbar), fb->edit_path_btn);

	fb->new_folder_btn = gtk_button_new_from_icon_name("folder-new-symbolic");
	gtk_widget_set_tooltip_text(fb->new_folder_btn, _("Create folder"));
	gtk_widget_add_css_class(fb->new_folder_btn, "flat");
	gtk_widget_set_margin_start(fb->new_folder_btn, 6);
	gtk_widget_set_visible(fb->new_folder_btn, FALSE);
	gtk_box_append(GTK_BOX(toolbar), fb->new_folder_btn);

	fb->new_folder_popover = GTK_POPOVER(gtk_popover_new());
	gtk_widget_set_parent(GTK_WIDGET(fb->new_folder_popover), fb->new_folder_btn);
	{
		GtkWidget *nf_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		fb->new_folder_entry = GTK_ENTRY(gtk_entry_new());
		gtk_entry_set_placeholder_text(fb->new_folder_entry, _("Folder name"));
		gtk_entry_set_activates_default(fb->new_folder_entry, FALSE);
		gtk_widget_set_hexpand(GTK_WIDGET(fb->new_folder_entry), TRUE);
		g_signal_connect(fb->new_folder_entry, "activate",
		                 G_CALLBACK(on_new_folder_entry_activate), fb);
		GtkWidget *nf_create = gtk_button_new_with_label(_("Create"));
		gtk_widget_add_css_class(nf_create, "suggested-action");
		g_signal_connect(nf_create, "clicked",
		                 G_CALLBACK(on_new_folder_confirm_clicked), fb);
		gtk_box_append(GTK_BOX(nf_box), GTK_WIDGET(fb->new_folder_entry));
		gtk_box_append(GTK_BOX(nf_box), nf_create);
		gtk_popover_set_child(fb->new_folder_popover, nf_box);
	}
	g_signal_connect(fb->new_folder_popover, "show",
	                 G_CALLBACK(on_new_folder_popover_show), fb);
	g_signal_connect_swapped(fb->new_folder_btn, "clicked",
	                         G_CALLBACK(gtk_popover_popup), fb->new_folder_popover);

	/* Search toggle — sits to the left of Open in the top toolbar.  Bound
	 * bidirectionally to the search bar's `search-mode-enabled` so clicking
	 * the button reveals/hides the bar, and Escape (handled by the bar
	 * itself) flips the button back off. */
	GtkWidget *search_button = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(search_button), "system-search-symbolic");
	gtk_widget_set_tooltip_text(search_button, _("Filter the list by filename"));
	gtk_widget_set_margin_start(search_button, 6);
	gtk_box_append(GTK_BOX(toolbar), search_button);

	fb->open_button = gtk_button_new_with_label(_("Open"));
	gtk_widget_add_css_class(fb->open_button, "suggested-action");
	gtk_widget_set_sensitive(fb->open_button, FALSE);
	gtk_widget_set_margin_start(fb->open_button, 8);
	g_signal_connect(fb->open_button, "clicked", G_CALLBACK(on_open_clicked), fb);
	gtk_box_append(GTK_BOX(toolbar), fb->open_button);

	/* Search bar (initially collapsed).  Lives between the toolbar and the
	 * paned area so revealing it slides the file list down without
	 * affecting the sidebar.  The bar handles Escape-to-close natively. */
	fb->search_entry = GTK_SEARCH_ENTRY(gtk_search_entry_new());
	gtk_widget_set_hexpand(GTK_WIDGET(fb->search_entry), TRUE);
	/* GtkSearchEntry is NOT a GtkEntry in GTK4 — it implements GtkEditable
	 * directly without inheriting from GtkEntry — so the entry placeholder
	 * setter trips a GTK_IS_ENTRY assertion.  Use the search-entry-specific
	 * placeholder setter instead. */
	gtk_search_entry_set_placeholder_text(fb->search_entry,
	                                      _("Filter by name…"));
	g_signal_connect(fb->search_entry, "search-changed",
	                 G_CALLBACK(on_search_changed), fb);

	fb->search_bar = GTK_SEARCH_BAR(gtk_search_bar_new());
	gtk_search_bar_set_child(fb->search_bar, GTK_WIDGET(fb->search_entry));
	gtk_search_bar_connect_entry(fb->search_bar,
	                             GTK_EDITABLE(fb->search_entry));
	/* show_close_button=FALSE: the toolbar toggle is our open/close. */
	gtk_search_bar_set_show_close_button(fb->search_bar, FALSE);
	g_signal_connect(fb->search_bar, "notify::search-mode-enabled",
	                 G_CALLBACK(on_search_mode_changed), fb);
	g_object_bind_property(search_button, "active",
	                       fb->search_bar, "search-mode-enabled",
	                       G_BINDING_BIDIRECTIONAL | G_BINDING_SYNC_CREATE);

	/* Models: DirectoryList → FilterListModel → SortListModel →
	 *         SingleSelection → ColumnView.
	 *
	 * Reference ownership: every chain step takes a transfer-full model
	 * arg, and we want to keep our own ref too — so we g_object_ref each
	 * field that we hand off, and we g_object_ref the selection again
	 * when handing it to gtk_column_view_new (also transfer-full).
	 *
	 * Sort is driven by a GtkMultiSorter: first child pins directories to
	 * the top, second child is the GtkColumnView's own sorter — which
	 * updates automatically when the user clicks a column header.
	 */
	/* Request the modified-time attribute too — needed for the Modified
	 * column and date sort. */
	fb->dir_list = gtk_directory_list_new("standard::*,time::modified", NULL);
	/* Live-update the listing when files are added/modified/removed in
	 * the current folder.  Default is already TRUE on construction, but
	 * we set it explicitly so the intent is visible. */
	gtk_directory_list_set_monitored(fb->dir_list, TRUE);
	fb->file_filter  = gtk_custom_filter_new(fileinfo_filter_func, fb, NULL);
	fb->filter_model = gtk_filter_list_model_new(
		G_LIST_MODEL(g_object_ref(fb->dir_list)),
		GTK_FILTER(g_object_ref(fb->file_filter)));
	/* Sort AFTER filter so we don't waste cycles ordering hidden items.
	 * The actual sorter is installed below once the column view exists. */
	fb->sort_model   = gtk_sort_list_model_new(
		G_LIST_MODEL(g_object_ref(fb->filter_model)),
		NULL);
	{
		GtkSingleSelection *ss = gtk_single_selection_new(
			G_LIST_MODEL(g_object_ref(fb->sort_model)));
		gtk_single_selection_set_autoselect(ss, FALSE);
		gtk_single_selection_set_can_unselect(ss, TRUE);
		fb->selection = GTK_SELECTION_MODEL(ss);
	}

	/* ColumnView with four sortable columns: Name | Modified | Type | Size.
	 * The dense-rows CSS class trims default Adwaita row padding so the
	 * picker can show ~50% more files in the same window height. */
	fb->columnview = GTK_COLUMN_VIEW(
		gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(fb->selection))));
	gtk_column_view_set_show_column_separators(fb->columnview, FALSE);
	gtk_column_view_set_show_row_separators(fb->columnview, FALSE);
	gtk_widget_add_css_class(GTK_WIDGET(fb->columnview), "siril-dense-rows");
	/* Column order: Name | Size | Type | Modified. */
	{
		GtkSignalListItemFactory *fname = GTK_SIGNAL_LIST_ITEM_FACTORY(
			gtk_signal_list_item_factory_new());
		g_signal_connect(fname, "setup", G_CALLBACK(name_setup_cb), NULL);
		g_signal_connect(fname, "bind",  G_CALLBACK(name_bind_cb),  fb);
		GtkColumnViewColumn *cn = gtk_column_view_column_new(
			_("Name"), GTK_LIST_ITEM_FACTORY(fname));
		gtk_column_view_column_set_resizable(cn, TRUE);
		gtk_column_view_column_set_expand(cn, TRUE);
		gtk_column_view_column_set_sorter(cn,
			GTK_SORTER(gtk_custom_sorter_new(name_column_compare, NULL, NULL)));
		gtk_column_view_append_column(fb->columnview, cn);
		g_object_unref(cn);
	}
	{
		GtkSignalListItemFactory *fsize = GTK_SIGNAL_LIST_ITEM_FACTORY(
			gtk_signal_list_item_factory_new());
		g_signal_connect(fsize, "setup", G_CALLBACK(size_setup_cb), NULL);
		g_signal_connect(fsize, "bind",  G_CALLBACK(size_bind_cb),  fb);
		GtkColumnViewColumn *cs = gtk_column_view_column_new(
			_("Size"), GTK_LIST_ITEM_FACTORY(fsize));
		gtk_column_view_column_set_resizable(cs, TRUE);
		gtk_column_view_column_set_sorter(cs,
			GTK_SORTER(gtk_custom_sorter_new(size_column_compare, NULL, NULL)));
		gtk_column_view_append_column(fb->columnview, cs);
		g_object_unref(cs);
	}
	{
		GtkSignalListItemFactory *ftype = GTK_SIGNAL_LIST_ITEM_FACTORY(
			gtk_signal_list_item_factory_new());
		g_signal_connect(ftype, "setup", G_CALLBACK(type_setup_cb), NULL);
		g_signal_connect(ftype, "bind",  G_CALLBACK(type_bind_cb),  fb);
		GtkColumnViewColumn *ct = gtk_column_view_column_new(
			_("Type"), GTK_LIST_ITEM_FACTORY(ftype));
		gtk_column_view_column_set_resizable(ct, TRUE);
		gtk_column_view_column_set_sorter(ct,
			GTK_SORTER(gtk_custom_sorter_new(type_column_compare, NULL, NULL)));
		gtk_column_view_append_column(fb->columnview, ct);
		g_object_unref(ct);
	}
	{
		GtkSignalListItemFactory *fmod = GTK_SIGNAL_LIST_ITEM_FACTORY(
			gtk_signal_list_item_factory_new());
		g_signal_connect(fmod, "setup", G_CALLBACK(modified_setup_cb), NULL);
		g_signal_connect(fmod, "bind",  G_CALLBACK(modified_bind_cb),  fb);
		GtkColumnViewColumn *cm = gtk_column_view_column_new(
			_("Modified"), GTK_LIST_ITEM_FACTORY(fmod));
		gtk_column_view_column_set_resizable(cm, TRUE);
		gtk_column_view_column_set_sorter(cm,
			GTK_SORTER(gtk_custom_sorter_new(modified_column_compare, NULL, NULL)));
		gtk_column_view_append_column(fb->columnview, cm);
		g_object_unref(cm);
	}

	/* Wire the column view's sorter into the sort list model, behind a
	 * dirs-first primary sorter — so directories always group at the top,
	 * regardless of which column the user picks to sort by. */
	{
		GtkMultiSorter *multi = gtk_multi_sorter_new();
		gtk_multi_sorter_append(GTK_MULTI_SORTER(multi),
			GTK_SORTER(gtk_custom_sorter_new(dirs_first_compare, NULL, NULL)));
		GtkSorter *cv_sorter = gtk_column_view_get_sorter(fb->columnview);
		if (cv_sorter)
			gtk_multi_sorter_append(GTK_MULTI_SORTER(multi),
				g_object_ref(cv_sorter));
		fb->combined_sorter = GTK_SORTER(multi);
		gtk_sort_list_model_set_sorter(fb->sort_model, fb->combined_sorter);
	}

	/* Default sort: ascending name.  Clicking the same header again
	 * inverts the direction; clicking another header switches the key. */
	{
		GtkColumnViewColumn *name_col =
			g_list_model_get_item(G_LIST_MODEL(
				gtk_column_view_get_columns(fb->columnview)), 0);
		if (name_col) {
			gtk_column_view_sort_by_column(fb->columnview, name_col,
				GTK_SORT_ASCENDING);
			g_object_unref(name_col);
		}
	}

	gtk_column_view_set_single_click_activate(fb->columnview, FALSE);
	g_signal_connect(fb->columnview, "activate",
	                 G_CALLBACK(on_row_activated), fb);
	g_signal_connect(fb->selection, "selection-changed",
	                 G_CALLBACK(on_selection_changed), fb);

	GtkWidget *list_scrolled = gtk_scrolled_window_new();
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(list_scrolled),
	                              GTK_WIDGET(fb->columnview));
	gtk_widget_set_hexpand(list_scrolled, TRUE);
	gtk_widget_set_vexpand(list_scrolled, TRUE);
	/* File list keeps a usable minimum width so the columns don't get
	 * crushed when the dialog narrows; once we hit this floor, further
	 * shrinkage propagates outward (and eventually hits the window's own
	 * min size set below). */
	gtk_widget_set_size_request(list_scrolled, 320, -1);

	/* Preview pane.  The .siril-zero-pad class clears the global
	 * `box { padding: 1px }` for these specific boxes — that padding
	 * adds 2 px to the box's reported size_request, and GtkPaned was
	 * intermittently allocating the resulting outer min 1 px short of
	 * what the children needed (producing the `Trying to measure GtkBox
	 * for width N, but it needs at least N+1` warning flood during
	 * `gtk_window_present`'s size negotiation).  Zeroing the padding
	 * makes box.measure() return exactly size_request so the paned's
	 * allocation matches the children's needs. */
	GtkWidget *preview_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_add_css_class(preview_box, "siril-zero-pad");

	{
		GtkWidget *heading = make_pane_heading(_("Preview"));
		gtk_widget_set_margin_start(heading, 6);
		gtk_widget_set_margin_end(heading, 6);
		gtk_box_append(GTK_BOX(preview_box), heading);
	}
	fb->preview = GTK_PICTURE(gtk_picture_new());
	/* The preview area is exactly the thumbnail size the user picked in
	 * Preferences (128 / 256 / 512).  Using it as the size_request makes it
	 * both the floor — so the list|preview divider can be dragged in until
	 * the pane is just thumbnail_size wide, no wider — and the reference
	 * size every preview is rendered against (see the load paths, which cap
	 * thumbnails, FITS and raster alike, at thumbnail_size). */
	int preview_min = com.pref.gui.thumbnail_size;  /* 128, 256 or 512 */
	gtk_widget_set_size_request(GTK_WIDGET(fb->preview), preview_min, preview_min);
	/* The picture FILLS the preview pane (so its allocation can never
	 * exceed the window edge → it is never clipped) and uses SCALE_DOWN so
	 * the thumbnail is drawn at its native size, centred inside the
	 * allocation.  SCALE_DOWN only ever scales *down*, never up — and since
	 * every preview texture is itself capped at thumbnail_size, the image
	 * never renders larger than thumbnail_size even when the pane is dragged
	 * wider (it just gains centring margin). */
	gtk_picture_set_can_shrink(fb->preview, TRUE);
	gtk_picture_set_content_fit(fb->preview, GTK_CONTENT_FIT_SCALE_DOWN);
	gtk_widget_set_hexpand(GTK_WIDGET(fb->preview), TRUE);
	gtk_widget_set_vexpand(GTK_WIDGET(fb->preview), TRUE);
	gtk_widget_set_halign(GTK_WIDGET(fb->preview), GTK_ALIGN_FILL);
	gtk_widget_set_valign(GTK_WIDGET(fb->preview), GTK_ALIGN_FILL);
	gtk_box_append(GTK_BOX(preview_box), GTK_WIDGET(fb->preview));

	fb->metadata_label = GTK_LABEL(gtk_label_new(""));
	gtk_label_set_wrap(fb->metadata_label, TRUE);
	gtk_label_set_wrap_mode(fb->metadata_label, PANGO_WRAP_WORD_CHAR);
	gtk_label_set_xalign(fb->metadata_label, 0.5);
	gtk_label_set_justify(fb->metadata_label, GTK_JUSTIFY_CENTER);
	gtk_widget_set_halign(GTK_WIDGET(fb->metadata_label), GTK_ALIGN_CENTER);
	gtk_widget_set_margin_start(GTK_WIDGET(fb->metadata_label), 6);
	gtk_widget_set_margin_end(GTK_WIDGET(fb->metadata_label), 6);
	gtk_label_set_selectable(fb->metadata_label, TRUE);
	gtk_label_set_use_markup(fb->metadata_label, TRUE);
	gtk_box_append(GTK_BOX(preview_box), GTK_WIDGET(fb->metadata_label));

	GtkWidget *list_preview_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	fb->list_preview_paned = list_preview_paned;
	gtk_paned_set_start_child(GTK_PANED(list_preview_paned), list_scrolled);
	gtk_paned_set_end_child(GTK_PANED(list_preview_paned), preview_box);
	/* Restore the user's last list|preview divider position (persisted in
	 * the init file), falling back to the built-in default. */
	gtk_paned_set_position(GTK_PANED(list_preview_paned),
	                       com.pref.gui.open_dialog_paned_pos > 0
	                       ? com.pref.gui.open_dialog_paned_pos
	                       : LIST_PREVIEW_POSITION);
	gtk_paned_set_resize_start_child(GTK_PANED(list_preview_paned), FALSE);
	gtk_paned_set_resize_end_child(GTK_PANED(list_preview_paned), TRUE);
	gtk_paned_set_shrink_start_child(GTK_PANED(list_preview_paned), FALSE);
	gtk_paned_set_shrink_end_child(GTK_PANED(list_preview_paned), FALSE);

	/* Outer paned: sidebar | (list | preview).  shrink_start_child=FALSE
	 * keeps the sidebar from being shrunk below its 180 px size_request. */
	GtkWidget *outer_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	fb->outer_paned = outer_paned;
	gtk_paned_set_start_child(GTK_PANED(outer_paned), build_sidebar(fb));
	gtk_paned_set_end_child(GTK_PANED(outer_paned), list_preview_paned);
	gtk_paned_set_position(GTK_PANED(outer_paned),
	                       com.pref.gui.open_dialog_sidebar_pos > 0
	                       ? com.pref.gui.open_dialog_sidebar_pos
	                       : 180);
	gtk_paned_set_resize_start_child(GTK_PANED(outer_paned), FALSE);
	gtk_paned_set_shrink_start_child(GTK_PANED(outer_paned), FALSE);
	gtk_widget_set_vexpand(outer_paned, TRUE);

	/* Filter dropdown (initially empty; populated as filters are added). */
	fb->filter_combo = GTK_DROP_DOWN(gtk_drop_down_new_from_strings((const char *[]){NULL}));
	gtk_widget_set_visible(GTK_WIDGET(fb->filter_combo), FALSE);
	g_signal_connect(fb->filter_combo, "notify::selected",
	                 G_CALLBACK(on_filter_changed), fb);

	/* Debayer toggle — hidden unless siril_file_browser_set_show_debayer_toggle
	 * is called.  Initial state and sibling-button binding happen there. */
	fb->debayer_check = GTK_CHECK_BUTTON(gtk_check_button_new_with_label(_("Debayer")));
	gtk_widget_set_tooltip_text(GTK_WIDGET(fb->debayer_check),
		_("Debayer CFA images on open.  Linked to the same setting as the "
		  "Conversion tab's Debayer toggle."));
	gtk_widget_set_visible(GTK_WIDGET(fb->debayer_check), FALSE);

	/* Action row (bottom).  Cancel / Open have moved to the top toolbar, so
	 * the bottom strip now carries only the secondary controls: the
	 * Debayer toggle (when shown) sits at the left and the filter dropdown
	 * anchors at the right. */
	GtkWidget *action_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start(action_row, 6);
	gtk_widget_set_margin_end(action_row, 6);
	gtk_widget_set_margin_bottom(action_row, 6);
	gtk_widget_set_margin_top(action_row, 3);
	gtk_box_append(GTK_BOX(action_row), GTK_WIDGET(fb->debayer_check));
	GtkWidget *spacer = gtk_label_new(NULL);
	gtk_widget_set_hexpand(spacer, TRUE);
	gtk_box_append(GTK_BOX(action_row), spacer);
	gtk_box_append(GTK_BOX(action_row), GTK_WIDGET(fb->filter_combo));

	/* Assemble.  Search bar sits between toolbar and the paned content so
	 * its slide-down animation pushes the file list, not the sidebar. */
	GtkWidget *root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_box_append(GTK_BOX(root), toolbar);
	gtk_box_append(GTK_BOX(root), GTK_WIDGET(fb->search_bar));
	gtk_box_append(GTK_BOX(root), outer_paned);
	gtk_box_append(GTK_BOX(root), action_row);
	gtk_window_set_child(fb->window, root);

	g_signal_connect(fb->window, "close-request",
	                 G_CALLBACK(on_close_request), fb);

	update_nav_sensitivity(fb);
}

/* Clear all per-open state on the long-lived singleton.  Called on
 * every siril_file_browser_new() after the first to make the next open
 * indistinguishable from a fresh dialog.  Whatever the previous caller
 * configured (filters, initial location, multi-select, debayer toggle,
 * search, sort, selection, navigation history) is wiped here. */
static void reset_browser_state(SirilFileBrowser *fb) {
	if (!fb) return;

	/* Cancel any in-flight breadcrumb-rebuild idle from the previous run. */
	if (fb->pending_breadcrumb_idle) {
		g_source_remove(fb->pending_breadcrumb_idle);
		fb->pending_breadcrumb_idle = 0;
	}

	/* Drop the per-open external-demosaic toggle wiring.  The Convert tab
	 * button is owned by the main UI builder so we just disconnect our
	 * handler; the next caller that wants it will rewire via
	 * siril_file_browser_set_show_debayer_toggle. */
	if (fb->external_demosaic_btn && fb->external_demosaic_handler) {
		g_signal_handler_disconnect(fb->external_demosaic_btn,
		                            fb->external_demosaic_handler);
		fb->external_demosaic_handler = 0;
	}
	fb->external_demosaic_btn = NULL;
	fb->show_debayer_toggle = FALSE;
	if (fb->debayer_check)
		gtk_widget_set_visible(GTK_WIDGET(fb->debayer_check), FALSE);

	/* Folder-picker mode is per-open; clear it so the next caller starts as
	 * a normal file picker unless it opts back in. */
	fb->directory_only = FALSE;
	if (fb->new_folder_btn)
		gtk_widget_set_visible(fb->new_folder_btn, FALSE);

	/* Revert multi-select back to single — set_select_multiple() rebuilds
	 * the selection model and rewires the columnview for us. */
	if (fb->select_multiple)
		siril_file_browser_set_select_multiple(fb, FALSE);
	else if (fb->selection)
		gtk_selection_model_unselect_all(fb->selection);

	/* Filters: drop all, hide the dropdown until a new caller adds some. */
	if (fb->filters) g_ptr_array_set_size(fb->filters, 0);
	fb->active_filter_idx = -1;
	if (fb->filter_combo) {
		gtk_drop_down_set_model(fb->filter_combo, NULL);
		gtk_widget_set_visible(GTK_WIDGET(fb->filter_combo), FALSE);
	}
	if (fb->file_filter)
		gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);

	/* Result + response. */
	g_clear_pointer(&fb->result_path, g_free);
	g_slist_free_full(fb->result_paths, g_free);
	fb->result_paths = NULL;
	fb->response = GTK_RESPONSE_NONE;

	/* Navigation: drop current folder, history, breadcrumb trail.
	 * The caller will set initial_folder/initial_file before run(),
	 * or run() itself will fall back to com.wd. */
	g_clear_object(&fb->current_folder);
	g_clear_object(&fb->initial_file);
	g_clear_object(&fb->breadcrumb_deepest);
	if (fb->back_history)    g_ptr_array_set_size(fb->back_history, 0);
	if (fb->forward_history) g_ptr_array_set_size(fb->forward_history, 0);

	/* Search: clear text and collapse the bar. */
	g_clear_pointer(&fb->search_text, g_free);
	if (fb->search_entry)
		gtk_editable_set_text(GTK_EDITABLE(fb->search_entry), "");
	if (fb->search_bar)
		gtk_search_bar_set_search_mode(fb->search_bar, FALSE);

	/* Recent-files list: drop the cached store (will be rebuilt on demand). */
	g_clear_object(&fb->recent_store);
	fb->in_recent_mode = FALSE;

	/* Pref-driven settings: re-read in case the user changed them between
	 * dialog runs. */
	fb->show_hidden_files = browser_initial_show_hidden();

	/* Open button starts insensitive until a selection lands. */
	if (fb->open_button)
		gtk_widget_set_sensitive(fb->open_button, FALSE);

	update_nav_sensitivity(fb);
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
			G_LIST_MODEL(g_object_ref(fb->sort_model))));
	} else {
		GtkSingleSelection *ss = gtk_single_selection_new(
			G_LIST_MODEL(g_object_ref(fb->sort_model)));
		gtk_single_selection_set_autoselect(ss, FALSE);
		gtk_single_selection_set_can_unselect(ss, TRUE);
		fresh = GTK_SELECTION_MODEL(ss);
	}
	g_clear_object(&fb->selection);
	fb->selection = fresh;  /* takes the +1 ref we got from _new */
	gtk_column_view_set_model(fb->columnview, fb->selection);
	g_signal_connect(fb->selection, "selection-changed",
	                 G_CALLBACK(on_selection_changed), fb);
}

void siril_file_browser_set_directory_only(SirilFileBrowser *fb, gboolean dir_only) {
	if (!fb) return;
	fb->directory_only = dir_only;
	/* Re-run the filter so any already-bound rows rebind and pick up the
	 * greyed / non-selectable state (no-op for a not-yet-shown dialog). */
	if (fb->file_filter)
		gtk_filter_changed(GTK_FILTER(fb->file_filter), GTK_FILTER_CHANGE_DIFFERENT);
	/* Open targets the current folder in this mode, so it's usable from the
	 * start; in normal mode it stays gated on an actual selection. */
	if (fb->open_button)
		gtk_widget_set_sensitive(fb->open_button, dir_only || any_selected(fb));
	/* The create-folder button only makes sense when picking a directory. */
	if (fb->new_folder_btn)
		gtk_widget_set_visible(fb->new_folder_btn, dir_only);
}

gint siril_file_browser_run(SirilFileBrowser *fb) {
	if (!fb) return GTK_RESPONSE_NONE;
	if (!fb->current_folder) {
		/* Default to Siril's working directory — that's where the user's
		 * current project lives, so opening a new file is almost always
		 * relative to it.  Fall back to $HOME only when wd is unset (e.g.
		 * the first dialog before any "cd" command has been issued). */
		if (com.wd && *com.wd && g_file_test(com.wd, G_FILE_TEST_IS_DIR)) {
			siril_file_browser_set_initial_folder(fb, com.wd);
		} else {
			const gchar *home = g_get_home_dir();
			if (home) siril_file_browser_set_initial_folder(fb, home);
		}
	}
	fb->response = GTK_RESPONSE_NONE;
	g_clear_pointer(&fb->result_path, g_free);
	g_slist_free_full(fb->result_paths, g_free);
	fb->result_paths = NULL;

	/* Force every row to re-bind so the per-item directory-only state
	 * (files greyed + non-selectable) reflects THIS run.  The browser is a
	 * process-wide singleton that recycles its row widgets, so reshowing
	 * the same folder would otherwise reuse rows still carrying the
	 * previous run's state — leaving files greyed after a folder pick, or
	 * un-greyed when switching back into folder mode.  Detaching and
	 * re-attaching the model drops all realized rows and rebuilds them
	 * (setup + bind) against the current fb->directory_only.  Park focus
	 * first: disposing focused rows mid-walk trips a stale-widget assert
	 * (same precaution as apply_current_folder). */
	if (fb->columnview && fb->selection) {
		if (fb->window)
			gtk_window_set_focus(fb->window, NULL);
		gtk_column_view_set_model(fb->columnview, NULL);
		gtk_column_view_set_model(fb->columnview, fb->selection);
	}

	fb->loop = g_main_loop_new(NULL, FALSE);
	gtk_window_present(fb->window);
	g_main_loop_run(fb->loop);
	g_main_loop_unref(fb->loop);
	fb->loop = NULL;
	/* Persist the current size and divider positions before hiding (the
	 * allocation is still valid here; once hidden, gtk_widget_get_width/
	 * height read 0).  Only write the init file when something actually
	 * changed to avoid needless disk churn on every open. */
	if (com.pref.gui.remember_windows) {
		int w = gtk_widget_get_width(GTK_WIDGET(fb->window));
		int h = gtk_widget_get_height(GTK_WIDGET(fb->window));
		int sidebar = fb->outer_paned
			? gtk_paned_get_position(GTK_PANED(fb->outer_paned)) : 0;
		int divider = fb->list_preview_paned
			? gtk_paned_get_position(GTK_PANED(fb->list_preview_paned)) : 0;
		gboolean changed = FALSE;
		if (w > 0 && h > 0 &&
		    (w != com.pref.gui.open_dialog_w || h != com.pref.gui.open_dialog_h)) {
			com.pref.gui.open_dialog_w = w;
			com.pref.gui.open_dialog_h = h;
			changed = TRUE;
		}
		if (sidebar > 0 && sidebar != com.pref.gui.open_dialog_sidebar_pos) {
			com.pref.gui.open_dialog_sidebar_pos = sidebar;
			changed = TRUE;
		}
		if (divider > 0 && divider != com.pref.gui.open_dialog_paned_pos) {
			com.pref.gui.open_dialog_paned_pos = divider;
			changed = TRUE;
		}
		if (changed)
			writeinitfile();
	}
	gtk_widget_set_visible(GTK_WIDGET(fb->window), FALSE);
	/* Re-enable the parent we disabled in _new (Plan C: substitute for
	 * the modal flag).  Done after hiding so the visual transition is
	 * "dialog disappears → parent regains input".  parent is cleared so
	 * the next _new call can rebind to whatever new parent it gets. */
	if (fb->parent) {
		gtk_widget_set_sensitive(GTK_WIDGET(fb->parent), TRUE);
		fb->parent = NULL;
	}
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

/* siril_file_browser_destroy was removed in the singleton conversion —
 * the browser is now a long-lived process-singleton (see g_singleton
 * above), and per-open state is wiped by reset_browser_state() on each
 * subsequent _new() call.  Callers no longer pair _new with _destroy. */
