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
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/utils.h"
#include "core/proto.h"
#include "io/annotation_catalogues.h"
#include "algos/siril_wcs.h"
#include "gui-gtk4/image_display.h"
#include <gtk/gtk.h>

#include "annotations_pref.h"

static gchar *astro_catalogue[] = {
		N_("Messier Catalogue (M)"),
		N_("New General Catalogue (NGC)"),
		N_("Index Catalogue (IC)"),
		N_("Lynds Catalogue of Dark Nebulae (LdN)"),
		N_("Sharpless Catalogue (Sh2)"),
		N_("Star Catalogue"),
		N_("IAU Constellations Catalogue"),
		N_("IAU Constellations Names"),
		N_("User Deep Sky Objects Catalogue"),
		N_("User Solar System Objects Catalogue"),
		N_("Velocity vectors for Solar System Objects")
};
// update the size of gui_config.catalog if changed

/* ---- Row GObject ------------------------------------------------- */

#define SIRIL_TYPE_CATALOGUE_ROW (siril_catalogue_row_get_type())
G_DECLARE_FINAL_TYPE(SirilCatalogueRow, siril_catalogue_row, SIRIL, CATALOGUE_ROW, GObject)

struct _SirilCatalogueRow {
	GObject parent_instance;
	gboolean use_flag;
	gchar   *name;
	guint    idx;            /* index into com.pref.gui.catalog[] */
	gboolean trigger_refresh; /* whether toggle should call refresh_annotation_visibility */
};

G_DEFINE_TYPE(SirilCatalogueRow, siril_catalogue_row, G_TYPE_OBJECT)

static void siril_catalogue_row_finalize(GObject *obj) {
	SirilCatalogueRow *self = SIRIL_CATALOGUE_ROW(obj);
	g_clear_pointer(&self->name, g_free);
	G_OBJECT_CLASS(siril_catalogue_row_parent_class)->finalize(obj);
}
static void siril_catalogue_row_class_init(SirilCatalogueRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_catalogue_row_finalize;
}
static void siril_catalogue_row_init(SirilCatalogueRow *self) { (void)self; }

/* Single shared store; both the settings-dialog view and the popover view
 * point at the same GListStore so toggling one updates the other. */
static GListStore *catalogue_store = NULL;

static void ensure_catalogue_store(gboolean *catalog) {
	if (catalogue_store) return;
	catalogue_store = g_list_store_new(SIRIL_TYPE_CATALOGUE_ROW);
	for (guint i = 0; i < G_N_ELEMENTS(astro_catalogue); i++) {
		SirilCatalogueRow *row = g_object_new(SIRIL_TYPE_CATALOGUE_ROW, NULL);
		row->use_flag = catalog ? catalog[i] : FALSE;
		row->name     = g_strdup(_(astro_catalogue[i]));
		row->idx      = i;
		row->trigger_refresh = FALSE;
		g_list_store_append(catalogue_store, row);
		g_object_unref(row);
	}
}

/* Forward decl: defined further below; needed by get_astrometry_catalogue_values. */
static void on_annotate_secondary_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);

/* ---- Cell factories --------------------------------------------- */

static void on_use_toggled(GtkCheckButton *btn, gpointer user_data) {
	(void)user_data;
	GtkListItem *li = g_object_get_data(G_OBJECT(btn), "siril-listitem");
	if (!li) return;
	SirilCatalogueRow *row = SIRIL_CATALOGUE_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gboolean new_val = siril_toggle_get_active(GTK_WIDGET(btn));
	row->use_flag = new_val;
	if (row->idx < G_N_ELEMENTS(astro_catalogue))
		com.pref.gui.catalog[row->idx] = new_val;
	if (row->trigger_refresh && has_wcs(gfit) &&
			siril_toggle_get_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "annotate_button")))) {
		refresh_annotation_visibility();
		refresh_found_objects();
		redraw(REDRAW_OVERLAY);
	}
}

static void use_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *btn = gtk_check_button_new();
	gtk_widget_set_halign(btn, GTK_ALIGN_CENTER);
	gtk_list_item_set_child(li, btn);
	g_signal_connect(btn, "toggled", G_CALLBACK(on_use_toggled), NULL);
}

static void use_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	gboolean trigger = (gboolean)GPOINTER_TO_INT(u);
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	SirilCatalogueRow *row = SIRIL_CATALOGUE_ROW(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", li);
	row->trigger_refresh = trigger;
	g_signal_handlers_block_by_func(btn, on_use_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(btn), row->use_flag);
	g_signal_handlers_unblock_by_func(btn, on_use_toggled, NULL);
}

static void use_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", NULL);
}

static void name_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}

static void name_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilCatalogueRow *row = SIRIL_CATALOGUE_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row->name);
}

static GtkWidget *build_catalogue_columnview(gboolean trigger_refresh) {
	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(catalogue_store)));
	GtkColumnView *cv = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));
	gtk_widget_add_css_class(GTK_WIDGET(cv), "siril-dense-rows");

	GtkSignalListItemFactory *fuse = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fuse, "setup",  G_CALLBACK(use_setup_cb),  NULL);
	g_signal_connect(fuse, "bind",   G_CALLBACK(use_bind_cb),   GINT_TO_POINTER(trigger_refresh));
	g_signal_connect(fuse, "unbind", G_CALLBACK(use_unbind_cb), NULL);

	GtkSignalListItemFactory *fname = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fname, "setup", G_CALLBACK(name_setup_cb), NULL);
	g_signal_connect(fname, "bind",  G_CALLBACK(name_bind_cb),  NULL);

	GtkColumnViewColumn *cu = gtk_column_view_column_new("", GTK_LIST_ITEM_FACTORY(fuse));
	gtk_column_view_column_set_resizable(cu, FALSE);
	gtk_column_view_append_column(cv, cu);
	g_object_unref(cu);

	GtkColumnViewColumn *cn = gtk_column_view_column_new(N_("Catalogue"), GTK_LIST_ITEM_FACTORY(fname));
	gtk_column_view_column_set_resizable(cn, TRUE);
	gtk_column_view_column_set_expand(cn, TRUE);
	gtk_column_view_append_column(cv, cn);
	g_object_unref(cn);

	/* No g_object_unref(sel): gtk_column_view_new() is transfer-full for
	 * the model, so cv consumed the ref returned by single_selection_new. */
	return GTK_WIDGET(cv);
}

void fill_astrometry_catalogue(gboolean *catalog) {
	ensure_catalogue_store(catalog);
	/* Sync store to current values (in case prefs changed after first build). */
	guint n = g_list_model_get_n_items(G_LIST_MODEL(catalogue_store));
	for (guint i = 0; i < n; i++) {
		SirilCatalogueRow *row = SIRIL_CATALOGUE_ROW(g_list_model_get_item(G_LIST_MODEL(catalogue_store), i));
		if (row && i < G_N_ELEMENTS(astro_catalogue))
			row->use_flag = catalog[i];
		g_object_unref(row);
	}

	/* Inject the column view into each scrolled window placeholder if not
	 * already populated.  The settings_window's "scrolled_astro_catalogue"
	 * and the main popover's "scrolled_popover_catalogue" hold these. */
	GtkScrolledWindow *sw1 = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolled_astro_catalogue"));
	if (sw1 && !gtk_scrolled_window_get_child(sw1))
		gtk_scrolled_window_set_child(sw1, build_catalogue_columnview(FALSE));
	GtkScrolledWindow *sw2 = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolled_popover_catalogue"));
	if (sw2 && !gtk_scrolled_window_get_child(sw2))
		gtk_scrolled_window_set_child(sw2, build_catalogue_columnview(TRUE));
	gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(sw2), TRUE);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(sw2), TRUE);

	/* GTK4: secondary-click on the annotate toggle button opens the
	 * catalog selector popover.  Replaces the GTK3 .ui-bound
	 * "button-press-event" handler. */
	GtkWidget *annotate_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "annotate_button"));
	if (annotate_button && !g_object_get_data(G_OBJECT(annotate_button), "siril-rmb-controller")) {
		GtkGesture *click = gtk_gesture_click_new();
		gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click),
		                              GDK_BUTTON_SECONDARY);
		g_signal_connect(click, "pressed",
		                 G_CALLBACK(on_annotate_secondary_pressed), NULL);
		gtk_widget_add_controller(annotate_button, GTK_EVENT_CONTROLLER(click));
		g_object_set_data(G_OBJECT(annotate_button), "siril-rmb-controller", click);
	}
}

/* GTK4 secondary-click handler for the annotate toggle button: opens the
 * catalog-selection popover.  The popover itself is a GtkPopover defined
 * in the .ui as `popover_catalogs`; here we just point it at the button
 * and pop it up. */
static void on_annotate_secondary_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)n_press; (void)x; (void)y; (void)user_data;
	GtkWidget *popover = GTK_WIDGET(gtk_builder_get_object(gui.builder, "popover_catalogs"));
	if (!popover) return;
	/* The popover lives in the .ui as a top-level object; re-parent it to
	 * the annotate button (its visual anchor) the first time we pop it. */
	GtkWidget *anchor = gtk_event_controller_get_widget(
	    GTK_EVENT_CONTROLLER(gesture));
	GtkWidget *cur_parent = gtk_widget_get_parent(popover);
	if (cur_parent != anchor) {
		if (cur_parent)
			gtk_widget_unparent(popover);
		gtk_widget_set_parent(popover, anchor);
	}
	gtk_popover_popup(GTK_POPOVER(popover));
}


void get_astrometry_catalogue_values() {
	if (!catalogue_store) return;
	guint n = g_list_model_get_n_items(G_LIST_MODEL(catalogue_store));
	for (guint i = 0; i < n && i < G_N_ELEMENTS(astro_catalogue); i++) {
		SirilCatalogueRow *row = SIRIL_CATALOGUE_ROW(g_list_model_get_item(G_LIST_MODEL(catalogue_store), i));
		if (row) com.pref.gui.catalog[i] = row->use_flag;
		g_object_unref(row);
	}
}

/* The .ui-wired GtkCellRendererToggle handlers (signature
 * GtkCellRendererToggle*, gchar*, gpointer) are obsolete now that the
 * cell renderers are gone.  The .ui files no longer reference them, but
 * keep stubs so any stale .ui or partial-build cache that still wires the
 * handler links cleanly. */
void on_cellrendeur_catalog_use_toggled(GtkWidget *cell_renderer,
		gchar *path, gpointer user_data) {
	(void)cell_renderer; (void)path; (void)user_data;
}

void on_cellrendeur2_catalog_use_toggled(GtkWidget *cell_renderer,
		gchar *path, gpointer user_data) {
	(void)cell_renderer; (void)path; (void)user_data;
}
