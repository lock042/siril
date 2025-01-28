/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "core/proto.h"
#include "io/annotation_catalogues.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"

#include "annotations_pref.h"

static GtkListStore *list_store_catalogue = NULL;

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
		N_("User Solar System Objects Catalogue")
};
// update the size of gui_config.catalog if changed

static void get_astrolist_store() {
	if (list_store_catalogue == NULL)
		list_store_catalogue = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststore_astrometry"));
}

void fill_astrometry_catalogue(gboolean *catalog) {
	GtkTreeIter iter;
	get_astrolist_store();
	gtk_list_store_clear(list_store_catalogue);

	for(int i = 0; i < G_N_ELEMENTS(astro_catalogue); i++) {
		gtk_list_store_append(list_store_catalogue, &iter);
		gtk_list_store_set(list_store_catalogue, &iter,
				0, catalog[i],
				1, _(astro_catalogue[i]),
				-1);
	}
}

gboolean on_annotate_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	if (event->button == GDK_BUTTON_SECONDARY) {
#if GTK_CHECK_VERSION(3, 22, 0)
		gtk_popover_popup(GTK_POPOVER(gtk_builder_get_object(gui.builder, "popover_catalogs")));
#else
		gtk_widget_show(gtk_builder_get_object(gui.builder, "popover_catalogs"));
#endif
	}


	return FALSE;
}

void get_astrometry_catalogue_values() {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store_catalogue), &iter);

	int i = 0;
	while (valid) {
		gboolean value = TRUE;
		gtk_tree_model_get(GTK_TREE_MODEL(list_store_catalogue), &iter, 0, &value, -1);
		com.pref.gui.catalog[i] = value;
		valid = gtk_tree_model_iter_next (GTK_TREE_MODEL(list_store_catalogue), &iter);
		i++;
	}
}

void on_cellrendeur_catalog_use_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *path, gpointer user_data) {
	GtkTreeIter iter;
	gboolean value;

	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(list_store_catalogue), &iter, path);
	gtk_tree_model_get(GTK_TREE_MODEL(list_store_catalogue), &iter, 0, &value, -1);
	gtk_list_store_set(list_store_catalogue, &iter, 0, !value, -1);
}

void on_cellrendeur2_catalog_use_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *path, gpointer user_data) {
    GtkTreeIter iter;
    gboolean value;
    GtkTreePath *tree_path;
    gint *indices;

    gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(list_store_catalogue), &iter, path);
    gtk_tree_model_get(GTK_TREE_MODEL(list_store_catalogue), &iter, 0, &value, -1);
    gtk_list_store_set(list_store_catalogue, &iter, 0, !value, -1);

    tree_path = gtk_tree_model_get_path(GTK_TREE_MODEL(list_store_catalogue), &iter);
    indices = gtk_tree_path_get_indices(tree_path);
    if (indices) {
        com.pref.gui.catalog[indices[0]] = !value;
		if (has_wcs(&gfit) && gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "annotate_button")))) {
			refresh_annotation_visibility();
			refresh_found_objects();
			redraw(REDRAW_OVERLAY);
		}
    }

    gtk_tree_path_free(tree_path);
}

