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

#include "compositing/hdr.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/file_browser.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/open_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "io/image_format_fits.h"

static GtkListBox *hdr_list = NULL;
static GtkSpinButton *hdr_feather = NULL;


static double row_exposure(GtkListBoxRow *row) {
	double *exp = g_object_get_data(G_OBJECT(row), "exposure");
	return exp ? *exp : 0.0;
}

static int sort_by_exposure(GtkListBoxRow *a, GtkListBoxRow *b, gpointer user_data) {
	double x = row_exposure(a), y = row_exposure(b);
	return (x < y) - (x > y);
}

void hdr_dialog_init_statics(void) {
	if (hdr_list) return;
	hdr_list = GTK_LIST_BOX(gtk_builder_get_object(gui.builder, "hdr_file_list"));
	hdr_feather = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "hdr_feather"));
	gtk_list_box_set_sort_func(hdr_list, sort_by_exposure, NULL, NULL);
}

static gboolean list_has_path(const gchar *path) {
	for (GtkWidget *w = gtk_widget_get_first_child(GTK_WIDGET(hdr_list)); w; w = gtk_widget_get_next_sibling(w)) {
		const gchar *p = g_object_get_data(G_OBJECT(w), "path");
		if (p && !g_strcmp0(p, path))
			return TRUE;
	}
	return FALSE;
}

static void add_row(const gchar *path, fits *f) {
	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
	gtk_widget_set_margin_start(box, 6);
	gtk_widget_set_margin_end(box, 6);
	gtk_widget_set_margin_top(box, 3);
	gtk_widget_set_margin_bottom(box, 3);

	gchar *base = g_path_get_basename(path);
	GtkWidget *name = gtk_label_new(base);
	g_free(base);
	gtk_label_set_ellipsize(GTK_LABEL(name), PANGO_ELLIPSIZE_MIDDLE);
	gtk_label_set_xalign(GTK_LABEL(name), 0.f);
	gtk_widget_set_hexpand(name, TRUE);
	gtk_widget_set_tooltip_text(name, path);
	gtk_box_append(GTK_BOX(box), name);

	gchar *size = g_strdup_printf("%ld×%ld, %s", f->naxes[0], f->naxes[1],
			f->naxes[2] == 3 ? _("colour") : _("mono"));
	GtkWidget *dims = gtk_label_new(size);
	g_free(size);
	gtk_widget_add_css_class(dims, "dim-label");
	gtk_box_append(GTK_BOX(box), dims);

	double exposure = f->keywords.exposure > 0.0 ? f->keywords.exposure : 0.0;
	gchar *exptxt = exposure > 0.0 ? g_strdup_printf("%g s", exposure) : g_strdup(_("no exposure"));
	GtkWidget *exp = gtk_label_new(exptxt);
	g_free(exptxt);
	gtk_label_set_width_chars(GTK_LABEL(exp), 10);
	gtk_label_set_xalign(GTK_LABEL(exp), 1.f);
	if (exposure <= 0.0)
		gtk_widget_set_tooltip_text(exp, _("The exposure time is missing from the header: "
					"the exposure ratios will be measured on the images"));
	gtk_box_append(GTK_BOX(box), exp);

	GtkWidget *row = gtk_list_box_row_new();
	gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), box);
	g_object_set_data_full(G_OBJECT(row), "path", g_strdup(path), g_free);
	double *pexp = g_new(double, 1);
	*pexp = exposure;
	g_object_set_data_full(G_OBJECT(row), "exposure", pexp, g_free);
	gtk_list_box_append(hdr_list, row);
}

void on_hdr_add_clicked(GtkButton *button, gpointer user_data) {
	hdr_dialog_init_statics();
	GtkWindow *parent = GTK_WINDOW(gtk_builder_get_object(gui.builder, "hdr_dialog"));
	SirilFileBrowser *fb = siril_file_browser_new(parent, _("Add images"));
	siril_file_browser_set_select_multiple(fb, TRUE);
	if (com.wd && *com.wd)
		siril_file_browser_set_initial_folder(fb, com.wd);
	siril_file_browser_add_filter_pattern(fb,
			_("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"),
			FITS_EXTENSIONS, TRUE);
	GSList *paths = NULL;
	if (siril_file_browser_run(fb) == GTK_RESPONSE_ACCEPT)
		paths = siril_file_browser_get_paths(fb);

	for (GSList *l = paths; l; l = l->next) {
		const gchar *path = l->data;
		if (!path || list_has_path(path))
			continue;
		fits f = { 0 };
		if (read_fits_metadata_from_path(path, &f)) {
			siril_log_error(_("Could not open file: %s\n"), path);
			continue;
		}
		add_row(path, &f);
		clearfits(&f);
	}
	g_slist_free_full(paths, g_free);
	gtk_list_box_invalidate_sort(hdr_list);
}

void on_hdr_remove_clicked(GtkButton *button, gpointer user_data) {
	hdr_dialog_init_statics();
	GList *rows = gtk_list_box_get_selected_rows(hdr_list);
	for (GList *l = rows; l; l = l->next)
		gtk_list_box_remove(hdr_list, GTK_WIDGET(l->data));
	g_list_free(rows);
}

void on_hdr_clear_clicked(GtkButton *button, gpointer user_data) {
	hdr_dialog_init_statics();
	gtk_list_box_remove_all(hdr_list);
}

void on_hdr_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("hdr_dialog");
}

gboolean hdr_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("hdr_dialog");
	return TRUE;
}

void on_hdr_apply_clicked(GtkButton *button, gpointer user_data) {
	hdr_dialog_init_statics();
	struct hdr_args *args = calloc(1, sizeof(struct hdr_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	for (GtkWidget *w = gtk_widget_get_first_child(GTK_WIDGET(hdr_list)); w; w = gtk_widget_get_next_sibling(w)) {
		const gchar *p = g_object_get_data(G_OBJECT(w), "path");
		if (p)
			args->files = g_slist_append(args->files, g_strdup(p));
	}
	if (g_slist_length(args->files) < 2) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Not enough images"),
				_("HDR composition needs at least two images of different exposures."));
		free_hdr_args(args);
		return;
	}
	hdr_params_init(&args->params);
	args->params.feather = gtk_spin_button_get_value_as_int(hdr_feather);

	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_close_dialog("hdr_dialog");
	if (!start_in_new_thread(hdr_worker, args))
		free_hdr_args(args);
}
