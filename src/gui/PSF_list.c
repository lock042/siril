/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/PSF_list.h"
#include "gui/progress_and_log.h"
#include "algos/siril_wcs.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/ccd-inspector.h"
#include "algos/astrometry_solver.h"

static GtkListStore *liststore_stars = NULL;

enum {
	COLUMN_CHANNEL,		// int
	COLUMN_B,			// gdouble
	COLUMN_A,			// gdouble
	COLUMN_X0,			// gdouble
	COLUMN_Y0,			// gdouble
	COLUMN_FWHMX,		// gdouble
	COLUMN_FWHMY,		// gdouble
	COLUMN_MAG,		    // gdouble
	COLUMN_BETA,		// gdouble
	COLUMN_ROUNDNESS,	// gdouble
	COLUMN_ANGLE,		// gdouble
	COLUMN_RMSE,		// gdouble
	COLUMN_INDEX,       // gint
	N_COLUMNS
};

enum {					//different context_id of the GtkStatusBar
	COUNT_STATE
};

static gchar *units = "";

static void gdouble_fwhmx_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_FWHMX, &var, -1);
	buf = g_strdup_printf("%.2f%s", var, units);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_fwhmy_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_FWHMY, &var, -1);
	buf = g_strdup_printf("%.2f%s", var, units);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_x0_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_X0, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_y0_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_Y0, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_mag_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_MAG, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_r_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_ROUNDNESS, &var, -1);
	buf = g_strdup_printf("%.3f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_angle_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_ANGLE, &var, -1);
	if (var == 0.0)
		buf = g_strdup_printf("%s", "N/A");
	else
		buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_beta_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_BETA, &var, -1);
	if (var == -1.0)
		buf = g_strdup_printf("%s", "N/A");
	else
		buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_rmse_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_RMSE, &var, -1);
	buf = g_strdup_printf("%.2e", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void get_stars_list_store() {
	if (liststore_stars == NULL)
		liststore_stars = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststore_stars"));

	GtkTreeViewColumn *col;
	GtkCellRenderer *cell;

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn7"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_x0"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_x0_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn8"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_y0"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_y0_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn9"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_fwhmx"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_fwhmx_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn10"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_fwhmy"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_fwhmy_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn_mag"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_mag"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_mag_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn29"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_beta"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_beta_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn14"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_r"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_r_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn6"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_angle"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_angle_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn15"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_rmse"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_rmse_cell_data_function, NULL, NULL);
}

static void display_PSF(psf_star **result) {
	if (result) {
		gchar *msg;
		int i = 0;
		double FWHMx = 0.0, FWHMy = 0.0, B = 0.0, A = 0.0, beta = 0.0, r = 0.0, angle = 0.0,
				rmse = 0.0;
		gboolean unit_is_arcsec;
		int n = 0, layer;
		starprofile profiletype;

		while (result[i]) {
			double fwhmx, fwhmy;
			char *unit;
			gboolean is_as = get_fwhm_as_arcsec_if_possible(result[i], &fwhmx, &fwhmy, &unit);
			if (i == 0) {
				unit_is_arcsec = is_as;
				layer = result[i]->layer;
				profiletype = result[i]->profile;
			}
			else if (is_as != unit_is_arcsec) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars FWHM must have the same units."));
				return;
			}
			else if (layer != result[i]->layer ) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars properties must all be computed on the same layer"));
				return;
			}
			else if (profiletype != result[i]->profile) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars must all be modeled with the same profile type"));
				return;
			}
			if (!result[i]->has_saturated) {
				B += result[i]->B;
				A += result[i]->A;
				FWHMx += fwhmx;
				beta += result[i]->beta;

				FWHMy += fwhmy;
				angle += result[i]->angle;
				rmse += result[i]->rmse;
				n++;
				}
			i++;
		}
		if (i <= 0 || n <= 0) return;
		/* compute average */
		B = B / (double)n;
		A = A / (double)n;
		FWHMx = FWHMx / (double)n;
		FWHMy = FWHMy / (double)n;
		beta =beta / (double)n;
		r = FWHMy / FWHMx;
		angle = angle / (double)n;
		rmse = rmse / (double)n;

		if (profiletype == PSF_GAUSSIAN) {
		msg = g_strdup_printf(_("Average Gaussian PSF\n\n"
				"N:\t%d stars (%d saturated and excluded)\nB:\t%.6f\nA:\t%.6f\nFWHMx:\t%.2f%s\n"
				"FWHMy:\t%.2f%s\nr:\t%.3f\nAngle:\t%.2f deg\nrmse:\t%.3e\n"),
				i, i - n, B, A, FWHMx, result[0]->units, FWHMy,
				result[0]->units, r, angle, rmse);
		} else {
		msg = g_strdup_printf(_("Average Moffat PSF\n\n"
				"N:\t%d stars (%d saturated and excluded)\nB:\t%.6f\nA:\t%.6f\nBeta:\t%.2f\nFWHMx:\t%.2f%s\n"
				"FWHMy:\t%.2f%s\nr:\t%.3f\nAngle:\t%.2f deg\nrmse:\t%.3e\n"),
				i, i - n, B, A, beta, FWHMx, result[0]->units, FWHMy,
				result[0]->units, r, angle, rmse);
		}
		show_data_dialog(msg, _("Average Star Data"), "stars_list_window", NULL);
		g_free(msg);
	}
}


static gint get_index_of_selected_star(gdouble x, gdouble y) {
	int i = 0;

	while (com.stars && com.stars[i]) {
		if ((com.stars[i]->xpos == x) && (com.stars[i]->ypos == y)) {
			return i;
		}
		i++;
	}
	return -1;
}

static void display_status() {
	gchar *text;
	int i = 0;
	GtkStatusbar *statusbar;

	statusbar = GTK_STATUSBAR(lookup_widget("statusbar_PSF"));

	while (com.stars && com.stars[i])
		i++;
	if (gui.selected_star == -1) {
		if (i > 0) {
			text = ngettext("%d star", "%d stars", i);
			text = g_strdup_printf(text, i);
		} else {
			text = g_strdup(" ");
		}
	} else {
		text = g_strdup_printf(_("Star %d of %d"), gui.selected_star + 1, i);
	}
	gtk_statusbar_push(statusbar, COUNT_STATE, text);
	g_free(text);
}

void set_iter_of_clicked_psf(double x, double y) {
	GtkTreeView *treeview = GTK_TREE_VIEW(lookup_widget("Stars_stored"));
	GtkTreeModel *model = gtk_tree_view_get_model(treeview);
	GtkTreeIter iter;
	gboolean valid;
	gboolean is_as;
	const double radian_conversion = ((3600.0 * 180.0) / M_PI) / 1.0E3;
	double invpixscalex = 1.0;
	double bin_X = com.pref.binning_update ? (double) gfit.binning_x : 1.0;
	if (com.stars && com.stars[0]) {// If the first star has units of arcsec, all should have
		is_as = (strcmp(com.stars[0]->units,"px"));
	} else {
		return; // If com.stars is empty there is no point carrying on
	}
	if (is_as) {
		invpixscalex = 1.0 / (radian_conversion * (double) gfit.pixel_size_x / gfit.focal_length) * bin_X;
	}
	valid = gtk_tree_model_get_iter_first(model, &iter);
	while (valid) {
		gdouble xpos, ypos, fwhmx;
		gtk_tree_model_get(model, &iter, COLUMN_X0, &xpos, -1);
		gtk_tree_model_get(model, &iter, COLUMN_Y0, &ypos, -1);
		gtk_tree_model_get(model, &iter, COLUMN_FWHMX, &fwhmx, -1);
		fwhmx *= invpixscalex;
		gdouble distsq = (xpos - x) * (xpos - x) + (ypos - y) * (ypos - y);
		gdouble psflimsq = 6. * fwhmx * fwhmx;
		if (distsq < psflimsq) {
			GtkTreeSelection *selection = gtk_tree_view_get_selection(treeview);
			GtkTreePath *path = gtk_tree_model_get_path(model, &iter);
			if (!path) return;
			gtk_tree_selection_select_path(selection, path);
			gtk_tree_view_scroll_to_cell(treeview, path, NULL, TRUE, 0.5, 0.0);
			gtk_tree_path_free(path);
			gui.selected_star = get_index_of_selected_star(xpos, ypos);
			gtk_window_present(GTK_WINDOW(lookup_widget("stars_list_window")));
			display_status();
			redraw(REDRAW_OVERLAY);
			return;
		}
		valid = gtk_tree_model_iter_next(model, &iter);
	}
	siril_debug_print("Point clicked does not correspond to a known star\n");
	return;
}

static int compare(void const *a, void const *b) {
	guint const *pa = a;
	guint const *pb = b;

	return *pa - *pb;
}

static void update_column_index(GtkTreeModel *treeModel, guint *sel, guint size) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_first(treeModel, &iter);

	while (valid) {
		gint idx;
		gtk_tree_model_get(treeModel, &iter, COLUMN_INDEX, &idx, -1);

		/* find the index we need to update */
		int i;
		for (i = size - 1; i >= 0 && idx < sel[i]; i--) ;

		gtk_list_store_set(liststore_stars, &iter, COLUMN_INDEX, idx - i - 1, -1);
		valid = gtk_tree_model_iter_next (treeModel, &iter);
	}
}

static void remove_selected_star() {
	GtkTreeSelection *selection;
	GList *references, *list;

	GtkTreeView *treeView = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "Stars_stored"));
	GtkTreeModel *treeModel = gtk_tree_view_get_model(treeView);

	selection = gtk_tree_view_get_selection(treeView);
	references = get_row_references_of_selected_rows(selection, treeModel);

	guint size = g_list_length(references);
	guint *sel = calloc(size, sizeof(guint));

	int i = 0;
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(treeModel, &iter, path)) {
				GValue g_idx = G_VALUE_INIT;
				gtk_tree_model_get_value(treeModel, &iter, COLUMN_INDEX, &g_idx);
				int idx = g_value_get_int(&g_idx);
				sel[i++] = idx;

				gtk_list_store_remove(liststore_stars, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	qsort (sel, size, sizeof *sel, compare);

	for(int i = size - 1; i >= 0; i--) {
		remove_star(sel[i] - 1);
	}

	update_column_index(treeModel, sel, size);

	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
	free(sel);
	gui.selected_star = -1;
	display_status();
}

static void remove_all_stars(){
	clear_stars_list(TRUE);
	gui.selected_star = -1;
	display_status();
	redraw(REDRAW_OVERLAY);
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("Star list file (*.lst)"));
	gtk_file_filter_add_pattern(f, "*.lst");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

static void save_stars_dialog() {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	GtkWindow *parent = GTK_WINDOW(lookup_widget("stars_list_window"));
	gint res;

	widgetdialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, "stars.lst");
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = gtk_file_chooser_get_filename(dialog);
		save_list(file, MAX_STARS, com.stars, 0, &com.pref.starfinder_conf, -1, TRUE); // passing layer as -1 as we are not sure all stars have been detected on same layer

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

static void add_star_to_list(psf_star *star, int i) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;

	get_stars_list_store();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "treeview-selection"));
	if (star == NULL) {
		gtk_list_store_clear(liststore_stars);
		return;		// just clear the list
	}

	double fwhmx = star->fwhmx_arcsec < 0 ? star->fwhmx : star->fwhmx_arcsec;
	double fwhmy = star->fwhmy_arcsec < 0 ? star->fwhmy : star->fwhmy_arcsec;

	gtk_list_store_append (liststore_stars, &iter);
	gtk_list_store_set (liststore_stars, &iter,
			COLUMN_CHANNEL, star->layer,
			COLUMN_B, star->B,
			COLUMN_A, star->A,
			COLUMN_X0, star->xpos,
			COLUMN_Y0, star->ypos,
			COLUMN_FWHMX, fwhmx,
			COLUMN_FWHMY, fwhmy,
			COLUMN_MAG, star->mag + com.magOffset,
			COLUMN_BETA, star->beta,
			COLUMN_ROUNDNESS, fwhmy / fwhmx,
			COLUMN_ANGLE, -star->angle, // the image had been flipped before the findstar
			COLUMN_RMSE, star->rmse,
			COLUMN_INDEX, i + 1,
			-1);

	units = star->units;
}

static void fill_stars_list(fits *fit, psf_star **stars) {
	int i = 0;
	if (stars == NULL)
		return;
	add_star_to_list(NULL, 0);	// clear

	while (stars[i]) {
		/* update units if needed */
		fwhm_to_arcsec_if_needed(fit, stars[i]);
		add_star_to_list(stars[i], i);
		i++;
	}
	gui.selected_star = -1;
	display_status();
}

/********************* public ***********************/

void refresh_star_list(psf_star **star){
	get_stars_list_store();
	gtk_list_store_clear(liststore_stars);
	fill_stars_list(&gfit, com.stars);
	redraw(REDRAW_OVERLAY);
}

void clear_stars_list(gboolean refresh_GUI) {
	if (com.stars) {
		if (refresh_GUI && !com.headless) {
			get_stars_list_store();
			gtk_list_store_clear(liststore_stars);
		}
		g_mutex_lock(&com.mutex);

		if (com.stars[0]) {
			/* freeing found stars. It must not be done when the only star in
			 * com.stars is the same as com.seq.imgparam[xxx].fwhm, as set in
			 * set_fwhm_star_as_star_list(), because it will be reused */
			if (com.stars[1] || !com.star_is_seqdata) {
				int i = 0;
				while (i < MAX_STARS && com.stars[i])
					free_psf(com.stars[i++]);
			}
		}
		free(com.stars);
		com.stars = NULL;
		g_mutex_unlock(&com.mutex);
	}
	com.star_is_seqdata = FALSE;
	if (refresh_GUI && !com.headless)
		display_status();
}

static int get_comstar_count() {
	int i = 0;
	while (com.stars[i])
		i++;
	return i;
}

void pick_a_star() {
	int layer = match_drawing_area_widget(gui.view[select_vport(gui.cvport)].drawarea, FALSE);
	if (layer != -1) {
		if (!(com.selection.h && com.selection.w))
			return;
		if (com.selection.w > 300 || com.selection.h > 300) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
					_("To determine the PSF, please make a selection around a star."));
			return;
		}
		int new_index;
		psf_star *new_star = add_star(&gfit, layer, &new_index);
		if (new_star) {
			add_star_to_list(new_star, get_comstar_count() - 1);
			display_status();
			siril_open_dialog("stars_list_window");
		} else
			return;
	}
	redraw(REDRAW_OVERLAY);
}

static gchar *build_wcs_url(gchar *ra, gchar *dec) {
	if (!has_wcs(&gfit)) return NULL;

	double resolution = get_wcs_image_resolution(&gfit);

	gchar *tol = g_strdup_printf("%lf", resolution * 3600 * 15);

	GString *url = g_string_new("https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=");
	url = g_string_append(url, ra);
	url = g_string_append(url, dec);
	url = g_string_append(url, "&Radius=");
	url = g_string_append(url, tol);
	url = g_string_append(url, "&Radius.unit=arcsec");
	url = g_string_append(url, "#lab_basic");

	gchar *simbad_url = g_string_free(url, FALSE);
	gchar *cleaned_url = url_cleanup(simbad_url);

	g_free(tol);
	g_free(simbad_url);

	return cleaned_url;
}

static const char *SNR_quality(double SNR) {
	if (SNR > 40.0) return _("Excellent");
	if (SNR > 25.0) return _("Good");
	if (SNR > 15.0) return _("Fair");
	if (SNR > 10.0) return _("Poor");
	if (SNR > 0.0) return _("Bad");
	else return _("N/A");
}

void popup_psf_result(psf_star *result, rectangle *area, fits *fit) {
	gchar *msg, *coordinates, *url = NULL;
	char buffer2[50];
	const char *str;
	if (com.magOffset > 0.0)
		str = _("true reduced");
	else
		str = _("relative");

	double x = result->x0 + area->x;
	double y = area->y + area->h - result->y0;
	if (has_wcs(&gfit)) {
		double world_x, world_y;
		SirilWorldCS *world_cs;

		pix2wcs(&gfit, x, (double) gfit.ry - y, &world_x, &world_y);
		world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
		if (world_cs) {
			gchar *ra = siril_world_cs_alpha_format(world_cs, "%02d %02d %.3lf");
			gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d %02d %.3lf");

			url = build_wcs_url(ra, dec);
			// TODO: change with vizier
			// TODO: use box size as radius

			g_free(ra);
			g_free(dec);

			ra = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
			dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");

			coordinates = g_strdup_printf("x0=%.2fpx\t%s J2000\n\t\ty0=%.2fpx\t%s J2000", x, ra, y, dec);

			g_free(ra);
			g_free(dec);
			siril_world_cs_unref(world_cs);
		} else {
			coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
		}
	} else {
		coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
	}

	double fwhmx, fwhmy;
	char *unts;
	get_fwhm_as_arcsec_if_possible(result, &fwhmx, &fwhmy, &unts);
	const gchar *chan = isrgb(fit) ? channel_number_to_name(result->layer) : _("monochrome");
	if (result->beta > 0.0) {
		g_snprintf(buffer2, 50, ", beta=%0.1f, %s channel", result->beta, chan);
	}
	else {
		g_snprintf(buffer2, 50, "%s, %s channel", "", chan);
	}
	msg = g_strdup_printf(_("PSF fit Result (%s%s):\n\n"
				"Centroid Coordinates:\n\t\t%s\n\n"
				"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\t\tr=%.2f\n"
				"Angle:\n\t\t%0.2fdeg\n\n"
				"Background Value:\n\t\tB=%.6f\n\n"
				"Maximal Intensity:\n\t\tA=%.6f\n\n"
				"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
				"Signal-to-noise ratio:\n\t\tSNR=%.1fdB (%s)\n\n"
				"RMSE:\n\t\tRMSE=%.3e"),
			(result->profile == PSF_GAUSSIAN) ? "Gaussian" : "Moffat", buffer2,
			coordinates, fwhmx, unts, fwhmy, unts, fwhmy / fwhmx,
			result->angle, result->B, result->A, str,
			result->mag + com.magOffset, result->s_mag, result->SNR,
			SNR_quality(result->SNR), result->rmse);
	g_free(coordinates);
	show_data_dialog(msg, _("PSF and quick photometry results"), NULL, url);
	g_free(msg);
	g_free(url);
}

/***************** callbacks ****************/

void on_treeview_selection_changed(GtkTreeSelection *selection, gpointer user_data) {
	GtkTreeIter iter;
	GValue value_x = G_VALUE_INIT;
	GValue value_y = G_VALUE_INIT;
	const gchar *area[] = {"drawingarear", "drawingareag", "drawingareab", "drawingareargb" };
	GtkWidget *widget = lookup_widget(area[gui.cvport]);

	GtkTreeView *treeView = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "Stars_stored"));
	GtkTreeModel *treeModel = gtk_tree_view_get_model(treeView);

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE)
		return;	//The tree is empty

	GList *list = gtk_tree_selection_get_selected_rows(selection, &treeModel);
	if (g_list_length(list) == 1) {
		gtk_tree_model_get_iter(treeModel, &iter, (GtkTreePath*) list->data);
		gdouble x0, y0;

		gtk_tree_model_get_value(treeModel, &iter, COLUMN_X0, &value_x);
		x0 = g_value_get_double(&value_x);
		gtk_tree_model_get_value(treeModel, &iter, COLUMN_Y0, &value_y);
		y0 = g_value_get_double(&value_y);

		g_value_unset(&value_x);
		g_value_unset(&value_y);

		// Set this to draw blue crosshairs
		gui.selected_star = get_index_of_selected_star(x0, y0);
		// Centre selected star
		GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("toggle_star_centered"));
		if (gtk_toggle_button_get_active(toggle)) {
			double z = get_zoom_val();
			gui.display_offset.x = (gtk_widget_get_allocated_width(widget) / 2 - x0 * z);
			gui.display_offset.y = (gtk_widget_get_allocated_height(widget) / 2 - y0 * z);
			adjust_vport_size_to_image();
		}

		display_status();
		redraw(REDRAW_OVERLAY);
	}
}

void on_Stars_stored_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {

		remove_selected_star();
	}
}

void on_stars_list_window_hide(GtkWidget *object, gpointer user_data) {
	gui.selected_star = -1;
}

void on_sum_button_clicked(GtkButton *button, gpointer user_data) {
	display_PSF(com.stars);
}

void on_add_button_clicked(GtkButton *button, gpointer user_data) {
	int layer = match_drawing_area_widget(gui.view[gui.cvport].drawarea, FALSE);
	int index;
	add_star(&gfit, layer, &index);
	if (index > -1)
		add_star_to_list(com.stars[index], index);
	display_status();
	refresh_star_list(com.stars);
}

void on_remove_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_star();
}

void on_remove_all_button_clicked(GtkButton *button, gpointer user_data) {
	remove_all_stars();
	clear_sensor_tilt();
}

void on_export_button_clicked(GtkButton *button, gpointer user_data) {
	if (com.stars) {
		save_stars_dialog();
	} else {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Nothing to export"),
				_("There are no stars in the list."));
	}
}

void on_stars_list_window_show(GtkWidget *widget, gpointer user_data) {
	update_peaker_GUI();
}

void on_button_stars_list_ok_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("stars_list_window");
}

