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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "core/processing.h"
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

static GtkListStore *liststore_stars = NULL;

enum {
	COLUMN_CHANNEL,		// int
	COLUMN_B,			// gdouble
	COLUMN_A,			// gdouble
	COLUMN_X0,			// gdouble
	COLUMN_Y0,			// gdouble
	COLUMN_RA,          // gdouble
	COLUMN_DEC,			// gdouble
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

static void gdouble_ra_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_RA, &var, -1);
	if (var > 9E9) {
		buf = g_strdup("N/A");
	} else {
		SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(var, 0);

		if (com.pref.gui.show_deciasec) {
			buf = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%04.1lfs");
		} else {
			buf = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
		}
		siril_world_cs_unref(world_cs);

	}
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_dec_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_DEC, &var, -1);
	if (var > 9E9) {
		buf = g_strdup("N/A");
	} else {
		SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(0, var);

		if (com.pref.gui.show_deciasec) {
			buf = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%04.1lf\"");
		} else {
			buf = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
		}
		siril_world_cs_unref(world_cs);
	}
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

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn30"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_ra"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_ra_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn31"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cell_dec"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_dec_cell_data_function, NULL, NULL);

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
	double bin_X = com.pref.binning_update ? (double) gfit.keywords.binning_x : 1.0;
	if (com.stars && com.stars[0]) {// If the first star has units of arcsec, all should have
		is_as = (strcmp(com.stars[0]->units, "px"));
	} else {
		return; // If com.stars is empty there is no point carrying on
	}
	if (is_as) {
		invpixscalex = 1.0 / (radian_conversion * (double) gfit.keywords.pixel_size_x / gfit.keywords.focal_length) * bin_X;
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
			gtk_tree_selection_unselect_all(selection);
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

static void update_column_index(GtkTreeModel *treeModel, const guint *sel, guint size) {
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

	for (int j = size - 1; j >= 0; j--) {
		remove_star(sel[j] - 1);
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

static void export_to_csv(GtkTreeView *treeview, const char *filename) {
	GError *error = NULL;
	gboolean ret = TRUE;
	GFile *file = g_file_new_for_path(filename);
	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "Cannot export CSV\n");
		}
		g_object_unref(file);
		return;
	}

	GDataOutputStream *data_stream = g_data_output_stream_new(
			G_OUTPUT_STREAM(output_stream));

	// Write header
	gint n_columns = gtk_tree_view_get_n_columns(treeview);
	for (gint i = 0; i < n_columns; i++) {
		if (gtk_tree_view_column_get_visible(gtk_tree_view_get_column(treeview, i))) {
			const gchar *column_name = gtk_tree_view_column_get_title(gtk_tree_view_get_column(treeview, i));
			ret &= g_data_output_stream_put_string(data_stream, column_name, NULL, NULL);
			if (i < n_columns - 1) {
				gint next_column = i + 1;
				while (!gtk_tree_view_column_get_visible(gtk_tree_view_get_column(treeview, next_column))) {
					next_column++;
					if (next_column >= n_columns) {
						break;
					}
				}
				if (next_column < n_columns) {
					ret &= g_data_output_stream_put_string(data_stream, ",", NULL, NULL);
				}
			}
		}
	}
	ret &= g_data_output_stream_put_string(data_stream, "\n", NULL, NULL);

	// Write row data
	GtkTreeModel *model = gtk_tree_view_get_model(treeview);
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);
	while (valid) {
		for (gint i = 0; i < n_columns; i++) {
			if (gtk_tree_view_column_get_visible(
					gtk_tree_view_get_column(treeview, i))) {
				GValue value = G_VALUE_INIT;
				gtk_tree_model_get_value(model, &iter, i, &value);
				if (G_VALUE_TYPE(&value) == G_TYPE_DOUBLE) {
					gdouble dbl_value = g_value_get_double(&value);
					gchar *str_value = g_strdup_printf("%g", dbl_value);
					ret &= g_data_output_stream_put_string(data_stream, str_value, NULL, NULL);
					g_free(str_value);
				} else if (G_VALUE_TYPE(&value) == G_TYPE_INT) {
					gint int_value = g_value_get_int(&value);
					gchar *str_value = g_strdup_printf("%d", int_value);
					ret &= g_data_output_stream_put_string(data_stream, str_value, NULL, NULL);
					g_free(str_value);
				}
				if (i < n_columns - 1) {
					gint next_column = i + 1;
					while (!gtk_tree_view_column_get_visible(gtk_tree_view_get_column(treeview, next_column))) {
						next_column++;
						if (next_column >= n_columns) {
							break;
						}
					}
					if (next_column < n_columns && gtk_tree_view_column_get_visible(gtk_tree_view_get_column(treeview, next_column))) {
						ret &= g_data_output_stream_put_string(data_stream, ",", NULL, NULL);
					}
				}
				g_value_unset(&value);
			}
		}
		ret &= g_data_output_stream_put_string(data_stream, "\n", NULL, NULL);

		valid = gtk_tree_model_iter_next(model, &iter);
	}
	if (!ret)
		siril_log_color_message(_("Error: error writing the CSV.\n"), "red");
	g_object_unref(data_stream);
	g_object_unref(output_stream);
	g_object_unref(file);
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
	gtk_file_chooser_set_local_only(dialog, FALSE);
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_file_chooser_get_filename(dialog);
		export_to_csv(GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "Stars_stored")), file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

int get_ra_and_dec_from_star_pos(psf_star *star, gdouble *alpha, gdouble *delta) {
	int ret = 1;
	if (has_wcs(&gfit)) {
		// coordinates of the star in FITS/WCS coordinates
		double fx, fy;
		display_to_siril(star->xpos, star->ypos, &fx, &fy, gfit.ry);

		double ra, dec;
		pix2wcs(&gfit, fx, fy, &ra, &dec);
		// *alpha = ra would work too instead of all this?
		SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(ra, dec);
		if (world_cs) {
			double a = siril_world_cs_get_alpha(world_cs);
			double d = siril_world_cs_get_delta(world_cs);

			siril_world_cs_unref(world_cs);

			*alpha = a;
			*delta = d;
			ret = 0;
		}
	}
	return ret;
}

static void add_star_to_list(psf_star *star, int i) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;
	double ra, dec;

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
			COLUMN_RA, star->ra,
			COLUMN_DEC, star->dec,
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
	if (!stars) return;
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

// consider using update_star_list instead
void refresh_star_list(){
	get_stars_list_store();
	gtk_list_store_clear(liststore_stars);
	fill_stars_list(&gfit, com.stars);
}

/* this can be called from any thread as long as refresh_GUI is false, it's
 * synchronized with the main thread with a mutex */
void clear_stars_list(gboolean refresh_GUI) {
	g_mutex_lock(&com.mutex); // Lock at the beginning to protect the check and modification
	psf_star **stars = com.stars;

	if (stars) {
		com.stars = NULL; // Set com.stars to NULL while holding the lock
		g_mutex_unlock(&com.mutex);

		if (refresh_GUI && !com.headless) {
			get_stars_list_store();
			gtk_list_store_clear(liststore_stars);
		}

		if (stars[0]) {
			/* freeing found stars. It must not be done when the only star in
			* com.stars is the same as com.seq.imgparam[xxx].fwhm, as set in
			* set_fwhm_star_as_star_list(), because it will be reused */
			if (stars[1] || !com.star_is_seqdata) {
				int i = 0;
				while (i < MAX_STARS && stars[i])
					free_psf(stars[i++]);
			}
			free(stars);
		}
	} else {
		g_mutex_unlock(&com.mutex); // Unlock if com.stars is NULL
	}

	com.star_is_seqdata = FALSE;
	gui.selected_star = -1;
	if (refresh_GUI && !com.headless)
		display_status();
}

struct star_update_s {
	psf_star **stars;
	gboolean update_GUI;
};

static gboolean update_stars_idle(gpointer p) {
	struct star_update_s *args = (struct star_update_s *)p;
	clear_stars_list(TRUE);
	com.stars = args->stars;
	if (args->update_GUI && !com.headless)
		fill_stars_list(&gfit, com.stars);
	redraw(REDRAW_OVERLAY);
	free(args);
	return FALSE;
}

/* Update the list of stars (com.stars) in a safe way.
 * assuming stars to not be shared with a sequence (star_is_seqdata false)
 * the PSF window's list will be cleared, only refilled if update_PSF_list is true
 */
void update_star_list(psf_star **new_stars, gboolean update_PSF_list, gboolean wait_for_update) {
	struct star_update_s *args = calloc(1, sizeof(struct star_update_s));
	args->stars = new_stars;
	args->update_GUI = update_PSF_list;
	if (!com.headless) {
		if (wait_for_update)
			execute_idle_and_wait_for_it(update_stars_idle, args);
		else
			siril_add_idle(update_stars_idle, args);
	}
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

void popup_psf_result(psf_star *result, rectangle *area, fits *fit) {
	gchar *url = NULL;
	gchar *msg = format_psf_result(result, area, fit, &url);
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
	if (layer == -1)
		layer = 1;
	int index;
	add_star(&gfit, layer, &index);
	if (index > -1)
		add_star_to_list(com.stars[index], index);
	display_status();
	refresh_star_list();
}

void on_remove_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_star();
}

void on_remove_all_button_clicked(GtkButton *button, gpointer user_data) {
	remove_all_stars();
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

