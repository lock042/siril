/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/dialog_preview.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/open_dialog.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "tinyexpr.h"
#include "pixel_math_runner.h"

static const gchar *variables[] = {"I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10"};
#define MAX_IMAGES G_N_ELEMENTS(variables)

enum {
	COLUMN_IMAGE_NUM,		// string
	COLUMN_IMAGE_PATH,		// string
	N_COLUMNS
};

static fits var_fit[MAX_IMAGES] = { 0 };

static GtkListStore *pixel_math_list_store = NULL;
static GtkTreeView *pixel_math_tree_view = NULL;
static GtkTreeModel *pixel_math_tree_model = NULL;
static GtkLabel *pixel_math_status_bar = NULL;
static GtkTextView *pixel_math_text_view = NULL;

static void init_widgets() {
	if (!pixel_math_tree_view) {
		pixel_math_tree_view = GTK_TREE_VIEW(gtk_builder_get_object(builder, "pixel_math_treeview"));
		pixel_math_tree_model = gtk_tree_view_get_model(pixel_math_tree_view);
		pixel_math_list_store = GTK_LIST_STORE(gtk_builder_get_object(builder, "pixel_math_liststore"));
		pixel_math_status_bar = GTK_LABEL(lookup_widget("pixel_math_status"));
		pixel_math_text_view = GTK_TEXT_VIEW(lookup_widget("pixel_math_textview"));

	}
	g_assert(pixel_math_tree_view);
	g_assert(pixel_math_tree_model);
	g_assert(pixel_math_list_store);
	g_assert(pixel_math_status_bar);
	g_assert(pixel_math_text_view);
}

void remove_spaces_from_str(gchar *s) {
	gchar *d = s;
	do {
		while (g_ascii_isspace(*d)) {
			++d;
		}
	} while((*s++ = *d++));
}

static gchar* get_pixel_math_expression() {
	GtkTextIter start, end;
	init_widgets();

	GtkTextBuffer *buf = gtk_text_view_get_buffer(pixel_math_text_view);

	gtk_text_buffer_get_bounds(buf, &start, &end);
	return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}

static void output_status_bar(int status) {
	switch (status) {
	case 0:
		gtk_label_set_text(pixel_math_status_bar, "");
		break;
	default:
		gtk_label_set_text(pixel_math_status_bar, _("Syntax error"));
	}
}

static gboolean end_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;
	stop_processing_thread();// can it be done here in case there is no thread?

	if (!args->ret) {
		//	/* Create new image */
		com.seq.current = UNRELATED_IMAGE;
		com.uniq = calloc(1, sizeof(single));
		com.uniq->filename = strdup(_("new empty image"));
		com.uniq->fileexist = FALSE;
		com.uniq->nb_layers = args->fit->naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = args->fit;
		clearfits(&gfit);
		copyfits(args->fit, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		open_single_image_from_gfit();
	}
	set_cursor_waiting(FALSE);
	g_free(args->expression);
	output_status_bar(args->ret);

	free(args);
	return FALSE;
}

static gpointer apply_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;

	fits *fit = args->fit;
	int nb_rows = args->nb_rows;
	const gchar *expression = args->expression;

	int err;
	gboolean failed = FALSE;

	size_t nbdata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		te_variable *vars = malloc(nb_rows * sizeof(te_variable));
		double *x = malloc(nb_rows * sizeof(double));

		for (int i = 0; i < nb_rows; i++) {
			vars[i].name = variables[i];
			vars[i].address = &x[i];
			vars[i].context = NULL;
			vars[i].type = 0;
		}
		te_expr *n = te_compile(expression, vars, nb_rows, &err);
		if (!n) {
			failed = TRUE;
		} else {
			size_t px;
#ifdef _OPENMP
#pragma omp for private(px) schedule(static)
#endif
			for (px = 0; px < nbdata; px++) {
				/* The variables can be changed here, and eval can be called as many
				 * times as you like. This is fairly efficient because the parsing has
				 * already been done. */
				for (int i = 0; i < nb_rows; i++) {
					x[i] = var_fit[i].fdata[px];
				}

				fit->fdata[px] = (float) te_eval(n);
			}

			te_free(n);
			args->ret = 0;
		}
		free(vars);
		free(x);
	}
	if (failed) {
		/* Show the user where the error is at. */
		printf("%s\n", expression);
		printf("\t%*s^\nError near here\n", err - 1, "");
		args->ret = err;
	}

	siril_add_idle(end_pixel_math_operation, args);
	return GINT_TO_POINTER(args->ret);
}

static const gchar *get_pixel_math_var_paths(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;

	init_widgets();

	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path);
	gtk_tree_model_get_value(pixel_math_tree_model, &iter, COLUMN_IMAGE_PATH, &value);

	return g_value_get_string(&value);
}

static int get_pixel_math_number_of_rows(){
	if (GTK_IS_TREE_MODEL(pixel_math_list_store))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store), NULL);
	else return 0;
}

static int pixel_math_evaluate(gchar *expression) {
	int nb_rows = 0;

	int width = -1;
	int height = -1;
	int channel = -1;

	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < MAX_IMAGES) {
		const gchar *path = get_pixel_math_var_paths(nb_rows);
		if (readfits(path, &var_fit[nb_rows], NULL, TRUE)) return -1;
		if (width == - 1) {
			width = var_fit[nb_rows].rx;
			height = var_fit[nb_rows].ry;
			channel = var_fit[nb_rows].naxes[2];
		} else {
			int w = var_fit[nb_rows].rx;
			int h = var_fit[nb_rows].ry;
			int c = var_fit[nb_rows].naxes[2];

			if ((width != w) || (height != h) || (channel != c)) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Image size must be identical"),
						_("The images used in the Pixel Math tool must have the same size."));
				return 1;
			}
		}
		nb_rows++;
	}

	if (!nb_rows) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No images loaded"), _("You must load images first."));
		return 1;
	}

	fits *fit = NULL;
	if (new_fit_image(&fit, width, height, channel, DATA_FLOAT))
		return 1;

	struct pixel_math_data *args = malloc(sizeof(struct pixel_math_data));

	args->expression = expression;
	args->fit = fit;
	args->nb_rows = nb_rows;
	args->ret = 0;

	start_in_new_thread(apply_pixel_math_operation, args);

	return 0;
}

static const gchar *get_pixel_math_var_name(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;

	init_widgets();

	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path);
	gtk_tree_model_get_value(pixel_math_tree_model, &iter, COLUMN_IMAGE_NUM, &value);

	return g_value_get_string(&value);
}

/* Add an image to the list. */
static void add_image_to_variable_list(const gchar *path, int id) {
	GtkTreeIter iter;

	gtk_list_store_append(pixel_math_list_store, &iter);
	gtk_list_store_set(pixel_math_list_store, &iter,
			COLUMN_IMAGE_NUM, variables[id],
			COLUMN_IMAGE_PATH, path,
			-1);

}

static void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern) {
	gchar **patterns;
	gint i;

	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	patterns = g_strsplit(pattern, ";", -1);
	for (i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	gtk_file_chooser_set_filter(file_chooser, f);
}

static void select_image(int id) {
	GtkWidget *dialog;
	fileChooserPreview *preview = NULL;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
	gint res;

	dialog = gtk_file_chooser_dialog_new("Open File", GTK_WINDOW(lookup_widget("pixel_math_dialog")), action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,	NULL);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), com.wd);
	gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(dialog), FALSE);
	gtk_filter_add(GTK_FILE_CHOOSER(dialog), _("FITS Files (*.fit, *.fits, *.fts, *.fits.fz)"), FITS_EXTENSIONS);
	siril_file_chooser_add_preview(GTK_FILE_CHOOSER(dialog), preview);

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		char *filename;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		filename = gtk_file_chooser_get_filename(chooser);

		add_image_to_variable_list(filename, id);
		g_free(filename);
	}

	gtk_widget_destroy(dialog);
}

void on_pixel_math_add_var_button_clicked(GtkButton *button, gpointer user_data) {
	init_widgets();

	int id = get_pixel_math_number_of_rows();
	if (id == MAX_IMAGES) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Cannot load new image"),
				_("You've reached the maximum of loaded image file."));
	} else {
		select_image(id);
	}
}

void on_pixel_math_remove_var_button_clicked(GtkButton *button, gpointer user_data) {
	GtkTreeSelection *selection;
	GList *references, *list;

	init_widgets();
	selection = gtk_tree_view_get_selection(pixel_math_tree_view);
	references = get_row_references_of_selected_rows(selection, pixel_math_tree_model);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path)) {
				gtk_list_store_remove(pixel_math_list_store, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
}

void on_apply_pixel_math_clicked(GtkButton *button, gpointer user_data) {
	gchar *expression = get_pixel_math_expression();
	remove_spaces_from_str(expression);

	pixel_math_evaluate(expression);
}

void on_pixel_math_treeview_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	GtkTextIter iter;

	init_widgets();

	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(pixel_math_text_view);
	gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_pixel_math_var_name(i[0]);

	gtk_text_buffer_get_iter_at_mark(tbuf, &iter, gtk_text_buffer_get_insert(tbuf));

	gtk_text_buffer_insert(tbuf, &iter, str, strlen(str));
	gtk_widget_grab_focus(GTK_WIDGET(pixel_math_text_view));
}
