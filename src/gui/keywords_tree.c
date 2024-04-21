/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/dialogs.h"
#include "gui/utils.h"

#include "keywords_tree.h"

static GtkListStore *key_liststore = NULL;

static int ffs2c(const char *instr, /* I - null terminated input string  */
char *outstr, /* O - null terminated quoted output string */
const int *status) /* IO - error status */
/*
 convert an input string to a quoted string. Leading spaces
 are significant.  FITS string keyword values must be at least
 8 chars long so pad out string with spaces if necessary.
 Example:   km/s ==> 'km/s    '
 Single quote characters in the input string will be replace by
 two single quote characters. e.g., o'brian ==> 'o''brian'
 */
{
	size_t len, ii, jj;

	if (*status > 0) /* inherit input status value if > 0 */
		return (*status);

	if (!instr) /* a null input pointer?? */
	{
		strcpy(outstr, "''"); /* a null FITS string */
		return (*status);
	}

	outstr[0] = '\''; /* start output string with a quote */

	len = strlen(instr);
	if (len > 68)
		len = 68; /* limit input string to 68 chars */

	for (ii = 0, jj = 1; ii < len && jj < 69; ii++, jj++) {
		outstr[jj] = instr[ii]; /* copy each char from input to output */
		if (instr[ii] == '\'') {
			jj++;
			outstr[jj] = '\''; /* duplicate any apostrophies in the input */
		}
	}

	for (; jj < 9; jj++) /* pad string so it is at least 8 chars long */
		outstr[jj] = ' ';

	if (jj == 70) /* only occurs if the last char of string was a quote */
		outstr[69] = '\0';
	else {
		outstr[jj] = '\''; /* append closing quote character */
		outstr[jj + 1] = '\0'; /* terminate the string */
	}

	return (*status);
}


enum {
	COLUMN_KEY,		// string
	COLUMN_VALUE,		// string
	COLUMN_COMMENT,		// string
	COLUMN_DTYPE,       // char
	COLUMN_PROTECTED, // gboolean
	COLUMN_COLOR, // string
	COLUMN_EDITABLE, //gboolean
	N_COLUMNS
};

static void get_keylist_store() {
	if (key_liststore == NULL)
		key_liststore = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "key_liststore"));
}


static void add_key_to_tree(const gchar *key, const gchar *value, const gchar *comment, const char dtype, gboolean protected, gboolean editable) {
	GtkTreeIter iter;

	gtk_list_store_append(key_liststore, &iter);
	gtk_list_store_set(key_liststore, &iter,
			COLUMN_KEY, key,
			COLUMN_VALUE, value,
			COLUMN_COMMENT, comment,
			COLUMN_DTYPE, dtype,
			COLUMN_PROTECTED, protected,
			COLUMN_COLOR, protected ? "salmon" : NULL,
			COLUMN_EDITABLE, !protected && editable,
			-1);
}

static void init_dialog() {
	static GtkTreeSelection *selection = NULL;
	get_keylist_store();
	if (!selection) {
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "key_selection"));
	}
	gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);
	gtk_list_store_clear(key_liststore);
}

static int listFITSKeywords(fits *fit, gboolean editable) {
	char card[FLEN_CARD];
	void *memptr;
	size_t memsize = IOBUFLEN;
	int status = 0;
    int nkeys, ii;
	fitsfile *fptr = NULL;
	memptr = malloc(memsize);
	if (!memptr) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fits_create_memfile(&fptr, &memptr, &memsize, IOBUFLEN, realloc, &status);
	if (status) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return 1;
	}

	if (fits_create_img(fptr, fit->bitpix, fit->naxis, fit->naxes, &status)) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return 1;
	}

	fits tmpfit = { 0 };
	copy_fits_metadata(fit, &tmpfit);
	tmpfit.fptr = fptr;
	save_fits_header(&tmpfit);

	fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */

	for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */
		char keyname[FLEN_KEYWORD];
		char value[FLEN_VALUE];
		char comment[FLEN_COMMENT];
		char dtype;
		int length = 0;

		if (fits_read_record(fptr, ii, card, &status))
			break;
		fits_get_keyname(card, keyname, &length, &status);
		fits_parse_value(card, value, comment, &status);
		fits_get_keytype(value, &dtype, &status);
		status = 0;
		add_key_to_tree(keyname, value, comment, dtype, fits_get_keyclass(card) == TYP_STRUC_KEY, editable);
	}


	if (status == END_OF_FILE)
		status = 0; /* Reset after normal error */

	fits_movabs_hdu(tmpfit.fptr, 1, NULL, &status);

	if (status)
		fits_report_error(stderr, status); /* print any error message */

	fits_movabs_hdu(fptr, 0, NULL, &status);
	fits_close_file(fptr, &status);
	clearfits(&tmpfit);

	free(memptr);

	return (status);
}

void on_keywords_dialog_show(GtkWidget *dialog, gpointer user_data) {
	refresh_keywords_dialog();
}


void on_val_edited(GtkCellRendererText *renderer, char *path, char *new_val, gpointer user_data) {
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(lookup_widget("key_treeview")));
	GtkTreeIter iter;
	gchar *FITS_key;
	gboolean protected;
	char dtype;

	gtk_tree_model_get_iter_from_string(model, &iter, path);
	gtk_tree_model_get(model, &iter, COLUMN_KEY, &FITS_key, COLUMN_DTYPE, &dtype, COLUMN_PROTECTED, &protected, -1);
	if (!protected) {
		char valstring[FLEN_VALUE];
		int status = 0;
		/* update FITS */
		if (dtype == 'C' && (new_val[0] != '\'' || new_val[strlen(new_val) - 1] != '\'')) {
			ffs2c(new_val, valstring, &status);
		} else {
			strcpy(valstring, new_val);
		}
		if (!updateFITSKeyword(&gfit, FITS_key, valstring)) {
			gtk_list_store_set(key_liststore, &iter, COLUMN_VALUE, valstring, -1);
		}
	}
}

void on_key_close_btn_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("keywords_dialog");
}

static gchar *list_all_keywords() {
	GtkTreeSelection *selection;
	GList *references, *list;
	GtkTreeView *tree_view = GTK_TREE_VIEW(lookup_widget("key_treeview"));
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GString *string = g_string_new("");

	selection = gtk_tree_view_get_selection(tree_view);
	references = get_row_references_of_selected_rows(selection, model);

	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*) list->data);
		if (path) {
			if (gtk_tree_model_get_iter(model, &iter, path)) {
				gchar *key, *value, *comment;
				int status = 0;
				char card[FLEN_CARD];
				gtk_tree_model_get(model, &iter, COLUMN_KEY, &key, COLUMN_VALUE, &value, COLUMN_COMMENT, &comment, -1);

				fits_make_key(key, value, comment, card, &status);

				g_string_prepend(string, "\n");
				g_string_prepend(string, card);
				g_free(key);
				g_free(value);
				g_free(comment);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);

	return g_string_free(string, FALSE);
}

void on_key_selection_changed(GtkTreeSelection *selection, gpointer user_data) {
	GtkTreeIter iter;
	GtkWidget *widget;
	GList *list;

	GtkTreeView *tree_view = GTK_TREE_VIEW(lookup_widget("key_treeview"));
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	gboolean is_empty = gtk_tree_model_get_iter_first(model, &iter) == FALSE;
	gboolean are_selected = FALSE;

	if (!is_empty) {
		list = gtk_tree_selection_get_selected_rows(selection, &model);
		are_selected = g_list_length(list)  > 0;
	}

	widget = lookup_widget("export_keywords_button");
	gtk_widget_set_sensitive(widget, !is_empty && are_selected);
}



static void save_key_to_clipboard() {
	/* Get the clipboard object */
	GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);

	gchar *list = list_all_keywords();
	gtk_clipboard_set_text(clipboard, list, -1);
	g_free(list);
}

void on_export_keywords_button_clicked(GtkButton *button, gpointer user_data) {
	save_key_to_clipboard();
}

void refresh_keywords_dialog() {
	init_dialog();
	gboolean is_a_single_image_loaded = single_image_is_loaded() &&
			(!sequence_is_loaded() || (sequence_is_loaded() &&
			(com.seq.current == RESULT_IMAGE || com.seq.current == SCALED_IMAGE)));
	listFITSKeywords(&gfit, is_a_single_image_loaded);
}
