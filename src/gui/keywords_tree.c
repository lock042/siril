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
#include "io/fits_keywords.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/utils.h"

#include "keywords_tree.h"

static GtkListStore *key_liststore = NULL;
static GtkTreeView *key_treeview = NULL;
static GtkTextView *key_textview = NULL;
static GtkNotebook *key_notebook = NULL;
static GtkTreeSelection *key_selection = NULL;
static GtkWidget *key_export_button = NULL;


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

static void on_mark_set(GtkTextBuffer *buffer, GtkTextIter *iter, GtkTextMark *mark, gpointer user_data) {
    GtkTextIter start, end;

    GtkWidget *widget = lookup_widget("export_keywords_button");
    gtk_widget_set_sensitive(widget, gtk_text_buffer_get_selection_bounds(buffer, &start, &end));
}

static void init_static_ui() {
	if (key_liststore == NULL) {
		key_liststore = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "key_liststore"));
		key_treeview = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "key_treeview"));
		key_textview = GTK_TEXT_VIEW(lookup_widget("FITS_header_txt"));
		key_notebook = GTK_NOTEBOOK(lookup_widget("notebook-keywords"));
		key_selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "key_selection"));
		key_export_button = lookup_widget("export_keywords_button");

	    g_signal_connect(gtk_text_view_get_buffer(GTK_TEXT_VIEW(key_textview)), "mark-set", G_CALLBACK(on_mark_set), NULL);
	}
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
	init_static_ui();

	gtk_tree_selection_set_mode(key_selection, sequence_is_loaded() ? GTK_SELECTION_SINGLE : GTK_SELECTION_MULTIPLE);
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
	tmpfit.bitpix = fit->bitpix;
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
		add_key_to_tree(keyname, value, comment, dtype, keyword_is_protected(card), editable);
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

static void remove_selected_keys () {
	GtkTreeSelection *selection;
	GList *references;

	GtkTreeModel *treeModel = gtk_tree_view_get_model(key_treeview);

	selection = gtk_tree_view_get_selection(key_treeview);

	if (sequence_is_loaded()) {
		GtkTreeIter iter;
		GValue g_key = G_VALUE_INIT;
		gboolean valid_selection = gtk_tree_selection_get_selected(selection, &treeModel, &iter);

		if (valid_selection) {
			gtk_tree_model_get_value(treeModel, &iter, COLUMN_KEY, &g_key);
			if (G_VALUE_HOLDS_STRING(&g_key)) {
				struct keywords_data *kargs = calloc(1, sizeof(struct keywords_data));

				kargs->FITS_key = g_strdup((gchar *)g_value_get_string(&g_key));
				kargs->value = NULL;
				kargs->comment = NULL;

				if (siril_confirm_dialog(_("Operation on the sequence"),
						_("These keywords will be deleted from each image of "
						"the entire sequence. Are you sure?”"), _("Proceed"))) {
					gtk_list_store_remove(GTK_LIST_STORE(treeModel), &iter);

					start_sequence_keywords(&com.seq, kargs);
				} else {
					free(kargs);
				}
		        g_value_unset(&g_key);
			}
		}
	} else {
		references = get_row_references_of_selected_rows(selection, treeModel);

		if (g_list_length(references) <= 0) return;

		for (GList *list = references; list; list = list->next) {
			GtkTreeIter iter;
			GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
			if (path) {
				if (gtk_tree_model_get_iter(treeModel, &iter, path)) {
					GValue g_key = G_VALUE_INIT;
					gtk_tree_model_get_value(treeModel, &iter, COLUMN_KEY, &g_key);
				    if (G_VALUE_HOLDS_STRING(&g_key)) {
				        gchar *FITS_key = (gchar *)g_value_get_string(&g_key);
						updateFITSKeyword(&gfit, FITS_key, NULL, NULL, NULL, TRUE, FALSE);
						gtk_list_store_remove(GTK_LIST_STORE(treeModel), &iter);

				        g_value_unset(&g_key);
				    }
				}
				gtk_tree_path_free(path);
			}
		}
		gtk_tree_selection_unselect_all(selection);
	    g_list_free_full(references, (GDestroyNotify)gtk_tree_row_reference_free);
	}
}

void on_key_treeview_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {

		remove_selected_keys();
	}
}

void on_cell_editing_started(GtkCellRenderer *renderer, GtkCellEditable *editable, const gchar *path, gpointer user_data) {
	GtkWidget *treeview = GTK_WIDGET(user_data);
	g_signal_handlers_block_by_func(treeview, G_CALLBACK(on_key_treeview_key_release_event), NULL);
}

void on_cell_editing_canceled(GtkCellRenderer *renderer, gpointer user_data) {
	GtkWidget *treeview = GTK_WIDGET(user_data);
	g_signal_handlers_unblock_by_func(treeview, G_CALLBACK(on_key_treeview_key_release_event), NULL);
}

void on_key_edited(GtkCellRendererText *renderer, char *path, char *new_val, gpointer user_data) {
	GtkWidget *treeview = GTK_WIDGET(user_data);
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(treeview));
	GtkTreeIter iter;
	gchar *old_keyname;
	gboolean protected;

	gtk_tree_model_get_iter_from_string(model, &iter, path);
	gtk_tree_model_get(model, &iter, COLUMN_KEY, &old_keyname, COLUMN_PROTECTED, &protected, -1);
	if (!protected) {
		/* update FITS keyname */
		if (g_strcmp0(old_keyname, new_val)) {
			if (strlen(new_val) > 8) {
				siril_log_color_message(_("Keyname can contain a maximum of 8 characters.\n"), "red");
			} else {
				if (!updateFITSKeyword(&gfit, old_keyname, new_val, NULL, NULL, TRUE, FALSE)) {
					gtk_list_store_set(key_liststore, &iter, COLUMN_KEY, new_val, -1);
				}
			}
		}
	}
	g_signal_handlers_unblock_by_func(treeview, G_CALLBACK(on_key_treeview_key_release_event), NULL);
}


void on_val_edited(GtkCellRendererText *renderer, char *path, char *new_val, gpointer user_data) {
	GtkWidget *treeview = GTK_WIDGET(user_data);
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(treeview));
	GtkTreeIter iter;
	gchar *FITS_key, *FITS_comment, *original_val;
	gboolean protected;
	char dtype;

	gtk_tree_model_get_iter_from_string(model, &iter, path);
	gtk_tree_model_get(model, &iter, COLUMN_KEY, &FITS_key, COLUMN_VALUE, &original_val, COLUMN_COMMENT, &FITS_comment, COLUMN_DTYPE, &dtype, COLUMN_PROTECTED, &protected, -1);
	if (!protected) {
		char valstring[FLEN_VALUE];
		/* update FITS key */
		process_keyword_string_value(new_val, valstring, dtype == 'C' && (new_val[0] != '\'' || new_val[strlen(new_val) - 1] != '\''));
		if (g_strcmp0(original_val, valstring)) {
			if (!updateFITSKeyword(&gfit, FITS_key, NULL, valstring, FITS_comment, TRUE, FALSE)) {
				gtk_list_store_set(key_liststore, &iter, COLUMN_VALUE, valstring, -1);
			}
		}
	}
	g_signal_handlers_unblock_by_func(treeview, G_CALLBACK(on_key_treeview_key_release_event), NULL);
}

void on_comment_edited(GtkCellRendererText *renderer, char *path, char *new_comment, gpointer user_data) {
	GtkWidget *treeview = GTK_WIDGET(user_data);
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(treeview));
	GtkTreeIter iter;
	gchar *FITS_key, *original_comment, *valstring;
	gboolean protected;
	char dtype;

	gtk_tree_model_get_iter_from_string(model, &iter, path);
	gtk_tree_model_get(model, &iter, COLUMN_KEY, &FITS_key, COLUMN_VALUE, &valstring, COLUMN_COMMENT, &original_comment, COLUMN_DTYPE, &dtype, COLUMN_PROTECTED, &protected, -1);
	if (!protected) {
		char commentstring[FLEN_COMMENT];
		/* update FITS comment */
		gsize len = g_strlcpy(commentstring, new_comment, FLEN_COMMENT);
		if (len >= FLEN_COMMENT) {
			siril_debug_print("Exceeded FITS COMMENT length\n");
		}
		if (g_strcmp0(original_comment, new_comment)) {
			if (!updateFITSKeyword(&gfit, FITS_key, NULL, valstring, commentstring, TRUE, FALSE)) {
				gtk_list_store_set(key_liststore, &iter, COLUMN_COMMENT, commentstring, -1);
			}
		}
	}
	g_signal_handlers_unblock_by_func(treeview, G_CALLBACK(on_key_treeview_key_release_event), NULL);
}

void on_key_close_btn_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("keywords_dialog");
}

static gchar *list_all_selected_keywords() {
	GtkTreeSelection *selection;
	GList *references, *list;
	GtkTreeModel *model = gtk_tree_view_get_model(key_treeview);
	GString *string = g_string_new("");

	selection = gtk_tree_view_get_selection(key_treeview);
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
	gtk_tree_selection_unselect_all(selection);
    g_list_free_full(references, (GDestroyNotify)gtk_tree_row_reference_free);

	return g_string_free(string, FALSE);
}

void on_key_selection_changed(GtkTreeSelection *selection, gpointer user_data) {
	GtkTreeIter iter;
	GList *list;

	GtkTreeModel *model = gtk_tree_view_get_model(key_treeview);
	gboolean is_empty = gtk_tree_model_get_iter_first(model, &iter) == FALSE;
	gboolean are_selected = FALSE;

	if (!is_empty) {
		list = gtk_tree_selection_get_selected_rows(selection, &model);
		are_selected = g_list_length(list)  > 0;
	}

	gtk_widget_set_sensitive(key_export_button, !is_empty && are_selected);
}

static void on_entry_comment_changed(GtkEntry *entry, gpointer user_data) {
	GtkEntry *entry_value = GTK_ENTRY(user_data);
	GtkEntry *entry_comment = GTK_ENTRY(entry);
	gint max_length = 67; // FLEN_CARD - 8 - 2 - 3 - 1

	const gchar *value_text = gtk_entry_get_text(entry_value);
	const gchar *comment_text = gtk_entry_get_text(entry_comment);

	/* In fact, there may be an error of 2 for strings in value,
	 when the two ' are automatically added.
	 We leave it like that, I think it's not too serious because we have
	 a mechanism that will truncate the character string during recording.
	 */
	gint total_length = strlen(value_text) + strlen(comment_text);

	if (total_length >= max_length) {
		gint allowed_length = max_length - strlen(value_text);
		gchar *new_text = g_strndup(comment_text, allowed_length);
		g_signal_handlers_block_by_func(entry, G_CALLBACK(on_entry_comment_changed), user_data);
		gtk_entry_set_text(entry_comment, new_text);
		g_signal_handlers_unblock_by_func(entry, G_CALLBACK(on_entry_comment_changed), user_data);
		g_free(new_text);
	}
}

static void on_entry_value_changed(GtkEntry *entry, gpointer user_data) {
	GtkEntry *entry_comment = GTK_ENTRY(user_data);
	GtkEntry *entry_value = GTK_ENTRY(entry);
	gint max_length = 67; // FLEN_CARD - 8 - 2 - 3 - 1

	const gchar *value_text = gtk_entry_get_text(entry_value);
	const gchar *comment_text = gtk_entry_get_text(entry_comment);

	/* In fact, there may be an error of 2 for strings in value,
	 when the two ' are automatically added.
	 We leave it like that, I think it's not too serious because we have
	 a mechanism that will truncate the character string during recording.
	 */
	gint total_length = strlen(value_text) + strlen(comment_text);

	if (total_length >= max_length) {
		gint allowed_length = max_length - strlen(comment_text);
		gchar *new_text = g_strndup(value_text, allowed_length);
		g_signal_handlers_block_by_func(entry, G_CALLBACK(on_entry_value_changed), user_data);
		gtk_entry_set_text(entry_value, new_text);
		g_signal_handlers_unblock_by_func(entry, G_CALLBACK(on_entry_value_changed), user_data);
		g_free(new_text);
	}
}

static void insert_text_handler(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	gchar *result = replace_wide_char(text);

	g_signal_handlers_block_by_func(G_OBJECT (editable), G_CALLBACK (insert_text_handler), data);
	gtk_editable_insert_text(editable, result, strlen(result), position);
	g_signal_handlers_unblock_by_func(G_OBJECT (editable), G_CALLBACK (insert_text_handler), data);

	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

static void insert_text_handler_key(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);

	GString *filtered_text = g_string_new(NULL);
	for (gint i = 0; i < length; i++) {
		// Append only characters that are not spaces
		if (!g_unichar_isspace(g_utf8_get_char(text + i))) {
			g_string_append_c(filtered_text, text[i]);
		}
	}

	gchar *result = replace_wide_char(filtered_text->str);

	g_signal_handlers_block_by_func(G_OBJECT(editable), G_CALLBACK(insert_text_handler_key), data);
	gtk_editable_insert_text(editable, result, strlen(result), position);
	g_signal_handlers_unblock_by_func(G_OBJECT(editable), G_CALLBACK(insert_text_handler_key), data);

	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
	g_string_free(filtered_text, TRUE);
}

static void on_entry_activate(GtkEntry *entry, gpointer user_data) {
	GtkWidget *add_button = GTK_WIDGET(user_data);
	gtk_button_clicked(GTK_BUTTON(add_button));
}

static void scroll_to_end() {
	GtkTreePath *path;
	GtkTreeIter iter;

	GtkTreeModel *model = gtk_tree_view_get_model(key_treeview);
	gint rows = gtk_tree_model_iter_n_children(model, NULL);

	if (rows > 0) {
		path = gtk_tree_path_new_from_indices(rows - 1, -1);

		gtk_tree_model_get_iter(model, &iter, path);
		gtk_tree_view_scroll_to_cell(key_treeview, path, NULL, FALSE, 0.0, 0.0);

		gtk_tree_path_free(path);
	}
}

void on_add_keyword_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *grid;
	GtkWidget *label_name;
	GtkWidget *label_value;
	GtkWidget *label_comment;
	GtkWidget *entry_name;
	GtkWidget *entry_value;
	GtkWidget *entry_comment;
	GtkWidget *add_button;
	GtkDialogFlags flags = GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT;

	// Create the dialog window with buttons
	dialog = gtk_dialog_new_with_buttons(_("Add New Keyword"),
			GTK_WINDOW(user_data), flags, _("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Add"), GTK_RESPONSE_OK,
			NULL);

	add_button = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_OK);
	gtk_widget_grab_focus(add_button);

	// Add the suggested-action style class to the Add button
	GtkStyleContext *context = gtk_widget_get_style_context(add_button);
	gtk_style_context_add_class(context, "suggested-action");

	// Set the dialog to be non-resizable
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);

	// Get the content area of the dialog
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

	// Create a grid to hold the labels and entries
	grid = gtk_grid_new();
	gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
	gtk_grid_set_row_spacing(GTK_GRID(grid), 5); // Set row spacing to 5
	gtk_grid_set_column_spacing(GTK_GRID(grid), 5); // Set column spacing to 5
	gtk_widget_set_margin_bottom(grid, 10); // Set bottom margin to 10 pixels
	gtk_container_add(GTK_CONTAINER(content_area), grid);

	// Create the Name label and entry
	label_name = gtk_label_new("Name:");
	gtk_widget_set_halign(label_name, GTK_ALIGN_START); // Align label to the start (left)
	entry_name = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(entry_name), 8); // Set the max length to 8. Don't want HIERARCH convention
	gtk_grid_attach(GTK_GRID(grid), label_name, 0, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_name, 1, 0, 1, 1);

	// Create the Value label and entry
	label_value = gtk_label_new("Value:");
	gtk_widget_set_halign(label_value, GTK_ALIGN_START); // Align label to the start (left)
	entry_value = gtk_entry_new();
	gtk_grid_attach(GTK_GRID(grid), label_value, 0, 1, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_value, 1, 1, 1, 1);

	// Create the Comment label and entry
	label_comment = gtk_label_new("Comment:");
	gtk_widget_set_halign(label_comment, GTK_ALIGN_START); // Align label to the start (left)
	entry_comment = gtk_entry_new();
	gtk_grid_attach(GTK_GRID(grid), label_comment, 0, 2, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_comment, 1, 2, 1, 1);

	// Set tooltips for each entry
	gtk_widget_set_tooltip_text(entry_name, "Enter the name of the keyword. Maximum 8 characters. Only ASCII characters are accepted.");
	gtk_widget_set_tooltip_text(entry_value, "Enter the value for the keyword. Only ASCII characters are accepted.");
	gtk_widget_set_tooltip_text(entry_comment, "Enter a comment or description for the keyword. Only ASCII characters are accepted.");

	// Connect the changed signal for both entry_value and entry_comment
	g_signal_connect(entry_value, "changed", G_CALLBACK(on_entry_value_changed), entry_comment);
	g_signal_connect(entry_comment, "changed", G_CALLBACK(on_entry_comment_changed), entry_value);

	// Connect the activate signal of each entry to activate the Add button
	g_signal_connect(entry_name, "activate", G_CALLBACK(on_entry_activate), add_button);
	g_signal_connect(entry_value, "activate", G_CALLBACK(on_entry_activate), add_button);
	g_signal_connect(entry_comment, "activate", G_CALLBACK(on_entry_activate), add_button);

	// Connect the insert-text signal of each entry
	g_signal_connect(entry_name, "insert-text", G_CALLBACK(insert_text_handler_key), add_button);
	g_signal_connect(entry_value, "insert-text", G_CALLBACK(insert_text_handler), add_button);
	g_signal_connect(entry_comment, "insert-text", G_CALLBACK(insert_text_handler), add_button);

	gtk_widget_show_all(dialog);

	gint result;
	do {
		result = gtk_dialog_run(GTK_DIALOG(dialog));
		if (result == GTK_RESPONSE_OK) {
			const gchar *key = gtk_entry_get_text(GTK_ENTRY(entry_name));
			const gchar *value = gtk_entry_get_text(GTK_ENTRY(entry_value));
			const gchar *comment = gtk_entry_get_text(GTK_ENTRY(entry_comment));

			if (g_strcmp0(key, "") == 0) key = NULL;
			if (g_strcmp0(value, "") == 0) value = NULL;
			if (g_strcmp0(comment, "") == 0) comment = NULL;

			if (comment || (value && key)) {
				char valstring[FLEN_VALUE] = { 0 };
				process_keyword_string_value(value, valstring, string_has_space(value));

				if (sequence_is_loaded()) {
					struct keywords_data *kargs = calloc(1, sizeof(struct keywords_data));

					kargs->FITS_key = g_strdup(key);
					kargs->value = valstring[0] == '\0' ? NULL : g_strdup(valstring);
					kargs->comment = g_strdup(comment);

					if (siril_confirm_dialog(_("Operation on the sequence"),
							_("These keywords will be added / modified in each image of "
							  "the entire sequence. Are you sure?"), _("Proceed"))) {
						start_sequence_keywords(&com.seq, kargs);
						break;
					} else {
						free(kargs);
						continue;
					}
				} else {
					updateFITSKeyword(&gfit, key, NULL, valstring[0] == '\0' ? NULL : valstring, comment, TRUE, FALSE);
					refresh_keywords_dialog();
					scroll_to_end();
					break;
				}
			}
		} else {
			break;
		}
	} while (TRUE);

	gtk_widget_destroy(dialog);
}

void on_delete_keyword_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_keys();
}

static void show_header_text(char *text) {
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(key_textview);
	GtkTextIter itDebut;
	GtkTextIter itFin;

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, text, strlen(text));
}

static void save_key_to_clipboard() {
	/* Get the clipboard object */
	GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
	GtkTextBuffer *buffer;
	GtkTextIter start, end;
	gchar *list = NULL;

	switch (gtk_notebook_get_current_page(key_notebook)) {
	case 0:
		list = list_all_selected_keywords();
		break;
	case 1:
		buffer = gtk_text_view_get_buffer(key_textview);
		if (gtk_text_buffer_get_selection_bounds(buffer, &start, &end)) {
			list = gtk_text_buffer_get_text(buffer, &start, &end, FALSE);
		}
		break;
	}
	if (list) {
		gtk_clipboard_set_text(clipboard, list, -1);
		g_free(list);
	}
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
	if (gfit.header)
		show_header_text(gfit.header);
}

void on_notebook_keywords_switch_page (GtkNotebook* self, GtkWidget* page, guint page_num, gpointer user_data) {
	GtkWidget *button = GTK_WIDGET(user_data);
    GtkTreeSelection *selection;
    GtkTextBuffer *buffer;
    GtkTextIter start, end;

	switch(page_num) {
	case 0:
	    selection = gtk_tree_view_get_selection(key_treeview);

		gtk_widget_set_sensitive(button, gtk_tree_selection_count_selected_rows(selection) > 0);
		break;
	case 1:
	    buffer = gtk_text_view_get_buffer(key_textview);

		gtk_widget_set_sensitive(button, gtk_text_buffer_get_selection_bounds(buffer, &start, &end));
		break;
	}
}
