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
#include "core/proto.h"

#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "io/fits_keywords.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/utils.h"

#include "keywords_tree.h"

/* ---- Per-row GObject for the keywords list ----------------------- */
#define SIRIL_TYPE_KEYWORD_ROW (siril_keyword_row_get_type())
G_DECLARE_FINAL_TYPE(SirilKeywordRow, siril_keyword_row, SIRIL, KEYWORD_ROW, GObject)
struct _SirilKeywordRow {
	GObject parent_instance;
	gchar  *key;
	gchar  *value;
	gchar  *comment;
	char    dtype;
	gboolean protected_flag;
	gboolean editable;
};
G_DEFINE_TYPE(SirilKeywordRow, siril_keyword_row, G_TYPE_OBJECT)
static void siril_keyword_row_finalize(GObject *obj) {
	SirilKeywordRow *self = SIRIL_KEYWORD_ROW(obj);
	g_clear_pointer(&self->key,     g_free);
	g_clear_pointer(&self->value,   g_free);
	g_clear_pointer(&self->comment, g_free);
	G_OBJECT_CLASS(siril_keyword_row_parent_class)->finalize(obj);
}
static void siril_keyword_row_class_init(SirilKeywordRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_keyword_row_finalize;
}
static void siril_keyword_row_init(SirilKeywordRow *self) { (void)self; }

static GListStore *key_store = NULL;
static GtkSortListModel *key_sortmodel = NULL;
static GtkSelectionModel *key_selection_model = NULL;
static GtkColumnView *key_columnview = NULL;
static GtkScrolledWindow *key_scrolled = NULL;
static GtkTextView *key_textview = NULL;
static GtkNotebook *key_notebook = NULL;
static GtkWidget *key_export_button = NULL;

/* Forward decl: defined further down, referenced by the columnview build path. */
static gboolean on_key_columnview_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data);

/* The sequence_is_loaded ? SINGLE : MULTI selection mode is dynamic.  We keep
 * a pointer to whichever is currently in use; ensure_keyword_view rebuilds
 * the selection when the mode changes. */
static gboolean key_use_single_selection = FALSE;

static void on_mark_set(GtkTextBuffer *buffer, GtkTextIter *iter, GtkTextMark *mark, gpointer user_data) {
	(void)iter; (void)mark; (void)user_data;
	GtkTextIter start, end;
	GtkWidget *widget = key_export_button;
	gtk_widget_set_sensitive(widget, gtk_text_buffer_get_selection_bounds(buffer, &start, &end));
}

/* Every column in this view is a GtkEditableLabel, which is focusable and
 * swallows the click that would otherwise select the row.  Because no cell
 * leaves a non-editable surface to click, the selection model stayed
 * permanently empty, so the Delete button and "Copy Selection" (both of which
 * read the selection) silently did nothing.  Selecting the row whenever any of
 * its cells gains focus restores both, and makes keyboard navigation select
 * too. */
static void on_key_cell_focus_enter(GtkEventControllerFocus *ctrl, gpointer user_data) {
	(void)ctrl;
	GtkListItem *li = GTK_LIST_ITEM(user_data);
	if (!key_selection_model || !li)
		return;
	guint pos = gtk_list_item_get_position(li);
	if (pos == GTK_INVALID_LIST_POSITION)
		return;
	gtk_selection_model_select_item(key_selection_model, pos, TRUE);
}

static void key_string_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_editable_label_new("");
	gtk_widget_set_halign(lbl, GTK_ALIGN_FILL);
	gtk_list_item_set_child(li, lbl);

	GtkEventController *foc = gtk_event_controller_focus_new();
	g_signal_connect(foc, "enter", G_CALLBACK(on_key_cell_focus_enter), li);
	gtk_widget_add_controller(lbl, foc);
}

typedef enum { KEY_COL_KEY, KEY_COL_VALUE, KEY_COL_COMMENT } KeyStringCol;

/* Commit on focus-out / activate of the GtkEditableLabel: we read the current
 * text and if it changed call updateFITSKeyword.  GtkEditableLabel emits the
 * "changed" signal continuously while the label is being edited — too noisy.
 * Instead listen for "editing" property going from TRUE to FALSE. */
typedef struct {
	GtkListItem *li;
	KeyStringCol kind;
	gulong editing_handler;
} KeyEditCtx;

static void on_key_edit_done(GtkEditableLabel *lbl, KeyEditCtx *ctx);

static void key_edit_ctx_destroy(gpointer data) {
	KeyEditCtx *ctx = (KeyEditCtx *)data;
	if (ctx->editing_handler && ctx->li) {
		GtkWidget *w = gtk_list_item_get_child(ctx->li);
		if (w) g_signal_handler_disconnect(w, ctx->editing_handler);
	}
	g_free(ctx);
}

static void on_key_edit_done(GtkEditableLabel *lbl, KeyEditCtx *ctx) {
	if (gtk_editable_label_get_editing(lbl))
		return;
	SirilKeywordRow *row = SIRIL_KEYWORD_ROW(gtk_list_item_get_item(ctx->li));
	if (!row || row->protected_flag || !row->editable) return;

	const gchar *new_val = gtk_editable_get_text(GTK_EDITABLE(lbl));
	switch (ctx->kind) {
		case KEY_COL_KEY:
			if (g_strcmp0(row->key, new_val) != 0) {
				if (strlen(new_val) > 8) {
					siril_log_error(_("Keyname can contain a maximum of 8 characters.\n"));
					gtk_editable_set_text(GTK_EDITABLE(lbl), row->key);
				} else if (!updateFITSKeyword(gfit, row->key, (gchar*)new_val, NULL, NULL, TRUE, FALSE)) {
					g_free(row->key);
					row->key = g_strdup(new_val);
				} else {
					gtk_editable_set_text(GTK_EDITABLE(lbl), row->key);
				}
			}
			break;
		case KEY_COL_VALUE: {
			char valstring[FLEN_VALUE];
			process_keyword_string_value(new_val, valstring,
					row->dtype == 'C' && (new_val[0] != '\'' || new_val[strlen(new_val)-1] != '\''));
			if (g_strcmp0(row->value, valstring) != 0) {
				if (!updateFITSKeyword(gfit, row->key, NULL, valstring, row->comment, TRUE, FALSE)) {
					g_free(row->value);
					row->value = g_strdup(valstring);
					gtk_editable_set_text(GTK_EDITABLE(lbl), valstring);
				} else {
					gtk_editable_set_text(GTK_EDITABLE(lbl), row->value);
				}
			}
			break;
		}
		case KEY_COL_COMMENT: {
			char commentstring[FLEN_COMMENT];
			gsize len = g_strlcpy(commentstring, new_val, FLEN_COMMENT);
			if (len >= FLEN_COMMENT)
				siril_log_debug("Exceeded FITS COMMENT length\n");
			if (g_strcmp0(row->comment, commentstring) != 0) {
				if (!updateFITSKeyword(gfit, row->key, NULL, row->value, commentstring, TRUE, FALSE)) {
					g_free(row->comment);
					row->comment = g_strdup(commentstring);
				} else {
					gtk_editable_set_text(GTK_EDITABLE(lbl), row->comment);
				}
			}
			break;
		}
	}
}

static void on_editing_property_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)pspec;
	on_key_edit_done(GTK_EDITABLE_LABEL(obj), (KeyEditCtx *)user_data);
}

static void key_string_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	KeyStringCol kind = (KeyStringCol)GPOINTER_TO_INT(u);
	GtkEditableLabel *lbl = GTK_EDITABLE_LABEL(gtk_list_item_get_child(li));
	SirilKeywordRow *row = SIRIL_KEYWORD_ROW(gtk_list_item_get_item(li));
	const gchar *text = "";
	if (row) {
		text = (kind == KEY_COL_KEY) ? row->key
		     : (kind == KEY_COL_VALUE) ? row->value
		     : row->comment;
		text = text ? text : "";
	}
	gtk_editable_set_text(GTK_EDITABLE(lbl), text);

	gboolean edit_ok = row && !row->protected_flag && row->editable;
	g_object_set(lbl, "editable", edit_ok, NULL);

	/* Protected cards are shown with their keyword name in salmon, as in the
	 * GTK3 build (which bound COLUMN_COLOR="salmon" to the Keyword column's
	 * foreground).  Labels are recycled by the factory, so toggle the class on
	 * every bind.  Only the Keyword column carries the colour. */
	if (kind == KEY_COL_KEY && row && row->protected_flag)
		gtk_widget_add_css_class(GTK_WIDGET(lbl), "siril-protected-keyword");
	else
		gtk_widget_remove_css_class(GTK_WIDGET(lbl), "siril-protected-keyword");

	KeyEditCtx *ctx = g_new0(KeyEditCtx, 1);
	ctx->li = li;
	ctx->kind = kind;
	ctx->editing_handler = g_signal_connect(lbl, "notify::editing",
			G_CALLBACK(on_editing_property_changed), ctx);
	g_object_set_data_full(G_OBJECT(lbl), "siril-edit-ctx", ctx, key_edit_ctx_destroy);
}

static void key_string_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkEditableLabel *lbl = GTK_EDITABLE_LABEL(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(lbl), "siril-edit-ctx", NULL);
}

static GtkListItemFactory *make_key_string_factory(KeyStringCol c) {
	GtkListItemFactory *f = gtk_signal_list_item_factory_new();
	g_signal_connect(f, "setup",  G_CALLBACK(key_string_setup_cb),  NULL);
	g_signal_connect(f, "bind",   G_CALLBACK(key_string_bind_cb),   GINT_TO_POINTER(c));
	g_signal_connect(f, "unbind", G_CALLBACK(key_string_unbind_cb), NULL);
	return f;
}

/* Sort comparators */
static gint keystr(const SirilKeywordRow *r, KeyStringCol k) {
	(void)r; (void)k;
	return 0;  /* unused: helper to silence warning if not needed */
}
static gint cmp_key(gconstpointer a, gconstpointer b, gpointer u) { (void)u; (void)keystr;
	return g_strcmp0(((const SirilKeywordRow*)a)->key, ((const SirilKeywordRow*)b)->key);
}
static gint cmp_val(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return g_strcmp0(((const SirilKeywordRow*)a)->value, ((const SirilKeywordRow*)b)->value);
}
static gint cmp_comment(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	return g_strcmp0(((const SirilKeywordRow*)a)->comment, ((const SirilKeywordRow*)b)->comment);
}

/* GtkSelectionModel selection-changed → toggle Export button */
static void on_key_selection_model_changed(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data) {
	(void)sm; (void)position; (void)n_items; (void)user_data;
	gboolean is_empty = (g_list_model_get_n_items(G_LIST_MODEL(key_store)) == 0);
	gboolean any_selected = FALSE;
	if (!is_empty && key_selection_model) {
		GtkBitset *bs = gtk_selection_model_get_selection(key_selection_model);
		any_selected = (gtk_bitset_get_size(bs) > 0);
		gtk_bitset_unref(bs);
	}
	gtk_widget_set_sensitive(key_export_button, !is_empty && any_selected);
}

static void rebuild_selection_model(void) {
	if (!key_columnview) return;
	if (key_selection_model) {
		g_signal_handlers_disconnect_by_func(key_selection_model, on_key_selection_model_changed, NULL);
		g_clear_object(&key_selection_model);
	}
	/* wrap the sorted model so header sorting is preserved across rebuilds */
	GListModel *base = G_LIST_MODEL(g_object_ref(key_sortmodel));
	if (key_use_single_selection) {
		key_selection_model = GTK_SELECTION_MODEL(gtk_single_selection_new(base));
		gtk_single_selection_set_can_unselect(GTK_SINGLE_SELECTION(key_selection_model), TRUE);
		gtk_single_selection_set_autoselect(GTK_SINGLE_SELECTION(key_selection_model), FALSE);
	} else {
		key_selection_model = GTK_SELECTION_MODEL(gtk_multi_selection_new(base));
	}
	g_signal_connect(key_selection_model, "selection-changed", G_CALLBACK(on_key_selection_model_changed), NULL);
	gtk_column_view_set_model(key_columnview, key_selection_model);
}

static void ensure_keyword_view(void) {
	if (key_columnview) {
		gboolean want_single = sequence_is_loaded();
		if (want_single != key_use_single_selection) {
			key_use_single_selection = want_single;
			rebuild_selection_model();
		}
		return;
	}

	if (!key_scrolled)
		key_scrolled = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolled_keywords"));
	if (!key_textview)
		key_textview = GTK_TEXT_VIEW(gtk_builder_get_object(gui.builder, "FITS_header_txt"));
	if (!key_notebook)
		key_notebook = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook-keywords"));
	if (!key_export_button)
		key_export_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "export_keywords_button"));

	if (!key_store)
		key_store = g_list_store_new(SIRIL_TYPE_KEYWORD_ROW);

	key_use_single_selection = sequence_is_loaded();

	key_columnview = GTK_COLUMN_VIEW(gtk_column_view_new(NULL));
	gtk_column_view_set_show_column_separators(key_columnview, TRUE);
	gtk_widget_add_css_class(GTK_WIDGET(key_columnview), "siril-dense-rows");

	GtkColumnViewColumn *c;
	c = gtk_column_view_column_new(N_("Keyword"), make_key_string_factory(KEY_COL_KEY));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_key, NULL, NULL)));
	gtk_column_view_append_column(key_columnview, c); g_object_unref(c);

	c = gtk_column_view_column_new(N_("Value"), make_key_string_factory(KEY_COL_VALUE));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_expand(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_val, NULL, NULL)));
	gtk_column_view_append_column(key_columnview, c); g_object_unref(c);

	c = gtk_column_view_column_new(N_("Comment"), make_key_string_factory(KEY_COL_COMMENT));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_expand(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp_comment, NULL, NULL)));
	gtk_column_view_append_column(key_columnview, c); g_object_unref(c);

	GtkSorter *view_sorter = gtk_column_view_get_sorter(key_columnview);
	key_sortmodel = gtk_sort_list_model_new(G_LIST_MODEL(g_object_ref(key_store)),
			view_sorter ? g_object_ref(view_sorter) : NULL);
	rebuild_selection_model();

	gtk_scrolled_window_set_child(key_scrolled, GTK_WIDGET(key_columnview));

	/* GTK4: replaces the GTK3 "key-release-event" signal that fired on
	 * the legacy GtkTreeView. */
	GtkEventController *key_ctl = gtk_event_controller_key_new();
	g_signal_connect(key_ctl, "key-pressed",
	                 G_CALLBACK(on_key_columnview_key_pressed), NULL);
	gtk_widget_add_controller(GTK_WIDGET(key_columnview), key_ctl);

	g_signal_connect(gtk_text_view_get_buffer(key_textview), "mark-set", G_CALLBACK(on_mark_set), NULL);
}

static void add_key_to_tree(const gchar *key, const gchar *value, const gchar *comment,
		const char dtype, gboolean protected_flag, gboolean editable) {
	SirilKeywordRow *row = g_object_new(SIRIL_TYPE_KEYWORD_ROW, NULL);
	row->key      = g_strdup(key);
	row->value    = g_strdup(value);
	row->comment  = g_strdup(comment);
	row->dtype    = dtype;
	row->protected_flag = protected_flag;
	row->editable = !protected_flag && editable;
	g_list_store_append(key_store, row);
	g_object_unref(row);
}

static void init_dialog() {
	ensure_keyword_view();
	g_list_store_remove_all(key_store);
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

	fits_get_hdrspace(fptr, &nkeys, NULL, &status);

	for (ii = 1; ii <= nkeys; ii++) {
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
		add_key_to_tree(keyname, value, comment, dtype, keyword_is_protected(card, fit), editable);
	}

	if (status == END_OF_FILE)
		status = 0;

	fits_movabs_hdu(tmpfit.fptr, 1, NULL, &status);

	if (status)
		fits_report_error(stderr, status);

	fits_movabs_hdu(fptr, 0, NULL, &status);
	fits_close_file(fptr, &status);
	clearfits(&tmpfit);

	free(memptr);

	return (status);
}

void on_keywords_dialog_show(GtkWidget *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	gui_function(refresh_keywords_dialog, NULL);
}

static void remove_selected_keys() {
	if (!key_selection_model || !key_store) return;

	if (sequence_is_loaded()) {
		guint pos = gtk_single_selection_get_selected(GTK_SINGLE_SELECTION(key_selection_model));
		if (pos == GTK_INVALID_LIST_POSITION) return;
		SirilKeywordRow *row = SIRIL_KEYWORD_ROW(g_list_model_get_item(G_LIST_MODEL(key_sortmodel), pos));
		if (!row) return;
		struct keywords_data *kargs = calloc(1, sizeof(struct keywords_data));
		kargs->FITS_key = g_strdup(row->key);
		kargs->value = NULL;
		kargs->comment = NULL;
		if (siril_confirm_dialog(_("Operation on the sequence"),
				_("These keywords will be deleted from each image of "
				  "the entire sequence. Are you sure?"), _("Proceed"))) {
			guint store_pos;
			if (g_list_store_find(key_store, row, &store_pos))
				g_list_store_remove(key_store, store_pos);
			start_sequence_keywords(&com.seq, kargs);
		} else {
			g_free(kargs->FITS_key);
			free(kargs);
		}
		g_object_unref(row);
	} else {
		GtkBitset *bs = gtk_selection_model_get_selection(key_selection_model);
		guint64 n = gtk_bitset_get_size(bs);
		if (n == 0) { gtk_bitset_unref(bs); return; }
		/* Selection positions index the sorted model: collect the row objects
		 * first, then remove each from the store by identity. */
		SirilKeywordRow **rows = calloc(n, sizeof(SirilKeywordRow *));
		for (guint64 i = 0; i < n; i++) {
			guint pos = gtk_bitset_get_nth(bs, (guint)i);
			rows[i] = SIRIL_KEYWORD_ROW(g_list_model_get_item(G_LIST_MODEL(key_sortmodel), pos));
		}
		gtk_bitset_unref(bs);
		for (guint64 i = 0; i < n; i++) {
			if (!rows[i]) continue;
			updateFITSKeyword(gfit, rows[i]->key, NULL, NULL, NULL, TRUE, FALSE);
			guint store_pos;
			if (g_list_store_find(key_store, rows[i], &store_pos))
				g_list_store_remove(key_store, store_pos);
			g_object_unref(rows[i]);
		}
		free(rows);
		gtk_selection_model_unselect_all(key_selection_model);
	}
}

/* GTK4 key handler.  Body: Delete / BackSpace removes the selected rows.
 * The GTK3 version was wired by .ui as a "key-release-event" signal;
 * GTK4 uses GtkEventControllerKey, attached to the GtkColumnView in
 * key_columnview_init_statics() above (the function is forward-declared
 * near the top of the file so the controller wiring can reference it). */
static gboolean on_key_columnview_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data) {
	(void)controller; (void)keycode; (void)state; (void)user_data;
	if (keyval == GDK_KEY_Delete || keyval == GDK_KEY_KP_Delete
			|| keyval == GDK_KEY_BackSpace) {
		remove_selected_keys();
		return TRUE;
	}
	return FALSE;
}

/* Stub legacy GtkCellRendererText handlers — superseded by the inline
 * GtkEditableLabel in the column view.  Kept as no-ops so any stale .ui
 * cache still binds successfully. */
void on_cell_editing_started(GtkWidget *renderer, GtkEditable *editable, const gchar *path, gpointer user_data) {
	(void)renderer; (void)editable; (void)path; (void)user_data;
}
void on_cell_editing_canceled(GtkWidget *renderer, gpointer user_data) {
	(void)renderer; (void)user_data;
}
void on_key_edited(GtkWidget *renderer, char *path, char *new_val, gpointer user_data) {
	(void)renderer; (void)path; (void)new_val; (void)user_data;
}
void on_val_edited(GtkWidget *renderer, char *path, char *new_val, gpointer user_data) {
	(void)renderer; (void)path; (void)new_val; (void)user_data;
}
void on_comment_edited(GtkWidget *renderer, char *path, char *new_val, gpointer user_data) {
	(void)renderer; (void)path; (void)new_val; (void)user_data;
}

void on_key_close_btn_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	siril_close_dialog("keywords_dialog");
}

static gchar *list_all_selected_keywords() {
	if (!key_selection_model || !key_store) return g_strdup("");
	GString *string = g_string_new("");
	GtkBitset *bs = gtk_selection_model_get_selection(key_selection_model);
	guint64 n = gtk_bitset_get_size(bs);
	for (guint64 i = 0; i < n; i++) {
		guint pos = gtk_bitset_get_nth(bs, (guint)i);
		SirilKeywordRow *row = SIRIL_KEYWORD_ROW(g_list_model_get_item(G_LIST_MODEL(key_sortmodel), pos));
		if (row) {
			int status = 0;
			char card[FLEN_CARD];
			fits_make_key(row->key, row->value, row->comment, card, &status);
			g_string_prepend(string, "\n");
			g_string_prepend(string, card);
			g_object_unref(row);
		}
	}
	gtk_bitset_unref(bs);
	gtk_selection_model_unselect_all(key_selection_model);
	return g_string_free(string, FALSE);
}

/* Legacy stub: the GtkTreeSelection signal no longer fires; the GtkSelectionModel
 * "selection-changed" handler on_key_selection_model_changed handles this now. */
void on_key_selection_changed(GtkWidget *selection, gpointer user_data) {
	(void)selection; (void)user_data;
}

static void on_entry_comment_changed(GtkEntry *entry, gpointer user_data) {
	GtkEntry *entry_value = GTK_ENTRY(user_data);
	GtkEntry *entry_comment = GTK_ENTRY(entry);
	gint max_length = 67;

	const gchar *value_text = gtk_editable_get_text(GTK_EDITABLE(entry_value));
	const gchar *comment_text = gtk_editable_get_text(GTK_EDITABLE(entry_comment));

	gint total_length = strlen(value_text) + strlen(comment_text);

	if (total_length >= max_length) {
		gint allowed_length = max_length - strlen(value_text);
		gchar *new_text = g_strndup(comment_text, allowed_length);
		g_signal_handlers_block_by_func(entry, G_CALLBACK(on_entry_comment_changed), user_data);
		gtk_editable_set_text(GTK_EDITABLE(entry_comment), new_text);
		g_signal_handlers_unblock_by_func(entry, G_CALLBACK(on_entry_comment_changed), user_data);
		g_free(new_text);
	}
}

static void on_entry_value_changed(GtkEntry *entry, gpointer user_data) {
	GtkEntry *entry_comment = GTK_ENTRY(user_data);
	GtkEntry *entry_value = GTK_ENTRY(entry);
	gint max_length = 67;

	const gchar *value_text = gtk_editable_get_text(GTK_EDITABLE(entry_value));
	const gchar *comment_text = gtk_editable_get_text(GTK_EDITABLE(entry_comment));

	gint total_length = strlen(value_text) + strlen(comment_text);

	if (total_length >= max_length) {
		gint allowed_length = max_length - strlen(comment_text);
		gchar *new_text = g_strndup(value_text, allowed_length);
		g_signal_handlers_block_by_func(entry, G_CALLBACK(on_entry_value_changed), user_data);
		gtk_editable_set_text(GTK_EDITABLE(entry_value), new_text);
		g_signal_handlers_unblock_by_func(entry, G_CALLBACK(on_entry_value_changed), user_data);
		g_free(new_text);
	}
}

static void insert_text_handler(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	(void)entry; (void)length;
	GtkEditable *editable = GTK_EDITABLE(entry);
	gchar *result = replace_wide_char(text);

	g_signal_handlers_block_by_func(G_OBJECT(editable), G_CALLBACK(insert_text_handler), data);
	gtk_editable_insert_text(editable, result, strlen(result), position);
	g_signal_handlers_unblock_by_func(G_OBJECT(editable), G_CALLBACK(insert_text_handler), data);

	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

static void insert_text_handler_key(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	(void)entry;
	GtkEditable *editable = GTK_EDITABLE(entry);

	GString *filtered_text = g_string_new(NULL);
	for (gint i = 0; i < length; i++) {
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
	(void)entry;
	GtkWidget *add_button = GTK_WIDGET(user_data);
	g_signal_emit_by_name(GTK_BUTTON(add_button), "clicked");
}

static void scroll_to_end() {
	if (!key_columnview || !key_store) return;
	guint n = g_list_model_get_n_items(G_LIST_MODEL(key_store));
	if (n > 0)
		gtk_column_view_scroll_to(key_columnview, n - 1, NULL,
				GTK_LIST_SCROLL_NONE, NULL);
}

/* Phase 14G.4 helpers: connect button "clicked" to a GMainLoop quit, and
 * record the response code so the synchronous caller can branch on it. */
struct kw_dlg_ctx { GMainLoop *loop; gint result; };
static void kw_dlg_set_result_ok(GtkButton *b, gpointer ud)     { (void)b; struct kw_dlg_ctx *c = ud; c->result = GTK_RESPONSE_OK;     if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); }
static void kw_dlg_set_result_cancel(GtkButton *b, gpointer ud) { (void)b; struct kw_dlg_ctx *c = ud; c->result = GTK_RESPONSE_CANCEL; if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); }
static gboolean kw_dlg_close_request(GtkWindow *w, gpointer ud) { (void)w; struct kw_dlg_ctx *c = ud; c->result = GTK_RESPONSE_CANCEL; if (c->loop && g_main_loop_is_running(c->loop)) g_main_loop_quit(c->loop); return FALSE; }

void on_add_keyword_button_clicked(GtkButton *button, gpointer user_data) {
	(void)user_data;
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

	/* Phase 14G.4: GtkDialog → GtkWindow.  We build a vertical content
	 * area + explicit Cancel/Add button row, and bridge the synchronous
	 * gtk_dialog_run loop with a nested GMainLoop. */
	dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Add New Keyword"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(dialog), TRUE);
	/* The signal in keywords_dialog.ui carries no "object" attribute, so the
	 * handler's user_data is NULL under GTK4 (there is no
	 * gtk_builder_connect_signals to supply a default).  Deriving the parent
	 * from the button's toplevel gives the FITS Header window as a proper
	 * transient parent; without one, a modal GtkWindow grabs input with no
	 * anchor and wedges the whole session on Wayland/Xorg. */
	gtk_window_set_transient_for(GTK_WINDOW(dialog),
			GTK_WINDOW(gtk_widget_get_root(GTK_WIDGET(button))));
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);

	content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_window_set_child(GTK_WINDOW(dialog), content_area);

	grid = gtk_grid_new();
	gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
	gtk_grid_set_row_spacing(GTK_GRID(grid), 5);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 5);
	gtk_widget_set_margin_bottom(grid, 10);
	gtk_box_append(GTK_BOX(content_area), grid);

	label_name = gtk_label_new("Name:");
	gtk_widget_set_halign(label_name, GTK_ALIGN_START);
	entry_name = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(entry_name), 8);
	gtk_grid_attach(GTK_GRID(grid), label_name, 0, 0, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_name, 1, 0, 1, 1);

	label_value = gtk_label_new("Value:");
	gtk_widget_set_halign(label_value, GTK_ALIGN_START);
	entry_value = gtk_entry_new();
	gtk_grid_attach(GTK_GRID(grid), label_value, 0, 1, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_value, 1, 1, 1, 1);

	label_comment = gtk_label_new("Comment:");
	gtk_widget_set_halign(label_comment, GTK_ALIGN_START);
	entry_comment = gtk_entry_new();
	gtk_grid_attach(GTK_GRID(grid), label_comment, 0, 2, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), entry_comment, 1, 2, 1, 1);

	gtk_widget_set_tooltip_text(entry_name, "Enter the name of the keyword. Maximum 8 characters. Only ASCII characters are accepted.");
	gtk_widget_set_tooltip_text(entry_value, "Enter the value for the keyword. Only ASCII characters are accepted.");
	gtk_widget_set_tooltip_text(entry_comment, "Enter a comment or description for the keyword. Only ASCII characters are accepted.");

	g_signal_connect(entry_value, "changed", G_CALLBACK(on_entry_value_changed), entry_comment);
	g_signal_connect(entry_comment, "changed", G_CALLBACK(on_entry_comment_changed), entry_value);

	/* Action row — create add_button first so the entry "activate" signals
	 * below can pass it as user_data without triggering uninitialized-use. */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	gtk_widget_set_margin_top(bbox, 6);
	GtkWidget *cancel_button = gtk_button_new_with_mnemonic(_("_Cancel"));
	add_button = gtk_button_new_with_mnemonic(_("_Add"));
	gtk_widget_add_css_class(add_button, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), cancel_button);
	gtk_box_append(GTK_BOX(bbox), add_button);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), add_button);
	gtk_widget_grab_focus(add_button);

	g_signal_connect(entry_name, "activate", G_CALLBACK(on_entry_activate), add_button);
	g_signal_connect(entry_value, "activate", G_CALLBACK(on_entry_activate), add_button);
	g_signal_connect(entry_comment, "activate", G_CALLBACK(on_entry_activate), add_button);

	g_signal_connect(entry_name, "insert-text", G_CALLBACK(insert_text_handler_key), add_button);
	g_signal_connect(entry_value, "insert-text", G_CALLBACK(insert_text_handler), add_button);
	g_signal_connect(entry_comment, "insert-text", G_CALLBACK(insert_text_handler), add_button);

	/* Sync wait via nested GMainLoop. */
	struct kw_dlg_ctx ctx = { g_main_loop_new(NULL, FALSE), GTK_RESPONSE_CANCEL };
	g_signal_connect(cancel_button, "clicked",        G_CALLBACK(kw_dlg_set_result_cancel), &ctx);
	g_signal_connect(add_button,    "clicked",        G_CALLBACK(kw_dlg_set_result_ok),     &ctx);
	g_signal_connect(dialog,        "close-request",  G_CALLBACK(kw_dlg_close_request),     &ctx);

	gtk_widget_set_visible(dialog, TRUE);

	gint result;
	do {
		ctx.result = GTK_RESPONSE_CANCEL;
		g_main_loop_run(ctx.loop);
		result = ctx.result;
		if (result == GTK_RESPONSE_OK) {
			const gchar *key = gtk_editable_get_text(GTK_EDITABLE(entry_name));
			const gchar *value = gtk_editable_get_text(GTK_EDITABLE(entry_value));
			const gchar *comment = gtk_editable_get_text(GTK_EDITABLE(entry_comment));

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
						g_free(kargs->FITS_key);
						g_free(kargs->comment);
						g_free(kargs->value);
						free(kargs);
						continue;
					}
				} else {
					updateFITSKeyword(gfit, key, NULL, valstring[0] == '\0' ? NULL : valstring, comment, TRUE, FALSE);
					gui_function(refresh_keywords_dialog, NULL);
					scroll_to_end();
					break;
				}
			}
		} else {
			break;
		}
	} while (TRUE);

	g_main_loop_unref(ctx.loop);
	gtk_window_destroy(GTK_WINDOW(dialog));
}

void on_delete_keyword_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
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
	GdkClipboard *clipboard = gdk_display_get_clipboard(gdk_display_get_default());
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
		gdk_clipboard_set_text(clipboard, list);
		g_free(list);
	}
}

void on_export_keywords_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	save_key_to_clipboard();
}

gboolean refresh_keywords_dialog(gpointer user_data) {
	(void)user_data;
	init_dialog();
	gboolean is_a_single_image_loaded = single_image_is_loaded() &&
			(!sequence_is_loaded() || (sequence_is_loaded() &&
			(com.seq.current == RESULT_IMAGE || com.seq.current == SCALED_IMAGE)));
	listFITSKeywords(gfit, is_a_single_image_loaded);
	if (gfit->header)
		show_header_text(gfit->header);
	return FALSE;
}

void on_notebook_keywords_switch_page(GtkNotebook* self, GtkWidget* page, guint page_num, gpointer user_data) {
	(void)self; (void)page;
	GtkWidget *button = GTK_WIDGET(user_data);
	GtkTextBuffer *buffer;
	GtkTextIter start, end;

	switch(page_num) {
	case 0: {
		gboolean any = FALSE;
		if (key_selection_model) {
			GtkBitset *bs = gtk_selection_model_get_selection(key_selection_model);
			any = (gtk_bitset_get_size(bs) > 0);
			gtk_bitset_unref(bs);
		}
		gtk_widget_set_sensitive(button, any);
		break;
	}
	case 1:
		buffer = gtk_text_view_get_buffer(key_textview);
		gtk_widget_set_sensitive(button, gtk_text_buffer_get_selection_bounds(buffer, &start, &end));
		break;
	}
}
