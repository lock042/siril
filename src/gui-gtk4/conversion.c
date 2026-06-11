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

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <string.h>

#include "gui-gtk4/conversion.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "algos/sorting.h"
#include "registration/mpp/mpp_config.h"   /* enum mpp_avi_bayer */

static char *destroot = NULL;

/* ---- Per-row GObject for the convert list ----------------------- */
#define SIRIL_TYPE_CONVERT_ROW (siril_convert_row_get_type())
G_DECLARE_FINAL_TYPE(SirilConvertRow, siril_convert_row, SIRIL, CONVERT_ROW, GObject)
struct _SirilConvertRow {
	GObject parent_instance;
	gchar  *filename;       /* basename */
	gchar  *filename_full;  /* path */
	gchar  *size_str;       /* formatted size */
	gint64  size_int64;
	gchar  *date_str;       /* formatted date */
	guint64 date_unix;
};
G_DEFINE_TYPE(SirilConvertRow, siril_convert_row, G_TYPE_OBJECT)
static void siril_convert_row_finalize(GObject *obj) {
	SirilConvertRow *self = SIRIL_CONVERT_ROW(obj);
	g_clear_pointer(&self->filename,      g_free);
	g_clear_pointer(&self->filename_full, g_free);
	g_clear_pointer(&self->size_str,      g_free);
	g_clear_pointer(&self->date_str,      g_free);
	G_OBJECT_CLASS(siril_convert_row_parent_class)->finalize(obj);
}
static void siril_convert_row_class_init(SirilConvertRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_convert_row_finalize;
}
static void siril_convert_row_init(SirilConvertRow *self) { (void)self; }

/* GTK4: model, view, selection - all built programmatically. */
static GListStore *liststore_convert = NULL;
static GtkSortListModel *sortmodel_convert = NULL;
static GtkColumnView *columnview_convert = NULL;

/* Forward decls: defined further down; needed by ensure_convert_view(). */
static gboolean on_convert_drop(GtkDropTarget *target, const GValue *value,
		double x, double y, gpointer user_data);
static gboolean on_convert_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data);
static GtkMultiSelection *selection_convert = NULL;
static GtkScrolledWindow *scrolled_convert  = NULL;

static gboolean warning_is_displayed = FALSE;
static GtkCheckButton *conv_multiple_seq = NULL;
static GtkCheckButton *conv_demosaicing_btn = NULL;
static GtkCheckButton *conv_symlink = NULL;
static GtkEntry *conv_start_entry = NULL;
static GtkWidget *conv_go_button = NULL;
static GtkLabel *conv_status_label = NULL;
static GtkEntry *conv_root_entry = NULL;
static GtkComboBox *conv_avi_bayer_combo = NULL;
static GtkWidget *conv_avi_bayer_label = NULL;

static void conversion_init_statics(void) {
	if (conv_multiple_seq) return;
	conv_multiple_seq = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "multiple_seq"));
	conv_demosaicing_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "demosaicingButton"));
	conv_symlink = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "convert_symlink"));
	conv_start_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "startIndiceEntry"));
	conv_go_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "convert_button"));
	conv_status_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "statuslabel_convert"));
	conv_root_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "convroot_entry"));
	conv_avi_bayer_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_convert_avi_bayer"));
	conv_avi_bayer_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_convert_avi_bayer"));
	/* GTK4 GtkComboBoxText applies the .ui "active" before its <items> are
	 * added, leaving it blank; set the default (Auto) in code. */
	if (conv_avi_bayer_combo)
		gtk_combo_box_set_active(conv_avi_bayer_combo, 0);
}

static void check_for_conversion_form_completeness();
static void on_input_files_change();
static sequence_type get_activated_output_type();

/* ---- Cell factories --------------------------------------------- */

typedef enum {
	CONV_COL_NAME, CONV_COL_SIZE, CONV_COL_DATE
} ConvColKind;

static const gchar *convert_row_get_col(SirilConvertRow *r, ConvColKind k) {
	switch (k) {
		case CONV_COL_NAME: return r->filename;
		case CONV_COL_SIZE: return r->size_str;
		case CONV_COL_DATE: return r->date_str;
	}
	return "";
}

static void convert_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	/* Zero vertical margins on the cell child so the row hugs the text
	 * baseline; combined with the .siril-dense-rows CSS class on the
	 * columnview this matches the compactness of the file-open dialog. */
	gtk_widget_set_margin_top(lbl, 0);
	gtk_widget_set_margin_bottom(lbl, 0);
	gtk_list_item_set_child(li, lbl);
}

static void convert_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	ConvColKind kind = (ConvColKind)GPOINTER_TO_INT(u);
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilConvertRow *row = SIRIL_CONVERT_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, convert_row_get_col(row, kind));
}

static GtkListItemFactory *make_convert_factory(ConvColKind k) {
	GtkListItemFactory *f = gtk_signal_list_item_factory_new();
	g_signal_connect(f, "setup", G_CALLBACK(convert_setup_cb), NULL);
	g_signal_connect(f, "bind",  G_CALLBACK(convert_bind_cb),  GINT_TO_POINTER(k));
	return f;
}

/* ---- Sort / compare functions used by GtkSorter ----------------- */

static gint name_compare(gconstpointer a, gconstpointer b, gpointer user_data) {
	(void)user_data;
	const SirilConvertRow *ra = a;
	const SirilConvertRow *rb = b;
	gchar *k1 = g_utf8_collate_key_for_filename(ra->filename, -1);
	gchar *k2 = g_utf8_collate_key_for_filename(rb->filename, -1);
	gint r = g_strcmp0(k1, k2);
	g_free(k1); g_free(k2);
	return r;
}
static gint size_compare(gconstpointer a, gconstpointer b, gpointer user_data) {
	(void)user_data;
	const SirilConvertRow *ra = a;
	const SirilConvertRow *rb = b;
	if (ra->size_int64 < rb->size_int64) return -1;
	if (ra->size_int64 > rb->size_int64) return 1;
	return 0;
}
static gint date_compare(gconstpointer a, gconstpointer b, gpointer user_data) {
	(void)user_data;
	const SirilConvertRow *ra = a;
	const SirilConvertRow *rb = b;
	if (ra->date_unix < rb->date_unix) return -1;
	if (ra->date_unix > rb->date_unix) return 1;
	return 0;
}

/* ---- View creation ---------------------------------------------- */

static void on_selection_changed(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data) {
	(void)sm; (void)position; (void)n_items; (void)user_data;
	void update_statusbar_convert(void);
	update_statusbar_convert();
}

static void ensure_convert_view(void) {
	if (columnview_convert) return;
	if (!scrolled_convert)
		scrolled_convert = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolledwindow8"));

	liststore_convert = g_list_store_new(SIRIL_TYPE_CONVERT_ROW);
	columnview_convert = GTK_COLUMN_VIEW(gtk_column_view_new(NULL));
	/* Match the file-open dialog's tight row spacing (see siril.css). */
	gtk_widget_add_css_class(GTK_WIDGET(columnview_convert), "siril-dense-rows");

	GtkColumnViewColumn *cn = gtk_column_view_column_new(N_("Name"), make_convert_factory(CONV_COL_NAME));
	gtk_column_view_column_set_resizable(cn, TRUE);
	gtk_column_view_column_set_expand(cn, TRUE);
	gtk_column_view_column_set_sorter(cn,
			GTK_SORTER(gtk_custom_sorter_new(name_compare, NULL, NULL)));
	gtk_column_view_append_column(columnview_convert, cn);
	g_object_unref(cn);

	GtkColumnViewColumn *cs = gtk_column_view_column_new(N_("Size"), make_convert_factory(CONV_COL_SIZE));
	gtk_column_view_column_set_resizable(cs, TRUE);
	gtk_column_view_column_set_sorter(cs,
			GTK_SORTER(gtk_custom_sorter_new(size_compare, NULL, NULL)));
	gtk_column_view_append_column(columnview_convert, cs);
	g_object_unref(cs);

	GtkColumnViewColumn *cd = gtk_column_view_column_new(N_("Date"), make_convert_factory(CONV_COL_DATE));
	gtk_column_view_column_set_resizable(cd, TRUE);
	gtk_column_view_column_set_sorter(cd,
			GTK_SORTER(gtk_custom_sorter_new(date_compare, NULL, NULL)));
	gtk_column_view_append_column(columnview_convert, cd);
	g_object_unref(cd);

	GtkSorter *view_sorter = gtk_column_view_get_sorter(columnview_convert);
	sortmodel_convert = gtk_sort_list_model_new(G_LIST_MODEL(g_object_ref(liststore_convert)),
			view_sorter ? g_object_ref(view_sorter) : NULL);
	selection_convert = gtk_multi_selection_new(G_LIST_MODEL(g_object_ref(sortmodel_convert)));
	g_signal_connect(selection_convert, "selection-changed", G_CALLBACK(on_selection_changed), NULL);
	gtk_column_view_set_model(columnview_convert, GTK_SELECTION_MODEL(selection_convert));

	gtk_scrolled_window_set_child(scrolled_convert, GTK_WIDGET(columnview_convert));

	/* GTK4 event controllers replace the GTK3 "key-release-event" and
	 * "drag-data-received" signals on the legacy GtkTreeView. */
	GtkEventController *key_ctl = gtk_event_controller_key_new();
	g_signal_connect(key_ctl, "key-pressed",
	                 G_CALLBACK(on_convert_key_pressed), NULL);
	gtk_widget_add_controller(GTK_WIDGET(columnview_convert), key_ctl);

	GtkDropTarget *drop = gtk_drop_target_new(GDK_TYPE_FILE_LIST, GDK_ACTION_COPY);
	g_signal_connect(drop, "drop", G_CALLBACK(on_convert_drop), NULL);
	gtk_widget_add_controller(GTK_WIDGET(columnview_convert),
	                          GTK_EVENT_CONTROLLER(drop));
}

static void init_widgets() {
	ensure_convert_view();
	g_assert(liststore_convert);
}

/* Called once from the main-window startup path so the empty Convert
 * tree view is present from the start.  Previously ensure_convert_view
 * only ran when the user added a file, leaving the scrolled window
 * with no child until then. */
void conversion_tab_setup(void) {
	ensure_convert_view();
}

static void format_index_convert(GtkEntry *entry) {
	int idx = g_ascii_strtoull(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL, 10);
	gchar *str = NULL;

	if (idx < 1) {
		str = g_strdup_printf("00001");
	} else if (idx < 10000) {
		str = g_strdup_printf("%05d", idx);
	} else if (idx > INDEX_MAX) {
		str = g_strdup_printf("65535");
	}
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
}

int count_converted_files() {
	init_widgets();
	return (int)g_list_model_get_n_items(G_LIST_MODEL(liststore_convert));
}

int count_selected_files() {
	init_widgets();
	if (!selection_convert) return 0;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(selection_convert));
	int n = (int)gtk_bitset_get_size(bs);
	gtk_bitset_unref(bs);
	return n;
}

static void initialize_convert() {
	gchar *file_data;
	GList *list = NULL;

	init_widgets();

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (file_name_has_invalid_chars(destroot)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid char"), _("Please remove invalid characters in the sequence name "
				"before trying to convert images into a new sequence again."));
		return;
	}

	if (g_file_test(destroot, G_FILE_TEST_EXISTS)) {
		char *title = siril_log_message(_("A file named %s already exists. "
				"Do you want to replace it?\n"), destroot);
		gboolean replace = siril_confirm_dialog(title, _("The file already exists. "
				"Replacing it will overwrite its contents."), _("Replace File"));
		if (!replace) return;
	}

	/* iterate in displayed (sorted) order so the output sequence numbering
	 * follows what the user sees */
	guint nrows = g_list_model_get_n_items(G_LIST_MODEL(sortmodel_convert));
	if (nrows == 0) return;	//The list is empty

	gboolean no_sequence_to_convert = TRUE;
	gboolean there_is_an_image = FALSE;
	gboolean there_is_an_xtrans = FALSE;
	gboolean there_is_a_film = FALSE;
	int count = 0;
	for (guint i = 0; i < nrows; i++) {
		SirilConvertRow *row = SIRIL_CONVERT_ROW(g_list_model_get_item(G_LIST_MODEL(sortmodel_convert), i));
		file_data = g_strdup(row->filename_full);
		g_object_unref(row);

		list = g_list_prepend(list, file_data);

		const char *src_ext = get_filename_ext(file_data);

		image_type type = get_type_for_extension(src_ext);
		if (type == TYPEAVI || type == TYPESER) {
			no_sequence_to_convert = FALSE;
			if (type == TYPEAVI)
				there_is_a_film = TRUE;
		}
		else if (type == TYPEUNDEF) {
			char *title = siril_log_message(_("Filetype is not supported, cannot convert: %s\n"), src_ext);
			gchar *msg = g_strdup_printf(_("File extension '%s' is not supported.\n"
				"Verify that you typed the extension correctly.\n"
				"If so, you may need to install third-party software to enable "
				"this file type conversion, look at the README file.\n"
				"If the file type you are trying to load is listed in supported "
				"formats, you may notify the developers that the extension you are "
				"trying to use should be recognized for this type."), src_ext);
			siril_message_dialog(GTK_MESSAGE_ERROR, title, msg);
			g_free(msg);
			g_list_free_full(list, g_free);
			return;
		}
		else if (type == TYPERAW && !g_ascii_strcasecmp(src_ext, "raf")) {
			there_is_an_xtrans = TRUE;
			there_is_an_image = TRUE;
		}
		// because of fitseq, we can't use this check for FITS
		else if (type != TYPEFITS)
			there_is_an_image = TRUE;
		count++;
	}

	sequence_type output_type = get_activated_output_type();
	int nb_allowed;
	if (!allow_to_open_files(count, &nb_allowed) && output_type == SEQ_REGULAR) {
		gboolean confirm = siril_confirm_dialog(_("Too many files are being converted."),
				_("You are about to convert a large amount of files into standard FITS files. "
				"However, your OS limits the number of files that will be processed in the same time."
				"You may want to convert your input files into a FITS sequence."), _("Convert to FITS Sequence"));
		if (!confirm) return;
		output_type = SEQ_FITSEQ;
	}

	conversion_init_statics();
	gboolean multiple = siril_toggle_get_active(GTK_WIDGET(conv_multiple_seq));
	gboolean debayer = siril_toggle_get_active(GTK_WIDGET(conv_demosaicing_btn));
	gboolean symbolic_link = siril_toggle_get_active(GTK_WIDGET(conv_symlink));

	if (output_type == SEQ_REGULAR && debayer && symbolic_link) {
		siril_log_message(_("Symbolic links cannot be used when demosaicing the images, new images will be created\n"));
		symbolic_link = FALSE;
	}
	if (multiple && there_is_an_image) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("Creating multiple sequences can only be done with only sequences as input."));
		g_list_free_full(list, g_free);
		return;
	}
	if (output_type == SEQ_SER && there_is_an_xtrans && !debayer) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("FujiFilm XTRANS sensors are not supported by SER v2 (CFA-style) standard. You may use FITS sequences instead."));
		g_list_free_full(list, g_free);
		return;
	}

	siril_log_info(_("Conversion: processing %d files...\n"), count);

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	list = g_list_reverse(list);
	gchar **files_to_convert = glist_to_array(list, &count);

	struct _convert_data *args = calloc(1, sizeof(struct _convert_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		g_strfreev(files_to_convert);
		return;
	}
	if (output_type == SEQ_REGULAR) {
		conversion_init_statics();
		GtkEntry *startEntry = conv_start_entry;
		format_index_convert(startEntry);
		const gchar *index = gtk_editable_get_text(GTK_EDITABLE(startEntry));
		args->start = (g_ascii_strtoll(index, NULL, 10) <= 0
						|| g_ascii_strtoll(index, NULL, 10) >= INDEX_MAX) ?	1 : g_ascii_strtoll(index, NULL, 10);
	}
	else args->start = 0;
	args->list = files_to_convert;
	args->total = count;
	args->nb_converted_files = 0;
	args->input_has_a_seq = !no_sequence_to_convert;
	args->input_has_a_film = there_is_a_film;
	args->destroot = strdup(destroot);
	args->debayer = debayer;
	args->make_link = symbolic_link;
	args->output_type = output_type;
	args->multiple_output = multiple;
	{
		const int ab = (there_is_a_film && conv_avi_bayer_combo)
		             ? gtk_combo_box_get_active(conv_avi_bayer_combo)
		             : MPP_AVI_BAYER_AUTO;
		args->avi_bayer_pattern = (ab >= MPP_AVI_BAYER_AUTO && ab <= MPP_AVI_BAYER_GRBG)
		                       ? ab : MPP_AVI_BAYER_AUTO;
	}
	gettimeofday(&(args->t_start), NULL);
	if (!start_in_new_thread(convert_thread_worker, args)) {
		g_strfreev(args->list);
		g_free(args->destroot);
	}
	return;
}

void on_conv_root_entry_activate(GtkEntry *entry, gpointer user_data) {
	(void)entry; (void)user_data;
	initialize_convert();
}

void on_convert_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	initialize_convert();
}

static void add_file_to_list(GFile *file) {
	GFileInfo *info = g_file_query_info(file, G_FILE_ATTRIBUTE_TIME_MODIFIED ","
			G_FILE_ATTRIBUTE_STANDARD_SIZE, 0, NULL, NULL);
	guint64 mtime = g_file_info_get_attribute_uint64(info, G_FILE_ATTRIBUTE_TIME_MODIFIED);
	GDateTime *dt = g_date_time_new_from_unix_local(mtime);
	gchar *date = g_date_time_format(dt, "%c");
	const gchar *filename = g_file_peek_path(file);
	gchar *bname = g_file_get_basename(file);
	gchar *size = g_format_size(g_file_info_get_size(info));

	SirilConvertRow *row = g_object_new(SIRIL_TYPE_CONVERT_ROW, NULL);
	row->filename      = bname;          /* takes ownership */
	row->filename_full = g_strdup(filename);
	row->size_str      = size;           /* takes ownership */
	row->size_int64    = g_file_info_get_size(info);
	row->date_str      = date;           /* takes ownership */
	row->date_unix     = mtime;
	g_list_store_append(liststore_convert, row);
	g_object_unref(row);

	g_object_unref(info);
	g_date_time_unref(dt);
}

static void remove_selected_files_from_list() {
	init_widgets();
	if (!selection_convert) return;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(selection_convert));
	guint64 n = gtk_bitset_get_size(bs);
	if (n == 0) { gtk_bitset_unref(bs); return; }

	SirilConvertRow **rows = calloc(n, sizeof(SirilConvertRow *));
	for (guint64 i = 0; i < n; i++) {
		guint pos = gtk_bitset_get_nth(bs, (guint)i);
		rows[i] = SIRIL_CONVERT_ROW(g_list_model_get_item(G_LIST_MODEL(sortmodel_convert), pos));
	}
	gtk_bitset_unref(bs);
	for (guint64 i = 0; i < n; i++) {
		guint store_pos;
		if (rows[i] && g_list_store_find(liststore_convert, rows[i], &store_pos))
			g_list_store_remove(liststore_convert, store_pos);
		if (rows[i])
			g_object_unref(rows[i]);
	}
	free(rows);
	gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(selection_convert));
}

void fill_convert_list(GSList *list) {
	GSList *l;
	init_widgets();

	for (l = list; l; l = l->next) {
		char *filename;

		filename = (char *) l->data;
		GFile *file = g_file_new_for_path(filename);

		add_file_to_list(file);
		g_free(filename);
		g_object_unref(file);
	}
	check_for_conversion_form_completeness();
	on_input_files_change();
}

void on_clear_convert_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	init_widgets();
	g_list_store_remove_all(liststore_convert);
	check_for_conversion_form_completeness();
	on_input_files_change();
}

void on_remove_convert_button_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	init_widgets();
	remove_selected_files_from_list();
	check_for_conversion_form_completeness();
	on_input_files_change();
}

/* GTK4 drag-and-drop receiver: replaces the GTK3 "drag-data-received"
 * signal.  GtkDropTarget delivers a typed value (GdkFileList for file
 * drops); we extract GFile entries, filter by extension, and feed them
 * into fill_convert_list() exactly like the legacy handler did. */
static gboolean on_convert_drop(GtkDropTarget *target, const GValue *value,
		double x, double y, gpointer user_data) {
	(void)target; (void)x; (void)y; (void)user_data;
	if (!G_VALUE_HOLDS(value, GDK_TYPE_FILE_LIST))
		return FALSE;
	const GSList *files = g_value_get_boxed(value);
	GSList *list = NULL;
	gint bad_files = 0;
	for (const GSList *l = files; l; l = l->next) {
		GFile *f = l->data;
		gchar *path = g_file_get_path(f);
		if (!path) { bad_files++; continue; }
		const char *src_ext = get_filename_ext(path);
		if (!src_ext || get_type_for_extension(src_ext) == TYPEUNDEF) {
			bad_files++;
			g_free(path);
		} else {
			list = g_slist_prepend(list, path);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	fill_convert_list(list);
	if (bad_files) {
		gchar *loc_str = ngettext("%d file was ignored while drag and drop\n",
				"%d files were ignored while drag and drop\n", bad_files);
		loc_str = g_strdup_printf(loc_str, bad_files);
		char *msg = siril_log_message(loc_str);
		siril_message_dialog(GTK_MESSAGE_INFO, msg,
				_("Files with unknown extension cannot be dropped in this area. "
						"Therefore they are ignored."));
		g_free(loc_str);
	}
	g_slist_free(list);
	on_input_files_change();
	return TRUE;
}

/* GTK4 key handler: replaces the GTK3 "key-release-event".  Delete /
 * BackSpace removes selected rows. */
static gboolean on_convert_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data) {
	(void)controller; (void)keycode; (void)state; (void)user_data;
	if (keyval == GDK_KEY_Delete || keyval == GDK_KEY_KP_Delete
			|| keyval == GDK_KEY_BackSpace) {
		remove_selected_files_from_list();
		check_for_conversion_form_completeness();
		on_input_files_change();
		return TRUE;
	}
	return FALSE;
}


static void check_for_conversion_form_completeness() {
	conversion_init_statics();
	init_widgets();
	gboolean valid = (g_list_model_get_n_items(G_LIST_MODEL(liststore_convert)) > 0);
	gtk_widget_set_sensitive(conv_go_button, destroot && destroot[0] != '\0' && valid);
}

/* Scan the input list for AVI / film entries so the AVI-Bayer combo can be
 * shown only when relevant. We can't cheaply distinguish mono-encoded from
 * RGB-encoded AVIs without opening each file, so the visibility heuristic is
 * just "any AVI present"; the per-frame apply step in io/conversion.c skips
 * the pattern when the autodetected fit happens to be 3-channel colour. */
static gboolean input_list_has_any_film(void) {
	if (!liststore_convert) return FALSE;
	const guint n = g_list_model_get_n_items(G_LIST_MODEL(liststore_convert));
	for (guint i = 0; i < n; i++) {
		SirilConvertRow *row = SIRIL_CONVERT_ROW(g_list_model_get_item(G_LIST_MODEL(liststore_convert), i));
		const char *src_ext = row->filename_full ? get_filename_ext(row->filename_full) : NULL;
		const image_type type = src_ext ? get_type_for_extension(src_ext) : TYPEUNDEF;
		g_object_unref(row);
		if (type == TYPEAVI) return TRUE;
	}
	return FALSE;
}

static void update_avi_bayer_combo_visibility(void) {
	conversion_init_statics();
	if (!conv_avi_bayer_combo || !conv_avi_bayer_label) return;
	const gboolean show = input_list_has_any_film();
	gtk_widget_set_visible(GTK_WIDGET(conv_avi_bayer_combo), show);
	gtk_widget_set_visible(conv_avi_bayer_label, show);
}

static void on_input_files_change() {
	/* Sort functions are now per-column GtkSorters wired in
	 * ensure_convert_view(); nothing else to do here besides
	 * refreshing the status. */
	update_statusbar_convert();
	update_avi_bayer_combo_visibility();
}

/******************************* Callback functions ***********************************/

// TODO: put a red lining around the entry instead of removing bad chars
void insert_text_handler(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (is_forbiden_in_filename(text[i]))
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
		widget_set_class(GTK_WIDGET(entry), "", "warning");
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

void update_statusbar_convert() {
	conversion_init_statics();

	int nb_files = count_converted_files();
	if (nb_files == 0)
		gtk_label_set_text(conv_status_label, " ");
	else {
		int selected = count_selected_files();
		gchar *str1, *total;

		str1 = ngettext("%d file loaded", "%d files loaded", nb_files);
		str1 = g_strdup_printf(str1, nb_files);
		if (selected == 0) {
			total = g_strdup(str1);
		} else {
			gchar *str2 = ngettext("%d file selected", "%d files selected", selected);
			str2 = g_strdup_printf(str2, selected);
			total = g_strdup_printf("%s, %s", str1, str2);
			g_free(str2);
		}
		gtk_label_set_text(conv_status_label, total);
		g_free(str1);
		g_free(total);
	}
}

void on_treeview_selection_convert_changed(GObject *sm,
		gpointer user_data) {
	(void)sm; (void)user_data;
	update_statusbar_convert();
}

void process_destroot(sequence_type output_type) {
	conversion_init_statics();
	const gchar *name = gtk_editable_get_text(GTK_EDITABLE(conv_root_entry));
	if (*name == '\0') {
		free(destroot);
		destroot = NULL;
		return;
	}
	if (destroot) {
		free(destroot);
		destroot = NULL;
	}

	destroot = strdup(name); // Need to ensure destroot is always allocated
							 // from the stdlib memory pool not the glib one
	gboolean seq_exists = FALSE;
	if (output_type == SEQ_SER) {
		if (!g_str_has_suffix(destroot, ".ser")) {
			destroot = str_append(&destroot, ".ser");
		}
		seq_exists = check_if_seq_exist(destroot, FALSE);
	}
	else if (output_type == SEQ_FITSEQ) {
		if (!g_str_has_suffix(destroot, com.pref.ext)) {
			destroot = str_append(&destroot, com.pref.ext);
		}
		seq_exists = check_if_seq_exist(destroot, FALSE);
	}
	else {
		if (g_str_has_suffix(destroot, ".seq")) {
			char* temp = remove_ext_from_filename(destroot);
			if (!temp) {
				PRINT_ALLOC_ERR;
				return;
			}
			free(destroot);
			destroot = temp;
		}
		seq_exists = check_if_seq_exist(destroot, TRUE);
	}

	if (seq_exists && !warning_is_displayed) {
		set_icon_entry(conv_root_entry, "dialog-warning");
		warning_is_displayed = TRUE;
	}
	else if (!seq_exists && warning_is_displayed) {
		set_icon_entry(conv_root_entry, NULL);
		warning_is_displayed = FALSE;
	}
}

/* 0: FITS images
 * 1: SER sequence
 * 2: FITS sequence
 */
static sequence_type get_activated_output_type() {
	static GtkDropDown *combo = NULL;
	if (!combo)
		combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "prepro_output_type_combo1"));
	return (sequence_type)gtk_drop_down_get_selected(combo);
}

void on_convtoroot_changed(GtkEditable *editable, gpointer user_data){
	(void)editable; (void)user_data;
	process_destroot(get_activated_output_type());
	check_for_conversion_form_completeness();
}

void on_demosaicing_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	(void)user_data;
	com.pref.debayer.open_debayer = siril_toggle_get_active(GTK_WIDGET(togglebutton));
}

void on_prepro_output_type_combo1_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec; (void)user_data;
	conversion_init_statics();
	sequence_type output = gtk_drop_down_get_selected(combo);
	gboolean seqfile_output = output == SEQ_SER || output == SEQ_FITSEQ;
	gtk_widget_set_visible(GTK_WIDGET(conv_multiple_seq), seqfile_output);
	gtk_widget_set_visible(GTK_WIDGET(conv_start_entry), !seqfile_output);
	if (!seqfile_output)
		siril_toggle_set_active(GTK_WIDGET(conv_multiple_seq), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(conv_symlink), !seqfile_output);
	process_destroot(output);
	check_for_conversion_form_completeness();
}

void on_startIndiceEntry_activate(GtkEntry *entry, gpointer user_data) {
	(void)user_data;
	format_index_convert(entry);
}

