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

#ifdef _WIN32
#include <windows.h>
#endif
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/plot.h"
#include "gui-gtk4/registration.h"	// for update_reg_interface
#include "gui-gtk4/stacking.h"	// for update_stack_interface
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/PSF.h"
#include "registration/registration.h"
#include <dirent.h>

#include "sequence_list.h"

static gboolean fill_sequence_list_idle(gpointer p);

static const char *bg_colour[] = { "WhiteSmoke", "#1B1B1B" };
static const char *ref_bg_colour[] = { "Beige", "#4A4A39" };
static gchar *qualfmt = "%.3f";

static GtkWidget *combo = NULL, *arcsec = NULL, *sourceCombo = NULL;
static GtkDropDown *seqlist_dialog_combo = NULL;
static GtkHeaderBar *seqlist_headerbar = NULL;
static GtkWidget *seqlist_buttonbar = NULL;
static GtkWidget *seqlist_refframe2 = NULL;
/* GTK4: GtkSearchEntry no longer derives from GtkEntry — both implement
 * GtkEditable instead.  Hold a GtkEditable* and use the editable API. */
static GtkEditable *seqlist_search_entry = NULL;
static GtkScrolledWindow *seqlist_scrolled = NULL;
/* Phase 11: GTK4 model + view replacing liststore1 + GtkTreeView treeview1. */
static GListStore        *seq_store      = NULL;
static GtkColumnView     *seq_columnview = NULL;
static GtkMultiSelection *seq_selection  = NULL;
static GtkColumnViewColumn *seq_quality_col = NULL; /* dynamic-titled column */
static GtkWidget *seqlist_drawframe_check = NULL;
static GtkDropDown *seqlist_sequence_combobox = NULL;
static GtkLabel *seqlist_statusbar = NULL;

static void seq_build_columnview(void);

static void sequence_list_init_statics(void) {
	if (seq_columnview) return;
	seqlist_dialog_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "seqlist_dialog_combo"));
	seqlist_headerbar = GTK_HEADER_BAR(gtk_builder_get_object(gui.builder, "seqlistbar"));
#if defined(OS_OSX) && GTK_CHECK_VERSION(4, 18, 0)
	gtk_header_bar_set_use_native_controls(seqlist_headerbar, TRUE);
#endif
	seqlist_buttonbar = GTK_WIDGET(gtk_builder_get_object(gui.builder, "seqlist_buttonbar"));
	seqlist_refframe2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "refframe2"));
	seqlist_search_entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "seqlistsearch"));
	seqlist_scrolled = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolledwindow1"));
	seqlist_drawframe_check = GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawframe_check"));
	seqlist_sequence_combobox = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "sequence_list_combobox"));
	seqlist_statusbar = GTK_LABEL(gtk_builder_get_object(gui.builder, "seqlist_statusbar"));
	combo = GTK_WIDGET(gtk_builder_get_object(gui.builder, "plotCombo"));
	sourceCombo = GTK_WIDGET(gtk_builder_get_object(gui.builder, "plotSourceCombo"));
	arcsec = GTK_WIDGET(gtk_builder_get_object(gui.builder, "arcsecPhotometry"));
	seq_build_columnview();
}

static int selected_source = -1;
static gboolean is_arcsec = FALSE;
static gboolean use_photometry = FALSE;

static int cached_real_index = -1;
static gchar *cached_original_filename = NULL;

struct _seq_list {
	gpointer tview;  /* Phase 11: unused (legacy field, kept for ABI safety in case any other compilation unit happens to share this layout). */
	sequence *seq;
	int layer;
};

enum {
	COLUMN_IMNAME,		// string
	COLUMN_SHIFTX,		// int
	COLUMN_SHIFTY,		// int
	COLUMN_SELECTED,	// gboolean
	COLUMN_FWHM,		// gdouble
	COLUMN_CURRENT,		// int weight, current file loaded, display IMNAME in bold
	COLUMN_REFERENCE,	// background color depending on the image being reference
	COLUMN_INDEX,		// int
	N_COLUMNS
};

/* ── Phase 11: per-row GObject for the sequence list ─────────────────────
 * Mirrors the legacy 8-column GtkListStore.  Index is 1-based to match
 * the original wire format; real_index = index - 1 in the underlying
 * com.seq.imgparam[] array.
 */
#define SIRIL_TYPE_SEQ_ROW (siril_seq_row_get_type())
G_DECLARE_FINAL_TYPE(SirilSeqRow, siril_seq_row, SIRIL, SEQ_ROW, GObject)
struct _SirilSeqRow {
	GObject parent_instance;
	gchar  *imname;       /* basename */
	gint    shiftx;
	gint    shifty;
	gboolean selected;
	gdouble fwhm;         /* compared against -DBL_MAX for "N/A" */
	gint    weight;       /* Pango weight: 400 normal, 800 bold */
	gchar  *ref_bg;       /* background colour name */
	gint    index;        /* 1-based */
};
G_DEFINE_TYPE(SirilSeqRow, siril_seq_row, G_TYPE_OBJECT)
static void siril_seq_row_finalize(GObject *obj) {
	SirilSeqRow *self = SIRIL_SEQ_ROW(obj);
	g_clear_pointer(&self->imname, g_free);
	g_clear_pointer(&self->ref_bg, g_free);
	G_OBJECT_CLASS(siril_seq_row_parent_class)->finalize(obj);
}
static void siril_seq_row_class_init(SirilSeqRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_seq_row_finalize;
}
static void siril_seq_row_init(SirilSeqRow *self) { (void)self; }

/* Phase 12: GtkStatusbar replaced by GtkLabel; context-id no longer used. */

void on_ref_frame_toggled(GtkCheckButton *togglebutton, gpointer user_data);

/******* Static functions **************/

/* ── Phase 11: GtkColumnView cell factories ─────────────────────────────
 * Each visible column has a setup/bind pair.  "weight" (bold for current
 * row) is applied to the File column label via PangoAttrList.  "ref_bg"
 * (alternate background colour for the reference row) is applied via
 * a CSS class on the cell label using gtk_widget_add_css_class.
 */

typedef enum {
	SEQ_COL_INDEX, SEQ_COL_FILE, SEQ_COL_SHIFTX, SEQ_COL_SHIFTY,
	SEQ_COL_QUALITY
} SeqColKind;

static GtkCssProvider *seq_css_provider = NULL;

static void ensure_seq_css(void) {
	if (seq_css_provider) return;
	seq_css_provider = gtk_css_provider_new();
	gtk_css_provider_load_from_string(seq_css_provider,
			"label.seq-ref { background-color: alpha(@warning_color, 0.25); }"
	);
	gtk_style_context_add_provider_for_display(gdk_display_get_default(),
			GTK_STYLE_PROVIDER(seq_css_provider),
			GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
}

static void seq_setup_label_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}

static void seq_apply_ref_bg(GtkLabel *lbl, SirilSeqRow *row) {
	GtkWidget *w = GTK_WIDGET(lbl);
	gboolean is_ref = (row && row->index > 0 && (row->index - 1) == com.seq.reference_image);
	if (is_ref) {
		if (!gtk_widget_has_css_class(w, "seq-ref"))
			gtk_widget_add_css_class(w, "seq-ref");
	} else {
		gtk_widget_remove_css_class(w, "seq-ref");
	}
}

static void seq_bind_index_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gchar *txt = g_strdup_printf("%d", row->index);
	gtk_label_set_text(lbl, txt);
	g_free(txt);
	seq_apply_ref_bg(lbl, row);
}

static void seq_bind_file_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gtk_label_set_text(lbl, row->imname ? row->imname : "");
	PangoAttrList *attrs = pango_attr_list_new();
	pango_attr_list_insert(attrs, pango_attr_weight_new((PangoWeight)row->weight));
	gtk_label_set_attributes(lbl, attrs);
	pango_attr_list_unref(attrs);
	seq_apply_ref_bg(lbl, row);
}

static void seq_bind_shiftx_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gchar *txt = g_strdup_printf("%d", row->shiftx);
	gtk_label_set_text(lbl, txt);
	g_free(txt);
	seq_apply_ref_bg(lbl, row);
}

static void seq_bind_shifty_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gchar *txt = g_strdup_printf("%d", row->shifty);
	gtk_label_set_text(lbl, txt);
	g_free(txt);
	seq_apply_ref_bg(lbl, row);
}

static void seq_bind_quality_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	if (row->fwhm > -DBL_MAX) {
		gchar buf[20];
		g_snprintf(buf, sizeof(buf), qualfmt, row->fwhm);
		gtk_label_set_text(lbl, buf);
	} else {
		gtk_label_set_text(lbl, "N/A");
	}
	seq_apply_ref_bg(lbl, row);
}

/* Sel toggle column: inline GtkCheckButton bound to row->selected. */
static void on_seq_selected_toggled(GtkCheckButton *btn, gpointer user_data) {
	(void)user_data;
	GtkListItem *li = g_object_get_data(G_OBJECT(btn), "siril-listitem");
	if (!li) return;
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gboolean new_val = siril_toggle_get_active(GTK_WIDGET(btn));
	if (new_val == row->selected) return;
	row->selected = new_val;
	int real_index = row->index - 1;
	if (real_index < 0 || real_index >= com.seq.number) return;
	/* Mirror the original on_seqlist_image_selection_toggled flow. */
	void exclude_include_single_frame(int);
	exclude_include_single_frame(real_index);
}

static void seq_setup_sel_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *btn = gtk_check_button_new();
	gtk_widget_set_halign(btn, GTK_ALIGN_CENTER);
	gtk_list_item_set_child(li, btn);
	g_signal_connect(btn, "toggled", G_CALLBACK(on_seq_selected_toggled), NULL);
}
static void seq_bind_sel_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	SirilSeqRow *row = SIRIL_SEQ_ROW(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", li);
	g_signal_handlers_block_by_func(btn, on_seq_selected_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(btn), row && row->selected);
	g_signal_handlers_unblock_by_func(btn, on_seq_selected_toggled, NULL);
}
static void seq_unbind_sel_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkCheckButton *btn = GTK_CHECK_BUTTON(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", NULL);
}

/* Column comparators */
static gint cmp_seq_index(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	const SirilSeqRow *ra = a, *rb = b;
	return ra->index - rb->index;
}
static gint cmp_seq_quality(gconstpointer a, gconstpointer b, gpointer u) { (void)u;
	const SirilSeqRow *ra = a, *rb = b;
	if (ra->fwhm < rb->fwhm) return -1;
	if (ra->fwhm > rb->fwhm) return  1;
	return 0;
}

/* Forward declaration: the original on_treeview1_cursor_changed pulled from
 * the GtkTreeSelection.  In GTK4 we wire to GtkSelectionModel
 * "selection-changed" instead. */
static void on_seq_selection_changed_model(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data);

/* set_sensitive_layers removed: GtkDropDown handles per-item sensitivity
 * directly via factories rather than the GtkCellLayout cell-data-func that
 * the GTK3 version needed. */

static void update_seqlist_dialog_combo(int layer) {
	if (!sequence_is_loaded()) return;
	sequence_list_init_statics();
	int activelayer;

	GtkDropDown *seqcombo = seqlist_dialog_combo;
	g_signal_handlers_block_by_func(GTK_DROP_DOWN(seqcombo), on_seqlist_dialog_combo_changed, NULL);
	siril_drop_down_clear_strings(seqcombo);

	if (com.seq.nb_layers == 1) {
		siril_drop_down_append_text(seqcombo, _("B&W channel"));
		activelayer = 0;
		gtk_widget_set_sensitive(GTK_WIDGET(seqcombo), FALSE);
	} else {
		siril_drop_down_append_text(seqcombo, _("Red channel"));
		siril_drop_down_append_text(seqcombo, _("Green channel"));
		siril_drop_down_append_text(seqcombo, _("Blue channel"));
		activelayer = layer >= 0 ? layer : 0;
		gtk_widget_set_sensitive(GTK_WIDGET(seqcombo), (use_photometry) ? FALSE : TRUE);
	}
	/* GTK3-era cell-layout dance for greyed-out invalid layers no longer
	 * applies to GtkDropDown.  Sensitivity of the dropdown overall is
	 * handled by gtk_widget_set_sensitive() above.  Greying-out
	 * individual entries can be reintroduced via a GtkSignalListItemFactory
	 * in a follow-up pass. */
	g_signal_handlers_unblock_by_func(GTK_DROP_DOWN(seqcombo), on_seqlist_dialog_combo_changed, NULL);

	gtk_drop_down_set_selected(GTK_DROP_DOWN(seqcombo), activelayer);

	// siril_log_debug("Resetting cached filename and index\n");
	g_free(cached_original_filename);
	cached_original_filename = NULL;
	cached_real_index = -1;
}

static void initialize_title() {
	sequence_list_init_statics();
	gchar *seq_basename;

	if (sequence_is_loaded())
		seq_basename = g_path_get_basename(com.seq.seqname);
	else
		seq_basename = g_strdup(_("No sequence loaded"));

	/* GTK4: gtk_header_bar_set_title removed */;
	/* GTK4: gtk_header_bar_set_subtitle removed */;
	gtk_widget_set_sensitive(seqlist_buttonbar, sequence_is_loaded());
	gtk_widget_set_sensitive(GTK_WIDGET(seqlist_dialog_combo), sequence_is_loaded());
	gtk_widget_set_sensitive(seqlist_refframe2, sequence_is_loaded());

	g_free(seq_basename);
}

/* GTK4: GtkColumnView has no equivalent of gtk_tree_view_set_search_entry().
 * Re-implement the same behaviour by hand: typing a number locates the row
 * with that 1-based index, selects it, and scrolls to it. */
static void on_seqlist_search_changed(GtkSearchEntry *entry, gpointer user_data) {
	(void)user_data;
	if (!seq_store || !seq_selection) return;
	const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
	if (!text || !*text) return;

	gchar *end = NULL;
	long want = strtol(text, &end, 10);
	if (end == text || want <= 0) return;

	guint n = g_list_model_get_n_items(G_LIST_MODEL(seq_store));
	for (guint i = 0; i < n; i++) {
		SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), i));
		gboolean match = (row->index == (gint)want);
		g_object_unref(row);
		if (match) {
			gtk_selection_model_select_item(GTK_SELECTION_MODEL(seq_selection), i, TRUE);
			if (seq_columnview)
				gtk_column_view_scroll_to(seq_columnview, i, NULL,
				                          GTK_LIST_SCROLL_FOCUS | GTK_LIST_SCROLL_SELECT,
				                          NULL);
			break;
		}
	}
}

static void initialize_search_entry() {
	sequence_list_init_statics();
	gtk_editable_set_text(seqlist_search_entry, "");
	if (!g_object_get_data(G_OBJECT(seqlist_search_entry), "siril-search-wired")) {
		g_signal_connect(seqlist_search_entry, "search-changed",
		                 G_CALLBACK(on_seqlist_search_changed), NULL);
		g_object_set_data(G_OBJECT(seqlist_search_entry), "siril-search-wired",
		                  GINT_TO_POINTER(1));
	}
}

/* Phase 11: build the GtkColumnView (idempotent — called once). */
static void seq_build_columnview(void) {
	if (seq_columnview) return;
	if (!seqlist_scrolled) return;

	ensure_seq_css();
	seq_store = g_list_store_new(SIRIL_TYPE_SEQ_ROW);
	seq_selection = gtk_multi_selection_new(G_LIST_MODEL(g_object_ref(seq_store)));
	g_signal_connect(seq_selection, "selection-changed",
			G_CALLBACK(on_seq_selection_changed_model), NULL);
	seq_columnview = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(seq_selection))));

	GtkColumnViewColumn *c;
	#define ADD_LBL_COLUMN(title, bind_cb, with_sort) do { \
		GtkSignalListItemFactory *_f = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new()); \
		g_signal_connect(_f, "setup", G_CALLBACK(seq_setup_label_cb), NULL); \
		g_signal_connect(_f, "bind",  G_CALLBACK(bind_cb),            NULL); \
		c = gtk_column_view_column_new(title, GTK_LIST_ITEM_FACTORY(_f)); \
		gtk_column_view_column_set_resizable(c, TRUE); \
	} while (0)

	ADD_LBL_COLUMN("#", seq_bind_index_cb, TRUE);
	gtk_column_view_column_set_sorter(c,
			GTK_SORTER(gtk_custom_sorter_new(cmp_seq_index, NULL, NULL)));
	gtk_column_view_append_column(seq_columnview, c); g_object_unref(c);

	ADD_LBL_COLUMN(N_("File"), seq_bind_file_cb, FALSE);
	gtk_column_view_column_set_expand(c, TRUE);
	gtk_column_view_append_column(seq_columnview, c); g_object_unref(c);

	ADD_LBL_COLUMN(N_("X"), seq_bind_shiftx_cb, FALSE);
	gtk_column_view_append_column(seq_columnview, c); g_object_unref(c);

	ADD_LBL_COLUMN(N_("Y"), seq_bind_shifty_cb, FALSE);
	gtk_column_view_append_column(seq_columnview, c); g_object_unref(c);

	/* Sel toggle column */
	GtkSignalListItemFactory *fsel = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fsel, "setup",  G_CALLBACK(seq_setup_sel_cb),  NULL);
	g_signal_connect(fsel, "bind",   G_CALLBACK(seq_bind_sel_cb),   NULL);
	g_signal_connect(fsel, "unbind", G_CALLBACK(seq_unbind_sel_cb), NULL);
	c = gtk_column_view_column_new(N_("Sel"), GTK_LIST_ITEM_FACTORY(fsel));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_append_column(seq_columnview, c); g_object_unref(c);

	ADD_LBL_COLUMN(N_("Quality"), seq_bind_quality_cb, TRUE);
	gtk_column_view_column_set_sorter(c,
			GTK_SORTER(gtk_custom_sorter_new(cmp_seq_quality, NULL, NULL)));
	gtk_column_view_append_column(seq_columnview, c);
	seq_quality_col = c;  /* keep ref for dynamic title set */
	g_object_ref(seq_quality_col);
	g_object_unref(c);

	#undef ADD_LBL_COLUMN

	gtk_scrolled_window_set_child(seqlist_scrolled, GTK_WIDGET(seq_columnview));
}

/* Forward stub for back-compat; the new code does not need it. */
static void get_list_store() {
	sequence_list_init_statics();
}

/* Add an image to the list. If seq is NULL, the list is cleared. */
/* Build a SirilSeqRow for `seq[index]` *without* inserting it into the
 * GListStore.  The bulk-fill path batches every row into a single
 * g_list_store_splice — populating a 44 k-frame sequence one append at
 * a time emits 44 k items-changed signals, each of which the
 * GtkColumnView reacts to, and the result is multi-second main-loop
 * stalls.  One splice = one signal. */
static SirilSeqRow *build_seq_row(sequence *seq, int index, int layer) {
	char imname[256];
	char *basename;
	int shiftx = -1, shifty = -1;
	double fwhm = -DBL_MAX;
	int color;
	double bin;

	if (!use_photometry) { // reporting registration data
		if (seq->regparam && seq->regparam[layer]) {
			double dx, dy;
			translation_from_H(seq->regparam[layer][index].H, &dx, &dy);
			shiftx = round_to_int(dx);
			shifty = round_to_int(dy);
			switch (selected_source) {
				case r_FWHM:
					if (is_arcsec) {
						bin = com.pref.binning_update ? (double) gfit->keywords.binning_x : 1.0;
						convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][index].fwhm, bin, (double) gfit->keywords.pixel_size_x, gfit->keywords.focal_length, &fwhm);
					} else {
						fwhm = seq->regparam[layer][index].fwhm;
					}
					break;
				case r_WFWHM:
					if (is_arcsec) {
						bin = com.pref.binning_update ? (double) gfit->keywords.binning_x : 1.0;
						convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][index].weighted_fwhm, bin, (double) gfit->keywords.pixel_size_x, gfit->keywords.focal_length, &fwhm);
					} else {
						fwhm = seq->regparam[layer][index].weighted_fwhm;
					}
					break;
				case r_ROUNDNESS:
					fwhm = seq->regparam[layer][index].roundness;
					break;
				case r_QUALITY:
					fwhm = seq->regparam[layer][index].quality;
					break;
				case r_BACKGROUND:
					fwhm = seq->regparam[layer][index].background_lvl;
					break;
				case r_NBSTARS:
					fwhm = seq->regparam[layer][index].number_of_stars;
					break;
				default:
					break;
			}
		}
	} else { //reporting photometry data for the reference star
		psf_star **psfs = seq->photometry[0];
		if (psfs && psfs[index]) {
			shiftx = roundf_to_int(psfs[index]->xpos);
			shifty = roundf_to_int(psfs[index]->ypos);
			switch (selected_source) {
				case ROUNDNESS:
					fwhm = psfs[index]->fwhmy / psfs[index]->fwhmx;
					break;
				case FWHM:
					if (is_arcsec) {
						fwhm_to_arcsec_if_needed(gfit, psfs[index]);
						fwhm = psfs[index]->fwhmx_arcsec < 0 ? psfs[index]->fwhmx : psfs[index]->fwhmx_arcsec;
					} else {
						fwhm = psfs[index]->fwhmx;
					}
					break;
				case AMPLITUDE:
					fwhm = psfs[index]->A;
					break;
				case MAGNITUDE:
					fwhm = psfs[index]->mag;
					break;
				case BACKGROUND:
					fwhm = psfs[index]->B;
					break;
				case X_POSITION:
					fwhm = psfs[index]->xpos;
					break;
				case Y_POSITION:
					fwhm = psfs[index]->ypos;
					break;
				case SNR:
					fwhm = psfs[index]->SNR;
					break;
				default:
					break;
			}
		}
	}

	color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
	basename = g_path_get_basename(seq_get_image_filename(seq, index, imname));
	SirilSeqRow *row = g_object_new(SIRIL_TYPE_SEQ_ROW, NULL);
	row->imname   = basename;     /* takes ownership */
	row->shiftx   = shiftx;
	row->shifty   = shifty;
	row->selected = seq->imgparam[index].incl;
	row->fwhm     = fwhm;
	row->weight   = (index == seq->current) ? 800 : 400;
	row->ref_bg   = g_strdup(index == seq->reference_image ? ref_bg_colour[color] : bg_colour[color]);
	row->index    = index + 1;
	return row;
}

/* Phase 11: change the included flag of the row whose 1-based index is
 * encoded in the path string (a stringified row number).  Walks the
 * GListStore to find the matching row. */
static void sequence_list_change_selection(gchar *path, gboolean new_value) {
	get_list_store();
	if (!seq_store) return;
	guint pos = (guint)g_ascii_strtoull(path, NULL, 10);
	if (pos >= g_list_model_get_n_items(G_LIST_MODEL(seq_store))) return;
	SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), pos));
	if (row) {
		row->selected = new_value;
		g_object_unref(row);
		g_list_model_items_changed(G_LIST_MODEL(seq_store), pos, 1, 1);
	}
}

/* Phase 11: read the real (0-based) image index from a row position in
 * the GListStore. */
static int get_image_index_from_position(guint pos) {
	get_list_store();
	if (!seq_store) return -1;
	if (pos >= g_list_model_get_n_items(G_LIST_MODEL(seq_store))) return -1;
	SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), pos));
	if (!row) return -1;
	int idx = row->index - 1;
	g_object_unref(row);
	return idx;
}

static void unselect_select_frame_from_list(gpointer unused) {
	(void)unused;
	get_list_store();
	if (!seq_selection || !seq_store) return;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(seq_selection));
	gboolean isfirst = TRUE;
	gboolean initvalue = TRUE;
	guint64 n = gtk_bitset_get_size(bs);

	/* iterate from end to start so removals (none here) keep positions stable */
	for (guint64 i = 0; i < n; i++) {
		guint pos = gtk_bitset_get_nth(bs, (guint)i);
		SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), pos));
		if (!row) continue;
		gint real_index = row->index - 1;
		g_object_unref(row);
		if (real_index < 0 || real_index >= com.seq.number) continue;
		if (isfirst) {
			initvalue = com.seq.imgparam[real_index].incl;
			isfirst = FALSE;
		}
		toggle_image_selection((int)pos, real_index, initvalue);
	}
	gtk_bitset_unref(bs);
	if (com.seq.reference_image >= 0 && com.seq.reference_image < com.seq.number
			&& !com.seq.imgparam[com.seq.reference_image].incl) {
		com.seq.reference_image = -1; // reinit to -1 in order to find new reference
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	}
	update_reg_interface(FALSE);
	update_stack_interface(FALSE);
	adjust_sellabel();
	sequence_list_change_reference();
	writeseqfile(&com.seq);
	if (!com.script) {
		redraw(REDRAW_OVERLAY);
		drawPlot();
	}
}

static void display_status() {
	sequence_list_init_statics();
	GtkLabel *statusbar = seqlist_statusbar;
	gchar *text;
	if (sequence_is_loaded()) {
		text = g_strdup_printf("%d/%d", com.seq.current + 1, com.seq.number);
		gtk_label_set_text(statusbar, text);
		g_free(text);
	} else {
		gtk_label_set_text(statusbar, "");
	}
}

/* method handling all include or all exclude from a sequence */
static void sequence_setselect_all(gboolean include_all) {
	int i;

	if (!com.seq.imgparam)
		return;
	for (i = 0; i < com.seq.number; ++i) {
		/* TODO: we include all and exclude all without checking if frame is already
		 * included of already excluded.
		 * Indeed, in the case of a list sorted by quality, index in lisrt and
		 * real index are different.
		 * We should find a way to compute the real index without using UI
		 * like we do with get_real_index_from_index_in_list()
		 *
		 */
//		if (com.seq.imgparam[i].incl != include_all) {
			com.seq.imgparam[i].incl = include_all;
			sequence_list_change_selection_index(i, i);
//		}
	}
	if (include_all) {
		com.seq.selnum = com.seq.number;
		siril_log_message(_("Selected all images from sequence\n"));
	} else {
		com.seq.selnum = 0;
		com.seq.reference_image = -1;
		siril_log_message(_("Unselected all images from sequence\n"));
		sequence_list_change_reference();
		adjust_refimage(com.seq.current);
	}
	update_reg_interface(FALSE);
	update_stack_interface(TRUE);
	writeseqfile(&com.seq);
	redraw(REDRAW_OVERLAY);
	drawPlot();
	adjust_sellabel();
}

/**** Callbacks *****/

/* Phase 11: superseded by on_seq_selection_changed_model wired to the
 * GtkSelectionModel "selection-changed" signal.  The .ui binding for
 * cursor-changed no longer exists; stub kept for any stale .ui cache. */
void on_treeview1_cursor_changed(GtkWidget *tree_view, gpointer user_data) {
	(void)tree_view; (void)user_data;
}

/* Phase 11: GTK4 replacement for the GTK3 cursor-changed handler.  When
 * exactly one row is selected, load that image (mirrors the original
 * single-click-load semantics). */
static void on_seq_selection_changed_model(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data) {
	(void)position; (void)n_items; (void)user_data;
	GtkBitset *bs = gtk_selection_model_get_selection(sm);
	if (gtk_bitset_get_size(bs) == 1) {
		guint pos = gtk_bitset_get_nth(bs, 0);
		int idx = get_image_index_from_position(pos);
		if (idx >= 0 && idx != com.seq.current) {
			fprintf(stdout, "loading image %d\n", idx);
			if (seq_load_image(&com.seq, idx, TRUE))
				seq_load_image(&com.seq, com.seq.reference_image, TRUE);
		}
	}
	gtk_bitset_unref(bs);
	display_status();
	update_reg_interface(TRUE);
}

void on_seqlist_dialog_combo_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *widget = GTK_DROP_DOWN(obj);
	(void)pspec;
	int active = gtk_drop_down_get_selected(GTK_DROP_DOWN(widget));
	if (active >= 0) {
		fill_sequence_list(&com.seq, active, FALSE);
		g_signal_handlers_block_by_func(GTK_DROP_DOWN(widget), on_seqlist_dialog_combo_changed, NULL);
		drawPlot();
		g_signal_handlers_unblock_by_func(GTK_DROP_DOWN(widget), on_seqlist_dialog_combo_changed, NULL);
	}
}

void on_column_clicked(GtkDropDown *widget, gpointer user_data) {
	int active = com.seq.current;
	sequence_list_select_row_from_index(active, FALSE);
}

void update_seqlist(int layer) {
	initialize_title();
	update_seqlist_dialog_combo(layer);
	initialize_search_entry();
	display_status();
	if (sequence_is_loaded()) sequence_list_select_row_from_index(com.seq.current, FALSE);
}

/* called on sequence loading (set_seq), on layer tab change and on registration data update.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_sequence_list(sequence *seq, int layer, gboolean as_idle) {
	struct _seq_list *args;
	if (seq == NULL || layer >= seq->nb_layers) return;

	args = calloc(1, sizeof(struct _seq_list));
	args->seq = seq;
	args->layer = layer;
	sequence_list_init_statics();
	args->tview = NULL;  /* unused under GTK4 */

	if (as_idle)
		g_idle_add(fill_sequence_list_idle, args);
	else fill_sequence_list_idle(args);
}

static gboolean fill_sequence_list_idle(gpointer p) {
	int i;
	struct _seq_list *args = (struct _seq_list *)p;
	if (!args)
		return FALSE;
	if (!args->seq) {
		g_free(args);
		return FALSE;
	}
	selected_source = gtk_drop_down_get_selected(GTK_DROP_DOWN(combo));
	is_arcsec = siril_toggle_get_active(GTK_WIDGET(arcsec));
	use_photometry = (gboolean)gtk_drop_down_get_selected(GTK_DROP_DOWN(sourceCombo));
	qualfmt = (args->seq && ((use_photometry && (selected_source == BACKGROUND)) || (!use_photometry && (selected_source == r_BACKGROUND))) && (get_data_type(args->seq->bitpix) == DATA_FLOAT)) ? ("%.5f") : ("%.3f");

	get_list_store();
	if (seq_store) g_list_store_remove_all(seq_store);

	const gchar *qual_title = N_("FWHM");
	if (!use_photometry) { // reporting registration data
		if (args->seq->regparam && args->seq->regparam[args->layer]) {
			switch (selected_source) {
				case r_FWHM:       qual_title = N_("FWHM");       break;
				case r_WFWHM:      qual_title = N_("wFWHM");      break;
				case r_ROUNDNESS:  qual_title = N_("Roundness");  break;
				case r_QUALITY:    qual_title = N_("Quality");    break;
				case r_BACKGROUND: qual_title = N_("Background"); break;
				case r_NBSTARS:    qual_title = N_("#Stars");     break;
				default: break;
			}
		}
	} else { //reporting photometry data for the reference star
		psf_star **psfs = args->seq->photometry[0];
		if (psfs && psfs[0]) {
			switch (selected_source) {
				case ROUNDNESS:  qual_title = N_("Roundness");  break;
				case FWHM:       qual_title = N_("FWHM");       break;
				case AMPLITUDE:  qual_title = N_("Amplitude");  break;
				case MAGNITUDE:  qual_title = N_("Magnitude");  break;
				case BACKGROUND: qual_title = N_("Background"); break;
				case X_POSITION: qual_title = N_("X Position"); break;
				case Y_POSITION: qual_title = N_("Y Position"); break;
				case SNR:        qual_title = N_("SNR");        break;
				default: break;
			}
		}
	}
	if (seq_quality_col)
		gtk_column_view_column_set_title(seq_quality_col, qual_title);

	/* Block selection-changed during bulk fill so the cursor-changed
	 * cascade does not run for every row added. */
	if (seq_selection)
		g_signal_handlers_block_by_func(seq_selection, on_seq_selection_changed_model, NULL);

	/* Build every row up-front, then push them into the GListStore in a
	 * single splice call.  For multi-thousand-frame SER / FITS sequences
	 * this turns a multi-second stall (one items-changed signal per row,
	 * each provoking work in the column view) into a single update. */
	if (args->seq->number > 0) {
		guint n = (guint) args->seq->number;
		gpointer *items = g_new0(gpointer, n);
		guint built = 0;
		for (i = 0; i < args->seq->number; i++) {
			SirilSeqRow *r = build_seq_row(args->seq, i, args->layer);
			if (r) items[built++] = r;
		}
		g_list_store_splice(seq_store, 0, 0, items, built);
		for (guint k = 0; k < built; k++)
			g_object_unref(items[k]);
		g_free(items);
	}

	if (seq_selection)
		g_signal_handlers_unblock_by_func(seq_selection, on_seq_selection_changed_model, NULL);

	//select and scroll to image already loaded as gfit
	sequence_list_select_row_from_index(args->seq->current, FALSE);

	free(args);
	return FALSE;
}

void exclude_include_single_frame(int index) {
	siril_log_message(_("%s image %d in sequence %s\n"),
			com.seq.imgparam[index].incl ? _("Excluding") : _("Including"),
			index + 1, com.seq.seqname);

	com.seq.imgparam[index].incl = !com.seq.imgparam[index].incl;
	if (com.seq.imgparam[index].incl)
		com.seq.selnum++;
	else {
		com.seq.selnum--;
	}
	if (com.seq.selnum == 1 || com.seq.reference_image == index) {
		if (com.seq.reference_image == index) {
			com.seq.reference_image = -1;
			com.seq.reference_image = sequence_find_refimage(&com.seq);
		} else {
			com.seq.reference_image = index;
		}
		adjust_refimage(index);
		sequence_list_change_reference();
	}
	update_reg_interface(FALSE);
	update_stack_interface(FALSE);
	redraw(REDRAW_OVERLAY);
	drawPlot();
	adjust_sellabel();
	writeseqfile(&com.seq);
}

void select_unselect_frames_from_list(gboolean *selected, gboolean keep) {
	gboolean need_ref_update = FALSE;

	for (int i = 0; i < com.seq.number; i++) { // use real index
		com.seq.imgparam[i].incl = (keep) ? selected[i] : !selected[i] && com.seq.imgparam[i].incl;
		if (!com.seq.imgparam[i].incl && com.seq.reference_image == i) { // reference image is not included anymore
			need_ref_update = TRUE;
		}
	}
	if (need_ref_update) {
		com.seq.reference_image = -1;
		com.seq.reference_image = sequence_find_refimage(&com.seq);
		adjust_refimage(com.seq.reference_image);
		sequence_list_change_reference();
	}

	redraw(REDRAW_OVERLAY);
	fix_selnum(&com.seq, FALSE);
	update_reg_interface(FALSE);
	update_stack_interface(FALSE);
	adjust_sellabel();
	writeseqfile(&com.seq);
	siril_log_message(_("Selection update finished, %d images are selected in the sequence\n"), com.seq.selnum);
}

/* Phase 11: superseded by on_seq_selected_toggled (inline GtkCheckButton).
 * Stub kept for any stale .ui binding. */
void on_seqlist_image_selection_toggled(GtkWidget *cell_renderer,
		gchar *char_path, gpointer user_data) {
	(void)cell_renderer; (void)char_path; (void)user_data;
}

void toggle_image_selection(int index_in_list, int real_index, gboolean initvalue) {
	gchar *msg = NULL;
	gboolean before_change = com.seq.imgparam[real_index].incl;
	if (initvalue) {
		com.seq.imgparam[real_index].incl = FALSE;
		if (before_change != com.seq.imgparam[real_index].incl) { // decrement only on value change
			--com.seq.selnum;
			msg = g_strdup_printf(_("Image %d has been unselected from sequence\n"), real_index + 1);
			if (com.seq.reference_image == real_index) {
				com.seq.reference_image = -1;  // invalidate to trigger new reference search if ref frame is deselected
				GtkToggleButton *refframebutton = GTK_TOGGLE_BUTTON(seqlist_refframe2); // invalidate toggle button, only effective if the frame is the first one in the selection
				g_signal_handlers_block_by_func(refframebutton, on_ref_frame_toggled, NULL);
				siril_toggle_set_active(GTK_WIDGET(refframebutton), FALSE);
				g_signal_handlers_unblock_by_func(refframebutton, on_ref_frame_toggled, NULL);
			}
		}
	} else {
		com.seq.imgparam[real_index].incl = TRUE;
		if (before_change != com.seq.imgparam[real_index].incl) { // increment only on value change
			++com.seq.selnum;
			msg = g_strdup_printf(_("Image %d has been selected from sequence\n"), real_index + 1);
		}
	}
	if (msg){
		siril_log_message(msg);
		g_free(msg);
	}
	com.seq.reference_image = sequence_find_refimage(&com.seq);
	if (!com.script) {
		sequence_list_change_selection_index(index_in_list, real_index);
	}
}

void on_selected_frames_select(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	unselect_select_frame_from_list(NULL);
}

void on_seqexcludeall_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean exclude_all = siril_confirm_dialog(_("Exclude all images?"),
			_("This erases previous image selection and there's no possible undo."),
			_("Exclude All"));
	if (exclude_all)
		sequence_setselect_all(FALSE);
}

void on_seqselectall_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean select_all = siril_confirm_dialog(_("Include all images?"),
			_("This erases previous image selection and there's no possible undo."),
			_("Include All"));
	if (select_all)
		sequence_setselect_all(TRUE);
}

void on_ref_frame_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	if (!sequence_is_loaded())
		return;

	(void)user_data;
	free_reference_image();
	get_list_store();

	/* Phase 11: this fired when the user clicked the "ref frame" toggle
	 * button.  Behaviour is unchanged from GTK3: toggle ref-status of the
	 * currently-loaded image (com.seq.current).  The original code walked
	 * the tree-selection just to check that there was at least one row
	 * selected; under GTK4 we early-out on an empty model instead. */
	if (!seq_store || g_list_model_get_n_items(G_LIST_MODEL(seq_store)) == 0)
		return;

	if (!siril_toggle_get_active(GTK_WIDGET(togglebutton))) {
		if (com.seq.reference_image == com.seq.current) {
			com.seq.reference_image = -1;
			com.seq.reference_image = sequence_find_refimage(&com.seq);
		}
	} else {
		com.seq.reference_image = com.seq.current;
		test_and_allocate_reference_image(-1);
		if (!com.seq.imgparam[com.seq.current].incl) {
			toggle_image_selection(com.seq.current, com.seq.current, FALSE);
		}
	}

	redraw(REDRAW_OVERLAY);
	sequence_list_change_reference();
	update_stack_interface(FALSE);
	adjust_sellabel();
	writeseqfile(&com.seq);
	drawPlot();
}

void on_draw_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (!sequence_is_loaded())
		return;
	redraw(REDRAW_OVERLAY);
}

/****************** modification of the list store (tree model) ******************/

void sequence_list_change_selection_index(int index_in_list, int real_index) {
	gchar *path = g_strdup_printf("%d", index_in_list);
	sequence_list_change_selection(path, com.seq.imgparam[real_index].incl);
	g_free(path);
}

void sequence_list_change_current() {
	get_list_store();
	if (!seq_store) { drawPlot(); return; }
	guint n = g_list_model_get_n_items(G_LIST_MODEL(seq_store));
	for (guint i = 0; i < n; i++) {
		SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), i));
		if (!row) continue;
		gint real_index = row->index - 1;
		gint new_weight = (real_index == com.seq.current) ? 800 : 400;
		if (row->weight != new_weight) {
			row->weight = new_weight;
			g_object_unref(row);
			g_list_model_items_changed(G_LIST_MODEL(seq_store), i, 1, 1);
			continue;
		}
		g_object_unref(row);
	}
	drawPlot();
}

void sequence_list_change_reference() {
	int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;
	get_list_store();
	if (!seq_store) return;
	guint n = g_list_model_get_n_items(G_LIST_MODEL(seq_store));
	for (guint i = 0; i < n; i++) {
		SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), i));
		if (!row) continue;
		gint real_index = row->index - 1;
		const gchar *want = (real_index == com.seq.reference_image) ? ref_bg_colour[color] : bg_colour[color];
		if (g_strcmp0(row->ref_bg, want) != 0) {
			g_free(row->ref_bg);
			row->ref_bg = g_strdup(want);
			g_object_unref(row);
			g_list_model_items_changed(G_LIST_MODEL(seq_store), i, 1, 1);
			continue;
		}
		g_object_unref(row);
	}
}

void clear_sequence_list() {
	sequence_list_init_statics();
	if (!seq_store) return;
	if (seq_selection)
		g_signal_handlers_block_by_func(seq_selection, on_seq_selection_changed_model, NULL);
	g_list_store_remove_all(seq_store);
	if (seq_selection)
		g_signal_handlers_unblock_by_func(seq_selection, on_seq_selection_changed_model, NULL);
}

void adjust_refimage(int n) {
	sequence_list_init_statics();
	g_signal_handlers_block_by_func(seqlist_refframe2, on_ref_frame_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(seqlist_refframe2), com.seq.reference_image == n);
	g_signal_handlers_unblock_by_func(seqlist_refframe2, on_ref_frame_toggled, NULL);
}

void sequence_list_select_row_from_index(int index, gboolean do_load_image) {
	sequence_list_init_statics();
	if (!seq_store || !seq_selection) return;
	guint n = g_list_model_get_n_items(G_LIST_MODEL(seq_store));
	guint found_pos = (guint)-1;
	for (guint i = 0; i < n; i++) {
		SirilSeqRow *row = SIRIL_SEQ_ROW(g_list_model_get_item(G_LIST_MODEL(seq_store), i));
		if (!row) continue;
		if (row->index - 1 == index) { found_pos = i; g_object_unref(row); break; }
		g_object_unref(row);
	}
	if (found_pos == (guint)-1) goto try_load;

	gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(seq_selection));
	gtk_selection_model_select_item(GTK_SELECTION_MODEL(seq_selection), found_pos, TRUE);
	if (seq_columnview)
		gtk_column_view_scroll_to(seq_columnview, found_pos, NULL, GTK_LIST_SCROLL_NONE, NULL);

try_load: ;

	if (do_load_image) {
		if (seq_load_image(&com.seq, index, TRUE)) // if loading fails, we fall back reloading the reference image
			seq_load_image(&com.seq, com.seq.reference_image, TRUE);
		update_reg_interface(FALSE);
		redraw(REDRAW_OVERLAY);
	}
}

void update_icons_sequence_list(gboolean is_dark) {
	GtkWidget *w;
	if (is_dark) {
		w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/frame_dark.svg");
	} else {
		w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/frame.svg");
	}
	sequence_list_init_statics();
	gtk_button_set_child(GTK_BUTTON(seqlist_drawframe_check), w);
	gtk_widget_set_visible(w, TRUE);
}

/*****************************************************************************
 *    Managing the list of sequence (the combo box) from the sequence tab
 *****************************************************************************/

#ifdef _WIN32
static int ListSequences(const gchar *sDir, const char *sequence_name_to_select,
		GtkDropDown *seqcombo, int *index_of_seq_to_load, gboolean *found) {
	WIN32_FIND_DATAW fdFile;
	HANDLE hFind = NULL;
	char sPath[2048];
	char filename[256];
	int number_of_loaded_sequences = 0;
	wchar_t *wpath;

	//Specify a file mask. *.seq = We want only seq file!
	sprintf(sPath, "%s\\*.seq", sDir);

	wpath = g_utf8_to_utf16(sPath, -1, NULL, NULL, NULL);
	if (wpath == NULL)
		return 1;

	if ((hFind = FindFirstFileW(wpath, &fdFile)) == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Path not found: [%s]\n", sDir);
		g_free(wpath);
		return 1;
	}
	g_free(wpath);

	do {
		//Is the entity a File or Folder?
		if (!(fdFile.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			gchar *cFileName = g_utf16_to_utf8(fdFile.cFileName, -1, NULL, NULL, NULL);
			if (cFileName == NULL) {
				return 1;
			}
			sequence *seq = readseqfile(cFileName);
			if (seq != NULL) {
				strncpy(filename, cFileName, 255);
				free_sequence(seq, TRUE);
				siril_drop_down_append_text(seqcombo, filename);
				if (sequence_name_to_select
						&& !strncmp(filename, sequence_name_to_select,
							strlen(filename))) {
					*index_of_seq_to_load = number_of_loaded_sequences;
					*found = TRUE;
				}
				++number_of_loaded_sequences;
			}
			g_free(cFileName);
		}
	} while (FindNextFileW(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	return number_of_loaded_sequences;
}
#endif

/** This method populates the sequence combo box with the sequences found in the CWD.
 *  If only one sequence is found, or if a sequence whose name matches the
 *  possibly NULL argument is found, it is automatically selected, which triggers
 *  its loading
 *  @param sequence_name_to_select the name of the input sequence
 *  @return 0 if success
 */
int update_sequences_list(const char *sequence_name_to_select) {
	GtkDropDown *seqcombo;
	int number_of_loaded_sequences = 0;
	int index_of_seq_to_load = -1;
	char *seqname = NULL;
	gboolean found = FALSE;

	// clear the previous list
	sequence_list_init_statics();
	seqcombo = seqlist_sequence_combobox;
	siril_drop_down_clear_strings(seqcombo);

	if (sequence_name_to_select) {
		if (g_str_has_suffix(sequence_name_to_select, ".seq"))
			seqname = strdup(sequence_name_to_select);
		else {
			seqname = malloc(strlen(sequence_name_to_select) + 5);
			sprintf(seqname, "%s.seq", sequence_name_to_select);
		}
	}

#ifdef _WIN32
	number_of_loaded_sequences = ListSequences(com.wd, seqname, seqcombo, &index_of_seq_to_load, &found);
#else
	struct dirent **list = NULL;
	int i, n;


	n = scandir(com.wd, &list, 0, alphasort);
	if (n < 0)
		perror("scandir");

	for (i = 0; i < n; ++i) {
		char *suf;

		if ((suf = strstr(list[i]->d_name, ".seq")) && strlen(suf) == 4) {
			sequence *seq = readseqfile(list[i]->d_name);
			if (seq != NULL) {
				free_sequence(seq, TRUE);
				char *filename = list[i]->d_name;
				siril_drop_down_append_text(seqcombo, filename);
				if (seqname && !strcmp(filename, seqname)) {
					index_of_seq_to_load = number_of_loaded_sequences;
					found = TRUE;
				}
				++number_of_loaded_sequences;
			}
		}
	}
	for (i = 0; i < n; i++)
		free(list[i]);
	free(list);
#endif

	if (!number_of_loaded_sequences) {
		fprintf(stderr, "No valid sequence found in CWD.\n");
		if (seqname) free(seqname);
		return -1;
	} else if (!seqname || found) {
		fprintf(stdout, "Loaded %d %s\n", number_of_loaded_sequences,
				ngettext("sequence", "sequences", number_of_loaded_sequences));
	} else return -1;

	if (seqname) free(seqname);

	if (number_of_loaded_sequences > 1 && index_of_seq_to_load < 0) {
		/* GTK4 has no public gtk_drop_down_get_popover().  GtkDropDown
		 * registers a class-level "dropdown.popup" action (used by its
		 * default key bindings); invoke it directly to open the menu. */
		gtk_widget_grab_focus(GTK_WIDGET(seqcombo));
		gtk_widget_activate_action(GTK_WIDGET(seqcombo), "dropdown.popup", NULL);
	} else if (index_of_seq_to_load >= 0)
		gtk_drop_down_set_selected(GTK_DROP_DOWN(seqcombo), index_of_seq_to_load);
	else
		gtk_drop_down_set_selected(GTK_DROP_DOWN(seqcombo), 0);
	return 0;
}

/* Phase 11: per-row "original filename" tooltip on the seq list.  In GTK4,
 * row tooltips on a GtkColumnView are set on the cell label inside the
 * bind callback rather than via a centralised query-tooltip handler.
 * This stub keeps the symbol present for any stale .ui cache; the
 * cached_real_index / cached_original_filename machinery is reused here
 * so we still maintain the lookup state, but the actual tooltip is left
 * to a follow-up pass that wires it into seq_bind_file_cb. */
gboolean on_treeview1_query_tooltip(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	(void)widget; (void)x; (void)y; (void)keyboard_tip; (void)tooltip; (void)data;
	return FALSE;
}
