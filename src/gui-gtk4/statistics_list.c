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
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "io/sequence.h"
#include "io/single_image.h"

static GtkScrolledWindow *stat_scrolled_window = NULL;
static GtkCheckButton *stat_checkbutton = NULL;
static GtkCheckButton *stat_cfa_checkbutton = NULL;
static GtkLabel *stat_name_label = NULL;
static GtkLabel *stat_selec_label = NULL;

/* GTK4: model and view live entirely in C — the .ui placeholder is filled
 * in init_dialog().  Each row is a SirilStatRow. */
static GListStore    *stat_store      = NULL;
static GtkColumnView *stat_columnview = NULL;

static void statistics_list_init_statics(void) {
	if (stat_checkbutton) return;
	stat_checkbutton     = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "statCheckButton"));
	stat_cfa_checkbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "cfastatsCheckButton"));
	stat_name_label      = GTK_LABEL(gtk_builder_get_object(gui.builder, "statNameLabel"));
	stat_selec_label     = GTK_LABEL(gtk_builder_get_object(gui.builder, "statSelecLabel"));
	stat_scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolledwindow2"));
}

char *statName[] = {
		N_("mean"),
		N_("median"),
		N_("sigma"),
		N_("avgDev"),
		N_("MAD"),
		N_("sqrt(BWMV)"),
		N_("min"),
		N_("max"),
};

/* ---- Row GObject -------------------------------------------------- */

#define SIRIL_TYPE_STAT_ROW (siril_stat_row_get_type())
G_DECLARE_FINAL_TYPE(SirilStatRow, siril_stat_row, SIRIL, STAT_ROW, GObject)

struct _SirilStatRow {
	GObject parent_instance;
	gchar *name;
	gchar *r;
	gchar *g;
	gchar *b;
};

G_DEFINE_TYPE(SirilStatRow, siril_stat_row, G_TYPE_OBJECT)

static void siril_stat_row_finalize(GObject *obj) {
	SirilStatRow *self = SIRIL_STAT_ROW(obj);
	g_clear_pointer(&self->name, g_free);
	g_clear_pointer(&self->r,    g_free);
	g_clear_pointer(&self->g,    g_free);
	g_clear_pointer(&self->b,    g_free);
	G_OBJECT_CLASS(siril_stat_row_parent_class)->finalize(obj);
}
static void siril_stat_row_class_init(SirilStatRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_stat_row_finalize;
}
static void siril_stat_row_init(SirilStatRow *self) { (void)self; }

static SirilStatRow *siril_stat_row_new(const gchar *name, const gchar *r, const gchar *g, const gchar *b) {
	SirilStatRow *self = g_object_new(SIRIL_TYPE_STAT_ROW, NULL);
	self->name = g_strdup(name);
	self->r    = g_strdup(r);
	self->g    = g_strdup(g);
	self->b    = g_strdup(b);
	return self;
}

/* ---- Cell factory: plain GtkLabel pulled from row property -------- */

typedef enum { STAT_COL_NAME, STAT_COL_R, STAT_COL_G, STAT_COL_B } StatColKind;

static const gchar *stat_row_get_col(SirilStatRow *row, StatColKind k) {
	switch (k) {
		case STAT_COL_NAME: return row->name;
		case STAT_COL_R:    return row->r;
		case STAT_COL_G:    return row->g;
		case STAT_COL_B:    return row->b;
	}
	return "";
}

static void stat_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}

static void stat_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	StatColKind kind = (StatColKind)GPOINTER_TO_INT(u);
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilStatRow *row = SIRIL_STAT_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, stat_row_get_col(row, kind));
}

static GtkListItemFactory *make_stat_factory(StatColKind kind) {
	GtkListItemFactory *f = gtk_signal_list_item_factory_new();
	g_signal_connect(f, "setup", G_CALLBACK(stat_setup_cb), NULL);
	g_signal_connect(f, "bind",  G_CALLBACK(stat_bind_cb),  GINT_TO_POINTER(kind));
	return f;
}

static void ensure_stat_view(void) {
	if (stat_columnview && stat_store) return;
	statistics_list_init_statics();

	stat_store = g_list_store_new(SIRIL_TYPE_STAT_ROW);
	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(stat_store)));
	GtkColumnView *cv = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));
	gtk_widget_set_name(GTK_WIDGET(cv), "statTreeView");

	GtkColumnViewColumn *c;
	c = gtk_column_view_column_new(N_("Name"),  make_stat_factory(STAT_COL_NAME)); gtk_column_view_column_set_resizable(c, TRUE); gtk_column_view_column_set_fixed_width(c, 100); gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Red"),   make_stat_factory(STAT_COL_R));    gtk_column_view_column_set_resizable(c, TRUE); gtk_column_view_column_set_fixed_width(c, 100); gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Green"), make_stat_factory(STAT_COL_G));    gtk_column_view_column_set_resizable(c, TRUE); gtk_column_view_column_set_fixed_width(c, 100); gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Blue"),  make_stat_factory(STAT_COL_B));    gtk_column_view_column_set_resizable(c, TRUE); gtk_column_view_column_set_fixed_width(c, 100); gtk_column_view_append_column(cv, c); g_object_unref(c);

	gtk_scrolled_window_set_child(stat_scrolled_window, GTK_WIDGET(cv));
	stat_columnview = cv;
	/* No g_object_unref(sel): gtk_column_view_new() is transfer-full for
	 * the model, so cv consumed the ref returned by single_selection_new. */
}

static void display_stat(const double *value, const double *normalization, char *format, int nblayer, int i, data_type type) {
	(void)type;
	char rvalue[20], gvalue[20], bvalue[20];

	sprintf(rvalue, format, value[RLAYER] / normalization[RLAYER]);

	if (nblayer > 1) {
		sprintf(gvalue, format, value[GLAYER] / normalization[GLAYER]);
		sprintf(bvalue, format, value[BLAYER] / normalization[BLAYER]);
	} else {
		sprintf(gvalue, "--");
		sprintf(bvalue, "--");
	}

	SirilStatRow *row = siril_stat_row_new(_(statName[i]), rvalue, gvalue, bvalue);
	g_list_store_append(stat_store, row);
	g_object_unref(row);
}

static double get_value_from_stat(imstats *stat, int index) {
	if (!stat) return 0.0;
	switch (index) {
	case 0: return stat->mean;
	case 1: return stat->median;
	case 2: return stat->sigma;
	case 3: return stat->avgDev;
	case 4: return stat->mad;
	case 5: return stat->sqrtbwmv;
	case 6: return stat->min;
	case 7: return stat->max;
	}
	return 0.0;
}

static void init_dialog() {
	statistics_list_init_statics();
	ensure_stat_view();
	g_list_store_remove_all(stat_store);
	gtk_widget_set_sensitive(GTK_WIDGET(stat_cfa_checkbutton),
			gfit->rx > 0 && gfit->naxes[2] == 1 && gfit->keywords.bayer_pattern[0] != '\0');
}

static void add_chan_stats_to_list(imstats **stat, int nblayer, data_type type, gboolean normalized) {
	double normalization[] = { 1.0, 1.0, 1.0 };
	char format[6];

	for (int i = 0; i < (int)G_N_ELEMENTS(statName); i++) {
		double value[3] = { 0.0 };

		value[RLAYER] = get_value_from_stat(stat[RLAYER], i);
		if (nblayer > 1) {
			value[GLAYER] = get_value_from_stat(stat[GLAYER], i);
			value[BLAYER] = get_value_from_stat(stat[BLAYER], i);
		}

		if (normalized) {
			normalization[RLAYER] = stat[RLAYER] ? stat[RLAYER]->normValue : 0.0;
			if (nblayer > 1) {
				normalization[GLAYER] = stat[GLAYER] ? stat[GLAYER]->normValue : 0.0;
				normalization[BLAYER] = stat[BLAYER] ? stat[BLAYER]->normValue : 0.0;
			}
			if (value[RLAYER] < 1E-5 && value[RLAYER] > 0.0) {
				sprintf(format, "%%.5e");
			} else {
				sprintf(format, "%%.7lf");
			}
		} else {
			if (type == DATA_FLOAT) {
				/* by default it is shown in ushort mode */
				normalization[RLAYER] = 1.0 / USHRT_MAX_DOUBLE;
				if (nblayer > 1) {
					normalization[GLAYER] = 1.0 / USHRT_MAX_DOUBLE;
					normalization[BLAYER] = 1.0 / USHRT_MAX_DOUBLE;
				}
			}
			sprintf(format, "%%.1lf");
		}
		display_stat(value, normalization, format, nblayer, i, type);
	}
}

void on_statButtonClose_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	siril_close_dialog("StatWindow");
}

void computeStat() {
	gboolean normalized, use_cfa;
	int channel;
	gchar *name, *selection;

	statistics_list_init_statics();
	normalized = siril_toggle_get_active(GTK_WIDGET(stat_checkbutton));
	use_cfa = siril_toggle_get_active(GTK_WIDGET(stat_cfa_checkbutton));

	if (single_image_is_loaded())
		name = g_strdup_printf("%s", com.uniq->filename);
	else if (sequence_is_loaded())
		name = g_strdup_printf(_("Image %d/%d from the sequence %s"),
				com.seq.current, com.seq.number + 1, com.seq.seqname);
	else
		name = g_strdup_printf(_("unknown image"));

	gtk_label_set_text(stat_name_label, name);
	g_free(name);

	if (com.selection.h && com.selection.w) {
		selection = g_strdup_printf(_("Size of selection in pixel: (%d,%d)"),
				com.selection.w, com.selection.h);
	} else {
		selection = g_strdup_printf(_("No selection"));
	}

	gtk_label_set_text(stat_selec_label, selection);
	g_free(selection);

	init_dialog();

	int nb_channels = (int)gfit->naxes[2];
	if (use_cfa) {
		if (nb_channels == 1 && gfit->keywords.bayer_pattern[0] != '\0')
			if ((com.selection.h || com.selection.w) && (com.selection.w < 2 || com.selection.h < 2))
				use_cfa = FALSE;
			else
				nb_channels = 3;
		else {
			use_cfa = FALSE;
		}
	}
	imstats *stat[3] = { NULL, NULL, NULL };
	for (channel = 0; channel < nb_channels; channel++) {
		int super_chan = channel;
		if (use_cfa)
			super_chan = -channel - 1;
		stat[channel] = statistics(NULL, -1, gfit, super_chan, &com.selection, STATS_MAIN, MULTI_THREADED);
		if (!stat[channel]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
		}
	}
	add_chan_stats_to_list(stat, nb_channels, gfit->type, normalized);

	for (channel = 0; channel < nb_channels; channel++) {
		free_stats(stat[channel]);
	}
}

void on_statCheckButton_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	(void)togglebutton; (void)user_data;
	set_cursor_waiting(TRUE);
	computeStat();
	set_cursor_waiting(FALSE);
}

void on_statButtonRun_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	set_cursor_waiting(TRUE);
	computeStat();
	set_cursor_waiting(FALSE);
}
