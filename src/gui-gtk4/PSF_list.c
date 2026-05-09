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
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/PSF_list.h"
#include "gui-gtk4/progress_and_log.h"
#include "algos/siril_wcs.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"

static GtkLabel *psf_statusbar = NULL;
static GtkScrolledWindow *psf_scrolled = NULL;
static GtkWindow *psf_stars_list_window = NULL;
static GtkCheckButton *psf_toggle_star_centered = NULL;

/* Phase 11: GTK4 model + view replacing liststore_stars + Stars_stored. */
static GListStore        *psf_store      = NULL;
static GtkColumnView     *psf_columnview = NULL;
static GtkMultiSelection *psf_selection  = NULL;

static void psf_build_columnview(void);

static void psf_list_init_statics(void) {
	if (psf_statusbar) return;
	psf_statusbar = GTK_LABEL(gtk_builder_get_object(gui.builder, "statusbar_PSF"));
	psf_scrolled = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolledwindow5"));
	psf_stars_list_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "stars_list_window"));
	psf_toggle_star_centered = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "toggle_star_centered"));
	psf_build_columnview();
}

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

/* ── Phase 11: per-row GObject for the PSF list ─────────────────────────
 * Mirrors the legacy 15-column GtkListStore.  Index is 1-based to match
 * the original on-screen indexing.  All double fields use NaN as the
 * "no value" sentinel; the bind callbacks substitute a blank string.
 */
#define SIRIL_TYPE_PSF_ROW (siril_psf_row_get_type())
G_DECLARE_FINAL_TYPE(SirilPSFRow, siril_psf_row, SIRIL, PSF_ROW, GObject)
struct _SirilPSFRow {
	GObject parent_instance;
	gint    channel;
	gdouble B;
	gdouble A;
	gdouble x0;
	gdouble y0;
	gdouble ra;
	gdouble dec;
	gdouble fwhmx;
	gdouble fwhmy;
	gdouble mag;
	gdouble beta;
	gdouble roundness;
	gdouble angle;
	gdouble rmse;
	gint    index;
};
G_DEFINE_TYPE(SirilPSFRow, siril_psf_row, G_TYPE_OBJECT)
static void siril_psf_row_class_init(SirilPSFRowClass *klass) { (void)klass; }
static void siril_psf_row_init(SirilPSFRow *self) { (void)self; }

/* Phase 12: GtkStatusbar replaced by GtkLabel; context-id no longer used. */

static gchar *units = "";

/* ── Phase 11: GtkColumnView cell factories ─────────────────────────────
 * Each visible column has setup/bind callbacks reading from the row.
 * setup creates a GtkLabel with right alignment for numeric columns and
 * left for the channel column.  bind formats the value.  unbind is
 * unused (no per-bind state). */

typedef enum {
	PSF_COL_CHANNEL, PSF_COL_B, PSF_COL_A,
	PSF_COL_X0, PSF_COL_Y0,
	PSF_COL_RA, PSF_COL_DEC,
	PSF_COL_FWHMX, PSF_COL_FWHMY,
	PSF_COL_MAG, PSF_COL_BETA,
	PSF_COL_ROUNDNESS, PSF_COL_ANGLE, PSF_COL_RMSE
} PSFColKind;

static void psf_setup_label_left_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void psf_setup_label_right_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 1.0);
	gtk_list_item_set_child(li, lbl);
}

static void psf_bind_double_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f;
	PSFColKind kind = (PSFColKind)GPOINTER_TO_INT(u);
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilPSFRow *row = SIRIL_PSF_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	gchar *txt = NULL;
	switch (kind) {
		case PSF_COL_CHANNEL: txt = g_strdup_printf("%d", row->channel); break;
		case PSF_COL_B:       txt = g_strdup_printf("%g", row->B);       break;
		case PSF_COL_A:       txt = g_strdup_printf("%g", row->A);       break;
		case PSF_COL_X0:      txt = g_strdup_printf("%.2f", row->x0);    break;
		case PSF_COL_Y0:      txt = g_strdup_printf("%.2f", row->y0);    break;
		case PSF_COL_RA:      txt = g_strdup_printf("%g", row->ra);      break;
		case PSF_COL_DEC:     txt = g_strdup_printf("%g", row->dec);     break;
		case PSF_COL_FWHMX:   txt = g_strdup_printf("%.2f%s", row->fwhmx, units); break;
		case PSF_COL_FWHMY:   txt = g_strdup_printf("%.2f%s", row->fwhmy, units); break;
		case PSF_COL_MAG:     txt = g_strdup_printf("%.6f", row->mag);   break;
		case PSF_COL_BETA:    txt = g_strdup_printf("%g", row->beta);    break;
		case PSF_COL_ROUNDNESS: txt = g_strdup_printf("%g", row->roundness); break;
		case PSF_COL_ANGLE:   txt = g_strdup_printf("%g", row->angle);   break;
		case PSF_COL_RMSE:    txt = g_strdup_printf("%.3e", row->rmse);  break;
	}
	gtk_label_set_text(lbl, txt ? txt : "");
	g_free(txt);
}

#define DEFINE_PSF_CMP(name, field) \
	static gint cmp_psf_##name(gconstpointer a, gconstpointer b, gpointer u) { (void)u; \
		const SirilPSFRow *ra = a, *rb = b; \
		if (ra->field < rb->field) return -1; \
		if (ra->field > rb->field) return  1; \
		return 0; \
	}
DEFINE_PSF_CMP(B, B)
DEFINE_PSF_CMP(A, A)
DEFINE_PSF_CMP(x0, x0)
DEFINE_PSF_CMP(y0, y0)
DEFINE_PSF_CMP(ra, ra)
DEFINE_PSF_CMP(dec, dec)
DEFINE_PSF_CMP(fwhmx, fwhmx)
DEFINE_PSF_CMP(fwhmy, fwhmy)
DEFINE_PSF_CMP(mag, mag)
DEFINE_PSF_CMP(beta, beta)
DEFINE_PSF_CMP(roundness, roundness)
DEFINE_PSF_CMP(angle, angle)
DEFINE_PSF_CMP(rmse, rmse)
#undef DEFINE_PSF_CMP

/* Forward declaration: defined further down (replaces the GTK3
 * GtkTreeSelection "changed" handler). */
static void on_psf_selection_changed_model(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data);
static gboolean on_psf_columnview_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data);


/* Phase 11: build the GtkColumnView (idempotent — called once). */
static void psf_build_columnview(void) {
	if (psf_columnview) return;
	if (!psf_scrolled) return;

	psf_store = g_list_store_new(SIRIL_TYPE_PSF_ROW);
	psf_selection = gtk_multi_selection_new(G_LIST_MODEL(g_object_ref(psf_store)));
	g_signal_connect(psf_selection, "selection-changed",
			G_CALLBACK(on_psf_selection_changed_model), NULL);
	psf_columnview = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(psf_selection))));

	GtkColumnViewColumn *c;

	#define ADD_NUM_COLUMN(title, kind, cmpfn) do { \
		GtkSignalListItemFactory *_f = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new()); \
		g_signal_connect(_f, "setup", G_CALLBACK(psf_setup_label_right_cb), NULL); \
		g_signal_connect(_f, "bind",  G_CALLBACK(psf_bind_double_cb),       GINT_TO_POINTER(kind)); \
		c = gtk_column_view_column_new(title, GTK_LIST_ITEM_FACTORY(_f)); \
		gtk_column_view_column_set_resizable(c, TRUE); \
		gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmpfn, NULL, NULL))); \
		gtk_column_view_append_column(psf_columnview, c); g_object_unref(c); \
	} while (0)

	/* Channel column - left aligned */
	GtkSignalListItemFactory *fch = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fch, "setup", G_CALLBACK(psf_setup_label_left_cb), NULL);
	g_signal_connect(fch, "bind",  G_CALLBACK(psf_bind_double_cb),      GINT_TO_POINTER(PSF_COL_CHANNEL));
	c = gtk_column_view_column_new("Ch.", GTK_LIST_ITEM_FACTORY(fch));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_append_column(psf_columnview, c); g_object_unref(c);

	ADD_NUM_COLUMN("B",     PSF_COL_B,     cmp_psf_B);
	ADD_NUM_COLUMN("A",     PSF_COL_A,     cmp_psf_A);
	ADD_NUM_COLUMN("x0",    PSF_COL_X0,    cmp_psf_x0);
	ADD_NUM_COLUMN("y0",    PSF_COL_Y0,    cmp_psf_y0);
	ADD_NUM_COLUMN("Ra",    PSF_COL_RA,    cmp_psf_ra);
	ADD_NUM_COLUMN("Dec",   PSF_COL_DEC,   cmp_psf_dec);
	ADD_NUM_COLUMN("FWHMx", PSF_COL_FWHMX, cmp_psf_fwhmx);
	ADD_NUM_COLUMN("FWHMy", PSF_COL_FWHMY, cmp_psf_fwhmy);
	ADD_NUM_COLUMN("Mag",   PSF_COL_MAG,   cmp_psf_mag);
	ADD_NUM_COLUMN("Beta",  PSF_COL_BETA,  cmp_psf_beta);
	ADD_NUM_COLUMN("r",     PSF_COL_ROUNDNESS, cmp_psf_roundness);
	ADD_NUM_COLUMN("Angle", PSF_COL_ANGLE, cmp_psf_angle);
	ADD_NUM_COLUMN("RMSE",  PSF_COL_RMSE,  cmp_psf_rmse);

	#undef ADD_NUM_COLUMN

	gtk_scrolled_window_set_child(psf_scrolled, GTK_WIDGET(psf_columnview));

	/* GTK4: Delete-key handling on the column view (replaces the GTK3
	 * "key-release-event" on Stars_stored). */
	GtkEventController *key_ctl = gtk_event_controller_key_new();
	g_signal_connect(key_ctl, "key-pressed",
	                 G_CALLBACK(on_psf_columnview_key_pressed), NULL);
	gtk_widget_add_controller(GTK_WIDGET(psf_columnview), key_ctl);
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
	gint result = -1;

	g_rw_lock_reader_lock(&com.stars_lock);
	while (com.stars && com.stars[i]) {
		if ((com.stars[i]->xpos == x) && (com.stars[i]->ypos == y)) {
			result = i;
			break;
		}
		i++;
	}
	g_rw_lock_reader_unlock(&com.stars_lock);
	return result;
}

static void display_status() {
	gchar *text;
	int i = 0;

	psf_list_init_statics();

	g_rw_lock_reader_lock(&com.stars_lock);
	while (com.stars && com.stars[i])
		i++;
	g_rw_lock_reader_unlock(&com.stars_lock);
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
	gtk_label_set_text(psf_statusbar, text);
	g_free(text);
}

void set_iter_of_clicked_psf(double x, double y) {
	psf_list_init_statics();
	if (!psf_store || !psf_selection) return;
	gboolean is_as;
	const double radian_conversion = ((3600.0 * 180.0) / M_PI) / 1.0E3;
	double invpixscalex = 1.0;
	double bin_X = com.pref.binning_update ? (double) gfit->keywords.binning_x : 1.0;
	g_rw_lock_reader_lock(&com.stars_lock);
	if (com.stars && com.stars[0]) {
		is_as = (strcmp(com.stars[0]->units, "px"));
	} else {
		g_rw_lock_reader_unlock(&com.stars_lock);
		return;
	}
	g_rw_lock_reader_unlock(&com.stars_lock);
	if (is_as) {
		invpixscalex = 1.0 / (radian_conversion * (double) gfit->keywords.pixel_size_x / gfit->keywords.focal_length) * bin_X;
	}
	guint n = g_list_model_get_n_items(G_LIST_MODEL(psf_store));
	for (guint i = 0; i < n; i++) {
		SirilPSFRow *row = SIRIL_PSF_ROW(g_list_model_get_item(G_LIST_MODEL(psf_store), i));
		if (!row) continue;
		gdouble xpos = row->x0, ypos = row->y0;
		gdouble fwhmx = row->fwhmx * invpixscalex;
		g_object_unref(row);
		gdouble distsq = (xpos - x) * (xpos - x) + (ypos - y) * (ypos - y);
		gdouble psflimsq = 6. * fwhmx * fwhmx;
		if (distsq < psflimsq) {
			gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(psf_selection));
			gtk_selection_model_select_item(GTK_SELECTION_MODEL(psf_selection), i, TRUE);
			if (psf_columnview)
				gtk_column_view_scroll_to(psf_columnview, i, NULL, GTK_LIST_SCROLL_NONE, NULL);
			gui.selected_star = get_index_of_selected_star(xpos, ypos);
			gtk_window_present(psf_stars_list_window);
			display_status();
			redraw(REDRAW_OVERLAY);
			return;
		}
	}
	siril_debug_print("Point clicked does not correspond to a known star\n");
}

static int compare(void const *a, void const *b) {
	guint const *pa = a;
	guint const *pb = b;

	return *pa - *pb;
}

/* Phase 11: walk the GListStore and decrement each row's index by the
 * number of removed entries with smaller index.  sel must be sorted
 * ascending. */
static void update_column_index_glist(const guint *sel, guint size) {
	if (!psf_store) return;
	guint n = g_list_model_get_n_items(G_LIST_MODEL(psf_store));
	for (guint k = 0; k < n; k++) {
		SirilPSFRow *row = SIRIL_PSF_ROW(g_list_model_get_item(G_LIST_MODEL(psf_store), k));
		if (!row) continue;
		int i;
		for (i = size - 1; i >= 0 && row->index < (gint)sel[i]; i--) ;
		gint new_idx = row->index - i - 1;
		if (new_idx != row->index) {
			row->index = new_idx;
			g_object_unref(row);
			g_list_model_items_changed(G_LIST_MODEL(psf_store), k, 1, 1);
			continue;
		}
		g_object_unref(row);
	}
}

static void remove_selected_star() {
	psf_list_init_statics();
	if (!psf_store || !psf_selection) return;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(psf_selection));
	guint64 nsel = gtk_bitset_get_size(bs);
	if (nsel == 0) { gtk_bitset_unref(bs); return; }

	guint *sel = calloc(nsel, sizeof(guint));

	/* Read the row indices BEFORE removing anything (positions stable). */
	for (guint64 i = 0; i < nsel; i++) {
		guint pos = gtk_bitset_get_nth(bs, (guint)i);
		SirilPSFRow *row = SIRIL_PSF_ROW(g_list_model_get_item(G_LIST_MODEL(psf_store), pos));
		if (row) {
			sel[i] = (guint)row->index;
			g_object_unref(row);
		}
	}

	/* Now remove from highest position downward so positions stay stable. */
	for (guint64 i = nsel; i > 0; i--) {
		guint pos = gtk_bitset_get_nth(bs, (guint)(i - 1));
		g_list_store_remove(psf_store, pos);
	}
	gtk_bitset_unref(bs);

	qsort(sel, nsel, sizeof *sel, compare);

	for (gint j = (gint)nsel - 1; j >= 0; j--) {
		remove_star(sel[j] - 1);
	}

	update_column_index_glist(sel, (guint)nsel);

	gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(psf_selection));
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
	G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
	gtk_file_filter_add_pattern(f, "*.lst");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
	G_GNUC_END_IGNORE_DEPRECATIONS
}

/* Phase 11: walk the GListStore and write a fixed-format CSV with the 14
 * visible columns in display order.  Header names match the column
 * titles built in psf_build_columnview(). */
static void psf_export_csv(const char *filename) {
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

	GDataOutputStream *data_stream = g_data_output_stream_new(G_OUTPUT_STREAM(output_stream));

	ret &= g_data_output_stream_put_string(data_stream,
			"Ch.,B,A,x0,y0,Ra,Dec,FWHMx,FWHMy,Mag,Beta,r,Angle,RMSE\n",
			NULL, NULL);

	if (psf_store) {
		guint n = g_list_model_get_n_items(G_LIST_MODEL(psf_store));
		for (guint i = 0; i < n; i++) {
			SirilPSFRow *row = SIRIL_PSF_ROW(g_list_model_get_item(G_LIST_MODEL(psf_store), i));
			if (!row) continue;
			gchar *line = g_strdup_printf("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
					row->channel, row->B, row->A, row->x0, row->y0,
					row->ra, row->dec, row->fwhmx, row->fwhmy,
					row->mag, row->beta, row->roundness, row->angle, row->rmse);
			ret &= g_data_output_stream_put_string(data_stream, line, NULL, NULL);
			g_free(line);
			g_object_unref(row);
		}
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
	psf_list_init_statics();
	GtkWindow *parent = psf_stars_list_window;
	gint res;

	widgetdialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
	siril_file_chooser_set_current_folder_path(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	/* GTK4: gtk_file_chooser_set_do_overwrite_confirmation removed */;
	gtk_file_chooser_set_current_name(dialog, "stars.lst");
	G_GNUC_END_IGNORE_DEPRECATIONS
	/* GTK4: gtk_file_chooser_set_local_only removed */;
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_file_chooser_get_filename(dialog);
		psf_export_csv(file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

int get_ra_and_dec_from_star_pos(psf_star *star, gdouble *alpha, gdouble *delta) {
	int ret = 1;
	if (has_wcs(gfit)) {
		// coordinates of the star in FITS/WCS coordinates
		double fx, fy;
		display_to_siril(star->xpos, star->ypos, &fx, &fy, gfit->ry);

		double ra, dec;
		pix2wcs(gfit, fx, fy, &ra, &dec);
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
	psf_list_init_statics();
	if (!psf_store) return;
	if (star == NULL) {
		g_list_store_remove_all(psf_store);
		return;
	}

	double fwhmx = star->fwhmx_arcsec < 0 ? star->fwhmx : star->fwhmx_arcsec;
	double fwhmy = star->fwhmy_arcsec < 0 ? star->fwhmy : star->fwhmy_arcsec;

	SirilPSFRow *row = g_object_new(SIRIL_TYPE_PSF_ROW, NULL);
	row->channel   = star->layer;
	row->B         = star->B;
	row->A         = star->A;
	row->x0        = star->xpos;
	row->y0        = star->ypos;
	row->ra        = star->ra;
	row->dec       = star->dec;
	row->fwhmx     = fwhmx;
	row->fwhmy     = fwhmy;
	row->mag       = star->mag + com.magOffset;
	row->beta      = star->beta;
	row->roundness = fwhmy / fwhmx;
	row->angle     = -star->angle;
	row->rmse      = star->rmse;
	row->index     = i + 1;
	g_list_store_append(psf_store, row);
	g_object_unref(row);

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
	psf_list_init_statics();
	if (psf_store) g_list_store_remove_all(psf_store);
	g_rw_lock_reader_lock(&com.stars_lock);
	fill_stars_list(gfit, com.stars);
	g_rw_lock_reader_unlock(&com.stars_lock);
}

/* this can be called from any thread; com.stars_lock (writer) serialises it */
/* clear_stars_list and clear_stars_list_as_idle moved to algos/PSF.c */

void clear_psf_list_display(void) {
	psf_list_init_statics();
	gui.selected_star = -1;
	if (psf_store) g_list_store_remove_all(psf_store);
	display_status();
}


struct star_update_s {
	psf_star **stars;
	gboolean update_GUI;
};

static gboolean update_stars_idle(gpointer p) {
	struct star_update_s *args = (struct star_update_s *)p;
	clear_stars_list(TRUE);
	g_rw_lock_writer_lock(&com.stars_lock);
	com.stars = args->stars;
	g_rw_lock_writer_unlock(&com.stars_lock);
	if (args->update_GUI && !com.headless) {
		g_rw_lock_reader_lock(&gfit->rwlock);
		g_rw_lock_reader_lock(&com.stars_lock);
		fill_stars_list(gfit, com.stars);
		g_rw_lock_reader_unlock(&com.stars_lock);
		g_rw_lock_reader_unlock(&gfit->rwlock);
	}
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
	g_rw_lock_reader_lock(&com.stars_lock);
	int i = 0;
	while (com.stars[i])
		i++;
	g_rw_lock_reader_unlock(&com.stars_lock);
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
		psf_star *new_star = add_star(gfit, layer, &new_index);
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

/* Phase 11: superseded by on_psf_selection_changed_model wired to the
 * GtkSelectionModel "selection-changed" signal.  Stub kept for any
 * stale .ui binding. */
void on_treeview_selection_changed(GObject *selection, gpointer user_data) {
	(void)selection; (void)user_data;
}

/* Phase 11: GTK4 replacement for the old GtkTreeSelection "changed"
 * handler.  When exactly one row is selected, set the crosshairs and
 * (optionally) centre the viewport on the star's position. */
static void on_psf_selection_changed_model(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data) {
	(void)position; (void)n_items; (void)user_data;
	psf_list_init_statics();
	if (!psf_store) return;
	GtkBitset *bs = gtk_selection_model_get_selection(sm);
	if (gtk_bitset_get_size(bs) != 1) {
		gtk_bitset_unref(bs);
		return;
	}
	guint pos = gtk_bitset_get_nth(bs, 0);
	gtk_bitset_unref(bs);
	SirilPSFRow *row = SIRIL_PSF_ROW(g_list_model_get_item(G_LIST_MODEL(psf_store), pos));
	if (!row) return;
	gdouble x0 = row->x0, y0 = row->y0;
	g_object_unref(row);

	const gchar *area[] = {"drawingarear", "drawingareag", "drawingareab", "drawingareargb"};
	GtkWidget *widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, area[gui.cvport]));

	gui.selected_star = get_index_of_selected_star(x0, y0);
	GtkCheckButton *toggle = psf_toggle_star_centered;
	if (siril_toggle_get_active(GTK_WIDGET(toggle))) {
		double z = get_zoom_val();
		gui.display_offset.x = (gtk_widget_get_width(widget) / 2 - x0 * z);
		gui.display_offset.y = (gtk_widget_get_height(widget) / 2 - y0 * z);
		adjust_vport_size_to_image();
	}
	display_status();
	redraw(REDRAW_OVERLAY);
}

/* GTK4 key handler: Delete / BackSpace removes the selected star.
 * Replaces the GTK3 "key-release-event" signal on the GtkTreeView. */
static gboolean on_psf_columnview_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data) {
	(void)controller; (void)keycode; (void)state; (void)user_data;
	if (keyval == GDK_KEY_Delete || keyval == GDK_KEY_KP_Delete
			|| keyval == GDK_KEY_BackSpace) {
		remove_selected_star();
		return TRUE;
	}
	return FALSE;
}


void on_stars_list_window_hide(GtkWidget *object, gpointer user_data) {
	gui.selected_star = -1;
}

void on_sum_button_clicked(GtkButton *button, gpointer user_data) {
	g_rw_lock_reader_lock(&com.stars_lock);
	psf_star **stars_snap = com.stars;
	g_rw_lock_reader_unlock(&com.stars_lock);
	display_PSF(stars_snap);
}

void on_add_button_clicked(GtkButton *button, gpointer user_data) {
	int layer = match_drawing_area_widget(gui.view[gui.cvport].drawarea, FALSE);
	if (layer == -1)
		layer = 1;
	int index;
	psf_star *new_star = add_star(gfit, layer, &index);
	if (index > -1)
		add_star_to_list(new_star, index);
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
	g_rw_lock_reader_lock(&com.stars_lock);
	gboolean have_stars_export = (com.stars != NULL);
	g_rw_lock_reader_unlock(&com.stars_lock);
	if (have_stars_export) {
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
