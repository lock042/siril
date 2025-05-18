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

#ifdef _WIN32
#include <windows.h>
#endif
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/plot.h"
#include "gui/registration.h"	// for update_reg_interface
#include "gui/stacking.h"	// for update_stack_interface
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

static  GtkWidget *combo = NULL, *arcsec = NULL, *sourceCombo = NULL;
static int selected_source = -1;
static gboolean is_arcsec = FALSE;
static gboolean use_photometry = FALSE;

struct _seq_list {
	GtkTreeView *tview;
	sequence *seq;
	int layer;
};

static GtkListStore *list_store = NULL;

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

enum {					//different context_id of the GtkStatusBar
	COUNT_STATE
};

void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data);

/******* Static functions **************/
static void fwhm_quality_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble quality;
	gchar buf[20];
	gtk_tree_model_get(model, iter, COLUMN_FWHM, &quality, -1);
	if (quality > -DBL_MAX)
		g_snprintf(buf, sizeof(buf), qualfmt, quality);
	else
		g_strlcpy(buf, "N/A", sizeof(buf));
	g_object_set(renderer, "text", buf, NULL);
}

static void set_sensitive_layers(GtkCellLayout *cell_layout,
		GtkCellRenderer *cell, GtkTreeModel *tree_model, GtkTreeIter *iter,
		gpointer data) {
	gboolean sensitive = TRUE;

	if (sequence_is_loaded() && !use_photometry) {
		GtkTreePath *path = gtk_tree_model_get_path(tree_model, iter);
		if (!path) return;
		const gint *index = gtk_tree_path_get_indices(path); // search by index to avoid translation problems
		if (!(com.seq.regparam))
			return;
		if (com.seq.regparam[*index] == NULL)
			sensitive = FALSE;
		gtk_tree_path_free(path);
	}
	g_object_set(cell, "sensitive", sensitive, NULL);
}

static void update_seqlist_dialog_combo(int layer) {
	if (!sequence_is_loaded()) return;
	int activelayer;

	GtkComboBoxText *seqcombo = GTK_COMBO_BOX_TEXT(lookup_widget("seqlist_dialog_combo"));
	g_signal_handlers_block_by_func(GTK_COMBO_BOX(seqcombo), on_seqlist_dialog_combo_changed, NULL);
	gtk_combo_box_text_remove_all(seqcombo);

	if (com.seq.nb_layers == 1) {
		gtk_combo_box_text_append_text(seqcombo, _("B&W channel"));
		activelayer = 0;
		gtk_widget_set_sensitive(GTK_WIDGET(seqcombo), FALSE);
	} else {
		gtk_combo_box_text_append_text(seqcombo, _("Red channel"));
		gtk_combo_box_text_append_text(seqcombo, _("Green channel"));
		gtk_combo_box_text_append_text(seqcombo, _("Blue channel"));
		activelayer = layer >= 0 ? layer : 0;
		gtk_widget_set_sensitive(GTK_WIDGET(seqcombo), (use_photometry) ? FALSE : TRUE);
	}
	//setting sensitvity of the different layers in the layer combo
	gtk_cell_layout_clear(GTK_CELL_LAYOUT(seqcombo));
	GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
	gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(seqcombo), renderer, TRUE);
	gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(seqcombo), renderer, "text", 0, NULL);
	gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(seqcombo), renderer, set_sensitive_layers, NULL, NULL);
	g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(seqcombo), on_seqlist_dialog_combo_changed, NULL);

	gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), activelayer);
}

static void initialize_title() {
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("seqlistbar"));
	gchar *seq_basename;

	if (sequence_is_loaded())
		seq_basename = g_path_get_basename(com.seq.seqname);
	else
		seq_basename = g_strdup(_("No sequence loaded"));

	gtk_header_bar_set_title(bar, _("Frame List"));
	gtk_header_bar_set_subtitle(bar, seq_basename);
	gtk_widget_set_sensitive(lookup_widget("seqlist_buttonbar"), sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("seqlist_dialog_combo"), sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("refframe2"), sequence_is_loaded());

	g_free(seq_basename);
}

static void initialize_search_entry() {
	GtkTreeView *tree_view = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "treeview1"));
	GtkEntry *entry = GTK_ENTRY(lookup_widget("seqlistsearch"));

	gtk_entry_set_text (entry, "");
	gtk_tree_view_set_search_entry(tree_view, entry);
}

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststore1"));

		GtkTreeViewColumn *col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5"));
		GtkCellRenderer *cell = GTK_CELL_RENDERER(gtk_builder_get_object(gui.builder, "cellrenderertext5"));
		gtk_tree_view_column_set_cell_data_func(col, cell, fwhm_quality_cell_data_function, NULL, NULL);
	}
}

/* Add an image to the list. If seq is NULL, the list is cleared. */
static void add_image_to_sequence_list(sequence *seq, int index, int layer) {
	GtkTreeIter iter;
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
						bin = com.pref.binning_update ? (double) gfit.keywords.binning_x : 1.0;
						convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][index].fwhm, bin, (double) gfit.keywords.pixel_size_x, gfit.keywords.focal_length, &fwhm);
					} else {
						fwhm = seq->regparam[layer][index].fwhm;
					}
					break;
				case r_WFWHM:
					if (is_arcsec) {
						bin = com.pref.binning_update ? (double) gfit.keywords.binning_x : 1.0;
						convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][index].weighted_fwhm, bin, (double) gfit.keywords.pixel_size_x, gfit.keywords.focal_length, &fwhm);
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
						fwhm_to_arcsec_if_needed(&gfit, psfs[index]);
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
	gtk_list_store_append (list_store, &iter);
	gtk_list_store_set (list_store, &iter,
			COLUMN_IMNAME, basename,
			COLUMN_SHIFTX, shiftx,
			COLUMN_SHIFTY, shifty,
			COLUMN_SELECTED, seq->imgparam[index].incl,
			COLUMN_FWHM, fwhm,
			COLUMN_CURRENT, index == seq->current ? 800 : 400,
			// weight value is 400 by default "normal":
			// http://developer.gnome.org/gtk3/stable/GtkCellRendererText.html#GtkCellRendererText--weight
			COLUMN_REFERENCE, index == seq->reference_image ?
			ref_bg_colour[color] : bg_colour[color],
			COLUMN_INDEX, (index + 1),
			-1);
	/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
	g_free(basename);
}

static void sequence_list_change_selection(gchar *path, gboolean new_value) {
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_list_store_set(list_store, &iter, COLUMN_SELECTED, new_value, -1);
}

static int get_image_index_from_path(GtkTreePath *path) {
	GValue value = G_VALUE_INIT;
	gint index;
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_tree_model_get_value(GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
	index = g_value_get_int(&value) - 1;
	g_value_unset(&value);
	return index;
}

static gint get_real_index_from_index_in_list(GtkTreeModel *model, GtkTreeIter *iter) {
	gint real_index;
	gtk_tree_model_get(model, iter, COLUMN_INDEX, &real_index, -1);
	return real_index - 1;
}

static void unselect_select_frame_from_list(GtkTreeView *tree_view) {
	GtkTreeSelection *selection;
	GtkTreeModel *model;
	GList *references, *list;
	gboolean isfirst = TRUE;
	gboolean initvalue;

	get_list_store();

	model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection(tree_view);
	references = get_row_references_of_selected_rows(selection, model);
	references = g_list_reverse(references); // for some reason, refs are returned in the reverse order of selection

	for (list = references; list; list = list->next) {
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			const gint *new_index = gtk_tree_path_get_indices(path);
			GtkTreeIter iter;
			gtk_tree_model_get_iter(GTK_TREE_MODEL(list_store), &iter, path);
			gint real_index = get_real_index_from_index_in_list(model, &iter);
			if (isfirst) { // replicate behavior of first selected row for the subsequent rows
				initvalue = com.seq.imgparam[real_index].incl;
				isfirst = FALSE;
			}
			toggle_image_selection(new_index[0], real_index, initvalue);
			gtk_tree_path_free(path);
		}
	}
	if (!com.seq.imgparam[com.seq.reference_image].incl) {
		com.seq.reference_image = -1; // reinit to -1 in order to find new reference
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	}
	g_list_free(references);
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
	GtkStatusbar *statusbar = GTK_STATUSBAR(lookup_widget("seqlist_statusbar"));
	gchar *text;
	if (sequence_is_loaded()) {
		text = g_strdup_printf("%d/%d", com.seq.current + 1, com.seq.number);
		gtk_statusbar_push(statusbar, COUNT_STATE, text);
		g_free(text);
	} else {
		gtk_statusbar_push(statusbar, COUNT_STATE, "");
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

void on_treeview1_cursor_changed(GtkTreeView *tree_view, gpointer user_data) {
	GtkTreeModel *tree_model;
	GtkTreeSelection *selection;
	GtkTreeIter iter;
	GList *list;

	tree_model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection (tree_view);

	list = gtk_tree_selection_get_selected_rows(selection, &tree_model);
	if (g_list_length(list) == 1) {
		gint idx;
		GValue value = G_VALUE_INIT;
		gtk_tree_model_get_iter(tree_model, &iter, (GtkTreePath *)list->data);

		gtk_tree_model_get_value(tree_model, &iter, COLUMN_INDEX, &value);
		idx = g_value_get_int(&value) - 1;
		if (idx != com.seq.current) {
			fprintf(stdout, "loading image %d\n", idx);
			if (seq_load_image(&com.seq, idx, TRUE)) // if loading fails, we fall back reloading the reference image
				seq_load_image(&com.seq, com.seq.reference_image, TRUE);
		}
		g_value_unset(&value);
	}
	g_list_free_full(list, (GDestroyNotify) gtk_tree_path_free);
	display_status();
	update_reg_interface(TRUE);
}

void on_seqlist_dialog_combo_changed(GtkComboBoxText *widget, gpointer user_data) {
	int active = gtk_combo_box_get_active(GTK_COMBO_BOX(widget));
	if (active >= 0) {
		fill_sequence_list(&com.seq, active, FALSE);
		g_signal_handlers_block_by_func(GTK_COMBO_BOX(widget), on_seqlist_dialog_combo_changed, NULL);
		drawPlot();
		g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(widget), on_seqlist_dialog_combo_changed, NULL);
	}
}

void on_column_clicked(GtkComboBoxText *widget, gpointer user_data) {
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
	args->tview = GTK_TREE_VIEW(lookup_widget("treeview1"));

	g_signal_handlers_block_by_func(args->tview, on_treeview1_cursor_changed, NULL);

	if (as_idle)
		gdk_threads_add_idle(fill_sequence_list_idle, args);
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
	if (combo == NULL) combo = lookup_widget("plotCombo");
	if (sourceCombo == NULL) sourceCombo = lookup_widget("plotSourceCombo");
	if (arcsec == NULL) arcsec = lookup_widget("arcsecPhotometry");
	selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	is_arcsec = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(arcsec));
	use_photometry = (gboolean)gtk_combo_box_get_active(GTK_COMBO_BOX(sourceCombo));
	qualfmt = (args->seq && ((use_photometry && (selected_source == BACKGROUND)) || (!use_photometry && (selected_source == r_BACKGROUND))) && (get_data_type(args->seq->bitpix) == DATA_FLOAT)) ? ("%.5f") : ("%.3f");

	if (list_store) gtk_list_store_clear(list_store);
	get_list_store();
	if (!use_photometry) { // reporting registration data
		if (args->seq->regparam && args->seq->regparam[args->layer]) {
			switch (selected_source) {
				case r_FWHM:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("FWHM"));
					break;
				case r_WFWHM:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("wFWHM"));
					break;
				case r_ROUNDNESS:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Roundness"));
					break;
				case r_QUALITY:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Quality"));
					break;
				case r_BACKGROUND:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Background"));
					break;
				case r_NBSTARS:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("#Stars"));
					break;
				default:
					break;
			}
		} else {
			gtk_tree_view_column_set_title (GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("FWHM"));
		}
	} else { //reporting photometry data for the reference star
		psf_star **psfs = args->seq->photometry[0];
		if (psfs && psfs[0]) {
			switch (selected_source) {
				case ROUNDNESS:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Roundness"));
					break;
				case FWHM:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("FWHM"));
					break;
				case AMPLITUDE:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Amplitude"));
					break;
				case MAGNITUDE:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Magnitude"));
					break;
				case BACKGROUND:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Background"));
					break;
				case X_POSITION:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("X Position"));
					break;
				case Y_POSITION:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("Y Position"));
					break;
				case SNR:
					gtk_tree_view_column_set_title(GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("SNR"));
					break;
				default:
					break;
			}
		} else {
			gtk_tree_view_column_set_title (GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(gui.builder, "treeviewcolumn5")), _("FWHM"));
		}
	}
	gint sort_column_id;
	GtkSortType order;
	// store sorted state of list_store, disable sorting, disconnect from the view, fill, reconnect and re-apply sort
	gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(list_store), &sort_column_id, &order);
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID, GTK_SORT_ASCENDING);
	gtk_tree_view_set_model(args->tview, NULL);
	if (args->seq->number > 0) {
		for (i = 0; i < args->seq->number; i++) {
			add_image_to_sequence_list(args->seq, i, args->layer);
		}
	}
	gtk_tree_view_set_model(args->tview, GTK_TREE_MODEL(list_store));
	gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list_store), sort_column_id, order);

	//select and scroll to image already loaded as gfit
	sequence_list_select_row_from_index(args->seq->current, FALSE);
	g_signal_handlers_unblock_by_func(args->tview, on_treeview1_cursor_changed, NULL);

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

void on_seqlist_image_selection_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *char_path, gpointer user_data) {
	GtkTreePath *path = gtk_tree_path_new_from_string(char_path);
	gint index = get_image_index_from_path(path);
	gtk_tree_path_free(path);
	if (index < 0 || index >= com.seq.number) return;
	fprintf(stdout, "toggle selection index = %d\n", index);

	sequence_list_change_selection(char_path, !com.seq.imgparam[index].incl);
	exclude_include_single_frame(index);
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
				GtkToggleButton *refframebutton = GTK_TOGGLE_BUTTON(lookup_widget("refframe2")); // invalidate toggle button, only effective if the frame is the first one in the selection
				g_signal_handlers_block_by_func(refframebutton, on_ref_frame_toggled, NULL);
				gtk_toggle_button_set_active(refframebutton, FALSE);
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
	unselect_select_frame_from_list((GtkTreeView *)user_data);
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

void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (!sequence_is_loaded())
		return;

	free_reference_image();
	get_list_store();
	GtkTreeModel *model = gtk_tree_view_get_model((GtkTreeView* ) user_data);
	GtkTreeSelection *selection = gtk_tree_view_get_selection((GtkTreeView* ) user_data);
	GList *references = get_row_references_of_selected_rows(selection, model);

	references = g_list_first(references);
	if (!references) return;

	GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)references->data);
	if (path) {
		if (!gtk_toggle_button_get_active(togglebutton)) {
			if (com.seq.reference_image == com.seq.current) {
				com.seq.reference_image = -1;
				com.seq.reference_image = sequence_find_refimage(&com.seq);
			}
		} else {
			com.seq.reference_image = com.seq.current;
			test_and_allocate_reference_image(-1);
			// a reference image should not be excluded to avoid confusion
			if (!com.seq.imgparam[com.seq.current].incl) {
				toggle_image_selection(com.seq.current, com.seq.current, FALSE);
			}
		}
		gtk_tree_path_free(path);
	}

	g_list_free(references);
	redraw(REDRAW_OVERLAY);
	sequence_list_change_reference();
	update_stack_interface(FALSE);// get stacking info and enable the Go button
	adjust_sellabel();	// reference image is named in the label
	writeseqfile(&com.seq);
	drawPlot();		// update plots
}

void on_draw_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (!sequence_is_loaded())
		return;
	redraw(REDRAW_OVERLAY);
}

/****************** modification of the list store (tree model) ******************/

void sequence_list_change_selection_index(int index_in_list, int real_index) {
	GtkTreePath *path = gtk_tree_path_new_from_indices(index_in_list, -1);
	if (path) {
		sequence_list_change_selection(gtk_tree_path_to_string(path), com.seq.imgparam[real_index].incl);
		gtk_tree_path_free(path);
	}
}

void sequence_list_change_current() {
	GtkTreeIter iter;
	gboolean valid;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gint index = 0;
		GValue value = G_VALUE_INIT;

		gtk_tree_model_get_value (GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
		index = g_value_get_int(&value) - 1;
		g_value_unset(&value);
		gtk_list_store_set(list_store, &iter, COLUMN_CURRENT,
				(index == com.seq.current) ? 800 : 400, -1);
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
	drawPlot();
}

void sequence_list_change_reference() {
	GtkTreeIter iter;
	gboolean valid;
	int color;

	color = (com.pref.gui.combo_theme == 0) ? 1 : 0;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gint index = 0;
		GValue value = G_VALUE_INIT;

		gtk_tree_model_get_value (GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
		index = g_value_get_int(&value) - 1;
		g_value_unset(&value);
		gtk_list_store_set(list_store, &iter,
						COLUMN_REFERENCE,
						(index == com.seq.reference_image) ?
						ref_bg_colour[color] : bg_colour[color], -1);
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
}

void clear_sequence_list() {
	g_signal_handlers_block_by_func(GTK_TREE_VIEW(lookup_widget("treeview1")), on_treeview1_cursor_changed, NULL);
	get_list_store();
	gtk_list_store_clear(list_store);
	g_signal_handlers_unblock_by_func(GTK_TREE_VIEW(lookup_widget("treeview1")), on_treeview1_cursor_changed, NULL);
}

void adjust_refimage(int n) {
	static GtkWidget *ref_butt2 = NULL, *treeview1 = NULL;
	if (ref_butt2 == NULL) {
		ref_butt2 = lookup_widget("refframe2");
		treeview1 = lookup_widget("treeview1");
	}

	g_signal_handlers_block_by_func(ref_butt2, on_ref_frame_toggled, treeview1);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ref_butt2), com.seq.reference_image == n);
	g_signal_handlers_unblock_by_func(ref_butt2, on_ref_frame_toggled, treeview1);
}

void sequence_list_select_row_from_index(int index, gboolean do_load_image) {
	GtkWidget *tree_view = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model;
	gboolean valid, found = FALSE;

	tree_view = lookup_widget("treeview1");
	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid && !found)
	{
		gint int_data;
		gtk_tree_model_get(GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &int_data, -1);
		if (int_data - 1 == index) {
			found = TRUE;
		} else {
			valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
		}
	}
	model = gtk_tree_view_get_model(GTK_TREE_VIEW(tree_view));
	if (!model) return;
	GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree_view)); //clear current selection
	gtk_tree_selection_unselect_all(selection);
	GtkTreePath *path = gtk_tree_model_get_path(model, &iter);
	if (!path) return;
	gtk_tree_selection_select_path(selection, path);
	gtk_tree_view_scroll_to_cell (GTK_TREE_VIEW(tree_view), path, NULL, FALSE, FALSE, FALSE);
	gtk_tree_path_free(path);

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
	gtk_button_set_image(GTK_BUTTON(GTK_TOGGLE_BUTTON(lookup_widget("drawframe_check"))), w);
	gtk_widget_show(w);
}

/*****************************************************************************
 *    Managing the list of sequence (the combo box) from the sequence tab
 *****************************************************************************/

#ifdef _WIN32
static int ListSequences(const gchar *sDir, const char *sequence_name_to_select,
		GtkComboBoxText *seqcombo, int *index_of_seq_to_load, gboolean *found) {
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
				gtk_combo_box_text_append_text(seqcombo, filename);
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
	GtkComboBoxText *seqcombo;
	int number_of_loaded_sequences = 0;
	int index_of_seq_to_load = -1;
	char *seqname = NULL;
	gboolean found = FALSE;

	// clear the previous list
	seqcombo = GTK_COMBO_BOX_TEXT(lookup_widget("sequence_list_combobox"));
	gtk_combo_box_text_remove_all(seqcombo);

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
	struct dirent **list;
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
				gtk_combo_box_text_append_text(seqcombo, filename);
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
		gtk_combo_box_popup(GTK_COMBO_BOX(seqcombo));
	} else if (index_of_seq_to_load >= 0)
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), index_of_seq_to_load);
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), 0);
	return 0;
}

gboolean on_treeview1_query_tooltip(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	GtkTreeIter iter;
	GtkTreeView *tree_view = GTK_TREE_VIEW(widget);
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreePath *path = NULL;

	if (!sequence_is_loaded()) return FALSE;

	if (!gtk_tree_view_get_tooltip_context(tree_view, &x, &y, keyboard_tip, &model, &path, &iter)) {
		return FALSE;
	}
	char buffer[512];
	fits fit = { 0 };
	gint real_index = get_real_index_from_index_in_list(model, &iter);
	if (seq_read_frame_metadata(&com.seq, real_index, &fit)) {
		return FALSE;
	}

	if (fit.keywords.filename[0] == '\0') return FALSE;

	g_snprintf(buffer, 511, "<b>Original filename:</b> %s", fit.keywords.filename);
	clearfits(&fit);
	gtk_tooltip_set_markup(tooltip, buffer);

	gtk_tree_view_set_tooltip_row(tree_view, tooltip, path);
	gtk_tree_path_free(path);

	return TRUE;
}
