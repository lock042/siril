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
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"


static GtkCheckButton *split_cfa_seq_btn = NULL;
static GtkEntry *split_cfa_entry = NULL;
static GtkDropDown *split_cfa_method_combo = NULL;
static GtkDropDown *split_haoiii_scaling_combo = NULL;
static GtkWidget *split_label10 = NULL;
static GtkWidget *split_labelhaoiiiname = NULL;

static void split_cfa_init_statics(void) {
	if (split_cfa_seq_btn) return;
	split_cfa_seq_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkSplitCFASeq"));
	split_cfa_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entrySplitCFA"));
	split_cfa_method_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_split_cfa_method"));
	split_haoiii_scaling_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_haoiii_scaling"));
	split_label10 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label10"));
	split_labelhaoiiiname = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelhaoiiiscaling"));
}

void on_split_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("split_cfa_dialog");
}

void on_split_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	split_cfa_init_statics();
	GtkCheckButton *seq = split_cfa_seq_btn;
	GtkEntry *entrySplitCFA = split_cfa_entry;
	gint method = gtk_drop_down_get_selected(split_cfa_method_combo);

	if (siril_toggle_get_active(GTK_WIDGET(seq)) && sequence_is_loaded()) {

		if (method == 0 || method == 2) {
			struct multi_output_data *args = calloc(1, sizeof(struct multi_output_data));
			set_cursor_waiting(TRUE);
			args->seq = &com.seq;
			args->user_data = calloc(1, sizeof(extraction_scaling));
			*(extraction_scaling *)args->user_data = gtk_drop_down_get_selected(split_haoiii_scaling_combo);
			args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(entrySplitCFA)));

			switch (method) {
				case 0:
					args->n = 4;
					args->prefixes = calloc(5, sizeof(const char*));
					if (!args->seqEntry || (args->seqEntry && args->seqEntry[0] == '\0')) {
						free(args->seqEntry);
						args->seqEntry = strdup("CFA");
					}
					size_t len = strlen(args->seqEntry);
					// Strip any trailing '_' in order to insert the frame number
					if (len > 0 && args->seqEntry[len - 1] == '_') {
						args->seqEntry[len - 1] = '\0';  // Replace trailing '_' with null terminator
					}
					for (int i = 0 ; i < 4 ; i++) {
						args->prefixes[i] = g_strdup_printf("%s%d_", args->seqEntry, i);
					}
					apply_split_cfa_to_sequence(args);
					break;
				case 2:
					args->n = 2;
					args->prefixes = calloc(3, sizeof(const char*));
					args->prefixes[0] = g_strdup("Ha_");
					args->prefixes[1] = g_strdup("OIII_");
					apply_extractHaOIII_to_sequence(args);
					break;
			}
		} else {
			struct simple_extract_data *args = calloc(1, sizeof(struct simple_extract_data));
			set_cursor_waiting(TRUE);
			args->seq = &com.seq;
			args->scaling = gtk_drop_down_get_selected(split_haoiii_scaling_combo);
			args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(entrySplitCFA)));

			switch (method) {
				case 1:
					if (args->seqEntry && args->seqEntry[0] == '\0') {
						free(args->seqEntry);
						args->seqEntry = strdup("Ha_");
					}
					apply_extractHa_to_sequence(args);
					break;
				case 3:
					if (args->seqEntry && args->seqEntry[0] == '\0') {
						free(args->seqEntry);
						args->seqEntry = strdup("Green_");
					}
					apply_extractGreen_to_sequence(args);
					break;
				default:
					siril_log_debug("unhandled case!\n");
					free(args->seqEntry);
					free(args);
			}
		}
	} else {
		int scaling = gtk_drop_down_get_selected(split_haoiii_scaling_combo);
		siril_log_debug("Scaling %d\n", scaling);

		/* Compute base filename before starting the thread */
		gchar *filename = NULL;
		if (sequence_is_loaded() && !single_image_is_loaded()) {
			filename = g_path_get_basename(com.seq.seqname);
		} else if (com.uniq && com.uniq->filename) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}

		struct cfa_extract_args *cfa_args = calloc(1, sizeof(struct cfa_extract_args));
		cfa_args->destroy_fn = free_cfa_extract_args;
		cfa_args->operation = method;
		cfa_args->scaling = scaling;

		switch (method) {
			case 0: /* split_cfa */
				cfa_args->channel[0] = g_strdup_printf("CFA0_%s%s", filename, com.pref.ext);
				cfa_args->channel[1] = g_strdup_printf("CFA1_%s%s", filename, com.pref.ext);
				cfa_args->channel[2] = g_strdup_printf("CFA2_%s%s", filename, com.pref.ext);
				cfa_args->channel[3] = g_strdup_printf("CFA3_%s%s", filename, com.pref.ext);
				break;
			case 1: /* extractHa */
				cfa_args->pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
				if (cfa_args->pattern < BAYER_FILTER_MIN || cfa_args->pattern > BAYER_FILTER_MAX) {
					siril_log_error(_("This image does not have a Bayer CFA pattern, cannot extract Ha.\n"));
					free_cfa_extract_args(cfa_args);
					g_free(filename);
					return;
				}
				cfa_args->channel[0] = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
				break;
			case 2: /* extractHaOIII */
				cfa_args->pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
				if (cfa_args->pattern < BAYER_FILTER_MIN || cfa_args->pattern > BAYER_FILTER_MAX) {
					siril_log_error(_("This image does not have a Bayer CFA pattern, cannot extract Ha/OIII channels.\n"));
					free_cfa_extract_args(cfa_args);
					g_free(filename);
					return;
				}
				cfa_args->channel[0] = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
				cfa_args->channel[1] = g_strdup_printf("OIII_%s%s", filename, com.pref.ext);
				break;
			case 3: /* extractGreen */
				cfa_args->pattern = get_validated_cfa_pattern(gfit, FALSE, FALSE);
				if (cfa_args->pattern < BAYER_FILTER_MIN || cfa_args->pattern > BAYER_FILTER_MAX) {
					siril_log_error(_("This image does not have a Bayer CFA pattern, cannot extract green channel.\n"));
					free_cfa_extract_args(cfa_args);
					g_free(filename);
					return;
				}
				cfa_args->channel[0] = g_strdup_printf("Green_%s%s", filename, com.pref.ext);
				break;
			default:
				siril_log_debug("unhandled case!\n");
				free_cfa_extract_args(cfa_args);
				g_free(filename);
				return;
		}
		g_free(filename);

		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		args->fit = gfit;
		args->image_hook = cfa_extract_image_hook;
		args->idle_function = end_generic_image_reset_cursor;
		args->description = _("CFA Extraction");
		args->verbose = TRUE;
		args->custom_undo = TRUE;
		args->user = cfa_args;
		set_cursor_waiting(TRUE);
		if (!start_in_new_thread(generic_image_worker, args)) {
			free_generic_img_args(args);
			set_cursor_waiting(FALSE);
		}
	}
}

void on_combo_split_cfa_method_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	split_cfa_init_statics();
	GtkWidget *w = split_label10;
	GtkWidget *cb = GTK_WIDGET(split_haoiii_scaling_combo);
	GtkWidget *cbl = split_labelhaoiiiname;
	GtkWidget *txt = GTK_WIDGET(split_cfa_entry);
	gint method = gtk_drop_down_get_selected(box);

	gtk_widget_set_sensitive(cb, method == 2);
	gtk_widget_set_sensitive(cbl, method == 2);
	gtk_widget_set_sensitive(w, method != 2);
	gtk_widget_set_sensitive(txt, method != 2);
	switch (method) {
		case 0:
			gtk_editable_set_text(GTK_EDITABLE(txt), "CFA_");
			break;
		case 1:
			gtk_editable_set_text(GTK_EDITABLE(txt), "Ha_");
			break;
		case 3:
			gtk_editable_set_text(GTK_EDITABLE(txt), "Green_");
			break;
	}
}

