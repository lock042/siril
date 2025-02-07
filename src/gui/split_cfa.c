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

#include "core/siril.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "algos/extraction.h"
#include "io/sequence.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

void on_split_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("split_cfa_dialog");
}

void on_split_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq = GTK_TOGGLE_BUTTON(lookup_widget("checkSplitCFASeq"));
	GtkEntry *entrySplitCFA = GTK_ENTRY(lookup_widget("entrySplitCFA"));
	gint method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_split_cfa_method")));

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {

		if (method == 0 || method == 2) {
			struct multi_output_data *args = calloc(1, sizeof(struct multi_output_data));
			set_cursor_waiting(TRUE);
			args->seq = &com.seq;
			args->user_data = calloc(1, sizeof(extraction_scaling));
			*(extraction_scaling *)args->user_data = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_haoiii_scaling")));
			args->seqEntry = strdup(gtk_entry_get_text(entrySplitCFA));

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
			args->scaling = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_haoiii_scaling")));
			args->seqEntry = strdup(gtk_entry_get_text(entrySplitCFA));

			switch (method) {
				case 1:
					if (args->seqEntry && args->seqEntry[0] == '\0')
						args->seqEntry = strdup("Ha_");
					apply_extractHa_to_sequence(args);
					break;
				case 3:
					if (args->seqEntry && args->seqEntry[0] == '\0')
						args->seqEntry = strdup("Green_");
					apply_extractGreen_to_sequence(args);
					break;
				default:
					siril_debug_print("unhandled case!\n");
					free(args->seqEntry);
					free(args);
			}
		}
	} else {
		int scaling = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_haoiii_scaling")));
		siril_debug_print("Scaling %d\n",scaling);
		switch (method) {
			case 0:
				process_split_cfa(0);
				break;
			case 1:
				extract_Ha(scaling);
				break;
			case 2:
				extract_HaOIII(scaling);
				break;
			case 3:
				process_extractGreen(0);
				break;
			default:
				fprintf(stderr, "unhandled case!\n");
		}
	}
}

void on_combo_split_cfa_method_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *w = lookup_widget("label10");
	GtkWidget *cb = lookup_widget("combo_haoiii_scaling");
	GtkWidget *cbl = lookup_widget("labelhaoiiiscaling");
	GtkWidget *txt = lookup_widget("entrySplitCFA");
	gint method = gtk_combo_box_get_active(box);

	gtk_widget_set_sensitive(cb, method == 2);
	gtk_widget_set_sensitive(cbl, method == 2);
	gtk_widget_set_sensitive(w, method != 2);
	gtk_widget_set_sensitive(txt, method != 2);
	switch (method) {
		case 0:
			gtk_entry_set_text(GTK_ENTRY(txt), "CFA_");
			break;
		case 1:
			gtk_entry_set_text(GTK_ENTRY(txt), "Ha_");
			break;
		case 3:
			gtk_entry_set_text(GTK_ENTRY(txt), "Green_");
			break;
	}
}

