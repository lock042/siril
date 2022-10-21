/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "algos/extraction.h"
#include "io/sequence.h"
#include "gui/dialogs.h"
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
		struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

		set_cursor_waiting(TRUE);
		args->seq = &com.seq;
		args->scaling = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_haoiii_scaling")));
		args->scaling = 0; // Temporary until the glade element is added (TODO)
		args->seqEntry = gtk_entry_get_text(entrySplitCFA);
		switch (method) {
			case 0:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "CFA_";
				apply_split_cfa_to_sequence(args);
				break;
			case 1:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "Ha_";
				apply_extractHa_to_sequence(args);
				break;
			case 2:
				apply_extractHaOIII_to_sequence(args);
				break;
			case 3:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "Green_";
				apply_extractGreen_to_sequence(args);
				break;
			default:
				fprintf(stderr, "unhandled case!\n");
		}
	} else {
		switch (method) {
			case 0:
				process_split_cfa(0);
				break;
			case 1:
				process_extractHa(0);
				break;
			case 2:
				process_extractHaOIII(0);
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
	GtkWidget *txt = lookup_widget("entrySplitCFA");
	gint method = gtk_combo_box_get_active(box);

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

