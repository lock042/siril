/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"

void on_menuitem_arithm_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	siril_open_dialog("arithm_dialog");
}

static int get_oper() {
	return gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_arithm_oper")));
}

static gchar *get_operand_filename() {
	GtkFileChooser *operand = GTK_FILE_CHOOSER(lookup_widget("arithm_file_chooser"));
	return gtk_file_chooser_get_filename(operand);
}

void on_arithm_apply_clicked(GtkButton *button, gpointer user_data) {
	fits *operand;
	gchar *filename = get_operand_filename();
	if (!filename) return;

	operand = calloc(1, sizeof(fits));
	if (readfits(filename, operand, NULL, gfit.type == DATA_FLOAT)) {
		g_free(filename);
		return;
	}
	g_free(filename);
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkArithmSeq")))
			&& sequence_is_loaded()) {
		struct arithm_data *args = malloc(sizeof(struct arithm_data));

		args->oper = get_oper();
		args->operand = operand;

		set_cursor_waiting(TRUE);

		args->seqEntry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryArithm")));
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "arithm_";
		args->seq = &com.seq;
		apply_arithm_extraction_to_sequence(args);
	} else {
		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, "Processing: Arithmetic");

		imoper(&gfit, operand, get_oper(), !com.pref.force_to_16bit);

		clearfits(operand);
		invalidate_stats_from_fit(&gfit);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);

		set_cursor_waiting(FALSE);
	}
}

void on_arithm_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("arithm_dialog");
}
