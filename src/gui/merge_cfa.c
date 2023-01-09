/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include "algos/demosaicing.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/sequence_list.h"
#include "gui/PSF_list.h"
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

static gchar *f_cfa0 = NULL, *f_cfa1 = NULL, *f_cfa2 = NULL, *f_cfa3 = NULL;
fits cfa0 = { 0 }, cfa1 = { 0 }, cfa2 = { 0 }, cfa3 = { 0 };
static gboolean cfa0_loaded = FALSE, cfa1_loaded = FALSE, cfa2_loaded = FALSE, cfa3_loaded = FALSE;

void reset_controls() {
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa0")), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa1")), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa2")), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa3")), com.wd);
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa0")));
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa1")));
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa2")));
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa3")));
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("merge_cfa_pattern")), 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("merge_cfa_seqapply")), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("merge_cfa_seq_controls")), FALSE);

}

void on_merge_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	reset_controls();
	clearfits(&cfa0);
	clearfits(&cfa1);
	clearfits(&cfa2);
	clearfits(&cfa3);
	siril_close_dialog("merge_cfa_dialog");
}

void on_merge_cfa_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_controls();
}

void on_merge_cfa_show(GtkWidget *widget, gpointer user_data) {
	reset_controls();
}

void on_merge_cfa_seqapply_toggled(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_button_get_active(button))
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("merge_cfa_seq_controls")), TRUE);
	else
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("merge_cfa_seq_controls")), FALSE);
	gtk_window_resize(GTK_WINDOW(lookup_widget("merge_cfa_dialog")), 1, 1);
}

void on_merge_cfa_filechooser_CFA0_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	clearfits(&cfa0);
	f_cfa0 = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(f_cfa0, &cfa0, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		cfa0_loaded = FALSE;
		return;
	} else {
		cfa0_loaded = TRUE;
	}
}

void on_merge_cfa_filechooser_CFA1_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	clearfits(&cfa1);
	f_cfa1 = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(f_cfa1, &cfa1, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		cfa1_loaded = FALSE;
		return;
	} else {
		cfa1_loaded = TRUE;
	}
}
void on_merge_cfa_filechooser_CFA2_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	clearfits(&cfa2);
	f_cfa2 = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(f_cfa2, &cfa2, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		cfa2_loaded = FALSE;
		return;
	} else {
		cfa2_loaded = TRUE;
	}
}
void on_merge_cfa_filechooser_CFA3_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	clearfits(&cfa3);
	f_cfa3 = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(f_cfa3, &cfa3, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		cfa3_loaded = FALSE;
		return;
	} else {
		cfa3_loaded = TRUE;
	}
}

void apply_to_seq() {
	if (com.seq.type == SEQ_SER) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: sequence is SER"),
				_("Only FITS format is supported for CFA merging"));
		return;
	}
	GtkEntry *entryMergeCFAin = GTK_ENTRY(lookup_widget("entryMergeCFAin"));
	GtkEntry *entryMergeCFAout = GTK_ENTRY(lookup_widget("entryMergeCFAout"));
	struct merge_cfa_data *args = calloc(1, sizeof(struct merge_cfa_data));
	args->seqEntryIn = gtk_entry_get_text(entryMergeCFAin);
	if (args->seqEntryIn && args->seqEntryIn[0] == '\0')
		args->seqEntryIn = "CFA_";
	args->seqEntryOut = gtk_entry_get_text(entryMergeCFAout);
	if (args->seqEntryOut && args->seqEntryOut[0] == '\0')
		args->seqEntryOut = "mCFA_";
	GtkComboBox *combo_pattern = GTK_COMBO_BOX(lookup_widget("merge_cfa_pattern"));
	gint p = gtk_combo_box_get_active(combo_pattern);
	args->pattern = (sensor_pattern) p;
	args->seq = &com.seq;

	apply_mergecfa_to_sequence(args);
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_close_dialog("merge_cfa_dialog");
}

void apply_to_img() {
	GtkComboBox *combo_pattern = GTK_COMBO_BOX(lookup_widget("merge_cfa_pattern"));
	gint p = gtk_combo_box_get_active(combo_pattern);
	sensor_pattern pattern = (sensor_pattern) p;
	fits *out = NULL;
	if (!(cfa0_loaded && cfa1_loaded && cfa2_loaded && cfa3_loaded)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: images are not all loaded"),
				_("Merge CFA cannot proceed"));
	} else {
		gboolean x_compat = (cfa0.naxes[0] == cfa1.naxes[0] && cfa1.naxes[0] == cfa2.naxes[0] && cfa2.naxes[0] == cfa3.naxes[0]);
		gboolean y_compat = (cfa0.naxes[1] == cfa1.naxes[1] && cfa1.naxes[1] == cfa2.naxes[1] && cfa2.naxes[1] == cfa3.naxes[1]);
		gboolean c_compat = (cfa0.naxes[2] == cfa1.naxes[2] && cfa1.naxes[2] == cfa2.naxes[2] && cfa2.naxes[2] == cfa3.naxes[2] && cfa3.naxes[2] == 1);
		gboolean t_compat = (cfa0.type == cfa1.type && cfa1.type == cfa2.type && cfa2.type == cfa3.type);
		if (!(x_compat && y_compat && c_compat && t_compat)) {
			siril_log_color_message(_("Input files are incompatible (all must be mono with the same size and bit depth). Aborting...\n"), "red");
			if(!x_compat)
				siril_log_message(_("X dimensions incompatible\n"));
			if(!y_compat)
				siril_log_message(_("Y dimensions incompatible\n"));
			if(!c_compat)
				siril_log_message(_("Channels not all mono\n"));
			if(!t_compat)
				siril_log_message(_("Input files not all the same bit depth\n"));
			return;
		} else {
			set_cursor_waiting(TRUE);
			out = merge_cfa(&cfa0, &cfa1, &cfa2, &cfa3, pattern);
			siril_log_message("Bayer pattern produced: 1 layer, %dx%d pixels\n", out->rx, out->ry);
			close_single_image();
			copyfits(out, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
			clearfits(out);
			clear_stars_list(TRUE);
			com.seq.current = UNRELATED_IMAGE;
			if (!create_uniq_from_gfit(strdup(_("Unsaved Bayer pattern merge")), FALSE))
				com.uniq->comment = strdup(_("Bayer pattern merge"));
			open_single_image_from_gfit();
			initialize_display_mode();
			update_zoom_label();
			display_filename();
			set_precision_switch();
			sliders_mode_set_state(gui.sliders);
			init_layers_hi_and_lo_values(MIPSLOHI);
			set_cutoff_sliders_max_values();
			set_cutoff_sliders_values();
			set_display_mode();
			redraw(REMAP_ALL);
			sequence_list_change_current();
			set_cursor_waiting(FALSE);
			reset_controls();
			control_window_switch_to_tab(OUTPUT_LOGS);
			siril_close_dialog("merge_cfa_dialog");
		}
	}
}

void on_merge_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *to_seq = GTK_TOGGLE_BUTTON(lookup_widget("merge_cfa_seqapply"));
	if (gtk_toggle_button_get_active(to_seq) && sequence_is_loaded())
		apply_to_seq();
	else
		apply_to_img();
}


