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
#include "core/command_line_processor.h" // for load_sequence - TODO: move this to sequence.c
#include "core/processing.h"
#include "core/processing_thread.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"
#include "algos/siril_wcs.h"
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

static GtkFileChooser *mcfa_chooser0 = NULL;
static GtkFileChooser *mcfa_chooser1 = NULL;
static GtkFileChooser *mcfa_chooser2 = NULL;
static GtkFileChooser *mcfa_chooser3 = NULL;
static GtkComboBox *mcfa_pattern_combo = NULL;
static GtkToggleButton *mcfa_seqapply_btn = NULL;
static GtkWidget *mcfa_seq_controls = NULL;
static GtkEntry *mcfa_entry_in = NULL;
static GtkEntry *mcfa_entry_out = NULL;
static GtkWindow *mcfa_dialog_window = NULL;

static void merge_cfa_init_statics(void) {
	if (mcfa_chooser0) return;
	mcfa_chooser0 = GTK_FILE_CHOOSER(gtk_builder_get_object(gui.builder, "filechooser_cfa0"));
	mcfa_chooser1 = GTK_FILE_CHOOSER(gtk_builder_get_object(gui.builder, "filechooser_cfa1"));
	mcfa_chooser2 = GTK_FILE_CHOOSER(gtk_builder_get_object(gui.builder, "filechooser_cfa2"));
	mcfa_chooser3 = GTK_FILE_CHOOSER(gtk_builder_get_object(gui.builder, "filechooser_cfa3"));
	mcfa_pattern_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "merge_cfa_pattern"));
	mcfa_seqapply_btn = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "merge_cfa_seqapply"));
	mcfa_seq_controls = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_cfa_seq_controls"));
	mcfa_entry_in = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryMergeCFAin"));
	mcfa_entry_out = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryMergeCFAout"));
	mcfa_dialog_window = GTK_WINDOW(gtk_builder_get_object(gui.builder, "merge_cfa_dialog"));
}

void reset_controls() {
	merge_cfa_init_statics();
	gtk_file_chooser_set_current_folder(mcfa_chooser0, com.wd);
	gtk_file_chooser_set_current_folder(mcfa_chooser1, com.wd);
	gtk_file_chooser_set_current_folder(mcfa_chooser2, com.wd);
	gtk_file_chooser_set_current_folder(mcfa_chooser3, com.wd);
	gtk_file_chooser_unselect_all(mcfa_chooser0);
	gtk_file_chooser_unselect_all(mcfa_chooser1);
	gtk_file_chooser_unselect_all(mcfa_chooser2);
	gtk_file_chooser_unselect_all(mcfa_chooser3);
	gtk_combo_box_set_active(mcfa_pattern_combo, 0);
	gtk_toggle_button_set_active(mcfa_seqapply_btn, FALSE);
	gtk_widget_set_visible(mcfa_seq_controls, FALSE);
	gtk_entry_set_text(mcfa_entry_in, "CFA");
	gtk_entry_set_text(mcfa_entry_out, "mCFA_");
}

static void close_everything() {
	reset_controls();
	clearfits(&cfa0);
	clearfits(&cfa1);
	clearfits(&cfa2);
	clearfits(&cfa3);
	siril_close_dialog("merge_cfa_dialog");
}
void on_merge_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	close_everything();
}

gboolean merge_cfa_hide_on_delete(GtkWidget *widget) {
	close_everything();
	return TRUE;
}

void on_merge_cfa_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_controls();
}

void on_merge_cfa_show(GtkWidget *widget, gpointer user_data) {
	reset_controls();
}

void on_merge_cfa_seqapply_toggled(GtkToggleButton *button, gpointer user_data) {
	merge_cfa_init_statics();
	gtk_widget_set_visible(mcfa_seq_controls, gtk_toggle_button_get_active(button));
	gtk_window_resize(mcfa_dialog_window, 1, 1);
}

void on_merge_cfa_filechooser_CFA0_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	clearfits(&cfa0);
	f_cfa0 = siril_file_chooser_get_filename(filechooser);
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
	f_cfa1 = siril_file_chooser_get_filename(filechooser);
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
	f_cfa2 = siril_file_chooser_get_filename(filechooser);
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
	f_cfa3 = siril_file_chooser_get_filename(filechooser);
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

gchar* replace_string_number(const gchar* input, const gchar* string, guint n) {
	gchar* result;
	GRegex* regex;
	gchar* pattern;
	gchar* replacement;

	// Escape the string to handle special regex characters
	gchar* escaped_string = g_regex_escape_string(string, -1);

	// Create a regex pattern to match string followed by any number
	pattern = g_strdup_printf("%s(\\d+)", escaped_string);
	regex = g_regex_new(pattern, 0, 0, NULL);

	// Create the replacement string
	replacement = g_strdup_printf("%s%u", string, n);

	// Replace the first occurrence of "stringm" with "stringn"
	result = g_regex_replace(regex, input, -1, 0, replacement, 0, NULL);

	// Clean up
	g_regex_unref(regex);
	g_free(pattern);
	g_free(replacement);
	g_free(escaped_string);

	return result;
}

gint find_substring_number(const gchar* input, const gchar* substring) {
	GRegex* regex;
	GMatchInfo* match_info;
	gint result = -1;
	gchar* pattern;

	// Escape the substring to handle special regex characters
	gchar* escaped_substring = g_regex_escape_string(substring, -1);

	// Create a regex pattern to match substring followed by a digit
	pattern = g_strdup_printf("%s(\\d)", escaped_substring);
	regex = g_regex_new(pattern, 0, 0, NULL);

	if (g_regex_match(regex, input, 0, &match_info)) {
		gchar* match_string = g_match_info_fetch(match_info, 1);
		if (match_string) {
			result = g_ascii_strtoll(match_string, NULL, 10);
			g_free(match_string);
		}
	}

	g_match_info_free(match_info);
	g_regex_unref(regex);
	g_free(pattern);
	g_free(escaped_substring);

	return result;
}

void apply_to_seq() {
	sequence *seq0 = NULL, *seq1 = NULL, *seq2 = NULL, *seq3 = NULL;
	struct merge_cfa_data *args = NULL;
	gchar *seqname0 = NULL, *seqname1 = NULL, *seqname2 = NULL, *seqname3 = NULL;

	// Get input sequence marker
	merge_cfa_init_statics();
	GtkEntry *entryMergeCFAin = mcfa_entry_in;
	GtkEntry *entryMergeCFAout = mcfa_entry_out;
	const gchar *seqmarker = gtk_entry_get_text(entryMergeCFAin);

	int comseqnumber = find_substring_number(com.seq.seqname, seqmarker);
	if (comseqnumber == -1) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
							 _("Error: input sequence CFA marker not found"));
		return;
	}

	siril_log_message(_("Identified the loaded sequence as CFA%d\n"), comseqnumber);

	// Allocate sequence names
	seqname0 = replace_string_number(com.seq.seqname, seqmarker, 0);
	seqname1 = replace_string_number(com.seq.seqname, seqmarker, 1);
	seqname2 = replace_string_number(com.seq.seqname, seqmarker, 2);
	seqname3 = replace_string_number(com.seq.seqname, seqmarker, 3);

	// Allocate merge data structure
	args = calloc(1, sizeof(struct merge_cfa_data));
	if (!args) {
		goto cleanup;
	}

	const gchar *outText = gtk_entry_get_text(entryMergeCFAout);
	args->seqEntryOut = strdup(outText && outText[0] ? outText : "mCFA_");
	if (!args->seqEntryOut) {
		goto cleanup;
	}

	GtkComboBox *combo_pattern = mcfa_pattern_combo;
	args->pattern = (sensor_pattern)gtk_combo_box_get_active(combo_pattern);

	// Load sequences
	if (comseqnumber != 0) {
		seq0 = load_sequence(seqname0, NULL);
		if (!seq0 || seq0->nb_layers > 1) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
								 _("Error: sequence has more than one channel, this is for mono (CFA subchannel) sequences"));
			goto cleanup;
		}
	} else {
		seq0 = &com.seq;
	}

	if (comseqnumber != 1) {
		seq1 = load_sequence(seqname1, NULL);
		if (!seq1 || seq1->nb_layers > 1) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
								 _("Error: sequence has more than one channel, this is for mono (CFA subchannel) sequences"));
			goto cleanup;
		}
	} else {
		seq1 = &com.seq;
	}

	if (comseqnumber != 2) {
		seq2 = load_sequence(seqname2, NULL);
		if (!seq2 || seq2->nb_layers > 1) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
								 _("Error: sequence has more than one channel, this is for mono (CFA subchannel) sequences"));
			goto cleanup;
		}
	} else {
		seq2 = &com.seq;
	}

	if (comseqnumber != 3) {
		seq3 = load_sequence(seqname3, NULL);
		if (!seq3 || seq3->nb_layers > 1) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
								 _("Error: sequence has more than one channel, this is for mono (CFA subchannel) sequences"));
			goto cleanup;
		}
	} else {
		seq3 = &com.seq;
	}

	// Verify sequence compatibility
	if ((seq0->rx != seq1->rx || seq0->rx != seq2->rx || seq0->rx != seq3->rx) ||
				(seq0->ry != seq1->ry || seq0->ry != seq2->ry || seq0->ry != seq3->ry) ||
				(seq0->nb_layers != seq1->nb_layers || seq0->nb_layers != seq2->nb_layers ||
				seq0->nb_layers != seq3->nb_layers) ||
				(seq0->bitpix != seq1->bitpix || seq0->bitpix != seq2->bitpix ||
				seq0->bitpix != seq3->bitpix) ||
				(seq0->number != seq1->number || seq0->number != seq2->number ||
				seq0->number != seq3->number)) {
		siril_log_error(_("Error: sequences don't match (dimensions, bitdepth, "
		"number of images must all be the same)\n"));
	goto cleanup;
	}

	// Set sequences in args
	args->seq0 = seq0;
	args->seq1 = seq1;
	args->seq2 = seq2;
	args->seq3 = seq3;

	// Apply merge and cleanup
	apply_mergecfa_to_sequence(args);
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_close_dialog("merge_cfa_dialog");

	// Normal exit - don't free the sequences as they're now owned by apply_mergecfa_to_sequence
	g_free(seqname0);
	g_free(seqname1);
	g_free(seqname2);
	g_free(seqname3);
	return;

	cleanup:
	// Free all allocated resources
	g_free(seqname0);
	g_free(seqname1);
	g_free(seqname2);
	g_free(seqname3);

	if (args) {
		free(args->seqEntryOut);
		free(args);
	}

	if (seq0 && !check_seq_is_comseq(seq0))
		free_sequence(seq0, TRUE);
	if (seq1 && !check_seq_is_comseq(seq1))
		free_sequence(seq1, TRUE);
	if (seq2 && !check_seq_is_comseq(seq2))
		free_sequence(seq2, TRUE);
	if (seq3 && !check_seq_is_comseq(seq3))
		free_sequence(seq3, TRUE);
}

/* Worker data for single-image merge */
struct merge_cfa_img_data {
	gchar *f_cfa0, *f_cfa1, *f_cfa2, *f_cfa3;
	sensor_pattern pattern;
	gboolean success;
};

static void free_merge_cfa_img_data(struct merge_cfa_img_data *data) {
	g_free(data->f_cfa0);
	g_free(data->f_cfa1);
	g_free(data->f_cfa2);
	g_free(data->f_cfa3);
	free(data);
}

static gboolean merge_cfa_img_idle(gpointer p) {
	struct merge_cfa_img_data *data = (struct merge_cfa_img_data *)p;
	stop_processing_thread();
	if (data->success) {
		g_rw_lock_reader_lock(&gfit->rwlock);
		gui_function(open_single_image_from_gfit, NULL);
		initialize_display_mode();
		update_zoom_label();
		display_filename();
		gui_function(set_precision_switch, NULL);
		sliders_mode_set_state(gui.sliders);
		set_cutoff_sliders_max_values();
		g_rw_lock_reader_unlock(&gfit->rwlock);
		set_cutoff_sliders_values();
		set_display_mode();
		redraw(REMAP_ALL);
		sequence_list_change_current();
		reset_controls();
		control_window_switch_to_tab(OUTPUT_LOGS);
		siril_close_dialog("merge_cfa_dialog");
	}
	set_cursor_waiting(FALSE);
	free_merge_cfa_img_data(data);
	return FALSE;
}

static gpointer merge_cfa_img_worker(gpointer p) {
	struct merge_cfa_img_data *data = (struct merge_cfa_img_data *)p;
	fits c0 = { 0 }, c1 = { 0 }, c2 = { 0 }, c3 = { 0 };
	int retval = readfits(data->f_cfa0, &c0, NULL, FALSE);
	retval += readfits(data->f_cfa1, &c1, NULL, FALSE);
	retval += readfits(data->f_cfa2, &c2, NULL, FALSE);
	retval += readfits(data->f_cfa3, &c3, NULL, FALSE);
	if (retval) {
		siril_log_error(_("Error loading files!\n"));
		clearfits(&c0); clearfits(&c1); clearfits(&c2); clearfits(&c3);
		data->success = FALSE;
		siril_add_idle(merge_cfa_img_idle, data);
		return GINT_TO_POINTER(1);
	}
	fits *out = merge_cfa(&c0, &c1, &c2, &c3, data->pattern);
	clearfits(&c0); clearfits(&c1); clearfits(&c2); clearfits(&c3);
	if (!out) {
		siril_log_error(_("Merge CFA failed.\n"));
		data->success = FALSE;
		siril_add_idle(merge_cfa_img_idle, data);
		return GINT_TO_POINTER(1);
	}
	siril_log_message("Bayer pattern produced: 1 layer, %dx%d pixels\n", out->rx, out->ry);
	g_rw_lock_writer_lock(&gfit->rwlock);
	close_single_image();
	copyfits(out, gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copy_fits_metadata(out, gfit);
	update_sampling_information(gfit, 0.5f);
	update_bayer_pattern_information(gfit, data->pattern);
	free_wcs(gfit);
	update_fits_header(gfit);
	g_rw_lock_writer_unlock(&gfit->rwlock);
	clearfits(out);
	free(out);
	clear_stars_list(TRUE);
	com.seq.current = UNRELATED_IMAGE;
	if (!create_uniq_from_gfit(strdup(_("Unsaved Bayer pattern merge")), FALSE))
		com.uniq->comment = strdup(_("Bayer pattern merge"));
	if (!com.headless) {
		notify_gfit_data_modified();
		init_layers_hi_and_lo_values(MIPSLOHI);
	}
	data->success = TRUE;
	siril_add_idle(merge_cfa_img_idle, data);
	return GINT_TO_POINTER(0);
}

void apply_to_img() {
	GtkComboBox *combo_pattern = mcfa_pattern_combo;
	gint p = gtk_combo_box_get_active(combo_pattern);
	sensor_pattern pattern = (sensor_pattern) p;

	if (!(cfa0_loaded && cfa1_loaded && cfa2_loaded && cfa3_loaded)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error: images are not all loaded"),
				_("Merge CFA cannot proceed"));
		return;
	}

	gboolean x_compat = (cfa0.naxes[0] == cfa1.naxes[0] && cfa1.naxes[0] == cfa2.naxes[0] && cfa2.naxes[0] == cfa3.naxes[0]);
	gboolean y_compat = (cfa0.naxes[1] == cfa1.naxes[1] && cfa1.naxes[1] == cfa2.naxes[1] && cfa2.naxes[1] == cfa3.naxes[1]);
	gboolean c_compat = (cfa0.naxes[2] == cfa1.naxes[2] && cfa1.naxes[2] == cfa2.naxes[2] && cfa2.naxes[2] == cfa3.naxes[2] && cfa3.naxes[2] == 1);
	gboolean t_compat = (cfa0.type == cfa1.type && cfa1.type == cfa2.type && cfa2.type == cfa3.type);
	if (!(x_compat && y_compat && c_compat && t_compat)) {
		siril_log_error(_("Input files are incompatible (all must be mono with the same size and bit depth). Aborting...\n"));
		if (!x_compat)
			siril_log_message(_("X dimensions incompatible\n"));
		if (!y_compat)
			siril_log_message(_("Y dimensions incompatible\n"));
		if (!c_compat)
			siril_log_message(_("Channels not all mono\n"));
		if (!t_compat)
			siril_log_message(_("Input files not all the same bit depth\n"));
		return;
	}

	/* Snapshot filenames before handing off to the worker */
	struct merge_cfa_img_data *args = calloc(1, sizeof(struct merge_cfa_img_data));
	if (!args)
		return;
	args->f_cfa0 = g_strdup(f_cfa0);
	args->f_cfa1 = g_strdup(f_cfa1);
	args->f_cfa2 = g_strdup(f_cfa2);
	args->f_cfa3 = g_strdup(f_cfa3);
	args->pattern = pattern;

	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(merge_cfa_img_worker, args)) {
		free_merge_cfa_img_data(args);
		set_cursor_waiting(FALSE);
	}
}

void on_merge_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	merge_cfa_init_statics();
	GtkToggleButton *to_seq = mcfa_seqapply_btn;
	if (gtk_toggle_button_get_active(to_seq) && sequence_is_loaded())
		apply_to_seq();
	else
		apply_to_img();
}


