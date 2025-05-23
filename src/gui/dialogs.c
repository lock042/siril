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
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/message_dialog.h"
#include "gui/curves.h"
#include "algos/background_extraction.h"
#include "filters/asinh.h"
#include "filters/epf.h"
#include "filters/clahe.h"
#include "filters/median.h"
#include "filters/saturation.h"
#include "filters/unpurple.h"
#include "filters/wavelets.h"

#include "gui/newdeconv.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "gui/dialogs.h"
#include "gui/dialog_preview.h"
#include "nina_light_curve.h"
#include "compstars.h"

static gboolean dialog_is_opened = FALSE;
static gboolean processing_dialog_is_opened = FALSE;

static const SirilDialogEntry entries[] =
{
	{"aavso_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"annotate_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"asinh_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_asinh_cancel},
	{"astrometry_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"denoise_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, close_denoise},
	{"background_extraction_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_background_cancel},
	{"binxy_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"canon_fixbanding_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"ccm_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"CLAHE_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_clahe_cancel},
	{"composition_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"compstars", get_compstars_dialog, OTHER_DIALOG, FALSE, NULL},
	{"color_calibration", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"cosmetic_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"crop_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"curves_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_curves_cancel},
	{"cut_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"cut_coords_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"cut_spectroscopy_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"bdeconv_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, close_deconv},
	{"dialog_FFT", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"dialog_star_remix", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"edge_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"epf_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_epf_cancel},
	{"extract_channel_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"extract_wavelets_layers_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"file_information", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"histogram_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_histo_cancel},
	{"keywords_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"icc_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"astrometry_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"linearmatch_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"Median_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, median_close},
	{"merge_cfa_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"mouse_actions_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"nina_light_curve", get_nina_lc_dialog, OTHER_DIALOG, FALSE, NULL},
	{"pixel_math_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"resample_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"rgradient_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"rotation_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"s_pcc_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"satu_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_satu_cancel},
	{"SCNR_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"settings_window", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"seqlist_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"split_cfa_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"starnet_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"stars_list_window", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"StatWindow", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"unpurple_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_unpurple_cancel},
	{"wavelets_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_wavelets_cancel}
};

static SirilDialogEntry get_entry_by_id(gchar *id) {
	int i;
	SirilDialogEntry empty = {NULL, NULL, NO_DIALOG, FALSE, NULL};

	for (i = 0; i < G_N_ELEMENTS(entries); i++) {
		if (!g_ascii_strcasecmp(id, entries[i].identifier)) {
			return entries[i];
		}
	}
	return empty;
}

static GtkWidget *get_widget_by_id(gchar *id) {
	SirilDialogEntry entry = get_entry_by_id(id);
	if (entry.get_window)
		return entry.get_window();
	return lookup_widget(entry.identifier);
}

static GtkWidget *get_widget_by_index(int index) {
	SirilDialogEntry entry = entries[index];
	if (entry.get_window)
		return entry.get_window();
	return lookup_widget(entry.identifier);
}

static gboolean check_escape(GtkWidget *widget, GdkEventKey *event,
		gpointer data) {
	if (event->keyval == GDK_KEY_Escape) {
		gtk_widget_hide(widget);
		return TRUE;
	}
	return FALSE;
}

void siril_open_dialog(gchar *id) {
	gint x, y;
	gboolean win_already_shown = FALSE;
	SirilDialogEntry entry = get_entry_by_id(id);
	GtkWindow *win = GTK_WINDOW(get_widget_by_id(id));

	// We cannot open the dialog if a python script claims the thread, to prevent conflict over
	// updating gfit
	if (entry.type == IMAGE_PROCESSING_DIALOG && com.python_claims_thread) {
		queue_error_message_dialog(_("Error"), _("Error: cannot open an image processing dialog while a Python script claims the processing thread. "
									"Wait for the Python script to finish processing first."));
		return;
	}

	if (entry.type != INFORMATION_DIALOG) {
		for (int i = 0; i < G_N_ELEMENTS(entries); i++) {
			GtkWidget *w = get_widget_by_index(i);
			if (gtk_widget_get_visible(w) && entries[i].type != INFORMATION_DIALOG) {
				win_already_shown = TRUE;
				gtk_window_get_position(GTK_WINDOW(w), &x, &y);
				if (entries[i].has_preview)
					entries[i].apply_function();
				gtk_widget_hide(w);
				break;
			}
		}
	}
	if (entry.type == SEARCH_ENTRY_DIALOG) {
		g_signal_connect(win, "key_press_event", G_CALLBACK(check_escape), NULL);
	}

	if (win_already_shown && x >= 0 && y >= 0) {
		gtk_window_move(win, x, y);
	} else {
		gtk_window_set_position(win, GTK_WIN_POS_CENTER);
	}
	gtk_window_set_type_hint(win, GDK_WINDOW_TYPE_HINT_DIALOG);
	gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget("control_window")));
	gtk_window_present_with_time(win, GDK_CURRENT_TIME);
	dialog_is_opened = TRUE;
	if (entry.type == IMAGE_PROCESSING_DIALOG) {
		siril_debug_print("### Opening imgproc dialog: %s\n", entry.identifier);
		processing_dialog_is_opened = TRUE;
	}
}

void siril_close_dialog(gchar *id) {
	gtk_widget_hide(get_widget_by_id(id));
	dialog_is_opened = FALSE;
	SirilDialogEntry entry = get_entry_by_id(id);
	if (entry.type == IMAGE_PROCESSING_DIALOG) {
		siril_debug_print("### Closing imgproc dialog: %s\n", entry.identifier);
		processing_dialog_is_opened = FALSE;
	}
}

void siril_close_preview_dialogs() {
	for (int i = 0; i < G_N_ELEMENTS(entries); i++) {
		GtkWidget *w = get_widget_by_index(i);
		if (gtk_widget_get_visible(w) && (entries[i].has_preview)) {
			entries[i].apply_function();
			gtk_widget_hide(w);
		}
	}
}

gboolean siril_widget_hide_on_delete(GtkWidget *widget) {
    dialog_is_opened = FALSE;
    gtk_widget_hide(widget);
    return TRUE;
}

gboolean is_a_dialog_opened() {
	return dialog_is_opened;
}

gboolean is_an_image_processing_dialog_opened() {
	return processing_dialog_is_opened;
}

/************ file chooser ************/

void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern, gboolean set_default) {
	gchar **patterns;
	gint i;

	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	patterns = g_strsplit(pattern, ";", -1);
	for (i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	if (set_default)
		gtk_file_chooser_set_filter(file_chooser, f);
}

SirilWidget *siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action) {
	gchar *title;
	SirilWidget *w;
	if (action == GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER) {
		title = g_strdup(_("Select Folder"));
	} else {
		title = g_strdup(_("Open File"));
	}
	w = gtk_file_chooser_dialog_new(title, parent, action, _("_Cancel"),
			GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);
	g_free(title);
	return w;
}

SirilWidget *siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Add Files"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Add"), GTK_RESPONSE_ACCEPT,
			NULL);
}

SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Save File"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);
}

gint siril_dialog_run(SirilWidget *widgetdialog) {
	return gtk_dialog_run(GTK_DIALOG(GTK_FILE_CHOOSER(widgetdialog)));
}

void siril_widget_destroy(SirilWidget *widgetdialog) {
    while (!is_callback_called()) {
        g_main_context_iteration(NULL, FALSE);
    }
	gtk_widget_destroy(widgetdialog);
}
