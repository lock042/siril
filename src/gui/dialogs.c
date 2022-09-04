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
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "algos/background_extraction.h"
#include "filters/asinh.h"
#include "filters/deconv.h"
#include "filters/clahe.h"
#include "filters/saturation.h"
#include "filters/wavelets.h"

#include "gui/dialogs.h"
#include "nina_light_curve.h"

static const SirilDialogEntry entries[] =
{
	{"asinh_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_asinh_cancel},
	{"background_extraction_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_background_cancel},
	{"canon_fixbanding_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"CLAHE_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_clahe_cancel},
	{"composition_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"color_calibration", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"cosmetic_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"crop_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"deconvolution_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_deconv_cancel},
	{"dialog_FFT", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"dialog_star_remix", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"edge_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"extract_channel_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"extract_wavelets_layers_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"file_information", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"histogram_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_histo_cancel},
	{"ImagePlateSolver_Dial", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"linearmatch_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"Median_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"nina_light_curve", get_nina_lc_dialog, OTHER_DIALOG, FALSE, NULL},
	{"pixel_math_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"resample_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"rgradient_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"rotation_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"satu_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_satu_cancel},
	{"SCNR_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"search_objects", NULL, SEARCH_ENTRY_DIALOG, FALSE, NULL},
	{"settings_window", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"seqlist_dialog", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"split_cfa_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"starnet_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
	{"stars_list_window", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"StatWindow", NULL, INFORMATION_DIALOG, FALSE, NULL},
	{"synthstar_dialog", NULL, IMAGE_PROCESSING_DIALOG, FALSE, NULL},
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
}

void siril_close_dialog(gchar *id) {
	gtk_widget_hide(get_widget_by_id(id));
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
