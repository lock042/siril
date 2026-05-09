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
#include "core/processing_thread.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/histogram.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/curves.h"
#include "algos/background_extraction.h"
#include "gui-gtk4/asinh.h"
#include "filters/epf.h"
#include "gui-gtk4/clahe.h"
#include "gui-gtk4/median.h"
#include "gui-gtk4/saturation.h"
#include "filters/unpurple.h"
#include "gui-gtk4/wavelets.h"

#include "gui-gtk4/newdeconv.h"
#include "filters/nlbayes/call_nlbayes.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/dialog_preview.h"
#include "nina_light_curve.h"
#include "compstars.h"
#include "catmag.h"

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
	{"catmag", get_catmag_dialog, OTHER_DIALOG, FALSE, NULL},
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
	{"wavelets_dialog", NULL, IMAGE_PROCESSING_DIALOG, TRUE, apply_wavelets_cancel},
	{"mask_from_color_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_from_image_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_from_stars_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_blur_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_feather_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_scale_dialog", NULL, OTHER_DIALOG, FALSE, NULL},
	{"mask_thresholds_dialog", NULL, OTHER_DIALOG, FALSE, NULL}
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

int number_of_dialogs() {
	return G_N_ELEMENTS(entries);
}

static GtkWidget *get_widget_by_id(gchar *id) {
	SirilDialogEntry entry = get_entry_by_id(id);
	if (entry.get_window)
		return entry.get_window();
	return GTK_WIDGET(gtk_builder_get_object(gui.builder, entry.identifier));
}

static GtkWidget *get_widget_by_index(int index) {
	SirilDialogEntry entry = entries[index];
	if (entry.get_window)
		return entry.get_window();
	return GTK_WIDGET(gtk_builder_get_object(gui.builder, entry.identifier));
}

/* Phase 14: GTK4 routes key input through GtkEventControllerKey instead
 * of the legacy "key-press-event" signal.  The escape-to-close behaviour
 * for the search-entry dialog is hooked via gtk_event_controller_key,
 * with the controller attached to the window once on first open. */
static gboolean check_escape(GtkEventControllerKey *ctrl, guint keyval,
		guint keycode, GdkModifierType state, gpointer data) {
	(void)ctrl; (void)keycode; (void)state;
	if (keyval == GDK_KEY_Escape) {
		GtkWidget *widget = GTK_WIDGET(data);
		gtk_widget_set_visible(widget, FALSE);
		return TRUE;
	}
	return FALSE;
}


void siril_open_dialog(gchar *id) {
	(void)0; /* Phase 14: gtk_window_get_position / gtk_window_move are gone
	          * in GTK4; window placement is now compositor-driven.  We no
	          * longer record/restore positions across reopens. */
	SirilDialogEntry entry = get_entry_by_id(id);
	GtkWindow *win = GTK_WINDOW(get_widget_by_id(id));

	// We cannot open the dialog if a python script claims the thread, to prevent conflict over
	// updating gfit
	if (entry.type == IMAGE_PROCESSING_DIALOG && processing_is_reserved_for_python()) {
		queue_error_message_dialog(_("Error"), _("Error: cannot open an image processing dialog while a Python script claims the processing thread. "
									"Wait for the Python script to finish processing first."));
		return;
	}

	if (entry.type != INFORMATION_DIALOG) {
		for (int i = 0; i < G_N_ELEMENTS(entries); i++) {
			GtkWidget *w = get_widget_by_index(i);
			if (gtk_widget_get_visible(w) && entries[i].type != INFORMATION_DIALOG) {
				if (entries[i].has_preview)
					entries[i].apply_function();
				gtk_widget_set_visible(w, FALSE);
				break;
			}
		}
	}
	if (entry.type == SEARCH_ENTRY_DIALOG) {
		/* Attach a key controller once per window; subsequent opens reuse it. */
		if (!g_object_get_data(G_OBJECT(win), "siril-escape-ctrl")) {
			GtkEventController *ctrl = gtk_event_controller_key_new();
			g_signal_connect(ctrl, "key-pressed", G_CALLBACK(check_escape), win);
			gtk_widget_add_controller(GTK_WIDGET(win), ctrl);
			g_object_set_data(G_OBJECT(win), "siril-escape-ctrl", ctrl);
		}
	}

	/* GTK4: GTK_WIN_POS_CENTER and GDK_WINDOW_TYPE_HINT_DIALOG are gone.
	 * The window manager centres modal/transient windows automatically
	 * when set_transient_for is set. */
	gtk_window_set_transient_for(win, GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	gtk_window_present(win);
	dialog_is_opened = TRUE;
	if (entry.type == IMAGE_PROCESSING_DIALOG) {
		siril_debug_print("### Opening imgproc dialog: %s\n", entry.identifier);
		processing_dialog_is_opened = TRUE;
	}
}

void siril_close_dialog(gchar *id) {
	gtk_widget_set_visible(get_widget_by_id(id), FALSE);
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
			gtk_widget_set_visible(w, FALSE);
		}
	}
}

// WARNING: do not use siril_widget_hide_on_delete() for IMAGE_PROCESSING_DIALOGs. These
// must call siril_close_dialog(builder_id) and therefore must have a custom handler
// as the GtkBuilder ID does not have a reverse lookup function.

gboolean siril_widget_hide_on_delete(GtkWidget *widget) {
    dialog_is_opened = FALSE;
    gtk_widget_set_visible(widget, FALSE);
    return TRUE;
}

gboolean is_a_dialog_opened() {
	return dialog_is_opened;
}

gboolean is_an_image_processing_dialog_opened() {
	return processing_dialog_is_opened;
}

void mark_imgproc_dialog_closed() {
	siril_debug_print("### Closing imgproc dialog via custom hide_on_delete callback\n");
	dialog_is_opened = FALSE;
	processing_dialog_is_opened = FALSE;
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
	G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	if (set_default)
		gtk_file_chooser_set_filter(file_chooser, f);
	G_GNUC_END_IGNORE_DEPRECATIONS
}

/* Phase 14G.2: dialogs.c hub — full GtkFileDialog migration is a 13-file
 * refactor (every caller does GTK_FILE_CHOOSER(widgetdialog) and calls
 * methods on it after construction).  For now we keep the public API
 * shape (SirilWidget* + sync siril_dialog_run + siril_widget_destroy)
 * by continuing to use GtkFileChooserDialog + gtk_dialog_run, both of
 * which are deprecated-but-functional in GTK 4.10.  The deprecation
 * warnings are silenced here instead of leaking through ~13 callers.
 *
 * The standalone file-chooser dialogs in masks_gui.c, python_gui.c and
 * pixelmath.c have already been migrated to GtkFileDialog in this same
 * Phase 14G.2 effort. */
SirilWidget *siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action) {
	gchar *title;
	SirilWidget *w;
	if (action == GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER) {
		title = g_strdup(_("Select Folder"));
	} else {
		title = g_strdup(_("Open File"));
	G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
	}
	w = gtk_file_chooser_dialog_new(title, parent, action, _("_Cancel"),
			GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);
	G_GNUC_END_IGNORE_DEPRECATIONS
	g_free(title);
	return w;
}

G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
SirilWidget *siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Add Files"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Add"), GTK_RESPONSE_ACCEPT,
			NULL);
G_GNUC_END_IGNORE_DEPRECATIONS
}

G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Save File"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);
G_GNUC_END_IGNORE_DEPRECATIONS
}

struct siril_dialog_ctx { gint response; GMainLoop *loop; };

static void siril_dialog_response_cb(GtkDialog *d, gint r, gpointer p) {
	(void) d;
	struct siril_dialog_ctx *c = p;
	c->response = r;
	g_main_loop_quit(c->loop);
}

gint siril_dialog_run(SirilWidget *widgetdialog) {
	/* GTK4 removed gtk_dialog_run.  Spin a local GMainLoop until the dialog
	 * emits "response", then return the response code.  Callers continue to
	 * use this synchronous shape; per-call-site async migration is tracked
	 * elsewhere in the plan. */
	struct siril_dialog_ctx ctx = { GTK_RESPONSE_NONE, g_main_loop_new(NULL, FALSE) };
	gulong id = g_signal_connect(widgetdialog, "response",
		G_CALLBACK(siril_dialog_response_cb), &ctx);
	gtk_window_present(GTK_WINDOW(widgetdialog));
	g_main_loop_run(ctx.loop);
	g_signal_handler_disconnect(widgetdialog, id);
	g_main_loop_unref(ctx.loop);
	return ctx.response;
}

void siril_widget_destroy(SirilWidget *widgetdialog) {
    while (!is_callback_called()) {
        g_main_context_iteration(NULL, FALSE);
    }
	gtk_window_destroy(GTK_WINDOW(widgetdialog));
}
