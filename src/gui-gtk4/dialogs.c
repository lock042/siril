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

/* When a transient dialog is hidden, some compositors (notably Wayland) do
 * not automatically return keyboard focus to the parent window.  The main
 * window then stops receiving key events, so the app/win accelerators
 * (Ctrl+P, Ctrl+O, ...) silently stop working until the user clicks back
 * into it.  Re-presenting the control window restores it as the active
 * toplevel so the shortcuts keep working. */
static void reactivate_main_window() {
	GtkWidget *main_win = GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"));
	if (main_win && gtk_widget_get_visible(main_win))
		gtk_window_present(GTK_WINDOW(main_win));
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
		reactivate_main_window();
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
		siril_log_debug("### Opening imgproc dialog: %s\n", entry.identifier);
		processing_dialog_is_opened = TRUE;
	}
}

void siril_close_dialog(gchar *id) {
	gtk_widget_set_visible(get_widget_by_id(id), FALSE);
	reactivate_main_window();
	dialog_is_opened = FALSE;
	SirilDialogEntry entry = get_entry_by_id(id);
	if (entry.type == IMAGE_PROCESSING_DIALOG) {
		siril_log_debug("### Closing imgproc dialog: %s\n", entry.identifier);
		processing_dialog_is_opened = FALSE;
	}
}

void siril_close_preview_dialogs() {
	gboolean closed_any = FALSE;
	for (int i = 0; i < G_N_ELEMENTS(entries); i++) {
		GtkWidget *w = get_widget_by_index(i);
		if (gtk_widget_get_visible(w) && (entries[i].has_preview)) {
			entries[i].apply_function();
			gtk_widget_set_visible(w, FALSE);
			closed_any = TRUE;
		}
	}
	if (closed_any)
		reactivate_main_window();
}

// WARNING: do not use siril_widget_hide_on_delete() for IMAGE_PROCESSING_DIALOGs. These
// must call siril_close_dialog(builder_id) and therefore must have a custom handler
// as the GtkBuilder ID does not have a reverse lookup function.

gboolean siril_widget_hide_on_delete(GtkWidget *widget) {
    dialog_is_opened = FALSE;
    gtk_widget_set_visible(widget, FALSE);
    reactivate_main_window();
    return TRUE;
}

gboolean is_a_dialog_opened() {
	return dialog_is_opened;
}

gboolean is_an_image_processing_dialog_opened() {
	return processing_dialog_is_opened;
}

void mark_imgproc_dialog_closed() {
	siril_log_debug("### Closing imgproc dialog via custom hide_on_delete callback\n");
	dialog_is_opened = FALSE;
	processing_dialog_is_opened = FALSE;
}

/* ─── SirilFileChooser: GtkFileDialog wrapper ───────────────────────────
 *
 * GTK 4.10 deprecated GtkFileChooserDialog in favour of GtkFileDialog,
 * which is async-only.  SirilFileChooser preserves the synchronous
 * create → configure → run → read → destroy shape that all callers
 * use, by spinning a local GMainLoop in siril_fc_run() until the
 * underlying _open/_save/_select_folder/_open_multiple finish
 * callback fires.
 *
 * open_dialog.c stays on GtkFileChooserDialog (legacy path is local to
 * that file) because GtkFileDialog has no slot for a custom preview
 * widget.
 */

typedef enum {
	SFC_OPEN,            /* gtk_file_dialog_open */
	SFC_OPEN_MULTIPLE,   /* gtk_file_dialog_open_multiple */
	SFC_SAVE,            /* gtk_file_dialog_save */
	SFC_SELECT_FOLDER    /* gtk_file_dialog_select_folder */
} SirilFCMode;

struct _SirilFileChooser {
	SirilFCMode      mode;
	GtkFileDialog   *dialog;          /* owned */
	GtkWindow       *parent;          /* not owned */
	GListStore      *filters;         /* of GTK_TYPE_FILE_FILTER, lazy */
	GtkFileFilter   *default_filter;  /* in `filters`, not extra-ref'd */
	gchar           *current_name;
	gchar           *initial_folder;
	gchar           *initial_file;
	gboolean         select_multiple;
	/* result populated by finish() */
	GFile           *result_file;
	GListModel      *result_files;
	gint             response;
	GMainLoop       *loop;
};

static SirilFileChooser *siril_fc_new(GtkWindow *parent,
                                      SirilFCMode mode,
                                      const gchar *title) {
	SirilFileChooser *fc = g_new0(SirilFileChooser, 1);
	fc->mode    = mode;
	fc->parent  = parent;
	fc->dialog  = gtk_file_dialog_new();
	fc->response = GTK_RESPONSE_NONE;
	if (title)
		gtk_file_dialog_set_title(fc->dialog, title);
	return fc;
}

static SirilFCMode mode_from_action(GtkFileChooserAction action,
                                    gboolean prefer_multi) {
	switch (action) {
	case GTK_FILE_CHOOSER_ACTION_SAVE:
		return SFC_SAVE;
	case GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER:
		return SFC_SELECT_FOLDER;
	case GTK_FILE_CHOOSER_ACTION_OPEN:
	default:
		return prefer_multi ? SFC_OPEN_MULTIPLE : SFC_OPEN;
	}
}

SirilFileChooser *siril_fc_open(GtkWindow *parent, GtkFileChooserAction action) {
	const gchar *title = (action == GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER)
	                     ? _("Select Folder") : _("Open File");
	return siril_fc_new(parent, mode_from_action(action, FALSE), title);
}

SirilFileChooser *siril_fc_add(GtkWindow *parent, GtkFileChooserAction action) {
	return siril_fc_new(parent, mode_from_action(action, TRUE), _("Add Files"));
}

SirilFileChooser *siril_fc_save(GtkWindow *parent, GtkFileChooserAction action) {
	return siril_fc_new(parent, mode_from_action(action, FALSE), _("Save File"));
}

void siril_fc_set_current_name(SirilFileChooser *fc, const gchar *name) {
	if (!fc) return;
	g_free(fc->current_name);
	fc->current_name = name ? g_strdup(name) : NULL;
}

void siril_fc_set_current_folder_path(SirilFileChooser *fc, const gchar *path) {
	if (!fc) return;
	g_free(fc->initial_folder);
	fc->initial_folder = path ? g_strdup(path) : NULL;
}

void siril_fc_set_filename(SirilFileChooser *fc, const gchar *path) {
	if (!fc) return;
	g_free(fc->initial_file);
	fc->initial_file = path ? g_strdup(path) : NULL;
}

void siril_fc_set_select_multiple(SirilFileChooser *fc, gboolean multi) {
	if (!fc) return;
	fc->select_multiple = multi;
	if (multi && fc->mode == SFC_OPEN)
		fc->mode = SFC_OPEN_MULTIPLE;
	else if (!multi && fc->mode == SFC_OPEN_MULTIPLE)
		fc->mode = SFC_OPEN;
}

void siril_fc_add_filter(SirilFileChooser *fc, GtkFileFilter *filter,
                         gboolean set_default) {
	if (!fc || !filter) return;
	if (!fc->filters)
		fc->filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
	g_list_store_append(fc->filters, filter);
	if (set_default || !fc->default_filter)
		fc->default_filter = filter;
}

void siril_fc_add_filter_pattern(SirilFileChooser *fc, const gchar *title,
                                 const gchar *pattern, gboolean set_default) {
	if (!fc || !title || !pattern) return;
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	gchar **patterns = g_strsplit(pattern, ";", -1);
	for (gint i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	g_strfreev(patterns);
	siril_fc_add_filter(fc, f, set_default);
	g_object_unref(f);
}

GtkFileFilter *siril_fc_get_filter(SirilFileChooser *fc) {
	return fc ? fc->default_filter : NULL;
}

static void siril_fc_apply_pending(SirilFileChooser *fc) {
	if (fc->filters) {
		gtk_file_dialog_set_filters(fc->dialog, G_LIST_MODEL(fc->filters));
		if (fc->default_filter)
			gtk_file_dialog_set_default_filter(fc->dialog, fc->default_filter);
	}
	if (fc->initial_folder) {
		GFile *gf = g_file_new_for_path(fc->initial_folder);
		gtk_file_dialog_set_initial_folder(fc->dialog, gf);
		g_object_unref(gf);
	}
	if (fc->initial_file) {
		GFile *gf = g_file_new_for_path(fc->initial_file);
		gtk_file_dialog_set_initial_file(fc->dialog, gf);
		g_object_unref(gf);
	}
	if (fc->current_name && fc->mode == SFC_SAVE)
		gtk_file_dialog_set_initial_name(fc->dialog, fc->current_name);
}

static void siril_fc_finish_cb(GObject *src, GAsyncResult *res, gpointer user_data) {
	SirilFileChooser *fc = user_data;
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *err = NULL;

	switch (fc->mode) {
	case SFC_OPEN:
		fc->result_file = gtk_file_dialog_open_finish(fd, res, &err);
		break;
	case SFC_OPEN_MULTIPLE:
		fc->result_files = gtk_file_dialog_open_multiple_finish(fd, res, &err);
		break;
	case SFC_SAVE:
		fc->result_file = gtk_file_dialog_save_finish(fd, res, &err);
		break;
	case SFC_SELECT_FOLDER:
		fc->result_file = gtk_file_dialog_select_folder_finish(fd, res, &err);
		break;
	}

	if (fc->result_file || fc->result_files)
		fc->response = GTK_RESPONSE_ACCEPT;
	else if (err && g_error_matches(err, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_DISMISSED))
		fc->response = GTK_RESPONSE_CANCEL;
	else
		fc->response = GTK_RESPONSE_NONE;
	g_clear_error(&err);

	if (fc->loop)
		g_main_loop_quit(fc->loop);
}

gint siril_fc_run(SirilFileChooser *fc) {
	if (!fc) return GTK_RESPONSE_NONE;
	/* Reset any leftover state from a previous run on this chooser. */
	g_clear_object(&fc->result_file);
	g_clear_object(&fc->result_files);
	fc->response = GTK_RESPONSE_NONE;
	siril_fc_apply_pending(fc);
	fc->loop = g_main_loop_new(NULL, FALSE);
	switch (fc->mode) {
	case SFC_OPEN:
		gtk_file_dialog_open(fc->dialog, fc->parent, NULL,
		                     siril_fc_finish_cb, fc);
		break;
	case SFC_OPEN_MULTIPLE:
		gtk_file_dialog_open_multiple(fc->dialog, fc->parent, NULL,
		                              siril_fc_finish_cb, fc);
		break;
	case SFC_SAVE:
		gtk_file_dialog_save(fc->dialog, fc->parent, NULL,
		                     siril_fc_finish_cb, fc);
		break;
	case SFC_SELECT_FOLDER:
		gtk_file_dialog_select_folder(fc->dialog, fc->parent, NULL,
		                              siril_fc_finish_cb, fc);
		break;
	}
	g_main_loop_run(fc->loop);
	g_main_loop_unref(fc->loop);
	fc->loop = NULL;
	return fc->response;
}

gchar *siril_fc_get_filename(SirilFileChooser *fc) {
	if (!fc) return NULL;
	if (fc->result_file)
		return g_file_get_path(fc->result_file);
	if (fc->mode == SFC_SAVE && fc->current_name)
		return g_strdup(fc->current_name);
	return NULL;
}

GSList *siril_fc_get_filenames(SirilFileChooser *fc) {
	GSList *out = NULL;
	if (!fc) return NULL;
	if (fc->result_files) {
		guint n = g_list_model_get_n_items(fc->result_files);
		for (guint i = 0; i < n; ++i) {
			GFile *gf = G_FILE(g_list_model_get_item(fc->result_files, i));
			if (gf) {
				gchar *p = g_file_get_path(gf);
				if (p) out = g_slist_append(out, p);
				g_object_unref(gf);
			}
		}
	} else if (fc->result_file) {
		gchar *p = g_file_get_path(fc->result_file);
		if (p) out = g_slist_append(out, p);
	}
	return out;
}

void siril_fc_destroy(SirilFileChooser *fc) {
	if (!fc) return;
	g_clear_object(&fc->result_file);
	g_clear_object(&fc->result_files);
	g_clear_object(&fc->filters);
	g_clear_object(&fc->dialog);
	g_free(fc->current_name);
	g_free(fc->initial_folder);
	g_free(fc->initial_file);
	g_free(fc);
}
