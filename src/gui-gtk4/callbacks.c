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

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>

#include "core/siril.h"
#include "core/processing_thread.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "core/masks.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/siril_language.h"
#include "core/siril_networking.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "compositing/compositing.h"
#include "gui-gtk4/cut.h"
#include "gui-gtk4/open_dialog.h"
#include "gui-gtk4/icc_profile.h"
#include "gui-gtk4/keywords_tree.h"
#include "gui-gtk4/registration.h"
#include "gui-gtk4/photometric_cc.h"
#include "gui-gtk4/stacking.h"
#include "gui-gtk4/python_gui.h"
#include "gui-gtk4/undo_gui.h"
/* Forward declaration: defined in gui/sequence_export.c */
void update_export_crop_label();
#include "algos/siril_wcs.h"
#include "io/annotation_catalogues.h"
#include "io/films.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/siril_git.h"
#include "io/siril_pythonmodule.h"
#include "annotations_pref.h"
#include "image_display.h"
#include "image_interactions.h"
#include "single_image.h"
#include "sequence_list.h"
#include "callbacks.h"

#include "algos/astrometry_solver.h"
#include "utils.h"
#include "plot.h"
#include "message_dialog.h"
#include "PSF_list.h"
#include "histogram.h"
#include "remixer.h"
#include "script_menu.h"
#include "progress_and_log.h"
#include "dialogs.h"
#include "siril_intro.h"
#include "siril_preview.h"
#include "siril-window.h"
#include "registration_preview.h"
#include "io/healpix/fluxcache_cat.h"

static GList *roi_callbacks = NULL;
static gchar *display_item_name[] = { "linear_item", "log_item", "square_root_item", "squared_item", "asinh_item", "auto_item", "histo_item", "softproof_item"};

void set_viewer_mode_widgets_sensitive(gboolean sensitive) {
	GtkWidget *scalemax = lookup_widget("scalemax");
	GtkWidget *scalemin = lookup_widget("scalemin");
	GtkWidget *entrymin = lookup_widget("min_entry");
	GtkWidget *entrymax = lookup_widget("max_entry");
	GtkWidget *minmax = lookup_widget("radiobutton_minmax");
	GtkWidget *hilo = lookup_widget("radiobutton_hilo");
	GtkWidget *user = lookup_widget("radiobutton_user");

	gtk_widget_set_sensitive(scalemax, sensitive);
	gtk_widget_set_sensitive(scalemin, sensitive);
	gtk_widget_set_sensitive(entrymin, sensitive);
	gtk_widget_set_sensitive(entrymax, sensitive);
	gtk_widget_set_sensitive(minmax, sensitive);
	gtk_widget_set_sensitive(hilo, sensitive);
	gtk_widget_set_sensitive(user, sensitive);
}

static void update_theme_button(const gchar *button_name, const gchar *path) {
	gchar *image;
	GtkWidget *w_image;

	image = g_strdup_printf("/org/siril/ui/pixmaps/%s", path);
	w_image = gtk_image_new_from_resource(image);
	gtk_widget_set_visible(w_image, TRUE);
	gtk_button_set_child(GTK_BUTTON(lookup_widget(button_name)), w_image);
	gtk_widget_set_visible(lookup_widget(button_name), TRUE);

	g_free(image);
}

/* Phase 13: GTK4 GdkClipboard text reads are asynchronous.  The original
 * GTK3 code called gtk_clipboard_wait_for_text() synchronously inside the
 * "owner-change" handler.  We split that into the signal callback (which
 * just kicks off the async read) and a continuation that does the actual
 * label-colour update once the text is available. */

static void clipboard_text_ready_cb(GObject *source, GAsyncResult *result, gpointer user_data) {
	(void)user_data;
	GdkClipboard *clipboard = GDK_CLIPBOARD(source);
	GError *error = NULL;
	char *clipboard_content = gdk_clipboard_read_text_finish(clipboard, result, &error);
	if (error) {
		g_clear_error(&error);
	}
	const char *content = clipboard_content ? clipboard_content : "NoRules";

	GtkLabel *label_name_of_seq = GTK_LABEL(lookup_widget("label_name_of_seq"));
	const char *format_green = "<span foreground=\"green\">%s</span>";
	const char *format_white = "<span foreground=\"white\">%s</span>";
	char *markup;

	if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		if (strcmp(filename, content) == 0) {
			markup = g_markup_printf_escaped(format_green, _("Image: "));
		} else {
			markup = g_markup_printf_escaped(format_white, _("Image: "));
		}
		gtk_label_set_markup(label_name_of_seq, markup);
		g_free(markup);
		g_free(filename);
	}

	if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);
		if (strcmp(seq_basename, content) == 0) {
			markup = g_markup_printf_escaped(format_green, _("Sequence: "));
		} else {
			markup = g_markup_printf_escaped(format_white, _("Sequence: "));
		}
		gtk_label_set_markup(label_name_of_seq, markup);
		g_free(markup);
		g_free(seq_basename);
	}

	g_free(clipboard_content);
}

void handle_owner_change(GdkClipboard *clipboard, gpointer data) {
	(void)data;
	/* Only fire the async text read when the clipboard actually has
	 * text on offer.  Without this guard, e.g. KDE Spectacle's
	 * "Copy" — which puts an image/png on the clipboard — triggers
	 * GTK4's text-read code path on data that can't be coerced to
	 * a string, and that path has crashed reliably for us (the
	 * fault surfaces inside libgtk / libgio's async dispatch with
	 * no Siril frames on the stack).  Checking formats up-front
	 * is also cheaper than firing a doomed async read. */
	GdkContentFormats *formats = gdk_clipboard_get_formats(clipboard);
	if (!formats) return;
	if (!gdk_content_formats_contain_gtype(formats, G_TYPE_STRING)) {
		/* Clipboard now holds something non-textual (image, files,
		 * custom blob).  The label-colour cue isn't applicable — just
		 * skip this update. */
		return;
	}
	gdk_clipboard_read_text_async(clipboard, NULL, clipboard_text_ready_cb, NULL);
}

gboolean launch_clipboard_survey(gpointer user_data) {
	(void)user_data;
	if (com.script)
		return FALSE;

	/* launch_clipboard_survey is wired to several entry points (startup,
	 * open action, click-on-seq-label, gui_iface impl) — without this
	 * guard each one stacks an additional handler on the clipboard,
	 * leaking and firing the read N times per owner change. */
	static gboolean connected = FALSE;
	if (connected) return FALSE;

	GdkClipboard *clipboard = gdk_display_get_clipboard(gdk_display_get_default());
	if (!clipboard) return FALSE;
	g_signal_connect(clipboard, "changed", G_CALLBACK(handle_owner_change), NULL);
	connected = TRUE;
	return FALSE;
}

void on_press_seq_field() {
	GdkClipboard *clipboard = gdk_display_get_clipboard(gdk_display_get_default());

	if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		gdk_clipboard_set_text(clipboard, filename);
		g_free(filename);
	} else if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);
		gdk_clipboard_set_text(clipboard, seq_basename);
		g_free(seq_basename);
	}
}

static void update_icons_to_theme(gboolean is_dark) {
	siril_log_debug("Loading %s theme...\n", is_dark ? "dark" : "light");
	if (is_dark) {
		update_theme_button("annotate_button", "astrometry_dark.svg");
		update_theme_button("wcs_grid_button", "wcs-grid_dark.svg");
		update_theme_button("photometry_button", "photometry_dark.svg");
		update_theme_button("cut_button", "cut_dark.svg");

		update_theme_button("rotate90_anticlock_button", "rotate-acw_dark.svg");
		update_theme_button("rotate90_clock_button", "rotate-cw_dark.svg");
		update_theme_button("mirrorx_button", "mirrorx_dark.svg");
		update_theme_button("mirrory_button", "mirrory_dark.svg");

		update_theme_button("process_starfinder_button", "starfinder_dark.svg");
		update_theme_button("sum_button", "sum_dark.svg");
		update_theme_button("export_button", "export_dark.svg");

		update_theme_button("histoToolAutoStretch", "mtf_dark.svg");

		update_theme_button("histo_button", "histo_dark.svg");
	} else {
		update_theme_button("annotate_button", "astrometry.svg");
		update_theme_button("wcs_grid_button", "wcs-grid.svg");
		update_theme_button("photometry_button", "photometry.svg");
		update_theme_button("cut_button", "cut.svg");

		update_theme_button("rotate90_anticlock_button", "rotate-acw.svg");
		update_theme_button("rotate90_clock_button", "rotate-cw.svg");
		update_theme_button("mirrorx_button", "mirrorx.svg");
		update_theme_button("mirrory_button", "mirrory.svg");

		update_theme_button("process_starfinder_button", "starfinder.svg");
		update_theme_button("sum_button", "sum.svg");
		update_theme_button("export_button", "export.svg");

		update_theme_button("histoToolAutoStretch", "mtf.svg");

		update_theme_button("histo_button", "histo.svg");
	}
}

void siril_set_theme(int active) {
	GtkSettings *settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", active == 0, NULL);
	update_icons_to_theme(active == 0);
	update_icons_sequence_list(active == 0);
}

void on_combo_theme_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	int active = gtk_drop_down_get_selected(box);

	siril_set_theme(active);
}

int populate_roi() {
	if (com.python_command)
		return 1;
	if (gui.roi.selection.w == 0 || gui.roi.selection.h == 0)
		return 1;
	int retval = 0;
	size_t npixels_roi = gui.roi.selection.w * gui.roi.selection.h;
	size_t npixels_gfit = gfit->rx * gfit->ry;
	size_t nchans = gfit->naxes[2];
	g_assert(nchans == 1 || nchans == 3);
	gboolean rgb = (nchans == 3);

	clearfits(&gui.roi.fit);
	copyfits(gfit, &gui.roi.fit, CP_FORMAT, -1);

	gui.roi.fit.rx = gui.roi.fit.naxes[0] = gui.roi.selection.w;
	gui.roi.fit.ry = gui.roi.fit.naxes[1] = gui.roi.selection.h;
	gui.roi.fit.naxes[2] = nchans;
	gui.roi.fit.naxis = (nchans == 1 ? 2 : 3);

	/* ---------------- IMAGE ROI COPY ---------------- */
	if (gui.roi.fit.type == DATA_FLOAT) {
		gui.roi.fit.fdata = malloc(npixels_roi * nchans * sizeof(float));
		if (!gui.roi.fit.fdata)
			retval = 1;

		gui.roi.fit.fpdata[0] = gui.roi.fit.fdata;
		gui.roi.fit.fpdata[1] = rgb ? gui.roi.fit.fdata + npixels_roi : gui.roi.fit.fdata;
		gui.roi.fit.fpdata[2] = rgb ? gui.roi.fit.fdata + 2 * npixels_roi : gui.roi.fit.fdata;

		for (uint32_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h; y++) {
				float *srcindex =
					gfit->fdata +
					(npixels_gfit * c) +
					((gfit->ry - y - (gui.roi.selection.y + 1)) * gfit->rx) +
					gui.roi.selection.x;

				float *destindex =
					gui.roi.fit.fdata +
					(npixels_roi * c) +
					(gui.roi.fit.rx * y);

				memcpy(destindex, srcindex, gui.roi.selection.w * sizeof(float));
			}
		}
	} else {
		gui.roi.fit.data = malloc(npixels_roi * nchans * sizeof(WORD));
		if (!gui.roi.fit.data)
			retval = 1;

		gui.roi.fit.pdata[0] = gui.roi.fit.data;
		gui.roi.fit.pdata[1] = rgb ? gui.roi.fit.data + npixels_roi : gui.roi.fit.data;
		gui.roi.fit.pdata[2] = rgb ? gui.roi.fit.data + 2 * npixels_roi : gui.roi.fit.data;

		for (uint32_t c = 0; c < nchans; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h; y++) {
				WORD *srcindex =
					gfit->data +
					(npixels_gfit * c) +
					((gfit->ry - y - (gui.roi.selection.y + 1)) * gfit->rx) +
					gui.roi.selection.x;

				WORD *destindex =
					gui.roi.fit.data +
					(npixels_roi * c) +
					(y * gui.roi.fit.rx);

				memcpy(destindex, srcindex, gui.roi.selection.w * sizeof(WORD));
			}
		}
	}

	/* ---------------- MASK ROI COPY ---------------- */
	gui.roi.fit.mask_active = FALSE;
	gui.roi.fit.mask = NULL;

	if (gfit->mask && gfit->mask->data) {
		size_t n = gui.roi.selection.w * gui.roi.selection.h;
		size_t src_w = gfit->rx;
		size_t src_h = gfit->ry;

		gui.roi.fit.mask = malloc(sizeof(mask_t));
		if (!gui.roi.fit.mask) {
			retval = 1;
			goto finalize;
		}

		gui.roi.fit.mask->bitpix = gfit->mask->bitpix;

		size_t elem_size = 0;
		switch (gfit->mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			free(gui.roi.fit.mask);
			gui.roi.fit.mask = NULL;
			goto finalize;
		}

		gui.roi.fit.mask->data = malloc(n * elem_size);
		if (!gui.roi.fit.mask->data) {
			free(gui.roi.fit.mask);
			gui.roi.fit.mask = NULL;
			retval = 1;
			goto finalize;
		}

		/* Copy mask ROI with same flip logic as image */
		for (uint32_t y = 0; y < gui.roi.selection.h; y++) {
			size_t src_y = src_h - y - (gui.roi.selection.y + 1);
			size_t src_x = gui.roi.selection.x;

			void *src = (uint8_t*)gfit->mask->data +
						(src_y * src_w + src_x) * elem_size;

			void *dst = (uint8_t*)gui.roi.fit.mask->data +
						(y * gui.roi.selection.w) * elem_size;

			memcpy(dst, src, gui.roi.selection.w * elem_size);
		}

		gui.roi.fit.mask_active = gfit->mask_active;
	}

finalize:
	backup_roi();
	gui.roi.active = TRUE;
	return retval;
}

static void call_roi_callbacks() {
	if (com.python_command)
		return;
	GList *current = roi_callbacks;
	while (current != NULL) {
		ROICallback func = (ROICallback)(current->data);
		func();
		current = current->next;
	}
}

void add_roi_callback(ROICallback func) {
	roi_callbacks = g_list_append(roi_callbacks, func);
}

void remove_roi_callback(ROICallback func) {
	roi_callbacks = g_list_remove_all(roi_callbacks, func);
}

static void roi_info_message_if_needed() {
	/* Phase 14G.1: GtkMessageDialog removed; route through the central
	 * "confirm + remember" wrapper which uses a custom GtkWindow under
	 * the hood (GtkAlertDialog cannot host a checkbox). */
	gboolean dont_show_again = FALSE;
	siril_confirm_dialog_and_remember(
		_("Region of Interest supported"),
		_("This image operation supports ROI processing. This is "
		  "indicated by the green border around the ROI. In image "
		  "operation dialogs that do not support ROI processing, and "
		  "when no dialog is open, the ROI will have a red outline.\n\n"
		  "While a ROI is active, processing will only preview the ROI. "
		  "When the Apply button is clicked, the operation will apply "
		  "to the whole image."),
		_("_OK"),
		&dont_show_again);
	com.pref.gui.enable_roi_warning = !dont_show_again;
}

void roi_supported(gboolean state) {
	gui.roi.operation_supports_roi = state;
	// Provide a one-time notification about ROI support
	if (state && gui.roi.active && com.pref.gui.enable_roi_warning)
		roi_info_message_if_needed();
	redraw(REDRAW_OVERLAY);
}

static GMutex roi_mutex;

void lock_roi_mutex() {
	g_mutex_lock(&roi_mutex);
}

void unlock_roi_mutex() {
	g_mutex_unlock(&roi_mutex);
}

gpointer on_set_roi() {
	if (com.selection.w == 0 && com.selection.h == 0)
		return GINT_TO_POINTER(0);
	rectangle sel = com.selection;
	g_mutex_lock(&roi_mutex); // Wait until any thread previews are finished
	if (com.python_command) {
		g_mutex_unlock(&roi_mutex);
		return GINT_TO_POINTER(0);
	}
	cancel_pending_update();
	if (gui.roi.operation_supports_roi && com.pref.gui.enable_roi_warning)
		roi_info_message_if_needed();
	// Ensure any pending ROI changes are overwritten by the backup
	// Must copy the whole backup to gfit and gui.roi.fit to account
	// for switching between full image and ROI
	memcpy(&gui.roi.selection, &sel, sizeof(rectangle));
	gui.roi.active = TRUE;
	if (gui.roi.selection.w > 0 && gui.roi.selection.h > 0 && is_preview_active())
		copy_backup_to_gfit();
	populate_roi();
	backup_roi();
	// Call any callbacks that need calling
	call_roi_callbacks();
	g_mutex_unlock(&roi_mutex);
	return GINT_TO_POINTER(0);
}

gpointer on_clear_roi() {
	if (gui.roi.active) {
		g_mutex_lock(&roi_mutex); // Wait until any thread previews are finished
		cancel_pending_update();
		copy_backup_to_gfit();
		clearfits(&gui.roi.fit);
		memset(&gui.roi, 0, sizeof(roi_t));
		redraw(REDRAW_OVERLAY);
		// Call any callbacks that need calling
		call_roi_callbacks();
		g_mutex_unlock(&roi_mutex);
	}
	return GINT_TO_POINTER(0);
}

static void initialize_theme_GUI() {
	GtkDropDown *box;

	box = GTK_DROP_DOWN(lookup_widget("combo_theme"));

	g_signal_handlers_block_by_func(box, on_combo_theme_changed, NULL);
	gtk_drop_down_set_selected(box, com.pref.gui.combo_theme);
	g_signal_handlers_unblock_by_func(box, on_combo_theme_changed, NULL);
	update_icons_to_theme(com.pref.gui.combo_theme == 0);
	update_icons_sequence_list(com.pref.gui.combo_theme == 0);
}

void load_prefered_theme(gint theme) {
	GtkSettings *settings;

	settings = gtk_settings_get_default();

	g_object_set(settings, "gtk-application-prefer-dark-theme", com.pref.gui.combo_theme == 0, NULL);
}

void set_sliders_value_to_gfit() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;

	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemax"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemin"));// scalemin
	}

	gfit->keywords.hi = gtk_adjustment_get_value(adj1);
	gfit->keywords.lo = gtk_adjustment_get_value(adj2);
}

/* Sets maximum value for contrast scales. Minimum is always 0.
 * Should be done on first image load of a sequence and when single images are loaded.
 * Max value is taken from gfit->maxi, recomputed if not present.
 */
void set_cutoff_sliders_max_values() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;
	gdouble max_val;
	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemax"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemin"));// scalemin
	}
	/* set max value for range according to number of bits of original image
	 */

	max_val = (gfit->type == DATA_FLOAT ? USHRT_MAX_DOUBLE : (gfit->orig_bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE));
	siril_log_debug(_("Setting MAX value for cutoff sliders adjustments (%f)\n"), max_val);
	gtk_adjustment_set_lower(adj1, 0.0);
	gtk_adjustment_set_lower(adj2, 0.0);
	gtk_adjustment_set_upper(adj1, max_val);
	gtk_adjustment_set_upper(adj2, max_val);
}

/* Sets the value scalemin and scalemax sliders so that they match the hi and
 * lo values stored for the current viewport.
 * It also sets the 'cut' checkbox status.
 * The function should be called on every tab change and file opening */
void set_cutoff_sliders_values() {
	gchar buffer[10];
	static GtkAdjustment *adjmin = NULL, *adjmax = NULL;
	static GtkEntry *maxentry = NULL, *minentry = NULL;
	static GtkCheckButton *cutmax = NULL;
	if (adjmin == NULL) {
		adjmax = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemax")); // scalemax
		adjmin = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemin")); // scalemin
		maxentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "max_entry"));
		minentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "min_entry"));
		cutmax = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkcut_max"));
	}

	/* Snapshot gui.hi/gui.lo under com.mutex: notify_gfit_data_modified() may
	 * write them from the worker thread concurrently. */
	g_mutex_lock(&com.mutex);
	WORD hi = gui.hi;
	WORD lo = gui.lo;
	g_mutex_unlock(&com.mutex);

	siril_log_debug(_("Setting ranges scalemin=%d, scalemax=%d\n"), lo, hi);
	gtk_adjustment_set_value(adjmin, (gdouble)lo);
	gtk_adjustment_set_value(adjmax, (gdouble)hi);
	g_snprintf(buffer, 6, "%u", hi);
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(maxentry), buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_snprintf(buffer, 6, "%u", lo);
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(minentry), buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	siril_toggle_set_active(GTK_WIDGET(cutmax), gui.cut_over);
}

gboolean set_cutoff_sliders_values_idle(gpointer p) {
	set_cutoff_sliders_values();
	return FALSE;
}

/* Remap gfit with the new display mode during a sequence operation.
 * Uses trylock so it never blocks the GTK thread; if the lock fails the
 * completion idle will repaint when the job finishes. */
static gboolean try_remap_for_mode_change_idle(gpointer p) {
	if (g_rw_lock_reader_trylock(&gfit->rwlock)) {
		remap_all();
		g_rw_lock_reader_unlock(&gfit->rwlock);
		redraw(REMAP_ALL);
	}
	return FALSE;
}

void on_display_item_toggled(GtkCheckButton *checkmenuitem, gpointer user_data) {
	if (!siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(checkmenuitem)))) return;

	static GtkLabel *label_display_menu = NULL;
	if (!label_display_menu) {
		label_display_menu = GTK_LABEL(lookup_widget("display_button_name"));
	}

	gui.rendering_mode = get_display_mode_from_menu();
	siril_log_debug("Display mode %d\n", gui.rendering_mode);
	gboolean override_label = FALSE;

	if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type != DATA_FLOAT)
		siril_log_message(_("Current image is not 32 bit. Standard 16 bit AutoStretch will be used.\n"));
	if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap) {
		allocate_hd_remap_indices();
		siril_log_message(_("The AutoStretch display mode will use a %d bit LUT\n"), (int) log2(gui.hd_remap_max));
	} else {
		hd_remap_indices_cleanup();
		if (gui.rendering_mode == STF_DISPLAY)
			siril_log_message(_("The AutoStretch display mode will use a 16 bit LUT\n"));
	}
	if (!override_label)
		gtk_label_set_text(label_display_menu, gtk_check_button_get_label(GTK_CHECK_BUTTON(checkmenuitem)));

	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	siril_window_autostretch_actions(app_win, gui.rendering_mode == STF_DISPLAY && gfit->naxes[2] == 3);

	com.gui_icc.same_primaries = same_primaries(current_icc_profile(), com.gui_icc.monitor, com.gui_icc.soft_proof ? com.gui_icc.soft_proof : NULL);

	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (processing_is_job_active()) {
			/* A sequence operation is running — gfit is not being written by
			 * the worker, but notify_gfit_data_modified() has a blanket guard
			 * that skips the remap when any job is active.  Post a small idle
			 * that tries a reader lock and remaps directly if it succeeds.
			 * (For single-image ops the menu is insensitive so we won't reach
			 * here in that case.) */
			siril_add_idle(try_remap_for_mode_change_idle, NULL);
		} else {
			notify_gfit_data_modified();
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}

	/* The display-mode menu is a radio group; once the user picks a mode
	 * the menu has done its job and should dismiss itself. */
	GtkWidget *popover = gtk_widget_get_ancestor(GTK_WIDGET(checkmenuitem),
	                                              GTK_TYPE_POPOVER);
	if (popover)
		gtk_popover_popdown(GTK_POPOVER(popover));
}

void on_mask_enable_toggled(GtkCheckButton *button, gpointer user_data) {
	gboolean state = siril_toggle_get_active(GTK_WIDGET(button));
	gfit->mask_active = state;
	if (com.pref.gui.mask_tints_vports) {
		notify_gfit_data_modified();
		redraw(REMAP_ALL); // draw or remove the red tint from the image
	}
}

void on_mask_show_toggled(GtkCheckButton *button, gpointer user_data) {
	gboolean state = siril_toggle_get_active(GTK_WIDGET(button));
	com.pref.gui.mask_tints_vports = state;
	siril_log_message(state ? _("Mask visibility enabled\n") : _("Mask visibility disabled\n"));
	notify_gfit_data_modified();
	redraw(REMAP_ALL);
}

void on_mask_clear_clicked(GtkButton *button, gpointer user_data) {
	if (!gfit || !gfit->mask) {
		return;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->mask_hook = mask_clear_hook;
	args->description = _("Clear mask");
	args->verbose = TRUE;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

// Callback function to switch tabs when button is clicked (GTK4 GtkGestureClick)
static void on_mask_button_clicked(GtkGestureClick *gesture, int n_press, double x, double y, gpointer data) {
	GtkWidget *button = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	guint btn = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	if (btn == GDK_BUTTON_PRIMARY) { // Left click
		GtkNotebook *notebook = GTK_NOTEBOOK(data);
		int page = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "page-num"));
		gtk_notebook_set_current_page(notebook, page);
	} else if (btn == GDK_BUTTON_SECONDARY) { // Right click
		GtkWidget *popover = GTK_WIDGET(gtk_menu_button_get_popover(GTK_MENU_BUTTON(button)));
		gtk_popover_popup(GTK_POPOVER(popover));
		gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
	}
}

static void initialize_mask_tab_label() {
	// Get the notebook
	GtkNotebook *notebook = GTK_NOTEBOOK(lookup_widget("notebook1"));

	// Find the mask tab (position 4 based on your XML)
	int mask_tab_position = 4;

	// Create the tab label container
	GtkWidget *tab_label_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);

	// Create the menu button with a down arrow.  GTK4 GtkMenuButton is no
	// longer a GtkButton subclass — use the menu-button setter directly.
	GtkWidget *menu_button = gtk_menu_button_new();
	gtk_menu_button_set_label(GTK_MENU_BUTTON(menu_button), "Mask");

	// Create the popover (Phase 15.C: gtk_popover_new takes no args in GTK4)
	GtkWidget *popover = gtk_popover_new();

	// Create a vbox to hold menu items inside the popover
	GtkWidget *menu_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_margin_start(menu_box, 6);
	gtk_widget_set_margin_end(menu_box, 6);
	gtk_widget_set_margin_top(menu_box, 6);
	gtk_widget_set_margin_bottom(menu_box, 6);

	// Create menu items as check buttons
	GtkWidget *enable_check = gtk_check_button_new_with_label(_("Mask active (GUI)"));
	gtk_widget_set_tooltip_text(enable_check, _("Toggle mask on or off for operations carried out using the GUI. Does not apply to Siril commands or python scripts"));
	GtkWidget *show_check = gtk_check_button_new_with_label(_("Show mask as tint"));
	gtk_widget_set_tooltip_text(show_check, _("Show the mask as a red tinted overlay in the image viewports. (Note: the mask can always be viewed in the mask tab)"));

	// Create separator
	GtkWidget *separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_margin_top(separator, 6);
	gtk_widget_set_margin_bottom(separator, 6);

	// Create clear button
	GtkWidget *clear_button = gtk_button_new_with_label(_("Clear Mask"));
	gtk_widget_set_tooltip_text(clear_button, _("Clear the current image mask"));

	// Add items to menu box
	gtk_widget_set_halign(enable_check, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(enable_check, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(menu_box), enable_check);
	gtk_widget_set_halign(show_check, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(show_check, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(menu_box), show_check);
	gtk_widget_set_halign(separator, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(separator, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(menu_box), separator);
	gtk_widget_set_halign(clear_button, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(clear_button, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(menu_box), clear_button);

	// Show all widgets in menu box
	gtk_widget_set_visible(menu_box, TRUE);

	// Add menu box to popover
	gtk_popover_set_child(GTK_POPOVER(popover), menu_box);

	// Set the popover to the menu button
	gtk_menu_button_set_popover(GTK_MENU_BUTTON(menu_button), popover);

	// Pack the menu button into the tab label box
	gtk_widget_set_halign(menu_button, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(menu_button, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(tab_label_box), menu_button);

	// Show all widgets
	gtk_widget_set_visible(menu_button, TRUE);
	gtk_widget_set_visible(tab_label_box, TRUE);

	// Expose widgets to the builder with IDs
	gtk_builder_expose_object(gui.builder, "mask_tab_label_box", G_OBJECT(tab_label_box));
	gtk_builder_expose_object(gui.builder, "mask_menu_button", G_OBJECT(menu_button));
	gtk_builder_expose_object(gui.builder, "mask_enable_check", G_OBJECT(enable_check));
	gtk_builder_expose_object(gui.builder, "mask_show_check", G_OBJECT(show_check));
	gtk_builder_expose_object(gui.builder, "mask_clear_button", G_OBJECT(clear_button));

	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(show_check)), com.pref.gui.mask_tints_vports);
	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(enable_check)), gfit->mask_active);

	// Get the tab content (the vbox_mask)
	GtkWidget *tab_content = gtk_notebook_get_nth_page(notebook, mask_tab_position);

	// Replace the tab label
	gtk_notebook_set_tab_label(notebook, tab_content, tab_label_box);

	// Store the page number in the menu button
	g_object_set_data(G_OBJECT(menu_button), "page-num", GINT_TO_POINTER(mask_tab_position));

	// GTK4: switch tabs on left-click via GtkGestureClick
	{
		GtkGesture *gesture = gtk_gesture_click_new();
		gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gesture), 0); // any button
		g_signal_connect(gesture, "pressed", G_CALLBACK(on_mask_button_clicked), notebook);
		gtk_widget_add_controller(menu_button, GTK_EVENT_CONTROLLER(gesture));
	}

	// Connect signals
	g_signal_connect(enable_check, "toggled", G_CALLBACK(on_mask_enable_toggled), NULL);
	g_signal_connect(show_check, "toggled", G_CALLBACK(on_mask_show_toggled), NULL);
	g_signal_connect(clear_button, "clicked", G_CALLBACK(on_mask_clear_clicked), NULL);
}

void on_autohd_item_toggled(GtkCheckButton *menuitem, gpointer user_data) {
	gui.use_hd_remap = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(menuitem)));
	if (gui.rendering_mode == STF_DISPLAY) {
		if (gui.use_hd_remap) {
			if (gfit->type == DATA_FLOAT)
				allocate_hd_remap_indices();
			siril_log_message(_("The AutoStretch display mode will use a %d bit LUT\n"), (int) log2(gui.hd_remap_max));
		} else {
			hd_remap_indices_cleanup();
			siril_log_message(_("The AutoStretch display mode will use a 16 bit LUT\n"));
		}
		notify_gfit_data_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_button_apply_hd_bitdepth_clicked(GtkSpinButton *button, gpointer user_data) {
	int bitdepth = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_hd_bitdepth")));
	siril_log_debug("bitdepth: %d\n", bitdepth);
	if (gui.hd_remap_max != 1 << bitdepth) {
		siril_log_message(_("Setting HD AutoStretch display mode bit depth to %d...\n"), bitdepth);
		com.pref.hd_bitdepth = bitdepth;
		gui.hd_remap_max = 1 << bitdepth;
		if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type == DATA_FLOAT) {
			allocate_hd_remap_indices();
			notify_gfit_data_modified();
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}
}

/* Sets the display mode combo box to the value stored in the relevant struct.
 * The operation is purely graphical. */
void set_display_mode() {
	/* Phase 15: GtkMenu / GTK_MENU removed in GTK4.  The display_menu
	 * widget tracking is no longer needed — we directly query each
	 * GtkCheckButton by id. */
	GtkCheckButton *button;
	static GtkLabel *label_display_menu = NULL;
	if (!label_display_menu) {
		label_display_menu = GTK_LABEL(lookup_widget("display_button_name"));
	}

	switch (gui.rendering_mode) {
	default:
	case LINEAR_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[LINEAR_DISPLAY]));
		break;
	case LOG_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[LOG_DISPLAY]));
	break;
	case SQRT_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[SQRT_DISPLAY]));
	break;
	case SQUARED_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[SQUARED_DISPLAY]));
		break;
	case ASINH_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[ASINH_DISPLAY]));
		break;
	case STF_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[STF_DISPLAY]));
	break;
	case HISTEQ_DISPLAY:
		button = GTK_CHECK_BUTTON(lookup_widget(display_item_name[HISTEQ_DISPLAY]));
		break;
	}
	g_signal_handlers_block_by_func(button, on_display_item_toggled, NULL);
	siril_toggle_set_active(GTK_WIDGET(button), TRUE);
	g_signal_handlers_unblock_by_func(button, on_display_item_toggled, NULL);

	gtk_label_set_text(label_display_menu, gtk_check_button_get_label(GTK_CHECK_BUTTON(button)));

}

gboolean set_display_mode_idle(gpointer user_data) {
	set_display_mode();
	return FALSE;
}

void set_unlink_channels(gboolean unlinked) {
	siril_log_debug("channels unlinked: %d\n", unlinked);
	gui.unlink_channels = unlinked;
	notify_gfit_data_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}

/* fill the label indicating how many images are selected in the gray and
 * which one is the reference image, at the bottom of the main window */
void adjust_sellabel() {
	static GtkLabel *global_label = NULL, *label_name_of_seq = NULL;
	gchar *buffer_global, *buffer_title;

	if (global_label == NULL) {
		global_label = GTK_LABEL(lookup_widget("labelseq"));
		label_name_of_seq = GTK_LABEL(lookup_widget("label_name_of_seq"));
	}

	if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);

		buffer_global = g_strdup_printf(_("%s, %d/%d images selected"),
				seq_basename, com.seq.selnum, com.seq.number);
		buffer_title = g_strdup(_("Sequence: "));

		g_free(seq_basename);
	} else if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);

		buffer_global = g_strdup_printf("%s", filename);
		buffer_title = g_strdup(_("Image: "));

		g_free(filename);
	} else {
		buffer_global = g_strdup(_("- none -"));
		buffer_title = NULL;

		gtk_widget_set_sensitive(lookup_widget("goregister_button"), FALSE);
	}

	gtk_label_set_text(label_name_of_seq, buffer_title);
	gtk_label_set_text(global_label, buffer_global);

	g_free(buffer_global);
	g_free(buffer_title);
}

static gboolean update_seq_gui_idle(gpointer p) {
	adjust_sellabel();
	update_reg_interface(FALSE);
	update_stack_interface(FALSE);
	redraw(REDRAW_OVERLAY);
	drawPlot();
	writeseqfile(&com.seq);
	return FALSE;
}

gpointer update_seq_gui_idle_thread_func(gpointer data) {
	execute_idle_and_wait_for_it(update_seq_gui_idle, NULL);
	return FALSE;
}

void set_icon_entry(GtkEntry *entry, gchar *string) {

	gtk_entry_set_icon_from_icon_name(entry, GTK_ENTRY_ICON_SECONDARY, string);
	if (string) {
		gchar *text = g_strdup(_("This sequence name already exists!! "
				"Please change the name before converting."));
		gtk_entry_set_icon_tooltip_text (entry, GTK_ENTRY_ICON_SECONDARY, text);

		g_free(text);
	}
}

gboolean update_MenuItem(gpointer user_data) {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	gboolean is_a_single_image_loaded;	/* An image is loaded. Not a sequence or only the result of stacking process */
	gboolean is_a_singleRGB_image_loaded;	/* A RGB image is loaded. Not a sequence or only the result of stacking process */
	gboolean any_image_is_loaded;		/* Something is loaded. Single image or Sequence */

	is_a_single_image_loaded = single_image_is_loaded() && (!sequence_is_loaded() || (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE || com.seq.current == SCALED_IMAGE)));
	is_a_singleRGB_image_loaded = isrgb(gfit) && single_image_is_loaded();
	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	/* some toolbar buttons */
	gtk_widget_set_sensitive(lookup_widget("toolbarbox"), any_image_is_loaded);

	gboolean wcs_ok = any_image_is_loaded && has_wcs(gfit);
	GAction *action_annotate = g_action_map_lookup_action(G_ACTION_MAP(app_win), "annotate-object");
	GAction *action_grid = g_action_map_lookup_action(G_ACTION_MAP(app_win), "wcs-grid");

	/* any image with wcs information is needed */
	siril_window_enable_wcs_proc_actions(app_win, wcs_ok);
	siril_window_enable_wcs_disto_proc_actions(app_win, wcs_ok && gfit->keywords.wcslib->lin.dispre);

	/* untoggle if disabled */
	if (!wcs_ok) {
		GVariant *state = g_action_get_state(action_annotate);
		if (g_variant_get_boolean(g_action_get_state(action_annotate))) {
			g_action_change_state(action_annotate, g_variant_new_boolean(FALSE));
		}
		g_variant_unref(state);
		state = g_action_get_state(action_grid);
		if (g_variant_get_boolean(g_action_get_state(action_grid))) {
			g_action_change_state(action_grid, g_variant_new_boolean(FALSE));
		}
		g_variant_unref(state);
	}

	/* undo and redo */
	GAction *action_undo = g_action_map_lookup_action(G_ACTION_MAP(app_win), "undo");
	GAction *action_redo = g_action_map_lookup_action(G_ACTION_MAP(app_win), "redo");
	g_simple_action_set_enabled(G_SIMPLE_ACTION (action_undo), is_undo_available());
	g_simple_action_set_enabled(G_SIMPLE_ACTION (action_redo), is_redo_available());

	/* update undo/redo button tooltips */
	if (is_undo_available()) {
		historic *h = (historic *) com.undo_stack->data;
		gchar *str = g_strdup_printf(_("Undo: \"%s\""), h->history);
		gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), str);
		g_free(str);
	} else {
		gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), _("Nothing to undo"));
	}
	if (is_redo_available()) {
		historic *h = (historic *) com.redo_stack->data;
		gchar *str = g_strdup_printf(_("Redo: \"%s\""), h->history);
		gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), str);
		g_free(str);
	} else {
		gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), _("Nothing to redo"));
	}

	/* save and save as */
	GAction *action_save = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save");
	GAction *action_save_as = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save-as");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save), is_a_single_image_loaded && com.uniq->fileexist);
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save_as), any_image_is_loaded);

	/* popup menu */

	/* fits header */
	GAction *action_fits_header = g_action_map_lookup_action (G_ACTION_MAP(app_win), "fits-header");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_fits_header), any_image_is_loaded && gfit->header != NULL);
	/* annotate dialog */
	GAction *action_annotate_dialog = g_action_map_lookup_action (G_ACTION_MAP(app_win), "annotate-dialog");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_annotate_dialog), any_image_is_loaded && has_wcs(gfit));
	/* Lightcurve process */
	GAction *action_nina_light_curve = g_action_map_lookup_action (G_ACTION_MAP(app_win), "nina_light_curve");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_nina_light_curve), sequence_is_loaded() && has_wcs(gfit));
	/* selection is needed */
	siril_window_enable_if_selection_actions(app_win, com.selection.w && com.selection.h);
	/* selection and sequence is needed */
	siril_window_enable_if_selection_sequence_actions(app_win, com.selection.w && com.selection.h && (sequence_is_loaded() || valid_rgbcomp_seq()));
	/* selectoin and isrgb is needed */
	siril_window_enable_if_selection_rgb_actions(app_win, com.selection.w && com.selection.h && isrgb(gfit));

	/* Image processing Menu */

	/* single RGB image is needed */
	siril_window_enable_rgb_proc_actions(app_win, is_a_singleRGB_image_loaded);
	/* single RGB image with wcs information is needed */
	siril_window_enable_rgb_wcs_proc_actions(app_win, is_a_singleRGB_image_loaded && has_wcs(gfit));
	/* single or sequence RGB is needed */
	siril_window_enable_any_rgb_proc_actions(app_win, is_a_singleRGB_image_loaded|| sequence_is_loaded());
	/* any image is needed */
	siril_window_enable_any_proc_actions(app_win, any_image_is_loaded);
	/* any mono image is needed */
	siril_window_enable_any_mono_proc_actions(app_win, any_image_is_loaded && !isrgb(gfit));
	/* only single image needed */
	siril_window_enable_single_proc_actions(app_win, is_a_single_image_loaded);
	/* no images needed */
	siril_window_enable_none_proc_actions(app_win, TRUE);

	/* Image information menu */
	siril_window_enable_image_actions(app_win, any_image_is_loaded);

	/* auto-stretch actions */
	siril_window_autostretch_actions(app_win, gui.rendering_mode == STF_DISPLAY && gfit->naxes[2] == 3);

	/* keywords list */
	if (gtk_widget_is_visible(lookup_widget("keywords_dialog")))
		refresh_keywords_dialog(NULL);
	return FALSE;
}

void sliders_mode_set_state(sliders_mode sliders) {
	/* GTK4: the three radiobutton_* widgets are GtkCheckButton in the .ui;
	 * casting through GTK_TOGGLE_BUTTON triggers a runtime
	 * "invalid cast from 'GtkCheckButton' to 'GtkToggleButton'" critical.
	 * Hold a generic GtkWidget* and route through the polymorphic shim. */
	GtkWidget *radiobutton;	// Must not be static
	gchar *str[] = { "radiobutton_hilo", "radiobutton_minmax", "radiobutton_user" };
	void *func[] = { on_radiobutton_hilo_toggled, on_radiobutton_minmax_toggled, on_radiobutton_user_toggled };

	radiobutton = lookup_widget(str[sliders]);

	g_signal_handlers_block_by_func(radiobutton, func[sliders], NULL);
	siril_toggle_set_active(radiobutton, TRUE);
	g_signal_handlers_unblock_by_func(radiobutton, func[sliders], NULL);
}

gboolean sliders_mode_set_state_idle(gpointer p) {
	sliders_mode sliders = *(sliders_mode*) p;
	sliders_mode_set_state(sliders);
	return FALSE;
}

display_mode get_display_mode_from_menu() {
	static GtkWidget *menu = NULL;
	if (!menu)
		menu = lookup_widget("menu_display");
	for (int i = 0; i < G_N_ELEMENTS(display_item_name); i++) {
		gboolean is_active = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget(display_item_name[i]))));
		if (is_active) return (display_mode) i;
	}
	return LINEAR_DISPLAY;
}

void set_initial_display_mode(display_mode x) {
	switch(x) {
		case LINEAR_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("linear_item"))), TRUE);
			break;
		case LOG_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("log_item"))), TRUE);
			break;
		case SQRT_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("square_root_item"))), TRUE);
			break;
		case SQUARED_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("squared_item"))), TRUE);
			break;
		case ASINH_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("asinh_item"))), TRUE);
			break;
		case STF_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("auto_item"))), TRUE);
			break;
		case HISTEQ_DISPLAY:
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("histo_item"))), TRUE);
			break;
	}
}

static void set_initial_histogram_display_mode(int x) {
	GtkCheckButton *button = GTK_CHECK_BUTTON(lookup_widget("HistoCheckLogButton"));
	g_signal_handlers_block_by_func(button, on_histo_toggled, NULL);
	switch(x) {
	case LINEAR_DISPLAY:
		siril_toggle_set_active(GTK_WIDGET(button), FALSE);
		break;
	case LOG_DISPLAY:
		siril_toggle_set_active(GTK_WIDGET(button), TRUE);
	}
	g_signal_handlers_unblock_by_func(button, on_histo_toggled, NULL);
}


void update_prepro_interface(gboolean allow_debayer) {
	static GtkCheckButton *udark = NULL, *uoffset = NULL, *uflat = NULL,
			       *checkAutoEvaluate = NULL;
	static GtkWidget *prepro_button = NULL, *cosme_grid = NULL, *dark_optim = NULL;
       	static GtkWidget *equalize = NULL, *auto_eval = NULL, *flat_norm = NULL;
       	static GtkWidget *debayer = NULL, *fix_xtrans = NULL;
	static GtkDropDown *output_type = NULL;
	if (udark == NULL) {
		udark = GTK_CHECK_BUTTON(lookup_widget("usedark_button"));
		uoffset = GTK_CHECK_BUTTON(lookup_widget("useoffset_button"));
		uflat = GTK_CHECK_BUTTON(lookup_widget("useflat_button"));
		checkAutoEvaluate = GTK_CHECK_BUTTON(lookup_widget("checkbutton_auto_evaluate"));
		output_type = GTK_DROP_DOWN(lookup_widget("prepro_output_type_combo"));
		prepro_button = lookup_widget("prepro_button");
		cosme_grid = lookup_widget("grid24");
		dark_optim = lookup_widget("comboDarkOptimize");
		equalize = lookup_widget("checkbutton_equalize_cfa");
		auto_eval = lookup_widget("checkbutton_auto_evaluate");
		flat_norm = lookup_widget("entry_flat_norm");
		debayer = lookup_widget("checkButton_pp_dem");
		fix_xtrans = lookup_widget("fix_xtrans_af");
	}

	gtk_widget_set_sensitive(prepro_button,
			siril_toggle_get_active(GTK_WIDGET(debayer)) ||
			((sequence_is_loaded() || single_image_is_loaded()) &&
			(siril_toggle_get_active(GTK_WIDGET(udark)) ||
			siril_toggle_get_active(GTK_WIDGET(uoffset)) ||
			siril_toggle_get_active(GTK_WIDGET(uflat)))));

	gtk_widget_set_sensitive(cosme_grid, siril_toggle_get_active(GTK_WIDGET(udark)));

	gtk_widget_set_sensitive(dark_optim, siril_toggle_get_active(GTK_WIDGET(udark)));

	gtk_widget_set_sensitive(equalize, siril_toggle_get_active(GTK_WIDGET(uflat)));

	gtk_widget_set_sensitive(auto_eval, siril_toggle_get_active(GTK_WIDGET(uflat)));

	gtk_widget_set_sensitive(flat_norm,
			siril_toggle_get_active(GTK_WIDGET(uflat)) &&
			!siril_toggle_get_active(GTK_WIDGET(checkAutoEvaluate)));

	gtk_widget_set_sensitive(debayer, allow_debayer);
	gtk_widget_set_sensitive(fix_xtrans, siril_toggle_get_active(GTK_WIDGET(udark)) || siril_toggle_get_active(GTK_WIDGET(uoffset)));

	gtk_widget_set_sensitive(GTK_WIDGET(output_type), sequence_is_loaded());
	int type = com.seq.type;
	if (com.seq.type < 0 || com.seq.type > SEQ_FITSEQ)
		type = SEQ_REGULAR;
	gtk_drop_down_set_selected(output_type, type);
}

void clear_sampling_setting_box() {
	GtkDropDown *binning = GTK_DROP_DOWN(
			gtk_builder_get_object(gui.builder, "combobinning"));
	GtkEntry* focal_entry = GTK_ENTRY(lookup_widget("focal_entry"));
	GtkEntry* pitchX_entry = GTK_ENTRY(lookup_widget("pitchX_entry"));
	GtkEntry* pitchY_entry = GTK_ENTRY(lookup_widget("pitchY_entry"));

	gtk_editable_set_text(GTK_EDITABLE(focal_entry), "");
	gtk_editable_set_text(GTK_EDITABLE(pitchX_entry), "");
	gtk_editable_set_text(GTK_EDITABLE(pitchY_entry), "");
	gtk_drop_down_set_selected(binning, 0);
}

static const char *untranslated_vport_number_to_name(int vport) {
	switch (vport) {
		case RED_VPORT:
			return ("red");
		case GREEN_VPORT:
			return ("green");
		case BLUE_VPORT:
			return ("blue");
		case RGB_VPORT:
			return ("rgb");
	}
	return NULL;
}

const gchar *layer_name_for_gfit(int layer) {
	if (gfit->naxes[2] == 1) {
		if (layer != 0) {
			siril_log_debug("unloaded layer name requested %d\n", layer);
			return "unset";
		}
		return _("Luminance");
	}
	switch (layer) {
		case 0:
			return _("Red");
		case 1:
			return _("Green");
		case 2:
			return _("Blue");
		default:
			siril_log_debug("unloaded layer name requested %d\n", layer);
			return "unset";
	}
}

int match_drawing_area_widget(const GtkWidget *drawing_area, gboolean allow_rgb) {
	/* could be done with a for i=0 loop, to get rid of these defines */
	if (drawing_area == gui.view[RED_VPORT].drawarea)
		return RED_VPORT;
	if (drawing_area == gui.view[GREEN_VPORT].drawarea)
		return GREEN_VPORT;
	else if (drawing_area == gui.view[BLUE_VPORT].drawarea)
		return BLUE_VPORT;
	else if (allow_rgb && drawing_area == gui.view[RGB_VPORT].drawarea)
		return RGB_VPORT;
	else if (drawing_area == gui.view[MASK_VPORT].drawarea)
		return MASK_VPORT;
	return -1;
}

void update_display_selection() {
	static const gchar *label_selection[] = { "labelselection_red", "labelselection_green", "labelselection_blue", "labelselection_rgb", "labelselection_mask" };
	static gchar selection_buffer[256] = { 0 };
	if (com.selection.w && com.selection.h) {
		g_sprintf(selection_buffer, _("W: %dpx H: %dpx ratio: %.4f"), com.selection.w, com.selection.h,
			(double)com.selection.w / (double)com.selection.h);
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_selection[gui.cvport])), selection_buffer);
	} else {
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_selection[gui.cvport])), "");
	}
}

static void update_roi_from_selection() {
	if (gui.roi.active) {
		restore_roi();
		copy_roi_into_gfit();
	}
	memcpy(&gui.roi.selection, &com.selection, sizeof(rectangle));
	gboolean active = (gui.roi.selection.w > 0 && gui.roi.selection.h > 0);
	if (active)
		on_set_roi();
	else
		on_clear_roi();
}

void update_roi_config() {
	if (com.pref.gui.roi_mode == ROI_AUTO)
		register_selection_update_callback(update_roi_from_selection);
	else
		unregister_selection_update_callback(update_roi_from_selection);

	gtk_widget_set_visible(lookup_widget("submenu_gray_roi"), com.pref.gui.roi_mode == ROI_MANUAL);
}

void update_display_fwhm() {
	if (gui.cvport == MASK_VPORT) return;
	static const gchar *label_fwhm[] = { "labelfwhm_red", "labelfwhm_green", "labelfwhm_blue", "labelfwhm_rgb" };
	static gchar fwhm_buffer[256] = { 0 };

	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		g_sprintf(fwhm_buffer, " ");
	} else if (com.selection.w && com.selection.h) {// Now we don't care about the size of the sample. Minimization checks that
		if (com.selection.w < 300 && com.selection.h < 300 && com.selection.w > 5 && com.selection.h > 5) {
			double roundness;
			double fwhm_val = psf_get_fwhm(gfit, select_vport(gui.cvport), &com.selection, &roundness);
			g_sprintf(fwhm_buffer, _("fwhm: %.2f px, r: %.2f"), fwhm_val, roundness);
		} else
			g_sprintf(fwhm_buffer, _("fwhm: N/A"));
	} else {
		g_sprintf(fwhm_buffer, _("fwhm: N/A"));
	}
	gtk_label_set_text(GTK_LABEL(lookup_widget(label_fwhm[gui.cvport])), fwhm_buffer);
}

/* displays the opened image file name in the layers window.
 * if a unique file is loaded, its details are used instead of any sequence data
 */
void display_filename() {
	gboolean local_filename = FALSE;
	int nb_layers, vport;
	char *filename;
	gchar *orig_filename = NULL;
	GError *error = NULL;
	if (single_image_is_loaded() && com.uniq->filename && com.uniq->filename[0] != '\0') {	// unique image
		filename = com.uniq->filename;
		orig_filename = g_file_read_link(filename, &error);
		nb_layers = com.uniq->chans;
	} else if (sequence_is_loaded()) {	// sequence
		filename = malloc(256);
		strcpy(filename, " ");
		local_filename = TRUE;
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		orig_filename = g_file_read_link(filename, &error);
		nb_layers = com.seq.nb_layers;
	} else {
		filename = " ";
		nb_layers = 3;
	}
	vport = nb_layers > 1 ? nb_layers + 1 : nb_layers;

	gchar *base_name = g_path_get_basename(filename);
	if (strlen(base_name) > 4096) // Prevent manifestly excessive lengths
		base_name[4095] = 0;
	gchar *orig_base_name = NULL;
	GString *concat_base_name = NULL;
	if (orig_filename && orig_filename[0] != '\0') {
		orig_base_name = g_path_get_basename(orig_filename); // taints orig_filename
		concat_base_name = g_string_new(orig_base_name);
		if (concat_base_name->len > 4096) // Sanitize length of tainted string
			g_string_truncate(concat_base_name, 4096);
		concat_base_name = g_string_append(concat_base_name, " <==> ");
		concat_base_name = g_string_append(concat_base_name, base_name);
	}
	for (int channel = 0; channel < vport; channel++) {
		gchar *c = g_strdup_printf("labelfilename_%s", untranslated_vport_number_to_name(channel));
		GtkLabel *fn_label = GTK_LABEL(lookup_widget(c));
		if (orig_base_name) {
			gtk_label_set_text(fn_label, concat_base_name->str);
			gtk_label_set_ellipsize(fn_label, PANGO_ELLIPSIZE_START);
		} else {
			gtk_label_set_text(fn_label, base_name);
			gtk_label_set_ellipsize(fn_label, PANGO_ELLIPSIZE_NONE);
		}
		g_free(c);
	}
	GtkLabel *mask_label = GTK_LABEL(lookup_widget("labelfilename_mask"));
	if (orig_base_name) {
		gtk_label_set_text(mask_label, concat_base_name->str);
		gtk_label_set_ellipsize(mask_label, PANGO_ELLIPSIZE_START);
	} else {
		gtk_label_set_text(mask_label, base_name);
		gtk_label_set_ellipsize(mask_label, PANGO_ELLIPSIZE_NONE);
	}

	if (local_filename)
		free(filename);
	g_free(base_name);
	g_free(orig_filename);
	g_free(orig_base_name);
	if (error) g_error_free(error);
	if (concat_base_name) g_string_free(concat_base_name, TRUE);
}

void on_precision_item_toggled(GtkCheckButton *checkmenuitem, gpointer user_data) {
	if (!single_image_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Cannot convert a sequence file"),
				_("A sequence file cannot be converted to 32 bits. This operation can only be done on a single file."));
	} else {
		size_t ndata = gfit->naxes[0] * gfit->naxes[1] * gfit->naxes[2];

		if (gfit->type == DATA_FLOAT) {
			gboolean convert = siril_confirm_dialog(_("Precision loss"),
					_("Converting the image from 32 bits to 16 bits may lead to a loss of numerical accuracy. "
							"Getting back to 32 bits will not recover this loss.\n"
							"Are you sure you want to convert your data?"), _("Convert to 16 bits"));
			if (convert) {
				if (is_preview_active()) {
					fits *backup_gfit = get_preview_gfit_backup();
					fit_replace_buffer(backup_gfit, float_buffer_to_ushort(backup_gfit->fdata, ndata), DATA_USHORT);
				}
				fit_replace_buffer(gfit, float_buffer_to_ushort(gfit->fdata, ndata), DATA_USHORT);
				invalidate_gfit_histogram();
				update_gfit_histogram_if_needed();
				notify_gfit_data_modified();
				redraw(REMAP_ALL);
			}
		} else if (gfit->type == DATA_USHORT) {
			if (is_preview_active()) {
				fits *backup_gfit = get_preview_gfit_backup();
				fit_replace_buffer(backup_gfit, ushort_buffer_to_float(backup_gfit->data, ndata), DATA_FLOAT);
			}
			fit_replace_buffer(gfit, ushort_buffer_to_float(gfit->data, ndata), DATA_FLOAT);
			invalidate_gfit_histogram();
			update_gfit_histogram_if_needed();
			notify_gfit_data_modified();
			redraw(REMAP_ALL);
		}
	}
	gui_function(set_precision_switch, NULL);
}

gboolean set_precision_switch(gpointer user_data) {
	if (!com.script) {
		GtkLabel *label = GTK_LABEL(lookup_widget("precision_button_name"));
		GtkCheckButton *float_button = GTK_CHECK_BUTTON(lookup_widget("32bits_item"));
		GtkCheckButton *ushort_button = GTK_CHECK_BUTTON(lookup_widget("16bits_item"));

		gtk_label_set_text(label, gfit->type == DATA_USHORT ? _("16 bits") : _("32 bits"));
		g_signal_handlers_block_by_func(float_button, on_precision_item_toggled, NULL);
		siril_toggle_set_active(GTK_WIDGET(float_button), gfit->type == DATA_FLOAT);
		siril_toggle_set_active(GTK_WIDGET(ushort_button), gfit->type == DATA_USHORT);
		g_signal_handlers_unblock_by_func(float_button,	on_precision_item_toggled, NULL);
	}
	return FALSE;
}

/* updates the combo box of registration layers to reflect data availability */
int set_layers_for_registration() {
	static GtkDropDown *cbbt_layers = NULL;
	int i;
	int reminder;

	if (cbbt_layers == NULL)
		cbbt_layers = GTK_DROP_DOWN(
				gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	reminder = gtk_drop_down_get_selected(GTK_DROP_DOWN(cbbt_layers));

	siril_drop_down_clear_strings(cbbt_layers);
	for (i = 0; i < com.seq.nb_layers; i++) {
		gchar *layer;
		const gchar *layer_name = layer_name_for_gfit(i);
		layer = g_strdup_printf("%d: %s", i, layer_name);
		if (com.seq.regparam && com.seq.regparam[i]) {
			str_append(&layer,  " (*)");
			if (reminder == -1 || ((reminder >= 0) && !(com.seq.regparam[reminder]))) // set as default selection
				reminder = i;
		}
		siril_drop_down_append_text(cbbt_layers, layer);

		g_free(layer);
	}
	/* First initialization */
	if (reminder == -1) {
		if (com.seq.nb_layers == 3)
			gtk_drop_down_set_selected(GTK_DROP_DOWN(cbbt_layers), 1);
		else
			gtk_drop_down_set_selected(GTK_DROP_DOWN(cbbt_layers), 0);
	}
	/* Already initialized or default selection to channel with data */
	else
		gtk_drop_down_set_selected(GTK_DROP_DOWN(cbbt_layers), reminder);

	return reminder;
}

void show_data_dialog(char *text, char *title, gchar *parent, gchar *url) {
	GtkTextView *tv = GTK_TEXT_VIEW(lookup_widget("data_txt"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tv);
	GtkTextIter itDebut;
	GtkTextIter itFin;
	GtkWindow *win = GTK_WINDOW(lookup_widget("data_dialog"));

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, text, strlen(text));
	gtk_window_set_title(win, title);

	gtk_widget_set_visible(lookup_widget("data_dialog"), TRUE);
	gtk_widget_set_visible(lookup_widget("data_extra_button"), url != NULL);
	if (url) {
		gtk_link_button_set_uri((GtkLinkButton *)lookup_widget("data_extra_button"), url);
	}
	if (parent) {
		gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget(parent)));
	} else {
		gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget("control_window")));
	}
}

/**
 * Get the active window on toplevels
 * @return the GtkWindow activated
 */
GtkWindow *siril_get_active_window() {
	GtkWindow *win = NULL;
	GList *list = gtk_window_list_toplevels();

	for (GList *l = list; l; l = l->next) {
		if (gtk_window_is_active((GtkWindow *) l->data)) {
			win = (GtkWindow *) l->data;
			break;
		}
	}
	g_list_free(list);
	return win;
}

void initialize_FITS_name_entries() {
	GtkEntry *mbias, *mdark, *mflat, *final_stack;
	gchar *str[4] = { NULL };

	mbias = GTK_ENTRY(lookup_widget("offsetname_entry"));
	mdark = GTK_ENTRY(lookup_widget("darkname_entry"));
	mflat = GTK_ENTRY(lookup_widget("flatname_entry"));
	final_stack = GTK_ENTRY(lookup_widget("entryresultfile"));

	if (com.pref.prepro.bias_lib) {
		if (com.pref.prepro.use_bias_lib) {
			str[0] = g_strdup_printf("%s", com.pref.prepro.bias_lib);
		}
	}

	if (com.pref.prepro.dark_lib) {
		if (com.pref.prepro.use_dark_lib) {
			str[1] = g_strdup_printf("%s", com.pref.prepro.dark_lib);
		}
	}

	if (com.pref.prepro.flat_lib) {
		if (com.pref.prepro.use_flat_lib) {
			str[2] = g_strdup_printf("%s", com.pref.prepro.flat_lib);
		}

	}

	if (com.pref.prepro.stack_default) {
		if (com.pref.prepro.use_stack_default) {
			str[3] = g_strdup_printf("%s", com.pref.prepro.stack_default);
		}

	}

	const char *t_bias = gtk_editable_get_text(GTK_EDITABLE(mbias));
	if (!str[0] && t_bias[0] == '\0') {
		str[0] = g_strdup_printf("master-bias%s", com.pref.ext);
	}
	if (str[0]) {
		gtk_editable_set_text(GTK_EDITABLE(mbias), str[0]);
	}

	const char *t_dark = gtk_editable_get_text(GTK_EDITABLE(mdark));
	if (!str[1] && t_dark[0] == '\0') {
		str[1] = g_strdup_printf("master-dark%s", com.pref.ext);
	}
	if (str[1]) {
		gtk_editable_set_text(GTK_EDITABLE(mdark), str[1]);
	}

	const char *t_flat = gtk_editable_get_text(GTK_EDITABLE(mflat));
	if (!str[2] && t_flat[0] == '\0') {
		str[2] = g_strdup_printf("master-flat%s", com.pref.ext);
	}
	if (str[2]) {
		gtk_editable_set_text(GTK_EDITABLE(mflat), str[2]);
	}

	const char *t_stack = gtk_editable_get_text(GTK_EDITABLE(final_stack));
	if (!str[3] && t_stack[0] == '\0') {
		if (sequence_is_loaded())
			str[3] = g_strdup_printf("%s_stacked%s", com.seq.seqname, com.pref.ext);
		else
			str[3] = g_strdup_printf("stack_result%s", com.pref.ext);
	}
	if (str[3]) {
		gtk_editable_set_text(GTK_EDITABLE(final_stack), str[3]);
	}

	for (int i = 0; i < 4; i++)
		g_free(str[i]);
}

/* when a sequence is loaded, the processing (stacking) output file name is
 * modified to include the name of the sequence */
void set_output_filename_to_sequence_name() {
	static GtkEntry *output_file = NULL;
	gchar *msg;
	if (!com.seq.seqname || *com.seq.seqname == '\0')
		return;
	output_file = GTK_ENTRY(lookup_widget("entryresultfile"));
	if (com.pref.prepro.stack_default) {
		if (com.pref.prepro.use_stack_default) {
			gtk_editable_set_text(GTK_EDITABLE(output_file), com.pref.prepro.stack_default);
			return;
		}
	}
	const gchar* com_ext = get_com_ext(com.seq.fz);

	msg = g_strdup_printf("%s%sstacked%s", com.seq.seqname,
			g_str_has_suffix(com.seq.seqname, "_") ?
			"" : (g_str_has_suffix(com.seq.seqname, "-") ? "" : "_"), com_ext);
	gtk_editable_set_text(GTK_EDITABLE(output_file), msg);

	g_free(msg);
}

gboolean show_or_hide_mask_tab_idle(gpointer p) {
	if (com.headless) return FALSE;
	/* Use trylock: if a processing job holds the write lock (e.g. during
	 * application shutdown) we skip the update rather than deadlocking. */
	if (!g_rw_lock_reader_trylock(&gfit->rwlock))
		return FALSE;
	gboolean has_mask = (gfit->mask != NULL);
	g_rw_lock_reader_unlock(&gfit->rwlock);
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget *page = gtk_notebook_get_nth_page(Color_Layers, MASK_VPORT);
	if (has_mask) {
		gtk_widget_set_visible(page, TRUE);
	} else {
		gtk_widget_set_visible(page, FALSE);
	}
	return FALSE;
}

void show_or_hide_mask_tab() {
	siril_add_pythonsafe_idle(show_or_hide_mask_tab_idle, NULL);
	return;
}

gboolean close_tab(gpointer user_data) {
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget* page;

	/* FLIS: a mono base layer (gfit->naxes[2] == 1) does NOT mean the
	 * displayed image is mono — the composite is chromatic if any layer
	 * is RGB or any mono layer carries a non-greyscale LAYER_COLOR tint
	 * (e.g. a 1-mono-base + 3-tinted-mono SHO stack).  Consult the layer
	 * stack instead of the active layer for FLIS images.  Plain FITS
	 * keep the original gfit-based check. */
	gboolean is_mono = is_current_image_flis()
	    ? !flis_composite_is_chromatic()
	    : (com.seq.nb_layers == 1 || gfit->naxes[2] == 1);

	if (is_mono) {
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_set_visible(page, FALSE);
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_set_visible(page, FALSE);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_set_visible(page, FALSE);
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("B&W"));
	} else {
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("Red"));
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_set_visible(page, TRUE);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_set_visible(page, TRUE);
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_set_visible(page, TRUE);
	}
	show_or_hide_mask_tab_idle(NULL);
	return FALSE;
}

void activate_tab(int vport) {
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook1"));
	if (gtk_notebook_get_current_page(notebook) != vport)
		gtk_notebook_set_current_page(notebook, vport);
	// com.cvport is set in the event handler for changed page
}

gboolean init_right_tab(gpointer user_data) {
	/* FLIS: gfit is the active layer, which may be mono even when the
	 * composite is chromatic (e.g. a mono-base FLIS with a tinted mono
	 * narrowband layer on top).  Read the composite's chromaticity from
	 * the layer stack so the user gets RGB viewports for anything that
	 * displays in colour. */
	gboolean rgb = is_current_image_flis()
	    ? flis_composite_is_chromatic()
	    : isrgb(gfit);
	activate_tab(rgb ? RGB_VPORT : RED_VPORT);
	return FALSE;
}

gboolean update_spinCPU(gpointer user_data) {
	int max = GPOINTER_TO_INT(user_data);
	static GtkSpinButton *spin_cpu = NULL;

	if (spin_cpu == NULL) {
		spin_cpu = GTK_SPIN_BUTTON(lookup_widget("spinCPU"));
	}
	if (max > 0) {
		gtk_spin_button_set_range (spin_cpu, 1, (gdouble) max);
	}
	gtk_spin_button_set_value (spin_cpu, (gdouble) com.max_thread);
	return FALSE;
}

/*****************************************************************************
 *             I N I T I A L I S A T I O N      F U N C T I O N S            *
 ****************************************************************************/

/* Activate handler for "<ns>.page" — sets the popover's GtkStack to a
 * named page.  user_data is the GtkStack widget itself. */
static void on_menu_page_activate(GSimpleAction *action, GVariant *parameter,
                                  gpointer user_data) {
	(void) action;
	GtkStack *stack = GTK_STACK(user_data);
	if (!parameter) return;
	const char *name = g_variant_get_string(parameter, NULL);
	gtk_stack_set_visible_child_name(stack, name);
}

/* When a stacked menu popover closes, jump its stack back to the "root"
 * page so the next time the user opens the menu they land on the top
 * level, not the submenu they were last on.  Transition is suppressed
 * during the reset so the slide-back isn't visible if the popover is
 * still mid-fade. */
static void on_menu_popover_closed(GtkPopover *popover, gpointer user_data) {
	(void) popover;
	GtkStack *stack = GTK_STACK(user_data);
	if (!GTK_IS_STACK(stack)) return;
	GtkStackTransitionType t = gtk_stack_get_transition_type(stack);
	gtk_stack_set_transition_type(stack, GTK_STACK_TRANSITION_TYPE_NONE);
	gtk_stack_set_visible_child_name(stack, "root");
	gtk_stack_set_transition_type(stack, t);
}

/* The menu popovers in siril.ui are GtkBox + GtkButton (rather than
 * GMenu / GtkPopoverMenu).  A button created with <property name="label">
 * holds a single GtkLabel child whose xalign defaults to 0.5, so the
 * text sits centred — which makes the rows look unlike menu items.
 * CSS can't reach GtkLabel::xalign, so we walk the tree at startup and
 * left-align every plain-label button.  Submenu-nav buttons (whose
 * child is a GtkBox containing label + arrow) and check buttons
 * (already laid out left-to-right by GTK) are skipped automatically. */
static void align_button_label_left_recursive(GtkWidget *widget) {
	if (GTK_IS_BUTTON(widget)
	    && !GTK_IS_CHECK_BUTTON(widget)
	    && !GTK_IS_TOGGLE_BUTTON(widget)
	    && !GTK_IS_MENU_BUTTON(widget)) {
		GtkWidget *child = gtk_button_get_child(GTK_BUTTON(widget));
		if (GTK_IS_LABEL(child)) {
			gtk_label_set_xalign(GTK_LABEL(child), 0.0f);
			gtk_widget_set_halign(child, GTK_ALIGN_FILL);
			gtk_widget_set_hexpand(child, TRUE);
		}
	}
	for (GtkWidget *child = gtk_widget_get_first_child(widget); child;
	     child = gtk_widget_get_next_sibling(child)) {
		align_button_label_left_recursive(child);
	}
}

static void align_menu_popover_labels_left(void) {
	static const char *const menu_popovers[] = {
		"main_menu", "menu_display", "menu_plot", "menu_seq_clean",
		"menugray", "menumask_rmb", "processing_menu", "tools_menu",
		"proofing_menu", "snapshot_menu", "precision_menu",
		"clear_menu", "display_plot_menu", "output_plot_menu",
		NULL
	};
	for (int i = 0; menu_popovers[i]; i++) {
		GtkWidget *p = lookup_widget(menu_popovers[i]);
		if (p) align_button_label_left_recursive(p);
	}
}

/* Install a "<action_namespace>.page(s)" action onto the named popover so
 * that GtkModelButtons inside the popover can drive the submenu stack via
 *   action-name="<ns>.page"   action-target="'<page-name>'"
 * The popover must have been built with a GtkStack id="<popover>_stack"
 * as its child (see uifiles/siril.ui menu rewrite).
 *
 * Also walks every descendant GtkButton in the popover and, for those
 * whose action-name is *not* the page-navigation action, connects a
 * "clicked" handler that pops the popover down — this restores the
 * "click an item, the menu closes" behaviour that GtkModelButton
 * provided in GTK3 (we replaced GtkModelButton with GtkButton on the
 * navigation buttons to suppress its auto-popdown; that side effect
 * also affected the action items, which is the bug we patch here). */
static void on_menu_button_popdown_clicked(GtkButton *btn, gpointer user_data) {
	(void) btn;
	GtkPopover *popover = GTK_POPOVER(user_data);
	if (popover)
		gtk_popover_popdown(popover);
}

static void connect_popdown_handlers_recursive(GtkWidget *widget,
                                               GtkPopover *popover,
                                               const char *nav_action_name) {
	if (GTK_IS_BUTTON(widget)) {
		const char *action = NULL;
		if (GTK_IS_ACTIONABLE(widget))
			action = gtk_actionable_get_action_name(GTK_ACTIONABLE(widget));
		/* Navigation buttons drive the stack (action-name == "<ns>.page");
		 * those must NOT pop the popover down.  Action items (anything
		 * else, including no action) should dismiss the menu when clicked. */
		if (!action || g_strcmp0(action, nav_action_name) != 0) {
			g_signal_connect_after(widget, "clicked",
			                       G_CALLBACK(on_menu_button_popdown_clicked),
			                       popover);
		}
	}
	for (GtkWidget *child = gtk_widget_get_first_child(widget);
	     child; child = gtk_widget_get_next_sibling(child)) {
		connect_popdown_handlers_recursive(child, popover, nav_action_name);
	}
}

void attach_display_scale_release_handlers(void);  /* forward — defined further down */

static void wire_popover_menu_stack(const char *popover_id, const char *action_namespace) {
	GtkPopover *popover = GTK_POPOVER(lookup_widget(popover_id));
	if (!popover) {
		g_warning("wire_popover_menu_stack: popover '%s' not found", popover_id);
		return;
	}
	gchar *stack_id = g_strdup_printf("%s_stack", popover_id);
	GtkStack *stack = GTK_STACK(lookup_widget(stack_id));
	g_free(stack_id);
	if (!stack) {
		g_warning("wire_popover_menu_stack: %s_stack not found inside %s",
		         popover_id, popover_id);
		return;
	}
	GSimpleAction *page_action = g_simple_action_new("page", G_VARIANT_TYPE_STRING);
	g_signal_connect(page_action, "activate",
	                 G_CALLBACK(on_menu_page_activate), stack);
	GSimpleActionGroup *group = g_simple_action_group_new();
	g_action_map_add_action(G_ACTION_MAP(group), G_ACTION(page_action));
	g_object_unref(page_action);
	gtk_widget_insert_action_group(GTK_WIDGET(popover), action_namespace,
	                               G_ACTION_GROUP(group));
	g_object_unref(group);

	/* After-the-fact popdown wiring for action buttons. */
	gchar *nav = g_strdup_printf("%s.page", action_namespace);
	connect_popdown_handlers_recursive(GTK_WIDGET(popover), popover, nav);
	g_free(nav);

	/* Reset the stack to the "root" page whenever this popover closes,
	 * so re-opening the menu always starts at the top level. */
	g_signal_connect(popover, "closed",
	                 G_CALLBACK(on_menu_popover_closed), stack);
}

static void load_accels() {
	static const gchar * const accelmap[] = {
		"app.quit",                   "<Primary>q", NULL,
		"app.preferences",            "<Primary>p", NULL,
		"app.open",                   "<Primary>o", NULL,
		"app.save",                   "<Primary>s", NULL,
		"app.save-as",                "<Primary><Shift>s", NULL,

		"win.undo",                   "<Primary>z", NULL,
		"win.redo",                   "<Primary><Shift>z", NULL,
		"win.close",                  "<Primary>w", NULL,
		"win.cwd",                    "<Primary>d", NULL,
		"win.full-screen",            "<Primary>f", NULL,

		"win.conversion",             "F1", NULL,
		"win.sequence",               "F2", NULL,
		"win.prepro",                 "F3", NULL,
		"win.registration",           "F4", NULL,
		"win.plot",                   "F5", NULL,
		"win.stacking",               "F6", NULL,
		"win.logs",                   "F7", NULL,

		"win.hide-show-toolbar",      "<Primary>T", NULL,

		"win.zoom-out",               "<Primary>minus", "<Primary>KP_Subtract", NULL,
		"win.zoom-in",                "<Primary>plus", "<Primary>KP_Add", NULL,
		"win.zoom-fit",               "<Primary>0", NULL,
		"win.zoom-one",               "<Primary>1", NULL,

		"win.annotate-dialog",        "<Primary>slash", NULL,
		"win.astrometry",             "<Primary><Shift>a", NULL,
		"win.histo_display",          "<Primary>h", NULL,
		"win.pcc-processing",         "<Primary><Shift>p", NULL,
		"win.spcc-processing",        "<Primary><Shift>c", NULL,
		"win.compstars",              "<Primary><Shift>b", NULL,
		"win.catmag",                 "<Primary><Shift>m", NULL,
		"win.nina_light_curve",       "<Primary><Shift>n", NULL,
		"win.pickstar",               "<Primary>space", NULL,
		"win.dyn-psf",                "<Primary>F6", NULL,
		"win.clipboard",              "<Primary><Shift>x", NULL,
		"win.statistics",             "<alt>s", NULL,

		"win.negative-processing",    "<Primary>i", NULL,
		"win.rotation90-processing",  "<Primary>Right", NULL,
		"win.rotation270-processing", "<Primary>Left", NULL,
		"win.mirrorx-processing",     "<Primary>Up", NULL,
		"win.mirrory-processing",     "<Primary>Down", NULL,

		"win.regframe",               "<Primary>r", NULL,
		"win.show-layers",            "<Primary>l", NULL,

		NULL /* Terminating NULL */
	};

	set_accel_map(accelmap);
}

void set_accel_map(const gchar * const *accelmap) {
	GApplication *application = g_application_get_default();

	for (const gchar *const *it = accelmap; it[0]; it += g_strv_length((gchar**) it) + 1) {
		gtk_application_set_accels_for_action(GTK_APPLICATION(application), it[0], &it[1]);
	}
}

/* Initialize the rendering mode from the GUI */
void initialize_display_mode() {
	gui.rendering_mode = get_display_mode_from_menu();
}

/* initialize non-preferences GUI from settings */
static void init_GUI_from_settings() {
	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("demosaicingButton"))), com.pref.debayer.open_debayer);
}

gboolean set_GUI_CWD(gpointer user_data) {
	if (!com.wd)
		return FALSE;
	/* GTK4: gtk_header_bar_set_title / set_subtitle were removed.  The
	 * replacement is set_title_widget with a custom GtkLabel containing
	 * markup for the two lines (title + dimmed subtitle).
	 *
	 * For a FLIS image the subtitle becomes "name.flis — layer 'X' [n/N]"
	 * so the active layer and stack position are visible at a glance.
	 * For sequence / plain-FITS / no-image we keep the working directory. */
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("headerbar"));
	gchar *subtitle = NULL;
	if (single_image_is_loaded() && is_current_image_flis() && com.uniq) {
		guint n_layers = g_slist_length(com.uniq->layers);
		flis_layer_t *active = flis_active_layer();
		const gchar *basename = NULL;
		gchar *bn = NULL;
		if (com.uniq->filename && *com.uniq->filename) {
			bn = g_path_get_basename(com.uniq->filename);
			basename = bn;
		} else {
			basename = _("(untitled FLIS)");
		}
		if (active && n_layers > 0) {
			subtitle = g_strdup_printf(
				"%s — layer '%s' [%d/%u]",
				basename,
				active->layer_name ? active->layer_name : "?",
				com.uniq->active_layer + 1, n_layers);
		} else {
			subtitle = g_strdup_printf("%s — FLIS", basename);
		}
		g_free(bn);
	}
	if (!subtitle)
		subtitle = siril_truncate_str(com.wd, 50);
	gchar *markup = g_markup_printf_escaped(
		"<b>Siril-%s</b>\n<small>%s</small>", VERSION, subtitle);
	GtkWidget *title = gtk_label_new(NULL);
	gtk_label_set_markup(GTK_LABEL(title), markup);
	gtk_label_set_justify(GTK_LABEL(title), GTK_JUSTIFY_CENTER);
	gtk_header_bar_set_title_widget(bar, title);
	g_free(markup);
	g_free(subtitle);
	return FALSE;
}

static void initialize_preprocessing() {
	GtkCheckButton *cfaButton, *eqButton, *xtransButton;

	cfaButton = GTK_CHECK_BUTTON(lookup_widget("cosmCFACheck"));
	siril_toggle_set_active(GTK_WIDGET(cfaButton), com.pref.prepro.cfa);
	eqButton = GTK_CHECK_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	siril_toggle_set_active(GTK_WIDGET(eqButton), com.pref.prepro.equalize_cfa);
	xtransButton = GTK_CHECK_BUTTON(lookup_widget("fix_xtrans_af"));
	siril_toggle_set_active(GTK_WIDGET(xtransButton), com.pref.prepro.fix_xtrans);

	update_prepro_interface(FALSE);
}

void update_sampling_in_information() {
	GtkLabel *sampling_label = GTK_LABEL(lookup_widget("info_sampling_label"));
	// display sampling based on information displayed in the window
	// otherwise, we could use get_wcs_image_resolution(fits *fit)
	if (gfit->keywords.focal_length > 0.0 && gfit->keywords.pixel_size_x > 0.0) {
		float pixel_size = gfit->keywords.pixel_size_x;
		if (gfit->keywords.binning_x > 1 && com.pref.binning_update)
			pixel_size *= gfit->keywords.binning_x;
		double scale = get_resolution(gfit->keywords.focal_length, pixel_size);
		char buf[20];
		sprintf(buf, "%.4f", scale);
		gtk_label_set_text(sampling_label, buf);
	}
	else gtk_label_set_text(sampling_label, _("unknown"));
}

/* update the Information window from gfit */
void set_GUI_CAMERA() {
	/* information will be taken from gfit or from the preferences and overwritten in
	 * gfit by the callbacks */
	char buf[20];
	double fl = gfit->keywords.focal_length > 0.0 ? gfit->keywords.focal_length : com.pref.starfinder_conf.focal_length;
	sprintf(buf, "%.3f", fl);
	gtk_editable_set_text(GTK_EDITABLE(lookup_widget("focal_entry")), buf);

	double pixsz = gfit->keywords.pixel_size_x > 0.0 ? gfit->keywords.pixel_size_x : com.pref.starfinder_conf.pixel_size_x;
	sprintf(buf, "%.2f", pixsz);
	gtk_editable_set_text(GTK_EDITABLE(lookup_widget("pitchX_entry")), buf);

	pixsz = gfit->keywords.pixel_size_y > 0.0 ? gfit->keywords.pixel_size_y : com.pref.starfinder_conf.pixel_size_x;
	sprintf(buf, "%.2f", pixsz);
	gtk_editable_set_text(GTK_EDITABLE(lookup_widget("pitchY_entry")), buf);

	GtkDropDown *binning = GTK_DROP_DOWN(lookup_widget("combobinning"));
	if (!gfit->keywords.binning_x || !gfit->keywords.binning_y) {
		gtk_drop_down_set_selected(binning, 0);
	}
	else if (gfit->keywords.binning_x == gfit->keywords.binning_y) {
		/* squared binning */
		gtk_drop_down_set_selected(binning, (gint) gfit->keywords.binning_x - 1);
	} else {
		short coeff =
			gfit->keywords.binning_x > gfit->keywords.binning_y ?
			gfit->keywords.binning_x / gfit->keywords.binning_y :
			gfit->keywords.binning_y / gfit->keywords.binning_x;
		switch (coeff) {
			case 2:
				gtk_drop_down_set_selected(binning, 4);
				break;
			case 3:
				gtk_drop_down_set_selected(binning, 5);
				break;
			default:
				siril_log_message(_("Binning %d x %d is not supported\n"),
						gfit->keywords.binning_x, gfit->keywords.binning_y);
		}
	}

	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("toggleButtonUnbinned"))), com.pref.binning_update);
	//siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("saveinfo_toggle"))), com.pref.astrometry.update_default_scale);
}

/* Phase 20 / 8B: GTK4 GtkDropTarget no longer uses GtkTargetEntry; the
 * drawing-area drop handler is configured directly in
 * single_image.c::siril_drag_single_image_set_dest with GDK_TYPE_FILE_LIST.
 * The orphan declaration that lived here is removed. */

static gboolean on_control_window_configure_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	gui_function(save_main_window_state, NULL);
	return FALSE;
}

static gboolean on_control_window_window_state_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	gui_function(save_main_window_state, NULL);
	return FALSE;
}

static gboolean paned_first_resize = TRUE;

static void pane_notify_position_cb(GtkPaned *paned, gpointer user_data) {
	if (paned_first_resize) {
		if (com.pref.gui.remember_windows && com.pref.gui.pan_position > 0) {
			gtk_paned_set_position(paned, com.pref.gui.pan_position);
		}
		paned_first_resize = FALSE;
	}
}

/* Force paned position restore (called after window is shown) */
void force_paned_restore() {
	if (!com.script && com.pref.gui.remember_windows && com.pref.gui.pan_position > 0) {
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		gtk_paned_set_position(paned, com.pref.gui.pan_position);
	}
}

gboolean on_main_panel_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	int position = gtk_paned_get_position((GtkPaned*) widget);
	if (com.pref.gui.remember_windows) {
		com.pref.gui.pan_position = position;
	}
	int max_position;
	g_object_get(G_OBJECT(widget), "max-position", &max_position, NULL);
	if (position == max_position) {
		GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
		com.pref.gui.pan_position = -1;
		// hide it
		GAction *action_panel = g_action_map_lookup_action(G_ACTION_MAP(app_win), "panel");
		g_action_activate(action_panel, NULL);
		gtk_paned_set_position((GtkPaned*) widget, -1);	// reset to default
	}

	return FALSE;
}

static gboolean gui_ready = FALSE;

gboolean is_gui_ready() {
	return gui_ready;
}

gpointer update_spcc(gpointer user_data) {
	// 1. Update the repository
	if (is_online()) {
#ifdef HAVE_LIBGIT2
		auto_update_gitspcc(TRUE);
#endif
		// 2. Update the SPCC combos
		// populate_spcc_combos_async runs in this thread but the actual call to
		// populate_spcc_combos is run from it in an idle in the GTK thread
		populate_spcc_combos_async(NULL);
	}
	return GINT_TO_POINTER(0);
}

// Updates the repository and then refreshes the scripts menu.
gpointer update_scripts(gpointer user_data) {
	if (is_online()) {
#ifdef HAVE_LIBGIT2
		// 1. Update the repository
//		if (com.python_init_thread)
//			g_thread_join(com.python_init_thread);
//		com.python_init_thread = NULL;
		auto_update_gitscripts(TRUE);
#endif
		// 2. Update the menu (not verbose)
		gui_mutex_lock();
		// refresh_script_menu runs in an idle in the GTK thread
		execute_idle_and_wait_for_it(refresh_script_menu_idle, GINT_TO_POINTER(0));
		gui_mutex_unlock();
	}
	return GINT_TO_POINTER(0);
}

gpointer initialize_script_menu_and_spcc_widgets_serially(gpointer user_data) {
	// First fill the SPCC combos and script menu so they are immediately available
	execute_idle_and_wait_for_it(initialize_script_menu_idle, GINT_TO_POINTER(1)); // script menu first as it's quick based on an already loaded list
	populate_spcc_combos_async(NULL); // this after, as it takes a little time to load the actual files
	// Now, in threads, update the repositories and then update the combos and menu again to reflect any changes
	if (com.pref.auto_script_update && is_online())
		g_thread_unref(g_thread_new("update_scripts", update_scripts, NULL)); // this is slow as will require repository syncing
	if (com.pref.spcc.auto_spcc_update && is_online())
		g_thread_unref(g_thread_new("update_spcc", update_spcc, NULL)); // this is slow as will require repository syncing
	return FALSE;
}

// Only called from preferences when reactivating the SPCC repository
gpointer initialize_spcc(gpointer user_data) {
	// 1. Initialize the SPCC combos
	// populate_spcc_combos_async runs in this thread but the actual call to
	// populate_spcc_combos is run from it in an idle in the GTK thread
	populate_spcc_combos_async(NULL); // controls the GUI mutex
	// 2. Update the repository
	if (com.pref.spcc.auto_spcc_update && is_online()) {
#ifdef HAVE_LIBGIT2
		auto_update_gitspcc(TRUE);
#endif
		// 3. Update the SPCC combos
		// populate_spcc_combos_async runs in this thread but the actual call to
		// populate_spcc_combos is run from it in an idle in the GTK thread
		populate_spcc_combos_async(NULL); // controls the GUI mutex again
	}
	return GINT_TO_POINTER(0);
}

// Initializes the scripts menu, updates the repository and then refreshes the scripts menu
gpointer initialize_scripts(gpointer user_data) {
	// 1. Initialize the menu (verbose)
	// initialize_script_menu_idle runs as an idle in the GTK thread
	gui_mutex_lock();
	execute_idle_and_wait_for_it(initialize_script_menu_idle, GINT_TO_POINTER(1));
	gui_mutex_unlock();
	// 2. Update the repository
	if (com.pref.auto_script_update && is_online()) {
#ifdef HAVE_LIBGIT2
		auto_update_gitscripts(TRUE);
#endif
		// 3. Update the menu (not verbose)
		// refresh_script_menu runs in an idle in the GTK thread
		gui_mutex_lock();
		execute_idle_and_wait_for_it(refresh_script_menu_idle, GINT_TO_POINTER(0));
		gui_mutex_unlock();
	}
	return GINT_TO_POINTER(0);
}
gboolean first_start_cb(gpointer user_data) {
	if (g_strcmp0(com.pref.gui.first_start, PACKAGE_VERSION)) {
		com.pref.gui.first_start = g_strdup(PACKAGE_VERSION);
		writeinitfile();

		gchar *ver = g_strdup_printf(_("Welcome to %s"), PACKAGE_STRING);
		int ret = siril_confirm_dialog(ver,
				_("Hello, this is the first time you use this new version of Siril. Please, have a seat and take the time "
				  "to watch the short introduction we have prepared for you. "
				  "Be aware you can replay this introduction at any times in the Miscellaneous tab of the preferences dialog box."),
				_("See Introduction"));
		if (ret)
			start_intro_script();
		g_free(ver);
	}
	return G_SOURCE_REMOVE;
}

/* GTK4-only GApplication-level actions.  Registered from
 * siril_app_startup via the per-toolkit hook so main.c never has to
 * reference the open_recent_action_activate symbol (which has no GTK3
 * counterpart). */
void register_toolkit_app_actions(GApplication *app) {
	static const GActionEntry entries[] = {
		{ "open-recent", open_recent_action_activate, "s", NULL, NULL },
	};
	g_action_map_add_action_entries(G_ACTION_MAP(app), entries,
			G_N_ELEMENTS(entries), app);
}

void initialize_all_GUI(gchar *supported_files) {
	/* initializing internal structures with widgets (drawing areas) */
	gui.view[RED_VPORT].drawarea  = lookup_widget("drawingarear");
	gui.view[GREEN_VPORT].drawarea= lookup_widget("drawingareag");
	gui.view[BLUE_VPORT].drawarea = lookup_widget("drawingareab");
	gui.view[RGB_VPORT].drawarea  = lookup_widget("drawingareargb");
	gui.view[MASK_VPORT].drawarea = lookup_widget("drawingareamask");
	/* GTK4: register the draw function for each viewport drawing area.  The
	 * old GTK3 "draw" signal connection (in siril.ui) is replaced by an
	 * explicit gtk_drawing_area_set_draw_func() per viewport. */
	install_drawarea_draw_funcs();
	gui.preview_area[0] = lookup_widget("drawingarea_reg_manual_preview1");
	gui.preview_area[1] = lookup_widget("drawingarea_reg_manual_preview2");
	memset(&gui.roi, 0, sizeof(roi_t)); // Clear the ROI
	initialize_image_display();
	init_mouse();

	/* GTK4: attach event controllers (click, motion, scroll) to each main
	 * viewport drawing area.  Replaces the GTK3 "scroll-event" /
	 * "button-press-event" / "motion-notify-event" / "enter-notify-event" /
	 * "leave-notify-event" signal wiring that used to live in siril.ui. */
	attach_drawingarea_event_controllers(gui.view[RED_VPORT].drawarea);
	attach_drawingarea_event_controllers(gui.view[GREEN_VPORT].drawarea);
	attach_drawingarea_event_controllers(gui.view[BLUE_VPORT].drawarea);
	attach_drawingarea_event_controllers(gui.view[RGB_VPORT].drawarea);
	attach_drawingarea_event_controllers(gui.view[MASK_VPORT].drawarea);

	/* set a transient */
	gtk_window_set_transient_for(GTK_WINDOW(lookup_widget("edge_dialog")), GTK_WINDOW(lookup_widget("control_window")));

	/* populate language combo */
	siril_language_fill_combo(com.pref.lang);

	/* Keybord Shortcuts */
	load_accels();

	/* Mouse event mask */
	set_mouse_event_mask();

	/* Wire up the menu-page actions for the Image Processing and Tools
	 * popover menus.  Each popover holds a GtkStack whose pages are the
	 * root menu and the various submenus; clicking a submenu trigger or
	 * back-button activates "<ns>.page" with the page name as parameter.
	 * (See uifiles/siril.ui after the Phase 8B/18 menu restructuring.) */
	wire_popover_menu_stack("processing_menu", "processing");
	wire_popover_menu_stack("tools_menu", "tools");

	/* Left-align the label of every plain-label GtkButton inside our
	 * menu popovers (siril.ui builds menus as GtkBox + GtkButton, so the
	 * labels default to xalign=0.5 which looks unlike a menu).  Submenu
	 * navigation buttons and check buttons are already aligned and
	 * skipped by the walker. */
	align_menu_popover_labels_left();

	/* Populate the Recent Files dropdown next to the Open button.  GTK4
	 * removed GtkRecentChooserMenu; we build a GMenuModel from
	 * GtkRecentManager (still present in GTK4) and attach it to the
	 * recent_menu_button, which was otherwise insensitive. */
	populate_recent_files_menu();

	/* Wire the GtkGestureClick that drives the icc_main_window_button.
	 * GTK3 used a button-press-event on a GtkButton with custom
	 * left/right dispatch; that signal doesn't exist in GTK4. */
	icc_main_window_button_attach_gesture();

	/* Wire GtkGestureClick "released" on scalemin / scalemax so the
	 * end-of-drag REMAP_ALL redraw fires (Phase 18 stripped the
	 * button-release-event binding from the .ui). */
	attach_display_scale_release_handlers();

	/* Select combo boxes that trigger some text display or other things */
	gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("comboboxstack_methods")), 0);

	GtkLabel *label_supported = GTK_LABEL(lookup_widget("label_supported_types"));
	gtk_label_set_text(label_supported, supported_files);

	adjust_sellabel();

	/* initialize theme */
	initialize_theme_GUI();

	/* map all actions for main window */
	siril_window_map_actions(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	siril_window_map_actions(GTK_APPLICATION_WINDOW(lookup_widget("seqlist_dialog")));

	/* secondary-click popover for undo/redo history navigation */
	setup_undo_redo_long_press();

	/* initialize menu gui */
	gui_function(update_MenuItem, NULL);

	/* initialize complex Mask tab label */
	initialize_mask_tab_label();
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget *page = gtk_notebook_get_nth_page(Color_Layers, MASK_VPORT);
	gtk_widget_set_visible(page, FALSE);


	/* initialize scripts and SPCC in threads:
	 * 1) initialize the scripts menu / SPCC widgets
	 * 2) sync the scripts repositories
	 * 3) update the scripts menu / SPCC widgets
	 */
	if (!is_online()) {
		siril_log_message(_("Siril started in offline mode. Will not attempt to update siril-scripts or siril-spcc-database...\n"));
	}
	g_thread_unref(g_thread_new("init_scripts_and_spcc", initialize_script_menu_and_spcc_widgets_serially, NULL));

	/* initialize command processor */
	init_command();

	/* initialize preprocessing */
	initialize_preprocessing();

	/* initialize registration methods */
	initialize_registration_methods();

	/* initialize stacking methods */
	initialize_stacking_methods();

	/* set focal and pixel pitch */
	init_astrometry();

	initialize_FITS_name_entries();

	initialize_log_tags();  // Can be called multiple times safely (idempotent)

	/* Initialize ROI settings */
	update_roi_config();

	/* register some callbacks */
	register_selection_update_callback(update_export_crop_label);
	register_selection_update_callback(update_display_selection);
	register_selection_update_callback(update_display_fwhm);

	/* Drag-and-drop targets are now installed per-widget by the GtkDropTarget
	 * sites in single_image.c, conversion.c and compositing.c.
	 * Sequence-list multi-selection lives on the GtkMultiSelection model in
	 * sequence_list.c.  The old GtkBuilder-based wiring that lived here is
	 * gone with the GtkTreeView removal. */

	g_signal_connect(lookup_widget("dialog_star_remix"), "close-request", G_CALLBACK(on_remix_close_clicked), NULL);
	g_signal_connect(lookup_widget("dialog_star_remix"), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	g_signal_connect(lookup_widget("histogram_dialog"), "close-request", G_CALLBACK(on_button_histo_close_clicked), NULL);
	g_signal_connect(lookup_widget("histogram_dialog"), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	/* Phase 11.11: pref_astro_tree_view + liststore_astrometry are gone.
	 * The astrometry catalogue panel is now built by annotations_pref.c. */

	siril_drag_single_image_set_dest();

	gui_function(set_GUI_CWD, NULL);
	siril_log_message(_("Default FITS extension is set to %s\n"), com.pref.ext);

	set_initial_display_mode((display_mode) com.pref.gui.default_rendering_mode);
	set_initial_histogram_display_mode(com.pref.gui.display_histogram_mode);

	update_spinCPU(GINT_TO_POINTER(com.max_thread));

	fill_astrometry_catalogue(com.pref.gui.catalog);
	init_GUI_from_settings();

	// init the Plot tab
	drawPlot();

	// init the cut structure
	initialize_cut_struct(&gui.cut);
	gui.cut.cut_measure = TRUE;
	gui.measure_start = (point){ -1., -1. };
	gui.measure_end = (point){ -1., -1. };

	/* every 0.5sec update memory display */
	g_timeout_add(500, update_displayed_memory, NULL);

	/* now that everything is loaded we can connect these signals
	 * Doing it in the glade file is a bad idea because they are called too many times during loading */
	/* GTK4: configure-event and window-state-event are gone.  Use GtkWindow
	 * property notifications.  notify::default-width and ::default-height
	 * fire when the user resizes the window; notify::maximized / ::fullscreened
	 * cover the state changes that window-state-event used to handle. */
	g_signal_connect(lookup_widget("control_window"), "notify::default-width",
	                 G_CALLBACK(on_control_window_configure_event), NULL);
	g_signal_connect(lookup_widget("control_window"), "notify::default-height",
	                 G_CALLBACK(on_control_window_configure_event), NULL);
	g_signal_connect(lookup_widget("control_window"), "notify::maximized",
	                 G_CALLBACK(on_control_window_window_state_event), NULL);
	g_signal_connect(lookup_widget("control_window"), "notify::fullscreened",
	                 G_CALLBACK(on_control_window_window_state_event), NULL);
	g_signal_connect(lookup_widget("main_panel"), "notify::position", G_CALLBACK(pane_notify_position_cb), NULL );

	/* GTK4 GtkExpander state-flag quirk: scan every expander in the
	 * builder (main UI and all dialog .ui files) and clear any stuck
	 * GTK_STATE_FLAG_INSENSITIVE bits in their child subtrees. */
	clear_stuck_insensitive_in_builder(gui.builder);

	gui_ready = TRUE;
}

/*****************************************************************************
 *      P U B L I C      C A L L B A C K      F U N C T I O N S              *
 ****************************************************************************/

void on_register_all_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_minscale_changed(GtkRange *range, gpointer user_data) {
	GtkEntry *minentry = (GtkEntry *)user_data;
	gchar buffer[12];
	int value = (int)gtk_range_get_value(range);
	if (value == gui.lo)
		return;
	gui.lo = value;
	g_sprintf(buffer, "%u", value);

	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(minentry), buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_maxscale_changed(GtkRange *range, gpointer user_data) {
	GtkEntry *maxentry = (GtkEntry *)user_data;
	gchar buffer[12];
	int value = (int)gtk_range_get_value(range);
	if (value == gui.hi)
		return;
	gui.hi = value;
	g_sprintf(buffer, "%u", value);

	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(maxentry), buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
}

/* End-of-drag redraw for the display sliders.  Phase 18 stripped the
 * GTK3 <signal name="button-release-event"/> binding, so the user-visible
 * effect was "sliders don't do anything".  GtkScale internally claims
 * its press/release events via its built-in GtkGestureClick, so adding a
 * sibling click gesture in BUBBLE phase never sees the release.  Using
 * GtkEventControllerLegacy in CAPTURE phase gets the GdkEvent before the
 * scale's internal gesture has a chance to consume it. */
static gboolean on_display_scale_legacy_event(GtkEventControllerLegacy *ctrl,
                                              GdkEvent *event, gpointer user_data) {
	(void) ctrl; (void) user_data;
	if (gdk_event_get_event_type(event) != GDK_BUTTON_RELEASE)
		return FALSE;   /* let the scale handle it normally */
	if (gui.sliders != USER) {
		gui.sliders = USER;
		sliders_mode_set_state(gui.sliders);
	}
	notify_gfit_data_modified();
	gfit_modified_update_gui();
	return FALSE;
}

void attach_display_scale_release_handlers(void) {
	const char *ids[] = { "scalemin", "scalemax" };
	for (size_t i = 0; i < G_N_ELEMENTS(ids); i++) {
		GtkWidget *scale = lookup_widget(ids[i]);
		if (!scale) continue;
		GtkEventController *legacy = gtk_event_controller_legacy_new();
		gtk_event_controller_set_propagation_phase(legacy, GTK_PHASE_CAPTURE);
		g_signal_connect(legacy, "event",
		                 G_CALLBACK(on_display_scale_legacy_event), NULL);
		gtk_widget_add_controller(scale, legacy);
	}
}

/* Legacy entry points retained so any external builder bindings still
 * resolve.  The actual release-handling is done via the legacy event
 * controller above. */
gboolean on_minscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	(void) widget; (void) event; (void) user_data;
	if (gui.sliders != USER) {
		gui.sliders = USER;
		sliders_mode_set_state(gui.sliders);
	}
	notify_gfit_data_modified();
	gfit_modified_update_gui();
	return FALSE;
}

gboolean on_maxscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	(void) widget; (void) event; (void) user_data;
	if (gui.sliders != USER) {
		gui.sliders = USER;
		sliders_mode_set_state(gui.sliders);
	}
	notify_gfit_data_modified();
	gfit_modified_update_gui();
	return FALSE;
}

/* a checkcut checkbox was toggled. Update the layer_info and others if chained. */
void on_checkcut_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	gui.cut_over = siril_toggle_get_active(GTK_WIDGET(togglebutton));
	notify_gfit_data_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}
/* gamut check was toggled. */
void on_gamutcheck_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	if (current_image_color_managed()) {
		lock_display_transform();
		if (com.gui_icc.proofing_transform)
			cmsDeleteTransform(com.gui_icc.proofing_transform);
		com.gui_icc.proofing_transform = initialize_proofing_transform();
		unlock_display_transform();
	}
	notify_gfit_data_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}

void on_cosmEnabledCheck_toggled(GtkCheckButton *button, gpointer user_data) {
	GtkWidget *CFA, *cc_frame_dark, *cc_frame_bad_pixel_map;
	gboolean is_active;

	CFA = lookup_widget("cosmCFACheck");
	cc_frame_dark = lookup_widget("cc_frame_dark");
	cc_frame_bad_pixel_map = lookup_widget("cc_frame_bad_pixel_map");

	is_active = siril_toggle_get_active(GTK_WIDGET(button));

	gtk_widget_set_sensitive(CFA, is_active);
	gtk_widget_set_sensitive(cc_frame_dark, is_active);
	gtk_widget_set_sensitive(cc_frame_bad_pixel_map, is_active);
}

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* focal_entry = gtk_editable_get_text(GTK_EDITABLE(editable));
	gfit->keywords.focal_length = g_ascii_strtod(focal_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchX_entry = gtk_editable_get_text(GTK_EDITABLE(editable));
	gfit->keywords.pixel_size_x = (float) g_ascii_strtod(pitchX_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchY_entry = gtk_editable_get_text(GTK_EDITABLE(editable));
	gfit->keywords.pixel_size_y = (float) g_ascii_strtod(pitchY_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_combobinning_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	gint index = gtk_drop_down_get_selected(box);

	switch (index) {
		case 0:
		case 1:
		case 2:
		case 3:
			gfit->keywords.binning_x = gfit->keywords.binning_y = (unsigned int) index + 1;
			break;
		case 4:
			gfit->keywords.binning_x = 1;
			gfit->keywords.binning_y = 2;
			break;
		case 5:
			gfit->keywords.binning_x = 1;
			gfit->keywords.binning_y = 3;
			break;
		case -1:
			gfit->keywords.binning_x = gfit->keywords.binning_y = 1;
			gtk_drop_down_set_selected(box, 0);
			break;
		default:
			fprintf(stderr, "Should not happen\n");
	}
	update_sampling_in_information();
	drawPlot();
}

void on_file_information_close_clicked(GtkButton *button, gpointer user_data) {
	GtkCheckButton *save_as_prefs_button = GTK_CHECK_BUTTON(lookup_widget("saveinfo_toggle"));
	if (siril_toggle_get_active(GTK_WIDGET(save_as_prefs_button))) {
		com.pref.starfinder_conf.focal_length = gfit->keywords.focal_length;
		com.pref.starfinder_conf.pixel_size_x = gfit->keywords.pixel_size_x;
		siril_log_message(_("Saved focal length %.2f and pixel size %.2f as default values\n"), gfit->keywords.focal_length, gfit->keywords.pixel_size_x);
	}
	gtk_widget_set_visible(lookup_widget("file_information"), FALSE);
	refresh_star_list();
}

void on_toggleButtonUnbinned_toggled(GtkCheckButton *button, gpointer user_data) {
	com.pref.binning_update = siril_toggle_get_active(GTK_WIDGET(button));
	update_sampling_in_information();
}

void on_button_clear_sample_clicked(GtkButton *button, gpointer user_data) {
	clear_sampling_setting_box();
}

static rectangle get_window_position(GtkWindow *window) {
	rectangle rec = { 0 };

	/* GTK4: gtk_window_get_position and gtk_window_get_size were removed —
	 * window position is owned by the compositor and cannot be queried.
	 * Use the default size as a best-effort fallback for w/h; x/y stay 0
	 * (the saved-state restore path centres the window when (x,y) == (0,0)). */
	(void)window;
	gtk_window_get_default_size(window, &rec.w, &rec.h);
	return rec;
}

gboolean save_main_window_state(gpointer user_data) {
	if (!com.script && com.pref.gui.remember_windows) {
		static GtkWindow *main_w = NULL;

		if (!main_w)
			main_w = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
		com.pref.gui.main_w_pos = get_window_position(main_w);
		com.pref.gui.is_maximized = gtk_window_is_maximized(main_w);
	}
	return FALSE;
}

gboolean load_main_window_state(gpointer user_data) {
	if (!com.script && com.pref.gui.remember_windows) {
		GtkWidget *win = lookup_widget("control_window");
		GdkRectangle workarea = { 0 };
		int scale_factor = 1;
		GdkSurface *surface = GTK_IS_NATIVE(win) ? gtk_native_get_surface(GTK_NATIVE(win)) : NULL;
		GdkMonitor *monitor = NULL;
		if (surface) {
			monitor = gdk_display_get_monitor_at_surface(gdk_surface_get_display(surface), surface);
		}
		if (monitor) {
			gdk_monitor_get_geometry(monitor, &workarea);
			scale_factor = gdk_monitor_get_scale_factor(monitor);
		}

		int w = com.pref.gui.main_w_pos.w;
		int h = com.pref.gui.main_w_pos.h;

		int x = CLAMP(com.pref.gui.main_w_pos.x, 0, workarea.width - w) * scale_factor;
		int y = CLAMP(com.pref.gui.main_w_pos.y, 0, workarea.height - h) * scale_factor;

		if (w > 0 && h > 0) {
			if (com.pref.gui.is_maximized) {
				gtk_window_maximize(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)));
			} else {
				/* Phase 17.5: gtk_window_move and gtk_window_resize are
				 * gone in GTK4 (compositor-managed placement / sizing).
				 * Default size is set instead; the WM is free to ignore
				 * the (x, y) request. */
				(void)x; (void)y;
				gtk_window_set_default_size(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), w, h);
			}
		}

		/* Now we handle the main panel */
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		GtkImage *image = GTK_IMAGE(gtk_button_get_child(GTK_BUTTON(lookup_widget("button_paned"))));
		GtkWidget *widget = gtk_paned_get_end_child(paned);

		gtk_widget_set_visible(widget, com.pref.gui.is_extended);
		if (com.pref.gui.is_extended) {
			gtk_image_set_from_icon_name(image, "pan-end-symbolic");
		} else {
			gtk_image_set_from_icon_name(image, "pan-start-symbolic");
		}
	}
	return FALSE;
}

/* Restore paned position after window is visible */
gboolean restore_paned_position(gpointer user_data) {
	if (!com.script && com.pref.gui.remember_windows && com.pref.gui.pan_position > 0) {
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		gtk_paned_set_position(paned, com.pref.gui.pan_position);
	}
	return FALSE;
}

void gtk_main_quit() {
	/* Pump the GTK main loop until any active processing job finishes or a
	 * timeout expires.  Cancellation is requested by siril_quit() before
	 * calling us.  Pumping lets any pending execute_idle_and_wait_for_it
	 * callbacks dispatch so the worker can unblock, see the cancel flag, and
	 * release gfit->rwlock before we attempt cleanup that also needs it
	 * (e.g. close_single_image → show_or_hide_mask_tab_idle).  Without this,
	 * the main thread deadlocks: worker holds the write lock waiting for an
	 * idle, main thread waits for the write lock. */
	gint64 deadline = g_get_monotonic_time() + 2 * G_TIME_SPAN_SECOND;
	siril_log_status(_("### Application Quit ###\n\nShutting down the processing subsystem..."));
	while (processing_is_job_active() && g_get_monotonic_time() < deadline) {
		if (g_main_context_pending(g_main_context_default()))
			g_main_context_iteration(g_main_context_default(), FALSE);
		else
			g_usleep(1000);
	}
	close_sequence(FALSE);	// save unfinished business
	close_single_image();	// close the previous image and free resources
	kill_child_process((GPid) -1, TRUE); // kill running child processes if any
	writeinitfile();		// save settings (like window positions)
	cmsUnregisterPlugins(); // unregister any lcms2 plugins
	g_slist_free_full(com.pref.gui.script_path, g_free);
	cleanup_common_profiles(); // close lcms2 data structures
	exit(EXIT_SUCCESS);
}

gboolean on_siril_window_delete(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	return !siril_quit();
}

// Returns TRUE if quitting should proceed, FALSE to cancel
gboolean siril_quit(void) {
	if (script_editor_has_unsaved_changes()) {
		gboolean quit_anyway = siril_confirm_dialog(
			_("Unsaved script changes"),
			_("There are unsaved changes in the script editor. Quit anyway?"),
			_("Quit")
		);
		if (!quit_anyway) {
			on_open_pythonpad(NULL, NULL);
			return FALSE; // Cancel quit
		}
	}

	if (!com.pref.gui.silent_quit) {
		gboolean quit = siril_confirm_dialog_and_remember(
			_("Closing application"),
			_("Are you sure you want to quit?"), _("Exit"),
			&com.pref.gui.silent_quit
		);
		if (!quit) {
			fprintf(stdout, "Staying on the application.\n");
			return FALSE; // Cancel quit
		}
	}

	/* Signal any running processing job to abort before the GTK main loop
	 * stops — jobs check processing_should_continue() in their frame loops. */
	processing_request_cancel();
	gtk_main_quit();
	return TRUE; // Proceed with quit
}

/* We give one signal event by toggle button to fix a bug. Without this solution
 * the toggled function was called 2 times
 * one call to select new button, one call to unselect previous one */
void on_radiobutton_minmax_toggled(GtkCheckButton *togglebutton,
		gpointer user_data) {
	WORD oldhi = gui.hi, oldlo = gui.lo;
	if (siril_toggle_get_active(GTK_WIDGET(togglebutton))) {
		gui.sliders = MINMAX;
		init_layers_hi_and_lo_values(gui.sliders);
		set_cutoff_sliders_values();
		if (gui.hi != oldhi || gui.lo != oldlo) {
			notify_gfit_data_modified();
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}
}

void on_radiobutton_hilo_toggled(GtkCheckButton *togglebutton,
		gpointer user_data) {
	WORD oldhi = gui.hi, oldlo = gui.lo;
	if (siril_toggle_get_active(GTK_WIDGET(togglebutton))) {
		gui.sliders = MIPSLOHI;
		init_layers_hi_and_lo_values(gui.sliders);
		set_cutoff_sliders_values();
		if (gui.hi != oldhi || gui.lo != oldlo) {
			notify_gfit_data_modified();
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}
}

void on_radiobutton_user_toggled(GtkCheckButton *togglebutton,
		gpointer user_data) {
	if (siril_toggle_get_active(GTK_WIDGET(togglebutton))) {
		gui.sliders = USER;
		init_layers_hi_and_lo_values(gui.sliders);
		set_cutoff_sliders_values();
// No redraw needed here as the slider values don't actually change until the user changes the values
	}
}

// Not a callback as such but this is called from multiple stretch
// functions in different files, so here is a good place to put it

void setup_stretch_sliders() {
	gboolean changed = FALSE;
	GtkCheckButton *button = GTK_CHECK_BUTTON(lookup_widget("radiobutton_user"));
	if (!siril_toggle_get_active(GTK_WIDGET(button))) {
		siril_toggle_set_active(GTK_WIDGET(button), TRUE);
		gui.sliders = USER;
	}
	if (gui.hi != USHRT_MAX) {
		gui.hi = USHRT_MAX;
		changed = TRUE;
	}
	if (gui.lo != 0) {
		gui.lo = 0;
		changed = TRUE;
	}
	if (changed) {
		set_cutoff_sliders_values();
		notify_gfit_data_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_max_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

		if (gui.sliders != USER) {
			gui.sliders = USER;
			sliders_mode_set_state(gui.sliders);
		}
		gui.hi = value;

		set_cutoff_sliders_values();

		notify_gfit_data_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

gboolean on_max_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_editable_set_text(GTK_EDITABLE(widget), "65535");
	return FALSE;
}

gboolean on_min_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_editable_set_text(GTK_EDITABLE(widget), "0");
	return FALSE;
}

void on_min_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

		if (gui.sliders != USER) {
			gui.sliders = USER;
			sliders_mode_set_state(gui.sliders);
		}
		gui.lo = value;

		set_cutoff_sliders_values();

		notify_gfit_data_modified();
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_seqproc_entry_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *widget = GTK_DROP_DOWN(obj);
	(void)pspec;
	gchar *name = siril_drop_down_get_active_text(GTK_DROP_DOWN(widget));
	if (name && name[0] != '\0') {
		gchar *type;

		set_cursor_waiting(TRUE);
		const char *ext = get_filename_ext(name);
		if (!strcmp(ext, "ser")) {
			name[strlen(name) - 1] = 'q';
			type = " SER";
#ifdef HAVE_FFMS2
		} else if (!check_for_film_extensions(ext)) {
			int len = strlen(ext);
			strncpy(name + strlen(name) - len - 1, "seq", len + 1);
			type = " AVI";
#endif
		} else
			type = "";
		gchar *msg = g_strdup_printf(_("Selected %s sequence %s..."), type, name);
		set_progress_bar_data(msg, PROGRESS_DONE);
		set_seq(name);
		set_cursor_waiting(FALSE);
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
		g_free(msg);
	}
	g_free(name);
	gui_function(launch_clipboard_survey, NULL);
}

/* signal handler for the gray window layer change */
void on_notebook1_switch_page(GtkNotebook *notebook, GtkWidget *page,
		guint page_num, gpointer user_data) {
	gui.cvport = page_num;
	redraw(REDRAW_OVERLAY);
	update_display_selection();	// update the dimensions of the selection when switching page
	if (page_num < MASK_VPORT)
		update_display_fwhm();
}

struct checkSeq_filter_data {
	int retvalue;
};

static gboolean end_checkSeq(gpointer p) {
	struct checkSeq_filter_data *args = (struct checkSeq_filter_data *) p;
	end_generic(NULL);
	if (!args->retvalue)
		update_sequences_list(NULL);
	else control_window_switch_to_tab(OUTPUT_LOGS);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	free(args);
	return FALSE;
}

static gpointer checkSeq(gpointer p) {
	struct checkSeq_filter_data *args = (struct checkSeq_filter_data *) p;
	args->retvalue = check_seq();
	siril_add_idle(end_checkSeq, args);
	return GINT_TO_POINTER(0);
}

void on_checkseqbutton_clicked(GtkButton *button, gpointer user_data) {
	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	set_cursor_waiting(TRUE);
	set_progress_bar_data(_("Searching for sequences in "
				"the current working directory..."), PROGRESS_NONE);

	struct checkSeq_filter_data *args = calloc(1, sizeof(struct checkSeq_filter_data));

	args->retvalue = 0;
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(checkSeq, args)) {
		free(args);
	}
}

void on_button_data_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_set_visible(lookup_widget("data_dialog"), FALSE);
}

void on_comboboxreglayer_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *widget = GTK_DROP_DOWN(obj);
	(void)pspec;
	if (gtk_drop_down_get_selected(widget) == -1)
		return;
	free_reference_image();
	update_stack_interface(TRUE);
	update_reg_interface(FALSE);
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = gtk_spin_button_get_value_as_int(spinbutton);
}

void on_menu_rgb_align_select(GtkWidget *menuitem, gpointer user_data) {
	gboolean sel_is_drawn = ((com.selection.w > 0.0) && (com.selection.h > 0.0));

	gtk_widget_set_sensitive(lookup_widget("rgb_align_dft"), sel_is_drawn);
	gtk_widget_set_sensitive(lookup_widget("rgb_align_psf"), sel_is_drawn);
}

void on_gotoStacking_button_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void on_checkbutton_auto_evaluate_toggled(GtkCheckButton *button,
		gpointer user_data) {
	GtkWidget *entry = (GtkWidget *)user_data;
	gtk_widget_set_sensitive(entry, !siril_toggle_get_active(GTK_WIDGET(button)));
}

void on_clean_sequence_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean cleanreg = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("seq_clean_reg"))));
	gboolean cleanstat = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("seq_clean_stat"))));
	gboolean cleansel = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("seq_clean_sel"))));

	if ((cleanreg || cleanstat || cleansel) && sequence_is_loaded()) {
		GString *warning = g_string_new(_("This erases the following data, and there's no possible undo:\n"));
		if (cleanreg) warning = g_string_append(warning, _("\n- Registration"));
		if (cleanstat) warning = g_string_append(warning, _("\n- Statistics"));
		if (cleansel) warning = g_string_append(warning, _("\n- Selection"));

		gchar *str = g_string_free(warning, FALSE);

		gboolean clear = siril_confirm_dialog(_("Clear Sequence Data?"), str, _("Clear Data"));
		g_free(str);

		if (clear) {
			clean_sequence(&com.seq, cleanreg, cleanstat, cleansel);
			update_stack_interface(TRUE);
			update_reg_interface(FALSE);
			adjust_sellabel();
			reset_plot();
			set_layers_for_registration();
			if (cleanstat)
				gfit_modified_update_gui();
			siril_message_dialog(GTK_MESSAGE_INFO, _("Sequence"), _("The requested data of the sequence has been cleaned."));
		}
	}
}

void on_clean_gaia_cache_clicked(GtkButton *button, gpointer user_data) {
	gchar *msg = g_strdup_printf(_("You are about to clean your Gaia online cache to remove any "
				"entries that have not been accessed for %d days. "
				"This operation cannot be undone."), com.pref.astrometry.gaia_cache_duration);

	int confirm = siril_confirm_dialog(_("Cache Cleaner"), msg, _("Proceed"));
	g_free(msg);
	if (!confirm) {
		return;
	}

	flux_cache_purge(com.pref.astrometry.gaia_cache_duration);
}

void on_purge_user_catalogue_clicked(GtkButton *button, gpointer user_data) {
	const gchar *caller = gtk_buildable_get_buildable_id(GTK_BUILDABLE(button));
	gboolean is_sso = g_strrstr(caller, "SSO") != NULL;
	gchar *msg = g_strdup_printf(_("You are about to purge your %s user catalog. This means the "
				"file containing the manually added objects will be deleted. "
				"This operation cannot be undone."), (is_sso) ? _("SOLAR SYSTEM") : _("DEEP SKY"));

	int confirm = siril_confirm_dialog(_("Catalog deletion"), msg, _("Purge Catalogue"));
	g_free(msg);
	if (!confirm) {
		return;
	}
	siril_cat_index cat_index = (is_sso) ? CAT_AN_USER_SSO : CAT_AN_USER_DSO;
	gchar *filename = get_annotation_catalog_filename(cat_index, FALSE);
	GFile *file = g_file_new_for_path(filename);

	g_autoptr(GError) local_error = NULL;
	if (g_file_query_exists(file, NULL) && !g_file_delete(file, NULL, &local_error)) {
		g_warning("Failed to delete %s: %s", g_file_peek_path(file), local_error->message);
	} else {
		gboolean was_displayed = com.found_object != NULL;
		purge_user_catalogue(cat_index);
		if (was_displayed) {
			redraw(REDRAW_OVERLAY);
		}
	}
	g_free(filename);
	g_object_unref(file);
}

/* Per-row GObject for the running-tasks list.  Owns the child_info pointer
 * (does not free it: kill_running_processes() owns the child list).  */
#define SIRIL_TYPE_CHILD_TASK_ROW (siril_child_task_row_get_type())
G_DECLARE_FINAL_TYPE(SirilChildTaskRow, siril_child_task_row, SIRIL, CHILD_TASK_ROW, GObject)

struct _SirilChildTaskRow {
	GObject parent_instance;
	child_info *child;
	gchar      *name;
	gchar      *time_str;
};
G_DEFINE_TYPE(SirilChildTaskRow, siril_child_task_row, G_TYPE_OBJECT)

static void siril_child_task_row_finalize(GObject *obj) {
	SirilChildTaskRow *self = SIRIL_CHILD_TASK_ROW(obj);
	g_clear_pointer(&self->name,     g_free);
	g_clear_pointer(&self->time_str, g_free);
	G_OBJECT_CLASS(siril_child_task_row_parent_class)->finalize(obj);
}
static void siril_child_task_row_class_init(SirilChildTaskRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_child_task_row_finalize;
}
static void siril_child_task_row_init(SirilChildTaskRow *self) { (void)self; }

static void child_task_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void child_task_bind_name_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilChildTaskRow *row = SIRIL_CHILD_TASK_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row->name);
}
static void child_task_bind_time_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilChildTaskRow *row = SIRIL_CHILD_TASK_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row->time_str);
}

/* Phase 14G.4: GtkWindow surrogate for the synchronous-style
 * child-process selection dialog.  GTK4 removed gtk_dialog_run; we run
 * our own nested GMainLoop and let the explicit Cancel/OK button
 * "clicked" handlers feed the result to it. */
struct child_dlg_state { gint response; GMainLoop *loop; };


static void child_dlg_cancel_clicked(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct child_dlg_state *st = (struct child_dlg_state *)ud;
	st->response = GTK_RESPONSE_CANCEL;
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
}
static void child_dlg_ok_clicked(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct child_dlg_state *st = (struct child_dlg_state *)ud;
	st->response = GTK_RESPONSE_OK;
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
}
static gboolean child_dlg_close_request(GtkWindow *w, gpointer ud) {
	(void)w;
	struct child_dlg_state *st = (struct child_dlg_state *)ud;
	st->response = GTK_RESPONSE_CANCEL;
	if (st->loop && g_main_loop_is_running(st->loop))
		g_main_loop_quit(st->loop);
	return FALSE;
}

GPid show_child_process_selection_dialog(GSList *children) {
	/* Phase 14G.4: GtkDialog → GtkWindow.  Sync semantics preserved via
	 * a nested GMainLoop fed by per-button "clicked" handlers. */
	GtkWidget *dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Select task to stop"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dialog),
		GTK_WINDOW(lookup_widget("control_window")));

	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 400);
	/* GTK4: GTK_WIN_POS_CENTER is gone — the WM centres modal/transient
	 * windows automatically when a transient parent is set. */

	GtkWidget *scrolled_window = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window),
			GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_has_frame(GTK_SCROLLED_WINDOW(scrolled_window), TRUE);

	/* Build the GtkColumnView model */
	GListStore *store = g_list_store_new(SIRIL_TYPE_CHILD_TASK_ROW);
	child_mutex_lock();
	for (GSList *l = children; l != NULL; l = l->next) {
		child_info *child = (child_info *)l->data;
		gchar *time_str = g_date_time_format(child->datetime, "%H:%M:%S");
		SirilChildTaskRow *row = g_object_new(SIRIL_TYPE_CHILD_TASK_ROW, NULL);
		row->child    = child;
		row->name     = g_strdup(child->name);
		row->time_str = time_str;  /* takes ownership */
		g_list_store_append(store, row);
		g_object_unref(row);
	}
	child_mutex_unlock();

	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(store)));
	GtkColumnView *cv = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));

	GtkSignalListItemFactory *fname = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fname, "setup", G_CALLBACK(child_task_setup_cb),     NULL);
	g_signal_connect(fname, "bind",  G_CALLBACK(child_task_bind_name_cb), NULL);
	GtkColumnViewColumn *cname = gtk_column_view_column_new(_("Process Name"), GTK_LIST_ITEM_FACTORY(fname));
	gtk_column_view_column_set_expand(cname, TRUE);
	gtk_column_view_append_column(cv, cname);
	g_object_unref(cname);

	GtkSignalListItemFactory *ftime = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(ftime, "setup", G_CALLBACK(child_task_setup_cb),     NULL);
	g_signal_connect(ftime, "bind",  G_CALLBACK(child_task_bind_time_cb), NULL);
	GtkColumnViewColumn *ctime = gtk_column_view_column_new(_("Start Time"), GTK_LIST_ITEM_FACTORY(ftime));
	gtk_column_view_column_set_expand(ctime, FALSE);
	gtk_column_view_append_column(cv, ctime);
	g_object_unref(ctime);

	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(cv));

	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_widget_set_hexpand(scrolled_window, TRUE);
	gtk_widget_set_vexpand(scrolled_window, TRUE);
	gtk_box_append(GTK_BOX(content_area), scrolled_window);

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *btn_cancel = gtk_button_new_with_mnemonic(_("_Cancel"));
	GtkWidget *btn_ok     = gtk_button_new_with_mnemonic(_("_Stop Process"));
	gtk_widget_add_css_class(btn_ok, "destructive-action");
	gtk_box_append(GTK_BOX(bbox), btn_cancel);
	gtk_box_append(GTK_BOX(bbox), btn_ok);
	gtk_box_append(GTK_BOX(content_area), bbox);

	gtk_window_set_child(GTK_WINDOW(dialog), content_area);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), btn_ok);

	struct child_dlg_state state = { GTK_RESPONSE_NONE, g_main_loop_new(NULL, FALSE) };
	g_signal_connect(btn_cancel, "clicked",        G_CALLBACK(child_dlg_cancel_clicked), &state);
	g_signal_connect(btn_ok,     "clicked",        G_CALLBACK(child_dlg_ok_clicked),     &state);
	g_signal_connect(dialog,     "close-request",  G_CALLBACK(child_dlg_close_request),  &state);
	gtk_widget_set_visible(dialog, TRUE);
	g_main_loop_run(state.loop);
	g_main_loop_unref(state.loop);

	GPid selected_pid = (GPid) -1;
	if (state.response == GTK_RESPONSE_OK) {
		guint pos = gtk_single_selection_get_selected(sel);
		if (pos != GTK_INVALID_LIST_POSITION) {
			SirilChildTaskRow *row = SIRIL_CHILD_TASK_ROW(g_list_model_get_item(G_LIST_MODEL(store), pos));
			if (row && row->child)
				selected_pid = row->child->childpid;
			g_clear_object(&row);
		}
	}

	g_object_unref(sel);
	g_object_unref(store);
	gtk_window_destroy(GTK_WINDOW(dialog));

	return selected_pid;
}

gboolean set_seq_browser_active(gpointer user_data) {
	gboolean state = (gboolean) GPOINTER_TO_INT(user_data);
	if (!state) { // Deactivating: check if the dialog is already open, if so close it
		if (gtk_widget_is_visible(lookup_widget("seqlist_dialog"))) {
			siril_close_dialog("seqlist_dialog");
		}
	}
	GtkWidget *widget = lookup_widget("seqlist_button");
	gtk_widget_set_sensitive(widget, state);
	return FALSE;
}

int seq_qphot(sequence *seq, int layer) {
	framing_mode framing = REGISTERED_FRAME;
	if (framing == REGISTERED_FRAME && !seq->regparam[layer])
		framing = ORIGINAL_FRAME;
	if (framing == ORIGINAL_FRAME) {
		GtkCheckButton *follow = GTK_CHECK_BUTTON(lookup_widget("followStarCheckButton"));
		if (siril_toggle_get_active(GTK_WIDGET(follow)))
			framing = FOLLOW_STAR_FRAME;
	}
	siril_log_message(_("Running the PSF on the sequence, layer %d\n"), layer);
	return seqpsf(seq, layer, FALSE, TRUE, FALSE, framing, TRUE, com.script);
}

gboolean in_gtk_thread(void) {
    return g_main_context_is_owner(g_main_context_default());
}

gboolean ensure_seqlist_dialog_closed_idle(gpointer user_data) {
	if (gtk_widget_is_visible(lookup_widget("seqlist_dialog"))) {
		siril_close_dialog("seqlist_dialog");
	}
	return FALSE;
}

void ensure_seqlist_dialog_closed() {
	if (com.headless) return; // nothing to do
	if (in_gtk_thread())
		ensure_seqlist_dialog_closed_idle(NULL);
	else
		execute_idle_and_wait_for_it(ensure_seqlist_dialog_closed_idle, NULL);
}
