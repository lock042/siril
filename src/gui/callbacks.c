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

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/siril_app_dirs.h"
#include "core/siril_language.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "filters/graxpert.h"
#include "gui/cut.h"
#include "gui/graxpert.h"
#include "gui/keywords_tree.h"
#include "gui/registration.h"
#include "gui/photometric_cc.h"
#include "gui/stacking.h"
#include "algos/siril_wcs.h"
#include "algos/star_finder.h"
#include "io/annotation_catalogues.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/siril_pythonmodule.h"
#include "annotations_pref.h"
#include "compositing/align_rgb.h"
#include "preferences.h"
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
#include "fix_xtrans_af.h"
#include "siril_intro.h"
#include "siril_preview.h"
#include "siril-window.h"
#include "registration_preview.h"

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
	gtk_widget_show(w_image);
	gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget(button_name)), w_image);
	gtk_widget_show(lookup_widget(button_name));

	g_free(image);
}

void handle_owner_change(GtkClipboard *clipboard, GdkEvent *event, gpointer data) {
	/*Only surveys the name of the opened item vs the clipboard content and change the color accoringly*/

	GtkLabel *label_name_of_seq = NULL;
	const char *format_green = "<span foreground=\"green\">%s</span>";
	const char *format_white = "<span foreground=\"white\">%s</span>";
	char *markup;

	label_name_of_seq = GTK_LABEL(lookup_widget("label_name_of_seq"));

	/* Get the clipboard object */
	clipboard = gtk_clipboard_get (GDK_SELECTION_CLIPBOARD);
	/* Get the clipboard content */
	char *clipboard_content = gtk_clipboard_wait_for_text(clipboard);
	/* Initialize for non-NULL */
	if (clipboard_content == NULL) {clipboard_content = "NoRules";}

	/* Set the right color*/
	if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		if ((strcmp(filename, clipboard_content) == 0)) {
			markup = g_markup_printf_escaped (format_green, _("Image: "));
			gtk_label_set_markup(label_name_of_seq, markup);
		} else {
			markup = g_markup_printf_escaped (format_white, _("Image: "));
			gtk_label_set_markup(label_name_of_seq, markup);
		}
		g_free(filename);
		filename = NULL;
	}


	if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);
		if ((strcmp(seq_basename, clipboard_content) == 0)) {
			markup = g_markup_printf_escaped (format_green, _("Sequence: "));
			gtk_label_set_markup(label_name_of_seq, markup);
		} else {
			markup = g_markup_printf_escaped (format_white, _("Sequence: "));
			gtk_label_set_markup(label_name_of_seq, markup);
		}
		g_free(seq_basename);
	}

}

gboolean launch_clipboard_survey(gpointer user_data) {
	if (com.script)
		return FALSE;
	GtkClipboard *clipboard = NULL;

	/* Get the clipboard object */
	clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);

	/* To launch the Handle*/
	g_signal_connect(clipboard, "owner-change", G_CALLBACK(handle_owner_change), NULL);
	return FALSE;
}

void on_press_seq_field() {
	GtkClipboard *clipboard = NULL;

	/* Get the clipboard object */
	clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);

	/* For an Image */
	if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		/* Set clipboard text */
		gtk_clipboard_set_text(clipboard, filename, -1);
		g_free(filename);
		/* For a Sequence */
	} else if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);
		/* Set clipboard text */
		gtk_clipboard_set_text(clipboard, seq_basename, -1);
		g_free(seq_basename);
	}
}

static void update_icons_to_theme(gboolean is_dark) {
	siril_debug_print("Loading %s theme...\n", is_dark ? "dark" : "light");
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
	}
}

void siril_set_theme(int active) {
	GtkSettings *settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", active == 0, NULL);
	update_icons_to_theme(active == 0);
	update_icons_sequence_list(active == 0);
}

void on_combo_theme_changed(GtkComboBox *box, gpointer user_data) {
	int active = gtk_combo_box_get_active(box);

	siril_set_theme(active);
}

int populate_roi() {
	if (com.python_command)
		return 1;
	if (gui.roi.selection.w == 0 || gui.roi.selection.h == 0)
		return 1;
	int retval = 0;
	size_t npixels_roi = gui.roi.selection.w * gui.roi.selection.h;
	size_t npixels_gfit = gfit.rx * gfit.ry;
	size_t nchans = gfit.naxes[2];
	g_assert(nchans ==1 || nchans == 3);
	gboolean rgb = (nchans == 3);
	clearfits(&gui.roi.fit);
	copyfits(&gfit, &gui.roi.fit, CP_FORMAT, -1);
	gui.roi.fit.rx = gui.roi.fit.naxes[0] = gui.roi.selection.w;
	gui.roi.fit.ry = gui.roi.fit.naxes[1] = gui.roi.selection.h;
	gui.roi.fit.naxes[2] = nchans;
	gui.roi.fit.naxis = nchans == 1 ? 2 : 3;
	if (gui.roi.fit.type == DATA_FLOAT) {
		gui.roi.fit.fdata = malloc(npixels_roi * nchans * sizeof(float));
		if (!gui.roi.fit.fdata)
			retval = 1;
		gui.roi.fit.fpdata[0] = gui.roi.fit.fdata;
		gui.roi.fit.fpdata[1] = rgb? gui.roi.fit.fdata + npixels_roi : gui.roi.fit.fdata;
		gui.roi.fit.fpdata[2] = rgb? gui.roi.fit.fdata + 2 * npixels_roi : gui.roi.fit.fdata;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				float *srcindex = gfit.fdata + (npixels_gfit * c) + ((gfit.ry - y - (gui.roi.selection.y + 1)) * gfit.rx) + gui.roi.selection.x;
				float *destindex = gui.roi.fit.fdata + (npixels_roi * c) + (gui.roi.fit.rx * y);
				memcpy(destindex, srcindex, (gui.roi.selection.w) * sizeof(float));
			}
		}
	} else {
		gui.roi.fit.data = malloc(npixels_roi * nchans * sizeof(WORD));
		if (!gui.roi.fit.data)
			retval = 1;
		gui.roi.fit.pdata[0] = gui.roi.fit.data;
		gui.roi.fit.pdata[1] = rgb? gui.roi.fit.data + npixels_roi : gui.roi.fit.data;
		gui.roi.fit.pdata[2] = rgb? gui.roi.fit.data + 2 * npixels_roi : gui.roi.fit.data;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				WORD *srcindex = gfit.data + (npixels_gfit * c) + ((gfit.ry - y - (gui.roi.selection.y + 1)) * gfit.rx) + gui.roi.selection.x;
				WORD *destindex = gui.roi.fit.data + (npixels_roi * c) + y * gui.roi.fit.rx;
				memcpy(destindex, srcindex, gui.roi.selection.w * sizeof(WORD));
			}
		}
	}
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
		GtkWidget *check = NULL;
		GtkWindow *parent = siril_get_active_window();
		if (!GTK_IS_WINDOW(parent))
			parent = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
		GtkWidget *dialog = gtk_message_dialog_new(parent, GTK_DIALOG_MODAL, GTK_MESSAGE_INFO, GTK_BUTTONS_OK, _("Region of Interest supported"));
		gtk_message_dialog_format_secondary_text(GTK_MESSAGE_DIALOG(dialog),
					_("This image operation supports ROI processing. This is "
					"indicated by the green border around the ROI. In image "
					"operation dialogs that do not support ROI processing, and "
					"when no dialog is open, the ROI will have a red outline.\n\n"
					"While a ROI is active, processing will only preview the ROI. "
					"When the Apply button is clicked, the operation will apply "
					"to the whole image."));
		check = gtk_check_button_new_with_mnemonic(_("_Do not show this dialog again"));
		gtk_box_pack_end(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG (dialog))),
				check, FALSE, FALSE, 0);
		gtk_widget_set_halign(check, GTK_ALIGN_START);
		gtk_widget_set_margin_start(check, 6);
		gtk_widget_show(check);
		gtk_dialog_run(GTK_DIALOG(dialog));
		com.pref.gui.enable_roi_warning = !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check));
		gtk_widget_destroy(dialog);
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
	if (gui.roi.selection.w == 0 && gui.roi.selection.h == 0)
		return GINT_TO_POINTER(0);
	g_mutex_lock(&roi_mutex); // Wait until any thread previews are finished
	if (com.python_command) {
//		on_clear_roi();
		g_mutex_unlock(&roi_mutex);
		return GINT_TO_POINTER(0);
	}
	cancel_pending_update();
	if (gui.roi.operation_supports_roi && com.pref.gui.enable_roi_warning)
		roi_info_message_if_needed();
	// Ensure any pending ROI changes are overwritten by the backup
	// Must copy the whole backup to gfit and gui.roi.fit to account
	// for switching between full image and ROI
	memcpy(&gui.roi.selection, &com.selection, sizeof(rectangle));
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
	GtkComboBox *box;

	box = GTK_COMBO_BOX(lookup_widget("combo_theme"));

	g_signal_handlers_block_by_func(box, on_combo_theme_changed, NULL);
	gtk_combo_box_set_active(box, com.pref.gui.combo_theme);
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

	gfit.keywords.hi = gtk_adjustment_get_value(adj1);
	gfit.keywords.lo = gtk_adjustment_get_value(adj2);
}

/* Sets maximum value for contrast scales. Minimum is always 0.
 * Should be done on first image load of a sequence and when single images are loaded.
 * Max value is taken from gfit.maxi, recomputed if not present.
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

	max_val = (gfit.type == DATA_FLOAT ? USHRT_MAX_DOUBLE : (gfit.orig_bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE));
	siril_debug_print(_("Setting MAX value for cutoff sliders adjustments (%f)\n"), max_val);
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
	static GtkToggleButton *cutmax = NULL;
	if (adjmin == NULL) {
		adjmax = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemax")); // scalemax
		adjmin = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustmentscalemin")); // scalemin
		maxentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "max_entry"));
		minentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "min_entry"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(gui.builder, "checkcut_max"));
	}

	siril_debug_print(_("Setting ranges scalemin=%d, scalemax=%d\n"), gui.lo, gui.hi);
	gtk_adjustment_set_value(adjmin, (gdouble)gui.lo);
	gtk_adjustment_set_value(adjmax, (gdouble)gui.hi);
	g_snprintf(buffer, 6, "%u", gui.hi);
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_snprintf(buffer, 6, "%u", gui.lo);
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	gtk_toggle_button_set_active(cutmax, gui.cut_over);
}

void on_display_item_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	if (!gtk_check_menu_item_get_active(checkmenuitem)) return;

	static GtkLabel *label_display_menu = NULL;
	if (!label_display_menu) {
		label_display_menu = GTK_LABEL(lookup_widget("display_button_name"));
	}

	gui.rendering_mode = get_display_mode_from_menu();
	siril_debug_print("Display mode %d\n", gui.rendering_mode);
	gboolean override_label = FALSE;

	if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit.type != DATA_FLOAT)
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
		gtk_label_set_text(label_display_menu, gtk_menu_item_get_label(GTK_MENU_ITEM(checkmenuitem)));

	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	siril_window_autostretch_actions(app_win, gui.rendering_mode == STF_DISPLAY && gfit.naxes[2] == 3);

	gui.icc.same_primaries = same_primaries(gfit.icc_profile, gui.icc.monitor, gui.icc.soft_proof ? gui.icc.soft_proof : NULL);

	if (single_image_is_loaded() || sequence_is_loaded()) {
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_autohd_item_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	gui.use_hd_remap = gtk_check_menu_item_get_active(menuitem);
	if (gui.rendering_mode == STF_DISPLAY) {
		if (gui.use_hd_remap) {
			if (gfit.type == DATA_FLOAT)
				allocate_hd_remap_indices();
			siril_log_message(_("The AutoStretch display mode will use a %d bit LUT\n"), (int) log2(gui.hd_remap_max));
		} else {
			hd_remap_indices_cleanup();
			siril_log_message(_("The AutoStretch display mode will use a 16 bit LUT\n"));
		}
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_button_apply_hd_bitdepth_clicked(GtkSpinButton *button, gpointer user_data) {
	int bitdepth = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_hd_bitdepth")));
	siril_debug_print("bitdepth: %d\n", bitdepth);
	if (gui.hd_remap_max != 1 << bitdepth) {
		siril_log_message(_("Setting HD AutoStretch display mode bit depth to %d...\n"), bitdepth);
//		set_cursor_waiting(TRUE);

		com.pref.hd_bitdepth = bitdepth;
		gui.hd_remap_max = 1 << bitdepth;
		if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit.type == DATA_FLOAT) {
			allocate_hd_remap_indices();
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
//		set_cursor_waiting(FALSE);
	}
}


/* Sets the display mode combo box to the value stored in the relevant struct.
 * The operation is purely graphical. */
void set_display_mode() {
	GtkCheckMenuItem *button;
	static GtkMenu *display_menu = NULL;
	static GtkLabel *label_display_menu = NULL;
	if (!display_menu) {
		display_menu = GTK_MENU(lookup_widget("menu_display"));
		label_display_menu = GTK_LABEL(lookup_widget("display_button_name"));
	}

	switch (gui.rendering_mode) {
	default:
	case LINEAR_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[LINEAR_DISPLAY]));
		break;
	case LOG_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[LOG_DISPLAY]));
	break;
	case SQRT_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[SQRT_DISPLAY]));
	break;
	case SQUARED_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[SQUARED_DISPLAY]));
		break;
	case ASINH_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[ASINH_DISPLAY]));
		break;
	case STF_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[STF_DISPLAY]));
	break;
	case HISTEQ_DISPLAY:
		button = GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[HISTEQ_DISPLAY]));
		break;
	}
	g_signal_handlers_block_by_func(button, on_display_item_toggled, NULL);
	gtk_check_menu_item_set_active(button, TRUE);
	g_signal_handlers_unblock_by_func(button, on_display_item_toggled, NULL);

	gtk_label_set_text(label_display_menu, gtk_menu_item_get_label(GTK_MENU_ITEM(button)));

}

void set_unlink_channels(gboolean unlinked) {
	siril_debug_print("channels unlinked: %d\n", unlinked);
	gui.unlink_channels = unlinked;
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
	is_a_singleRGB_image_loaded = isrgb(&gfit) && single_image_is_loaded();
	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	/* some toolbar buttons */
	gtk_widget_set_sensitive(lookup_widget("toolbarbox"), any_image_is_loaded);

	gboolean wcs_ok = any_image_is_loaded && has_wcs(&gfit);
	GAction *action_annotate = g_action_map_lookup_action(G_ACTION_MAP(app_win), "annotate-object");
	GAction *action_grid = g_action_map_lookup_action(G_ACTION_MAP(app_win), "wcs-grid");

	/* any image with wcs information is needed */
	siril_window_enable_wcs_proc_actions(app_win, wcs_ok);
	siril_window_enable_wcs_disto_proc_actions(app_win, wcs_ok && gfit.keywords.wcslib->lin.dispre);

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

	set_undo_redo_tooltip();

	/* save and save as */
	GAction *action_save = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save");
	GAction *action_save_as = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save-as");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save), is_a_single_image_loaded && com.uniq->fileexist);
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save_as), any_image_is_loaded);

	/* popup menu */

	/* fits header */
	GAction *action_fits_header = g_action_map_lookup_action (G_ACTION_MAP(app_win), "fits-header");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_fits_header), any_image_is_loaded && gfit.header != NULL);
	/* annotate dialog */
	GAction *action_annotate_dialog = g_action_map_lookup_action (G_ACTION_MAP(app_win), "annotate-dialog");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_annotate_dialog), any_image_is_loaded && has_wcs(&gfit));
	/* Lightcurve process */
	GAction *action_nina_light_curve = g_action_map_lookup_action (G_ACTION_MAP(app_win), "nina_light_curve");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_nina_light_curve), sequence_is_loaded() && has_wcs(&gfit));
	/* selection is needed */
	siril_window_enable_if_selection_actions(app_win, com.selection.w && com.selection.h);
	/* selection and sequence is needed */
	siril_window_enable_if_selection_sequence_actions(app_win, com.selection.w && com.selection.h && sequence_is_loaded());
	/* selectoin and isrgb is needed */
	siril_window_enable_if_selection_rgb_actions(app_win, com.selection.w && com.selection.h && isrgb(&gfit));

	/* Image processing Menu */

	/* single RGB image is needed */
	siril_window_enable_rgb_proc_actions(app_win, is_a_singleRGB_image_loaded);
	/* single RGB image with wcs information is needed */
	siril_window_enable_rgb_wcs_proc_actions(app_win, is_a_singleRGB_image_loaded && has_wcs(&gfit));
	/* single or sequence RGB is needed */
	siril_window_enable_any_rgb_proc_actions(app_win, is_a_singleRGB_image_loaded|| sequence_is_loaded());
	/* any image is needed */
	siril_window_enable_any_proc_actions(app_win, any_image_is_loaded);
	/* any mono image is needed */
	siril_window_enable_any_mono_proc_actions(app_win, any_image_is_loaded && !isrgb(&gfit));
	/* only single image needed */
	siril_window_enable_single_proc_actions(app_win, is_a_single_image_loaded);
	/* no images needed */
	siril_window_enable_none_proc_actions(app_win, TRUE);

	/* Image information menu */
	siril_window_enable_image_actions(app_win, any_image_is_loaded);

	/* auto-stretch actions */
	siril_window_autostretch_actions(app_win, gui.rendering_mode == STF_DISPLAY && gfit.naxes[2] == 3);

	/* keywords list */
	if (gtk_widget_is_visible(lookup_widget("keywords_dialog")))
		refresh_keywords_dialog(NULL);
	return FALSE;
}

void sliders_mode_set_state(sliders_mode sliders) {
	GtkToggleButton *radiobutton;	// Must not be static
	gchar *str[] = { "radiobutton_hilo", "radiobutton_minmax", "radiobutton_user" };
	void *func[] = { on_radiobutton_hilo_toggled, on_radiobutton_minmax_toggled, on_radiobutton_user_toggled };

	radiobutton = GTK_TOGGLE_BUTTON(lookup_widget(str[sliders]));

	g_signal_handlers_block_by_func(radiobutton, func[sliders], NULL);
	gtk_toggle_button_set_active(radiobutton, TRUE);
	g_signal_handlers_unblock_by_func(radiobutton, func[sliders], NULL);
}

display_mode get_display_mode_from_menu() {
	static GtkWidget *menu = NULL;
	if (!menu)
		menu = lookup_widget("menu_display");
	for (int i = 0; i < G_N_ELEMENTS(display_item_name); i++) {
		gboolean is_active = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(lookup_widget(display_item_name[i])));
		if (is_active) return (display_mode) i;
	}
	return LINEAR_DISPLAY;
}

void set_initial_display_mode(display_mode x) {
	switch(x) {
		case LINEAR_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("linear_item")), TRUE);
			break;
		case LOG_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("log_item")), TRUE);
			break;
		case SQRT_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("square_root_item")), TRUE);
			break;
		case SQUARED_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("squared_item")), TRUE);
			break;
		case ASINH_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("asinh_item")), TRUE);
			break;
		case STF_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("auto_item")), TRUE);
			break;
		case HISTEQ_DISPLAY:
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("histo_item")), TRUE);
			break;
	}
}

static void set_initial_histogram_display_mode(int x) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckLogButton"));
	g_signal_handlers_block_by_func(button, on_histo_toggled, NULL);
	switch(x) {
	case LINEAR_DISPLAY:
		gtk_toggle_button_set_active(button, FALSE);
		break;
	case LOG_DISPLAY:
		gtk_toggle_button_set_active(button, TRUE);
	}
	g_signal_handlers_unblock_by_func(button, on_histo_toggled, NULL);
}


void update_prepro_interface(gboolean allow_debayer) {
	static GtkToggleButton *udark = NULL, *uoffset = NULL, *uflat = NULL,
			       *checkAutoEvaluate = NULL;
	static GtkWidget *prepro_button = NULL, *cosme_grid = NULL, *dark_optim = NULL;
       	static GtkWidget *equalize = NULL, *auto_eval = NULL, *flat_norm = NULL;
       	static GtkWidget *debayer = NULL, *fix_xtrans = NULL;
	static GtkComboBox *output_type = NULL;
	if (udark == NULL) {
		udark = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
		uoffset = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));
		uflat = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
		checkAutoEvaluate = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto_evaluate"));
		output_type = GTK_COMBO_BOX(lookup_widget("prepro_output_type_combo"));
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
			(sequence_is_loaded() || single_image_is_loaded())
			&& (gtk_toggle_button_get_active(udark)
				|| gtk_toggle_button_get_active(uoffset)
				|| gtk_toggle_button_get_active(uflat)));
	gtk_widget_set_sensitive(cosme_grid, gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(dark_optim, gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(equalize, gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(auto_eval, gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(flat_norm,
			gtk_toggle_button_get_active(uflat) &&
			!gtk_toggle_button_get_active(checkAutoEvaluate));

	gtk_widget_set_sensitive(debayer, allow_debayer && gtk_widget_get_sensitive(prepro_button));
	gtk_widget_set_sensitive(fix_xtrans, gtk_toggle_button_get_active(udark) || gtk_toggle_button_get_active(uoffset));

	gtk_widget_set_sensitive(GTK_WIDGET(output_type), sequence_is_loaded());
	int type = com.seq.type;
	if (com.seq.type < 0 || com.seq.type > SEQ_FITSEQ)
		type = SEQ_REGULAR;
	gtk_combo_box_set_active(output_type, type);
}

void clear_sampling_setting_box() {
	GtkComboBox *binning = GTK_COMBO_BOX(
			gtk_builder_get_object(gui.builder, "combobinning"));
	GtkEntry* focal_entry = GTK_ENTRY(lookup_widget("focal_entry"));
	GtkEntry* pitchX_entry = GTK_ENTRY(lookup_widget("pitchX_entry"));
	GtkEntry* pitchY_entry = GTK_ENTRY(lookup_widget("pitchY_entry"));

	gtk_entry_set_text(focal_entry, "");
	gtk_entry_set_text(pitchX_entry, "");
	gtk_entry_set_text(pitchY_entry, "");
	gtk_combo_box_set_active(binning, 0);
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
	if (gfit.naxes[2] == 1) {
		if (layer != 0) {
			siril_debug_print("unloaded layer name requested %d\n", layer);
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
			siril_debug_print("unloaded layer name requested %d\n", layer);
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
	return -1;
}

void update_display_selection() {
	static const gchar *label_selection[] = { "labelselection_red", "labelselection_green", "labelselection_blue", "labelselection_rgb" };
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
	static const gchar *label_fwhm[] = { "labelfwhm_red", "labelfwhm_green", "labelfwhm_blue", "labelfwhm_rgb" };
	static gchar fwhm_buffer[256] = { 0 };

	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		g_sprintf(fwhm_buffer, " ");
	} else if (com.selection.w && com.selection.h) {// Now we don't care about the size of the sample. Minimization checks that
		if (com.selection.w < 300 && com.selection.h < 300 && com.selection.w > 5 && com.selection.h > 5) {
			double roundness;
			double fwhm_val = psf_get_fwhm(&gfit, select_vport(gui.cvport), &com.selection, &roundness);
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
	if (single_image_is_loaded()) {	// unique image
		filename = com.uniq->filename;
		orig_filename = g_file_read_link(filename, &error);
		nb_layers = com.uniq->nb_layers;
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

	if (local_filename)
		free(filename);
	g_free(base_name);
	g_free(orig_filename);
	g_free(orig_base_name);
	if (error) g_error_free(error);
	if (concat_base_name) g_string_free(concat_base_name, TRUE);
}

void on_precision_item_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	if (!single_image_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Cannot convert a sequence file"),
				_("A sequence file cannot be converted to 32 bits. This operation can only be done on a single file."));
	} else {
		size_t ndata = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];

		if (gfit.type == DATA_FLOAT) {
			gboolean convert = siril_confirm_dialog(_("Precision loss"),
					_("Converting the image from 32 bits to 16 bits may lead to a loss of numerical accuracy. "
							"Getting back to 32 bits will not recover this loss.\n"
							"Are you sure you want to convert your data?"), _("Convert to 16 bits"));
			if (convert) {
				if (is_preview_active()) {
					fits *backup_gfit = get_preview_gfit_backup();
					fit_replace_buffer(backup_gfit, float_buffer_to_ushort(backup_gfit->fdata, ndata), DATA_USHORT);
				}
				fit_replace_buffer(&gfit, float_buffer_to_ushort(gfit.fdata, ndata), DATA_USHORT);
				invalidate_gfit_histogram();
				update_gfit_histogram_if_needed();
				redraw(REMAP_ALL);
			}
		} else if (gfit.type == DATA_USHORT) {
			if (is_preview_active()) {
				fits *backup_gfit = get_preview_gfit_backup();
				fit_replace_buffer(backup_gfit, ushort_buffer_to_float(backup_gfit->data, ndata), DATA_FLOAT);
			}
			fit_replace_buffer(&gfit, ushort_buffer_to_float(gfit.data, ndata), DATA_FLOAT);
			invalidate_gfit_histogram();
			update_gfit_histogram_if_needed();
			redraw(REMAP_ALL);
		}
	}
	gui_function(set_precision_switch, NULL);
}

gboolean set_precision_switch(gpointer user_data) {
	if (!com.script) {
		GtkLabel *label = GTK_LABEL(lookup_widget("precision_button_name"));
		GtkCheckMenuItem *float_button = GTK_CHECK_MENU_ITEM(lookup_widget("32bits_item"));
		GtkCheckMenuItem *ushort_button = GTK_CHECK_MENU_ITEM(lookup_widget("16bits_item"));

		gtk_label_set_text(label, gfit.type == DATA_USHORT ? _("16 bits") : _("32 bits"));
		g_signal_handlers_block_by_func(float_button, on_precision_item_toggled, NULL);
		gtk_check_menu_item_set_active(float_button, gfit.type == DATA_FLOAT);
		gtk_check_menu_item_set_active(ushort_button, gfit.type == DATA_USHORT);
		g_signal_handlers_unblock_by_func(float_button,	on_precision_item_toggled, NULL);
	}
	return FALSE;
}

/* updates the combo box of registration layers to reflect data availability */
int set_layers_for_registration() {
	static GtkComboBoxText *cbbt_layers = NULL;
	int i;
	int reminder;

	if (cbbt_layers == NULL)
		cbbt_layers = GTK_COMBO_BOX_TEXT(
				gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	reminder = gtk_combo_box_get_active(GTK_COMBO_BOX(cbbt_layers));

	gtk_combo_box_text_remove_all(cbbt_layers);
	for (i = 0; i < com.seq.nb_layers; i++) {
		gchar *layer;
		const gchar *layer_name = layer_name_for_gfit(i);
		layer = g_strdup_printf("%d: %s", i, layer_name);
		if (com.seq.regparam && com.seq.reglayer > -1) {
			str_append(&layer,  " (*)");
			if (reminder == -1 || ((reminder >= 0) && (com.seq.reglayer != reminder))) // set as default selection
				reminder = i;
		}
		gtk_combo_box_text_append_text(cbbt_layers, layer);

		g_free(layer);
	}
	/* First initialization */
	if (reminder == -1) {
		if (com.seq.nb_layers == 3)
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 1);
		else
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 0);
	}
	/* Already initialized or default selection to channel with data */
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), reminder);

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

	gtk_widget_show(lookup_widget("data_dialog"));
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

	const char *t_bias = gtk_entry_get_text(mbias);
	if (!str[0] && t_bias[0] == '\0') {
		str[0] = g_strdup_printf("master-bias%s", com.pref.ext);
	}
	if (str[0]) {
		gtk_entry_set_text(mbias, str[0]);
	}

	const char *t_dark = gtk_entry_get_text(mdark);
	if (!str[1] && t_dark[0] == '\0') {
		str[1] = g_strdup_printf("master-dark%s", com.pref.ext);
	}
	if (str[1]) {
		gtk_entry_set_text(mdark, str[1]);
	}

	const char *t_flat = gtk_entry_get_text(mflat);
	if (!str[2] && t_flat[0] == '\0') {
		str[2] = g_strdup_printf("master-flat%s", com.pref.ext);
	}
	if (str[2]) {
		gtk_entry_set_text(mflat, str[2]);
	}

	const char *t_stack = gtk_entry_get_text(final_stack);
	if (!str[3] && t_stack[0] == '\0') {
		if (sequence_is_loaded())
			str[3] = g_strdup_printf("%s_stacked%s", com.seq.seqname, com.pref.ext);
		else
			str[3] = g_strdup_printf("stack_result%s", com.pref.ext);
	}
	if (str[3]) {
		gtk_entry_set_text(final_stack, str[3]);
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
			gtk_entry_set_text(output_file, com.pref.prepro.stack_default);
			return;
		}
	}
	const gchar* com_ext = get_com_ext(com.seq.fz);

	msg = g_strdup_printf("%s%sstacked%s", com.seq.seqname,
			g_str_has_suffix(com.seq.seqname, "_") ?
			"" : (g_str_has_suffix(com.seq.seqname, "-") ? "" : "_"), com_ext);
	gtk_entry_set_text(output_file, msg);

	g_free(msg);
}

gboolean close_tab(gpointer user_data) {
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget* page;

	if (com.seq.nb_layers == 1 || gfit.naxes[2] == 1) {
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("B&W"));
	} else {
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("Red"));
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_show(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_show(page);
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_show(page);
	}
	return FALSE;
}

void activate_tab(int vport) {
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook1"));
	if (gtk_notebook_get_current_page(notebook) != vport)
		gtk_notebook_set_current_page(notebook, vport);
	// com.cvport is set in the event handler for changed page
}

gboolean init_right_tab(gpointer user_data) {
	activate_tab(isrgb(&gfit) ? RGB_VPORT : RED_VPORT);
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
		"win.pcc-processing",         "<Primary><Shift>p", NULL,
		"win.spcc-processing",        "<Primary><Shift>c", NULL,
		"win.compstars",              "<Primary><Shift>b", NULL,
		"win.nina_light_curve",       "<Primary><Shift>n", NULL,
		"win.pickstar",               "<Primary>space", NULL,
		"win.dyn-psf",                "<Primary>F6", NULL,
		"win.clipboard",              "<Primary><Shift>x", NULL,

		"win.negative-processing",    "<Primary>i", NULL,
		"win.rotation90-processing",  "<Primary>Right", NULL,
		"win.rotation270-processing", "<Primary>Left", NULL,
		"win.mirrorx-processing",     "<Primary>Up", NULL,
		"win.mirrory-processing",     "<Primary>Down", NULL,

		"win.regframe",               "<Primary>r", NULL,

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
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton")), com.pref.debayer.open_debayer);
}

gboolean set_GUI_CWD(gpointer user_data) {
	if (!com.wd)
		return FALSE;
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("headerbar"));

	gchar *str = g_strdup_printf("Siril-%s", VERSION);
	gtk_header_bar_set_title(bar , str);

	gchar *truncated_wd = siril_truncate_str(com.wd, 50);
	gtk_header_bar_set_subtitle(bar, truncated_wd);

	g_free(str);
	g_free(truncated_wd);
	return FALSE;
}

static void initialize_preprocessing() {
	GtkToggleButton *cfaButton, *eqButton, *xtransButton;

	cfaButton = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	gtk_toggle_button_set_active(cfaButton, com.pref.prepro.cfa);
	eqButton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	gtk_toggle_button_set_active(eqButton, com.pref.prepro.equalize_cfa);
	xtransButton = GTK_TOGGLE_BUTTON(lookup_widget("fix_xtrans_af"));
	gtk_toggle_button_set_active(xtransButton, com.pref.prepro.fix_xtrans);

	update_prepro_interface(FALSE);
}

void update_sampling_in_information() {
	GtkLabel *sampling_label = GTK_LABEL(lookup_widget("info_sampling_label"));
	// display sampling based on information displayed in the window
	// otherwise, we could use get_wcs_image_resolution(fits *fit)
	if (gfit.keywords.focal_length > 0.0 && gfit.keywords.pixel_size_x > 0.0) {
		float pixel_size = gfit.keywords.pixel_size_x;
		if (gfit.keywords.binning_x > 1 && com.pref.binning_update)
			pixel_size *= gfit.keywords.binning_x;
		double scale = get_resolution(gfit.keywords.focal_length, pixel_size);
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
	double fl = gfit.keywords.focal_length > 0.0 ? gfit.keywords.focal_length : com.pref.starfinder_conf.focal_length;
	sprintf(buf, "%.3f", fl);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("focal_entry")), buf);

	double pixsz = gfit.keywords.pixel_size_x > 0.0 ? gfit.keywords.pixel_size_x : com.pref.starfinder_conf.pixel_size_x;
	sprintf(buf, "%.2f", pixsz);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchX_entry")), buf);

	pixsz = gfit.keywords.pixel_size_y > 0.0 ? gfit.keywords.pixel_size_y : com.pref.starfinder_conf.pixel_size_x;
	sprintf(buf, "%.2f", pixsz);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchY_entry")), buf);

	GtkComboBox *binning = GTK_COMBO_BOX(lookup_widget("combobinning"));
	if (!gfit.keywords.binning_x || !gfit.keywords.binning_y) {
		gtk_combo_box_set_active(binning, 0);
	}
	else if (gfit.keywords.binning_x == gfit.keywords.binning_y) {
		/* squared binning */
		gtk_combo_box_set_active(binning, (gint) gfit.keywords.binning_x - 1);
	} else {
		short coeff =
			gfit.keywords.binning_x > gfit.keywords.binning_y ?
			gfit.keywords.binning_x / gfit.keywords.binning_y :
			gfit.keywords.binning_y / gfit.keywords.binning_x;
		switch (coeff) {
			case 2:
				gtk_combo_box_set_active(binning, 4);
				break;
			case 3:
				gtk_combo_box_set_active(binning, 5);
				break;
			default:
				siril_log_message(_("Binning %d x %d is not supported\n"),
						gfit.keywords.binning_x, gfit.keywords.binning_y);
		}
	}

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggleButtonUnbinned")), com.pref.binning_update);
	//gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("saveinfo_toggle")), com.pref.astrometry.update_default_scale);
}

static GtkTargetEntry drop_types[] = {
	{ "text/uri-list", 0, 0 }
};

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

static void pane_notify_position_cb(GtkPaned *paned, gpointer user_data) {
	static gboolean first_resize = TRUE;
	if (first_resize) {
		if (com.pref.gui.remember_windows && com.pref.gui.pan_position > 0) {
			gtk_paned_set_position(paned, com.pref.gui.pan_position);
		}
		first_resize = FALSE;
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

void initialize_all_GUI(gchar *supported_files) {
	/* pre-check the Gaia archive status */
	check_gaia_archive_status();
	/* initializing internal structures with widgets (drawing areas) */
	gui.view[RED_VPORT].drawarea  = lookup_widget("drawingarear");
	gui.view[GREEN_VPORT].drawarea= lookup_widget("drawingareag");
	gui.view[BLUE_VPORT].drawarea = lookup_widget("drawingareab");
	gui.view[RGB_VPORT].drawarea  = lookup_widget("drawingareargb");
	gui.preview_area[0] = lookup_widget("drawingarea_reg_manual_preview1");
	gui.preview_area[1] = lookup_widget("drawingarea_reg_manual_preview2");
	memset(&gui.roi, 0, sizeof(roi_t)); // Clear the ROI
	initialize_image_display();
	init_mouse();

	/* set a transient */
	gtk_window_set_transient_for(GTK_WINDOW(lookup_widget("edge_dialog")), GTK_WINDOW(lookup_widget("control_window")));

	/* populate language combo */
	siril_language_fill_combo(com.pref.lang);

	/* Keybord Shortcuts */
	load_accels();

	/* Mouse event mask */
	set_mouse_event_mask();

	/* Select combo boxes that trigger some text display or other things */
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboboxstack_methods")), 0);

	GtkLabel *label_supported = GTK_LABEL(lookup_widget("label_supported_types"));
	gtk_label_set_text(label_supported, supported_files);

	adjust_sellabel();

	/* initialize theme */
	initialize_theme_GUI();

	/* map all actions for main window */
	siril_window_map_actions(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	siril_window_map_actions(GTK_APPLICATION_WINDOW(lookup_widget("seqlist_dialog")));

	/* initialize menu gui */
	gui_function(update_MenuItem, NULL);
	initialize_script_menu(TRUE);

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

	initialize_log_tags();

	/* Initialize ROI settings */
	update_roi_config();

	/* register some callbacks */
	register_selection_update_callback(update_export_crop_label);
	register_selection_update_callback(update_display_selection);
	register_selection_update_callback(update_display_fwhm);

	/* support for converting files by dragging onto the GtkTreeView */
	gtk_drag_dest_set(lookup_widget("treeview_convert"),
			GTK_DEST_DEFAULT_MOTION | GTK_DEST_DEFAULT_HIGHLIGHT, drop_types, G_N_ELEMENTS(drop_types),
			GDK_ACTION_COPY);

	/* Due to a glade bug, this property is often removed, lets code it */
	GtkTreeSelection *selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "treeview_selection_convert"));
	gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);
	g_signal_connect(selection, "changed", G_CALLBACK(on_treeview_selection_convert_changed), NULL);

	/* Due to another bug in glade we write these callbacks here
	 * In glade there are not sorted as we want */
	g_signal_connect(lookup_widget("dialog_star_remix"), "delete-event", G_CALLBACK(on_remix_close_clicked), NULL);
	g_signal_connect(lookup_widget("dialog_star_remix"), "delete-event", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	g_signal_connect(lookup_widget("histogram_dialog"), "delete-event", G_CALLBACK(on_button_histo_close_clicked), NULL);
	g_signal_connect(lookup_widget("histogram_dialog"), "delete-event", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "treeview-selection"));
	gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);

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

	/* every 0.5sec update memory display */
	g_timeout_add(500, update_displayed_memory, NULL);

#ifndef HAVE_LIBCURL
	// SPCC is not available if compiled without networking
	gtk_widget_set_visible(lookup_widget("proc_spcc"), FALSE);
	// Remove it from the RGB composition color calibration methods list too
	gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(lookup_widget("rgbcomp_cc_method")), 2);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("rgbcomp_cc_method")), 1);
#endif

	/* now that everything is loaded we can connect these signals
	 * Doing it in the glade file is a bad idea because they are called too many times during loading */
	g_signal_connect(lookup_widget("control_window"), "configure-event", G_CALLBACK(on_control_window_configure_event), NULL);
	g_signal_connect(lookup_widget("control_window"), "window-state-event", G_CALLBACK(on_control_window_window_state_event), NULL);
	g_signal_connect(lookup_widget("main_panel"), "notify::position", G_CALLBACK(pane_notify_position_cb), NULL );
	/* populate SPCC combos in a thread */
	g_thread_unref(g_thread_new("spcc_combos", populate_spcc_combos_async, NULL));
	/* GraXpert checks, if required */
	g_thread_unref(g_thread_new("graxpert_checks", graxpert_setup_async, NULL));
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
	gtk_entry_set_text(minentry, buffer);
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
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
}

gboolean on_minscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (gui.sliders != USER) {
		gui.sliders = USER;
		sliders_mode_set_state(gui.sliders);
	}
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	return FALSE;
}

gboolean on_maxscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (gui.sliders != USER) {
		gui.sliders = USER;
		sliders_mode_set_state(gui.sliders);
	}
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	return FALSE;
}

/* a checkcut checkbox was toggled. Update the layer_info and others if chained. */
void on_checkcut_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	gui.cut_over = gtk_toggle_button_get_active(togglebutton);
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}
/* gamut check was toggled. */
void on_gamutcheck_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (gfit.color_managed) {
		lock_display_transform();
		if (gui.icc.proofing_transform)
			cmsDeleteTransform(gui.icc.proofing_transform);
		gui.icc.proofing_transform = initialize_proofing_transform();
		unlock_display_transform();
	}
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
}

void on_cosmEnabledCheck_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *CFA, *cc_frame_dark, *cc_frame_bad_pixel_map;
	gboolean is_active;

	CFA = lookup_widget("cosmCFACheck");
	cc_frame_dark = lookup_widget("cc_frame_dark");
	cc_frame_bad_pixel_map = lookup_widget("cc_frame_bad_pixel_map");

	is_active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(CFA, is_active);
	gtk_widget_set_sensitive(cc_frame_dark, is_active);
	gtk_widget_set_sensitive(cc_frame_bad_pixel_map, is_active);
}

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* focal_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.keywords.focal_length = g_ascii_strtod(focal_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchX_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.keywords.pixel_size_x = (float) g_ascii_strtod(pitchX_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchY_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.keywords.pixel_size_y = (float) g_ascii_strtod(pitchY_entry, NULL);
	update_sampling_in_information();
	drawPlot();
}

void on_combobinning_changed(GtkComboBox *box, gpointer user_data) {
	gint index = gtk_combo_box_get_active(box);

	switch (index) {
		case 0:
		case 1:
		case 2:
		case 3:
			gfit.keywords.binning_x = gfit.keywords.binning_y = (unsigned int) index + 1;
			break;
		case 4:
			gfit.keywords.binning_x = 1;
			gfit.keywords.binning_y = 2;
			break;
		case 5:
			gfit.keywords.binning_x = 1;
			gfit.keywords.binning_y = 3;
			break;
		case -1:
			gfit.keywords.binning_x = gfit.keywords.binning_y = 1;
			gtk_combo_box_set_active(box, 0);
			break;
		default:
			fprintf(stderr, "Should not happen\n");
	}
	update_sampling_in_information();
	drawPlot();
}

void on_file_information_close_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *save_as_prefs_button = GTK_TOGGLE_BUTTON(lookup_widget("saveinfo_toggle"));
	if (gtk_toggle_button_get_active(save_as_prefs_button)) {
		com.pref.starfinder_conf.focal_length = gfit.keywords.focal_length;
		com.pref.starfinder_conf.pixel_size_x = gfit.keywords.pixel_size_x;
		siril_log_message(_("Saved focal length %.2f and pixel size %.2f as default values\n"), gfit.keywords.focal_length, gfit.keywords.pixel_size_x);
	}
	gtk_widget_hide(lookup_widget("file_information"));
	refresh_star_list();
}

void on_toggleButtonUnbinned_toggled(GtkToggleButton *button, gpointer user_data) {
	com.pref.binning_update = gtk_toggle_button_get_active(button);
	update_sampling_in_information();
}

void on_button_clear_sample_clicked(GtkButton *button, gpointer user_data) {
	clear_sampling_setting_box();
}

static rectangle get_window_position(GtkWindow *window) {
	gint x, y, w, h;
	rectangle rec = { 0 };

	gtk_window_get_position(window, &x, &y);
	gtk_window_get_size(window, &w, &h);
	rec.x = x;
	rec.y = y;
	rec.w = w;
	rec.h = h;
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
		GdkWindow *window = gtk_widget_get_window(win);
		GdkMonitor *monitor = gdk_display_get_monitor_at_window(gdk_window_get_display(window), window);

		gdk_monitor_get_workarea(monitor, &workarea);
		int scale_factor = gdk_monitor_get_scale_factor(monitor);

		int w = com.pref.gui.main_w_pos.w;
		int h = com.pref.gui.main_w_pos.h;

		int x = CLAMP(com.pref.gui.main_w_pos.x, 0, workarea.width - w) * scale_factor;
		int y = CLAMP(com.pref.gui.main_w_pos.y, 0, workarea.height - h) * scale_factor;

		if (w > 0 && h > 0) {
			if (com.pref.gui.is_maximized) {
				gtk_window_maximize(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)));
			} else {
				gtk_window_move(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), x, y);
				gtk_window_resize(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), w, h);
			}
		}

		/* Now we handle the main panel */
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
		GtkWidget *widget = gtk_paned_get_child2(paned);

		gtk_widget_set_visible(widget, com.pref.gui.is_extended);
		if (com.pref.gui.is_extended) {
			gtk_image_set_from_icon_name(image, "pan-end-symbolic", GTK_ICON_SIZE_BUTTON);
		} else {
			gtk_image_set_from_icon_name(image, "pan-start-symbolic", GTK_ICON_SIZE_BUTTON);
		}
	}
	return FALSE;
}

void gtk_main_quit() {
	writeinitfile();		// save settings (like window positions)
	close_sequence(FALSE);	// save unfinished business
	close_single_image();	// close the previous image and free resources
	kill_child_process((GPid) -1, TRUE); // kill running child processes if any
	cmsUnregisterPlugins(); // unregister any lcms2 plugins
	g_slist_free_full(com.pref.gui.script_path, g_free);
	cleanup_common_profiles(); // close lcms2 data structures
	exit(EXIT_SUCCESS);
}

void siril_quit() {
	if (com.pref.gui.silent_quit) {
		gtk_main_quit();
	}
	gboolean quit = siril_confirm_dialog_and_remember(_("Closing application"),
			_("Are you sure you want to quit?"), _("Exit"), &com.pref.gui.silent_quit);
	if (quit) {
		gtk_main_quit();
	} else {
		fprintf(stdout, "Staying on the application.\n");
	}
}

/* We give one signal event by toggle button to fix a bug. Without this solution
 * the toggled function was called 2 times
 * one call to select new button, one call to unselect previous one */
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	WORD oldhi = gui.hi, oldlo = gui.lo;
	if (gtk_toggle_button_get_active(togglebutton)) {
		gui.sliders = MINMAX;
		init_layers_hi_and_lo_values(gui.sliders);
		set_cutoff_sliders_values();
		if (gui.hi != oldhi || gui.lo != oldlo) {
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}
}

void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	WORD oldhi = gui.hi, oldlo = gui.lo;
	if (gtk_toggle_button_get_active(togglebutton)) {
		gui.sliders = MIPSLOHI;
		init_layers_hi_and_lo_values(gui.sliders);
		set_cutoff_sliders_values();
		if (gui.hi != oldhi || gui.lo != oldlo) {
			redraw(REMAP_ALL);
			gui_function(redraw_previews, NULL);
		}
	}
}

void on_radiobutton_user_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton)) {
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
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_user"));
	if (!gtk_toggle_button_get_active(button)) {
		gtk_toggle_button_set_active(button, TRUE);
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
		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_max_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

		if (gui.sliders != USER) {
			gui.sliders = USER;
			sliders_mode_set_state(gui.sliders);
		}
		gui.hi = value;

		set_cutoff_sliders_values();

		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

gboolean on_max_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "65535");
	return FALSE;
}

gboolean on_min_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "0");
	return FALSE;
}

void on_min_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

		if (gui.sliders != USER) {
			gui.sliders = USER;
			sliders_mode_set_state(gui.sliders);
		}
		gui.lo = value;

		set_cutoff_sliders_values();

		redraw(REMAP_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void on_seqproc_entry_changed(GtkComboBox *widget, gpointer user_data) {
	gchar *name = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(widget));
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
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	set_cursor_waiting(TRUE);
	set_progress_bar_data(_("Searching for sequences in "
				"the current working directory..."), PROGRESS_NONE);

	struct checkSeq_filter_data *args = malloc(sizeof(struct checkSeq_filter_data));

	args->retvalue = 0;
	set_cursor_waiting(TRUE);
	start_in_new_thread(checkSeq, args);
}

void on_button_data_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("data_dialog"));
}

void on_comboboxreglayer_changed(GtkComboBox *widget, gpointer user_data) {
	if (gtk_combo_box_get_active(widget) == -1)
		return;
	free_reference_image();
	update_stack_interface(TRUE);
	update_reg_interface(FALSE);
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = gtk_spin_button_get_value_as_int(spinbutton);
}

void on_menu_rgb_align_select(GtkMenuItem *menuitem, gpointer user_data) {
	gboolean sel_is_drawn = ((com.selection.w > 0.0) && (com.selection.h > 0.0));

	gtk_widget_set_sensitive(lookup_widget("rgb_align_dft"), sel_is_drawn);
	gtk_widget_set_sensitive(lookup_widget("rgb_align_psf"), sel_is_drawn);
}

void on_gotoStacking_button_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void on_checkbutton_auto_evaluate_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *entry = (GtkWidget *)user_data;
	gtk_widget_set_sensitive(entry, !gtk_toggle_button_get_active(button));
}

void on_clean_sequence_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean cleanreg = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("seq_clean_reg")));
	gboolean cleanstat = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("seq_clean_stat")));
	gboolean cleansel = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("seq_clean_sel")));

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
				notify_gfit_modified();
			siril_message_dialog(GTK_MESSAGE_INFO, _("Sequence"), _("The requested data of the sequence has been cleaned."));
		}
	}
}

void on_purge_user_catalogue_clicked(GtkButton *button, gpointer user_data) {
	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE(button));
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

GPid show_child_process_selection_dialog(GSList *children) {
	// Create the dialog
	GtkWidget *dialog = gtk_dialog_new_with_buttons(
		"Select task to stop",
		NULL,  // No parent window
		GTK_DIALOG_MODAL,
		"_Cancel", GTK_RESPONSE_CANCEL,
		"_Stop Process", GTK_RESPONSE_OK,
		NULL
	);
	gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_OK);

	// Set minimum dialog size
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 400);
	gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_CENTER);

	// Create a scrolled window to hold the list
	GtkWidget *scrolled_window = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(
		GTK_SCROLLED_WINDOW(scrolled_window),
		GTK_POLICY_AUTOMATIC,  // Horizontal scrollbar as needed
		GTK_POLICY_AUTOMATIC   // Vertical scrollbar as needed
	);
	gtk_scrolled_window_set_shadow_type(
		GTK_SCROLLED_WINDOW(scrolled_window),
		GTK_SHADOW_IN  // Add a border around the scrolled window
	);

	// Create a list store and tree view
	GtkListStore *list_store = gtk_list_store_new(
		3,  // 3 columns
		G_TYPE_POINTER,  // Pointer to child_info
		G_TYPE_STRING,   // Process name
		G_TYPE_STRING    // Formatted start time
	);

	GtkWidget *tree_view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(list_store));

	// Enable multiple selection
	GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree_view));
	gtk_tree_selection_set_mode(selection, GTK_SELECTION_SINGLE);

	// Add tree view to the scrolled window
	gtk_container_add(GTK_CONTAINER(scrolled_window), tree_view);

	// Add columns to the tree view
	GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
	GtkTreeViewColumn *column;

	// Process Name Column
	column = gtk_tree_view_column_new_with_attributes(
		"Process Name",
		renderer,
		"text", 1,
		NULL
	);
	gtk_tree_view_column_set_expand(column, TRUE);  // Allow column to expand
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Start Time Column
	column = gtk_tree_view_column_new_with_attributes(
		"Start Time",
		renderer,
		"text", 2,
		NULL
	);
	gtk_tree_view_column_set_expand(column, FALSE);  // Keep start time column compact
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Populate the list store
	GtkTreeIter iter;
	for (GSList *l = children; l != NULL; l = l->next) {
		child_info *child = (child_info *)l->data;

		// Format time to 24-hour HH:MM
		gchar *time_str = g_date_time_format(child->datetime, "%H:%M:%S");

		gtk_list_store_append(list_store, &iter);
		gtk_list_store_set(list_store, &iter,
			0, child,         // Pointer to child_info
			1, child->name,   // Process name
			2, time_str,      // Formatted start time
			-1
		);

		// Free the formatted time string
		g_free(time_str);
	}

	// Add scrolled window to dialog
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_pack_start(GTK_BOX(content_area), scrolled_window, TRUE, TRUE, 0);
	gtk_widget_show_all(dialog);

	// Run the dialog
	GPid selected_pid = (GPid) -1;
	if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK) {
		GtkTreeModel *model;
		if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
			child_info *selected_child;
			gtk_tree_model_get(model, &iter, 0, &selected_child, -1);
			selected_pid = selected_child->childpid;
		}
	}

	// Clean up
	gtk_widget_destroy(dialog);

	return selected_pid;
}
