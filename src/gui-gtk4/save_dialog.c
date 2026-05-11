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

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/dialog_preview.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/icc_profile.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/callbacks.h"
#include "io/Astro-TIFF.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "save_dialog.h"

static image_type type_of_image = TYPEUNDEF;
static SirilFileChooser *saveDialog = NULL;

static GtkNotebook *sd_notebook_format = NULL;
static GtkWidget *sd_savepopup = NULL;
static GtkWidget *sd_filenameframe = NULL;
static GtkWidget *sd_button_savepopup = NULL;
static GtkWindow *sd_control_window = NULL;
static GtkEntry *sd_savetxt_entry = NULL;
static GtkCheckButton *sd_fits_8 = NULL, *sd_fits_16s = NULL, *sd_fits_16 = NULL;
static GtkCheckButton *sd_fits_32f = NULL;
static GtkCheckButton *sd_update_hilo = NULL, *sd_checksum = NULL;
static GtkEntry *sd_size_estimate_entry = NULL;
static GtkWidget *sd_header_snapshot_btn = NULL;
static GtkScrolledWindow *sd_scrolledwindow3 = NULL;

static void save_dialog_init_statics(void) {
	if (sd_notebook_format) return;
	sd_notebook_format = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebookFormat"));
	sd_savepopup = GTK_WIDGET(gtk_builder_get_object(gui.builder, "savepopup"));
	sd_filenameframe = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filenameframe"));
	sd_button_savepopup = GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_savepopup"));
	sd_control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	sd_savetxt_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "savetxt"));
	sd_fits_8 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton_save_fit8"));
	sd_fits_16s = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton_save_fit16s"));
	sd_fits_16 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton_save_fit16"));
	sd_fits_32f = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton_save_fit32f"));
	sd_update_hilo = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_update_hilo"));
	sd_checksum = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_chksum"));
	sd_size_estimate_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "size_estimate_entry"));
	sd_header_snapshot_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "header_snapshot_button"));
	sd_scrolledwindow3 = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scrolledwindow3"));
}

static void set_filters_save_dialog(SirilFileChooser *fc) {
	GString *all_filter = NULL;
	const gchar *fits_filter = "*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.FIT.fz;*.fits.fz;*.FITS.fz;*.fts.fz;*.FTS.fz";
	const gchar *bmp_filter = "*.bmp;*.BMP";
	const gchar *netpbm_filter = "*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM";

	all_filter = g_string_new(fits_filter);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, netpbm_filter);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, bmp_filter);

	siril_fc_add_filter_pattern(fc, _("FITS Files (*.fit, *.fits, *.fts)"), fits_filter, FALSE);
	siril_fc_add_filter_pattern(fc, _("BMP Files (*.bmp)"), bmp_filter, FALSE);
#ifdef HAVE_LIBJPEG
	const gchar *jpg_filter = "*.jpg;*.JPG;*.jpeg;*.JPEG";
	siril_fc_add_filter_pattern(fc, _("JPEG Files (*.jpg, *.jpeg)"), jpg_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, jpg_filter);
#endif

#ifdef HAVE_LIBJXL
	const gchar *jxl_filter = "*.jxl;*.JXL";
	siril_fc_add_filter_pattern(fc, _("JPEG XL Files (*.jxl)"), jxl_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, jxl_filter);
#endif

#ifdef HAVE_LIBPNG
	const gchar *png_filter = "*.png;*.PNG";
	siril_fc_add_filter_pattern(fc, _("PNG Files (*.png)"), png_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, png_filter);
#endif

#ifdef HAVE_LIBTIFF
	const gchar *tif_filter = "*.tif;*.TIF;*.tiff;*.TIFF";
	siril_fc_add_filter_pattern(fc, _("TIFF Files (*.tif, *.tiff)"), tif_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, tif_filter);
#endif

	/* NETPBM FILES */
	siril_fc_add_filter_pattern(fc, _("Netpbm Files (*.ppm, *.pnm, *.pgm)"), netpbm_filter, FALSE);

	gchar *supported_file = g_string_free(all_filter, FALSE);
	siril_fc_add_filter_pattern(fc, _("Supported Image Files"), supported_file, TRUE);

	g_free(supported_file);
}

static image_type get_filetype(const gchar *filter) {
	image_type type = TYPEUNDEF;
	int i = 0;

	gchar **string = g_strsplit_set(filter, "*(),.", -1);

	while (string[i]) {
		if (!g_strcmp0(string[i], "fit")) {
			type = TYPEFITS;
			break;
		} else if (!g_strcmp0(string[i], "bmp")) {
			type = TYPEBMP;
			break;
		} else if (!g_strcmp0(string[i], "jpg")) {
			type = TYPEJPG;
			break;
		} else if (!g_strcmp0(string[i], "png")) {
			type = TYPEPNG;
			break;
		} else if (!g_strcmp0(string[i], "tif")) {
			type = TYPETIFF;
			break;
		} else if (!g_strcmp0(string[i], "jxl")) {
			type = TYPEJXL;
			break;
		} else if (!g_strcmp0(string[i], "ppm")) {
			type = TYPEPNM;
			break;
		}
		i++;
	}
	g_strfreev(string);

	return type;
}

static void set_copyright_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextBuffer *tbuf;
	GtkTextIter itStart, itEnd;
	gchar *copyright;

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(gtk_builder_get_object(gui.builder, "Copyright_txt"));

	tbuf = gtk_text_view_get_buffer(TIFF_txt);

	if (com.pref.copyright && (com.pref.copyright[0] != '\0')) {
		copyright = g_strdup(com.pref.copyright);
	} else {
		copyright = g_strdup_printf("%c%s v%s", toupper(PACKAGE[0]), (char *)PACKAGE + 1, VERSION);
	}

	gtk_text_buffer_get_bounds(tbuf, &itStart, &itEnd);
	gtk_text_buffer_delete(tbuf, &itStart, &itEnd);
	gtk_text_buffer_set_text(tbuf, copyright, strlen(copyright));

	g_free(copyright);
}

static void set_description_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextIter itStart, itEnd;

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(gtk_builder_get_object(gui.builder, "Description_txt"));

	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(TIFF_txt);

	gtk_text_buffer_get_bounds(tbuf, &itStart, &itEnd);
	gtk_text_buffer_delete(tbuf, &itStart, &itEnd);

	if (gtk_drop_down_get_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_type_of_tiff"))) == 0) {
		gchar *astro_tiff = AstroTiff_build_header(gfit);

		gtk_text_buffer_get_end_iter(tbuf, &itEnd);
		gtk_text_buffer_insert(tbuf, &itEnd, astro_tiff, -1);
		gtk_text_buffer_get_end_iter(tbuf, &itEnd);
		gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
		g_free(astro_tiff);
	} else {
		/* History already written in header */
		if (gfit->history) {
			GSList *list;
			for (list = gfit->history; list; list = list->next) {
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, (gchar *)list->data, -1);
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
			}
		}
		/* New history.  com.undo_stack is head=most-recent; walk it
		 * back-to-front so the TIFF description reads in chronological
		 * order. */
		for (GList *l = g_list_last(com.undo_stack); l; l = l->prev) {
			historic *h = (historic *)l->data;
			if (h->history[0] != '\0') {
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, h->history, strlen(h->history));
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
			}
		}
	}
}

static void prepare_savepopup() {
	save_dialog_init_statics();
	GtkWidget *savepopup = sd_savepopup;
	GtkWidget *savetxt = sd_filenameframe;
	GtkWidget *button_savepopup = sd_button_savepopup;
	int tab;

	GtkWindow *parent = sd_control_window;

	switch (type_of_image) {
	case TYPEBMP:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		tab = PAGE_MISC;
		break;
	case TYPEPNG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving PNG"));
		tab = PAGE_MISC;
		break;
	case TYPEPNM:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		tab = PAGE_MISC;
		break;
	case TYPEJPG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPG"));
		tab = PAGE_JPG;
		break;
	case TYPEJXL:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPEG XL"));
		tab = PAGE_JXL;
		break;
	case TYPETIFF:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		set_copyright_in_TIFF();
		set_description_in_TIFF();
		set_icc_description_in_TIFF();
		tab = PAGE_TIFF;
		break;
	default:
	case TYPEFITS:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving FITS"));
		tab = PAGE_FITS;
	}

	gtk_window_set_transient_for(GTK_WINDOW(savepopup), parent);

	gtk_widget_set_visible(savetxt, FALSE);
	gtk_notebook_set_current_page(sd_notebook_format, tab);
	gtk_widget_grab_focus(button_savepopup);
}

static void init_dialog() {
	if (saveDialog == NULL) {
		GtkWindow *parent = siril_get_active_window();
		saveDialog = siril_fc_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
		set_filters_save_dialog(saveDialog);
	}
}

static void close_dialog() {
	if (saveDialog != NULL) {
		siril_fc_destroy(saveDialog);
		saveDialog = NULL;
	}
}

static gchar *get_filename_and_replace_ext() {
	gchar *basename;

	if (sequence_is_loaded() && com.seq.current != RESULT_IMAGE) {
		char fname[256];
		/* set the output file name default as the current image.jpg */
		seq_get_image_filename(&com.seq, com.seq.current, fname);
		basename = g_path_get_basename(fname);
	} else {
		basename = g_path_get_basename(com.uniq->filename);
	}

	gboolean is_format_valid = get_type_from_filename(basename)
			& (TYPEFITS | TYPEBMP | TYPETIFF | TYPEPNG | TYPEJPG | TYPEPNM);
	if (!is_format_valid) {
		char *file_no_ext = remove_ext_from_filename(basename);
		g_free(basename);
		basename = g_strdup_printf("%s%s", file_no_ext, com.pref.ext);
		free(file_no_ext);
	}

	return basename;
}

static image_type get_image_type_from_filter(GtkFileFilter *filter) {
	return get_filetype(gtk_file_filter_get_name(filter));
}

/* Adjust the basename's extension to match the filter's expected format.
 * GtkFileDialog has no notify::filter signal, so we reconcile only after
 * the user accepts: return a possibly-corrected newly-allocated filename
 * (caller must g_free) or NULL if no change is needed. */
static gchar *reconcile_filename_with_filter(const gchar *filename, GtkFileFilter *filter) {
	if (!filename || !filter) return NULL;
	char *file_no_ext = remove_ext_from_filename(filename);
	image_type format = get_image_type_from_filter(filter);
	gchar *new_filename = NULL;

	switch (format) {
	default:
	case TYPEFITS:
		new_filename = g_strdup_printf("%s%s", file_no_ext, com.pref.ext);
		break;
	case TYPEBMP:
		new_filename = g_strdup_printf("%s.bmp", file_no_ext);
		break;
	case TYPEPNM:
		new_filename = g_strdup_printf("%s%s", file_no_ext,
		                              (gfit->naxes[2] == 1) ? ".pgm" : ".ppm");
		break;
#ifdef HAVE_LIBTIFF
	case TYPETIFF:
		new_filename = g_strdup_printf("%s.tif", file_no_ext);
		break;
#endif
#ifdef HAVE_LIBJPEG
	case TYPEJPG:
		new_filename = g_strdup_printf("%s.jpg", file_no_ext);
		break;
#endif
#ifdef HAVE_LIBJXL
	case TYPEJXL:
		new_filename = g_strdup_printf("%s.jxl", file_no_ext);
		break;
#endif
#ifdef HAVE_LIBPNG
	case TYPEPNG:
		new_filename = g_strdup_printf("%s.png", file_no_ext);
		break;
#endif
	}
	free(file_no_ext);
	return new_filename;
}

gchar* extract_extension_from_filter(const gchar *filter_name) {
	// Regex pattern to match first extension inside parentheses
	GRegex *regex;
	GError *error = NULL;
	gchar *extension = NULL;

	// Compile regex to match first extension after *. inside parentheses
	regex = g_regex_new("\\*\\.([a-zA-Z0-9]+)", G_REGEX_MULTILINE, 0, &error);

	if (error != NULL) {
		g_warning("Regex compilation error: %s", error->message);
		g_error_free(error);
		return NULL;
	}

	// Try to match the regex
	GMatchInfo *match_info;
	if (g_regex_match(regex, filter_name, 0, &match_info)) {
		// Get the first captured group (extension)
		extension = g_match_info_fetch(match_info, 1);
	}

	// Free resources
	g_match_info_free(match_info);
	g_regex_unref(regex);

	return extension;
}

static int save_dialog() {
	init_dialog();

	gchar *fname = get_filename_and_replace_ext();
	siril_fc_set_current_name(saveDialog, fname);
	g_free(fname);

	while (TRUE) {
		gint res = siril_fc_run(saveDialog);

		if (res != GTK_RESPONSE_ACCEPT) {
			close_dialog();
			return res;
		}

		gchar *filename = siril_fc_get_filename(saveDialog);
		if (!filename) {
			close_dialog();
			return GTK_RESPONSE_CANCEL;
		}

		GtkFileFilter *current_filter = siril_fc_get_filter(saveDialog);
		const gchar *filter_name = current_filter ? gtk_file_filter_get_name(current_filter) : NULL;
		gchar *filter_extension = filter_name ? extract_extension_from_filter(filter_name) : NULL;

		image_type file_type = get_type_from_filename(filename);
		image_type filter_format = filter_extension ? get_type_for_extension(filter_extension) : TYPEUNDEF;
		image_type selected_type = (filter_format != TYPEUNDEF) ? filter_format : file_type;

		if (file_type != TYPEUNDEF && filter_format != TYPEUNDEF && file_type != filter_format) {
			gboolean confirm = siril_confirm_dialog(_("Extension Mismatch"),
							_("The given file extension does not match the chosen file type. Do you want to save the image using this name anyway?"),
							_("Save Anyway"));

			if (!confirm) {
				/* User declined: reconcile the basename to the filter's
				 * expected extension and re-show the dialog. */
				gchar *fixed = reconcile_filename_with_filter(filename, current_filter);
				if (fixed) {
					gchar *bname = g_path_get_basename(fixed);
					siril_fc_set_current_name(saveDialog, bname);
					g_free(bname);
					g_free(fixed);
				}
				g_free(filename);
				g_free(filter_extension);
				continue;
			}

			selected_type = file_type;
		}

		type_of_image = selected_type;

		GtkEntry *savetext = sd_savetxt_entry;
		gtk_editable_set_text(GTK_EDITABLE(savetext), filename);

		prepare_savepopup();
		g_free(filename);
		g_free(filter_extension);

		return res;
	}
}

// idle function executed at the end of the Save Data processing
static gboolean end_save(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
	stop_processing_thread();

	if (args->retval)
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("File saving failed. Check the logs for more info."));

	gtk_editable_set_text(GTK_EDITABLE(args->entry), "");
	gtk_widget_set_visible(sd_savepopup, FALSE);
	display_filename(); // update filename display
	adjust_sellabel();
	gui_function(set_precision_switch, NULL);
	set_cursor_waiting(FALSE);
	close_dialog();
	gui_function(update_MenuItem, NULL);

	g_free(args->copyright);
	g_free(args->description);
	free(args);
	return FALSE;
}

static gboolean test_for_viewer_mode() {
	gboolean confirm = TRUE;
	if (get_display_mode_from_menu() != LINEAR_DISPLAY && !com.pref.gui.silent_linear) {
		confirm = siril_confirm_dialog_and_remember(_("Viewer mode is not linear"),
				_("You are saving an image that is displayed in a non linear mode. "
						"What you see is not what you will get. Switch the viewer mode to linear to save what you see."), _("Save Anyway"), &com.pref.gui.silent_linear);
	}
	return confirm;
}

static void initialize_data(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;

	save_dialog_init_statics();
	GtkCheckButton *fits_8 = sd_fits_8;
	GtkCheckButton *fits_16s = sd_fits_16s;
	GtkCheckButton *fits_16 = sd_fits_16;
	GtkCheckButton *update_hilo = sd_update_hilo;
	GtkCheckButton *checksum = sd_checksum;
#ifdef HAVE_LIBJPEG
	GtkSpinButton *qlty_spin_button = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "quality_spinbutton"));
	args->quality = gtk_spin_button_get_value_as_int(qlty_spin_button);
#endif
#ifdef HAVE_LIBJXL
	GtkSpinButton *quality_spin_button = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "jxl_quality_spinbutton"));
	args->jxl_quality = gtk_spin_button_get_value_as_int(quality_spin_button);
	GtkSpinButton *effort_spin_button = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "jxl_effort_spinbutton"));
	args->jxl_effort = gtk_spin_button_get_value_as_int(effort_spin_button);
	GtkCheckButton *toggle_button_8bit = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "jxl_force_8bit"));
	args->jxl_force_8bit = siril_toggle_get_active(GTK_WIDGET(toggle_button_8bit));
#endif
#ifdef HAVE_LIBTIFF
	GtkCheckButton *button_8 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton8bits"));
	GtkCheckButton *button_32 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton32bits"));
	args->bitspersamples = siril_toggle_get_active(GTK_WIDGET(button_8)) ? 8 : siril_toggle_get_active(GTK_WIDGET(button_32)) ? 32 : 16;
	get_tif_data_from_ui(gfit, &args->description, &args->copyright);
	args->tiff_compression = get_tiff_compression();
#endif
	args->entry = sd_savetxt_entry;
	args->filename = gtk_editable_get_text(GTK_EDITABLE(args->entry));

	if (siril_toggle_get_active(GTK_WIDGET(fits_8)))
		args->bitpix = BYTE_IMG;
	else if (siril_toggle_get_active(GTK_WIDGET(fits_16s)))
		args->bitpix = SHORT_IMG;
	else if (siril_toggle_get_active(GTK_WIDGET(fits_16)))
		args->bitpix = USHORT_IMG;
	else
		args->bitpix = FLOAT_IMG;

	args->update_hilo = siril_toggle_get_active(GTK_WIDGET(update_hilo));
	args->checksum = siril_toggle_get_active(GTK_WIDGET(checksum));
}

#ifdef HAVE_LIBJPEG
static long calculate_jpeg_size(struct savedial_data *args) {
    gchar *tmp_template;
    int fd;
    long file_size = -1;

    /* Build a secure temporary filename template */
    tmp_template = g_build_filename(
        g_get_tmp_dir(),
        "siril-preview-XXXXXX.jpg",
        NULL
    );

    fd = g_mkstemp(tmp_template);
    if (fd == -1) {
        g_warning("Failed to create temporary file: %s", g_strerror(errno));
        g_free(tmp_template);
        return -1;
    }

    /* We only need the filename; savejpg() will reopen it */
    close(fd);

    /* Save JPEG */
    if (savejpg(tmp_template, gfit, args->quality, FALSE)) {
        siril_debug_print("Failed to save JPEG to temporary file");
        goto cleanup;
    }

    /* Get file size */
    GStatBuf stat_buf;
    if (g_stat(tmp_template, &stat_buf) == 0) {
        file_size = (long) stat_buf.st_size;
    } else {
        g_warning("Unable to stat '%s': %s", tmp_template, g_strerror(errno));
    }

cleanup:
    if (g_unlink(tmp_template) != 0) {
        siril_debug_print("g_unlink() failed for temp file\n");
    }

    g_free(tmp_template);
    return file_size;
}

struct jpeg_size_result {
    struct savedial_data *args;
    gchar *txt;
};

static gboolean end_calculate_jpeg_size(gpointer p) {
    struct jpeg_size_result *res = p;

    gtk_editable_set_text(GTK_EDITABLE(sd_size_estimate_entry), res->txt);

    stop_processing_thread();
    set_cursor_waiting(FALSE);

    g_free(res->txt);
    g_free(res->args->copyright);
    g_free(res->args->description);
    free(res->args);
    g_free(res);

    return FALSE;
}

static gpointer calculate_jpeg_size_thread(gpointer p) {
    struct savedial_data *args = p;
    long file_size = calculate_jpeg_size(args);

    struct jpeg_size_result *res = g_new0(struct jpeg_size_result, 1);
    res->args = args;
    res->txt = g_format_size_full(file_size, G_FORMAT_SIZE_IEC_UNITS);

    siril_add_idle(end_calculate_jpeg_size, res);
    return NULL;
}

#else

static gboolean end_calculate_jpeg_size(gpointer p) {
    stop_processing_thread();
    set_cursor_waiting(FALSE);
    return FALSE;
}

static gpointer calculate_jpeg_size_thread(gpointer p) {
	printf("Should not happen\n");
	siril_add_idle(end_calculate_jpeg_size, NULL);
	return NULL;
}

#endif

static gpointer mini_save_dialog(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
#ifdef HAVE_LIBPNG
	uint32_t bytes_per_sample;
#endif
	args->retval = 0;

	if (args->filename[0] != '\0') {
		switch (type_of_image) {
		case TYPEBMP:
			args->retval = savebmp(args->filename, gfit);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:;
			args->retval = savejpg(args->filename, gfit, args->quality, TRUE);
			break;
#endif
#ifdef HAVE_LIBJXL
		case TYPEJXL:
			args->retval = savejxl(args->filename, gfit, args->jxl_effort, args->jxl_quality, args->jxl_force_8bit);
			break;
#endif
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			args->retval = savetif(args->filename, gfit, args->bitspersamples, args->description, args->copyright, args->tiff_compression, TRUE, TRUE);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			bytes_per_sample = gfit->orig_bitpix != BYTE_IMG ? 2 : 1;
			args->retval = savepng(args->filename, gfit, bytes_per_sample, gfit->naxes[2] == 3);
			break;
#endif
		default:
			type_of_image = TYPEFITS;
			/* no break */
		case TYPEFITS:
			gfit->bitpix = args->bitpix;
			/* Check if MIPS-HI and MIPS-LO must be updated. If yes,
			 * Values are taken from the layer 0 */
			if (args->update_hilo) {
				gfit->keywords.hi = gui.hi;
				gfit->keywords.lo = gui.lo;
				if (gfit->orig_bitpix == BYTE_IMG
						&& (gfit->keywords.hi > UCHAR_MAX || gfit->keywords.lo > UCHAR_MAX)) {
					gfit->keywords.hi = UCHAR_MAX;
					gfit->keywords.lo = 0;
				} else if (gfit->orig_bitpix == SHORT_IMG
						&& (gfit->keywords.hi > SHRT_MAX || gfit->keywords.lo > SHRT_MAX)) {
					gfit->keywords.hi = UCHAR_MAX;
					gfit->keywords.lo = 0;
				}
				if (gfit->orig_bitpix == BYTE_IMG && gfit->bitpix != BYTE_IMG) {
					gfit->keywords.hi = USHRT_MAX;
					gfit->keywords.lo = 0;
				}
			} else {
				gfit->keywords.hi = 0;
				gfit->keywords.lo = 0;
			}
			gfit->checksum = args->checksum;
			args->retval = savefits(args->filename, gfit);
			if (!args->retval && single_image_is_loaded()) {
				/* Use the same extension logic as savefits() so the tracked
				 * filename always matches the file actually written on disk. */
				gchar *actual = set_right_extension(args->filename);
				free(com.uniq->filename);
				com.uniq->filename = strdup(actual);
				g_free(actual);
				com.uniq->fileexist = TRUE;
			}
			break;
		case TYPEPNM:
			args->retval = saveNetPBM(args->filename, gfit);
			break;
		}
	}
	siril_add_idle(end_save, args);
	return NULL;
}

void set_entry_filename() {
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		char filename[256];
		/* set the output file name default as the current image.jpg */
		GtkEntry *entry = sd_savetxt_entry;
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		char *file_no_ext = remove_ext_from_filename(filename);
		gtk_editable_set_text(GTK_EDITABLE(entry), file_no_ext);
		free(file_no_ext);
	}
}

void on_savetxt_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(editable);
	GtkWidget *button = sd_button_savepopup;

	const gchar *name = gtk_editable_get_text(GTK_EDITABLE(entry));
	gtk_widget_set_sensitive(button, (*name != '\0'));
}

void on_size_estimate_toggle_toggled(GtkCheckButton *button, gpointer user_data) {
	struct savedial_data *args = calloc(1, sizeof(struct savedial_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	initialize_data(args);

	if (!processing_is_job_active() && siril_toggle_get_active(GTK_WIDGET(button))) {
		if (!processing_is_job_active()) {
			if (!start_in_new_thread(calculate_jpeg_size_thread, args)) {
				g_free(args->copyright);
				g_free(args->description);
				free(args);
			}
			return;
		}
	}
	// Clear preview size
	gtk_editable_set_text(GTK_EDITABLE(sd_size_estimate_entry), "");

	g_free(args->description);
	g_free(args->copyright);
	free (args);
}

void on_quality_spinbutton_value_changed(GtkSpinButton *button, gpointer user_data) {
	if (siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "size_estimate_toggle"))))) {
		struct savedial_data *args = calloc(1, sizeof(struct savedial_data));
		if (!args) {
			PRINT_ALLOC_ERR;
			return;
		}
		initialize_data(args);
		if (!processing_is_job_active()) {
			if (!start_in_new_thread(calculate_jpeg_size_thread, args)) {
				g_free(args->copyright);
				g_free(args->description);
				free(args);
			}
		} else {
			g_free(args->description);
			g_free(args->copyright);
			free (args);
			return;
		}
	}
}

void on_button_savepopup_clicked(GtkButton *button, gpointer user_data) {
	struct savedial_data *args = calloc(1, sizeof(struct savedial_data));

	set_cursor_waiting(TRUE);
	initialize_data(args);
	if (test_for_viewer_mode()) {
		if (!start_in_new_thread(mini_save_dialog, args)) {
			g_free(args->copyright);
			g_free(args->description);
			free(args);
		}
	} else {
		g_free(args->copyright);
		g_free(args->description);
		free(args);
		siril_add_idle(end_generic, NULL);
	}
}

void on_savetxt_activate(GtkEntry *entry, gpointer user_data) {
	struct savedial_data *args = calloc(1, sizeof(struct savedial_data));

	set_cursor_waiting(TRUE);
	initialize_data(args);
	if (test_for_viewer_mode()) {
		if (!start_in_new_thread(mini_save_dialog, args)) {
			g_free(args->copyright);
			g_free(args->description);
			free(args);
		}
	} else {
		g_free(args->copyright);
		g_free(args->description);
		free(args);
		siril_add_idle(end_generic, NULL);
	}
}

void on_button_cancelpopup_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_set_visible(sd_savepopup, FALSE);
}

void on_header_save_as_button_clicked() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		GtkWidget *savepopup = sd_savepopup;

		if (save_dialog() == GTK_RESPONSE_ACCEPT) {
			/* now it is not needed for some formats */
			if (type_of_image & (TYPEBMP | TYPEPNG | TYPEPNM)) {
				struct savedial_data *args = calloc(1, sizeof(struct savedial_data));

				set_cursor_waiting(TRUE);
				initialize_data(args);
				if (test_for_viewer_mode()) {
					if (!start_in_new_thread(mini_save_dialog, args)) {
						g_free(args->description);
						g_free(args->copyright);
						free(args);
					}
				} else {
					g_free(args->copyright);
					g_free(args->description);
					free(args);
					siril_add_idle(end_generic, NULL);
				}
			} else {
				close_dialog();
				GtkWindow *parent = siril_get_active_window();
				if (!GTK_IS_WINDOW(parent)) {
					parent = sd_control_window;
				}
				gtk_window_set_transient_for(GTK_WINDOW(savepopup), parent);
				gtk_widget_set_visible(savepopup, TRUE);
				gtk_window_present(GTK_WINDOW(savepopup));
			}
		}
	}
}

static gboolean snapshot_notification_close(gpointer user_data) {
	GtkWidget *popover = GTK_WIDGET(user_data);
	/* gtk_widget_set_visible(FALSE) hides the popover widget but does
	 * NOT release the GTK4 popover input grab — that grab keeps
	 * swallowing every click in the application until the popover is
	 * finalised.  gtk_popover_popdown() ends the grab cleanly; we then
	 * unparent so the popover is freed once its remaining refs (the
	 * paintable / Cairo surface chain) drop. */
	if (GTK_IS_POPOVER(popover)) {
		gtk_popover_popdown(GTK_POPOVER(popover));
		gtk_widget_unparent(popover);
	}
	return FALSE;
}

static GtkWidget *snapshot_notification(GtkWidget *widget, const gchar *filename, GdkPaintable *paintable) {
	gchar *text;
	if (filename)
		text = g_strdup_printf(_("Snapshot <b>%s</b> was saved into the working directory."), filename);
	else
		text = g_strdup((_("Snapshot was saved into the clipboard.")));
	GtkWidget *popover = popover_new_with_image(widget, text, paintable);
	gtk_popover_popup(GTK_POPOVER(popover));
	g_free(text);
	return popover;
}

void on_header_snapshot_button_clicked(gboolean clipboard) {
	gchar *timestamp, *filename;
	GtkWidget *widget;
	const gchar *area[] = {"drawingarear", "drawingareag", "drawingareab", "drawingareargb" };

	/* Dismiss the snapshot menu popover (it doesn't auto-popdown on
	 * action-button clicks the way GtkPopoverMenu would).  Without this
	 * the popover's input grab swallows every subsequent click and the
	 * application appears frozen until the menu times out — which it
	 * won't, since it has no timeout.  Safe to call when invoked from
	 * a keyboard shortcut (the popover is hidden/popped-down already). */
	GObject *snap_menu = gtk_builder_get_object(gui.builder, "snapshot_menu");
	if (snap_menu && GTK_IS_POPOVER(snap_menu))
		gtk_popover_popdown(GTK_POPOVER(snap_menu));

	widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, area[gui.cvport]));
	timestamp = build_timestamp_filename();
	filename = g_strdup_printf("%s.png", timestamp);

	g_free(timestamp);
	/* create cr from the surface */
	int width = max(gfit->rx, gtk_widget_get_width(widget));
	int height = max(gfit->ry, gtk_widget_get_height(widget));
	/* Phase 2 removed view->full_surface; create a fresh image surface
	 * (cairo_surface_create_similar needed a non-NULL template). */
	cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
	cairo_t *cr = cairo_create(surface);

	/* add image and all annotations if available */
	add_image_and_label_to_cairo(cr, gui.cvport);

	cairo_destroy(cr);

	double z = get_zoom_val();

	gint x1 = max(0, (int) gui.display_offset.x);
	gint y1 = max(0, (int) gui.display_offset.y);

	gint x2 = min(width * z + (int) gui.display_offset.x, gtk_widget_get_width(widget));
	gint y2 = min(height * z + (int) gui.display_offset.y, gtk_widget_get_height(widget));

	/* Crop into a fresh ARGB32 image surface, then build a GdkTexture
	 * from it.  Cairo sub-surfaces aren't image surfaces and don't have
	 * a get_data accessor, so we paint the cropped region into a new
	 * image surface that the memory-texture helper can wrap. */
	cairo_surface_t *sub = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32, x2 - x1, y2 - y1);
	cairo_t *scr = cairo_create(sub);
	cairo_set_source_surface(scr, surface, -x1, -y1);
	cairo_paint(scr);
	cairo_destroy(scr);
	cairo_surface_destroy(surface);

	GdkTexture *tex = siril_texture_from_cairo_surface(sub);
	cairo_surface_destroy(sub);

	if (!tex) {
		g_free(filename);
		return;
	}

	if (clipboard) {
		GdkClipboard *cb = gdk_display_get_clipboard(gdk_display_get_default());
		gdk_clipboard_set_texture(cb, tex);
		GtkWidget *popover = snapshot_notification(sd_header_snapshot_btn, NULL,
		                                           GDK_PAINTABLE(tex));
		g_timeout_add(5000, (GSourceFunc) snapshot_notification_close, (gpointer) popover);
	} else {
		GFile *file = g_file_new_build_filename(com.wd, filename, NULL);
		gchar *path = g_file_get_path(file);
		/* gdk_texture_save_to_png is synchronous, but PNG-encoding a
		 * screen-sized snapshot is fast enough (single-digit ms) that
		 * the async dance with gdk_pixbuf_save_to_stream_async isn't
		 * worth keeping just to dodge a brief main-loop hitch. */
		gboolean ok = path && gdk_texture_save_to_png(tex, path);
		if (ok) {
			gchar *basename = g_file_get_basename(file);
			GtkWidget *popover = snapshot_notification(sd_header_snapshot_btn,
			                                           basename, GDK_PAINTABLE(tex));
			g_timeout_add(5000, (GSourceFunc) snapshot_notification_close, (gpointer) popover);
			g_free(basename);
		} else {
			siril_log_message(_("Cannot take snapshot\n"));
		}
		g_free(path);
		g_object_unref(file);
	}
	g_object_unref(tex);
	g_free(filename);
}

#undef NEW_SIZE

void on_header_save_button_clicked() {
	if (single_image_is_loaded() && com.uniq->fileexist) {
		savefits(com.uniq->filename, gfit);
	}
}

void on_savepopup_show(GtkWidget *widget, gpointer user_data) {
	save_dialog_init_statics();
	GtkScrolledWindow *scrolled_window = sd_scrolledwindow3;
	gint height, width;

	if (type_of_image & TYPETIFF) {
		GtkCheckButton *b16 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton16bits"));
		GtkCheckButton *b32 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "radiobutton32bits"));
		siril_toggle_set_active(GTK_WIDGET(b32), gfit->type == DATA_FLOAT);
		siril_toggle_set_active(GTK_WIDGET(b16), gfit->type == DATA_USHORT);
		width = 400;
		height = 100;
	} else {
		width = 100;
		height = 50;
	}
	/* we need to handle TYPEUNDEF in the case the extension
	 * is not explicitly written. In this case the FITS format will be chosen
	 */
	if (type_of_image & (TYPEFITS | TYPEUNDEF)) {
		GtkCheckButton *b16bitu = sd_fits_16;
		GtkCheckButton *b32bits = sd_fits_32f;
		siril_toggle_set_active(GTK_WIDGET(b32bits), gfit->type == DATA_FLOAT);
		siril_toggle_set_active(GTK_WIDGET(b16bitu), gfit->type == DATA_USHORT);
	}

	gtk_scrolled_window_set_min_content_height(scrolled_window, height);
	gtk_scrolled_window_set_min_content_width(scrolled_window, width);
}

void on_combo_type_of_tiff_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	set_description_in_TIFF();
}

gboolean get_tiff_compression(void) {
	if (!com.headless) {
		GtkCheckButton *button = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "radiobuttonCompDeflate")));
		if (siril_toggle_get_active(GTK_WIDGET(button)))
			return TRUE;
	}
	return FALSE;
}

void get_tif_data_from_ui(fits *fit, gchar **description, gchar **copyright) {
	if (!com.script && !com.headless) {
		/*******************************************************************
		 * If the user saves a tif from the graphical menu, he can set
		 * the Description and the Copyright of the Image
		 ******************************************************************/

		GtkTextIter itDebut;
		GtkTextIter itFin;

		GtkTextView *description_txt_view = GTK_TEXT_VIEW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "Description_txt")));
		GtkTextBuffer *desbuf = gtk_text_view_get_buffer(description_txt_view);
		gtk_text_buffer_get_start_iter(desbuf, &itDebut);
		gtk_text_buffer_get_end_iter(desbuf, &itFin);
		*description = gtk_text_buffer_get_text(desbuf, &itDebut, &itFin, TRUE);

		GtkTextView *copyright_txt_view = GTK_TEXT_VIEW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "Copyright_txt")));
		GtkTextBuffer *copybuf = gtk_text_view_get_buffer(copyright_txt_view);
		gtk_text_buffer_get_start_iter(copybuf, &itDebut);
		gtk_text_buffer_get_end_iter(copybuf, &itFin);
		*copyright = gtk_text_buffer_get_text(copybuf, &itDebut, &itFin, TRUE);

	}
}
