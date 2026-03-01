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

/*
 * GTK File Chooser static functions
 */

#include <stdio.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/initfile.h"
#include "core/icc_profile.h"
#include "algos/sorting.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/fits_sequence.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialog_preview.h"
#include "gui/conversion.h"
#include "gui/dialogs.h"

#include "open_dialog.h"

static void set_single_filter_dialog(GtkFileChooser *chooser, const gchar *name, const gchar *filter) {
	gtk_filter_add(chooser, name, filter, TRUE);
}

static void set_filters_dialog(GtkFileChooser *chooser, int whichdial) {
	GString *all_filter = NULL;
	gchar *fits_filter = FITS_EXTENSIONS;
	gchar *netpbm_filter = "*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM";
	gchar *pic_filter = "*.pic;*.PIC";
	gchar *ser_filter = "*.ser;*.SER";
	if (whichdial != OD_CONVERT && whichdial != OD_OPEN) {
		gtk_filter_add(chooser, _("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"),
				fits_filter, gui.file_ext_filter == TYPEFITS);
	} else {
		all_filter = g_string_new(fits_filter);
	}
	if (whichdial == OD_OPEN || whichdial == OD_CONVERT) {
#ifdef HAVE_LIBRAW
		/* RAW FILES */
		int nb_raw;
		char *raw;
		int i;

		nb_raw = get_nb_raw_supported();
		raw = calloc(sizeof(char), nb_raw * 12 + 1); // we assume the extension size of 3 char "*.xxx;*.XXX;" = 12
		for (i = 0; i < nb_raw; i++) {
			char *ext;
			gchar *upcase;

			upcase = g_ascii_strup(supported_raw[i].extension, strlen(supported_raw[i].extension));
			ext = g_strdup_printf("*.%s;*.%s;", supported_raw[i].extension,
					upcase);
			strcat(raw, ext);

			g_free(ext);
			g_free(upcase);
		}
		raw[strlen(raw) - 1] = '\0';

		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, raw);


		free(raw);

#endif
		/* GRAPHICS FILES */
		GString *s_supported_graph, *s_pattern;
		gchar *graphics_supported, *graphics_filter;

		s_supported_graph = g_string_new(_("Graphics Files (*.bmp"));
		s_pattern = g_string_new("*.bmp;*.BMP;");
#ifdef HAVE_LIBJPEG
		s_supported_graph = g_string_append(s_supported_graph, ", *.jpg, *.jpeg");
		s_pattern = g_string_append(s_pattern, "*.jpg;*.JPG;*.jpeg;*.JPEG;");
#endif

#ifdef HAVE_LIBJXL
		s_supported_graph = g_string_append(s_supported_graph, ", *.jxl");
		s_pattern = g_string_append(s_pattern, "*.jxl;*.JXL;");
#endif

#ifdef HAVE_LIBHEIF
		s_supported_graph = g_string_append(s_supported_graph, ", *.heic, *.heif, *.avif");
		s_pattern = g_string_append(s_pattern, "*.heic;*.HEIC;*.heif;*.HEIF;*.avif;*.AVIF;");
#endif

#ifdef HAVE_LIBPNG
		s_supported_graph = g_string_append(s_supported_graph, ", *.png");
		s_pattern = g_string_append(s_pattern, "*.png;*.PNG;");
#endif

#ifdef HAVE_LIBTIFF
		s_supported_graph = g_string_append(s_supported_graph, ", *.tif, *.tiff");
		s_pattern = g_string_append(s_pattern, "*.tif;*.TIF;*.tiff;*.TIFF;");
#endif

#ifdef HAVE_LIBXISF
		s_supported_graph = g_string_append(s_supported_graph, ", *.xisf");
		s_pattern = g_string_append(s_pattern, "*.xisf;*.XISF");
#endif
		s_supported_graph = g_string_append(s_supported_graph, ")");

		graphics_supported = g_string_free(s_supported_graph, FALSE);
		graphics_filter = g_string_free(s_pattern, FALSE);

		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, graphics_filter);
		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, netpbm_filter);
		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, pic_filter);
		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, ser_filter);

#ifdef HAVE_FFMS2
		/* FILM FILES */
		int nb_film;
		char *film_filter;
		int j;

		nb_film = get_nb_film_ext_supported();
		film_filter = calloc(sizeof(char), nb_film * 14 + 1);// we assume the extension size of 4 char "*.xxxx;*.XXXX;" = 14
		for (j = 0; j < nb_film; j++) {
			char *ext;
			gchar *upcase;

			upcase = g_ascii_strup(supported_film[j].extension,
					strlen(supported_film[j].extension));
			ext = g_strdup_printf("*.%s;*.%s;", supported_film[j].extension,
					upcase);
			strcat(film_filter, ext);

			g_free(ext);
			g_free(upcase);
		}
		if (strlen(film_filter) > 0)
			film_filter[strlen(film_filter) - 1] = '\0';

		all_filter = g_string_append(all_filter, ";");
		all_filter = g_string_append(all_filter, film_filter);

		free(film_filter);
#endif
		g_free(graphics_supported);
		g_free(graphics_filter);

		if (whichdial == OD_CONVERT || whichdial == OD_OPEN) {
			gchar *filter = g_string_free(all_filter, FALSE);

			gtk_filter_add(chooser, _("All supported files"), filter, TRUE);
			g_free(filter);
		}
	}
}

static void on_debayer_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	gtk_toggle_button_set_active((GtkToggleButton *)user_data, gtk_toggle_button_get_active(togglebutton));
}

static void siril_add_debayer_toggle_button(GtkFileChooser *dialog) {
	GtkToggleButton *main_debayer_button = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	GtkWidget *toggle_debayer = gtk_check_button_new_with_label(_("Debayer"));

	gtk_widget_show(toggle_debayer);
	gtk_file_chooser_set_extra_widget(dialog, toggle_debayer);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_debayer), gtk_toggle_button_get_active(main_debayer_button));
	g_signal_connect(GTK_TOGGLE_BUTTON(toggle_debayer), "toggled", G_CALLBACK(on_debayer_toggled), (gpointer) main_debayer_button);
}

static gchar* get_calibration_file_directory(GtkWidget *entry) {
    if (!GTK_IS_ENTRY(GTK_ENTRY(entry))) {
        return NULL;
    }

    // Get the filename from the file chooser button (returns absolute path)
    const gchar *filename = gtk_entry_get_text(GTK_ENTRY(entry));

    if (!filename || filename[0] == '\0') {
        // No file selected
        return g_strdup(com.wd);
    }

    // Since gtk_file_chooser_get_filename() always returns absolute paths,
    // we always extract the directory portion
    gchar *dirname = g_path_get_dirname(filename);
    return dirname;
}

static void opendial(int whichdial) {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	fileChooserPreview *preview = NULL;
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	gint res;
	int retval;

	if (!com.wd)
		return;

	GtkWidget *entry = NULL;
	switch (whichdial) {
	case OD_NULL:
		fprintf(stderr, "whichdial undefined, should not happen\n");
		return;
	case OD_FLAT:
		if (whichdial == OD_FLAT)
			entry = lookup_widget("flatname_entry");
	case OD_DARK:
		if (whichdial == OD_DARK)
			entry = lookup_widget("darkname_entry");
	case OD_OFFSET:
		if (whichdial == OD_OFFSET)
			entry = lookup_widget("offsetname_entry");
	case OD_FLATLIB:
		if (whichdial == OD_FLATLIB)
			entry = lookup_widget("flatlib_entry");
	case OD_DARKLIB:
		if (whichdial == OD_DARKLIB)
			entry = lookup_widget("darklib_entry");
	case OD_OFFSETLIB:
		if (whichdial == OD_OFFSETLIB)
			entry = lookup_widget("biaslib_entry");
	case OD_DISTOLIB:
		if (whichdial == OD_DISTOLIB)
			entry = lookup_widget("distolib_entry");
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gchar *lastdir = get_calibration_file_directory(entry);
		gtk_file_chooser_set_current_folder(dialog, lastdir);
		g_free(lastdir);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		if (whichdial == OD_FLATLIB) {
			if (com.pref.prepro.flat_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.flat_lib);
				gtk_file_chooser_set_current_folder(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_DARKLIB) {
			if (com.pref.prepro.dark_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.dark_lib);
				gtk_file_chooser_set_current_folder(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_OFFSETLIB) {
			if (com.pref.prepro.bias_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.bias_lib);
				gtk_file_chooser_set_current_folder(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_DISTOLIB) {
			if (com.pref.prepro.disto_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.disto_lib);
				gtk_file_chooser_set_current_folder(dialog, path);
				g_free(path);
			}
		}
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		if (whichdial == OD_DISTOLIB) {
			set_single_filter_dialog(dialog, _("Master-distortion: file (*.wcs)"), "*.wcs;*.WCS");
		} else {
			set_filters_dialog(dialog, whichdial);
		}
		siril_file_chooser_add_preview(dialog, preview);
		break;
	case OD_CWD:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		break;
	case OD_OPEN:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_filters_dialog(dialog, whichdial);
		siril_file_chooser_add_preview(dialog, preview);
		siril_add_debayer_toggle_button(dialog);
		break;
	case OD_CONVERT:
		widgetdialog = siril_file_chooser_add(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		gtk_file_chooser_set_select_multiple(dialog, TRUE);
		set_filters_dialog(dialog, whichdial);
		break;
	case OD_BADPIXEL:
		widgetdialog = siril_file_chooser_add(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_single_filter_dialog(dialog, _("Cosmetic correction file (*.lst)"), "*.lst;*.LST");
	}

	if (!dialog)
		return;

	wait:
	res = siril_dialog_run(widgetdialog);

	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *list = NULL;
		gchar *filename, *err;
		gboolean anything_loaded;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		GtkEntry *flat_entry, *dark_entry, *bias_entry, *bad_pixel_entry;
		GtkEntry *flatlib_entry, *darklib_entry, *biaslib_entry, *distolib_entry;
		GtkToggleButton *flat_button, *dark_button, *bias_button;
		GtkWidget *pbutton;

		filename = siril_file_chooser_get_filename(chooser);
		if (!filename)
			goto wait;

		if (fitseq_is_fitseq(filename, NULL)
				&& ((whichdial == OD_FLAT) || (whichdial == OD_DARK) || (whichdial == OD_OFFSET))) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Format not supported."),
					_("FITS sequences are not supported for master files. Please select a single FITS file."));
			goto wait;
		}

		pbutton  = lookup_widget("prepro_button");
		flat_entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		dark_entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		bias_entry = GTK_ENTRY(lookup_widget("offsetname_entry"));
		flatlib_entry = GTK_ENTRY(lookup_widget("flatlib_entry"));
		darklib_entry = GTK_ENTRY(lookup_widget("darklib_entry"));
		biaslib_entry = GTK_ENTRY(lookup_widget("biaslib_entry"));
		distolib_entry = GTK_ENTRY(lookup_widget("distolib_entry"));
		bad_pixel_entry = GTK_ENTRY(lookup_widget("pixelmap_entry"));

		flat_button = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
		dark_button = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
		bias_button = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));

		anything_loaded = sequence_is_loaded() || single_image_is_loaded();

		switch (whichdial) {
		case OD_FLAT:
			gtk_entry_set_text(flat_entry, filename);
			gtk_toggle_button_set_active(flat_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_DARK:
			gtk_entry_set_text(dark_entry, filename);
			gtk_toggle_button_set_active(dark_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_OFFSET:
			gtk_entry_set_text(bias_entry, filename);
			gtk_toggle_button_set_active(bias_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_FLATLIB:
			gtk_entry_set_text(flatlib_entry, filename);
			break;

		case OD_DARKLIB:
			gtk_entry_set_text(darklib_entry, filename);
			break;

		case OD_OFFSETLIB:
			gtk_entry_set_text(biaslib_entry, filename);
			break;

		case OD_DISTOLIB:
			gtk_entry_set_text(distolib_entry, filename);
			break;

		case OD_CWD:
			if (!siril_change_dir(filename, &err)) {
				if (com.pref.wd)
					g_free(com.pref.wd);
				com.pref.wd = g_strdup(com.wd);
				writeinitfile();
				gui_function(set_GUI_CWD, NULL);
			} else {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), err);
			}
			break;

		case OD_OPEN:
			set_cursor_waiting(TRUE);
			retval = open_single_image(filename);
			icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_LOAD);
			set_cursor_waiting(FALSE);
			if (retval == OPEN_IMAGE_CANCEL) goto wait;
			break;

		case OD_CONVERT:
			list = siril_file_chooser_get_filenames(chooser);
			list = g_slist_sort(list, (GCompareFunc) strcompare);
			fill_convert_list(list);
			g_slist_free(list);
			break;

		case OD_BADPIXEL:
			gtk_entry_set_text(bad_pixel_entry, filename);
			break;

		}
		g_free(filename);
	}
	siril_preview_free(preview);
	siril_widget_destroy(widgetdialog);
}

void cwd_btton_clicked() {
	opendial(OD_CWD);
}

void on_offsetfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_OFFSET);
}

void on_darkfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_DARK);
}

void on_flatfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_FLAT);
}

void on_offsetlibfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_OFFSETLIB);
}

void on_darklibfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_DARKLIB);
}

void on_flatlibfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_FLATLIB);
}

void on_distolibfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_DISTOLIB);
}

void header_open_button_clicked() {
	opendial(OD_OPEN);
}

void on_select_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	opendial(OD_CONVERT);
}

void on_pixelmap_button_clicked(GtkToggleButton *button, gpointer user_data) {
	opendial(OD_BADPIXEL);
}

/** callback function for recent document opened
 *
 * @param chooser
 * @param user_data
 */
void on_open_recent_action_item_activated(GtkRecentChooser *chooser,
		gpointer user_data) {
	gchar *uri, *path;
	GError *error = NULL;

	uri = gtk_recent_chooser_get_current_uri(chooser);

	path = g_filename_from_uri(uri, NULL, &error);
	if (error) {
		g_warning("Could not convert uri \"%s\" to a local path: %s", uri,
				error->message);
		g_clear_error(&error);
		g_free(uri);
		return;
	}

	open_single_image(path);
	icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_LOAD);

	g_free(uri);
	g_free(path);
}

