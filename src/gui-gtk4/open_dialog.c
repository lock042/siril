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
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "core/initfile.h"
#include "core/icc_profile.h"
#include "algos/sorting.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/fits_sequence.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/dialog_preview.h"
#include "gui-gtk4/conversion.h"
#include "gui-gtk4/dialogs.h"

#include "open_dialog.h"

static GtkCheckButton *od_demosaicing_btn = NULL;
static GtkWindow *od_control_window = NULL;
static GtkEntry *od_flatname_entry = NULL;
static GtkEntry *od_darkname_entry = NULL;
static GtkEntry *od_offsetname_entry = NULL;
static GtkEntry *od_flatlib_entry = NULL;
static GtkEntry *od_darklib_entry = NULL;
static GtkEntry *od_biaslib_entry = NULL;
static GtkEntry *od_distolib_entry = NULL;
static GtkWidget *od_prepro_button = NULL;
static GtkEntry *od_pixelmap_entry = NULL;
static GtkCheckButton *od_useflat_btn = NULL;
static GtkCheckButton *od_usedark_btn = NULL;
static GtkCheckButton *od_useoffset_btn = NULL;

static void open_dialog_init_statics(void) {
	if (od_control_window) return;
	od_demosaicing_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "demosaicingButton"));
	od_control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	od_flatname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "flatname_entry"));
	od_darkname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "darkname_entry"));
	od_offsetname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "offsetname_entry"));
	od_flatlib_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "flatlib_entry"));
	od_darklib_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "darklib_entry"));
	od_biaslib_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "biaslib_entry"));
	od_distolib_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "distolib_entry"));
	od_prepro_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "prepro_button"));
	od_pixelmap_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "pixelmap_entry"));
	od_useflat_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "useflat_button"));
	od_usedark_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "usedark_button"));
	od_useoffset_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "useoffset_button"));
}

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
	siril_toggle_set_active(GTK_WIDGET((GtkToggleButton *)user_data), siril_toggle_get_active(GTK_WIDGET(togglebutton)));
}

static void siril_add_debayer_toggle_button(GtkFileChooser *dialog) {
	open_dialog_init_statics();
	GtkCheckButton *main_debayer_button = od_demosaicing_btn;
	GtkWidget *toggle_debayer = gtk_check_button_new_with_label(_("Debayer"));

	gtk_widget_set_visible(toggle_debayer, TRUE);
	/* GTK4: gtk_file_chooser_set_extra_widget removed */;
	siril_toggle_set_active(GTK_WIDGET(toggle_debayer), siril_toggle_get_active(GTK_WIDGET(main_debayer_button)));
	g_signal_connect(GTK_TOGGLE_BUTTON(toggle_debayer), "toggled", G_CALLBACK(on_debayer_toggled), (gpointer) main_debayer_button);
}

static gchar* get_calibration_file_directory(GtkWidget *entry) {
    if (!GTK_IS_ENTRY(GTK_ENTRY(entry))) {
        return NULL;
    }

    // Get the filename from the file chooser button (returns absolute path)
    const gchar *filename = gtk_editable_get_text(GTK_EDITABLE(entry));

    if (!filename || filename[0] == '\0') {
        // No file selected
        return g_strdup(com.wd);
    }

    // Since siril_file_chooser_get_filename() always returns absolute paths,
    // we always extract the directory portion
    gchar *dirname = g_path_get_dirname(filename);
    return dirname;
}

static void opendial(int whichdial) {
	open_dialog_init_statics();
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	fileChooserPreview *preview = NULL;
	GtkWindow *control_window = od_control_window;
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
			entry = GTK_WIDGET(od_flatname_entry);
	case OD_DARK:
		if (whichdial == OD_DARK)
			entry = GTK_WIDGET(od_darkname_entry);
	case OD_OFFSET:
		if (whichdial == OD_OFFSET)
			entry = GTK_WIDGET(od_offsetname_entry);
	case OD_FLATLIB:
		if (whichdial == OD_FLATLIB)
			entry = GTK_WIDGET(od_flatlib_entry);
	case OD_DARKLIB:
		if (whichdial == OD_DARKLIB)
			entry = GTK_WIDGET(od_darklib_entry);
	case OD_OFFSETLIB:
		if (whichdial == OD_OFFSETLIB)
			entry = GTK_WIDGET(od_biaslib_entry);
	case OD_DISTOLIB:
		if (whichdial == OD_DISTOLIB)
			entry = GTK_WIDGET(od_distolib_entry);
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gchar *lastdir = get_calibration_file_directory(entry);
		siril_file_chooser_set_current_folder_path(dialog, lastdir);
		g_free(lastdir);
		/* GTK4: gtk_file_chooser_set_local_only removed */;
		if (whichdial == OD_FLATLIB) {
			if (com.pref.prepro.flat_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.flat_lib);
				siril_file_chooser_set_current_folder_path(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_DARKLIB) {
			if (com.pref.prepro.dark_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.dark_lib);
				siril_file_chooser_set_current_folder_path(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_OFFSETLIB) {
			if (com.pref.prepro.bias_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.bias_lib);
				siril_file_chooser_set_current_folder_path(dialog, path);
				g_free(path);
			}
		} else if (whichdial == OD_DISTOLIB) {
			if (com.pref.prepro.disto_lib != NULL) {
				gchar *path = g_path_get_dirname(com.pref.prepro.disto_lib);
				siril_file_chooser_set_current_folder_path(dialog, path);
				g_free(path);
			}
		G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
		}
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		G_GNUC_END_IGNORE_DEPRECATIONS
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
		siril_file_chooser_set_current_folder_path(dialog, com.wd);
		G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
		/* GTK4: gtk_file_chooser_set_local_only removed */;
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		G_GNUC_END_IGNORE_DEPRECATIONS
		break;
	case OD_OPEN:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		siril_file_chooser_set_current_folder_path(dialog, com.wd);
		G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
		/* GTK4: gtk_file_chooser_set_local_only removed */;
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		G_GNUC_END_IGNORE_DEPRECATIONS
		set_filters_dialog(dialog, whichdial);
		siril_file_chooser_add_preview(dialog, preview);
		siril_add_debayer_toggle_button(dialog);
		break;
	case OD_CONVERT:
		widgetdialog = siril_file_chooser_add(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		siril_file_chooser_set_current_folder_path(dialog, com.wd);
		G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
		/* GTK4: gtk_file_chooser_set_local_only removed */;
		gtk_file_chooser_set_select_multiple(dialog, TRUE);
		G_GNUC_END_IGNORE_DEPRECATIONS
		set_filters_dialog(dialog, whichdial);
		break;
	case OD_BADPIXEL:
		widgetdialog = siril_file_chooser_add(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		siril_file_chooser_set_current_folder_path(dialog, com.wd);
		G_GNUC_BEGIN_IGNORE_DEPRECATIONS  /* Batch C/D pending — see /tmp/deprecation-migration-plan.md */
		/* GTK4: gtk_file_chooser_set_local_only removed */;
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		G_GNUC_END_IGNORE_DEPRECATIONS
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

		filename = siril_file_chooser_get_filename(chooser);
		if (!filename)
			goto wait;

		if (fitseq_is_fitseq(filename, NULL)
				&& ((whichdial == OD_FLAT) || (whichdial == OD_DARK) || (whichdial == OD_OFFSET))) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Format not supported."),
					_("FITS sequences are not supported for master files. Please select a single FITS file."));
			goto wait;
		}

		GtkEntry *flat_entry = od_flatname_entry;
		GtkEntry *dark_entry = od_darkname_entry;
		GtkEntry *bias_entry = od_offsetname_entry;
		GtkEntry *flatlib_entry = od_flatlib_entry;
		GtkEntry *darklib_entry = od_darklib_entry;
		GtkEntry *biaslib_entry = od_biaslib_entry;
		GtkEntry *distolib_entry = od_distolib_entry;
		GtkEntry *bad_pixel_entry = od_pixelmap_entry;
		GtkCheckButton *flat_button = od_useflat_btn;
		GtkCheckButton *dark_button = od_usedark_btn;
		GtkCheckButton *bias_button = od_useoffset_btn;
		GtkWidget *pbutton = od_prepro_button;

		anything_loaded = sequence_is_loaded() || single_image_is_loaded();

		switch (whichdial) {
		case OD_FLAT:
			gtk_editable_set_text(GTK_EDITABLE(flat_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(flat_button), TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_DARK:
			gtk_editable_set_text(GTK_EDITABLE(dark_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(dark_button), TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_OFFSET:
			gtk_editable_set_text(GTK_EDITABLE(bias_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(bias_button), TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_FLATLIB:
			gtk_editable_set_text(GTK_EDITABLE(flatlib_entry), filename);
			break;

		case OD_DARKLIB:
			gtk_editable_set_text(GTK_EDITABLE(darklib_entry), filename);
			break;

		case OD_OFFSETLIB:
			gtk_editable_set_text(GTK_EDITABLE(biaslib_entry), filename);
			break;

		case OD_DISTOLIB:
			gtk_editable_set_text(GTK_EDITABLE(distolib_entry), filename);
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
			gtk_editable_set_text(GTK_EDITABLE(bad_pixel_entry), filename);
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

void on_select_convert_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_CONVERT);
}

void on_pixelmap_button_clicked(GtkToggleButton *button, gpointer user_data) {
	opendial(OD_BADPIXEL);
}

/* Phase 14: GtkRecentChooserMenu is gone, but GtkRecentManager itself
 * still exists in GTK4.  We populate a GMenuModel from it and attach it
 * to the recent-files GtkMenuButton (`recent_menu_button`).  Each item
 * activates `app.open-recent('<path>')`. */
void on_open_recent_action_item_activated(GObject *chooser,
		gpointer user_data) {
	(void)chooser; (void)user_data;
	/* Legacy .ui binding; the GMenu items below drive the open path now. */
}

void open_recent_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	(void)action; (void)user_data;
	if (!parameter) return;
	const gchar *path = g_variant_get_string(parameter, NULL);
	if (!path || !*path) return;
	if (!g_file_test(path, G_FILE_TEST_EXISTS)) {
		siril_log_color_message(_("Recent file no longer exists: %s\n"),
		                        "salmon", path);
		/* Best-effort cleanup of the stale recent entry. */
		gchar *uri = g_filename_to_uri(path, NULL, NULL);
		if (uri) {
			gtk_recent_manager_remove_item(
				gtk_recent_manager_get_default(), uri, NULL);
			g_free(uri);
		}
		return;
	}
	gchar *image_dir = g_path_get_dirname(path);
	if (image_dir) {
		siril_change_dir(image_dir, NULL);
		g_free(image_dir);
	}
	if (!com.script)
		gui_function(set_GUI_CWD, NULL);
	open_single_image(path);
}

static int recent_info_cmp_visited_desc(gconstpointer a, gconstpointer b) {
	GDateTime *ta = gtk_recent_info_get_visited((GtkRecentInfo *) a);
	GDateTime *tb = gtk_recent_info_get_visited((GtkRecentInfo *) b);
	if (!ta && !tb) return 0;
	if (!ta) return 1;
	if (!tb) return -1;
	return g_date_time_compare(tb, ta);  /* descending */
}

/* Build (or rebuild) the recent-files GMenuModel and attach it to the
 * recent_menu_button.  Call once at startup; GtkRecentManager auto-
 * updates its store but we don't subscribe to "changed" yet, so the
 * menu reflects the state at startup. */
void populate_recent_files_menu(void) {
	GtkMenuButton *btn = GTK_MENU_BUTTON(lookup_widget("recent_menu_button"));
	if (!btn) return;
	GtkRecentManager *mgr = gtk_recent_manager_get_default();
	GList *items = g_list_sort(gtk_recent_manager_get_items(mgr),
	                           recent_info_cmp_visited_desc);

	GMenu *menu = g_menu_new();
	int count = 0;
	const int max_items = 12;
	for (GList *l = items; l && count < max_items; l = l->next) {
		GtkRecentInfo *info = l->data;
		const gchar *uri = gtk_recent_info_get_uri(info);
		if (!uri) continue;
		gchar *path = g_filename_from_uri(uri, NULL, NULL);
		if (!path) continue;
		if (!g_file_test(path, G_FILE_TEST_EXISTS)) {
			g_free(path);
			continue;
		}
		/* Restrict to siril-openable image / sequence types so the menu
		 * doesn't fill with unrelated recent files (PDFs, scripts, etc.)
		 * tracked by GtkRecentManager from other apps. */
		if (get_type_from_filename(path) == TYPEUNDEF) {
			g_free(path);
			continue;
		}
		gchar *display = g_path_get_basename(path);
		GMenuItem *item = g_menu_item_new(display, NULL);
		g_menu_item_set_action_and_target_value(item, "app.open-recent",
			g_variant_new_string(path));
		g_menu_append_item(menu, item);
		g_object_unref(item);
		g_free(display);
		g_free(path);
		count++;
	}
	g_list_free_full(items, (GDestroyNotify) gtk_recent_info_unref);

	if (count > 0) {
		gtk_menu_button_set_menu_model(btn, G_MENU_MODEL(menu));
		gtk_widget_set_sensitive(GTK_WIDGET(btn), TRUE);
	} else {
		gtk_menu_button_set_menu_model(btn, NULL);
		gtk_widget_set_sensitive(GTK_WIDGET(btn), FALSE);
	}
	g_object_unref(menu);
}
