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
#include "gui-gtk4/conversion.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/file_browser.h"

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

static void set_filters_browser(SirilFileBrowser *fb, int whichdial) {
	GString *all_filter = NULL;
	gchar *fits_filter = FITS_EXTENSIONS;
	gchar *netpbm_filter = "*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM";
	gchar *pic_filter = "*.pic;*.PIC";
	gchar *ser_filter = "*.ser;*.SER";
	if (whichdial != OD_CONVERT && whichdial != OD_OPEN) {
		siril_file_browser_add_filter_pattern(fb,
				_("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"),
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
			siril_file_browser_add_filter_pattern(fb, _("All supported files"), filter, TRUE);
			g_free(filter);
		}
	}
}

/* GTK4: gtk_file_chooser_set_extra_widget is gone, and SirilFileBrowser
 * has no built-in extra-widget slot — the debayer toggle that used to
 * sit alongside the Open dialog has been dropped here.  Toggling debayer
 * on a per-open basis can still be done from the main demosaic control. */

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

/* Resolve the initial folder for a given calibration master case.  Falls
 * back to com.wd if no specific prior path is recorded. */
static gchar *initial_folder_for(int whichdial, GtkWidget *entry) {
	switch (whichdial) {
	case OD_FLATLIB:
		if (com.pref.prepro.flat_lib && *com.pref.prepro.flat_lib)
			return g_path_get_dirname(com.pref.prepro.flat_lib);
		break;
	case OD_DARKLIB:
		if (com.pref.prepro.dark_lib && *com.pref.prepro.dark_lib)
			return g_path_get_dirname(com.pref.prepro.dark_lib);
		break;
	case OD_OFFSETLIB:
		if (com.pref.prepro.bias_lib && *com.pref.prepro.bias_lib)
			return g_path_get_dirname(com.pref.prepro.bias_lib);
		break;
	case OD_DISTOLIB:
		if (com.pref.prepro.disto_lib && *com.pref.prepro.disto_lib)
			return g_path_get_dirname(com.pref.prepro.disto_lib);
		break;
	}
	return get_calibration_file_directory(entry);
}

/* Image cases — use SirilFileBrowser so the user gets a thumbnail
 * preview of the selected frame.  Returns the picked path (single) or
 * NULL on cancel; populates `out_paths` when whichdial == OD_CONVERT. */
static gchar *pick_image(int whichdial, GtkWindow *parent,
                         GtkWidget *entry, GSList **out_paths) {
	const gchar *title;
	switch (whichdial) {
	case OD_OPEN:      title = _("Open Image"); break;
	case OD_CONVERT:   title = _("Add Files");  break;
	case OD_FLAT:      title = _("Select Flat Frame");      break;
	case OD_DARK:      title = _("Select Dark Frame");      break;
	case OD_OFFSET:    title = _("Select Bias / Offset Frame"); break;
	case OD_FLATLIB:   title = _("Select Master Flat from Library");   break;
	case OD_DARKLIB:   title = _("Select Master Dark from Library");   break;
	case OD_OFFSETLIB: title = _("Select Master Bias from Library");   break;
	default:           title = _("Open Image"); break;
	}

	SirilFileBrowser *fb = siril_file_browser_new(parent, title);
	if (whichdial == OD_CONVERT)
		siril_file_browser_set_select_multiple(fb, TRUE);
	/* Show the Debayer toggle on the cases where the user might actually
	 * want to flip debayer-on-open: single open and Add Files (convert).
	 * The toggle is linked to com.pref.debayer.open_debayer and the
	 * Convert tab's demosaicingButton, so flipping it here changes the
	 * setting consistently across the app. */
	if (whichdial == OD_OPEN || whichdial == OD_CONVERT)
		siril_file_browser_set_show_debayer_toggle(fb, TRUE);
	gchar *initial = initial_folder_for(whichdial, entry);
	if (initial && *initial)
		siril_file_browser_set_initial_folder(fb, initial);
	g_free(initial);
	set_filters_browser(fb, whichdial);

	gchar *picked = NULL;
	if (siril_file_browser_run(fb) == GTK_RESPONSE_ACCEPT) {
		if (whichdial == OD_CONVERT && out_paths)
			*out_paths = siril_file_browser_get_paths(fb);
		picked = siril_file_browser_get_path(fb);
	}
	siril_file_browser_destroy(fb);
	/* Drain pending main-loop events so the NSWindow teardown for the
	 * (modal) browser completes before opendial blocks on the synchronous
	 * load.  On macOS the AppKit run loop defers NSWindow destruction
	 * (and the modal-session release) to the next loop iteration; if we
	 * dive into open_single_image immediately, the dead browser's modal
	 * session stays alive across the whole load and mouse-down events
	 * vanish into nowhere — motion still flows, so X/Y/RGB labels keep
	 * updating, but the rest of the UI looks frozen until opendial
	 * returns.  Linux is unaffected (no CFRunLoop bridge) but the drain
	 * is harmless there. */
	while (g_main_context_pending(NULL))
		g_main_context_iteration(NULL, FALSE);
	return picked;
}

/* Non-image cases — use GtkFileDialog via the SirilFileChooser wrapper. */
static gchar *pick_non_image(int whichdial, GtkWindow *parent) {
	GtkFileChooserAction action =
		(whichdial == OD_CWD) ? GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER
		                      : GTK_FILE_CHOOSER_ACTION_OPEN;
	SirilFileChooser *fc = siril_fc_open(parent, action);
	siril_fc_set_current_folder_path(fc, com.wd);
	if (whichdial == OD_DISTOLIB) {
		gchar *initial = initial_folder_for(OD_DISTOLIB, NULL);
		if (initial) {
			siril_fc_set_current_folder_path(fc, initial);
			g_free(initial);
		}
		siril_fc_add_filter_pattern(fc,
			_("Master-distortion: file (*.wcs)"),
			"*.wcs;*.WCS", TRUE);
	} else if (whichdial == OD_BADPIXEL) {
		siril_fc_add_filter_pattern(fc,
			_("Cosmetic correction file (*.lst)"),
			"*.lst;*.LST", TRUE);
	}
	gchar *picked = NULL;
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT)
		picked = siril_fc_get_filename(fc);
	siril_fc_destroy(fc);
	return picked;
}

static void opendial(int whichdial) {
	open_dialog_init_statics();
	GtkWindow *control_window = od_control_window;

	if (!com.wd)
		return;

	if (whichdial == OD_NULL) {
		fprintf(stderr, "whichdial undefined, should not happen\n");
		return;
	}

	/* Map each calibration case to the entry widget that holds the
	 * "previous selection" path (used to compute the initial folder). */
	GtkWidget *entry = NULL;
	switch (whichdial) {
	case OD_FLAT:      entry = GTK_WIDGET(od_flatname_entry);   break;
	case OD_DARK:      entry = GTK_WIDGET(od_darkname_entry);   break;
	case OD_OFFSET:    entry = GTK_WIDGET(od_offsetname_entry); break;
	case OD_FLATLIB:   entry = GTK_WIDGET(od_flatlib_entry);    break;
	case OD_DARKLIB:   entry = GTK_WIDGET(od_darklib_entry);    break;
	case OD_OFFSETLIB: entry = GTK_WIDGET(od_biaslib_entry);    break;
	case OD_DISTOLIB:  entry = GTK_WIDGET(od_distolib_entry);   break;
	}

	gchar *filename = NULL;
	GSList *multi = NULL;

	for (;;) {
		g_free(filename);
		filename = NULL;
		g_slist_free_full(multi, g_free);
		multi = NULL;

		switch (whichdial) {
		case OD_FLAT:
		case OD_DARK:
		case OD_OFFSET:
		case OD_FLATLIB:
		case OD_DARKLIB:
		case OD_OFFSETLIB:
		case OD_OPEN:
		case OD_CONVERT:
			filename = pick_image(whichdial, control_window, entry,
			                      whichdial == OD_CONVERT ? &multi : NULL);
			break;
		case OD_CWD:
		case OD_DISTOLIB:
		case OD_BADPIXEL:
			filename = pick_non_image(whichdial, control_window);
			break;
		}

		if (!filename) {
			g_slist_free_full(multi, g_free);
			return;  /* user cancelled */
		}

		/* Master-frame validation: reject FITS sequences. */
		if (fitseq_is_fitseq(filename, NULL)
		    && (whichdial == OD_FLAT || whichdial == OD_DARK || whichdial == OD_OFFSET)) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Format not supported."),
				_("FITS sequences are not supported for master files. Please select a single FITS file."));
			continue;
		}

		gchar *err = NULL;
		gboolean anything_loaded = sequence_is_loaded() || single_image_is_loaded();
		int retval;

		switch (whichdial) {
		case OD_FLAT:
			gtk_editable_set_text(GTK_EDITABLE(od_flatname_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(od_useflat_btn), TRUE);
			gtk_widget_set_sensitive(od_prepro_button, anything_loaded);
			break;
		case OD_DARK:
			gtk_editable_set_text(GTK_EDITABLE(od_darkname_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(od_usedark_btn), TRUE);
			gtk_widget_set_sensitive(od_prepro_button, anything_loaded);
			break;
		case OD_OFFSET:
			gtk_editable_set_text(GTK_EDITABLE(od_offsetname_entry), filename);
			siril_toggle_set_active(GTK_WIDGET(od_useoffset_btn), TRUE);
			gtk_widget_set_sensitive(od_prepro_button, anything_loaded);
			break;
		case OD_FLATLIB:
			gtk_editable_set_text(GTK_EDITABLE(od_flatlib_entry), filename);
			break;
		case OD_DARKLIB:
			gtk_editable_set_text(GTK_EDITABLE(od_darklib_entry), filename);
			break;
		case OD_OFFSETLIB:
			gtk_editable_set_text(GTK_EDITABLE(od_biaslib_entry), filename);
			break;
		case OD_DISTOLIB:
			gtk_editable_set_text(GTK_EDITABLE(od_distolib_entry), filename);
			break;
		case OD_CWD:
			if (!siril_change_dir(filename, &err)) {
				if (com.pref.wd) g_free(com.pref.wd);
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
			if (retval == OPEN_IMAGE_CANCEL) continue;  /* re-prompt */
			break;
		case OD_CONVERT:
			multi = g_slist_sort(multi, (GCompareFunc) strcompare);
			fill_convert_list(multi);  /* frees each entry's gchar*; we just free the spine */
			g_slist_free(multi);
			multi = NULL;
			break;
		case OD_BADPIXEL:
			gtk_editable_set_text(GTK_EDITABLE(od_pixelmap_entry), filename);
			break;
		}

		g_free(filename);
		g_slist_free_full(multi, g_free);
		return;
	}
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
		siril_log_warning(_("Recent file no longer exists: %s\n"), path);
		/* Best-effort cleanup of the stale recent entry. */
		gchar *uri = g_filename_to_uri(path, NULL, NULL);
		if (uri) {
			gtk_recent_manager_remove_item(
				gtk_recent_manager_get_default(), uri, NULL);
			g_free(uri);
		}
		return;
	}
	/* The GVariant target is owned by the menu item and will outlive
	 * this handler, but downstream code passes the pointer into idles
	 * that may run after we return — copy onto the heap to be safe. */
	gchar *path_copy = g_strdup(path);
	gchar *image_dir = g_path_get_dirname(path_copy);
	if (image_dir) {
		siril_change_dir(image_dir, NULL);
		g_free(image_dir);
	}
	if (!com.script)
		gui_function(set_GUI_CWD, NULL);
	open_single_image(path_copy);
	g_free(path_copy);
}

static int recent_info_cmp_modified_desc(gconstpointer a, gconstpointer b) {
	/* Sort by "modified" time, not "visited" - gtk_recent_manager_add_full()
	 * stamps modified but leaves visited NULL, so sorting by visited would
	 * push freshly-added entries to the bottom of the 1000-item list. */
	GDateTime *ta = gtk_recent_info_get_modified((GtkRecentInfo *) a);
	GDateTime *tb = gtk_recent_info_get_modified((GtkRecentInfo *) b);
	if (!ta && !tb) return 0;
	if (!ta) return 1;
	if (!tb) return -1;
	return g_date_time_compare(tb, ta);  /* descending */
}

static void recent_manager_changed_cb(GtkRecentManager *mgr, gpointer user_data) {
	(void)mgr; (void)user_data;
	populate_recent_files_menu();
}

/* Build (or rebuild) the recent-files GMenuModel and attach it to the
 * recent_menu_button.  First call subscribes to GtkRecentManager's
 * "changed" signal so subsequent updates (including the async flush
 * after gtk_recent_manager_add_item) auto-refresh the dropdown — the
 * direct populate call from impl_on_image_loaded runs before the store
 * has flushed, so the just-opened file would otherwise not appear until
 * the next manual rebuild. */
void populate_recent_files_menu(void) {
	GtkMenuButton *btn = GTK_MENU_BUTTON(lookup_widget("recent_menu_button"));
	if (!btn) return;
	GtkRecentManager *mgr = gtk_recent_manager_get_default();
	static gboolean changed_subscribed = FALSE;
	if (!changed_subscribed) {
		changed_subscribed = TRUE;
		g_signal_connect(mgr, "changed",
		                 G_CALLBACK(recent_manager_changed_cb), NULL);
	}
	GList *items = g_list_sort(gtk_recent_manager_get_items(mgr),
	                           recent_info_cmp_modified_desc);
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
		/* Restrict to FITS and TIFF only — the formats users actually want
		 * to reopen quickly, and the ones Siril adds via add_full(). */
		image_type t = get_type_from_filename(path);
		if (t != TYPEFITS && t != TYPETIFF) {
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
