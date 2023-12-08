/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/dialog_preview.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/Astro-TIFF.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "save_dialog.h"

static image_type type_of_image = TYPEUNDEF;
static SirilWidget *saveDialog = NULL;

static void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern, gboolean set_default) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	gchar **patterns = g_strsplit(pattern, ";", -1);
	for (gint i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	if (set_default)
		gtk_file_chooser_set_filter(file_chooser, f);
}

static void set_filters_save_dialog(GtkFileChooser *chooser) {
	GString *all_filter = NULL;
	const gchar *fits_filter = "*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.FIT.fz;*.fits.fz;*.FITS.fz;*.fts.fz;*.FTS.fz";
	const gchar *bmp_filter = "*.bmp;*.BMP";
	const gchar *netpbm_filter = "*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM";

	all_filter = g_string_new(fits_filter);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, netpbm_filter);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, bmp_filter);

	gtk_filter_add(chooser, _("FITS Files (*.fit, *.fits, *.fts)"), fits_filter, FALSE);

	/* GRAPHICS FILES */
	/* BMP FILES */
	gtk_filter_add(chooser, _("BMP Files (*.bmp)"), bmp_filter, FALSE);
#ifdef HAVE_LIBJPEG
	const gchar *jpg_filter = "*.jpg;*.JPG;*.jpeg;*.JPEG";

	gtk_filter_add(chooser, _("JPEG Files (*.jpg, *.jpeg)"), jpg_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, jpg_filter);
#endif

#ifdef HAVE_LIBJXL
	const gchar *jxl_filter = "*.jxl;*.JXL";

	gtk_filter_add(chooser, _("JPEG XL Files (*.jxl)"), jxl_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, jxl_filter);
#endif
#ifdef HAVE_LIBHEIF
	const gchar *heifavif_filter = "*.heif;*.HEIF;*.heic;*.HEIC;*.avif;*.AVIF";

	gtk_filter_add(chooser, _("AVIF & HEIF Files (*.avif, *.heif, *.heic)"), heifavif_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, heifavif_filter);
#endif

#ifdef HAVE_LIBPNG
	const gchar *png_filter = "*.png;*.PNG";

	gtk_filter_add(chooser, _("PNG Files (*.png)"), png_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, png_filter);
#endif

#ifdef HAVE_LIBTIFF
	const gchar *tif_filter = "*.tif;*.TIF;*.tiff;*.TIFF";

	gtk_filter_add(chooser, _("TIFF Files (*.tif, *.tiff)"), tif_filter, FALSE);
	all_filter = g_string_append(all_filter, ";");
	all_filter = g_string_append(all_filter, tif_filter);
#endif

	/* NETPBM FILES */
	gtk_filter_add(chooser, _("Netpbm Files (*.ppm, *.pnm, *.pgm)"), netpbm_filter, FALSE);

	gchar *supported_file = g_string_free(all_filter, FALSE);
	gtk_filter_add(chooser, _("Supported Image Files"), supported_file, TRUE);

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
		} else if (!g_strcmp0(string[i], "heif") || !g_strcmp0(string[i], "heic")) {
			type = TYPEHEIF;
			break;
		} else if (!g_strcmp0(string[i], "avif")) {
			type = TYPEAVIF;
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
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));

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
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Description_txt"));

	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(TIFF_txt);

	gtk_text_buffer_get_bounds(tbuf, &itStart, &itEnd);
	gtk_text_buffer_delete(tbuf, &itStart, &itEnd);

	if (gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_type_of_tiff"))) == 0) {
		gchar *astro_tiff = AstroTiff_build_header(&gfit);

		gtk_text_buffer_get_end_iter(tbuf, &itEnd);
		gtk_text_buffer_insert(tbuf, &itEnd, astro_tiff, -1);
		gtk_text_buffer_get_end_iter(tbuf, &itEnd);
		gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
		g_free(astro_tiff);
	} else {
		/* History already written in header */
		if (gfit.history) {
			GSList *list;
			for (list = gfit.history; list; list = list->next) {
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, (gchar *)list->data, -1);
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
			}
		}
		/* New history */
		if (com.history) {
			for (int i = 0; i < com.hist_display; i++) {
				if (com.history[i].history[0] != '\0') {
					gtk_text_buffer_get_end_iter(tbuf, &itEnd);
					gtk_text_buffer_insert(tbuf, &itEnd, com.history[i].history, strlen(com.history[i].history));
					gtk_text_buffer_get_end_iter(tbuf, &itEnd);
					gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
				}
			}
		}
	}
}

static void prepare_savepopup() {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	static GtkWidget *savetxt = NULL;
	int tab;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(lookup_widget("notebookFormat"));
		savepopup = lookup_widget("savepopup");
		savetxt = lookup_widget("filenameframe");
	}

	GtkWindow *parent = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));

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
	case TYPEAVIF:
	case TYPEHEIF:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving AVIF / HEIF"));
		tab = PAGE_AVIF;
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
	gtk_notebook_set_current_page(notebookFormat, tab);
}

static void init_dialog() {
	if (saveDialog == NULL) {
		GtkWindow *parent = siril_get_active_window();
		saveDialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	}
}

static void close_dialog() {
	if (saveDialog != NULL) {
		siril_widget_destroy(saveDialog);
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

static void filter_changed(gpointer user_data) {
	GtkFileChooser *chooser = GTK_FILE_CHOOSER(user_data);
	GtkFileFilter *filter = gtk_file_chooser_get_filter(chooser);
	gchar *filename = siril_file_chooser_get_filename(chooser);
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
		if (gfit.naxes[2] == 1) {
			new_filename = g_strdup_printf("%s.pgm", file_no_ext);
		} else {
			new_filename = g_strdup_printf("%s.ppm", file_no_ext);
		}
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
#ifdef HAVE_LIBHEIF
	case TYPEAVIF:
		new_filename = g_strdup_printf("%s.avif", file_no_ext);
		break;
	case TYPEHEIF:
		new_filename = g_strdup_printf("%s.heif", file_no_ext);
		break;
#endif
#ifdef HAVE_LIBPNG
	case TYPEPNG:
		new_filename = g_strdup_printf("%s.png", file_no_ext);
		break;
#endif
	}
	if (new_filename) {
		gchar *bname = g_path_get_basename(new_filename);
		gtk_file_chooser_set_current_name(chooser, bname);

		g_free(bname);
		g_free(new_filename);
	}

	g_free(file_no_ext);
	g_free(filename);
}

static int save_dialog() {
	init_dialog();

	GtkFileChooser *chooser = GTK_FILE_CHOOSER(saveDialog);
	fileChooserPreview *preview = NULL;
	gchar *fname = get_filename_and_replace_ext();

	gtk_file_chooser_set_current_name(chooser, fname);
	gtk_file_chooser_set_do_overwrite_confirmation(chooser, TRUE);
	gtk_file_chooser_set_local_only(chooser, FALSE);
	set_filters_save_dialog(chooser);
	siril_file_chooser_add_preview(chooser, preview);
	gtk_file_chooser_unselect_all(chooser);
	g_signal_connect(chooser, "notify::filter", G_CALLBACK(filter_changed), (gpointer) chooser);
	g_free(fname);

	gint res = siril_dialog_run(saveDialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *filename = siril_file_chooser_get_filename(chooser);
		type_of_image = get_type_from_filename(filename);

		GtkEntry *savetext = GTK_ENTRY(lookup_widget("savetxt"));
		gtk_entry_set_text(savetext, filename);

		prepare_savepopup();

		g_free(filename);
		return res;
	}
	close_dialog();
	siril_preview_free(preview);

	return res;
}

// idle function executed at the end of the Save Data processing
static gboolean end_save(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
	stop_processing_thread();

	if (args->retval)
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("File saving failed. Check the logs for more info."));

	gtk_entry_set_text(args->entry, "");
	gtk_widget_hide(lookup_widget("savepopup"));
	display_filename(); // update filename display
	set_precision_switch();
	set_cursor_waiting(FALSE);
	close_dialog();	// is this different from the hide above?
	update_MenuItem();

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

static gboolean initialize_data(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;

	GtkToggleButton *fits_8 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *fits_16s = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16s"));
	GtkToggleButton *fits_16 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16"));
	GtkToggleButton *update_hilo = (GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_update_hilo")));
#ifdef HAVE_LIBJPEG
	GtkSpinButton *qlty_spin_button = GTK_SPIN_BUTTON(lookup_widget("quality_spinbutton"));
	args->quality = gtk_spin_button_get_value_as_int(qlty_spin_button);
#endif
#ifdef HAVE_LIBJXL
	GtkSpinButton *quality_spin_button = GTK_SPIN_BUTTON(lookup_widget("jxl_quality_spinbutton"));
	args->jxl_quality = gtk_spin_button_get_value_as_int(quality_spin_button);
	GtkSpinButton *effort_spin_button = GTK_SPIN_BUTTON(lookup_widget("jxl_effort_spinbutton"));
	args->jxl_effort = gtk_spin_button_get_value_as_int(effort_spin_button);
	GtkToggleButton *toggle_button_8bit = GTK_TOGGLE_BUTTON(lookup_widget("jxl_force_8bit"));
	args->jxl_force_8bit = gtk_toggle_button_get_active(toggle_button_8bit);
#endif
#ifdef HAVE_LIBHEIF
	GtkSpinButton *avif_quality_spin_button = GTK_SPIN_BUTTON(lookup_widget("avif_quality_spinbutton"));
	args->heif_quality = gtk_spin_button_get_value_as_int(avif_quality_spin_button);
	GtkComboBox *avif_bitdepth = GTK_COMBO_BOX(lookup_widget("avif_bitdepth"));
	int depth = gtk_combo_box_get_active(avif_bitdepth);
	args->avif_bitdepth = depth == 0 ? 8 : depth == 1 ? 10 : 12;
	GtkToggleButton *toggle_button_lossless = GTK_TOGGLE_BUTTON(lookup_widget("avif_lossless"));
	args->lossless = gtk_toggle_button_get_active(toggle_button_lossless);
#endif
#ifdef HAVE_LIBTIFF
	GtkToggleButton *button_8 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton8bits"));
	GtkToggleButton *button_32 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton32bits"));
	args->bitspersamples = gtk_toggle_button_get_active(button_8) ? 8 : gtk_toggle_button_get_active(button_32) ? 32 : 16;
	get_tif_data_from_ui(&gfit, &args->description, &args->copyright);
	args->tiff_compression = get_tiff_compression();
#endif
	args->entry = GTK_ENTRY(lookup_widget("savetxt"));
	args->filename = gtk_entry_get_text(args->entry);

	if (gtk_toggle_button_get_active(fits_8))
		args->bitpix = BYTE_IMG;
	else if (gtk_toggle_button_get_active(fits_16s))
		args->bitpix = SHORT_IMG;
	else if (gtk_toggle_button_get_active(fits_16))
		args->bitpix = USHORT_IMG;
	else
		args->bitpix = FLOAT_IMG;

	args->update_hilo = gtk_toggle_button_get_active(update_hilo);

	return test_for_viewer_mode();
}

static gpointer mini_save_dialog(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
#ifdef HAVE_LIBPNG
	uint32_t bytes_per_sample;
#endif
	args->retval = 0;

	if (args->filename[0] != '\0') {
		switch (type_of_image) {
		case TYPEBMP:
			args->retval = savebmp(args->filename, &gfit);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
			args->retval = savejpg(args->filename, &gfit, args->quality);
			break;
#endif
#ifdef HAVE_LIBJXL
		case TYPEJXL:
			args->retval = savejxl(args->filename, &gfit, args->jxl_effort, args->jxl_quality, args->jxl_force_8bit);
			break;
#endif
#ifdef HAVE_LIBHEIF
		case TYPEAVIF:
			args->retval = saveheifavif(args->filename, &gfit, args->heif_quality, args->lossless, TRUE, args->avif_bitdepth);
			break;
		case TYPEHEIF:
			args->retval = saveheifavif(args->filename, &gfit, args->heif_quality, args->lossless, FALSE, args->avif_bitdepth);
			break;
#endif
			#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			args->retval = savetif(args->filename, &gfit, args->bitspersamples, args->description, args->copyright, args->tiff_compression, TRUE, TRUE);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
			args->retval = savepng(args->filename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
			break;
#endif
		default:
			type_of_image = TYPEFITS;
			/* no break */
		case TYPEFITS:
			gfit.bitpix = args->bitpix;
			/* Check if MIPS-HI and MIPS-LO must be updated. If yes,
			 * Values are taken from the layer 0 */
			if (args->update_hilo) {
				gfit.hi = gui.hi;
				gfit.lo = gui.lo;
				if (gfit.orig_bitpix == BYTE_IMG
						&& (gfit.hi > UCHAR_MAX || gfit.lo > UCHAR_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				} else if (gfit.orig_bitpix == SHORT_IMG
						&& (gfit.hi > SHRT_MAX || gfit.lo > SHRT_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				}
				if (gfit.orig_bitpix == BYTE_IMG && gfit.bitpix != BYTE_IMG) {
					gfit.hi = USHRT_MAX;
					gfit.lo = 0;
				}
			} else {
				gfit.hi = 0;
				gfit.lo = 0;
			}
			args->retval = savefits(args->filename, &gfit);
			if (!args->retval && single_image_is_loaded()) {
				com.uniq->filename = strdup(args->filename);
				com.uniq->fileexist = TRUE;
			}
			break;
		case TYPEPNM:
			args->retval = saveNetPBM(args->filename, &gfit);
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
		GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		char *file_no_ext = remove_ext_from_filename(filename);
		gtk_entry_set_text(entry, file_no_ext);
		free(file_no_ext);
	}
}

void on_savetxt_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(editable);
	GtkWidget *button = lookup_widget("button_savepopup");

	const gchar *name = gtk_entry_get_text(entry);
	gtk_widget_set_sensitive(button, (*name != '\0'));
}

void on_button_savepopup_clicked(GtkButton *button, gpointer user_data) {
	struct savedial_data *args = calloc(1, sizeof(struct savedial_data));

	set_cursor_waiting(TRUE);
	if (initialize_data(args)) {
		start_in_new_thread(mini_save_dialog, args);
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
	if (initialize_data(args)) {
		start_in_new_thread(mini_save_dialog, args);
	} else {
		g_free(args->copyright);
		g_free(args->description);
		free(args);
		siril_add_idle(end_generic, NULL);
	}
}

void on_button_cancelpopup_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("savepopup"));
}

void on_header_save_as_button_clicked() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		GtkWidget *savepopup = lookup_widget("savepopup");

		if (save_dialog() == GTK_RESPONSE_ACCEPT) {
			/* now it is not needed for some formats */
			if (type_of_image & (TYPEBMP | TYPEPNG | TYPEPNM)) {
				struct savedial_data *args = calloc(1, sizeof(struct savedial_data));

				set_cursor_waiting(TRUE);
				if (initialize_data(args)) {
					start_in_new_thread(mini_save_dialog, args);
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
					parent = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
				}
				gtk_window_set_transient_for(GTK_WINDOW(savepopup), parent);
				gtk_widget_show(savepopup);
				gtk_window_present_with_time(GTK_WINDOW(savepopup),	GDK_CURRENT_TIME);
			}
		}
	}
}

static gboolean snapshot_notification_close(gpointer user_data) {
	gtk_widget_hide(GTK_WIDGET(user_data));
	return FALSE;
}

static GtkWidget *snapshot_notification(GtkWidget *widget, const gchar *filename, GdkPixbuf *pixbuf) {
	gchar *text;
	if (filename)
		text = g_strdup_printf(_("Snapshot <b>%s</b> was saved into the working directory."), filename);
	else
		text = g_strdup((_("Snapshot was saved into the clipboard.")));
	GtkWidget *popover = popover_new_with_image(widget, text, pixbuf);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	g_free(text);
	return popover;
}

static void snapshot_callback(GObject *source_object, GAsyncResult *result,
		gpointer user_data) {
	GError *error = NULL;

	if (!gdk_pixbuf_save_to_stream_finish(result, &error)) {
		siril_log_message(_("Cannot take snapshot: %s\n"), error->message);
		g_clear_error(&error);
	} else {
		gchar *filename = (gchar *)user_data;
		GtkWidget *widget = lookup_widget("header_snapshot_button");
		GtkWidget *popover = snapshot_notification(widget, filename, (GdkPixbuf *)source_object);
		g_timeout_add(5000, (GSourceFunc) snapshot_notification_close, (gpointer) popover);

		g_free(filename);
	}
}

void on_header_snapshot_button_clicked(gboolean clipboard) {
	GError *error = NULL;
	gchar *timestamp, *filename;
	GtkWidget *widget;
	const gchar *area[] = {"drawingarear", "drawingareag", "drawingareab", "drawingareargb" };

	widget = lookup_widget(area[gui.cvport]);
	timestamp = build_timestamp_filename();
	filename = g_strdup_printf("%s.png", timestamp);

	g_free(timestamp);
	/* create cr from the surface */
	int width = max(gfit.rx, gtk_widget_get_allocated_width(widget));
	int height = max(gfit.ry, gtk_widget_get_allocated_height(widget));
	cairo_surface_t *surface = cairo_surface_create_similar(gui.view[gui.cvport].full_surface, CAIRO_CONTENT_COLOR_ALPHA, width, height);
	cairo_t *cr = cairo_create(surface);

	/* add image and all annotations if available */
	add_image_and_label_to_cairo(cr, gui.cvport);

	cairo_destroy(cr);

	double z = get_zoom_val();

	gint x1 = max(0, (int) gui.display_offset.x);
	gint y1 = max(0, (int) gui.display_offset.y);

	gint x2 = min(width * z + (int) gui.display_offset.x, gtk_widget_get_allocated_width(widget));
	gint y2 = min(height * z + (int) gui.display_offset.y, gtk_widget_get_allocated_height(widget));

	GdkPixbuf *pixbuf = gdk_pixbuf_get_from_surface(surface, x1, y1, x2 - x1, y2 - y1);
	if (pixbuf) {
		if (clipboard) {
			GtkClipboard *cb = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
			gtk_clipboard_set_image(cb, pixbuf);
			gtk_clipboard_store(cb);

			GtkWidget *w = lookup_widget("header_snapshot_button");
			GtkWidget *popover = snapshot_notification(w, NULL, pixbuf);
			g_timeout_add(5000, (GSourceFunc) snapshot_notification_close, (gpointer) popover);
		} else {
			GFile *file = g_file_new_build_filename(com.wd, filename, NULL);

			GOutputStream *stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
			if (stream == NULL) {
				if (error != NULL) {
					g_warning("%s\n", error->message);
					g_clear_error(&error);
				}
				cairo_surface_destroy(surface);
				g_free(filename);
				g_object_unref(pixbuf);
				g_object_unref(file);
				return;
			}
			gdk_pixbuf_save_to_stream_async(pixbuf, stream, "png", NULL, snapshot_callback, (gpointer) g_file_get_basename(file), NULL);

			g_object_unref(stream);
			g_object_unref(file);
		}
		g_object_unref(pixbuf);
	}
	cairo_surface_destroy(surface);
	g_free(filename);
}

#undef NEW_SIZE

void on_header_save_button_clicked() {
	if (single_image_is_loaded() && com.uniq->fileexist) {
		savefits(com.uniq->filename, &gfit);
	}
}

void on_savepopup_show(GtkWidget *widget, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scrolledwindow3"));
	gint height, width;

	if (type_of_image & TYPETIFF) {
		GtkToggleButton *b16 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton16bits"));
		GtkToggleButton *b32 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton32bits"));
		gtk_toggle_button_set_active(b32, gfit.type == DATA_FLOAT);
		gtk_toggle_button_set_active(b16, gfit.type == DATA_USHORT);
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
		GtkToggleButton *b16bitu = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16"));
		GtkToggleButton *b32bits = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit32f"));
		gtk_toggle_button_set_active(b32bits, gfit.type == DATA_FLOAT);
		gtk_toggle_button_set_active(b16bitu, gfit.type == DATA_USHORT);
	}

	gtk_scrolled_window_set_min_content_height(scrolled_window, height);
	gtk_scrolled_window_set_min_content_width(scrolled_window, width);
}

void on_combo_type_of_tiff_changed(GtkComboBox* box, gpointer user_data) {
	set_description_in_TIFF();
}
