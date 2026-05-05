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
#include "core/processing.h"
#include "core/exif.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "io/ser.h"
#include "io/image_format_fits.h"

#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif
#ifdef HAVE_LIBXISF
#include "io/SirilXISFWraper.h"
#endif
#ifdef HAVE_LIBJXL
#include "io/SirilJpegXLWrapper.h"
#endif
#include "dialog_preview.h"
#include "core/proto.h"
#include "filters/mtf.h"

/* Helper and function moved here from io/ser.c because it creates a
 * GdkPixbuf for GUI preview purposes only (sole caller: this file). */
static GdkPixbufDestroyNotify free_preview_data(guchar *pixels, gpointer data) {
	free(pixels);
	free(data);
	return FALSE;
}

GdkPixbuf* get_thumbnail_from_ser(const char *filename, gchar **descr) {
	GdkPixbuf *pixbuf = NULL;
	int MAX_SIZE = com.pref.gui.thumbnail_size;
	gchar *description = NULL;
	int i, j, k, l, N, M;
	int w, h, pixScale, Ws, Hs, n_channels, n_frames, bit;
	int sz;
	struct ser_struct ser;
	fits fit = { 0 };

	ser_init_struct(&ser);
	if (ser_open_file(filename, &ser)) {
		return NULL;
	}
	float *pix = malloc(MAX_SIZE * sizeof(float));
	float *ima_data = NULL, *ptr = NULL, byte, n, max, min, wd, avr;
	guchar *pixbuf_data = NULL;

	w = ser.image_width;
	h = ser.image_height;
	sz = w * h;

	switch (ser.color_id) {
	case SER_MONO:
		n_channels = 1;
		break;
	default:
		n_channels = 3;
	}

	ima_data = malloc(sz * n_channels * sizeof(float));
	pixbuf_data = malloc(3 * MAX_SIZE * MAX_SIZE * sizeof(guchar));

	ser_read_frame(&ser, 0, &fit, FALSE, FALSE);

	if (n_channels == 1) {
		for (i = 0; i < sz; i++) {
			ima_data[i] = (float)fit.pdata[RLAYER][i];
		}
	} else {
		for (i = 0; i < sz; i++) {
			ima_data[i * 3 + 0] = (float)fit.pdata[RLAYER][i];
			ima_data[i * 3 + 1] = (float)fit.pdata[GLAYER][i];
			ima_data[i * 3 + 2] = (float)fit.pdata[BLAYER][i];
		}
	}
	clearfits(&fit);

	i = (int) ceil((float) w / MAX_SIZE);
	j = (int) ceil((float) h / MAX_SIZE);
	pixScale = (i > j) ? i : j;
	if (pixScale == 0) {
		free(ima_data);
		free(pixbuf_data);
		free(pix);
		return NULL;
	}
	Ws = w / pixScale;
	Hs = h / pixScale;

	n_frames = ser.frame_count;
	bit = ser.bit_pixel_depth;

	if (n_channels == 1) {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), bit, n_frames,
				ngettext("frame", "frames", n_frames));
	} else {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), bit, n_frames,
				ngettext("frame", "frames", n_frames));
	}

	float *pix_r = malloc(MAX_SIZE * sizeof(float));
	float *pix_g = malloc(MAX_SIZE * sizeof(float));
	float *pix_b = malloc(MAX_SIZE * sizeof(float));

	M = 0;
	for (i = 0; i < Hs; i++) {
		for (j = 0; j < MAX_SIZE; j++) {
			if (n_channels == 1) {
				pix[j] = 0;
			} else {
				pix_r[j] = 0;
				pix_g[j] = 0;
				pix_b[j] = 0;
			}
		}
		float m = 0.f;
		for (l = 0; l < pixScale; l++, m++) {
			if (n_channels == 1) {
				ptr = &ima_data[M * w];
			} else {
				ptr = &ima_data[M * w * 3];
			}
			N = 0;
			for (j = 0; j < Ws; j++) {
				n = 0.;
				if (n_channels == 1) {
					byte = 0.;
					for (k = 0; k < pixScale; k++, n++) {
						if (N++ < w)
							byte += *ptr++;
						else
							break;
					}
					if (n == 0)
						n = 1.f;
					pix[j] += byte / n;
				} else {
					float byte_r = 0., byte_g = 0., byte_b = 0.;
					for (k = 0; k < pixScale; k++, n++) {
						if (N++ < w) {
							byte_r += *ptr++;
							byte_g += *ptr++;
							byte_b += *ptr++;
						} else {
							break;
						}
					}
					if (n == 0)
						n = 1.f;
					pix_r[j] += byte_r / n;
					pix_g[j] += byte_g / n;
					pix_b[j] += byte_b / n;
				}
			}
			if (++M >= h)
				break;
		}
		if (n_channels == 1) {
			ptr = &ima_data[i * Ws];
			for (l = 0; l < Ws; l++)
				*ptr++ = pix[l] / m;
		} else {
			for (l = 0; l < Ws; l++) {
				ima_data[(i * Ws + l) * 3 + 0] = pix_r[l] / m;
				ima_data[(i * Ws + l) * 3 + 1] = pix_g[l] / m;
				ima_data[(i * Ws + l) * 3 + 2] = pix_b[l] / m;
			}
		}
	}

	if (bit > 8) {
		int prev_size = Ws * Hs;
		float *ima_ch = NULL;
		if (n_channels > 1) {
			ima_ch = malloc(prev_size * n_channels * sizeof(float));
			if (!ima_ch) {
				ima_ch = NULL;
			} else {
				for (i = 0; i < prev_size; i++) {
					ima_ch[0 * prev_size + i] = ima_data[i * 3 + 0];
					ima_ch[1 * prev_size + i] = ima_data[i * 3 + 1];
					ima_ch[2 * prev_size + i] = ima_data[i * 3 + 2];
				}
			}
		} else {
			ima_ch = ima_data;
		}

		if (ima_ch != NULL) {
			float maxmax = 65535.f;
			float minmin = 0.f;
			float scale = 1.f / (maxmax - minmin);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
			for (int idx = 0 ; idx < prev_size * n_channels; idx++) {
				ima_ch[idx] = (ima_ch[idx] - minmin) * scale;
			}
			fits *tmp = NULL;
			new_fit_image_with_data(&tmp, Ws, Hs, n_channels, DATA_FLOAT, ima_ch);
			struct mtf_params mtfp[3] = {
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE }
			};
			find_unlinked_midtones_balance_default(tmp, mtfp);
			apply_unlinked_mtf_to_fits(tmp, tmp, mtfp);
			if (n_channels > 1) {
				for (i = 0; i < prev_size; i++) {
					ima_data[i * 3 + 0] = ima_ch[0 * prev_size + i];
					ima_data[i * 3 + 1] = ima_ch[1 * prev_size + i];
					ima_data[i * 3 + 2] = ima_ch[2 * prev_size + i];
				}
			}
			tmp->fdata = NULL;
			tmp->fpdata[0] = NULL;
			tmp->fpdata[1] = NULL;
			tmp->fpdata[2] = NULL;
			clearfits(tmp);
			free(tmp);
			if (n_channels > 1) {
				free(ima_ch);
			}
		}
	}

	if (n_channels == 1) {
		ptr = ima_data;
		sz = Ws * Hs;
		max = min = *ptr;
		avr = 0;
		for (i = 0; i < sz; i++, ptr++) {
			float tmp = *ptr;
			if (tmp > max) max = tmp;
			else if (tmp < min) min = tmp;
			avr += tmp;
		}
		avr /= (float) sz;
		wd = max - min;
		if (wd == 0.f) wd = 1.f;
		avr = (avr - min) / wd;
		if (avr > 1.)
			wd /= avr;
		ptr = ima_data;
		for (i = Hs - 1; i > -1; i--) {
			guchar *pptr = &pixbuf_data[Ws * i * 3];
			for (j = 0; j < Ws; j++) {
				guchar val = (guchar) roundf_to_BYTE(255.f * (*ptr - min) / wd);
				*pptr++ = val;
				*pptr++ = val;
				*pptr++ = val;
				ptr++;
			}
		}
	} else {
		float max_r = ima_data[0], min_r = ima_data[0];
		float max_g = ima_data[1], min_g = ima_data[1];
		float max_b = ima_data[2], min_b = ima_data[2];
		sz = Ws * Hs;
		for (i = 0; i < sz; i++) {
			float r = ima_data[i * 3 + 0];
			float g = ima_data[i * 3 + 1];
			float b = ima_data[i * 3 + 2];
			if (r > max_r) max_r = r; else if (r < min_r) min_r = r;
			if (g > max_g) max_g = g; else if (g < min_g) min_g = g;
			if (b > max_b) max_b = b; else if (b < min_b) min_b = b;
		}
		float wd_r = max_r - min_r;
		float wd_g = max_g - min_g;
		float wd_b = max_b - min_b;
		if (wd_r == 0.f) wd_r = 1.f;
		if (wd_g == 0.f) wd_g = 1.f;
		if (wd_b == 0.f) wd_b = 1.f;
		for (i = Hs - 1; i > -1; i--) {
			guchar *pptr = &pixbuf_data[Ws * i * 3];
			for (j = 0; j < Ws; j++) {
				int idx = ((Hs - 1 - i) * Ws + j) * 3;
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 0] - min_r) / wd_r);
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 1] - min_g) / wd_g);
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 2] - min_b) / wd_b);
			}
		}
	}

	ser_close_file(&ser);
	pixbuf = gdk_pixbuf_new_from_data(pixbuf_data,
			GDK_COLORSPACE_RGB, FALSE, 8,
			Ws, Hs, Ws * 3,
			(GdkPixbufDestroyNotify) free_preview_data, NULL);
	free(ima_data);
	free(pix);
	free(pix_r);
	free(pix_g);
	free(pix_b);
	*descr = description;
	return pixbuf;
}

static gboolean preview_allocated = FALSE; // flag needed when user load image before preview was displayed.
static gboolean callback_is_called = TRUE;

struct _updta_preview_data {
	GtkFileChooser *file_chooser;
	gchar *filename;
	gchar *description;
	GdkPixbuf *pixbuf;
	GFileInfo *file_info;
	fileChooserPreview *preview;
};

struct _fileChooserPreview {
	GtkWidget *image;
	GtkWidget *name_label;
	GtkWidget *dim_label;
	GtkWidget *size_label;
};

static fileChooserPreview *new_preview_object() {
	fileChooserPreview *object = g_new(fileChooserPreview, 1);
	object->image = gtk_image_new();
	object->name_label = gtk_label_new(NULL);
	object->dim_label = gtk_label_new(NULL);
	object->size_label = gtk_label_new(NULL);
	return object;
}

static gboolean end_update_preview_cb(gpointer p) {
	struct _updta_preview_data *args = (struct _updta_preview_data *) p;

	const char *bytes_str;
	char *size_str = NULL;
	char *name_str = NULL;
	char *info_str = NULL;
	GFileType type;

	fileChooserPreview *preview = args->preview;
	siril_debug_print("preview idle\n");

	if (!preview_allocated || !preview || !(GTK_IS_IMAGE(preview->image))) {
		set_cursor_waiting(FALSE);
		callback_is_called = TRUE;
		return FALSE;
	}

	name_str = g_path_get_basename(args->filename);

	if (!args->file_info) {
		g_free(name_str);
		set_cursor_waiting(FALSE);
		callback_is_called = TRUE;
		return FALSE;
	}
	type = g_file_info_get_file_type(args->file_info);

	/* try to read file size */
	if (args->pixbuf && (bytes_str = gdk_pixbuf_get_option(args->pixbuf, "tEXt::Thumb::Size")) != NULL) {
		int bytes = g_ascii_strtoll(bytes_str, NULL, 10);
		size_str = g_format_size(bytes);
	} else {
		if (type == G_FILE_TYPE_REGULAR) {
			size_str = g_format_size(g_file_info_get_size(args->file_info));
		}
	}

	/* load icon */
	if (type == G_FILE_TYPE_REGULAR && args->pixbuf) {
		gtk_image_set_from_pixbuf(GTK_IMAGE(preview->image), args->pixbuf);
		info_str = args->description;
	} else if (type == G_FILE_TYPE_DIRECTORY) {
		gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "folder", GTK_ICON_SIZE_DIALOG);
		gtk_image_set_pixel_size(GTK_IMAGE(preview->image), com.pref.gui.thumbnail_size);
		info_str = g_strdup(_("Folder"));
	} else {
		image_type im_type = get_type_from_filename(args->filename);
		if (im_type == TYPEAVI || im_type == TYPESER ||
				(im_type == TYPEFITS && fitseq_is_fitseq(args->filename, NULL)))
			gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "video", GTK_ICON_SIZE_DIALOG);
		else {
			gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "image", GTK_ICON_SIZE_DIALOG);
			if (args->description)
				info_str = args->description;
		}
		gtk_image_set_pixel_size(GTK_IMAGE(preview->image), com.pref.gui.thumbnail_size);
	}

	/* information strings */
	const char *format = "<span style=\"italic\">%s</span>";
	char *markup = g_markup_printf_escaped(format, name_str);
	gtk_label_set_markup(GTK_LABEL(preview->name_label), markup);
	gtk_label_set_ellipsize(GTK_LABEL(preview->name_label), PANGO_ELLIPSIZE_MIDDLE);
	gtk_label_set_width_chars(GTK_LABEL(preview->name_label), 25);
	gtk_label_set_max_width_chars(GTK_LABEL(preview->name_label), 25);

	gtk_label_set_text(GTK_LABEL(preview->dim_label), info_str);
	gtk_label_set_text(GTK_LABEL(preview->size_label), size_str);

	if (args->pixbuf)
		g_object_unref(args->pixbuf);
	g_free(markup);
	g_free(name_str);
	g_free(info_str);
	g_free(size_str);

	g_object_unref(args->file_info);
	g_free(args->filename);

	free(args);
	args = NULL;
	set_cursor_waiting(FALSE);
	callback_is_called = TRUE;
	return FALSE;
}

static gpointer update_preview(gpointer p) {
	uint8_t *buffer = NULL;
	size_t size;
	char *mime_type = NULL;
	GdkPixbuf *pixbuf = NULL;
	image_type im_type;
	gboolean libheif_is_ok = FALSE;

	callback_is_called = FALSE;

	struct _updta_preview_data *args = (struct _updta_preview_data *) p;

	args->description = NULL;

	im_type = get_type_from_filename(args->filename);

	if (im_type == TYPEFITS) {
		/* try FITS file */
		pixbuf = get_thumbnail_from_fits(args->filename, &args->description);
	}
#ifdef HAVE_LIBXISF
	else if (im_type == TYPEXISF) {
		GdkPixbuf *pixtmp = get_thumbnail_from_xisf(args->filename, &args->description);

		if (pixtmp) {
			/* The size of the XISF thumbnails is different. We want to resize it */
			const int MAX_SIZE = com.pref.gui.thumbnail_size;

			int w = gdk_pixbuf_get_width(pixtmp);
			int h = gdk_pixbuf_get_height(pixtmp);

			const int x = (int) ceil((float) w / MAX_SIZE);
			const int y = (int) ceil((float) h / MAX_SIZE);
			const int pixScale = (x > y) ? x : y;	// picture scale factor
			const int Ws = w / pixScale; 			// picture width in pixScale blocks
			const int Hs = h / pixScale; 			// -//- height pixScale

			pixbuf = gdk_pixbuf_scale_simple(pixtmp, Ws, Hs, GDK_INTERP_BILINEAR);

			g_object_unref(pixtmp);
		}
	}
#endif
#ifdef HAVE_LIBJXL
	else if (im_type == TYPEJXL) {
		size_t jxl_size = 0;
		uint8_t* jxl_data = NULL;
		GError *error = NULL;
		gboolean success = g_file_get_contents(args->filename, (gchar**) &jxl_data, &jxl_size, &error);
		if (!success) {
			g_error_free(error);
			goto cleanup2;
		}
		siril_debug_print("Generating JXL preview, filesize %lu\n", jxl_size);
		GdkPixbuf *pixtmp = get_thumbnail_from_jxl(jxl_data, &args->description, jxl_size);
		if (pixtmp) {
		/* The size of the JXL thumbnails is different. We want to resize it */
			const int MAX_SIZE = com.pref.gui.thumbnail_size;

			int w = gdk_pixbuf_get_width(pixtmp);
			int h = gdk_pixbuf_get_height(pixtmp);

			const int x = (int) ceil((float) w / MAX_SIZE);
			const int y = (int) ceil((float) h / MAX_SIZE);
			const int pixScale = (x > y) ? x : y;	// picture scale factor
			const int Ws = w / pixScale; 			// picture width in pixScale blocks
			const int Hs = h / pixScale; 			// -//- height pixScale

			pixbuf = gdk_pixbuf_scale_simple(pixtmp, Ws, Hs, GDK_INTERP_BILINEAR);

			g_object_unref(pixtmp);
		}
	}
#endif
	else if (im_type == TYPESER) {
		pixbuf = get_thumbnail_from_ser(args->filename, &args->description);
	} else {
		if (im_type != TYPEUNDEF && !siril_get_thumbnail_exiv(args->filename, &buffer, &size,
				&mime_type)) {
			// Scale the image to the correct size
			GdkPixbuf *tmp;
				GdkPixbufLoader *loader = gdk_pixbuf_loader_new();
			if (!gdk_pixbuf_loader_write(loader, buffer, size, NULL))
				goto cleanup;
			// Calling gdk_pixbuf_loader_close forces the data to be parsed by the
			// loader. We must do this before calling gdk_pixbuf_loader_get_pixbuf.
			if (!gdk_pixbuf_loader_close(loader, NULL))
				goto cleanup;
			if (!(tmp = gdk_pixbuf_loader_get_pixbuf(loader)))
				goto cleanup;
			float ratio = 1.0 * gdk_pixbuf_get_height(tmp) / gdk_pixbuf_get_width(tmp);
			int width = com.pref.gui.thumbnail_size, height = com.pref.gui.thumbnail_size * ratio;
			pixbuf = gdk_pixbuf_scale_simple(tmp, width, height, GDK_INTERP_BILINEAR);
			args->description = siril_get_file_info(args->filename, pixbuf);

			cleanup: gdk_pixbuf_loader_close(loader, NULL);
			free(buffer);
			g_object_unref(loader); // This should clean up tmp as well
		}

		/* if no pixbuf created try to directly read the file */
		/* libheif < 1.6.2 has a bug, therefore we can't open preview if libheif is too old
		 * bug fixed in https://github.com/strukturag/libheif/commit/fbd6d28e8604ecb53a2eb33b522a664b6bcabd0b*/
#ifdef HAVE_LIBHEIF
		libheif_is_ok = LIBHEIF_HAVE_VERSION(1, 6, 2);
#endif

		if (!pixbuf && (im_type != TYPEHEIF || libheif_is_ok)) {
			pixbuf = gdk_pixbuf_new_from_file_at_size(args->filename,
					com.pref.gui.thumbnail_size, com.pref.gui.thumbnail_size, NULL);
			args->description = siril_get_file_info(args->filename, pixbuf);
		}
	}
#ifdef HAVE_LIBJXL
	cleanup2:
#endif
	free(mime_type);
	args->pixbuf = pixbuf;
	gdk_threads_add_idle(end_update_preview_cb, args);
	return GINT_TO_POINTER(0);
}

static void update_preview_cb(GtkFileChooser *file_chooser, gpointer p) {
	gchar *uri;
	GFile *file;
	GFileInfo *file_info;
	fileChooserPreview *preview = (fileChooserPreview *)p;

	uri = gtk_file_chooser_get_preview_uri(file_chooser);
	if (uri == NULL) {
		gtk_file_chooser_set_preview_widget_active(file_chooser, FALSE);
		return;
	}

	file = g_file_new_for_uri(uri);
	file_info = g_file_query_info(file,
				       G_FILE_ATTRIBUTE_TIME_MODIFIED ","
				       G_FILE_ATTRIBUTE_STANDARD_TYPE ","
				       G_FILE_ATTRIBUTE_STANDARD_SIZE ","
				       G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
				       0, NULL, NULL);

	gtk_file_chooser_set_preview_widget_active(file_chooser, TRUE);

	struct _updta_preview_data *data = calloc(1, sizeof(struct _updta_preview_data));
	data->filename = g_file_get_path(file);
	data->file_info = file_info;
	data->file_chooser = file_chooser;
	data->preview = preview;

	g_free(uri);
	g_object_unref(file);

	// this is a graphical operation, we don't use the main processing thread for it, it could block file opening
	g_thread_new("thumbnail", update_preview, data);
}

void siril_preview_free(fileChooserPreview *preview) {
	g_free(preview);
	preview_allocated = FALSE;
}

void siril_file_chooser_add_preview(GtkFileChooser *dialog, fileChooserPreview *preview) {
	if (com.pref.gui.show_thumbnails) {
		GtkWidget *vbox;

		vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
		gtk_container_set_border_width(GTK_CONTAINER(vbox), 12);

		preview = new_preview_object();
		preview_allocated = TRUE;

		gtk_label_set_justify(GTK_LABEL(preview->name_label), GTK_JUSTIFY_CENTER);
		gtk_label_set_justify(GTK_LABEL(preview->dim_label), GTK_JUSTIFY_CENTER);
		gtk_label_set_justify(GTK_LABEL(preview->dim_label), GTK_JUSTIFY_CENTER);

		gtk_widget_set_size_request(preview->image, com.pref.gui.thumbnail_size, com.pref.gui.thumbnail_size);

		gtk_box_pack_start(GTK_BOX(vbox), preview->image, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview->name_label, FALSE, TRUE, 10);
		gtk_box_pack_start(GTK_BOX(vbox), preview->size_label, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview->dim_label, FALSE, TRUE, 0);

		gtk_widget_show_all(vbox);

		gtk_file_chooser_set_preview_widget(dialog, vbox);
		gtk_file_chooser_set_use_preview_label(dialog, FALSE);
		gtk_file_chooser_set_preview_widget_active(dialog, FALSE);

		g_signal_connect(dialog, "update-preview", G_CALLBACK(update_preview_cb), (gpointer)preview);
	}
}

gboolean is_callback_called() {
	return callback_is_called;
}
