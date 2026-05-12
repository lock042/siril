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
#include <glib.h>

#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "core/siril_log.h"
#include "utils.h"
#include "message_dialog.h"

struct _label_data {
	const char *label_name;
	char *text;
	const char *class_to_add;
	const char *class_to_remove;
};

static GMutex gui_mutex = { 0 };

void gui_mutex_lock() {
	g_mutex_lock(&gui_mutex);
}

void gui_mutex_unlock() {
	g_mutex_unlock(&gui_mutex);
}

static gboolean set_label_text_idle(gpointer p) {
	struct _label_data *args = (struct _label_data *) p;
	GtkWidget *widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, args->label_name));
	GtkLabel *label = GTK_LABEL(widget);

	if (args->class_to_add || args->class_to_remove) {
		GtkStyleContext *context = gtk_widget_get_style_context(widget);
		if (args->class_to_add)
			gtk_style_context_add_class(context, args->class_to_add);
		if (args->class_to_remove)
			gtk_style_context_remove_class(context, args->class_to_remove);
	}
	if (strcmp(gtk_label_get_text(label), args->text))
		gtk_label_set_text(label, args->text);
	free(args->text);
	free(args);
	return FALSE;
}

static void set_label_text_from_main_thread(const char *label_name, const char *text, const char *class_to_add, const char *class_to_remove) {
	struct _label_data *data = malloc(sizeof(struct _label_data));
	data->label_name = label_name;
	data->text = strdup(text);
	data->class_to_add = class_to_add;
	data->class_to_remove = class_to_remove;
	gdk_threads_add_idle(set_label_text_idle, data);
}

void widget_set_class(GtkWidget *entry, const char *class_to_add, const char *class_to_remove) {
	GtkStyleContext *context = gtk_widget_get_style_context(entry);
	if (class_to_add)
		gtk_style_context_add_class(context, class_to_add);
	if (class_to_remove)
		gtk_style_context_remove_class(context, class_to_remove);
}

GtkWidget* lookup_widget(const gchar *widget_name) {
	return GTK_WIDGET(gtk_builder_get_object(gui.builder, widget_name));
}

GObject* lookup_gobject(const gchar *gobject_name) {
	return ((GObject*) gtk_builder_get_object(gui.builder, gobject_name));
}

GtkAdjustment* lookup_adjustment(const gchar *adjustment_name) {
	return GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, adjustment_name));
}

static GtkNotebook *utils_notebook_center_box = NULL;

static gboolean switch_tab(gpointer user_data) {
	main_tabs tab = (main_tabs) GPOINTER_TO_INT(user_data);
	if (!utils_notebook_center_box)
		utils_notebook_center_box = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_center_box"));
	gtk_notebook_set_current_page(utils_notebook_center_box, tab);
	return FALSE;
}

void control_window_switch_to_tab(main_tabs tab) {
	if (com.script || com.headless)
		return;
	gui_function(switch_tab, GINT_TO_POINTER(tab));
}

/**
 * Create a popover with icon and text
 * @param widget is the parent widget where the popover arises from
 * @param text will be shown in the popover
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new(GtkWidget *widget, const gchar *text) {
	return popover_new_with_image(widget, text, NULL);
}

/**
 * Create a popover with icon and text
 * @param widget is the parent widget where the popover arises from
 * @param text will be shown in the popover
 * @param pixbuf will be shown in the popover
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPixbuf *pixbuf) {
	GtkWidget *popover, *box, *image, *label;

	popover = gtk_popover_new(widget);
	label = gtk_label_new(NULL);
	gtk_label_set_line_wrap_mode(GTK_LABEL(label), PANGO_WRAP_WORD);
	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);

	if (pixbuf) {
		float ratio = 1.0 * gdk_pixbuf_get_height(pixbuf) / gdk_pixbuf_get_width(pixbuf);
		int width = 128, height = 128 * ratio;
		GdkPixbuf *new_pixbuf = gdk_pixbuf_scale_simple(pixbuf, width, height, GDK_INTERP_BILINEAR);
		image = gtk_image_new_from_pixbuf(new_pixbuf);
		g_object_unref(new_pixbuf);
	} else {
		image = gtk_image_new_from_icon_name("dialog-information-symbolic", GTK_ICON_SIZE_DIALOG);
	}

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 100);

	gtk_box_pack_start(GTK_BOX(box), image, FALSE, FALSE, 10);
	gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 10);
	gtk_container_add(GTK_CONTAINER(popover), box);

	/* make all sensitive in case where parent is not */
	gtk_widget_set_sensitive(label, TRUE);
	gtk_widget_set_sensitive(box, TRUE);
	gtk_widget_set_sensitive(popover, TRUE);

	gtk_widget_show_all(box);

	return popover;
}

GList *get_row_references_of_selected_rows(GtkTreeSelection *selection,
		GtkTreeModel *model) {
	GList *ref = NULL;
	GList *sel, *s;

	sel = gtk_tree_selection_get_selected_rows(selection, &model);

	for (s = sel; s; s = s->next) {
		GtkTreeRowReference *rowref = gtk_tree_row_reference_new(model,	(GtkTreePath *) s->data);
		ref = g_list_prepend(ref, rowref);
	}
	g_list_free_full(sel, (GDestroyNotify) gtk_tree_path_free);
	return ref;
}

/* size is in Bytes */
void set_GUI_MEM(guint64 used, const gchar *label) {
	if (com.headless)
		return;
	gchar *str;
	if (used > 0) {
		gchar *mem = g_format_size_full(used, G_FORMAT_SIZE_IEC_UNITS);
		str = g_strdup_printf(_("Mem: %s"), mem);
		g_free(mem);
	} else {
		str = g_strdup(_("Mem: N/A"));
	}
	set_label_text_from_main_thread(label, str, NULL, NULL);
	g_free(str);
}

void set_GUI_DiskSpace(gint64 space, const gchar *label) {
	if (com.headless)
		return;
	gchar *str;
	gboolean set_class = FALSE;

	if (space > 0) {
		if (space < 1073741824) { // we want to warn user if space is less than 1GiB
			set_class = TRUE;
		}
		gchar *mem = g_format_size_full(space, G_FORMAT_SIZE_IEC_UNITS);
		str = g_strdup_printf(_("Disk Space: %s"), mem);
		g_free(mem);
	} else {
		set_class = TRUE;
		str = g_strdup(_("Disk Space: N/A"));
	}
	set_label_text_from_main_thread(label, str,
			set_class ? "label-info" : NULL,
			set_class ? NULL : "label-info");
	g_free(str);
}

void set_suggested(GtkWidget *widget) {
	gtk_style_context_add_class(gtk_widget_get_style_context(widget),
			GTK_STYLE_CLASS_SUGGESTED_ACTION);
}

void unset_suggested(GtkWidget *widget) {
	gtk_style_context_remove_class(gtk_widget_get_style_context(widget),
			GTK_STYLE_CLASS_SUGGESTED_ACTION);
}

// sets the nth button of a container highlighted (useful for stack switchers)
// the Glist can be populated with gtk_container_get_children()
void set_switcher_buttons_colors(GList *list, int n) {
	GList *l = list;
	int i = 0;
	while (l != NULL) {
		GtkWidget *button = (GtkWidget *)l->data;
		if (i == n) {
			set_suggested(button);
		} else {
			unset_suggested(button);
		}
		l = l->next;
		i++;
	}
}

/* Managing GUI calls *
 * We try to separate core concerns from GUI (GTK+) concerns to keep the core
 * code efficient and clean, and avoid calling GTK+ functions in threads other
 * than the main GTK+ thread.
 * GTK+ functions calls should be grouped in idle functions: these are
 * functions called by the GTK+ main thread whenever it's idle. The problem is
 * that is we queue several core functions, like in a script execution, the
 * idle functions may not execute for a purpose corresponding to the current
 * core data. To fix this issue, when siril is executed in GUI mode, we
 * synchronize the GUI calls with the core calls with the following functions.
 * This will slow the execution of some commands down, but will ensure that all
 * is fine in the GUI.
 */

struct idle_data {
	gboolean idle_finished;
	GMutex mutex;
	GCond cond;
	gboolean (*idle)(gpointer);
	gpointer user;
};

static gboolean wrapping_idle(gpointer arg) {
	struct idle_data *data = (struct idle_data *)arg;
	fprintf(stderr, "Entering wrapping_idle for function %p\n", data->idle);
	data->idle(data->user);
	siril_debug_print("idle %p signaling end\n", data->idle);
	g_mutex_lock(&data->mutex);
	data->idle_finished = TRUE;
	g_cond_signal(&data->cond);
	g_mutex_unlock(&data->mutex);
	return FALSE;
}

void execute_idle_and_wait_for_it(gboolean (* idle)(gpointer), gpointer arg) {
	if (com.headless) {
		siril_debug_print("execute_idle_and_wait_for_it called headless, this should not happen!\n");
		return;
	}

	if (g_main_context_is_owner(g_main_context_default())) {
		idle(arg);
		return;
	}

	struct idle_data *data = g_malloc(sizeof(struct idle_data));
	data->idle_finished = FALSE;
	data->idle = idle;
	data->user = arg;
	g_mutex_init(&data->mutex);
	g_cond_init(&data->cond);

	siril_debug_print("queueing idle %p\n", idle);
	gdk_threads_add_idle(wrapping_idle, data);
	siril_debug_print("waiting for idle %p\n", idle);

	/* Wait for the idle to run.  Simple g_cond_wait() deadlocks when the
	 * calling thread is the main thread but is not currently inside a main
	 * loop dispatch (e.g. during shutdown cleanup): nobody pumps the main
	 * context so wrapping_idle never executes.
	 *
	 * Solution: poll with a short timeout.  On each timeout, try to acquire
	 * the main context ourselves.  If we succeed (nobody else is iterating),
	 * pump it once so the pending idle dispatches, then release.  We must
	 * drop our own mutex before pumping, because wrapping_idle will try to
	 * lock it to signal completion. */
	g_mutex_lock(&data->mutex);
	while (!data->idle_finished) {
		gint64 deadline = g_get_monotonic_time() + G_TIME_SPAN_MILLISECOND;
		if (!g_cond_wait_until(&data->cond, &data->mutex, deadline)) {
			/* Timed out — nobody delivered the idle yet.  Try to pump. */
			if (g_main_context_acquire(g_main_context_default())) {
				g_mutex_unlock(&data->mutex);
				g_main_context_iteration(g_main_context_default(), FALSE);
				g_main_context_release(g_main_context_default());
				g_mutex_lock(&data->mutex);
			}
		}
	}
	g_mutex_unlock(&data->mutex);

	g_mutex_clear(&data->mutex);
	g_cond_clear(&data->cond);
	g_free(data);
	siril_debug_print("idle %p wait is over\n", idle);
}

int select_vport(int vport) {
	return vport == RGB_VPORT ? GREEN_VPORT : vport;
}

gboolean check_ok_if_cfa() {
	gboolean retval;
	if (gfit->naxes[2] == 1 && gfit->keywords.bayer_pattern[0] != '\0') {
		int confirm = siril_confirm_dialog(_("Undebayered CFA image loaded"),
				_("You are about to apply a function that is not intended for use on an undebayered CFA image. Are you sure you wish to proceed?"), _("Proceed"));
		retval = confirm ? TRUE : FALSE;
	} else {
		retval = TRUE;
	}
	return retval;
}

point closest_point_on_line(point in, point p1, point p2) {
	point out = { 0 };
	if (p1.x == p2.x) {
		out.x = p1.x;
		out.y = in.y;
	} else if (p1.y == p2.y) {
		out.x = in.x;
		out.y = p1.y;
	} else {
		double a = in.x;
		double b = in.y;
		double x1 = p1.x;
		double y1 = p1.y;
		double x2 = p2.x;
		double y2 = p2.y;
		double m1 = (y2-y1)/(x2-x1);
		double m2 = -1.0 / m1;
		double x = (m1*x1-m2*a+b-y1)/(m1-m2);
		double y = m2*(x-a)+b;
		out.x = x;
		out.y = y;
	}
	return out;
}

// Add the GtkFileFilter filter_name to the GtkFileChooser widget_name
// If the widget already obeys the filter then it is not added a second time
void siril_set_file_filter(GtkFileChooser* chooser, const gchar* filter_name, gchar *filter_display_name) {
	GtkFileFilter* filter = GTK_FILE_FILTER(lookup_gobject(filter_name));
	GSList* filter_list = gtk_file_chooser_list_filters(chooser);
	gtk_file_filter_set_name(filter, filter_display_name);
	GSList* iterator = filter_list;
	gboolean add_filter = TRUE;
	while (iterator) {
		if (iterator->data == filter) {
			add_filter = FALSE;
			break;
		}
		iterator++;
	}
	g_slist_free(filter_list);
	if (add_filter)
		gtk_file_chooser_add_filter(chooser, filter);
}

/* apply_limits moved to core/utils.c */

gboolean value_check(fits *fit) {
	if (fit->type == DATA_USHORT)
		return TRUE;

	double maxval, minval;
	int retval = quick_minmax(fit, &minval, &maxval);
	if (retval)
		return TRUE;
	OverrangeResponse result;
	if (maxval > 10) {
		// All FITS files are scaled to [0,1] on opening, so if we have extreme values at this stage then these are highly likely hot pixels and we should clip
		result = RESPONSE_CLIP;
	// Now we know we need a scaling response
	} else if (minval < -0.1) {
		// As above, images are opened in [0,1] so significant negative values are outliers and should be cropped, so
		result = RESPONSE_RESCALE_CLIPNEG;
	} else {
		// Small negative values may be legitimate as the result of noise, and in this case we should rescale all
		result = RESPONSE_RESCALE_ALL;
	}
	result = apply_limits(fit, minval, maxval, result);
	return result == RESPONSE_CANCEL ? FALSE : TRUE;
}

GdkRGBA uint32_to_gdk_rgba(uint32_t packed_rgba) {
    GdkRGBA rgba;

    // Extract each 8-bit component (assuming RGBA order)
    rgba.red   = ((packed_rgba >> 24) & 0xFF) / 255.0;
    rgba.green = ((packed_rgba >> 16) & 0xFF) / 255.0;
    rgba.blue  = ((packed_rgba >> 8)  & 0xFF) / 255.0;
    rgba.alpha = (packed_rgba & 0xFF) / 255.0;

    return rgba;
}

/* ── HEIF image-selector dialog ─────────────────────────────────────────── */
#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#define MAX_THUMBNAIL_SIZE com.pref.gui.thumbnail_size

struct HeifImage {
	uint32_t ID;
	char caption[100];
	struct heif_image *thumbnail;
	int width, height;
};

static gboolean load_thumbnails(struct heif_context *heif, struct HeifImage *images) {
	int numImages = heif_context_get_number_of_top_level_images(heif);
	uint32_t *IDs = malloc(numImages * sizeof(uint32_t));
	heif_context_get_list_of_top_level_image_IDs(heif, IDs, numImages);

	for (int i = 0; i < numImages; i++) {
		images[i].ID = IDs[i];
		images[i].caption[0] = 0;
		images[i].thumbnail = NULL;

		struct heif_image_handle *handle;
		struct heif_error err = heif_context_get_image_handle(heif, IDs[i], &handle);
		if (err.code) { g_printf("%s\n", err.message); continue; }

		int width  = heif_image_handle_get_width(handle);
		int height = heif_image_handle_get_height(handle);
		if (heif_image_handle_is_primary_image(handle))
			sprintf(images[i].caption, "%dx%d (%s)", width, height, _("primary"));
		else
			sprintf(images[i].caption, "%dx%d", width, height);

		struct heif_image_handle *thumbnail_handle;
		heif_item_id thumbnail_ID;
		int nThumbnails = heif_image_handle_get_list_of_thumbnail_IDs(handle, &thumbnail_ID, 1);
		if (nThumbnails > 0) {
			err = heif_image_handle_get_thumbnail(handle, thumbnail_ID, &thumbnail_handle);
			if (err.code) { g_printf("%s\n", err.message); continue; }
		} else {
			err = heif_context_get_image_handle(heif, IDs[i], &thumbnail_handle);
			if (err.code) { g_printf("%s\n", err.message); continue; }
		}

		struct heif_image *thumbnail_img;
		err = heif_decode_image(thumbnail_handle, &thumbnail_img,
				heif_colorspace_RGB, heif_chroma_interleaved_RGB, NULL);
		if (err.code) { g_printf("%s\n", err.message); continue; }

		int thumbnail_width  = heif_image_handle_get_width(thumbnail_handle);
		int thumbnail_height = heif_image_handle_get_height(thumbnail_handle);
		if (thumbnail_width > MAX_THUMBNAIL_SIZE || thumbnail_height > MAX_THUMBNAIL_SIZE) {
			float factor_h = thumbnail_width  / (float)MAX_THUMBNAIL_SIZE;
			float factor_v = thumbnail_height / (float)MAX_THUMBNAIL_SIZE;
			int new_width, new_height;
			if (factor_v > factor_h) {
				new_height = MAX_THUMBNAIL_SIZE;
				new_width  = thumbnail_width / factor_v;
			} else {
				new_height = thumbnail_height / factor_h;
				new_width  = MAX_THUMBNAIL_SIZE;
			}
			struct heif_image *scaled_img = NULL;
			err = heif_image_scale_image(thumbnail_img, &scaled_img, new_width, new_height, NULL);
			if (err.code) { g_printf("%s\n", err.message); continue; }
			heif_image_release(thumbnail_img);
			thumbnail_img    = scaled_img;
			thumbnail_width  = new_width;
			thumbnail_height = new_height;
		}
		heif_image_handle_release(thumbnail_handle);
		heif_image_handle_release(handle);
		images[i].thumbnail = thumbnail_img;
		images[i].width     = thumbnail_width;
		images[i].height    = thumbnail_height;
	}
	return TRUE;
}

gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image) {
	int numImages = heif_context_get_number_of_top_level_images(heif);
	struct HeifImage *heif_images = malloc(numImages * sizeof(struct HeifImage));
	gboolean success = load_thumbnails(heif, heif_images);
	if (!success) { free(heif_images); return FALSE; }

	GtkWidget *dlg = gtk_dialog_new_with_buttons(_("Load HEIF image content"),
			GTK_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"))),
			GTK_DIALOG_MODAL,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_OK"), GTK_RESPONSE_OK, NULL);
	gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_OK);

	GtkContainer *content_area = GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(dlg)));
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 12);

	GtkWidget *frame = gtk_frame_new(_("Select image"));
	gtk_container_add(content_area, GTK_WIDGET(frame));
	gtk_widget_show(frame);

	GtkListStore *liststore;
	GtkTreeIter iter;
	liststore = gtk_list_store_new(2, G_TYPE_STRING, GDK_TYPE_PIXBUF);
	for (int i = 0; i < numImages; i++) {
		gtk_list_store_append(liststore, &iter);
		gtk_list_store_set(liststore, &iter, 0, heif_images[i].caption, -1);
		int stride;
		const uint8_t *data = heif_image_get_plane_readonly(
				heif_images[i].thumbnail, heif_channel_interleaved, &stride);
		GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB,
				FALSE, 8, heif_images[i].width, heif_images[i].height, stride,
				NULL, NULL);
		gtk_list_store_set(liststore, &iter, 1, pixbuf, -1);
	}

	GtkWidget *iconview = gtk_icon_view_new();
	gtk_icon_view_set_model((GtkIconView *)iconview, (GtkTreeModel *)liststore);
	gtk_icon_view_set_text_column((GtkIconView *)iconview, 0);
	gtk_icon_view_set_pixbuf_column((GtkIconView *)iconview, 1);
	gtk_icon_view_set_item_width((GtkIconView *)iconview, MAX_THUMBNAIL_SIZE);

	GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_widget_set_size_request(scroll, -1, 400);
	g_object_set(scroll, "expand", TRUE, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
			GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
	gtk_container_add(GTK_CONTAINER(frame), scroll);
	gtk_container_add(GTK_CONTAINER(scroll), iconview);
	gtk_widget_show(scroll);
	gtk_widget_show(iconview);

	/* pre-select the primary image */
	int selected_idx = -1;
	for (int i = 0; i < numImages; i++) {
		if (heif_images[i].ID == *selected_image) { selected_idx = i; break; }
	}
	if (selected_idx != -1) {
		GtkTreePath *path = gtk_tree_path_new_from_indices(selected_idx, -1);
		gtk_icon_view_select_path((GtkIconView *)iconview, path);
		gtk_tree_path_free(path);
	}

	gtk_widget_show(dlg);
	gboolean run = (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_OK);
	if (run) {
		GList *selected_items = gtk_icon_view_get_selected_items((GtkIconView *)iconview);
		if (selected_items) {
			GtkTreePath *path = (GtkTreePath *)(selected_items->data);
			const gint *indices = gtk_tree_path_get_indices(path);
			*selected_image = heif_images[indices[0]].ID;
			g_list_free_full(selected_items, (GDestroyNotify)gtk_tree_path_free);
		}
	}
	gtk_widget_destroy(dlg);

	for (int i = 0; i < numImages; i++)
		heif_image_release(heif_images[i].thumbnail);
	free(heif_images);
	return run;
}
#endif /* HAVE_LIBHEIF */

/* ── File-chooser filename helpers ────────────────────────────────────────── */

gchar *siril_file_chooser_get_filename(GtkFileChooser *chooser) {
	gchar *filename = NULL;
	gchar *uri = gtk_file_chooser_get_uri(GTK_FILE_CHOOSER(chooser));

	if (uri != NULL) {
		filename = g_filename_from_uri(uri, NULL, NULL);
		if (filename != NULL) {
			char *scheme = g_uri_parse_scheme(uri);
			if (g_strcmp0(scheme, "file") == 0) {
				printf("The URI points to a local file.\n");
			} else {
				printf("The URI is non-local (scheme: %s).\n", uri);
			}
			g_free(scheme);
		}
		g_free(uri);
	}
	/* Fallback for Save dialogs where the typed name has no URI yet
	 * (gtk_file_chooser_get_uri returns NULL before the file exists). */
	if (!filename)
		filename = gtk_file_chooser_get_current_name(chooser);
	return filename;
}

GSList *siril_file_chooser_get_filenames(GtkFileChooser *chooser) {
	GSList *filenames = NULL;
	GSList *uris = gtk_file_chooser_get_uris(GTK_FILE_CHOOSER(chooser));

	for (GSList *iter = uris; iter != NULL; iter = g_slist_next(iter)) {
		const gchar *uri = (const gchar *)iter->data;
		gchar *filename = g_filename_from_uri(uri, NULL, NULL);

		if (filename != NULL) {
			printf("filename=%s\n", filename);
			filenames = g_slist_append(filenames, filename);
		}
	}

	g_slist_free(uris);

	return filenames;
}
