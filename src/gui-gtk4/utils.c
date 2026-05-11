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
#include "gui-gtk4/callbacks.h"
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
		if (args->class_to_add)
			gtk_widget_add_css_class(widget, args->class_to_add);
		if (args->class_to_remove)
			gtk_widget_remove_css_class(widget, args->class_to_remove);
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
	g_idle_add(set_label_text_idle, data);
}

void widget_set_class(GtkWidget *entry, const char *class_to_add, const char *class_to_remove) {
	if (class_to_add)
		gtk_widget_add_css_class(entry, class_to_add);
	if (class_to_remove)
		gtk_widget_remove_css_class(entry, class_to_remove);
}

/* Force a sensitivity propagation cycle on the given subtree.  In GTK4,
 * GtkExpander's user-child subtree can keep GTK_STATE_FLAG_INSENSITIVE
 * stuck after an initial set_sensitive(FALSE)/set_sensitive(TRUE) cycle.
 * The flag sits in the descendants' state_flags but gtk_widget_set_sensitive(TRUE)
 * is a no-op when sensitive is already TRUE — so no re-propagation runs.
 * Workaround: toggle the local sensitive property FALSE then TRUE on each
 * affected widget; the FALSE→TRUE transition is non-trivial and forces GTK
 * to run a fresh propagation, clearing the stuck INSENSITIVE bit. */
static void force_sensitivity_repropagate(GtkWidget *w) {
	if (!w) return;
	if (gtk_widget_get_sensitive(w)) {
		gtk_widget_set_sensitive(w, FALSE);
		gtk_widget_set_sensitive(w, TRUE);
	}
	for (GtkWidget *c = gtk_widget_get_first_child(w); c; c = gtk_widget_get_next_sibling(c))
		force_sensitivity_repropagate(c);
}

void expander_clear_stuck_insensitive(GtkWidget *expander) {
	if (!expander) return;
	for (GtkWidget *c = gtk_widget_get_first_child(expander); c; c = gtk_widget_get_next_sibling(c))
		force_sensitivity_repropagate(c);
}

void clear_stuck_insensitive_in_tree(GtkWidget *root) {
	if (!root) return;
	if (GTK_IS_EXPANDER(root))
		expander_clear_stuck_insensitive(root);
	for (GtkWidget *c = gtk_widget_get_first_child(root); c; c = gtk_widget_get_next_sibling(c))
		clear_stuck_insensitive_in_tree(c);
}

void clear_stuck_insensitive_in_builder(GtkBuilder *builder) {
	if (!builder) return;
	GSList *objects = gtk_builder_get_objects(builder);
	for (GSList *l = objects; l; l = l->next) {
		if (GTK_IS_EXPANDER(l->data))
			expander_clear_stuck_insensitive(GTK_WIDGET(l->data));
	}
	g_slist_free(objects);
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
 * @param paintable will be shown in the popover (NULL for an info icon)
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPaintable *paintable) {
	GtkWidget *popover, *box, *image, *label;

	popover = gtk_popover_new();
	/* The supplied `widget` is sometimes a GtkMenuButton (e.g. the
	 * snapshot/save buttons on the headerbar).  GtkMenuButton owns its
	 * own popover machinery; calling gtk_widget_set_parent on a second
	 * popover with the menu button as parent crashes inside realize
	 * because the menu button's surface isn't a valid popup parent.
	 *
	 * Walking up via gtk_widget_get_root also produced a non-window
	 * GtkRoot in some cases (notably popovers inside popovers), so we
	 * just look up the main control_window from the builder.  It's
	 * always realized by the time any of these notifications fire and
	 * always provides a valid native surface. */
	GtkWidget *parent = NULL;
	if (gui.builder) {
		GObject *obj = gtk_builder_get_object(gui.builder, "control_window");
		if (GTK_IS_WIDGET(obj))
			parent = GTK_WIDGET(obj);
	}
	if (!parent && widget) {
		GtkRoot *r = gtk_widget_get_root(widget);
		if (r && GTK_IS_WIDGET(r)) parent = GTK_WIDGET(r);
	}
	if (parent)
		gtk_widget_set_parent(popover, parent);
	if (widget && parent && parent != widget) {
		/* Point the popover at the originating widget so the
		 * notification visually anchors on the snapshot button. */
		graphene_rect_t bounds;
		if (gtk_widget_compute_bounds(widget, parent, &bounds)) {
			GdkRectangle rect = {
				(int) bounds.origin.x,
				(int) bounds.origin.y,
				(int) bounds.size.width,
				(int) bounds.size.height,
			};
			gtk_popover_set_pointing_to(GTK_POPOVER(popover), &rect);
		}
	}
	label = gtk_label_new(NULL);
	gtk_label_set_wrap_mode(GTK_LABEL(label), PANGO_WRAP_WORD);
	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);

	if (paintable) {
		int iw = gdk_paintable_get_intrinsic_width(paintable);
		int ih = gdk_paintable_get_intrinsic_height(paintable);
		float ratio = iw > 0 ? (float)ih / (float)iw : 1.0f;
		int width = 128, height = (int)(128 * ratio);
		/* GtkPicture handles arbitrary-size paintables: it scales the
		 * paintable into the widget's allocation with preserved aspect
		 * ratio.  No need to pre-scale the source. */
		image = gtk_picture_new_for_paintable(paintable);
		gtk_picture_set_content_fit(GTK_PICTURE(image), GTK_CONTENT_FIT_CONTAIN);
		gtk_picture_set_can_shrink(GTK_PICTURE(image), TRUE);
		gtk_widget_set_size_request(image, width, height);
	} else {
		image = gtk_image_new_from_icon_name("dialog-information-symbolic");
	}

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 100);

	gtk_widget_set_halign(image, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(image, GTK_ALIGN_CENTER);
	gtk_widget_set_margin_start(image, 10);
	gtk_widget_set_margin_end(image, 10);
	gtk_widget_set_margin_top(image, 10);
	gtk_widget_set_margin_bottom(image, 10);
	gtk_box_append(GTK_BOX(box), image);
	gtk_widget_set_halign(label, GTK_ALIGN_CENTER);
	gtk_widget_set_valign(label, GTK_ALIGN_CENTER);
	gtk_widget_set_margin_start(label, 10);
	gtk_widget_set_margin_end(label, 10);
	gtk_widget_set_margin_top(label, 10);
	gtk_widget_set_margin_bottom(label, 10);
	gtk_box_append(GTK_BOX(box), label);
	gtk_popover_set_child(GTK_POPOVER(popover), box);

	/* make all sensitive in case where parent is not */
	gtk_widget_set_sensitive(label, TRUE);
	gtk_widget_set_sensitive(box, TRUE);
	gtk_widget_set_sensitive(popover, TRUE);

	gtk_widget_set_visible(box, TRUE);

	return popover;
}

/* GTK4 replacement for gtk_container_get_children: returns a newly-allocated
 * GList of direct children in source order. Caller frees with g_list_free(). */
GList *widget_list_children(GtkWidget *parent) {
	GList *out = NULL;
	for (GtkWidget *c = gtk_widget_get_last_child(parent); c; c = gtk_widget_get_prev_sibling(c))
		out = g_list_prepend(out, c);
	return out;
}

/* Phase 11 transition: this helper still serves the GtkTreeView callers
 * (PSF_list.c, keywords_tree.c, sequence_list.c, pixelmath.c, conversion.c)
 * that have not yet been migrated to GListStore + GtkSelectionModel.
 * Once those files migrate to GtkSingleSelection / GtkMultiSelection, the
 * caller-side patterns become "iterate the selection model's GtkBitset and
 * read items by position" — at that point this function can be deleted.
 * Until then, silence the unavoidable deprecation warnings here. */
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
	/* Phase 15: GTK_STYLE_CLASS_SUGGESTED_ACTION removed from GTK4
	 * headers; use the CSS class name directly. */
	gtk_widget_add_css_class(widget, "suggested-action");
}

void unset_suggested(GtkWidget *widget) {
	gtk_widget_remove_css_class(widget, "suggested-action");
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
	g_idle_add(wrapping_idle, data);
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

/* Phase 14G.7: full GTK4 rewrite.  The legacy version used a custom
 * GtkDialog with a GtkIconView fed by a GtkListStore — all of which
 * are gone (GtkIconView removed) or deprecated (GtkDialog).  GTK4
 * replacement: a GtkWindow hosting a GtkGridView fed by a GListStore
 * of per-row GObjects.  Sync semantics preserved via nested GMainLoop.
 *
 * Per-row object: SirilHeifRow holds caption + GdkTexture (cached). */
#define SIRIL_TYPE_HEIF_ROW (siril_heif_row_get_type())
G_DECLARE_FINAL_TYPE(SirilHeifRow, siril_heif_row, SIRIL, HEIF_ROW, GObject)
struct _SirilHeifRow {
	GObject     parent_instance;
	gchar      *caption;
	GdkTexture *texture;
	guint       array_index;  /* offset into the heif_images[] array */
};
G_DEFINE_TYPE(SirilHeifRow, siril_heif_row, G_TYPE_OBJECT)
static void siril_heif_row_finalize(GObject *obj) {
	SirilHeifRow *self = SIRIL_HEIF_ROW(obj);
	g_clear_pointer(&self->caption, g_free);
	g_clear_object(&self->texture);
	G_OBJECT_CLASS(siril_heif_row_parent_class)->finalize(obj);
}
static void siril_heif_row_class_init(SirilHeifRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_heif_row_finalize;
}
static void siril_heif_row_init(SirilHeifRow *self) { (void)self; }

/* GtkSignalListItemFactory callbacks for the GtkGridView cells. */
static void heif_factory_setup(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	GtkWidget *image = gtk_image_new();
	gtk_widget_set_size_request(image, MAX_THUMBNAIL_SIZE, MAX_THUMBNAIL_SIZE);
	GtkWidget *label = gtk_label_new(NULL);
	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 18);
	gtk_box_append(GTK_BOX(box), image);
	gtk_box_append(GTK_BOX(box), label);
	gtk_list_item_set_child(li, box);
}
static void heif_factory_bind(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *box = gtk_list_item_get_child(li);
	if (!box) return;
	GtkWidget *image = gtk_widget_get_first_child(box);
	GtkWidget *label = gtk_widget_get_last_child(box);
	SirilHeifRow *row = SIRIL_HEIF_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	if (row->texture)
		gtk_image_set_from_paintable(GTK_IMAGE(image), GDK_PAINTABLE(row->texture));
	gtk_label_set_text(GTK_LABEL(label), row->caption ? row->caption : "");
}

struct heif_dlg_ctx {
	GMainLoop          *loop;
	gboolean            accepted;
	GtkSingleSelection *sel;
	GListStore         *store;
};

static void heif_dlg_cancel_clicked(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct heif_dlg_ctx *c = ud;
	c->accepted = FALSE;
	if (c->loop && g_main_loop_is_running(c->loop))
		g_main_loop_quit(c->loop);
}
static void heif_dlg_ok_clicked(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct heif_dlg_ctx *c = ud;
	c->accepted = TRUE;
	if (c->loop && g_main_loop_is_running(c->loop))
		g_main_loop_quit(c->loop);
}
static gboolean heif_dlg_close_request(GtkWindow *w, gpointer ud) {
	(void)w;
	struct heif_dlg_ctx *c = ud;
	c->accepted = FALSE;
	if (c->loop && g_main_loop_is_running(c->loop))
		g_main_loop_quit(c->loop);
	return FALSE;
}
static void heif_dlg_grid_activate(GtkGridView *gv, guint position, gpointer ud) {
	(void)gv; (void)position;
	/* double-click on a thumbnail accepts the dialog */
	heif_dlg_ok_clicked(NULL, ud);
}

gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image) {
	int numImages = heif_context_get_number_of_top_level_images(heif);
	struct HeifImage *heif_images = malloc(numImages * sizeof(struct HeifImage));
	gboolean success = load_thumbnails(heif, heif_images);
	if (!success) { free(heif_images); return FALSE; }

	GtkWidget *dlg = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dlg), _("Load HEIF image content"));
	gtk_window_set_modal(GTK_WINDOW(dlg), TRUE);
	gtk_window_set_destroy_with_parent(GTK_WINDOW(dlg), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dlg),
		GTK_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"))));
	gtk_window_set_default_size(GTK_WINDOW(dlg), 600, 500);

	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_window_set_child(GTK_WINDOW(dlg), content_area);

	GtkWidget *frame = gtk_frame_new(_("Select image"));
	gtk_widget_set_hexpand(frame, TRUE);
	gtk_widget_set_vexpand(frame, TRUE);
	gtk_box_append(GTK_BOX(content_area), frame);

	GListStore *store = g_list_store_new(SIRIL_TYPE_HEIF_ROW);
	int selected_idx = -1;
	for (int i = 0; i < numImages; i++) {
		int stride;
		const uint8_t *data = heif_image_get_plane_readonly(
				heif_images[i].thumbnail, heif_channel_interleaved, &stride);
		SirilHeifRow *row = g_object_new(SIRIL_TYPE_HEIF_ROW, NULL);
		row->caption = g_strdup(heif_images[i].caption);
		/* The HEIF data is borrowed from heif_images[].thumbnail and
		 * stays valid until heif_image_release() runs below at function
		 * end.  The store + rows (which own the textures) are dropped
		 * before then, so a static GBytes wrapper is safe. */
		row->texture = data
			? siril_texture_from_rgb_bytes(data,
					(gsize)stride * heif_images[i].height,
					heif_images[i].width, heif_images[i].height,
					stride, FALSE, NULL, NULL)
			: NULL;
		row->array_index = (guint)i;
		g_list_store_append(store, row);
		g_object_unref(row);
		if (heif_images[i].ID == *selected_image)
			selected_idx = i;
	}

	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(store)));
	GtkSignalListItemFactory *factory = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(factory, "setup", G_CALLBACK(heif_factory_setup), NULL);
	g_signal_connect(factory, "bind",  G_CALLBACK(heif_factory_bind),  NULL);

	GtkWidget *gridview = gtk_grid_view_new(GTK_SELECTION_MODEL(sel),
		GTK_LIST_ITEM_FACTORY(factory));
	gtk_grid_view_set_max_columns(GTK_GRID_VIEW(gridview), 6);
	gtk_grid_view_set_min_columns(GTK_GRID_VIEW(gridview), 1);

	GtkWidget *scroll = gtk_scrolled_window_new();
	gtk_widget_set_size_request(scroll, -1, 400);
	gtk_widget_set_hexpand(scroll, TRUE);
	gtk_widget_set_vexpand(scroll, TRUE);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
			GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
	gtk_frame_set_child(GTK_FRAME(frame), scroll);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroll), gridview);

	if (selected_idx >= 0)
		gtk_single_selection_set_selected(sel, (guint)selected_idx);

	/* Action area */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *btn_cancel = gtk_button_new_with_mnemonic(_("_Cancel"));
	GtkWidget *btn_ok     = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(btn_ok, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), btn_cancel);
	gtk_box_append(GTK_BOX(bbox), btn_ok);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dlg), btn_ok);

	struct heif_dlg_ctx ctx = {
		.loop = g_main_loop_new(NULL, FALSE),
		.accepted = FALSE,
		.sel = sel,
		.store = store,
	};
	g_signal_connect(btn_cancel, "clicked",        G_CALLBACK(heif_dlg_cancel_clicked), &ctx);
	g_signal_connect(btn_ok,     "clicked",        G_CALLBACK(heif_dlg_ok_clicked),     &ctx);
	g_signal_connect(dlg,        "close-request",  G_CALLBACK(heif_dlg_close_request),  &ctx);
	g_signal_connect(gridview,   "activate",       G_CALLBACK(heif_dlg_grid_activate),  &ctx);

	gtk_window_present(GTK_WINDOW(dlg));
	g_main_loop_run(ctx.loop);
	g_main_loop_unref(ctx.loop);

	gboolean run = ctx.accepted;
	if (run) {
		guint pos = gtk_single_selection_get_selected(sel);
		if (pos != GTK_INVALID_LIST_POSITION) {
			SirilHeifRow *row = SIRIL_HEIF_ROW(g_list_model_get_item(G_LIST_MODEL(store), pos));
			if (row) {
				*selected_image = heif_images[row->array_index].ID;
				g_object_unref(row);
			}
		}
	}

	g_object_unref(factory);
	g_object_unref(sel);
	g_object_unref(store);
	gtk_window_destroy(GTK_WINDOW(dlg));

	for (int i = 0; i < numImages; i++)
		heif_image_release(heif_images[i].thumbnail);
	free(heif_images);
	return run;
}
#endif /* HAVE_LIBHEIF */

/* ── Path-button accessors ─────────────────────────────────────────────────
 *
 * Phase 18 replaced every GtkFileChooserButton in the .ui files with a
 * plain GtkButton wired through siril_path_button_init().  The selected
 * file path is stashed on the button via g_object_set_data so the rest
 * of the codebase can keep using "char *path" semantics.  These four
 * helpers read/write that stashed data.
 *
 * For actual dialog-based open/save flows, callers should use
 * SirilFileChooser (siril_fc_*) which wraps GtkFileDialog
 * synchronously; the deprecated GtkFileChooser interface is no longer
 * touched here. */

gchar *siril_file_chooser_get_filename(GtkWidget *button) {
	if (!button || !GTK_IS_BUTTON(button)) return NULL;
	const gchar *p = g_object_get_data(G_OBJECT(button), "siril-path");
	return p ? g_strdup(p) : NULL;
}

GSList *siril_file_chooser_get_filenames(GtkWidget *button) {
	/* Path buttons hold a single path; multi-select isn't supported
	 * on this widget.  Return a single-element list to keep the
	 * caller's iteration shape. */
	gchar *p = siril_file_chooser_get_filename(button);
	return p ? g_slist_append(NULL, p) : NULL;
}

gboolean siril_file_chooser_set_filename(GtkWidget *button, const char *path) {
	if (!button || !path || !GTK_IS_BUTTON(button)) return FALSE;
	g_object_set_data_full(G_OBJECT(button), "siril-path",
	                       g_strdup(path), g_free);
	gchar *base = g_path_get_basename(path);
	gtk_button_set_label(GTK_BUTTON(button), base && *base ? base : "(None)");
	g_free(base);
	return TRUE;
}

gboolean siril_file_chooser_set_current_folder_path(GtkWidget *button, const char *path) {
	if (!button || !path || !GTK_IS_BUTTON(button)) return FALSE;
	g_object_set_data_full(G_OBJECT(button), "siril-initial-folder",
	                       g_strdup(path), g_free);
	return TRUE;
}

/* ---- GtkFileChooserButton replacement ----
 * GTK4 removed GtkFileChooserButton.  Phase 18 turned those .ui widgets
 * into plain GtkButtons; siril_path_button_init() wires a click handler
 * that opens a GtkFileDialog and updates the button's stored path. */

typedef struct {
	gchar    *title;
	gchar    *filter_id;
	gboolean  is_folder;
} SirilPathButtonSpec;

static void siril_path_button_spec_free(gpointer p) {
	SirilPathButtonSpec *s = p;
	if (!s) return;
	g_free(s->title);
	g_free(s->filter_id);
	g_free(s);
}

static void on_siril_path_button_dialog_done(GObject *src, GAsyncResult *res, gpointer ud) {
	GtkButton *btn = GTK_BUTTON(ud);
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	const SirilPathButtonSpec *spec = g_object_get_data(G_OBJECT(btn), "siril-spec");
	GError *err = NULL;
	GFile *file = NULL;
	if (spec && spec->is_folder)
		file = gtk_file_dialog_select_folder_finish(fd, res, &err);
	else
		file = gtk_file_dialog_open_finish(fd, res, &err);
	if (file) {
		gchar *path = g_file_get_path(file);
		if (path) {
			g_object_set_data_full(G_OBJECT(btn), "siril-path", path, g_free);
			gchar *base = g_path_get_basename(path);
			gtk_button_set_label(btn, base && *base ? base : "(None)");
			g_free(base);
		}
		g_object_unref(file);
	}
	if (err)
		g_clear_error(&err);
}

static void on_siril_path_button_clicked(GtkButton *btn, gpointer user_data) {
	const SirilPathButtonSpec *spec = user_data;
	GtkFileDialog *fd = gtk_file_dialog_new();
	if (spec && spec->title)
		gtk_file_dialog_set_title(fd, spec->title);
	if (spec && spec->filter_id) {
		GtkFileFilter *filter = GTK_FILE_FILTER(
			gtk_builder_get_object(gui.builder, spec->filter_id));
		if (filter) {
			GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
			g_list_store_append(filters, filter);
			gtk_file_dialog_set_default_filter(fd, filter);
			gtk_file_dialog_set_filters(fd, G_LIST_MODEL(filters));
			g_object_unref(filters);
		}
	}
	/* Pre-populate from any previously selected path or initial-folder hint. */
	const gchar *prev = g_object_get_data(G_OBJECT(btn), "siril-path");
	if (prev && *prev) {
		GFile *gf = g_file_new_for_path(prev);
		gtk_file_dialog_set_initial_file(fd, gf);
		g_object_unref(gf);
	} else {
		const gchar *folder = g_object_get_data(G_OBJECT(btn), "siril-initial-folder");
		if (folder && *folder) {
			GFile *gf = g_file_new_for_path(folder);
			gtk_file_dialog_set_initial_folder(fd, gf);
			g_object_unref(gf);
		}
	}
	GtkRoot *root = gtk_widget_get_root(GTK_WIDGET(btn));
	GtkWindow *parent = GTK_IS_WINDOW(root) ? GTK_WINDOW(root) : NULL;
	if (spec && spec->is_folder)
		gtk_file_dialog_select_folder(fd, parent, NULL,
		                              on_siril_path_button_dialog_done, btn);
	else
		gtk_file_dialog_open(fd, parent, NULL,
		                     on_siril_path_button_dialog_done, btn);
	g_object_unref(fd);
}

/* GTK 4.10+ deprecated the per-widget gtk_style_context_add_provider in
 * favour of gtk_style_context_add_provider_for_display.  This helper
 * registers the given CSS source string once per provider_id at the
 * default display level — call sites then use gtk_widget_add_css_class
 * (and matching `.<class>` selectors in the CSS) to scope the rules. */
void siril_register_css_for_display(const char *provider_id, const char *css) {
	static GHashTable *registered = NULL;
	if (!registered)
		registered = g_hash_table_new(g_str_hash, g_str_equal);
	if (g_hash_table_contains(registered, provider_id))
		return;
	GtkCssProvider *provider = gtk_css_provider_new();
	gtk_css_provider_load_from_string(provider, css);
	GdkDisplay *display = gdk_display_get_default();
	if (display)
		gtk_style_context_add_provider_for_display(display,
			GTK_STYLE_PROVIDER(provider), GTK_STYLE_PROVIDER_PRIORITY_USER);
	g_object_unref(provider);
	g_hash_table_insert(registered, (gpointer) provider_id, GINT_TO_POINTER(1));
}

/* Load an SVG/PNG resource as a GdkPaintable via GTK4's native pipeline,
 * which handles SVG without depending on the librsvg GdkPixbuf loader.
 * Tries gtk_icon_paintable_new_for_file with a resource:// URI first
 * (this is what gtk_image_set_from_resource uses internally for vector
 * formats); falls back to a GdkTexture for raster resources. */
GdkPaintable *siril_paintable_from_resource(const char *resource_path, int size) {
	if (!resource_path) return NULL;
	gchar *uri = g_strconcat("resource://", resource_path, NULL);
	GFile *file = g_file_new_for_uri(uri);
	g_free(uri);
	if (!file) return NULL;

	/* GtkIconPaintable handles SVG via GTK's builtin renderer. */
	GtkIconPaintable *icon = gtk_icon_paintable_new_for_file(file, size, 1);
	if (icon) {
		g_object_unref(file);
		return GDK_PAINTABLE(icon);
	}

	/* Fall back to GdkTexture for raster formats. */
	GError *err = NULL;
	GdkTexture *tex = gdk_texture_new_from_file(file, &err);
	g_object_unref(file);
	if (err) {
		g_clear_error(&err);
		return NULL;
	}
	return tex ? GDK_PAINTABLE(tex) : NULL;
}

gboolean siril_cairo_paint_resource(cairo_t *cr, const char *resource_path,
                                    double x, double y, int width, int height) {
	if (!cr || !resource_path) return FALSE;
	GdkPaintable *paintable = siril_paintable_from_resource(resource_path, width);
	if (!paintable) return FALSE;
	GtkSnapshot *snap = gtk_snapshot_new();
	gdk_paintable_snapshot(paintable, GDK_SNAPSHOT(snap), width, height);
	GskRenderNode *node = gtk_snapshot_free_to_node(snap);
	g_object_unref(paintable);
	if (!node) return FALSE;
	cairo_save(cr);
	cairo_translate(cr, x, y);
	gsk_render_node_draw(node, cr);
	cairo_restore(cr);
	gsk_render_node_unref(node);
	return TRUE;
}

GdkTexture *siril_texture_from_cairo_surface(cairo_surface_t *surface) {
	if (!surface) return NULL;
	if (cairo_image_surface_get_format(surface) != CAIRO_FORMAT_ARGB32)
		return NULL;
	int w = cairo_image_surface_get_width(surface);
	int h = cairo_image_surface_get_height(surface);
	int stride = cairo_image_surface_get_stride(surface);
	cairo_surface_flush(surface);
	unsigned char *data = cairo_image_surface_get_data(surface);
	if (!data || w <= 0 || h <= 0) return NULL;
	cairo_surface_reference(surface);
	GBytes *bytes = g_bytes_new_with_free_func(
		data, (gsize)stride * h,
		(GDestroyNotify)cairo_surface_destroy, surface);
	GdkTexture *tex = gdk_memory_texture_new(
		w, h, GDK_MEMORY_B8G8R8A8_PREMULTIPLIED, bytes, stride);
	g_bytes_unref(bytes);
	return tex;
}

GdkTexture *siril_texture_from_pixbuf(GdkPixbuf *pb) {
	if (!pb || gdk_pixbuf_get_bits_per_sample(pb) != 8
	    || gdk_pixbuf_get_colorspace(pb) != GDK_COLORSPACE_RGB)
		return NULL;
	int w = gdk_pixbuf_get_width(pb);
	int h = gdk_pixbuf_get_height(pb);
	int rs = gdk_pixbuf_get_rowstride(pb);
	gboolean alpha = gdk_pixbuf_get_has_alpha(pb);
	int nch = gdk_pixbuf_get_n_channels(pb);
	if (w <= 0 || h <= 0 || rs <= 0 || (nch != 3 && nch != 4))
		return NULL;
	guchar *pixels = gdk_pixbuf_get_pixels(pb);
	g_object_ref(pb);
	GBytes *bytes = g_bytes_new_with_free_func(pixels, (gsize)rs * h,
	                                            g_object_unref, pb);
	GdkTexture *tex = gdk_memory_texture_new(
		w, h, alpha ? GDK_MEMORY_R8G8B8A8 : GDK_MEMORY_R8G8B8,
		bytes, rs);
	g_bytes_unref(bytes);
	return tex;
}

GdkTexture *siril_texture_from_rgb_bytes(const guchar *data, gsize len,
                                          int width, int height, int stride,
                                          gboolean has_alpha,
                                          GDestroyNotify notify,
                                          gpointer notify_data) {
	if (!data || width <= 0 || height <= 0 || stride <= 0)
		return NULL;
	GBytes *bytes = notify
		? g_bytes_new_with_free_func((gconstpointer)data, len, notify, notify_data)
		: g_bytes_new_static((gconstpointer)data, len);
	GdkTexture *tex = gdk_memory_texture_new(
		width, height,
		has_alpha ? GDK_MEMORY_R8G8B8A8 : GDK_MEMORY_R8G8B8,
		bytes, stride);
	g_bytes_unref(bytes);
	return tex;
}

void siril_path_button_init(GtkWidget *button, const gchar *title,
                            const gchar *filter_id, gboolean is_folder) {
	if (!button || !GTK_IS_BUTTON(button)) return;
	/* If already wired, skip. */
	if (g_object_get_data(G_OBJECT(button), "siril-spec")) return;
	SirilPathButtonSpec *spec = g_new0(SirilPathButtonSpec, 1);
	spec->title    = title ? g_strdup(title) : NULL;
	spec->filter_id = filter_id ? g_strdup(filter_id) : NULL;
	spec->is_folder = is_folder;
	g_object_set_data_full(G_OBJECT(button), "siril-spec", spec,
	                       siril_path_button_spec_free);
	g_signal_connect(button, "clicked",
	                 G_CALLBACK(on_siril_path_button_clicked), spec);
}

gchar *siril_drop_down_get_active_text(GtkDropDown *dd) {
	if (!dd) return NULL;
	GObject *item = gtk_drop_down_get_selected_item(dd);
	if (!item || !GTK_IS_STRING_OBJECT(item)) return NULL;
	const char *s = gtk_string_object_get_string(GTK_STRING_OBJECT(item));
	return s ? g_strdup(s) : NULL;
}

void siril_drop_down_clear_strings(GtkDropDown *dd) {
	if (!dd) return;
	/* If the dropdown already has a GtkStringList, empty it in-place
	 * with gtk_string_list_splice rather than swapping in a new model.
	 *
	 * Why: this routine is sometimes called from inside a notify::selected
	 * handler that fires while GtkDropDown's popup is still processing
	 * the click (its internal GtkListView::activate signal is mid-emit).
	 * Replacing the model via gtk_drop_down_set_model would drop the only
	 * strong ref to the old model and free its items while the popup is
	 * still walking them, causing a use-after-free in
	 * g_object_notify_by_pspec when the GTK4 stack unwinds.
	 *
	 * Splicing out all items keeps the same model object alive, so any
	 * still-running internal GtkListView dispatch sees an empty list
	 * instead of a freed one. */
	GListModel *m = gtk_drop_down_get_model(dd);
	if (m && GTK_IS_STRING_LIST(m)) {
		guint n = g_list_model_get_n_items(m);
		if (n > 0)
			gtk_string_list_splice(GTK_STRING_LIST(m), 0, n, NULL);
		return;
	}
	GtkStringList *sl = gtk_string_list_new(NULL);
	gtk_drop_down_set_model(dd, G_LIST_MODEL(sl));
	g_object_unref(sl);
}

void siril_drop_down_append_text(GtkDropDown *dd, const gchar *text) {
	if (!dd || !text) return;
	GListModel *m = gtk_drop_down_get_model(dd);
	GtkStringList *sl = (m && GTK_IS_STRING_LIST(m)) ? GTK_STRING_LIST(m) : NULL;
	if (!sl) {
		sl = gtk_string_list_new(NULL);
		gtk_drop_down_set_model(dd, G_LIST_MODEL(sl));
		g_object_unref(sl);
		m = gtk_drop_down_get_model(dd);
		sl = GTK_STRING_LIST(m);
	}
	gtk_string_list_append(sl, text);
}

gboolean siril_toggle_get_active(GtkWidget *w) {
	if (!w) return FALSE;
	if (GTK_IS_CHECK_BUTTON(w)) return gtk_check_button_get_active(GTK_CHECK_BUTTON(w));
	if (GTK_IS_TOGGLE_BUTTON(w)) return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
	return FALSE;
}

void siril_toggle_set_active(GtkWidget *w, gboolean active) {
	if (!w) return;
	if (GTK_IS_CHECK_BUTTON(w)) gtk_check_button_set_active(GTK_CHECK_BUTTON(w), active);
	else if (GTK_IS_TOGGLE_BUTTON(w)) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), active);
}
