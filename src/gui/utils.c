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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include <glib.h>

#include "algos/demosaicing.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "utils.h"
#include "core/siril_log.h"
#include "message_dialog.h"

struct _label_data {
	const char *label_name;
	char *text;
	const char *class_to_add;
	const char *class_to_remove;
};

static gboolean set_label_text_idle(gpointer p) {
	struct _label_data *args = (struct _label_data *) p;
	GtkWidget *widget = lookup_widget(args->label_name);
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

static gboolean switch_tab(gpointer user_data) {
	main_tabs tab = (main_tabs) GPOINTER_TO_INT(user_data);
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook_center_box"));
	gtk_notebook_set_current_page(notebook, tab);
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
	// Return immediately if headless, idles should only be used for GUI operations
	// that must be run in the GTK context.
	if (com.headless) {
		siril_debug_print("execute_idle_and_wait_for_it called headless, this should not happen!\n");
		return;
	}

	// Check if we're already in the main thread (because if we are,
	// the mutex control below will fail and will cause a hang)
	if (g_main_context_is_owner(g_main_context_default())) {
		// If in main thread, just call the function directly
		idle(arg);
		return;
	}

	struct idle_data data = { .idle_finished = FALSE, .idle = idle, .user = arg };
	g_mutex_init(&data.mutex);
	g_cond_init(&data.cond);
	siril_debug_print("queueing idle %p\n", idle);
	gdk_threads_add_idle(wrapping_idle, &data);

	siril_debug_print("waiting for idle %p\n", idle);
	g_mutex_lock (&data.mutex);
	while (!data.idle_finished)
		g_cond_wait(&data.cond, &data.mutex);
	g_mutex_unlock (&data.mutex);

	g_mutex_clear(&data.mutex);
	g_cond_clear(&data.cond);
	siril_debug_print("idle %p wait is over\n", idle);
}

int select_vport(int vport) {
	return vport == RGB_VPORT ? GREEN_VPORT : vport;
}

gboolean check_ok_if_cfa() {
	gboolean retval;
	if (gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0') {
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

static GtkWidget* create_overrange_dialog(GtkWindow *parent, const gchar *title, const gchar *message) {
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *hbox;
	GtkWidget *image;
	GtkWidget *label;

	dialog = gtk_dialog_new_with_buttons(
		title,
		parent,
		GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
		NULL, NULL  // We'll add buttons manually
	);

	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

	// Create a horizontal box to hold the icon and message
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), 12);

	// Add warning icon
	image = gtk_image_new_from_icon_name("dialog-warning", GTK_ICON_SIZE_DIALOG);
	gtk_box_pack_start(GTK_BOX(hbox), image, FALSE, FALSE, 0);

	// Create label with the message
	label = gtk_label_new(message);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 50);  // Constrain maximum width
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);

	// Add the hbox to the content area
	gtk_box_pack_start(GTK_BOX(content_area), hbox, TRUE, TRUE, 0);

	// Create and add buttons vertically
	const char* button_labels[] = {
		_("Cancel"), _("Clip"), _("Rescale\n(values > 0 only)"), _("Rescale\n(all values)")
	};
	const char* button_tooltips[] = {
		_("Cancel without making any changes"),
		_("Clip pixel values to the range 0.0 - 1.0"),
		_("Rescale pixel values to fit the range 0.0 - 1.0, clipping negative pixel values"),
		_("Rescale pixel values to fit the range 0.0 to 1.0, applying an offset to bring negatve pixel values into range")
	};
	OverrangeResponse responses[] = {
		RESPONSE_CANCEL, RESPONSE_CLIP, RESPONSE_RESCALE_CLIPNEG, RESPONSE_RESCALE_ALL
	};

	for (int i = 0; i < G_N_ELEMENTS(responses); i++) {
		GtkWidget *button = gtk_dialog_add_button(GTK_DIALOG(dialog), button_labels[i], responses[i]);
		// Get the label widget from the button
		GtkWidget *label = gtk_bin_get_child(GTK_BIN(button));

		// Check if the child is a label
		if (GTK_IS_LABEL(label)) {
			// Set the alignment to center
			gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
			gtk_label_set_xalign(GTK_LABEL(label), 0.5);
			gtk_widget_set_valign(label, GTK_ALIGN_CENTER);
		}
		gtk_widget_set_tooltip_text(button, button_tooltips[i]);
	}

	// Set dialog to be not resizable
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);

	gtk_widget_show_all(dialog);

	return dialog;
}

// Function to apply limits based on the chosen method
OverrangeResponse apply_limits(fits *fit, double minval, double maxval, OverrangeResponse method) {
	switch (method) {
		case RESPONSE_CLIP:;
			clip(fit);
			break;
		case RESPONSE_RESCALE_CLIPNEG:;
			clipneg(fit);
			if (maxval > 1.0)
				soper(fit, (1.0 / maxval), OPER_MUL, TRUE);
			break;
		case RESPONSE_RESCALE_ALL:;
			double range = maxval - minval;
			if (minval < 0.0)
				soper(fit, minval, OPER_SUB, TRUE);
			if (range > 1.0)
				soper(fit, (1.0 / range), OPER_MUL, TRUE);
			break;
		default:
			// Covers RESPONSE_CANCEL. We have removed the Proceed button
			// Indeed, we should not use algorithm that are not intended to be used with negative pixels
			// (default is trigged when the dialog is closed with the cross, we need to initialize method to RESPONSE_CANCEL)
			// Do nothing, no need to notify gfit modified so just return
			method = RESPONSE_CANCEL;
	}

	if (fit == &gfit)
		notify_gfit_modified();

	return method;
}

gboolean value_check(fits *fit) {
	if (fit->type == DATA_USHORT)
		return TRUE;

	double maxval, minval;
	int retval = quick_minmax(fit, &minval, &maxval);
	if (retval)
		return TRUE;

	if (maxval > 1.0 || minval < 0.0) {
		gchar *msg = g_strdup_printf(_("This image contains pixel values outside the range 0.0 - 1.0 (min = %.3f, max = %.3f). This can cause unwanted behaviour. Choose how to handle this.\n"), minval, maxval);
		GtkWidget *dialog = create_overrange_dialog(siril_get_active_window(), _("Warning"), msg);
		OverrangeResponse result = (OverrangeResponse) gtk_dialog_run(GTK_DIALOG(dialog));
		gtk_widget_destroy(dialog);
		g_free(msg);

		if (result == RESPONSE_CANCEL)
			return FALSE;

		result = apply_limits(fit, minval, maxval, result);
		if (result == RESPONSE_CANCEL)
			return FALSE;
	}
	return TRUE;
}
