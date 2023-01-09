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

#include "core/siril.h"
#include "core/proto.h"
#include <glib.h>

#include "utils.h"

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

void control_window_switch_to_tab(main_tabs tab) {
	if (com.script)
		return;
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook_center_box"));
	gtk_notebook_set_current_page(notebook, tab);
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
		GdkPixbuf *new_pixbuf = gdk_pixbuf_scale_simple(pixbuf, width, height,
				GDK_INTERP_BILINEAR);
		image = gtk_image_new_from_pixbuf(new_pixbuf);
		g_object_unref(new_pixbuf);
	} else {
		image = gtk_image_new_from_icon_name("dialog-information-symbolic",
					GTK_ICON_SIZE_DIALOG);
	}

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 80);

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
		if (space < 1073741824) { // we want to warn user of space is less than 1GiB
			set_class = TRUE;
		}
		gchar *mem = g_format_size_full(space, G_FORMAT_SIZE_IEC_UNITS);
		str = g_strdup_printf(_("Disk Space: %s"), mem);
		g_free(mem);
	} else {
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

/* Managing GUI calls *
 * We try to separate core concerns from GUI (GTK+) concerns to keep the core
 * code efficient and clean, and avoir calling GTK+ functions in threads other
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
	data->idle(data->user);

	siril_debug_print("idle %p signaling end\n", data->idle);
	g_mutex_lock(&data->mutex);
	data->idle_finished = TRUE;
	g_cond_signal(&data->cond);
	g_mutex_unlock(&data->mutex);

	return FALSE;
}

void execute_idle_and_wait_for_it(gboolean (* idle)(gpointer), gpointer arg) {
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
