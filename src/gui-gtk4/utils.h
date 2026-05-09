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
#ifndef SRC_GUI_UTILS_H_
#define SRC_GUI_UTILS_H_
#include <gtk/gtk.h>
#include "core/siril.h"
#include "gui-gtk4/gtk3_event_compat.h"
#include "gui-gtk4/gui_state.h"

/* main_tabs enum moved to core/siril.h */

void gui_mutex_lock();
void gui_mutex_unlock();
GtkWidget* lookup_widget (const gchar *widget_name);
GObject* lookup_gobject(const gchar *gobject_name);
GtkAdjustment* lookup_adjustment(const gchar *adjustment_name);
void control_window_switch_to_tab(main_tabs tab);
GtkWidget* popover_new(GtkWidget *widget, const gchar *text);
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPixbuf *pixbuf);
GList *get_row_references_of_selected_rows(GtkTreeSelection *selection, GtkTreeModel *model);
GList *widget_list_children(GtkWidget *parent);
void set_GUI_MEM(guint64 used, const gchar *label);
void set_GUI_DiskSpace(gint64 mem, const gchar *label);
void set_suggested(GtkWidget *widget);
void unset_suggested(GtkWidget *widget);
void set_switcher_buttons_colors(GList *list, int n);

void widget_set_class(GtkWidget *entry, const char *class_to_add, const char *class_to_remove);

void execute_idle_and_wait_for_it(gboolean (* idle)(gpointer), gpointer arg);
int select_vport(int vport);
gboolean check_ok_if_cfa();
point closest_point_on_line(point in, point p1, point p2);
void siril_set_file_filter(GtkFileChooser* chooser, const gchar* filter_name, gchar *filter_display_name);
gboolean value_check(fits *fit); // checks for pixel values outside [0.0, 1.0]
gchar* get_control_window_id();
GdkRGBA uint32_to_gdk_rgba(uint32_t packed_rgba);

gchar    *siril_file_chooser_get_filename (GtkFileChooser *chooser);
GSList   *siril_file_chooser_get_filenames(GtkFileChooser *chooser);
gboolean  siril_file_chooser_set_filename (GtkFileChooser *chooser, const char *path);
gboolean  siril_file_chooser_set_current_folder_path(GtkFileChooser *chooser, const char *path);

/* GTK4 removed GtkFileChooserButton.  Phase 18 converted those widgets to
 * plain GtkButtons with label "(None)"; this helper wires up a click
 * handler that opens a GtkFileDialog (or select-folder) and stores the
 * chosen path on the button as g_object_data "siril-path", which the
 * siril_file_chooser_get_filename / set_filename shims understand.
 * Pass filter_id=NULL for no filter; is_folder=TRUE selects a directory. */
void siril_path_button_init(GtkWidget *button, const gchar *title,
                            const gchar *filter_id, gboolean is_folder);

/* Register a one-shot CSS provider for the default GdkDisplay containing
 * the given CSS source string.  Idempotent — a provider is only created
 * the first time a given (provider_id) key is seen; subsequent calls are
 * no-ops.  Use this together with gtk_widget_add_css_class() to scope
 * style rules in GTK 4.10+ where per-widget gtk_style_context_add_provider
 * is deprecated. */
void siril_register_css_for_display(const char *provider_id, const char *css);

/* Load an SVG (or any GdkPaintable-supported format) from a GResource via
 * GTK4's paintable infrastructure.  This works even when the librsvg
 * GdkPixbuf loader isn't installed, which is why the gdk_pixbuf_new_*
 * calls fail at runtime with "Unrecognized image file format".  Returns
 * NULL on failure.  Caller must g_object_unref the returned paintable. */
GdkPaintable *siril_paintable_from_resource(const char *resource_path, int size);

/* Convenience: render the resource at (x, y) on the cairo context
 * scaled to (width, height), using siril_paintable_from_resource()
 * internally.  Returns FALSE if the resource can't be loaded. */
gboolean siril_cairo_paint_resource(cairo_t *cr, const char *resource_path,
                                    double x, double y, int width, int height);

/* GtkDropDown helpers (replace deprecated GtkComboBoxText paths).
 * siril_drop_down_get_active_text returns a freshly-allocated copy of the
 * currently-selected string-list entry (caller g_free), or NULL if no
 * selection. */
gchar *siril_drop_down_get_active_text(GtkDropDown *dd);
/* Empty the dropdown's underlying GtkStringList (creates a new empty one if
 * the model isn't a GtkStringList). */
void   siril_drop_down_clear_strings  (GtkDropDown *dd);
/* Append a string to the dropdown's GtkStringList (creates one if needed). */
void   siril_drop_down_append_text    (GtkDropDown *dd, const gchar *text);

/* GTK4 split GtkCheckButton out from GtkToggleButton — they are
 * unrelated types now.  Many siril call sites carry a generic
 * GtkWidget* and don't know which they have.  These shims dispatch on
 * actual type.  Pass the underlying GtkWidget* (or any subclass);
 * passing NULL or a non-button widget returns FALSE/no-op. */
gboolean siril_toggle_get_active(GtkWidget *w);
void     siril_toggle_set_active(GtkWidget *w, gboolean active);

#ifdef HAVE_LIBHEIF
struct heif_context; /* forward declaration — avoids pulling in <libheif/heif.h> */
/* HEIF image-selector dialog: lets the user pick one image from a
 * multi-image HEIF file.  Returns TRUE if a selection was confirmed. */
gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image);
#endif

#endif /* SRC_GUI_UTILS_H_ */
