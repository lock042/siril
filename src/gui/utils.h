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
#ifndef SRC_GUI_UTILS_H_
#define SRC_GUI_UTILS_H_
#include "core/siril.h"

typedef enum {
	FILE_CONVERSION,
	IMAGE_SEQ,
	PRE_PROC,
	REGISTRATION,
	PLOT,
	STACKING,
	OUTPUT_LOGS
} main_tabs;


GtkWidget* lookup_widget (const gchar *widget_name);
GObject* lookup_gobject(const gchar *gobject_name);
GtkAdjustment* lookup_adjustment(const gchar *adjustment_name);
void control_window_switch_to_tab(main_tabs tab);
GtkWidget* popover_new(GtkWidget *widget, const gchar *text);
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPixbuf *pixbuf);
GList *get_row_references_of_selected_rows(GtkTreeSelection *selection, GtkTreeModel *model);
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
const char* get_cfa_from_pattern(sensor_pattern pattern);
void interpolate_nongreen(fits *fit);
OverrangeResponse apply_limits(fits *fit, double minval, double maxval, OverrangeResponse method);
gboolean value_check(fits *fit); // checks for pixel values outside [0.0, 1.0]
gchar* get_control_window_id();

#endif /* SRC_GUI_UTILS_H_ */
