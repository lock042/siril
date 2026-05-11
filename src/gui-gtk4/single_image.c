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
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/message_dialog.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "io/sequence.h"

static const char *drawing_area[] = { "drawingarear", "drawingareag", "drawingareab", "drawingareargb"};

/* Phase 20 / 8B: GTK4 replaced gtk_drag_dest_set / GtkSelectionData /
 * GdkDragContext with GtkDropTarget.  GdkFileList is the native payload
 * type for file/URI drops. */
static gboolean on_drawingarea_drop(GtkDropTarget *target, const GValue *value,
		double x, double y, gpointer user_data) {
	(void) target;
	(void) x;
	(void) y;
	(void) user_data;

	if (!G_VALUE_HOLDS(value, GDK_TYPE_FILE_LIST))
		return FALSE;

	GdkFileList *file_list = g_value_get_boxed(value);
	GSList *files = gdk_file_list_get_files(file_list);
	if (!files)
		return FALSE;

	GFile *file = files->data;
	gchar *filename = g_file_get_path(file);
	if (!filename)
		return FALSE;

	const char *src_ext = get_filename_ext(filename);
	gboolean confirm = TRUE;
	if (src_ext) {
		if ((!strncmp(src_ext, "seq", 4)) || (get_type_for_extension(src_ext) != TYPEUNDEF)) {
			if (single_image_is_loaded() || sequence_is_loaded()) {
				confirm = siril_confirm_dialog(_("An image (or sequence) is already loaded"),
						_("Are you sure you want to close everything and open the new image?"), _("Open"));
			}
			if (confirm) {
				if (!strncmp(src_ext, "seq", 4)) {
					gchar *sequence_dir = g_path_get_dirname(filename);
					if (!siril_change_dir(sequence_dir, NULL)) {
						if (check_seq()) {
							siril_log_message(_("No sequence `%s' found.\n"), filename);
						} else {
							set_seq(filename);
							if (com.pref.wd)
								g_free(com.pref.wd);
							com.pref.wd = g_strdup(com.wd);
							if (!com.script) {
								gui_iface.populate_seq_combo(filename);
								gui_iface.set_gui_cwd();
							}
						}
						g_free(sequence_dir);
					}
				} else {
					if (get_type_for_extension(src_ext) != TYPEUNDEF) {
						open_single_image(filename);
					}
				}
			}
		}
	}
	g_free(filename);
	return TRUE;
}

static GdkDragAction on_drop_enter(GtkDropTarget *target, double x, double y,
                                   gpointer user_data) {
	(void)target; (void)x; (void)y; (void)user_data;
	return GDK_ACTION_COPY;
}

static GdkDragAction on_drop_motion(GtkDropTarget *target, double x, double y,
                                    gpointer user_data) {
	(void)target; (void)x; (void)y; (void)user_data;
	return GDK_ACTION_COPY;
}

static void install_one_drop_target(GtkWidget *w) {
	if (!w) return;
	GtkEventController *old = g_object_get_data(G_OBJECT(w), "siril-drop-target");
	if (old) {
		gtk_widget_remove_controller(w, old);
		g_object_set_data(G_OBJECT(w), "siril-drop-target", NULL);
	}
	GtkDropTarget *target = gtk_drop_target_new(GDK_TYPE_FILE_LIST, GDK_ACTION_COPY);
	g_signal_connect(target, "enter",  G_CALLBACK(on_drop_enter),       NULL);
	g_signal_connect(target, "motion", G_CALLBACK(on_drop_motion),      NULL);
	g_signal_connect(target, "drop",   G_CALLBACK(on_drawingarea_drop), NULL);
	gtk_widget_add_controller(w, GTK_EVENT_CONTROLLER(target));
	g_object_set_data(G_OBJECT(w), "siril-drop-target", target);
}

static void install_all_drop_targets(void) {
	for (int i = 0; i < G_N_ELEMENTS(drawing_area); i++) {
		GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(gui.builder, drawing_area[i]));
		install_one_drop_target(w);
	}
}

/* Re-install drop targets and toggle can-target on the toplevel.  GTK
 * 4.14 has a surface-state bug where, after an xdg_toplevel.configure
 * (resize / maximise change), GDK_DRAG_ENTER events stop being
 * dispatched to widget DnD controllers.  Tearing down and re-creating
 * the GtkDropTarget plus a can-target toggle is the only workaround
 * that consistently re-establishes the dispatch path. */
static gboolean refresh_drop_targets_idle(gpointer ud) {
	(void)ud;
	GtkWidget *win = GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"));
	if (win) {
		gtk_widget_set_can_target(win, FALSE);
		gtk_widget_set_can_target(win, TRUE);
	}
	install_all_drop_targets();
	return G_SOURCE_REMOVE;
}

static void on_window_state_changed(GObject *obj, GParamSpec *pspec, gpointer ud) {
	(void)obj; (void)pspec; (void)ud;
	g_idle_add(refresh_drop_targets_idle, NULL);
}

void siril_drag_single_image_set_dest() {
	install_all_drop_targets();

	GtkWidget *win = GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"));
	if (win && !g_object_get_data(G_OBJECT(win), "siril-drop-state-hooks")) {
		g_signal_connect(win, "notify::maximized",
		                 G_CALLBACK(on_window_state_changed), NULL);
		g_signal_connect(win, "notify::fullscreened",
		                 G_CALLBACK(on_window_state_changed), NULL);
		g_signal_connect(win, "notify::default-width",
		                 G_CALLBACK(on_window_state_changed), NULL);
		g_signal_connect(win, "notify::default-height",
		                 G_CALLBACK(on_window_state_changed), NULL);
		g_object_set_data(G_OBJECT(win), "siril-drop-state-hooks",
		                  GINT_TO_POINTER(1));
	}
}
