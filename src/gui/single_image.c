/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "io/sequence.h"

static const char *drawing_area[] = { "drawingarear", "drawingareag", "drawingareab", "drawingareargb"};

void siril_drag_single_image_set_dest() {
	for (int i = 0; i < G_N_ELEMENTS(drawing_area); i++) {
		gtk_drag_dest_set(lookup_widget(drawing_area[i]), GTK_DEST_DEFAULT_MOTION | GTK_DEST_DEFAULT_DROP | GTK_DEST_DEFAULT_HIGHLIGHT, NULL, 0, GDK_ACTION_COPY | GDK_ACTION_ASK);
		gtk_drag_dest_add_uri_targets(lookup_widget(drawing_area[i]));
	}
}

void on_drawingarea_drag_data_received(GtkWidget *widget,
		GdkDragContext *context, gint x, gint y,
		GtkSelectionData *selection_data, guint info, guint time,
		gpointer user_data) {

	GdkAtom target;
	GtkWidget *src;

	target = gtk_selection_data_get_target(selection_data);

	if (!gtk_targets_include_uri(&target, 1))
		return;

	/* if the request is from another process this will return NULL */
	src = gtk_drag_get_source_widget(context);

	/* if the drag request originates from the current siril instance, ignore
	 the request if the source window is the same as the dest window */
	if (src && gtk_widget_get_toplevel(src) == gtk_widget_get_toplevel(widget)) {
		gdk_drag_status(context, 0, time);
		return;
	}

	if (gdk_drag_context_get_suggested_action(context) & GDK_ACTION_COPY) {
		GError *error = NULL;
		gboolean confirm = TRUE;
		const guchar *data = gtk_selection_data_get_data(selection_data);
		gchar **uris = g_uri_list_extract_uris((gchar *) data);
		/* we can open only one image */
		gchar *filename = g_filename_from_uri(uris[0], NULL, &error);
		const char *src_ext = get_filename_ext(filename);
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
									populate_seqcombo(filename);
									gui_function(set_GUI_CWD, NULL);
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
		g_strfreev(uris);
	}
}
