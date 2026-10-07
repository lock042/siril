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

#if HAVE_X11_DND_FALLBACK
#include <gdk/x11/gdkx.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#endif

static const char *drawing_area[] = { "drawingarear", "drawingareag", "drawingareab", "drawingareargb"};

/* Apply a dropped file (FITS/seq/whatever) — extension-dispatched to
 * open_single_image() or set_seq() with the usual confirm-on-replace
 * prompt.  Shared by the Wayland (GtkDropTarget) path and the X11
 * fallback. */
static void process_dropped_filename(gchar *filename) {
	const char *src_ext = get_filename_ext(filename);
	if (!src_ext)
		return;
	if (strncmp(src_ext, "seq", 4) != 0 && get_type_for_extension(src_ext) == TYPEUNDEF)
		return;

	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (!siril_confirm_dialog(_("An image (or sequence) is already loaded"),
				_("Are you sure you want to close everything and open the new image?"), _("Open")))
			return;
	}

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
		}
		g_free(sequence_dir);
	} else {
		open_single_image(filename);
	}
}

/* Wayland / portable path: GtkDropTarget delivers a GdkFileList GValue. */
static gboolean on_drawingarea_drop(GtkDropTarget *target, const GValue *value,
		double x, double y, gpointer user_data) {
	(void) target; (void) x; (void) y; (void) user_data;

	if (!G_VALUE_HOLDS(value, GDK_TYPE_FILE_LIST))
		return FALSE;
	GSList *files = gdk_file_list_get_files(g_value_get_boxed(value));
	if (!files)
		return FALSE;
	gchar *filename = g_file_get_path(files->data);
	if (!filename)
		return FALSE;
	process_dropped_filename(filename);
	g_free(filename);
	return TRUE;
}

#if HAVE_X11_DND_FALLBACK
/* X11 fallback for the upstream GTK4 bug whereby gdk_drop_read_value_async
 * never completes on the X11 backend, so GtkDropTarget's "drop" signal is
 * never emitted (see dnd-investigation.md).  We let GtkDropTarget continue
 * to answer the XDND protocol (so the source sees a valid target), but
 * intercept GDK_DROP_START in capture phase and fetch the text/uri-list
 * selection directly via Xlib, bypassing the broken value-read entirely.
 *
 * Synchronous: blocks the GUI thread for at most ~500 ms while waiting
 * for SelectionNotify.  Real responses arrive in microseconds. */

#define X11_SELECTION_TIMEOUT_US (500 * 1000)

static gchar **x11_fetch_uri_list(GtkWidget *widget, guint32 time) {
	GtkNative *native = gtk_widget_get_native(widget);
	if (!native) return NULL;
	GdkSurface *surface = gtk_native_get_surface(native);
	if (!surface || !GDK_IS_X11_SURFACE(surface)) return NULL;

	GdkDisplay *display = gdk_surface_get_display(surface);
	/* GDK_DISPLAY_XDISPLAY / GDK_SURFACE_XID are deprecated in GTK 4.18
	 * with no replacement: GTK is moving away from X11 specifics, but this
	 * workaround requires direct Xlib access and there is no GTK4 equivalent. */
	Display *xdpy = GDK_DISPLAY_XDISPLAY(display);
	Window xwin = GDK_SURFACE_XID(surface);

	/* Use XInternAtom directly instead of the deprecated
	 * gdk_x11_get_xatom_by_name_for_display() wrapper. */
	Atom xdnd_sel = XInternAtom(xdpy, "XdndSelection", False);
	Atom uri_list = XInternAtom(xdpy, "text/uri-list", False);
	Atom property = XInternAtom(xdpy, "SIRIL_XDND_DATA", False);

	XConvertSelection(xdpy, xdnd_sel, uri_list, property, xwin, time);
	XFlush(xdpy);

	XEvent xev;
	gboolean got = FALSE;
	gint64 deadline = g_get_monotonic_time() + X11_SELECTION_TIMEOUT_US;
	while (g_get_monotonic_time() < deadline) {
		if (XCheckTypedWindowEvent(xdpy, xwin, SelectionNotify, &xev)
		    && xev.xselection.selection == xdnd_sel) {
			got = TRUE;
			break;
		}
		g_usleep(1000);
	}
	if (!got || xev.xselection.property == None)
		return NULL;

	Atom actual_type;
	int actual_format;
	unsigned long nitems, bytes_after;
	unsigned char *data = NULL;
	if (XGetWindowProperty(xdpy, xwin, property, 0, 1 << 20, True,
				AnyPropertyType, &actual_type, &actual_format,
				&nitems, &bytes_after, &data) != Success || !data)
		return NULL;

	gchar *text = g_strndup((const gchar *) data, nitems);
	XFree(data);
	gchar **uris = g_uri_list_extract_uris(text);
	g_free(text);
	return uris;
}

static gboolean on_x11_drop_start(GtkEventControllerLegacy *controller,
		GdkEvent *event, gpointer user_data) {
	(void) user_data;
	if (gdk_event_get_event_type(event) != GDK_DROP_START)
		return FALSE;
	GdkDisplay *display = gdk_event_get_display(event);
	if (!GDK_IS_X11_DISPLAY(display))
		return FALSE;
	GdkDrop *drop = gdk_dnd_event_get_drop(event);
	if (!drop)
		return FALSE;

	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(controller));
	gchar **uris = x11_fetch_uri_list(widget, gdk_event_get_time(event));
	if (uris && uris[0]) {
		gchar *filename = g_filename_from_uri(uris[0], NULL, NULL);
		if (filename) {
			process_dropped_filename(filename);
			g_free(filename);
		}
		gdk_drop_finish(drop, GDK_ACTION_COPY);
	} else {
		gdk_drop_finish(drop, 0);
	}
	g_strfreev(uris);
	return TRUE; /* swallow so GtkDropTarget never tries the broken async read */
}
#endif /* HAVE_X11_DND_FALLBACK */

static void install_one_drop_target(GtkWidget *w) {
	if (!w) return;
	GtkDropTarget *target = gtk_drop_target_new(GDK_TYPE_FILE_LIST, GDK_ACTION_COPY);
	g_signal_connect(target, "drop", G_CALLBACK(on_drawingarea_drop), NULL);
	gtk_widget_add_controller(w, GTK_EVENT_CONTROLLER(target));

#if HAVE_X11_DND_FALLBACK
	if (GDK_IS_X11_DISPLAY(gtk_widget_get_display(w))) {
		GtkEventController *legacy = gtk_event_controller_legacy_new();
		gtk_event_controller_set_propagation_phase(legacy, GTK_PHASE_CAPTURE);
		g_signal_connect(legacy, "event", G_CALLBACK(on_x11_drop_start), NULL);
		gtk_widget_add_controller(w, legacy);
		/* Mark the widget so siril.css can suppress the :drop(active)
		 * outline that would otherwise stick after our intercepted drop. */
		gtk_widget_add_css_class(w, "siril-x11-dnd-fallback");
	}
#endif
}

void siril_drag_single_image_set_dest() {
	for (int i = 0; i < G_N_ELEMENTS(drawing_area); i++) {
		install_one_drop_target(
			GTK_WIDGET(gtk_builder_get_object(gui.builder, drawing_area[i])));
	}
}
