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

/*
 * GTK-backed implementations of the SirilGuiInterface vtable.
 *
 * Compiled only in GUI builds (listed under src_files_gui in meson.build).
 * main.c calls siril_register_gui_iface() after load_ui_files() to replace
 * the no-op stubs from gui_iface_stubs.c with these implementations.
 */

#include <gtk/gtk.h>
#include "core/gui_iface.h"
#include "core/proto.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/callbacks.h"
#include "gui/geometry.h"
#include "gui/utils.h"

/* ── Group A: Progress ───────────────────────────────────────────────────── */

static void impl_set_progress(double fraction, const char *msg) {
	set_progress_bar_data(msg, fraction);
}

static void impl_set_busy(gboolean busy) {
	set_cursor_waiting(busy);
}

/* ── Group B: Logging ────────────────────────────────────────────────────── */

static void impl_log_message(const char *msg, const char *color) {
	gui_log_message(msg, color);
}

/* ── Group C: Dialogs ────────────────────────────────────────────────────── */

static GtkMessageType to_gtk_msg_type(SirilMessageType type) {
	switch (type) {
		case SIRIL_MSG_WARNING:  return GTK_MESSAGE_WARNING;
		case SIRIL_MSG_ERROR:    return GTK_MESSAGE_ERROR;
		case SIRIL_MSG_QUESTION: return GTK_MESSAGE_QUESTION;
		default:                 return GTK_MESSAGE_INFO;
	}
}

static void impl_message_dialog(SirilMessageType type, const char *title,
                                const char *text) {
	siril_message_dialog(to_gtk_msg_type(type), (char *)title, (char *)text);
}

static gboolean impl_confirm_dialog(const char *title, const char *msg,
                                    const char *button_accept) {
	return siril_confirm_dialog((gchar *)title, (gchar *)msg,
	                            (gchar *)button_accept);
}

static void impl_open_dialog(const char *id) {
	siril_open_dialog((gchar *)id);
}

static void impl_close_dialog(const char *id) {
	siril_close_dialog((gchar *)id);
}

/* ── Group D: Image display ──────────────────────────────────────────────── */

static void impl_redraw_image(SirilRedrawType remap) {
	redraw((remap_type)remap);
}

static void impl_redraw_image_async(SirilRedrawType remap) {
	queue_redraw((remap_type)remap);
}

static void impl_redraw_image_sync(SirilRedrawType remap) {
	queue_redraw_and_wait_for_it((remap_type)remap);
}

static void impl_delete_selection(void) {
	delete_selected_area();
}

/* ── Groups E, F: no-op placeholders until phases 4–5 wire them up ── */

static void impl_on_sequence_opened(void) {}
static void impl_on_image_loaded(void) {}
static void impl_on_image_closed(void) {}

/* ── Group F: Panel / tab switching ─────────────────────────────────────── */

static void impl_show_panel(const char *panel_name, gboolean visible) {
	if (g_strcmp0(panel_name, "output_logs") == 0) {
		if (visible)
			control_window_switch_to_tab(OUTPUT_LOGS);
		return;
	}
	(void)visible;
}

/* ── Group G: Misc GUI state ─────────────────────────────────────────────── */

static void impl_update_status_bar(void) {
	gui_function(update_zoom_label_idle, NULL);
}

static void impl_update_menu_state(void) {
	gui_function(update_MenuItem, NULL);
}

/* ── Group H: Geometry / ROI / Mask state ────────────────────────────────── */

static void impl_on_geometry_changed(void) {
	on_clear_roi();
}

static void impl_on_mask_state_changed(void) {
	show_or_hide_mask_tab();
}

static void impl_on_crop_complete(void) {
	gui_function(crop_gui_updates, NULL);
}

/* ── Registration ────────────────────────────────────────────────────────── */

void siril_register_gui_iface(void) {
	gui_iface.set_progress          = impl_set_progress;
	gui_iface.set_busy              = impl_set_busy;
	gui_iface.log_message           = impl_log_message;
	gui_iface.message_dialog        = impl_message_dialog;
	gui_iface.confirm_dialog        = impl_confirm_dialog;
	gui_iface.open_dialog           = impl_open_dialog;
	gui_iface.close_dialog          = impl_close_dialog;
	gui_iface.redraw_image          = impl_redraw_image;
	gui_iface.redraw_image_async    = impl_redraw_image_async;
	gui_iface.redraw_image_sync     = impl_redraw_image_sync;
	gui_iface.delete_selection      = impl_delete_selection;
	gui_iface.on_sequence_opened    = impl_on_sequence_opened;
	gui_iface.on_image_loaded       = impl_on_image_loaded;
	gui_iface.on_image_closed       = impl_on_image_closed;
	gui_iface.show_panel            = impl_show_panel;
	gui_iface.update_status_bar     = impl_update_status_bar;
	gui_iface.update_menu_state     = impl_update_menu_state;
	gui_iface.on_geometry_changed   = impl_on_geometry_changed;
	gui_iface.on_mask_state_changed = impl_on_mask_state_changed;
	gui_iface.on_crop_complete      = impl_on_crop_complete;
}
