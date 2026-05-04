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
#include "core/processing.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/callbacks.h"
#include "gui/geometry.h"
#include "gui/gui_state.h"
#include "gui/plot.h"
#include "gui/PSF_list.h"
#include "gui/registration_preview.h"
#include "gui/sequence_list.h"
#include "gui/siril_plot.h"
#include "gui/stacking.h"
#include "gui/utils.h"
#include "io/single_image.h"

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
	queue_message_dialog(to_gtk_msg_type(type), title, text);
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

/* Called after stacking completes successfully on the GTK main thread.
 * Consolidates all display-state update calls that were previously scattered
 * through end_stacking() in stacking.c. */
static void impl_on_stack_complete(void) {
	clear_stars_list(TRUE);
	initialize_display_mode();
	sliders_mode_set_state(gui.sliders);
	/* Reader lock guards set_cutoff_sliders_max_values() which reads
	 * gfit->type/orig_bitpix.  The hi/lo assignment writes keyword fields
	 * that were set on the worker thread by notify_gfit_data_modified(). */
	g_rw_lock_reader_lock(&gfit->rwlock);
	display_filename();
	gui_function(set_precision_switch, NULL);
	set_cutoff_sliders_max_values();
	gfit->keywords.hi = gui.hi;
	gfit->keywords.lo = gui.lo;
	g_rw_lock_writer_unlock(&gfit->rwlock);
	gfit_modified_update_gui();
	set_display_mode();
	gui_function(update_MenuItem, NULL);
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	sequence_list_change_current();
	update_stack_interface(TRUE);
}

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

/* ── Group I: Statistics ─────────────────────────────────────────────────── */

static void impl_on_stats_ready(void) {
	computeStat();
	siril_open_dialog("StatWindow");
}

/* ── Group J: Photometry ─────────────────────────────────────────────────── */

static void impl_on_photometry_changed(void) {
	drawPlot();
	notify_new_photometry();
	redraw(REDRAW_OVERLAY);
}

static void impl_show_siril_plot(gpointer spl_data) {
	siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
}

/* ── Group K: Star list ──────────────────────────────────────────────────── */

static void impl_update_star_list(psf_star **stars, gboolean update_psf_list,
                                  gboolean wait) {
	update_star_list(stars, update_psf_list, wait);
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
	gui_iface.on_stack_complete     = impl_on_stack_complete;
	gui_iface.show_panel            = impl_show_panel;
	gui_iface.update_status_bar     = impl_update_status_bar;
	gui_iface.update_menu_state     = impl_update_menu_state;
	gui_iface.on_geometry_changed   = impl_on_geometry_changed;
	gui_iface.on_mask_state_changed = impl_on_mask_state_changed;
	gui_iface.on_crop_complete      = impl_on_crop_complete;
	gui_iface.on_stats_ready        = impl_on_stats_ready;
	gui_iface.on_photometry_changed = impl_on_photometry_changed;
	gui_iface.show_siril_plot       = impl_show_siril_plot;
	gui_iface.update_star_list      = impl_update_star_list;
}
