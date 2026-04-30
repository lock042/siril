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
 * Stub (no-op / stderr) implementations for the SirilGuiInterface vtable.
 *
 * This file defines the gui_iface global and is always compiled — in both
 * GUI and CLI builds.  The GUI build additionally compiles
 * gui/gui_iface_impl.c, which calls siril_register_gui_iface() from main.c
 * to replace these stubs with GTK-backed implementations.
 *
 * The CLI build never calls siril_register_gui_iface(), so these stubs
 * remain active for the lifetime of the process.
 */

#include <stdio.h>
#include "core/gui_iface.h"

/* ── Stub implementations ─────────────────────────────────────────────────── */

static void stub_set_progress(double fraction, const char *msg) {
	(void)fraction; (void)msg;
}

static void stub_set_busy(gboolean busy) {
	(void)busy;
}

static void stub_log_message(const char *msg, const char *color) {
	(void)color;
	/* The normal siril_log path already writes to stdout; suppress duplicate
	 * output here.  Enable if you need GUI-widget log echoed to stderr. */
	(void)msg;
}

static void stub_message_dialog(SirilMessageType type, const char *title,
                                const char *text) {
	const char *prefix;
	switch (type) {
		case SIRIL_MSG_WARNING:  prefix = "WARNING"; break;
		case SIRIL_MSG_ERROR:    prefix = "ERROR";   break;
		case SIRIL_MSG_QUESTION: prefix = "QUESTION"; break;
		default:                 prefix = "INFO";    break;
	}
	fprintf(stderr, "siril [%s] %s: %s\n", prefix,
	        title ? title : "", text ? text : "");
}

static gboolean stub_confirm_dialog(const char *title, const char *msg,
                                    const char *button_accept) {
	(void)title; (void)msg; (void)button_accept;
	/* Auto-accept in headless/CLI mode so scripts are not blocked. */
	return TRUE;
}

static void stub_open_dialog(const char *id)  { (void)id; }
static void stub_close_dialog(const char *id) { (void)id; }

static void stub_redraw_image(SirilRedrawType remap)       { (void)remap; }
static void stub_redraw_image_async(SirilRedrawType remap) { (void)remap; }
static void stub_redraw_image_sync(SirilRedrawType remap)  { (void)remap; }
static void stub_delete_selection(void) {}

static void stub_on_sequence_opened(void) {}
static void stub_on_image_loaded(void) {}
static void stub_on_image_closed(void) {}

static void stub_show_panel(const char *panel_name, gboolean visible) {
	(void)panel_name; (void)visible;
}

static void stub_update_status_bar(void) {}
static void stub_update_menu_state(void) {}
static void stub_on_geometry_changed(void) {}
static void stub_on_mask_state_changed(void) {}
static void stub_on_crop_complete(void) {}
static void stub_on_stats_ready(void) {}
static void stub_on_photometry_changed(void) {}
static void stub_show_siril_plot(gpointer spl_data) { (void)spl_data; }
static void stub_update_star_list(psf_star **stars, gboolean update_psf_list,
                                  gboolean wait) {
	(void)stars; (void)update_psf_list; (void)wait;
}

/* ── Global instance ──────────────────────────────────────────────────────── */

SirilGuiInterface gui_iface = {
	.set_progress         = stub_set_progress,
	.set_busy             = stub_set_busy,
	.log_message          = stub_log_message,
	.message_dialog       = stub_message_dialog,
	.confirm_dialog       = stub_confirm_dialog,
	.open_dialog          = stub_open_dialog,
	.close_dialog         = stub_close_dialog,
	.redraw_image         = stub_redraw_image,
	.redraw_image_async   = stub_redraw_image_async,
	.redraw_image_sync    = stub_redraw_image_sync,
	.delete_selection     = stub_delete_selection,
	.on_sequence_opened   = stub_on_sequence_opened,
	.on_image_loaded      = stub_on_image_loaded,
	.on_image_closed      = stub_on_image_closed,
	.show_panel           = stub_show_panel,
	.update_status_bar    = stub_update_status_bar,
	.update_menu_state    = stub_update_menu_state,
	.on_geometry_changed  = stub_on_geometry_changed,
	.on_mask_state_changed= stub_on_mask_state_changed,
	.on_crop_complete     = stub_on_crop_complete,
	.on_stats_ready       = stub_on_stats_ready,
	.on_photometry_changed= stub_on_photometry_changed,
	.show_siril_plot      = stub_show_siril_plot,
	.update_star_list     = stub_update_star_list,
};
