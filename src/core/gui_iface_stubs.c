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
static gboolean stub_is_dialog_open(void) { return FALSE; }

static void stub_redraw_image(SirilRedrawType remap)       { (void)remap; }
static void stub_redraw_image_async(SirilRedrawType remap) { (void)remap; }
static void stub_redraw_image_sync(SirilRedrawType remap)  { (void)remap; }
static void stub_delete_selection(void) {}
static void stub_queue_redraw_mask(void) {}

static void stub_on_sequence_opened(void) {}
static void stub_on_image_loaded(void) {}
static void stub_on_image_closed(void) {}
static void stub_on_stack_complete(void) {}
static void stub_update_sequences_list(const char *seqname) { (void)seqname; }

static void stub_show_panel(const char *panel_name, gboolean visible) {
	(void)panel_name; (void)visible;
}
static void stub_script_widgets_enable(gboolean enable) { (void)enable; }
static void stub_set_seq_browser_active(gboolean active) { (void)active; }

static void stub_update_status_bar(void) {}
static void stub_update_menu_state(void) {}
static void stub_set_suppress_redraws(gboolean suppress) { (void)suppress; }
static void stub_populate_roi(void) {}
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
static void stub_clear_star_list(void) {}
static int stub_get_reg_layer(void) { return -1; }

static void stub_execute_idle_sync(GSourceFunc func, gpointer data) { func(data); }
static GPid stub_select_child_process(GSList *children) { (void)children; return (GPid)0; }

/* D additions */
static void stub_invalidate_histogram(void) {}
static void stub_update_histogram(void) {}
static void stub_redraw_mask_idle(void) {}

/* G additions */
static void stub_on_channel_count_changed(void) {}
static void stub_on_precision_changed(void) {}

/* H additions */
static gboolean stub_roi_is_active(void) { return FALSE; }
static void stub_get_roi_selection(rectangle *rect) { (void)rect; }
static void stub_clear_roi(void) {}
static void stub_restore_roi(const rectangle *rect) { (void)rect; }
static void stub_reset_display_transform(void) {}

/* N – Preview */
static gboolean stub_is_preview_active(void) { return FALSE; }
static void stub_hide_preview(void) {}
static void stub_copy_gfit_to_backup(void) {}
static void stub_copy_gfit_icc_to_backup(void) {}

/* O – ICC info */
static void stub_check_icc_identical_to_monitor(void) {}
static void stub_set_source_information(void) {}

/* Additional slots (5.6–5.10) */
static void stub_data_dialog(SirilMessageType type, const char *title,
                             const char *text, const char *data) {
	stub_message_dialog(type, title, text);
	if (data) fprintf(stderr, "%s\n", data);
}
static void stub_redraw_previews(void) {}
static void stub_open_single_image_from_gfit(void) {}
static void stub_update_mem_usage(guint64 used_bytes) { (void)used_bytes; }
static void stub_update_disk_space(gint64 space_bytes, const char *label_id) {
	(void)space_bytes; (void)label_id;
}
static void stub_update_mask_enable(gboolean state) { (void)state; }
static void stub_set_display_range(int lo, int hi) { (void)lo; (void)hi; }
static void stub_check_gaia_status(void) {}
static void stub_trigger_gaia_check(void) {}

/* Steps 5.11–5.15 */
static void stub_remap_all_vports(void) {}
static void stub_quit_application(void) { exit(0); }
static void stub_refresh_script_menu(void) {}
static void stub_clear_backup(void) {}
static void stub_livestacking_setup_gui(gboolean has_dark, gboolean has_flat,
                                        int reg_type) {
	(void)has_dark; (void)has_flat; (void)reg_type;
}
static void stub_livestacking_teardown_gui(void) {}

/* ── Global instance ──────────────────────────────────────────────────────── */

SirilGuiInterface gui_iface = {
	.set_progress         = stub_set_progress,
	.set_busy             = stub_set_busy,
	.log_message          = stub_log_message,
	.message_dialog       = stub_message_dialog,
	.confirm_dialog       = stub_confirm_dialog,
	.open_dialog            = stub_open_dialog,
	.close_dialog           = stub_close_dialog,
	.is_dialog_open         = stub_is_dialog_open,
	.redraw_image           = stub_redraw_image,
	.redraw_image_async     = stub_redraw_image_async,
	.redraw_image_sync      = stub_redraw_image_sync,
	.delete_selection       = stub_delete_selection,
	.queue_redraw_mask      = stub_queue_redraw_mask,
	.on_sequence_opened     = stub_on_sequence_opened,
	.on_image_loaded        = stub_on_image_loaded,
	.on_image_closed        = stub_on_image_closed,
	.on_stack_complete      = stub_on_stack_complete,
	.update_sequences_list  = stub_update_sequences_list,
	.show_panel             = stub_show_panel,
	.script_widgets_enable  = stub_script_widgets_enable,
	.set_seq_browser_active = stub_set_seq_browser_active,
	.update_status_bar      = stub_update_status_bar,
	.update_menu_state      = stub_update_menu_state,
	.set_suppress_redraws   = stub_set_suppress_redraws,
	.populate_roi           = stub_populate_roi,
	.on_geometry_changed    = stub_on_geometry_changed,
	.on_mask_state_changed  = stub_on_mask_state_changed,
	.on_crop_complete       = stub_on_crop_complete,
	.on_stats_ready         = stub_on_stats_ready,
	.on_photometry_changed  = stub_on_photometry_changed,
	.show_siril_plot        = stub_show_siril_plot,
	.update_star_list       = stub_update_star_list,
	.clear_star_list        = stub_clear_star_list,
	.get_reg_layer          = stub_get_reg_layer,
	.execute_idle_sync           = stub_execute_idle_sync,
	.select_child_process        = stub_select_child_process,
	.invalidate_histogram        = stub_invalidate_histogram,
	.update_histogram            = stub_update_histogram,
	.redraw_mask_idle            = stub_redraw_mask_idle,
	.on_channel_count_changed    = stub_on_channel_count_changed,
	.on_precision_changed        = stub_on_precision_changed,
	.roi_is_active               = stub_roi_is_active,
	.get_roi_selection           = stub_get_roi_selection,
	.clear_roi                   = stub_clear_roi,
	.restore_roi                 = stub_restore_roi,
	.reset_display_transform     = stub_reset_display_transform,
	.is_preview_active           = stub_is_preview_active,
	.hide_preview                = stub_hide_preview,
	.copy_gfit_to_backup         = stub_copy_gfit_to_backup,
	.copy_gfit_icc_to_backup     = stub_copy_gfit_icc_to_backup,
	.check_icc_identical_to_monitor = stub_check_icc_identical_to_monitor,
	.set_source_information      = stub_set_source_information,
	.data_dialog                 = stub_data_dialog,
	.redraw_previews             = stub_redraw_previews,
	.open_single_image_from_gfit = stub_open_single_image_from_gfit,
	.update_mem_usage            = stub_update_mem_usage,
	.update_disk_space           = stub_update_disk_space,
	.update_mask_enable          = stub_update_mask_enable,
	.set_display_range           = stub_set_display_range,
	.check_gaia_status           = stub_check_gaia_status,
	.trigger_gaia_check          = stub_trigger_gaia_check,
	.remap_all_vports            = stub_remap_all_vports,
	.quit_application            = stub_quit_application,
	.refresh_script_menu         = stub_refresh_script_menu,
	.clear_backup                = stub_clear_backup,
	.livestacking_setup_gui      = stub_livestacking_setup_gui,
	.livestacking_teardown_gui   = stub_livestacking_teardown_gui,
};
