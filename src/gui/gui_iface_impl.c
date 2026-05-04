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
#include "core/icc_profile.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/histogram.h"
#include "gui/icc_profile.h"
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
#include "gui/siril_preview.h"
#include "gui/registration.h"
#include "gui/script_menu.h"
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

static gboolean impl_is_dialog_open(void) {
	return is_an_image_processing_dialog_opened();
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

static void impl_queue_redraw_mask(void) {
	queue_redraw_mask();
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

static void impl_update_sequences_list(const char *seqname) {
	update_sequences_list(seqname);
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

static void impl_script_widgets_enable(gboolean enable) {
	script_widgets_enable(enable);
}

static void impl_set_seq_browser_active(gboolean active) {
	gui_function(set_seq_browser_active, GINT_TO_POINTER(active));
}

/* ── Group G: Misc GUI state ─────────────────────────────────────────────── */

static void impl_update_status_bar(void) {
	gui_function(update_zoom_label_idle, NULL);
}

static void impl_update_menu_state(void) {
	gui_function(update_MenuItem, NULL);
}

static gboolean set_display_mode_menu_sensitive_idle(gpointer p) {
	gtk_widget_set_sensitive(
		GTK_WIDGET(gtk_builder_get_object(gui.builder, "menu_display_button")),
		GPOINTER_TO_INT(p));
	return FALSE;
}

static void impl_set_suppress_redraws(gboolean suppress) {
	if (suppress) {
		g_atomic_int_set(&gui.suppress_drawarea_redraw, 1);
		siril_add_idle(set_display_mode_menu_sensitive_idle, GINT_TO_POINTER(FALSE));
	} else {
		g_atomic_int_set(&gui.suppress_drawarea_redraw, 0);
	}
}

static void impl_populate_roi(void) {
	populate_roi();
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

/* ── Group M: Thread utilities ───────────────────────────────────────────── */

static void impl_execute_idle_sync(GSourceFunc func, gpointer data) {
	execute_idle_and_wait_for_it(func, data);
}

static GPid impl_select_child_process(GSList *children) {
	return show_child_process_selection_dialog(children);
}

/* ── Group D additions: Histogram / image modification state ─────────────── */

static void impl_invalidate_histogram(void) {
	invalidate_gfit_histogram();
}

static void impl_update_histogram(void) {
	update_gfit_histogram_if_needed();
}

static void impl_redraw_mask_idle(void) {
	redraw_mask_idle(NULL);
}

/* ── Group G additions: Channel / precision display state ────────────────── */

static void impl_on_channel_count_changed(void) {
	gui_function(close_tab, NULL);
}

static void impl_on_precision_changed(void) {
	gui_function(set_precision_switch, NULL);
}

/* ── Group H additions: ROI state getters/setters ────────────────────────── */

static gboolean impl_roi_is_active(void) {
	return gui.roi.active;
}

static void impl_get_roi_selection(rectangle *rect) {
	memcpy(rect, &gui.roi.selection, sizeof(rectangle));
}

static void impl_clear_roi(void) {
	on_clear_roi();
}

static void impl_restore_roi(const rectangle *rect) {
	memcpy(&com.selection, rect, sizeof(rectangle));
	on_set_roi();
}

static void impl_reset_display_transform(void) {
	lock_display_transform();
	if (gui.icc.proofing_transform) {
		cmsDeleteTransform(gui.icc.proofing_transform);
		gui.icc.proofing_transform = NULL;
	}
	unlock_display_transform();
}

/* ── Group N: Preview / backup state ─────────────────────────────────────── */

static gboolean impl_is_preview_active(void) {
	return is_preview_active();
}

static void impl_hide_preview(void) {
	siril_preview_hide();
}

static void impl_copy_gfit_to_backup(void) {
	copy_gfit_to_backup();
}

static void impl_copy_gfit_icc_to_backup(void) {
	copy_gfit_icc_to_backup();
}

/* ── Group O: ICC information callbacks ──────────────────────────────────── */

static void impl_check_icc_identical_to_monitor(void) {
	check_gfit_profile_identical_to_monitor();
}

static void impl_set_source_information(void) {
	set_source_information();
}

/* ── Groups K, L: Star list / Registration state ─────────────────────────── */

static void impl_update_star_list(psf_star **stars, gboolean update_psf_list,
                                  gboolean wait) {
	update_star_list(stars, update_psf_list, wait);
}

static void impl_clear_star_list(void) {
	clear_stars_list(FALSE);
}

static int impl_get_reg_layer(void) {
	return get_registration_layer_from_GUI(&com.seq);
}

/* ── Registration ────────────────────────────────────────────────────────── */

void siril_register_gui_iface(void) {
	gui_iface.set_progress          = impl_set_progress;
	gui_iface.set_busy              = impl_set_busy;
	gui_iface.log_message           = impl_log_message;
	gui_iface.message_dialog        = impl_message_dialog;
	gui_iface.confirm_dialog        = impl_confirm_dialog;
	gui_iface.open_dialog            = impl_open_dialog;
	gui_iface.close_dialog           = impl_close_dialog;
	gui_iface.is_dialog_open         = impl_is_dialog_open;
	gui_iface.redraw_image           = impl_redraw_image;
	gui_iface.redraw_image_async     = impl_redraw_image_async;
	gui_iface.redraw_image_sync      = impl_redraw_image_sync;
	gui_iface.delete_selection       = impl_delete_selection;
	gui_iface.queue_redraw_mask      = impl_queue_redraw_mask;
	gui_iface.on_sequence_opened     = impl_on_sequence_opened;
	gui_iface.on_image_loaded        = impl_on_image_loaded;
	gui_iface.on_image_closed        = impl_on_image_closed;
	gui_iface.on_stack_complete      = impl_on_stack_complete;
	gui_iface.update_sequences_list  = impl_update_sequences_list;
	gui_iface.show_panel             = impl_show_panel;
	gui_iface.script_widgets_enable  = impl_script_widgets_enable;
	gui_iface.set_seq_browser_active = impl_set_seq_browser_active;
	gui_iface.update_status_bar      = impl_update_status_bar;
	gui_iface.update_menu_state      = impl_update_menu_state;
	gui_iface.set_suppress_redraws   = impl_set_suppress_redraws;
	gui_iface.populate_roi           = impl_populate_roi;
	gui_iface.on_geometry_changed    = impl_on_geometry_changed;
	gui_iface.on_mask_state_changed  = impl_on_mask_state_changed;
	gui_iface.on_crop_complete       = impl_on_crop_complete;
	gui_iface.on_stats_ready         = impl_on_stats_ready;
	gui_iface.on_photometry_changed  = impl_on_photometry_changed;
	gui_iface.show_siril_plot        = impl_show_siril_plot;
	gui_iface.update_star_list       = impl_update_star_list;
	gui_iface.clear_star_list        = impl_clear_star_list;
	gui_iface.get_reg_layer          = impl_get_reg_layer;
	gui_iface.execute_idle_sync           = impl_execute_idle_sync;
	gui_iface.select_child_process        = impl_select_child_process;
	gui_iface.invalidate_histogram        = impl_invalidate_histogram;
	gui_iface.update_histogram            = impl_update_histogram;
	gui_iface.redraw_mask_idle            = impl_redraw_mask_idle;
	gui_iface.on_channel_count_changed    = impl_on_channel_count_changed;
	gui_iface.on_precision_changed        = impl_on_precision_changed;
	gui_iface.roi_is_active               = impl_roi_is_active;
	gui_iface.get_roi_selection           = impl_get_roi_selection;
	gui_iface.clear_roi                   = impl_clear_roi;
	gui_iface.restore_roi                 = impl_restore_roi;
	gui_iface.reset_display_transform     = impl_reset_display_transform;
	gui_iface.is_preview_active           = impl_is_preview_active;
	gui_iface.hide_preview                = impl_hide_preview;
	gui_iface.copy_gfit_to_backup         = impl_copy_gfit_to_backup;
	gui_iface.copy_gfit_icc_to_backup     = impl_copy_gfit_icc_to_backup;
	gui_iface.check_icc_identical_to_monitor = impl_check_icc_identical_to_monitor;
	gui_iface.set_source_information      = impl_set_source_information;
}
