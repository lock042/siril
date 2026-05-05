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
static gboolean stub_roi_operation_supports(void) { return FALSE; }
static gpointer stub_get_roi_fit(void) { return NULL; }
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
static void stub_apply_display_icc_compensation(gpointer fit) { (void)fit; }

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

/* 6.1–6.4 additions */
static void stub_on_sequence_closed(gboolean loading_next) { (void)loading_next; }
static void stub_populate_seq_combo(const char *realname) { (void)realname; }
static void stub_set_gui_cwd(void) {}
static void stub_refresh_keywords_dialog(void) {}
static void stub_set_wcs_overlay(gboolean show) { (void)show; }
static void stub_launch_clipboard_survey(void) {}
static void stub_clear_log_buffer(void) {}
static void stub_update_spin_cpu(void) {}
static void stub_new_selection_zone(void) {}
static int  stub_get_active_vport(void) { return 0; }
static void stub_console_set_status(const char *msg, int line) { (void)msg; (void)line; }
static void stub_console_clear_status(void) {}
static void stub_end_script_gui(void) {}
static gboolean stub_get_star_follow_state(void) { return FALSE; }
static void stub_show_command_help(void) {}

/* S – Pixel-math status */
static void stub_update_pixel_math_status(int ret) { (void)ret; }

/* SC – Annotation display */
static void stub_activate_annotation_display(void) {}

/* U – ICC status */
static void stub_update_icc_status_icon(gpointer fit, gboolean active) { (void)fit; (void)active; }
static gboolean stub_get_gamut_check_active(void) { return FALSE; }

/* V – Registration panel status */
static void stub_update_registration_status(const gchar *msg) { (void)msg; }

/* R – Python bridge UI */
static void stub_update_single_image_display(void) {}
static void stub_seq_redisplay_frame(int index) { (void)index; }
static void stub_set_poly_drawing(guint32 color, gboolean fill) { (void)color; (void)fill; }
static void stub_set_rendering_mode(int mode) { (void)mode; }
static void stub_set_channels_linked(gboolean state) { (void)state; }
static void stub_set_sliders_mode(int mode) { (void)mode; }
static void stub_set_cutoff_values(int lo, int hi) { (void)lo; (void)hi; }
static void stub_update_zoom_label(void) {}
static double stub_get_zoom_value(void) { return 1.0; }
static int stub_activate_action(const char *name, gboolean appmap) { (void)name; (void)appmap; return ACTION_SUCCESS; }
static void stub_reset_display_offset(void) {}

static void stub_enable_display_mode_menu(void) {}
static void stub_switch_to_tab(int tab) { (void)tab; }

/* SG – Miscellaneous single-file accesses */
static void stub_set_last_opened_filetype(int type) { (void)type; }
static void stub_free_reference_image_display(void) {}
static psf_star *stub_get_qphot_result(void) { return NULL; }

/* SF – Python IPC display state */
static int stub_get_channel_for_vport(void) { return 0; }
static int stub_get_rendering_mode(void) { return 0; }  /* LINEAR */
static gboolean stub_get_channels_linked(void) { return TRUE; }
static void stub_get_display_offset(double *x, double *y) {
	if (x) *x = 0.0;
	if (y) *y = 0.0;
}
static void stub_set_display_offset(double x, double y) { (void)x; (void)y; }
static void stub_set_zoom_value(double zoom) { (void)zoom; }
static GSList *stub_get_user_polygons(void) { return NULL; }
static void stub_add_user_polygon_to_list(gpointer polygon) { (void)polygon; }

/* SD – Display range state */
static void stub_get_display_lo_hi(int *lo, int *hi) {
	if (lo) *lo = 0;
	if (hi) *hi = 0xFFFF;
}
static int stub_get_sliders_mode(void) { return 0; } /* MIPSLOHI */
static void stub_update_display_range_after_load(int sliders, int lo, int hi) {
	(void)sliders; (void)lo; (void)hi;
}

/* Phase 3 – Display / slider state */
static void stub_sliders_mode_set_state(int mode) { (void)mode; }
static void stub_set_cutoff_sliders_max_values(void) {}
static void stub_set_cutoff_sliders_values(void) {}
static void stub_update_display_mode_state(void) {}
static void stub_compute_histo_for_fit(gpointer fit) { (void)fit; }
static void stub_refresh_histogram_if_visible(void) {}
static void stub_fill_sequence_list(gpointer seq, int layer, gboolean as_idle) {
	(void)seq; (void)layer; (void)as_idle;
}
static void stub_sequence_list_change_current(void) {}
static void stub_enable_view_reference_checkbox(gboolean status) { (void)status; }
static void stub_close_tab(void) {}
static void stub_init_right_tab(void) {}
static void stub_initialize_display_mode(void) {}
static void stub_display_filename(void) {}
static void stub_update_display_fwhm(void) {}
static void stub_update_prepro_interface(gboolean allow_debayer) { (void)allow_debayer; }
static void stub_adjust_sellabel(void) {}
static void stub_adjust_reginfo(void) {}
static void stub_adjust_refimage(int n) { (void)n; }
static int  stub_set_layers_for_registration(void) { return 0; }
static void stub_set_precision_switch(void) {}
static void stub_set_GUI_CAMERA(void) {}
static void stub_update_menu_item(void) {}
static void stub_update_seqlist(int layer) { (void)layer; }
static void stub_update_sequence_overlay_async(void) {}
static void stub_ensure_seqlist_dialog_closed(void) {}
static void stub_copy_roi_into_gfit(void) {}
static void stub_lock_roi_mutex(void) {}
static void stub_unlock_roi_mutex(void) {}
static void stub_show_or_hide_mask_tab(void) {}
static void stub_show_or_hide_mask_tab_async(void) {}
static int  stub_number_of_dialogs(void) { return 0; }
static void stub_clear_previews(void) {}
static int  stub_toggle_remixer_window_visibility(int invocation,
                                                   gpointer fit_left,
                                                   gpointer fit_right) {
	(void)invocation; (void)fit_left; (void)fit_right; return 0;
}
static gboolean stub_heif_dialog(gpointer heif, uint32_t *selected_image) {
	(void)heif; (void)selected_image; return FALSE;
}

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
	.roi_operation_supports      = stub_roi_operation_supports,
	.get_roi_fit                 = stub_get_roi_fit,
	.get_roi_selection           = stub_get_roi_selection,
	.clear_roi                   = stub_clear_roi,
	.restore_roi                 = stub_restore_roi,
	.reset_display_transform     = stub_reset_display_transform,
	.is_preview_active           = stub_is_preview_active,
	.hide_preview                = stub_hide_preview,
	.copy_gfit_to_backup         = stub_copy_gfit_to_backup,
	.copy_gfit_icc_to_backup     = stub_copy_gfit_icc_to_backup,
	.check_icc_identical_to_monitor = stub_check_icc_identical_to_monitor,
	.set_source_information            = stub_set_source_information,
	.apply_display_icc_compensation    = stub_apply_display_icc_compensation,
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
	.on_sequence_closed          = stub_on_sequence_closed,
	.populate_seq_combo          = stub_populate_seq_combo,
	.set_gui_cwd                 = stub_set_gui_cwd,
	.refresh_keywords_dialog     = stub_refresh_keywords_dialog,
	.set_wcs_overlay             = stub_set_wcs_overlay,
	.launch_clipboard_survey     = stub_launch_clipboard_survey,
	.clear_log_buffer            = stub_clear_log_buffer,
	.update_spin_cpu             = stub_update_spin_cpu,
	.new_selection_zone          = stub_new_selection_zone,
	.get_active_vport            = stub_get_active_vport,
	.console_set_status          = stub_console_set_status,
	.console_clear_status        = stub_console_clear_status,
	.end_script_gui              = stub_end_script_gui,
	.get_star_follow_state       = stub_get_star_follow_state,
	.show_command_help           = stub_show_command_help,
	.update_pixel_math_status    = stub_update_pixel_math_status,
	.activate_annotation_display = stub_activate_annotation_display,
	.update_icc_status_icon      = stub_update_icc_status_icon,
	.get_gamut_check_active      = stub_get_gamut_check_active,
	.update_registration_status  = stub_update_registration_status,
	.update_single_image_display = stub_update_single_image_display,
	.seq_redisplay_frame         = stub_seq_redisplay_frame,
	.set_poly_drawing            = stub_set_poly_drawing,
	.set_rendering_mode          = stub_set_rendering_mode,
	.set_channels_linked         = stub_set_channels_linked,
	.set_sliders_mode            = stub_set_sliders_mode,
	.set_cutoff_values           = stub_set_cutoff_values,
	.update_zoom_label           = stub_update_zoom_label,
	.get_zoom_value              = stub_get_zoom_value,
	.activate_action             = stub_activate_action,
	.reset_display_offset        = stub_reset_display_offset,
	.enable_display_mode_menu        = stub_enable_display_mode_menu,
	.switch_to_tab                   = stub_switch_to_tab,
	.set_last_opened_filetype        = stub_set_last_opened_filetype,
	.free_reference_image_display    = stub_free_reference_image_display,
	.get_qphot_result                = stub_get_qphot_result,
	.get_channel_for_vport           = stub_get_channel_for_vport,
	.get_rendering_mode              = stub_get_rendering_mode,
	.get_channels_linked             = stub_get_channels_linked,
	.get_display_offset              = stub_get_display_offset,
	.set_display_offset              = stub_set_display_offset,
	.set_zoom_value                  = stub_set_zoom_value,
	.get_user_polygons               = stub_get_user_polygons,
	.add_user_polygon_to_list        = stub_add_user_polygon_to_list,
	.get_display_lo_hi               = stub_get_display_lo_hi,
	.get_sliders_mode                = stub_get_sliders_mode,
	.update_display_range_after_load = stub_update_display_range_after_load,

	/* Phase 3 */
	.sliders_mode_set_state          = stub_sliders_mode_set_state,
	.set_cutoff_sliders_max_values   = stub_set_cutoff_sliders_max_values,
	.set_cutoff_sliders_values       = stub_set_cutoff_sliders_values,
	.update_display_mode_state       = stub_update_display_mode_state,
	.compute_histo_for_fit           = stub_compute_histo_for_fit,
	.refresh_histogram_if_visible    = stub_refresh_histogram_if_visible,
	.fill_sequence_list              = stub_fill_sequence_list,
	.sequence_list_change_current    = stub_sequence_list_change_current,
	.enable_view_reference_checkbox  = stub_enable_view_reference_checkbox,
	.close_tab                       = stub_close_tab,
	.init_right_tab                  = stub_init_right_tab,
	.initialize_display_mode         = stub_initialize_display_mode,
	.display_filename                = stub_display_filename,
	.update_display_fwhm             = stub_update_display_fwhm,
	.update_prepro_interface         = stub_update_prepro_interface,
	.adjust_sellabel                 = stub_adjust_sellabel,
	.adjust_reginfo                  = stub_adjust_reginfo,
	.adjust_refimage                 = stub_adjust_refimage,
	.set_layers_for_registration     = stub_set_layers_for_registration,
	.set_precision_switch            = stub_set_precision_switch,
	.set_GUI_CAMERA                  = stub_set_GUI_CAMERA,
	.update_menu_item                = stub_update_menu_item,
	.update_seqlist                  = stub_update_seqlist,
	.update_sequence_overlay_async   = stub_update_sequence_overlay_async,
	.ensure_seqlist_dialog_closed    = stub_ensure_seqlist_dialog_closed,
	.copy_roi_into_gfit              = stub_copy_roi_into_gfit,
	.lock_roi_mutex                  = stub_lock_roi_mutex,
	.unlock_roi_mutex                = stub_unlock_roi_mutex,
	.show_or_hide_mask_tab           = stub_show_or_hide_mask_tab,
	.show_or_hide_mask_tab_async     = stub_show_or_hide_mask_tab_async,
	.number_of_dialogs               = stub_number_of_dialogs,
	.clear_previews                  = stub_clear_previews,
	.toggle_remixer_window_visibility = stub_toggle_remixer_window_visibility,
	.heif_dialog                     = stub_heif_dialog,
};
