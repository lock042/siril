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
#include "algos/background_extraction.h"
#include "algos/ccd-inspector.h"
#include "gui/cut.h"
#include "gui/dialogs.h"
#include "gui/histogram.h"
#include "gui/icc_profile.h"
#include "gui/keywords_tree.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/callbacks.h"
#include "gui/geometry.h"
#include "gui/gui_state.h"
#include "gui/masks_gui.h"
#include "gui/photometric_cc.h"
#include "gui/plot.h"
#include "gui/PSF_list.h"
#include "gui/registration_preview.h"
#include "gui/sequence_list.h"
#include "gui/siril_plot.h"
#include "gui/siril_preview.h"
#include "gui/script_console.h"
#include "gui/user_polygons.h"
#include "io/annotation_catalogues.h"
#include "io/sequence.h"
#include "livestacking/gui.h"
#include "gui/registration.h"
#include "gui/script_menu.h"
#include "gui/siril_actions.h"
#include "gui/siril-window.h"
#include "gui/stacking.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

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

/* ── Groups E, F: no-op placeholders (filled in by 6.1/6.2 impls below) ── */

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

/* ── Steps 6.1–6.4: sequence / image / command / console slots ──────────── */

/* Helper: free the sequence comboboxreglayer before closing a sequence. */
static void free_cbbt_layers(void) {
	GtkComboBoxText *cbbt = GTK_COMBO_BOX_TEXT(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "comboboxreglayer")));
	gtk_combo_box_text_remove_all(cbbt);
}

/* Helper: fill the AVI-export size fields from the current sequence. */
static void fillSeqAviExport(void) {
	char width[6], height[6];
	GtkEntry *heightEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviHeight")));
	GtkEntry *widthEntry  = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviWidth")));
	g_snprintf(width,  sizeof(width),  "%d", com.seq.rx);
	g_snprintf(height, sizeof(height), "%d", com.seq.ry);
	gtk_entry_set_text(widthEntry,  width);
	gtk_entry_set_text(heightEntry, height);
	if (com.seq.type == SEQ_SER && com.seq.ser_file) {
		GtkEntry *fpsEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviFps")));
		char fps[7];
		if (com.seq.ser_file->fps <= 0.0)
			g_snprintf(fps, sizeof(fps), "25.000");
		else
			g_snprintf(fps, sizeof(fps), "%2.3lf", com.seq.ser_file->fps);
		gtk_entry_set_text(fpsEntry, fps);
	}
}

/* GTK idle: update all GUI widgets after a sequence is opened. */
static gboolean set_seq_gui(gpointer user_data) {
	sequence *seq = (sequence *) user_data;
	init_layers_hi_and_lo_values(MIPSLOHI);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();
	int layer = set_layers_for_registration();
	update_seqlist(layer);
	fill_sequence_list(seq, max(layer, 0), FALSE);
	set_output_filename_to_sequence_name();
	sliders_mode_set_state(gui.sliders);
	initialize_display_mode();
	update_zoom_label();
	reset_plot();
	reset_3stars();
	set_display_mode();
	display_filename();
	gui_function(set_precision_switch, NULL);
	adjust_refimage(seq->current);
	update_prepro_interface(seq->type == SEQ_REGULAR ||
	                        seq->type == SEQ_FITSEQ ||
	                        seq->type == SEQ_SER);
	update_reg_interface(FALSE);
	update_stack_interface(FALSE);
	adjust_reginfo();
	update_gfit_histogram_if_needed();
	adjust_sellabel();
	fillSeqAviExport();
	update_MenuItem(NULL);
	set_GUI_CAMERA();
	gui_function(close_tab, NULL);
	gui_function(init_right_tab, NULL);
	notify_gfit_data_modified();
	gui_iface.redraw_image(REMAP_ALL);
	drawPlot();
	return FALSE;
}

/* GTK idle: clean up GUI state when a sequence is closed. */
static gboolean close_sequence_idle(gpointer data) {
	free_cbbt_layers();
	clear_sequence_list();
	clear_stars_list(TRUE);
	reset_3stars();
	clear_previews();
	free_reference_image();
	update_stack_interface(TRUE);
	adjust_sellabel();
	update_seqlist(-1);
	initialize_cut_struct(&gui.cut);
	if (!data) {
		GtkComboBox *seqcombo = GTK_COMBO_BOX(GTK_WIDGET(
			gtk_builder_get_object(gui.builder, "sequence_list_combobox")));
		gtk_combo_box_set_active(seqcombo, -1);
	}
	return FALSE;
}

/* GTK idle: populate the sequence selector with a single entry. */
static gboolean populate_seqcombo(gpointer user_data) {
	const gchar *realname = (const gchar *) user_data;
	control_window_switch_to_tab(IMAGE_SEQ);
	GtkComboBoxText *combo = GTK_COMBO_BOX_TEXT(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "sequence_list_combobox")));
	gtk_combo_box_text_remove_all(combo);
	gchar *rname = g_path_get_basename(realname);
	gtk_combo_box_text_append(combo, 0, rname);
	g_signal_handlers_block_by_func(GTK_COMBO_BOX(combo),
		on_seqproc_entry_changed, NULL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(combo),
		on_seqproc_entry_changed, NULL);
	g_free(rname);
	return FALSE;
}

static void impl_on_sequence_opened(void) {
	if (g_main_context_is_owner(g_main_context_default()))
		set_seq_gui(&com.seq);
	else
		execute_idle_and_wait_for_it(set_seq_gui, &com.seq);
}

static void impl_on_sequence_closed(gboolean loading_next) {
	if (g_main_context_is_owner(g_main_context_default()))
		close_sequence_idle(GINT_TO_POINTER(loading_next));
	else
		execute_idle_and_wait_for_it(close_sequence_idle,
		                             GINT_TO_POINTER(loading_next));
}

static void impl_populate_seq_combo(const char *realname) {
	gui_function(populate_seqcombo, (gpointer)realname);
}

/* GTK idle: free_image_data_gui — moved from io/single_image.c */
static gboolean free_image_data_gui(gpointer p) {
	disable_iso12646_conditions(TRUE, FALSE, FALSE);
	clear_user_polygons();
	delete_selected_area();
	reset_plot();
	siril_close_preview_dialogs();
	display_filename();
	update_zoom_label();
	update_display_fwhm();
	adjust_sellabel();
	gui_function(update_MenuItem, NULL);
	reset_3stars();
	gui_function(close_tab, NULL);

	GtkComboBox *binning = GTK_COMBO_BOX(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "combobinning")));
	GtkEntry *focal_entry  = GTK_ENTRY(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "focal_entry")));
	GtkEntry *pitchX_entry = GTK_ENTRY(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "pitchX_entry")));
	GtkEntry *pitchY_entry = GTK_ENTRY(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "pitchY_entry")));
	g_signal_handlers_block_by_func(focal_entry,  on_focal_entry_changed,  NULL);
	g_signal_handlers_block_by_func(pitchX_entry, on_pitchX_entry_changed, NULL);
	g_signal_handlers_block_by_func(pitchY_entry, on_pitchY_entry_changed, NULL);
	g_signal_handlers_block_by_func(binning,       on_combobinning_changed, NULL);

	clear_stars_list(TRUE);
	clear_backup();
	clear_sampling_setting_box();
	sample_mutex_lock();
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	sample_mutex_unlock();
	cleanup_annotation_catalogues(TRUE);
	reset_display_offset();
	reset_menu_toggle_button();
	reset_zoom_default();
	free(gui.qphot);
	gui.qphot = NULL;
	gui.show_wcs_disto = FALSE;
	clear_sensor_tilt();

	g_signal_handlers_unblock_by_func(focal_entry,  on_focal_entry_changed,  NULL);
	g_signal_handlers_unblock_by_func(pitchX_entry, on_pitchX_entry_changed, NULL);
	g_signal_handlers_unblock_by_func(pitchY_entry, on_pitchY_entry_changed, NULL);
	g_signal_handlers_unblock_by_func(binning,       on_combobinning_changed, NULL);
	siril_debug_print("free_image_data_idle() complete\n");

	for (int vport = 0; vport < MAXVPORT; vport++) {
		struct image_view *view = &gui.view[vport];
		if (view->buf) { free(view->buf); view->buf = NULL; }
		if (view->full_surface) {
			cairo_surface_destroy(view->full_surface);
			view->full_surface = NULL;
		}
		view->full_surface_stride = 0;
		view->full_surface_height = 0;
		if (view->disp_surface) {
			cairo_surface_destroy(view->disp_surface);
			view->disp_surface = NULL;
		}
		view->view_width = -1;
		view->view_height = -1;
	}
	clear_previews();
	free_reference_image();
	siril_debug_print("free_image_data_gui() complete\n");
	return FALSE;
}

static void impl_on_image_closed(void) {
	if (com.script || com.python_command)
		execute_idle_and_wait_for_it(free_image_data_gui, NULL);
	else if (!g_main_context_is_owner(g_main_context_default()))
		siril_add_idle(free_image_data_gui, NULL);
	else
		free_image_data_gui(NULL);
}

static void impl_on_image_loaded(void) {
	if (g_main_context_is_owner(g_main_context_default()))
		open_single_image_from_gfit(NULL);
	else
		execute_idle_and_wait_for_it(open_single_image_from_gfit, NULL);
}

/* ── Clear-log-buffer idle (moved from command.c) ────────────────────────── */
static gboolean clear_log_buffer_idle(gpointer user_data) {
	GtkTextView *text = GTK_TEXT_VIEW(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "output")));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
	return FALSE;
}

/* ── 6.1–6.4 slot implementations ────────────────────────────────────────── */

static void impl_set_gui_cwd(void) {
	gui_function(set_GUI_CWD, NULL);
}

static void impl_refresh_keywords_dialog(void) {
	gui_function(refresh_keywords_dialog, NULL);
}

static void impl_set_wcs_overlay(gboolean show) {
	gui.show_wcs_disto = show;
}

static void impl_launch_clipboard_survey(void) {
	gui_function(launch_clipboard_survey, NULL);
}

static void impl_clear_log_buffer(void) {
	gui_function(clear_log_buffer_idle, NULL);
}

static void impl_update_spin_cpu(void) {
	gui_function(update_spinCPU, GINT_TO_POINTER(0));
}

static void impl_new_selection_zone(void) {
	gui_function(new_selection_zone, NULL);
}

static int impl_get_active_vport(void) {
	return select_vport(gui.cvport);
}

static void impl_console_set_status(const char *msg, int line) {
	console_log_status(msg, line);
}

static void impl_console_clear_status(void) {
	gui_function((GSourceFunc)console_clear_status_bar, NULL);
}

static void impl_end_script_gui(void) {
	gui_function(end_script, NULL);
}

static gboolean impl_get_star_follow_state(void) {
	GtkToggleButton *follow = GTK_TOGGLE_BUTTON(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "followStarCheckButton")));
	return gtk_toggle_button_get_active(follow);
}

static void impl_show_command_help(void) {
	/* show_command_help_popup is in command_line_processor.c; route via gui_function
	 * by calling it with the "command" entry widget as the argument. */
	GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "command")));
	gui_function((GSourceFunc)show_command_help_popup, entry);
}

/* ── Steps 5.6–5.10: additional slots ───────────────────────────────────── */

/* Group C addition */
static void impl_data_dialog(SirilMessageType type, const char *title,
                             const char *text, const char *data) {
	siril_data_dialog(to_gtk_msg_type(type), (char *)title, (char *)text,
	                  (gchar *)data);
}

/* Group D addition */
static void impl_redraw_previews(void) {
	gui_function(redraw_previews, NULL);
}

/* Group E addition */
static void impl_open_single_image_from_gfit(void) {
	gui_function(open_single_image_from_gfit, NULL);
}

/* Group G additions */
static void impl_update_mem_usage(guint64 used_bytes) {
	set_GUI_MEM(used_bytes, "labelmem");
}

static void impl_update_disk_space(gint64 space_bytes, const char *label_id) {
	set_GUI_DiskSpace(space_bytes, label_id);
}

static void impl_update_mask_enable(gboolean state) {
	gui_function(set_mask_active_idle, GINT_TO_POINTER(state));
}

static void impl_set_display_range(int lo, int hi) {
	gui.lo = lo;
	gui.hi = hi;
	set_cutoff_sliders_values();
}

static void impl_check_gaia_status(void) {
	check_gaia_archive_status();
}

static void impl_trigger_gaia_check(void) {
	gaia_check(NULL);
}

/* ── Steps 5.11–5.15: additional slots ───────────────────────────────────── */

static void impl_remap_all_vports(void) {
	remap_all();
}

static void impl_quit_application(void) {
	gtk_main_quit();
}

static void impl_refresh_script_menu(void) {
	gui_mutex_lock();
	execute_idle_and_wait_for_it(refresh_script_menu_idle, NULL);
	gui_mutex_unlock();
}

static void impl_clear_backup(void) {
	clear_backup();
}

static void impl_livestacking_setup_gui(gboolean has_dark, gboolean has_flat,
                                        int reg_type) {
	gui.rendering_mode = STF_DISPLAY;
	set_display_mode();
	force_unlinked_channels();
	GtkWidget *toolbar = GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkToolMainBar"));
	if (gtk_widget_is_visible(toolbar))
		show_hide_toolbox();
	livestacking_display_config(has_dark, has_flat, (transformation_type)reg_type);
}

static void impl_livestacking_teardown_gui(void) {
	GtkWidget *toolbar = GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkToolMainBar"));
	if (!gtk_widget_is_visible(toolbar))
		show_hide_toolbox();
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

static gboolean impl_roi_operation_supports(void) {
	return gui.roi.operation_supports_roi;
}

static gpointer impl_get_roi_fit(void) {
	return (gpointer)&gui.roi.fit;
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
	if (com.gui_icc.proofing_transform) {
		cmsDeleteTransform(com.gui_icc.proofing_transform);
		com.gui_icc.proofing_transform = NULL;
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
	clear_psf_list_display();
}

static int impl_get_reg_layer(void) {
	return get_registration_layer_from_GUI(&com.seq);
}

/* ── SC: Annotation display ──────────────────────────────────────────────── */

static void impl_activate_annotation_display(void) {
	GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(
		gtk_builder_get_object(gui.builder, "annotate_button"));
	refresh_annotation_visibility();
	if (!gtk_toggle_tool_button_get_active(button)) {
		gtk_toggle_tool_button_set_active(button, TRUE);
	} else {
		refresh_found_objects();
		gui_iface.redraw_image(REDRAW_OVERLAY);
	}
}


/* ── U: ICC status icon ──────────────────────────────────────────────────── */

/* cm_struct and cm_worker moved here from core/icc_profile.c */
struct cm_struct {
	fits *fit;
	gboolean active;
};

static gboolean cm_worker(gpointer user_data) {
	struct cm_struct *data = (struct cm_struct *) user_data;
	fits *fit = data->fit;
	gboolean active = data->active;
	gchar *buffer = NULL, *monitor = NULL, *proof = NULL;
	gchar *name = g_build_filename("/org/siril/ui/", "pixmaps",
		active ? "color_management.svg" : "color_management_off.svg", NULL);
	gchar *tooltip = NULL;
	if (active) {
		if (fit->icc_profile) {
			buffer = siril_color_profile_get_description(fit->icc_profile);
			monitor = siril_color_profile_get_description(com.gui_icc.monitor);
		}
		if (com.gui_icc.soft_proof)
			proof = siril_color_profile_get_description(com.gui_icc.soft_proof);
		else
			proof = g_strdup(monitor);
		tooltip = g_strdup_printf(_("Image is color managed\nImage profile: %s\nMonitor profile: %s\nSoft proofing profile: %s"), buffer, monitor, proof);
		if (!tooltip)
			tooltip = g_strdup(_("Image is color managed\n\nLeft click: Color management dialog\nRight click: toggle ISO12646 color assessment mode"));
	} else {
		tooltip = g_strdup(_("Image is not color managed\n\nLeft click: Color management dialog\nRight click: toggle ISO12646 color assessment mode"));
	}
	GtkWidget *image = GTK_WIDGET(gtk_builder_get_object(gui.builder, "color_managed_icon"));
	GtkWidget *button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "icc_main_window_button"));
	gtk_image_set_from_resource((GtkImage *) image, name);
	gtk_widget_set_tooltip_text(button, tooltip);
	g_free(name);
	g_free(buffer);
	g_free(monitor);
	g_free(proof);
	g_free(tooltip);
	return FALSE;
}

static void impl_update_icc_status_icon(gpointer p, gboolean active) {
	fits *fit = (fits *) p;
	struct cm_struct data = { fit, active };
	if (g_main_context_is_owner(g_main_context_default())) {
		cm_worker(&data);
	} else {
		execute_idle_and_wait_for_it(cm_worker, &data);
	}
}

static gboolean impl_get_gamut_check_active(void) {
	GtkToggleButton *checkgamut = GTK_TOGGLE_BUTTON(
		gtk_builder_get_object(gui.builder, "checkgamut"));
	return gtk_toggle_button_get_active(checkgamut);
}

/* ── V: Registration panel status ────────────────────────────────────────── */

static void impl_update_registration_status(const gchar *msg) {
	registration_update_label(msg);
}

/* ── Group S: Pixel-math status ──────────────────────────────────────────── */

static void impl_update_pixel_math_status(int ret) {
	GtkLabel *label = GTK_LABEL(GTK_WIDGET(
		gtk_builder_get_object(gui.builder, "pixel_math_status")));
	if (!label) return;
	gtk_label_set_text(label, ret == 0 ? "" : _("Syntax error"));
}

/* ── Group R: Python bridge UI ───────────────────────────────────────────── */

static void impl_update_single_image_display(void) {
	if (g_main_context_is_owner(g_main_context_default()))
		update_single_image_from_gfit(NULL);
	else
		execute_idle_and_wait_for_it(update_single_image_from_gfit, NULL);
}

static void impl_seq_redisplay_frame(int index) {
	int idx = index;
	execute_idle_and_wait_for_it(seq_load_image_in_thread, &idx);
}

static void impl_set_poly_drawing(guint32 color, gboolean fill) {
	gui.poly_fill = fill;
	gui.poly_ink  = uint32_to_gdk_rgba(color);
	init_draw_poly();
}

static void impl_set_rendering_mode(int mode) {
	gui.rendering_mode = (display_mode)mode;
	execute_idle_and_wait_for_it(set_display_mode_idle, NULL);
}

static void impl_set_channels_linked(gboolean state) {
	gui.unlink_channels = !state;
	execute_idle_and_wait_for_it(chain_channels_idle_callback, GINT_TO_POINTER(state));
}

static void impl_set_sliders_mode(int mode) {
	sliders_mode sliders = (sliders_mode)mode;
	execute_idle_and_wait_for_it(sliders_mode_set_state_idle, &sliders);
}

static void impl_set_cutoff_values(int lo, int hi) {
	gui.lo = lo;
	gui.hi = hi;
	execute_idle_and_wait_for_it(set_cutoff_sliders_values_idle, NULL);
}

static void impl_update_zoom_label(void) {
	execute_idle_and_wait_for_it(update_zoom_label_idle, NULL);
}

static double impl_get_zoom_value(void) {
	return get_zoom_val();
}

static int impl_activate_action(const char *name, gboolean appmap) {
	return (int)queue_activate_action_if_enabled(name, appmap);
}

static void impl_reset_display_offset(void) {
	reset_display_offset();
}

/* ── SF: Python IPC display state ───────────────────────────────────────── */

static int impl_get_channel_for_vport(void) {
	return match_drawing_area_widget(gui.view[select_vport(gui.cvport)].drawarea, FALSE);
}

static int impl_get_rendering_mode(void) {
	return (int)gui.rendering_mode;
}

static gboolean impl_get_channels_linked(void) {
	return !gui.unlink_channels;
}

static void impl_get_display_offset(double *x, double *y) {
	if (x) *x = gui.display_offset.x;
	if (y) *y = gui.display_offset.y;
}

static void impl_set_display_offset(double x, double y) {
	gui.display_offset.x = x;
	gui.display_offset.y = y;
}

static void impl_set_zoom_value(double zoom) {
	gui.zoom_value = zoom;
}

static GSList *impl_get_user_polygons(void) {
	return gui.user_polygons;
}

static void impl_add_user_polygon_to_list(gpointer polygon) {
	gui.user_polygons = g_slist_append(gui.user_polygons, polygon);
}

/* ── SD: Display range state ─────────────────────────────────────────────── */

static void impl_get_display_lo_hi(int *lo, int *hi) {
	if (lo) *lo = (int)gui.lo;
	if (hi) *hi = (int)gui.hi;
}

static int impl_get_sliders_mode(void) {
	return (int)gui.sliders;
}

static void impl_update_display_range_after_load(int sliders, int lo, int hi) {
	gui.sliders = (sliders_mode)sliders;
	gui.lo = (WORD)lo;
	gui.hi = (WORD)hi;
}

/* ── Registration ────────────────────────────────────────────────────────── */

static void impl_apply_display_icc_compensation(gpointer p);

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
	gui_iface.on_sequence_opened      = impl_on_sequence_opened;
	gui_iface.on_sequence_closed      = impl_on_sequence_closed;
	gui_iface.populate_seq_combo      = impl_populate_seq_combo;
	gui_iface.on_image_loaded         = impl_on_image_loaded;
	gui_iface.on_image_closed         = impl_on_image_closed;
	gui_iface.set_gui_cwd             = impl_set_gui_cwd;
	gui_iface.refresh_keywords_dialog = impl_refresh_keywords_dialog;
	gui_iface.set_wcs_overlay         = impl_set_wcs_overlay;
	gui_iface.launch_clipboard_survey = impl_launch_clipboard_survey;
	gui_iface.clear_log_buffer        = impl_clear_log_buffer;
	gui_iface.update_spin_cpu         = impl_update_spin_cpu;
	gui_iface.new_selection_zone      = impl_new_selection_zone;
	gui_iface.get_active_vport        = impl_get_active_vport;
	gui_iface.console_set_status      = impl_console_set_status;
	gui_iface.console_clear_status    = impl_console_clear_status;
	gui_iface.end_script_gui          = impl_end_script_gui;
	/* on_image_loaded and on_image_closed are wired above in the 6.1/6.2 block */
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
	gui_iface.roi_operation_supports      = impl_roi_operation_supports;
	gui_iface.get_roi_fit                 = impl_get_roi_fit;
	gui_iface.get_roi_selection           = impl_get_roi_selection;
	gui_iface.clear_roi                   = impl_clear_roi;
	gui_iface.restore_roi                 = impl_restore_roi;
	gui_iface.reset_display_transform     = impl_reset_display_transform;
	gui_iface.is_preview_active           = impl_is_preview_active;
	gui_iface.hide_preview                = impl_hide_preview;
	gui_iface.copy_gfit_to_backup         = impl_copy_gfit_to_backup;
	gui_iface.copy_gfit_icc_to_backup     = impl_copy_gfit_icc_to_backup;
	gui_iface.check_icc_identical_to_monitor = impl_check_icc_identical_to_monitor;
	gui_iface.set_source_information         = impl_set_source_information;
	gui_iface.apply_display_icc_compensation = impl_apply_display_icc_compensation;
	gui_iface.data_dialog                 = impl_data_dialog;
	gui_iface.redraw_previews             = impl_redraw_previews;
	gui_iface.open_single_image_from_gfit = impl_open_single_image_from_gfit;
	gui_iface.update_mem_usage            = impl_update_mem_usage;
	gui_iface.update_disk_space           = impl_update_disk_space;
	gui_iface.update_mask_enable          = impl_update_mask_enable;
	gui_iface.set_display_range           = impl_set_display_range;
	gui_iface.check_gaia_status           = impl_check_gaia_status;
	gui_iface.trigger_gaia_check          = impl_trigger_gaia_check;
	gui_iface.remap_all_vports            = impl_remap_all_vports;
	gui_iface.quit_application            = impl_quit_application;
	gui_iface.refresh_script_menu         = impl_refresh_script_menu;
	gui_iface.clear_backup                = impl_clear_backup;
	gui_iface.livestacking_setup_gui      = impl_livestacking_setup_gui;
	gui_iface.livestacking_teardown_gui   = impl_livestacking_teardown_gui;
	gui_iface.get_star_follow_state       = impl_get_star_follow_state;
	gui_iface.show_command_help           = impl_show_command_help;
	gui_iface.update_pixel_math_status    = impl_update_pixel_math_status;
	gui_iface.activate_annotation_display = impl_activate_annotation_display;
	gui_iface.update_icc_status_icon      = impl_update_icc_status_icon;
	gui_iface.get_gamut_check_active      = impl_get_gamut_check_active;
	gui_iface.update_registration_status  = impl_update_registration_status;
	gui_iface.update_single_image_display = impl_update_single_image_display;
	gui_iface.seq_redisplay_frame         = impl_seq_redisplay_frame;
	gui_iface.set_poly_drawing            = impl_set_poly_drawing;
	gui_iface.set_rendering_mode          = impl_set_rendering_mode;
	gui_iface.set_channels_linked         = impl_set_channels_linked;
	gui_iface.set_sliders_mode            = impl_set_sliders_mode;
	gui_iface.set_cutoff_values           = impl_set_cutoff_values;
	gui_iface.update_zoom_label           = impl_update_zoom_label;
	gui_iface.get_zoom_value              = impl_get_zoom_value;
	gui_iface.activate_action             = impl_activate_action;
	gui_iface.reset_display_offset        = impl_reset_display_offset;
	gui_iface.get_channel_for_vport           = impl_get_channel_for_vport;
	gui_iface.get_rendering_mode              = impl_get_rendering_mode;
	gui_iface.get_channels_linked             = impl_get_channels_linked;
	gui_iface.get_display_offset              = impl_get_display_offset;
	gui_iface.set_display_offset              = impl_set_display_offset;
	gui_iface.set_zoom_value                  = impl_set_zoom_value;
	gui_iface.get_user_polygons               = impl_get_user_polygons;
	gui_iface.add_user_polygon_to_list        = impl_add_user_polygon_to_list;
	gui_iface.get_display_lo_hi               = impl_get_display_lo_hi;
	gui_iface.get_sliders_mode                = impl_get_sliders_mode;
	gui_iface.update_display_range_after_load = impl_update_display_range_after_load;
}

static void impl_apply_display_icc_compensation(gpointer p) {
	fits *fit = (fits *)p;
	if (!fit || !com.gui_icc.monitor) return;
	int depth = fit->naxes[2];
	if (depth == 1) {
		fits_change_depth(fit, 3);
		if (fit->type == DATA_FLOAT) {
			memcpy(fit->fpdata[1], fit->fdata, fit->rx * fit->ry * sizeof(float));
			memcpy(fit->fpdata[2], fit->fdata, fit->rx * fit->ry * sizeof(float));
		} else {
			memcpy(fit->pdata[1], fit->data, fit->rx * fit->ry * sizeof(WORD));
			memcpy(fit->pdata[2], fit->data, fit->rx * fit->ry * sizeof(WORD));
		}
	}
	cmsHPROFILE temp = copyICCProfile(fit->icc_profile);
	cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = copyICCProfile(com.gui_icc.monitor);
	siril_colorspace_transform(fit, temp);
	cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = copyICCProfile(temp);
	cmsCloseProfile(temp);
	if (depth == 1) {
		size_t npixels = fit->rx * fit->ry;
		if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++)
				fit->fdata[i] = (fit->fpdata[0][i] + fit->fpdata[1][i] + fit->fpdata[2][i]) / 3.f;
		} else {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++)
				fit->data[i] = (fit->pdata[0][i] + fit->pdata[1][i] + fit->pdata[2][i]) / 3;
		}
		fits_change_depth(fit, 1);
	}
}
