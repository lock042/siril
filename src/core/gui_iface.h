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
 * Abstract GUI interface (vtable) for the Siril processing core.
 *
 * This header is GTK-free and may be included from any translation unit.
 * It defines the residual "bridge" between the processing core and the
 * GUI toolkit, covering the interface groups described in the
 * GUI-separation plan (groups A–D initially; E–G added incrementally in
 * phase 2.4 onwards as later phases demand them).
 *
 * All slots default to no-op stub implementations defined in
 * src/core/gui_iface_stubs.c.  The GUI build overrides them by calling
 * siril_register_gui_iface() from main.c after GTK initialisation.
 *
 * Usage in processing code:
 *   #include "core/gui_iface.h"
 *   gui_iface.set_progress(0.5, _("Half done"));
 *   if (gui_iface.confirm_dialog(...)) { ... }
 */

#ifndef SRC_CORE_GUI_IFACE_H_
#define SRC_CORE_GUI_IFACE_H_

#include <glib.h>          /* gboolean, gchar — GLib only, no GTK */
#include "core/settings.h" /* rectangle — used by Group H ROI slots */

/* Forward declaration for Group K — avoids including siril.h or PSF.h */
struct fwhm_struct;
typedef struct fwhm_struct psf_star;

#ifdef __cplusplus
extern "C" {
#endif

/* ── Group A: Progress bar constants ─────────────────────────────────────── */
/* Moved here from gui/progress_and_log.h so processing code can use them
 * without pulling in GTK headers. */
#define PROGRESS_NONE     -2.0  /* do not update the progress bar fraction */
#define PROGRESS_PULSATE  -1.0  /* pulsate (indeterminate) mode */
#define PROGRESS_RESET     0.0  /* reset progress bar to zero */
#define PROGRESS_DONE      1.0  /* fill progress bar to 100 % */
#define PROGRESS_TEXT_RESET ""  /* reset the progress bar text to "Ready." */

/* ── Group D: Redraw type ────────────────────────────────────────────────── */
/*
 * Defines what must be regenerated before the next paint.  Moved here from
 * gui/image_display.h so non-GUI code can name these values without
 * including GTK headers.  gui/image_display.h aliases this as remap_type
 * for backward compatibility.
 */
typedef enum {
	REDRAW_OVERLAY, /* only annotation/overlay layers changed */
	REDRAW_IMAGE,   /* image pixel values changed; re-render but no remap */
	REMAP_ALL,      /* image data or display LUT changed; remap then render */
} SirilRedrawType;

/* ── ActionResult — returned by activate_action() ───────────────────────── */
typedef enum {
    ACTION_SUCCESS        = 0,
    ACTION_NOT_FOUND      = 1,
    ACTION_DISABLED       = 2,
    ACTION_WINDOW_MISSING = 3,
    ACTION_NULL_DATA      = 4,
} ActionResult;

/* ── Group C: Message type ───────────────────────────────────────────────── */
/* Toolkit-neutral replacement for GtkMessageType used in the dialog slot. */
typedef enum {
	SIRIL_MSG_INFO,
	SIRIL_MSG_WARNING,
	SIRIL_MSG_ERROR,
	SIRIL_MSG_QUESTION,
} SirilMessageType;

/* ── GUI interface vtable ────────────────────────────────────────────────── */
/*
 * All slots are initialised to no-op stubs (gui_iface_stubs.c).
 * The GUI build replaces them with GTK implementations via
 * siril_register_gui_iface() (gui/gui_iface_impl.c).
 *
 * Convention for callers: just call the slot directly — the stubs are
 * always valid function pointers, so no NULL guard is needed.
 */
typedef struct {
	/* A – Progress / status ------------------------------------------------ */
	/* fraction: PROGRESS_* constant or [0.0, 1.0]; msg may be NULL */
	void     (*set_progress)(double fraction, const char *msg);
	/* busy=TRUE: show "waiting" cursor; busy=FALSE: restore normal cursor */
	void     (*set_busy)(gboolean busy);

	/* B – Logging ---------------------------------------------------------- */
	/* Append msg to the GUI script-log widget in the given colour.
	 * color is a Pango colour name string (e.g. "red"), or NULL for default. */
	void     (*log_message)(const char *msg, const char *color);

	/* C – Modal dialogs ---------------------------------------------------- */
	void     (*message_dialog)(SirilMessageType type, const char *title,
	                           const char *text);
	/* Returns TRUE if the user clicked the accept button. */
	gboolean (*confirm_dialog)(const char *title, const char *msg,
	                           const char *button_accept);
	/* id is the GtkBuilder identifier of the dialog widget. */
	void     (*open_dialog)(const char *id);
	void     (*close_dialog)(const char *id);
	/* Returns TRUE if an image-processing dialog is currently open. */
	gboolean (*is_dialog_open)(void);
	/* Rich-data dialog (e.g. update notes); data is optional supplementary
	 * text shown below the main message (may be NULL). */
	void     (*data_dialog)(SirilMessageType type, const char *title,
	                        const char *text, const char *data);

	/* D – Image display ---------------------------------------------------- */
	/* Synchronous redraw (call only from the GTK main thread). */
	void     (*redraw_image)(SirilRedrawType remap);
	/* Queue a redraw from any thread; returns immediately. */
	void     (*redraw_image_async)(SirilRedrawType remap);
	/* Queue a redraw from any thread and block until it completes. */
	void     (*redraw_image_sync)(SirilRedrawType remap);
	/* Clear any active selection rectangle from the image display. */
	void     (*delete_selection)(void);
	/* Queue a redraw of the mask overlay only (from any thread). */
	void     (*queue_redraw_mask)(void);
	/* Refresh all preview windows (registration / filter previews). */
	void     (*redraw_previews)(void);
	/* Remap all display viewports (recalculate display LUT/buffers). */
	void     (*remap_all_vports)(void);

	/* E – Sequence / image state notifications ----------------------------- */
	/* Called after a sequence is fully opened and ready for use. */
	void     (*on_sequence_opened)(void);
	/* Called when the current sequence is closed/unloaded.
	 * loading_next: TRUE when immediately loading another sequence (skip
	 * deselecting the combo box). */
	void     (*on_sequence_closed)(gboolean loading_next);
	/* Called after a single image is loaded into gfit and displayed. */
	void     (*on_image_loaded)(void);
	/* Called when the current single image is closed/unloaded. */
	void     (*on_image_closed)(void);
	/* Called after stacking completes successfully; gfit holds the result. */
	void     (*on_stack_complete)(void);
	/* Update the sequence list display; seqname is the basename to select,
	 * or NULL to refresh without changing selection. */
	void     (*update_sequences_list)(const char *seqname);
	/* Open/display gfit as a single image in the main window. */
	void     (*open_single_image_from_gfit)(void);
	/* Populate the sequence-selector combo with realname as the only entry. */
	void     (*populate_seq_combo)(const char *realname);

	/* F – Panel / tab switching -------------------------------------------- */
	/* Show or hide a named UI panel (sidebar tab, floating window, etc.). */
	void     (*show_panel)(const char *panel_name, gboolean visible);
	/* Enable or disable script-related UI widgets (run/stop buttons, etc.). */
	void     (*script_widgets_enable)(gboolean enable);
	/* FFMS2/AVI: enable or disable the sequence browser widget. */
	void     (*set_seq_browser_active)(gboolean active);

	/* G – Misc GUI state --------------------------------------------------- */
	/* Refresh the main-window status bar (disk space, memory, zoom, …). */
	void     (*update_status_bar)(void);
	/* Sync menu/toolbar enable-state to current application state. */
	void     (*update_menu_state)(void);
	/* Suppress (TRUE) or restore (FALSE) drawarea redraws during processing.
	 * When suppressing, also disables the display-mode menu button. */
	void     (*set_suppress_redraws)(gboolean suppress);
	/* Repopulate the ROI display from current gfit data (call while holding
	 * the gfit read lock). */
	void     (*populate_roi)(void);

	/* D additions – Histogram / image modification state ------------------ */
	/* Mark the histogram as stale (call while holding at least a read lock). */
	void     (*invalidate_histogram)(void);
	/* Recompute histogram if stale (call while holding at least a read lock). */
	void     (*update_histogram)(void);
	/* Queue an idle mask redraw (may also trigger a full image redraw). */
	void     (*redraw_mask_idle)(void);

	/* F additions – Application lifecycle -------------------------------- */
	/* Quit the application's main event loop. */
	void     (*quit_application)(void);
	/* Update the CWD display label in the main window. */
	void     (*set_gui_cwd)(void);

	/* G additions – Channel / precision / command state ----------------- */
	/* Called when the channel count of gfit may have changed. */
	void     (*on_channel_count_changed)(void);
	/* Called when the data type (float/ushort) of gfit has changed. */
	void     (*on_precision_changed)(void);
	/* Refresh the script menu on the GTK main thread (serialised). */
	void     (*refresh_script_menu)(void);
	/* Refresh the keywords-tree dialog from current gfit keywords. */
	void     (*refresh_keywords_dialog)(void);
	/* Show/hide the WCS distortion overlay on the image canvas. */
	void     (*set_wcs_overlay)(gboolean show);
	/* Launch the photometric/SPCC catalogue survey dialog. */
	void     (*launch_clipboard_survey)(void);
	/* Clear the main log text buffer. */
	void     (*clear_log_buffer)(void);
	/* Update the CPU-thread-count spin button. */
	void     (*update_spin_cpu)(void);
	/* Sync the selection-rectangle UI state from com.selection. */
	void     (*new_selection_zone)(void);
	/* Return the currently active display viewport index. */
	int      (*get_active_vport)(void);
	/* Return TRUE if the star-following toggle is active. */
	gboolean (*get_star_follow_state)(void);
	/* Show the command help popup for the current console command. */
	void     (*show_command_help)(void);
	/* Update the memory-usage label (used_bytes = current RSS). */
	void     (*update_mem_usage)(guint64 used_bytes);
	/* Update a disk-space label (space_bytes = free bytes, label_id = widget name). */
	void     (*update_disk_space)(gint64 space_bytes, const char *label_id);
	/* Sync the mask-enable toggle button to state (blocks its signal handler). */
	void     (*update_mask_enable)(gboolean state);
	/* Set the display lo/hi cutoff values and refresh slider widgets. */
	void     (*set_display_range)(int lo, int hi);
	/* Trigger a Gaia archive connectivity check. */
	void     (*check_gaia_status)(void);
	/* Initiate an on-demand Gaia catalogue check. */
	void     (*trigger_gaia_check)(void);

	/* H – Geometry / ROI / Mask state ------------------------------------- */
	/* Called before a geometry-altering operation; clears ROI if active. */
	void     (*on_geometry_changed)(void);
	/* Called when a mask is removed or its visibility state changes. */
	void     (*on_mask_state_changed)(void);
	/* Called after a crop completes; clears stars, selection, display offset. */
	void     (*on_crop_complete)(void);
	/* Returns TRUE if an ROI selection is currently active. */
	gboolean (*roi_is_active)(void);
	/* Copy the current ROI selection rectangle into *rect. */
	void     (*get_roi_selection)(rectangle *rect);
	/* Clear any active ROI selection. */
	void     (*clear_roi)(void);
	/* Set *rect as the new ROI selection and notify the GUI. */
	void     (*restore_roi)(const rectangle *rect);
	/* Reset (delete) the ICC proofing display transform under its lock. */
	void     (*reset_display_transform)(void);

	/* I – Statistics ------------------------------------------------------- */
	/* Called to populate and open the statistics dialog for the current image. */
	void     (*on_stats_ready)(void);

	/* J – Photometry ------------------------------------------------------- */
	/* Called after photometry data is updated; refreshes the plot panel. */
	void     (*on_photometry_changed)(void);
	/* Open a new siril-plot window for the given siril_plot_data pointer. */
	void     (*show_siril_plot)(gpointer spl_data);

	/* K – Star list -------------------------------------------------------- */
	/* Update the star list display and optionally the PSF list panel. */
	void     (*update_star_list)(psf_star **stars, gboolean update_psf_list,
	                             gboolean wait);
	/* Clear the star list display (equivalent to clear_stars_list(FALSE)). */
	void     (*clear_star_list)(void);

	/* L – Registration state ----------------------------------------------- */
	/* Returns the active registration layer from the GUI combo box, or -1 if
	 * no GUI is available or the selected layer has no registration data. */
	int      (*get_reg_layer)(void);

	/* M – Thread utilities ------------------------------------------------- */
	/* Run func(data) on the GTK main thread and block until it completes.
	 * The stub calls func(data) directly on the calling thread.
	 * MUST NOT be called from the GTK main thread (deadlock). */
	void     (*execute_idle_sync)(GSourceFunc func, gpointer data);
	/* Open a dialog to select which child process to kill; returns the chosen
	 * GPid, or 0 if unavailable (single child or headless). */
	GPid     (*select_child_process)(GSList *children);

	/* N – Preview / backup state ------------------------------------------ */
	/* Returns TRUE if a preview effect is currently applied to gfit. */
	gboolean (*is_preview_active)(void);
	/* Cancel the active preview and restore gfit from the backup buffer. */
	void     (*hide_preview)(void);
	/* Copy gfit pixel data into the preview backup buffer. */
	void     (*copy_gfit_to_backup)(void);
	/* Copy the gfit ICC profile into the preview backup (after ICC ops). */
	void     (*copy_gfit_icc_to_backup)(void);
	/* Discard the preview backup buffer entirely. */
	void     (*clear_backup)(void);

	/* O – ICC information callbacks --------------------------------------- */
	/* Check whether gfit's ICC profile matches the monitor profile and update
	 * the color-management icon accordingly. */
	void     (*check_icc_identical_to_monitor)(void);
	/* Update the image source/profile information labels in the main window. */
	void     (*set_source_information)(void);
	/* Apply the monitor ICC profile to fit so that a stretched preview looks
	 * correct on the current display.  Callers cast fits* to gpointer.
	 * No-op in headless mode. */
	void     (*apply_display_icc_compensation)(gpointer fit);

	/* Q – Script console -------------------------------------------------- */
	/* Append msg to the script-log status bar; line is the script line number. */
	void     (*console_set_status)(const char *msg, int line);
	/* Clear the script-log status bar. */
	void     (*console_clear_status)(void);
	/* Run all post-script GUI cleanup on the GTK main thread. */
	void     (*end_script_gui)(void);

	/* P – Livestacking GUI lifecycle --------------------------------------- */
	/* Called when livestacking starts; switches display mode and hides toolbar. */
	void     (*livestacking_setup_gui)(gboolean has_dark, gboolean has_flat,
	                                   int reg_type);
	/* Called when livestacking stops; restores toolbar visibility. */
	void     (*livestacking_teardown_gui)(void);

	/* S – Pixel-math status ------------------------------------------------ */
	/* Called when a pixel-math operation completes; updates the status label.
	 * ret is the pixel_math_data.ret value (0 = success, non-zero = error). */
	void (*update_pixel_math_status)(int ret);

	/* SC – Annotation display ----------------------------------------------- */
	/* Ensure the annotation overlay is visible and up-to-date after a
	 * catalogue search.  Activates the annotate toggle if needed. */
	void (*activate_annotation_display)(void);

	/* T – Window handle ----------------------------------------------------- */
	/* Return the main application window (GtkWindow*) as void* so the header
	 * stays GTK-free.  Callers in GUI code cast to GtkWindow*. */
	void *(*get_main_window)(void);

	/* U – ICC status --------------------------------------------------------- */
	/* Update the ICC color-management status icon and tooltip in the toolbar.
	 * fit is cast to fits* in the implementation; use gpointer to keep the
	 * header GTK/siril.h-free. */
	void     (*update_icc_status_icon)(gpointer fit, gboolean active);
	/* Return TRUE if the gamut check toggle in the ICC dialog is active. */
	gboolean (*get_gamut_check_active)(void);

	/* V – Registration panel status ----------------------------------------- */
	/* Set the info label text in the 3-star registration panel. */
	void (*update_registration_status)(const gchar *msg);

	/* R – Python bridge UI ------------------------------------------------- */
	/* Refresh the single-image display (wraps update_single_image_from_gfit).
	 * Threading-aware: dispatches to main thread if called from worker. */
	void     (*update_single_image_display)(void);
	/* Reload and display sequence frame at index (wraps seq_load_image_in_thread). */
	void     (*seq_redisplay_frame)(int index);
	/* Set the polygon drawing ink colour (packed RGBA) and fill flag, then
	 * initialise the draw-polygon state machine. */
	void     (*set_poly_drawing)(guint32 color, gboolean fill);
	/* Set the display rendering mode and refresh the display. */
	void     (*set_rendering_mode)(int mode);
	/* Set whether autostretch channels are linked/unlinked and refresh. */
	void     (*set_channels_linked)(gboolean state);
	/* Set the slider (display) mode and refresh the display.
	 * mode is cast to sliders_mode in the implementation. */
	void     (*set_sliders_mode)(int mode);
	/* Set the lo/hi cutoff slider values and refresh the display. */
	void     (*set_cutoff_values)(int lo, int hi);
	/* Update the zoom label widget (wraps update_zoom_label_idle). */
	void     (*update_zoom_label)(void);
	/* Return the current effective zoom value (handles ZOOM_FIT). */
	double   (*get_zoom_value)(void);
	/* Activate a named GAction on the application or window map.
	 * Returns an ActionResult int (ACTION_SUCCESS=0 etc.). */
	int      (*activate_action)(const char *name, gboolean appmap);
	/* Reset the image display pan offset to (0, 0). */
	void     (*reset_display_offset)(void);
} SirilGuiInterface;

/* The single global GUI interface instance.  Defined in gui_iface_stubs.c. */
extern SirilGuiInterface gui_iface;

/*
 * Called once from main.c (GUI build only) after the GtkBuilder UI is
 * loaded.  Replaces stub implementations with GTK-backed ones.
 * Defined in gui/gui_iface_impl.c (compiled only in GUI builds).
 */
void siril_register_gui_iface(void);

#ifdef __cplusplus
}
#endif

#endif /* SRC_CORE_GUI_IFACE_H_ */
