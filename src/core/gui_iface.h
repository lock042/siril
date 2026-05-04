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

#include <glib.h>   /* gboolean, gchar — GLib only, no GTK */

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

	/* D – Image display ---------------------------------------------------- */
	/* Synchronous redraw (call only from the GTK main thread). */
	void     (*redraw_image)(SirilRedrawType remap);
	/* Queue a redraw from any thread; returns immediately. */
	void     (*redraw_image_async)(SirilRedrawType remap);
	/* Queue a redraw from any thread and block until it completes. */
	void     (*redraw_image_sync)(SirilRedrawType remap);
	/* Clear any active selection rectangle from the image display. */
	void     (*delete_selection)(void);

	/* E – Sequence / image state notifications ----------------------------- */
	/* Called after a sequence is fully opened and ready for use. */
	void     (*on_sequence_opened)(void);
	/* Called after a single image is loaded into gfit and displayed. */
	void     (*on_image_loaded)(void);
	/* Called when the current single image is closed/unloaded. */
	void     (*on_image_closed)(void);
	/* Called after stacking completes successfully; gfit holds the result. */
	void     (*on_stack_complete)(void);

	/* F – Panel / tab switching -------------------------------------------- */
	/* Show or hide a named UI panel (sidebar tab, floating window, etc.). */
	void     (*show_panel)(const char *panel_name, gboolean visible);

	/* G – Misc GUI state --------------------------------------------------- */
	/* Refresh the main-window status bar (disk space, memory, zoom, …). */
	void     (*update_status_bar)(void);
	/* Sync menu/toolbar enable-state to current application state. */
	void     (*update_menu_state)(void);

	/* H – Geometry / ROI / Mask state ------------------------------------- */
	/* Called before a geometry-altering operation; clears ROI if active. */
	void     (*on_geometry_changed)(void);
	/* Called when a mask is removed or its visibility state changes. */
	void     (*on_mask_state_changed)(void);
	/* Called after a crop completes; clears stars, selection, display offset. */
	void     (*on_crop_complete)(void);

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
