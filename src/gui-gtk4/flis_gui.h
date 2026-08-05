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
 */

#ifndef FLIS_GUI_H
#define FLIS_GUI_H

#include <gtk/gtk.h>
#include <glib.h>

/*
 * flis_gui — FLIS layers panel (stage 4).
 *
 * Modeless top-level GtkWindow that owns the layer list, the per-layer
 * property panel, the layer-mask sub-frame, the toolbar, and the
 * context menu.  Built entirely in procedural C (no .ui file beyond
 * the small flis_layer_row.ui row template) because the panel morphs
 * heavily with selection state and the cost of expressing that morph
 * in a static GtkBuilder template exceeds the value of the
 * declaration.
 *
 * Public surface:
 *   • flis_gui_toggle_visible() — called by the win.show-layers
 *     action (Ctrl+L).  Allocates the panel on first show and
 *     thereafter shows/hides it.  Safe in any GTK4 build, no-op in
 *     headless contexts (the action isn't registered there).
 *   • flis_gui_update_from_idle() — called by image_format_flis
 *     primitives (via gui_iface.flis_gui_update) after a mutation
 *     to schedule a refresh on the main loop.  Safe from worker
 *     threads.  The first call after a non-FLIS → FLIS image open
 *     also populates the panel from scratch.
 */

void flis_gui_toggle_visible(void);
void flis_gui_update_from_idle(void);

/* Nondestructive History popover, anchored below the header-bar clock button
 * (nde_history_button).  _init builds and attaches it at startup
 * (initialize_all_GUI); _toggle_visible is the programmatic toggle behind
 * win.show-nde-history.  Main thread only. */
void flis_gui_history_popover_init(void);
void flis_gui_history_toggle_visible(void);

/* Dismiss the pinned history popover when no single image is loaded any
 * more; fed from update_MenuItem.  Main thread only. */
void flis_gui_history_notify_image_state(gboolean single_image_loaded);

/* Coalesced refresh of the History window only (NDE provenance mirror).
 * Routed from gui_iface.nde_history_changed; callable from any thread. */
void flis_gui_history_update_from_idle(void);

/* Coalesced refresh for per-motion drag ticks: immediate paint (GPU
 * compose path reads live values) + one idle-priority full pipeline
 * rebuild once the motion-event burst drains. */
void flis_drag_tick_refresh(void);

/* Idempotent show: present the panel (building it if needed) only when
 * a FLIS is currently loaded.  No-op for plain FITS / sequences / no
 * image.  Wired from single_image.c::open_single_image_from_gfit via
 * gui_iface.flis_gui_present_if_flis so opening a FLIS automatically
 * reveals the layers panel without stealing focus when the user is
 * working with non-FLIS images. */
void flis_gui_present_if_flis(void);

/* TRUE when the layers-panel selection currently sits on a GROUP row.
 * Reads the panel's cached selection state (plain ints), so it is safe
 * to call from worker threads via gui_iface.flis_group_is_selected. */
gboolean flis_panel_group_is_selected(void);

/* item_id of the selected group row, 0 when the selection is not a group. */
gint flis_panel_selected_group_id(void);

/* TRUE when the selected group plausibly composites to colour. */
gboolean flis_panel_selected_colour_group(void);

/* TRUE when more than one panel row is selected (GtkMultiSelection). */
gboolean flis_panel_multi_is_selected(void);

#endif /* FLIS_GUI_H */
