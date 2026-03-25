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

#include <glib.h>
#include <gtk/gtk.h>

/**
 * flis_gui_init:
 *
 * One-time initialisation called once after the main UI has been loaded at
 * startup.  Connects signals that cannot be expressed in Glade (specifically
 * the GtkAdjustment value-changed signal shared by the opacity scale and
 * spinbutton).  Call this from wherever Siril performs its post-load GUI
 * setup.
 */
void flis_gui_init(void);

/**
 * flis_gui_open:
 *
 * Creates (first call) and shows the FLIS layers panel window.
 * Safe to call when no image is loaded or when a plain FITS is loaded;
 * the panel shows an appropriate placeholder in those states.
 * Idempotent: calling it when the window is already visible just raises it.
 */
void flis_gui_open(void);

/**
 * flis_gui_update:
 *
 * Rebuilds the layer list and updates all property widgets to reflect the
 * current state of com.uniq.  Must be called from the GTK main thread.
 * Call this after any operation that changes the layer stack or the
 * properties of any layer.
 */
void flis_gui_update(void);

/**
 * flis_gui_update_from_idle:
 *
 * Thread-safe variant of flis_gui_update().  Schedules the update to run
 * on the GTK main thread via siril_add_idle().  Safe to call from a
 * processing thread.
 */
void flis_gui_update_from_idle(void);

/**
 * flis_remove_selected_lmask:
 *
 * Removes the layer mask from the currently selected layer, prompting for
 * confirmation.  Saves an undo state before removal.  No-op if no layer is
 * selected or it has no layer mask.  Must be called from the GTK main thread.
 */
void flis_remove_selected_lmask(void);

/**
 * flis_get_selected_layer_id:
 *
 * Returns the item_id of the currently selected layer in the FLIS panel,
 * or 0 if no layer is selected or the current image is not a FLIS.
 */
gint flis_get_selected_layer_id(void);

/**
 * flis_populate_layer_combo:
 * @combo: a GtkComboBoxText to populate
 * @filter_by_canvas_size: if TRUE, only layers whose pixel dimensions match
 *   gfit (the composited canvas) are included.  Pass FALSE when the mask
 *   source dimensions are not yet known (e.g. file-mode dialogs before a
 *   file has been selected).
 *
 * Fills @combo with the layer names from the current FLIS, in layer order.
 * No-op if the current image is not a FLIS.
 */
void flis_populate_layer_combo(GtkComboBoxText *combo, gboolean filter_by_canvas_size);

/**
 * flis_combo_select_active_layer:
 * @combo: a GtkComboBoxText previously filled by flis_populate_layer_combo()
 *
 * Selects the entry whose text matches the currently active layer's name.
 * Falls back to the first entry if the active layer is not present in the
 * combo (e.g. filtered out by canvas-size).  No-op if @combo is NULL.
 */
void flis_combo_select_active_layer(GtkComboBoxText *combo);

#endif /* FLIS_GUI_H */
