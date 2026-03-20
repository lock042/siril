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

#endif /* FLIS_GUI_H */
