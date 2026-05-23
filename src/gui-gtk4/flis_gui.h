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

#endif /* FLIS_GUI_H */
