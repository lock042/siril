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
 * GUI-side helpers for geometry operations (crop, rotate, mirror, etc.).
 * Compiled only in GUI builds.  Pure processing stays in algos/geometry.c.
 */

#include <gtk/gtk.h>
#include "core/siril.h"
#include "gui-gtk4/PSF_list.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/geometry.h"

/* Idle callback that refreshes GUI state after a crop.
 * Scheduled via gui_iface.on_crop_complete() from algos/geometry.c and
 * core/command.c so that the GUI is updated after every crop path. */
gboolean crop_gui_updates(gpointer user) {
	clear_stars_list(TRUE);
	delete_selected_area();
	reset_display_offset();
	update_zoom_label();
	return FALSE;
}
