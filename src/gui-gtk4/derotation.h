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

#ifndef SRC_GUI_GTK4_DEROTATION_H_
#define SRC_GUI_GTK4_DEROTATION_H_

#include <gtk/gtk.h>
#include <glib.h>

/* Launch handler for the registration tab's "Derotation…" button. */
void on_seqmpp_derotation_button_clicked(GtkButton *button, gpointer user_data);

/* Overlay support: true while the window is open, and the current disk fit
 * (full-frame pixels) so the image canvas can draw a live outline. */
gboolean derotation_is_open(void);

/* Called when a sequence is loaded/displayed: if the derotation tool is open,
 * refresh its per-sequence metadata (reference epoch, frame rate) and re-fit the
 * disk for the newly loaded sequence. No-op when the tool is closed. */
void derotation_sequence_changed(void);

gboolean derotation_get_disk(double *cx, double *cy, double *radius, double *rpol,
                             double *pa_deg, gboolean *mirrored);
/* Selected target body (0 = Jupiter, 1 = Saturn, 2 = Mars), or -1 if unknown.
 * Lets the overlay draw body-specific fitting guides. */
int derotation_get_body(void);

/* Interactive fit (mouse). Coordinates are display/image pixels (y-down).
 * hit_test returns a drag mode (0 none, 1 move centre, 2 equatorial radius,
 * 3 polar radius); set/get_drag hold the active drag; drag_to applies it. */
int  derotation_hit_test(double dx, double dy);
void derotation_set_drag(int mode);
int  derotation_get_drag(void);
void derotation_drag_to(double dx, double dy);

/* Rotation grab-handle centre in display/image pixels (y-down). Normally beyond
 * the north-pole tip; when that falls outside the image it is pulled to the
 * midpoint of the north axis (inside the disc) and *inside is set TRUE. Used by
 * both the overlay and the hit-test so they stay in sync. Returns FALSE if no
 * fit is active. */
gboolean derotation_rotate_handle(double *hx, double *hy, gboolean *inside);

#endif /* SRC_GUI_GTK4_DEROTATION_H_ */
