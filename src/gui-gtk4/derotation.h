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
gboolean derotation_get_disk(double *cx, double *cy, double *radius,
                             double *pa_deg, gboolean *mirrored);

#endif /* SRC_GUI_GTK4_DEROTATION_H_ */
