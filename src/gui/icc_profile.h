/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#ifndef SRC_GUI_ICC_PROFILE_H_
#define SRC_GUI_ICC_PROFILE_H_
#include <stdint.h>
#include <lcms2.h>

void set_icc_description_in_TIFF();
void initialize_icc_preferences_widgets();
void set_source_information();
gboolean on_icc_main_window_button_clicked(GtkWidget *btn, GdkEventButton *event, gpointer userdata);
void enable_iso12646_conditions();
void disable_iso12646_conditions(gboolean revert_zoom, gboolean revert_panel, gboolean revert_rendering_mode);

#endif
