/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_GUI_SIRIL_WINDOW_H_
#define SRC_GUI_SIRIL_WINDOW_H_

void siril_window_enable_image_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_wcs_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_wcs_disto_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_autostretch_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_rgb_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_any_rgb_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_rgb_wcs_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_any_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_any_mono_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_single_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_none_proc_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_if_selection_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_if_selection_rgb_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_enable_if_selection_sequence_actions(GtkApplicationWindow *window, gboolean enable);
void siril_window_map_actions(GtkApplicationWindow *window);

#endif /* SRC_GUI_SIRIL_WINDOW_H_ */
