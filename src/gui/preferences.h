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
#ifndef SRC_GUI_PREFERENCES_H_
#define SRC_GUI_PREFERENCES_H_

void notify_script_update();
void update_libraw_and_debayer_interface();
void update_photometry_interface();
void initialize_path_directory(const gchar *path);
void set_libraw_settings_menu_available(gboolean activate);
void initialize_compression_param();
void initialize_default_preferences();
void update_preferences_from_model();
gchar *get_swap_dir();

#endif /* SRC_GUI_PREFERENCES_H_ */
