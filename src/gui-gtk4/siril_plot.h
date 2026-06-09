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

#ifndef SRC_GUI_SIRIL_PLOT_H_
#define SRC_GUI_SIRIL_PLOT_H_

#include "core/siril.h"
#include "io/siril_plot.h"
gboolean create_new_siril_plot_window(gpointer p);
gchar* build_save_filename(gchar *prepend, gchar *ext, gboolean forsequence, gboolean add_time_stamp);
gboolean save_siril_plot_to_clipboard(siril_plot_data *spl_data, int width, int height);

#endif /* SRC_GUI_PLOT_H_ */
