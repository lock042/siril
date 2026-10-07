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
#ifndef SRC_GUI_DIALOG_PREVIEW_H_
#include <gtk/gtk.h>
#define SRC_GUI_DIALOG_PREVIEW_H_

/* The GtkFileChooser-based preview pipeline was removed: GTK4 has no
 * file-chooser preview API, and previews are now rendered by the custom
 * file browser (siril_file_browser_default_preview in file_browser.c, which
 * reuses extract_thumbnail_from_ser/_fits/_avi directly). */

#endif /* SRC_GUI_DIALOG_PREVIEW_H_ */
