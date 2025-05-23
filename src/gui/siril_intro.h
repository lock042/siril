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
#ifndef SRC_GUI_SIRIL_INTRO_H_
#define SRC_GUI_SIRIL_INTRO_H_

#include <gtk/gtk.h>

typedef enum {
	SIRIL_INTRO_GTK_POPOVER,
	SIRIL_INTRO_WORKAROUND_POPOVER // this workaround is used when popover must be displayed on a child. popover is cropped, here we use a GtkDialog.
} type_of_popover;


typedef struct {
	GtkWidget *popover;
	GtkWidget *widget;
} SirilUIIntro;

typedef struct {
	const gchar *widget;
	const gchar *tip;
	const gint delay;
	type_of_popover type;
} SirilTipIntro;

void start_intro_script();

#endif /* SRC_GUI_SIRIL_INTRO_H_ */
