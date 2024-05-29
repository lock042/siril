/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#ifndef SRC_GUI_ANNOTATE_H_
#define SRC_GUI_ANNOTATE_H_

typedef struct {
	float limit_mag;
	gboolean photometric;
	super_bool display_tag;
	super_bool display_log;
	siril_cat_index cat;
	gchar *obscode;
	gboolean default_obscode_used;
	int trixel;
	gchar *outfilename;
} conesearch_params;

typedef struct {
	gboolean clear;
	char *file;
	char *name;
	SirilWorldCS *coords;
	gboolean display_log;
	gboolean display_tag;
} show_params;

int execute_conesearch(conesearch_params *params);
int execute_show_command(show_params *params);

#endif /* SRC_GUI_ANNOTATE_H_ */
