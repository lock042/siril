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

#ifndef SRC_GUI_CUT_H_
#define SRC_GUI_CUT_H_

gboolean reset_cut_gui_filedependent(gpointer user_data);
double get_conversion_factor(fits *fit);
void measure_line(fits* fit, point start, point finish, gboolean pref_as);
void initialize_cut_struct(cut_struct *arg);
gboolean cut_struct_is_valid(cut_struct *arg);
void free_cut_args(cut_struct *arg);
/* Launch a single-image cut profile via generic_image_worker (read-only path).
 * Selects the correct hook (profile / tri / cfa) from arg->tri / arg->cfa.
 * Ownership of arg transfers to the framework on success; caller frees on failure. */
gboolean launch_cut_single(cut_struct *arg);
/* Low-level workers used by the sequence path (cut_image_hook in apply_cut_to_sequence). */
gpointer cut_profile(gpointer p);
gpointer tri_cut(gpointer p);
gpointer cfa_cut(gpointer p);
void apply_cut_to_sequence(cut_struct *cut_args);
void update_spectro_labels();

#endif
