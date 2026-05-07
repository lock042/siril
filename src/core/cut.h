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

/*
 * Non-GUI lifecycle helpers for cut_struct (intensity profile / spectrogram).
 * The GTK dialog callbacks remain in gui/cut.c.
 */

#ifndef SRC_CORE_CUT_H_
#define SRC_CORE_CUT_H_

#include "core/siril.h"

/* Initialise all fields of a cut_struct to safe defaults.
 * Frees any existing heap strings inside *arg before zeroing them. */
void initialize_cut_struct(cut_struct *arg);

/* Free all heap-allocated strings inside *arg, then free the struct itself.
 * Do NOT call with &gui.cut; use initialize_cut_struct for that instead. */
void free_cut_args(cut_struct *arg);

/* Validate the configuration in *arg; returns FALSE and logs a message on
 * the first constraint violation, TRUE if everything is in order. */
gboolean cut_struct_is_valid(cut_struct *arg);

#endif /* SRC_CORE_CUT_H_ */
