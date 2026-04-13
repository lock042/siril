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
#ifndef SRC_GUI_LINEAR_MATCH_H_
#define SRC_GUI_LINEAR_MATCH_H_

#include "core/siril.h"
#include "core/processing.h"

struct linear_match_data {
	void (*destroy_fn)(void *); /* Must be first member */
	fits ref;                   /* Reference image (owned by this struct) */
	double low;
	double high;
};

struct linear_match_data *new_linear_match_data(fits *ref_fit, double low, double high);
void free_linear_match_data(void *p);
int linear_match_image_hook(struct generic_img_args *args, fits *fit, int threads);
void apply_linear_to_fits(fits *fit, double *a, double *b);

#endif /* SRC_GUI_LINEAR_MATCH_H_ */
