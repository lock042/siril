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
#ifndef SRC_PIXELMATH_PIXEL_MATH_RUNNER_H_
#define SRC_PIXELMATH_PIXEL_MATH_RUNNER_H_

#include "core/siril.h"

struct pixel_math_data {
	fits *fit;
	int nb_rows;
	gchar **varname;
	int ret;
	gboolean from_ui;
	gboolean has_gfit;
	gboolean do_sum;
	gchar *expression1;
	gchar *expression2;
	gchar *expression3;
	gboolean single_rgb;
	gboolean rescale;
	float min, max;
};

int load_pm_var(const gchar *var, int index, int *w, int *h, int *c);
void free_pm_var(int nb);
gpointer apply_pixel_math_operation(gpointer p);

#endif /* SRC_PIXELMATH_PIXEL_MATH_RUNNER_H_ */
