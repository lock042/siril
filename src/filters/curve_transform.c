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

#include <glib.h>
#include "curve_transform.h"
#include "core/proto.h"
#include <math.h>

void apply_curve(fits *from, fits *to, struct curve_params params, gboolean multithreaded) {
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	g_assert(from->type == to->type);
	const size_t layersize = from->naxes[0] * from->naxes[1];

	cubic_spline_data cspline_data;
	// Linear data
	double slopes[MAX_POINTS - 1];

	if (params.algorithm == CUBIC_SPLINE)
		cubic_spline_fit(params.points, &cspline_data);
	else if (params.algorithm == LINEAR)
		linear_fit(params.points, slopes);


	if (from->type == DATA_USHORT) {
		float norm = (float) get_normalized_value(from);
		float inv_norm = 1.f / norm;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t i = 0; i < from->naxes[2]; i++) {
			if (params.do_channel[i]) {
				for (size_t j = 0; j < layersize; j++) {
					float pixel_value = from->pdata[i][j] * inv_norm;
					if (params.algorithm == LINEAR)
						to->pdata[i][j] = roundf_to_WORD(linear_interpolate(pixel_value, params.points, slopes) * norm);
					else if (params.algorithm == CUBIC_SPLINE)
						to->pdata[i][j] = roundf_to_WORD(cubic_spline_interpolate(pixel_value, &cspline_data) * norm);
				}
			} else
				memcpy(to->pdata[i], from->pdata[i], layersize * sizeof(WORD));
		}
	} else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t i = 0; i < from->naxes[2]; i++) {
			if (params.do_channel[i]) {
				for (size_t j = 0; j < layersize; j++) {
					float pixel_value = from->fpdata[i][j];
					if (params.algorithm == LINEAR)
						to->fpdata[i][j] = linear_interpolate(pixel_value, params.points, slopes);
					else if (params.algorithm == CUBIC_SPLINE)
						to->fpdata[i][j] = cubic_spline_interpolate(pixel_value, &cspline_data);
				}
			} else
				memcpy(to->fpdata[i], from->fpdata[i], layersize * sizeof(float));
		}
	}
	return;
}

void linear_fit(GList *points, double *slopes) {
	g_assert(g_list_length(points) >= 2);

	// Precalculate the slope between each pair of points
	GList *current = points;
	for (int i = 0; i < g_list_length(points) - 1; i++) {
		point *point1 = (point *) current->data;
		point *point2 = (point *) current->next->data;
		slopes[i] = (point2->y - point1->y) / (point2->x - point1->x);
		current = current->next;
	}
}

float linear_interpolate(float x, GList *points, double *slopes) {
	if (x > ((point *) g_list_last(points)->data)->x)
		return ((point *) g_list_last(points)->data)->y;
	else if (x < ((point *) points->data)->x)
		return ((point *) points->data)->y;

	x = fmax(0, fmin(1, x));

	// Find the point that is closest to x with its x value < x
	GList *point1 = points;
	int point_index = 0;
	while (point1->next->next != NULL && ((point *) point1->next->data)->x < x) {
		point_index++;
		point1 = point1->next;
	}

	float x1 = ((point *) point1->data)->x;
	float y1 = ((point *) point1->data)->y;
	return y1 + slopes[point_index] * (x - x1);
}

void cubic_spline_fit(GList *points, cubic_spline_data *cspline_data) {
	g_assert(g_list_length(points) >= 2);

	cspline_data->n = g_list_length(points);
	double *h = malloc((cspline_data->n - 1) * sizeof(double));
	double *alpha = malloc((cspline_data->n - 1) * sizeof(double));
	double *l = malloc(cspline_data->n * sizeof(double));
	double *mu = malloc(cspline_data->n * sizeof(double));
	double *z = malloc(cspline_data->n * sizeof(double));

	l[0] = 1;
	mu[0] = 0;
	z[0] = 0;

	GList *current = points;
	for (int i = 0; i < cspline_data->n; i++) {
		point *p = (point *) current->data;
		cspline_data->x_values[i] = p->x;
		cspline_data->y_values[i] = p->y;
		current = current->next;
	}

	for (int i = 0; i < cspline_data->n - 1; i++)
		h[i] = cspline_data->x_values[i + 1] - cspline_data->x_values[i];

	for (int i = 1; i < cspline_data->n - 1; i++)
		alpha[i] = (3 / h[i]) * (cspline_data->y_values[i + 1] - cspline_data->y_values[i]) -
				   (3 / h[i - 1]) * (cspline_data->y_values[i] - cspline_data->y_values[i - 1]);

	for (int i = 1; i < cspline_data->n - 1; i++) {
		l[i] = 2 * (cspline_data->x_values[i + 1] - cspline_data->x_values[i - 1]) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[cspline_data->n - 1] = 1;
	z[cspline_data->n - 1] = 0;
	cspline_data->c[cspline_data->n - 1] = 0;

	double factor = 1.0 / 3;
	for (int j = cspline_data->n - 2; j >= 0; j--) {
		cspline_data->c[j] = z[j] - mu[j] * cspline_data->c[j + 1];
		cspline_data->b[j] = (cspline_data->y_values[j + 1] - cspline_data->y_values[j]) / h[j] -
							 h[j] * (cspline_data->c[j + 1] + 2 * cspline_data->c[j]) * factor;
		cspline_data->d[j] = (cspline_data->c[j + 1] - cspline_data->c[j]) / (3 * h[j]);
	}

	free(h);
	free(alpha);
	free(l);
	free(mu);
	free(z);
}

float cubic_spline_interpolate(float x, cubic_spline_data *cspline_data) {
	if (x > cspline_data->x_values[cspline_data->n - 1])
		return cspline_data->y_values[cspline_data->n - 1];
	else if (x < cspline_data->x_values[0])
		return cspline_data->y_values[0];


	x = fmax(0, fmin(1, x));
	int i = 0;
	while (i < cspline_data->n - 1 && x > cspline_data->x_values[i + 1])
		i++;

	double diff = x - cspline_data->x_values[i];
	double diff_sq = diff * diff;
	double interpolated_value =
			cspline_data->y_values[i] + cspline_data->b[i] * diff + cspline_data->c[i] * diff_sq +
			cspline_data->d[i] * diff * diff_sq;
	return fmax(0, fmin(1, interpolated_value));
}
