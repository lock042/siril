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

#include <string.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/fitting.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "linear_match.h"

static void apply_linear_to_fits_ushort(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(fit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->pdata[channel][i] = round_to_WORD(fit->pdata[channel][i] * a[channel] + b[channel] * USHRT_MAX_DOUBLE);
		}
	}
}

static void apply_linear_to_fits_float(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(fit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->fpdata[channel][i] = fit->fpdata[channel][i] * a[channel] + b[channel];
		}
	}
}

void apply_linear_to_fits(fits *fit, double *a, double *b) {
	if (fit->type == DATA_USHORT) {
		apply_linear_to_fits_ushort(fit, a, b);
	} else if (fit->type == DATA_FLOAT) {
		apply_linear_to_fits_float(fit, a, b);
	}
}

struct linear_match_data *new_linear_match_data(fits *ref_fit, double low, double high) {
	struct linear_match_data *data = calloc(1, sizeof(struct linear_match_data));
	if (!data)
		return NULL;
	data->destroy_fn = free_linear_match_data;
	data->low = low;
	data->high = high;
	/* Caller must have already loaded ref_fit; we take ownership of its data */
	data->ref = *ref_fit;
	/* Zero out the caller's copy so clearfits() won't free the same buffers */
	memset(ref_fit, 0, sizeof(fits));
	return data;
}

void free_linear_match_data(void *p) {
	struct linear_match_data *data = (struct linear_match_data *)p;
	if (!data)
		return;
	clearfits(&data->ref);
	free(data);
}

int linear_match_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct linear_match_data *data = (struct linear_match_data *)args->user;
	double a[3] = { 0.0 }, b[3] = { 0.0 };
	if (find_linear_coeff(fit, &data->ref, data->low, data->high, a, b, NULL))
		return 1;
	apply_linear_to_fits(fit, a, b);
	return 0;
}

gchar *linear_match_log_hook(gpointer p, log_hook_detail detail) {
	struct linear_match_data *data = (struct linear_match_data *)p;
	return g_strdup_printf(_("Linear match (low: %.3f, high: %.3f)"), data->low, data->high);
}
