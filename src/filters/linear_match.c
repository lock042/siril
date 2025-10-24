/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"

static void apply_linear_to_fits_ushort(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(&gfit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->pdata[channel][i] = round_to_WORD(fit->pdata[channel][i] * a[channel] + b[channel] * USHRT_MAX_DOUBLE);
		}
	}
}

static void apply_linear_to_fits_float(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(&gfit);
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
