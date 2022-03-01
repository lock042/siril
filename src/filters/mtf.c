/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "mtf.h"
#include "core/proto.h"
#include "algos/statistics.h"


void apply_linked_mtf_to_fits(fits *from, fits *to, struct mtf_params params) {

	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1] * from->naxes[2];
	g_assert(from->type == to->type);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			float pxl = (float)from->data[i] / norm;
			float mtf = MTFp(pxl, params);
			to->data[i] = round_to_WORD(mtf * norm);
		}
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			to->fdata[i] = MTFp(from->fdata[i], params);
		}
	}
	else return;

	invalidate_stats_from_fit(to);
}

float MTF(float x, float m, float lo, float hi) {
	if (x <= lo)
		return 0.f;
	if (x >= hi)
		return 1.f;

	float xp = (x - lo) / (hi - lo);

	return ((m - 1.f) * xp) / (((2.f * m - 1.f) * xp) - m);
}

float MTFp(float x, struct mtf_params params) {
	return MTF(x, params.midtones, params.shadows, params.highlights);
}

int find_linked_midtones_balance(fits *fit, struct mtf_params *result) {
	float c0 = 0.0, c1 = 0.0;
	float m = 0.0;
	int i, n, invertedChannels = 0;
	imstats *stat[3];

	n = fit->naxes[2];

	int retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_MAD, MULTI_THREADED, -1, stat);
	for (i = 0; i < n; ++i) {
		if (retval && stat[i])
			free_stats(stat[i]);
		else if (stat[i]->median / stat[i]->normValue > 0.5)
			++invertedChannels;
	}
	if (retval)
		return -1;

	if (invertedChannels < n) {
		for (i = 0; i < n; ++i) {
			float median, mad, normValue;

			normValue = (float)stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			c0 += median + shadowsClipping * mad;
			m += median;
		}
		c0 /= (float) n;
		if (c0 < 0.f) c0 = 0.f;
		float m2 = m / (float) n - c0;
		result->midtones = MTF(m2, targetBackground, 0.f, 1.f);
		result->shadows = c0;
		result->highlights = 1.0;
	} else {
		for (i = 0; i < n; ++i) {
			float median, mad, normValue;

			normValue = (float) stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			m += median;
			c1 += median - shadowsClipping * mad;
		}
		c1 /= (float) n;
		if (c1 > 1.f) c1 = 1.f;
		float m2 = c1 - m / (float) n;
		result->midtones = 1.f - MTF(m2, targetBackground, 0.f, 1.f);
		result->shadows = 0.f;
		result->highlights = c1;

	}
	for (i = 0; i < n; ++i)
		free_stats(stat[i]);
	return 0;
}

