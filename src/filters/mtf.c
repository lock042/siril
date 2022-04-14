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
	int i, invertedChannels = 0;
	imstats *stat[3];

	int nb_channels = (int)fit->naxes[2];

	int retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
	for (i = 0; i < nb_channels; ++i) {
		if (stat[i]) {
			if (retval)
				free_stats(stat[i]);
			else if (stat[i]->median / stat[i]->normValue > 0.5)
				++invertedChannels;
		}
	}
	if (retval) {
		result->shadows = 0.0f;
		result->midtones = 0.2f;
		result->highlights = 1.0f;
		return -1;
	}

	if (invertedChannels < nb_channels) {
		for (i = 0; i < nb_channels; ++i) {
			float median, mad, normValue;

			normValue = (float)stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			c0 += median + SHADOWS_CLIPPING * mad;
			m += median;
		}
		c0 /= (float) nb_channels;
		if (c0 < 0.f) c0 = 0.f;
		float m2 = m / (float) nb_channels - c0;
		result->midtones = MTF(m2, TARGET_BACKGROUND, 0.f, 1.f);
		result->shadows = c0;
		result->highlights = 1.0f;

		siril_debug_print("autostretch: (%f, %f, %f)\n",
				result->shadows, result->midtones, result->highlights);
	} else {
		for (i = 0; i < nb_channels; ++i) {
			float median, mad, normValue;

			normValue = (float) stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			m += median;
			c1 += median - SHADOWS_CLIPPING * mad;
		}
		c1 /= (float) nb_channels;
		if (c1 > 1.f) c1 = 1.f;
		float m2 = c1 - m / (float) nb_channels;
		result->midtones = 1.f - MTF(m2, TARGET_BACKGROUND, 0.f, 1.f);
		result->shadows = 0.f;
		result->highlights = c1;

	}
	for (i = 0; i < nb_channels; ++i)
		free_stats(stat[i]);
	return 0;
}

void apply_unlinked_mtf_to_fits(fits *from, fits *to, struct mtf_params *params) {
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
#ifdef _OPENMP
		int threads = com.max_thread >= 3 ? 3 : com.max_thread;
#pragma omp parallel for num_threads(threads) schedule(static) if (threads> 1)
#endif
		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
			for (size_t i = 0; i < ndata; i++) {
				float pxl = (float)from->pdata[chan][i] / norm;
				float mtf = MTFp(pxl, params[chan]);
				to->pdata[chan][i] = round_to_WORD(mtf * norm);
			}
		}
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
		int threads = com.max_thread >= 3 ? 3 : com.max_thread;
#pragma omp parallel for num_threads(threads) schedule(static) if (threads> 1)
#endif
		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
			for (size_t i = 0; i < ndata; i++) {
				to->fpdata[chan][i] = MTFp(from->fpdata[chan][i], params[chan]);
			}
		}
	}
	else return;

	invalidate_stats_from_fit(to);
}


int find_unlinked_midtones_balance(fits *fit, float shadows_clipping, float target_bg, struct mtf_params *results) {
	int i, invertedChannels = 0;
	imstats *stat[3];

	int nb_channels = (int)fit->naxes[2];

	int retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
	for (i = 0; i < nb_channels; ++i) {
		if (retval) {
			results[i].shadows = 0.0f;
			results[i].midtones = 0.2f;
			results[i].highlights = 1.0f;
			if (stat[i])
				free_stats(stat[i]);
		}
		else if (stat[i]->median / stat[i]->normValue > 0.5)
				++invertedChannels;
	}
	if (retval)
		return -1;

	if (invertedChannels < nb_channels) {
		for (i = 0; i < nb_channels; ++i) {
			float normValue = (float)stat[i]->normValue;
			float median = (float) stat[i]->median / normValue;
			float mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			float c0 = median + shadows_clipping * mad;
			if (c0 < 0.f) c0 = 0.f;
			float m2 = median - c0;
			results[i].midtones = MTF(m2, target_bg, 0.f, 1.f);
			results[i].shadows = c0;
			results[i].highlights = 1.0;
			siril_debug_print("autostretch for channel %d: (%f, %f, %f)\n",
					i, results[i].shadows, results[i].midtones, results[i].highlights);
		}
	} else {
		for (i = 0; i < nb_channels; ++i) {
			float normValue = (float) stat[i]->normValue;
			float median = (float) stat[i]->median / normValue;
			float mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			float c1 = median - shadows_clipping * mad;
			if (c1 > 1.f) c1 = 1.f;
			float m2 = c1 - median;
			results[i].midtones = 1.f - MTF(m2, target_bg, 0.f, 1.f);
			results[i].shadows = 0.f;
			results[i].highlights = c1;
			siril_debug_print("autostretch for channel %d: (%f, %f, %f)\n",
					i, results[i].shadows, results[i].midtones, results[i].highlights);
		}

	}
	for (i = 0; i < nb_channels; ++i)
		free_stats(stat[i]);
	return 0;
}

