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

#include <glib.h>
#include "mtf.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"

void apply_linked_mtf_to_fits(fits *from, fits *to, struct mtf_params params, gboolean multithreaded) {

	const gboolean do_channel[3] = { params.do_red, params.do_green, params.do_blue };

	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);
#ifdef _OPENMP
	int threads = min(com.max_thread, 2); // not worth using many threads here
#endif
	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
		// Set up a LUT
		WORD *lut = siril_malloc((USHRT_MAX + 1) * sizeof(WORD));

		// This is only a small loop: 8 threads seems to be about as many as is worthwhile
		// because of the thread startup cost
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(threads) schedule(static) if (threads > 1)
#endif
		for (int i = 0 ; i <= USHRT_MAX ; i++) { // Fill LUT
			lut[i] = roundf_to_WORD(USHRT_MAX_SINGLE * MTFp(i * invnorm, params));
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				WORD *tolayer = to->pdata[j], *fromlayer = from->pdata[j];
				for (size_t i = 0; i < layersize; i++) {
					tolayer[i] = lut[fromlayer[i]];
				}
			} else
				memcpy(to->pdata[j], from->pdata[j], layersize * sizeof(WORD));
		}
		siril_free(lut);
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				float *tolayer = to->fpdata[j], *fromlayer = from->fpdata[j];
				for (size_t i = 0; i < layersize; i++) {
					tolayer[i] = MTFp(fromlayer[i], params);
				}
			} else
				memcpy(to->fpdata[j], from->fpdata[j], layersize * sizeof(float));
		}
	}
	else return;

	invalidate_stats_from_fit(to);
}

// In the general case this cannot return the inverse of MTF() because MTF() may clip.
// However when autostretch parameters are used the amount of clipping is negligible
// and this pseudoinverse gives a good approximation.
float MTF_pseudoinverse(float y, struct mtf_params params) {
	return ((((params.shadows + params.highlights) * params.midtones
			- params.shadows) * y - params.shadows * params.midtones
			+ params.shadows)
			/ ((2 * params.midtones - 1) * y - params.midtones + 1));
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

void apply_linked_pseudoinverse_mtf_to_fits(fits *from, fits *to, struct mtf_params params, gboolean multithreaded) {
// This is for use in reversing the pre-stretch applied to linear images for starnet++ input.
// It does not support selected channels.
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);

	const gboolean do_channel[3] = { params.do_red, params.do_green, params.do_blue };
#ifdef _OPENMP
	int threads = min(com.max_thread, 2); // not worth using many threads here
#endif
	siril_log_message(_("Applying inverse MTF with values %f, %f, %f\n"),
			params.shadows, params.midtones, params.highlights);
	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
		// Set up a LUT
		WORD *lut = siril_malloc((USHRT_MAX + 1) * sizeof(WORD));
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(multithreaded)
#endif
		for (int i = 0 ; i <= USHRT_MAX ; i++) { // Fill LUT
			lut[i] = roundf_to_WORD(USHRT_MAX_SINGLE * MTF_pseudoinverse(i * invnorm, params));
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				for (size_t i = 0; i < layersize; i++) {
					to->pdata[j][i] = lut[from->pdata[j][i]];
				}
			} else
				memcpy(to->pdata[j], from->pdata[j], layersize * sizeof(WORD));
		}
		siril_free(lut);
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				for (size_t i = 0; i < layersize; i++) {
					to->fpdata[j][i] = MTF_pseudoinverse(from->fpdata[j][i], params);
				}
			} else
				memcpy(to->fpdata[j], from->fpdata[j], layersize * sizeof(float));
		}

	}
	else return;

	invalidate_stats_from_fit(to);
}

int find_linked_midtones_balance(fits *fit, float shadows_clipping, float target_bg, struct mtf_params *result) {
	float c0 = 0.0, c1 = 0.0;
	float m = 0.0;
	int i, invertedChannels = 0;
	imstats *stat[3];

	int nb_channels = (int)fit->naxes[2];
	int retval = compute_all_channels_statistics_single_image(fit,
			STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
	for (i = 0; i < nb_channels; i++) {
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
			float normValue = (float)stat[i]->normValue;
			float median = (float) stat[i]->median / normValue;
			float mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			c0 += median + shadows_clipping * mad;
			m += median;
		}
		c0 /= (float) nb_channels;
		if (c0 < 0.f) c0 = 0.f;
		float m2 = m / (float) nb_channels - c0;
		result->shadows = c0;
		result->midtones = MTF(m2, target_bg, 0.f, 1.f);
		result->highlights = 1.0f;

		siril_debug_print("autostretch: (%f, %f, %f)\n",
				result->shadows, result->midtones, result->highlights);
	} else {
		for (i = 0; i < nb_channels; ++i) {
			float normValue = (float)stat[i]->normValue;
			float median = (float) stat[i]->median / normValue;
			float mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			m += median;
			c1 += median - shadows_clipping * mad;
		}
		c1 /= (float) nb_channels;
		if (c1 > 1.f) c1 = 1.f;
		float m2 = c1 - m / (float) nb_channels;
		result->midtones = 1.f - MTF(m2, target_bg, 0.f, 1.f);
		result->shadows = 0.f;
		result->highlights = c1;

	}
	for (i = 0; i < nb_channels; i++)
		free_stats(stat[i]);
	return 0;
}

int find_linked_midtones_balance_default(fits *fit, struct mtf_params *result) {
	return find_linked_midtones_balance(fit,
			AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, result);
}

void apply_unlinked_mtf_to_fits(fits *from, fits *to, struct mtf_params *params) {
	int threads = com.max_thread;
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
		// Set up a LUT
		WORD *lut = siril_malloc((USHRT_MAX + 1) * sizeof(WORD));

		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(threads) schedule(static) if (threads > 1)
#endif
			for (int i = 0 ; i <= USHRT_MAX ; i++) { // Fill LUT
				lut[i] = roundf_to_WORD(USHRT_MAX_SINGLE * MTFp(i * invnorm, params[chan]));
			}
			siril_log_message(_("Applying MTF to channel %d with values %f, %f, %f\n"), chan,
					params[chan].shadows, params[chan].midtones, params[chan].highlights);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < ndata; i++) {
				to->pdata[chan][i] = lut[from->pdata[chan][i]];
			}
		}
		siril_free(lut);
	}
	else if (from->type == DATA_FLOAT) {
		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
			siril_log_message(_("Applying MTF to channel %d with values %f, %f, %f\n"), chan,
					params[chan].shadows, params[chan].midtones, params[chan].highlights);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (threads > 1)
#endif
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

	int retval = compute_all_channels_statistics_single_image(fit,
			STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
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
			siril_debug_print("autostretch for channel %d: (%f, %f, %f)\n", i,
					results[i].shadows, results[i].midtones, results[i].highlights);
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
			siril_debug_print("autostretch for channel %d: (%f, %f, %f)\n", i,
					results[i].shadows, results[i].midtones, results[i].highlights);
		}

	}
	for (i = 0; i < nb_channels; ++i)
		free_stats(stat[i]);
	return 0;
}

int find_unlinked_midtones_balance_default(fits *fit, struct mtf_params *results) {
	return find_unlinked_midtones_balance(fit,
			AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, results);
}

