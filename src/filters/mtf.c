/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include <math.h>
#include "mtf.h"
#include "algos/sorting.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "algos/statistics.h"

void apply_linked_mtf_to_fits(fits *from, fits *to, struct mtf_params params, gboolean multithreaded) {

	const gboolean do_channel[3] = { params.do_red, params.do_green, params.do_blue };

	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);

	siril_log_message(_("Applying MTF with values %f, %f, %f\n"),
			params.shadows, params.midtones, params.highlights);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				for (size_t i = 0; i < layersize; i++) {
					float pxl = from->pdata[j][i] * invnorm;
					float mtf = MTFp(pxl, params);
					to->pdata[j][i] = roundf_to_WORD(mtf * norm);
				}
			} else
				memcpy(to->pdata[j], from->pdata[j], layersize * sizeof(WORD));
		}
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				for (size_t i = 0; i < layersize; i++) {
					to->fpdata[j][i] = MTFp(from->fpdata[j][i], params);
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

	siril_log_message(_("Applying inverse MTF with values %f, %f, %f\n"),
			params.shadows, params.midtones, params.highlights);
	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
#endif
		for (size_t j = 0; j < from->naxes[2] ; j++) {
			if (do_channel[j]) {
				for (size_t i = 0; i < layersize; i++) {
					float pxl = from->pdata[j][i] * invnorm;
					float mtf = MTF_pseudoinverse(pxl, params);
					to->pdata[j][i] = roundf_to_WORD(mtf * norm);
				}
			} else
				memcpy(to->pdata[j], from->pdata[j], layersize * sizeof(WORD));
		}
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
	cmsHTRANSFORM transform = NULL;

	int nb_channels = (int)fit->naxes[2];

	if (com.icc.available) {
		fits_check_icc(fit);
		if (nb_channels == 1)
			// This is a bit of a fudge: it should really be using the
			// display profile and TYPE_RGB_FLT_PLANAR, but only bothering
			// to use the first channel
			transform = cmsCreateTransform(fit->icc_profile, TYPE_GRAY_FLT, com.icc.mono_standard, TYPE_GRAY_FLT, gui.icc.rendering_intent, 0);
		else
			transform = cmsCreateTransform(fit->icc_profile, TYPE_RGB_FLT_PLANAR, gui.icc.monitor, TYPE_RGB_FLT_PLANAR, gui.icc.rendering_intent, 0);
	}
	int retval = 0;
	if (gui.icc.available)
		retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_CDF, MULTI_THREADED, stat);
	else
		retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);

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

	float *invnormValue = malloc(nb_channels * sizeof(float));
	float *median = malloc(nb_channels * sizeof(float));
	float *mad = malloc(nb_channels * sizeof(float));
	float *percentiles = malloc(nb_channels * NBUCKETS * sizeof(float));
	for (i = 0; i < nb_channels; ++i) {
		invnormValue[i] = 1.f / (float)stat[i]->normValue;
		median[i] = (float) stat[i]->median * invnormValue[i];
		// if color management is active we need to transform the median
		// and calculate the MAD in the display colorspace
		if (gui.icc.available) {
			int start = i * NBUCKETS;
			for (int j = 0; j < NBUCKETS; j++)
				percentiles[start + j] = stat[i]->cdf[j] * invnormValue[i]; 

		} else {
			mad[i] = stat[i]->mad * invnormValue[i]; 
		}
	}
	if (gui.icc.available) {
		cmsDoTransform(transform, (void *) median, (void *) median, 1);
		cmsDoTransform(transform, (void*) percentiles, (void*) percentiles, NBUCKETS);
	}
	for (i = 0 ; i < nb_channels ; i++) {
		if (gui.icc.available) {
			float *devs = malloc(NBUCKETS * sizeof(float));
			for (int j = 0 ; j < NBUCKETS ; j++) {
				devs[j] = fabsf(*(percentiles + (i * NBUCKETS) + j) - median[i]);
			}
			mad[i] = quickmedian_float(devs, NBUCKETS);
			free(devs);
		}
		/* this is a guard to avoid breakdown point */
		if (mad[i] == 0.f) mad[i] = 0.001f;
		printf("layer:%d - median: %8.6f - mad: %8.6f\n", i, median[i], mad[i]);
	}
	if (invertedChannels < nb_channels) {
		for (i = 0 ; i < nb_channels ; i++) {
			c0 += median[i] + shadows_clipping * mad[i];
			m += median[i];
		}
		c0 /= (float) nb_channels;
		if (c0 < 0.f) c0 = 0.f;
		float m2 = m / (float) nb_channels - c0;
		result->midtones = MTF(m2, target_bg, 0.f, 1.f);
		result->shadows = c0;
		result->highlights = 1.0f;
	} else {
		for (i = 0 ; i < nb_channels ; ++i) {
			c1 += median[i] - shadows_clipping * mad[i];
			m += median[i];
		}
		c1 /= (float) nb_channels;
		if (c1 > 1.f) c1 = 1.f;
		float m2 = c1 - m / (float) nb_channels;
		result->midtones = 1.f - MTF(m2, target_bg, 0.f, 1.f);
		result->shadows = 0.f;
		result->highlights = c1;
	}
	free(percentiles);
	free(median);
	free(invnormValue);
	free(mad);

	for (i = 0; i < nb_channels; i++)
		free_stats(stat[i]);
	return 0;

}

int find_linked_midtones_balance_default(fits *fit, struct mtf_params *result) {
	return find_linked_midtones_balance(fit,
			AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, result);
}

void apply_unlinked_mtf_to_fits(fits *from, fits *to, struct mtf_params *params) {
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
		float invnorm = 1.0f / norm;
#ifdef _OPENMP
		int threads = com.max_thread >= 3 ? 3 : com.max_thread;
#endif
		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
			siril_log_message(_("Applying MTF to channel %d with values %f, %f, %f\n"), chan,
					params[chan].shadows, params[chan].midtones, params[chan].highlights);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
			for (size_t i = 0; i < ndata; i++) {
				float pxl = (float)from->pdata[chan][i] * invnorm;
				float mtf = MTFp(pxl, params[chan]);
				to->pdata[chan][i] = roundf_to_WORD(mtf * norm);
			}
		}
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
		int threads = com.max_thread >= 3 ? 3 : com.max_thread;
#endif
		for (int chan = 0; chan < (int)from->naxes[2]; chan++) {
			siril_log_message(_("Applying MTF to channel %d with values %f, %f, %f\n"), chan,
					params[chan].shadows, params[chan].midtones, params[chan].highlights);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
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
	float c0 = 0.0, c1 = 0.0;
	float m = 0.0;
	int i, invertedChannels = 0;
	imstats *stat[3];
	cmsHTRANSFORM transform = NULL;

	int nb_channels = (int)fit->naxes[2];

	if (com.icc.available) {
		fits_check_icc(fit);
		if (nb_channels == 1)
			// Fudge, see the linked version above
			transform = cmsCreateTransform(fit->icc_profile, TYPE_GRAY_FLT, com.icc.mono_standard, TYPE_GRAY_FLT, gui.icc.rendering_intent, 0);
		else
			transform = cmsCreateTransform(fit->icc_profile, TYPE_RGB_FLT_PLANAR, gui.icc.monitor, TYPE_RGB_FLT_PLANAR, gui.icc.rendering_intent, 0);
	}

	int retval = 0;
	if (gui.icc.available)
		retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_CDF, MULTI_THREADED, stat);
	else
		retval = compute_all_channels_statistics_single_image(fit, STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
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

	float *invnormValue = malloc(nb_channels * sizeof(float));
	float* median = malloc(nb_channels * sizeof(float));
	float* mad = malloc(nb_channels * sizeof(float));
	float *percentiles = malloc(nb_channels * NBUCKETS * sizeof(float));
	for (i = 0; i < nb_channels; ++i) {
		invnormValue[i] = 1.f / (float)stat[i]->normValue;
		median[i] = (float) stat[i]->median * invnormValue[i];
		// if color management is active we need to transform the median
		// and calculate the MAD in the display colorspace
		if (gui.icc.available) {
			int start = i * NBUCKETS;
			for (int j = 0; j < NBUCKETS; j++)
				percentiles[start + j] = stat[i]->cdf[j] * invnormValue[i]; 

		} else {
			mad[i] = stat[i]->mad * invnormValue[i]; 
		}
	}
	if (gui.icc.available) {
		cmsDoTransform(transform, (void *) median, (void *) median, 1);
		cmsDoTransform(transform, (void*) percentiles, (void*) percentiles, NBUCKETS);
	}
	for (i = 0 ; i < nb_channels ; ++i) {
		if (gui.icc.available) {
			float *devs = malloc(NBUCKETS * sizeof(float));
			for (int j = 0 ; j < NBUCKETS ; j++) {
				devs[j] = fabsf(*(percentiles + (i * NBUCKETS) + j) - median[i]);
			}
			mad[i] = quickmedian_float(devs, NBUCKETS);
			free(devs);
		}
		/* this is a guard to avoid breakdown point */
		if (mad[i] == 0.f) mad[i] = 0.001f;
		printf("layer:%d - median: %8.6f - mad: %8.6f\n", i, median[i], mad[i]);
		if (invertedChannels < nb_channels) {
			c0 = median[i] + shadows_clipping * mad[i];
			m = median[i];
			c0 /= (float) nb_channels;
			if (c0 < 0.f) c0 = 0.f;
			float m2 = m / (float) nb_channels - c0;
			results[i].midtones = MTF(m2, target_bg, 0.f, 1.f);
			results[i].shadows = c0;
			results[i].highlights = 1.0f;
		} else {
			c1 = median[i] - shadows_clipping * mad[i];
			m = median[i];
			c1 /= (float) nb_channels;
			if (c1 > 1.f) c1 = 1.f;
			float m2 = c1 - m / (float) nb_channels;
			results[i].midtones = 1.f - MTF(m2, target_bg, 0.f, 1.f);
			results[i].shadows = 0.f;
			results[i].highlights = c1;
		}
	}
	free(percentiles);
	free(median);
	free(invnormValue);
	free(mad);

	for (i = 0; i < nb_channels; i++)
		free_stats(stat[i]);
	return 0;
}

int find_unlinked_midtones_balance_default(fits *fit, struct mtf_params *results) {
	return find_unlinked_midtones_balance(fit,
			AS_DEFAULT_SHADOWS_CLIPPING, AS_DEFAULT_TARGET_BACKGROUND, results);
}

