/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics_float.h>

#include "core/siril.h"
#include "stacking/siril_fit_linear.h"
#include "stacking/stacking.h"
#include "algos/sorting.h"
#include "algos/statistics.h"

static int percentile_clipping(float pixel, float sig[], float median,
		guint64 rej[]) {
	float plow = sig[0];
	float phigh = sig[1];

	if (median - pixel > median * plow) {
		rej[0]++;
		return -1;
	} else if (pixel - median > median * phigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

/* Rejection of pixels, following sigma_(high/low) * sigma.
 * The function returns 0 if no rejections are required, 1 if it's a high
 * rejection and -1 for a low-rejection */
static int sigma_clipping_float(float pixel, float sigma, float sigmalow,
		float sigmahigh, float median, guint64 rej[]) {

	if (median - pixel > sigma * sigmalow) {
		rej[0]++;
		return -1;
	} else if (pixel - median > sigma * sigmahigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

static int line_clipping(float pixel, float sig[], float sigma, int i, float a,
		float b, guint64 rej[]) {
	float sigmalow = sig[0];
	float sigmahigh = sig[1];

	if (a * i + b - pixel> sigma * sigmalow) {
		rej[0]++;
		return -1;
	} else if (pixel - a * i - b > sigma * sigmahigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

int apply_rejection_float(struct _data_block *data, int nb_frames,
		struct stacking_args *args, guint64 crej[2]) {
	int N = nb_frames;	// N is the number of pixels kept from the current stack
	double median = 0.0;
	int pixel, output, changed, n, r = 0;
	int firstloop = 1;

	float *stack = (float*) data->stack;
	float *w_stack = (float*) data->w_stack;
	int *rejected = (int*) data->rejected;
	const float siglow = args->sig[0];
	const float sighigh = args->sig[1];

	/* prepare median and check that the stack is not mostly zero */
	switch (args->type_of_rejection) {
	case PERCENTILE:
	case SIGMA:
		median = quickmedian_float(stack, N);
		if (median == 0.0)
			return 0;
		break;
	case SIGMEDIAN:
	case WINSORIZED:
	default:
		break;
	}

	switch (args->type_of_rejection) {
	case PERCENTILE:
		for (int frame = 0; frame < N; frame++) {
			rejected[frame] = percentile_clipping(stack[frame], args->sig, (float) median, crej);
		}
		for (pixel = 0, output = 0; pixel < N; pixel++) {
			if (!rejected[pixel]) {
				// copy only if there was a rejection
				if (pixel != output)
					stack[output] = stack[pixel];
				output++;
			}
		}
		N = output;
		break;
	case SIGMA:
		do {
			const float sigma = siril_stats_float_sd(stack, N);
			if (!firstloop)
				median = quickmedian_float(stack, N);
			else
				firstloop = 0;
			for (int frame = 0; frame < N; frame++) {
				if (N - r <= 4) {
					// no more rejections
					rejected[frame] = 0;
				} else {
					rejected[frame] = sigma_clipping_float(stack[frame], sigma,
							siglow, sighigh, (float) median, crej);
					if (rejected[frame])
						r++;
				}
			}
			for (pixel = 0, output = 0; pixel < N; pixel++) {
				if (!rejected[pixel]) {
					// copy only if there was a rejection
					if (pixel != output)
						stack[output] = stack[pixel];
					output++;
				}
			}
			changed = N != output;
			N = output;
		} while (changed && N > 3);
		break;
	case SIGMEDIAN:
		do {
			const float sigma = siril_stats_float_sd(stack, N);
			const float medianf = quickmedian_float(stack, N);
			n = 0;
			for (int frame = 0; frame < N; frame++) {
				if (sigma_clipping_float(stack[frame], sigma, siglow, sighigh, medianf, crej)) {
					stack[frame] = medianf;
					n++;
				}
			}
		} while (n > 0);
		break;
	case WINSORIZED:
		do {
			float sigma0;
			float sigma = siril_stats_float_sd(stack, N);
			const float medianf = quickmedian_float(stack, N);
			memcpy(w_stack, stack, N * sizeof(float));
			do {
				const float m0 = medianf - 1.5f * sigma;
				const float m1 = medianf + 1.5f * sigma;
				for (int jj = 0; jj < N; jj++) {
					w_stack[jj] = min(m1, max(m0, w_stack[jj]));
				}
				sigma0 = sigma;
				sigma = 1.134f * siril_stats_float_sd(w_stack, N);
			} while (fabsf(sigma - sigma0) > sigma0 * 0.0005f);
			for (int frame = 0; frame < N; frame++) {
				if (N - r <= 4) {
					// no more rejections
					rejected[frame] = 0;
				} else {
					rejected[frame] = sigma_clipping_float(stack[frame], sigma, siglow, sighigh, medianf, crej);
					if (rejected[frame] != 0)
						r++;
				}

			}
			for (pixel = 0, output = 0; pixel < N; pixel++) {
				if (!rejected[pixel]) {
					// copy only if there was a rejection
					stack[output] = stack[pixel];
					output++;
				}
			}
			changed = N != output;
			N = output;
		} while (changed && N > 3);
		break;
	case LINEARFIT:
		do {
			quicksort_f(stack, N);
			for (int frame = 0; frame < N; frame++) {
				data->yf[frame] = stack[frame];
			}
			float a, b;
			siril_fit_linear(data->xf, data->yf, data->m_x, data->m_dx2, N, &b,	&a);
			float sigma = 0.f;
			for (int frame = 0; frame < N; frame++)
				sigma += fabsf(stack[frame] - (a * frame + b));
			sigma /= (float) N;
			for (int frame = 0; frame < N; frame++) {
				if (N - r <= 4) {
					// no more rejections
					rejected[frame] = 0;
				} else {
					rejected[frame] = line_clipping(stack[frame], args->sig,
							sigma, frame, a, b, crej);
					if (rejected[frame] != 0)
						r++;
				}
			}
			for (pixel = 0, output = 0; pixel < N; pixel++) {
				if (!rejected[pixel]) {
					// copy only if there was a rejection
					if (pixel != output)
						stack[output] = stack[pixel];
					output++;
				}
			}
			changed = N != output;
			N = output;
		} while (changed && N > 3);
		break;
	default:
	case NO_REJEC:
		;		// Nothing to do, no rejection
	}
	return N;
}

