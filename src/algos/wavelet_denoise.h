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

#ifndef SRC_ALGOS_WAVELET_DENOISE_H_
#define SRC_ALGOS_WAVELET_DENOISE_H_

#include <stddef.h>
#include <glib.h>

/* Per-scale wavelet denoising — noise foundation.
 *
 * Siril's wavelet tool uses the isotropic undecimated starlet (a trous)
 * transform (see src/algos/pave.c). For white Gaussian noise of unit variance
 * the standard deviation of detail scale j is a fixed "noise propagation
 * factor" e_j that depends only on the scaling function and the scale, not on
 * the image. Combined with a single global noise estimate sigma, the per-scale
 * noise std is sigma_j = sigma * e_j, which is the basis of k-sigma
 * thresholding and bivariate shrinkage on this transform. */

/* Maximum number of planes supported by the factor cache. */
#define WD_MAX_PLAN 16

/* Compute the per-scale noise propagation factors e_j for the given transform
 * (TO_PAVE_LINEAR or TO_PAVE_BSPLINE) and number of planes. The factors are the
 * L2 norms of the equivalent detail filters, obtained exactly by transforming a
 * unit impulse through the actual kernels. e_out must hold at least
 * (nbr_plan - 1) doubles; it is filled for the detail scales 0..nbr_plan-2 (the
 * coarsest plane is the residual and has no factor). Results are cached per
 * (type, nbr_plan). Returns 0 on success, non-zero on bad arguments/allocation. */
int wavelet_noise_factors(int type, int nbr_plan, double *e_out);

/* Robust (MAD-based) Gaussian noise std of a coefficient band of n samples:
 * 1.4826 * median(|x - median(x)|). Returns < 0 on allocation failure. The
 * input band is not modified. */
double wavelet_mad_sigma_float(const float *band, size_t n);

/* Estimate the global image noise sigma from the finest detail plane (band0,
 * n samples) and that scale's propagation factor e1: sigma = MAD_sigma / e1.
 * Returns < 0 on error (e1 <= 0 or allocation failure). */
double wavelet_estimate_noise_float(const float *band0, size_t n, double e1);

/* Estimate the global noise sigma of an on-disk a trous transform (.wave file)
 * from its finest detail plane. Writes the result to *sigma_out. Returns 0 on
 * success, non-zero on read/estimation failure. */
int wavelet_sigma_from_file(const char *filename, double *sigma_out);

/* Shrinkage method. WD_THRESHOLD is implemented now; WD_BISHRINK (the intended
 * default) and WD_GSM are added in later phases. */
enum wd_method {
	WD_THRESHOLD = 0,
	WD_BISHRINK  = 1,
	WD_GSM       = 2,
};

/* Source of the per-scale noise standard deviation. PROPAGATED uses a single
 * global sigma scaled by the per-scale propagation factors e_j (the default,
 * "obey the propagation factors"); PER_BAND measures each band's own MAD. */
enum wd_sigma_source {
	WD_SIGMA_PROPAGATED = 0,
	WD_SIGMA_PER_BAND   = 1,
};

/* Per-scale denoising parameters. The effective per-scale threshold is
 * t_j = k * f[j] * sigma_j, with sigma_j = sigma_global * e_j by default.
 * Default f[j] = 1 obeys the propagation factors; raise f[0] to crush the
 * finest scale and lower f[j] to be gentle on coarser scales. */
struct denoise_params {
	gboolean enabled;     /* master switch */
	int method;           /* enum wd_method */
	float k;              /* global k-sigma multiplier (default 3) */
	float f[WD_MAX_PLAN]; /* per detail-scale factor, default 1 */
	int sigma_source;     /* enum wd_sigma_source */
	gboolean soft;        /* soft vs hard thresholding (WD_THRESHOLD) */
	gboolean anscombe;    /* the transform was built in the Anscombe VST domain
	                       * (set at decomposition); reconstruction must invert it */
};

/* Anscombe variance-stabilising transform. Applied to the image (in count-like
 * units, hence the scale) before the linear wavelet decomposition so that
 * signal-dependent (Poisson) noise becomes approximately Gaussian with constant
 * variance, which is what the shrinkage assumes. anscombe_inverse is the exact
 * algebraic inverse, so VST with no denoising round-trips to the original. */
#define ANSCOMBE_USHORT_SCALE 1.0    /* USHORT wavelet buffer is already in ADU */
#define ANSCOMBE_FLOAT_SCALE  65535.0 /* map normalised [0,1] data to ADU-like */
void anscombe_forward(float *data, size_t n, double scale);
void anscombe_inverse(float *data, size_t n, double scale);

/* Initialise dp to safe defaults: disabled, threshold/soft, k = 3, all
 * per-scale factors = 1, propagated sigma. Callers tweak fields afterwards. */
void denoise_params_init(struct denoise_params *dp);

/* Shrink the detail planes of an in-memory a trous transform in place. pave_data
 * is the contiguous plane buffer (plane p at offset p*Nl*Nc); type is
 * TO_PAVE_LINEAR/TO_PAVE_BSPLINE; the coarsest plane (nbr_plan-1, the residual)
 * is left untouched. No-op when dp is NULL or disabled. Returns 0 on success. */
int wavelet_denoise_planes(float *pave_data, int type, int nbr_plan, int Nl,
		int Nc, const struct denoise_params *dp);

#endif /* SRC_ALGOS_WAVELET_DENOISE_H_ */
