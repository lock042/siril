/* @(#)pave.c	19.1 (ES0-DMD) 02/25/03 13:34:39 */
/*===========================================================================
  Copyright (C) 1995 European Southern Observatory (ESO)

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public
  License along with this program;
  If not, see <http://www.gnu.org/licenses/>.

  Corresponding concerning ESO-MIDAS should be addressed as follows:
	Internet e-mail: midas@eso.org
	Postal address: European Southern Observatory
			Data Management Division
			Karl-Schwarzschild-Strasse 2
			D 85748 Garching bei Muenchen
			GERMANY
===========================================================================*/

/******************************************************************************
**                   Copyright (C) 1993 by European Southern Observatory
*******************************************************************************
**
**    UNIT
**
**    Version: 19.1
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/25
**
**    File:  pave.c
**
*******************************************************************************
**
**    DESCRIPTION  routines for the a trous algorithme
**    ----------
*******************************************************************************
**
** pave_2d_tfo (Imag, Pave, Nl, Nc, Nbr_Plan, Type_To)
** float *Imag, *Pave;
** int Nl, Nc, Nbr_Plan;
** int Type_To;
**
** computes the wavelet transform without reduction of the sampling
**
** Type_To = TO_PAVE_LINEAR for a linear scaling function
** Type_To = TO_PAVE_BSPLINE for a b3-spline scaling function
**
*******************************************************************************
**
** pave_2d_build (Pave, Imag, Nl, Nc, Nbr_Plan)
** float *Imag, *Pave;
** int Nl, Nc, Nbr_Plan;
**
** reconstruction of the image from its wavelet transform
**
*******************************************************************************
**
** pave_2d_extract_plan (Pave, Imag, Nl, Nc, Num_Plan)
** float *Imag, *Pave;
** int Nl, Nc, Num_Plan;
**
** extracts a plan from the wavelet transform
**
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "algos/Def_Math.h"
#include "algos/Def_Mem.h"
#include "algos/Def_Wavelet.h"

/* Half-widths and 1D taps of the à trous scaling functions, given from the
 * centre outwards (t[0] is the centre, t[k] the weight of the ±k*Step
 * neighbours). The 2D kernels are the outer products of these with themselves:
 *
 *   linear   [1 2 1]/4        ->  the 3x3 kernel  1/4, 1/8, 1/16
 *   b3-spline [1 4 6 4 1]/16  ->  the 5x5 kernel  0.140625, 0.09375, 0.0625,
 *                                 0.0234375, 0.015625, 0.00390625
 *
 * so smoothing separably (one horizontal then one vertical pass) computes
 * exactly the same filter as a direct 2D convolution — 6 taps per pixel
 * instead of 9, and 10 instead of 25 — differing only in the order of the
 * float additions. Replicated borders stay exact under separation because the
 * row and column clamps are independent of each other. */
#define PAVE_LINEAR_RADIUS 1
#define PAVE_BSPLINE_RADIUS 2
static const float pave_linear_taps[PAVE_LINEAR_RADIUS + 1] = { 0.5f, 0.25f };
static const float pave_bspline_taps[PAVE_BSPLINE_RADIUS + 1] = { 0.375f, 0.25f,
		0.0625f };

/* Border handling: replicate the edge pixel. */
static inline int clamp_ind(int ind, int N) {
	if (ind < 0)
		return 0;
	if (ind >= N)
		return N - 1;
	return ind;
}

/****************************************************************************/

/* Horizontal half of the separable smoothing. Rows are independent, and the
 * middle of each row needs no index clamping, so that part is a plain
 * unit-stride loop the compiler can vectorise. */
static void pave_hpass(const float *src, float *dst, int Nl, int Nc, int Step,
		int r, const float *t, int threads) {
	int lo = r * Step;
	int hi = Nc - r * Step;
	if (lo > Nc)
		lo = Nc;
	if (hi < lo)
		hi = lo;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for (int i = 0; i < Nl; i++) {
		const float *s = src + (size_t) i * Nc;
		float *d = dst + (size_t) i * Nc;

		/* left and right borders: the only pixels that need the clamp */
		for (int j = 0; j < lo; j++) {
			float v = t[0] * s[j];
			for (int k = 1; k <= r; k++)
				v += t[k] * (s[clamp_ind(j - k * Step, Nc)]
						+ s[clamp_ind(j + k * Step, Nc)]);
			d[j] = v;
		}
		for (int j = hi; j < Nc; j++) {
			float v = t[0] * s[j];
			for (int k = 1; k <= r; k++)
				v += t[k] * (s[clamp_ind(j - k * Step, Nc)]
						+ s[clamp_ind(j + k * Step, Nc)]);
			d[j] = v;
		}

		/* interior */
		if (r == PAVE_BSPLINE_RADIUS) {
			const float t0 = t[0], t1 = t[1], t2 = t[2];
			const int s1 = Step, s2 = 2 * Step;
#ifdef _OPENMP
#pragma omp simd
#endif
			for (int j = lo; j < hi; j++)
				d[j] = t0 * s[j] + t1 * (s[j - s1] + s[j + s1])
						+ t2 * (s[j - s2] + s[j + s2]);
		} else {
			const float t0 = t[0], t1 = t[1];
			const int s1 = Step;
#ifdef _OPENMP
#pragma omp simd
#endif
			for (int j = lo; j < hi; j++)
				d[j] = t0 * s[j] + t1 * (s[j - s1] + s[j + s1]);
		}
	}
}

/* Vertical half of the separable smoothing. Each output row reads 2r+1 source
 * rows; those row pointers are hoisted out of the pixel loop, which is then
 * unit-stride over the whole row. */
static void pave_vpass(const float *src, float *dst, int Nl, int Nc, int Step,
		int r, const float *t, int threads) {
	int lo = r * Step;
	int hi = Nl - r * Step;
	if (lo > Nl)
		lo = Nl;
	if (hi < lo)
		hi = lo;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for (int i = 0; i < Nl; i++) {
		float *d = dst + (size_t) i * Nc;
		const float *c = src + (size_t) i * Nc;

		if (i >= lo && i < hi) {
			if (r == PAVE_BSPLINE_RADIUS) {
				const float *m1 = src + (size_t) (i - Step) * Nc;
				const float *p1 = src + (size_t) (i + Step) * Nc;
				const float *m2 = src + (size_t) (i - 2 * Step) * Nc;
				const float *p2 = src + (size_t) (i + 2 * Step) * Nc;
				const float t0 = t[0], t1 = t[1], t2 = t[2];
#ifdef _OPENMP
#pragma omp simd
#endif
				for (int j = 0; j < Nc; j++)
					d[j] = t0 * c[j] + t1 * (m1[j] + p1[j])
							+ t2 * (m2[j] + p2[j]);
			} else {
				const float *m1 = src + (size_t) (i - Step) * Nc;
				const float *p1 = src + (size_t) (i + Step) * Nc;
				const float t0 = t[0], t1 = t[1];
#ifdef _OPENMP
#pragma omp simd
#endif
				for (int j = 0; j < Nc; j++)
					d[j] = t0 * c[j] + t1 * (m1[j] + p1[j]);
			}
		} else {
			/* top and bottom bands: clamp the row indices once per row */
			const float *rm[PAVE_BSPLINE_RADIUS], *rp[PAVE_BSPLINE_RADIUS];
			for (int k = 1; k <= r; k++) {
				rm[k - 1] = src + (size_t) clamp_ind(i - k * Step, Nl) * Nc;
				rp[k - 1] = src + (size_t) clamp_ind(i + k * Step, Nl) * Nc;
			}
			for (int j = 0; j < Nc; j++) {
				float v = t[0] * c[j];
				for (int k = 1; k <= r; k++)
					v += t[k] * (rm[k - 1][j] + rp[k - 1][j]);
				d[j] = v;
			}
		}
	}
}

/* Smooth src into dst at scale Num_Plan, using scratch (Nl*Nc floats) for the
 * intermediate of the separable pass. scratch must not alias src or dst. */
static int pave_2d_smooth(const float *src, float *dst, float *scratch, int Nl,
		int Nc, int Num_Plan, int Type_To, int threads) {
	const int Step = 1 << Num_Plan;
	int r;
	const float *t;

	switch (Type_To) {
	case TO_PAVE_LINEAR:
		r = PAVE_LINEAR_RADIUS;
		t = pave_linear_taps;
		break;
	case TO_PAVE_BSPLINE:
		r = PAVE_BSPLINE_RADIUS;
		t = pave_bspline_taps;
		break;
	default:
		siril_log_message(_("pave_2d_smooth: unknown transform\n"));
		return 1;
	}

	pave_hpass(src, scratch, Nl, Nc, Step, r, t, threads);
	pave_vpass(scratch, dst, Nl, Nc, Step, r, t, threads);
	return 0;
}

/***************************************************************************/

int pave_2d_tfo(float *Pict, float *Pave, int Nl, int Nc, int Nbr_Plan,
		int Type_To, int threads) {
	threads = wavelet_threads(threads);
	if (Type_To != TO_PAVE_LINEAR && Type_To != TO_PAVE_BSPLINE) {
		siril_log_message(_("pave_2d_tfo: unknown transform\n"));
		return 1;
	}

	const size_t N = (size_t) Nl * (size_t) Nc;

	/* Pict belongs to the caller — for float images it is the live image
	 * buffer — so the recursion runs on our own copy. Two scratch planes are
	 * enough: the detail plane being produced doubles as the intermediate of
	 * the separable smoothing, and is overwritten with the detail
	 * coefficients immediately afterwards. */
	float *cur = malloc(N * sizeof(float));
	float *next = malloc(N * sizeof(float));
	if (!cur || !next) {
		free(cur);
		free(next);
		PRINT_ALLOC_ERR;
		return 1;
	}
	memcpy(cur, Pict, N * sizeof(float));

	for (int Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++) {
		float *Plan = Pave + N * (size_t) Num_Plan;

		if (pave_2d_smooth(cur, next, Plan, Nl, Nc, Num_Plan, Type_To, threads)) {
			free(cur);
			free(next);
			return 1;
		}

		/* the wavelet plane is the difference between two resolutions */
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(threads) schedule(static)
#endif
		for (size_t i = 0; i < N; i++)
			Plan[i] = cur[i] - next[i];

		float *swap = cur;
		cur = next;
		next = swap;
	}

	/* copy the low resolution image in the transform */
	memcpy(Pave + N * (size_t) (Nbr_Plan - 1), cur, N * sizeof(float));

	free(cur);
	free(next);
	return 0;
}

/***************************************************************************/

/* Number of pixels reconstructed per block. The running sum for a block stays
 * in L1 while every plane contributes to it, so the output is written to
 * memory once instead of being read-modify-written once per plane. */
#define PAVE_BUILD_BLOCK 4096

int pave_2d_build(float *Pave, float *Imag, int Nl, int Nc, int Nbr_Plan,
		const float *coef, int threads) {
	threads = wavelet_threads(threads);
	const size_t N = (size_t) Nl * (size_t) Nc;
	const size_t nblocks = (N + PAVE_BUILD_BLOCK - 1) / PAVE_BUILD_BLOCK;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for (size_t b = 0; b < nblocks; b++) {
		const size_t beg = b * PAVE_BUILD_BLOCK;
		size_t end = beg + PAVE_BUILD_BLOCK;
		if (end > N)
			end = N;
		float *out = Imag + beg;
		const size_t len = end - beg;
		int started = 0;

		for (int Num_Plan = 0; Num_Plan < Nbr_Plan; Num_Plan++) {
			const float c = coef[Num_Plan];
			if (c == 0.f)
				continue; /* a zeroed scale contributes nothing */
			const float *plan = Pave + N * (size_t) Num_Plan + beg;
			if (!started) {
#ifdef _OPENMP
#pragma omp simd
#endif
				for (size_t i = 0; i < len; i++)
					out[i] = c * plan[i];
				started = 1;
			} else {
#ifdef _OPENMP
#pragma omp simd
#endif
				for (size_t i = 0; i < len; i++)
					out[i] += c * plan[i];
			}
		}
		if (!started)
			memset(out, 0, len * sizeof(float));
	}
	return 0;
}

/***************************************************************************/

int pave_2d_extract_plan(float *Pave, float *Imag, int Nl, int Nc, int Num_Plan) {
	const size_t N = (size_t) Nl * (size_t) Nc;

	memcpy(Imag, Pave + N * (size_t) Num_Plan, N * sizeof(float));

	return 0;
}

/****************************************************************************/
