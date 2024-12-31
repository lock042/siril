/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_world_cs.h"
#include "algos/photometry.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "algos/star_finder.h"
#include "filters/median.h"

#include "PSF.h"


#define MAX_ITER_NO_ANGLE  20		//Number of iterations in the minimization with no angle
#define MAX_ITER_ANGLE     20		//Number of iterations in the minimization with angle
#define MIN_HALF_RADIUS     1		// Minimum radius around center pixel to initialize FWHM
#define MIN_LARGE_SAMPLING     3	// Minimum FWHM to use a more sophisticated PSF initialization
#define XTOL 1e-3
#define GTOL 1e-3
#define FTOL 1e-3

#define DEBUG_PSF 0 // flag to show progress of fitting process - may flood output if numerous stars

const double radian_conversion = ((3600.0 * 180.0) / M_PI) / 1.0E3;

// we also zero at bg level so that we don't have to bother substracting bg
// in all the subsequent operations
static gsl_matrix *prepare_init_matrix(gsl_matrix *in, double bg, gboolean frompeaker) {
	size_t width = in->size2;
	size_t height = in->size1;
	size_t x, y;
	gsl_matrix *out = gsl_matrix_alloc (in->size1, in->size2);
	if (!out) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	gsl_matrix_memcpy (out, in);
	if (frompeaker) {
		gsl_matrix_add_constant(out, -bg);
		return out;
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			double a = get_median_gsl(in, x, y, width, height, 1, FALSE, FALSE);
			gsl_matrix_set(out, y, x, a - bg);
		}
	}
	return out;
}

static double S_from_FWHM(double FWHM, double beta, starprofile profile) {
	if (profile == PSF_GAUSSIAN)
		return SQR(FWHM) * INV_4_LOG2;
	return SQR(FWHM) * 0.25 / (pow(2., 1. / beta) - 1.);
}

#if DEBUG_PSF
static double FWHM_from_S(double S, double beta, starprofile profile) {
	if (profile == PSF_GAUSSIAN)
		return 2. * sqrt(S * log(2.));
	return 2. * sqrt(S * (pow(2., 1./ beta) - 1.));
}
#endif

// static double s_from_FWHM(double FWHM, double beta, starprofile profile) {
// 	if (profile == PSF_GAUSSIAN)
// 		return FWHM / _2_SQRT_2_LOG2;
// 	return FWHM * 0.5 /  sqrt(pow(2., 1./ beta) - 1.);
// }

static double FWHM_from_s(double s, double beta, starprofile profile) {
	if (profile == PSF_GAUSSIAN)
		return s * _2_SQRT_2_LOG2;
	return 2. * s * sqrt(pow(2., 1./ beta) - 1.);
}

/* Compute initial values for the algorithm from data in the pixel value matrix */
static gsl_vector* psf_init_data(gsl_matrix* z, double bg, gboolean frompeaker) {
	gsl_vector *MaxV = gsl_vector_alloc(6);
	if (!MaxV) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	double max, halfA;
	double S = 0., Ixx = 0., Iyy = 0., Ixy = 0., Sx = 0., Sy = 0.;
	size_t NbRows = z->size1; //y
	size_t NbCols = z->size2; //x
	size_t i, j;

	/* find maximum */
	/* first we remove hot pixels in the matrix if not coming from the peaker
	   to make sure we get a real maximum. Otherwise we just remove the bakcground
	*/
	gsl_matrix *m_tmp = prepare_init_matrix(z, bg, frompeaker);
	if (!m_tmp) return NULL;
	// if the data is coming from the peaker, we know the central pixel is already
	// the max intensity pixel. We don't bother searching it as it may catch another
	// brighter star located in the box
	if (!frompeaker) {
		max = gsl_matrix_max(m_tmp);
		gsl_matrix_max_index(m_tmp, &i, &j); // i=y , j=x
	} else {
		max = gsl_matrix_get(m_tmp, NbRows / 2, NbCols / 2);
		i = NbRows / 2;
		j = NbCols / 2;
	}
	halfA = max * 0.5; // half amplitude
	gsl_vector_set(MaxV, 0, max); // A

	int yc = (int)i;
	int xc = (int)j;
	int ii1 = yc + MIN_HALF_RADIUS;
	int ii2 = yc - MIN_HALF_RADIUS;
	int jj1 = xc + MIN_HALF_RADIUS;
	int jj2 = xc - MIN_HALF_RADIUS;

	while (ii1 < NbRows - 1 && gsl_matrix_get(m_tmp, ii1 + 1, xc) > halfA) {
		ii1++;
	}
	while (ii2 > 0 && gsl_matrix_get(m_tmp, ii2 - 1, xc) > halfA) {
		ii2--;
	}

	while (jj1 < NbCols - 1 && gsl_matrix_get(m_tmp, yc, jj1 + 1) > halfA) {
		jj1++;
	}
	while (jj2 > 0 && gsl_matrix_get(m_tmp, yc, jj2 - 1) > halfA) {
		jj2--;
	}
	double FWHMx = jj1 - jj2;
	double FWHMy = ii1 - ii2;

	if (min(FWHMx, FWHMy) < MIN_LARGE_SAMPLING) {
		gsl_vector_set(MaxV, 1, (double)(jj1 + jj2 + 1) / 2.0); //x0
		gsl_vector_set(MaxV, 2, (double)(ii1 + ii2 + 1) / 2.0); //y0
		gsl_vector_set(MaxV, 3, max(FWHMx, FWHMy)); //FWHM x
		gsl_vector_set(MaxV, 4, min(FWHMx, FWHMy)); //FWHM y
		gsl_vector_set(MaxV, 5, 0.); //angle
		gsl_matrix_free(m_tmp);
		return MaxV;
	}

	// If the star is large enough, we will make a smarter init for the fitting:
	// This is based on computing centroid, fwhm, roundness and angle of the trace
	// of the star above the half amplitude. We use the trace instead of the actual
	// data in order to be resilient to outlier values.
	// In order to only account for the star (and not data from an adjacent star
	// also in the box), we use the central row of the trace first, and then move
	// downwards. We only select pixels which are adjacent to the previous row by moving
	// left and right cursors. For the right cursor, this means:
	// - If the pixel one down the previous right boundary is >halfA, we search if we need
	//  to move right (and stop when we don't find a new pixel > halfA or hit the end of the row)
	// -If not, we move left until we find again a pixel > halfA or we hit the left cursor
	// We then do same for the left cursor on same row until the cursors are same and
	// we are not anymore above halfA
	// And again, starting from the central row but upwards
	// For all the pixels above halfA, we compute their zero, first and second moments
	// wrt origin.
	// Centroid position wrt origin is then simply M1/M0
	// Inertia matrix about the centroid is assembled by correcting each M2 by
	// M0.xx, M0.yy and M0.xy
	// Then assuming that the trace is an ellipse, we use Singular Value Decomposition
	// to find the principal inertia moments and the rotation (angle)
	// From the inertia moments, we can then calculate main FWHM and roundness
	// As positions are 0 at the center of first pixel, 0.5 is added to centroids to
	// keep consistency with the convention.

	// There is still a shortfall with this init, which is that for saturated stars,
	// the amplitude is unknown and we will definitely init the optimization
	// with too small a value. But there's no quick and esay way to estimate A,
	// so we'll have to live with that for now.

	// computing moments
	gboolean goon = TRUE;
	int y = yc, cr = jj1, cl = jj2;

	// we will scan each row downwards to check pixel above halfA and adjacent to the previous row
	while (goon) {
		double y2 = SQR((double)y);
		if (cr < NbCols - 1 && gsl_matrix_get(m_tmp, y, cr) > halfA) {// we continue expanding right
			while (cr < NbCols - 1 && gsl_matrix_get(m_tmp, y, cr + 1) > halfA) {
				cr++;
			}
		} else {
			cr--;
			while (cr > 0 && cr > cl && gsl_matrix_get(m_tmp, y, cr - 1) <= halfA) {
				cr--;
			}
		}
		if (cl >= 0 && gsl_matrix_get(m_tmp, y, cl) > halfA) {// we continue expanding left
			while (cl > 1 && gsl_matrix_get(m_tmp, y, cl - 1) > halfA) {
				cl--;
			}
		} else {
			if (cl < cr) {
				cl++;
				while (cl < cr && gsl_matrix_get(m_tmp, y, cl + 1) <= halfA) {
					cl++;
				}
			}
		}
		if (cr == cl && gsl_matrix_get(m_tmp, y, cr) <= halfA) {
//			goon = FALSE; // Not needed, after the break goon is overwritten immediately
			break;
		}
		// we've found all the valid pixels, we can compute all the useful quantities
		for (int x = cl; x <= cr; x++) {
			double x2 = SQR((double)x);
			S += 1.;
			Sx += (double)x;
			Sy += (double)y;
			Ixx += y2;
			Iyy += x2;
			Ixy += (double)(x * y);
		}
		y++;
		if (y == NbRows)
			break;
	}

	goon = TRUE;
	y = yc - 1, cr = jj1, cl = jj2;
	if (y < 0)
		goon = FALSE;

	// we will scan each row upwards to check pixel above halfA and adjacent to the previous row
	while (goon) {
		double y2 = SQR((double)y);
		if (cr < NbCols - 1 && gsl_matrix_get(m_tmp, y, cr) > halfA) {// we continue expanding right
			while (cr < NbCols - 1 && gsl_matrix_get(m_tmp, y, cr + 1) > halfA) {
				cr++;
			}
		} else {
			cr--;
			while (cr > 0 && cr > cl && gsl_matrix_get(m_tmp, y, cr - 1) <= halfA) {
				cr--;
			}
		}
		if (cl >= 0 && gsl_matrix_get(m_tmp, y, cl) > halfA) {// we continue expanding left
			while (cl > 1 && gsl_matrix_get(m_tmp, y, cl - 1) > halfA) {
				cl--;
			}
		} else {
			if (cl < cr) {
				cl++;
				while (cl < cr && gsl_matrix_get(m_tmp, y, cl + 1) <= halfA) {
					cl++;
				}
			}
		}
		if (cr == cl && gsl_matrix_get(m_tmp, y, cr) <= halfA) {
			goon = FALSE;
			break;
		}
		// we've found all the valid pixels, we can compute all the useful quantities
		for (int x = cl; x <= cr; x++) {
			double x2 = SQR((double)x);
			S += 1.;
			Sx += (double)x;
			Sy += (double)y;
			Ixx += y2;
			Iyy += x2;
			Ixy += (double)(x * y);
		}
		y--;
		if (y < 0)
			break;
	}

	gsl_matrix_free(m_tmp);
	// centroid coordinates
	double x0 = Sx / S;
	double y0 = Sy / S;

	// inertia moments about the x/y axes at centroid
	Ixx -= S * y0 * y0;
	Iyy -= S * x0 * x0;
	Ixy -= S * x0 * y0;

	// SVD terms from https://lucidar.me/en/mathematics/singular-value-decomposition-of-a-2x2-matrix/
	double Su00 = Ixx * Ixx + Ixy * Ixy;
	double Su01 = (Ixx + Iyy) * Ixy;
	double Su11 = Iyy * Iyy + Ixy * Ixy;

	double ang = 90 + 0.5 * atan2( 2 * Su01, Su00 - Su11) * 180. / M_PI; //(eq 4)
	double SUsum = Su00 + Su11; // (eq7)
	double SUdif = sqrt(SQR(Su00 - Su11) + 4 * SQR(Su01)); // (eq7)

	// Inertia moment in the principal axes
	Iyy = sqrt((SUsum + SUdif) * 0.5); // (eq6)
	Ixx = sqrt((SUsum - SUdif) * 0.5); // (eq6)

	double r = sqrt(Ixx / Iyy);
	double FWHM = 2 * sqrt(S / M_PI / r);

	// vector init
	gsl_vector_set(MaxV, 1, x0 + 0.5); //x0
	gsl_vector_set(MaxV, 2, y0 + 0.5); //y0
	gsl_vector_set(MaxV, 3, FWHM); //FWHM x
	gsl_vector_set(MaxV, 4, FWHM * r); //FWHM y
	gsl_vector_set(MaxV, 5, ang); //angle
	return MaxV;
}

/* Basic magnitude computation. This is not really accurate, all pixels are
 * taken into account. But this fast function is used if the other one
 * failed and for star detection when magnitude is not needed.
 */
static double psf_get_mag(gsl_matrix* z, double B) {
	double intensity = 1.0;
	size_t NbRows = z->size1;
	size_t NbCols = z->size2;

	for (size_t i = 0; i < NbRows; i++) {
		for (size_t j = 0; j < NbCols; j++)
			intensity += gsl_matrix_get(z, i, j) - B;
	}
	return -2.5 * log10(intensity);
}

/* Gaussian */
static int psf_Gaussian_f_ang(const gsl_vector * x, void *PSF_data,
		gsl_vector * f) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t n = ((struct PSF_data *) PSF_data)->n;
	size_t i, j, k = 0;
	double *y = ((struct PSF_data *) PSF_data)->y;
	gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
	double B = gsl_vector_get(x, 0);
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = fabs(gsl_vector_get(x, 4));
	double r = 0.5 * (cos(gsl_vector_get(x, 5)) + 1.);
	double SY = SQR(r) * SX;
	double alpha = gsl_vector_get(x, 6);
	double tmpx, tmpy, tmpc, sumres = 0.;
	double ca = cos(alpha);
	double sa = sin(alpha);

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			if (mask[NbCols * i + j]) {
				tmpx = ca * (j + 0.5 - x0) - sa * (i + 0.5 - y0);
				tmpy = sa * (j + 0.5 - x0) + ca * (i + 0.5 - y0);
				tmpc = exp(-(SQR(tmpx) / SX + SQR(tmpy) / SY));
				gsl_vector_set(f, k,
						(B + A * tmpc - y[k]));
				sumres += (B + A * tmpc - y[k])
						* (B + A * tmpc - y[k]);
				k++;
			}
		}
	}
	((struct PSF_data *) PSF_data)->rmse = sqrt(sumres / n);
	return GSL_SUCCESS;
}

static int psf_Gaussian_df_ang(const gsl_vector * x, void *PSF_data,
		gsl_matrix * J) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t i, j, k = 0;
	const gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = fabs(gsl_vector_get(x, 4));
	double r = 0.5 * (cos(gsl_vector_get(x, 5)) + 1.);
	double SY = SQR(r) * SX;
	double alpha = gsl_vector_get(x, 6);
	double ca = cos(alpha);
	double sa = sin(alpha);
	double sc = sin(gsl_vector_get(x, 5));
	double tmpx, tmpy, tmpc, tmpd;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			if (mask[NbCols * i + j]) {
				tmpx = ca * (j + 0.5 - x0) - sa * (i + 0.5 - y0);
				tmpy = sa * (j + 0.5 - x0) + ca * (i + 0.5 - y0);
				tmpc = exp(-(SQR(tmpx) / SX + SQR(tmpy) / SY));
				gsl_matrix_set(J, k, 0, 1.); //dB
				gsl_matrix_set(J, k, 1, tmpc); //dA
				tmpd = 2 * A * tmpc * (  tmpx / SX * ca + tmpy / SY * sa); // dx0
				gsl_matrix_set(J, k, 2, tmpd);
				tmpd = 2 * A * tmpc * ( -tmpx / SX * sa + tmpy / SY * ca); // dy0
				gsl_matrix_set(J, k, 3, tmpd);
				tmpd = tmpc * A * (SQR(tmpx / SX) + SQR(tmpy / SX / r)); // dSX
				gsl_matrix_set(J, k, 4, tmpd);
				tmpd = -A * tmpc * sc * SQR(tmpy) / SY / r; // dfc
				gsl_matrix_set(J, k, 5, tmpd);
				tmpd = 2 * A * tmpc * tmpx * tmpy *(1. / SX - 1. / SY); // dalpha
				gsl_matrix_set(J, k, 6, tmpd);
				k++;
			}
		}
	}
	return GSL_SUCCESS;
}

/* Moffat */
static int psf_Moffat_f_ang(const gsl_vector * x, void *PSF_data, gsl_vector * f) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t n = ((struct PSF_data *) PSF_data)->n;
	size_t i, j, k = 0;
	const double *y = ((struct PSF_data *) PSF_data)->y;
	const gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
	double B = gsl_vector_get(x, 0);
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = fabs(gsl_vector_get(x, 4)); // Ro_x ^ 2
	double r = 0.5 * (cos(gsl_vector_get(x, 5)) + 1.);
	double SY = SQR(r) * SX; // Ro_y ^ 2
	double alpha = gsl_vector_get(x, 6);
	double ca = cos(alpha);
	double sa = sin(alpha);
	double beta = MOFFAT_BETA_UBOUND * 0.5 * (cos(gsl_vector_get(x, 7)) + 1.);
	// // bounding beta to MOFFAT_BETA_UBOUND if free, otherwise setting to PSF_data->beta
	// double beta = (((struct PSF_data *) PSF_data)->betafree ?
	// 		min(gsl_vector_get(x, 6), MOFFAT_BETA_UBOUND) :
	// 		((struct PSF_data *) PSF_data)->beta);
	double tmpx, tmpy, tmpa, sumres = 0.;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			if (mask[NbCols * i + j]) {
				tmpx = ca * (j + 0.5 - x0) - sa * (i + 0.5 - y0);
				tmpy = sa * (j + 0.5 - x0) + ca * (i + 0.5 - y0);
				tmpa = pow(1 + SQR(tmpx) / SX + SQR(tmpy) / SY, -beta);
				gsl_vector_set(f, k,
						(B + A * tmpa - y[k]));
				sumres += (B + A * tmpa - y[k])
						* (B + A * tmpa - y[k]);
				k++;
			}
		}
	}
	((struct PSF_data *) PSF_data)->rmse = sqrt(sumres / n);
	return GSL_SUCCESS;
}

static int psf_Moffat_df_ang(const gsl_vector * x, void *PSF_data, gsl_matrix * J) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t i, j, k = 0;
	const gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = fabs(gsl_vector_get(x, 4)); // Ro_x ^ 2
	double r = 0.5 * (cos(gsl_vector_get(x, 5)) + 1.);
	double SY = SQR(r) * SX; // Ro_y ^ 2
	double sc = sin(gsl_vector_get(x, 5));
	double alpha = gsl_vector_get(x, 6);
	double ca = cos(alpha);
	double sa = sin(alpha);
	double beta = MOFFAT_BETA_UBOUND * 0.5 * (cos(gsl_vector_get(x, 7)) + 1.);
	double sfbeta = sin(gsl_vector_get(x, 7));
	// // bounding beta to MOFFAT_BETA_UBOUND if free, otherwise setting to PSF_data->beta
	// double beta = (((struct PSF_data *) PSF_data)->betafree ?
	// 		min(gsl_vector_get(x, 6), MOFFAT_BETA_UBOUND) :
	// 		((struct PSF_data *) PSF_data)->beta);
	double tmpx, tmpy, tmpa, tmpb, tmpc, tmpd;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			if (mask[NbCols * i + j]) {
				tmpx = ca * (j + 0.5 - x0) - sa * (i + 0.5 - y0);
				tmpy = sa * (j + 0.5 - x0) + ca * (i + 0.5 - y0);
				tmpa = 1 + SQR(tmpx) / SX + SQR(tmpy) / SY;
				tmpb = pow(tmpa, -beta);
				tmpc = A * beta * pow(tmpa, -beta - 1.0);
				gsl_matrix_set(J, k, 0, 1.); // dB
				gsl_matrix_set(J, k, 1, tmpb); // dA
				tmpd = 2 * tmpc * (  tmpx / SX * ca + tmpy / SY * sa); // dx0
				gsl_matrix_set(J, k, 2, tmpd);
				tmpd = 2 * tmpc * ( -tmpx / SX * sa + tmpy / SY * ca); // dy0
				gsl_matrix_set(J, k, 3, tmpd);
				tmpd = tmpc * (SQR(tmpx / SX) + SQR(tmpy / SX / r));
				gsl_matrix_set(J, k, 4, tmpd); // dSX
				tmpd = -tmpc * sc * SQR(tmpy) / SY / r;
				gsl_matrix_set(J, k, 5, tmpd); // dfc
				tmpd = 2 * tmpc * tmpx * tmpy *(1. / SX - 1. / SY); // dalpha
				gsl_matrix_set(J, k, 6, tmpd);
				tmpd = 0.5 * A * MOFFAT_BETA_UBOUND * sfbeta * log(tmpa) * tmpb; // dfbeta
				gsl_matrix_set(J, k, 7, tmpd);
				k++;
			}
		}
	}
	return GSL_SUCCESS;
}

#if DEBUG_PSF
struct callback_params{
	starprofile profile;		// star profile type
};

static void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
	gsl_vector *f = gsl_multifit_nlinear_residual(w);
	gsl_vector *x = gsl_multifit_nlinear_position(w);
	struct callback_params *c_params = (struct callback_params *) params;
	gboolean ismoffat = c_params->profile != PSF_GAUSSIAN;

	if (ismoffat) {
		if (iter == 0) {
			fprintf(stdout, "i\tB\tA\tx0\ty0\tFWHM\tr\tangle\tbeta\t|f|\n");
		}
		fprintf(stdout, "%2zu\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",
						iter,
						gsl_vector_get(x, 0),
						gsl_vector_get(x, 1),
						gsl_vector_get(x, 2),
						gsl_vector_get(x, 3),
						FWHM_from_S(fabs(gsl_vector_get(x, 4)), 0.5 * MOFFAT_BETA_UBOUND * (cos(gsl_vector_get(x, 7)) + 1.), PSF_MOFFAT_BFREE), // FWHM
						0.5 * (cos(gsl_vector_get(x, 5)) + 1.), // roundness
						gsl_vector_get(x, 6) * 180. / M_PI,
						0.5 * MOFFAT_BETA_UBOUND * (cos(gsl_vector_get(x, 7)) + 1.), // beta
						gsl_blas_dnrm2(f));
	} else {
		if (iter == 0) {
			fprintf(stdout, "i\tB\tA\tx0\ty0\tFWHM\tr\tangle\t|f|\n");
		}
		fprintf(stdout, "%2zu\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",
						iter,
						gsl_vector_get(x, 0),
						gsl_vector_get(x, 1),
						gsl_vector_get(x, 2),
						gsl_vector_get(x, 3),
						FWHM_from_S(fabs(gsl_vector_get(x, 4)), 0., PSF_GAUSSIAN), // FWHM
						0.5 * (cos(gsl_vector_get(x, 5)) + 1.), // roundness
						gsl_vector_get(x, 6) * 180. / M_PI,
						gsl_blas_dnrm2(f));

	}
}
#endif

/* The function returns the fitted parameters with angle. However it returns
 * NULL if the number of parameters is => to the pixel number.
 */
static psf_star *psf_minimiz_angle(gsl_matrix* z, double background, double sat, int convergence, gboolean from_peaker, gboolean for_photometry, struct phot_config *phot_set, gboolean verbose, starprofile profile, psf_error *error) {
	size_t i, j, k = 0;
	size_t NbRows = z->size1; //characteristics of the selection : height and width
	size_t NbCols = z->size2;
	const size_t p = (profile == PSF_GAUSSIAN) ? 7 : 8;	// Number of parameters fitted
	int status;
	gboolean *mask = NULL;
	gsl_vector *MaxV = NULL;
	gsl_matrix *covar = NULL;
	psf_star *psf = NULL;
	double *y = NULL;
	int max_iter;
	gsl_multifit_nlinear_workspace *work = NULL;

	if (error) *error = PSF_NO_ERR;
	// computing the mask to discard clipped values
	mask = malloc(NbRows * NbCols * sizeof(gboolean));
	if (!mask) {
		PRINT_ALLOC_ERR;
		if (error) *error = PSF_ERR_ALLOC;
		goto free_and_exit;
	}
	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			mask[NbCols * i + j] = gsl_matrix_get(z, i, j) < sat;
			if (mask[NbCols * i + j]) k++;
		}
	}

	const size_t n = k;
	if (n <= p) { // could happen if star is mostly saturated (hand-selection case)
		if (error) *error = PSF_ERR_WINDOW_TOO_SMALL;
		goto free_and_exit;
	}

	psf = new_psf_star();
	covar = gsl_matrix_alloc(p, p);
	y = malloc(n * sizeof(double));
	if (!psf || !covar || !y) {
		PRINT_ALLOC_ERR;
		if (error) *error = PSF_ERR_ALLOC;
		if (psf) free_psf(psf);
		psf = NULL;
		goto free_and_exit;
	}

	max_iter = MAX_ITER_ANGLE * ((k < NbRows * NbCols) ? 3 : 1) * convergence * ((profile == PSF_GAUSSIAN) ? 1 : 2);

	MaxV = psf_init_data(z, background, from_peaker);
	if (!MaxV) {
		PRINT_ALLOC_ERR;
		if (error) *error = PSF_ERR_ALLOC;
		free_psf(psf);
		psf = NULL;
		goto free_and_exit;
	}

	double beta = (profile == PSF_GAUSSIAN) ? -1. : 2; // TODO: to be changed if we implement PSF_MOFFAT_BFIXED
	double fbeta = (profile == PSF_GAUSSIAN) ? 0. : acos(2. * beta / MOFFAT_BETA_UBOUND - 1.);

	struct PSF_data d = { n, y, NbRows, NbCols, 0. , mask };
	double FWHM = gsl_vector_get(MaxV, 3);
	double roundness = gsl_vector_get(MaxV, 4) / gsl_vector_get(MaxV, 3);
	double a_init = gsl_vector_get(MaxV, 5) * M_PI / 180.; // angle in radians
	// if roundness is 1., we decrease it a bit so as not to be stuck on the boundary
	// as it is messes up the initial gradient calcs
	if (roundness == 1.) {
		roundness = 0.9;
		a_init = 0.;
	}
	double fr = acos(2. * roundness - 1.); // r = 0.5 *(cos(fc)+1) to bound it between 0 and 1
	double x0_init = gsl_vector_get(MaxV, 1);
	double y0_init = gsl_vector_get(MaxV, 2);
	double x_init[8] = { background, // B
						gsl_vector_get(MaxV, 0), // A
						x0_init, // x0
						y0_init, // y0
						S_from_FWHM(FWHM, beta, profile), // SX
						fr, //
						a_init, // angle
						fbeta}; // beta = betamax * 0.5 * (cos(fbeta) + 1)
	gsl_vector_view x = gsl_vector_view_array(x_init, p);

	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	fdf_params.trs = gsl_multifit_nlinear_trs_lm; // levenberg-marquardt
	gsl_multifit_nlinear_fdf fdf;
	if (profile == PSF_GAUSSIAN) {
		fdf.f = &psf_Gaussian_f_ang;
		fdf.df = &psf_Gaussian_df_ang;
	} else {
		fdf.f = &psf_Moffat_f_ang;
		fdf.df = &psf_Moffat_df_ang;
	}
	fdf.fvv = NULL;
	fdf.n = n;
	fdf.p = p;
	fdf.params = &d;

	k = 0;
	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			if (mask[NbCols * i + j]) {
				y[k] = gsl_matrix_get(z, i, j);
				k++;
			}
		}
	}
	g_assert(k == n);

	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	work = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);
	int info;

	/* initialize solver */
	gsl_multifit_nlinear_init(&x.vector, &fdf, work);

	/* iterate until convergence */
#if DEBUG_PSF
	struct callback_params c_params = { profile };
	status = gsl_multifit_nlinear_driver(max_iter, XTOL, GTOL, FTOL,
	callback, &c_params, &info, work);
#else
	status = gsl_multifit_nlinear_driver(max_iter, XTOL, GTOL, FTOL,
	NULL, NULL, &info, work);
#endif

	if (status != GSL_SUCCESS) {
		if (error) *error = PSF_ERR_DIVERGED;
	}
#if DEBUG_PSF
	siril_debug_print("Successful criterion#:%d\n",info);
#endif

	/* computing the covariance to estimate the errors*/
	gsl_matrix * J;

	J = gsl_multifit_nlinear_jac(work);
	gsl_multifit_nlinear_covar (J, 0.0, covar);

#define FIT(i) gsl_vector_get(work->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))	//for now, errors are not displayed

	/*Output structure with parameters fitted */
	psf->profile = profile;
	psf->B = FIT(0);
	psf->A = FIT(1);
	psf->x0 = FIT(2);
	psf->y0 = FIT(3);
	psf->beta = (profile == PSF_GAUSSIAN) ? -1.0 : 0.5 * MOFFAT_BETA_UBOUND * (cos(FIT(7)) + 1.);
	psf->sx = (profile == PSF_GAUSSIAN) ? sqrt(fabs(FIT(4)) * 0.5) : sqrt(fabs(FIT(4))); // Gaussian: sigma, Moffat: Ro
	double r = 0.5 * (cos(FIT(5)) + 1.);
	psf->sy = psf->sx * r;
	psf->fwhmx = FWHM_from_s(psf->sx, psf->beta, profile);	//Set the real FWHMx with regards to the sx parameter
	psf->fwhmy = FWHM_from_s(psf->sy, psf->beta, profile);	//Set the real FWHMy with regards to the Sy parameter
	psf->angle = -FIT(6) * 180.0 / M_PI;

	/* In some cases convergence give crazy values
	 * very high. Here we add a sanity check to avoid
	 * pseudo infinite loop with the while.
	 */
	if (fabs(psf->angle) > 10000) {
		free_psf(psf);
		psf = NULL;
		goto free_and_exit;
	}
	/* The angle must be => -90 and <= 90
	 * Otherwise, the solution may be degenerate
	 * and produce an angle > 90. So we're
	 * looking for the solution between
	 * the interval we want */
	while (fabs(psf->angle) > 90.0) {
		if (psf->angle > 0.0)
			psf->angle -= 180.0;
		else
			psf->angle += 180.0;
	}
	//Units
	psf->units = "px";
	// Photometry
	if (for_photometry)
		psf->phot = getPhotometryData(z, psf, phot_set, verbose, error);
	else {
		psf->phot = NULL;
		psf->phot_is_valid = FALSE;
	}
	// Magnitude
	if (psf->phot && psf->phot->valid) {
		psf->mag = psf->phot->mag;
		psf->s_mag = psf->phot->s_mag;
		psf->SNR = psf->phot->SNR;
		psf->phot_is_valid = psf->phot->valid;
	} else {
		psf->mag = psf_get_mag(z, psf->B);
		psf->s_mag = 9.999;
		psf->SNR = 0;
		psf->phot_is_valid = FALSE;
	}
	//RMSE
	psf->rmse = d.rmse;
	// absolute uncertainties
	// TODO: this will need to be revisited if of use as we are using intermediate variables
	psf->B_err = ERR(0) / FIT(0);
	psf->A_err = ERR(1) / FIT(1);
	psf->x_err = ERR(2) / FIT(2);
	psf->y_err = ERR(3) / FIT(3);
	psf->sx_err = ERR(4) / FIT(4);
	psf->sy_err = ERR(5) / FIT(5);
	psf->ang_err = ERR(6) / FIT(6);
	psf->beta_err = (profile == PSF_GAUSSIAN) ? 0. : ERR(7) / FIT(7);

	//we free the memory
free_and_exit:
	if (y) free(y);
	if (mask) free(mask);
	if (MaxV) gsl_vector_free(MaxV);
	if(work) gsl_multifit_nlinear_free(work);
	if (covar) gsl_matrix_free(covar);
	return psf;
}
/******************************************************************************/

/* Returns the largest FWHM in pixels
 * The optional output parameter roundness is the ratio between the two axis FWHM */
double psf_get_fwhm(fits *fit, int layer, rectangle *selection, double *roundness) {
	psf_star *result = psf_get_minimisation(fit, layer, selection, FALSE, NULL, TRUE, com.pref.starfinder_conf.profile, NULL);
	if (result == NULL) {
		*roundness = 0.0;
		return 0.0;
	}
	double retval;
	retval = result->fwhmx;
	if (roundness)
		*roundness = result->fwhmy / result->fwhmx;
	free_psf(result);
	return retval;
}

/* Computes the FWHM on data in the selection rectangle of image fit.
 * Selection rectangle is passed as third argument.
 * Return value is a structure, type psf_star, that has to be freed after use.
 * verbose is used in photometry only, to inform that inner is too small for example
 */
psf_star *psf_get_minimisation(fits *fit, int layer, rectangle *area,
		gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error) {
	int stridefrom, i, j;
	psf_star *result;
	if (error) *error = PSF_NO_ERR;
	double bg = background(fit, layer, area, SINGLE_THREADED);
	double sat;
	if (bg == -1.0) {
		if (error) *error = PSF_ERR_INVALID_IMAGE;
		return NULL;
	}

	// fprintf(stdout, "background: %g\n", bg);
	gsl_matrix *z = gsl_matrix_alloc(area->h, area->w);
	stridefrom = fit->rx - area->w;

	// create the matrix with values from the selected rectangle
	// area coordinates are in display coordinates, the matrix is read top-down but in
	// FITS coordinates
	if (fit->type == DATA_USHORT) {
		WORD *from = fit->pdata[layer] +
			(fit->ry - area->y - area->h) * fit->rx + area->x;
		sat = (fit->orig_bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;

		for (i = 0; i < area->h; i++) {
			for (j = 0; j < area->w; j++) {
				gsl_matrix_set(z, i, j, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *from = fit->fpdata[layer] +
			(fit->ry - area->y - area->h) * fit->rx + area->x;
		sat = 1.;

		for (i = 0; i < area->h; i++) {
			for (j = 0; j < area->w; j++) {
				gsl_matrix_set(z, i, j, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else {
		gsl_matrix_free(z);
		if (error) *error = PSF_ERR_UNSUPPORTED;
		return NULL;
	}

	result = psf_global_minimisation(z, bg, sat, com.pref.starfinder_conf.convergence, FALSE, for_photometry, phot_set, verbose, profile, error);

	if (result) {
		fwhm_to_arcsec_if_needed(fit, result);
		result->layer = layer;
	}
	gsl_matrix_free(z);
	return result;
}

/* This function is the global minimisation. Every call to the minimisation
 * must come over here.
 * If fit_angle, it will check if the difference between Sx and Sy is larger
 * than or equal to 0.01 pixel. In this case, Dynamic PSF fits additional angle
 * parameter which is the rotation angle of the X axis with respect to the
 * centroid coordinates: so, by design we set Sx>Sy.
 * If the difference is smaller OR if fit_Angle is equal to FALSE (in the case
 * of the star_finder algorithm), no angle parameter is fitted.
 * The function returns NULL if values look bizarre.
 * The photometry config is only required if for_photometry is TRUE.
 * Error can be reported if error is provided, giving a reason for the failure
 * of the minimisation.
 */
psf_star *psf_global_minimisation(gsl_matrix* z, double bg, double sat, int convergence,
		gboolean from_peaker, gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error) {
	if (error) *error = PSF_NO_ERR;
//	gboolean photometry_computed = FALSE; // This is never used except in the dead code commented out later

	psf_star *psf = NULL;
	if (!(psf = psf_minimiz_angle(z, bg, sat, convergence, from_peaker, for_photometry, phot_set, verbose, profile, error))) {
		return NULL;
	}
//	photometry_computed = TRUE;

	/* We quickly test the result. If it is bad we return NULL */
	if (!isfinite(psf->fwhmx) || !isfinite(psf->fwhmy) ||
			psf->fwhmx <= 0.0 || psf->fwhmy <= 0.0) {
		free_psf(psf);
		if (error && *error == PSF_NO_ERR)
			*error = PSF_ERR_DIVERGED;
		return NULL;
	}

/* This code is logically dead. Commenting out prior to removal.
 *	// Photometry
	if (for_photometry && !photometry_computed &&
			(!error || *error == PSF_NO_ERR || *error == PSF_ERR_DIVERGED)) {
		psf->phot = getPhotometryData(z, psf, phot_set, verbose, error);
		if (psf->phot) {
			psf->mag = psf->phot->mag;
			psf->s_mag = psf->phot->s_mag;
			psf->SNR = psf->phot->SNR;
			psf->phot_is_valid = psf->phot->valid;
		}
		else {
			psf->phot_is_valid = FALSE;
			psf->s_mag = 9.999;
			psf->SNR = 0;
		}
	}
*/
	return psf;
}

static gchar *build_wcs_url(gchar *ra, gchar *dec) {
	if (!has_wcs(&gfit)) return NULL;

	double resolution = get_wcs_image_resolution(&gfit);

	gchar *tol = g_strdup_printf("%lf", resolution * 3600 * 15);

	GString *url = g_string_new("https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=");
	url = g_string_append(url, ra);
	url = g_string_append(url, dec);
	url = g_string_append(url, "&Radius=");
	url = g_string_append(url, tol);
	url = g_string_append(url, "&Radius.unit=arcsec");
	url = g_string_append(url, "#lab_basic");

	gchar *simbad_url = g_string_free(url, FALSE);
	gchar *cleaned_url = url_cleanup(simbad_url);

	g_free(tol);
	g_free(simbad_url);

	return cleaned_url;
}

static const char *SNR_quality(double SNR) {
	if (SNR > 40.0) return _("Excellent");
	if (SNR > 25.0) return _("Good");
	if (SNR > 15.0) return _("Fair");
	if (SNR > 10.0) return _("Poor");
	if (SNR > 0.0) return _("Bad");
	else return _("N/A");
}

gchar *format_psf_result(psf_star *result, const rectangle *area, fits *fit, gchar **url) {
	gchar *msg, *coordinates;
	char buffer2[50];
	const char *str;
	if (com.magOffset > 0.0)
		str = _("true reduced");
	else
		str = _("relative");

	// coordinates of the star in the displayed image
	double xpos = result->x0 + area->x;
	double ypos = area->y + area->h - result->y0;

	if (has_wcs(&gfit)) {
		// coordinates of the star in FITS/WCS coordinates
		double fx, fy;
		display_to_siril(xpos, ypos, &fx, &fy, gfit.ry);

		double ra, dec;
		pix2wcs(&gfit, fx, fy, &ra, &dec);
		SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(ra, dec);
		if (world_cs) {
			gchar *strra, *strdec;
			if (url) {
				strra = siril_world_cs_alpha_format(world_cs, "%02d %02d %.3lf");
				strdec = siril_world_cs_delta_format(world_cs, "%c%02d %02d %.3lf");
				*url = build_wcs_url(strra, strdec);
				// TODO: change with vizier
				// TODO: use box size as radius
				g_free(strra);
				g_free(strdec);
			}

			if (com.pref.gui.show_deciasec) {
				strra = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%04.1lfs");
				strdec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%04.1lf\"");
			} else {
				strra = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
				strdec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
			}

			coordinates = g_strdup_printf("x0=%.2fpx\t%s J2000\n\t\ty0=%.2fpx\t%s J2000", xpos, strra, ypos, strdec);

			g_free(strra);
			g_free(strdec);
			siril_world_cs_unref(world_cs);
		} else {
			coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", xpos, ypos);
		}
	} else {
		coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", xpos, ypos);
	}

	double fwhmx, fwhmy;
	char *unts;
	get_fwhm_as_arcsec_if_possible(result, &fwhmx, &fwhmy, &unts);
	const gchar *chan = isrgb(fit) ? channel_number_to_name(result->layer) : _("monochrome");
	if (result->beta > 0.0) {
		g_snprintf(buffer2, 50, _(", beta=%0.1f, %s channel"), result->beta, chan);
	}
	else {
		g_snprintf(buffer2, 50, _(", %s channel"), chan);
	}
	msg = g_strdup_printf(_("PSF fit Result (%s%s):\n\n"
				"Centroid Coordinates:\n\t\t%s\n\n"
				"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\t\tr=%.2f\n"
				"Angle:\n\t\t%0.2fdeg\n\n"
				"Background Value:\n\t\tB=%.6f\n\n"
				"Maximal Intensity:\n\t\tA=%.6f\n\n"
				"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
				"Signal-to-noise ratio:\n\t\tSNR=%.1fdB (%s)\n\n"
				"RMSE:\n\t\tRMSE=%.3e"),
			(result->profile == PSF_GAUSSIAN) ? _("Gaussian") : _("Moffat"), buffer2,
			coordinates, fwhmx, unts, fwhmy, unts, fwhmy / fwhmx,
			result->angle, result->B, result->A, str,
			result->mag + com.magOffset, result->s_mag, result->SNR,
			SNR_quality(result->SNR), result->rmse);
	g_free(coordinates);
	return msg;
}

/* If the pixel pitch and the focal length are known and filled in the
 * setting box, we convert FWHM in pixel to arcsec by multiplying
 * the FWHM value with the sampling value */
void fwhm_to_arcsec_if_needed(fits* fit, psf_star *result) {

	if (!result) return;
	if (fit->keywords.focal_length <= 0.0 || fit->keywords.pixel_size_x <= 0.f
			|| fit->keywords.pixel_size_y <= 0.f || fit->keywords.binning_x <= 0
			|| fit->keywords.binning_y <= 0) {
		result->fwhmx_arcsec = -1.0;
		result->fwhmy_arcsec = -1.0;
		return;
	}

	double bin_X, bin_Y;

	bin_X = com.pref.binning_update ? (double) fit->keywords.binning_x : 1.0;
	bin_Y = com.pref.binning_update ? (double) fit->keywords.binning_y : 1.0;

	result->fwhmx_arcsec = result->fwhmx * (radian_conversion * (double)fit->keywords.pixel_size_x / fit->keywords.focal_length) * bin_X;
	result->fwhmy_arcsec = result->fwhmy * (radian_conversion * (double)fit->keywords.pixel_size_y / fit->keywords.focal_length) * bin_Y;
	result->units = "\"";
}

// returns boolean if it was possible (true if arcsec)
gboolean get_fwhm_as_arcsec_if_possible(psf_star *star, double *fwhmx, double *fwhmy, char **unit) {
	if (!strcmp(star->units, "px")) {
		*fwhmx = star->fwhmx;
		*fwhmy = star->fwhmy;
		*unit = star->units;
		return FALSE;
	}
	if (star->fwhmx_arcsec <= 0.0) {
		fprintf(stderr, "FWHM wrongly stored as arcsec\n");
		star->units = "px";
		return get_fwhm_as_arcsec_if_possible(star, fwhmx, fwhmy, unit);
	}
	*fwhmx = star->fwhmx_arcsec;
	*fwhmy = star->fwhmy_arcsec;
	*unit = star->units;
	return TRUE;
}

gboolean convert_single_fwhm_to_arcsec_if_possible(double fwhm, double bin, double px_size, double flength, double *result) {
	double arcsec = fwhm * (radian_conversion * px_size / flength) * bin;
	if (arcsec <= 0.0 || isnan(arcsec) || !isfinite(arcsec)) {
		*result = 0;
		return FALSE;
	}
	*result = arcsec;
	return TRUE;
}

void psf_star_init(psf_star *s) {
	memset(s, 0, sizeof(psf_star));
	s->units = "px";
}

psf_star *new_psf_star() {
	psf_star *star = calloc(1, sizeof(psf_star));
	star->phot = NULL;
	star->units = "px";
	return star;
}

psf_star *duplicate_psf(psf_star *psf) {
	if (!psf)
		return NULL;
	psf_star *new_psf = new_psf_star();
	memcpy(new_psf, psf, sizeof(psf_star));
	if (psf->phot) {
		new_psf->phot = malloc(sizeof(photometry));
		memcpy(new_psf->phot, psf->phot, sizeof(photometry));
	} else {
		new_psf->phot = NULL;
	}
	if (psf->star_name)
		new_psf->star_name = g_strdup(psf->star_name);
	return new_psf;
}

void free_psf(psf_star *psf) {
	if (psf->phot) free(psf->phot);
	if (psf->star_name) g_free(psf->star_name);
	free(psf);
}
