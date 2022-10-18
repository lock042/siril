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
#define MOFFAT_BETA_UBOUND 10. // Max allowable value for Moffat beta
#define EPSILON            0.001
#define XTOL 1e-3
#define GTOL 1e-3
#define FTOL 1e-3

#define DEBUG_PSF 0 // flag to show progress of fitting process - may flood output if numerous stars

const double radian_conversion = ((3600.0 * 180.0) / M_PI) / 1.0E3;

static gsl_matrix *removeHotPixels(gsl_matrix *in) {
	size_t width = in->size2;
	size_t height = in->size1;
	size_t x, y;
	gsl_matrix *out = gsl_matrix_alloc (in->size1, in->size2);
	if (!out) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	gsl_matrix_memcpy (out, in);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			double a = get_median_gsl(in, x, y, width, height, 1, FALSE, FALSE);
			gsl_matrix_set(out, y, x, a);
		}
	}
	return out;
}

static double S_from_FWHM(double FWHM, double beta, starprofile profile) {
	if (profile == GAUSSIAN)
		return SQR(FWHM) * INV_4_LOG2;
	return SQR(FWHM) * 0.25 / (pow(2., 1. / beta) - 1.);
}

#if DEBUG_PSF
static double FWHM_from_S(double S, double beta, starprofile profile) {
	if (profile == GAUSSIAN)
		return 2. * sqrt(S * log(2.));
	return 2. * sqrt(S * (pow(2., 1./ beta) - 1.));
}
#endif

// static double s_from_FWHM(double FWHM, double beta, starprofile profile) {
// 	if (profile == GAUSSIAN)
// 		return FWHM / _2_SQRT_2_LOG2;
// 	return FWHM * 0.5 /  sqrt(pow(2., 1./ beta) - 1.);
// }

static double FWHM_from_s(double s, double beta, starprofile profile) {
	if (profile == GAUSSIAN)
		return s * _2_SQRT_2_LOG2;
	return 2. * s * sqrt(pow(2., 1./ beta) - 1.);
}

/* Compute initial values for the algorithm from data in the pixel value matrix */
static gsl_vector* psf_init_data(gsl_matrix* z, double bg) {
	gsl_vector *MaxV = gsl_vector_alloc(5);
	if (!MaxV) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	double max, halfA;
	size_t NbRows = z->size1; //y
	size_t NbCols = z->size2; //x
	size_t i, j;

	/* find maximum */
	/* first we remove hot pixels in the matrix */
	gsl_matrix *m_tmp = removeHotPixels(z);
	if (!m_tmp) return NULL;
	max = gsl_matrix_max(m_tmp);
	gsl_matrix_max_index(m_tmp, &i, &j); // i=y , j=x
	gsl_matrix_free(m_tmp);
	gsl_vector_set(MaxV, 1, j); // x0
	gsl_vector_set(MaxV, 2, i); // y0
	halfA = (max - bg) * 0.5; // half amplitude
	gsl_vector_set(MaxV, 0, max - bg); // A

	size_t ii1 = (size_t) gsl_vector_get(MaxV, 2);
	size_t ii2 = (size_t) gsl_vector_get(MaxV, 2);
	size_t jj1 = (size_t) gsl_vector_get(MaxV, 1);
	size_t jj2 = (size_t) gsl_vector_get(MaxV, 1);
	size_t perm1 = (size_t) gsl_vector_get(MaxV, 2);
	size_t perm2 = (size_t) gsl_vector_get(MaxV, 1);

	while (gsl_matrix_get(z, ii1, perm2) - bg > halfA && (ii1 < NbRows - 1.0)) {
		ii1++;
	}
	while (gsl_matrix_get(z, ii2, perm2) - bg > halfA && (ii2 > 0)) {
		ii2--;
	}

	while (gsl_matrix_get(z, perm1, jj1) - bg > halfA && (jj1 < NbCols - 1)) {
		jj1++;
	}
	while (gsl_matrix_get(z, perm1, jj2) - bg > halfA && (jj2 > 0)) {
		jj2--;
	}
	gsl_vector_set(MaxV, 1, (jj1 + jj2 + 1) / 2.0); //x0
	gsl_vector_set(MaxV, 2, (ii1 + ii2 + 1) / 2.0); //y0
	gsl_vector_set(MaxV, 3, jj1 - jj2); //FWHM x
	gsl_vector_set(MaxV, 4, ii1 - ii2); //FWHM y

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
	gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
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
	double *y = ((struct PSF_data *) PSF_data)->y;
	gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
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
	gboolean *mask = ((struct PSF_data *) PSF_data)->mask;
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
	gboolean ismoffat = c_params->profile != GAUSSIAN;

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
						FWHM_from_S(fabs(gsl_vector_get(x, 4)), 0.5 * MOFFAT_BETA_UBOUND * (cos(gsl_vector_get(x, 7)) + 1.), MOFFAT_BFREE), // FWHM
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
						FWHM_from_S(fabs(gsl_vector_get(x, 4)), 0., GAUSSIAN), // FWHM
						0.5 * (cos(gsl_vector_get(x, 5)) + 1.), // roundness
						gsl_vector_get(x, 6) * 180. / M_PI,
						gsl_blas_dnrm2(f));

	}
}
#endif

/* The function returns the fitted parameters with angle. However it returns
 * NULL if the number of parameters is => to the pixel number.
 */
static psf_star *psf_minimiz_angle(gsl_matrix* z, double background, double sat, int convergence, gboolean for_photometry, struct phot_config *phot_set, gboolean verbose, starprofile profile, psf_error *error) {
	size_t i, j, k = 0;
	size_t NbRows = z->size1; //characteristics of the selection : height and width
	size_t NbCols = z->size2;
	const size_t p = (profile == GAUSSIAN) ? 7 : 8;	// Number of parameters fitted
	int status;
	gboolean *mask = NULL;
	gsl_vector *MaxV = NULL;
	gsl_matrix *covar = NULL;
	psf_star *psf;
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

	max_iter = MAX_ITER_ANGLE * ((k < NbRows * NbCols) ? 3 : 1) * convergence;

	MaxV = psf_init_data(z, background);
	if (!MaxV) {
		PRINT_ALLOC_ERR;
		if (error) *error = PSF_ERR_ALLOC;
		if (psf) free_psf(psf);
		psf = NULL;
		goto free_and_exit;
	}

	double beta = (profile == GAUSSIAN) ? -1. : 2; // TODO: to be changed if we implement MOFFAT_BFIXED
	double fbeta = acos(2. * beta / MOFFAT_BETA_UBOUND - 1.);

	struct PSF_data d = { n, y, NbRows, NbCols, 0. , mask };
	double FWHM = max(gsl_vector_get(MaxV, 3), gsl_vector_get(MaxV, 4));
	double a_init = (gsl_vector_get(MaxV, 3) > gsl_vector_get(MaxV, 4)) ? 0. : M_PI / 2;
	double x_init[8] = { background, // B
						gsl_vector_get(MaxV, 0), // A
						gsl_vector_get(MaxV, 1), // x0
						gsl_vector_get(MaxV, 2), // y0
						S_from_FWHM(FWHM, beta, profile), // SX
						M_PI / 2., // r = 0.5*(cos(fc)+1) to bound it between 0 and 1 - we init at r = 0.5 to be able to start in a place where the direction of variation is well defined
						a_init, // angle
						fbeta}; // beta = betamax * 0.5 * (cos(fbeta) + 1)
	gsl_vector_view x = gsl_vector_view_array(x_init, p);

	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	fdf_params.trs = gsl_multifit_nlinear_trs_lm; // levenberg-marquardt
	gsl_multifit_nlinear_fdf fdf;
	if (profile == GAUSSIAN) {
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
	covar = gsl_matrix_alloc (p, p);

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
	psf->beta = (profile == GAUSSIAN) ? -1.0 : 0.5 * MOFFAT_BETA_UBOUND * (cos(FIT(7)) + 1.);
	psf->sx = (profile == GAUSSIAN) ? sqrt(fabs(FIT(4)) * 0.5) : sqrt(fabs(FIT(4))); // Gaussian: sigma, Moffat: Ro
	double r = 0.5 * (cos(FIT(5)) + 1.);
	psf->sy = psf->sx * r;
	psf->fwhmx = FWHM_from_s(psf->sx, psf->beta, profile);	//Set the real FWHMx with regards to the sx parameter
	psf->fwhmy = FWHM_from_s(psf->sy, psf->beta, profile);	//Set the real FWHMy with regards to the Sy parameter
	psf->angle = -FIT(6) * 180.0 / M_PI;

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
	// TODO: this will need to be revisited if of use as we are using interemdiate variables
	psf->B_err = ERR(0) / FIT(0);
	psf->A_err = ERR(1) / FIT(1);
	psf->x_err = ERR(2) / FIT(2);
	psf->y_err = ERR(3) / FIT(3);
	psf->sx_err = ERR(4) / FIT(4);
	psf->sy_err = ERR(5) / FIT(5);
	psf->ang_err = ERR(6) / FIT(6);
	psf->beta_err = (profile == GAUSSIAN) ? 0. : ERR(7) / FIT(7);

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
	psf_star *result = psf_get_minimisation(fit, layer, selection, FALSE, FALSE, NULL, TRUE, GAUSSIAN, NULL);
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
psf_star *psf_get_minimisation(fits *fit, int layer, rectangle *area, gboolean fit_angle,
		gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error) {
	int stridefrom, i, j;
	psf_star *result;
	if (error) *error = PSF_NO_ERR;
	double bg = background(fit, layer, area, SINGLE_THREADED);
	double sat;
	if (bg == 1.0) {
		if (error) *error = PSF_ERR_INVALID_IMAGE;
		return NULL;
	}

	// fprintf(stdout, "background: %g\n", bg);
	gsl_matrix *z = gsl_matrix_alloc(area->h, area->w);
	stridefrom = fit->rx - area->w;

	// create the matrix with values from the selected rectangle
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

	result = psf_global_minimisation(z, bg, sat, com.pref.starfinder_conf.convergence, fit_angle, for_photometry, phot_set, verbose, profile, error);

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
psf_star *psf_global_minimisation(gsl_matrix* z, double bg, double sat, int convergence, gboolean fit_angle,
		gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error) {
	if (error) *error = PSF_NO_ERR;
	gboolean photometry_computed = FALSE;

	psf_star *psf = NULL;
	if (!(psf = psf_minimiz_angle(z, bg, sat, convergence, for_photometry, phot_set, verbose, profile, error))) {
		return NULL;
	}
	photometry_computed = TRUE;

	/* We quickly test the result. If it is bad we return NULL */
	if (!isfinite(psf->fwhmx) || !isfinite(psf->fwhmy) ||
			psf->fwhmx <= 0.0 || psf->fwhmy <= 0.0) {
		free_psf(psf);
		if (error && *error == PSF_NO_ERR)
			*error = PSF_ERR_DIVERGED;
		return NULL;
	}

	// Photometry
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

	return psf;
}

void psf_display_result(psf_star *result, rectangle *area) {
	char *buffer, *buffer2, *coordinates;
	char *str;
	if (com.magOffset > 0.0)
		str = _("true reduced");
	else
		str = _("relative");

	double x = result->x0 + area->x;
	double y = area->y + area->h - result->y0;

	if (has_wcs(&gfit)) {
		double world_x, world_y;
		SirilWorldCS *world_cs;
		pix2wcs(&gfit, x, (double) gfit.ry - y, &world_x, &world_y);
		world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
		if (world_cs) {
			gchar *ra = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
			gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
			coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px (%s , %s)", x, y, ra, dec);

			siril_world_cs_unref(world_cs);
			g_free(ra);
			g_free(dec);
		} else {
			coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px", x, y);
		}
	} else {
		coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px", x, y);
	}

	double fwhmx, fwhmy;
	char *unts;
	get_fwhm_as_arcsec_if_possible(result, &fwhmx, &fwhmy, &unts);

	buffer = g_strdup_printf(_("PSF fit Result:\n"
			"%s\n"
			"FWHM X=%0.2f%s, FWHM Y=%0.2f%s\n"
			"Angle=%0.2f deg\n"
			"Background value=%0.6f\n"
			"Maximal intensity=%0.6f\n"
			"Magnitude (%s)=%0.2f\n"
			"SNR=%.1fdB\n"
			"RMSE=%.3e"),
			coordinates,
			fwhmx, unts, fwhmy, unts,
			result->angle,
			result->B,
			result->A,
			str,
			result->mag + com.magOffset,
			result->SNR,
			result->rmse);
	if (result->beta > 0.0) {
		buffer2 = g_strdup_printf(_("\nMoffat fitting:\n"
									"Beta=%0.2f\n"), result->beta);
	}
	else {
		buffer2 = g_strdup_printf("\n");
	}
	buffer = g_strdup_printf("%s%s", buffer, buffer2);
	siril_log_message(buffer);
	g_free(buffer);
	g_free(buffer2);
	g_free(coordinates);
}

/* If the pixel pitch and the focal length are known and filled in the
 * setting box, we convert FWHM in pixel to arcsec by multiplying
 * the FWHM value with the sampling value */
void fwhm_to_arcsec_if_needed(fits* fit, psf_star *result) {

	if (!result) return;
	if (fit->focal_length <= 0.0 || fit->pixel_size_x <= 0.f
			|| fit->pixel_size_y <= 0.f || fit->binning_x <= 0
			|| fit->binning_y <= 0) {
		result->fwhmx_arcsec = -1.0;
		result->fwhmy_arcsec = -1.0;
		return;
	}

	double bin_X, bin_Y;

	bin_X = fit->unbinned ? (double) fit->binning_x : 1.0;
	bin_Y = fit->unbinned ? (double) fit->binning_y : 1.0;

	result->fwhmx_arcsec = result->fwhmx * (radian_conversion * (double)fit->pixel_size_x / fit->focal_length) * bin_X;
	result->fwhmy_arcsec = result->fwhmy * (radian_conversion * (double)fit->pixel_size_y / fit->focal_length) * bin_Y;
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

psf_star *new_psf_star() {
	psf_star *star = calloc(1, sizeof(psf_star));
	star->phot = NULL;

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
	return new_psf;
}

void free_psf(psf_star *psf) {
	if (psf->phot) free(psf->phot);
	free(psf);
}
