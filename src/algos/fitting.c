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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include "core/siril.h"
#include "core/proto.h"
#include "algos/noise.h"
#include "core/processing.h"
#include "core/siril_log.h"

/************* robust linear fit *****************/

static int dofit(const gsl_multifit_robust_type *T, const gsl_matrix *X, const gsl_vector *y, gsl_vector *c, gsl_matrix *cov, double *sigma, gboolean *mask) {
	gsl_multifit_robust_workspace *work = gsl_multifit_robust_alloc (T, X->size1, X->size2);
	int s = gsl_multifit_robust (X, y, c, cov, work);
	gsl_multifit_robust_stats fit_stats = gsl_multifit_robust_statistics(work);
	gsl_vector *weights = fit_stats.weights;
	*sigma = fit_stats.sigma_rob;
	for (size_t i = 0; i < weights->size; i++) {
		if (gsl_vector_get(weights, i) > 0.) // weights are null if gsl_multifit_robust_type = gsl_multifit_robust_bisquare
			mask[i] = TRUE;
	}
	gsl_multifit_robust_free (work);
	return s;
}

// Fits y = a + bx to double *xdata, double *ydata with each array of length n
int robust_linear_fit(double *xdata, double *ydata, int n, double *a, double *b, double *sigma, gboolean *mask) {
	const size_t p = 2; /* linear fit */
	gsl_matrix *X = NULL, *cov = NULL;
	gsl_vector *y = NULL, *c = NULL;
	X = gsl_matrix_alloc (n, p);
	y = gsl_vector_alloc (n);
	c = gsl_vector_alloc (p);
	cov = gsl_matrix_alloc (p, p);
	if (!X || !y || !c || !cov) {
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(c);
		gsl_matrix_free(cov);
		return 1;
	}
	for (int i = 0 ; i < n ; i++) {
		gsl_vector_set(y, i, ydata[i]);
	}
	/* construct design matrix X for linear fit */
	for (int i = 0; i < n; ++i) {
		gsl_matrix_set (X, i, 0, 1.0);
		gsl_matrix_set (X, i, 1, xdata[i]);
	}

	/* perform robust and OLS fit */
	int retval = dofit(gsl_multifit_robust_bisquare, X, y, c, cov, sigma, mask);
	*a = gsl_vector_get(c,0);
	*b = gsl_vector_get(c,1);
	gsl_matrix_free (X);
	gsl_vector_free (y);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
	return retval;
}

// Function to compute threshold based on FITS noise sigma. Threshold is 1 * sigma,
// which is strict but provides good fit accuracy.

#define THRESHOLD_SIGMA_MULTIPLIER 1.0

static double compute_threshold(fits *fit) {
	double retval;
	struct noise_data *data = calloc(1, sizeof(struct noise_data));
    if (data == NULL) {
        PRINT_ALLOC_ERR;
        return NAN;
    }
	data->fit = fit;
	// All other members are FALSE / 0.0 anyway because of calloc()
	noise_worker(data);
	if (gfit.naxes[2] == 1)
		retval = data->bgnoise[0];
	else
		retval = sqrt(	data->bgnoise[0] * data->bgnoise[0] +
						data->bgnoise[1] * data->bgnoise[1] +
						data->bgnoise[2] * data->bgnoise[2]);
	free(data);
	return THRESHOLD_SIGMA_MULTIPLIER * retval;
}

// Function to perform polynomial fitting using GSL
static void gsl_polynomial_fit(double *x, double *y, int n, int degree, double *coeffs, double *uncertainties) {
	gsl_matrix *X, *cov;
	gsl_vector *y_vec, *c;
	double chisq;

	X = gsl_matrix_alloc(n, degree + 1);
	y_vec = gsl_vector_alloc(n);
	c = gsl_vector_alloc(degree + 1);
	cov = gsl_matrix_alloc(degree + 1, degree + 1);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= degree; j++) {
			gsl_matrix_set(X, i, j, pow(x[i], j));
		}
		gsl_vector_set(y_vec, i, y[i]);
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, degree + 1);
	gsl_multifit_linear(X, y_vec, c, cov, &chisq, work);
	gsl_multifit_linear_free(work);

	for (int i = 0; i <= degree; i++) {
		coeffs[i] = gsl_vector_get(c, i);
		if (uncertainties != NULL) {
			uncertainties[i] = sqrt(gsl_matrix_get(cov, i, i)); // Standard deviation (uncertainty)
		}
	}

	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y_vec);
	gsl_vector_free(c);
}

// Function to evaluate the polynomial at a given point
double evaluate_polynomial(double *coeffs, int degree, double x) {
	double result = 0.0;
	for (int i = degree; i >= 0; i--) {
		result = result * x + coeffs[i];
	}
	return result;
}

// RANSAC for polynomial fit
void ransac_polynomial_fit(double *x, double *y, int n, int degree, double *best_coeffs, double *threshold, int max_iters, fits *fit, double *best_uncertainties) {
	int max_inliers = 0;
	double best_error = INFINITY;
	double *temp_coeffs = (double *)malloc((degree + 1) * sizeof(double));
	int *inliers = (int *)malloc(n * sizeof(int));
	if (!temp_coeffs || !inliers) {
		PRINT_ALLOC_ERR;
		free(temp_coeffs);
		free(inliers);
		return;
	}
	*threshold = compute_threshold(fit);
	for (int iter = 0; iter < max_iters; iter++) {
		// Randomly select a subset of points
		int subset_size = 3 * degree + 1;
		int subset_indices[subset_size];
		for (int i = 0; i < subset_size; i++) {
			subset_indices[i] = rand() % n;
		}

		double subset_x[subset_size];
		double subset_y[subset_size];
		for (int i = 0; i < subset_size; i++) {
			subset_x[i] = x[subset_indices[i]];
			subset_y[i] = y[subset_indices[i]];
		}

		// Fit polynomial to subset
		gsl_polynomial_fit(subset_x, subset_y, subset_size, degree, temp_coeffs, NULL);

		// Count inliers
		int num_inliers = 0;
		double error_sum = 0.0;
		for (int i = 0; i < n; i++) {
			double y_pred = evaluate_polynomial(temp_coeffs, degree, x[i]);
			double error = fabs(y[i] - y_pred);
			if (error < *threshold) {
				inliers[num_inliers++] = i;
				error_sum += error;
			}
		}

		// Update best model if better
		if (num_inliers > max_inliers || (num_inliers == max_inliers && error_sum < best_error)) {
			max_inliers = num_inliers;
			best_error = error_sum;
			for (int i = 0; i <= degree; i++) {
				best_coeffs[i] = temp_coeffs[i];
			}
		}
	}

	// Compute uncertainties for the best model
	double *best_fit_x = (double *)malloc(max_inliers * sizeof(double));
	double *best_fit_y = (double *)malloc(max_inliers * sizeof(double));
	if (best_fit_x == NULL || best_fit_y == NULL) {
		PRINT_ALLOC_ERR;
		goto cleanup;
	}
	for (int i = 0; i < max_inliers; i++) {
		best_fit_x[i] = x[inliers[i]];
		best_fit_y[i] = y[inliers[i]];
	}
	gsl_polynomial_fit(best_fit_x, best_fit_y, max_inliers, degree, best_coeffs, best_uncertainties);

cleanup:

	free(temp_coeffs);
	free(inliers);
	free(best_fit_x);
	free(best_fit_y);
}
