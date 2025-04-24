/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "algos/sorting.h"
#include "core/optimize_utils.h"
#include "core/proto.h"
#include "core/siril_log.h"

/************* robust 1D polynomial fit *****************/

// Helper function to perform robust fitting
static int dofit(const gsl_multifit_robust_type *T, const gsl_matrix *X, const gsl_vector *y, gsl_vector *c, gsl_matrix *cov, double *sigma, gboolean *mask) {
	gsl_multifit_robust_workspace *work = gsl_multifit_robust_alloc(T, X->size1, X->size2);
	if (!work) {
		return GSL_ENOMEM; // Return error code if memory allocation fails
	}
	int s = gsl_multifit_robust(X, y, c, cov, work);
	if (s != GSL_SUCCESS) {
		goto end_dofit;
	}
	gsl_multifit_robust_stats fit_stats = gsl_multifit_robust_statistics(work);
	gsl_vector *weights = fit_stats.weights;
	*sigma = fit_stats.sigma_rob;
	if (mask != NULL) {
		for (size_t i = 0; i < weights->size; i++) {
			mask[i] = gsl_vector_get(weights, i) > 0.0 ? TRUE : FALSE;
		}
	}

end_dofit:

	gsl_multifit_robust_free(work);
	return s;
}

// Combined function for robust polynomial fitting
int robust_polynomial_fit(double *xdata, double *ydata, int n, int degree, double *coeffs, double *uncertainties, gboolean *mask, double *sigma) {
	gsl_matrix *X = gsl_matrix_alloc(n, degree + 1);
	gsl_matrix *cov = gsl_matrix_alloc(degree + 1, degree + 1);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *c = gsl_vector_alloc(degree + 1);
	int retval;
	if (!X || !y || !c || !cov) {
		retval = GSL_ENOMEM;
		goto cleanup_and_return;
	}

	// Construct the design matrix X for polynomial fitting
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= degree; j++) {
			gsl_matrix_set(X, i, j, pow(xdata[i], j));
		}
		gsl_vector_set(y, i, ydata[i]);
	}

	// Perform robust fitting
	retval = dofit(gsl_multifit_robust_bisquare, X, y, c, cov, sigma, mask);
	if (retval != GSL_SUCCESS) {
		goto cleanup_and_return;
	}

	// Retrieve coefficients and optionally uncertainties
	for (int i = 0; i <= degree; i++) {
		coeffs[i] = gsl_vector_get(c, i);
		if (uncertainties != NULL) {
			uncertainties[i] = sqrt(gsl_matrix_get(cov, i, i));
		}
	}

cleanup_and_return:

	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);

	return retval;
}

/**
 * Siegel's Repeated Median Regression for linear fit
 * Extremely robust method that can handle up to 50% contamination
 */
int repeated_median_fit(double *xdata, double *ydata, int n, double *a,
		double *b, double *sigma, gboolean *mask) {

	if (n < 2) {
		return -1; // Need at least 2 points for a line
	}

	// Arrays to hold all possible pairwise slopes
	double *all_slopes = g_malloc(n * (n - 1) / 2 * sizeof(double));
	double *point_medians = g_malloc(n * sizeof(double)); // For each point, median of slopes with all other points

	// For each point, calculate slopes with all other points
	for (int i = 0; i < n; i++) {
		int slope_count = 0;
		double *slopes_for_point = g_malloc((n - 1) * sizeof(double));

		for (int j = 0; j < n; j++) {
			if (i != j) {
				// Avoid division by zero
				if (xdata[i] != xdata[j]) {
					slopes_for_point[slope_count] = (ydata[j] - ydata[i])
							/ (xdata[j] - xdata[i]);
					slope_count++;
				}
			}
		}

		// Calculate median of slopes for this point
		point_medians[i] = quickmedian_double(slopes_for_point, slope_count);

		g_free(slopes_for_point);
	}

	// The final slope is the median of all the point medians
	double slope = quickmedian_double(point_medians, n);

	// Calculate intercepts for each point
	double *intercepts = g_malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		intercepts[i] = ydata[i] - slope * xdata[i];
	}

	// The final intercept is the median of all intercepts
	double intercept = quickmedian_double(intercepts, n);

	// Set output parameters
	*b = slope;
	*a = intercept;

	// Calculate residuals and identify outliers
	double residual_sum = 0.0;
	int inlier_count = 0;

	double *absolute_residuals = g_malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		double predicted = intercept + slope * xdata[i];
		double residual = ydata[i] - predicted;
		absolute_residuals[i] = fabs(residual);
		residual_sum += residual * residual;
	}

	// Calculate MAD (Median Absolute Deviation) for residuals
	double mad = quickmedian_double(absolute_residuals, n);

	// Normalize by 1.4826 to make MAD equivalent to standard deviation for normal distribution
	mad *= MAD_NORM;

	// Mark points as inliers/outliers using MAD
	for (int i = 0; i < n; i++) {
		if (mask != NULL) {
			// Points with absolute residual > 3 * MAD are considered outliers
			// (TODO: This threshold can be adjusted (common values: 2.5-3.5))
			mask[i] = (absolute_residuals[i] <= 3.0 * mad);
		}

		if (absolute_residuals[i] <= 3.0 * mad) {
			inlier_count++;
		}
	}

	// Calculate sigma (standard deviation of residuals)
	if (sigma != NULL) {
		if (inlier_count > 2) {
			*sigma = sqrt(residual_sum / inlier_count);
		} else {
			*sigma = sqrt(residual_sum / n);
		}
	}

	// Cleanup
	g_free(all_slopes);
	g_free(point_medians);
	g_free(intercepts);
	g_free(absolute_residuals);

	return 0;
}

// Function to evaluate a 1D polynomial at a given point
double evaluate_polynomial(double *coeffs, int degree, double x) {
	double result = 0.0;
	for (int i = degree; i >= 0; i--) {
		result = result * x + coeffs[i];
	}
	return result;
}

// Functions used to find linear fit coefficients for fitting an entire fit to a reference fit
// These are used by linear match
// Robust fitting is not used here because it is far too slow
static int find_linear_coeff_ushort(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	double c0, c1, cov00, cov01, cov11, sumsq;
	size_t ref_size = reference_fit->rx * reference_fit->ry;

	if (memcmp(target_fit->naxes, reference_fit->naxes, sizeof target_fit->naxes)) {
		gchar *err = siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		if (error) {
			*error = err;
		}
		return -1;
	}

	low *= USHRT_MAX_DOUBLE;
	high *= USHRT_MAX_DOUBLE;

	siril_log_color_message(_("Linear fit functions:\n"), "green");
	for (int channel = 0; channel < reference_fit->naxes[2]; channel++) {
		ssize_t j = 0;
		double *x = malloc(ref_size * sizeof(double));
		double *y = malloc(ref_size * sizeof(double));
		for (size_t i = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->pdata[channel][i], low, high) && (reference_fit->pdata[channel][i] != 0)) {
				if (target_fit->type == DATA_FLOAT && (target_fit->fpdata[channel][i] != 0)) {
					x[j] = (double) target_fit->fpdata[channel][i];
				} else if (target_fit->type != DATA_FLOAT && target_fit->pdata[channel][i] != 0) {
					x[j] = (double) target_fit->pdata[channel][i] * INV_USHRT_MAX_DOUBLE;
				} else 	continue;
				y[j] = (double) reference_fit->pdata[channel][i] * INV_USHRT_MAX_DOUBLE;
				j++;
			}
		}
		j--;
		if (j > 1) {
			gsl_fit_linear(x, 1, y, 1, (int)j, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
			siril_log_color_message("y_0 = %e + %e*x_0 (%d)\n", "blue", c0, c1, (int)j);
			free(x);
			free(y);
			a[channel] = c1;
			b[channel] = c0;
		} else {
			gchar *err = siril_log_color_message(_("Error! Need at least 2 points...\n"), "red");
			if (error) {
				*error = err;
			}
			free(x);
			free(y);
			return -1;
		}
	}
	return 0;
}

static int find_linear_coeff_float(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	double c0, c1, cov00, cov01, cov11, sumsq;
	size_t ref_size = reference_fit->rx * reference_fit->ry;

	if (memcmp(target_fit->naxes, reference_fit->naxes, sizeof target_fit->naxes)) {
		gchar *err = siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		if (error) {
			*error = err;
		}
		return -1;
	}

	siril_log_color_message(_("Linear fit functions:\n"), "green");
	for (int channel = 0; channel < reference_fit->naxes[2]; channel++) {
		ssize_t j = 0;
		double *x = malloc(ref_size * sizeof(double));
		double *y = malloc(ref_size * sizeof(double));
		for (size_t i = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->fpdata[channel][i], low, high) && (reference_fit->fpdata[channel][i] != 0)) {
				if (target_fit->type == DATA_FLOAT && (target_fit->fpdata[channel][i] != 0)) {
					x[j] = (double) target_fit->fpdata[channel][i];
				} else if (target_fit->type != DATA_FLOAT && target_fit->pdata[channel][i] != 0) {
					x[j] = (double) target_fit->pdata[channel][i] * INV_USHRT_MAX_DOUBLE;
				} else 	continue;
				y[j] = (double) reference_fit->fpdata[channel][i];
				j++;
			}
		}
		j--;
		if (j > 1) {
			gsl_fit_linear(x, 1, y, 1, (int)j, &c0, &c1, &cov00, &cov01, &cov11,	&sumsq);
			siril_log_color_message("y_0 = %e + %e*x_0 (%d)\n", "blue", c0, c1, (int)j);
			free(x);
			free(y);
			a[channel] = c1;
			b[channel] = c0;
		} else {
			gchar *err = siril_log_color_message(_("Error! Need at least 2 points...\n"), "red");
			if (error) {
				*error = err;
			}
			free(x);
			free(y);
			return -1;
		}
	}
	return 0;
}

int find_linear_coeff(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	if (reference_fit->type == DATA_USHORT) {
		return find_linear_coeff_ushort(target_fit, reference_fit, low, high, a, b, error);
	} else if (reference_fit->type == DATA_FLOAT) {
		return find_linear_coeff_float(target_fit, reference_fit, low, high, a, b, error);
	}
	return 1;
}
