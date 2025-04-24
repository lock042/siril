#ifndef _SIRIL_FITTING_H
#define _SIRIL_FITTING_H

int robust_polynomial_fit(double *xdata, double *ydata, int n, int degree, double *coeffs, double *uncertainties, gboolean *mask, double *sigma);
int robust_linear_fit(double *xdata, double *ydata, int n, double *a, double *b, double *sigma, gboolean *mask);
int repeated_median_fit(double *xdata, double *ydata, int n, double *a, double *b, double *sigma, gboolean *mask);
double evaluate_polynomial(double *coeffs, int degree, double x);

int find_linear_coeff(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error);

#endif
