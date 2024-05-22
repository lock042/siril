#ifndef _SIRIL_FITTING_H
#define _SIRIL_FITTING_H

int robust_linear_fit(double *xdata, double *ydata, int n, double *a, double *b, double *sigma, gboolean *mask);
double evaluate_polynomial(double *coeffs, int degree, double x);
void ransac_polynomial_fit(double *x, double *y, int n, int degree, double *best_coeffs);

#endif
