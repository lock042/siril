#ifndef _SIRIL_STATS_FLOAT_H
#define _SIRIL_STATS_FLOAT_H

#include "core/siril.h"

imstats* statistics_internal_float(fits *fit, int layer, rectangle *selection,
		int option, imstats *stats, int bitpix, threading_type threads);
int compute_means_from_flat_cfa_float(fits *fit, double mean[36]);
#if 0
int IKSS(float *data, size_t n, double *location, double *scale, gboolean multithread);
#endif
int IKSSlite(float *data, size_t n, const float median, float mad, double *location, double *scale, threading_type threads);

float siril_stats_robust_mean(const float sorted_data[],
		const size_t stride, const size_t size, double *deviation);

int robustmean(int n, const double *x, double *mean, double *stdev);
int robust_linear_fit(double *xdata, double *ydata, int n, double *a, double *b, double *sigma, gboolean *mask);

#endif
