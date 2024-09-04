#ifndef SIRIL_CURVE_TRANSFORM_H
#define SIRIL_CURVE_TRANSFORM_H

#include "core/siril.h"
#include <glib.h>

#define MAX_POINTS 100

// Order must match the drop-down menu in the curve window
enum curve_algorithm {
	CUBIC_SPLINE,
	LINEAR
};

struct curve_params {
	GList *points;
	enum curve_algorithm algorithm;
	gboolean do_channel[3];
};

typedef struct {
	double a[MAX_POINTS], b[MAX_POINTS], c[MAX_POINTS], d[MAX_POINTS];
	double x_values[MAX_POINTS], y_values[MAX_POINTS];
	int n;
} cubic_spline_data;

void apply_curve(fits *from, fits *to, struct curve_params params, gboolean mutlithreaded);

void linear_fit(GList *points, double *slopes);

float linear_interpolate(float x, GList *points, double *slopes);

void cubic_spline_fit(GList *points, cubic_spline_data *cspline_data);

float cubic_spline_interpolate(float x, cubic_spline_data *cspline_data);

#endif