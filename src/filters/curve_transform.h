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
	destructor destroy_fn;  // Must be first member
	GList *points;
	enum curve_algorithm algorithm;
	gboolean do_channel[3];
	fits *fit; // just a reference, not freed
	gboolean verbose;
	gboolean for_preview;
};

typedef struct {
	double a[MAX_POINTS], b[MAX_POINTS], c[MAX_POINTS], d[MAX_POINTS];
	double x_values[MAX_POINTS], y_values[MAX_POINTS];
	int n;
} cubic_spline_data;

/* Allocator and destructor functions */
struct curve_params *new_curve_params();
void free_curve_params(void *args);

/* Core curve application function */
void apply_curve(fits *from, fits *to, struct curve_params *params, gboolean multithreaded);

/* Image processing hook for generic_image_worker */
int curve_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle functions */
gchar *curves_log_hook(gpointer p, log_hook_detail detail);
gboolean curve_preview_idle(gpointer p);
gboolean curve_apply_idle(gpointer p);

void linear_fit(GList *points, double *slopes);

float linear_interpolate(float x, GList *points, double *slopes);

void cubic_spline_fit(GList *points, cubic_spline_data *cspline_data);

float cubic_spline_interpolate(float x, cubic_spline_data *cspline_data);

#endif
