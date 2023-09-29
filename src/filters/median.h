#ifndef SRC_FILTERS_MEDIAN_H_
#define SRC_FILTERS_MEDIAN_H_

#include <glib.h>
#include <gsl/gsl_matrix.h>

/* median filter data from GUI */
struct median_filter_data {
	fits *fit;
	int ksize;
	double amount;
	int iterations;
	gboolean previewing;
};

gpointer median_filter(gpointer p);

void median_roi_callback();
double get_median_ushort(const WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);
double get_median_float(const float *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);
double get_median_gsl(gsl_matrix *mat, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);

#endif /* SRC_FILTERS_MEDIAN_H_ */
