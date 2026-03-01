#ifndef SRC_FILTERS_MEDIAN_H_
#define SRC_FILTERS_MEDIAN_H_

#include <glib.h>
#include <gsl/gsl_matrix.h>

/* median filter data from GUI */
struct median_filter_data {
	destructor destroy_fn;
	fits *fit;
	int ksize;
	double amount;
	int iterations;
	gboolean previewing;
};

gpointer median_filter(gpointer p);
int median_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
gchar* median_log_hook(gpointer p, log_hook_detail detail);
void median_roi_callback();
void median_close();
double get_median_ushort(const WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);
double get_median_float(const float *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);
double get_median_gsl(gsl_matrix *mat, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self);

#endif /* SRC_FILTERS_MEDIAN_H_ */
