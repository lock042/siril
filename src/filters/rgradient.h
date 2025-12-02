#ifndef SRC_ALGOS_RGRADIENT_H_
#define SRC_ALGOS_RGRADIENT_H_

#include "core/siril.h"
#include "core/processing.h"

/* rgradient filter data */
struct rgradient_data {
	destructor destroy_fn;  // Must be first member
	fits *fit;  // just a reference, not freed
	double xc, yc, dR, da;
	gboolean verbose;
};

/* Allocator and destructor functions */
struct rgradient_data *new_rgradient_data();
void free_rgradient_data(void *args);

/* Image processing hook */
int rgradient_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle function */
gboolean rgradient_idle(gpointer p);
gchar *rgradient_log_hook(struct generic_img_args *p);

#endif /* SRC_ALGOS_RGRADIENT_H_ */
