#ifndef SRC_FILTERS_SATURATION_H_
#define SRC_FILTERS_SATURATION_H_

#include <glib.h>
#include "core/processing.h"

typedef struct {
	void (*free)(gpointer); /* Destructor required by generic_img_args */
	double coeff;
	double background_factor;
	double h_min;
	double h_max;
} saturation_params;

gchar* satu_log_hook(gpointer p, log_hook_detail detail);
void satu_set_hues_from_types(saturation_params *args, int type);
int saturation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
/* GUI functions (apply_satu_cancel, satu_change_between_roi_and_image) declared in gui/saturation.h */

#endif /* SRC_FILTERS_SATURATION_H_ */
