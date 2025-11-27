#ifndef SRC_GUI_SATURATION_H_
#define SRC_GUI_SATURATION_H_

#include <glib.h>

typedef struct {
	void (*free)(gpointer); // Destructor required by generic_img_args
	double coeff;
	double background_factor;
	double h_min;
	double h_max;
} saturation_params;

void satu_change_between_roi_and_image();
void apply_satu_cancel();
void satu_set_hues_from_types(saturation_params *args, int type);
int saturation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

#endif /* SRC_GUI_SATURATION_H_ */
