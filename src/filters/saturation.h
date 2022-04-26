#ifndef SRC_GUI_SATURATION_H_
#define SRC_GUI_SATURATION_H_

#include <glib.h>

/* color saturation data from GUI */
struct enhance_saturation_data {
	fits *input, *output;
	double coeff, h_min, h_max, background_factor;
	gboolean for_preview;
};

void apply_satu_cancel();
gpointer enhance_saturation(gpointer p);
void satu_set_hues_from_types(struct enhance_saturation_data *args, int type);

#endif /* SRC_GUI_SATURATION_H_ */
