#ifndef SRC_GUI_BILAT_H_
#define SRC_GUI_BILAT_H_

typedef enum {
	EP_BILATERAL,
	EP_GUIDED
} ep_filter_t;

int edge_preserving_filter(fits *fit, fits *guide, double d, double sigma_col, double sigma_space, ep_filter_t filter_type, gboolean verbose);
void bilat_change_between_roi_and_image();
void apply_bilat_cancel();

#endif /* SRC_GUI_ASINH_H_ */
