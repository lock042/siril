#ifndef SRC_GUI_PAYNE_H_
#define SRC_GUI_PAYNE_H_

int paynelut(fits *fit, double beta, double intensity, double lower, double shoulder, double headrom, double offset, gboolean rgb_space, gboolean inverse);
void apply_payne_cancel();

#define STRETCH_PAYNE_NORMAL 1
#define STRETCH_PAYNE_INVERSE 2
#define STRETCH_ASINH_NORMAL 3
#define STRETCH_ASINH_INVERSE 4

#endif /* SRC_GUI_PAYNE_H_ */
