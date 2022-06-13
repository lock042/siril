#ifndef SRC_GUI_PAYNE_H_
#define SRC_GUI_PAYNE_H_

int paynelut(fits *fit, double beta, double intensity, double lower, double shoulder, double headroom, double offset, int colourmodel, int stretchtype);
void apply_payne_cancel();

#define STRETCH_LINEAR 0
#define STRETCH_PAYNE_NORMAL 1
#define STRETCH_PAYNE_INVERSE 2
#define STRETCH_ASINH 3
#define STRETCH_INVASINH 4

#define COL_HUMANLUM 0
#define COL_EVENLUM 1
#define COL_INDEP 2

#endif /* SRC_GUI_PAYNE_H_ */
