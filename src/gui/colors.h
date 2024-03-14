#ifndef GUI_COLORS_H_
#define GUI_COLORS_H_

#include "algos/colors.h"

struct ccm_data {
	ccm matrix;
	float power;
	fits *fit;
	sequence *seq;
	char *seqEntry;
};

void initialize_calibration_interface();
void negative_processing();

#endif /* !GUI_COLORS_H */
