#ifndef _BLENDING_H
#define _BLENDING_H

#include "stacking.h"

void init_ramp();
float get_ramped_value(float val);

void compute_downscaled_mask_size(int rx, int ry, int *rx_out, int *ry_out, double *fx, double *fy);

int compute_masks(struct stacking_args *args);

#endif
