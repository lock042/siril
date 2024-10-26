#ifndef _BLENDING_H
#define _BLENDING_H

#include "stacking.h"

double get_mask_scale();
void compute_downscaled_mask_size(int rx, int ry, int *rx_out, int *ry_out, double *fx, double *fy);

gchar *get_mask_filename(sequence *seq, int index);
int compute_masks(struct stacking_args *args);

#endif
