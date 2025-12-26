#ifndef SIRIL_MASKS_H
#define SIRIL_MASKS_H

int mask_create_test(fits *fit);
int mask_create_ones_like(fits *fit);
int mask_create_zeroes_like(fits *fit);
int mask_invert(mask_t *mask);

#endif
