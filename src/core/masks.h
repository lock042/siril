#ifndef SIRIL_MASKS_H
#define SIRIL_MASKS_H

void free_mask(mask_t* mask);
int mask_create_test(fits *fit, uint8_t bitpix);
int mask_create_ones_like(fits *fit, uint8_t bitpix);
int mask_create_zeroes_like(fits *fit, uint8_t bitpix);
int mask_invert(mask_t *mask);

#endif
