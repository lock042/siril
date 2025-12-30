#ifndef SIRIL_MASKS_H
#define SIRIL_MASKS_H

void free_mask(mask_t* mask);

int mask_create_test(fits *fit, uint8_t bitpix);
int mask_create_ones_like(fits *fit, uint8_t bitpix);
int mask_create_zeroes_like(fits *fit, uint8_t bitpix);
int mask_create_from_channel(fits *fit, int chan, uint8_t bitpix);
int mask_create_from_luminance(fits *fit, float rw, float gw, float bw, uint8_t bitpix);
int mask_create_from_luminance_even(fits *fit, uint8_t bitpix);
int mask_create_from_luminance_human(fits *fit, uint8_t bitpix);

int mask_apply_gaussian_blur(fits *fit, float radius);
int mask_binarize(fits *fit, float min_val, float max_val);
int mask_invert(fits *fit);

#endif
