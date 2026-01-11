#ifndef SIRIL_MASKS_H
#define SIRIL_MASKS_H

void mask_from_image_dialog_set_file_mode(gboolean file_mode);
int get_default_mask_bitpix();

void set_mask_active(fits *fit, gboolean state);
void free_mask(mask_t* mask);

gpointer clear_mask_worker(gpointer p);
gpointer autostretch_mask_worker(gpointer p);
gpointer binarize_mask_from_gui_worker(gpointer p);
gpointer invert_mask_worker(gpointer p);

mask_t *fits_to_mask(fits *mfit);
fits *mask_to_fits(fits *fit);
void set_poly_in_mask(UserPolygon *poly, fits *fit, gboolean state);

int mask_create_test(fits *fit, uint8_t bitpix);
int mask_create_ones_like(fits *fit, uint8_t bitpix);
int mask_create_zeroes_like(fits *fit, uint8_t bitpix);
int mask_create_from_channel(fits *fit, fits *source, int chan, uint8_t bitpix);
int mask_create_from_luminance(fits *fit, fits *source, float rw, float gw, float bw, uint8_t bitpix);
int mask_create_from_luminance_even(fits *fit, fits *source, uint8_t bitpix);
int mask_create_from_luminance_human(fits *fit, fits *source, uint8_t bitpix);
int mask_create_from_image(fits *fit, gchar *filename, int chan, uint8_t bitpix,
                           double weight_r, double weight_g, double weight_b);
int mask_create_from_stars(fits *fit, float n_fwhm, uint8_t bitpix);
int mask_create_from_chromaticity_luminance(fits *fit, fits *source,
                                            float chrom_center_r, float chrom_center_g, float chrom_center_b,
                                            float chrom_tolerance,
                                            float lum_min, float lum_max,
                                            int feather_radius, gboolean invert,
                                            uint8_t bitpix);
int mask_create_from_color_hsv(fits *fit, fits *source,
                                float h_min, float h_max,
                                float s_min, float s_max,
                                float v_min, float v_max,
                                int feather_radius, gboolean invert,
                                uint8_t bitpix);

int mask_autostretch(fits *fit);
int mask_apply_gaussian_blur(fits *fit, float radius);
int mask_binarize(fits *fit, float min_val, float max_val);
int mask_invert(fits *fit);
int mask_feather(fits *fit, float feather_dist, feather_mode mode);
int mask_scale(fits *fit, float f);
int mask_change_bitpix(fits* fit, uint8_t new_bitpix);

#endif
