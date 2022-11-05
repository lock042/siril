#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

extern const char *filter_pattern[];
extern const size_t num_filter_patterns;
int get_cfa_pattern_index_from_string(const char *bayer);

WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int bit_depth);
int debayer(fits*, interpolation_method, sensor_pattern pattern);

#ifdef __cplusplus
extern "C" {
#endif
int adjust_Bayer_pattern(fits *fit, sensor_pattern *pattern);
WORD *debayer_buffer_superpixel_ushort(WORD *buf, int *width, int *height, sensor_pattern pattern);
float *debayer_buffer_superpixel_float(float *buf, int *width, int *height, sensor_pattern pattern);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer);
#ifdef __cplusplus
}
#endif

sensor_pattern get_bayer_pattern(fits *fit);
void clear_Bayer_information(fits *fit);

void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);

#ifdef __cplusplus
extern "C" {
#endif
/* from demosaicing_rtp.cpp */
WORD *debayer_buffer_new_ushort(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6], int bit_depth);

float *debayer_buffer_new_float(float *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]);
#ifdef __cplusplus
}
#endif

WORD *extract_CFA_buffer_ushort(fits *fit, int layer, size_t *newsize);
WORD *extract_CFA_buffer_area_ushort(fits *fit, int layer, rectangle *bounds, size_t *newsize);
float *extract_CFA_buffer_float(fits *fit, int layer, size_t *newsize);

#endif
