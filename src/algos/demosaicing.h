#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

extern const char *filter_pattern[];
extern const size_t num_filter_patterns;
sensor_pattern get_cfa_pattern_index_from_string(const char *bayer);
sensor_pattern get_validated_cfa_pattern(fits *fit, gboolean force_debayer);
int get_compiled_pattern(fits *fit, BYTE pattern[36], int *pattern_size);
int FC_array(int row, int col, BYTE* bpattern, int size);
gboolean compare_compiled_pattern(BYTE *refpattern, BYTE *pattern, int pattern_size);

int debayer(fits*, interpolation_method, sensor_pattern pattern);

#ifdef __cplusplus
extern "C" {
#endif
WORD *debayer_buffer_superpixel_ushort(WORD *buf, int *width, int *height, sensor_pattern pattern);
float *debayer_buffer_superpixel_float(float *buf, int *width, int *height, sensor_pattern pattern);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer);
#ifdef __cplusplus
}
#endif

struct merge_cfa_data {
	sequence *seq0;
	sequence *seq1;
	sequence *seq2;
	sequence *seq3;
	char *seqEntryOut;
	sensor_pattern pattern;
};
void update_bayer_pattern_information(fits *fit, sensor_pattern pattern);
void apply_mergecfa_to_sequence(struct merge_cfa_data *merge_cfa_args);

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

#ifdef __cplusplus
extern "C" {
#endif
fits* merge_cfa (fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3, sensor_pattern pattern);
#ifdef __cplusplus
}
#endif
WORD *extract_CFA_buffer_ushort(fits *fit, int layer, size_t *newsize);
WORD *extract_CFA_buffer_area_ushort(fits *fit, int layer, rectangle *bounds, size_t *newsize);
float *extract_CFA_buffer_float(fits *fit, int layer, size_t *newsize);
float *extract_CFA_buffer_area_float(fits *fit, int layer, rectangle *bounds, size_t *newsize);

#endif
