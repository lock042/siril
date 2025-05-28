#ifndef _GPS_PARSER_H
#define _GPS_PARSER_H

#include "core/processing.h"

struct _qhy_struct {
	// native endianness
	uint32_t sequence_number;
	uint16_t image_width;
	uint16_t image_height;
	double latitude;	// format is SDDMMMMMMM (sign, degrees, minutes)
	double longitude;
	BYTE start_flag;	// GPS status, bits 0 (PPS), 1 (second lock) and 4+5 (recv flag)
	BYTE end_flag;
	BYTE now_flag;
	GDateTime *start;
	GDateTime *end;
	GDateTime *now;
	uint32_t count_of_PPS;
	gboolean flags_are_shifted;	// in the metadata they're not, in NINA headers they are
};

enum timestamp_type { EXP_START, EXP_MIDDLE, EXP_END };

int parse_gps_image(fits *fit, struct _qhy_struct *qhy_header);
int parse_gps_from_header(fits *fit, const char *filename, struct _qhy_struct *qhy_header);
int update_fit_from_qhy_header(fits *fit, struct _qhy_struct *qhy_header);

void print_qhy_data(struct _qhy_struct *qhy);
void release_qhy_struct(struct _qhy_struct *data);

int get_rs_gps_data(fits *fit, struct gps_rs_data *data);
GDateTime *get_timestamp_for_pixel(struct gps_rs_data *data, enum timestamp_type type, int x, int y);
struct gps_rs_data *clone_gps_data(struct gps_rs_data *in);

void apply_crop_to_gps_data(fits *fit, rectangle *bounds);
void apply_flip_to_gps_data(fits *fit);
void apply_binning_to_gps_data(fits *fit);

int gps_extract_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
                rectangle *_, int threads);
#endif
